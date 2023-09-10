import os
import sys
import time
from datetime import datetime
import networkx as nx

from sklearn.svm import SVC

import grakel as gko 
import grakel.kernels as gkr 

import numpy as np

from sklearn.base import clone
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

from kernel_grid_search import KernelGridSearchCV

workflow_path = os.environ.get('path_workflow_screendop')
if(workflow_path==None):
    raise Exception('You forgot to setup the path_workflow environment variable with the path to the workflow')
    
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]
sys.path.insert(0, workflow_path)
from utils import process_memory

# Step 1
class Generate_personalised_networks:
    def __init__(self, cutoff_up, cutoff_down):
        self.cutoff_up = cutoff_up
        self.cutoff_down = cutoff_down

    def load_interactome(self):
        interactome=[]

        f=open( f"{workflow_path}/interactome_hint.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            interactome.append([l[0], l[1]])
        f.close()

        return interactome

    def load_expression_table(self, folder, expr_table_file):
        expr_data={}
        c=0
        f=open(folder+'/'+expr_table_file,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(c==0):
                for sample in l[1:]:
                    expr_data[sample]={}
            else:
                j=1
                for sample in expr_data.keys():
                    value=0.0
                    if(j<len(l)):
                        value=float(l[j])
                    expr_data[sample][l[0]]=value
                    j+=1
            c+=1
        f.close()

        return expr_data

    def _get_node_label(self, exp_value):
        cutup = self.cutoff_up
        cutdown = self.cutoff_down
        
        node_label='neutral'
        node_nlabel = 0
        if ( exp_value >= cutup ):
            node_label='up'
            node_nlabel = 1
        elif ( exp_value <= cutdown ):
            node_label='down'
            node_nlabel = 2 
        return node_label, node_nlabel
    
    def _get_edge_label(self, exp_value_p1, exp_value_p2):
        cutup = self.cutoff_up
        cutdown = self.cutoff_down
        
        edge_label = 'no_gde'
        edge_nlabel = 0
        if( ( exp_value_p1 >= cutup and ( exp_value_p2 > cutdown and exp_value_p2 < cutup ) ) or ( exp_value_p1 > cutdown and exp_value_p1 < cutup and exp_value_p1 >= cutup ) ) : 
            edge_label = 'mixed'
            edge_nlabel = 1
        elif( exp_value_p1 >= cutup and exp_value_p2 >= cutup ):
            edge_label = 'complex_up'
            edge_nlabel = 2
        elif( exp_value_p1 <= cutdown and exp_value_p2 <= cutdown ):
            edge_label = 'complex_down'
            edge_nlabel = 3
        return edge_label, edge_nlabel
    
    def make_patients_graph(self, folder, mask_expr_table_file, raw_expr_table_file):
        print("Loading interactome...")
        positive_pairs=self.load_interactome()
        before = len(positive_pairs)
        
        print("Loading mask expression table...")
        expr_data=self.load_expression_table(folder, mask_expr_table_file)
        
        print("Loading raw expression table...")
        raw_expr_data=self.load_expression_table(folder, raw_expr_table_file)
        
        samples = list(expr_data.keys())
        proteins = set( expr_data[samples[0]] )
        positive_pairs = list( filter( lambda x: ( x[0] in proteins and x[1] in proteins ), positive_pairs ) )
        print( '\tfiltered interactome from ', before, 'to', len(positive_pairs) )

        print("Generating the patients graphs and exporting to graphml...")
        if(not os.path.isdir(folder+"/patient_graphs")):
            os.system("mkdir "+folder+"/patient_graphs")
        #else:
        #    os.system("rm "+folder+"/patient_graphs/*")

        oneorother=0
        both=0
        i=1
        for k in samples:
            #print(i, '/', len(samples) )
            if(k!=""):
                #print('proteins under -1:', len( list( filter( lambda x: x<-1, list(expr_data[k].values())) ) ) )
                #print('proteins above 1:', len( list( filter( lambda x: x>1, list(expr_data[k].values())) ) ) )
                   
                gr=nx.Graph()
                nodes=set()
                edges=set()
                for p in positive_pairs:
                    p1=p[0]
                    p2=p[1]
                    if(p1 in expr_data[k].keys() and p2 in expr_data[k].keys()):
                        edge_label, edge_nlabel = self._get_edge_label(expr_data[k][p1], expr_data[k][p2])
                        #if( not p1+'-'+p2 in edges and edge_nlabel!=0 ):
                        if( not p1+'-'+p2 in edges and edge_nlabel > 1 ):
                            node_label, node_nlabel = self._get_node_label( raw_expr_data[k][p1] )
                            if( not p1 in nodes ):
                                gr.add_node(p1, label=node_nlabel, ntype=node_label, expression=raw_expr_data[k][p1] )
                                #gr.add_node(p1, label=node_nlabel )
                                nodes.add(p1)
                            
                            node_label, node_nlabel = self._get_node_label( raw_expr_data[k][p2] )
                            if( not p2 in nodes ):
                                gr.add_node(p2, label=node_nlabel, ntype=node_label, expression=raw_expr_data[k][p2] )
                                #gr.add_node(p2, label=node_nlabel )
                                nodes.add(p2)
                        
                            w = ( raw_expr_data[k][p1] + raw_expr_data[k][p2]) / 2
                            if( edge_nlabel == 2 ):
                                w = abs(w)
                            elif( edge_nlabel == 3 ):
                                w = abs(w)*(-1)
                            gr.add_edge(p1, p2, label=edge_nlabel, etype=edge_label, weight=w )
                            #gr.add_edge(p1, p2, label=edge_nlabel )
                            edges.add( p1+'-'+p2 )
                            
                if( len(gr.edges)>0 ):
                    nx.write_graphml(gr, folder+"/patient_graphs/graph_"+k+".graphml")
                        
                i+=1
        #print('both: ', both)
        #print('one or the other:', oneorother)
        
# Step 2
class KernelMatrixBasedPrediction:
    def __init__(self, folder, nodes_enrich, edges_enrich, identifier):
        self.modeNode = nodes_enrich
        self.modeEdge = edges_enrich
        fout = f"{folder}/{identifier}_nEnrich-{nodes_enrich}_eEnrich-{edges_enrich}"
        self.fout = fout
        if(not os.path.isdir( self.fout )):
            os.system( f"mkdir {self.fout}")
        
    def _get_attr_node(self, g, id_):
        attrs = list( filter( lambda y: y!=None, [ x if( not isinstance(x, str) ) else None for x in list(g.nodes[id_].values()) ] ) )
        return attrs
        
    def _get_attr_edge(self, g, id_):
        attrs = list( filter( lambda y: y!=None, [ x if( not isinstance(x, str) ) else None for x in list(g.edges[id_].values()) ] ) )
        return attrs
    
    def load_dataset_networks_grakel(self, folder, labels_file):
        data_directory = folder+'/patient_graphs'
        y_real={}
        with open(folder+'/'+labels_file) as f:
            for line in f:
                l=line.replace("\n","").split("\t")
                y_real[l[0]]=float(l[1])
        
        #Xm1 = []
        #Xm2 = []
        
        X = []
        y = []
        samples = []
        for filename in os.listdir(data_directory):
            full_path=os.path.join(data_directory, filename)
            
            sample=filename.replace("graph_","").replace(".graphml","")
            samples.append(sample)
            
            y.append(y_real[sample])
            #print(full_path)
            
            #edgesm1 = {}
            #edgesm2 = {}
            gr = nx.read_graphml( full_path )
            nodesAttr = {}
            edgesAttr = {}
            for ed in list(gr.edges):
                p1 = ed[0]
                exp = gr.nodes[p1]['expression']
                lb = gr.nodes[p1]['label']
                
                if( self.modeNode=='label' ):
                    nodesAttr[p1] = [lb]
                elif( self.modeNode=='weight' ):
                    nodesAttr[p1] = [exp]
                else:
                    nodesAttr[p1] = [lb, exp]
                 
                p2 = ed[1]
                exp = gr.nodes[p2]['expression']
                lb = gr.nodes[p2]['label']
                
                if( self.modeNode=='label' ):
                    nodesAttr[p2] = [lb]
                elif( self.modeNode=='weight' ):
                    nodesAttr[p2] = [exp]
                else:
                    nodesAttr[p2] = [lb, exp]
                
                w = gr.edges[ ed ]['weight']
                lb = gr.edges[ ed ]['label']
                #if( not p1 in edgesm1 ):
                #     edgesm1[p1]=[]
                #edgesm1[p1].append(p2)
                #edgesm2[ed] = lb
                
                if( self.modeEdge=='label' ):
                    edgesAttr[ed] = [lb]
                elif( self.modeEdge=='weight' ):
                    edgesAttr[ed] = [w]
                else:
                    edgesAttr[ed] = [lb, w]
                
            #Xm1.append( gko.Graph(edgesm1, node_labels = nodesAttr) )
            #Xm2.append( gko.Graph(edgesm2, edge_labels = edgesAttr) )
            #X.append( [ set(gr.edges), { n: self._get_attr_node(gr, n) for n in gr.nodes }, { e: self._get_attr_edge(gr, e) for e in gr.edges } ] )
            X.append( [ set(gr.edges), nodesAttr, edgesAttr ] )
            
        y = np.array(y)
        
        #return [Xm1, Xm2, X, y]
        return [X, y, samples]

    def balance_dataset(self, x, y):
        i=0
        pos=[]
        neg=[]
        for el in y:
            if(el == 0):
                neg.append(i)
            else:
                pos.append(i)
            i+=1
        
        n = len( pos )
        shrink = neg
        rest = pos
        if( len(pos) > len(neg) ):
            n = len(neg)
            shrink = pos
            rest = neg
        print('before downsampling', len(matrix) )
        index = np.random.choice( shrink, size=n, replace=False)
        index = list(index)+rest
        Xbal = []
        ybal=[]
        i=0
        for i in index:
            ybal.append(y[i])
            Xbal.append(x[i, index])
            print('after downsampling', len(Xbal) )
            Xbal = np.array(Xbal)
            ybal = np.array(ybal)
        return Xbal, ybal

    def execute_classification(self, folder, labels_file):
        print("Loading Patients graphs...")
        data = self.load_dataset_networks_grakel(folder, labels_file)
        
        if(not os.path.isdir( f"{self.fout}/feature_matrices")):
            os.system( f"mkdir {self.fout}/feature_matrices" )

        X=data[0]
        y=data[1]
        samples = data[2]

        #kernels = ['VertexHistogram', 'EdgeHistogram', 'WeisfeilerLehman', 'PropagationAttr','RandomWalk']
        kernels = ['VertexHistogram', 'EdgeHistogram', 'WeisfeilerLehman']
        results={}
        for kn in kernels:
            print('- Kernel matrix for', kn)
            gk = eval( f'gkr.{kn}(normalize=True)' )
            
            print('\tCalculating feature matrix ', kn)
            matrix = gk.fit_transform(X)
            indnan = np.isnan(matrix)
            matrix[indnan] = 0
            f=open( f"{self.fout}/feature_matrices/matrix_{kn}.tsv","w")
            head = [ f"f{i}" for i in range(1, len(matrix[0]) ) ]
            f.write( '\t'.join(head)+'\n' )
            for l in matrix:
                svals = [ str(v) for v in l]
                f.write( '\t'.join(svals)+'\n' )
            f.close()
            
            print('\tExecuting tunning ', kn)
            grid = { 'C': [0.01, 0.1, 1, 10, 100], 'gamma': [1,0.1,0.01,0.001] }
            clf = SVC(kernel='precomputed')
            #grid_search = KernelGridSearchCV( clf, param_grid = grid, cv = 5,  random_state = 0 )
            
            #matrix, y = self.balance_dataset(matrix, y)
            
            grid_search = KernelGridSearchCV( clf, param_grid = grid, cv = 5,  random_state = 42 )
            grid_search.fit(matrix, y)
            params = grid_search.best_params_
            clf = grid_search.best_estimator_
            results[kn] = grid_search.complete_score
            #results[kn]['accuracy'] = grid_search.best_score_
            results[kn]['C'] = params['C']
            results[kn]['gamma'] = params['gamma']
            
            """
            G_train, G_test, y_train, y_test = train_test_split(Xbal, ybal, test_size=0.1, random_state=42)
            parameters = {'kernel':('linear', 'rbf'), 'C': 10. ** np.arange(-2,3) }
            clfg = GridSearchCV(clf, parameters)
            clfg.fit(G_train, y_train)
            y_pred = clfg.predict(G_test)
            acc = accuracy_score(y_test, y_pred)
            """
            
            #grid_search.fit( Xbal, ybal)
            """
            ncv = 3
            if( len(ybal) > 2000 and len(ybal) < 20000 ):
                ncv = 5
            if( len(ybal) > 20000 ):
                ncv = 10
            cv = StratifiedKFold( n_splits = ncv, shuffle = True, random_state = 0 )
            best_ac = -1
            best_clf = None
            scores = []
            i=1
            for train, test in cv.split(np.zeros(len(ybal)), ybal):
                X_train = Xbal[train][:, train]
                #print(y)
                y_train = ybal[train]
                X_test  = Xbal[test][:, train]
                y_test  = ybal[test]

                clf.fit(X_train, y_train)

                # The class only supports the accuracy score for now.
                ac = accuracy_score(y_test, clf.predict(X_test))
                scores.append(ac)
                print('ac', 'round {i}', ac )
                if(ac > best_ac):
                    best_ac = ac
                    best_clf = clf
                    
                i+=1

            score = np.mean(scores)
            print('mean', score)
            #clf = grid_search.best_estimator_
            #results[kn] = grid_search.best_score_ * 100.0
            results[kn] = best_ac
            """
            #results[kn] = acc
            print('acc', results[kn] )

        return samples, results

# Step 3
class Export_results:
    def __init__(self, fout):
        self.fout = fout
    
    def _compute_similarity_in_samples(self, nmsamples, folder, labels_file):
        y_real={}
        with open(folder+'/'+labels_file) as f:
            for line in f:
                l=line.replace("\n","").split("\t")
                y_real[l[0]]=float(l[1])
                
        union = []
        for h in ['all_similar','same_class_similar', 'opposite_class_similar']:
            for s in ['_qty', '_samples', '_values']:
                union.append(h+s)
        
        gf = open( f"{self.fout}/comparison_classes_graphSim.tsv", "w")
        gf.write("kernel\ttarget_sample\t%s\n" %( '\t'.join(union) ) )
        for f in os.listdir( f"{self.fout}/feature_matrices" ):
            kn = f.split('.')[0].split('_')[1]
            
            i=-1
            g = open( f"{self.fout}/feature_matrices/{f}", 'r')
            for line in g:
                if(i>=0):
                    ls = line.replace('\n','').split('\t')
                    s1 = nmsamples[i]
                    class1 = y_real[s1]
                    
                    j = 0
                    samples = { 'all': [], 'same_class': [], 'opposite_class': [] }
                    values = { 'all': [], 'same_class': [], 'opposite_class': [] }
                    for v in ls:
                        if(i!=j):
                            s2 = nmsamples[j]
                            class2 = y_real[s2]
                            
                            v = float(v)
                            if( v > 0.7 ):
                                v = str(v)
                                samples['all'].append(s2)
                                values['all'].append(v)
                                if( class1==class2 ):
                                    samples['same_class'].append(s2)
                                    values['same_class'].append(v)
                                else:
                                    samples['opposite_class'].append(s2)
                                    values['opposite_class'].append(v)
                        j+=1
                    
                    line = [kn, s1]
                    for k in samples.keys():
                        vals = [ str(x) for x in values[k] ]
                        line += [ str( len(samples[k]) ), (','.join(samples[k])), (','.join( vals )) ]
                    line = '\t'.join(line)
                    gf.write(line+'\n')
                    
                i+=1
            g.close()
        gf.close()
    
    def make_report(self, folder, samples, results, labels_file):
        print("Exporting report...")
        
        self._compute_similarity_in_samples(samples, folder, labels_file)
        
        kernels = list(results.keys())
        keys = '\t'.join( list( results[ kernels[0] ].keys() ) )
        
        f=open( f"{self.fout}/outcome_prediction_report.tsv","w")
        f.write( f"kernel\t{keys}\n")
        for k in kernels:
            vals = [ str(x) for x in results[k].values() ]
            vals = '\t'.join(vals)
            f.write( "%s\t%s\n" %(k, vals) )
        f.close()
        

class Pipeline:
    def _profile(self, folder, obj, description, call, mask_expr_table_file, raw_expr_table_file, labels_file, cutoff_up, cutoff_down, nodes_enrich, edges_enrich, identifier, samples, results):
        
        tbefore = time.time()
        membefore = process_memory()
        
        result = eval( f"obj.{call}" )
        
        diffMem = process_memory() - membefore
        diffTime = time.time() - tbefore
        timeExec = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
        
        with open( f'{folder}/log_execution_details.tsv', 'a') as f:
            f.write( f"NetworkBased\t{description}\t{timeExec}\t{diffTime}\t{diffMem}\n" )
        
        return result
        
    def run(self, folder, mask_expr_table_file, raw_expr_table_file, labels_file, cutoff_up, cutoff_down, nodes_enrich, edges_enrich, identifier):
        samples = None 
        results = None
        
        step1 = Generate_personalised_networks(cutoff_up, cutoff_down)
        self._profile(folder, step1, 'Step 1 - Personalized sample networks generation', 'make_patients_graph(folder, mask_expr_table_file, raw_expr_table_file )', mask_expr_table_file, raw_expr_table_file, labels_file, cutoff_up, cutoff_down, nodes_enrich, edges_enrich, identifier, samples, results)
        # step1.make_patients_graph(folder, mask_expr_table_file, raw_expr_table_file )
        
        for ne in nodes_enrich:
            for ee in edges_enrich:
                print( f'--- Enrichment {ne} and {ee}')
                step2 = KernelMatrixBasedPrediction(folder, ne, ee, identifier)
                fout = step2.fout
                samples, results = self._profile(folder, step2, f'Step 2 - Execution classification: enrichment nodes with {ne} and edges with {ee}', 'execute_classification(folder, labels_file)', mask_expr_table_file, raw_expr_table_file, labels_file, cutoff_up, cutoff_down, nodes_enrich, edges_enrich, identifier, samples, results)
                # samples, results=step2.execute_classification(folder, labels_file)

                step3 = Export_results(fout)
                self._profile(folder, step3, f'Step 3 - Exporting results: enrichment nodes with {ne} and edges with {ee}', 'make_report( folder, samples, results, labels_file)', mask_expr_table_file, raw_expr_table_file, labels_file, cutoff_up, cutoff_down, nodes_enrich, edges_enrich, identifier, samples, results)
                # step3.make_report( folder, samples, results, labels_file)
                
