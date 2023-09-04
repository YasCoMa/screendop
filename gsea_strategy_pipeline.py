import os
import pandas as pd
import statistics as st
import numpy as np
import gseapy as gp

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
from imblearn.over_sampling import SMOTE

workflow_path = os.environ.get('path_workflow_screendop')
if(workflow_path==None):
    raise Exception('You forgot to setup the path_workflow environment variable with the path to the workflow')
    
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

class Gsea_strategy:
    def __init__(self, folder, flag_balance):
        self.folder = folder
        self.flag_balance = flag_balance
        ext=''
        if( flag_balance ):
            ext='_with_balancing'
        
        self.fout = f'{folder}/gsea_strategy{ext}'
        if( not os.path.isdir(self.fout) ):
            os.system( f"mkdir {self.fout}" )

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
    
    # step 1      
    def generate_features_gsea(self, raw_expr_table_file, patients_labels_file):
        folder = self.folder
        ffout = self.fout
        
        mp={}
        df = pd.read_csv( f"{workflow_path}/mapping_hgnc_uniprot.txt", sep='\t', header=None)
        for i in df.index:
            mp[ df.iloc[i,1] ] = df.iloc[i,0]
        genesokuni = set(mp.keys())

        gmt_lib='KEGG_2021_HUMAN'
        gmt = gp.parser.download_library(gmt_lib, 'Human')
        genesok = set()
        for k in gmt.keys():
            genesok = genesok.union( set(gmt[k]) )
        
        lb = pd.read_csv( f'{folder}/{patients_labels_file}', sep='\t', header=None)
        lb.columns = ['sample','class']
        lb['sample'] = lb['sample'].apply( lambda x: str(x) )
        
        df = pd.read_csv( f'{folder}/{raw_expr_table_file}', sep='\t')
        df = df[ df['id'].isin(genesokuni) ]

        cols = df.columns.tolist()
        cols[0] = 'Gene'
        df['id'] = df['id'].apply(lambda x: mp[x])
        df = df[ df['id'].isin(genesok) ]
        #df['id'] = [ mp[v] for v in df['id'].values ]
        df.columns = cols
                    
        df = df[ df['Gene'].isin(genesok) ]
        ss = gp.ssgsea(data = df, gene_sets=gmt, outdir=None, sample_norm_method='rank', no_plot=True)
        rd = ss.res2d
        
        rd['Name'] = [ str(v) for v in rd['Name'] ]
        samples = rd['Name'].unique()
        paths = rd['Term'].unique()
        f = open( f'{ffout}/gsea_features.tsv', 'w')
        f.write( '\t'.join( ['sample', 'y']+list(paths) )+'\n' )
        for s in samples:
            y = lb[ lb['sample'] == s ]['class'].values[0]
            es = rd[ rd['Name'] == s ]['ES'].values
            vals = [ str(v) for v in es ]
            f.write( '\t'.join( [s, str(y)]+list(vals) )+'\n' )
        f.close()
    
    # step 2
    def train_evaluate(self):
        ffout = self.fout
        
        df = pd.read_csv( f'{ffout}/gsea_features.tsv', sep='\t')
        y = list( df['y'].values )
        X = np.array( df.iloc[ :, 2: ].values )
        print( 'Features (n):', len(X[0]) )
        
        if( self.flag_balance ):
            oversample = SMOTE()
            X, y = oversample.fit_resample(X, y)
        
        model = AdaBoostClassifier()
        cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
        
        # define the grid of values to search
        grid = dict()
        grid['n_estimators'] = [10, 50, 100, 500]
        grid['learning_rate'] = [0.0001, 0.001, 0.01, 0.1, 1.0]
        grid_search = GridSearchCV( estimator=model, param_grid=grid, n_jobs=-1, cv=cv, scoring='accuracy')
        grid_result = grid_search.fit(X, y)
        
        lr = 0
        ne = 0
        ant = 0
        means = grid_result.cv_results_['mean_test_score']
        stds = grid_result.cv_results_['std_test_score']
        params = grid_result.cv_results_['params']
        for mean, stdev, param in zip(means, stds, params):
            if( mean > ant ):
                ant=mean
                lr = param['learning_rate']
                ne = param['n_estimators']
        
        # get best parameters from grid search and fit one model
        clf = AdaBoostClassifier( learning_rate=lr, n_estimators=ne, random_state=0)
        metrics=['accuracy','precision', 'recall','f1']
        
        list_metrics = {}
        for m in metrics:
            crs = cross_val_score(clf, X, y, cv=10, scoring=m)
            list_metrics[m] = crs
        
        return list_metrics, grid_result

    # step3
    def _report_performance(self, list_metrics, grid_result):
        ffout = self.fout
        
        means = grid_result.cv_results_['mean_test_score']
        stds = grid_result.cv_results_['std_test_score']
        params = grid_result.cv_results_['params']
        
        f = open( f'{ffout}/grid_search_parameters.tsv', 'w')
        f.write( 'mean_accuracy\tstd_accuracy\tlearning_rate\tn_estimator\n' )
        for mean, stdev, param in zip(means, stds, params):
            f.write( '%.6f\t%.6f\t%.6f\t%.6f\n' %(mean, stdev, param['learning_rate'], param['n_estimators']) )
        f.close()
        
        metrics=['accuracy','precision', 'recall','f1']
        f = open( f'{ffout}/cross_val_evaluation.tsv', 'w')
        f.write( 'metric\tmean\tstdev\tvalues\n' )
        for m in metrics:
            crs = list_metrics[m]
            vm = [ str(v) for v in crs ]
            vm = ','.join(vm)
            mean = st.mean(crs)
            stdev = st.stdev(crs)
            f.write( '%s\t%.6f\t%.6f\t%s\n' %(m, mean, stdev, vm) )
        f.close()  
    
    def report_performance_sim_Matrix_ppi(self, mask_expr_table_file, raw_expr_table_file, patients_labels_file, list_metrics, grid_results):
        folder = self.folder
        fout = self.fout
        
        self._report_performance(list_metrics, grid_results)
        
        # Check and discuss the impact of both diseases in the ppi network
        # check the intersection of highly correlated genes with protein interactions in the disease, non disease and all samples
        ppi = pd.read_csv( f"{workflow_path}/interactome_hint.tsv", sep='\t', header=None)
        ppi.columns=['p1','p2','class']
        ppis = {}
        for i in ppi.index:
            k1 = ppi.loc[i, 'p1']
            k2 = ppi.loc[i, 'p2']
            
            if(not k1 in ppis):
                ppis[k1]=set()
            ppis[ k1 ].add( k2 )
            
            if(not k2 in ppis):
                ppis[k2]=set()
            ppis[ k2 ].add( k1 )
        
        conds = ['positive', 'negative']
        values=[ '> 0.7', '< -0.7']
        
        f = open( f'{fout}/report_overlapping_corr_ppis.tsv', 'w')
        f.write("sampleset\ttype_correlation\tprotein\tqty_neighbors_corr\tqty_neighbors_ppi\tqty_neighbors_intersection\tneighbors_corr\tneighbors_corr_values\n")
        f.close()
        
        lb = pd.read_csv( f'{folder}/{patients_labels_file}', sep='\t', header=None)
        lb.columns = ['sample','class']
        samplesets = { 'dead': list( lb[ lb['class']==1 ]['sample'].values ), 'alive': list( lb[ lb['class']==0 ]['sample'].values ), 'all': list(lb['sample'].values) }
        for g in samplesets:
            mask = pd.read_csv( f'{folder}/{mask_expr_table_file}', sep='\t')
            genesok = mask['id'].values # significant dges
            
            cols = [str(v) for v in samplesets[g]]
            df = pd.read_csv( f'{folder}/{raw_expr_table_file}', sep='\t', index_col=0)
            #tr = df[ df.index.isin(genesok) ][ cols ].transpose().corr()
            tr = df[ cols ].transpose().corr()
            n = len(tr.index)
            
            for j in [0,1]:
                cutoff = values[j]
                clas = conds[j]
                
                cnt={}
                corr={}
                for c in tr.columns:
                    cnt[c] = list( tr[ eval(f'tr[c] {cutoff}') ].index )
                    corr[c] = [ str(v) for v in list( tr[ eval(f'tr[c] {cutoff}') ][c] ) ]
                    aux1=[]
                    aux2=[]
                    for a,b in zip(cnt[c], corr[c]):
                        if(a!=c):  
                            aux1.append(a)
                            aux2.append(b)
                    cnt[c] = aux1
                    corr[c] = aux2        
                    
                    nn = ','.join( cnt[c] )
                    nv = ','.join( corr[c] )
                    intersec = []
                    ni = -1
                    np = -1
                    if(c in ppis):
                        np = len(ppis[c])
                        intersec = list( ppis[c].intersection( set(cnt[c]) ) )
                        ni=len(intersec)
                    intersec = ','.join(intersec)    
                    if( len(cnt[c]) > 0):
                        with open( f'{fout}/report_overlapping_corr_ppis.tsv', 'a') as gf:
                             gf.write( f"{g}\t{clas}\t{c}\t{ len(cnt[c]) }\t{ np }\t{ ni }\t{ nn }\t{ nv }\n" )
                 
                b = list( cnt.items() )
                c = list( filter(lambda x: len(x[1])>0, b))
                print( f"{g} -- {clas} -- Number of proteins that have more than one neighbor {clas}ly correlated: ", len(c))
    
class Pipeline:
    def run(self, folder, mask_expr_table_file, raw_expr_table_file, patient_labels_file, flag_balance):
        p = Gsea_strategy(folder, flag_balance)
        p.generate_features_gsea(raw_expr_table_file, patient_labels_file)
        list_metrics, grid_results = p.train_evaluate()
        p.report_performance_sim_Matrix_ppi(mask_expr_table_file, raw_expr_table_file, patient_labels_file, list_metrics, grid_results)

