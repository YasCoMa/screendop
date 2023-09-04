
# Application on outcome disease prediction on leukaemia dataset
# http://psb.stanford.edu/psb-online/proceedings/psb07/borgwardt.pdf
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE425

import os
import statistics as st
import urllib.request
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt 
import json
from rnanorm import UQ

class Outcome_experiment_leukaemia:
    def get_gpl_mapping(self):
        for i in range(317,320):
            f=open("data_preparation/leukaemia/map_gpl_"+str(i)+".tsv","w")
            f.close()

            link="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL"+str(i)+"&db=GeoDb_blob03"
            f = urllib.request.urlopen(link)
            file = f.read()
            a=str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace("\\t","\t").replace('>"','>')
            f=open("temp.txt","w")
            f.write(a)
            f.close()

            ok=False
            f=open("temp.txt","r")
            for line in f:
                if(ok):
                    l=line.replace("\n","").split("\t")
                    if(len(l)>3):
                        with open("data_preparation/leukaemia/map_gpl_"+str(i)+".tsv","a") as g:
                            g.write("%s\t%s\n" %(l[0], l[2]))
                if(line.find("<strong>ID")!=-1):
                    ok=True
            f.close()

    def build_dictionary(self):
        gpl={}
        f=open("data_preparation/leukaemia/gpl_samples.csv","r")
        for line in f:
            l=line.replace("\n","").split(",")
            gpl[l[0]]=l[1:]
        f.close()

        mapping_gpl={}
        for i in range(317,320):
            temp={}
            f=open("data_preparation/leukaemia/map_gpl_"+str(i)+".tsv","r")
            for line in f:
                l=line.replace("\n","").split("\t")
                temp[l[0]]=l[1]
            f.close()
            mapping_gpl[str(i)] = temp
        return [gpl, mapping_gpl]

    def getting_data(self):
        mapinfo=self.build_dictionary()
        gpl=mapinfo[0]
        mapping=mapinfo[1]

        header={}
        findgpl={}
        for k in gpl.keys():
            header[k]=[]
            findgpl[k]=0
            f=open("data_preparation/leukaemia/gpl-"+k+"_samples_laukaemia.tsv","w")
            f.close()

        samples={}
        c=0
        for i in range(6228,6347):
            key=""
            for k in gpl.keys():
                if("GSM"+str(i) in gpl[k]):
                    findgpl[k]+=1
                    key=k
            #print(key)
            values=[]
            link="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM"+str(i)+"&db=GeoDb_blob02"
            f = urllib.request.urlopen(link)
            file = f.read()
            a=str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace("\\t","\t").replace('>"','>')
            f=open("temp.txt","w")
            f.write(a)
            f.close()

            ok=False
            f=open("temp.txt","r")
            for line in f:
                if(ok):
                    l=line.replace("\n","").split("\t")
                    if(len(l)>3 ):
                        #print(l[0],mapping[key][l[0]])
                        if(mapping[key][l[0]]!=""):
                            values.append(l[-1])
                            if(findgpl[key]==1):
                                header[key].append(mapping[key][l[0]])

                if(line.find("<strong>ID_REF")!=-1):
                    ok=True
            f.close()

            samples[str(i)+"-"+key]=values
            c+=1

        gplid="317"
        for h in header.keys():
            text=["gene"]
            for m in samples.keys():
                info=m.split("-")
                gplid=info[1]
                if(gplid==h):
                    text.append(info[0])
            #print(text)
            with open("data_preparation/leukaemia/gpl-"+h+"_samples_laukaemia.tsv","a") as g:
                g.write("\t".join(text)+"\n")
            
            for j in range(len(header[h])):
                text=[header[h][j]]
                for m in samples.keys():
                    info=m.split("-")
                    gplid=info[1]
                    if(gplid==h):
                        text.append(samples[m][j])

                #print(text)
                with open("data_preparation/leukaemia/gpl-"+h+"_samples_laukaemia.tsv","a") as g:
                    g.write("\t".join(text)+"\n") 

    def do_mapping_hgnc_uniprot(self):
        mapping={}
        f=open("data_preparation/leukaemia/mapping_hgnc_uniprot.txt", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            mapping[l[0]]=l[1]
        f.close()

        return mapping

    def do_mapping_samples(self):
        mapping={}
        f=open("data_preparation/leukaemia/names_samples.tsv", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            mapping[l[1]]=l[0].replace("GSM","")
        f.close()

        return mapping

    def prepare_expression_table_and_labels(self):
        print("Loading mapping between hgnc and uniprot")
        mapping=self.do_mapping_hgnc_uniprot()

        print("Loading expression data by sample")
        expr_data={}
        for i in range(317,320):
            samples=[]
            c=0
            f=open("data_preparation/leukaemia/gpl-"+str(i)+"_samples_laukaemia.tsv", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(c==0):
                    samples=l[1:]
                    for k in samples:
                        if(not k in expr_data.keys()):
                            expr_data[k]={}
                else:
                    j=1
                    if(l[0] in mapping.keys()):
                        for k in samples:
                            value=0.0
                            if(l[j]!=""):
                                value=float(l[j])
                            if(not l[0] in expr_data[k].keys()):
                                expr_data[k][mapping[l[0]]]=[]
                            expr_data[k][mapping[l[0]]].append(value)

                            j+=1
                c+=1
            f.close()

        print("Loading mapping between accesion number for sample and the name given in the experimental design")
        mapping_samples=self.do_mapping_samples()

        print("Loading outcome of the disease")
        y_real={}
        f=open("data_preparation/leukaemia/y_real_deadOrAlive.tsv", "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            outcome=0
            if(l[-3]=="dead"):
                outcome=1
            y_real[mapping_samples[l[0]]]=outcome
        f.close()
        
        print("Generating normalized raw expression table...")
        f=open("data_preparation/leukaemia/patients_labels_leukaemia.tsv","w")
        f.close()

        f=open("data_preparation/leukaemia/raw_expression_table.tsv","w")
        #f.write("\t")
        f.close()
        ysamples = list(y_real.keys())
        idsok = list( filter( lambda x: (x in ysamples), list(expr_data.keys())) )
        
        first=""
        for k in idsok:
            with open("data_preparation/leukaemia/patients_labels_leukaemia.tsv","a") as gf:
                gf.write("%s\t%i\n" %(k,y_real[k]))

        with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
            gf.write("id\t"+('\t'.join(idsok))+'\n')

        for gene in expr_data[ idsok[0] ].keys():
            i=0
            dat=[]
            for k in idsok:
                v = 0
                if( gene in expr_data[k].keys() ):
                    v = st.mean( expr_data[k][gene] )
                    
                if(i==0):
                    dat.append(gene)
                    #with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
                    #    gf.write("%s" %( gene ) )
                dat.append( str(v) )
                #with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
                #    gf.write("\t%.4f" %( v ) )
                i+=1

            with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
                gf.write( ( '\t'.join(dat) )+"\n" )

        print("Generating mask expression table...")
        df = pd.read_csv('data_preparation/leukaemia/differential_gene_expression.csv', sep=',', index_col=0)
        df = df[ df['P.Value']<=0.05 ]
        mp={}
        for i in df.index:
            mp[i] = df.loc[i, 'logFC']

        with open("data_preparation/leukaemia/mask_expression_table.tsv","w") as gf:
            gf.write("id\t"+('\t'.join(idsok))+'\n')

        f=open("data_preparation/leukaemia/mask_expression_table.tsv","a")
        f.close()
        for gene in expr_data[ idsok[0] ].keys():
            i=0
            dat=[]
            for k in idsok:
                v = 0
                if( gene in mp.keys() ):
                    v = mp[gene]
                    
                if(i==0):
                    dat.append(gene)
                    #with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
                    #    gf.write("%s" %( gene ) )
                dat.append( str(v) )
                #with open("data_preparation/leukaemia/raw_expression_table.tsv","a") as gf:
                #    gf.write("\t%.4f" %( v ) )
                i+=1

            with open("data_preparation/leukaemia/mask_expression_table.tsv","a") as gf:
                gf.write( ( '\t'.join(dat) )+"\n" )
                
        print("Data available in data_preparation/leukaemia/")

# do the networks for each patient using 0 as cutoff
# export to graphml and use the graph kernel to separate the patients

class Outcome_experiment_ovarian:

    def mapping_ilumina_uniprot(self):
        mapping={}
        f=open("data_preparation/ovarian/uniprot_to_illumina.txt","r")
        i=0
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[0]!="" and l[1]!=""):
                mapping[l[0]]=l[1]
            i+=1 
        f.close()
        return mapping

    def prepare_files_for_deseq(self):
        f=open("data_preparation/ovarian/GSE140082_series_matrix.txt","r")
        for line in f:
            l=line.replace("\n","").replace('\"',"")
            #if(l.find("Sample_geo_accession")!=-1):
            #    samples_names=l.split("\t")[1:]
            if(l.find("!Series_sample_id")!=-1 ):
                samples=l.split("\t")[1].split(" ")[:-1]
                #print(len(samples), samples)
            
            if(l.find("!Sample_characteristics_ch1")!=-1 and l.find("treatment:")!=-1):
                conditions=l.replace("treatment: ","").split("\t")[1:]
                #print(len(conditions), conditions)

            if(l.find("!Sample_characteristics_ch1")!=-1 and l.find("t1_cluster_name:")!=-1):
                moltypes=l.replace("t1_cluster_name: ","").split("\t")[1:]
                #print(len(moltypes), moltypes)
        f.close()

        subconditions=["immunoreactive","mesenchymal", "proliferative", "differentiated"]
        i=0
        map_={}
        for s in samples:
            map_[s]=moltypes[i]
            i+=1

        for s in subconditions:
            print(s)

            subsamples=[]
            for s2 in samples:
                if map_[s2]==s:
                    subsamples.append(s2)

            f=open("data_preparation/ovarian/"+s+"_count_matrix.txt","w")
            f.write("id,"+(",".join(subsamples))+"\n")
            f.close()
            print(len(subsamples))

            co=0
            f=open("data_preparation/ovarian/GSE140082_geo.rawdata.csv","r") 
            for line in f:
                if(co!=0):
                    l=line.replace("\n","").split(",")
                    txt=[l[0]]
                    c=1
                    a=0
                    for value in l:
                        if(c%2==0):
                            if (map_[samples[a]]==s):
                                if(int(float(value))<0):
                                    value=str(-1*int(float(value)))
                                    value=str(0)
                                else:
                                    value=str(int(float(value)))

                                txt.append(value)
                            a+=1
                        c+=1
                    
                    with open(s+"_count_matrix.txt","a") as gf:
                        gf.write(",".join(txt)+"\n")
                co+=1
            f.close()
            print(len(txt))

            f=open("data_preparation/ovarian/"+s+"_design_matrix.txt","w")
            f.write(",treatment\n")
            i=0
            for c in conditions:
                if(map_[samples[i]]==s):
                    f.write(samples[i]+","+c+"_"+s+"\n")
                i+=1
            f.close()

            # template deseq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
            f=open("data_preparation/ovarian/"+s+"_rscript.R","w")
            f.write("library(DESeq2) \n")
            f.write("data=as.matrix(read.csv('data_preparation/ovarian/"+s+"_count_matrix.txt', sep=',', row.names='id', header=TRUE)) \n")
            f.write("design=read.csv('data_preparation/ovarian/"+s+"_design_matrix.txt', sep=',', row.names=1, header=TRUE) \n")
            f.write("data[is.na(data)] <- 0\n")
            f.write("dataset <- DESeqDataSetFromMatrix(countData = data, colData = design, design = ~treatment) \n")
            f.write("dds <- DESeq(dataset) \n")
            f.write("res <- results(dds) \n")
            f.write("data_ = as.data.frame(res) \n")
            f.write("data_ = data_[!is.na(data_$log2FoldChange),] \n")
            f.write("data_ = data_[!is.na(data_$padj),] \n")
            #f.write("data_ = data_[data_$padj<=0.05,] \n")
            f.write("write.table(data_, file='data_preparation/ovarian/"+s+"_result.csv', sep=',', quote=FALSE) \n")
            f.close()

            os.system("Rscript data_preparation/ovarian/"+s+"_rscript.R")

    def _save_raw_expression_table_normalized(self, subconditions, mapping):
        df = pd.DataFrame()
        genes = []
        allcols = []
        for s in subconditions:
            dff = pd.read_csv("data_preparation/ovarian/"+s+"_count_matrix.txt", index_col=0)
            passed = set()
            genes = []
            i=0
            inds=[]
            for y in dff.index:
                if ( y in mapping.keys() ):
                    passed.add(y)
                    if(not mapping[y] in genes):
                        genes.append(mapping[y])
                        inds.append(i)
                i+=1
            cols = list(dff.columns)
            dff = dff.iloc[ inds, : ]
            #dff.index = list(set(genes))
            allcols += cols
            df = pd.concat([df, dff[ cols ] ], axis=1, ignore_index=True)
    
            #for c in dff.columns:
            #    df[c] = dff[c]
        df.index = genes
        df.columns = allcols
        df = df.transpose()
        """
        head = df.columns
        index = df.index
        df = pd.DataFrame(df.values[1:,:], columns=head)
        df.index = index
        df = df.fillna(0)
        print( len(df), len(df.columns), df.head())
           """     
        transformer = UQ()
        transformer.fit(df)
        index = df.index
        cols = df.columns
        df = transformer.transform(df)
        
        df = pd.DataFrame(df, columns = cols, index = index)
        df = df.transpose()
        df.reset_index(inplace=True)
        cols = df.columns.tolist()
        cols[0]='id'
        df.columns=cols
        df = df[ ['id']+cols[1:] ]
        df.to_csv('data_preparation/ovarian/raw_expression_table.tsv', sep='\t', index=None)
    
    def prepare_expression_table_and_labels(self):
        mapping=self.mapping_ilumina_uniprot()

        #samples_names=[]
        y=[]
        y_real={}
        expr_data={}
        start_table=False
        samples=[]

        c=0
        f=open("data_preparation/ovarian/GSE140082_series_matrix.txt","r")
        for line in f:
            l=line.replace("\n","").replace('\"',"")
            #if(l.find("Sample_geo_accession")!=-1):
            #    samples_names=l.split("\t")[1:]
            if(l.find("Sample_characteristics_ch1")!=-1 and l.find("final_osid")!=-1):
                y=l.replace("final_osid: ","").split("\t")[1:]

            if(l.find("!Series_sample_id")!=-1 ):
                samples=l.split("\t")[1].split(" ")[:-1]
                #print(len(samples), samples)

            if(l.find("!Sample_characteristics_ch1")!=-1 and l.find("t1_cluster_name:")!=-1):
                moltypes=l.replace("t1_cluster_name: ","").split("\t")[1:]
                #print(len(moltypes), moltypes)

            if(l.find("ID_REF")!=-1):
                start_table=True

            if(start_table):
                l=l.split("\t")
                if(c==0):
                    nf=0
                    samples=l[1:]
                    for k in samples:
                        expr_data[k]={}
                        y_real[k]=int(y[nf])
                        nf+=1
                c+=1

                break
        f.close()

        subconditions=["immunoreactive","mesenchymal", "proliferative", "differentiated"]
        
        print("Generating normalized raw expression table...")
        self._save_raw_expression_table_normalized(subconditions, mapping)
            
        i=0
        map_={}
        edata={}
        for s in samples:
            map_[s]=moltypes[i]
            i+=1
            c=0

        for s in subconditions:
            edata[s]={}
            f=open("data_preparation/ovarian/"+s+"_result.csv","r")
            for line in f:
                l=line.replace("\n","").split(",")
                if(c!=0):
                    if(l[0] in mapping.keys() and float(l[5])<=0.05  ):
                        if(not mapping[l[0]] in edata[s].keys()):
                            edata[s][mapping[l[0]]]=0
                        edata[s][mapping[l[0]]]=float(l[2])
                c+=1
            f.close()

        expr_data={}
        for s in samples:
            expr_data[s]={}
            for p in edata[map_[s]].keys():
                expr_data[s][p]=edata[map_[s]][p]
        
        print( len(samples), len(expr_data[samples[0]]) )
        
        print("Generating mask expression table...")
        f=open("data_preparation/ovarian/patients_labels_ovarian.tsv","w")
        f.close()

        f=open("data_preparation/ovarian/mask_expression_table.tsv","w")
        #f.write("id")
        f.close()
        
        first=""
        dat=[]
        for k in expr_data.keys():
            if(k in y_real.keys()):
                first=k
                with open("data_preparation/ovarian/patients_labels_ovarian.tsv","a") as gf:
                    gf.write("%s\t%i\n" %(k,y_real[k]))
                dat.append(k)
                #with open("data_preparation/ovarian/mask_expression_table.tsv","a") as gf:
                #    gf.write("\t%s" %( k ))
        with open("data_preparation/ovarian/mask_expression_table.tsv","a") as gf:
            gf.write("id\t"+('\t'.join( dat ))+'\n' )
        
        for gene in expr_data[first].keys():
            i=0
            dat=[gene]
            for k in expr_data.keys():
                v = 0
                if(gene in expr_data[k].keys()):
                    v = expr_data[k][gene]
                    #if(i==0):
                    #    dat.append(gene)
                    
                        #with open("data_preparation/ovarian/mask_expression_table.tsv","a") as gf:
                        #    gf.write("%s" %( p ) )

                dat.append( str( v ) )
                    #with open("data_preparation/ovarian/mask_expression_table.tsv","a") as gf:
                    #    gf.write("\t%.4f" %( expr_data[k][p] ) )
                i+=1
            with open("data_preparation/ovarian/mask_expression_table.tsv","a") as gf:
                gf.write( ('\t'.join( dat ))+'\n' )


        print("Data available in data_preparation/ovarian/")

class Running_config:
    def run(self, args):
        if(args.running_type in [1,2]):
            if(args.running_type==1):
                a=Outcome_experiment_leukaemia()
            if(args.running_type==2):
                a=Outcome_experiment_ovarian()

            a.prepare_expression_table_and_labels()
        else:
            print("Error: invalid option")

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Pipeline to prepare the files for leukaemia and ovarian cancer diseases for later execution in the outcome disease prediction pipeline', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-rt", "--running_type", action="store", help="1 - Run with Leukaemia data\n2 - Run with ovarian cancer data", type=int)
    args = parser.parse_args()
    r=Running_config()
    r.run(args)
