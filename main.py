
import os
import json
import sys
from collections.abc import Iterable

workflow_path = os.environ.get('path_workflow_screendop')
if(workflow_path==None):
    raise Exception('You forgot to setup the path_workflow environment variable with the path to the workflow')
    
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

sys.path.insert(0, workflow_path)

from vanilla_kernelFunc_pipeline import Pipeline as PipeStrategyA
from gsea_strategy_pipeline import Pipeline as PipeStrategyB

class Pipeline_ScreenDOP:
    
    def _validate_input(self, config):
        flag = False
        if( os.path.isfile(config) ):
            try:
                with open(config, 'r') as g:
                    cfg = json.load(g)
                    
                aux = True
                if('experiments' in cfg):
                    for e in cfg['experiments']:
                        if( ('folder' in e) and ('mask_expression_table' in e) and ('raw_expression_table' in e) ):
                            if ( (e['folder']!='' and e['folder']!=None) and (e['mask_expression_table']!='' and e['mask_expression_table']!=None) and (e['raw_expression_table']!='' and e['raw_expression_table']!=None) ):
                                if( os.path.isdir(e['folder']) ):
                                    folder = e['folder']
                                    
                                    aux = aux and self._validate_file(e, 'Labels file', 'labels_file')
                                    aux = aux and self._validate_file(e, 'Raw Expression file', 'raw_expression_table')
                                    aux = aux and self._validate_file(e, 'Mask Expression file', 'mask_expression_table')
                                    
                                    aux = aux and True
                                else:
                                    aux = False
                                    print(f'Error - {ide}: The directory was not found in {folder}')
                            else:
                                aux = False
                                print("Error - Mandatory fields are empty")
                        else:
                            aux = False
                            print("Error - Configuration file does not have the mandatory fields")
                else:
                    aux = False
                    print("Error - Configuration file in a bad format")
            except:
                aux = False
                print("Error - Configuration file cannot not be loaded and parsed")
        else:
            print("Error - Configuration file not found")
        
        flag = aux
           
        return flag
    
    def _validate_file(self, e, name, ide):
        flag = False
        if( ide in e ):
            if ( (e[ ide ]!='' and e[ ide ]!=None) ):
                folder = e['folder']
                if(folder[-1]=='/'):
                    folder = folder[:-1]
                if( os.path.isfile(folder+'/'+e[ ide ]) ):
                    flag = True
                else:
                    print (f'Error: {name} was not found in {folder}')
            else:
                print( f"Error - {name} field is empty")
        else:
            print( f"Error - {name} field is missing from configuration")
            
        return flag
    
    def check_cutoffs(self, e):
        cuts = {}
        
        i=0
        vals = [1, -1]
        for op in ['up', 'down']:
            cuts[op] = vals[i]
            
            ide = f'deg_cutoff_{op}'
            name = f'Cutoff {op} regulated'
            if( ide in e ):
                if ( (e[ ide ]!='' and e[ ide ]!=None) ):
                    cuts[op] = e[ide]
                else:
                    print( f"Information - {name} field is empty - switching to default value")
            else:
                print( f"Information - {name} field is missing from configuration - switching to default value")
            i+=1
            
        return cuts
    
    def check_flagBalance(self, e):
        flag = False
        ide = 'flag_balance'
        name = 'Balancing dataset flag'
        if( ide in e ):
            if ( (e[ ide ]!='' and e[ ide ]!=None) ):
                flag = e[ ide ]
            else:
                print( f"Information - {name} field is empty - switching to default value")
        else:
            print( f"Information - {name} field is missing from configuration - switching to default value")
            i+=1
            
        return flag
    
    def check_net_enrichment(self, e):
        net_ops = {}
        
        options = ["label", "weight", "all"]
        for op in ['nodes', 'edges']:
            net_ops[op]=['all']
            ide = f'{op}_enrichment'
            name = f'Enrichment {op}'
            if( ide in e ):
                if ( (e[ ide ]!='' and e[ ide ]!=None) ):
                    flag=True
                    if( isinstance(e[ide], Iterable) ):
                        e[ide] = set(e[ide])
                    elif( isinstance(e[ide], str) ):
                        e[ide] = set( [ e[ide] ] )
                    else:
                        flag=False
                        
                    if( flag ):
                        intersec = e[ide].intersection( set(options) )
                        if( len(intersec) == len(e[ide]) ):
                            net_ops[op] = list( e[ide] )
                        else:
                            print( f"Information - {name} invalid value - switching to default value")
                    else:
                        print( f"Information - {name} invalid value type - switching to default value")
                else:
                    print( f"Information - {name} field is empty - switching to default value")
            else:
                print( f"Information - {name} field is missing from configuration - switching to default value")
        
        return net_ops
    
    def run(self, option, config):
        if( self._validate_input(config) ):
            if( option in range(4) ):
                with open(config, 'r') as g:
                    cfg = json.load(g)
                
                i=1
                for e in cfg['experiments']:
                    ide = f'exp-{i}'
                    if( 'identifier' in e ):
                        if( e['identifier']!='' and e['identifier']!=None ) :
                            ide = e['identifier']
                            
                    print(f'Running Experiment [{ide}]')
                    
                    folder = e['folder']
                    if(folder[-1]=='/'):
                        folder = folder[:-1]
                        
                    if( not os.path.isfile(f'{folder}/log_execution_details.tsv') ):
                        with open( f'{folder}/log_execution_details.tsv', 'w') as f:
                            f.write( f"strategy\tdescription\tmoment\texecution_time\tmemory_usage\n" )
                    
                    flagop = ( option==0 ) 
                    if( flagop or option==1): # vanilla implementation
                        print('\tExecuting Vanilla implementation')
                        
                        cuts = self.check_cutoffs(e)
                        up = cuts['up']
                        down = cuts['down']
                        
                        net_enrich = self.check_net_enrichment(e)
                        enodes = net_enrich['nodes']
                        eedges = net_enrich['edges']
                        
                        a = PipeStrategyA(  )
                        a.run( folder, e['mask_expression_table'], e['raw_expression_table'], e['labels_file'], up, down, enodes, eedges, ide )
                        
                    if( flagop or option==2): # enriched network implementation
                        print('\tExecuting Enriched network implementation')
                        
                        flagBal = self.check_flagBalance(e)
                        
                        a = PipeStrategyB(  )
                        a.run( folder, e['mask_expression_table'], e['raw_expression_table'], e['labels_file'], flagBal )
                    i+=1
            else:
                print('Error - Invalid option')
            
import sys

#op = int(sys.argv[1])
#config = sys.argv[2]

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' ScreenDOP - Screening of strategies for disease outcome prediction', formatter_class=RawTextHelpFormatter)

parser.add_argument("-cf", "--configuration_file", action="store", help="(For both modes) Folder to store the files (use the folder where the required files can be found, ex.: /home/user/experiment/ )\n")

parser.add_argument("-rt", "--running_step", action="store", help="0 - Run all steps\n\
1 - Run strategy 1: Graph kernel implementation\n\
1 - Run strategy 2: Enriched network implementation\n\
", type=int)

args = parser.parse_args()

a = Pipeline_ScreenDOP()
a.run( args.running_step, args.configuration_file)
