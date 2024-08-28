#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:37:38 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Xcms_base(pipeline_tools.FeatureFinderPipeline):
    def __init__(self):
        self.name = NotImplemented
        self.cores = 1
        self.memory = 8
        self.timeout = '08:00:00'
        self.algorithm = NotImplemented
        self.get_params()
    
    def get_params(self):
        params = super().get_params()
        
        merge_params = {}
        self.merge_param_set = set(merge_params.keys())
        params.update(merge_params)
        
        return params

    def setup_workspace(self):
        import subprocess
        import os
        import shutil
        
        #copy tool script to working directory
        self.tool_path = os.path.split(os.path.abspath(pipeline_tools.__file__))[0]
        self.tool_path = os.path.join(os.path.split(self.tool_path)[0], 'tools')
        tool = os.path.join(self.tool_path, 'xcms_quantify_features.R')
        shutil.copy2(tool, 'pms_quantify_peptides.py')
        
        #make conda environment for running tool
        if not os.path.isdir('~/.conda/envs/xcms_env'):
            subprocess.run(' '.join(['conda create -n xcms_env', 
                                     'r=4.3',
                                     'r-optparse=1.7.5', 
                                     'r-configr=0.3.5', 
                                     'bioconductor-xcms=4.0.0',
                                     'bioconductor-msexperiment=1.4.0', 
                                     '-c bioconda', 
                                     '-c conda-forge']),
                           shell = True)
        
    def run_job(self, job):
        super().run_job(job)
        
        import os
        import subprocess
        import shutil
        from time import time
        import traceback

        import pandas as pd
        import numpy as np
        
        #set up temporary workspace
        tmpdir = str(os.getpid())
        print(tmpdir, flush = True)
        os.mkdir(tmpdir)
        mzmls = [f for f in os.listdir() if f.endswith('.mzML') and f.startswith('20210827')]
        base_names = [f[:-5] for f in mzmls]
        psms = [f for f in os.listdir() if f.endswith('_PSMs.txt')]
        os.chdir(tmpdir)
        try:
            for file in mzmls + psms:
                os.link(f'../{file}', file)
                
            #make params file
            # THIS WILL REQUIRE READING IN THE DEFAULT PARAMS FILE IN ORDER TO GET THE RIGHT DATA TYPES

            start = time()
            peptide_results = []
            for base_name in base_names:
                mzml_file = next(mzml for mzml in mzmls if mzml.startswith(base_name))
                command = ' '.join(['Rscript ../xcms_quantify_features.R',
                                    f'--mzml {mzml_file}',
                                    f'--output {base_name}.results',
                                    '--xcms_params xcms_params',
                                    '--peakmerge_params merge_params',
                                    f'--algorithm {self.algorithm}'])
                subprocess.run(command, shell = True)
                
                features = pd.read_csv(f'{base_name}.results', sep = '\t').replace(0, np.nan)
                features.columns = ['mz',
                                    'mzmin',
                                    'mzmax',
                                    'rt',
                                    'rt_start',
                                    'rt_end',
                                    'intensity',
                                    'intb',
                                    'maxo',
                                    'sn',
                                    'sample']
                features = features[['rt_start', 'rt_end', 'mz', 'intensity']]
                
                psms = pd.read_csv(f'{base_name}_PSMs.txt', sep = '\t')
                psms['mass'] = psms['Theo. MH+ [Da]'] - pipeline_tools.H
                psms['rt'] = psms['RT [min]']
                psms['sequence'] = psms['Annotated Sequence']
                psms = psms[['mass', 'rt', 'sequence']]  
                
                peptide_results.append(self.peptide_rollup(features, psms))
                
            end = time()
            
            #process results
            quant_depth, mre = self.calc_metrics(peptide_results[0], peptide_results[1])
            runtime = end - start
            
            result_line = list(job.values()) + [quant_depth, mre, runtime]
            with open('../outcomes.tsv', 'a') as tsv:
                tsv.write('\t'.join(str(r) for r in result_line) + '\n')
            
        except Exception as e:
            traceback.print_exc(e)
        finally:
            #clean up temporary files
            os.chdir('..')
            shutil.rmtree(tmpdir)




