#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:20:57 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Pyopenms(pipeline_tools.PepQuantPipeline):
    def __init__(self):
        self.name = 'Pyopenms'
        self.cores = 2
        self.memory = 8
        self.timeout = '08:00:00'
    
    def get_params(self):
        params = {'onMultiMatch':['drop','sum','max'],
                  'intensity:bins':[5, 1, 3, 8, 10, 15, 20],
                  'mass_trace:mz_tolerance':[0.03, 0.01, 0.02, 0.04, 0.06],
                  'mass_trace:min_spectra':[10, 8, 5, 3, 15],
                  'mass_trace:max_missing':[1, 0, 2, 3],
                  'mass_trace:slope_bound':[0.1, 0.01, 0.2, 0.4, 0.8, 1, 1.2],
                  'isotopic_pattern:charge_low':[1, 2],
                  'isotopic_pattern:charge_high':[4, 3],
                  'isotopic_pattern:mz_tolerance':[0.03, 0.01, 0.02, 0.04, 0.06],
                  'isotopic_pattern:intensity_percentage':[10.0, 5, 1, 20],
                  'isotopic_pattern:intensity_percentage_optional':[0.1, 1, 5, 10],
                  'isotopic_pattern:optional_fit_improvement':[2.0, 0.1, 5, 10],
                  'isotopic_pattern:mass_window_width':[25.0, 1, 10, 30, 50, 100],
                  'seed:min_score':[0.8, 0, 0.2, 0.4, 1],
                  'feature:min_score':[0.7, 0, 0.2, 0.4, 1],
                  'feature:min_isotope_fit':[0.8, 0, 0.2, 0.4, 1],
                  'feature:min_trace_score':[0.5, 0, 0.2, 0.4, 0.8, 1],
                  'feature:min_rt_span':[0.333, 0, 0.2, 0.4, 0.8, 1],
                  'feature:max_rt_span':[2.5, 1, 1.5, 3, 4],
                  'feature:rt_shape':["symmetric","asymmetric"],
                  'feature:max_intersection':[0.35, 0, 0.2, 0.4, 0.8, 1]}
        return params

    def setup_workspace(self):
        import subprocess
        import os
        import shutil
        
        #copy tool script to working directory
        tool = os.path.split(os.path.abspath(pipeline_tools.__file__))[0]
        tool = os.path.join(os.path.split(tool)[0], 'tools/pms_quantify_peptides.py')
        shutil.copy2(tool, 'pms_quantify_peptides.py')
        
        #make conda environment for running tool
        if not os.path.exists('~/.conda/envs/pyopenms_env'):
            subprocess.run(' '.join(['conda create -n pyopenms_env', 
                                     'python=3.11.9',
                                     'pyopenms=3.1.0',
                                     'pandas=2.2.2', 
                                     'numpy=1.23.5', 
                                     'sortedcontainers=2.4.0', 
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
            with open('params', 'w') as params:
                params.write('\n'.join([f'{k}\t{v}' for k,v in job.items()]))
            
            start = time()
            peptide_results = []
            for base_name in base_names:
                mzml_file = next(mzml for mzml in mzmls if mzml.startswith(base_name))
                psm_file = next(psm for psm in psms if psm.startswith(base_name))
                command = ' '.join(['conda run -n pyopenms_env',
                                    'python pms_quantify_peptides.py',
                                    f'--mzml {mzml_file}',
                                    f'--psms {psm_file}',
                                    '--params params',
                                    f'--output {base_name}.results'])
                subprocess.run(command, shell = True)
                peptide_results.append(pd.read_csv(f'{base_name}.results', sep = '\t'))
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

