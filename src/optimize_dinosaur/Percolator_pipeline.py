#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:29:31 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Percolator(pipeline_tools.Pipeline):
    def __init__(self):
        self.name = 'Percolator'
        self.cores = 3
        self.timeout = '01:00:00'
        self.memory = 12

    def get_params(self):
        choices = {'default-direction':['Xcorr', 
                                        'enzC', 
                                        '-lnExpect', 
                                        'deltCn', 
                                        '-PepLen', 
                                        '-lnrSp', 
                                        'lnNumSP', 
                                        '-enzInt', 
                                        'IonFrac', 
                                        '-Sp', 
                                        'enzN', 
                                        '-absdM'],
                   'testFDR':[0.01, 0.001, 0.05],
                   'trainFDR':[0.01, 0.001, 0.05],
                   'maxiter':[10,5,15,20],
                   'init-weights':[False, 'weights.tsv'],
                   'unitnorm':[False, True],
                   'nested-xval-bins':[1,3,5],
                   'tab-in':['-y /data/S1_N1.pin', 
                             '-y /data/S1_N5.pin', 
                             '-Y /data/S2_N1.pin', 
                             '-Y /data/S2_N5.pin']}
        return choices
    
    def get_metrics(self):
        self.metrics = {'depth':1,
                        'FDR':-1}
        return self.metrics
    
    def setup_workspace(self):
        import subprocess
        import os
        import shutil
        
        #make a percolator singularity container
        subprocess.run('git clone https://github.com/stavis1/proteomics_cluster_submission', shell = True)
        os.chdir('proteomics_cluster_submission/exes/')
        subprocess.run('wget https://github.com/percolator/percolator/releases/download/rel-3-06-05/percolator-noxml-v3-06-linux-amd64.deb', shell = True)
        subprocess.run('singularity build --fakeroot percolator.sif percolator.def', shell = True)
        shutil.move('percolator.sif', '../../')
        os.chdir('../../')
        shutil.rmtree('proteomics_cluster_submission')
        
        #download comet
        subprocess.run('wget https://github.com/UWPR/Comet/releases/download/v2024.01.1/comet.linux.exe', shell = True)
        subprocess.run('chmod +x comet.linux.exe', shell = True)
        
        #run comet jobs
        for search_type in [1,2]:
            for n_outputs in [1,5]:
                with open('comet.params', 'r') as params_in:
                    with open('tmp.params', 'w') as params_out:
                        for line in params_in:
                            if line.startswith('decoy_search'):
                                line = f'decoy_search = {search_type}\n'
                            elif line.startswith('num_output_lines'):
                                line = f'num_output_lines = {n_outputs}'
                            params_out.write(line)
                
                command = ' '.join(['./comet.linux.exe',
                                    '-Ptmp.params',
                                    '-Doptimize_percolator.faa',
                                    f'-NS{search_type}_N{n_outputs}',
                                    'optimize_percolator.mzML'])
                subprocess.run(command, shell = True)
        
    def run_job(self, job):
        super().run_job(job)
        import os
        import shutil
        import subprocess
        import pandas as pd
        import numpy as np
        import time
        import re

        #set up temporary workspace
        temp_dir = str(os.getpid())
        files = [f for f in os.listdir() if os.path.isfile(f)]
        os.mkdir(temp_dir)
        os.chdir(temp_dir)
        for file in files:
            os.link(f'../{file}', file)
        try:
            #run percolator
            singularity_params = '--fakeroot --containall --bind ./:/data/ -w --unsquash'
            perc_params = '-U -m /data/results.pout'
            job_params = ' '.join(f'--{k} {v}' if v != 'True' else f'--{k}' for k,v in job.items() if v != 'False')
            job_params = re.sub(r'--tab-in ', '', job_params)
            command = f'singularity run {singularity_params} percolator.sif percolator {perc_params} {job_params}'
            print(command)
            
            start = time.time()
            subprocess.run(command, shell = True)
            end = time.time()
            
            #parse results
            psms = pd.read_csv('results.pout', sep = '\t')
            psms = psms[psms['q-value'] < 0.01]
            psms['fp'] = [all('ecoli' in p for p in ps.split(';')) for ps in psms['proteinIds']]
            fp = np.sum(psms['fp'])*16.6
            tp = np.sum(np.logical_not(psms['fp']))
            fdr = fp/(tp+fp)
            depth = tp
            runtime = end - start
            
            result_line = list(job.values()) + [depth, fdr, runtime]
            with open('../outcomes.tsv', 'a') as tsv:
                tsv.write('\t'.join(str(r) for r in result_line) + '\n')

        except Exception as e:
            print(e)
        #clean up temporary files
        os.chdir('..')
        shutil.rmtree(temp_dir)


