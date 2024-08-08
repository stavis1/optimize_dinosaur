#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 17:11:13 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

H = 1.007276

class Dinosaur(pipeline_tools.FeatureFinderPipeline):
    def __init__(self):
        self.name = 'Dinosaur'
        self.cores = 4
        self.timeout = '01:00:00'
        self.get_params()
    
    def get_params(self):
        super().get_params()
        corr_set = ['0.6', '0.4', '0.8', '0.2', '0.9', '0.95', '0.0', '0.01']
        dinosaur_params = {'averagineCorr' : corr_set,
                           'averagineExplained' : ['0.5', '0.45', '0.48', '0.51', '0.53', '0.55'],
                           'chargePairPPM' : ['7.0', '2.0', '8.0', '10.0', '15.0', '20.0'],
                           'deisoCorr' : corr_set,
                           'hillMaxMissing' : ['1', '2', '4', '6', '12'],
                           'hillMinLength' : ['3', '1', '2', '6', '12'],
                           'hillPPM' : ['8.0', '5.0', '3.0', '10.0', '20.0'],
                           'hillPeakFactor' : ['2', '1', '3', '4', '5'],
                           'hillPeakFactorMinLength' : ['40', '30', '35', '20', '50', '60', '70'],
                           'hillSmoothMeanWindow' : ['1', '2', '3', '5', '7'],
                           'hillSmoothMedianWindow' : ['1', '2', '3', '5', '7'],
                           'hillValleyFactor' : ['1.3', '1.1', '1.8', '1.01', '2.0', '1.4'],
                           'noHillSplit' : ['false','true']}
        self.dinosaur_param_set = set(dinosaur_params.keys())

        self.param_choices.update(dinosaur_params)
        return self.param_choices

    def setup_workspace(self):
        import os
        import subprocess 
        
        if not os.path.exists('Dinosaur.jar'):
            subprocess.run('wget https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar -O Dinosaur.jar',
                           shell = True)
    
    def run_job(self, job):
        super().run_job(job)
        self.set_params(job)
        
        import os
        import subprocess
        import shutil
        from time import time

        import pandas as pd
        
        #set up temporary workspace
        tmpdir = str(os.getpid())
        os.mkdir(tmpdir)
        mzmls = [f for f in os.listdir() if f.endswith('.mzML')]
        base_names = [f[:-5] for f in mzmls]
        psms = [f for f in os.listdir() if f.endswith('_PSMs.txt')]
        os.chdir(tmpdir)
        try:
            for file in mzmls + psms:
                os.link(f'../{file}', file)
            start = time()
            
            #run Dinosaur    
            with open('dinosaur.params', 'w') as params:
                params.write('\n'.join(f'{k}={v}' for k,v in job.items() if k in self.dinosaur_param_set))
            
            for mzml in mzmls:
                subprocess.run(f'java -jar ../Dinosaur.jar --advParams={os.path.abspath("dinosaur.params")} --concurrency={self.cores} {mzml}', shell = True)
            
            #run peptide rollup
            peptide_results = []
            for base_name in base_names:
                features = pd.read_csv(f'{base_name}.features.tsv', sep = '\t')
                features['rt_start'] = features['rtStart']
                features['rt_end'] = features['rtEnd']
                features['mz'] = (features['mass']/features['charge']) + H
                features['intensity'] = features['intensitySum']/features['charge']
                features = features[['rt_start', 'rt_end', 'mz', 'intensity']]
                
                psms = pd.read_csv(f'{base_name}_PSMs.txt', sep = '\t')
                psms['mass'] = psms['Theo. MH+ [Da]'] - H
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
            print(e)
        #clean up temporary files
        os.chdir('..')
        shutil.rmtree(tmpdir)