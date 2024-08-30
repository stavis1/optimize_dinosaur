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
        shutil.copy2(tool, 'xcms_quantify_features.R')
        
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
    
    def write_toml(self, data, out_path):
        '''
        writes a TOML file formatted for XCMS parameters
        takes:
            data: a dictionary of parameters
            out_path: the file name for the output toml
        returns:
            None
        '''
        toml = [f'{key} = [ {value},]' for key, value in data.items()]
        with open(out_path, 'w') as toml_file:
            toml_file.write('\n'.join(toml))        
    
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
                
            #make params files
            self.write_toml(dict(i for i in job.items() if i[0] in self.xcms_param_set),
                            'xcms_params')
            self.write_toml(dict(i for i in job.items() if i[0] in self.merge_param_set),
                            'merge_params')

            start = time()
            peptide_results = []
            for base_name in base_names:
                mzml_file = next(mzml for mzml in mzmls if mzml.startswith(base_name))
                command = ' '.join(['conda run -n xcms_env',
                                    'Rscript ../xcms_quantify_features.R',
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
                features['rt_start'] = features['rt_start']/60
                features['rt_end'] = features['rt_end']/60
                
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

class Xcms_cw(Xcms_base):
    def __init__(self):
        super().__init__()
        self.name = 'Xcms_cw'
        self.algorithm = 'xcms_cw'
    
    def get_params(self):
        params = super().get_params()
        xcms_params = {'ppm':[25,20,15,10,8,5,4],
                       'peakwidth':[f'{l}, {u}' for l in (20, 30, 40, 15, 10) for u in (50, 40, 70, 90, 120, 150, 300)],
                       'snthresh':[10,15,8,5,4,3,2],
                       'prefilter':[f'{l}, {u}' for l in (3, 2, 1, 4, 5, 6) for u in (100, 1000, 10000, 10)],
                       'mzCenterFun':['"wMean"', '"mean"', '"apex"', '"wMeanApex3"', '"meanApex3"'],
                       'integrate':[1,2],
                       'mzdiff':[-0.001, -0.005, -0.01, 0.001, 0.005, 0.01, 0],
                       'fitgauss':['false','true'],
                       'noise':[0, 10, 100, 1000, 10000],
                       'firstBaselineCheck':['true','false'],
                       'extendLengthMSW':['false','true']}
        self.xcms_param_set = set(xcms_params.keys())
        params.update(xcms_params)
        return params

class Xcms_cwip(Xcms_base):
    def __init__(self):
        super().__init__()
        self.name = 'Xcms_cwip'
        self.algorithm = 'xcms_cwip'
    
    def get_params(self):
        params = super().get_params()
        xcms_params = {'ppm':[25,20,15,10,8,5,4],
                       'peakwidth':[f'{l}, {u}' for l in (20, 30, 40, 15, 10) for u in (50, 40, 70, 90, 120, 150, 300)],
                       'snthresh':[10,15,8,5,4,3,2],
                       'prefilter':[f'{l}, {u}' for l in (3, 2, 1, 4, 5, 6) for u in (100, 1000, 10000, 10)],
                       'mzCenterFun':['"wMean"', '"mean"', '"apex"', '"wMeanApex3"', '"meanApex3"'],
                       'integrate':[1,2],
                       'mzdiff':[-0.001, -0.005, -0.01, 0.001, 0.005, 0.01, 0],
                       'fitgauss':['false','true'],
                       'noise':[0, 10, 100, 1000, 10000],
                       'firstBaselineCheck':['true','false'],
                       'extendLengthMSW':['false','true'],
                       'snthreshIsoROIs':[6.25, 5, 3, 2, 1.5, 7, 8],
                       'maxCharge':[3, 4, 5],
                       'maxIso':[5, 4, 3, 2, 6],
                       'mzIntervalExtension':['true', 'false']}
        self.xcms_param_set = set(xcms_params.keys())
        params.update(xcms_params)
        return params

class Xcms_mf(Xcms_base):
    def __init__(self):
        super().__init__()
        self.name = 'Xcms_mf'
        self.algorithm = 'xcms_mf'
    
    def get_params(self):
        params = super().get_params()
        xcms_params = {'binSize':[0.1,0.05,0.01,0.15,0.2],
                       'impute':['"none"','"lin"','"linbase"','"intlin"'],
                       'baseValue':[0, 1, 10, 100, 1000],
                       'distance':[1,2,3,4,5],
                       'sigma':[12.73994, 10, 8, 6, 15, 18, 20],
                       'max':[5, 4, 3, 6, 7, 8, 9, 10],
                       'snthres':[10,15,8,5,4,3,2],
                       'steps':[2, 1, 3, 4, 5],
                       'mzdiff':[0.6, 0.4, 0.2, 0.05, 0.8, 1]}
        self.xcms_param_set = set(xcms_params.keys())
        params.update(xcms_params)
        return params

class Xcms_kalman(Xcms_base):
    def __init__(self):
        super().__init__()
        self.name = 'Xcms_kalman'
        self.algorithm = 'xcms_kalman'
    
    def get_params(self):
        params = super().get_params()
        xcms_params = {'ppm':[25,20,15,10,8,5,4],
                       'peakwidth':[f'{l}, {u}' for l in (20, 30, 40, 15, 10) for u in (50, 40, 70, 90, 120, 150, 300)],
                       'snthresh':[10,15,8,5,4,3,2],
                       'prefilter':[f'{l}, {u}' for l in (3, 2, 1, 4, 5, 6) for u in (100, 1000, 10000, 10)],
                       'mzCenterFun':['"wMean"', '"mean"', '"apex"', '"wMeanApex3"', '"meanApex3"'],
                       'integrate':[1,2],
                       'mzdiff':[-0.001, -0.005, -0.01, 0.001, 0.005, 0.01, 0],
                       'fitgauss':['false','true'],
                       'noise':[0, 10, 100, 1000, 10000],
                       'criticalValue':[1.125, 0.1, 0.5, 1.5, 2, 2.5, 3],
                       'consecMissedLimit':[2,1,3],
                       'unions':[1,0],
                       'checkBack':[0,1],
                       'withWave':['false','true']}
        self.xcms_param_set = set(xcms_params.keys())
        params.update(xcms_params)
        return params


