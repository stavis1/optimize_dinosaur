#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 13:07:45 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Asari(pipeline_tools.FeatureFinderPipeline):
    def __init__(self):
        self.name = 'Asari'
        self.cores = 1
        self.memory = 95
        self.timeout = '01:00:00'
        self.get_params()

    def get_params(self):
        super().get_params()
        self.run_params = {'mode': 'pos',
                           'project_name': 'asari_optimization',
                           'multicores': self.cores,
                           'outdir': 'output',
                           # 'reference': 'null',
                           'database_mode': 'ondisk'}
        
        asari_params = {'mz_tolerance_ppm': [5, 4, 8, 10, 15],              # ppm, high selectivity meaning no overlap neighbors to confuse; 
                                                            # Low selectivity regions will be still inspected to determine the true number of features
                    
                        # chromatogram and peak parameters
                        'min_timepoints': [6, 5, 4, 3, 7],                # minimal number of data points in elution profile. scipy find_peaks treat `width` as FWHM, thus half of this value.
                        'signal_noise_ratio': [2, 1, 3, 4],            # peak height at least x fold over local noise
                        'min_intensity_threshold': [1000, 10000, 100000],    # minimal intensity for mass track extraction, filtering baseline
                        'min_peak_height': [100000, 10000, 1000000],           # minimal peak height.
                        'min_peak_ratio': [0.001, 0.0001, 0.01, 0.00001],             # minimal ratio of a peak of the max height of its ROI, relevant to small peaks next to big ones.
                        'wlen': [25, 10, 50, 100],                         # window size for evaluating prominence in peaks. Important to resolve clustered narrow peaks.
                        'autoheight': ['false', 'true'],                # min_peak_height can be estimated automatically by setting autoheight on in CLI  
                        'gaussian_shape': [0.5, 0.2, 0.1, 0.8, 0.9],              # min cutoff of goodness of fitting to Gauss model
                        'peak_area': ['sum', 'auc', 'gauss'],                 # `sum` for simple sum, `auc` for area under the curve, `gauss` for gaussian
                        
                        # retention time alignment
                        'rt_align_method': ['lowess', 'tolerance'],        # 'lowess', 'tolerance', or to implement           
                        'rt_align_on': ['true', 'false'],                # False to bypass retention time alignment
                        'rtime_tolerance': [50, 100, 200, 20],              # feature rtime shift threshold under 10 seconds; or 10% of rtime   
                        'cal_min_peak_height': [100000, 10000, 1000000],      # minimal peak height required for peaks used for RT calibration
                        'max_retention_shift': ['null',10,20,50,100],        # landmark peak pairs with a scan number delta greater than this are not used for RT calibration
                        'num_lowess_iterations': [3,1]}
        self.asari_param_set = set(asari_params.keys())
        self.param_choices.update(asari_params)
        
        self.asari_param_dtypes = {'mz_tolerance_ppm': 'int',
                                   'min_timepoints': 'int',
                                   'signal_noise_ratio': 'float',
                                   'min_intensity_threshold': 'float',
                                   'min_peak_height': 'float',
                                   'min_peak_ratio': 'float',
                                   'wlen': 'int',
                                   'autoheight': 'bool',
                                   'gaussian_shape': 'float',
                                   'peak_area': 'str',
                                   'rt_align_method': 'str',
                                   'rt_align_on': 'bool',
                                   'rtime_tolerance': 'float',
                                   'cal_min_peak_height': 'float',
                                   'max_retention_shift': 'float',
                                   'num_lowess_iterations': 'int'}
        
        return self.param_choices

    def setup_workspace(self):
        import subprocess
        subprocess.run('conda create -y -n asari_env "python=3.12.5" pip -c conda-forge', shell = True)
        subprocess.run('conda run -n asari_env pip install "asari-metabolomics==1.13.1"', shell = True)

    def run_job(self, job):
        super().run_job(job)
        
        import subprocess
        from time import time
        import os
        import pandas as pd
        import shutil
        
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
                
            with open('asari.params', 'w') as yaml:
                for param, value in self.run_params.items():
                    yaml.write(f'{param}: {value}\n')
                for param in self.asari_param_set:
                    if not self.params[param] == 'null':
                        yaml.write(f'{param}: !!{self.asari_param_dtypes[param]} {self.params[param]}\n')
            
            start = time()
            subprocess.run('conda run -n asari_env asari process -p asari.params -i ./', shell = True)
            
            #run peptide rollup
            peptide_results = []
            outdir = next(f for f in os.listdir() if f.startswith('output_'))
            features = pd.read_csv(os.path.join(outdir,'export/full_Feature_table.tsv'), sep = '\t')
            features = features[['rtim_left_base', 'rtime_right_base', 'mz'] + base_names]
            features.columns = ['rt_start', 'rt_end', 'mz'] + base_names
            for base_name in base_names:
                feature_subset = features[features[base_name] > 0]
                feature_subset['intensity'] = feature_subset[base_name]
                
                psms = pd.read_csv(f'{base_name}_PSMs.txt', sep = '\t')
                psms['mass'] = psms['Theo. MH+ [Da]'] - pipeline_tools.H
                psms['rt'] = psms['RT [min]']
                psms['sequence'] = psms['Annotated Sequence']
                psms = psms[['mass', 'rt', 'sequence']]
                
                peptide_results.append(self.peptide_rollup(feature_subset, psms))
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
