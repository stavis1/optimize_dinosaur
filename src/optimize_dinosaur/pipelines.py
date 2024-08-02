#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:45:57 2024

@author: 4vt
"""
from collections import defaultdict
H = 1.007276

class keydefaultdict(defaultdict):
    '''
    subclass of defaultdict that passes the key to the first
    argument of the default function
    '''
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

class Pipeline():
    def __init__(self):
        self.name = NotImplemented
        self.cores = NotImplemented
        self.timeout = NotImplemented
    
    def get_params(self):
        '''
        takes no arguments
        returns a dictionary of {parameter:[value options]}
        the first element of the options list is used as the default
        '''
        raise NotImplementedError()
    
    def set_params(self, params):
        '''
        takes a dictionary of {parameter:value}
        '''
        self.params = params
    
    def get_metrics(self):
        '''
        takes no arguments
        returns a dictionary of {metric:1 if larger is good else -1}
        '''
        raise NotImplementedError()
    
    def setup_workspace(self):
        '''
        takes no arguments
        runs any pipeline specific code for setting up the workspace
        this includes installing software
        '''
        pass
    
    def run_job(self, job):
        '''
        the only argument is a dictionary of parameter choices
        runs the pipeline and appends the outcome to outcomes.tsv
        '''
        with open('attempted_solutions.tsv', 'a') as tsv:
            tsv.write('\t'.join(str(v) for v in job.values()) + '\n')

class PepQuantPipeline(Pipeline):
    def get_metrics(self):
        return {'quant_depth':1,
                'mean_relative_error':-1}

    def calc_metrics(self, quant1, quant2):
        '''
        arguments:
            takes two dataframes that are the output of self.peptide_rollup()
        returns:
            (quant_depth, mean_relative_error)
        '''
        from collections import defaultdict
        import numpy as np
        
        #build array of intensity values
        quants = defaultdict(lambda:[np.nan]*2)
        for i,table in enumerate((quant1, quant2)):
            for sequence, intensity in zip(table['sequence'], table['intensity']):
                quants[sequence][i] = intensity
        quants = np.array(list(quants.values()))
        
        #calculate the number of quantified peptides
        quant_depth = np.sum(np.isfinite(quants))
        
        #calculate mean relative error
        quants = quants[np.all(np.isfinite(quants), axis = 1)]
        means = np.mean(quants, axis = 1)
        diffs = np.abs(quants[:,0] - quants[:,1])
        mre = np.mean(diffs/means)

        return (quant_depth, mre)
    
    def get_params(self):
        self.param_choices = {'ppm':[5, 2, 8, 10, 15, 20],
                              'rt_wiggle':[0, 0.01, 0.05, 0.1,]}
        self.pep_rollup_param_set = set(self.param_choices.keys())

class FeatureFinderPipeline(PepQuantPipeline):
    def map_feature(self, psm_idx):
        feature_set = set([])
        for charge in range(1,6):
            rt = psm_rt[psm_idx]
            mz = psm_mass[psm_idx]/charge + H
            ppm = (mz/1e6)*self.params['ppm']
            rtstart_set = set((i[1] for i in rtstart_idx.irange((rt-max_rt_width,), (rt,))))
            rtend_set = set((i[1] for i in rtend_idx.irange((rt,), (rt+max_rt_width,))))
            rt_set = rtstart_set.intersection(rtend_set)
            mz_set = set((i[1] for i in mz_idx.irange((mz-ppm,),(mz+ppm,))))
            feature_set.update(rt_set.intersection(mz_set))
        return feature_set
    
    def peptide_rollup(self, features, psms):
        '''
        arguments:
            features (a dataframe) must have columns: rt_start, rt_end, mz, intensity
            psms (a dataframe) must have columns: mass, rt, sequence
        returns:
            a dataframe with columns: sequence, intensity
        note that retention time should be in minutes
        '''
        from multiprocessing import Pool

        from sortedcontainers import SortedList
        import numpy as np
        import pandas as pd
        
        #set up indexes as globals to use in parallelized map_feature() calls
        global psm_rt
        psm_rt = {i:rt for i,rt in zip(psms.index, psms['rt'])}
        global psm_mass
        psm_mass = {i:mass for i,mass in zip(psms.index, psms['mass'])}
        global rtstart_idx
        rtstart_idx = SortedList(zip(features['rt_start'] - self.params['rt_wiggle'], features.index))
        global rtend_idx
        rtend_idx = SortedList(zip(features['rt_end'] + self.params['rt_wiggle'], features.index))
        global max_rt_width
        max_rt_width = max(features['rt_end'] - features['rt_start'])
        global mz_idx
        mz_idx = SortedList(zip(features['mz'], features.index))
        intensity_map = {idx:intensity for idx, intensity in zip(features.index, features['intensity'])}
        
        #connect features to PSMs
        with Pool(self.cores) as p:
            feature_map = p.map(self.map_feature, psms.index)
        
        class peptide():
            def __init__(self, seq):
                self.seq = seq
                self.psm_indices = []
                self.features = set([])
            
            def add_psm(self, psm_index, features):
                self.psm_indices.append(psm_index)
                self.features.update(features)
            
            def remove_bad_features(self, bad_features):
                self.features = [f for f in self.features if f not in bad_features]
                
            def calculate_intensity(self, intensity_map):
                self.intensity = np.sum([intensity_map[f] for f in self.features])
            
            def report(self):
                return (self.seq, 
                        self.intensity)
        
        #initialize peptide objects
        peptides = keydefaultdict(peptide)
        for seq, psm, feature_set in zip(psms['sequence'], psms.index, feature_map):
            if feature_set:
                peptides[seq].add_psm(psm, feature_set)

        #remove degenerate features
        feature_peptides = defaultdict(lambda:[])
        for peptide in peptides.values():
            for feature in peptide.features:
                feature_peptides[feature].append(peptide.seq)
        bad_features = set(f for f,p in feature_peptides.items() if len(p) > 1)

        for peptide in peptides.values():
            peptide.remove_bad_features(bad_features)
        peptide_list = [pep for pep in peptides.values() if pep.features]

        #calculate intensity
        for peptide in peptide_list:
            peptide.calculate_intensity(intensity_map)
        
        #make results dataframe
        peptide_data = pd.DataFrame([p.report() for p in peptide_list],
                                    columns = ('sequence', 'intensity'))
        return peptide_data

class DinosaurRunner(FeatureFinderPipeline):
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

