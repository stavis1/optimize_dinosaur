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
        self.set_params(job)

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


class FeatureFinderPipeline(PepQuantPipeline):    
    def get_params(self):
        self.param_choices = {'ppm':[5, 2, 8, 10, 15, 20],
                              'rt_wiggle':[0, 0.01, 0.05, 0.1,]}
        self.pep_rollup_param_set = set(self.param_choices.keys())
        
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


