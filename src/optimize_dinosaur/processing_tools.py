#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:35:44 2024

@author: 4vt
"""
def map_feature(psm_idx):
    rt = psm_rt[psm_idx]
    mass = psm_mass[psm_idx]
    ppm = (mass/1e6)*float(job_params['ppm'])
    rtstart_set = set((i[1] for i in rtstart_idx.irange((rt-max_rt_width,), (rt,))))
    rtend_set = set((i[1] for i in rtend_idx.irange((rt,), (rt+max_rt_width,))))
    rt_set = rtstart_set.intersection(rtend_set)
    if job_params['add_H']:
        mass_set = set((i[1] for i in mass_H_idx.irange((mass-ppm,),(mass+ppm,))))
    else:
        mass_set = set((i[1] for i in mass_idx.irange((mass-ppm,),(mass+ppm,))))
    feature_set = rt_set.intersection(mass_set)
    return feature_set

def peptide_rollup(base_name, job):
    global job_params
    job_params = job
    from multiprocessing import Pool
    from collections import defaultdict
    from sortedcontainers import SortedList
    import pandas as pd
    import numpy as np
    
    H = 1.007276

    #read in data
    charges = [int(c) for c in job['charges'].split(',')]
    features = pd.read_csv(f'{base_name}.features.tsv', sep = '\t')
    features = features[[c in charges for c in features['charge']]]
    psms = pd.read_csv(f'{base_name}_PSMs.txt', sep = '\t')

    #build psm index dictionaries for fast lookup
    global psm_rt
    psm_rt = {i:r for i,r in zip(psms.index, psms['RT [min]'])}
    global psm_mass
    psm_mass = {i:m for i,m in zip(psms.index, psms['Theo. MH+ [Da]'] - H)}
    seq_prots = {s:p for s,p in zip(psms['Annotated Sequence'], psms['Protein Accessions'])}

    #build feature indices for fast lookup
    global rtstart_idx
    rtstart_idx = SortedList(zip(features['rtStart'], features.index))
    global rtend_idx
    rtend_idx = SortedList(zip(features['rtEnd'], features.index))
    global mass_idx
    mass_idx = SortedList(zip(features['mass'], features.index))
    global mass_H_idx
    mass_H_idx = SortedList(zip(features['mass'] - (H*features['charge']), features.index))
    features['intensityCorr'] = features['intensitySum']/features['charge']
    feature_intensity = {idx:i for idx,i in zip(features.index, features['intensityCorr'])}
    
    global max_rt_width
    max_rt_width = max(features['rtEnd'] - features['rtStart'])

    #parallelize map_feature() calls
    with Pool(8) as p:
        feature_map = p.map(map_feature, range(psms.shape[0]))

    #peptide rollup 
    class peptide():
        def __init__(self, seq):
            self.seq = seq
            self.psm_indices = []
            self.features = set([])
            self.proteins = seq_prots[self.seq]
        
        def add_psm(self, psm_index, features):
            self.psm_indices.append(psm_index)
            self.features.update(features)
        
        def remove_bad_features(self, bad_features):
            self.features = [f for f in self.features if f not in bad_features]
            
        def calculate_intensity(self, intensity_map):
            self.intensity = np.sum([intensity_map[f] for f in self.features])
        
        def report(self):
            return (self.seq, 
                    self.intensity, 
                    ';'.join((str(i) for i in self.psm_indices)), 
                    ';'.join((str(i) for i in self.features)),
                    self.proteins)

    class keydefaultdict(defaultdict):
        def __missing__(self, key):
            if self.default_factory is None:
                raise KeyError( key )
            else:
                ret = self[key] = self.default_factory(key)
                return ret

    peptides = keydefaultdict(peptide)
    for seq, psm, feature_set in zip(psms['Annotated Sequence'], psms.index, feature_map):
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
        peptide.calculate_intensity(feature_intensity)

    #report
    peptide_data = pd.DataFrame(np.array([p.report() for p in peptide_list]),
                                columns = ('sequence', 'intensity', 'psm_indices', 'feature_indices', 'proteins'))
    peptide_data.to_csv(f'{base_name}.peptides.txt', sep = '\t', index = False)

def process_results(base_names, job):
    from collections import defaultdict
    import numpy as np
    import pandas as pd
    
    peptide_quants = defaultdict(lambda:[np.nan]*2)
    allseqs = set([])
    for i,base_name in enumerate(base_names):
        PSM_data = pd.read_csv(f'{base_name}_PSMs.txt', sep = '\t')
        allseqs.update(PSM_data['Annotated Sequence']) 
        
        peptide_data = pd.read_csv(f'{base_name}.peptides.txt', sep = '\t')
        for seq, intensity in zip(peptide_data['sequence'], peptide_data['intensity']):
            peptide_quants[seq][i] = intensity
    peptide_quants = np.array(list(peptide_quants.values()))
    depth = np.sum(np.isnan(peptide_quants))/len(allseqs)
    jaccard = np.sum(np.all(np.isfinite(peptide_quants), axis = 1))/peptide_quants.shape[0]
    peptide_quants = peptide_quants[np.all(np.isfinite(peptide_quants), axis = 1),:]
    means = np.mean(peptide_quants, axis = 1)
    diffs = np.abs(peptide_quants[:,0] - peptide_quants[:,1])
    mre = np.mean(diffs/means)
    result_line = list(job.values()) + [mre, jaccard, depth]
    with open('../outcomes.tsv', 'a') as tsv:
        tsv.write('\t'.join(str(r) for r in result_line) + '\n')

def run_job(job):
    import os
    import subprocess
    import shutil
    
    from optimize_dinosaur.shared_data import dinosaur_param_set
    
    #write job paramters
    with open('attempted_solutions.tsv', 'a') as tsv:
        tsv.write('\t'.join(str(v) for v in job.values()) + '\n')
    
    #set up temporary workspace
    tmpdir = str(os.getpid())
    os.mkdir(tmpdir)
    mzmls = [f for f in os.listdir() if f.endswith('.mzML')]
    base_names = [f[:-5] for f in mzmls]
    psms = [f for f in os.listdir() if f.endswith('_PSMs.txt')]
    os.chdir(tmpdir)
    for file in mzmls + psms:
        os.link(f'../{file}', file)
    
    
    #run Dinosaur    
    with open('dinosaur.params', 'w') as params:
        params.write('\n'.join(f'{k}={v}' for k,v in job.items() if k in dinosaur_param_set))
    
    for mzml in mzmls:
        subprocess.run(f'java -jar ../Dinosaur.jar --advParams={os.path.abspath("dinosaur.params")} --concurrency=8 {mzml}', shell = True)
    
    #run peptide rollup
    for name in base_names:
        peptide_rollup(name, job)
    
    #process results
    process_results(base_names, job)
    
    #clean up temporary files
    os.chdir('..')
    shutil.rmtree(tmpdir)




