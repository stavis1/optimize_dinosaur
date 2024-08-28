#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 17:07:22 2024

@author: 4vt
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--mzml', action = 'store', required = True,
                    help = 'the .mzML spetrum file')
parser.add_argument('--psms', action = 'store', required = True,
                    help = 'the psm identification file from Proteome Discoverer')
parser.add_argument('--output', action = 'store', required = True,
                    help = 'the output file name')
parser.add_argument('--params', action = 'store', required = True,
                    help = 'the tsv of feature detection parameters')
args = parser.parse_args()

import re
import pyopenms as oms
import pandas as pd
import numpy as np
from sortedcontainers import SortedList
from collections import defaultdict

psms = pd.read_csv(args.psms, sep = '\t')
params = pd.read_csv(args.params, sep = '\t', header = None)
omm_idx = next(i for i,p in zip(params.index, params.iloc[:,0]) if p == 'onMultiMatch')
args.onMultiMatch = params.iloc[omm_idx, 1]
params = params.drop(omm_idx)

# Prepare data loading (save memory by only
# loading MS1 spectra into memory)
options = oms.PeakFileOptions()
options.setMSLevels([1])
fh = oms.MzMLFile()
fh.setOptions(options)

# Load data
input_map = oms.MSExperiment()
fh.load(args.mzml, input_map)
input_map.updateRanges()

ff = oms.FeatureFinder()
ff.setLogType(oms.LogType.CMD)

# Run the feature finder
name = "centroided"
features = oms.FeatureMap()
seeds = oms.FeatureMap()
ff_params = oms.FeatureFinder().getParameters(name)
for key, val in zip(params.iloc[:,0],params.iloc[:,1]):
    if re.search(r'[^\d\.]', val) is None:
        if '.' in val:            
            val = float(val)
        else:
            val = int(val)
    ff_params.setValue(key, val)
ff.run(name, input_map, features, ff_params, seeds)
features.setUniqueIds()

def get_data(feature, idx):
    windows = []
    for hull in feature.getConvexHulls():
        hull = hull.getBoundingBox2D()
        #rtstart, rtend, mzlow, mzhigh, intensity, index
        windows.append((hull[0][0]/60, 
                        hull[1][0]/60,
                        hull[0][1], 
                        hull[1][1],
                        feature.getIntensity(),
                        idx))
    return windows

feature_data = [w for i,f in enumerate(features) for w in get_data(f, i)]

#set up data structures for faster lookups
intensity = {f[5]:f[4] for f in feature_data}
rt_start = SortedList((f[0], f[5]) for f in feature_data)
rt_end = SortedList((f[1], f[5]) for f in feature_data)
max_Δrt = max(f[1]-f[0] for f in feature_data)
mz_low = SortedList((f[2], f[5]) for f in feature_data)
mz_high = SortedList((f[3], f[5]) for f in feature_data)
max_Δmz = max(f[3]-f[2] for f in feature_data)

def match_feature(mz, rt):
    rt_starts = set(f[1] for f in rt_start.irange((rt-max_Δrt,), (rt,)))
    rt_ends = set(f[1] for f in rt_end.irange((rt,), (rt+max_Δrt,)))
    mz_starts = set(f[1] for f in mz_low.irange((mz-max_Δmz,), (mz,)))
    mz_ends = set(f[1] for f in mz_high.irange((mz,), (mz+max_Δmz,)))
    
    matches = rt_starts
    for match_set in (rt_ends, mz_starts, mz_ends):
        matches = matches.intersection(match_set)
    
    return matches

def peptide_rollup(feature_sets):
    all_features = []
    for feature_set in feature_sets:
        if len(feature_set) > 1:
            if args.onMultiMatch == 'drop':
                continue
            elif args.onMultiMatch == 'sum':
                all_features.extend(feature_set)
            elif args.onMultiMatch == 'max':
                all_features.append(max(feature_set, key = lambda x: intensity[x]))
        elif feature_set:
            all_features.extend(feature_set)
    all_features = set(all_features)
    pep_intensity = np.sum([intensity[f] for f in all_features])
    return pep_intensity

feature_matches = defaultdict(lambda:[])
for seq, mz, rt in zip(psms['Annotated Sequence'], psms['m/z [Da]'], psms['RT [min]']):
    feature_matches[seq].append(match_feature(mz, rt))
    
peptides = list(set(psms['Annotated Sequence']))
intensities = [peptide_rollup(feature_matches[s]) for s in peptides]

report = pd.DataFrame({'sequence':peptides,
                       'intensity':intensities})
report.to_csv(args.output, sep = '\t', index = False)

