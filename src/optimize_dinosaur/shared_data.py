#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:59:19 2024

@author: 4vt
"""
from copy import copy

corr_set = ['0.6', '0.4', '0.8', '0.2', '0.9', '0.95', '0.0']
dinosaur_params = {'averagineCorr' : corr_set,
                   'averagineExplained' : ['0.5', '0.45', '0.48', '0.51', '0.53', '0.55'],
                   'chargePairPPM' : ['7.0', '1.0', '20.0'],
                   'deisoCorr' : corr_set,
                   'hillMaxMissing' : ['1', '6', '12'],
                   'hillMinLength' : ['3', '1', '2', '6', '12'],
                   'hillPPM' : ['8.0', '5.0', '3.0', '10.0', '20.0'],
                   'hillPeakFactor' : ['2', '1', '3', '4', '5'],
                   'hillPeakFactorMinLength' : ['40', '30', '35', '20', '50', '60', '70'],
                   'hillSmoothMeanWindow' : ['1', '2', '3', '5', '7'],
                   'hillSmoothMedianWindow' : ['1', '2', '3', '5', '7'],
                   'hillValleyFactor' : ['1.3', '1.1', '1.8', '1.01', '2.0', '1.4'],
                   'noHillSplit' : ['false','true']}
dinosaur_param_set = set(dinosaur_params.keys())

pep_rollup_params = {'ppm':['5', '2', '8', '10', '15', '20']}
pep_rollup_param_set = set(pep_rollup_params.keys())

params = copy(dinosaur_params)
params.update(pep_rollup_params)

