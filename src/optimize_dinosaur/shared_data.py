#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:59:19 2024

@author: 4vt
"""

corr_set = [0.6, 0.5, 0.4, 0.8]
params = {'averagineCorr' : corr_set,
          'averagineExplained' : [0.5, 0.3, 0.7],
          'backgroundQuantile' : [0.0, 0.01, 0.05],
          'chargePairCorr' : corr_set,
          'chargePairPPM' : [7.0, 5, 10],
          'deisoCorr' : corr_set,
          'deisoDecompMaxSeeds' : [100, 50, 150],
          'deisoOverlap' : [3, 2, 5],
          'deisoSigmas' : [5.0, 3, 6],
          'deisoValleyFactor' : [1.3, 1.1, 1.8],
          'hillMaxMissing' : [1, 3, 5, 7],
          'hillMinLength' : [3, 5, 7],
          'hillMzGuessLength' : [3, 5, 7],
          'hillNBoots' : [150, 100, 200],
          'hillPPM' : [8.0, 5, 10],
          'hillPeakFactor' : [2, 1.5, 3],
          'hillPeakFactorMinLength' : [40, 30, 50],
          'hillSmoothMeanWindow' : [1, 3, 5],
          'hillSmoothMedianWindow' : [1, 3, 5],
          'hillValleyFactor' : [1.3, 1.1, 1.8],
          'massCalcNBoots' : [150, 100, 200],
          'massEstPoints' : [3, 5, 7],
          'maxBootSize' : [300, 200, 400],
          'maxIntensity' : ['false', 'true'],
          'noHillSplit' : ['false','true']}