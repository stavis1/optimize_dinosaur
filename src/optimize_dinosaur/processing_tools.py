#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:35:44 2024

@author: 4vt
"""

def peptide_rollup(base_name, params):
    pass

def process_results(base_names):
    pass

def run_job(job):
    import os
    import subprocess
    import shutil
    
    from shared_data import dinosaur_param_set
    
    #write job paramters
    with open('attempted_solutions.tsv', 'a') as tsv:
        tsv.write('\t'.join(str(v) for v in job.values()) + '\n')
    
    #set up temporary workspace
    tmpdir = str(os.getpid())
    os.mkdir(tmpdir)
    mzmls = [f for f in os.listdir() if f.endswith('.mzML')]
    base_names = [f[:-5] for f in mzmls]
    psms = [f for f in os.listdir() if f.endswith('.pout')]
    os.chdir(tmpdir)
    for file in mzmls + psms:
        os.link(f'../{file}', file)
    
    
    #run Dinosaur
    def flag(k,v):
        if type(v) == bool:
            return f'--{k}' if v else ''
        else:
            return f'--{k} {v}'
    
    dinosaur_param_string = ' '.join(flag(k,v) for k,v in job.items() if k in dinosaur_param_set)
    
    for mzml in mzmls:
        subprocess.run(f'java -jar ../Dinosaur.jar {dinosaur_param_string} {mzml}', shell = True)
    
    #run peptide rollup
    for name in base_names:
        peptide_rollup(name, job)
    
    #process results
    process_results(base_names)
    
    #clean up temporary files
    os.chdir('..')
    shutil.rmtree(tmpdir)
    pass




