#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:34:42 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Flashlfq(pipeline_tools.PepQuantPipeline):
    def __init__(self):
        self.name = 'Flashlfq'
        self.cores = 2
        self.timeout = '03:00:00'
        self.memory = 8
    
    def get_params(self):
        param_choices = {'nor':['false', 'true'],
                         'ppm':[10, 8, 5, 3, 15, 20],
                         'iso':[5, 3, 8, 10, 15, 20],
                         'int':['false', 'true'],
                         'nis':[2,1,3,4],
                         'chg':['false', 'true'],
                         'mbr':['false', 'true'],
                         'mrt':[2.5, 2, 1, 3, 4, 5, 7, 10]}
        return param_choices
    
    def set_params(self, params):
        super().set_params(params)
        self.params.update({'idt':'/data/psms.tsv',
                            'rep':'/data/',
                            'out':'/data/',
                            'thr':2})
    
    def setup_workspace(self):
        import subprocess
        import os
        import pandas as pd
        
        subprocess.run('singularity build flashlfq.sif docker://smithchemwisc/flashlfq:latest', shell = True)
        
        flashlfq_cols = ['File Name', 
                         'Base Sequence',
                         'Full Sequence',
                         'Peptide Monoisotopic Mass',
                         'Scan Retention Time',
                         'Precursor Charge',
                         'Protein Accession']
        psm_files = [f for f in os.listdir() if f.endswith('_PSMs.txt')]
        all_psms = []
        for psm_file in psm_files:
            psms = pd.read_csv(psm_file, sep = '\t')
            psms['File Name'] = [f[:-4] + '.mzML' for f in psms['Spectrum File']]
            psms['Base Sequence'] = psms['Sequence']
            psms['Full Sequence'] = psms['Annotated Sequence']
            psms['Peptide Monoisotopic Mass'] = psms['Theo. MH+ [Da]']# - pipeline_tools.H
            psms['Scan Retention Time'] = psms['RT [min]']
            psms['Precursor Charge'] = psms['Charge']
            psms['Protein Accession'] = psms['Protein Accessions']
            all_psms.append(psms[flashlfq_cols])
        all_psms = pd.concat(all_psms)
        all_psms.to_csv('psms.tsv', sep = '\t', index = False)
            
    def run_job(self, job):
        super().run_job(job)

        import os
        import subprocess
        import shutil
        from time import time
        import traceback

        import pandas as pd
        
        #set up temporary workspace
        tmpdir = str(os.getpid())
        print(tmpdir, flush = True)
        os.mkdir(tmpdir)
        mzmls = [f for f in os.listdir() if f.endswith('.mzML') and f.startswith('20210827')]
        psms = [f for f in os.listdir() if f.endswith('_PSMs.txt')]
        os.chdir(tmpdir)
        try:
            for file in mzmls + psms:
                os.link(f'../{file}', file)
            os.link('../psms.tsv', 'psms.tsv')
            
            flfq_params = ' '.join(f'--{k}={v}' if v != 'true' else '--{k}' for k,v in self.params.items() if v != 'false')
            command = f'singularity run --bind ./:/data/ --containall ../flashlfq.sif {flfq_params}'
            print(command, flush = True)
            start = time()
            subprocess.run(command, shell = True)
            end = time()
            
            #process results
            peptides = pd.read_csv('QuantifiedPeptides.tsv', sep = '\t')
            peptide_results = []
            for mzml in mzmls:
                mzml_col = f'Intensity_{mzml[:-5]}'
                subset = peptides[mzml_col] > 0
                subset['sequence'] = subset['Sequence']
                subset['intensity'] = subset[mzml_col]
                peptide_results.append(subset[['sequence', 'intensity']])
            quant_depth, mre = self.calc_metrics(peptide_results[0], peptide_results[1])
            runtime = end - start
            
            result_line = list(job.values()) + [quant_depth, mre, runtime]
            with open('../outcomes.tsv', 'a') as tsv:
                tsv.write('\t'.join(str(r) for r in result_line) + '\n')
            
        except Exception as e:
            traceback.print_exc(e)
        #clean up temporary files
        os.chdir('..')
        shutil.rmtree(tmpdir)
