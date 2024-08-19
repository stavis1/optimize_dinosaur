#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:24:12 2024

@author: 4vt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 13:07:45 2024

@author: 4vt
"""

from optimize_dinosaur import pipeline_tools

class Osfd(pipeline_tools.FeatureFinderPipeline):
    def __init__(self):
        self.name = 'Osfd'
        self.cores = 1
        self.memory = 4
        self.timeout = '04:00:00'
        self.get_params()
    
    def get_params(self):
        super().get_params()
        osfd_params = {'mz_step':[0.02, 0.01, 0.03, 0.05, 0.001, 0.005],
                       'SN':[3,2,1.5,1,4,5],
                       'int_threshold':[1,1000,10000],
                       'peakwidth_min':[5,1,10,15,20],
                       'peakwidth_max':[300, 200, 100, 60, 500, 700],
                       'maxPeaksPerSignal':[10,8,5,3,1,15,20],
                       'precursormzTol':[20,15,10,8,5,4],
                       'linear_binning':['TRUE', 'FALSE']}
        self.osfd_param_set = set(osfd_params.keys())

        self.param_choices.update(osfd_params)
        return self.param_choices

    def setup_workspace(self):
        import subprocess
        subprocess.run('singularity build --sandbox --fakeroot osfd.sif docker://stavisvols/osfd',
                       shell = True)
        
    def run_job(self, job):
        super().run_job(job)
        
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
            
            #run OSFD                
            args = ' '.join(f'--{k} {v}' for k,v in job.items() if k in self.osfd_param_set)
            singularity_command = 'singularity run --containall --fakeroot --bind ./:/data/ ../osfd.sif'
            for mzml, base_name in zip(mzmls, base_names):
                osfd_command = 'Rscript /osfd/peakpicking.R {args} -i /data/{mzml} -o /data/{base_name}.features'
                subprocess.run(f'{singularity_command} {osfd_command}',
                               shell = True)
            
            #run peptide rollup
            peptide_results = []
            for base_name in base_names:
                features = pd.read_csv(f'{base_name}.features', sep = '\t')
                
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
            print(e)
        #clean up temporary files
        os.chdir('..')
        shutil.rmtree(tmpdir)

