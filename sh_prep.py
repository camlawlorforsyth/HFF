
import os
import glob

from astropy.table import Table

# clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']

def prep_for_running_on_astro_and_quixote() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    outfile = 'astro1.sh'
    
    with open(outfile, 'a') as slurm :
        slurm.write('#!/usr/bin/env bash\n')
        slurm.write('\n')
    
    for cluster in clusters :
        phot_paths = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                    cluster)
        phots = glob.glob(phot_paths)
        
        for file in phots :
            file = file.replace(os.sep, '/') # compatibility for Windows
            ID = int(file.split('_')[2]) # the galaxy ID to fit the bins for
            
            table = Table.read(file)
            bins = table['bin'] # get a list of bin values
            for binNum in bins : # loop over all the bins in the table
                outLog = '{}/logs/{}_ID_{}_bin_{}.log'.format(cluster, cluster, ID, binNum)
                cmd = ('nohup python params.py' +
                       ' --infile {} '.format(file) +
                       ' --binNum {} >> {} 2>&1\n'.format(binNum, outLog))
                
                with open(outfile, 'a') as slurm :
                    slurm.write(cmd)
    
    return

def integrated_prep(population) : # population == 'Q' | 'SF'
    
    outfile = 'output/sh_runs/integrated_complete_sample_{}.sh'.format(population)
    
    table = Table.read('output/tables/sample_final.fits')
    table = table[table['pop'] == population]
    
    with open(outfile, 'a') as slurm :
        slurm.write('#!/usr/bin/env bash\n')
        slurm.write('\n')
    
    for cluster, ID in zip(table['cluster'], table['id']) :
        file = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster, cluster, ID)
        outLog = '{}/logs/{}_ID_{}.log'.format(cluster, cluster, ID)
        
        cmd = ('nohup python params_integrated_{}.py'.format(population) +
               ' --infile {} >> {} 2>&1\n'.format(file, outLog))
        
        with open(outfile, 'a') as slurm :
            slurm.write(cmd)
    
    return

"""
for cluster in clusters :
    os.makedirs('{}/h5'.format(cluster), # ensure the output directory for the
                exist_ok=True)           # results is available
    
    os.makedirs('{}/logs'.format(cluster), # ensure the output directory for
                exist_ok=True)             # log files is available
    
    os.makedirs('{}/sh'.format(cluster), # ensure the run directory for the
                exist_ok=True)           # slurm files is available
    
    # get a list of fits files containing photometric data for all bins for
    # a given galaxy, as denoted by ID
    photometries = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                  cluster)
    phot_files = glob.glob(photometries)
    
    phots = []
    for file in phot_files :
        file = file.replace(os.sep, '/') # compatibility for Windows
        phots.append(file)
    
    # loop over all the fits files in the directory
    for file in phots :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
        
        table = Table.read(file)
        bins = table['bin'] # get a list of bin values
        
        for binNum in bins : # loop over all the bins in the table
            
            outfile = '{}/h5/ID_{}_bin_{}'.format(cluster, ID, binNum)
            outLog = '{}/logs/ID_{}_bin_{}.log'.format(cluster, ID, binNum)
            errLog = '{}/logs/ID_{}_bin_{}_err.log'.format(cluster, ID, binNum)
            
            cmd = ('python params.py' +
                   ' --infile {} --outfile {}'.format(file, outfile) +
                   ' --binNum {} --dynesty'.format(binNum) +
                   ' >> {} 2>> {}'.format(outLog, errLog))
            
            outRun = '{}/sh/ID_{}_bin_{}_run.sh'.format(cluster, ID, binNum)
            
            with open(outRun, 'a') as slurm :
                slurm.write(cmd)
            
            '''
            with open(outRun, 'a') as slurm :
                slurm.write('#!/bin/bash\n')
                slurm.write('#SBATCH --account=clawlorf\n')
                slurm.write('#SBATCH --time=03:00:00\n')
                slurm.write('#SBATCH --job-name={}_{}_{}\n'.format(cluster, ID, binNum))
                slurm.write('#SBATCH --output={}\n'.format(outLog))
                slurm.write('#SBATCH --error={}\n'.format(errLog))
                slurm.write('source ../env/bin/activate\n')
                slurm.write('module load python/3.8.10\n')
                slurm.write(cmd)
            '''
"""
