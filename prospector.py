
import os
import glob
import subprocess

from astropy.table import Table

clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
for cluster in clusters :
    outDir = cluster + '/' + 'h5'
    filterFile = cluster + '/' + 'filters.txt'
    photometries = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                  cluster)
    
    os.makedirs(outDir, exist_ok=True) # ensure the output directory for the
                                       # results is available
    
    # get a list of fits files containing photometric data for all bins for
    # a given galaxy, as denoted by ID
    tables = glob.glob(photometries)
    
    # loop over all the fits files in the directory
    for file in tables :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
        table = Table.read(file)
        bins = table['bin'] # get a list of bin values
        
        # create the output directory for the given galaxy
        outGal = '{}/{}'.format(outDir, ID)
        os.makedirs(outGal, exist_ok=True)
        
        for binNum in bins : # loop over all the bins in the table
            
            # parameters necessary for fitting and writing output
            redshift = table['z'][binNum]
            lumDist = table['lumDist'][binNum]
            outfile = '{}/{}_ID_{}_BIN_{}'.format(outGal, cluster, ID, binNum)
            
            # create argument list to pass to params.py
            args = ['python', 'params.py',
                    '--object_redshift', str(redshift),
                    '--luminosity_distance', str(lumDist),
                    '--infile', str(file),
                    '--filterFile', str(filterFile),
                    '--binNum', str(int(binNum)),
                    '--verbose', str(int(0)), # False
                    '--dynesty',
                    '--outfile', outfile]
            subprocess.run(args) # run the subprocess
