
import os
import glob
import subprocess
import time

from astropy.table import Table

begin = time.time()

clusters = ['m416'] # ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
for cluster in clusters :
    outDir = '{}/h5'.format(cluster)
    filterFile = '{}/{}_filters.txt'.format(cluster, cluster)
    photometries = '{}/photometry/{}_ID_2876_photometry.fits'.format(cluster,
                                                                  cluster)
    
    os.makedirs(outDir, exist_ok=True) # ensure the output directory for the
                                       # results is available
    
    # get a list of fits files containing photometric data for all bins for
    # a given galaxy, as denoted by ID
    phot_files = glob.glob(photometries)
    
    # loop over all the fits files in the directory
    for file in phot_files :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
        table = Table.read(file)
        bins = table['bin'] # get a list of bin values
        
        # create the output directory for the given galaxy
        outGal = '{}/{}'.format(outDir, ID)
        os.makedirs(outGal, exist_ok=True)
        
        for binNum in [10] : #bins : # loop over all the bins in the table
            
            # parameters necessary for fitting and writing output
            redshift = table['z'][binNum]
            # lumDist = table['lumDist'][binNum]
            outfile = '{}/{}_ID_{}_BIN_{}_npTT1'.format(outGal, cluster, ID, binNum)
            
            # create argument list to pass to params.py
            args = ['python', 'params.py', #'mpiexec', '-np', '10', this goes before 'python'
                    '--object_redshift', str(redshift),
                    # '--luminosity_distance', str(lumDist),
                    '--fixed_metallicity', str(0.02),
                    '--infile', str(file),
                    '--filterFile', str(filterFile),
                    '--binNum', str(int(binNum)),
                    '--verbose', str(int(1)), # True
                    '--emcee',
                    # '--dynesty',
                    # '--nested_dlogz_init=500', # see prospect/utils/prospect_args.py
                    # '--nested_posterior_thresh=50', # for more information
                    '--outfile', outfile]
            try :
                subprocess.run(args) # run the subprocess
            except :
                pass

end = time.time()

print(end - begin)
