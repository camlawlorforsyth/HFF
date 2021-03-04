
import glob
import subprocess

from astropy.table import Table

cluster = 'a2744'
outDir = cluster + '/h5/'
filterFile = cluster + '/filters.txt'
photometries = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster, cluster)

# get a list of fits files containing photometric data for all vorbins for
# a given galaxy, as denoted by ID
tables = glob.glob(photometries)

# loop over all the fits files in the directory
for file in tables :
    ID = file.split('_')[2] # the galaxy ID to fit the vorbins for
    table = Table.read(file)
    vorbins = table['vorbin'] # get a list of vorbin values
    for vorbin in vorbins : # loop over all the vorbins in the table
        if table['use'][vorbin] : # complete SED fitting for useable vorbins
            
            # parameters necessary for fitting and writing output
            redshift = table['z'][vorbin]
            lumDist = table['lumDist'][vorbin]
            outfile = '{}{}_ID_{}_BIN_{}'.format(outDir, cluster, ID, vorbin)
            
            # create argument list to pass to params.py
            args = ['python', 'params.py',
                    '--object_redshift', str(redshift),
                    '--luminosity_distance', str(lumDist),
                    '--infile', str(file),
                    '--filterFile', str(filterFile),
                    '--vorbinNum', str(int(vorbin)),
                    '--verbose', str(int(0)), # False
                    '--dynesty',
                    '--outfile', outfile]
            subprocess.run(args) # run the subprocess
