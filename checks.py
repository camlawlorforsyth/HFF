
import os

from astropy.table import Table

import plotting as plt

def check_bins(cluster, plot=True) :
    
    photDir = '{}/photometry'.format(cluster)
    files = os.listdir(photDir)
    
    lengths = []
    for file in files :
        table = Table.read('{}/{}'.format(photDir, file))
        length = len(table)
        lengths.append(length)
    
    if plot :
        plt.histogram(lengths, 'Number of Annuli per Galaxy', title=cluster,
                      bins=19)
    
    return
