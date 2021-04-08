
import os
import numpy as np

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
        numBins = int(np.ceil(1.3*np.sqrt(len(lengths))))
        plt.histogram(lengths, 'Number of Annuli per Galaxy', title=cluster,
                      bins=numBins)
    
    return
