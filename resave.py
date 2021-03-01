
from astropy.table import Table

def save_as_fits(file, outfile, rename_cols=False) :
    
    table = Table.read(file, format='ascii')
    
    if rename_cols :
        names = tuple(table.colnames)
        new_names = ('id','z', 'ltau', 'metal', 'lage', 'Av',
                     'lmass', 'lsfr', 'lssfr', 'la2t', 'chi2')
        table.rename_columns(names, new_names)
    
    table.write(outfile)
    
    return

'''
# PHOTOMETRY - CLUSTERS
save_as_fits('catalogs/abell370clu_catalogs/hffds_abell370clu_v3.9.cat',
             'catalogs/abell_370_cluster_photometry.fits')
save_as_fits('catalogs/abell1063clu_catalogs/hffds_abell1063clu_v3.9.cat',
             'catalogs/abell_S1063_cluster_photometry.fits')
save_as_fits('catalogs/abell2744clu_catalogs/hffds_abell2744clu_v3.9.cat',
             'catalogs/abell_2744_cluster_photometry.fits')
save_as_fits('catalogs/macs0416clu_catalogs/hffds_macs0416clu_v3.9.cat',
             'catalogs/macs_0416_cluster_photometry.fits')
save_as_fits('catalogs/macs0717clu_catalogs/hffds_macs0717clu_v3.9.cat',
             'catalogs/macs_0717_cluster_photometry.fits')
save_as_fits('catalogs/macs1149clu_catalogs/hffds_macs1149clu_v3.9.cat',
             'catalogs/macs_1149_cluster_photometry.fits')

# PHOTOMETRY - PARALLELS
save_as_fits('catalogs/abell370par_catalogs/hffds_abell370par_v3.9.cat',
             'catalogs/abell_370_parallel_photometry.fits')
save_as_fits('catalogs/abell1063par_catalogs/hffds_abell1063par_v3.9.cat',
             'catalogs/abell_S1063_parallel_photometry.fits')
save_as_fits('catalogs/abell2744par_catalogs/hffds_abell2744par_v3.9.cat',
             'catalogs/abell_2744_parallel_photometry.fits')
save_as_fits('catalogs/macs0416par_catalogs/hffds_macs0416par_v3.9.cat',
             'catalogs/macs_0416_parallel_photometry.fits')
save_as_fits('catalogs/macs0717par_catalogs/hffds_macs0717par_v3.9.cat',
             'catalogs/macs_0717_parallel_photometry.fits')
save_as_fits('catalogs/macs1149par_catalogs/hffds_macs1149par_v3.9.cat',
             'catalogs/macs_1149_parallel_photometry.fits')

# FAST OUTPUT - CLUSTERS
save_as_fits('catalogs/abell370clu_catalogs/fast/abell370clu_v3.9/abell370clu_v3.9.fout',
             'catalogs/abell_370_cluster_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/abell1063clu_catalogs/fast/abell1063clu_v3.9/abell1063clu_v3.9.fout',
             'catalogs/abell_S1063_cluster_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/abell2744clu_catalogs/fast/abell2744clu_v3.9/abell2744clu_v3.9.fout',
             'catalogs/abell_2744_cluster_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs0416clu_catalogs/fast/macs0416clu_v3.9/macs0416clu_v3.9.fout',
             'catalogs/macs_0416_cluster_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs0717clu_catalogs/fast/macs0717clu_v3.9/macs0717clu_v3.9.fout',
             'catalogs/macs_0717_cluster_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs1149clu_catalogs/fast/macs1149clu_v3.9/macs1149clu_v3.9.fout',
             'catalogs/macs_1149_cluster_fastoutput.fits', rename_cols=True)

# FAST OUTPUT - PARALLELS
save_as_fits('catalogs/abell370par_catalogs/fast/abell370par_v3.9/abell370par_v3.9.fout',
             'catalogs/abell_370_parallel_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/abell1063par_catalogs/fast/abell1063par_v3.9/abell1063par_v3.9.fout',
             'catalogs/abell_S1063_parallel_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/abell2744par_catalogs/fast/abell2744par_v3.9/abell2744par_v3.9.fout',
             'catalogs/abell_2744_parallel_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs0416par_catalogs/fast/macs0416par_v3.9/macs0416par_v3.9.fout',
             'catalogs/macs_0416_parallel_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs0717par_catalogs/fast/macs0717par_v3.9/macs0717par_v3.9.fout',
             'catalogs/macs_0717_parallel_fastoutput.fits', rename_cols=True)
save_as_fits('catalogs/macs1149par_catalogs/fast/macs1149par_v3.9/macs1149par_v3.9.fout',
             'catalogs/macs_1149_parallel_fastoutput.fits', rename_cols=True)
'''
