
import numpy as np

from astropy.io import fits
from astropy.table import Table

def add_gaussian_background(infile, outfile, rms) :
    
    with fits.open(infile) as hdu :
        hdr = hdu[0].header
        data = hdu[0].data
    
    background = np.random.normal(loc=0.0, scale=rms, size=data.shape)
    
    hdu = fits.PrimaryHDU(data + background)
    hdu.header = hdr
    hdu.writeto(outfile)
    
    return

def replace(cluster, old_ids, new_ids, print_ids=False) :
    
    if print_ids :
        temp = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
        temp = temp[temp['id'] > 20000]
        print(list(temp['id']))
    
    old_segmap = 'sextractor/output/{}_bcg_segmap.fits.gz'.format(cluster)
    new_segmap = 'sextractor/{}_model_segmap.fits.gz'.format(cluster)
    
    with fits.open(old_segmap) as hdu :
        hdr = hdu[0].header
        data = hdu[0].data
    
    if len(old_ids) == len(new_ids) :
        for old_id, new_id in zip(old_ids, new_ids) :
            data[data == old_id] = new_id
        
        hdu = fits.PrimaryHDU(data)
        hdu.header = hdr
        hdu.writeto(new_segmap)
    else :
        print('The number of elements is not the same.')
    
    return

def replace_all() :
    
    replace('a370',
            [ 4, 32, 43, 44, 18, 14, 37, 42, 24, 38,
              6, 12, 10, 27, 45, 39, 41, 50, 22, 11,
             17, 25, 13, 19],
            [20001, 20003, 20004, 20007, 20020, 20021, 20023, 20024, 20026, 20027,
             20031, 20034, 20035, 20042, 20043, 20044, 20045, 20046, 20051, 20053,
             20055, 20058, 20061, 20062])
    
    replace('a1063',
            [26, 30, 34, 36,  8, 11, 40, 38, 28, 31,
             14, 17, 15, 16, 19, 22, 23,  6, 39, 44,
             37,  4, 41, 43, 50, 53, 46, 48, 33, 56,
             51, 32, 10,  9, 35, 20, 29, 27, 20060], # 20060 not detected in segmap
            [20002, 20003, 20005, 20008, 20010, 20012, 20014, 20015, 20016, 20017,
             20018, 20019, 20020, 20022, 20023, 20024, 20025, 20028, 20029, 20030,
             20031, 20033, 20035, 20036, 20040, 20042, 20044, 20045, 20046, 20049,
             20050, 20053, 20054, 20055, 20056, 20057, 20058, 20059, 20060])
    
    replace('a2744',
            [11,  6, 14, 22, 17, 27, 31, 36, 19, 29,
             33, 30, 40, 41,  7, 16, 42, 47, 49, 38,
             45, 25, 28, 23,  9, 39, 37, 21, 26, 32,
             12, 13, 15, 18, 34, 24, 54, 55, 53, 58,
             59, 48,  8,  3,  1,  4, 56, 64, 62, 61,
             60, 10,  5, 51, 43, 50, 52, 44, 20, 35],
            [20001, 20002, 20003, 20004, 20005, 20006, 20007, 20008, 20009, 20010,
             20011, 20012, 20013, 20014, 20015, 20016, 20017, 20018, 20019, 20020,
             20021, 20022, 20023, 20024, 20025, 20026, 20027, 20028, 20029, 20030,
             20031, 20032, 20033, 20034, 20035, 20036, 20037, 20038, 20039, 20040,
             20041, 20042, 20043, 20044, 20045, 20046, 20048, 20049, 20050, 20051,
             20053, 20054, 20055, 20056, 20057, 20058, 20059, 20061, 20062, 20063])
    
    replace('m416',
            [11,  6,  5, 23, 12, 10,  9, 25, 31, 34,
              7, 20, 22, 19, 21, 18,  8, 16, 36, 35,
             38, 24, 26, 27, 29, 30, 33, 37, 39, 28,
             32, 41, 42, 43, 14, 13, 15, 17,  2, 40,
              4, 3],
            [20001, 20002, 20003, 20005, 20006, 20007, 20008, 20009, 20010, 20011,
             20012, 20013, 20014, 20015, 20016, 20017, 20018, 20019, 20020, 20021,
             20022, 20023, 20024, 20025, 20026, 20027, 20028, 20029, 20030, 20031,
             20032, 20033, 20034, 20035, 20036, 20037, 20038, 20039, 20040, 20043,
             20045, 20046])
    
    replace('m717',
            [ 9,  1, 19, 22, 18, 13, 11, 14, 27, 25,
             29, 32, 21, 23, 10, 31,  6,  8,  3,  4,
             24],
            [20004, 20005, 20006, 20007, 20008, 20009, 20010, 20011, 20015, 20016,
             20017, 20018, 20020, 20021, 20022, 20023, 20027, 20028, 20029, 20030,
             20031])
    
    replace('m1149',
            [17, 23, 21, 11,  6,  3,  8,  9, 16, 44,
             48, 30, 37, 10, 51, 49, 52, 43, 47, 39,
             40, 27, 18, 28, 26, 38, 22, 35, 14,  5,
              4, 13, 50, 45, 41, 54, 33, 19, 31, 29,
             24, 34,  1, 32, 25, 36],
            [20002, 20003, 20004, 20005, 20006, 20007, 20008, 20009, 20011, 20012,
             20014, 20016, 20017, 20019, 20020, 20021, 20022, 20023, 20024, 20025,
             20026, 20027, 20028, 20029, 20030, 20031, 20032, 20033, 20034, 20035,
             20036, 20038, 20040, 20041, 20042, 20043, 20044, 20045, 20046, 20047,
             20048, 20049, 20050, 20052, 20053, 20054])
    
    replace('a370par',
            [5, 1],
            [20004, 20005])
    
    replace('a1063par',
            [2, 3, 4, 6, 13],
            [20005, 20006, 20010, 20011, 20014])
    
    replace('a2744par',
            [ 5,  7,  9, 11, 10,  2,  3, 14, 13, 15,
             18, 12,  6,  4, 16],
            [20001, 20002, 20003, 20004, 20006, 20007, 20008, 20009, 20010, 20012,
             20015, 20016, 20017, 20018, 20019])
    
    replace('m416par',
            [5, 6, 2, 1, 3, 7],
            [20001, 20002, 20004, 20005, 20006, 20007])
    
    replace('m717par',
            [4, 2, 1, 3, 5],
            [20001, 20002, 20003, 20004, 20007])
    
    replace('m1149par',
            [5, 4, 6],
            [20001, 20002, 20006])
    
    return

def sextractor_prep() :    
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149',
                'a370par', 'a1063par', 'a2744par', 'm416par', 'm717par', 'm1149par']
    
    models = ['misc/abell370clu_misc/images/bcgs_models/abell370_f160w_bcgs_model.fits.gz',
              'misc/abell1063clu_misc/images/bcgs_models/abell1063_f160w_bcgs_model.fits.gz',
              'misc/abell2744clu_misc/images/bcgs_models/abell2744_f160w_bcgs_model.fits.gz',
              'misc/macs0416clu_misc/images/bcgs_models/macs0416_f160w_bcgs_model.fits.gz',
              'misc/macs0717clu_misc/images/bcgs_models/macs0717_f160w_bcgs_model.fits.gz',
              'misc/macs1149clu_misc/images/bcgs_models/macs1149_f160w_bcgs_model.fits.gz',
              'misc/abell370par_misc/images/bcgs_models/abell370par_f160w_bcgs_model.fits.gz',
              'misc/abell1063par_misc/images/bcgs_models/abell1063par_f160w_bcgs_model.fits.gz',
              'misc/abell2744par_misc/images/bcgs_models/abell2744par_f160w_bcgs_model.fits.gz',
              'misc/macs0416par_misc/images/bcgs_models/macs0416par_f160w_bcgs_model.fits.gz',
              'misc/macs0717par_misc/images/bcgs_models/macs0717par_f160w_bcgs_model.fits.gz',
              'misc/macs1149par_misc/images/bcgs_models/macs1149par_f160w_bcgs_model.fits.gz']
    
    RMSes = [0.008631396, 0.023791686, 0.010081889, 0.0073933504, 0.005936407, 0.008023425,
             0.0042280676, 0.0049573397, 0.0036131504, 0.012034563, 0.005247639, 0.0035911163]
    
    for cluster, model, RMS in zip(clusters, models, RMSes) :
        outfile = 'sextractor/images/{}_f160w_model_sextractor.fits.gz'.format(cluster)
        add_gaussian_background(model, outfile, RMS)
    
    return
