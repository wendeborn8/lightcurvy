import numpy as np

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u

from tess_stars2px import tess_stars2px_function_entry

###

def get_radec(obj):
    results = np.array(Simbad.query_object(obj))
    c = SkyCoord(str(results['RA'][0]) + str(results['DEC'][0]), unit = (u.hourangle, u.deg))
    ra, dec = c.ra.deg, c.dec.deg
    return ra, dec

###

def clean_lc(lc, sigma = 5):
    
    for filt, lc_filt in lc.groupby('filter'):
        for source, lc_source in lc_filt.groupby('source'): 
            med, std = np.nanmedian(lc_source['mag']), np.nanstd(lc_source['mag'])
            good_inds = np.array(lc_source.index)[np.logical_and(np.array(lc_source['mag']) > med - sigma*std, np.array(lc_source['mag']) < med + sigma*std)]
            bad_inds = np.array(lc_source.index)[~np.logical_and(lc_source['mag'] > med - sigma*std, lc_source['mag'] < med + sigma*std)]
            lc.loc[good_inds, 'flag'] = True
            lc.loc[bad_inds, 'flag'] = False
            
    return lc

###

def get_tic_num(obj):
    IDs = Simbad.query_objectids(obj)
    for ID in IDs:
        ID = ID['ID']
        if ID.startswith('TIC '):
            tic_num = int(ID.split()[1])
            break
    return tic_num

###

def find_tess_sectors(obj):
    tic_num = get_tic_num(obj)
    ra, dec = get_radec(obj)
    results = np.array(Simbad.query_object(obj))
    c = SkyCoord(str(results['RA'][0]) + str(results['DEC'][0]), unit = (u.hourangle, u.deg))
    ra, dec = c.ra.deg, c.dec.deg
    _, _, _, sectors, _, _, _, _, _ = tess_stars2px_function_entry(tic_num, ra, dec)
    return sectors