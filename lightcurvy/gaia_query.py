import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import os
import pandas as pd
import numpy as np

from .mag_to_flux import mag_to_flux

###

def gaia_query(obj = None, ra = None, dec = None,
               search_radius = 0.0015,
               save = True, overwrite = False, obj_save = None, save_dir = None, **kwargs):
    
    print('Grabbing Gaia light curve')
    
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        sae_dir = os.getcwd()
    savename = os.path.join(save_dir, 'gaia', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'gaia'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'Gaia light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        obj_gaia_table = Gaia.query_object(obj, radius = f'{search_radius} deg').to_pandas()
        gaia_id = obj_gaia_table.iloc[0, :]['SOURCE_ID']
        
        print(f'\tObject: {obj}\n\tGaia ID: {gaia_id}\n\tSearch Radius: {search_radius}')
        
        datalink = Gaia.load_data(ids = [gaia_id], 
                                  retrieval_type = 'EPOCH_PHOTOMETRY',
                                  data_structure = 'INDIVIDUAL',
                                  data_release = 'Gaia DR3',
                                  format = 'csv')
        
        if len(datalink) == 0:
            print('\t- - - - - - - - - -')
            print(f'\tFound no Gaia epoch photometry for {obj} (Gaia ID = {gaia_id}).')
            return pd.DataFrame()
        
        else:
    
            lc = next(iter(datalink.values()))[0].to_pandas()
            
            lc['mjd'] = lc['time'] + 57000 - 0.5
            lc['mag_err'] = 1.08574 * (lc['flux_error'] / lc['flux'])
            flag_mask = np.logical_or(lc['rejected_by_photometry'] == 'true', lc['rejected_by_variability'] == 'true')
            lc['flag'] = [1 if f else 0 for f in flag_mask]
            
            # lc.drop(columns = ['source_id', 'transit_id', 'solution_id', 'flux_over_error', 'time', 'flux', 'flux_error', 'rejected_by_photometry', 'rejected_by_variability', 'other_flags'], inplace = True)
            lc.rename(columns = {'band' : 'filter'}, inplace = True)
            lc['source'] = 'gaia'
            lc = mag_to_flux(lc)
            lc.sort_values(by = ['mjd'], inplace = True)
            lc.reset_index(drop = True, inplace = True)
            
            N_filt = len(set(lc['filter']))
            print('\t- - - - - - - - - -')
            print(f'\tFound {len(lc.index)} ZTF photometry points in {N_filt} bands:')
            for filt in set(lc['filter']):
                lc_filt = lc.query('filter == @filt')
                N = len(lc_filt.index)
                fraction_flagged = np.sum(lc_filt['flag']) / N
                print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
                        
            if save:
                lc.to_csv(savename, index = False)
            
            return lc