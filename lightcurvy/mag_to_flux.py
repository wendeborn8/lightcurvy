import pandas as pd
import numpy as np
from tqdm import tqdm
import os

###

def mag_to_flux(df, verbose = True):
    
    filter_infos_file = os.path.join(os.path.dirname(__file__), r'filter_infos.csv')
    filter_infos = pd.read_csv(filter_infos_file, index_col = 0)
    
    filters = np.array(df['filter'])
    filters_set = set(filters)
    mags = np.array(df['mag']).astype(float)
    if 'mag_err' in df.columns:
        err = True
        mag_errs = np.array(df['mag_err']).astype(float)
    else:
        if verbose:
            print(r'\tColumn \'mag_err\' not found. Not calculating flux uncertainties.')
        err = False
        
    if verbose:
        bad_filters = [filt for filt in filters_set if filt not in list(filter_infos.index)]
        print('\t- - - - - - - - - -')
        print(f'\tConverting {len(mags)} magnitudes to fluxes in {len(filters_set)} filters: {filters_set}.')
        if len(bad_filters) > 0:
            print(f'\t\t{bad_filters} have no entry in Filter Info file, so no fluxes will be calculated.')
    
    wls_eff, systems, zps_Jy, zps_Flam = [], [], [], []
    for filt in tqdm(filters, total = len(filters), desc = '\tConverting magnitudes to fluxes'):
        if filt in set(list(filter_infos.index)):
            wl_eff, system = filter_infos.loc[filt, 'wl_eff'], filter_infos.loc[filt, 'system']
            wls_eff.append(wl_eff)
            systems.append(system)
            
            zps_Jy.append(filter_infos.loc[filt, f'zp_{system}_Jy'])
            zps_Flam.append(filter_infos.loc[filt, f'zp_{system}_erg'])
            
        else:
            wls_eff.append(np.nan)
            systems.append(np.nan)
            zps_Jy.append(np.nan)
            zps_Flam.append(np.nan)
            
    wls_eff, zps_Jy, zps_Flam = np.array(wls_eff), np.array(zps_Jy), np.array(zps_Flam)
    
    # Convert mags to fluxes using zero points
    df['flux_Jy'] = zps_Jy * 10**(-mags / 2.5)
    df['flux_Flam'] = zps_Flam * 10**(-mags / 2.5)
    df['flux_lamFlam'] = wls_eff * df['flux_Flam']  
    
    if err:
        df['err_Jy'] = 0.921034 * zps_Jy * mag_errs * np.sqrt(np.exp(-1.84204 * mags))
        df['err_Flam'] = 0.921034 * zps_Flam * mag_errs * np.sqrt(np.exp(-1.84204 * mags))
        df['err_lamFlam'] = wls_eff * df['err_Flam']
        
    return df
    