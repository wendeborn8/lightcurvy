import os
import pandas as pd
import numpy as np
from sklearn import linear_model

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import conf as mast_config
from astropy.io import fits

from tglc.quick_lc import tglc_lc
import lightkurve as lk
import eleanor

from . import util
from .mag_to_flux import mag_to_flux

###

def tglc_query(tic_id, obj, root, sectors = None, ffi_cutout_size = 50, limit_mag = 16, timeout = 6000, prior = 0.5, flux_col = 'flux_psf', err_col = 'flux_psf_err'):
    
    if sectors is None:
        sectors = util.find_tess_sectors(obj = obj)
        
    if tic_id is None:
        tic_num = util.get_tic_num(obj = obj)
        tic_id = f'TIC {tic_num}'
        
    with mast_config.set_temp('timeout', timeout):
        for sector in sectors:
            try:
                tglc_lc(target = tic_id, 
                        local_directory = root, 
                        size = ffi_cutout_size,
                        save_aper = True,
                        limit_mag = limit_mag,
                        get_all_lc = False, 
                        first_sector_only = False,
                        sector = int(sector),
                        prior = 0.9)
            except Exception as e:
                print(e)
            
    lc_dir = os.path.join(root, 'lc')
    files = [os.path.join(lc_dir, file) for file in os.listdir(lc_dir) if file.endswith('.fits') and 'spoc' not in file and 'eleanor' not in file]
    lcs = []
    for file in files:
        
        try:         
            hdus = fits.open(file)
            data = hdus[1].data
            header = hdus[1].header
            sector = hdus[0].header['SECTOR']
        
            lc = pd.DataFrame()
            for hdr_key, my_key in zip(['time', 'psf_flux', 'PSF_ERR', 'aperture_flux', 'APER_ERR', 'TESS_flags', 'TGLC_flags'], ['time', 'flux_psf', 'flux_psf_err', 'flux_aper', 'flux_aper_err', 'flagged_TESS', 'flagged_TGLC']):
                if hdr_key in ['PSF_ERR', 'APER_ERR']:
                    values = [header[hdr_key] for _ in data['time']]
                else:
                    values = data[hdr_key]
                lc[my_key] = values
                
            lc['mjd'] = lc['time'] + 57000.0 - 0.5
            lc['sector'] = sector
                          
            lcs.append(lc)
            hdus.close()
            
        except:
            continue
        
    # Combine the light curves for each sector into one and save
    lc = pd.concat(lcs, axis = 0)
    lc['mag'] = -2.5 * np.log10(lc[flux_col]) + 20.44
    lc['mag_err'] = 1.08574 * (lc[err_col] / lc[flux_col])
    lc['filter'] = 'TESS'
    lc['source'] = 'TESS_tglc'
    lc['ffi_cutout_size'] = ffi_cutout_size
    lc['limit_mag'] = limit_mag
    lc['prior'] = prior
    lc['flag'] = [1 if (flag_tess or flag_tglc) else 0 for flag_tess, flag_tglc in zip(lc['flagged_TESS'], lc['flagged_TGLC'])]
    lc = mag_to_flux(lc)
    lc.dropna(subset = ['mjd', flux_col], inplace = True)
    lc.sort_values(by = 'mjd', inplace = True)
    lc.reset_index(drop = True, inplace = True)
    
    return lc

###

def spoc_query(obj, obj_save, download_dir, tic_id = None, search_radius = 0.0015, flux_col = 'pdcsap_flux', err_col = 'pdcsap_flux_err'):
    
    print(f'\tFetching SPOC TESS light curves for {obj} and saving to {download_dir}...')
    os.makedirs(download_dir, exist_ok = True)
    spoc_dir = os.path.join(download_dir, 'mastDownload', 'TESS')
    os.makedirs(spoc_dir, exist_ok = True)
    
    # Downlaod all the available SPOC light curves
    if tic_id is None:
        search_results = lk.search_lightcurve(obj, radius = search_radius)
    else:
        search_results = lk.search_lightcurve(tic_id)
    df = search_results.table.to_pandas()
        
    
    if len(search_results) == 0:
        print(f'\t\tFound no SPOC light curves for {obj}...')
        lc = pd.DataFrame()
        lc.to_csv(os.path.join(download_dir, f'{obj_save}_spoc.csv'), index = False)
        return
        
    df.query('provenance_name == \'SPOC\'', inplace = True)
    df.query('not obs_id.str.endswith(\'fast\')', inplace = True)
    if tic_id is not None:
        tic_num = str(tic_id.split()[-1])
        df.query('target_name == @tic_num', inplace = True)
    inds = df.index
    search_results = search_results[list(inds)]
    if len(search_results) == 0:
        print(f'\t\tFound no SPOC light curves for {obj}...')
        lc = pd.DataFrame()
        lc.to_csv(os.path.join(download_dir, f'{obj_save}_spoc.csv'), index = False)
        return
    
    else:
        search_results.download_all(download_dir = download_dir, quality_bitmask = 'none')
    
    # Convert those FITS light curves files to CSV files
    files = []
    for dirr in os.listdir(spoc_dir):
        for file in os.listdir(os.path.join(spoc_dir, dirr)):
            files.append(os.path.join(spoc_dir, dirr, file))

    lcs, sectors = [], []
    for file in files:
        basename = os.path.basename(file)
        sector = int(basename.split('-')[1].lstrip('s'))
        df = lk.io.tess.read_tess_lightcurve(file).to_pandas()
        df['sector'] = sector
        sectors.append(sector)
        lcs.append(df)
        
    print(f'\t\tFound SPOC light curves for sectors {sectors}...')
        
    lc = pd.concat(lcs)
    lc.index.names = ['index']
    lc['time'] = list(lc.index)
    lc['flag'] = [0 if qual == 0 else 1 for qual in lc['quality']]
    lc['mjd'] = lc['time'] + 57000 - 0.5  
    lc['mag'] = -2.5 * np.log10(lc[flux_col]) + 20.44
    lc['mag_err'] = 1.08574 * (lc[err_col] / lc[flux_col])
    lc['filter'] = 'TESS'
    lc['source'] = 'TESS_spoc'
    lc['search_radius'] = search_radius
    lc = mag_to_flux(lc)
    lc.sort_values(by = 'mjd', inplace = True)
    lc.reset_index(inplace = True, drop = True)
          
    return lc
    
###

def eleanor_query(obj, download_dir, flux_col = 'flux_pca', sectors = None, tic_num = None,
                  height = 15, width = 15, bkg_size = 31, do_psf = False, do_pca = True, regressors = 'corner'):
    
    if sectors is None:
        sectors = util.find_tess_sectors(obj = obj)
    if tic_num is None:
        tic_num = util.get_tic_num(obj = obj)
    
    lcs = []
    sectors = sectors[sectors < 75]
    star = eleanor.multi_sectors(tic = tic_num, sectors = list(sectors))
    for s, sector in zip(star, sectors):
        try:
            data = eleanor.TargetData(s, height = height, width = width, bkg_size = bkg_size, do_psf = do_psf, do_pca = do_pca, regressors = regressors)
        except Exception as e:
            print(f'\t{e}')
            print(f'Issue retrieving eleanor light curve for Sector {sector}')
                
        lc = pd.DataFrame({'mjd' : data.time + 57000 - 0.5, 
                           'flux_raw' : data.raw_flux, 'flux_corr' : data.corr_flux, 'flux_pca' : data.pca_flux, 
                           'flux_err' : data.flux_err, 'quality' : data.quality})
        
        lc['flag'] = [0 if qual == 0 else 1 for qual in lc['quality']]
        lc['mag'] = -2.5 * np.log10(lc[flux_col]) + 20.44
        lc['mag_err'] = 1.08574 * (lc['flux_err'] / lc[flux_col])
        lc['source'] = 'TESS_eleanor'
        lc['filter'] = 'TESS'
        lc['sector'] = sector
        
        lcs.append(lc)
        
    lc = pd.concat(lcs)
    lc = mag_to_flux(lc)
    lc.sort_values(by = 'mjd', inplace = True)
    lc.reset_index(inplace = True, drop = True)
    
    return lc

###
    
def tess_query(obj = None, ra = None, dec = None,
               save = True, overwrite = False, obj_save = None, save_dir = None,
               sectors = None, ffi_cutout_size = 50, limit_mag = 15, timeout = 6000,
               prior = 0.5, flux_col = 'flux_psf', err_col = 'flux_psf_err',
               eleanor_size = 15, eleanor_bkg_size = 31,
               tess_sources = ['tglc', 'spoc', 'eleanor'],
               **kwargs):
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        root = os.path.join(os.getcwd(), 'TESS', obj_save) + os.path.sep
    else:
        root = os.path.join(save_dir, 'TESS', obj_save) + os.path.sep
    tglc_savename = os.path.join(root, rf'{obj_save}_tglc.csv')
    spoc_savename = os.path.join(root, f'{obj_save}_spoc.csv')
    eleanor_savename = os.path.join(root, f'{obj_save}_eleanor.csv')
    savename = os.path.join(root, f'{obj_save}.csv')
    os.makedirs(root, exist_ok = True)
    
    # Find the sectors, using TESS-point
    sectors = util.find_tess_sectors(obj = obj)
    
    # Find the TIC ID/Num
    tic_num = util.get_tic_num(obj = obj)
    tic_id = f'TIC {tic_num}'

    print(f'Fetching TESS light curves for:')
    print(f'\tObject: {obj}')
    print(f'\tObject Savename: {obj_save}')
    print(f'\tTIC Number: {tic_num}')
    print(f'\tSectors: {sectors}')
    print(f'\tSaving results to: {root}')
    print('\t' + 10 * '- ')
    
    
    lcs = []
    
    # Get SPOC Light Curve(s)
    if overwrite or not os.path.exists(spoc_savename):
        print(f'\tFetching TESS SPOC light curves for sectors {sectors}...')
        spoc_lc = spoc_query(obj = obj, obj_save = obj_save, download_dir = root, tic_id = tic_id, )
        if spoc_lc is not None:
            lcs.append(spoc_lc)
            if save:
                spoc_lc.to_csv(spoc_savename, index = False)
            print(f'\n\t\t\tSuccessfully fetched TESS SPOC light curve(s) for Sector(s) {sectors}.')
        else:
            spoc_lc = pd.DataFrame()
            if save:
                spoc_lc.to_csv(spoc_savename, index = False)
    else:
        print(f'\tTESS SPOC light curve for {obj} already exists at {spoc_savename}. Will not overwrite.')
        lcs.append(pd.read_csv(spoc_savename))
    print('\t' + 10 * '- ')
        
        
    # Get Eleanor Light Curve(s)
    if overwrite or not os.path.exists(eleanor_savename):
        print(f'\tFetching TESS eleanor light curves for sectors {sectors}...')
        eleanor_lc = eleanor_query(obj = obj, download_dir = root, flux_col = 'flux_pca', sectors = sectors, tic_num = None, height = eleanor_size, width = eleanor_size, bkg_size = eleanor_bkg_size, do_psf = False, do_pca = True, regressors = 'corner')
        if eleanor_lc is not None:
            lcs.append(eleanor_lc)
            if save:
                eleanor_lc.to_csv(eleanor_savename, index = False)
            print(f'\n\t\t\tSuccessfully fetched TESS tglc light curve(s) for Sector(s) {sectors}.')
        else:
            eleanor_lc = pd.DataFrame()
            if save:
                eleanor_lc.to_csv(eleanor_savename, index = False)
    else:
        print(f'\tTESS eleanor light curve for {obj} already exists at {eleanor_savename}. Will not overwrite.')
        lcs.append(pd.read_csv(eleanor_savename))
    print('\t' + 10 * '- ')
        
        
    # Get tglc Light Curve(s)
    if overwrite or not os.path.exists(tglc_savename):
        print(f'\tFetching TESS tglc light curves for Sectors {sectors}...')
        tglc_lc = tglc_query(tic_id = tic_id, obj = obj, root = root, sectors = sectors, ffi_cutout_size = ffi_cutout_size, limit_mag = limit_mag, timeout = timeout, prior = prior, flux_col = flux_col, err_col = err_col)
        if tglc_lc is not None:
            lcs.append(tglc_lc)
            if save:
                tglc_lc.to_csv(tglc_savename, index = False)
            print(f'\n\t\t\tSuccessfully fetched TESS tglc light curve(s) for Sector(s) {sectors}.')
        else:
            tglc_lc = pd.DataFrame()
            if save:
                tglc_lc.to_csv(tglc_savename, index = False)
    else:
        print(f'\tTESS tglc light curve for {obj} already exists at {tglc_savename}. Will not overwrite.')
        lcs.append(pd.read_csv(tglc_savename))
    
    
    lc = pd.concat(lcs)
    lc.sort_values(by = 'mjd', inplace = True)
    lc.reset_index(drop = True, inplace = True)
    if save:
        lc.to_csv(savename, index = False)
                
    return lc


