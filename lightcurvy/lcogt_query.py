import io
import re
import sys
import time
import requests
import os
import pandas as pd
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
from tqdm import tqdm
import json
import copy
from astropy.io import fits

from . import util
from .mag_to_flux import mag_to_flux

###

def get_token(username = None, password = None, token = None):
    
    if (username is None or password is None) and token is None:
        raise NameError('Either username AND password, or token must be given')
    elif (username is None or password is None) and token is not None:
        token = token
    elif (username is not None and password is not None) and token is None:
        try:
            response = requests.post('https://observe.lco.global/api/api-token-auth/',
                                     data = {'username': username,
                                             'password': password}).json()
            token = response['token']
        except:
            print('Hmm. That didnt work. Username/password probably not accepted. Proceeding as anonymous user.')
    elif username is not None and password is not None and token is not None:
        try:
            response = requests.post('https://observe.lco.global/api/api-token-auth/',
                                     data = {'username': username,
                                             'password': password}).json()
            token_retrieved = response['token']
            if token != token_retrieved:
                print('Retrieved token from username/password does not match given token. Using the token retrieved using the given username/password.')
                return token_retrieved
        except:
            print('Username, password, AND token all given, but username/password not accepted. Using the given token only.')

    return token

def get_url(proposal_id,
            start,
            end,
            reduction_level,
            primary_optical_element,
            telescope_id,
            target_name,
            covers,
            limit = 1000):
    
    url = 'https://archive-api.lco.global/frames/?'\
        f'proposal_id={proposal_id}&'\
        f'start={start}&'\
        f'end={end}&'\
        f'covers={covers}&'\
        'public=true&'\
        f'limit={limit}&'\
        f'target_name={target_name}&'\
        f'reduction_level={reduction_level}&'\
        f'primary_optical_element={primary_optical_element}&'\
        f'telescope_id={telescope_id}&'\
        f'exclude_calibrations=true'
    
    return url

def get_frames(url, username, password, token, verbose = False):
    
    print(f'Fetching frames for the url: {url}')
    
    token = get_token(username, password, token)
    if token is None:
        if (username is not None or password is not None):
            print('Token could not be retrieved. Proceeding as anonymous user.')
        response = requests.get(url)
    else:
        if verbose:
            print(f'Token retrieved for user {username}: {token}.')
        response = requests.get(url, headers = {'Authorization' : f'token {token}'})
    
    response_json = json.loads(response.text)
    frames = response_json['results']
    if response_json['next']:
        more = True
        response_new = copy.deepcopy(response_json)
        while more:
            response_new = requests.get(response_new['next'], headers = {'Authorization' : f'token {token}'}).json()
            frames += response_new['results']
            if response_new['next']:
                more = True
            else:
                more = False
    
    return frames

###

def lcogt_query(obj = None, ra = None, dec = None,
                save = True, overwrite = False, overwrite_images = False, obj_save = None, image_savedir = None,
                lcogt_username = 'wendebo2@bu.edu', lcogt_password = 'Qzectbum13579@LCO', lcogt_token = '',
                start = '2014-01-01', end = '', verbose = False,
                telescopes = 'all', reduction_level = 91, limit = 1000, filters = 'all', proposal_id = '', **kwargs
                ):
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if image_savedir is None:
        image_savedir = os.path.join(r'D:\My Drive\Data', 'LCOGT', obj_save)
    savename = os.path.join(r'D:\My Drive\Data', 'LCOGT', rf'{obj_save}.csv')
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'LCOGT light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:

        covers = f'POINT%28{ra}%20{dec}%29'
        target_name = ''
            
        if isinstance(filters, str):
            if filters in ['all', '', '[all]', '[\'all\']']:
                filters = ['']
            else:
                if ',' in filters:
                    filters = filters.split(',')
                else:
                    filters = [filters]
        
        if telescopes in ['all', '']:
            telescopes = ['']
        else:
            if ',' in telescopes:
                telescopes = telescopes.split(',')
            else:
                telescopes = [telescopes]
        
        # If the output directory (outdir) doesn't exist, make it
        os.makedirs(image_savedir, exist_ok = True)
        
        print(f'Grabbing LCOGT light curves\n\tObject: {obj}\n\tRA: {ra:.4f}\n\tDec: {dec:.4f}\n\tStart/End Dates: {start}/{end}\n\tTelescopes: {telescopes}\n\tFilters: {filters}\n\tLCOGT Username / Password / Token: {lcogt_username} / {lcogt_password} / {lcogt_token}')
    
        
        # Collect all the valid frames
        all_frames = []
        for filt in filters:
            for telescope_id in telescopes:
                request_url = get_url(proposal_id = proposal_id, start = start, end = end, 
                                      reduction_level = reduction_level, primary_optical_element = filt,
                                      target_name = target_name, telescope_id = telescope_id,
                                      limit = limit, covers = covers)
                
                frames = get_frames(request_url, username = lcogt_username, password = lcogt_password, token = lcogt_token, verbose = verbose)
                all_frames += frames
        
        # Try to download/write each frame
        existing_files, new_files, failed_files = 0, 0, 0
        print('\t- - - - - - - - - -')
        for frame in tqdm(all_frames, desc = f'\tFound {len(all_frames)} frames. Downloading to {image_savedir}'):
            
            # Get the url and filename for each frame
            url, filename = frame['url'], frame['filename']
            if verbose:
                print(filename)
            
            # Determine the output filename
            outfile = os.path.join(image_savedir, filename)
            
            # Write each image to disk, checking if it exists, how big the file size is (if it exists), and whether existing files should be overwritten       
            if os.path.exists(savename):
                lc = pd.read_csv(savename, index_col = 0)
            else:
                lc = pd.DataFrame()
            if os.path.exists(outfile):
                if overwrite_images:
                    if verbose:
                        print(f'{filename} already exists. Overwriting.')
                elif os.path.getsize(outfile) < 1:
                    if verbose:
                        print('File exists, but it is 0 bytes. Overwriting.')
                elif filename in list(lc.index):
                    if verbose:
                        print('Photometry for file has already been extracted. Skipping.')
                    existing_files += 1
                    continue
                else:
                    if verbose:
                        print(f'{filename} already exists. Will not overwrite.')
                    existing_files += 1
                    continue
            with open(outfile, 'wb+') as f:
                try:
                    f.write(requests.get(url).content)
                    new_files += 1
                except:
                    failed_files += 1
                    
        print(f'\tTotal Files Found: {len(frames)}')
        print(f'\t\t{new_files} New')
        print(f'\t\t{existing_files} Already Exist')
        print(f'\t\t{failed_files} Failed to Download')
        
        
        # Extract the light curve
        lc = pd.DataFrame()
        fits_files = [os.path.join(image_savedir, file) for file in os.listdir(image_savedir) if file.endswith('.fits') or file.endswith('.fits.fz')]
        
        print('\t- - - - - - - - - -')
        for file in tqdm(fits_files, total = len(fits_files), desc = f'\tGrabbing photometry for {obj} from reduced LCOGT images.'):

            with fits.open(file) as hdus:
                if 'CAT' not in hdus:
                    # print(f'{os.path.basename(file)} has no CAT extension. Skipping.')
                    hdr = dict(hdus['SCI'].header)
                    try:
                        filt, mjd = hdr['FILTER'], float(hdr['MJD-OBS'])
                    except:
                        filt, mjd = hdr['FILTER'], np.nan
                    lc.loc[file, ['mjd', 'filter']] = [mjd, filt]
                    continue
                else:
                    hdr, cat = dict(hdus['SCI'].header), pd.DataFrame(hdus['CAT'].data)
                    try:
                        filt, mjd = hdr['FILTER'], float(hdr['MJD-OBS'])
                    except:
                        filt, mjd = hdr['FILTER'], np.nan
                    
                    if 'ra' not in cat.columns:
                        lc.loc[file, ['mjd', 'filter', 'mag', 'mag_err']] = [mjd, filt, np.nan, np.nan]
                        
                    else:
                        dras, ddecs = abs(ra - cat['ra']), abs(dec - cat['dec'])
                        dists = np.sqrt(dras**2 + ddecs**2)
                        i_close = np.argmin(dists)
                        row = cat.iloc[i_close, :]
                        
                        flag = row['flag']
                        
                        if 'mag' in cat.columns:                
                            mag, mag_err = row['mag'], row['magerr']
                        else:
                            mag, mag_err = np.nan, np.nan
                            
                        lc.loc[file, ['mjd', 'filter', 'mag', 'mag_err', 'flag']] = [mjd, filt, mag, mag_err, flag]
                        
        lc['source'] = 'lcogt'
                        
        lc.dropna(subset = ['mjd', 'mag'], inplace = True)
        lc = mag_to_flux(lc)
        lc.sort_values(by = ['mjd'], inplace = True)
        lc.reset_index(drop = True, inplace = True)
        
        # Print some stats about the queried data
        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} LCOGT photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
                        
        if save:
            lc.to_csv(savename, index = True)        
        
        return lc
                
