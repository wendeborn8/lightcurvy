import io
import re
import sys
import time
import requests
import os
import pandas as pd
import numpy as np

from lightcurvy import util
from mag_to_flux import mag_to_flux

###

def atlas_query(atlas_username = None, atlas_password = None, atlas_token = None,
                obj = None, ra = None, dec = None,
                save = True, overwrite = False, obj_save = None, save_dir = None,
                use_reduced = True, **kwargs):
    
    print('Grabbing ATLAS light curve')
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'atlas', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'atlas'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'ATLAS light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        BASEURL = "https://fallingstar-data.com/forcedphot/"
        if atlas_token is None:
            if atlas_username is None or atlas_password is None:
                print('ATLAS: If no token given, must provide username and password.')
                return
            else:
                resp = requests.post(url=f"{BASEURL}/api-token-auth/", data={'username': f'{atlas_username}', 'password': f'{atlas_password}'})
    
                if resp.status_code == 200:
                    token = resp.json()['token']
                    print(f'Your token is {token}')
                else:
                    print(f'ERROR {resp.status_code}')
                    print(resp.json())
                    return
                
        headers = {'Authorization': f'Token {atlas_token}', 'Accept': 'application/json'}
        
        print(f'\tObject: {obj}\n\tRA: {ra}\n\tDec: {dec}\n\tATLAS Username / Password / Token: {atlas_username} / {atlas_password} / {atlas_token}')
        print('\t- - - - - - - - - -')
        
        task_url = None
        ncheck = 1
        while not task_url:
            with requests.Session() as s:
                resp = s.post(f"{BASEURL}/queue/", headers=headers, data={
                    'ra': ra, 'dec': dec, 'use_reduced' : use_reduced, 'mjd_min' : 55000})
    
                if resp.status_code == 201:  # successfully queued
                    task_url = resp.json()['url']
                    print(f'\tThe task URL is {task_url}')
                elif resp.status_code == 429:  # throttled
                    message = resp.json()["detail"]
                    print(f'\t{resp.status_code} {message}')
                    t_sec = re.findall(r'available in (\d+) seconds', message)
                    t_min = re.findall(r'available in (\d+) minutes', message)
                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0]) * 60
                    else:
                        waittime = 30 * ncheck
                    ncheck += 1
                    print(f'\tWaiting {waittime} seconds')
                    time.sleep(waittime)
                else:
                    print(f'\tERROR {resp.status_code}')
                    print(resp.json())
                    sys.exit()
                    
        result_url = None
        ncheck = 1
        while not result_url:
            with requests.Session() as s:
                resp = s.get(task_url, headers=headers)
    
                if resp.status_code == 200:  # HTTP OK
                    waittime = ncheck * 30
                    if resp.json()['finishtimestamp']:
                        result_url = resp.json()['result_url']
                        print(f"\tTask is complete with results available at {result_url}")
                        break
                    elif resp.json()['starttimestamp']:
                        print(f"\tTask is running (started at {resp.json()['starttimestamp']}). Checking again in {waittime} seconds...")
                    else:
                        print(f"\tWaiting for job to start. Checking again in {waittime} seconds...")
                    time.sleep(waittime)
                    ncheck += 1
                else:
                    print(f'\tERROR {resp.status_code}')
                    print(resp.json())
                    sys.exit()
    
        with requests.Session() as s:
            textdata = s.get(result_url, headers=headers).text
            
        lc = pd.read_csv(io.StringIO(textdata.replace("###", "")), delim_whitespace=True)
        lc.rename(columns = {'MJD' : 'mjd', 'm' : 'mag', 'dm' : 'mag_err', 
                             'uJy' : 'flux_uJy', 'duJy' : 'flux_err_uJy', 'RA' : 'ra', 'Dec' : 'dec',
                             'Obs' : 'image_id', 'F' : 'filter', 'err' : 'flag'}, inplace = True)

        lc['source'] = 'atlas'
        # lc.query('mag > 0 and mag < 25', inplace = True)
        lc = mag_to_flux(lc)
        lc.sort_values(by = ['mjd'], inplace = True)
        lc.reset_index(drop = True, inplace = True)

        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} ATLAS photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
                        
        if save:
            lc.to_csv(savename, index = False)
        
        return lc