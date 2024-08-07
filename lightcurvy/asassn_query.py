import os
import pandas as pd
import numpy as np
import os
import requests
import io
from base64 import encodebytes

from pyasassn.client import SkyPatrolClient
from pyasassn.utils import Vcams, gcams
client = SkyPatrolClient()

from . import util
from .mag_to_flux import mag_to_flux

#%%

def _deserialize(buffer):
    # Deserialize from bytes
    return pd.read_parquet(io.BytesIO(buffer))

def asassn_query(obj = None, ra = None, dec = None,
                 save = True, overwrite = False, obj_save = None, save_dir = None,
                 search_radius = 0.0015,
                 **kwargs):
    
    print('Grabbing ASAS-SN light curve')
        
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Make the paths and filenames
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'asassn', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'asassn'), exist_ok = True)
    
    # Check if this light curve has already been generated
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'ASAS-SN light curve already exists for {obj} ({savename}), but will not overwrite. Grabbing light curve from {savename}.')
        
        return lc
    
    else:
        print(f'\tObject: {obj}\n\tRA: {ra}\n\tDec: {dec}\n\tSearch Radius: {search_radius}')

        
        # The following ~15 lines are taken more-or-less directly from SkyPatrol. By default it uses multiprocessing to fetch data, but the overheads here can be *enormous*, so I needed to create a synchronous process instead. For large datasets this may end up being slower, but for every target I've tried it's effectively infinitely faster.
        catalog = 'master_list'
        cols = ["asas_sn_id", "ra_deg", "dec_deg", 'catalog_sources']
        download = True
        file_format = 'parquet'
        
        url = f"http://asassn-lb01.ifa.hawaii.edu:9006/lookup_cone/radius{search_radius}_ra{ra}_dec{dec}"
        response = requests.post(url, json={"catalog": catalog, "cols": cols, "format": "arrow", "download": download})
        
        tar_df = _deserialize(response.content)
        
        query_id = (f"conectr-{search_radius}_conera-{ra}_conedec-{dec}|catalog-{catalog}|cols-" + "/".join(cols))
        query_hash = encodebytes(bytes(query_id, encoding="utf-8")).decode()

        # Get lightcurve ids to pull
        tar_ids = list(tar_df["asas_sn_id"])

        n_chunks = int(np.ceil(len(tar_ids) / 1000))
        results = [client._get_lightcurve_chunk(query_hash, idx, catalog, None, file_format) for idx in range(n_chunks)]
        
        chunks = []
        for data, n in results:
            chunks.append(data)
            
        lc = pd.concat(chunks)
        
        # Add/adjust some columns
        lc['filter'] = ['V' if camera in Vcams else 'g' if camera in gcams else None for camera in lc['camera']]
        lc['flag'] = [0 if quality == 'G' else 1 if quality == 'B' else None for quality in lc['quality']]
        lc.rename(columns = {'asas_sn_id' : 'object_id', 'limit' : 'mag_limit'}, inplace = True)
        # lc.drop(columns = ['object_id', 'mag_limit', 'flux', 'flux_err', 'fwhm'], inplace = True)
        lc['source'] = 'asassn'
        lc['mjd'] = lc['jd'] - 2400000.5
        
        #
        lc = mag_to_flux(lc)
        
        # Print some stats about the queried data
        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} ASAS-SN photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
        
        # Save, if desired
        if save:
            lc.to_csv(savename, index = False)
        
        return lc