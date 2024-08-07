A small set of scripts to fetch astronomical light curve data from a variety of sources, including:
- [ASAS-SN](http://asas-sn.ifa.hawaii.edu/skypatrol/), [ZTF](https://www.ztf.caltech.edu/), [AAVSO](https://www.aavso.org/LCGv2/), [Atlas](https://fallingstar-data.com/forcedphot/), [WISE](https://irsa.ipac.caltech.edu/Missions/wise.html), [Gaia](https://gaia.aip.de/cms/data/gdr3/), [Pan-STARRS](https://catalogs.mast.stsci.edu/panstarrs/), [LCOGT](https://archive.lco.global/), and TESS.

with plans to include SkyMapper light curves, as well.

---
Using lightcurvy is pretty straightforward. After defining some basics (such as your target's name), you create a lightcurvy object:
```
obj = 'GM Aur'                   # This is the SIMBAD name of your target
obj_save = 'gmaur'               # This is how any resulting files will be saved. For example, I like to avoid spaces and any special characters in my filenames
ra, dec = None, None             # If your target does not have a SIMBAD name, you can provide RA/Dec instead, though some function may not work as intended
save = True                      # Whether to save the resulting light curves
overwrite = False                # Whether to overwrite existing light curve files
save_dir = r'D:\My Drive\Data'   # Where all of the resulting data will be saved
lcy = lightcurvy.lightcurvy(obj = obj, obj_save = obj_save, ra = ra, dec = dec,
                            save = save, overwrite = overwrite, save_dir = save_dir)
```

On its own, this won't do much. The real meat of `lightcurvy` comes from `lightcurvy.query_all`, which actually fetches the light curves. Let's again first declare some relevant parameters, then use `query_all`.

```
sources = ['asassn', 'ztf', 'atlas', 'aavso', 'wise', 'panstarrs', 'gaia', 'lcogt']    # Which light curves to search for data. This is currently all of them
search_radius = 0.0015                                                                 # Several sources require you to include some radius to search within, in degrees (about 5 arcseconds)
lcogt_start, lcogt_end = 2014-01-01                                                    # If fetching LCOGT data, the start/end dates to search for
overwrite_images = False                                                               # If fetching LCOGT images, whether to overwrite (and therefor re-download) any existing images
lc = lcy.query_all(sources = sources, search_radius = search_radius,
                   lcogt_start = lcogt_start, lcogt_end = lcogt_end, overwrite_images = overwrite_images)
```

lightcurvy will now query each source in `sources`, fetching light curves and saving them to their respective sub-directory in `save_dir`, updating you on progress as it goes. The source-specific light curves include *all* columns that are included from the original database, along with some extras that are created. These extra columns are primarily related to the source's flux, converting magnitude to flux (Jy, erg/s/cm2/A, and erg/s/cm2) using effective wavelengths and zero-points provided by [SVO](http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse).

Data from each source is fetched in different ways, so the progress/updates will look different. That also means each can take wildly different amounts of time.

AAVSO data doesn't appear to be located in a central, query-able database, but is instead spread across many pages that each needs to be scraped. If anyone knows about a better way to fetch AAVSO data, please let me know!
![image](https://github.com/user-attachments/assets/032bc7e6-0086-404c-ab8e-e31e6786e59c)

ATLAS photometry, for example, is fetched via forcing photometry on ATLAS' image database on their servers. for this reason it can take a while to initiate then run the job. If you request many in a short period opf time, ATLAS may throttle your connection/account, but I've never had this happen.
![2024-07-25 08_39_15-JWST UX Tau A - Online LaTeX Editor Overleaf](https://github.com/user-attachments/assets/bdfdf9fd-e537-42ec-9227-9796b61b2ef6)

ASAS-SN, ZTF, Gaia, Pan-STARRS, and WISE data are pretty quick to fetch, taking no more than a minute or so each.
![image](https://github.com/user-attachments/assets/44575f79-7106-43d9-81a1-d7abf6942880) ![image](https://github.com/user-attachments/assets/0663d23e-99c6-44e3-9d0d-1ec66da2c0cb)

LCOGT data is probably the most time-consuming of all. LCOGT's archive is an image archive, not a photometry/light curve archive, so no centralized cache of light curves exists. Their image reduction/calibration pipeline BANZAI does, however perform source extraction and magnitude calibration for images taken using a griz filter. So here `lightcurvy` must download the respective image and grab the photometry data from the header. If, like for GM Aur, there are thousands of images, this can take hours and consume many gigabytes of storage space. This is why I include extra parameters for LCOGT, to hopefully minimize the time as much as possible.
![image](https://github.com/user-attachments/assets/42cf7bd1-2ba6-418c-91db-1443c63ae6bd)

**It's important to note that not all targets have photometry in all sources. For example, AAVSO doesn't accept/use coordinates in their searches, so a target without a common name will probably not have any AAVSO data. Others are restricted to certain parts of the sky: ZTF is northern hemisphere only, Pan-STARRS only covers ~75% of the sky, etc, so if you don't get data from a specific source, check whether you target is even observable.**

---

Now that we've fetched our light curves, the heavy lifting of `lightcurvy` has been completed. Because the use-cases for each user will differ, deeper analyses of the resulting data is primarily left up to the user. That said, some basic tools are included for convenience. 

We can plot the light curves, though with up to seven sources of photometry in many different bands, it can get quite chaotic.
```
lcy.plot_lc()
lcy.plot_sed()
```
![image](https://github.com/user-attachments/assets/d7d4fb8e-340d-4591-8157-52fd7194ba60)
![image](https://github.com/user-attachments/assets/6ee6f537-c49f-4cba-9254-6a40ba32f258)


The resulting plots (particularly the light curves) aren't great. There are some really spurious points in the ATLAS light curves which wreck the flux scaling in the light curve plot and introduce large errorbars in the SED plot. Let's re-plot this, but remove flagged data and restrict the light curve to a more reasonable date range where the data is fairly dense.
```
lcy.plot_lc(mjd_min = 59000, mjd_max = 60000, flag = True)
lcy.plot_sed(mjd_min = None, mjd_max = None, flag = True)
```
![image](https://github.com/user-attachments/assets/49865185-3532-4385-aebc-9dbbb7fd7a0a)
![image](https://github.com/user-attachments/assets/3fd16b7a-aa04-4d78-93ad-094e329fdf17)


Much better! The SED looks (though ATLAS and Pan-STARRS are still far from perfect) and the light curve shows some structure. Manual tuning/pruning of each respective source is still probably necessary; the calibration between sources with like bands (like ASAS-SN *g* and LCOGT *gp*) isn't perfect.

---

Another useful part of `lightcurvy` is the `auto_period` module that allows for identification and removal of periodic signals from light curves.

Let's take the multi-source, multi-band light curve from above and restrict it to LCOGT *gp* data. These are data I obtained at high SNR and that are fairly trustworthy and delivers reliable periodic signatures. You can see clear structure in this light curve that might appear periodic, but it's hard to tell. 
```
lc = lcy.lc
lc.query('source == \'lcogt\'', inplace = True)
lc.query('mjd > 59400 and mjd < 59600')
g = lc.query('filter == \'gp\'')
g.query('flux_Flam > 2', inplace = True) # Remove some spruious data
plt.scatter(g['mjd'], g['flux_Flam'], s = 15, alpha = 0.75)
plt.xlabel('Date (MJD)')
plt.ylabel('Flux')
plt.show()
```

Let's do a quick run through `auto_period` with basic parameters.
```
from lightcurvy.auto_period import auto_period
x_col = 'mjd'                 # Column name for x values
y_col = 'flux_Flam'           # Column name for y values
ye_col =  'err_Flam'          # Column name for uncertainty
subtract_polyfit = 0          # What degree polynomial to fit and subtract from the light curve
n_passes = 1                  # How many times to run through the data (see below)
n_mc = 1                      # Whether to use a Monte Carlo/bootstrapping approach in period-finding. If n_mc is in [0, 1, False, None], will not use Monte Carlo approach
rand_frac = 1                 # What fraction of the data to use (only used if n_mc > 1, see below)
period_range = [0.1, 15.1]    # The range of periods to search through. Large (>~1/2 the time range) and small (very close to 0) periods can cause issues
results = auto_period.auto_period(data = g, x_col = x_col, y_col = y_col, ye_col = ye_col,
                                  subtract_polyfit = subtract_polyfit, n_passes = n_passes, 
                                  n_mc = n_mc, rand_frac = rand_frac, period_range = period_range)
```
There is a clear ~6-day period in this light curve, identified by both the Lomb-Scargle and Phase Dispersion Minimization algorithms. The top panel is the light curve with fitted Gaussian Processes, the middle is the periodograms with the peak noted, and the bottom panel is the phase-folded light curve with the Gaussian Process fit. It's rather sinusoidal, pointing to this periodic signature originating from rotation.
![image](https://github.com/user-attachments/assets/bd6547f4-30ad-4cc9-9f77-2aa279a2bbf6)

That will suffice for most, but maybe you're looking for secondary, or even tertiary periods. To do so, we can increase `n_passes` to something like 3. This will run the periodogram, subtract off that fitted Gaussian Process and re-run the periodogram 2 more times, allowing you to quickly investigate secondary periods. It also allows you to remove spurious periodic signatures that might be aliases of the primary period.
```
x_col = 'mjd'                 # 
y_col = 'flux_Flam'           # 
ye_col =  'err_Flam'          # 
subtract_polyfit = 1          # Subtract a 1st degree polynomial (linear fit) from the light curve. Long term trends can impact the periodogram
n_passes = 3                  # How many times to compute the periodogram and fit/subtract a Gaussian Process
n_mc = 100                    # Sample 100 light curves. This helps smooth over small bumps in the periodogram
rand_frac = 0.75              # Use 75% of the data for each Monte Carlo iteration. This helps reduce the the impact of outlying points that might be having an outsize impact on the periodogram
period_range = [0.1, 15.1]    # 
results = auto_period.auto_period(data = g, x_col = x_col, y_col = y_col, ye_col = ye_col,
                                  subtract_polyfit = subtract_polyfit, n_passes = n_passes, 
                                  n_mc = n_mc, rand_frac = rand_frac, period_range = period_range)
```
The first pass is much the same as before, nothing to see there, really. But now on the second pass, it subtracts that fitted GP (the solid black line in the top panel) from the light curve and recalculates the period. The 6-day period is completely gone, leaving behind a broad signal near 7.9 days. In this particular case, it doesn't complete a 3rd pass because it wasn't able to find a significant period. This behavior can be controlled by `power_min`. Perhaps with a cleaner light curve a more obvious secondary period could be found.
![image](https://github.com/user-attachments/assets/eaffe0ac-b842-4f23-91ed-3ddb648a54bc)

The orignal `lightcurvy` object also has the ability to perform `auto_period` on each filter in the light curve. Here, I've created a `lightcurvy` object from *just* LCOGT data (again, it tends to be cleaner) and run `auto_period`.
```
lcy = lightcurvy.Lightcurvy(obj = obj, obj_save = obj_save, save_dir = save_dir)
lcy.query_all(sources = ['lcogt'], overwrite = False, save = False)
lcy.plot_lc(mjd_min = 59400, mjd_max = 59600, norm = True)
lcy.multi_period(mjd_min = 59400, mjd_max = 59600)
```
Here the light curves have been normalized by the median flux and you can see the similar behaviors of each filter.
![image](https://github.com/user-attachments/assets/1064f568-880e-44de-ac5a-4178de64b6c2)

After running `multi_period`, we're left with periodograms for each filter. *gp*, *rp*, and *ip* show that clear 6-day period. *zs* is much noisier and does show signal near 6 days, but the peak happens to be near 3 days.
![image](https://github.com/user-attachments/assets/3f0de574-b6e4-4538-a9f3-0a0181057a98)

---
