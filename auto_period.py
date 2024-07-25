import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy

from scipy.ndimage import median_filter as MDF
from scipy.signal import savgol_filter as SGF
from scipy import signal
from scipy.interpolate import splev, splrep

from astropy.timeseries import LombScargle as LS
from pdmpy import pdm as PDM
from PyAstronomy.pyTiming import pyPDM

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, RationalQuadratic, ExpSineSquared

#%%

def auto_period(data = None, x_col = 'time', y_col = 'value', ye_col = 'error',
                file = None, sep = ',',
                xs = None, ys = None, yes = None, 
                methods = ['ls', 'pdm'],
                plot = True, normalize = False, subtract_polyfit = 1, n_passes = 5, 
                n_mc = 1, rand_frac = 1, combine_method = 'multiply',
                n_points = 1000, period_range = [0.1, 19.9], power_min = 0.1,
                fap = None, n_terms = 1,
                skip = 1
                ):
    
    
    ##### Read Data #####
    # If both data and file are not provded
    if data is None and file is None:
        if xs is not None and ys is not None:
            if not all((len(xs) == len(yes), len(xs) == len(ys))):
                print('Lengths of \"xs\" ({len(xs)}), \"ys\" ({len(ys)}), and \"yes" ({len(yes)}) are not all the same.')
                return
        else:
            print('If \"data\" and \"file\" are not provided, must at least provide \"xs\" and \"ys\".')
            return
        
        
    # If data is not provided but file is
    else:
        if data is None and file is not None:
            if isinstance(file, str):
                try:
                    data = pd.read_csv(file, sep = sep)
                    if len(data.columns) < 3:
                        print(f'Only {len(data.columns)} columns read from {file=} with {sep=}. Make sure file and/or sep are correct.')
                        return
                except:
                    print(f'Couldn\'t read {file=} as a Pandas DataFrame with {sep=}')
                    return
            else:
                print('If neither \"data\" nor \"xs\" and \"ys\" are provided, must provide \"file\" as a filename.')
                return
            
        # Make columns are correct
        if x_col not in data.columns:
            print(f'{x_col=} not in columns of resulting data. Assuming instead it is the 1st column.')
            x_col = list(data.columns)[0]
        if y_col not in data.columns:
            print(f'{y_col=} not in columns of file. Assuming instead it is the 2nd column.')
            y_col = list(data.columns)[1]
        if ye_col not in data.columns:
            print(f'{ye_col=} not in columns of file. Assuming instead it is the 3rd column.')
            ye_col = data.columns[2]
            
        # Remove nans from "data" DataFrame
        data.dropna(subset = [x_col, y_col, ye_col], inplace = True)
            
        # Grab xs, ys, and yes from data
        xs = data[x_col]
        ys = data[y_col]
        yes = data[ye_col]
    
    # Make sure they are numpy arrays
    # Pandas Series don't always play well with certain functions
    xs, ys, yes = np.array(xs).astype(float), np.array(ys).astype(float), np.array(yes).astype(float)

                    
    # Process the uncertainties
    
    if subtract_polyfit:
        ys_polyfit = np.poly1d(np.polyfit(xs, ys, subtract_polyfit))(xs)
        ys = ys - ys_polyfit + np.nanmedian(ys)
    
    LS_methods = ['ls', 'lombscargle', 'lomb-scargle', 'lomb']
    PDM_methods = ['pdm', 'phasedispersion', 'phasedispersionminimization']
    AC_methods = ['ac', 'auto', 'autocorrelate', 'autocorrelation', 'autocorr']
    
    results_all = {}
    i = 1
    while True:
        print(f'Pass {i}')
        results = {}
        for method in methods:
            method.lower()
            method.replace(' ', '')
            if method.lower() in LS_methods:
                results['ls'] = auto_ls(xs, ys, yes, period_range, n_points, n_mc, rand_frac, fap, n_terms)
            if method.lower() in PDM_methods:
                results['pdm'] = auto_pdm(xs, ys, yes, period_range, n_points, n_mc, rand_frac)
            if method.lower() in AC_methods:
                results['ac'] = auto_ac(xs, ys, yes, period_range, n_points, n_mc, rand_frac)
        # return results
        if len(methods) > 1:
            results['combined'] = combine(results, method = combine_method)
        else:
            results['combined'] = {}
            for key, val in results[methods[0]].items():
                results['combined'][key] = val
            
        results['peaks'] = find_peaks(results['combined'], power_min = power_min)
        
        if len(results['peaks']['period_peaks']):
            P = results['peaks']['period_peaks'][0]
            results_all[f'Pass {i}'] = results
            print(f'    - - - - - - - - - -\n    Period = {P:.4f}')
        else:
            P = None
        
        if P is not None:
            gp_fit = GP_fit(xs, ys, yes, P = P, skip = skip)
            gp, ys_gp_sub = gp_fit['gp'], gp_fit['ys_gp_sub']
        
        if plot:
            fig, (ax_lc, ax_pgrams, ax_fold) = plt.subplots(3, 1, figsize = (7, 7))
            
            # Light curve
            ax_lc.scatter(xs, ys, s = 3, alpha = 0.5, color = 'grey')
            # ax_lc.invert_yaxis()
            ax_lc.set_ylabel('Value')
            ax_lc.set_xlabel('Time')
            ax_lc.set_title(f'Light Curve - Pass {i}')
            
            # Periodogram(s)
            for method in methods:
                ax_pgrams.plot(results[method]['period'], results[method]['power'], linewidth = 1, label = method.upper())
            
            if len(methods) > 1:
                ax_pgrams.plot(results['combined']['period'], results['combined']['power'], color = 'black', linewidth = 1.5, label = 'Combined')
            ax_pgrams.set_ylabel('Power')
            ax_pgrams.set_xlabel('Period [days]')
            ax_pgrams.set_title(f'Periodograms')
            ax_pgrams.set_xlim(period_range[0], period_range[1])
            ax_pgrams.legend()
            ax_pgrams.set_ylim(-0.01, 1.01)
            ax_pgrams.axhline(power_min, color = 'grey', linestyle = ':', linewidth = 1, alpha = 0.5, zorder = 0)
            if P is not None:
                ax_pgrams.axvline(P, color = 'grey', linestyle = '--', linewidth = 1, zorder = 0, alpha = 0.5)
            
            # Phase-folded light curve
            if P is not None:
                xs_interp = np.linspace(min(xs), max(xs), 5000)
                phase = (xs % P) / P
                phase_interp = (xs_interp % P) / P
                ax_lc.plot(xs_interp, gp.predict(phase_interp.reshape(-1, 1)), color = 'black', linewidth = 1)
                ax_fold.scatter(phase, ys, c = 'grey', s = 3, alpha = 0.5)
                ax_fold.scatter(phase_interp, gp.predict(phase_interp.reshape(-1, 1)), s = 0.5, color = 'black', alpha = 0.25)
                ax_fold.set_title(f'Light Curve Folded to {P:.4f}-day Period')
                # ax_fold.invert_yaxis()
            else:
                ax_fold.set_title(f'No Suitable Period Detected - Cannot Phase-Fold')
            ax_fold.set_ylabel('Magnitude')
            ax_fold.set_xlabel('Phase')
            
            fig.tight_layout()
            plt.show()
        
        results['period'] = P
        results_all[f'Pass {i}'] = results
        i += 1
        if i > n_passes:
            break
        if P is None:
            break
        else:
            ys = copy.deepcopy(ys_gp_sub)
            
    results_all['n_mc'] = n_mc
    results_all['n_points'] = n_points
    results_all['rand_frac'] = rand_frac
    results_all['skip'] = skip
    results_all['normalize'] = normalize
    results_all['fap'] = fap
    results_all['period_range'] = period_range
    results_all['methods'] = methods
    results_all['subtract_polyfit'] = subtract_polyfit
                    
    return results_all

###

# Automatically perform the Lomb-Scargle periodogram on a set of x and y data
# Optionally provide uncertainties and use a bootstrap method
def auto_ls(xs, ys, yes, 
            period_range = [0.1, 19.9], n_points = 1000,
            n_mc = 1, rand_frac = 1, power_min = 0.1,
            fap = 0.01, n_terms = 1):
    
    results = {}
    
    xs, ys, yes = np.array(xs), np.array(ys), np.array(yes)
    periods = np.linspace(min(period_range), max(period_range), n_points)
    frequencies = 1 / periods
    results['period'] = periods
    results['frequency'] = frequencies
    
    if n_mc in [0, 1, False, None, np.nan]:
        model = LS(xs, ys, nterms = n_terms)
        powers = model.power(frequencies)
        results['power'] = powers
        results['power_high'] = None
        results['power_low'] = None
        
        if fap in [0, 100, None, np.nan, False]:
            results['power_fap'] = None
            results['power_fap_high'] = None
            results['power_fap)low'] = None
        else:
            power_fap = model.false_alarm_level(fap)
            results['power_fap'] = power_fap
            results['power_fap_high'] = None
            results['power_fap)low'] = None

    else:
        ys_rand_all = np.random.normal(loc = ys, scale = yes, size = (n_mc, len(ys)))
        inds_rand = [np.random.choice(range(len(xs)), size = int(rand_frac*len(xs)), replace = False) for _ in range(n_mc)]
        
        powers_all, power_fap_all = [], []
        for i, inds in tqdm(enumerate(inds_rand), total = n_mc, desc = '    Lomb-Scargle'):
            xs_rand = xs[inds]
            ys_rand = ys_rand_all[i][inds]
            model = LS(xs_rand, ys_rand, nterms = n_terms)
            powers = model.power(frequencies)
            powers_all.append(powers)
            
            if fap in [0, 100, None, np.nan, False]:
                power_fap_all.append(None)
            else:
                power_fap = model.false_alarm_level(fap)
                power_fap_all.append(power_fap)
                
        powers, power_high, power_low = np.nanpercentile(powers_all, (50, 84.1, 15.9), axis = 0)
        power_fap, power_fap_high, power_fap_low = np.nanpercentile(power_fap_all, (50, 84.1, 15.9))
        
        results['power'] = powers
        results['power_high'] = power_high
        results['power_low'] = power_low
        results['power_fap'] = power_fap
        results['power_fap_high'] = power_fap_high
        results['power_fap_low'] = power_fap_low
        
    
    # Window Power
    ones = np.ones(len(ys))
    powers_window = LS(xs, ones, fit_mean = False, center_data = False).power(frequencies)
    results['power_window'] = powers_window
    
    # Peaks
    results = find_peaks(results, power_min = power_min)
    
    results['n_mc'] = n_mc
    results['rand_frac'] = rand_frac
    results['fap'] = fap
    results['n_terms'] = n_terms
    
    return results

###

#

def auto_pdm(xs, ys, yes,
             period_range = [0.1, 19.9], n_points = 1000,
             n_mc = 1, rand_frac = 1, power_min = 0.1):
    
    results = {}
    
    xs, ys, yes = np.array(xs), np.array(ys), np.array(yes)
    # periods = np.linspace(min(period_range), max(period_range), n_p)
    # frequencies = 1 / periods
    fmin, fmax = 1/max(period_range), 1/min(period_range)
    dfreq = (fmax - fmin) / n_points
    
    dval = (max(period_range) - min(period_range)) / n_points
    scanner = pyPDM.Scanner(minVal=min(period_range), maxVal=max(period_range), dVal=dval, mode="period")
        
    if n_mc in [0, 1, False, None, np.nan]:
        # frequencies, thetas = PDM(xs, ys, f_min = fmin, f_max = fmax, delf = dfreq, nbin = 25)
        # periods = 1 / frequencies
        # periods, frequencies, thetas = periods[::-1], frequencies[::-1], thetas[::-1]
        periods, thetas = pyPDM.PyPDM(xs, ys).pdmEquiBinCover(10, 3, scanner)
        frequencies = 1 / periods
        
        powers = abs(1 - thetas)
        powers = MDF(powers, size = 51)
        powers = SGF(powers, 21, 3)
        
        results['period'] = periods
        results['frequency'] = frequencies
        results['power'] = powers
        results['power_high'] = None
        results['power_low'] = None
        results['power_fap'] = None
        results['power_fap_high'] = None
        results['power_fap_low'] = None
        
    else:   
        ys_rand_all = np.random.normal(loc = ys, scale = yes, size = (n_mc, len(ys)))
        inds_rand = [np.random.choice(range(len(xs)), size = int(rand_frac*len(xs)), replace = False) for _ in range(n_mc)]
        
        powers_all, power_fap_all = [], []
        for i, inds in tqdm(enumerate(inds_rand), total = n_mc, desc = '    Phase Dispersion Minimization'):
            xs_rand = xs[inds]
            ys_rand = ys_rand_all[i][inds]
            # frequencies, thetas = PDM(xs_rand, ys_rand, f_min = fmin, f_max = fmax, delf = dfreq, nbin = 25)
            # periods = 1 / frequencies
            # periods, frequencies, thetas = periods[::-1], frequencies[::-1], thetas[::-1]
            periods, thetas = pyPDM.PyPDM(xs_rand, ys_rand).pdmEquiBinCover(10, 3, scanner)
            frequencies = 1 / periods
            
            powers = abs(1 - thetas)
            powers_all.append(powers)
            
            if i == 0:
                results['period'] = periods
                results['frequency'] = frequencies
                
        powers, power_high, power_low = np.nanpercentile(powers_all, (50, 84.1, 15.9), axis = 0)
        powers = SGF(powers, 21, 3)
        power_high = SGF(power_high, 21, 3)
        power_low = SGF(power_low, 21, 3)
        
        results['power'] = powers
        results['power_high'] = power_high
        results['power_low'] = power_low
        results['power_fap'] = None
        results['power_fap_high'] = None
        results['power_fap_low'] = None
        
    # Peaks
    results = find_peaks(results, power_min = power_min)
    
    results['n_mc'] = n_mc
    results['rand_frac'] = rand_frac
    results['fap'] = None
    results['n_terms'] = None
    
    return results

###

def auto_ac(xs, ys, yes,
            period_range = [0.1, 19.9], n_points = 1000,
            n_mc = 1, rand_frac = 1):
    
    results = {}
    xs, ys, yes = np.array(xs), np.array(ys), np.array(yes)
    period_min, period_max = min(period_range), max(period_range)
    
    # Resample to an even grid
    diff = np.nanmedian(abs(np.diff(xs))) / 4
    diff = (max(xs) - min(xs)) / n_points

    x_interp = np.arange(min(xs), max(xs), diff)
    # x_interp = np.linspace(min(xs), max(xs), n_points)
    y_interp = np.interp(x_interp, xs, ys)
    # plt.scatter(x_interp, y_interp, s = 1)
    # y_interp = y_interp / np.linalg.norm(y_interp) # normalize y_interp
    
    if n_mc in [0, 1, False, None, np.nan]:
        corr = signal.correlate(y_interp, y_interp)
        lags = signal.correlation_lags(len(x_interp), len(x_interp))
        # corr /= corr[list(lags).index(0)]
        
        corr = pd.DataFrame({'lag' : lags*diff, 'power' : corr})
        corr.query('lag > @period_min and lag < @period_max', inplace = True)
        
        yfit = np.poly1d(np.polyfit(corr['lag'], corr['power'], deg = 1))(corr['lag'])
        corr['power'] -= yfit
        
        corr.loc[:, 'power'][corr['power'] < 0] = 0
        corr.reset_index(inplace = True, drop = True)

        corr['power'] /= max(corr['power'])
        i_first0 = np.argmin(np.array(corr['power']))
        lag_i_first0 = np.array(corr['lag'])[i_first0]
        corr.loc[:, 'power'][corr['lag'] < lag_i_first0] = 0
        corr.sort_values('lag', inplace = True)
                
        results['period'] = np.array(corr['lag'])
        results['frequency'] = 1 / np.array(corr['lag'])
        results['power'] = np.array(corr['power'])
        results['power_high'] = None
        results['power_low'] = None
        results['power_fap'] = None
        results['power_fap_high'] = None
        results['power_fap_low'] = None
    
    return results

###

def auto_dcf(xs, ys, yes,
            period_range = [0.1, 19.9], n_points = 1000,
            n_mc = 1, rand_frac = 1):
    
    results = {}
    
    xs, ys, yes = np.array(xs), np.array(ys), np.array(yes)
    # periods = np.linspace(min(period_range), max(period_range), n_p)
    # frequencies = 1 / periods
    # fmin, fmax = 1/max(period_range), 1/min(period_range)
    # dfreq = (fmax - fmin) / n_points
    
    lags = np.linspace(min(period_range), max(period_range), n_points)
    dt = (max(period_range) - min(period_range)) / n_points
        
    if n_mc in [0, 1, False, None, np.nan]:
        lc = np.array([xs, ys, yes]).T
        dcf, dcf_err = GDCF(lc, lc, lags, dt)
        
        width = int(0.05 * n_points)
        # dcf = SGF(dcf, window_length = 31, polyorder = 3)
        dcf = MDF(dcf, size = width)
        
        results['period'] = lags
        results['frequency'] = 1 / lags
        results['power'] = dcf
        results['power_high'] = None
        results['power_low'] = None
        
        # t1, t2 = lc1.iloc[:, 0], lc2.iloc[:, 0]
        # y1, y2 = lc1.iloc[:, 1], lc2.iloc[:, 1]
        # e1, e2 = lc1.iloc[:, 2], lc2.iloc[:, 2]
        
        # if dt is None:
        #     dt1, dt2 = np.nanmean(np.diff(t1)), np.nanmean(np.diff(t2))
        #     dt = np.mean([dt1, dt2]) / 2
            
        # if lag_min is None and lag_max is None:
        #     lag_max = min([max(t1) - min(t1), max(t2) - min(t2)])
        #     lag_min= -lag_max
        
        
        # N = np.around((lag_max - lag_min) / float(dt))
        # lags = np.linspace(lag_min + (dt/2.0), lag_max - (dt/2.0), int(N))
        
        # y1s = np.random.normal(loc = y1, scale = 0.1*y1, size = (n_mc, len(y1)))
        # y2s = np.random.normal(loc = y2, scale = 0.1*y2, size = (n_mc, len(y2)))
        
        # dcfs, dcf_errs = [], []
        # for i in range(n_mc):
            
        #     _lc1 = np.array([t1, y1s[i], e1]).T
        #     _lc2 = np.array([t2, y2s[i], e2]).T
        # dcf, dcf_err = gdcf(_lc1, _lc2, lags, dt)
        
    else:
        pass
    
    return results

###

def combine(results, normalize = True, method = 'x', power_min = 0.1):
    
    periods, powers = [], []
    powers_high, powers_low = [], []
    for key, value in results.items():
        periods.append(value['period'])
        powers.append(value['power'])
        powers_high.append(value['power_high'])
        powers_low.append(value['power_low'])
        
    results_combined = {}
        
    # Get the period array
    # 1. Use all points
    periods_all = []
    for period in periods:
        periods_all += list(period)
    periods_all = np.array(sorted(periods_all))
    
    if normalize:
        powers_combined = np.ones(len(periods_all))
        powers_high_combined = np.ones(len(periods_all))
        powers_low_combined = np.ones(len(periods_all))
    
    powers_all, powers_high_all, powers_low_all = [], [], []
    norm = 1
    for period, power, power_high, power_low in zip(periods, powers, powers_high, powers_low):
        power_interp = np.interp(periods_all, period, power)
        powers_all.append(power_interp)
            
        if power_high is not None:
            power_high_interp = np.interp(periods_all, period, power_high)
            powers_high_all.append(power_high_interp)
                
        if power_low is not None:
            power_low_interp = np.interp(periods_all, period, power_low)
            powers_low_all.append(power_low_interp)
    
    if method == 'mean':
        powers_combined = np.nanmean(powers_all, axis = 0)
        if power_high is not None:
            powers_high_combined = np.nanmean(powers_high_all, axis = 0) 
        if power_low is not None:
            powers_low_combined = np.nanmean(powers_low_all, axis = 0)
            
    elif method in ['x', 'X', 'mult', 'multiply', 'times']:
        powers_combined = np.prod(powers_all, axis = 0) / np.nanmedian(powers_all, axis = 0)
        if power_high is not None:
            powers_high_combined = np.prod(powers_high_all, axis = 0) / np.nanmedian(powers_high_all, axis = 0)
        if power_low is not None:
            powers_low_combined = np.prod(powers_low_all, axis = 0) / np.nanmedian(powers_low_all, axis = 0)
            
    results_combined['period'] = periods_all
    results_combined['power'] = powers_combined
    if power_high is not None:
        results_combined['power_high'] = powers_high_combined
        results_combined['power_low'] = powers_low_combined
        
    results = find_peaks(results_combined, power_min = power_min)
    
    return results_combined

###

def GDCF(ts1, ts2, t, dt):

    '''
        Subroutine - gdcf
          DCF algorithm with gaussian weighting
    '''

    h = dt/4.0
    gkrn = lambda x: np.exp(-1.0 * np.abs(x)**2 / (2.0 * h**2)) \
           / np.sqrt(2.0 * np.pi * h)
    cntrbt = gkrn(3.290527*h)

    dcf = np.zeros(t.shape[0])
    dcferr = np.zeros(t.shape[0])
    n = np.zeros(t.shape[0])

    dst = np.empty((ts1.shape[0], ts2.shape[0]))
    for i in range(ts1.shape[0]):
        for j in range(ts2.shape[0]):
            dst[i,j] = ts2[j,0] - ts1[i,0]

    for k in range(t.shape[0]):
        gdst = gkrn(dst - t[k])
        ts1idx, ts2idx = np.where(gdst >= cntrbt)

        mts2 = np.mean(ts2[ts2idx,1])
        mts1 = np.mean(ts1[ts1idx,1])
        n[k] = ts1idx.shape[0]

        dcfdnm = np.sqrt((np.var(ts1[ts1idx,1]) - np.mean(ts1[ts1idx,2])**2) \
                         * (np.var(ts2[ts2idx,1]) - np.mean(ts2[ts2idx,2])**2))

        dcfs = (ts2[ts2idx,1] - mts2) * (ts1[ts1idx,1] - mts1) / dcfdnm
        dcf[k] = np.sum(dcfs) / float(n[k])
        dcferr[k] = np.sqrt(np.sum((dcfs - dcf[k])**2)) / float(n[k] - 1)

    return dcf, dcferr

###

def GP_fit(xs, ys, yes, P = 1,
           skip = 1):
    
    phase = np.array((xs % P) / P)
    xs_skip, ys_skip, yes_skip, phase_skip = np.array(xs[::skip]), np.array(ys[::skip]), np.array(yes[::skip]), np.array(phase[::skip])
    
    y_mean, y_std = np.nanmedian(ys), np.nanmedian(abs(np.diff(ys)))
    
    kernel = y_mean**2 * ExpSineSquared(length_scale = 100, periodicity = 1, 
                                        length_scale_bounds = 'fixed', periodicity_bounds = 'fixed')
    gp = GaussianProcessRegressor(kernel = kernel, alpha = y_std, n_restarts_optimizer = 9)
    gp.fit(phase_skip.reshape(-1, 1), ys_skip.reshape(-1, 1))
    
    ys_gp_sub = ys - gp.predict(phase.reshape(-1, 1)) + np.nanmedian(ys)
    
    results = {}
    results['gp'] = gp
    results['ys_gp_sub'] = ys_gp_sub
    
    return results

###

def find_peaks(results, n_peaks = 5, 
               power_min = 0.05, distance_frac = 0.05):
    
    period, power = results['period'], results['power']
    
    peaks, props = signal.find_peaks(x = power, 
                                     height = power_min,
                                     distance = distance_frac*len(power), 
                                     )
    peak_inds_sorted = np.argsort([power[ind] / np.sqrt(1 + 2*period[ind]) for ind in peaks])[::-1]
    peak_inds_sorted = np.argsort([power[ind] for ind in peaks])[::-1]
    peak_inds = peaks[peak_inds_sorted][0:n_peaks]
    
    results['peaks_inds'] = peak_inds
    results['power_peaks'] = power[peak_inds]
    results['period_peaks'] = period[peak_inds]
    
    return results


