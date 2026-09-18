## Imports
import os
import gc
import h5py
import numpy as np
from tqdm import tqdm
from astropy.io import fits
from multiprocess import Pool
import logging
import resource
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (4096, hard))

## Initialization ######################################################
# The only part you have to modify unless you know what you are doing ;) 

cwd = "/Volumes/dataDRS/spirou/lbl"
outdir = '/change/me' # Directory of output files

temps = [3000,3500,4000,4500,5000,5500,6000]

bands = {
    "Y": [980,1080],
    "J": [1100,1400],
    "H": [1450,1800],
    "K": [2000,2500]
}

targets = ['AUMIC', 'EV_LAC', 'V1298TAU']
templates = targets

########################################################################

## Logs 

logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

## Functions 

def odd_ratio_mean(value, err, odd_ratio=1e-4, nmax=10):
    # Vectorized implementation of odd_ratio_mean

    # Mask NaNs
    mask = np.isfinite(value) & np.isfinite(err)
    if not np.any(mask):
        return np.nan, np.nan

    # Apply mask
    value = value[mask]
    err = err[mask]

    # Initial guess
    guess = np.nanmedian(value)

    for nite in range(nmax):
        nsig = (value - guess) / err
        gg = np.exp(-0.5 * nsig**2)
        odd_bad = odd_ratio / (gg + odd_ratio)
        odd_good = 1 - odd_bad
        w = odd_good / err**2
        guess = np.nansum(value * w) / np.nansum(w)

    bulk_error = np.sqrt(1 / np.nansum(odd_good / err**2))

    return guess, bulk_error
    
def night_bin(times, rv, drv=None, binsize=0.5):
    
    """
    Bin data by night, using a moving window of size `binsize`.
    
    Parameters
    ----------
    times : numpy.ndarray
        The time array.
    rv : numpy.ndarray
        The RV array.
    drv : numpy.ndarray, optional
        The error on the RV array. If provided, the weighted mean and variance of the RV values will be
        calculated. If not provided, the unweighted mean and variance will be calculated.
    binsize : float, optional
        The size of the moving window, in days.
    
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray
       The time, RV, and error on RV arrays for the binned data.
    """
    
    n = len(rv)
    res_times = np.empty(n)
    res_rv = np.empty(n)
    res_drv = np.empty(n)

    times_temp = []
    rv_temp = []
    drv_temp = []

    time_0 = times[0]

    res_index = 0
    for index in range(n):
        if np.abs(times[index] - time_0) <= binsize:
            times_temp.append(times[index])
            rv_temp.append(rv[index])
            if drv is not None:
                drv_temp.append(drv[index])

        else:
            
            times_temp = np.array(times_temp)
            res_times[res_index] = times_temp.mean()

            rv_temp = np.array(rv_temp)
            if drv is not None:
                drv_temp = np.array(drv_temp)

            mask = ~np.isnan(rv_temp)
        
            if mask.sum() > 0:
                if drv is not None:
                    weights = 1 / (drv_temp[mask]) ** 2
                    weights /= weights.sum()
                    average = (weights * rv_temp[mask]).sum()
                    var = (weights ** 2 * (drv_temp[mask]) ** 2).sum()
                    res_rv[res_index] = average
                    res_drv[res_index] = np.sqrt(var)
                else:
                    res_rv[res_index] = rv_temp[mask].mean()
                    res_drv[res_index] = rv_temp[mask].std()

                res_index += 1
            else:
                res_rv[res_index] = np.nan
                res_drv[res_index] = np.nan
                res_index += 1
                

            time_0 = times[index]
            times_temp = [time_0]
            rv_temp = [rv[index]]
            if drv is not None:
                drv_temp = [drv[index]]
    
    if np.abs(times[index] - time_0) <= binsize:
        times_temp = np.array(times_temp)
        res_times[res_index] = times_temp.mean()

        rv_temp = np.array(rv_temp)
        if drv is not None:
            drv_temp = np.array(drv_temp)

        mask = ~np.isnan(rv_temp)

        if mask.sum() > 0:
            if drv is not None:
                weights = 1 / (drv_temp[mask]) ** 2
                weights /= weights.sum()
                average = (weights * rv_temp[mask]).sum()
                var = (weights ** 2 * (drv_temp[mask]) ** 2).sum()
                res_rv[res_index] = average
                res_drv[res_index] = np.sqrt(var)
            else:
                res_rv[res_index] = rv_temp[mask].mean()
                res_drv[res_index] = rv_temp[mask].std()

            res_index += 1
        else:
            res_rv[res_index] = np.nan
            res_drv[res_index] = np.nan
            res_index += 1


        time_0 = times[index]
        times_temp = [time_0]
        rv_temp = [rv[index]]
        if drv is not None:
            drv_temp = [drv[index]]

    if drv is not None:
        return res_times[:res_index], res_rv[:res_index], res_drv[:res_index]
    else:
        return res_times[:res_index], res_rv[:res_index]

def get_data(target, template):

    path = f"{cwd}/lblrv/{target}_{template}"

    times_lbl = []

    dtemp_all = {T: [] for T in temps}
    sdtemp_all = {T: [] for T in temps}

    wave_start = None
    wave_end = None
    order = None

    for root, dirs, files in os.walk(path):

        for f in tqdm(sorted(files), desc="Loading FITS files"):

            if "lbl.fits" not in f:
                continue

            with fits.open(os.path.join(root,f), memmap=True) as hdul:

                times_lbl.append(
                    hdul[0].header["BJD"]
                )

                for T in temps:

                    dtemp_all[T].append(
                        hdul[1].data[f"DTEMP{T}"]
                    )

                    sdtemp_all[T].append(
                        hdul[1].data[f"SDTEMP{T}"]
                    )

                if wave_start is None:

                    wave_start = hdul[1].data["WAVE_START"]
                    wave_end   = hdul[1].data["WAVE_END"]
                    order      = hdul[1].data["ORDER"]

        gc.collect()

    times_lbl = np.array(times_lbl)

    dtemp_cube = np.stack(
        [np.array(dtemp_all[T]) for T in temps],
        axis=2
    )

    sdtemp_cube = np.stack(
        [np.array(sdtemp_all[T]) for T in temps],
        axis=2
    )

    return (
        times_lbl,
        order,
        wave_start,
        wave_end,
        dtemp_cube,
        sdtemp_cube
    )



def nightly_bin_cube(
        times_lbl,
        dtemp_cube,
        sdtemp_cube,
        nproc=None,
        chunksize=1,
                    ):
    n_temp = dtemp_cube.shape[2]
    n_line = dtemp_cube.shape[1]

    def _worker(args):
        it, il, y, dy = args
        tb, yb, dyb = night_bin(times_lbl, y, dy)
        return it, il, tb, yb, dyb

    jobs = [
        (it, il, dtemp_cube[:, il, it], sdtemp_cube[:, il, it])
        for it in range(n_temp)
        for il in range(n_line)
    ]

    with Pool(nproc) as pool:
        results = list(tqdm(
            pool.imap(_worker, jobs, chunksize=chunksize),
            total=len(jobs),
            desc="Night binning"
        ))

    # grab time_binned + its length from the first result to preallocate
    time_binned = results[0][2]
    n_bin = len(time_binned)

    dtemp_binned = np.empty((n_bin, n_line, n_temp))
    sdtemp_binned = np.empty((n_bin, n_line, n_temp))

    for it, il, tb, y, dy in results:
        dtemp_binned[:, il, it] = y
        sdtemp_binned[:, il, it] = dy

    return time_binned, dtemp_binned, sdtemp_binned

def compute_band_averages(
        wave_start,
        dtemp_binned,
        sdtemp_binned
):

    n_night = dtemp_binned.shape[0]
    n_temp  = dtemp_binned.shape[2]
    n_band  = len(bands)

    dtemp_band = np.full(
        (n_night,n_band,n_temp),
        np.nan
    )

    sdtemp_band = np.full(
        (n_night,n_band,n_temp),
        np.nan
    )

    for ib,(band_name,limits) in enumerate(
            bands.items()
    ):

        logging.info(f"Computing {band_name}")

        mask = (
            (wave_start >= limits[0]) &
            (wave_start <= limits[1])
        )

        for it in range(n_temp):

            for inight in range(n_night):

                y,dy = odd_ratio_mean(
                    dtemp_binned[inight,mask,it],
                    sdtemp_binned[inight,mask,it]
                )

                dtemp_band[inight,ib,it] = y
                sdtemp_band[inight,ib,it] = dy

    return (
        dtemp_band,
        sdtemp_band
    )


def save_hdf5(
        filename,
        time_binned,
        dtemp_band,
        sdtemp_band
):

    with h5py.File(filename,"w") as h5:

        h5["time"] = time_binned

        h5["dtemp"] = dtemp_band
        h5["sdtemp"] = sdtemp_band

        h5["temperatures"] = np.array(
            temps
        )

        h5["bands"] = np.array(
            list(bands.keys()),
            dtype="S"
        )


def main_loop(target,template):

    (
        times_lbl,
        order,
        wave_start,
        wave_end,
        dtemp_cube,
        sdtemp_cube
    ) = get_data(target,template)

    ## Night binning

    (
        time_binned,
        dtemp_binned,
        sdtemp_binned
    ) = nightly_bin_cube(
        times_lbl,
        dtemp_cube,
        sdtemp_cube
    )

    ## Computing band averages

    (
        dtemp_band,
        sdtemp_band
    ) = compute_band_averages(
        wave_start,
        dtemp_binned,
        sdtemp_binned
    )

    outfile = f"{outdir}/{target}_DTEMP_bands.h5"

    logging.info(
        f"Saving {outfile}"
    )

    save_hdf5(
        outfile,
        time_binned,
        dtemp_band,
        sdtemp_band
    )

## Main 

for i in range(len(targets)):
    try :
        logging.info('## ' + targets[i])
        main_loop(targets[i], templates[i])
        gc.collect()
    except Exception as e:
        logging.info(e)


