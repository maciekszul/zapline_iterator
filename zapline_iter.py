import numpy as np
from scipy.signal import welch
from meegkit.dss import dss_line
import matplotlib.pylab as plt
import numpy as np
import copy


def nan_basic_interp(array):
    nans, ix = np.isnan(array), lambda x: x.nonzero()[0]
    array[nans] = np.interp(ix(nans), ix(~nans), array[~nans])
    return array



def zapline_until_gone(data, fline, sfreq, win_sz=10, spot_sz=2.5, nfft=512, show=False, prefix="dss_iter"):
    """Iterative application of the dss_line power line artifact removal

    Parameters
    ----------
    data: data, shape=(n_samples, n_chans, n_trials)
        Input data.
    fline: float
        Line frequency
    sfreq:  float
        Sampling frequency
    win_sz: float
        Half of the width of the window around the target frequency that is
        used to fit the polynomial (default=10)
    spot_sz: float
        Half of the width of the window around the target frequency that is
        used to remove the peak and interpolate (default=2.5)
    nfft: int
        FFT size for the internal PSD calculation (default=512)
    viz: bool
        Produce a visual output of each iteration. (default=False)
    prefix: str
        Path and first part of the visualisation output file 
        "{prefix}_{iteration number}.png" (default="dss_iter")

    Returns
    -------
    data: array, shape=(n_samples, n_chans, n_trials)
        Denoised data.
    iterations: number of iterations
    """

    iterations = 0
    aggr_resid = []

    freq_rn = [fline - win_sz, fline + win_sz]
    freq_sp = [fline - spot_sz, fline + spot_sz]
    freq, psd = welch(data, fs=sfreq, nfft=nfft, axis=0)

    freq_rn_ix = [
        np.where(freq >= freq_rn[0])[0][0], 
        np.where(freq <= freq_rn[1])[0][-1]
    ]
    freq_used = freq[freq_rn_ix[0]:freq_rn_ix[1]]
    freq_sp_ix = [
        np.where(freq_used >= freq_sp[0])[0][0], 
        np.where(freq_used <= freq_sp[1])[0][-1]
    ]

    if len(psd.shape) == 3:
        mean_psd = np.mean(psd, axis=(1,2))[freq_rn_ix[0]:freq_rn_ix[1]]
    if len(psd.shape) == 2:
        mean_psd = np.mean(psd, axis=(1))[freq_rn_ix[0]:freq_rn_ix[1]]

    mean_psd_wospot = copy.copy(mean_psd)
    mean_psd_wospot[freq_sp_ix[0]: freq_sp_ix[1]] = np.nan
    mean_psd_tf = nan_basic_interp(mean_psd_wospot)
    pf = np.polyfit(freq_used, mean_psd_tf, 3)
    p = np.poly1d(pf)
    clean_fit_line = p(freq_used)

    max_psd = []
    max_mean = []
    max_resid = []

    norm_vals=[]
    resid_lims=[]
    while True:
        iterations += 1
        data, art = dss_line(data, fline, sfreq, nremove=1)
        del art
        freq, psd = welch(data, fs=sfreq, nfft=nfft, axis=0)
        if len(psd.shape) == 3:
            mean_psd = np.mean(psd, axis=(1,2))[freq_rn_ix[0]:freq_rn_ix[1]]
        if len(psd.shape) == 2:
            mean_psd = np.mean(psd, axis=(1))[freq_rn_ix[0]:freq_rn_ix[1]]

        residuals = mean_psd - clean_fit_line
        mean_score = np.mean(residuals[freq_sp_ix[0]: freq_sp_ix[1]])
        aggr_resid.append(mean_score)
        
        print("Iteration {} score: {}".format(iterations, mean_score))

        if show:
            f, (ax1, ax2, ax3, ax4) = plt.subplots(
                1, 4, figsize=(12,6), 
                facecolor="white", 
                gridspec_kw={"wspace":0.2}
            )

            if len(psd.shape) == 3:
                mean_sens = np.mean(psd, axis=2)

            if len(psd.shape) == 2:
                mean_sens = psd

            for sensor in range(mean_sens.shape[1]):
                psd_s = mean_sens[:,sensor]
                y = psd_s[freq_rn_ix[0]:freq_rn_ix[1]]
                ax1.plot(freq_used, y)
            ax1.set_title("Mean PSD \nacross trials\n")
            max_psd.append(mean_sens[freq_rn_ix[0]:freq_rn_ix[1],:].max())
            ax1.set_ylim([0, max_psd[0]])

            ax2.plot(freq_used, mean_psd_tf, c="gray")
            ax2.plot(freq_used, mean_psd, c="blue")
            ax2.plot(freq_used, clean_fit_line, c="red")
            ax2.set_title("Mean PSD across \ntrials and sensors\n")
            max_mean.append(mean_psd.max())
            ax2.set_ylim([0, max_mean[0]])


            tf_ix = np.where(freq_used <= fline)[0][-1]
            ax3.plot(residuals, freq_used)
            scat_color = "green"
            if mean_score <= 0:
                scat_color = "red"
            ax3.scatter(residuals[tf_ix], freq_used[tf_ix], c=scat_color)
            ax3.set_title("Residuals\n")
            max_resid.append(residuals.max())
            ax3.set_xlim([-max_resid[0], max_resid[0]])

            ax4.scatter(np.arange(iterations), aggr_resid)
            ax4.set_title("Iterations\n")
            ax4.set_ylim([0 - (aggr_resid[0]*0.1), aggr_resid[0] + (aggr_resid[0]*0.1)])
            
            plt.savefig("{}_{}.png".format(prefix, str(iterations).zfill(3)))
            plt.close("all")
        
        if mean_score <= 0:
            break

        iterations += 1

    return [data, iterations]

