# Zapline Iterator

    zapline_until_gone(data, target_freq, sfreq, win_sz=10, spot_sz=2.5, viz=False, prefix="zapline_iter")

Returns: clean data, number of iterations

Function iteratively applies the Zapline algorithm.

data: assumed that the function is a part of the MNE-Python pipeline,
the input should be an output of {MNE object}.get_data() function. The shape 
should be Trials x Sensors x Time for epochs.

target_freq: frequency + harmonics that comb-like approach will be applied
with Zapline

sfreq: sampling frequency, the output of {MNE object}.info["sfreq"]

win_sz: 2x win_sz = window around the target frequency


spot_sz: 2x spot_sz = width of the frequency peak to remove
viz: produce a visual output of each iteration,

prefix: provide a path and first part of the file 
"{prefix}_{iteration number}.png"

Requires: meegkit, mne, matplotlib, numpy