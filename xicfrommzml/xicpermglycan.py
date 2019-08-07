"""Idea is to do extract XIC (MS1) AND averaged MS2 for a list of ions
over a mzml file.

MS1 XIC:
    use the skeleton from extractxiccanx.py
    recode it to open once and get many ions at the same time.
    output as a df of columns: [neutral mass, pairs of (rt,intensity) from MS1]
    NB: the RT units must be resampled to all match.

MS2 averaged:
    cycle within a RT window
    collect MS2 matching the precursor of interest
    output as a dictionary (neutral mass, rt, [pairs of (mz,intensity) from MS2])

"""

import pymzml
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.models import Label, LabelSet, TeeHead, Arrow, Range1d
import pickle
from collections import OrderedDict
import os
import numpy as np

"""
Proposed code structure:
    Create the output np.array for the MS1 -> precursor_array [neutral_mass_of_interest, pairs_of_(rt, intensity)]
    Create a dict of np.arrays for the MS2 final outputs
    
    Create a np.array to store (mz_of_interest, ms2_observed_precursor_mz, ms2 spectrum)
    Use an np.array because I wish to store the SPECTRUM, because pymzml has a addition/division function,
        see section 1.3.2 in the manual for pymzml 2.2.5.
        
    Into a target_dataframe, Get list of mz of interest and their expected RT start and end
    
    cycle through mzml:
        slice the target_dataframe by RT to get only the hits of interest.
        cycle through target_dataframe:
            If has_peak,
            collect [mz_of_interest, rt_observed, mz_observed, intensity] ms1 value into precursor_array
            
"""