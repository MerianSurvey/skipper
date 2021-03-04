#!/usr/bin/env python
"""
IO related functions.
"""
__author__ = "Shany Danieli"

import numpy as np
import pandas as pd
import astropy.units as u

# __all__ = []


def read_night_summary(date, dirc):
    '''Read the nightly summary .dat file generated at the end
    of each night by the nightly observer 
    (inq on the observer 2 machine)
    Parameters
    ----------
    date: str
         Observation date
         format is: YYYYMMDD_DD (e.g. 20210305_06)
    dirc: str
         The path to the directory where the nightly summary file is saved.
    Returns
    -------
    observed_dates: list
        List of dates that were observed during the input date.
    Notes
    -----
    '''
    file = dirc+date+'.dat'
    datfile = open(file,'r')
    count = 0
    for line in open(file): count += 1
    index = np.arange(count)

    df = pd.DataFrame(index=index, columns = ['expnum', 'ra', 'dec', 'ut', 'filter', 'exptime','secz', 'type', 'object'])

    i = 0
    for line in datfile:
        line_lst = line.split()
        df.loc[df.index[i], 'expnum'] = line_lst[0]
        df.loc[df.index[i], 'ra'] = line_lst[1]
        df.loc[df.index[i], 'dec'] = line_lst[2]
        df.loc[df.index[i], 'ut'] = line_lst[3]
        df.loc[df.index[i], 'filter'] = line_lst[4]
        df.loc[df.index[i], 'exptime'] = line_lst[5]
        df.loc[df.index[i], 'secz'] = line_lst[6]
        df.loc[df.index[i], 'type'] = line_lst[7]
        df.loc[df.index[i], 'object'] = " ".join(line_lst[8:])
        i += 1

    observed_dates = df['object'].tolist()
    return observed_dates

