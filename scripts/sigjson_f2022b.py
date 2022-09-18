# \\ observing script generation for semester 2022B

#import datetime
#import sys
import numpy as np
import pandas as pd
from skipper import planner

#import make_pointings
import our_pointings

####
#### PARAMETERS
####
obsdates = np.genfromtxt ( '../scripts/f2022_dates.txt', comments='#', dtype=int)
obskeys = [ f'{x[0]:02d}-{x[1]:02d}-{x[2]:02d}' for x in obsdates[:,:3]]
obsfilters = np.genfromtxt ( '../scripts/f2022_dates.txt', comments='#', dtype=str)[:,4]
_field_priorities = {'VVDSearly':0, 'VVDSlate':0, 'VVDS':0, 'XMM':1, 'btwnXV':2, 'XMMhigh':4, 'newRAbtwnXV':5}
####
####

#\\ years = obsdates[:,0]
#\\ months = obsdates[:,1]
#\\ days = obsdates[:,2]
#\\ slots = obsdates[:,3]


#'VVDSearly':0, 'VVDSlate':1, 'VVDS':1, 'btwnXV':2, 'XMM':3}   

def whichfield ( year, month, day ):
    key = f'{year:02d}-{month:02d}-{day:02d}'
    index = obskeys.index(key)
        
    mfilt = obsfilters[index]
    nightslot = obsdates[index, 3]
    return mfilt, None, nightslot

def plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mfilt=None, 
                    pointings=None, priorities=None, **kwargs ):  
    if pointings is None:
        halpha_pointings, oiii_pointings =  our_pointings.load_fallfields()
    if priorities is None:
        priorities = _field_priorities
    
    if pointings is not None:
        mastercat = pointings
    elif mfilt == 'N540':
        mastercat = oiii_pointings
    elif mfilt == 'N708':
        mastercat = halpha_pointings
    else:
        raise ValueError (f"Filter {mfilt} not recognized.")
    
    is_queued = planner.plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mastercat,
                                       whichfield=whichfield, priorities=priorities, **kwargs )
    return is_queued
