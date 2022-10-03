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

def plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mfilt=None, slot=None,
                    pointings=None, priorities=None, **kwargs ):
    print(f'DAY:       {day}')
    print(f'MONTH:     {month}')
    print(f'YEAR:      {year}')
    print(f'TELEFILE:  {tele_fname}')
    print(f'COPILOT:   {copilot_fname}')      
    w = open('../resources/aart.txt','r').read()
    print(w)
            
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
                                       current_slot=slot,
                                       current_filter = mfilt,
                                       whichfield=None, priorities=priorities, **kwargs )
    return is_queued

def duplicate_obsscript (halpha_pointings, json_input, output_filename ):
    from astropy import coordinates
    from astropy import units as u
    from skipper import observe 
    jscript = pd.read_json(json_input)
    
    planned_oiii_obscoords = coordinates.SkyCoord( jscript['RA'], 
                                                   jscript['dec'], 
                                                   unit=('deg','deg'))
    # \\ match via coordinates
    # \\ load N708 pointings
    halpha_locations = coordinates.SkyCoord( halpha_pointings['RA'], halpha_pointings['dec'], 
                                            unit=('deg','deg'))
    matchids, sep, _ = planned_oiii_obscoords.match_to_catalog_sky ( halpha_locations )
    matches = sep < 0.05*u.arcsec

    assert matches.sum() == len(planned_oiii_obscoords)

    match_indices = halpha_pointings.index[matchids[matches]]

    dummy_obscat = observe.ObsCatalog ( propid='2020B-0288', seqid='2022B')

    dummy_obscat.to_json ( halpha_pointings.reindex(match_indices), fp=output_filename, 
                          insert_checksky_exposures=False)    