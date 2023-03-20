# \\ observing script generation for semester 2022B

#import datetime
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skipper import planner, observe

#import make_pointings
import our_pointings

####
#### PARAMETERS
####
obsdates = np.genfromtxt ( '../scripts/s2023_dates.txt', comments='#', dtype=int)
obskeys = [ f'{x[0]:02d}-{x[1]:02d}-{x[2]:02d}' for x in obsdates[:,:3]]
obsfilters = np.genfromtxt ( '../scripts/s2023_dates.txt', comments='#', dtype=str)[:,4]
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
                    pointings=None, priorities=None, verbose=True, **kwargs ):
    print(f'DAY:       {day}')
    print(f'MONTH:     {month}')
    print(f'YEAR:      {year}')
    if verbose:
        print(f'TELEFILE:  {tele_fname}')
        print(f'COPILOT:   {copilot_fname}')      
        w = open('../resources/aart.txt','r').read()
        print(w)
            
    if pointings is None:
        halpha_pointings, oiii_pointings =  our_pointings.load_springfields()
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
                                       whichfield=None, priorities=priorities,
                                       verbose=verbose, **kwargs )
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
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser ( prog='sigjson_2023a.py', description='Merian obs planner for 2023A')
    parser.add_argument ( '--year', '-Y', action='store', default=2023, help='year of observations', type=int )
    parser.add_argument ( '--month', '-M', action='store', help='month of observations', type=int )
    parser.add_argument ( '--day', '-D', action='store', help='day of observations', type=int )    
    parser.add_argument ( '--make_figure', action='store_true', help='makes a figure that shows queued observations for tonight as queued.png')
    parser.add_argument ( '--telemetry_file', '-t', action='store' )
    parser.add_argument ( '--copilot_file', '-c', default=f'{os.environ["HOME"]}/Downloads/db_merian.fits', action='store' )    
    parser.add_argument ( '--dryrun', action='store_true' )
    args = parser.parse_args ()
    
    #plan_tomorrow ( args.day, args.month, args.year, args.telemetry_file, args.copilot_file)

    #fields = oiii_pointings['object'].str.extract(r'(.*?(?=_))')[0]
    
    night_index = np.where([(night[:3] == [args.year, args.month, args.day]).all() for idx,night in enumerate(obsdates)])[0][0]
    night = obsdates[night_index]
    
    #"../json/"$year$fmt_month$fmt_day"/output.log"
    jsondir=f'../json/{args.year}{args.month:02d}{args.day:02d}/'
    if not os.path.exists(jsondir) and not args.dryrun:
        os.makedirs(jsondir)
        
    plan_args = [night[2], night[1], night[0], args.telemetry_file, args.copilot_file,]
    plan_kwargs = dict( slot = night[3],
                        mfilt=obsfilters[night_index], 
                        is_queued=None, 
                        pad_last_hour=True, 
                        maxairmass=1.8, 
                        save=not args.dryrun, 
                        verbose=True
    )
    if args.dryrun:
        print(f'[sigjson] We are doing a dry run of {args.year}/{args.month}/{args.day}')
        is_queued = plan_tomorrow(  *plan_args, **plan_kwargs)
    else:        
        with open(f'{jsondir}output.log', 'w') as sys.stdout:
            is_queued = plan_tomorrow( *plan_args, **plan_kwargs )
        sys.stdout = sys.__stdout__ # \\ gotta reset stdout 
        
        print(open(f'{jsondir}output.log','r').read())
    
    if args.make_figure:
        halpha_pointings, oiii_pointings =      our_pointings.load_springfields()
        pointings_d = {'N540':oiii_pointings,'N708':halpha_pointings} 
        pointings = pointings_d[obsfilters[night_index]]
        coo = observe.CopilotOutput ( args.copilot_file, pointings )

        to_obs = is_queued.loc[~is_queued['qstamp'].isna()]
        
        fig, ax = plt.subplots(1,1, figsize=(20,3))

        plt.scatter ( pointings['RA'], pointings['dec'], 
                    facecolor='None', edgecolor='lightgrey', s=30**2, lw=1, label='planned' )
        plt.scatter ( coo.merian_sidecar['racenter'], coo.merian_sidecar['deccenter'], s=30**2, color='lightgrey', label='executed')
        pointings = pointings_d[obsfilters[night_index]]
        plt.scatter ( pointings.reindex(to_obs.index)['RA'], pointings.reindex(to_obs.index)['dec'], 
                    facecolor='None', edgecolor='r', s=30**2, lw=1, label='queued [w. padding]' )
        plt.title ( obsfilters[night_index] )
        plt.ylim(-2.2,5)
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')     
        plt.tight_layout ()  
        plt.legend () 
        plt.savefig('./queued.png')
        print('Made queued plot at queued.png')
    print('Finished planning!')