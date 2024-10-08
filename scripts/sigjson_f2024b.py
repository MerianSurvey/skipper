# \\ observing script generation for semester 2023B

#import datetime
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates, table
from astropy.io import fits
import pandas as pd
from skipper import planner, observe

#import make_pointings
import our_pointings

####
#### PARAMETERS
####
date_file = '../scripts/f2024_dates.txt'
obsdates = np.genfromtxt ( date_file, comments='#', dtype=int)
obskeys = [ f'{x[0]:02d}-{x[1]:02d}-{x[2]:02d}' for x in obsdates[:,:3]]
obsfilters = np.genfromtxt ( date_file, comments='#', dtype=str)[:,4]
_field_priorities = {'VVDSearly':100, 'VVDSlate':100, 'VVDS':100, 'XMM':100, 'btwnXV':2, 'XMMhigh':2, 'newRAbtwnXV':2, 'coobserved':1, 'newEarlyRA':3, 'mustget':0}
####
####

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
    
    # -----------------------------------------------------------------
    # -- Prioritize pointings that already have single band coverage --
    # and super-prioritize ones from Yifei's lists
    # -----------------------------------------------------------------
    # -- Set up co filter
    if mfilt == 'N540':
        cofilt = 'N708'
        co_skySB = 21.
        co_teffmin = 300.
        co_pointings = halpha_pointings
    else:
        cofilt = 'N540'
        co_skySB = 22.1
        co_teffmin = 300. # XXX
        co_pointings = oiii_pointings
    
    co_coo = observe.CopilotOutput ( copilot_fname, pointings=co_pointings,  skySB_0 = co_skySB )   
    co_finished = co_coo.identify_completed_pointings ( co_teffmin )
    observed_cofilter = coordinates.SkyCoord ( co_finished['rabore'], co_finished['decbore'], unit=('deg','deg') )
    pointings_mfilter = coordinates.SkyCoord ( mastercat['RA'], mastercat['dec'], unit=('deg','deg') )
    idx, d2d, _ = pointings_mfilter.match_to_catalog_sky ( observed_cofilter )
    has_observed_match = d2d.to('arcsec').value < 2.
    
    # \\ all observed pointings should have a match in the second filter pointing    
    assert has_observed_match.sum() == observed_cofilter.shape[0]
    mastercat['priority_name'] = mastercat['object'].str.extract(r'(.*?(?=_))')[0]
    mastercat.loc[has_observed_match, 'priority_name'] = 'coobserved'
    # \\ super prioritize Yifei's list
    highest_priority_pointings = table.Table(fits.getdata(f'../pointings/high_priority/{mfilt}_pointing_list_fall_2024.fits',1))
    is_hp = np.in1d(mastercat['object'], highest_priority_pointings['object'])
    mastercat.loc[is_hp, 'priority_name'] = 'mustget'

    
    is_queued, emptyhours = planner.plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mastercat,
                                       current_slot=slot,
                                       current_filter = mfilt,
                                       whichfield=None, priorities=priorities,
                                       verbose=verbose, 
                                       flag_emptyhours=True,                                       
                                       **kwargs )
    if len(emptyhours) > 0:
        print(f'Warning! Unqueued hours: {emptyhours}')
    return is_queued, emptyhours

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
    
def wrapra ( ra, pivot = 320. ):
    wrapped_ra = np.where(ra>pivot, ra-360, ra )
    return wrapped_ra
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser ( prog='sigjson_f2023b.py', description='Merian obs planner for 2023B')
    parser.add_argument ( '--year', '-Y', action='store', default=2023, help='year of observations', type=int )
    parser.add_argument ( '--month', '-M', action='store', help='month of observations', type=int )
    parser.add_argument ( '--day', '-D', action='store', help='day of observations', type=int )    
    parser.add_argument ( '--make_figure', action='store_true', help='[DEFUNCT] makes a figure that shows queued observations for tonight as queued.png')
    parser.add_argument ( '--telemetry_file', '-t', action='store', help='path to telemetry file.' )
    parser.add_argument ( '--copilot_file', '-c', default=f'{os.environ["HOME"]}/Downloads/db_merian.fits', action='store', help='path to copilot FITS db' )    
    parser.add_argument ( '--dryrun', action='store_true', help='do not save results' )
    parser.add_argument ( '--maxairmass', action='store', default=1.9, help='the maximum airmass at which exposures will be queued')
    parser.add_argument ( '--ignoresync', action='store_true', help='do not verify that telemetry and copilot files are synced')
    args = parser.parse_args ()
    
    if args.ignoresync:
        print('WARNING! You are not checking to make sure that the copilot and telemetry files are synced.')
    
    print(f'Running with airmass cut of {float(args.maxairmass)}')

    night_index = np.where([(night[:3] == [args.year, args.month, args.day]).all() for idx,night in enumerate(obsdates)])[0][0]
    night = obsdates[night_index]
        
    jsondir=f'../json/{args.year}{args.month:02d}{args.day:02d}/'
    if not os.path.exists(jsondir) and not args.dryrun:
        os.makedirs(jsondir)
        
    plan_args = [night[2], night[1], night[0], args.telemetry_file, args.copilot_file,]
    plan_kwargs = dict( slot = night[3],
                        mfilt=obsfilters[night_index], 
                        is_queued=None, 
                        pad_last_hour=True, 
                        maxairmass=float(args.maxairmass),
                        save=not args.dryrun, 
                        verbose=True,
                        ignore_synchronicity=args.ignoresync
    )
    if args.dryrun:
        print(f'[sigjson] We are doing a dry run of {args.year}/{args.month}/{args.day}')
        is_queued, emptyhours = plan_tomorrow(  *plan_args, **plan_kwargs)
    else:        
        with open(f'{jsondir}output.log', 'w') as sys.stdout:
            is_queued, emptyhours = plan_tomorrow( *plan_args, **plan_kwargs )
        sys.stdout = sys.__stdout__ # \\ gotta reset stdout 
        print(open(f'{jsondir}output.log','r').read())
    
    if args.make_figure:
        halpha_pointings, oiii_pointings =  our_pointings.load_fallfields()
        pointings_d = {'N540':oiii_pointings,'N708':halpha_pointings} 
        teff_min_d = {'N540':300, 'N708':300}
        fig, axarr = plt.subplots(2,1, figsize=(20,6))
        
        for idx, cfilter in enumerate(['N540','N708']):
            ax = axarr[idx]
            teff_min = teff_min_d[cfilter]
            pointings = pointings_d[cfilter]
            coo = observe.CopilotOutput ( args.copilot_file, pointings )
            completed = coo.identify_completed_pointings ( teff_min )
            to_obs = is_queued.loc[~is_queued['qstamp'].isna()]                        
            
            ax.scatter ( wrapra(pointings['RA']), pointings['dec'], 
                        facecolor='None', edgecolor='lightgrey', s=30**2, lw=1, label='planned' )
            ax.scatter ( wrapra(completed['racenter']), completed['deccenter'], s=30**2, color='#26b7f0', label='executed', alpha=0.2)
            pointings = pointings_d[cfilter]
            if cfilter == obsfilters[night_index]:
                print( f'We are observing {cfilter} ({idx})')                

                ax.scatter ( wrapra(pointings.reindex(to_obs.index)['RA']), pointings.reindex(to_obs.index)['dec'], 
                            facecolor='None', edgecolor='r', s=30**2, lw=1, label='queued [w. padding]' )
            ax.set_title ( cfilter )     
            ax.legend ()  
        
           
        for ax in axarr:
            ax.set_ylim(-2.7,5)            
            ax.set_xlabel('RA (deg)')
            ax.set_ylabel('Dec (deg)')     
            
        
        
        plt.tight_layout ()  
        
        figname = f'queued_{args.year}{args.month:02d}{args.day:02d}.png'
        if not args.dryrun:
            plt.savefig(jsondir + figname)
            print(f'Made queued plot at {jsondir}/{figname}')
        else:
            plt.show ()
    print('Finished planning!')