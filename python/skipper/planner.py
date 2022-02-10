from asyncio.format_helpers import extract_stack
import os
import pytz
import datetime
import numpy as np
import pandas as pd
from astropy import coordinates
from . import observe

# \\ Define relevant timezones
fmt = '%Y/%m/%d %I:%M %p'
et = pytz.timezone("America/New_York")
pt = pytz.timezone("America/Los_Angeles")
ct = pytz.timezone("Asia/Chongqing")


_BACKUP_FIELDS = ['SXDS','COSMOS']
_backup_centers = {'SXDS':coordinates.SkyCoord ("35.739030633438745 -4.7489828727193775", unit='deg'),
                    'COSMOS':coordinates.SkyCoord ("10h00m28.6s+02d12m21.0s") }

def load_telemetry ( fname ):
    return pd.read_csv(fname)


def print_backupaltitudes (obs_start, obs_end, backup_fields=None):
    '''
    Plot backup field altitudes over the course of the (half) night
    '''

    if backup_fields is None:
        backup_fields = _BACKUP_FIELDS
    backup_centers = [ _backup_centers[name] for name in backup_fields ]
    
    ctio = observe.ObservingSite ()
    obsframe = ctio.define_obsframe ( obs_start=obs_start, obs_end=obs_end )
    alt_l = [ ctio.get_altitude(cc, obsframe) for cc in backup_centers]
    
    dtime = [ alt_l[0][ix].obstime.datetime for ix in range(len(alt_l[0]))]
    hd = 'time (UTC)\t\t'
    for iw in range(len(backup_centers)):
        hd = f'{hd}{backup_fields[iw]}\t'
    bfly = open('../resources/bfly.txt','r').read()
    print(bfly)
    print(hd)    
    for iv in range(len(dtime)):
        dt = dtime[iv]
        st = f'{dt.strftime(fmt)}\t'
        for iw in range(len(backup_centers)):            
            airmass = alt_l[iw].secz[iv]
            st = f'{st}{airmass:.2f}\t'
        print(st) 
 
def nextbackupscript ( tele, backup_fields=None ):
    '''
    Check which backup scripts have already been observed, and print out the next one
    '''
    if backup_fields is None:
        backup_fields = _BACKUP_FIELDS
    for name in backup_fields:
        for filt in ['g','r']:
            iu=0
            while True:
                fname = f'../json/backup_scripts/{name}/5min/{name}_5minAGN_{filt}_{iu:02d}.json' 
                if not os.path.exists ( fname ):
                    iu += 1
                    continue
                elif iu > 10:
                    raise ValueError ("No back-up scripts available! Need to regenerate")
                else:
                    json = pd.read_json ( fname )
                    has_observed = np.in1d(json['object'].iloc[1:], tele['object']).any ()
                    if has_observed:
                        iu+=1
                    else:
                        print('------')
                        print(f'Next script for {name} [{filt}] is {fname}')
                        break
                    
def plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mastercat, 
                   whichfield=None,
                   current_filter = None,
                   current_field = None,
                   current_slot = None, 
                   is_queued=None,                  
                   cut_at_contract=True,
                   priorities=None,
                   extra = None, # manual extra hours                    
                   **kwargs ):
    '''
    Plan "tomorrow" in the observing time frame (i.e. the next night that we will be observing)
    
    mastercat: pd.DataFrame
        The catalog of all pointings in the field(s) to be observed 
    cut_at_contract: bool, default=True
        If True, the night end will occur, at latest, at 6:05 AM (Chilean work contract rule)
    '''
    tele = load_telemetry ( tele_fname )

    # \\ figure out which field and filter we're going to be observing in,
    if whichfield is not None:
        mfilt, field, slot = whichfield (year,month,day)
    else: 
        mfilt = current_filter
        field = current_field
        slot = current_slot
        
    
    # \\ QUICKFIX for 2021B TODO: fix in observe.py
    if mfilt == 'N708':
        skySB_0 = 21.
        teff_min = 200.
    elif mfilt == 'N540':
        skySB_0 = 22.1
        teff_min = 300.
    else:
        raise ValueError (f"the filter {mfilt} is not recognized!")

    print(f"On {year}/{month}/{day}, we are observing {field} in {mfilt}")

    if slot==0:
        tag = 'full'
    elif slot==1:
        tag='first half'
    elif slot==2:
        tag = 'second half'
    elif slot==3:
        tag ='second 3/4'
    elif slot==4:
        tag = f'second half and the last {extra} hours of the first half'
        
    print(f'We are observing the {tag} of the night')
    #exp_exposures = tele.query('(exptime>599.)&(object!="G09")').shape[0]
    has_observed = np.in1d(mastercat['object'], tele['object'])
    
    # \\ also check for exposures that need to be reobserved
    coo = observe.CopilotOutput ( copilot_fname, pointings=mastercat,  skySB_0 = skySB_0 )
    reobs = coo.flag_for_reobservation ( min_teff = teff_min )
    needs_reobservation = np.in1d(mastercat['object'], reobs)
    print(f'{needs_reobservation.sum()} pointings in this catalog need reobservation!')
    has_observed = has_observed & ~needs_reobservation
    
    
    #assert has_observed.sum() == exp_exposures, "We have observed exposures that aren't in the master catalog?!"
    ocat = observe.ObsCatalog(comment='--', proposer='Leathaud', propid='2020B-0288', seqid='F2021B')

    if is_queued is None:
        # \\ build is_queued <- previously observed objects
        is_queued = pd.DataFrame ( index=mastercat.index,
                                columns=['is_queued','qstamp','has_priority'])
        is_queued['is_queued'] = False
        is_queued.loc[has_observed, 'is_queued'] = True

    # \\ Define the observatory site -- default is CTIO
    ctio = observe.ObservingSite ()
    night_start, night_end = ctio.get_sunriseset ( year, month, day, cut_at_contract=cut_at_contract )
    if cut_at_contract:    
        ten_start, ten_end = ctio.get_sunriseset ( year, month, day, alt=-10, cut_at_contract=False)
        midpoint = ten_start + 0.5*(night_end-ten_start)
    else:
        midpoint = night_start + 0.5*(night_end-night_start)
    if slot == 0:
        print('[predict] night slot: Full night')
        obs_start = night_start
        obs_end = night_end
    elif slot == 1:
        print('[predict] night slot: First half')
        obs_start = night_start
        obs_end = midpoint
    elif slot==2:
        print('[predict] night slot: Second half')
        obs_start = midpoint
        obs_end = night_end 
    elif slot==3:
        print('[predict] night slot: Second 3/4')
        obs_start = night_start + 0.5*(midpoint - night_start)
        obs_end = night_end
    elif slot==4:
        if extra == 1:
            s = ''
        else:
            s = 's'
        print (f'[predict] night slot: Second half and the last {extra} hour{s} of the first half')
        obs_start = midpoint - datetime.timedelta(hours=extra)
        obs_end = night_end
        

    print(f"obsStart: {obs_start.astimezone(ctio.timezone).strftime(fmt)} Santiago")
    print(f"          {obs_start.astimezone(et).strftime(fmt)} ET")
    print(f"          {obs_start.astimezone(pt).strftime(fmt)} PT")
    print(f"          {obs_start.strftime(fmt)} UTC")
    print(f"obsEnd:   {obs_end.astimezone(ctio.timezone).strftime(fmt)} Santiago")
    print(f"          {obs_end.astimezone(et).strftime(fmt)} ET")
    print(f"          {obs_end.astimezone(pt).strftime(fmt)} PT")
    print(f"          {obs_end.strftime(fmt)} UTC")
    
    moon_cillum, moon_altreport = ctio.track_moon ( obs_start, obs_end)
    print(f'Moon illumination is: {moon_cillum:.2f}')
    print(f'Moon max altitude during observation is: {moon_altreport:.2f}')

    is_queued_tmrw = ocat.plan_night ( obs_start, ctio, catalog=mastercat, obs_end=obs_end,
                                     is_queued=is_queued.copy(),
                                     maxairmass=1.5, 
                                     object_priority=priorities,
                                     **kwargs )

    print_backupaltitudes (obs_start, obs_end )
    nextbackupscript ( tele )
    return is_queued_tmrw