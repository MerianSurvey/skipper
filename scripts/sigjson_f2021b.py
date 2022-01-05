import sys
import pytz
import os
#import datetime
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
from astropy import coordinates
#from astropy import units as u
#from astropy.io import fits
from skipper import observe#, qa
import make_pointings

fmt = '%Y/%m/%d %I:%M %p'
et = pytz.timezone("America/New_York")
pt = pytz.timezone("America/Los_Angeles")
_BACKUP_FIELDS = ['SXDS','COSMOS']

######################### ==>
# \\ Filter and field assignments for F2021B
datelist_vvdsxmm_N540 = [(2021,11,ix) for ix in np.arange(24,31)]
datelist_vvdsxmm_N540 += [(2021,12,ix) for ix in np.arange(1,5)]
nightslot_vvdsxmm_N540 = [1 for ix in np.arange(24,31)]
nightslot_vvdsxmm_N540 += [1 for ix in np.arange(1,5)]
priorities_N540 = {'VVDSearly':0, 'VVDSlate':1, 'VVDS':1, 'XMM':2}    

datelist_vvdsxmm_n708 =  [(2021,9,ix) for ix in np.arange(10, 14)] 
#datelist_vvdsxmm_n708 += [(2021,11,ix) for ix in np.arange(24,31)]
nightslot_vvdsxmm_n708 = [2 for ix in np.arange(10,14)]
#nightslot_vvdsxmm_n708 += [1 for ix in np.arange(24,31)]
priorities_n708 = {'VVDSearly':0, 'VVDSlate':1,'VVDS':1, 'XMM':2}    

datelist_cosmosgama_N540 = [(2021,12,31), (2022,1,1)]
datelist_cosmosgama_N540 += [(2022,1,ix) for ix in np.arange(2,5)]
datelist_cosmosgama_N540 += [(2022,1,ix) for ix in np.arange(6,8)]
datelist_cosmosgama_N540 += [(2022,1,ix) for ix in np.arange(9,12)]
datelist_cosmosgama_N540 += [(2022,1,ix) for ix in np.arange(25, 32)]
nightslot_cosmosgama_N540 = [2,2]
nightslot_cosmosgama_N540 += [2 for ix in np.arange(2,5)]
nightslot_cosmosgama_N540 += [2 for ix in np.arange(6,8)]
nightslot_cosmosgama_N540 += [2 for ix in np.arange(9,12)]
nightslot_cosmosgama_N540 += [2 for ix in np.arange(25,32)]
priorities_cosmosgama_N540 = {'COSMOS':0, 'GAMA':1} 

# \\ total list    
datelist = datelist_vvdsxmm_N540 + datelist_vvdsxmm_n708 + datelist_cosmosgama_N540
nightslot = nightslot_vvdsxmm_N540 + nightslot_vvdsxmm_n708 + nightslot_cosmosgama_N540
nightslot_d = dict ( [(key,val) for key,val in zip(datelist,nightslot)])
priorities = {('VVDSXMM','N540'):priorities_N540, ('VVDSXMM','n708'):priorities_n708, ('COSMOSGAMA','N540'):priorities_cosmosgama_N540}
filter_l = len(datelist_vvdsxmm_N540) * ['N540'] + len(datelist_vvdsxmm_n708) *['n708'] + len(datelist_cosmosgama_N540)*['N540']
filter_d = dict ( [ (key,val) for key, val in zip(datelist, filter_l )])
field_l = len(datelist_vvdsxmm_N540) * ['VVDSXMM'] + len(datelist_vvdsxmm_n708) *['VVDSXMM'] + len(datelist_cosmosgama_N540)*['COSMOSGAMA']
field_d = dict ( [ (key,val) for key, val in zip(datelist, field_l )])
######################### <==

def whichfield ( year, month, day ):
    tpl = (year,month,day)
    field = field_d[tpl]
    mfilt = filter_d[tpl]
    nightslot = nightslot_d[tpl]
    return mfilt, field, nightslot

def load_mastercat ( filter_name, early_vvds=True ):
    vvds = pd.read_csv ( f'../pointings/vvds_{filter_name}.csv', index_col='object.1')
    if early_vvds:
        mask = vvds.RA < 345.
        vvds.loc[mask, 'object'] = vvds.loc[mask, 'object'].apply ( lambda x: x.split('_')[0] +'early' + '_' + '_'.join(x.split('_')[1:]) )
        vvds.loc[~mask, 'object'] = vvds.loc[~mask, 'object'].apply ( lambda x: x.split('_')[0] +'late' + '_' +  '_'.join(x.split('_')[1:]))
    xmm = pd.read_csv ( f'../pointings/xmm_{filter_name}.csv', index_col='object.1')

    mastercat = pd.concat([vvds, xmm])
    mastercat['proposer'] = 'Leauthaud'
    return mastercat

def load_mastercat_cosmos ( fname = '../pointings/cosmosgama_N540.csv' ):
    mastercat = pd.read_csv ( fname )
    mastercat = mastercat.set_index('object')
    mastercat['wait'] = "False"
    mastercat['proposer'] = 'Leauthaud'
    mastercat = mastercat.drop('object.1', axis=1)
    mastercat['object'] = mastercat.index
    return mastercat


def write_backupjson ():
    raise NotImplementedError

    for filt in ['g','r']:
        catalog_l, ocat, frd = make_pointings.build_backup (filter=filt)

        for ij in np.arange(10):
            fp = f'../json/backup_scripts/COSMOS_5minAGN_{filt}_{ij+1:02d}.json'
            print(fp)
            ocat.to_json ( catalog_l[ij], fp=fp, 
                           end_with_onemin=False )

def load_telemetry ( fname ):
    return pd.read_csv(fname)

def lazy_ephem ( day ):
    sunset, sunrise = ctio.get_sunriseset ( 2021, 3, day, alt=-10. )
    middle = sunset + 0.5*(sunrise-sunset)

    print(f'START: {sunset.strftime(fmt)} UTC; {sunset.astimezone(et).strftime(fmt)} ET')
    print(f'MID:   {middle.strftime(fmt)} UTC; {middle.astimezone(et).strftime(fmt)} ET')
    print(f'END:   {sunrise.strftime(fmt)} UTC; {sunrise.astimezone(et).strftime(fmt)} ET')
 
def predict_f2021b ( filter_name, datelist, nightslot_l, priorities=None, field='VVDSXMM', **kwargs ):
    if priorities is None:
        priorities = {}
        
    if field == 'VVDSXMM':
        mastercat = load_mastercat (filter_name)
    elif field == 'COSMOSGAMA':
        mastercat = load_mastercat_cosmos ()
    
    ocat = observe.ObsCatalog(comment='--', proposer='Leathaud', propid='2020B-0288', seqid='S2021B')

    # \\ build is_queued 
    is_queued = pd.DataFrame ( index=mastercat.index,
                               columns=['is_queued','qstamp','has_priority'])
    is_queued['is_queued'] = False

    # \\ Define the observatory site -- default is CTIO
    ctio = observe.ObservingSite ()
    
    
    for ix,date in enumerate(datelist):
        year, month, day = date        
        night_start, night_end = ctio.get_sunriseset ( year, month, day )
        if nightslot_l[ix] == 0:
            print('[predict] night slot: Full night')
            obs_start = night_start
            obs_end = night_end
        elif nightslot_l[ix] == 1:
            print('[predict] night slot: First half')
            obs_start = night_start
            obs_end = obs_start + 0.5*(night_end-obs_start)
        else:
            print('[predict] night slot: Second half')
            obs_start = night_start + 0.5*(night_end-night_start)
            obs_end = night_end 

        print(f"obsStart: {obs_start.astimezone(ctio.timezone).strftime(fmt)} Santiago")
        print(f"          {obs_start.astimezone(et).strftime(fmt)} ET")
        print(f"          {obs_start.strftime(fmt)} UTC")
        print(f"obsEnd:   {obs_end.astimezone(ctio.timezone).strftime(fmt)} Santiago")
        print(f"          {obs_end.astimezone(et).strftime(fmt)} ET")
        print(f"          {obs_end.strftime(fmt)} UTC")

        is_queued = ocat.plan_night ( obs_start, ctio, catalog=mastercat, obs_end=obs_end,
                                        is_queued=is_queued.copy(),
                                        save=False,
                                        maxairmass=1.5, object_priority=priorities,**kwargs )
        
    return is_queued

def write_backupjson ():
    from astropy import coordinates
    for filt in ['g','r']:
        center = coordinates.SkyCoord ("35.739030633438745 -4.7489828727193775", unit='deg')
        catalog_l, ocat, frd = make_pointings.build_backup (filter=filt, center=center, name='SXDS')

        for ij in np.arange(10):
            fp = f'../json/backup_scripts/SXDS_5minAGN_{filt}_{ij+1:02d}.json'
            print(fp)
            ocat.to_json ( catalog_l[ij], fp=fp, 
                           end_with_onemin=False )

def print_backupaltitudes (obs_start, obs_end, backup_fields=None):
    '''
    Plot backup field altitudes over the course of the (half) night
    '''
    _backup_centers = {'SXDS':coordinates.SkyCoord ("35.739030633438745 -4.7489828727193775", unit='deg'),
                      'COSMOS':coordinates.SkyCoord ("10h00m28.6s+02d12m21.0s") }
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
    if backup_fields is None:
        backup_fields = _BACKUP_FIELDS
    for name in backup_fields:
        for filt in ['g','r']:
            iu=0
            while True:
                fname = f'../json/backup_scripts/{name}_5minAGN_{filt}_{iu:02d}.json' 
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
                        
        
def plan_tomorrow ( day, month, year, tele_fname, copilot_fname, cut_at_contract=True, **kwargs ):
    '''
    Plan tomorrow in F2021B

    TODO: port this function to observe.py
    
    cut_at_contract: bool, default=True
        If True, the night end will occur, at latest, at 6:05 AM (Chilean work contract rule)
    '''
    tele = load_telemetry ( tele_fname )

    # \\ figure out which field and filter we're going to be observing in,
    # \\ TODO : manual override
    mfilt, field, slot = whichfield (year,month,day)
    
    # \\ QUICKFIX for 2021B TODO: fix in observe.py
    if mfilt == 'N708':
        skySB_0 = 21.
        teff_min = 200.
    elif mfilt == 'N540':
        skySB_0 = 22.1
        teff_min = 300.
    
    if field == 'VVDSXMM':
        mastercat = load_mastercat (mfilt)
    elif field == 'COSMOSGAMA':
        mastercat = load_mastercat_cosmos () # \\ only N540
    print(f"On {year}/{month}/{day}, we are observing {field} in {mfilt}")

    if slot==0:
        tag = 'full'
    elif slot==1:
        tag='first half'
    elif slot==2:
        tag = 'second half'
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
    else:
        print('[predict] night slot: Second half')
        obs_start = midpoint
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
                                     maxairmass=1.5, object_priority=priorities[(field,mfilt)],**kwargs )

    print_backupaltitudes (obs_start, obs_end )
    nextbackupscript ( tele )
    return is_queued_tmrw

kw_types = {'save':lambda x: x.lower == 'true'}

if __name__ == '__main__':
    print('='*27)
    print('== Merian Planner F2021B ==')
    print('='*27)


    if len(sys.argv)==1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('F2021B JSON file generation.')
        print('\n(usage) python sigjson_2021b.py [run day] [run month] [run year] [path/to/telemetry/output] [path/to/copilot/output]')
        print('  [run day] integer, the date of the beginning of the night')
        print('  [path/to/telemetry/output] path to the CSV containing the telemetry output from SISPI')
        print('(contact kadofong at princeton dot edu for help)\n')
    else:
        kwarg_keys = [ x[2:] for x in sys.argv[5::2] ]
        kwarg_cont = sys.argv[6::2]
        kwargs = dict(zip(kwarg_keys,kwarg_cont))
        for key in kwargs:
            kwargs[key] = kw_types[key](kwargs[key])
        
        print(f'DAY:       {sys.argv[1]}')
        print(f'MONTH:     {sys.argv[2]}')
        print(f'YEAR:      {sys.argv[3]}')
        print(f'TELEFILE:  {sys.argv[4]}')
        print(f'COPILOT:   {sys.argv[5]}')
        w = open('../resources/aart.txt','r').read()
        print(w)
        plan_tomorrow ( int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],**kwargs )
