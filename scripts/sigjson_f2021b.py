import sys
import pytz
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import table
from astropy import units as u
from astropy.io import fits
from skipper import observe, qa
import make_pointings

fmt = '%Y/%m/%d %I:%M %p'
et = pytz.timezone("America/New_York")

def load_mastercat ( filter_name ):
    vvds = pd.read_csv ( f'../pointings/vvds_{filter_name}.csv', index_col='object.1')
    xmm = pd.read_csv ( f'../pointings/xmm_{filter_name}.csv', index_col='object.1')

    mastercat = pd.concat([vvds, xmm])
    mastercat['proposer'] = 'Leauthaud'
    return mastercat

def load_mastercat_cosmos ( fname = '../pointings/cosmosgama_n536.csv' ):
    mastercat = pd.read_csv ( fname )
    mastercat = mastercat.set_index('object.1')
    mastercat['wait'] = "False"
    mastercat['proposer'] = 'Leauthaud'
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

def plan_tomorrow ( day, month, year, tele_fname, add_extracosmos=False,  **kwargs ):
    '''
    Plan tomorrow in F2021B

    TODO: port this function to observe.py
    '''
    mastercat = load_mastercat ()
    tele = load_telemetry ( tele_fname )

    exp_exposures = tele.query('(exptime>599.)&(object!="G09")').shape[0]
    has_observed = np.in1d(mastercat['object'], tele['object'])

    #assert has_observed.sum() == exp_exposures, "We have observed exposures that aren't in the master catalog?!"

    ocat = observe.ObsCatalog(comment='--', proposer='Leathaud', propid='2020B-0288', seqid='S2021A')

    # \\ build is_queued <- previously observed objects
    is_queued = pd.DataFrame ( index=mastercat.index,
                               columns=['is_queued','qstamp','has_priority'])
    is_queued['is_queued'] = False
    is_queued.loc[has_observed, 'is_queued'] = True

    # \\ Define the observatory site -- default is CTIO
    ctio = observe.ObservingSite ()
    priorities = {'COSMOS':0, 'GAMA':1}

    obs_start, twibeg = ctio.get_sunriseset ( year, month, day )
    obs_end = obs_start + 0.5*(twibeg-obs_start)

    print(f"obsStart: {obs_start.astimezone(ctio.timezone).strftime(fmt)} Santiago")
    print(f"          {obs_start.astimezone(et).strftime(fmt)} ET")
    print(f"          {obs_start.strftime(fmt)} UTC")
    print(f"obsEnd:   {obs_end.astimezone(ctio.timezone).strftime(fmt)} Santiago")
    print(f"          {obs_end.astimezone(et).strftime(fmt)} ET")
    print(f"          {obs_end.strftime(fmt)} UTC")

    is_queued_tmrw = ocat.plan_night ( obs_start, ctio, catalog=mastercat, obs_end=obs_end,
                                     is_queued=is_queued.copy(),
                                     maxairmass=1.5, object_priority=priorities,**kwargs )

    # \\ add extra COSMOS great seeing queue
    if add_extracosmos:
        print('')

        print('-'*31)
        print('-- COSMOS great seeing queue --')
        print('-'*31)
        print('\n(only trigger if seeing is <.75" for more than 30min prior to the hour)\n')
        extra_cosmos = load_extracosmos ()
        is_queued_ec2 = ocat.plan_night ( obs_start, ctio, catalog=extra_cosmos,
                                          obs_end=obs_end,
                                          checksky_at_start=False,
                                          maxairmass=1.2, prefix='extracosmos_',
                                          pointingdb_fname=tele_fname, **kwargs )
    
    return is_queued_tmrw

kw_types = {'save':lambda x: x.lower == 'true'}

if __name__ == '__main__':
    print('='*27)
    print('== Merian Planner S2021A ==')
    print('='*27)


    if len(sys.argv)==1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('S2021A JSON file generation.')
        print('\n(usage) python sigjson_2021a.py [run day] [path/to/telemetry/output]')
        print('  [run day] integer, the date of the beginning of the night')
        print('  [path/to/telemetry/output] path to the CSV containing the telemetry output from SISPI')
        print('(contact kadofong at princeton dot edu for help)\n')
    else:
        kwarg_keys = [ x[2:] for x in sys.argv[3::2] ]
        kwarg_cont = sys.argv[4::2]
        kwargs = dict(zip(kwarg_keys,kwarg_cont))
        for key in kwargs:
            kwargs[key] = kw_types[key](kwargs[key])

        plan_tomorrow ( int(sys.argv[1]), sys.argv[2],**kwargs )
