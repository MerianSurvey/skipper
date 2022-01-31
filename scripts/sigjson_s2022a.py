import sys
import numpy as np
import pandas as pd
from skipper import planner

import make_pointings

obsdir="/Users/kadofong/Google Drive/My Drive/MerianSurvey/observing/DataLogs/"
priorities = {'COSMOS':0,'GAMA':3, 'GAMAhigh':2, 'GAMAearly':1, 'GAMAlate':3}

######################### ==>
dates = [ (2022,2,idx) for idx in range(2,11)]
dates += [ (2022,2,28) ]
dates += [ (2022,3,idx) for idx in range(1,10)]
#dates = np.asarray(dates)

slots = [ 4 for idx in range(2,6) ] # first Feb nights, 2/28, and 3/1
slots += [ 2 for idx in range(6,13)]
slots += [ 3 for idx in range(2,10) ] # second 3/4 for rest of run
slots = np.asarray(slots)

filters = [ 'N540' for idx in range(2,16) ]
filters += [ 'N708' for idx in range(5,10) ]
filters = np.asarray(filters)


assert len(filters) == len(slots)
assert len(slots) == len(dates)

nightslot_d = dict ( [(key,val) for key,val in zip(dates,slots) ] )
filter_d = dict ( [ (key,val) for key, val in zip(dates,filters) ] )
######################### <==

def whichfield ( year, month, day ):
    tpl = (year,month,day)
    field = 'COSMOSGAMA'
    mfilt = filter_d[tpl]
    nightslot = nightslot_d[tpl]
    return mfilt, field, nightslot

def _load_mastercat_cosmos ( fname = '../pointings/cosmosgama_N540.csv' ):
    mastercat = pd.read_csv ( fname )
    mastercat = mastercat.set_index('object')
    mastercat['wait'] = "False"
    mastercat['proposer'] = 'Leauthaud'
    mastercat = mastercat.drop('object.1', axis=1)
    mastercat['object'] = mastercat.index
    return mastercat

def load_mastercat ():
    # \\ load Shany's pointings (expanded for S21A HSC coverage)
    halpha_s2022a = pd.read_csv('../pointings/gama_2022A.csv', index_col='object.1')

    # \\ load Shany's old pointings
    oiii_pointings = _load_mastercat_cosmos ()
    halpha_pointings = _load_mastercat_cosmos ( '../pointings/S2021A.csv')
    
    # \\ add priorities to GAMA field
    is_high = (halpha_s2022a['dec'] > 1.5)
    is_early = (halpha_s2022a['RA'] < 160.)
    is_late = (halpha_s2022a['RA'] > 210.)

    halpha_s2022a['priority_name'] = 'GAMA'
    halpha_s2022a.loc[is_high&~is_early, 'priority_name'] = 'GAMAhigh'
    halpha_s2022a.loc[is_early, 'priority_name'] = 'GAMAearly'
    halpha_s2022a.loc[is_late, 'priority_name'] = 'GAMAlate'    
    
    # \\ copy & modify Halpha catalog to be OIII catalog
    get_catalog_objects = lambda x: x['object'].str.extract(r'(.*?(?=_))')[0]
    oiii_s2022a = halpha_s2022a.copy()
    oiii_s2022a['filter'] = 'N540'
    oiii_s2022a['object'] = [ xo.replace('N708','N540') for xo in oiii_s2022a['object'] ]
    oiii_s2022a.index = oiii_s2022a['object']
    oiii_s2022a['expTime'] = 900.

    # \\ add COSMOS pointings
    halpha_cosmos = halpha_pointings.loc[get_catalog_objects(halpha_pointings)=='COSMOS'].copy()
    halpha_cosmos['priority_name'] = 'COSMOS'
    halpha_s2022a = pd.concat([halpha_s2022a, halpha_cosmos], sort=False)                                   

    oiii_cosmos = oiii_pointings.loc[get_catalog_objects(oiii_pointings)=='COSMOS'].copy()
    oiii_cosmos['priority_name'] = 'COSMOS'
    oiii_s2022a = pd.concat([oiii_s2022a, oiii_cosmos], sort=False)    
    
    return halpha_s2022a, oiii_s2022a

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

 
def plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mfilt=None, **kwargs ):        
    halpha_s2022a, oiii_s2022a = load_mastercat ()
    if mfilt is None:
        mfilt, _, slot = whichfield ( year, month, day )
        if slot == 4:
            extra = 1. # ONE extra hour from Soo
        else:
            extra= None
    
    if mfilt == 'N540':
        mastercat = oiii_s2022a
    elif mfilt == 'N708':
        mastercat = halpha_s2022a
    else:
        raise ValueError (f"Filter {mfilt} not recognized.")
    
    is_queued = planner.plan_tomorrow ( day, month, year, tele_fname, copilot_fname, mastercat, whichfield=whichfield,
                   priorities=priorities, extra=extra, **kwargs )
    return is_queued

                   
def predict_all ( tele_fname=None, copilot_fname=None ):
    # \\\ Remaining nights of F2021B run
    rdates = [(2022,1,idx) for idx in range(25,32) ]
    rdates = np.asarray(rdates)    
        
    # \\ Load pointings, observation logs 
    _, oiii_s2022a = load_mastercat ()   
    if tele_fname is None:
        tele_fname = f'{obsdir}/tele20220112.csv'
    if copilot_fname is None:
        copilot_fname = f'{obsdir}/db_merian.fits'
    
    # \\ Predict remaining nights in F2021B
    is_queued_oiii = None
    is_queued_halpha = None
    for date in rdates:
        is_queued_oiii = planner.plan_tomorrow ( date[2],date[1],date[0], tele_fname, copilot_fname, oiii_s2022a, 
                                                current_field='COSMOSGAMA',
                                                current_filter='N540',                                                
                                                current_slot = 2,
                                                priorities=priorities,
                                                is_queued=is_queued_oiii,
                                                pad_last_hour=False,
                                                save=False)
        print('\n<(-.-)>\n')
    
    for date in dates:
        mfilt, _, _ = whichfield ( date[0],date[1],date[2] )
        if mfilt == 'N708':
            is_queued_halpha = plan_tomorrow ( date[2], date[1], date[0], tele_fname, copilot_fname, 
                                              is_queued=is_queued_halpha, pad_last_hour=False, save=False )
        elif mfilt == 'N540':
            is_queued_oiii = plan_tomorrow ( date[2], date[1], date[0], tele_fname, copilot_fname, 
                                              is_queued=is_queued_oiii, pad_last_hour=False, save=False )
        print('\n<(-.-)>\n')
    
    return is_queued_oiii, is_queued_halpha
                
    

kw_types = {'save':lambda x: x.lower == 'true'}

if __name__ == '__main__':
    print('='*27)
    print('== Merian Planner S2022A ==')
    print('='*27)


    if len(sys.argv)==1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('S2022A JSON file generation.')
        print('\n(usage) python sigjson_s2022a.py [run day] [run month] [run year] [path/to/telemetry/output] [path/to/copilot/output]')
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
