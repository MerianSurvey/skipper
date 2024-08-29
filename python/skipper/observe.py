#!/usr/bin/env python
"""
Base functions to generate observing JSONs 
"""
__author__ = "Erin Kado-Fong"

import os
import datetime
import warnings
import numpy as np
import pandas as pd
import pytz
from astropy import coordinates, table
from astropy import units as u
from astropy.time import Time
from astropy.io import fits

from . import utils, qa

dim_d = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
fmt = '%Y/%m/%d %I:%M %p'

class ObsCatalog (object):
    def __init__ (self,
                  comment='',
                  proposer='Leauthaud',
                  program='Merian',
                  propid='2020B-XXXX',
                  seqid='TEST',
                  count=1,
                  wait="False"):
        self.comment = comment
        self.proposer = proposer
        self.program = program
        self.propid = propid
        self.wait = wait
        self.count = count
        self.seqid = seqid
        self.columns = np.asarray( ['comment','filter','seqtot','seqnum','seqid','expType',
                        'object','proposer','program','RA','propid','dec',
                        'expTime', 'count', 'wait'] )
        self._singular = [True,False,True,False,True,False,False,
                           True,True,False,True,False,False,True,True]

    def objnamer_coordinates ( self, cat_row ):
        '''
        Define an object name in the SDSS convention using its coordinates 
        '''
        coord = coordinates.SkyCoord(cat_row['RA'],cat_row['dec'],unit=u.deg)
        coorstr = coord.to_string('hmsdms',sep='').split()
        ra_truncated = utils.truncate(coorstr[0],2)
        dec_truncated = utils.truncate(coorstr[1],2)
        filter_name = cat_row['filter']
        field_name = cat_row['object']
        return f'{field_name}_J{ra_truncated}{dec_truncated}_{filter_name}'

    def objnamer_healpix ( self, cat_row ):
        import healpy as hp
        pixel_index = hp.ang2pix(80000, cat_row['RA'],
                                   cat_row['dec'], lonlat=True)
        filter_name = cat_row['filter']
        field_name = cat_row['object']
        return f'{field_name}_HP{pixel_index}_{filter_name}'
        
    def build_object_names (self, catalog, naming_scheme='coord'):
        obj_df = pd.Series(index=catalog.index, dtype=str)

        namer_d = {'healpix':self.objnamer_healpix,
                   'coord':self.objnamer_coordinates}
        for ix in range(catalog.shape[0]):
            row = catalog.iloc[ix]            
            obj_name = namer_d[naming_scheme] (row )
            obj_df.loc[catalog.index[ix]] = obj_name
        return obj_df

    def build_catalog ( self, ra_l, dec_l, object_l, filter_l, expType_l, expTime_l, build_object_names=True, **kwargs):
        '''
        Generate pandas DataFrame that has our observing info

        If build_object_names is True, then we will build the object names
        according to self.get_object_name ()
        '''
        self.seqtot = len(ra_l)
        catalog = pd.DataFrame ( index=np.arange(self.seqtot), columns=self.columns )
        for col in self.columns[self._singular]:
            catalog[col] = getattr(self, col)
            
        catalog['RA'] = ra_l
        catalog['dec'] = dec_l
        catalog['seqnum'] = catalog.index + 1
        catalog['filter'] = filter_l
        catalog['expType'] = expType_l
        catalog['expTime'] = expTime_l


        # \\ add OBJECT names
        catalog['object'] = object_l
        if build_object_names:
            catalog['object'] = self.build_object_names(catalog, **kwargs)

        self.catalog = catalog
        return catalog.set_index('object', drop=False)

    def write_jsonlog ( self, fp, logname='../json/json.log', user=None ):
        if user is None:
            try:
                user=os.environ['USER']
            except KeyError:
                raise ValueError ('No user was specified and cannot be read!')
        
        dtime = datetime.datetime.now().strftime('%Y/%m/%d %I:%M %p')
        with open(logname,'a') as ff:
            print(f'{fp} written on {dtime} by {user}', file=ff)
        
    def to_json (self, catalog=None, fp='../json/obsscript.json',
                 insert_onemin_exposures=True,
                 insert_random_ome=False,
                 slew_scale=8.*u.deg,
                 insert_checksky_exposures=False,
                 previous_position=None,
                 end_with_onemin=False,
                 verbose=True):
        '''
        Format to JSON with small tweaks to enhance readability

        insert_onemin_exposures:
          Adds a 1 minute exposure in between each science exposure (ODIN)

        insert_random_ome:
          Adds random 1 minute exposures in between science exposure of slew size < slew_scale
        insert_checksky_exposures:
          Adds 3 1 minute exposures at the start of the script to check sky brightness
        '''
        nonsci_count = 0
        if catalog is None:
            catalog = self.catalog
        if insert_onemin_exposures:
            if verbose:
                print('[to_json] Inserting 60s focus exposures')
            catalog = catalog.reset_index(drop=True)
            for name, row in catalog.iterrows():
                throw = row.copy()
                if name == 0:
                    prev_row = previous_position
                else:
                    prev_row = catalog.loc[name-1]

                if prev_row is not None:
                    coord_next = coordinates.SkyCoord(throw['RA'], throw['dec'], unit=u.deg)
                    coord_prev = coordinates.SkyCoord(prev_row['RA'], prev_row['dec'], unit=u.deg)
                    slew_size = coord_next.separation(coord_prev)
                    if slew_size > slew_scale:
                        insert_throw=True
                        if verbose:
                            print(f'[to_json] Big slew from {prev_row["object"]} to {throw["object"]}')
                    elif insert_random_ome:
                        pull = np.random.uniform(0.,5)
                        
                        if pull > 3.:
                            insert_throw=True
                            if verbose:
                                print(f'[to_json] Adding random OME before {throw["object"]}')
                    else:
                        insert_throw=False
                else:
                    insert_throw=True

                if insert_throw:
                    throw['expTime'] = 60.
                    throw['object'] =  f'1minexp_{nonsci_count:03d}'
                    throw['comment'] = 'OneMinuteFocusExposure'
                    throw.name = throw.name - 0.5
                    catalog.loc[throw.name] = throw
                    nonsci_count += 1

                if end_with_onemin and name >= catalog.index.max():
                    # \\ if last and we are supposed to end with a
                    # \\ one min exposure, pad there
                    # \\ can do gteq because IF we are doing this THEN
                    # \\ OME insertion must be active, so the index MUST be
                    # \\ float type
                    print(f'[to_json] Ending with OME')
                    throw['expTime'] = 60.
                    throw['object'] =  f'1minexp_{nonsci_count:03d}'
                    throw['comment'] = 'OneMinuteFocusExposure'
                    throw.name = throw.name + 0.5
                    catalog.loc[throw.name] = throw
                    nonsci_count += 1
                   

            catalog = catalog.sort_index().reset_index(drop=True)
        if insert_checksky_exposures:
            if verbose:
                print('[to_json] Inserting CheckSky 60s exposures')
            catalog = catalog.reset_index(drop=True)
            checksky = catalog.iloc[0].copy ()
            checksky['expTime'] = 60.
            checksky['comment'] = 'CheckSkyExposure'
            for ix in np.arange(-3,0,1):
                checksky['object'] = f'checksky_{nonsci_count:03d}'
                nonsci_count += 1
                catalog.loc[ix] = checksky
            catalog = catalog.sort_index().reset_index(drop=True)

        catalog.loc[:,'seqnum'] = np.arange(1,catalog.shape[0]+1)
        catalog.loc[:,'seqtot'] = catalog.shape[0]
            
        fx = {',':',\n', '[':'[\n', ']':'\n]', '{':'{\n', '}':'\n}' }
        json_str = catalog.to_json(orient='records')

        for key, repkey in fx.items():
            json_str = json_str.replace(key, repkey)

        with open(fp,'w') as infp:
            infp.write(json_str)
        
        #self.write_jsonlog ( fp ) # \\ this is annoying on merge and not useful in practice


    def observing_time ( self ):
        return self.catalog.expTime.sum() / 3600.

    def as_skycoord ( self, catalog=None):
        if catalog is None:
            catalog = self.catalog
        coords = coordinates.SkyCoord ( catalog['RA'],
                                        catalog['dec'], unit='deg')
        return coords

    def plan_night ( self, obs_start,
                     obssite,
                     catalog=None,
                     maxairmass=1.3,
                     obs_end=None, 
                     is_queued=None,
                     object_priority=None,
                     exclude_hour_indices=None,
                     save=True,
                     checksky_at_start=True,
                     pad_first_hour=False,
                     pad_last_hour=True,
                     prefix='',
                     pointingdb_fname=None,
                     verbose=True,
                     flag_emptyhours=False,
                     ):
        '''
        Using obstime and obssite (CTIO), generate a plan from the night
        via airmass optimization.

        TODO: cut time at beginning at end based on ephem
        
        args:
        =====
        obs_start (datetime.Datetime): date and time for observing to start. Can be in either
          local time or UTC, as long as time zone is specified.
        obssite (observe.ObservingSite): Observatory object where we'll
          observer.
        catalog (observe.ObsCatalog): Catalog from which to generate
          observing plan. If None, use self.
        maxairmass (float): maximum airmass at which we will observe
        obs_end (datetime.Datetime): date and time for observing to end. If no obs_end is
          specified, a plan will be generated for 6 hours.
        is_queued (pandas.DataFrame): list of pointings that have
          already been observed
        object_priority (dict-like): if given, a list of objects in the
          observing catalog that will take priority over the rest of the
          catalog. Other pointings will only be queued if the priority
          objects are all already queued OR above the maximum airmass
          needed for observation.
          object priority should be given as {OBJECT_NAME:OBJECT_PRIORITY},
          where 0 is highest priority.
        exclude_hour_indices (list-like): if given, will exclude hours by index 
        save (bool, default=True): if True, save output files as JSON observing scripts
        checksky_at_start (bool, default=True): if True, add 3 60s exposures to the start
          of the first script in order to check sky brightness
        '''
        if catalog is None:
            catalog = self.catalog
        
        if 'priority_name' not in catalog.columns:
            if verbose:
                print('[plan_night] Inferring priorities from object names')
                catalog_objects = catalog['object'].str.extract(r'(.*?(?=_))')[0]
        else:
            print('[plan_night] Using "priority_name" column as object priorities')
            catalog_objects = catalog['priority_name']
        
        if exclude_hour_indices is None:
            exclude_hour_indices = []     
            
        if flag_emptyhours:
            emptyhours = []   

        # \\ if time < 12, we've started past midnight.
        start_hour = obs_start.astimezone(obssite.timezone).hour
        if start_hour > 12:
            dstr = obs_start.astimezone(obssite.timezone).strftime('%Y%m%d')
        else:
            tstart = obs_start.astimezone(obssite.timezone) - datetime.timedelta(days=1.)
            dstr = tstart.strftime('%Y%m%d')
            
        dpath = f'../json/{dstr}'
        if (not os.path.exists(dpath)) and save:         
            os.mkdir(dpath)

        # \\ Define Observing Frame from obstime
        obsframe = obssite.define_obsframe ( obs_start=obs_start, obs_end=obs_end)

        # \\ Get altitudes for the night
        alt_l = [ obssite.get_altitude(cc, obsframe) for cc in self.as_skycoord(catalog)]

        # \\ If we've not already got a queue, make one
        if is_queued is None:
            is_queued = pd.DataFrame (index=catalog.index,
                                      columns=['is_queued'])
            is_queued['is_queued'] = False
            is_queued['qstamp'] = ''

        # \\ If we've got object priorities, set those
        if object_priority is not None: # \\ allow overwrite if priorities change
            is_queued['has_priority'] = np.inf
            for key,val in object_priority.items():                
                is_that_object = catalog_objects==key
                is_queued.loc[is_that_object, 'has_priority'] = val
        elif 'has_priority' not in is_queued.columns: # \\ don't overwrite if already there
            is_queued['has_priority'] = np.inf
            
        # \\ START planning the night
        if checksky_at_start: # \\ the first hour might not have any exposures so add this flag
            has_checkedsky = False
        else:
            has_checkedsky = True

        previous_position = None # \\ initialize as null; no starting focus exp
        nexcluded = 0
        while True:
            hr = len(alt_l[0]) - 1 - nexcluded
            if hr in exclude_hour_indices:
                nexcluded+=1
            else:
                break
        
        for ix in range(len(alt_l[0])): # \\ for each hour,
            if ix in exclude_hour_indices:
                continue
            htime = alt_l[0][ix].obstime.datetime     
            #print(f'[observe l348] Reference time is {htime}. Obs start is {obs_start}')
            
            cmaxairmass = maxairmass
            #print(f'[observe l351] {cmaxairmass}')       
            hstr = htime.strftime('%Y%m%d_%H')

            cmass = pd.DataFrame(index=catalog.index)

            cmass['airmass'] = [ ai.secz[ix] for ai in alt_l]
            cmass['is_possible'] = True            
            cmass.loc[cmass.airmass<0,'is_possible'] = False
            cmass.loc[is_queued.is_queued, 'is_possible'] = False
            print(f'[observe] (Minimum airmass available: {cmass.loc[cmass.is_possible, "airmass"].min()})')
            ntargets = (cmass.loc[cmass.is_possible].airmass<cmaxairmass).sum()
            #print(ntargets)
            while ntargets == 0:                
                ntargets = (cmass.loc[cmass.is_possible].airmass<cmaxairmass).sum()
                #print(ntargets)
                cmaxairmass += 0.1
                if cmaxairmass > 2.2:
                    print('Preparing to reobserve pointings:')
                    cmass.loc[(cmass['airmass']<cmaxairmass)&(cmass['airmass']>0.), 'is_possible'] = True
            if cmaxairmass > maxairmass:
                print (f'[observe] No targets found at {maxairmass}. Raising max airmass allowable to {cmaxairmass}...')
            cmass.loc[cmass.airmass>cmaxairmass, 'is_possible'] = False
            #print(cmass['is_possible'].sum())
            cmass['going_to_queue'] = False

            total_queued_time = 0.
            if ix == 0:                
                total_available_time = (obsframe.obstime[ix] + 1.*u.hr - Time(obs_start) - 0.5*u.hour).to(u.second).value
                if (3600. - total_available_time) < (catalog['expTime'].mean()*.1): # if we're within ten pct of an exposure just round up
                    total_available_time = 3600.

                if pad_first_hour and total_available_time % 600. > 300:
                    
                    # \\ if we've got time for at least half an exposure, put one int
                    print(f'\n[plan_night] padding the first hour script with an exposure even though we only have {total_available_time % 600.:.1f}s left')
                    total_available_time += 600. - total_available_time% 600.                    
            elif ix==(len(alt_l[0])-1-nexcluded):
                #print(obsframe.obstime[ix]+1.*u.hr)
                #print(Time(obs_end))
                if pad_last_hour:
                    pad_time = catalog['expTime'].mean() * 6                    
                else:
                    pad_time = 0.
                total_available_time = (Time(obs_end) - obsframe.obstime[ix]+1.*u.hr - 0.5*u.hour).to(u.second).value + pad_time
                #print(total_available_time)
                # \\ Add two extra exposures at the end just in case.
                # \\ We definitely don't want
                # \\ to run out of queued objects!
            else:
                total_available_time = 3600.
                
            if ix != (len(alt_l[0])-1-nexcluded): # \\ don't do this when we intentionally add time
                assert total_available_time <= 3600.1, f'{obs_end}, {total_available_time:.0f}s'
            if total_available_time < catalog.expTime.mean():
                print(f'({total_available_time:.0f}s) Not enough time for an exposure. Skipping...')
                continue
            print(f'==> {hstr}, {total_available_time}s available')
            if pad_last_hour and ix==(len(alt_l[0])-1):
                print(f'[plan_night] padding the last hour script by {pad_time}s') # \\ say this here to make 
                # \\ log more readable
                
            for cprior in np.sort(is_queued.has_priority.unique()):                
                # \\ go through each object priority level
                avail_queue_time = total_available_time - total_queued_time                
                is_this_priority = is_queued.reindex(cmass.index)['has_priority'] == cprior                
                targets = cmass.is_possible & is_this_priority
                pidx = cmass.loc[targets].sort_values('airmass').index
                going_to_queue = catalog.reindex(pidx)['expTime'].cumsum() <= avail_queue_time

                # \\ update remaining available queue time
                new_qtime = catalog.reindex(pidx).loc[going_to_queue, 'expTime'].sum()
                total_queued_time += new_qtime
                print(f'{new_qtime}s filled by priority={cprior} objects')
                cmass.loc[going_to_queue.loc[going_to_queue].index, 'going_to_queue'] = True

                amass = cmass.reindex(pidx).loc[going_to_queue, 'airmass']
                cmass.loc[going_to_queue.loc[going_to_queue].index, 'airmass'] = amass

            avail_queue_time = total_available_time - total_queued_time 
            #pidx = cmass.loc[cmass.is_possible].sort_values('airmass').index
            #g2q = catalog.reindex(pidx)['expTime'].cumsum () <= 3600.

            hfile = catalog.loc[cmass.going_to_queue]
            pindex = is_queued.reindex(hfile.index).sort_values('has_priority').index
            hfile = hfile.reindex(pindex) # \\ sort by priority so that we finish
            # \\ high priority targets before slewing
            if hfile.shape[0]==0:
                print('!!! Nothing to queue !!!')
                warnings.warn (f'Queue empty at {hstr}')
                if flag_emptyhours:
                    emptyhours.append( hstr )
            elif avail_queue_time > catalog['expTime'].mean():
                print(f'!!! Cannot fill queue !!! {avail_queue_time}, {catalog["expTime"].mean()}')
                warnings.warn(f'Queue unfilled at {hstr}')
                if flag_emptyhours:
                    emptyhours.append( hstr )
            if hfile.shape[0]>0:
                if save:
                    fname = f'../json/{dstr}/{prefix}{hstr}.json'
                    self.to_json(hfile, fp=fname,
                                 insert_checksky_exposures=not has_checkedsky,
                                 previous_position=previous_position)
                    has_checkedsky=True
                    qa.validate_json(fname, 
                                     ix >0 and htime or obs_start, 
                                     pointingdb_fname)
                previous_position = hfile.iloc[-1]

            is_queued.loc[cmass.index[cmass.going_to_queue], 'is_queued'] = True
            is_queued.loc[cmass.index[cmass.going_to_queue], 'qstamp'] = hstr
            is_queued.loc[cmass.index[cmass.going_to_queue], 'airmass'] = cmass.loc[cmass.going_to_queue, 'airmass']

        if flag_emptyhours:
            return is_queued, emptyhours
        return is_queued
    

class ObservingSite ( object ):
    '''
    [Determine airmass as a function of time for a given observing catalog
    and date.]
    Defines an observing site from which we may calculate observing conditions,
    given an target list and date.
    '''
    def __init__ ( self, site='CTIO', timezone=None):
        if site=='CTIO':
            self.site = coordinates.EarthLocation ( lat='-30d10m10.78s',
                                                    lon='-70d48m23.49s',
                                                    height=2200.*u.m )
            self.timezone= pytz.timezone ( 'America/Santiago' )
        elif site=='palomar':
            self.site = coordinates.EarthLocation ( lat='+33d21m22.7s',
                                                    lon='-116d51m53.6s',
                                                    height=1712.*u.m )
            self.timezone= pytz.timezone ( 'America/Los_Angeles' )            
        else:
            # // if not CTIO, trust the user to put in an EarthLocation
            # // and pytz.timezone ()
            self.site = site
            if type(timezone) == str:
                self.timezone = pytz.timezone(timezone)
            else:
                self.timezone = timezone

    def get_sunriseset ( self, year, month, day, alt=-14., cut_at_contract=False, contract_time=(6,5), return_track=False ):
        '''
        DECam observing begins and ends at Sun altitude=-14 deg.
        
        cut_at_contract: bool, default=False
            if True, cut at <contract_time> in accordance to Chilean labor laws
        contract_time: tuple, default=(6,5)
            hour,min in Chilean local time at which observations must cease
        '''
        print (f'[observe] Computing sunrise and sunset on {year}/{month}/{day} at altitude = {alt}' )
        utc_midnight = pytz.utc.localize ( datetime.datetime ( year, month, day, 0, 0 ) )
        utc_offset = int(self.get_utcoffset (utc_midnight))

        utc_start = pytz.utc.localize ( datetime.datetime ( year, month, day, 12-utc_offset, 0))
        
        if (day == dim_d[month]) & (month==12):
            print('happy new years')
            utc_end = pytz.utc.localize ( datetime.datetime ( year+1, 1, 1, 12-utc_offset,0) )
        elif day == dim_d[month]:
            utc_end = pytz.utc.localize ( datetime.datetime ( year, month+1, 1, 12-utc_offset,0) )
        else:
            utc_end = pytz.utc.localize ( datetime.datetime ( year, month, day+1, 12-utc_offset,0) )
        


        grid = np.arange(Time(utc_start), Time(utc_end),10.*u.min)
        fgrid = np.arange(Time(utc_start), Time(utc_end), 1.*u.min)
        sun_alt = []
        for ts in grid:
            sun_coord = coordinates.get_sun ( ts )
            obsframe = coordinates.AltAz ( obstime=ts, location=self.site )
            sun_alt.append( sun_coord.transform_to(obsframe).alt )

        sun_alt = np.asarray( [ sa.value for sa in sun_alt ] )
        if return_track:
            return grid, sun_alt
        
        fgrid_unix = np.asarray([ gg.unix for gg in fgrid ])
        grid_unix = np.asarray([ gg.unix for gg in grid ])
        
        sun_falt = np.interp(fgrid_unix, grid_unix, sun_alt)
        
        observable = sun_falt <= alt
        obs_can_start, obs_must_end = fgrid[observable][[0,-1]]
        
        lcl = lambda x: pytz.utc.localize(x.to_datetime())
        night_start = lcl(obs_can_start)
        night_end = lcl(obs_must_end)
        
        
        if cut_at_contract:            
            utcoff = self.get_utcoffset ( night_end )
            contract_end = pytz.utc.localize(datetime.datetime ( night_end.year, night_end.month, night_end.day, contract_time[0]-int(utcoff), contract_time[1],))
            mat_end = min ( night_end, contract_end )
            
            if contract_end < night_end:
                print ( 'True night end is:')                
                print(f"obsEnd:   {night_end.astimezone(self.timezone).strftime(fmt)} Santiago")        
                print(f"          {night_end.strftime(fmt)} UTC")
                                
                print ( 'Updated night end is:')
                print(f"obsEnd:   {contract_end.astimezone(self.timezone).strftime(fmt)} Santiago")        
                print(f"          {contract_end.strftime(fmt)} UTC")    
            return night_start, mat_end
        else:
            return night_start, night_end
    def track_moon ( self, start_time, end_time, alt_returntype='max' ):
        '''
        Over a starting and ending time, return the illumination (fractional) and 
        maximum altitude of the moon.
        '''
        from .astroplan_moon import moon_illumination
        
        time_start_time = Time(start_time)
        time_end_time = Time(end_time)
        
        grid = np.arange(time_start_time, time_end_time,10.*u.min)
        fgrid = np.arange(time_start_time, time_end_time, 1.*u.min)
        
        moon_alt = []
        for ts in grid:
            moon_coord = coordinates.get_moon ( ts )
            obsframe = coordinates.AltAz ( obstime=ts, location=self.site )
            moon_alt.append( moon_coord.transform_to(obsframe).alt )
        moon_alt = np.asarray ( [ ma.value for ma in moon_alt ])
        
        fgrid_unix = np.asarray([ gg.unix for gg in fgrid ])
        grid_unix = np.asarray([ gg.unix for gg in grid ])
        moon_falt = np.interp(fgrid_unix, grid_unix, moon_alt)        
        
        moon_altreport = getattr(np, alt_returntype)(moon_falt)
        moon_altreport = max(0., moon_altreport)
        
        moon_cillum = np.mean([moon_illumination (time_start_time), 
                               moon_illumination(time_end_time)])
        
        return moon_cillum, moon_altreport
        

    def get_utcoffset ( self, date ):
        '''
        From a datetime object, return the UTC offset at the observing site
        at that date/time.
        '''
        # // because of daylight savings time,
        # // we need to calculate the UTC offset *after* we know the
        # // observation date.
        # // WARNING ) pytz-2018 has an outdated handling of Chilean
        # // daylight savings time. Make sure your pytz version is up-to-date (2020).
        return date.astimezone ( self.timezone ).utcoffset().total_seconds()/3600.

    def define_obsframe ( self, obs_start=None, nstep=1., lim=6.,
                          obs_end=None ):
        '''
        For an individual night of an observing run, generate the
        observing frame. The observing frame specifies observatory location
        + time for +/- lim hours about the fiducial obs_datetime in steps of nstep
        hours.
        '''
        #utcoffset = self.get_utcoffset ( obs_datetime )

        utc_start = Time(obs_start) - obs_start.astimezone(pytz.utc).minute*u.minute - obs_start.astimezone(pytz.utc).second*u.second
        if obs_end is not None:
            utc_end = Time(obs_end) - obs_end.astimezone(pytz.utc).minute*u.minute -obs_end.astimezone(pytz.utc).second*u.second
        else:
            utc_end = utc_start + lim*u.hour

        
        #frame = np.arange ( -lim, lim+nstep/2., nstep) * u.hour
        #frame = np.arange(0, (utc_end-utc_start).to_value(u.hour)+1,1)*u.hour
        timeframe = np.arange(utc_start+0.5*u.hour, utc_end+1*u.hour, 1.*u.hour) #utc_start + frame
        #timeframe += 0.5*u.hour # calculate airmass in middle of hour
        obsframe = coordinates.AltAz ( obstime = timeframe, location=self.site)
        return obsframe

    def get_altitude ( self, target_coord, obsframe ):
        '''
        For an observing frame and target coordinates, get altitude as a
        function of time
        '''
        alt = target_coord.transform_to ( obsframe )
        return alt

class CopilotOutput ( object ):
    def __init__ ( self, filename, pointings,  transparency_0=1., seeing_0=1., skySB_0=21.):
        self.sidecar = self.load ( filename )
        self.sidecar['t_eff'] = self.t_effective ( transparency_0, seeing_0, skySB_0)
        self.merian_sidecar = self.sidecar.loc[np.in1d(self.sidecar['object'], pointings)]
        
    def load ( self, filename ):
        _sidecar = fits.getdata ( filename, 1 )
        sidecar = table.Table ( _sidecar ).to_pandas ()        
        return sidecar
    
    def t_effective ( self, transparency_0, seeing_0, skySB_0 ):
        t0 = ( self.sidecar['transparency'] / transparency_0 )**2
        t1 = ( self.sidecar['seeing']/seeing_0 )**-2
        t2 = 10.**( (self.sidecar['sky']-skySB_0)/2.5 )
        return t0*t1*t2*self.sidecar['exptime']
        
    def identify_completed_pointings ( self, min_teff ):
        merian_sidecar = self.merian_sidecar        
        sorter = merian_sidecar.sort_values('t_eff')['object'].duplicated(keep='last')
        merian_sidecar = merian_sidecar.loc[~sorter] # \\ remove frames that        
        # \\ have already been reobserved; only consider the highest t_eff exposure
        return merian_sidecar.loc[merian_sidecar['t_eff']>=min_teff]
    
    def flag_for_reobservation ( self, min_teff=200., ):
        '''
        Flag exposures that need reobservation
        '''
        merian_sidecar = self.merian_sidecar
        #return merian_sidecar.loc[merian_sidecar['t_eff']<min_teff, 'object'].values
        sorter = merian_sidecar.sort_values('t_eff')['object'].duplicated(keep='last')
        merian_sidecar = merian_sidecar.loc[~sorter] # \\ remove frames that
        # \\ have already been reobserved; only consider the highest t_eff exposure
        return merian_sidecar.loc[merian_sidecar['t_eff']<min_teff, 'object'].values
        
    def get_observed_catalog ( self, min_teff ):
        '''
        Return a copilot database that only has finished pointings (i.e. existing and do not need
        reobservation).
        
        Default min t_eff:
            Halpha : 200 s
            OIII : 300 s
        i.e. 1/3 t_exp
        '''
        #coo_o_spring.merian_sidecar[~np.in1d(coo_o_spring.merian_sidecar['object'], needs_reobservation)]
        needs_reobservation = self.flag_for_reobservation (min_teff)
        return self.merian_sidecar[~np.in1d(self.merian_sidecar['object'], needs_reobservation)]