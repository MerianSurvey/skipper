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
from astropy import coordinates
from astropy import units as u
from astropy.time import Time
from astropy.io import fits

class ObsCatalog (object):
    def __init__ (self,
                  comment='',
                  proposer='Leauthaud',
                  program='Merian',
                  propid='2020B-XXXX',
                  seqid='TEST',
                  count=1,
                  wait=False):
        self.comment = comment
        self.proposer = proposer
        self.program = program
        self.propid = propid
        self.wait = wait
        self.count = count
        self.seqid = seqid
        self.columns = np.asarray( ['comment','filter','seqtot','seqnum','seqid','expType',
                        'object','proposer','program','RA','propid','dec',
                        'expTime', 'count'] )
        self._singular = [True,False,True,False,True,False,False,
                           True,True,False,True,False,False,True]


    def build_catalog ( self, ra_l, dec_l, object_l, filter_l, expType_l, expTime_l):
        '''
        Generate pandas DataFrame that has our observing info
        '''
        self.seqtot = len(ra_l)
        catalog = pd.DataFrame ( index=np.arange(self.seqtot), columns=self.columns )
        for col in self.columns[self._singular]:
            catalog[col] = getattr(self, col)

        catalog['RA'] = ra_l
        catalog['dec'] = dec_l
        catalog['seqnum'] = catalog.index + 1
        catalog['object'] = object_l
        catalog['filter'] = filter_l
        catalog['expType'] = expType_l
        catalog['expTime'] = expTime_l

        self.catalog = catalog
        return catalog

    def write_jsonlog ( self, fp, logname='../json/json.log', user=None ):
        if user is None:
            try:
                user=os.environ['USER']
            except KeyError:
                raise ValueError ('No user was specified and cannot be read!')
        
        dtime = datetime.datetime.now().strftime('%Y/%m/%d %I:%M %p')
        with open(logname,'a') as ff:
            print(f'{fp} written on {dtime} by {user}', file=ff)
        
    def to_json (self, catalog=None, fp='../json/obsscript.json'):
        '''
        Format to JSON with small tweaks to enhance readability
        '''
        if catalog is None:
            catalog = self.catalog

        catalog.loc[:,'seqnum'] = np.arange(1,catalog.shape[0]+1)
        catalog.loc[:,'seqtot'] = catalog.shape[0]
            
        fx = {',':',\n', '[':'[\n', ']':'\n]', '{':'{\n', '}':'\n}' }
        json_str = catalog.to_json(orient='records')

        for key, repkey in fx.items():
            json_str = json_str.replace(key, repkey)

        open(fp, 'w').write(json_str)
        self.write_jsonlog ( fp )

    def observing_time ( self ):
        return self.catalog.expTime.sum() / 3600.

    def as_skycoord ( self, catalog=None):
        if catalog is None:
            catalog = self.catalog
        coords = coordinates.SkyCoord ( catalog['RA'],
                                        catalog['dec'], unit='deg')
        return coords

    def plan_night ( self, obs_start, obssite, catalog=None, maxairmass=1.3,
                     obs_end=None, 
                     is_queued=None, object_priority=None):
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
        '''
        if catalog is None:
            catalog = self.catalog
        dstr = obs_start.strftime('%Y%m%d')
        dpath = f'../json/{dstr}'
        if not os.path.exists(dpath):
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
                is_that_object = catalog['object']==key
                is_queued.loc[is_that_object, 'has_priority'] = val
        elif 'has_priority' not in is_queued.columns: # \\ don't overwrite if already there
            is_queued['has_priority'] = np.inf
            
        # \\ START planning the night
        for ix in range(len(alt_l[0])): # \\ for each hour,
            htime = alt_l[0][ix].obstime.datetime            
            hstr = htime.strftime('%Y%m%d_%H')

            cmass = pd.DataFrame(index=catalog.index)

            cmass['airmass'] = [ ai.secz[ix] for ai in alt_l]
            cmass['is_possible'] = True
            cmass.loc[cmass.airmass>maxairmass, 'is_possible'] = False
            cmass.loc[cmass.airmass<0,'is_possible'] = False
            cmass.loc[is_queued.is_queued, 'is_possible'] = False
            cmass['going_to_queue'] = False

            total_queued_time = 0.
            if ix == 0:
                total_available_time = (obsframe.obstime[ix] + 1.*u.hr - Time(obs_start)-0.5*u.hour).to(u.second).value

            elif ix==(len(alt_l[0])-1):
                total_available_time = (obsframe.obstime[ix]+1.*u.hr - Time(obs_end)-0.5*u.hour).to(u.second).value
            else:
                total_available_time = 3600.
            assert total_available_time <= 3600.1, f'{obs_end}, {total_available_time:.0f}s'
            if total_available_time < catalog.expTime.mean():
                print(f'({total_available_time:.0f}s) Not enough time for an exposure. Skipping...')
                continue
            print(f'==> {hstr}, {total_available_time}s available')
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

            hfile = catalog.loc[cmass.going_to_queue] #catalog.reindex(pidx).loc[g2q]
            if hfile.shape[0]==0:
                print('Nothing to queue!!')
            elif avail_queue_time > catalog['expTime'].mean():
                print(f'Cannot fill queue!! {avail_queue_time}, {catalog["expTime"].mean()}')
                warnings.warn('Queue unfilled at {hstr}')
            if hfile.shape[0]>0:
                self.to_json(hfile, fp=f'../json/{dstr}/{hstr}.json')

            is_queued.loc[cmass.index[cmass.going_to_queue], 'is_queued'] = True
            is_queued.loc[cmass.index[cmass.going_to_queue], 'qstamp'] = hstr
            is_queued.loc[cmass.index[cmass.going_to_queue], 'airmass'] = cmass.loc[cmass.going_to_queue, 'airmass']

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
                                                    height=2241.*u.m )
            self.timezone= pytz.timezone ( 'America/Santiago' )
        else:
            # // if not CTIO, trust the user to put in an EarthLocation
            # // and pytz.timezone ()
            self.site = site
            if type(timezone) == str:
                self.timezone = pytz.timezone(timezone)
            else:
                self.timezone = timezone

    def get_sunriseset ( self, year, month, day, alt=-10. ):
        '''
        DECam observing begins and ends at Sun altitude=-10 deg.
        '''
        utc_midnight = pytz.utc.localize ( datetime.datetime ( year, month, day, 0, 0 ) )
        utc_offset = int(self.get_utcoffset (utc_midnight))

        utc_start = pytz.utc.localize ( datetime.datetime ( year, month, day, 12-utc_offset, 0))
        utc_end = pytz.utc.localize ( datetime.datetime ( year, month, day+1, 12-utc_offset,0) )

        grid = np.arange(Time(utc_start), Time(utc_end),0.1*u.hour)
        sun_alt = []
        for ts in grid:
            sun_coord = coordinates.get_sun ( ts )
            obsframe = coordinates.AltAz ( obstime=ts, location=self.site )
            sun_alt.append( sun_coord.transform_to(obsframe).alt )

        sun_alt = np.asarray( [ sa.value for sa in sun_alt ] )
        observable = sun_alt <= alt
        obs_can_start, obs_must_end = grid[observable][[0,-1]]
        lcl = lambda x: pytz.utc.localize(x.to_datetime())
        return lcl(obs_can_start), lcl(obs_must_end)

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

        print(utc_end)
        #frame = np.arange ( -lim, lim+nstep/2., nstep) * u.hour
        frame = np.arange(0, (utc_end-utc_start).to_value(u.hour)+1,1)*u.hour
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
