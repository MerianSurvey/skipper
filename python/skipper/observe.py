#!/usr/bin/env python
"""
Base functions to generate observing JSONs 
"""
__author__ = "Erin Kado-Fong"

import os
import datetime
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
                  propid='2020B-XXXX' ):
        self.comment = comment
        self.proposer = proposer
        self.program = program
        self.propid = propid

        self.columns = np.asarray( ['comment','filter','seqtot','seqnum','expType',
                        'object','proposer','program','RA','propid','dec',
                        'expTime'] )
        self._singular = [True,False,True,False,False,False,
                           True,True,False,True,False,False]


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

    def observing_time ( self ):
        return self.catalog.expTime.sum() / 3600.

    def as_skycoord ( self, catalog=None):
        if catalog is None:
            catalog = self.catalog
        coords = coordinates.SkyCoord ( catalog['RA'],
                                        catalog['dec'], unit='deg')
        return coords

    def plan_night ( self, obstime, obssite, catalog=None, maxairmass=1.3,
                     obsstart=None, obsend=None,
                     is_queued=None, object_priority=None):
        '''
        Using obstime and obssite (CTIO), generate a plan from the night
        via airmass optimization.

        TODO: cut time at beginning at end based on ephem
        
        args:
        =====
        obstime (datetime.Datetime): date and central time for observing.
          We will plan the night for +/-6 hours around this central time by 
          default; if obsstart and obsend are supplied, we will plan
          between obsstart and obsend.
        obssite (observe.ObservingSite): Observatory object where we'll
          observer.
        catalog (observe.ObsCatalog): Catalog from which to generate
          observing plan. If None, use self.
        maxairmass (float): maximum airmass at which we will observe
        !!NOTIMPLEMENTED) obsstart (datetime.Datetime): date and time for observing to start.
        !!NOTIMPLEMENTED) obsend (datetime.Datetime): date and time for observing to end. 
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
        dstr = obstime.strftime('%Y%m%d')
        dpath = f'../json/{dstr}'
        if not os.path.exists(dpath):
            os.mkdir(dpath)

        # \\ Define Observing Frame from obstime
        obsframe = obssite.define_obsframe ( obstime)

        # \\ Get altitudes for the night
        alt_l = [ obssite.get_altitude(cc, obsframe) for cc in self.as_skycoord(catalog)]

        # \\ If we've not already got a queue, make one
        if is_queued is None:
            is_queued = pd.DataFrame (index=catalog.index,
                                      columns=['is_queued'])
            is_queued['is_queued'] = False

        # \\ If we've got object priorities, set those
        if object_priority is not None:
            is_queued['has_priority'] = np.inf
            for key,val in object_priority.items():
                is_that_object = catalog['object']==key
                is_queued.loc[is_that_object, 'has_priority'] = val
        else:
            is_queued['has_priority'] = np.inf

        # \\ START planning the night
        for ix in range(len(alt_l[0])): # \\ for each hour,
            htime = alt_l[0][ix].obstime.datetime            
            cmass = pd.DataFrame(index=catalog.index)

            cmass['airmass'] = [ ai.secz[ix] for ai in alt_l]
            cmass['is_possible'] = True
            cmass.loc[cmass.airmass>maxairmass, 'is_possible'] = False
            cmass.loc[cmass.airmass<0,'is_possible'] = False
            cmass.loc[is_queued.is_queued, 'is_possible'] = False
            cmass['going_to_queue'] = False
            total_queued_time = 0.
            for cprior in is_queued.has_priority.unique():                
                # \\ go through each object priority level
                avail_queue_time = 3600. - total_queued_time                
                is_this_priority = is_queued.reindex(cmass.index)['has_priority'] == cprior                
                targets = cmass.is_possible & is_this_priority
                pidx = cmass.loc[targets].sort_values('airmass').index
                going_to_queue = catalog.reindex(pidx)['expTime'].cumsum() <= avail_queue_time

                # \\ update remaining available queue time
                total_queued_time += catalog.reindex(pidx).loc[going_to_queue, 'expTime'].sum()
                print(f'{total_queued_time}s filled by priority={cprior} objects')
                cmass.loc[going_to_queue.loc[going_to_queue].index, 'going_to_queue'] = True

            #pidx = cmass.loc[cmass.is_possible].sort_values('airmass').index
            #g2q = catalog.reindex(pidx)['expTime'].cumsum () <= 3600.

            hfile = catalog.loc[cmass.going_to_queue] #catalog.reindex(pidx).loc[g2q]
            hstr = htime.strftime('%Y%m%d_%H')
            if hfile.shape[0]>0:
                self.to_json(hfile, fp=f'../json/{dstr}/{hstr}.json')

            is_queued.loc[cmass.index[cmass.going_to_queue], 'is_queued'] = True

        return catalog, is_queued
    

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

    def define_obsframe ( self, obs_datetime, nstep=1., lim=6. ):
        '''
        For an individual night of an observing run, generate the
        observing frame. The observing frame specifies observatory location
        + time for +/- lim hours about the fiducial obs_datetime in steps of nstep
        hours.
        '''
        #utcoffset = self.get_utcoffset ( obs_datetime )
        frame = np.arange ( -lim, lim+nstep/2., nstep) * u.hour
        timeframe = Time ( obs_datetime ) + frame
        obsframe = coordinates.AltAz ( obstime = timeframe, location=self.site)
        return obsframe#, utcoffset

    def get_altitude ( self, target_coord, obsframe ):
        '''
        For an observing frame and target coordinates, get altitude as a
        function of time
        '''
        alt = target_coord.transform_to ( obsframe )
        return alt
