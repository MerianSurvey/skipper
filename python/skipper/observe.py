#!/usr/bin/env python
"""
Base functions to generate observing JSONs 
"""
__author__ = "Erin Kado-Fong"

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

        return catalog
        
    def to_json (self, catalog, fp='../json/obsscript.json'):
        '''
        Format to JSON with small tweaks to enhance readability
        '''
        fx = {',':',\n', '[':'[\n', ']':'\n]', '{':'{\n', '}':'\n}' }
        json_str = catalog.to_json(orient='records')

        for key, repkey in fx.items():
            json_str = json_str.replace(key, repkey)

        open(fp, 'w').write(json_str)

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
        # // observation date. So here I'll set up that machinery.
        # // (We're actually observing during the daylight savings time
        # // switch in 2021A, so we need to make sure this is robust)
        return date.astimezone ( self.timezone ).utcoffset().total_seconds()/3600.

    def define_obsframe ( self, obs_datetime, nstep=0.5, lim=6. ):
        '''
        For an individual night of an observing run, generate the
        observing frame. The observing frame specifies observatory location
        + time for +/- lim hours about the fiducial obs_datetime in steps of nstep
        hours.
        '''
        utcoffset = self.get_utcoffset ( obs_datetime )
        frame = np.arange ( -lim, lim+nstep/2., nstep) * u.hour
        timeframe = Time ( obs_datetime ) + frame
        obsframe = coordinates.AltAz ( obstime = timeframe, location=self.site)
        return obsframe, utcoffset

    def get_altitude ( self, target_coord, obsframe ):
        '''
        For an observing frame and target coordinates, get altitude as a
        function of time
        '''
        alt = target_coord.transform_to ( obsframe )
        return alt

    def myownformatter ( self, dtobj ):
        '''
        There's a bug in matplotlib 2.0.0 that forcibly converts all datetime
        objects to UTC when plotting, so we need to get around that with this hack
        '''
        return f'{dtobj.hour:02d}:00'
    
    def plot_altitude ( self, alt ):
        fig, ax = plt.subplots(1,1,figsize=(8,6))

        ax.plot ( alt.obstime.datetime, alt.alt, lw=3 )
        ax.set_xlabel ("UTC")

        twiny = ax.twiny ()
        twiny.set_xlim ( ax.get_xlim() )
        twiny.set_xticks ( ax.get_xticks() )
        twiny.set_xticklabels ([ myownformatter(mdates.num2date(ct).astimezone(ctio.timezone))
                                for ct in ax.get_xticks()] )
        twiny.set_xlabel('local time')

        xformatter = mdates.DateFormatter('%H:%M')
        plt.gcf().axes[0].xaxis.set_major_formatter(xformatter)
        #plt.gcf().axes[1].xaxis.set_major_formatter(xformatter)
        plt.grid ()
