#!/usr/bin/env python
"""
Base functions to generate observing JSONs 
"""
__author__ = "Erin Kado-Fong"

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
    def __init__ ( self, site='CTIO', timezone=None):
        if site=='CTIO':
            self.site = coordinates.EarthLocation ( lat='-30d10m10.78s',
                                                    lon='-70d48m23.49s',
                                                    height=2241.*u.m )
            timezone='America/Santiago'

            
        else:
            self.site = site


            

    def get_utcoffset ( self, zone_name):
        from pytz import timezone
        
        
