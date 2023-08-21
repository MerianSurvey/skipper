# qa.py
# Joseph Wick
#
# Validation for merian JSON observing scripts
# Ensures:
#   1. File is in proper JSON format
#   2. Observing target is visible and within threshold airmax at observing time

import json

import os

import astropy
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import pandas as pd

from datetime import datetime, timedelta
import pytz

exptime_d = {'N540':10.*u.minute, 'N708':15.*u.minute}

def tmp():
    print('Hello, world')

def validate_json(file, obs_start, obssite, csv=None, maxairmass=1.5, logname='../json/json.log'):
    '''
    Validation for JSON observing scripts
    Verifies:
        1. File is in proper JSON format
        2. Observing target is within threshold airmax at observing time

        args:
        ===========
        - file (string): path to json file to validate
        - obs_start (datetime.Datetime): date and time for observing to start. Can be in either
            local time or UTC, as long as time zone is specified.
        - obssite (observe.ObsCatalog): Observatory object where we'll be observing, if None defaults to Cerro Tololo
        - maxairmass (float): maximum airmass at which we will observe
        - csv (string): path to the csv file with all already observed objects
        ===========
    '''

    valid_json = True
    valid_vis = True

    if __checkJSON__(file, csv):
        valid_json = True
        valid_vis = 1
    else:
        valid_json = False
        valid_vis = -1

    if valid_json:
        if __checkVis__(file, obs_start, obssite, maxairmass, csv) == False:
            valid_vis=False



    if valid_json and valid_vis:
        print('QA COMPLETE: All tests passed')

        try:
            user=os.environ['USER']
        except KeyError:
            raise ValueError ('No user was specified and cannot be read!')

        dtime = datetime.now().strftime('%Y/%m/%d %I:%M %p')

        with open(logname,'a') as ff:
            print(f'{file} VALIDATED at {dtime} by {user}', file=ff)




def __checkJSON__(f, csv):
    '''
    Checks if a file is in valid json format, makes sure important keys are present, checks value format
    '''
    reqdKeys = ['filter', 'RA', 'dec', 'expTime']

    content = ''
    with open(f, 'r') as content_file:
        content = content_file.read()

    try:
        j = json.loads(content)
    except:
        print('JSON Warning: Error Loading JSON in file:', f)
        return False

    for group in j:
        # get present keys
        keys=[]
        for key in group:
            keys.append(key)
        for i,key in enumerate(reqdKeys):
            if key not in keys:
                print('JSON Warning: missing key "' + key + '" in file ' + f + ' block ' + str(i))
                return False
        # check that values are valid
        for key in keys:
            val = group[key]
            if not __validJSONVal__(val, key):
                print('JSON Warning: "' + val + '" is invalid value format for key "' + key + '" in file ' + f)
                return False
    return True

# validJSONVal
# return boolean
def __validJSONVal__(val, forKey):
    '''
    Checks if JSON value is acceptable for a given key:
        - makes sure there's no spaces
        - makes sure numerical values are numerical
    '''
    if hasattr(val, '__len__') and ' ' in val:
        return False

    if forKey == 'expTime' or forKey == 'RA' or forKey == 'dec': #must be float
        try:
            float(val)
        except ValueError:
            return False

    elif forKey == 'count': #must be integer
        try:
            int(val)
        except ValueError:
            return False

    elif forKey == 'focus': #float or list of floats
        for f in __listify__(val):
            try:
                float(val)
            except ValueError:
                if val!= 'None': return False

    return True

# Checks is airmass is within acceptable range
# using input start time and site
def __checkVis__(file, obs_start, obssite, maxairmass, csv):
    '''
    Checks airmass against upper limit
    '''

    # load into json object
    content = ''
    with open(file, 'r') as content_file:
        content = content_file.read()
    j = json.loads(content)
    print(f'[qa l155] Reference time: {obs_start}')

    #site = obssite #EarthLocation.of_site('Cerro Tololo')
    site=None
    if obssite == None:
        site = EarthLocation( lat='-30d10m10.78s', lon='-70d48m23.49s', height=2241.*u.m )
    else:
        site = obssite
    timezone= pytz.timezone ( 'America/Santiago' )
    utc_start = Time(obs_start) #- obs_start.astimezone(pytz.utc).minute*u.minute - obs_start.astimezone(pytz.utc).second*u.second
    print(f'[qa] UTC start time: {utc_start}')

    #already observed objects
    if csv is not None:
        done = pd.read_csv(csv)['object'].to_list()

    ok=True

    for i,exposure in enumerate(j):
        # ra/dec coordinates
        ra = exposure['RA']
        dec = exposure['dec']
        #mediumband = exposure['filter']
        exptime = exposure['expTime']
        ra_dec = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

        # get airmass
        utc_start += exptime * u.second #exptime_d[mediumband]
        obsframe = astropy.coordinates.AltAz( obstime = utc_start, location=site)
        wt = (ra_dec.transform_to(obsframe)).secz
        #print(wt)

        # check airmass against our threshold
        if wt > maxairmass or wt<0:
            print(f'Warning: Airmass of {wt:.2f} in exposure ' + str(i+1) + ' of file ' + file)
            ok=False

        # check object repeat
        if csv is not None:
            __checkObsdPrior__(exposure['object'], done, file)

    return ok

# CheckObsdPrior
# checks if a given object has already been observed
# ie is present in the csv
def __checkObsdPrior__(object, done, file):
    ignore = ['G09', '1minexp', 'checksky']

    for i in ignore:
        if i in object:
            return True

    if object in done and object not in ignore:
        m = 'Repeat Object: ' + object + ' in file ' + file
        raise PointingError(m)

# listify
# takes a string in form [a, b, c] and returns a list [a, b, c]
# if string is just a single value x, returns [x]
def __listify__(s):
    if '[' in s:
        return s[1:len(s)-1].split(', ')
    else:
        return [s]

class PointingError(Exception):
    '''Exception raised for a pointing object that has already been observec'''
    def __init__(self, message='pointing error'):
        self.message=message
        super().__init__(self.message)

class JSONError(Exception):
    '''Exception raised for incorrect JSON format'''

    def __init__(self, message='json error'):
        self.message = message
        super().__init__(self.message)
