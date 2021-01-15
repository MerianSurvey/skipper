#!/usr/bin/env python
"""
Create various tiling schemes for the Merian Survey.
"""
__author__ = "Shany Danieli"

import os
from os.path import join,dirname,abspath
import fitsio
import logging
import numpy as np


def get_datadir():
    return join(dirname(dirname(dirname(abspath(__file__)))),'data')

def get_datafile(filename):
    dirname = get_datadir()
    filepath = os.path.join(dirname,filename)

    if not os.path.exists(filepath):
        msg = "File does not exists: %s"%filepath
        logging.warn(msg)
    else:
        return filepath    

DECALS = fitsio.read(get_datafile('decam-tiles_obstatus.fits'))[:47616]

N = int(len(DECALS)/3)
DTYPE = [
    ('TILEID','>i4'),('PASS','>i2'),('RA','>f8'),('DEC','>f8'),
]

def decals_rotate(ra,dec,dx,dy):
    """Perform a Euler angle rotation of the celesitial coordinates and
    return the rotated position of ra,dec in the original
    coordinate system. dx,dy specify the Z and Y Euler rotation
    angles in decimal degrees, respectively.
    It seems that the origin 3 rotation angles used for DECaLS in
    'decam-tiles_obstatus.fits' were:
    dx,dy = (0.0,0.0), (-0.2917, 0.0833), and (-0.5861, 0.1333)
    For a 4th rotation, I've used
    dx,dy = (-0.8805,0.1833)
    Parameters:
    -----------
    ra : Input right ascension
    dec: Input declination
    dx : Rotation in x-dimension/Euler Z/ra (deg)
    dy : Rotation in y-dimension/Euler Y/dec (deg)
    Returns:
    --------
    ra, dec : Dithered ra,dec tuple
    """
    if float(dx) == 0.0 and float(dy) == 0.0: return ra,dec
    from astropy.modeling.rotations import EulerAngleRotation
    R = EulerAngleRotation(dx,dy,0,'zyx')
    ra1,dec1 = R(ra,dec)
    ra1 += 360 * (ra1 < 0)
    return ra1,dec1


def create_fields(base, offsets, dither):
    base = np.copy(base)
    fields =[]
    for i,off in enumerate(offsets):
        data = np.zeros(len(base),dtype=DTYPE)
        data['TILEID'] = base['TILEID']
        data['PASS'] = i+1
        ra_dither,dec_dither = dither(base['RA'],base['DEC'],off[0],off[1])
        data['RA']  = ra_dither
        data['DEC'] = dec_dither
        fields.append(data)
    return np.concatenate(fields)

def create_decals_fields():
    """ Replicate first 3 DECaLS tilings from
    decam-tiles_obstatus.fits and add a 4th tiling.
    """
    offsets = [(0.,0.),(-0.2917, 0.0833),
               (-0.5861, 0.1333),
               (-0.8805,0.1833)
    ]
    base = DECALS[:N]
    dither = decals_rotate
    return create_fields(base, offsets, dither)

if __name__ == "__main__":

    decals = create_decals_fields()
    datadir = join(dirname(dirname(dirname(abspath(__file__)))),'data')
    outfile = datadir+'/decam-tiles-decals-merian.fits.gz'
    fitsio.write(outfile,decals,clobber=True)





