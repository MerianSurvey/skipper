#!/usr/bin/env python
"""
Calculate sky coverage for a decam tiling.
"""
__author__ = "Shany Danieli"


# from __future__ import print_function
import sys

if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
from os.path import join,dirname,abspath
import pylab as plt
import numpy as np
import fitsio
from astropy.coordinates import SkyCoord

try:
    import astrometry
except ModuleNotFoundError:
    from astrometry.util.fits import *
    from astrometry.util.util import wcs_pv2sip_hdr, Tan
    from astrometry.util.resample import *
    from astrometry.libkd.spherematch import match_radec
    from astrometry.util.plotutils import *
    from astrometry.plot.plotstuff import *
    from astrometry.util.util import anwcs_new_sip

from shapely import geometry
from shapely import affinity
from shapely.ops import unary_union
from shapely.geometry import Point
from descartes.patch import PolygonPatch


fig_dir = '/Users/shanydanieli/projects/merian/skipper/figures/skycoverage/'


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

def decam():
    
    plt.figure(figsize=(4,3))
    plt.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.15)
    
    datadir = join(dirname(dirname(dirname(abspath(__file__)))),'data')

    # T = fits_table(datadir+'/decam-tiles_obstatus.fits')
    T = fits_table(datadir+'/decam-tiles-decals-merian.fits')
    T.rename('pass', 'passnum')
    ra,dec = 0.933, 0.
    I,J,d = match_radec(T.ra, T.dec, ra, dec, radius_in_deg=5.) #2.8) 
    print(len(I), 'tiles near 0,0')
    T.cut(I)
    T.dist = d
    print('dists:', d)
    print('Passes:', T.passnum)

    F = fitsio.FITS(os.path.join(datadir+'/decam.wcs')) 
    wcs = []
    for i in range(1, len(F)):
        hdr = F[i].read_header()
        wcs.append(wcs_pv2sip_hdr(hdr, W=2046, H=4094))

    W,H = 5000, 5000
    pixsc = 4./3600.
    targetwcs = Tan(ra, dec, W/2.+0.5, H/2.+0.5, -pixsc, 0., 0., pixsc,
                    float(W), float(H))
    II = np.lexsort((T.dist, T.passnum))
    
    # This is for making the (vector) PDF format tiling images.
    for maxit in [0, 6, 30, 31, 37, 61, 62, 68, 90, 91, 92, 119]:
        #mx = { 1: 2, 2: 4, 3: 6 }[t.passnum]
        plot = Plotstuff(outformat='pdf', ra=ra, dec=dec, width=W*pixsc,
                         size=(W,H), outfn=fig_dir+'tile-%02i.pdf' % maxit)
        plot.color = 'white'
        plot.alpha = 1.
        plot.plot('fill')
    
        out = plot.outline
        out.fill = True
        out.stepsize = 1024.
        print('maxit: '+str(maxit))
        if maxit >100:
            print('YES!!!!!')
            plot.color = 'pink'
        plot.color = 'black'
        plot.alpha = 0.4
        plot.apply_settings()
    
        for it,t in enumerate(T[II]):
            print('Tile', it, 'pass', t.passnum)
            for w in wcs:
                w.set_crval((t.ra, t.dec))
                out.wcs = anwcs_new_sip(w)
                plot.plot('outline')
            if it == maxit:
                print('Writing', it)
                plot.write()
                break
    
    # And this is for PNG-format tiling images and histograms.
    cov = np.zeros((H,W), np.uint8)
    for it,t in enumerate(T[II]):
        print('Tile', it, 'pass', t.passnum)
        for w in wcs:
            w.set_crval((t.ra, t.dec))
            #print('WCS:', w)
            try:
                Yo,Xo,Yi,Xi,nil = resample_with_wcs(targetwcs, w)
            except:
                #import traceback
                #traceback.print_exc()
                continue
            cov[Yo,Xo] += 1

    
        if it in [0, 6, 30, 31, 37, 61, 62, 68, 90]:
            mx = { 1: 2, 2: 4, 3: 6, 4:8 }[t.passnum]

            # plt.clf()
            # plt.imshow(cov, interpolation='nearest', origin='lower', vmin=0, vmax=mx)
            # plt.colorbar()
            # plt.savefig('tile-%02i.png' % it) 
    
            plt.imsave(fig_dir+'tile-%02i.png' % it, cov, origin='lower', vmin=0, vmax=mx, cmap=antigray)
            #plt.imsave('tile-%02i.pdf' % it, cov, origin='lower', vmin=0, vmax=mx, cmap=antigray, format='pdf')
    
        if it in [30, 61, 90, 119]:
            from collections import Counter
            if it == 30:
                print(cov.tolist())

            print('Coverage counts:', Counter(cov.ravel()).most_common())
            bins = -0.5 + np.arange(8)
            plt.clf()
            # n,b,p = plt.hist(cov.ravel(), bins=bins, normed=True)
            n,b,p = plt.hist(cov.ravel(), bins=bins, density=True)
            #plt.hist(cov.ravel(), bins=bins, normed=True, cumulative=True, histtype='step')
            # Cumulative histogram from the right...
            xx,yy = [],[]
            for blo,bhi,ni in reversed(list(zip(bins, bins[1:], n))):
                nc = float(np.sum(cov.ravel() > blo)) / len(cov.ravel())
                yy.extend([nc,nc])
                xx.extend([bhi,blo])
                if ni > 0:
                    if nc != ni:
                        if nc > ni+0.03:
                            # If there's room, label the histogram bin above, else below
                            plt.text((blo+bhi)/2., ni, '%.1f \%%' % (100.*ni), ha='center', va='bottom', color='k')
                        else:
                            plt.text((blo+bhi)/2., ni-0.01, '%.1f \%%' % (100.*ni), ha='center', va='top', color='k')
                    plt.text((blo+bhi)/2., nc, '%.1f \%%' % (100.*nc), ha='center', va='bottom', color='k')
    
            plt.plot(xx, yy, 'k-')
    
            plt.xlim(bins.min(), bins.max())
            plt.ylim(0., 1.1)
            plt.xlabel('Number of exposures')
            plt.ylabel('Fraction of sky')
            #plt.title('DECaLS tiling, %i pass%s' % (t.passnum, t.passnum > 1 and 'es' or ''))
            #plt.savefig('hist-%02i.png' % it)
            plt.savefig(fig_dir+'hist-%02i.pdf' % it)


class FocusedRandomDither ( object ):
    '''
    Create a focused random dither pattern around a central
    point. Similar dither pattern to HSC UDeep, but using evenly 
    spaced grid in azimuth instead of fixed 5 pointings.
    '''
    def __init__ ( self, center,
                   offset_radius=0.08,
                   random_max=0.125,
                   start_at_center=True,
                   ndither=40,
                   fov_radius = (3.18/np.pi)**0.5,
                   ):
        '''
        center: [array-like] (RA, Dec) in degrees or SkyCoord object
        '''
        if type(center) == SkyCoord:
            self.center = (center.ra.deg, center.dec.deg )
        else:
            assert len(center)==2, "Expected coordinates of length 2: (RA,Dec)."
            self.center = center

        
        self.offset_radius = offset_radius
        self.random_max = random_max
        self.start_at_center = start_at_center
        self.ndither = ndither
        assert self.ndither >= 1, "Cannot have <1 exposure!"
        self.fov_radius = fov_radius

    def _make_grid ( self, ngridstep=300, extent=None ):
        if extent is None:
            extent = self.random_max * 20
            
        gra_a = np.linspace(self.center[0] - extent,
                            self.center[0]  + extent, ngridstep)
        gdec_a = np.linspace(self.center[1] - extent,
                             self.center[1] + extent, ngridstep)

        step = np.diff(gra_a)[0]

        pl_l = []
        for cra in gra_a:
            for cdec in gdec_a:
                pl_l.append ( Point ( cra, cdec ).buffer ( step ) )

        grid = np.zeros([self.ndither,ngridstep,ngridstep])  
        self.grid = grid
        self.gridpoints = pl_l
        self.grid_ra = gra_a
        self.grid_dec = gdec_a   

    def get_dcenter ( self, theta ):
        '''
        Find dither center for a given angle. Randomly draw offsets
        in RA and Dec from a uniform distribution [-random_max, random_max]
        '''
        dely = self.offset_radius * np.sin ( theta )
        delx = self.offset_radius * np.cos ( theta )

        pull = lambda: np.random.uniform ( -self.random_max,
                                              self.random_max )
        cra = self.center[0] + delx + pull ()
        cdec = self.center[1] + dely + pull ()
        return cra, cdec

    def circular_fov ( self, cra, cdec ):
        dfov = Point ( cra, cdec ).buffer ( self.fov_radius )
        return dfov

    def decam_fp ( self, cra, cdec, rotation=0. ):
        from skymap.instrument.decam import DECamFocalPlane
        
        decam = DECamFocalPlane ()
        decam_arr = decam.rotate ( cra, cdec )
        ccd_l = []
        for xi in range(decam_arr.shape[0]):
            #corners[:,0] += cra # \\ treating as Cartesian b.c.
            #corners[:,1] += cdec # \\ SMALL AREA!!
            ccd = geometry.Polygon ( decam_arr[xi] )
            #ccd = affinity.rotate ( ccd, rotation, 'center' )
            ccd_l.append(ccd)

        ccd_mpoly = unary_union ( ccd_l )
        ccd_mpoly = affinity.rotate(ccd_mpoly, rotation, 'center' )
        return ccd_mpoly
    def get_centers (self):
        theta_a = np.linspace(0, np.pi*2, self.ndither+1)[:-1]
        centers = np.zeros ( [self.ndither, 2] )
        for ii,ita in enumerate(theta_a):
            if self.start_at_center and ii==0:
                centers[ii,0] = self.center[0]
                centers[ii,1] = self.center[1]
            else:
                cra, cdec = self.get_dcenter ( ita )
                centers[ii,0] = cra
                centers[ii,1] = cdec
        return centers
        
    def compute_coverage ( self, target_area, footprint=None,
                           rotate=False ):
        if footprint is None:
            footprint = self.circular_fov
            
        #grid = self.grid.copy ()

        theta_a = np.linspace(0, np.pi*2, self.ndither+1)[:-1]
        if rotate:
            centers = np.zeros ( [self.ndither, 3] )
        else:
            centers = np.zeros ( [self.ndither, 2] )
            
        area_a = np.zeros(self.ndither)
        poly_l = []
        ii=0
        pl_l = self.gridpoints

        for ii,ita in enumerate(theta_a):
            if self.start_at_center and ii==0:
                centers[ii,0] = self.center[0]
                centers[ii,1] = self.center[1]
            else:
                cra, cdec = self.get_dcenter ( ita )
                centers[ii,0] = cra
                centers[ii,1] = cdec

            if rotate:
                centers[ii,2] = np.random.uniform(0,359)
                dfov = footprint ( *centers[ii] )
            else:
                dfov = footprint ( centers[ii,0], centers[ii,1] )
            
            iarea = target_area.intersection ( dfov )
            area_a[ii] = iarea.area
            poly_l.append(dfov)

            #aa = [ dfov.contains ( pl_l[ix] ) for ix in range(len(pl_l)) ]
            #grid[ii, np.asarray(aa).reshape(grid.shape[1:])] = 1

        #self.grid = grid
        #self.poly_l = poly_l
        #self.area_a = area_a
        #self.centers = centers
        return poly_l, area_a, centers

    def evaluate_coverage ( self ):
        total_area_covered = unary_union ( self.poly_l ).area
        cy = (self.grid[in_cosmos] > cnexp).sum(axis=1)/in_cosmos.sum()
        

if __name__ == "__main__":

    decam()







