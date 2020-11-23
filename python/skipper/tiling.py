#!/usr/bin/env python
"""
Calculate sky coverage for a decam tiling.
"""
__author__ = "Shany Danieli"


# from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')
from os.path import join,dirname,abspath
import pylab as plt
import numpy as np
import fitsio

from astrometry.util.fits import *
from astrometry.util.util import wcs_pv2sip_hdr, Tan
from astrometry.util.resample import *
from astrometry.libkd.spherematch import match_radec
from astrometry.util.plotutils import *
from astrometry.plot.plotstuff import *
from astrometry.util.util import anwcs_new_sip

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
                         size=(W,H), outfn='tile-%02i.pdf' % maxit)
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
    
            plt.imsave('tile-%02i.png' % it, cov, origin='lower', vmin=0, vmax=mx, cmap=antigray)
            #plt.imsave('tile-%02i.pdf' % it, cov, origin='lower', vmin=0, vmax=mx, cmap=antigray, format='pdf')
    
        if it in [30, 61, 90, 119]:
            from collections import Counter
            if it == 30:
                print(cov.tolist())
            quit()

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
            plt.savefig('hist-%02i.pdf' % it)



if __name__ == "__main__":

    decam()







