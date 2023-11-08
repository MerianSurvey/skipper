import numpy as np
import matplotlib.pyplot as plt 

from ekfplot import plot as ek
from ekfplot import colors as ec

from skipper import observe

coo = observe.CopilotOutput('/Users/kadofong/Downloads/db_merian.fits', None )

fig, axarr = plt.subplots(1,2, figsize=(12,5))
exptimes = {'N708':600., 'N540':900.}
sb_bounds = {'N708':(19.5,21.5), 'N540':(20.5,22.5)}
assumed = {'N708':21., 'N540':22.1}
colors = {'N540': ec.ColorBase('#0d7ab5'), 'N708':ec.ColorBase('#b51b0d')}
for axindex, filt in enumerate(['N708','N540']):
    ax = axarr[axindex]
    ek.text ( 0.025, 0.975, filt, color=colors[filt].base, fontsize=20, ax=ax )
    
    is_thisfilter = coo.sidecar.apply( lambda x: filt in x['object'], axis=1 )
    teff_cut = exptimes[filt] / 3.
    is_goodexposure = is_thisfilter&(coo.sidecar['t_eff']>teff_cut)    
    bins = np.linspace(*sb_bounds[filt],30)
    for midx,mask in enumerate([is_thisfilter, is_goodexposure]):
        cc = colors[filt].modulate(0.3*midx).base
        
        med = coo.sidecar.loc[mask,'sky'].median()
        ax.hist ( coo.sidecar.loc[mask, 'sky'], bins=bins, color=cc, alpha=0.7, 
                 histtype=(midx==1) and 'step' or 'bar',
                 lw=3)

        ax.axvline ( coo.sidecar.loc[mask,'sky'].median(), color=cc, lw=3, ls='--' )
        if midx == 1:
            ek.text ( 0.025, 0.9, r'$\langle \mu_{\rm sky} \rangle=%.1f$ mag arcsec$^{-2}$' % med, 
                      ax=ax )
            ax.axvline(assumed[filt], color='k', lw=3, ls='--')

for ax in axarr:    
    ax.set_xlabel(r'Sky Surface Brightness [mag arcsec$^{-2}$]')
    ax.set_ylabel(r'Number of Exposures')
ek.text ( 0.025, 0.82, '''Filled histogram = all exposures.
Unfilled histogram = only exposures 
above min t_eff.''',
         ax=axarr[0],
         fontsize=8
        )
plt.savefig('/Users/kadofong/Downloads/SBdiagnostic.png')