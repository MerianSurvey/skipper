import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import coordinates
from skipper import observe,tiling,visualize

fmt = '%Y/%m/%d %I:%M %p'

def build_cosmos ():
    '''
    Build COSMOS dithering pattern from FocusedRandomDither
    '''
    center = coordinates.SkyCoord ("10h00m28.6s+02d12m21.0s")
    size =  (1.4, 1.4)
    np.random.seed(1003)
    edges = [ (center.ra.deg-size[0]/2., center.dec.deg-size[0]/2.),
              (center.ra.deg+size[0]/2., center.dec.deg-size[0]/2.),
              (center.ra.deg+size[0]/2., center.dec.deg+size[0]/2.),
              (center.ra.deg-size[0]/2., center.dec.deg+size[0]/2.) ]    

    frd = tiling.FocusedRandomDither (center, random_max=0.1,
                                      offset_radius=0.01, ndither=40)
    centers = frd.get_centers ()
    ocat = observe.ObsCatalog(comment='--', proposer='LeathaudGreene', 
                          propid='2020B-0288', seqid='S2021A')
    catalog = ocat.build_catalog(centers[:,0], centers[:,1],
                                 'COSMOS', 'N708', 'object', 10.*60)
    assert not catalog.object.duplicated().any()
    return catalog, ocat

def build_gama ():
    '''
    TO-DO) fill in real GAMA pointings
    '''
    ocat = observe.ObsCatalog(comment='--', proposer='LeathaudGreene',
                              propid='2020B-0288', seqid='S2021A')
    catalog = pd.read_csv("../json/gama.csv", index_col=0)
    catalog['filter'] = 'N708'
    catalog['object'] = ocat.build_object_names (catalog)
    return catalog, ocat

def plan_tomorrow ( mastercat, day, month=3,
                    year=2021, observed_objects=None ):
    cosmos,ocat = build_cosmos ()
    gama,_ = build_gama ()
    
    mastercat = pd.concat([gama,cosmos])
    mastercat.index = mastercat['object']
    assert not mastercat.object.duplicated().any()

    # \\ Define the observatory site -- default is CTIO
    ctio = observe.ObservingSite ()
    priorities = {'COSMOS':0, 'GAMA':1, 'earlyGAMA':2}

    # \\ We need to deprioritize early targets in order to preserve
    # \\ them for the end of the run
    obs_start, twibeg = ctio.get_sunriseset ( 2021, 3, 17 )
    obs_end = obs_start + 0.5*(twibeg-obs_start)

    obsframe = ctio.define_obsframe ( obs_start=obs_start, obs_end=obs_end )
    alt_l = [ ctio.get_altitude(cc, obsframe) for \
                  cc in ocat.as_skycoord(mastercat)] 

    airmass = np.array([ai.secz[0] for ai in alt_l])
    mask =(airmass>0)&(airmass<1.41)&(mastercat['object']!='COSMOS')
    print (f'Deprioritizing N={mask.sum()} targets')
    mastercat.loc[mask, 'object'] = mastercat.loc[mask,'object'].apply(lambda x: 'early' + x )


    

def predict_s2021a ():
    cosmos,ocat = build_cosmos ()
    gama,_ = build_gama ()
    
    mastercat = pd.concat([gama,cosmos])
    mastercat.index = mastercat['object']
    assert not mastercat.object.duplicated().any()

    # \\ Define the observatory site -- default is CTIO
    ctio = observe.ObservingSite ()
    priorities = {'COSMOS':0, 'GAMA':1, 'earlyGAMA':2}

    # \\ We need to deprioritize early targets in order to preserve
    # \\ them for the end of the run
    obs_start, twibeg = ctio.get_sunriseset ( 2021, 3, 17 )
    obs_end = obs_start + 0.5*(twibeg-obs_start)

    obsframe = ctio.define_obsframe ( obs_start=obs_start, obs_end=obs_end )
    alt_l = [ ctio.get_altitude(cc, obsframe) for \
                  cc in ocat.as_skycoord(mastercat)] 

    airmass = np.array([ai.secz[0] for ai in alt_l])
    mask =(airmass>0)&(airmass<1.41)&(mastercat['object']!='COSMOS')
    print (f'Deprioritizing N={mask.sum()} targets')
    mastercat.loc[mask, 'object'] = mastercat.loc[mask,'object'].apply(lambda x: 'early' + x )

    # \\ Nights for S2021A
    nights = np.arange(5,18,1)
    nights = np.concatenate([nights[:7], nights[8:]])

    obs_start, twibeg = ctio.get_sunriseset ( 2021, 3, nights[0] )
    obs_end = obs_start + 0.5*(twibeg-obs_start)

    print(obs_start.astimezone(ctio.timezone).strftime(fmt))
    print(obs_end.astimezone(ctio.timezone).strftime(fmt))

    is_queued = ocat.plan_night ( obs_start, ctio, catalog=mastercat,
                                  obs_end=obs_end,
                                  maxairmass=1.5,
                                  object_priority=priorities )

    for night in nights[1:]:
        obs_start, twibeg = ctio.get_sunriseset ( 2021, 3, night)
        obs_end = obs_start + 0.5*(twibeg-obs_start)

        print(obs_start.astimezone(ctio.timezone).strftime(fmt))
        print(obs_end.astimezone(ctio.timezone).strftime(fmt))

        is_queued = ocat.plan_night ( obs_start, ctio,
                                      catalog=mastercat, obs_end=obs_end,
                                      maxairmass=1.5, is_queued=is_queued,
                                      save=False)
    return mastercat, is_queued
