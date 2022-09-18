# \\ load catalog of master pointings
import re
import pandas as pd
from astropy import coordinates
from astropy import units as u

def _load_mastercat_cosmos ( fname = '../pointings/cosmosgama_N540.csv' ):
    mastercat = pd.read_csv ( fname )
    mastercat = mastercat.set_index('object')
    mastercat['wait'] = "False"
    mastercat['proposer'] = 'Leauthaud'
    mastercat = mastercat.drop('object.1', axis=1)
    mastercat['object'] = mastercat.index
    return mastercat


def load_springfields ():
    '''
    Loads Halpha and OIII catalogs for the springs fields
    '''
    # \\ load Shany's pointings (expanded for S21A HSC coverage)
    halpha_s2022a = pd.read_csv('../pointings/gama_2022A.csv', index_col='object.1')

    # \\ load Shany's old pointings
    oiii_pointings = _load_mastercat_cosmos ()
    halpha_pointings = _load_mastercat_cosmos ( '../pointings/S2021A.csv')
    
    # \\ add priorities to GAMA field
    is_high = (halpha_s2022a['dec'] > 1.5)
    is_early = (halpha_s2022a['RA'] < 160.)
    is_late = (halpha_s2022a['RA'] > 210.)

    halpha_s2022a['priority_name'] = 'GAMA'
    halpha_s2022a.loc[is_high&~is_early, 'priority_name'] = 'GAMAhigh'
    halpha_s2022a.loc[is_early, 'priority_name'] = 'GAMAearly'
    halpha_s2022a.loc[is_late, 'priority_name'] = 'GAMAlate'    
    
    # \\ copy & modify Halpha catalog to be OIII catalog
    get_catalog_objects = lambda x: x['object'].str.extract(r'(.*?(?=_))')[0]
    oiii_s2022a = halpha_s2022a.copy()
    oiii_s2022a['filter'] = 'N540'
    oiii_s2022a['object'] = [ xo.replace('N708','N540') for xo in oiii_s2022a['object'] ]
    oiii_s2022a.index = oiii_s2022a['object']
    oiii_s2022a['expTime'] = 900.

    # \\ add COSMOS pointings
    halpha_cosmos = halpha_pointings.loc[get_catalog_objects(halpha_pointings)=='COSMOS'].copy()
    halpha_cosmos['priority_name'] = 'COSMOS'
    halpha_s2022a = pd.concat([halpha_s2022a, halpha_cosmos], sort=False)                                   

    oiii_cosmos = oiii_pointings.loc[get_catalog_objects(oiii_pointings)=='COSMOS'].copy()
    oiii_cosmos['priority_name'] = 'COSMOS'
    oiii_s2022a = pd.concat([oiii_s2022a, oiii_cosmos], sort=False)    
    
    return halpha_s2022a, oiii_s2022a

def _load_fallfields_perfilter ( filter_name, early_vvds=True ):
    ''' 
    Loads Halpha or OIII catalog for the fall fields
    '''
    vvds = pd.read_csv ( f'../pointings/vvds_{filter_name}.csv', index_col='object.1')
    if early_vvds:
        mask = vvds.RA < 345.
        vvds.loc[mask, 'object'] = vvds.loc[mask, 'object'].apply ( lambda x: x.split('_')[0] +'early' + '_' + '_'.join(x.split('_')[1:]) )
        vvds.loc[~mask, 'object'] = vvds.loc[~mask, 'object'].apply ( lambda x: x.split('_')[0] +'late' + '_' +  '_'.join(x.split('_')[1:]))
    xmm = pd.read_csv ( f'../pointings/xmm_{filter_name}.csv', index_col='object.1')

    mastercat = pd.concat([vvds, xmm])
    mastercat['proposer'] = 'Leauthaud'
    return mastercat    

def load_fallfields_f2021 ( *args, **kwargs ):
    '''
    Loads Halpha and OIII catalogs for the fall fields -- the OLD pointings
    from F2021B. These do not include the expanded coverage from S21A
    '''
    hcat = _load_fallfields_perfilter ( 'N708', *args, **kwargs )
    ocat = _load_fallfields_perfilter ( 'N540', *args, **kwargs )
    return hcat, ocat

def construct_fall_superset ( newp, oldp ):
    '''
    Construct the superset of the new fall pointings & old fall pointings without
    duplicating objects.
    
    Also preserves old object names
    '''
    new_coords = coordinates.SkyCoord ( newp['RA'], newp['dec'], unit=('deg','deg'))
    old_coords = coordinates.SkyCoord ( oldp['RA'], oldp['dec'], unit=('deg','deg'))

    matchid, sep, _ = old_coords.match_to_catalog_sky ( new_coords )
    match_limit = 0.05 * u.arcsec

    oldnames_wmatch = oldp.loc[sep<=match_limit].index
    newnames_wmatch = newp.index[matchid[sep<=match_limit]]

    new_unique_pointings = newp.index.difference(newnames_wmatch)
    pointings_to_add = newp.reindex(new_unique_pointings)
    pointings_to_add.index = [ x.replace('XMM_VVDS','btwnXV') for x in pointings_to_add.index ]
    
    # \\ object.1 -> object
    pointings_to_add = pointings_to_add.rename ( {'object.1':'object'}, axis=1)
    pointings_to_add['object'] = pointings_to_add['object'].apply(lambda x:x.replace('XMM_VVDS','btwnXV') )
    
    superset = pd.concat([oldp, pointings_to_add])
    
    # \\ cut off early VVDS fields -< not for now
    superset = superset.loc[(superset['RA']<175.)|(superset['RA']>339.)]
    # \\ deprioritize middle where we have no DEC coverage
    is_mid = (superset['RA']>(360.-8.))|(superset['RA']<28.)
    is_mid |= (superset['RA']>(360-13.))&(superset['RA']<(360.-5.))&(superset['dec']>3.)
    superset.loc[is_mid, 'object'] = superset.loc[is_mid, 'object'].apply(lambda x: x.replace('btwnXV', 'newRAbtwnXV'))
    above_xmm = (superset['RA']>38.)&(superset['RA']<=40.)
    superset.loc[above_xmm, 'object'] = superset.loc[above_xmm, 'object'].apply(lambda x: x.replace('btwnXV', 'XMMhigh'))
    #assert superset.shape[0] == (newp.shape[0] + oldp.shape[0] - oldnames_wmatch.shape[0])
    return superset

def load_newfallpointings(csv):
    x = pd.read_csv(csv, index_col=0)
    x.loc[x['RA']<0.,'RA'] = x.loc[x['RA']<0.,'RA'] + 360.
    return x

def load_fallfields ( *args, **kwargs):
    '''
    Loads Halpha & OIII pointings for the fall fields -- these DO include
    new pointings from S21A, and we should use these catalogs
    as the fall field catalogs.    
    '''
    # a mapping hiccup: copilot output saves 'object' which has
    # the VVDS early/late designations.
    hcat, ocat = load_fallfields_f2021 ()
    oiii_pointings = construct_fall_superset ( load_newfallpointings('../pointings/xmm_vvds_N540.csv'),
                                               ocat )
    halpha_pointings = construct_fall_superset ( load_newfallpointings('../pointings/xmm_vvds_N708.csv'),
                                                 hcat )
    return halpha_pointings, oiii_pointings