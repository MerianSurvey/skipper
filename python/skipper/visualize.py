import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates as mdates


def myownformatter ( dtobj ):
    '''
    There's a bug in matplotlib 2.0.0 that forcibly converts all datetime
    objects to UTC when plotting, so we need to get around that with this hack
    '''
    return f'{dtobj.hour:02d}'



def plot_altitude ( alt, obssite, ax=None, format=False, **kwargs ):
    if ax is None:        
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        format=True

    ax.plot ( alt.obstime.datetime, alt.alt, lw=3, **kwargs )

    if format:
        ax.set_xlabel ("UTC")
        ax.set_ylim(0,90)
        twiny = ax.twiny ()
        twiny.set_xlim ( ax.get_xlim() )
        twiny.set_xticks ( ax.get_xticks() )
        twiny.set_xticklabels ([ myownformatter(mdates.num2date(ct).astimezone(obssite.timezone))
                                for ct in ax.get_xticks()] )
        twiny.set_xlabel('local time (La Serena)')

        xformatter = mdates.DateFormatter('%H')
        plt.gcf().axes[0].xaxis.set_major_formatter(xformatter)
        #plt.gcf().axes[1].xaxis.set_major_formatter(xformatter)
        plt.grid ()

        twinx = ax.twinx ()
        twinx.set_ylim ( ax.get_ylim() )

        mask = ax.get_yticks () > 0
        twinx.set_yticks ( ax.get_yticks()[mask] )
        twinx.set_yticklabels([ f'{np.sin(np.deg2rad(cy))**-1:.2f}'
                                for cy in ax.get_yticks()[mask] ])

        twinx.set_ylabel('airmass')
        ax.set_ylabel('altitude (deg)')
    return ax
