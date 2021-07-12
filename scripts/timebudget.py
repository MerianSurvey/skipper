# \\ calculate how many exposures I might have left,
# \\ based on the time and sequence ... yes, lazy, I know

import datetime

def timetoendexp ( dt, timeleft, exptime=10. ):
    '''
    Time to end of exposure
    '''
    return dt + datetime.timedelta(seconds=timeleft)

