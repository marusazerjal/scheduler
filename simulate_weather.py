#python script that takes the weather data from HATSOUTH and tries to work out how much time was lost due to only high humidity.

import numpy as np
import astropy.time as time
import pdb
import pylab as pl

import params_simulator

N=30000

add_8yr = time.TimeDelta((3600.0*24.0)*(6.0*365.0+2.0*366.0), format='sec') # 2009 --> 2017
add_5yr = time.TimeDelta((3600.0*24.0)*(4.0*365.0+2.0*366.0), format='sec') # 2009 --> 2015

try:
    #~ clouddata=CLOUDDATA[:N]
    weather_dict=WEATHER_DICT
    dates=np.array([k for k in weather_dict.iterkeys()])
except NameError:
    print 'Reading weather data...'
    CLOUDDATA=np.loadtxt(params_simulator.params['weather_sso_database'],dtype='str',delimiter=',')
    #~ clouddata=CLOUDDATA[:N]
    clouddata=CLOUDDATA
    #~ WEATHER_DICT={(time.Time(x[0],format='iso')+add_5yr).value[:-7]: {'skydiff': float(x[2]), 'wind': float(x[3]), 'dew': float(x[4]), 'rain': float(x[5]), 'wet': float(x[6]), 'light': float(x[7]), 'wind_average': float(x[9]), 'wind_max': float(x[10]), 'humidity': float(x[11])} for x in clouddata} # no seconds
    #~ WEATHER_DICT={(time.Time(x[0],format='iso')+add_5yr).value[:-7]: {'skydiff': float(x[2]), 'light': float(x[7]), 'humidity': float(x[11])} for x in clouddata} # no seconds
    WEATHER_DICT={(time.Time(x[0],format='iso')+add_5yr).mjd: {'skydiff': float(x[2]), 'light': float(x[7]), 'humidity': float(x[11])} for x in clouddata} # no seconds
    weather_dict=WEATHER_DICT
    dates=np.array([k for k in weather_dict.iterkeys()])
    
    '''
    Maybe filter out daytime here. You do it once here and then you save time later. Maybe filter out rain and overcast weather as well.
    '''
    
    print 'Finished with reading the weather data.'

def find_weather_data(t):
    idx = (np.abs(dates-t)).argmin()
    k=dates[idx]
    diff=np.abs(k-t)
    if diff<0.042: # 1h [JD]
        return weather_dict[k]
    else:
        return None

def is_weather_good(t):
    wd=find_weather_data(t)
    if wd is None:
        return True # For simulation purposes
    elif wd['light']<400 and wd['skydiff']<-20 and wd['humidity']<90:
        return True
    else:
        return False
    # cond=(av_wind<12)&(max_wind<18)&(skydiff<-20)&(dew>-2)&(light<400)
    # add seeing

    
def is_weather_acceptable_for_a_bright_tile(t):
    wd=find_weather_data(t)
    if wd is None:
        return True # For simulation purposes
    elif wd['light']<400 and wd['skydiff']<-15 and wd['humidity']<90:
        return True
    else:
        return False      

def is_it_night(t):
    wd=find_weather_data(t)
    if wd is None:
        return True # For simulation purposes
    elif wd['light']>400:
        return False
    else:
        return True

def exclude_bad_weather_from_the_weather_database():
    for date, wd in weather_dict.iteritems():
        if wd['light']<100 and wd['skydiff']<-15 and wd['humidity']<90:
            # clear night
            
    

if __name__ == "__main__":
    t=(time.Time('2014-10-24 11:20', format='iso')).mjd
    
    #~ find_weather_data(t)
    
    print is_weather_good(t)
    #~ print is_weather_acceptable_for_a_bright_tile('2014-10-24 11:20')
    print is_weather_acceptable_for_a_bright_tile(t)
