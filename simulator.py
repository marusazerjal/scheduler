'''
Scheduler simulation over a few years.
Use one tiling file (it includes most of the stars from the priority list) and simulate observing nights through years. Include weather, etc.
'''

import numpy as np
import random
from datetime import date, timedelta
from os import listdir
from os.path import isfile, join
import copy
import simplejson as json
import matplotlib.pylab as plt

from astropy.time import Time
from astroplan.moon import moon_illumination

import params_simulator
reload(params_simulator)
import scheduler
reload(scheduler)


# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.
# It is better to change data_output_folder = 'data_output/' to something else (in the params.py file). Don't forget to add nearest_heighbours file.


def find_dates_to_observe():
    '''
    Determine nights to observe. Binary outcome: YES or NO.
    Include:
    (1) Moon phase
    (2) Overcast
    (3) Dome closed due to technical and other issues
    Poor conditions (but still acceptable for bright tiles) are determined later.
    '''
    
    # Fraction of cloudy days per month
    clouds_stats=np.loadtxt(params_simulator.params['weather_clouds_stats'], comments='#')
    clouds_stats=clouds_stats[:,1]
    
    DELTA = params_simulator.params['date_finish']-params_simulator.params['date_start']
    
    dates=[]
    for i in range(DELTA.days + 1):
        t=params_simulator.params['date_start'] + timedelta(days=i)
        ts=str(t)
        
        # Moon phase
        moon_phase=moon_illumination(Time(ts))
        if moon_phase>params_simulator.params['MIN_MOON_PHASE_TO_OBSERVE']:
            p1=True
        else:
            p1=False
        
        # Overcast
        ti=int(t.month)
        print '***', t
        print t.month
        print ti
        print
        if random.random()>clouds_stats[ti]: # TODO: include more detailed monthly weather statistics
            p2=True
        else:
            p2=False
        
        # Dome closed due to technical and other issues
        if random.random()>params_simulator.params['technical_issues_rate']:
            p3=True
        else:
            p3=False
        
        p=p1*p2*p3
        if p:
            dates.append(ts)

    print 'Number of nights expected in %d days between'%DELTA.days, params_simulator.params['date_start'], 'and', params_simulator.params['date_finish'], 'is %d.'%len(dates)
    return dates    

def simulate():
    '''
    Simulate observations.
    '''
    dates=find_dates_to_observe()
    
    s=scheduler.Scheduler()
    f=open(params_simulator.params['simulate_dates_file'], 'wb')
    for date in dates:
        print 'START NIGHT:', date
        f.write(date.replace('-', '')+'\n')
        
        # skip if already exists (was computed during the previous simulations)
        #~ if os.path.isdir("/home/el")):
            #~ continue
        
        s.observing_plan(date=date, remove_twilight=True)

        # Seeing
        #~ seeing=random.gauss(2.0, 1.0) # Maybe gauss is not the best distribution
    f.close()
    
    
if __name__ == "__main__":
    simulate()
