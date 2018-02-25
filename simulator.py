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
import pickle

from astropy.time import Time, TimeDelta
from astroplan.moon import moon_illumination
from astropy.coordinates import get_moon

import params_simulator
reload(params_simulator)
import scheduler
reload(scheduler)
import simulate_weather
reload(simulate_weather)

class Dictlist(dict):
	def __setitem__(self, key, value):
		try:
			self[key]
		except KeyError:
			super(Dictlist, self).__setitem__(key, [])
		self[key].append(value)

#~ print 'Start'

# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.
# It is better to change data_output_folder = 'data_output/' to something else (in the params.py file). Don't forget to add nearest_heighbours file.

# TODO
'''
Print ranking for each tile, or for all the stars observed.
Or maybe just print IDs and then read tiling file to get percentages of stars observed for each priority level.
'''

print 'Start simulation'

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
        
        # Moon position
        # Make smaller time steps. In hours.
        #~ utc = Time(Time.now(), format='iso', scale='utc')
        #~ moon = get_moon(utc)
        
        # Moon phase
        #~ moon_phase=moon_illumination(Time(ts))
        #~ if moon_phase>params_simulator.params['MIN_MOON_PHASE_TO_OBSERVE']:
            #~ p1=True
        #~ else:
            #~ p1=False
        p1=True
        
        # Overcast
        ti=int(t.month)-1 # start with 0 in array)
        if random.random()>clouds_stats[ti]: # TODO: include more detailed monthly weather statistics
            p2=True
        else:
            p2=False
        
        print p2
        
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
        print
        print 'START NIGHT:', date
        f.write(date.replace('-', '')+'\n')
        
        # skip if already exists (was computed during the previous simulations)
        #~ if os.path.isdir("/home/el")):
            #~ continue
        
        # Thin clouds. Observe only bright stars --> magnitude limit.
        mag_limit=True
        
        s.observing_plan(date=date, remove_twilight=True, bright_time=True)

        # Seeing
        #~ seeing=random.gauss(2.0, 1.0) # Maybe gauss is not the best distribution
    f.close()

def simulate2():
    s=scheduler.Scheduler()
    
    ra_current=None
    dec_current=None
    
    unique_targets=set()
    repeats=Dictlist()
    total_number_cumulative_unique = 0 # including repeated observations
    total_number_cumulative = 0 # including repeated observations

    # Print output
    f=open(params_simulator.params['simulator_statistics_output'], 'wb')
    f.write('# time; tile_id; Ntargets_in_this_tile; NUnique_targets_in_this_tile; NTotal_number_cumulative_unique; NTotal_number_cumulative; priority; weight; mag_max; json_filename \n')
    
    def time_per_tile(mag):
        if mag<12.5:
            exposure = 60.0
        elif 12.5<=mag<=14.5: # todo
            exposure = 5.0*60.0
        else:
            exposure = 10.0*60.0

        #~ calibration_time=60.0 # todo; done hourly
        reconfig_time=3.5*60.0 # between 2 and 5 minutes
        
        slew_time = 60.0 # 60 seconds. is this realistic?
        exposure += slew_time
        
        # Assumption: slewing during reconfig time
        
        total_time = np.max([exposure, reconfig_time])
        
        #~ total_time = exposure#+calibration_time# +reconfig_time: reconfig of the next tile is done during the exposure of the current tile
        
        return total_time

    def find_next_tile(ts=None, ra_current=None, dec_current=None, limiting_magnitude=None):
        best_tile, json_filename = s.next_tile(date=ts, ra_current=ra_current, dec_current=dec_current, bright_time=True, limiting_magnitude=limiting_magnitude, simulation_nickname=params_simulator.params['simulation_nickname'])
        if best_tile is not None:
            print best_tile
            
            # Unique targets
            science_targets=best_tile.TaipanTile.get_assigned_targets_science()
            ids=[x.idn for x in science_targets]
            unique_targets.update(ids)
            unique_tmp=set(ids).difference(unique_targets) # unique targets in the best_tile tile
            
            total_number_cumulative_unique += len(unique_tmp)
            total_number_cumulative += len(ids)
            
            # Repeated observations: just make list of all observations of all stars for now
            for x in ids:
                repeats[x]=ts
            
            f.write('%s; %d; %d; %d; %d; %d; %d; %f; %.1f; %s \n'%(ts, best_tile.TaipanTile.field_id, len(ids), len(unique_tmp), total_number_cumulative_unique, total_number_cumulative, best_tile.TaipanTile.priority, best_tile.weight*1000.0, best_tile.TaipanTile.mag_max, json_filename))
            
            mag=best_tile.TaipanTile.mag_max # TODO: check if this is upper or lower mag in the tile
            dt=time_per_tile(mag)
            ra_current=best_tile.TaipanTile.ra
            dec_current=best_tile.TaipanTile.dec
        else:
            dt=10*60.0 # jump for 10 minutes and hope weather gets better
        dt=TimeDelta(dt, format='sec')
        return ra_current, dec_current, dt
            
    i=0
    t=params_simulator.params['date_start']
    hourly_calibration_frames=TimeDelta(0, format='sec')
    while t < params_simulator.params['date_finish']:
        ts=str(t)[:-7]

        # Weather
        # TODO: sky illumination is included, so only night time is considered (when Sun below horizon). What about the Moon?
        if simulate_weather.is_weather_good(t.mjd):
            msg = 'Weather good.'
            ra_current, dec_current, dt = find_next_tile(ts=ts, ra_current=ra_current, dec_current=dec_current)
        
        elif simulate_weather.is_weather_acceptable_for_a_bright_tile(t.mjd):
            msg = 'Weather acceptable for a bright tile.'
            # TODO: determine limiting_magnitude
            ra_current, dec_current, dt = find_next_tile(ts=ts, ra_current=ra_current, dec_current=dec_current, limiting_magnitude=9.0)
        
        else:
            msg = 'Weather bad.'
            dt=TimeDelta(10*60.0, format='sec') # jump for 10 minutes

        t+=dt
        print t, 'dt', dt, params_simulator.params['date_finish'], msg
        
        hourly_calibration_frames+=dt
        if hourly_calibration_frames>TimeDelta(3600.0, format='sec'):
            hourly_calibration_frames=TimeDelta(0, format='sec')
            t+=TimeDelta(2*60.0, format='sec') # 2 minutes for hourly arcs and flats
        #~ i+=1
        #~ if i>10:
            #~ break
            
    f.close()
    
    with open('repeats_with_priorities.pkl', 'wb') as f:
        pickle.dump(repeats, f)
    
if __name__ == "__main__":
    #~ simulate()
    simulate2()
