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
from astropy.coordinates import get_sun#get_moon
from astropy.coordinates import AltAz
from astroplan import Observer

import params
reload(params)
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


# Count tiles. When you run out of all the tiles, terminate the simulator.
# time with no observing: ruined by bad weather.

# Remote twilight time from time available to observe.

print 'Start simulation'

class Simulator():
    def __init__(self):
        self.scheduler = scheduler.Scheduler()
        self.observatory = Observer.at_site("Anglo-Australian Observatory") # TODO: enter LAT and LON coordinates
        
        self.ra_current=None
        self.dec_current=None

        self.number_of_tiles = self.scheduler.number_of_all_tiles() # total number of all the tiles in this tiling run
        self.number_of_tiles_observed=0

        self.unique_targets=set()
        self.repeats=Dictlist()
        self.total_number_cumulative_unique = 0 # EXcluding repeated observations
        self.total_number_cumulative = 0 # INcluding repeated observations
    
        self.number_of_tiles_observed = 0
        self.number_of_all_tiles = scheduler.number_of_all_tiles()

        self.dt = TimeDelta(0, format='sec')
        self.time_with_no_observing = TimeDelta(0, format='sec')

        # print output
        self.f = open(params_simulator.params['simulator_statistics_output'], 'wb')
        self.f.write('# time; tile_id; Ntargets_in_this_tile; Nunique_targets_in_this_tile; Ntotal_number_cumulative_unique; Ntotal_number_cumulative; priority; weight; mag_max; json_filename \n')
        
        # Print calibration files (times)
        self.fc = open(params_simulator.params['simulator_statistics_output_calibration'], 'wb')

    def __finish__(self):
        # AT THE END
        self.f.close()
        self.fc.close()


    def time_per_tile(self, mag):
        #~ if mag<12.49:
            #~ exposure = 60.0
        #~ elif 12.49<mag<14.49: # todo
            #~ exposure = 5.0*60.0
        #~ else:
            #~ exposure = 10.0*60.0

        exposure=params_simulator.params['exposure_time'][mag]

        #~ calibration_time=60.0 # todo; done hourly
        reconfig_time=params_simulator.params['reconfig_time'] # between 2 and 5 minutes
        
        # We slew during reconfig time
        #~ slew_time = 90.0 # 60 seconds. is this realistic?
        #~ exposure += slew_time
        
        # Assumption: slewing during reconfig time
        
        #~ total_time = np.max([exposure, reconfig_time]) # wrong. we dont have two plates
        total_time = exposure + reconfig_time
        
        #~ total_time = exposure#+calibration_time# +reconfig_time: reconfig of the next tile is done during the exposure of the current tile
        
        return total_time

    def find_next_tile(self, ts=None, limiting_magnitude=None):
        '''
        Calls scheduler to find the next tile.
        And statistics.
        '''
        
        # Next tile
        best_tile, json_filename = self.scheduler.next_tile(date=ts, ra_current=self.ra_current, dec_current=self.dec_current, bright_time=True, limiting_magnitude=limiting_magnitude)
        
        
        # Statistics and technical stuff
        if best_tile is not None:
            print best_tile
            self.number_of_tiles_observed+=1
            self.print_statistics(best_tile=best_tile, json_filename=json_filename, ts=ts)

            self.ra_current=best_tile.TaipanTile.ra
            self.dec_current=best_tile.TaipanTile.dec  
            
            mag=best_tile.TaipanTile.mag_max # TODO: check if this is upper or lower mag in the tile
            dt=self.time_per_tile(mag)
            self.time_with_no_observing = TimeDelta(0, format='sec')
            self.bad_weather_duration = TimeDelta(0, format='sec')

        else:
            dt=10*60.0 # jump for 10 minutes and hope weather gets better
            self.time_with_no_observing += TimeDelta(10.0*60.0, format='sec')
            print 'time_with_no_observing:', self.time_with_no_observing

        dt=TimeDelta(dt, format='sec')
        self.dt=dt

    def print_statistics(self, best_tile=None, json_filename=None, ts=None):
        # Unique targets
        #~ science_targets=best_tile.TaipanTile.get_assigned_targets_science()
        #~ ids=[x.idn for x in science_targets]
        #~ self.unique_targets.update(ids)
        #~ unique_tmp=set(ids).difference(self.unique_targets) # unique targets in the best_tile tile
        
        #~ self.total_number_cumulative_unique += len(unique_tmp)
        #~ self.total_number_cumulative += len(ids)
        
        self.number_of_tiles_observed += 1
        
        #~ # Repeated observations: just make list of all observations of all stars for now
        #~ for x in ids:
            #~ self.repeats[x]=ts
        
        #~ self.f.write('%s; %d; %d; %d; %d; %d; %d; %f; %.1f; %s \n'%(ts, best_tile.TaipanTile.field_id, len(ids), len(unique_tmp), self.total_number_cumulative_unique, self.total_number_cumulative, best_tile.TaipanTile.priority, best_tile.weight*1000.0, best_tile.TaipanTile.mag_max, json_filename))
        self.f.write('%s; %d; %d; %f; %.1f; %s \n'%(ts, best_tile.TaipanTile.field_id, best_tile.TaipanTile.priority, best_tile.weight*1000.0, best_tile.TaipanTile.mag_max, json_filename))
        
        #~ print len(self.unique_targets), len(unique_tmp), self.total_number_cumulative_unique, self.total_number_cumulative

    def calibration_frames(self, t=None, ts=None):
        dt=TimeDelta(0, format='sec')
        
        # Biases every evening
        # TODO: what if FW has no obs time and only Taipan is observing?
        is_night2 = simulate_weather.is_it_night(t.mjd)
        if self.is_night is False and is_night2 is True: # start of the night
            self.is_night=True
            self.fc.write('%s; %d; bias \n'%(ts, 15*2)) # 15 frames per arm
        if self.is_night is True and is_night2 is False: # change to day in the morning
            self.is_night=False     
        
        # Arcs and flats every hour       
        if self.bad_weather_duration<TimeDelta(3600.0, format='sec'): # no arcs and flats if there is more than 1 hour of bad weather in a row
            self.hourly_calibration_frames+=self.dt
            if self.hourly_calibration_frames>TimeDelta(3600.0, format='sec'):
                self.hourly_calibration_frames=TimeDelta(0, format='sec')
                self.fc.write('%s; %d; arc flat \n'%(ts, 4)) # 4 frames: 1 arc and 1 flat per arm
                dt=TimeDelta(2*60.0, format='sec') # 2 minutes for hourly arcs and flats
        
        return dt
    
    def run(self):
        t=params_simulator.params['date_start']
        self.hourly_calibration_frames=TimeDelta(0, format='sec')

        self.is_night = False
        self.bad_weather_duration=TimeDelta(0, format='sec')

        while t < params_simulator.params['date_finish'] and self.number_of_tiles_observed < self.number_of_tiles and self.time_with_no_observing < TimeDelta(3600.0*24.0*10.0, format='sec'):
            ts=str(t)[:-7]
            date=ts[:10]

            # NIGHT
#~ times = midnight + delta_midnight
#~ altazframe = AltAz(obstime=times, location=bear_mountain)
#~ sunaltazs = get_sun(times).transform_to(altazframe)
#~ sunaltazs.alt
            
            
            sun_set = self.observatory.sun_set_time(Time(date) + TimeDelta(3600.0, format='sec')).datetime     
            sun_rise = self.observatory.sun_rise_time(Time(date) - TimeDelta(3600.0, format='sec') + TimeDelta(1.0, format='jd')).datetime # NEXT DAY
            #~ print date, sun_set, sun_rise
            if t>=sun_set and t<sun_rise:
                pass # it is night
            else: # daytime
                t+=TimeDelta(10*60.0, format='sec') # TODO: check if this increment is consistent with other things
                continue
        

            # Weather
            # TODO: sky illumination is included, so only night time is considered (when Sun below horizon). What about the Moon?
            if simulate_weather.is_weather_good(t.mjd):
                msg = 'Weather good.'
                self.find_next_tile(ts=ts)
            
            elif simulate_weather.is_weather_acceptable_for_a_bright_tile(t.mjd):
                msg = 'Weather acceptable for a bright tile.'
                # TODO: determine limiting_magnitude
                self.find_next_tile(ts=ts)
            
            else:
                msg = 'Weather bad.'
                self.dt=TimeDelta(10*60.0, format='sec') # jump for 10 minutes and hope weather gets better
                self.bad_weather_duration+=TimeDelta(10*60.0, format='sec')

            # If hourly arcs and flats are taken, dt_cal is not 0.
            dt_cal=self.calibration_frames(t=t, ts=ts)

            t+=self.dt+dt_cal
            print t, 'dt', self.dt, 'finish:', params_simulator.params['date_finish'], msg
            
            #~ hourly_calibration_frames+=self.dt
            #~ if hourly_calibration_frames>TimeDelta(3600.0, format='sec'):
                #~ hourly_calibration_frames=TimeDelta(0, format='sec')
                #~ t+=TimeDelta(2*60.0, format='sec') # 2 minutes for hourly arcs and flats
                #~ self.fc.write('%s; %d; arc flat \n'%(ts, 4)) # 2 minutes

            print ''

        #~ with open('%srepeats.pkl'%params.params['data_output_folder'], 'wb') as f:
            #~ pickle.dump(self.repeats, f)
        
            if self.time_with_no_observing > TimeDelta(3600.0*24.0*10.0, format='sec'):
                break
        
        self.__finish__()
    
if __name__ == "__main__":
    sim = Simulator()
    sim.run()
