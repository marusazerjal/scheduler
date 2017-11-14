"""
FunnelWeb Scheduler: Decide which tile to observe next, based on sidereal time and observing conditions (weather, seeing etc.).
"""

import numpy as np
import math
import pickle
import json
#~ import urllib2
import datetime
from collections import defaultdict

from astropy.time import Time
from astropy.time import TimeDelta
from astropy.coordinates import get_moon
from astroplan import Observer

import taipan.core as tp

# FunnelWeb Scheduler
import manage_list_of_observed_tiles
import params
import obstile
import visualization

"""
Load the tiling data
"""
reload(manage_list_of_observed_tiles)
reload(params)
reload(obstile)
reload(visualization)

print 'Input tiling file', params.params['input_tiling_filename']
try:
    TILES = data[0]
    SETTINGS = data[2]
except NameError:
    print 'Reading data...'
    fl = open(params.params['input_tiling_filename'], 'rb')
    data = pickle.load(fl)
    fl.close()				

    TILES = data[0]
    SETTINGS = data[2]

print 'Number of input tiles:', len(TILES)



print 'Nearest neighbours file', params.params['nearest_neighbours_filename']
try:
    NN = nndata
except NameError:
    print 'Reading nearest neighbours...'
    fl = open(params.params['nearest_neighbours_filename'], 'rb')
    nndata = pickle.load(fl)
    fl.close()				

    NN = nndata


print 'List of observed tiles', params.params['observed_tiles_internal_filename']


class Scheduler():
    """
    Schedule FunnelWeb tile observations.
    """
    def __init__(self):
        """
        Parameters
        ----------
        ra_current, dec_current: float, degrees
            Current position of telescope. This position is needed to determine slew time to the next target. This data comes directly from Jeeves.
            For the start of the night it is best to start at meridian (telescope is not parked there). Later on take current telescope position into account.
            
        Data needed from Jeeves:
            ra_current, dec_current
            weather
        """
        self.tiles = TILES # TODO: do I now have 2 copies of TILES?
        self.determine_internal_tile_id_and_priorities()
        self.observed_tiles = manage_list_of_observed_tiles.load_internal_list_of_observed_tile_ids()
        self.group_tiles_per_magnitude_range() # Need this for magnitude range weight determination

        self.observatory = Observer.at_site("Anglo-Australian Observatory") # TODO: enter LAT and LON coordinates

    def determine_internal_tile_id_and_priorities(self):
        """
        Determine internal field_id and priority.
        """
        tiles2=[]
        priorities=[]
        for tile_id, x in enumerate(self.tiles):
            x.field_id=tile_id
            prior=x.calculate_tile_score(method=SETTINGS['ranking_method'], disqualify_below_min=SETTINGS['disqualify_below_min'], combined_weight=SETTINGS['combined_weight'], exp_base=SETTINGS['exp_base'])
            x.priority=prior
            priorities.append(prior)
            tiles2.append(x)
        self.tiles=tiles2
        self.max_priority=np.max(priorities)

    # TODO
    def next_tile(self, ra_current=None, dec_current=None):
        """
        Find the next tile to be observed. This is what Jeeves is going to call each time.
        Add ra_current and dec_current (and weather). This is how Jeeves is going to call this method.
        """
        # TODO: Now or in 10 minutes?
        self.utc = Time(Time.now(), format='iso', scale='utc')
        self.local_sidereal_time = self.utc.sidereal_time('mean', longitude=params.params['LON']).value

        self.moon = get_moon(self.utc)

        #~ self.weather_conditions()
        #~ if weather_data_filename:
            #~ self.weather_data=np.loadtxt(weather_data_filename) # STRING!!
        #~ else:
            #~ self.weather_data=None

        # Telescope position from Jeeves
        self.ra_current=ra_current
        self.dec_current=dec_current
        if ra_current is None:
            self.ra_current=self.local_sidereal_time
        if dec_current is None:
            self.dec_current=-30.0 # TODO

        best_tile=self.find_best_tile()

        # Update list of observed tiles
        manage_list_of_observed_tiles.add_tile_id_internal_to_the_list({best_tile.TaipanTile.field_id})
        
        # TODO: what should be a format for Jeeves?
        return best_tile

    def observing_plan(self, date=None, time=None):
        """
        Make a list of best tiles observable through the night.
        
        Output: observing_plan.txt
        # observing_ideal observing_start observing_stop fieldID path_to_config_file
        1030 1000 1100 fieldID /observers_files/taipan/YYYYMMDD/tile_pk_HHMMSS.obs_config.json
        
        observing_ideal = meridian transit time
        observing_start, observing_stop: observing_ideal +/- dt, dt is determined using weights
        
        Should times be given in UT?
        What is the json file (e.g. list of stars with coordinates etc.?)
        
        TODO: What is ObsConfig file?
        """     
        if date is None:
            datenow=datetime.datetime.now().date()
            date='%d-%02d-%02d'%(datenow.year, datenow.month, datenow.day)
            date = Time('%s'%date)

        sun_set = self.observatory.sun_set_time(Time(date)).datetime
        sunset=self.observatory.datetime_to_astropy_time(sun_set) # converts into different format (datetime to astropy)
        if time is not None:
            todo=True
            # If time > sunset: continue from 'time' on.
            # In case observing plan is re-run during the night.
        
        sun_rise = self.observatory.sun_rise_time(Time(date) + TimeDelta(1.0, format='jd')).datetime
        
        f=open(params.params['observing_plan_filename'], 'wb')

        # UTC times
        times=[]
        for i in range(100):
            dt = datetime.timedelta(minutes=i*params.params['TIME_PER_TILE'])
            t = sun_set + dt
            if t < sun_rise:
                times.append(t)
            else:
                break
        
        print 'Selecting tiles for', date
        timezone_correction=TimeDelta(3600.0*11.0, format='sec')
        sunset_lt=self.observatory.datetime_to_astropy_time(sun_set) + timezone_correction
        sunrise_lt=self.observatory.datetime_to_astropy_time(sun_rise) + timezone_correction
        dark_time = sun_rise-sun_set
        dark_time = dark_time.seconds
        dark_time_hours = int(dark_time/3600.0)
        dark_time_minutes = (dark_time - float(dark_time_hours)*3600.0)/60.0
        print 'Sunset LT %d-%02d-%02d %02d:%02d'%(sunset_lt.value.year, sunset_lt.value.month, sunset_lt.value.day, sunset_lt.value.hour, sunset_lt.value.minute)
        print 'Sunrise LT %d-%02d-%02d %02d:%02d'%(sunrise_lt.value.year, sunrise_lt.value.month, sunrise_lt.value.day, sunrise_lt.value.hour, sunrise_lt.value.minute)
        print 'Nighttime duration %02d:%02d'%(dark_time_hours, dark_time_minutes)
        print 'Number of tiles this night:', len(times)
        print
        
        # Convert UTC times to LST
        times_lst = [Time(t).utc.sidereal_time('mean', longitude=params.params['LON']).value for t in times]

        time_efficiency=[]
        telescope_positions=[] # for testing purposes

        selected_tiles=[]
        count=1
        new_tile_ids_to_be_added_to_list_of_all_observed_stars=set()
        for t, lst in zip(times, times_lst):
            t_start=datetime.datetime.now()      
            utc = self.observatory.datetime_to_astropy_time(t)
            self.moon = get_moon(utc)

            try:
                self.ra_current=best_tile.TaipanTile.ra
                self.dec_current=best_tile.TaipanTile.dec
            except:
                self.ra_current=lst
                self.dec_current=-70.0 # TODO

            self.local_sidereal_time=lst # TODO: check if this violates any other things. Why cant I insert LST to init_best_tiles...??
  
            best_tile = self.find_best_tile()

            # Update list of observed tiles
            self.observed_tiles.add(best_tile.TaipanTile.field_id)
            new_tile_ids_to_be_added_to_list_of_all_observed_stars.add(best_tile.TaipanTile.field_id)
            
            """
            Code from this point on is only output formatting
            """

            local_time=t+datetime.timedelta(hours=11.0)
            
            print '%02d/%02d'%(count, len(times)), 'LT=%d-%02d-%02d %02d:%02d:%02d'%(local_time.year, local_time.month, local_time.day, local_time.hour, local_time.minute, local_time.second), 'UT=%d-%02d-%02d %02d:%02d:%02d'%(times[count-1].year, times[count-1].month, times[count-1].day, times[count-1].hour, times[count-1].minute, times[count-1].second), 'LST=%s'%(('%02.2f'%lst).rjust(5)), best_tile #, best_tile.meridian_transit_time.datetime

            telescope_positions.append([best_tile.TaipanTile.ra, best_tile.TaipanTile.dec, self.moon.ra.value, self.moon.dec.value, best_tile.angular_moon_distance])
            

            # Print out the data
            H_amp = best_tile.estimate_best_time_interval_to_observe_tile()
            if H_amp is None:
                todo=True
                
            # TODO: check if this times are during the night, not e.g. just before sunset or just after sunrise, otherwise limit them within sunset and sunrise time
            observing_start = t - datetime.timedelta(hours=H_amp)
            observing_stop = t + datetime.timedelta(hours=H_amp)
            observing_ideal = Time(t - datetime.timedelta(hours=best_tile.hour_angle)).utc.sidereal_time('mean', longitude=params.params['LON']).value # Do I insert time NOW or time of meridian crossing?
            
            
            """
            Print output
            """
            line='%02d%02d %02d%02d %02d%02d %05d /observers_files/funnelweb/YYYYMMDD/%05d_HHMMSS.obs_config.json'%(t.hour, t.minute, observing_start.hour, observing_start.minute, observing_stop.hour, observing_stop.minute, best_tile.TaipanTile.field_id, best_tile.TaipanTile.field_id)
           
            # TODO: what happens with observing_ideal for tiles at ALT=90? Because there is a limit at 85 degrees.
            f.write(line+'\n')
        
            count+=1
            
            t_end=datetime.datetime.now() 
            time_efficiency.append((t_end-t_start).seconds)
            
            # PLOT
            #~ visualization.plot_selected_tile_with_neighbourhood(moon=self.moon, lst=lst, best_tiles=self.best_tiles_to_observe_now, tiles=self.tiles, best_tile=best_tile, i=count-1, ra_current=self.ra_current, dec_current=self.dec_current, telescope_positions=telescope_positions, observed_tile_ids=self.observed_tiles)            
            
        f.close()
        
        manage_list_of_observed_tiles.add_tile_id_internal_to_the_list(new_tile_ids_to_be_added_to_list_of_all_observed_stars)
        
        print 'Average time to find the next tile: [seconds]', np.mean(time_efficiency)
        
        #~ telescope_positions=np.array(telescope_positions)
        #~ return telescope_positions

    def find_best_tile(self):
        """
        This is a master function. Use multiple methods to find a best tile to be observed now.
        First, all tiles are assumed to be appropriate to observe. Then we exclude tiles that are not appropriate at the moment.
        """
        # Exclude: global
        self.select_only_tiles_within_hour_angle_amp(S=self.local_sidereal_time)
        #~ self.observed_tiles=manage_list_of_observed_tiles.load_internal_list_of_observed_tile_ids()
        self.exclude_tiles_already_observed()

        # Initialize
        self.best_tiles_to_observe_now = [obstile.ObsTile(tp=x, local_sidereal_time=self.local_sidereal_time, moon=self.moon, max_priority=self.max_priority) for x in self.tiles_tmp]
        
        # Exclude: local
        # Position restrictions
        self.lower_altitude_limit() # Tiles visible on this day
        self.exclude_zenith()
        self.exclude_tiles_too_close_to_the_moon()
        
        # Magnitude restrictions
        #~ self.exclude_tiles_below_limiting_magnitude() # add this later when weather is taken into account
        
        # Determine weights for candidate tiles:
        for x in self.best_tiles_to_observe_now:
            x.weighting(nearest_neighbours=NN, observed_tile_id_internal=self.observed_tiles, tiles_mag_range=self.tiles_mag_range, ra_current=self.ra_current, dec_current=self.dec_current)
        
        b = sorted(self.best_tiles_to_observe_now, key=lambda y: y.weight, reverse=True)
        self.best_tiles_to_observe_now=b
        
        # Get data
        best_tile=self.best_tiles_to_observe_now[0]
        return best_tile        

    def group_tiles_per_magnitude_range(self):
        """
        Distribute tiles into magnitude ranges.
        """
        tiles_mag_range = defaultdict(list)
        for x in self.tiles:
            tiles_mag_range[tuple([float(x.mag_min), float(x.mag_max)])].append(x.field_id)
        self.tiles_mag_range={k: set(v) for k, v in tiles_mag_range.iteritems()}        

    def lower_altitude_limit(self):
        """
        Find tiles with altitude high enough to be observable.
        """
        visible_tiles=[]
        for x in self.best_tiles_to_observe_now:
            if x.alt>params.params['ALT_MIN']:
                visible_tiles.append(x)
        self.best_tiles_to_observe_now=visible_tiles

    def exclude_zenith(self):
        """
        Exclude zenith due to technical obstacles.
        """
        visible_tiles=[]
        for x in self.best_tiles_to_observe_now:
            if x.alt<params.params['ALT_MAX']:
                visible_tiles.append(x)
        self.best_tiles_to_observe_now=visible_tiles

    def select_only_tiles_within_hour_angle_amp(self, S=None):
        """
        Exclude tiles that are not visible right now. Keep only those close to the local meridian (H_amp away from the local meridian at most).
        """
        H_amp = params.params['HOUR_ANGLE_AMP']
        
        ra_min = S - H_amp
        if ra_min<0.0:
            ra_min=24.0+ra_min
        ra_min *= 15.0
        
        ra_max = S + H_amp
        if ra_max>24.0:
            ra_max=ra_max-24.0
        ra_max *= 15.0

        if ra_min<ra_max:
            b=[x for x in self.tiles if x.ra>ra_min and x.ra<ra_max]
        elif ra_min>ra_max:
            b1=[x for x in self.tiles if x.ra>ra_min]
            b2=[x for x in self.tiles if x.ra<ra_max]
            b=b1+b2

        self.tiles_tmp=b
            
    def exclude_tiles_too_close_to_the_moon(self):
        """
        Exclude tiles closer than self.moon_angdist_min from the Moon.
        """
        good_tiles=[x for x in self.best_tiles_to_observe_now if x.angular_moon_distance>params.params['MOON_ANGDIST_MIN']]
        self.best_tiles_to_observe_now=good_tiles    
    
    def exclude_tiles_below_limiting_magnitude(self):
        """
        Exclude tiles within magnitude ranges that include stars fainter than the limiting magnitude.
        """
        if self.limiting_magnitude:
            good_tiles=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.mag_max<=self.limiting_magnitude]
            self.best_tiles_to_observe_now=good_tiles
        else:
            pass

    def exclude_tiles_already_observed(self):
        """
        Exclude tiles that have already been observed
        """
        # TODO: could be faster if I used sets instead of lists. But is is not possible to use set for self.tiles_tmp in this moment.
        tiles2=[x for x in self.tiles_tmp if x.field_id not in self.observed_tiles]

        self.tiles_tmp=tiles2
                
    def weather_conditions(self):
        """
        Define limiting magnitude regarding the weather conditions.

        Parameters
        ----------
        seeing:
        
        Return:
        ----------        
        limiting_magnitude: float
            Limiting magnitude for given weather conditions.

        WHAT is a relation between seeing and limiting magnitude?
        
        sky transparency (--> limiting magnitude) 
        
        
        Sky transparency: -10 is probably rain, -20 is clear (with possible cirruses). But this is the average.       
        
        
        """

        #~ if self.weather_data:
            # seeing=?
                        
            #~ self.limiting_magnitude=TODO

        # TODO: read this file only once in the __init__, not for each tile.
        data = urllib2.urlopen("http://site.aao.gov.au/met/metdata.dat").read(500) # read only 20 000 chars
        data = data.split("\n") # then split it into lines
        
        
        # Check if date and time are close enough to our observations
        
        date=data[0].replace('"', '').replace('.', '').replace(' ', '')
        date_spl=date.split('-')
        month=int(date_spl[0])
        day=int(date_spl[1])
        year=int(date_spl[2])
        
        spl=data[1].split('\t')
        
        time=spl[0].split(':')
        hour=int(time[0])
        minute=int(time[1])
        
        time_weather = datetime.datetime(year, month, day, hour, minute)
        time_now = datetime.datetime.now()
        
        if time_weather>time_now:
            dt=time_weather-time_now
        else:
            dt=time_now-time_weather
        
        dt=str(dt).split(':')
        dt=float(dt[0])*60.0+float(dt[1])
        print 'weather dt', dt
        # Note: There is 1 hour difference between Stromlo and SSO. Why??!?
        
        #~ seeing=float(spl[13]) # there is no seeing info in the data
        sky_transparency=float(spl[13])

        #~ print 'seeing, sky_transparency', seeing, sky_transparency
        print 'sky_transparency', sky_transparency

    def print_selected_tile_to_json(self):
        """
        Generate dictionary of output data and write to a json file.
        """
        # Generate dictionary
        t = self.best_tiles_to_observe_now[0]
        data={
        'tile_id': t.TaipanTile.field_id
        }
        #~ 'TaipanTile': t.TaipanTile # Class is not serializable
        
        # Write
        # Each time add timestamp to json file. Is this a good idea? Because then Jeeves has to search for the latest filename each time.
        timestamp=self.utc.value.replace(' ', '--')
        filename='selected_tile_%s.json'%timestamp
        with open(filename, 'w') as outfile:  
            json.dump(data, outfile)    
        print '%s created.'%filename


"""
-------------------
Tests
-------------------
All the possible situations:
- only part of the sky is clear

"""

def next_tile():
    """
    Run scheduler for testing purposes.
    """
    s=Scheduler()

    best_tile = s.next_tile()
    print best_tile  

def observing_plan(date=None, time=None):
    """
    Run observing plan for testing purposes.
    """

    t_start=datetime.datetime.now()

    s=Scheduler()
    s.observing_plan(date=date, time=time)
    
    t_end=datetime.datetime.now()
    dt=t_end-t_start
    print 'Total time: %d s (%.2f min)'%(dt.seconds, dt.seconds/60.0)
            
if __name__ == "__main__":

    """
    Run Scheduler
    """
    next_tile()
    
    #~ observing_plan(date='2017-11-04', time='02:42:42')
    #~ observing_plan(date='2017-11-03') # default: tonight
    #~ observing_plan() # default: tonight
