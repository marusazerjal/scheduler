"""
FunnelWeb: Decide which tile to observe next, based on sidereal time and observing conditions (weather, seeing etc.).
"""

# Should this code select only one tile to be observe next, or a few best candidates?

# Assumption: none of the tiles from the input file have not been observed yet. If Scheduler is run muptiple times during the night make sure you exclude tiles that have already been observed this night.
# Tiling will be run about once per month anyway so there should be a list of observed tiles.


import numpy as np
import math
import random
import pickle
import json
import datetime

from astropy.time import Time
from astropy.coordinates import get_moon
from astropy.coordinates import SkyCoord
from astropy.coordinates.earth import EarthLocation
import astropy.units as u
from astroplan import Observer
from astroplan import FixedTarget

import line_profiler

# Taipan
import taipan.core as tp
# Have to run the code in python 2. Astropy requires numpy>1.9.1. How to do that? It doesn't work with pip

"""
Constants
"""
# Minimal and maximal allowed altitude of the center of the tile
ALT_MIN=27.0
ALT_MAX=85.0

# Minimal angular Moon distance
MOON_ANGDIST_MIN=20.0

TIME_PER_TILE=10 # [minutes], 10 minutes altogether per tile
SLEW_TIME_MIN=60.0 # [seconds], take slew time into consideration above this limit


"""
Observer's coordinates
"""
# Siding Spring Observatory
# TODO: this data is automatic result from Google. Check if it is correct.
LAT=-31.2749 # deg South
LON=149.0685 # deg East


tiling_filename='171308_1647_fw_tiling.pkl'
try:
    tiles = data[0]#[:10]
    settings = data[2]
except NameError:
    print 'Reading data...'
    fl = open(tiling_filename,'rb')
    data = pickle.load(fl)
    fl.close()				

    tiles = data[0]#[:10]
    settings = data[2]

print 'Number of tiles:', len(tiles)



class Scheduler():
    """
    Schedule FunnelWeb tile observations.
    """
    def __init__(self, tiling_filename=None, observed_tiles_filename=None, weather_data_filename=None, alt_min=None, alt_max=None, ra_min=None, ra_max=None, moon_angdist_min=None, lat=None, lon=None, desired_priority=None, minimal_priority=None, time_per_tile=TIME_PER_TILE, limiting_magnitude=None, slew_time_min=SLEW_TIME_MIN, ra_current=None, dec_current=None):
        """
        Parameters
        ----------
        tiling_filename: string
            A pickle data filename of tiles generated by the tiling code.
        observed_tiles_filename: string
            A file with a list of tile IDs (a list of integers) that have already been observed. TODO: MAKE SURE THAT TILE IDs ARE UNIQUE.
        weather_data_filename: string
            Filename of the weather data.
        alt_min: float, degrees
            Minimal allowed altitude of a tile above the horizon. Default: 30.0.
        alt_max: float, degrees
            Maximal allowed altitude of a tile. (Avoid observing zenith). Default: 70.0.
        ra_min, ra_max: float, degrees
            Limit tiles within certain RA limits (e.g. when only part of the sky is clear).
        moon_angdist_min: float, degrees
            Angular distance between the Moon and the tile. Default: 20.0.
        lat, lon: float, degrees
            Observer's latitude and longitude. Default: lat=-31.2749, lon=149.0685 (Siding Spring).
        desired_priority: float
            Desired priority to be observed.
        minimal_priority: float
            User defined minimal tile priority to be observed.
        processing_time: float, seconds
            Total time spend per tile (slewing, exposure, calibration etc.). Default: 300.0.
        limiting_magnitude: float
            Magnitude of a faintest star to be observed. Default: None.
        ra_current, dec_current: float, degrees
            Current position of telescope. This position is needed to determine slew time to the next target.
            
            
        Data needed from Jeeves:
            ra_current, dec_current
            weather
        
        """
        
        """
        Read the data
        """
        #~ if tiling_filename:
            #~ try:
                #~ self.tiles = data[0]#[:10]
                #~ self.settings = data[2]
            #~ except NameError:
                #~ print 'Reading data...'
                #~ fl = open(tiling_filename,'rb')
                #~ data = pickle.load(fl)
                #~ fl.close()				
            
                #~ self.tiles = data[0]#[:10]
                #~ self.settings = data[2]
            
            #~ print 'Number of tiles:', len(self.tiles)
        #~ else:
            #~ print 'tiling_filename is None.'
            #~ exit(0)

        self.tiles = tiles
        self.settings = settings
        
        if observed_tiles_filename:
            try:
                observed_tile_ids=np.loadtxt(observed_tiles_filename)
                self.observed_tile_ids=observed_tile_ids # a list of integers (IDs)
            except:
                self.observed_tile_ids=[]
                print 'WARNING: File with tile IDs already observed cannot be read.'
        else:
            self.observed_tile_ids=[]
        
        
        if weather_data_filename:
            self.weather_data=np.loadtxt(weather_data_filename) # STRING!!
        else:
            self.weather_data=None
        
        """
        Determine field_id and priority.
        Calculate_the_score should remain here.
        """
        # TODO: Later remove field_id as tile_id will be already available in the TaipanTiles.
        tiles2=[]
        priorities=[]
        #~ print 'WARNING: CORRECT TILE PRIORITIES!'
        for tile_id, x in enumerate(self.tiles):
            x.field_id=tile_id
            prior=x.calculate_tile_score(method=self.settings['ranking_method'], disqualify_below_min=self.settings['disqualify_below_min'], combined_weight=self.settings['combined_weight'], exp_base=self.settings['exp_base'])
            x.priority=prior
            #~ x.priority=tile_id
            priorities.append(prior)
            tiles2.append(x)
        self.tiles=tiles2
        self.max_priority=np.max(priorities)
        #~ print priorities
        
        #~ prior=sorted(priorities)
        #~ print prior
        
        #~ exit(0)

        """
        Observer's coordinates
        """
        if lat:
            self.lat=lat
        else:
            self.lat=LAT
        if lon:
            self.lon=lon
        else:
            self.lon=LON

        #observatory=EarthLocation('Anglo-Australian Observatory')
        #~ observatory = Observer(at_site="Anglo-Australian Observatory")
        #~ sun_set = apo.sun_set_time(time, which="next")



        """
        Minimal altitude of a tile above the horizon.
        """
        if alt_min:
            if alt_min>ALT_MIN:
                self.alt_min=alt_min
            else: # user is not allowed to do funny things
                self.alt_min=ALT_MIN
                print 'WARNING: Minimal altitude too low. Set to %g degrees.'%ALT_MIN
        else:
            self.alt_min = ALT_MIN # degrees


        """
        Maximal tile altitude valid to observe - due to technical obstacles (cannot observe zenith). ???? DO WE NEED THIS?
        """
        if alt_max:
            if alt_max<ALT_MAX:
                self.alt_max=alt_max
            else: # user is not allowed to do funny things
                self.alt_max=ALT_MAX
                print 'WARNING: Maximal altitude limit too high. Set to %g degrees.'%ALT_MAX
        else:
            self.alt_max = ALT_MAX # degrees


        """
        Limit RA in case of partly cloudy sky.
        """
        if ra_min:
            if ra_min>1e-10 and ra_min<360.0+1e-10:
                self.ra_min=ra_min
            else:
                self.ra_min=None
                print 'WARNING: You must enter 0 < ra_min < 360.'
        else:
            self.ra_min=None
            
        if ra_max:
            if ra_max>1e-10 and ra_max<360.0+1e-10:
                self.ra_max=ra_max
            else:
                self.ra_max=None
                print 'WARNING: You must enter 0 < ra_max < 360.'
        else:
            self.ra_max=None

        if self.ra_min and self.ra_max:
            if self.ra_min>self.ra_max or np.abs(self.ra_min-self.ra_max)<1e-6:
                print 'WARNING: ra_min and ra_max are not within valid limits or ra_min>ra_max.'



        """
        Limiting magnitude
        """
        if limiting_magnitude:
            self.limiting_magnitude=limiting_magnitude
        else:
            self.limiting_magnitude=None

        """
        Minimal angular distance to the Moon
        """
        if moon_angdist_min is not None:
            if moon_angdist_min>MOON_ANGDIST_MIN: # user is not allowed to do funny things
                self.moon_angdist_min=moon_angdist_min
        else:
            self.moon_angdist_min = MOON_ANGDIST_MIN # degrees    


        """
        ------------------------------------
        Time dependent variables
        ------------------------------------
        """

        """
        Determine local sidereal time
        """
        self.utc = Time(Time.now(), format='iso', scale='utc')
        self.local_sidereal_time = self.utc.sidereal_time('mean', longitude=self.lon).value


        """
        Current telescope position
        """
        if ra_current and dec_current:
            self.ra_current=ra_current
            self.dec_current=dec_current
        else:
            # todo: at what coordinates is telescope parked?
            self.ra_current=self.local_sidereal_time*15.0
            self.dec_current=-30.0


        """
        Moon
        """
        self.moon = get_moon(self.utc)
        # Astropy Moon coordinates agree with Stellarium within a degree. This is good enough for the Scheduler.
        #~ exit(0)


        # User defined settings
        # Manually enter desired priority of the field. If None, tiles with highest priorities will be selected automatically.
        self.desired_priority=desired_priority
        self.minimal_priority=minimal_priority
        # TODO: Error check if value entered is valid
        
        # Processing time (time per tile)
        if time_per_tile:
            self.time_per_tile=time_per_tile
        else:
            self.time_per_tile=TIME_PER_TILE


        # At the beginning all tiles are assumed good to observe, then we throw away those not appropriate at the moment
        #~ self.best_tiles_to_observe_now=[ObsTile(TaipanTile=x) for x in self.tiles]
        self.best_tiles_to_observe_now=[ObsTile(tp=x, lat=self.lat, local_sidereal_time=self.local_sidereal_time, moon=self.moon, ra_current=self.ra_current, dec_current=self.dec_current, max_priority=self.max_priority) for x in self.tiles]


    def __str__(self):
        string="""
        __str__:
        TODO=True
        """
        return string

    def __repr__(self):
        string="""
        __repr__:
        TODO=True
        """
        return string

    #~ @profile
    def find_the_best_tile_to_observe_now(self): # A combination of different methods
        """
        This is a master function. Use multiple methods to find a best tile to be observed now.
        First, all tiles are assumed to be appropriate to observe. Then we exclude tiles that are not appropriate at the moment.
        """
        
        # Exclude tiles already observed. Only tiles from the input tiling file are considered.
        self.exclude_tiles_already_observed()
        
        # Position restrictions
        self.lower_RA_limit()
        self.upper_RA_limit()
        self.lower_altitude_limit()
        self.exclude_zenith()
        self.exclude_tiles_too_close_to_the_moon()
        
        # Magnitude restrictions
        self.exclude_tiles_below_limiting_magnitude() # user defined
        
        # Priorities
        # Determine weights for candidate tiles:
        for x in self.best_tiles_to_observe_now:
            x.weighting()
        
        b = sorted(self.best_tiles_to_observe_now, key=lambda y: y.weight, reverse=True)
        self.best_tiles_to_observe_now=b
        
        #~ self.find_tiles_closest_to_local_meridian() # TODO: tile_priority overrides this sorting
        #~ self.tile_priority()
        
        # Get data
        best_tile=self.best_tiles_to_observe_now
        return best_tile

    def make_list_of_best_tiles_through_the_night(self): # FIX THIS!!!
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
        todo=True
      
      
        # Observatory
        observatory = Observer.at_site("Anglo-Australian Observatory") # TODO: enter LAT and LON coordinates
        
        # Determine sunset and sunrise times
        datenow=datetime.datetime.now().date()
        date='%d-%02d-%02d'%(datenow.year, datenow.month, datenow.day)
        time = Time('%s 23:30:00'%date)
        
        sun_set = observatory.sun_set_time(time, which="previous").datetime
        sun_rise = observatory.sun_rise_time(time, which="previous").datetime # previous. This is correct.
        
        sunset=observatory.datetime_to_astropy_time(sun_set)

        # Start of the night, end of the night. Time step: self.processing_time
        #~ self.utc = Time(Time.now(), format='iso', scale='utc')
        #~ self.utc.sidereal_time('mean', longitude=self.lon).value
        
        #~ dt=datetime.timedelta(minutes=self.time_per_tile) # in minutes

        # UTC times
        times = [sun_set + datetime.timedelta(minutes=i*self.time_per_tile) for i in range(1, 100) if sun_set + datetime.timedelta(minutes=i*self.time_per_tile)<=sun_rise]
        
        # Convert UTC times to LST
        times_lst = [Time(t).utc.sidereal_time('mean', longitude=self.lon).value for t in times]

        selected_tiles=[]
        for t, lst in zip(times, times_lst):            
            #~ utc = Time('%d-%02d-%02d %02d:%02d:%02d'%(t.year, t.month, t.day, t.hour, t.minute, t.second), format='iso', scale='utc')
            utc = observatory.datetime_to_astropy_time(t)
            moon = get_moon(utc)

            self.best_tiles_to_observe_now=[ObsTile(tp=x, lat=self.lat, local_sidereal_time=lst, moon=moon, observatory=observatory, sunset=sunset, ra_current=self.ra_current, dec_current=self.dec_current, max_priority=self.max_priority) for x in self.tiles]

            
            # Find best tile for time t
            tls=self.find_the_best_tile_to_observe_now()
            
            # Tiles cannot be selected twice
            do_selection=True
            i=0
            while do_selection:
                best_tile=tls[i]
                if best_tile.TaipanTile.field_id not in selected_tiles:
                    selected_tiles.append(best_tile.TaipanTile.field_id)
                    break
                else:
                    i+=1
            
            print t, 'LST=%s'%(('%.11f'%lst).rjust(14)), best_tile#, best_tile.meridian_transit_time.datetime
            
            # Provisional data. UPDATE!!
            observing_start=sun_set # This could be perhaps when a tile is at 70% of its maximum altitude. Should be larger fraction for more northern tiles.
            observing_stop=sun_rise
            observing_ideal = 'meridian transition time'
            #~ print '%02d%02d %02d%02d %02d%02d %05d /observers_files/funnelweb/YYYYMMDD/%05d_HHMMSS.obs_config.json'%(t.hour, t.minute, observing_start.hour, observing_start.minute, observing_stop.hour, observing_stop.minute, best_tile.TaipanTile.field_id, best_tile.TaipanTile.field_id)
        
        # Order by hour angle, starting from the evening.
        

    def exclude_tiles_already_observed(self):
        """
        Exclude tiles with tile IDs in the list of tiles already observed.
        """
        if len(self.observed_tile_ids)>0:
            good=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.field_id not in self.observed_tile_ids]
            self.best_tiles_to_observe_now=good
        else:
            pass
    

    def lower_altitude_limit(self):
        """
        Find tiles with altitude high enough to be observable.
        """
        visible_tiles=[]
        for x in self.best_tiles_to_observe_now:
            if x.alt>self.alt_min:
                visible_tiles.append(x)
        self.best_tiles_to_observe_now=visible_tiles

    def exclude_zenith(self):
        """
        Exclude zenith due to technical obstacles.
        """
        visible_tiles=[]
        for x in self.best_tiles_to_observe_now:
            if x.alt<self.alt_max:
                visible_tiles.append(x)
        self.best_tiles_to_observe_now=visible_tiles

    def lower_RA_limit(self):
        """
        If self.ra_min is not None, take only tiles with centre RA>=ra_min
        """
        if self.ra_min:
            good=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.ra>=self.ra_min]
            self.best_tiles_to_observe_now=good
            if len(good)<1:
                print 'WARNING: There is no tiles with RA>=%f.'%self.ra_min
        else:
            pass

    def upper_RA_limit(self):
        """
        If self.ra_max is not None, take only tiles with centre RA<=ra_min
        """
        if self.ra_max:
            good=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.ra<=self.ra_max]
            self.best_tiles_to_observe_now=good
            if len(good)<1:
                print 'WARNING: There is no tiles with RA<=%f.'%self.ra_max
        else:
            pass

        
    def find_tiles_closest_to_local_meridian(self):
        """
        Sort tiles according to their hour_angle. Tiles closest to local meridian are listed first.
        """
        tiles=sorted(self.best_tiles_to_observe_now, key=lambda x: np.abs(x.hour_angle))
        self.best_tiles_to_observe_now=tiles
            
    def exclude_tiles_too_close_to_the_moon(self):
        """
        Exclude tiles closer than self.moon_angdist_min from the Moon.
        """
        good_tiles=[x for x in self.best_tiles_to_observe_now if x.angular_moon_distance>self.moon_angdist_min]
        self.best_tiles_to_observe_now=good_tiles
        #~ print 'Number of tiles far enough from the Moon', len(good_tiles)
    
    
    def exclude_tiles_below_limiting_magnitude(self):
        """
        Exclude tiles within magnitude ranges that include stars fainter than the limiting magnitude.
        """
        if self.limiting_magnitude:
            good_tiles=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.mag_max<=self.limiting_magnitude]
            self.best_tiles_to_observe_now=good_tiles
        else:
            pass
        
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
        
        
        """
        todo=True
        
        
        #~ if self.weather_data:
            # seeing=?
            
            
            #~ self.limiting_magnitude=TODO
        
        
            
    def tile_priority(self):
        # Tile priority: summing up the priorities of the TaipanTargets within the tile.
        if self.minimal_priority:
            good=[x for x in self.best_tiles_to_observe_now if x.TaipanTile.priority>=self.minimal_priority]
            if len(good)<1:
                print 'WARNING: There is no tiles with desired minimal_priority.'
                self.best_tiles_to_observe_now=good
            else:
                good=sorted(good, key=lambda x: x.TaipanTile.priority, reverse=True)
        else:
            # Choose the highest priority tiles
            good=sorted(self.best_tiles_to_observe_now, key=lambda x: x.TaipanTile.priority, reverse=True)
        self.best_tiles_to_observe_now=good    
        # If there is no tiles with minimal_priority, select the next one


        # TODO: DIFFICULTY
        # x.difficulty Difficulty: summing up the difficulties of the TaipanTargets within the tile. ##### WHAT IS THIS? it is probably relevant


    def weight_field_density(self):
        """
        Weighting for field density. I.e. the number of fields above a cutoff priority in each region of the sky (or each declination).
        """
        todo=True
    

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

class ObsTile():
    """
    Observational parameters of a tile. This class makes no observational decisions, it only determines parameters.
    """
    #~ @profile
    def __init__(self, tp=None, lat=None, local_sidereal_time=None, moon=None, observatory=None, sunset=None, ra_current=None, dec_current=None, max_priority=None):
        """
        Parameters
        ----------        
        TaipanTile:
        ra, dec: tile coordinates, degrees
        lat: observer's latitude, degrees
        local_sidereal_time: Local sidereal time, hours
        moon: astropy moon object
        observatory: astroplan observatory object
        sunset: astropy Time object: sunset time
        
        hour_angle: Hour angle of a tile at time local_sidereal_time, hours
        alt: altitude of a center of a tile for a given hour_angle, degrees
        alt_max: maximal altitude of a center of a tile (in upper culmination), degrees
        alt_diff: 
        angular_moon_distance: Angular distance between the tile and the Moon, degrees
        """
        self.TaipanTile=tp
        self.ra=tp.ra
        self.dec=tp.dec
        self.lat=lat
        self.ra_current=ra_current
        self.dec_current=dec_current
        self.max_priority=max_priority
        
        self.local_sidereal_time=local_sidereal_time
        self.moon=moon
        self.observatory=observatory
        self.sunset=sunset

        """
        Compute observational parameters
        """
        self.hour_angle=self.determine_hour_angle()
        self.alt=self.altitude_of_an_object_in_the_sky()
        self.alt_max=self.max_altitude_of_an_object_in_the_sky()
        self.alt_diff=self.alt_max-self.alt
        self.angular_moon_distance=self.distance_between_two_points_in_the_sky(alpha1=self.ra, delta1=self.dec, alpha2=self.moon.ra.value, delta2=self.moon.dec.value)
        
        #~ self.meridian_transit_time=None
        #~ try:
            #~ self.next_meridian_transit_time()
        #~ except:
            #~ pass
        #~ self.weight=self.weighting()

    @property
    def local_sidereal_time(self):
        return self.local_sidereal_time

    @local_sidereal_time.setter
    def local_sidereal_time(self, value):    
        self.local_sidereal_time = value
        
    def __str__(self):
        string = 'TP TILE %s: RA=%s, Dec=%s, Ranking=%s, w=%s, Altitude=%s, H=%s, Moon_dist=%s, mag_max=%s' % (('%d'%self.TaipanTile.field_id).rjust(5), ('%3.1f'%self.TaipanTile.ra).rjust(5), ('%2.1f'%self.TaipanTile.dec).rjust(5), ('%d'%self.TaipanTile.priority).rjust(5), ('%2.4f'%self.weight).rjust(8), ('%d'%self.alt).rjust(2), ('%.2f'%self.hour_angle).rjust(6), ('%d'%self.angular_moon_distance).rjust(3), ('%d'%self.TaipanTile.mag_max).rjust(2))
        return string

    def __repr__(self):
        string = 'TP TILE %d: RA=%3.1f, Dec=%2.1f, Ranking=%d, Altitude=%d, %d, %d, H=%.2f, Moon_dist=%d, mag_max=%d' % (self.TaipanTile.field_id, self.TaipanTile.ra, self.TaipanTile.dec, self.TaipanTile.priority, self.alt, self.alt_max, self.alt_diff, self.hour_angle, self.angular_moon_distance, self.TaipanTile.mag_max)
        return string

    #~ @profile
    def distance_between_two_points_in_the_sky(self, alpha1=None, delta1=None, alpha2=None, delta2=None):
        """
        Determine the distance between two points in the sky.
        
        Parameters
        ----------
        alpha1, delta1: float, degrees
            RA and Dec for object 1
        alpha2, delta2: float, degrees
            RA and Dec for object 2
        
        Output
        ----------
        A: float, degrees
            Distance between two points in the sky.
        """
        a1=np.deg2rad(alpha1)
        d1=np.deg2rad(delta1)
        a2=np.deg2rad(alpha2)
        d2=np.deg2rad(delta2)
        
        cos_A = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(a1-a2)
        A=np.arccos(cos_A) # [0, 180] in radians
        A=np.rad2deg(A)
        return A # degrees

    #~ @profile
    def determine_azimuth(self, ra=None, dec=None, alt=None, H=None):
        #~ alpha=np.deg2rad(ra)
        delta=np.deg2rad(dec)
        h=np.deg2rad(alt)
        phi=np.deg2rad(self.lat)
        hour_angle = np.deg2rad(H*15.0)
        
        cosA = (np.sin(delta) - np.sin(h)*np.sin(phi)) / (np.cos(h)*np.cos(phi))
        sinA = - np.sin(hour_angle) * np.cos(delta) / np.cos(h)
        
        A = math.atan2(sinA, cosA)
        A = np.rad2deg(A)
        
        #~ A=np.arccos(cosA)
        #~ A=np.rad2deg(A)
        
        #~ H=np.deg2rad(H*15.0)
        #~ if np.sin(H)>0:
            #~ A=360.0-A
    
            
        return A

    def determine_hour_angle(self, ra=None):
        """
        Determine hour angle of a tile.
        
        Output
        ----------
        H: float, hours
            Hour angle of the tile.        
        """
        if ra:
            H = self.local_sidereal_time - ra/15.0
        else:
            H = self.local_sidereal_time - self.ra/15.0
        #~ print 'hour angle', self.local_sidereal_time
        if H>12:
            H=24.0-H
        return H

    def altitude_of_an_object_in_the_sky(self, hour_angle=None, dec=None):
        """
        Determine altitude of a tile above horizon.
        
        Output
        ----------
        h: float, degrees
            Altitute of the tile above horizon.        
        """
        # No need for atmospheric refraction correction.
        # Disagrees with Stellarium for about 5 degrees (Stellarium is using apparent altitude, but 5 degrees difference was for a star in zenith!!)
        
        if hour_angle and dec:
            hour_angle=np.deg2rad(hour_angle*15.0)
            dec=np.deg2rad(dec)
        else:
            hour_angle=np.deg2rad(self.hour_angle*15.0)
            dec=np.deg2rad(self.dec)
        
        lat=np.deg2rad(self.lat)
        sin_h = np.cos(hour_angle)*np.cos(dec)*np.cos(lat) + np.sin(dec)*np.sin(lat)
        h=np.arcsin(sin_h)
        h=np.rad2deg(h)
        return h

    def max_altitude_of_an_object_in_the_sky(self):
        """
        Determine maximal altitude of a tile above horizon.
        
        Output
        ----------
        h: float, degrees
            Maximal altitute of the tile above horizon.        
        """
        dec=np.deg2rad(self.dec)
        lat=np.deg2rad(self.lat)
        sin_h = np.cos(dec-lat)
        h=np.arcsin(sin_h)
        h=np.rad2deg(h)
        return h


    #~ def next_meridian_transit_time(self):
        #~ """
        #~ Next local meridian transit time.
        #~ """

        #~ coord = SkyCoord(ra=self.TaipanTile.ra*u.deg, dec=self.TaipanTile.dec*u.deg)
        #~ target = FixedTarget(coord=coord, name='TP TILE %d'%self.TaipanTile.field_id)
        #~ meridian_transit_time = self.observatory.target_meridian_transit_time(self.sunset, target, which="next")            
        #~ self.meridian_transit_time=meridian_transit_time


            
    #~ @profile
    def weighting(self):
        """
        Weighting between H (hour angle), Ranking (priority) and slew time.
        """

        w_altitude = self.weighting_altitude() # [0, 1]
        w_slew_time = self.weighting_slew_time() # [0, 1]
        w_ranking = float(self.TaipanTile.priority) / float(self.max_priority) # [0, 1]
        
        w = w_ranking * w_altitude * w_slew_time * 100.0
        #~ ww=[1.0, 1.0, 3.0]
        #~ w = (w_ranking*ww[1] + w_altitude*ww[1] + w_slew_time*ww[2]) / np.sum(ww) * 100.0
        
        #~ print w_altitude, w_slew_time, w_ranking, w
        
        self.weight = w
        
        #~ return weight
        
    def weighting_altitude(self):
        H=self.hour_angle
        dec=self.TaipanTile.dec
        alt=self.alt
        alt_max=self.alt_max # altitude at local meridian
        
        w = float(alt)/float(alt_max)
        
        return w

    #~ @profile
    def weighting_slew_time(self):
        """
        Slew time function: a maximum of angular distance and difference in azimuth.
        """       
        # Current position:
        ra1=self.ra_current
        dec1=self.dec_current
        
        # Next position:
        ra2=self.TaipanTile.ra
        dec2=self.TaipanTile.dec
        
        angle=self.distance_between_two_points_in_the_sky(alpha1=ra1, delta1=dec1, alpha2=ra2, delta2=dec2)
        if angle>180.0:
			angle=360.0-angle
        
        # Azimuth      
        H1=self.determine_hour_angle(ra=ra1)
        alt1=self.altitude_of_an_object_in_the_sky(hour_angle=H1, dec=dec1)
        az1=self.determine_azimuth(ra=ra1, dec=dec1, alt=alt1, H=H1)
        
        az2=self.determine_azimuth(ra=ra2, dec=dec2, alt=self.alt, H=self.hour_angle)
        
        az=np.abs(az2-az1)
        if az>180.0:
			az=360.0-az
        
          
        # If slew time is less than e.g. 60s (approximate readout time and minimum starbug reconfiguration) then ignore slew time.
        # TODO: how do I know relation between the angle and slew time?

        m=np.max([angle, az])
        result=(180.0 - m) / 180.0
        #~ result=((180.0 - m) / 180.0)**4
        #~ result=np.cos(m)
        result=np.exp(-result)
        #~ result=np.exp(-m)
        
        #~ print 'azimuth', ra1, ra2, dec1, dec2, az1, az2, az, m, result

        return result




                

# Why is altitude determined here, and not in Scheduler (where we can automatically skip tiles too low in the sky so we don't have to determine moon distance etc. of this file: if we are going to determine tile list for the entire night, then altitude changes overnight anyway. But...!

    

"""
-------------------
Tests
-------------------
All the possible situations:
- only part of the sky is clear

"""

def run_scheduler(tiling_filename, observed_tiles_filename):
    """
    Run scheduler for testing purposes.
    """
    # Basic example with no limitations other than altitude and Moon distance
    s=Scheduler(tiling_filename=tiling_filename, observed_tiles_filename=observed_tiles_filename)

    # Example with limitations
    #~ s=Scheduler(tiling_filename=tiling_filename, observed_tiles_filename=observed_tiles_filename, alt_min=20.0, alt_max=50.0, moon_angdist_min=20.0, ra_min=200, ra_max=250, limiting_magnitude=9)
    

    t=s.find_the_best_tile_to_observe_now()
    t0=t[0]
    #~ print 'Best tile:', t.alt, t.hour_angle, t.angular_moon_distance, t.TaipanTile.priority
    print 'Best tile:', t0
    s.print_selected_tile_to_json()
    
    for x in t[:10]:
        print x    


def observing_plan(tiling_filename, observed_tiles_filename):
    """
    Run observing plan for testing purposes.
    """

    s=Scheduler(tiling_filename=tiling_filename, observed_tiles_filename=observed_tiles_filename)
    s.make_list_of_best_tiles_through_the_night()


def test_tile_distribution_in_the_sky(filename):
    """
    Plot RA, Dec of all the candidate tiles from the tiling file.
    """
    fl = open(filename,'rb')
    tiles = pickle.load(fl)[0]
    fl.close()
    print 'Number of tiles:', len(tiles)
    
    c=np.array([[x.ra/15.0, x.dec] for x in tiles])
    
    import matplotlib.pylab as plt
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.scatter(c[:,0], c[:,1], s=5)
    plt.show()

def test_read_json_file(filename):
    """
    Check if json file is readable.
    """
    with open(filename) as json_file:  
        data = json.load(json_file)
    print data




def construct_weighting():
    lat=np.deg2rad(-31.0)
    
    dec=np.linspace(0, -90, 91)
    decr=[np.deg2rad(x) for x in dec]

    h_max = [np.rad2deg(np.arcsin(np.cos(x-lat))) for x in decr]
    
    import random
    htmp = [x*random.random() for x in h_max]
    
    p = [x/y for x, y, in zip(htmp, h_max)]
    # If I take only h_tmp / h_max, then the highest fields are observed very high, but in the same time also lower altitudes for such fields would be good enough.
    
    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dec, h_max)
    #~ ax.plot(dec, htmp)
    
    ax2=ax.twinx()
    ax2.plot(dec, p)
    ax2.axhline(y=0.9, color='k')
    
    w=[x*np.cos(np.deg2rad(d-lat)) for x, d in zip(p, h_max)]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dec, w)
    
    
    plt.show()


        
            
if __name__ == "__main__":
    #~ dir = '/Users/mireland/Google Drive/FunnelWeb/TargetSelection/tiling_results/'
    #~ tiling_filename=dir + '171308_1647_fw_tiling.pkl'
    
    tiling_filename='171308_1647_fw_tiling.pkl'
    #~ tiling_filename='171308_1647_fw_tiling_short.pkl'
    #~ tiling_filename='171308_1647_fw_tiling_short2.pkl'
    observed_tiles_filename='%s_observed_tiles.txt'%tiling_filename[:-4]

    """
    Run Scheduler
    """
    #~ run_scheduler(tiling_filename, observed_tiles_filename)
    
    observing_plan(tiling_filename, observed_tiles_filename)

    
    """
    Tests
    """
    #~ test_tile_distribution_in_the_sky(tiling_filename)
    #~ test_read_json_file('selected_tile_2017-10-09--01:42:30.379.json')
    #~ construct_weighting()
    


    #~ print 'Reading data...'
    #~ fl = open(tiling_filename,'rb')
    #~ data = pickle.load(fl)
    #~ fl.close()
    
    #~ n=100
    #~ d = data[0][:n]
    #~ s = data[2]
    #~ t = data[1][:n]
    #~ data=[d, t, s]
    
    #~ o=open('171308_1647_fw_tiling_short2.pkl', 'wb')
    #~ pickle.dump(data, o)
    #~ o.close()
    
    
    
    #~ import pdb
    #~ pdb.set_trace()
    
    #~ reload(fwtl)
    #~ fwtiler = FWTiler(...)
    #~ %run -i script.py
