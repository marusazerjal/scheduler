import numpy as np
import math
import random
import pickle
import json
import urllib2
import datetime
from collections import defaultdict

from astropy.time import Time
from astropy.coordinates import get_moon
from astropy.coordinates import SkyCoord
from astropy.coordinates.earth import EarthLocation
import astropy.units as u
from astroplan import Observer
from astroplan import FixedTarget

import line_profiler
import pdb

import taipan.core as tp

# FunnelWeb Scheduler
import manage_list_of_observed_tiles
import params

reload(manage_list_of_observed_tiles)
reload(params)


ALT_MIN = params.params['ALT_MIN']
ALT_MAX = params.params['ALT_MAX']
HOUR_ANGLE_AMP = params.params['HOUR_ANGLE_AMP']
MOON_ANGDIST_MIN = params.params['MOON_ANGDIST_MIN']
MOON_ANGDIST_OK = params.params['MOON_ANGDIST_OK']
SLEW_TIME_MIN = params.params['SLEW_TIME_MIN'] # Do I need this?

class ObsTile():
    """
    Observational parameters of a tile. This class makes no observational decisions, it only determines parameters.
    """
    def __init__(self, tp=None, local_sidereal_time=None, moon=None, max_priority=None):
        """
        Parameters
        ----------        
        TaipanTile:
        local_sidereal_time: Local sidereal time, hours
        moon: astropy moon object
        
        max_priority: Maximal tile priority in the current tiling. Used for normalization.
        
        hour_angle: Hour angle of a tile at time local_sidereal_time, hours
        alt: altitude of a center of a tile for a given hour_angle, degrees
        alt_max: maximal altitude of a center of a tile (in upper culmination), degrees
        alt_diff: 
        angular_moon_distance: Angular distance between the tile and the Moon, degrees
        """
        self.TaipanTile=tp
        self.lat=params.params['LAT']
        self.max_priority=max_priority # To normalize ranking
        self.local_sidereal_time=local_sidereal_time
        self.moon=moon
        #~ self.observatory=observatory

        """
        Compute observational parameters
        """
        self.hour_angle=self.determine_hour_angle()
        self.alt=self.altitude_of_an_object_in_the_sky()
        
        # These params could be set later, only for visible tiles
        self.alt_max=self.max_altitude_of_an_object_in_the_sky()
        self.angular_moon_distance=self.distance_between_two_points_in_the_sky(alpha1=self.TaipanTile.ra, delta1=self.TaipanTile.dec, alpha2=self.moon.ra.value, delta2=self.moon.dec.value)


    @property
    def local_sidereal_time(self):
        return self.local_sidereal_time

    @local_sidereal_time.setter
    def local_sidereal_time(self, value):    
        self.local_sidereal_time = value
        
    def __str__(self):
        string = 'TP TILE %s: RA=%s, Dec=%s, Ranking=%s, w=%s, Altitude=%s, H=%s, Moon_dist=%s, mag_max=%s' % (('%d'%self.TaipanTile.field_id).rjust(5), ('%3.1f'%self.TaipanTile.ra).rjust(5), ('%2.1f'%self.TaipanTile.dec).rjust(5), ('%d'%self.TaipanTile.priority).rjust(5), ('%2.4f'%(self.weight*1000.0)).rjust(8), ('%d'%self.alt).rjust(2), ('%.2f'%self.hour_angle).rjust(6), ('%d'%self.angular_moon_distance).rjust(3), ('%d'%self.TaipanTile.mag_max).rjust(2)) #+ ' surface density %f %d %d'%(self.surface_density[0], self.surface_density[1], self.surface_density[2])
        return string

    #~ def __repr__(self):
        #~ string = 'TP TILE %d: RA=%3.1f, Dec=%2.1f, Ranking=%d, Altitude=%d, %d, %d, H=%.2f, Moon_dist=%d, mag_max=%d' % (self.TaipanTile.field_id, self.TaipanTile.ra, self.TaipanTile.dec, self.TaipanTile.priority, self.alt, self.alt_max, self.alt_diff, self.hour_angle, self.angular_moon_distance, self.TaipanTile.mag_max)
        #~ return string

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

    def determine_azimuth(self, ra=None, dec=None, alt=None, H=None):
        """
        Azimuth of the tile.
        """
        #~ alpha=np.deg2rad(ra)
        # TODO: remove ra.
        delta=np.deg2rad(dec)
        h=np.deg2rad(alt)
        phi=np.deg2rad(self.lat)
        hour_angle = np.deg2rad(H*15.0)
        
        #~ cosA = (np.sin(delta) - np.sin(h)*np.sin(phi)) / (np.cos(h)*np.cos(phi))
        #~ sinA = - np.sin(hour_angle) * np.cos(delta) / np.cos(h)

        cosA = (np.sin(h)*np.sin(phi) - np.sin(delta)) / (np.cos(h)*np.cos(phi))
        sinA = np.sin(hour_angle) * np.cos(delta) / np.cos(h)
        
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
            H = self.local_sidereal_time - self.TaipanTile.ra/15.0
        if H>12:
            H=24.0-H
        if H<-12:
            H=24.0+H
        # TODO: check again why does this happen
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
            dec=np.deg2rad(self.TaipanTile.dec)
        
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
        dec=np.deg2rad(self.TaipanTile.dec)
        lat=np.deg2rad(self.lat)
        sin_h = np.cos(dec-lat)
        h=np.arcsin(sin_h)
        h=np.rad2deg(h)
        return h

        

    def estimate_best_time_interval_to_observe_tile(self):
        """
        Determine hour angle amplitude to estimate time when observing is still acceptable
        """
        dec=np.deg2rad(self.TaipanTile.dec)
        lat=np.deg2rad(self.lat)
        
        # If star is above the horizon (our 'horizon' is at ALT_MIN) at all times:
        if self.lat+self.TaipanTile.dec+ALT_MIN<-90.0:
            distance_from_the_pole=(90.0+self.TaipanTile.dec)# self.dec<0
            alt_low=self.alt_max-2.0*distance_from_the_pole
            h_good_low = alt_low + 2.0*distance_from_the_pole*params.params['ALTITUDE_LOW_FRACTION']
        else:
            h_good_low = self.alt_max * params.params['ALTITUDE_LOW_FRACTION']

        if h_good_low<ALT_MIN:
            h_good_low=ALT_MIN            
        h_good_low = np.deg2rad(h_good_low)
        
        cosH = (np.sin(h_good_low) - np.sin(lat)*np.sin(dec)) / (np.cos(lat)*np.cos(dec))

        try:
            H = np.arccos(cosH)
            H = np.rad2deg(H)
            H = np.abs(H)
            H = H / 15.0 # hours
        except:
            H=None # TODO: how to treat |cosH|>1?? This is for circumpolar stars.
        return H

    def weighting(self, nearest_neighbours=None, observed_tile_id_internal=None, tiles_mag_range=None, ra_current=None, dec_current=None):
        """
        Weighting between H (hour angle), Ranking (priority) and slew time.
        """
        w_altitude = self.weighting_altitude() # [0, 1]
        w_slew_time = self.weighting_slew_time(ra_current=ra_current, dec_current=dec_current) # [0, 1]
        w_moon = self.weighting_moon_distance() # [0, 1]
        w_density = self.weighting_field_density_knn(nearest_neighbours=nearest_neighbours, observed_tile_id_internal=observed_tile_id_internal) # [0, 1]
        w_mag_range = self.weighting_magnitude_range(observed_tile_id_internal=observed_tile_id_internal, tiles_mag_range=tiles_mag_range)
        
        w_ranking = float(self.TaipanTile.priority) / float(self.max_priority) # [0, 1]
        # Some tiles have ranking equal to 0. So we set probability to 0.5:
        #~ if w_ranking<1e-12:
            #~ w_ranking=0.5
        
        w = w_ranking * w_altitude * w_slew_time * w_moon * w_density * w_mag_range

        self.weight = w

        
    def weighting_altitude(self):
        """
        Return altitude_current / altitude_at_meridian.
        """
        alt=self.alt
        alt_max=self.alt_max # altitude at local meridian

        # If star is above the horizon (our 'horizon' is at ALT_MIN) at all times:
        if self.lat+self.TaipanTile.dec+ALT_MIN<-90.0:
            distance_from_the_pole=(90.0+self.TaipanTile.dec)# self.dec<0
            alt_min=alt_max-2.0*distance_from_the_pole
            w = float(alt-alt_min)/float(2.0*distance_from_the_pole)
        else:
            w = float(alt)/float(alt_max)
        
        return w

    def weighting_moon_distance(self):
        d=self.angular_moon_distance
        if d<MOON_ANGDIST_MIN:
            w=0.0
        elif d>MOON_ANGDIST_OK:
            w=1
        else:
            k=(1.0-0.5)/(MOON_ANGDIST_OK-MOON_ANGDIST_MIN)
            n=0.5-k*MOON_ANGDIST_MIN
            w=k*d+n
        return w 

    def weighting_slew_time(self, ra_current=None, dec_current=None):
        """
        Slew time function: a maximum of angular distance and difference in azimuth.
        """       
        # Current position:
        ra1=ra_current
        dec1=dec_current
        
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
        #~ result=np.exp(-result)

        return result

    def weighting_magnitude_range(self, observed_tile_id_internal=None, tiles_mag_range=None):
        """
        Fainter magnitude ranges contain more tiles and are thus selected more likely. This makes all magnitude ranges to be selected with equal probability.
        """
        mag_range=(self.TaipanTile.mag_min, self.TaipanTile.mag_max)
        tls=tiles_mag_range[mag_range]

        n_observed=float(len(tls.intersection(observed_tile_id_internal)))        
        n_total=float(len(tls))

        w = 1.0 - n_observed / n_total

        return w

    def weighting_field_density(self, tiles_mag_range=None, internal_observed_tiles_mag_range=None):
        """
        Determine local field density (for a particular magnitude range)
        """
        ra=self.TaipanTile.ra # deg
        dec=self.TaipanTile.dec # deg
        
        mag_range=(float(self.TaipanTile.mag_min), float(self.TaipanTile.mag_max))

        try:
            observed_tiles=internal_observed_tiles_mag_range[mag_range]
        except:
            observed_tiles=[]
            
        n=0
        for x in observed_tiles: # TODO: make hour_angle exclusion loop (because after time number of observed tiles is going to grow (or perhaps not because new tiling will be generated each month)
            if np.abs(ra-x.ra)<params.params['TILE_DENSITY_RADIUS'] and np.abs(dec-x.dec)<params.params['TILE_DENSITY_RADIUS']:
                d=self.distance_between_two_points_in_the_sky(alpha1=ra, delta1=dec, alpha2=x.ra, delta2=x.dec)
                if d<params.params['TILE_DENSITY_RADIUS']:
                    n+=1
        n_observed=float(n)


        # All candidate tiles from the tiling code within tile_density_radius
        n=0
        for x in tiles_mag_range[mag_range]:
            if np.abs(ra-x.ra)<params.params['TILE_DENSITY_RADIUS'] and np.abs(dec-x.dec)<params.params['TILE_DENSITY_RADIUS']:
                d=self.distance_between_two_points_in_the_sky(alpha1=ra, delta1=dec, alpha2=x.ra, delta2=x.dec)
                if d<params.params['TILE_DENSITY_RADIUS']:
                    #~ print ra, x.ra, dec, x.dec
                    n+=1
        n_total=float(n)

        
        if n<0.5:
            #~ test=np.array(sorted([[x.ra, x.dec] for x in tiles_mag_range[mag_range]], key=lambda x: x[0]))
            #~ print test
            for x in tiles_mag_range[mag_range]:
                print x.ra, x.dec
            pdb.set_trace()
            
            return 0.0 # apparently out of the range
        
        weight = 1.0 - n_observed / n_total

        self.surface_density=[weight, n_observed, n_total]


        return weight

    def weighting_field_density_knn(self, nearest_neighbours=None, observed_tile_id_internal=None):
        """
        Determine local field density (for a particular magnitude range).
        List of nearest neighbours contains only neighbours within magnitude range of a tile in question
        """
        
        nn=nearest_neighbours[self.TaipanTile.field_id]
        nn=nn[:params.params['K_NEAREST_NEIGHBOURS']]
        nn=set([int(x[0]) for x in nn]) # get rid of distances, keep only field_id
        
        number_of_observed_tiles=len(nn.intersection(observed_tile_id_internal))

        n_observed=float(number_of_observed_tiles)
        n_total=float(len(nn))
        
        weight = 1.0 - n_observed / n_total

        return weight
    
