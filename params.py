"""
Params file for FunnelWeb Scheduler
"""

# Essential params

# Input tiling file (including path to this tile):
input_tiling_filename_folder = 'data_input/'
input_tiling_filename = '171308_1647_fw_tiling.pkl'
data_output_folder = 'data_output/'

#~ dir = '/Users/mireland/Google Drive/FunnelWeb/TargetSelection/tiling_results/'
#~ tiling_filename=dir + '171308_1647_fw_tiling.pkl'


#~ # Supplementary params (if not given, default values will be assigned)
args={'date': '2017-11-10', # Default: today
        'time': '02:42:42' # Default: Now
        }

# todo: assign time here


"""
CONSTANTS
"""
params={'input_tiling_filename': '%s%s'%(input_tiling_filename_folder, input_tiling_filename),
        'observed_tiles_internal_filename': '%s%s_observed_tiles_internal.dat'%(data_output_folder, input_tiling_filename[:-11]),
        'observed_tiles_external_filename': '%s%s_observed_tiles_external.dat'%(data_output_folder, input_tiling_filename[:-11]),
        'observing_plan_filename': '%sobserving_plan_%s.dat'%(data_output_folder, args['date']),
        'nearest_neighbours_filename': '%snearest_neighbours_%s.pkl'%(data_output_folder, input_tiling_filename[:-11]),
        
        'ALT_MIN': 27.0, # Minimal altitude of the tile good to observe
        'ALT_MAX': 85.0, # Maximal altitude of the tile good to observe
        'ALTITUDE_LOW_FRACTION': 0.8, # Fraction = Current altitude / altitude at local meridian. Altitude_low_fraction is thus the lowest acceptable fraction for a tile to be observed.

        'HOUR_ANGLE_AMP': 2,#3, # hours, consider only tiles with H = +/- HOUR_ANGLE_AMP
        # But be careful: when moon is on the local meridian, what happens?

        # TODO: Should I include magnitude dependent Moon distances?
        'MOON_ANGDIST_MIN': 20.0, # Minimal acceptable distance
        'MOON_ANGDIST_OK': 30.0, # Distance where sky illumination by the Moon is negligible

        'TIME_PER_TILE': 10,#10 # [minutes], 10 minutes altogether per tile
        'SLEW_TIME_MIN': 60.0, # [seconds], take slew time into consideration above this limit

        'TILE_DENSITY_RADIUS': 6.0, #deg, OBSOLETE; radius to compute observed tile surface density
        'K_NEAREST_NEIGHBOURS': 15, # Number of nearest neighbours to consider instead of TILE_DENSITY_RADIUS

        # Siding Spring Observatory
        # TODO: this data is automatic result from Google. Check if it is correct.
        'LAT': -31.2749, # deg South
        'LON': 149.0685 # deg East
        }


