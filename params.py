"""
Params file for FunnelWeb Scheduler
"""

import os

# Input data
input_tiling_filename_folder = 'data_input/'
input_tiling_filename = '180222_2311_21_fw_tiling.pkl'
#~ input_tiling_filename = '180409_0157_25_fw_tiling.pkl' # just two different priorities

data_output_folder = 'data_output/data_paper2/'

# Create this folder
if not os.path.exists(data_output_folder):
    os.makedirs(data_output_folder)

"""
CONSTANTS
"""
params={'input_tiling_filename': '%s%s'%(input_tiling_filename_folder, input_tiling_filename),
        'nearest_neighbours_filename': '%snearest_neighbours_%s.pkl'%(input_tiling_filename_folder, input_tiling_filename[:-11]),
        'observed_tiles_internal_filename': '%s%s_observed_tiles_internal.dat'%(data_output_folder, input_tiling_filename[:-11]),
        'observed_tiles_external_filename': '%s%s_observed_tiles_external.dat'%(data_output_folder, input_tiling_filename[:-11]),
        'observing_plan_filename': '%sobserving_plan'%(data_output_folder),
        #~ 'obs_config_json_folder': 'observers_files/funnelweb/',
        'obs_config_json_folder': '%sjson/'%data_output_folder,
        
        'exponent_base_add': 4.0, # add this number to the default exponent base of 3
        'highest_tile_score_if_any_priority_5_star': True,
        'highest_tile_score_if_any_priority_5_star_value': 1e+8,
        
        'ALT_MIN': 27.0, # Minimal altitude of the tile good to observe
        'ALT_MAX': 85.0, # Maximal altitude of the tile good to observe
        'ALTITUDE_LOW_FRACTION': 0.8, # Fraction = Current altitude / altitude at local meridian. Altitude_low_fraction is thus the lowest acceptable fraction for a tile to be observed.

        'HOUR_ANGLE_AMP': 2,#3, # hours, consider only tiles with H = +/- HOUR_ANGLE_AMP
        # But be careful: when moon is on the local meridian, what happens?

        # TODO: Should I include magnitude dependent Moon distances?
        'MOON_ANGDIST_MIN': 20.0, # Minimal acceptable distance
        'MOON_ANGDIST_OK': 30.0, # Distance where sky illumination by the Moon is negligible

        'TIME_PER_TILE': 10, # [minutes], observing plan, 10 minutes altogether per tile
        #~ 'SLEW_TIME_MIN': 60.0, # [seconds], take slew time into consideration above this limit

        #~ 'TILE_DENSITY_RADIUS': 6.0, #deg, OBSOLETE; radius to compute observed tile surface density
        'K_NEAREST_NEIGHBOURS': 15, # Number of nearest neighbours to consider instead of TILE_DENSITY_RADIUS
        
        #~ 'N_BEST_TILES_MERIDIAN': 10, # Select tile closest to the meridian among the first N_BEST_TILES_MERIDIAN tiles.
        'CONSIDER_TILES_ABOVE_THIS_WEIGHT': 0.95,

        # Siding Spring Observatory
        # TODO: this data is automatic result from Google. Check if it is correct.
        'LAT': -31.2749, # deg South
        'LON': 149.0685, # deg East
        'EL': 1164.0, # m, elevation
        
        
        # SIMULATOR
        'data_output_folder': data_output_folder, # simulator needs this
        
        # if ranking is 0
        #~ 'weighting_ranking_0': False,
        #~ 'w_ranking_0_value': 1e-6
        }
