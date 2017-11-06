"""
Params file for FunnelWeb Scheduler
"""

# Essential params

# Input tiling file (including path to this tile):
input_tiling_filename = '171308_1647_fw_tiling.pkl'

#~ dir = '/Users/mireland/Google Drive/FunnelWeb/TargetSelection/tiling_results/'
#~ tiling_filename=dir + '171308_1647_fw_tiling.pkl'

# No need to edit
params={'input_tiling_filename': input_tiling_filename,
        'observed_tiles_internal_filename': '%s_observed_tiles_internal.dat'%(input_tiling_filename[:-11]),
        'observed_tiles_external_filename': '%s_observed_tiles_external.dat'%(input_tiling_filename[:-11])
        }


# Supplementary params (if not given, default values will be assigned)
args={'date': '2017-11-04', # Default: today
        'time': '02:42:42' # Default: Now
        }

# todo: assign time here
