"""
Params file for FunnelWeb Scheduler Simulator
"""
from astropy.time import Time
import params
reload(params)

data_output_folder = params.params['data_output_folder']

params={   'simulator_statistics_output': '%ssimulator_statistics.dat'%data_output_folder,
           'simulator_statistics_output_calibration': '%ssimulator_statistics_calibration.dat'%data_output_folder,
           
           'date_start': Time('2018-06-22 11:15'), # UT
           'date_finish': Time('2021-12-22 12:00'),
           #~ 'date_start': Time('2018-07-03 11:15'), # UT
           #~ 'date_finish': Time('2018-07-04 12:00'),
           
           'reconfig_time': 5.0*60.0, # seconds
           'exposure_time': {8.5: 60.0, 10.5: 60.0, 12.5: 300.0, 14.5: 600.0}, # mag is mag_max of the tile, exposure time is in seconds
           'slew_time': 60.0, # we are not sure if we can configure during slewing
           'readout_time': 21.0, # long discussion, check notes if this is correct
           
           'weather_sso_database': 'simulator/fulldb_cadence_5mins.csv',
           'volume_per_exposure': 68 # 2 x 34 MB raw data files per exposure
        }
