"""
Params file for FunnelWeb Scheduler
"""
from astropy.time import Time
import params
reload(params)

data_output_folder = params.params['data_output_folder']

#~ simulation_nickname = 'priorities11'

params={   
           #~ 'simulation_nickname': simulation_nickname,
           'simulator_statistics_output': '%ssimulator_statistics.dat'%data_output_folder,
           'simulator_statistics_output_calibration': '%ssimulator_statistics_calibration.dat'%data_output_folder,
           #~ 'observing_plan_folder': 'data_output_simulator/observing_plans_%s/'%simulation_nickname,
           #~ 'obs_config_json_folder': 'observers_files/funnelweb_simulation_%s/'%simulation_nickname,
           
           'date_start': Time('2018-06-22 11:15'), # UT
           'date_finish': Time('2018-06-23 08:00'),
           
           'reconfig_time': 5.0*60.0, # seconds
           'exposure_time': {8.5: 60.0, 10.5: 60.0, 12.5: 300.0, 14.5: 600.0}, # mag is mag_max of the tile, exposure time is in seconds
           
           'weather_sso_database': 'simulator/fulldb_cadence_5mins.csv',
           'volume_per_exposure': 68 # 2 x 34 MB raw data files per exposure
        }
