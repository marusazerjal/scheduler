"""
Params file for FunnelWeb Scheduler
"""
from astropy.time import Time

simulation_nickname = 'priorities11'

params={   'simulation_nickname': simulation_nickname,
           'simulator_statistics_output': 'simulator/simulator_statistics_%s.dat'%simulation_nickname,
           'simulator_statistics_output_calibration': 'simulator/simulator_statistics_calibration_%s.dat'%simulation_nickname,
           'observing_plan_folder': 'data_output_simulator/observing_plans_%s/'%simulation_nickname,
           'obs_config_json_folder': 'observers_files/funnelweb_simulation_%s/'%simulation_nickname,
           'date_start': Time('2018-06-22 11:15'), # UT
           'date_finish': Time('2021-06-25 08:00'),
           'weather_sso_database': 'simulator/fulldb_cadence_5mins.csv',
           'volume_per_exposure': 68 # 2 x 34 MB raw data files per exposure
        }
