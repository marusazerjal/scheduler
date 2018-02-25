"""
Params file for FunnelWeb Scheduler
"""
#~ from datetime import date
from astropy.time import Time

#~ data_output_folder = 'data_output/'
#~ data_output_folder = 'data_output_simulator/'

simulation_nickname = 'v1'

params={   'simulation_nickname': simulation_nickname,
           'simulator_statistics_output': 'simulator/simulator_statistics_%s.dat'%simulation_nickname,
           'observing_plan_folder': 'data_output_simulator/observing_plans_%s/'%simulation_nickname,
           'obs_config_json_folder': 'observers_files/funnelweb_simulation_%s/'%simulation_nickname,
           'date_start': Time('2018-06-22 11:15'), # UT
           'date_finish': Time('2030-10-25 08:00'),
           'weather_sso_database': 'simulator/fulldb_cadence_5mins.csv',
           'volume_per_exposure': 68 # 2 x 34 MB raw data files per exposure
        }




#~ params={'simulation_nickname': simulation_nickname,
           #~ 'simulate_dates_file': 'simulator/simulator_dates_bright_time.dat',
           #~ 'simulator_statistics_output': 'simulator/simulator_statistics_%s.dat'%simulation_nickname,
           #~ 'observing_plan_folder': 'data_output_simulator/observing_plans_%s/'%simulation_nickname,
           #~ 'MIN_MOON_PHASE_TO_OBSERVE': 0.5, # FunnelWeb gets bright time. This is around 16 nights per month with phase>0.5.
           #~ 'date_start': date(2014, 10, 24),  # start date
           #~ 'date_finish': date(2014, 10, 25),  # end date
           #~ 'weather_clouds_stats': 'simulator/coonabarabran_cloudy_days_stats.dat',
           #~ 'weather_sso_database': 'simulator/fulldb_cadence_5mins.csv',
           #~ 'technical_issues_rate': 0.05, # How many times observations are not possible due to the technical issues
           
           #~ 'obs_config_json_folder': 'observers_files/funnelweb_simulation_%s/'%simulation_nickname
        #~ }

