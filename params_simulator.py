"""
Params file for FunnelWeb Scheduler
"""
from datetime import date

# Essential params

# Input tiling file (including path to this tile):
input_tiling_filename_folder = 'data_input/'
#~ input_tiling_filename = '171308_1647_fw_tiling.pkl'
input_tiling_filename = '170812_1742_31_fw_tiling.pkl' # with FW input catalogue
#~ input_tiling_filename = '170509_1943_fw_tiling.pkl' # 600mb
#~ data_output_folder = 'data_output/'
data_output_folder = 'data_output_simulator/'



params={'simulate_dates_file': 'simulator/simulator_dates3.dat',
           'MIN_MOON_PHASE_TO_OBSERVE': 0.5, # FunnelWeb gets bright time. This is around 16 nights per month with phase>0.5.
           'date_start': date(2018, 4, 14),  # start date
           'date_finish': date(2018, 4, 24),  # end date
           'weather_clouds_stats': 'simulator/coonabarabran_cloudy_days_stats.dat',
           'technical_issues_rate': 0.05 # How many times observations are not possible due to the technical issues
        }
