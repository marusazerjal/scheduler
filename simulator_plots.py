'''
Scheduler simulation over a few years.
Use one tiling file (it includes most of the stars from the priority list) and simulate observing nights through years. Include weather, etc.
'''

import numpy as np

from datetime import date, timedelta
from os import listdir
from os.path import isfile, join
import simplejson as json
import matplotlib.pylab as plt

from astropy.time import Time


import params_simulator

# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.
# It is better to change data_output_folder = 'data_output/' to something else (in the params.py file). Don't forget to add nearest_heighbours file.

#~ '''
#~ Params
#~ '''
#~ MIN_MOON_PHASE_TO_OBSERVE=0.5 # FunnelWeb gets bright time. This is around 16 nights per month with phase>0.5.
#~ D1 = date(2018, 3, 14)  # start date
#~ D1 = date(2018, 4, 14)  # start date
#~ D2 = date(2018, 5, 14)  # end date
#~ D2 = date(2019, 5, 25)  # end date
#~ DELTA = D2 - D1         # timedelta
#~ SIMULATE_DATES_FILE='simulator/simulator_dates2.dat'

    
def statistics():
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    #~ print dates
    
    # Read json files for each date
    nights={}
    for date in dates:
        path=params.params['obs_config_json_folder']+date
        json_files = [f for f in listdir(path)] # if isfile(join(mypath, f))
        #~ print date, json_files
        #~ print
        
        # Merge all json data together. We want data for the entire night.
        mags_tmp=[]
        coo_tmp=[]
        coo_tile_tmp=[]
        for j in json_files:
            path2=path+'/'+j
            #~ print path2
            with open(path2) as json_data:
                d = json.load(json_data)
                #~ print [k for k in d.iterkeys()]
                targets=d['targets']
                #~ print len(targets)
                #~ print targets
                for tg in targets:
                    mags_tmp.append(tg['mag'])
                    coo_tmp.append([tg['ra'], tg['dec']])
                
                tl=d['fieldCentre']
                coo_tile_tmp=[tl['ra'], tl['dec']] # ra and dec are strings!
            #~ break
        
        nights[date]={'mag': mags_tmp, 'coo': coo_tmp, 'coo_tiles': coo_tile_tmp}
        # coo_tiles: because we want percentage of the field done
        #~ break
    
    # TODO: ORDER nights for time simulation
    #~ print nights
    keys=[k for k in nights.iterkeys()]
    keys=sorted(keys)
    
    # magnitude distributions
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for k, v in nights.iteritems():
        print k
        mag=v['mag']
        ax.hist(mag, histtype='step')
    
    
    # Number of stars over time. Cumulative
    number_of_stars_per_night=np.array([[int(k), len(nights[k]['coo'])] for k in keys])
    number_of_tiles_per_night=np.array([[int(k), len(nights[k]['coo_tiles'])] for k in keys])
    
    print number_of_stars_per_night
    n=0
    number_of_stars_per_night2=copy.deepcopy(number_of_stars_per_night)
    #~ day0=number_of_stars_per_night2[0,0] # wrong. it is not continuous.
    for i in range(len(number_of_stars_per_night)):
        print number_of_stars_per_night2[i,1], n
        #~ number_of_stars_per_night2[i,0]-=day0
        number_of_stars_per_night2[i,0]=i
        n1=number_of_stars_per_night2[i,1]
        number_of_stars_per_night2[i,1]+=n
        n+=n1
    number_of_stars_per_night2=np.array(number_of_stars_per_night2)
    print number_of_stars_per_night2
    
    fig=plt.figure()
    ax=fig.add_subplot(211)
    ax.plot(number_of_stars_per_night2[:,0], number_of_stars_per_night2[:,1])    
    #~ ax=fig.add_subplot(212)
    #~ ax.hist(number_of_tiles_per_night, histtype='step', cumulative=True)    
    
    plt.show()

    
if __name__ == "__main__":
    statistics()
