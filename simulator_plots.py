'''
Scheduler simulation over a few years.
Use one tiling file (it includes most of the stars from the priority list) and simulate observing nights through years. Include weather, etc.
'''

import numpy as np

#~ from datetime import date, timedelta
import datetime
from os import listdir
from os.path import isfile, join
import copy
import simplejson as json
import matplotlib.pylab as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook

from astropy.time import Time

import params
import params_simulator

# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.
    
def how_number_of_observed_stars_grows_with_time(path=None, simulator_dates_file=None):
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    if path is None: # Take the recent results
        path0=params.params['obs_config_json_folder']
        dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    else: # results from a motley simulation in a different folder
        path0=path
        dates=np.loadtxt(simulator_dates_file, dtype='str')
    
    def read_the_data():
        # Read json files for each date
        nights={}
        count=0
        for date in dates:
            path=path0+date
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
                    coo_tile_tmp.append([tl['ra'], tl['dec']]) # ra and dec are strings!
                #~ break
                    count+=1
            nights[date]={'mag': mags_tmp, 'coo': coo_tmp, 'coo_tiles': coo_tile_tmp}
            # coo_tiles: because we want percentage of the field done
            #~ break
        print 'Counted number of tiles', count
        return nights
    
    
    nights=read_the_data() 
    # TODO: ORDER nights for time simulation
    #~ print nights


    def cumulative_numbers():
        keys=[k for k in nights.iterkeys()]
        keys=sorted(keys)
        
        # magnitude distributions; for animation
        #~ fig=plt.figure()
        #~ ax=fig.add_subplot(111)
        #~ for k, v in nights.iteritems():
            #~ print k
            #~ mag=v['mag']
            #~ ax.hist(mag, histtype='step')
        
        
        # Number of stars over time. Cumulative
        number_of_stars_per_night=np.array([[int(k), len(nights[k]['coo'])] for k in keys])
        number_of_tiles_per_night=np.array([[int(k), len(nights[k]['coo_tiles'])] for k in keys])
        
        #~ for k in keys:
            #~ print k, nights[k]['coo_tiles']
        
        #~ print 'number_of_tiles_per_night', number_of_tiles_per_night
        
        #~ print number_of_stars_per_night
        n=0
        number_of_stars_per_night2=copy.deepcopy(number_of_stars_per_night)
        #~ day0=number_of_stars_per_night2[0,0] # wrong. it is not continuous.
        for i in range(len(number_of_stars_per_night)):
            #~ print number_of_stars_per_night2[i,1], n
            number_of_stars_per_night2[i,0]=i
            n1=number_of_stars_per_night2[i,1]
            number_of_stars_per_night2[i,1]+=n
            n+=n1
        number_of_stars_per_night2=np.array(number_of_stars_per_night2)
        #~ print number_of_stars_per_night2
        

        # TILES
        n=0
        number_of_tiles_per_night2=copy.deepcopy(number_of_tiles_per_night)
        for i in range(len(number_of_tiles_per_night)):
            number_of_tiles_per_night2[i,0]=i
            n1=number_of_tiles_per_night2[i,1]
            number_of_tiles_per_night2[i,1]+=n
            n+=n1
        number_of_tiles_per_night2=np.array(number_of_tiles_per_night2)

        return number_of_stars_per_night2, number_of_tiles_per_night2, number_of_stars_per_night, number_of_tiles_per_night


    number_of_stars_per_night2, number_of_tiles_per_night2, number_of_stars_per_night, number_of_tiles_per_night = cumulative_numbers()
    
    '''
    PLOT
    '''
    #~ years = mdates.YearLocator()   # every year
    #~ months = mdates.MonthLocator()  # every month
    #~ yearsFmt = mdates.DateFormatter('%Y%m')
    dates_datetime=[datetime.date(year=int(d[:4]), month=int(d[4:6]), day=int(d[6:])) for d in dates]
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    #~ ax.plot(number_of_stars_per_night2[:,0], number_of_stars_per_night2[:,1])
    ax.plot(dates_datetime, number_of_stars_per_night2[:,1])
    
    # twinxy. xlabels
    ax2=ax.twiny()
    #~ ax2.set_xticks(dates_datetime)
    #~ ax2.set_xticklabels(dates)
    ax2.set_xticks(number_of_stars_per_night2[:,0])
    #~ ax2.xaxis.set_major_locator(50)
    #~ ax2.xaxis.set_major_formatter(yearsFmt)
    #~ ax2.xaxis.set_minor_locator(10)

    #~ ax2.set_xticklabels(dates)
    #~ ax2.xaxis.set_major_locator(years)
    #~ ax2.xaxis.set_major_formatter(yearsFmt)
    #~ ax2.xaxis.set_minor_locator(months)
    
    # Add percentage of the input catalog, cumulative plot
    num_targets=15505354.0
    ax3=ax.twinx()
    ax3.plot(dates_datetime, number_of_stars_per_night2[:,1]/num_targets, c='r')
    ax3.set_ylim(0, 1.0)
    print number_of_stars_per_night2[-1,1], number_of_stars_per_night2[-1,1]/num_targets
    
    #~ ax=fig.add_subplot(212)
    #~ ax.hist(number_of_tiles_per_night, histtype='step', cumulative=True)    



    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dates_datetime, number_of_tiles_per_night2[:,1])
    ax2=ax.twiny()
    ax2.set_xticks(number_of_tiles_per_night2[:,0])
    num_tiles=14201.0
    ax3=ax.twinx()
    ax3.plot(dates_datetime, number_of_tiles_per_night2[:,1]/num_tiles, c='r')
    ax3.set_ylim(0, 1.0)
    print number_of_tiles_per_night2[-1,1], number_of_tiles_per_night2[-1,1]/num_tiles


    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dates_datetime, number_of_stars_per_night[:,1])
    ax2=ax.twinx()
    ax2.plot(dates_datetime, number_of_tiles_per_night[:,1])

    
    plt.show()

def number_of_nights_observed_each_bright_time():
    todo=True
    # histogram for each observing run

def how_priorities_drop_over_time():
    todo=True
    # read tiling file
    # compute ranking
    # read all json files and get fieldID

def how_magnitude_distribution_is_changing_over_time():
    todo=True
    # but it depends on priority list
    
if __name__ == "__main__":
    how_number_of_observed_stars_grows_with_time(path=params_simulator.params['obs_config_json_folder'], simulator_dates_file='simulator/simulator_dates_motley.dat')
