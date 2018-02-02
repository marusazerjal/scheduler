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
from astropy import units as u
from astropy.coordinates import SkyCoord

import params
import params_simulator

# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.


print 'Input tiling file', params.params['input_tiling_filename']
#~ try:
    #~ TILES = data[0]
    #~ SETTINGS = data[2]
#~ except NameError:
    #~ print 'Reading data...'
    #~ fl = open(params.params['input_tiling_filename'], 'rb')
    #~ data = pickle.load(fl)
    #~ fl.close()				

    #~ TILES = data[0]
    #~ SETTINGS = data[2]
#~ MAX_PRIORITY=0
#~ print 'Number of input tiles:', len(TILES)


#~ def determine_internal_tile_id_and_priorities():
    #~ """
    #~ Determine internal field_id and priority.
    #~ """
    #~ tiles2=[]
    #~ priorities=[]
    #~ for tile_id, x in enumerate(TILES):
        #~ x.field_id=tile_id
        #~ prior=x.calculate_tile_score(method=SETTINGS['ranking_method'], disqualify_below_min=SETTINGS['disqualify_below_min'], combined_weight=SETTINGS['combined_weight'], exp_base=SETTINGS['exp_base'])
        #~ x.priority=prior
        #~ priorities.append(prior)
        #~ tiles2.append(x)
    #~ TILES=tiles2
    #~ MAX_PRIORITY=np.max(priorities)
#~ determine_internal_tile_id_and_priorities()
    
def how_number_of_observed_stars_grows_with_time(path=None, simulator_dates_file=None):
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    if path is None: # Take the recent results
        path0=params_simulator.params['obs_config_json_folder']
        dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    else: # results from a motley simulation in a different folder
        path0=path
        dates=np.loadtxt(simulator_dates_file, dtype='str')
    
    #~ print path0, dates
    
    def read_the_data():
        # Read json files for each date
        nights={}
        count=0
        actual_dates=[]
        for date in dates:
            path=path0+date
            try: # why does it happen that nothing is observed this night?
                json_files = [f for f in listdir(path)] # if isfile(join(mypath, f))
                actual_dates.append(date)
            except:
                continue
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
        print 'Number of nights missing in the data', len(dates)-len(actual_dates)
        return nights, actual_dates
    
    print 'Reading data...'
    nights, dates = read_the_data() 
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
    plt.title('Number of stars')
    ax=fig.add_subplot(111)
    #~ ax.plot(number_of_stars_per_night2[:,0], number_of_stars_per_night2[:,1])
    ax.plot(dates_datetime, number_of_stars_per_night2[:,1])
    print number_of_stars_per_night2[:,1]
    
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
    num_targets=2071223.0
    ax.set_ylim(0, num_targets)
    ax3=ax.twinx()
    ax3.plot(dates_datetime, number_of_stars_per_night2[:,1]/num_targets, c='r', label='percentage')
    #~ ax3.set_ylim(0, 1.0)
    print number_of_stars_per_night2[-1,1], number_of_stars_per_night2[-1,1]/num_targets
    
    #~ ax=fig.add_subplot(212)
    #~ ax.hist(number_of_tiles_per_night, histtype='step', cumulative=True)    



    #~ fig=plt.figure()
    #~ plt.title('Number of tiles')
    #~ ax=fig.add_subplot(111)
    #~ ax.plot(dates_datetime, number_of_tiles_per_night2[:,1])
    #~ ax2=ax.twiny()
    #~ ax2.set_xticks(number_of_tiles_per_night2[:,0])
    num_tiles=16517.0
    #~ ax3=ax.twinx()
    #~ ax3.plot(dates_datetime, number_of_tiles_per_night2[:,1]/num_tiles, c='r')
    #~ ax3.set_ylim(0, 1.0)
    print number_of_tiles_per_night2[-1,1], number_of_tiles_per_night2[-1,1]/num_tiles


    fig=plt.figure()
    plt.title('Number of stars and tiles per night')
    ax=fig.add_subplot(111)
    #~ ax.plot(dates_datetime, number_of_stars_per_night[:,1], label='Stars', c='k')
    ax.bar(dates_datetime, number_of_stars_per_night[:,1], align='center', color='green', ecolor='black')
    #~ ax2=ax.twinx()
    #~ ax2.plot(dates_datetime, number_of_tiles_per_night[:,1], label='Tiles', c='r')
    plt.legend()

    
    plt.show()

def survey_completeness_after_one_year_of_observations():
    todo=True

def sky_coverage_after_one_year(path=None, simulator_dates_file=None): # animation for the webpage
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    '''
    Add local meridian at midnight for each date.
    Add Moon for each date.
    '''
    
    if path is None: # Take the recent results
        path0=params_simulator.params['obs_config_json_folder']
        dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    else: # results from a motley simulation in a different folder
        path0=path
        dates=np.loadtxt(simulator_dates_file, dtype='str')
    
    def read_the_data():
        # Read json files for each date
        nights={}
        count=0
        actual_dates=[]
        
        def coordinates_string_to_float(ra, dec):
            ra0=ra.split(':')
            #~ print ra0
            ra1=float(ra0[0])
            ra2=float(ra0[1])
            ra3=float(ra0[2])
            RA=ra1+ra2/60.0+ra3/3600.0
            RA=RA*15.0
            
            de0=dec.split(':')
            de1=np.abs(float(de0[0]))
            de2=float(de0[1])
            de3=float(de0[2])
            
            #~ if de1<1.0
                
            if '-' in dec or dec.replace(' ', '')[0]=='-':
            #~ if de1<0.0:
                sign=-1.0
            else:
                sign=1.0
            DE=de1+de2/60.0+de3/3600.0
            DE=DE*sign
            
            #~ if DE>-2.0:
                #~ print dec, DE
            
            return RA, DE
        
        for i, date in enumerate(dates):
            #~ print i
            if i%10==0:
                print i, len(dates)
            path=path0+date
            try: # why does it happen that nothing is observed this night?
                json_files = [f for f in listdir(path)] # if isfile(join(mypath, f))
                actual_dates.append(date)
            except:
                continue
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
                        #~ c=SkyCoord('%s %s'%(tg['ra'], tg['dec']), unit=(u.hourangle, u.deg))
                        #~ coo_tmp.append([c.ra.value, c.dec.value])
                        #~ print c
                        #~ coo_tmp.append([tg['ra'], tg['dec']])
                        
                        ra, dec = coordinates_string_to_float(tg['ra'], tg['dec'])
                        coo_tmp.append([ra, dec])
                    
                    tl=d['fieldCentre']
                    coo_tile_tmp.append([tl['ra'], tl['dec']]) # ra and dec are strings!
                #~ break
                    count+=1
            nights[date]={'mag': mags_tmp, 'coo': np.array(coo_tmp), 'coo_tiles': np.array(coo_tile_tmp)}
            # coo_tiles: because we want percentage of the field done
            #~ break
        print 'Counted number of tiles', count
        print 'Number of nights missing in the data', len(dates)-len(actual_dates)
        return nights, actual_dates
    
    print 'Reading data...'
    nights, dates = read_the_data() 
    # TODO: ORDER nights for time simulation
    #~ print nights

    # SORT DATES
    keys=[k for k in nights.iterkeys()]
    keys=sorted(keys)
    
    # Number of stars over time. Cumulative
    data=[]
    for k in keys:
        c=nights[k]['coo']
        data.append([k, c])
    
    # All
    fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
    ax=fig.add_subplot(111)
    for x in data:
        d=x[1]
        ax.scatter(d[:,0], d[:,1], s=1, alpha=0.1, c='k')
    plt.savefig('simulator/all.png')
    
    # Animation
    for i in range(len(data)):
        fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
        ax=fig.add_subplot(111)
        for x in data[:i+1]:
            d=x[1]
            #~ date=d[0]
            ax.scatter(d[:,0], d[:,1], s=1, alpha=0.1, c='k')
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 0)
        plt.savefig('simulator/sky_coverage_over_time/nights/%d.png'%int(data[i][0]))
        print i, len(data), data[i][0]
    #~ plt.show() 
 
def allsky_animation_of_how_fields_are_observed_through_time(path=None, simulator_dates_file=None): # animation for the webpage
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    if path is None: # Take the recent results
        path0=params_simulator.params['obs_config_json_folder']
        dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    else: # results from a motley simulation in a different folder
        path0=path
        dates=np.loadtxt(simulator_dates_file, dtype='str')
    
    def read_the_data():
        # Read json files for each date
        nights={}
        count=0
        actual_dates=[]
        
        def coordinates_string_to_float(ra, dec):
            ra0=ra.split(':')
            #~ print ra0
            ra1=float(ra0[0])
            ra2=float(ra0[1])
            ra3=float(ra0[2])
            RA=ra1+ra2/60.0+ra3/3600.0
            RA=RA*15.0
            
            de0=dec.split(':')
            de1=np.abs(float(de0[0]))
            de2=float(de0[1])
            de3=float(de0[2])
            
            #~ if de1<1.0
                
            if '-' in dec or dec.replace(' ', '')[0]=='-':
            #~ if de1<0.0:
                sign=-1.0
            else:
                sign=1.0
            DE=de1+de2/60.0+de3/3600.0
            DE=DE*sign
            
            #~ if DE>-2.0:
                #~ print dec, DE
            
            return RA, DE
        
        
        folder=params_simulator.params['observing_plan_folder']
        plans=[f for f in listdir(folder)]
        plans=sorted(plans) # sorted by date
        #~ print plans
        
        json_files=[]
        for i, path in enumerate(plans):
            if i%100==0:
                print i, len(plans)
            plan=np.loadtxt(folder+path, dtype='string')
            try:
                for p in plan[:,4]:
                    p=p.replace('observers_files/funnelweb/', '')
                    json_path=params_simulator.params['obs_config_json_folder']+p
                    #~ print json_path
                    json_files.append(params_simulator.params['obs_config_json_folder']+p)
            except:
                pass
        
        print 'json files ready'


        #~ mags_tmp=[]
        coo_tmp=[]
        #~ coo_tile_tmp=[]
        coo=[]
        for i, path2 in enumerate(json_files):
            if i%100==0:
                print i, len(json_files), path2
            coo_tmp=[]
            with open(path2) as json_data:
                d = json.load(json_data)
                targets=d['targets']
                for tg in targets:
                    #~ mags_tmp.append(tg['mag'])                    
                    ra, dec = coordinates_string_to_float(tg['ra'], tg['dec'])
                    coo_tmp.append([ra, dec])
                
                #~ tl=d['fieldCentre']
                #~ coo_tile_tmp.append([tl['ra'], tl['dec']]) # ra and dec are strings!
                
            date=path2.split('/')[-2]
            coo.append([date, np.array(coo_tmp)])

        return coo
    
    print 'Reading data...'
    data = read_the_data() 
    # TODO: ORDER nights for time simulation
    #~ print nights

    # SORT DATES
    #~ keys=[k for k in nights.iterkeys()]
    #~ keys=sorted(keys)
    
    print 'Plotting...'
    
    #~ print 'Plotting all...'
    # All
    #~ fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
    #~ ax=fig.add_subplot(111)
    #~ for x in data:
        #~ d=x[1]
        #~ ax.scatter(d[:,0], d[:,1], s=1, alpha=0.1, c='k')
    #~ plt.savefig('simulator/all3.png')
    
    print 'Plotting snapshots...'
    # Animation
    for i in range(len(data)):
        fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
        ax=fig.add_subplot(111)
        for x in data[:i+1]:
            d=x[1]
            ax.scatter(d[:,0], d[:,1], s=1, alpha=0.1, c='k')
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 0)
        plt.savefig('simulator/sky_coverage_over_time/snapshot_%d_%s.png'%(i, data[i][0]))
        print i, len(data), data[i][0]

 
def ra_deg_histograms_over_time(path=None, simulator_dates_file=None):
    '''
    Statistics of observed tiles.
    (1) Cumulative plot: How number of tiles (stars) observed grows with time.
    (2) Magnitude histograms over time.
    (3) Sky coverage: histograms (similar to Taipan): number of stars per RA-DEC bins. Completeness.
    (4) Animation: observed tiles in the sky throught time
    '''
    
    if path is None: # Take the recent results
        path0=params_simulator.params['obs_config_json_folder']
        dates=np.loadtxt(params_simulator.params['simulate_dates_file'], dtype='str')
    else: # results from a motley simulation in a different folder
        path0=path
        dates=np.loadtxt(simulator_dates_file, dtype='str')
    
    def read_the_data():
        # Read json files for each date
        nights={}
        count=0
        actual_dates=[]
        
        def coordinates_string_to_float(ra, dec):
            ra0=ra.split(':')
            #~ print ra0
            ra1=float(ra0[0])
            ra2=float(ra0[1])
            ra3=float(ra0[2])
            RA=ra1+ra2/60.0+ra3/3600.0
            RA=RA*15.0
            
            de0=dec.split(':')
            de1=np.abs(float(de0[0]))
            de2=float(de0[1])
            de3=float(de0[2])
            
            #~ if de1<1.0
                
            if '-' in dec or dec.replace(' ', '')[0]=='-':
            #~ if de1<0.0:
                sign=-1.0
            else:
                sign=1.0
            DE=de1+de2/60.0+de3/3600.0
            DE=DE*sign
            
            #~ if DE>-2.0:
                #~ print dec, DE
            
            return RA, DE
        
        for i, date in enumerate(dates):
            #~ print i
            if i%10==0:
                print i, len(dates)
            path=path0+date
            try: # why does it happen that nothing is observed this night?
                json_files = [f for f in listdir(path)] # if isfile(join(mypath, f))
                actual_dates.append(date)
            except:
                continue
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
                        #~ c=SkyCoord('%s %s'%(tg['ra'], tg['dec']), unit=(u.hourangle, u.deg))
                        #~ coo_tmp.append([c.ra.value, c.dec.value])
                        #~ print c
                        #~ coo_tmp.append([tg['ra'], tg['dec']])
                        
                        ra, dec = coordinates_string_to_float(tg['ra'], tg['dec'])
                        coo_tmp.append([ra, dec])
                    
                    tl=d['fieldCentre']
                    coo_tile_tmp.append([tl['ra'], tl['dec']]) # ra and dec are strings!
                #~ break
                    count+=1
            nights[date]={'mag': mags_tmp, 'coo': np.array(coo_tmp), 'coo_tiles': np.array(coo_tile_tmp)}
            # coo_tiles: because we want percentage of the field done
            #~ break
        print 'Counted number of tiles', count
        print 'Number of nights missing in the data', len(dates)-len(actual_dates)
        return nights, actual_dates
    
    print 'Reading data...'
    nights, dates = read_the_data() 
    # TODO: ORDER nights for time simulation
    #~ print nights

    # SORT DATES
    keys=[k for k in nights.iterkeys()]
    keys=sorted(keys)
    
    # Number of stars over time. Cumulative
    data=[]
    for k in keys:
        c=nights[k]['coo']
        data.append([k, c])
    
    #~ print 'Reading data...'
    #~ data = read_the_data() 
    # TODO: ORDER nights for time simulation
    #~ print nights

    # SORT DATES
    #~ keys=[k for k in nights.iterkeys()]
    #~ keys=sorted(keys)
    
    print 'Plotting...'
    
    #~ print 'Plotting all...'
    # All
    #~ fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
    #~ ax=fig.add_subplot(111)
    #~ for x in data:
        #~ d=x[1]
        #~ ax.scatter(d[:,0], d[:,1], s=1, alpha=0.1, c='k')
    #~ plt.savefig('simulator/all3.png')
    
    print 'Plotting snapshots...'
    # Animation
    fig=plt.figure()
    for i in range(len(data)):
        #~ fig=plt.figure()
        fig.clf()
        axra=fig.add_subplot(211)
        axde=fig.add_subplot(212)
        tmp=[]
        for x in data[:i+1]:
            d=x[1]
            for y in d:
                tmp.append(y)
        tmp=np.array(tmp)
        #~ print tmp
        try:
            axra.hist(tmp[:,0], bins=np.linspace(0, 360, 100))
            axde.hist(tmp[:,1], bins=np.linspace(-90, 10, 100))
            axra.set_xlim(0, 360)
            axde.set_xlim(-90, 10)
            plt.savefig('simulator/sky_coverage_over_time/radec/snapshot_radec_%d_%s.png'%(i, data[i][0]))
            print i, len(data), data[i][0]
        except:
            pass
 

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
    #~ how_number_of_observed_stars_grows_with_time(path=params_simulator.params['obs_config_json_folder'], simulator_dates_file='simulator/simulator_dates_motley.dat')
    how_number_of_observed_stars_grows_with_time()
    
    #~ sky_coverage_after_one_year() # this one
    #~ allsky_animation_of_how_fields_are_observed_through_time()
    
    #~ ra_deg_histograms_over_time()
