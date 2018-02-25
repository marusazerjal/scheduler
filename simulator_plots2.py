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
import pickle

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord

import params
import params_simulator

simulation_nickname = 'v1'

filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat'

# CLEAN LIST OF TILES ALREADY OBSERVED BEFORE EACH SIMULATION!
# Remote twilight time from time available to observe.

class Dictlist(dict):
	def __setitem__(self, key, value):
		try:
			self[key]
		except KeyError:
			super(Dictlist, self).__setitem__(key, [])
		self[key].append(value)


'''
TESTS OF SIMULATION
'''


'''
Repeated observations: standards are relatively short list, compared to the entired input catalog. So standard stars are going to be observed many times. And standard stars are just the usual stars (plus color cut etc.). So we can use standard stars as science targets and study their variability.

'''

try:
    TILES = data[0]
    SETTINGS = data[2]
except NameError:
    print 'Reading tiling data...'
    fl = open(params.params['input_tiling_filename'], 'rb')
    data = pickle.load(fl)
    fl.close()				

    TILES = data[0]
    SETTINGS = data[2]

print 'Number of input tiles:', len(TILES)



def determine_internal_tile_id_and_priorities():
    tiles2=[]
    for tile_id, x in enumerate(TILES):
        x.field_id=tile_id
        prior=x.calculate_tile_score(method=SETTINGS['ranking_method'], disqualify_below_min=SETTINGS['disqualify_below_min'], combined_weight=SETTINGS['combined_weight'], exp_base=SETTINGS['exp_base'])
        x.priority=prior
        tiles2.append(x)
    return tiles2


def test_if_stars_in_tiling_code_are_assigned_only_once():
    tiles=determine_internal_tile_id_and_priorities()
    
    d=Dictlist()
    t=Dictlist()
    for tile in tiles:
        for x in tile.get_assigned_targets_science(include_science_standards=False):
        #~ for x in tile.get_assigned_targets_standard():
            d[x.idn]=tile.field_id
    
    stats=[]
    for k, v in d.iteritems():
        if len(v)>1: # Gaia_ids that are assigned to more than 1 tile
            print k, v
            stats.append(len(v))

    for k, v in d.iteritems():
        if len(v)>1:
            for i in v:
                t[i]=k
    
    # number of stars per tile that were assigned to some other tile as well
    s=[]
    for k, v in t.iteritems():
        print k, len(v)
        s.append(len(v))
    print len(t)
    
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax1.hist(stats)

    ax2=fig.add_subplot(212)
    ax2.hist(s)
    plt.show()
    
    
    

def test_if_any_of_the_stars_was_observed_more_than_once():
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')

    ids=set()
    jf=set()
    
    c=0
    a=0
    repeated=Dictlist()
    dc=Dictlist()
    for x in data:
        json_filename=x[-1].replace(' ', '')
        if json_filename not in jf:
            jf.update([json_filename])
        
            with open(json_filename) as json_data:
                d = json.load(json_data)
                targets=d['targets']
                for tg in targets:
                    idt=tg['targetID']
                    dc[idt]=[x[0], x[1], tg['dec'], tg['ra']]
                    a+=1
                    if idt not in ids:
                        c+=1
                        ids.update([idt])
                    else:
                        #~ print idt
                        repeated[idt]=x[0]
                    if len(ids)%1000==0:
                        print len(ids)
    print 'Total len', len(ids)
    #~ print ids
    print c, a
    
    #~ for k, v in repeated.iteritems():
        #~ print k, v
    
    for k, v in dc.iteritems():
        if len(v)>1:
            print k, v
    

def test_if_any_of_the_fields_was_observed_more_than_once(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')

    r=Dictlist()

    for x in data:
        r[x[1]]=x[0]
    for k, v in r.iteritems():
        if len(v)>1:
            print k, v


def read_the_data(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')
    
    nights=dict()

    for x in data[:10]:
        
    
    # time; tile_id; Ntargets; priority; weight; mag_max; json_filename 
    #~ 2019-10-24 17:35; 3549; 140; 1260; 14.257791; 10.0; observers_files/funnelweb_simulation_v1/20191024/3549_173500.obs_config.json 
    #~ 2019-10-24 17:40; 3542; 140; 1260; 30.743440; 10.0; observers_files/funnelweb_simulation_v1/20191024/3542_174000.obs_config.json 
    #~ 2019-10-24 17:46; 10099; 143; 1287; 48.924326; 14.0; observers_files/funnelweb_simulation_v1/20191024/10099_174600.obs_config.json 
    #~ 2019-10-24 17:55; 15001; 143; 1287; 38.175544; 14.0; observers_files/funnelweb_simulation_v1/20191024/15001_175500.obs_config.json

        json_filename=x[-1].replace(' ', '')

        d=x[0]
        date = datetime.date(year=int(d[:4]), month=int(d[5:7]), day=int(d[8:10]))

        mags_tmp=[]
        coo_tmp=[]
        coo_tile_tmp=[]

        with open(json_filename) as json_data:
            d = json.load(json_data)
            targets=d['targets']
            for tg in targets:
                mags_tmp.append(tg['mag'])
                coo_tmp.append([tg['ra'], tg['dec']])
            
            tl=d['fieldCentre']
            coo_tile_tmp.append([tl['ra'], tl['dec']]) # ra and dec are strings!

        nights[date]={'mag': mags_tmp, 'coo': coo_tmp, 'coo_tiles': coo_tile_tmp}

 

def number_of_targets_per_night(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')

    r=Dictlist()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        r[date]=int(x[2])
        
    dates=[k for k in r.iterkeys()]
    #~ Ntiles=[len(v) for v in r.itervalues()]
    numbers=[np.sum(v) for v in r.itervalues()]
    
    print dates
    print numbers
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.scatter(dates, numbers)
    plt.show() 

def cumulative_number_of_targets_over_time(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')

    # clear out doubles
    ids=set()
    data2=[]
    for x in data:
        if x[0] in ids:
            continue
        else:
            data2.append(x)
            ids.update([x[0]])
    print len(data), len(data2)
    data=data2

    r=Dictlist()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        r[date]=int(x[2])
    dates=sorted(r)
    numbers=[np.sum(r[date]) for date in dates]
    cumulative=np.cumsum(numbers)
    ntotal=2071223.0
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dates, cumulative)
    ax2=ax.twinx()
    ax2.plot(dates, cumulative/ntotal)
    plt.show()     

def check_if_any_tile_is_observed_more_than_once(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')

    r=Dictlist()

    for x in data:
        r[x[1]]=x[0]
    for k, v in r.iteritems():
        if len(v)>1:
            print k, v
    
def are_any_stars_assigned_more_than_once():
    todo=True    

def priorities_over_time(filename_with_json_filenames=None):
    data=np.loadtxt(filename_with_json_filenames, dtype='string', comments='#', delimiter=';')
    print len(data)
    r=Dictlist()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        ax.scatter(date, int(x[3]), c='k', s=1)

    plt.show() 

    
if __name__ == "__main__":
    #~ read_the_data(filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat')
    #~ number_of_targets_per_night(filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat')
    #~ priorities_over_time(filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat')
    #~ cumulative_number_of_targets_over_time(filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat')
    
    # TESTS
    test_if_stars_in_tiling_code_are_assigned_only_once()
    #~ test_if_any_of_the_stars_was_observed_more_than_once()
    #~ check_if_any_tile_is_observed_more_than_once(filename_with_json_filenames='simulator/simulator_statistics_motley/simulator_statistics_all4.dat')
