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
import collections

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord

import params
reload(params)
import params_simulator
reload(params_simulator)
import scheduler
reload(scheduler)

class Dictlist(dict):
	def __setitem__(self, key, value):
		try:
			self[key]
		except KeyError:
			super(Dictlist, self).__setitem__(key, [])
		self[key].append(value)

################################

root='data8'
root='test3_julij4'

print 'root:', root

filename_with_json_filenames='%s/simulator_statistics.dat'%root
filename_data_volume='%s/simulator_statistics_calibration.dat'%root

################################

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

def data_volume_with_time():
    '''
    How data volume increases with time. Science shots, arc, flat, bias. 2 arms.
    '''
    
    fits_size = 34.0 # MB
    
    d=np.loadtxt(filename_data_volume, dtype='string', delimiter=';')
    calib=Dictlist()
    for x in d:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        n=int(x[1])
        calib[date]=n
 
    tls=Dictlist()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tls[date]=2 # 2 arms --> 2 fits files per tile
    
    dates=sorted(tls)
    
    #~ print tls[dates[0]]
    
    #~ science=[np.sum(tls[date]) for date in dates]
    #~ calib=[np.sum(calib[date]) for date in dates]
    
    #~ print science
    
    total=[np.sum(tls[date])+np.sum(calib[date]) for date in dates]    
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.bar(dates, total)
    ax2=ax.twinx()
    ax2.plot(dates, np.cumsum(total)*fits_size / 1024.0 / 1024.0, c='k') # [TB], 34 MB per fits file
    mean=np.mean(total)
    ax.axhline(y=mean, c='r', label='Mean = %d fits files (%d GB)'%(mean, mean*fits_size / 1024.0))
    ax.legend(loc=2)
    ax.set_ylabel('Number of fits files per night')
    ax2.set_ylabel('Data volume [TB]')
    plt.show() 

def number_of_targets_per_night():
    '''
    Total number over time, plus number of stars with high priorities over time.
    '''
    '''
    Duration of the survey.
    Number of UNIQUE targets.
    '''
    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    r=Dictlist()
    tls=Dictlist()
    
    unique=set()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
        idn=set([y.idn for y in tile.get_assigned_targets_science()])
        new_unique_targets=idn.difference(unique)
        unique.update(idn)       
        r[date]=len(new_unique_targets)
        tls[date]=tileid
    dates=sorted(r)
    numbers=[np.sum(r[date]) for date in dates]

    numbers_tiles=[len(tls[date]) for date in dates]
        
    fig=plt.figure()
    ax=fig.add_subplot(111)
    #~ ax.scatter(dates, numbers)
    ax2=ax.twinx()
    ax2.plot(dates, numbers_tiles, c='k')
    ax.bar(dates, numbers)
    mean=np.mean(numbers)
    ax.axhline(y=mean, c='r', label='Mean = %d'%mean)
    ax.legend(loc=1)
    plt.show() 

def percentage_of_local_area_covered_magnitude_ranges():
    '''
    Snapshots for each magnitude bin.
    '''

    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    # tiles per magnitude range
    tiles_mags=Dictlist()
    for x in TILES:
        tiles_mags[x.mag_max]=x.field_id
    
    mags=sorted(tiles_mags, reverse=True)

    mag_ranges=s.magnitude_ranges()
    magranges={8.5: mag_ranges[0], 10.5: mag_ranges[1], 12.5: mag_ranges[2], 14.5: mag_ranges[3]}
    
    fig2=plt.figure(figsize=(16.0, 10.0), dpi=600)
    ncols=4
    nrows=4
    
    #~ ii=0
    for ii, mag in enumerate(mags):
        tls=set(tiles_mags[mag])
        
        # add magnitude range!
        COO=Dictlist()
        for x in TILES:
            if x.mag_max==mag:
                COO[(x.ra, x.dec)]=1
        
        density_date={}
        observed=Dictlist()

        nights=Dictlist()
        for x in data:
            date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
            tileid=int(x[1])
            if tileid in tls: 
                nights[date]=tileid

        dates=sorted(nights)


        #~ i=0
        for date in dates:
            tileids=nights[date]
        #~ for date, tileids in nights.iteritems(): # all tiles on this date
            #~ i+=1
            #~ if i>100:
                #~ break
            
            for tileid in tileids:
                tile=tiles[tileid]
                c=(tile.ra, tile.dec)
                
                observed[c]=1
            
            d=[]
            for k, v in observed.iteritems():
                n=float(len(v))
                N=float(len(COO[k]))
                
                dens=n/N
                d.append([k, dens])
                
            density_date[date]=d

        dates=sorted(density_date)

        i=0
        #~ for k in dates:
        if ii==0:
            l=float(len(dates))
            #~ snapshot_dates=[dates[int(l/4.0)], dates[int(l/2.0)], dates[int(3.0*l/4.0)], dates[-1]]
            snapshot_dates=[dates[int(l/16.0)], dates[int(l/8.0)], dates[int(l/4.0)], dates[int(l/2.0)]]
            print snapshot_dates
        #~ ii=+1
        
        for kk, k in enumerate(snapshot_dates):
            try:
                v=density_date[k]
            except:
                last_date=dates[-1]
                if k>last_date:
                #~ print mag, kk, k, dates[-1]
                    v=density_date[dates[-1]] # some magnitude ranges are done in the first half a year and dont have later dates
                else:
                    index=np.abs([k-d for d in dates]).argmin()
                    v=density_date[dates[index]]
                    print mag, kk, k, dates[-1], dates[index]
            #~ fig=plt.figure()
            #ax=fig.add_subplot(111, projection='aitoff')
            #~ ax=fig.add_subplot(111)
            #~ ax.set_title(str(k)+',  mag '+str(magranges[mag]))
            #p=np.array([[np.deg2rad(x[0][0]), np.deg2rad(x[0][1]), x[1]] for x in v])
            #~ p=np.array([[x[0][0], x[0][1], x[1]] for x in v])
            #~ cb=ax.scatter(p[:,0], p[:,1], c=p[:,2], vmin=0, vmax=1)
            #~ plt.colorbar(cb)
            #plt.clim(0, 1)
            #~ ax.set_xlim(0, 360)
            #~ ax.set_ylim(-90, 10)
            #~ plt.savefig('density_mag/fig_%d_%d.png'%(mag, i))


            #~ ax=fig.add_subplot(111, projection='aitoff')
            idx=nrows*ncols-((ii+1)*4-kk)+1
            ax=fig2.add_subplot(ncols, nrows, idx)
            if ii==3:
                ax.set_title(str(k))
            if kk==0:
                ax.set_ylabel(str(magranges[mag]))
            #~ p=np.array([[np.deg2rad(x[0][0]), np.deg2rad(x[0][1]), x[1]] for x in v])
            p=np.array([[x[0][0], x[0][1], x[1]] for x in v])


            #~ RA=p[:,0]
            
            #~ org=0
            #~ x = np.remainder(RA+360-org,360) # shift RA values
            #~ ind = x>180
            #~ x[ind] -=360    # scale conversion to [-180, 180]
            #~ x=-x    # reverse the scale: East to the left
            #~ RA=x
            #~ tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
            #~ tick_labels = np.remainder(tick_labels+360+org,360)
            #~ ax.scatter(np.radians(RA),np.radians(p[:,1]), c=p[:,2], vmin=0, vmax=1, s=5)  # convert degrees to radians
            
            
            
            
            cb=ax.scatter(p[:,0], p[:,1], c=p[:,2], vmin=0, vmax=1, s=5)
            #~ plt.colorbar(cb)
            #~ plt.clim(0, 1)
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 10)
            

            i+=1

    fig2.savefig('density_mag/all.png')
def percentage_of_local_area_covered_magnitude_ranges_aitoff():
    '''
    Snapshots for each magnitude bin.
    '''

    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    # tiles per magnitude range
    tiles_mags=Dictlist()
    for x in TILES:
        tiles_mags[x.mag_max]=x.field_id
    
    mags=sorted(tiles_mags, reverse=True)

    mag_ranges=s.magnitude_ranges()
    magranges={8.5: mag_ranges[0], 10.5: mag_ranges[1], 12.5: mag_ranges[2], 14.5: mag_ranges[3]}
    
    fig2=plt.figure(figsize=(16.0, 10.0), dpi=600)
    ncols=4
    nrows=4
    
    for ii, mag in enumerate(mags):
        tls=set(tiles_mags[mag])
        
        # add magnitude range!
        COO=Dictlist()
        for x in TILES:
            if x.mag_max==mag:
                COO[(x.ra, x.dec)]=1
        
        density_date={}
        observed=Dictlist()

        nights=Dictlist()
        for x in data:
            date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
            tileid=int(x[1])
            if tileid in tls: 
                nights[date]=tileid

        dates=sorted(nights)

        for date in dates:
            tileids=nights[date]
           
            for tileid in tileids:
                tile=tiles[tileid]
                c=(tile.ra, tile.dec)
                
                observed[c]=1
            
            d=[]
            for k, v in observed.iteritems():
                n=float(len(v))
                N=float(len(COO[k]))
                
                dens=n/N
                d.append([k, dens])
                
            density_date[date]=d

        dates=sorted(density_date)

        i=0
        if ii==0:
            l=float(len(dates))
            snapshot_dates=[dates[int(l/16.0)], dates[int(l/8.0)], dates[int(l/4.0)], dates[int(l/2.0)]]
            print snapshot_dates
        
        for kk, k in enumerate(snapshot_dates):
            try:
                v=density_date[k]
            except:
                last_date=dates[-1]
                if k>last_date:
                #~ print mag, kk, k, dates[-1]
                    v=density_date[dates[-1]] # some magnitude ranges are done in the first half a year and dont have later dates
                else:
                    index=np.abs([k-d for d in dates]).argmin()
                    v=density_date[dates[index]]
                    print mag, kk, k, dates[-1], dates[index]


            idx=nrows*ncols-((ii+1)*4-kk)+1
            ax=fig2.add_subplot(ncols, nrows, idx, projection='aitoff')
            if ii==3:
                ax.set_title(str(k))
            if kk==0:
                ax.set_ylabel(str(magranges[mag]))
            #~ p=np.array([[np.deg2rad(x[0][0]), np.deg2rad(x[0][1]), x[1]] for x in v])
            p=np.array([[x[0][0], x[0][1], x[1]] for x in v])

            cb=ax.scatter(np.radians(p[:,0]-180.0), np.radians(p[:,1]), c=p[:,2], vmin=0, vmax=1, s=5)
            #~ plt.colorbar(cb)
            #~ plt.clim(0, 1)
          

            i+=1

    fig2.savefig('figs/testall.png')
def percentage_of_local_area_covered_magnitude_ranges_aitoff_all_stars_animation():
    '''
    Snapshots for each night. For all stars together
    '''

    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}

    COO0=Dictlist()
    for x in TILES:
        COO0[(x.ra, x.dec)]=1
    
    COO={k: float(len(v)) for k, v in COO0.iteritems()}

    nights=Dictlist()
    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        nights[date]=tiles[tileid]
    
    dates=sorted(nights)

    test=[]
    for date in dates:
        v=nights[date]
        test.append([x.ra/15.0 for x in v])

    fig=plt.figure()
    ax=fig.add_subplot(111)
    for date, x in zip(dates, test):
        ax.scatter([date for i in range(len(x))], x, c='k', alpha=0.05)

    plt.show()
        
    #~ coo={}
    #~ observed=Dictlist()
    #~ for date in dates:
        #~ v=nights[date]
        #~ # update night with observed tiles
        #~ for x in v:
            #~ observed[(x.ra, x.dec)]=1
        
        #~ # density in the morning
        #~ tmp=[]
        #~ for c, l in observed.iteritems():
            #~ n_total=COO[c]
            #~ n_observed=float(len(l))
            #~ d = n_observed / n_total
            #~ tmp.append([c[0], c[1], d])
        #~ coo[date]=np.array(tmp)
    
    #~ print len(dates)
    
    #~ for i, date in enumerate(dates):
        #~ fig=plt.figure(figsize=(16.0, 10.0), dpi=600)
        #~ ax=fig.add_subplot(111, projection='aitoff')
        #~ ax.set_title(date)
        #~ p=coo[date]
        #~ cb=ax.scatter(np.radians(p[:,0]-180.0), np.radians(p[:,1]), c=p[:,2], vmin=0, vmax=1, s=10)
        #~ plt.colorbar(cb)
        #~ fig.savefig('figs/animation/density_%d.png'%i)
        
        #~ if i>10:
            #~ break


def repeated_observations():
    '''
    Standard stars are observed more than once.
    '''

    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    r=Dictlist()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
                
        for y in tile.get_assigned_targets_science():
            r[y.idn]=date


    repeats=[len(v) for v in r.itervalues()]

    
    c=collections.Counter(repeats)
    print c

    p=np.array([[k, v] for k, v in c.iteritems()])



    #~ r=[]
    #~ m=0
    #~ for k, v in d.iteritems():
        #~ l=len(v)
        #~ if l>1:
            #~ r.append(l)
        #~ if l>m:
            #~ m=l
    #~ print 'Number of stars with more than 1 observation:', len(r)
    #~ print 'Number of observations: max', m

    fig=plt.figure()
    ax=fig.add_subplot(111)
    #~ ax.hist(r, bins=np.arange(22)-0.5)
    ax.bar(p[:,0], p[:,1])
    ax.set_yscale('log')
    #~ ax.set_xticks(range(22))
    #~ ax.set_xlim([1, 16])
    #~ ax.set_ylim([0.5, 10000])
    plt.show()

def priorities_and_magnitudes_over_time():
    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    priorities=Dictlist()
    for x in data:
        tileid=int(x[1])
        tile=tiles[tileid]
        p=[y.priority for y in tile.get_assigned_targets_science(include_science_standards=False)]
        for y in p:
            priorities[y]=1
    
    priorities=Dictlist()
    for x in TILES:
        p=[y.priority for y in x.get_assigned_targets_science(include_science_standards=False)]
        for y in p:
            priorities[y]=1
    priorities_total={k: len(v) for k, v in priorities.iteritems()}
    print priorities_total
    
    r=Dictlist()
    #~ pr=[[] for x in range(6)]
    for x in data:
        pr=[[] for xx in range(6)]
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
        p=[y.priority for y in tile.get_assigned_targets_science(include_science_standards=False)]
        
        for z in p:
            pr[z].append(1)
        
        summary=[len(z) for z in pr]
        #~ print summary
        r[date]=summary
    
    final={}    
    for k, v in r.iteritems():
        v=np.array(v)
        vv=v.reshape(len(v), 6)
        #~ print k, vv, np.sum(vv, axis=0)
        final[k]=np.sum(vv, axis=0)
        
    dates=sorted(final)
    f=[]
    for date in dates:
        num=final[date]
        f.append(num)
    f=np.array(f)

    
    fig=plt.figure()
    ax=fig.add_subplot(211)
    ax.axvline(x=365.0/2.0, lw=0.5, c='grey')
    ax.axvline(x=365.0, lw=0.5, c='grey')
    ax.axvline(x=365.0*1.5, lw=0.5, c='grey')
    ax.axvline(x=365.0*2.0, lw=0.5, c='grey')    
    ax.axhline(y=0.5, lw=0.5, c='grey') 
    colors=['k', 'r', 'blue', 'green', 'orange', 'purple']
    for i in range(6):
        try:
            cumulative=np.cumsum(f[:,i])
            n=float(priorities_total[i])
            #~ print cumulative
            ax.plot(dates, cumulative/n, c=colors[i], label=i)
        except:
            print 'NOT POSSIBLE', i
            pass
    ax.legend(loc=4)
    ax.set_xlim(dates[0], dates[-1])
    ax.set_ylabel('Fraction of stars observed')
    ax.xaxis.tick_top()
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    #~ plt.show()  
    
    
    '''
    How magnitude ranges get completed over time. It seems that the faintest get completed later (probably because they have the highest number of tiles).
    '''
    
    mag_ranges=s.magnitude_ranges()
    
    # number of all stars per each magnitude bin
    MAGS0=Dictlist()
    unique=set()
    for x in TILES:
        idn=set([y.idn for y in x.get_assigned_targets_science()])
        new_unique_targets=idn.difference(unique)
        unique.update(idn)

        for y in x.get_assigned_targets_science():
            if y.idn in new_unique_targets:
                m=y.mag
                for i, mr in enumerate(mag_ranges):
                    if m>=mr[0] and m<mr[1]:
                        k=mr
                        MAGS0[k]=m

    mranges=sorted(MAGS0)
    MAGS=[float(len(MAGS0[k])) for k in mranges]
    
    print 'MAGS', MAGS

    
    r=Dictlist()
    r_all=Dictlist()
    unique=set()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
        idn=set([y.idn for y in tile.get_assigned_targets_science()])
        new_unique_targets=idn.difference(unique)
        unique.update(idn)
        
        r_all[date]=len(new_unique_targets)
        
        mags=[[] for i in range(len(mag_ranges))]
        for y in tile.get_assigned_targets_science():
            if y.idn in new_unique_targets:
                m=y.mag
                for i, mr in enumerate(mag_ranges):
                    if m>=mr[0] and m<mr[1]:
                        mags[i].append(m)

        mags_summary=[float(len(m)) for m in mags]
        r[date]=mags_summary
    
    dates=sorted(r)
    
    t0=dates[0]
    days=[int((x-t0).days) for x in dates]
    
    numbers_all=[np.sum(r_all[date]) for date in dates]
    cumulative=np.cumsum(numbers_all)
    ntotal=1963635.0

    
    # Convert lists to np.arrays for every night
    nmr=len(mag_ranges)
    r2=[[] for i in range(nmr)]
    for date in dates:
        v=r[date]
        v=np.array(v)
        d=v.reshape(len(v), nmr)
        for i in range(nmr):
            s=np.sum(d[:,i])
            r2[i].append(s)
          
    # PLOT AX2
    colors=['blue', 'r', 'green', 'purple']
    ax=fig.add_subplot(212)
    ax.axvline(x=365.0/2.0, lw=0.5, c='grey')
    ax.axvline(x=365.0, lw=0.5, c='grey')
    ax.axvline(x=365.0*1.5, lw=0.5, c='grey')
    ax.axvline(x=365.0*2.0, lw=0.5, c='grey')    
    ax.axhline(y=0.5, lw=0.5, c='grey')    
    for i in range(len(mag_ranges)):
        numbers=np.cumsum(r2[i])/MAGS[i]
        ax.plot(days, numbers, c=colors[i], label=mag_ranges[i])
    ax.plot(days, cumulative/ntotal, c='k', label='All stars')
    l=ax.legend(loc=4, frameon=False)
    l.get_frame().set_facecolor('white')
    ax.set_xlim(0, days[-1])
    ax.set_xlabel('Duration of the survey [days]') # These are calendar days (including days with dark time when Taipan is observing)
    ax.set_ylabel('Fraction of stars observed')
    print 'now show:'
    plt.show()      
    
def obstile_weights_versus_time():
    '''
    The last tiles to be observed (what fraction) have weights equal to 0.
    '''
    r=Dictlist()
    for x in data:
        pr=[[] for xx in range(6)]
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        w=float(x[3])
        r[date]=w
    
    dates=sorted(r)
    #~ weights=[r[x] for x in dates]
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for date in dates:
        v=r[date]
        ax.scatter([date for i in range(len(v))], v, c='k', alpha=0.01)
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1e+10)
    plt.show()
    
def ranking_versus_time():
    '''
    The last tiles to be observed (what fraction) have weights equal to 0.
    '''
    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}

    r=Dictlist()
    for x in data:
        pr=[[] for xx in range(6)]
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
        r[date]=tile.priority
    
    dates=sorted(r)
    #~ weights=[r[x] for x in dates]
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for date in dates:
        v=r[date]
        ax.scatter([date for i in range(len(v))], v)
    ax.set_yscale('log')
    plt.show()

def ra_every_day():
    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    r=Dictlist()
    tls=Dictlist()
    
    unique=set()

    for x in data:
        date = datetime.date(year=int(x[0][:4]), month=int(x[0][5:7]), day=int(x[0][8:10]))
        tileid=int(x[1])
        tile=tiles[tileid]
        r[date]=[tile.ra, tile.dec]

    dates=sorted(r)
    p=[]
    for date in dates:
        d=r[date]
        d=np.array(d)
        #~ p.append([np.min(d[:,0]), np.max(d[:,0])])
        p.append([d[0,0], d[-1,0]])
    p=np.array(p)
         
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(dates, p[:,0]/15.0, c='k')
    ax.plot(dates, p[:,1]/15.0, c='r')
    ax.scatter(dates, p[:,0]/15.0, c='k')
    ax.scatter(dates, p[:,1]/15.0, c='r')
    plt.show()     

def ra_tonight():
    '''
    Total number over time, plus number of stars with high priorities over time.
    '''
    '''
    Duration of the survey.
    Number of UNIQUE targets.
    '''
    s=scheduler.Scheduler()
    TILES=s.tiles # with priorities and tile_ids
    tiles={x.field_id: x for x in TILES}
    
    #~ data=np.loadtxt('test1/observing_plan_20180304.dat', dtype='string')
    data=np.loadtxt('test2_julij4/observing_plan_20180304.dat', dtype='string')
    
    d=[]
    for x in data:
        print x[3]
        tileid=int(x[3])
        tile=tiles[tileid]
        
        d.append([tile.ra, tile.dec])

    d=np.array(d)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(d[:,0]/15.0, d[:,1], c='k')
    plt.show() 



if __name__ == "__main__":
    #~ number_of_targets_per_night()
    #~ repeated_observations()
    #~ priorities_and_magnitudes_over_time()
    #~ percentage_of_local_area_covered_magnitude_ranges()
    #~ percentage_of_local_area_covered_magnitude_ranges_aitoff()
    #~ percentage_of_local_area_covered_magnitude_ranges_aitoff_all_stars_animation()
    
    #~ ra_tonight()
    ra_every_day()
    
    #~ obstile_weights_versus_time()
    #~ ranking_versus_time()
    
    #~ data_volume_with_time()
