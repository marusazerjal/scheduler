"""
FunnelWeb Scheduler: Decide which tile to observe next, based on sidereal time and observing conditions (weather, seeing etc.).
"""

import numpy as np
import pickle
import datetime
from collections import defaultdict
import matplotlib.pyplot as plt

from sklearn.neighbors import KDTree

import taipan.core as tp

# FunnelWeb Scheduler
import manage_list_of_observed_tiles
import params
import obstile
import visualization




"""
Load the tiling data
"""
reload(manage_list_of_observed_tiles)
reload(params)
reload(obstile)
reload(visualization)

print 'Input tiling file', params.params['input_tiling_filename']
try:
    TILES = data[0]
    SETTINGS = data[2]
except NameError:
    print 'Reading data...'
    fl = open(params.params['input_tiling_filename'],'rb')
    data = pickle.load(fl)
    fl.close()				

    TILES = data[0]
    SETTINGS = data[2]

print 'Number of input tiles:', len(TILES)

#~ TILES=TILES[:100]

# Distribute into magnitude ranges
from collections import defaultdict
tiles_mag_range = defaultdict(list)
for i, x in enumerate(TILES):
    tiles_mag_range[tuple([float(x.mag_min), float(x.mag_max)])].append([i, np.deg2rad(x.ra), np.deg2rad(x.dec)])
tiles_mag_range={k: v for k, v in tiles_mag_range.iteritems()}

#~ for k, v in tiles_mag_range.iteritems():
    #~ print k, v

#~ X=np.array([[np.deg2rad(x.ra), np.deg2rad(x.dec), i] for i, x in enumerate(TILES)])
#~ X=np.array([[np.deg2rad(x.ra), np.deg2rad(x.dec)] for x in TILES])
X=np.array(tiles_mag_range[(5.0, 8.0)])
print 'len X', len(X)

def distance_between_two_points_in_the_sky(v1, v2):
    """
    Determine the distance between two points in the sky.
    
    Parameters
    ----------
    alpha1, delta1: float, degrees
        RA and Dec for object 1
    alpha2, delta2: float, degrees
        RA and Dec for object 2
    
    Output
    ----------
    A: float, degrees
        Distance between two points in the sky.
    """
    a1=v1[0]
    d1=v1[1]
    a2=v2[0]
    d2=v2[1]
    
    cos_A = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(a1-a2)

    #~ if cos_A<-1.0 or cos_A>1.0:
        #~ A=np.arccos(cos_A) # [0, 180] in radians
        #~ print 'invalid', cos_A, v1, v2, A
    if cos_A>-1.0 and cos_A<1.0: # There are troubles with cos_A==1
        A=np.arccos(cos_A) # [0, 180] in radians
    else:
        A=0.0
    #~ A=np.rad2deg(A)
    return A # radians

#~ tree = KDTree(X[:2], metric='haversine', leaf_size=50)              
#~ tree = KDTree(X[:2], metric=lambda v1, v2: distance_between_two_points_in_the_sky(v1, v2), leaf_size=50)              
#~ dist, ind = tree.query(X[0,:2], k=3)       
#~ print dist
#~ print ind

def knn():
    import sklearn
    from sklearn.neighbors import NearestNeighbors


    knn=NearestNeighbors(n_neighbors=50,
                     leaf_size=50,
                     algorithm='auto',
                     metric=lambda v1, v2: distance_between_two_points_in_the_sky(v1, v2),
                     n_jobs=8
                     )
    knn.fit(X)
    distances, indices = knn.kneighbors(X)
    print 'done'

    N=50

    c=np.array([[TILES[i].ra, TILES[i].dec] for i in indices[N]])
    f=open('tt.dat', 'wb')
    for x in c:
        #~ print x[0], x[1]
        f.write('%g %g\n'%(x[0], x[1]))
    f.write('%g %g\n'%(TILES[N].ra, TILES[N].dec))
    f.close()
    #~ print
    print [TILES[N].ra], [TILES[N].dec]
    #~ fig=plt.figure()
    #~ ax=fig.add_subplot(111)
    #~ ax.scatter(c[:,0], c[:,1], c='k')
    #~ ax.scatter([TILES[0].ra], [TILES[0].dec], c='red')
    #~ plt.show()

    #~ params.params['nearest_neighbours_filename']

    #~ for x in indices:
        #~ print x[0]

def brute_force(N=100):
    """
    Save first N nearest neighbours.
    """    
    result={}

    for mag_range, X in tiles_mag_range.iteritems():
        X=np.array(X)
        print 'len X', len(X)
        
        t_start=datetime.datetime.now()
        d=np.array([[x[0], y[0], distance_between_two_points_in_the_sky([x[1], x[2]], [y[1], y[2]])] for x in X for y in X if np.abs(y[0]-x[0])>1])
        t_end=datetime.datetime.now()
        dt=t_end-t_start
        print mag_range, 'Total time: %d s (%.2f min)'%(dt.seconds, dt.seconds/60.0)

        dist = defaultdict(list)
        for x in d:
            dist[int(x[0])].append([int(x[1]), np.rad2deg(x[2])])
        dist_mag_range={k: sorted(v, key=lambda x: x[1]) for k, v in dist.iteritems()}

        for k, v in dist_mag_range.iteritems():
            #~ print k, v[:30]
            result[k]=v[:N] # Save only first N neighbours

    f = open(params.params['nearest_neighbours_filename'], 'wb')
    pickle.dump(result, f)
    f.close()

if __name__ == "__main__":
    brute_force()
