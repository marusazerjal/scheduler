"""
Test tiling
"""

import numpy as np
import matplotlib.pyplot as plt

import taipan.core as tp

print 'Input tiling file', params.params['input_tiling_filename']
try:
    TILES = data[0]
    SETTINGS = data[2]
except NameError:
    print 'Reading data...'
    fl = open(params.params['input_tiling_filename'], 'rb')
    data = pickle.load(fl)
    fl.close()				

    TILES = data[0]
    SETTINGS = data[2]

print 'Number of input tiles:', len(TILES)

#~ targets=TILES[5240].get_assigned_targets_science()

t=[]
for x in TILES:
    targets=x.get_assigned_targets_science()
    r=[y.ra for y in targets]
    d=[y.dec for y in targets]
    t.append([np.max(r)-np.min(r), np.max(d)-np.max(d), x.ra, x.dec])
#~ targets=TILES[5673].get_assigned_targets_science()
#~ d=np.array([[x.ra, x.dec] for x in targets])

t=np.array(t)

fig=plt.figure()
ax=fig.add_subplot(111)
#~ ax.scatter(d[:,0], d[:,1])
ax.scatter(t[:,2], t[:,3], c=t[:,0], vmin=5, vmax=15)
#~ ax.hist(t[:,0], bins=100)

plt.show()
