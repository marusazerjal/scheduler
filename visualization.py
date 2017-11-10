"""
FunnelWeb Scheduler: Decide which tile to observe next, based on sidereal time and observing conditions (weather, seeing etc.).
"""

import numpy as np
import matplotlib.pyplot as plt

# FunnelWeb
import params

reload(params)

def plot_selected_tile_with_neighbourhood(moon=None, lst=None, best_tiles=None, tiles=None, best_tile=None, i=None, ra_current=None, dec_current=None, telescope_positions=None, observed_tile_ids=None):
    lst=lst*15.0
    if lst>180:
        lst=lst-360.0
    plot_ra_amp=50

    #~ print 'plot LST', lst

    w=np.array([[x.TaipanTile.ra, x.TaipanTile.dec, x.weight] for x in best_tiles if x.hour_angle<6.0])
    w=np.array([x if x[0]<180 else [x[0]-360.0, x[1], x[2]] for x in w])
    w=np.array(sorted(w, key=lambda y: y[2]))
    
    
    telescope_positions=np.array([[x[0] if x[0]<180 else x[0]-360.0, x[1]] for x in telescope_positions])

    # Tiles already observed
    observed=np.array([[x.ra, x.dec] for x in tiles if x.field_id in observed_tile_ids])
    observed=np.array([x if x[0]<180 else [x[0]-360.0, x[1]] for x in observed])

    tiles=np.array([[x.ra, x.dec] for x in tiles])
    tiles=np.array([x if x[0]<180 else [x[0]-360.0, x[1]] for x in tiles])


   
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.axvline(x=lst)
    ax.axvline(x=lst-15.0, alpha=0.5)
    ax.axvline(x=lst+15.0, alpha=0.5)
    
    # All candidate tiles
    ax.scatter(tiles[:,0], tiles[:,1], c='grey', alpha=0.1)
    
    m=[moon.ra.value if moon.ra.value<180.0 else moon.ra.value-360.0, moon.dec.value]
    ax.scatter(m[0], m[1], c='red', s=200)
    
    # All considered tiles
    cb=ax.scatter(w[:,0], w[:,1], c=w[:,2], vmin=0, vmax=1, cmap=plt.cm.RdYlGn_r)
    plt.colorbar(cb)
    
    # Plot tiles already observed
    ax.scatter(observed[:,0], observed[:,1], c='k')
    
    # Current telescope position
    ctp=[ra_current if ra_current<180.0 else ra_current-360.0, dec_current]
    #~ ax.scatter(ctp[0], ctp[1], c='None', edgecolor='k', s=150)
    ax.scatter(telescope_positions[:,0], telescope_positions[:,1], c='None', edgecolor='k', s=150)
    ax.plot(telescope_positions[:,0], telescope_positions[:,1], c='k')

    
    bt=[best_tile.TaipanTile.ra if best_tile.TaipanTile.ra<180.0 else best_tile.TaipanTile.ra-360.0, best_tile.TaipanTile.dec]
    #~ ax.scatter([best_tile.TaipanTile.ra], [best_tile.TaipanTile.dec], c='None', edgecolor='r', s=150)
    ax.scatter(bt[0], bt[1], c='None', edgecolor='r', s=150)  
    
    
    ax.set_xlim([lst-plot_ra_amp, lst+plot_ra_amp])
    ax.set_xlabel('RA [deg, -180:180]')
    ax.set_ylabel('DEC [deg]')
    #~ plt.show()
    plt.savefig('seq/%3d.png'%i)
    plt.close()  

def test_tile_distribution_in_the_sky():
    """
    Plot RA, Dec of all the candidate tiles from the tiling file.
    """
    fl = open(params.params['input_tiling_filename'],'rb')
    import pickle
    tiles = pickle.load(fl)[0]
    fl.close()
    print 'Number of tiles:', len(tiles)
    
    c=np.array([[x.ra/15.0, x.dec] for x in tiles])

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.scatter(c[:,0], c[:,1], s=5)
    plt.show()
            
if __name__ == "__main__":
    todo=True
