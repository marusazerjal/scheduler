"""
FunnelWeb: Manage list of observed tiles.
Talk to Jeeves.

There are two separate lists:
(1) Temporary list with numbered tile ids internal to Scheduler
(2) External list with TILE IDs

We need (1) because while the last tile is still being observed we already need to select the next one. So we assume that the last tile is going to be observed successfully and add tile_id of the last tile on the temporary list (1) of observed tiles.
If later Jeeves confirms that this tile was observed successfully, generate TILE ID and add it to the (2) External list of observed TILE IDs. 
If observation was UNsuccessfull, remove this tile from the temporary list. (But I need to know its internal TILE ID!!)
"""


"""
(1) From Scheduler: Receive internal TILE_ID_internal that was selected
(2) Write TILE_ID_internal to a file
(3) From Jeeves: Receive observation success info.
(4) Generate TILE_ID.
(5) Update file with TILE_IDs.
(6) If observation was NOT successfull, remove TILE_ID_internal from the temporary list.
"""


import numpy as np
import datetime
from collections import defaultdict

# FunnelWeb
import params

reload(params)


def news_from_Jeeves(success=None, tile_id=None):
    """
    If observation was successfull, assign external tile_id and add it on the list of observed tiles.
    Add this tile to pickle file of observed tiles.
    Do nothing to internal list because observation was successfull.
    """
    todo=True
    
    if success:
        todo=True
    
def load_internal_list_of_observed_tile_ids_and_divide_into_mag_ranges(tiles):
    """
    tiles: tiling code pickle file. Need this for mag ranges.
    """
    td={x.field_id: x for x in tiles}
    
    filename=params.params['observed_tiles_internal_filename']

    try:
        data=np.loadtxt(filename, ndmin=1, comments='#')
        data=[int(x) for x in data]

        d = defaultdict(list)

        # Split into magnitude ranges
        for t in data:
            x=td[t]
            d[(float(x.mag_min), float(x.mag_max))].append(x)
        d2={k: v for k, v in d.iteritems()} # because obs_tile_coo is defaultdict (creates new list in the dict if a key that does not exist yet is called)
        
        return d2

    except:
        return {}  
    
def load_internal_list_of_observed_tile_ids():
    """
    Load list of observed tile ids
    """
    filename=params.params['observed_tiles_internal_filename']

    try:
        data=np.loadtxt(filename, ndmin=1, comments='#')
        data=[int(x) for x in data]
        return data

    except:
        return []


# Internal (temporary) list of observed tiles
def add_tile_id_internal_to_the_list(tile_id_internal=None):
    filename=params.params['observed_tiles_internal_filename']

    try:
        data=np.loadtxt(filename, ndmin=1)
        data=[int(x) for x in data]
        data.append(tile_id_internal)
    except:
        data=[tile_id_internal]
        
    # TODO: later, when number of observed tiles becomes big, add only the last line. Do not overwrite file each time.
    f=open(filename, 'wb')
    for x in data:
        #~ f.write('%d %s %s \n'%(x[0], str(x[1]), str(x[2])))
        f.write('%d \n'%x)
    f.close()

def remove_tile_id_internal_from_the_list(tile_id_internal=None, filename=None):
    try:
        data=np.loadtxt(filename, ndmin=1)
        data=[int(x) for x in data if int(x)!=tile_id_internal] # TODO: 3 columns now
        
        f=open(filename, 'wb')
        for x in data:
            f.write(str(int(x))+'\n')
        f.close()
    except:
        # List of observed tiles was not created yet. Nothing to remove.
        pass


def assign_external_tile_id(best_tile=None): # TODO: should this be kept here or in file manager?
    tile_id_time = datetime.datetime.now()
    # TODO: time should be in UTC (more universal + there is no problems with daylight savings time)
    sign=0
    if best_tile.TaipanTile.dec>0.0:
        sign=1
    tile_id = '%03d%01d%02d%04d%02d%02d%02d%02d'%(best_tile.TaipanTile.ra, sign, np.abs(best_tile.TaipanTile.dec), tile_id_time.year, tile_id_time.month, tile_id_time.day, tile_id_time.hour, tile_id_time.minute)
    
    return int(tile_id) 

# External (eternal) list of observed tiles    
def add_tile_id_to_the_list(tile_id=None, filename=None):
    # We need to know some info about this tile. At least coordinates. Observing time as well.
    # We will probably need expusure time, weather info etc. as well.

   # Assign TILE ID
    tile_id=assign_external_tile_id(best_tile=None)


    try:
        data=np.loadtxt(filename, ndmin=1)
        data=[int(x) for x in data]
        data.append(tile_id)
    except:
        data=[tile_id]
        
    # TODO: later, when number of observed tiles becomes big, add only the last line. Do not overwrite file each time.
    f=open(filename, 'wb')
    for x in data:
        f.write(str(int(x))+'\n')
    f.close()


if __name__ == "__main__":
    manage_temporary_list_of_observed_tiles()
