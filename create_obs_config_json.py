"""
Create FunnelWeb ObsConfig.json file for a selected tile.
"""
import os
import errno
import numpy as np
import simplejson as json

from astropy import units as u
from astropy.coordinates import SkyCoord

import random # just for testing. delete!

#~ from json import encoder
#~ encoder.FLOAT_REPR = lambda o: format(o, '.4f')
from decimal import Decimal

import params

reload(params)

#~ class DecimalEncoder(json.JSONEncoder):
    #~ def _iterencode(self, o, markers=None):
        #~ if isinstance(o, decimal.Decimal):
            #~ # wanted a simple yield str(o) in the next line,
            #~ # but that would mean a yield on the line with super(...),
            #~ # which wouldn't work (see my comment below), so...
            #~ return (str(o) for o in [o])
        #~ return super(DecimalEncoder, self)._iterencode(o, markers)

def create_ObsConfig_json(tile=None, utc=None):
    tile=tile.TaipanTile
    
    fibres={v: k for k, v in tile.fibres.iteritems()}

    # targets: including standards?
    science = [{"type": "science", 'bugLemoID': fibres[x], 'ra': x.ra, "dec": x.dec, "pmRA": 50.0, "pmDec": -10.0, "mag": Decimal(str(x.mag)), "xMicrons": -50400.0, "yMicrons": -54500.0, "targetID": str(x.idn)} for x in tile.get_assigned_targets_science()] # FIXME
    guides = [{"type": "guide", 'bugLemoID': fibres[x], 'ra': x.ra, "dec": x.dec, "pmRA": 50.0, "pmDec": -10.0, "xMicrons": -50400.0, "yMicrons": -54500.0, "targetID": str(x.idn), "mag": Decimal(str(x.mag))} for x in tile.get_assigned_targets_guide()] # FIXME
    
    
    # not perfect, but it is good enough
    n=len(tile.get_assigned_targets_guide())
    PA = [random.uniform(-np.pi, np.pi) for i in range(n)]
    R = [random.uniform(0, 2.9) for i in range(n)]
    
    #~ sky = [{"type": "sky", 'bugLemoID': fibres[x], 'ra': tile.ra+random.uniform(-2.7, 2.7), "dec": tile.dec+random.uniform(-2.7, 2.7), "xMicrons": -50400.0, "yMicrons": -54500.0} for x in tile.get_assigned_targets_guide()] # FIXME
    sky = [{"type": "sky", 'bugLemoID': fibres[x], 'ra': tile.ra+r*np.sin(pa), "dec": tile.dec+r*np.cos(pa), "xMicrons": -50400.0, "yMicrons": -54500.0} for x, r, pa in zip(tile.get_assigned_targets_guide(), R, PA)] # FIXME

    objects = science + guides + sky
    
    # is utc time right now or when this tile is supposed to be observed?
    utc=str(utc)
    utc=utc.replace(' ', 'T').split('.')[0]
    
    #"fieldID": For now it is internal fieldID. It is unique only within the current tiling. Should I change it?

    data={
        "schemaID": 1,
        "instrumentName": "AAO.Taipan", 
        "origin": [
        {
           "name": "FunnelWebSurvey.Tiler",
           "software": "tiler_executable_name",
           "version": "v2.34_20170109", # FIXME
           "execDate": "2017-01-23T16:33:58+11:00" # FIXME
        }
        ], 
        "routable": "unknown", 
        "fieldID": tile.field_id,
        "filePurpose": "CoD.Testing", # FIXME
        "configFormatVersion": 0.4, # FIXME
        "fieldCentre": {
        "ra": tile.ra,
        "dec": tile.dec
        }, 

        "telModel":{ # FIXME
          "cenWave": 6000.0,
          "obsWave": 5000.0,
          "UT": "2018-10-22 12:00:00"
        },

        "targets": objects

        }

    t=utc.split('T')[1].replace(':', '')
    date=utc.split('T')[0].replace('-', '')
    
    #~ if simulation_nickname is not None:
        #~ folder = params.params['obs_config_json_folder'].replace('funnelweb', 'funnelweb_simulation_%s'%simulation_nickname) + date + '/'
    #~ else:
    folder = params.params['obs_config_json_folder'] + date + '/'
    
    try:
        os.makedirs(folder)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #~ print 'Printing json data'
    #~ print data
    
    filename = folder + '%d_%s.obs_config.json'%(tile.field_id, t)
    with open(filename, 'w') as outfile:  
        #~ print data
        json.dump(data, outfile, indent=4, sort_keys=True)    # , cls=DecimalEncoder
    print '%s created.'%filename
    
    return filename
