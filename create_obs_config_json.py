"""
Create FunnelWeb ObsConfig.json file for a selected tile.
"""
import os
import errno
import numpy as np
import json

from astropy import units as u
from astropy.coordinates import SkyCoord

import params

reload(params)

def create_ObsConfig_json(tile=None, utc=None):
    tile=tile.TaipanTile
    
    ra_tile, dec_tile = format_coordinates(tile.ra, tile.dec)
    
    fibres={v: k for k, v in tile.fibres.iteritems()}
    
    # targets: including standards?
    targets=[]
    for x in tile.get_assigned_targets_science():
        ra, dec = format_coordinates(x.ra, x.dec)
        targets.append({'sbID': fibres[x], 'ra': ra, "dec": dec, "xMicrons": -50400.0, "yMicrons": -54500.0, "targetID": x.idn, "mag": x.mag})

    guides=[]
    for x in tile.get_assigned_targets_guide():
        ra, dec = format_coordinates(x.ra, x.dec)
        guides.append({'sbID': fibres[x], 'ra': ra, "dec": dec, "xMicrons": -50400.0, "yMicrons": -54500.0, "targetID": x.idn, "mag": x.mag})
    
    # What does sky need 'mag'?
    sky=[{'sbID': fibres[x], 'ra': "02:41:04.8", "dec": "-08:15:21", "xMicrons": -50400.0, "yMicrons": -54500.0, "targetID": x.idn, "mag": x.mag} for x in tile.get_assigned_targets_guide()]
    
    # is utc time right now or when this tile is supposed to be observed?
    utc=str(utc)
    utc=utc.replace(' ', 'T').split('.')[0]
    
    #"fieldID": For now it is internal fieldID. It is unique only within the current tiling. Should I change it?
    
    data={
        "configFormatVersion": 0.1,
        "instrumentName": "AAO.TAIPAN",
        "filePurpose": "Instrument configuration",
        "notes": "Any notes that the tiler and or router want to add. We may or may not actually need this, but it seemed a good idea to include a field for any explanatory text that we may wish to include.",
        "origin":
        [
        {"name": "FunnelWebSurvey.Tiler", "software": "tiler_executable", "version": "v2.34_20170109", "execDate": "2017-01-23T16:33:58+11:00"},
        {"name": "TAIPAN.Router", "software": "findRoute", "version": "v1.23_20170101", "execDate": "2017-01-24T16:33:58+11:00"}
        ],

        "configurationPhase": "tiling",
        "survey": "FunnelWeb",
        "fieldID": tile.field_id,
        "fieldCentre":
        {
        "ra": ra_tile,
        "dec": dec_tile,
        "UT": utc
        },
        
        "targets": targets,
        "guideStars": guides,
        "sky": sky,

        "additionalFITSkeywords":
        {
        "KEYWORD1": "These FITS header items are passed straight to the instrument software...",
        "ANOTHER": "...which writes the FITS file. The instrument doesn't act on them in any other way.",
        "YET_MORE": "someValue"
        },
        
        "routable": 'false',
        "starbugRoute":
        [
        {
            "tickID": 0,
            "tickType": "initialPositions",
            "data":
            [
            {"sbID": 1, "xMicrons": -50400.0, "yMicrons": -74500.0, "thetaDeg": 0.0},
            {"sbID": 2, "xMicrons": -26400.0, "yMicrons": -98500.0, "thetaDeg": 10.0},
            {"sbID": 3, "xMicrons": -54800.0, "yMicrons": -90900.0, "thetaDeg": 20.0}
            ]
        },
        {
            "tickID": 1, 
            "tickType": "rotation",
            "data":
            [
            {"sbID": 1, "thetaDeg": 0.0, "deltaTheta": 0.0},
            {"sbID": 2, "thetaDeg": 0.0, "deltaTheta": -10.0},
            {"sbID": 3, "thetaDeg": 20.0, "deltaTheta": 0.0}
            ]
        },
        {
            "tickID": 2,
            "tickType": "translation",
            "data":
            [
            {"sbID": 1, "xMicrons": -50400.0, "yMicrons": -54500.0, "deltaX": 0.0, "deltaY": 20000.0},
            {"sbID": 2, "xMicrons": -26400.0, "yMicrons": -78500.0, "deltaX": 0.0, "deltaY": 20000.0},
            {"sbID": 3, "xMicrons": -54800.0, "yMicrons": 0.0, "deltaX": 0.0, "deltaY": 0.0},
            {"sbID": 4, "xMicrons": 10000.0, "yMicrons": 0.0, "deltaX": 0.0, "deltaY": 0.0}
            ]
        }
        ]
    }

    t=utc.split('T')[1].replace(':', '')
    date=utc.split('T')[0].replace('-', '')
    
    folder = params.params['obs_config_json_folder'] + date + '/'
    try:
        os.makedirs(folder)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    
    filename = folder + '%d_%s.obs_config.json'%(tile.field_id, t)
    with open(filename, 'w') as outfile:  
        json.dump(data, outfile, indent=4, sort_keys=True)    
    print '%s created.'%filename
    
    return filename

def format_coordinates(ra, dec):
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    ra1='%02d:%02d:%02.1f'%(c.ra.hms.h, c.ra.hms.m, c.ra.hms.s)
    dec1='%02d:%02d:%02.0f'%(c.dec.dms.d, np.abs(c.dec.dms.m), np.abs(c.dec.dms.s))
    # TODO: round coordinates??
    return ra1, dec1   
