




def locs_from_dataless(dataless_path):

    from obspy.xseed import Parser
    import numpy as np
    sp = Parser(dataless_path)
    info = {}
    
    metadata = Parser.getInventory(sp)

#    print(metadata['channels'][0])

    lats = np.asarray([float(i['latitude']) for i in metadata['channels']])
    lons = np.asarray([float(i['longitude']) for i in metadata['channels']])
    stations = np.asarray([i['channel_id'].split('..')[0] for i in metadata['channels']])
    net_stats = np.asarray(['{}.{}'.format(i.split('.')[0], i.split('.')[1]) for i in stations])
    #lats = np.unique(lats)    
    #lons = np.unique(lons)         
    #stations = np.unique(stations)
    
    for i, stat in enumerate(net_stats):    
        info[stat] = [lons[i], lats[i]]
        
    return info

#dataless_path = 'UOM.dataless'
#info = locs_from_dataless(dataless_path)

#print(info); print(len(info))
