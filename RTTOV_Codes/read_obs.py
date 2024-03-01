import numpy as np
import netCDF4
from netCDF4 import Dataset

def read_obs(PATH, filename, latmin, latmax, lonmin, lonmax):
    ncfile = Dataset(PATH + filename)
    lat    = (ncfile.variables['latitude'])[:]
    lon    = (ncfile.variables['longitude'])[:]
    SAZ    = (ncfile.variables['SAZ'])[:]
    SAA    = (ncfile.variables['SAA'])[:]
    SOZ    = (ncfile.variables['SOZ'])[:]
    SOA    = (ncfile.variables['SOA'])[:]
    tb_11  = (ncfile.variables['tbb_11'])[:]
    tb_12  = (ncfile.variables['tbb_12'])[:]
    tb_13  = (ncfile.variables['tbb_13'])[:]
    tb_14  = (ncfile.variables['tbb_14'])[:]
    tb_15  = (ncfile.variables['tbb_15'])[:]
    tb_16  = (ncfile.variables['tbb_16'])[:]
    ncfile.close()
    print(filename)

    time = filename[-28:-24]
    hh= time [0:2]
    mm =time [2:4]

    hh = int(hh)
    mm = int(mm)

    ind_lat = np.where ((lat >= latmin) * (lat<=latmax))
    lat     = lat [ind_lat]
    ind_lon = np.where ((lon >=lonmin)  * (lon<=lonmax))
    lon     = lon [ind_lon]
    
    SAZ     = SAZ   [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    SOZ     = SOZ   [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    SOA     = SOA   [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    SAA     = SAA   [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_11   = tb_11 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_12   = tb_12 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_13   = tb_13 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_14   = tb_14 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_15   = tb_15 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    tb_16   = tb_16 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]

    
    lat    = lat [::5]
    lon    = lon [::5]

    SAZ    = SAZ [::5, ::5]
    SAA    = SAA [::5, ::5]
    SOZ    = SOZ [::5, ::5]
    SOA    = SOA [::5, ::5]
    tb_11  = tb_11 [::5, ::5]
    tb_12  = tb_12 [::5, ::5]
    tb_13  = tb_13 [::5, ::5]
    tb_14  = tb_14 [::5, ::5]
    tb_15  = tb_15 [::5, ::5]
    tb_16  = tb_16 [::5, ::5]

    nprofiles = lat.shape[0] * lon.shape[0]
    SAZ       = SAZ.reshape(nprofiles)
    SAA       = SAA.reshape(nprofiles)
    SOZ       = SOZ.reshape(nprofiles)
    SOA       = SOA.reshape(nprofiles)
    tb_11     = tb_11.reshape(nprofiles)
    tb_12     = tb_12.reshape(nprofiles)
    tb_13     = tb_13.reshape(nprofiles)
    tb_14     = tb_14.reshape(nprofiles)
    tb_15     = tb_15.reshape(nprofiles)
    tb_16     = tb_16.reshape(nprofiles)

    return lat, lon, SAZ, SAA, SOZ, SOA, tb_11, tb_12, tb_13,\
            tb_14, tb_15, tb_16, hh, mm



    


