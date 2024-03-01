import numpy as np
import netCDF4
from netCDF4 import Dataset
def read_era (PATH, file_profile, file_surface, hh, mm, latmin, latmax, lonmin, lonmax):

    # extract surface parameters
    ncfile = Dataset(PATH + file_surface)
    lat    = (ncfile.variables['latitude'])[:]
    lon    = (ncfile.variables['longitude'])[:]
    u10    = (ncfile.variables['u10'])[:]
    v10    = (ncfile.variables['v10'])[:]
    d2m    = (ncfile.variables['d2m'])[:]
    t2m    = (ncfile.variables['t2m'])[:]
    mdww   = (ncfile.variables['mdww'])[:]
    skt    = (ncfile.variables['skt'])[:]
    sp     = (ncfile.variables['sp'])[:] / 100  # convert into hPa
    z_surf = (ncfile.variables['z'])[:]
    lsm    = (ncfile.variables['lsm'])[:]
    ncfile.close()

    # extract atmospheric profiles
    ncfile = Dataset(PATH+ file_profile)
    T       = (ncfile.variables['t'])[:] # Temperature profile
    Q       = (ncfile.variables['q'])[:]  # specific humidity
    O3      = (ncfile.variables['o3'])[:] # Ozone mass mixing ratio
    time    = (ncfile.variables['time'])[:] # time
    cfrac   = (ncfile.variables['cc'])[:] # fraction of cloud cover
    z       = (ncfile.variables['z'])[:] # geopotential
    ciwc    = (ncfile.variables['ciwc'])[:] # cloud ice water
    clwc    = (ncfile.variables['clwc'])[:] # cloud liquid water
    crwc    = (ncfile.variables['crwc'])[:] # cloud rain water
    cswc    = (ncfile.variables['cswc'])[:] # cloud snow water
    ncfile.close()

    # sptaial collocation:
    ind_lat = np.where ((lat >= latmin) * (lat<= latmax))
    lat_era = lat [ind_lat]

    ind_lon = np.where ((lon >= lonmin) * (lon<= lonmax))
    lon_era = lon [ind_lon]

    # collocated_surface:
    u10  =u10  [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    v10  =v10  [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    d2m  =d2m  [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    t2m  =t2m  [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    mdww =mdww [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    skt  =skt  [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    sp   =sp   [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    z_surf =z_surf [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    lsm  =lsm   [:, ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]

    # collocated_profiles
    T = T [:, :,  ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    Q = Q [:, :,  ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    O3= O3[:, :,  ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    cfrac= cfrac [:, :,  ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    z= z [:, :,  ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
    ciwc = ciwc [:, :, ind_lat[0][0]:ind_lat[0][-1] + 1, ind_lon[0][0]:ind_lon[0][-1] + 1]
    clwc = clwc [:, :, ind_lat[0][0]:ind_lat[0][-1] + 1, ind_lon[0][0]:ind_lon[0][-1] + 1]
    crwc = crwc[:, :, ind_lat[0][0]:ind_lat[0][-1] + 1, ind_lon[0][0]:ind_lon[0][-1] + 1]
    cswc = cswc[:, :, ind_lat[0][0]:ind_lat[0][-1] + 1, ind_lon[0][0]:ind_lon[0][-1] + 1]

    # time index:
    # extract hour from ERA:
    #dates = netCDF4.num2date(time, units='hours since 1900-01-01 00:00:00')
    # hh could be time index
    nlevels = 37
    time    = 24
    nprofiles = lat_era.shape[0] * lon_era.shape[0]
    T = T.reshape (time, nlevels, nprofiles)
    Q = Q.reshape (time, nlevels, nprofiles)
    O3 = O3.reshape (time, nlevels, nprofiles)
    cfrac = cfrac.reshape (time, nlevels, nprofiles)
    z = z.reshape (time, nlevels, nprofiles)
    ciwc = ciwc.reshape (time, nlevels, nprofiles)
    clwc = clwc.reshape (time, nlevels, nprofiles)
    crwc = crwc.reshape (time, nlevels, nprofiles)
    cswc = cswc.reshape (time, nlevels, nprofiles)

    u10 = u10.reshape (time, nprofiles)
    v10 = v10.reshape (time, nprofiles)
    d2m = d2m.reshape (time, nprofiles)
    t2m = t2m.reshape (time, nprofiles)
    mdww = mdww.reshape (time, nprofiles)
    skt  = skt.reshape (time, nprofiles)
    sp   = sp.reshape (time, nprofiles)
    z_surf = z_surf.reshape (time, nprofiles)
    lsm    = lsm.reshape (time, nprofiles)

    u10 = u10 [hh, :]
    v10 = v10 [hh, :]
    d2m = d2m [hh, :]
    t2m = t2m [hh, :]
    mdww = mdww [hh, :]
    skt  = skt [hh, :]
    sp  = sp [hh, :]
    z_surf  = z_surf [hh, :]
    lsm  = lsm [hh, :]

    lsm [lsm <= 0.5]= 0.5     # indicate ocean
    lsm [lsm > 0.5] = 0.1     # indicate land
    
    lsm [lsm==0.5] = 1       # in according to RTTOV structure
    lsm [lsm==0.1] = 0.0     # in accordance with RTTOV structure

    # subset for obs:
    T  = np.transpose (T[hh, :, :])
    Q  = np.transpose (Q[hh, :, :])
    O3  = np.transpose (O3[hh, :, :])
    cfrac  = np.transpose (cfrac[hh, :, :])
    z      = np.transpose (z[hh, :, :])
    ciwc   = np.transpose (ciwc[hh, :, :])
    clwc   = np.transpose (clwc[hh, :, :])
    crwc   = np.transpose (crwc[hh, :, :])
    cswc   = np.transpose (cswc[hh, :, :])

    return lat_era, lon_era, u10, v10, d2m, t2m, mdww, skt, sp, z_surf, \
    lsm, T, Q, O3, cfrac, z, ciwc, clwc, crwc, cswc, nprofiles, nlevels

    


