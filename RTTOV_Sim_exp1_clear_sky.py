import numpy as np
import read_era
import read_obs
import pyrttov
import scipy.constants
import sys
g = scipy.constants.g
from netCDF4 import Dataset
import datetime
import glob

# Step 1
# Read observation
PATH_OBS_Folder  = '/home/rohit/Rohit/NUS_Project/Volcano_ash/'
PATH_ERA_Folder  = '/home/rohit/Rohit/NUS_Project/Volcano_ash/ERA/'
PATH_SAVE_Folder  = '/home/rohit/Rohit/NUS_Project/Volcano_ash/'

year     = 2015
monthdeb = 11
daydeb   = 4 
hourdeb  = 0
ncycles  = 1

# Assign domain for simulation
latmin = -12
latmax = -5
lonmin = 112
lonmax= 120

base = datetime.datetime(year, monthdeb, daydeb, hourdeb)
timedelta = datetime.timedelta(1,0,0,0,0,0) #1day
datetime_series = [base+i*timedelta for i in range(ncycles)]

file_surface = 'era_surface.nc'
file_profile = 'era_profiles.nc'

for date in datetime_series:
    file_OBS  = glob.glob(PATH_OBS_Folder+'%4.4i/%2.2i/%2.2i/*_02401.nc' %(date.year, date.month, date.day))
    for ifile_OBS in file_OBS:
        if '.nc' not in ifile_OBS:continue
        if 'rttov13' in ifile_OBS:continue
        print(ifile_OBS)

        PATH_OBS = PATH_OBS_Folder + '%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)
        file_OBS = ifile_OBS[-44:]
        print(ifile_OBS)
        FILE_OUT = ifile_OBS[-44:]+'.rttov13.'+'_AHI_Clear.NCDF'
        print(FILE_OUT)
             
        # reading observations
        lat_obs, lon_obs, SAZ, SAA, SOZ, SOA, tb_11, tb_12, tb_13,\
        tb_14, tb_15, tb_16, hh,mm = read_obs.read_obs (PATH_OBS, file_OBS, latmin, latmax, lonmin, lonmax)
 
        # reading ERA reanalysis data and collocated to observatio lat-lon grid
        PATH_ERA   = PATH_ERA_Folder + '%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day) 
        lat_era, lon_era, u10, v10, d2m, t2m, mdww, skt, sp, z_surf, \
            lsm, T, Q, O3, cfrac, z, ciwc, clwc, crwc, cswc, nprofiles, nlevels = read_era.read_era(PATH_ERA,file_profile,\
                                                                                                    file_surface, hh, mm,\
                                                                                                    latmin, latmax, lonmin,lonmax)
        # create an lat-lon matrix
        lon, lat = np.meshgrid(lon_era, lat_era)
        lon = lon.reshape(nprofiles)
        lat = lat.reshape(nprofiles)

        # Computing surface humidity at 2m 
        q2m = np.zeros(nprofiles, dtype= np.float64)
        for iprofile in range (nprofiles):
            e= 6.1078 * np.exp(17.2693882 * (d2m [iprofile] - 273.16)/(d2m[iprofile]-35.86))
            q2m[iprofile] = 0.622 * e /(sp[iprofile] - 0.378 *e)

        #Arrays for pyrttov (RTTOV structure)
        dates         = np.zeros([nprofiles,6],          dtype=np.float64)
        s2m           = np.zeros([nprofiles,6],          dtype=np.float64)
        skin          = np.zeros([nprofiles,9],          dtype=np.float64)
        surfgeom      = np.zeros([nprofiles,3],          dtype=np.float64)
        surftype      = np.zeros([nprofiles, 2],         dtype=np.float64)
        angles        = np.zeros([nprofiles,4],          dtype=np.float64)
        icecloud      = np.zeros([nprofiles,2],          dtype=np.float64)

        # Assigning inputs to RTTOV arrays
        s2m[:,0] = sp[:]
        s2m[:,1] = t2m[:]
        s2m[:,2] = q2m[:]
        s2m[:,3] = u10[:]
        s2m[:,4] = v10[:]
        s2m[:,5] = 100000


        skin[:,0] = skt[:]
        skin[:,1] = 35.
        skin[:,2] = 0.
        skin[:,3] = 0.
        skin[:,4:] = [3.,5.,15.,0.1,0.3]

        surfgeom[:,0] = lat[:]
        surfgeom[:,1] = lon[:]
        surfgeom[:,2] = z_surf[:]/(g*1000)

        surftype[:, 0] = lsm[:]

        angles[:, 0] =SAZ   # satellite zenith angle
        angles[:, 1] =SAA   # satellite azimuth angle
        angles[:, 2] =SOZ   # solar zenith angle
        angles[:, 3] =SOA   # Solar Azimuth angle

        # ERA Pressure levels
        P = np.array([1, 2, 3, 5,7,10,20,30,50,70,100,125,150, 175,200,225,250,300,350,\
                      400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000])

        # Construct pressure level matrix for n profiles:
        P_profiles = np.zeros(( nprofiles, nlevels), dtype=np.float64)
        for ilevel in range(nlevels):
            P_profiles[:,ilevel] = P[ilevel]

        # defined an RTTOV profile structrue:
        myCloudProfiles = pyrttov.Profiles(nprofiles, nlevels)
        myCloudProfiles.DateTimes = dates

        myCloudProfiles.GasUnits = 1  # kg_per_kg
        myCloudProfiles.P        = P_profiles
        myCloudProfiles.T        = T
        myCloudProfiles.Angles   = angles
        myCloudProfiles.S2m      = s2m
        myCloudProfiles.Skin     = skin
        myCloudProfiles.SurfType = surftype
        myCloudProfiles.SurfGeom = surfgeom
        


        myCloudProfiles.Gases = np.zeros((2, nprofiles, nlevels), dtype=np.float32)
        myCloudProfiles.Gases[0, :, :]= Q            # humidity
        myCloudProfiles.Gases[1, :, :]= O3            # Ozone
        myCloudProfiles.GasId = [1, 2]

        # Run RTTOV

        AHIRttov = pyrttov.Rttov()

        AHIRttov.Profiles = myCloudProfiles
        nchan_ahi = 6
        chan_list_ahi = (11, 12, 13, 14, 15, 16)
        AHIRttov.FileCoef='/home/rohit/Rohit/rttov/rttov132/rtcoef_rttov13/rttov13pred54L/rtcoef_himawari_8_ahi_o3.dat'
        AHIRttov.Options.OzoneData = True               # Ozone data flag is True
        AHIRttov.Options.AddInterp = True               # Use the RTTOV interpolator
        AHIRttov.Options.AddAerosl = False              # Activate aerosol simulations
        AHIRttov.Options.AddClouds = False              # Activate cloud simulations
        AHIRttov.Options.VerboseWrapper = True          # turn on verbose wrapper output
        AHIRttov.Options.Opdep13GasClip= False           
        AHIRttov.Options.IrScattModel = False           # Use chou-scaling for thermal radiated radiation
        AHIRttov.Options.UserAerOptParam = False        # Use explicit cloud optical properties
        AHIRttov.Options.StoreRad = True                # Store all radiance outputs
        AHIRttov.Options.AddSolar = False               # Enable solar calcuations for solar affected channels
        AHIRttov.Options.RayleighSingleScatt = False    # the Rayleigh single-scattering calculation for visible channels
        AHIRttov.Options.IrSeaEmisModel =2              # IR sea-surface emissivity model (IREMIS)
              


        ## Load the instrument (reads all channels)
        try:
            AHIRttov.loadInst(chan_list_ahi)
        except pyrttov.RttovError as e:
            sys.stderr.write("Error loading instrument(s): {!s}".format(e))
            sys.exit(1)

        
        print ('Run RTTOV')
        try:
            AHIRttov.runDirect()
    
        except pyrttov.RttovError as e:
            sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
            sys.exit(1)

        # All-sky Tb (Clear-sky + cloudy sky)
        Tb_all       = AHIRttov.Bt
        Tb_clear     = AHIRttov.BtClear

        PATH_SAVE = PATH_SAVE_Folder + '%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)
        
        ncfile = Dataset(PATH_SAVE + FILE_OUT,'w')
        ncfile.createDimension('nprofiles',size=nprofiles)                                                         
        ncfile.createDimension('nchannels',size=6)
       
        datanc          = ncfile.createVariable('Tb_all',np.float64,('nprofiles','nchannels'))
        datanc[:,:]     = np.float64(Tb_all)

        datanc          = ncfile.createVariable('Tb_clear',np.float64,('nprofiles','nchannels'))
        datanc[:,:]     = np.float64(Tb_clear)

        ncfile.close()
       


