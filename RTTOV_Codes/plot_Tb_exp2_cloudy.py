import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.animation as animation
import matplotlib
import nclcmaps
import datetime
import cartopy.crs as ccrs
from cartopy import config
import cartopy.mpl.ticker as cticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.feature as cfeature
from matplotlib import gridspec

# Step 2
# Read Himawari OBS Data

PATH_OBS_Folder  = '/home/rohit/Rohit/NUS_Project/Volcano_ash/'
PATH_SAVE        = '/home/rohit/Rohit/NUS_Project/Volcano_ash/'
year     = 2015
monthdeb = 11
daydeb   = 4 
hourdeb  = 0
ncycles  = 1

latmin = -12
latmax = -5
lonmin = 112
lonmax = 120

base = datetime.datetime(year, monthdeb, daydeb, hourdeb)
timedelta = datetime.timedelta(1,0,0,0,0,0) #1day
datetime_series = [base+i*timedelta for i in range(ncycles)]

for date in datetime_series:
    file_OBS  = glob.glob(PATH_OBS_Folder + '%4.4i/%2.2i/%2.2i/*_02401.nc' %(date.year, date.month, date.day))
    for ifile_OBS in file_OBS:
        if '.nc' not in ifile_OBS:continue
        if 'rttov13' in ifile_OBS:continue
        print(ifile_OBS)

        PATH_OBS = PATH_OBS_Folder + '%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)
        ifile_OBS = ifile_OBS[-44:]
        
        FILE_OUT = ifile_OBS[-44:-3]+'.rttov13.'+'_AHI_Cloudy.NCDF'
                
        time = ifile_OBS[-28:-24]
        hh = time[0:2]
        mm = time[2:4]

        day = date.day
        month = date.month
        
        ncfile = Dataset (PATH_OBS + ifile_OBS)
        lat = (ncfile.variables['latitude'])[:]
        lon = (ncfile.variables['longitude'])[:]
        tb_11 = (ncfile.variables['tbb_11'])[:]
        tb_12 = (ncfile.variables['tbb_11'])[:]
        tb_13 = (ncfile.variables['tbb_13'])[:]
        tb_14 = (ncfile.variables['tbb_14'])[:]
        tb_15 = (ncfile.variables['tbb_15'])[:]
        tb_16 = (ncfile.variables['tbb_16'])[:]
        ncfile.close()
        
        lat    = lat [::5]
        lon    = lon [::5]
        tb_11  = tb_11[::5, ::5]
        tb_12  = tb_12[::5, ::5]
        tb_13  = tb_13[::5, ::5]
        tb_14  = tb_14[::5, ::5]
        tb_15  = tb_15[::5, ::5]
        tb_16  = tb_16[::5, ::5]
        
        ind_lat = np.where ((lat >=latmin) * (lat<=latmax))
        lat = lat [ind_lat]
        ind_lon = np.where ((lon >=lonmin) * (lon<=lonmax))
        lon = lon [ind_lon]
        
        tb_11   = tb_11 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
        tb_12   = tb_12 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
        tb_13   = tb_13 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
        tb_14   = tb_14 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
        tb_15   = tb_15 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]
        tb_16   = tb_16 [ind_lat[0][0]:ind_lat[0][-1]+1, ind_lon[0][0]:ind_lon[0][-1]+1]

        BTD = tb_13- tb_15
   
        BTD [BTD > 5]  =  np.nan
        BTD [BTD < -6] =  np.nan
     
        SWD = tb_14 - tb_15
                
        R = tb_15 - tb_13
        G = tb_13- tb_11
        B = tb_13

        ind = np.where((R<-4) | (R>2))
        R [ind] =np.nan

        ind = np.where((G < 0) | (G > 6))
        G[ind] = np.nan

        ind = np.where((B < 248) | (B > 303))
        B[ind] = np.nan
        
        rgb = np.dstack([R, G, B])
        rgb =np.flipud(rgb)
        rgb_masked = np.ma.masked_invalid(rgb)

        
        # RTTOV Output
        ncfile = Dataset (PATH_OBS + FILE_OUT)
        Tb_all = (ncfile.variables['Tb_all'])[:]
        ncfile.close()

        Tb_all   = Tb_all.reshape  (np.size(lat), np.size(lon), 6)

        BTD_sim = Tb_all [:,:,2]-Tb_all[:,:,4]

        BTD_sim [BTD_sim > 5]  =  np.nan
        BTD_sim [BTD_sim < -6] =  np.nan

        SWD_sim = Tb_all [:,:,3]-Tb_all[:,:,4]

        R_sim = Tb_all[:, :, 4] - Tb_all[:, :, 2]
        G_sim = Tb_all[:, :, 2] - Tb_all[:, :, 0]
        B_sim = Tb_all[:, :,2]

        ind = np.where((R_sim<-4) | (R_sim>2))
        R_sim [ind] =np.nan

        ind = np.where((G_sim < 0) | (G_sim > 6))
        G_sim[ind] = np.nan

        ind = np.where((B_sim < 248) | (B_sim > 303))
        B_sim[ind] = np.nan
        
        rgb_sim = np.dstack([R_sim, G_sim, B_sim])
        rgb_sim =np.flipud(rgb_sim)
        rgb_sim_masked = np.ma.masked_invalid(rgb_sim)

        # O-B
        diff_11 = tb_11 - Tb_all [:, :, 0]
        diff_12 = tb_12 - Tb_all [:, :, 1]
        diff_13 = tb_13 - Tb_all [:, :, 2]
        diff_14 = tb_14 - Tb_all [:, :, 3]
        diff_15 = tb_15 - Tb_all [:, :, 4]
        diff_16 = tb_16 - Tb_all [:, :, 5]
        
        lat0 = -8.4
        lon0 =  116.45

        vmin = 200
        vmax  = 315
        res_colorbar = 5
        clevs = np.arange(0, (vmax-vmin)/res_colorbar+1)*res_colorbar+vmin
        xv, yv = np.meshgrid(lon, lat)

        #cmap = plt.cm.viridis
        cmap = nclcmaps.cmap('amwg256')

        img_extent = (lonmin, lonmax, latmin, latmax)

        # Plot Tbs (obs vs Sim)
        nrows = 2
        ncols = 6
        fig, ax = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': ccrs.PlateCarree()},
                figsize=(15,4), gridspec_kw={'height_ratios': [1, 1], 'width_ratios':[1,1,1,1,1,1]})

        for irow in range(nrows):
            for icol in range (ncols):
                # subplot 1 (tb_11, OBS)
                if ((irow==0)*(icol==0)):
                    title = '8.6\u03bcm,OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_11, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 2 (tb_11, SIM)
                elif ((irow==1)*(icol==0)):
                    title = '8.6\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,0], cmap=cmap, vmin=vmin, vmax=vmax)   
                # subplot 3 (tb_12, OBS)
                elif ((irow==0)*(icol==1)):
                    title = '9.6\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_12, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 4 (tb_12, SIM)
                elif ((irow==1)*(icol==1)):
                    title = '9.6\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,1], cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 5 (tb_13, OBS)
                elif ((irow==0)*(icol==2)):
                    title = '10.4\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_13, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 6 (tb_13, SIM)
                elif ((irow==1)*(icol==2)):
                    title = '10.4\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,2], cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 7 (tb_14, OBS)
                elif ((irow==0)*(icol==3)):
                    title = '11.2\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_14, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 8 (tb_14, SIM)
                elif ((irow==1)*(icol==3)):
                    title = '11.2\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,3], cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 9 (tb_15, OBS)
                elif ((irow==0)*(icol==4)):
                    title = '12.4\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_15, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 10 (tb_15, SIM)
                elif ((irow==1)*(icol==4)):
                    title = '12.4\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,4], cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 11 (tb_16, OBS)
                elif ((irow==0)*(icol==5)):
                    title = '13.3\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv, tb_16, cmap=cmap, vmin=vmin, vmax=vmax)
                # subplot 12 (tb_16, SIM)
                elif ((irow==1)*(icol==5)):
                    title = '13.3\u03bcm, SIM'
                    cs= ax[irow][icol].pcolormesh(xv, yv, Tb_all[:,:,5], cmap=cmap, vmin=vmin, vmax=vmax)
        
                ax[irow][icol].set_title(title, fontsize= 8)
                ax[irow][icol].set_xticks(np.arange(lonmin, lonmax+3,3), crs=ccrs.PlateCarree())
                lon_formatter = cticker.LongitudeFormatter()
                ax[irow][icol].xaxis.set_major_formatter(lon_formatter)
                ax[irow][icol].set_yticks(np.arange(latmin,latmax+2, 2), crs=ccrs.PlateCarree())
                lat_formatter = cticker.LatitudeFormatter()
                ax[irow][icol].yaxis.set_major_formatter(lat_formatter)
                ax[irow][icol].coastlines()
                ax[irow][icol].plot(lon0, lat0,   marker='x', color='grey')    
                plt.colorbar(cs, ax= ax[irow][icol], orientation='vertical', shrink= 0.6)
        

        fig.subplots_adjust(right=0.9)
        plt.tight_layout()      
                
        plt.savefig(PATH_SAVE +'%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)+ 'Exp2_Cld_Tbs'+ str(day)+str(month)+str(year)+'_'+str(time)+'.png', dpi=300, bbox_inches = 'tight')
        plt.close()

        # histogram diff:


        fig = plt.figure()

        

        nrows = 2
        ncols = 3
        fig.set_figheight(8)
        fig.set_figwidth(15)
        spec = gridspec.GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1, 1], wspace=0.3, hspace= 0.3,
                         height_ratios=[1,1] )
        nbins=100

        for ifigure in range (6):
            ax = fig.add_subplot(spec[ifigure])
            ax.set_xlim(-10, 10)
            if ifigure ==0:
                title = 'O-B, 8.6\u03bcm, \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                ax.hist(diff_11.flatten(), bins=nbins)
            elif ifigure==1:
                title = 'O-B, 9.6\u03bcm'
                ax.hist(diff_12.flatten(), bins= nbins)
                ax.set_xlim(-5,30)
            elif ifigure==2:
                title = 'O-B, 10.4\u03bcm'
                ax.hist(diff_13.flatten(), bins = nbins)
            elif ifigure==3:
                title = 'O-B, 11.2\u03bcm'
                ax.hist(diff_14.flatten(), bins= nbins)
            elif ifigure==4:
                title = 'O-B, 12.4\u03bcm' 
                ax.hist(diff_15.flatten(), bins= nbins)
            else:
                title = 'O-B, 13.3\u03bcm'
                ax.hist(diff_16.flatten(), bins=nbins)

            ax.set_title (title)
            ax.set_xlabel('O-B [K]')
                         
                #ax[irow][icol].set_yticks(np.arange(latmin,latmax+2,2), crs=ccrs.PlateCarree())
        fig.tight_layout()           
                
        plt.savefig(PATH_SAVE +'%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)+ 'Exp2_Cld_hist'+ str(day)+str(month)+str(year)+'_'+str(time)+'.png', dpi=300, bbox_inches = 'tight')
        plt.close()





        

        vmin = -5
        vmax  = 5
        res_colorbar = 0.5
        clevs = np.arange(0, (vmax-vmin)/res_colorbar+1)*res_colorbar+vmin
        cmap = nclcmaps.cmap('GMT_relief_oceanonly')
        nrows = 2
        ncols = 3
        
        fig, ax = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': ccrs.PlateCarree()},
                figsize=(8,4), gridspec_kw={'height_ratios': [1, 1], 'width_ratios':[1,1,0.8]})

        for irow in range(nrows):
            for icol in range (ncols):
                # subplot 1 (OBS)
                if ((irow==0)*(icol==0)):
                    title = '10.4\u03bcm-12.4\u03bcm, OBS \n'+ str(day) +'/' + str(month)+'/'+ str(year)+','+\
                        str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv,BTD, cmap=cmap, vmin=vmin, vmax=vmax)
                    plt.colorbar(cs, ax= ax[irow][icol], orientation='vertical', shrink= 0.6)
                    
                elif ((irow==1)*(icol==0)):
                    title = '10.4\u03bcm-12.4\u03bcm, SIM '
                    cs= ax[irow][icol].pcolormesh(xv, yv,BTD_sim, cmap=cmap, vmin=vmin, vmax=vmax)
                    plt.colorbar(cs, ax= ax[irow][icol], orientation='vertical', shrink= 0.6)
                elif ((irow==0)*(icol==1)):
                    title = '11.3\u03bcm-12.4\u03bcm,OBS \n'+ str(day) + '/' +str(month)+ '/'+ str(year)+','+\
                str(hh)+'h'+ str(mm)+'min'
                    cs= ax[irow][icol].pcolormesh(xv, yv,SWD, cmap=cmap, vmin=vmin, vmax=vmax)
                    plt.colorbar(cs, ax= ax[irow][icol], orientation='vertical', shrink= 0.6)
                elif ((irow==1)*(icol==1)):
                    title = '11.3\u03bcm-12.4\u03bcm,SIM '
                    cs= ax[irow][icol].pcolormesh(xv, yv,SWD_sim, cmap=cmap, vmin=vmin, vmax=vmax)
                    plt.colorbar(cs, ax= ax[irow][icol], orientation='vertical', shrink= 0.6)

                elif ((irow==0)*(icol==2)):
                    title= 'RGB, OBS'
                    cs= ax[irow][icol].imshow(rgb_masked, aspect= 'equal',interpolation='none', origin='lower',
                            extent= img_extent, transform =ccrs.PlateCarree())
                    
                elif ((irow==1)*(icol==2)):
                    title= 'RGB, SIM'
                    cs= ax[irow][icol].imshow(rgb_sim_masked, aspect= 'equal',interpolation='none', origin='lower',
                         extent= img_extent, transform =ccrs.PlateCarree())
                    
                    
                ax[irow][icol].set_title(title, fontsize= 8)
                ax[irow][icol].set_xticks(np.arange(lonmin, lonmax+3,3), crs=ccrs.PlateCarree())
                lon_formatter = cticker.LongitudeFormatter()
                ax[irow][icol].xaxis.set_major_formatter(lon_formatter)
                ax[irow][icol].set_yticks(np.arange(latmin,latmax+2,2), crs=ccrs.PlateCarree())
                lat_formatter = cticker.LatitudeFormatter()
                ax[irow][icol].yaxis.set_major_formatter(lat_formatter)
                ax[irow][icol].coastlines()
                ax[irow][icol].plot(lon0, lat0,   marker='x', color='grey')
                
        fig.subplots_adjust(right=0.9)
        plt.tight_layout()      
                
        plt.savefig(PATH_SAVE +'%4.4i/%2.2i/%2.2i/' %(date.year, date.month, date.day)+ 'Exp2_Cld_RGB'+ str(day)+str(month)+str(year)+'_'+str(time)+'.png', dpi=300, bbox_inches = 'tight')
        plt.close()

        

                

                







