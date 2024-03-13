#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:16:08 2024

@author: arturo
"""
#https://polar.ncep.noaa.gov/waves/examples/usingpython.shtml
#http://www.himpactwxlab.com/home/how-to-wiki/write-grib2-data-with-pygrib
#https://jswhit.github.io/pygrib/docs/
#https://pyhogs.github.io/intro_netcdf4.html
import numpy as np
import pygrib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans
import netCDF4 as nc4
from datetime import datetime,timedelta
from netCDF4 import Dataset
from numpy import load,save,loadtxt,savetxt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import math
from scipy.interpolate import Rbf
from scipy.io import loadmat
from matplotlib import path
import csv

#Paleta de Colores 
colo=np.array([
[255, 255, 255],
[255, 170, 0  ],
[163, 255, 114],
[112, 168, 0  ],
[0  , 169, 230],
[0  , 77 , 168],
[169, 0  , 230]])/255

colo1=ListedColormap(colo)
anomcolo=np.array([
[168, 112, 0  ],
[255, 255, 114],
[255, 255, 255],
[209, 255, 114],
[ 76, 230,   0],
[0  , 168 , 132]])/255
colo1=ListedColormap(colo)
anomcolo1=ListedColormap(anomcolo)
res=100
Nr=1
Nc=1

#ENTRADA DE DATOS  MAX-MIN
infiles1 = list(open('/home/arturo/Desktop/era5/infilesTmax', 'r'))
infiles2 = list(open('/home/arturo/Desktop/era5/infilesTmin', 'r'))
length=np.shape(infiles1)

#Ciclo de rango
for ff in range(0,length[0]):
 print(infiles1[ff][:-1])
 f1=infiles1[ff][:-1]
 f2=infiles2[ff][:-1] 
 nc1=Dataset(f1)
 nc2=Dataset(f2)

#Declaracion de Variables 
 lon1 = nc1.variables['longitude'][:]
 lat1 = nc1.variables['latitude'][:]
 loon = nc1.variables['longitude'].size
 laat = nc1.variables['latitude'].size
 lon,lat = np.meshgrid(lon1,lat1)
 tmaxd=nc1.variables['t2m'][:]
 tmind=nc2.variables['t2m'][:]
#Umbral
 Tu=39;
 Tl=5;
 ll=np.shape(tmaxd)
 UC_1 = np.empty([ll[0],ll[1],ll[2]])
#ciclo de formula para unidades calor. "no se toca" 
 for d in range(0, ll[0]):
   alpha = (tmaxd[d,:,:]-tmind[d,:,:])/2
   theta2=np.arcsin((Tu-((tmaxd[d,:,:]+tmind[d,:,:])/2))/alpha[:,:])
   theta1=np.arcsin((Tl-((tmaxd[d,:,:]+tmind[d,:,:])/2))/alpha[:,:])
   for i in range(0, ll[1]):
     for j in range(0,ll[2]):
       if (tmaxd[d,i,j]>Tu and tmind[d,i,j]>Tu): #Caso 1
         UC_1[d,i,j]=Tu-Tl
       elif (tmaxd[d,i,j]<Tl and tmind[d,i,j]<Tl): #Caso 2
         UC_1[d,i,j]=0
       elif (tmaxd[d,i,j]<Tu and tmind[d,i,j]>Tl): #Caso 3
         UC_1[d,i,j]=6*(tmaxd[d,i,j]+tmind[d,i,j]-2*Tl)/12;
       elif (tmaxd[d,i,j]<Tu and tmind[d,i,j]<Tl): #Caso 4
#     disp('Caso 4. Interceptado por el umbral más bajo.')
         UC_1[d,i,j]=(1/np.pi)*(((tmaxd[d,i,j]+tmind[d,i,j])/2-Tl)*(np.pi/2-theta1[i,j])+(alpha[i,j]*np.cos(theta1[i,j])))
       elif (tmaxd[d,i,j]>Tu and tmind[d,i,j]>Tl): #Caso 5
#     disp('Caso 5')
         UC_1[d,i,j]=1/np.pi*((((tmaxd[d,i,j]+tmind[d,i,j])/2)-Tl)*(theta2[i,j]+np.pi/2)+((Tu-Tl)*(np.pi/2-theta2[i,j]))-(alpha[i,j]*np.cos(theta2[i,j])))
       elif (tmaxd[d,i,j]>Tu and tmind[d,i,j]<Tl): #Caso6
#     disp('Caso 6')
         UC_1[d,i,j]=1/np.pi*((((tmaxd[d,i,j]+tmind[d,i,j])/2)-Tl)*(theta2[i,j]-theta1[i,j])+(alpha[i,j]*(np.cos(theta1[i,j])-np.cos(theta2[i,j])))+((Tu-Tl)*(np.pi/2-theta2[i,j])))        
       else:
         UC_1[d,i,j]=0#np.nan #0 #Aqui deja en nan para identificar donde si son valores diferentes a 0
UC = np.sum(UC_1,axis=0)
#Trazar la inofmracionn en un mapa 
# fig,ax2 =plt.subplots(Nr,Nc,figsize=(12,12),sharex=True) #,constrained_layout=True
# ax=plt.subplot(Nr,Nc,1)
# m = Basemap(projection='merc',llcrnrlon=-118,urcrnrlon=-86,llcrnrlat=14,urcrnrlat=33,resolution='l')
# m.readshapefile('/home/arturo/Documents/gadm36_MEX_shp/gadm36_MEX_1',name='FIRST_ESTA',drawbounds=True,color='k',linewidth=1)

# x, y = m(lon,lat)
# parallels=np.arange(10,41,10.)
# clevs2=[0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250] 
# cmap2=plt.cm.get_cmap('YlOrRd',50) 
# norm = BoundaryNorm(clevs2, ncolors=cmap2.N, clip=True)

# UCmask=maskoceans(lon,lat,UC)

# cs = m.contourf(x,y,UCmask,clevs2,cmap=cmap2, extend='max',norm=norm)#
# m.drawparallels(np.arange(-90,90,5),labels=[True,False,True,False],fontsize=16,linewidth=0.0001) #para escribir latitudes
# m.drawmeridians(np.arange(-180,180,10),labels=[1,0,0,1],latmax=33,fontsize=16,linewidth=0.0001)#para escribir longitudes


# plt.title('INIFAP. FV3GFS (0.25°)',fontsize=18)
# cb_ax = fig.add_axes([.92, 0.27, 0.02, 0.46])
# cbar = fig.colorbar(cs, cax=cb_ax,pad="12%",drawedges=False, ticks=[0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
# cbar.ax.tick_params(labelsize=16)

# fig.savefig("./test.png", format='png', bbox_inches='tight',dpi=res)

#parte de aqui modificasmos mañana 23 de febrero 2024Enero_1964.nc 
les =np.empty((laat*loon,4))
G2=np.reshape(UC,(1,laat*loon))
E2=np.reshape(lon,(1,laat*loon))
C2=np.reshape(lat,(1,laat*loon))


# ###############################
#Save Document SCV
for f in range(0,laat*loon):
  les[f,:] = list([f+1,E2[0,f],C2[0,f],G2[0,f]])  
savetxt('../../Desktop/era5/'+str(infiles1[0][35:38])+'_'+str(infiles1[0][26:30])+'.csv',les, delimiter=',')

with open('../../Desktop/era5/'+str(infiles1[0][35:38])+'_'+str(infiles1[0][26:30])+'.csv','r') as csvread:
    rd=csv.reader(csvread)
    lines=list(rd)
    lines.insert(0,['punto','lon','lat','t2m'])
    
with open('../../Desktop/era5/'+str(infiles1[0][35:38])+'_'+str(infiles1[0][26:30])+'.csv','w',newline='') as csvfile:
    writer=csv.writer(csvfile)
    writer.writerows(lines)

