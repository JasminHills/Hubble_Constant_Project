# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 14:51:42 2019

@author: hills
Estimate of Hubble's constant: 70.74646 +/- 4.8241 km/s/Mpc

"""
#Estimate of Hubble's constant: 70.74646 +/- 4.8241 km/s/Mpc
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import glob as glob
def readfilez():
  count =1
  good_wv= sp.empty((30,999))
  good_int= sp.empty((30, 999))
  good_obs= sp.empty((30))
  filelist= glob.glob( 'Data\Halpha_spectral_data\Halpha_spectral_data-?*.csv')

  line1split=""
  for i in filelist:
    val=pd.read_csv(i, delimiter=',', header=2)
    val=pd.DataFrame.to_numpy(val)
    with open(i, 'r') as csv:
        line1=csv.readline()
        line1splt=line1.split(',')
        resp=line1splt[3].strip()
        obs=line1splt[2].strip()
        obsv= obs.split(' ')
        obsv[1].strip()
        num=int(obsv[1])
        if resp=="Instrument Response: Good":
          wvlnth= val.T[0]
          intensity=val.T[1]
          plt.scatter(wvlnth,intensity)
          good_wv[count-1]= wvlnth
          good_int[count-1]= intensity
          good_obs[count-1]=num
          count=count+1
  return good_wv, good_int, good_obs    
def findmax(arr):
    arrsz=len(arr)
    print(arrsz)
    val, biggest= 0,0
    for x in range(0, arrsz-1):
       if arr[x]>biggest:
              val=x
              biggest=arr[x]
    return val
def initialguess(wvlnth, intensity):
  length = len(intensity)-1
  slope= (intensity[length]-intensity[0])/(wvlnth[length]-wvlnth[0])
  intercept= intensity[0]-wvlnth[0]*slope
  distance= intensity-((wvlnth)*slope+intercept)
  val=findmax(distance)
  gass= wvlnth[val]
  a=intensity[val]
  sig=sp.std(intensity) 
  init=[a,gass,sig,slope,intercept]
  return init

def maxiTester(m,c,y,x):
    lenth=len(x)
    totaldist=sp.empty(lenth-3)
    for count in range(6,lenth-6):
        for a in range(0,7):
            totaldist[count]=totaldist[count]+y[count-a]-(x[count-a]*m-c)
            totaldist[count]=y[count+a]-(x[count+a]*m-c)
    return x[findmax(totaldist)]

def fitfunc(x, a,mu,sig,m,c):
  gausie=a*sp.exp(-(x-mu)**2/(2*sig**2))
  line=m*x+c
  return gausie+line

def graphingp1(wvlnth, intensity):
  initial_guess=initialguess(wvlnth, intensity)
  wvlnth=sp.nan_to_num(wvlnth)
  intensity=sp.nan_to_num(intensity)
  
  po, po_cov=sp.optimize.curve_fit(fitfunc,wvlnth,intensity,initial_guess)
  print(po[0], po[1], po[2], po[3], po[4])
  
  plt.plot(wvlnth,intensity)
  plt.plot(wvlnth,fitfunc(wvlnth, po[0], po[1], po[2], po[3], po[4]) )
  plt.show()
  maxival= findmax(intensity)
  maxival=maxi(po[3], po[4],intensity, wvlnth)
  return maxival

readfilez()
v=sp.empty(30)
d= sp.empty(30)
good_wv,good_int, good_obs =readfilez()
print(len(good_obs))
for count in range(0,30):
     if not good_obs[count]== None:
        number=graphingp1(good_wv[count], good_int[count])
        print("num:",number)
        if number==759.5:
           number=guesses[count]
        ratsq= (number/656.28)**2
        v[count]=3*(10**8)*((ratsq-1)/(1+ratsq))
        print v[count]
        obsv=good_obs[count]
        print("obsv",obsv)
        distances=pd.read_csv('Data\Distance_Mpc.csv', delimiter=None, header=1)
        distances=pd.DataFrame.to_numpy(distances)
        observation=distances.T[0]
        distance=distances.T[1]
        length=observation.size
        for x in range(0, length):
            if obsv==observation[x]:
                d[count]=distance[x]

v=v/1000
fit_Hubble, cov_Hubble = sp.polyfit(d, v, 1, cov = 1)     
vGalaxy = sp.poly1d(fit_Hubble)
plt.plot(d, vGalaxy(d), color = 'blue')
plt.title('Redshift Speed vs Distance')
plt.xlabel('Distances (Mpc)')
plt.ylabel('Velocities (Km/s)')
print("Equation: ", fit_Hubble," Uncertainty: ",sp.sqrt(cov_Hubble[0,0]) )
#plt.plot(d,pfit(d), label='Distance vs Velocity')      
 DONEEEEEEEEEEEEEEE
  
#%%fit_1,cov_1=sp.polyfit(d,v,1,cov=True)

#%%sig_0=sp.sqrt(cov_1(0,0))
#%%sig_1=sp.sqrt(cov_1(1,1))


