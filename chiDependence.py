# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 23:19:58 2019

@author: Markus
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def offsetBoxplot(filename,coordinate,zerocoordinate,particletype):
    offsetarray = []
    dirname = filename+zerocoordinate
    df_raw = pd.read_csv(dirname+"/SPOTS-ch0.csv")
    df = df_raw[np.logical_and(df_raw["Good?"] == 1,df_raw["Above threshold?"] == "yes")]
    df = df[df['streak or spot?']=='spot']
    offsetarray.append(df['Correlation Coefficient'])
    
    for offset in range(1,4,1):
        dirname = filename+coordinate+"{:1d}".format(offset)
        
        df_raw = pd.read_csv(dirname+"/SPOTS-ch0.csv")
        
        df = df_raw[np.logical_and(df_raw["Good?"] == 1,df_raw["Above threshold?"] == "yes")]
        
        df = df[df['streak or spot?']=='{:}'.format(particletype)]
        df["offsetx"] = offset
    
        offsetarray.append(df['Correlation Coefficient'])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.boxplot(offsetarray,positions = np.arange(4),patch_artist=True,notch=False,showmeans=True)
    ax.set_xlabel(coordinate+" offset [px]")
    ax.set_ylabel("Cross-correlation coefficient")
    ax.set_ylim([-0.4,1.2])
    ax.grid()
    plt.savefig(coordinate+"-Offset-{:}.png".format(particletype))
    plt.show()
    
offsetBoxplot("ExpPositiveXCorr-","X","X0","spot")
offsetBoxplot("ExpPositiveXCorr-","Y","X0","spot")
offsetBoxplot("ExpStreaksPositiveXCorr-","X","Y0","streak")
offsetBoxplot("ExpStreaksPositiveXCorr-","Y","Y0","streak")

offsetarray = []
dirnames = ["ExpSpotsRandomNoise","ExpSpotsRandomSpots","ExpStreaksRandomNoise","ExpStreaksRandomStreaks"]
for elem in dirnames:
    df_raw = pd.read_csv(elem+"/SPOTS-ch0.csv")
    df = df_raw[np.logical_and(df_raw["Good?"] == 1,df_raw["Above threshold?"] == "yes")]
    df = df[df['streak or spot?']=='spot']
    offsetarray.append(df['Correlation Coefficient'])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.boxplot(offsetarray[0:2],positions = np.arange(2),patch_artist=True,notch=False,showmeans=True,labels=["vs Noise","vs Spots"])
ax.set_xlabel("")
ax.set_ylabel("Cross-correlation coefficient")
ax.set_ylim([-1,1.2])
ax.grid()
plt.savefig("SpotsNoise.png")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.boxplot(offsetarray[2:4],positions = np.arange(2),patch_artist=True,notch=False,showmeans=True,labels=["vs Noise","vs Streaks"])
ax.set_xlabel("")
ax.set_ylabel("Cross-correlation coefficient")
ax.set_ylim([-1,1.2])
ax.grid()
plt.savefig("StreaksNoise.png")
plt.show()
