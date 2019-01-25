import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


def beamint(pos1,pos2,amp=1,w0=1):
    return amp * np.exp(0 - (((pos1-pos2)**2).sum())/w0**2)

def createMobileImage(numpixel=100,pixelsize=0.1,tau=0.001,w0=0.3,N=10,D=1e-9):
    #all units in um and s    
    i_map = np.zeros((numpixel,numpixel)) #pixelmap/ image
    i_list = i_map.flatten()    #flattened pixelmap (list of intensities)

    #track coordinates in um
    tracks = []
    for i in range(N):
        #format: t = [[xcoords],[ycoords]]
        t = np.array([np.random.normal(loc=0.0,scale=np.sqrt(2*D*tau),size=(numpixel*numpixel)).cumsum(),
                np.random.normal(loc=0.0,scale=np.sqrt(2*D*tau),size=(numpixel*numpixel)).cumsum()])
        #print(t)
        startpos = numpixel*pixelsize*np.random.rand(2)
        #print(startpos)
        t[0] += startpos[0]
        t[1] += startpos[1]
        #print(t)
        tracks.append(t)

    for t in range(len(i_list)):
        xcoord = (t % numpixel)
        ycoord = (t-xcoord)/numpixel * pixelsize
        xcoord *= pixelsize
        
        for track in tracks:
            i_list[t] += beamint(np.array([xcoord,ycoord]),track[:,t],w0=w0)

    i_map = i_list.reshape((numpixel,numpixel))
    #print(i_map)
    return i_map

def createStationaryImage(numpixel=100,pixelsize=0.1,w0=0.3,N=10):
    i_map = np.zeros((numpixel,numpixel)) #pixelmap/ image
    i_list = i_map.flatten()    #flattened pixelmap (list of intensities)

    #track coordinates in um
    tracks = []
    for i in range(N):
        #format: t = [[xcoords],[ycoords]]
        startpos = numpixel*pixelsize*np.random.rand(2)
        #print(startpos)
        tracks.append(startpos)

    for t in range(len(i_list)):
        xcoord = (t % numpixel)
        ycoord = (t-xcoord)/numpixel * pixelsize
        xcoord *= pixelsize
        
        for track in tracks:
            i_list[t] += beamint(np.array([xcoord,ycoord]),track,w0=w0)

    i_map = i_list.reshape((numpixel,numpixel))
    #print(i_map)
    return i_map

def createNoiseImage(numpixel=100):
    i_map = np.abs(np.random.normal(loc=1.,scale=10.,size=[100,100]))
    return i_map


if __name__=="__main__":
    numpixel = 100
    pixelsize = 0.1 #um
    tau = 0.001 #s
    w0 = 0.2 #um
    N = 0 #number of particles
    D = 1e-9 #um^2/s
    
    #i_map = createImage(numpixel,pixelsize,tau,w0,N,D)
    #i_map = createNoise(numpixel)
    i_map = createStationaryImage()
    
    
    plt.imshow(i_map)
    plt.show()

