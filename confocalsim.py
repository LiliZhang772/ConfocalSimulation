import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from skimage import io


def beamint(pos1,pos2,amp=1,w0=1):
    return amp * np.exp(0 - (((pos1-pos2)**2).sum())/w0**2)

def createMobileImage(numpixel=100,pixelsize=0.1,tau=0.001,w0=0.3,N=10,D=1e-9,subtau=10):
    #all units in um and s    
    i_map = np.zeros((numpixel,numpixel)) #pixelmap/ image
    i_list = i_map.flatten()    #flattened pixelmap (list of intensities)

    #track coordinates in um
    tracks = []
    for i in range(N):
        #format: t = [[xcoords],[ycoords]]
        t = np.array([np.random.normal(loc=0.0,scale=np.sqrt(2*D*tau/subtau),size=(numpixel*numpixel*subtau)).cumsum(),
                np.random.normal(loc=0.0,scale=np.sqrt(2*D*tau/subtau),size=(numpixel*numpixel*subtau)).cumsum()])
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
            saver = 0
            for j in range(subtau):
                 saver += beamint(np.array([xcoord,ycoord]),track[:,t*subtau+j],w0=w0,amp=1/subtau)
                 #saver /= subtau
            i_list[t] += saver

    i_map = i_list.reshape((numpixel,numpixel))
    #print(i_map)
    return i_map, tracks

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
            i_list[t] += beamint(np.array([xcoord,ycoord]),track,w0=w0,amp=1)

    i_map = i_list.reshape((numpixel,numpixel))
    #print(i_map)
    return i_map, tracks

def createNoiseImage(numpixel=100,noise=10.):
    i_map = np.abs(np.random.normal(loc=0.1,scale=noise,size=[100,100]))
    return i_map


if __name__=="__main__":
    numpixel = 100
    pixelsize = 0.1 #um
    tau = 0.001 #s
    w0 = 0.3 #um
    Nmob = 5 #number of mobile particles
    Nstat = 5 #number of stationary particles
    Nnoise = 0.2 #noise amount
    D = 0.7 #um^2/s
    subtau = 10

    cpp = 300 #kHz
    
    i_map_mob, tracks = createMobileImage(numpixel,pixelsize,tau,w0,Nmob,D,subtau)
    i_map_stat, parts = createStationaryImage(numpixel,pixelsize,w0,Nstat)
    i_map_noise = createNoiseImage(numpixel,Nnoise)
    
    i_map = i_map_mob + i_map_stat + i_map_noise

    i_map = i_map * cpp 

    i_map[i_map > 0] += np.random.normal(loc=0.0,scale=np.sqrt(i_map[i_map > 0]))

    i_map[i_map < 0] = 0
    i_map[i_map >= 2**16] = 2**16 - 1

    i_map = i_map.astype(np.uint16)

    
    io.imsave("confocalImage.tif",i_map)
    
    fig1 = plt.figure(figsize=(12,7))
    ax1 = fig1.add_subplot(121)
    ax1.imshow(i_map)
    ax1.set_xlabel(r"x [$\mu$m]")
    ax1.set_ylabel(r"y [$\mu$m]")
    ax2 = fig1.add_subplot(122)
    for track in tracks:
        ax2.plot(track[0],track[1],'y-')
    for track in tracks:
        ax2.plot(track[0,0],track[1,0],'go')
        ax2.plot(track[0,-1],track[1,-1],'ro')
    for part in parts:
        ax2.plot(part[0],part[1],'ko')
    ax2.plot([-1,-2],[-1,-2],'y-',label="mobile tracks")
    ax2.plot([-1],[-1],'ko',label="immobile")
    ax2.plot([-1],[-1],'go',label="mobile start")
    ax2.plot([-1],[-1],'ro',label="mobile end")
    ax2.axis('scaled')
    ax2.set_xlim((0,int(numpixel*pixelsize)))
    ax2.set_ylim((int(numpixel*pixelsize),0))
    ax2.set_xlabel(r"x [$\mu$m]")
    ax2.set_ylabel(r"y [$\mu$m]")
    ax2.legend(bbox_to_anchor=(0,0,1.0,-0.13),loc=2,ncol=2,borderaxespad=0.,mode='expand')
    plt.savefig("confocalPlot.png",dpi=70)
    plt.show()

