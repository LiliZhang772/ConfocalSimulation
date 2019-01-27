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

    '''
    def gaussian(pos,w,amp=1):
        xcent = 5 + pos[0] - int(pos[0])
        ycent = 5 + pos[0] - int(pos[0])
        gmap = np.zeros((numpixel,numpixel))
        for i in range(len(gmap)):
            for j in range(len(gmap[i])):
                gmap[i,j] = amp * np.exp(0 - (xcent-j-5)**2/w[1]**2+(ycent-i-5)**2/w[0]**2)
        return gmap
    '''

    #track coordinates in um
    tracks = []
    for i in range(N):
        #format: t = [[xcoords],[ycoords]]
        startpos = numpixel*pixelsize*np.random.rand(2)
        #print(startpos)
        tracks.append(startpos)

    i_list = i_map.flatten()    #flattened pixelmap (list of intensities)
    for t in range(len(i_list)):
        xcoord = (t % numpixel)
        ycoord = (t-xcoord)/numpixel * pixelsize
        xcoord *= pixelsize
        
        for track in tracks:
            i_list[t] += beamint(np.array([xcoord,ycoord]),track,w0=w0,amp=1)

    i_map = i_list.reshape((numpixel,numpixel))

    return i_map, tracks

def createNoiseImage(numpixel=100,noise=10.):
    i_map = np.abs(np.random.normal(loc=0.1,scale=noise,size=[100,100]))
    return i_map


def simulateConfocal(nmobile,nimmobile,fnoise,numpix=100,diffconst=1.0,width=0.3,noise=0.2,filename="confocalImage"):
    numpixel = numpix
    pixelsize = 0.1 #um
    tau = 0.001 #s
    w0 = width #um
    Nmob = nmobile #number of mobile particles
    Nstat = nimmobile #number of stationary particles
    Nnoise = noise #noise amount
    D = diffconst #um^2/s
    subtau = 10

    cpp = 300 #kHz
    
    i_map_mob, tracks = createMobileImage(numpixel,pixelsize,tau,w0,Nmob,D,subtau)
    i_map_stat, parts = createStationaryImage(numpixel,pixelsize,w0,Nstat)
    if fnoise:
        i_map_noise = createNoiseImage(numpixel,Nnoise)
    else:
        i_map_noise = np.zeros((numpixel,numpixel))
    
    i_map = i_map_mob + i_map_stat + i_map_noise

    i_map = i_map * cpp 

    i_map[i_map > 0] += np.random.normal(loc=0.0,scale=np.sqrt(i_map[i_map > 0]))

    i_map[i_map < 0] = 0
    i_map[i_map >= 2**16] = 2**16 - 1

    i_map = i_map.astype(np.uint16)

    return i_map,tracks,parts
    
def saveTIFF(i_map,filename):
    io.imsave(filename+".tif",i_map)
    return
    
def savePlot(i_map,tracks,parts,pixelsize=0.1,filename="exampleImage"):
    numpixel = len(i_map)
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
    plt.savefig(filename+".png",dpi=70)
    #plt.show()
    plt.close("all")
    return

if __name__=="__main__":
    imap,tracks,parts = simulateConfocal(nmobile=0,nimmobile=3,fnoise=False)
    savePlot(imap,tracks,parts)

