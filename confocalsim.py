import numpy as np
import matplotlib.pyplot as plt
from skimage import io

#Gaussian Beam intensity
def beamint(pos1,pos2,amp=1,w0=1):
    return amp * np.exp(0 - 2*(((pos1-pos2)**2).sum())/w0**2)

def makeTracks(num_particles, diff_const, numpixel, pixelsize, tau, subtau):
    #track coordinates in um
    tracks = []
    for i in range(num_particles):
        #format: t = [[xcoords],[ycoords]]
        if diff_const > 0:
            t = np.array([np.random.normal(loc=0.0,scale=np.sqrt(2*diff_const*tau/subtau),size=(numpixel*numpixel*subtau)).cumsum(),
                np.random.normal(loc=0.0,scale=np.sqrt(2*diff_const*tau/subtau),size=(numpixel*numpixel*subtau)).cumsum()])
        else:
            t = np.zeros((2,numpixel*numpixel*subtau))
        startpos = numpixel*pixelsize*np.random.rand(2)
        #print(startpos)
        t[0] += startpos[0]
        t[1] += startpos[1]
        #print(t)
        tracks.append(t)
    return tracks

def makeIntensityMap(tracks, numpixel, pixelsize, width):
    #all units in um and s    
    i_map = np.zeros((numpixel,numpixel)) #pixelmap/ image
    i_list = i_map.flatten()    #flattened pixelmap (list of intensities)
    
    #subtau: subdivisions for pixel dwell time
    if len(tracks) > 0:
        subtau = int(len(tracks[0][0])/(numpixel*numpixel))
        #make the intensity map
        for t in range(len(i_list)):
            xcoord = (t % numpixel)
            ycoord = (t-xcoord)/numpixel * pixelsize
            xcoord *= pixelsize
            
            for track in tracks:
                if subtau >= 1:
                    saver = 0
                    for j in range(subtau):
                         saver += beamint(np.array([xcoord,ycoord]),track[:,t*subtau+j],w0=width,amp=1/subtau)
                    i_list[t] += saver
                else:
                    i_list[t] += beamint(np.array([xcoord,ycoord]),track[:,0],w0=width,amp=1)

    i_map = i_list.reshape((numpixel,numpixel))
    return i_map


#create a noise image
def createNoiseImage(numpixel,stdNoise,avNoise=0):
    if stdNoise > 0:
        i_map = np.abs(np.random.normal(loc=avNoise,scale=stdNoise,size=[numpixel,numpixel]))
    else:
        i_map = np.zeros((numpixel,numpixel))+avNoise
    return i_map



def createFullImageFromTracks(tracks_mobile,tracks_immobile, numpixel,pixelsize,w0,cpp):
    #Additive White Gaussian Noise
    #Make the image
    #create mobile image map
    i_map_mobile = makeIntensityMap(tracks_mobile,numpixel,pixelsize,w0)
    #create immobile image map
    i_map_immobile = makeIntensityMap(tracks_immobile,numpixel,pixelsize,w0)
    #Add particle images together:
    i_map = i_map_mobile + i_map_immobile
    #Create the background:
    #create background noise

    #scale to CPP value
    i_map = i_map * cpp

    i_map += 1.2917
    #Ensure positive values
    i_map[i_map < 0] = 0
    #Add poisson noise
    i_map = np.random.poisson(i_map)
    #Have maximum pixel value
    i_map[i_map >= 2**16] = 2**16 - 1
    #make a 16-bit integer array
    i_map = i_map.astype(np.uint16)
    return i_map

#main function: create the full confocal image
def simSingleChannelConfocal(nmobile, nimmobile, diffconst, filename="confocalImage", cpp=10):
    pixelsize = 0.1 #um
    numpixel=100
    tau = 0.001 #s
    w0 = 0.3 #um
    subtau = 10

    #Build the particles:
    #create mobile tracks:
    tracks_mobile = makeTracks(nmobile, diffconst ,numpixel, pixelsize, tau, subtau)
    #create stationary particles:
    tracks_immobile = makeTracks(nimmobile, 0, numpixel, pixelsize, tau, 1)
    #Make the images from tracks
    i_map = createFullImageFromTracks(tracks_mobile,tracks_immobile,numpixel,pixelsize,w0,cpp)
    #save the image and the track plot in a PNG
    savePlot(i_map,tracks_mobile,tracks_immobile)
    saveTIFF(i_map,"single_TIFF")
    return 

def simBoundDualChannelConfocal(nmobile, nimmobile, diffconst, offset=0, filenames=["CH0","CH1"]):
    pixelsize = 0.1 #um
    numpixel=100
    tau = 0.001 #s
    w0green = 0.32 #um
    w0red = 0.37 #um
    cppgreen = 20 #kHz * 1ms (pixel dwell time)
    cppred = 20 #kHz * 1ms (pixel dwell time)
    subtau = 10

    if len(filenames) != 2:
        print("Must give 2 filenames")
        raise IndexError

    #Build the particles:
    #create mobile tracks:
    ch0_tracks_mobile = makeTracks(nmobile, diffconst ,numpixel, pixelsize, tau, subtau)
    #create stationary particles:
    ch0_tracks_immobile = makeTracks(nimmobile, 0, numpixel, pixelsize, tau, 1)
    #Make the images from tracks
    ch0_i_map = createFullImageFromTracks(ch0_tracks_mobile,ch0_tracks_immobile,numpixel,pixelsize,w0green,cppgreen)
    #save the image and the track plot in a PNG
    savePlot(ch0_i_map,ch0_tracks_mobile,ch0_tracks_immobile,filename=filenames[0]+"_PLOT")
    saveTIFF(ch0_i_map,filenames[0])
    #offset mobile and immobile tracks by "offset"
    ch1_tracks_mobile = []
    for track in ch0_tracks_mobile:
        rand_orientation = np.random.rand
        randorient = np.random.uniform(low=-np.pi,high=np.pi)
        randx = offset*np.cos(randorient)
        randy = offset*np.sin(randorient)
        saver1 = track[0] + randx*pixelsize
        saver2 = track[1] + randy*pixelsize
        ch1_tracks_mobile.append(np.array([saver1,saver2]))
    ch1_tracks_immobile = []
    for track in ch0_tracks_immobile:
        randorient = np.random.uniform(low=-np.pi,high=np.pi)
        randx = offset*np.cos(randorient)
        randy = offset*np.sin(randorient)
        saver1 = track[0] + randx*pixelsize
        saver2 = track[1] + randy*pixelsize
        ch1_tracks_immobile.append(np.array([saver1,saver2]))
    #Make the images from tracks
    ch1_i_map = createFullImageFromTracks(ch1_tracks_mobile,ch1_tracks_immobile,numpixel,pixelsize,w0red,cppred)
    #save the image and the track plot in a PNG
    savePlot(ch1_i_map,ch1_tracks_mobile,ch1_tracks_immobile,filename=filenames[1]+"_PLOT")
    saveTIFF(ch1_i_map,filenames[1])
    
    return [ch0_i_map, ch1_i_map]

def simFreeDualChannelConfocal(nmob, nimmob, diffconst, filenames=["CH0","CH1"],cpp=20):
    pixelsize = 0.1 #um
    numpixel=100
    tau = 0.001 #s
    w0 = [0.32,0.37] #um
    subtau = 10

    nmobile = nmob
    nimmobile = nimmob
    i_maps = []

    for i in range(2):
        #Build the particles:
        #create mobile tracks:
        tracks_mobile = makeTracks(nmobile[i], diffconst ,numpixel, pixelsize, tau, subtau)
        #create stationary particles:
        tracks_immobile = makeTracks(nimmobile[i], 0, numpixel, pixelsize, tau, 1)
        #Make the images from tracks
        i_map = createFullImageFromTracks(tracks_mobile,tracks_immobile,numpixel,pixelsize,w0[i],cpp)
        #save the image and the track plot in a PNG
        savePlot(i_map,tracks_mobile,tracks_immobile,filename=filenames[i]+"_Plot".format(i))
        saveTIFF(i_map,filenames[i])
        i_maps.append(i_map)
    return i_maps

def simParticleAndNoiseDualChannelConfocal(nmobile, nimmobile, diffconst, offset=[0,0], filenames=["CH0","CH1"],cpp=10):
    pixelsize = 0.1 #um
    numpixel=100
    tau = 0.001 #s
    w0 = 0.3 #um
    subtau = 10

    if len(filenames) != 2:
        print("Must give 2 filenames")
        raise IndexError

    #Build the particles:
    #create mobile tracks:
    ch0_tracks_mobile = makeTracks(nmobile, diffconst ,numpixel, pixelsize, tau, subtau)
    #create stationary particles:
    ch0_tracks_immobile = makeTracks(nimmobile, 0, numpixel, pixelsize, tau, 1)
    #Make the images from tracks
    ch0_i_map = createFullImageFromTracks(ch0_tracks_mobile,ch0_tracks_immobile,numpixel,pixelsize,w0,cpp)
    #save the image and the track plot in a PNG
    savePlot(ch0_i_map,ch0_tracks_mobile,ch0_tracks_immobile,filename=filenames[0]+"_PLOT")
    saveTIFF(ch0_i_map,filenames[0])
    #Noise Channel
    ch1_tracks_mobile = []
    ch1_tracks_immobile = []
    #Make the images from tracks
    ch1_i_map = createFullImageFromTracks(ch1_tracks_mobile,ch1_tracks_immobile,numpixel,pixelsize,w0,cpp)
    #save the image and the track plot in a PNG
    savePlot(ch1_i_map,ch1_tracks_mobile,ch1_tracks_immobile,filename=filenames[1]+"_PLOT")
    saveTIFF(ch1_i_map,filenames[1])

    return [ch0_i_map, ch1_i_map]

#save the image
def saveTIFF(i_map,filename):
    io.imsave(filename+".tif",i_map)
    return

#save the plot file
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
        ax2.plot(part[0,0],part[1,0],'ko')
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
    print("Testing the simulation")
    simSingleChannelConfocal(nmobile=20,nimmobile=5,diffconst=2.0)
    #simBoundDualChannelConfocal(nmobile=4,nimmobile=10,diffconst=2.0,offset=[10,0])
    #simFreeDualChannelConfocal(nmobile=0,nimmobile=10,diffconst=2.0)
