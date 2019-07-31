import confocalsim as cs
import datetime as dt
import os
import sys
import numpy as np


time = dt.datetime.now().strftime("%Y%m%d")

workingDir = os.path.abspath('.')


xoffset = 3
serialnum = 14000

expDir = "ExpStreaksPositiveXCorr-X{:1d}".format(xoffset)
if not os.path.isdir(expDir):
    os.mkdir(expDir)
os.chdir(expDir)
ofile = open("Logfile.log",'w')

#====== Immobile Particles =====
print("immobile particles")
ofile.write("\n---Immobile Particles---\n")
number = "{:06d}".format(serialnum)
ofile.write("TIME="+time+"\n")
ofile.write("START_SAMPLE="+number+"\n")
for samplenum in range(0,100,1):
    print(samplenum)
    number = "{:06d}".format(samplenum+serialnum)

    dirname = time+number
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        print("Directory exists! Overwriting now.")
    sys.stdout.flush()
    
    fnames = []
    for channelnum in range(2):
        ending = ".bin.counts_ch{:01d}.stack".format(channelnum)
        fnames.append(dirname+'/'+time+number+ending)
    cs.simBoundDualChannelConfocal(nmobile=15,nimmobile=0,diffconst=2.0,offset=[xoffset,0],filenames=fnames)
        
ofile.write("END_SAMPLE="+number+"\n")
ofile.write("XOFFSET={:02.04f}".format(xoffset)+"\n")
ofile.close()
