import confocalsim as cs
import datetime as dt
import os
import sys
import numpy as np


time = dt.datetime.now().strftime("%Y%m%d")

workingDir = os.path.abspath('.')


serialnum = 42000
diff_const = 10.0

expDir = "RandomSpotstoStreaks_20to20"
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
    cs.simFreeDualChannelConfocal(nmob=[0,20], nimmob=[20,0], diffconst=diff_const, filenames=fnames,cpp=20)
        
ofile.write("END_SAMPLE="+number+"\n")
ofile.write("DIFF_CONST={:1.01f}".format(diff_const)+"\n")
ofile.write("Noise output\n")
ofile.write("Offset=[0,3]\n")
ofile.write("Done"+"\n")
ofile.close()
