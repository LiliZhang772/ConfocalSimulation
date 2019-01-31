import confocalsim as cs
import datetime as dt
import os
import sys
import numpy as np


ofile = open("Logfile.log",'w')
time = dt.datetime.now().strftime("%Y%m%d")

#====== Mobile Particles ======
'''
print("mobile particles")
ofile.write("\n---Mobile Particles---\n")
dcounter = -1
for samplenum in range(100):
    dconsts = np.arange(10)*0.2 + 0.1
    number = "{:06d}".format(samplenum)
    if samplenum % 10 == 0:
        dcounter += 1
        print(dcounter)
    ofile.write("TIME="+time+"\n")
    ofile.write("SAMPLE="+number+"\n")
    ofile.write("D={:02.04f}".format(dconsts[dcounter])+"\n")
    dirname = time+number

    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        print("Directory exists! Overwriting now.")
        #continue

    sys.stdout.flush()

    i_map,tracks,parts = cs.simulateConfocal(nmobile=6,nimmobile=0,diffconst=dconsts[dcounter],fnoise=False,filename="testimage")
    for channelnum in range(2):
        ending = ".bin.counts_ch{:01d}.stack".format(channelnum)
        cs.saveTIFF(i_map,dirname+'/'+time+number+ending)
    cs.savePlot(i_map,tracks,parts,filename=time+number)
    sys.stdout.flush()


#====== Immobile Particles =====

print("immobile particles")
ofile.write("\n---Immobile Particles---\n")
wcounter = -1
for samplenum in range(100):
    widths = np.arange(10)*0.05 + 0.05
    number = "{:06d}".format(samplenum+10000)
    if samplenum % 10 == 0:
        wcounter += 1
        print(wcounter)
    ofile.write("TIME="+time+"\n")
    ofile.write("SAMPLE="+number+"\n")
    ofile.write("W0={:02.04f}".format(widths[dcounter])+"\n")
    dirname = time+number

    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        print("Directory exists! Overwriting now.")
        #continue

    sys.stdout.flush()

    i_map,tracks,parts = cs.simulateConfocal(nmobile=0,nimmobile=6,width=widths[wcounter],fnoise=False,filename="testimage")
    for channelnum in range(2):
        ending = ".bin.counts_ch{:01d}.stack".format(channelnum)
        cs.saveTIFF(i_map,dirname+'/'+time+number+ending)
    cs.savePlot(i_map,tracks,parts,filename=time+number)
    sys.stdout.flush()
'''

#===== Noise ========

print("Noise")
ofile.write("\n---Noise Particles---\n")
ncounter = -1
for samplenum in range(100):
    number = "{:06d}".format(samplenum+20000)
    noise = np.arange(10)*0.05+0.05

    if samplenum % 10 == 0:
        ncounter += 1
        print(ncounter)

    ofile.write("TIME="+time+"\n")
    ofile.write("SAMPLE="+number+"\n")
    ofile.write("NoiseSTD={:02.04f}".format(noise[ncounter])+"\n")
    dirname = time+number

    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        print("Directory exists! Overwriting now.")
        #continue

    sys.stdout.flush()

    i_map,tracks,parts = cs.simulateConfocal(nmobile=0,nimmobile=0,noise=noise[ncounter],fnoise=True,filename="testimage")
    for channelnum in range(2):
        ending = ".bin.counts_ch{:01d}.stack".format(channelnum)
        cs.saveTIFF(i_map,dirname+'/'+time+number+ending)
    cs.savePlot(i_map,tracks,parts,filename=time+number)
    sys.stdout.flush()


ofile.close()
