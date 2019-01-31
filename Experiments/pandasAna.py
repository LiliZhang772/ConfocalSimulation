import os
import sys
import re
import numpy as np
import pandas as pd

#BName = "STREAKS"
#BName = "SPOTS"
BName = "NOISE"

savepath = "C:/Users/Markus/LittleHelpers/ConfocalSimulation/Analysis"
if not os.path.isdir(savepath):
    os.mkdir(savepath)
ch0savecsv = "{:}-ch0.csv".format(BName)
ch1savecsv = "{:}-ch1.csv".format(BName)
saveinfo = "{:}-ch1-info.txt".format(BName)

of = open(os.path.join(savepath,saveinfo),'w')

if BName == "STREAKS":
    dirpattern = re.compile(r"\d{8}00\d{4}")
if BName == "SPOTS":
    dirpattern = re.compile(r"\d{8}01\d{4}")
if BName == "NOISE":
    dirpattern = re.compile(r"\d{8}02\d{4}")
ch0df = pd.DataFrame()
dirindex = []
ch0dflist = []
ch1df = pd.DataFrame()
ch1dflist = []
dircount = 0
filecount = 0
for (dirpath,dirnames,filenames) in os.walk('.'):
    if dircount == 0:
        dirindex = [i for i,word in enumerate(dirnames) if dirpattern.match(word)]
        dircount = len(dirindex)
        print(dircount)
        continue

    print(dirpath)
    of.write("PATH={:}\n".format(dirpath))
    dirname = os.path.split(dirpath)[1]
    if not dirpattern.match(dirname):
        print(dirname + " doesn't Match!")
        continue
    fp0 = re.compile(r"{:}.bin.counts_ch0.stack.tifslice#\d.txt".format(dirname))
    fp1 = re.compile(r"{:}.bin.counts_ch1.stack.tifslice#\d.txt".format(dirname))
    ch0idx = [i for i, word in enumerate(filenames) if fp0.match(word)]
    ch1idx = [i for i, word in enumerate(filenames) if fp1.match(word)]

    if len(ch0idx) == 1 and len(ch1idx) == 1:
        filecount += 1
        #print(filecount)
        ch0save= pd.read_csv(os.path.join(dirpath,filenames[ch0idx[0]]),sep='\t')#,index_col=' ')
        ch1save= pd.read_csv(os.path.join(dirpath,filenames[ch1idx[0]]),sep='\t')#,index_col=' ')
        ch0save['date'] = [int(dirname[:8])]*len(ch0save)
        ch1save['date'] = [int(dirname[:8])]*len(ch1save)
        ch0df = ch0df.append(ch0save)
        ch1df = ch1df.append(ch1save)
        print(filenames[ch0idx[0]])
        print(filenames[ch1idx[0]])
        of.write("CHANNEL0={:}\n".format(filenames[ch0idx[0]]))
        of.write("CHANNEL1={:}\n".format(filenames[ch1idx[0]]))
    else:
        print("False")
    print("="*72)
    of.write("="*72+"\n")
print(filecount)
of.write("\n\n")
of.write("SUMMARY\n=======\n")
of.write("FOUND_FOLDERS={:}\n".format(dircount))
of.write("FOUND_FILES={:}\n".format(filecount))

ch0df.to_csv(os.path.join(savepath,ch0savecsv),sep=',',header=True)
ch1df.to_csv(os.path.join(savepath,ch1savecsv),sep=',',header=True)
of.write("OUTPUT-CH0={:}\n".format(os.path.join(savepath,ch0savecsv)))
of.write("OUTPUT-CH1={:}\n".format(os.path.join(savepath,ch1savecsv)))

