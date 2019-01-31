
import os
import sys
import re
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import itertools

doublesize=20
singlesize=16

rc('font', family='serif')
rc('font', serif=['Computer Modern'])
rc('font', size=doublesize)
#rc('text', usetex=True)


def read_all_particles(conc):
    allParticles = pd.DataFrame()

    for baxconc in conc:
        filename0 = "{:}-ch0.csv".format(baxconc)
        filename1 = "{:}-ch1.csv".format(baxconc)

        ch0save = pd.read_csv(filename0,sep=',')
        ch1save = pd.read_csv(filename1,sep=',')

        ch0save['Baxconc'] = baxconc
        ch1save['Baxconc'] = baxconc

        ch0save['Channel'] = 0
        ch1save['Channel'] = 1
        
        condition0 = [True] * len(ch0save)
        condition1 = [True] * len(ch1save)
        condition0 = np.logical_and(ch0save['Good?']==1,condition0)
        condition1 = np.logical_and(ch1save['Good?']==1,condition1)
        condition0 = np.logical_and(ch0save['Above threshold?']=='yes',condition0)
        condition1 = np.logical_and(ch1save['Above threshold?']=='yes',condition1)
        
        '''
        if baxconc == 500:
            condition0 = np.logical_and(condition0,ch0save['date'] != 20130808)
            condition1 = np.logical_and(condition1,ch1save['date'] != 20130808)
        '''
        allParticles = allParticles.append(ch0save[condition0])
        allParticles = allParticles.append(ch1save[condition1])

    return allParticles

def norm_factors(conc):

    critlist = []
    countlist = []
    conckeep = []


    for baxconc in conc:
        finfo = open("Bax{:}-ch1-info.txt".format(baxconc),'r')
        contents = finfo.read()
        finfo.close()

        saver = re.findall(r'CHANNEL0=[0-9]{,8}',contents)
        saverset = set(saver)
            
        for elem in set(saver):
            dateint = int(elem[9:])
            if baxconc == 500 and dateint == 20130808:
                #print("YES HELLO HELLO ++++++++++++++++++++++++++++++++++++++++++++++++++")
                saverset.remove(elem)
                continue
            countlist.append(contents.count(elem))
            critlist.append(int(elem[9:]))
            conckeep.append(baxconc)
        #print("Bax Concentration: {:}".format(baxconc))
        #print("Dates:")
        #print(saverset)

    return critlist,countlist,conckeep



def read_data(Bax):
    filename0 = "Bax{:}-ch0.csv".format(Bax)
    filename1 = "Bax{:}-ch1.csv".format(Bax)

    ch0df = pd.read_csv(filename0,sep=',')
    ch1df = pd.read_csv(filename1,sep=',')

    anadf1 = ch0df[np.logical_and(ch0df['streak or spot?'] == 'spot',ch0df['Good?']==1)]
    anadf2 = ch0df[np.logical_and(ch0df['streak or spot?'] == 'streak',ch0df['Good?']==1)]
    anadf3 = ch1df[np.logical_and(ch1df['streak or spot?'] == 'spot',ch1df['Good?']==1)]
    anadf4 = ch1df[np.logical_and(ch1df['streak or spot?'] == 'streak',ch1df['Good?']==1)]
    return anadf1,anadf2,anadf3,anadf4



def fithistogram(histo,savename,start=4,dist=4.5,sigmamax=2,path='./',fitall=True):
    xvals = (histo[1][1:] + histo[1][:-1])/2   
    xvrange = xvals < 30

    def gaussian(x,a,mu,sig):
        return a*np.exp(0 - (x - mu)**2/(2.*sig**2))

    def fitcurve(x,start,dist,amplitudes,sig):
        fitc = np.zeros(len(x))
        for i in range(len(amplitudes)):
            fitc += gaussian(x,amplitudes[i],(i)*dist+start,sig)
        return fitc


    def fitthis(x,start,dist,a1,a2,a3,a4,a5,sig):
        amps = [a1,a2,a3,a4,a5]
        return fitcurve(x,start,dist,amps,sig)

    def fitthis2(x,a1,a2,a3,a4,a5,sig):
        amps = [a1,a2,a3,a4,a5]
        return fitcurve(x,start,dist,amps,sig)

    try:
        if fitall:
            popt,pcov = curve_fit(fitthis,xvals[xvrange],histo[0][xvrange],
                    bounds=[[4,3,1,1,1,0.1,0.1,0.1],[10,7,np.inf,np.inf,np.inf,np.inf,10,sigmamax]])
        else:
            popt,pcov = curve_fit(fitthis2,xvals[xvrange],histo[0][xvrange],
                    bounds=[[0.1,0.1,0.1,0.1,0.1,0.1],[np.inf,np.inf,np.inf,np.inf,10,sigmamax]])
        sigs = [popt[-1]] * 3
        cols = ["Start","Distance","A1","A2","A3","A4","A5","sigma"]
        if fitall:
            pcov = list(np.sqrt(np.diagonal(pcov)))
        else:
            popt = [start,dist] + list(popt)
            pcov = [0,0] + list(np.sqrt(np.diagonal(pcov)))
    except RuntimeError:
        print("Fit did not converge")
        amplis = [23,5,7,5,5]
        popt = [5,9] + amplis + [1]
        pcov = [0]*8

    widths = 0.9*(histo[1][1]-histo[1][0])
    xvalsalt = np.arange(1000)/1000*(xvals[-1] - xvals[0])+xvals[0]
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111)
    ax.bar(xvals-(xvals[1]-xvals[0])/2,histo[0],width=widths,align='center',alpha=0.5,color='blue')
    ax.plot(xvalsalt,fitthis(xvalsalt,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7]),'k-')
    sumfun = np.zeros(len(histo[0]))
    for i in range(len(popt)-3):
        gg = gaussian(xvalsalt,popt[2+i],(i)*popt[1]+popt[0],popt[-1])
        ax.plot(xvalsalt,gg)
    ax.set_xlabel("Intensity [kHz]")
    ax.set_ylabel("Count")
    ax.set_xlim((0,50))
    plt.savefig(path + savename+"-fit.png")
    plt.clf()
    return popt,pcov


def make_stocheo(data,savename,bincount=81,normfactor=1,path='./',fitall=True,start=4,dist=5):
    if not os.path.isdir(path):
        os.mkdir(path)
    histo = np.histogram(data['a'].values,bins=bincount,range=(0,70))
    if fitall:
        popt,pcov = fithistogram(histo,savename,path=path,fitall=True)
    else:
        popt,pcov = fithistogram(histo,savename,path=path,fitall=False,start=start,dist=dist)
    xvals = (histo[1][1:] + histo[1][:-1])/2   
    bins = []
    binsize = popt[1]/2.
    val = popt[0] - binsize
    while val < xvals.max():
        bins.append(val)
        val += 2*binsize
    bins.append(val)
    thishist = np.histogram(data['a'].values,bins=bins)
    newxvals = (thishist[1][1:]+thishist[1][:-1])/2.
    #newxvals /= popt[1]
    #xvals /= popt[1]
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111)
    ax.bar(newxvals-(newxvals[1]-newxvals[0])/2,thishist[0]/normfactor,width=0.9*(newxvals[1]-newxvals[0]),color='blue',yerr=np.sqrt(thishist[0]-thishist[0]**2/thishist[0].sum())/normfactor)
    ax.bar(xvals-(xvals[1]-xvals[0])/2,histo[0]/normfactor,width=0.9*(xvals[1]-xvals[0]),color='green')
    #plt.axis([0,8,0,thishist[0].max()+5])
    ax.set_xlabel("Intensity [kHz]")
    ax.set_ylabel("Count")
    ax.set_xlim((0,50))
    plt.savefig(path + savename+"-stocheo.png")
    plt.clf()
    return list(popt) + list(pcov) + list(thishist[0][0:9]) + list(np.sqrt(thishist[0][0:9])) #list(np.sqrt(thishist[0][0:5] - thishist[0][0:5]**2/thishist[0].sum())/normfactor)


def crosscorrelation(data,savename,cutoff=0.3,bincount=81,normfactor=1,path='./'):
    histo = np.histogram(data['Correlation Coefficient'],bins=bincount,range=(-1,1))
    histoerr = np.sqrt(histo[0] - histo[0]**2/histo[0].sum())
    xvals = (histo[1][1:] + histo[1][:-1])/2   
    xdispl = xvals - (xvals[1]-xvals[0]) / 2
    #matplotlib.rcParams.update({'font.size': 30})
    figcross = plt.figure(figsize=(9,7))
    ax1 = figcross.add_subplot(111)
    ax1.bar(xdispl,histo[0]/normfactor,width=0.9*(xvals[1]-xvals[0]),yerr=histoerr/normfactor)
    ax1.axvline(cutoff-(xvals[1]-xvals[0])/2,color='red')
    ax1.axvline(0,color='black')
    ax1.set_xlabel(r"Cross-correlation coefficient $\chi$")
    ax1.set_ylabel("Particle count")
    ax1.set_xlim([-1,1])
    plt.savefig(path + savename + "-crosscorr.png")
    plt.clf()
    highcount = len(data[data["Correlation Coefficient"] >= cutoff])
    allcount = len(data) 
    return allcount,highcount, np.sqrt(highcount) #np.sqrt(highcount - highcount ** 2/allcount)


def main():
    immobile_dir = "./ImmobileStocheo/"
    mobile_dir = "./MobileStocheo/"
    corr_dir = "./Cross-Correlation/"
    if not os.path.isdir("./ImmobileStocheo/"):
        os.mkdir(immobile_dir)
    if not os.path.isdir(mobile_dir):
        os.mkdir(mobile_dir)
    if not os.path.isdir(corr_dir):
        os.mkdir(corr_dir)

    columns = ['Ioffset','Imono','A1','A2','A3','A4','A5','sig','err_Ioffset','err_Imono','err_A1','err_A2','err_A3','err_A4','err_A5','err_sig','Binsum1','Binsum2','Binsum3','Binsum4','Binsum5','Binsum6','Binsum7','Binsum8','Binsum9','err_Binsum1','err_Binsum2','err_Binsum3','err_Binsum4','err_Binsum5','err_Binsum6','err_Binsum7','err_Binsum8','err_Binsum9','all_correlation','high_correlation','all_correlation_alt','high_correlation_alt','err_high_corr','err_high_corr_alt','Ntotal','normfactor','Itotal','Baxconc','Bidconc','Baxints','Bidints']

    #conc = ['NOISE','STREAKS','SPOTS']
    conc = ['NOISE']
    allparts = read_all_particles(conc)
    mobileparts = allparts[allparts['streak or spot?']=='streak']
    immobileparts = allparts[allparts['streak or spot?']=='spot']
    undefinedparts = allparts[allparts['streak or spot?']=='Undefined']

    rc('font',size=singlesize)
    figex = plt.figure(figsize=(9,7))
    axex = figex.add_subplot(111)
    axex.plot(mobileparts['wx'],mobileparts['wy'],marker='o',linestyle='',color='orange',label="mobile particles")
    axex.plot(immobileparts['wx'],immobileparts['wy'],marker='o',linestyle='',color='green',label="immobile particles")
    axex.plot(undefinedparts['wx'],undefinedparts['wy'],marker='o',linestyle='',color='grey',label="undefined particles")
    axex.set_xlabel(r"$w_x$ [$nm$]")
    axex.set_ylabel(r"$w_x$ [$nm$]")
    axex.set_xlim((0,500))
    axex.set_ylim((0,500))
    axex.legend(loc='best')
    figex.savefig("eccentricityMap-NOISE.png")
    plt.clf()
    rc('font',size=doublesize)
    
    conc = ['STREAKS']
    allparts = read_all_particles(conc)
    mobileparts = allparts[allparts['streak or spot?']=='streak']
    immobileparts = allparts[allparts['streak or spot?']=='spot']
    undefinedparts = allparts[allparts['streak or spot?']=='Undefined']

    rc('font',size=singlesize)
    figex = plt.figure(figsize=(9,7))
    axex = figex.add_subplot(111)
    axex.plot(mobileparts['wx'],mobileparts['wy'],marker='o',linestyle='',color='orange',label="mobile particles")
    axex.plot(immobileparts['wx'],immobileparts['wy'],marker='o',linestyle='',color='green',label="immobile particles")
    axex.plot(undefinedparts['wx'],undefinedparts['wy'],marker='o',linestyle='',color='grey',label="undefined particles")
    axex.set_xlabel(r"$w_x$ [$nm$]")
    axex.set_ylabel(r"$w_x$ [$nm$]")
    axex.set_xlim((0,500))
    axex.set_ylim((0,500))
    axex.legend(loc='best')
    figex.savefig("eccentricityMap-STREAKS.png")
    plt.clf()
    rc('font',size=doublesize)
    
    conc = ['SPOTS']
    allparts = read_all_particles(conc)
    mobileparts = allparts[allparts['streak or spot?']=='streak']
    immobileparts = allparts[allparts['streak or spot?']=='spot']
    undefinedparts = allparts[allparts['streak or spot?']=='Undefined']

    rc('font',size=singlesize)
    figex = plt.figure(figsize=(9,7))
    axex = figex.add_subplot(111)
    axex.plot(mobileparts['wx'],mobileparts['wy'],marker='o',linestyle='',color='orange',label="mobile particles")
    axex.plot(immobileparts['wx'],immobileparts['wy'],marker='o',linestyle='',color='green',label="immobile particles")
    axex.plot(undefinedparts['wx'],undefinedparts['wy'],marker='o',linestyle='',color='grey',label="undefined particles")
    axex.set_xlabel(r"$w_x$ [$nm$]")
    axex.set_ylabel(r"$w_x$ [$nm$]")
    axex.set_xlim((0,500))
    axex.set_ylim((0,500))
    axex.legend(loc='best')
    figex.savefig("eccentricityMap-SPOTS.png")
    plt.clf()
    rc('font',size=doublesize)

    '''
    #nf[0] = Date; nf[1] = Number of Images; nf[2] = BaxConcentration
    #nf = norm_factors(conc)

    #print(nf)

    concdnf = [np.array(conc),np.zeros(len(conc)),np.array(conc)]
    for i in range(len(nf[0])):
        for j in range(len(conc)):
            if nf[2][i] == concdnf[2][j]:
                concdnf[1][j] += nf[1][i]
    #print(concdnf)

    conc_and_intens = []
    for i in range(len(nf[0])):
        c0 = len(allparts[np.logical_and(np.logical_and(allparts['date'] == nf[0][i],allparts['Channel'] == 0),allparts['Baxconc'] == nf[2][i])])/nf[1][i]/100
        i0 = allparts[np.logical_and(np.logical_and(allparts['date'] == nf[0][i],allparts['Channel'] == 0),allparts['Baxconc'] == nf[2][i])]['a'].sum()/nf[1][i]/100
        c1 = len(allparts[np.logical_and(np.logical_and(allparts['date'] == nf[0][i],allparts['Channel'] == 1),allparts['Baxconc'] == nf[2][i])])/nf[1][i]/100
        i1 = allparts[np.logical_and(np.logical_and(allparts['date'] == nf[0][i],allparts['Channel'] == 1),allparts['Baxconc'] == nf[2][i])]['a'].sum()/nf[1][i]/100
        conc_and_intens.append([c0,c1,i0,i1])

    cni = np.transpose(conc_and_intens)
    names = []
    for i in range(len(nf[0])):
        names.append(str(nf[0][i]) + " conc: " + str(nf[2][i])+" pM")
    fig = plt.figure(figsize=(13,7))
    ax1 = fig.add_subplot(121)
    ax1.bar(np.arange(len(nf[0]))-0.2,cni[0],color='blue',width=0.4,align='center',label="Bax")
    ax1.bar(np.arange(len(nf[0]))+0.2,cni[1],color='orange',width=0.4,align='center',label="Bid")
    ax1.set_title("Particles / um^2")
    ax1.set_xticks(np.arange(len(nf[0])))
    ax1.set_xticklabels(names,rotation=35,fontdict={'horizontalalignment':'right'})
    ax1.set_ylabel('Particle Conc [um^(-2)]')
    ax1.legend()
    ax2 = fig.add_subplot(122)
    ax2.bar(np.arange(len(nf[0]))-0.2,cni[2],color='blue',width=0.4,align='center',label="Bax")
    ax2.bar(np.arange(len(nf[0]))+0.2,cni[3],color='orange',width=0.4,align='center',label="Bid")
    ax2.set_title("Intensity /um^2")
    ax2.set_xticks(np.arange(len(nf[0])))
    ax2.set_xticklabels(names,rotation=35,fontdict={'horizontalalignment':'right'})
    ax2.set_ylabel('Intensity [kHz/um^2]')
    ax2.legend()
    fig.tight_layout()
    fig.savefig("Daily_Concentration.png",dpi=200)
    plt.clf()


    nf = concdnf
    conc_and_intens = []
    for i in range(len(nf[0])):
        c0 = len(allparts[np.logical_and(allparts['Channel'] == 0,allparts['Baxconc'] == nf[2][i])])/nf[1][i]/100
        i0 = allparts[np.logical_and(allparts['Channel'] == 0,allparts['Baxconc'] == nf[2][i])]['a'].sum()/nf[1][i]/100
        c1 = len(allparts[np.logical_and(allparts['Channel'] == 1,allparts['Baxconc'] == nf[2][i])])/nf[1][i]/100
        i1 = allparts[np.logical_and(allparts['Channel'] == 1,allparts['Baxconc'] == nf[2][i])]['a'].sum()/nf[1][i]/100
        conc_and_intens.append([c0,c1,i0,i1])

    cni = np.transpose(conc_and_intens)
    names = []
    for i in range(len(nf[0])):
        names.append("Conc: " + str(nf[2][i])+" pM")
    fig = plt.figure(figsize=(13,7))
    ax1 = fig.add_subplot(121)
    ax1.bar(np.arange(len(nf[0]))-0.2,cni[0],color='blue',width=0.4,align='center',label="Bax")
    ax1.bar(np.arange(len(nf[0]))+0.2,cni[1],color='orange',width=0.4,align='center',label="Bid")
    ax1.set_title("Particles / um^2")
    ax1.set_xticks(np.arange(len(nf[0])))
    ax1.set_xticklabels(names,rotation=35,fontdict={'horizontalalignment':'right'})
    ax1.set_ylabel('Particle Conc [um^(-2)]')
    ax1.legend()
    ax2 = fig.add_subplot(122)
    ax2.bar(np.arange(len(nf[0]))-0.2,cni[2],color='blue',width=0.4,align='center',label="Bax")
    ax2.bar(np.arange(len(nf[0]))+0.2,cni[3],color='orange',width=0.4,align='center',label="Bid")
    ax2.set_title("Intensity /um^2")
    ax2.set_xticks(np.arange(len(nf[0])))
    ax2.set_xticklabels(names,rotation=35,fontdict={'horizontalalignment':'right'})
    ax2.set_ylabel('Intensity [kHz/um^2]')
    ax2.legend()
    fig.tight_layout()
    fig.savefig("Concentration-byConc.png",dpi=200)
    plt.clf()




    criteria = ['streak or spot?','Channel']
    permuts = [r for r in itertools.product(['spot','streak'],[0,1])]

    #Data sorted by channel and spot/streak
    alldata = []
    #Data output data
    colData = []


    daysort = False
    if daysort == True:
        category = nf[0]
        normcount = nf[1]
        concentration = nf[2]
    else:
        category = np.unique(nf[2])
        normcount = []
        for cat in category:
            sav = 0
            for i in range(len(nf[1])):
                if nf[2][i] == cat:
                    sav += nf[1][i]
            normcount.append(sav)
        concentration = np.array(category)

    
    category=["NOISE"]
    normcount=[1]
    concentration=[1]
    print("Categories are: {:}".format(category))
    print("Normalization Counts are: {:}".format(normcount))
    print("Concentrations for the Categories are: {:}".format(concentration))
    
    corrcoeff = pd.DataFrame(columns=['date','NumImages','Baxconc','streak or spot?','Channel','NumParts','NumCorrParts'],index=np.arange(len(category)*len(permuts)))

    for p in permuts:
        alldata.append(allparts[np.logical_and(allparts[criteria[0]]==p[0],allparts[criteria[1]]==p[1])])
        colData.append(pd.DataFrame(columns=columns,index=np.arange(len(category))))

            
    bcount = 101
    #totalnumexp = []

    histo = np.histogram(allparts[np.logical_and(allparts['streak or spot?']=='spot',allparts["Channel"]==0)]['a'].values,bins=bcount,range=(0,70))
    popt,pcov = fithistogram(histo,"Baxfithistogram",sigmamax=2)
    print(popt)
    baxstart = popt[0]
    baxdist = popt[1]

    histo = np.histogram(allparts[np.logical_and(allparts['streak or spot?']=='spot',allparts["Channel"]==1)]['a'].values,bins=bcount,range=(0,70))
    popt,pcov = fithistogram(histo,"tBidfithistogram",sigmamax=2)
    print(popt)
    bidstart = popt[0]
    biddist = popt[1]

    corrcount = 30
    KDimmob = []
    KDimmob_err = []
    KDmob = []
    KDmob_err = []
    #corrsavefile = open(corr_dir+"corrfile.csv",'w')
    #corrsavefile.write("Cval,KDMobile,KDMobileSTD,KDImmobile,KDImmobileSTD\n")
    for corrcount in range(30,35,5):
        cval= corrcount * 0.01
        print()
        print("******")
        print("Correlation Cutoff: {:}".format(cval))
        sys.stdout.flush()
        

        c_nf = 0
        for elem in category:
            print()
            print("Category: {:}".format(elem))
            normfactor = normcount[c_nf] * 100 # um^2
            
            for i in range(len(permuts)):
                print("Permutation: {:}".format(i))
                if daysort == True:
                    sortcondition = np.logical_and(alldata[i]['date']==elem,alldata[i]['Baxconc']==concentration[c_nf])
                else:
                    sortcondition = alldata[i]['Baxconc']==concentration[c_nf]
                colData[i].loc[c_nf,columns[:34]] = pd.Series(make_stocheo(alldata[i][sortcondition],('Bax{:}-Bax{:}pM' if i % 2 == 0 else 'Bax{:}-Bid{:}').format(elem,concentration[c_nf]),normfactor=1,bincount=bcount,path=('ImmobileStocheo/' if i < 2 else 'MobileStocheo/'),fitall=False,start=(baxstart if i % 2 == 0 else bidstart),dist=(baxdist if i % 2 == 0 else biddist)),index=columns[:34],name=c_nf)
                colData[i].at[c_nf,'Ntotal'] = len(alldata[i][sortcondition])*1.0
                colData[i].at[c_nf,'normfactor'] = normfactor
                colData[i].at[c_nf,'Itotal'] = alldata[i][sortcondition]['a'].sum()
                colData[i].loc[c_nf,'Baxconc'] = concentration[c_nf]
                colData[i].loc[c_nf,('all_correlation','high_correlation','err_high_corr')] = crosscorrelation(alldata[i][sortcondition],('Bax{:}-Bax{:}' if i % 2 == 0 else 'Bax{:}-Bid{:}').format(elem,concentration[c_nf]),cutoff=cval,normfactor=1,bincount=bcount,path=('ImmobileStocheo/' if i < 2 else 'MobileStocheo/'))
                corrcoeff.at[c_nf*len(permuts)+i,['date','NumImages','Baxconc','streak or spot?','Channel']] = [category[c_nf],normcount[c_nf],concentration[c_nf],permuts[i][0],permuts[i][1]]
                corrcoeff.at[c_nf*len(permuts)+i,'NumParts'] = len(alldata[i][np.logical_and(alldata[i]['date']==elem,alldata[i]['Baxconc']==concentration[c_nf])])
                corrcoeff.at[c_nf*len(permuts)+i,'NumCorrParts'] = len(alldata[i][np.logical_and(np.logical_and(alldata[i]['date']==elem,alldata[i]['Baxconc']==concentration[c_nf]),alldata[i]['Correlation Coefficient']>=cval)])
                #print(corrcoeff)

            c_nf += 1

        
        #Plot the total number of experiments for each concentration used.
        corrcoeff.to_csv(corr_dir+"correlation_counts_day-{:}.csv".format(corrcount),sep=',')


        print("Cross-Correlation for cval = {:}".format(corrcount*.01))
        fig2 = plt.figure(figsize=(9,7))
        ax1 = fig2.add_subplot(111)
        ax1.errorbar(colData[0]['Baxconc'],colData[0]['high_correlation']/colData[0]['normfactor'],label="Bax",marker='o',linestyle='',yerr=colData[0]['err_high_corr']/colData[0]['normfactor'])
        ax1.errorbar(colData[1]['Baxconc'],colData[1]['high_correlation']/colData[1]['normfactor'],label="Bid",marker='o',linestyle='',yerr=colData[1]['err_high_corr']/colData[1]['normfactor'])
        ax1.set_xlabel(r"Bax concentration [pM]")
        ax1.set_ylabel(r"Correlated Particles [counts/$\mu$m$^2$]")
        #ax1.set_ylim([0,0.1])
        ax1.set_xlim([0,2100])
        ax1.legend()
        fig2.tight_layout()
        fig2.savefig(corr_dir+"cross-correlation-Immobile-{:}.png".format(corrcount))
        plt.clf()
        fig2 = plt.figure(figsize=(9,7))
        ax2 = fig2.add_subplot(111)
        ax2.errorbar(colData[2]['Baxconc'],colData[2]['high_correlation']/colData[2]['normfactor'],label="Bax",marker='o',linestyle='',yerr=colData[2]['err_high_corr']/colData[2]['normfactor'])
        ax2.errorbar(colData[3]['Baxconc'],colData[3]['high_correlation']/colData[3]['normfactor'],label="Bid",marker='o',linestyle='',yerr=colData[3]['err_high_corr']/colData[3]['normfactor'])
        ax2.set_xlabel(r"Bax concentration [pM]")
        ax2.set_ylabel("Correlated Particles [counts/$\mu$m$^2$]")
        ax2.set_ylim([0,0.014])
        ax2.set_xlim([0,2100])
        ax2.legend()
        fig2.tight_layout()
        fig2.savefig(corr_dir+"cross-correlation-Mobile-{:}.png".format(corrcount))
        plt.clf()

        #print(colData[0][colData[0]['high_correlation']==0]['Baxconc'])

        #============
        #KD value calculations
        yconc = colData[0]['high_correlation']/colData[0]['normfactor']
        xconc1 = (colData[0]['all_correlation']-colData[0]['high_correlation'])/colData[0]['normfactor']
        xconc2 = (colData[1]['all_correlation']-colData[1]['high_correlation'])/colData[1]['normfactor']
        yerr = np.sqrt(colData[0]['high_correlation'].astype('float64'))/colData[0]['normfactor']
        xerr1 = (np.sqrt(colData[0]['all_correlation'].astype('float64'))+np.sqrt(colData[0]['high_correlation'].astype('float64')))/colData[0]['normfactor']
        xerr2 = (np.sqrt(colData[1]['all_correlation'].astype('float64'))+np.sqrt(colData[1]['high_correlation'].astype('float64')))/colData[1]['normfactor']
        xerr = np.sqrt(xerr1.astype('float64')**2 * xconc2.astype('float64')**2 + xconc1.astype('float64')**2 * xerr2.astype('float64')**2)
        yconcm = colData[2]['high_correlation']/colData[2]['normfactor']
        xconc1m = (colData[2]['all_correlation']-colData[2]['high_correlation'])/colData[2]['normfactor']
        xconc2m = (colData[3]['all_correlation']-colData[3]['high_correlation'])/colData[3]['normfactor']
        yerrm = np.sqrt(colData[2]['high_correlation'].astype('float64'))/colData[2]['normfactor']
        xerr1m = (np.sqrt(colData[2]['all_correlation'].astype('float64'))+np.sqrt(colData[2]['high_correlation'].astype('float64')))/colData[2]['normfactor']
        xerr2m = (np.sqrt(colData[3]['all_correlation'].astype('float64'))+np.sqrt(colData[3]['high_correlation'].astype('float64')))/colData[3]['normfactor']
        xerrm = np.sqrt(xerr1m.astype('float64')**2 * xconc2m.astype('float64')**2 + xconc1m.astype('float64')**2 * xerr2m.astype('float64')**2)

        #print(permuts)
        #print(yconc,xconc1,xconc2)
        #print(yconcm,xconc1m,xconc2m)

        KDimmobile = (xconc1*xconc2/yconc).mean()
        KDimmob.append(KDimmobile)
        KDimmobile_err = (xconc1*xconc2/yconc).std()
        KDimmob_err.append(KDimmobile_err)
        KDmobile = (xconc1m*xconc2m/yconcm).mean()
        KDmob.append(KDmobile)
        KDmobile_err = (xconc1m*xconc2m/yconcm).std()
        KDmob_err.append(KDmobile_err)
        #corrsavefile.write("{:},{:},{:},{:},{:}\n".format(cval,KDmobile,KDmobile_err,KDimmobile,KDimmobile_err))
        #============


        #Plot KD Values
        xvals = np.arange(1,100001,100)/100000
        fig3 = plt.figure(figsize=(9,7))
        ax = fig3.add_subplot(111)
        ax.errorbar(xconc1*xconc2,yconc,xerr=xerr,yerr=yerr,marker='o',linestyle='',label="Immobile",color='green')
        ax.errorbar(xconc1m*xconc2m,yconcm,xerr=xerrm,yerr=yerrm,marker='o',linestyle='',label="Mobile",color='orange')
        ax.plot(xvals,xvals/KDimmobile,'-',color='green',label=r'{:4.2f} $\mu$mi$^-2$'.format(KDimmobile))
        ax.plot(xvals,xvals/KDmobile,'-',color='orange',label=r'{:4.2f} $\mu$m$^-2$'.format(KDmobile))
        ax.set_xlabel(r'c[Bax]*c[Bid] in [counts$^2/\mu$m$^4$]')
        ax.set_ylabel(r'c[Bid-Bax] in [counts$/\mu$m$^2$]')
        ax.set_ylim(0.0001,1)
        ax.set_xlim(0.00001,1)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend(loc='best')
        fig3.savefig(corr_dir+"KD-value-{:}.png".format(corrcount))
        plt.clf()
        plt.close('all')

    #corrsavefile.close()




    fig_totalexp = plt.figure(figsize=(8,7))
    axte = fig_totalexp.add_subplot(111)
    ypos = np.arange(len(totalnumexp))
    axte.bar(ypos,totalnumexp,align='center',alpha=0.5,label="# of Experiments")
    xticknames = []
    for elem in conc:
        xticknames.append("{:} pM".format(elem))
    axte.set_xticklabels(xticknames)
    axte.set_xlabel("Bax concentration [pM]")
    axte.set_ylabel("Number of Images")
    fig_totalexp.savefig("TotalNumExperiments.png")

#CONTINUE HERE! colData1 not defined!
    fig_numparts = plt.figure(figsize=(13, 7))
    axnp = fig_numparts.add_subplot(121)
    axnp.errorbar(conc,colData1['Ntotal'],yerr=np.sqrt(colData1['Ntotal'].astype('float64')),label="Immobile Bax")
    axnp.errorbar(conc,colData2['Ntotal'],yerr=np.sqrt(colData2['Ntotal'].astype('float64')),label="Mobile Bax"  )
    axnp.errorbar(conc,colData3['Ntotal'],yerr=np.sqrt(colData3['Ntotal'].astype('float64')),label="Immobile bid")
    axnp.errorbar(conc,colData4['Ntotal'],yerr=np.sqrt(colData4['Ntotal'].astype('float64')),label="Mobile bid"  )
    axnp.set_xlabel("Bax concentration [pM]")
    axnp.set_ylabel("Number of Particles")
    axnp2 = fig_numparts.add_subplot(122)
    axnp2.errorbar(conc,colData1['Ntotal']/colData1['normfactor'],label="Immobile Bax",yerr=np.sqrt((colData1['Ntotal']).astype('float64'))/colData1['normfactor'])
    axnp2.errorbar(conc,colData2['Ntotal']/colData2['normfactor'],label="Mobile Bax"  ,yerr=np.sqrt((colData2['Ntotal']).astype('float64'))/colData2['normfactor'])
    axnp2.errorbar(conc,colData3['Ntotal']/colData3['normfactor'],label="Immobile bid",yerr=np.sqrt((colData3['Ntotal']).astype('float64'))/colData3['normfactor'])
    axnp2.errorbar(conc,colData4['Ntotal']/colData4['normfactor'],label="Mobile bid"  ,yerr=np.sqrt((colData4['Ntotal']).astype('float64'))/colData4['normfactor'])
    axnp2.set_xlabel("Bax concentration [pM]")
    axnp2.set_ylabel("Particles [counts/um^2]")
    axnp.legend()
    axnp2.legend()
    fig_numparts.tight_layout()
    fig_numparts.savefig("TotalParticleNumbers.png")

    fig_Iparts = plt.figure(figsize=(13, 7))
    axIp = fig_Iparts.add_subplot(121)
    axIp.errorbar(conc,colData1['Itotal'],yerr=np.sqrt(colData1['Itotal'].astype('float64')),label="Immobile Bax")
    axIp.errorbar(conc,colData2['Itotal'],yerr=np.sqrt(colData2['Itotal'].astype('float64')),label="Mobile Bax"  )
    axIp.errorbar(conc,colData3['Itotal'],yerr=np.sqrt(colData3['Itotal'].astype('float64')),label="Immobile bid")
    axIp.errorbar(conc,colData4['Itotal'],yerr=np.sqrt(colData4['Itotal'].astype('float64')),label="Mobile bid"  )
    axIp.set_xlabel("Bax concentration [pM]")
    axIp.set_ylabel("Total Intensity [kHz]")
    axIp2 = fig_Iparts.add_subplot(122)
    axIp2.errorbar(conc,colData1['Itotal']/colData1['normfactor'],label="Immobile Bax",yerr=np.sqrt((colData1['Itotal']).astype('float64'))/colData1['normfactor'])
    axIp2.errorbar(conc,colData2['Itotal']/colData2['normfactor'],label="Mobile Bax"  ,yerr=np.sqrt((colData2['Itotal']).astype('float64'))/colData2['normfactor'])
    axIp2.errorbar(conc,colData3['Itotal']/colData3['normfactor'],label="Immobile bid",yerr=np.sqrt((colData3['Itotal']).astype('float64'))/colData3['normfactor'])
    axIp2.errorbar(conc,colData4['Itotal']/colData4['normfactor'],label="Mobile bid"  ,yerr=np.sqrt((colData4['Itotal']).astype('float64'))/colData4['normfactor'])
    axIp2.set_xlabel("Bax concentration [pM]")
    axIp2.set_ylabel("Total Intensity per area [kHz/um^2] ")
    axIp.legend()
    axIp2.legend()
    fig_Iparts.tight_layout()
    fig_Iparts.savefig("TotalIntenity.png")

    print("Intensity Plots")
    #Immobile peak plot
    figs = []
    axs = []
    count = 0
    titles=['Immobile Bax','Immobile bid']
    plotnames = ['immobileBaxPeaks.png','immobileBidPeaks.png']
    for elem in [colData[0],colData[1]]:
        figs.append(plt.figure(figsize=(9,7)))
        axs.append(figs[-1].add_subplot(111))
        for dat in ['A1','A2','A3','A4']:
            axs[-1].errorbar(elem['Baxconc'],elem[dat]/elem['normfactor'],marker='o',linestyle='',yerr=elem['err_'+dat]/elem['normfactor'],label=dat)
        axs[-1].set_xlabel("Bax concentration [pM]")
        axs[-1].set_ylabel("Particles [counts/um^2]")
        axs[-1].set_xlim((0,2100))
        axs[-1].set_ylim((0,0.003))
        axs[-1].legend(loc='best')
        figs[-1].tight_layout()
        figs[-1].savefig("fit-"+plotnames[count])
        figs.append(plt.figure(figsize=(9,7)))
        axs.append(figs[-1].add_subplot(111))
        labels = ['monomer','dimer','trimer','tetramer','pentamer','hexamer','heptamer']
        lcount = 0
        for dat in ['Binsum1','Binsum2','Binsum3','Binsum4','Binsum5','Binsum6','Binsum7']:
            axs[-1].errorbar(elem['Baxconc'],elem[dat]/elem['normfactor'],marker='o',linestyle='',yerr=elem['err_'+dat]/elem['normfactor'],label=labels[lcount])
            lcount += 1
            
        axs[-1].set_xlabel("Bax concentration [pM]")
        axs[-1].set_ylabel(r"Particles [counts/$\mu$m$^2$]")
        axs[-1].set_xlim((0,2100))
        axs[-1].set_ylim((0,0.013))
        axs[-1].legend(loc='best')
        figs[-1].tight_layout()
        figs[-1].savefig("abs-"+plotnames[count])
        figs.append(plt.figure(figsize=(9,7)))
        axs.append(figs[-1].add_subplot(111))
        sumbins = elem[['Binsum1','Binsum2','Binsum3','Binsum4','Binsum5','Binsum6']]
        sumbins = sumbins.sum(axis=1)
        lcount = 0
        for dat in ['Binsum1','Binsum2','Binsum3','Binsum4','Binsum5','Binsum6','Binsum7']:
            axs[-1].errorbar(elem['Baxconc'],elem[dat]/sumbins,marker='o',linestyle=':',yerr=elem['err_'+dat]/sumbins,label=labels[lcount])
            lcount += 1
        axs[-1].set_xlabel("Bax concentration [pM]")
        axs[-1].set_ylabel(r"particle percentage")
        axs[-1].set_xlim((0,2100))
        axs[-1].set_ylim((0,1.))
        axs[-1].legend(loc='best')
        figs[-1].tight_layout()
        figs[-1].savefig("rel-"+plotnames[count])
        count += 1
    plt.clf()


    #Fit offset and Monomer intensity for Bax and Bid
    fig1 = plt.figure(figsize=(13,7))
    fig1.suptitle("Intensities")
    ax1 = fig1.add_subplot(221)
    ax1.errorbar(colData[0]['Baxconc'],colData[0]['Ioffset'],label="Bax",marker='o',linestyle='',yerr=colData[0]['err_Ioffset'])
    ax1.errorbar(colData[1]['Baxconc'],colData[1]['Ioffset'],label="bid",marker='o',linestyle='',yerr=colData[1]['err_Ioffset'])
    ax1.set_title("Offset Immobile")
    ax1.set_ylim([0,10])
    ax1.set_xlabel("Bax conc [counts/um^2]")
    ax1.set_ylabel("Intensity [kHz]")
    ax1.legend(loc='best')
    ax2 = fig1.add_subplot(223)
    ax2.errorbar(colData[0]['Baxconc'],colData[0]['Imono'],label="Bax",marker='o',linestyle='',yerr=colData[0]['err_Imono'])
    ax2.errorbar(colData[1]['Baxconc'],colData[1]['Imono'],label="bid",marker='o',linestyle='',yerr=colData[1]['err_Imono'])
    ax2.set_title("Monomer Immobile")
    ax2.set_ylim([0,10])
    ax2.set_xlabel("Bax conc [counts/um^2]")
    ax2.set_ylabel("Intensity [kHz]")
    ax2.legend(loc='best')
    ax3 = fig1.add_subplot(222)
    ax3.errorbar(colData[2]['Baxconc'],colData[2]['Ioffset'],label="Bax",marker='o',linestyle='',yerr=colData[2]['err_Ioffset'])
    ax3.errorbar(colData[3]['Baxconc'],colData[3]['Ioffset'],label="bid",marker='o',linestyle='',yerr=colData[3]['err_Ioffset'])
    ax3.set_title("Offset Mobile")
    ax3.set_ylim([0,10])
    ax3.set_xlabel("Bax conc [counts/um^2]")
    ax3.set_ylabel("Intensity [kHz]")
    ax3.legend(loc='best')
    ax4 = fig1.add_subplot(224)
    ax4.errorbar(colData[2]['Baxconc'],colData[2]['Imono'],label="Bax",marker='o',linestyle='',yerr=colData[2]['err_Imono'])
    ax4.errorbar(colData[3]['Baxconc'],colData[3]['Imono'],label="bid",marker='o',linestyle='',yerr=colData[3]['err_Imono'])
    ax4.set_title("Monomer Mobile")
    ax4.set_ylim([0,10])
    ax4.set_xlabel("Bax conc [counts/um^2]")
    ax4.set_ylabel("Intensity [kHz]")
    ax4.legend(loc='best')
    fig1.tight_layout()
    fig1.savefig("Intensitychange-Counts.png")
    plt.clf()
    '''




if __name__=="__main__":
    main()
