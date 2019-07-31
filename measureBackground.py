# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 14:47:10 2019

@author: Markus
"""

#Measure Background

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skimage import io



def giveBackground(image):
    image = image.flatten()
    print("Image Mean:   {:}".format(image.mean()))
    print("Image STD:    {:}".format(image.std()))
    print("Image Median: {:}".format(np.median(image)))
    print("Image Max:    {:}".format(image.max()))
    print("Image Min:    {:}".format(image.min()))
    print()
    return


def plotImages(img1,img2):
    
    fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    counter = 0
    for img in (img1,img2):
        axes[0,counter].imshow(img)
        axes[1,counter].boxplot(img.flatten(),patch_artist=True,widths=0.3)
        axes[1,counter].set_xlim([0,2])
        counter += 1
    df = pd.DataFrame(image.flatten(),columns=["pixVal"])
    df['pix_quantiles'] = pd.qcut(df['pixVal'],4,duplicates='drop')
    return

def main():
    image1 = io.imread("single_TIFF.tif")
    image2 = io.imread("background_exp.tif")
    
    print("Simulation")
    giveBackground(image1)
    print("Experiment")
    giveBackground(image2)
    
    plotImages(image1,image2)
    return
    
    
if __name__=="__main__":
    main()
    

    
