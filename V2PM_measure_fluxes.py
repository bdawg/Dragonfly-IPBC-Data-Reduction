# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 15:36:15 2015

@author: bnorris

This takes the fits files from the Xenics camera and puts the measured waveguide
fluxes into appropriate variables.
"""


import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import pdb
from scipy.io.idl import readsav
from circles import circles
from numpy.polynomial.polynomial import polyval


global p1Orp2, nWgsOP, nWgsIN, nBLs, nPists, segmentOrder, showPlots, nonLinCoeffs
p1Orp2 = 1 #Which side of the baseline to piston? (1 or 2)
nWgsOP = 120
nWgsIN = 8
nBLs = 28
nPists = 61
wavelength = 1.55 #microns
showPlots = True
saveFigs = False

nonLinCoeffsFile = 'coeffs_detailed_6.idlvar'

segmentOrder = np.array([23, 25, 27, 9, 29, 36, 31, 33])

outFile = 'extractedFluxes_ipbc_20140930_3_xtr2-nonlin6'

############################################################################################################################################################
# LAB DATA
############################################################################################################################################################

#dataDir='./dragonfly20140930data/'
dataDir='/Volumes/mojito2/snert/dragonfly/labtests/20140930/data_fits/'

darkFrame="ipbc_20140930_3_dark.mat.fits"

singleChannels=['ipbc_20140930_3_SingleWGScan_WG_23.mat.fits','ipbc_20140930_3_SingleWGScan_WG_25.mat.fits',
       'ipbc_20140930_3_SingleWGScan_WG_27.mat.fits','ipbc_20140930_3_SingleWGScan_WG_9.mat.fits',
       'ipbc_20140930_3_SingleWGScan_WG_29.mat.fits','ipbc_20140930_3_SingleWGScan_WG_36.mat.fits',
       'ipbc_20140930_3_SingleWGScan_WG_31.mat.fits','ipbc_20140930_3_SingleWGScan_WG_33.mat.fits']


baselineChannels1=[
'ipbc_20140930_3_BaselineScan_WGs_23_25_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_27_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_23_9_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_29_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_23_36_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_31_p1.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_23_33_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_27_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_25_9_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_29_p1.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_25_36_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_31_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_25_33_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_9_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_27_29_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_36_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_27_31_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_33_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_9_29_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_9_36_p1.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_9_31_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_9_33_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_29_36_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_29_31_p1.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_29_33_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_36_31_p1.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_36_33_p1.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_31_33_p1.mat.fits']



baselineChannels2=[
'ipbc_20140930_3_BaselineScan_WGs_23_25_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_27_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_23_9_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_29_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_23_36_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_23_31_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_23_33_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_27_p2.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_25_9_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_29_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_25_36_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_25_31_p2.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_25_33_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_9_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_27_29_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_36_p2.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_27_31_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_27_33_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_9_29_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_9_36_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_9_31_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_9_33_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_29_36_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_29_31_p2.mat.fits', 
'ipbc_20140930_3_BaselineScan_WGs_29_33_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_36_31_p2.mat.fits',
'ipbc_20140930_3_BaselineScan_WGs_36_33_p2.mat.fits', 'ipbc_20140930_3_BaselineScan_WGs_31_33_p2.mat.fits']


allOnScans=["ipbc_20140930_3b_AllOnScan_WG_23.mat.fits","ipbc_20140930_3b_AllOnScan_WG_25.mat.fits",
"ipbc_20140930_3b_AllOnScan_WG_27.mat.fits","ipbc_20140930_3b_AllOnScan_WG_9.mat.fits",
"ipbc_20140930_3b_AllOnScan_WG_29.mat.fits","ipbc_20140930_3b_AllOnScan_WG_36.mat.fits",
"ipbc_20140930_3b_AllOnScan_WG_31.mat.fits","ipbc_20140930_3b_AllOnScan_WG_33.mat.fits"
]




#==============================================================================
# Dark Frame, etc.
#==============================================================================
dark=pf.open(dataDir+darkFrame)
dark = dark['PRIMARY'].data
dark=np.average(dark,axis=0)
nonLinCoeffs = readsav(nonLinCoeffsFile).coeffs

if showPlots:
    plt.figure(1)
    plt.pause(5)

#==============================================================================
# Function to extract fluxes from a given frame
#==============================================================================
# TODO - Rigth now you have to manually put in the following numbers. Make
# an interactive GUI for this.
def fluxesFromIm(image, xstart=86., xend=540., ycoord=234., radius=1.5):
    ycents = np.concatenate((np.repeat(ycoord+1,50),np.repeat(ycoord,70)))
    xcents = np.linspace(xstart, xend, nWgsOP)
    xpx=range(len(image[0,:]))    
    ypx=range(len(image[:,0]))
    xinds,yinds=np.meshgrid(xpx,ypx)
    #pdb.set_trace()

    image = np.squeeze(polyval(image, nonLinCoeffs)) #Non-linear correction
    
    if False: #showPlots: <- This is too slow.
    #if showPlots: #<- This is too slow.
        plt.figure(2)
        #plt.clf()
        plt.imshow(image, interpolation='nearest')
        plt.plot(xcents,ycents,'+r',linewidth=1)
        circles(xcents, ycents, radius, facecolor='none',edgecolor='r')
        plt.pause(0.01)
        
    fluxes = np.zeros(nWgsOP)
    for wg in range(nWgsOP):  
        circleInds = np.sqrt( (xinds-xcents[wg])**2 + (yinds-ycents[wg])**2 ) < radius
        fluxes[wg]=np.average(image[circleInds])

    return fluxes




#==============================================================================
# Extract Single Waveguide Scans
#==============================================================================
wgInInds = range(nWgsIN)
singleWGscan_allChipOuts = np.zeros((nPists,nWgsOP,nWgsIN))
print '######## Doing Single Waveguide Scans ########'

for i in wgInInds:
    wg2pist = i
    print 'Pistoning waveguide %d' % i
    filename = dataDir+singleChannels[i]
    hdulist=pf.open(filename)    
    imCube=hdulist['PRIMARY'].data

    for pist in range(nPists):
        print 'Piston step %d of %d' % (pist+1, nPists)
        im=imCube[pist,:,:] - dark
        fluxes = fluxesFromIm(im)
        singleWGscan_allChipOuts[pist,:,i] = fluxes

    if showPlots:
        plt.figure(1)
        plt.clf()
        im2show = singleWGscan_allChipOuts[:,:,i].T
        plt.imshow(im2show, interpolation='nearest')
        plt.colorbar(orientation='vertical')
        plt.pause(0.01)        
        
        if saveFigs:
            fname="fig_ExtractSingleWGScan_%d.pdf" % i       
            plt.savefig(fname)

singleWGscan_allChipOuts_avgd = np.average(singleWGscan_allChipOuts, axis=0)





#==============================================================================
# Extract Baseline Scans
#==============================================================================
bl2wg = np.array([  [23,25 ],
                    [23,27],
                    [23,9 ],
                    [23,29],
                    [23,36 ],
                    [23,31],
                    [23,33 ],
                    [25,27],
                    [25,9 ],
                    [25,29],
                    [25,36 ],
                    [25,31],
                    [25,33 ],
                    [27,9],
                    [27,29 ],
                    [27,36],
                    [27,31 ],
                    [27,33],
                    [9,29 ],
                    [9,36],
                    [9,31 ],
                    [9,33],
                    [29,36 ],
                    [29,31],
                    [29,33 ],
                    [36,31],
                    [36,33 ],
                    [31,33] ])

blInds = range(nBLs)
BLscan_allChipOuts = np.zeros((nPists,nWgsOP,nBLs))
BLscan_pistVals = np.zeros((nPists,nBLs))
print '######## Doing Baseline Scans ########'

for i in blInds:
    wg2pist = bl2wg[i, p1Orp2-1]
    print 'Baseline %d, pistoning waveguide %d' % (i, wg2pist)
    if p1Orp2 == 1:
        filename = dataDir+baselineChannels1[i]
    else:
        filename = dataDir+baselineChannels2[i]
    hdulist=pf.open(filename)    
    imCube=hdulist['PRIMARY'].data 
    
    # Put pistvals in radians. NB *2 is because OPL = 2*MEMS_piston    
    BLscan_pistVals[:,i] = (imCube[:,0,0]*2. / wavelength) * 2.*np.pi
    print BLscan_pistVals[:,i]
    
    for pist in range(nPists):
        print 'Piston step %d of %d' % (pist+1, nPists)
        im=imCube[pist,:,:] - dark
        fluxes = fluxesFromIm(im)
        BLscan_allChipOuts[pist,:,i] = fluxes
        
    if showPlots:
        plt.figure(1)
        plt.clf()
        im2show = BLscan_allChipOuts[:,:,i].T
        plt.imshow(im2show, interpolation='nearest')
        plt.colorbar(orientation='vertical')
        plt.pause(0.01)     

        if saveFigs:
            fname="fig_ExtractBLScan_%d.pdf" % i       
            plt.savefig(fname)


#==============================================================================
# Extract AllOn Scans (for later testing)
#==============================================================================

wgInInds = range(nWgsIN)
allOnWGscan_allChipOuts = np.zeros((nPists,nWgsOP,nWgsIN))
print '######## Doing AllOn Scans ########'

for i in wgInInds:
    wg2pist = i
    print 'Pistoning waveguide %d' % i
    filename = dataDir+allOnScans[i]
    hdulist=pf.open(filename)    
    imCube=hdulist['PRIMARY'].data
        
    for pist in range(nPists):
        print 'Piston step %d of %d' % (pist+1, nPists)
        im=imCube[pist,:,:] - dark
        fluxes = fluxesFromIm(im)
        allOnWGscan_allChipOuts[pist,:,i] = fluxes
        
    if showPlots:
        plt.figure(1)
        plt.clf()
        im2show = allOnWGscan_allChipOuts[:,:,i].T
        plt.imshow(im2show, interpolation='nearest')
        plt.colorbar(orientation='vertical')
        plt.pause(0.01)        
        
        if saveFigs:
            fname="fig_ExtractAllOnScan_%d.pdf" % i       
            plt.savefig(fname)



#==============================================================================
# Save everything
#==============================================================================

#FIXME
np.savez(outFile, singleWGscan_allChipOuts_avgd=singleWGscan_allChipOuts_avgd, 
                  singleWGscan_allChipOuts=singleWGscan_allChipOuts, 
                  BLscan_allChipOuts=BLscan_allChipOuts,
                  BLscan_pistVals=BLscan_pistVals,
                  allOnWGscan_allChipOuts=allOnWGscan_allChipOuts)


#==============================================================================
# To restore, you'd do:
#file='extractedFluxes_ipbc_20140930_3_xtr1.npz'
#npzfile=np.load(file)
#singleWGscan_allChipOuts=npzfile['singleWGscan_allChipOuts']
#etc.
#npzfile.files will tell you what variables are in there (by name)
#==============================================================================










