# -*- coding: utf-8 -*-
"""
Make a V2PM and P2VM from fluxes in a npz file.
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import curve_fit
from time import sleep

nWgsOP = 120
nWgsIN = 8
nBLs = 28
nPists = 61
showPlots = True

segmentOrder = np.array([23, 25, 27, 9, 29, 36, 31, 33])

#inFile = 'extractedFluxes_ipbc_20140930_3_xtr1'
inFile = 'extractedFluxes_ipbc_20140930_3_xtr2-nonlin6'


fitPhases = True # If false, restore phases from phaseFile. 
                  # If true, fit phases and save to phaseFile.
phaseFile = 'phases_ipbc_20140930_3_xtr2-nonlin6'

analyseTestData = True #Use P2VM to analyse data in allOnWGscan_allChipOuts

#==============================================================================
# Restore flux values from file
#==============================================================================
npzfile=np.load(inFile+'.npz')
singleWGscan_allChipOuts_avgd = npzfile['singleWGscan_allChipOuts_avgd']
allOnWGscan_allChipOuts = npzfile['allOnWGscan_allChipOuts']
BLscan_allChipOuts = npzfile['BLscan_allChipOuts']
BLscan_pistVals=npzfile['BLscan_pistVals']



#==============================================================================
# Extract the phases and amplitudes from the baseline scans...
# Note: 'pistons' are always measured in radians. They have been converted from
# microns in V2PM_measure_fluxes.py
#==============================================================================
#blInds = [0]#range(nBLs)
blInds = range(nBLs)
wgInds = range(nWgsIN)
opInds = range(nWgsOP)
BLphases = np.zeros((nWgsOP,nBLs))
BLamps = np.zeros((nWgsOP,nBLs))


def sinFunc(x, sinParam0, sinParam1 ,sinParam2 ,sinParam3):
    # sinParams = [ampl, freq, phase, offset]
    y = sinParam0 * np.sin(sinParam1*x + sinParam2) + sinParam3
    return y
    
    
if fitPhases:    
    count=1
    plotlim = 1 #Plot OPs with amplitudes > this.
    
    for bl in blInds:
        pistAmts = BLscan_pistVals[:,bl]
        print "Fitting baseline %d" % bl
        if showPlots:
            plt.pause(0.001)
            plt.figure(5)  
            plt.clf()
     
        for op in opInds:
            current = BLscan_allChipOuts[:,op,bl]
            # Estimate starting params:
            offsetGuess = np.mean(current)
            amplGuess = np.std(current)*2./np.sqrt(2) #Check this is right...
            p0 = [amplGuess, 0.9, np.pi, offsetGuess]
            fitp, fitpcov = curve_fit(sinFunc, pistAmts, current, p0)
            BLphases[op,bl] = fitp[2]
            BLamps[op,bl] = fitp[0]
    
            if showPlots:
                plt.subplot(11,11,op+1)
                plt.plot(pistAmts,current,'+')
                plt.plot(pistAmts, sinFunc(pistAmts,fitp[0],fitp[1],fitp[2],fitp[3]), linewidth=1)
                #plt.ylim(0,10)
                
    #            if np.abs(fitp[0]) > plotlim:
    #                print fitp
    #                #print pistAmts
    #                plt.subplot(6,6,count)
    #                plt.plot(pistAmts,current,'+')
    #                plt.plot(pistAmts, sinFunc(pistAmts,fitp[0],fitp[1],fitp[2],fitp[3]), linewidth=2)
    #                #plt.ylim(0,10)
    #                count = count + 1
    #                plt.pause(0.001)
    #            
                    
    # Make amplitudes positive
    BLphases[BLamps < 0] = BLphases[BLamps < 0] + np.pi
    BLamps[BLamps < 0] = -1.*BLamps[BLamps < 0]
    
    # A -pi/2 offset to everything makes the mean closure phase zero (i.e. it
    # removes the offset.) But its still quite the same as the artifical v2pm. 
    # But thus doesn't affect CP stability, just offset. 
    # BLphases = BLphases - np.pi/2

    np.savez(phaseFile, BLphases=BLphases, BLamps=BLamps)

else:
    npzPfile=np.load(phaseFile+'.npz')
    BLphases = npzPfile['BLphases']
    BLamps = npzPfile['BLamps']
    
BLcomplx = BLamps * np.exp(1j*BLphases)


#==============================================================================
# Assemble V2PM and P2VM
#==============================================================================

V2PM = np.complex128(np.zeros((nWgsOP,nWgsIN+nBLs)))
V2PM[:,0:nWgsIN] = singleWGscan_allChipOuts_avgd
V2PM[:,nWgsIN:nWgsIN+nBLs] = BLcomplx
V2PMabs = np.abs(V2PM)
V2PMphi = np.angle(V2PM)
V2PMphi[V2PMabs < 50] = 0. #Just for nice plotting

if showPlots:
    plt.figure(3)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('V2PM abs')
    plt.imshow(V2PMabs,interpolation='nearest')
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.title('V2PM phi')
    plt.imshow(V2PMphi,interpolation='nearest')
    plt.colorbar()

P2VM=np.linalg.pinv(V2PM)
newP2VM = P2VM
P2VMabs = np.abs(P2VM)
P2VMphi = np.angle(P2VM)
#P2VMphi[P2VMabs < 50] = 0. #Just for nice plotting

if showPlots:
    plt.figure(4)
    plt.clf()
    plt.subplot(1,2,1)
    plt.title('P2VM abs')
    plt.imshow(P2VMabs,interpolation='nearest')
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.title('P2VM phi')
    plt.imshow(P2VMphi,interpolation='nearest')
    plt.colorbar()



#==============================================================================
# Analyse test data
#==============================================================================
if analyseTestData:
    pistonedWG = 1 # Use data set with which waveguide pistoned?
    
    cVis = np.complex128(np.zeros((nWgsIN+nBLs,nPists)))
    for i in range(nPists):
        cVis[:,i] = np.dot(P2VM, allOnWGscan_allChipOuts[i,:,pistonedWG])
        
    V2s = np.abs(cVis)**2
    phis = np.angle(cVis)
    
    if showPlots:
        plt.figure(5)
        plt.clf()
        plt.subplot(2,1,1)
        plt.title('Recovered $V^{2}$s')
        plt.imshow(V2s,interpolation='nearest')#,cmap='hot')
        plt.colorbar(orientation='vertical')
        plt.subplot(2,1,2)
        plt.title('Recovered phases')
        plt.imshow(phis,interpolation='nearest')#,cmap='hot')
        plt.colorbar(orientation='vertical')
    
    
    
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
 

    #==============================================================================
    # Dodgy manual CP testing...
    #==============================================================================
    chosenBLs = [0,11,5]
      
    p1 = phis[nWgsIN+chosenBLs[0]]
    p2 = phis[nWgsIN+chosenBLs[1]]
    p3 = -phis[nWgsIN+chosenBLs[2]]
    
#    p1 = np.repeat(0.,61)
#    p2 = phis[18]
#    p3 = -phis[20]
     
    testCPs = p1 + p2 + p3


    cv1 = cVis[nWgsIN+chosenBLs[0]]
    cv2 = cVis[nWgsIN+chosenBLs[1]]
    cv3 = cVis[nWgsIN+chosenBLs[2]]
    
    cv1=np.repeat(1j,61)
    cv1 = cVis[23]
    cv2 = cVis[18]
    cv3 = cVis[20]
    cv1 = np.roll(cv1,1)
    cv3 = np.roll(cv3,1)
     
    testBS = cv1 * cv2 * np.conjugate(cv3)
    testCPs = np.angle(testBS)
    p1 = np.angle(cv1)
    p2 = np.angle(cv2)
    p3 = np.angle(cv3)
    
    plt.figure(6)
    plt.clf()
    plt.subplot(4,1,1)
    plt.plot(pistAmts,p1)
    plt.subplot(4,1,2)
    plt.plot(pistAmts,p2)
    plt.subplot(4,1,3)
    plt.plot(pistAmts,p3)
    plt.subplot(4,1,4)
    plt.plot(pistAmts,testCPs)
    
    
    
    
    
    
    
    
    

