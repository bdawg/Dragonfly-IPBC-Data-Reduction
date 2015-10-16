3# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 12:13:06 2015

@author: bnorris

Analyse the fitted phase data and calciulate closure phases.

"""

import numpy as np
import matplotlib.pyplot as plt
import pdb

inputFittedPhaseFile = 'extractedPhases_amplGT1e-9_ipbc_20140930_3_xtr2-nonlin6_2'
nWgsOP = 120
nWgsIN = 8
nBLs = 28
nPists = 61
nCPs = 56
showPlots = True

#np.savez(outputFittedPhaseFile, BLphaseABCD = BLphaseABCD, 
#                                BLphaseABCDchisqr = BLphaseABCDchisqr, 
#                                BLintensABCD = BLintensABCD, 
#                                BLphases = BLphases, 
#                                BLamps = BLamps, 
#                                BL2wgsOP = BL2wgsOP, 
#                                wgsIN2WgsOP = wgsIN2WgsOP, 
#                                pistAmts = pistAmts, 
#                                allOnPhases = allOnPhases, 
#                                allOnAmps = allOnAmps, 
#                                relBLph = relBLph )
                                
                                
inFile=np.load(inputFittedPhaseFile+'.npz')
BLphaseABCD = inFile['BLphaseABCD']
BLphaseABCDchisqr = inFile['BLphaseABCDchisqr']
wgsIN2WgsOP = inFile['wgsIN2WgsOP']
BL2wgsOP = inFile['BL2wgsOP']
pistAmts = inFile['pistAmts']


if showPlots:
    plt.figure(1)
    plt.clf()
    for i in range(nWgsIN):
        plt.subplot(4, 4, 2*i+1)    
        plt.imshow(BLphaseABCD[:,:,i],interpolation='nearest')
        plt.title('Pistoning waveguide %d' % i)
        #plt.colorbar()
        plt.subplot(4, 4, 2*i+2)  
        plt.imshow(BLphaseABCDchisqr[:,:,i],interpolation='nearest')


#==============================================================================
# Find closure phase triangles
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
                    
segmentOrder = np.array([23, 25, 27, 9, 29, 36, 31, 33])

for i in range(nWgsIN):
    bl2wg[bl2wg == segmentOrder[i]] = i


triangles = np.zeros((nCPs,3))
count=0
for i in range(0, nWgsIN-2):
    for j in range(i+1, nWgsIN-1):
        for k in range(j+1, nWgsIN):
            w=np.asscalar(np.where((bl2wg == [i,j]).all(axis=1))[0])  
            triangles[count, 0] = w
            w=np.asscalar(np.where((bl2wg == [j,k]).all(axis=1))[0])
            triangles[count, 1] = w
            w=np.asscalar(np.where((bl2wg == [i,k]).all(axis=1))[0])
            triangles[count, 2] = w
            
            print triangles[count,:]                        
            count += 1




#==============================================================================
# Calculate closure phases
#==============================================================================
# BLphaseABCD = np.zeros((nBLs, nPists, nWgsIN))

chisqCut = 0.1 #Mask out chisqs above this
allCPs = np.zeros((nCPs, nPists, nWgsIN))
allCPMask = np.zeros((nCPs, nPists, nWgsIN),dtype=np.bool)

for scannedWG in range(0,nWgsIN):
    for cpInd in range(0,nCPs):
        for pistInd in range(nPists):
            i0 = triangles[cpInd, 0]
            i1 = triangles[cpInd, 1]
            i2 = triangles[cpInd, 2]
            
#            allCPs[cpInd, pistInd, scannedWG] = BLphaseABCD[i0,pistInd, scannedWG] \
#                + BLphaseABCD[i1,pistInd, scannedWG] - BLphaseABCD[i2,pistInd, scannedWG]
#            if allCPs[cpInd, pistInd, scannedWG] < 0:
#                allCPs[cpInd, pistInd, scannedWG] = allCPs[cpInd, pistInd, scannedWG] + 2.*np.pi

            # Do as exps to avoid wrapping problems
            tripProd = np.exp(1j*BLphaseABCD[i0,pistInd, scannedWG]) * \
                        np.exp(1j*BLphaseABCD[i1,pistInd, scannedWG]) * \
                        np.conj(np.exp(1j*BLphaseABCD[i2,pistInd, scannedWG]))
            allCPs[cpInd, pistInd, scannedWG]=np.angle(tripProd)


            if (BLphaseABCDchisqr[i0,pistInd, scannedWG] > chisqCut) or \
                (BLphaseABCDchisqr[i1,pistInd, scannedWG] > chisqCut) or \
                (BLphaseABCDchisqr[i2,pistInd, scannedWG] > chisqCut):
                   allCPMask[cpInd, pistInd, scannedWG] = True

#        if showPlots:
#            plt.figure(2)
#            plt.clf()
#            plt.plot(pistAmts, BLphaseABCD[i0,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, BLphaseABCD[i1,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, BLphaseABCD[i2,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, allCPs[cpInd, :, scannedWG])
#            print np.std(allCPs[cpInd, :, scannedWG])/np.pi*180.
#            plt.pause(0.5)


allCPsMasked = np.ma.MaskedArray(allCPs, mask = allCPMask)




#==============================================================================
# Plot CPs
#==============================================================================

## Show ALL the CPs
#if showPlots:
#    for scannedWG in range(0,nWgsIN):
#        plt.figure(2)
#        plt.clf()
#        for cpInd in range(0,nCPs):
#            i0 = triangles[cpInd, 0]
#            i1 = triangles[cpInd, 1]
#            i2 = triangles[cpInd, 2]
#
#            plt.subplot(7,8,cpInd+1)
#            plt.plot(pistAmts, BLphaseABCD[i0,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, BLphaseABCD[i1,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, BLphaseABCD[i2,:, scannedWG], linestyle='--')
#            plt.plot(pistAmts, allCPsMasked[cpInd, :, scannedWG] - np.mean(allCPsMasked[cpInd, :, scannedWG]))
#            sd = np.std(allCPsMasked[cpInd, :, scannedWG])/np.pi*180.
#            print sd
#            plt.title('SD (deg): %f' % sd)
#            #plt.pause(0.001)
#        plt.pause(1)
        
        
# Just show pistoned triangles
if showPlots:
    saveFigs = False    
    savePref='iobc_fringefitting_CPplots_'
    
    ### nPistonedTriangles = 56###21
    nPistonedTriangles = 21
    allSDs = np.zeros((nPistonedTriangles, nWgsIN))

    for scannedWG in range(0,nWgsIN):
        plt.figure(2)
        plt.clf()
        count=0
        for cpInd in range(0,nCPs):
                        
            i0 = triangles[cpInd, 0]
            i1 = triangles[cpInd, 1]
            i2 = triangles[cpInd, 2]

            curWgsIN = [bl2wg[i0,:], bl2wg[i1,:], bl2wg[i2,:]]
            curWgsIN = np.asarray(curWgsIN)
            if np.sum(curWgsIN == scannedWG):
            ### if True:#np.sum(curWgsIN == scannedWG):
                plt.subplot(6,4,count+1)
                ### plt.subplot(7,8,count+1)
                plt.plot(pistAmts, BLphaseABCD[i0,:, scannedWG], linestyle='--')
                plt.plot(pistAmts, BLphaseABCD[i1,:, scannedWG], linestyle='--')
                plt.plot(pistAmts, BLphaseABCD[i2,:, scannedWG], linestyle='--')
                
                nmCP = allCPsMasked[cpInd, :, scannedWG] - np.mean(allCPsMasked[cpInd, :, scannedWG])
                if np.sum((nmCP > np.pi) + (nmCP < -np.pi)):
                    nmCP[nmCP > np.pi] = nmCP[nmCP > np.pi]-2.*np.pi
                    nmCP[nmCP < -np.pi] = nmCP[nmCP < -np.pi]+2.*np.pi
                    nmCP = nmCP - np.mean(nmCP)
                plt.plot(pistAmts, nmCP, linewidth=1, color='k')
                
                #sd = np.std(allCPsMasked[cpInd, :, scannedWG])/np.pi*180.
                sd = np.std(nmCP)/np.pi*180.
                allSDs[count, scannedWG] = sd
                #print sd
                plt.title('SD (deg): %f' % sd)
                #plt.pause(0.001)
                count += 1 
        plt.suptitle('Pistoning waveguide %d' % scannedWG, verticalalignment='bottom')
        plt.pause(.01)
        
        if saveFigs:
            fname=savePref+"%d.png" % scannedWG       
            plt.savefig(fname)
        
    plt.figure(3)
    plt.clf()
    plt.hist(np.reshape(allSDs, -1))
    plt.xlabel('Closure phase SD (degrees)')
    plt.ylabel('Frequency (across all triangles and piston choices')
        
  
        
        


#==============================================================================
# This will plot *all* triangles...
#==============================================================================
#if showPlots:
#    saveFigs = False    
#    savePref='iobc_fringefitting_CPplots_'
#    
#    ### nPistonedTriangles = 56###21
#    nPistonedTriangles = 56
#    allSDs = np.zeros((nPistonedTriangles, nWgsIN))
#
#    for scannedWG in [5]:#range(0,nWgsIN):
#        plt.figure(12)
#        plt.clf()
#        count=0
#        for cpInd in range(0,nCPs):
#                        
#            i0 = triangles[cpInd, 0]
#            i1 = triangles[cpInd, 1]
#            i2 = triangles[cpInd, 2]
#
#            curWgsIN = [bl2wg[i0,:], bl2wg[i1,:], bl2wg[i2,:]]
#            curWgsIN = np.asarray(curWgsIN)
#            if True:#np.sum(curWgsIN == scannedWG):
#                plt.subplot(7,8,count+1)
#                ### plt.subplot(7,8,count+1)
#                plt.plot(pistAmts, BLphaseABCD[i0,:, scannedWG], linestyle='--')
#                plt.plot(pistAmts, BLphaseABCD[i1,:, scannedWG], linestyle='--')
#                plt.plot(pistAmts, BLphaseABCD[i2,:, scannedWG], linestyle='--')
#                
#                nmCP = allCPsMasked[cpInd, :, scannedWG] - np.mean(allCPsMasked[cpInd, :, scannedWG])
#                if np.sum((nmCP > np.pi) + (nmCP < -np.pi)):
#                    nmCP[nmCP > np.pi] = nmCP[nmCP > np.pi]-2.*np.pi
#                    nmCP[nmCP < -np.pi] = nmCP[nmCP < -np.pi]+2.*np.pi
#                    nmCP = nmCP - np.mean(nmCP)
#                plt.plot(pistAmts, nmCP, linewidth=1, color='k')
#                
#                #sd = np.std(allCPsMasked[cpInd, :, scannedWG])/np.pi*180.
#                sd = np.std(nmCP)/np.pi*180.
#                allSDs[count, scannedWG] = sd
#                #print sd
#                plt.title('SD (deg): %f' % sd)
#                #plt.pause(0.001)
#                count += 1 
#        plt.suptitle('Pistoning waveguide %d' % scannedWG, verticalalignment='bottom')
#        plt.pause(.01)
#        
#        if saveFigs:
#            fname=savePref+"%d.png" % scannedWG       
#            plt.savefig(fname)
#        
#    plt.figure(3)
#    plt.clf()
#    plt.hist(np.reshape(allSDs, -1))
#    plt.xlabel('Closure phase SD (degrees)')
#    plt.ylabel('Frequency (across all triangles and piston choices')
#        
        
