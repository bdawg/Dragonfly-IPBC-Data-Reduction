# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 17:57:21 2015

@author: bnorris
"""

import numpy as np
import matplotlib.pyplot as plt

#allCPsMasked = np.ma.MaskedArray(allCPs, mask = allCPMask)


scannedWG = 5
cpInd = 11#27

count=0
plt.figure(2)
plt.clf()

dl = pistAmts/(2.*np.pi)*1.55

                
i0 = triangles[cpInd, 0]
i1 = triangles[cpInd, 1]
i2 = triangles[cpInd, 2]

curWgsIN = [bl2wg[i0,:], bl2wg[i1,:], bl2wg[i2,:]]
curWgsIN = np.asarray(curWgsIN)
#if np.sum(curWgsIN == scannedWG):
if True:#np.sum(curWgsIN == scannedWG):

    ### plt.subplot(7,8,count+1)
    plt.plot(dl, BLphaseABCD[i0,:, scannedWG]/np.pi*180., linestyle='--', label='Baseline 1')
    plt.plot(dl, BLphaseABCD[i1,:, scannedWG]/np.pi*180., linestyle='--', label='Baseline 2')
    plt.plot(dl, BLphaseABCD[i2,:, scannedWG]/np.pi*180., linestyle='--', label='Baseline 3')
    meanCP = np.mean(allCPsMasked[cpInd, :, scannedWG]) /np.pi*180.
    nmCP = allCPsMasked[cpInd, :, scannedWG]/np.pi*180. - np.mean(allCPsMasked[cpInd, :, scannedWG])
    if np.sum((nmCP > np.pi) + (nmCP < -np.pi)):
        nmCP[nmCP > np.pi] = nmCP[nmCP > np.pi]-2.*np.pi
        nmCP[nmCP < -np.pi] = nmCP[nmCP < -np.pi]+2.*np.pi
        nmCP = nmCP - np.mean(nmCP)
    plt.plot(dl, nmCP, linewidth=1, color='k', label='Closure Phase')
    
    #sd = np.std(allCPsMasked[cpInd, :, scannedWG])/np.pi*180.
    sd = np.std(nmCP)
    allSDs[count, scannedWG] = sd
    #print sd
    plt.title('SD (deg): %1.2f, mean (deg): %3.1f' % (sd, meanCP) )
    plt.ylabel('Phase (Degrees)')
    plt.xlabel(r'Phase delay applied $(\mu m)$')

    #plt.ylim(-10,10)
    plt.ylim(-200,200)
    plt.legend()    
    #plt.pause(0.001)
    count += 1 
    
#plt.suptitle('Pistoning waveguide %d' % scannedWG, verticalalignment='bottom')
    
  
plt.pause(.01)



plt.clf()
for i in range(8):
    plt.semilogy(singleWGscan_allChipOuts_avgd[:,i])
plt.ylim(1,1000)




