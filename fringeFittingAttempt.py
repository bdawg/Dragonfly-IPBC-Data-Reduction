# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 10:55:05 2015

@author: bnorris


Try a new approach: using fringe-fitting rather than V2PM.
This might suck for V2s, but perhaps the cross-terms are less crucial for
phases. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import curve_fit
import lmfit

nWgsOP = 120
nWgsIN = 8
nBLs = 28
nPists = 61
showPlots = True

inFile = 'extractedFluxes_ipbc_20140930_3_xtr2-nonlin6'
phaseFile = 'phases_ipbc_20140930_3_xtr2-nonlin6'

fitAllOnPhases = True # If false, restore allOnphases from allOnPhaseFile. 
                       # If true, fit phases and save to allOnPhaseFile.
allOnPhaseFile = 'AllOnPhases_ipbc_20140930_3_xtr2-nonlin6'


outputFittedPhaseFile = 'extractedPhases_amplGT1e-9_ipbc_20140930_3_xtr2-nonlin6_2'

#==============================================================================
# Restore flux values from file
#==============================================================================
npzfile=np.load(inFile+'.npz')
singleWGscan_allChipOuts_avgd = npzfile['singleWGscan_allChipOuts_avgd']
singleWGscan_allChipOuts = npzfile['singleWGscan_allChipOuts']
allOnWGscan_allChipOuts = npzfile['allOnWGscan_allChipOuts']
BLscan_allChipOuts = npzfile['BLscan_allChipOuts']
BLscan_pistVals=npzfile['BLscan_pistVals']

npzPfile=np.load(phaseFile+'.npz')
BLphases = npzPfile['BLphases']
BLamps = npzPfile['BLamps']


#==============================================================================
# First, match baselines to outputs by looking at BLamps of waveguide scans
#==============================================================================
cutoff = 50
BL2wgsOPBool = np.zeros((nWgsOP,nBLs),dtype=bool)
BL2wgsOP=[]
for i in range(nBLs):
    BL2wgsOPBool[:,i] = BLamps[:,i] > cutoff
    BL2wgsOP.append(np.nonzero(BL2wgsOPBool[:,i]))
BL2wgsOP = np.squeeze(np.asarray(BL2wgsOP))



#==============================================================================
# Now match input waveguides to outputs using singleWGscan_allChipOuts_avgd
# This the brightest 29 OP wgs (28 plus photometric)
#==============================================================================
wgsIN2WgsOP = []
nHighest=29
photLimLow = 3
photLimHi = 116
opNums=np.array(range(nWgsOP))

for i in range(nWgsIN): 
    cur = singleWGscan_allChipOuts_avgd[:,i]
    hiInds = (-cur).argsort()[:nHighest]
    
    # Remove photometric chans:
    cut = np.nonzero(hiInds <= photLimLow)
    hiInds = np.delete(hiInds, cut)
    cut = np.nonzero(hiInds >= photLimHi)
    hiInds = np.delete(hiInds, cut)
    
    wgsIN2WgsOP.append(hiInds)

    if showPlots:
        plt.figure(1)
        plt.clf()
        plt.plot(opNums,cur,'o')
        plt.plot(opNums[hiInds],cur[hiInds],'or')
        plt.plot(opNums[hiInds],cur[hiInds],'x',markersize=10)
        plt.pause(0.0001)

wgsIN2WgsOP = np.squeeze(np.asarray(wgsIN2WgsOP))



#==============================================================================
# Get the phases from AllOnScans, and/or save/restore to a file.
# Note: 'pistons' are always measured in radians. They have been converted from
# microns in V2PM_measure_fluxes.py
# XXX Actually, is this necessary? XXX
#==============================================================================
opInds = range(nWgsOP)
ipInds = range(nWgsIN)
pistAmts = BLscan_pistVals[:,0] #This assumes they're all the same! Fix this.
allOnPhases = np.zeros((nWgsOP,nWgsIN))
allOnAmps = np.zeros((nWgsOP,nWgsIN))
            
            
def sinFunc(x, sinParam0, sinParam1 ,sinParam2 ,sinParam3):
    # sinParams = [ampl, freq, phase, offset]
    y = sinParam0 * np.sin(sinParam1*x + sinParam2) + sinParam3
    return y
    
if fitAllOnPhases:
    count=1
    
    for ip in ipInds:
        print "Fitting waveguide scan %d" % ip
        if showPlots:
            plt.figure(2)  
            plt.clf()
            plt.imshow(allOnWGscan_allChipOuts[:,:,ip],interpolation='nearest')
            plt.pause(0.001)
            plt.figure(3)  
            plt.clf()
     
        for op in opInds:
            current = allOnWGscan_allChipOuts[:,op,ip]
            # Estimate starting params:
            offsetGuess = np.mean(current)
            amplGuess = np.std(current)*2./np.sqrt(2) #Check this is right...
            p0 = [amplGuess, 0.9, np.pi, offsetGuess]
            try:
                fitp, fitpcov = curve_fit(sinFunc, pistAmts, current, p0)
            except RuntimeError:
                print "Fit failed on waveguide %d, output %d." % (ip, op)
            allOnPhases[op,ip] = fitp[2]
            allOnAmps[op,ip] = fitp[0]
            #print fitp[1]
    
            if showPlots:
                plt.subplot(11,11,op+1)
                plt.plot(pistAmts,current,'+')
                plt.plot(pistAmts, sinFunc(pistAmts,fitp[0],fitp[1],fitp[2],fitp[3]), linewidth=1)
                #plt.ylim(0,10)
                    
    # Make amplitudes positive
    allOnPhases[allOnAmps < 0] = allOnPhases[allOnAmps < 0] + np.pi
    allOnAmps[allOnAmps < 0] = -1.*allOnAmps[allOnAmps < 0]
    
    # A -pi/2 offset to everything makes the mean closure phase zero (i.e. it
    # removes the offset.) But its still quite the same as the artifical v2pm. 
    # But thus doesn't affect CP stability, just offset. 
    # BLphases = BLphases - np.pi/2

    np.savez(allOnPhaseFile, allOnPhases=allOnPhases, allOnAmps=allOnAmps)

else:
    npzAOPfile=np.load(allOnPhaseFile+'.npz')
    allOnPhases = npzAOPfile['allOnPhases']
    allOnAmps = npzAOPfile['allOnAmps']
 


#==============================================================================
# TODO!!
# Obtain relative phase offsets between OP waveguides for each baseline, using
# allOnPhases data.
# TODO: Comapre this to phases from singleBLscan...
#==============================================================================



#==============================================================================
# Look at relative phase offsets between OP waveguides for each baseline, using
# singleWgScan data, and put in new variable.
#==============================================================================
relBLph = np.zeros((4,nBLs))

if showPlots:
    plt.figure(4)
    plt.clf()
        
for bl in range(nBLs):
    curInds=BL2wgsOP[bl,:]
    curPh = BLphases[curInds,bl]
    #print curPh/(np.pi)
    newPh = curPh - curPh[0]
    relBLph[:,bl] = newPh
    
    if showPlots:
        plt.figure(4)  
        #plt.clf()
        plt.subplot(5,6,bl+1)
        curPh=np.sort(curPh)
        plt.plot((curPh))
        plt.plot((curPh),'o')
        plt.xlim(-0.5, 3.5)
        plt.xlabel("Waveguide no. (sorted)")
        plt.ylabel("Phase")
        #plt.pause(0.5)
    
#pdb.set_trace()

#==============================================================================
# Find baseline phase for each piston state using ABCD fit
#==============================================================================

#nMeasBLs = 7 #This is the number of baselines to measure for each set of data, 
#    # i.e. those that correpond to the pistoned wavdeguide.


# This is (phase of this baseline, for this piston position, for this wg scan)
BLphaseABCD = np.zeros((nBLs, nPists, nWgsIN))
BLphaseABCDchisqr = np.zeros((nBLs, nPists, nWgsIN))
BLintensABCD = np.zeros((4, nBLs, nPists, nWgsIN))
xinds = np.arange(-2*np.pi,2.*np.pi,0.1) # For plotting



#def sinFunc(x, sinParam0, sinParam1 ,sinParam2 ,sinParam3):
#    # sinParams = [ampl, freq, phase, offset]
#    y = sinParam0 * np.sin(sinParam1*x + sinParam2) + sinParam3
#    return y
#


def sinFuncLMfit(params, x, data):
    # sinParams = [ampl, freq, phase, offset]
    parvals = params.valuesdict()
    y = parvals['ampl'] * np.sin(parvals['freq']*x + parvals['phas']) + parvals['offs']
    return y-data


for scannedWG in range(nWgsIN):
    for bl in range(nBLs):
        print "WG scan %d, Baseline %d" % (scannedWG, bl)
        #print "\a"
        curOPs = BL2wgsOP[bl,:]
        
        for pistInd in range(nPists):
            
            curWGouts = [] #The four waveguides corresponding to this baseline                
            curIntens = allOnWGscan_allChipOuts[pistInd, curOPs, scannedWG]
            BLintensABCD[:, bl, pistInd, scannedWG] = curIntens
            
            normdIntens = curIntens/np.mean(
                allOnWGscan_allChipOuts[:, curOPs, scannedWG],axis=0)       
            curIntens = normdIntens
            
            curPhOf = relBLph[:,bl]#+2.*np.pi #Make positive
#            print curPhOf
            
            offsetGuess = np.mean(curIntens)
            amplGuess = np.std(curIntens)*2./np.sqrt(2) #Check this is right...
#            print "Offsetguess: %f, Amplguess: %f " % (offsetGuess, amplGuess)


             # Do with curve_fit (numpy)
#            p0 = [amplGuess, 1., np.pi, offsetGuess]
#            fitp, fitpcov = curve_fit(sinFunc, curPhOf, curIntens, p0)
#            print fitp


            # Do with lmfit:
            x = curPhOf
            data = curIntens
            params=lmfit.Parameters()
            #                 (Name,  Value,     Vary, Min,  Max,  Expr)
            params.add_many( ('ampl', amplGuess, True, 0.3, None, None),
                             ('freq', 1.,        True, 0.9, 1.1, None),
                             ('phas', 0.,        True, None, None, None),
                             ('offs', offsetGuess, True, None, None, None) )
            # XXX: Since ampl>0.3, unmodulated BLs will have wrong amplitude.
            # If you care about these, perhaps run fit again with amp > 1e-9 
                             
            result = lmfit.minimize(sinFuncLMfit, params, args=(x,data))
            bestFit = data + result.residual
            print params.valuesdict()
            print result.chisqr
            print ''
            
            BLphaseABCD[bl, pistInd, scannedWG] = params.valuesdict()['phas']
            BLphaseABCDchisqr[bl, pistInd, scannedWG] = result.chisqr

            if showPlots:
                plt.figure(5)
                plt.clf()            
                plt.plot(curPhOf, curIntens,'o',markersize=20)
                plt.ylim(0,2)
                
                # For curve_fit:                
                #plt.plot(xinds, sinFunc(xinds, fitp[0], fitp[1], fitp[2], fitp[3]))    
                
                # For lmfit:
                emptyData = np.repeat(0., xinds.size)
                plotRes = sinFuncLMfit(params, xinds, emptyData)
                plt.plot(xinds, plotRes)
                plt.pause(0.001)
                
    plt.pause(1)


np.savez(outputFittedPhaseFile, BLphaseABCD = BLphaseABCD, 
                                BLphaseABCDchisqr = BLphaseABCDchisqr, 
                                BLintensABCD = BLintensABCD, 
                                BLphases = BLphases, 
                                BLamps = BLamps, 
                                BL2wgsOP = BL2wgsOP, 
                                wgsIN2WgsOP = wgsIN2WgsOP, 
                                pistAmts = pistAmts, 
                                allOnPhases = allOnPhases, 
                                allOnAmps = allOnAmps, 
                                relBLph = relBLph )
                          











