#!/usr/bin/env python

import sys
import re
import UncertMath
import BinCluster

def probAdjustEquil( binProb, rates, uncert ):  # this function adjusts bin pops in binProb using rates and uncert matrices

    ### SWITCH BETWEEN SIMPLE CALCULATION AND WEIGHTED AVERAGE CALC
    fullCalcBins = 0  # 1 for weighted avg, 0 for simple calc
    fullCalcClust = 0 # 1 for weighted avg, 0 for simple calc
    threshhold = 0.0  # minimum weight (relative to max) for another value to be averaged
                      # only matters if fullCalcBins = 1 (or later perhaps if fullCalcClust = 1)
    
    # check rate matrix is square
    Ni = len(rates)
    Nj = len(rates[0])
    if Ni != Nj:
        print '\nWARNING: Not a square matrix!\n'
    
    # create new object: rate matrix of scalars with uncertainties
    ratesUncert = []
    for i in range(Ni):
        tmpRow = []
        for j in range(Nj):
            number = rates[i][j]
            min = rates[i][j]-uncert[i][j]
            max = rates[i][j]+uncert[i][j]
            tmpRow.append(   UncertMath.ScalarUncert( [ number, min, max ] )   )
        ratesUncert.append(tmpRow)
    
    # STEP 1a: Create matrix of ratios of probabilities based on DIRECT estimates
    # that is, ij element is p_i / p_j = k_ji / k_ij
    #print'#######  STEP 1a'
    ratiosDirect = [] 
    for i in range(Ni):
        tmpRow = []
        for j in range(Nj):
            ratio =  ratesUncert[j][i] / ratesUncert[i][j]  # ij element is p_i / p_j = k_ji / k_ij
            tmpRow.append( ratio )
        ratiosDirect.append(tmpRow)
    
    
    # STEP 1b: Create averaged matrix of ratios of probabilities based on both direct and indirect estimates
    # Indirect means '3rd bin' estimates: p_i / p_j = ( k_ki / k_ik ) ( k_jk / k_kj )
    # Turns out this is not helpful, so generally set fullCalcBins = 0 
    #print'#######  STEP 1b'
    ratiosAvg = [] 
    for i in range(Ni):
        tmpRow = []
        for j in range(Nj):
            #print 'i =', i, '  j =', j
            ijTmp = []
            maxWt = ratiosDirect[i][j].wt
            if maxWt=='UnDef':
                maxWt = 0.0
            if fullCalcBins==1:
                for k in range(Ni):
                    ratioIndirect = ( ratesUncert[k][i] / ratesUncert[i][k] )   *   ( ratesUncert[j][k] / ratesUncert[k][j] )
                    if not ratioIndirect.wt=='UnDef': 
                        #print '\nindirect:', ratioIndirect.wt, '  vs  max*thold:', maxWt*threshhold
                        if ratioIndirect.wt > (maxWt*threshhold): 
                            #print '... bins', i, j, 'indirect above threshhold'
                            #print '... compare values: indirect:', ratioIndirect.num, '  direct:', ratiosDirect[i][j].num
                            ijTmp.append(   ratioIndirect   ) # add indirect
            ijTmp.append(  ratiosDirect[i][j] ) # add direct estimate to the list
            tmpVec = UncertMath.VectorUncert(ijTmp)  # create vector object to use for weighted average
            ijAvg = UncertMath.ScalarUncert( tmpVec.WtAvg() )  # averaged scalar object
            #ijAvg = UncertMath.ScalarUncert( tmpVec.WtAvgCent() )  # averaged scalar object
            #print '\nFinal Average: i  j', i, j, ijAvg.num
            tmpRow.append( ijAvg )
        ratiosAvg.append( tmpRow )
    
    
    # STEP 2: Form clusters
    
    #print'#######  STEP 2'
    # STEP 2a: Sort probability ratios based on uncertainty
    # create list for sorting
    ijSort = []  # elements will be [uncertainty of ij or ji prob ratio, i or j, j or i]
                 # order ij chosen so that p_i < p_j
    for i in range(Ni):
        for j in range(i+1,Nj): # exclude ii pairs
            ijVal = ratiosAvg[i][j].num
            ijUncert = ratiosAvg[i][j].uncert
            jiVal = ratiosAvg[j][i].num
            jiUncert = ratiosAvg[j][i].uncert
            #print i, j, ijUncert, jiUncert
            #ratiosAvg[i][j].display()
            if ijVal < jiVal: # convention for sorting: only use ij pair for which ratio < 1
                if not ijUncert == 'UnDef':
                    ijSort.append( [ijUncert, i, j] )
            else:
                if not jiUncert == 'UnDef':
                    ijSort.append( [jiUncert, j, i] )
    
    ### sort ratios from lowest to highest uncertainty
    #print 'Before sort:'
    #print ijSort
    ijSort.sort()
    #print 'After sort:'
    #print ijSort
    
    # STEP 2b: Create initial ClusterList object for clustering
    tmpList = []  # object used for clustering
    ifClust = []  # track whether a bin has been clustered/group with other bins or not
    for i in range(Ni):
        tmpList.append( [ i, None, UncertMath.ScalarUncert( ['UnDef', 'UnDef', 'UnDef' ] ) ] )  # dummy value to start  
            # [ bin, cluster index, relative probability of bin in cluster]
            # last item is a place-holder ScalarUncert object which will get updated
        ifClust.append( 0 )  # zero indicates all bins initially unclustered
    #print '\nInitial Cluster List:'
    clusters = BinCluster.ClusterList( tmpList, ratiosAvg )
    #clusters.display()

    # note any bins with zero prob for omission from cluster process - 8/17/2011
    zeroes = []
    b = 0  # bin index
    for p in binProb:
        #print 'zero? bin',b,'  prob=',p
        if p==0:
            zeroes.append( b )
        b = b+1
    #print '... zeroes =', zeroes
    
    # join the pair with lowest uncertainty
    #print '\nJoining'
    flag = 'NotDone'
    while len(ijSort) > 0 and not flag=='Done':
        #print '\n', ijSort[0]
        ii = ijSort[0][1] # ii, jj are bins to be clustered
        jj = ijSort[0][2]
        # omit bin from clustering if zero prob -- 8/17/2011
        yesno = 'yes'
        for bb in zeroes:
            if ii==bb or jj==bb:
                yesno = 'no'
        #print 'yesno clustering? ii, jj, yesno', ii, jj, yesno
        if yesno=='yes':
            if fullCalcClust==1:
                flag = clusters.join( ii, jj )
            else:
                flag = clusters.joinSimple( ii, jj )
            #print 'flag=',flag
            #clusters.display()
            ifClust[ ii ] = 1  # note that bin has been clustered (8/17/2011: does not seem to be used)
            ifClust[ jj ] = 1  # note that bin has been clustered
        #print ifClust
        del ijSort[0]
    
    # Bin populations must be re-normalized -- each cluster separately to equal sum originally in cluster
    clustList = clusters.getList()
    #print clustList
    probSum = 0.0
    probInClust = 0.0
    for binList in clustList:  # loop over each cluster - most will be empty after full cluster combination process
        probSum = 0.0  # for summing prob in each cluster
        #print 'Next cluster:', binList
        if len(binList)>0:
            #print 'list of bins:', binList
            for bin in binList:  # calculate old sum of bin probabilities
                probSum = probSum + binProb[bin]
            probInClust = probInClust + probSum
            for bin in binList:  # re-normalize to maintain old sum of bin probabilities
                newBinPop = clusters[bin][2]  # new prob, but not normalized
                #print '... bin', bin, '  old prob =', binProb[bin], '  probSum = ', probSum
                binProb[bin] = probSum * newBinPop  # normalized using prob in cluster
                #print '... ... old new prob =', newBinPop, '   new new prob =', binProb[bin]
    print '... Prob in clusters =', probInClust
        
    # 7/18/2011 checking after Carsen's bug report - probabilities not normalized
    print '\nChecking bin probabilities'
    binCount = 0
    probSum = 0.0
    for p in binProb:
       print 'bin', binCount, '  binProb = ', p
       probSum = probSum + p
       binCount = binCount + 1
    print '..... Total Prob is', probSum
