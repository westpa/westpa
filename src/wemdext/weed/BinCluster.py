import UncertMath

class ClusterList(): # list of lists, with each element = [bin, cluster index, relative prob in cluster]
    # must pass prob ratios as second argument
    # assumes that no clusters exist initially
    def __init__(self, value1, value2):
        self.data = value1
        self.ratios = value2
        self.Nbin = len( self.data )
        #print 'For clustering, Nbin=',self.Nbin
        self.clusterContents = []  # a list of the existing clusters
        self.count = 0
        self.flag = 'NotDone'
    def __getitem__(self,index):
        return [ self.data[index][0], self.data[index][1], self.data[index][2].num, self.data[index][2].min, self.data[index][2].max ]
    def display(self):
        print 'ClusterList:'
        for x in self.data:
            print 'bin', x[0], '  cluster', x[1]
            x[2].display()
        print 'Clusters so far:', self.clusterContents
    def getList(self):
        return self.clusterContents
    def join(self, i, j):
        print 'In join function: joining bins', i, j
        # TO-DO: need to stop it if everything joined - could happen before all ratios considered
        if self.data[i][1]==None and self.data[j][1]==None: # both bins previously unjoined
            print '\n... New cluster being made'
            self.data[i][1] = self.count  # create new cluster
            self.data[j][1] = self.count
            rij = self.ratios[i][j]  # ScalarUncert object
            denom = rij.plusOne()
            self.data[i][2] = rij / denom  # relative probability for bin i
            self.data[j][2] = denom.recip()  # relative probability for bin j
            self.count = self.count + 1
            tmpList = [i,j]  # for sorting
            tmpList.sort()
            self.clusterContents.append( tmpList )
        elif self.data[i][1]==None or self.data[j][1]==None: # only one bin previously unjoined 
            print '\n... Singlet bin being joined to existing cluster'
            if self.data[i][1]==None:  # determine which bin is unjoined
                idum = i  # join idum to cluster containing jdum, by convention
                jdum = j
            else:
                idum = j  # join idum to cluster containing jdum, by convention
                jdum = i
            print 'adding bin', idum
            jClust = self.data[jdum][1]  # cluster to which bin idum will be joined
            vecTmp = []
            print 'Here1'
            for k in self.clusterContents[jClust]:  # loop over bins already in cluster
                rik = self.ratios[idum][k]
                pk = self.data[k][2]
                piTmp = rik * pk  # estimate for p_idum / P_cluster based on 'path' through bin k
                                  # Note that here P_cluster is value before addition of bin idum
                vecTmp.append( piTmp )
            print 'Here2'
            vecTmpUnc = UncertMath.VectorUncert( vecTmp )
            print 'Here3'
            #piTmpAvg = UncertMath.ScalarUncert( vecTmpUnc.WtAvgCent() )
            piTmpAvg = UncertMath.ScalarUncert( vecTmpUnc.WtAvg() )
            # now, compute relative prob of each in bin in *new* cluster (including bin idum)
            denom = piTmpAvg.plusOne()
            self.data[idum][2] = piTmpAvg / denom 
            for k in self.clusterContents[jClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / denom
            self.data[idum][1] = jClust # bookkeeping: joins cluster that jdum is already in
            self.clusterContents[jClust].append( idum )  # continue bookkeeping for joining process
            tmpList = self.clusterContents[jClust]  # for sorting
            tmpList.sort()
            self.clusterContents[jClust] = tmpList
        elif not self.data[i][1]==self.data[j][1]:  # if both bins previously in DIFFERENT clusters
            print '\n... Joining cluster to cluster ..........'
            iClust = self.data[i][1]  # cluster identities
            jClust = self.data[j][1]
            # for now use a single value for ijClustRatio; NEED TO AVERAGE THIS
            #rij = self.ratios[i][j]   # bin prob ratio
            #pi =  self.data[i][2]     # current fractional pop of bin in cluster
            #pj =  self.data[j][2]
            #ijClustRatio = rij * pj / pi  # prob ratio of i to j cluster
            # estimate ratio of cluster populations using all available bin pairs
            ratioList = []
            for k in self.clusterContents[iClust]:  # bins in iClust
                pk = self.data[k][2]
                for m in self.clusterContents[jClust]:  # bins in jClust
                    pm = self.data[m][2]
                    rkm = self.ratios[k][m]   # bin prob ratio
                    ijClustRatioTmp = rkm * pm / pk
                    print 'Cluster ratio estimate for', iClust, jClust, '=', ijClustRatioTmp.num, '   wt=', ijClustRatioTmp.wt
                    ratioList.append( ijClustRatioTmp )
            tmpVec = UncertMath.VectorUncert( ratioList )  # create vector object to use for weighted average
            #ijClustRatio = UncertMath.ScalarUncert( tmpVec.WtAvgCent() )  # averaged scalar object
            ijClustRatio = UncertMath.ScalarUncert( tmpVec.WtAvg() )  # averaged scalar object
            print '.... Averaged estimate:', ijClustRatio.num, '   wt=', ijClustRatio.wt
            # reassign relative prob of bins in each cluster 
            iDenom = ( ijClustRatio.recip() ).plusOne()
            for k in self.clusterContents[iClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / iDenom
            jDenom = ijClustRatio.plusOne()
            for k in self.clusterContents[jClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / jDenom
            # bookkeeping: convention jClust bins join iClust cluster
            for k in self.clusterContents[jClust]:  # loop over bins in jClust cluster
                self.data[k][1] = iClust
                self.clusterContents[iClust].append( k )
            tmpList = self.clusterContents[iClust]  # for sorting
            tmpList.sort()
            self.clusterContents[iClust] = tmpList
            length = len( self.clusterContents[iClust] )
            print 'length=',length
            if length == self.Nbin:
                self.flag = 'Done'
            self.clusterContents[jClust] = []       # 'zero out' jClust; removing it would change other indices
        else:
            print '\n... No joining.  Bins already in same cluster'
            print 'i in', self.data[i][1], '   j in', self.data[j][1]
        return self.flag
    def joinSimple(self, i, j):
        # SIMPLE VERSION: ONLY USES IJ PAIR WHEN TWO CLUSTERS ARE JOINED
        #print 'In join function: joining bins', i, j
        # TO-DO: need to stop it if everything joined - could happen before all ratios considered
        if self.data[i][1]==None and self.data[j][1]==None: # both bins previously unjoined
            #print '\n... New cluster being made'
            self.data[i][1] = self.count  # create new cluster
            self.data[j][1] = self.count
            rij = self.ratios[i][j]  # ScalarUncert object
            denom = rij.plusOne()
            self.data[i][2] = rij / denom  # relative probability for bin i
            self.data[j][2] = denom.recip()  # relative probability for bin j
            self.count = self.count + 1
            tmpList = [i,j]  # for sorting
            tmpList.sort()
            self.clusterContents.append( tmpList )
        elif self.data[i][1]==None or self.data[j][1]==None: # only one bin previously unjoined 
            #print '\n... Singlet bin being joined to existing cluster'
            if self.data[i][1]==None:  # determine which bin is unjoined
                idum = i  # join idum to cluster containing jdum, by convention
                jdum = j
            else:
                idum = j  # join idum to cluster containing jdum, by convention
                jdum = i
            #print 'adding bin', idum
            jClust = self.data[jdum][1]  # cluster to which bin idum will be joined
            #print 'Here1'
            # use only one 'k' value: k=jdum
            rik = self.ratios[idum][jdum]
            pk = self.data[jdum][2]
            piTmp = rik * pk  # estimate for p_idum / P_cluster based on 'path' through bin k
                              # Note that here P_cluster is value before addition of bin idum
            #print 'Here3'
            # now, compute relative prob of each bin in *new* cluster (including bin idum)
            denom = piTmp.plusOne()
            self.data[idum][2] = piTmp / denom 
            for k in self.clusterContents[jClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / denom
            self.data[idum][1] = jClust # bookkeeping: joins cluster that jdum is already in
            self.clusterContents[jClust].append( idum )  # continue bookkeeping for joining process
            tmpList = self.clusterContents[jClust]  # for sorting
            tmpList.sort()
            self.clusterContents[jClust] = tmpList
        elif not self.data[i][1]==self.data[j][1]:  # if both bins previously in DIFFERENT clusters
            #print '\n... Joining cluster to cluster ..........'
            iClust = self.data[i][1]  # cluster identities
            jClust = self.data[j][1]
            # THE SIMPLE PART: just use a single value for ijClustRatio
            rij = self.ratios[i][j]   # bin prob ratio
            pi =  self.data[i][2]     # current fractional pop of bin in cluster
            pj =  self.data[j][2]
            ijClustRatio = rij * pj / pi  # prob ratio of i to j cluster
            # reassign relative prob of bins in each cluster 
            iDenom = ( ijClustRatio.recip() ).plusOne()
            for k in self.clusterContents[iClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / iDenom
            jDenom = ijClustRatio.plusOne()
            for k in self.clusterContents[jClust]:  # loop over bins already in cluster
                self.data[k][2] = self.data[k][2] / jDenom
            # bookkeeping: convention jClust bins join iClust cluster
            for k in self.clusterContents[jClust]:  # loop over bins in jClust cluster
                self.data[k][1] = iClust
                self.clusterContents[iClust].append( k )
            tmpList = self.clusterContents[iClust]  # for sorting
            tmpList.sort()
            self.clusterContents[iClust] = tmpList
            length = len( self.clusterContents[iClust] )
            #print 'length=',length
            if length == self.Nbin:
                self.flag = 'Done'
            self.clusterContents[jClust] = []       # 'zero out' jClust; removing it would change other indices
        #else:
            #print '\n... No joining.  Bins already in same cluster'
            #print 'i in', self.data[i][1], '   j in', self.data[j][1]
        return self.flag

