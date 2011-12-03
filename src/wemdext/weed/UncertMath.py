from math import sqrt

class ScalarUncert():  # a 3-element array: ( x, x - dx(-), x + dx(+) )
    # NOTE: ASYMMETRIC UNCERTAINTY BUILT IN
    # ASSUMES ALL VALUES POSITIVE; WILL RE-SET MIN TO BE ABOVE ZERO
    def __init__(self,value):
        self.data = value
        self.num = self.data[0]
        self.min = self.data[1]
        self.max = self.data[2]
        self.neg = 0
        self.tiny = 1e-6  # used for setting near-zero minimum values
        if self.min <= 0.0:
            self.data[1] = self.tiny * self.num
            self.min = self.data[1]
        if self.num<0 or self.max<0:
            print 'Problem: Negative value!'
            print 'num =', self.num, '   max =', self.max
        if not self.min == self.max and not self.num==0:
            self.wt = (  self.num / ( self.max - self.min )  )**2 # weight based on fractional uncertainty
            # self.wt = 1.0 / ( self.max - self.min )**2
            self.uncert = ( self.max - self.min ) / self.num
        else: # if uncertainty not finite, assume data is worthless 
            #print 'Problem in ScalarUncert: num=',self.num,'  min=',self.min,'  max=',self.max
            self.wt = 'UnDef'
            #self.wt = 0.0
            self.uncert = 'UnDef'
            self.data[0] = 'UnDef'  # explicitly un-define to avoid future math problems -- e.g., division by zero
            self.data[1] = 'UnDef'
            self.data[2] = 'UnDef'
    # def __getitem__(self,i):
        # return self.data[i]
    def display(self):
        print 'the ScalarUncert value is ', self.data
    def __add__(self,other):
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef' or other.num=='UnDef' or other.min=='UnDef' or other.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = self.num + other.num  # simple addition
            lower = self.min + other.min   # simple minimum value
            upper = self.max + other.max   # simple maximum value
        if self.neg==0:
            if number<0 or lower<0 or upper<0:
                self.neg = 1
                print 'Add: Negative result found:', number, lower, upper
                print 'Self:', self.num, self.min, self.max
                print 'Other:', other.num, other.min, other.max
        return ScalarUncert([number,lower,upper])
    def __sub__(self,other):
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef' or other.num=='UnDef' or other.min=='UnDef' or other.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = self.num - other.num  # simple subtraction
            lower = self.min - other.max   # simple minimum value
            upper = self.max - other.min   # simple maximum value
        if self.neg==0:
            if number<0 or lower<0 or upper<0:
                self.neg = 1
                print 'Sub: Negative result found:', number, lower, upper
                print 'Self:', self.num, self.min, self.max
                print 'Other:', other.num, other.min, other.max
        return ScalarUncert([number,lower,upper])
    def __mul__(self,other):
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef' or other.num=='UnDef' or other.min=='UnDef' or other.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = self.num * other.num  # simple multiplication
            lower = self.min * other.min   # simple minimum value
            upper = self.max * other.max   # simple maximum value
        if self.neg==0:
            if number<0 or lower<0 or upper<0:
                self.neg = 1
                print 'Mul: Negative result found:', number, lower, upper
                print 'Self:', self.num, self.min, self.max
                print 'Other:', other.num, other.min, other.max
        return ScalarUncert([number,lower,upper])
    def __div__(self,other):
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef' or other.num=='UnDef' or other.min=='UnDef' or other.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        elif other.num==0 or other.min==0 or other.max==0:
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = self.num / other.num  # simple division
            lower = self.min / other.max   # simple minimum value
            upper = self.max / other.min   # simple maximum value
        if self.neg==0:
            if number<0 or lower<0 or upper<0:
                self.neg = 1
                print 'Div: Negative result found:', number, lower, upper
                print 'Self:', self.num, self.min, self.max
                print 'Other:', other.num, other.min, other.max
        return ScalarUncert([number,lower,upper])
    def plusOne(self):
        # adds one to number and uncertainties -- i.e., adds unity that has no uncertainty
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = self.num + 1.0  # simple addition
            lower = self.min + 1.0   # simple minimum value
            upper = self.max + 1.0   # simple maximum value
        return ScalarUncert([number,lower,upper])
    def recip(self):
        # takes reciprocal
        if self.num=='UnDef' or self.min=='UnDef' or self.max=='UnDef':
            number = 'UnDef'
            lower = 'UnDef'
            upper = 'UnDef'
        else:
            number = 1.0 / self.num  # simple addition
            lower = 1.0 / self.max   # simple minimum value
            upper = 1.0 / self.min   # simple maximum value
        return ScalarUncert([number,lower,upper])

class VectorUncert():  # an array of numbers and uncertainties: ( x1, x1 - dx1(-), x1 + dx1(+)  ,  x2, x2 - dx2(-), x2 + dx2(+)  , ...   )
    # NOTE: ASYMMETRIC UNCERTAINTY BUILT IN
    def __init__(self,value):
        self.data = value
    def WtAvg(self):
        #print 'in WtAvg function ....'
        norm = 0
        sum = 0
        wtmax = 0
        dsum = 0
        count = 0
        self.tiny = 1e-6

        # first obtain average and also check that at least one value was defined
        for x in self.data:               # simple weighted averages: number, lower, upper
            # MATT: I'm not sure if I used if syntax as intended below:
            if not x.wt=='UnDef' and not x.num=='UnDef' and not x.min=='UnDef' and not x.max=='UnDef':
                #x.display()
                #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                count = count + 1
                norm = norm + x.wt            # normalization for weighted average
                sum = sum + x.wt * x.num      # sum for weighted average
        #print 'count of nonUnDef =', count
        if count > 0: # be sure at least one value was defined
            #print 'sum=',sum,'  norm=',norm
            avg = sum / norm
            #print 'avg=', avg
            # average is used to calculate weighted uncertainty
            for x in self.data:
                #x.display()
                #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                if not x.wt=='UnDef' and not x.num=='UnDef' and not x.min=='UnDef' and not x.max=='UnDef':
                    #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                    if x.num > avg:                # seeking max deviation from average
                        term = ( x.max - avg )**2  # if above average, use max
                    else:
                        term = ( avg - x.min )**2  # if below average, use min
                    #print 'wt =', x.wt, '   term =', term
                    dsum = dsum + x.wt * term      # sum for weighted average of deviations
                    if x.wt > wtmax:                 # maximum weight for effective n
                        wtmax = x.wt
            #print 'wtmax=', wtmax
            neff = norm / wtmax  # effective number of samples based on uncertainties
            #print 'avg=', avg, '   norm=', norm, '   neff =', neff
            dev = 0.5 * sqrt( dsum / ( norm * neff ) )
            tmpMin = avg-dev
            if tmpMin < 0:
                tmpMin = self.tiny * avg  # a tiny fraction of average, so as to be > 0
            tmpMax = avg+dev
            return [avg, tmpMin, tmpMax ]  # just a list, not a ScalarUncert as written
        else:
            #print 'Problem in WtAvg due to count=0'
            return ['UnDef', 'UnDef', 'UnDef']
    def WtAvgCent(self):
        # modified weighted average which accounts not only for weights but range of possible values
        # Initially, a 'centralized' weighted average is calculated, which emphasizes high-weight values.
        # Then, modified values are re-averaged: values are moved closer to initial avg, within uncertainty range
        #print 'in WtAvg function ....'
        wtmax = 0
        dsum = 0
        count = 0
        self.tiny = 1e-6

        norm = 0
        sum = 0
        normCent = 0
        sumCent = 0
        # first obtain centralized average and also check that at least one value was defined
        wtList = []
        for x in self.data:               # simple weighted averages: number, lower, upper
            if (not x.wt=='UnDef') and (not x.num=='UnDef') and (not x.min=='UnDef') and (not x.max=='UnDef'):
                #x.display()
                #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                count = count + 1
                norm = norm + x.wt            # normalization for weighted average
                sum = sum + x.wt * x.num      # sum for weighted average
                wtsq = x.wt**8
                normCent = normCent + wtsq            # normalization for weighted average
                sumCent = sumCent + wtsq * x.num      # sum for weighted average
                wtList.append( [ x.wt, wtsq, x.wt, x.num, x.min, x.max, 0.0 ] )
        #print 'count of nonUnDef =', count
        if count > 0: # be sure at least one value was defined
            #print '\nNormalized list: [ wtNorm, wtsqNorm, wt, num, min, max ]:'
            wtList.sort()
            for z in wtList:
                z[0] = z[0] / norm
                z[1] = z[1] / normCent
                #print 'after:', z 
            tmpList = []
            #print 'sum=',sum,'  norm=',norm
            avg = sum / norm
            avgCent = sumCent / normCent
            #print 'In weighted avg function: avg=', avg, '   avgCent = ', avgCent
            # values are now modified and average is re-calculated
            norm = 0
            sum = 0
            #print '\n... Now re-averaging'
            for x in self.data:               # go through values again
                if (not x.wt=='UnDef') and (not x.num=='UnDef') and (not x.min=='UnDef') and (not x.max=='UnDef'):
                    #x.display()
                    #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                    norm = norm + x.wt            # normalization for weighted average
                    # change values based on avgInit
                    if x.max > avgCent and x.min < avgCent:  # centralized avg in range
                        numCent = avgCent           # use centralized value as representative value
                    elif x.max <= avgCent:          # centralized avg bigger than max
                        numCent =  x.max            # use max
                    elif x.min >= avgCent:          # centralized avg smaller than min
                        numCent = x.min             # use min
                    else:
                        print 'DANGER: unknown case'
                    tmpList.append( [ x.wt, x.num, x.min, x.max, numCent ] )
                    sum = sum + x.wt * numCent
            tmpList.sort()
            #print '... sorted items [ wt, num, min, max, numCent]:'
            #for y in tmpList:
                #print '*', y
            avg = sum / norm  # the new average based on values shifted centrally, within their ranges
            #print '... Final centralized average =', avg, '\n'
            # average is used to calculate weighted uncertainty
            for x in self.data:
                #x.display()
                #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                if not x.wt=='UnDef' and not x.num=='UnDef' and not x.min=='UnDef' and not x.max=='UnDef':
                    #print 'count=',count,'  wt=',x.wt,'  num=',x.num,'  min=',x.min,'  max=',x.max
                    if x.num > avg:                # seeking max deviation from average
                        term = ( x.max - avg )**2  # if above average, use max
                    else:
                        term = ( avg - x.min )**2  # if below average, use min
                    #print 'wt =', x.wt, '   term =', term
                    dsum = dsum + x.wt * term      # sum for weighted average of deviations
                    if x.wt > wtmax:                 # maximum weight for effective n
                        wtmax = x.wt
            #print 'wtmax=', wtmax
            neff = norm / wtmax  # effective number of samples based on uncertainties
            #print 'avg=', avg, '   norm=', norm, '   neff =', neff
            dev = 0.5 * sqrt( dsum / ( norm * neff ) )
            tmpMin = avg-dev
            if tmpMin < 0:
                tmpMin = self.tiny * avg  # a tiny fraction of average, so as to be > 0
            tmpMax = avg+dev
            return [avg, tmpMin, tmpMax ]  # just a list, not a ScalarUncert as written
        else:
            #print 'Problem in WtAvg due to count=0'
            return ['UnDef', 'UnDef', 'UnDef']

