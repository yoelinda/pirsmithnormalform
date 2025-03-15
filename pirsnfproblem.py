from pirsmithnormalform import matrix
class pirSNFProblem:
    def __init__(self, A, debug=False):
        self.A = A.copy()
        self.elementT = type(A.get(0, 0))
        self.J = A.copy()
        self.S = matrix.pirMatrix.id(A.h,type(A.get(0,0)))
        self.T = matrix.pirMatrix.id(A.w,type(A.get(0,0)))
        self.Sinv=matrix.pirMatrix.id(A.h,type(A.get(0,0)))
        self.Tinv=matrix.pirMatrix.id(A.w,type(A.get(0,0)))
        self.debug = debug
    def isValid(self):
        if not self.S * self.A * self.T == self.J:
            return False
        zero = self.elementT.getZero()
        for i in range(self.J.h):
            for j in range(self.J.w):
                if i != j and self.J.get(i, j) != zero:
                    return False
        diagLength = min(self.J.h, self.J.w)
        if diagLength > 1:
            lastDiag = self.J.get(0, 0)
            seenZero = (lastDiag == zero)
            for i in range(1, diagLength):
                currentDiag = self.J.get(i, i)
                if currentDiag == zero:
                    seenZero = True
                if seenZero and currentDiag != zero: 
                    return False 
                if not seenZero and not *\ 
                  lastDiag.istotaldivisorof(currentDiag)[0]:
                    return False
        return True
    def cSwap(self, i, j): 
        if self.debug:
            print("cSwap call")
        if i == j:
            return
        for k in range(self.J.h):
            temp = self.J.get(k, i)
            self.J.set(k, i, self.J.get(k, j))
            self.J.set(k, j, temp)
        adjustment=matrix.pirMatrix.id(self.T.h,*\
          self.elementT)
        adjustment.set(i, i, self.elementT.getZero())
        adjustment.set(j, j, self.elementT.getZero())
        adjustment.set(i, j, self.elementT.getOne())
        adjustment.set(j, i, self.elementT.getOne())
        self.T = self.T * adjustment  
        self.Tinv =  adjustment * self.Tinv 
    def cLC(self,i,j,a,b,glcd=None,matrixM=None,Minv=None):
        if self.debug:
            print("cLC call")
        if glcd is None or a.isUnit():
            c = self.elementT.getZero()
            d = self.elementT.getOne()
            Minv = matrix.pirMatrix.id(2, type(a))
            Minv.set(1, 0, -b)
        else:
            c = matrixM.get(0,1)
            d = matrixM.get(1,1)
        for k in range(self.J.h):
            temp = self.J.get(k, i)
            self.J.set(k, i, self.J.get(k, i)*a *\
              + self.J.get(k, j)*b )
            self.J.set(k,j,temp*c+self.J.get(k,j)*d)
        adjustment=matrix.pirMatrix.id(self.T.h,*\
          self.elementT)
        adjustmentinv=matrix.pirMatrix.id(self.T.h,*\
          self.elementT)
        adjustment.set(i, i, a)
        adjustmentinv.set(i, i, Minv.get(0,0))
        if i != j:
            adjustment.set(j, i, b)
            adjustment.set(i, j, c)
            adjustment.set(j, j, d)
            adjustmentinv.set(j, i, Minv.get(1,0) )
            adjustmentinv.set(i, j, Minv.get(0,1) )
            adjustmentinv.set(j, j, Minv.get(1,1) )
        self.T = self.T * adjustment
        self.Tinv =  adjustmentinv * self.Tinv 
    def rSwap(self, i, j):
        if self.debug:
            print("rSwap call")
        if i == j:
            return
        for k in range(self.J.w):
            temp = self.J.get(i, k)
            self.J.set(i, k, self.J.get(j, k))
            self.J.set(j, k, temp)
        adjustment = matrix.pirMatrix.id(self.S.h, *\
          self.elementT)
        adjustment.set(i, j, self.elementT.getOne())
        adjustment.set(j, i, self.elementT.getOne())
        adjustment.set(i, i, self.elementT.getZero())
        adjustment.set(j, j, self.elementT.getZero())
        self.S = adjustment * self.S
        self.Sinv = self.Sinv * adjustment
    def rLC(self,i,j,a,b,grcd=None,matrixM=None,Minv=None):
        if self.debug:
            print("rLC call")
        if grcd is None or a.isUnit():
            c = self.elementT.getZero()
            d = self.elementT.getOne()
            Minv=matrix.pirMatrix.id(2,type(self.A.get(0,0)))
            Minv.set(0, 1, -b)
        else:
            c = matrixM.get(1,0)
            d = matrixM.get(1,1)
        for k in range(self.J.w):
            temp = self.J.get(i, k)
            self.J.set(i,k,a*self.J.get(i,k)+ *\
            b*self.J.get(j,k))
            self.J.set(j,k,c*temp+d*self.J.get(j, k))
        adjustment = matrix.pirMatrix.id(self.S.h, *\
          self.elementT)
        adjustmentinv = matrix.pirMatrix.id(self.S.h, *\
          self.elementT)
        adjustment.set(i, i, a)
        adjustmentinv.set(i, i, Minv.get(0,0))
        if i != j:
            adjustment.set(i, j, b)
            adjustment.set(j, i, c)
            adjustment.set(j, j, d)
            adjustmentinv.set(j, i, Minv.get(1,0) )
            adjustmentinv.set(i, j, Minv.get(0,1) )
            adjustmentinv.set(j, j, Minv.get(1,1) )
        self.S = adjustment * self.S
        self.Sinv = self.Sinv * adjustmentinv
    def makeitdiagonal(self):
        for i in range(min(self.J.h, self.J.w)):
            if self.J.get(i, i) == self.elementT.getZero():
                foundReplacement = False
                j = i
                k = i
                for j in range(i, self.J.h):
                    if foundReplacement:
                        j = j-1
                        break
                    for k in range(i, self.J.w):
                        if self.J.get(j,k)!= *\
                          self.elementT.getZero():
                            foundReplacement = True
                            break
                if not foundReplacement:
                    break
                else:
                    self.rSwap(i, j)
                    self.cSwap(i, k)
            doneIteration = False
            while not doneIteration:
                if self.J.get(i, i).isUnit():
                    break
                doneIteration = True
                for j in range(i + 1, self.J.h):
                    a= self.J.get(i, i)
                    b= self.J.get(j, i)
                    grcd, matrixM, Minv = self.J.get(i, *\
                      i).extended_grcd(self.J.get(j, i))
                    x = matrixM.get(0,0)
                    y = matrixM.get(0,1)
                    if self.J.get(i, *\
                      i).isrightUnitMultipleOf(grcd):
                        pass
                    elif self.J.get(j, *\
                      i).isrightUnitMultipleOf(grcd):
                        self.rSwap(i, j)
                        doneIteration = False
                    else:
                        self.rLC(i, j, x, y, grcd, *\
                          matrixM, Minv)
                        doneIteration = False
                for j in range(i + 1, self.J.w):
                    a= self.J.get(i, i)
                    b= self.J.get(i, j)
                    glcd, matrixM, Minv = self.J.get(i, *\
                      i).extended_glcd(self.J.get(i, j))
                    x = matrixM.get(0,0)
                    y = matrixM.get(1,0)
                    if self.J.get(i, *\
                      i).isleftUnitMultipleOf(glcd):
                        pass
                    elif self.J.get(i, \
                      j).isleftUnitMultipleOf(glcd):
                        self.cSwap(i, j)
                        doneIteration = False
                    else:
                        self.cLC(i, j, x, y, glcd, *\
                         matrixM, Minv)
                        doneIteration = False
            doneZeroing = False
            while not doneZeroing:
                doneZeroing = True
                for j in range(i + 1, self.J.h):
                    if self.J.get(j, i) != *\
                      self.elementT.getZero():
                        self.rLC(j,i,self.elementT.getOne(),
                                 -self.J.get(j, 
                                  i).__rightfloordiv__(
                                    self.J.get(i, i)))
                        if self.J.get(j, i) != 
                          self.elementT.getZero():
                            doneZeroing = False
                for j in range(i + 1, self.J.w):
                    if self.J.get(i, j) != 
                      self.elementT.getZero():
                        self.cLC(j,i,self.elementT.getOne(),
                                 -self.J.get(i, 
                                  j).__leftfloordiv__(
                                    self.J.get(i, i)))
                        if self.J.get(i, j) != 
                          self.elementT.getZero():
                            doneZeroing = False
        diagLength = min(self.J.h, self.J.w)
        for i in range(diagLength-1):
            ii = diagLength-i -1
            if self.J.get(ii,ii).a.degree<self.J.get(ii-1,
              ii-1).a.degree :
                self.cSwap(ii-1, ii)
                self.rSwap(ii-1, ii)
    def attempt_totaldivisors(self):    
        diagSize = min(self.J.w, self.J.h) - 1
        for ii in range(diagSize):
            i = diagSize - ii -1
            if self.J.get(i + 1, i + 1) == 
              self.elementT.getZero():
                return
            VEC0 = self.J.get(i,i).istotaldivisorof(
              self.J.get(i+1,i+1))
            if not VEC0[0]:
                if VEC0[2]=="left":
                    self.rLC(i,i+1,self.elementT.getOne(), 
                      self.elementT( VEC0[1] ) )
                    glcd, matrixM, Minv = self.J.get(i, 
                      i).extended_glcd(self.J.get(i, i+1))
                    x = matrixM.get(0,0)
                    y = matrixM.get(1,0)
                    self.cLC(i,i+1,x,y,glcd,matrixM,Minv)
                    grcd, matrixM, Minv = self.J.get(i, 
                      i).extended_grcd(self.J.get(i+1, i))
                    x = matrixM.get(0,0)
                    y = matrixM.get(0,1)
                    self.rLC(i,i+1,x,y,grcd,matrixM,Minv)
                    self.cLC(i+1,i,self.elementT.getOne(), 
                      -self.J.get(i,i+1).__leftfloordiv__(
                        self.J.get(i, i)))
                    self.rLC(i+1,i,self.elementT.getOne(), 
                      -self.J.get(i+1,i).__rightfloordiv__(
                        self.J.get(i, i)))
                if VEC0[2]=="right":
                    self.cLC(i,i+1,self.elementT.getOne(),
                      self.elementT(VEC0[1]))
                    grcd, matrixM, Minv = self.J.get(i, 
                      i).extended_grcd(self.J.get(i+1, i))
                    x = matrixM.get(0,0)
                    y = matrixM.get(0,1)
                    self.rLC(i,i+1,x,y,grcd,matrixM,Minv)
                    glcd, matrixM, Minv = self.J.get(i, 
                      i).extended_glcd(self.J.get(i, i+1))
                    x = matrixM.get(0,0)
                    y = matrixM.get(1,0)
                    self.cLC(i,i+1,x,y,glcd,matrixM,Minv)
                    self.rLC(i+1,i,self.elementT.getOne(), 
                      -self.J.get(i+1,i).__rightfloordiv__(
                        self.J.get(i, i)))
                    self.cLC(i+1,i,self.elementT.getOne(), 
                      -self.J.get(i,i+1).__leftfloordiv__(
                        self.J.get(i, i)))
    def computeSNF(self):    
        self.makeitdiagonal()
        while not self.isValid():
            self.attempt_totaldivisors()
        self.makemonic()