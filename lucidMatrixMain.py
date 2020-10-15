import numpy as np
import sympy as sy

λ = sy.symbols('λ')


class MatrixCalcs():
    def __init__(self, m1, md = "nan"):
        
        #main matrixes#
        self.mat1 = np.array(m1, dtype=sy.core.add.Add)
        self.mat2 = np.array(m1, dtype=sy.core.add.Add)

        #dimentions# (n x m1) * (m2 x k) = (n x k)
        self.n = 0
        self.m1 = 0
        self.m2 = 0
        self.k  = 0
         
        #preps#
        self.gaus = np.array(self.mat1.copy(), dtype=sy.core.add.Add)
        
        self.prod = self.mat1.copy()
        self.det = None
        
        
        self.sub = None
        self.row = 0
       
        self.mode = md
        self.__setShape()
        
    
     ############overlasting av operatorer###################   
    
    def __str__(self):    
        s = "%s" % self.mode
        s +="\n--------------------\n"
        for i in range(self.n):
            s += "( "
            for j in range(self.m1):     
                s+= "["+str(self.mat1[i][j]) +"] "               
            s += ") \n"   
        s +="--------------------\n"
        return s
            
 
    def __mul__(self, b):       
        self.mat2 = b.mat1.copy() 
        self.__setShape()
        return MatrixCalcs(self.matMult(), "A * B")
               
    def __sub__(self, b):
        self.mat2 = b.mat1.copy() 
        self.__setShape()
        return MatrixCalcs(self.matSub(), "A - B")
                
    def __add__(self, b):
        self.mat2 = b.mat1.copy()    
        self.__setShape()
        return MatrixCalcs(self.matAdd(), "A + B")
    
   ##############################################################
    
        

     
    def __genMat(self, d): 
        return np.array([[0]*d[1]]*d[0], dtype=sy.core.add.Add) #return NxM array containing zeros.
    
    
    def __setShape(self):
        
        a = np.shape(self.mat1)
        b = np.shape(self.mat2)
        
        if(len(a) <= 1):
            a = np.shape([self.mat1])
        if(len(b) <= 1):
            b = np.shape([self.mat2])
       
        self.n = a[0]
        self.m1 = a[1]
        self.m2 = b[0]
        self.k = b[1]
        self.prod = np.array([[0]*self.k]*self.n)
        self.row = np.array([0]*self.m2) 
              
    
    def __getShape(self, m):
        a = np.shape(m)
        if(len(a) <= 1):
            a = np.shape([m])
        return a
  
    
    def __rowXcol(self, r, c):
        x = r*c
        return sum(x)
               
    
    def __genRow(self,c):                                                    
        for i in range(self.m2):
            self.row[i] = self.mat2[i][c]
        return self.row
    
 
    
    def __errorCheck(self, mode):
        if(mode == "mult"):
            if(self.m1 != self.m2):                                               
                print("ERROR: Wrong matrix dimentions.")
                return 1
            
        if(mode == "add"):
            if(self.__getShape(self.mat1) != self.__getShape(self.mat2)):                                               
                print("ERROR: Both matrixes must have the same dimentions.")
                return 1
            
        if(mode == "sub"):
            if(self.__getShape(self.mat1) != self.__getShape(self.mat2)):                                               
                print("ERROR: Both matrixes must have the same dimentions.")
                return 1   
         
        if(mode == "det"):
            s = self.__getShape(self.mat1)
            if(s[0] != s[1]):                                               
                print("ERROR: Not a square matrix.")
                return 1
            
        return 0
            
    
    
    
    def matMult(self):  
        if(self.__errorCheck("mult")): 
            return 0

        for i in range(self.n):
            for j in range(self.k):
                if(self.n <= 1):
                    self.prod[i][j] = self.__rowXcol(self.mat1, self.__genRow(j))
                else:
                    self.prod[i][j] = self.__rowXcol(self.mat1[i], self.__genRow(j))

        return self.prod
   
    
  
    def matAdd(self):    
        if(self.__errorCheck("add")):
            return 0
    
        matSum = self.__genMat(self.__getShape(self.mat1))
  
        for i in range(self.n):
            for j in range(self.m1):
                if(self.n <= 1):
                    matSum[i][j] = (self.mat1[j] + self.mat2[j])
                else:
                    matSum[i][j] = (self.mat1[i][j] + self.mat2[i][j])
         
        return matSum
    
    
    
    def matSub(self):        
        if(self.__errorCheck("sub")): 
            return 0
        
        matSum = self.__genMat(self.__getShape(self.mat1))
        for i in range(self.n):
            for j in range(self.m1):
                if(self.n <= 1):
                    matSum[i][j] = (self.mat1[j] - self.mat2[j] )
                else:
                    matSum[i][j] = (self.mat1[i][j] - self.mat2[i][j])
        
        return matSum
        
    
        
    def matDet(self):
        if(self.__errorCheck("det")): 
            return 0
        shape = self.__getShape(self.mat1)[0]
        
        det = 0
        if(shape > 2):
            for i in range(shape):
                x = self.mat1[0][i] * ((-1) ** i)
                det += x * (self.__shrinkDet(self.mat1, i, shape-1))
            return det
        self.det = self.__detCalc(self.mat1)
        return self.det
    
    
    
    
    def gausJordan(self):
        
        # error check here.
        rSize = len(self.gaus)      
        cSize = len(self.gaus[0])      
        self.__gausAlg(rSize, cSize)
        self.__redAlg(rSize, cSize)
     
        return MatrixCalcs(self.gaus, "Gaus-Jordan")

     

    



    def __detCalc(self, m): 
        a = (m[0][0]*m[1][1]) - (m[0][1]*m[1][0])
        return a
        
    
    def __fetchNextMat(self, mx, n, shape):
        shrkMx = np.zeros((shape, shape))
        for j in range(shape):
            for k in range(shape):
                if(k < n):
                    shrkMx[j][k] = mx[j+1][k]
                else:
                    shrkMx[j][k] = mx[j+1][k+1]
                    
        return shrkMx
        
        
    def __shrinkDet(self, mx, n, shape):

        shrkMx = self.__fetchNextMat(mx, n, shape)
        det = 0
        if(shape > 2):
            for i in range(shape):
                if(mx == self.mat1):
                    x = shrkMx[0][i] * ((-1) ** i)
                else:
                    x = mx[0][i] * ((-1) ** i)
           
                det += x * (self.__shrinkDet(shrkMx, i, shape-1))
            
            return det        
        return self.__detCalc(shrkMx)
        
               
    
            
        
    
    def __gausLoop(self, cSize, k, mode, r1, r2 = 0):
        
        if(mode == "mult"):     
           # print(self.gaus, "r%d * %d" %(r1,k))
            for c in range(cSize):                
                self.gaus[r1][c] *= k
           # print(self.gaus, "\n")
            return
        
        if(mode == "add"):    
           # print(self.gaus, "r%d  + r%d" % (r2, r1))
            for c in range(cSize):
                tmp = self.gaus[r1][c] * k
                self.gaus[r2][c] += tmp
           # print(self.gaus, "\n")
            return
 

    def __checkIfAllZero(self, rSize, cSize):
        for i in range(1,rSize):
            for j in range(cSize):
                if(self.gaus[i][j] != 0):
                    return 0
        return 1
                    
            
    def __checkIfGaused(self):
       
        if (self.n < self.m1):
            n = self.n
        else:
            n = self.m1
            
        for i in range(n):
            if(self.gaus[i][i] != 1 and self.gaus[i][i] != 0):
                return 0
        if(self.gaus[i][i] == abs(0)):
            return 2
        return 1
            
            
        
                
    def __gausAlg(self, rSize, cSize):
        gaused = False  
        r1 = 0
        r2 = 1
        c = 0
        k = 0
        
        while(not gaused): 
            if(self.gaus[r1][c] == 0):
                # gaused = self.__checkIfGaused()
                if(self.__checkIfAllZero(rSize, cSize)):
                    gaused = True
                
            
            elif(self.gaus[r1][c] != 1):
                k = (1/self.gaus[r1][c])
                self.__gausLoop(cSize, k, "mult", r1)
                
                
            elif(self.gaus[r1][c] == 1):    
                for n in range(r2, rSize):
                    k = self.gaus[n][c]*(-1)
                    self.__gausLoop(cSize, k, "add", r1, n)
                
                r1 += 1
                if(c+1 < cSize):
                    c += 1
                r2 += 1
                
                if(self.__checkIfGaused() == 1):
                    gaused = True
                elif(self.__checkIfGaused() == 2):
                    if(c+1 < cSize):
                        c += 1
                    else:
                        gaused = True
                    
                    
                
        
                
        
     
    def __redAlg(self, rSize, cSize): #redusert
        if(self.n < self.m1):
            cIts = rSize
        else:
            cIts = cSize
    
        for r in range(1, rSize):
            r2 = r
            cs = r
            for c in range(cs, cIts):
                if(self.gaus[r-1][c] == 0):
                    c = cSize
                    
                elif(self.gaus[r-1][c] != 0):
                    k = self.gaus[r-1][c]*(-1)
                    self.__gausLoop(cSize, k, "add", r2, r-1)
                if(r2+1 < rSize):
                    r2 +=1
            
                    
            
            
      
            
      
        
      
        
      
        
      
        
class EgenVe():
    def __init__(self, m1):
        r, c = np.shape(m1)
        self.eigVals = None
          
        self.lamMat = self.__fillMat(m1, λ)  
        
          
        
        self.a = MatrixCalcs(m1, "%dx%d matrix" %(r,c))
        self.b = MatrixCalcs(self.lamMat, "λ*I")
        self.c = None
        
        self.mat1 = self.a - self.b
        self.eigVects = []
        self.__solveDet()
       
        
    
        
 

    # def __str__(self):
        
    def __findEigVects(self, x):
        tmp = self.__fillMat(self.lamMat, x)
        self.c = MatrixCalcs(tmp, "λx*I")       
        self.mat1 = self.a - self.c
        vkt = self.mat1.gausJordan()
        print(vkt)
        
       
        
        
        
            
     
    
    def __fillMat(self, m1, x):
        mat = np.zeros(np.shape(m1), dtype=sy.core.add.Add)
        for i in range(len(m1)):
            mat[i][i] = x
        return mat
       
    def __solveDet(self):
         expr = sy.sympify(str(self.mat1.matDet()))
         expr = sy.expand(expr)
         res = sy.solve(expr)
         self.eigVals = res
         for i in range(len(self.eigVals)):
             self.__findEigVects(self.eigVals[i])
       
         
        
        
                
                
        
        
        
        
        
        
      
        
        
            
                
        
            
           
    
          
                    

 
        
        
        
    

        
        
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
        
    







    
    



















