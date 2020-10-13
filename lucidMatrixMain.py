import numpy as np

class MatrixCalcs():
    def __init__(self):
        
        #main matrixes#
        self.mat1 = 0
        self.mat2 = 0

        #dimentions# (n x m1) * (m2 x k) = (n x k)
        self.n = 0
        self.m1 = 0
        self.m2 = 0
        self.k  = 0
         
        #preps#
        self.prod = 0
        self.row = 0
        self.det = 0
        
        
    def __setMats(self, a, b = [0,0]):
        self.mat1 = a
        self.mat2 = b
        self.__setShape()
     
    def __genMat(self, d): 
        return np.array([[0]*d[1]]*d[0]) #return NxM array containing zeros.
    
    
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
    
    def __returnMatrix(self, m):
        if (self.__getShape(m)[0] == 1):
            return m[0]  
        return m
    
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
            
    
    
    
    def matMult(self, a, b):     
        self.__setMats(a, b)
        if(self.__errorCheck("mult")): 
            return 0
   
        for i in range(self.n):
            for j in range(self.k):
                if(self.n <= 1):
                    self.prod[i][j] = self.__rowXcol(self.mat1, self.__genRow(j))
                else:
                    self.prod[i][j] = self.__rowXcol(self.mat1[i], self.__genRow(j))

        return self.__returnMatrix(self.prod)
    
  
    def matAdd(self, a, b):
        self.__setMats(a, b)
        if(self.__errorCheck("add")):
            return 0
    
        matSum = self.__genMat(self.__getShape(self.mat1))
        
        for i in range(self.n):
            for j in range(self.m1):
                if(self.n <= 1):
                    matSum[i][j] = (self.mat1[j] + self.mat2[j])
                else:
                    matSum[i][j] = (self.mat1[i][j] + self.mat2[i][j])
                
        return self.__returnMatrix(matSum)
    
    
    
    def matSub(self, a, b):
        self.__setMats(a, b)
        if(self.__errorCheck("sub")): 
            return 0
        
        matSum = self.__genMat(self.__getShape(self.mat1))
        
        for i in range(self.n):
            for j in range(self.m1):
                if(self.n <= 1):
                    matSum[i][j] = (self.mat1[j] - self.mat2[j] )
                else:
                    matSum[i][j] = (self.mat1[i][j] - self.mat2[i][j])
                
        return self.__returnMatrix(matSum)
        
        

        
    def __detCalc(self, m): 
        a = (m[0][0]*m[1][1]) - (m[0][1]*m[1][0])
        # print("        (%d * %d) - (%d * %d) = %d |  \n\n" % (m[0][0],m[1][1], m[0][1],m[1][0], a))
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
        
               
    def matDet(self, a):
        self.__setMats(a)
        if(self.__errorCheck("det")): 
            return 0
        shape = self.__getShape(self.mat1)[0]
        det = 0
        if(shape > 2):
            for i in range(shape):
                x = self.mat1[0][i] * ((-1) ** i)
                det += x * (self.__shrinkDet(self.mat1, i, shape-1))
            return det
        return self.__detCalc(self.mat1)
            
        
    
    
    def redAlg(self, r, c):
        reduced = False
        tMat1 = np.transpose(self.mat1)
        
     
                    
    
        
    
    def gausJordan(self, a):
        # stage 1: find trappeform
        # stage 2: reduser
        self.__setMats(a)
        s = self.__getShape(self.mat1)
        # error check here.
        self.redAlg(s[0], s[1])
        
        
        
    


            
        
        
class lucidMath():
    def __init__(self):
        self.obj = MatrixCalcs()
        
    def product(self, a, b):
        return self.obj.matMult(a, b)
        
   
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
         # if(i%2 == 0):
         #            x = self.mat1[0][i] 
         #            print("x =", x)
         #            self.det += x * (self.__shrinkDet(self.mat1, i, shape-1))
         #        else:
         #            y = self.mat1[0][i] 
         #            print("y =", y)
         #            self.det -= y * (self.__shrinkDet(self.mat1, i, shape-1))
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
        
    







    
    



















