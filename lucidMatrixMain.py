import numpy as np

class MatrixCalcs():
    def __init__(self, m1, md = "nan"):
        
        #main matrixes#
        self.mat1 = m1
        self.mat2 = m1
        self.mat1org = self.mat1.copy()
       

        

        #dimentions# (n x m1) * (m2 x k) = (n x k)
        self.n = 0
        self.m1 = 0
        self.m2 = 0
        self.k  = 0
         
        #preps#
        self.prod = 0
        self.sub = None
        self.row = 0
        self.det = 0
        self.gaus = None
        self.mode = md
        
    
     ############overlasting av operatorer###################   
    
    def __str__(self):   
        s = "mode: %s" % self.mode
        s +="\n--------------------\n"
        for i in range(len(self.mat1)):
            s += str(self.mat1[i])
            s += "\n"
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
    
        
    def __setMats(self):
        self.mat1 = self.mat1org.copy()
        self.mat2 = self.mat1.copy()
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

        self.mat1 = self.mat1org
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
        self.mat1 = self.mat1org      
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
        
         
        self.mat1 = self.mat1org
        return matSum
        
    
        
    def matDet(self, a=None):
        # if(a != None):
        #     self.__setMats(a)
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
    
    
    
    
    def gausJordan(self):
        self.__setMats()
        print(self.mat1org)
        # error check here.
        rSize = len(self.mat1)      
        cSize = len(self.mat1[0])  
        
        self.__redAlgOne(rSize, cSize)
        print(self.mat1org)
        self.__redAlgTwo(rSize, cSize)
        print(self.mat1org)
        self.gaus = self.mat1.copy() 
        self.mat1 = self.mat1org
     
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
            for c in range(cSize):
                self.mat1[r1][c] *= k
            return
        
        if(mode == "add"):
            for c in range(cSize):
                tmp = self.mat1[r1][c] * k
                self.mat1[r2][c] += tmp
            return
 

    def __checkIfAllZero(self, rSize, cSize):
        for i in range(1,rSize):
            for j in range(cSize):
                if(self.mat1[i][j] != 0):
                    return 0
        return 1
                    
            
        
                
    def __redAlgOne(self, rSize, cSize):
        reduced = False  
        r1 = 0
        r2 = 1
        c = 0
        k = 0
        
        while(not reduced): 
            if(self.mat1[r1][c] == 0):
                if(self.__checkIfAllZero(rSize, cSize)):
                    reduced = True
                
            
            elif(self.mat1[r1][c] != 1):
                k = (1/self.mat1[r1][c])
                self.__gausLoop(rSize, k, "mult", r1)
                
            elif(self.mat1[r1][c] == 1):    
                for n in range(r2, rSize):
                    k = self.mat1[n][c]*(-1)
                    self.__gausLoop(cSize, k, "add", r1, n)
                
                r1 += 1
                c += 1
                r2 += 1
                if(r2 == rSize):
                    reduced = True
        
        if(self.mat1[rSize-1][cSize-1] != 1):
            if(self.mat1[rSize-1][cSize-1] != 0):
                k =  (1/self.mat1[rSize-1][cSize-1])
                self.__gausLoop(cSize, k, "mult", rSize-1)
                
            
            
                    
     
    def __redAlgTwo(self, rSize, cSize):
        reduced = False
        r1 = 0 
        r2 = 1
        c = 1 
   
        while(not reduced):
            if(self.mat1[r1][c] != 0):
                for n in range(r2, rSize):
                    if(self.mat1[n][c] != 0):
                        k = self.mat1[r1][c]*(-1)
                        self.__gausLoop(cSize, k, "add", n, r1)
                         
            r1 += 1
            c += 1
            r2 += 1
            if(r1 == rSize-2):
                reduced = True
            elif(r2 == rSize):
                reduced = True
    
          
                    

 
        
        
        
    

        
        
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
        
    







    
    



















