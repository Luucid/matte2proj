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
        
        
        
     
    def setMats(self, a, b):
        self.mat1 = a
        self.mat2 = b
        self.setShape()
        
        if(self.m1 != self.m2):                                               
            print("wrong matrix dimentions.")
            return 0

    def setShape(self):
        
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
              
    
    def getProdShape(self):
        a = np.shape(self.prod)
        if(len(a) <= 1):
            a = np.shape([self.prod])
        return a
  
    
    def __rowXcol(self, r, c):
        x = r*c
        return sum(x)
               
    
    def __genRow(self,c):                                                    
        for i in range(self.m2):
            self.row[i] = self.mat2[i][c]
        return self.row
    
    def __returnMatrix(self):
        if (self.getProdShape()[0] == 1):
            return self.prod[0]  
        return self.prod
    
    def matMult(self, a, b):     
        self.setMats(a, b)

        for i in range(self.n):
            for j in range(self.k):
                if(self.n <= 1):
                    self.prod[i][j] = self.__rowXcol(self.mat1, self.__genRow(j))
                else:
                    self.prod[i][j] = self.__rowXcol(self.mat1[i], self.__genRow(j))

        return self.__returnMatrix()
    
  
       
        
    
    
    def matPrint(self):   
        print(self.mat1, "*", self.mat2)
        print("(%dx%d)*(%dx%d)" % (self.n, self.m1, self.m2, self.k))
        



























