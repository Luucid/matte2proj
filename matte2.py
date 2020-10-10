import numpy as np

# Ta imot to matriser og gi produktet av dem som output. CHECK
# Ta imot to matriser og gi summen av dem.
# Ta imot to matriser og finne differansen mellom dem
# Finne determinanten til en matrise
# Finne den redusert trappeform til en matrise (Gauss-Jordan eliminasjon)


class MatrixCalcs():
    def __init__(self, a, b):
        self.mat1 = np.array(a)
        self.mat2 = np.array(b)
        
        #dimentions# (n x m1) * (m2 x k) = (n x k)
        self.n = len(self.mat1)
        self.m1 = len(self.mat1[0])
        self.m2 = len(self.mat2)
        self.k = len(self.mat2[0])
        
        #preps#
        self.prod = np.array([[0]*self.k]*self.n)
        self.row = np.array([0]*self.m2) 
   



    def matPrint(self, a):   
        for i in range(len(a)):
            print(a[i])
        print("\n")
    
  
    
    def __rowXcol(self, r, c):
        x = r*c
        return sum(x)
               
    
    def __genRow(self,c):                                                    
        for i in range(self.m2):
            self.row[i] = self.mat2[i][c]
        return self.row
    
    
    def matMult(self):                                                    
        if(self.m1 != self.m2):                                               
            print("wrong matrix dimentions.")
            return 0

        for i in range(self.n):
            for j in range(self.k):
                self.prod[i][j] = self.__rowXcol(self.mat1[i], self.__genRow(j))
        
        
        return self.prod
       
        
    
    
   
        

        





m1 = np.array([[1,2], [4,5]])
m2 = np.array([[1,2,3],[4,5,6]])
    
mx = MatrixCalcs(m1, m2)

prod = mx.matMult()

print(prod)


    
    
    
    
    
