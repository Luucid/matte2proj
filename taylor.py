import numpy as np
import sympy as sy
import matplotlib.pyplot as plt


x = sy.Symbol('x')


            
class TaylorRekke():
    def __init__(self, f, n):
        self.f = str(f)
        self.n = n+1
        self.tr = np.zeros(self.n, dtype=sy.core.add.Add)
        self.__calcTayl()
        

    def plotTayl(self,a=-4, b=4, res=0.1):
        xVals = np.arange(a, b, res, dtype=sy.core.add.Add)
        yVals = np.zeros(len(xVals), dtype=sy.core.add.Add)
        for i in range(len(yVals)):
            yVals[i] = self.__funk(xVals[i])
            yVals[i] = sy.sympify(yVals[i])
            yVals[i] = yVals[i].subs('e', np.e)
        plt.grid(True)
        plt.plot(xVals, yVals, "r-")
        plt.show()
        plt.pause(10)
        plt.close()
       
    def __funk(self, v):
        return self.f.subs(x, v)     
            
    def __calcTayl(self):
        nfak = 1
        f = sy.sympify(self.f)
        g = f.subs(x, 0)
        self.tr[0] = (g/nfak)*(x**0)
        
        for n in range(1, self.n):
            nfak *= n
            f = sy.diff(f, x)
            g = f.subs(x, 0)
            self.tr[n] = (g/nfak)*(x**n)
            
        fs = str(self.tr[0])
        for n in range(1, len(self.tr)):
            if(self.tr[n] != 0):
                fs += " + " + str(self.tr[n])   
         
    
       
        print("\n\n")
        sy.pprint(sy.sympify(self.f))
        print("---------------\n\n")
        self.f = sy.sympify(fs) 
        sy.pprint(self.f)
        print("---------------\n\n")
        
        cstm = input("customize plot settings? (y/n): ")
        if(cstm == "y" or cstm == "Y"):
            a = float(input("enter start x: "))
            b = float(input("enter end x: "))
            c = float(input("enter resolution: "))
            self.plotTayl(a, b, c)
        else:
            self.plotTayl()
        


        







    
    