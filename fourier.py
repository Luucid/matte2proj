import numpy as np
import sympy as sy
import matplotlib.pyplot as plt



x = sy.Symbol('x')
π = sy.Symbol('π')
N = sy.Symbol('N')
pi = np.pi



class Fourier():
    '''
    (f(x), period-start, period-end, n-terms, customize settings = False)
    
    '''
    def __init__(self, fx, pStart, pEnd, n, cust = False):
        ################################
        self.n = n
        self.fx = sy.sympify(fx)
        self.c0 = pStart
        self.c1 = pEnd
        self.T = abs(pStart) + abs(pEnd)
        self.s = 0
        self.fr = None
        self.cust = cust
        ################################
        
        
        self.__four()
        self.__prnt(fx)
        
        
    

    def __prnt(self, lbl):      
        sy.pprint(self.fr)
        n = 1000  
        a = -8
        b = 8
        animate = False
        if(self.cust):
            a = int(input("enter xMin: "))
            b = int(input("enter xMax: "))
            n = int(input("enter resolution: "))
            animate = bool(input("turn on animation?(1/0): "))
        
        
        xv = np.linspace(a, b, n)
        yv = np.zeros(n)    
        f = sy.lambdify(x, self.fr)     
        for i in range(n):
            yv[i] = f(xv[i])
                
        
        
        lbl = sy.trigsimp(lbl)
        plt.close()
        plt.pause(3)
        plt.xlabel("Fourier-rekke for f(x) = " + str(lbl))
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.subplot().spines["left"].set_visible(0)
        plt.subplot().spines["right"].set_visible(0)
        plt.subplot().spines["top"].set_visible(0)
        plt.yticks([])
    
        col = ["b-", "c-", "g-", "y-", "r-"] 
        cl = len(col)
        colCtr = 0
        curCol = col[0]
        lPops = n//(8*np.pi)
        teiting = 0
        
        for i in range(1, n):
            j = i-1
         
            plt.plot([xv[j], xv[i]], [yv[j], yv[i]], "%s"%(curCol), linewidth=0.8)
            if(animate):
                if(j % 10 == 0):
                    plt.pause(0.05)
            
            if(i % lPops == 0):
                teiting = (teiting + 1) % 5
                if(teiting == 0):
                    colCtr  += 1
                    curCol = col[colCtr % cl] 
                
        
        plt.show()
        
    
    def __four(self):
        n = self.n
        c0 = self.c0
        c1 = self.c1
        T = self.T
        fx = self.fx
        
        
        a0 = (1/T) * sy.integrate(fx, (x, c0, c1)) 
        an = (2/T) * sy.integrate(fx*sy.cos((2*pi*N*x)/T), (x, c0, c1)) 
        bn = (2/T) * sy.integrate(fx*sy.sin((2*pi*N*x)/T), (x, c0, c1))    
        
        for i in range(1, n+1):        
            self.s += (an.subs(N, i)*sy.cos((2*pi*i*x)/T) + bn.subs(N, i)*sy.sin((2*pi*i*x)/T))
    
        self.fr = a0 + self.s
            
        
    
    
    
    

   
   



      

    
    
    
    
    
    
    
        
