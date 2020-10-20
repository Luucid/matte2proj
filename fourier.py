import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

x = sy.Symbol('x')
π = sy.Symbol('π')
pi = np.pi



def prnt():
    

    plt.close()
    plt.pause(3)
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.axis([-8, 8, -3, 3])
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
    
    
    for i in range(1,n):
        j = i-1
     
        plt.plot([xv[j], xv[i]], [yv[j], yv[i]], "%s"%(curCol), linewidth=0.8)
        if(j % 5 == 0):
            plt.pause(0.05)
        
        if(i % lPops == 0):
            teiting = (teiting + 1) % 5
            if(teiting == 0):
                colCtr  += 1
                curCol = col[colCtr % cl] 
            
            
           
        
     
        
    # plt.figure(dpi=1200) 
    # plt.grid(True)
    # plt.plot(xv, yv, "r-", linewidth=0.5)
    plt.show()
    





def four(fx, pStart, pEnd, n):
    
    c0 = pStart
    c1 = pEnd
    T = abs(c0) + abs(c1)
    
    
    a0 = (1/T) * sy.integrate(fx, (x, c0, c1))     
    a0 = sy.sympify(a0)
    ffx = 0
    
    for i in range(1, n+1):
        an = (2/T) * sy.integrate(fx*sy.cos((2*pi*i*x)/T), (x, c0, c1)) 
        bn = (2/T) * sy.integrate(fx*sy.sin((2*pi*i*x)/T), (x, c0, c1))
        
        s = (an*sy.cos((2*pi*i*x)/T) + bn*sy.sin((2*pi*i*x)/T))

        ffx += s

    return a0 + ffx
        
    


def f(x):
    return x

fr = four(f(x), -1*(pi/2), pi/2, 4)


fr  = sy.sympify(fr)

print(fr)

n = 1000

xv = np.linspace(-8, 8, n)
yv = np.zeros(n)



for i in range(n):
    yv[i] = fr.subs(x, xv[i])


prnt()        








    
