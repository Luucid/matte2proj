import lucidMatrixMain as lm
import taylor
import fourier
from fourier import x, pi


a = [[0,2,2],
     [0,1,1],
     [-1,2,1]]


b = [[1,2],
     [-1,3]]

c = [[1,-2,3],
     [1,0,1],
     [1,3,-2]]


# mB = lm.MatrixCalcs(a, "3x3 matrix")
# print(mB)

# mBgaus = mB.gausJordan()
# print(mBgaus)

# mBdet = mB.matDet()
# print("determinant: ", mBdet)

# print(mB * mB)
# print(mB + mB)
# print(mB - mB)

# print("EIGEN VECTOR \n\n\n")
# lm.EgenVe(a)

# taylor.TaylorRekke("sin(x)", 9)

fourier.Fourier(x**2, -1*(pi/2), pi/2, 1000, False)









    
    
    
    
    
