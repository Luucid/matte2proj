import lucidMatrixMain as lm
import taylor



a = [[0,2,2],
     [0,1,1],
     [-1,2,1]]


b = [[1,2],
     [-1,3]]

c = [[1,-2,3],
     [1,0,1],
     [1,3,-2]]


mB = lm.MatrixCalcs(a, "3x3 matrix")
print(mB)

mBgaus = mB.gausJordan()
print(mBgaus)

mBdet = mB.matDet()
print("determinant: ", mBdet)

print(mB * mB)
print(mB + mB)
print(mB - mB)

print("EIGEN VECTOR \n\n\n")
lm.EgenVe(c)

taylor.TaylorRekke("sin(x)", 9)











    
    
    
    
    
