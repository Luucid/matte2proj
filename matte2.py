import lucidMatrixMain as lm

# Ta imot to matriser og gi produktet av dem som output. CHECK
# Ta imot to matriser og gi summen av dem. CHECK
# Ta imot to matriser og finne differansen mellom dem CHECK
# Finne determinanten til en matrise CHECK
# Finne den redusert trappeform til en matrise (Gauss-Jordan eliminasjon) CHECK


m0 = [[-2,-1,2,1],
      [5,5,-5,5],
      [6,8,2,4], 
      [3,6,9,12]]

a = [[1,4],
     [4,5],
     [7,8]]

b = [[2,4],
     [4,5]]





mA = lm.MatrixCalcs(a, "3x2 matrix")
print(mA)
mAgaus = mA.gausJordan()
print(mAgaus)


mB = lm.MatrixCalcs(b, "2x2 matrix")
print(mB)
mBgaus = mB.gausJordan()
print(mBgaus)


# k = lm.MatrixCalcs(m0, "4x4 matrix")
# print(k)
# kGaus = k.gausJordan()
# print(kGaus)









    
    
    
    
    
