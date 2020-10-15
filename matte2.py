import lucidMatrixMain as lm

# Ta imot to matriser og gi produktet av dem som output. CHECK
# Ta imot to matriser og gi summen av dem. CHECK
# Ta imot to matriser og finne differansen mellom dem CHECK
# Finne determinanten til en matrise CHECK
# Finne den redusert trappeform til en matrise (Gauss-Jordan eliminasjon) CHECK


# egenverdier CHECK
# egenvektorer


m0 = [[-2,-1,2,1],
      [5,5,-5,5],
      [6,8,2,4], 
      [3,6,9,12]]

# a = [[1,4,7,9],
#      [4,5,6,7],
#      [7,8,9,10]]









b = [[3,4],
     [2,1]]


mB = lm.MatrixCalcs(b, "2x2 matrix")
print(mB)
mBgaus = mB.gausJordan()
print(mBgaus)


x = lm.EgenVe(b)









    
    
    
    
    
