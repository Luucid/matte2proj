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

test = [[2,4,6],
        [4,5,6],
        [7,8,9]]

m1 = [[1,2],[3,4]]   
m2 = [[5,6],[7,8]]    
m3 = [1,2]
m4 = [3,4]  

ms = [[4, 4],
      [2, 2]]

ma = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]
mb = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]



mA = lm.MatrixCalcs(test, "3x3 matrix")
mB = lm.MatrixCalcs(test, "3x3 matrix")
print(mA)
print(mB)

product = mA*mB
print(product)



add = mA+mB
print(add)

sub = mA-mB
print(sub)



print(mA)
print(mB)

gaus = mA.gausJordan()
print(gaus)


print(mA)
print(mB)



    
    
    
    
    
