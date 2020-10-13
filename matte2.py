import lucidMatrixMain as lm

# Ta imot to matriser og gi produktet av dem som output. CHECK
# Ta imot to matriser og gi summen av dem. CHECK
# Ta imot to matriser og finne differansen mellom dem CHECK
# Finne determinanten til en matrise CHECK

# Finne den redusert trappeform til en matrise (Gauss-Jordan eliminasjon)




mx = lm.MatrixCalcs()

m0 = [[-2,-1,2,1],
      [5,5,-5,5],
      [6,8,2,4], 
      [3,6,9,12]]

test = [[-1,2,3],
        [4,5,6],
        [7,8,9]]

m1 = [[1,2],[3,4]]   
m2 = [[5,6],[7,8]]    
m3 = [1,2]
m4 = [3,4]  

ms = [[1, 1],
      [2, 3]]

ma = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]
mb = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]


prod = mx.matMult(m2, m1)
add = mx.matAdd(ma, mb)
sub = mx.matSub(m3, m4)
det = mx.matDet(m0)

print(det)
# print(prod)
# print("\n")
# print(add)
# print("\n")
# print(sub)


    
    
    
    
    
