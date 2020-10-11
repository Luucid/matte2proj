
import lucidMatrixMain as lm

# Ta imot to matriser og gi produktet av dem som output. CHECK
# Ta imot to matriser og gi summen av dem. CHECK
# Ta imot to matriser og finne differansen mellom dem

# Finne determinanten til en matrise
# Finne den redusert trappeform til en matrise (Gauss-Jordan eliminasjon)




mx = lm.MatrixCalcs()

m0 = [[1,2,3,4],[5,6,7, 8],[9, 10, 11, 12], [13,14,15,16]]
m1 = [[1,2],[3,4]]   
m2 = [[5,6],[7,8]]    
m3 = [1,2]
m4 = [3,4]  

ma = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]
mb = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]












prod = mx.matMult(m2, m1)
add = mx.matAdd(ma, mb)
sub = mx.matSub(m3, m4)
det = mx.matDet(m0)

# print(prod)
# print("\n")
# print(add)
# print("\n")
# print(sub)


    
    
    
    
    
