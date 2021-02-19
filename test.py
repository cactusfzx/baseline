import cvxpy as cp
import numpy as np

import mosek, sys

env = mosek.Env()
env.set_Stream(mosek.streamtype.log, lambda x: sys.stdout.write(x))
env.echointro(1)

x = cp.Variable ()

y = cp.Variable ()

t = cp.Variable(pos=True)

objective_fn = cp.ceil ( x + y + cp.abs(x-y))

objective = cp.Minimize ( objective_fn+1 )

#constraint = cp.ceil ( x + y) <= t
print((cp.ceil ( x + y + cp.abs(x-y))+1).curvature)

problem = cp.Problem ( objective, [x<=2,y<=3,x>=0.5,y>=0.5] )

problem.solve ( qcp=True )

print ( "Optimal value: ", problem.value )
print ( "x: ", x.value )
print ( "y: ", y.value )

# Generate a random problem
np.random.seed(0)
m, n= 40, 25

A = np.random.rand(m, n)
b = np.random.randn(m)
# Construct a CVXPY problem
x = cp.Variable(n, integer=True)
objective = cp.Minimize(cp.sum_squares(A @ x - b))
prob = cp.Problem(objective)
prob.solve(solver=cp.MOSEK)

print("Status: ", prob.status)
print("The optimal value is", prob.value)
print("A solution x is")
print(x.value)

a = 1
print("test max",max(a,0))
print("test min",min(a,0))

one_cell_AoI = []
single_user_AoI = []

one_cell_AoI.append( np.random.randint(5,10,[3,2]))
print(one_cell_AoI[0][2,1])



