#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import cvxpy as cp
import numpy as np

x0 = np.array([1,1,0])
x = cp.Variable(3)

def norm_p(y):
    return np.sqrt((y[0] +y[1]+y[2])**2 + (y[0] -y[1]+y[2])**2  + (-2*y[0] +y[1]+y[2])**2 + (-y[0] -2*y[1]+y[2])**2 )

constraints = [
      x[0]   +x[1]+x[2]   >= 0, 
      x[0]   -x[1]+x[2]   >= 0,
   -2*x[0]   +x[1]+x[2]   >= 0,
     -x[0] -2*x[1]+x[2]   >= 0
]
prob = cp.Problem(cp.Minimize(cp.norm(x - x0,2)),constraints)
prob.solve()
print("f1 = ", prob.value)

Q = np.zeros([3,3])
for i in range(3):
    ei = np.array([0,0,0])
    ei[i] = 1
    for j in range(3):
        ej = np.array([0,0,0])
        ej[j] = 1
        Q[i,j] = (norm_p(ei+ej)**2)/4.0-(norm_p(ei-ej)**2)/4.0
    
prob = cp.Problem(cp.Minimize(cp.quad_form(x-x0, Q)),constraints)
prob.solve()
print("f2 = ", np.sqrt(prob.value))

print("f3 = ", np.sqrt(min(x0[0]+x0[1]+x0[2],0)**2 + min(x0[0]-x0[1]+x0[2],0)**2  + min(-2*x0[0]+x0[1]+x0[2],0)**2 + min(-x0[0]-2*x0[1]+x0[2],0)**2))
