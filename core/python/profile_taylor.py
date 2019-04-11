# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:22:26 2017

@author: blake
"""


# IMPORTS

import numpy as np
from mpmath import iv

# PRECISION CONTROL
iv.dps = 16

# PDE PARAMETERS
p = {'gamma':iv.mpf(1)/9}
p['alpha'] = 1/p['gamma']




# METHOD PARAMETERS

# truncation bound on parametrization method p_{MN}
M = 10
N = 10

# PARAMETRIZATION METHOD
Qu = np.empty((4,M,N), dtype = iv.mpf)

# fixed point
Qu[0][0][0] = iv.mpf(1)
Qu[1][0][0] = iv.mpf(0)
Qu[2][0][0] = iv.mpf(0)
Qu[3][0][0] = iv.mpf(0)

mu1_u = iv.sqrt(p['alpha'])
mu2_u = 1/iv.sqrt(p['gamma'])

# first eigenvector
Qu[0][1][0] = iv.mpf(1)
Qu[1][1][0] = mu1_u
Qu[2][1][0] = iv.mpf(0)
Qu[3][1][0] = iv.mpf(0)

# second eigenvector
Qu[0][0][1] = iv.mpf(0)
Qu[1][0][1] = iv.mpf(0)
Qu[2][0][1] = iv.mpf(1)
Qu[3][0][1] = mu2_u



Jmat = iv.zeros(4)
Jmat[0,1] = iv.mpf(1)
Jmat[1,0] = p['alpha']+Qu[2][0][0]**2
Jmat[1,2] = 2*Qu[0][0][0]*Qu[2][0][0] 
Jmat[2,3] = iv.mpf(1)
Jmat[3,0] = -Qu[2][0][0]**2/p['gamma']
Jmat[3,2] = (1-2*Qu[0][0][0]*Qu[2][0][0])/p['gamma']


I4 = iv.zeros(4)
for j in range(0,4):
    I4[j,j] = iv.mpf(1)







'''


a = iv.matrix([['0.1','0.3','1.0'],
                ['7.1','5.5','4.8'],
                ['3.2','4.4','5.6']])

b = iv.matrix([iv.mpc(4,1),iv.mpc(0.6,0),iv.mpc(0.5,0)])


c = iv.lu_solve(a, b)


print(a*c-b)
'''
