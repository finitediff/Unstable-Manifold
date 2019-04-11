"""This driver is designed to solve for the parameterization
functions for the beta model using Sympy.  It exports to 
a file the functions it produces"""

import sympy
from sympy import Function, Symbol
from sympy import diff

# Create the symbols x,y,z and t, so we can take derivatives.
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
t = Symbol('t')

# Create the functions for the pde.
u = Function('u')(x,y,z,t)
rho = Function('rho')(x,y,z,t)
h = Function('h')(x,y,z,t)

# Create additional operation functions.
div = Function('div')
grad = Function('grad')
tensor = Function('tensor')

# Define the derivative functions.
_x = lambda inVal: (diff(inVal, x))
_y = lambda inVal: (diff(inVal, y))
_z = lambda inVal: (diff(inVal, z))
_t = lambda inVal: (diff(inVal, t))

# Define the Beta function.
eq1 = _t(rho) + div(rho*u)
print(eq1)





# Solve for P{n}

""" We can assume we're given the information for P_n-1.
We wish to derive the equation that lets us use this
to solve for the next values."""


# Given P{n} as a dictionary.
# Given 