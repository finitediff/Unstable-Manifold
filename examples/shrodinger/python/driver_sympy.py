"""This driver is designed to solve for the parameterization
functions for the beta model using Sympy.  It exports to 
a file the functions it produces"""

import sympy
from sympy import Function, Symbol
from sympy import diff
from sympy import Sum

"""___________________________________________________________
Define the class expSum."""

# Define the class expSum to be used for infinite summations.
class expSum(Function):

    # Multiplication with another expSum
    def __mul__(self, other):

        # If the other object is also an expSum
        if isinstance(other, expSum):
            i = Symbol('i')
            funct = self.args[1]
            iterVar = self.args[2]
            coeff1 = self.args[0].subs(iterVar, i)
            coeff2 = other.args[0].subs(iterVar, iterVar-i)

            mySum = Sum(coeff1 * coeff2, (i, 0, iterVar))
            return expSum(mySum, funct, iterVar)

        # Otherwise, follow the original function rules
        else:
            return super().__mul__(other)

    # expSum raised to a power.
    def __pow__(self, pow):
        out = 1
        for i in range(pow):
            out *= self
        return out

    # Differentiation, hard coded to work with Theta = sigma*e**(lambda*t)
    def diff(self, inVar):
        coeff1 = self.args[0]
        funct = self.args[1]
        iterVar = self.args[2]

        # If the derivative variable is t in P_n*theta(t)**n
        if (inVar in self.args[1].args):
            lamda = Symbol('lambda')
            return expSum(iterVar * coeff1 * lamda, funct, iterVar)

        # Otherwise, we're just priming the coefficients.
        else:
            return expSum(coeff1.diff(inVar), funct, iterVar)

    # Differentiation, hard coded to work with Theta = sigma*e**(lambda*t)
    def _eval_derivative(self, x):
        return self.diff(x)

"""___________________________________________________________
Set up for symbolic solving."""

# Define the derivative functions.
_x = lambda inVal: (inVal.diff(x))
_t = lambda inVal: (inVal.diff(t))
_xx = lambda inVal: _x(_x(inVal))

"""___________________________________________________________
Define the PDE."""

# Define spatial variables.
x = Symbol('x')
t = Symbol('t')

# Define system parameters
mu = Symbol('mu')
nu = Symbol('nu')

# Create system functions
u = Function('u')(x,t)
v = Function('v')(x,t)

# Create the PDE.
eq1 = _t(u) + _xx(v) - mu*v + (u**2 + v**2)*v
eq2 = -_t(v) + _xx(u) - u + (u**2 + v**2)*u - 2*nu*v
equations = (eq1, eq2)

"""___________________________________________________________
Substitute in parameterized functions."""

# Create the parameterized functions
theta = Function('theta')(t)
n = Symbol('n')
U_n = Function('U_n')(n, x)
V_n = Function('V_n')(n, x)

# Create the expSums with coefficient 
# "U_n", a function of "n" and "x", and "theta" a function of "t"
U_subs = expSum(U_n, theta, n)
V_subs = expSum(V_n, theta, n)

# Substitute in the expSums.
eq1 = eq1.subs(u, U_subs)
eq2 = eq2.subs(u, U_subs)
eq1 = eq1.subs(v, V_subs)
eq2 = eq2.subs(v, V_subs)

# The profile solution has been verified.
# U(x) = sqrt(2) sech(x)
# V(x) = 0

# Solve for P{n}

""" We can assume we're given the information for P_n-1.
We wish to derive the equation that lets us use this
to solve for the next values."""


# Given P{n} as a dictionary.
# Given 

