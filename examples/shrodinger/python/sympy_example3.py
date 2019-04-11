from sympy import Sum, Function, Symbol, oo
import sympy
from sympy import Function, Symbol
from sympy import diff
from sympy import Sum

"""___________________________________________________________
Define the class expSum."""
class Sum2(Sum):

    # Overriding the __mul__ method.
    def __mul__(self, other):

        if isinstance(other, Sum2):
            i = Symbol('i')
            n = Symbol('n')
            return Sum2(Sum(self.args[0].subs(n, i)*other.args[0].subs(n, n-i), (i,0,n)), (n,0,oo))

        else:
            super().__mul__(other)


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
u = Function('U')(x,t)
v = Function('V')(x,t)

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
U_subs = Sum2(U_n*theta**n, (n, 0, oo))
V_subs = Sum2(V_n*theta**n, (n, 0, oo))

# Substitute in the expSums.
eq1 = eq1.subs(u, U_subs)
eq2 = eq2.subs(u, U_subs)
eq1 = eq1.subs(v, V_subs)
eq2 = eq2.subs(v, V_subs)