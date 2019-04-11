import sympy
from sympy import expand, Function, Symbol, oo, diff, Sum, Derivative

"""A sum class overriding the multiplication operator to allow for double sums to be formed."""
class Sum2(Sum):
    
    #Class variables
    curr_i = 0

    # Overriding the __mul__ method.
    def __mul__(self, other):

        if isinstance(other, Sum2):
            i = Symbol('i_' + str(Sum2.curr_i))
            Sum2.curr_i += 1
            n = Symbol('n')
            return Sum2(Sum2(self.args[0].subs(n, i)*other.args[0].subs(n, n-i), (i,0,n)), (n,0,oo))

        else:
            super().__mul__(other)

"""A function that when called with a sympy expression will recursively look at each item
and multiply out all multpilications that can be done."""

def multiplyItOut(inEq):

    # Load the multiplication class.
    from sympy.core.mul import Mul as Mul

    # Recursively look at each subcomponent.
    items = [multiplyItOut(item) for item in inEq.args]

    # If it's a multiplyer object.
    if type(inEq) == Mul:
        out = 1
        for item in items:
            out *= item
        return out
    
    # If no arguments
    elif len(list(items)) == 0:
        return inEq

    # If it's not,
    else:
        return inEq.func(*list(items))


    

"""___________________________________________________
The Code ______________________________________"""
x = Symbol('x')
t = Symbol('t')
n = Symbol('n')

f = Function('f')(n, x)
a = Sum2(f*t**n, (n,0,oo))

# Works
print(a*a)

# Doesn't work.
c = (Derivative(a,x)*a).doit()
print("c =",c)
#print(c.doit())
#print(expand(c))

# Print the arguments.
print("Multiplied out", multiplyItOut(c))


