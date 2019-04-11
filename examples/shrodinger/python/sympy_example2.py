from sympy import Sum, Function, Symbol, oo

class Sum2(Sum):

    # Overriding the __mul__ method.
    def __mul__(self, other):

        if isinstance(other, Sum2):
            i = Symbol('i')
            n = Symbol('n')
            return Sum2(Sum(self.args[0].subs(n, i)*other.args[0].subs(n, n-i), (i,0,n)), (n,0,oo))

        else:
            super().__mul__(other)

