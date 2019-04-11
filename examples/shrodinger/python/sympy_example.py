def multiplyThatSum(inputVal):
    return inputVal

from sympy import Sum, Function, Symbol, oo

# Define variables
n = Symbol('n')
x = Symbol('x')
t = Symbol('t')

# Define functions
theta = Function('theta')(t)
p = Function('p')(n,x)
q = Function('q')(n,x)

# Create Summations
pSum = Sum(p*theta**n, (n,0,oo))
qSum = Sum(q*theta**n, (n,0,oo))

# Multiply
out = pSum * qSum
print(out)

out = out.multiplyThatSum()
print(out)


Sum(Sum((p(i,x)*q(n-i,x))*theta**n, (i,0,n)), (n,0,oo))
