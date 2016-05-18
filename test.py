#!/usr/bin/env python

from sympy import *

x, y, z ,alpha= symbols('x y z alpha')

rsqrt = sqrt(x*x + y*y + z*z)

f = exp(-alpha * rsqrt)

dx = diff(f,x)
ddx = diff(diff(f,x),x)


print dx
print ddx
