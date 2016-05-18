
from sympy import (symbols,
                   sqrt,
                   diff,
                   trigsimp,
                   assoc_legendre,
                   sin,
                   cos,
                   exp,
                   sympify,
                   printing,
                   Symbol,
                   Wild,
                   srepr)

from math import ceil,pow

from sympy.physics.hydrogen import R_nl

x, y, z, theta, phi, one = symbols('x y z theta phi one', real=True)

r2_2d = x**2 + y**2
r2_3d = r2_2d + z**2

r2d = sqrt(r2_2d)
r3d = sqrt(r2_3d)

r, r2, r_2d, k, x2, y2, z2 = symbols('r r^2 r_2d k x2 y2 z2', real=True, positive = True)

class orbitalGenerator(object):
    
    def __init__(self,nOrbitals):

        self.nOrbitals = nOrbitals

        self.nlm = {}

   

    def getOrbitals(self):

        orbitals = []

        self.nShells = int(ceil(pow(self.nOrbitals,1.0/3.0)) + 1)

        o = 0
        for n in range(1,self.nShells+1):
            for l in range(n):
                for m in range(-l,l+1):
                    self.nlm[o] = [n,l,m]
                    o += 1
    

        o = 0

        for i in range(len(self.nlm)):
            n,l,m = self.nlm[i]
       
            func = self.getSphericalFunc(l, m)

            func = R_nl(n,l,r3d,Z=k) * func
            
            func = func.simplify().collect(x).collect(y).collect(z)

            orbitals.append(self.removeConstants(func).collect(k))

            if o >= self.nOrbitals:
                break

            o += 1

        return orbitals

    def removeConstants(self, func):

        #Stips off constants
        newfunc = 1

        p = Wild('p')

        for arg in func.args:
            
            if ((not arg.is_constant()) and (not arg.match(k**p))):

                newfunc *= arg

        return newfunc
          
    def sphere2Cart(self, func):
          
        func = func.subs(sin(2*phi), sympify(2)*z*r2d/r3d**2)
        func = func.subs(cos(2*phi), 1 - r2d**2/r3d**2)
        func = func.subs(cos(2*theta), (x**2 - y**2)/r2d**2)
        func = func.subs(sin(2*theta), sympify(2)*x*y/r2d**2)
          
        func = func.subs(cos(theta), x/r2d)
        func = func.subs(sin(theta), y/r2d)
        func = func.subs(cos(phi), z/r3d)
        func = func.subs(sin(phi), r2d/(r3d))
          
        return func  
        
        
        
    def getSphericalFunc(self, l, m):
        if m < 0:
            m = abs(m)
            fac = sin(m*theta)
        else:
            fac = cos(m*theta)

        P = assoc_legendre(l, m, cos(phi))
        res = fac*P

        res = self.sphere2Cart(trigsimp(res))

        #Takes care of Abs when the argument is real..
        res = res.subs(r2d, r_2d).subs(r3d, r).subs(r, r3d).subs(r_2d, r2d)
        
        return res; 
