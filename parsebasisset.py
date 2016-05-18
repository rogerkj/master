#!/usr/bin/env python

import re
import sys

import math

from sympy import (symbols,Symbol,exp,sympify,collect,rcollect,simplify,factor,Rational)

x, y, z , k = symbols('x y z k', real=True)


r2 , one = symbols('r2 one', real = True)


r2_2d = x**2 + y**2
r2_3d = r2_2d + z**2


class parsebasisset(object):

    def __init__(self):
        self.alpha_s = {}
        self.coeff_s = {}
        
        self.alpha_p = {}
        self.coeff_p = {}
    
        self.alpha_d = {}
        self.coeff_d = {}



    def normalize(self,alpha,expo):

        i = expo[0]
        j = expo[1]
        k = expo[2]

        N = (2.0*alpha/math.pi)**(0.75)*math.sqrt((8.0*alpha)**(i+j+k))
                                                  
        N *= self.factorial(i)*self.factorial(j)*self.factorial(k)/(self.factorial(2*i)*self.factorial(2*j)*self.factorial(2*k)) 
        
        return N

        

    def parse(self,filename):


        f = open(filename,"r")

        line = f.readline()

        while (line != "$basis\n"):
            line = f.readline();
        
        for i in range(4):
            line = f.readline()

        s = 0
        p = 0
        d = 0

        while (line != "*\n") :

            matchObj = re.match("\\s*(\d)\\s*(\w)\\s*",line)

            if (not matchObj):
                print "Orbital not found"
                sys.exit()
            
            orbital = matchObj.group(2)
            nr_el   = int(matchObj.group(1))

            if (orbital == "s"):
                self.alpha_s[s] = {}
                self.coeff_s[s] = {}

            elif (orbital == "p"):
                self.alpha_p[p] = {}
                self.coeff_p[p] = {}

            elif (orbital == "d"):
                self.alpha_d[d] = {}
                self.coeff_d[d] = {}

            else:
                print "orbital: " + orbital + " not supported!"
                sys.exit()


            for i in range(nr_el):
    
                line = f.readline()

                values = re.match("\\s*(-?\d+\\.?\d+)\\s*(-?\d+\\.?\d+)\\s*",line)

                if values:
        
                    if (orbital == "s"):
                        self.alpha_s[s][i] = float(values.group(1))
                        self.coeff_s[s][i] = float(values.group(2))

                    if (orbital == "p"):
                        self.alpha_p[p][i] = float(values.group(1))
                        self.coeff_p[p][i] = float(values.group(2))

                    if (orbital == "d"):
                        self.alpha_d[d][i] = float(values.group(1))
                        self.coeff_d[d][i] = float(values.group(2))


            line = f.readline()

            if (orbital == "s"):
                s += 1

            if (orbital == "p"):
                p += 1

            if (orbital == "d"):
                d += 1

        f.close()


        sum = []

        index = 0

        for i in range(len(self.alpha_s)):
            
            sum.append(0)

            for j in range(len(self.alpha_s[i])):
                N = (2.0*self.alpha_s[i][j]/math.pi)**(3.0/4.0)
                sum[index] += simplify((N*self.coeff_s[i][j] * exp(-(one*self.alpha_s[i][j])*r2_3d)))

            index += 1
        
        for i in range(len(self.alpha_p)):

            for t in range(3):

                sum.append(0)

                for j in range(len(self.alpha_p[i])):
                    N = (2.0*self.alpha_p[i][j]/math.pi)**(3.0/4.0)*2.0*math.sqrt(self.alpha_p[i][j])
                
                    m = 0
                    n = 0
                    o = 0

                    if t == 0:
                        m = 1

                    if t == 1:
                        n = 1

                    if t == 2:
                        o = 1

                    sum[index] += simplify(x**m*y**n*z**o*(N*self.coeff_p[i][j] * exp(-(one*self.alpha_p[i][j])*r2_3d)))

                index += 1


        expo = [[2,0,0],[0,2,0],[0,0,2],[1,1,0],[1,0,1],[0,1,1]]

        print str(len(expo))
        

        for i in range(len(self.alpha_d)):
                          
            for t in range(len(expo)):

                sum.append(0)

                for j in range(len(self.alpha_d[i])):
                
                    N = self.normalize(self.alpha_d[i][j],expo[t])
                
                    exponent = x**expo[t][0]*y**expo[t][1]*z**expo[t][2]

                    sum[index] += simplify(exponent*(N*self.coeff_d[i][j] * exp(-(one*self.alpha_d[i][j])*r2_3d)))

                index += 1


        print sum

        return sum




    def factorial(self,n):

        value = 1.0
        i = 1.0

        while (i < n):
            i += 1
            value *= i

        return value
    
