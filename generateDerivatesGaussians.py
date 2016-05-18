#!/usr/bin/env python

from sympy import *
from sympy.printing import ccode


nDimentions = 3


x0, y0, z0 = symbols('x0 y0 z0')

from parsebasisset import (x,y,z,k,one)

argument = symbols('argument')
rSingleParticle = symbols('rSingleParticle')


rsqrt = sqrt(x0*x0 + y0*y0 + z0*z0)

r = x0*x0 + y0*y0 + z0*z0

alpha = symbols('alpha')

    
class generateDerivatesGaussians(object):
    
    def __init__(self,nF):

        self.nFunctions = nF


    def generate(self,fx):

        print "Length" + str(len(fx))


        for i in range(len(fx)):
            fx[i] = fx[i].subs(x,x0,simultaneous=True).subs(y,y0,simultaneous=True).subs(z,z0,simultaneous=True).subs(k,alpha,simultaneous=True)

        

        cords = ["x0","y0","z0"]

        #dfx = zeros(self.nFunctions,nDimentions)
        
        dfx = [[0 for i in range(nDimentions)] for j in range(self.nFunctions)]
        d2fx = [0 for j in range(self.nFunctions)]
        
        for i in range(self.nFunctions):
            for j in range(nDimentions):
                exec "dfx[i][j] = diff(fx[i],"+cords[j]+")"

        for i in range(self.nFunctions):
            for j in range(nDimentions):
                exec "d2fx[i] += diff(diff(fx[i],"+cords[j]+"),"+cords[j]+")"


        for i in range(self.nFunctions):
            fx[i] = fx[i].subs(rsqrt,argument,simultaneous=True)
            fx[i] = fx[i].subs(r,rSingleParticle,simultaneous=True)
            fx[i] = fx[i].subs(one,sympify(1))

            d2fx[i] = d2fx[i].subs(rsqrt,argument,simultaneous=True)
            d2fx[i] = d2fx[i].subs(r,rSingleParticle,simultaneous=True)
            d2fx[i] = d2fx[i].subs(one,sympify(1))


        for i in range(self.nFunctions):
            for j in range(nDimentions):
                dfx[i][j] = dfx[i][j].subs(rsqrt,argument,simultaneous=True)
                dfx[i][j] = dfx[i][j].subs(r,rSingleParticle,simultaneous=True)
                dfx[i][j] = dfx[i][j].subs(one,sympify(1))
              

        common_dfx = [1 for j in range(self.nFunctions)]
        common_d2fx = [1 for j in range(self.nFunctions)]

        #dfx,common_dfx = self.find_common(dfx)
        #d2fx,common_d2fx = self.find_common(d2fx)

       #print d2fx
        fxcode =  ['' for j in range(self.nFunctions)]
        dfxcode_common = ['' for i in range(self.nFunctions)]
        d2fxcode_common = ['' for i in range(self.nFunctions)]
        
        d2fxcode = ['' for j in range(self.nFunctions)]


        w = Wild('w')

        for i in range(self.nFunctions):
            fxcode[i] = ccode(simplify(fx[i]))
            dfxcode_common[i] = ccode(simplify(common_dfx[i]))
            d2fxcode_common[i] = ccode(simplify(common_d2fx[i]))
            
            d2fxcode[i] = ccode(collect(d2fx[i],exp(w)))


        dfxcode = [['' for i in range(nDimentions)] for j in range(self.nFunctions)]
    

        for i in range(self.nFunctions):
            for j in range(nDimentions):
                dfxcode[i][j] = ccode(simplify(dfx[i][j]))
           

        for i in range(self.nFunctions):
            for j in range(nDimentions):
                dfxcode[i][j] = self.remove_pows(dfxcode[i][j])
               
                
        for i in range(self.nFunctions):
            fxcode[i] = self.remove_pows(fxcode[i])
            dfxcode_common[i] = self.remove_pows(dfxcode_common[i])
            d2fxcode_common[i] = self.remove_pows(d2fxcode_common[i])
            d2fxcode[i] = self.remove_pows(d2fxcode[i])


        for i in range(self.nFunctions):
            for n in range(nDimentions):
                exec "fxcode[i] = fxcode[i].replace('"+cords[n]+"','r(i,"+str(n)+")')"
                exec "dfxcode_common[i] = dfxcode_common[i].replace('"+cords[n]+"','r(i,"+str(n)+")')"
                exec "d2fxcode_common[i] = d2fxcode_common[i].replace('"+cords[n]+"','r(i,"+str(n)+")')"

                exec "d2fxcode[i] = d2fxcode[i].replace('"+cords[n]+"','r(i,"+str(n)+")')"


        for i in range(self.nFunctions):
            for j in range(nDimentions):
                for n in range(nDimentions):
                    exec "dfxcode[i][j] = dfxcode[i][j].replace('"+cords[n]+"','r(i,"+str(n)+")')"
                   


        

        print "dfxcode"

        for i in range(self.nFunctions):
            for j in range(nDimentions):
                print dfxcode[i][j]

        print "d2fxcode"

        for i in range(self.nFunctions):
            print d2fxcode[i]

        

                
        print "Printing to file"

        f = open('Gaussians.cpp','w')

        f.write("#include <armadillo>\n")
        f.write("#include <math.h>\n")
        f.write('#include "Gaussians.h"\n')
        f.write("using namespace arma;\n")
        f.write("using namespace std;\n")
        
        f.write("\n\n")
        

        f.write("void Gaussians::setR(double R){\n")
        f.write("}\n\n")
        


        f.write("Gaussians::Gaussians(int nDim,int nPart,int ch,double a,double b) {\n")
    
        f.write("\tnDimensions = nDim;\n")
        f.write("\tnParticles = nPart;\n")

        f.write("\tcharge = ch;\n")

        f.write("\talpha = a;\n")
        f.write("\tbeta = b;\n\n")

        f.write("}\n")
        
        f.write("\n\n");

        f.write("double Gaussians::waveFunction(const mat &r,int nParticle,int orbital) {\n\n")

        f.write("\tint i = nParticle;\n")
        f.write("\tdouble rSingleParticle = 0;\n")
        f.write("\tdouble argument = 0.0;\n\n")
        
        f.write("\tfor(int j = 0; j < nDimensions; j++) {\n")
        f.write("\t\trSingleParticle += r(i,j) * r(i,j);\n")
        f.write("\t}\n\n")

 #       f.write("\targument += sqrt(rSingleParticle);\n")
        
        f.write("\tswitch (orbital) {\n")

        for i in range(self.nFunctions):
            f.write("\tcase "+str(i)+":\n")
            f.write("\t\treturn "+fxcode[i] + ";\n")

        f.write("\t}\n}\n\n")

        f.write("rowvec Gaussians::gradient(const mat &r,int nParticle,int orbital) {\n\n")

        f.write("\tint i = nParticle;\n")
        f.write("\tdouble rSingleParticle = 0;\n")
        f.write("\tdouble argument = 0.0;\n\n")
  
        f.write("\tfor(int j = 0; j < nDimensions; j++) {\n")
        f.write("\t\trSingleParticle += r(i,j) * r(i,j);\n")
        f.write("\t}\n\n")

 #       f.write("\targument += sqrt(rSingleParticle);\n")

        f.write("\trowvec retvec(3);\n\n")

        f.write("\tswitch (orbital) {\n")

        for i in range(self.nFunctions):
            f.write("\tcase "+str(i)+":\n")
            for j in range(nDimentions):
                f.write("\t\tretvec("+str(j)+") = " + dfxcode[i][j]+ ";\n")
            f.write("\t\treturn retvec;\n\n")

        f.write("\t}\n}\n")

        f.write("double Gaussians::laplacian(const mat &r,int nParticle,int orbital) {\n\n")

        f.write("\tint i = nParticle;\n")
        f.write("\tdouble rSingleParticle = 0;\n")
        f.write("\tdouble argument = 0.0;\n\n")
  
        f.write("\tfor(int j = 0; j < nDimensions; j++) {\n")
        f.write("\t\trSingleParticle += r(i,j) * r(i,j);\n")
        f.write("\t}\n\n")

#        f.write("\targument += sqrt(rSingleParticle);\n")

        f.write("\trowvec retvec(3);\n\n")

        f.write("\tswitch (orbital) {\n")

        for i in range(self.nFunctions):
            f.write("\tcase "+str(i)+":\n")
            f.write("\t\treturn " + d2fxcode[i] + ";\n")
       

        f.write("\t}\n\n}\n")


        f.close()



    def make_list(self,infunc):

        outarr = []

        for i in range(len(infunc.args)):
            argstmp = infunc.args[i]
            outarr.append(argstmp)

        return outarr
        

    def find_common(self,dfxin):
        
        common_dfx = zeros(self.nFunctions)

#        dfx = zeros(self.nFunctions,nDimentions)
        dfx = [[0 for i in range(nDimentions)] for j in range(self.nFunctions)]

        for f in range(self.nFunctions):

            common = 1
           
            dfxtmp0 = self.make_list(dfxin[f][0].simplify())
            dfxtmp1 = self.make_list(dfxin[f][1].simplify())
            dfxtmp2 = self.make_list(dfxin[f][2].simplify())


            #       print ""
            #       print str(len(dfxtmp0))
        
            #      for i in range(len(dfxtmp0)):
            #          print dfxtmp0[i]
    

            i = 0

            while (i < len(dfxtmp0)):
                
                found_common = False
                
                j = 0

                while (j < len(dfxtmp1)):
                    
                    l = 0

                    while (l < len(dfxtmp2)):
                        if (dfxtmp0[i].equals(dfxtmp1[j]) and dfxtmp0[i].equals(dfxtmp2[l])):
                            if (not found_common):
                                
                              
                                mul = dfxtmp0.pop(i)
                                common *= mul
                            
                                dfxtmp1.pop(j)
                                dfxtmp2.pop(l)
                        
                                i = -1
                                j = len(dfxtmp1)
                                l = len(dfxtmp2)
                            
                                
                            
                            found_common = True
                            
                    
                        l += 1
                    j += 1
                i += 1

            common_dfx[f] = common

            print common
        
            tmp = 1
            for i in range(len(dfxtmp0)):
                tmp *= dfxtmp0[i]
                
            dfx[f][0] = tmp

            tmp = 1
            for i in range(len(dfxtmp1)):
                tmp *= dfxtmp1[i]
                
            dfx[f][1] = tmp

            tmp = 1
            for i in range(len(dfxtmp2)):
                tmp *= dfxtmp2[i]
                
            dfx[f][2] = tmp

            

        return dfx,common_dfx


    def remove_pows(self,string):
        
        cords = ['x0','y0','z0','alpha','argument']

        for i in range(len(cords)):
            for n in range(1,5):
                
                powstr = "(" + cords[i]
                for m in range(n-1):
                    powstr += "*" + cords[i]

                powstr = powstr + ")"

                string = string.replace("pow(" + cords[i] + ", "+str(n) + ")",powstr) 


        return string


                                    
