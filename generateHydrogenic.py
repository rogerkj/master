#!/usr/bin/env python

from sympy import *

from generateDerivatesHydrogenic import generateDerivatesHydrogenic

from orbitalGenerator import orbitalGenerator

def main():

    og = orbitalGenerator(5)

    gd = generateDerivatesHydrogenic()

    orbitals = og.getOrbitals()

    #  print "Orbitals generated"

    gd.generate(orbitals)

  #  print "Done!"

   
if __name__ == "__main__":
    main()
