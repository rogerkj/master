#!/usr/bin/env python

from sympy import *

from generateDerivatesGaussians import generateDerivatesGaussians

from parsebasisset import parsebasisset

def main():

    parser = parsebasisset()

    orbitals = parser.parse("ne.dat")

    gd = generateDerivatesGaussians(5)

    gd.generate(orbitals)

if __name__ == "__main__":
    main()
