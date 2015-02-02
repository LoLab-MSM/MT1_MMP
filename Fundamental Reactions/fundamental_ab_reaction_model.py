"""This is reaction between molecule a and b to produce ab"""
#a + c <-> ac

from pysb import * #calling all pysb classes
Model() # calling Model class. We don't have to put arguments inside since it all has default in it (in core file)

Monomer('a',['b']) #calling monomer class in core file
Monomer('c',['b'])

##Diffusion
#default: diffusion constant is zero#
#Parameter('Da', 1)
#Parameter('Dc', 1)
#Diffusion('a(b=None), Da') #(species, value or parameter)
#Diffusion('a(c=None), Dc')

""""Solve RDE using fipy"""
#diffusion constants as a list
#D=[1,1,0]
Parameter('k',1)
Parameter('l',1e-10)

Rule('d', a(b=None) + c(b=None) <> a(b=1)%c(b=1), k, l)

Parameter('ao',100)
Parameter('co',200)
Initial(a(b=None),ao)
Initial(c(b=None),co)

Observable('dt', a(b=1)%c(b=1))
Observable('at', a(b=None))
Observable('ct', c(b=None))

