from pysb import *


Model()

Monomer('C8', ['b'])
Monomer('Bid', ['b','s'],{'s':['u','p']})

Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)

Rule('C8Bid_bind', C8(b=None) + Bid(b=None, s='u') <> C8(b=1) % Bid(b=1, s='u'), kf, kr)
Rule('pBid', C8(b=1) % Bid(b=1, s='u') >> C8(b=None) + Bid(b=None, s='p'), kc)

Parameter('C8_O', 1000)
Parameter('Bid_O', 10000)
Initial(C8(b=None), C8_O)
Initial(Bid(b=None, s='u'), Bid_O)

Observable('uBid', Bid(b=None, s='u'))
Observable('obspBid', Bid(b=None, s='p'))
Observable('obsC8', C8(b=None))
Observable('C8Bid', C8(b=1) % Bid(b=1, s='u'))


