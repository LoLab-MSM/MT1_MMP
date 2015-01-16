from pysb import *

Model()

def monomer_abc_model():
    """Let a, b, and c be MMP2, TIMP2, and MT1-MMP, respectivelly.
    Monomer a has only one binding site. Each of monomer b and c has two sites"""
    Monomer('a',['a1'])
    Monomer('b',['b1','b2'])
    Monomer('c',['c1','c2'])

def rate_constant_abc_model():
    #default rate constants
    Parameter('kab', 2.1e7)
    Parameter('kbc', 2.74e6)
    Parameter('lbc', 2e-4)
    Parameter('kcc', 2*2e6)
    Parameter('lcc', 1e-2)

def rule_original_abc_model():
    """Monomer a can bind b. Monomer b can bind to monomer a and c on each sites.
    Monomer c can form dimer and bind b."""
    #binding criteria : (ab) b1 with a1, (bc) b2 with c1,(cc) c2 with itself
    Rule('ab', a(a1=None) + b(b1=None) >> a(a1=1)%b(b1=1), kab)
    Rule('bc', b(b2=None) + c(c1=None) <> b(b2=1)%c(c1=1), kbc, lbc)
    Rule('cc', c(c2=None) + c(c2=None) <> c(c2=1)%c(c2=1), kcc, lcc)


##model knockout 1
##remove ab's supplies
##remove forward reaction of ab + {c, cc, bcc, abcc}
##keep rule bc and cc
#Rule('ab', a(a1=None) + b(b1=None, b2=None) >> a(a1=1)%b(b1=1), kab)
#Rule('abc', a(a1=1)%b(b1=1, b2=None) + c(c1=None) <> a(a1=1)%b(b1=1, b2=2)%c(c1=2), kbc, lbc)


##model knockout 2
##remove abc's supplies
##remove forward reaction of abc + {c, bc, abc} and 
##keep rule bc and cc


def initial_condition_abc_model():
    #from the data
    Initial(a(a1=None), Parameter('ao', 1e-6))
    Initial(b(b1=None, b2=None), Parameter('bo', 1.57e-7))
    Initial(c(c1=None, c2=None), Parameter('co', 1e-6))

def observe_abc_model():
    """From the rules we have, in total, 12 components"""
    Observable('ta', a(a1=None))
    Observable('tb', b(b1=None, b2=None))
    Observable('tc', c(c1=None, c2=None))
    Observable('tab', a(a1=1)%b(b1=1, b2=None))
    Observable('tbc', b(b1=None, b2=1)%c(c1=1, c2=None))
    Observable('tcc', c(c1=None, c2=1)%c(c1=None, c2=1))
    Observable('abc', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=None))
    Observable('bcc', b(b1=None, b2=1)%c(c1=1, c2=2)%c(c1=None, c2=2))
    Observable('abcc', a(a1=1)%b(b1=1,b2=2)%c(c1=2,c2=3)%c(c1=None,c2=3))
    Observable('bccb', b(b1=None, b2=1)%c(c1=1, c2=2)%c(c2=2, c1=3)%b(b2=3, b1=None))
    Observable('abccb', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=3)%c(c2=3, c1=4)%b(b2=4, b1=None))
    Observable('abccba', a(a1=1)%b(b1=1, b2=2)%c(c1=2, c2=3)%c(c2=3, c1=4)%b(b2=4, b1=5)%a(a1=5))
