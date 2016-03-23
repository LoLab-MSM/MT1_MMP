from pysb import *
from pysb.bng import generate_equations
Model()

def monomer_abc_model():
    """Let a, b, and c be MMP2, TIMP2, and MT1-MMP, respectivelly.
    Monomer a has only one binding site. Each of monomer b and c has two sites"""
    Monomer('A',['b'])
    Monomer('B',['a','c'])
    Monomer('C',['b','c'])

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
    Rule('AB', A(b=None) + B(a=None) >> A(b=1)%B(a=1), kab)
    Rule('BC', B(c=None) + C(b=None) <> B(c=1)%C(b=1), kbc, lbc)
    Rule('CC', C(c=None) + C(c=None) <> C(c=1)%C(c=1), kcc, lcc)
    
def initial_condition_abc_model():
    #from the data
    Initial(A(b=None), Parameter('ao', 1e-6))
    Initial(B(a=None, c=None), Parameter('bo', 1.57e-7))
    Initial(C(b=None, c=None), Parameter('co', 1e-6))

monomer_abc_model()
rate_constant_abc_model()
rule_original_abc_model()
initial_condition_abc_model()
generate_equations(model)
print model.species

