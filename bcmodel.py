#simplest model from MT1-MMP model (abc model)
from pysb import *
Model()

def monomer_bc_model():
    """Let b and c be TIMP2 and MT1-MMP. Monomer b has one binding site while
    monomer c has two binding sites"""    
    Monomer('b', ['b1'])
    Monomer('c', ['c1', 'c2'])

def rate_constant_bc_model():
    #default rate constants
    Parameter('kbc', 2.74e6)
    Parameter('lbc', 2e-4)
    Parameter('kcc', 2*2e6)
    Parameter('lcc', 1e-2)

def rule_original_bc_model():
    """monomer b has only 1 site that can bind monomer c.
    Monomer c has two sites where the one site can bind b and another one can 
    form dimer"""
    Rule('bc', b(b1=None) + c(c1=None) <> b(b1=1)%c(c1=1), kbc, lbc)
    Rule('cc', c(c2=None) + c(c2=None) <> c(c2=1)%c(c2=1), kcc, lcc)
    

def initial_condition_bc_model():
    """the real initial data for monomer be varies from 0.5e-6 to 1e-6.
    For monomer c we use the default initial data"""
    Initial(b(b1=None), Parameter('bo', 1.57e-7))
    Initial(c(c1=None, c2=None), Parameter('co', 1e-6))

def observe_bc_model():
    """By the rules, we get 6 species"""
    Observable('tb', b(b1=None))
    Observable('tc', c(c1=None, c2=None))
    Observable('tbc', b(b1=1)%c(c1=1,c2=None))
    Observable('tcc', c(c1=None,c2=1)%c(c1=None,c2=1))
    Observable('bcc', b(b1=1)%c(c1=1, c2=2)%c(c1=None,c2=2))
    Observable('bccb', b(b1=1)%c(c1=1, c2=2)%c(c2=2, c1=3)%b(b1=3))