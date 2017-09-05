from pysb import *


def model_pysb():
    Monomer('A1',['a1','a2'])
    Monomer('A2',['a1','a3'])
    Monomer('A3',['a2'])
    
    Parameter('k_A1_A1', 2*2e6)
    Parameter('l_A1_A1', 1e-2)
    Parameter('k_A1_A2', 2.74e6)
    Parameter('l_A1_A2', 2e-4)
    Parameter('k_A2_A3', 2.1e7)
    
    
    Rule('A1_A1', A1(a1=None) + A1(a1=None) <> A1(a1=1)%A1(a1=1), k_A1_A1, l_A1_A1)
    Rule('A1_A2', A1(a2=None) + A2(a1=None) <> A1(a2=1)%A2(a1=1), k_A1_A2, l_A1_A2)
    Rule('A2_A3', A2(a3=None) + A3(a2=None) >> A2(a3=1)%A3(a2=1), k_A2_A3)
    
    Initial(A1(a1=None, a2=None), Parameter('A1_0', 1e-6))
    Initial(A2(a1=None, a3=None), Parameter('A2_0', 1.57e-7))
    Initial(A3(a2=None), Parameter('A3_0', 1e-6))
    
def model_old():
    Monomer('A1',['a1','a2'])
    Monomer('A2',['a1','a3'])
    Monomer('A3',['a2'])
    
    Parameter('k_A1_A1', 2e6)
    Parameter('l_A1_A1', 1e-2)
    Parameter('k_A1_A2', 2.74e6)
    Parameter('l_A1_A2', 2e-4)
    Parameter('k_A2_A3', 2.1e7)
    
    Rule('A1_A1', A1(a1=None) + A1(a1=None) <> A1(a1=1)%A1(a1=1), k_A1_A1, l_A1_A1)
    Rule('A1_A2', A1(a2=None) + A2(a1=None) <> A1(a2=1)%A2(a1=1), k_A1_A2, l_A1_A2)
    Rule('A2_A3', A2(a3=None) + A3(a2=None) >> A2(a3=1)%A3(a2=1), k_A2_A3)
    
    Initial(A1(a1=None, a2=None), Parameter('A1_0', 1e-6))
    Initial(A2(a1=None, a3=None), Parameter('A2_0', 1.57e-7))
    Initial(A3(a2=None), Parameter('A3_0', 1e-6))

def observe_abc_model():
    Observable('ta', A3(a2=None))
    Observable('tb', A2(a1=None,a3=None))
    Observable('tc', A1(a1=None,a2=None))
    Observable('tab', A3(a2=1) % A2(a3=1,a1=None))
    Observable('tbc', A2(a3=None,a1=1) % A1(a1=None,a2=1))
    Observable('tcc', A1(a1=1,a2=None) % A1(a1=1,a2=None))
    Observable('tabc', A3(a2=2) % A2(a1=1,a3=2) % A1(a1=None,a2=1))
    Observable('tbcc',A2(a1=1,a3=None) % A1(a1=3,a2=1) % A1(a1=3,a2=None))
    Observable('tabcc', A3(a2=1) % A2(a1=2,a3=1) % A1(a1=3,a2=2) % A1(a1=3, a2=None))
    Observable('tbccb', A2(a1=2,a3=None) % A1(a1=3,a2=2) % A1(a1=3, a2=1) % A2(a1=1,a3=None))
    Observable('tabccb', A3(a2=4) % A2(a1=2,a3=4) % A1(a1=3,a2=2) % A1(a1=3, a2=1) % A2(a1=1,a3=None))
    Observable('tabccba', A3(a2=4) % A2(a1=2,a3=4) % A1(a1=3,a2=2) % A1(a1=3, a2=1) % A2(a1=1,a3=5) % A3(a2=5))

def return_model(model_type):
    Model()
    if model_type=='pysb':
        model_pysb()
        observe_abc_model()
    if model_type=='old':
        model_old()
        observe_abc_model()
    return model