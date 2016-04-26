from pysb import *

def test():
    Monomer('A1',['a2'])
    Monomer('A2',['a1'])
    
    Parameter('k_A1_A2', 2.74e6)
    Parameter('l_A1_A2', 2e-4)
    
    Rule('A1_A2', A1(a2=None) + A2(a1=None) <> A1(a2=1)%A2(a1=1), k_A1_A2, l_A1_A2)
    
    Initial(A1(a2=None), Parameter('A1_0', 1e-6))
    Initial(A2(a1=None), Parameter('A2_0', 1.57e-7))

'''A_(n)-A_(n-1)-...-A_2-A_1>'''

def two_monomers():
    Monomer('A1',['a1','a2'])
    Monomer('A2',['a1'])
    
    Parameter('k_A1_A1', 2*2e6)
    Parameter('l_A1_A1', 1e-2)
    Parameter('k_A1_A2', 2.74e6)
    Parameter('l_A1_A2', 2e-4)
    
    Rule('A1_A1', A1(a1=None) + A1(a1=None) <> A1(a1=1)%A1(a1=1), k_A1_A1, l_A1_A1)
    Rule('A1_A2', A1(a2=None) + A2(a1=None) <> A1(a2=1)%A2(a1=1), k_A1_A2, l_A1_A2)
    
    Initial(A1(a1=None, a2=None), Parameter('A1_0', 1e-6))
    Initial(A2(a1=None), Parameter('A2_0', 1.57e-7))

def three_monomers():
    Monomer('A1',['a1','a2'])
    Monomer('A2',['a1','a3'])
    Monomer('A3',['a2'])
    
    Parameter('k_A1_A1', 2*2e6)
    Parameter('l_A1_A1', 1e-2)
    Parameter('k_A1_A2', 2.74e6)
    Parameter('l_A1_A2', 2e-4)
    Parameter('k_A2_A3', 2.1e7)
    Parameter('l_A2_A3', 1e-2)
    
    Rule('A1_A1', A1(a1=None) + A1(a1=None) <> A1(a1=1)%A1(a1=1), k_A1_A1, l_A1_A1)
    Rule('A1_A2', A1(a2=None) + A2(a1=None) <> A1(a2=1)%A2(a1=1), k_A1_A2, l_A1_A2)
    Rule('A2_A3', A2(a3=None) + A3(a2=None) <> A2(a3=1)%A3(a2=1), k_A2_A3, l_A2_A3)
    
    Initial(A1(a1=None, a2=None), Parameter('A1_0', 1e-6))
    Initial(A2(a1=None, a3=None), Parameter('A2_0', 1.57e-7))
    Initial(A3(a2=None), Parameter('A3_0', 1e-6))    
    
    #Observable('X3_X5_X7', A1(a1=None,a2=None) + A2(a1=1,a3=None)%A1(a2=1) + A1(a1=None,a2=1)%A2(a1=1,a3=1)%A3(a2=1)) ????
    Observable('X2_X4', A2(a1=None,a3=None) + A3(a2=1)%A2(a3=1,a1=None))
    #Observable('X3_2X6_X8_X9', A1(a1=None,a2=None) + 2*???)
    Observable('X1', A3(a2=0))
    #Observable('X2_X5_X8_2X10_X11',???
    
def return_model(model_type):
    Model()
    if model_type=='2_monomers':
        two_monomers()
    if model_type=='3_monomers':
        three_monomers()
    if model_type=='testt':
        test()
    return model