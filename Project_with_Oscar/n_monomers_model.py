from pysb import *
import re
from pysb.bng import generate_equations
import pandas as pd
import numpy as np

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def n_monomers_model(n):
    monomers = ['b_{0}'.format(i) for i in range(1, n+1)]
    mons = [0] * n
    for i, j in enumerate(monomers):
        if i == len(monomers)-1:
            mons[i] = Monomer(j, ["s{0}".format(1)])
        else:
            mons[i] = Monomer(j, ["s{0}".format(1), "s{0}".format(2)])
    return mons


def n_rate_constants(n):
    n_constants = (n-1) * 2 + 1
    params = [0] * n_constants
    for i in range(n_constants):
        params[i] = Parameter('k{0}'.format(i), 1)
    return params


def n_rules(n):
    mon = n_monomers_model(n)
    # Initial conditiond
    for i in mon:
        state_sites = {j: None for j in i.sites}
        Initial(i(state_sites), Parameter(i.name+'_0'))

    pars = n_rate_constants(n)
    for i in range(n):
        if i == 0:
            m1_site = mon[i].sites[1]
            Rule('rule_{0}'.format(i), mon[i]({m1_site: None}) + mon[i]({m1_site: None}) <> mon[i]({m1_site: 1}) %
                 mon[i]({m1_site: 1}), pars[i], pars[i+1])

        elif i == n-1:
            m1_site = mon[i].sites[0]
            m2_site = mon[i-1].sites[0]
            Rule('rule_{0}'.format(i), mon[i]({m1_site: None}) + mon[i-1]({m2_site: None}) >> mon[i]({m1_site: 1}) %
                 mon[i-1]({m2_site: 1}), pars[-1])

        else:
            m1_site = mon[i].sites[1]
            m2_site = mon[i-1].sites[0]
            Rule('rule_{0}'.format(i), mon[i]({m1_site: None}) + mon[i-1]({m2_site: None}) <> mon[i]({m1_site: 1}) %
                 mon[i-1]({m2_site: 1}), pars[i+1], pars[i+2])
    return


def generate_model(n):
    Model()
    n_rules(n)
    return model


def species_mon_names(species):
    """

    :param species: list of pysb species
    :return: dictionary of species with their monomers names
    """
    sp_dict = {}
    for i in species:
        bla = tuple([j.monomer.name for j in i.monomer_patterns])
        sp_dict[bla] = i
    return sp_dict


def get_number(string):
    bla = [int(s) for s in re.findall(r'\d+', string)]
    if len(bla) == 1:
        return bla[0]
    else:
        return bla


def get_b_lu_ld(df, group='lu', group_idx=1):
    if group_idx > len(df.columns):
        raise ValueError('group_idx larger than number of monomers')

    group_fixed_idx = group_idx - 1
    if group == 'ld':
        df = df.iloc[::-1]
        return np.diag(df, group_fixed_idx)
    elif group == 'lu':
        return np.diag(df, group_fixed_idx)
    elif group == 'b':
        return df[group_idx]
    else:
        raise ValueError('Parameter value not valid')


def group_species(model):
    if not model.species:
        generate_equations(model)
    sp_mon = species_mon_names(model.species)
    mons_polymer = {i: [] for i in range(1, len(model.monomers)+1)}
    for i, j in sp_mon.items():
        mon_idx = get_number(i[-1])
        mons_polymer[mon_idx].append(j)

    # sorting lists by length
    for i in mons_polymer:
        mons_polymer[i].sort(key=len)
        species = [sp_mon[j] for j in mons_polymer[i]]
        mons_polymer[i] = species

    # dataframe to group species
    df = pd.DataFrame(np.nan, index=range(1, (2*len(model.monomers))+1), columns=range(1, len(model.monomers)+1))
    for i, j in mons_polymer.items():
        df[i].iloc[:len(j)] = j