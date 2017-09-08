from pysb import *
import re
from pysb.bng import generate_equations
import pandas as pd
import numpy as np


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split('(\d+)', text) ]


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
    values = range(1, n_constants, 2)
    for i, j in enumerate(values):
        params[j-1] = Parameter('k{0}'.format(i+1), 1)
        params[j] = Parameter('l{0}'.format(i+1), 1)

    params[n_constants - 1] = Parameter('k{0}'.format(n), 1)
    return params


def n_rules(n):
    mon = n_monomers_model(n)
    # Initial conditions
    for i in mon:
        state_sites = {j: None for j in i.sites}
        Initial(i(state_sites), Parameter(i.name+'_0', 5))

    pars = n_rate_constants(n)
    for i in range(n):
        if i == 0:
            m1_site = mon[i].sites[1]
            Rule('rule_{0}'.format(i+1), mon[i]({m1_site: None}) + mon[i]({m1_site: None}) <> mon[i]({m1_site: 1}) %
                 mon[i]({m1_site: 1}), pars[i], pars[i+1])

        elif i == n-1:
            m1_site = mon[i].sites[0]
            m2_site = mon[i-1].sites[0]
            Rule('rule_{0}'.format(i+1), mon[i]({m1_site: None}) + mon[i-1]({m2_site: None}) >> mon[i]({m1_site: 1}) %
                 mon[i-1]({m2_site: 1}), pars[-1])

        else:
            m1_site = mon[i].sites[1]
            m2_site = mon[i-1].sites[0]
            Rule('rule_{0}'.format(i+1), mon[i]({m1_site: None}) + mon[i-1]({m2_site: None}) <> mon[i]({m1_site: 1}) %
                 mon[i-1]({m2_site: 1}), pars[2*i], pars[(2*i)+1])
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


def get_groups(df, group, group_idx):
    """

    :param df: Dataframe with species grouped
    :param group: str, can be 'lu' or 'ld'
    :param group_idx: group index to get
    :return:
    """
    if group_idx > len(df.columns):
        raise ValueError('group_idx larger than number of monomers')

    group_fixed_idx = group_idx - 1
    if group == 'ld':
        ld_idx = np.negative(group_idx)
        diagonal = pd.Series(np.diag(df, ld_idx))
        #  removing nans
        # ld_sum = np.sum(diagonal.dropna())
        return diagonal.dropna()
    elif group == 'lu':
        # lu_sum = np.sum(np.diag(df, group_fixed_idx))
        return np.diag(df, group_fixed_idx)
    elif group == 'b':
        # b_sum = np.sum(df[group_idx].dropna())
        return df[group_idx].dropna()
    else:
        raise ValueError('Parameter value not valid')


def get_complex_pattern_ic(model, cp):
    ic = [i[0] for i in model.initial_conditions]
    #  TODO: might be necessary to use is_equivalent()
    idx = [i for i, x in enumerate(ic) if x.is_equivalent_to(cp)]
    if not idx:
        ic_value = 0
    else:
        ic_value = model.initial_conditions[idx[0]][1].value
    return ic_value


def get_total_monomer(df, monomer_idx, model):
    n = df.columns[-1]
    b_ld_total = np.array([np.concatenate([get_groups(df, group='b', group_idx=i),
                                           get_groups(df, group='ld', group_idx=i)]) for i in range(monomer_idx, n+1)])
    if len(b_ld_total) > 1:
        b_ld_total = np.concatenate(b_ld_total)
    elif len(b_ld_total) == 1:
        b_ld_total = b_ld_total[0]

    lu_m1_total = np.array([get_groups(df, group='lu', group_idx=i + 1) for i in range(monomer_idx, n)])
    if len(lu_m1_total) > 1:
        lu_m1_total = np.concatenate(lu_m1_total)
    elif len(lu_m1_total) == 1:
        lu_m1_total = lu_m1_total[0]

    b_ld_value = 0
    if not b_ld_total.size:
        pass
    else:
        for ld in b_ld_total:
            b_ld_value += get_complex_pattern_ic(model, ld)

    lu_value = 0
    if not lu_m1_total.size:
        pass
    else:
        for lu in lu_m1_total:
            lu_value += get_complex_pattern_ic(model, lu)

    total_monomer = [b_ld_total, lu_m1_total]
    total_value = b_ld_value - lu_value
    return total_value


def group_species(model):
    if not model.species:
        generate_equations(model)
    sp_mon = species_mon_names(model.species)
    mons_polymer = {i: [] for i in range(1, len(model.monomers)+1)}
    for i, j in sp_mon.items():
        mon_idx = get_number(i[-1])
        mons_polymer[mon_idx].append(i)

    # sorting lists by length
    for i in mons_polymer:
        mons_polymer[i].sort(key=len)
        species = [sp_mon[j] for j in mons_polymer[i]]
        mons_polymer[i] = species

    # dataframe to group species
    df = pd.DataFrame(np.nan, index=range(1, (2*len(model.monomers))+1), columns=range(1, len(model.monomers)+1))
    for i, j in mons_polymer.items():
        df[i].iloc[:len(j)] = j

    return df


def get_lu_mplus1_discriminant(model, mplus1):
    n = len(model.monomers)
    if mplus1 == 1 or mplus1 >= n:
        raise ValueError('only values between 1 and n-1 are valid')
    k_mplus1 = model.rules[mplus1-1].rate_forward
    l_mplus1 = model.rules[mplus1-1].rate_reverse
