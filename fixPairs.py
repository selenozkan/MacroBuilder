from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import asarray, array, dot, set_printoptions, append
import numpy
import sys

parser = PDBParser(PERMISSIVE=1, QUIET=1)

# we have to build a second prot_chains_o with consistency and no redundancies.

def checkredcon(chain_names):
    """
    Checks inconsistency (chain-identifier of same structure different in different pairs) or redundancy
    (chain name appears for more than one structure) of chain-identifiers for the input structures. Returns
    a string with 'True' or 'False' value.
    """
    chains = []
    fix = False
    for prot in chain_names:
        if len(chain_names[prot]) != 1:
            fix = True
        for chain in chain_names[prot][0]:
            if chain not in chains:
                chains.append(chain)
            else:
                fix = True
    return (fix)


def get_chain_name(chain, chains):
    """
    Renames inconsistent/redundant chain-identifiers with a letter that has not been used by other pairs.
    If all the letters has been used, exits the program.
    """
    alf = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'
    for let in alf:
        if let != chain and let not in chains:
            return (let)
    for let in alf:
        for let2 in alf:
            if let+let2 != chain and let+let2 not in chains:
                return let+let2
    print('The limit of chain names has been reached for re-naming inconsistent / redundant chain-identifiers. Program terminated.')
    exit()

def build_fixed_prot(prot_chains_o, multi, chain_n):
    """Used for sets of pairs showing inconsistency/redundancy in their chain-identifiers.Assigns chain-identifiers
    to structures in such way that the whole set of pairs with fully consistent and nonredundant chain-identifiers.
    Returns a dictionary in which keys are the list of pairs; and values are dictionaries with elements of each pair
    as keys and corresponding chain-identifiers as values.
    """
    fixed_prot_chains = {}
    new_chain_names = {}
    all_chains = []
    for pair in prot_chains_o:
        chains_pair = {}
        for prot in prot_chains_o[pair]:
            chains = prot_chains_o[pair][prot]
            if prot not in new_chain_names:
                n_chains = []
                for chain in chains:
                    if multi == True:
                        chain = get_chain_name(chain, chain_n)
                        chain_n.append(chain)
                        n_chains.append(chain)
                    else:
                        if chain not in all_chains:
                            all_chains.append(chain)
                            n_chains.append(chain)
                        else:
                            chain = get_chain_name(chain, all_chains)
                            all_chains.append(chain)
                            n_chains.append(chain)
                new_chain_names[prot] = n_chains
            else:
                n_chains = []
                for chain in prot_chains_o[pair][prot]:
                    if chain not in new_chain_names[prot]:
                        n_chains = new_chain_names[prot]
                    else:
                        n_chains = new_chain_names[prot]
            chains_pair[prot] = n_chains
        fixed_prot_chains[pair] = chains_pair
    return (fixed_prot_chains)
