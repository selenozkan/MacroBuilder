from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import asarray, array, dot, set_printoptions, append
import numpy
from elements import ELEMENTS
import math
parser = PDBParser(PERMISSIVE=1, QUIET=1)

def get_model(structure_id):
    """Returns structure object from .pdb file"""
    filename = structure_id + ".pdb"
    structure = parser.get_structure(structure_id, filename)
    return (structure[0])

def get_coordinates(model, chains, mol_type):
    """Returns the coordinates of the backbone of structure. Chooses atoms according to the molecule type."""
    coord = []
    backbone_ids = {"N":["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
                    "P":["CA"],
                    "B": ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "CA"]}
    for chain in model.get_list():
        if chain.get_id() in chains:
            for residue in chain.get_list():
                for atom in residue:
                    if mol_type != '-':
                        if atom.get_id() in backbone_ids[mol_type]:
                            coord.append(list(atom.get_coord()))
                        else:
                            continue
                    else:
                        coord.append(list(atom.get_coord()))
    coord = numpy.array(coord)
    return coord

def get_coordinates_info(model_distance, chains, mol_type):
    """Returns the coordinates of the backbone, end coordiantes of the chains and
    chain ids. The function is used for calculating the distance matrix between
    the previous and newly added structure"""
    coord = []
    chains_span = []
    chains_id = []
    N = 0
    backbone_ids = {"N":["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
                    "P":["CA"],
                    "B": ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "CA"]}
    for chain in model_distance.get_list():
        if chain.get_id() in chains:
            for residue in chain.get_list():
                for atom in residue:
                    if mol_type != '-':
                        if atom.get_id() in backbone_ids[mol_type]:
                            coord.append(list(atom.get_coord()))
                            N = N + 1
                        else:
                            continue
                    else:
                        coord.append(list(atom.get_coord()))
                        N = N + 1
            chains_span.append(N)
            chains_id.append(chain.get_id())
    coord = numpy.array(coord)
    return coord, chains_span, chains_id


def adjust_space(string, n_space, order):
    """Adjusts spacing between columns; used for generating new file in .pdb format"""
    space = ' ' * (n_space - len(string))
    if order == '1':
        return (space + string)
    else:
        return (string + space)

def write_line(resname, residue_id, chain_id, atom_id, atom_ele, coordinates, occupancy, bfactor, n, fileout):
    """Returns line in the format required for .pdb file"""
    residue_id = adjust_space(residue_id, 7, '2')
    atom_id = adjust_space(atom_id, 4, '2')
    coords1 = adjust_space(coordinates[0], 8, '1')
    coords2 = adjust_space(coordinates[1], 8, '1')
    coords3 = adjust_space(coordinates[2], 8, '1')
    coords_string = coords1 + coords2 + coords3
    bfactor = adjust_space(('%.2f')%(bfactor), 6, '1')
    atom_ele = adjust_space(atom_ele, 12, '1')
    string = adjust_space(str(n), 7, '1')
    chain_name = adjust_space(chain_id, 2, '2')
    pdb_line = ('ATOM%s  %s%s %s%s%s  %.2f%s%s  \n') % (
        string, atom_id, resname, chain_name, residue_id, coords_string, occupancy, bfactor, atom_ele)
    fileout.write(pdb_line)

def fix_files(prot_chains_o, fixed_chains, pairs_f, multi, prot_mol_o):
    """Corrects unmatching chain-identifiers with same ones if introduced differently for the same
    structure in the input file. Returns a dictionary of keys as list of pairs and values as correct
    chain-identifiers."""
    message = ""
    Messages = []
    m = 0
    for pair in prot_chains_o:
        m = m + 1
        for prot in prot_chains_o[pair]:
            if prot_chains_o[pair][prot] != fixed_chains[pair][prot]:
                #if multi == True:
                #    outputname = pairs_f[str(pair)].split('.pdb')[0] + str(m) + '.fixed.pdb'
                #else:
                outputname = pairs_f[str(pair)].split('.pdb')[0] + '.fixed.pdb'
                chain_pre = prot_chains_o[pair][prot]
                chain_post = fixed_chains[pair][prot]
                structure_id = pairs_f[str(pair)].split('.pdb')[0]
                structure = parser.get_structure(structure_id, pairs_f[str(pair)])
                model = structure[0]
                output = open(outputname, 'w')
                n = 1
                for chain in model.get_list():
                    for residue in chain.get_list():
                        if chain.id in chain_pre:
                            chain_id = chain_post[chain_pre.index(chain.id)]
                        else:
                            chain_id = chain.id
                        for atom in residue.get_list():
                            coords = atom.get_coord()
                            coords_string_aux = [("%.3f") % (coords[0]), ("%.3f") % (coords[1]), ("%.3f") % (coords[2])]
                            write_line(residue.resname, str(residue.get_id()[1]), chain_id, atom.get_id(),
                                           atom.get_id()[0], coords_string_aux, atom.get_occupancy(),
                                           atom.get_bfactor(), n, output)
                            n = n + 1
                    output.write('TER\n')
                message = ("%s\t\t%s\t\t%s\t\t%s\t\t%s\n") % (
                    (pairs_f[str(pair)].split('.pdb')[0]).split('/')[-1], str(prot), (str(chain_pre)).split('\'')[1],
                    (str(chain_post)).split('\'')[1], (str(outputname)).split('/')[-1])
                output.close()
                pairs_f[str(pair)] = outputname
                Messages.append(message)
                output.close()
    return (pairs_f, Messages)



def RTadd(modely, modelx, rot, tran, outputname, prot_chains_o, pair2, super_prot, prot_mol_o):
    """Creates new pdb file of the base structure; using all the information related to
    superimposition of the new structure. Returns list of chains for new chain id(s)
    unique to newly added structure, the ones in the base and the ones newly added structure.
    """
    output = open(outputname, 'w')
    n = 1
    for chain in modely.get_list():
        for residue in chain.get_list():
            for atom in residue.get_list():
                coords = atom.get_coord()
                new_coords = dot(atom.get_coord(), rot) + tran
                coords_string_aux = [("%.3f") % (new_coords[0]), ("%.3f") % (new_coords[1]), ("%.3f") % (new_coords[2])]
                write_line(residue.resname, str(residue.get_id()[1]), chain.id, atom.get_id(), atom.get_id()[0],
                               coords_string_aux, atom.get_occupancy(), atom.get_bfactor(), n, output)
                n = n + 1
        output.write('TER\n')

    for chain in modelx.get_list():
        if chain.get_id() not in prot_chains_o[str(pair2)][super_prot]:
            for residue in chain.get_list():
                for atom in residue.get_list():
                    coords = atom.get_coord()
                    coords_string_aux = [("%.3f") % (coords[0]), ("%.3f") % (coords[1]), ("%.3f") % (coords[2])]
                    write_line(residue.resname, str(residue.get_id()[1]), chain.id, atom.get_id(), atom.get_id()[0],
                                   coords_string_aux, atom.get_occupancy(), atom.get_bfactor(), n, output)
                n = n + 1
            output.write('TER\n')
    output.write('END\n')
    output.close()
    list_base, list_newpair, chainstotest = [],[],[]
    for chain in modelx.get_list():
        list_base.append(chain.get_id())  # chains of base
    for chain in modely.get_list():
        list_newpair.append(chain.get_id())  # chains of new pair
    for a in list_newpair:
        if a not in list_base:
            chainstotest.append(a)
    return chainstotest, list_base, list_newpair


def get_clashes(model1,threshold=1.0):
    """Checks for clashes by measuring the distance between radius of two atoms"""
    Messages2=[]
    for chain in model1.get_list():
        for residue in chain.get_list():
            for atom in residue.get_list():
                a_x = float(atom.get_coord()[0])
                a_y = float(atom.get_coord()[1])
                a_z = float(atom.get_coord()[2])
                a_elem = ELEMENTS[atom.get_id()[0]]
                a_r = a_elem.vdwrad
                for chain2 in model1.get_list():
                    if chain2.get_id() != chain.get_id():
                        for residue2 in chain2.get_list():
                            for atom2 in residue2.get_list():
                                b_x = float(atom2.get_coord()[0])
                                b_y = float(atom2.get_coord()[1])
                                b_z = float(atom2.get_coord()[2])
                                b_elem = ELEMENTS[atom2.get_id()[0]]
                                b_r = b_elem.vdwrad
                                dist = math.sqrt((a_x-b_x) ** 2 + (a_y-b_y) ** 2 + (a_z-b_z) ** 2)
                                sum_rad = a_r + b_r
                                if dist < sum_rad and abs(dist - sum_rad)>threshold:
                                    message2 = ("%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n") % (
                                        atom.get_id(), residue.resname, chain.get_id(),
                                        atom2.get_id(), residue2.resname, chain2.get_id(),
                                        dist, sum_rad)
                                    Messages2.append(message2)
    return (Messages2)
