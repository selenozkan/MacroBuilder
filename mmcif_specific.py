from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.MMCIFParser
from Bio.PDB import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import asarray, array, dot, set_printoptions, append
import numpy
parser = PDBParser(PERMISSIVE=1, QUIET=1)
from scipy.spatial.distance import cdist
from elements import ELEMENTS
import math

# Creation of parser of mmCIF files (similar to PDB parser)

# This parser returns a Structure Object
parser = Bio.PDB.MMCIFParser()

def get_model(structure_id):
    """Returns a dictionary with all metatags from a file with format PDBx/mmCIF.
    If the metatag has more than one item, it returns a list."""
    filename=structure_id + '.cif'
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(filename)
    return mmcif_dict

def get_coordinates(dict_, chain_name, mol_type):
    """Returns the coordinates of backbone atoms for a given model"""
    coord = []
    backbone_ids = {"N":["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
                    "P":["CA"],
                    "B": ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "CA"]}
    for a in range(len(dict_['_atom_site.Cartn_x'])):
        if dict_['_atom_site.auth_asym_id'][a] in chain_name:
            if mol_type != '-':
                if dict_['_atom_site.label_atom_id'][a] in backbone_ids[mol_type]:
                    coords1 = float(dict_['_atom_site.Cartn_x'][a])
                    coords2 = float(dict_['_atom_site.Cartn_y'][a])
                    coords3 = float(dict_['_atom_site.Cartn_z'][a])
                    coord.append([coords1, coords2, coords3])
                else:
                    continue
            else:
                coords1 = float(dict_['_atom_site.Cartn_x'][a])
                coords2 = float(dict_['_atom_site.Cartn_y'][a])
                coords3 = float(dict_['_atom_site.Cartn_z'][a])
                coord.append([coords1, coords2, coords3])
    coord = numpy.array(coord)
    return coord

def get_coordinates_info(model_distance, chains, mol_type):
    """Returns the coordinates of the backbone, end coordiantes of the chains and
    chain ids. The function is used for calculating the distance matrix between
    the previous and newly added structure"""
    coord = []
    chains_span = []
    chains_id = []
    M = 0
    T = 0
    backbone_ids = {"N":["P", "O5'", "C5'", "C4'", "C3'", "O3'"],
                    "P":["CA"],
                    "B": ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "CA"]}
    for a in range(len(model_distance['_atom_site.Cartn_x'])):
        if model_distance['_atom_site.auth_asym_id'][a] in chains:
            if model_distance['_atom_site.auth_asym_id'][a] not in chains_id:
                chains_id.append(model_distance['_atom_site.auth_asym_id'][a])
                T=T+M
                if M != 0:
                    chains_span.append(T)
                M=0
            if mol_type != '-':
                if model_distance['_atom_site.label_atom_id'][a] in backbone_ids[mol_type]:
                    coords1 = float(model_distance['_atom_site.Cartn_x'][a])
                    coords2 = float(model_distance['_atom_site.Cartn_y'][a])
                    coords3 = float(model_distance['_atom_site.Cartn_z'][a])
                    coord.append([coords1, coords2, coords3])
                    M = M + 1
                else:
                    continue
            else:
                coords1 = float(model_distance['_atom_site.Cartn_x'][a])
                coords2 = float(model_distance['_atom_site.Cartn_y'][a])
                coords3 = float(model_distance['_atom_site.Cartn_z'][a])
                coord.append([coords1, coords2, coords3])
                M = M + 1
    chains_span.append(T+M)
    coord = numpy.array(coord)
    return coord, chains_span, chains_id



def adjust_space(string, n_space, order):
    """Adjusts spacing between columns; used for generating new file in .mmCIF format"""

    space = ' ' * (n_space - len(string))
    if order == '1':
        return (space + string)
    else:
        return (string + space)

def add_header(output):
    """Returns header in the format required for .mmCIF file"""
    output.write("""data_
#
loop_
_atom_type.symbol
C
MG
N
O
P
S
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
""")

def write_line(group_PDB_atom, id_atom, type_symbol_atom, label_atom_id_atom, label_comp_id_atom,
                     label_asym_id_atom, label_entity_id_atom, label_seq_atom, occupancy_atom, B_iso_equiv_atom,
                     auth_seq_id_atom, auth_comp_id_atom, auth_asym_id_atom, auth_atom_id_atom, PDB_model_atom,
                     coords_string_aux, n, fileout):
    """Returns line in the format required for .mmCIF file"""
    group_PDB = adjust_space(group_PDB_atom, 7, '2')
    string = adjust_space(str(n), 6, '2')
    symbol_atom = adjust_space(type_symbol_atom, 3, '2')
    quotids = ["O5'", "C5'", "C4'", "C3'", "O3'", "O4'", "C2'", "O2'", "C1'"]
    if label_atom_id_atom in quotids:
        labelatomid = adjust_space('"'+label_atom_id_atom+ '"', 6, '2')
        authatomidatom = adjust_space('"'+auth_atom_id_atom+ '"', 6, '2')
    else:
        labelatomid = adjust_space(label_atom_id_atom, 6, '2')
        authatomidatom = adjust_space(auth_atom_id_atom, 6, '2')
    labelcompid = adjust_space(label_comp_id_atom, 4, '2')
    labelasymid = adjust_space(label_asym_id_atom, 3, '2')
    labelentityid = adjust_space(label_entity_id_atom, 2, '2')
    labelseq = adjust_space(label_seq_atom, 4, '2')
    coords1 = adjust_space(coords_string_aux[0], 8, '2')
    coords2 = adjust_space(coords_string_aux[1], 8, '2')
    coords3 = adjust_space(coords_string_aux[2], 8, '2')
    coords_string = coords1 + coords2 + coords3
    occupancy_atom = adjust_space(str(occupancy_atom), 5, '2')
    B_iso_equiv_atom = "{0:.2f}".format(float(B_iso_equiv_atom))
    B_isoequivatom = adjust_space(str(B_iso_equiv_atom), 7, '2')
    authseqid = adjust_space(auth_seq_id_atom, 5, '2')
    authcompid = adjust_space(auth_comp_id_atom, 4, '2')
    authasym = adjust_space(auth_asym_id_atom, 3, '2')
    PDBmodelatom = adjust_space(PDB_model_atom, 2, '2')
    mmcif_line = ('%s%s%s%s. %s%s%s%s? %s%s%s? %s%s%s%s%s\n') % (
        group_PDB, string, symbol_atom, labelatomid, labelcompid, labelasymid, labelentityid,
        labelseq, coords_string, occupancy_atom, B_isoequivatom, authseqid, authcompid, authasym,
        authatomidatom, PDBmodelatom)
    fileout.write(mmcif_line)

def fix_files(prot_chains_o, fixed_chains, pairs_f, multi, prot_mol_o):
    """Corrects unmatching chain-identifiers with same ones if introduced differently for the same
    structure in the input file. Depending on the molecule types, adds info to the header.Returns
    a dictionary of keys as list of pairs and values as correct chain-identifiers."""
    message = ""
    Messages = []
    m = 0
    for pair in prot_chains_o:
        m = m + 1
        for prot in prot_chains_o[pair]:
            if prot_chains_o[pair][prot] != fixed_chains[pair][prot]:
                chain_pre = prot_chains_o[pair][prot]
                #if multi == True:
                #    outputname = pairs_f[str(pair)].split('.cif')[0] + str(m) + '.fixed.cif'
                #else:
                outputname = pairs_f[str(pair)].split('.cif')[0] + '.fixed.cif'
                chain_post = fixed_chains[pair][prot]
                structure_id = pairs_f[str(pair)].split('.cif')[0]
                model = get_model(structure_id)
                output = open(outputname, 'w')
                if prot_mol_o[str(pair)][prot] == 'P':
                    add_header(output)
                elif prot_mol_o[str(pair)][prot] =='N':
                    nuc.add_header(output, model)
                n = 1
                for a in range(len(model['_atom_site.Cartn_x'])):
                    if model['_atom_site.auth_asym_id'][a] in chain_pre:
                        chain_id = chain_post[chain_pre.index(model['_atom_site.auth_asym_id'][a])]
                    else:
                        chain_id = model['_atom_site.auth_asym_id'][a]
                    coords1 = ("%.3f") % (float(model['_atom_site.Cartn_x'][a]))
                    coords2 = ("%.3f") % (float(model['_atom_site.Cartn_y'][a]))
                    coords3 = ("%.3f") % (float(model['_atom_site.Cartn_z'][a]))
                    coords_string_aux = [coords1, coords2, coords3]
                    write_line(model['_atom_site.group_PDB'][a], model['_atom_site.id'][a],
                         model['_atom_site.type_symbol'][a], model['_atom_site.label_atom_id'][a],
                         model['_atom_site.label_comp_id'][a], model['_atom_site.label_asym_id'][a],
                         model['_atom_site.label_entity_id'][a], model['_atom_site.label_seq_id'][a],
                         model['_atom_site.occupancy'][a], model['_atom_site.B_iso_or_equiv'][a],
                         model['_atom_site.auth_seq_id'][a], model['_atom_site.auth_comp_id'][a],
                         chain_id, model['_atom_site.auth_atom_id'][a],
                         model['_atom_site.pdbx_PDB_model_num'][a], coords_string_aux, n, output)
                    n = n + 1
                message = ("%s\t\t%s\t\t%s\t\t%s\t\t%s\n") % (
                    (pairs_f[str(pair)].split('.pdb')[0]).split('/')[-1], str(prot), (str(chain_pre)).split('\'')[1],
                    (str(chain_post)).split('\'')[1], (str(outputname)).split('/')[-1])
                output.close()
                pairs_f[str(pair)] = outputname
                Messages.append(message)
                output.close()
    return (pairs_f, Messages)


def RTadd(dict_y, dict_x, rot, tran, outputname, prot_chains_o, pair2, super_prot, prot_mol_o):
    """Creates new mmcif file of the base structure; using all the information related to
    superimposition of the new structure. Returns list of chains for new chain id(s)
    unique to newly added structure, the ones in the base and the ones newly added structure.
    """
    output = open(outputname, 'w')
    if prot_mol_o[str(pair2)][super_prot] == 'P':
        add_header(output)
    elif prot_mol_o[str(pair2)][super_prot] =='N':
        nuc.add_header(output, dict_x)
    list_base, list_newpair, chainstotest = [],[],[]     # chains of base, new pairs and unique to newpair
    n = 1
    for a in range(len(dict_y['_atom_site.Cartn_x'])):
        coords_string = [float(dict_y['_atom_site.Cartn_x'][a]), float(dict_y['_atom_site.Cartn_y'][a]) , float(dict_y['_atom_site.Cartn_z'][a])]
        new_coords = dot(coords_string, rot) + tran
        coords1 = ("%.3f") % (new_coords[0])
        coords2 = ("%.3f") % (new_coords[1])
        coords3 = ("%.3f") % (new_coords[2])
        coords_string_aux = [coords1, coords2, coords3]
        write_line(dict_y['_atom_site.group_PDB'][a], dict_y['_atom_site.id'][a],
                     dict_y['_atom_site.type_symbol'][a], dict_y['_atom_site.label_atom_id'][a],
                     dict_y['_atom_site.label_comp_id'][a], dict_y['_atom_site.label_asym_id'][a],
                     dict_y['_atom_site.label_entity_id'][a], dict_y['_atom_site.label_seq_id'][a],
                     dict_y['_atom_site.occupancy'][a], dict_y['_atom_site.B_iso_or_equiv'][a],
                     dict_y['_atom_site.auth_seq_id'][a], dict_y['_atom_site.auth_comp_id'][a],
                     dict_y['_atom_site.auth_asym_id'][a], dict_y['_atom_site.auth_atom_id'][a],
                     dict_y['_atom_site.pdbx_PDB_model_num'][a], coords_string_aux, n, output)
        if dict_y['_atom_site.auth_asym_id'][a] not in list_newpair:
            list_newpair.append(dict_y['_atom_site.auth_asym_id'][a])
        n = n + 1

    for b in range(len(dict_x['_atom_site.Cartn_x'])):
        if dict_x['_atom_site.auth_asym_id'][b] not in prot_chains_o[str(pair2)][super_prot]:
            coords1 = ("%.3f") % (float(dict_x['_atom_site.Cartn_x'][b]))
            coords2 = ("%.3f") % (float(dict_x['_atom_site.Cartn_y'][b]))
            coords3 = ("%.3f") % (float(dict_x['_atom_site.Cartn_z'][b]))
            coords_string_aux = [coords1, coords2, coords3]
            write_line(dict_x['_atom_site.group_PDB'][b], dict_x['_atom_site.id'][b],
                         dict_x['_atom_site.type_symbol'][b], dict_x['_atom_site.label_atom_id'][b],
                         dict_x['_atom_site.label_comp_id'][b], dict_x['_atom_site.label_asym_id'][b],
                         dict_x['_atom_site.label_entity_id'][b], dict_x['_atom_site.label_seq_id'][b],
                         dict_x['_atom_site.occupancy'][b], dict_x['_atom_site.B_iso_or_equiv'][b],
                         dict_x['_atom_site.auth_seq_id'][b], dict_x['_atom_site.auth_comp_id'][b],
                         dict_x['_atom_site.auth_asym_id'][b], dict_x['_atom_site.auth_atom_id'][b],
                         dict_x['_atom_site.pdbx_PDB_model_num'][b], coords_string_aux, n, output)
            n = n + 1
        if dict_x['_atom_site.auth_asym_id'][b] not in list_base:
            list_base.append(dict_x['_atom_site.auth_asym_id'][b])

    output.close()
    for c in list_newpair:
        if c not in list_base:
            chainstotest.append(c)
    return chainstotest, list_base, list_newpair

def get_clashes(dict_1, threshold=1.0):
    """Checks for clashes by measuring the distance between radius of two atoms"""
    Message2 = []
    for a in range(len(dict_1['_atom_site.Cartn_x'])):
        a_x = float(dict_1['_atom_site.Cartn_x'][a])
        a_y = float(dict_1['_atom_site.Cartn_y'][a])
        a_z = float(dict_1['_atom_site.Cartn_z'][a])
        a_elem = ELEMENTS[dict_1['_atom_site.type_symbol'][a]]
        a_r = a_elem.vdwrad
        for b in range(len(dict_1['_atom_site.Cartn_x'])):
            if dict_1['_atom_site.auth_asym_id'][a] != dict_1['_atom_site.auth_asym_id'][b]:
                b_x = float(dict_1['_atom_site.Cartn_x'][b])
                b_y = float(dict_1['_atom_site.Cartn_y'][b])
                b_z = float(dict_1['_atom_site.Cartn_z'][b])
                b_elem = ELEMENTS[dict_1['_atom_site.type_symbol'][b]]
                b_r = b_elem.vdwrad
                dist = math.sqrt((a_x-b_x)**2 + (a_y-b_y)**2 + (a_z-b_z)**2)
                sum_rad = a_r + b_r
                if dist < sum_rad and abs(dist - sum_rad)>threshold:
                    message2 = ("%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n") % (dict_1['_atom_site.type_symbol'][a],
                                                                           dict_1['_atom_site.label_comp_id'][a],
                                                                           dict_1['_atom_site.auth_asym_id'][a],
                                                                           dict_1['_atom_site.type_symbol'][b],
                                                                           dict_1['_atom_site.label_comp_id'][b],
                                                                           dict_1['_atom_site.auth_asym_id'][b],
                                                                           dist, sum_rad)
                    Message2.append(message2)
    return (Message2)
