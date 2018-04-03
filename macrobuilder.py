from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import asarray, array, dot, set_printoptions, append
import numpy
import sys
import argparse
import fixPairs as fix
from scipy.spatial.distance import cdist
import pylab
import html_templs as ht
import os
from pylab import rcParams
import pkg_resources

def test():
    print('You are using a MacroBuilder function.')

def parsePairFile(filename):
    """
    Parses the pairs input file and returns the following data structures:
    - pairs_files: dictionary in which keys are the list of pairs and values are the file names.
    - prot_chains_o: dictionary in which keys are list of pairs; values are dictionary of each element of the pair as keys and their chains as values.
    - tuple of pairs.
    - elements of the pairs.
    - chain_names: dictionary in which keys are elements of pairs and values are list of chains.
    """
    try:
        pairfile = open(filename, 'r')
    except IOError as e:
        print('Input Error message:', os.strerror(e.errno))
    pairs_files = {}
    prot_chains_o = {}
    prot_mol_o = {}
    pairs = []
    prots = []
    chain_names = {}
    for line in pairfile:
        fields = line.split('\t')
        pairs.append([fields[0], fields[1]])
        chains_pair = {}
        moltype_pair = {}
        for prot in pairs[-1]:
            if prot not in chain_names:
                chain_names[prot] = []
            chains = fields[2 + pairs[-1].index(prot)].split(',')
            chains_pair[prot] = chains
            if chains not in chain_names[prot]:
                chain_names[prot].append(chains)
            if prot not in prots:
                prots.append(prot)
            mol = fields[4 + pairs[-1].index(prot)]
            moltype_pair[prot] = mol
        prot_chains_o[str(pairs[-1])] = chains_pair
        prot_mol_o[str(pairs[-1])] = moltype_pair
        pairs_files[str(pairs[-1])] = fields[6].split('\n')[0]
    pairfile.close()
    return (pairs_files, prot_chains_o, tuple(pairs), prots, chain_names, prot_mol_o)

def get_ordered(pairs):
    """Returns a list in which elements are list of pairs in the order to be superimposed"""
    pairs_ = list(pairs)
    current = pairs_[0]
    ordered_pairs = [current]
    Pairs = pairs_
    Pairs.remove(current)
    base_all = []
    base_all.append(current[0])
    base_all.append(current[1])
    while(len(Pairs)!=0):
        for pair in Pairs:
            for prot in pair:
                if prot in base_all:
                    base_all.append(pair[0])
                    base_all.append(pair[1])
                    Pairs.remove(pair)
                    ordered_pairs.append(pair)
                    break
    return (ordered_pairs)


def get_pair_prot(pairs, Prots):
    for pair in pairs:
        for prot in pair:
            if prot in Prots:
                return(pair)

def get_rot_tran(y, x):
    """Returns rotation, translation and RMDS values of the superimposed atoms."""
    sup = SVDSuperimposer()
    sup.set(x, y)  # AC over AD
    sup.run()
    rms = sup.get_rms()
    rot, tran = sup.get_rotran()
    return (rot, tran, rms)

def distance_matrix(distmodel, chainstotest, list_base, list_newpair, mol_type):
    """Returns a distance matrix that is constructed for the base and newly added structure."""

    list_uniq = []
    for each in list_newpair:
        if each in chainstotest:
            list_uniq.append(each)

    coords1, c_info_base, c_info_base_id = ftype.get_coordinates_info(distmodel, list_base, mol_type)
    coords2, c_info_uniq, c_info_uniq_id = ftype.get_coordinates_info(distmodel, list_uniq, mol_type)
    distances = cdist(coords1, coords2)
    return distances, c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq

def plot_distances(distance, c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq, base_id, idy, out_id):
    """Plots the distance map for the chains of two structures."""
    pylab.matshow(numpy.transpose(distance))
    pylab.colorbar()
    pylab.gca().xaxis.tick_bottom()
    pylab.xlabel(base_id)
    pylab.ylabel(idy)
    # base axis
    x_0 = 0.01
    y_0 = c_info_uniq[-1] + c_info_uniq[-1]/15.0
    m = len(c_info_base_id)
    last = 0
    l = 0
    for a in range(m):
        if c_info_base[a] == last:
            continue
        else:
            x_l = float(c_info_base[a])
            if a == 0:
                x_t = c_info_base[a]/2.3
            else:
                x_t = c_info_base[a-1] + (c_info_base[a] - c_info_base[a-1])/2.3
            y_t = y_0
            pylab.axvline(linewidth=1, color='black', x = x_l)
            pylab.text(x_t, y_t,str(c_info_base_id[a]), fontsize=12)
            last = c_info_base[a]
    n = len(c_info_uniq_id)
    last = 0
    x_0 = c_info_base[-1] + c_info_base[-1]/40.0
    y_0 = 1
    for e in range(n):
        if c_info_uniq == last:
            continue
        else:
            if e == 0:
                y_t = c_info_uniq[e]/2.6
            else:
                y_t = c_info_uniq[e-1] + (c_info_uniq[e] - c_info_uniq[e-1])/2.6
            y_l = float(c_info_uniq[e])
            pylab.axhline(linewidth=1, color='black', y = y_l)
            pylab.text(x_0, y_t,str(c_info_uniq_id[e]), fontsize=12)
            last = c_info_uniq[e]
    pylab.xlim(0, c_info_base[-1])
    pylab.ylim(0, c_info_uniq[-1])
    rcParams['figure.figsize'] = 5, 10
    pylab.savefig(out_id+'.png', dpi=300)
    pylab.close()

def get_file(filein):
    string = ''
    f_pairs = open(filein, "r")
    for line in f_pairs:
        string = string + line
    return string

parser = PDBParser(PERMISSIVE=1, QUIET=1)

if __name__ == "__main__":
    parser2 = argparse.ArgumentParser()
    parser2.add_argument('-i', '--input', dest="infile", action="store", nargs="?", default='.',
                         help="Input FASTA file", required = True)
    parser2.add_argument('-o', '--output', dest="outfile", action="store", default=None, required = True,
                        help="Ouput file. If not defined, it prints output to standard output.")
    parser2.add_argument('-f', '--format', dest="Format", action="store", required = True, type=str,
                         help="Input and output files format")
    parser2.add_argument('-multi', '--multimere', dest = "multi", action="store_true", default = False,
                         help="To build a multimer from a dimer (one file). See documentation for help on how to prepare the --input for this purpose.")
    parser2.add_argument('-c', '--clashes', dest="clash_thres", action="store", default=False,
                         help='''To compute the distance of all the atoms in the itermediary
                         complexes and return pairs of atoms whose distance is smaller
                         than the sum of their covalent radius by a difference bigger
                         than the threshold. -c followed by <threshold>.''')
    parser2.add_argument('-mx', '--matrix', dest="matrix", action="store_true", default=False,
                         help="Rotation and translation matrices")

    options = parser2.parse_args()
    # depending on the format; import related modules
    global ftype
    if options.Format=='pdb':
        import pdb_specific as ftype
        ter = '.pdb'
    elif options.Format=='mmcif':
        import mmcif_specific as ftype
        ter = '.cif'
    else:
       print('Format type '+ options.Format + ' not supported. Please see --help for format options.')
       exit()

    output = open(options.outfile, "w")

    output.write("""###################################################\n
SUPERIMPOSING By PROTEIN STRUCTURE version 1.0\n
###################################################
by Antía Fernández Pintos, Eva Martín del Pico, Selen Ozkan
Universitat Pompeu Fabra, Bioinformatics For Health Sciences, 2018\n\n""")
    output.write("File Type:\t%s\n" % ter)

    pairsfilename = options.infile
    (pairs_f, prot_chains_o, pairs, prots, chain_names, prot_mol_o) = parsePairFile(pairsfilename)

    # check if molecule types are consistent for each element of the pair
    check_all = {}
    for value in prot_mol_o.values():
        for k, v in value.items():
            if k in check_all.keys():
                if not v == check_all[k]:
                    print(
                        "\n!!! Inconsistency detected for the molecule type of %s !!!\nPlease check input file to make sure molecule type is same in different pairs.\nProgram terminated.\n" % k)
                    exit()
            else:
                check_all[k] = v

    num_pairs = len(pairs_f)
    output.write("No. of Pairs:\t%s\n\n" %num_pairs)

    fixing = fix.checkredcon(chain_names)
    if fixing == True:
        chain_n = []
        for a in chain_names:
            for e in chain_names[a]:
                for chain in e:
                    chain_n.append(chain)
        fixed_chains = fix.build_fixed_prot(prot_chains_o, options.multi, chain_n)
        if options.Format == 'pdb':
            for value in fixed_chains.values():
                for k,v in value.items():
                    for elements in v:
                        if len(elements) > 1:
                            print('The limit of chain names has been reached for re-naming inconsistent / redundant chain-identifiers. Program terminated.')
                            exit()
        (pairs_f, Messages) = ftype.fix_files(prot_chains_o, fixed_chains, pairs_f, options.multi, prot_mol_o)
        prot_chains_o = fixed_chains

        output.write("# FIXED PAIRS\n\nFILE\t\t\tFIXED_PAIR\tCHAIN_ID\tNEW_CHAIN_ID\tNEWFILE\n")
        for m in Messages:
            output.write(m)
        print(
            "\n!!! Warning: inconsistency / redundancy detected. Chain names are adjusted to overcome the problem. See the output file for details.\n\n")

    output.write("\n# REPORT\n\nFILE\t\tRMSD\t\tMIN.DIS\t\tMAX.DIS\t\tSUPERIMPOSED\n")

    # html ouput file:
    htmlout = open('results.html', 'w')
    ngl_p = pkg_resources.resource_filename('macrobuilder', 'ngl.js')
    ngl = get_file(ngl_p)
    htmlout.write(ht.header_html.format(ngl=ngl))
    args = ''
    for i in sys.argv[1:]:
        args = args + ' ' + i
    job = ht.job_html.format(arguments = args, pairs_file=get_file(options.infile), n_pairs = num_pairs)
    htmlout.write(job)
    html_bases = ''
    html_table = ht.table_header

    if options.matrix:
        matrix_output = open("matrix.trans", "w")
        matrix_output.write("ROTATION and TRANSLATION MATRICES of SUPERIMPOSED PROTEINS\n")

    # iterating the superimposition for all pairs
    orpairs = get_ordered(pairs)
    pair1 = orpairs[0]
    base = pair1
    base_files = {}
    filex = pairs_f[str(pair1)]
    idx = filex.split(ter)[0]
    base_id = idx
    n = 1
    for pair2 in orpairs[1:]:
        if n>1:
            pair1 = get_pair_prot(orpairs, prot_chains_o[str(pair2)])
        print("Starting on new pairs...\n")
        super_prot = set(pair1).intersection(pair2).pop()
        # check if the superimposed pairs are the same molecule type
        if prot_mol_o[str(pair1)][super_prot] != prot_mol_o[str(pair2)][super_prot]:
            print("!!! Cannot perform superimposition - superimposed elements must be same molecule type. Please check input file to make sure molecule types are written correctly.\nProgram terminated.\n")
            exit()
        filex = pairs_f[str(pair1)]
        idx = filex.split(ter)[0]
        filey = pairs_f[str(pair2)]
        idy = filey.split(ter)[0]
        modely = ftype.get_model(idy)
        modelbase = ftype.get_model(base_id)
        mol_type = prot_mol_o[str(pair1)][super_prot]
        base_id = 'base' + str(n)
        (rot, tran, rms) = get_rot_tran(ftype.get_coordinates(modely,prot_chains_o[str(pair1)][super_prot], mol_type), ftype.get_coordinates(modelbase, prot_chains_o[str(pair2)][super_prot], mol_type))
        chainstotest, list_base, list_newpair = ftype.RTadd(modely, modelbase, rot, tran, base_id + ter, prot_chains_o, pair2, super_prot, prot_mol_o)
        print("Rotation, translation calculations completed.\n")
        print("RMSD value for the pair %s - %s is calculated.\n" % (str(idx).split('/')[-1], str(idy).split('/')[-1]))

        # Distance between the atoms of base and newly introduced pair
        model_distance = ftype.get_model(base_id)
        distance,c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq = distance_matrix(model_distance, chainstotest, list_base, list_newpair, "B")
        # Plotting the distance matrix
        plot_distances(distance, c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq, base_id, idy, base_id)

        output.write("%s\t\t%.3f\t\t%.3f\t\t%.3f\t\t%s - %s\n" % (
        base_id, rms, numpy.min(distance), numpy.max(distance), str(idx).split('/')[-1], str(idy).split('/')[-1]))

        print("Minimum and maximum distance between the previous structure and the added structure are calculated.\n\n")
        print(
            "Superimposition for the pair is complete. Files related to %s are created. Go to current directory to see the file.\n" % base_id)
        print("#########\n")

        if options.matrix:
            matrix_output.write("\n\n# %s \n" % (str(idy).split('/')[-1]))
            for each in rot:
                matrix_output.write("%s\n" % each)
            matrix_output.write("\n%s\n\n" % tran)

        html_base_file = open(base_id + '.html', 'w')
        html_base_file.write(ht.header_html.format(ngl=ngl))
        html_base_file.write(ht.base_html.format(pairs="<b>%s</b> (in base) and <b>%s</b>"%(pair1, pair2),complex_name=base_id, RMSD= "%.8f"%rms, min_dist= "%.7f"%numpy.min(distance),
                                                      max_dist="%.7f"%numpy.max(distance),
                                                      n=n,
                                                      structure_f= 'base'+str(n)+ter,
                                                      plot_f= 'base'+str(n)+ '.png',
                                                      chain_names_base=c_info_base_id))
        html_base_file.write(ht.footer_html)
        html_base_file.close()
        entry = ht.table_entry.format(complex_name=base_id, pair_1 = pair1, pair_2 = pair2, RMSD = "%.8f"%rms, min_dist= "%.7f"%numpy.min(distance), max_dist= "%.7f"%numpy.max(distance))
        html_table = html_table + entry


        n = n + 1

    if options.clash_thres:
        print("Checking for clashes..\n\n")
        output.write(("\n\n# CLASHES\nThreshold: %s\n\nATOM\tRES\tCHAIN\tATOM\tRES\tCHAIN\tDIST\tSUM_OF_RADII\n") % options.clash_thres)
        Messages2 = ftype.get_clashes(ftype.get_model(base_id), float(options.clash_thres))
        for m2 in Messages2:
            output.write(m2)
    modely = ftype.get_model(idy)
    all_chains = list_base + list_newpair
    if mol_type in ['N','P']:
        mol_type_d = 'B'
    else:
        mol_type_d = '-'
    distance, c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq = distance_matrix(model_distance, all_chains, all_chains, all_chains, mol_type_d)
    plot_distances(distance, c_info_base_id, c_info_base, c_info_uniq_id, c_info_uniq, base_id, base_id, 'finalplot')
    html_result = ht.final_html.format(structure_file=base_id+ter, plot_file='finalplot.png')
    htmlout.write(html_result)
    html_table = html_table + ht.table_bot
    htmlout.write(html_table)
    htmlout.write(ht.footer_html)

    output.close()
    if options.matrix:
        matrix_output.close()
    print(
        "To see a more advanced version of the report, please check \'results.html\' file that is created in your working directory.\n")
