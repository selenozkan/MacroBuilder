###################################################
##                                               ##
##         Readme for MacroBuilder               ##
##                                               ##
###################################################
Last update: 3rd April 2018


Macromolecular complexes builder tool developed by:
- Eva Martín del Pico
- Selen Özkan
- Antía Fernández Pintos
MSc Bioinformatics for Health Sciences, Universitat Pompeu Fabra (UPF), Barcelona

This documentation provides the following information about the program MacroBuilder:
   -Requirements
   -Installation
   -Usage
   -Outputs

*** FOR FURTHER INFORMATION AND EXAMPLES ABOUT HOW TO PREPARE INPUT FILES AND USAGE OF THIS PROGRAM, CONSULT THE USER MANUAL OF MACROBUILDER AT: https://github.com/EvaMart/MacroBuilder/blob/master/manual.mkd ***

-------------------------------------------------------------------------
-------------------------------------------------------------------------

                   **********************
                   **   REQUIREMENTS   **
                   **********************

REQUIREMENTS:

1. Python version 3.x 
   Last Python version (and others), as well as its documentation are available at: https://docs.python.org/3/ .


2. Python packages:

   2.1. Biopython (version 1.70)
   2.2. NumPy (version 1.14.1)
   2.3. SciPy (version 1.0.0)
   2.4. Matplotlib (version 1.5.0)


-------------------------------------------------------------------------
-------------------------------------------------------------------------

                   **********************
                   **   INSTALLATION   **
                   **********************

Unix based systems:

A) As a program:
   Download source code (MacroBuilder_1.0.tar.gz) and uncompress it in the desired directory.
   It is ready to use.
   
A) As a package:
   Download macrobuilder.tar.gz and uncompress the package in the desired directory.
   Within the same directory, open the terminal and run:
   
   -----------------------------------------------------------
   pip install .
   -----------------------------------------------------------
   
   To check if it is correctly installed, open python3 in the terminal and type:
   
   -----------------------------------------------------------
   >>> import macrobuilder
   >>> macrobuilder.test()
   -----------------------------------------------------------
   
   Should return:
   
   -----------------------------------------------------------
   You are using a MacroBuilder function.
   -----------------------------------------------------------

-------------------------------------------------------------------------
-------------------------------------------------------------------------

                   ***************
                   **   USAGE   **
                   ***************
                   
                   
COMMAND LINE ARGUMENTS

MacroBuilder is a command line program. The command for its execution is of the following form:

-------------------------------------------------------------------------------------------
python3 macrobuilder.py -i <pairs_input_file> -f <format> -o <output_name> [<Optional>]
-------------------------------------------------------------------------------------------

where:
<pairs_input_file> is a required argument. Plain text containing details about the pairs from which to build the complex (see examples below).
<format> is a required argument. Format of the 3D molecular structures specified inside the <input file> that will be used to build the complex. MacroBuilder supports PDB (-f pdb) and mmCIF (-f mmcif).
<output_name> is a required argument. Name of the output report. 

Optional arguments [<Optional>]:

   -multi, --multimer (no argument needed): default configuration is false. Select this option if you want to build an homomultimer from one identical homodimer. 

   -mx, --matrix (no argument needed): default value is false. If true, the program will return an output file with rotation and translation matrices. 

   -c, --clashes (argument: <threshold>): by default is false. If true, MacroBuilder will compute the pairwise distances of all atoms in the intermediates and the final complex. If distance between two atoms is greater than the sum of their covalent radii by more than the <threshold>, the "clash" will be reported. <threshold> unit: Amstrongs. 
   


PAIRS INPUT FILE

The input file is a plain text in which each line contains information about one pair to be used in the construction of the complex. Each line is of the form:
-----------------------------------------------------------------------------------------------------------
<struct1_id>  <struct2_id>  <struct1_chain_ids> <struct2_chain_ids> <type_struct1>  <type_struct2>   <file>
-----------------------------------------------------------------------------------------------------------

    <struct1_id>: identifier of one of the structures in the pair.
    <struct2_id>: identifier of the other structure in the pair. 
    <struct1_chain_ids>: identifiers of chain(s) in struct1. If more than 1 chain, must be separated with commas (e.g. A,B,C,K).
    <struct2_chain_ids>: identifiers of chain(s) in struct2. If more than 1 chain, must be separated with commas (e.g. A,B,C,K).
    <type_struct1>: molecular type of struct1. Three posibilities: "P" for peptidic, "N" for nucleotic (DNA, RNA) and "-" for others or mixtures.  
    <type_struct2>: molecular type of struct1. Three posibilities: "P" for peptidic, "N" for nucleotic (DNA, RNA) and "-" for others or mixtures.
    <file>: path of the 3D atomic structure file of the pair. 

   
*** FOR FURTHER INFORMATION AND EXAMPLES ABOUT HOW TO PREPARE INPUT FILES AND USAGE OF THIS PROGRAM, CONSULT THE USER MANUAL OF MACROBUILDER AT: https://github.com/EvaMart/MacroBuilder/blob/master/manual.mkd ***



ATOMIC 3D STRUCTURES OF PAIRS
The user must provide the 3D structure of each pair used in the construction of the complex, whose paths are specified in the pairs input file (explained above). MacroBuilder supports PDB and mmCIF formats.

-------------------------------------------------------------------------
-------------------------------------------------------------------------

                   *****************
                   **   OUTPUTS   **
                   *****************

   1) Report:
      Plain text consisting of at least one block of relevant information about the macrocomplex building process. Report of structural information including RMSD values, min/max distances is the main block of the report. If there are inconsistencies and/or redundancies in the chain identifiers and thus files need to be fixed, a second block containing information relative to this proccess is added before the mentioned one. If the user introduced the -c option, found clashes will be shown as an additional block at the end of the file.
      
    2) Final macrocomplex and intermediates:
       The final macrocomplex as well as all the intermediary complexes generated in its construction (bases) are returned in the same format as the 3D atomic structures provided by the user. The files are called base<N>.<ter>, being <N> the intermediate to which it belong and <ter> the appropriate extension (.pdb or .cif).
       
    3) Fixed atomic 3D structure files:
       3D atomic structure files generated when redundancies and/or inconsistencies are detected. They contain the exact same structure as the provided files, except for the chain identifiers, which are changed to ensure consistency and non-redundancy.
       
    4) Distance Plots for intermediates and final complexes (finalplot.png, baseN.png)
       In each step of the construction, the distance between the previously generated base and the newly added structure is computed and a plot of them is generated. Distance plots are called baseN.png, and N being the intermediate to which it belongs. 
       A plot of the distances between all the atoms in the final macrocomplex,*finalplot.png*, is also generated.
    
    5) HTML results (results.html):
       This html file contains the results of a MacroBuilder run in a more visual form. It includes the details about the job, 3D atomic representations of the final complex and its distance plot.
       In addition, a list of relevant parameters in each pair addition (RMSD of common chain in pair and minimum and maximum distances between old and newly added atoms) and a link to the rpresentations of 3D structure and plot of the resulting intermediate. 
      
   6) Translation and rotation matrices (matrix.trans):
      This file contains, for each base, the rotation and translation matrices that the main algorithm used to position the last added structure. And excerpt from example data showing the rotation (top) and translation (bottom) matrices for one base. 
