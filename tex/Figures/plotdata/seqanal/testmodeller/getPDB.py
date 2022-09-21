import Bio
from Bio.PDB import PDBList
'''Selecting structures from PDB'''
pdbl = PDBList()
PDBlist2=['6nz7','4we6','4we8','4uo0','6e4x','1qu1']
for i in PDBlist2:
    pdbl.retrieve_pdb_file(i,pdir='PDB',file_format='pdb')
