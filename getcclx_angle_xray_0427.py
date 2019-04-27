import sys,os
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB import NeighborSearch, PDBParser, Selection
import Bio.KDTree
import numpy as np
from Bio.PDB import *
import warnings

import Bio.PDB.StructureBuilder as StructureBuilder
from Bio.PDB.Residue import DisorderedResidue
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)
## Fang-Yu Lin Feb, 2017
## For PDB survey of halogen C-X..Y angle
## usage: python name.py

def calangle(c, x, p):
	coord_c = c
	coord_x = x
	coord_p = p
	a = np.array(coord_c)
	b = np.array(coord_x)
	c = np.array(coord_p)
	ba = a - b
	bc = c - b
	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle_pre = np.arccos(cosine_angle)
        return np.degrees(angle_pre)

def caldist(c, x, p):
	coord_c = c
	coord_x = x
	coord_p = p
	a = np.array(coord_c)
	b = np.array(coord_x)
	c = np.array(coord_p)
	bc = c - b
	return  np.linalg.norm(bc)

surr_residue_dist = 3.5
cx_bond_dist = 2.00  ## this is just for quickly identifying the C-X bond, should be larger than what we've known. 
halogen='CL'   ## Need to be capital
## open the log
sys.stdout = open('ccl_angle.log', 'w')
## open data
dfoh = open('./results/detail_cl_oh_angle.dat', 'w')
dfsh = open('./results/detail_cl_hs_angle.dat', 'w')
dfn_bb = open('./results/detail_cl_nbb_angle.dat', 'w')
dfnh = open('./results/detail_cl_nh_angle.dat', 'w')
dfcarbonyl = open('./results/detail_cl_carbonyl_angle.dat', 'w')
dfcarbonyl2_sidechain = open('./results/detail_cl_carbonylsidechaine_angle.dat', 'w')
dfcx_info = open('info_cx.dat', 'w')
dfcount_info = open('info_count.dat', 'w')
dffail_info = open('info_failed.dat', 'w')
## record the number of hits
number_all=0
number_cxo=0         ## o on side chaine
number_cxns=0        ## n on side chaine
number_cxnbb=0       ## nh on bacbone
number_cxs=0         ## sulfor
number_cxco=0        ## carbonyl
number_cxsideco=0    ## carbonyl on ASN, GLN 
number_fail=0

pdbfile = open('./input/ccl_ethane_benzene_res3.0.sheet', 'r')     
temp = pdbfile.readline()
list1 = temp.split()

for name in list1:
	pdbl = PDBList()
	name2=name.lower()
	pdb1 ='pdb'+name2+".ent"
	number_all=number_all+1
	try:
		pdbl.retrieve_pdb_file(name, pdir=".")
		parser=PDBParser(PERMISSIVE=1)
		structure = parser.get_structure(name,pdb1)
		for model in structure:
			for chain in model:
				for residue in chain:
					if residue.id[0] != ' ' and residue.id[0] != 'W':  # non-water HETATM record
						for atom in residue:
							atom_list_ligand = Selection.unfold_entities(residue, 'A')
							atom_list_protein = Selection.unfold_entities(model, 'A')
							if atom.element == halogen:
								center = atom.get_coord()  # center or halogen atom
								ns_ligand = NeighborSearch(atom_list_ligand)
								ns_protein = NeighborSearch(atom_list_protein)
								dfcx_info.write(str(ns_ligand.search(center, cx_bond_dist, "A"))+"\n")
								if len(ns_ligand.search(center, cx_bond_dist, "A"))==2:
									neighbors_temp = ns_ligand.search(center, cx_bond_dist, "A")[0] #  Cl , C-Cl,Br <=1.93
									if neighbors_temp.element == 'C':
										neighbors_c = ns_ligand.search(center, cx_bond_dist, "A")[0] #  Cl , C-Cl,Br <=1.93
										neighbors_x = ns_ligand.search(center, cx_bond_dist, "A")[1] #  C  , C-Cl,Br <=1.93
									else:
										if neighbors_temp.element == halogen:
											neighbors_x = ns_ligand.search(center, cx_bond_dist, "A")[0] #Cl , C-Cl,Br <=1.93
											neighbors_c = ns_ligand.search(center, cx_bond_dist, "A")[1] #C , C-Cl,Br <=1.93
								
									neighbors_p_array = ns_protein.search(center, surr_residue_dist, "A") #  C  , C-Cl,Br <=1.93
									for neighbors_p in neighbors_p_array:

										##  OH related -  SER, THR, TYR 
										if neighbors_p.name == 'OG' or neighbors_p.name== 'OG1' or neighbors_p.name== 'OH':
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfoh.write( str(angle)+ "\t"+ str(dist_xp)+"\t" +str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")
											number_cxo=number_cxo+1 

										##  HS related -  CYS 
										if neighbors_p.name == 'SG' :
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfsh.write( str(angle)+ "\t"+ str(dist_xp)+"\t"+str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")

											number_cxs=number_cxs+1 
										##  H2N residue HN - ASN , GLN, LYS, ARG, HIS, TRP
										if neighbors_p.name == 'ND2' or neighbors_p.name == 'NE2' or neighbors_p.name == 'NZ' or  neighbors_p.name == 'NE' or neighbors_p.name == 'NH1' or neighbors_p.name == 'NH2'or neighbors_p.name == 'ND1' or neighbors_p.name == 'NE1' :
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfnh.write( str(angle)+ "\t"+ str(dist_xp)+"\t"+str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")

											number_cxns=number_cxns+1 

										##  HN backbone related - N
										if neighbors_p.name == 'N' :
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfn_bb.write( str(angle)+ "\t"+ str(dist_xp)+"\t"+str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")
											number_cxnbb=number_cxnbb+1 

										##  C=O backbone related 
										if neighbors_p.name == 'O' :
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfcarbonyl.write( str(angle)+ "\t"+ str(dist_xp)+"\t"+str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")
											number_cxco=number_cxco+1 


										##  C=O on side chains related, ASN, GLN, ASP, GLU
										if neighbors_p.name == 'OD1' or neighbors_p.name == 'OE1' or neighbors_p.name == 'OD2' or neighbors_p.name == 'OE2':
											angle = calangle(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dist_xp = caldist(neighbors_c.get_coord(), neighbors_x.get_coord(), neighbors_p.get_coord())
											dfcarbonyl2_sidechain.write( str(angle)+ "\t"+ str(dist_xp)+"\t"+str(neighbors_c.get_name())+"\t"+ str(neighbors_x.get_name())+ "\t"+ str(neighbors_p.get_name())+"\t"+str(neighbors_p.get_parent())+"\t\t"+str(pdb1)+"\n")
											number_cxsideco=number_cxsideco+1 
		os.system('rm %(pdb1)s' % locals())
	except:
		number_fail=number_fail+1 
		dffail_info.write(str(name)+" ")
		pass
	
dfcount_info.write("all: "+str(number_all)+"\n"+"failed: "+str(number_fail)+"\n")
dfcount_info.close()
dffail_info.close()
dfoh.close()
dfsh.close()
dfnh.close()
dfn_bb.close()
dfcarbonyl.close()
dfcarbonyl2_sidechain.close()
dfcx_info.close()
