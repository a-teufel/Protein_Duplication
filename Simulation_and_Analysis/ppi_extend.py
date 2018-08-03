#!/usr/bin/python
# program to calculate a bunch of metrics from pdbs from the forward simulation, also calculates the ability of these evolved proteins to bind ancestral partners

import re, shutil, random, os, math, sys, subprocess, glob
import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

#from toolbox import mutate_residue
from toolbox import cleanATOM
from rosetta.protocols.analysis import *
from rosetta.protocols.vip import *
from rosetta.protocols.grafting import *
from time import time

import numpy as np
from Bio import *
from Bio.PDB import *
from interface_analysis import *

def main():
  args  =  sys.argv
  in_file = args[1]
  out_file = args[2]
  distance_cutoff = float(sys.argv[3])
  init(extra_options='-mute basic -mute core -mute protocols -mute Warning')

  all_lines = (open(in_file, 'r')).readlines()
  print(len(all_lines))
  
  #get the protein used to initalize the forward simulation
  initial_pose = pose_from_pdb(str('burn1ABC_renumb.pdb'))
  
  #save each one of its changes 
  chains=initial_pose.split_by_chain()

  ancestral1 = chains[1]
  ancestral2 = chains[2]
  ancestral3 = chains[3]

  ancestral1.dump_pdb("Ans_A.pdb")
  ancestral2.dump_pdb("Ans_B.pdb")
  ancestral3.dump_pdb("Ans_C.pdb")  


  ancestral_structure1=capture_pdb_one("Ans_A_cap.pdb","Ans_A.pdb")
  ancestral_structure2=capture_pdb_one("Ans_B_cap.pdb","Ans_B.pdb")
  ancestral_structure3=capture_pdb_one("Ans_C_cap.pdb","Ans_C.pdb")

  all_data = []
  i=0
  for a_line in all_lines:
    split = a_line.split(',')
    if split[0] == 'Variant':
      continue
    if split[0] == 'WT':
      continue
    else:
      print(split[0])
      pos=re.sub("[^0-9^.]", "", split[0])
	  #figure out if a position is in A B or C
      print(pos)
      if int(pos) <= ancestral1.total_residue():
        all_data.append([i,pos,split[0],split[1],split[2],split[3],split[4],split[5],split[6],split[7],split[8], 'A'])
		i=i+1
      if int(pos) > ancestral1.total_residue() and int(pos) <= ancestral1.total_residue()+ancestral2.total_residue():
		all_data.append([i,pos,split[0],split[1],split[2],split[3],split[4],split[5],split[6],split[7],split[8], 'B'])
		i=i+1
      if int(pos) > ancestral1.total_residue()+ancestral2.total_residue():
		all_data.append([i,pos,split[0],split[1],split[2],split[3],split[4],split[5],split[6],split[7],split[8], 'C'])
		i=i+1

  #make the output datafile
  output = open(out_file , 'w')
  to_file = 'name\tcount\tidentityA\tidentityB\tindentityC\tindentityBC\tinterface_identity_AB\tinterface_identity_BA\tinterface_identity_AC\tinterface_identity_CA\tinterface_identity_BC\tnoninterface_identity_AB\tnoninterface_identity_BC\tnoninterface_identity_AC\tnoninterface_identity_CA\tnoninterface_identity_BC\tevolved_interaction_AB\tevolved_interaction_AB_recalc\tevolved_interaction_AC\tevolved_interaction_AC_recal\tancestral_interaction_AB\tancestral_interaction_BA\tancestral_interaction_AC\tancestral_interaction_CA\tstabilityA\tstabilityB\tstabilityC\n'
  output.write(to_file)
  output.close()

  #Make binding comparison using ancestral chain A
  for i in range(0, len(all_data)):
    interval = 10
    #score_ob = pyros.Scores()

    a_mutant = all_data[i]

    print(a_mutant)
    if i % interval == 0:

        print(a_mutant[0])
		#split each chain up
        current_struct = pose_from_pdb(str(a_mutant[0])+'.pdb')
		chains=current_struct.split_by_chain()
		A=chains[1]
		B=chains[2]
		C=chains[3]
        A.dump_pdb("tempA.pdb")
		B.dump_pdb("tempB.pdb")
		C.dump_pdb("tempC.pdb")
	
		new_a=capture_pdb_one("tempA_cap.pdb","tempA.pdb")
		new_b=capture_pdb_one("tempB_cap.pdb","tempB.pdb")
		new_c=capture_pdb_one("tempC_cap.pdb","tempC.pdb")
		newAB = "A_B_combined_current.pdb"
		newcombined_fileAB = make_combined_file([new_a,new_b], newAB)

		#calc evolved interactions
		newAB_pose = pose_from_pdb(newAB)
		inter_Evo_AB = InterfaceEnergy_split(newAB_pose)

		newAC = "A_C_combined_current.pdb"
		newcombined_fileAC = make_combined_file([new_a,new_c], newAC)

		newAC_pose = pose_from_pdb(newAC)
		inter_Evo_AC = InterfaceEnergy_split(newAC_pose)

	
		#make combinded file for A and ancest B
        new_combo1 = new_a[0:-4] + '_B_combined.pdb'
        combined_fileAB = make_combined_file([new_a,ancestral_structure2], new_combo1)
        new_combo1_rev = 'A_'+new_b[0:-4]+'_combinded.pdb'
		combined_fileBA =make_combined_file([ancestral_structure1,new_b],new_combo1_rev)

		#calc ancestral binding
		combindedAB = pose_from_pdb(new_combo1)
		combindedBA = pose_from_pdb(new_combo1_rev)
 
        inter_AB = InterfaceEnergy_split(combindedAB)
		inter_BA = InterfaceEnergy_split(combindedBA)
		
       
		#make combinded file for A and ancest C
        new_combo2 = new_a[0:-4] + '_C_combined.pdb'
        combined_fileAC = make_combined_file([new_a,ancestral_structure3], new_combo2)
		new_combo2_rev = 'A_'+new_c[0:-4]+'_combinded.pdb'
		combined_fileCA = make_combined_file([ancestral_structure1, new_c], new_combo2_rev)	
	

		combindedAC = pose_from_pdb(new_combo2)
		combindedCA = pose_from_pdb(new_combo2_rev)

        inter_AC = InterfaceEnergy_split(combindedAC)
		inter_CA = InterfaceEnergy_split(combindedCA)
       

		#calc id from new A and ancest A
        identityA = calc_identity(ancestral_structure1, new_a)

		#calc id from new B and ancest B
        identityB = calc_identity(ancestral_structure2, new_b)
        
		#calc id from new C and ancest C
        identityC = calc_identity(ancestral_structure3, new_c)
       
		#calc B C diver
		identityBC = calc_identity(new_b, new_c)
		print(identityBC)


		#calc the distance between current A and B
        p = PDBParser()
        structureAB = p.get_structure('temp_AB.pdb', newcombined_fileAB)
        distances = calculate_ca_distance(structureAB, 'A')
  
		#get inferface site between A and B
        interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)
	
        interface_identityAB = calc_identity_sites(ancestral_structure1, new_a, interface_sites)

		#get non inferface site between A and B
        non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
        non_interface_identityAB = calc_identity_sites(ancestral_structure1, new_a, non_interface_sites)


		#calc the distance between current B and A
        p = PDBParser()
        structureBA = p.get_structure('temp_BA.pdb', newcombined_fileAB)
        distances = calculate_ca_distance(structureBA, 'B')

		#get inferface site between B and A
        interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)
        interface_identityBA = calc_identity_sites(ancestral_structure2, new_b, interface_sites)
		interface_identityBC = calc_identity_sites(new_b, new_c, interface_sites)

		#get non inferface site between B and A
        non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
        
        non_interface_identityBA = calc_identity_sites(ancestral_structure2, new_b, non_interface_sites)
		non_interface_identityBC = calc_identity_sites(new_b, new_c, non_interface_sites)



		#calc the distance between current A and C (this should always be the same a dist between A and B) this just does the same stuff as above
		p = PDBParser()
        structureAC = p.get_structure('temp_AC.pdb', newcombined_fileAC)
        distances = calculate_ca_distance(structureAC, 'A')
       
        interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)
        interface_identityAC = calc_identity_sites(ancestral_structure1, new_a, interface_sites)
 
        non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
        non_interface_identityAC = calc_identity_sites(ancestral_structure1, new_a, non_interface_sites)

	
		p = PDBParser()
        structureCA = p.get_structure('temp_CA.pdb', newcombined_fileAC)
        distances = calculate_ca_distance(structureCA, 'C')
      
        interface_sites = (distances[distances[:, 1] <= distance_cutoff][:, 0]).astype(int)

        interface_identityCA = calc_identity_sites(ancestral_structure3, new_c, interface_sites)
        
        non_interface_sites = (distances[distances[:, 1] > distance_cutoff][:, 0]).astype(int)
     
        non_interface_identityCA = calc_identity_sites(ancestral_structure3, new_c, non_interface_sites)

        #get stabilities of each chain
		stabA=str(a_mutant[3]) 
		stabB=str(a_mutant[4])
		stabC=str(a_mutant[5])  
		
		output = open(out_file , 'a')
    	to_file = str(a_mutant[1]) + '_' + str(a_mutant[11]) +'\t' +str(i) +'\t' + str(identityA) + '\t' + str(identityB)+ '\t'  + str(identityC) + '\t'  + str(identityBC)+'\t' +str(interface_identityAB)+'\t' +str(interface_identityBA) + '\t'  +str(interface_identityAC) + '\t'  +str(interface_identityCA)+ '\t'+ str(interface_identityBC) +'\t' +str(non_interface_identityAB) +'\t' +str(non_interface_identityBC) +'\t' +str(non_interface_identityAC) +'\t' +str(non_interface_identityCA) +'\t' +str(non_interface_identityBC)+'\t' + str(a_mutant[6]) +'\t'+ str(inter_Evo_AB) +'\t' + str(a_mutant[7]) + '\t'+ str(inter_Evo_AC) +'\t'+ str(inter_AB) + '\t'+ str(inter_BA)  + '\t'+ str(inter_AC) + '\t'+ str(inter_CA)  + '\t' + str(stabA) + '\t' + str(stabB) + '\t' + str(stabC) +'\n'
    	output.write(to_file)
    	output.close()

  print("done")

def capture_pdb(out_name, fn, chain_letter):
  parser = PDB.PDBParser()
  structure = parser.get_structure('working_pdb', fn)
  
  writer = PDB.PDBIO()
  writer.set_structure(structure)
  writer.save(out_name, select=SelectChains(chain_letter))
  return(out_name)


def capture_pdb_one(out_name, fn):
  parser = PDBParser()
  structure = parser.get_structure('working_pdb', fn)
  
  writer = PDBIO()
  writer.set_structure(structure)
  writer.save(out_name)
  return(out_name)
    
def calc_identity(fn1, fn2):
  parser = PDBParser()
  structure1 = parser.get_structure('working_pdb1', fn1)
  structure2 = parser.get_structure('working_pdb2', fn2)

  count = 0
  identical = 0

  pp1 = ''
  ppb1 = PPBuilder()
  for pp in ppb1.build_peptides(structure1):
	pp1 += pp.get_sequence()

  pp2 = ''
  ppb2 = PPBuilder()
  for pp in ppb2.build_peptides(structure2):
	pp2 += pp.get_sequence()  

  for i in range(0, len(pp1)):
     if pp1[i] == pp2[i]:
       identical += 1
     count += 1

  return(float(identical)/float(count))

def calc_identity_sites(fn1, fn2, sites):
  parser = PDBParser()
  structure1 = parser.get_structure('working_pdb1', fn1)
  structure2 = parser.get_structure('working_pdb2', fn2)

  count = 0
  identical = 0

  pp1 = ''
  ppb1 = PPBuilder()
  for pp in ppb1.build_peptides(structure1):
	pp1 += pp.get_sequence()

  pp2 = ''
  ppb2 = PPBuilder()
  for pp in ppb2.build_peptides(structure2):
	pp2 += pp.get_sequence()  

  for i in sites:
     if pp1[int(i)] == pp2[int(i)]:
       identical += 1
     count += 1

  return(float(identical)/float(count))



def make_combined_file(files, path):
  with open(path, 'w') as outfile:
    for fname in files:
      with open(fname) as infile:
        for line in infile:
          outfile.write(line)
  subprocess.call("sed -i 's/END//g' " + path, shell=True)
  return(path)

def clean_up(delete_files):
  import glob
  for name in delete_files:
    for files in glob.glob(name):
      if os.path.isfile(files):
        os.remove(files)

def rmWS(string):
    return("".join(string.split()))

def InterfaceEnergy_split(pdb):

  sf = get_fa_scorefxn()
  interface_mover = InterfaceAnalyzerMover(1, False, sf, False, True, True, False )
  interface_mover.apply(pdb)
  return(interface_mover.get_interface_dG())


#Run main program
if __name__ == '__main__':
   main()
