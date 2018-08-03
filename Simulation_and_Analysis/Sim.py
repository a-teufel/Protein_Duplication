#!/share/apps/python-2.7.2/bin/python
# Program to simulate the evolution of a protein when one chain of a hetrodimer protein is duplicated
# For a more basic version of this code please see
# Teufel AI, Wilke CO. Accelerated simulation of evolutionary trajectories in origin-fixation models. Journal of The Royal Society Interface. 2017 Feb 1;14(127):20160906.
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# run this program like: python Sim.py mypdb
# where mypdb is the name of the pdb file you want to use without the ".pdb" part
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 

# Imports
import sys
import argparse
import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover


from toolbox import cleanATOM
from rosetta.protocols.analysis import *
from rosetta.protocols.vip import *
from rosetta.protocols.grafting import *
from time import time


def main():
    #takes name of pdb file without the extention
    args =  sys.argv	
    pdb_file = args[1]
    out_file = args[2]
    score_type = int(args[3])
    #set up timer to figure out how long the code took to run
    t0=time()

    # Initialize Rosetta.
    init(extra_options='-mute basic -mute core -mute protocol -mute warn')

    # Constants
    PACK_RADIUS = 5
    #Amino acids, notice there is no C
    AAs = ("A","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    #Number of mutations to accept
    max_accept_mut = 2000
    #Population size
    N = 1
    #Beta (temp term)
    beta = 1

    #Prepare data headers
    data = ['Variant,ChainA,ChainB,ChainC,InterfaceAB,InterfaceAC,"delta-delta-G",Probability,Generation\n']

    initial_pose = pose_from_pdb(pdb_file)

    #Set up ScoreFunction
    sf = get_fa_scorefxn()
       
    #Set up MoveMap This is where you turn the bb and side chain flexibility on and off
    mm = MoveMap()
    mm.set_bb(False)

    #Get the init score of the struct to calc the threshold
    pre_pre_packing_score = sf(initial_pose)
    print(pre_pre_packing_score)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('dfpmin_armijo_nonmonotone')

    cp_init_pdb = Pose()
    cp_init_pdb.assign(initial_pose)
    chains=cp_init_pdb.split_by_chain()

    #split up AB inter and AC inter 
    initial_poseAB = Pose()
    initial_poseAB.assign(initial_pose)
    initial_poseAC = Pose()
    initial_poseAC.assign(initial_pose)

    init_chain_moverAB = SwitchChainOrderMover()
    init_chain_moverAB.chain_order("12")
    init_chain_moverAB.apply(initial_poseAB)

    init_chain_moverAC = SwitchChainOrderMover()
    init_chain_moverAC.chain_order("13")
    init_chain_moverAC.apply(initial_poseAC)

    #score the inital stabs of each chain
    wt_a=sf(chains[1])

    wt_b=sf(chains[2])

    wt_c=sf(chains[3])

    #score the intial interfaces 
    inter_AB=InterfaceEnergy_split(initial_poseAB)

    inter_AC=InterfaceEnergy_split(initial_poseAC)

    #init thresholds set to half of the init stabilities, if you want to do a different protein change these
    threshold_a=-138.41754752
    threshold_b=-61.378619136
    threshold_c=-61.378619136
    threshold_inter_ab=-10.3726691079
    threshold_inter_ac=-10.3726691079

    data.append('WT,' + str(wt_a)+','+str(wt_b)+','+str(wt_c)+','+str(inter_AB)+','+str(inter_AC)+',0.0,0.0,0\n')

	#check the inital starting score
    init_score=score_all(initial_pose,sf,min_mover,beta,threshold_a, threshold_b, threshold_c,threshold_inter_ab,threshold_inter_ac,score_type)
    print(init_score)

    #number of residues to select from
    n_res = initial_pose.total_residue()
    print(n_res)
  
    #start sim
    i=0
    gen=0
    while i < max_accept_mut:
            #update the number of generations that have pased
            gen+=1

	    print 'accepts:', i 

	    #pick a place to mutate
	    mut_location = random.randint(1, n_res)
	    #mut_location = random.randint(1, 10)

	    #get the amino acid at that position
	    res = initial_pose.residue(mut_location)

	    #don't mess with C, just choose again
	    while(res.name1() == 'C'):
			mut_location = random.randint(1, n_res)
	    	#get the amino acid at that position
	    	res = initial_pose.residue(mut_location)


	    #choose the amino acid to mutate to
	    toname = res.name1()
	    new_mut_key = random.randint(0,len(AAs)-1)
	    proposed_res = AAs[new_mut_key]
	  
	    #don't bother mutating to the same amino acid it just takes more time
	    while(proposed_res == res.name1()):
			new_mut_key = random.randint(0,len(AAs)-1)
	        proposed_res = AAs[new_mut_key]

	    #init mutant with current 
	    mutant_pose = Pose()
	    mutant_pose.assign(initial_pose)
		
		#mutate 
	    mutant_pose=mutate_residue_chain(mutant_pose, mut_location, proposed_res, PACK_RADIUS, sf)
		
	    #score mutant
	     mut_score=score_all(mutant_pose,sf,min_mover,beta,threshold_a, threshold_b, threshold_c,threshold_inter_ab,threshold_inter_ac,score_type)

	    #get the probability that the mutation will be accepted
	    probability = calc_prob_scores(mut_score['score'], init_score['score'], N)
		
	    rand = random.random()

	    #test to see if mutation is accepted
	    if float(rand) < float(probability):
			print "accepted" 	
		
			#make a name for the new mutant
			variant_name = str(toname) + str(initial_pose.pdb_info().number(mut_location)) + str(proposed_res)


			# Assuming some burn in phase, make this zero if you want to store everything
			if i>=0:
				#save name and energy change
				data.append(variant_name +',' + str(mut_score['a'])+','+str(mut_score['b'])+','+str(mut_score['c'])+','+str(mut_score['ab'])+','+str(mut_score['ac'])+',' + str(mut_score['score'] - init_score['score']) + "," + str(probability) + "," + str(gen) + "\n")

				#save the new accepted mutation	
				pdb_name=str(i)+".pdb"	
				mutant_pose.dump_pdb(pdb_name)

			#update the wildtype 
			initial_pose = mutant_pose
			init_score = mut_score

			#update number of accepts
	    	i+=1

    #end of sim
    print '\nMutations and scoring complete.'
    t1 = time()

    # Output results
    data_filename = out_file
    with open(data_filename, "w") as f:
        f.writelines(data)

    print 'Data written to:', data_filename
    print 'program takes %f' %(t1-t0)


###assorted functions that have to do with scoring and prob of acceptance ####

#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):
 
  xi = calc_x(stab_org, beta, thresholds)
  xj = calc_x(stab_mut, beta, thresholds)

  if xj >= xi:
    return(float(1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))


#score functions for met-hastings selection
def calc_prob_scores(stab_mut, stab_org,N):
 
  xi = stab_org
  xj = stab_mut

  if xj >= xi:
    return(float(1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))


def calc_x(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total += -math.log(safe_calc(exponent) + 1)
  return(total)


def calc_x_list(data, beta, thresholds):
  total = 0
  print(data)
  for i in range(0, len(data)):
    #Need to make sure you check numbers that are too big for the math library
    exponent = float(beta) * (float(data[i]) - float(thresholds[i]))
    total += -math.log(safe_calc(exponent) + 1)
  return(total)

def calc_x_list_2(data, beta, thresholds):
  total = 0
  print(data)
  for i in range(0, len(data)-1):
    print(data[i])
    #Need to make sure you check numbers that are too big for the math library
    exponent = float(beta) * (float(data[i]) - float(thresholds[i]))
    total += -math.log(safe_calc(exponent) + 1)
  #subtract off the delterious binding
  exponent = float(beta) * (float(data[len(data)-1]) - float(thresholds[len(data)-1]))
  total -= -math.log(safe_calc(exponent) + 1)
  return(total)


def calc_x_list_3(data, beta, thresholds):
  total = 0
  print(data)
  for i in range(0, len(data)-2):
    #Need to make sure you check numbers that are too big for the math library
    exponent = float(beta) * (float(data[i]) - float(thresholds[i]))
    total += -math.log(safe_calc(exponent) + 1)

  exponent = float(beta) * (float(data[len(data)-2]) - float(thresholds[len(data)-2]))
  t1= -math.log(safe_calc(exponent) + 1)

  exponent = float(beta) * (float(data[len(data)-1]) - float(thresholds[len(data)-1]))
  t2= -math.log(safe_calc(exponent) + 1)
  #just add the best one to the total
  total+=max(t1,t2)
  return(total)

def safe_calc(exponent):
  if exponent > 700:
    print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))


def InterfaceEnergy_split(pdb):

  sf = get_fa_scorefxn()
  interface_mover = InterfaceAnalyzerMover(1, False, sf, False, True, True, False )
  interface_mover.apply(pdb)
  return(interface_mover.get_interface_dG())



def score_all(pdb,sf,min_mover,beta,threshold_a, threshold_b, threshold_c,threshold_inter_ab,threshold_inter_ac,score_type):

  
  cp_init_pdb = Pose()
  cp_init_pdb.assign(pdb)
  chains = cp_init_pdb.split_by_chain()

  initial_poseAB = Pose()
  initial_poseAB.assign(pdb)

  initial_poseAC = Pose()
  initial_poseAC.assign(pdb)

  init_chain_moverAB = SwitchChainOrderMover()
  init_chain_moverAB.chain_order("12")
  init_chain_moverAB.apply(initial_poseAB)

  init_chain_moverAC = SwitchChainOrderMover()
  init_chain_moverAC.chain_order("13")
  init_chain_moverAC.apply(initial_poseAC)

  init_a=sf(chains[1])

  init_b=sf(chains[2])

  init_c=sf(chains[3])

  init_inter_AB=InterfaceEnergy_split(initial_poseAB)

  init_inter_AC=InterfaceEnergy_split(initial_poseAC)


  data=[init_a,init_b,init_c,init_inter_AB,init_inter_AC]
  thresh=[threshold_a,threshold_b,threshold_c,threshold_inter_ab,threshold_inter_ac]
  

  init_score=0
  
  #score everything (bind both)  
  if score_type == 1:
	init_score=calc_x_list(data,beta,thresh)
	
  #pentalize binding B' (bind B and not B')
  if score_type == 2:
	init_score=calc_x_list_2(data,beta,thresh)

  #bind just the best duplicate (bind max)
  if score_type == 3:
	init_score=calc_x_list_3(data,beta,thresh)
   
  #just stabilitys and no binding (no bind) 
  if score_type == 4:
	data=[init_a,init_b,init_c]
    thresh=[threshold_a,threshold_b,threshold_c]
	init_score=calc_x_list(data,beta,thresh)
	
  #stabs and binding B (bind B)
  if score_type == 5:
	data=[init_a,init_b,init_c,init_inter_AB]
    thresh=[threshold_a,threshold_b,threshold_c,threshold_inter_ab]
	init_score=calc_x_list(data,beta,thresh)

  dic = {'score':init_score, 'a':init_a, 'b':init_b, 'c':init_c, 'ab':init_inter_AB, 'ac':init_inter_AC}


  return(dic)


def mutate_residue_chain(pose, mutant_position, mutant_aa,
        pack_radius = 0.0, pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose, 30, A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue(mutant_position, mutant_aa)
    #mutator.apply(test_pose)
    #test_pose = Pose()
    #test_pose.assign( pose )

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )


    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = get_fa_scorefxn() #  create_score_function('standard')

    task = standard_packer_task(pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
   
    # prevent residues from packing if they are in a different chain then the mutant
    task2=restrict_non_nbrs_from_repacking_chain(pose, mutant_position, task, pack_radius)
    task2.nonconst_residue_task(mutant_position
        ).restrict_absent_canonical_aas(aa_bool)

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover(pack_scorefxn, task2)
    packer.apply(pose)
   
    return(pose)

def restrict_non_nbrs_from_repacking_chain(pose, res, task, pack_radius):
    #print(task)
    center = pose.residue( res ).xyz( pose.residue( res ).nbr_atom() );
    #print "Res: pack radius: "+repr(pack_radius)
    for i in range(1, pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius and only repack on the same chain that is being mutated

            nbr = pose.residue( i ).xyz( pose.residue( i ).nbr_atom() );
            dist = nbr.distance(center)
        
            if dist > pack_radius:
                task.nonconst_residue_task(i).prevent_repacking()
            if dist > 0.0 and dist < pack_radius:
                task.nonconst_residue_task(i).restrict_to_repacking()
	    if not pose.pdb_info().chain(i) == pose.pdb_info().chain(res):
		task.nonconst_residue_task(i).prevent_repacking()

    return task


#Run main program
if __name__ == '__main__':
   main()

