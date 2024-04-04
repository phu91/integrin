# Updated on 11/29/2022
import MDAnalysis as mda
import pandas as pd
import numpy as np
import math, sys
import argparse
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import rms, align
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# FUNCTIONS
def rmsf_calculation(selected_str,skipping):
    # ref = u.select_atoms(selected_str)      ## FIRST FRAME
    # Aligning the traj to the REF frame (CRYOEM frame). ALIGNMENT ON EACH SEGMENT
    # print(len(u.atoms))
    average = align.AverageStructure(u, u, select=selected_str,
                                 ref_frame=0).run(step=skipping)  ## THE REFERENCE is the the first frame = cryoEM frame. 
    print("Aligning TRAJ: %s  || to the AVERAGED structure || %s"%(selected_str,selected_str))
    ref = average.results.universe
    aligner = align.AlignTraj(u, ref,
                          select=selected_str,
                          in_memory=True).run(step=skipping)

    selected_segment = u.select_atoms(selected_str)
    R = rms.RMSF(selected_segment, verbose=True).run(step=skipping)
    return R

def extract_frames(startFrame,endFrame):
    all = u.select_atoms("all")
    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.pdb"%(startFrame,endFrame)) as pdb:
        pdb.write(all)

    with mda.Writer("SEGMENT_FRAME_%s_TO_%s.xtc"%(startFrame,endFrame), all.n_atoms) as W:
        for ts in u.trajectory[startFrame:endFrame:traj_skip]:
            W.write(all)


# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Required positional argument
parser.add_argument('--top', type=str,
                    help='A required topology (PSF, PDB, etc.) file')

parser.add_argument('--traj', type=str,
                    help='A required topology (XTC, DCD, etc.) file')

parser.add_argument('--begin', type=int, default=1,
                    help='Starting Frame. Default = 1 FRAME 1')

parser.add_argument('--end', type=int, default=-1,
                    help='Ending Frame. Default = -1 ALL FRAME')

parser.add_argument('--skip', type=int, default=1,
                    help='Skipping rate for Radial Distribution Function calculations. Default = 1 frames')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Add a system name to output file')


args = parser.parse_args()


top_file =  args.top
traj_file = args.traj
traj_skip = args.skip
traj_begin = args.begin
traj_end = args.end
systemname = args.system

u = mda.Universe(top_file,traj_file,in_memory=True)
u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
n_atom_origin = len(u.atoms)

## selected strings
chainA = 'protein and segid PROA and name CA'
chainB = 'protein and segid PROB and name CA'
chainE = 'protein and segid PROE and name CA'
chainF = 'protein and segid PROF and name CA'
chainI = 'protein and segid PROI and name CA'

## REFERENCE FRAME 0
chainA_frame_0 = u.select_atoms('protein and segid PROA and name CA')
chainB_frame_0 = u.select_atoms('protein and segid PROB and name CA')
chainE_frame_0 = u.select_atoms('protein and segid PROE and name CA')
chainF_frame_0 = u.select_atoms('protein and segid PROF and name CA')
chainI_frame_0 = u.select_atoms('protein and segid PROI and name CA')
# print(chainA.residues)
ref_list = [chainA_frame_0,chainB_frame_0,chainE_frame_0,chainF_frame_0,chainI_frame_0]
chain_str = [chainA,chainB,chainE,chainF,chainI]
chain_list = ['ITGAV','ITGB8','LTGFb1A','LTGFb1B','GARP']
chain_name_list = ['A','B','E','F','I']

if traj_end != -1:
    extract_frames(traj_begin,traj_end)
    u = mda.Universe("SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end),
                     "SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end),
                     in_memory=True)
    print("\n########################################################")
    print("TOP : SEGMENT_FRAME_%s_TO_%s.pdb"%(traj_begin,traj_end))
    print("TRAJ: SEGMENT_FRAME_%s_TO_%s.xtc"%(traj_begin,traj_end))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("########################################################\n")
else:
    print("\n########################################################")
    print("TOP : %s"%(top_file))
    print("TRAJ: %s"%(traj_file))
    print("NUMBER OF ATOMS (ORIGINAL): %s"%(n_atom_origin))
    print("NUMBER OF ATOMS (CURRENT) : %s"%(len(u.atoms)))
    print("########################################################\n")
    pass

with open("RMSF_%s.dat"%(systemname),"w+") as rmsf_out:
    for ind, (chain) in enumerate(chain_str):
        # print(ind,chain)
        rmsf_out.write("#chain resname resid rmsf sysname\n")
        RMSF = rmsf_calculation(chain,traj_skip)
        
        # print()
        # print(len(RMSF.results.rmsf))
        # print(len(ref_list[ind].residues))
        for RES,rmsf in tqdm(zip(ref_list[ind].residues,RMSF.results.rmsf),total=len(ref_list[ind].residues),desc=chain):
            print(RES.resid,chain_name_list[ind],rmsf)
            # if ind ==0:
                # print("A",RES.resid)
                # rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
                # rmsf_out.flush()
                # RES.atoms.tempfactors = rmsf
                # RES.chain==chain_name_list[ind]
                # ref.atoms.write('rmsf_tempfactors_%s_%s.pdb'%(chain_list[ind],systemname))
#             if ind ==1:
#                 # print("B",RES.resid)
#                 rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
#                 rmsf_out.flush()
#                 RES.atoms.tempfactors = rmsf
#                 ref.atoms.write('rmsf_tempfactors_%s_%s.pdb'%(chain_list[ind],systemname))
#             if ind ==2:
#                 # print("I",RES.resid)
#                 rmsf_out.write("%s\t%s\t%s\t%s\t%s\n"%(chain_list[ind],RES.resname,RES.resid,rmsf,systemname))
#                 rmsf_out.flush()
#                 RES.atoms.tempfactors = rmsf
#                 ref.atoms.write('rmsf_tempfactors_%s_%s.pdb'%(chain_list[ind],systemname))