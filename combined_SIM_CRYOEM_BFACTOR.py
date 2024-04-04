import MDAnalysis as mda
import pandas as pd
import argparse



parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--sim_data', type=str, default='',
                    help='INPUT to RMSF profile extracted from SIMULATIONS')
parser.add_argument('--cryo_pdb', type=str, default='',
                    help='INPUT to RMSF profile extracted from CRYOEM')

args = parser.parse_args()

ifile1 =  args.sim_data
ifile2 =  args.cryo_pdb


sim_rmsf = pd.read_csv(ifile1,
                names=['chain','resname','resid','rmsf','system'],
                delim_whitespace=True)

#v = sim.loc[(sim['resid']=='1') & (sim['resname']=='LEU')]
sim_rmsf.loc[(sim_rmsf['chain']=='LTGFb1A'),'chain']='A'
sim_rmsf.loc[(sim_rmsf['chain']=='LTGFb1B'),'chain']='B'
sim_rmsf.loc[(sim_rmsf['chain']=='GARP'),'chain']='I'

# print(sim.loc[(sim['chain']=='GARP')])
sim_rmsf_copy = sim_rmsf.copy()
sim_rmsf_copy = sim_rmsf_copy.query("system=='REP1'")

cryo_universe = mda.Universe(ifile2,ifile2)
ca_cryo = cryo_universe.select_atoms("name CA")

## CHANGE ALL CHAIN GARP to CHAIN I
for id,name,chain,bfactor in zip(ca_cryo.resids,ca_cryo.resnames,ca_cryo.segids,ca_cryo.bfactors):
    if chain!='I':
        rmsf_per_residue = sim_rmsf_copy.loc[(sim_rmsf_copy['resid']=='%s'%(id)) & (sim_rmsf_copy['resname']=='%s'%(name)) & (sim_rmsf_copy['chain']=='%s'%(chain))]
    else:
        rmsf_per_residue = sim_rmsf_copy.loc[(sim_rmsf_copy['resid']=='%s'%(id-19)) & (sim_rmsf_copy['resname']=='%s'%(name)) & (sim_rmsf_copy['chain']=='%s'%(chain))]
        # print(rmsf_per_residue)
    if rmsf_per_residue is not None:
        print(id,name,chain,bfactor,*rmsf_per_residue['rmsf'])