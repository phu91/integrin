import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse,os
from matplotlib import rc
from statannot import add_stat_annotation
import MDAnalysis as mda 

# sns.set_context("poster")

def inputPart(whichPart):
    if whichPart=='full':
        chain_list=['A','B','E','F','I']
        # xrange_list = [[1,361],[1,361],[1,567]]
        yrange_list =[[0.5,6.5],[0,5],[0,6],[0,6],[0,6]]
        y2range_list =[[0,90],[0,120],[0,100],[0,100],[0,100]]
        chain_name_list=['ITGAV','ITGB8','LTGFb1 Dimer 1','LTGFb1 Dimer 2','GARP']
        fig,axes = plt.subplots(5,1)
    else:
        chain_list=['A','B','I']
        # xrange_list = [[1,361],[1,361],[1,567]]
        yrange_list =[[0,6],[0,4.5],[0,6]]
        y2range_list =[[0,130],[0,120],[0,130]]
        chain_name_list=['LTGFb1 Dimer 1','LTGFb1 Dimer 2','GARP']
        fig,axes = plt.subplots(3,1)
    return chain_list,yrange_list,y2range_list,chain_name_list,fig,axes

def generate_PDB_with_new_RMSF(selectedPDB_RMSF):
    if selectedPDB_RMSF=='yes':
        template_pdb = mda.Universe(template,template)
        select_template = template_pdb.select_atoms("protein")
        select_template.tempfactors=0
        select_template.write("template_ca_only.pdb")
        template_pdb_ca = mda.Universe("template_ca_only.pdb","template_ca_only.pdb")
        all = template_pdb_ca.select_atoms("resid 1:585")
        # print(all)
        for ind, (RESID,RESNAME,CHAIN,RMSF) in enumerate(zip(data.resid,data.resname,data.chain,data.rmsf)):
            # SEGNAME='X'
            # print(CHAIN)
            # if CHAIN=='A':
            #     SEGNAME='PROA'
            # if CHAIN=='B':
            #     SEGNAME='PROB'
            # if CHAIN=='I':
            #     SEGNAME='PROI'
            # print(RESID,RESNAME,SEGNAME)
            u = template_pdb_ca.select_atoms("resname %s and resid %s and segid %s"%(RESNAME,RESID,CHAIN))
            # print(u)
            u.tempfactors=RMSF
            # print(u.tempfactors)

        all.write("RMSF_model_%s.pdb"%(systemName))
        print("\n===> RMSF_model_%s.pdb is generated"%(systemName))
    else:
        print("\n===> No PDB with new RMSF generate!")


parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

parser.add_argument('--cryo', type=str, default='',
                    help='cryoEM PDB as template.')

parser.add_argument('--part', type=str, default='yes',
                    help='select "full" for FULL COMPLEX and "yes" for SMALL COMPLEX. Default: "yes"')

parser.add_argument('--nrep', type=int, default='1',
                    help='Number of replicas in data. Default: 1')

parser.add_argument('--system', type=str, default='UNKNOWN',
                    help='Output name')

parser.add_argument('--pdbout', type=str, default='no',
                    help='Output PDB with NEW RMSF from input')
args = parser.parse_args()

ifile =  args.input
part = args.part
template=args.cryo
systemName=args.system
nrep=args.nrep
pdbout=args.pdbout

data = pd.read_csv(ifile,comment='#',
                delim_whitespace=True,
                names=['resname','resid','chain','bfactor','rmsf','system']
                )
                #resname resid chain bfactor rmsf system
print(data)
# print(data.iloc[:1114]) REP1
# print(data.iloc[1114:]) REP2
chain_list,yrange_list,y2range_list,chain_name_list,fig,axes=inputPart(part)
generate_PDB_with_new_RMSF(pdbout)

### PLOT
weight=20

if nrep==1:
    hue=None
else:
    hue='system'

for ind,ax in enumerate(axes):
    ax1 = sns.lineplot(data=data.query("chain=='%s'"%(chain_list[ind])),
                    x='resid',
                    y='rmsf',
                    color='tab:blue',
                    hue=hue,
                    legend='auto',
                    marker='o',
                    ax=ax
                    )
    ax1.set_title("%s"%(chain_name_list[ind]),fontweight='bold')
    ax1.set_ylabel("RMSF",color='tab:blue',fontweight='bold')
    ax1.set_ylim(yrange_list[ind][0],yrange_list[ind][1])
    ax1.set_xlabel("RESID",fontweight='bold')
    ax2 = ax1.twinx()
    # ax2.set_ylim([0,130])
    ax2.set_ylim(y2range_list[ind][0],y2range_list[ind][1])
    ax2.set_ylabel("B-Factor",color='#d800a2',fontweight='bold')
    sns.lineplot(data=data.query("chain=='%s'"%(chain_list[ind])),  ### THE END OF REP 1
                    x="resid", 
                    y="bfactor",
                    color="#d800a2",
                    ax=ax2
    )

    # ax2.legend(loc='upper left')
# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42

plt.gcf().set_size_inches(7.5,8.5)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s.png"%(ifile[:-4]),dpi=700)
plt.savefig("%s.eps"%(ifile[:-4]),dpi=700)
# os.remove("template_ca_only.pdb")
plt.show()