import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("paper")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

parser.add_argument('--part', type=str, default='yes',
                    help='select "full" for FULL COMPLEX and "yes" for SMALL COMPLEX. Default: "yes"')

args = parser.parse_args()

ifile =  args.input
part = args.part

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['resid','resname','chain','bfactor','rmsf'])
# print(data)
if part=='full':
    chain_list=['A','B','E','F','I']
    yrange_list =[[0.5,6.5],[0,5],[0,6],[0,6],[0,6]]
    chain_name_list=['ITGAV','ITGB8','LTGFb1 Dimer 1','LTGFb1 Dimer 2','GARP']
    fig,axes = plt.subplots(5,1)
else:
    chain_list=['A','B','I']
    yrange_list =[[0.5,6.5],[0,5],[0,6]]
    chain_name_list=['LTGFb1 Dimer 1','LTGFb1 Dimer 2','GARP']
    fig,axes = plt.subplots(3,1)

for ind,ax in enumerate(axes):
    # print(chain_list[ind])
    ax1 = sns.lineplot(data=data.query("chain=='%s'"%(chain_list[ind])),
                    x='resid',
                    y='rmsf',
                    color='tab:blue',
                    ax=ax
                    # hue='chain'
                    )
    ax1.set_title("%s"%(chain_name_list[ind]))
    ax1.set_ylabel("RMSF",color='tab:blue')
    ax1.set_ylim(yrange_list[ind][0],yrange_list[ind][1])
    
    ax2 = ax1.twinx()
    # ax2.set_ylim([0,130])
    ax2.set_ylabel("B-Factor",color='tab:red')
    data.query("chain=='%s'"%(chain_list[ind])).plot(x="resid", y="bfactor", ax=ax2, legend=False, color="tab:red")
# ax1.figure.legend()

# ### MISCELLANEOUS ###
plt.suptitle("%s"%(ifile[:-4]),va='top')
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
plt.gcf().set_size_inches(7.5,6.5)   ## Wide x Height
# plt.locator_params(axis='both', nbins=5)
plt.tight_layout()
plt.savefig("%s"%(ifile[:-3]),dpi=700)
plt.show()