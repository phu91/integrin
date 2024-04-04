import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import argparse
from matplotlib import rc
from statannot import add_stat_annotation

sns.set_context("notebook")

parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--input', type=str, default='',
                    help='INPUT to RMSF profile')

args = parser.parse_args()

ifile =  args.input

data = pd.read_csv(ifile,comment='#',
                   delim_whitespace=True,
                   names=['CHAIN','RESNAME','RESID','RMSF','SYSTEM'])
# print(data['RMSF'].max())
# u=data.query("SYSTEM=='first100ns'").reset_index()
# v=data.query("SYSTEM=='last100ns'").reset_index()

# # print(u,v)
# delta=(v['RMSF']-u['RMSF'])

# delta_df = pd.DataFrame({'RESID':u['RESID'],
#                       'CHAIN':u['CHAIN'],
#                       'DELTA_RMSF':delta})

# # print(delta_df)

# order=['LTGFb1A','LTGFb1B','GARP']
# box_pairs = [
#     (("LTGFb1A", "first100ns"), ("LTGFb1A", "last100ns")),
#     (("LTGFb1B", "first100ns"), ("LTGFb1B", "last100ns")),
#     (("GARP", "first100ns"), ("GARP", "last100ns")),
# ]

# ##
# # hue='SYSTEM'

# system_list=['first100ns','last100ns']

# g = sns.boxplot(data=data,
#              x='CHAIN',
#              y='RMSF',
#              style='SYSTEM',
#              hue='SYSTEM'
#             #  split=True,
#             #  showfliers=False,
#             #  lw=3,
#              )

# # g.ylim([0,6])
# # fig,axes = plt.subplots(2,1,sharex=False,sharey=True,layout="constrained")
# # axes = axes.flatten()

# # for ind,ax in enumerate(axes):
# #     # print(ax)
# #     sns.boxplot(data=data.query("SYSTEM=='%s' and CHAIN !='GARP'"%(system_list[ind])),
# #                  x='CHAIN',
# #                  y='RMSF',
# #                 #  hue='CHAIN',
# #                 #  lw=3,
# #                  ax=ax)

order=['LTGFb1A','LTGFb1B','GARP']

fig,axes = plt.subplots(3,1,sharex=True,sharey=True,layout="constrained")
axes = axes.flatten()
for ind,ax in enumerate(axes):
    g = sns.lineplot(data=data.query("CHAIN=='%s'"%(order[ind])),
                 x='RESID',
                 y='RMSF',
                #  hue='CHAIN',
                #  size='RMSF',
                 hue='SYSTEM',
                 lw=3,
                 ax=ax
                 )
    g.set_title("CHAIN %s"%(order[ind]))
# test_results = add_stat_annotation(g,data=data, x='CHAIN', y='RMSF',hue='SYSTEM',
#                                    box_pairs=box_pairs,
#                                    test='Mann-Whitney', text_format='star',
#                                    loc='inside', verbose=2)

# # ### MISCELLANEOUS ###
# plt.suptitle("%s"%(ifile[:-4]),va='top')
# plt.rcParams['ps.useafm'] = True
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rcParams['pdf.fonttype'] = 42
# plt.gcf().set_size_inches(7.5,5)   ## Wide x Height
# # plt.locator_params(axis='both', nbins=5)
# plt.tight_layout()
# plt.savefig("%s"%(ifile[:-3]))
plt.show()