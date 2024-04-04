#!/bin/bash

name='RMSF_visualization_REP2.pdb'
grep PROA rmsf_tempfactors_chain_LTGFb1A.pdb > rmsf_tempfactors_chain_LTGFb1A.pdb2
grep PROB rmsf_tempfactors_chain_LTGFb1B.pdb > rmsf_tempfactors_chain_LTGFb1B.pdb2
grep PROI rmsf_tempfactors_chain_GARP.pdb > rmsf_tempfactors_chain_GARP.pdb2

cat rmsf_tempfactors_chain_LTGFb1A.pdb2 rmsf_tempfactors_chain_LTGFb1B.pdb2 rmsf_tempfactors_chain_GARP.pdb2 > $name

