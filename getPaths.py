"""
TITLE : getPaths
AUTHOR : Hernansaiz Ballesteros, Rosa.
            rosa.hernansaiz@bioquant.uni-heidelberg.de
DESCRIPTION : Example of how to use topo_fun.py applied to CARNIVAL's results
                of the Cell lines analysed in CARNIVAL1000
                (https://github.com/saezlab/CARNIVAL1000).
             Get pickle with the dictionary structure:
             {sample:[icar, nodI, nodE, isoNode, subpaths]}
             Get netsum.csv with some statistics of the pickles
             (sample, cancer_tissue, network, #vertices, #edges,
              #iniciators, #effectors, #isolated, #subpaths)
LICENSE : GPL-v3
"""

#==========#
# Packages #
#==========#

import os
import time
from tqdm import tqdm
import pandas as pd
import igraph as ig
import pickle
import topo_fun as tp

#===================================#
# Parameters and required variables #
#===================================#

#set working directory
path = "/Users/rosherbal/Projects/network_tools"
os.chdir(path)

net_type =  ['omnipathNOsign']#["omnipath", "GWpkn"]
save_folder = 'CMP_noSign_subpaths/' #CMP_subpaths

carnival_folder = "../00-Networks/footprints_tissue_specific/CARNIVAL"
smpl_res = os.listdir(carnival_folder + "/" + net_type[0])
#smpl_res.remove('SIDM00319')

meta_samples = pd.read_csv('../00-CMP/metadata/model_list_2019-04-04_0916.csv',
                            sep = ",")
meta_samples = meta_samples[meta_samples['model_id'].isin(smpl_res)]

cancer_type = list(set(meta_samples['cancer_type']))

#=======
# Main #
#=======

# Create the structure
# {sample:[icar, nodI, nodE, isoNode, subpaths]}
# also a pd.df with the metrics of all the previous objects
for i in tqdm(range(len(cancer_type))):
    ctype = cancer_type[i]
    print(ctype)
    samples = meta_samples[meta_samples['cancer_type'] == ctype]["model_id"].tolist()
    for network in net_type:
        lstats = []
        netWord = {}
        for sample in samples:

            # read information per sample
            pcar = carnival_folder + "/" + network + "/" + sample

            # add information to a dictionary
            netWord[sample] = tp.getTopology(Cfile = pcar)

            # get some metrics
            aux = {"sample": sample, "cancer_tissue": ctype,
                     "network": network,
                     "vertices": netWord[sample][0].vcount(),
                     "edges": netWord[sample][0].ecount(),
                     "iniciators": len(netWord[sample][1]),
                     "effectors": len(netWord[sample][2]),
                     "isolated": len(netWord[sample][3]),
                     "subpaths": len(netWord[sample][4])}
            lstats.append(aux)

        with open(save_folder + ctype + '_netextend-' + net_type[0] + '.pkl', 'wb') \
                    as ofile:
            pickle.dump(netWord, ofile)
        ofile.close()

        stats = pd.DataFrame(lstats)
        stats.to_csv(save_folder + ctype + '_netsum-' + net_type[0] + '.csv',
                        header=True, index=None)
