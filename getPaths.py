"""
TITLE : getPaths
AUTHOR : Hernansaiz Ballesteros, Rosa.
            rosa.hernansaiz@bioquant.uni-heidelberg.de
DESCRIPTION : Example of how to use topo_fun.py applied to CARNIVAL's results
                of the Cell lines analysed in CARNIVAL1000
                (https://github.com/saezlab/CARNIVAL1000).
             Get pickle with the structure:
             {c_type:{network:sample:[icar, nodI, nodE, isoNode, subpaths]}}
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

net_type =  ["omnipath", "GWpkn"]

carnival_folder = "../00-Networks/footprints_tissue_specific/CARNIVAL"
smpl_res = os.listdir(carnival_folder + "/" + net_type[0])

meta_samples = pd.read_csv('../00-CMP/metadata/model_list_2019-04-04_0916.csv',
                            sep = ",")
meta_samples = meta_samples[meta_samples['model_id'].isin(smpl_res)]

cancer_type = list(set(meta_samples['cancer_type']))

cancer_tissue = {}
lstats = []

#=======
# Main #
#=======

# Create the structure
# {c_type:{network:sample:[icar, nodI, nodE, isoNode, subpaths]}}
# also a pd.df with the metrics of all the previous objects
for i in tqdm(range(len(cancer_type))):
    ctype = cancer_type[i]
    samples = meta_samples[meta_samples['cancer_type'] == ctype]["model_id"].tolist()
    netWord = {}
    print(ctype)
    for network in net_type:
        netWord[network] = {}
        print(network)
        for sample in samples:
            # read information per sample
            pcar = carnival_folder + "/" + network + "/" + sample
            smpl_pd = pd.read_csv(pcar + "/weightedModel_1.txt", sep = "\t")
            smpl_nAtt = pd.read_csv(pcar + "/nodesAttributes_1.txt", sep = "\t")

            # create igraph object with weights and node attributes
            nodI, nodE, icar = tp.CARNIVAL2igraph_convert(smpl_pd)
            icar, isoNode = tp.CARNIVALnode_attribute(icar, smpl_nAtt)

            # finde all subpaths from initial nodes to effectors
            paths = tp.pathfinder(icar, nodI, nodE)

            # add information to a dictionary
            netWord[network][sample] = [icar, nodI, nodE, isoNode, paths]

            # get some metrics
            aux = {"sample": sample, "cancer_tissue": ctype,
                     "network": network, "vertices": icar.vcount(),
                     "edges": icar.ecount(), "iniciators": len(nodI),
                     "effectors": len(nodE), "isolated": len(isoNode),
                     "subpaths": len(paths)}
            lstats.append(aux)

    cancer_tissue[ctype] = netWord
    with open('CMP_subpaths/' + ctype + '_netextend-omnipath_GWpkn.pkl', 'wb') \
                as ofile:
        pickle.dump(cancer_tissue[ctype], ofile)
    ofile.close()

stats = pd.DataFrame(lstats)
stats.to_csv('CMP_subpaths/netsum-omnipath_GWpkn.csv', header=True, index=None)
