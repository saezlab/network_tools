"""
TITLE : Topological functions for network analysis
AUTHOR : Hernansaiz Ballesteros, Rosa. rosa.hernansaiz@bioquant.uni-heidelberg.de
DESCRIPTION : Functions to deal with comparisons between networks based on toplogy
LICENSE : GPL-v3
"""

def CARNIVAL2igraph_convert(dfSample):
    """
    Transform CARNIVAL's output read as pd.df into igraph.
    It also gives the origin and ending nodes of the graph.

    :param dfSample dict: {sample:pd.df} weightedModel_1.txt from CARNIVAL
    :return: nodAux (initial nodes), endAux (effectors), igrAux (igraph)
    :rtype: dict, dict, igraph object
    """

    import pandas as pd
    import igraph as ig
    import numpy as np

    nodAux = {} #sample-list origin Nodes
    endAux = {} #sample-list end Nodes
    igrAux = {} #sample-igraph
    for sample in dfSample.keys():
        Onodes = []
        mgraph = []
        for x in dfSample[sample].values:
            if len(x[np.where(x == "Perturbation")]) > 0 :
            #if x contains the 'Perturbation' edge, it's a root vertex
                if x[2] in dfSample[sample]['Node1'].values:
                #if that node is not isolated, it is an initial node
                    Onodes.append(x[2])
            else:
            #it is part of the graph
                mgraph.append({'source':x[0],'target':x[2],'weight':x[1]*x[3]})

        nodAux[sample] = Onodes
        dfaux = dfSample[sample][dfSample[sample].Node1 != 'Perturbation']
        endAux[sample] = [j for j in dfaux[~dfaux['Node2'].isin(dfaux['Node1'])].Node2.drop_duplicates()]
        mgraph = pd.DataFrame(mgraph)
        igrAux[sample] = ig.Graph.TupleList(mgraph.itertuples(index=False),
                                            weights=True, directed = True)

    return(nodAux, endAux, igrAux)

def pathfinder(igra, ori, fin):
    """
    Use pypath to find all paths between all initial (ori) and
    final (fin) nodes in a graph (igra)

    :param igra dict: {sample:igrapht} nodesAttributes_1.txt from CARNIVAL
    :param ori dict: ori{sample:list of initial nodes} nodesAttributes_1.txt from CARNIVAL
    :param fin dict: in{sample:list of final nodes} nodesAttributes_1.txt from CARNIVAL
    :return: topaux {sample:list of lists (subpaths)}
    :rtype: dict
    """
    import igraph as ig
    from pypath import main, data_formats
    import pickle

    topaux = {}
    pa = main.PyPath()
    for sample in igra.keys():
        topaux[sample] = pa.find_all_paths(start=ori[sample],
                            end=fin[sample],
                            attr = 'name',
                            maxlen=len(igra[sample].get_edgelist()),
                            graph = igra[sample])
        #       topaux[sample] = set(map(tuple,topaux[sample]))
        uniquelist = []
        for i in topaux[sample]:
            if i not in uniquelist:
                uniquelist.append(i)
        topaux[sample] = uniquelist

    return(topaux)
