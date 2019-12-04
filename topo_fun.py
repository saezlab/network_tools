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

    :param dfSample pd.df: weightedModel_1.txt from CARNIVAL
    :return: nodAux (initial nodes), endAux (effectors), igrAux (igraph)
    :rtype: list, list, igraph object
    """

    import pandas as pd
    import igraph as ig
    import numpy as np

    Enodes = [] #list effectors
    Onodes = [] #list initial Nodes
    igr = [] #sample-igraph
    mgraph = []
    for line in dfSample.values:
        if len(line[np.where(line == "Perturbation")]) > 0 :
        #if x contains the 'Perturbation' edge, it's a root vertex
            if line[2] in dfSample['Node1'].values:
            #if that node is not isolated, it is an initial node
                Onodes.append(line[2])
        else:
        #it is part of the graph
            mgraph.append({'source':line[0],
                           'target':line[2],
                           'weight':line[1]*line[3]})

    dfaux = dfSample[dfSample.Node1 != 'Perturbation']
    Enodes = [j for j in dfaux[~dfaux['Node2'].isin(dfaux['Node1'])].Node2.drop_duplicates()]
    mgraph = pd.DataFrame(mgraph)
    igr = ig.Graph.TupleList(mgraph.itertuples(index=False),
                                        weights=True, directed = True)

    return(Onodes, Enodes, igr)

def CARNIVALnode_attribute(igra, nodeSample):
    """
    Add node status to vertices in the graph as attributes.
    It get up and down information.

    :param igra igraph: igraph object to add node activity.
    :param nodeSample pd.df: pd.df with node activities
    :return: igra with node attributes, isolated_nodes
    :rtype: igraph object, list
    """
    import numpy as np

    isolated_nodes = []
    for line in nodeSample.values:
        #if it is not information about the perturbation
        no_pert = len(line[np.where(line == "Perturbation")])
        if  no_pert <= 0 and line[1] != 100:
            try:
                igra.vs.find(line[0])
            except ValueError as e:
                isolated_nodes.append(line[0])
                continue
            v = int(str(igra.vs.find(line[0])).split(",")[1]) #id Vertex
            igra.vs[v]["up"] = line[2]
            igra.vs[v]["down"] = line[3]

    return igra, isolated_nodes

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
