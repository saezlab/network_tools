"""
TITLE : topo_fun
AUTHOR : Hernansaiz Ballesteros, Rosa.
            rosa.hernansaiz@bioquant.uni-heidelberg.de
DESCRIPTION : topological functions to deal with CARNIVAL's results.

LICENSE : GPL-v3
"""

def carnival2directGraph(weigthSample, inverse=True, verbose=False):
    """
    Creates a dirct graph from CARNIVAL's weight's file.
    :param weithSample string or pd.df: weightedModel_1.txt from CARNIVAL
    :param inverse boolean: type of carnival analysis. inverse=True (default)
    :param verbose boolean: give printed information about repeated edges.
                            False (default)
    :return: DG (direct graph)
    :rtype: networkX graph
    """
    import pandas as pd
    import networkx as nx

    if type(weigthSample) == str:
        wm = pd.read_csv(weigthSample, sep='\t') # read file
    else:
        wm = weigthSample

    if inverse:
        wm = wm[wm.Node1 != 'Perturbation'] # remove perturbation node
        wm.reset_index(inplace = True, drop=True) # reset indexes

    dupEdge = wm[wm[['Node1', 'Node2']].duplicated(keep=False)] # get duplicates

    if dupEdge.shape[0] > 0:

        idxDE = dupEdge[['Node1', 'Node2']].groupby(dupEdge[['Node1', 'Node2']].columns.tolist()).apply(lambda x: tuple(x.index)).tolist()

        if verbose:
            print('Number of repeated edges:', dupEdge.shape[0]/2)

        # get the highest value for the repeated edge (keeping sign)
        index2drop = []
        for dupIndex in idxDE:
            weight = [wm.iloc[dupIndex[0]]['Weight'], wm.iloc[dupIndex[1]]['Weight']]
            i = weight.index(min(weight))
            index2drop.append(dupIndex[i])

        wm.drop(labels=[dupIndex[i]], axis='index', inplace=True)

    # create direct graph
    DG = nx.DiGraph()

    auxL = []
    for row in wm.values:
    	auxL.append((row[0], row[2], int(row[1])*int(row[3])))

    DG.add_weighted_edges_from(auxL)

    return(DG)

def get_initiators_effectors(weigthSample, inverse=True):
    """
    Get initiators (initial nodes) and effectors (ending nodes)
    from CARNIVAL weightedModel_1
    :param weithSample string or pd.df: weightedModel_1.txt from CARNIVAL
    :param inverse boolean: type of carnival analysis. inverse=True (default)
    :return: iNodes (initial nodes), eNodes (effectors)
    :rtype: set
    """

    import pandas as pd

    if type(weigthSample) == str:
        wm = pd.read_csv(weigthSample, sep='\t') # read file
    else:
        wm = weigthSample

    if inverse:
        wm = wm[wm.Node1 != 'Perturbation'] # remove perturbation node
        wm.reset_index(inplace = True) # reset indexes

    iNodes = set(wm['Node1']) - set(wm['Node2'])
    eNodes = set(wm['Node2']) - set(wm['Node1'])

    return(iNodes, eNodes)

def get_measurments(DG, extended=False):
    """
    Get network and node-based measurments for a network as networkX.DiGraph
    :param DG networkx.DiGraph: directed graph frim networkX
    :param extended boolean: False default. get extra measures (diameter,
                                                 avg closeness centrality,
                                                 avg eccentricity,
                                                 avg eigenvector centrality)
    :return: measure with values: number nodes, edges, network's density,
                                  avg betweenness centrality,
                                  avg degree centrality
                                  and avg eigenvector centrality
    :rtype: dict
    """
    import networkx as nx
    import numpy as np

    measure = {}
    measure['nNodes'] = nx.number_of_nodes(DG)
    measure['nEdges'] = nx.number_of_edges(DG)
    measure['density'] = nx.density(DG)
    bc = nx.betweenness_centrality(DG)
    measure['avg betweenness centrality'] = np.mean(list(bc.values()))
    dc = nx.degree_centrality(DG)
    measure['avg degree centrality'] = np.mean(list(dc.values()))
    if extended:
        measure['diameter'] = nx.diameter(DG)
        cc = nx.closeness_centrality(DG)
        measure['avg closeness centrality'] = np.mean(list(cc.values()))
        ec = nx.eccentricity(DG)
        measure['avg eccentricity'] = np.mean(list(ec.values()))
        ev = nx.eigenvector_centrality(DG)
        measure['avg eigenvector centrality'] = np.mean(list(ev.values()))
    return(measure)

def CARNIVAL2igraph_convert(weigthSample):
    """
    Transform CARNIVAL's output read as pd.df into igraph.
    It also gives the origin and ending nodes of the graph.

    :param weithSample pd.df: weightedModel_1.txt from CARNIVAL
    :return: Onodes (initial nodes), Enodes (effectors), igra (igraph)
    :rtype: list, list, igraph object
    """

    import pandas as pd
    import igraph as ig
    import numpy as np

    Enodes = [] #list effectors
    Onodes = [] #list initial Nodes
    igra = [] #sample-igraph
    mgraph = []
    for line in weigthSample.values:
        if len(line[np.where(line == "Perturbation")]) > 0 :
        #if x contains the 'Perturbation' edge, it's a root vertex
            if line[2] in weigthSample['Node1'].values:
            #if that node is not isolated, it is an initial node
                Onodes.append(line[2])
        else:
        #it is part of the graph
            mgraph.append({'source':line[0],
                           'target':line[2],
                           'weight':line[1]*line[3]})

    dfaux = weigthSample[weigthSample.Node1 != 'Perturbation']
    Enodes = [j for j in dfaux[~dfaux['Node2'].isin(dfaux['Node1'])].Node2.drop_duplicates()]
    mgraph = pd.DataFrame(mgraph)
    igra = ig.Graph.TupleList(mgraph.itertuples(index=False),
                                        weights=True, directed = True)

    return(Onodes, Enodes, igra)

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

def pathfinder(igra, Onodes, Enodes):
    """
    Use pypath to find all paths between all initial (Onodes) and
    final (Enodes) nodes in a graph (igra)

    :param igra igraph: igraph object
    :param Onodes list: initial nodes
    :param Enodes list: effector nodes
    :return: uniquelist contains all subpath form initial to effector nodes
    :rtype: list of lists
    """
    import igraph as ig
    from pypath import main, data_formats

    pa = main.PyPath()

    topaux = pa.find_all_paths(start=Onodes,
                    end=Enodes,
                    attr = 'name',
                    maxlen=len(igra.get_edgelist()),
                    graph = igra)
    #check to remove duplicates
    uniquelist = []
    for i in topaux:
        if i not in uniquelist:
            uniquelist.append(i)

    return(uniquelist)

def getTopology(Cfile):
    """
    Get all information form CARNIVAL:
      - network in igraph format (icar)
      - unique initial nodes (nodI) (for inverse CARNIVAL)
      - unique effectors (nodE)
      - isolated nodes (isoNode)
      - subpaths (paths)
    :param Cfile str: directory to CARNIVAL results
    :return: contains icar, nodI, nodE, isoNode, paths
    :rtype: list of lists
    """
    import pandas as pd

    # carnival outputs
    smpl_pd = pd.read_csv(Cfile + "/weightedModel_1.txt", sep = "\t")
    smpl_nAtt = pd.read_csv(Cfile + "/nodesAttributes_1.txt", sep = "\t")

    # create igraph object with weights and node attributes
    nodI, nodE, icar = CARNIVAL2igraph_convert(smpl_pd)
    icar, isoNode = CARNIVALnode_attribute(icar, smpl_nAtt)

    # finde all subpaths from initial nodes to effectors
    paths = pathfinder(icar, nodI, nodE)

    return [icar, nodI, nodE, isoNode, paths]
