"""
TITLE : topo_fun
AUTHOR : Hernansaiz Ballesteros, Rosa.
            rosa.hernansaiz@bioquant.uni-heidelberg.de
DESCRIPTION : topological functions to deal with CARNIVAL's results.

LICENSE : GPL-v3
"""

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
