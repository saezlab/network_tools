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
    :param weithSample string or pd.df: interaction file with weights
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

    # rename columns
    wm.columns = ['Node1', 'Sign', 'Node2', 'Weight']

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
    :param weithSample string or pd.df: interaction file with weights
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

    # rename columns
    wm.columns = ['Node1', 'Sign', 'Node2', 'Weight']

    iNodes = set(wm['Node1']) - set(wm['Node2'])
    eNodes = set(wm['Node2']) - set(wm['Node1'])

    return(iNodes, eNodes)

def get_measurments(DG, extended=False):
    """
    Get network and node-based measurments for a network as networkX.DiGraph
    :param DG networkx.DiGraph: directed graph from networkX
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

def adjacency2DG(adjaMTX):
    """
    Creates a networkX direct graph based on the adjacency matrix
    :param adjaMTX pd.df: pandas data frame with the network interactions
                            rows are the sources and columns the targets

    :return: DG: directed graph from networkX
    :rtype: networkx.DiGraph
    """
    import networkx as nx
    import pandas as pd

    sources = adjaMTX.index.tolist()
    targets = adjaMTX.columns.tolist()

    DG = nx.DiGraph()

    for source in sources:
        for target in targets:
            w = adjaMTX.loc[source, target]
            if w != 0:
                DG.add_edge(source, target, weight=w)

    return(DG)

def annotate_shared_nodes(nodAtt, adjaMTX, value="S"):
    """
    :param adjaMTX pd.df: pandas data frame with the network interactions
                            rows are the sources and columns the targets
    :param nodAtt pd.df: node atributes data frame
    :param value str: value to annotate the nodes in adjaMTX into nodAtt
    :return: nodAtt: same nodAtt pd.df with an exta column: core
    :rtype: pd.df
    """

    import pandas as pd

    sources = adjaMTX.index.tolist()
    targets = adjaMTX.columns.tolist()

    #rename nodAtt DataFrame
    nodAtt.columns = ['Node', 'ZeroAct', 'UpAct', 'DownAct', 'AvgAct', 'NodeType']

    #Add empty column to annotate
    nodAtt["core"] = ""

    for source in sources:
        for target in targets:
            w = adjaMTX.loc[source, target]
            if w != 0:
                nodAtt.loc[nodAtt['Node']==source, "core"] = value
                nodAtt.loc[nodAtt['Node']==target, "core"] = value

    return(nodAtt)

def annotate_shared_edges(edgAtt, adjaMTX, value="S"):
    """
    :param adjaMTX pd.df: pandas data frame with the network interactions
                            rows are the sources and columns the targets
    :param edgAtt pd.df: node atributes data frame
    :param value str: value to annotate the nodes in adjaMTX into edgAtt
    :return: edgAtt: same edgAtt pd.df with an exta column: core
    :rtype: pd.df
    """

    import pandas as pd

    sources = adjaMTX.index.tolist()
    targets = adjaMTX.columns.tolist()

    #rename edgAtt DataFrame
    edgAtt.columns = ['Source', 'Sign', 'Target', 'Weight']

    #Add empty column to annotate
    edgAtt["core"] = ""

    for source in sources:
        for target in targets:
            w = adjaMTX.loc[source, target]
            if w != 0:
                edgAtt.loc[(edgAtt['Source']==source) &
                           (edgAtt['Target']==target), "core"] = value

    return(edgAtt)
