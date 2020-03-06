network_tools
##########

Collection of Python and R functions to run network-based analysis
in signed and directed networks.
Follow the jupyter notebook tutorial to see what you can do.

About the code
================================

Developed on:

- python 3.7.3 installed through ANACONDA [Clang 4.0.1 (tags/RELEASE_401/final)] in macOS Catalina version 10.15.1.

**Packages (Version)**
 * pandas (1.0.1)
 * networkx (2.4.0)
 * numpy (1.18.1)


- R 3.6.0 (2019-04-26) -- "Planting of a Tree" Platform: x86_64-apple-darwin15.6.0 (64-bit)

Notes
================================
1. Weight and sign are collapse in the same value: weight.
Negative weights indicate inhibition, while positive one are activations.

2. We use directed graphs to calculate network and node-based measurements.
This type of graphs only allow ONE interaction between two nodes.
As the weight contains also de sign,
if an interaction is repeated with opposite signs between 2 nodes
(e.g. A -| B; A -> B), only one weight can be kept.
Thus, the weight that is kept is the highest (keeping the sign).

Measurements
================================
The function __get_measurments__ calculate the follow network's features:

* Number of edges
* Number of nodes
* Density
* Avg betweenness centrality
* Avg degree centrality
* Avg eigenvector centrality
