This folder contains the instances proposed in "Vidal, T. (2017). Node, edge, arc routing and turn penalties: Multiple problems -- One neighborhood extension. Operations Research, 65(4), 992–1010".
It is available at: https://github.com/vidalt/HGS-CARP

These benchmark instances have been derived from original test problems from:
[1] Corberán, R. Martí, E. Martínez, and D. Soler, “The rural postman problem on mixed graphs with turn penalties,” J. Oper. Res. Soc., vol. 29, no. 7, pp. 887–903, 2002.
[2] L. Bach, G. Hasle, and S. Wøhlk, “A lower bound for the node, edge, and arc routing problem,” Comput. Oper. Res., vol. 40, no. 4, pp. 943–952, 2013.

The fleet is assumed to be limited to a number of vehicles specified in the instance definition. Each vehicle has a capacity of Q. Nodes are indexed from 1 to the number of nodes.
The instance file specifies, in turn, general information about the instance, the list of NODES, EDGES, ARCS, and TURNS.
For each node, edge, and arc, a boolean value “IS_REQUIRED” states whether a service is required at this location. Turns are identified as triplets of node indices I,J and K.
The specific format of the instance is described below.

Name: <Instance name>
#Vehicles: <Max. number of vehicles, -1 if unconstrained> 
Capacity: <Vehicle capacity Q>
Depot: <Index of depot node>
#Nodes: <number of nodes>
#Edges: <number of edges>
#Arcs: <number of arcs>
#Required N: <number of required nodes>
#Required E: <number of required edges>
#Required A: <number of required arcs>
#Nb-Turns: <number of turns>

----------NODES----------
INDEX QTY IS-REQUIRED X Y
<list of nodes>

----------EDGES----------
INDEX-I INDEX-J QTY IS-REQUIRED TR-COST
<list of edges>

-----------ARCS----------
INDEX-I INDEX-J QTY IS-REQUIRED TR-COST
<list of arcs>

----------TURNS----------
INDEX-I INDEX-J INDEX-K COST TYPE
<list of turns>

