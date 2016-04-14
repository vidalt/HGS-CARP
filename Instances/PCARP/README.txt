Instance Format:

Name = Name of the instance
nbnodes= Number of nodes in the graph
nbarcs= Number of arcs
nbtasks= Number or required arcs
nbnon_req= Number of non-required arcs
nbtrucks= Number of vehicles
capacity= Capacity of the vehicles
total_cost= Sum of costs of required arcs (value from the original CARP instances, not adapted to the PCARP)
Depot= Index of the depot

!List of internal arcs:
arc= index of the arc,
From= initial extremity,
to= final extremity,
Qty= demand,
Trav= cost for traversing,
Col= cost for collecting,
inv= index of the opposed arc,
succ= indices of successor arcs, 
freq= frequency of services,
minspac= can be ignored for the PCARP,
maxspac= can be ignored for the PCARP.

For each instance, the allowed patterns (day combinations) and customer demands should be 
set as specified in "Chu, F., Labadi, N. & Prins, C., 2006. A Scatter Search for the periodic capacitated arc routing problem. European Journal of Operational Research, 169(2), pp.586–605"

