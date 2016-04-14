The Instance Definition Format and cost definition:

Below is a description of the format of the text file that defines each problem instance. A fleet of identical vehicles is assumed. The number of vehicles is given (-1 if unconstrained). Each vehicle has a capacity of Q. Nodes are indexed from 1 to the number of nodes. A required node i has demand q_i. A required edge or arc (i,j) has demand q_ij. Each edge is indexed from 1 to the number of edges. Required edges are listed separately before the non-required edges are listed. The same goes for arcs. The index of the depot node is explicitly given.

By 'traversal cost' c_ij of an edge or arc (i,j) we understand the cost for just traversing it (no service). Nodes have no traversal cost. The deadheading cost , i.e., the cost for pure movement between required items must be calculated by a shortest path algorithm, as there may be multiple arcs/edges between two given nodes, and the triangle inequality may not hold.

Service costs are given for all required nodes (s_i), edges and arcs (s_ij). These are the cost of servicing the demand in the required items.

The optimal value given here (if known) is the total traversal and deadheading cost of the solution, not including total service cost, which is constant. Thus for required nodes, no cost is added, and for required edges and arcs, only the traversal cost is added. Deadheading costs are determined as explained above. Upper and lower bounds are calculated according to this definition. This is also the reporting convention that has been and should be used in reports and papers. The instance definitions and best known bounds can all be found on the individual benchmark pages.

_____________________________________________________________________
 

Name:          <Instance name>
Optimal value: <Optimal value, -1 if unknown>
#Vehicles:     <Max. number of vehicles, -1 if unconstrained>
Capacity:      <Vehicle capacity Q>
Depot:         <Index of depot node>
#Nodes:        <number of nodes>
#Edges:        <number of edges>
#Arcs:         <number of arcs> 
#Required N:   <number of required nodes>
#Required E:   <number of required edges>
#Required A:   <number of required arcs>

% Required nodes:  Ni q_i s_ij
NODE INDEX, DEMAND, SERVICE COST

% Required edges: Ek i j q_ij c_ij s_ij  
EDGE INDEX, FROM NODE, TO NODE, TRAVERSAL COST, DEMAND, SERVICE COST 

% Non-required edges: NrEl i j c_ij 
EDGE INDEX, FROM NODE, TO NODE, TRAVERSAL COST

% Required arcs: Ar i j q_ij c_ij
ARC INDEX, FROM NODE, TO NODE, TRAVERSAL COST, DEMAND, SERVICE COST 

% Non-required arcs: NrAs i j c_ij
ARC INDEX, FROM NODE, TO NODE, TRAVERSAL COST