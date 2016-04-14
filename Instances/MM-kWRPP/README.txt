
The files included in this webpage contain instances for Arc Routing Problems defined on 
Undirected, Mixed and Windy graphs.

For each type of problem there is a .zip file with the name of the problem to which the 
instances contained in it belong. For example, MGRP.ZIP contains instances of the Mixed 
General Routing Problem (MGRP). Also a .txt file with the same name is included, containing 
the optimal value (if known) for each instance. The first column in these files corresponds 
to the name of the instance. If the optimal value is known, it is shown in column 2 followed 
by the word "optimal". If not, a lower bound is given in column 2, followed by an upper 
bound (or, if there is no known upper bound, a "-"). A .txt file with the characteristics 
of the instances (numer of nodes, arcs, ....) is also included. These results have been 
obtained with a Branch and Cut algorithm for the WGRP, presented in "A Branch & Cut Algorithm 
for the Windy General Routing Problem" by Corberán, Plana and Sanchis (Networks, 2007).

With respect to the Min-Max K-vehicles WRPP, a zip file with the same name is included, 
containing 4 .txt files corresponding to the solutions with 2, 3, 4 and 5 vehicles. Each 
file contains the lower bound obtained at the root node (column 2) and the optimal value 
(if known) for each instance. If the optimal value is known, it is shown in column 3 
followed by the word "optimal". If not, an upper bound is given in column 4 and the final lower 
bound obtained after 1 hour CPU time is shown in column 5. A .txt 
file with the characteristics of the instances (numer of nodes, arcs, ....) is also 
included. These results have been obtained with an enhanced Branch and Cut algorithm for the 
Min Max K-WRPP, presented in "New facets and an enhanced Branch-And-Cut for the Min-Max 
K-vehicles Windy Rural Postman Problem" by Benavent, Corberán, Plana and Sanchis (2009).

Note that all the instances are represented as "windy instances", i.e. every link 
(arc or edge) has two associated costs (see below).

The files containing the instances have the following format:

NOMBRE : <name of instance>
COMENTARIO : <comment>
VERTICES : <number of vertices>
ARISTAS_REQ : <number of required edges>
ARISTAS_NOREQ : <number of non-required edges>
LISTA_ARISTAS_REQ :
(<u1>,<u2>) coste <cu1> <cu2>
.
.
.
LISTA_ARISTAS_NOREQ :
(<v1>,<v2>) coste <cv1> <cv2>
.
.
.



Entries <name of instance> and <comment> may be empty.
For each edge between nodes u1 and u2 there is a line 

	(<u1>,<u2>) coste <cu1> <cu2>

where cu1 is the cost of traversing the edge from u1 to u2 and cu2 is the cost of traversing it 
from u2 to u1. If the graph is not windy, both costs have the same value. If the edge is required, 
it wil be listed under "LISTA_ARISTAS_REQ :", and if it is non-required, it will be listed under 
"LISTA_ARISTAS_NOREQ :".

If the graph is mixed, an arc from node u1 to node u2 with cost cu will be represented by

	(<u1>,<u2>) coste <cu> 9999

or

	(<u2>,<u1>) coste 9999 <cu>