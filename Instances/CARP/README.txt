========================================================================
|                                                                      |
|                           CARPLIB                                    |
|                                                                      |
|                     Version of November 5, 2005                      |
|                                                                      |
|                     Dpto. Estadistica e I.O.                         |
|                     Universitat de Valencia                          |
|                                                                      |
========================================================================


In the current version all problems are defined on a graph not complete. 
We consider two different costs on the edges with positive demand: 
servicing costs and traversing costs. Servicing cost are the cost of traversing
and servicing the corresponding edge. Traversing cost denotes the cost of
traversing, without servicing, the edge. 

==========================================================================
|                                                                        |
| 1. THE FILE FORMAT                                                     |
|                                                                        |
==========================================================================

Basically, each file consists of two parts: a specification part and a 
data part. The specification part contains information on the file format 
and on its contents. The data part contains explicit data.

==========================================================================
| 1.1. The specification part                                            |
==========================================================================

All entries in this section consist of lines of the form

-----------------------------------------------
<keyword> : <value>
-----------------------------------------------

where <keyword> denotes an alphanumerical keyword and <value> denotes 
alphanumerical or numerical data. The terms <string>, <integer> and <real>
denote character string, integer or real data, respectively. Integer and
real numbers are given in free format. All the keywords are in spanish language.

Below we give a list of all available keywords.

-----------------------------------------------
NOMBRE : <string>
-----------------------------------------------

Used as an identification of the data file (name of the instance).

---------------------------------------------------
COMENTARIO : <integer> (<string>)
---------------------------------------------------

An additional comment on the data. Usually, an upper bound.

-----------------------------------------------
VERTICES : <integer>
-----------------------------------------------

Specifies the number of nodes.

-----------------------------------------------
ARISTAS_REQ : <integer>
-----------------------------------------------

Specifies the number of edges with positive demand (required edges).

-----------------------------------------------
ARISTAS_NOREQ : <integer>
-----------------------------------------------

Specifies the number of edges with zero demand (non required edges).

-----------------------------------------------
VEHICULOS : <integer>
-----------------------------------------------

Specifies the number of vehicles available.

-----------------------------------------------
CAPACIDAD : <integer>
-----------------------------------------------

Specifies the capacity of a truck.

----------------------------------------------------
TIPO_COSTES_ARISTAS : EXPLICITOS
----------------------------------------------------

Edge costs are listed explicitly in the corresponding 
section.

-----------------------------------------------
COSTE_TOTAL_REQ : <integer>
-----------------------------------------------

Specifies the sum of the servicing costs of all the required edges.



==========================================================================
| 1.2. The data part                                                     |
==========================================================================

The data are given in the corresponding
data sections which follow the specification part. Each data section is
started with a corresponding keyword. The length of the section is either 
implicitly known from the format specification or the section is 
terminated by special end-of-section terminators.

-----------------------------------------------
 LISTA_ARISTAS_REQ : 
-----------------------------------------------

Edges with positive demand are given in this section. Each line is of the form

( <integer> , <integer> )  coste  <integer>   demanda  <integer> 

The pair of integers between parenthesis give the node indices of the endpoints of the edge. The other two integers specify the traversing cost and the demand, respectively.

-----------------------------------------------
DEPOSITO :    <integer>
-----------------------------------------------

Contains the index of depot node (always node 1). 



==========================================================================
|                                                                        |
| 2.  FILES                                                              |
|                                                                        |
==========================================================================

The data files that are available in CARP directory in the 
above format correspond to those used in the papers:

 1.  BCCM directory:

     "The Capacitatted Arc Routing Problem. Lower Bounds", E. Benavent, V. Campos,
     A. Corberan and M. Mota, Networks 22, (1992) pp. 669-690.

     Thus instance val1A.dat correspond to problem 1.A in that paper, and so on.  

 2.  GDB directory:
     
     "Computational Experiments with Algorithms for a class of Routing Problem",  
     B. L. Golden, J. S. DeArmon and E. K. Baker, Computers and Operations Research 10,
     (1983) pp. 47-59.

     This reference contains 25 instances, but instances 8 and 9 were removed because 
     they contained inconsistences. 

     Also, in some instances, we have interchanged the original depot with node 1   
     (instances 8, 9, 10, 11, and 13).

 3.  KSHS directory:
     
     "An exact algorithm for the Capacitatted Arc Routing Problem using Parallel  
     Branch and Bound method", M. Kiuchi, Y. Shinano, R. Hirabayashi and Y. Saruwatari,
     Abstracts of the "1995 Spring National Conference of the Oper. Res. Soc. of Japan"
     pp. 28-29 (in Japanese). 
 
 4.  EGLESE directory:
     
     These instances are based on some instances used in the following two works:

     "Vehicle Routeing for Winter Gritting", L.Y.O. Li (1992). PH. D. Thesis, Dept. of 
     Management Science, Lancaster University.

     "An Interactive Algorithm for Vehicle Routeing for Winter-Gritting", L.Y.O. Li and 
     R.W. Eglese, Journal of the Operational Research Society 47, (1996) pp. 217-228.

     A more detailed explanation on how the instances has been generated can be found in:

     "A Cutting Plane Algorithm for the Capacitated Arc Routing Problem", J.M. Belenguer 
     and E. Benavent, Computers and Operations Research 30 (5), (2003) pp. 705-728.