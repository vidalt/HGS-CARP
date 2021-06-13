/*  ---------------------------------------------------------------------- //
    Hybrid Genetic Search for Arc Routing Problems -- HGS-CARP
    Copyright (C) 2016 Thibaut VIDAL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

// EPSILON GAP in the local search, to avoid numerical issues when accepting improvement
#define EPSILON_LS 0.001

#include "Noeud.h"
#include "SeqData.h"
#include <stdlib.h>
#include <stdio.h> 
#include <vector>
#include <list>
#include <math.h>
using namespace std ;

// Little structure for Ejection Chains
struct EC_element {
	
	// Node that is represented
	Noeud * myNode ;

	// Position for insertion (case where the node is a depot)
	Noeud * bestInsertionPlace ;

	// Cost of the labels arriving here
	double cost ;

	// Predecessor in the best path
	EC_element * pred ;

	// Number of customers nodes in the chain (in case of equality, long chains are priviledged).
	int nbCustNodesInChain ;

	// ID of the route to which the node is associated
	int routeID ;
};

// Local Search class
class LocalSearch
{

private:
	
	// acces to the data of the problem
	Params * params ;

	// is the data structure already initialized
	bool allAttributesSet ;

	// number of days in the problem
	int nbDays ;

	// pointer towards the first entry of the auxiliary data structures on nodes
	SeqData * seqdeb ;

	// Is the search finished
	bool rechercheTerminee ;

	// Access to the associated individual
	Individu * individu ;

	// Small pointers and variables which are used during the LS
	Noeud * noeudU ;
	Noeud * x ;
    Noeud * noeudV ;
	Noeud * y ;
	Route * routeU ;
	Route * routeV ;
	Vehicle * vehicleU ;
	Vehicle * vehicleV ;
	bool testingIncumbentPattern ; // little variable used during the local search (PI) to see if we are testing the pattern for which the customer is currently inserted
	bool deplacementIntraJour ; // flag raised if its possible to improve the location of a customer in its own day (PVRP and PCARP), used in PI
	bool firstLoop ; // are we in the first loop (for PI mutations)
	vector <SeqData *> myseqs ; // temporary vector to keep some pointers towards preprocessed data

	// shuffling procedure (for RI)
	void melangeParcours () ;

	// updating the list of moves (for RI)
	void updateMoves ();

	// little table which is used to sum up the costs of the moves
	vector < vector < double > > resultMoves ;

	// used during moves preprocessing, using the lower bounds
	vector < vector < bool > > shouldBeTested ;

	// say that some moves need to be tested again
	void reinitSingleDayMoves(Route * r);
	void reinitAllSingleDayMoves();

	// Functions for the evaluation of the classical moves (in RI)
	int interRouteGeneralInsert(); // Inter-Route General Insert (testing all moves together for a customer pair [i,j] enables to gain a few evaluations)
	int interRoute2Opt(); // 2-OPT* (without inversion of the routes)
	int interRoute2OptInv() ; // 2-OPT* (with inversion of the routes)
	int intraRouteGeneralInsertDroite(); // Intra-Route General Insert (testing all moves together for a customer pair [i,j] enables to gain a few evaluations)
    int intraRoute2Opt (); // Intra-Route 2-Opt

	// change the pattern of "client" if its possible to find a better pattern (PI)
	int searchBetterPattern (int client);
	void computeCoutInsertion(Noeud * client, int pattern) ; // subprocedures for PI
	void evalInsertClient (Route * R, Noeud * U, int pattern) ; // subprocedures for PI

public:

	// nbTotal Moves, for printouts
	int nbTotalRISinceBeginning ;
	int nbTotalPISinceBeginning ;

	// little data structure to choose the order of the exploration of the moves.
	// updated to contain only the visits in their respective days (PVRP and PCARP), and regularly shuffled
	vector < vector < int > > ordreParcours ;
	void addOP (int day, int client) ; // functions for updating the ordreParcours data structure
	void removeOP (int day, int client) ; // functions for updating the ordreParcours data structure

	// counters on the moves (for printing)
	int nbInterSwap ;
	int nbIntraSwap ;
	int nbInter2Opt ;
	int nbIntra2Opt ;
	int nbEjectionChains ;
	int nbEjectionChainsNodes ;

	// for each day, keeping a pointer routeVide[day] on the first empty route (to avoid iterating many identical empty routes) 
	// initially set to NULL
	vector < Route * > routeVide ;
	// setting the pointer for a given day
	void setRouteVide(int day);

	// Linked list of the customer (used to store the solution)
	// there is a pointer for the successor and the predecessor, as well as for the route and pre-processed data structures
	Noeud ** clients ; // Elements representing the customers
	Noeud ** depots ; // Elements representing the depots
	Noeud ** depotsFin ; // Sentinels at the end of the routes
	Route ** routes ; // Elements representing the routes

	// running the complete local search (RI-PI-RI)
	void runSearchTotal ();

	// main function for the RI procedure
	int mutationSameDay (int day);

	// main function for the PI procedure
	int mutationDifferentDay ();

	// says that all moves involving the customer "cli" in the day "day" have been tested.
	// to avoid testing them again as long as nothing has been modified
	void nodeTestedForEachRoute (int cli, int day);

	/* O(n²) VERSION OF EJECTION CHAINS */
	// Picks a random order for the routes
	// relies on some auxiliary graph structures to solve a shortest path sub-problem */
	int ec_nbRoutes ;
	vector < int > ec_nbElements ;
	vector < vector < EC_element > > ejectionGraph ;
	bool ejectionChains (int day); 	// returns true in case of success of ejection

	/* ROUTINES USED TO MODIFY THE SOLUTION, WHEN APPLYING A MOVE */

	void insertNoeud(Noeud * U, Noeud * V) ; // Inserts client U after V
	void removeNoeud(Noeud * U) ; // Removes client U
	void addNoeud(Noeud * U, Noeud * V) ; // Adds client U after V

	// To call intra-route moves on smaller sequences
	// if a limit on the size of the subsequence is set, it will cut the seqdata in smaller pieces
	// By default, there we use an upper bound,  such that all subsequences have been preprocessed (cf Params.h)
	inline void addSeqDataInPieces (Noeud * node, int length, int day);
	inline void addReverseSeqDataInPieces (Noeud * node, int length, int day); // the same with the reverse sequence. For 2-opt moves.

	/* PROCEDURE OF INSERTION OF THE MISSING NODES IN THE PIX CROSSOVER */
	void placeManquants ();

	// Constructor
	LocalSearch();

    // Constructor with all data structures
	LocalSearch(Params * params, Individu * individu);

	~LocalSearch(void);
};

#endif
