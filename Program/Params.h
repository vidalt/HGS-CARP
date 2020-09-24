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

#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include "math.h"
#include <time.h>
#include <algorithm>
#include "Client.h"
#include "Vehicle.h"
using namespace std ;

// little function used to clear some arrays
template <class C> void FreeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); 
              it != cntr.end(); ++it ) {
        delete * it;
    }
    cntr.clear();
}

// Pre-definition, to allow compilation with self-references
class Params ;
class Vehicle ;

// A little auxiliary structure
struct pairB 
{
	int myInt ;
	int iCour ; 
	Params * myparams ;
};

class Params
{
public:

	/* ------------------------- PROBLEM DATA -------------------------- */

	// random seed
	int seed ;

	// path to the instance
	string pathToInstance ;

	// path to the solution
	string pathToSolution ;

	// path to the BKS (just to read the value and overwrite if needed)
	string pathToBKS ;

	// Problem type
	/*
	type =     // This lists the problems which can be solved with this algorithm
			   30 CARP (Capacitated Arc Routing Problem)
			   31 NEARP (Edges, Arcs and Nodes), also called MCGRP in the literature
			   32 PCARP (Periodic CARP)
			   33 MDCARP (Multi-depot CARP)
			   34 NEARP-TP (NEARP with Turn penalties) --> also set the flag "TURN_PENALTIES" when compiling
			   35 MM-kWRPP (Min-Max Windy Rural Postman Problem)
			   */
	int type ;
	bool multiDepot ; // is there multiple depots in the problem
	bool periodique ; // is there multiple periods in the problem
	bool isTurnPenalties ; // is there turn penalties

	int ar_NodesRequired ; // number of nodes which require a visit
	int ar_NodesNonRequired ; // number of other nodes
	int ar_EdgesRequired ; // number of edges which require a visit
	int ar_EdgesNonRequired ; // number of other edges
	int ar_ArcsRequired ; // number of arcs which require a visit
	int ar_ArcsNonRequired ; // number of other arcs

	// shortest paths between nodes of the original network
	// to keep the things clearer we don't use the "0 node", every index starts from one
	// the depot is among these nodes.
	void ar_InitializeDistanceNodes() ;
	vector < vector < double > > ar_distanceNodes ;

	/* CASE OF VEHICLE ROUTING PROBLEMS WITH TURN PENALTIES */
	// paths between arcs in the original network
	// the arc 0 is a fake arc which corresponds to the depot to himself
	// the distance between an arc i and arc j is the sum of distances of these arcs with the turn penalties on the way.
	int ar_nbArcsDistance ;
	int ar_nbTurns ;
	int ar_maxNbModes ;
	vector < Arc > ar_Arcs ;
	vector < vector < Arc * > > ar_correspondingArc ; // correspondingArc[i][j] returns a pointer to the associated arc between node i and j, otherwise, if non existing an error.
	vector < vector < double > > ar_distanceArcs ;
	void ar_computeDistancesArcs(); // computes the distance in the line graph

	// number of customers/services considered in the vehicle routing problem
	int nbClients ;

	// sum of customers x visit frequency (multi-period problems)
	int nbTotalServices ;

	// sum of demand quantity of the customers
	double totalDemand ;

	// number of days
	int nbDays ;
	int ancienNbDays ; // copy of the number of days

	// number of vehicles per depot
	int nbVehiculesPerDep ;

	// number of depots
	int nbDepots ;

	// for each day, sum of the capacity of the vehicles
	vector <double> dayCapacity ;

	// list of vehicles available for each day
	vector < vector < Vehicle > > ordreVehicules ;

	// number of vehicles available for each day
	vector <int> nombreVehicules ;

	// array containing the information of each separate client/service
	vector<Client> cli ;

	// travel time (was used for the CVRP) now its mainly used as an intermediate structure to compute the granular search proximity
	vector<vector<double > > timeCost ;

	// isCorrelated[i][j] returns true if and only if j is considered to be among the closest customers to i (granular search parameter)
	vector < vector <bool> > isCorrelated ;

	/* ----------------- METHOD PARAMETERS & OTHER DATA ---------------- */

	// Are we running an Iterated Local Search ?
	// In this case the behavior of the method is changed at several points
	bool isILS_general ;

	// Are we running a feasibility problem ?
	// In this case we would stop the search as soon as a feasible solution is found
	bool isSearchingFeasible ;

	// population size parameters
	int mu ; // Default 25
	int lambda ; // Default 40
	int el ; // Default 8

	// number of close individuals, taken into account in the distance measure (diversity management)
	int nbCountDistMeasure ; // Default 3

	// penalty coefficients (are adapted during the search)
	double penalityCapa ;
	double penalityLength ;

	// target of feasible individuals following LS
	double minValides ; // Default 0.25
	double maxValides ; // Default 0.35

	// number of close customers considered in RI (granular search)
	int granularity ; // Default 40

	// how much additional capacity consumption (multiplicator) allowed in Split
	double borne ; // Default 2

	// max size of a SeqData
	// The preprocessing effort could also be limited to O(n^4/3) using hierarchies as in Irnich 2008 (JOC)
	int sizeSD ; // Default 10

	/* ------------------------  PARSING ROUTINES  -------------------- */

	// incoming data stream
	ifstream fichier ;

	// setting the parameters of the method
	void setMethodParams () ;

	// get the data from the stream
	void preleveDonnees (string nomInstance) ;
	void ar_parseOtherLinesCARP(); // some sub-procedures when reading the various instance formats
	void ar_parseOtherLinesNEARP();
	void ar_computeDistancesNodes();
	void ar_parseOtherLinesNEARP_TP();
	int ar_tempIndexDepot ;
	int parsing_courNbArcs ;

	// reading a customer from the stream
	void getClient (int i);

	// builds the other data structures (granular search etc...)
	void calculeStructures () ;

	// sets the good patterns for a customer
	// part of the instance definition in the PCARP
	void setPatterns_PCARP(Client& myCli);

	// Used to initialize the data structures for multi-depot problems
	// Each depot is considered as a day (it works in the same way in the local search and all components of the method)
	void processDataStructuresMD () ;
	
	// shuffle the lists of closest customers
	void shuffleProches () ;

	// constructor
	Params(string nomInstance, string nomSolution, string nomBKS, int seedRNG, int type, int nbVeh, int nbDep, bool isSearchingFeasible);

	// destructor
	~Params(void);
};
#endif

