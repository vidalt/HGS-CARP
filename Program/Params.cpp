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

#include "Params.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

void Params::setMethodParams()
{
	/* MAIN PARAMETERS OF THE METHOD */

	isILS_general = false ; // Are we running the ILS version of the code
	mu = 25 ; // Population size
	lambda = 40 ; // Number of individuals per generation
	el = 12 ; // Number of elite

	nbCountDistMeasure = 3 ; // Number of close individuals considered in the distance measure (diversity management)
	granularity = 40 ; // Restriction of the LS moves to 40 closest nodes
	minValides = 0.15 ; // Target proportion of feasible solution
	maxValides = 0.20 ; // Target proportion of feasible solution
	penalityCapa = 50 ; // Initial penalties (will evolve during the search)
	penalityLength = 50; // Initial penalties (will evolve during the search)

	// The ELS/ILS requires slightly different parameter setting to get the right number of children and solutions, as specified in Prins 2009
	if (isILS_general) 
	{ mu = 5 ; lambda = 5 ; el = 1 ; minValides = 0.6 ; maxValides = 0.7 ; }
}

void Params::preleveDonnees (string nomInstance)
{
	// Main method to read a problem instance
	vector <Vehicle> tempI ;
	double vc ;
	string contenu, useless2 ;
	nbTotalServices = 0 ;
	totalDemand = 0 ;
	multiDepot = false ;
	periodique = false ;
	isTurnPenalties = false ;
	char * myChars = new char[100]; 

	if (type == 30 || type == 33) // This is a standard CARP (can also be an experiment for MDCARP, when the number of depots is defined to be greater than 1)
	{
		if (type == 33)
			multiDepot = true ;

		// Reading all lines one by one
		getline(fichier, contenu); 
		getline(fichier, contenu);
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_NodesNonRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesNonRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> nbVehiculesPerDep ;

		// The instance only provides a lower bound on the necessary number of vehicles, per definition of the CARP, more vehicles are allowed (here we put one more)
		// Still, in all solutions the minimum number of vehicles turned out to be used
		nbVehiculesPerDep ++ ; 
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> vc ;
		getline(fichier, contenu);
		getline(fichier, contenu);		
		getline(fichier, contenu);
		getline(fichier, contenu);

		ar_NodesRequired = 0 ;
		ar_ArcsRequired = 0 ;
		ar_ArcsNonRequired = 0 ;

		ar_InitializeDistanceNodes();
	}
	else if (type == 31) // This is the standard NEARP (also called MCGRP)
	{
		getline(fichier, contenu);
		getline(fichier, contenu);
		int nodesTotal;
		int edgesTotal;
		int arcsTotal;
		fichier >> useless2 ;

		// If the number of vehicles per depot has been specified in the commandline
		// Otherwise we read it from the file (not always available, depending on the benchmark set)
		if (nbVehiculesPerDep == -1)
			fichier >> nbVehiculesPerDep ;
		else 
			fichier >> useless2 ;

		fichier >> useless2 ;
		fichier >> vc ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_tempIndexDepot ;
		fichier >> useless2 ;
		fichier >> nodesTotal ;
		fichier >> useless2 ;
		fichier >> edgesTotal ;
		fichier >> useless2 ;
		fichier >> arcsTotal ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_NodesRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_ArcsRequired ;
		ar_NodesNonRequired = nodesTotal - ar_NodesRequired ;
		ar_EdgesNonRequired = edgesTotal - ar_EdgesRequired ;
		ar_ArcsNonRequired  = arcsTotal - ar_ArcsRequired ;  
		getline(fichier, contenu);
		getline(fichier, contenu);

		ar_InitializeDistanceNodes();
	}
	else if (type == 32) // This is the PCARP (file format provided by N. Labadi)
	{	
		periodique = true ;
		// Reading all lines one by one
		getline(fichier, contenu); // gets a complete line
		fichier.getline(myChars,1000,'='); // This gets all the next characters until finding "="
		fichier >> ar_NodesNonRequired ;
		fichier.getline(myChars,1000,'=');
		fichier.getline(myChars,1000,'=');
		fichier >> ar_EdgesRequired ;
		fichier.getline(myChars,1000,'=');
		fichier.getline(myChars,1000,'=');
		if (nbVehiculesPerDep == -1) // Getting the number of vehicles per depot (in case was not specific in commandline)
			fichier >> nbVehiculesPerDep ;
		fichier.getline(myChars,1000,'=');
		fichier >> vc ;
		getline(fichier, contenu);
		getline(fichier, contenu);
		getline(fichier, contenu);
		getline(fichier, contenu);

		ar_NodesRequired = 0 ;
		ar_ArcsRequired = 0 ;
		ar_ArcsNonRequired = 0 ;
		ar_EdgesNonRequired = 0 ;

		ar_InitializeDistanceNodes();
	}
	else if (type == 34) // This is the NEARP-TP
	{
		isTurnPenalties = true ;

		getline(fichier, contenu);
		fichier >> useless2 ;
		fichier >> nbVehiculesPerDep ;
		fichier >> useless2 ;
		fichier >> vc ;
		fichier >> useless2 ;
		fichier >> ar_tempIndexDepot ;

		int totalNodes ;
		int totalEdges ;
		int totalArcs ;
		fichier >> useless2 ;
		fichier >> totalNodes ;
		fichier >> useless2 ;
		fichier >> totalEdges ;
		fichier >> useless2 ;
		fichier >> totalArcs ;
		fichier >> useless2 ;
		fichier >> ar_NodesRequired ;
		fichier >> useless2 ;
		fichier >> ar_EdgesRequired ;
		fichier >> useless2 ;
		fichier >> ar_ArcsRequired ;

		ar_NodesNonRequired = totalNodes - ar_NodesRequired ;
		ar_EdgesNonRequired = totalEdges - ar_EdgesRequired ;
		ar_ArcsNonRequired = totalArcs - ar_ArcsRequired ;

		fichier >> useless2 ;
		fichier >> ar_nbTurns ;

		getline(fichier, contenu);
		getline(fichier, contenu);

		// Trying to detect if something went wrong when reading the instance
		// These things could easily happen when specifying the wrong problem type for a given input data
		if (ar_nbArcsDistance < 0 || ar_nbArcsDistance > 1000000)
		throw string("PARSING ERROR : Incorrect number of arcs. A likely cause is the use of the wrong problem type");

		// initializing list of arcs (used for the modes)
		ar_nbArcsDistance = 1 + totalArcs + 2*totalEdges ;
		ar_Arcs = vector <Arc> (ar_nbArcsDistance) ;
		vector <Arc*> myTempVect = vector <Arc*> (this->ar_NodesNonRequired + this->ar_NodesRequired + 1, NULL) ;
		ar_correspondingArc = vector < vector <Arc*> > (this->ar_NodesNonRequired + this->ar_NodesRequired + 1, myTempVect);
	}
	else if (type == 35) // Min-Max WRPP, instance format of Angel Corberan
	{
		getline(fichier, contenu);
		getline(fichier, contenu);
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_NodesNonRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesRequired ;
		fichier >> useless2 ;
		fichier >> useless2 ;
		fichier >> ar_EdgesNonRequired ;
		getline(fichier, contenu);
		getline(fichier, contenu);
		vc = 10000 ; // There is no capacity limit in this problem.
		ar_NodesRequired = 0 ;
		ar_ArcsRequired = 0 ;
		ar_ArcsNonRequired = 0 ;
		ar_InitializeDistanceNodes();
	}
	else
		throw string ("Incorrect problem type");

	// Setting the other instance parameters
	if (type == 32) nbDays = 5 ; // PCARP instances are defined on 5 days
	else nbDays = 1 ; // Other Instances
	ancienNbDays = nbDays ;
	if (type != 33) nbDepots = 1 ; // for MDCARP instances, the number of depots was already specified in commandline
	nbClients = ar_ArcsRequired + ar_EdgesRequired + ar_NodesRequired ;

	// Printing some warning in case of wrong use of the instances and parameters
	if (nbVehiculesPerDep == -1) throw string ("WARNING : some type of instances do not specify a fleet size, please specify an upper bound manually using -veh XX in the commandline");
	#ifdef TURN_PENALTIES
	if (type != 34) throw string("When solving problem instances without turn penalties, please set the flag TURN_PENALTIES to false");
	#else
	if (type == 34) throw string("When solving problem instances with turn penalties, please set the flag TURN_PENALTIES to true");
	#endif

	// Trying to detect if something went wrong when reading the instance
	// These things could easily happen when specifying the wrong problem type for a given input data
	if (nbClients < 0 || nbClients > 1000000)
		throw string("ERROR WHEN READING : Number of services has not been correctly read. A very likely cause is the use of the wrong problem type for a given problem instance");

	// Building the list of vehicles
	ordreVehicules.push_back(tempI) ;
	nombreVehicules.push_back(0);
	dayCapacity.push_back(0);
	for (int kk=1 ; kk <= nbDays; kk ++)
	{
		ordreVehicules.push_back(tempI) ;
		dayCapacity.push_back(0);
		nombreVehicules.push_back(nbDepots*nbVehiculesPerDep);
		for (int i=0 ; i < nbDepots ; i++)
		{
			for (int j=0 ; j < nbVehiculesPerDep ; j++)
			{
				ordreVehicules.at(kk).push_back(Vehicle(i,1000000,vc)); // Duration constraint set to a high value
				dayCapacity.at(kk) += vc ;
			}
		}
	}

	// Reading the list of customers 
	// Not all instance formats follow this convention
	// Such that some dedicated parsing procedures are sometimes needed.
	cli = vector<Client> (nbDepots + nbClients);
	for (int i = 0 ; i < nbClients + nbDepots ; i++) 
		getClient (i);

	// The file formats for CARP, NEARP and NEARP-TP may require some distinct parsing routines due to different formats
	if (type == 30 || type == 33 || type == 35) // CARP or MM-kWRPP
	{
		ar_parseOtherLinesCARP();
		ar_computeDistancesNodes(); // We only need the distance between endpoints

		bool has_components = false;
		for (int i=0; i < ar_distanceNodes.size(); ++i) {
			for(int j=0; j < ar_distanceNodes.at(i).size(); ++j) {
				if (ar_distanceNodes.at(i).at(j) == 1.e20) {
					has_components = true;
				}
			} 
			if (has_components) {
				throw string("The Provided Network has components, please make sure the network is connected");
			}
		}
	}
	else if (type == 31) // NEARP
	{
		ar_parseOtherLinesNEARP();
		ar_computeDistancesNodes();  // We only need the distance between endpoints
	}
	else if (type == 32) // PCARP
	{
		ar_computeDistancesNodes(); // We only need the distance between endpoints
	}
	else if (type == 34) // NEARP-TP
	{
		ar_parseOtherLinesNEARP_TP();
		ar_computeDistancesArcs();  // For the NEARP-TP, we rely on the distance between arcs (computed in the line graph)
	}
}

bool compPredicate(pairB i, pairB j) 
{ 
	return (i.myparams->timeCost.at(i.myInt).at(i.iCour) < i.myparams->timeCost.at(j.myInt).at(i.iCour)) ;
}

void Params::calculeStructures () 
{
	// Initializing other structures of the search
	vector < vector < bool > > tempB ;
	vector <bool> tempB2 ;
	vector <pairB> myVector ;
	pairB myPair ;
	myPair.myparams = this ;

	// Initialization of the table of "correlation" (granular search restriction)
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		isCorrelated.push_back(tempB2);
		for (int j=0 ; j < nbClients + nbDepots ; j++)
			isCorrelated.at(i).push_back(i < nbDepots || j < nbDepots);
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		cli.at(i).sommetsVoisins.clear();
		cli.at(i).sommetsVoisinsAvant.clear();
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		myVector.clear();
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myPair.myInt = j ;
			myPair.iCour = i ;
			if ( i != j ) myVector.push_back(myPair);
		}

		// For each customer, sorting the list of other customers
		std::sort(myVector.begin(),myVector.end(),compPredicate);

		// And keeping only the closest ones
		for (int j=0 ; j < min(nbClients,granularity) ; j++)
		{
			cli.at(i).sommetsVoisinsAvant.push_back(myVector.at(j).myInt);
			cli.at(myVector.at(j).myInt).sommetsVoisins.push_back(i);
			isCorrelated.at(myVector.at(j).myInt).at(i) = true ;
		}
	}

	for ( int i=0 ; i < nbDepots + nbClients ; i++ )
	{
		cli.at(i).computeVisitsDyn(nbDays,ancienNbDays);
		cli.at(i).computeJourSuiv(nbDays,ancienNbDays);
	}
}

void Params::getClient (int i)
{
	// Reading/Initializing the data for each customer
	string tempstring ;
	pattern p ;
	char * myChars = new char[100]; 
	int t_arc, t_from, t_to, t_qty, t_trav, t_col, t_inv, t_freq ;
	bool isNewService ;

	Client& myCli = cli.at(i);

	myCli.custNum = i ;
	myCli.freq = 1 ;
	myCli.serviceDuration = 0 ;

	if (type == 30 || type == 33) // CARP instances (also used when running MDCARP experiments, which only involve to define more depots)
	{
		myCli.ar_nodesExtr0 = -1 ;
		myCli.ar_nodesExtr1 = -1 ;
		myCli.serviceDuration = 0 ;

		if (i < nbDepots)
		{
			myCli.demand = 0 ;
			myCli.ar_nodeType = AR_DEPOT ;
			myCli.ar_serviceCost01 = 0. ;
			myCli.ar_serviceCost10 = 0. ;
			myCli.freq = 0 ;
		}
		else
		{
			fichier >> tempstring ;
			fichier >> myCli.ar_nodesExtr0 ;
			fichier >> tempstring ;
			fichier >> myCli.ar_nodesExtr1 ;
			fichier >> tempstring ;
			fichier >> tempstring ;
			fichier >> myCli.ar_serviceCost01 ;
			myCli.ar_serviceCost10 = myCli.ar_serviceCost01 ;

			ar_distanceNodes.at(myCli.ar_nodesExtr0).at(myCli.ar_nodesExtr1) = myCli.ar_serviceCost01 ;
			ar_distanceNodes.at(myCli.ar_nodesExtr1).at(myCli.ar_nodesExtr0) = myCli.ar_serviceCost01 ;

			fichier >> tempstring ;
			fichier >> myCli.demand ;
			totalDemand += myCli.demand ;
			myCli.ar_nodeType = AR_CLIENT_EDGE ;
		}
	}
	else if (type == 31) // NEARP instances (we don't really read the file here but at least initialize the structures)
	{
		myCli.ar_nodesExtr0 = -1 ;
		myCli.ar_nodesExtr1 = -1 ;
		myCli.demand = 0 ;
		if (i == 0)
		{
			myCli.ar_nodeType = AR_DEPOT ;
			myCli.ar_nodesExtr0 = ar_tempIndexDepot ;
			myCli.ar_nodesExtr1 = ar_tempIndexDepot ;
			myCli.freq = 0 ;
		}
		else if (i <= ar_NodesRequired)
			myCli.ar_nodeType = AR_CLIENT_NODE ;
		else if (i <= ar_NodesRequired + ar_EdgesRequired)
			myCli.ar_nodeType = AR_CLIENT_EDGE ;
		else if (i <= ar_NodesRequired + ar_EdgesRequired + ar_ArcsRequired)
			myCli.ar_nodeType = AR_CLIENT_ARC ;
		else 
			cout << "Error in parsing nbNodes" << endl ;

		myCli.ar_serviceCost01 = 0. ;
		myCli.ar_serviceCost10 = 0. ;
	}
	else if (type == 32) // PCARP instances
		// Here we need to make a slight trick, because the lines in the file contain both arc directions
		// (two lines per service)
		// and not necessarily consecutive
		// Thus we will test if the arc corresponds to an edge that has been already encountered or not
	{
		if (i == 0)
		{
			getline(fichier, tempstring); // removing the line
			myCli.demand = 0 ;
			myCli.ar_serviceCost01 = 0. ;
			myCli.ar_serviceCost10 = 0. ;
			myCli.ar_nodeType = AR_DEPOT ;
			myCli.ar_nodesExtr0 = 1 ; // The depot is always index 1 in these instances
			myCli.ar_nodesExtr1 = 1 ;
			myCli.freq = 0 ;
		}
		else
		{
			isNewService = false ;
			while (!isNewService) // looking for the next new service
			{
				fichier.getline(myChars,1000,'=');
				fichier >> t_arc ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_from ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_to ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_qty ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_trav ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_col ;
				fichier.getline(myChars,1000,'=');
				fichier >> t_inv ;
				fichier.getline(myChars,1000,';');
				fichier.getline(myChars,1000,'=');
				fichier >> t_freq ;
				getline(fichier, tempstring); // finishing this line
				if (t_inv > t_arc) isNewService = true ; // We found a new service
			}

			myCli.ar_nodesExtr0 = t_from ;
			myCli.ar_nodesExtr1 = t_to ;
			myCli.ar_serviceCost01 = t_col ;
			myCli.ar_serviceCost10 = t_col ;
			ar_distanceNodes.at(myCli.ar_nodesExtr0).at(myCli.ar_nodesExtr1) = t_trav ;
			ar_distanceNodes.at(myCli.ar_nodesExtr1).at(myCli.ar_nodesExtr0) = t_trav ;
			myCli.demand = t_qty ;
			totalDemand += 7.0 * myCli.demand ;
			myCli.ar_nodeType = AR_CLIENT_EDGE ;
			myCli.freq = t_freq ;
		}
	}
	else if (type == 34) // NEARP-TP instances (we don't really read the file here but at least initialize the structures)
	{
		myCli.demand = 0 ;
		if (i == 0)
		{
			myCli.ar_nodeType = AR_DEPOT ;
			myCli.freq = 0 ;
		}
		else if (i <= ar_NodesRequired)
			myCli.ar_nodeType = AR_CLIENT_NODE ;
		else if (i <= ar_NodesRequired + ar_EdgesRequired)
			myCli.ar_nodeType = AR_CLIENT_EDGE ;
		else if (i <= ar_NodesRequired + ar_EdgesRequired + ar_ArcsRequired)
			myCli.ar_nodeType = AR_CLIENT_ARC ;
		else 
			cout << "Error when parsing nodes" << endl ;
	}
	else if (type == 35) // Min-Max WRPP
	{
		myCli.serviceDuration = 0 ;
		myCli.demand = 0 ; // No capacity constraints in this problem

		if (i == 0)
		{
			myCli.ar_nodesExtr0 = 1 ;
			myCli.ar_nodesExtr1 = 1 ;
			myCli.ar_nodeType = AR_DEPOT ;
			myCli.ar_serviceCost01 = 0. ;
			myCli.ar_serviceCost10 = 0. ;
			myCli.freq = 0 ;
		}
		else
		{
			fichier >> tempstring ;
			fichier >> myCli.ar_nodesExtr0 ;
			fichier >> tempstring ;
			fichier >> myCli.ar_nodesExtr1 ;
			fichier >> tempstring ;
			fichier >> tempstring ;
			fichier >> myCli.ar_serviceCost01 ;
			fichier >> myCli.ar_serviceCost10 ;

			// setting the distance between the nodes (here its symmetric)
			ar_distanceNodes.at(myCli.ar_nodesExtr0).at(myCli.ar_nodesExtr1) = myCli.ar_serviceCost01 ;
			ar_distanceNodes.at(myCli.ar_nodesExtr1).at(myCli.ar_nodesExtr0) = myCli.ar_serviceCost10 ;

			myCli.ar_nodeType = AR_CLIENT_EDGE ;
		}
	}

	// Setting the pattern information (on which day the customer can be visited)
	if (type != 32 || i == 0)
	{
		p.dep = 0 ;
		p.cost = 0 ;
		p.pat = 1 ;
		myCli.visits.push_back(p);
		myCli.visitsOrigin.push_back(p);
	}
	else 
	{
		// Need to define the visit patterns for the PCARP
		// The conventions are defined in Chu et al 2006 
		setPatterns_PCARP(myCli);
	}

	if (i >= 1) nbTotalServices += myCli.freq ;
	delete [] myChars ;
}

void Params::setPatterns_PCARP(Client& myCli)
{
	vector<int> listPatternsInstance = {16,8,4,2,1,20,18,17,10,9,5,26,25,22,21,19,13,11,30,29,27,23,15,31};
	pattern p ;
	p.dep = 0 ;
	p.cost = 0 ;
	int deb, end ;

	if (myCli.freq == 1)
	{deb = 0 ; end = 4 ;}
	else if (myCli.freq == 2)
	{deb = 5 ; end = 10 ;}
	else if (myCli.freq == 3)
	{deb = 11 ; end = 17 ;}
	else if (myCli.freq == 4)
	{deb = 18 ; end = 22 ;}
	else if (myCli.freq == 5)
	{deb = 23 ; end = 23 ;}
	else 
		throw ("ERROR FREQ PCARP");

	for (int i=deb ; i <= end ; i++) // The patterns [0..4] are acceptable per instance definition
	{
		p.pat = listPatternsInstance.at(i) ;
		myCli.visits.push_back(p);
		myCli.visitsOrigin.push_back(p);
	}
}

Params::Params(string nomInstance, string nomSolution, string nomBKS, string nomDistMat, int seedRNG, int type, int nbVeh, int nbDep, bool isSearchingFeasible):type(type), nbVehiculesPerDep(nbVeh), nbDepots(nbDep), isSearchingFeasible(isSearchingFeasible)
{
	// Main constructor of Params
	pathToInstance = nomInstance ;
	pathToSolution = nomSolution ;
	pathToBKS = nomBKS ;
    pathToDistanceMatrix = nomDistMat;
	borne = 2.0 ;
	sizeSD = 10 ;

	seed = seedRNG;
	if (seed == 0) // using the time to generate a seed when seed = 0 
		srand((unsigned int)time(NULL));
	else 
		srand(seed);

	// Opening the instance file
	fichier.open(nomInstance.c_str());

	// Reading the instance file
	if (fichier.is_open())
		preleveDonnees (nomInstance);
	else 
		throw string(" Impossible to find instance file ");

	// Setting the method parameters
	setMethodParams();

	// If its a problem with multiple, we consider it as an equivalent problem with multiple periods
	// one for each depot
	if (multiDepot) processDataStructuresMD();

	// Computing the other data structures
	calculeStructures();	
}

Params::~Params(void)
{
}

void Params::shuffleProches () 
{
	int temp,temp2 ;

	// Shuffling the list of close customers for each customer
	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli.at(i).sommetsVoisins.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli.at(i).sommetsVoisins.size() - a1) ;
			temp =  cli.at(i).sommetsVoisins.at(a1) ;
			cli.at(i).sommetsVoisins.at(a1) = cli.at(i).sommetsVoisins.at(temp2);
			cli.at(i).sommetsVoisins.at(temp2) = temp ;
		}
	}

	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli.at(i).sommetsVoisinsAvant.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli.at(i).sommetsVoisinsAvant.size() - a1) ;
			temp =  cli.at(i).sommetsVoisinsAvant.at(a1) ;
			cli.at(i).sommetsVoisinsAvant.at(a1) = cli.at(i).sommetsVoisinsAvant.at(temp2);
			cli.at(i).sommetsVoisinsAvant.at(temp2) = temp ;
		}
	}
}

void Params::ar_InitializeDistanceNodes()
{
	if (ar_NodesNonRequired+ar_NodesRequired  < 0 || ar_NodesNonRequired+ar_NodesRequired  > 1000000)
		throw string("ERROR WHEN READING : Number of nodes has not been correctly read. A very likely cause is the use of the wrong problem type for a given problem instance");

	// Build the ar_distanceNodes data structures
	vector <double> myTemp = vector <double> (ar_NodesNonRequired+ar_NodesRequired+1) ; 
	for (int i=0 ; i <= ar_NodesNonRequired+ar_NodesRequired ; i++)
		myTemp.at(i) = 1.e20 ;

	ar_distanceNodes.clear();
	for (int i=0 ; i <= ar_NodesNonRequired+ar_NodesRequired ; i++)
		ar_distanceNodes.push_back(myTemp);
}

void Params::ar_computeDistancesNodes()
{
	cout << "Computing distance Nodes" << endl;
	for (int ii = 1 ; ii <= ar_NodesNonRequired + ar_NodesRequired ; ii++)
		ar_distanceNodes.at(ii).at(ii) = 0 ;

	cout << "Progress: 0.000%\r" << flush;

	// simple application of the Floyd Warshall algorithm
	double last_perc = 0.0;

	for (int k=1 ; k <= ar_NodesNonRequired + ar_NodesRequired ; k++)
	{
		double current_perc = (float)k*100/(ar_NodesNonRequired+ar_NodesRequired);
		if (current_perc - last_perc > 1.0) {
			cout << "Progress: " << current_perc << "%                        \r" << flush;
			last_perc = current_perc;
		}
		for (int i=1 ; i <= ar_NodesNonRequired + ar_NodesRequired ; i++)
		{
			for (int j=1 ; j <= ar_NodesNonRequired + ar_NodesRequired ; j++)
			{
				if (ar_distanceNodes.at(i).at(k) + ar_distanceNodes.at(k).at(j) < ar_distanceNodes.at(i).at(j))
					ar_distanceNodes.at(i).at(j) = ar_distanceNodes.at(i).at(k) + ar_distanceNodes.at(k).at(j) ;
			}
		}
	}

	// Then, we would still like to include some distance information between services.
	// The distance between two services is the minimum distance between the closest endpoints of the edge
	// This is used by the granular search
	double d ;
	timeCost = vector<vector<double> >(nbClients + nbDepots+1);
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		timeCost.at(i) = vector<double> (nbClients + nbDepots+1);
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			d = min(min(ar_distanceNodes.at(cli.at(i).ar_nodesExtr0).at(cli.at(j).ar_nodesExtr0),
				ar_distanceNodes.at(cli.at(i).ar_nodesExtr0).at(cli.at(j).ar_nodesExtr1)),
				min(ar_distanceNodes.at(cli.at(i).ar_nodesExtr1).at(cli.at(j).ar_nodesExtr0),
				ar_distanceNodes.at(cli.at(i).ar_nodesExtr1).at(cli.at(j).ar_nodesExtr1)));

			timeCost.at(i).at(j) = d ;
		}
	}
	cout << "Computing distance nodes finished" << endl;
}

void Params::ar_computeDistancesArcs()
{
    if (!Params::readDistanceMatrix(ar_distanceArcs,pathToDistanceMatrix)) {
        // computes the distance in the line graph (for the turn penalties)
        // simple application of the Floyd Warshall algorithm
        cout << "Computing distance args" << endl;
        cout << "Progress: 0.000%\r" << flush;
        double last_perc = 0.0;
        auto start = chrono::system_clock::now();
    
        for (int k=0 ; k < ar_nbArcsDistance ; k++)
        {
            double current_perc = (float)k*100/(ar_nbArcsDistance);
            if (current_perc-last_perc > 1.0) {
                chrono::duration<double> elapsed = chrono::system_clock::now() - start;
                double remainig_seconds = elapsed.count()/current_perc*(100.0-current_perc);
                
                cout << "Progress: " << current_perc << "%    [" << (int) remainig_seconds/60 << ":" << setfill('0') << setw(2) << (int) remainig_seconds % 60 << "] remaining                       \r" << flush;
                last_perc = current_perc;
            }
            for (int i=0 ; i < ar_nbArcsDistance ; i++)
            {
                for (int j=0 ; j < ar_nbArcsDistance  ; j++)
                {
                    if (ar_distanceArcs.at(i).at(k) + ar_distanceArcs.at(k).at(j) < ar_distanceArcs.at(i).at(j))
                        ar_distanceArcs.at(i).at(j) = ar_distanceArcs.at(i).at(k) + ar_distanceArcs.at(k).at(j) ;
                }
            }
        }
    
        Params::writeDistanceMatrix(ar_distanceArcs,pathToDistanceMatrix);
    }
    

	// Then, we would still like to include some distance information between services.
	// The distance between two services is the minimum distance between one mode of each service
	// This is used by the granular search
	double myDistanceMin ;
	double myTemp ;

	timeCost = vector<vector<double> > (nbClients + nbDepots);
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		timeCost.at(i) = vector<double> (nbClients + nbDepots);
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myDistanceMin = 1.e20 ;
			for (int ii=0 ; ii < (int)cli.at(i).ar_Modes.size() ; ii++)
			{
				for (int jj=0 ; jj < (int)cli.at(j).ar_Modes.size() ; jj++)
				{
					myTemp = ar_distanceArcs.at(cli.at(i).ar_Modes.at(ii)->indexArc).at(cli.at(j).ar_Modes.at(jj)->indexArc) ;
					if (myTemp < myDistanceMin)
						myDistanceMin = myTemp ;
				}
			}
			timeCost.at(i).at(j) = myDistanceMin ;
		}
	}
}

void Params::ar_parseOtherLinesCARP()
{
	// Parsing routine for CARP
	string contenu;
	string useless;
	int startNode ;
	int endNode ;
	double myCost ;

	getline(fichier, contenu);
	if (ar_EdgesNonRequired  > 0) getline(fichier, contenu);

	for (int k=0 ; k < ar_EdgesNonRequired ; k++)
	{
		fichier >> useless ;
		fichier >> startNode ;
		fichier >> useless ;
		fichier >> endNode ;
		fichier >> useless ;
		fichier >> useless ;
		fichier >> myCost ;
		ar_distanceNodes.at(startNode).at(endNode) = myCost ;
		if (type != 35) // not a "windy problem"
			ar_distanceNodes.at(endNode).at(startNode) = myCost ;
		else
			fichier >> ar_distanceNodes.at(endNode).at(startNode) ;
	}

	// at the end we need to set the depot locations 
	// (except for the WRPP, where it was already defined as node 1 per convention of the instances)
	if (type == 30)
	{
		fichier >> useless ;
		fichier >> useless ;
		fichier >> startNode ; // in CARP instances, this information is included in the file

		cli.at(0).ar_nodesExtr0 = startNode ;
		cli.at(0).ar_nodesExtr1 = startNode ;
	}
	else if (type == 33)
	{
		cli.at(0).ar_nodesExtr0 = 1 ;
		cli.at(0).ar_nodesExtr1 = 1 ;
		// Case with two depots (see instance specs in Kansou and Yassine 2009)
		if (nbDepots == 2)
		{
			cli.at(1).ar_nodesExtr0 = ar_NodesNonRequired ;
			cli.at(1).ar_nodesExtr1 = ar_NodesNonRequired ;
		}
		else if (nbDepots >= 3) // Case with three or more depots (see instance specs in Kansou and Yassine 2009)
		{
			for (int d=1 ; d < nbDepots ; d++)
			{
				cli.at(d).ar_nodesExtr0 = d*(int)(ar_NodesNonRequired/(nbDepots-1)) ;
				cli.at(d).ar_nodesExtr1 = d*(int)(ar_NodesNonRequired/(nbDepots-1)) ;
			}
		}
	}
}

void Params::ar_parseOtherLinesNEARP()
{
	// Parsing routine for NEARP
	string contenu;
	string useless;
	double myTravelCost ;
	int iCour = 1 ;
	int idNode ;
	string myString ;
	string myStringTemp ;
	int side0 ;
	int side1 ;

	//  REQUIRED NODES
	getline(fichier, contenu);
	for (int k=0 ; k < ar_NodesRequired ; k++)
	{
		fichier >> myStringTemp ;
		myString = myStringTemp.substr(1,myStringTemp.size());
		idNode = atoi(myString.c_str());
		cli.at(iCour).ar_nodesExtr0 = idNode ;
		cli.at(iCour).ar_nodesExtr1 = idNode ;
		fichier >> cli.at(iCour).demand ;
		totalDemand += cli.at(iCour).demand ;
		fichier >> cli.at(iCour).ar_serviceCost01 ;
		cli.at(iCour).ar_serviceCost01 = 0 ; // to remove the base cost from the results.
		cli.at(iCour).ar_serviceCost10 = 0 ; // to remove the base cost from the results.
		iCour ++ ;
	}
	getline(fichier, contenu);
	if (ar_NodesRequired > 0) getline(fichier, contenu);

	// REQUIRED EDGES
	getline(fichier, contenu);
	for (int k=0 ; k < ar_EdgesRequired ; k++)
	{
		fichier >> useless ;
		fichier >> cli.at(iCour).ar_nodesExtr0 ;
		fichier >> cli.at(iCour).ar_nodesExtr1 ;

		// setting the distance between the nodes (here its symmetric because it's an edge)
		fichier >> myTravelCost ;
		ar_distanceNodes.at(cli.at(iCour).ar_nodesExtr0).at(cli.at(iCour).ar_nodesExtr1) = myTravelCost ;
		ar_distanceNodes.at(cli.at(iCour).ar_nodesExtr1).at(cli.at(iCour).ar_nodesExtr0) = myTravelCost ;

		fichier >> cli.at(iCour).demand ;
		totalDemand += cli.at(iCour).demand ;
		fichier >> cli.at(iCour).ar_serviceCost01 ;
		cli.at(iCour).ar_serviceCost01 = myTravelCost ;
		cli.at(iCour).ar_serviceCost10 = myTravelCost ;
		iCour ++ ;
	}
	getline(fichier, contenu);
	if (ar_EdgesRequired > 0) getline(fichier, contenu);

	// NON-REQUIRED EDGES
	getline(fichier, contenu);
	for (int k=0 ; k < ar_EdgesNonRequired ; k++)
	{
		fichier >> useless ;
		fichier >> side0 ;
		fichier >> side1 ;
		fichier >> myTravelCost ;
		ar_distanceNodes.at(side0).at(side1) = myTravelCost ;
		ar_distanceNodes.at(side1).at(side0) = myTravelCost ;
	}
	getline(fichier, contenu);
	if (ar_EdgesNonRequired > 0) getline(fichier, contenu);

	// REQUIRED ARCS
	getline(fichier, contenu);
	for (int k=0 ; k < ar_ArcsRequired ; k++)
	{
		fichier >> useless ;
		fichier >> cli.at(iCour).ar_nodesExtr0 ;
		fichier >> cli.at(iCour).ar_nodesExtr1 ;

		// setting the distance between the nodes
		fichier >> myTravelCost ;
		ar_distanceNodes.at(cli.at(iCour).ar_nodesExtr0).at(cli.at(iCour).ar_nodesExtr1) = myTravelCost ;

		fichier >> cli.at(iCour).demand ;
		totalDemand += cli.at(iCour).demand ;
		fichier >> cli.at(iCour).ar_serviceCost01 ;
		cli.at(iCour).ar_serviceCost01 = myTravelCost ;
		cli.at(iCour).ar_serviceCost10 = 1.e20 ; // This is an arc
		iCour ++ ;
	}
	getline(fichier, contenu);
	if (ar_ArcsRequired > 0) getline(fichier, contenu);

	// NON-REQUIRED ARCS
	getline(fichier, contenu);
	for (int k=0 ; k < ar_ArcsNonRequired ; k++)
	{
		fichier >> useless ;
		fichier >> side0 ;
		fichier >> side1 ;
		fichier >> myTravelCost ;
		ar_distanceNodes.at(side0).at(side1) = myTravelCost ;
	}

	if (iCour != ar_ArcsRequired + ar_EdgesRequired + ar_NodesRequired + 1)
		cout << "Problem number of arcs, edges and nodes when parsing" << endl ;
}

void Params::ar_parseOtherLinesNEARP_TP()
{
	// Parsing routine for NEARP-TP
	string contenu ;
	double uselessDbl ;

	int p_qty ;
	bool p_isReq ;
	int p_indexI ;
	int p_indexJ ;
	int p_indexK ;
	double p_traversal ;
	double p_costTurn ;

	parsing_courNbArcs = 0 ;
	int parsing_courService = 0 ;


	/* FIRST SPECIFY THE DEPOT CHARACTERISTICS */
	ar_Arcs.at(0).indexArc = 0 ;
	ar_Arcs.at(0).cost = 0 ;
	ar_Arcs.at(0).nb_Turns = 0 ;
	ar_Arcs.at(0).nodeBegin = ar_tempIndexDepot ;
	ar_Arcs.at(0).nodeEnd = ar_tempIndexDepot ;
	cli.at(0).ar_Modes.push_back(&ar_Arcs.at(0));
	cli.at(0).ar_nbModes = 1 ;
	parsing_courNbArcs ++ ;
	parsing_courService ++ ;

	/* THEN READ THE CHARACTERISTICS OF THE NODE SERVICES */
	getline(fichier, contenu);
	getline(fichier, contenu);

	for (int i=0 ; i < this->ar_NodesNonRequired + this->ar_NodesRequired ; i++)
	{
		fichier >> p_indexI ;
		fichier >> p_qty ;
		fichier >> p_isReq ;
		if (p_isReq)
		{
			cli.at(parsing_courService).ar_nodesExtr0 = p_indexI ;
			cli.at(parsing_courService).demand = p_qty ;
			totalDemand += cli.at(parsing_courService).demand ;
			parsing_courService ++ ;
		}
		fichier >> uselessDbl ;
		fichier >> uselessDbl ;
	}

	/* THEN READ THE CHARACTERISTICS OF THE EDGE SERVICES */
	getline(fichier, contenu);
	getline(fichier, contenu);
	getline(fichier, contenu);
	if (this->ar_NodesNonRequired + this->ar_NodesRequired > 0)
		getline(fichier, contenu);

	for (int i=0 ; i < this->ar_EdgesNonRequired + this->ar_EdgesRequired ; i++)
	{
		fichier >> p_indexI ;
		fichier >> p_indexJ ;
		fichier >> p_qty ;
		fichier >> p_isReq ;
		fichier >> p_traversal ;
		ar_Arcs.at(parsing_courNbArcs).indexArc = parsing_courNbArcs ;
		ar_Arcs.at(parsing_courNbArcs).nodeBegin = p_indexI ;
		ar_Arcs.at(parsing_courNbArcs).nodeEnd = p_indexJ ;
		ar_Arcs.at(parsing_courNbArcs).cost = p_traversal ;
		ar_correspondingArc.at(ar_Arcs.at(parsing_courNbArcs).nodeBegin).at(ar_Arcs.at(parsing_courNbArcs).nodeEnd) = &ar_Arcs.at(parsing_courNbArcs) ;
		parsing_courNbArcs ++ ;
		ar_Arcs.at(parsing_courNbArcs).indexArc = parsing_courNbArcs ;
		ar_Arcs.at(parsing_courNbArcs).nodeBegin = p_indexJ ;
		ar_Arcs.at(parsing_courNbArcs).nodeEnd = p_indexI ;
		ar_Arcs.at(parsing_courNbArcs).cost = p_traversal ;
		ar_correspondingArc.at(ar_Arcs.at(parsing_courNbArcs).nodeBegin).at(ar_Arcs.at(parsing_courNbArcs).nodeEnd) = &ar_Arcs.at(parsing_courNbArcs) ;
		parsing_courNbArcs ++ ;

		if (p_isReq)
		{
			cli.at(parsing_courService).demand = p_qty ;
			totalDemand += cli.at(parsing_courService).demand ;
			cli.at(parsing_courService).ar_Modes.push_back(&ar_Arcs.at(parsing_courNbArcs-2));
			cli.at(parsing_courService).ar_Modes.push_back(&ar_Arcs.at(parsing_courNbArcs-1));
			cli.at(parsing_courService).ar_nbModes = 2 ;
			parsing_courService ++ ;
		}
	}

	/* THEN READ THE CHARACTERISTICS OF THE ARCS SERVICES */
	getline(fichier, contenu);
	getline(fichier, contenu);
	getline(fichier, contenu);
	if (this->ar_EdgesNonRequired + this->ar_EdgesRequired > 0)
		getline(fichier, contenu);

	for (int i=0 ; i < this->ar_ArcsNonRequired + this->ar_ArcsRequired ; i++)
	{
		fichier >> p_indexI ;
		fichier >> p_indexJ ;
		fichier >> p_qty ;
		fichier >> p_isReq ;
		fichier >> p_traversal ;
		ar_Arcs.at(parsing_courNbArcs).indexArc = parsing_courNbArcs ;
		ar_Arcs.at(parsing_courNbArcs).nodeBegin = p_indexI ;
		ar_Arcs.at(parsing_courNbArcs).nodeEnd = p_indexJ ;
		ar_Arcs.at(parsing_courNbArcs).cost = p_traversal ;
		ar_correspondingArc.at(ar_Arcs.at(parsing_courNbArcs).nodeBegin).at(ar_Arcs.at(parsing_courNbArcs).nodeEnd) = &ar_Arcs.at(parsing_courNbArcs) ;
		parsing_courNbArcs ++ ;

		if (p_isReq)
		{
			cli.at(parsing_courService).demand = p_qty ;
			totalDemand += cli.at(parsing_courService).demand ;
			cli.at(parsing_courService).ar_Modes.push_back(&ar_Arcs.at(parsing_courNbArcs-1)); // HERE is -1
			cli.at(parsing_courService).ar_nbModes = 1 ;
			parsing_courService ++ ;
		}
	}


	/* BUILDING THE ARC DISTANCE STRUCTURES */

	vector <double> myTemp = vector <double> (ar_nbArcsDistance) ; 
	for (int i=0 ; i < ar_nbArcsDistance ; i++)
		myTemp.at(i) = 1.e20 ;

	ar_distanceArcs.clear();
	for (int i=0 ; i < ar_nbArcsDistance ; i++)
		ar_distanceArcs.push_back(myTemp);

	for (int i=0 ; i < ar_nbArcsDistance ; i++)
		ar_distanceArcs.at(i).at(i) = 0 ; 

	// distance of any connected pair of edges with a connected at the depot location should have a turn penalty of 0
	for (int i=0 ; i < ar_nbArcsDistance ; i++)
	{
		for (int j=0 ; j < ar_nbArcsDistance ; j++)
		{
			if (ar_Arcs.at(i).nodeEnd == ar_tempIndexDepot && ar_Arcs.at(j).nodeBegin == ar_tempIndexDepot)
				ar_distanceArcs.at(i).at(j) = ar_Arcs.at(i).cost ;
		}
	}
    
    // distance should be without turn penalty as default. will be overwritten below when parsing the turns
    for (int i=0 ; i < ar_nbArcsDistance ; i++)
    {
        for (int j=0; j < ar_nbArcsDistance ; j++) {
            if (ar_Arcs.at(i).nodeEnd == ar_Arcs.at(j).nodeBegin)
                ar_distanceArcs.at(i).at(j) = ar_Arcs.at(i).cost;
        }
    }

	/* PARSING THE TURNS */
	getline(fichier, contenu);
	getline(fichier, contenu);
	getline(fichier, contenu);
	if (this->ar_ArcsNonRequired + this->ar_ArcsRequired > 0)
		getline(fichier, contenu);

	for (int i=0 ; i < this->ar_nbTurns ; i++)
	{
		fichier >> p_indexI ;
		fichier >> p_indexJ ;
		fichier >> p_indexK ;
		fichier >> p_costTurn ;
		fichier >> contenu ;

		int corrArc1 = ar_correspondingArc.at(p_indexI).at(p_indexJ)->indexArc ;
		int corrArc2 = ar_correspondingArc.at(p_indexJ).at(p_indexK)->indexArc ;

		ar_distanceArcs.at(corrArc1).at(corrArc2) = ar_correspondingArc.at(p_indexI).at(p_indexJ)->cost + p_costTurn ;
	}

	/* SPECIFYING THE MODES FOR THE NODE DELIVERIES */
	for (int i=1 ; i <= this->ar_NodesRequired ; i++)
	{
		cli.at(i).ar_nbModes = 0 ;
		for (int k=0 ; k < (int)ar_Arcs.size() ; k++)
		{
			if (ar_Arcs.at(k).nodeEnd == cli.at(i).ar_nodesExtr0)
			{
				cli.at(i).ar_Modes.push_back(&ar_Arcs.at(k));
				cli.at(i).ar_nbModes ++ ;
			}
		}
	}

	/* IDENTIFYING THE MAXIMUM NUMBER OF MODES FOR A DELIVERY */
	ar_maxNbModes = 0 ;
	for (int i=0 ; i <= this->nbClients ; i++)
		if (cli.at(i).ar_nbModes > ar_maxNbModes)
			ar_maxNbModes = cli.at(i).ar_nbModes ;
}


// effectue la conversion de MDPVRP en PVRP
void Params::processDataStructuresMD () 
{
	int nbPat ;
	vector <Vehicle> temp ;
	pattern p ;
	vector < vector <Vehicle> > ordreVehiculesAncien ;

	// we change the number of days
	// j1 - j2 - j3
	// da - db - dc
	// --> j1a - j2a - j3a - j1b - j2b - j3b - j1c - j2c - j3c
	nbDays = nbDays*nbDepots ;

	ordreVehiculesAncien = ordreVehicules ;
	ordreVehicules.clear();
	nombreVehicules.clear();
	dayCapacity.clear();

	ordreVehicules.push_back(temp);
	nombreVehicules.push_back(0);
	dayCapacity.push_back(0);
	for (int k=1 ; k <= nbDays ; k++ )
	{
		ordreVehicules.push_back(temp);
		nombreVehicules.push_back(nbVehiculesPerDep) ;
		dayCapacity.push_back(0);
		for (int d=0 ; d<nbVehiculesPerDep ; d++ )
		{
			ordreVehicules.at(k).push_back(ordreVehiculesAncien.at( (k-1)%ancienNbDays + 1 ).at(((k - 1)/ancienNbDays) * nbVehiculesPerDep + d ));
			dayCapacity.at(k) += ordreVehicules.at(k).at(d).vehicleCapacity ;
		}
	}

	// we update the patterns
	for (int i = nbDepots ;  i < nbDepots + nbClients ; i++ )
	{
		nbPat = (int)cli.at(i).visitsOrigin.size();
		cli.at(i).visits.clear();

		for (int j=0 ; j < nbDepots ; j++)
		{
			p.dep = j ;
			p.cost = 0 ;
			for ( int pat = 0 ; pat < nbPat ; pat ++ )
			{
				p.pat = cli.at(i).visitsOrigin.at(pat).pat ;
				cli.at(i).visits.push_back(p) ;
			}
		}
		cli.at(i).visitsOrigin = cli.at(i).visits ;
	}	
}

void Params::writeDistanceMatrix (const vector<vector<double> >& matrix, string filename) {
    if (filename == "")
        return;
        
    ofstream dfile;
    try {
        dfile.open(filename.c_str());
        if (dfile.good()) {
            for (int k=0 ; k < matrix.size() ; k++) {
                for (int i=0; i < matrix.at(k).size() ; i++) {
                    dfile << matrix.at(k).at(i) << " ";
                }
                dfile << endl;
            }
            dfile.close();
        }
    } catch (ofstream::failure e) {
        cout << "Could not write distance matrix" << endl;
    }
}

bool Params::readDistanceMatrix (vector<vector<double> >& matrix, string filename) {
    if (filename == "")
        return false;
        
    ifstream rmat;
    string continu;
    bool success = false;
    
    try {
        rmat.open(filename.c_str());
        if (rmat.good()) {
            cout << "Found matrix. restoring" << endl;
            cout << "Progress: 0.000%\r" << flush;
            double last_perc = 0.0;
            auto start = chrono::system_clock::now();
        
            for (int k=0 ; k < matrix.size() ; k++) {
                double current_perc = (float)k*100/(matrix.size());
                if (current_perc-last_perc > 1.0) {
                    chrono::duration<double> elapsed = chrono::system_clock::now() - start;
                    double remainig_seconds = elapsed.count()/current_perc*(100.0-current_perc);
                
                    cout << "Progress: " << current_perc << "%    [" << (int) remainig_seconds/60 << ":" << setfill('0') << setw(2) << (int) remainig_seconds % 60 << "] remaining                       \r" << flush;
                    last_perc = current_perc;
                }
            
                for (int i=0; i < matrix.at(k).size() ; i++) {
                    rmat >> matrix.at(k).at(i);
                }
                getline(rmat,continu);
            }
            rmat.close();
            success = true;
        }
    } catch (ifstream::failure e) {
        cout << "Could not load distance matrix" << endl;
    }
    return success;
}
