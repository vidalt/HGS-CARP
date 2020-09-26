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

#include <stdlib.h>
#include <stdio.h> 
#include <string>
#include "Genetic.h"
#include "commandline.h"

using namespace std;

int main (int argc, char *argv[])
{
	bool minFleetSize ;
	bool minMaxTour ;
	int nbpop = 0;
	vector < Population * > populationTab ;
	vector < Params * > mesParametresTab ;
	Population * lastPop = NULL ;
	Population * population2 ; 
	Population * population ; 
	Params * mesParametres ;
	Params * mesParametres2 ;
	clock_t nb_ticks_allowed;
	int veh ;
	double distConstraint ;
	int nbOverallLoop = 0 ;
	cout << endl ;

	try
	{
		// Reading the commandline
		commandline c(argc, argv);

		if (!c.is_valid())
			throw string("Commandline could not be read, Usage : gencarp instance -type problemType [-t cpu-time] [-sol solutionPath]  [-dmat distanceMatrix] [-s seed] [-veh nbVehicles] [-dep nbDepots]");

		minFleetSize = (c.get_type() == 32) ; // For the PCARP, we need to minimize fleet size as first objective, then minimize distance as second objective
		minMaxTour = (c.get_type() == 35) ; // For the MM-kWRPP, we need to minimize the length of the maximum route

		/* CLASSIC CASE OF OPTIMIZATION, BASED ON DISTANCE : for the CVRP, CARP, NEARP, MDCARP... */
		/* THIS IS THE MAIN START OF THE PROGRAM */
		if (!minFleetSize && !minMaxTour)
		{
			// Number of clock ticks allowed for the program
			nb_ticks_allowed = c.get_cpu_time() * CLOCKS_PER_SEC;

			// initialisation of the Parameters
			mesParametres = new Params(c.get_path_to_instance(),c.get_path_to_solution(),c.get_path_to_BKS(),c.get_path_to_distmat(),c.get_seed(),c.get_type(),c.get_nbVeh(),c.get_nbDep(),false) ;

			// Running the algorithm
			population = new Population(mesParametres) ;
			Genetic solver(mesParametres,population,nb_ticks_allowed,true);

			solver.evolve(20000,1); // First parameter (20000) controls the number of iterations without improvement before termination

			// Printing the solution
			population->ExportBest(c.get_path_to_solution());
			population->ExportBKS(c.get_path_to_BKS());

			delete population ;
			delete mesParametres ;
			cout << endl ;
			return 0 ;
		}


		/* SOME PROBLEMS CONSIDERED IN THE PAPER INVOLVE ANOTHER OBJECTIVE, such as fleet size minimization, or minimization of the maximum tour */
		/* THIS IS DONE HERE BY RUNNING ITERATIVELY THE ALGORITHM with a decreasing fleet or distance constraint */
		// fleet size minimization (minFleetSize = true) -- applying the algorithm with a decreasing fleet size, as long as a feasible solution is found
		// or minimization of the maximum tour (minMaxTour = true) -- applying the algorithm with a decreasing tour duration constraint
		else
		{
			veh = c.get_nbVeh(); // start with an upper bound on the number of vehicles
			nb_ticks_allowed = c.get_cpu_time() * CLOCKS_PER_SEC;
			distConstraint = 1.e30 ; // or with a permissive distance constraint
			bool validExist = true;
			while (validExist) // A feasible solution has been found, we can continue to decrease (either the number of vehicles or the distance constraint, depending on the case)
			{	
				// Setting the parameters of the next problem
				mesParametresTab.push_back(new Params(c.get_path_to_instance(),c.get_path_to_solution(),c.get_path_to_BKS(),c.get_path_to_distmat(),c.get_seed(),c.get_type(),veh,c.get_nbDep(),true)) ;
				nbpop = (int)mesParametresTab.size() ;
				nbOverallLoop ++ ; // counting the number of subproblems which have been resolved

				// For safety, to evacuate any chance of infinite loop and printout.
				// No considered instances should lead to more than 10000 overall modifications of the fleet size or distance constraint
				if (nbOverallLoop >= 10000)
					throw string ("Fleet or distance minimization, too many overall loops, there must be a problem, aborting the run");

				// Setting the distance constraint (only effective for the MM-kWRPP)
				for (int v=0 ; v < mesParametresTab[nbpop-1]->nbVehiculesPerDep ; v++) 
					mesParametresTab[nbpop-1]->ordreVehicules[1][v].maxRouteTime = distConstraint ;
				
				if (minMaxTour && minFleetSize) throw string("This program was not designed to optimize jointly the fleet size and length of the maximum tour");

				// Keeping the current penalty values
				if (lastPop != NULL) 
				{
					mesParametresTab[nbpop-1]->penalityCapa = lastPop->params->penalityCapa ; 
					mesParametresTab[nbpop-1]->penalityLength = lastPop->params->penalityLength ;
				}

				// Constructing the new population
				populationTab.push_back (new Population(mesParametresTab[nbpop-1])) ;

				// Adding the individuals found in previous iterations to help the search to start
				if (nbpop >= 2) populationTab[nbpop-1]->addAllIndividus(populationTab[nbpop-2]);
				if (nbpop >= 3) populationTab[nbpop-1]->addAllIndividus(populationTab[nbpop-3]);
				// Solving
				Genetic solver(mesParametresTab[nbpop-1],populationTab[nbpop-1],nb_ticks_allowed,true);
				cout << "######### GA evolution ######### : " << "| FLEET SIZE : " << veh << " | DIST CONSTRAINT : " << mesParametresTab[nbpop-1]->ordreVehicules[1][0].maxRouteTime <<  endl ;
				solver.evolve(2000,1);
				
				// Checking if we need to go to the next fleet or distance constraint value
				if (populationTab[nbpop-1]->getIndividuBestValide () != NULL)
				{
					if (minFleetSize) 
						veh -- ; // reducing the fleet size (PCARP)
					else if (minMaxTour) 
						distConstraint = populationTab[nbpop-1]->getIndividuBestValide()->maxRoute -1 ; // or reducing the distance below the best current solution (MM-kWRPP)
				}
				else 
					validExist = false ;
				
				// in the case of fleet size minimization for the PCARP (type == 32), we can test to see if there is enough capacity left to service all customers (trivial lower bound on fleet size)
				// in this case, don't need to pursue the search further
				if (c.get_type() == 32 && mesParametresTab[nbpop-1]->totalDemand > veh*mesParametresTab[nbpop-1]->ordreVehicules[1][0].vehicleCapacity*mesParametresTab[nbpop-1]->nbDays) 
				{
					cout << "Insufficient capacity -- we can stop decreasing the fleet size" << endl ;
					validExist = false ;
				}

				cout << "  " << endl ;
			}

			// Case of the minimization of the max route length
			// At the end of the process, the search is finished, we return the solution
			if (minMaxTour)
			{
				populationTab[nbpop-2]->ExportBest(c.get_path_to_solution());
				populationTab[nbpop-2]->ExportBKS(c.get_path_to_BKS());
			}
			else
			// Case of the minimization of the fleet size
			// We should not forget the secondary objective in the hierarchy, which is now to minimize the distance for the resulting fleet size
			// Thus, a last optimization run is done
			{
				veh ++ ;
				cout << "######### Second phase : minimizing Distance with " << veh << " vehicles" << endl ;
				mesParametres2 = new Params(c.get_path_to_instance(),c.get_path_to_solution(),c.get_path_to_BKS(),c.get_path_to_distmat(),c.get_seed(),c.get_type(),veh,c.get_nbDep(),false) ;
				population2 = new Population(mesParametres2) ;
				if (nbpop >= 1 && populationTab[nbpop-1]->getIndividuBestValide() != NULL) 
					population2->addAllIndividus(populationTab[nbpop-1]);
				else if (nbpop >= 2) 
					population2->addAllIndividus(populationTab[nbpop-2]);
				Genetic solver(mesParametres2,population2,nb_ticks_allowed,true);
				solver.evolve(10000,1);

				// Returning the final solution
				population2->ExportBest(c.get_path_to_solution());
				population2->ExportBKS(c.get_path_to_BKS());
				
				// Clearing the data structures
				delete population2 ;
				delete mesParametres2 ;
			}

			// Clearing the data structures
			FreeClear (populationTab) ;
			FreeClear (mesParametresTab) ;
			cout << endl ;
			return 0;
		}
	}
	catch(const string& e)
	{
		cout << e << endl ;
		cout << endl ;
		return 0 ;
	}
}
