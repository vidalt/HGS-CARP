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

#include "Genetic.h"

void Genetic::evolve (int maxIterNonProd, int nbRec)
{
	// Should we run either ILS or and HGA 
	if (params->isILS_general) 
		evolveILS();
	else 
		evolveHGA(maxIterNonProd, nbRec);
}

void Genetic::evolveHGA (int maxIterNonProd, int nbRec)
{
	// Main code of the HGA 
	Individu * parent1 ;
	Individu * parent2 ;
	int place, place2 ;
	place2 = 10000 ;
	nbIterNonProd = 1 ;
	nbIter = 0 ;
	int resultCross ;
	string temp ;
	double fitBeforeRepair ;
	CoutSol bestSolFeasibility ;
	clock_t debut = clock() ; // When iterating several time the HGA (e.g. PCARP, the time limit applies to one iteration -- fleet size or max distance value)

	if (population->getIndividuBestValide() != NULL) bestSolFeasibility = population->getIndividuBestValide()->coutSol ;
	else bestSolFeasibility = population->getIndividuBestInvalide()->coutSol ;
	for (int i=0 ; i<population->invalides->nbIndiv ; i++)
		if (population->invalides->individus[i]->coutSol.isBetterFeas(bestSolFeasibility)) bestSolFeasibility = population->invalides->individus[i]->coutSol ;

	rejeton->localSearch->nbTotalRISinceBeginning = 0 ;
	rejeton->localSearch->nbTotalPISinceBeginning = 0 ;

	cout << "| Start of GA | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	while (nbIterNonProd < maxIterNonProd && (clock() - debut <= ticks) && (!params->isSearchingFeasible || population->getIndividuBestValide() == NULL))
	{
		// CROSSOVER
		parent1 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		parent2 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		rejeton->recopieIndividu(rejeton,parent1); // Put them in adequate data structures
		rejeton2->recopieIndividu(rejeton2,parent2); // Put them in adequate data structures

		if (!params->periodique && !params->multiDepot) 
			resultCross = crossOX(); // Pick OX crossover if its a single-period problem
		else 
			resultCross = crossPIX() ; // Otherwise PIX (see Vidal et al 2012 -- OR)

		// SPLIT
		rejeton->generalSplit();

		// LOCAL SEARCH
		rejeton->updateLS();
		rejeton->localSearch->runSearchTotal();
		rejeton->updateIndiv();
		population->updateNbValides(rejeton);
		place = population->addIndividu(rejeton) ;

		// POSSIBLE REPAIR
		if (!rejeton->estValide) 
		{
			fitBeforeRepair = rejeton->coutSol.evaluation ;
			if (rand() % 2 == 0) // 50% chance to do repair on an infeasible individual
			{
				reparer();
				if (rejeton->coutSol.evaluation < fitBeforeRepair - 0.01 || rejeton->coutSol.evaluation > fitBeforeRepair + 0.01 || rejeton->estValide) 
					place2 = population->addIndividu(rejeton) ;
				if (rejeton->estValide)
					place = place2 ;
				else 
					place = min(place,place2);
			}
		}

		// SOME TRACES
		if ( (rejeton->estValide && place == 0) || (rejeton->coutSol.isBetterFeas(bestSolFeasibility) && population->valides->nbIndiv == 0))
		{	
			if (traces && population->valides->nbIndiv > 0) 
				cout << "NEW BEST FEASIBLE " << place << " " << population->getIndividuBestValide()->coutSol.evaluation << " distance : " << rejeton->coutSol.distance << " nbRoutes : " << rejeton->coutSol.routes << " capaViol : " << rejeton->coutSol.capacityViol << " lengthViol : " << rejeton->coutSol.lengthViol << endl << endl ;
			if (traces && population->valides->nbIndiv == 0 ) 
				cout << "NEW BEST INFEASIBLE "<< place << " " << rejeton->coutSol.evaluation                            << " distance : " << rejeton->coutSol.distance << " nbRoutes : " << rejeton->coutSol.routes << " capaViol : " << rejeton->coutSol.capacityViol << " lengthViol : " << rejeton->coutSol.lengthViol << endl << endl ;
			if (rejeton->coutSol.isBetterFeas(bestSolFeasibility)) 
				bestSolFeasibility = rejeton->coutSol ;
			nbIterNonProd = 1 ; 
		}
		else nbIterNonProd ++ ;


		// DIVERSIFICATION
		if (nbRec > 0 && nbIterNonProd % (maxIterNonProd/3+1) == maxIterNonProd/3) 
		{	
			if (traces) cout << "Diversification" << endl ;
			population->diversify();
		}

		// PENALTY MANAGEMENT
		if (nbIter % 30 == 0) 
			gererPenalites () ;

		// MORE TRACES
		if (traces && nbIter % 500 == 0)
		{
			population->afficheEtat(nbIter);
			cout << " | NbTotalMovesLS : (RI) " << rejeton->localSearch->nbTotalRISinceBeginning << " | (PI) " << rejeton->localSearch->nbTotalPISinceBeginning << " | " << endl ;
			cout << " | interSwap " << rejeton->localSearch->nbInterSwap ;
			cout << " | intraSwap " << rejeton->localSearch->nbIntraSwap ;
			cout << " | inter2opt " << rejeton->localSearch->nbInter2Opt ;
			cout << " | intra2opt " << rejeton->localSearch->nbIntra2Opt ;
			cout << " | " << endl ;
			cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
			cout << endl ;
		}
		nbIter ++ ;
	}

	// END OF THE ALGORITHM
	if (traces)
	{
		cout << "Time Elapsed : " << clock() << endl ;
		cout << "Number of Iterations : " << nbIter << endl ;
	}
}

void Genetic::evolveILS ()
{
	int nbGRASP = 5 ;
	int nbILS = 100 ;
	int nbCHILD = 50 ;
	bool isFirstLoop ;
	Individu * parent ;
	clock_t timeBest2 ;
	rejeton->localSearch->nbTotalRISinceBeginning = 0 ;
	rejeton->localSearch->nbTotalPISinceBeginning = 0 ;
	nbIter = 0 ;
	clock_t debut = clock();
	rejetonBestFoundAll->coutSol.evaluation = 1.e30 ;

	cout << "| Debut evolution ILS | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	for (int phaseGrasp = 0 ; phaseGrasp < nbGRASP ; phaseGrasp ++)
	{
		// NEW RANDOM START
		cout << endl  << "------------ RESTART ---------" << endl << endl ;
		params->penalityCapa = 50 ;
		params->penalityLength = 50;
		Individu * myIndiv = new Individu(params, 1.0) ; 
		rejetonP1->recopieIndividu(rejetonP1,myIndiv);
		delete myIndiv ;
		rejetonBestFound->coutSol.evaluation = 1.e30 ;
		isFirstLoop = true ;
		
		for (int nbGenerationNonProd = 0 ; nbGenerationNonProd < nbILS && (clock() - debut <= ticks) ; nbGenerationNonProd ++)
		{
			if (!isFirstLoop)
			{
				parent = population->getIndividuBestValide();
				if (parent == NULL) parent = population->getIndividuBestInvalide();
				rejetonP1->recopieIndividu(rejetonP1,parent);
			}
			else 
				isFirstLoop = false ;

			// Clear the population
			population->clear();

			for (int i=0 ; i < nbCHILD ; i++)
			{
				// SHAKING
				rejeton->recopieIndividu(rejeton,rejetonP1);
				rejeton->shakingSwap(2 + (int)params->nbClients/200);
				rejeton->generalSplit();

				// LOCAL SEARCH
				rejeton->updateLS();
				rejeton->localSearch->runSearchTotal();
				rejeton->updateIndiv();
				population->updateNbValides(rejeton);
				// If the solution is infeasible, do a LS with higher penalty to restore feasibiliy
				if (!rejeton->estValide) reparer();
				population->addIndividu(rejeton) ;

				// Checking if its a solution improvement
				if (rejeton->estValide && rejeton->coutSol.evaluation < rejetonBestFound->coutSol.evaluation - 0.001) 
				{	
					nbGenerationNonProd = -1 ;
					rejetonBestFound->recopieIndividu(rejetonBestFound,rejeton);
					if (rejetonBestFound->coutSol.evaluation < rejetonBestFoundAll->coutSol.evaluation - 0.001)
					{
						cout << "NEW BEST EVER : " << rejetonBestFound->coutSol.evaluation << endl ;
						rejetonBestFoundAll->recopieIndividu(rejetonBestFoundAll,rejetonBestFound);
						timeBest2 = population->timeBest ;
					}
					else
						cout << "NEW BEST      : " << rejeton->coutSol.evaluation << endl ;
				}

				// Regular adaptation of penalty parameters (every 30 LS)
				if (nbIter % 30 == 0) 
					gererPenalites () ;
				nbIter ++ ;
			}

			if (nbIter % 500 == 0)
			{
				population->afficheEtat(nbIter);
				cout << " | NbTotalMovesLS : (RI) " << rejeton->localSearch->nbTotalRISinceBeginning << " | (PI) " << rejeton->localSearch->nbTotalPISinceBeginning << " | " << endl ;
				cout << " | interSwap " << rejeton->localSearch->nbInterSwap ;
				cout << " | intraSwap " << rejeton->localSearch->nbIntraSwap ;
				cout << " | inter2opt " << rejeton->localSearch->nbInter2Opt ;
				cout << " | intra2opt " << rejeton->localSearch->nbIntra2Opt ;
				cout << " | " << endl ;
				cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
				cout << endl ;
			}
		}
	}

	// fin de l'algorithme , diverses informations affichées
	if (traces) 
	{
		cout << "temps passe : " << clock() << endl ;
		cout << "fin evolution ILS, nombre d'iterations : " << nbIter << endl ;
	}

	// ajouter la meilleure solution trouvée dans la population pour l'écriture en fin d'algorithme
	population->addIndividu(rejetonBestFoundAll);
	population->timeBest = timeBest2 ; // need to correct the time to best solution.
}

void Genetic::reparer ()
{
	double temp, temp2  ;

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;

	// First Tentative of Repair
	params->penalityCapa *= 10 ;
	params->penalityLength *= 10 ;
	rejeton->updateLS();
	rejeton->localSearch->runSearchTotal();
	rejeton->updateIndiv();

	// If the first tentative failed, second tentative with higher penalty
	if (!rejeton->estValide) 
	{
		params->penalityCapa *= 10 ;
		params->penalityLength *= 10 ;
		rejeton->generalSplit();
		rejeton->updateLS();
		rejeton->localSearch->runSearchTotal();
		rejeton->updateIndiv();
	}

	params->penalityCapa = temp ;
	params->penalityLength = temp2 ;
	rejeton->measureSol();
}

void Genetic::gererPenalites ()
{
	bool changeDone = false ;
	double fractionCharge = population->fractionValidesCharge() ;
	double fractionTemps = population->fractionValidesTemps() ;

	// if there are not enough feasible solutions
	if ( fractionCharge < params->minValides && params->penalityCapa < 5000)
	{
		params->penalityCapa = (double)((float)params->penalityCapa * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionCharge > params->maxValides && params->penalityCapa > 0.01)
	{
		params->penalityCapa =  (double)((float)params->penalityCapa * 0.85) ;
		changeDone = true ;
	}

	// if there are too many feasible solutions
	if ( fractionTemps < params->minValides && params->penalityLength < 5000)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionTemps > params->maxValides && params->penalityLength > 0.01)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 0.85) ;
		changeDone = true ;
	}

	if (changeDone) 
		population->validatePen(population->invalides);
}

int Genetic::crossOX ()
{
	int temp, tempSuiv ;

	// We pick the beginning and end of the crossover zone
	int debut = rand() % params->nbClients ;
	int fin = rand() % params->nbClients ;
	while (fin == debut && params->nbClients > 1)
		fin = rand() % params->nbClients ;

	// We initialize a little frequency table to know if each customer was placed or not
	for (int i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient[i] = 1 ;

	int j = debut ;
	// we keep the elements from "debut" to "end"
	while ((j % params->nbClients) != ((fin + 1) % params->nbClients))
	{ 
		freqClient[rejeton->chromT[1][j % params->nbClients]] = 0 ;
		j ++ ;
	}

	// We fill the rest of the elements in the order of the second parent
	for (int i=1 ; i <= params->nbClients ; i++)
	{
		temp = rejeton2->chromT[1][(fin + i) % params->nbClients] ;
		tempSuiv = rejeton2->chromT[1][(fin + i + 1) % params->nbClients] ;
		if (freqClient[temp] == 1)
		{
			rejeton->chromT[1][j % params->nbClients] = temp ;
			j ++ ;
		}
	}

	return 0 ;
}

int Genetic::crossPIX ()
{
	vector < int > vide, garder, joursPerturb, tableauFin, tableauEtat ;
	vector < vector <int> > garder2 ;
	int jj,i,ii,j,temp,size ;
	int debut, fin, day ;
	int j1,j2 ;
	int placeInsertion ;

	// We initialize some little data structures
	for (i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		freqClient[i] = params->cli[i].freq ;
		rejeton->chromP[i].pat = 0 ;
		rejeton->chromP[i].dep = -1 ;
	}


	// We take the days in random order
	for (int k=1 ; k <= params->nbDays ; k++ )
		joursPerturb.push_back(k) ;
	for (i = 0 ; i < (int)joursPerturb.size() ; i++)
	{
		jj = i + rand() % ((int)joursPerturb.size() - i) ;
		temp = joursPerturb[i] ;
		joursPerturb[i] = joursPerturb[jj] ;
		joursPerturb[jj] = temp ;
	}

	// We pick j1 and j2
	j1 = rand() % params->nbDays ;
	j2 = rand() % params->nbDays ;
	if (j1 > j2)
	{
		temp = j2 ;
		j2 = j1 ;
		j1 = temp ;
	}

	for (int k = 0 ; k <= params->nbDays ; k ++ ) 
		garder2.push_back(vide);

	// Inheriting the visits
	for (int k = 0 ; k < params->nbDays ; k ++ )
	{
		day = joursPerturb[k];

		// First case, we copy a segment (these visits will be temporarily kept in the data structure "garder2")
		if (k < j1 && !rejeton->chromT[day].empty())
		{
			debut = (int)(rand() % rejeton->chromT[day].size()) ;
			fin = (int)(rand() % rejeton->chromT[day].size()) ;
			tableauFin.push_back(fin);
			j = debut ;
			while ( j != ((fin + 1) % rejeton->chromT[day].size()) )
			{
				freqClient[rejeton->chromT[day][j]] -= 1 ;
				rejeton->chromP[rejeton->chromT[day][j]].pat += (int)pow((float)2,(int)((params->nbDays-day)%params->ancienNbDays)) ;
				rejeton->chromP[rejeton->chromT[day][j]].dep = (day-1)/params->ancienNbDays ;
				garder2[day].push_back(rejeton->chromT[day][j]) ;
				j = (j+1) % rejeton->chromT[day].size() ;
			}
			rejeton->chromT[day].clear();
		}
		else if (k < j2) // Second case, we copy nothing on this day
		{
			rejeton->chromT[day].clear() ;
			tableauFin.push_back(-1);
		}
		else // Third case, we copy everything on this day
		{
			tableauFin.push_back(0);
			for (j=0 ; j < (int)rejeton->chromT[day].size() ; j++)
			{
				freqClient[rejeton->chromT[day][j]] -= 1 ;
				rejeton->chromP[rejeton->chromT[day][j]].pat += (int)pow((float)2,(int)((params->nbDays-day)%params->ancienNbDays)) ;
				rejeton->chromP[rejeton->chromT[day][j]].dep = (day-1)/params->ancienNbDays ;
			}
		}
	}

	// We complete with the second parent
	for (int k = 0 ; k < params->nbDays ; k ++ )
	{
		day = joursPerturb[k] ;
		fin = tableauFin[k] ;
		if (k < j2)
		{
			for (i=0 ; i < (int)rejeton2->chromT[day].size() ; i++)
			{
				ii = rejeton2->chromT[day][ (i + fin + 1) % (int)rejeton2->chromT[day].size() ] ;
				if (freqClient[ii] != 0
					&& params->cli[ii].jourSuiv[rejeton->chromP[ii].pat][(day-1)%params->ancienNbDays+1] == (int)((day-1)%params->ancienNbDays+1) 
					&& ( rejeton->chromP[ii].dep == -1 || rejeton->chromP[ii].dep == (day-1)/params->ancienNbDays ))
				{
					rejeton->chromT[day].push_back(ii);
					freqClient[ii] -= 1 ;
					rejeton->chromP[ii].pat += (int)pow((float)2,(int)((params->nbDays-day)%params->ancienNbDays)) ;
					rejeton->chromP[ii].dep = (day-1)/params->ancienNbDays ;
				}
			}
		}
	}

	// we complete with the elements of garder2 (which come from the first parent for the days where the parents are mixed)
	for (int k=1 ; k<=params->nbDays ; k++)
	{
		garder.clear();
		// Choose a random place of insertion
		size = (int)rejeton->chromT[k].size() ;
		if (size != 0) placeInsertion = rand() % size ; 
		else placeInsertion = 0 ;
		for (int iii=placeInsertion ; iii <  size ; iii ++)
			garder.push_back(rejeton->chromT[k][iii]);
		for (int iii=placeInsertion ; iii <  size ; iii ++)
			rejeton->chromT[k].pop_back();
		for (int iii=0 ; iii < (int)garder2[k].size() ; iii ++)
			rejeton->chromT[k].push_back(garder2[k][iii]);
		for (int iii=0 ; iii < (int)garder.size() ; iii ++)
			rejeton->chromT[k].push_back(garder[iii]);
	}

	rejeton->toPlace.clear();
	// We gather in "toPlace" those elements with missing visits
	for (i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		if (freqClient[i] > 0)
		{
			for (j=0 ; j< freqClient[i] ; j++)
				rejeton->toPlace.push_back(i);
		}
	}

	// We randomize toPlace
	for (i = 0 ; i < (int)rejeton->toPlace.size() ; i++)
	{
		jj = i + rand() % ((int)rejeton->toPlace.size() - i) ;
		temp = rejeton->toPlace[i] ;
		rejeton->toPlace[i] = rejeton->toPlace[jj] ;
		rejeton->toPlace[jj] = temp ;
	}

	// We call Split to obtain a full solution
	rejeton->generalSplit();
	rejeton->updateLS();
	rejeton->localSearch->placeManquants(); // and we perform a least-cost insertion of the missing visits
	rejeton->updateIndiv();

	return 0 ;
}

Genetic::Genetic(Params * params,Population * population, clock_t ticks, bool traces) : 
ticks(ticks), traces(traces), population(population), params(params)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient.push_back(params->cli[i].freq);

	// Creating the Individuals that serve to perform the Local Search and other operations
	rejeton = new Individu (params, true) ; 
	rejeton2 = new Individu(params, true) ; 
	rejetonP1 = new Individu(params, true) ; 
	rejetonP2 = new Individu(params, true) ; 
	rejetonBestFound = new Individu(params, true) ; 
	rejetonBestFoundAll = new Individu(params, true) ; 
	rejeton->localSearch = new LocalSearch(params,rejeton) ;
} 

Genetic::~Genetic(void)
{ 
	delete rejeton ;
	delete rejeton2 ;
	delete rejetonP1 ;
	delete rejetonP2 ;
	delete rejetonBestFound ;
	delete rejetonBestFoundAll ;
}

