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

#include "LocalSearch.h"
#include "Individu.h"

void LocalSearch::runSearchTotal ()
{
	// Shuffling the order of move evaluations
	params->shuffleProches();
	melangeParcours();
	int nbMoves = 0 ;

	/* RUNNING THE LS */
	// RI -- Route improvement
	updateMoves ();
	reinitAllSingleDayMoves();
	for (int day = 1 ; day <= params->nbDays ; day++)
		nbMoves += mutationSameDay (day) ;
	nbTotalRISinceBeginning += nbMoves ;

	// PI and RI local search (only for PCARP and MDCARP, see Vidal et al 2012 (OR)
	if (params->periodique || params->multiDepot)
	{
		// PI -- Pattern improvement
		nbMoves = mutationDifferentDay () ; 
		nbTotalPISinceBeginning += nbMoves ;

		if (nbMoves > 0)
		{
			// RI -- Route improvement
			updateMoves (); 
			reinitAllSingleDayMoves();
			nbMoves = 0 ;
			for (int day = 1 ; day <= params->nbDays ; day++)
				nbMoves += mutationSameDay (day) ;
			nbTotalRISinceBeginning += nbMoves ;
		}
	}
}

int LocalSearch::mutationSameDay (int day)
{
	// Local Search for one given day
	// Based on the classic moves (Relocate, Swap, CROSS, 2-Opt and 2-Opt*)
	rechercheTerminee = false ;
	double costBeforeRemoval, costAfterRemoval;
	bool gainWhenRemoving ;
	int moveEffectue = 0;
	int nbMoves = 0;
	bool routeVideTestee ;
	Noeud * tempNoeud ;
	int size2 ;

	// We search and apply moves until a local minimum is attained
	while (!rechercheTerminee)
	{
		rechercheTerminee = true ;
		moveEffectue = 0 ;

		// For every node U in random order
		for ( int posU = 0 ; posU < (int)ordreParcours[day].size() ; posU ++ )
		{
			posU -= moveEffectue ; // We return on the last node if there has been a move
			nbMoves += moveEffectue ;
			moveEffectue = 0 ;
			routeVideTestee = false ;

			noeudU = &clients[day][ordreParcours[day][posU]];
			routeU = noeudU->route ;
			vehicleU = routeU->vehicle ;
			x = noeudU->suiv ;

			// In the CARP, some service removal do not reduce the cost of the route
			// In this case, very few moves can improve the solution (only 2-opt variants), and thus SWAP and RELOCATE variants do not need to be tested
			costBeforeRemoval = noeudU->seq0_i->evaluation(routeU->depot->pred->seq0_i,routeU->vehicle);
			costAfterRemoval = noeudU->seq0_i->evaluation(noeudU->pred->seq0_i,noeudU->suiv->seqi_n,routeU->vehicle);
			gainWhenRemoving = (costAfterRemoval < costBeforeRemoval - 0.1) ;

			// For every node V in random order
			size2 = (int)noeudU->moves.size();
			for ( int posV = 0 ; posV < size2 && moveEffectue == 0 ; posV ++ )
			{
				noeudV = &clients[day][noeudU->moves[posV]] ;
				routeV = noeudV->route ;
				vehicleV = routeV->vehicle ;

				// If we have not yet tested the moves involving the node U and the route of node V
				// (This flag is reset to false as soon as there is a modification in the route)
				if (!noeudV->route->nodeAndRouteTested[noeudU->cour])
				{
					y = noeudV->suiv ;
					if (routeV->cour != routeU->cour)
					{
						if (moveEffectue != 1)
						{
							tempNoeud = noeudV ;
							noeudV = noeudV->suiv ;
							y = noeudV->suiv ;
							// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
							// Case where they are in different routes
							if (gainWhenRemoving) moveEffectue = interRouteGeneralInsert(); 
							noeudV = tempNoeud ;
							y = noeudV->suiv ;
						}

						// 2-Opt*
						if (moveEffectue != 1) 
							moveEffectue = interRoute2Opt ();

						// 2-Opt* (second type, where the routes can be reversed)
						if (moveEffectue != 1) 
							moveEffectue = interRoute2OptInv ();
					}
					else
					{
						tempNoeud = noeudV ;
						noeudV = noeudV->suiv ;
						y = noeudV->suiv ;

						// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
						// Case where they are in the same route
						if (moveEffectue != 1 && gainWhenRemoving) 
							moveEffectue = intraRouteGeneralInsertDroite();
						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
				}
			}

			noeudV = noeudU->suiv ;
			routeV = noeudV->route ;
			vehicleV = routeV->vehicle ;
			y = noeudV->suiv ;
			if (!noeudV->route->nodeAndRouteTested[noeudU->cour])
			{
				while (moveEffectue != 1 && !noeudV->estUnDepot)
				{ 
					// Testing 2-Opt between U and V (if the restriction of the granular search allows) 
					if (params->isCorrelated[noeudU->pred->cour][noeudV->cour] || params->isCorrelated[noeudU->cour][noeudV->suiv->cour]) 
						moveEffectue = intraRoute2Opt ();
					noeudV = noeudV->suiv ;
					y = noeudV->suiv ;
				}
			}

			// Special cases : testing the insertions behind the depot, and the empty routes
			for (int route = 0 ; route < params->nbVehiculesPerDep && moveEffectue == 0 ; route ++)
			{
				noeudV = &depots[day][route] ;
				routeV = noeudV->route ;
				y = noeudV->suiv ;
				if ( (!noeudV->route->nodeAndRouteTested[noeudU->cour]) && (!y->estUnDepot || !routeVideTestee))
				{
					if (y->estUnDepot) routeVideTestee = true ;
					if (routeV != routeU)
					{
						tempNoeud = noeudV ;
						noeudV = depots[day][route].suiv ;
						y = noeudV->suiv ;

						// Insertion after the depot, in a different route
						if (gainWhenRemoving && (params->isCorrelated[noeudU->cour][noeudV->cour] || params->isCorrelated[noeudU->cour][y->cour]) && moveEffectue != 1 ) 
							moveEffectue = interRouteGeneralInsert();

						noeudV = noeudV->route->depot ;
						y = noeudV->suiv ;

						// 2-Opt* after the depot
						if (params->isCorrelated[noeudU->pred->cour][noeudV->cour] && moveEffectue != 1) 
							moveEffectue = interRoute2Opt ();

						// 2-Opt* after the depot
						if ((params->isCorrelated[x->cour][y->cour] || params->isCorrelated[y->cour][x->cour]) && moveEffectue != 1) 
							moveEffectue = interRoute2OptInv ();

						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
					else
					{
						tempNoeud = noeudV ;
						noeudV = depots[day][route].suiv ;
						y = noeudV->suiv ;

						// Insertion after the depot, in the same route
						if ((params->isCorrelated[noeudU->cour][noeudV->cour] || params->isCorrelated[noeudU->cour][y->cour]) && moveEffectue != 1) 
							moveEffectue = intraRouteGeneralInsertDroite();

						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
				}
			}

			// Say that we have tested the node U with all routes
			if (moveEffectue == 0)
				nodeTestedForEachRoute(noeudU->cour,day);
		}
	}
	// Calling the ejection chains at the end of the LS
	nbMoves += ejectionChains(day);
	return nbMoves ;
}

int LocalSearch::mutationDifferentDay ()
{
	// Local Search to improve the pattern choices for customers (PCARP and MDCARP)
	// Only a single move, which is a relocate of all occurences of a customer in the best combination of days
	// This move is still used for problems with a single period, in this case it only does a relocate without granular search restriction
	rechercheTerminee = false ;
	int nbMoves = 0 ;
	firstLoop = true ;

	for (int day = 1 ; day <= params->nbDays ; day++)
		for (int r=0 ; r < params->nombreVehicules[day] ; r++)
			routes[day][r].initiateInsertions() ;

	while ( !rechercheTerminee )
	{
		rechercheTerminee = true ;
		// Searching for a better insertion place for all customers
		for ( int posU = 0 ; posU < params->nbClients ; posU ++ )
			nbMoves += searchBetterPattern(ordreParcours[0][posU]);
		firstLoop = false ;
	}
	return nbMoves ;
}

int LocalSearch::interRouteGeneralInsert()
{
	// For a pair of nodes U, V, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	// Some route evaluations can be gained (about 40%) by doing these computations is a combined manner
	int ibest, jbest ;
	double moveMin ;
	double temp ;
	SeqData * seq = noeudU->seq0_i ;
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++)
			resultMoves[i][j] = 1.e29 ; 

	// This keeps the current cost
	resultMoves[0][0] = costZero ;

	/* EVALUATION OF MOVE LOWER BOUNDS */
	// We start to compute lower bounds on move evaluations, so that we can discard any moves which are not promising.

	// Can send something if U is not a depot
	if (!noeudU->estUnDepot)
	{
		// We evaluate the result of the move as a combination of existing subsequences
		resultMoves[1][0] = seq->evaluationLB(noeudU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		// Can we receive something (if V is not a depot)
		// In the following, all the conditionals are set up to avoid moving depots
		if (!noeudV->estUnDepot)
			resultMoves[1][1] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq1,x->seqi_n,routeU->vehicle);

		// And so on...
		if (!x->estUnDepot)
		{
			resultMoves[2][0] = seq->evaluationLB(noeudU->pred->seq0_i,x->suiv->seqi_n,routeU->vehicle);
			resultMoves[3][0] = resultMoves[2][0] ;
			if (!noeudV->estUnDepot)
			{
				resultMoves[2][1] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq1,x->suiv->seqi_n,routeU->vehicle);
				resultMoves[3][1] = resultMoves[2][1] ;
				if (!y->estUnDepot)
				{
					resultMoves[2][2] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq12,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[2][3] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq21,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[3][2] = resultMoves[2][2];
					resultMoves[3][3] = resultMoves[2][3];
				}
			}
		}
	}

	if (!noeudU->estUnDepot)
	{
		resultMoves[1][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq1,noeudV->seqi_n,routeV->vehicle);
		if (!x->estUnDepot)
		{
			resultMoves[2][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,noeudV->seqi_n,routeV->vehicle);
			resultMoves[3][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeV->vehicle);
		}
	}

	if (!noeudV->estUnDepot)
	{
		if (!noeudU->estUnDepot)
		{
			resultMoves[1][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->estUnDepot)
			{
				resultMoves[2][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,y->seqi_n,routeV->vehicle);
				resultMoves[3][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->estUnDepot && !noeudU->estUnDepot && !x->estUnDepot)
		{
			temp = seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,y->suiv->seqi_n,routeV->vehicle);
			resultMoves[2][2] += temp ;
			resultMoves[2][3] += temp ;
			temp = seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,y->suiv->seqi_n,routeV->vehicle);
			resultMoves[3][2] += temp ;
			resultMoves[3][3] += temp ;
		}
	}

	/* WE IDENTIFY WHICH MOVES CAN BE FILTERED OUT USING THE LOWER BOUND */

	for (int i=0 ; i<4; i++)
	{
		for (int j=0 ; j < 4 ; j++)
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true ;
	resultMoves[0][0] = costZero ;

	/* AND NOW WE TEST THE MOVES THAT HAVE A CHANCE TO BE IMPROVING */
	// Exactly the same code as previously, but using "seq->evaluation" instead of "seq->evaluationLB"

	if (!noeudU->estUnDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(noeudU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		if (!noeudV->estUnDepot && shouldBeTested[1][1])
			resultMoves[1][1] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq1,x->seqi_n,routeU->vehicle);

		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0] || shouldBeTested[3][0]) 
			{
				resultMoves[2][0] = seq->evaluation(noeudU->pred->seq0_i,x->suiv->seqi_n,routeU->vehicle);
				resultMoves[3][0] = resultMoves[2][0] ;
			}
			if (!noeudV->estUnDepot)
			{
				if (shouldBeTested[2][1] || shouldBeTested[3][1])
				{
					resultMoves[2][1] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq1,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[3][1] = resultMoves[2][1] ;
				}
				if (!y->estUnDepot)
				{
					if (shouldBeTested[2][2] || shouldBeTested[3][2])
					{
						resultMoves[2][2] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq12,x->suiv->seqi_n,routeU->vehicle);
						resultMoves[3][2] = resultMoves[2][2];
					}
					if (shouldBeTested[2][3] || shouldBeTested[3][3])
					{
						resultMoves[2][3] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq21,x->suiv->seqi_n,routeU->vehicle);
						resultMoves[3][3] = resultMoves[2][3];
					}
				}
			}
		}
	}

	if (!noeudU->estUnDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq1,noeudV->seqi_n,routeV->vehicle);
		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0]) resultMoves[2][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,noeudV->seqi_n,routeV->vehicle);
			if (shouldBeTested[3][0]) resultMoves[3][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeV->vehicle);
		}
	}

	if (!noeudV->estUnDepot)
	{
		if (!noeudU->estUnDepot)
		{
			if (shouldBeTested[1][1]) resultMoves[1][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->estUnDepot)
			{
				if (shouldBeTested[2][1]) resultMoves[2][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,y->seqi_n,routeV->vehicle);
				if (shouldBeTested[3][1]) resultMoves[3][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->estUnDepot && !noeudU->estUnDepot && !x->estUnDepot)
		{
			if (shouldBeTested[2][2] || shouldBeTested[2][3])
			{
				temp = seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,y->suiv->seqi_n,routeV->vehicle);
				resultMoves[2][2] += temp ;
				resultMoves[2][3] += temp ;
			}
			if (shouldBeTested[3][2] || shouldBeTested[3][3])
			{
				temp = seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,y->suiv->seqi_n,routeV->vehicle);
				resultMoves[3][2] += temp ;
				resultMoves[3][3] += temp ;
			}
		}
	}

	// We identify the best move among relocate, swap, CROSS and I-CROSS between the node pair U,V
	ibest = 0 ; jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4; i++ )
	{
		for (int j=0 ; j < 4 ; j++)
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// If no improving move between U and V, we return
	if ( ibest == 0 && jbest == 0) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	Noeud * placeV = noeudV->pred ;
	Noeud * placeU = noeudU->pred ;
	reinitSingleDayMoves(placeU->route); // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(placeV->route); // Say that all moves involving this route must be tested again

	// Moving the nodes in the solution structure
	if ( ibest == 1 || ibest == 2) insertNoeud(noeudU,placeV);
	if ( ibest == 2) insertNoeud(x,noeudU);
	if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(noeudU,x); }

	if ( jbest == 1 || jbest == 2) insertNoeud(noeudV,placeU);
	if ( jbest == 2) insertNoeud(y,noeudV);
	if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(noeudV,y); }

	// Update the pre-processed data on the subsequences of the route
	placeU->route->updateRouteData(false);
	placeV->route->updateRouteData(false);
	setRouteVide(noeudU->jour); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInterSwap ++ ;
	return 1 ; // Return Success
}

int LocalSearch::interRoute2Opt ()
{
	// Testing 2-Opt* between U and V
	double cost ;
	SeqData * seq = noeudU->seq0_i ;

	if  (routeU->depot->cour != routeV->depot->cour || (noeudU->pred->estUnDepot && noeudV->estUnDepot)) 
		return 0 ; // Cannot do a 2-Opt* if the depot is placed in the wrong way

	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ; // Cost before the move

	// Lower bound on the move value
	cost = seq->evaluationLB(noeudU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluationLB(noeudV->seq0_i,noeudU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; // Exit if no chance of improvement

	// Exact move evaluation
	cost = seq->evaluation(noeudU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluation(noeudV->seq0_i,noeudU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; 

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitSingleDayMoves(routeU);  // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(routeV);  // Say that all moves involving this route must be tested again
	Noeud * tempU = noeudU ;
	noeudU = noeudU->pred ;
	x = noeudU->suiv ;

	// Updating the solution
	Noeud * count ;
	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = depotU->pred ;
	Noeud * depotVFin = depotV->pred ;
	Noeud * depotUpred = depotUFin->pred ;

	count = y ;
	while ( !count->estUnDepot )
	{
		count->route = routeU ;
		count = count->suiv ;
	}

	count = x ;
	while ( !count->estUnDepot )
	{
		count->route = routeV ;
		count = count->suiv ;
	}

	noeudU->suiv = y ;
	y->pred = noeudU ;
	noeudV->suiv = x ;
	x->pred = noeudV ;

	if (x->estUnDepot)
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		noeudV->suiv = depotVFin ;
		depotVFin->pred = noeudV ;
	}
	else
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		depotVFin->pred = depotUpred ;
		depotVFin->pred->suiv = depotVFin ;
	}

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
	setRouteVide(noeudU->jour); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInter2Opt ++ ;
	noeudU = tempU ;
	x = noeudU->suiv ;
	return 1 ; // Return Success
}


int LocalSearch::interRoute2OptInv()
{
	// 2-Opt* with route inversions
	SeqData * seq = noeudU->seq0_i ;
	double cost ;
	double costTemp ;
	double costTempReverse = 1.e30 ;
	bool reverseRouteU, reverseRouteV ; 
	if (noeudU->route == noeudV->route) { return 0 ; }
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// Compute a lower bound on the value of the move
	cost = - costZero ;
	costTemp = seq->evaluationLB(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluationLB(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluationLB(noeudV->seq0_i,noeudU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluationLB(noeudU->seq0_i,noeudV->seqi_0,routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;  // Exit if no chance of improvement

	// Test the real move cost
	cost = - costZero ;
	costTemp = seq->evaluation(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluation(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluation(noeudV->seq0_i,noeudU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluation(noeudU->seq0_i,noeudV->seqi_0,routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitSingleDayMoves(noeudU->route);  // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(noeudV->route);  // Say that all moves involving this route must be tested again
	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = routeU->depot->pred ;
	Noeud * depotVFin = routeV->depot->pred ;
	Noeud * depotUSuiv = depotU->suiv ;

	// Update the solution
	Noeud * temp ;
	Noeud * yy = y ;
	Noeud * uu = noeudU ;

	while ( !yy->estUnDepot )
	{
		temp = yy->suiv ;
		yy->suiv = yy->pred ;
		yy->pred = temp ;
		yy->route = routeU ;
		yy = temp ;
	}

	while ( !uu->estUnDepot )
	{
		temp = uu->pred ;
		uu->pred = uu->suiv ;
		uu->suiv = temp ;
		uu->route = routeV ;
		uu = temp ;
	}

	noeudV->suiv = noeudU ;
	noeudU->pred = noeudV ;
	y->suiv = x ;
	x->pred = y ;

	if (y->estUnDepot)
	{
		depotVFin->suiv = depotV ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->suiv = depotVFin ;
		depotU->suiv = x ;
		x->pred = depotU ;
	}
	else if (noeudU->estUnDepot )
	{
		depotU->suiv = depotVFin->pred ;
		depotU->suiv->pred = depotU ;
		depotU->pred = depotUFin ;
		depotVFin->pred = noeudV ;
		noeudV->suiv = depotVFin ;
	}
	else
	{
		depotU->suiv = depotVFin->pred ;
		depotU->suiv->pred = depotU ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->suiv = depotVFin ;
	}

	// Reverse if needed
	if (reverseRouteU) routeU->reverse();
	if (reverseRouteV) routeV->reverse();

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
	setRouteVide(noeudU->jour); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInter2Opt ++ ;
	return 1 ; // Return Success
}

int LocalSearch::intraRouteGeneralInsertDroite ()
{
	// For a pair of nodes U, V, IN THE SAME ROUTE, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	Noeud * tempU = noeudU ;
	Noeud * tempV = noeudV ;
	bool turned = false ;

	int decalage = noeudV->place - noeudU->place ; // Decalage is the difference of index between U and V
	if (decalage >= -1 && decalage <= 1) return 0 ; // This means that U and V are already consecutive, testing these moves is useless

	// If decalage < 0, this means that V is on the left of U in the route
	// We don't want to deal with this case, so we simply inverse the notations of U and V 
	// And we deal in the following with only the case where V is on the right
	if ( decalage < 0 )
	{
		noeudU = tempV ;
		noeudV = tempU ;
		x = noeudU->suiv ;
		y = noeudV->suiv ;
		decalage = -decalage ;
		turned = true ;
	}

	double moveMin ;
	int ibest, jbest ;
	SeqData * seq = noeudU->seq0_i ;
	Noeud * Upred = noeudU->pred ;
	Noeud * placeU = noeudU->pred ;
	Noeud * placeV = noeudV->pred ;
	Noeud * xsuiv = x->suiv ;
	Noeud * ysuiv = y->suiv ;
	double costZero = routeU->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++ )
			resultMoves[i][j] = 1.e29 ;

	resultMoves[0][0] = costZero ;

	if (decalage >= 3) // General case
	{
		if (!noeudV->estUnDepot && turned) // relocate of V before U (V cannot be a depot, but U can be)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(noeudV->seq1);
			// The intra-route moves are the only moves which can involve general subsequences (i,j) in the middle of the route
			// To reduce a bit the pre-processing, we keep only sequences of limited size (<= 15), even if this involves to concatenate more pieces
			// The concatenation of the good number of pieces is done by "addSeqDataInPieces"
			addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
			myseqs.push_back(y->seqi_n);
			resultMoves[0][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->estUnDepot) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq12);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[0][2] = seq->evaluationLB(myseqs,routeU->vehicle);
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq21);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[0][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!turned)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(noeudV->seqi_n);
			resultMoves[1][0] = seq->evaluationLB(myseqs,routeU->vehicle);
		}
		if (!noeudV->estUnDepot) // exchange U and V
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(noeudV->seq1);
			addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(y->seqi_n);
			resultMoves[1][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->estUnDepot && turned) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq12);
				addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[1][2] = seq->evaluationLB(myseqs,routeU->vehicle);

				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq21);
				addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[1][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!x->estUnDepot)
		{
			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
				myseqs.push_back(noeudU->seq12);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[2][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}

			if (!noeudV->estUnDepot) // exchange U,X and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange U,X and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[2][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[2][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}

			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
				myseqs.push_back(noeudU->seq21);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[3][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
			if (!noeudV->estUnDepot) // exchange X,U and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange X,U and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[3][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[3][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}
		}
	}
	else if (decalage == 2) // U and V are almost consecutive, we are in this situation : UXVY. Useful moves must be tested as special cases
	{
		// XUVY
		resultMoves[1][0] = seq->evaluationLB(Upred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeU->vehicle);
		if (!noeudV->estUnDepot)
		{
			// VXUY
			resultMoves[1][1] = seq->evaluationLB(Upred->seq0_i,noeudV->seq1,noeudU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			resultMoves[2][1] = seq->evaluationLB(Upred->seq0_i,noeudV->seq1,noeudU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->estUnDepot)
			{
				// VYXU
				resultMoves[1][2] = seq->evaluationLB(Upred->seq0_i,noeudV->seq12,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// YVXU
				resultMoves[1][3] = seq->evaluationLB(Upred->seq0_i,noeudV->seq21,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// VYUX
				resultMoves[2][2] = seq->evaluationLB(Upred->seq0_i,noeudV->seq12,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
				// YVUX
				resultMoves[2][3] = seq->evaluationLB(Upred->seq0_i,noeudV->seq21,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// HERE STARTS THE SECOND PHASE AFTER PREPROCESSING
	// ONLY NON-FILTERED MOVES THAT HAVE A CHANCE TO BE IMPROVING
	/////////////////////////////////////////////////////////////

	for (int i=0 ; i<4 ; i++ )
	{
		for (int j=0 ; j<4 ; j++ )
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true  ;
	resultMoves[0][0] = costZero ;


	// Same procedure as previously, but with "seq->evaluation" instead of "seq->evaluationLB"
	if (decalage >= 3)
	{
		if (!noeudV->estUnDepot && turned)
		{

			if (shouldBeTested[0][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq1);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
				myseqs.push_back(y->seqi_n);
				resultMoves[0][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->estUnDepot) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[0][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[0][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[0][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place,noeudU->jour);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[0][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!turned && shouldBeTested[1][0])
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(noeudV->seqi_n);
			resultMoves[1][0] = seq->evaluation(myseqs,routeU->vehicle);
		}
		if (!noeudV->estUnDepot) // exchange U and V
		{
			if (shouldBeTested[1][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq1);
				addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(y->seqi_n);
				resultMoves[1][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->estUnDepot && turned) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[1][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
					myseqs.push_back(noeudU->seq1);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[1][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[1][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(x,noeudV->place-1-x->place,noeudU->jour);
					myseqs.push_back(noeudU->seq1);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[1][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
				myseqs.push_back(noeudU->seq12);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[2][0] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!noeudV->estUnDepot) // exchange U,X and V
			{
				if (shouldBeTested[2][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange U,X and (V,Y or Y,V)
				{
					if (shouldBeTested[2][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq12);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
						myseqs.push_back(noeudU->seq12);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[2][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[2][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq21);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
						myseqs.push_back(noeudU->seq12);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[2][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}

			if (shouldBeTested[3][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
				myseqs.push_back(noeudU->seq21);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[3][0] = seq->evaluation(myseqs,routeU->vehicle);
			}
			if (!noeudV->estUnDepot) // exchange X,U and V
			{
				if (shouldBeTested[3][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange X,U and (V,Y or Y,V)
				{
					if (shouldBeTested[3][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq12);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
						myseqs.push_back(noeudU->seq21);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[3][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[3][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq21);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place,noeudU->jour);
						myseqs.push_back(noeudU->seq21);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[3][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}
		}
	}
	else if (decalage == 2)
	{
		// XUVY
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(Upred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeU->vehicle);
		if (!noeudV->estUnDepot)
		{
			// VXUY
			if (shouldBeTested[1][1]) resultMoves[1][1] = seq->evaluation(Upred->seq0_i,noeudV->seq1,noeudU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			if (shouldBeTested[2][1]) resultMoves[2][1] = seq->evaluation(Upred->seq0_i,noeudV->seq1,noeudU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->estUnDepot)
			{
				// VYXU
				if (shouldBeTested[1][2]) resultMoves[1][2] = seq->evaluation(Upred->seq0_i,noeudV->seq12,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// YVXU
				if (shouldBeTested[1][3]) resultMoves[1][3] = seq->evaluation(Upred->seq0_i,noeudV->seq21,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// VYUX
				if (shouldBeTested[2][2]) resultMoves[2][2] = seq->evaluation(Upred->seq0_i,noeudV->seq12,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
				// YVUX
				if (shouldBeTested[2][3]) resultMoves[2][3] = seq->evaluation(Upred->seq0_i,noeudV->seq21,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
			}
		}
	}

	// Identify the best move involving U,V and apply
	ibest = 0 ;
	jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4 ; i++)
	{
		for (int j=0 ; j<4 ; j++ )
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// No improving move has been found
	if ( ibest == 0 && jbest == 0 ) 
	{ 
		noeudU = tempU ;
		noeudV = tempV ;
		x = noeudU->suiv ;
		y = noeudV->suiv ;
		return 0 ;
	}


	///////////////////////////////////////////
	// HERE STARTS THE MOVE APPLICATION PROCESS
	///////////////////////////////////////////

	reinitSingleDayMoves(routeU);  // Say that all moves involving this route must be tested again

	// Update the solution
	if (decalage >= 3)
	{
		if ( ibest == 1 || ibest == 2) insertNoeud(noeudU,placeV);
		if ( ibest == 2 ) insertNoeud(x,noeudU);
		if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(noeudU,x); }

		if ( jbest == 1 || jbest == 2) insertNoeud(noeudV,placeU);
		if ( jbest == 2) insertNoeud(y,noeudV);
		if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(noeudV,y); }
	}

	// Special cases of decalage == 2
	else if (decalage == 2)
	{
		if (ibest == 1 && jbest == 0) { insertNoeud(noeudU,x); }
		else if (ibest == 1 && jbest == 1) { insertNoeud(x,noeudV); insertNoeud(noeudU,x);  }
		else if (ibest == 1 && jbest == 2) { insertNoeud(x,y); insertNoeud(noeudU,x);  }
		else if (ibest == 1 && jbest == 3) { insertNoeud(noeudV,y); insertNoeud(x,noeudV); insertNoeud(noeudU,x);  }
		else if (ibest == 2 && jbest == 1) { insertNoeud(noeudV,noeudU->pred) ;  }
		else if (ibest == 2 && jbest == 2) { insertNoeud(noeudV,noeudU->pred) ;  insertNoeud(y,noeudV) ; }
		else if (ibest == 2 && jbest == 3) { insertNoeud(y,noeudU->pred) ;  insertNoeud(noeudV,y) ; }
		else throw string ("ERROR move intra-route") ;
	} 

	routeU->updateRouteData(false); // Update the pre-processed data on the subsequences of the route
	setRouteVide(noeudU->jour); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntraSwap ++ ;
	noeudU = tempU ;
	noeudV = tempV ;
	x = noeudU->suiv ;
	y = noeudV->suiv ;
	return 1 ; // Return Success
}


int LocalSearch::intraRoute2Opt ()
{
	// Evaluation procedure for 2-Opt
	Noeud * nodeNum = noeudU->suiv ;
	Noeud * nodeUpred = noeudU->pred ;
	Noeud * temp ;
	SeqData * seq = noeudU->seq0_i ;

	double cost ;
	double costZero = routeU->currentRouteCost ;

	myseqs.clear();
	myseqs.push_back(noeudU->pred->seq0_i);
	addReverseSeqDataInPieces(noeudU,noeudV->place-noeudU->place,noeudU->jour);
	myseqs.push_back(noeudV->suiv->seqi_n);

	// Compute the lower bound on move value and exits if no possible improvement
	cost = seq->evaluationLB(myseqs,routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Compute the real move value
	cost = seq->evaluation(myseqs,routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Apply the move and update the solution
	reinitSingleDayMoves(routeU);
	noeudU->pred = nodeNum ;
	noeudU->suiv = y ;

	while ( nodeNum != noeudV )
	{
		temp = nodeNum->suiv ;
		nodeNum->suiv = nodeNum->pred ;
		nodeNum->pred = temp ;
		nodeNum = temp ;
	}

	noeudV->suiv = noeudV->pred ;
	noeudV->pred = nodeUpred ;
	nodeUpred->suiv = noeudV ;
	y->pred = noeudU ;
	routeU->updateRouteData(false);  // Update the pre-processed data on the subsequences of the route
	setRouteVide(noeudU->jour); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntra2Opt ++ ;
	return 1 ; // Return Success
}

int LocalSearch::searchBetterPattern (int client)
{
	pattern pattern1 = individu->chromP[client] ;
	pattern pattern2, meilleurPattern ;
	int indexMeilleur = -1 ;
	meilleurPattern.pat = -1000000 ;
	int temp, calcul, depot ;
	double depense = 0 ;
	double meilleureDepense = 1.e30 ;
	Noeud * noeudTravail ;
	deplacementIntraJour = false ; 	// this flag is only raised in case there is a better insertion in the same day

	for (int pat = 0 ; pat < (int)params->cli[client].visits.size() ; pat ++)
	{
		// testing a new pattern
		pattern2 = params->cli[client].visits[pat] ;
		testingIncumbentPattern = (individu->chromP[client].pat == pattern2.pat && individu->chromP[client].dep == pattern2.dep) ;
		calcul = pattern2.pat ;
		depot = pattern2.dep ;
		depense = pattern2.cost ;
		for (int k = 0 ; k < params->ancienNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
			{
				noeudTravail = &clients[params->ancienNbDays-k + depot*params->ancienNbDays][client] ;
				// If the insertion in this day has not been computed yet 
				// (Note that it would be possible to do way faster by considering the fact that only the demand of 
				// a single delivery changes... and sometimes its even the same for different patterns)
				// Still, this is a simplistic implementation for the tests PCARP, no need for the highest performance
				if (noeudTravail->coutInsertion[pat] > 1.e29 || firstLoop)
					computeCoutInsertion(noeudTravail,pat) ;
				depense += noeudTravail->coutInsertion[pat] ;
			}
		}

		if (depense < meilleureDepense - EPSILON_LS
			|| (depense < meilleureDepense + EPSILON_LS && pattern2.dep == pattern1.dep && pattern2.pat == pattern1.pat) )
		{
			meilleureDepense = depense ;
			meilleurPattern = pattern2 ;
			indexMeilleur = pat ;
		}
	}

	if (meilleurPattern.pat == -1000000) 
		throw string ("ERROR when computing the best pattern !") ;

	// Applying the move if a better pattern has been found
	if ( meilleurPattern.pat != pattern1.pat || meilleurPattern.dep != pattern1.dep || deplacementIntraJour)
	{
		// removing the current occurences of this customer
		calcul = pattern1.pat ;
		for (int k = 0 ; k < params->ancienNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
				removeNoeud(&clients[params->ancienNbDays-k+ pattern1.dep*params->ancienNbDays][client]);
		}

		// (PCARP) Updating the chromP (necessary to do now, not later, otherwise the wrong data is set to pre-process the SeqData)
		// When adding the nodes in the next block of instructions
		individu->chromP[client] = meilleurPattern ;

		// Inserting in the new locations
		calcul = meilleurPattern.pat ;
		depot = meilleurPattern.dep ;
		for (int k = 0 ; k < params->ancienNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
			{
				Noeud * hereCli = &clients[params->ancienNbDays-k+ meilleurPattern.dep*params->ancienNbDays][client] ;
				Noeud * thereCli = hereCli->placeInsertion[indexMeilleur] ;
				addNoeud(hereCli,thereCli);
			}	
		}
		//cout << "Inserting Node " << client << " With pattern " << meilleurPattern.pat << " Flag is : " << deplacementIntraJour << endl ;

		rechercheTerminee = false ;
		return 1 ;
	}
	else return 0 ;
}

void LocalSearch::computeCoutInsertion(Noeud * client, int pattern) 
{
	// Computing the best cost for inserting a client in its day
	Route * myRoute ;
	client->coutInsertion[pattern] = 1.e30 ;
	client->placeInsertion[pattern] = NULL ;

	noeudU = client ;
	x = noeudU->suiv ;
	routeU = noeudU->route ;

	// Find the best insertion for each route
	for (int r=0 ; r < params->nombreVehicules[client->jour] ; r++)
	{
		myRoute = &routes[client->jour][r] ;
		if ( myRoute->coutInsertionClient[client->cour][pattern] > 1.e29 || firstLoop ) 
			evalInsertClient(myRoute,client,pattern) ;

		if ( myRoute->coutInsertionClient[client->cour][pattern] < client->coutInsertion[pattern] - EPSILON_LS)
		{
			client->coutInsertion[pattern] = myRoute->coutInsertionClient[client->cour][pattern] ;
			client->placeInsertion[pattern] = myRoute->placeInsertionClient[client->cour][pattern] ;
		}
	}

	// If its possible to improve the placement of a customer in a day where its already placed, and according to its current pattern, then we raise this flag
	if (client->estPresent // its placed here
		&& testingIncumbentPattern // with the same pattern
		&& client->route->coutInsertionClient[client->cour][pattern] > client->coutInsertion[pattern] + EPSILON_LS) // but we can do better
		deplacementIntraJour = true ;
}

void LocalSearch::evalInsertClient (Route * R,Noeud * U,int pattern) 
{
	// Computing the least cost for inserting a client U in a route R
	SeqData * seq = U->seq0_i ;
	Noeud * courNoeud ;
	double leastCost = 1.e30 ;
	double cost ;
	bool firstLoopDep = true ;

	// Tweaking the code to work with the PCARP
	// We need to account for the fact that the demand may change based on the pattern choice
	// Here doing something a bit ugly, which is to store and modify from outside the pre-processed Seqdata "U->seq1"
	// to include the good delivery quantity
	double tempDemand = U->seq1->load ;
	U->seq1->load = params->cli[U->cour].demandPatDay[params->cli[U->cour].visits[pattern].pat][U->jour];

	// Some memory structures to avoid recomputing these things again and again
	R->coutInsertionClient[U->cour][pattern] = 1.e30 ;
	R->placeInsertionClient[U->cour][pattern] = NULL ;

	if (!(U->route == R) || !U->estPresent)
	{
		// Case 1 : U is not in the route
		courNoeud = R->depot->suiv ;
		while (!courNoeud->estUnDepot || firstLoopDep)
		{
			if (courNoeud->estUnDepot) firstLoopDep = false ;
			cost = seq->evaluation(courNoeud->pred->seq0_i,U->seq1,courNoeud->seqi_n,R->vehicle);
			if ( cost < leastCost )
			{
				leastCost = cost ;
				R->placeInsertionClient[U->cour][pattern] = courNoeud->pred ;
			}
			courNoeud = courNoeud->suiv ;
		}
		R->coutInsertionClient[U->cour][pattern] = leastCost - seq->evaluation(R->depot->pred->seq0_i,R->vehicle);
	}
	else
	{
		// Case 2 : U is already in the route R
		leastCost = seq->evaluation(U->pred->seq0_i,U->seq1,U->suiv->seqi_n,R->vehicle);
		R->placeInsertionClient[U->cour][pattern] = U->pred ;
		courNoeud = R->depot->suiv ;
		while (!courNoeud->estUnDepot || firstLoopDep)
		{
			if (courNoeud->estUnDepot) firstLoopDep = false ;

			if (courNoeud->place < U->place)
			{
				myseqs.clear();
				myseqs.push_back(courNoeud->pred->seq0_i);
				myseqs.push_back(U->seq1);
				addSeqDataInPieces(courNoeud,U->place-1-courNoeud->place,courNoeud->jour);
				myseqs.push_back(U->suiv->seqi_n);
				cost = seq->evaluation(myseqs,R->vehicle);
			}
			else if (courNoeud->place > U->place + 1)
			{
				myseqs.clear();
				myseqs.push_back(U->pred->seq0_i);
				addSeqDataInPieces(U->suiv,courNoeud->place-1-U->suiv->place,courNoeud->jour);
				myseqs.push_back(U->seq1);
				myseqs.push_back(courNoeud->seqi_n);
				cost = seq->evaluation(myseqs,R->vehicle);
			}
			else cost = 1.e30 ;

			if ( cost < leastCost - EPSILON_LS )
			{
				leastCost = cost ;
				R->placeInsertionClient[U->cour][pattern] = courNoeud->pred ;
			}
			courNoeud = courNoeud->suiv ;
		}
		R->coutInsertionClient[U->cour][pattern] = leastCost - seq->evaluation(U->pred->seq0_i,U->suiv->seqi_n,R->vehicle);
	}

	// PCARP, setting back the pre-processed demand to its correct value
	U->seq1->load = tempDemand ;
}

void LocalSearch::melangeParcours ()
{
	// Shuffling the ordreParcours vector
	int j,temp ;
	for (int k = 0 ; k <= params->nbDays ; k++)
	{
		for (int i = 0 ; i < (int)ordreParcours[k].size() - 1 ; i++)
		{
			j = i + rand() % ((int)ordreParcours[k].size() - i) ;
			temp = ordreParcours[k][i] ;
			ordreParcours[k][i] = ordreParcours[k][j] ;
			ordreParcours[k][j] = temp ;
		}
	}
}

void LocalSearch::updateMoves ()
{
	int client, client2, size ;
	for (int k=1 ; k<=params->nbDays ; k++)
	{
		for (int i=0 ; i<( int)ordreParcours[k].size() ; i++)
		{
			client = ordreParcours[k][i] ;
			clients[k][client].moves.clear();
			size = (int)params->cli[client].sommetsVoisinsAvant.size();
			for (int a1 = 0 ; a1 < size ; a1++ )
			{
				client2 = params->cli[client].sommetsVoisinsAvant[a1] ;
				if (client2 >= params->nbDepots && clients[k][client2].estPresent) clients[k][client].moves.push_back(client2);
			}
		}
	}
}

void LocalSearch::setRouteVide(int day)
{
	routeVide[day] = NULL ;
	int route = 0 ;
	while (routeVide[day] == NULL && route < params->nbVehiculesPerDep )
	{
		if (depots[day][route].suiv->estUnDepot) routeVide[day] = depots[day][route].route ;
		route ++ ;
	}
}

void LocalSearch::nodeTestedForEachRoute (int cli, int day)
{
	for (int route = 0 ; route < (int)params->nombreVehicules[day] ; route ++)
		routes[day][route].nodeAndRouteTested[cli]=true ;	
}

void LocalSearch::placeManquants ()
{
	// Insertion procedure for the crossover PIX
	int k ;
	pattern meilleurPattern ;
	double depense, meilleureDepense ;
	int indexMeilleur = -1;
	Noeud * noeudTravail ;
	int calcul1, calcul2 ;
	pattern pattern1, pattern2 ;

	firstLoop = true ;

	for (int day = 1 ; day <= params->nbDays ; day++)
		for (int r=0 ; r < params->nombreVehicules[day] ; r++)
			routes[day][r].initiateInsertions() ;

	// We iterate on missing visits
	for (int i=0 ; i < (int)individu->toPlace.size() ; i++ )
	{
		meilleureDepense = 1.e30 ;
		k = individu->toPlace[i] ;
		pattern1 = individu->chromP[k] ;

		for (int pat = 0 ; pat < (int)params->cli[k].visits.size() ; pat ++)
		{
			// testing a new pattern
			pattern2 = params->cli[k].visits[pat] ;
			if (pattern::isSubset(pattern1,pattern2))
			{
				calcul2 = pattern2.pat ;
				depense = pattern2.cost ;
				for (int kk = 0 ; kk < params->ancienNbDays ; kk++)
				{
					if (calcul2 % 2 == 1)
					{
						noeudTravail = &clients[params->ancienNbDays-kk + pattern2.dep*params->ancienNbDays][k] ;
						computeCoutInsertion(noeudTravail,pat) ;
						depense += noeudTravail->coutInsertion[pat] ;
					}
					calcul2 = calcul2/2 ;
				}

				if (depense < meilleureDepense)
				{
					meilleureDepense = depense ;
					meilleurPattern = pattern2 ;
					indexMeilleur = pat ;
				}
			}
		}

		// Updating the chromP with the chosen pattern
		individu->chromP[k] = meilleurPattern ;

		// Inserting in the new locations
		calcul1 = pattern1.pat ;
		calcul2 = meilleurPattern.pat ;
		for (int kk = 0 ; kk < params->ancienNbDays ; kk++)
		{
			if (calcul2 % 2 == 1 && calcul1 % 2 == 0)
			{
				Noeud * hereCli = &clients[params->ancienNbDays-kk+ meilleurPattern.dep*params->ancienNbDays][k] ;
				Noeud * thereCli = hereCli->placeInsertion[indexMeilleur] ;
				addNoeud(hereCli,thereCli);
			}
			calcul1 = calcul1/2 ;
			calcul2 = calcul2/2 ;
		}
	}
}

void LocalSearch::insertNoeud(Noeud * U, Noeud * V)
{
	if (U->pred != V && U != V)
	{
		U->pred->suiv = U->suiv ;
		U->suiv->pred = U->pred ;
		V->suiv->pred = U ;
		U->pred = V ;
		U->suiv = V->suiv ;
		V->suiv = U ;
		U->route = V->route ;
	}
}

void LocalSearch::removeNoeud(Noeud * U)
{
	Noeud * temp =  U->suiv ;
	reinitSingleDayMoves(U->route);
	U->pred->suiv = U->suiv ;
	U->suiv->pred = U->pred ;
	U->route = NULL ;
	U->pred = NULL ;
	U->suiv = NULL ;
	temp->route->updateRouteData(false);
	setRouteVide(U->jour);

	// Managing the other data structures
	individu->chromT[U->jour].pop_back();
	removeOP(U->jour,U->cour);
	U->estPresent = false ;

	// Say that the insertions on this day need to be computed again
	temp->route->initiateInsertions();
	for (int i= params->nbDepots ; i< params->nbDepots + params->nbClients ; i++)
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			clients[U->jour][i].coutInsertion[p] = 1.e30 ;
}

void LocalSearch::addNoeud(Noeud * U, Noeud * V)
{
	reinitSingleDayMoves(V->route);
	V->suiv->pred = U ;
	U->pred = V;
	U->suiv = V->suiv ;
	V->suiv = U ;

	// Update the routes
	U->estPresent = true ;
	U->route = V->route ;
	U->route->updateRouteData(false);
	setRouteVide(U->jour);

	// Manage the other data structures
	individu->chromT[U->jour].push_back(0);
	addOP(U->jour,U->cour);

	// Say that the insertions on this day need to be computed again
	U->route->initiateInsertions();
	for (int i= params->nbDepots ; i< params->nbDepots + params->nbClients ; i++)
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			clients[U->jour][i].coutInsertion[p] = 1.e30 ;
}

void LocalSearch::removeOP (int day, int client) 
{
	// Remove one occurence in the ordreParcours structure
	int it = 0  ;
	while (ordreParcours[day][it] != client) { it ++ ; }
	ordreParcours[day][it] = ordreParcours[day][(int)ordreParcours[day].size() - 1];
	ordreParcours[day].pop_back();
}

void LocalSearch::addOP (int day, int client) 
{
	// Add one element in the ordreParcours structure
	int it, temp2 ;
	if (ordreParcours[day].size() != 0)
	{
		it = (int)rand() % ordreParcours[day].size() ;
		temp2 = ordreParcours[day][it] ;
		ordreParcours[day][it] = client ;
		ordreParcours[day].push_back(temp2);
	}
	else
		ordreParcours[day].push_back(client);
}

void LocalSearch::reinitSingleDayMoves(Route * r)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		r->nodeAndRouteTested[i] = false ;

	for (Noeud * tempNoeud = r->depot->suiv ; !tempNoeud->estUnDepot ; tempNoeud = tempNoeud->suiv)
		for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
			depots[tempNoeud->jour][route].route->nodeAndRouteTested[tempNoeud->cour] = false ;
}

void LocalSearch::reinitAllSingleDayMoves()
{
	for (int k=1 ; k<= params->nbDays ; k++)
		for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
			for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
				depots[k][route].route->nodeAndRouteTested[i] = false ;
}

bool compPredicateEC(EC_element * i,EC_element * j)
{
	return (i->cost < j->cost - 0.0001 || (i->cost < j->cost + 0.0001 && i->nbCustNodesInChain > j->nbCustNodesInChain));
}

bool LocalSearch::ejectionChains (int day)
{
	SeqData * seq = depots[day][0].seq0_i ;
	int myNodeIndex ;
	int myRouteIndex ;
	Noeud * myNoeud ;
	Vehicle * myVehicle ;
	Noeud * myDepot ;
	double myPrevCost ;
	double myDeltaCost ;

	// 1) CHOOSE AN ORDER FOR THE BINS
	vector < Noeud * > ordreBins ;
	for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
		ordreBins.push_back(&depots[day][route]);
	std::random_shuffle(ordreBins.begin(),ordreBins.end());

	// 2) UPDATE THE DATA STRUCTURE WILL ALL NECESSARY INFORMATIONS AND THE GOOD SIZE
	myRouteIndex = 0 ;
	for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
	{
		myNoeud = ordreBins[route] ;
		// if the route is not empty, this makes a new layer in the ejection chains graph
		if (!myNoeud->suiv->estUnDepot)
		{
			ejectionGraph[myRouteIndex][0].myNode = myNoeud ;
			ejectionGraph[myRouteIndex][0].cost = 1.e20 ;
			ejectionGraph[myRouteIndex][0].pred = NULL ;
			ejectionGraph[myRouteIndex][0].nbCustNodesInChain = 0 ;
			ejectionGraph[myRouteIndex][0].routeID = route ;
			myNoeud = myNoeud->suiv ;
			myNodeIndex = 1 ;
			while (!myNoeud->estUnDepot)
			{
				ejectionGraph[myRouteIndex][myNodeIndex].myNode = myNoeud ;
				ejectionGraph[myRouteIndex][myNodeIndex].cost = 1.e20 ;
				ejectionGraph[myRouteIndex][myNodeIndex].pred = NULL ;
				ejectionGraph[myRouteIndex][myNodeIndex].nbCustNodesInChain = 0 ;
				ejectionGraph[myRouteIndex][myNodeIndex].routeID = route ;
				myNoeud = myNoeud->suiv ;
				myNodeIndex ++ ;
			}
			ec_nbElements[myRouteIndex] = myNodeIndex ;
			myRouteIndex ++ ;
		}
	}
	ec_nbRoutes = myRouteIndex ;

	if (ec_nbRoutes < 2)
	{
		//cout << "Only one route for EC, abort" << endl ;
		return false ;
	}

	// 3) SOLVE THE SHORTEST PATH PROBLEM
	// for each bin in order
	for (int r=0 ; r < ec_nbRoutes ; r++)
	{
		myDepot = ejectionGraph[r][0].myNode ;
		myVehicle = myDepot->route->vehicle ;
		myPrevCost = myDepot->pred->seq0_i->evaluation(myDepot->pred->seq0_i,myVehicle); 

		for (int rPrev=0 ; rPrev < r  ; rPrev ++)
		{
			// For the depot node, propagate from another previous depot node
			if (ejectionGraph[rPrev][0].cost < ejectionGraph[r][0].cost - 0.0001 || 
				(ejectionGraph[rPrev][0].cost < ejectionGraph[r][0].cost + 0.0001 && ejectionGraph[rPrev][0].nbCustNodesInChain > ejectionGraph[r][0].nbCustNodesInChain))
			{
				ejectionGraph[r][0].cost = ejectionGraph[rPrev][0].cost ;
				ejectionGraph[r][0].pred = &ejectionGraph[rPrev][0] ;
				ejectionGraph[r][0].bestInsertionPlace = NULL ;
				ejectionGraph[r][0].nbCustNodesInChain = ejectionGraph[rPrev][0].nbCustNodesInChain ;
			}

			// for the depot node, propagate from any predecessor node
			for (int iiPrev=1 ; iiPrev < ec_nbElements[rPrev] ; iiPrev++)
			{
				// for any such predecessor, test the best place of insertion in the route
				for (int ii=0 ; ii < ec_nbElements[r] ; ii++)
				{
					myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->seq0_i,
						ejectionGraph[rPrev][iiPrev].myNode->seq1,
						ejectionGraph[r][ii].myNode->suiv->seqi_n,
						myVehicle) - myPrevCost ; 
					if (myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][0].cost - 0.0001 || 
						(myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][0].cost + 0.0001 && ejectionGraph[rPrev][iiPrev].nbCustNodesInChain > ejectionGraph[r][0].nbCustNodesInChain))
					{
						ejectionGraph[r][0].cost = myDeltaCost + ejectionGraph[rPrev][iiPrev].cost ;
						ejectionGraph[r][0].pred = &ejectionGraph[rPrev][iiPrev] ;
						ejectionGraph[r][0].bestInsertionPlace = ejectionGraph[r][ii].myNode ;
						ejectionGraph[r][0].nbCustNodesInChain = ejectionGraph[rPrev][iiPrev].nbCustNodesInChain ;
					}
				}
			}
		}

		for (int ii=1 ; ii < ec_nbElements[r] ; ii++)
		{
			// for each customer node in the bin, propagate from the general 0 node
			myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->pred->seq0_i,
				ejectionGraph[r][ii].myNode->suiv->seqi_n,
				myVehicle) - myPrevCost ;
			if (myDeltaCost < ejectionGraph[r][ii].cost - 0.0001)
			{
				ejectionGraph[r][ii].cost = myDeltaCost ;
				ejectionGraph[r][ii].pred = NULL ;
				ejectionGraph[r][ii].bestInsertionPlace = NULL ;
				ejectionGraph[r][ii].nbCustNodesInChain = 1 ;
			}

			// for each customer node in the bin, propagate from any previous 0 (depot) node.
			for (int rPrev=0 ; rPrev < r  ; rPrev ++)
			{
				if (myDeltaCost + ejectionGraph[rPrev][0].cost < ejectionGraph[r][ii].cost - 0.0001 ||
					(myDeltaCost + ejectionGraph[rPrev][0].cost < ejectionGraph[r][ii].cost + 0.0001 && ejectionGraph[rPrev][0].nbCustNodesInChain >= ejectionGraph[r][ii].nbCustNodesInChain))
				{
					ejectionGraph[r][ii].cost = myDeltaCost + ejectionGraph[rPrev][0].cost ;
					ejectionGraph[r][ii].pred = &ejectionGraph[rPrev][0] ;
					ejectionGraph[r][ii].bestInsertionPlace = NULL ;
					ejectionGraph[r][ii].nbCustNodesInChain = ejectionGraph[rPrev][0].nbCustNodesInChain + 1 ;
				}
			}

			// for each customer node in the bin, propagate from a previous customer node.
			for (int rPrev=0 ; rPrev < r  ; rPrev ++)
			{
				for (int iiPrev=1 ; iiPrev < ec_nbElements[rPrev] ; iiPrev++)
				{
					myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->pred->seq0_i,
						ejectionGraph[rPrev][iiPrev].myNode->seq1,
						ejectionGraph[r][ii].myNode->suiv->seqi_n,
						myVehicle) - myPrevCost ;

					if (myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][ii].cost - 0.0001 ||
						(myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][ii].cost + 0.0001 &&  ejectionGraph[rPrev][iiPrev].nbCustNodesInChain >= ejectionGraph[r][ii].nbCustNodesInChain))
					{
						ejectionGraph[r][ii].cost = myDeltaCost + ejectionGraph[rPrev][iiPrev].cost ;
						ejectionGraph[r][ii].pred = &ejectionGraph[rPrev][iiPrev] ;
						ejectionGraph[r][ii].bestInsertionPlace = NULL ;
						ejectionGraph[r][ii].nbCustNodesInChain = ejectionGraph[rPrev][iiPrev].nbCustNodesInChain + 1 ;
					}
				}
			}
		}
	}

	// 4) KEEP THE BEST SHORTEST PATHS in increasing order of costs.
	vector <EC_element *> orderEnds ;
	for (int r=1 ; r < ec_nbRoutes ; r++)
		orderEnds.push_back(&ejectionGraph[r][0]);
	std::sort(orderEnds.begin(),orderEnds.end(),compPredicateEC);

	// 5) APPLY THE MOVES TO UPDATE THE SOLUTION (track back the nodes and apply successive relocates).
	// Apply the best chain
	if (orderEnds[0]->cost < -0.0001)
	{
		EC_element * elementCour ;
		EC_element * elementPred ;
		//cout << "----------------- EJECTION CHAINS PHASE ---------------" << endl ;
		for (int myID = 0 ; myID < 1 ; myID++)
		{
			nbEjectionChains ++ ;
			nbEjectionChainsNodes += orderEnds[myID]->nbCustNodesInChain ;
			elementCour = orderEnds[myID];
			elementPred = elementCour->pred ;
			Noeud * insertionPosition = NULL;
			Noeud * insertionPositionTemp = NULL;
			while (elementCour != NULL)
			{
				if (!(elementPred == NULL || elementPred->myNode->estUnDepot)) 
				{
					//before is a node, after is a depot (insert in the best location, and keep track of where he was.)
					if (elementCour->myNode->estUnDepot)
					{
						insertionPosition = elementPred->myNode->pred ;		
						insertNoeud(elementPred->myNode,elementCour->bestInsertionPlace);
						reinitSingleDayMoves(elementCour->bestInsertionPlace->route);
						reinitSingleDayMoves(insertionPosition->route);
						elementCour->bestInsertionPlace->route->updateRouteData(false);
						insertionPosition->route->updateRouteData(false);
					}
					else // before is a node, after is a node
					{
						insertionPositionTemp = elementPred->myNode->pred ;
						insertNoeud(elementPred->myNode,insertionPosition);
						reinitSingleDayMoves(insertionPosition->route);
						reinitSingleDayMoves(insertionPositionTemp->route);
						insertionPosition->route->updateRouteData(false);
						insertionPositionTemp->route->updateRouteData(false);
						insertionPosition = insertionPositionTemp ;
					}
				}
				elementCour = elementPred ;
				if (elementPred != NULL) elementPred = elementCour->pred ;
			}
		}
		setRouteVide(day);
		return true ;
	}
	else
		return false ;
}

void LocalSearch::addSeqDataInPieces (Noeud * node, int length, int day)
{
	Noeud * cour = node ;
	SeqData * courSeq ;
	int courLenght = length + 1 ;
	int temp ;
	while (courLenght > 0)
	{
		temp = min(courLenght,params->sizeSD) ;
		courSeq = cour->seqi_j[temp-1];
		myseqs.push_back(courSeq);
		courLenght -= temp ;
		if (courLenght != 0) cour = clients[day][courSeq->lastNode].suiv ;
	}
}

void LocalSearch::addReverseSeqDataInPieces (Noeud * node, int length, int day)
{
	Noeud * cour = node ;
	SeqData * courSeq ;
	list <SeqData*> myseqTemp ;
	int courLenght = length + 1 ;
	int temp ;
	while (courLenght > 0)
	{
		temp = min(courLenght,params->sizeSD) ;
		courSeq = cour->seqj_i[temp-1];
		myseqTemp.push_back(courSeq);
		courLenght -= temp ;
		if (courLenght != 0) cour = clients[day][courSeq->firstNode].suiv ;
	}

	for ( list<SeqData*>::reverse_iterator it=myseqTemp.rbegin() ; it != myseqTemp.rend() ; it++)
		myseqs.push_back(*it);
}

LocalSearch::LocalSearch(void)
{
	allAttributesSet = false ;
	seqdeb = NULL ;
}

LocalSearch::LocalSearch(Params * params,Individu * individu) : params (params),individu(individu)
{
	nbDays = params->nbDays ;
	int nbVeh ;
	allAttributesSet = true ;
	nbInterSwap = 0 ;
	nbIntraSwap = 0 ;
	nbInter2Opt = 0 ;
	nbIntra2Opt = 0 ;
	nbEjectionChains = 0 ;
	nbEjectionChainsNodes = 0 ;
	seqdeb = NULL ;
	nbTotalRISinceBeginning = 0 ;
	nbTotalPISinceBeginning = 0 ;
	vector < Noeud * > tempNoeud ;
	vector <int> temp2 ;

	clients = new Noeud * [params->nbDays+1] ;
	depots = new Noeud * [params->nbDays+1] ;
	depotsFin = new Noeud * [params->nbDays+1] ;
	routes = new Route * [params->nbDays+1] ;

	for (int kk = 1 ; kk <= params->nbDays ; kk++)
	{
		nbVeh = params->nombreVehicules[kk] ;
		clients[kk] = new Noeud [params->nbClients + params->nbDepots + 1] ;
		depots[kk] = new Noeud [nbVeh] ;
		depotsFin[kk] = new Noeud [nbVeh] ;
		routes[kk] = new Route [nbVeh] ; 

		for (int i = 0 ; i <  params->nbClients + params->nbDepots ; i ++ )
			clients[kk][i] = Noeud(false,i,kk,false,NULL,NULL,NULL,params);

		for (int i = 0 ; i < nbVeh ; i++ )
		{
			depots[kk][i] = Noeud(true,params->ordreVehicules[kk][i].depotNumber,kk,false,NULL,NULL,NULL,params);
			depotsFin[kk][i] = Noeud(true,params->ordreVehicules[kk][i].depotNumber,kk,false,NULL,NULL,NULL,params);
			routes[kk][i] = Route(i,0,&params->ordreVehicules[kk][i],params,individu,kk);
			depots[kk][i].route = &routes[kk][i] ;
			depotsFin[kk][i].route = &routes[kk][i] ;
			routes[kk][i].depot = &depots[kk][i] ;
		}
	}

	for (int i=0 ; i < 4 ; i++)
	{
		resultMoves.push_back(vector<double>(4));
		shouldBeTested.push_back(vector<bool>(4));
	}

	for (int day = 0 ; day <= params->nbDays ; day++)
	{
		ordreParcours.push_back(temp2);
		routeVide.push_back(NULL);
	}

	for (int i=params->nbDepots ; i < params->nbDepots + params->nbClients ; i++)
		ordreParcours[0].push_back(i);

	// Initialization for ejection chains
	for (int v=0 ; v < params->nbVehiculesPerDep ; v ++)
		ejectionGraph.push_back(vector <EC_element> (params->nbClients+1));
	ec_nbElements = vector <int> (params->nbVehiculesPerDep);

	// Initialization of the SeqData structures, for all nodes
	// These structures were designed to stay contiguous in memory (however it had only a negligible impact on performance)
	int nbSeqsSet = 0 ;
	int taillemyseqDatas = (params->sizeSD + params->sizeSD + 4)*(params->nbClients + params->nbDepots + 2*params->nbVehiculesPerDep)*params->nbDays ;
	SeqData * myseqDatas = new SeqData [taillemyseqDatas] ;
	seqdeb = myseqDatas ;
	for (int k=1 ; k <= params->nbDays ; k++)
	{
		for (int i=0 ; i < params->nbClients + params->nbDepots ; i++) 
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(clients[k][i].cour,params,individu,k,false);  

			clients[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			clients[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			clients[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			clients[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				clients[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				clients[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			clients[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}

		for (int i=0 ; i < params->nbVehiculesPerDep  ; i++) 
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(depots[k][i].cour,params,individu,k,false);

			depots[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			depots[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			depots[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			depots[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depots[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depots[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			depots[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}

		for (int i=0 ; i < params->nbVehiculesPerDep ; i++) 
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(depotsFin[k][i].cour,params,individu,k,false);

			depotsFin[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			depotsFin[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			depotsFin[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			depotsFin[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depotsFin[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depotsFin[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			depotsFin[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}
	}
}

LocalSearch::~LocalSearch(void)
{	
	SeqData * seqdeb2 = (SeqData*) seqdeb ;
	delete [] seqdeb2 ;

	if (allAttributesSet)
	{
		for (int kk = 1 ; kk <= nbDays ; kk++)
		{
			delete [] clients[kk] ;
			delete [] depots[kk] ;
			delete [] depotsFin[kk] ;
			delete [] routes[kk] ;
		}
		delete [] clients ;
		delete [] depots ;
		delete [] depotsFin ;
		delete [] routes ;
	}
}

