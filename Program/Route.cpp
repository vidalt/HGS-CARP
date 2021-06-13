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

#include "Route.h"

Route::Route(void){}

Route::Route(int cour, Noeud * depot, Vehicle * vehicle, Params * params, Individu * indiv, int day) : params(params), individu(indiv), cour(cour), day(day), depot(depot), vehicle(vehicle) 
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		coutInsertionClient.push_back(vector <double> ()) ;
		placeInsertionClient.push_back(vector <Noeud *> ());
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
		{
			coutInsertionClient[i].push_back(1.e30);
			placeInsertionClient[i].push_back(NULL);
		}
		nodeAndRouteTested.push_back(false);
	}
}

Route::~Route(void){}

void Route::reverse ()
{
	// Reversing the order of the nodes in the sequence
	if (!depot->suiv->estUnDepot)
	{
		Noeud * temp ;
		Noeud * myDepot = depot ;
		Noeud * myDepotFin = depot->pred ;
		Noeud * current = depot->suiv ;

		while ( !current->estUnDepot )
		{
			temp = current->suiv ;
			current->suiv = current->pred ;
			current->pred = temp ;
			current = temp ;
		}

		temp = myDepot->suiv ;
		myDepot->suiv = myDepotFin->pred ;
		myDepotFin->pred = temp ;
		myDepot->suiv->pred = myDepot ;
		myDepotFin->pred->suiv = myDepotFin ;
	}
}

void Route::updateRouteData (bool isForPrinting)
{
	bool firstIt ;
	int place = 0 ;
	double Xvalue = 0 ;
	double Yvalue = 0 ;
	double nbNodes = 1 ;

	// Computing the auxiliary data on any subsequence (0..i), using forward recursion
	Noeud * noeud = depot ;
	noeud->place = place ;
	noeud->seq0_i->initialisation(noeud->cour,params,individu,day,isForPrinting);
	noeud->seqi_0->initialisation(noeud->cour,params,individu,day,false);
	Xvalue += params->cli[noeud->cour].coord.x ;
	Yvalue += params->cli[noeud->cour].coord.y ;

	firstIt = true ;
	while (!noeud->estUnDepot || firstIt)
	{
		firstIt = false ;
		noeud = noeud->suiv ;
		Xvalue += params->cli[noeud->cour].coord.x ;
		Yvalue += params->cli[noeud->cour].coord.y ;
		nbNodes ++ ;
		place ++ ;
		noeud->place = place ;
		if (isForPrinting)
			noeud->seq0_i->concatOneAfterWithPathTracking(noeud->pred->seq0_i,noeud->cour,individu,day);
		else
			noeud->seq0_i->concatOneAfter(noeud->pred->seq0_i,noeud->cour,individu,day);
		noeud->seqi_0->concatOneBefore(noeud->pred->seqi_0,noeud->cour,individu,day);
	}

	// Computing the auxiliary data on any subsequence (i..n), using backward recursion
	noeud = depot->pred ;
	noeud->seqi_n->initialisation(noeud->cour,params,individu,day,false);
	noeud->seqn_i->initialisation(noeud->cour,params,individu,day,false);

	firstIt = true ;
	while ( !noeud->estUnDepot || firstIt )
	{
		firstIt = false ;
		noeud = noeud->pred ;
		noeud->seqi_n->concatOneBefore(noeud->suiv->seqi_n,noeud->cour,individu,day);
		noeud->seqn_i->concatOneAfter(noeud->suiv->seqn_i,noeud->cour,individu,day);
	}

	// Computing the auxiliary data on any subsequence (i..j), using forward recursion
	// To gain a bit of time, we limit this preprocessing to subsequences such that i..j does not contain more than "sizeSD" elements
	// (More intelligent strategies could be used (e.g., the hierarchical approach of Irnich))
	Noeud * noeudi ;
	Noeud * noeudj ;
	noeudi = depot ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		noeudi->seqi_j[0]->initialisation(noeudi->cour,params,individu,day,false);
		noeudj = noeudi->suiv ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			noeudi->seqi_j[j]->concatOneAfter(noeudi->seqi_j[j-1],noeudj->cour,individu,day);
			noeudj = noeudj->suiv ;
		}
		noeudi = noeudi->suiv ;
	}

	noeudi = depot->pred ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		noeudi->seqj_i[0]->initialisation(noeudi->cour,params,individu,day,false);
		noeudj = noeudi->pred ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			noeudj->seqj_i[j]->concatOneAfter(noeudj->suiv->seqj_i[j-1],noeudj->cour,individu,day);
			noeudj = noeudj->pred ;
		}
		noeudi = noeudi->pred ;
	}

	// Checking the route feasibility
	double dist ;
	double violLoad ;
	double violDuration ;
	currentRouteCost = depot->pred->seq0_i->evaluation(depot->pred->seq0_i,depot->seq1,this->vehicle,dist,violDuration,violLoad);
	
	if (violDuration < 0.001 && violLoad < 0.001) 
		isFeasible = true ;
	else 
		isFeasible = false ;
}

// no insertion are computed
void Route::initiateInsertions()
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			coutInsertionClient[i][p] = 1.e30 ;
}

