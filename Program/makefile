all : gencarp

CCC = g++
CCFLAGS = -O3 -Wall  
LIBS= -lm
TARGETDIR=.
CPPFLAGS += \
         -I.
        
OBJS2 = \
        $(TARGETDIR)/Client.o \
        $(TARGETDIR)/Vehicle.o \
        $(TARGETDIR)/commandline.o \
	    $(TARGETDIR)/Genetic.o \
        $(TARGETDIR)/Individu.o \
        $(TARGETDIR)/LocalSearch.o \
        $(TARGETDIR)/main.o \
        $(TARGETDIR)/Noeud.o \
        $(TARGETDIR)/SeqData.o \
        $(TARGETDIR)/Params.o \
        $(TARGETDIR)/Population.o \
        $(TARGETDIR)/Route.o
        
	
$(TARGETDIR)/gencarp: $(OBJS2)
	$(CCC)  $(CCFLAGS) $(LIBS) -o $(TARGETDIR)/gencarp $(OBJS2)

$(TARGETDIR)/Client.o: Client.h Client.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Client.cpp -o $(TARGETDIR)/Client.o

$(TARGETDIR)/Vehicle.o: Vehicle.h Vehicle.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Vehicle.cpp -o $(TARGETDIR)/Vehicle.o
	
$(TARGETDIR)/commandline.o: commandline.h commandline.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c commandline.cpp -o $(TARGETDIR)/commandline.o
	
$(TARGETDIR)/Genetic.o: Genetic.h Genetic.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Genetic.cpp -o $(TARGETDIR)/Genetic.o

$(TARGETDIR)/Individu.o: Individu.h Individu.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Individu.cpp -o $(TARGETDIR)/Individu.o

$(TARGETDIR)/LocalSearch.o: LocalSearch.h LocalSearch.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c LocalSearch.cpp -o $(TARGETDIR)/LocalSearch.o
	
$(TARGETDIR)/main.o: main.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c main.cpp -o $(TARGETDIR)/main.o
	
$(TARGETDIR)/Noeud.o: Noeud.h Noeud.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Noeud.cpp -o $(TARGETDIR)/Noeud.o

$(TARGETDIR)/SeqData.o: SeqData.h
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c SeqData.cpp -o $(TARGETDIR)/SeqData.o

$(TARGETDIR)/Params.o: Params.h Params.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Params.cpp -o $(TARGETDIR)/Params.o

$(TARGETDIR)/Population.o: Population.h Population.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Population.cpp -o $(TARGETDIR)/Population.o

$(TARGETDIR)/Route.o: Route.h Route.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Route.cpp -o $(TARGETDIR)/Route.o

test: gencarp
	./gencarp ../Instances/CARP/gdb1.dat -type 30 -sol testCARP.sol -seed 1
	./gencarp ../Instances/MCGRP/mgval_0.25_1A.dat -type 31 -sol testMCGRP.sol -seed 1
	./gencarp ../Instances/PCARP/gdb10.txt -type 32  -sol testPCARP.sol -seed 1 -veh 10
	./gencarp ../Instances/CARP/egl-e1-A.dat -type 33 -sol testMDCARP.sol -seed 1 -dep 4
	./gencarp ../Instances/MM-kWRPP/P202100 -type 35 -sol testMM-kWRPP.sol -seed 1 -veh 3

clean:
	 $(RM) \
	$(TARGETDIR)/gencarp \
	$(TARGETDIR)/Client.o \
        $(TARGETDIR)/Vehicle.o \
       $(TARGETDIR)/main.o \
     $(TARGETDIR)/commandline.o \
	 $(TARGETDIR)/Genetic.o \
     $(TARGETDIR)/Individu.o \
     $(TARGETDIR)/LocalSearch.o \
     $(TARGETDIR)/Noeud.o \
     $(TARGETDIR)/SeqData.o \
     $(TARGETDIR)/Params.o \
     $(TARGETDIR)/Population.o \
     $(TARGETDIR)/Route.o



