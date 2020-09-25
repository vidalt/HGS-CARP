all : gencarp

CCC = g++
CCFLAGS =-DTURN_PENALTIES=1 -Wall -std=c++14 -g -Ofast -I/usr/local/opt/openssl/include/
LDFLAGS =-Wall -std=c++14 -g -O0 -L/usr/local/opt/openssl/lib/
LIBS= -lm -lcrypto
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
	$(CCC)  $(CCFLAGS) $(LDFLAGS) $(LIBS) -o $(TARGETDIR)/gencarp $(OBJS2)

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



