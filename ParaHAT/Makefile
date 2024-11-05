CC=mpic++
#CFLAGS=-O0 -g3 -ggdb -c -Wall
CFLAGS=-O3 -c -Wall
DFLAGS=-g -Wall
LDFLAGS= -lz -lm -lpthread

KSWSOURCE=ksw.c 
ALGNSOURCES=main.cpp aligner.cpp graph.cpp  readfl.cpp hash.cpp form.cpp utils.c
INDSOURCES=whash.cpp genhash.cpp formh.cpp

INDOBJECTS=$(INDSOURCES:.cpp=.o) $(KSWSOURCE:.c=.o)
ALGNOBJECTS=$(ALGNSOURCES:.cpp=.o) $(KSWSOURCE:.c=.o)

INDEXER=ParaHAT-indexer
ALIGNER=ParaHAT-aligner

all: $(INDSOURCES) $(ALGNSOURCES) $(KSWSOURCE) $(ALIGNER) $(INDEXER)
	
$(ALIGNER): $(ALGNOBJECTS)
	$(CC)  $(ALGNOBJECTS) -o $@ $(LDFLAGS)

$(INDEXER): $(INDOBJECTS)
	$(CC)  $(INDOBJECTS) readfl.o -o $@ $(LDFLAGS)

debug:
	$(CC) $(DFLAGS) $(ALGNSOURCES) $(KSWSOURCE) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.o $(ALIGNER) $(INDEXER) a.out

 
