# Makefile for cvrp_aco
VERSION=1.0

OPTIM_FLAGS=-O2
WARN_FLAGS=-Wall -ansi -pedantic
CPPFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS) -lpthread
CC=g++
# To change the default timer implementation, uncomment the line below
# or call 'make TIMER=unix'
TIMER=unix
LDLIBS=-lm

OBJS= antColony.o io.o localSearch.o main.o $(TIMER)_timer.o move.o neighbourSearch.o parallelAco.o problem.o simulatedAnnealing.o utilities.o vrpHelper.o
EXE=main

all: clean cvrp_aco

clean:
	@$(RM) *.o main

cvrp_aco: $(OBJS)
	$(CC) $(CPPFLAGS) $(OBJS) -o $(EXE)

antColony.o: antColony.cpp antColony.h

io.o: io.cpp io.h

localSearch.o: localSearch.cpp localSearch.h

main.o: main.cpp

$(TIMER)_timer.o: $(TIMER)_timer.cpp timer.h

move.o: move.cpp move.h

neighbourSearch.o: neighbourSearch.cpp neighbourSearch.h

parallelAco.o: parallelAco.cpp parallelAco.h

problem.o: problem.cpp problem.h

simulatedAnnealing.o: simulatedAnnealing.cpp simulatedAnnealing.h

utilities.o: utilities.cpp utilities.h

vrpHelper.o: vrpHelper.cpp vrpHelper.h


dist : DIST_SRC_FILES=*.cpp *.h Makefile
dist : all
	@(mkdir -p ../cvrp_aco-$(VERSION)			\
	&& rsync -rlpC --exclude=.svn $(DIST_SRC_FILES) ../cvrp_aco-$(VERSION)/ \
        && cd .. 	\
	&& tar cf - cvrp_aco-$(VERSION) | gzip -f9 > cvrp_aco-$(VERSION).tar.gz \
	&& rm -rf ./cvrp_aco-$(VERSION)					\
	&& echo "cvrp_aco-$(VERSION).tar.gz created." && cd $(CWD) )
