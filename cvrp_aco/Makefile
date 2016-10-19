# Makefile for acovrp
VERSION=1.03

OPTIM_FLAGS=-O
WARN_FLAGS=-Wall -ansi -pedantic
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS)
# To change the default timer implementation, uncomment the line below
# or call 'make TIMER=unix'
TIMER=unix
LDLIBS=-lm

cvrp_aco: vrpHelper.o utilities.o problem.o main.o $(TIMER)_timer.o parallelAco.o antColony.o io.o localSearch.o

all: clean cvrp_aco

clean:
	@$(RM) *.o cvrp_aco

cvrp_aco: vrpHelper.o utilities.o problem.o main.o $(TIMER)_timer.o parallelAco.o antColony.o io.o localSearch.o

parallelAco.o: parallelAco.cpp parallelAco.h

antColony.o: antColony.cpp antColony.h

acovrp.o: acovrp.c

vrp.o: vrp.c vrp.h

ants.o: ants.c ants.h

InOut.o: InOut.c InOut.h

utilities.o: utilities.c utilities.h

ls.o: ls.c ls.h

$(TIMER)_timer.o: $(TIMER)_timer.c timer.h

dist : DIST_SRC_FILES=*.c *.h *.vrp Makefile gpl.txt
dist : all
	@(mkdir -p ../acovrp-$(VERSION)			\
	&& rsync -rlpC --exclude=.svn $(DIST_SRC_FILES) ../acovrp-$(VERSION)/ \
        && cd .. 	\
	&& tar cf - acovrp-$(VERSION) | gzip -f9 > acovrp-$(VERSION).tar.gz \
	&& rm -rf ./acovrp-$(VERSION)					\
	&& echo "acovrp-$(VERSION).tar.gz created." && cd $(CWD) )
