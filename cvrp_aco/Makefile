# Makefile for acovrp
VERSION=1.03

OPTIM_FLAGS=-O
WARN_FLAGS=-Wall -ansi -pedantic
CFLAGS=$(WARN_FLAGS) $(OPTIM_FLAGS)
# To change the default timer implementation, uncomment the line below
# or call 'make TIMER=unix'
TIMER=dos
#TIMER=unix
LDLIBS=-lm

acovrp: acovrp.o vrp.o utilities.o ants.o InOut.o $(TIMER)_timer.o ls.o

all: clean acovrp

clean:
	@$(RM) *.o acovrp

acovrp: acovrp.o vrp.o utilities.o ants.o InOut.o $(TIMER)_timer.o ls.o

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
