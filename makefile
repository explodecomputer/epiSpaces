CC = gcc
CFLAGS = -c -O2 -Wall
LDFLAGS = -lm -Wall -O2
FCC = gfortran
SOURCES1 = main.c
OBJECTS1 = $(SOURCES1:.c=.o)
EXE1 = epiFitness

all: $(SOURCES1) $(EXE1)

$(EXE1): $(OBJECTS1)
	$(CC) $(LDFLAGS) $(OBJECTS1) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

all:
	rm $(OBJECTS1)

#ifeq ($(UNAME),Darwin)
#all:
#	$(FCC) -dynamiclib $(SOURCES2) -o $(SO)
#	rm $(OBJECTS1)
#else
#all:
#	$(FCC) -fPIC -c $(SOURCES2)
#	$(FCC) -shared $(OBJECTS2) -o $(SO)
#	rm $(OBJECTS1) $(OBJECTS2)
#endif
#clean:
#	rm $(OBJECTS1) $(EXE1)

