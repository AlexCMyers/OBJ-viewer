## Compiler
CC=g++
## Linker
LD=$(CC)
## Flags
CPPFLAGS = -std=c++11 -Wall -g -DLINUX
LFLAGS = -lglut -L/usr/lib -L/usr/X11R6/lib -lXi -lXmu -lGL -lGLU -lm 

TARGETS = $(PROGFILES:.cpp=)

PROGFILES = \
        assn2.cpp \
        $(NULL)

targets default: $(TARGETS)

$(PROGFILES:.cpp=): assn2.o 
	$(CC) -o assn2 assn2.o ${LFLAGS}

depend :
	makedepend ${PROGFILES}
# DO NOT DELETE
