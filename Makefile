CC=g++
SRCDIR=src
SRCEXT=cpp
BUILDDIR=build
EXAMPLESDIR=examples
BENCHMARKSDIR=benchmark
INCLUDEDIR=include
TOOLSDIR=tools

#CFLAGS=-g -W -Wall -O0 -DDEBUG -std=c++11
CFLAGS=-W -Wall -O3 -std=c++11 -pg
#CFLAGS=-W -Wall -O0 -std=c++11  -g

SOURCES=$(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS=$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
EXAMPLES=$(patsubst $(EXAMPLESDIR)/%.$(SRCEXT),%,$(shell find $(EXAMPLESDIR) -type f -name *.$(SRCEXT)))
BENCHMARKS=$(patsubst $(BENCHMARKSDIR)/%.$(SRCEXT),%,$(shell find $(BENCHMARKSDIR) -type f -name *.$(SRCEXT)))
TOOLS=$(patsubst $(TOOLSDIR)/%.$(SRCEXT),%,$(shell find $(TOOLSDIR) -type f -name *.$(SRCEXT)))
LIB=-lsdsl -ldivsufsort -ldivsufsort64 -lboost_filesystem -lboost_system -lboost_iostreams
INC=-I include

all: $(OBJECTS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(INCLUDEDIR)/%.hpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) bin/*"; $(RM) -r $(BUILDDIR) bin/*

examples: $(EXAMPLES)

benchmarks: $(BENCHMARKS)

tools: $(TOOLS)

$(EXAMPLES): $(OBJECTS)
	@echo " Building examples...";
	@echo " $(CC) $(CFLAGS) $^ -o bin/$@ $(EXAMPLESDIR)/$@.$(SRCEXT) $(INC) $(LIB)"; $(CC) $(CFLAGS) $^ -o bin/$@ $(EXAMPLESDIR)/$@.$(SRCEXT) $(INC) $(LIB)

$(BENCHMARKS): $(OBJECTS)
	@echo " Building benchmarks...";
	@echo " $(CC) $(CFLAGS) $^ -o bin/$@ $(BENCHMARKSDIR)/$@.$(SRCEXT) $(INC) $(LIB) -lGkArrays -lz"; $(CC) $(CFLAGS) $^ -o bin/$@ $(BENCHMARKSDIR)/$@.$(SRCEXT) $(INC) $(LIB) -lGkArrays -lz

$(TOOLS): $(OBJECTS)
	@echo " $(CC) $(CFLAGS) $^ -o bin/$@ $(TOOLSDIR)/$@.$(SRCEXT) $(INC) $(LIB)"; $(CC) $(CFLAGS) $^ -o bin/$@ $(TOOLSDIR)/$@.$(SRCEXT) $(INC) $(LIB)
