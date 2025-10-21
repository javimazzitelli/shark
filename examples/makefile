LIBS = -lgsl -lgslcblas -lm -lstdc++ -fno-math-errno
CXX = gcc #-std=gnu++11
SHARKDIR = ..

GENERICMAIN = example

all: $(TTBARMAIN) $(GENERICMAIN)

example.o: example.cpp
	$(CXX) -I$(SHARKDIR)/inc -c example.cpp

$(GENERICMAIN): $(LIBshark) example.o
	$(CXX) example.o -L$(SHARKDIR)/lib -Wl,-rpath,$(SHARKDIR)/lib -L$(SHARKDIR)/lib -lshark $(LIBS) -o $(GENERICMAIN)

clean:
	rm *.o example
