CPP=g++
CXXFLAGS=-O3 -march=native
LIBS= -lntl -lgmp

%.o : %.cpp
	$(CPP) $(CXXFLAGS) -c -o $@ $<
 
bench : bench.o relax.o middle_product.o short_product.o
	libtool --mode=link $(CPP) $(CXXFLAGS) -o $@ $^ $(LIBS)

test : test.o relax.o middle_product.o short_product.o
	libtool --mode=link $(CPP) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o bench test

all: bench test
