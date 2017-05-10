CXX=g++
INC=../include
SRC=../src
SRC_EXE=../src_exe

all: objects lib tester

lib: objects lib
	$(CXX) -shared -o libAtmosphere.so Atmosphere.o

objects: $(INC)/Atmosphere.h $(SRC)/Atmosphere.cc
	$(CXX) -c -fPIC $(SRC)/Atmosphere.cc -I$(INC) -o Atmosphere.o

install: all

clean:
	rm -f *.o *.so a.out tester

tester: $(SRC_EXE)/tester.cc
	$(CXX) $(SRC_EXE)/tester.cc -o tester -m64 -L./ -lAtmosphere -I$(INC)
