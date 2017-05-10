objects: Atmosphere.h Atmosphere.cc
	clang -c -fPIC Atmosphere.cc -o Atmosphere.o

all: objects
	clang -shared -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-install_name,/Users/harmscho/OwnLibs/Atmosphere/lib/libAtmosphere.so -o libAtmosphere.so Atmosphere.o

install: all
	mv *.so ../lib/
	cp *.h ../include/

clean:
	rm -f *.o *.dylib a.out tester

tester: tester.cc install
	g++ tester.cc -o tester -m64 -L/Users/harmscho/OwnLibs/Atmosphere/lib -lAtmosphere -I/Users/harmscho/OwnLibs/Atmosphere/include
