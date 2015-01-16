CXXFLAGS=$(shell root-config --cflags) -fPIC 
ROOTLIBS=$(shell root-config --libs)

BASE=TreeHandler

SONAME=$(addsuffix .so,${BASE}) 
OBJ=$(addsuffix .o,${BASE}) 
DICTOBJ=$(addsuffix Dict.o,${BASE}) 
SRC=$(addsuffix .cxx,${BASE}) 
DICTSRC=$(addsuffix Dict.cxx,${BASE}) 
HDR=$(addsuffix .h,${BASE}) 
DICTHDR=$(addsuffix Dict.h,${BASE}) 

$(SONAME) : $(OBJ)  
	g++ $(ROOTLIBS) -shared -o $(SONAME) $(OBJ) $(DICTOBJ)
	cp $(HDR) $(BASEDIR)/include 
	cp $(SONAME) $(BASEDIR)/lib 
$(OBJ) : $(DICTSRC) $(SRC)
	g++ $(CXXFLAGS) -c $(DICTSRC) $(SRC)
$(DICTSRC) : $(HDR)
	rootcint -f $(DICTSRC) -c $(HDR) LinkDef.h

clean :
	rm -f *~ *.o *.so ${DICTSRC} ${DICTHDR} ${BASEDIR}/lib/${SONAME} ${BASEDIR}/include/${HDR}
