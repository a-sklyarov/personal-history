
CC		= g++ -O2 -g -pedantic -Wall
INCLUDE		= -I$(CAV_INC) -I$(GSL_INC) -I$(FFTW_INC) 
LIBS		= $(CAV_LIB) $(GSL_LIB) $(FFTW_LIB)
RM		= /bin/rm -f

SOURCES_1	= CoreTask1.cc aperture.cc screen.cc
SOURCES_2	= CoreTask2.cc aperture.cc screen.cc
SOURCES_3	= SupplTask1a.cc aperture.cc screen.cc
SOURCES_4	= SupplTask1b.cc aperture.cc screen.cc

OBJECTS_1	= $(SOURCES_1:.cc=.o)
OBJECTS_2	= $(SOURCES_2:.cc=.o)
OBJECTS_3	= $(SOURCES_3:.cc=.o)
OBJECTS_4	= $(SOURCES_4:.cc=.o)

EXECUTABLE_1 	= CoreTask1.exe
EXECUTABLE_2	= CoreTask2.exe
EXECUTABLE_3 	= SupplTask1a.exe
EXECUTABLE_4	= SupplTask1b.exe

all: CoreTask1 CoreTask2 SupplTask1a SupplTask1b

CoreTask1: $(SOURCES_1) $(EXECUTABLE_1)
	
CoreTask2: $(SOURCES_2) $(EXECUTABLE_2)
	
SupplTask1a: $(SOURCES_3) $(EXECUTABLE_3)
	
SupplTask1b: $(SOURCES_4) $(EXECUTABLE_4)
	
$(EXECUTABLE_1): $(OBJECTS_1) 
	$(CC) $(INCLUDE) $(OBJECTS_1) -o $@ $(LIBS)

$(EXECUTABLE_2): $(OBJECTS_2) 
	$(CC) $(INCLUDE) $(OBJECTS_2) -o $@ $(LIBS)

$(EXECUTABLE_3): $(OBJECTS_3) 
	$(CC) $(INCLUDE) $(OBJECTS_3) -o $@ $(LIBS)

$(EXECUTABLE_4): $(OBJECTS_4) 
	$(CC) $(INCLUDE) $(OBJECTS_4) -o $@ $(LIBS)

%.o : %.cc
	$(CC) $(INCLUDE) -c $< -o $@

clean: 
	$(RM) *.o *.exe

