# MAKEFILE FOR simple C++ programming
# by David MacKay

# if you want the debugger kdbg to work nicely, REMOVE the -O2 flag
# if you want to get all warnings enabled, INCLUDE the -O2 flag

CFLAGS = $(INCDIRS)  \
	-ansi -pedantic -g -O2\
	-Wall -Wconversion\
	-Wformat  -Wshadow\
	-Wpointer-arith -Wcast-qual -Wwrite-strings\
	-D__USE_FIXED_PROTOTYPES__

LIBS = -l stdc++ -lm

CXX = g++

.PHONY: clean
clean:
	rm *.o

%: %.cc
	$(CXX) $(CFLAGS) $(LIBS) $< -o $@

