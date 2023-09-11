# C source code
CSRC	= io_png/io_png.c \
			mt19937ar.c

# C++ source code
CXXSRC	= flutter.cpp \
	  borders.cpp\
	  fftw_routines.cpp\
	  standard_routines.cpp\
	  midway.cpp\
	  codes_flutter.cpp\
		demo_flutter.cpp


# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= demo_flutter

default	: $(BIN)

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	= -lpng -lm -lfftw3


# use DEBUG
ifdef DEBUG
CFLAGS	+= -g
CXXFLAGS	+= -g
LDFLAGS += -g
endif

# build the local png library
.PHONY	: libpng
libpng	:
	$(MAKE) -C io_png/libs libpng

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the opject code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
	$(MAKE) -C ./io_png/libs $@
distclean	: clean
	$(RM) $(BIN)
	$(MAKE) -C ./io_png/libs $@
