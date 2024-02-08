#!/bin/make -f
# disable implicit rules and variables
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# name of your program
PROG := scf

SRCS := main.f90

# sources provided by your lab assistents
SRCS += integrals.f90\
        slater.f90\
        print_matrix.f90\
        input_reader.f90\
        prog.f90\
        tools.f90\
        linear_algebra.f90\
        io_tools.f90\
        input_reader.f90

# make object files from source names
OBJS := $(patsubst %, build/%.o, $(SRCS))

# configuration
FC := gfortran
LD := $(FC)
RM := rm -f
BUILD := build
vpath %.f90 src:app
# use -O0 for debugging, -O2 is fine for production runs
FFLAGS := -O
FFLAGS += -ffree-form
# the nasty GNU compiler used to cut off lines, this is over now
FFLAGS += -ffree-line-length-none
# please, write Fortran2008, we have left the era of F77 already
FFLAGS += -pedantic
FFLAGS += -std=f2008
# put all the modules in a directory
FFLAGS += -J$(BUILD)

# local Lapack installation present?
HAS_LAPACK := no
# in case your system has no LAPACK installation present, we provide
# a minimal LAPACK implementation that will be used instead

.PHONY: all
all: setup $(PROG)

# libraries you want to use
ifeq ($(HAS_LAPACK),yes)
LIBS := -llapack -lblas
else
LAPACK_OBJS := $(BUILD)/lapack.f90.o
$(BUILD)/lapack.f90.o: lapack.f90
	$(FC) -O2 -c $< -o $@
endif

# this defines how to build your program
$(PROG): $(OBJS) $(LAPACK_OBJS)
	$(FC) $^ $(LIBS) -o $@

$(OBJS): $(BUILD)/%.o: %
	$(FC) $(FFLAGS) -c $< -o $@

# dependencies must be defined explicitly
$(BUILD)/main.f90.o:    $(BUILD)/integrals.f90.o\
                        $(BUILD)/slater.f90.o\
                        $(BUILD)/print_matrix.f90.o\
                        $(BUILD)/input_reader.f90.o\
                        $(BUILD)/prog.f90.o\
                        $(BUILD)/tools.f90.o\
                        $(BUILD)/linear_algebra.f90.o\
                        $(BUILD)/io_tools.f90.o\
                        $(BUILD)/input_reader.f90.o

.PHONY: setup clean veryclean
# make a build directory
setup: $(BUILD) $(BUILD)/lib
$(BUILD) $(BUILD)/lib:
	mkdir -p $@
# remove all object files generated by this Makfile
# filter first to avoid cleaning up (=deleting) student codez.
clean:
	$(RM) $(filter %.o, $(OBJS) $(LAPACK_OBJS))

# clean up, rigorously remove *all* object files, and modules, and binaries
veryclean:
	$(RM) $(wildcard $(BUILD)/*.o) $(wildcard $(BUILD)/lib/*.o) $(wildcard $(BUILD)/*.mod) $(PROG)