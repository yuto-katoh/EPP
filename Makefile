FC = gfortran
CC = /usr/local/Cellar/gcc/11.1.0_1/bin/gcc-11

FFLAGS = -m64
CFLAGS =

PROGRAM = run.exe
EXE0 = p__ehy.o e__ehy.o
PRM0 = c__eprm.o v__in.o
OBJ0 = p__initmg.o p__initcs.o p__init.o p__passv.o p__passx.o
RAND0 = mt19937arm.o

all : $(PROGRAM)

$(PRM0): %.o : %.F90
	$(FC) $(FFLAGS) -c $<
$(RAND0): %.o : %.c
	$(CC) $(CFLAGS) -c $<
$(OBJ0): %.o : %.F90
	$(FC) $(FFLAGS) -c $<
$(EXE0): %.o : %.F90
	$(FC) $(FFLAGS) -c $<

$(PROGRAM): $(PRM0) $(RAND0) $(OBJ0) $(EXE0)
	$(FC) -o $@ $(FFLAGS) $(EXE0) $(PRM0) $(RAND0) $(OBJ0)

clean:
	rm -f *.o *.mod
