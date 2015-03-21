# PreDef = -DGENERATIONMIXING  -DONLYDOUBLE -DSEESAWIII
PreDef = -DGENERATIONMIXING  -DONLYDOUBLE
# setting various paths
InDir = ../include
Mdir = ${InDir}
name = ../lib/libSPheno.a

#
# options for various compilers
#

# Intels ifort, default in optimized mode
F90 = ifort
comp = -c -O -module ${Mdir} -I${InDir} 
LFlagsB = -O  

# Intels ifort, debug modus
ifeq (${F90},ifortg)
 F90 = ifort
 comp = -c -g -module ${Mdir} -I${InDir} 
 LFlagsB = -g  
endif

# gfortran
ifeq (${F90},gfortran)
 comp = -c -O -J${Mdir} -I${InDir}
 LFlagsB = -O  
endif

# g95 
ifeq (${F90},g95)
 comp = -c -O -fmod=${Mdir} -I${InDir}
 LFlagsB = -O  
endif

# Lahey F95 compiler
ifeq (${F90},lf95)
 comp = -c -O -M ${Mdir} -I${InDir}
 LFlagsB = -O  
endif
 
# NAG f95/2003
ifeq (${F90},nagfor)
 comp = -c -O  -DONLYDOUBLE -mdir ${Mdir} -I${InDir}   
 LFlagsB = -O
endif
 
.SUFFIXES : .o .ps .f90 .F90 .a
bin/SPheno: ${name} SPheno3.o
	${F90} -o SPheno ${LFlagsB} SPheno3.o ../lib/${name}
	mv SPheno ../bin
${name}: ${name}(Control.o)  ${name}(Mathematics.o)  ${name}(RGEs.o)       \
    ${name}(MathematicsQP.o)  ${name}(LoopFunctions.o) ${name}(StandardModel.o) \
  ${name}(Model_Data.o) ${name}(Couplings.o)  ${name}(SusyMasses.o) \
   ${name}(LoopCouplings.o)  ${name}(DecayFunctions.o)  \
  ${name}(SusyDecays.o) ${name}(ThreeBodyPhaseSpace.o)               \
  ${name}(ThreeBodyPhaseSpaceS.o) ${name}(Chargino3.o)               \
  ${name}(Gluino3.o) ${name}(Neutralino3.o)                          \
  ${name}(Stop3BodyDecays.o) ${name}(Slepton3Body.o) ${name}(BranchingRatios.o)              \
  ${name}(EplusEminusProduction.o) ${name}(TwoLoopHiggsMass.o) \
  ${name}(LoopMasses.o) ${name}(SugraRuns.o) ${name}(Experiment.o) \
  ${name}(LowEnergy.o) ${name}(NMSSM_tools.o) ${name}(RPtools.o) \
  ${name}(LHC_observables.o) ${name}(InputOutput.o) 
clean:
	rm -f *.o *~ */*.o */*~
cleanall:
	rm -f bin/SPheno3 lib/*.a *~ */*.o */*~ include/*
#
# Suffix rules
#
.f90.a:
	${F90} ${comp} $<
	ar -ruc $@ $*.o
	rm -f $*.o
.F90.a:
	${F90} ${comp}  ${PreDef} $<
	ar -ruc $@ $*.o
	rm -f $*.o
.f90.o: 
	${F90} ${comp}  $< 
.f90.ps:
	a2ps  -o $*.ps $<
.h.ps:
	a2ps -o  $*.ps $<

