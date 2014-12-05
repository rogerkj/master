# Standard makefile

CXX= c++
CPPFLAGS = -O3 -llapack -lblas -larmadillo
PROG = main
LIB1 = lib
LIB2 = mcintegrator
LIB3 = LocalEnergy
WAVE = WaveFunction
LEN = LocalEnergyNumeric
HYD = Hydrogenic
LEO = LocalEnergyOpt
QF = QuantumForce
QFN = QuantumForceNumeric
DI = Diatomic
ORB = Orbitals
GAUSS = Gaussians

${PROG} : 	${PROG}.o ${LIB1}.o  ${LIB2}.o ${LIB3}.o ${WAVE}.o ${LEN}.o ${HYD}.o ${LEO}.o ${QF}.o ${QFN}.o ${DI}.o ${ORB}.o ${GAUSS}.o
		${CXX} ${PROG}.o ${LIB1}.o ${LIB2}.o ${LIB3}.o ${WAVE}.o ${LEN}.o ${HYD}.o ${LEO}.o ${QF}.o ${QFN}.o ${ORB}.o ${DI}.o ${GAUSS}.o -o ${PROG} ${CPPFLAGS}


clean:
	rm -f ${PROG} ${PROG}.o
	rm -f ${LIB1} ${LIB1}.o
	rm -f ${LIB2} ${LIB2}.o
	rm -f ${LIB3} ${LIB3}.o
	rm -f ${WAVE} ${WAVE}.o
	rm -f ${LEN} ${LEN}.o
	rm -f ${HYD} ${HYD}.o
	rm -f ${LEO} ${LEO}.o
	rm -f ${QF} ${QF}.o
	rm -f ${QFN} ${QFN}.o
	rm -f ${DI} ${DI}.o
	rm -f ${ORB} ${ORB}.o
	rm -f ${GAUSS} ${GAUSS}.o
