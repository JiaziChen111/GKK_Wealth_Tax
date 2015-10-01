
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a
simple_2: Sergio_Simple_2.a

NRTYPE.o: NRTYPE.f90
	gfortran -c NRTYPE.f90

NRUTIL.o: NRUTIL.f90 NRTYPE.o
	gfortran -c NRUTIL.f90 NRTYPE.o

Toolbox.o: Toolbox.f90 NRTYPE.o NRUTIL.o
	gfortran Toolbox.f90 NRTYPE.o NRUTIL.o

Sergio_Simple_2.a: GKK_simple_V2.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple_V2.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple_2.a
	./Sergio_Simple_2.a

Sergio_Simple.a: GKK_simple.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple.a
	./Sergio_Simple.a

Sergio.a: GKK_Wealth_Tax_Sergio.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_Wealth_Tax_Sergio.f95 NRUTIL.o NRTYPE.o -o Sergio.a
	./Sergio.a