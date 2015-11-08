
GKK_Main: GKK_Main.a
GKK_Main_Server: GKK_Main_Server.a
GKK_Opt_Taxes: GKK_Opt_Taxes.a
GKK_Wealth_Tax: Sergio.a
simple: Sergio_Simple.a
simple_2: Sergio_Simple_2.a

NRTYPE.o: NRTYPE.F90
	gfortran -c NRTYPE.F90

NRUTIL.o: NRUTIL.F90 NRTYPE.o
	gfortran -c NRUTIL.F90 NRTYPE.o	

Toolbox.o: Toolbox.f90 NRTYPE.o NRUTIL.o
	gfortran -c Toolbox.f90 NRTYPE.o NRUTIL.o

parameters.o: parameters.f90 NRTYPE.o
	gfortran -c parameters.f90 NRTYPE.o

global.o: global.f90 parameters.o
	gfortran -c global.f90 parameters.o

programfunctions.o: programfunctions.f90 parameters.o global.o Toolbox.o
	gfortran -c programfunctions.f90 parameters.o global.o Toolbox.o

Sergio_Simple_2.a: GKK_simple_V2.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple_V2.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple_2.a
	./Sergio_Simple_2.a

Sergio_Simple.a: GKK_simple.f95 NRTYPE.o NRUTIL.o
	gfortran GKK_simple.f95 NRUTIL.o NRTYPE.o -o Sergio_Simple.a
	./Sergio_Simple.a

Sergio.a: GKK_Wealth_Tax_Sergio.f95 NRTYPE.o NRUTIL.o Toolbox.o
	gfortran GKK_Wealth_Tax_Sergio.f95 NRUTIL.o NRTYPE.o Toolbox.o -o Sergio.a
	./Sergio.a

GKK_Main.a: GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Main.a
	./GKK_Main.a

GKK_Main_Server.a: GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Main.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Main_Server.a

GKK_Opt_Taxes.a: GKK_Optimal_Taxes.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o
	gfortran GKK_Optimal_Taxes.f95 NRUTIL.o NRTYPE.o Toolbox.o parameters.o global.o programfunctions.o -o GKK_Opt_Taxes.a
	./GKK_Opt_Taxes.a