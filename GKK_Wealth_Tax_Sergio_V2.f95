

! This program solves for the wealth taxation model of 
! Guvenen, Kuruscu and Kambourov

! The code allows for:
	! Progressive labor income taxation
	! Threshold on wealth for wealth taxation
	! Non-Separable utility

! When utility is separable it is:
	! U(c,h) = log(c) + phi*log(1-h)

! When utility is non-separable it is:
	! U(c,h) = (c^(gamma)(1-l)^(1-gamma))^(1-sigma) / (1-sigma)


!========================================================================================
!========================================================================================
!========================================================================================




PROGRAM main
	USE parameters
	USE GLOBAL
	use programfunctions
	use Toolbox

	IMPLICIT NONE
	! Variables to measure running time
		REAL(DP) :: start_time, finish_time
	! Values of parameters for different parametrizations (Low and High)
		REAL(DP) :: betaL, betaH, phiL, phiH, sigmalambL,sigmalambH, sigmazL, sigmazH, rhozL, rhozH
	! Number of different values to be considered for each parameter 
		INTEGER  ::nbeta, nphi, nsigmalambda, nrhoz, nsigmaz
	! Counters for runs of the model
		INTEGER  :: parindx1,  parindx2, parindx3, parindx4, parindx5
	! Compute benchmark or load results
		INTEGER  :: read_write_bench

	! Resutls Folder
		write(Result_Folder,'(f4.2)') Threshold_Factor
		Result_Folder = './V2_Results/Factor_'//trim(Result_Folder)//'/'
		! call execute_command_line( 'mkdir -p ' // trim(Result_Folder) )
		call system( 'mkdir -p ' // trim(Result_Folder) )
		print*, "Results are stored in directory: ", Result_Folder
	
	! Unused values of parameters
		! the following solves equilibrium of capital tax economy
! 		params=[  0.9500,    0.7000,    0.2996,    0.5673,    1.2280]
! 		params=[  0.9510,    0.5000,    0.4060,    0.5680,    1.2250]
! 		params=[  0.9503,    0.6500,    0.3288,    0.5667,    1.2280]
! 		params=[  0.9510,    0.5250,    0.3942,    0.5680,    1.2247]
! 		params=[  0.9511,    0.4200,    0.4400,    0.5680,    1.2250]
! 		params=[  0.9522,    0.3400,    0.4740,    0.5690,    1.2240]
! 		params=[  0.9506,    0.6000,    0.3564,    0.5667,    1.2280]

! 		! NEW PARAMETERS

! 		params=[ 0.947  ,  0.4 , 0.490 , 0.340 , 1.01 ]
! 		params=[ 0.9455 ,  0.6 , 0.381 , 0.335 , 1.00 ]
! 		params=[ 0.9455 ,  0.8 , 0.255 , 0.34  , 1.00 ]
! 		params=[ 0.948  ,  0.2 , 0.56  , 0.34  , 1.02 ]

! 		! Latests parameters
! 		params=[0.9489479959, 0.40, 0.514595031738281, 0.34, 1.01]


	!print*,'---------------------------       PSI    NOT  ADJUSTED   ---------------------------'
	!print*,'------------------------- RETIREMENT BENEFITS ADJUSTED - DBN ADJUSTED ------------------------'
	!print*,'------------------------- CONS TAX SIMPLE ------------------------'
	print*,'na=',na,'update_period=',update_period


		
! 		beta             = params(1)
! 		rho_z            = params(2)
! 		sigma_z_eps      = params(3)
! 		sigma_lambda_eps = params(4)
! 		phi              = params(5)
! 		!print*,  beta, rho_z, sigma_z_eps, sigma_lambda_eps,  phi

! 		sigma            = 4.0_dp
! 		gamma            = 1.0_dp/(1.0_dp+phi)


		Params =[ 0.9436, 0.00, 0.50, 0.70444445, 0.34, 0.4494 ] ! tauL=0.224, tauC=0.075 calibration
		beta   = params(1)
		mu_z   = params(2) ! this is just shifting the z grids. it is zero now.
		rho_z  = params(3) 
		sigma_z_eps      =params(4)
		sigma_lambda_eps = params(5)
		gamma  = params(6)
		sigma  = 4.0_dp

	! Set parameters to be used in all simulations economy
		OPEN   (UNIT=3, FILE=trim(Result_Folder)//'params', STATUS='replace')
		WRITE(unit=3, FMT=*) params
		CLOSE (unit=3)


	! Start timining of  the process
		call cpu_time(start_time) 

		
	! Set initia lvalues of R, Wage, Ebar to find equilibrium
		! ------- DO NOT REMOVE THE LINES BELOW

		rr    =  4.906133597851297E-002 
		!! rr    =  0.10_dp
		wage  =  1.97429920063330 
		Ebar  =  1.82928004963637  
		Ebar_bench = Ebar
		 

		! ------- DO NOT REMOVE THE LINES ABOVE

	!--------------------------------------------------------------------------------------------------------------------------

		!call cpu_time(start_time) 
		! 
		!
		!
		!betaL =0.947
		!betaH=0.948
		!rhozL = 0.4
		!rhozH = 0.4
		!sigmazL = 0.48
		!sigmazH = 0.49
		!sigmalambL = 0.34
		!sigmalambH =0.35
		!phiL    = 0.99
		!phiH   = 1.01
		!
		!
		!!betaL =0.9455
		!!betaH=0.946
		!!rhozL = 0.6
		!!rhozH = 0.6
		!!sigmazL = 0.381
		!!sigmazH = 0.385
		!!sigmalambL = 0.335
		!!sigmalambH =0.34
		!!phiL    = 0.99
		!!phiH   = 1.00
		!
		!
		!betaL =0.9455
		!betaH=0.946
		!rhozL = 0.8
		!rhozH = 0.8
		!sigmazL = 0.253
		!sigmazH = 0.255
		!sigmalambL = 0.335
		!sigmalambH =0.34
		!phiL    = 1.00
		!phiH   = 1.02
		!
		!betaL =0.948
		!betaH=0.949
		!rhozL = 0.2
		!rhozH = 0.2
		!sigmazL = 0.56
		!sigmazH = 0.57
		!sigmalambL = 0.34
		!sigmalambH =0.35
		!phiL    = 1.02
		!phiH   = 1.04
		!
		!
		!nbeta =2
		!nrhoz=1
		!nsigmaz=2
		!nsigmalambda=2
		!nphi=2
		!
		!Min_SSE_Moments=1000.0_DP
		!
		!DO parindx3=1,nsigmalambda
		!DO parindx2=1,nphi
		!DO parindx1=1,nbeta
		!DO parindx4=1,nrhoz
		!DO parindx5=1,nsigmaz
		!
		!    beta = betaL + real(parindx1-1,8) *(betaH-betaL)/max(real(nbeta-1,8),1.0_DP)
		!    phi   = phiL   + real(parindx2-1,8) *(phiH-phiL)/max(real(nphi-1,8),1.0_DP)
		!    sigma_lambda_eps = sigmalambL + real(parindx3-1,8)*(sigmalambH -sigmalambL) / max(real(nsigmalambda-1,8),1.0_DP)
		!    rho_z= rhozL   +  real(parindx4-1,8)*(rhozH-rhozL) / max(real(nrhoz-1,8),1.0_DP)
		!    sigma_z_eps = sigmazL +  real(parindx5-1,8)*(sigmazH-sigmazL) / max(real(nsigmaz-1,8),1.0_DP)
		!
		!    CALL  INITIALIZE
		!    CALL FIND_DBN_EQ
		!    CALL COMPUTE_STATS
		!    IF (SSE_Moments .lt. Min_SSE_Moments ) THEN
		!        Min_SSE_Moments =SSE_Moments
		!        params= [ beta, rho_z, sigma_z_eps, sigma_lambda_eps, phi ]
		!        Min_Moments = [  Wealth_Output, prct1_wealth , prct10_wealth, Std_Log_Earnings_25_60, meanhours_25_60  ]
		!    ENDIF
		!    !CALL WRITE_TO_FILE
		!
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!ENDDO
		!print*, params
		!print*, Min_Moments
	!

	PRINT*,''
	Print*,'--------------- SOLVING BENCHMARK WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'CAPITAL TAX ECONOMY'

	! Benchmark economy
		solving_bench=1

	! Set taxes for benchmark economy
		tauK = 0.25_DP
		tauL = 0.30_DP
		tauW_bt = 0.00_DP
		tauW_at = 0.00_DP
		Y_a_threshold = 0.00_DP 

	! Solve for the model and compute stats
	read_write_bench = 0
	print*,"	Initializing program"
		CALL INITIALIZE
	if (read_write_bench.eq.0) then
		print*,"	Computing equilibrium distribution"
		CALL FIND_DBN_EQ
		print*,"	Computing satitics"
		CALL COMPUTE_STATS
		print*,"	Computing government spending"
		CALL GOVNT_BUDGET
		print*,"	Writing variables"
		CALL WRITE_VARIABLES(1)
		print*,"	Computing Value Function"
		CALL COMPUTE_VALUE_FUNCTION_SPLINE 
		print*,"	Saving results in text files to be read later"
		CALL Write_Benchmark_Results(read_write_bench)
	else
		print*,"	Reading benchmark results from files"
		CALL Write_Benchmark_Results(read_write_bench)
	end if 

	! Aggregate variables in benchmark economy
		GBAR_bench  = GBAR
		QBAR_bench  = QBAR 
		NBAR_bench  = NBAR 
		Ebar_bench  = EBAR
		rr_bench    = rr
		wage_bench  = wage
		Y_bench     = YBAR
		tauK_bench  = tauK
		tauL_bench  = tauL
		DBN_bench   = DBN1
		tauw_bt_bench = tauW_bt
		tauw_at_bench = tauW_at
		Y_a_threshold_bench = Y_a_threshold

		write(*,*) "Benchmark variables"
		write(*,*) "GBAR=",GBAR,"EBAR=",EBAR,"NBAR=",NBAR,"QBAR=",QBAR,"rr=",rr,"wage=",wage

	PRINT*,''
	Print*,'--------------- SOLVING EXPERIMENT WITH BEST PARAMETERS -----------------'
	PRINT*,''
	print*,'Wealth Tax Economy'
	
	! Experiment economy
		solving_bench=0
	! Set capital taxes to zero
		tauK = 0.0_DP
	! Set Y_a_threshold
		write(*,*) "Y_a threshold is set to a proportion of the mean wealth under current distribution"
		!Y_a_threshold = 0.0_dp ! big_p   !8.1812138704441200
		call Find_TauW_Threshold(DBN_bench,W_bench)  
		Y_a_threshold = Threshold_Factor*Ebar_bench !0.75_dp
		Wealth_factor = Y_a_threshold/W_bench

	! Find wealth taxes that balances budget
	print*, "	Computing Wealth Tax to balance the budget"
		! Set initial value for G in experimental economy and for wealth taxes
		GBAR_exp = 0.0_DP
		tauW_bt  = tauWmin_bt
		tauW_at  = tauWmin_at
		tauWindx = 0.0_DP
		! Solve for the model increasing wealth taxes until revenue is enough to finance G_benchamark
		DO WHILE (GBAR_exp .lt. GBAR_bench)
			! Set old G and new value of tauW
			GBAR_exp_old = GBAR_exp
			tauW_bt = tauWmin_bt + tauWindx * tauWinc_bt
			tauW_at = tauWmin_at + tauWindx * tauWinc_at
			! Solve the model
			CALL FIND_DBN_EQ
			CALL GOVNT_BUDGET

			! Get new G
			GBAR_exp = GBAR 
			! Iteratioins  
			tauWindx = tauWindx + 1.0_DP   
			write(*,*) "Bracketing GBAR: tauW_bt=", tauW_bt*100, "And tauW_at=", tauW_at*100
			print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
		ENDDO

		! Set tauW as weighted average of point in  the grid to balance budget more precisely
			tauW_up_bt  = tauW_bt
			tauW_low_bt = tauW_bt  -  tauWinc_bt
			tauW_bt     = tauW_low_bt + tauWinc_bt * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
			tauW_up_at  = tauW_at
			tauW_low_at = tauW_at  -  tauWinc_at  
			tauW_at     = tauW_low_at + tauWinc_at * (GBAR_bench - GBAR_exp_old )/(GBAR_exp - GBAR_exp_old)
			print*,''
			print*,'GBAR bracketed by taxes:'
			print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
			print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
			print*,''

		! Solve (again) experimental economy
			CALL FIND_DBN_EQ
			CALL GOVNT_BUDGET

		! Find tauW that exactly balances the budget (up to precisioin 0.1) using bisection
			GBAR_exp = GBAR
			print*,"Gbar at midpoint of bracket and GBAR at benchmark"
			print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			print*,''
			print*,'Bisection for TauW:'
			DO WHILE (  abs(100.0_DP*(1.0_DP-GBAR_exp/GBAR_bench)) .gt. 0.1 ) ! as long as the difference is greater than 0.1% continue
			    if (GBAR_exp .gt. GBAR_bench ) then
			        tauW_up_bt  = tauW_bt 
			        tauW_up_at  = tauW_at 
			    else
			        tauW_low_bt = tauW_bt
			        tauW_low_at = tauW_at
			    endif
			    tauW_bt = (tauW_low_bt + tauW_up_bt)/2.0_DP
			    tauW_at = (tauW_low_at + tauW_up_at)/2.0_DP
			    CALL FIND_DBN_EQ
			    CALL GOVNT_BUDGET
			    GBAR_exp = GBAR
			    print*,'tauW_low_bt =', tauW_low_bt*100, '% tauW_up_bt=', tauW_up_bt*100, '% tauW_bt=', tauW_bt*100, "%"
				print*,'tauW_low_at =', tauW_low_at*100, '% tauW_up_at=', tauW_up_at*100, '% tauW_at=', tauW_at*100, "%"
				print*, "Current Threshold for wealth taxes", Y_a_threshold, "Share above threshold=", Threshold_Share
				print*,'GBAR_exp =', GBAR_exp,'GBAR_bench=',GBAR_bench
			ENDDO

	! AGgregate variable in experimental economy
		GBAR_exp  = GBAR
		QBAR_exp  = QBAR 
		NBAR_exp  = NBAR  
		Y_exp 	  = YBAR
		Ebar_exp  = EBAR
		rr_exp    = rr
		wage_exp  = wage
		tauK_exp  = tauK
		tauL_exp  = tauL
		DBN_exp   = DBN1
		tauw_bt_exp = tauW_bt
		tauw_at_exp = tauW_at
		Y_a_threshold_exp = Y_a_threshold


	!!!	CALL COMPUTE_VALUE_FUNCTION_SPLINE 
	!!! CALL Write_Experimental_Results()
	CALL COMPUTE_STATS
	CALL WRITE_VARIABLES(0)
	! Compute welfare gain between economies
		CALL COMPUTE_WELFARE_GAIN

	! Write in files some stats
		OPEN (UNIT=19, FILE=trim(Result_Folder)//'Stats_Resuls', STATUS='replace') 
			WRITE(UNIT=19, FMT=*) "Threshold_Factor"     	, Threshold_Factor
			WRITE(UNIT=19, FMT=*) "Wealth_Factor"		  	, Wealth_Factor
			WRITE(UNIT=19, FMT=*) "Threshold"			  	, Y_a_Threshold
			WRITE(UNIT=19, FMT=*) "Wealth_Tax_Above"	    , TauW_at
			WRITE(UNIT=19, FMT=*) "Welfare_Gain_Pop(bench)" , Welfare_Gain_Pop_bench
			WRITE(UNIT=19, FMT=*) "Welfare_Gain_Pop(exp)"   , Welfare_Gain_Pop_exp
			WRITE(UNIT=19, FMT=*) "Welfare_Gain_NB(bench)"  , Welfare_Gain_NB_bench
			WRITE(UNIT=19, FMT=*) "Welfare_Gain_NB(exp)"    , Welfare_Gain_NB_exp
			WRITE(UNIT=19, FMT=*) "Output_Gain(prct)"	  	, 100.0_DP*(Y_exp/Y_bench-1.0) 
			WRITE(UNIT=19, FMT=*) "W/GDP"				  	, Wealth_Output
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_1%' 	, prct1_wealth
			WRITE(UNIT=19, FMT=*) 'Wealth_held_by_Top_10%'	, prct10_wealth
			WRITE(UNIT=19, FMT=*) 'STD_Labor_Earnings'	  	, Std_Log_Earnings_25_60
			WRITE(UNIT=19, FMT=*) 'Mean_Labor_Earnings'   	, meanhours_25_60
			WRITE(UNIT=19, FMT=*) 'Moments'				  	, SSE_Moments 
			WRITE(UNIT=19, FMT=*) 'Share_Above_Threshold'	, Threshold_Share
			WRITE(UNIT=19, FMT=*) ' '
			WRITE(UNIT=19, FMT=*) 'GBAR_bench'				, GBAR_bench                   , 'GBAR_exp=' , GBAR_exp
			WRITE(UNIT=19, FMT=*) 'QBAR_bench'				, QBAR_bench                   , 'QBAR_exp=' , QBAR_exp
			WRITE(UNIT=19, FMT=*) 'NBAR_bench'				, NBAR_bench                   , 'NBAR_exp=' , NBAR_exp
			WRITE(UNIT=19, FMT=*) 'YBAR_bench'				, Y_bench                  	   , 'YBAR_exp=' , Y_exp
			WRITE(UNIT=19, FMT=*) 'EBAR_bench'				, EBAR_bench                   , 'EBAR_exp=' , EBAR_exp
			WRITE(UNIT=19, FMT=*) 'rr_bench'				, rr_bench                     , 'rr_exp='   , rr_exp
			WRITE(UNIT=19, FMT=*) 'wage_bench'				, wage_bench                   , 'wage_exp=' , wage_exp
		CLOSE(Unit=19)
	
	print*,'---------------------------'
	print*,''
	print*,'Output Gain Prct=', 100.0_DP*(Y_exp/Y_bench-1.0) 
	print*,''
	print*,'---------------------------'

	print*," "
	print*,"Wealth_factor=",Wealth_factor
	print*," "

	call cpu_time(finish_time)
	print*,'Total time =',finish_time-start_time

! 	print*, " "
! 	print*, "Asset Function"
! 	do ai=1,na 
! 		print*,agrid(ai), Aprime(20,ai,:,3,3)
! 	end do 

END PROGRAM main

!========================================================================================
!========================================================================================
!========================================================================================

