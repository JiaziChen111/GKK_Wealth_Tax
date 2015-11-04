MODULE parameters

    use nrtype

    ! Huge number
    real(dp), parameter         ::  big_p   = HUGE(1.0_dp)

    ! Switch for computing social security with benchmark earnings
    	! If KeepSSatBench=1 then E_bar is kept at E_bar_bench for experiments
    INTEGER(I4B),  PARAMETER :: KeepSSatBench=1

    ! Switch for seprable vs non-separable utiltiy
    	! If Utility_Type=1 then use separable utility
    INTEGER(I4B),  PARAMETER :: Utility_Type=0

	! Labor efficiency shocks
		! log(y)=  lambda + kappa + e 
		! lambda: inidividual fixed effect (fixed within generation)
		! kappa: life-cycle component
		! e: indiveidual transitory component (varies across time)
    ! Permanent labor earnings componenet (lambda)
    REAL(DP), PARAMETER      :: rho_lambda=0.5_DP 
    REAL(DP)             	 :: sigma_lambda_eps
    INTEGER(I4B),  PARAMETER :: nlambda=5         ! Number of grid points

    ! Transitory labor earnings componenet (e)
    REAL(DP), PARAMETER  	 :: rho_e=0.9_DP, sigma_e_eps=0.20_DP
    INTEGER(I4B),  PARAMETER :: ne=5              ! Number of grid points

    ! Entrepreneurial ability (z)
    REAL(DP)         	     :: rho_z, sigma_z_eps, mu_z
    INTEGER(I4B),  PARAMETER :: nz=7              ! Number of grid points

 

    ! Utility: Discount factor (beta) and Disutility from labor (phi) 
	REAL(DP)                 :: beta, phi, sigma, gamma 
    
	! Production 
		! Final good producer
		REAL(DP), PARAMETER  :: alpha=0.33_DP, Aprod=1.0_DP
		! Intermediate good (or home production)
		REAL(DP), PARAMETER  :: mu=0.9_DP
		! Depreciation rate
		REAL(DP), PARAMETER  :: DepRate=0.0_DP
	

	! Life cycle: retirement age, maximum age
	INTEGER(I4B), PARAMETER  :: MaxAge=81, RetAge=45


	! Asset grid: nodes (na,fine_na), min (amin), max (amax), curvature (a_curv) 
	INTEGER(I4B), PARAMETER  :: na=201, fine_na=801
	REAL(DP)    , PARAMETER  :: a_theta=4.0_DP , amax=100000.0_DP, amin=0.0001_DP

	

	! Control for updates on stationary distribution
		! Every "update_period" iterations policy functions are updated
	INTEGER(I4B), PARAMETER :: update_period=5
		! The distribution is iterated until convergence or until "MaxSimuTime" iterations
	INTEGER(I4B), PARAMETER :: MaxSimuTime=500 
		! Age categories are established
	INTEGER(I4B), PARAMETER :: max_age_category=7

	! Parameters for external functions and procedures
		! Number of std std away from mean in tauchen
		REAL(DP), PARAMETER  :: mtauchen=3.0_DP
		! Tolerance for brent algorithm
		REAL(DP), PARAMETER  :: brent_tol=0.00000001_DP
		

	! Taxes
		! Wealth tax: minimum wealth tax to consider and increments for balancing budget
		REAL(DP), PARAMETER  :: tauWmin_bt=0.00_DP, tauWinc_bt=0.000_DP ! Minimum tax below threshold and increments
		REAL(DP), PARAMETER  :: tauWmin_at=0.01_DP, tauWinc_at=0.002_DP ! Minimum tax above threshold and increments
		REAL(DP), PARAMETER  :: Threshold_Factor = 0.00_dp 
		! Consumption tax
		REAL(DP), PARAMETER  :: tauC=0.075_DP
		! Labor income tax: This is a progresive tax.
			! 1-psi controls the level of tax, and tauPL controls progressivity
		REAL(DP), PARAMETER  :: tauPL=0.185_DP, psi_PL=0.77_DP  
		!REAL(DP), PARAMETER  :: tauPL=0.0_DP, psi_PL=0.776_DP  

END MODULE parameters