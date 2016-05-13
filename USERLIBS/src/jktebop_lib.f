!=======================================================================
!     PROGRAM JKTEBOP             John Taylor   j.k.taylor@warwick.ac.uk
!                                 Astrophysics Group  Warwick University
!-----------------------------------------------------------------------
! This is the library version of jktebop
!-----------------------------------------------------------------------
! V(1) = surface brightness ratio    V(15) = third light
! V(2) = sum of fractional radii     V(16) = phase correction
! V(3) = ratio of stellar radii      V(17) = light scaling factor
! V(4) = linear LD for star A        V(18) = integration ring size (deg)
! V(5) = linear LD for star B        V(19) = orbital period (days)
! V(6) = orbital inclination         V(20) = ephemeris timebase (days)
! V(7) = e cos(omega) OR ecentricity V(21) = nonlinear LD for star A
! V(8) = e sin(omega) OR omega       V(22) = nonlinear LD for star B
! V(9) = gravity darkening 1         VEXTRA(1) = primary star radius
! V(10)= gravity darkening 2         VEXTRA(2) = secondary star radius
! V(11) = primary reflected light    VEXTRA(3) = stellar light ratio
! V(12) = secondary reflected light  VEXTRA(4) = eccentricity
! V(13) = stellar mass ratio         VEXTRA(5) = periastron longitude
! V(14) = tidal lead/lag angle (deg) VEXTRA(6) = reduced chi-squared
!-----------------------------------------------------------------------
! Version 1: Simplex minimisation algorithm  and new input / output used
! Version 2: Monte Carlo simulation and parameter perturbation algorithm
! Version 3: Adjustments to Monte Carlo LD coeffs and input/output files
! Version 4: Now solves for sum of radii, convergence criterion modified
! Version 5: Added TASK0 to find LD and GD coeffs.   Minor modifications
! Version 6: Reflection and  scale  factor  can all be fixed or adjusted
! Version 7: Can use either (e,w) or (ecosw,esinw); SFACT modified to be
!            in magnitudes; observ'l errors found; spherical star option
! Version 8: Bootstrapping error analysis algorithm added and output mod
! Version 9: Command-line arguments allowed, Monte Carlo without parame-
!            eter kicking option, and fitting for period and Tzero added
! Version 10: Can now have 9999 datapoints. Whole code now in magnitudes
! Version 11: Bug fixes, tasks renumbered,  added sigma clipping and the
!             global fit procedures, but not thoroughly tested these yet
! Version 12: Nonlinear limb darkening law, fitting for times of minimum
!             light, and FMAX corrections included (from Alvaro Gimenez)
! Version 13: Removed  BILINEAR  and modified  TASK1  to just call JKTLD
!             Modified input file and arrays  for the diff types of data
! Version 14: Fixed the requirement for inputtd INTRING to be an integer
!             Fixed formal errors (when observational ones not supplied)
! Last modified: 10th February 2006  (working on beta Aurigae)
!-----------------------------------------------------------------------
! Possible modifications in future:
! 1) extend to WD2003 and WINK
! 4) port to F90 or F95 to use long lines, modules, improved output
! 5) JVC suggests that formal fitting errors are x2 too small
! 6) incorporate change of omega (apsidal motion)
! 8) go through subroutine LIGHT to understand it all
!-----------------------------------------------------------------------
! Miscellaneous notes:
! 1) Phase shift has been redefined compared to original EBOP so that
!    it corresponds directly to the phase of the primary minimum.
! 2) MRQMIN adjusts coeffs only if VARY (called 'ia' in MRQMIN) is 1.
! 3) If VARY=2 then the parameter is fixed during the initial fit but is
!    perturbed by a set amount (flat distribution) for later analyses.
! 4) If VARY(11) and/or  VARY(12) are -1 then V(11) and/or V(12) are
!    calculated from the system geometry; if 0 they are fixed at the
!    input value and if 1 are freely adjusted to best fit.
! 5) If the mass ratio is <= 0 then both stars are assumed to be spheres
! 6) If ecosw > 5.0 then (ecosw,esinw) will be taken to be  (10+e,omega)
!    and fitting will occur using e and omega as parameters. e and omega
!    can be strongly correlated,  but this is useful if  e  is known but
!    omnega is not;   this can happen for EBs exhibiting apsidal motion.
! 7) Observational errors are looked for in the input light curve file.
!    If they are not found then equal weight is given to each point.
! 8) Nonlinear LD is now supported for the two-coefficient logarithmic,
!    quadratic and square-root laws. The type of law must be specified.
!    on input. Star B can also be forced to the same coeffs as star A.
!    BUT: normalisation for logarithmic not possible (not got equations)
! 9) Fitting for times of minimum light is directly possible. The cycle
!    numbers and times are inputted on lines immediately below all the
!    parameter lines in the input file.
!-----------------------------------------------------------------------
! Task numbers and purposes:
! (1) This outputs LD coefficients for a given Teff, logg, [M/H], Vmicro
! (2) This outputs a model light curve for fixed input parameters.
! (3) This fits a model to an observed light curve and outputs results.
! (4) This fits a model, rejects discrepant observations, and refits.
! (5) This does a pseudo-global minimisation by perturbing input params.
! (6) This investigates how different parameters vary around best fit.
! (7) This conducts bootstrapping simulations to find robust errors.
! (8) This conducts Monte Carlo simulations to find robust errors.
!=======================================================================
! Language: JKTEBOP is written in FORTRAN 77 using standard F77 syntax
! with a few straightforward extras, and compiled using g77.
! Possible non-F77 bits: == <= < > >= /= ! and part-line comments
!=======================================================================
!=======================================================================
      SUBROUTINE TASK1 ()           ! This task outputs  limb  darkening
            ! coefficients for given Teff and logg and [M/H] and Vmicro.
            ! It simply interfaces with the  JKTLD  code, which performs
            ! bilinear interpolation in Teff,logg for given [M/H],Vmicro
            ! Usage:   jktld  <Teff>  <logg>  <M/H>  <Vmicro>  <outfile>
      implicit none
      real*8 TEFF,LOGG              ! IN: Parameters  to  interpolate to
      real*8 MOH,VMICRO             ! IN: Other parameters forthe tables
      character*20 CTEFF,CLOGG      ! LOCAL: character version of values
      character*20 CMOH,CMICRO      ! LOCAL: character version of values
      character*30 OUTFILE          ! LOCAL: name of output file to make

      write (*,'(A40,$)') "Enter the effective temperature (K)  >> "
      read (*,*) TEFF
      write (*,'(A40,$)') "Enter the surface gravity (log cm/s) >> "
      read (*,*) LOGG
      write (*,'(A40,$)') "Enter the metal abundance  ([M/H])   >> "
      read (*,*) MOH
      write (*,'(A40,$)') "Enter the microturbulence velocity   >> "
      read (*,*) VMICRO
      write (*,'(A40,$)') "Enter the output file name to create >> "
      read (*,*) OUTFILE

      if ( TEFF < 3500.0d0 ) then
        write (*,'(A39,A41)') "### Warning: a Teff below 3500 K is out",
     &                      " of range of most of the LD coeff tables."
      else if ( TEFF < 2000.0d0 ) then
        write (*,'(A39,A41)') "### Warning: a Teff below 2000 K is out",
     &                      " of range of the LD coefficient tables.  "
      else if ( TEFF > 50000.0d0 ) then
        write (*,'(A39,A41)') "### Warning: a Teff above 50000 K is ou",
     &                      "t of range of the LD coefficient tables. "
      end if
      if ( LOGG > 5.0d0 .or. LOGG < 0.0d0 )
     &  write (*,'(A39,A41)') "### Warning: log(g)a outside the range ",
     &                      "0.0 to 5.0 are not covered in the tables."
      if ( MOH /= 0.0 .and. VMICRO /= 2.0 )
     &  write (*,'(A39,A41)') "### Warning: for [M/H] /= 0 and Vmicro ",
     &                      "/= 2 you will probably get no results.   "

      write (CTEFF,'(F20.8)') TEFF
      write (CLOGG,'(F20.8)') LOGG
      write (CMOH,'(F20.8)') MOH
      write (CMICRO,'(F20.8)') VMICRO

      CALL SYSTEM ( "jktld " // CTEFF // " " //  CLOGG // " " // CMOH
     &                              // " " // CMICRO // " " // OUTFILE )

      END SUBROUTINE TASK1
!=======================================================================
      SUBROUTINE TASK2 (V,LDTYPE)         ! Produces a model light curve
      implicit none
      real*8 V(22)                        ! IN: light  curve  parameters
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      integer i,ERROR                     ! LOCAL: counters & error flag
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      real*8 PHASE                        ! LOCAL:  phase for evaluation
      real*8 GETMODEL                     ! FUNCTION: evaluate the model

      V(19) = 1.0d0           ! Set period to 1.0
      V(20) = 0.0d0           ! Set Tzero to 0.0

      MAG = GETMODEL (V,LDTYPE,0.0d0,1,LP,LS)
      V(11) = 0.4d0 * LS * (V(2)/(1.0d0+V(3)))**2
      V(12) = 0.4d0 * LP * (V(2)/(1.0d0+(1.0d0/V(3))))**2

      write (*,'(A40,A40)') ">> The reflection coefficients come from ",
     &                      "the system geometry, not the input file. "

      write (62,'(A20)') "# PHASE   MAGNITUDE "
      do i = 0,10000
        PHASE = i / 10000.0d0
        write (62,100) PHASE,GETMODEL(V,LDTYPE,PHASE,1,LP,LS)
      end do
      close (62)

100   FORMAT (F8.6,1X,F12.8)

      END SUBROUTINE TASK2
!=======================================================================
      SUBROUTINE TASK34 (TASK,V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,
     &                                                            SIGMA)
            ! Finds the best-fitting light curve params for the dataset.
            ! If TASK=4 then it rejects all datapoints which are SIGMA
            ! sigma away from the best fit and refits the remainder.
      implicit none
      integer TASK                        ! IN: which task to do
      real*8 V(22)                        ! IN: light  curve  parameters
      integer VARY(22)                    ! IN: parameters vary or fixed
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATA(3,99999)                ! IN: time, magnitude, error
      integer DTYPE(99999)                ! IN: type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: number of  types of data
      real*8 SIGMA                        ! IN:  std. dev. for  clipping
      real*8 CHISQ                        ! SUB: chi-square of model fit
      real*8 VERR(22)                     ! SUB: parameter formal errors
      real*8 OMC(99999)                   ! LOCAL: (O-C) residual values
      real*8 SIG                          ! LOCAL: rms of the O-C values
      real*8 RESIDSQSUM                   ! LOCAL: sum of resid. squares
      integer ITER,IFAIL                  ! LOCAL: iter number & success
      real*8 MAG,LP,LS                    ! LOCAL: EBOP/GETMODEL  output
      integer KEEP(99999),COUNT           ! LOCAL: Datapoint bookkeeping
      integer i                           ! LOCAL: loop counter
      real*8 GETMODEL                     ! FUNCTION: evaluate the model

      do i = 1,99999
        OMC(i) = 0.0d0
        KEEP(i) = 0
      end do
      SIG = 0.0d0

            ! Output initial values and then find the best fit.
            ! If TASK=3 then output various results and finish.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,0)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                       VERR,IFAIL)
      if ( TASK == 3 ) then
        CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,1)
        write(*,'(A32,I3,A45)') ">> Best fit has been found from ",ITER,
     &                   " iterations and written to the parameter file"
        RETURN
      end if

            ! If TASK=4 then output the best fit, then remove discrepant
            ! observations  (those with values more than SIGMA away from
            ! the best-fitting model value) and refit and output results
            ! If individual observational errors are available  then use
            ! these. If not then calculate the standard deviation of the
            ! light curve datapoints and use that instead. Note that the
            ! light ratios and minimum times always have uncertainties.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,2)
      write (*,'(A40,A40)') ">> Best fit to all data points has been ",
     &                      "found and written to the parameter file."

      RESIDSQSUM = 0.0d0
      do i = 1,NDATA
        if ( DTYPE(i) == 1 ) then
          MAG = GETMODEL (V,LDTYPE,DATA(1,i),DTYPE(i),LP,LS)
          OMC(i) = MAG - DATA(2,i)
          RESIDSQSUM = RESIDSQSUM + OMC(i)**2
        end if
      end do
      if (DATA(3,1) < 0.0d0) SIG = sqrt(RESIDSQSUM/real(NDATA-NLR-NMIN))

          ! Now put the array indices of the good datapoints into KEEP
      COUNT = 0
      do i = 1,NDATA                      ! if no observational errors
        if ( DTYPE(i) == 1 .and. DATA(3,i) < 0.0d0 ) then
          if ( abs(OMC(i)) <= SIG*SIGMA ) then
            COUNT = COUNT + 1
            KEEP(COUNT) = i
          end if
        else                             ! if ob'l errors were supplied
          if ( abs(OMC(i)/DATA(3,i)) <= SIGMA ) then
            COUNT = COUNT + 1
            KEEP(COUNT) = i
          end if
        end if
      end do

          ! Now keep only those datapoints which are specified by KEEP
          ! Have to recount the number of light ratios and minimum times

      do i = 1,COUNT
        DATA(1,i) = DATA(1,KEEP(i))
        DATA(2,i) = DATA(2,KEEP(i))
        DATA(3,i) = DATA(3,KEEP(i))
        DTYPE(i) = DTYPE(KEEP(i))
      end do

      write (*,'(A2,I5,A8,I5,A20,A40)') ">>",NDATA-COUNT," out of ",
     &                                 NDATA," datapoints have bee",
     &                     "n rejected from the set of observations."

      NDATA = COUNT
      NLR = 0
      NMIN = 0
      do i = 1,NDATA
        if ( DTYPE(i) == 2 ) NLR = NLR + 1
        if ( DTYPE(i) == 3 ) NMIN = NMIN + 1
      end do

      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                       VERR,IFAIL)
      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,1)
      write(*,'(A32,I3,A45)') ">> Best fit has been found from ",ITER,
     &                 " iterations and written to the parameter file"

      END SUBROUTINE TASK34
!=======================================================================
      SUBROUTINE TASK6 (V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN)
            ! Fits model parameters to an observed light curve. Then
            ! investigates each adjustable parameter by fixing it at a
            ! range of values round the best fit and seeing what the
            ! effect is on the other parameters.
      implicit none
      real*8 V(22)                        ! IN: The light curve params
      integer VARY(22)                    ! IN: Par adjustment integers
      integer LDTYPE(2)                   ! IN: LD law type foreach star
      real*8 DATA(3,99999)                ! IN: Observational data
      integer DTYPE(99999)                ! IN: Type of each datapoint
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      integer NUMVARY                     ! LOCAL: Number of vary params
      integer VWHERE(22)                  ! LOCAL: Which parameters vary
      real*8 VSTORE(22)                   ! LOCAL: Store best-fit params
      integer VARYFLAG                    ! LOCAL: Store a VARY integer
      real*8 VERR(22),CHISQ               ! LOCAL: Output from FITEBOP
      integer i,j,k,IFAIL,ITER            ! LOCAL: Loop counters etc

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,0)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                       VERR,IFAIL)
      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,2)
      write (*,'(A40,A40)')  ">> Fitting process completed. Best fit i",
     &                      "s found and outputted to parameter file. "

      NUMVARY = 0
      j = 1
      do i = 1,22
        if ( i /= 11 .and. i /= 12 ) then
          if ( VARY(i) /= 0 ) then
            NUMVARY = NUMVARY + 1
            VWHERE(j) = i
            j = j + 1
          end if
        end if
      end do

      do i = 1,22
        VSTORE(i) = V(i)
      end do

            ! Params with VARY=1 are adjustable and those with VARY=2
            ! are fixed when finding the best fit but are perturbed like
            ! adjustable parameters here.
            ! Foreach perturbed parameter store its value and adjustment
            ! integer, then step through various values whilst adjusting
            ! all adjustable parameters to find the best fit.   Start at
            ! the best-fit value  and gradually step away from it in one
            ! direction. Then do the same for the other direction.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,0)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                       VERR,IFAIL)

      do i = 1,NUMVARY
        j = VWHERE(i)
        do k = 1,22
          V(k) = VSTORE(k)
        end do
        j = VWHERE(i)
        if ( j /= 16 .and. j /= 17 ) then
          VARYFLAG = VARY(j)
          VARY(j) = 0
          CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                                                 CHISQ,VERR,100+j)

          do k = 0,20
            if ( j == 6 ) V(j) = VSTORE(j) - 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 - k/40.0d0)
            CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,
     &                                            ITER,CHISQ,VERR,IFAIL)
            CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                                                 CHISQ,VERR,200+j)
          end do

          do k = 1,22
            V(k) = VSTORE(k)
          end do

          do k = 1,20
            if ( j == 6 ) V(j) = VSTORE(j) + 0.1d0*k
            if ( j /= 6 ) V(j) = VSTORE(j) * (1.0d0 + k/40.0d0)
            if ( j /= 6 .or. V(6) < 90.1d0 ) then
              CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,
     &                                            ITER,CHISQ,VERR,IFAIL)
              CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,
     &                                            ITER,CHISQ,VERR,200+j)
            end if
          end do

          VARY(j) = VARYFLAG
        end if

      end do

      END SUBROUTINE TASK6
!=======================================================================
      SUBROUTINE TASK578 (TASK,V,VARY,LDTYPE,DATA,DTYPE,NDATA,NLR,NMIN,
     &                                                             NSIM)
!   This big subroutine executes three of the JKTEBOP tasks:
! TASK 5:  this has some characteristics of a global search. Once a best
!   fit has been found the params are perturbed and refitted many times.
!   Each result is outputted  and the best one is given as the solution.
! TASK 7:  this performs a bootstrapping error analysis. Once a best fit
!   is found many new datasets are made by copying (with replacement) at
!   random from the actual light curve. Each one is fitted to find 1-sig
!   errors. All light ratios and times of minimum are kept each time.
! TASK 8:  this does a Monte Carlo simulation error analysis.  Once best
!   fit has been found, the model is evaluated  at the phases of the obs
!   to make a synthetic light curve. Then, for each simulation, Gaussian
!   noise is added and the data are refitted to find the 1-sigma spread.
! For Tasks 7 and 8, the start params are optionally perturbed each time
!   to avoid finding unrealistically small errors in a local minimum.
! Monte Carlo simulations  are an excellent error indicator,  but do not
!   fully account for systematic errors.   Bootstrapping is better,  but
!   the correct solution to this problem is to remove the systematics.
! For information on bootstrapping, Monte Carlo, and minimisation algor-
! ithms, read Numerical Recipes in Fortran 77 (Press et al 1993) chap.15
      implicit none
      integer TASK                  ! IN: The task to undertake
      real*8 V(22)                  ! IN: The photometric parameters
      integer VARY(22)              ! IN: Parameter adjustment integers
      integer LDTYPE(2)             ! IN: Type of LD law for each star
      real*8 DATA(3,99999)          ! IN: Observational data
      integer DTYPE(99999)          ! IN: Types of observational data
      integer NDATA,NLR,NMIN        ! IN: Numbers of different datatypes
      integer NSIM                  ! IN: Number of  simulations  to  do
      character PERTURB*1           ! LOCAL: Perturb params ('y' or 'n')
      real*8 PFACTOR                ! LOCAL: Scale factor for perturbatn
      real*8 VERR(22)               ! SUB:   Param  formal uncertainties
      real*8 DV(22)                 ! LOCAL: Perturbations for paramters
      real*8 VEXTRA(6)              ! SUB:   Extra (dependent) parametrs
      real*8 VALL(28,100000)        ! LOCAL: All the best-fitting params
      integer NVARY                 ! LOCAL: Number of varying parametrs
      integer ITER                  ! SUB:   Numberof fitting iterations
      real*8 INDATA(3,99999)        ! LOCAL: Synthetic obs  to be fitted
      real*8 INV(22)                ! SUB:   Parameters of  synthetic LC
      real*8 CHISQ,ERRSIZE          ! SUB:   Chi-squared of the best fit
      integer SEEDSTART,SEED        ! FUNCTION: Start randomnumber maker
      real*8 RANDOMG                ! FUNCTION: Gaussian  random  number
      real*8 RANDOM                 ! FUNCTION: Flat-distrib  random num
      real*8 GETMODEL               ! FUNCTION: Gets  model  predictions
      real*8 LP,LS,MAG              ! LOCAL: EBOP/GETMODEL output values
      real*8 HELP1                  ! LOCAL: Useful variable storage
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters + error flags

            ! Define parameter labels for output and the amounts to vary
            ! each by if they have not been adjusted to the best fit but
            ! are required to be perturbed during Monte Carlo analysis.
            ! A positive number indicates that it is a multiplicative
            ! factor whilst a negative number that it is additive.

      DATA DV / 0.1d0, 0.1d0, 0.2d0,-0.1d0,-0.1d0,-1.0d0, 0.1d0, 0.1d0,
     &          0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.0d0,-0.05d0,0.1d0,
     &          0.1d0, 0.0d0, 1.d-6,-1.d-4,-0.1d0,-0.1d0 /

            ! Amount to perturb parameters by before each refit. For
            ! Tasks 7 and 8 this should be the amounts specified in DV
            ! whereas for Task 5 it is four times the amounts.
            ! Then determine whether perturbations should occur at all.

      PFACTOR = 1.0d0
      if ( TASK == 5 ) PFACTOR = 4.0d0

      if ( NSIM < 0 ) then
        PERTURB = "n"
        NSIM = abs(NSIM)
      else
        PERTURB = "y"
      end if

            ! Now find the best fit to the light curve and output it.

      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,0)
      CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                       VERR,ERROR)
      if ( TASK == 5 ) then
        write(62,'(A44)') "---------------------------------------     "
        write(62,'(A44)') "Best fit found from initial parameters      "
      end if
      CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,2)

      write (*,'(A40,A40)')  ">> Fitting process completed. Best fit i",
     &                      "s found and outputted to parameter file. "

            ! Store some dependent quantities for output later on.
            ! Write column headings for simulation results output file.

      VEXTRA(1) = V(2) / (1.0d0 + V(3))                      ! r_1
      VEXTRA(2) = V(2) / (1.0d0 + (1.0d0/V(3)))              ! r_2
      VEXTRA(3) = GETMODEL (V,LDTYPE,0.25d0+V(16),2,LP,LS)   ! L_2 / L_1
      VEXTRA(4) = sqrt(V(7)**2 + V(8)**2)                    ! e
      VEXTRA(5) = atan2(V(8),V(7))*57.2957795d0              ! omega
      if (VEXTRA(5)<0.0d0)  VEXTRA(5)=VEXTRA(5)+360.0d0
      VEXTRA(6) = CHISQ

      write (63,'(A)')  "#   N  ITER     SB2                r_1+r_2  "//
     &     "           k                 LDu1               LDu2     "//
     &     "           i                 e(cos)w            e(sin)w  "//
     &     "          GD1                GD2                refl1    "//
     &     "          refl2               q                 tidalangl"//
     &     "e         L_3                phase_corr         sfact    "//
     &     "          int_ring            P               T_zero     "//
     &     "          LDn1               LDn2               r_1      "//
     &     "          r_2                L_2 / L_1           e       "//
     &     "         omega               (chi)^2"

            ! Start the random number generator. Calculate the number of
            ! variable params.  Store original data and best-fit params.

      SEED = SEEDSTART ()

      NVARY = 0
      do i = 1,22
        if ( VARY(i) == 1 ) NVARY = NVARY + 1
      end do

      do i = 1,22
        INV(i) = V(i)
      end do

      do i = 1,NDATA
        INDATA(1,i) = DATA(1,i)
        INDATA(2,i) = DATA(2,i)
        INDATA(3,i) = DATA(3,i)
      end do

            ! If the inputted NSIM is below zero no perturbations should
            ! be applied to the initial  light curve parameter estimates
            ! prior to each MC simulation. This is useful if convergence
            ! only happens for a  very narrow range of parameter values.

      if ( TASK /= 5 ) then
        if (PERTURB == "n") write (62,'(A18,A60)') ">> Number of simul",
     &   "ations is below zero: the parameters will not be perturbed. "
        if (PERTURB == "y") write (62,'(A18,A60)') ">> Number of simul",
     &   "ations is above zero: the parameters will be perturbed.     "
        write (62,*) " "
      end if

      if ( TASK == 5 ) write(*,'(A46,I6)')
     &             ">> Number of perturbed parameter refits to do:",NSIM
      if ( TASK == 7 ) write(*,'(A45,I6)')
     &              ">> Number of bootstrapping simulations to do:",NSIM
      if ( TASK == 8 ) write(*,'(A43,I6)')
     &                ">> Number of Monte Carlo simulations to do:",NSIM
      write(*,'(A13,$)') ">> Completed:"

            ! For Monte Carlo simulations must evaluate the best-fitting
            ! model light curve at the observed phases. Also, if there
            ! were no errors in the input light curve, must calculate
            ! the rms of the residuals of the fit and use this to scale
            ! the Gaussian noise added to the simulated light curve.

      if ( TASK == 8 ) then

        HELP1 = 0.0d0
        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) then
            DATA(2,i) = GETMODEL (V,LDTYPE,DATA(1,i),DTYPE(i),LP,LS)
            HELP1 = HELP1 + (DATA(2,i) - INDATA(2,i))**2
          end if
        end do

        if ( DATA(3,1) <= 0.0d0 ) then
          ERRSIZE = sqrt(HELP1 / (NDATA-NVARY-NLR-NMIN) )
          do i = 1,NDATA
            if ( DTYPE(i) == 1 ) INDATA(3,i) = ERRSIZE
          end do
          write (62,'(A36,A44)') "No observational errors were supplie",
     &                   "d with the input light curve.  Noise of size"
          write (62,'(F6.3,A14,A60)')abs(ERRSIZE)*1.d3,"mmag (standard",
     &   " error of best fit) will be added to the synthetic datasets."
        else
          write (62,'(A36,A44)') "Observational errors were supplied w",
     &                   "ith the input light curve.  These have been "
          write (62,'(A36,A44)') "assumed to be correct and used to se",
     &                   "t the size of the simulated Gaussian noise. "
        end if
      end if

!-----------------------------------------------------------------------

      do i = 1,NSIM                  ! Loop over number  of  simulations
500     continue                     ! Enter here if previous one failed

            ! For bootstrapping, randomly sample the observations (with
            ! replacement) to create a new light curve to fit. Don't do
            ! this with the light ratios or minimum times - they should
            ! be preserved as they are.

        if ( TASK == 7 ) then
          do j = 1,NDATA
            if ( DTYPE(j) == 1 ) then
              do m = 1,100000
                k = 1 + int( random(SEED) * (real(NDATA)-0.000001) )
                if ( DTYPE(k) == 1 ) exit
              end do
              INDATA(1,j) = DATA(1,k)
              INDATA(2,j) = DATA(2,k)
              INDATA(3,j) = DATA(3,k)
            else
              INDATA(1,j)= DATA(1,j)
              INDATA(2,j)= DATA(2,j)
              INDATA(3,j)= DATA(3,j)
            end if
          end do
        end if

            ! For Monte Carlo, create a simulated light curve by adding
            ! Gaussian noise to the best-fit model light curve (which
            ! has been put in the DATA array, replacing the actual obs).

        if ( TASK == 8 ) then
          do j = 1,NDATA
            INDATA(2,j) = DATA(2,j) + randomg(SEED,0.0d0,INDATA(3,j))
          end do
        end if

            ! Now perturb the initial values of the fitted parameters.

        do j = 1,22
          if ( (VARY(j)==1 .and. PERTURB=="y") .or. VARY(j)==2 ) then
            if ( DV(j) > 0.0d0 ) then
              INV(j) = V(j) + PFACTOR*(RANDOM(SEED)-0.5)*V(j)*DV(j)
            else if ( DV(j) <= 0.0d0 ) then
              INV(j) = V(j) + PFACTOR*(RANDOM(SEED)-0.5)*abs(DV(j))
            end if
          end if
        end do

        if ( INV(6) > 90.0d0 ) GOTO 500         ! Avoid inclination > 90

            ! Now fit the simulated light curve and store the results

        CALL FITEBOP (INDATA,DTYPE,NDATA,NLR,NMIN,INV,VARY,LDTYPE,ITER,
     &                                                 CHISQ,VERR,ERROR)
        if ( ERROR /= 0 )  GOTO 500

        do j = 1,22
          VALL(j,i) = INV(j)
        end do
        VALL(23,i) = INV(2) / (1.0d0 + INV(3))                 ! r_1
        VALL(24,i) = INV(2) / (1.0d0 + (1.0d0/INV(3)))         ! r_2
        VALL(25,i) = GETMODEL(INV,LDTYPE,0.25d0+V(16),2,LP,LS) ! L2/L1
        VALL(26,i) = sqrt(INV(7)**2 + INV(8)**2)               ! e
        VALL(27,i) = atan2(INV(8),INV(7))*57.2957795d0         ! omega
        if (VALL(27,i)<0.0d0) VALL(27,i)=VALL(27,i)+360.0d0
        VALL(28,i) = CHISQ

            ! Write result to the big output file. Output the simulation
            ! number to screen if required. Check that a halt has not
            ! been requested. Then finish the iteration loop.

        write(63,'(I5,1X,I3,30(1X,f18.10))')
     &                        i,ITER,(INV(j),j=1,22),(VALL(j,i),j=23,28)

        if ( i < 10 ) write (*,'(1X,I1,$)') i
        if ( (real(i)/10.0) == int(real(i)/10.0) ) then
          if ( i >= 10 .and. i < 100 )   write (*,'(1X,I2,$)') i
          if ( i >= 100 .and. i < 1000 )  write (*,'(1X,I3,$)') i
          if ( i >= 1000 .and. i < 10000 ) write (*,'(1X,I4,$)') i
          if ( i >= 10000 )                 write (*,'(1X,I5,$)') i
        end if

        CALL STOPCHECK (ERROR)
        if ( ERROR /= 0 ) then
          NSIM = i
          exit
        end if
      end do
      write (*,*) " "

!-----------------------------------------------------------------------

            ! For Task 5, pick out the result with the lowest chi-square
            ! value, put its parameters into V, refit the result, and
            ! print the final quantities to the output file.

      if ( TASK == 5 ) then
        j = 0
        HELP1 = 1.0d8
        do i = 1,NSIM
          if ( VALL(28,i) < HELP1 ) then
            j = i
            HELP1 = VALL(28,i)
          end if
        end do

        do i = 1,22
          V(i) = VALL(i,j)
        end do
        CHISQ = HELP1
        CALL FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ
     &                                                      ,VERR,ERROR)

        write(62,'(A44)') "----------------------------------------    "
        write(62,'(A44)') "Overall best fit from all parameter sets    "
        CALL OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,CHISQ,
     &                                                           VERR,2)
      end if

            ! For Task 7 and 8, call OUTPUTSIM to write out 1-sigma unc-
            ! retainties in each of the fitted and dependent parameters.

      if ( TASK == 7 .or. TASK == 8 )
     &                     CALL OUTPUTSIM (TASK,V,VARY,VEXTRA,VALL,NSIM)

      END SUBROUTINE TASK578
!=======================================================================
!=======================================================================
      SUBROUTINE STOPCHECK (STATUS)  ! Checks if the file "jktebop.stop"
                                     ! exists and if so puts STATUS=1000
      implicit none
      integer ERROR,STATUS

      STATUS = 0
      ERROR = 100
      open (99,file="jktebop.stop",status="old",iostat=ERROR)
      close (99,status="delete")
      if ( ERROR == 0 ) STATUS = 1000

      END SUBROUTINE STOPCHECK
!=======================================================================
      SUBROUTINE OUTPUTSIM (TASK,V,VARY,VEXTRA,VALL,NSIM)
      implicit none
      integer TASK                  ! IN: Which task for  output results
      real*8 V(22)                  ! IN: The photometric model paramtrs
      integer VARY(22)              ! IN: Which parametrs fixed/adjusted
      real*8 VEXTRA(6)              ! IN: Extra  (dependent)  parameters
      real*8 VALL(28,100000)        ! IN: All the best-fitting parametrs
      integer NSIM                  ! IN: The number of simulations done
      real*8 ARRAY(NSIM)            ! LOCAL: All results for one paramtr
      real*8 HELP1,HELP2,HELP3      ! LOCAL: Useful storage for variabls
      real*8 PERCENT                ! LOCAL: Percentage error for a par.
      character*5 NAME(28)          ! LOCAL: Names of the variables etc.
      integer i,j,k                 ! LOCAL: Loop counters
      real*8 SELECT,SIGMA           ! FUNCTIONS: Give standard deviation

      DATA NAME/"    J","r1+r2","    k","LD_u1","LD_u2","  inc","ecosw",
     &          "esinw","  GD1","  GD2","refl1","refl2","    q","T_lag",
     &          "  L_3"," d(P)","sfact","iring","P_orb"," T_0 ","LD_n1",
     &          "LD_n2","  r_1","  r_2","L2/L1","    e","omega","Rchi2"/

            ! Now for each adjusted parameter (and the other calculated
            ! quantities, put the NSIM parameter evaluations into an
            ! array and select the median and 1_sigma bounds (using the
            ! median +/- 0.5_sigma) and output to the parameter file.

      write (62,*) " "
      write (62,'(A40,A54)') "========================================",
     &         "======================================================"
      if ( TASK == 7 ) write (62,'(A13,I6,A27)')
     &          "Results from ",abs(NSIM)," bootstrapping simulations:"
      if ( TASK == 8 ) write (62,'(A13,I6,A25)')
     &                 "Results from ",NSIM," Monte Carlo simulations:"
      write (62,'(A40,A54)') "========================================",
     &         "======================================================"
      write (62,*) " "
      write (62,'(A40,A54)') "              Best fit     std.deviation",
     &         "       median         +68.3%        -68.3%        %   "

      do j = 1,28
        do k = 1,NSIM
          ARRAY(k) = mod(VALL(j,k),1.0d3)
        end do
        HELP1=SELECT(ARRAY,NSIM,int(NSIM*0.1585d0))     ! Get median and
        HELP2=SELECT(ARRAY,NSIM,int(NSIM*0.5d0))        ! 68.3% (1sigma)
        HELP3=SELECT(ARRAY,NSIM,int(NSIM*0.8415d0))     ! intervals.
        PERCENT = 50.d0 * abs((abs(HELP3-HELP2)+abs(HELP1-HELP2))/HELP2)

        if ( j <= 22 .and. VARY(j) /= 0 ) then
          if ( j == 1 .or. j == 2 .or. j == 3 ) then
            write (62,100) j,NAME(j),mod(V(j),1.0d3),SIGMA(ARRAY,NSIM)
     &                 ,mod(HELP2,1.0d3),HELP2-HELP1,HELP3-HELP2,PERCENT
          else if ( j /= 14 .and. j /= 17 ) then
            write (62,101) j,NAME(j),mod(V(j),1.0d3),SIGMA(ARRAY,NSIM)
     &                         ,mod(HELP2,1.0d3),HELP2-HELP1,HELP3-HELP2
          end if
        end if

        if ( j == 23 .or. j == 24 .or. j == 25 ) then
          write (62,102) NAME(j),VEXTRA(j-22),SIGMA(ARRAY,NSIM),
     &                             HELP2,HELP2-HELP1,HELP3-HELP2,PERCENT
        else if ( j == 26 .or. j == 27 .or. j == 28) then
          write (62,103) NAME(j),VEXTRA(j-22),SIGMA(ARRAY,NSIM),
     &                                     HELP2,HELP2-HELP1,HELP3-HELP2
        end if
      end do
      write (62,*) " "

100   FORMAT (I2,1X,A5,1X,F14.10," +/-",F13.10,2X,
     &                    F14.10," +",F13.10," -",F13.10," (",F5.2,"%)")
101   FORMAT (I2,1X,A5,1X,F14.10," +/-",F13.10,2X,
     &                    F14.10," +",F13.10," -",F13.10)
102   FORMAT    (3X,A5,1X,F14.10," +/-",F13.10,2X,
     &                    F14.10," +",F13.10," -",F13.10," (",F5.2,"%)")
103   FORMAT    (3X,A5,1X,F14.10," +/-",F13.10,2X,
     &                    F14.10," +",F13.10," -",F13.10)

      END SUBROUTINE OUTPUTSIM
!=======================================================================
      SUBROUTINE OUTPUT (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                                                  CHISQ,VERR,WHAT)
            ! This subroutine writes information to the output files.
            ! Precisely what to write is given by parameter WHAT.
            ! If WHAT=0, output initial parameters and useful quantities
            ! If WHAT=1, then output the final fitted parameters, useful
            !   quantities, residual curve, and best-fitting model curve
            ! If WHAT=2, output best fit parametrs and useful quantities
            ! If WHAT=101-116 then this has been invoked from TASK6 to
            !   output the name of the parameter=(WHAT-100)
            ! If WHAT=201-216 then this has come from TASK6 to output
            !   fit results for a certain value of parameter=(WHAT-200)
      implicit none
      integer WHAT                        ! IN:  Indicate what to output
      real*8 DATA(3,99999)                ! IN:  Observational data
      integer DTYPE(99999)                ! IN:  Type of each data point
      integer NDATA,NLR,NMIN              ! IN:  Numbers of  data points
      real*8 V(22),VERR(22)               ! IN:  Parameters  and  errors
      integer VARY(22)                    ! IN:  Par adjustment integers
      integer LDTYPE(2)                   ! IN:  LDlaw type foreach star
      integer ITER                        ! IN:  numbr of fit iterations
      real*8 CHISQ                        ! IN:   chi-squared of the fit
      integer NVARY,VWHERE(22)            ! LOCAL: which parameters vary
      character VNAME(22)*18,VNAM(26)*5   ! LOCAL: names  of  parameters
      real*8 ECC,OMEGA,ECOSW,ESINW        ! LOCAL: orbital    parameters
      real*8 RP,RS                        ! LOCAL: fractional star radii
      real*8 RESIDSUM,RESIDABSSUM         ! LOCAL: sums of the residuals
      real*8 RESIDSQSUM,SE                ! LOCAL: resid, standard error
      real*8 LTOTAL                       ! LOCAL: total light of system
      real A,B,EPSILON1,EPSILON2,EPSILON  ! LOCAL: stellar shape  params
      real*8 MAG,LP,LS                    ! LOCAL: LIGHT subroutine pars
      real*8 PHASE,HELP                   ! LOCAL: some useful variables
      integer i,j,k,ERROR,STATUS          ! LOCAL: counters & errorflags
      real*8 GETPHASE                     ! FUNCTION: calc orbital phase
      real*8 GETMODEL                     ! FUNCTION: calcs model output

      DATA VNAME/             "Surface br ratio  ","Sum of the radii  ",
     &   "Ratio of the radii","Limb darkening u1 ","Limb darkening u2 ",
     &   "Orbit inclination ","ecc cos(omega)    ","ecc sin(omega)    ",
     &   "Gravity dark 1    ","Gravity dark 2    ","Reflected light 1 ",
     &   "Reflected light 2 ","Mass ratio q      ","Tidal lead/lag ang",
     &   "Third light       ","Phase correction  ","Light scale factor",
     &   "Integration ring  ","Orbital period    ","Ephemeris timebase",
     &   "Limb darkening n1 ","Limb darkening n2 "/
      DATA VNAM/"  J  ","r1+r2","  k  ","LD_u1","LD_u2"," inc ","ecosw",
     &          "esinw"," GD1 "," GD2 ","refl1","refl2","  q  ","T_lag",
     &          " L_3 "," d(P)","sfact","iring","P_orb"," T_0 ","LD_n1",
     &          "LD_n2"," r_1 "," r_2 ","L2/L1"," rms "/

            ! Find stellar radii, reflection coeffs, total system light
            ! and the number of adjustable parameters.

      RP = V(2) / (1.0d0 + V(3))
      RS = V(2) / (1.0d0 + (1.0d0/V(3)))

      MAG = GETMODEL(V,LDTYPE,0.0d0,1,LP,LS)
      LTOTAL = LP + LS + V(15)
      if ( VARY(11) == -1 )  V(11) = 0.4d0 * LS * RP**2
      if ( VARY(12) == -1 )  V(12) = 0.4d0 * LP * RS**2

      NVARY = 0
      do i = 1,22
        if ( VARY(i) == 1 ) NVARY = NVARY + 1
      end do

            ! Output the values of the parameters and various other info

      if ( WHAT == 0 ) then
        write(62,107) "Initial values of the parameters:               "
      else if ( WHAT == 1 .or. WHAT == 2 ) then
        write(62,*)   " "
        write(62,107) "---------------------------------------         "
        write(62,110) ITER," iterations of EBOP completed     "
        write(62,107) "---------------------------------------         "
        write(62,*)   " "
        write(62,105)   "Warning: the uncertainties outputted imm",
     &                  "ediately below are formal               "
        write(62,105)   "errors and cannot be fully trusted. Use ",
     &                  "Task 6,7,8 for good errors.             "
        if ( DATA(3,1) < 0.0d0 ) then
          write(62,105) "The parameter errors have been scaled by",
     &                  " sqrt{reduced chi^2}"
        end if
        write(62,*)   " "
        write(62,107) "Final values of the parameters:                 "
      end if

      if ( WHAT == 0 .or. WHAT == 1  .or. WHAT == 2 ) then
        if ( V(6) > 90.0d0 )  V(6) = 180.0d0 - V(6)
        if ( V(6) < 0.0d0 ) V(6) = V(6) + 180.0d0
        if ( LDTYPE(2) == 0 ) then
          V(5) = V(4)
          VERR(5) = VERR(4)
        end if

        if ( WHAT == 0 ) then
          do i = 1,22
            if ( i == 7 .and. V(7) > 5.0d0 ) then
              if (VARY(i)/=1) write(62,100) i,"Eccentricity      ",
     &                                        V(i)-10.0d0
              if (VARY(i)==1) write(62,101) i,"Eccentricity      ",
     &                                        V(i)-10.0d0,"  (adjusted)"
            else if ( i == 8 .and. V(7) > 5.0d0  ) then
              if (VARY(i)/=1) write(62,100) i,"Periast. longitude",V(i)
              if (VARY(i)==1) write(62,101) i,"Periast. longitude",V(i),
     &                                                    "  (adjusted)"
            else
              if (VARY(i)/=1) write(62,100) i,VNAME(i),V(i)
              if (VARY(i)==1)write(62,101)i,VNAME(i),V(i),"  (adjusted)"
            end if
          end do
        else
          do i = 1,22
            if ( i == 7 .and. V(7) > 5.0d0 ) then
              if ( VARY(i) /= 1 ) write (62,113) i,"Eccentricity      ",
     &                                           V(i)-10.0d0,"  (fixed)"
              if ( VARY(i) == 1 ) write (62,114) i,"Eccentricity      ",
     &                                           V(i)-10.0d0,VERR(i)
            else if ( i == 8 .and. V(7) > 5.0d0  ) then
              if ( VARY(i) /= 1 ) write (62,113) i,"Periast. longitude",
     &                                           V(i)-10.0d0,"  (fixed)"
              if ( VARY(i) == 1 ) write (62,114) i,"Periast. longitude",
     &                                           V(i)-10.0d0,VERR(i)
            else
              if ( i == 11 .and. VARY(i) == -1 ) then
                write(62,111) i,VNAME(i),V(i),"  (from geometry)"
              else if ( i == 12 .and. VARY(i) == -1 ) then
                write(62,111) i,VNAME(i),V(i),"  (from geometry)"
              else
                if(VARY(i)/=1) write(62,113)i,VNAME(i),V(i),"  (fixed)"
                if(VARY(i)==1) write(62,114)i,VNAME(i),V(i),VERR(i)
              end if
            end if
          end do
        end if

        write (62,*) " "

116     FORMAT ("Limb darkening type for primary star:    ",A16)
117     FORMAT ("Limb darkening type for secondary star:  ",A16)
        if ( LDTYPE(1) == 1 ) write (62,116) "linear          "
        if ( LDTYPE(1) == 2 ) write (62,116) "logarithmic     "
        if ( LDTYPE(1) == 3 ) write (62,116) "square-root     "
        if ( LDTYPE(1) == 4 ) write (62,116) "quadratic       "
        if ( LDTYPE(2) == 1 ) write (62,117) "linear          "
        if ( LDTYPE(2) == 2 ) write (62,117) "logarithmic     "
        if ( LDTYPE(2) == 3 ) write (62,117) "square-root     "
        if ( LDTYPE(2) == 4 ) write (62,117) "quadratic       "
        if ( LDTYPE(2) == 0 ) write (62,117) "same as primary "

        if ( WHAT == 0 .and. (LDTYPE(1) == 2 .or. LDTYPE(2) == 2) ) then
          write (62,*) " "
          write (62,105) "### WARNING: the logarithmic LD law flux",
     &                   " normalisation is only approximate.     "
          write (62,105) "Do not use this law if you need a very p",
     &                   "recise (better than 1%) light ratio.    "
        end if

        write (62,*) " "

        write(62,102)"Fractional primary radius:          ",RP
        write(62,102)"Fractional secondary radius:        ",RS
        if ( V(2) >= 0.6d0 ) then
          write (62,*) " "
          write (62,105) "### WARNING: average radius is greater t",
     &                   "han 0.3 so the radii may be wrong by 5%!"
          write (62,105) "### See North & Zahn (2004, New.Ast.Rev.",
     &                   ", 47, 741) for details. Do not use EBOP!"
          write (62,*) " "
        else if ( V(2) >= 0.5d0 )  then
          write (62,*) " "
          write (62,105) "### WARNING: average radius is greater t",
     &                   "han 0.25 so the radii may be wrong by 1%"
          write (62,105) "### See North & Zahn (2004, New Astronom",
     &                   "y Review, 47, 741) for further details. "
          write (62,*) " "
        end if

        write(62,102)"Stellar light ratio (at phase 0.25):",LS/LP
        write(62,102)"Primary contribut'n to system light:",LP/LTOTAL
        write(62,102)"Secondary contrib'n to system light:",LS/LTOTAL
        write(62,*)  " "

        if ( V(7) > 5.0d0 ) then
          ECOSW = (V(7)-10.0d0) * cos(V(8)/57.2957795d0)
          ESINW = (V(7)-10.0d0) * sin(V(8)/57.2957795d0)
          write (62,102) "Eccentricity*cos(omega):            ",ECOSW
          write (62,102) "Eccentricity*sin(omega):            ",ESINW
        else if ( V(7) /= 0.0d0 .and. V(8) /= 0.0d0 ) then
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
          write (62,102) "Eccentricity:                       ",ECC
          write (62,102) "Omega (degrees):                    ",OMEGA
        end if

        if ( VARY(11) /= 1 .or. VARY(12) /= -1 ) then
          HELP = 0.4d0 * LS * RP**2
          write (62,102) "Geometric reflection coeff (star A):",HELP
          HELP = 0.4d0 * LP * RS**2
          write (62,102) "Geometric reflection coeff (star B):",HELP
        end if

        CALL BIAX (real(RP),real(abs(V(13))),A,B,EPSILON1)
        CALL BIAX (real(RS),real(1.0/abs(V(13))),A,B,EPSILON2)
        if ( V(13) <= 0.0d0 ) EPSILON1 = 0.0d0
        if ( V(13) <= 0.0d0 ) EPSILON2 = 0.0d0

        write (62,102) "Oblateness of the primary star:     ",EPSILON1
        write (62,102) "Oblateness of the secondary star:   ",EPSILON2
        if ( (WHAT == 1 .or. WHAT == 2) .and. V(13) <= 0.0d0 ) then
          CALL BIAX(real(RP),real(abs(V(13))),A,B,EPSILON)
          write(62,102)"Expected oblateness of primary:     ",EPSILON
          CALL BIAX(real(RS),real(1./abs(V(13))),A,B,EPSILON)
          write(62,102)"Expected oblateness of secondary:   ",EPSILON
        end if

        if ( EPSILON1 > 0.04 .or. EPSILON2 > 0.04 )  then
          write (62,*) " "
          write (62,105) "### WARNING: oblateness is above the rec",
     &                   "ommended maximum value for EBOP of 0.04."
          write (62,105) "See Popper & Etzel (1981, AJ, 86, 102) f",
     &                   "or a justification and further details. "
        end if
        write (62,*) " "
      end if

!-----------------------------------------------------------------------
            ! Now output the observations and the calculated residuals.

      if ( WHAT == 1 .or. WHAT == 2 ) then
        if ( WHAT == 1 ) write (63,'(A21,A44)') "#     TIME      MAGNI",
     &                   "TUDE    ERROR      PHASE      MODEL    (O-C)"

        RESIDSUM = 0.0d0
        RESIDABSSUM = 0.0d0
        RESIDSQSUM = 0.0d0
        do i = 1,NDATA
          if ( DTYPE(i) == 1 ) then
            MAG = GETMODEL (V,LDTYPE,DATA(1,i),DTYPE(i),LP,LS)
            if ( WHAT == 1 ) write(63,104)DATA(1,i),DATA(2,i),DATA(3,i),
     &               GETPHASE(DATA(1,i),V(19),V(20)),MAG,(DATA(2,i)-MAG)
            RESIDSUM = RESIDSUM + DATA(2,i) - MAG
            RESIDABSSUM = RESIDABSSUM + abs(DATA(2,i)-MAG)
            RESIDSQSUM = RESIDSQSUM + (DATA(2,i)-MAG)**2
          end if
        end do

        write(62,102) "Sum of residuals for LC data (mag): ",   RESIDSUM
        write(62,102) "Sum of absolute values of LC resids:",RESIDABSSUM
        write(62,102) "Sum of the squares of the LC resids:", RESIDSQSUM
        write(62,115) "Total datapoints / fitting params:  ",NDATA,NVARY
        if ( WHAT == 1 .or. WHAT == 2 ) then
          write (62,102) "Standard error (one LC obs, mmag):  ",
     &              sqrt(RESIDSQSUM / (NDATA-NVARY-NLR-NMIN)) * 1000.0d0
          write (62,102) "rms of residuals (one LC obs, mmag):",
     &                    sqrt(RESIDSQSUM / (NDATA-NLR-NMIN)) * 1000.0d0
          if ( CHISQ > 0.0d0 )
     &    write (62,102) "Reduced chi-squared from errorbars: ", CHISQ
        end if
        write (62,*) " "
      end if

      if ( WHAT == 1 ) then
        write (64,'(A22)') "#  PHASE    MAGNITUDE "
        do i = 0,10000
          PHASE = i / 10000.0d0
          write(64,12)PHASE,GETMODEL(V,LDTYPE,V(20)+V(19)*PHASE,1,LP,LS)
12        FORMAT (F8.6,1X,F12.8)
        end do
      end if

!-----------------------------------------------------------------------
      ! Now output the results and residuals for the times of minimum/a
      ! and the spectroscopic light ratio(s).

118   FORMAT (F9.1,1X,F13.5,1X,F8.5,2X,F13.5,1X,F6.2)
119   FORMAT (F13.5,1X,F7.4,1X,F6.4,2X,F7.4,F7.2)

      if ( WHAT == 1 .or. WHAT == 2 ) then
        if ( NLR > 0 ) then
          write (62,1071) "--------------------------------------------"
          write (62,1071) "Results for the spectroscopic light ratios  "
          write (62,1071) "--------------------------------------------"
          write (62,1071) "    Time(HJD)  Ratio   Error   Model   Sigma"
          do i = 1,NDATA
            if ( DTYPE(i) == 2 ) then
              HELP = GETMODEL (V,LDTYPE,DATA(1,i),DTYPE(i),LP,LS)
              write (62,119) DATA(1,i),DATA(2,i),DATA(3,i),HELP,
     &                                        (HELP-DATA(2,i))/DATA(3,i)
            end if
          end do
!           write (62,1071) "--------------------------------------------"
          write (62,*) " "
        end if

        if ( NMIN > 0 ) then
          write (62,105) "----------------------------------------",
     &                   "--------------                          "
          write (62,105) "Results for the times of minimum light  ",
     &                   "                                        "
          write (62,105) "----------------------------------------",
     &                   "--------------                          "
          write (62,105) "    Cycle    Time(HJD)    Error      Mod",
     &                   "el(HJD)  Sigma                          "
          do i = 1,NDATA
            if ( DTYPE(i) == 3 ) then
              HELP = GETMODEL (V,LDTYPE,DATA(1,i),DTYPE(i),LP,LS)
              write (62,118) DATA(1,i),DATA(2,i),DATA(3,i),HELP,
     &                                        (HELP-DATA(2,i))/DATA(3,i)
            end if
          end do
!           write (62,1071) "--------------------------------------------"
          write (62,*) " "
          write (62,*) " "
        end if
      end if

      if ( WHAT > 100 .and. WHAT < 117 ) then
        j = 0
        do i = 1,22
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
        write (62,*) " "
        write (62,'(A41,A20)')
     &       "Now investigating the fitting parameter: ",VNAME(WHAT-100)
        write (62,'(A4,2X,24(A5,6X))')   "ITER",VNAM(WHAT-100),VNAM(26),
     &                (VNAM(VWHERE(k)),k=1,j),VNAM(23),VNAM(24),VNAM(25)
        write(*,'(A44,A20)')
     &   ">> Now investigating the fitting parameter:  ",VNAME(WHAT-100)
      end if

!-----------------------------------------------------------------------

      if ( WHAT > 200 .and. WHAT < 217 ) then
        RESIDSQSUM = 0.0d0
        do i = 1,NDATA
          MAG = GETMODEL (V,LDTYPE,DATA(1,i),1,LP,LS)
          RESIDSQSUM = RESIDSQSUM + (DATA(2,i)-MAG)**2
        end do
        SE = sqrt(RESIDSQSUM / NDATA) * 1000.0d0
        j = 0
        do i = 1,22
          if ( VARY(i) == 1 ) then
            j = j + 1
            VWHERE(j) = i
          end if
        end do
        write(62,'(I3,23(F11.7))')    ITER,V(WHAT-200),SE,(V(VWHERE(k)),
     &          k=1,j),V(2)/(1.0d0+V(3)),V(2)/(1.0d0+(1.0d0/V(3))),
     &         GETMODEL (V,LDTYPE,0.25d0+V(16),2,LP,LS)
      end if
!-----------------------------------------------------------------------

100   FORMAT (I2,2X,A18,2X,F18.10)
101   FORMAT (I2,2X,A18,2X,F18.10,A12)
102   FORMAT (A36,2X,F16.10)
1021  FORMAT (A36,2X,E16.9)
104   FORMAT (F14.6,1X,F9.4,1X,F9.5,3X,F10.8,1X,F9.5,1X,F8.5)
105   FORMAT (A40,A40)
107   FORMAT (A48)
1071  FORMAT (A44)
108   FORMAT (F14.6,F10.6)
110   FORMAT (I6,A34)
111   FORMAT (I2,2X,A18,2X,F18.10,A17)
112   FORMAT (A39,F7.4)
113   FORMAT (I2,2X,A18,2X,F18.10,A9)
114   FORMAT (I2,2X,A18,2X,F18.10," +/- ",F14.10)
115   FORMAT (A36,7X,I6,1X,I4)

      END SUBROUTINE OUTPUT
!=======================================================================
!=======================================================================
      DOUBLEPRECISION FUNCTION GETPHASE (HJD,PERIOD,TZERO)
            ! Returns phase from given time and orbital ephemeris
      implicit none
      real*8 HJD,PERIOD,TZERO

      GETPHASE = (HJD - TZERO) / PERIOD
      GETPHASE = GETPHASE - int(GETPHASE)
      if ( GETPHASE < 0.0d0 ) GETPHASE = GETPHASE + 1.0d0

      END FUNCTION GETPHASE
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMIN (TZERO,PERIOD,ECC,OMEGA,CYCLE)
            ! Returns time of minimum for given cycle and ephemeris.  If
            ! the orbit is circular thenthe cycle number is used without
            ! restriction so can refer to any phase. If the orbit is ec-
            ! centric then the cycle number should be integer (indicates
            ! primary minimum) or half-integer (secondary minimum).
      implicit none
      real*8 TZERO,PERIOD           ! IN: reference time, orbital period
      real*8 ECC,OMEGA              ! IN: orbital (e,w) or (ecosw,esinw)
      real*8 CYCLE                  ! IN: cycle number of minimum to use
      real*8 E,W                    ! LOCAL: eccentricity and peri.long.
      real*8 THETA,EE,TANEE         ! LOCAL: true and eccentric anomaly
      real*8 PSEP                   ! LOCAL: phase diff between minima
      real*8 PI,DEG2RAD,TINY        ! LOCAL: useful variables

      PI = atan(1.0d0) * 4.0d0
      DEG2RAD = 45.0d0 / atan(1.0d0)
      TINY = 1.0d-6

            ! First must deal with the possibility that e and omega are
            ! actually e*cos(omega) and e*sin(omega)

      if ( ECC > 9.0d0 ) then
        E = ECC - 10.0d0
        W = OMEGA
      else
        E = sqrt( OMEGA*OMEGA + ECC*ECC )
        W = atan2( ECC,OMEGA )
      end if

            ! If orbit is circular then simply use the orbital ephemeris
            ! If orbit is eccentric then must calculate the phase diffe-
            ! rence between two successive primary and secondary minima.

      if ( abs(E) < TINY ) then
        GETMIN = TZERO  +  PERIOD * CYCLE
      else
        THETA = 3.0d0*PI/2.0d0 - W
        TANEE = sqrt( (1.0d0-E) / (1.0d0+E) ) * tan(THETA/2.0d0)
        EE = atan(TANEE) * 2.0d0
        PSEP = (EE - E*sin(EE)) / (2.0d0 * PI)
        if ( PSEP < 0.0d0 ) PSEP = PSEP + 1.0d0

        if ( mod(abs(CYCLE),1.0d0) < TINY ) then       ! primary minimum
          GETMIN = TZERO + PERIOD * CYCLE
        else                                         ! secondary minimum
          GETMIN = TZERO + PERIOD * int(CYCLE)
          if ( CYCLE < 1.0d0 ) GETMIN = GETMIN - PERIOD
          GETMIN = GETMIN + PERIOD*PSEP
        end if
      end if

      END FUNCTION GETMIN
!=======================================================================
      DOUBLEPRECISION FUNCTION GETMODEL (V,LDTYPE,TIME,TYPE,LA,LB)
            ! Output a predicted model value according to the parameters
            ! in array V. Precise meaning of the value depends on DTYPE.
            ! DTYPE=1  it outputs an EBOP magnitude for given time
            ! DTYPE=2  it outputs a light ratio for the given time
            ! DTYPE=3  outputs a time of eclipse for the given =CYCLE=
      implicit none
      real*8 V(22)                  ! IN: Photometric parameters
      integer LDTYPE(2)             ! IN: LD law type for the two stars
      real*8 TIME                   ! IN: The given TIME, PHASE or CYCLE
      integer TYPE                  ! IN: 1, 2 or 3 dep on wanted result
      real*8 LA,LB                  ! OUT: Light produced by each star
      real FMAG,LP,LS               ! LOCAL: LIGHT subroutine output
      real*8 ECC,OMEGA              ! LOCAL: eccentricity, perilongitude
      real*8 GETMIN,GETPHASE          ! FUNCTIONS

      if ( TYPE == 1 ) then
        CALL LIGHT(V,LDTYPE,real(GETPHASE(TIME,V(19),V(20))),FMAG,LP,LS)
        LA = dble(LP)
        LB = dble(LS)
        GETMODEL = dble(FMAG)

      else if ( TYPE == 2 ) then
        CALL LIGHT(V,LDTYPE,real(GETPHASE(TIME,V(19),V(20))),FMAG,LP,LS)
        GETMODEL = dble(LS) / dble(LP)

      else if ( TYPE == 3 ) then
        if ( V(7) > 5.0d0 ) then
          ECC = V(7) - 10.0d0
          OMEGA = V(8)
        else
          ECC = sqrt(V(7)**2 + V(8)**2)
          OMEGA = atan2(V(8),V(7)) * 45.0d0 / atan(1.0d0)
          if ( OMEGA < 0.0d0 ) OMEGA = OMEGA + 360.0d0
        end if
        GETMODEL = GETMIN (V(20),V(19),ECC,OMEGA,TIME)

      else
        GETMODEL = -100.0d0
        print*,"### ERROR: wrong datatype asked for in GETMODEL:",TYPE
        STOP
      end if

      END FUNCTION GETMODEL
!=======================================================================
      SUBROUTINE FITEBOP (DATA,DTYPE,NDATA,NLR,NMIN,V,VARY,LDTYPE,ITER,
     &                                                 CHISQ,VERR,ERROR)
            ! This subroutine calls the  MRQMIN algorithm  (Press et al,
            ! 1992, Numerical recipes in FORTRAN 77, p.678)  which finds
            ! the best-fitting EBOP model for the data using the
            ! Levenberg-Marquardt optimisation method.
            ! Unfortunately I have had to use a COMMON block here to
            ! avoid passing parameters through MRQMIN and related code.
      implicit none
      real*8 DATA(3,99999)                ! IN: Observational data
      integer DTYPE(99999)                ! IN: Types of datapoints
      integer NDATA,NLR,NMIN              ! IN: Numbers of datapoints
      real*8 V(22)                        ! IN: EBOP parameter values
      integer VARY(22)                    ! IN: Whether parameters vary
      integer LDTYPE(2)                   ! IN: Type of LD law for stars
      integer ITER                        ! OUT: Number of iterations
      real*8 CHISQ                        ! OUT: Reduced chi-squared
      real*8 VERR(22)                     ! OUT: Formal parameter errors
      integer ERROR                       ! OUT: Whether fit successful
      real*8 X(99999),Y(99999),SIG(99999) ! LOCAL: data to go intoMRQMIN
      real*8 LAMBDA                       ! LOCAL: Marquardt lambda
      real*8 COVAR(22,22),ALPHA(22,22)    ! LOCAL: Covariance etc matrix
      real*8 OCHISQ                       ! LOCAL: Previous chi-squared
      integer i,j,NVARY                   ! LOCAL: helpful integers
      integer CLDTYPE(2),CTYPE(99999)     ! COMMON: for EBOP subroutine

      common / FOREBOP / CLDTYPE,CTYPE

      CLDTYPE(1) = LDTYPE(1)
      CLDTYPE(2) = LDTYPE(2)
      do i = 1,NDATA
        CTYPE(i) = DTYPE(i)
      end do

      do i = 1,NDATA
        X(i)   = DATA(1,i)
        Y(i)   = DATA(2,i)
        SIG(i) = abs(DATA(3,i))
      end do

      NVARY = 0
      do i = 1,22
        if ( VARY(i) == 1 ) NVARY = NVARY + 1
      end do

            ! Now find the best fit using MRQMIN. This requires an init-
            ! ial call with LAMBDA less than zero to initialise things.
            ! Then iterate to the best fit and assume this has been got
            ! once LAMBDA > 10^7. If no observational errors have been
            ! supplied then calculate them and send in once more. And
            ! finally set LAMBDA = 0.0 to get useful parameters out.

105      LAMBDA = -1.0d0
      CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,22,COVAR,ALPHA,22,CHISQ,
     &                                                     LAMBDA,ERROR)
      OCHISQ = 1.0d10
      do i = 1,200
        if ( V(6) > 89.9d0 .and. VARY(6) /= 0) V(6) = 89.9d0
        CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,22,COVAR,ALPHA,22,CHISQ,
     &                                                     LAMBDA,ERROR)
        if ( LAMBDA >= 1.0d8 .and. abs(OCHISQ/CHISQ) < 1.001 ) exit
        if ( ERROR /= 0 ) exit
        OCHISQ = CHISQ
      end do
      ITER = i

      LAMBDA = 0.0d0
      CALL MRQMIN (X,Y,SIG,NDATA,V,VARY,22,COVAR,ALPHA,22,CHISQ,
     &                                                     LAMBDA,ERROR)

            ! Now record the formal errors outputted by MRQMIN

      do i = 1,22
        VERR(i) = sqrt(COVAR(i,i))
        if (DATA(3,1) < 0.0d0) VERR(i)=VERR(i)*sqrt(CHISQ/(NDATA-NVARY))
      end do

      if ( V(6) > 90.0d0 ) V(6) = 180.0d0 - V(6)

      if ( LDTYPE(2) == 0 ) then
        V(5) = V(4)
        VERR(5) = VERR(4)
        V(22) = V(21)
        VERR(22) = VERR(21)
      end if

      if ( DATA(3,1) > 0.0d0 )  CHISQ = CHISQ / (NDATA-NVARY)
      if ( DATA(3,1) <= 0.0d0 ) CHISQ = -1.0d0

      END SUBROUTINE FITEBOP
!=======================================================================
      SUBROUTINE EBOP (INDEX,X,V,Y,DYDA,NCOEFFS,VARY)
            ! This evaluates the goodness of fit and the numerical deri-
            ! vatives for each coefficient for one datapoint using EBOP.
      implicit none
      integer INDEX                 ! IN: Which data point being studied
      real*8 X                      ! IN: Time to  calculate results for
      real*8 V(22)                  ! IN: The   photometric   parameters
      integer NCOEFFS               ! IN: Total number of adjustd params
      integer VARY(22)              ! IN:  Which params  being  adjusted
      real*8 Y                      ! OUT:  Output result for input time
      real*8 DYDA(22)               ! OUT:   The numerical differentials
      real*8 LP,LS                  ! LOCAL: Light produced by each star
      real*8 OUT1,OUT2              ! LOCAL: Help in finding derivatives
      real*8 STORE                  ! LOCAL: Help in finding derivatives
      real*8 DV(22),DVV             ! LOCAL: Amount to perturb params by
      integer i,j,k,ERROR           ! LOCAL: Loop counters and errorflag
      real*8 GETMODEL               ! FUNCTION: Returns model evaluation
      integer LDTYPE(2)             ! IN/COMMON: LD law types  for stars
      integer DTYPE(99999)          ! IN/COMMON: number of LC datapoints

      common / FOREBOP / LDTYPE,DTYPE

      data DV/ 0.01d0,0.01d0,0.01d0,-0.05d0,-0.05d0,-0.1d0,-0.001d0,
     &         -0.001d0,0.1d0,0.1d0,0.01d0,0.01d0,0.01d0,-1.0d0,-0.01d0,
     &           -0.0001d0,-0.01d0,0.5d0,0.1d-6,-0.001d0,0.01d0,0.01d0 /

            ! First get the model prediction for this datapoint. And use
            ! this call to GETMODEL to get the  reflection coefficients.

      Y = GETMODEL (V,LDTYPE,X,DTYPE(INDEX),LP,LS)
      if (VARY(11)==-1)  V(11) = 0.4d0*LS*(V(2)/(1.0d0+V(3)))**2
      if (VARY(12)==-1)  V(12) = 0.4d0*LP*(V(2)/(1.0d0+(1.0d0/V(3))))**2

            ! Now for each adjustable parameter work out the adjustment
            ! not make for calculating the partial derivative

      do i = 1,NCOEFFS
        if ( VARY(i) == 1) then

          if ( DV(i) >= 0.0d0 ) then
            DVV = abs(DV(i) * V(i))
          else
            DVV = abs(DV(i))
          end if

          STORE = V(i)
          V(i) = STORE + DVV
          OUT1 = GETMODEL (V,LDTYPE,X,DTYPE(INDEX),LP,LS)
          V(i) = STORE - DVV
          OUT2 = GETMODEL (V,LDTYPE,X,DTYPE(INDEX),LP,LS)
          V(i) = STORE
          DYDA(i) = (OUT1 - OUT2) / (2.0d0 * DVV)

        else
          DYDA(i) = 0.0d0
        end if
      end do

      if ( V(6) > 90.0d0) then
        DYDA(6) = DYDA(6) + sign( (exp(V(6)-90.0d0))**20 , DYDA(6))
!        DYDA(6) = DYDA(6)*exp(V(6)-90.0d0))**20
!        DYDA(6) = DYDA(6)*exp(V(6)-90.0d0))**2
      end if

      END SUBROUTINE EBOP
!=======================================================================
! !=======================================================================
!=======================================================================
      SUBROUTINE BIAX (R,Q,A,B,EPS)
            ! EBOP subroutine to calculate biaxial ellipsoid dimensions
            ! and oblateness for each star after Chandrasekhar (1933).

      if ( Q <= 0.0 )  then
        A = R
        B = R
        EPS = 0.0
      else
        A = R * ( 1.0 + (1.0 + 7.0*Q)/6.0 * R**3.0)
        B = R * ( 1.0 + (1.0 - 2.0*Q)/6.0 * R**3.0)
        EPS = (A - B) / A
        B=( (1.0 - EPS) * R**3.0) ** (1.0/3.0)
        A = B / (1.0 - EPS)
      end if

      END SUBROUTINE BIAX
!=======================================================================
      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
            ! EBOP subroutine to calculate e and w from e(cos)w e(sin)w

      if ( ECOSW == 0.0  .and.  ESINW == 0.0 ) then
        E = 0.0
        W = 0.0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0 / 3.1415926536
      end if

      END SUBROUTINE GETEW
!=======================================================================
      SUBROUTINE LIGHT (V,LDTYPE,PHASE,FMAG,LP,LS)
      real*8 V(22)
      REAL LP,LS,LECL,LE
      real LD1U,LD2U          ! linear LD coeff for each star
      real LD1Q,LD2Q          ! quadratic LD coeff for each star
      real LD1S,LD2S          ! square-root LD coeff for each star
      real LD1L,LD2L          ! logarithmic LD coeff for each star
      real LDU,LDQ,LDS,LDL    ! LD coeffs for the star in question
      integer LDTYPE(2)       ! LD law type for both stars
      integer GIMENEZ         ! 1 to use original FMAX calculations
                              ! 2 to use Gimenez' modified calcs
                              ! 3 to use Gimenez' calcs for nonlinear LD
      data GIMENEZ / 3 /

      data PI,TWOPI,RAD / 3.1415926536E0,6.28318531E0,0.0174532925E0 /
C
C        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
C        USE SPERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
C        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
C        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
C
      BS     = real(V(1))
      RP     = real(V(2)/(1.0d0+V(3)))
      RATIO  = real(V(3))
      FI     = real(V(6))

      if ( V(7) > 5.0d0 ) then
        ECOSW = real( (V(7)-10.0d0) * cos(V(8)/57.2957795d0) )
        ESINW = real( (V(7)-10.0d0) * sin(V(8)/57.2957795d0) )
      else
        ECOSW  = real(V( 7))
        ESINW  = real(V( 8))
      end if

      YP     = real(V( 9))
      YS     = real(V(10))
      SP     = real(V(11))
      SS     = real(V(12))
      Q      = real(V(13))
      TANGL  = real(V(14))
      EL     = real(V(15))
      DPH    = 1.0 - real(V(16))
      SFACT  = real(V(17))!10.0d0**(-0.4*V(17)))
      DGAM   = real(V(18))

      if ( Q <= 0.0 ) then
        CALL BIAX (RP,0.0,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,0.0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,1.0/Q,RSA,RSB,ES)
      end if

      LD1U = real(V(4))         ! linear term
      LD2U = real(V(5))
      LD1L = 0.0                ! log term
      LD2L = 0.0
      LD1S = 0.0                ! sqrt term
      LD2S = 0.0
      LD1Q = 0.                 ! quadratic term
      LD2Q = 0.0
      if ( LDTYPE(1) == 2 ) LD1L = real(V(21))
      if ( LDTYPE(1) == 3 ) LD1S = real(V(21))
      if ( LDTYPE(1) == 4 ) LD1Q = real(V(21))
      if ( LDTYPE(2) == 2 ) LD2L = real(V(22))
      if ( LDTYPE(2) == 3 ) LD2S = real(V(22))
      if ( LDTYPE(2) == 4 ) LD2Q = real(V(22))
      if ( LDTYPE(2) == 0 ) then
        LD2U = LD1U
        LD2L = LD1L
        LD2S = LD1S
        LD2Q = LD1Q
      end if

!        write(23,*),ldtype
!        write(23,*),ld1u,ld2u
!        write(23,*),ld1l,ld2l
!        write(23,*),ld1s,ld2s
!        write(23,*),ld1q,ld2q
!        write(23,*)," "

C
C        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
C
      THETA=PHASE+DPH
C
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0E0  - SINI2
C
C        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
C
C     EQUATION 9
C        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
C
C        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
C
C        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0E0
      GO TO 25
C
C        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
C        (NON-ZERO ECCENTRICITY ONLY . . . )
C
C     EQUATION 6
C
   17 OMEGA = 450.0E0  - W
   23 IF (OMEGA - 360.0E0)         22,21,21
   21 OMEGA = OMEGA - 360.0E0
      GO TO 23
   22 OMEGA = OMEGA*RAD
C        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
C
C        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
C        CORRESPONDING TO V=OMEGA=90-W
C        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0E0-E*E)*SINW,COSW+E)
C
C        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
C
C        MEAN ANOMALY
      FMA=FMN+FMA0
C     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
C
      DO 10 J=1,15
C        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0E0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
C        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0E-05)     15,15,10
   10 CONTINUE
C
C
C        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0E0-E*E)/DENOM
C
C        RADIUS VECTOR
      RV = (1.0E0-E*E)/(1.0E0+E*COSV)
C
C        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
C        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
C
   25 COS2=COSVW*COSVW
      SIN2=1.0E0-COS2
C
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
C
C
C        PHOTOMETRIC EFFECTS
C
C

!-----------------------------------------------------------------------
! Alvaro Gimenez has corrected the treatment of stellar shapes.
! This newer treatment can be used by putting GIMENEZ=2
! His teatment for nonlinear LD can be used by putting GIMENEZ=3
!-----------------------------------------------------------------------
! This whole thing affects only the brightness normalisation of the two
! eclipsing stars: any problems here wil affect the radiative parameters
! but not the geometric parameters (radii, inclination etc).
!-----------------------------------------------------------------------
      FMAXP = 0.0
      FMAXS = 0.0
      DELTP = 0.0
      DELTS = 0.0
      SHORT = 0.0
      if ( GIMENEZ == 1 ) then

!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))  ! Original
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)                 ! lines
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))  ! if
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)                 ! stars
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP     ! are
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES     ! oblate
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0              ! Original
!       FMAXS=1.0E0-US/3.0E0              ! lines if
!       DELTP=0.0E0                       ! stars
!       DELTS=0.0E0                       ! are
!       SHORT=0.0                         ! spherical

        if ( Q >= 0.0 ) then
          FMAXP=((1.0E0-LD1U)+0.666666667E0*LD1U*(1.0E0+0.2E0*EP))
     1        *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
          FMAXS=((1.0E0-LD2U)+0.666666667E0*LD2U*(1.0E0+0.2E0*ES))
     1        *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
          DELTP=(15.0E0+LD1U)/(15.0E0-5.0E0*LD1U)*(1.0E0+YP)*EP
          DELTS=(15.0E0+LD2U)/(15.0E0-5.0E0*LD2U)*(1.0E0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0
          FMAXS=1.0E0-LD2U/3.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ == 2 ) then

!       FMAXP=(1.0E0-UP*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP
!      1      *(3.0E0-13.0E0/15.0E0*UP))/(1.0E0-EP)
!       FMAXS=(1.0E0-US*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES
!      1      *(3.0E0-13.0E0/15.0E0*US))/(1.0E0-ES)
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
!       SHORT=SINI2*CSVWT*CSVWT

!    26 FMAXP=1.0E0-UP/3.0E0              ! Original
!       FMAXS=1.0E0-US/3.0E0              ! lines if
!       DELTP=0.0E0                       ! stars
!       DELTS=0.0E0                       ! are
!       SHORT=0.0                         ! spherical

        if ( Q >= 0.0 ) then
          FMAXP=(1.0E0-LD1U*(1.0E0-2.0E0/5.0E0*EP)/3.0E0+YP*EP
     1          *(3.0E0-13.0E0/15.0E0*LD1U))/(1.0E0-EP)
          FMAXS=(1.0E0-LD2U*(1.0E0-2.0E0/5.0E0*ES)/3.0E0+YS*ES
     1          *(3.0E0-13.0E0/15.0E0*LD2U))/(1.0E0-ES)
          DELTP=(15.0E0+LD1U)/(15.0E0-5.0E0*LD1U)*(1.0E0+YP)*EP
          DELTS=(15.0E0+LD2U)/(15.0E0-5.0E0*LD2U)*(1.0E0+YS)*ES
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0
          FMAXS=1.0E0-LD2U/3.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!-----------------------------------------------------------------------

      else if ( GIMENEZ == 3 ) then
! And this is Gimenez's code for including nonlinear LD. He includes
! the linear (UP), quadratic (UP, U2P) and square-root (UP, U3P) laws.

!      FMAXP=1.0E0-UP*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.5E0-13.0E0*UP/30.0E0-U2P/5.0E0-23.0E0*U3P/90.0E0)
!      FMAXP=FMAXP/(1.0E0-EP)
!      FMINP=1.0E0-UP*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-
!     1      U2P*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-
!     1      U3P*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
!     1      *(1.0E0-7.0E0*UP/15.0E0-4.0E0*U2P/15.0E0-13.0E0*U3P/45.0E0)
!      FMINS=1.0E0-US*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.0E0-7.0E0*US/15.0E0-4.0E0*U2S/15.0E0-13.0E0*U3S/45.0E0)
!      FMAXS=1.0E0-US*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-
!     1      U2S*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-
!     1      U3S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
!     1      *(1.5E0-13.0E0*US/30.0E0-U2S/5.0E0-23.0E0*U3S/90.0E0)
!      FMAXS=FMAXS/(1.0E0-ES)
!      DELTP=1.0E0-FMINP/FMAXP
!      DELTS=1.0E0-FMINS/FMAXS
!      SHORT=SINI2*CSVWT*CSVWT

!   26 FMAXP=1.0E0-UP/3.0E0-U2P/6.0E0-U3P/5.0E0
!      FMAXS=1.0E0-US/3.0E0-U2S/6.0E0-U3S/5.0E0
!      DELTP=0.0E0
!      DELTS=0.0E0
!      SHORT=0.0

        if ( Q >= 0.0d0 ) then
          FMAXP=1.0E0-LD1U*(1.0E0-2.0E0*EP/5.0E0)/3.0E0-
     1          LD1Q*(1.0E0-3.0E0*EP/5.0E0)/6.0E0-
     1          LD1S*(1.0E0-4.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
     1         *(1.5E0-13.0E0*LD1U/30.0E0-LD1Q/5.0E0-23.0E0*LD1S/90.0E0)
          FMAXP=FMAXP/(1.0E0-EP)
          FMINP=1.0E0-LD1U*(1.0E0+4.0E0*EP/5.0E0)/3.0E0-
     1          LD1Q*(1.0E0+6.0E0*EP/5.0E0)/6.0E0-
     1          LD1S*(1.0E0+8.0E0*EP/9.0E0)/5.0E0+2.0E0*YP*EP
     1   *(1.0E0-7.0E0*LD1U/15.0E0-4.0E0*LD1Q/15.0E0-13.0E0*LD1S/45.0E0)
          FMINS=1.0E0-LD2U*(1.0E0+4.0E0*ES/5.0E0)/3.0E0-
     1          LD2Q*(1.0E0+6.0E0*ES/5.0E0)/6.0E0-
     1          LD2S*(1.0E0+8.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
     1   *(1.0E0-7.0E0*LD2U/15.0E0-4.0E0*LD2Q/15.0E0-13.0E0*LD2S/45.0E0)
          FMAXS=1.0E0-LD2U*(1.0E0-2.0E0*ES/5.0E0)/3.0E0-
     1          LD2Q*(1.0E0-3.0E0*ES/5.0E0)/6.0E0-
     1          LD2S*(1.0E0-4.0E0*ES/9.0E0)/5.0E0+2.0E0*YS*ES
     1         *(1.5E0-13.0E0*LD2U/30.0E0-LD2Q/5.0E0-23.0E0*LD2S/90.0E0)
          FMAXS=FMAXS/(1.0E0-ES)
          DELTP=1.0E0-FMINP/FMAXP
          DELTS=1.0E0-FMINS/FMAXS
          SHORT=SINI2*CSVWT*CSVWT
        else
          FMAXP=1.0E0-LD1U/3.0E0-LD1Q/6.0E0-LD1S/5.0E0
          FMAXS=1.0E0-LD2U/3.0E0-LD2Q/6.0E0-LD2S/5.0E0
          DELTP=0.0E0
          DELTS=0.0E0
          SHORT=0.0
        end if
!----------------------------------------------------------------------
      end if
!----------------------------------------------------------------------
! Complete original code before the above messing:
! C
! C
! C        PHOTOMETRIC EFFECTS
! C
! C
! C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
!       IF (EP .EQ. 0.  .AND.  ES .EQ. 0.)   GO TO 26
! C
! C        EITHER OR BOTH STARS ARE OBLATE
! C
!       FMAXP=((1.0E0-UP)+0.666666667E0*UP*(1.0E0+0.2E0*EP))
!      1      *(1.0E0+3.0E0*YP*EP)/(1.0E0-EP)
!       FMAXS=((1.0E0-US)+0.666666667E0*US*(1.0E0+0.2E0*ES))
!      1      *(1.0E0+3.0E0*YS*ES)/(1.0E0-ES)
! C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
! C        FROM QUADRATURE TO MINIMUM
! C        FACE ON TO END ON
!       DELTP=(15.0E0+UP)/(15.0E0-5.0E0*UP)*(1.0E0+YP)*EP
!       DELTS=(15.0E0+US)/(15.0E0-5.0E0*US)*(1.0E0+YS)*ES
! C        FORE-SHORTENING FUNCTION OF OBLATENESS
!       SHORT=SINI2*CSVWT*CSVWT
!       GO TO 27
! C
! C        BOTH STARS ARE SPHERICAL
! C
!    26 FMAXP=1.0E0-UP/3.0E0
!       FMAXS=1.0E0-US/3.0E0
!       DELTP=0.0E0
!       DELTS=0.0E0
!       SHORT=0.0
!----------------------------------------------------------------------

C
C        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
C        THE NORMALIZING FACTOR
      OTOT=OP+OS
C        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0E0-DELTP*SHORT)
      LS=OS/OTOT*(1.0E0-DELTS*SHORT)
C
C        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0E0  .AND.  SS .EQ. 0.0E0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5E0+0.5E0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0E0
      DLS=0.0E0
C
C        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
C
C     PRIMARY ECLIPSE
C
   30 R1 = RP
      R2 = RS
!----------------------------------------------------------------------!
! JKT mod (10/8/2006): the line these replaced was      UU = UP        !
!----------------------------------------------------------------------!
      LDU = LD1U                                                       !
      LDL = LD1L                                                       !
      LDS = LD1S                                                       !
      LDQ = LD1Q                                                       !
!----------------------------------------------------------------------!
      LE=LP
      DLE=DLP
      GO TO 60
C
C
C     SECONDARY ECLIPSE
C
   40 R1 = RS
      R2 = RP
!-----------------------------------------------------------------------
! JKT mod (10/8/2006): the line these replaced was      UU = US        !
!----------------------------------------------------------------------!
      LDU = LD2U                                                       !
      LDL = LD2L                                                       !
      LDS = LD2S                                                       !
      LDQ = LD2Q                                                       !
!----------------------------------------------------------------------!
      LE=LS
      DLE=DLS
C
   60 SUM = 0.0E0
      ALAST = 0.0E0
      AREA=0.0E0
C
C     EQUATION  5
C
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
      IF (DD .LE. 1.0E-06)  DD=0.0
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
C
C     EQUATION 17
C
      GAMN = 90.01E0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0E0
      RK = 0.0E0
      GAM = 0.0E0
   50 GAM = GAM + DGAMA
C        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
C
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
C
      AA = 0.0E0
C        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
C        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUM = 0.0E0
      GO TO 49
C        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
C
C     EQUATION 12
C
  270 S1 = ABS((R12 - R22 - DD)*0.5E0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)
     1   - (R2-A1)*SQRT(2.0E0*R2*A1-A1*A1))
     2   +R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0E0*RR*B2-B2*B2)
      GO TO 215
C
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0E0
      GO TO 260
C
  230 AA = PI*R12
      GO TO 215
C
C     EQUATION 10
C
  205 S = ABS((R12 - R22 + DD)*0.5E0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0E0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0E0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
C
C     EQUATION 1
C
  210 S = ABS((R12 - R22 + DD)*0.5E0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0E0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0E0*R2*B - B*B)
      AA = A1 + AA1
C
  215 DAREA = AA - ALAST
!----------------------------------------------------------------------!
! JKT modification (10/9/2006). The removed line was:                  !
!     SUM = SUM + DAREA*(1.0E0  - UU + UU*COS(GAM-DGM))                !
!----------------------------------------------------------------------!
      COSGAM = cos(GAM-DGM)                                            !
      SUM = SUM + DAREA*(1.0 - LDU*(1.0-COSGAM) - LDL*COSGAM*log(COSGAM)!
     &          - LDS*(1.0-sqrt(COSGAM))  - LDQ*(1.0-COSGAM)**2)       !
!----------------------------------------------------------------------!
      ALAST = AA
      AREA = AREA + DAREA
C
      RK = RR
      GO TO 50
C
C        LIGHT LOSS FROM ECLIPSE
C
   49 ADISK = PI*R1*R1
!----------------------------------------------------------------------!
! JKT modification (10/9/2006).  See 1992A+A...259..227D for more info.!
! The factor 10.3616329185 was calculated independently by JKT         !
! The removed line was:           ALPHA = SUM/(ADISK*(1.0E0-UU/3.0E0)) !
!----------------------------------------------------------------------!
      ALPHA = 1.0 - LDU/3.0 + LDL/10.3616329185 - LDS/5.0 - LDQ/6.0    !
      ALPHA = SUM/(ADISK*ALPHA)                                        !
!----------------------------------------------------------------------!
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
C
C        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
C        SCALE FACTOR APPLIED
C
C      FLITE = ((LP+LS-LECL+REFL)*(1.0E0-EL)+EL)*SFACT

      FLITE = ((LP+LS-LECL+REFL)*(1.0E0-EL)+EL)
      FMAG = -2.5d0 * log10(FLITE) + SFACT
C
      RETURN
      END
!=======================================================================
!=======================================================================
!=================     NUMERICAL RECIPES SUBROUTINES     ===============
!=======================================================================
!=======================================================================
      SUBROUTINE MRQMIN (x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     &                                                     alamda,ifail)
      implicit none
      integer NDATA,MA,NCA          ! IN: NDATA, numcoeffs, maxcoeffs
      real*8 X(ndata),Y(ndata)      ! IN: data to be fitted
      real*8 SIG(ndata)             ! IN: data errors in y
      real*8 A(ma)                  ! IN: coefficients
      integer IA(ma)                ! IN: adjust (1) or fix (0) coeffs
      real*8 COVAR(nca,nca)         ! OUT: curvature matrix
      real*8 ALPHA(nca,nca)         ! OUT: covariance matrix
      real*8 ALAMDA                 ! IN/OUT: Marquardt lambda factor
      real*8 CHISQ                  ! OUT: chi-squared of the fit

      integer MMAX,j,k,l,m,mfit,ifail
      parameter (MMAX = 22)
      real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.0d0)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).eq.1) mfit=mfit+1
11      continue
        alamda=0.0010d0
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).eq.1) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.0d0+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1,ifail)
      if ( ifail /= 0 ) return
      if(alamda.eq.0.0d0)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).eq.1) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
      if(chisq.lt.ochisq)then
        alamda=0.1d0*alamda
        ochisq=chisq
        j=0
        do 17 l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.0d0*alamda
        chisq=ochisq
      endif
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq)
      IMPLICIT NONE
      integer ma,nalp,ndata,ia(ma),MMAX,mfit,i,j,k,l,m
      parameter ( MMAX = 22 )
      real*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata)
      real*8 y(ndata),dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).eq.1) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.0d0
12      continue
        beta(j)=0.0d0
13    continue
      chisq=0.0d0
      do 16 i=1,ndata

        CALL EBOP (i,x(i),a,ymod,dyda,ma,ia)

        sig2i=1.0d0/(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).eq.1) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).eq.1) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE GAUSSJ (a,n,np,b,m,mp,ifail)
      implicit none
      integer m,mp,n,np,NMAX,ifail
      real*8 a(np,np),b(np,mp),big,dum,pivinv
      parameter (NMAX=22)
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      ifail=0
      irow=0
      icol=0
      do j=1,n
        ipiv(j)=0
      end do
      do i=1,n
        big=0.0d0
        do j=1,n
          if(ipiv(j).ne.1)then
            do k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                print*,"Singular matrix in gaussj ",k
                ifail = 1
                return
              endif
            end do
          endif
        end do
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          end do
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          end do
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.0d0) then
          print*,"singular matrix in gaussj ",ICOL
          ifail = 1
          return
        end if
        pivinv=1.0d0/a(icol,icol)
        a(icol,icol)=1.0d0
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        end do
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        end do
        do ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            end do
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            end do
          endif
        end do
      end do
      do l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          end do
        endif
      end do
      return
      END
!-----------------------------------------------------------------------
      SUBROUTINE COVSRT (covar,npc,ma,ia,mfit)
      implicit none
      integer ma,mfit,npc,ia(ma),i,j,k
      real*8 covar(npc,npc),swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.0d0
          covar(j,i)=0.0d0
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).eq.1)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
!=======================================================================
!=======================================================================
!======================================================================
      INTEGER FUNCTION SEEDSTART ()
            ! This uses the Fortran-intrinsic function SYSTEM_CLOCK to
            ! generate an integer SEED dependent on real time.
            ! SYSTEM_CLOCK returns the present time, the number of
            ! times a second this changes (100 on my computer) and the
            ! maximum value of present time (after which present time
            ! starts again from zero) (2147483647 on my computer).
            ! SEED is outputted as a four-figure integer.
      implicit none
      integer a,b                   ! unused constant values

      CALL SYSTEM_CLOCK (SEEDSTART,a,b)
                        ! SEEDSTART becomes the current clock count

      END FUNCTION SEEDSTART
!=====================================================================
      DOUBLEPRECISION FUNCTION RANDOM (SEED)
            ! Minimal random number generator from Park and Miller.
            ! Return a random deviate between 0.0 and 1.0 with a flat
            ! ditribution. Its period is (2**31 - 1) = 2.15e9
            ! The constants are the best found by Park and Miller.
            ! Set and reset SEED to any integer value to inititalize
            ! the sequence and do not change it thereafter - the code
            ! updates SEED itself.
            ! Does not work if SEED=0, unless it bit-swaps.

      implicit none
      integer SEED                  ! All-important seed value
      integer IA,IM,IR,IQ,MASK      ! Various constants
      real*8 AM                     ! Another constant
      integer K                     ! Yet another constant

      IA = 16807
      IM = 2147483647               ! 2**31 ie largest integer the ..
      AM = 1.0d0 / IM               ! .. computer can handle.
      IQ = 127773
      IR = 2836
      MASK = 123459876

      SEED = IEOR(SEED,MASK)        ! IEOR is integer exclusive-or
      K = SEED / IQ                 ! .. bit-swapping, the non-
      SEED = IA * (SEED - K*IQ) - K*IR    ! .. numeric operation
      if (SEED < 0.0d0) SEED = SEED + IM  ! .. needed by all random
      RANDOM = AM * SEED                  ! .. number generstors.
      SEED = IEOR(SEED,MASK)        ! Seed is updated here

      END FUNCTION RANDOM
!=======================================================================
      DOUBLEPRECISION FUNCTION RANDOMG (SEED,MEAN,SD)
            ! This produces a random number with a Gaussian distribut-
            ! ion. It uses RANDG to generate a random number distrib
            ! with zero mean and unity variance then scales the result
            ! with SD and offsets it with MEAN to produce the effect.
      implicit none
      integer SEED                  ! Seeding integer
      real*8 MEAN,SD                ! Desired mean and S.D. (variance)
      integer ISET                  ! LOCAL: integer variable
      real*8 FAC,GSET,RSQ,V1,V2     ! LOCAL: real variables
      real*8 RANDOM                 ! FUNCTION: flat-distrib random number generator
      SAVE ISET,GSET

      if ( SEED < 0 ) ISET = 0      ! Reinitialise
      if ( ISET == 0 ) then
        do
          V1 = 2.0d0 * RANDOM(SEED) - 1.0d0
          V2 = 2.0d0 * RANDOM(SEED) - 1.0d0
          RSQ = V1**2 + V2**2
          if ( RSQ /= 0.0d0 .and. RSQ <= 1.0d0 ) exit
        end do
        FAC = SQRT( -2.0d0 * log(RSQ) / RSQ )
        GSET = V1 * FAC
        RANDOMG = V2 * FAC
        ISET = 1
      else
        RANDOMG = GSET
        ISET = 0
      end if
      RANDOMG = ( RANDOMG * SD ) + MEAN

      END FUNCTION RANDOMG
!=====================================================================
      DOUBLEPRECISION FUNCTION SELECT (ARRAY,NUM,K)
            ! Returns the Kth smallest value in ARRAY(1:NUM).
            ! ARRAY is simply sorted during this procedure.
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: Input array
      integer K                     ! OUT: value to find
      real*8 STORE                  ! LOCAL: storage variable
      integer i,j,TAG               ! LOCAL: loop counters and a flag

      TAG = 0
      do i = 1,NUM
        STORE = 10.0**10.0
        do j = i,NUM
            if ( ARRAY(j) < STORE )  then
            STORE = ARRAY(j)
            TAG = j
          end if
        end do
        ARRAY(TAG) = ARRAY(i)
        ARRAY(i) = STORE
      end do
      SELECT = ARRAY(K)

      END FUNCTION SELECT
!=======================================================================
      DOUBLEPRECISION FUNCTION SIGMA (ARRAY,NUM)
            ! Returns the standard deviation of array ARRAY
      implicit none
      integer NUM                   ! IN: size of the array
      real*8 ARRAY(NUM)             ! IN: array
      real*8 MEAN,VAR,SUM,SUMSQ     ! LOCAL: properties of ARRAY values
      integer i                     ! LOCAL: loop counter

      SIGMA = 0.0d0
      SUM = 0.0d0
      SUMSQ = 0.0d0
      do i = 1,NUM
        SUM = SUM + ARRAY(i)
        SUMSQ = SUMSQ + ARRAY(i)**2
      end do

      MEAN = SUM / NUM
      VAR = (SUMSQ / NUM) - MEAN**2
      if ( VAR > 0.0d0 )  SIGMA = sqrt( VAR )

      END FUNCTION SIGMA
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
