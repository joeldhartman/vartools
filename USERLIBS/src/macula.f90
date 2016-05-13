! ==============================================================================
! 				MACULA version 1.0			       !
! ==============================================================================
!
! AUTHOR: David Kipping
!         Harvard-Smithsonian Center for Astrophysics
!         Please report any problems to: dkipping@cfa.harvard.edu
!         This code is distributed under the GNU General Public License.
!
! CITATION: If using this code, please cite:
!            Kipping, D. M., 2012, 'An Analytic Model for Rotational Modulations
!            in the Photometry of Spotted Stars', MNRAS, accepted
!
! DESCRIPTION: MACULA computes the photometric flux variations due to starspots
!              on the surface of a rotating star. Star is assumed to be
!              spherical and exhibit four-coefficient non-linear limb darkening.
!              Starspots are circular, gray and exhibit four-coefficient
!              non-linear limb darkening too. MACULA allows for multiple 
!              starspots which are assumed to not overlap. Differential
!              rotation and starspot evolution is included. Partial derivatives
!              are provided of the model flux with respect to the model inputs,
!              as well as time. Finally, expected transit depth variations due
!              to unocculted spots is evaluated.
!              
! INPUTS:
!
! ------------------------------------------------------------------------------
! Theta_star(j) = Parameter vector for the star's intrinsic parameters
! ------------------------------------------------------------------------------
! Istar 	= Theta_star(1)		! Inclination of the star [rads]
! Peq 		= Theta_star(2)		! Rot'n period of the star's equator [d]
! kappa2 	= Theta_star(3)		! Quadratic differential rotation coeff
! kappa4 	= Theta_star(4)		! Quartic differential rotation coeff
! c1 		= Theta_star(5)		! 1st of four-coeff stellar LD terms
! c2 		= Theta_star(6)		! 2nd of four-coeff stellar LD terms
! c3 		= Theta_star(7)		! 3rd of four-coeff stellar LD terms
! c4 		= Theta_star(8)		! 4th of four-coeff stellar LD terms
! d1 		= Theta_star(9)		! 1st of four-coeff spot LD terms
! d2	 	= Theta_star(10)	! 2nd of four-coeff spot LD terms
! d3	 	= Theta_star(11)	! 3rd of four-coeff spot LD terms
! d4	 	= Theta_star(12)	! 4th of four-coeff spot LD terms
! ------------------------------------------------------------------------------
! Theta_spot(j,k) = Parameters of the k^th spot
! ------------------------------------------------------------------------------
! Lambda0(k) 	= Theta_spot(1,k)	! Longitude of spot at time tref(k)
! Phi0(k) 	= Theta_spot(2,k)	! Latitude of spot at time tref(k)
! alphamax(k)	= Theta_spot(3,k)	! Angular spot size at time tmax(k)
! fspot(k)	= Theta_spot(4,k)	! Spot-to-star flux contrast of spot k
! tmax(k)	= Theta_spot(5,k)	! Time at which spot k is largest
! life(k)	= Theta_spot(6,k)	! Lifetime of spot k (FWFM) [days]
! ingress(k)	= Theta_spot(7,k)	! Ingress duration of spot k [days]
! egress(k)	= Theta_spot(8,k)	! Egress duration of spot k  [days]
! ------------------------------------------------------------------------------
! Theta_inst(j,m) = Instrumental/nuisance parameters
! ------------------------------------------------------------------------------
! U(m) 		= Theta_inst(1,m)	! Baseline flux level for m^th data set
! B(m) 		= Theta_inst(2,m)	! Blend factor for m^th data set
! ------------------------------------------------------------------------------
! Non-fitted input parameters
! ------------------------------------------------------------------------------
! Tstart(m)				! Time stamp @ start of m^th data series 
! Tend(m)				! Time stamp @ end of m^th data series
! t(i)					! Time stamp of i^th data point [days]
! ndata					! Number of data points
! mmax					! Number of data sets
! Nspot					! Number of starspots (over all times)
! derivatives				! Compute partial derivatives?
! temporal				! Compute temporal derivatives?
! TdeltaV				! Compute transit depth variations?
! ------------------------------------------------------------------------------
!
! OUTPUTS:
!
! Fmod(i)				! Model flux
! deltaratio(i)				! Transit depth variations prediction
! dFmod_star(i,j)			! Partial derivs of Fmod wrt Theta_star
! dFmod_spot(i,j,k)			! Partial derivs of Fmod wrt Theta_spot
! dFmod_inst(i,j,m)			! Partial derivs of Fmod wrt Theta_inst
! dFmoddt(i)				! Partial derivs of Fmod wrt time
!
! ==============================================================================

MODULE maculamod

implicit none

CONTAINS

! ==============================================================================
! ============================= SUBROUTINE: MACULA =============================
! ==============================================================================

SUBROUTINE macula(t,ndata,Nspot,mmax,derivatives,temporal,TdeltaV,& !Controls
                  Theta_star,Theta_spot,Theta_inst,&	  !System parameters
                  Tstart,Tend,&				  !Data start/end times
                  Fmod,dFmod_star,dFmod_spot,dFmod_inst,dFmoddt,deltaratio) BIND(C, name="MACULAFUNC")!Outs
use ISO_C_BINDING
implicit none

 ! === PRE-AMBLE ===
 ! Note - There should be no need to ever change these three parameters.
 INTEGER, PARAMETER :: pstar = 12
 INTEGER, PARAMETER :: pspot = 8
 INTEGER, PARAMETER :: pinst = 2
 ! === INPUTS ===
 REAL(8), DIMENSION(ndata), INTENT(IN) :: t
 REAL(8), DIMENSION(pstar), INTENT(IN) :: Theta_star
 REAL(8), DIMENSION(pspot,Nspot), INTENT(IN) :: Theta_spot
 REAL(8), DIMENSION(pinst,mmax), INTENT(IN) :: Theta_inst
 REAL(8), DIMENSION(mmax), INTENT(IN) :: Tstart, Tend
 LOGICAL, INTENT(IN) :: derivatives ! Calculate partial derivatives?
 LOGICAL, INTENT(IN) :: temporal    ! Calculate temporal derivatives?
 LOGICAL, INTENT(IN) :: TdeltaV     ! Calculate transit depth variations?
 ! === VARIABLES ===
 INTEGER :: i, j, k, l, m, n
 INTEGER, INTENT(IN) :: ndata, mmax, Nspot
 INTEGER :: jmax
 REAL(8), DIMENSION(Nspot) :: tref ! By default, macula will set tref(k)=tmax(k)
 REAL(8), PARAMETER :: pi = 3.141592653589793D0
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: piI = 0.3183098861837907D0
 REAL(8), PARAMETER :: tol = 1.0D-4 ! alpha values below this will be ignored
 COMPLEX(8), PARAMETER :: minusone = -1.0D0
 ! === SECTION 1 ===
 REAL(8), DIMENSION(pstar+pspot*Nspot+pinst*mmax) :: Theta
 ! === SECTION 2 ===
 REAL(8), DIMENSION(5) :: c, d
 REAL(8), DIMENSION(mmax) :: U, B
 REAL(8), DIMENSION(mmax,ndata) :: Box
 REAL(8), DIMENSION(Nspot) :: Phi0, Prot
 REAL(8), DIMENSION(Nspot,ndata) :: Lambda, beta, alpha
 REAL(8), DIMENSION(Nspot,ndata) :: tdel1, tdel2, tdel3, tdel4
 REAL(8), DIMENSION(Nspot) :: alphamax, fspot, tmax, life, ingress, egress
 ! === SECTION 3 ===
 REAL(8), DIMENSION(Nspot,ndata) :: zetaneg, zetapos
 REAL(8), DIMENSION(5,Nspot,ndata) :: Upsilon, w
 COMPLEX(8), DIMENSION(Nspot,ndata) :: Psi, Xi, Acomplex
 REAL(8), DIMENSION(Nspot,ndata) :: q, A
 REAL(8) :: Fab0
 REAL(8), DIMENSION(ndata) :: Fab
 REAL(8), DIMENSION(ndata), INTENT(OUT) :: Fmod, deltaratio
 ! === SECTION 4 ===
 REAL(8), DIMENSION(5,pstar+pspot*Nspot+pinst*mmax) :: dc, dd
 REAL(8), DIMENSION(mmax,pstar+pspot*Nspot+pinst*mmax) :: dU, dB
 REAL(8), DIMENSION(Nspot,pstar+pspot*Nspot+pinst*mmax) :: dfspot
 ! === SECTION 5 ===
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dalpha
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dbeta
 ! === SECTION 6 ===
 COMPLEX(8), DIMENSION(Nspot,ndata) :: epsil, dAdacomplex, dAdbcomplex
 REAL(8), DIMENSION(Nspot,ndata) :: dAda, dAdb
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) ::  dA
 ! === SECTION 7 ===
 REAL(8), DIMENSION(5,Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dUpsilon, dw
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dzetaneg
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dzetapos
 REAL(8), DIMENSION(Nspot,ndata,pstar+pspot*Nspot+pinst*mmax) :: dq
 REAL(8), DIMENSION(mmax,ndata,pstar+pspot*Nspot+pinst*mmax) :: dFtilde
 REAL(8), DIMENSION(ndata,pstar+pspot*Nspot+pinst*mmax) :: dFab, dFmod
 REAL(8), DIMENSION(Nspot,ndata) :: dzetanegda, dzetaposda
 REAL(8), DIMENSION(pstar+pspot*Nspot+pinst*mmax) :: dFab0
 ! === SECTION 8 ===
 REAL(8), DIMENSION(Nspot,ndata) :: dalphadt, dbetadt, dzetanegdt, dzetaposdt
 REAL(8), DIMENSION(Nspot,ndata) :: dqdt, dAdt
 REAL(8), DIMENSION(5,Nspot,ndata) :: dwdt, dUpsilondt
 REAL(8), DIMENSION(mmax,ndata) :: dFtildedt
 REAL(8), DIMENSION(ndata) :: dFabdt
 REAL(8), DIMENSION(ndata), INTENT(OUT) :: dFmoddt
 ! === SECTION 9 ===
 REAL(8), DIMENSION(ndata,pstar), INTENT(OUT) :: dFmod_star
 REAL(8), DIMENSION(ndata,pspot,Nspot), INTENT(OUT) :: dFmod_spot
 REAL(8), DIMENSION(ndata,pinst,mmax), INTENT(OUT) :: dFmod_inst

 ! ===================================
 ! === SECTION 1: THETA ASSIGNMENT ===
 ! ===================================

 jmax = pstar + pspot*Nspot + pinst*mmax

 l = 0
 DO j=1,pstar
   l = l + 1
   Theta(l) = Theta_star(j)
 END DO
 DO k=1,Nspot
   DO j=1,pspot
     l = l + 1
     Theta(l) = Theta_spot(j,k)
   END DO
 END DO
 DO j=1,pinst
   DO m=1,mmax
     l = l + 1
     Theta(l) = Theta_inst(j,m)
   END DO
 END DO

 ! Thus we have...
 !Theta_star(j) = Theta(j)
 !Theta_spot(j,k) = Theta(pstar + pspot*(k-1) + j)
 !Theta_inst(m) = Theta(pstar + pspot*Nspot + m)

 ! ===================================
 ! === SECTION 2: BASIC PARAMETERS ===
 ! ===================================

 ! c and d assignment
 c(2) = Theta(5); c(3) = Theta(6); c(4) = Theta(7); c(5) = Theta(8)
 d(2) = Theta(9); d(3) =Theta(10); d(4)= Theta(11); d(5)= Theta(12)
 c(1) = 1.0D0 - c(2) - c(3) - c(4) - c(5) !c0
 d(1) = 1.0D0 - d(2) - d(3) - d(4) - d(5) !d0

 ! U and B assignment
 DO m=1,mmax
   U(m) = Theta(pstar+pspot*Nspot+m)
   B(m) = Theta(pstar+pspot*Nspot+mmax+m)
 END DO

 ! Box-car function (labelled as Pi_m in the paper)
 DO i=1,ndata
   DO m=1,mmax
     Box(m,i) = heaviside(t(i)-Tstart(m)) - heaviside(t(i)-Tend(m))
   END DO
 END DO

 ! Phi0, Prot & Lambda calculation
 DO k=1,Nspot
   Phi0(k) = Theta(pstar+pspot*(k-1)+2)
 END DO
 DO k=1,Nspot
   Prot(k) =Theta(2)/(1.0D0-Theta(3)*DSIN(Phi0(k))**2-Theta(4)*DSIN(Phi0(k))**4)
 END DO
 DO i=1,ndata
   DO k=1,Nspot
     Lambda(k,i) = Theta(pstar+pspot*(k-1)+1) + 2.0D0*pi*(t(i)-tref(k))/Prot(k)
   END DO
 END DO

 ! beta calculation
 DO i=1,ndata
   DO k=1,Nspot
     beta(k,i) = DACOS( DCOS(Theta(1))*DSIN(Phi0(k)) &
                 + DSIN(Theta(1))*DCOS(Phi0(k))*DCOS(Lambda(k,i)) )
   END DO
 END DO

 ! alpha calculation
 DO k=1,Nspot
   alphamax(k)	= Theta(pstar + pspot*(k-1) + 3)
   fspot(k)	= Theta(pstar + pspot*(k-1) + 4)
   tmax(k)	= Theta(pstar + pspot*(k-1) + 5)
   life(k)	= Theta(pstar + pspot*(k-1) + 6)
   ingress(k)	= Theta(pstar + pspot*(k-1) + 7)
   egress(k)	= Theta(pstar + pspot*(k-1) + 8)
   IF( ingress(k) .LT. 0.0416667D0 ) THEN ! 1 hour; minimum ingress time
      ingress = 0.0416667D0
   END IF
   IF( egress(k) .LT. 0.0416667D0 ) THEN ! 1 hour; minimum egress time
      egress = 0.0416667D0
   END IF
 END DO
 ! macula defines the reference time = maximum spot-size time
 ! However, one can change the line below to whatever they wish.
 tref(:) = tmax(:)
 DO i=1,ndata
   DO k=1,Nspot
     tdel1(k,i) = t(i)-tmax(k)+0.5D0*life(k)+ingress(k)
     tdel2(k,i) = t(i)-tmax(k)+0.5D0*life(k)
     tdel3(k,i) = t(i)-tmax(k)-0.5D0*life(k)
     tdel4(k,i) = t(i)-tmax(k)-0.5D0*life(k)-egress(k)
     alpha(k,i) = (( tdel1(k,i)*heaviside(tdel1(k,i)) &
                  - tdel2(k,i)*heaviside(tdel2(k,i)) )/ingress(k)) - &
                  (( tdel3(k,i)*heaviside(tdel3(k,i)) &
                  - tdel4(k,i)*heaviside(tdel4(k,i)) )/egress(k))
     alpha(k,i) = alphamax(k)*alpha(k,i)
     !write(*,*) life(k),ingress(k),egress(k)
   END DO
 END DO

 ! =================================
 ! === SECTION 3: COMPUTING FMOD ===
 ! =================================

 ! zetapos and zetaneg
 DO i=1,ndata
   DO k=1,Nspot
     zetapos(k,i) = zeta(beta(k,i)+alpha(k,i))
     zetaneg(k,i) = zeta(beta(k,i)-alpha(k,i))
   END DO
 END DO

 ! Upsilon
 DO i=1,ndata
   DO k=1,Nspot
     DO n=0,4
       Upsilon(n+1,k,i) = zetaneg(k,i)**2 - zetapos(k,i)**2 &
                          + kronecker(zetapos(k,i),zetaneg(k,i))
       Upsilon(n+1,k,i) = (DSQRT(zetaneg(k,i)**(n+4)) &
                          - DSQRT(zetapos(k,i)**(n+4)))&
                          /Upsilon(n+1,k,i)
     END DO
   END DO
 END DO
        
 ! w
 DO i=1,ndata
   DO k=1,Nspot
     DO n=0,4
       w(n+1,k,i) = ( 4.0D0*(c(n+1)-d(n+1)*fspot(k)) )/( n+4.0D0 )
       w(n+1,k,i) = w(n+1,k,i)*Upsilon(n+1,k,i)
     END DO
   END DO
 END DO

 ! Area A
 DO i=1,ndata
   DO k=1,Nspot
     ! Psi, Xi & Acomplex are all complex parameters
     Psi(k,i) = SQRT(-minusone - (DCOS(alpha(k,i))/DSIN(beta(k,i)))**2 )
     IF( alpha(k,i) .GT. tol ) THEN
       Xi(k,i) = DSIN(alpha(k,i))*&
                 CACOS(-1.0D0/(DTAN(alpha(k,i))*DTAN(beta(k,i))))
     ELSE
       Xi(k,i) = 0.0D0 ! In the limit of alpha->0, Xi->0
     END IF
     IF( beta(k,i) .LT. (halfpi+alpha(k,i)) .AND. alpha(k,i) .GT. tol ) THEN
       Acomplex(k,i) = CACOS( DCOS(alpha(k,i))/DSIN(beta(k,i)) ) &
                     + DCOS(beta(k,i))*DSIN(alpha(k,i))*Xi(k,i) &
                     - DCOS(alpha(k,i))*DSIN(beta(k,i))*Psi(k,i)
     ELSE
       Acomplex(k,i) = 0.0D0
     END IF
     A(k,i) = REAL(Acomplex(k,i))
     ! This should not be needed, but just incase...
     IF( IsNan(A(k,i)) ) THEN
       A(k,i) = 0.0D0
       !write(*,*) 'ERROR: ',A(k,i),alpha(k,i)
     END IF
   END DO
 END DO

 ! q
 DO i=1,ndata
   DO k=1,Nspot
     q(k,i) = (A(k,i)*piI)*SUM(w(:,k,i))
   END DO
 END DO

 ! Fab0
 Fab0 = 0.0D0
 DO n=0,4
   Fab0 = Fab0 + (n*c(n+1))/(n+4.0D0)
 END DO
 Fab0 = 1.0D0 - Fab0

 ! Fab
 DO i=1,ndata
   Fab(i) = Fab0 - SUM(q(:,i))
 END DO

 ! Fmod
 DO i=1,ndata
   Fmod(i) = 0.0D0
   DO m=1,mmax
     Fmod(i) = Fmod(i) + U(m)*Box(m,i)*( Fab(i)/(Fab0*B(m))+(B(m)-1.0D0)/B(m) )
   END DO
 END DO

 ! delta_{obs}/delta
 IF( TdeltaV ) THEN
   DO i=1,ndata
     deltaratio(i) = 0.0D0
     DO m=1,mmax
       deltaratio(i) = deltaratio(i) + U(m)*Box(m,i)
     END DO
     deltaratio(i) = 2.0D0 - (deltaratio(i)/Fmod(i))
     deltaratio(i) = deltaratio(i)/( 1.0D0 - piI*SUM(A(:,i)) )
   END DO
 ELSE
   deltaratio(:) = 1.0D0
 END IF

 ! ====================================
 ! === SECTION 4: BASIS DERIVATIVES ===
 ! ====================================
 
 ! Master if-loop
 IF( derivatives ) THEN

 ! Derivatives of c & d. Note that c_n = c(n+1) due to memory constraints.
 dc(:,:) = 0.0D0 !dc/dTheta=0 for all n and j unless Theta_j=c_n...
 dc(2,5) = 1.0D0;!d(c1)/d(Theta5) = d(c1)/(c1) = 1
 dc(3,6) = 1.0D0;!d(c2)/d(Theta6) = d(c2)/(c2) = 1
 dc(4,7) = 1.0D0;!d(c3)/d(Theta7) = d(c3)/(c3) = 1
 dc(5,8) = 1.0D0;!d(c4)/d(Theta8) = d(c4)/(c4) = 1
 dc(1,5) =-1.0D0;!d(c0)/d(Theta5) = d(c0)/(c1) =-1
 dc(1,6) =-1.0D0;!d(c0)/d(Theta6) = d(c0)/(c2) =-1
 dc(1,7) =-1.0D0;!d(c0)/d(Theta7) = d(c0)/(c3) =-1
 dc(1,8) =-1.0D0;!d(c0)/d(Theta8) = d(c0)/(c4) =-1
 dd(:,:) = 0.0D0 !dd/dTheta=0 for all n and j unless Theta_j=c_n...
 dd(2,9) = 1.0D0;!d(d1)/d(Theta9) = d(d1)/(d1) = 1
 dd(3,10)= 1.0D0;!d(d2)/d(Theta10)= d(d2)/(d2) = 1
 dd(4,11)= 1.0D0;!d(d3)/d(Theta11)= d(d3)/(d3) = 1
 dd(5,12)= 1.0D0;!d(d4)/d(Theta12)= d(d4)/(d4) = 1
 dd(1,9) =-1.0D0;!d(c0)/d(Theta9) = d(d0)/(d1) =-1
 dd(1,10)=-1.0D0;!d(c0)/d(Theta10)= d(d0)/(d2) =-1
 dd(1,11)=-1.0D0;!d(c0)/d(Theta11)= d(d0)/(d3) =-1
 dd(1,12)=-1.0D0;!d(c0)/d(Theta12)= d(d0)/(d4) =-1

 ! Derivatives of U
 dU(:,:) = 0.0D0
 DO m=1,mmax
   ! d(U_m)/d(Theta_j) for j=pstar+pspot*Nspot+m => Theta_j=U_m
   dU(m,pstar+pspot*Nspot+m) = 1.0D0 
 END DO

 ! Derivatives of B
 dB(:,:) = 0.0D0
 DO m=1,mmax
   ! d(B_m)/d(Theta_j) for j=pstar+pspot*Nspot+m+mmax => Theta_j=B_m
   dB(m,pstar+pspot*Nspot+mmax+m) = 1.0D0
 END DO

 ! Derivatives of fspot (4th spot parameter)
 DO k=1,Nspot
   dfspot(k,:) = 0.0D0
   !d(fspot_k)/d(Theta_j) for j=pstar+pspot*(k-1)+4 => Theta_j=fspot_k
   dfspot(k,pstar+pspot*(k-1)+4) = 1.0D0
 END DO

 ! ===========================================
 ! === SECTION 5: ALPHA & BETA DERIVATIVES ===
 ! ===========================================

 ! Derivatives of alpha(alphamax,tmax,life,ingress,egress) 
 ! [function of 4*Nspot parameters]
 DO i=1,ndata
   DO k=1,Nspot
     dalpha(k,i,:) = 0.0D0 
     ! wrt alphamax (3rd spot parameter)
     dalpha(k,i,pstar+pspot*(k-1)+3) = alpha(k,i)/alphamax(k)
     ! wrt tmax (5th spot parameter)
     dalpha(k,i,pstar+pspot*(k-1)+5) = -(alphamax(k)/ingress(k))*( &
                                            heaviside(tdel1(k,i)) &
                                          - heaviside(tdel2(k,i)) ) &
                                       + (alphamax(k)/egress(k))*( &
                                            heaviside(tdel3(k,i)) &
                                          - heaviside(tdel4(k,i)) )
     ! wrt life (6th spot parameter)
     dalpha(k,i,pstar+pspot*(k-1)+6) = 0.5D0*(alphamax(k)/ingress(k))&
                                       *( heaviside(tdel1(k,i)) &
                                          - heaviside(tdel2(k,i)) )&
                                       + 0.5D0*(alphamax(k)/egress(k))&
                                       *( heaviside(tdel3(k,i)) &
                                          - heaviside(tdel4(k,i)) )
     ! wrt ingress (7th spot parameter)
     dalpha(k,i,pstar+pspot*(k-1)+7) = -(alphamax(k)/ingress(k)**2)&
                                       *(0.50D0*(tdel1(k,i)+tdel2(k,i)))&
                                       *( heaviside(tdel1(k,i)) - &
                                          heaviside(tdel2(k,i)) )
     ! wrt egress (8th spot parameter)
     dalpha(k,i,pstar+pspot*(k-1)+8) = (alphamax(k)/egress(k)**2)&
                                       *(0.50D0*(tdel3(k,i)+tdel4(k,i)))&
                                       *( heaviside(tdel3(k,i)) - &
                                          heaviside(tdel4(k,i)) )
   END DO
 END DO

 ! Stellar-derivatives of beta(Istar,Phi0,Lambda0,Peq,kappa2,kappa4) 
 ! [function of 4+2*Nspot parameters]
 DO i=1,ndata
   DO k=1,Nspot
     ! pre-amble
     dbeta(k,i,:) = 0.0D0
     ! wrt Istar (1st star parameter)
     dbeta(k,i,1) = DSIN(Phi0(k))*DSIN(Theta(1)) &
                    - DCOS(Lambda(k,i))*DCOS(Phi0(k))*DCOS(Theta(1)) 
     dbeta(k,i,1) = dbeta(k,i,1)/DSIN(beta(k,i))
     ! wrt Peq (2nd star parameter)
     dbeta(k,i,2) = DCOS(Phi0(k))*DSIN(Lambda(k,i))*DSIN(Theta(1))&
                    /DSIN(beta(k,i))
     dbeta(k,i,2) = dbeta(k,i,2)*2.0D0*pi*(t(i)-tref(k))/(Theta(2)) ! Temporary
     ! wrt kappa2 (3rd star parameter)
     dbeta(k,i,3) = -dbeta(k,i,2)*DSIN(Phi0(k))**2
     ! wrt kappa4 (4th star parameter)
     dbeta(k,i,4) = -dbeta(k,i,2)*DSIN(Phi0(k))**4
     ! wrt Peq continued
     dbeta(k,i,2) = -dbeta(k,i,2)/Prot(k)
   END DO
 END DO

 ! Spot-derivatives of beta
 DO i=1,ndata
   DO k=1,Nspot
     ! wrt Lambda [1st spot parameter]
     dbeta(k,i,pstar+pspot*(k-1)+1) = DSIN(Theta(1))*DCOS(Phi0(k))*&
                                      DSIN(Lambda(k,i))/DSIN(beta(k,i))
     ! wrt Phi0 [2nd spot parameter]
     dbeta(k,i,pstar+pspot*(k-1)+2) = 2.0D0*Theta(3)*DCOS(Phi0(k))**2 &
                                      + Theta(4)*DSIN(2.0D0*Phi0(k))**2
     dbeta(k,i,pstar+pspot*(k-1)+2) = dbeta(k,i,pstar+pspot*(k-1)+2)*2.0D0*pi*&
                                      (t(i)-tref(k))/Theta(2)
     dbeta(k,i,pstar+pspot*(k-1)+2) = DCOS(Lambda(k,i)) &
                                      - dbeta(k,i,pstar+pspot*(k-1)+2)
     dbeta(k,i,pstar+pspot*(k-1)+2) = dbeta(k,i,pstar+pspot*(k-1)+2)*&
                                      DSIN(Theta(1))*DSIN(Phi0(k))/&
                                      DSIN(beta(k,i))
   END DO
 END DO

 ! ================================
 ! === SECTION 6: A DERIVATIVES ===
 ! ================================

 ! Semi-derivatives of A
 DO i=1,ndata
   DO k=1,Nspot
     ! epsil, dAdacomplex & dAdacomplex are all complex parameters
     epsil(k,i) = ( DCOS(2.0D0*alpha(k,i)) + DCOS(2.0D0*beta(k,i)) )/( &
                  DSIN(beta(k,i))**2*Psi(k,i) )
     dAdacomplex(k,i) = -DSIN(alpha(k,i))*DSIN(beta(k,i))*epsil(k,i) &
                        +2.0D0*DCOS(alpha(k,i))*DCOS(beta(k,i))*Xi(k,i)
     dAdbcomplex(k,i) = 0.5D0*DCOS(alpha(k,i))*DCOS(beta(k,i))*epsil(k,i) &
                        - DSIN(alpha(k,i))*DSIN(beta(k,i))*Xi(k,i)
     dAda(k,i) = REAL(dAdacomplex(k,i))
     dAdb(k,i) = REAL(dAdbcomplex(k,i))
   END DO
 END DO

 ! Derivatives of A
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       dA(k,i,j) = dAda(k,i)*dalpha(k,i,j) + dAdb(k,i)*dbeta(k,i,j)
     END DO
   END DO
 END DO

 ! ====================================
 ! === SECTION 7: FINAL DERIVATIVES ===
 ! ====================================

 ! Derivatives of zeta wrt alpha (and implicitly beta)
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       dzetanegda(k,i) = heaviside( halfpi - (beta(k,i)-alpha(k,i)) )*&
                         heaviside( beta(k,i) - alpha(k,i) )
       dzetanegda(k,i) = dzetanegda(k,i)*DSIN(beta(k,i)-alpha(k,i))
       dzetaposda(k,i) = -heaviside( halfpi - (beta(k,i)+alpha(k,i)) )*&
                         heaviside( beta(k,i) + alpha(k,i) )
       dzetaposda(k,i) = dzetaposda(k,i)*DSIN(beta(k,i)+alpha(k,i))
     END DO
   END DO
 END DO

 ! Derivatives of zeta
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       dzetaneg(k,i,j) = dzetanegda(k,i)*( dalpha(k,i,j) - dbeta(k,i,j) )
       dzetapos(k,i,j) = dzetaposda(k,i)*( dalpha(k,i,j) + dbeta(k,i,j) )
     END DO
   END DO
 END DO

 ! Derivatives of Upsilon
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       DO n=0,4
         dUpsilon(n+1,k,i,j) = DSQRT( zetaneg(k,i)**(n+2) )*dzetaneg(k,i,j) &
                           - DSQRT( zetapos(k,i)**(n+2) )*dzetapos(k,i,j)
         dUpsilon(n+1,k,i,j) = 0.5D0*(n+4.0D0)*dUpsilon(n+1,k,i,j) &
                               - 2.0D0*Upsilon(n+1,k,i)*&
                               ( dzetaneg(k,i,j)-dzetapos(k,i,j) )
         dUpsilon(n+1,k,i,j) = dUpsilon(n+1,k,i,j)/( zetaneg(k,i)**2 &
                               - zetapos(k,i)**2 &
                               + kronecker(zetapos(k,i),zetaneg(k,i)) )
       END DO
    END DO
   END DO
 END DO

 ! Derivatives of w
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       DO n=0,4
         dw(n+1,k,i,j) = Upsilon(n+1,k,i)*dc(n+1,j) + (c(n+1)-d(n+1)*fspot(k))*&
                         dUpsilon(n+1,k,i,j) - d(n+1)*Upsilon(n+1,k,i)*&
                         dfspot(k,j) - fspot(k)*Upsilon(n+1,k,i)*dd(n+1,j)
         dw(n+1,k,i,j) = (4.0D0*dw(n+1,k,i,j))/(n+4.0D0)
       END DO
     END DO
   END DO
 END DO

 ! Derivatives of q
 DO j=1,jmax
   DO i=1,ndata
     DO k=1,Nspot
       dq(k,i,j) = ( A(k,i)*SUM(dw(:,k,i,j)) + dA(k,i,j)*SUM(w(:,k,i)) )*piI
     END DO
   END DO
 END DO

 ! Derivatives of Fab
 DO j=1,jmax
   DO i=1,ndata
     dFab(i,j) = dFab0(j) - SUM(dq(:,i,j))
   END DO
 END DO

 ! Derivatives of Fab0
 dFab0(:) = 0.0D0               ! wrt everything else
 dFab0(5) = -0.2D0		! wrt c1
 dFab0(6) = -0.33333333D0	! wrt c2
 dFab0(7) = -0.42857143D0	! wrt c3
 dFab0(8) = -0.5D0 		! wrt c4

 ! Derivatives of Ftilde
 DO j=1,jmax
   DO i=1,ndata
     DO m=1,mmax
       dFtilde(m,i,j) = Fab0*B(m)*(Fab(i)+Fab0*(B(m)-1.0D0))*dU(m,j) &
                        + U(m)*( B(m)*Fab0*dFab(i,j) - B(m)*Fab(i)*dFab0(j) &
		        + Fab0*(Fab0-Fab(i))*dB(m,j) )
       dFtilde(m,i,j) = ( Box(m,i)/(Fab0**2*B(m)**2) )*dFtilde(m,i,j)
     END DO
   END DO
 END DO

 ! Derivatives of Fmod
 DO j=1,jmax
   DO i=1,ndata
     dFmod(i,j) = SUM(dFtilde(:,i,j))
   END DO
 END DO

 ! =======================================
 ! === SECTION 8: TEMPORAL DERIVATIVES ===
 ! =======================================

 IF( temporal ) THEN

 ! Temporal derivatives of alpha and beta
 DO i=1,ndata
   DO k=1,Nspot
     dalphadt(k,i) = (alphamax(k)/ingress(k))&
                     *( heaviside(tdel1(k,i)) - heaviside(tdel2(k,i)) ) &
                     - (alphamax(k)/egress(k))&
                     *( heaviside(tdel3(k,i)) - heaviside(tdel4(k,i)) )
     dbetadt(k,i) = 1.0D0 - ( DSIN(Theta(1))*DCOS(Phi0(k))*DCOS(Lambda(k,i)) + &
                    DCOS(Theta(1))*DSIN(Phi0(k)) )**2
     dbetadt(k,i) = (2.0D0*pi*DSIN(Theta(1))*DCOS(Phi0(k))*DSIN(Lambda(k,i)) )/&
                    ( Prot(k)*DSQRT(dbetadt(k,i)) )
   END DO
 END DO

 ! Temporal derivatives of zeta
 DO i=1,ndata
   DO k=1,Nspot
     dzetanegdt(k,i) = dzetanegda(k,i)*( dalphadt(k,i) - dbetadt(k,i) )
     dzetaposdt(k,i) = dzetaposda(k,i)*( dalphadt(k,i) + dbetadt(k,i) )
   END DO
 END DO

 ! Temporal derivatives of Upsilon
 DO i=1,ndata
   DO k=1,Nspot
     DO n=0,4
       dUpsilondt(n+1,k,i) = DSQRT( zetaneg(k,i)**(n+2) )*dzetanegdt(k,i) &
                         - DSQRT( zetapos(k,i)**(n+2) )*dzetaposdt(k,i)
       dUpsilondt(n+1,k,i) = 0.5D0*(n+4.0D0)*dUpsilondt(n+1,k,i) &
                             - 2.0D0*Upsilon(n+1,k,i)*&
                             ( dzetanegdt(k,i)-dzetaposdt(k,i) )
       dUpsilondt(n+1,k,i) = dUpsilondt(n+1,k,i)/( zetaneg(k,i)**2 &
                             - zetapos(k,i)**2 &
                             + kronecker(zetapos(k,i),zetaneg(k,i)) )
     END DO
   END DO
 END DO

 ! Temporal derivatives of w
 DO i=1,ndata
   DO k=1,Nspot
     DO n=0,4
       dwdt(n+1,k,i) = dUpsilondt(n+1,k,i)*&
                       ( ( 4.0D0*( c(n+1) - d(n+1)*fspot(k) ) )/(n+4.0D0) )
     END DO
   END DO
 END DO

 ! Temporal derivatives of A
 DO i=1,ndata
   DO k=1,Nspot
     dAdt(k,i) = dAda(k,i)*dalphadt(k,i) + dAdb(k,i)*dbetadt(k,i)
   END DO
 END DO

 ! Temporal derivatives of q
 DO i=1,ndata
   DO k=1,Nspot
     dqdt(k,i) = ( A(k,i)*SUM(dwdt(:,k,i)) + dAdt(k,i)*SUM(w(:,k,i)) )*piI
   END DO
 END DO

 ! Temporal derivatives of Fab
 DO i=1,ndata
   dFabdt(i) = - SUM(dqdt(:,i))
 END DO

 ! Temporal derivatives of Ftilde
 DO i=1,ndata
   DO m=1,mmax
     dFtildedt(m,i) = ( (U(m)*Box(m,i))/(B(m)*Fab0) )*dFabdt(i)
   END DO
 END DO

 ! Temporal derivatives of Fmod
 DO i=1,ndata
   dFmoddt(i) = SUM(dFtildedt(:,i))
 END DO

 ELSE
   dFmoddt(:) = 0.0D0
 END IF

 ! =======================================
 ! === SECTION 9: RE-SPLIT DERIVATIVES ===
 ! =======================================

 ! Derivatives provided for Theta_star, Theta_inst, Theta_spot discretely
 l=0
 DO j=1,pstar
   l = l + 1
   dFmod_star(:,j) = dFmod(:,l)
 END DO
 DO k=1,Nspot
   DO j=1,pspot
     l = l + 1
     dFmod_spot(:,j,k) = dFmod(:,l)
   END DO
 END DO
 DO m=1,mmax
   DO j=1,pinst
     l = l + 1
     dFmod_inst(:,j,m) = dFmod(:,l)
   END DO
 END DO
 
 ! Master if-loop
 ELSE
   dFmod_star(:,:)   = 0.0D0
   dFmod_spot(:,:,:) = 0.0D0
   dFmod_inst(:,:,:) = 0.0D0
   dFmoddt(:) = 0.0D0
 END IF

END SUBROUTINE macula
! ==============================================================================

! ==============================================================================
! ================================= FUNCTION: ZETA =============================
! ==============================================================================
REAL(8) FUNCTION zeta(x)

implicit none

 REAL(8), INTENT(IN) :: x
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0

 zeta = DCOS(x)*heaviside(x)*heaviside(halfpi-x) + heaviside(-x)

 END FUNCTION zeta
! ==============================================================================

! ==============================================================================
! =========================== FUNCTION: KRONECKER ==============================
! ==============================================================================
REAL(8) FUNCTION kronecker(x,y)

implicit none

 REAL(8), INTENT(IN) :: x,y

 IF( x .EQ. y ) THEN
   kronecker = 1.0D0
 ELSE
  kronecker = 0.0D0
 END IF

 END FUNCTION kronecker
! ==============================================================================

! ==============================================================================
! ============================ FUNCTION: HEAVISIDE =============================
! ==============================================================================
REAL(8) FUNCTION heaviside(x)

implicit none

 REAL(8), INTENT(IN) :: x

 IF( x .LT. 0.0D0 ) THEN
   heaviside = 0.0D0
 ELSE IF( x .EQ. 0.0D0 ) THEN
   heaviside = 0.5D0
 ELSE IF( x .GT. 0.0D0 ) THEN
   heaviside = 1.0D0
 END IF

 END FUNCTION heaviside
! ==============================================================================

! ==============================================================================
! ========================== FUNCTION: COMPLEX ARCCOS ==========================
! ==============================================================================
COMPLEX(8) FUNCTION CACOS(y)

implicit none

 REAL(8), INTENT(IN) :: y
 COMPLEX(8) :: ycomplex
 COMPLEX(8), PARAMETER :: minusone = -1.0D0
 COMPLEX(8), PARAMETER :: halfpi = 1.5707963267948966D0

 ycomplex = y
 CACOS = halfpi + SQRT(minusone)*LOG( SQRT(minusone)*ycomplex &
         + SQRT(1.0D0 - ycomplex**2) )

 END FUNCTION CACOS
! ==============================================================================

END MODULE maculamod
