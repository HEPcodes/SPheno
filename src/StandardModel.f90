Module StandardModel

! load modules
Use Control
Use Mathematics
! load modules


! global variables
 ! Z-boson
 Real(dp), Save :: mZ,mZ2,gamZ,gamZ2,gmZ,gmZ2,BrZqq(5),BrZll(3),BrZinv 
 ! W-boson
 Real(dp), Save :: mW,mW2,gamW,gamW2,gmW,gmW2,BrWqq(2),BrWln(3)
 ! fermion masses
 Real(dp), Save :: mf_l(3), mf_l2(3), mf_u(3), mf_u2(3), mf_d(3), mf_d2(3)  &
   & , mf_nu(3), mf_l_mZ(3), mf_d_mZ(3), mf_u_mZ(3), mf_d_mt(3), mf_u_mt(3) &
   & , mf_l_mt(3), mf_d_125(3), mf_u_125(3), mf_l_125(3)
 ! scale at which the light quark masses for u, d, c, s are defined
 Real(dp), Save :: Q_light_quarks
 ! fermion decay widths
 Real(dp), Save :: g_T = 2.0_dp
 ! meson masses
 Real(dp), Parameter :: m_pi0=0.135_dp, m_pip=0.140_dp 
 Real(dp), Parameter :: m_pi02=m_pi0**2, m_pip2=m_pip**2
 ! decay constants
 Real(dp), Parameter :: f_pi=0.092_dp
 ! couplings
 Real(dp), Save :: Alpha, Alpha_mZ, AlphaS_mZ, G_F, Alpha_mZ_MS, AlphaS_mB &
   & , AlphaS_mT, Alpha_mT, AlphaS_125, Alpha_125
 Real(dp), Save :: Delta_Alpha_Lepton, Delta_Alpha_Hadron
! pseudo observables
 Real(dp), Save :: Rho_parameter
 ! partial widht of mu and tau
 Real(dp), Save :: GammaMu, GammaTau

 Real(dp), Save :: KFactorLee  ! for ISR corrections in e+ e- 

 Complex(dp), Save :: CKM(3,3)
! update to PDG 2013
 Real(dp), Save :: lam_wolf=0.22535_dp, A_wolf=0.811_dp, rho_wolf=0.131_dp &
                 & , eta_wolf=0.345_dp

 Complex(dp), Save :: Unu(3,3)
 Real(dp), Save :: theta_12, theta_13, theta_23, delta_nu, alpha_nu1, alpha_nu2
 ! h-bar
 Real(dp), Save :: hbar = 6.58211928e-25_dp
!-----------------------------------------------------------------
! constants for B-physics, can be changed using the FLHA
! taking data from arXiv:0808.1297v3 and arXiv:0910.2928
! update: http://krone.physik.unizh.ch/~lunghi/webpage/LatAves/page7/page7.html
! update arXiv:1204.0791
! experimental data from PDG 2013
!-----------------------------------------------------------------
  Real(dp) :: MassBq(2) = (/ 5.2796_dp, 5.3667_dp /), etaB = 0.55_dp &
    & , FB(2) = (/ 0.193_dp, 0.239_dp /)                             &
    & , FBhatB(2) = (/ 0.216_dp, 0.275_dp /)                         &
    & , TauBq(2) = (/ 1.519e-12_dp, 1.516e-12_dp/) ! in seconds
  Real(dp) :: MassBm(2) = (/ 5.2793_dp, 6.2745_dp /)  &
    & , TauBm(2) = (/ 1.641e-12_dp, 0.452e-12_dp /) ! in seconds
!-----------------------------------------------------------
! constants for K-physics, can be changed using the FLHA
! experimental data from PDG 2013
!-----------------------------------------------------------
  Real(dp) :: MassK0 = 0.497616_dp 
  Real(dp) :: DeltaMK = 3.483e-15_dp 
  Real(dp) :: FK = 0.1558_dp
!-----------------------------------------------------------
! constants for edms
!-----------------------------------------------------------
  Real(dp) :: ecmfactor = 5.975e-15_dp, elimit = 1.6e-27_dp               &
    & , nlimit = 2.9e-26_dp, qcdCorrection = 1.53_dp                      &
    & , CqcdCorrection = 3.4_dp, chiralMass = 1.19_dp, deltaUp = 0.746_dp &
    & , deltaDown = -0.508_dp, deltaStrange = -0.226_dp
! global variables

 Private :: RGE10_SMa

Contains

 Subroutine CalculateRunningMasses(mf_l_in, mf_d_in, mf_u_in, Qlow, alpha &
     &  , alphas, Qhigh, mf_l_out, mf_d_out, mf_u_out, kont)
 !-----------------------------------------------------------------------
 ! calculates running masses except the top in the MSbar scheme
 ! at the scale Qout. The formulas for the decoupling procedure can
 ! be found in K.G.Chetyrkin et al., hep-ph/0004189
 ! input:
 !  mf_l_in ......... lepton onshell masses
 !  mf_d_in ......... d-quark masses, it is assumed the d and s are given
 !                    at the scale Qlow and that mb(mb)
 !  mf_u_in ......... u-quark masses, it is assumed the u and c are given
 !                    at the scale Qlow 
 !  Qlow ............ low energy scale, must be smaller or equal mb(mb)
 !  Qhigh ........... high energy scale, must be larger than Qlow and mb(mb)
 !  alpha ........... alpha_magnetic at Qhigh
 !  alphas .......... alpha_strong at Qhigh
 ! output:
 !  mf_l_out ........ mf_l(Qhigh)
 !  mf_d_out ........ mf_d(Qhigh)
 !  mf_u_out ........ mf_u(Qhigh) [except top quark ]
 !  kont ............ contains the error, =0 if everything is fine
 ! written by Werner Porod, 28.8.99
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  mf_l_in(3), mf_d_in(3), mf_u_in(3), Qlow, alpha &
     &  , alphas, Qhigh
  Real(dp), Intent(out) :: mf_l_out(3), mf_d_out(3), mf_u_out(3)
  Integer, Intent(inout) :: kont

  Real(dp) :: g9(9), g10(10), as_nf, as_nf_minus_1, as_nf_d_pi, aem , tz, dt
  Real(dp), Parameter :: zeta3 = 1.202056903159594285399738161511449990765_dp &
     & , zeta4 = 1.082323233711138191516003696541167902775_dp                 &
     & , B4 = -1.762800087073770864061897634679818807215_dp                   &
     & , c_as(2) = (/ 11._dp / 72._dp                                         &
     &             , 58067._dp / 13824._dp - 82043._dp * zeta3 / 27648._dp /) &
     & , c_m(2) = (/ 89._dp / 432._dp                                         &
     &            , 713._dp/486._dp - B4 / 36._dp - 221._dp * zeta3 / 288._dp &
     &          + 1.25_dp * zeta4  /)
     
  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateRunningMasses"

  kont = 0 
  !-------------------------------------------------------
  ! check if conditions on Qlow and Qhigh are fullfilled
  !-------------------------------------------------------
  If (Qlow.Gt.mf_d_in(3)) Then
   Iname = Iname - 1
   kont = -101
   Call AddError(101)
   Return
  End If
  If ((Qlow.Gt.Qhigh).Or.(Qhigh.Lt.mf_d_in(3))) Then
   Iname = Iname - 1
   kont = - 102
   Call AddError(102)
   Return
  End If

  !-------------------------------------------------------------------------
  ! calculate first coupling at low scales, putting abritary fermion masses
  ! because they are not needed at this stage
  !-------------------------------------------------------------------------
  g10(1) = Sqrt(4._dp * Pi * alphas)
  g10(2) = Sqrt(4._dp * Pi * alpha)
  g10(3:10)= 2._dp

  If (Qhigh.Gt.mf_d_in(3)) Then ! allow that Qhigh==mb(mb)
   tz = Log(mf_d_in(3)/Qhigh) ! at mb(mb) we have to decouple the b-quark
   dt = tz / 50._dp
   Call odeint(g10, 10, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  End If
  g9 = g10(1:9)
  !-------------------------------------------
  ! at mb(mb) we have to decouple the b-quark
  !--------------------------------------------
  If (Qlow.Lt.mf_d_in(3)) Then
   as_nf = g9(1)**2 / (4._dp * Pi)
   as_nf_d_pi = as_nf / Pi
   as_nf_minus_1 = as_nf *(1._dp+ as_nf_d_pi**2 *(c_as(1) +c_as(2)*as_nf_d_pi))
   g9(1) = Sqrt(4._dp * Pi * as_nf_minus_1)
   tz = Log(Qlow/mf_d_in(3)) 
   dt = tz / 50._dp
   Call odeint(g9, 9, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  End If
  !--------------------------------------------------
  ! lepton masses at Qlow, note that aem is alpha/pi
  !--------------------------------------------------
  aem = g9(2)**2 / (4._dp * Pi**2)
  g9(3:5) = mf_l_in * (1._dp - aem)
  !-----------
  ! m_u, m_c
  !-----------
  g9(6:7) = mf_u_in(1:2)
  !-----------
  ! m_d, m_s
  !-----------
  g9(8:9) = mf_d_in(1:2)

  !-------------------------
  ! running back to mb(mb)
  !-------------------------
  If (Qlow.Lt.mf_d_in(3)) Then
   tz = Log(Qlow/mf_d_in(3))
   dt = - tz / 50._dp
   Call odeint(g9, 9, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
   !-------------
   ! thresholds
   !-------------
   AlphaS_mB = as_nf
   g9(1) = Sqrt(4._dp * Pi * as_nf)
   g9(6:9) = g9(6:9) / (1._dp + as_nf_d_Pi**2 * (c_m(1) + c_m(2)*as_nf_d_Pi))  
  End If
  g10(1:9) = g9
  g10(10) = mf_d_in(3)
  !-------------------------
  ! running back to Qhigh
  !-------------------------
  tz = Log(mf_d_in(3)/Qhigh)
  dt = - tz / 50._dp
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)

  mf_l_out = g10(3:5)
  mf_u_out(1:2) = g10(6:7)
  mf_u_out(3) = mf_u_in(3)
  mf_d_out = g10(8:10)

  Iname = Iname - 1

 End Subroutine CalculateRunningMasses

 Subroutine FermionMass(Yuk,vev,mF,U_L,U_R,kont)
 !-----------------------------------------------------------------
 ! calculates fermion masses for a given 3*3 yukawa matrix and a vev
 ! input
 !  Yuk ....... 3*3 Yukawa matrix
 !  vev ....... veve
 ! output
 !  mF(i) ..... fermion masses
 !  U_L(i,j) .. left mixing matrix
 !  U_R(i,j) .. right mixing matrix
 ! written by Werner Porod, 30.07.01
 !----------------------------------------------------------------------
 Implicit None

  Complex(Dp), Intent(in) :: Yuk(3,3)
  Real(Dp), Intent(in) :: vev
  Complex(Dp), Intent(out) :: U_L(3,3), U_R(3,3)
  Real(Dp), Intent(out) :: mF(3)
  Integer, Intent(inout) :: kont

  Complex(Dp) :: mat3(3,3)=0, mat32(3,3)=0, U3(3,3), V3(3,3), phaseM
  Real(dp) :: v1(3,3), u1(3,3), mC2(3), test(2)
  Integer :: ierr, i1

  Iname = Iname + 1
  NameOfUnit(Iname) = "FermionMass"

  mF = 0
  U_L = 0
  U_R = 0

  mat3 = oosqrt2 * vev * Yuk
  mat32 = mat3

  mat32 = Matmul( Transpose( Conjg( mat3 ) ), mat3 )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, v1, ierr, test)
   v3 = v1
  Else
   Call EigenSystem(mat32, mC2, v3, ierr, test)
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( mat3, Transpose( Conjg( mat3 ) ) )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, u1, ierr, test)
   u3 = u1
  Else
   Call EigenSystem(mat32, mC2, u3, ierr, test)
  End If
  u3 = Conjg(u3)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( Matmul( Conjg(u3), mat3), Transpose( Conjg(v3) ) )
  Do i1=1,3
   If ( Abs( mat32(i1,i1) ).Ne.0._dp) Then
    phaseM =   mat32(i1,i1)   / Abs( mat32(i1,i1) )
    v3(i1,:) = phaseM * v3(i1,:)
   End If
  End Do

  Do i1=1,3
   If (Abs( v3(i1,i1) ).Ne.0._dp) Then
    phaseM = v3(i1,i1)   / Abs( v3(i1,i1) )
    v3(i1,:) = Conjg(phaseM) * v3(i1,:)
    u3(i1,:) = phaseM * u3(i1,:)
   End If
  End Do

  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'Warning in subroutine FermionMass, ierr = ',ierr
   Write(ErrCan,*) 'Yuk ',Yuk
   Write(ErrCan,*) 'vev ',vev
   Write(ErrCan,*) ' '
   kont = ierr
   mF(1) = Abs(yuk(1,1)*vev)
   mF(2) = Abs(yuk(2,2)*vev)
   mF(3) = Abs(yuk(3,3)*vev)
   U_L = id3C
   U_R = U_L
   Iname = Iname - 1
   Return
  Endif

  Do i1=1,2
   If ((mC2(i1).Lt.0._dp).And.(Abs(mc2(i1)/mc2(3)).Lt.(10._dp*Epsilon(1._dp)))) &
          &  mc2(i1) = 0._dp
  End Do

  mF = Sqrt(MC2) 
  U_L = u3
  U_R = v3

  Iname = Iname - 1

 End Subroutine FermionMass

 Subroutine QuarkMasses_and_PhaseShifts(Yd, Yu, vSM, mf_d, uD_L, uD_R &
                                     & , mf_u, uU_L, uU_R, CKM)
 !-----------------------------------------------------------------------------
 ! takes the Yukawa couplings, calculates the quark masses and the CKM
 ! from given Yukawa couplings. In this process also phase shifts are
 ! preformed such that the CKM is in the standard form
 ! written by Werner Porod, 27.02.2014  
 !-----------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: Yu(3,3), Yd(3,3)
  Real(dp), Intent(in) :: vSM(2)
  Complex(dp), Dimension(3,3), Intent(out) :: uD_L, uD_R, uU_L, uU_R
  Real(dp), Dimension(3), Intent(out) :: mf_d, mf_u
  Complex(dp), Intent(out), Optional :: CKM(3,3)
  Complex(dp) :: CKM_Q(3,3), EPhi
  Integer :: kont
  Real(dp) :: s13, c13, s23, s12, aR, aI

  Call FermionMass(Yd,vSM(1),mf_d,uD_L,uD_R,kont)
  Call FermionMass(Yu,vSM(2),mf_u,uU_L,uU_R,kont)

  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )
  if (Abs(CKM_Q(1,1)).ne.0._dp) &
        &   uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  if (Abs(CKM_Q(1,2)).ne.0._dp) &
        &   uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  if (Abs(CKM_Q(2,3)).ne.0._dp) &
        &   uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  if (Abs(CKM_Q(3,3)).ne.0._dp) &
        &   uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  
  if (Abs(CKM_Q(1,1)).ne.0._dp) &
        &   uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  if (Abs(CKM_Q(1,2)).ne.0._dp) &
        &   uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  if (Abs(CKM_Q(2,3)).ne.0._dp) &
        &   uU_R(2,:) = uU_R(2,:) / Conjg(CKM_Q(2,3)) * Abs(CKM_Q(2,3))
  if (Abs(CKM_Q(3,3)).ne.0._dp) &
        &   uU_R(3,:) = uU_R(3,:) / Conjg(CKM_Q(3,3)) * Abs(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  !--------------------------------------------------------------
  ! one more freedom left
  !--------------------------------------------------------------
  s13 = Abs(CKM_Q(1,3))
  c13 = Sqrt(1._dp - s13**2)
  s23 = Abs(CKM_Q(2,3))/c13
  s12 = Abs(CKM_Q(1,2))/c13

  aR = Real(CKM_Q(2,2),dp) + s12 * s23 * Real(CKM_Q(1,3),dp)
  aI =  s12 * s23 * Aimag(CKM_Q(1,3)) - Aimag(CKM_Q(2,2))
  Ephi = Cmplx(aR/Sqrt(aR**2+aI**2),aI/Sqrt(aR**2+aI**2),dp)

  uU_L(2:3,:) = Ephi * uU_L(2:3,:)
  uD_L(3,:) = Ephi * uD_L(3,:)
  Ephi = Conjg(Ephi)
  uU_R(2:3,:) = Ephi * uU_R(2:3,:)
  uD_R(3,:) = Ephi * uD_R(3,:)

  If (Present(CKM)) CKM = Matmul(uU_L,Conjg(Transpose(uD_L)))

 End Subroutine QuarkMasses_and_PhaseShifts

 Subroutine InitializeStandardModel
 !-----------------------------------------------------------------
 ! reads in Standard Model Paramters from the file StandardModel.in
 ! In case that the file does not exist, default values are used
 ! defined after label 200
 ! written by Werner Porod, 07.07.02
 !------------------------------------------------------------------
 Implicit None

  Integer          :: i1, kont
  Real(dp) :: oo_alpha, s12, s23, s13, phase, c13, c23, c12 &
                &, TimeMu, TimeTau, g10(10), tz, dt 

  Iname = Iname + 1
  NameOfUnit(Iname) = "InitializeStandardModel"

  mf_nu = 0._dp ! default for neutrino masses

  !------------------------------------------------------------------------
  ! Contributions to alpha(m_Z), based on F. Jegerlehner, hep-ph/0310234
  ! and Fanchiotti, Kniehl, Sirlin PRD 48 (1993) 307
  !------------------------------------------------------------------------
  Delta_Alpha_Lepton = 0.04020_dp
  Delta_Alpha_Hadron = 0.027651_dp

  Open(99,file='StandardModel.in',status='old',err=200)

 !---------
 ! Z-boson  
 !---------
  Read (99,800) mZ      ! mass
  Read (99,800) gamZ    ! width
  Read (99,*) BrZqq     ! branching ratios in q \bar{q} 
  Read (99,*) BrZll     ! branching ratios in leptons
  Read (99,800) BrZinv  ! invisible branching ratio

  mZ2 = mZ**2
  gamZ2 = gamZ**2
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2
 !---------
 ! W-boson
 !---------
  Read (99,800) mW      ! mass
  Read (99,800) gamW    ! width
  Read (99,*) BrWqq     ! branching ratios in q \bar{q}
  Read (99,*) BrWln     ! branching ratios in leptons

  mW2 = mW**2
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2

 !-----------------------------
 ! lepton masses: e, muon, tau
 !-----------------------------
  Read (99,800) mf_l(1)
  Read (99,800) mf_l(2)
  Read (99,800) mf_l(3)

 !---------------------------------------------------------
 ! scale where masses of light quarks are defined [in GeV]
 !---------------------------------------------------------
  Read (99,*) Q_light_quarks
 !--------------------------
 ! up-quark masses: u, c, t
 !--------------------------
  Read (99,800) mf_u(1)
  Read (99,800) mf_u(2)
  Read (99,800) mf_u(3)

 !----------------------------
 ! down-quark masses: d, s, b
 !----------------------------
  Read (99,800) mf_d(1)
  Read (99,800) mf_d(2)
  Read (99,800) mf_d(3)

  Do i1=1,3
   mf_l2(i1) = mf_l(i1)**2
   mf_u2(i1) = mf_u(i1)**2
   mf_d2(i1) = mf_d(i1)**2
  Enddo
 !--------------------------------------------------------------------
 ! couplings: Alpha(Q=0), Alpha(mZ), Alpha_S(mZ), Fermi constant G_F
 !--------------------------------------------------------------------
  Read (99,800) oo_alpha
  Alpha = 1._dp / oo_alpha

  Read (99,800) oo_alpha
  Alpha_mZ = 1._dp / oo_alpha

  Read (99,800) AlphaS_mZ
  Read (99,800) G_F

 !-----------------------------------------
 ! for ISR correction in e+ e- annihilation
 !-----------------------------------------
  KFactorLee = 1._dp + (Pi/3._dp - 1/(2._dp* Pi) ) * Alpha

 !------------
 ! CKM matrix
 !------------
  Read (99,800) s12
  Read (99,800) s23
  Read (99,800) s13
  Read (99,800) phase

  c12 = Sqrt(1._dp-s12*s12)
  c23 = Sqrt(1._dp-s23*s23)
  c13 = Sqrt(1._dp-s13*s13)

  CKM(1,1) = c12 * c13
  CKM(1,2) = s12 * c13
  CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
  CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,3) = s23 * c13
  CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,3) = c23 * c13

  !-----------------------------
  ! width of mu- and tau-lepton
  !-----------------------------
  Read (99,800) TimeMu
  Read (99,800) TimeTau
  
  GammaMu = G_F**2 * mf_l(2)**5 * (1._dp-8._dp * mf_l(1)**2 / mf_l(2)**2 ) &
        & * (1._dp + 0.5_dp * Alpha * (6.25-Pi2)/Pi) / (192._dp*pi*pi2)
  GammaTau = hbar / TimeTau

  Close(99)

  !-------------------------------------------
  ! masses and couplings at m_Z
  !-------------------------------------------
  Call CalculateRunningMasses(mf_l, mf_d, mf_u, Q_light_quarks, alpha_mZ &
     &  , alphas_mZ, mZ, mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)
  !-------------------------------------------
  ! masses and couplings at m_h~125 GeV
  !-------------------------------------------
  tz = Log(mZ/125._dp)
  dt = - tz / 50._dp
  g10(1) = Sqrt(4._dp * Pi * alphas_mZ)
  g10(2) = Sqrt(4._dp * Pi * alpha_mZ)
  g10(3:5) = mf_l_mZ
  g10(6:7) = mf_u_mZ(1:2)
  g10(8:10) = mf_d_mZ
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  AlphaS_125 = oo4pi * g10(1)**2
  Alpha_125 = oo4pi * g10(2)**2
  mf_l_125 = g10(3:5)
  mf_d_125 = g10(8:10)
  mf_u_125(1:2) = g10(6:7)
  mf_u_125(3) = mf_u(3) * (1._dp - alphaS_125/Pi * 4._dp/3._dp)
  !-------------------------------------------
  ! masses and couplings at m_t(m_t)~160 GeV
  !-------------------------------------------
  tz = Log(mZ/160._dp)
  dt = - tz / 50._dp
  g10(1) = Sqrt(4._dp * Pi * alphas_mZ)
  g10(2) = Sqrt(4._dp * Pi * alpha_mZ)
  g10(3:5) = mf_l_mZ
  g10(6:7) = mf_u_mZ(1:2)
  g10(8:10) = mf_d_mZ
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  AlphaS_mt = oo4pi * g10(1)**2
  Alpha_mt = oo4pi * g10(2)**2
  mf_l_mt = g10(3:5)
  mf_d_mt = g10(8:10)
  mf_u_mt(1:2) = g10(6:7)
  mf_u_mt(3) = mf_u(3) * (1._dp - alphaS_mt/Pi * 4._dp/3._dp)

  Iname = Iname - 1
  Return
!  200 Write(*,*) "File StandardModel.in does not exist, using default values."
!  Write(ErrCan,*) "File StandardModel.in does not exist, using default values."

 !---------
 ! Z-boson  
 !---------
200  mZ = 91.1876_dp            ! mass
  gamZ = 2.4952_dp           ! width
  BrZll(1:2) = 0.0336_dp    ! branching ratios in leptons
  BrZll(3) = 0.0338_dp      ! tau+ tau-
  BrZinv = 0.2_dp ! invisible branching ratio
  BrZqq(1) = 0.156_dp       ! d \bar{d}
  BrZqq(2) = 0.156_dp       ! s \bar{s}
  BrZqq(3) = 0.151_dp       ! b \bar{b}
  BrZqq(4) = 0.116_dp       ! u \bar{u}
  BrZqq(5) = 0.12_dp        ! c \bar{c}

  mZ2 = mZ**2
  gamZ2 = gamZ**2
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2
 !---------
 ! W-boson
 !---------
  mW = 80.385_dp      ! mass
  gamW = 2.085_dp     ! width
  BrWqq = 0.35_dp    ! branching ratios in q \bar{q}
  BrWln = 0.1_dp     ! branching ratios in leptons

  mW2 = mW**2
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2

 !-----------------------------
 ! lepton masses: e, muon, tau
 !-----------------------------
  mf_l(1) = 0.51099893e-3_dp
  mf_l(2) = 0.1056583715_dp
  mf_l(3) = 1.77682_dp

 !---------------------------------------------------------
 ! scale where masses of light quarks are defined [in GeV]
 !---------------------------------------------------------
  Q_light_quarks = 2._dp
 !--------------------------
 ! up-quark masses: u, c, t
 !--------------------------
  mf_u(1) = 0.0025_dp 
  mf_u(2) = 1.27_dp
  mf_u(3) = 173.1_dp

 !----------------------------
 ! down-quark masses: d, s, b
 !----------------------------
  mf_d(1) = 0.005_dp
  mf_d(2) = 0.095_dp
  mf_d(3) = 4.18_dp

  Do i1=1,3
   mf_l2(i1) = mf_l(i1)**2
   mf_u2(i1) = mf_u(i1)**2
   mf_d2(i1) = mf_d(i1)**2
  Enddo

 !--------------------------------------------------------------------
 ! couplings: Alpha(Q=0), Alpha(mZ), Alpha_S(mZ), Fermi constant G_F
 !--------------------------------------------------------------------
  Alpha = 1._dp / 137.035999074_dp
  Alpha_mZ = 1._dp / 127.9_dp

  AlphaS_mZ = 0.1184_dp
  G_F = 1.1663787e-5_dp

 !-----------------------------------------
 ! for ISR correction in e+ e- annihilation
 !-----------------------------------------
  KFactorLee = 1._dp + (Pi/3._dp - 1/(2._dp* Pi) ) * Alpha

 !------------
 ! CKM matrix
 !------------
   s12 = lam_wolf
   s23 = s12**2 * A_wolf
   s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2) 
   phase = Atan(eta_wolf/rho_wolf)

  c12 = Sqrt(1._dp-s12*s12)
  c23 = Sqrt(1._dp-s23*s23)
  c13 = Sqrt(1._dp-s13*s13)

  CKM(1,1) = c12 * c13
  CKM(1,2) = s12 * c13
  CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
  CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,3) = s23 * c13
  CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,3) = c23 * c13

  !-----------------------------
  ! width of mu- and tau-lepton
  !-----------------------------
  TimeMu = 2.1969811e-6_dp
  TimeTau = 2.9e-13_dp
  
  GammaMu = hbar / TimeMu
  GammaTau = hbar / TimeTau

  !-------------------------------------------
  ! masses and couplings at m_Z
  !-------------------------------------------
  Call CalculateRunningMasses(mf_l, mf_d, mf_u, Q_light_quarks, alpha_mZ &
     &  , alphas_mZ, mZ, mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)
  !-------------------------------------------
  ! masses and couplings at m_h~125 GeV
  !-------------------------------------------
  tz = Log(mZ/125._dp)
  dt = - tz / 50._dp
  g10(1) = Sqrt(4._dp * Pi * alphas_mZ)
  g10(2) = Sqrt(4._dp * Pi * alpha_mZ)
  g10(3:5) = mf_l_mZ
  g10(6:7) = mf_u_mZ(1:2)
  g10(8:10) = mf_d_mZ
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  AlphaS_125 = oo4pi * g10(1)**2
  Alpha_125 = oo4pi * g10(2)**2
  mf_l_125 = g10(3:5)
  mf_d_125 = g10(8:10)
  mf_u_125(1:2) = g10(6:7)
  mf_u_125(3) = mf_u(3) * (1._dp - alphaS_125/Pi * 4._dp/3._dp)
  !-------------------------------------------
  ! masses and couplings at m_t(m_t)~160 GeV
  !-------------------------------------------
  tz = Log(mZ/160._dp)
  dt = - tz / 50._dp
  g10(1) = Sqrt(4._dp * Pi * alphas_mZ)
  g10(2) = Sqrt(4._dp * Pi * alpha_mZ)
  g10(3:5) = mf_l_mZ
  g10(6:7) = mf_u_mZ(1:2)
  g10(8:10) = mf_d_mZ
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SMa, kont)
  AlphaS_mt = oo4pi * g10(1)**2
  Alpha_mt = oo4pi * g10(2)**2
  mf_l_mt = g10(3:5)
  mf_d_mt = g10(8:10)
  mf_u_mt(1:2) = g10(6:7)
  mf_u_mt(3) = mf_u(3) * (1._dp - alphaS_mt/Pi * 4._dp/3._dp)

  800  Format(f16.7)

  Iname = Iname - 1

 End Subroutine InitializeStandardModel


 Subroutine NeutrinoMasses(MnuL5, mN, N, kont)
 !----------------------------------------------------------------------
 ! calculates neutrino masses + mixing matrix N from the dim 5 operator
 ! input:
 !  MnuL5 .......... dim 5 operator
 ! output 
 !  mN(i) .......... neutrino mass_i
 !  N(i,j) ......... neutrino mixing matrix
 ! written by Werner Porod, 09.01.2006
 ! Note the factor of 1/2 appears because in the RGE running the Higgs field
 ! is still the complex field
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Complex(Dp), Intent(in) :: MnuL5(3,3)
  Real(Dp), Intent(out) :: mN(3)
  Complex(Dp), Intent(out) :: N(3,3)

  Integer :: i1,i2,ierr
  Complex(Dp) :: mat32(3,3), E3(3), phaseM, mat3(3,3)
  Real(Dp) :: N3a(3,3), eig(3), test(2)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutrinoMasses'

  mat3 = 0.5_dp * MnuL5

  If (Maxval(Abs(Aimag(Mat3))).Lt.  &                           ! matrix is reel
     & (Maxval(Abs(Real(Mat3,dp)))*1.e-3_dp * Epsilon(1._dp))) Then

   Call EigenSystem(Real(Mat3,dp), Eig, N3a, ierr, test)

   Do i1=1,3
    If (Eig(i1).Lt.0._dp) Then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N3a(i1,:)
    Else
     mN(i1) = Eig(i1)
     N(i1,:) =N3a(i1,:)
    End If
   End Do

   Do i1=1,2
    Do i2=i1+1,3
     If (mN(i1).Gt.mN(i2)) Then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E3 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E3
     End If
    End Do
   End Do

  Else

   mat32 = Matmul( Transpose(Conjg( Mat3 ) ), Mat3 )
   Call EigenSystem(mat32, Eig, N, ierr, test)
   mat32 = Matmul(Conjg(N), Matmul( Mat3, Transpose( Conjg( N ) ) ) )
   Do i1=1,3
    phaseM =   Sqrt( mat32(i1,i1)   / Abs( mat32(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN = Sqrt( Abs(Eig) ) ! abs to avoid problems with numerical zeros showing up
                         ! as tiny negative numbers
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (Errcan,*) 'Warning in subroutine NeutrinoMasses, ierr =',ierr
   Write (Errcan,*) 'MnuL5:',Cmplx(MnuL5(1,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(2,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(3,:))
   kont = ierr
   Iname = Iname - 1
   Return
  Endif
  !-------------------------------------------
  ! note, that my mixing matrix is the
  !  Transpose( Conjg( Nnu ) )
  ! compared to standard convention in neutrino physics
  !-------------------------------------------

  Iname = Iname - 1

 End Subroutine NeutrinoMasses

 Subroutine RGE10_SMa(len,t,gy,f)
 !--------------------------------------------------------
 ! RGEs within the SM assuming the MSbar scheme
 ! 2-loop RGEs for e
 ! 4-loop RGEs for g_3
 ! 2-loop RGEs for lepton masses
 ! 4-loop QCD and 2-loop QED RGES for quark masses
 ! Assumption: the only threhold to be checked is m_b
 ! input: t = Log(Q^2)
 !        gy(i) ... i=1  -> e(Q)
 !                  i=2  -> g_3
 !                  i=3  -> m_e
 !                  i=4  -> m_mu
 !                  i=5  -> m_tau
 !                  i=6  -> m_u
 !                  i=7  -> m_c
 !                  i=8  -> m_d
 !                  i=9  -> m_s
 !                  i=10 -> m_b, is optional
 ! output:
 !   f = d(gy)/d(t)
 ! written by Werner Porod, 03.12.03
 !--------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: t, gy(len)
  Real(dp), Intent(out) :: f(len)

  Integer :: i1
  Real(dp) :: g32, g34, g36, g38, e2, e4, g32e2, q
  Real(dp), Parameter :: b_e1(2) = (/ 76._dp / 9._dp , 80._dp / 9._dp /)    &
       & , b_e2(2) = (/ 460._dp / 27._dp , 464._dp / 27._dp /)              & 
       & , b_e3(2) = (/ 160._dp / 9._dp , 176._dp / 9._dp  /)               & 
       & , b_g1(2) = (/ -25._dp / 3._dp, -23._dp/3._dp /)                   &
       & , b_g2(2) = (/ -154._dp / 3._dp, -116._dp/3._dp /)                 &
       & , b_g3(2) = (/ 20._dp / 3._dp, 22._dp/3._dp /)                     &
       & , b_g4(2) = (/ -21943._dp/54._dp, 9769._dp/54._dp /)               &
       & , b_g5(2) = (/ -4918247._dp/1458._dp-414140._dp*zeta3/81._dp       &
       &             , 598391._dp/1458._dp - 352864._dp*zeta3/81._dp /)     &
       & , g_el1(2) = (/ -6._dp, -6._dp /)                                  &
       & , g_el2(2) = (/ 353._dp / 9._dp,  373._dp / 9._dp /)               & 
       & , g_eu1(2) = (/ -8._dp/3._dp, -8._dp/3._dp /)                      &
       & , g_eu2(2) = (/ 1472._dp / 81._dp, 1552._dp / 81._dp/)             & 
       & , g_eu3(2) = (/ -32._dp / 9._dp,  -32._dp / 9._dp/)                & 
       & , g_ed1(2) = (/ -2._dp/3._dp, -2._dp/3._dp /)                      &
       & , g_ed2(2) = (/ 377._dp / 81._dp,  397._dp / 81._dp /)             & 
       & , g_ed3(2) = (/ -8._dp / 9._dp,  -8._dp / 9._dp /)                 & 
       & , g_q1(2) = (/ - 8._dp , -8._dp /)                                 &
       & , g_q2(2) = (/ -1052._dp / 9._dp ,  -1012._dp / 9._dp /)           &
       & , g_q3(2) = (/ -144674._dp/81._dp + 1280._dp * zeta3 / 3._dp       &
       &              , -128858._dp/81._dp + 1600._dp * zeta3 / 3._dp /)    &
       & , g_q4(2) = (/ -7330357._dp/243._dp + 51584._dp* zeta3/3._dp       &
       &                - 16000._dp*zeta4 / 3._dp + 11200._dp* zeta5 /9._dp &
       &             , -1911065._dp/81._dp + 618400._dp* zeta3/27._dp       &
       &                - 18400._dp*zeta4 / 3._dp - 25600._dp* zeta5 /9._dp  /)

       
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RGE10_SMa'

  q = t

  If (len.Eq.9) Then ! check which beta function (anomalous dimension) to use
   i1 = 1
  Else If (len.Eq.10) Then
   i1 = 2
  Else
   Write(ErrCan,*) "Error in routine "//Trim(NameOfUnit(Iname))
   Write(ErrCan,*) "Length of the vector gy = ",len
   Call TerminateProgram
  End If

  g32 = gy(1)**2
  g34 = gy(1)**4
  g36 = gy(1)**6
  g38 = gy(1)**8
  e2 = gy(2)**2
  e4 = gy(2)**4
  g32e2 = g32 * e2 
 !--------
 ! g_3
 !--------
  f(1) = oo16pi2 * gy(1) * ( b_g1(i1)*g32                                     &
       &                   + oo16pi2 * ( b_g2(i1)*g34 + b_g3(i1)*g32e2        &
       &                               + oo16pi2 * ( b_g4(i1)*g36             &
       &                                           + oo16pi2 * b_g5(i1)*g38 )))
 !--------
 ! e
 !--------
  f(2) = oo16pi2 * gy(2) * ( b_e1(i1) * e2                                &
       &                   + oo16pi2 * (b_e2(i1) * e4 + b_e3(i1) * g32e2 ))
 !-----------------
 ! m_l, l=e,mu,tau
 !-----------------
  f(3:5) =  oo16pi2 * gy(3:5) * (g_el1(i1) * e2 + oo16pi2 *g_el2(i1) * e4)
 !---------
 ! m_u, m_c
 !---------
  f(6:7) = oo16pi2 * gy(6:7) * (g_eu1(i1) * e2 + g_q1(i1) * g32              &
         &                     + oo16pi2 * (g_eu2(i1)*e4 + g_eu3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))
 !---------------
 ! m_d, m_s, m_b
 !---------------
  f(8:len) = oo16pi2 * gy(8:len) * (g_ed1(i1) * e2 + g_q1(i1) * g32          &
         &                     + oo16pi2 * (g_ed2(i1)*e4 + g_ed3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))

  Iname = Iname - 1

 End Subroutine RGE10_SMa


End Module StandardModel
