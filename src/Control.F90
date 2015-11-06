Module Control
!-------------------------------------------------------------------
! This module contains the variables and subroutines that control
! input, output and error handling. They need to be in one module,
! because some fatal errors require a closing of all open files to save
! as much information as possible.
! ErrorLevel:
!
!  -2 ... do not print severe warnings
!
!  -1 ... do not print warnings
!
!   0 ... print every warning
!
!   1 ... abort in case of a severe warning
!
!   2 ... abort even in case of a warning
! 13.09.2002: including the definition of dp for the determination
!             of the precision used
! 05.01.04: putting all numerical constants here as this file is used
!           in every other module
! 01.03.04: start to build a collection system on error messages,
!           -creating new array i_errors containing the information
!            how often a warning and/or error has been produced
!           -new routines to add 1 in array errors, to clean errors
!            and to write the corresponding numbers 
!-------------------------------------------------------------------

 Interface Is_NaN
  Module Procedure IsNaN_r, IsNaN_v, IsNaN_m
 End Interface

 Interface FindPosition
  Module Procedure FindPosition_R
 End Interface

! global variables
! double and quadrupole precision
#ifdef QUADRUPOLE
Integer, Parameter :: dp = Selected_real_kind(25,450)
#else
Integer, Parameter :: dp = Selected_real_kind(15,300)
#endif

#ifdef ONLYDOUBLE
Integer, Parameter :: qp = Selected_real_kind(15,300)
#else
Integer, Parameter :: qp = Selected_real_kind(25,450)
#endif

 !------------------------------------
 ! version number
 !------------------------------------
 Character(len=8), Save :: version="v3.3.8"
 !------------------------------------
 ! variables for spectrum calculation
 !------------------------------------
 Real(dp) :: delta_mass = 1.e-5_dp    ! precision of mass calculation
 Logical :: WriteOut = .False.        ! write out debug information
 Integer :: n_run = 40         ! maximal number of iterations in mass calculation
 Logical :: GenerationMixing=.False.  ! if (s)fermion generations mix
 Logical :: FermionMassResummation=.True.  ! if contributions to fermion
                                           ! mass matrix are resummed
 !--------------------------------------
 ! variables for process calculations
 !--------------------------------------
 Logical :: L_CS = .False. ! if cross sections should be calculated at all
 Logical :: L_CSrp = .False. ! if cross sections should be calculated in case of RP violation
 Logical :: L_BR=.True.    ! if branching ratios should be calculated at all
 !---------------------------------------------------
 ! for debugging, to get information on the screen
 !---------------------------------------------------
 Logical :: output_screen = .False.
 !--------------------------------------------------------------
 ! if Higgs masses and mixings from external programs are used
 !--------------------------------------------------------------
 Logical :: External_Higgs = .False.
 !--------------------------------------------------------------
 ! if mass spectrum and mixings from external programs are used
 !--------------------------------------------------------------
 Logical :: External_Spectrum = .False.
 !------------------------------------------
 ! variables for R-parity violation
 !------------------------------------------
 Logical :: l_fit_RP_parameters = .False. ! if true fit bilinear R-parity parameters
                                          ! such that neutrino data are fullfilled
 !-------------------------------------------------------------
 ! variables for super PMNS basis if only mixing is given
 !-------------------------------------------------------------
 Logical, save :: fake_m_nu = .True.
 !-------------------------------------------------------------------
 ! in case one still wants to check results using the old BoundaryEW
 ! use entry 9 in SPhenoInput to set value to 1 to use old version
 !-------------------------------------------------------------------
 Logical, save :: UseNewBoundaryEW = .True.
 !-------------------------------------------------------------------
 ! in case one still wants to check results using the on-shell stop
 ! masses to get the renormalisation scale, 
 ! use entry 10 in SPhenoInput to set value to 1 
 !-------------------------------------------------------------------
 Logical, save :: UseNewScale = .True.
 !------------------------------------------
 ! warning and error system
 !------------------------------------------
 Integer, Save :: Iname, ErrorLevel=0, NumberOfOpenFiles=1, ErrCan=10
 Character (Len=40), Save :: NameOfUnit(40)
 Logical, Save :: Non_Zero_Exit = .False.
 Logical, Save :: Write_warning_to_screen = .False.
 !-------------------------
 ! for the error system 
 !-------------------------
 Character(len=60) :: Math_Error(31) =                                &
  & (/ "Routine OdeInt, stepsize smaller than minimum:              " &
  &  , "Routine OdeInt, maximal value > 10^36:                      " &
  &  , "Routine OdeInt, too many steps:                             " &
  &  , "Routine OdeIntB, boundary condition not fullfilled:         " &
  &  , "Routine OdeIntB, stepsize smaller than minimum:             " &
  &  , "Routine OdeIntB, maximal value > 10^36:                     " &
  &  , "Routine OdeIntB, too many steps:                            " &
  &  , "Routine OdeIntC, boundary condition not fullfilled:         " &
  &  , "Routine OdeIntC, stepsize smaller than minimum:             " &
  &  , "Routine OdeIntC, maximal value > 10^36:                     " &
  &  , "Routine OdeIntC, too many steps:                            " &
  &  , "Routine rkqs, stepsize underflow:                           " &
  &  , "Routine ComplexEigenSystem, array dimensions do not match:  " &
  &  , "Routine ComplexEigenSystem, numerical precision not reached:" &
  &  , "Routine RealEigenSystem, array dimensions do not match:     " &
  &  , "Routine RealEigenSystem, numerical precision not reached:   " &
  &  , "Routine tqli, array dimensions do not match:                " &
  &  , "Routine tqli, too many iterations:                          " &
  &  , "Routine Dgauss, too high accuracy required:                 " &
  &  , "Routine DgaussInt, too high accuracy required:              " &
  &  , "Routine Kappa, precision problem:                           " &
  &  , "Routine IntRomb: step size to small:                        " &
  &  , "Routine IntRomb: too many steps:                            " &
  &  , "Routine GaussJ: singular matrix                             " &
  &  , "Routine InverseMatrix: inversion failed                     " &
  &  , "Routine InvMat3: inversion failed, det(M)=0                 " &
  &  , "Routine bsstep, stepsize underflow:                         " &
  &  , "Routine pzextr: probable misuse, too much extrapolation     " &
  &  , "Routine rzextr: probable misuse, too much extrapolation     " &
  &  , "Routine RealEigenSystem, matrix contains NaN                " &
  &  , "Routine ComplexEigenSystem, matrix contains NaN             " &
  & /)
 Character(len=60) :: MathQP_Error(10) =                              &
  & (/ "Routine ComplexEigenSystemDP, array dimensions do not match:" &
  &  , "Routine ComplexEigenSystemDP, numer. precision not reached: " &
  &  , "Routine ComplexEigenSystemQP, array dimensions do not match:" &
  &  , "Routine ComplexEigenSystemQP, numer. precision not reached: " &
  &  , "Routine RealEigenSystem_DP, array dimensions do not match:  " &
  &  , "Routine RealEigenSystem_DP, numerical precision not reached:" &
  &  , "Routine RealEigenSystem_QP, array dimensions do not match:  " &
  &  , "Routine Tqli_QP, array dimensions do not match:             " &
  &  , "Routine Tqli_QP, too many iterations:                       " &
  &  , "Routine Tql2_QP, too many iterations:                       " &
  & /)
 Character(len=60) :: SM_Error(2) =                                   &
  & (/ "Routine CalculateRunningMasses: Q_low > m_b(m_b)            " &
  & ,  "Routine CalculateRunningMasses: Max(Q_low,m_b(m_b) > Q_max  " &
  & /)
 Character(len=60) :: SusyM_Error(33) =                                  &
  & (/ "Routine ChargedScalarMassEps1nt: negative mass squared      "    &
  & ,  "Routine ChargedScalarMassEps3nt: negative mass squared      "    &
  & ,  "Routine ChargedScalarMassLam3nt: negative mass squared      "    &
  & ,  "Routine CharginoMass3: |y_tau|^2 < 0                        "    &
  & ,  "Routine CharginoMass5: |y_tau|^2 < 0                        "    &
  & ,  "Routine PseudoScalarMassEps1nt: negative mass squared       "    &
  & ,  "Routine PseudoScalarMassEps3nt: negative mass squared       "    &
  & ,  "Routine PseudoScalarMassMSSMnt: negative mass squared       "    &
  & ,  "Routine PseudoScalarMassSpon1Gen: negative mass squared     "    &
  & ,  "Routine ScalarMassEps1nt: negative mass squared             "    &
  & ,  "Routine ScalarMassEps3nt: negative mass squared             "    &
  & ,  "Routine ScalarMassMSSMeff: negative mass squared            "    &
  & ,  "Routine ScalarMassMSSMnt: negative mass squared             "    &
  & ,  "Routine ScalarMassNMSSMeff: L*k*tanbq*mu = 0                "    &
  & ,  "Routine ScalarMassNMSSMeff: mS0^2_1 < 0                     "    &
  & ,  "Routine ScalarMassNMSSMeff: mP0^2_1 < 0                     "    &
  & ,  "Routine ScalarMassNMSSMeff: mH+^2 < 0                       "    &
  & ,  "Routine ScalarMassSpon1Gen: negative mass squared           "    &
  & ,  "Routine SdownMass3Lam: negative mass squared                "    &
  & ,  "Routine SfermionMass1Eps1: negative mass squared            "    &
  & ,  "Routine SfermionMass1Eps3: negative mass squared            "    &
  & ,  "Routine SfermionMass1MSSM: negative mass squared            "    &
  & ,  "Routine SfermionMass3MSSM: negative mass squared            "    &
  & ,  "Routine SquarkMass3Eps: negative mass squared               "    &
  & ,  "Routine TreeMassesEps1: negative sneutrino mass squared     "    &
  & ,  "Routine TreeMassesMSSM: negative sneutrino mass squared     "    &
  & ,  "Routine TreeMassesMSSM: negative pseudoscalar mass squared  "    &
  & ,  "Routine TreeMassesMSSM: negative charged higgs mass squared "    &
  & ,  "Routine TreeMassesMSSM2: negative sneutrino mass squared    "    &
  & ,  "Routine TreeMassesMSSM2: negative pseudoscalar mass squared "    &
  & ,  "Routine TreeMassesMSSM2: negative charged higgs mass squared"    &
  & ,  "Routine TreeMassesMSSM3: negative sneutrino mass squared    "    &
  & ,  "Routine TreeMassesNMSSM: negative sneutrino mass squared    "    &
  & /)
 Character(len=60) :: InOut_Error(15) =                                  &
  & (/ "Routine LesHouches_Input: unknown error occured:            "    &
  & ,  "Routine LesHouches_Input: Unknown entry for Block MODSEL    "    &
  & ,  "LesHouches_Input: model must be specified before parameters "    &
  & ,  "Routine LesHouches_Input: Unknown entry for Block MINPAR    "    &
  & ,  "LesHouches_Input: model has not been specified completly    "    &
  & ,  "LesHouches_Input: a serious error has been part of the input"    &
  & ,  "LesHouches_Input: Higgs sector has not been fully specified "    &
  & ,  "ReadMatrixC: indices exceed the given boundaries            "    &
  & ,  "ReadMatrixR: indices exceed the given boundaries            "    &
  & ,  "ReadVectorC: index exceeds the given boundaries             "    &
  & ,  "ReadVectorR: index exceeds the given boundaries             "    &
  & ,  "ReadTensorC: indices exceed the given boundaries            "    &
  & ,  "LesHouches_Input, GMSB: Lambda < M_M                        "    &
  & ,  "LesHouches_Input: non-perturbative Yukawa coupling as input "    &
  & ,  "LesHouches_Input: non-diagonl entries in a diagonal matrix  "    &
  & /)

 Character(len=60) :: Sugra_Error(19) =                                  &
  & (/ "Routine BoundaryEW: negative scalar mass as input           "    &
  & ,  "Routine BoundaryEW: mZ^2(mZ) < 0                            "    &
  & ,  "Routine BoundaryEW: sin^2(theta_DR) < 0                     "    &
  & ,  "Routine BoundaryEW: mW^2 < 0                                "    &
  & ,  "Routine BoundaryEW: m_l_DR/m_l < 0.1 or m_l_DR/m_l > 10     "    &
  & ,  "Routine BoundaryEW: m_b_DR/m_b < 0.1 or m_b_DR/m_b > 10     "    &
  & ,  "Routine BoundaryEW: m_t_DR/m_t < 0.1 or m_t_DR/m_t > 10     "    &
  & ,  "Routine RunRGE: entering non-perturbative regime            "    &
  & ,  "Routine RunRGE: g1 and g2 do not meet at the high scale     "    &
  & ,  "Routine RunRGE: entering non-perturbative regime at M_GUT   "    &
  & ,  "Routine RunRGE: entering non-perturbative regime at M_H3    "    &
  & ,  "Routine Sugra: run did not converge                         "    &
  & ,  "Routine Calculate_Gi_Yi: mZ^2(mZ) < 0                       "    &
  & ,  "Routine Calculate_Gi_Yi: too many iterations of mb(mb)      "    &
  & ,  "Routine Sugra:  |mu|^2 < 0 at m_Z                           "    &
  & ,  "Routine RunRGE_2: entering non-perturbative regime          "    &
  & ,  "Routine RunRGE_2: g1 and g2 do not meet at the high scale   "    &
  & ,  "Routine RunRGE_2: entering non-perturbative regime at M_GUT "    &
  & ,  "Routine RunRGE_2: entering non-perturbative regime at M_H3  "    &
  & /)
 Character(len=60) :: LoopMass_Error(25) =                               &
  & (/ "SleptonMass_1L: encountered a negative mass squared         "    &
  & ,  "SleptonMass_1L: p^2 iteration did not converge              "    &
  & ,  "SneutrinoMass_1L: encountered a negative mass squared       "    &
  & ,  "SneutrinoMass_1L: p^2 iteration did not converge            "    &
  & ,  "SquarkMass_1L: encountered a negative mass squared          "    &
  & ,  "SquarkMass_1L: p^2 iteration did not converge               "    &
  & ,  "LoopMassesMSSM: m_h0^2 <0 in the 1-loop calculation         "    &
  & ,  "LoopMassesMSSM: m_H+^2 <0 in the 1-loop calculation         "    &
  & ,  "LoopMassesMSSM: m_A0^2 <0 in the 1-loop calculation         "    &
  & ,  "LoopMassesMSSM: |mu|^2 > 10^20 in the 1-loop calculation    "    &
  & ,  "LoopMassesMSSM: |mu|^2 < 0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM: mZ^2(mZ) < 0                                "    &
  & ,  "LoopMassesMSSM2: m_h0^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM2: m_H+^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM2: m_A0^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM2: |mu|^2 > 10^20 in the 1-loop calculation   "    &
  & ,  "LoopMassesMSSM2: |mu|^2 < 0 in the 1-loop calculation       "    &
  & ,  "LoopMassesMSSM2: mZ^2(mZ) < 0                               "    &
  & ,  "LoopMassesMSSM3: m_h0^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM3: m_H+^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM3: m_A0^2 <0 in the 1-loop calculation        "    &
  & ,  "LoopMassesMSSM3: |mu|^2 > 10^20 in the 1-loop calculation   "    &
  & ,  "LoopMassesMSSM3: |mu|^2 < 0 in the 1-loop calculation       "    &
  & ,  "LoopMassesMSSM3: mZ^2(mZ) < 0                               "    &
  & ,  "Sigma_SM_chirally_enhanced: negative mass squared           "    &
  & /)
 Character(len=60) :: TwoLoopHiggs_Error(9) =                            &
  & (/ "PiPseudoScalar2: encountered negative stop mass squared     "    &
  & ,  "PiPseudoScalar2: encountered negative sbottom mass squared  "    &
  & ,  "PiPseudoScalar2: encountered negative stau mass squared     "    &
  & ,  "PiScalar2: encountered negative stop mass squared           "    &
  & ,  "PiScalar2: encountered negative sbottom mass squared        "    &
  & ,  "PiScalar2: encountered negative stau mass squared           "    &
  & ,  "Two_Loop_Tadpoles: encountered negative stop mass squared   "    &
  & ,  "Two_Loop_Tadpoles: encountered negative sbottom mass squared"    &
  & ,  "Two_Loop_Tadpoles: encountered negative stau mass squared   "    &
  & /)

 Integer, Private :: i_errors(1100)=0

!----------------------------
! numerical constants
!----------------------------
! powers and multiples of Pi
!----------------------------
Real(dp), Parameter :: Pi = 3.1415926535897932384626433832795029_dp           &
   & , Pi2 = Pi**2, Pi3 = Pi2*Pi                                              &
   & , pi2o3 = Pi2 / 3._dp, pi2o6 =0.5_dp * Pi2o3, pi2o12 = 0.5_dp * Pi2o6    &
   & , oo3pi = 1._dp / (3._dp * pi), oo4Pi = 1._dp / (4._dp * Pi)             &
   & , oo8Pi = 0.5_dp * oo4Pi, oo16Pi = 0.5_dp * oo8Pi                        &
   & , oo32Pi = 0.5_dp * oo16Pi, oo64Pi = 0.5_dp * oo32Pi                     &
   & , oo48Pi = oo4Pi / 12._dp, fo3Pi = 4._dp / (3._dp * Pi)                  &
   & , oo4pi2 = 1._dp / (4._dp * Pi2), oo3pi2 = 1._dp / (3._dp * Pi2)         &
   & , oo6pi2 = 1._dp / (6._dp * Pi2), oo8pi2 = 1._dp / (8._dp * Pi2)         &
   & , oo16pi2 = 0.5_dp * oo8pi2, oo64pi2 = 0.25_dp * oo16pi2                 &
   & , oo32pi2 = 0.5_dp * oo16pi2, oo48pi2 = oo16pi2 / 3._dp                  &
   & , oo36pi3 = 1._dp / (36._dp * Pi3)                                       &
   & , oo128pi3 = 1._dp / (128._dp * Pi3)                                     &
   & , oo256pi3 = 1._dp / (256._dp * Pi3), oo512pi3 = 0.5_dp * oo256pi3    
!----------------------------------------------
! special values of Logs
!----------------------------------------------
Real(dp), Parameter :: Log2 = 0.6931471805599453094172321214581765680755_dp
!----------------------------------------------
! special values of the Rieman zeta function
!----------------------------------------------
Real(dp), Parameter :: Zeta2 = Pi2 / 6._dp                    &   
   & , zeta3 = 1.202056903159594285399738161511449990765_dp   &
   & , zeta4 = 1.082323233711138191516003696541167902775_dp   &
   & , zeta5 = 1.036927755143369926331365486457034168057_dp
!----------------------------------------------
! special values of the Gamma function
!----------------------------------------------
 Real(dp), Parameter :: &
   &   Gamma_1o3 = 2.678938534707747633655692940974677644129_dp         &
   & , Gamma_2o3 = 1.354117939426400416945288028154513785519_dp
!---------------------------------------------
! roots and powers of roots
!---------------------------------------------
Real(dp), Parameter :: &
    &   Sqrt2 = 1.41421356237309504880168872420969807856967187537694807317_dp &
    & , ooSqrt2 = 1._dp / Sqrt2                                               &
    & , Sqrt3 = 1.73205080756887729352744634150587236694280525381038062805_dp &
    & , ooSqrt3 = 1._dp / Sqrt3
!------------------------------
! special numbers
!------------------------------
 Complex(dp), Parameter :: ZeroC=(0._dp,0._dp), Ic = (0._dp, 1._dp)
 Real(dp), Parameter :: NearlyZero = 10._dp * Epsilon(1._dp)
!------------------------------
! matrices 
!------------------------------
 Real(dp), Parameter :: id2R(2,2) = &
    & Reshape( Source = (/ 1, 0, 0, 1 /), shape = (/2, 2/) )
 Complex(dp), Parameter :: id2C(2,2) = &
    & Reshape( Source = (/ 1, 0, 0, 1 /), shape = (/2, 2/) )
 Complex(dp), Parameter :: id3C(3,3) = &
    & Reshape( Source = (/ 1, 0, 0, & 
    &                      0, 1, 0, &
    &                      0, 0, 1 /), shape = (/3, 3/) )
 Complex(dp), Parameter :: id4C(4,4) = &
    & Reshape( Source = (/ 1, 0, 0, 0, &
    &                      0, 1, 0, 0, &
    &                      0, 0, 1, 0, &
    &                      0, 0, 0, 1 /), shape = (/4, 4/) )
 Complex(dp), Parameter :: id6C(6,6) = &
    & Reshape( Source = (/ 1, 0, 0, 0, 0, 0, &
    &                      0, 1, 0, 0, 0, 0, &
    &                      0, 0, 1, 0, 0, 0, &
    &                      0, 0, 0, 1, 0, 0, &
    &                      0, 0, 0, 0, 1, 0, &
    &                      0, 0, 0, 0, 0, 1 /), shape = (/6, 6/) )
 Complex(dp), Parameter :: Zero33C(3,3) = ZeroC

 !--------------------------------
 ! needed for the type definition
 !--------------------------------
 Integer, Parameter, Private :: n2=200, n3=600

 Type particle2
  Real(dp) :: m, m2 ! mass
  Real(dp) :: g  ! total width
  Integer  :: id ! internal particle identity code
  Integer  :: id2(n2,2) ! particle codes for the 2-body final states
  Real(dp) :: gi2(n2) ! partial widths for 2-body final states
  Real(dp) :: bi2(n2) ! branching ratios for 2--body decays
 End Type particle2 ! still missing are names

 Type particle23
  Real(dp) :: m, m2  ! mass + mass squared
  Real(dp) :: g  ! total width
  Integer  :: id ! internal particle identity code
  Integer  :: id2(n2,2) ! particle codes for the 2-body final states
  Integer  :: id3(n3,3) ! particle codes for the 2-body final states
  Real(dp) :: gi2(n2), gi3(n3) ! partial widths for 2- and 3-body final states
  Real(dp) :: bi2(n2), bi3(n3) ! branching ratios for 2- and 3-body decays
 End Type particle23 ! still missing are names

! global variables

Contains 


 Subroutine AddError(i)
 !--------------------------------------------------------------------
 ! is used to count number of errors
 ! written by Werner Porod, 01.03.04
 !--------------------------------------------------------------------
 Implicit None
  Integer, intent(in) :: i
  i_Errors(i) = i_Errors(i) + 1
 End Subroutine AddError


 Character Function Bu(n)
 !---------------------------------------------------------------
 ! changes integer n [0,9] into the corresponding character
 !---------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Bu = Char(48+n)
 End Function Bu


 Subroutine CheckIndex(i,range,name1,name2)
 !--------------------------------------------------------------------
 ! checks if the index i is between 1 and range 
 ! in case that it is out range, the following actions happen:
 !    ErrorLevel >= 0, termination of program after error message is
 !                     written
 !    ErroLevel = -1, only warning, this can lead to side effects  
 ! written by Werner Porod, 20.9.2000
 !--------------------------------------------------------------------
 Implicit None

  Integer, Intent(In)            :: i,range
  Character (Len=5), Intent(In)  :: name1,name2

  If ((i.Lt.1).Or.(i.Gt.range)) Then
   If (ErrorLevel.Ge.0) Then
    Write (ErrCan,*) 'Severe problem in unit ',NameOfUnit(Iname)
    Write (ErrCan,*) 'Index ',name1,'is out of range: (',name1,name2, &
   &             ') = ',i,range
    Call TerminateProgram
   Else
    Write (ErrCan,*) 'Severe warning from unit ',NameOfUnit(Iname)
    Write (ErrCan,*) 'Index ',name1,'is out of range: (',name1,name2, &
   &             ') = ',i,range
   Endif
  Endif

 End Subroutine CheckIndex


 Subroutine CleanErrors()
 !--------------------------------------------------------------------
 ! is used to reset the array i_errors
 ! written by Werner Porod, 01.03.04
 !--------------------------------------------------------------------
 Implicit None
  i_Errors = 0
 End Subroutine CleanErrors


 Subroutine Closing()
 Implicit None
 Integer :: i1

 Do i1=ErrCan,ErrCan+NumberOfOpenFiles-1
  Close(i1)
 End Do

 End Subroutine Closing


 Subroutine GetError(i)
 !--------------------------------------------------------------------
 ! is used to count number of errors
 ! written by Werner Porod, 01.03.04
 !--------------------------------------------------------------------
 Implicit None
  Integer, intent(out) :: i(size(i_errors))
  i = i_Errors
 End Subroutine GetError


 Subroutine FindPosition_R(n, mat, wert, ii, jj)
 !------------------------------------------------------------------------
 ! finds the position of wert in an n times n matrix
 ! written by Werner Porod, 21.05.2007
 !------------------------------------------------------------------------
 Implicit None
  !-------
  ! input
  !-------
  Integer, Intent(in) :: n         ! dimension of matrix
  Real(dp), Intent(in) :: mat(n,n) ! matrix to be checked
  Real(dp), Intent(in) :: wert     ! value to be looked for
  !--------
  ! output
  !--------
  Integer, Intent(out) :: ii, jj   ! position mat(ii,jj) corresponding to wert
 
  Integer :: i1, i2
 
  Do i1=1,n
   Do i2=1,n
    If (wert.Eq.mat(i1,i2)) Then
     ii = i1
     jj = i2
     Return
    End If
   End Do
  End Do

 End Subroutine FindPosition_R


 Subroutine InitializeControl(io, outfile, name)
 !---------------------------------------------------------------
 ! input:
 !  io ........ channel, to which normal output is written
 !  outfile ... file name connected to channel io
 !  name ...... name of the main program
 !---------------------------------------------------------------
 Implicit None
  Character(len=*), intent(in) :: outfile
  Character(len=*), intent(in) :: name
  Integer, intent(in) :: io
  Iname = 1
  NameOfUnit(Iname) = name

  Open(ErrCan,file="Messages.out",status="unknown")
  Open(io,file=outfile,status="unknown")
  NumberOfOpenFiles = NumberOfOpenFiles + 1


  ErrorLevel = 0
  GenerationMixing=.False.
  L_BR=.True.
  L_CS=.True.

  Open(99,file="Control.in",status="old",err=200)
   Read(99,*) ErrorLevel
   Read(99,*) GenerationMixing
   Read(99,*) L_BR
   Read(99,*) L_CS
  Close(99)

 200  Return

 End Subroutine InitializeControl


 Logical Function IsNaN_r(x)
 !---------------------------------------------------------------------
 ! tests if x is NaN. Comparison of NaN with any number gives FALSE
 !---------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  IsNaN_r = .Not. ((x.Gt.0._dp).Or.(x.Lt.0._dp).Or.(x.Eq.0._dp))

 End Function IsNaN_r


 Logical Function IsNaN_v(x)
 !---------------------------------------------------------------------
 ! tests if at least one element of x is NaN. 
 ! Comparison of NaN with any number gives FALSE
 !---------------------------------------------------------------------
 Real(dp), Intent(in) :: x(:)

 Integer :: l1, i1

 l1 = Size(x)

 IsNaN_v = .False.
 
 Do i1=1,l1
  IsNaN_v = .Not. ((x(i1).Gt.0._dp).Or.(x(i1).Lt.0._dp).Or.(x(i1).Eq.0._dp))
  If (IsNaN_v) Exit
 End Do

 End Function IsNaN_v

 Logical Function IsNaN_m(x)
 !---------------------------------------------------------------------
 ! tests if at least one element of x is NaN. 
 ! Comparison of NaN with any number gives FALSE
 !---------------------------------------------------------------------
 Real(dp), Intent(in) :: x(:,:)

 Integer :: l1, i1, i2, l2

 l1 = Size(x,dim=1)
 l2 = Size(x,dim=2)

 IsNaN_m = .False.
 
 Do i1=1,l1
  Do i2=1,l2
   IsNaN_m = .Not. ((x(i1,i2).Gt.0._dp).Or.(x(i1,i2).Lt.0._dp).Or.(x(i1,i2).Eq.0._dp))
   If (IsNaN_m) Exit
  End Do
 End Do

 End Function IsNaN_m

 Subroutine ModelNotIncluded(i1,i2,i3,i4)
 !-------------------------------------------------------------------
 ! States which models are not yet included in a specific subroutine
 ! written by Werner Porod, 20.9.2000
 !-------------------------------------------------------------------
 Implicit None
  
  Integer, Intent(In) :: i1,i2,i3,i4

  If (ErrorLevel.Le.0) Then
   Write (ErrCan,*) 'Warning from unit ',NameOfUnit(Iname)
   If (i1.Ne.0) Write(ErrCan,*) ' MSSM not yet included '
   If (i2.Ne.0) Write(ErrCan,*) ' 1-generation epsilon model not yet included '
   If (i3.Ne.0) Write(ErrCan,*) ' 3-generation epsilon model not yet included '
   If (i4.Ne.0) Write(ErrCan,*) ' full RP model not yet included '
  End If

 End Subroutine ModelNotIncluded


 Subroutine TerminateProgram
 !-----------------------------------------------------------------------
 ! This subroutine terminates a program if a fatal error occurs.
 ! Before doing this, it writes the tree of calling subroutines to
 ! the file which is connected to the channel ErrCan
 ! written by Werner Porod, 20.9.2000
 !-----------------------------------------------------------------------
 Implicit None

  Integer :: i1

  Write (ErrCan,*) "  "
  Write (ErrCan,*) "ErrorLevel, Iname:",ErrorLevel, Iname
  Write (ErrCan,*) &
    & "The error has occured in the following chain of subroutines:"
  Do i1=1,Size(NameOfUnit)
   Write (ErrCan,*) NameOfUnit(i1)
  End Do
  Write (ErrCan,*) "  "
  Write (ErrCan,*) "Hopefully you find the error soon"
  Do i1=ErrCan,ErrCan+NumberOfOpenFiles-1
   Close(i1)
  End Do
  Stop "Subroutine TerminateProgram"

 End Subroutine TerminateProgram


 Subroutine WriteNumberOfErrors(io)
 Implicit None
  Integer, Intent(in) :: io

  Integer :: i1
! module Mathematics
  Do i1=1,100
   If (i_errors(i1).Gt.0) Write(io,*) Trim(Math_Error(i1)),i_errors(i1)
  End Do
! module MathematicsQP
  Do i1=1,100
   If (i_errors(1000+i1).Gt.0) Write(io,*) Trim(MathQP_Error(i1)),i_errors(i1)
  End Do
! module StandardModel
  Do i1=1,100
   If (i_errors(100+i1).Gt.0) Write(io,*) Trim(SM_Error(i1)),i_errors(100+i1)
  End Do

! module SusyMasses
  Do i1=1,100
   If (i_errors(200+i1).Gt.0) Write(io,*) &
      &    Trim(SusyM_Error(i1)),i_errors(200+i1)
  End Do

! module InputOutput
  Do i1=1,100
   If (i_errors(300+i1).Gt.0) Write(io,*) &
      &    Trim(InOut_Error(i1)),i_errors(300+i1)
  End Do

! module Sugra
  Do i1=1,100
   If (i_errors(400+i1).Gt.0) Write(io,*) &
      &    Trim(Sugra_Error(i1)),i_errors(400+i1)
  End Do

! module LoopMasses
  Do i1=1,100
   If (i_errors(500+i1).Gt.0) Write(io,*) &
      &    Trim(LoopMass_Error(i1)),i_errors(500+i1)
  End Do

! module TwoLoopHiggs
  Do i1=1,100
   If (i_errors(600+i1).Gt.0) Write(io,*) &
      &    Trim(LoopMass_Error(i1)),i_errors(550+i1)
  End Do

 End Subroutine WriteNumberOfErrors


End Module Control

