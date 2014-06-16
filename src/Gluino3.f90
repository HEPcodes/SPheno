Module Gluino3Decays

! load modules
Use Control
Use DecayFunctions, Only : FermionToFermionScalar
Use StandardModel, Only: mW, mW2, mf_d , mf_d2, mf_u, mf_u2, mZ2, CKM
Use ThreeBodyPhaseSpace
Use ThreeBodyPhaseSpaceS
Use LoopFunctions, Only: Igamma, I2gamma, Jgamma, Kgamma &
    & , GetRenormalizationScale
! load modules

! private variables
 ! variables for gluino -> Stop_i W b
 Real(dp), Private :: mg2, mb2, mt2, mgT, coeffT(-2:2), mSb2(2), mgsb(2) &
     & , mgSb2(2), mst2, mgT2
 Complex(dp), Private :: coeffSbB(3,0:2), coeffSbT(2,0:2,0:1)
 Logical :: WStB_contri(3)
 ! for check, if there is a numerical problem in the 3-body decays
 Real(dp), Private :: p_test 
 Real(dp), Private, Parameter :: prec=100._dp*Epsilon(1._dp)

Contains


 Subroutine GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
    & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_t, mf_d, USquark, gSU3            &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, cpl_UGSu_L, cpl_UGSu_R &
    & , cpl_SdSuW, DSquark, cpl_DNSd_L, cpl_DNSd_R, cpl_CUSd_L, cpl_CUSd_R     &
    & , cpl_DGSd_L, cpl_DGSd_R, mW, epsI, deltaM, Check_Real_States)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a gluino
 ! input:
 !    mN ......... neutralino masses
 !    mC ......... chargino masses
 !
 ! output:
 ! written by Werner Porod, 16.05.2001
 !  - taking the code from the Routine NeutralinoDecays.f as basis 
 ! 15.07.02: adding gluino -> W Stop_i B
 !           the method used in case of generation mixing is only valid
 !           if the third generation does not mix to strongly with other
 !           generations
 ! 24.08.03: - adding gluino -> chi^0_i gluon
 !           - adding possiblity to calculate 3-body states only via
 !             virtual states -> new logical variable Check_Real_States
 ! 18.09.2010: adapting to new variable type for particles
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_d, id_d(:), n_u, id_u(:), n_n, n_c, id_gl, n_su &
     & , n_sd, n_W
  Real(dp), Intent(in) :: mf_u(n_u), mf_d(n_d), epsI, deltaM, g_T, gSU3, mW(:)
  Complex(dp), Intent(in) :: cpl_UNSu_L(:,:,:), cpl_UNSu_R(:,:,:)          &
     & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:), cpl_UGSu_L(:,:)             &
     & , cpl_UGSu_R(:,:), cpl_DNSd_L(:,:,:), cpl_DNSd_R(:,:,:)             &
     & , cpl_CUSd_L(:,:,:), cpl_CUSd_R(:,:,:), cpl_DGSd_L(:,:)             &
     & , cpl_DGSd_R(:,:), cpl_SdSuW(:,:,:)
  Logical, Intent(in) :: Check_Real_States

  Type(particle2), Intent(in) :: Dsquark(:)
  Type(particle23), Intent(in) :: Usquark(:), Chi0(:), ChiPm(:)
  Type(particle23), intent(inout) :: Glu

  Real(dp) :: mUSquark(n_su), gUSquark(n_su), mDSquark(n_sd), gDSquark(n_sd)
  Integer :: i1, i2, i3, n_Sf4, n_CSf4, n_Sf8, i_part, i_c
  Real(dp) :: factor(2), mUsquark2(6), mDsquark2(6), gNff(3,3), gSbottom(2) &
     & , gCffp(3,3), mStop(2), mStop2(2), mSbottom(2), mSbottom2(2), gG     &
     & , g_SU(6), g_Sd(6)
  Real(dp) :: mN(n_n), mGlu, mC(n_c)
  Complex(dp) :: c_BGSb_L(2), c_BGSb_R(2), c_WSbSt(2,2), c_TGSt_L(2)       &
     & , c_TGSt_R(2)
  Real(dp), Allocatable :: IntegralsSf4(:,:)
  Complex(dp), Allocatable :: IntegralsCSf4(:,:), IntegralsSf8(:,:)
  logical :: check
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluinoThreeBodyDecays'

  !--------------------
  ! checking for model 
  !--------------------
   Allocate( IntegralsSf4(2500,10) )
   Allocate( IntegralsCSf4(2500,12) )
   Allocate( IntegralsSf8(2500,16) )

   IntegralsSf4 = 0._dp
   IntegralsCSf4 = 0._dp
   IntegralsSf8 = 0._dp

   n_Sf4 = 0
   n_CSf4 = 0
   n_Sf8 = 0

   If (Check_Real_States) then
    g_Sd = 0._dp
    g_Su = 0._dp
   else
    g_Sd = gDSquark
    g_Su = gUSquark
   end if
   check = Check_Real_States

   mGlu = Glu%m
   Glu%gi3 = 0._dp

   mDsquark = Dsquark%m
   mUsquark = Usquark%m
   mDsquark2 = Dsquark%m2
   mUsquark2 = Usquark%m2
   gDsquark = Dsquark%g
   gUsquark = Usquark%g

   mN = Chi0%m
   mC = ChiPm%m

   factor(1) = oo256pi3 / Abs(mglu)**3   ! for 3-body decays
   !--------------------------------------
   ! decays into a neutralino + 2 quarks
   !--------------------------------------
   i_c = 1
   Do i1 = 1, n_n
    If (Abs(mGlu).Gt.Abs(mN(i1))) Then
      Call GluToChi0qq(mGlu, i1,' u u ', mN, mf_u, mUSquark, g_Su             &
         & , cpl_UGSu_L, cpl_UGSu_R, cpl_UNSu_L, cpl_UNSu_R                   &
         & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8  &
         & , deltaM, epsI, GenerationMixing, check, factor(1), gNff) 
      Do i2=1,n_u
       Do i3=i2,n_u
        If (i2.eq.i3) then
         Glu%gi3(i_c) = gNff(i2,i3)
         Glu%id3(i_c,1) = Chi0(i1)%id
         Glu%id3(i_c,2) = id_u(i2)
         Glu%id3(i_c,3) = id_u(i2) + 1
         i_c = i_c +1
        Else
         Glu%gi3(i_c) = gNff(i2,i3)
         Glu%id3(i_c,1) = Chi0(i1)%id
         Glu%id3(i_c,2) = id_u(i2)
         Glu%id3(i_c,3) = id_u(i3) + 1
         Glu%gi3(i_c+1) = gNff(i3,i2)
         Glu%id3(i_c+1,1) = Chi0(i1)%id
         Glu%id3(i_c+1,2) = id_u(i3)
         Glu%id3(i_c+1,3) = id_u(i2) + 1
         i_c = i_c +2
        End If
       End Do
      End Do

      Call GluToChi0qq(mGlu, i1,' d d ', mN, mf_d, mDSquark, g_Sd             &
         & , cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R                   &
         & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8  &
         & , deltaM, epsI, GenerationMixing, check, factor(1), gNff)
      Do i2=1,n_d
       Do i3=i2,n_d
        If (i2.eq.i3) then
         Glu%gi3(i_c) = gNff(i2,i3)
         Glu%id3(i_c,1) = Chi0(i1)%id
         Glu%id3(i_c,2) = id_d(i2)
         Glu%id3(i_c,3) = id_d(i2) + 1
         i_c = i_c +1
        Else
         Glu%gi3(i_c) = gNff(i2,i3)
         Glu%id3(i_c,1) = Chi0(i1)%id
         Glu%id3(i_c,2) = id_d(i2)
         Glu%id3(i_c,3) = id_d(i3) + 1
         Glu%gi3(i_c+1) = gNff(i3,i2)
         Glu%id3(i_c+1,1) = Chi0(i1)%id
         Glu%id3(i_c+1,2) = id_d(i3)
         Glu%id3(i_c+1,3) = id_d(i2) + 1
         i_c = i_c +2
        End If
       End Do
      End Do

    End If

   End Do

   !--------------------------------------
   ! decay into charginos + 2 quarks
   !--------------------------------------
   Do i1=1,n_c
    If (Abs(mGlu).Gt.Abs(mC(i1))) Then
     Call GluToChimqqp(mGlu, i1, mC, mf_d, mf_u, mDSquark, g_Sd            &
          & , cpl_DGSd_L, cpl_DGSd_R, cpl_CUSd_L, cpl_CUSd_R, mUSquark      &
          & , g_Su, cpl_UGSu_L, cpl_UGSu_R, cpl_CDSu_L, cpl_CDSu_R          &
          & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8      &
          & , n_Sf8, deltaM, epsI, GenerationMixing, check, factor(1), gCffp)
     Do i2=1,n_d
      Do i3=1,n_u
       Glu%gi3(i_c) = gCffp(i2,i3)
       Glu%id3(i_c,1) = ChiPm(i1)%id
       Glu%id3(i_c,2) = id_d(i2)
       Glu%id3(i_c,3) = id_u(i3) + 1
       Glu%gi3(i_c+1) = gCffp(i2,i3)
       Glu%id3(i_c+1,1) = ChiPm(i1)%id + 1
       Glu%id3(i_c+1,2) = id_d(i2) + 1
       Glu%id3(i_c+1,3) = id_u(i3)
       i_c = i_c +2
      End Do
     End Do
    End If
   End Do
   !-------------------------------
   ! decay into neutralino + gluon
   !-------------------------------  
   factor(1) = - gSU3 * oo16pi2
   factor(2) = 0.125_dp / ( Pi * Abs(mGlu)**3 ) 
   Do i1=1,n_n
    If (Abs(mGlu).Gt.Abs(mN(i1)) ) Then
     Call GluToChi0Gluon(mGlu, i1, mN, mf_u, mUsquark2, cpl_UNSu_L           &
        & , cpl_UNSu_R, cpl_UGSu_L, cpl_UGSu_R, mf_d, mDsquark2, cpl_DNSd_L  &
        & , cpl_DNSd_R, cpl_DGSd_L, cpl_DGSd_R, factor, gG)
      Glu%gi2(100+i1) = gG
      Glu%id2(100+i1,1) = Chi0(i1)%id
      Glu%id2(100+i1,2) = id_gl
    End If
   End Do
   !-------------------------------
   ! decay into Stop W B
   ! in case of generation mixing, the position of the stops is determined
   ! using the couplings stop-top-gluino and in the same way for sbottoms
   !-------------------------------  
    mStop = mUSquark(5:6)
    mStop2 = mUSquark2(5:6)
    mSbottom = mDSquark(5:6)
    mSbottom2 = mDSquark2(5:6)
    gSbottom = g_Sd(5:6)
    c_BGSb_L = cpl_DGSd_L(3,5:6)
    c_BGSb_R = cpl_DGSd_R(3,5:6)
    c_TGSt_L = cpl_UGSu_L(3,5:6)
    c_TGSt_R = cpl_UGSu_R(3,5:6)
    c_WSbSt = cpl_SdSuW(5:6,5:6,1)
   i_part = 0
!   Call GluinoToStopWB(mglu, mStop, mStop2, mSbottom, mSbottom2, gSbottom     &
!       & , c_BGSb_L, c_BGSb_R, c_WSbSt, c_TGSt_L, c_TGSt_R, CKM(3,3), g_t &
!       & , epsI, i_part, Check_Real_States, gStWB(5:6,3) )
   !---------------------------
   ! simplifies life sometimes
   !---------------------------
   Glu%g = Sum(Glu%gi2) + Sum(Glu%gi3)
   If (Glu%g.ne.0._dp) then
    Glu%bi2 = Glu%gi2 / Glu%g
    Glu%bi3 = Glu%gi3 / Glu%g
   End If

  Deallocate( IntegralsSf4, IntegralsCSf4, IntegralsSf8 )

  Iname = Iname - 1
 End Subroutine GluinoThreeBodyDecays



 Subroutine GluToChi0qq(mGlu, i_out, state, mN, mf, mSf, gSf, cpl_FGSf_L   &
    & , cpl_FGSf_R, cpl_FNSf_L, cpl_FNSf_R, IntegralsSf4, n_Sf4            &
    & , IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI           &
    & , GenerationMixing, check, fac, gNff, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a gluino to a Neutralino + fermion pair
 ! Written by Werner Porod, 28.06.2001
 ! 24.08.03: adding the possiblity to check for real intermediate states
 !           with the help of the variable check
 !           If check=.True. then the contribution of real intermediate
 !           states will not be calculated
 !--------------------------------------------------------------------------
 Implicit None
  Character(len=5), Intent(in) :: state 
  Integer, Intent(in) :: i_out
  Logical, Intent(in) :: check
  Integer, Intent(inout) :: n_Sf4, n_CSf4, n_Sf8

  Real(dp), Intent(in) :: mGlu, mN(:), mf(:), mSf(:), gSf(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:)    &
                         & , cpl_FGSf_L(:,:), cpl_FGSf_R(:,:)
  Complex(dp), Intent(inout) :: IntegralsCSf4(:,:), IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gNff(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4, i_gen
  Real(dp) :: Boson2(2), mass(4), resR, resRa, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8), resCa

  Real(dp) :: gNffSum(3,3,80)
  Character(len=20) :: Contribution(3,3,80)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluToChi0qq'

  mass(1) = mGlu

  gNffSum = 0._dp
  Contribution = ' '

  Isum = 0
  !-----------------------------
  ! sfermion sfermion, diagonal
  !-----------------------------
  If (GenerationMixing) Then
   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)
    Do i1=1,3
     Do i3=1,3
      !------------------------
      ! t-channel
      !------------------------
      coup1(1) = cpl_FGSf_L(i1,i2)
      coup1(2) = cpl_FGSf_R(i1,i2)
      coup1(3) = Conjg(cpl_FNSf_R(i3,i_out,i2))
      coup1(4) = Conjg(cpl_FNSf_L(i3,i_out,i2))
      mass(2) = mf(i1)
      mass(3) = - mf(i3)
      mass(4) = mN(i_out)        
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      !------------------------
      ! u-channel
      !------------------------
      coup1(1) = Conjg(cpl_FGSf_R(i3,i2))  ! u-channel
      coup1(2) = Conjg(cpl_FGSf_L(i3,i2))
      coup1(3) = cpl_FNSf_L(i1,i_out,i2)
      coup1(4) = cpl_FNSf_R(i1,i_out,i2)
      mass(2) = mf(i3)
      mass(3) = - mN(i_out)        
      mass(4) = mf(i1)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                            &, IntegralsSf4, n_Sf4, resRa, check)
      gNffSum(i1,i3,Isum) = resR + resRa
      Contribution(i1,i3,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i3)
     End Do ! i3 fermions
    End Do  ! i1 fermions
   End Do   ! i2 sfermions

  Else
   Do i2=1,6  ! fermion generation
    i1 = (i2+1)/2
     Isum = Isum + 1
     Boson2(1) = mSf(i2)
     Boson2(2) = gSf(i2)
     !------------------------
     ! t-channel
     !------------------------
     coup1(1) = cpl_FGSf_L(i1,i2)
     coup1(2) = cpl_FGSf_R(i1,i2)
     coup1(3) = Conjg(cpl_FNSf_R(i1,i_out,i2))
     coup1(4) = Conjg(cpl_FNSf_L(i1,i_out,i2))
     mass(2) = mf(i1)
     mass(3) = - mf(i1)
     mass(4) = mN(i_out) 
     !------------------------
     ! t-channel
     !------------------------
     coup1(1) = cpl_FGSf_L(i1,i2)
     coup1(2) = cpl_FGSf_R(i1,i2)
     coup1(3) = Conjg(cpl_FNSf_R(i1,i_out,i2))
     coup1(4) = Conjg(cpl_FNSf_L(i1,i_out,i2))
     mass(2) = mf(i1)
     mass(3) = - mf(i1)
     mass(4) = mN(i_out) 
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resR, check)
     !------------------------
     ! u-channel
     !------------------------
     coup1(1) = Conjg(cpl_FGSf_R(i1,i2))
     coup1(2) = Conjg(cpl_FGSf_L(i1,i2))
     coup1(3) = cpl_FNSf_L(i1,i_out,i2)
     coup1(4) = cpl_FNSf_R(i1,i_out,i2)
     mass(2) = mf(i1)
     mass(3) = - mN(i_out)        
     mass(4) = mf(i1)
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resRa, check)
     gNffSum(i1,i1,Isum) = resR + resRa
     Contribution(i1,i1,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do   ! i1 u-quarks
  End If    ! GenerationMixing

  !----------------------------------
  ! Sfermion - Sfermion, non-diagonal
  !----------------------------------
  If (GenerationMixing) Then
   Do i3=1,5
    Boson4(1) = mSf(i3)
    Boson4(2) = gSf(i3)
    Do i4=i3+1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i4)
     Boson4(4) = gSf(i4)
     Do i1 = 1,3 ! fermions
      Do i2 = 1,3
      Contribution(i1,i2,Isum) = &
         &  'tt Sf_'//Bu(i3)//' Sf_'//Bu(i4)//' f_'//Bu(i1)//' f_'//Bu(i2)
       !-------------
       ! t-channel
       !-------------
       mass(2) = mf(i1)
       mass(3) = -mf(i2)
       mass(4) = mN(i_out)
       coup2(1) = cpl_FGSf_L(i1,i3)
       coup2(2) = cpl_FGSf_R(i1,i3)
       coup2(3) = Conjg(cpl_FGSf_R(i1,i4))
       coup2(4) = Conjg(cpl_FGSf_L(i1,i4))
       coup2(5) = Conjg(cpl_FNSf_R(i2,i_out,i3))
       coup2(6) = Conjg(cpl_FNSf_L(i2,i_out,i3))
       coup2(7) = cpl_FNSf_L(i2,i_out,i4)
       coup2(8) = cpl_FNSf_R(i2,i_out,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resC, check)
       !-------------
       ! u-channel
       !-------------
       mass(2) = mf(i2)
       mass(3) = -mN(i_out)
       mass(4) = mf(i1)
       coup2(1) = Conjg(cpl_FGSf_R(i2,i3))
       coup2(2) = Conjg(cpl_FGSf_L(i2,i3))
       coup2(3) = cpl_FGSf_L(i2,i4)
       coup2(4) = cpl_FGSf_R(i2,i4)
       coup2(5) = cpl_FNSf_L(i1,i_out,i3)
       coup2(6) = cpl_FNSf_R(i1,i_out,i3)
       coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i4))
       coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i4))
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resCa, check)
       gNffSum(i1,i2,Isum) = 2._dp * Real(resC+resCa,dp)
      End Do
     End Do
    End Do
   End Do

  Else ! .no.GenerationMixing

   Do i1 = 1,3
    i2 = 2*i1 - 1
    i3 = 2*i1
    Isum = Isum + 1
    Contribution(i1,i1,Isum) = &
        &  'tt Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
    Boson4(1) = mSf(i2)
    Boson4(2) = gSf(i2)
    Boson4(3) = mSf(i3)
    Boson4(4) = gSf(i3)
    !-------------
    ! t-channel
    !-------------
    mass(2) = mf(i1)
    mass(3) = -mf(i1)
    mass(4) = mN(i_out)
    coup2(1) = cpl_FGSf_L(i1,i2)
    coup2(2) = cpl_FGSf_R(i1,i2)
    coup2(3) = Conjg(cpl_FGSf_R(i1,i3))
    coup2(4) = Conjg(cpl_FGSf_L(i1,i3))
    coup2(5) = Conjg(cpl_FNSf_R(i1,i_out,i2))
    coup2(6) = Conjg(cpl_FNSf_L(i1,i_out,i2))
    coup2(7) = cpl_FNSf_L(i1,i_out,i3)
    coup2(8) = cpl_FNSf_R(i1,i_out,i3)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resC, check)
    !-------------
    ! u-channel
    !-------------
    mass(2) = mf(i1)
    mass(3) = -mN(i_out)
    mass(4) = mf(i1)
    coup2(1) = Conjg(cpl_FGSf_R(i1,i2))
    coup2(2) = Conjg(cpl_FGSf_L(i1,i2))
    coup2(3) = cpl_FGSf_L(i1,i3)
    coup2(4) = cpl_FGSf_R(i1,i3)
    coup2(5) = cpl_FNSf_L(i1,i_out,i2)
    coup2(6) = cpl_FNSf_R(i1,i_out,i2)
    coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
    coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resCa, check)
    gNffSum(i1,i1,Isum) = 2._dp * Real(resC+resCa,dp)
   End Do

  End If

  !-------------------------------------------------------
  ! Sfermion_{xyz}  t-channel - Sfermion_{xyz}  u-channel
  !-------------------------------------------------------
  If (GenerationMixing) Then
   Do i2= 1,6
    Boson4(1) = mSf(i2)
    Boson4(2) = gSf(i2)
    Do i3 = 1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i3)
     Boson4(4) = gSf(i3)
     Do i1 = 1,3 ! fermion
      coup2(1) = cpl_FGSf_L(i1,i2)
      coup2(2) = cpl_FGSf_R(i1,i2)
      coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
      coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
      Do i4=1,3
       Contribution(i1,i4,Isum) = &
          &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i4)
       mass(2) = mf(i4)
       mass(3) = -mN(i_out)
       mass(4) = mf(i1)
       coup2(3) = cpl_FGSf_L(i4, i3)
       coup2(4) = cpl_FGSf_R(i4, i3)
       coup2(5) = Conjg(cpl_FNSf_R(i4,i_out,i2))
       coup2(6) = Conjg(cpl_FNSf_L(i4,i_out,i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSf8, n_Sf8, resC, check)
        gNffSum(i1,i4,Isum) = - 2._dp * Real(resC,dp)
      End Do
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,3 !u-quarks
    i_gen = 2*i1-1
    Do i2 = i_gen, i_gen+1
     Boson4(1) = mSf(i2)
     Boson4(2) = gSf(i2)
     coup2(1) = cpl_FGSf_L(i1,i2)
     coup2(2) = cpl_FGSf_R(i1,i2)
     coup2(5) = Conjg(cpl_FNSf_R(i1,i_out,i2))
     coup2(6) = Conjg(cpl_FNSf_L(i1,i_out,i2))
     Do i3 = i_gen, i_gen+1
      Isum = Isum + 1
      Contribution(i1,i1,Isum) = &
         &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
      Boson4(3) = mSf(i3)
      Boson4(4) = gSf(i3)
      mass(2) = mf(i1)
      mass(3) = -mN(i_out)
      mass(4) = mf(i1)
      coup2(3) = cpl_FGSf_L(i1, i3)
      coup2(4) = cpl_FGSf_R(i1, i3)
      coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
      coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSf8, n_Sf8, resC, check)
      gNffSum(i1,i1,Isum) = - 2._dp * Real(resC,dp)
     End Do
    End Do
   End Do

  End If

  !----------
  ! Summing
  !----------
  gNff = 0._dp
  Do i1=1,3
   If (GenerationMixing) Then
    Do i2=1,3
     gNff(i1,i2) = Sum( gNffSum(i1,i2,1:Isum) )
     If (gNff(i1,i2).Lt.0._dp) Then
      p_test = Abs(gNff(i1,i2)) / Maxval(Abs(gNffSum(i1,i2,1:Isum)))
      gNff(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero

      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma(G -> Chi_'//Bu(i_out)//state//') < 0 :' &
      & ,i1,i2,gNff(i1,i2)
      Write(ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gNffSum(i1,i2,i3).Ne.0._dp) &
        &      Write(ErrCan,*) Contribution(i1,i2,i3),gNffSum(i1,i2,i3)
      End Do
     End If
    End Do

   Else
    gNff(i1,i1) = Sum( gNffSum(i1,i1,1:Isum) )
    If (gNff(i1,i1).Lt.0._dp) Then
     p_test = Abs(gNff(i1,i1)) / Maxval(Abs(gNffSum(i1,i1,1:Isum)))
     gNff(i1,i1) = 0._dp
     If (p_test.le.prec) cycle ! this is a numerical zero
     Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma( G -> Chi_'//Bu(i_out)//state//') < 0 :' &
      & ,i1,i1,gNff(i1,i1)
     Write(ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gNffSum(i1,i1,i3).Ne.0._dp) &
        & Write(ErrCan,*) Contribution(i1,i1,i3),gNffSum(i1,i1,i3)
     End Do
    End If
   End If
  End Do

  gNff = fac * gNff

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gNffSum = gNffSum * fac

   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
     Write (n_out,*) &
      & 'G -> Chi_'//Bu(i_out)//state//') :' &
      & ,i1,i2,gNff(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gNffSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gNffSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,3
     Write (n_out,*) &
      & 'G -> Chi_'//Bu(i_out)//state//') :' &
      & ,i1,i1,gNff(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gNffSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gNffSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Iname = Iname - 1

 End Subroutine GluToChi0qq

 Subroutine GluToChimqqp(mGlu, i_out, mC, mf, mfp, mSf, gSf, cpl_FGSf_L     &
    & , cpl_FGSf_R, cpl_CFpSf_L, cpl_CFpSf_R, mSfp, gSfp, cpl_FpGSfp_L      &
    & , cpl_FpGSfp_R, cpl_CFSfp_L, cpl_CFSfp_R, IntegralsSf4, n_Sf4         &
    & , IntegralsSfC4, n_SfC4, IntegralsSf8, n_Sf8, deltaM, epsI            &
    & , GenerationMixing, check, fac, gCffp, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a gluino to a Chargino + fermion pair
 ! Written by Werner Porod, 28.06.2001
 ! 24.08.03: adding the possiblity to check for real intermediate states
 !           with the help of the variable check
 !           If check=.True. then the contribution of real intermediate
 !           states will not be calculated
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_out
  Logical, Intent(in) :: check
  Integer, Intent(inout) :: n_Sf4, n_SfC4, n_Sf8  
  Real(dp), Intent(in) :: mGlu, mC(:), mf(:), mfp(:), mSf(:), gSf(:), mSfp(:) &
      & , gSfp(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: cpl_FGSf_L(:,:), cpl_FGSf_R(:,:)             &
      & , cpl_CFpSf_L(:,:,:), cpl_CFpSf_R(:,:,:), cpl_FpGSfp_L(:,:)       &
      & , cpl_FpGSfp_R(:,:), cpl_CFSfp_L(:,:,:), cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsSfC4(:,:), IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gCffp(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Real(dp) :: gCffpSum(3,3,80)
  Character(len=20) :: Contribution(3,3,80)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluToChimffp'

  mass(1) = mGlu
  gCffpSum = 0._dp
  Contribution = ' '
  Isum = 0

  !-------------------
  ! Sfp Sfp, diagonal
  !-------------------  
  If (GenerationMixing) Then
   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i1=1,3
     coup1(1) = cpl_FpGSfp_L(i1,i2)
     coup1(2) = cpl_FpGSfp_R(i1,i2)
     Do i3=1,3
      mass(2) = mfp(i1)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)
      coup1(3) = Conjg(cpl_CFSfp_R(i_out,i3,i2))
      coup1(4) = Conjg(cpl_CFSfp_L(i_out,i3,i2))
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gCffpSum(i3,i1,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do

  Else

   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)

    i1 = (i2+1)/2

    mass(2) = mfp(i1)
    mass(3) = -mf(i1)
    mass(4) = mC(i_out)
    coup1(1) = cpl_FpGSfp_L(i1,i2)
    coup1(2) = cpl_FpGSfp_R(i1,i2)
    coup1(3) = Conjg(cpl_CFSfp_R(i_out,i1,i2))
    coup1(4) = Conjg(cpl_CFSfp_L(i_out,i1,i2))
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
    gCffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !-------------------------
  ! Sf Sf, diagonal
  !-------------------------
  If (GenerationMixing) Then
   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    Do i1=1,3
     coup1(1) = Conjg(cpl_FGSf_R(i1,i2))
     coup1(2) = Conjg(cpl_FGSf_L(i1,i2))
     Do i3=1,3
      mass(2) = mf(i1)
      mass(3) = - mC(i_out)
      mass(4) = mfp(i3)
      coup1(3) = cpl_CFpSf_L(i_out,i3,i2)
      coup1(4) = cpl_CFpSf_R(i_out,i3,i2)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                            &, IntegralsSf4, n_Sf4, resR, check)
      gCffpSum(i1,i3,Isum) = resR
      Contribution(i1,i3,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do
  Else

   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mf(i1)
    mass(3) = - mC(i_out)
    mass(4) = mfp(i1)
    coup1(1) = Conjg(cpl_FGSf_R(i1,i2))
    coup1(2) = Conjg(cpl_FGSf_L(i1,i2))
    coup1(3) = cpl_CFpSf_L(i_out,i1,i2)
    coup1(4) = cpl_CFpSf_R(i_out,i1,i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resR, check)
    gCffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do
  End If

  !-----------------------------------------------
  ! Sfp_{xyz}  t-channel - Sfp_{xyz}  t-channel
  !-----------------------------------------------
  If (GenerationMixing) Then
   Do i1=1,5
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2=i1+1,6
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     Do i3 = 1,3
      coup2(1) = cpl_FpGSfp_L(i3,i1)
      coup2(2) = cpl_FpGSfp_R(i3,i1)
      coup2(3) = Conjg(cpl_FpGSfp_R(i3,i2))
      coup2(4) = Conjg(cpl_FpGSfp_L(i3,i2))
      Do i4=1,3
       mass(2) = mfp(i3)
       mass(3) = -mf(i4)
       mass(4) = mC(i_out)
       coup2(5) = Conjg(cpl_CFSfp_R(i_out,i4,i1))
       coup2(6) = Conjg(cpl_CFSfp_L(i_out,i4,i1))
       coup2(7) = cpl_CFSfp_L(i_out,i4,i2)
       coup2(8) = cpl_CFSfp_R(i_out,i4,i2)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsSfC4, n_SfC4, resC, check)
       gCffpSum(i4,i3,Isum) = 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1=1,5,2
    i2 = i1+1
    Isum = Isum + 1

    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)

    i3 = (i1+1)/2

    mass(2) = mfp(i3)
    mass(3) = -mf(i3)
    mass(4) = mC(i_out)

    coup2(1) = cpl_FpGSfp_L(i3,i1)
    coup2(2) = cpl_FpGSfp_R(i3,i1)
    coup2(3) = Conjg(cpl_FpGSfp_R(i3,i2))
    coup2(4) = Conjg(cpl_FpGSfp_L(i3,i2))
    coup2(5) = Conjg(cpl_CFSfp_R(i_out,i3,i1))
    coup2(6) = Conjg(cpl_CFSfp_L(i_out,i3,i1))
    coup2(7) = cpl_CFSfp_L(i_out,i3,i2)
    coup2(8) = cpl_CFSfp_R(i_out,i3,i2)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSfC4, n_SfC4, resC, check)
    gCffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
     & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If

  !------------------------------------------------
  ! Sfp_{xyz}  t-channel - Sf_{xyz} u-channel
  !------------------------------------------------
  If (GenerationMixing) Then
   Do i1= 1,6
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2 =1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 = 1,3
      coup2(1) = cpl_FpGSfp_L(i3,i1)
      coup2(2) = cpl_FpGSfp_R(i3,i1)
      Do i4 =1, 3
       mass(2) = mf(i4)
       mass(3) = -mC(i_out)
       mass(4) = mfp(i3)
       coup2(3) = cpl_FGSf_L(i4,  i2)
       coup2(4) = cpl_FGSf_R(i4,  i2)
       coup2(5) = Conjg(cpl_CFSfp_R(i_out,i4,i1))
       coup2(6) = Conjg(cpl_CFSfp_L(i_out,i4,i1))
       coup2(7) = Conjg(cpl_CFpSf_R(i_out, i3, i2))
       coup2(8) = Conjg(cpl_CFpSf_L(i_out, i3, i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSf8, n_Sf8, resC, check)
       gCffpSum(i4,i3,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1= 1,6
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    i3 = (i1+1)/2

    coup2(1) = cpl_FpGSfp_L(i3,i1)
    coup2(2) = cpl_FpGSfp_R(i3,i1)

    Do i2 =1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     mass(2) = mf(i3)
     mass(3) = -mC(i_out)
     mass(4) = mfp(i3)

     coup2(3) = cpl_FGSf_L(i3,  i2)
     coup2(4) = cpl_FGSf_R(i3,  i2)
     coup2(5) = Conjg(cpl_CFSfp_R(i_out,i3,i1))
     coup2(6) = Conjg(cpl_CFSfp_L(i_out,i3,i1))
     coup2(7) = Conjg(cpl_CFpSf_R(i_out, i3, i2))
     coup2(8) = Conjg(cpl_CFpSf_L(i_out, i3, i2))
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSf8, n_Sf8, resC, check)
     gCffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !-----------------------------------------------
  ! Sf_{xyz}  u-channel - Sf_{xyz}  u-channel
  !-----------------------------------------------
  If (GenerationMixing) Then

   Do i1=1,5
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Do i2=i1+1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 = 1,3
      coup2(1) = Conjg(cpl_FGSf_R(i3,  i1))
      coup2(2) = Conjg(cpl_FGSf_L(i3,  i1))
      coup2(3) = cpl_FGSf_L(i3,  i2)
      coup2(4) = cpl_FGSf_R(i3,  i2)
      Do i4=1,3
       mass(2) = mf(i3)
       mass(3) = - mC(i_out)
       mass(4) = mfp(i4)
       coup2(5) = cpl_CFpSf_L(i_out, i4, i1)
       coup2(6) = cpl_CFpSf_R(i_out, i4, i1)
       coup2(7) = Conjg(cpl_CFpSf_R(i_out, i4, i2))
       coup2(8) = Conjg(cpl_CFpSf_L(i_out, i4, i2))
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsSfC4, n_SfC4, resC, check)
       gCffpSum(i3,i4,Isum) = 2._dp * Real(resC,dp)
       Contribution(i3,i4,Isum) = &
        & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1=1,5,2
    Isum = Isum + 1
    i2 = i1+1
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    i3 = (i1+1)/2

    mass(2) = mf(i3)
    mass(3) = - mC(i_out)
    mass(4) = mfp(i3)

    coup2(1) = Conjg(cpl_FGSf_R(i3,  i1))
    coup2(2) = Conjg(cpl_FGSf_L(i3,  i1))
    coup2(3) = cpl_FGSf_L(i3,  i2)
    coup2(4) = cpl_FGSf_R(i3,  i2)
    coup2(5) = cpl_CFpSf_L(i_out, i3, i1)
    coup2(6) = cpl_CFpSf_R(i_out, i3, i1)
    coup2(7) = Conjg(cpl_CFpSf_R(i_out, i3, i2))
    coup2(8) = Conjg(cpl_CFpSf_L(i_out, i3, i2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSfC4, n_SfC4, resC, check)
    gCffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
      & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If


  !----------
  ! Summing
  !----------
  gCffp = 0._dp
  Do i1=1,3
   If (GenerationMixing) Then
    Do i2=1,3
     gCffp(i1,i2) = Sum( gCffpSum(i1,i2,1:Isum) )
     If (gCffp(i1,i2).Lt.0._dp) Then
      p_test = Abs(gCffp(i1,i2)) / Maxval(Abs(gCffpSum(i1,i2,1:Isum)))
      gCffp(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero
      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma(G -> Chi^-_'//Bu(i_out)//') < 0 :' &
      & ,i1,i2,gCffp(i1,i2)
      Write(ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffpSum(i1,i2,i3).Ne.0._dp) &
        &      Write(ErrCan,*) Contribution(i1,i2,i3),gCffpSum(i1,i2,i3)
      End Do
     End If
    End Do

   Else
    gCffp(i1,i1) = Sum( gCffpSum(i1,i1,1:Isum) )
    If (gCffp(i1,i1).Lt.0._dp) Then
     p_test = Abs(gCffp(i1,i1)) / Maxval(Abs(gCffpSum(i1,i1,1:Isum)))
     gCffp(i1,i1) = 0._dp
     If (p_test.le.prec) cycle ! this is a numerical zero
     Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma(G -> Chi^-_'//Bu(i_out)//') < 0 :' &
      & ,i1,i1,gCffp(i1,i1)
     Write(ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffpSum(i1,i1,i3).Ne.0._dp) &
        & Write(ErrCan,*) Contribution(i1,i1,i3),gCffpSum(i1,i1,i3)
     End Do
    End If
   End If
  End Do

  gCffp = fac * gCffp

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gCffpSum = gCffpSum * fac

   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
     Write (n_out,*) &
      & 'Gamma(G -> Chi^-_'//Bu(i_out)//') :' &
      & ,i1,i2,gCffp(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffpSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gCffpSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,3
     Write (n_out,*) &
      & 'Gamma(G -> Chi^-_'//Bu(i_out)//') :' &
      & ,i1,i1,gCffp(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffpSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gCffpSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Iname = Iname - 1

 End Subroutine GluToChimqqp
 
 Subroutine GluinoToStopWB(mglu, mStop, mStop2, mSbottom, mSbottom2, gSbottom &
                        &, c_BGSb_L, c_BGSb_R, c_WSbSt, c_TGSt_L, c_TGSt_R    &
                        &, c_Wbt, g_t, prec, i_part, RealStates, width)
 !--------------------------------------------------------------------------
 ! calculates the width gluino -> stop_1 b W in the MSSM
 ! written by Werner Porod
 ! - 11.02.02: implementation without generation mixing
 ! - 15.07.02: addapting code for implementation
 !             in routine GluinoThreeBodyDecays
 ! 24.08.03: adding the possiblity to check for real intermediate states
 !           with the help of the variable RealStates
 !           If RealStates=.True. then the contribution of real intermediate
 !           states will not be calculated
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: i_part
  Real(dp), Intent(in) :: mglu, mStop(2), mStop2(2), mSbottom(2)   &
      & , mSbottom2(2) , gSbottom(2), prec, g_t
  Complex(dp), Intent(in) :: c_BGSb_L(2), c_BGSb_R(2), c_WSbSt(2,2) &
      & , c_TGSt_L(2), c_TGSt_R(2), c_Wbt
  Logical, Intent(in) :: RealStates
  Real(dp), Intent(out) ::  width(2)

  Integer :: i1, i2
  Real(dp) :: mBmG, AbsCL2, AbsCR2, mTmG, mb2omW2, ReClCr, mBmT, mst2omW2 &
      & , smin, smax

  !----------------------------
  ! first kinematics
  !----------------------------
  width = 0._dp
  If (mglu.Lt.(mStop(1)+mf_d(3)+mW)) Return
 
  Iname = Iname + 1
  NameOfUnit(Iname) = "GluinoToStopWB"

  mg2 = mglu**2
  mb2 = mf_d2(3)
  mBmG = mglu * mf_d(3)
  mTmG = mf_u(3) * mglu
  mb2omW2 = mb2 / mW2
  mBmT = mf_d(3) * mf_u(3)

  mt2 = mf_u2(3)
  mst2omW2 = mst2 / mW2
  msb2 = mSbottom2

  Do i2=1,2
   WStB_contri = .True.
   If (RealStates) Then
    mgt = 0._dp
    mgsb = 0._dp
    If (mGlu.Gt.(mStop(i2)+mf_u(3))) WStB_contri(1) = .False.
    If (mGlu.Gt.(mSbottom(1)+mf_d(3))) WStB_contri(2) = .False.
    If (mGlu.Gt.(mSbottom(2)+mf_d(3))) WStB_contri(3) = .False.
    If (      (.Not.WStB_contri(1)).And.(.Not.WStB_contri(2))   &
       & .And.(.Not.WStB_contri(3)) ) Cycle
    mgt = mf_u(3) * g_t
    mgsb = gSbottom * mSbottom
   End If
   mgt2 = mgt**2
   mgsb2 = mgsb**2
   If (mglu.Lt.(mStop(i2)+mf_d(3)+mW)) Return
   mst2 = mStop2(i2)
   !--------------------------------------
   ! coefficients for sbottom exchange
   !--------------------------------------
   coeffSbB = ZeroC
   Do i1=1,2
    AbsCL2 = Abs(c_BGSb_L(i1))**2 + Abs(c_BGSb_R(i1))**2
    ReClCr = Real(c_BGSb_L(i1)* Conjg(c_BGSb_R(i1)),dp )
    coeffSbB(i1,2) = - 2._dp * AbsCL2
    coeffSbB(i1,1) = ( 2._dp * (mb2 + mG2 - mst2) + mW2) * AbsCL2        &
                 & + 8._dp * mBmG * ReClCr
    coeffSbB(i1,0) = ( (mb2 + mG2) * AbsCL2 + 4._dp * mBmG * ReClCr )    &
                   & * ( 2._dp * mst2 - mW2)
    coeffSbB(i1,:) =  coeffSbB(i1,:) * Abs(c_WSbSt(i1,i2))**2
   End Do
   coeffSbB(3,2) = - 2._dp * ( Conjg(c_BGSb_L(1))*c_BGSb_L(2)            &
                &           + Conjg(c_BGSb_R(1))*c_BGSb_R(2) )
   coeffSbB(3,1) = ( Conjg(c_BGSb_L(1))*c_BGSb_L(2)                      &
              &     + Conjg(c_BGSb_R(1))*c_BGSb_R(2) )                   &
              &   * (2._dp * (mb2 + mG2 - mst2 ) + mW2 )                 &
              &   + 4._dp * mBmG * ( Conjg(c_BGSb_L(1))*c_BGSb_R(2)      &
              &                    + Conjg(c_BGSb_R(1))*c_BGSb_L(2) )
   coeffSbB(3,0) = ( ( Conjg(c_BGSb_L(1))*c_BGSb_L(2)                    &
              &     + Conjg(c_BGSb_R(1))*c_BGSb_R(2) ) * (mb2 + mg2)     &
              &   + 2._dp * mBmG * ( Conjg(c_BGSb_L(1))*c_BGSb_R(2)      &
              &                    + Conjg(c_BGSb_R(1))*c_BGSb_L(2) ) )  &
              & * ( 2._dp * mst2 - mW2)
   coeffSbB(3,:) = 2._dp * coeffSbB(3,:) * c_WSbSt(1,i2) *Conjg(c_WSbSt(2,i2)) 

   !--------------------------------------
   ! coefficients for top exchange
   !--------------------------------------
   coeffT = 0._dp
   AbsCR2 = Abs(c_TGSt_R(i2))**2
   AbsCL2 = Abs(c_TGSt_L(i2))**2
   ReClCr = 2._dp * Real( Conjg(c_TGSt_L(i2)) * c_TGSt_R(i2),dp )
   coeffT(-2) = AbsCL2 * mt2 * ( 0.5_dp * mb2 * (1._dp + mb2omW2) - mW2)     &
           &          * (mg2 - mst2) 
   coeffT(-1) = AbsCR2 * ( 0.5_dp * mb2 * ( 1._dp + mb2omW2) - mW2)          &
           &          * (mg2 - mst2)                                         & 
           & + AbsCL2 * mt2 * (0.5_dp * mb2 * (1._dp  + mb2omW2 ) - mW2      &
           &                  + (mg2 - mst2) * (0.5_dp - mb2omW2 )  )        &
           & + ReClCr * mTmG * ( mb2 * (1._dp + mb2omW2) - 2._dp * mW2 )
   coeffT(0) = AbsCR2 * (0.5_dp * mb2 * (1._dp + mb2omW2) - mW2              &
          &            + (mg2-mst2) * (0.5_dp - mb2omW2) )                   &
          & + AbsCL2 * mt2 * ( 0.5_dp - mb2omW2 + 0.5_dp * (mg2-mst2) / mW2) &
          & + ReClCr * mTmG * (1._dp - 2._dp * mb2omW2)
   coeffT(1) = AbsCR2 * ( 0.5_dp * (1._dp + mg2/mW2 - mst2omW2) - mb2omW2 )  &
          & + 0.5_dp * ReClCr * mTmG / mW2 + 0.5_dp * AbsCL2 * mt2 / mW2
   coeffT(2) = 0.5_dp * AbsCR2 / mW2 
   coeffT = Abs(c_Wbt)**2 * coeffT
   !-----------------------------------------------------
   ! coefficients for interference top-sbottom exchange
   !-----------------------------------------------------
   coeffSbT = ZeroC
   Do i1=1,2
    coeffSbT(i1,2,1) = - c_TGSt_R(i2) * Conjg(c_BGSb_R(i1)) / mW2
    coeffSbT(i1,2,0) = mBmT * c_TGSt_L(i2) * Conjg(c_BGSb_L(i1)) / mW2
    coeffSbT(i1,1,1) = mBmT * c_TGSt_L(i2) * Conjg(c_BGSb_L(i1)) / mW2 &
                  & + c_TGSt_R(i2) * Conjg(c_BGSb_R(i1))               &
                  &   * ( 1._dp + mst2omW2 + mg2 / mW2)                      &
                  & + mBmG * c_TGSt_R(i2) * Conjg(c_BGSb_L(i1)) / mW2  &
                  & + mTmG * c_TGSt_L(i2) * Conjg(c_BGSb_R(i1)) / mW2
    coeffSbT(i1,1,0) = - c_TGSt_L(i2) * Conjg(c_BGSb_L(i1))           &
                  &     * mBmT * (1 + mb2omW2 + 2._dp * mst2omW2)           &
                  & + c_TGSt_R(i2) * Conjg(c_BGSb_L(i1))              &  
                  &     * mBmG * (1._dp - mb2omW2)                          &
                  & + c_TGSt_L(i2) * Conjg(c_BGSb_R(i1))              &  
                  &     * mTmG * (1._dp - mb2omW2)                          &
                  & + c_TGSt_R(i2) * Conjg(c_BGSb_R(i1))              &
                  &   * (mg2 - mst2) * (2._dp - mb2omW2)
    coeffSbT(i1,0,1) = c_TGSt_L(i2) * Conjg(c_BGSb_L(i1))             &
                  &     * mBmT * (1 - mst2omW2)                             &
                  & + c_TGSt_R(i2) * Conjg(c_BGSb_L(i1))              &  
                  &     * mBmG  * (1 - mst2omW2)                            &
                  & + c_TGSt_L(i2) * Conjg(c_BGSb_R(i1))              &  
                  &     * mTmG  * (1 - mst2omW2)                            &
                  & + c_TGSt_R(i2) * Conjg(c_BGSb_R(i1))              &
                  &   * mg2 * (1 - mst2omW2)
    coeffSbT(i1,0,0) = c_TGSt_L(i2) * Conjg(c_BGSb_L(i1))* mBmT         &
                &  * (mb2 * (1._dp+mst2omW2)-2._dp*mg2-mst2*(1._dp-mst2omW2)) &
                & + c_TGSt_R(i2) * Conjg(c_BGSb_L(i1)) * mBmG           &
                &       * (mb2*(1 + mst2omW2)-2._dp * mg2 + mst2 - mW2)       &
                & + c_TGSt_L(i2) * Conjg(c_BGSb_R(i1)) * mTmG           &
                &       * (mb2*(1 + mst2omW2)-2._dp * mg2 + mst2 - mW2)       &
                & + c_TGSt_R(i2) * Conjg(c_BGSb_R(i1))                  &
                &   * ( mg2 * ( mb2 * (1 + mst2omW2)                          &
                &             - 2._dp * (mg2 + mW2 - mst2) )                  &
                &     + mb2 * mst2 * (1._dp - mst2omW2) )
    coeffSbT(i1,:,:) = c_WSbSt(i1,i2) * coeffSbT(i1,:,:)
   End Do

   coeffSbT = - 2._dp * coeffSbT * Conjg(c_Wbt)

   smin = (mf_d(3)+mW)**2
   smax = (mglu - mstop(i2))**2
   If (i_part.Eq.1) Then
    CoeffSbT = 0
    coeffT = 0
   Else If (i_part.Eq.2) Then
    coeffSbB = 0
    CoeffSbT = 0
   Else If (i_part.Eq.3) Then
    coeffSbB = 0
    coeffT = 0
   End If
   width(i2) = oo256pi3 * Dgauss(GluStoWB,smin,smax,prec) / Abs(mglu)**3
  End Do

  Iname = Iname - 1

 End Subroutine  GluinoToStopWB
!\end{verbatim}

 Real(dp) Function GluStoWB(s)
 Implicit None
  Real(dp), Intent(in) :: s

  Integer :: i1
  Real(dp) :: s2, wsr, tmin, tmax, teil(6), sumI

  s2 = s**2
  wsr = kappa(s,mst2,mg2) * kappa(s, mW2, mb2)
  sumI = mg2 + mW2 + mst2 + mb2 - s - (mg2-mst2)*(mb2-mW2) / s
  tmin = 0.5_dp * (sumI - wsr/s)
  tmax = 0.5_dp * (sumI + wsr/s)
  teil = 0._dp
  !------------
  ! sbottom 
  !------------
  Do i1=1,2
   If (WStB_contri(i1+1)) Then
    teil(i1) = coeffSbB(i1,2) * I2t2dtm1g1(tmin, tmax, mSb2(i1), mgSb2(i1))   &
           & + coeffSbB(i1,1) * I2tdtm1g1(tmin, tmax, mSb2(i1), mgSb2(i1))    &
           & + coeffSbB(i1,0) * I2dtm1g1(tmin, tmax, mSb2(i1), mgSb2(i1))
   End If
  End Do
  If (WStB_contri(2).And.WStB_contri(3)) Then
   teil(3) = Real( coeffSbB(3,2)                                              &
          &    *I2t2dtm12g12(tmin, tmax, mSb2(1), mgSb(1), mSb2(2), mgSb2(2)) &
          &   + coeffSbB(3,1)                                                 &
          &    * I2tdtm12g12(tmin, tmax, mSb2(1), mgSb(1), mSb2(2), mgSb2(2)) &
          &   + coeffSbB(3,0)                                                 &
          &     * I2dtm12g12(tmin, tmax, mSb2(1), mgSb(1), mSb2(2), mgSb2(2)) &
          &     ,dp)
  End If
  !------
  ! top
  !------
  If (WStB_contri(1)) Then
   teil(4) = wsr * ( coeffT(-2) / s2 + coeffT(-1) / s + coeffT(0) &
          &       + coeffT(2) * s2 + coeffT(1) * s ) / ((s-mt2)**2 + mgt2)
  End If
  !----------------
  ! top sbottom
  !----------------
  Do i1=1,2
   If (WStB_contri(1).And.WStB_contri(i1+1)) Then
    teil(4+i1) = (coeffSbT(i1,0,1) * s + coeffSbT(i1,0,0) )                   &
            &      * I2dsmg1tmg2(s, tmin, tmax, mt2, mgt, mSb2(i1), mgSb(i1)) &
            & + (coeffSbT(i1,1,1) * s + coeffSbT(i1,1,0) )                    &
            &     * I2tdsmg1tmg2(s, tmin, tmax, mt2, mgt, mSb2(i1), mgSb(i1)) &
            & + (coeffSbT(i1,2,1) * s + coeffSbT(i1,2,0) )                    &
            &     * I2t2dsmg1tmg2(s, tmin, tmax, mt2, mgt, mSb2(i1), mgSb(i1))
   End If
  End Do

  GluStoWB = Sum(teil)

  !----------------------------------------------------------------------
  ! parts of the amplitudes squared need to be positive
  !----------------------------------------------------------------------
  If (   (teil(1).Lt.0._dp).Or.(teil(2).Lt.0._dp).Or.(teil(4).Lt.0._dp)  &
   & .Or.(Sum(teil(1:3)).Lt.0._dp).Or.(GluStoWB.Lt.0._dp) ) Then
    If (ErrorLevel.Ge.-1) Then
     Write(ErrCan,*) "Error in Function GluStoWB:"
     Write(ErrCan,*) "teil(1),teil(2),teil(4)",teil(1:2),teil(4)
     Write(ErrCan,*) "Sum(teil(1:3))",Sum(teil(1:3))
     Write(ErrCan,*) "GluStoWB",GluStoWB
     Write(ErrCan,*) "s,mg2,mst2,mb2,mw2",s,mg2,mst2,mb2,mw2
     If (ErrorLevel.Ge.0) Call TerminateProgram
    End If
  End If

 End Function GluStoWB

 Real(dp) Function GluinoToStopC(mglu, mStop2, Rstop, mSbottom2, mP02 &
        & , tanb, yukB, mu, A_b, c_GUSu_R, kont)
 Implicit None
  Integer, Intent(inout) :: kont
  Complex(dp), Intent(in) :: yukB, mu, A_b, Rstop(6,6), c_GUSu_R(3,6)
  Real(dp), Intent(in) :: tanb, mSbottom2(2), mP02(:), mstop2(6), mglu

  Real(dp) :: sinb2, cosb2, cos2b, fakt16pi, mass2(3), test(2), Q2
  Complex(dp) :: deltaL, deltaR, mat3(3,3), Rsf(3,3), coupL

  cosb2 = 1._dp / (1._dp + tanb**2)
  sinb2 = tanb**2 * cosb2
  cos2b = cosb2 - sinb2

  Q2 = GetRenormalizationScale()
  fakt16pi = Abs(YukB)**2 * oo16pi2 * CKM(2,3) * CKM(3,3) * Log(1.e32_dp/Q2) &
          & / cosb2


  ! corrected version with electroweak symmetry breaking
   deltal = - fakt16pi * (Sum(mSbottom2) + sinb2 * mP02(2) - Abs(mu)**2   &
          &              - 0.5_dp * cos2b * Mz2 + Abs(A_b/yukB)**2 )
   deltar = fakt16pi * A_b * mf_u(3) / yukB ! other convention than Hikasa

   mat3(1,1) = mstop2(5)
   mat3(1,2) = 0._dp
   mat3(1,3) = deltal * Rstop(5,5) + deltar * Rstop(5,6)
   mat3(2,1) = Conjg(mat3(1,2))
   mat3(2,2) = mstop2(6)
   mat3(2,3) = - deltal * Rstop(6,5) + deltar * Rstop(6,6)
   mat3(3,1) = Conjg( mat3(1,3) )
   mat3(3,2) = Conjg( mat3(2,3) )
   mat3(3,3) = mstop2(4)

   Call EigenSystem(mat3, mass2, Rsf, kont,test)
   If (kont.Eq.-14) then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    If (ErrorLevel.Eq.2) Call TerminateProgram
    kont = 0
   end if
   If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
!   Call AddNOW() 
    Write(ErrCan,*) "Warning, in subroutine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Diagonalization of mat3 has failed",kont
    If (ErrorLevel.Eq.2) Call TerminateProgram
   End If

   test(1) = Maxval( Abs(Rsf(3,:) ) )
   If (test(1).Eq.Abs(Rsf(3,1))) Then
    coupL = Rsf(3,2) * c_GUSu_R(2,4)
   Else
    coupL = Rsf(3,1) *  c_GUSu_R(2,4)
   End If
   Call FermionToFermionScalar(mglu, mF_u(2), Sqrt(mStop2(5)), coupL, ZeroC &
                          &, GluinoToStopC )

 End Function GluinoToStopC
 Subroutine GluToChi0Gluon(mGlu, i_out, mN, mf_u, mUsquark2, cpl_UNSu_L     &
       & , cpl_UNSu_R, cpl_UGSu_L, cpl_UGSu_R, mf_d, mDsquark2, cpl_DNSd_L  &
       & , cpl_DNSd_R, cpl_DGSd_L, cpl_DGSd_R, factor, gGluon)
 !---------------------------------------------------------------------
 ! Calculates the decay of a gluino to a neutralino + gluon
 ! Written by Werner Porod, 28.06.2001
 !---------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) ::  i_out
  Real(dp), Intent(in) :: mN(:), factor(2), mf_u(3), mUsquark2(6), mf_d(3) &
     & , mDsquark2(6), mGlu
  Complex(dp), Intent(in) :: cpl_UNSu_L(:,:,:), cpl_UNSu_R(:,:,:)  &
     & , cpl_UGSu_L(:,:), cpl_UGSu_R(:,:), cpl_DNSd_L(:,:,:)       &
     & , cpl_DNSd_R(:,:,:), cpl_DGSd_L(:,:), cpl_DGSd_R(:,:)
  Real(dp), Intent(out) :: gGluon

  Integer :: i2, i3, i_gen
  Real(dp) :: mj2, mi2, m12, m22
  Complex(dp) :: Gcoup, Iinte, Jinte, Kinte, coup1, coup2, I2inte

  mj2 = mGlu**2
!  factor(1) = - gSU3 * oo16pi2
!  factor(2) = 0.5_dp / ( Pi * Abs(mGlu)**3 ) 

  mi2 = mN(i_out)**2
  Gcoup = ZeroC
  !--------------------
  ! up-squark up-quark
  !--------------------
  Do i2=1,3
   m12 = mf_u(i2)**2
   If (GenerationMixing) Then
    Do i3=1,6
     m22 = mUsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_UGSu_R(i2,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UGSu_L(i2,i3) ) * cpl_UNSu_L(i2,i_out,i3)
     coup2 = cpl_UGSu_L(i2,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UGSu_R(i2,i3) ) * cpl_UNSu_L(i2,i_out,i3) 
     Gcoup = Gcoup + coup1 * mGlu * (I2inte - Kinte)     &
         &         - Conjg(coup1) * mN(i_out) * Kinte    &
         &         + coup2 * mf_u(i2) * Iinte 
    End Do    

   Else
    i_gen = 2*i2-1
    Do i3=i_gen,i_gen+1
     m22 = mUsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_UGSu_R(i2,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UGSu_L(i2,i3) ) * cpl_UNSu_L(i2,i_out,i3)
     coup2 = cpl_UGSu_L(i2,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
         & - Conjg( cpl_UGSu_R(i2,i3) ) * cpl_UNSu_L(i2,i_out,i3) 
     Gcoup = Gcoup + coup1 * mGlu * (I2inte - Kinte)     &
         &         - Conjg(coup1) * mN(i_out) * Kinte    &
         &         + coup2 * mf_u(i2) * Iinte
    End Do    
   End If
  End Do    
  !------------------------
  ! down-squark down-quark
  !------------------------
  Do i2=1,3
   m12 = mf_d(i2)**2
   If (GenerationMixing) Then
    Do i3=1,6
     m22 = mDsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_DGSd_R(i2,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DGSd_L(i2,i3) ) * cpl_DNSd_L(i2,i_out,i3)
     coup2 = cpl_DGSd_L(i2,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DGSd_R(i2,i3) ) * cpl_DNSd_L(i2,i_out,i3) 
     Gcoup = Gcoup + coup1 * mGlu * (I2inte - Kinte)     &
         &         - Conjg(coup1) * mN(i_out) * Kinte    &
         &         + coup2 * mf_d(i2) * Iinte
    End Do    

   Else
    i_gen = 2*i2-1
    Do i3=i_gen,i_gen+1
     m22 = mDsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_DGSd_R(i2,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DGSd_L(i2,i3) ) * cpl_DNSd_L(i2,i_out,i3)
     coup2 = cpl_DGSd_L(i2,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DGSd_R(i2,i3) ) * cpl_DNSd_L(i2,i_out,i3) 
     Gcoup = Gcoup + coup1 * mGlu * (I2inte - Kinte)     &
         &         - Conjg(coup1) * mN(i_out) * Kinte    &
         &         + coup2 * mf_d(i2) * Iinte
    End Do    
   End If
  End Do    

  gGluon = factor(2) * (mj2 - mi2)**3 * Abs(factor(1)*Gcoup)**2   

 End Subroutine GluToChi0Gluon


End Module Gluino3Decays

