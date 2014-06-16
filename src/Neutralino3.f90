Module Neut3Decays

! load modules
Use Control
Use LoopFunctions
Use ThreeBodyPhaseSpace

! for check, if there is a numerical problem in the 3-body decays
 Real(dp), Private :: p_test 
 Real(dp), Private, Parameter :: prec=100._dp*Epsilon(1._dp)

! load modules
Contains


 Subroutine NeutralinoThreeBodyDecays(n_in, n_l, id_l, n_nu, id_nu, n_d, id_d  &
    & , n_u, id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su    &
    & , n_S0, n_P0, n_Spm, id_ph, Chi0, mZ, gZ, L_nu, R_nu, mf_l, L_e, R_e     &
    & , mf_u, L_u, R_u, mf_d, L_d, R_d, cpl_NNZ_L, cpl_NNZ_R, ChiPm, mW, gW    &
    & , cpl_NuLW, cpl_UDW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L, cpl_CNW_R, Spm     &
    & , cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L, cpl_SmpLNu_R, cpl_SmpDU_L      &
    & , cpl_SmpDU_R, S0, cpl_NNS0_L, cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R        &
    & , cpl_DDS0_L, cpl_DDS0_R, cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R &
    & , P0, cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
    & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, Sup        &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, Glu, cpl_UGSu_L        &
    & , cpl_UGSu_R, Sdown, cpl_DNSd_L, cpl_DNSd_R, cpl_CUSd_L, cpl_CUSd_R      &
    & , cpl_DGSd_L, cpl_DGSd_R, Sneut, cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L    &
    & , cpl_CLSn_R, Slept, cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R    &
    & , gSU2, sW2, GenerationMixing, OnlySM, epsI, deltaM, Check_Real_States   &
    & , M_D, n_int_diag, n_int_off_diag, is_NMSSM)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a neutralino into fermions +
 ! the decay into another neutralino + photon, in case of spontaneous R-parity
 ! violation also the decays into majoron are included. The idea is that the
 ! couplings are calcualted in another subroutine and that this routine
 ! can be used in a model independent way.
 ! input:
 !    n_in ....... index of the decaying neutralino, in case n_in < 0,
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !    mC ......... chargino masses
 !
 ! output:
 ! 2-body decay modes
 !  gPhoton(i) .......... neutralino_i + photon
 ! 3-body decay modes
!  gNll(i,j,k) ......... Neutralino_i + lepton_j + lepton_k (j,k=1 e
 !                                                            j,k=2 mu
 !                                                            j,k=3 tau)
 !  gNqq(i,j) ......... Neutralino_i + quark_j + quark_j (j,k=1 u
 !                                                        j,k=2 d
 !                                                        j,k=3 c
 !                                                        j,k=4 s
 !                                                        j,k=5 t
 !                                                        j,k=6 b)
 !  gCln(i,j) ........... Chargino_i + l_j nu (j=1 -> e nu
 !                                              j=2 -> mu nu
 !                                              j=3 -> tau nu )
 !  gCDU(i,j) .......... Chargino_i + q_j q'_j (j=1 -> u d
 !                                                j=2 -> c s
 !                                                j=3 -> t b )
 !  gP(i) ............... the above described partial widths in one
 !                        array, useful for further processing
 !  gT .................. total width
 !  BR(i) ............... branching ratios, optional
 ! written by Werner Porod, 16.05.2001
 !  - taking the code from the Routine NeutralinoDecays.f as basis 
 !  02.07.2001: new variable OnlySM controls if also decays into
 !              MSSM chargino/neutralinos should be calculated
 !  11.09.03: splitting up neutrino final states because this simplifies
 !            the Les Houches Interface
 !  12.09.03: including the possiblity to check for real intermediates states
 !  13.09.03: adding charged conjugated states to output
 !  06.01.06: - replacing unused gMajaron (is calculated in routine for
 !              2-body decays) by m_d. The later one is the minimum phase
 !              space for the calculation of 3-body decays to avoid numerical
 !              problems 
 !            - adding n_int_diag, n_int_off_diag, which give the size of
 !              arrays for the intermediate contrictuions
 !  19.09.2010: adjusting for new variable types
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_in, n_int_diag, n_int_off_diag, n_nu, n_l, n_d, n_u &
     & , n_Z, n_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_s0, n_p0, n_Spm, id_ph
  Integer, Intent(in) :: id_nu(:), id_l(:), id_d(:), id_u(:), id_Z(:), id_W(:)
  Real(dp), Intent(in) :: mZ, gZ, L_nu, R_nu, mf_l(:), L_e, R_e, mf_u(:), L_u  &
     & , R_u, mf_d(:), L_d, R_d, mW, gW, gSU2, sW2, epsI, deltaM, M_D
  Complex(dp), Intent(in) :: cpl_NNZ_L(:,:), cpl_NNZ_R(:,:), cpl_CCZ_L(:,:)  &
     & , cpl_CCZ_R(:,:), cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_UDW(3,3)        &
     & , cpl_SmpCN_L(:,:,:), cpl_SmpCN_R(:,:,:), cpl_SmpLNu_L(:,:,:)         &
     & , cpl_SmpLNu_R(:,:,:), cpl_SmpDU_L(:,:,:), cpl_SmpDU_R(:,:,:)         &
     & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_CCS0_L(:,:,:)             &
     & , cpl_CCS0_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)             &
     & , cpl_LLS0_L(:,:,:), cpl_LLS0_R(:,:,:), cpl_UUS0_L(:,:,:)             &
     & , cpl_UUS0_R(:,:,:), cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:)             &
     & , cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:), cpl_DDP0_L(:,:,:)             &
     & , cpl_DDP0_R(:,:,:), cpl_LLP0_L(:,:,:), cpl_LLP0_R(:,:,:)             &
     & , cpl_UUP0_L(:,:,:), cpl_UUP0_R(:,:,:), cpl_UNSu_L(:,:,:)             &
     & , cpl_UNSu_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)             &
     & , cpl_UGSu_L(:,:), cpl_UGSu_R(:,:), cpl_DNSd_L(:,:,:)                 &
     & , cpl_DNSd_R(:,:,:), cpl_CUSd_L(:,:,:), cpl_CUSd_R(:,:,:)             &
     & , cpl_DGSd_L(:,:), cpl_DGSd_R(:,:), cpl_NuNSn_L(:,:,:)                &
     & , cpl_NuNSn_R(:,:,:), cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)            &
     & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:), cpl_CNuSl_L(:,:,:)            &
     & , cpl_CNuSl_R(:,:,:), cpl_NuLW(:,:)
  Logical, Intent(in) :: GenerationMixing, OnlySM, Check_Real_States
  Logical, Optional , Intent(in) :: is_NMSSM
  Type(particle2), Intent(in) :: Sdown(:), Spm(:), P0(:)
  Type(particle23), Intent(in) :: Sneut(:), Slept(:), Sup(:), ChiPm(:), S0(:) &
     & , Glu
  Type(particle23), Intent(inout) :: Chi0(:)

  Integer :: i_start, i_end, i_run, i_c, i1, i2, i3, n_Z4, n_S04, n_Sf4   &
     & , n_ZS04, n_ZSf8, n_S0P04, n_S0Sf8, n_CSf4, n_Sf8, n_W4, n_WSpm    &
     & , n_ZS08, n_S0P08, n_W8, n_ZW8, n_WSpm8, n_WSf8, n_Z8
  Real(dp) :: factor(3), mW2a, m_nu(3), mUsquark2(n_Su), gP0_in(n_P0)       &
     & , mDsquark2(n_Sd), mSlepton2(n_sle), gNff(3,3), gCffp(3,3), mGlu     &
     & , mSneutrino(n_snu), mN(n_n), mC(n_c), mP0(n_P0), mUsquark(n_Su)     &
     & , mDsquark(n_Sd), mSlepton(n_sle), mS0(n_S0), mSpm(n_Spm), diffM     &
     & , gNNN, gNCC, gGqq(3,3), gW_in, gZ_in, mSpm2(n_Spm), gS0_in(n_S0)    &
     & , gSpm_in(n_Spm), g_Su(n_Su), g_Sd(n_Sd), g_Sl(n_Sle), g_Sn(n_Snu)
  Real(dp), Allocatable :: IntegralsZ4(:,:), IntegralsS04(:,:)   &
     & , IntegralsSf4(:,:), IntegralsW4(:,:)
  Complex(dp), Allocatable :: IntegralsZS04(:,:), IntegralsZSf8(:,:)     &
      & , IntegralsS0P04(:,:), IntegralsS0Sf8(:,:), IntegralsCSf4(:,:)   &
      & , IntegralsSf8(:,:), IntegralsWSpm4(:,:), IntegralsWSf8(:,:)     &
      & , IntegralsZ8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)      &
      & , IntegralsW8(:,:), IntegralsZW8(:,:), IntegralsWSpm8(:,:)
  Logical :: check, l_nmssm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoThreeBodyDecays'
  !--------------------
  ! checking for model 
  !--------------------

  If (Present(is_NMSSM)) Then
   l_nmssm = is_NMSSM
  Else
   l_nmssm = .False.
  End If

  Allocate( IntegralsZ4(n_int_diag,8) )
  Allocate( IntegralsZ8(n_int_diag,12) )
  Allocate( IntegralsW4(n_int_diag,8) )
  Allocate( IntegralsW8(n_int_diag,12) )
  Allocate( IntegralsZW8(n_int_diag,12) )
  Allocate( IntegralsS04(n_int_off_diag,10) )
  Allocate( IntegralsSf4(n_int_off_diag,10) )
  Allocate( IntegralsZS04(n_int_off_diag,12) )
  Allocate( IntegralsZS08(n_int_off_diag,16) )
  Allocate( IntegralsWSpm4(n_int_off_diag,12) )
  Allocate( IntegralsWSpm8(n_int_off_diag,16) )
  Allocate( IntegralsZSf8(n_int_off_diag,16) )
  Allocate( IntegralsWSf8(n_int_off_diag,16) )
  Allocate( IntegralsS0P04(n_int_off_diag,12) )
  Allocate( IntegralsS0P08(n_int_off_diag,16) )
  Allocate( IntegralsS0Sf8(n_int_off_diag,16) )
  Allocate( IntegralsCSf4(n_int_off_diag,12) )
  Allocate( IntegralsSf8(n_int_off_diag,16) )

  If (n_in.Lt.0) Then
   i_start = 1
   i_end = n_n
   Do i1=1,n_n
    Chi0(i1)%gi3 = 0._dp
    Chi0(i1)%bi3 = 0._dp
   End Do

  Else If ( (n_in.Ge.1).And.(n_in.Le.n_n) ) Then 
   i_start = n_in 
   i_end = n_in
   Chi0(n_in)%gi3 = 0._dp
   Chi0(n_in)%bi3 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of n_in out of range, (n_in,n_n) = ',n_in,n_n
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mC = ChiPm%m
  mN = Chi0%m
  mSlepton = Slept%m
  mSneutrino = Sneut%m
  mDsquark = Sdown%m
  mUSquark = Sup%m
  mP0 = P0%m
  mS0 = S0%m
  mSpm = Spm%m
  mGlu = Glu%m

  mSpm2 = mSpm**2
  mUsquark2 = mUsquark**2 
  mDsquark2 = mDsquark**2 
  mSlepton2 = mSlepton**2 

  m_nu = 0._dp

  mW2a = mW**2

  If (Check_Real_States) Then
   gZ_in = 0._dp
   gW_in = 0._dp
   gS0_in = 0._dp
   gP0_in = 0._dp
   gSpm_in = 0._dp
   g_Su = 0._dp
   g_Sd = 0._dp 
   g_Sl = 0._dp 
   g_Sn =  0._dp
  Else
   gZ_in = gZ
   gW_in = gW
   gS0_in = S0%g
   gP0_in = P0%g
   gSpm_in = Spm%g
   g_Su = Sup%g
   g_Sd = Sdown%g
   g_Sl = Slept%g
   g_Sn = Sneut%g
  End If
  check = Check_Real_States

  Do i_run = i_start, i_end
   !-----------------
   ! intitialisation
   !-----------------
   IntegralsZ4 = 0._dp
   IntegralsW4 = 0._dp
   IntegralsS04 = 0._dp
   IntegralsSf4 = 0._dp
   IntegralsZS04 = 0._dp
   IntegralsWSpm4 = 0._dp
   IntegralsZ8 = 0._dp
   IntegralsZS08 = 0._dp
   IntegralsZSf8 = 0._dp
   IntegralsWSf8 = 0._dp
   IntegralsS0P04 = 0._dp
   IntegralsS0P08 = 0._dp
   IntegralsS0Sf8 = 0._dp
   IntegralsCSf4 = 0._dp
   IntegralsSf8 = 0._dp

   factor(1) = oo512pi3 / Abs(mN(i_run))**3   ! for 3-body decays
   factor(2) = 3._dp * factor(1)            ! including color factor
   factor(3) = 4._dp * factor(1)            ! for decays into gluinos
   Chi0(i_run)%gi3 = 0._dp
   i_c = 1
   !--------------------------------------
   ! decays into a neutralino + 2 fermions
   !--------------------------------------
   Do i1 = 1, i_run - 1
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    If (Abs(mN(i_run)).Gt.(Abs(mN(i1))+M_D)) Then
     Call  Chi0ToChi0ff(i_run, i1, ' u u ', mN, mZ, gZ_in, Cpl_NNZ_L,Cpl_NNZ_R &
     & , n_u, mf_u, L_u, R_u, IntegralsZ4, n_Z4, mS0, gS0_in, cpl_NNS0_L       &
     & , cpl_NNS0_R, cpl_UUS0_L, cpl_UUS0_R, IntegralsS04, n_S04, mP0, gP0_in  &
     & , cpl_NNP0_L, cpl_NNP0_R, cpl_UUP0_L, cpl_UUP0_R, mUSquark, g_Su        &
     & , cpl_UNSu_L, cpl_UNSu_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04    &
     & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8        &
     & , n_S0Sf8, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI     &
     & , GenerationMixing, check, factor(2), gNff,0,ErrCan)

     Do i2=1,n_u
      Do i3=i2,n_u
       If (i2.Eq.i3) Then
        Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = id_u(i2)
        Chi0(i_run)%id3(i_c,3) = id_u(i2) + 1
        i_c = i_c +1
       Else
        Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = id_u(i2)
        Chi0(i_run)%id3(i_c,3) = id_u(i3) + 1
        Chi0(i_run)%gi3(i_c+1) = gNff(i3,i2)
        Chi0(i_run)%id3(i_c+1,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c+1,2) = id_u(i3)
        Chi0(i_run)%id3(i_c+1,3) = id_u(i2) + 1
        i_c = i_c +2
       End If
      End Do
     End Do

     Call  Chi0ToChi0ff(i_run, i1, ' d d ', mN, mZ, gZ_in, Cpl_NNZ_L,Cpl_NNZ_R &
     & , n_d, mf_d, L_d, R_d, IntegralsZ4, n_Z4, mS0, gS0_in, cpl_NNS0_L       &
     & , cpl_NNS0_R, cpl_DDS0_L, cpl_DDS0_R, IntegralsS04, n_S04, mP0, gP0_in  &
     & , cpl_NNP0_L, cpl_NNP0_R, cpl_DDP0_L, cpl_DDP0_R, mDSquark, g_Sd        &
     & , cpl_DNSd_L, cpl_DNSd_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04    &
     & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8        &
     & , n_S0Sf8, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI     &
     & , GenerationMixing, check, factor(2), gNff,0,ErrCan)

     Do i2=1,n_d
      Do i3=i2,n_d
       If (i2.Eq.i3) Then
        Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = id_d(i2)
        Chi0(i_run)%id3(i_c,3) = id_d(i2) + 1
        i_c = i_c +1
       Else
        Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = id_d(i2)
        Chi0(i_run)%id3(i_c,3) = id_d(i3) + 1
        Chi0(i_run)%gi3(i_c+1) = gNff(i3,i2)
        Chi0(i_run)%id3(i_c+1,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c+1,2) = id_d(i3)
        Chi0(i_run)%id3(i_c+1,3) = id_d(i2) + 1
        i_c = i_c +2
       End If
      End Do
     End Do


     If (n_l.Ge.1) Then
      Call  Chi0ToChi0ff(i_run, i1, ' l l ', mN, mZ,gZ_in, Cpl_NNZ_L,Cpl_NNZ_R &
      & , n_l, mf_l, L_e, R_e, IntegralsZ4, n_Z4, mS0, gS0_in, cpl_NNS0_L      &
      & , cpl_NNS0_R, cpl_LLS0_L, cpl_LLS0_R, IntegralsS04, n_S04, mP0, gP0_in &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_LLP0_L, cpl_LLP0_R, mSlepton, g_Sl       &
      & , cpl_LNSl_L, cpl_LNSl_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04   &
      & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8       &
      & , n_S0Sf8, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI    &
      & , GenerationMixing, check, factor(1), gNff,0,ErrCan)

      Do i2=1,n_l
       Do i3=i2,n_l
        If (i2.Eq.i3) Then
         Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = id_l(i2)
         Chi0(i_run)%id3(i_c,3) = id_l(i2) + 1
         i_c = i_c +1
        Else
         Chi0(i_run)%gi3(i_c) = gNff(i2,i3)
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = id_l(i2)
         Chi0(i_run)%id3(i_c,3) = id_l(i3) + 1
         Chi0(i_run)%gi3(i_c+1) = gNff(i3,i2)
         Chi0(i_run)%id3(i_c+1,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c+1,2) = id_l(i3)
         Chi0(i_run)%id3(i_c+1,3) = id_l(i2) + 1
         i_c = i_c +2
        End If
       End Do
      End Do
     End If


     If (n_nu.Ge.1) Then
      Call  Chi0ToChi0NuNu(i_run, i1, mN, mZ, gZ_in, Cpl_NNZ_L, Cpl_NNZ_R      &
        & , n_nu, L_nu, R_nu, IntegralsZ4, n_Z4, mSneutrino, g_Sn, cpl_NuNSn_L &
        & , cpl_NuNSn_R, IntegralsSf4, n_Sf4, IntegralsZSf8, n_ZSf8            &
        & , IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8                         &
        & , deltaM, epsI, GenerationMixing, check, factor(1), gNff,0,ErrCan)
      Chi0(i_run)%gi3(i_c) = Sum(gNff)
      Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
      Chi0(i_run)%id3(i_c,2) = id_nu(1)
      Chi0(i_run)%id3(i_c,3) = id_nu(1) + 1
      i_c = i_c +1
     End If
    End If
   End Do

   !--------------------------------------
   ! decay into charginos + 2 SM fermions
   !--------------------------------------
   Do i1=1,n_c
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    If (Abs(mN(i_run)).Gt.(Abs(mC(i1))+M_D)) Then
     Call Chi0ToChimffp(i_run, i1, mN, mC, 3, mf_d, mf_u, mW, gW_in, Cpl_CNW_L &
          & , Cpl_CNW_R, Cpl_UDW, mSpm, gSpm_in, Cpl_SmpCN_L, Cpl_SmpCN_R      &
          & , Cpl_SmpDU_L, Cpl_SmpDU_R, mDSquark, g_Sd, cpl_DNSd_L, cpl_DNSD_R &
          & , cpl_CUSd_L, cpl_CUSd_R, mUSquark, g_Su, cpl_UNSu_L, cpl_UNSu_R   &
          & , cpl_CDSu_L, cpl_CDSu_R, IntegralsW4, n_W4, IntegralsS04, n_S04   &
          & , IntegralsSf4, n_Sf4, IntegralsWSpm4, n_WSpm, IntegralsWSf8       &
          & , n_WSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8, n_S0Sf8         &
          & , IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI         &
          & , GenerationMixing, check, factor(2), gCffp, 0, ErrCan)

     Do i2=1,n_d
      Do i3=1,n_u
       Chi0(i_run)%gi3(i_c) = gCffp(i2,i3)
       Chi0(i_run)%id3(i_c,1) = ChiPm(i1)%id
       Chi0(i_run)%id3(i_c,2) = id_d(i2)
       Chi0(i_run)%id3(i_c,3) = id_u(i3) + 1
       i_c = i_c +1
       Chi0(i_run)%gi3(i_c) = gCffp(i2,i3)
       Chi0(i_run)%id3(i_c,1) = ChiPm(i1)%id + 1
       Chi0(i_run)%id3(i_c,2) = id_d(i2) + 1
       Chi0(i_run)%id3(i_c,3) = id_u(i3)
       i_c = i_c +1
      End Do
     End Do

     If (n_l.ge.1) Then
      Call Chi0ToChimffp(i_run, i1, mN, mC, n_l, mf_l, m_nu, mW, gW_in          &
         & , Cpl_CNW_L, Cpl_CNW_R, Cpl_NuLW, mSpm, gSpm_in, Cpl_SmpCN_L         &
         & , Cpl_SmpCN_R, Cpl_SmpLNu_L, Cpl_SmpLNu_R, mSlepton, g_Sl            &
         & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, mSneutrino, g_Sn &
         & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, IntegralsW4, n_W4&
         & , IntegralsS04, n_S04, IntegralsSf4, n_Sf4, IntegralsWSpm4, n_WSpm   &
         & , IntegralsWSf8, n_WSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8     &
         & , n_S0Sf8, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI  &
         & , GenerationMixing, check, factor(1), gCffp,0,ErrCan)
      Do i2=1,n_l
       Do i3=1,n_nu
        Chi0(i_run)%gi3(i_c) = gCffp(i2,i3)
        Chi0(i_run)%id3(i_c,1) = ChiPm(i1)%id
        Chi0(i_run)%id3(i_c,2) = id_l(i2)
        Chi0(i_run)%id3(i_c,3) = id_nu(i3) + 1
        i_c = i_c +1
        Chi0(i_run)%gi3(i_c) = gCffp(i2,i3)
        Chi0(i_run)%id3(i_c,1) = ChiPm(i1)%id + 1
        Chi0(i_run)%id3(i_c,2) = id_l(i2) + 1
        Chi0(i_run)%id3(i_c,3) = id_nu(i3)
        i_c = i_c +1
       End Do
      End Do
     End If
    End If
   End Do
   !-------------------------------
   ! decay into 3 neutralinos
   !-------------------------------  
   If (.Not.OnlySM) Then
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    Do i1=1,i_run-1
     Do i2=1,i1
      Do i3=1,i2
       diffM = Abs(mN(i_run)) - Abs(mN(i1)) - Abs(mN(i2)) - Abs(mN(i3)) - M_D
       If (diffM.Gt.0._dp) Then
        Call Chi0To3Chi0(i_run, i1, i2, i3, mN, mZ, gZ_in, cpl_NNZ_L, cpl_NNZ_R &
             & , mS0, gS0_in, cpl_NNS0_L, cpl_NNS0_R, mP0, gP0_in, cpl_NNP0_L   &
             & , cpl_NNP0_R, IntegralsZ4, n_Z4, IntegralsS04, n_S04             &
             & , IntegralsZ8, n_Z8, IntegralsZS04, n_ZS04, IntegralsZS08        &
             & , n_ZS08, IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08       &
             & , deltaM, epsI, check, factor(1), gNNN,0,ErrCan)
        Chi0(i_run)%gi3(i_c) = gNNN
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = Chi0(i2)%id
        Chi0(i_run)%id3(i_c,3) = Chi0(i3)%id
        i_c = i_c +1
       End If
      End Do
     End Do
    End Do
    !-------------------------------------
    ! decay into neutralino + 2 charginos
    !-------------------------------------
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    Do i1=1,i_run-1
     Do i2=1,n_c
      Do i3=1,i2
       diffM = Abs(mN(i_run)) - Abs(mN(i1)) - Abs(mC(i2)) - Abs(mC(i3)) - M_D
       If (diffM.Gt.0._dp) Then
        Call Chi0ToChi0ChimChip(i_run, i1, i2, i3, mN, mC, mZ, gZ_in, cpl_NNZ_L &
           & , cpl_NNZ_R, cpl_CCZ_L, cpl_CCZ_R, mW, gW_in, cpl_CNW_L, cpl_CNW_R &
           & , mSpm, gSpm_in, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0_in, cpl_NNS0_L &
           & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, mP0, gP0_in, cpl_NNP0_L      &
           & , cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, IntegralsZ4, n_Z4            &
           & , IntegralsW4, n_W4, IntegralsS04, n_S04, IntegralsZW8, n_ZW8      &
           & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsW8, n_W8  &
           & , IntegralsWSpm4, n_WSpm, IntegralsWSpm8, n_WSpm8, IntegralsS0P04  &
           & , n_S0P04, IntegralsS0P08, n_S0P08                                 &
           & , deltaM, epsI, check, factor(1), gNcc,0,ErrCan)
        If (i2.Eq.i3) Then
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
        Else
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i3)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i2)%id + 1
         i_c = i_c +1
        End If
       End If
      End Do
     End Do
    End Do
   Else If (OnlySM.And.(l_nmssm)) Then ! do nothing, is necessary because
                                       ! as also the 1-gen. R-parity model
                                       ! contains 5 neutralinos

   Else If (OnlySM.And.(n_n.Gt.4).And.(i_run.Gt.3)) Then
    !-------------------------------
    ! decay into 2 or 3 neutralinos
    !-------------------------------  
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    Do i1=1,i_run-1
     Do i2=1,Min(i1,3)
      Do i3=1,i2
       diffM = Abs(mN(i_run)) - Abs(mN(i1)) - Abs(mN(i2)) - Abs(mN(i3)) - M_D
       If (diffM.Gt.0._dp) Then
        Call Chi0To3Chi0(i_run, i1, i2, i3, mN, mZ, gZ_in, cpl_NNZ_L, cpl_NNZ_R &
           & , mS0, gS0_in, cpl_NNS0_L, cpl_NNS0_R, mP0, gP0_in, cpl_NNP0_L     &
           & , cpl_NNP0_R, IntegralsZ4, n_Z4, IntegralsS04, n_S04, IntegralsZ8  &
           & , n_Z8, IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08               &
           & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08                 &
           & , deltaM, epsI, check, factor(1), gNNN,0,ErrCan)
        Chi0(i_run)%gi3(i_c) = gNNN
        Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
        Chi0(i_run)%id3(i_c,2) = Chi0(i2)%id
        Chi0(i_run)%id3(i_c,3) = Chi0(i3)%id
        i_c = i_c +1
       End If
      End Do
     End Do
    End Do
    !-------------------------------------
    ! decay into neutralino + 2 charginos
    !-------------------------------------
    n_Z4 = 0
    n_W4 = 0
    n_S04 = 0
    n_Sf4 = 0
    n_ZS04 = 0
    n_WSpm = 0
    n_WSpm8 = 0
    n_Z8 = 0
    n_W8 = 0
    n_ZW8 = 0
    n_ZSf8 = 0
    n_ZS08 = 0
    n_WSf8 = 0
    n_S0P04 = 0
    n_S0P08 = 0
    n_S0Sf8 = 0
    n_CSf4 = 0
    n_Sf8 = 0
    Do i1=1,3
     Do i2=1,n_c
      Do i3=1,Min(i2,3)
       diffM = Abs(mN(i_run)) - Abs(mN(i1)) - Abs(mC(i2)) - Abs(mC(i3)) - M_D
       If (diffM.Gt.0._dp) Then
        Call Chi0ToChi0ChimChip(i_run, i1, i2, i3, mN, mC, mZ, gZ_in, cpl_NNZ_L &
           & , cpl_NNZ_R, cpl_CCZ_L, cpl_CCZ_R, mW, gW_in, cpl_CNW_L, cpl_CNW_R &
           & , mSpm, gSpm_in, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0_in, cpl_NNS0_L &
           & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, mP0, gP0_in, cpl_NNP0_L      &
           & , cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, IntegralsZ4, n_Z4            &
           & , IntegralsW4, n_W4, IntegralsS04, n_S04, IntegralsZW8, n_ZW8      &
           & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsW8        &
           & , n_W8, IntegralsWSpm4, n_WSpm, IntegralsWSpm8, n_WSpm8            &
           & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08                 &
           & , deltaM, epsI, check, factor(1), gNcc,0,ErrCan)

        If (i2.Eq.i3) Then
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
        Else
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i3)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i2)%id + 1
         i_c = i_c +1
        End If
       End If
      End Do
     End Do
    End Do
    Do i1=4,i_run-1
     Do i2=1,3
      Do i3=1,i2
       diffM = Abs(mN(i_run)) - Abs(mN(i1)) - Abs(mC(i2)) - Abs(mC(i3)) - M_D
       If (diffM.Gt.0._dp) Then
        Call Chi0ToChi0ChimChip(i_run, i1, i2, i3, mN, mC, mZ, gZ_in, cpl_NNZ_L &
           & , cpl_NNZ_R, cpl_CCZ_L, cpl_CCZ_R, mW, gW_in, cpl_CNW_L, cpl_CNW_R &
           & , mSpm, gSpm_in, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0_in, cpl_NNS0_L &
           & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, mP0, gP0_in, cpl_NNP0_L      &
           & , cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, IntegralsZ4, n_Z4            &
           & , IntegralsW4, n_W4, IntegralsS04, n_S04, IntegralsZW8, n_ZW8      &
           & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsW8, n_W8  &
           & , IntegralsWSpm4, n_WSpm, IntegralsWSpm8, n_WSpm8                  &
           & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08                 &
           & , deltaM, epsI, check, factor(1), gNcc,0,ErrCan)

        If (i2.Eq.i3) Then
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
        Else
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i2)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i3)%id + 1
         i_c = i_c +1
         Chi0(i_run)%gi3(i_c) = gNCC
         Chi0(i_run)%id3(i_c,1) = Chi0(i1)%id
         Chi0(i_run)%id3(i_c,2) = ChiPm(i3)%id
         Chi0(i_run)%id3(i_c,3) = ChiPm(i2)%id + 1
         i_c = i_c +1
        End If
       End If
      End Do
     End Do
    End Do
   End If
   !---------------------------------
   ! decays into a gluino + 2 quarks
   !---------------------------------
   n_Sf4 = 0
   n_CSf4 = 0
   n_Sf8 = 0
   If ( Abs(mN(i_run)).Gt.(mGlu + M_D) ) Then
    Call Chi0toGqq(i_run, ' u u ', mN, mGlu, mf_u, mUSquark, g_Su       &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_UGSu_L, cpl_UGSu_R                  &
    & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8 &
    & , deltaM, epsI, GenerationMixing, check, factor(3), gGqq,0,ErrCan)
    Do i2=1,n_u
     Do i3=i2,n_u
      If (i2.Eq.i3) Then
       Chi0(i_run)%gi3(i_c) = gGqq(i2,i3)
       Chi0(i_run)%id3(i_c,1) = Glu%id
       Chi0(i_run)%id3(i_c,2) = id_u(i2)
       Chi0(i_run)%id3(i_c,3) = id_u(i2) + 1
       i_c = i_c +1
      Else
       Chi0(i_run)%gi3(i_c) = gGqq(i2,i3)
       Chi0(i_run)%id3(i_c,1) = Glu%id
       Chi0(i_run)%id3(i_c,2) = id_u(i2)
       Chi0(i_run)%id3(i_c,3) = id_u(i3) + 1
       Chi0(i_run)%gi3(i_c+1) = gGqq(i3,i2)
       Chi0(i_run)%id3(i_c+1,1) = Glu%id
       Chi0(i_run)%id3(i_c+1,2) = id_u(i3)
       Chi0(i_run)%id3(i_c+1,3) = id_u(i2) + 1
       i_c = i_c +2
      End If
     End Do
    End Do

    Call Chi0toGqq(i_run, ' d d ', mN, mGlu, mf_d, mDSquark, g_Sd       &
    & , cpl_DNSd_L, cpl_DNSd_R, cpl_DGSd_L, cpl_DGSd_R                  &
    & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8 &
    & , deltaM, epsI, GenerationMixing, check, factor(3), gGqq ,0,ErrCan)

    Do i2=1,n_u
     Do i3=i2,n_u
      If (i2.Eq.i3) Then
       Chi0(i_run)%gi3(i_c) = gGqq(i2,i3)
       Chi0(i_run)%id3(i_c,1) = Glu%id
       Chi0(i_run)%id3(i_c,2) = id_d(i2)
       Chi0(i_run)%id3(i_c,3) = id_d(i2) + 1
       i_c = i_c +1
      Else
       Chi0(i_run)%gi3(i_c) = gGqq(i2,i3)
       Chi0(i_run)%id3(i_c,1) = Glu%id
       Chi0(i_run)%id3(i_c,2) = id_d(i2)
       Chi0(i_run)%id3(i_c,3) = id_d(i3) + 1
       Chi0(i_run)%gi3(i_c+1) = gGqq(i3,i2)
       Chi0(i_run)%id3(i_c+1,1) = Glu%id
       Chi0(i_run)%id3(i_c+1,2) = id_d(i3)
       Chi0(i_run)%id3(i_c+1,3) = id_d(i2) + 1
       i_c = i_c +2
      End If
     End Do
    End Do
   End If
   !-------------------------------
   ! decay into neutralino + photon
   !-------------------------------  
   factor(1) = - gSU2 * Sqrt(sW2) * oo8pi2
   factor(2) = 0.25_dp * factor(1)
   factor(3) = 0.125_dp / (Pi * Abs(mN(i_run))**3 ) 

   Do i1=1,i_run-1
    Call  Chi0ToChi0Photon(i_run, i1, mN, mW2a, mC, Cpl_CNW_L, Cpl_CNW_R     &
       & , mSpm2, Cpl_SmpCN_L, Cpl_SmpCN_R, mf_u, mUsquark2, cpl_UNSu_L      &
       & , cpl_UNSu_R, mf_d, mDsquark2, cpl_DNSd_L, cpl_DNSd_R, mf_l         &
       & , mSlepton2, cpl_LNSl_L, cpl_LNSl_R, factor, gNNN )

    Chi0(i_run)%gi2(200-i_run+i1) = gNNN
    Chi0(i_run)%id2(200-i_run+i1,1) = Chi0(i1)%id
    Chi0(i_run)%id2(200-i_run+i1,2) = id_ph
   End Do

   Chi0(i_run)%g = Sum(Chi0(i_run)%gi2) + Sum(Chi0(i_run)%gi3)
   If (Chi0(i_run)%g.Ne.0._dp) Then
    Chi0(i_run)%bi2 = Chi0(i_run)%gi2 / Chi0(i_run)%g
    Chi0(i_run)%bi3 = Chi0(i_run)%gi3 / Chi0(i_run)%g
   End If

  End Do ! i_run

  Deallocate( IntegralsZ4, IntegralsZ8, IntegralsW4, IntegralsW8    &
    & , IntegralsZW8, IntegralsS04, IntegralsSf4, IntegralsZS04            &
    & , IntegralsZS08, IntegralsWSpm4, IntegralsWSpm8, IntegralsZSf8       &
    & , IntegralsWSf8, IntegralsS0P04, IntegralsS0P08, IntegralsS0Sf8      &
    & , IntegralsCSf4, IntegralsSf8 )

  Iname = Iname - 1

 End Subroutine NeutralinoThreeBodyDecays


 Subroutine Chi0To3Chi0(i_in, i1, i2, i3, mN, mZ, gZ, cpl_NNZ_L, cpl_NNZ_R    &
    & , mS0, gS0, cpl_NNS0_L, cpl_NNS0_R, mP0, gP0, cpl_NNP0_L, cpl_NNP0_R    &
    & , IntegralsZ4, n_Z4, IntegralsS04, n_S04, IntegralsZ8, n_Z8             &
    & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsS0P04, n_S0P04 &
    & , IntegralsS0P08, n_S0P08                                               &
    & , deltaM, epsI, check, fac, gNNN, WriteContribution, n_out)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a neutralino into 3 neutralinos
 ! input:
 !    i_in ....... index of the decaying neutralino
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !
 ! output:
 ! written by Werner Porod, 26.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, i1, i2, i3
  Integer, Intent(inout) :: n_Z4, n_S04, n_ZS04, n_S0P04, n_Z8, n_ZS08, n_S0P08
  Logical, Intent(in) :: check

  Real(dp), Intent(in) :: mN(:), mZ, gZ, mS0(:), gS0(:), mP0(:), gP0(:) &
     & , epsI, deltaM, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsS04(:,:)

  Complex(dp), Intent(in) :: Cpl_NNZ_L(:,:), Cpl_NNZ_R(:,:)             &
      & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_NNP0_L(:,:,:)       &
      & , cpl_NNP0_R(:,:,:)

  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsS0P04(:,:) &
      & , IntegralsZ8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)

  Real(dp), Intent(out) :: gNNN

  Integer, Intent(in)  :: WriteContribution, n_out

  Integer :: Isum, n1, n2, n_S0, n_P0
  Real(dp) :: Boson2(2), Boson4(4), mass(4), resR, resRa, resRb
  Complex(dp) :: coup1(4), coup2(8), resC, resCa, resCb, resCc, resCd, resCe

  Real(dp), Allocatable :: gNNNSum(:)
  Character(len=15), Allocatable :: Contribution(:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0To3Chi0'

  mass(1) = mN(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )

  Isum = 9 * (1 + n_S0 + n_P0)**2
  Allocate( gNNNSum(Isum) )
  Allocate( Contribution(Isum) )
   
  gNNNSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ

   !--------------
   ! s-channel
   !--------------
   coup1(1) = Cpl_NNZ_L(i1,i_in)
   coup1(2) = Cpl_NNZ_R(i1,i_in)
   coup1(3) = Cpl_NNZ_L(i2,i3)
   coup1(4) = Cpl_NNZ_R(i2,i3)
   mass(2) = mN(i1)
   mass(3) = -mN(i2)
   mass(4) = mN(i3)       
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check)
   !--------------
   ! t-channel
   !--------------
   coup1(1) = Cpl_NNZ_L(i2,i_in)
   coup1(2) = Cpl_NNZ_R(i2,i_in)
   coup1(3) = Cpl_NNZ_L(i1,i3)
   coup1(4) = Cpl_NNZ_R(i1,i3)
   mass(2) = mN(i2)
   mass(3) = -mN(i1)
   mass(4) = mN(i3)       
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resRa, check)
   !--------------
   ! u-channel
   !--------------
   coup1(1) = Cpl_NNZ_L(i3,i_in)
   coup1(2) = Cpl_NNZ_R(i3,i_in)
   coup1(3) = Cpl_NNZ_L(i1,i2)
   coup1(4) = Cpl_NNZ_R(i1,i2)
   mass(2) = mN(i3)
   mass(3) = -mN(i1)
   mass(4) = mN(i2)       
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resRb, check)
   gNNNSum(Isum) = resR + resRa + resRb
   Contribution(Isum) = 'Z'

  !-------------------
  ! Scalar S_0
  !-------------------
  Do n1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(n1)
   Boson2(2) = gS0(n1)

   !--------------
   ! s-channel
   !--------------
   mass(2) = mN(i1)
   mass(3) = -mN(i2)
   mass(4) = mN(i3)
   coup1(1) = cpl_NNS0_L(i1,i_in,n1)
   coup1(2) = cpl_NNS0_R(i1,i_in,n1)
   coup1(3) = cpl_NNS0_L(i2,i3,n1)
   coup1(4) = cpl_NNS0_R(i2,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   !--------------
   ! t-channel
   !--------------
   mass(2) = mN(i2)
   mass(3) = -mN(i1)
   mass(4) = mN(i3)
   coup1(1) = cpl_NNS0_L(i2,i_in,n1)
   coup1(2) = cpl_NNS0_R(i2,i_in,n1)
   coup1(3) = cpl_NNS0_L(i1,i3,n1)
   coup1(4) = cpl_NNS0_R(i1,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resRa, check)
   !--------------
   ! u-channel
   !--------------
   mass(2) = mN(i3)
   mass(3) = -mN(i1)
   mass(4) = mN(i2)
   coup1(1) = cpl_NNS0_L(i3,i_in,n1)
   coup1(2) = cpl_NNS0_R(i3,i_in,n1)
   coup1(3) = cpl_NNS0_L(i1,i2,n1)
   coup1(4) = cpl_NNS0_R(i1,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, resRb, check)
   gNNNSum(Isum) = resR+resRa+resRb
   Contribution(Isum) = 'S0_'//Bu(n1)
  End Do ! n1 ,  S^0_n1

  !-------------------
  ! Pseudoscalar P_0
  !-------------------
  Do n1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(n1)
   Boson2(2) = gP0(n1)
    
   !--------------
   ! s-channel
   !--------------
   mass(2) = mN(i1)
   mass(3) = -mN(i2)
   mass(4) = mN(i3)

   coup1(1) = cpl_NNP0_L(i1,i_in,n1)
   coup1(2) = cpl_NNP0_R(i1,i_in,n1)
   coup1(3) = cpl_NNP0_L(i2,i3,n1)
   coup1(4) = cpl_NNP0_R(i2,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   !--------------
   ! t-channel
   !--------------
   mass(2) = mN(i2)
   mass(3) = -mN(i1)
   mass(4) = mN(i3)

   coup1(1) = cpl_NNP0_L(i2,i_in,n1)
   coup1(2) = cpl_NNP0_R(i2,i_in,n1)
   coup1(3) = cpl_NNP0_L(i1,i3,n1)
   coup1(4) = cpl_NNP0_R(i1,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resRa, check)
   !--------------
   ! u-channel
   !--------------
   mass(2) = mN(i3)
   mass(3) = -mN(i1)
   mass(4) = mN(i2)

   coup1(1) = cpl_NNP0_L(i3,i_in,n1)
   coup1(2) = cpl_NNP0_R(i3,i_in,n1)
   coup1(3) = cpl_NNP0_L(i1,i2,n1)
   coup1(4) = cpl_NNP0_R(i1,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resRb, check)
   gNNNSum(Isum) = resR+resRa+resRb
   Contribution(Isum) = 'P0_'//Bu(n1)
  End Do ! n1 ,  P^0_n1

  !-------------------------------
  ! Z boson - Z-boson
  !-------------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ
  Boson4(3) = mZ
  Boson4(4) = gZ
  Isum = Isum + 1

  !-----------
  ! s-t channel
  !-----------
  mass(2) = mN(i2)
  mass(3) = -mN(i3)
  mass(4) = mN(i1)
  coup2(1) = Cpl_NNZ_L(i1,i_in)
  coup2(2) = Cpl_NNZ_R(i1,i_in)
  coup2(3) = Conjg(Cpl_NNZ_L(i2,i_in))
  coup2(4) = Conjg(Cpl_NNZ_R(i2,i_in))
  coup2(5) = Cpl_NNZ_R(i3,i2)
  coup2(6) = Cpl_NNZ_L(i3,i2)
  coup2(7) = Conjg(Cpl_NNZ_R(i3,i1))
  coup2(8) = Conjg(Cpl_NNZ_L(i3,i1))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZ8, n_Z8, resC, check)
  !-----------
  ! s-u channel
  !-----------
  mass(2) = mN(i3)
  mass(3) = -mN(i2)
  mass(4) = mN(i1)
  coup2(1) = Cpl_NNZ_L(i1,i_in)
  coup2(2) = Cpl_NNZ_R(i1,i_in)
  coup2(3) = Conjg(Cpl_NNZ_L(i3,i_in))
  coup2(4) = Conjg(Cpl_NNZ_R(i3,i_in))
  coup2(5) = Cpl_NNZ_R(i2,i3)
  coup2(6) = Cpl_NNZ_L(i2,i3)
  coup2(7) = Conjg(Cpl_NNZ_R(i2,i1))
  coup2(8) = Conjg(Cpl_NNZ_L(i2,i1))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZ8, n_Z8, resCa, check)
  !-----------
  ! t-u channel
  !-----------
  mass(2) = mN(i3)
  mass(3) = -mN(i1)
  mass(4) = mN(i2)
  coup2(1) = Cpl_NNZ_L(i2,i_in)
  coup2(2) = Cpl_NNZ_R(i2,i_in)
  coup2(3) = Conjg(Cpl_NNZ_L(i3,i_in))
  coup2(4) = Conjg(Cpl_NNZ_R(i3,i_in))
  coup2(5) = Cpl_NNZ_R(i1,i3)
  coup2(6) = Cpl_NNZ_L(i1,i3)
  coup2(7) = Conjg(Cpl_NNZ_R(i1,i2))
  coup2(8) = Conjg(Cpl_NNZ_L(i1,i2))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZ8, n_Z8, resCb, check)
  gNNNSum(Isum) = -2._dp * Real(resC+resCa+resCb,dp)
  Contribution(Isum) = 'Z Z s-t'

  !-------------------------------
  ! Z boson - S^0 boson  s-channel
  !-------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_S0 
    Isum = Isum + 1
    Boson4(3) = mS0(n1)
    Boson4(4) = gS0(n1)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i2,i3)
    coup2(6) = Cpl_NNZ_R(i2,i3)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i1,i3)
    coup2(6) = Cpl_NNZ_R(i1,i3)
    coup2(7) = Conjg(cpl_NNS0_R(i1,i3,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i1,i3,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCa, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i1,i2)
    coup2(6) = Cpl_NNZ_R(i1,i2)
    coup2(7) = Conjg(cpl_NNS0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCb, check)
    gNNNSum(Isum) = 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'Z S0_'//Bu(n1)//' s-s'
   End Do  ! n1 n_S0

   !----------------------------------
   ! Z boson - S^0 boson  s-t channel
   !----------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_S0 
    Isum = Isum + 1
    Boson4(3) = mS0(n1)
    Boson4(4) = gS0(n1)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mN(i1)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i3,i2)
    coup2(6) = Cpl_NNZ_L(i3,i2)
    coup2(7) = Conjg(cpl_NNS0_R(i3,i1,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i3,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mN(i1)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i3,i2)
    coup2(6) = Cpl_NNZ_R(i3,i2)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i1,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCa, check)
    !-------------
    ! t-s channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mN(i3)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i3,i1)
    coup2(6) = Cpl_NNZ_L(i3,i1)
    coup2(7) = Conjg(cpl_NNS0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCb, check)
    !-------------
    ! t-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i3,i1)
    coup2(6) = Cpl_NNZ_R(i3,i1)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i1,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCc, check)
    !-------------
    ! u-s channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i2,i1)
    coup2(6) = Cpl_NNZ_L(i2,i1)
    coup2(7) = Conjg(cpl_NNS0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCd, check)
    !-------------
    ! u-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNS0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNS0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i2,i1)
    coup2(6) = Cpl_NNZ_R(i2,i1)
    coup2(7) = Conjg(cpl_NNS0_R(i3,i1,n1))
    coup2(8) = Conjg(cpl_NNS0_L(i3,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCe, check)
    gNNNSum(Isum) = 2._dp * Real(resC-resCa+resCb-resCc+resCd-resCe,dp)
    Contribution(Isum) = 'Z S0_'//Bu(n1)//' s-t'
   End Do  ! n1 n_S0

  !-------------------------------
  ! Z boson - P^0 boson  s-channel
  !-------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_P0 
    Isum = Isum + 1
    Boson4(3) = mP0(n1)
    Boson4(4) = gP0(n1)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i2,i3)
    coup2(6) = Cpl_NNZ_R(i2,i3)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i1,i3)
    coup2(6) = Cpl_NNZ_R(i1,i3)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i3,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i3,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCa, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i1,i2)
    coup2(6) = Cpl_NNZ_R(i1,i2)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCb, check)
    gNNNSum(Isum) = 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'Z P0_'//Bu(n1)//' s-s'
   End Do  ! n1 n_P0

   !----------------------------------
   ! Z boson - P^0 boson  s-t channel
   !----------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_P0 
    Isum = Isum + 1
    Boson4(3) = mP0(n1)
    Boson4(4) = gP0(n1)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mN(i1)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i3,i2)
    coup2(6) = Cpl_NNZ_L(i3,i2)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i1,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mN(i1)
    coup2(1) = Cpl_NNZ_L(i1,i_in)
    coup2(2) = Cpl_NNZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i3,i2)
    coup2(6) = Cpl_NNZ_R(i3,i2)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCa, check)
    !-------------
    ! t-s channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mN(i3)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i3,i1)
    coup2(6) = Cpl_NNZ_L(i3,i1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCb, check)
    !-------------
    ! t-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = Cpl_NNZ_L(i2,i_in)
    coup2(2) = Cpl_NNZ_R(i2,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i3,i1)
    coup2(6) = Cpl_NNZ_R(i3,i1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCc, check)
    !-------------
    ! u-s channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n1))
    coup2(5) = Cpl_NNZ_R(i2,i1)
    coup2(6) = Cpl_NNZ_L(i2,i1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCd, check)
    !-------------
    ! u-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = Cpl_NNZ_L(i3,i_in)
    coup2(2) = Cpl_NNZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n1))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n1))
    coup2(5) = Cpl_NNZ_L(i2,i1)
    coup2(6) = Cpl_NNZ_R(i2,i1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i1,n1))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i1,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCe, check)
    gNNNSum(Isum) = 2._dp * Real(resC-resCa+resCb-resCc+resCd-resCe,dp)
    Contribution(Isum) = 'Z P0_'//Bu(n1)//' s-t'
   End Do  ! n1 n_P0

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = n1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n2))
    coup2(5) = cpl_NNS0_L(i2,i3,n1)
    coup2(6) = cpl_NNS0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i2,i_in,n1)
    coup2(2) = cpl_NNS0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i2,i_in,n2))
    coup2(5) = cpl_NNS0_L(i1,i3,n1)
    coup2(6) = cpl_NNS0_R(i1,i3,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i1,i3,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i1,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNS0_L(i3,i_in,n1)
    coup2(2) = cpl_NNS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i1,i2,n1)
    coup2(6) = cpl_NNS0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCb, check)
    gNNNSum(Isum) = 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' S0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i2,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i2,n1)
    coup2(6) = cpl_NNS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i3,i1,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i3,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i2,n1)
    coup2(6) = cpl_NNS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    !-----------
    ! t-u channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNS0_L(i2,i_in,n1)
    coup2(2) = cpl_NNS0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i1,n1)
    coup2(6) = cpl_NNS0_R(i3,i1,n1)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCb, check)
    gNNNSum(Isum) = - 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' S0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_NNS0_L(i2,i3,n1)
    coup2(6) = cpl_NNS0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i2,i_in,n1)
    coup2(2) = cpl_NNS0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n2))
    coup2(5) = cpl_NNS0_L(i1,i3,n1)
    coup2(6) = cpl_NNS0_R(i1,i3,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNS0_L(i3,i_in,n1)
    coup2(2) = cpl_NNS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i1,i2,n1)
    coup2(6) = cpl_NNS0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCb, check)
    gNNNSum(Isum) = 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' P0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i2,n1)
    coup2(6) = cpl_NNS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i2,n1)
    coup2(6) = cpl_NNS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    !-----------
    ! t-s channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i3)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNS0_L(i2,i_in,n1)
    coup2(2) = cpl_NNS0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i1,n1)
    coup2(6) = cpl_NNS0_R(i3,i1,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCb, check)
    !-----------
    ! t-u channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNS0_L(i2,i_in,n1)
    coup2(2) = cpl_NNS0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNS0_L(i3,i1,n1)
    coup2(6) = cpl_NNS0_R(i3,i1,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCc, check)
    !-----------
    ! u-s channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i3,i_in,n1)
    coup2(2) = cpl_NNS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_NNS0_L(i2,i1,n1)
    coup2(6) = cpl_NNS0_R(i2,i1,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCd, check)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNS0_L(i3,i_in,n1)
    coup2(2) = cpl_NNS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n2))
    coup2(5) = cpl_NNS0_L(i2,i1,n1)
    coup2(6) = cpl_NNS0_R(i2,i1,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCe, check)
    gNNNSum(Isum) = - 2._dp * Real(resC+resCa+resCb+resCc+resCd+resCe,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' P0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   Do n2 = n1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mN(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNP0_L(i1,i_in,n1)
    coup2(2) = cpl_NNP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_NNP0_L(i2,i3,n1)
    coup2(6) = cpl_NNP0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i1)
    mass(4) = mN(i3)
    coup2(1) = cpl_NNP0_L(i2,i_in,n1)
    coup2(2) = cpl_NNP0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n2))
    coup2(5) = cpl_NNP0_L(i1,i3,n1)
    coup2(6) = cpl_NNP0_R(i1,i3,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNP0_L(i3,i_in,n1)
    coup2(2) = cpl_NNP0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNP0_L(i1,i2,n1)
    coup2(6) = cpl_NNP0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCb, check)
    gNNNSum(Isum) = 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//' P0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNP0_L(i1,i_in,n1)
    coup2(2) = cpl_NNP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i2,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i2,i_in,n2))
    coup2(5) = cpl_NNP0_L(i3,i2,n1)
    coup2(6) = cpl_NNP0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i3,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i3,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mN(i1)
    coup2(1) = cpl_NNP0_L(i1,i_in,n1)
    coup2(2) = cpl_NNP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNP0_L(i3,i2,n1)
    coup2(6) = cpl_NNP0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    !-----------
    ! t-u channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i1)
    mass(4) = mN(i2)
    coup2(1) = cpl_NNP0_L(i2,i_in,n1)
    coup2(2) = cpl_NNP0_R(i2,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i3,i_in,n2))
    coup2(5) = cpl_NNP0_L(i3,i1,n1)
    coup2(6) = cpl_NNP0_R(i3,i1,n1)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i1,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i1,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCb, check)
    gNNNSum(Isum) = - 2._dp * Real(resC+resCa+resCb,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//' P0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !----------
  ! Summing
  !----------
  gNNN = sum( gNNNSum(1:Isum) )
  If (gNNN.Lt.0._dp) Then
   p_test = Abs(gNNN) / Maxval(Abs(gNNNSum(1:Isum)))
   gNNN = 0._dp
   If (p_test.gt.prec) then
    Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//Bu(i2)//Bu(i3)//') < 0 :' &
      & ,i_in,i1,i2,i3,gNNN
    Write(ErrCan,*) 'The different contributions are :'
    Do n1=1,Isum
     If (gNNNSum(n1).Ne.0._dp)  Write(ErrCan,*) Contribution(n1),gNNNSum(n1)
    End Do
   End If
  End If

  !----------------------------------------------
  ! symmetry due to majorana nature
  !----------------------------------------------
  If ((i1.Eq.i2).And.(i2.Eq.i3)) Then
   gNNN = gNNN / 6._dp
  Else If ((i1.Eq.i2).Or.(i2.Eq.i3).Or.(i1.Eq.i3)) Then
   gNNN = 0.5_dp * gNNN 
  End If

  gNNN = fac * gNNN

  !---------------------------
  ! for detailed information
  !---------------------------

  If (WriteContribution.ne.0) Then

   gNNNSum = gNNNSum * fac
 
   Write (n_out,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//Bu(i2)//Bu(i3)//') :' &
     & ,i_in,i1,i2,i3
   Write (n_out,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gNNNSum(n1).Ne.0._dp)  Write (n_out,*) Contribution(n1),gNNNSum(n1)
   End Do
  End If

  Deallocate( gNNNSum, Contribution )
   
  Iname = Iname - 1

 End Subroutine Chi0To3Chi0

 Subroutine Chi0ToChi0ChimChip(i_in, i1, i2, i3, mN, mC, mZ, gZ, cpl_NNZ_L    &
    & , cpl_NNZ_R, cpl_CCZ_L, cpl_CCZ_R, mW, gW, cpl_CNW_L, cpl_CNW_R, mSpm   &
    & , gSpm, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0, cpl_NNS0_L, cpl_NNS0_R      &
    & , cpl_CCS0_L, cpl_CCS0_R, mP0, gP0, cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L  &
    & , cpl_CCP0_R, IntegralsZ4, n_Z4, IntegralsW4, n_W4, IntegralsS04, n_S04 &
    & , IntegralsZW8, n_ZW8, IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08     &
    & , IntegralsW8, n_W8, IntegralsWSpm4, n_WSpm4, IntegralsWSpm8, n_WSpm8   &
    & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08                      &
    & , deltaM, epsI, check, fac, gNcc, WriteContribution, n_out)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a neutralino into another neutralino +
 ! two charginos
 ! input:
 !    i_in ....... index of the decaying neutralino, in case n_in < 0,
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !    mC ......... chargino masses
 !
 ! output:
 ! written by Werner Porod, 26.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, i1, i2, i3
  Integer, Intent(inout) :: n_Z4, n_W4, n_S04, n_ZS04, n_S0P04, n_ZW8, n_ZS08 &
     & , n_W8, n_S0P08, n_WSpm4, n_WSpm8
  Logical, Intent(in) :: check

  Real(dp), Intent(in) :: mN(:), mZ, gZ, mC(:), mW, gW, mSpm(:), gSpm(:)   &
     & , mS0(:), gS0(:), mP0(:), gP0(:), epsI, deltaM, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsW4(:,:)   &
     & , IntegralsS04(:,:)

  Complex(dp), Intent(in) :: cpl_NNZ_L(:,:), cpl_NNZ_R(:,:), cpl_CCZ_L(:,:)  &
     & , cpl_CCZ_R(:,:), cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_SmpCN_L(:,:,:)  &
     & , cpl_SmpCN_R(:,:,:), cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:)            &
     & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:), cpl_NNP0_L(:,:,:)             &
     & , cpl_NNP0_R(:,:,:), cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:)

  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsS0P04(:,:) &
      & , IntegralsZW8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)    &
      & , IntegralsW8(:,:), IntegralsWSpm4(:,:), IntegralsWSpm8(:,:)

  Real(dp), Intent(out) :: gNCC

  Integer, Intent(in) :: n_out, WriteContribution

  Integer :: Isum, n1, n2, n_S0, n_P0, n_Spm
  Real(dp) :: Boson2(2), Boson4(4), mass(4), resR, resRa
  Complex(dp) :: coup1(4), coup2(8), resC, resCa

  Real(dp), Allocatable :: gNCCSum(:)
  Character(len=12), Allocatable :: Contribution(:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0ToChi0ChimChip'

  mass(1) = mN(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )
  n_Spm = Size( mSpm )

  Isum = (3 + 2 * n_Spm + n_S0 + n_P0)**2
  Allocate( gNCCSum(Isum) )
  Allocate( Contribution(Isum) )
   
  gNCCSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3) 
  coup1(1) = Cpl_NNZ_L(i1,i_in)
  coup1(2) = Cpl_NNZ_R(i1,i_in)
  coup1(3) = Cpl_CCZ_L(i2,i3)
  coup1(4) = Cpl_CCZ_R(i2,i3)
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check)
  gNCCSum(Isum) = resR
  Contribution(Isum) = 'Z'

  !--------------
  ! W W
  !--------------
  Isum = Isum + 1
  Boson2(1) = mW
  Boson2(2) = gW
  !---------------
  ! t-channel
  !---------------
  coup1(1) = Conjg(Cpl_CNW_R(i2,i_in))
  coup1(2) = Conjg(Cpl_CNW_L(i2,i_in))
  coup1(3) = Cpl_CNW_L(i3,i1)
  coup1(4) = Cpl_CNW_R(i3,i1)
  mass(2) = mC(i2)
  mass(3) = -mN(i1)
  mass(4) = mC(i3)                 
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                       &, n_W4, resR, check)
  !---------------
  ! u-channel
  !---------------
  coup1(1) = Cpl_CNW_L(i3,i_in)
  coup1(2) = Cpl_CNW_R(i3,i_in)
  coup1(3) = Conjg(Cpl_CNW_L(i2,i1))
  coup1(4) = Conjg(Cpl_CNW_R(i2,i1))
  mass(2) = mC(i3)
  mass(3) = -mN(i1)
  mass(4) = mC(i2)       
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                       &, n_W4, resRa, check)
  gNCCSum(Isum) = resR+resRa
  Contribution(Isum) = 'W'

  !---------------------------------------------
  ! decays via S^\pm, without interference part 
  !---------------------------------------------
  Do n1=1,n_Spm
   Isum = Isum + 1
   Boson2(1) = mSpm(n1)
   Boson2(2) = gSpm(n1)
   !------------------
   ! t-channel
   !------------------
   mass(2) = mC(i2)
   mass(3) = -mN(i1)
   mass(4) = mC(i3)
   coup1(1) = Conjg(Cpl_SmpCN_R(n1,i2,i_in))
   coup1(2) = Conjg(Cpl_SmpCN_L(n1,i2,i_in))
   coup1(3) = Cpl_SmpCN_L(n1,i3,i1)
   coup1(4) = Cpl_SmpCN_R(n1,i3,i1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04 ,n_S04, resRa, check)
   !------------------
   ! u-channel
   !------------------
   mass(2) = mC(i3)
   mass(3) = -mN(i1)
   mass(4) = mC(i2)
   coup1(1) = Cpl_SmpCN_L(n1,i3,i_in)
   coup1(2) = Cpl_SmpCN_R(n1,i3,i_in)
   coup1(3) = Conjg(Cpl_SmpCN_R(n1,i2,i1))
   coup1(4) = Conjg(Cpl_SmpCN_L(n1,i2,i1))
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04 ,n_S04, resR, check)
   gNCCSum(Isum) = resR+resRa
   Contribution(Isum) = 'S-_'//Bu(n1)
  End Do  ! n1 S^\pm

  !-------------------
  ! Scalar S_0
  !-------------------
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  Do n1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(n1)
   Boson2(2) = gS0(n1)
   coup1(1) = cpl_NNS0_L(i1,i_in,n1)
   coup1(2) = cpl_NNS0_R(i1,i_in,n1)
   coup1(3) = cpl_CCS0_L(i3,i2,n1)
   coup1(4) = cpl_CCS0_R(i3,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   gNCCSum(Isum) = resR
   Contribution(Isum) = 'S0_'//Bu(n1)
  End Do ! n1 ,  S^0_n1

  !-------------------
  ! Pseudoscalar P_0
  !-------------------
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  Do n1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(n1)
   Boson2(2) = gP0(n1)
   coup1(1) = cpl_NNP0_L(i1,i_in,n1)
   coup1(2) = cpl_NNP0_R(i1,i_in,n1)
   coup1(3) = cpl_CCP0_L(i2,i3,n1)
   coup1(4) = cpl_CCP0_R(i2,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   gNCCSum(Isum) = resR
   Contribution(Isum) = 'P0_'//Bu(n1)
  End Do ! n1 ,  P^0_n1

  !-------
  ! Z - W
  !-------
  Boson4(1) = mW
  Boson4(2) = gW
  Boson4(3) = mZ
  Boson4(4) = gZ
  Isum = Isum + 1
  !-----------
  ! t-s channel
  !-----------
  mass(2) = mN(i1)
  mass(3) = -mC(i3)
  mass(4) = mC(i2)
  coup2(1) = Conjg(Cpl_CNW_R(i2,i_in))
  coup2(2) = Conjg(Cpl_CNW_L(i2,i_in))
  coup2(3) = Conjg(Cpl_NNZ_L(i1,i_in))
  coup2(4) = Conjg(Cpl_NNZ_R(i1,i_in))
  coup2(5) = Cpl_CNW_R(i3,i1)
  coup2(6) = Cpl_CNW_L(i3,i1)
  coup2(7) = Conjg(Cpl_CCZ_R(i3,i2))
  coup2(8) = Conjg(Cpl_CCZ_L(i3,i2))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZW8, n_ZW8, resC, check)
  !-----------
  ! u-s channel
  !-----------
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  coup2(1) = Cpl_CNW_L(i3,i_in)
  coup2(2) = Cpl_CNW_R(i3,i_in)
  coup2(3) = Conjg(Cpl_NNZ_L(i1,i_in))
  coup2(4) = Conjg(Cpl_NNZ_R(i1,i_in))
  coup2(5) = Conjg(Cpl_CNW_L(i2,i1))
  coup2(6) = Conjg(Cpl_CNW_R(i2,i1))
  coup2(7) = Conjg(Cpl_CCZ_L(i2,i3))
  coup2(8) = Conjg(Cpl_CCZ_R(i2,i3))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZW8, n_ZW8, resCa, check)
  gNCCSum(Isum) = 2._dp * Real(resC-resCa,dp)
  Contribution(Isum) = 'Z W'

  !-------------------------------
  ! Z boson - S^0 boson  s-channel
  !-------------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  coup2(1) = Cpl_NNZ_L(i1,i_in)
  coup2(2) = Cpl_NNZ_R(i1,i_in)
  coup2(5) = Cpl_CCZ_L(i2,i3)
  coup2(6) = Cpl_CCZ_R(i2,i3)
  Do n1=1,n_S0 
   Isum = Isum + 1
   Boson4(3) = mS0(n1)
   Boson4(4) = gS0(n1)
   coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n1))
   coup2(7) = Conjg(cpl_CCS0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_CCS0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
   gNCCSum(Isum) = 2._dp * Real(resC,dp)
   Contribution(Isum) = 'Z S0_'//Bu(n1)
  End Do  ! n1 n_S0

  !-------------------------------
  ! Z boson - P^0 boson  s-channel
  !-------------------------------
  Do n1=1,n_P0 
   Isum = Isum + 1
   Boson4(3) = mP0(n1)
   Boson4(4) = gP0(n1)
   coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n1))
   coup2(7) = Conjg(cpl_CCP0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_CCP0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS04, n_ZS04, resC, check)
   gNCCSum(Isum) = 2._dp * Real(resC,dp)
   Contribution(Isum) = 'Z P0_'//Bu(n1)
  End Do  ! n1 n_P0

  !---------------------------------
  ! Z boson - S+ boson  t-u channel
  !---------------------------------
  coup2(1) = Cpl_NNZ_L(i1, i_in)
  coup2(2) = Cpl_NNZ_R(i1, i_in)
  Do n1=1,n_Spm  
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-------------
   ! s-t channel
   !-------------
   mass(2) = mC(i2)
   mass(3) = -mC(i3)
   mass(4) = mN(i1)
   coup2(3) = Cpl_SmpCN_L(n1,i2,i_in)
   coup2(4) = Cpl_SmpCN_R(n1,i2,i_in)
   coup2(5) = Cpl_CCZ_R(i3,i2)
   coup2(6) = Cpl_CCZ_L(i3,i2)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i3,i1))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i3,i1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS08, n_ZS08, resC, check)
   !-------------
   ! u-t channel
   !-------------
   mass(2) = mC(i3)
   mass(3) = -mC(i2)
   mass(4) = mN(i1)
   coup2(3) = Conjg(Cpl_SmpCN_R(n1,i3,i_in))
   coup2(4) = Conjg(Cpl_SmpCN_L(n1,i3,i_in))
   coup2(5) = Cpl_CCZ_L(i3,i2)
   coup2(6) = Cpl_CCZ_R(i3,i2)
   coup2(7) = Cpl_SmpCN_L(n1,i2,i1)
   coup2(8) = Cpl_SmpCN_R(n1,i2,i1)
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS08, n_ZS08, resCa, check)
   gNCCSum(Isum) = 2._dp * Real(resC-resCa,dp)
   Contribution(Isum) = 'Z S-_'//Bu(n1)
  End Do ! n1 S^+m

  !-------------------------------
  ! W boson - W-boson
  !-------------------------------
  Boson4(1) = mW
  Boson4(2) = gW
  Boson4(3) = mW
  Boson4(4) = gW
  Isum = Isum + 1

  mass(2) = mC(i3)
  mass(3) = -mN(i1)
  mass(4) = mC(i2)
  coup2(1) = Conjg(Cpl_CNW_R(i2,i_in))
  coup2(2) = Conjg(Cpl_CNW_L(i2,i_in))
  coup2(3) = Conjg(Cpl_CNW_L(i3,i_in))
  coup2(4) = Conjg(Cpl_CNW_R(i3,i_in))
  coup2(5) = Cpl_CNW_L(i3,i1)
  coup2(6) = Cpl_CNW_R(i3,i1)
  coup2(7) = Cpl_CNW_R(i2,i1)
  coup2(8) = Cpl_CNW_L(i2,i1)
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsW8, n_W8, resC, check)
  gNCCSum(Isum) = -2._dp * Real(resC,dp)
  Contribution(Isum) = 'W W st'

  !-------------------------------
  ! W boson - S+ boson  s-channel
  !-------------------------------
  Do n1=1,n_Spm 
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-----------------
   ! t-channel
   !-----------------
   mass(2) = mC(i2)
   mass(3) = -mN(i1)
   mass(4) = mC(i3)
   coup2(1) = Conjg(Cpl_CNW_R(i2,i_in))
   coup2(2) = Conjg(Cpl_CNW_L(i2,i_in))
   coup2(3) = Cpl_SmpCN_L(n1,i2,i_in)
   coup2(4) = Cpl_SmpCN_R(n1,i2,i_in)
   coup2(5) = Cpl_CNW_L(i3,i1)
   coup2(6) = Cpl_CNW_R(i3,i1)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i3,i1))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i3,i1))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm4, n_WSpm4, resC, check)
   !-----------------
   ! u-channel
   !-----------------
   mass(2) = mC(i3)
   mass(3) = -mC(i2)
   mass(4) = mN(i1)
   coup2(1) = Cpl_CNW_L(i3,i_in)
   coup2(2) = Cpl_CNW_R(i3,i_in)
   coup2(3) = Conjg(Cpl_SmpCN_R(n1,i3,i_in))
   coup2(4) = Conjg(Cpl_SmpCN_L(n1,i3,i_in))
   coup2(5) = Conjg(Cpl_CNW_L(i2,i1))
   coup2(6) = Conjg(Cpl_CNW_R(i2,i1))
   coup2(7) = Cpl_SmpCN_L(n1,i2,i1)
   coup2(8) = Cpl_SmpCN_R(n1,i2,i1)
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm4, n_WSpm4, resCa, check)
   gNCCSum(Isum) = 2._dp * Real(resCa-resC,dp)
   Contribution(Isum) = 'W S-_'//Bu(n1)//' s-s'
  End Do  ! n1 Spm(n1)

  !---------------------------------
  ! W boson - S+ boson  t-u channel
  !---------------------------------
  Do n1=1,n_Spm  
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-------------
   ! t-u channel
   !-------------
   mass(2) = mC(i3)
   mass(3) = -mN(i1)
   mass(4) = mC(i2)
   coup2(1) = Conjg(Cpl_CNW_R(i2, i_in))
   coup2(2) = Conjg(Cpl_CNW_L(i2, i_in))
   coup2(3) = Conjg(Cpl_SmpCN_R(n1,i3,i_in))
   coup2(4) = Conjg(Cpl_SmpCN_L(n1,i3,i_in))
   coup2(5) = Cpl_CNW_L(i3,i1)
   coup2(6) = Cpl_CNW_R(i3,i1)
   coup2(7) = Cpl_SmpCN_L(n1,i2,i1)
   coup2(8) = Cpl_SmpCN_R(n1,i2,i1)
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-t channel
   !-------------
   mass(2) = mC(i2)
   mass(3) = -mN(i1)
   mass(4) = mC(i3)
   coup2(1) = Cpl_CNW_L(i3, i_in)
   coup2(2) = Cpl_CNW_R(i3, i_in)
   coup2(3) = Cpl_SmpCN_L(n1,i2,i_in)
   coup2(4) = Cpl_SmpCN_R(n1,i2,i_in)
   coup2(5) = Conjg(Cpl_CNW_R(i2,i1))
   coup2(6) = Conjg(Cpl_CNW_L(i2,i1))
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i3,i1))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i3,i1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gNCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W S-_'//Bu(n1)//' s-t'
  End Do ! n1 S^+m

  !---------------------------------
  ! W boson - S0 boson  t-s channel
  !---------------------------------
  Do n1=1,n_S0  
   Isum = Isum + 1
   Boson4(3) = mS0(n1)
   Boson4(4) = gS0(n1)
   !-------------
   ! t-s channel
   !-------------
   mass(2) = mN(i1)
   mass(3) = -mC(i3)
   mass(4) = mC(i2)
   coup2(1) = Conjg(Cpl_CNW_R(i2, i_in))
   coup2(2) = Conjg(Cpl_CNW_L(i2, i_in))
   coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_R(i3,i1)
   coup2(6) = Cpl_CNW_L(i3,i1)
   coup2(7) = Conjg(cpl_CCS0_R(i3,i2,n1))
   coup2(8) = Conjg(cpl_CCS0_L(i3,i2,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-s channel
   !-------------
   mass(2) = mN(i1)
   mass(3) = -mC(i2)
   mass(4) = mC(i3)
   coup2(1) = Cpl_CNW_L(i3, i_in)
   coup2(2) = Cpl_CNW_R(i3, i_in)
   coup2(5) = Conjg(Cpl_CNW_L(i2,i1))
   coup2(6) = Conjg(Cpl_CNW_R(i2,i1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gNCCSum(Isum) = -2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W S0_'//Bu(n1)//' s-t'
  End Do ! n1 S^0

  !---------------------------------
  ! W boson - P0 boson  t-s channel
  !---------------------------------
  Do n1=1,n_P0  
   Isum = Isum + 1
   Boson4(3) = mP0(n1)
   Boson4(4) = gP0(n1)
   !-------------
   ! t-s channel
   !-------------
   mass(2) = mN(i1)
   mass(3) = -mC(i3)
   mass(4) = mC(i2)
   coup2(1) = Conjg(Cpl_CNW_R(i2, i_in))
   coup2(2) = Conjg(Cpl_CNW_L(i2, i_in))
   coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_R(i3,i1)
   coup2(6) = Cpl_CNW_L(i3,i1)
   coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n1))
   coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-s channel
   !-------------
   mass(2) = mN(i1)
   mass(3) = -mC(i2)
   mass(4) = mC(i3)
   coup2(1) = Cpl_CNW_L(i3, i_in)
   coup2(2) = Cpl_CNW_R(i3, i_in)
   coup2(5) = Conjg(Cpl_CNW_L(i2,i1))
   coup2(6) = Conjg(Cpl_CNW_R(i2,i1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gNCCSum(Isum) = -2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W P0_'//Bu(n1)//' s-t'
  End Do ! n1 P^0

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = n1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n2))
    coup2(5) = cpl_CCS0_L(i2,i3,n1)
    coup2(6) = cpl_CCS0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_CCS0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_CCS0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gNCCSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//'S0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    coup2(1) = cpl_NNS0_L(i1,i_in,n1)
    coup2(2) = cpl_NNS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_CCS0_L(i2,i3,n1)
    coup2(6) = cpl_CCS0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gNCCSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^+_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = 1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mC(i3)
    mass(4) = mC(i2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i2,i_in))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i2,i_in))
    coup2(3) = Conjg(cpl_NNS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNS0_L(i1,i_in,n2))
    coup2(5) = Cpl_SmpCN_L(n1,i3,i1)
    coup2(6) = Cpl_SmpCN_R(n1,i3,i1)
    coup2(7) = Conjg(cpl_CCS0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCS0_L(i3,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_SmpCN_L(n1,i3,i_in)
    coup2(2) = Cpl_SmpCN_R(n1,i3,i_in)
    coup2(5) = Conjg(Cpl_SmpCN_R(n1,i2,i1))
    coup2(6) = Conjg(Cpl_SmpCN_L(n1,i2,i1))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    gNCCSum(Isum) = - 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'S0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  mass(2) = mN(i1)
  mass(3) = -mC(i2)
  mass(4) = mC(i3)
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   Do n2 = n1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    coup2(1) = cpl_NNP0_L(i1,i_in,n1)
    coup2(2) = cpl_NNP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = cpl_CCP0_L(i2,i3,n1)
    coup2(6) = cpl_CCP0_R(i2,i3,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gNCCSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do


  !-------------------------------------------------
  ! P^0_i boson  s-channel - S^+_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mC(i3)
    mass(4) = mC(i2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i2,i_in))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i2,i_in))
    coup2(3) = Conjg(cpl_NNP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_NNP0_L(i1,i_in,n2))
    coup2(5) = Cpl_SmpCN_L(n1,i3,i1)
    coup2(6) = Cpl_SmpCN_R(n1,i3,i1)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mN(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_SmpCN_L(n1,i3,i_in)
    coup2(2) = Cpl_SmpCN_R(n1,i3,i_in)
    coup2(5) = Conjg(Cpl_SmpCN_R(n1,i2,i1))
    coup2(6) = Conjg(Cpl_SmpCN_L(n1,i2,i1))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    gNCCSum(Isum) = - 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------
  ! S+ boson  s-channel - S+ boson  s-channel
  !-------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = n1+1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(n2)
    Boson4(4) = gSpm(n2)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i2)
    mass(3) = -mN(i1)
    mass(4) = mC(i3)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i2,i_in))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i2,i_in))
    coup2(3) = Cpl_SmpCN_L(n2,i2,i_in)
    coup2(4) = Cpl_SmpCN_R(n2,i2,i_in)
    coup2(5) = Cpl_SmpCN_L(n1,i3,i1)
    coup2(6) = Cpl_SmpCN_R(n1,i3,i1)
    coup2(7) = Conjg(Cpl_SmpCN_R(n2,i3,i1))
    coup2(8) = Conjg(Cpl_SmpCN_L(n2,i3,i1))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mN(i1)
    coup2(1) = Cpl_SmpCN_L(n1,i3,i_in)
    coup2(2) = Cpl_SmpCN_R(n1,i3,i_in)
    coup2(3) = Conjg(Cpl_SmpCN_R(n2,i3,i_in))
    coup2(4) = Conjg(Cpl_SmpCN_L(n2,i3,i_in))
    coup2(5) = Conjg(Cpl_SmpCN_R(n1,i2,i1))
    coup2(6) = Conjg(Cpl_SmpCN_L(n1,i2,i1))
    coup2(7) = Cpl_SmpCN_L(n2,i2,i1)
    coup2(8) = Cpl_SmpCN_R(n2,i2,i1)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    gNCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'S-_'//Bu(n2)//' s-s'
   End Do
  End Do

  !---------------------------------------------
  ! S+ boson  s-channel - S+ boson  t-channel
  !---------------------------------------------
  mass(2) = mC(i3)
  mass(3) = -mN(i1)
  mass(4) = mC(i2)
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 =1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(n2)
    Boson4(4) = gSpm(n2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i2,i_in))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i2,i_in))
    coup2(3) = Conjg(Cpl_SmpCN_R(n2,i3,i_in))
    coup2(4) = Conjg(Cpl_SmpCN_L(n2,i3,i_in))
    coup2(5) = Cpl_SmpCN_L(n1,i3,i1)
    coup2(6) = Cpl_SmpCN_R(n1,i3,i1)
    coup2(7) = Cpl_SmpCN_L(n2,i2,i1)
    coup2(8) = Cpl_SmpCN_R(n2,i2,i1)
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    gNCCSum(Isum) = - 2._dp * Real(resC,dp) 
    Contribution(Isum) = 'S-_'//Bu(n1)//'S-_'//Bu(n2)//' s-t'
   End Do
  End Do

  !----------
  ! Summing
  !----------
  gNCC = sum( gNCCSum(1:Isum) )
  If (gNCC.Lt.0._dp) Then
   p_test = Abs(gNCC) / Maxval(Abs(gNCCSum(1:Isum)))
   gNCC = 0._dp
   If (p_test.gt.prec) then 
    Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//' Chi^-_'//Bu(i2)// &
     & ' Chi^-_'//Bu(i3)//') < 0 :',i_in,i1,i2,i3,gNCC
    Write(ErrCan,*) 'The different contributions are :'
    Do n1=1,Isum
     If (gNCCSum(n1).Ne.0._dp)  Write(ErrCan,*) Contribution(n1),gNCCSum(n1)
    End Do
   End If
  End If

  gNCC = fac * gNCC

  !---------------------------
  ! for detailed information
  !---------------------------
  If (WriteContribution.ne.0) Then

   gNCCSum = gNCCSum * fac
 
   Write (n_out,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//' Chi^-_'//Bu(i2)// &
     & ' Chi^-_'//Bu(i3)//') :',i_in,i1,i2,i3
   Write (n_out,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gNCCSum(n1).Ne.0._dp)  Write (n_out,*) Contribution(n1),gNCCSum(n1)
   End Do
  End If

  Deallocate( gNCCSum, Contribution )
   
  Iname = Iname - 1

 End Subroutine Chi0ToChi0ChimChip

 Subroutine Chi0ToChi0ff(i_in, i_out, state, mN, mZ, gZ, Cpl_NNZ_L, Cpl_NNZ_R &
    & , n_f, mf, L_f, R_f, IntegralsZ4, n_Z4, mS0, gS0, cpl_NNS0_L            &
    & , cpl_NNS0_R, cpl_FFS0_L, cpl_FFS0_R, IntegralsS04, n_S04, mP0, gP0     &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_FFP0_L, cpl_FFP0_R, mSf, gSf, cpl_FNSf_L  &
    & , cpl_FNSf_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04, IntegralsZSf8 &
    & , n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8, n_S0Sf8              &
    & , IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI              &
    & , GenerationMixing, check, fac, gNff, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a neutralino to another Neutralino + fermion pair
 ! Written by Werner Porod, 04.06.2001
 ! 12.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Character(len=5), Intent(in) :: state 
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_Z4, n_S04, n_Sf4, n_ZS04, n_ZSf8, n_S0P04 &
     & , n_S0Sf8, n_CSf4, n_Sf8
  Logical, Intent(in) :: check
  Real(dp), Intent(in) :: mN(:), mZ, gZ, mf(:), L_f, R_f, mS0(:), gS0(:) &
      & , mP0(:), gP0(:), mSf(:), gSf(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsS04(:,:)         &
      & , IntegralsSf4(:,:)
  Complex(dp), Intent(in) :: Cpl_NNZ_L(:,:), Cpl_NNZ_R(:,:)             &
      & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_FFS0_L(:,:,:)       &
      & , cpl_FFS0_R(:,:,:), cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:)       &
      & , cpl_FFP0_L(:,:,:), cpl_FFP0_R(:,:,:), cpl_FNSf_L(:,:,:)       &
      & , cpl_FNSf_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsZSf8(:,:)  &
      & , IntegralsS0P04(:,:), IntegralsS0Sf8(:,:), IntegralsCSf4(:,:)  &
      & , IntegralsSf8(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gNff(:,:)

  Integer, Intent(in) :: n_out, WriteContribution

  Integer :: Isum, i1, i2, i3, i4, n_S0, n_P0, i_gen
  Real(dp) :: Boson2(2), mass(4), resR, resRa, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8), resCa

  Real(dp), Allocatable :: gNffSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0ToChi0ff'

  mass(1) = mN(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )


  Isum = (1 + n_S0 + n_P0 + 2*n_f)**2

  Allocate( gNffSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )
   
  gNffSum = 0._dp
  Contribution = ' '
  !----- 
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  coup1(1) = Cpl_NNZ_L(i_out,i_in)
  coup1(2) = Cpl_NNZ_R(i_out,i_in)
  coup1(3) = L_f
  coup1(4) = R_f

  Do i1=1,n_f
   mass(2) = mN(i_out)
   mass(3) = -mf(i1)
   mass(4) = mf(i1)
   Contribution(i1,i1,1) = 'Z f_'//Bu(i1)//' f_'//Bu(i1)
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, gNffSum(i1,i1,Isum), check )
  End Do

  !-------------------
  ! S_0 S_0, diagonal
  !-------------------
  Do i1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(i1)
   Boson2(2) = gS0(i1)
   coup1(1) = cpl_NNS0_L(i_out,i_in,i1) 
   coup1(2) = cpl_NNS0_R(i_out,i_in,i1)
   Do i2=1,n_f
    coup1(3) = cpl_FFS0_L(i2,i2,i1)
    coup1(4) = cpl_FFS0_R(i2,i2,i1)
    mass(2) = mN(i_out)
    mass(3) = - mf(i2)
    mass(4) = mf(i2)
    Contribution(i2,i2,Isum) = 'S0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, gNffSum(i2, i2, Isum), check )
   End Do
  End Do

  !-------------------
  ! P_0 P_0, diagonal
  !-------------------
  Do i1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(i1)
   Boson2(2) = gP0(i1)
   coup1(1) = cpl_NNP0_L(i_out,i_in,i1) 
   coup1(2) = cpl_NNP0_R(i_out,i_in,i1)
   Do i2=1,n_f
    coup1(3) = cpl_FFP0_L(i2,i2,i1)
    coup1(4) = cpl_FFP0_R(i2,i2,i1)
    mass(2) = mN(i_out)
    mass(3) = - mf(i2)
    mass(4) = mf(i2)
    Contribution(i2,i2,Isum) = 'P0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, gNffSum(i2, i2, Isum), check )
   End Do
  End Do

  !-----------------------------
  ! sfermion sfermion, diagonal
  !-----------------------------
  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)
    Do i1=1,n_f
     Do i3=1,n_f
      !------------------------
      ! t-channel
      !------------------------
      coup1(1) = cpl_FNSf_L(i1,i_in,i2)
      coup1(2) = cpl_FNSf_R(i1,i_in,i2)
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
      coup1(1) = Conjg(cpl_FNSf_R(i3,i_in,i2))  ! u-channel
      coup1(2) = Conjg(cpl_FNSf_L(i3,i_in,i2))
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
   Do i2=1,2*n_f  ! fermion generation
    i1 = (i2+1)/2
     Isum = Isum + 1
     Boson2(1) = mSf(i2)
     Boson2(2) = gSf(i2)
     !------------------------
     ! t-channel
     !------------------------
     coup1(1) = cpl_FNSf_L(i1,i_in,i2)
     coup1(2) = cpl_FNSf_R(i1,i_in,i2)
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
     coup1(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
     coup1(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
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

  !--------
  ! Z - S0
  !--------
  Boson4(1) = mZ
  Boson4(2) = gZ
  Do i1=1,n_S0 
   Isum = Isum + 1
   Boson4(3) = mS0(i1)
   Boson4(4) = gS0(i1)
   coup2(1) = Cpl_NNZ_L(i_out,i_in)
   coup2(2) = Cpl_NNZ_R(i_out,i_in)
   coup2(3) = Conjg(cpl_NNS0_R(i_out,i_in,i1))
   coup2(4) = Conjg(cpl_NNS0_L(i_out,i_in,i1))
   coup2(5) = L_f
   coup2(6) = R_f
   Do i2=1,n_f
    mass(2) = mN(i_out)
    mass(3) = -mf(i2)
    mass(4) = mf(i2)
    coup2(7) = Conjg(cpl_FFS0_R(i2,i2,i1))
    coup2(8) = Conjg(cpl_FFS0_L(i2,i2,i1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsZS04, n_ZS04, resC, check)
    gNffSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
    Contribution(i2,i2,Isum) = 'Z S0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End Do

  !--------
  ! Z - P0
  !--------
  Do i1=1,n_P0 
   Isum = Isum + 1
   Boson4(3) = mP0(i1)
   Boson4(4) = gP0(i1)
   coup2(1) = Cpl_NNZ_L(i_out,i_in)
   coup2(2) = Cpl_NNZ_R(i_out,i_in)
   coup2(3) = Conjg(cpl_NNP0_R(i_out,i_in,i1))
   coup2(4) = Conjg(cpl_NNP0_L(i_out,i_in,i1))
   coup2(5) = L_f
   coup2(6) = R_f
   Do i2=1,n_f
    mass(2) = mN(i_out)
    mass(3) = -mf(i2)
    mass(4) = mf(i2)
    coup2(7) = Conjg(cpl_FFP0_R(i2,i2,i1))
    coup2(8) = Conjg(cpl_FFP0_L(i2,i2,i1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsZS04, n_ZS04, resC, check)
    gNffSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
    Contribution(i2,i2,Isum) = 'Z P0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End Do

  !--------------------------
  ! Z boson - Sfermion_{xyz} 
  !-------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)
    Do i1=1,n_f! generation
     !-----------
     ! t-channel
     !-----------
     mass(2) = mf(i1)
     mass(3) = -mf(i1)
     mass(4) = mN(i_out)
     coup2(1) = Cpl_NNZ_L(i_out, i_in)
     coup2(2) = Cpl_NNZ_R(i_out, i_in)
     coup2(3) = Conjg( cpl_FNSf_R(i1, i_in, i2) )
     coup2(4) = Conjg( cpl_FNSf_L(i1, i_in, i2) )
     coup2(5) = L_f
     coup2(6) = R_f
     coup2(7) = cpl_FNSf_L(i1, i_out, i2)
     coup2(8) = cpl_FNSf_R(i1, i_out, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     !-----------
     ! u-channel
     !-----------
     coup2(3) = cpl_FNSf_L(i1, i_in, i2)
     coup2(4) = cpl_FNSf_R(i1, i_in, i2)
     coup2(5) = R_f
     coup2(6) = L_f
     coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
     coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resCa, check)
     gNffSum(i1,i1,Isum) = 2._dp *Real(resCa-resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  Else

   Do i1=1,n_f! generation
    i_gen = i1*2-1
    Do i2=i_gen,i_gen+1
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     !-----------
     ! t-channel
     !-----------
     mass(2) = mf(i1)
     mass(3) = -mf(i1)
     mass(4) = mN(i_out)
     coup2(1) = Cpl_NNZ_L(i_out, i_in)
     coup2(2) = Cpl_NNZ_R(i_out, i_in)
     coup2(3) = Conjg( cpl_FNSf_R(i1, i_in, i2) )
     coup2(4) = Conjg( cpl_FNSf_L(i1, i_in, i2) )
     coup2(5) = L_f
     coup2(6) = R_f
     coup2(7) = cpl_FNSf_L(i1, i_out, i2)
     coup2(8) = cpl_FNSf_R(i1, i_out, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     !-----------
     ! u-channel
     !-----------
     coup2(3) = cpl_FNSf_L(i1, i_in, i2)
     coup2(4) = cpl_FNSf_R(i1, i_in, i2)
     coup2(5) = R_f
     coup2(6) = L_f
     coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
     coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resCa, check)
     gNffSum(i1,i1,Isum) = 2._dp *Real(resCa-resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  End If

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_S0
   Boson4(1) = mS0(i1)
   Boson4(2) = gS0(i1)
   Do i2 = i1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(i2)
    Boson4(4) = gS0(i2)
    coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
    coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_NNS0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_NNS0_L(i_out,i_in,i2))
    Do i3=1,n_f
     mass(2) = mN(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFS0_L(i3,i3,i1)
     coup2(6) = cpl_FFS0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFS0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFS0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gNffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
       &  'S0_'//Bu(i1)//' S0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
    End Do
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_S0
   Boson4(1) = mS0(i1)
   Boson4(2) = gS0(i1)
   Do i2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(i2)
    Boson4(4) = gP0(i2)
    coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
    coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_NNP0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_NNP0_L(i_out,i_in,i2))
    Do i3=1,n_f
    Contribution(i3,i3,Isum) = &
       &  'S0_'//Bu(i1)//' P0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     mass(2) = mN(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFS0_L(i3,i3,i1)
     coup2(6) = cpl_FFS0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFP0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFP0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gNffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    End Do
   End Do
  End Do

  !--------------------------------------
  ! S0 boson  s-channel - S-fermion_{xyz} 
  !--------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_S0
    Boson4(1) = mS0(i1)
    Boson4(2) = gS0(i1)
    Do i2 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 =1,n_f
      !-----------
      ! t-channel
      !-----------
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mN(i_out)
      coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
      coup2(3) = Conjg(cpl_FNSf_R(i3,i_in, i2))
      coup2(4) = Conjg(cpl_FNSf_L(i3,i_in, i2))
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = cpl_FNSf_L(i3,i_out, i2)
      coup2(8) = cpl_FNSf_R(i3,i_out, i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      !-----------
      ! u-channel
      !-----------
      coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
      coup2(3) = cpl_FNSf_L( i3,i_in, i2)
      coup2(4) = cpl_FNSf_R( i3,i_in, i2)
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = Conjg(cpl_FNSf_R(i3,i_out, i2))
      coup2(8) = Conjg(cpl_FNSf_L(i3,i_out, i2))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resCa, check)
      gNffSum(i3,i3,Isum) = - 2._dp * Real(resC+resCa,dp)
      Contribution(i3,i3,Isum) = &
        &  'S0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_S0
    Boson4(1) = mS0(i1)
    Boson4(2) = gS0(i1)
    Do i3 =1,n_f
     i_gen = 2*i3-1
     Do i2 =i_gen,i_gen+1
      Isum = Isum + 1
     Contribution(i3,i3,Isum) = &
        &  'S0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
      Boson4(3) = mSf(i2)
      Boson4(4) = gSf(i2)
      !-----------
      ! t-channel
      !-----------
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mN(i_out)
      coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
      coup2(3) = Conjg(cpl_FNSf_R(i3,i_in, i2))
      coup2(4) = Conjg(cpl_FNSf_L(i3,i_in, i2))
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = cpl_FNSf_L(i3,i_out, i2)
      coup2(8) = cpl_FNSf_R(i3,i_out, i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      !-----------
      ! u-channel
      !-----------
      coup2(1) = cpl_NNS0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNS0_R(i_out,i_in,i1)
      coup2(3) = cpl_FNSf_L( i3,i_in, i2)
      coup2(4) = cpl_FNSf_R( i3,i_in, i2)
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = Conjg(cpl_FNSf_R(i3,i_out, i2))
      coup2(8) = Conjg(cpl_FNSf_L(i3,i_out, i2))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resCa, check)
      gNffSum(i3,i3,Isum) = - 2._dp * Real(resC+resCa,dp)
     End Do
    End Do
   End Do

  End If

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_P0
   Boson4(1) = mP0(i1)
   Boson4(2) = gP0(i1)
   Do i2 = i1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(i2)
    Boson4(4) = gP0(i2)
    coup2(1) = cpl_NNP0_L(i_out,i_in,i1)
    coup2(2) = cpl_NNP0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_NNP0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_NNP0_L(i_out,i_in,i2))
    Do i3=1,n_f
    Contribution(i3,i3,Isum) = &
       &  'P0_'//Bu(i1)//' P0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     mass(2) = mN(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFP0_L(i3,i3,i1)
     coup2(6) = cpl_FFP0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFP0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFP0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gNffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    End Do
   End Do
  End Do

  !--------------------------------------
  ! P0 boson  s-channel - S-fermion_{xyz} 
  !--------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_P0
    Boson4(1) = mP0(i1)
    Boson4(2) = gP0(i1)
    Do i2 =1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 =1,n_f
      !-----------
      ! t-channel
      !-----------
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mN(i_out)
      coup2(1) = cpl_NNP0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNP0_R(i_out,i_in,i1)
      coup2(3) = Conjg(cpl_FNSf_R(i3,i_in, i2))
      coup2(4) = Conjg(cpl_FNSf_L(i3,i_in, i2))
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = cpl_FNSf_L(i3,i_out, i2)
      coup2(8) = cpl_FNSf_R(i3,i_out, i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      !-----------
      ! u-channel
      !-----------
      coup2(1) = cpl_NNP0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNP0_R(i_out,i_in,i1)
      coup2(3) = cpl_FNSf_L( i3,i_in, i2)
      coup2(4) = cpl_FNSf_R( i3,i_in, i2)
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = Conjg(cpl_FNSf_R(i3,i_out, i2))
      coup2(8) = Conjg(cpl_FNSf_L(i3,i_out, i2))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resCa, check)
      gNffSum(i3,i3,Isum) = - 2._dp * Real(resC+resCa,dp)
      Contribution(i3,i3,Isum) = &
         &  'P0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_P0
    Boson4(1) = mP0(i1)
    Boson4(2) = gP0(i1)
    Do i3 =1,n_f
     i_gen = 2*i3-1
     Do i2 =i_gen,i_gen+1
      Isum = Isum + 1
      Contribution(i3,i3,Isum) = &
         &  'P0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
      Boson4(3) = mSf(i2)
      Boson4(4) = gSf(i2)
      !-----------
      ! t-channel
      !-----------
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mN(i_out)
      coup2(1) = cpl_NNP0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNP0_R(i_out,i_in,i1)
      coup2(3) = Conjg(cpl_FNSf_R(i3,i_in, i2))
      coup2(4) = Conjg(cpl_FNSf_L(i3,i_in, i2))
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = cpl_FNSf_L(i3,i_out, i2)
      coup2(8) = cpl_FNSf_R(i3,i_out, i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      !-----------
      ! u-channel
      !-----------
      coup2(1) = cpl_NNP0_L(i_out,i_in,i1)
      coup2(2) = cpl_NNP0_R(i_out,i_in,i1)
      coup2(3) = cpl_FNSf_L( i3,i_in, i2)
      coup2(4) = cpl_FNSf_R( i3,i_in, i2)
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = Conjg(cpl_FNSf_R(i3,i_out, i2))
      coup2(8) = Conjg(cpl_FNSf_L(i3,i_out, i2))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resCa, check)
      gNffSum(i3,i3,Isum) = - 2._dp * Real(resC+resCa,dp)
     End Do
    End Do
   End Do

  End If

  !----------------------------------
  ! Sfermion - Sfermion, non-diagonal
  !----------------------------------
  If (GenerationMixing) Then
   Do i3=1,2*n_f-1
    Boson4(1) = mSf(i3)
    Boson4(2) = gSf(i3)
    Do i4=i3+1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i4)
     Boson4(4) = gSf(i4)
     Do i1 = 1,n_f ! fermions
      Do i2 = 1,n_f
      Contribution(i1,i2,Isum) = &
         &  'tt Sf_'//Bu(i3)//' Sf_'//Bu(i4)//' f_'//Bu(i1)//' f_'//Bu(i2)
       !-------------
       ! t-channel
       !-------------
       mass(2) = mf(i1)
       mass(3) = -mf(i2)
       mass(4) = mN(i_out)
       coup2(1) = cpl_FNSf_L(i1,i_in,i3)
       coup2(2) = cpl_FNSf_R(i1,i_in,i3)
       coup2(3) = Conjg(cpl_FNSf_R(i1,i_in,i4))
       coup2(4) = Conjg(cpl_FNSf_L(i1,i_in,i4))
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
       coup2(1) = Conjg(cpl_FNSf_R(i2,i_in,i3))
       coup2(2) = Conjg(cpl_FNSf_L(i2,i_in,i3))
       coup2(3) = cpl_FNSf_L(i2,i_in,i4)
       coup2(4) = cpl_FNSf_R(i2,i_in,i4)
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

   Do i1 = 1,n_f
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
    coup2(1) = cpl_FNSf_L(i1,i_in,i2)
    coup2(2) = cpl_FNSf_R(i1,i_in,i2)
    coup2(3) = Conjg(cpl_FNSf_R(i1,i_in,i3))
    coup2(4) = Conjg(cpl_FNSf_L(i1,i_in,i3))
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
    coup2(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
    coup2(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
    coup2(3) = cpl_FNSf_L(i1,i_in,i3)
    coup2(4) = cpl_FNSf_R(i1,i_in,i3)
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
   Do i2= 1,2*n_f
    Boson4(1) = mSf(i2)
    Boson4(2) = gSf(i2)
    Do i3 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i3)
     Boson4(4) = gSf(i3)
     Do i1 = 1,n_f ! fermion
      coup2(1) = cpl_FNSf_L(i1,i_in,i2)
      coup2(2) = cpl_FNSf_R(i1,i_in,i2)
      coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
      coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
      Do i4=1,n_f
       Contribution(i1,i4,Isum) = &
          &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i4)
       mass(2) = mf(i4)
       mass(3) = -mN(i_out)
       mass(4) = mf(i1)
       coup2(3) = cpl_FNSf_L(i4,i_in, i3)
       coup2(4) = cpl_FNSf_R(i4,i_in, i3)
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

   Do i1 = 1,n_f ! fermions
    i_gen = 2*i1-1
    Do i2 = i_gen, i_gen+1
     Boson4(1) = mSf(i2)
     Boson4(2) = gSf(i2)
     coup2(1) = cpl_FNSf_L(i1,i_in,i2)
     coup2(2) = cpl_FNSf_R(i1,i_in,i2)
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
      coup2(3) = cpl_FNSf_L(i1,i_in, i3)
      coup2(4) = cpl_FNSf_R(i1,i_in, i3)
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
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gNff(i1,i2) = Sum( gNffSum(i1,i2,1:Isum) )
     If (gNff(i1,i2).Lt.0._dp) Then
      p_test = Abs(gNff(i1,i2)) / Maxval(Abs(gNffSum(i1,i2,1:Isum)))
      gNff(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero
      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//state//') < 0 :' &
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
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//state//') < 0 :' &
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
  If (WriteContribution.Ne.0) Then

   gNffSum = gNffSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//state//') :' &
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
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//state//') :' &
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

  Deallocate( gNffSum, Contribution)

  Iname = Iname - 1

 End Subroutine Chi0ToChi0ff


 Subroutine Chi0ToChi0NuNu(i_in, i_out, mN, mZ, gZ, Cpl_NNZ_L                 &
    & , Cpl_NNZ_R, n_f, L_f, R_f, IntegralsZ4, n_Z4, mSf, gSf, cpl_FNSf_L     &
    & , cpl_FNSf_R, IntegralsSf4, n_Sf4, IntegralsZSf8, n_ZSf8, IntegralsCSf4 &
    & , n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI, GenerationMixing, check    &
    & , fac, gNff, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a neutralino to another Neutralino + neutrino pair
 ! Written by Werner Porod, 24.06.2001
 ! 08.07.2001: in the sfermion part: t- is equal u-channel 
 ! 12.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_Z4, n_Sf4, n_ZSf8, n_CSf4, n_Sf8
  Logical, Intent(in) :: check

  Real(dp), Intent(in) :: mN(:), mZ, gZ, L_f, R_f, mSf(:), gSf(:), deltaM &
      & , epsI, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsSf4(:,:)
  Complex(dp), Intent(in) :: Cpl_NNZ_L(:,:), Cpl_NNZ_R(:,:)             &
      & , cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsZSf8(:,:), IntegralsCSf4(:,:)  &
      & , IntegralsSf8(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gNff(:,:)

  Integer, Intent(in) :: n_out, WriteContribution

  Integer :: Isum, i1, i2, i3, i4
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8), resCa

  Real(dp), Allocatable :: gNffSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0ToChi0NuNu'

  mass(1) = mN(i_in)

  Isum = (1 + 2*n_f)**2

  Allocate( gNffSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )
   
  gNffSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  coup1(1) = Cpl_NNZ_L(i_out,i_in)
  coup1(2) = Cpl_NNZ_R(i_out,i_in)
  coup1(3) = L_f
  coup1(4) = R_f
  mass(2) = mN(i_out)
  mass(3) = 0._dp
  mass(4) = 0._dp
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check )
  Do i1=1,n_f
   gNffSum(i1,i1,Isum) = resR
   Contribution(i1,i1,1) = 'Z f_'//Bu(i1)//' f_'//Bu(i1)
  End Do

  !-----------------------------
  ! sfermion sfermion, diagonal
  !-----------------------------
  If (GenerationMixing) Then
   mass(2) = 0._dp
   Do i2=1,n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)
    Do i1=1,n_f
     Do i3=1,n_f
      !------------------------
      ! t-channel
      !------------------------
      coup1(1) = cpl_FNSf_L(i1,i_in,i2)
      coup1(2) = cpl_FNSf_R(i1,i_in,i2)
      coup1(3) = Conjg(cpl_FNSf_R(i3,i_out,i2))
      coup1(4) = Conjg(cpl_FNSf_L(i3,i_out,i2))
      mass(3) = 0._dp
      mass(4) = mN(i_out)        
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      !------------------------
      ! u-channel
      !------------------------
!      coup1(1) = Conjg(cpl_FNSf_R(i3,i_in,i2))  ! u-channel
!      coup1(2) = Conjg(cpl_FNSf_L(i3,i_in,i2))
!      coup1(3) = cpl_FNSf_L(i1,i_out,i2)
!      coup1(4) = cpl_FNSf_R(i1,i_out,i2)
!      mass(3) = - mN(i_out)        
!      mass(4) = 0._dp
!      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
!                           &, IntegralsSf4, n_Sf4, resRa, check)
      gNffSum(i1,i3,Isum) = 2._dp * resR ! + resRa
      Contribution(i1,i3,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i3)
     End Do ! i3 fermions
    End Do  ! i1 fermions
   End Do   ! i2 sfermions

  Else
   mass(2) = 0._dp
   Do i2=1,n_f  ! fermion generation
    i1 = i2
     Isum = Isum + 1
     Boson2(1) = mSf(i2)
     Boson2(2) = gSf(i2)
     !------------------------
     ! t-channel
     !------------------------
     coup1(1) = cpl_FNSf_L(i1,i_in,i2)
     coup1(2) = cpl_FNSf_R(i1,i_in,i2)
     coup1(3) = Conjg(cpl_FNSf_R(i1,i_out,i2))
     coup1(4) = Conjg(cpl_FNSf_L(i1,i_out,i2))
     mass(3) = 0._dp
     mass(4) = mN(i_out)        
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resR, check)
     !------------------------
     ! u-channel
     !------------------------
!     coup1(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
!     coup1(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
!     coup1(3) = cpl_FNSf_L(i1,i_out,i2)
!     coup1(4) = cpl_FNSf_R(i1,i_out,i2)
!     mass(3) = - mN(i_out)        
!     mass(4) = 0._dp
!     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
!                          &, IntegralsSf4, n_Sf4, resRa, check)
     gNffSum(i1,i1,Isum) = 2._dp * resR ! + resRa
     Contribution(i1,i1,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do   ! i1 u-quarks
  End If    ! GenerationMixing

  !--------------------------
  ! Z boson - Sfermion_{xyz} 
  !-------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  mass(2) = 0._dp
  mass(3) = 0._dp
  mass(4) = mN(i_out)
  If (GenerationMixing) Then
   Do i2=1,n_f! generation
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)
    Do i1=1,n_f
     !-----------
     ! t-channel
     !-----------
     coup2(1) = Cpl_NNZ_L(i_out, i_in)
     coup2(2) = Cpl_NNZ_R(i_out, i_in)
     coup2(3) = Conjg( cpl_FNSf_R(i1, i_in, i2) )
     coup2(4) = Conjg( cpl_FNSf_L(i1, i_in, i2) )
     coup2(5) = L_f
     coup2(6) = R_f
     coup2(7) = cpl_FNSf_L(i1, i_out, i2)
     coup2(8) = cpl_FNSf_R(i1, i_out, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     !-----------
     ! u-channel
     !-----------
     coup2(3) = cpl_FNSf_L(i1, i_in, i2)
     coup2(4) = cpl_FNSf_R(i1, i_in, i2)
     coup2(5) = R_f
     coup2(6) = L_f
     coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
     coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resCa, check)
     gNffSum(i1,i1,Isum) = 2._dp *Real(resCa-resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  Else
   Do i1=1,n_f! generation
    i2 = i1
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     !-----------
     ! t-channel
     !-----------
     coup2(1) = Cpl_NNZ_L(i_out, i_in)
     coup2(2) = Cpl_NNZ_R(i_out, i_in)
     coup2(3) = Conjg( cpl_FNSf_R(i1, i_in, i2) )
     coup2(4) = Conjg( cpl_FNSf_L(i1, i_in, i2) )
     coup2(5) = L_f
     coup2(6) = R_f
     coup2(7) = cpl_FNSf_L(i1, i_out, i2)
     coup2(8) = cpl_FNSf_R(i1, i_out, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     !-----------
     ! u-channel
     !-----------
     coup2(3) = cpl_FNSf_L(i1, i_in, i2)
     coup2(4) = cpl_FNSf_R(i1, i_in, i2)
     coup2(5) = R_f
     coup2(6) = L_f
     coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
     coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resCa, check)
     gNffSum(i1,i1,Isum) = 2._dp *Real(resCa-resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do

  End If

  !----------------------------------
  ! Sfermion - Sfermion, non-diagonal
  !----------------------------------
  If (GenerationMixing) Then
   mass(2) = 0._dp
   Do i3=1,n_f-1
    Boson4(1) = mSf(i3)
    Boson4(2) = gSf(i3)
    Do i4=i3+1,n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i4)
     Boson4(4) = gSf(i4)
     Do i1 = 1,n_f ! fermions
      Do i2 = 1,n_f
       !-------------
       ! t-channel
       !-------------
       mass(3) = 0._dp
       mass(4) = mN(i_out)
       coup2(1) = cpl_FNSf_L(i1,i_in,i3)
       coup2(2) = cpl_FNSf_R(i1,i_in,i3)
       coup2(3) = Conjg(cpl_FNSf_R(i1,i_in,i4))
       coup2(4) = Conjg(cpl_FNSf_L(i1,i_in,i4))
       coup2(5) = Conjg(cpl_FNSf_R(i2,i_out,i3))
       coup2(6) = Conjg(cpl_FNSf_L(i2,i_out,i3))
       coup2(7) = cpl_FNSf_L(i2,i_out,i4)
       coup2(8) = cpl_FNSf_R(i2,i_out,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resC, check)
       !-------------
       ! u-channel
       !-------------
!       mass(3) = -mN(i_out)
!       mass(4) = 0._dp
!       coup2(1) = Conjg(cpl_FNSf_R(i2,i_in,i3))
!       coup2(2) = Conjg(cpl_FNSf_L(i2,i_in,i3))
!       coup2(3) = cpl_FNSf_L(i2,i_in,i4)
!       coup2(4) = cpl_FNSf_R(i2,i_in,i4)
!       coup2(5) = cpl_FNSf_L(i1,i_out,i3)
!       coup2(6) = cpl_FNSf_R(i1,i_out,i3)
!       coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i4))
!       coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i4))
!       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
!                              &, IntegralsCSf4, n_CSf4, resCa, check)
       gNffSum(i1,i2,Isum) = 4._dp * Real(resC,dp) ! 2._dp * Real(resC+resCa,dp)
       Contribution(i1,i2,Isum) = &
         &  'tt Sf_'//Bu(i3)//' Sf_'//Bu(i4)//' f_'//Bu(i1)//' f_'//Bu(i2)
      End Do
     End Do
    End Do
   End Do

  End If

  !-------------------------------------------------------
  ! Sfermion_{xyz}  t-channel - Sfermion_{xyz}  u-channel
  !-------------------------------------------------------
  mass(2) = 0._dp
  mass(3) = -mN(i_out)
  mass(4) = 0._dp
  If (GenerationMixing) Then
   Do i2= 1,n_f
    Boson4(1) = mSf(i2)
    Boson4(2) = gSf(i2)
    Do i3 = 1,n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i3)
     Boson4(4) = gSf(i3)
     Do i1 = 1,n_f ! fermion
      coup2(1) = cpl_FNSf_L(i1,i_in,i2)
      coup2(2) = cpl_FNSf_R(i1,i_in,i2)
      coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
      coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
      Do i4=1,n_f
       coup2(3) = cpl_FNSf_L(i4,i_in, i3)
       coup2(4) = cpl_FNSf_R(i4,i_in, i3)
       coup2(5) = Conjg(cpl_FNSf_R(i4,i_out,i2))
       coup2(6) = Conjg(cpl_FNSf_L(i4,i_out,i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSf8, n_Sf8, resC, check)
       gNffSum(i1,i4,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i1,i4,Isum) = &
          &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,n_f !
    i2 = i1
     Boson4(1) = mSf(i2)
     Boson4(2) = gSf(i2)
     coup2(1) = cpl_FNSf_L(i1,i_in,i2)
     coup2(2) = cpl_FNSf_R(i1,i_in,i2)
     coup2(5) = Conjg(cpl_FNSf_R(i1,i_out,i2))
     coup2(6) = Conjg(cpl_FNSf_L(i1,i_out,i2))
     i3 = i1
      Isum = Isum + 1
      Contribution(i1,i1,Isum) = &
         &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
      Boson4(3) = mSf(i3)
      Boson4(4) = gSf(i3)
      coup2(3) = cpl_FNSf_L(i1,i_in, i3)
      coup2(4) = cpl_FNSf_R(i1,i_in, i3)
      coup2(7) = Conjg(cpl_FNSf_R(i1,i_out,i3))
      coup2(8) = Conjg(cpl_FNSf_L(i1,i_out,i3))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSf8, n_Sf8, resC, check)
      gNffSum(i1,i1,Isum) = - 2._dp * Real(resC,dp)
   End Do

  End If

  !----------
  ! Summing
  !----------
  gNff = 0._dp
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gNff(i1,i2) = Sum( gNffSum(i1,i2,1:Isum) )
     If (gNff(i1,i2).Lt.0._dp) Then
      p_test = Abs(gNff(i1,i2)) / Maxval(Abs(gNffSum(i1,i2,1:Isum)))
      gNff(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero
      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//' nu nu ) < 0 :' &
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
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//' nu nu ) < 0 :' &
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
  If (WriteContribution.ne.0) Then

   gNffSum = gNffSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//' nu nu ) :' &
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
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//' nu nu ) :' &
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

  Deallocate( gNffSum, Contribution)

  Iname = Iname - 1

 End Subroutine Chi0ToChi0NuNu

 Subroutine Chi0ToChi0Photon(i_in, i_out, mN, mW2, mC, Cpl_CNW_L, Cpl_CNW_R &
       & , mSpm2, Cpl_SmpCN_L, Cpl_SmpCN_R, mf_u, mUsquark2, cpl_UNSu_L     &
       & , cpl_UNSu_R, mf_d, mDsquark2, cpl_DNSd_L, cpl_DNSd_R              &
       & , mf_l, mSlepton2, cpl_LNSl_L, cpl_LNSl_R, factor, gPhoton)
 !---------------------------------------------------------------------
 ! Calculates the decay of a neutralino to another Neutralino + Photon
 ! Written by Werner Porod, 04.06.2001
 ! - 27.08.03: including e_d, e_u in case of (s)quarks
 !---------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in, i_out
  Real(dp), Intent(in) :: mN(:), mW2, factor(3), mC(:),  mSpm2(:), mf_u(3) &
     & , mUsquark2(6), mf_d(3), mDsquark2(6), mf_l(:), mSlepton2(:)
  Complex(dp), Intent(in) :: Cpl_CNW_L(:,:), Cpl_CNW_R(:,:)        &
     & , Cpl_SmpCN_L(:,:,:), Cpl_SmpCN_R(:,:,:), cpl_UNSu_L(:,:,:) &
     & , cpl_UNSu_R(:,:,:), cpl_DNSd_L(:,:,:), cpl_DNSd_R(:,:,:)   &
     & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)
  Real(dp), Intent(out) :: gPhoton

  Integer :: i2, i3, n_char, n_Spm, i_gen,i_count
  Real(dp) :: mj2, mi2, m12, m22
  Complex(dp) :: Gcoup(2), Iinte, Jinte, Kinte, coup1, coup2, I2inte

  n_char = Size( mC )
  n_Spm = Size( mSpm2 )
  mj2 = mN(i_in)**2
!  factor(1) = - gSU2 * Sqrt(sW2) * oo8pi2
!  factor(2) = 0.25_dp * factor(1)
!  factor(3) = 0.5_dp / (Pi * Abs(mN(i_in))**3 ) 

  mi2 = mN(i_out)**2
  Gcoup = 0._dp
  !--------------------
  ! W contribution
  !--------------------
  Do i2=1,n_char
   m12 = mC(i2)**2
   Iinte = Igamma(mj2,mi2,mW2,m12)
   I2inte = I2gamma(mj2,mi2,mW2,m12)
   Jinte = Jgamma(mj2,mi2,mW2,m12)
   Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + mW2 * Jinte ) / (mi2 - mj2)
   coup1 = Cpl_CNW_R(i2,i_in) * Cpl_CNW_R(i2,i_out)              &
       & - Conjg( Cpl_CNW_L(i2,i_in) * Cpl_CNW_L(i2,i_out) ) 
   coup2 = Cpl_CNW_L(i2,i_in) * Cpl_CNW_R(i2,i_out)              &
       & - Conjg( Cpl_CNW_R(i2,i_in) * Cpl_CNW_L(i2,i_out) ) 
   Gcoup(1) = Gcoup(1)                                    &
       &    + coup1 * mN(i_in) * (I2inte - Jinte - Kinte) &
       &    - Conjg(coup1) * mN(i_out) * (Kinte - Jinte)  &
       &    + 2._dp * coup2 * mC(i2) * Jinte   
!   Write(10,*) "W chargino",i2
!   Write(10,*) "I ",Iinte
!   Write(10,*) "I2",I2inte
!   Write(10,*) "J ",Jinte
!   Write(10,*) "K ",Kinte
!   Write(10,*) "c1",coup1
!   Write(10,*) "c2",coup2
  End Do    
  !--------------------
  ! S+ contribution
  !--------------------
i_count=3
  Do i2=1,n_char
   m12 = mC(i2)**2
   Do i3=1,n_Spm
    m22 = mSpm2(i3)
    Iinte = Igamma(mj2,mi2,m22,m12)
    I2inte = I2gamma(mj2,mi2,m22,m12)
    Jinte = Jgamma(mj2,mi2,m22,m12)
    Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte ) / (mi2 - mj2)
    coup1 = Cpl_SmpCN_R(i3,i2,i_in) * Conjg( Cpl_SmpCN_R(i3,i2,i_out) )  &
        & - Conjg( Cpl_SmpCN_L(i3,i2,i_in) ) * Cpl_SmpCN_L(i3,i2,i_out)
    coup2 = Cpl_SmpCN_L(i3,i2,i_in) * Conjg( Cpl_SmpCN_R(i3,i2,i_out) )  &
        & - Conjg( Cpl_SmpCN_R(i3,i2,i_in) ) * Cpl_SmpCN_L(i3,i2,i_out) 
    Gcoup(2) = Gcoup(2)                                    &
        &    + coup1 * mN(i_in) * (I2inte - Kinte)         &
        &    - Conjg(coup1) * mN(i_out) * Kinte            &
        &    + coup2 * mC(i2) * Iinte
   End Do    
  End Do    
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
     coup1 = cpl_UNSu_R(i2,i_in,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UNSu_L(i2,i_in,i3) ) * cpl_UNSu_L(i2,i_out,i3)
     coup2 = cpl_UNSu_L(i2,i_in,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UNSu_R(i2,i_in,i3) ) * cpl_UNSu_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                                   &
         & - 2._dp * ( coup1 * mN(i_in) * (I2inte - Kinte) &
         &           - Conjg(coup1) * mN(i_out) * Kinte    &
         &           + coup2 * mf_u(i2) * Iinte )  
    End Do    

   Else
    i_gen = 2*i2-1
    Do i3=i_gen,i_gen+1
     m22 = mUsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_UNSu_R(i2,i_in,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
          & - Conjg( cpl_UNSu_L(i2,i_in,i3) ) * cpl_UNSu_L(i2,i_out,i3)
     coup2 = cpl_UNSu_L(i2,i_in,i3) * Conjg( cpl_UNSu_R(i2,i_out,i3) )  &
         & - Conjg( cpl_UNSu_R(i2,i_in,i3) ) * cpl_UNSu_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                                   &
         & - 2._dp * ( coup1 * mN(i_in) * (I2inte - Kinte) &
         &           - Conjg(coup1) * mN(i_out) * Kinte    &
         &           + coup2 * mf_u(i2) * Iinte )  
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
     coup1 = cpl_DNSd_R(i2,i_in,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DNSd_L(i2,i_in,i3) ) * cpl_DNSd_L(i2,i_out,i3)
     coup2 = cpl_DNSd_L(i2,i_in,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DNSd_R(i2,i_in,i3) ) * cpl_DNSd_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                                   &
         & + ( coup1 * mN(i_in) * (I2inte - Kinte) &
         &           - Conjg(coup1) * mN(i_out) * Kinte    &
         &           + coup2 * mf_d(i2) * Iinte ) 
    End Do    

   Else
    i_gen = 2*i2-1
    Do i3=i_gen,i_gen+1
     m22 = mDsquark2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_DNSd_R(i2,i_in,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DNSd_L(i2,i_in,i3) ) * cpl_DNSd_L(i2,i_out,i3)
     coup2 = cpl_DNSd_L(i2,i_in,i3) * Conjg( cpl_DNSd_R(i2,i_out,i3) )   &
         & - Conjg( cpl_DNSd_R(i2,i_in,i3) ) * cpl_DNSd_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                                   &
         & + ( coup1 * mN(i_in) * (I2inte - Kinte) &
         &           - Conjg(coup1) * mN(i_out) * Kinte    &
         &           + coup2 * mf_d(i2) * Iinte )  
    End Do    
   End If
  End Do    
  !----------------
  ! Slepton Lepton
  !----------------
  Do i2=1,5-n_char
   m12 = mf_l(i2)**2
   If (GenerationMixing) Then
    Do i3=1,2*(5-n_char) 
     m22 = mSlepton2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_LNSl_R(i2,i_in,i3) * Conjg( cpl_LNSl_R(i2,i_out,i3) )  &
          & - Conjg( cpl_LNSl_L(i2,i_in,i3) ) * cpl_LNSl_L(i2,i_out,i3)
     coup2 = cpl_LNSl_L(i2,i_in,i3) * Conjg( cpl_LNSl_R(i2,i_out,i3) )  &
          & - Conjg( cpl_LNSl_R(i2,i_in,i3) ) * cpl_LNSl_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                              &
         &    + coup1 * mN(i_in) * (I2inte - Kinte)   &
         &    - Conjg(coup1) * mN(i_out) * Kinte      &
         &    + coup2 * mf_l(i2) * Iinte
    End Do    

   Else
    i_gen = 2*i2-1
    Do i3=i_gen,i_gen+1
     m22 = mSlepton2(i3)
     Iinte = Igamma(mj2,mi2,m22,m12)
     I2inte = I2gamma(mj2,mi2,m22,m12)
     Jinte = Jgamma(mj2,mi2,m22,m12)
     Kinte = (1._dp + m12 * Iinte - mj2 * I2inte + m22 * Jinte) / (mi2 - mj2)
     coup1 = cpl_LNSl_R(i2,i_in,i3) * Conjg( cpl_LNSl_R(i2,i_out,i3) )  &
          & - Conjg( cpl_LNSl_L(i2,i_in,i3) ) * cpl_LNSl_L(i2,i_out,i3)
     coup2 = cpl_LNSl_L(i2,i_in,i3) * Conjg( cpl_LNSl_R(i2,i_out,i3) )  &
          & - Conjg( cpl_LNSl_R(i2,i_in,i3) ) * cpl_LNSl_L(i2,i_out,i3) 
     Gcoup(2) = Gcoup(2)                                &
         &    + coup1 * mN(i_in) * (I2inte - Kinte)     &
         &    - Conjg(coup1) * mN(i_out) * Kinte        &
         &    + coup2 * mf_l(i2) * Iinte
    End Do    
   End If
  End Do    

  gPhoton = factor(3) * (mj2 - mi2)**3 &
          &           * Abs(factor(1)*Gcoup(1)+factor(2)*Gcoup(2))**2   

 End Subroutine Chi0ToChi0Photon



 Subroutine Chi0ToChimffp(i_in, i_out, mN, mC, n_f, mf, mfp, mW, gW          &
    & , Cpl_CNW_L, Cpl_CNW_R, Cpl_FpFW, mSpm, gSpm, Cpl_SmpCN_L, Cpl_SmpCN_R &
    & , Cpl_SmpFFp_L, Cpl_SmpFFp_R, mSf, gSf, cpl_FNSf_L, cpl_FNSf_R         &
    & , cpl_CFpSf_L, cpl_CFpSf_R, mSfp, gSfp, cpl_FpNSfp_L, cpl_FpNSfp_R     &
    & , cpl_CFSfp_L, cpl_CFSfp_R                                             &
    & , IntegralsW4, n_W4, IntegralsSpm4, n_Spm4, IntegralsSf4, n_Sf4        &
    & , IntegralsWSpm4, n_WSpm, IntegralsWSf8, n_WSf8, IntegralsSpmC4        &
    & , n_SpmC4, IntegralsSpmSf8, n_SpmSf8, IntegralsSfC4, n_SfC4            &
    & , IntegralsSf8, n_Sf8, deltaM, epsI, GenerationMixing, check, fac      &
    & , gCffp, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a neutralino to a Chargino + fermion pair
 ! Written by Werner Porod, 26.06.2001
 ! 12.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_W4, n_Spm4, n_Sf4, n_WSpm, n_WSf8, n_SpmC4 &
     & , n_SpmSf8, n_SfC4, n_Sf8
  Logical, Intent(in) :: check
  
  Real(dp), Intent(in) :: mN(:), mC(:), mf(:), mfp(:), mW, gW, mSpm(:)     &
      & , gSpm(:), mSf(:), gSf(:), mSfp(:), gSfp(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsW4(:,:), IntegralsSpm4(:,:)         &
      & , IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: Cpl_CNW_L(:,:), Cpl_CNW_R(:,:), Cpl_FpFW(:,:)  &
      & , Cpl_SmpCN_L(:,:,:), Cpl_SmpCN_R(:,:,:), Cpl_SmpFFp_L(:,:,:)       &
      & , Cpl_SmpFFp_R(:,:,:), cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:)         &
      & , cpl_CFpSf_L(:,:,:), cpl_CFpSf_R(:,:,:), cpl_FpNSfp_L(:,:,:)       &
      & , cpl_FpNSfp_R(:,:,:), cpl_CFSfp_L(:,:,:), cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsWSpm4(:,:), IntegralsWSf8(:,:)   &
      & , IntegralsSpmC4(:,:), IntegralsSpmSf8(:,:), IntegralsSfC4(:,:)   &
      & , IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gCffp(:,:)

  Integer, Intent(in) :: n_out, WriteContribution

  Integer :: Isum, i1, i2, i3, i4, n_Spm, n_Sfp
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Real(dp), Allocatable :: gCffpSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0ToChimffp'

  If (n_f.Le.0) Then ! useful for the case of R-parity violation
   gCffp = 0._dp
   Iname = Iname - 1
   Return
  End If

  mass(1) = mN(i_in)

  n_Spm = Size( mSpm )
  n_Sfp = Size( mSfp )

  !-------------------------------------------
  ! to cover the case of sneutrinos correctly
  !-------------------------------------------
  If (n_Sfp.Le.3) Then
   n_Sfp = n_f
  Else
   n_Sfp = 2 * n_f
  End If

  Isum = (1 + n_Spm + 4 * n_f)**2
  Allocate( gCffpSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )

  gCffpSum = 0._dp
  Contribution = ' '

  !-----
  ! W W
  !-----
  Isum = 1
  Boson2(1) = mW
  Boson2(2) = gW

  coup1(1) = Cpl_CNW_L(i_out,i_in)
  coup1(2) = Cpl_CNW_R(i_out,i_in)
  coup1(4) = 0._dp

  If (GenerationMixing) Then
   Do i2=1,n_f
    Do i3=1,n_f
     mass(2) = mC(i_out)
     mass(3) = -mf(i2)
     coup1(3) = cpl_FpFW(i3,i2) 
     mass(4) = mfp(i3)
     Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                          &, n_W4, resR, check)
     gCffpSum(i2, i3, Isum) = resR
     Contribution(i2,i3,1) = 'W f_'//Bu(i2)//' fp_'//Bu(i3)
    End Do
   End Do

  Else
   Do i2=1,n_f
    mass(2) = mC(i_out)
    mass(3) = -mf(i2)
    coup1(3) = cpl_FpFW(i2,i2) 
    mass(4) = mfp(i2)
    Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                         &, n_W4, resR, check)
    gCffpSum(i2, i2, Isum) = resR
    Contribution(i2,i2,1) = 'W f_'//Bu(i2)//' fp_'//Bu(i2)
   End Do
  End If

  !-------------------
  ! S^- S^-, diagonal
  !-------------------
  Do i1=1,n_Spm
   Isum = Isum + 1

   coup1(1) = Cpl_SmpCN_L(i1, i_out,i_in)
   coup1(2) = Cpl_SmpCN_R(i1, i_out,i_in) 
   Boson2(1) = mSpm(i1)
   Boson2(2) = gSpm(i1)

   If (GenerationMixing) Then
    Do i2=1,n_f
     Do i3=1,n_f
      mass(2) = mC(i_out)
      mass(3) = -mf(i2)
      mass(4) = mfp(i3)
      coup1(3) = Conjg(cpl_SmpFFp_L(i1,i2,i3) )
      coup1(4) = Conjg(cpl_SmpFFp_R(i1,i2,i3) )
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSpm4 ,n_Spm4, resR, check)
      gCffpSum(i2,i3,Isum) = resR
      Contribution(i2,i3,Isum) = 'S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i3)
     End Do
    End Do

   Else
    Do i2=1,n_f
     mass(2) = mC(i_out)
     mass(3) = -mf(i2)
     mass(4) = mfp(i2)
     coup1(3) = Conjg(cpl_SmpFFp_L(i1,i2,i2) )
     coup1(4) = Conjg(cpl_SmpFFp_R(i1,i2,i2) )
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSpm4 ,n_Spm4, resR, check)
     gCffpSum(i2,i2,Isum) = resR
     Contribution(i2,i2,Isum) = 'S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i2)
    End Do

   End If
  End Do

  !-------------------
  ! Sfp Sfp, diagonal
  !-------------------
  
  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i1=1,n_f
     coup1(1) = cpl_FpNSfp_L(i1,i_in,i2)
     coup1(2) = cpl_FpNSfp_R(i1,i_in,i2)
     Do i3=1,n_f
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

   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)

    If (n_Sfp.Le.3) Then
     i1 = i2
    Else
     i1 = (i2+1)/2
    End If

    mass(2) = mfp(i1)
    mass(3) = -mf(i1)
    mass(4) = mC(i_out)
    coup1(1) = cpl_FpNSfp_L(i1,i_in,i2)
    coup1(2) = cpl_FpNSfp_R(i1,i_in,i2)
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
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    Do i1=1,n_f
     coup1(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
     coup1(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
     Do i3=1,n_f
      mass(2) = mf(i1)
      mass(3) = - mC(i_out)
      mass(4) = mfp(i3)
      coup1(3) = cpl_CFpSf_L(i_out,i3,i2)
      coup1(4) = cpl_CFpSf_R(i_out,i3,i2)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gCffpSum(i1,i3,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do
  Else

   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mf(i1)
    mass(3) = - mC(i_out)
    mass(4) = mfp(i1)
    coup1(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
    coup1(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
    coup1(3) = cpl_CFpSf_L(i_out,i1,i2)
    coup1(4) = cpl_CFpSf_R(i_out,i1,i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gCffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do
  End If

  !-------------------------------
  ! W boson - S+ boson  s-channel
  !-------------------------------
  Boson4(1) = mW
  Boson4(2) = gW
  coup2(1) = Cpl_CNW_L(i_out,i_in)
  coup2(2) = Cpl_CNW_R(i_out,i_in)
  coup2(6) = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_Spm 
    Isum = Isum + 1
    Boson4(3) = mSpm(i1)
    Boson4(4) = gSpm(i1)

    coup2(3) = Conjg(Cpl_SmpCN_R(i1,i_out,i_in))
    coup2(4) = Conjg(Cpl_SmpCN_L(i1,i_out,i_in))
    Do i2=1,n_f
     Do i3=1,n_f
      mass(2) = mC(i_out)
      mass(3) = -mf(i2)
      mass(4) = mfp(i3)
      coup2(5) = cpl_FpFW(i3,i2)
      coup2(7) = cpl_SmpFFp_L(i1,i2,i3)
      coup2(8) = cpl_SmpFFp_R(i1,i2,i3)
      Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsWSpm4, n_WSpm, resC, check)
      gCffpSum(i2,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i2,i3,Isum) = &
        & 'W S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i1=1,n_Spm 
    Isum = Isum + 1
    Boson4(3) = mSpm(i1)
    Boson4(4) = gSpm(i1)

    coup2(3) = Conjg(Cpl_SmpCN_R(i1,i_out,i_in))
    coup2(4) = Conjg(Cpl_SmpCN_L(i1,i_out,i_in))
    Do i2=1,n_f
     mass(2) = mC(i_out)
     mass(3) = -mf(i2)
     mass(4) = mfp(i2)
     coup2(5) = cpl_FpFW(i2,i2)
     coup2(7) = cpl_SmpFFp_L(i1,i2,i2)
     coup2(8) = cpl_SmpFFp_R(i1,i2,i2)
     Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsWSpm4, n_WSpm, resC, check)
     gCffpSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
     Contribution(i2,i2,Isum) = &
        & 'W S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i2)
    End Do
   End Do

  End If

  !---------------------------------
  ! W boson - Sfp_{xyz}  t-channel
  !---------------------------------
  coup2(1) = Cpl_CNW_L(i_out, i_in)
  coup2(2) = Cpl_CNW_R(i_out, i_in)
  coup2(6) = 0._dp 

  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)
    Do i1=1,n_f
     coup2(3) = Conjg( cpl_FpNSfp_R(i1, i_in, i2) )
     coup2(4) = Conjg( cpl_FpNSfp_L(i1, i_in, i2) )

     Do i3=1,n_f
      mass(2) = mfp(i1)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)

      coup2(5) = cpl_FpFW(i1,i3) 
      coup2(7) = cpl_CFSfp_L(i_out, i3, i2)
      coup2(8) = cpl_CFSfp_R(i_out, i3, i2)
      Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                 &, IntegralsWSf8, n_WSf8, resC, check)
      gCffpSum(i3,i1,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i3,i1,Isum) = &
        & 'W Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do

  Else
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)

    If (n_Sfp.Le.3) Then
     i1 = i2
    Else
     i1 = (i2+1)/2
    End If

    mass(2) = mfp(i1)
    mass(3) = -mf(i1)
    mass(4) = mC(i_out)

    coup2(3) = Conjg( cpl_FpNSfp_R(i1, i_in, i2) )
    coup2(4) = Conjg( cpl_FpNSfp_L(i1, i_in, i2) )
    coup2(5) = cpl_FpFW(i1,i1) 
    coup2(7) = cpl_CFSfp_L(i_out, i1, i2)
    coup2(8) = cpl_CFSfp_R(i_out, i1, i2)
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsWSf8, n_WSf8, resC, check)
    gCffpSum(i1,i1,Isum) = - 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
      & 'W Sfp_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !---------------------------------
  ! W boson - Sf_{xyz}  u-channel
  !---------------------------------
  coup2(1) = Cpl_CNW_L(i_out, i_in)
  coup2(2) = Cpl_CNW_R(i_out, i_in)
  coup2(5) = 0._dp 

  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    Do i1=1,n_f
     coup2(3) = cpl_FNSf_L(i1, i_in, i2)
     coup2(4) = cpl_FNSf_R(i1, i_in, i2)

     Do i3=1,n_f
      mass(2) = mf(i1)
      mass(3) = -mfp(i3)
      mass(4) = mC(i_out)

      coup2(6) = cpl_FpFW(i3,i1) 
      coup2(7) = Conjg( cpl_CFpSf_R(i_out, i3, i2) )
      coup2(8) = Conjg( cpl_CFpSf_L(i_out, i3, i2) )
      Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsWsf8, n_WSf8, resC, check)
      gCffpSum(i1,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i1,i3,Isum) = &
        & 'W Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mf(i1)
    mass(3) = -mfp(i1)
    mass(4) = mC(i_out)

    coup2(3) = cpl_FNSf_L(i1, i_in, i2)
    coup2(4) = cpl_FNSf_R(i1, i_in, i2)
    coup2(6) = cpl_FpFW(i1,i1) 
    coup2(7) = Conjg( cpl_CFpSf_R(i_out, i1, i2) )
    coup2(8) = Conjg( cpl_CFpSf_L(i_out, i1, i2) )
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsWsf8, n_WSf8, resC, check)
    gCffpSum(i1,i1,Isum) = 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
        & 'W Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !-------------------------------------------
  ! S+ boson  s-channel - S+ boson  s-channel
  !-------------------------------------------
  Do i1 = 1,n_Spm
   Boson4(1) = mSpm(i1)
   Boson4(2) = gSpm(i1)
   coup2(1) = Cpl_SmpCN_L(i1,i_out,i_in)
   coup2(2) = Cpl_SmpCN_R(i1,i_out,i_in)
   Do i2 = i1+1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(i2)
    Boson4(4) = gSpm(i2)
    coup2(3) = Conjg(Cpl_SmpCN_R(i2,i_out,i_in))
    coup2(4) = Conjg(Cpl_SmpCN_L(i2,i_out,i_in))

    If (GenerationMixing) Then
     Do i3=1,n_f
      Do i4=1,n_f
       mass(2) = mC(i_out)
       mass(3) = -mf(i3)
       mass(4) = mfp(i4)

       coup2(5) = Conjg(cpl_SmpFFp_R(i1,i3,i4) )
       coup2(6) = Conjg(cpl_SmpFFp_L(i1,i3,i4) )
       coup2(7) = cpl_SmpFFp_L(i2,i3,i4)
       coup2(8) = cpl_SmpFFp_R(i2,i3,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSpmC4, n_SpmC4, resC, check)
       gCFFpSum(i3,i4,Isum) = 2._dp * Real(resC,dp)
       Contribution(i3,i4,Isum) = &
        & 'S^-_'//Bu(i1)//' S^-_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do

    Else
     Do i3=1,n_f
      mass(2) = mC(i_out)
      mass(3) = -mf(i3)
      mass(4) = mfp(i3)
      coup2(5) = Conjg(cpl_SmpFFp_R(i1,i3,i3) )
      coup2(6) = Conjg(cpl_SmpFFp_L(i1,i3,i3) )
      coup2(7) = cpl_SmpFFp_L(i2,i3,i3)
      coup2(8) = cpl_SmpFFp_R(i2,i3,i3)
      Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsSpmC4, n_SpmC4, resC, check)
      gCFFpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i3,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' S^-_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
     End Do

    End If
   End Do
  End Do

  !---------------------------------------------
  ! S+ boson  s-channel - Sfp_{xyz}  t-channel
  !---------------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Cpl_SmpCN_L(i1,i_out,i_in)
    coup2(2) = Cpl_SmpCN_R(i1,i_out,i_in)
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     Do i3 =1,n_f
      coup2(3) = Conjg(cpl_FpNSfp_R(i3,i_in, i2))
      coup2(4) = Conjg(cpl_FpNSfp_L(i3,i_in, i2))
      Do i4=1,n_f
       mass(2) = mfp(i3)
       mass(3) = -mf(i4)
       mass(4) = mC(i_out)
       coup2(5) = Conjg(cpl_SmpFFp_R(i1,i4,i3) )
       coup2(6) = Conjg(cpl_SmpFFp_L(i1,i4,i3) )
       coup2(7) = cpl_CFSfp_L(i_out,i4,i2)
       coup2(8) = cpl_CFSfp_R(i_out,i4,i2)
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSpmSf8, n_SpmSf8, resC, check)
       gCffpSum(i4,i3,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Cpl_SmpCN_L(i1,i_out,i_in)
    coup2(2) = Cpl_SmpCN_R(i1,i_out,i_in)
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     If (n_Sfp.Le.3) Then
      i3 = i2
     Else
      i3 = (i2+1)/2
     End If

     mass(2) = mfp(i3)
     mass(3) = -mf(i3)
     mass(4) = mC(i_out)

     coup2(3) = Conjg(cpl_FpNSfp_R(i3,i_in, i2))
     coup2(4) = Conjg(cpl_FpNSfp_L(i3,i_in, i2))
     coup2(5) = Conjg(cpl_SmpFFp_R(i1,i3,i3) )
     coup2(6) = Conjg(cpl_SmpFFp_L(i1,i3,i3) )
     coup2(7) = cpl_CFSfp_L(i_out,i3,i2)
     coup2(8) = cpl_CFSfp_R(i_out,i3,i2)
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSpmSf8, n_SpmSf8, resC, check)
     gCffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !---------------------------------------------
  ! S+ boson  s-channel - Sf_{xyz}  u-channel
  !---------------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Cpl_SmpCN_L(i1,i_out,i_in)
    coup2(2) = Cpl_SmpCN_R(i1,i_out,i_in)
    Do i2 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 =1,n_f
      coup2(3) = cpl_FNSf_L(i3, i_in, i2)
      coup2(4) = cpl_FNSf_R(i3, i_in, i2)
      Do i4=1,n_f
       mass(2) = mf(i3)
       mass(3) = -mfp(i4)
       mass(4) = mC(i_out)

       coup2(5) = Conjg(cpl_SmpFFp_R(i1,i3,i4))
       coup2(6) = Conjg(cpl_SmpFFp_L(i1,i3,i4))
       coup2(7) = Conjg(cpl_CFpSf_R(i_out,i4,i2))
       coup2(8) = Conjg(cpl_CFpSf_L(i_out,i4,i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSpmSf8, n_SpmSf8, resC, check)
       gCffpSum(i3,i4,Isum) = - 2._dp * Real(resC,dp)      
       Contribution(i3,i4,Isum) = &
        & 'S^-_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Cpl_SmpCN_L(i1,i_out,i_in)
    coup2(2) = Cpl_SmpCN_R(i1,i_out,i_in)
    Do i2 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     i3 = (i2+1)/2

     mass(2) = mf(i3)
     mass(3) = -mfp(i3)
     mass(4) = mC(i_out)

     coup2(3) = cpl_FNSf_L(i3, i_in, i2)
     coup2(4) = cpl_FNSf_R(i3, i_in, i2)
     coup2(5) = Conjg(cpl_SmpFFp_R(i1,i3,i3))
     coup2(6) = Conjg(cpl_SmpFFp_L(i1,i3,i3))
     coup2(7) = Conjg(cpl_CFpSf_R(i_out,i3,i2))
     coup2(8) = Conjg(cpl_CFpSf_L(i_out,i3,i2))
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSpmSf8, n_SpmSf8, resC, check)
     gCffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)      
     Contribution(i3,i3,Isum) = &
      & 'S^-_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !-----------------------------------------------
  ! Sfp_{xyz}  t-channel - Sfp_{xyz}  t-channel
  !-----------------------------------------------
  If (GenerationMixing) Then
   Do i1=1,n_Sfp-1
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2=i1+1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     Do i3 = 1,n_f
      coup2(1) = cpl_FpNSfp_L(i3,i_in,i1)
      coup2(2) = cpl_FpNSfp_R(i3,i_in,i1)
      coup2(3) = Conjg(cpl_FpNSfp_R(i3,i_in,i2))
      coup2(4) = Conjg(cpl_FpNSfp_L(i3,i_in,i2))
      Do i4=1,n_f
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

  Else If (n_Sfp.Gt.3) Then ! no  contribution in case of sneutrinos
   Do i1=1,n_Sfp-1,2
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

    coup2(1) = cpl_FpNSfp_L(i3,i_in,i1)
    coup2(2) = cpl_FpNSfp_R(i3,i_in,i1)
    coup2(3) = Conjg(cpl_FpNSfp_R(i3,i_in,i2))
    coup2(4) = Conjg(cpl_FpNSfp_L(i3,i_in,i2))
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
   Do i1= 1,n_Sfp
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2 =1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 = 1,n_f
      coup2(1) = cpl_FpNSfp_L(i3,i_in,i1)
      coup2(2) = cpl_FpNSfp_R(i3,i_in,i1)
      Do i4 =1, n_f
       mass(2) = mf(i4)
       mass(3) = -mC(i_out)
       mass(4) = mfp(i3)
       coup2(3) = cpl_FNSf_L(i4, i_in, i2)
       coup2(4) = cpl_FNSf_R(i4, i_in, i2)
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
   Do i1= 1,n_Sfp
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    If (n_Sfp.Le.3) Then
     i3 = i1
    Else
     i3 = (i1+1)/2
    End If
    coup2(1) = cpl_FpNSfp_L(i3,i_in,i1)
    coup2(2) = cpl_FpNSfp_R(i3,i_in,i1)

    Do i2 =1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     mass(2) = mf(i3)
     mass(3) = -mC(i_out)
     mass(4) = mfp(i3)

     coup2(3) = cpl_FNSf_L(i3, i_in, i2)
     coup2(4) = cpl_FNSf_R(i3, i_in, i2)
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

   Do i1=1,2*n_f-1
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Do i2=i1+1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i3 = 1,n_f
      coup2(1) = Conjg(cpl_FNSf_R(i3, i_in, i1))
      coup2(2) = Conjg(cpl_FNSf_L(i3, i_in, i1))
      coup2(3) = cpl_FNSf_L(i3, i_in, i2)
      coup2(4) = cpl_FNSf_R(i3, i_in, i2)
      Do i4=1,n_f
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
   Do i1=1,2*n_f-1,2
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

    coup2(1) = Conjg(cpl_FNSf_R(i3, i_in, i1))
    coup2(2) = Conjg(cpl_FNSf_L(i3, i_in, i1))
    coup2(3) = cpl_FNSf_L(i3, i_in, i2)
    coup2(4) = cpl_FNSf_R(i3, i_in, i2)
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
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gCffp(i1,i2) = Sum( gCffpSum(i1,i2,1:Isum) )
     If (gCffp(i1,i2).Lt.0._dp) Then
      p_test = Abs(gCffp(i1,i2)) / Maxval(Abs(gCffpSum(i1,i2,1:Isum)))
      gCffp(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero
      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//') < 0 :' &
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
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//') < 0 :' &
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
  If (WriteContribution.ne.0) Then

   gCffpSum = gCffpSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//') :' &
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
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//') :' &
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

  Deallocate( gCffpSum, Contribution)

  Iname = Iname - 1

 End Subroutine Chi0ToChimffp


 Subroutine Chi0toGqq(i_in, state, mN, mGlu, mf, mSf, gSf, cpl_FNSf_L      &
    & , cpl_FNSf_R, cpl_FGSf_L, cpl_FGSf_R, IntegralsSf4                   &
    & , n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI    &
    & , GenerationMixing, check, fac, gGff, WriteContribution, n_out)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a neutralino into gluino + quarks
 ! input:
 !    n_in ....... index of the decaying neutralino, in case n_in < 0,
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !
 ! output:
 ! written by Werner Porod, 27.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !------------------------------------------------------------------
 Implicit None
  Character(len=5), Intent(in) :: state 
  Integer, Intent(in) :: i_in
  Integer, Intent(inout) :: n_Sf4, n_CSf4, n_Sf8
  Logical, Intent(in) :: check

  Real(dp), Intent(in) :: mN(:), mGlu, mf(:), mSf(:), gSf(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:) &
                         & , cpl_FGSf_L(:,:), cpl_FGSf_R(:,:)
  Complex(dp), Intent(inout) :: IntegralsCSf4(:,:), IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gGff(:,:)

  Integer, Intent(in) :: n_out, WriteContribution

  Integer :: Isum, i1, i2, i3, i4, i_gen
  Real(dp) :: Boson2(2), mass(4), resR, resRa, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8), resCa

  Real(dp), Allocatable :: gGffSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Chi0ToGqq'

  mass(1) = mN(i_in)

  Allocate( gGffSum(3, 3, 80) )
  Allocate( Contribution(3, 3, 80) )
   
  gGffSum = 0._dp
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
      coup1(1) = cpl_FNSf_L(i1,i_in,i2)
      coup1(2) = cpl_FNSf_R(i1,i_in,i2)
      coup1(3) = Conjg(cpl_FGSf_R(i3,i2))
      coup1(4) = Conjg(cpl_FGSf_L(i3,i2))
      mass(2) = mf(i1)
      mass(3) = - mf(i3)
      mass(4) = mGlu        
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      !------------------------
      ! u-channel
      !------------------------
      coup1(1) = Conjg(cpl_FNSf_R(i3,i_in,i2))  ! u-channel
      coup1(2) = Conjg(cpl_FNSf_L(i3,i_in,i2))
      coup1(3) = cpl_FGSf_L(i1,i2)
      coup1(4) = cpl_FGSf_R(i1,i2)
      mass(2) = mf(i3)
      mass(3) = - mGlu        
      mass(4) = mf(i1)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resRa, check)
      gGffSum(i1,i3,Isum) = resR + resRa
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
     coup1(1) = cpl_FNSf_L(i1,i_in,i2)
     coup1(2) = cpl_FNSf_R(i1,i_in,i2)
     coup1(3) = Conjg(cpl_FGSf_R(i1,i2))
     coup1(4) = Conjg(cpl_FGSf_L(i1,i2))
     mass(2) = mf(i1)
     mass(3) = - mf(i1)
     mass(4) = mGlu        
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resR, check)
     !------------------------
     ! u-channel
     !------------------------
     coup1(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
     coup1(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
     coup1(3) = cpl_FGSf_L(i1,i2)
     coup1(4) = cpl_FGSf_R(i1,i2)
     mass(2) = mf(i1)
     mass(3) = - mGlu        
     mass(4) = mf(i1)
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resRa, check)
     gGffSum(i1,i1,Isum) = resR + resRa
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
       mass(4) = mGlu
       coup2(1) = cpl_FNSf_L(i1,i_in,i3)
       coup2(2) = cpl_FNSf_R(i1,i_in,i3)
       coup2(3) = Conjg(cpl_FNSf_R(i1,i_in,i4))
       coup2(4) = Conjg(cpl_FNSf_L(i1,i_in,i4))
       coup2(5) = Conjg(cpl_FGSf_R(i2,i3))
       coup2(6) = Conjg(cpl_FGSf_L(i2,i3))
       coup2(7) = cpl_FGSf_L(i2,i4)
       coup2(8) = cpl_FGSf_R(i2,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resC, check)
       !-------------
       ! u-channel
       !-------------
       mass(2) = mf(i2)
       mass(3) = -mGlu
       mass(4) = mf(i1)
       coup2(1) = Conjg(cpl_FNSf_R(i2,i_in,i3))
       coup2(2) = Conjg(cpl_FNSf_L(i2,i_in,i3))
       coup2(3) = cpl_FNSf_L(i2,i_in,i4)
       coup2(4) = cpl_FNSf_R(i2,i_in,i4)
       coup2(5) = cpl_FGSf_L(i1,i3)
       coup2(6) = cpl_FGSf_R(i1,i3)
       coup2(7) = Conjg(cpl_FGSf_R(i1,i4))
       coup2(8) = Conjg(cpl_FGSf_L(i1,i4))
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resCa, check)
       gGffSum(i1,i2,Isum) = 2._dp * Real(resC+resCa,dp)
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
    mass(4) = mGlu
    coup2(1) = cpl_FNSf_L(i1,i_in,i2)
    coup2(2) = cpl_FNSf_R(i1,i_in,i2)
    coup2(3) = Conjg(cpl_FNSf_R(i1,i_in,i3))
    coup2(4) = Conjg(cpl_FNSf_L(i1,i_in,i3))
    coup2(5) = Conjg(cpl_FGSf_R(i1,i2))
    coup2(6) = Conjg(cpl_FGSf_L(i1,i2))
    coup2(7) = cpl_FGSf_L(i1,i3)
    coup2(8) = cpl_FGSf_R(i1,i3)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resC, check)
    !-------------
    ! u-channel
    !-------------
    mass(2) = mf(i1)
    mass(3) = -mGlu
    mass(4) = mf(i1)
    coup2(1) = Conjg(cpl_FNSf_R(i1,i_in,i2))
    coup2(2) = Conjg(cpl_FNSf_L(i1,i_in,i2))
    coup2(3) = cpl_FNSf_L(i1,i_in,i3)
    coup2(4) = cpl_FNSf_R(i1,i_in,i3)
    coup2(5) = cpl_FGSf_L(i1,i2)
    coup2(6) = cpl_FGSf_R(i1,i2)
    coup2(7) = Conjg(cpl_FGSf_R(i1,i3))
    coup2(8) = Conjg(cpl_FGSf_L(i1,i3))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resCa, check)
    gGffSum(i1,i1,Isum) = 2._dp * Real(resC+resCa,dp)
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
      coup2(1) = cpl_FNSf_L(i1,i_in,i2)
      coup2(2) = cpl_FNSf_R(i1,i_in,i2)
      coup2(7) = Conjg(cpl_FGSf_R(i1,i3))
      coup2(8) = Conjg(cpl_FGSf_L(i1,i3))
      Do i4=1,3
       mass(2) = mf(i4)
       mass(3) = -mGlu
       mass(4) = mf(i1)
       coup2(3) = cpl_FNSf_L(i4,i_in, i3)
       coup2(4) = cpl_FNSf_R(i4,i_in, i3)
       coup2(5) = Conjg(cpl_FGSf_R(i4,i2))
       coup2(6) = Conjg(cpl_FGSf_L(i4,i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSf8, n_Sf8, resC, check)
       gGffSum(i1,i4,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i1,i4,Isum) = &
          &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i4)
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
     coup2(1) = cpl_FNSf_L(i1,i_in,i2)
     coup2(2) = cpl_FNSf_R(i1,i_in,i2)
     coup2(5) = Conjg(cpl_FGSf_R(i1,i2))
     coup2(6) = Conjg(cpl_FGSf_L(i1,i2))
     Do i3 = i_gen, i_gen+1
      Isum = Isum + 1
      Boson4(3) = mSf(i3)
      Boson4(4) = gSf(i3)
      mass(2) = mf(i1)
      mass(3) = -mGlu
      mass(4) = mf(i1)
      coup2(3) = cpl_FNSf_L(i1,i_in, i3)
      coup2(4) = cpl_FNSf_R(i1,i_in, i3)
      coup2(7) = Conjg(cpl_FGSf_R(i1,i3))
      coup2(8) = Conjg(cpl_FGSf_L(i1,i3))
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSf8, n_Sf8, resC, check)
      gGffSum(i1,i1,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i1,i1,Isum) = &
         &  'tu Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
     End Do
    End Do
   End Do

  End If

  !----------
  ! Summing
  !----------
  gGff = 0._dp
  Do i1=1,3
   If (GenerationMixing) Then
    Do i2=1,3
     gGff(i1,i2) = Sum( gGffSum(i1,i2,1:Isum) )
     If (gGff(i1,i2).Lt.0._dp) Then
      p_test = Abs(gGff(i1,i2)) / Maxval(Abs(gGffSum(i1,i2,1:Isum)))
      gGff(i1,i2) = 0._dp
      If (p_test.le.prec) cycle ! this is a numerical zero
      Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'Gamma(Chi_'//Bu(i_in)//' -> g '//state//') < 0 :' &
      & ,i1,i2,gGff(i1,i2)
      Write(ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gGffSum(i1,i2,i3).Ne.0._dp) &
        &      Write(ErrCan,*) Contribution(i1,i2,i3),gGffSum(i1,i2,i3)
      End Do
     End If
    End Do

   Else
    gGff(i1,i1) = Sum( gGffSum(i1,i1,1:Isum) )
    If (gGff(i1,i1).Lt.0._dp) Then
     p_test = Abs(gGff(i1,i1)) / Maxval(Abs(gGffSum(i1,i1,1:Isum)))
     gGff(i1,i1) = 0._dp
     If (p_test.le.prec) cycle ! this is a numerical zero
     Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) 'Gamma(Chi_'//Bu(i_in)//' -> g '//state//') < 0 :' &
      & ,i1,i1,gGff(i1,i1)
     Write(ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gGffSum(i1,i1,i3).Ne.0._dp) &
        & Write(ErrCan,*) Contribution(i1,i1,i3),gGffSum(i1,i1,i3)
     End Do
    End If
   End If
  End Do

  gGff = fac * gGff

  !---------------------------
  ! for detailed information
  !---------------------------
  If (WriteContribution.ne.0) Then

   gGffSum = gGffSum * fac

   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
     Write (n_out,*)  'Gamma(Chi_'//Bu(i_in)//' -> g '//state//') :' &
      & ,i1,i2,gGff(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gGffSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gGffSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,3
     Write (n_out,*)  'Gamma(Chi_'//Bu(i_in)//' -> g '//state//') :' &
      & ,i1,i1,gGff(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gGffSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gGffSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Deallocate( gGffSum, Contribution)

  Iname = Iname - 1

 End Subroutine Chi0toGqq


End Module Neut3Decays

