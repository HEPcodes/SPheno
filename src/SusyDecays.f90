Module SusyDecays
! comments
! In this module the routines for the decays of SUSY particles are
! stored. 

! load modules
Use Control
Use DecayFunctions
Use Mathematics, Only:  Li2
Use LoopCouplings
! load modules

Contains


 Subroutine ChargedscalarTwoBodyDecays(i_in, n_s0, n_nu, id_nu, n_l, id_l, n_d &
      & , id_d, n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n &
      & , n_c, n_p0, n_Spm, Spm, mf_l, cpl_SmpLNu_L, cpl_SmpLNu_R, mf_d, mf_u  &
      & , cpl_SmpDU_L, cpl_SmpDU_R, Slept, Sneut, cpl_SmpSlSn, Sdown, Sup      &
      & , cpl_SmpSdSu, Chi0, ChiPm, cpl_SmpCN_L, cpl_SmpCN_R, mW, mZ, cpl_SmpZ &
      & , P0, cpl_SmpP03, cpl_SmpP0W, S0, cpl_SmpS03, cpl_SmpS0W, k_neut)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of pseudoscalars
 ! input:
 !  i_in ................. specifies the decaying pseudoscalar. The decays of
 !                         all pseudoscalars are calculated if i_in < 0.
 !  mSpm(i) .............. masses of charged scalars
 !  mf_l(i) .............. lepton masses
 !  cpl_SmpLNu_L(i,j,k) .. left coupling charged scalar - lepton - neutrino
 !  cpl_SmpLNu_R(i,j,k) .. right coupling charged scalar - lepton - neutrino
 !  mf_d(i) ............. d-quark masses
 !  mf_u(i) ............. u-quark masses
 !  cpl_SmpDU_L(i,j,k) .. left coupling charged scalar d-quark u-quark
 !  cpl_SmpDU_R(i,j,k) .. right coupling charged scalar d-quark u-quark
 !  mSlepton(i) ......... slepton masses
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_SmpSlSn(i,j,k) .. coupling charged scalar slepton sneutrino
 !  mSdown(i) ........... d-squark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_SmpSdSu(i,j,k) .. coupling charged scalar d-squark u-squark
 !  mN(i) ............... neutralino masses
 !  mC(i) ............... chargino masses
 !  cpl_SmpCN_L(i,j,k) .. left charged scalar - chargino - neutralino coupling
 !  cpl_SmpCN_R(i,j,k) .. right charged scalar - chargino - neutralino coupling
 !  mW .................. mass of the W-boson
 !  mZ .................. mass of the Z-boson
 !  cpl_SmpZ(i,j) ....... coupling charged scalar - scalar - Z
 !  mP0(i) .............. pseudoscalar masses
 !  cpl_SmpP03(i,j,k) ... charged scalar - charged scalar - pseudo scalar 
 !  cpl_SmpP0W(i) ....... charged scalar - pseudo scalar - W coupling
 !  mS0(i) .............. scalar masses
 !  cpl_SmpS03(i,j,k) ... charged scalar - charged scalar - scalar coupling
 !  cpl_SmpS0W(i,j) ..... charged scalar - scalar - W coupling
 !  k_neut .............. summing over neutrinos if =1
 ! written by Werner Porod, 30.04.2001
 ! 19.09.2010: adjusting to new variable types
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut, n_s0, n_nu, id_nu(:), n_l, id_l(:), n_d  &
      & , id_d(:), n_u, id_u(:), n_Z, id_Z(:), n_W, id_W(:), n_snu, n_sle, n_Sd &
      & , n_su, n_n, n_c, n_p0, n_Spm
  Real(dp), Intent(in) :: mf_l(:), mf_d(:), mf_u(:), mW(:), mZ(:)
  Complex(dp), Intent(in) ::  cpl_SmpLNu_L(:,:,:), cpl_SmpLNu_R(:,:,:)    &
          & , cpl_SmpDU_L(:,:,:), cpl_SmpDU_R(:,:,:), cpl_SmpSlSn(:,:,:)  &
          & , cpl_SmpSdSu(:,:,:), cpl_SmpCN_L(:,:,:), cpl_SmpCN_R(:,:,:)  &
          & , cpl_SmpP03(:,:,:), cpl_SmpP0W(:,:,:), cpl_SmpZ(:,:,:)       &
          & , cpl_SmpS03(:,:,:), cpl_SmpS0W(:,:,:)
  Type(particle2), Intent(in) :: Sdown(:), P0(:)
  Type(particle23), Intent(in) :: Sneut(:), Slept(:), Sup(:), Chi0(:) &
          & , ChiPm(:), S0(:)
  Type(particle2), intent(inout) :: Spm(:)

  Integer :: i1, i2, i_start, i_end, i_count, i3
  Real(dp) :: gam, m_in, mS0(n_S0),  mSlepton(n_sle), mP0(n_P0), mSpm(n_Spm) &
    & , mSneutrino(n_snu), mSdown(n_sd), mSup(n_su), mN(n_n), mC(n_c)
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedscalarTwoBodyDecays'


  If (i_in.Lt.0) Then
   i_start = 2
   i_end = n_Spm
   Do i1=1,n_S0
    Spm(i1)%gi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_Spm) ) Then 
   i_start = i_in 
   i_end = i_in
   Spm(i_in)%gi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_Spm) = ',i_in,n_Spm
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mSPm = Spm%m
  mSlepton = Slept%m
  mSneutrino = Sneut%m
  mSdown = Sdown%m
  mSup = Sup%m
  mN = Chi0%m
  mC = ChiPm%m
  mP0 = P0%m
  mS0 = S0%m

  Do i1 = i_start, i_end
   m_in = mSpm(i1)
   Spm(i1)%gi2 = 0._dp
   i_count = 1
   
   If (m_in.Eq.0._dp) Cycle ! massless particle

   !------------------
   ! into leptons
   !------------------
   Do i2 = 1,n_l
    Do i3 = 1,n_nu
     If ((Abs(cpl_SmpLNu_L(i1,i2,i3))+Abs(cpl_SmpLNu_R(i1,i2,i3))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_l(i2), 0._dp, cpl_SmpLNu_L(i1,i2,i3) &
                             &, cpl_SmpLNu_R(i1,i2,i3), gam )
      If (k_neut.Eq.1) Then
       Spm(i1)%gi2(i_count) = gam
       Spm(i1)%id2(i_count,1) = id_l(i2)+1
       Spm(i1)%id2(i_count,2) = id_nu(1)
      Else
       Spm(i1)%gi2(i_count) = gam
       Spm(i1)%id2(i_count,1) = id_l(i2)+1
       Spm(i1)%id2(i_count,2) = id_nu(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
    If (k_neut.Eq.1) i_count = i_count + 1
   End Do
   !------------------
   ! into quarks
   !------------------
   Do i2 = 1, n_d
    Do i3 = 1,n_u
     If ((Abs(cpl_SmpDU_L(i1,i2,i3))+Abs(cpl_SmpDU_R(i1,i2,i3))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_u(i3), cpl_SmpDU_L(i1,i2,i3) &
                             &, cpl_SmpDU_R(i1,i2,i3), gam )
      gam = 3._dp * gam 
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = id_d(i2) + 1
      Spm(i1)%id2(i_count,2) = id_u(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do
   !------------------
   ! into sleptons
   !------------------
   Do i2 = 1,n_sle
    Do i3 = 1,n_snu
     If (Abs(cpl_SmpSlSn(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSneutrino(i3)   &
                             &, cpl_SmpSlSn(i1,i2,i3), gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = Slept(i2)%id + 1
      Spm(i1)%id2(i_count,2) = Sneut(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !------------------
   ! into squarks
   !------------------
   Do i2 = 1, n_sd
    Do i3 = 1, n_su
     If (Abs(cpl_SmpSdSu(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSup(i3)   &
                             &, cpl_SmpSdSu(i1,i2,i3), gam )
      Spm(i1)%gi2(i_count) = 3._dp * gam
      Spm(i1)%id2(i_count,1) = Sdown(i2)%id + 1
      Spm(i1)%id2(i_count,2) = Sup(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !------------------------
   ! Charginos + Neutralinos
   !------------------------
   Do i2 = 1, n_c
    Do i3 = 1, n_n
     Call ScalarToTwoFermions(m_in, mC(i2), mN(i3), cpl_SmpCN_L(i1,i2,i3) &
                             &, cpl_SmpCN_R(i1,i2,i3), gam )
     Spm(i1)%gi2(i_count) = gam
     Spm(i1)%id2(i_count,1) = ChiPm(i2)%id
     Spm(i1)%id2(i_count,2) = Chi0(i3)%id 
     i_count = i_count + 1
    End Do
   End Do


   !-----------------------
   ! pseudoscalar + W
   !-----------------------
   Do i2 = 2,n_P0
    Do i3=1,n_W
     coupC = cpl_SmpP0W(i1, i2,i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in, mP0(i2), mW(i3), coupC, gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = P0(i2)%id
      Spm(i1)%id2(i_count,2) = id_W(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !-------------------------------
   ! charged scalar + pseudoscalar
   !-------------------------------
   Do i2 = 2,i1-1
    Do i3 = 2,n_P0
     coupC = cpl_SmpP03(i1, i2, i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSpm(i2), mP0(i3), coupC, gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = Spm(i2)%id
      Spm(i1)%id2(i_count,2) = P0(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do

   !--------------
   ! scalar + W
   !--------------
   Do i2 = 1,n_S0
    Do i3=1,n_W
     coupC = cpl_SmpS0W(i1, i2,i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in, mS0(i2), mW(i3), coupC, gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = S0(i2)%id
      Spm(i1)%id2(i_count,2) = id_W(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !-------------------------
   ! charged scalar + scalar
   !-------------------------
   Do i2 = 2,i1-1
    Do i3 = 1,n_S0
     coupC = cpl_SmpS03(i1, i2, i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSpm(i2), mS0(i3), coupC, gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = Spm(i2)%id
      Spm(i1)%id2(i_count,2) = S0(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do

   !----------------------
   ! charged scalar + Z
   !---------------------
   Do i2 = 2,i1-1
    Do i3=1,n_Z
    coupC = cpl_SmpZ(i1, i2, i3)
    If (Abs(CoupC).ne.0._dp) then
     Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mZ(i3), coupC, gam )
      Spm(i1)%gi2(i_count) = gam
      Spm(i1)%id2(i_count,1) = Spm(i2)%id
      Spm(i1)%id2(i_count,2) = id_Z(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   Spm(i1)%g = Sum(Spm(i1)%gi2)
   If (Spm(i1)%g.Ne.0._dp) Spm(i1)%bi2 = Spm(i1)%gi2 / Spm(i1)%g


  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine ChargedscalarTwoBodyDecays


 Subroutine CharginoTwoBodyDecays(i_in, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u &
  & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_s0     &
  & , n_p0, n_Spm, id_grav, ChiPm, Slept, c_CNuSl_L, c_CNuSl_R, Sneut          &
  & , c_CLSn_L, c_CLSn_R, mf_l, Sdown, c_CUSd_L, c_CUSd_R, mf_u, Sup, c_CDSu_L &
  & , c_CDSu_R, mf_d, Chi0, mW, c_CNW_L, c_CNW_R, Spm, c_SmpCN_L, c_SmpCN_R    &
  & , mZ, c_CCZ_L, c_CCZ_R, P0, c_CCP0_L, c_CCP0_R, S0, c_CCS0_L, c_CCS0_R     &
  & , m32, c_CGW_L, c_CGW_R, k_neut )
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of charginos:
 ! input:
 !  i_in ................. specifies the decaying chargino. The decays of
 !                all charginos are calculated if n_in < 0. In the case of
 !                generation diagonal models, the charginos are ordered
 !                according to their generation
 !  ChiPm(i) ............ chargino masses + decay information
 !  Slept(i) ........... slepton masses + decay information
 !  c_CNuSl_L(i,j,k) .. left chargino neutrino slepton coupling
 !  c_CNuSl_R(i,j,k) .. right chargino neutrino slepton coupling
 !  mSneutrino(i) ....... sneutrino masses + decay information
 !  c_CLSn_L(i,j,k) ... left chargino lepton sneutrino coupling
 !  c_CLSn_R(i,j,k) ... right chargino lepton sneutrino coupling
 !  mf_l(i) ............. lepton masses
 !  Sdown(i) ............ d-squark masses + decay information
 !  c_CUSd_L(i,j,k) ... left chargino u-quark d-squark coupling
 !  c_CUSd_R(i,j,k) ... right chargino u-quark d-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  Sup(i) .............. u-squark masses + decay information
 !  c_CDSu_L(i,j,k) ... left chargino d-quark u-squark coupling
 !  c_CDSu_R(i,j,k) ... right chargino d-quark u-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mN(i) ............... neutralino masses + decay information
 !  mW .................. mass of the W-boson
 !  c_CNW_L(i,j) ...... left chargino neutralino W coupling
 !  c_CNW_R(i,j) ...... right chargino neutralino W coupling
 !  Spm(i) .............. masses of charged scalars + decay information
 !  c_SmpCNW_L(i,j,k) . left charged scalar - chargino -neutralino coupling
 !  c_SmpCNW_R(i,j,k) . right charged scalar - chargino - neutralino coupling
 !  mZ .................. mass of the Z-boson
 !  c_CCZ_L(i,j) ...... left chargino Z coupling
 !  c_CCZ_R(i,j) ...... right chargino Z coupling
 !  P0(i) ............... pseudoscalar masses + decay information
 !  c_CCP0_L(i,j,k) ... left chargino pseudoscalar coupling
 !  c_CCP0_R(i,j,k) ... right chargino pseudoscalar coupling
 !  S0(i) ............... scalar masses + decay information
 !  c_CCS0_L(i,j,k) ... left chargino scalar coupling
 !  c_CCS0_R(i,j,k) ... right chargino scalar coupling
 !  k_neut ................ if =1 .... summing over neutrinos 
 !                          if =2 .... summing over all SM-fermions 
 ! written by Werner Porod, 26.04.2001
 ! 19.09.2010: adjusting to new variable types
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut, n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu  &
       & , n_sle, n_Sd, n_su, n_n, n_c, n_s0, n_p0, n_Spm, id_grav
  Integer, Intent(in) :: id_nu(:), id_l(:), id_d(:), id_u(:), id_Z(:), id_W(:)
  Real(dp), Intent(in) ::  mf_l(:), mf_d(:), mf_u(:), mW(:), mZ(:), m32
  Complex(dp), Intent(in) :: c_CNuSl_L(:,:,:), c_CNuSl_R(:,:,:)     &
         & , c_CLSn_L(:,:,:), c_CLSn_R(:,:,:), c_CUSd_L(:,:,:)      &
         & , c_CUSd_R(:,:,:), c_CDSu_L(:,:,:), c_CDSu_R(:,:,:)      &
         & , c_CNW_L(:,:,:), c_CNW_R(:,:,:), c_SmpCN_L(:,:,:)       &
         & , c_SmpCN_R(:,:,:), c_CCZ_L(:,:,:), c_CCZ_R(:,:,:)       &
         & , c_CCP0_L(:,:,:), c_CCP0_R(:,:,:), c_CCS0_L(:,:,:)      &
         & , c_CCS0_R(:,:,:), c_CGW_L(2), c_CGW_R(2)
  Type(particle2), Intent(in) :: Sdown(:), Spm(:), P0(:)
  Type(particle23), Intent(in) :: Sneut(:), Slept(:), Sup(:), Chi0(:), S0(:)
  Type(particle23), Intent(inout) :: ChiPm(:)
 
 
  Integer :: i1, i2, i_start, i_end, i_count, i3
  Real(dp) :: gam, m_in, mN(n_n), mSlepton(n_sle), mSneut(n_snu), mSdown(n_sd) &
         & , mSup(n_su), mC(n_c), mS0(n_S0), mSpm(n_Spm), mP0(n_P0), x1, x2, sq1
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CharginoTwoBodyDecays'


  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_c

   ChiPm%g = 0._dp
   Do i1=1,n_n
    ChiPm(i1)%gi2 = 0._dp
    ChiPm(i1)%bi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_c) ) Then 
   i_start = i_in 
   i_end = i_in
   ChiPm(i_in)%g = 0._dp
   ChiPm(i_in)%gi2 = 0._dp
   ChiPm(i_in)%bi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_c) = ',i_in,n_c
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mC = ChiPm%m
  mN = Chi0%m
  mSlepton = Slept%m
  mSneut = Sneut%m
  mSdown = Sdown%m
  mSup = Sup%m
  mP0 = P0%m
  mS0 = S0%m
  mSpm = Spm%m

  Do i1 = i_start, i_end
   m_in = mC(i1)
   If (Abs(m_in).Eq.0._dp) Cycle ! massless particle

   i_count = 1
   !--------------------------------------------------------
   ! Slepton neutrino, summing over neutrinos if k_neut=1
   !--------------------------------------------------------
   Do i2 = 1,n_sle
    Do i3 = 1, n_nu
     If ((Abs(c_CNuSl_L(i1,i3,i2))+Abs(c_CNuSl_R(i1,i3,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, 0._dp, mSlepton(i2) &
             & , c_CNuSl_L(i1,i3,i2), c_CNuSl_R(i1,i3,i2), gam)
      If (k_neut.Eq.1) Then
       ChiPm(i1)%gi2(i_count) = ChiPm(i1)%gi2(i_count) + gam
       ChiPm(i1)%id2(i_count,1) = Slept(i2)%id+1
       ChiPm(i1)%id2(i_count,2) = id_nu( (i2+1)/2 )
      Else
       ChiPm(i1)%gi2(i_count) = gam
       ChiPm(i1)%id2(i_count,1) = Slept(i2)%id+1
       ChiPm(i1)%id2(i_count,2) = id_nu(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
    If (k_neut.Eq.1)  i_count = i_count + 1
   End Do
   !-----------------------------------------------------------
   ! Sneutrino lepton
   !-----------------------------------------------------------
   Do i2 = 1,n_snu
    Do i3 = 1,n_l
     If ((Abs(c_CLSn_L(i1,i3,i2))+Abs(c_CLSn_R(i1,i3,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_l(i3), mSneut(i2) &
             & , c_CLSn_L(i1,i3,i2), c_CLSn_R(i1,i3,i2), gam)
      ChiPm(i1)%gi2(i_count) = gam
      ChiPm(i1)%id2(i_count,1) = Sneut(i2)%id
      ChiPm(i1)%id2(i_count,2) = id_l(i3)+1
      i_count = i_count + 1
     End If
    End Do
   End Do
   !----------------------------------------------------
   ! u-Squark d-quark
   !----------------------------------------------------
   Do i2 = 1,n_su
    Do i3 = 1,n_u
     If ((Abs(c_CDSu_L(i1,i3,i2))+Abs(c_CDSu_R(i1,i3,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_d(i3), mSup(i2) &
             & , c_CDSu_L(i1,i3,i2), c_CDSu_R(i1,i3,i2), gam)

      ChiPm(i1)%gi2(i_count) = 3._dp * gam ! colour 
      ChiPm(i1)%id2(i_count,1) = Sup(i2)%id
      ChiPm(i1)%id2(i_count,2) = id_d(i3)+1
      i_count = i_count + 1
     End If
    End Do
   End Do
   !----------------------------------------------------
   ! d-Squark u-quark
   !----------------------------------------------------
   Do i2 = 1,n_Sd
    Do i3 = 1,n_d
     If ((Abs(c_CUSd_L(i1,i3,i2))+Abs(c_CUSd_R(i1,i3,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_u(i3), mSdown(i2) &
             & , c_CUSd_L(i1,i3,i2), c_CUSd_R(i1,i3,i2), gam)

      ChiPm(i1)%gi2(i_count) = 3._dp * gam ! colour 
      ChiPm(i1)%id2(i_count,1) = Sdown(i2)%id+1
      ChiPm(i1)%id2(i_count,2) = id_u(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !------------------
   ! neutralino W
   !------------------
   Do i2 =1, n_n
    Do i3 =1, n_W
     Call FermionToFermionVectorBoson(m_in, mN(i2), mW(i3) &
           & , c_CNW_L(i1,i2,i3), c_CNW_R(i1,i2,i3), gam)

     ChiPm(i1)%gi2(i_count) = gam
     ChiPm(i1)%id2(i_count,1) = Chi0(i2)%id
     ChiPm(i1)%id2(i_count,2) = id_W(i3)
     i_count = i_count + 1
    End Do
   End Do

   !---------------------------
   ! charged scalar neutralino
   !--------------------------
   Do i2 = 2, n_Spm
    Do i3 = 1, n_n
     Call FermionToFermionScalar(m_in, mN(i3), mSpm(i2) &
            & ,c_SmpCN_L(i2,i1,i3), c_SmpCN_R(i2,i1,i3), gam)

     ChiPm(i1)%gi2(i_count) = gam
     ChiPm(i1)%id2(i_count,1) = SPm(i2)%id
     ChiPm(i1)%id2(i_count,2) = Chi0(i3)%id
     i_count = i_count + 1
    End Do
   End Do
   !------------------
   ! chargino Z
   !------------------
   Do i2 =1, i1-1
    Do i3 = 1, n_Z
     Call FermionToFermionVectorBoson(m_in, mC(i2), mZ(i3) &
           & , c_CCZ_L(i1,i2,i3), c_CCZ_R(i1,i2,i3), gam)
     ChiPm(i1)%gi2(i_count) = gam
     ChiPm(i1)%id2(i_count,1) = ChiPm(i2)%id
     ChiPm(i1)%id2(i_count,2) = id_Z(i3)
     i_count = i_count + 1
    End Do
   End Do

   !------------------------
   ! pseudoscalar neutralino
   !------------------------
   Do i2 = 2, n_P0
    Do i3 = 1, i1-1
     Call FermionToFermionScalar(m_in, mC(i3), mP0(i2) &
            & , c_CCP0_L(i1,i3,i2), c_CCP0_R(i1,i3,i2), gam)
     ChiPm(i1)%gi2(i_count) = gam
     ChiPm(i1)%id2(i_count,1) = ChiPm(i3)%id
     ChiPm(i1)%id2(i_count,2) = P0(i2)%id
     i_count = i_count + 1
    End Do
   End Do

   !-----------------
   ! scalar neutralino
   !-----------------
   Do i2 = 1, n_S0
    Do i3 = 1, i1-1
     Call FermionToFermionScalar(m_in, mC(i3), mS0(i2) &
            & , c_CCS0_L(i1,i3,i2), c_CCS0_R(i1,i3,i2), gam)
     ChiPm(i1)%gi2(i_count) = gam
     ChiPm(i1)%id2(i_count,1) = ChiPm(i3)%id
     ChiPm(i1)%id2(i_count,2) = S0(i2)%id
     i_count = i_count + 1
    End Do
   End Do
   !-----------------------------------------
   ! gravitino W
   !-----------------------------------------
   If (Abs(mC(i1)).Gt.(m32+mW(1))) Then ! to be changed
     x1 = m32/mC(i1)
     sq1 = Sqrt(x1)
     x2 = (mW(1)/mC(i1))**2
     ChiPm(i1)%gi2(i_count) = oo16pi * Abs(mC(i1))**5                         &
        &      * Sqrt(1._dp-2._dp*(x1+x2)+(x1-x2)**2)                         &
        &      * ( (Abs(c_CGW_L(1))**2 + Abs(c_CGW_R(1))**2)                  &
        &          * ( (1._dp-x1)**2 * (1._dp + 3._dp*x1)                     &
        &                              - x2 * (3._dp + x1**2                  &
        &                                  - 12._dp * x1 * sq1                &
        &                                  - x2 * (3._dp-x1-x2) ) )           &
        &         + (Abs(c_CGW_L(2))**2 + Abs(c_CGW_R(2))**2)                 &
        &           * ( (1._dp-x1)**2 * (1._dp + sq1)**2                      &
        &             - x2 * ( (1-sq1)**2 * (3._dp + 2._dp*sq1 - 9._dp*x1)    &
        &                             - x2 * (3._dp-2._dp*sq1-9._dp*x1-x2) )) )

     ChiPm(i1)%id2(i_count,1) = id_grav
     ChiPm(i1)%id2(i_count,2) = id_W(1)

     i_count = i_count + 1
   End If

   ChiPm(i1)%g = Sum(ChiPm(i1)%gi2)
   If (ChiPm(i1)%g.Gt.0._dp) ChiPm(i1)%bi2 = ChiPm(i1)%gi2 / ChiPm(i1)%g

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine CharginoTwoBodyDecays


 Subroutine GluinoTwoBodyDecays(n_d, id_d, n_Sd, n_u, id_u, n_Su, id_grav, id_gl &
       & , Glu, Sdown, cpl_DGSd_L, cpl_DGSd_R, mf_d, Sup, cpl_UGSu_L, cpl_UGSu_R &
       & , mf_u, m32, Fgmsb, k_neut)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of gluinos at tree-level:
 ! input:
 !  mGlu ................ gluino masses
 !  mSdown(i) ........... d-squark masses
 !  cpl_DGSd_L(i,j,k) ... left d-quark gluino d-squark coupling
 !  cpl_DGSd_R(i,j,k) ... right d-quark gluino d-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_UGSu_L(i,j,k) ... left u-quark gluino u-squark coupling
 !  cpl_UGSu_R(i,j,k) ... right u-quark gluino u-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  k_neut .............. if =1 .... summing over all quarks
 !  GenerationMixing .... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output:
 !  depends on the values of k_neut and GenerationMixing.
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 26.04.2001
 ! 24.08.03: up to now there has been a sum over charged conjuagted states
 !           due to the need for the Les Houches Interface this will be
 !           removed
 ! 05.03.2015: adding decay into gravitino and gluon
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: k_neut, n_d, id_d(:), n_Sd, n_u, id_u(:), n_Su &
      & , id_grav, id_gl
  Real(dp), Intent(in) ::  mf_d(:), mf_u(:), m32, Fgmsb
  Complex(dp), Intent(in) :: cpl_DGSd_L(:,:), cpl_DGSd_R(:,:)   &
                          &, cpl_UGSu_L(:,:), cpl_UGSu_R(:,:)

  Type(particle2), Intent(in) :: Sdown(:)
  Type(particle23), Intent(in) :: Sup(:)
  Type(particle23), Intent(inout) :: Glu

  Integer :: i1, i2, i_count 
  Real(dp) :: mSdown(n_Sd), mGlu, gam, mSup(n_su), x1
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluinoTwoBodyDecays'

  mGlu = Glu%m

  If (mglu.Eq.0._dp) Then ! massless particle
   Iname = Iname - 1
   Return
  End If
  
  mSdown = Sdown%m
  mSup = Sup%m
  Glu%gi2 = 0._dp

  i_count = 1

  !----------------------------------------------------
  ! u-Squark u-quark, summing over quarks if k_neut=1
  !----------------------------------------------------
  Do i1 = 1,n_su
   Do i2 = 1,n_u
    If ( (Abs(cpl_UGSu_L(i2,i1))+Abs(cpl_UGSu_R(i2,i1))).Ne.0._dp) Then
     Call  FermionToFermionScalar(mglu, mf_u(i2), mSup(i1) &
             & , cpl_UGSu_L(i2,i1), cpl_UGSu_R(i2,i1), gam)
     gam = 2._dp * gam ! colour
     If (k_neut.Eq.1) Then
      Glu%gi2(i_count) = Glu%gi2(i_count) + gam
      Glu%id2(i_count,1) = Sup(i1)%id
      Glu%id2(i_count,2) = id_u(1)+1
      Glu%gi2(i_count+1) = Glu%gi2(i_count+1) + gam
      Glu%id2(i_count+1,1) = Sup(i1)%id+1
      Glu%id2(i_count+1,2) = id_u(1)
     Else
      Glu%gi2(i_count) = gam
      Glu%id2(i_count,1) = Sup(i1)%id
      Glu%id2(i_count,2) = id_u(i2)+1
      Glu%gi2(i_count+1) = gam
      Glu%id2(i_count+1,1) = Sup(i1)%id+1
      Glu%id2(i_count+1,2) = id_u(i2)
      i_count = i_count + 2
     End If
    End If
   End Do
   If (k_neut.Eq.1)  i_count = i_count + 2
  End Do

  !----------------------------------------------------
  ! d-Squark d-quark, summing over quarks if k_neut=2
  !----------------------------------------------------
  Do i1 = 1,n_sd
   Do i2 = 1,n_d
    If ( (Abs(cpl_UGSu_L(i2,i1))+Abs(cpl_UGSu_R(i2,i1))).Ne.0._dp) Then
     Call FermionToFermionScalar(mglu, mf_d(i2), mSdown(i1) &
            & , cpl_DGSd_L(i2,i1), cpl_DGSd_R(i2,i1), gam)
     gam = 2._dp * gam ! colour
     If (k_neut.Eq.1) Then
      Glu%gi2(i_count) = Glu%gi2(i_count) + gam
      Glu%id2(i_count,1) = Sdown(i1)%id
      Glu%id2(i_count,2) = id_d(1)+1
      Glu%gi2(i_count+1) = Glu%gi2(i_count+1) + gam
      Glu%id2(i_count+1,1) = Sdown(i1)%id+1
      Glu%id2(i_count+1,2) = id_d(1)
     Else
      Glu%gi2(i_count) = gam
      Glu%id2(i_count,1) = Sdown(i1)%id
      Glu%id2(i_count,2) = id_d(i2)+1
      Glu%gi2(i_count+1) = gam
      Glu%id2(i_count+1,1) = Sdown(i1)%id+1
      Glu%id2(i_count+1,2) = id_d(i2)
      i_count = i_count + 2
     End If
    End If
   End Do
   If (k_neut.Eq.1) i_count = i_count + 2
  End Do
     
 
  !-----------------------------------------
  ! gravitino gluon
  !-----------------------------------------
  If (mglu.Gt.m32) Then
   x1 = m32 / mglu
   Glu%gi2(i_count) = oo16pi * Abs(mGlu)**5 / Fgmsb**2 &
                   &        * (1._dp-x1)**3 * (1._dp + 3._dp*x1)
   Glu%id2(i_count,1) = id_grav
   Glu%id2(i_count,2) = id_gl

   i_count = i_count + 1
  End If

  Glu%g = Sum(Glu%gi2)
  If (Glu%g.Ne.0._dp) Then
   Glu%bi2 = Glu%gi2 / Glu%g
  Else
   Glu%bi2 = 0._dp
  End If

  Iname = Iname - 1

 End Subroutine GluinoTwoBodyDecays


 Subroutine GluinoTwoBodyDecays_old(mGlu, mSdown, cpl_DGSd_L, cpl_DGSd_R, mf_d  &  
                              &, mSup, cpl_UGSu_L, cpl_UGSu_R, mf_u         &  
                              &, k_neut, GenerationMixing                   &
                              &, gP, gT, BR )
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of gluinos at tree-level:
 ! input:
 !  mGlu ................ gluino masses
 !  mSdown(i) ........... d-squark masses
 !  cpl_DGSd_L(i,j,k) ... left d-quark gluino d-squark coupling
 !  cpl_DGSd_R(i,j,k) ... right d-quark gluino d-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_UGSu_L(i,j,k) ... left u-quark gluino u-squark coupling
 !  cpl_UGSu_R(i,j,k) ... right u-quark gluino u-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  k_neut .............. if =1 .... summing over all quarks
 !  GenerationMixing .... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output:
 !  depends on the values of k_neut and GenerationMixing.
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 26.04.2001
 ! 24.08.03: up to now there has been a sum over charged conjuagted states
 !           due to the need for the Les Houches Interface this will be
 !           removed
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: k_neut
  Real(dp), Intent(in) ::  mGlu, mSdown(6), mf_d(3), mSup(6), mf_u(3)
  Complex(dp), Intent(in) :: cpl_DGSd_L(:,:), cpl_DGSd_R(:,:)   &
                          &, cpl_UGSu_L(:,:), cpl_UGSu_R(:,:)
  Real(dp), Intent(inout) :: gP(:), gT
  Real(dp), Optional, Intent(inout) :: BR(:)
  Logical, Intent(in) :: GenerationMixing
 
  Integer :: i1, i2, i_gen, i_count 
  Real(dp) :: gam
  Complex(dp) :: coupLC, coupRC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluinoTwoBodyDecays_old'

  gT = 0._dp
  gP = 0._dp

  If (mglu.Eq.0._dp) Then ! massless particle
   Iname = Iname - 1
   Return
  End If

  i_count = 1

  If (GenerationMixing) Then
   !----------------------------------------------------
   ! u-Squark u-quark, summing over quarks if k_neut=1
   !----------------------------------------------------
   Do i1 = 1,6
    Do i2 = 1,3
     coupLC = cpl_UGSu_L(i2,i1)
     coupRC = cpl_UGSu_R(i2,i1)
     Call  FermionToFermionScalar(mglu, mf_u(i2), mSup(i1), coupLC, coupRC &
                                &, gam)
     If (k_neut.Eq.1) Then
      gP(i_count) = gP(i_count) +  gam
      gP(i_count+1) = gP(i_count+1) +  gam
     Else
      gP(i_count) = gam 
      gP(i_count+1) = gam 
      i_count = i_count + 2
     End If
    End Do
    If (k_neut.Eq.1)  i_count = i_count + 2
   End Do
   !----------------------------------------------------
   ! d-Squark d-quark, summing over quarks if k_neut=2
   !----------------------------------------------------
   Do i1 = 1,6
    Do i2 = 1,3
     coupLC = cpl_DGSd_L(i2,i1)
     coupRC = cpl_DGSd_R(i2,i1)
     Call  FermionToFermionScalar(mglu, mf_d(i2), mSdown(i1), coupLC, coupRC &
                                &, gam)
     If (k_neut.Eq.1) Then
      gP(i_count) = gP(i_count) + gam
      gP(i_count+1) = gP(i_count+1) +  gam
     Else
      gP(i_count) = gam
      gP(i_count+1) = gam 
      i_count = i_count + 2
     End If
    End Do
    If (k_neut.Eq.1) i_count = i_count + 2
   End Do
     
  Else ! GenerationMixing = .FALSE.
   !------------------
   ! u-Squark u-quark
   !------------------
   Do i1 = 1,3
    Do i2 = 1,2
     i_gen = (i1-1)*2 + i2
     coupLC = cpl_UGSu_L(i1,i_gen)
     coupRC = cpl_UGSu_R(i1,i_gen)
     Call  FermionToFermionScalar(mglu, mf_u(i1), mSup(i_gen), coupLC &
                                &, coupRC, gam)
     gP(i_count) = gam
     gP(i_count+1) = gam
     i_count = i_count + 2
    End Do
   End Do
   !------------------
   ! d-Squark d-quark
   !------------------
   Do i1 = 1,3
    Do i2 = 1,2
     i_gen = (i1-1)*2 + i2
     coupLC = cpl_DGSd_L(i1,i_gen)
     coupRC = cpl_DGSd_R(i1,i_gen)
     Call  FermionToFermionScalar(mglu, mf_d(i1), mSdown(i_gen), coupLC &
                                &, coupRC, gam)
     gP(i_count) = gam
     gP(i_count+1) = gam
     i_count = i_count + 2
    End Do
   End Do
  End If ! GenerationMixing

  !--------------------------------------------------------------
  ! summation over colour gives a factor 2 
  !-------------------------------------------------------------
  gP = 2._dp * gP
  gT = Sum(gP)

  If ((Present(BR)).And.(gT.Eq.0)) Then
   BR = 0._dp
  Else If (Present(BR)) Then
   BR = gP / gT
  End If

  Iname = Iname - 1

 End Subroutine GluinoTwoBodyDecays_old

 Subroutine NeutralinoTwoBodyDecays(i_in, Chi0, n_nu, id_nu, n_l, id_l, n_d    &
   & , id_d, n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n    &
   & , n_c, n_s0, n_p0, n_Spm, id_ph, id_grav, Slept, c_LNSl_L, c_LNSl_R, mf_l &
   & , Sneut, c_NuNSn_L, c_NuNSn_R, Sdown, c_DNSd_L, c_DNSd_R, mf_d, Sup       &
   & , c_UNSu_L, c_UNSu_R, mf_u, ChiPm, mW, c_CNW_L, c_CNW_R, Spm, c_SmpCN_L   &
   & , c_SmpCN_R, mZ, c_NNZ_L, c_NNZ_R, P0, c_NNP0_L, c_NNP0_R, S0, c_NNS0_L   &
   & , c_NNS0_R, m32, c_NGP, c_NGZ, c_NGH, k_neut)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of neutralinos:
 ! input:
 !  i_in ................. specifies the decaying neutralino. The decays of
 !                all neutralinos are calculated if n_in < 0. In the case of
 !                generation diagonal models, the neutralinos are ordered
 !                according to their generation
 !  Chi0(i) ............. neutralino masses + decay information
 !  Slept(i) ............ slepton masses + decay information
 !  c_LNSl_L(i,j,k) ... left lepton neutralino slepton coupling
 !  c_LNSl_R(i,j,k) ... right lepton neutralino slepton coupling
 !  mf_l(i) ............. lepton masses
 !  Sneut(i) ............ sneutrino masses + decay information
 !  c_NuNSn_L(i,j,k) .. left neutrino neutralino sneutrino coupling
 !  c_NuNSn_R(i,j,k) .. right neutrino neutralino sneutrino coupling
 !  Sdown(i) ............ d-squark masses + decay information
 !  c_DNSd_L(i,j,k) ... left d-quark neutralino d-squark coupling
 !  c_DNSd_R(i,j,k) ... right d-quark neutralino d-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  Sup(i) .............. u-squark masses + decay information
 !  c_UNSu_L(i,j,k) ... left u-quark neutralino u-squark coupling
 !  c_UNSu_R(i,j,k) ... right u-quark neutralino u-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  ChiPm(i) ............ chargino masses + decay information
 !  mW .................. mass of the W-boson
 !  c_CNW_L(i,j) ...... left chargino neutralino W coupling
 !  c_CNW_R(i,j) ...... right chargino neutralino W coupling
 !  Spm(i) .............. masses of charged scalars + decay information
 !  c_SmpCNW_L(i,j,k) . left charged scalar - chargino -neutralino coupling
 !  c_SmpCNW_R(i,j,k) . right charged scalar - chargino - neutralino coupling
 !  mZ .................. mass of the Z-boson
 !  c_NNZ_L(i,j) ...... left neutralino Z coupling
 !  c_NNZ_R(i,j) ...... right neutralino Z coupling
 !  P0(i) ............... pseudoscalar masses + decay information
 !  c_NNP0_L(i,j,k) ... left neutralino pseudoscalar coupling
 !  c_NNP0_R(i,j,k) ... right neutralino pseudoscalar coupling
 !  S0(i) ............... scalar masses + decay information
 !  c_NNS0_L(i,j,k) ... left neutralino scalar coupling
 !  c_NNS0_R(i,j,k) ... right neutralino scalar coupling
 !  m32 ................ gravitino mass
 !  c_NGP, c_NGZ, c_NGH .. neutralino coupling to gravitino + either photon
 !                         or Z or h^0
 !  k_neut ................ if =1 .... summing over neutrinos, otherwise
 !                          no summation over neutrinos
 ! written by Werner Porod, 26.04.2001
 ! 10.09.03: give now explicitly branching ratios of charge conjugated states
 ! 17.09.10: using new type structures 
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut, n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu  &
       & , n_sle, n_Sd, n_su, n_n, n_c, n_s0, n_p0, n_Spm, id_ph, id_grav
  Integer, Intent(in) :: id_nu(:), id_l(:), id_d(:), id_u(:), id_Z(:), id_W(:)
  Real(dp), Intent(in) ::  mf_l(3), mf_d(3), mf_u(3), mW(:), mZ(:), m32
  Complex(dp), Intent(in) :: c_LNSl_L(:,:,:), c_LNSl_R(:,:,:)                 &
     & , c_NuNSn_L(:,:,:), c_NuNSn_R(:,:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:) &
     & , c_UNSu_L(:,:,:), c_UNSu_R(:,:,:), c_CNW_L(:,:,:), c_CNW_R(:,:,:)     &
     & , c_SmpCN_L(:,:,:), c_SmpCN_R(:,:,:), c_NNZ_L(:,:,:), c_NNZ_R(:,:,:)   &
     & , c_NNP0_L(:,:,:), c_NNP0_R(:,:,:), c_NNS0_L(:,:,:), c_NNS0_R(:,:,:)   &
     & , c_NGP, c_NGZ(2), c_NGH
  Type(particle2), Intent(in) :: Sdown(:), Spm(:), P0(:)
  Type(particle23), Intent(in) :: Sneut(:), Slept(:), Sup(:), ChiPm(:), S0(:)
  Type(particle23), Intent(inout) :: Chi0(:)

  Integer :: i1, i2, i3, i_start, i_end, i_count
  Real(dp) :: gam, m_in, mN(n_n), mSlepton(n_sle), mSneut(n_snu)     &
       & , mSdown(n_sd), mSup(n_su), mC(n_c), mS0(n_S0), mSpm(n_Spm) &
       & , mP0(n_P0), x1, x2, sq1
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoTwoBodyDecays'

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_n

   Chi0%g = 0._dp
   Do i1=1,n_n
    Chi0(i1)%gi2 = 0._dp
    Chi0(i1)%bi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_n) ) Then 
   i_start = i_in 
   i_end = i_in

   Chi0(i_in)%g = 0._dp
   Chi0(i_in)%gi2 = 0._dp
   Chi0(i_in)%bi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_n) = ',i_in,n_n
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mN = Chi0%m
  mC = ChiPm%m
  mSlepton = Slept%m
  mSneut = Sneut%m
  mSdown = Sdown%m
  mSup = Sup%m
  mP0 = P0%m
  mS0 = S0%m
  mSpm = Spm%m

  Do i1 = i_start, i_end
   m_in = mN(i1)
   If (Abs(m_in).Eq.0._dp) Cycle ! massless particle
   i_count = 1

   !----------------------------
   ! Slepton lepton
   !----------------------------
   Do i2 = 1,n_sle
    Do i3 = 1, n_l
     If ((Abs(c_LNSl_L(i3,i1,i2))+Abs(c_LNSl_R(i3,i1,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_l(i3), mSlepton(i2) &
             & , c_LNSl_L(i3,i1,i2), c_LNSl_R(i3,i1,i2), gam)
      Chi0(i1)%gi2(i_count) = gam
      Chi0(i1)%id2(i_count,1) = Slept(i2)%id
      Chi0(i1)%id2(i_count,2) = id_l(i3)+1
      Chi0(i1)%gi2(i_count+1) = gam
      Chi0(i1)%id2(i_count+1,1) = Slept(i2)%id+1
      Chi0(i1)%id2(i_count+1,2) = id_l(i3)
      i_count = i_count + 2
     End If
    End Do
   End Do
   !-----------------------------------------------------------
   ! Sneutrino neutrino, summing over neutrinos if k_neut=1
   !-----------------------------------------------------------
   Do i2 = 1,n_snu
    Do i3 = 1,n_nu
     If ((Abs(c_NuNSn_L(i3,i1,i2))+Abs(c_NuNSn_R(i3,i1,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, 0._dp, mSneut(i2) &
             & , c_NuNSn_L(i3,i1,i2), c_NuNSn_R(i3,i1,i2), gam)
      If (k_neut.Eq.1) Then
       Chi0(i1)%gi2(i_count) = Chi0(i1)%gi2(i_count) + gam
       Chi0(i1)%id2(i_count,1) = Sneut(i2)%id
       Chi0(i1)%id2(i_count,2) = id_nu(i2)+1
       Chi0(i1)%gi2(i_count+1) = Chi0(i1)%gi2(i_count+1) + gam
       Chi0(i1)%id2(i_count+1,1) = Sneut(i2)%id+1
       Chi0(i1)%id2(i_count+1,2) = id_nu(i2)
      Else
       Chi0(i1)%gi2(i_count) = gam
       Chi0(i1)%id2(i_count,1) = Sneut(i2)%id
       Chi0(i1)%id2(i_count,2) = id_nu(i3)+1
       Chi0(i1)%gi2(i_count+1) = gam
       Chi0(i1)%id2(i_count+1,1) = Sneut(i2)%id+1
       Chi0(i1)%id2(i_count+1,2) = id_nu(i3)
       i_count = i_count + 2
      End If
     End If
    End Do
    If (k_neut.Eq.1)  i_count = i_count + 2
   End Do
   !----------------------------------------------------
   ! u-Squark u-quark
   !----------------------------------------------------
   Do i2 = 1,n_su
    Do i3 = 1,n_u
     If ((Abs(c_UNSu_L(i3,i1,i2))+Abs(c_UNSu_R(i3,i1,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_u(i3), mSup(i2) &
             & , c_UNSu_L(i3,i1,i2), c_UNSu_R(i3,i1,i2), gam)

      Chi0(i1)%gi2(i_count) = 3._dp * gam ! colour 
      Chi0(i1)%id2(i_count,1) = Sup(i2)%id
      Chi0(i1)%id2(i_count,2) = id_u(i3)+1
      Chi0(i1)%gi2(i_count+1) = 3._dp * gam ! colour 
      Chi0(i1)%id2(i_count+1,1) = Sup(i2)%id+1
      Chi0(i1)%id2(i_count+1,2) = id_u(i3)

      i_count = i_count + 2
     End If
    End Do
   End Do
   !----------------------------------------------------
   ! d-Squark d-quark
   !----------------------------------------------------
   Do i2 = 1,n_Sd
    Do i3 = 1,n_d
     If ((Abs(c_DNSd_L(i3,i1,i2))+Abs(c_DNSd_R(i3,i1,i2))).Gt.0._dp) Then
      Call FermionToFermionScalar(m_in, mf_d(i3), mSdown(i2) &
             & , c_DNSd_L(i3,i1,i2), c_DNSd_R(i3,i1,i2), gam)

      Chi0(i1)%gi2(i_count) = 3._dp * gam ! colour 
      Chi0(i1)%id2(i_count,1) = Sdown(i2)%id
      Chi0(i1)%id2(i_count,2) = id_d(i3)+1
      Chi0(i1)%gi2(i_count+1) = 3._dp * gam ! colour 
      Chi0(i1)%id2(i_count+1,1) = Sdown(i2)%id+1
      Chi0(i1)%id2(i_count+1,2) = id_d(i3)

      i_count = i_count + 2
     End If
    End Do
   End Do

   !------------------
   ! chargino W
   !------------------
   Do i2 =1, n_c
    Do i3 =1, n_W
     Call FermionToFermionVectorBoson(m_in, mC(i2), mW(i3) &
           & , c_CNW_L(i2,i1,i3), c_CNW_R(i2,i1,i3), gam)

     Chi0(i1)%gi2(i_count) = gam
     Chi0(i1)%id2(i_count,1) = ChiPm(i2)%id
     Chi0(i1)%id2(i_count,2) = id_W(i3)+1
     Chi0(i1)%gi2(i_count+1) = gam
     Chi0(i1)%id2(i_count+1,1) = ChiPm(i2)%id+1
     Chi0(i1)%id2(i_count+1,2) = id_W(i3)

     i_count = i_count + 2
    End Do
   End Do

   !---------------------------
   ! charged scalar chargino
   !--------------------------
   Do i2 = 2, n_Spm
    Do i3 = 1, n_c
     Call FermionToFermionScalar(m_in, mC(i3), mSpm(i2) &
            & ,c_SmpCN_L(i2,i3,i1), c_SmpCN_R(i2,i3,i1), gam)

     Chi0(i1)%gi2(i_count) = gam
     Chi0(i1)%id2(i_count,1) = SPm(i2)%id
     Chi0(i1)%id2(i_count,2) = ChiPm(i3)%id+1
     Chi0(i1)%gi2(i_count+1) = gam
     Chi0(i1)%id2(i_count+1,1) = SPm(i2)%id+1
     Chi0(i1)%id2(i_count+1,2) = ChiPm(i3)%id

     i_count = i_count + 2
    End Do
   End Do
   !------------------
   ! neutralino Z
   !------------------
   Do i2 =1, i1-1
    Do i3 = 1, n_Z
     Call FermionToFermionVectorBoson(m_in, mN(i2), mZ(i3) &
           & , c_NNZ_L(i1,i2,i3), c_NNZ_R(i1,i2,i3), gam)

     Chi0(i1)%gi2(i_count) = gam
     Chi0(i1)%id2(i_count,1) = Chi0(i2)%id
     Chi0(i1)%id2(i_count,2) = id_Z(i3)

     i_count = i_count + 1
    End Do
   End Do

   !------------------------
   ! pseudoscalar neutralino
   !------------------------
   Do i2 = 2, n_P0
    Do i3 = 1, i1-1
     Call FermionToFermionScalar(m_in, mN(i3), mP0(i2) &
            & , c_NNP0_L(i1,i3,i2), c_NNP0_R(i1,i3,i2), gam)
     Chi0(i1)%gi2(i_count) = gam
     Chi0(i1)%id2(i_count,1) = Chi0(i3)%id
     Chi0(i1)%id2(i_count,2) = P0(i2)%id

     i_count = i_count + 1
    End Do
   End Do

   !-----------------
   ! scalar neutralino
   !-----------------
   Do i2 = 1, n_S0
    Do i3 = 1, i1-1
     Call FermionToFermionScalar(m_in, mN(i3), mS0(i2) &
            & , c_NNS0_L(i1,i3,i2), c_NNS0_R(i1,i3,i2), gam)
     Chi0(i1)%gi2(i_count) = gam
     Chi0(i1)%id2(i_count,1) = Chi0(i3)%id
     Chi0(i1)%id2(i_count,2) = S0(i2)%id

     i_count = i_count + 1
    End Do
   End Do
   If (Abs(mN(i1)).Gt.m32) Then
    !-----------------------------------------
    ! gravitino photon
    !-----------------------------------------
    sq1 = Abs(m32/mN(i1))
    x1 = sq1**2
    Chi0(i1)%gi2(i_count) = oo16pi * Abs(c_NGP)**2 * Abs(mN(i1))**5 &
                          &        * (1._dp-x1)**3 * (1._dp + 3._dp*x1)
    Chi0(i1)%id2(i_count,1) = id_grav
    Chi0(i1)%id2(i_count,2) = id_ph

    i_count = i_count + 1
    !-----------------------------------------
    ! gravitino Z
    !-----------------------------------------
    If (Abs(mN(i1)).Gt.(m32+mZ(1))) Then ! to be changed
     x2 = (mZ(1)/mN(i1))**2
     Chi0(i1)%gi2(i_count) = oo16pi * Abs(mN(i1))**5                          &
        &      * Sqrt(1._dp-2._dp*(x1+x2)+(x1-x2)**2)                         &
        &      * ( Abs(c_NGZ(1))**2  * ( (1._dp-x1)**2 * (1._dp + 3._dp*x1)   &
        &                           - x2 * (3._dp + x1**2                     &
        &                                  - 12._dp * x1 * sq1                &
        &                                  - x2 * (3._dp-x1-x2) ) )           &
        &         + Abs(c_NGZ(2))**2 * ( (1._dp-x1)**2 * (1._dp + sq1)**2     &
        &                      - x2 * ( (1-sq1)**2                            &
        &                               * (3._dp + 2._dp*sq1 - 9._dp*x1)      &
        &                             - x2 * (3._dp-2._dp*sq1-9._dp*x1-x2) )) )
     Chi0(i1)%id2(i_count,1) = id_grav
     Chi0(i1)%id2(i_count,2) = id_Z(1)

     i_count = i_count + 1
    End If
   !-----------------------------------------
   ! gravitino h0
   !-----------------------------------------
    If (Abs(mN(i1)).Gt.(m32+mS0(1))) Then ! h0
     x2 = (mS0(1)/mN(i1))**2
     Chi0(i1)%gi2(i_count) = oo16pi * Abs(c_NGH)**2 * Abs(mN(i1))**5        &
                     &       * Sqrt(1._dp-2._dp*(x1+x2)+(x1-x2)**2)         &
                     &       * ( (1._dp-x1)**2 * (1._dp + sq1)**2           &
                     &         - x2 * ( (1+sq1)**2                          &
                     &                  * (3._dp - 2._dp*sq1 + 3._dp*x1)    & 
                     &                - x2 * (3._dp+2._dp*sq1+3._dp*x1-x2) ))
     Chi0(i1)%id2(i_count,1) = id_grav
     Chi0(i1)%id2(i_count,2) = S0(1)%id

    End If
   End If

   Chi0(i1)%g = Sum(Chi0(i1)%gi2)
   If (Chi0(i1)%g.Gt.0._dp) Chi0(i1)%bi2 = Chi0(i1)%gi2 / Chi0(i1)%g

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine NeutralinoTwoBodyDecays


 Subroutine PseudoscalarTwoBodyDecays(i_in, n_s0, n_nu, id_nu, n_l, id_l, n_d  &
       & , id_d, n_u, id_u, n_Z, id_Z, n_W, id_W, n_sle, n_Sd, n_su, n_n, n_c  &
       & , n_p0, n_Spm, id_ph, id_gl, P0, mf_l, cpl_LLP0_L, cpl_LLP0_R         &
       & , mf_d, cpl_DDP0_L, cpl_DDP0_R, mf_u, cpl_UUP0_L, cpl_UUP0_R          &
       & , Slept, cpl_P0SlSl, Sdown, cpl_P0SdSd, Sup, cpl_P0SuSu, Chi0         &
       & , cpl_NNP0_L, cpl_NNP0_R, ChiPm, cpl_CCP0_L, cpl_CCP0_R, Spm          &
       & , cpl_SmpP03, S0, cpl_P0S03, mZ, cpl_P0S0Z, mW, cpl_SmpP0W            &
       & , cpl_GlGlP0, cpl_GGP0, mglu)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of pseudoscalars
 ! input:
 !  i_in ................. specifies the decaying pseudoscalar. The decays of
 !                         all pseudoscalars are calculated if i_in < 0.
 !  mP0(i) .............. scalar masses
 !  mf_l(i) ............. lepton masses
 !  cpl_LLP0_L(i,j,k) ... left coupling lepton-lepton-pseudoscalar
 !  cpl_LLP0_R(i,j,k) ... right coupling lepton-lepton-pseudoscalar
 !  mf_d(i) ............. d-quark masses
 !  cpl_DDP0_L(i,j,k) ... left coupling d-quark d-quark pseudoscalar
 !  cpl_DDP0_R(i,j,k) ... right coupling d-quark d-quark pseudoscalar
 !  mf_u(i) ............. u-quark masses
 !  cpl_UUP0_L(i,j,k) ... left coupling u-quark u-quark pseudoscalar
 !  cpl_UUP0_R(i,j,k) ... right coupling u-quark u-quark pseudoscalar
 !  mSlepton(i) ......... slepton masses
 !  cpl_P0SlSl(i,j,k) ... coupling pseudoscalar slepton slepton
 !  mSdown(i) ........... d-squark masses
 !  cpl_P0SdSd(i,j,k) ... coupling pseudoscalar d-squark d-squark
 !  mSup(i) ............. u-squark masses
 !  cpl_P0SuSu(i,j,k) ... coupling pseudoscalar u-squark u-squark
 !  mN(i) ............... neutralino masses
 !  cpl_NNP0_L(i,j,k) ... left neutralino-neutralino-pseudoscalar coupling
 !  cpl_NNP0_R(i,j,k) ... right neutralino-neutralino-pseudoscalar coupling
 !  mC(i) ............... chargino masses
 !  cpl_CCP0_L(i,j,k) ... left chargino-chargino-pseudoscalar coupling
 !  cpl_CCP0_R(i,j,k) ... right chargino-chargino-pseudoscalar coupling
 !  mW .................. mass of the W-boson
 !  cpl_S0WW(i) ......... scalar-W-W coupling
 !  mZ ................... mass of the Z-boson
 !  cpl_S0ZZ(i) ......... scalar-Z-Z coupling
 !  mSpm(i) .............. masses of charged scalars
 !  cpl_SmpS03(i,j,k) .... charged scalar - charged scalar - scalar coupling
 !  mP0(i) ............... pseudoscalar masses
 !  cpl_P0S03(i,j,k) ..... pseudoscalar - pseudoscalar - scalar coupling
 !  cpl_P0S0Z(i,j) ....... pseudoscalar-scalar-Z coupling
 !  cpl_SmpS0W(i,j) ...... charged scalar - scalar - W coupling
 ! output: 
 ! written by Werner Porod, 30.04.2001
 ! 15.11.02: adding QCD corrections for decays into fermions
 ! 14.09.03: adding charge conjugated states to output
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, n_s0, n_nu, id_nu(:), n_l, id_l(:), n_d, id_d(:)  &
       & , n_u, id_u(:), n_Z, id_Z(:), n_W, id_W(:), n_sle, n_Sd, n_su, n_n, n_c &
       & , n_p0, n_Spm, id_ph, id_gl
  Real(dp), Intent(in) :: mf_l(3), mf_d(3), mf_u(3), mW(:), mZ(:), mglu
  Real(dp), Intent(in) ::  cpl_P0S03(:,:,:)
  Complex(dp), Intent(in) ::  cpl_LLP0_L(:,:,:), cpl_LLP0_R(:,:,:)       &
          & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:), cpl_UUP0_L(:,:,:)    &
          & , cpl_UUP0_R(:,:,:), cpl_P0SlSl(:,:,:), cpl_P0SdSd(:,:,:)    &
          & , cpl_P0SuSu(:,:,:), cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:)    &
          & , cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:), cpl_SmpP03(:,:,:)    &
          & , cpl_SmpP0W(:,:,:), cpl_P0S0Z(:,:,:), cpl_GlGlP0(:), cpl_GGP0(:)
  Type(particle2), Intent(in) :: Sdown(:), Spm(:)
  Type(particle23), Intent(in) :: Slept(:), Sup(:), Chi0(:), ChiPm(:), S0(:)
  Type(particle2), intent(inout) :: P0(:)

  Integer :: i1, i2, i_start, i_end, i_count, i3
  Real(dp) :: m_in, alpha_3, gam, mSlepton(n_sle), mSdown(n_sd)  &
    & , mSup(n_su), mN(n_n), mC(n_c), mP0(n_P0), mSpm(n_Spm), mS0(n_S0)
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoscalarTwoBodyDecays'

  If (i_in.Lt.0) Then
   i_start = 2
   i_end = n_P0
   Do i1=1,n_S0
    P0(i1)%gi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_P0) ) Then 
   i_start = i_in 
   i_end = i_in
   P0(i_in)%gi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_P0) = ',i_in,n_P0
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mP0 = P0%m
  mSlepton = Slept%m
  mSdown = Sdown%m
  mSup = Sup%m
  mN = Chi0%m
  mC = ChiPm%m
  mSPm = Spm%m
  mS0 = S0%m

  Do i1 = i_start, i_end
   m_in = mP0(i1)
   i_count = 1
   If (m_in.Eq.0._dp) Cycle ! massless particle

   alpha_3 = AlphaSDR(m_in, mGlu, mSup, mSdown) ! needed for QCD corrections
   !------------------
   ! into leptons
   !------------------
   Do i2 = 1,n_l
    Do i3 = i2,n_l
     If ((Abs(cpl_LLP0_L(i2,i3,i1))+Abs(cpl_LLP0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i3), cpl_LLP0_L(i2,i3,i1) &
                             &, cpl_LLP0_R(i2,i3,i1), gam )
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = id_l(i2)
      P0(i1)%id2(i_count,2) = id_l(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = gam
       P0(i1)%id2(i_count,1) = id_l(i2) + 1
       P0(i1)%id2(i_count,2) = id_l(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do

   !------------------
   ! into d-quarks
   !------------------
   Do i2 = 1, n_d
    Do i3 = i2,n_d
     If ((Abs(cpl_DDP0_L(i2,i3,i1))+Abs(cpl_DDP0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i3), cpl_DDP0_L(i2,i3,i1) &
                             &, cpl_DDP0_R(i2,i3,i1), gam )
      gam = 3._dp * gam * FFqcd(mf_d(i2), m_in, alpha_3)
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = id_d(i2)
      P0(i1)%id2(i_count,2) = id_d(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = gam
       P0(i1)%id2(i_count,1) = id_d(i2) + 1
       P0(i1)%id2(i_count,2) = id_d(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into u-quarks
   !------------------
   Do i2 = 1,n_u
    Do i3 = i2,n_u
     If ((Abs(cpl_UUP0_L(i2,i3,i1))+Abs(cpl_UUP0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i3), cpl_UUP0_L(i2,i3,i1) &
                             &, cpl_UUP0_R(i2,i3,i1), gam )
      gam = 3._dp * gam * FFqcd(mf_u(i2), m_in, alpha_3)
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = id_u(i2)
      P0(i1)%id2(i_count,2) = id_u(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = gam
       P0(i1)%id2(i_count,1) = id_u(i2) + 1
       P0(i1)%id2(i_count,2) = id_u(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into sleptons
   !------------------
   Do i2 = 1,n_sle
    Do i3 = i2,n_sle
     If (Abs(cpl_P0SlSl(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSlepton(i3)   &
                             &, cpl_P0SlSl(i1,i2,i3), gam )
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = Slept(i2)%id
      P0(i1)%id2(i_count,2) = Slept(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = gam
       P0(i1)%id2(i_count,1) = Slept(i2)%id + 1
       P0(i1)%id2(i_count,2) = Slept(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into d-squarks
   !------------------
   Do i2 = 1, n_sd
    Do i3 = i2, n_sd
     If (Abs(cpl_P0SdSd(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSdown(i3)   &
                             &, cpl_P0SdSd(i1,i2,i3), gam )
      P0(i1)%gi2(i_count) = 3._dp * gam
      P0(i1)%id2(i_count,1) = Sdown(i2)%id
      P0(i1)%id2(i_count,2) = Sdown(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = 3._dp * gam
       P0(i1)%id2(i_count,1) = Sdown(i2)%id + 1
       P0(i1)%id2(i_count,2) = Sdown(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into u-squarks
   !------------------
   Do i2 = 1, n_su
    Do i3 = i2, n_su
     If (Abs(cpl_P0SuSu(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSup(i2), mSup(i3)   &
                             &, cpl_P0SuSu(i1,i2,i3), gam )
      P0(i1)%gi2(i_count) = 3._dp * gam
      P0(i1)%id2(i_count,1) = Sup(i2)%id
      P0(i1)%id2(i_count,2) = Sup(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = 3._dp * gam
       P0(i1)%id2(i_count,1) = Sup(i2)%id + 1
       P0(i1)%id2(i_count,2) = Sup(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !-------------
   ! Neutralinos
   !-------------
   Do i2 = 1, n_n
    Do i3 = i2, n_n
     Call ScalarToTwoFermions(m_in, mN(i2), mN(i3), cpl_NNP0_L(i2,i3, i1) &
                             &, cpl_NNP0_R(i2,i3, i1), gam )
     If (i2.Eq.i3) gam = 0.5_dp * gam ! Majorana
     P0(i1)%gi2(i_count) = gam
     P0(i1)%id2(i_count,1) = Chi0(i2)%id
     P0(i1)%id2(i_count,2) = Chi0(i3)%id 
     i_count = i_count + 1
    End Do
   End Do
   !-------------
   ! Charginos
   !-------------
   Do i2 = 1, n_c
    Do i3 = i2, n_c
     Call ScalarToTwoFermions(m_in, mC(i2), mC(i3), cpl_CCP0_L(i2,i3, i1) &
                             &, cpl_CCP0_R(i2,i3, i1), gam )
     P0(i1)%gi2(i_count) = gam
     P0(i1)%id2(i_count,1) = ChiPm(i2)%id
     P0(i1)%id2(i_count,2) = ChiPm(i3)%id + 1
     i_count = i_count + 1
     If (i2.Ne.i3) Then
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = ChiPm(i2)%id + 1
      P0(i1)%id2(i_count,2) = ChiPm(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !--------------------
   ! charged scalar + W
   !--------------------
   Do i2 = 2,n_Spm
    Do i3=1,n_W
     coupC = cpl_SmpP0W(i2, i1, i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mW(i3), coupC, gam )
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = SPm(i2)%id
      P0(i1)%id2(i_count,2) = id_W(i3) + 1
      i_count = i_count + 1
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = SPm(i2)%id + 1
      P0(i1)%id2(i_count,2) = id_W(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !-------------------
   ! 2 charged scalars
   !-------------------
   Do i2 = 2,n_Spm
    Do i3 = i2,n_Spm
     coupC = cpl_SmpP03(i2, i3, i1)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSpm(i2), mSpm(i3), coupC, gam)
      P0(i1)%gi2(i_count) = gam
      P0(i1)%id2(i_count,1) = SPm(i2)%id
      P0(i1)%id2(i_count,2) = SPm(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       P0(i1)%gi2(i_count) = gam
       P0(i1)%id2(i_count,1) = SPm(i2)%id + 1
       P0(i1)%id2(i_count,2) = SPm(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !-------------------
   ! scalars + Z
   !-------------------
   Do i2 = 1,n_S0
    Do i3=1,n_Z
     coupC = cpl_P0S0Z(i1, i2, i3)
     Call ScalarToScalarVectorBoson(m_in, mS0(i2), mZ(i3), coupC, gam )
     P0(i1)%gi2(i_count) = gam
     P0(i1)%id2(i_count,1) = S0(i2)%id
     P0(i1)%id2(i_count,2) = id_Z(i3)
     i_count = i_count + 1
    End Do
   End Do
   !------------------------
   ! pseudoscalar + scalar
   !------------------------
   Do i2 = 2,i1-1
    Do i3 = 1,n_S0
     coupC = cpl_P0S03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mP0(i2), mS0(i3), coupC, gam )
     P0(i1)%gi2(i_count) = gam
     P0(i1)%id2(i_count,1) = P0(i2)%id
     P0(i1)%id2(i_count,2) = S0(i3)%id
     i_count = i_count + 1
    End Do
   End Do

   !--------------
   ! two gluons
   !--------------
   P0(i1)%gi2(i_count) = G_F * m_in**3 * oosqrt2 * oo36pi3 *Abs(cpl_GlGlP0(i1))**2
   P0(i1)%id2(i_count,1) = id_gl
   P0(i1)%id2(i_count,2) = id_gl
   i_count = i_count + 1
   !--------------
   ! two photons
   !--------------
   P0(i1)%gi2(i_count) = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cpl_GGP0(i1))**2
   P0(i1)%id2(i_count,1) = id_ph
   P0(i1)%id2(i_count,2) = id_ph
   i_count = i_count + 1

   P0(i1)%g = Sum(P0(i1)%gi2) 
   If (P0(i1)%g.Ne.0._dp) P0(i1)%bi2 = P0(i1)%gi2 / P0(i1)%g

  End Do ! i1
 
  Iname = Iname - 1
 contains

  Real(dp) Function FFqcd(mf, mA, alpha_s)
  implicit none
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   if (ratio.ge.0.5_dp) return ! decay is kinematically forbitten

   if (ratio.ge.0.495_dp) return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (19._dp + 2._dp * beta2 + 3._dp * beta**4) / (16._dp * beta)      &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta
    fac =  fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  end  Function FFqcd

 End Subroutine PseudoscalarTwoBodyDecays


 Subroutine ScalarTwoBodyDecays(i_in, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d  &
       & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c &
       & , n_p0, n_Spm, id_ph, id_gl, S0, cpl_S03, cpl_GlGlS0, cpl_GGS0, mf_l  &
       & , cpl_LLS0_L, cpl_LLS0_R, mf_d, cpl_DDS0_L, cpl_DDS0_R, mf_u          & 
       & , cpl_UUS0_L, cpl_UUS0_R, Slept, cpl_S0SlSl, Sneut, cpl_S0SnSn        &
       & , Sdown, cpl_S0SdSd, Sup, cpl_S0SuSu, Chi0, cpl_NNS0_L, cpl_NNS0_R    & 
       & , ChiPm, cpl_CCS0_L, cpl_CCS0_R, mW, cpl_S0WW, cpl_S0WWvirt, mZ       &
       & , cpl_S0ZZ, cpl_S0ZZvirt, Spm, cpl_SmpS03, P0, cpl_P0S03, cpl_P0S0Z   &
       & , cpl_SmpS0W, mglu)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of scalars
 ! input:
 !  i_in ................. specifies the decaying scalar. The decays of
 !                         all scalarss are calculated if i_in < 0.
 !  mS0(i) .............. scalar masses
 !  cpl_S03(i,j,k) ...... trilinear self interaction of scalars
 !  cpl_GlGlS0(i1) ...... gluon-gluon-scalar coupling
 !  cpl_GGS0(i1) ........ photon-photon-scalar coupling
 !  mf_l(i) ............. lepton masses
 !  cpl_LLS0_L(i,j,k) ... left coupling lepton-lepton-scalar
 !  cpl_LLS0_R(i,j,k) ... right coupling lepton-lepton-scalar
 !  mf_d(i) ............. d-quark masses
 !  cpl_DDS0_L(i,j,k) ... left coupling d-quark d-quark scalar
 !  cpl_DDS0_R(i,j,k) ... right coupling d-quark d-quark scalar
 !  mf_u(i) ............. u-quark masses
 !  cpl_UUS0_L(i,j,k) ... left coupling u-quark u-quark scalar
 !  cpl_UUS0_R(i,j,k) ... right coupling u-quark u-quark scalar
 !  mSlepton(i) ......... slepton masses
 !  cpl_S0SlSl(i,j,k) ... coupling scalar slepton slepton
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_S0SnSn(i,j,k) ... coupling scalar sneutrino sneutrino
 !  mSdown(i) ........... d-squark masses
 !  cpl_S0SdSd(i,j,k) ... coupling scalar d-squark d-squark
 !  mSup(i) ............. u-squark masses
 !  cpl_S0SuSu(i,j,k) ... coupling scalar u-squark u-squark
 !  mN(i) ............... neutralino masses
 !  cpl_NNS0_L(i,j,k) ... left neutralino-neutralino-scalar coupling
 !  cpl_NNS0_R(i,j,k) ... right neutralino-neutralino-scalar coupling
 !  mC(i) ............... chargino masses
 !  cpl_CCS0_L(i,j,k) ... left chargino-chargino-scalar coupling
 !  cpl_CCS0_R(i,j,k) ... right chargino-chargino-scalar coupling
 !  mW .................. mass of the W-boson
 !  cpl_S0WW(i) ......... scalar-W-W coupling
 !  mZ ................... mass of the Z-boson
 !  cpl_S0ZZ(i) ......... scalar-Z-Z coupling
 !  mSpm(i) .............. masses of charged scalars
 !  cpl_SmpS03(i,j,k) .... charged scalar - charged scalar - scalar coupling
 !  mP0(i) ............... pseudoscalar masses
 !  cpl_P0S03(i,j,k) ..... pseudoscalar - pseudoscalar - scalar coupling
 !  cpl_P0S0Z(i,j) ....... pseudoscalar-scalar-Z coupling
 !  cpl_SmpS0W(i,j) ...... charged scalar - scalar - W coupling
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output: 
 !  depends on the value of  GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 30.04.2001
 ! 15.11.02: adding QCD corrections for decays into fermions
 ! 14.09.03: adding charge conjugated states to output
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, n_s0, n_nu, id_nu(:), n_l, id_l(:), n_d, id_d(:) &
       & , n_u, id_u(:), n_Z, id_Z(:), n_W, id_W(:), n_snu, n_sle, n_Sd, n_su   &
       & , n_n, n_c, n_p0, n_Spm, id_ph, id_gl
  Real(dp), Intent(in) :: mf_l(3), mf_d(3), mf_u(3), mW(:), mZ(:), mglu
  Real(dp), Intent(in) :: cpl_S03(:,:,:), cpl_S0WW(:,:), cpl_S0ZZ(:,:)     &
          & , cpl_P0S03(:,:,:), cpl_S0WWvirt(:,:), cpl_S0ZZvirt(:,:)
  Complex(dp), Intent(in) :: cpl_GlGlS0(:), cpl_GGS0(:), cpl_LLS0_L(:,:,:) &
          & , cpl_LLS0_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)      &
          & , cpl_UUS0_L(:,:,:), cpl_UUS0_R(:,:,:), cpl_S0SlSl(:,:,:)      &
          & , cpl_S0SnSn(:,:,:), cpl_S0SdSd(:,:,:), cpl_S0SuSu(:,:,:)      &
          & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_CCS0_L(:,:,:)      &
          & , cpl_CCS0_R(:,:,:), cpl_SmpS03(:,:,:), cpl_SmpS0W(:,:,:)      &
          & , cpl_P0S0Z(:,:,:)
  Type(particle2), Intent(in) :: Sdown(:), P0(:), Spm(:)
  Type(particle23), Intent(in) :: Sneut(:), Slept(:), Sup(:), Chi0(:), ChiPm(:)
  Type(particle23), intent(inout) :: S0(:)

  Integer :: i1, i2, i3, i_start, i_end, i_count
  Real(dp) :: mS0(n_S0), m_in, alpha_3, gam, mSlepton(n_sle), mSpm(n_Spm)  &
    & , mSneutrino(n_snu), mSdown(n_sd), mSup(n_su), mN(n_n), mC(n_c)      &
    & , mP0(n_P0)
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarTwoBodyDecays'

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_S0
   Do i1=1,n_S0
    S0(i1)%gi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_S0) ) Then 
   i_start = i_in 
   i_end = i_in
   S0(i_in)%gi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_S0) = ',i_in,n_S0
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mS0 = S0%m
  mSlepton = Slept%m
  mSneutrino = Sneut%m
  mSdown = Sdown%m
  mSup = Sup%m
  mN = Chi0%m
  mC = ChiPm%m
  mSPm = Spm%m
  mP0 = P0%m

  Do i1 = i_start, i_end

   S0(i1)%gi2 = 0._dp
   m_in = mS0(i1)
   i_count = 1

   If (m_in.Eq.0._dp) Cycle ! massless particle

   alpha_3 = AlphaSDR(m_in, mGlu, mSup, mSdown) ! needed for QCD corrections
   !------------------
   ! into leptons
   !------------------
   Do i2 = 1,n_l
    Do i3 = i2,n_l
     If ((Abs(cpl_LLS0_L(i2,i3,i1))+Abs(cpl_LLS0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i3), cpl_LLS0_L(i2,i3,i1) &
                             &, cpl_LLS0_R(i2,i3,i1), gam )
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = id_l(i2)
      S0(i1)%id2(i_count,2) = id_l(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = id_l(i2) + 1
       S0(i1)%id2(i_count,2) = id_l(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do

   !------------------
   ! into d-quarks
   !------------------
   Do i2 = 1, n_d
    Do i3 = i2,n_d
     If ((Abs(cpl_DDS0_L(i2,i3,i1))+Abs(cpl_DDS0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i3), cpl_DDS0_L(i2,i3,i1) &
                             &, cpl_DDS0_R(i2,i3,i1), gam )
      gam = 3._dp * gam * FFqcd(mf_d(i2), m_in, alpha_3)
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = id_d(i2)
      S0(i1)%id2(i_count,2) = id_d(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = id_d(i2) + 1
       S0(i1)%id2(i_count,2) = id_d(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into u-quarks
   !------------------
   Do i2 = 1,n_u
    Do i3 = i2,n_u
     If ((Abs(cpl_UUS0_L(i2,i3,i1))+Abs(cpl_UUS0_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i3), cpl_UUS0_L(i2,i3,i1) &
                             &, cpl_UUS0_R(i2,i3,i1), gam )
      gam = 3._dp * gam * FFqcd(mf_u(i2), m_in, alpha_3)
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = id_u(i2)
      S0(i1)%id2(i_count,2) = id_u(i3) + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = id_u(i2) + 1
       S0(i1)%id2(i_count,2) = id_u(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into sleptons
   !------------------
   Do i2 = 1,n_sle
    Do i3 = i2,n_sle
     If (Abs(cpl_S0SlSl(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSlepton(i3)   &
                             &, cpl_S0SlSl(i1,i2,i3), gam )
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = Slept(i2)%id
      S0(i1)%id2(i_count,2) = Slept(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = Slept(i2)%id + 1
       S0(i1)%id2(i_count,2) = Slept(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into sneutrinos
   !------------------
   Do i2 = 1, n_snu
    Do i3 = i2, n_snu
     If (Abs(cpl_S0SnSn(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSneutrino(i2), mSneutrino(i3)   &
                             &, cpl_S0SnSn(i1,i2,i3), gam )
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = Sneut(i2)%id
      S0(i1)%id2(i_count,2) = Sneut(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = Sneut(i2)%id + 1
       S0(i1)%id2(i_count,2) = Sneut(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into d-squarks
   !------------------
   Do i2 = 1, n_sd
    Do i3 = i2, n_sd
     If (Abs(cpl_S0SdSd(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSdown(i3)   &
                             &, cpl_S0SdSd(i1,i2,i3), gam )
      S0(i1)%gi2(i_count) = 3._dp * gam
      S0(i1)%id2(i_count,1) = Sdown(i2)%id
      S0(i1)%id2(i_count,2) = Sdown(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = 3._dp * gam
       S0(i1)%id2(i_count,1) = Sdown(i2)%id + 1
       S0(i1)%id2(i_count,2) = Sdown(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !------------------
   ! into u-squarks
   !------------------
   Do i2 = 1, n_su
    Do i3 = i2, n_su
     If (Abs(cpl_S0SuSu(i1,i2,i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSup(i2), mSup(i3)   &
                             &, cpl_S0SuSu(i1,i2,i3), gam )
      S0(i1)%gi2(i_count) = 3._dp * gam
      S0(i1)%id2(i_count,1) = Sup(i2)%id
      S0(i1)%id2(i_count,2) = Sup(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = 3._dp * gam
       S0(i1)%id2(i_count,1) = Sup(i2)%id + 1
       S0(i1)%id2(i_count,2) = Sup(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do
   !-------------
   ! Neutralinos
   !-------------
   Do i2 = 1, n_n
    Do i3 = i2, n_n
     Call ScalarToTwoFermions(m_in, mN(i2), mN(i3), cpl_NNS0_L(i2,i3, i1) &
                             &, cpl_NNS0_R(i2,i3, i1), gam )
     If (i2.Eq.i3) gam = 0.5_dp * gam ! Majorana
     S0(i1)%gi2(i_count) = gam
     S0(i1)%id2(i_count,1) = Chi0(i2)%id
     S0(i1)%id2(i_count,2) = Chi0(i3)%id 
     i_count = i_count + 1
    End Do
   End Do
   !-------------
   ! Charginos
   !-------------
   Do i2 = 1, n_c
    Do i3 = i2, n_c
     Call ScalarToTwoFermions(m_in, mC(i2), mC(i3), cpl_CCS0_L(i2,i3, i1) &
                             &, cpl_CCS0_R(i2,i3, i1), gam )
     S0(i1)%gi2(i_count) = gam
     S0(i1)%id2(i_count,1) = ChiPm(i2)%id
     S0(i1)%id2(i_count,2) = ChiPm(i3)%id + 1
     i_count = i_count + 1
     If (i2.Ne.i3) Then
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = ChiPm(i2)%id + 1
      S0(i1)%id2(i_count,2) = ChiPm(i3)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !------
   ! Z Z
   !------
   Do i2=1,n_Z
    coupC = cpl_S0ZZ(i1,i2)
    Call ScalarToTwoVectorBosons(m_in, mZ(i2), mZ(i2), coupC, gam )
    gam = 0.5_dp * gam ! identical particles 
    S0(i1)%gi2(i_count) = gam
    S0(i1)%id2(i_count,1) = id_Z(i2)
    S0(i1)%id2(i_count,2) = id_Z(i2)
    i_count = i_count + 1
   End Do
   !------
   ! W W
   !------
   Do i2=1,n_W
    coupC = cpl_S0WW(i1,i2)
    Call ScalarToTwoVectorBosons(m_in, mW(i2), mW(i2), coupC, gam)
    S0(i1)%gi2(i_count) = gam
    S0(i1)%id2(i_count,1) = id_W(i2)
    S0(i1)%id2(i_count,2) = id_W(i2) + 1
    i_count = i_count + 1
   End Do
   !-------------------
   ! pseudoscalars + Z
   !-------------------
   Do i2 = 2,n_P0
    Do i3=1,n_Z
     coupC = cpl_P0S0Z(i2, i1, i3)
     Call ScalarToScalarVectorBoson(m_in, mP0(i2), mZ(i3), coupC, gam)
     S0(i1)%gi2(i_count) = gam
     S0(i1)%id2(i_count,1) = P0(i2)%id
     S0(i1)%id2(i_count,2) = id_Z(i3)
     i_count = i_count + 1
    End Do
   End Do
   !-------------------
   ! 2 pseudoscalars 
   !-------------------
   Do i2 = 2,n_P0
    Do i3 = i2,n_P0
     coupC = cpl_P0S03(i2, i3, i1)
     Call ScalarToTwoScalars(m_in, mP0(i2), mP0(i3), coupC, gam )
     If (i2.Eq.i3) gam = 0.5_dp * gam ! identical particles 
     S0(i1)%gi2(i_count) = gam
     S0(i1)%id2(i_count,1) = P0(i2)%id
     S0(i1)%id2(i_count,2) = P0(i3)%id
     i_count = i_count + 1
    End Do
   End Do
   !-------------------
   ! 2 scalars 
   !-------------------
   Do i2 = 1,i1-1
    Do i3 = i2,i1-1
     coupC = cpl_S03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mS0(i2), mS0(i3), coupC, gam )
     If (i2.Eq.i3) gam = 0.5_dp * gam ! identical particles 
     S0(i1)%gi2(i_count) = gam
     S0(i1)%id2(i_count,1) = S0(i2)%id
     S0(i1)%id2(i_count,2) = S0(i3)%id
     i_count = i_count + 1
    End Do
   End Do
   !--------------------
   ! charged scalar + W
   !--------------------
   Do i2 = 2,n_Spm
    Do i3=1,n_W
     coupC = cpl_SmpS0W(i2, i1, i3)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mW(i3), coupC, gam )
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = SPm(i2)%id
      S0(i1)%id2(i_count,2) = id_W(i3) + 1
      i_count = i_count + 1
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = SPm(i2)%id + 1
      S0(i1)%id2(i_count,2) = id_W(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !-------------------
   ! 2 charged scalars
   !-------------------
   Do i2 = 2,n_Spm
    Do i3 = i2,n_Spm
     coupC = cpl_SmpS03(i2, i3, i1)
     If (Abs(CoupC).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSpm(i2), mSpm(i3), coupC, gam)
      S0(i1)%gi2(i_count) = gam
      S0(i1)%id2(i_count,1) = SPm(i2)%id
      S0(i1)%id2(i_count,2) = SPm(i3)%id + 1
      i_count = i_count + 1
      If (i2.Ne.i3) Then
       S0(i1)%gi2(i_count) = gam
       S0(i1)%id2(i_count,1) = SPm(i2)%id + 1
       S0(i1)%id2(i_count,2) = SPm(i3)%id
       i_count = i_count + 1
      End If
     End If
    End Do
   End Do

   !--------------
   ! two gluons
   !--------------
   S0(i1)%gi2(i_count) = G_F * m_in**3 * oosqrt2 * oo36pi3 *Abs(cpl_GlGlS0(i1))**2
   S0(i1)%id2(i_count,1) = id_gl
   S0(i1)%id2(i_count,2) = id_gl
   i_count = i_count + 1
   !--------------
   ! two photons
   !--------------
   S0(i1)%gi2(i_count) = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cpl_GGS0(i1))**2
   S0(i1)%id2(i_count,1) = id_ph
   S0(i1)%id2(i_count,2) = id_ph
   i_count = i_count + 1

   !-----------------------
   ! W W^*, 3-body decays
   !-----------------------
   S0(i1)%gi3 = 0._dp
   i_count = 1
   Do i2=1,n_W
    Call  ScalarToVectorbosonsVR(m_in, mW(i2), cpl_S0WWvirt(i1,i2), gam )
    S0(i1)%gi3(i_count) = 0.5_dp * gam ! formula is for sum over charges
    S0(i1)%id3(i_count,1) = id_W(i2)
    S0(i1)%id3(i_count,2) = id_W(i2) + 1
    i_count = i_count + 1
    S0(i1)%gi3(i_count) = 0.5_dp * gam ! formula is for sum over charges
    S0(i1)%id3(i_count,1) = id_W(i2) + 1
    S0(i1)%id3(i_count,2) = id_W(i2)
    i_count = i_count + 1
   End Do   

   !-----------------------
   ! Z Z^*
   !-----------------------
   Do i2=1,n_W
    Call  ScalarToVectorbosonsVR(m_in, mZ(i2), cpl_S0ZZvirt(i1,i2), gam )
    S0(i1)%gi3(i_count) = gam 
    S0(i1)%id3(i_count,1) = id_Z(i2)
    S0(i1)%id3(i_count,2) = id_Z(i2)
    i_count = i_count + 1
   End Do   

   S0(i1)%g = Sum(S0(i1)%gi2) + Sum(S0(i1)%gi3)
   If (S0(i1)%g.Ne.0._dp) then
    S0(i1)%bi2 = S0(i1)%gi2 / S0(i1)%g
    S0(i1)%bi3 = S0(i1)%gi3 / S0(i1)%g
   End If

  End Do ! i1
 
  Iname = Iname - 1

 Contains

  Real(dp) Function FFqcd(mf, mA, alpha_s)
  Implicit None
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   If (ratio.Ge.0.5_dp) Return ! decay is kinematically forbitten

   If (ratio.Ge.0.495_dp) Return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (3._dp + 34._dp * beta2 - 13._dp * beta**4) / (16._dp * beta**3)  &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta

    fac = fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  End  Function FFqcd

 End Subroutine ScalarTwoBodyDecays


 Subroutine SfermionTwoBodyDecays(i_in, n_Sf, n_f, id_f, n_n, n_g, n_c, n_fp &
          & , id_fp, n_W, id_W, n_Sfp, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav &
          & , Sf, mf, mfp, Chi0, c_FNSf_L, c_FNSf_R, ChiPm, c_CFpSf_L        &
          & , c_CFpSf_R, Sfp, mW, c_SfSfpW, mZ, c_SfSfZ, Spm, c_SmpSfSfp     &
          & , P0, c_P0SfSf, S0, c_S0SfSf, m_grav, F_eff, c_GraFSf_L          &
          & , c_GraFSf_R, k_neut, Glu, c_GQSq_L, c_GQSq_R)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of sfermions:
 ! input:
 !  i_in ................. specifies the decaying sfermion. The decays of
 !                all sfermions are calculated if n_in < 0. In the case of
 !                generation diagonal models, the sfermions are ordered
 !                according to their generation
 !  Sf(i) ............... sfermion masses + decay information
 !  mf(i) ............... the corresponding fermion masses
 !  mfp(i) .............. the corresponding fermion' masses
 !  Chi0(i) ............. neutralino masses + decay information
 !  c_FNSf_L(i,j,k) ... left fermion-neutralino-sfermion coupling
 !  c_FNSf_R(i,j,k) ... right fermion-neutralino-sfermion coupling
 !  ChiPm(i) ............ chargino masses + decay information
 !  c_CFpSf_L(i,j,k) .. left chargino fermion' sfermion coupling
 !  c_CFpSf_R(i,j,k) .. right chargino fermion' sfermion coupling
 !  Sfp(i) .............. the corresponding sfermion' masses
 !  mW(i) ............... mass of the W-boson
 !  c_SfSfpW(i,j,k) ... coupling sfermion-sfermion'-W
 !  mZ(i) ............... mass of the Z-boson
 !  c_SfSfZ(i,j,k) .... coupling sfermion-sfermion-Z
 !  Spm(i) .............. masses of charged scalars + decay information
 !  c_SmpSfSfp(i,j,k) . charged scalar - sfermion - sfermion' coupling
 !  P0(i) ............... pseudoscalar masses + decay information
 !  c_P0SfSf(i,j,k) ... pseudoscalar - sfermion - sfermion coupling
 !  S0(i) ............... scalar masses + decay information
 !  c_S0SfSf(i,j,k) ... scalar - sfermion - sfermion coupling
 !  k_neut .............. if =1 .... summing over fermions in the neutralino
 !                                   final states
 !                        if =2 .... summing over fermions in the chargino
 !                                   final states
 !                        if =3 .... summing over fermions in the chargino
 !                                   and neutralino final states
 !  Glu ................. Gluino mass + decay information, optional
 !  c_GQSq_L(i,j) ... left gluino quark squark coupling, optional
 !  c_GQSq_R(i,j) ... right gluino quark squark coupling, optional
 ! output: 
 ! written by Werner Porod, 16.04.2001
 !  19.04.2001: adding interface for decay into gluinos
 !  18.09.2010: switiching to new variable types for SUSY + Higgs particles
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut, n_f, n_n, n_g, n_c, n_fp, n_W, n_Sfp   &
        & , n_Spm, n_Z, n_P0, n_S0, n_Sf, id_f(:), id_fp(:), id_W(:), id_Z(:) &
        & , id_grav 
  Real(dp), Intent(in) :: mf(:), mfp(:), mW(:), mZ(:), m_grav, F_eff
  Type(particle23), Intent(in), Optional :: Glu
  Complex(dp), Intent(in) :: c_FNSf_L(:,:,:), c_FNSf_R(:,:,:)    &
                              & , c_CFpSf_L(:,:,:), c_CFpSf_R(:,:,:)
  Complex(dp), Intent(in) :: c_SfSfpW(:,:,:), c_SmpSfSfp(:,:,:)    &
             & , c_SfSfZ(:,:,:), c_P0SfSf(:,:,:), c_S0SfSf(:,:,:)  &
             & , c_GraFSf_L(:,:), c_GraFSf_R(:,:)
  Complex(dp), Intent(in), Optional :: c_GQSq_L(:,:), c_GQSq_R(:,:)
  Type(particle2), intent(in) :: Sfp(:), Spm(:), P0(:)
  Type(particle23), intent(in) :: Chi0(:), ChiPm(:), S0(:)
  Type(particle2), intent(inout) :: Sf(:)

  Integer :: i1, i2, i3, i_gen, i_start, i_end, i_count
  Real(dp) :: mSf(n_Sf), gam, m_in, mN(n_n), mC(n_c), mG, mSfp(n_sfp)  &
      & , mSpm(n_Spm), mP0(n_P0), mS0(n_S0)
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SfermionTwoBodyDecays'

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_sf
   Do i1=1,n_sf
    Sf(i1)%gi2 = 0._dp
    Sf(i1)%bi2 = 0._dp
   End Do

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_sf) ) Then 
   i_start = i_in 
   i_end = i_in
   Sf(i_in)%gi2 = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_sf) = ',i_in,n_sf
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  mSf = Sf%m
  mSfp = Sfp%m
  mN = Chi0%m
  mC = ChiPm%m
  If (Present(Glu)) mG = Glu%m
  mSpm = Spm%m
  mS0 = S0%m
  mP0 = P0%m

  Do i1 = i_start, i_end
   m_in = msf(i1)
   i_count = 1
   Sf(i1)%gi2 = 0._dp
   If (n_sf.eq.3) i_gen = i1
   If (n_sf.Eq.6) i_gen = (i1+1)/2
   !---------------------------------------------------------------------
   ! into neutralinos, if k_neut=1 or k_neut=3 summing over all fermions
   !---------------------------------------------------------------------
   Do i2 = 1, n_n
    Do i3 = 1,n_f
     If ((Abs(c_FNSf_L(i3,i2,i1))+Abs(c_FNSf_R(i3,i2,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf(i3), mN(i2) &
            & , c_FNSf_L(i3,i2,i1), c_FNSf_R(i3,i2,i1), gam)
      If ((k_neut.Eq.1).Or.(k_neut.Eq.3)) Then
       Sf(i1)%gi2(i_count) = Sf(i1)%gi2(i_count) + gam
       Sf(i1)%id2(i_count,1) = Chi0(i2)%id
       Sf(i1)%id2(i_count,2) = id_f(i_gen)
      Else
       Sf(i1)%gi2(i_count) = gam
       Sf(i1)%id2(i_count,1) = Chi0(i2)%id
       Sf(i1)%id2(i_count,2) = id_f(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
    If ((k_neut.Eq.1).Or.(k_neut.Eq.3))  i_count = i_count + 1
   End Do
   !--------------------------------------------------------------------
   ! into charginos, if k_neut=2 or k_neut=3 summing over all fermions
   !--------------------------------------------------------------------
   Do i2 = 1, n_c
    Do i3 = 1,n_fp
     If ((Abs(c_CFpSf_L(i2,i3,i1))+Abs(c_CFpSf_R(i2,i3,i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mfp(i3), mC(i2) &
           & , c_CFpSf_L(i2, i3, i1), c_CFpSf_R(i2, i3, i1), gam)
      If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) Then
       Sf(i1)%gi2(i_count) = Sf(i1)%gi2(i_count) + gam
       Sf(i1)%id2(i_count,1) = ChiPm(i2)%id
       Sf(i1)%id2(i_count,2) = id_fp(i_gen)
      Else
       Sf(i1)%gi2(i_count) = gam
       Sf(i1)%id2(i_count,1) = ChiPm(i2)%id
       Sf(i1)%id2(i_count,2) = id_fp(i3)
       i_count = i_count + 1
      End If
     End If
    End Do
    If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) i_count = i_count + 1
   End Do
   !----------------------------------------------------
   ! into gluino, if k_neut=3 summing over all quarks
   !----------------------------------------------------
   If (Present(Glu).And.Present(c_GQSq_L).And.Present(c_GQSq_R) ) Then
    Do i_gen = 1,n_f
     If ((Abs(c_GQSq_L(i_gen, i1))+Abs(c_GQSq_R(i_gen, i1))).ne.0._dp) then
      Call ScalarToTwoFermions(m_in, mf(i_gen), mG, c_GQSq_L(i_gen, i1) &
             & , c_GQSq_R(i_gen, i1), gam)
      If (k_neut.Eq.3) Then
       Sf(i1)%gi2(i_count) = Sf(i1)%gi2(i_count) + 16._dp * gam / 3._dp ! Colour factor
       Sf(i1)%id2(i_count,1) = Glu%id
       Sf(i1)%id2(i_count,2) = id_f(1)
      Else
       Sf(i1)%gi2(i_count) = 16._dp * gam / 3._dp ! Colour factor 
       Sf(i1)%id2(i_count,1) = Glu%id
       Sf(i1)%id2(i_count,2) = id_f(i_gen)
       i_count = i_count + 1
      End If
     End If
    End Do
    If (k_neut.Eq.3) i_count = i_count + 1
   End If
   !-----------------
   ! W-boson
   !-----------------
   Do i2=1,n_sfp
    Do i3=1,n_W
     If (Abs(c_SfSfpW(i1, i2, i3)).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in,mSfp(i2),mW(i3) &
             & ,c_SfSfpW(i1, i2, i3), gam)
      Sf(i1)%gi2(i_count) = gam 
      Sf(i1)%id2(i_count,1) = Sfp(i2)%id
      Sf(i1)%id2(i_count,2) = id_W(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do
   !-----------------
   ! charged scalar
   !-----------------
   Do i2=2,n_Spm
    Do i3=1,n_sfp
     If (Abs(c_SmpSfSfp(i2, i1, i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSfp(i3), mSpm(i2) &
              & , c_SmpSfSfp(i2, i1, i3), gam)
      Sf(i1)%gi2(i_count) = gam 
      Sf(i1)%id2(i_count,1) = Sfp(i3)%id
      Sf(i1)%id2(i_count,2) = Spm(i2)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !-----------------
   ! Z-boson
   !-----------------
   Do i2=1,i1-1
    Do i3=1,n_Z
     If (Abs(c_SfSfZ(i1, i2, i3)).ne.0._dp) then
      Call ScalarToScalarVectorBoson(m_in,mSf(i2),mZ(i3) &
             & ,c_SfSfZ(i1, i2, i3), gam)
      Sf(i1)%gi2(i_count) = gam 
      Sf(i1)%id2(i_count,1) = Sf(i2)%id
      Sf(i1)%id2(i_count,2) = id_Z(i3)
      i_count = i_count + 1
     End If
    End Do
   End Do
   !-----------------
   ! pseudoscalar
   !-----------------
   Do i2=2,n_P0
    Do i3=1,i1-1
     If (Abs(c_P0SfSf(i2, i1, i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSf(i3), mP0(i2) &
             & , c_P0SfSf(i2, i1, i3), gam)
      Sf(i1)%gi2(i_count) = gam 
      Sf(i1)%id2(i_count,1) = Sf(i3)%id
      Sf(i1)%id2(i_count,2) = P0(i2)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !-----------------
   ! scalar
   !-----------------
   Do i2=1,n_S0
    Do i3=1,i1-1
     If (Abs(c_S0SfSf(i2, i1, i3)).ne.0._dp) then
      Call ScalarToTwoScalars(m_in, mSf(i3), mS0(i2) &
             & , c_S0SfSf(i2, i1, i3), gam)
      Sf(i1)%gi2(i_count) = gam 
      Sf(i1)%id2(i_count,1) = Sf(i3)%id
      Sf(i1)%id2(i_count,2) = S0(i2)%id
      i_count = i_count + 1
     End If
    End Do
   End Do
   !-----------------------
   ! gravitino
   !-----------------------
   Do i2=1,n_f
    If ((Abs(c_GraFSf_L(i2,i1))+Abs(c_GraFSf_R(i2,i1))).ne.0._dp) then
     Call ScalarToFermionGravitino(m_in,mf(i2),m_grav,F_eff &
             & , c_GraFSf_L(i2,i1), c_GraFSf_R(i2,i1), gam)
     Sf(i1)%gi2(i_count) = gam 
     Sf(i1)%id2(i_count,1) = id_grav
     Sf(i1)%id2(i_count,2) = id_f(i2)
     i_count = i_count + 1
    End If
   End Do

   Sf(i1)%g = Sum(Sf(i1)%gi2)
   If (Sf(i1)%g.Ne.0._dp) Sf(i1)%bi2 = Sf(i1)%gi2 / Sf(i1)%g

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine SfermionTwoBodyDecays



 Subroutine SfermionTwoBodyDecays_old(i_in, mSf, mf, mfp                    &
          &, mN, c_FNSf_L, c_FNSf_R, mC, c_CFpSf_L, c_CFpSf_R   &
          &, mSfp, mW, c_SfSfpW, mZ, c_SfSfZ, mSpm, c_SmpSfSfp    &
          &, mP0, c_P0SfSf, mS0, c_S0SfSf                           &
          &, k_neut, GenerationMixing                                   &
          &, gP, gT, BR                                                 &
          &, mG, c_GQSq_L, c_GQSq_R)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of sfermions:
 ! input:
 !  i_in ................. specifies the decaying sfermion. The decays of
 !                all sfermions are calculated if n_in < 0. In the case of
 !                generation diagonal models, the sfermions are ordered
 !                according to their generation
 !  mSf(i) .............. sfermion masses
 !  mf(i) ............... the corresponding fermion masses
 !  mfp(i) .............. the corresponding fermion' masses
 !  mN(i) ............... neutralino masses
 !  c_FNSf_L(i,j,k) ... left fermion-neutralino-sfermion coupling
 !  c_FNSf_R(i,j,k) ... right fermion-neutralino-sfermion coupling
 !  mC(i) ............... chargino masses
 !  c_CFpSf_L(i,j,k) ... left chargino fermion' sfermion coupling
 !  c_CFpSf_R(i,j,k) ... right chargino fermion' sfermion coupling
 !  mSfp(i) .............. the corresponding sfermion' masses
 !  mW ................... mass of the W-boson
 !  c_SfSfpW(i,j) ...... coupling sfermion-sfermion'-W
 !  mZ ................... mass of the Z-boson
 !  c_SfSfZ(i,j) ....... coupling sfermion-sfermion-Z
 !  mSpm(i) .............. masses of charged scalars
 !  c_SmpSfSfp(i,j,k) .. charged scalar - sfermion - sfermion' coupling
 !  mP0(i) ............... pseudoscalar masses
 !  c_P0SfSf(i,j,k) .... pseudoscalar - sfermion - sfermion coupling
 !  mS0(i) ............... scalar masses
 !  c_S0SfSf(i,j,k) .... scalar - sfermion - sfermion coupling
 !  k_neut ................ if =1 .... summing over fermions in the neutralino
 !                                     final states
 !                          if =2 .... summing over fermions in the chargino
 !                                     final states
 !                          if =3 .... summing over fermions in the chargino
 !                                     and neutralino final states
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 !  mG ................... Gluino mass, optional
 !  c_GQSq_L(i,j) ... left gluino quark squark coupling, optional
 !  c_GQSq_R(i,j) ... right gluino quark squark coupling, optional
 ! output: 
 !  depends on the values of k_neut and GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  In addition the variables n_sfer (n_sferp) give the lengths of the
 !  sfermion (sfermions') depending on the sfermion:
 !    n_sfer, n_sferp = 3 for sneutrinos and 6 otherwise 
 !  gP(i,j) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 16.04.2001
 !  19.04.2001: adding interface for decay into gluinos
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut
  Real(dp), Intent(in) :: mSf(:), mf(3), mfp(3), mN(:), mC(:), mSfp(:) &
               & , mW, mZ, mSpm(:), mP0(:), mS0(:)
  Real(dp), Intent(in), Optional :: mG
  Complex(dp), Intent(in) :: c_FNSf_L(:,:,:), c_FNSf_R(:,:,:)    &
                              & , c_CFpSf_L(:,:,:), c_CFpSf_R(:,:,:)
  Complex(dp), Intent(in) :: c_SfSfpW(:,:), c_SmpSfSfp(:,:,:)    &
             & , c_SfSfZ(:,:), c_P0SfSf(:,:,:), c_S0SfSf(:,:,:)
  Complex(dp), Intent(in), Optional :: c_GQSq_L(:,:), c_GQSq_R(:,:)
  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)
  Logical, Intent(in) :: GenerationMixing

  Integer :: i1, i2, n_sfer, n_neut, n_char, i_gen, i_start, i_end, i_count &
         & , n_sferp, n_Spm, i3, n_P0, n_S0
  Real(dp) :: gam, m_in 
  Complex(dp) :: coupLC, coupRC, coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SfermionTwoBodyDecays_old'

  n_sfer = Size(mSf)
  n_sferp = Size(mSfp)
  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_sfer
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_sfer) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_sfer) = ',i_in,n_sfer
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = msf(i1)
   i_count = 1
   If (GenerationMixing) Then
    !---------------------------------------------------------------------
    ! into neutralinos, if k_neut=1 or k_neut=3 summing over all fermions
    !---------------------------------------------------------------------
    Do i2 = 1, n_neut
     Do i_gen = 1,3
      coupLC = c_FNSf_L(i_gen,i2,i1)
      coupRC = c_FNSf_R(i_gen,i2,i1)
      Call ScalarToTwoFermions(m_in, mf(i_gen), mN(i2), coupLC, coupRC, gam)
      If ((k_neut.Eq.1).Or.(k_neut.Eq.3)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If ((k_neut.Eq.1).Or.(k_neut.Eq.3)) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !--------------------------------------------------------------------
    ! into charginos, if k_neut=2 or k_neut=3 summing over all fermions
    !--------------------------------------------------------------------
    Do i2 = 1, n_char
     Do i_gen = 1,3
      coupLC = c_CFpSf_L(i2, i_gen, i1)
      coupRC = c_CFpSf_R(i2, i_gen, i1)

      Call ScalarToTwoFermions(m_in, mfp(i_gen), mC(i2), coupLC, coupRC, gam)
      If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !----------------------------------------------------
    ! into gluino, if k_neut=3 summing over all quarks
    !----------------------------------------------------
    If (Present(mG).And.Present(c_GQSq_L).And.Present(c_GQSq_R) ) Then
     Do i_gen = 1,3
      coupLC = c_GQSq_L(i_gen, i1)
      coupRC = c_GQSq_R(i_gen, i1)
      Call ScalarToTwoFermions(m_in, mf(i_gen), mG, coupLC, coupRC, gam)
      If (k_neut.Eq.3) Then
       gP(i1, i_count) = gP(i1, i_count) + 16._dp * gam / 3._dp ! Colour factor
      Else
       gP(i1, i_count) = 16._dp * gam / 3._dp ! Colour factor 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.3) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End If
    !-----------------
    ! W-boson
    !-----------------
    Do i2=1,n_sferp
     coupC = c_SfSfpW(i1, i2)
     Call ScalarToScalarVectorBoson(m_in,mSfp(i2),mW,coupC,gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !-----------------
    ! charged scalar
    !-----------------
    Do i2=2,n_Spm
     Do i3=1,n_sferp
      coupC = c_SmpSfSfp(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSfp(i3), mSpm(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do
    !-----------------
    ! Z-boson
    !-----------------
    Do i2=1,i1-1
     coupC = c_SfSfZ(i1, i2)
     Call ScalarToScalarVectorBoson(m_in, mSf(i2), mZ, coupC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !-----------------
    ! pseudoscalar
    !-----------------
    Do i2=2,n_P0
     Do i3=1,i1-1
      coupC = c_P0SfSf(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSf(i3), mP0(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do
    !-----------------
    ! scalar
    !-----------------
    Do i2=1,n_S0
     Do i3=1,i1-1
      coupC = c_S0SfSf(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSf(i3), mS0(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    If (n_sfer.Eq.6) Then
     i_gen = (i1+1) / 2
    Else
     i_gen = i1
    Endif
    !--------------------------
    ! into neutralinos
    !--------------------------
    Do i2 = 1,n_neut
     coupLC = c_FNSf_L(i_gen,i2,i1)
     coupRC = c_FNSf_R(i_gen,i2,i1)
     Call ScalarToTwoFermions(m_in, mf(i_gen), mN(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !--------------------------
    ! into charginos
    !--------------------------
    Do i2 = 1, n_char
     coupLC = c_CFpSf_L(i2, i_gen, i1)
     coupRC = c_CFpSf_R(i2, i_gen, i1)
     Call ScalarToTwoFermions(m_in, mfp(i_gen), mC(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !--------------------------
    ! into gluino
    !--------------------------
    If (Present(mG).And.Present(c_GQSq_L).And.Present(c_GQSq_R) ) Then
     coupLC = c_GQSq_L(i_gen, i1)
     coupRC = c_GQSq_R(i_gen, i1)
     Call ScalarToTwoFermions(m_in, mf(i_gen), mG, coupLC, coupRC, gam)
     gP(i1, i_count) = 16._dp * gam / 3._dp ! Colour factor 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End If

    !-----------------
    ! W-boson
    !-----------------
    If (n_sferp.Eq.3) Then
     coupC = c_SfSfpW(i1, i_gen)
     Call ScalarToScalarVectorBoson(m_in,mSfp(i_gen),mW,coupC,gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    Else
     Do i2 =1,2
      coupC = c_SfSfpW(i1, (i_gen-1)*2 + i2)
      Call ScalarToScalarVectorBoson(m_in,mSfp((i_gen-1)*2 + i2),mW,coupC,gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End If   

    !-----------------
    ! charged scalar
    !-----------------
    If (n_sferp.Eq.3) Then
     Do i2 = 2, n_Spm
      coupC = c_SmpSfSfp(i2, i1, i_gen)
      Call ScalarToTwoScalars(m_in, mSfp(i_gen), mSpm(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    Else
     Do i2 =2,n_Spm
      Do i3 =1,2
       coupC = c_SmpSfSfp(i2, i1, (i_gen-1)*2 + i3)
       Call ScalarToTwoScalars(m_in, mSfp((i_gen-1)*2+i3), mSpm(i2), coupC &
                             &, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End Do
     End Do
    End If   

    !-----------------
    ! Z-boson
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     If (i1.Eq.2) Then
      coupC = c_SfSfZ(1, 2)
      Call ScalarToScalarVectorBoson(m_in, mSf(1), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     Else If (i1.Eq.4) Then
      coupC = c_SfSfZ(3, 4)
      Call ScalarToScalarVectorBoson(m_in, mSf(3), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     Else If (i1.Eq.6) Then
      coupC = c_SfSfZ(5, 6)
      Call ScalarToScalarVectorBoson(m_in, mSf(5), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End If   
    End If   

    !-----------------
    ! pseudoscalar
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     Do i2=2,n_P0
      If (i1.Eq.2) Then
       coupC = c_P0SfSf(i2, 1, 2)
       Call ScalarToTwoScalars(m_in, mSf(1), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.4) Then
       coupC = c_P0SfSf(i2, 3, 4)
       Call ScalarToTwoScalars(m_in, mSf(3), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.6) Then
       coupC = c_P0SfSf(i2, 5, 6)
       Call ScalarToTwoScalars(m_in, mSf(5), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If   
     End Do
    End If   

    !-----------------
    ! scalar
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     Do i2=1,n_S0
      If (i1.Eq.2) Then
       coupC = c_S0SfSf(i2, 1, 2)
       Call ScalarToTwoScalars(m_in, mSf(1), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.4) Then
       coupC = c_S0SfSf(i2, 3, 4)
       Call ScalarToTwoScalars(m_in, mSf(3), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.6) Then
       coupC = c_S0SfSf(i2, 5, 6)
       Call ScalarToTwoScalars(m_in, mSf(5), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If   
     End Do
    End If   

   End If    ! GenerationMixing


   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine SfermionTwoBodyDecays_old

 Subroutine NeutralVectorTwoBodyDecays(mV, mf_d, c_DDV_L, c_DDV_R, mf_u       &
    & , c_UUV_L, c_UUV_R, mf_l, c_LLV_L, c_LLV_R, mf_nu, c_NuNuV_L, c_NuNuV_R & 
    & , mSdown, c_SdSdV, mSup, c_SuSuV, mSlept, c_SlSlV, mSneut, c_SnSnV      &
    & , mN, c_NNV_L, c_NNV_R, mC, c_CCV_L, c_CCV_R, mP0, mS0, c_P0S0V         &
    & , mSpm, c_SmSpV, GenerationMixing, gP, gT, BR)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of a heavy neutral vectorboson 
 ! input:
 !  mV .................. mass of the vector boson
 !  mf_d(i) ............. d-quark masses
 !  c_DDV_L(i,j) ........ left coupling d-quark d-quark vectorboson
 !  c_DDV_R(i,j) ........ right coupling d-quark d-quark vectorboson
 !  mf_u(i) ............. u-quark masses
 !  c_UUV_L(i,j) ........ left coupling u-quark u-quark vectorboson
 !  c_UUV_R(i,j) ........ right coupling u-quark u-quark vectorboson
 !  mf_l(i) ............. lepton masses
 !  c_LLV_L(i,j) ........ left coupling lepton lepton vectorboson
 !  c_LLV_R(i,j) ........ right coupling lepton lepton vectorboson
 !  mf_nu(i) ............ neutrino masses
 !  c_NuNuV_L(i,j) ...... left coupling neutrino neutrino vectorboson
 !  c_NuNuV_R(i,j) ...... right coupling neutrino neutrino vectorboson
 !  mSdown(i) ........... d-squark masses
 !  c_SdSdV(i,j) ........ coupling vectorboson d-squark d-squark
 !  mSup(i) ............. u-squark masses
 !  c_SuSuV(i,j) ........ coupling vectorboson u-squark u-squark
 !  mSlept(i) ........... slepton masses
 !  c_SlSlV(i,j) ........ coupling vectorboson slepton slepton
 !  mSneut(i) ........... sneutrino masses
 !  c_SnSnV(i,j) ........ coupling vectorboson sneutrino sneutrino
 !  mN(i) ............... neutralino masses
 !  c_NNV_L(i,j) ........ left neutralino-neutralino-vectorboson coupling
 !  c_NNV_R(i,j) ........ right neutralino-neutralino-vectorboson coupling
 !  mC(i) ............... chargino masses
 !  c_CCV_L(i,j) ........ left chargino-chargino-vectorboson coupling
 !  c_CCV_R(i,j) ........ right chargino-chargino-vectorboson coupling
 !  mP0(i) .............. pseudoscalar masses
 !  mS0(i) .............. scalar masses
 !  c_P0S0V(i,j) ........ pseudoscalar - pseudoscalar - vectorboson coupling
 !  mSpm(i) ............. masses of charged scalars
 !  c_SmSpV(i,j) ........ charged scalar - charged scalar - vectorboson coupling
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output: 
 !  depends on the value of  GenerationMixing and also on the
 !  lengths of mNu, mSneutm mN, mC, mSpm, mP0, and mS0 which are measured 
 !  by n_nu, n_snu, n_neut, n_char, n_Spm, n_P0 and n_S0, respectively
 !   (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 05.09.2010
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mV, mf_l(3), mf_d(3), mf_u(3), mf_nu(:), mSlept(:) &
     & , mSneut(:), mSdown(6), mSup(6), mN(:), mC(:), mSpm(:), mP0(:), mS0(:)
  Real(dp), Intent(in) :: c_P0S0V(:,:)
  Complex(dp), Intent(in) :: c_LLV_L(:,:), c_LLV_R(:,:), c_DDV_L(:,:)          &
     & , c_DDV_R(:,:), c_UUV_L(:,:), c_UUV_R(:,:), c_SlSlV(:,:), c_SnSnV(:,:)  &
     & , c_SdSdV(:,:), c_SuSuV(:,:), c_NNV_L(:,:), c_NNV_R(:,:), c_CCV_L(:,:)  &
     & , c_CCV_R(:,:), c_SmSpV(:,:), c_NuNuV_L(:,:), c_NuNuV_R(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(inout) :: gP(:), gT
  Real(dp), Optional, Intent(inout) :: BR(:)

  Integer :: i2, i3, i4, n_neut, n_char, i_count, n_Spm, n_P0, n_S0, n_nu, n_snu
  Real(dp) :: m1, m2
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralVectorTwoBodyDecays'

  n_nu = Size(mf_nu)
  n_snu = Size(mSneut)
  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  i_count = 1

  If (GenerationMixing) Then
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call VectorbosonToTwoFermions(mV, mf_d(i2), mf_d(i3), 3, c_DDV_L(i2,i3) &
                             &, c_DDV_R(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call VectorbosonToTwoFermions(mV, mf_u(i2), mf_u(i3), 3, c_UUV_L(i2,i3) &
                             &, c_UUV_R(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1,5-n_char 
     Do i3 = i2,5-n_char
      Call VectorbosonToTwoFermions(mV, mf_l(i2), mf_l(i3), 1, c_LLV_L(i2,i3) &
                             &, c_LLV_R(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1,n_nu
     Do i3 = i2,n_nu
      Call VectorbosonToTwoFermions(mV, mf_nu(i2), mf_nu(i3), 1        &
                   & , c_NuNuV_L(i2,i3), c_NuNuV_R(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call VectorbosonToTwoScalars(mV, mSdown(i2), mSdown(i3), 3 &
                            & , c_SdSdV(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call VectorbosonToTwoScalars(mV, mSup(i2), mSup(i3), 3, c_SuSuV(i2,i3)   &
                             &, gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,2*(5-n_char)
     Do i3 = i2,2*(5-n_char)
      Call VectorbosonToTwoScalars(mV, mSlept(i2), mSlept(i3), 1   &
                             &, c_SlSlV(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into sneutrinos
    !------------------
    Do i2 = 1, n_snu
     Do i3 = i2,n_snu
      Call VectorbosonToTwoScalars(mV, mSneut(i2), mSneut(i3), 1   &
                             &, c_SnSnV(i2,i3), gP(i_count) )
      If (i2.Ne.i3) Then
       gP(i_count+1) = gP(i_count) ! charge
       gT = gT + 2._dp * gP(i_count)
       i_count = i_count + 2
      Else
       gT = gT + gP(i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Call VectorbosonToTwoFermions(mV, mf_d(i2), mf_d(i2), 3, c_DDV_L(i2,i2) &
                            &, c_DDV_R(i2,i2), gP(i_count) )
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Call VectorbosonToTwoFermions(mV, mf_u(i2), mf_u(i2), 3, c_UUV_L(i2,i2) &
                            &, c_UUV_R(i2,i2), gP(i_count) )
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1, 5-n_char
     Call VectorbosonToTwoFermions(mV, mf_l(i2), mf_l(i2), 1, c_LLV_L(i2,i2) &
                            &, c_LLV_R(i2,i2), gP(i_count) )
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into neutrinos
    !------------------
    Do i2 = 1, n_nu
     Call VectorbosonToTwoFermions(mV, mf_nu(i2), mf_nu(i2), 1        &
                   &, c_NuNuV_L(i2,i2), c_NuNUV_R(i2,i2), gP(i_count) )
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do

    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = c_SdSdV((i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSdown((i2-1)*2+i3)
       m2 = mSdown((i2-1)*2+i4)
       Call VectorbosonToTwoScalars(mV, m1, m2, 3, coupC, gP(i_count) )
       If (i3.Eq.i4) Then
        gT = gT + gP(i_count)
        i_count = i_count + 1
       Else
        gP(i_count+1) = gP(i_count) ! colour + charge
        gT = gT + 2._dp * gP(i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = c_SuSuV((i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSup((i2-1)*2+i3)
       m2 = mSup((i2-1)*2+i4)
       Call VectorbosonToTwoScalars(mV, m1, m2, 3, coupC, gP(i_count) )
       If (i3.Eq.i4) Then
        gT = gT + gP(i_count)
        i_count = i_count + 1
       Else
        gP(i_count+1) = gP(i_count) ! colour + charge
        gT = gT + 2._dp * gP(i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do

    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,5-n_char
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = c_SlSlV((i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSlept((i2-1)*2+i3)
       m2 = mSlept((i2-1)*2+i4)
       Call VectorbosonToTwoScalars(mV, m1, m2, 1, coupC, gP(i_count) )
       If (i3.Eq.i4) Then
        gT = gT + gP(i_count)
        i_count = i_count + 1
       Else
        gP(i_count+1) = gP(i_count) ! charge
        gT = gT + 2._dp * gP(i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into sneutrinos
    !------------------
    Do i2 = 1,n_snu
     Call VectorbosonToTwoScalars(mV, mSneut(i2), mSneut(i2), 1   &
                            &, c_SnSnV(i2,i2), gP(i_count) )
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do

   End If    ! GenerationMixing

   !-------------
   ! Neutralinos
   !-------------
   Do i2 = 1, n_neut
    Do i3 = i2, n_neut
     Call VectorbosonToTwoFermions(mV, mN(i2), mN(i3), 1, c_NNV_L(i2,i3) &
                             &, c_NNV_R(i2,i3), gP(i_count) )
     If (i2.Eq.i3) gP(i_count) = 0.5_dp * gP(i_count) ! Majorana
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do
   End Do
   
   !-------------
   ! Charginos
   !-------------
   Do i2 = 1, n_char
    Do i3 = i2, n_char
     Call VectorbosonToTwoFermions(mV, mC(i2), mC(i3), 1, c_CCV_L(i2,i3) &
                             &, c_CCV_R(i2,i3), gP(i_count) )
     If (i2.Ne.i3) Then
      gP(i_count+1) = gP(i_count) ! charge
      gT = gT + 2._dp * gP(i_count)
      i_count = i_count + 2
     Else
      gT = gT + gP(i_count)
      i_count = i_count + 1
     End If
    End Do
   End Do
   
   !-------------------------
   ! pseudoscalars + scalar 
   !-------------------------
   Do i2 = 2,n_P0
    Do i3 = 1,n_S0
     coupC = c_P0S0V(i2,i3)
     Call VectorbosonToTwoScalars(mV, mP0(i2), mS0(i3), 1, coupC, gP(i_count))
     gT = gT + gP(i_count)
     i_count = i_count + 1
    End Do
   End Do

   !-------------------
   ! 2 charged scalars
   !-------------------
   Do i2 = 2,n_Spm
    Do i3 = i2,n_Spm
     coupC = c_SmSpV(i2, i3)
     Call VectorbosonToTwoScalars(mV, mSpm(i2), mSpm(i3), 1, coupC, gP(i_count) )
     If (i2.Ne.i3) Then
      gP(i_count+1) = gP(i_count) ! charge
      gT = gT + 2._dp * gP(i_count)
      i_count = i_count + 2
     Else
      gT = gT + gP(i_count)
      i_count = i_count + 1
     End If
    End Do
   End Do

   If ((Present(BR)).And.(gT.Eq.0)) Then
    BR = 0._dp
   Else If (Present(BR)) Then
    BR = gP / gT
   End If

  Iname = Iname - 1

 End Subroutine NeutralVectorTwoBodyDecays

  Subroutine Zdecays(gp, g, mZ, mf_l, mf_d, mf_u, BR_ll, BR_inv, BR_dd, BR_uu &
                    &, gT, gP_ll, gP_inv, gP_dd, gP_uu)
  implicit none

   Real(dp), Intent(in) :: gp, g, mf_l(3), mf_d(3), mf_u(3), mZ
   Real(dp), intent(out) :: BR_ll(3), BR_inv, BR_dd(3), BR_uu(2), gT
   Real(dp), Intent(out), optional :: gP_ll(3), gP_inv, gP_dd(3), gP_uu(2)

   Real(dp) :: sinW2, cL, cR, gam, g_ll(3), g_inv, g_dd(3), g_uu(2)

   sinW2 = gp**2/(gp**2 + g**2)
   gT = 0._dp
Write(*,*) gp,g,sinw2
   Call CoupFermionZ(-0.5_dp,-1._dp, g,sinW2,cL,cR)
   Call VectorBosonToTwoFermions(mZ,0._dp,0._dp,1,cL,cR,gam)
   g_LL(1:2) = gam
   Call VectorBosonToTwoFermions(mZ,mf_l(3),mf_l(3),1,cL,cR,gam)
   g_LL(3) = gam
   gT = sum(g_LL)
         
   call CoupFermionZ(0.5_dp,0._dp, g,sinW2,cL,cR)
   Call VectorBosonToTwoFermions(mZ,0._dp,0._dp,1,cL,cR,gam)
   g_inv = 3._dp * gam
   gT = gT + g_inv

   Call CoupFermionZ(-0.5_dp,-1._dp/3._dp, g,sinW2,cL,cR)
   Call VectorBosonToTwoFermions(mZ,0._dp,0._dp,3,cL,cR,gam) 
   g_dd(1:2) = gam
   Call VectorBosonToTwoFermions(mZ,mf_d(3),mf_d(3),3,cL,cR,gam) 
   g_dd(3) = gam
   gT = gT + Sum(g_dd)

   Call CoupFermionZ(0.5_dp,2._dp/3._dp, g,sinW2,cL,cR)
   Call VectorBosonToTwoFermions(mZ,0._dp,0._dp,3,cL,cR,gam)
   g_uu(1) = gam
   Call VectorBosonToTwoFermions(mZ,mf_u(2),mf_u(2),3,cL,cR,gam)
   g_uu(2) = gam
   gT = gT + Sum(g_uu)
 
   BR_ll = g_ll / gT
   BR_inv = g_inv / gT
   BR_dd = g_dd / gT
   BR_uu = g_uu / gT

   Write(*,*) "total",gT
   Write(*,*) " L",g_lL
   Write(*,*) "Nu",g_inv
   Write(*,*) " d",g_dd
   Write(*,*) " u",g_uu
   Write(*,*) "BR"
   Write(*,*) " L",br_lL
   Write(*,*) "Nu",br_inv
   Write(*,*) " d",br_dd
   Write(*,*) " u",br_uu
  If (Present(gP_ll)) gP_ll = g_ll
  If (Present(gP_inv)) gP_inv = g_inv
  If (Present(gP_dd)) gP_dd = g_dd
  If (Present(gP_uu)) gP_uu = g_uu

 End Subroutine Zdecays

End Module  SusyDecays

