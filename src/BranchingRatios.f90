Module BranchingRatios 

! load modules
Use Control
Use StandardModel
Use Chargino3Decays
Use Couplings
Use Neut3Decays
Use Slepton3BodyDecays
Use LoopCouplings
Use SusyMasses
Use SusyDecays
Use Model_Data
Use MSSM_Data
Use Gluino3Decays
Use Stop3BodyDecays
Use LoopFunctions
Use RGEs
! load modules

! interfaces
Interface CalculateBR
 Module Procedure CalculateBR_MSSM, CalculateBR_NMSSM, CalculateBR_RPeps &
      & , CalculateBR_RPlam
End Interface
! interfaces

! private variables
!Real(dp), Private :: erg, factor(3), invfac(2)
!Integer, private :: i2
Contains


 Subroutine CalculateBR_MSSM(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u, n_Z &
     & , id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0, n_p0  &
     & , n_Spm, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, ChiPm, U, V, Chi0 &
     & , N, Sneut, RSneut, Slept, RSlepton, Sup, RSup, Sdown, RSdown, uL_L     &
     & , uL_R, uD_L, uD_R, uU_L, uU_R, S0, RS0, P0, RP0, Spm, RSpm, epsI       &
     & , deltaM, CTBD, fac3, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, vevSM, Fgmsb    &
     & , m32, grav_fac) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 ! 28.02.14: changing units of m32 to GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu, n_sle, n_Sd &
     & , n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm, id_grav, id_gl, id_ph 
  Integer, Intent(in), Dimension(1) :: id_Z, id_W
  Integer, Intent(in), Dimension(3) :: id_nu, id_l, id_d, id_u
  
  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: RP0(2,2), RS0(2,2)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(2,2), U(2,2), V(2,2), N(4,4) &
         & , RSlepton(6,6), Rsneut(3,3), RSup(6,6), RSdown(6,6), uL_L(3,3) &
         & , uL_R(3,3), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu
  Real(dp), Intent(in) :: vevSM(2), Fgmsb, m32, grav_fac
  Logical, Intent(in) :: CTBD
  Type(particle2), Intent(inout) :: Sdown(6), Spm(2), P0(2)
  Type(particle23), Intent(inout) :: Slept(6), Chi0(4), Sup(6), ChiPm(2), Glu &
         & , S0(2), Sneut(3)
   ! contains only the information on the 2-body
  Type(particle2) :: Sup2(6), Slept2(6), Sneut2(3)
                             ! decays of the scalar up
  Type(particle2) :: SMp(2) ! contains the identies of the negative charged ones
  Type(particle23) :: ChiM(2) ! contains the identies of the negative charged ones
                              
  Integer :: i1, i2, k_neut
  Real(dp) :: sinW

  Complex(dp) :: c_SmpSlSn(2,6,3), c_SmpSdSu(2,6,6)   &
      & , c_SmpSnSl(2,3,6), c_SmpSuSd(2,6,6), c_SmpP03(2,2,2)    &
      & , c_SmpP0W(2,2,1), c_SmpS03(2,2,2), c_SmpS0W(2,2,1)          &
      & , c_SmpLNu_L(2,3,3), c_SmpLNu_R(2,3,3), c_SmpDU_L(2,3,3) &
      & , c_SmpDU_R(8,3,3), c_SmpZ(2,2,1), c_DUW(3,3)
  Real(dp) :: c_LLZ_L, c_LLZ_R, c_DDZ_L, c_DDZ_R  &
      & , c_UUZ_L, c_UUZ_R, c_NuNuZ_L, c_NuNuZ_R
  Complex(dp) :: c_CCZ_L(2,2,1), c_CCZ_R(2,2,1)           &
      & , c_NNZ_L(4,4,1), c_NNZ_R(4,4,1), c_NNS0_L(4,4,2)            &
      & , c_NNS0_R(4,4,2), c_NNP0_L(4,4,2), c_NNP0_R(4,4,2) 
  Complex(dp) :: c_GDSd_L(3,6), c_GDSd_R(3,6)          &
      & , c_DNSd_L(3,4,6), c_DNSd_R(3,4,6), c_GUSu_L(3,6)         &
      & , c_GUSu_R(3,6), c_UNSu_L(3,4,6), c_UNSu_R(3,4,6)         &
      & , c_LNSl_L(3,4,6), c_LNSl_R(3,4,6), c_NuNSn_L(3,4,3)      & 
      & , c_NuNSn_R(3,4,3), c_DDP0_L(3,3,2), c_LLP0_L(3,3,2)      &
      & , c_UUP0_L(3,3,2), c_DDP0_R(3,3,2), c_LLP0_R(3,3,2)       &
      & , c_UUP0_R(3,3,2), c_DDS0_L(3,3,2), c_LLS0_L(3,3,2)       &
      & , c_UUS0_L(3,3,2), c_DDS0_R(3,3,2), c_LLS0_R(3,3,2)       &
      & , c_UUS0_R(3,3,2)
  Complex(dp) :: c_CUSd_L(2,3,6), c_CUSd_R(2,3,6)      &
      & , c_CDSu_L(2,3,6), c_CDSu_R(2,3,6), c_CLSn_L(2,3,3)       &
      & , c_CLSn_R(2,3,3), c_CNuSl_L(2,3,6), c_CNuSl_R(2,3,6)
  Complex(dp) :: c_GlGlS0(2), c_GGS0(2),  c_GlGlP0(2), c_GGP0(2)
  Complex(dp) :: c_P0SdSd(2,6,6), c_P0SuSu(2,6,6), c_P0SlSl(2,6,6) &
      & , c_P0SnSn(2,3,3), c_P0S0Z(2,2,1) 
  Real(dp) :: c_P0S03(2,2,2), c_S0WWvirt(2,1), c_S0ZZvirt(2,1), vev
  Complex(dp) :: c_S0SdSd(2,6,6), c_S0SuSu(2,6,6) &
      & , c_S0SlSl(2,6,6), c_S0SnSn(2,3,3), c_LNuW(3,3)
  Real(dp) :: c_S03(2,2,2), c_S0WW(2,1), c_S0ZZ(2,1), c_FFpW
  Complex(dp) :: c_SdSuW(6,6,1), c_SuSdW(6,6,1), c_SlSnW(6,3,1) &
      & , c_SnSlW(3,6,1), c_SdSdZ(6,6,1), c_SlSlZ(6,6,1), c_SnSnZ(3,3,1)     &
      & , c_SuSuZ(6,6,1)
  Complex(dp) :: c_CCP0_L(2,2,2), c_CCP0_R(2,2,2)    &
      & , c_CCS0_L(2,2,2), c_CCS0_R(2,2,2), c_CNW_L(2,4,1)        &
      & , c_CNW_R(2,4,1), c_SmpCN_L(2,2,4), c_SmpCN_R(2,2,4), c_NGP &
      & , c_NGZ(2), c_NGH, c_CGW_L(2), c_CGW_R(2)
  Complex(dp), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R, c_GraLSl_L, c_GraLSl_R
  Complex(dp), Dimension(3,3) :: c_GraNuSn_L, c_GraNuSn_R

  Real(dp) :: tanb, sinW2, cosW, M2_H(2), g_i(3), CosW2
  Complex(dp) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_Q(3,3), M2_U(3,3), Mi(3) &
      & , B, Y_l_h(3,3), Y_d_h(3,3), Y_u_h(3,3), A_l_h(3,3), A_d_h(3,3)       &
      & , A_u_h(3,3), mu_h, cR
  Integer :: kont
  Complex(dp) :: coup, g_u(3), g_d(3), g_l(3), g_c(2), g_sl(6), g_sd(6), g_su(6)
  Complex(dp), Dimension(3,3) :: uD_L_h, uD_R_h, uL_L_h, uL_R_h, uU_L_h, uU_R_h
  Real(dp) :: g_W, g_Hp, m_W(1), m_Z(1), F_eff, gam, mf(3), mf_d2_h(3)   &
      & , mf_u2_h(3), mf_l2_h(3), mGlu, mC(2), mC2(2), mN(4), mN2(4)     &
      & , mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6), mSdown(6)    &
      & , mSdown2(6), mSup(6), mSup2(6), mP0(2), mP02(2), RP0_h(2,2)     &
      & , mS0(2), mS02(2), RS0_h(2,2), mSpm(2), mSpm2(2), sinb, cosb
  Complex(dp) :: U_h(2,2), V_h(2,2), N_h(4,4), Rsneut_h(3,3), RSlepton_h(6,6)  &
              & , RSdown_h(6,6), RSup_h(6,6), RSpm_h(2,2), PhiGlu, yuk, A      &
              & , Rsl(2,2) , Rsd(2,2) , Rsu(2,2) 
  Real(dp) :: fakt, r_T, width
  Real(dp), Parameter :: mf_nu(3)=0._dp, e_d = -1._dp/3._dp, e_u = 2._dp/3._dp
  Complex(dp), Parameter :: g_sd0(6) = 0._dp, g_su0(6) = 0._dp &
      & , g_u0(3) = (1._dp,0._dp), g_d0(3)  = (1._dp,0._dp)
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_MSSM'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)

  m_W = mW
  m_Z = mZ
  Spm(1)%g = gamW
  P0(1)%g = gamZ

  Sup2%m = Sup%m 
  Sup2%id = Sup%id 
  Slept2%m = Slept%m 
  Slept2%id = Slept%id 
  Sneut2%m  = Sneut%m 
  Sneut2%id = Sneut%id 
  ChiM%m = ChiPm%m
  ChiM%id = ChiPm%id + 1
  SMp%m = SPm%m
  SMp%id = SPm%id + 1

  Call AllCouplingsMSSM(gauge, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R, Y_u           &
    & , uU_L, uU_R, vevSM, RSpm, RP0, RS0, U, V, N, mu, PhaseGlu, RSlepton     &
    & , A_l, Rsneut, RSup, A_u, RSdown, A_d                                    &
    & , c_SmpSlSn, c_SmpSdSu, c_SmpSnSl, c_SmpSuSd, c_SmpP03, c_SmpP0W         &
    & , c_SmpS03, c_SmpS0W, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R       &
    & , c_SmpZ(:,:,1), c_DUW, c_LLZ_L, c_LLZ_R, c_DDZ_L, c_DDZ_R               &
    & , c_UUZ_L, c_UUZ_R, c_NuNuZ_L, c_NuNuZ_R, c_CCZ_L(:,:,1), c_CCZ_R(:,:,1) &
    & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_NNS0_L, c_NNS0_R                     &
    & , c_NNP0_L, c_NNP0_R, c_GDSd_L, c_GDSd_R, c_DNSd_L                       &
    & , c_DNSd_R, c_GUSu_L, c_GUSu_R, c_UNSu_L, c_UNSu_R                       &
    & , c_LNSl_L, c_LNSl_R, c_NuNSn_L, c_NuNSn_R, c_DDP0_L                     &
    & , c_LLP0_L, c_UUP0_L, c_DDP0_R, c_LLP0_R, c_UUP0_R                       &
    & , c_DDS0_L, c_LLS0_L, c_UUS0_L, c_DDS0_R, c_LLS0_R                       &
    & , c_UUS0_R, c_CUSd_L, c_CUSd_R, c_CDSu_L, c_CDSu_R                       &
    & , c_CLSn_L, c_CLSn_R, c_CNuSl_L, c_CNuSl_R, c_GlGlS0                     &
    & , c_P0SdSd, c_P0SuSu, c_P0SlSl, c_P0SnSn, c_P0S0Z(:,:,1), c_P0S03        &
    & , c_S0SdSd, c_S0SuSu, c_S0SlSl, c_S0SnSn, c_S03, c_S0WW(:,1)             &
    & , c_S0ZZ(:,1), c_FFpW, c_LNuW, c_SdSuW(:,:,1), c_SuSdW(:,:,1)            &
    & , c_SlSnW(:,:,1), c_SnSlW(:,:,1), c_SdSdZ(:,:,1), c_SlSlZ(:,:,1)         &
    & , c_SnSnZ(:,:,1), c_SuSuZ(:,:,1), c_CCP0_L, c_CCP0_R                     &
    & , c_CCS0_L, c_CCS0_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)                     &
    & , c_SmpCN_L, c_SmpCN_R, c_GraDSd_L, c_GraDSd_R, c_GraUSu_L, c_GraUSu_R   &
    & , c_GraLSl_L, c_GraLSl_R, c_GraNuSn_L, c_GraNuSn_R, GenerationMixing)

  tanb = vevSM(2) / vevSM(1)
  cosb = 1._dp / Sqrt(1._dp + tanb**2)
  sinb = cosb * tanb 
  Call LoopCouplingsMSSM(GenerationMixing, Y_d(3,3), tanb, A_d(3,3), mu &
               & , Sdown%m, Sup%m2, RSup, P0%m2, CKM, c_UNSu_R, c_GUSu_R)

  If (grav_fac.Ne.0._dp) Then
   F_eff = grav_fac*Fgmsb ! effective SUSY breaking scale in case of GSMB
  Else 
   F_eff = Fgmsb
  End If
  !------------------------------------------------
  ! 2-body decays of sneutrinos
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Snu, n_nu, id_nu, n_n, 0, n_c, n_l, id_l   &
          & , n_W, id_W, n_Sle, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sneut2 &
          & , mf_nu, mf_l, Chi0, c_NuNSn_L, c_NuNSn_R, ChiPm, c_CLSn_L        &
          & , c_CLSn_R, Slept2, m_W, c_SnSlW, m_Z, c_SnSnZ, Spm               &
          & , c_SmpSnSl, P0, c_P0SnSn, S0, c_S0SnSn, m32, F_eff               &
          & , c_GraNuSn_L, c_GraNuSn_R, 1)

  Do i1=1,3
   Sneut(i1)%g = Sneut2(i1)%g
   Sneut(i1)%gi2 = Sneut2(i1)%gi2
   Sneut(i1)%bi2 = Sneut2(i1)%bi2
   Sneut(i1)%id2 = Sneut2(i1)%id2
  End Do 

  !------------------------------------------------
  ! 3-body decays of sneutrinos
  !------------------------------------------------
  Call Sneutrino_3body(-1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u, n_c   &
    & , n_n, n_W, id_W, n_Sle, n_Snu, n_Spm, Sneut, mW, mf_l, mf_u, mf_d       &
    & , Chi0, ChiPm, Spm, Slept, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn        &
    & , c_CNuSl_L, c_CNuSl_R, c_LNuW, c_DUW, c_CLSn_L, c_CLSn_R, c_NuNSn_L     &
    & , c_NuNSn_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R                &
    & , GenerationMixing, epsI)

  !------------------------------------------------
  ! 2-body decays of sleptons
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Sle, n_l, id_l, n_n, 0, n_c, n_nu, id_nu   &
          & , n_W, id_W+1, n_Snu, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav       &
          & , Slept2, mf_l, mf_nu, Chi0, c_LNSl_L, c_LNSl_R, ChiM, c_CNuSl_L  &
          & , c_CNuSl_R, Sneut2, m_W, c_SlSnW, m_Z, c_SlSlZ, Smp, c_SmpSlSn   &
          & , P0, c_P0SlSl, S0, c_S0SlSl, m32, F_eff, c_GraLSl_L              &
          & , c_GraLSl_R, 2)

  Do i1=1,6
   Slept(i1)%g = Slept2(i1)%g
   Slept(i1)%gi2 = Slept2(i1)%gi2
   Slept(i1)%bi2 = Slept2(i1)%bi2
   Slept(i1)%id2 = Slept2(i1)%id2
  End Do

  !------------------------------------------------
  ! 3-body decays of sleptons
  !------------------------------------------------
  Call Slepton_3body(-1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u, n_c     &
    & , n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_S0, n_P0, n_Spm             &
    & , Slept, mZ, mW, mf_l, mf_u, mf_d, Chi0, ChiPm, P0, S0, Spm, Sneut       &
    & , c_P0SlSl, c_S0SlSl, c_SlSlZ, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn    &
    & , c_CNuSl_L, c_CNuSl_R, c_NuNuZ_L, c_NuNuZ_R, c_LLZ_L, c_LLZ_R, c_UUZ_L  &
    & , c_UUZ_R, c_DDZ_L, c_DDZ_R, c_LNuW, c_DUW, c_LLP0_L, c_LLP0_R           &
    & , c_UUP0_L, c_UUP0_R, c_DDP0_L, c_DDP0_R, c_LLS0_L, c_LLS0_R, c_UUS0_L   &
    & , c_UUS0_R, c_DDS0_L, c_DDS0_R, c_CLSn_L, c_CLSn_R, c_NuNSn_L, c_NuNSn_R &
    & , c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R, GenerationMixing, epsI)

  !------------------------------------------------
  ! 2-body decays of u-squarks
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Su, n_u, id_u, n_n, n_g, n_c, n_d, id_d &
          & , n_W, id_W, n_Sd, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sup2 &
          & , mf_u, mf_d, Chi0, c_UNSu_L, c_UNSu_R, ChiPm, c_CDSu_L  &
          & , c_CDSu_R, Sdown, m_W, c_SuSdW, m_Z, c_SuSuZ, Spm       &
          & , c_SmpSuSd, P0, c_P0SuSu, S0, c_S0SuSu, m32, F_eff      &
          & , c_GraUSu_L, c_GraUSu_R, 0, Glu, c_GUSu_L, c_GUSu_R)
  Do i1=1,6
   Sup(i1)%g = Sup2(i1)%g
   Sup(i1)%gi2 = Sup2(i1)%gi2
   Sup(i1)%bi2 = Sup2(i1)%bi2
   Sup(i1)%id2 = Sup2(i1)%id2
  End Do

  !------------------------------------------------
  ! 2-body decays of d-squarks
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Sd, n_d, id_d, n_n, n_g, n_c, n_u, id_u    &
          & , n_W, id_W+1, n_Su, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sdown &
          & , mf_d, mf_u, Chi0, c_DNSd_L, c_DNSd_R, ChiM, c_CUSd_L      &
          & , c_CUSd_R, Sup2, m_W, c_SdSuW, m_Z, c_SdSdZ, Smp           &
          & , c_SmpSdSu, P0, c_P0SdSd, S0, c_S0SdSd, m32, F_eff         &
          & , c_GraDSd_L, c_GraDSd_R, 0, Glu, c_GDSd_L, c_GDSd_R)

  !------------------------------------------------
  ! gluino decays
  !------------------------------------------------
  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(n_d, id_d, n_Sd, n_u, id_u, n_Su, id_grav, id_gl &
     & , Glu, Sdown, c_GDSd_L, c_GDSd_R, mf_d, Sup, c_GUSu_L, c_GUSu_R, mf_u &
     & , m32, Fgmsb, 0 )

   If (Glu%g.Lt.fac3*Glu%m) Then
    Glu%gi2 = 0._dp
    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .False. )

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .True. )

   End If

  Else ! calculation of 3-body decay modes is enforced

   Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
      & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
      & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
      & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
      & , c_GDSd_R, m_W, epsI, deltaM, .True. )

  End If

  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  !------------------------------------------------------------------------
  ! couplings to a real and a virtual gauge boson
  !------------------------------------------------------------------------
  c_S0WWvirt = c_S0WW / vev
  c_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * c_S0ZZ / vev
  !------------------------------------------------------------------------
  ! loop induced couplings to photons and gluons
  !------------------------------------------------------------------------
  Do i1=2,2
   g_W = (vevSM(1) * RS0(i1,1) + vevSM(2) * RS0(i1,2) ) / vev
   g_u = RS0(i1,2) * vev / vevSM(2)
   g_d = RS0(i1,1) * vev / vevSM(1)
   g_l = RS0(i1,1) * vev / vevSM(1)
   Forall(i2=1:2) g_c(i2) = c_CCS0_L(i2,i2,i1) * vev / ChiPm(i2)%m
   g_Hp = c_SmpS03(2,2,i1) * vev /( 2._dp * Spm(2)%m2)
   Forall(i2=1:6) g_sl(i2) = - 0.5_dp * c_S0SlSl(i1,i2,i2) * vev &
                           &          / Slept(i2)%m2
   Forall(i2=1:6) g_sd(i2) = - 0.5_dp * c_S0SdSd(i1,i2,i2) * vev &
                           &          / Sdown(i2)%m2
   Forall(i2=1:6) g_su(i2) = - 0.5_dp * c_S0SuSu(i1,i2,i2) * vev &
                           &          / Sup(i2)%m2
   Call CoupScalarPhoton(S0(i1)%m2, mW2, g_W, mf_u2, g_u, 1._dp, mf_d2, g_d &
    & , mf_l2, g_l, ChiPm%m2, g_c, Spm(2)%m2, g_Hp, Sup%m2, g_su, Sdown%m2  &
    & , g_sd, 1._dp, Slept%m2, g_sl, c_GGS0(i1), coup )

   c_GGS0(i1) = c_GGS0(i1) * oo4pi * gauge(2)**2 * sinW2

   Call CoupScalarGluon(S0(i1)%m2, mf_u2, g_u, mf_d2, g_d, Sup%m2, g_su &
                      & , Sdown%m2, g_sd, c_GlGlS0(i1), coup )
   r_GlGlS0(i1) = Abs(c_GlGlS0(i1) / coup)**2
   c_GlGlS0(i1) = c_GlGlS0(i1) * oo4pi * gauge(3)**2
  End Do

!-----------------------------------------------
! include running of couplings to m_h for h^0
!-----------------------------------------------
  Call GetParameters_at_Q(S0(1)%m, g_i, Y_l_h, Y_d_h, Y_u_h,  Mi, A_l_h &
              & , A_d_h, A_u_h, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu_h, B)

  If (GenerationMixing) Then
   Call FermionMass(Y_l_h,vevSM(1),mF, uL_L_h,uL_R_h,kont)
   mf_l2_h = mf**2
   Call FermionMass(Y_d_h,vevSM(1),mF, uD_L_h,uD_R_h,kont)
   mf_d2_h = mf**2
   Call FermionMass(Y_u_h,vevSM(2),mF, uU_L_h,uU_R_h,kont)
   mf_u2_h = mf**2
  Else
   Forall(i1=1:3) mf_l2_h(i1) = (vevSM(1)*Abs(Y_l_h(i1,i1)))**2
   Forall(i1=1:3) mf_d2_h(i1) = (vevSM(1)*Abs(Y_d_h(i1,i1)))**2
   Forall(i1=1:3) mf_u2_h(i1) = (vevSM(2)*Abs(Y_u_h(i1,i1)))**2
  End If

  Call TreeMassesMSSM(g_i(1), g_i(2), vevSM, Mi(1), Mi(2), Mi(3), mu_h, B     &
              & , tanb_mZ, M2_E, M2_L, A_l_h, Y_l_h                           &
              & , M2_D, M2_U, M2_Q, A_d_h, A_u_h, Y_d_h, Y_u_h                &
              & , mGlu, PhiGlu, mC, mC2, U_h, V_h, mN, mN2, N_h               &
              & , mSneut, mSneut2, Rsneut_h, mSlepton, mSlepton2, RSlepton_h  &
              & , mSdown, mSdown2, RSdown_h, mSup, mSup2, RSup_h              &
              & , mP0, mP02, RP0_h, mS0, mS02, RS0_h, mSpm, mSpm2, RSpm_h     &
              & , GenerationMixing, kont, .False. )

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
      Call CoupFermionScalar(i1, i2, 1, -0.5_dp, Y_d_h, uD_L_h, uD_R_h, RS0 &
                           &, c_DDS0_L(i1,i2,1), c_DDS0_R(i1,i2,1))
      Call CoupFermionScalar(i1, i2, 1, -0.5_dp, Y_l_h, uL_L_h, uL_R_h, RS0 &
                           &, c_LLS0_L(i1,i2,1), c_LLS0_R(i1,i2,1))
      Call CoupFermionScalar(i1, i2, 1, 0.5_dp, Y_u_h, uU_L_h, uU_R_h, RS0 &
                           &, c_UUS0_L(i1,i2,1), c_UUS0_R(i1,i2,1))
    End Do
   End Do
  Else
   Do i1 = 1,3
     Call CoupFermionScalar(1, -0.5_dp, Y_d_h(i1,i1), RS0   &
                          &, c_DDS0_L(i1,i1,1), c_DDS0_R(i1,i1,1))
     Call CoupFermionScalar(1, -0.5_dp, Y_l_h(i1,i1), RS0   &
                          &, c_LLS0_L(i1,i1,1), c_LLS0_R(i1,i1,1))
     Call CoupFermionScalar(1, 0.5_dp, Y_u_h(i1,i1), RS0   &
                          &, c_UUS0_L(i1,i1,1), c_UUS0_R(i1,i1,1))
   End Do
  End If

  Call CoupScalarW(1, g_i(2), vevSM, RS0, c_S0WW(1,1) )
  sinW2 = g_i(1)**2 / (g_i(1)**2 + g_i(2)**2)
  cosW2 = 1._dp-sinW2
  Call CoupScalarZ(1, g_i(2), cosW2, vevSM, RS0, c_S0ZZ(1,1) )
  c_S0WWvirt(1,1) = c_S0WW(1,1) / vev
  c_S0ZZvirt(1,1) = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * c_S0ZZ(1,1) / vev
  Call CoupChargedScalarScalar3(2, 2, 1, RSpm_h, RS0, vevSM, g_i(1), g_i(2) &
                                   &, c_SmpS03(2,2,1) )
  g_Hp = c_SmpS03(2,2,1) * vev /( 2._dp * Spm(2)%m2)

  g_u = RS0(1,2) * vev / vevSM(2)
  g_d = RS0(1,1) * vev / vevSM(1)
  g_l = RS0(1,1) * vev / vevSM(1)
  Do i2=1,2
   Call CoupCharginoScalar(i2, i2, 1, U_h, V_h, RS0, g_i(2), g_c(i2), cR)
   g_c(i2) = g_c(i2) * vev / mC(i2)
  End Do

  If (GenerationMixing) Then
   Do i2=1,6
    Call CoupScalarSfermion3(1, i2, i2, RS0, -0.5_dp, e_d, Y_d_h, Rsdown_h   &
                           &, A_d_h, mu_h, vevSM, g_i(1), g_i(2), g_sd(i2) )
    g_sd(i2) =  - 0.5_dp * g_sd(i2) * vev / mSdown2(i2)
    Call CoupScalarSfermion3(1, i2, i2, RS0, 0.5_dp, e_u, Y_u_h, Rsup_h      &
                           &, A_u_h, mu_h, vevSM, g_i(1), g_i(2), g_su(i2) )
    g_su(i2) =  - 0.5_dp * g_su(i2) * vev / mSup2(i2)
    Call CoupScalarSfermion3(1, i2, i2, RS0, -0.5_dp, -1._dp, Y_l_h        &
                 &, Rslepton_h, A_l_h, mu_h, vevSM, g_i(1), g_i(2), g_sl(i2) )
    g_sl(i2) =  - 0.5_dp * g_sl(i2) * vev / mSlepton2(i2)
   End Do
  Else
   Do i2=1,3
    Rsl = RSlepton_h(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
    Rsd = RSdown_h(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
    Rsu = RSup_h(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

    yuk = Y_l_h(i2,i2)
    A = A_l_h(i2,i2)
    Do i1=1,2
     Call CoupScalarSfermion3(1, i1, i1, RS0, -0.5_dp, -1._dp, Yuk, Rsl  &
                               &, A, mu_h, vevSM, g_i(1), g_i(2), cR )
     g_sl( (i2-1)*2+i1 ) =  - 0.5_dp * cR * vev / mSlepton2( (i2-1)*2 + i1 )
    End Do

    yuk = Y_d_h(i2,i2)
    A = A_d_h(i2,i2)
    Do i1=1,2
     Call CoupScalarSfermion3(1, i1, i1, RS0, -0.5_dp, e_d, Yuk, Rsd  &
                               &, A, mu_h, vevSM, g_i(1), g_i(2), cR )
     g_sd( (i2-1)*2+i1 ) =  - 0.5_dp * cR * vev / mSdown2( (i2-1)*2 + i1 )
    End Do

    yuk = Y_u_h(i2,i2)
    A = A_u_h(i2,i2)
    Do i1=1,2
     Call CoupScalarSfermion3(1, i1, i1, RS0, 0.5_dp, e_u, Yuk, Rsu  &
                               &, A, mu_h, vevSM, g_i(1), g_i(2), cR )
     g_su( (i2-1)*2+i1 ) =  - 0.5_dp * cR * vev / mSup2( (i2-1)*2 + i1 )
    End Do

   End Do
  End If

  g_W = (vevSM(1) * RS0(1,1) + vevSM(2) * RS0(1,2) ) / vev
  !-----------------------------------------------------------------
  ! including approx. QCD corrections in case of lighter h^0
  ! based on Graudenz, Djouadi, Spira, Zerwas NPB (1995), for gamma 
  !-----------------------------------------------------------------
  r_T = 1._dp - oo4pi2 * g_i(3)**2
  fakt = 1._dp + 2._dp * oo3pi2 * g_i(3)**2
  Call CoupScalarPhoton(S0(1)%m2, mW2, g_W, mf_u2_h, g_u, r_T, mf_d2_h, g_d   &
    & , mf_l2_h, g_l, mC2, g_c, Spm(2)%m2, g_Hp, mSup2, g_su, mSdown2 &
    & , g_sd, fakt, mSlepton2, g_sl, c_GGS0(1), coup )

  ! effectively we have the SM outside
  c_GGS0(1) = c_GGS0(1) * alpha_125 ! oo4pi * g_i(2)**2 * sinW2
  !---------------------------------------------------------------------
  ! including approx. QCD corrections in case of lighter h^0
  ! based on Graudenz, Djouadi, Spira, Zerwas NPB453 (1995) 17, gluons
  ! taking out the gamma gamma part, adding gluon gluon part
  ! sqrt appears as I change here the coupling and not the width
  !---------------------------------------------------------------------
  g_u = g_u * fakt 
  g_d = g_d * fakt 
  g_su = g_su * fakt 
  g_sd = g_sd * fakt 
  Call CoupScalarGluon(S0(1)%m2, mf_u2, g_u, mf_d2, g_d, mSup2, g_su &
                      & , mSdown2, g_sd, c_GlGlS0(1), coup )
  r_GlGlS0(1) = Abs(c_GlGlS0(1) / coup)**2
  fakt = Sqrt(1._dp + 215._dp * oo48pi2 * g_i(3)**2)
  c_GlGlS0(1) = c_GlGlS0(1) * alphaS_125 * fakt ! * oo4pi * g_i(3)**2 * fakt

  Call ScalarTwoBodyDecays(-1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u    &
      & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_p0 &
      & , n_Spm, id_ph, id_gl, S0, c_S03, c_GlGlS0, c_GGS0, mf_l, c_LLS0_L     &
      & , c_LLS0_R, mf_d, c_DDS0_L, c_DDS0_R, mf_u, c_UUS0_L, c_UUS0_R, Slept  &
      & , c_S0SlSl, Sneut, c_S0SnSn, Sdown, c_S0SdSd, Sup, c_S0SuSu, Chi0      & 
      & , c_NNS0_L, c_NNS0_R, ChiPm, c_CCS0_L, c_CCS0_R, m_W, c_S0WW           &
      & , c_S0WWvirt, m_Z, c_S0ZZ, c_S0ZZvirt, Spm, c_SmpS03, P0, c_P0S03      &
      & , c_P0S0Z, c_SmpS0W, Glu%m)
  
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  !------------------------------------------------------------------------
  ! loop induced couplings to photons and gluons
  !------------------------------------------------------------------------
  c_GlGlP0 = 0._dp
  c_GGP0 = 0._dp
  Do i1=2,2
   g_u = RP0(i1,2) * vev / vevSM(2)
   g_d = RP0(i1,1) * vev / vevSM(1)
   g_l = RP0(i1,1) * vev / vevSM(1)
   Forall(i2=1:2) g_c(i2) = c_CCS0_L(i2,i2,i1) * mW /( gauge(2) * ChiPm(i2)%m)
   g_Hp = c_SmpS03(2,2,i1) * mW /( gauge(2) * Spm(2)%m2)
   Call CoupPseudoScalarPhoton(P0(i1)%m2, mf_u2, g_u, mf_d2, g_d, mf_l2, g_l &
          & , ChiPm%m2, g_c, c_GGP0(i1), coup )
   c_GGP0(i1) = c_GGP0(i1) * oo4pi * gauge(2)**2 * sinW2
   Call CoupPseudoScalarGluon(P0(i1)%m2, mf_u2, g_u, mf_d2, g_d, c_GlGlP0(i1) &
          & , coup )
   r_GlGlP0 = Abs(c_GlGlP0(i1) / coup)
   c_GlGlP0(i1) = c_GlGlP0(i1) * oo4pi * gauge(3)**2
  End Do

  Call PseudoscalarTwoBodyDecays(2, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d    &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_sle, n_Sd, n_su, n_n, n_c, n_p0   &
      & , n_Spm, id_ph, id_gl, P0, mf_l, c_LLP0_L, c_LLP0_R, mf_d, c_DDP0_L    &
      & , c_DDP0_R, mf_u, c_UUP0_L, c_UUP0_R, Slept, c_P0SlSl, Sdown, c_P0SdSd &
      & , Sup, c_P0SuSu, Chi0, c_NNP0_L, c_NNP0_R, ChiPm, c_CCP0_L, c_CCP0_R   &
      & , Spm, c_SmpP03, S0, c_P0S03, m_Z, c_P0S0Z, m_W, c_SmpP0W, c_GlGlP0    &
      & , c_GGP0, Glu%m )

  Call ChargedscalarTwoBodyDecays(2, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d  &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c &
      & , n_p0, n_Spm, Spm, mf_l, c_SmpLNu_L, c_SmpLNu_R, mf_d, mf_u          &
      & , c_SmpDU_L, c_SmpDU_R, Slept, Sneut, c_SmpSlSn, Sdown, Sup           &
      & , c_SmpSdSu, Chi0, ChiPm, c_SmpCN_L, c_SmpCN_R, m_W, m_Z, c_SmpZ, P0  &
      & , c_SmpP03, c_SmpP0W, S0, c_SmpS03, c_SmpS0W, 1)

  Do i1=1,2
   c_CGW_L(1) = U(i1,1) /F_eff 
   c_CGW_R(1) = V(i1,1) /F_eff 
   c_CGW_L(2) = oosqrt2 * U(i1,2) * cosb /F_eff 
   c_CGW_R(2) = oosqrt2 *  V(i1,2) * sinb /F_eff 
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u      &
       & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c      &
       & , n_s0, n_p0, n_Spm, id_grav, ChiPm, Slept, c_CNuSl_L, c_CNuSl_R      &
       & , Sneut, c_CLSn_L, c_CLSn_R, mf_l, Sdown, c_CUSd_L, c_CUSd_R, mf_u    &
       & , Sup, c_CDSu_L, c_CDSu_R, mf_d, Chi0, m_W, c_CNW_L, c_CNW_R, Spm     &
       & , c_SmpCN_L, c_SmpCN_R, m_Z, c_CCZ_L, c_CCZ_R, P0, c_CCP0_L, c_CCP0_R &
       & , S0, c_CCS0_L, c_CCS0_R, m32, c_CGW_L, c_CGW_R, 1)
    k_neut =0
    If (ChiPm(i1)%g.Lt. fac3*Abs(ChiPm(i1)%m)) Then
     ChiPm(i1)%gi2 = 0._dp
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R           &
       & , c_CLSn_L, c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

    Else ! calculate 3-body decays with virtual interemdiate states only

     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R           &
       & , c_CLSn_L, c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R &
       & , GenerationMixing, k_neut, epsI, deltaM, .True. )
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R           &
       & , c_CLSn_L, c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

   End If

   If (i1.Eq.1) Then ! then calculate decay into pion and check this
                     ! partial width with respect to the total width
    Call  FermionToFermionPi(ChiPm(1)%m,Chi0(1)%m,m_pip,f_pi, G_F &
           & , c_CNW_L(1,1,1),c_CNW_R(1,1,1), width)

    If (width.Gt.ChiPm(1)%g) Then ! PDG 111 pi0, 211 pi+
     ChiPm(1)%gi2(1) = width
     ChiPm(1)%id2(1,1) = Chi0(1)%id
     ChiPm(1)%id2(1,2) = id_piP
     ChiPm(1)%gi3(1:9) = 0._dp
     ChiPm(1)%g = Sum(ChiPm(1)%gi2) + Sum(ChiPm(1)%gi3)
     ChiPm(1)%bi2(1) = width / ChiPm(1)%g
     ChiPm(1)%bi3 = ChiPm(1)%gi3 / ChiPm(1)%g
    End If 
   End If
  End Do
  
  OnlySM = .True.

  Do i1=1,4
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    c_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / F_eff
    c_NGZ(1) = (- N(i1,1)*sinW+ N(i1,2)*cosW) / F_eff
    c_NGZ(2) = oosqrt2 * (cosb * N(i1,3) - sinb * N(i1,4) ) / F_eff
    c_NGH =  oosqrt2 * (N(i1,3) * RS0(1,1) + N(i1,4)*RS0(1,2)) / F_eff
    Call NeutralinoTwoBodyDecays(i1, Chi0, n_nu, id_nu, n_l, id_l, n_d, id_d   &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W , n_snu, n_sle, n_Sd, n_su, n_n, n_c &
      & , n_s0, n_p0, n_Spm, id_ph, id_grav, Slept, c_LNSl_L, c_LNSl_R, mf_l   &
      & , Sneut, c_NuNSn_L, c_NuNSn_R, Sdown, c_DNSd_L, c_DNSd_R, mf_d, Sup    &
      & , c_UNSu_L, c_UNSu_R, mf_u, ChiPm, m_W, c_CNW_L, c_CNW_R, Spm          &
      & , c_SmpCN_L, c_SmpCN_R, m_Z, c_NNZ_L, c_NNZ_R, P0, c_NNP0_L, c_NNP0_R  &
      & , S0, c_NNS0_L, c_NNS0_R, m32, c_NGP, c_NGZ, c_NGH, 1 )

    If (Chi0(i1)%g.Lt.fac3*Abs(Chi0(i1)%m)) Then
     Chi0(i1)%gi2 = 0._dp
     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 2500)

    Else ! calculate only 3-body via virtual particles
     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .True., 0._dp, 200    &
      & , 2500)

    End If

   Else ! calculation of 3-body decay modes is enforced

     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 2500)

    End If ! CTBD

   If (i1.Eq.2) Then ! then calculate decay into pion and check this
                     ! partial width with respect to the total width
    Call  FermionToFermionPi(Chi0(i1)%m,Chi0(1)%m,m_pi0,f_pi, G_F/cosW2 &
           & , c_NNZ_L(i1,1,1),c_NNZ_R(i1,1,1), width)

    If (width.Gt.Chi0(2)%g) Then ! PDG 111 pi0, 211 pi+
     Chi0(2)%gi2(1) = width
     Chi0(2)%id2(1,1) = Chi0(1)%id
     Chi0(2)%id2(1,2) = id_pi0
     Chi0(2)%gi3(1:18) = 0._dp
     Chi0(2)%g = Sum(Chi0(2)%gi2) + Sum(Chi0(2)%gi3)
     Chi0(2)%bi2(1) = width / Chi0(2)%g
     Chi0(2)%bi3 = Chi0(2)%gi3 / Chi0(2)%g
    End If 
   End If

  End Do

  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
  If (GenerationMixing) Then
   Write(ErrCan,*) "Warning 3-body decays of ~t_1 are not in case of"
   Write(ErrCan,*) "generation mixing"
   Do i1=1,6
    Sup(i1)%gi3 = 0
    Sup(i1)%bi3 = 0
   End Do
  Else 
   If ((Sup(5)%g.Lt.fac3*Sup(5)%m).Or.CTBD) Then ! calculation including widths
    Sup(5)%gi2 = 0._dp
    !--------------------------------------------
    ! check if flavour violating decays exist
    !--------------------------------------------
    Do i1=1,2
     Call ScalarToTwoFermions(Sup(5)%m, mf_u(2), Chi0(i1)%m &
            & , c_UNSu_L(2,i1,5), c_UNSu_R(2,i1,5), gam)
     Sup(5)%gi2(i1) = gam
     Sup(5)%id2(i1,1) = Chi0(i1)%id
     Sup(5)%id2(i1,2) = id_u(2)
    End Do
    Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d   &
     & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
     & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
     & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
     & , Slept, c_CNuSl_R, epsI, .False.)

   Else  ! calculation excluding widths

    Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d   &
     & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
     & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
     & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
     & , Slept, c_CNuSl_R, epsI, .True.)

   End If
  End If ! generation mixing

  Iname = Iname - 1

 Contains

  Subroutine LoopCouplingsMSSM(GenerationMixing, yukB, tanb, A_b, mu, mSbottom2 &
                  & , mStop2, RStop, mP02 &
                  & , CKM, c_UNSu_R, c_GUSu_R)
  Implicit None
   Logical, Intent(in) :: GenerationMixing
   Real(dp), Intent(in) :: tanb, mSbottom2(6), mStop2(6), mP02(2)
   Complex(dp), Intent(in) :: YukB, RStop(6,6), A_b, mu, CKM(3,3)
   Complex(dp), Intent(inout) :: c_UNSu_R(3,4,6), c_GUSu_R(3,6)

   Integer ::  kont
   
   Real(dp) :: cosb2, sinb2, cos2b, fakt16pi, deltal, deltar, mass2(3), test(2)
   Complex(dp) :: Rsf(3,3), mat3(3,3)

   If (.Not.GenerationMixing) Then ! some additional couplings, loop induced
    cosb2 = 1._dp / (1._dp + tanb**2)
    sinb2 = tanb**2 * cosb2
    cos2b = cosb2 - sinb2
    ! the log is approximately log(m^2_gut/ 1 TeV^2)
    fakt16pi = Abs(YukB)**2 *oo16pi2 * CKM(2,3) * CKM(3,3) * Log(1.e26_dp) &
           & / cosb2

    ! corrected version with electroweak symmetry breaking
    deltal = - fakt16pi * (mSbottom2(5)+mSbottom2(6) + sinb2 * mP02(2)  &
          &           - Abs(mu)**2- 0.5_dp * cos2b * mP02(1) + Abs(A_b/yukB)**2 )
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

    Call EigenSystem(mat3, mass2, Rsf, kont, test)
    If ((kont.Eq.-14).Or.(kont.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = 0
    End If
    If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then
!   call AddNOW() 
     Write(ErrCan,*) "Warning, in subroutine "//NameOfUnit(Iname)
     Write(ErrCan,*) "Diagonalization of mat3 has failed",kont
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If

    test(1) = Maxval( Abs(Rsf(3,:) ) )
    If (test(1).Eq.Abs(Rsf(3,1))) Then
     c_UNSu_R(2,1:2,5) = Rsf(3,2) * c_UNSu_R(2,1:2,4)
     c_GUSu_R(2,5) = Rsf(3,2) * c_GUSu_R(2,4)
    Else
     c_UNSu_R(2,1:2,5) = Rsf(3,1) * c_UNSu_R(2,1:2,4)
     c_GUSu_R(2,5) = Rsf(3,1) * c_GUSu_R(2,4)
    End If

   End If ! GenerationMixing


  End Subroutine LoopCouplingsMSSM

 End Subroutine CalculateBR_MSSM

 Subroutine CalculateBR_NMSSM(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u    &
    & , n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0   &
    & , n_p0, n_Spm, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, ChiPm, U, V &
    & , Chi0, N, Sneut, RSneut, Slept, RSlepton, Sup, RSup, Sdown, RSdown     &
    & , uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, S0, RS0, P0, RP0, Spm, RSpm       &
    & , epsI, deltaM, CTBD, fac3, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, h0, Ah0   &
    & , lam, Alam, vevSM, vP, Fgmsb, m32, grav_fac) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - kont ........ control variable, is 0 if everything is o.k.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 ! 28.02.14: changing units of m32 to GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu, n_sle, n_Sd &
     & , n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm, id_grav, id_gl, id_ph 
  Integer, Intent(in), Dimension(1) :: id_Z, id_W
  Integer, Intent(in), Dimension(3) :: id_nu, id_l, id_d, id_u

  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: RP0(3,3), RS0(3,3)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(2,2), U(2,2), V(2,2), N(5,5) &
         & , RSlepton(6,6), Rsneut(3,3), RSup(6,6), RSdown(6,6), uL_L(3,3) &
         & , uL_R(3,3), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu, h0, Ah0, lam, Alam
  Real(dp), Intent(in) :: vevSM(2), vP, Fgmsb, m32, grav_fac
  Logical, Intent(in) :: CTBD

  Type(particle2), Intent(inout) :: Sdown(6), Spm(2), P0(3)
  Type(particle23), Intent(inout) :: Slept(6), Chi0(5), Sup(6), ChiPm(2), Glu &
      & , S0(3), Sneut(3)

   ! contains only the information on the 2-body
  Type(particle2) :: Sup2(6), Slept2(6), Sneut2(3)
  Type(particle2) :: SMp(2) ! contains the identies of the negative charged ones
  Type(particle23) :: ChiM(2) ! contains the identies of the negative charged ones

  Integer :: i1, k_neut
  Real(dp) :: sinW

  Complex(dp) :: c_SmpSlSn(2,6,3), c_SmpSdSu(2,6,6)   &
      & , c_SmpSnSl(2,3,6), c_SmpSuSd(2,6,6), c_SmpP03(2,2,3)    &
      & , c_SmpP0W(2,3,1), c_SmpS03(2,2,3), c_SmpS0W(2,3,1)          &
      & , c_SmpLNu_L(2,3,3), c_SmpLNu_R(2,3,3), c_SmpDU_L(2,3,3) &
      & , c_SmpDU_R(8,3,3), c_SmpZ(2,2,1), c_DUW(3,3)
  Real(dp) :: c_LLZ_L, c_LLZ_R, c_DDZ_L, c_DDZ_R  &
      & , c_UUZ_L, c_UUZ_R, c_NuNuZ_L, c_NuNuZ_R
  Complex(dp) :: c_CCZ_L(2,2,1), c_CCZ_R(2,2,1)           &
      & , c_NNZ_L(5,5,1), c_NNZ_R(5,5,1), c_NNS0_L(5,5,3)            &
      & , c_NNS0_R(5,5,3), c_NNP0_L(5,5,3), c_NNP0_R(5,5,3) 
  Complex(dp) :: c_GDSd_L(3,6), c_GDSd_R(3,6)          &
      & , c_DNSd_L(3,5,6), c_DNSd_R(3,5,6), c_GUSu_L(3,6)         &
      & , c_GUSu_R(3,6), c_UNSu_L(3,5,6), c_UNSu_R(3,5,6)         &
      & , c_LNSl_L(3,5,6), c_LNSl_R(3,5,6), c_NuNSn_L(3,5,3)      & 
      & , c_NuNSn_R(3,5,3), c_DDP0_L(3,3,3), c_LLP0_L(3,3,3)      &
      & , c_UUP0_L(3,3,3), c_DDP0_R(3,3,3), c_LLP0_R(3,3,3)       &
      & , c_UUP0_R(3,3,3), c_DDS0_L(3,3,3), c_LLS0_L(3,3,3)       &
      & , c_UUS0_L(3,3,3), c_DDS0_R(3,3,3), c_LLS0_R(3,3,3)       &
      & , c_UUS0_R(3,3,3)
  Complex(dp) :: c_CUSd_L(2,3,6), c_CUSd_R(2,3,6)      &
      & , c_CDSu_L(2,3,6), c_CDSu_R(2,3,6), c_CLSn_L(2,3,3)       &
      & , c_CLSn_R(2,3,3), c_CNuSl_L(2,3,6), c_CNuSl_R(2,3,6)
  Complex(dp) :: c_GlGlS0(3), c_GGS0(3), c_GlGlP0(3), c_GGP0(3)
  Complex(dp) :: c_P0SdSd(3,6,6), c_P0SuSu(3,6,6), c_P0SlSl(3,6,6) &
      & , c_P0SnSn(3,3,3), c_P0S0Z(3,3,1) 
  Real(dp) :: c_P0S03(3,3,3), c_S0WWvirt(3,1), c_S0ZZvirt(3,1), vev
  Complex(dp) :: c_S0SdSd(3,6,6), c_S0SuSu(3,6,6) &
      & , c_S0SlSl(3,6,6), c_S0SnSn(3,3,3), c_LNuW(3,3)
  Real(dp) :: c_S03(3,3,3), c_S0WW(3,1), c_S0ZZ(3,1), c_FFpW, sinb, cosb
  Complex(dp) :: c_SdSuW(6,6,1), c_SuSdW(6,6,1), c_SlSnW(6,3,1) &
      & , c_SnSlW(3,6,1), c_SdSdZ(6,6,1), c_SlSlZ(6,6,1), c_SnSnZ(3,3,1)     &
      & , c_SuSuZ(6,6,1)
  Complex(dp), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R, c_GraLSl_L, c_GraLSl_R
  Complex(dp), Dimension(3,3) :: c_GraNuSn_L, c_GraNuSn_R
  Complex(dp) :: c_CCP0_L(2,2,3), c_CCP0_R(2,2,3)                     &
      & , c_CCS0_L(2,2,3), c_CCS0_R(2,2,3), c_CNW_L(2,5,1)            &
      & , c_CNW_R(2,5,1), c_SmpCN_L(2,2,5), c_SmpCN_R(2,2,5), c_NGP &
      & , c_NGZ(2), c_NGH, c_CGW_L(2), c_CGW_R(2)

  Real(dp) :: tanb, sinW2, cosW, m_W(1), m_Z(1), F_eff, gam
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_NMSSM'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  
  m_W = mW
  m_Z = mZ
  Spm(1)%g = gamW
  P0(1)%g = gamZ

  Sup2%m = Sup%m 
  Sup2%id = Sup%id 
  Slept2%m = Slept%m 
  Slept2%id = Slept%id 
  Sneut2%m  = Sneut%m 
  Sneut2%id = Sneut%id 
  ChiM%m = ChiPm%m
  ChiM%id = ChiPm%id + 1
  SMp%m = SPm%m
  SMp%id = SPm%id + 1

  Call AllCouplings(gauge, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R &
    & , vevSM, h0, Ah0, lam, Alam, vP, RSpm, RP0, RS0, U, V, N, mu, PhaseGlu &
    & , RSlepton, A_l, Rsneut, RSup, A_u, RSdown, A_d                        &
    & , c_SmpSlSn, c_SmpSdSu, c_SmpSnSl, c_SmpSuSd, c_SmpP03       &
    & , c_SmpP0W(:,:,1), c_SmpS03, c_SmpS0W(:,:,1), c_SmpLNu_L, c_SmpLNu_R       &
    & , c_SmpDU_L, c_SmpDU_R, c_SmpZ(:,:,1), c_DUW, c_LLZ_L, c_LLZ_R    &
    & , c_DDZ_L, c_DDZ_R, c_UUZ_L, c_UUZ_R, c_NuNuZ_L, c_NuNuZ_R &
    & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_NNS0_L, c_NNS0_R   &
    & , c_NNP0_L, c_NNP0_R, c_GDSd_L, c_GDSd_R, c_DNSd_L           &
    & , c_DNSd_R, c_GUSu_L, c_GUSu_R, c_UNSu_L, c_UNSu_R           &
    & , c_LNSl_L, c_LNSl_R, c_NuNSn_L, c_NuNSn_R, c_DDP0_L         &
    & , c_LLP0_L, c_UUP0_L, c_DDP0_R, c_LLP0_R, c_UUP0_R           &
    & , c_DDS0_L, c_LLS0_L, c_UUS0_L, c_DDS0_R, c_LLS0_R           &
    & , c_UUS0_R, c_CUSd_L, c_CUSd_R, c_CDSu_L, c_CDSu_R           &
    & , c_CLSn_L, c_CLSn_R, c_CNuSl_L, c_CNuSl_R, c_GlGlS0         &
    & , c_P0SdSd, c_P0SuSu, c_P0SlSl, c_P0SnSn, c_P0S0Z(:,:,1), c_P0S03 &
    & , c_S0SdSd, c_S0SuSu, c_S0SlSl, c_S0SnSn, c_S03, c_S0WW(:,1)    &
    & , c_S0ZZ(:,1), c_FFpW, c_LNuW, c_SdSuW(:,:,1), c_SuSdW(:,:,1), c_SlSnW(:,:,1)        &
    & , c_SnSlW(:,:,1), c_SdSdZ(:,:,1), c_SlSlZ(:,:,1), c_SnSnZ(:,:,1), c_SuSuZ(:,:,1), c_CCP0_L    &
    & , c_CCP0_R, c_CCS0_L, c_CCS0_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)             &
    & , c_SmpCN_L, c_SmpCN_R, c_GraDSd_L, c_GraDSd_R, c_GraUSu_L, c_GraUSu_R &
    & , c_GraLSl_L, c_GraLSl_R, c_GraNuSn_L, c_GraNuSn_R, GenerationMixing)


  tanb = vevSM(2) / vevSM(1)
  cosb = 1._dp / Sqrt(1._dp + tanb**2)
  sinb = cosb * tanb 
   
  If (grav_fac.Ne.0._dp) Then
   F_eff = grav_fac*Fgmsb ! effective SUSY breaking scale in case of GSMB
  Else 
   F_eff = Fgmsb
  End If

  !------------------------------------------------
  ! 2-body decays of sneutrinos
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Snu, n_nu, id_nu, n_n, 0, n_c, n_l, id_l   &
          & , n_W, id_W, n_Sle, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sneut2 &
          & , mf_nu, mf_l, Chi0, c_NuNSn_L, c_NuNSn_R, ChiPm, c_CLSn_L        &
          & , c_CLSn_R, Slept2, m_W, c_SnSlW, m_Z, c_SnSnZ, Spm, c_SmpSnSl    &
          & , P0, c_P0SnSn, S0, c_S0SnSn, m32, F_eff, c_GraNuSn_L          &
          & , c_GraNuSn_R, 1)

  Do i1=1,3
   Sneut(i1)%g = Sneut2(i1)%g
   Sneut(i1)%gi2 = Sneut2(i1)%gi2
   Sneut(i1)%bi2 = Sneut2(i1)%bi2
   Sneut(i1)%id2 = Sneut2(i1)%id2
  End Do 

  !------------------------------------------------
  ! 3-body decays of sneutrinos
  !------------------------------------------------
  Call Sneutrino_3body(-1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u, n_c   &
    & , n_n, n_W, id_W, n_Sle, n_Snu, n_Spm, Sneut, mW, mf_l, mf_u, mf_d       &
    & , Chi0, ChiPm, Spm, Slept, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn        &
    & , c_CNuSl_L, c_CNuSl_R, c_LNuW, c_DUW, c_CLSn_L, c_CLSn_R, c_NuNSn_L     &
    & , c_NuNSn_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R                &
    & , GenerationMixing, epsI)

  !------------------------------------------------
  ! 2-body decays of sleptons
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Sle, n_l, id_l, n_n, 0, n_c, n_nu, id_nu  &
          & , n_W, id_W+1, n_Snu, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav      &
          & , Slept2, mf_l, mf_nu, Chi0, c_LNSl_L, c_LNSl_R, ChiM, c_CNuSl_L &
          & , c_CNuSl_R, Sneut2, m_W, c_SlSnW, m_Z, c_SlSlZ, Smp, c_SmpSlSn  &
          & , P0, c_P0SlSl, S0, c_S0SlSl, m32, F_eff, c_GraLSl_L          &
          & , c_GraLSl_R, 2)

  Do i1=1,6
   Slept(i1)%g = Slept2(i1)%g
   Slept(i1)%gi2 = Slept2(i1)%gi2
   Slept(i1)%bi2 = Slept2(i1)%bi2
   Slept(i1)%id2 = Slept2(i1)%id2
  End Do

  !------------------------------------------------
  ! 3-body decays of sleptons
  !------------------------------------------------
  Call Slepton_3body(-1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u, n_c     &
    & , n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_S0, n_P0, n_Spm             &
    & , Slept, mZ, mW, mf_l, mf_u, mf_d, Chi0, ChiPm, P0, S0, Spm, Sneut       &
    & , c_P0SlSl, c_S0SlSl, c_SlSlZ, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn    &
    & , c_CNuSl_L, c_CNuSl_R, c_NuNuZ_L, c_NuNuZ_R, c_LLZ_L, c_LLZ_R, c_UUZ_L  &
    & , c_UUZ_R, c_DDZ_L, c_DDZ_R, c_LNuW, c_DUW, c_LLP0_L, c_LLP0_R           &
    & , c_UUP0_L, c_UUP0_R, c_DDP0_L, c_DDP0_R, c_LLS0_L, c_LLS0_R, c_UUS0_L   &
    & , c_UUS0_R, c_DDS0_L, c_DDS0_R, c_CLSn_L, c_CLSn_R, c_NuNSn_L, c_NuNSn_R &
    & , c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R, GenerationMixing, epsI)


  !------------------------------------------------
  ! 2-body decays of u-squarks
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Su, n_u, id_u, n_n, n_g, n_c, n_d, id_d &
          & , n_W, id_W, n_Sd, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sup2 &
          & , mf_u, mf_d, Chi0, c_UNSu_L, c_UNSu_R, ChiPm, c_CDSu_L  &
          & , c_CDSu_R, Sdown, m_W, c_SuSdW, m_Z, c_SuSuZ, Spm       &
          & , c_SmpSuSd, P0, c_P0SuSu, S0, c_S0SuSu, m32, F_eff   &
          & , c_GraUSu_L, c_GraUSu_R, 0, Glu, c_GUSu_L, c_GUSu_R)
  Do i1=1,6
   Sup(i1)%g = Sup2(i1)%g
   Sup(i1)%gi2 = Sup2(i1)%gi2
   Sup(i1)%bi2 = Sup2(i1)%bi2
   Sup(i1)%id2 = Sup2(i1)%id2
  End Do

  !------------------------------------------------
  ! 2-body decays of d-squarks
  !------------------------------------------------
  Call SfermionTwoBodyDecays(-1, n_Sd, n_d, id_d, n_n, n_g, n_c, n_u, id_u    &
          & , n_W, id_W+1, n_Su, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sdown &
          & , mf_d, mf_u, Chi0, c_DNSd_L, c_DNSd_R, ChiM, c_CUSd_L      &
          & , c_CUSd_R, Sup2, m_W, c_SdSuW, m_Z, c_SdSdZ, Smp           &
          & , c_SmpSdSu, P0, c_P0SdSd, S0, c_S0SdSd, m32, F_eff      &
          & , c_GraDSd_L, c_GraDSd_R, 0, Glu, c_GDSd_L, c_GDSd_R)

  !------------------------------------------------
  ! gluino decays
  !------------------------------------------------
  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(n_d, id_d, n_Sd, n_u, id_u, n_Su, id_grav, id_gl &
     & , Glu, Sdown, c_GDSd_L, c_GDSd_R, mf_d, Sup, c_GUSu_L, c_GUSu_R, mf_u &
     & , m32, Fgmsb, 0 )

   If (Glu%g.Lt.fac3*Glu%m) Then
    Glu%gi2 = 0._dp
    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .False. )

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .True. )

   End If

  Else ! calculation of 3-body decay modes is enforced

   Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
      & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
      & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
      & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
      & , c_GDSd_R, m_W, epsI, deltaM, .True. )

  End If

  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  c_S0WWvirt = c_S0WW / vev
  c_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * c_S0ZZ / vev
  c_GGS0 = 0._dp
  c_GlGlS0 = 0._dp

  Call ScalarTwoBodyDecays(-1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u    &
      & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_p0 &
      & , n_Spm, id_ph, id_gl, S0, c_S03, c_GlGlS0, c_GGS0, mf_l, c_LLS0_L     &
      & , c_LLS0_R, mf_d, c_DDS0_L, c_DDS0_R, mf_u, c_UUS0_L, c_UUS0_R, Slept  &
      & , c_S0SlSl, Sneut, c_S0SnSn, Sdown, c_S0SdSd, Sup, c_S0SuSu, Chi0      &
      & , c_NNS0_L, c_NNS0_R, ChiPm, c_CCS0_L, c_CCS0_R, m_W, c_S0WW           &
      & , c_S0WWvirt, m_Z, c_S0ZZ, c_S0ZZvirt, Spm, c_SmpS03, P0, c_P0S03      &
      & , c_P0S0Z, c_SmpS0W, Glu%m)

  c_GGP0 = 0._dp
  c_GlGlP0 = 0._dp

  Do i1=2,3
   Call PseudoscalarTwoBodyDecays(i1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d  &
       & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_sle, n_Sd, n_su, n_n, n_c, n_p0  &
       & , n_Spm, id_ph, id_gl, P0, mf_l, c_LLP0_L, c_LLP0_R, mf_d, c_DDP0_L   &
       & , c_DDP0_R, mf_u, c_UUP0_L, c_UUP0_R, Slept, c_P0SlSl, Sdown, c_P0SdSd&
       & , Sup, c_P0SuSu, Chi0, c_NNP0_L, c_NNP0_R, ChiPm, c_CCP0_L, c_CCP0_R  &
       & , Spm, c_SmpP03, S0, c_P0S03, m_Z, c_P0S0Z, m_W, c_SmpP0W, c_GlGlP0   &
       & , c_GGP0, Glu%m )
  End Do

  Call ChargedscalarTwoBodyDecays(2, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d   &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c  &
      & , n_p0, n_Spm, Spm, mf_l, c_SmpLNu_L, c_SmpLNu_R, mf_d, mf_u           &
      & , c_SmpDU_L, c_SmpDU_R, Slept, Sneut, c_SmpSlSn, Sdown, Sup, c_SmpSdSu &
      & , Chi0, ChiPm, c_SmpCN_L, c_SmpCN_R, m_W, m_Z, c_SmpZ, P0, c_SmpP03    &
      & , c_SmpP0W, S0, c_SmpS03, c_SmpS0W, 1)

  Do i1=1,2
   c_CGW_L(1) = U(i1,1) /F_eff 
   c_CGW_R(1) = V(i1,1) /F_eff 
   c_CGW_L(2) = oosqrt2 * U(i1,2) * cosb /F_eff 
   c_CGW_R(2) = oosqrt2 * V(i1,2) * sinb /F_eff 
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u      &
       & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c      &
       & , n_s0, n_p0, n_Spm, id_grav, ChiPm, Slept, c_CNuSl_L, c_CNuSl_R      &
       & , Sneut, c_CLSn_L, c_CLSn_R, mf_l, Sdown, c_CUSd_L, c_CUSd_R, mf_u    &
       & , Sup, c_CDSu_L, c_CDSu_R, mf_d, Chi0, m_W, c_CNW_L, c_CNW_R, Spm     &
       & , c_SmpCN_L, c_SmpCN_R, m_Z, c_CCZ_L, c_CCZ_R, P0, c_CCP0_L, c_CCP0_R &
       & , S0, c_CCS0_L, c_CCS0_R, m32, c_CGW_L, c_CGW_R, 1)
    k_neut =0
    
    If (ChiPm(i1)%g.Lt. fac3*Abs(ChiPm(i1)%m)) Then
     ChiPm(i1)%gi2 = 0._dp
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

    Else ! calculate 3-body decays with virtual interemdiate states only

     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .True. )
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

   End If

  End Do

  OnlySM = .False.

  Do i1=1,5
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    c_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / F_eff
    c_NGZ(1) = (- N(i1,1)*sinW+ N(i1,2)*cosW) / F_eff
    c_NGZ(2) = oosqrt2 * (cosb * N(i1,3) - sinb * N(i1,4) ) / F_eff
    c_NGH = oosqrt2 * (N(i1,3) * RS0(1,1) + N(i1,4)*RS0(1,2) + N(i1,5)*RS0(1,3) )  &
          & / F_eff
    Call NeutralinoTwoBodyDecays(i1, Chi0, n_nu, id_nu, n_l, id_l, n_d, id_d   &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W , n_snu, n_sle, n_Sd, n_su, n_n, n_c &
      & , n_s0, n_p0, n_Spm, id_ph, id_grav, Slept, c_LNSl_L, c_LNSl_R, mf_l   &
      & , Sneut, c_NuNSn_L, c_NuNSn_R, Sdown, c_DNSd_L, c_DNSd_R, mf_d, Sup    &
      & , c_UNSu_L, c_UNSu_R, mf_u, ChiPm, m_W, c_CNW_L, c_CNW_R, Spm          &
      & , c_SmpCN_L, c_SmpCN_R, m_Z, c_NNZ_L, c_NNZ_R, P0, c_NNP0_L, c_NNP0_R  &
      & , S0, c_NNS0_L, c_NNS0_R, m32, c_NGP, c_NGZ, c_NGH, 1 )

    If (Chi0(i1)%g.Lt. fac3*Abs(Chi0(i1)%m)) Then
     Chi0(i1)%gi2 = 0._dp
     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 2500, .True.)

    Else ! calculate only 3-body via virtual particles

     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .True., 0._dp, 200    &
      & , 2500, .True.)
    End If

   Else ! calculation of 3-body decay modes is enforced

    Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u  &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 2500, .True.)

    End If ! CTBD

  End Do
     
  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
  If (GenerationMixing) Then
   Write(ErrCan,*) "Warning 3-body decays of ~t_1 are not in case of"
   Write(ErrCan,*) "generation mixing"
   Do i1=1,6
    Sup(i1)%gi3 = 0
    Sup(i1)%bi3 = 0
   End Do
  Else 
   If ((Sup(5)%g.Lt.fac3*Sup(5)%m).Or.CTBD) Then ! calculation including widths
    Sup(5)%gi2 = 0._dp
    !--------------------------------------------
    ! check if flavour violating decays exist
    !--------------------------------------------
    Do i1=1,2
     Call ScalarToTwoFermions(Sup(5)%m, mf_u(2), Chi0(i1)%m &
            & , c_UNSu_L(2,i1,5), c_UNSu_R(2,i1,5), gam)
     Sup(5)%gi2(i1) = gam
     Sup(5)%id2(i1,1) = Chi0(i1)%id
     Sup(5)%id2(i1,2) = id_u(2)
    End Do
    Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d     &
       & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
       & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
       & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
       & , Slept, c_CNuSl_R, epsI, .False.)

   Else  ! calculation excluding widths

    Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d     &
       & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
       & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
       & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
       & , Slept, c_CNuSl_R, epsI, .True.)

   End If
  End If ! generation mixing

  Iname = Iname - 1

 End Subroutine CalculateBR_NMSSM

 Subroutine CalculateBR_RPeps(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u    &
    & , n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0   &
    & , n_p0, n_Spm, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, ChiPm, U, V &
    & , Chi0, N, Sup, RSup, Sdown, RSdown, uD_L, uD_R, uU_L, uU_R, S0, RS0    &
    & , P0, RP0, Spm, RSpm, epsI, deltaM, CTBD, fac3, Y_d, A_d, Y_l, A_l      &
    & , Y_u, A_u, mu, eps, vevSM, vL, Fgmsb, m32, grav_fac) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 ! 28.02.14: changing units of m32 to GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu, n_sle, n_Sd &
     & , n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm, id_grav, id_gl, id_ph 
  Integer, Intent(in), Dimension(1) :: id_Z, id_W
  Integer, Intent(in), Dimension(3) :: id_nu, id_l, id_d, id_u

  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: RP0(5,5), RS0(5,5)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(8,8), U(5,5), V(5,5), N(7,7) &
         & , RSup(6,6), RSdown(6,6), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu, eps(3)
  Real(dp), Intent(in) :: vevSM(2), Fgmsb, m32, vL(3), grav_fac
  Logical, Intent(in) :: CTBD

  Type(particle2), Intent(inout) :: Sdown(6), Spm(8), P0(5)
  Type(particle23), Intent(inout) :: Chi0(7), Sup(6), ChiPm(5), Glu, S0(5)

  Type(particle2) :: Sup2(6) ! contains only the information on the 2-body
                             ! decays of the scalar up
  Type(particle23) :: Sneut(0), Slept(0) ! dummy arguments
  Type(particle2) :: SMp(8) ! contains the identies of the negative charged ones
  Type(particle23) :: ChiM(5) ! contains the identies of the negative charged ones

  Integer :: i1, k_neut
  Real(dp) :: sinW2

 !----------------------------
 ! dummy couplings
 !----------------------------
 Complex(dp) :: c_P0SlSl(5,6,6)=0._dp, c_CNuSl_R(5,3,6)=0._dp &
      & , c_CNuSl_L(5,3,6)=0._dp, c_CLSn_R(5,3,3)=0._dp       &
      & , c_CLSn_L(5,3,3)=0._dp, c_LLP0_R(3,3,5)=0._dp        &
      & , c_LLP0_L(3,3,5)=0._dp, c_LLS0_R(3,3,5)=0._dp        &
      & , c_LLS0_L(3,3,5)=0._dp, c_NuNSn_R(3,7,3)=0._dp       &
      & , c_NuNSn_L(3,7,3)=0._dp, c_LNSl_R(3,7,6)=0._dp       &
      & , c_LNSl_L(3,7,6)=0._dp, c_SmpLNu_R(8,3,3)=0._dp      &
      & , c_SmpLNu_L(8,3,3)=0._dp, c_SmpSlSn(8,6,3)=0._dp     &
      & , c_S0SlSl(5,6,6)=0._dp, c_S0SnSn(5,3,3)=0._dp        &
      & , c_LNuW(3,3)=0._dp
 Real(dp) :: c_LLZ_L = 0._dp, c_LLZ_R = 0._dp, c_NuNuZ_L = 0._dp &
      & , c_NuNuZ_R = 0._dp
 !----------------------------
 ! couplings
 !----------------------------
  Complex(dp) :: c_SmpSdSu(8,6,6), c_SmpSuSd(8,6,6)       &
      & , c_SmpP03(8,8,5), c_SmpP0W(8,5,1), c_SmpS03(8,8,5)            &
      & , c_SmpS0W(8,5,1), c_SmpDU_L(8,3,3), c_SmpDU_R(8,3,3)          &
      & , c_SmpZ(8,8,1), c_DUW(3,3)
  Real(dp) :: c_DDZ_L = 0._dp, c_DDZ_R = 0._dp, c_UUZ_L = 0._dp      &
      & , c_UUZ_R = 0._dp
  Complex(dp) :: c_CCZ_L(5,5,1), c_CCZ_R(5,5,1)           &
      & , c_NNZ_L(7,7,1), c_NNZ_R(7,7,1), c_NNS0_L(7,7,5)            &
      & , c_NNS0_R(7,7,5), c_NNP0_L(7,7,5), c_NNP0_R(7,7,5) 
  Complex(dp) :: c_GDSd_L(3,6), c_GDSd_R(3,6)          &
      & , c_DNSd_L(3,7,6), c_DNSd_R(3,7,6), c_GUSu_L(3,6)         &
      & , c_GUSu_R(3,6), c_UNSu_L(3,7,6), c_UNSu_R(3,7,6)         &
      & , c_DDP0_L(3,3,5), c_UUP0_L(3,3,5), c_DDP0_R(3,3,5)       &
      & , c_UUP0_R(3,3,5), c_DDS0_L(3,3,5), c_UUS0_L(3,3,5)       &
      & , c_DDS0_R(3,3,5), c_UUS0_R(3,3,5)
  Complex(dp) :: c_CUSd_L(5,3,6), c_CUSd_R(5,3,6)      &
      & , c_CDSu_L(5,3,6), c_CDSu_R(5,3,6)
  Complex(dp) :: c_GlGlP0(5), c_GGP0(5), c_GlGlS0(5), c_GGS0(5)
  Complex(dp) :: c_P0SdSd(5,6,6), c_P0SuSu(5,6,6), c_P0S0Z(5,5,1) 
  Real(dp) :: c_P0S03(5,5,5)
  Complex(dp) :: c_S0SdSd(5,6,6), c_S0SuSu(5,6,6)
  Real(dp) :: c_S03(5,5,5), c_S0WW(5,1), c_S0ZZ(5,1), c_S0WWvirt(5,1) &
      & , c_S0ZZvirt(5,1)
  Complex(dp) :: c_SdSuW(6,6,1), c_SuSdW(6,6,1), c_SdSdZ(6,6,1) &
      & , c_SuSuZ(6,6,1), c_NGP, c_NGZ(2), c_NGH, c_CGW_L(2), c_CGW_R(2)
  Complex(dp) :: c_CCP0_L(5,5,5), c_CCP0_R(5,5,5)    &
      & , c_CCS0_L(5,5,5), c_CCS0_R(5,5,5), c_CNW_L(5,7,1)        &
      & , c_CNW_R(5,7,1), c_SmpCN_L(8,5,7), c_SmpCN_R(8,5,7)

  Complex(dp), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R

  Real(dp) :: tanb, sinW, cosW, vev, m_Z(1), m_W(1), F_eff, sinb, cosb
  Complex(dp) :: bi(4)
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_RPeps'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  bi(1) = mu
  bi(2:4) = eps

  m_W = mW
  m_Z = mZ
  Spm(1)%g = gamW
  P0(1)%g = gamZ

  Sup2%m = Sup%m 
  Sup2%id = Sup%id 
  ChiM%m = ChiPm%m
  ChiM%id = ChiPm%id + 1
  SMp%m = SPm%m
  SMp%id = SPm%id + 1

  tanb = vevSM(2) / vevSM(1)
  cosb = 1._dp / Sqrt(1._dp + tanb**2)
  sinb = cosb * tanb 
   
  If (grav_fac.Ne.0._dp) Then
   F_eff = grav_fac*Fgmsb ! effective SUSY breaking scale in case of GSMB
  Else 
   F_eff = Fgmsb
  End If

  Call AllCouplingsEps3(gauge, Y_l, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, vevSM  &
    &  , vL, RSpm, RP0, RS0, U, V, N, bi, PhaseGlu, A_l, RSup, A_u    &
    &  , RSdown, A_d, GenerationMixing                                  &
    & , c_SmpSdSu, c_SmpSuSd, c_SmpP03, c_SmpP0W(:,:,1), c_SmpS03          &
    & , c_SmpS0W(:,:,1), c_SmpDU_L, c_SmpDU_R, c_SmpZ(:,:,1), c_DUW, c_DDZ_L    &
    & , c_DDZ_R, c_UUZ_L, c_UUZ_R, c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_NNZ_L(:,:,1)      &
    & , c_NNZ_R(:,:,1), c_NNS0_L, c_NNS0_R, c_NNP0_L, c_NNP0_R, c_GDSd_L &
    & , c_GDSd_R, c_DNSd_L, c_DNSd_R, c_GUSu_L, c_GUSu_R            &
    & , c_UNSu_L, c_UNSu_R, c_DDP0_L, c_UUP0_L, c_DDP0_R            &
    & , c_UUP0_R, c_DDS0_L, c_UUS0_L, c_DDS0_R, c_UUS0_R            &
    & , c_CUSd_L, c_CUSd_R, c_CDSu_L, c_CDSu_R, c_GlGlS0            &
    & , c_P0SdSd, c_P0SuSu, c_P0S0Z, c_P0S03, c_S0SdSd, c_S0SuSu  &
    & , c_S03, c_S0WW(:,1), c_S0ZZ(:,1), c_SdSuW(:,:,1), c_SuSdW(:,:,1), c_SdSdZ(:,:,1)          &
    & , c_SuSuZ(:,:,1), c_CCP0_L, c_CCP0_R, c_CCS0_L, c_CCS0_R, c_CNW_L(:,:,1)  &
    & , c_CNW_R(:,:,1), c_SmpCN_L, c_SmpCN_R)

c_GraDSd_L = 0
c_GraDSd_R = 0
c_GraUSu_L = 0
c_GraUSu_R = 0

  Do i1=1,6
   Sup2(i1)%id2 = 0
  End Do
  Call SfermionTwoBodyDecays(-1, n_Su, n_u, id_u, n_n, n_g, n_c, n_d, id_d &
          & , n_W, id_W, n_Sd, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sup2 &
          & , mf_u, mf_d, Chi0, c_UNSu_L, c_UNSu_R, ChiPm, c_CDSu_L  &
          & , c_CDSu_R, Sdown, m_W, c_SuSdW, m_Z, c_SuSuZ, Spm       &
          & , c_SmpSuSd, P0, c_P0SuSu, S0, c_S0SuSu, m32, F_eff   &
          & , c_GraUSu_L, c_GraUSu_R, 0, Glu, c_GUSu_L, c_GUSu_R)
  Do i1=1,6
   Sup(i1)%g = Sup2(i1)%g
   Sup(i1)%gi2 = Sup2(i1)%gi2
   Sup(i1)%bi2 = Sup2(i1)%bi2
   Sup(i1)%id2 = Sup2(i1)%id2
   Call check_charge(Sup(i1)%id,Sup(i1)%id2)
  End Do

  Call SfermionTwoBodyDecays(-1, n_Sd, n_d, id_d, n_n, n_g, n_c, n_u, id_u    &
          & , n_W, id_W+1, n_Su, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sdown &
          & , mf_d, mf_u, Chi0, c_DNSd_L, c_DNSd_R, ChiM, c_CUSd_L      &
          & , c_CUSd_R, Sup2, m_W, c_SdSuW, m_Z, c_SdSdZ, Smp           &
          & , c_SmpSdSu, P0, c_P0SdSd, S0, c_S0SdSd, m32, F_eff      &
          & , c_GraDSd_L, c_GraDSd_R, 0, Glu, c_GDSd_L, c_GDSd_R)
  Do i1=1,6
   Call check_charge(Sdown(i1)%id,Sdown(i1)%id2)
  End Do


  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(n_d, id_d, n_Sd, n_u, id_u, n_Su, id_grav, id_gl &
     & , Glu, Sdown, c_GDSd_L, c_GDSd_R, mf_d, Sup, c_GUSu_L, c_GUSu_R, mf_u &
     & , m32, Fgmsb, 0 )
   Call check_charge(Glu%id,Glu%id2)

   If (Glu%g.Lt.fac3*Glu%m) Then
    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .False. )

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .True. )

   End If

  Else ! calculation of 3-body decay modes is enforced

   Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
      & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
      & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
      & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
      & , c_GDSd_R, m_W, epsI, deltaM, .True. )

  End If
  Call check_charge(Glu%id,Glu%id3)

  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  c_S0WWvirt = c_S0WW / vev
  c_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * c_S0ZZ / vev
  c_GGS0 = 0._dp
  c_GlGlS0 = 0._dp

  Call ScalarTwoBodyDecays(-1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u    &
      & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_p0 &
      & , n_Spm, id_ph, id_gl, S0, c_S03, c_GlGlS0, c_GGS0, mf_l, c_LLS0_L     &
      & , c_LLS0_R, mf_d, c_DDS0_L, c_DDS0_R, mf_u, c_UUS0_L, c_UUS0_R, Slept  &
      & , c_S0SlSl, Sneut, c_S0SnSn, Sdown, c_S0SdSd, Sup, c_S0SuSu, Chi0      & 
      & , c_NNS0_L, c_NNS0_R, ChiPm, c_CCS0_L, c_CCS0_R, m_W, c_S0WW           &
      & , c_S0WWvirt, m_Z, c_S0ZZ, c_S0ZZvirt, Spm, c_SmpS03, P0, c_P0S03      &
      & , c_P0S0Z, c_SmpS0W, Glu%m)

  Do i1=1,n_S0
   Call check_charge(S0(i1)%id,S0(i1)%id2)
  End Do

  c_GGP0 = 0._dp
  c_GlGlP0 = 0._dp

  Do i1=2,5
   Call PseudoscalarTwoBodyDecays(i1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d  &
       & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_sle, n_Sd, n_su, n_n, n_c, n_p0  &
       & , n_Spm, id_ph, id_gl, P0, mf_l, c_LLP0_L, c_LLP0_R, mf_d, c_DDP0_L   &
       & , c_DDP0_R, mf_u, c_UUP0_L, c_UUP0_R, Slept, c_P0SlSl, Sdown, c_P0SdSd&
       & , Sup, c_P0SuSu, Chi0, c_NNP0_L, c_NNP0_R, ChiPm, c_CCP0_L, c_CCP0_R  &
       & , Spm, c_SmpP03, S0, c_P0S03, m_Z, c_P0S0Z, m_W, c_SmpP0W, c_GlGlP0   &
       & , c_GGP0, Glu%m )
   Call check_charge(P0(i1)%id,P0(i1)%id2)
  End Do

  Do i1=2,8
   Call ChargedscalarTwoBodyDecays(i1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c  &
      & , n_p0, n_Spm, Spm, mf_l, c_SmpLNu_L, c_SmpLNu_R, mf_d, mf_u           &
      & , c_SmpDU_L, c_SmpDU_R, Slept, Sneut, c_SmpSlSn, Sdown, Sup, c_SmpSdSu &
      & , Chi0, ChiPm, c_SmpCN_L, c_SmpCN_R, m_W, m_Z, c_SmpZ, P0, c_SmpP03    &
      & , c_SmpP0W, S0, c_SmpS03, c_SmpS0W, 1)
   Call check_charge(Spm(i1)%id,SPm(i1)%id2)

!     if ((abs(Chi0(4)%m).gt.Spm(2)%m).and.(i1.le.4)) then
!      z1 = 181 
!      z2 = 191
!      Do i2=1,3
!       Do i3=i2,3
!        call Slepton_to_Slepton_ll(i1, 2, i2, i3, Spm%m, Chi0(4:7)%m &
!                & , c_SmpCN_L(:,:,4:7) , c_SmpCN_R(:,:,4:7) , 1.e-5_dp, gam)
!        gP_Spm(i1,z1) = gam(1)
!        z1 = z1 + 1
!        if (i2.ne.i3) gP_Spm(i1,z1) = gam(1)
!        if (i2.ne.i3) z1 = z1 + 1
!        gP_Spm(i1,z2) = gam(2)
!        z2 = z2 + 1
!       end do
!      end do
!      gT_Spm(i1) = Sum(gP_Spm(i1,:))
!      BR_Spm(i1,:) = gP_Spm(i1,:) / gT_Spm(i1)
!     end if

  End Do

  Do i1=4,5
   c_CGW_L(1) = U(i1,1) /F_eff 
   c_CGW_R(1) = V(i1,1) /F_eff 
   c_CGW_L(2) = oosqrt2 * U(i1,2) * cosb / F_eff
   c_CGW_R(2) = oosqrt2 * V(i1,2) * sinb / F_eff
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u      &
       & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c      &
       & , n_s0, n_p0, n_Spm, id_grav, ChiPm, Slept, c_CNuSl_L, c_CNuSl_R      &
       & , Sneut, c_CLSn_L, c_CLSn_R, mf_l, Sdown, c_CUSd_L, c_CUSd_R, mf_u    &
       & , Sup, c_CDSu_L, c_CDSu_R, mf_d, Chi0, m_W, c_CNW_L, c_CNW_R, Spm     &
       & , c_SmpCN_L, c_SmpCN_R, m_Z, c_CCZ_L, c_CCZ_R, P0, c_CCP0_L, c_CCP0_R &
       & , S0, c_CCS0_L, c_CCS0_R, m32, c_CGW_L, c_CGW_R, 1)
    k_neut =0
    Call check_charge(ChiPm(i1)%id,ChiPm(i1)%id2)

    If (ChiPm(i1)%g.Lt. fac3*Abs(ChiPm(i1)%m)) Then
     ChiPm(i1)%gi2 = 0._dp
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

    Else ! calculate 3-body decays with virtual interemdiate states only

     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .True. )
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

   End If
   Call check_charge(ChiPm(i1)%id,ChiPm(i1)%id3)

  End Do

  OnlySM = .True.

  Do i1=4,7
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    c_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / F_eff
    c_NGZ(1) = (- N(i1,1)*sinW+ N(i1,2)*cosW) / F_eff
    c_NGZ(2) = oosqrt2 * (cosb * N(i1,3) - sinb * N(i1,4) ) / F_eff
    c_NGH =  oosqrt2 * (N(i1,3) * RS0(1,1) + N(i1,4)*RS0(1,2)) / F_eff
    Call NeutralinoTwoBodyDecays(i1, Chi0, n_nu, id_nu, n_l, id_l, n_d, id_d   &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W , n_snu, n_sle, n_Sd, n_su, n_n, n_c &
      & , n_s0, n_p0, n_Spm, id_ph, id_grav, Slept, c_LNSl_L, c_LNSl_R, mf_l   &
      & , Sneut, c_NuNSn_L, c_NuNSn_R, Sdown, c_DNSd_L, c_DNSd_R, mf_d, Sup    &
      & , c_UNSu_L, c_UNSu_R, mf_u, ChiPm, m_W, c_CNW_L, c_CNW_R, Spm          &
      & , c_SmpCN_L, c_SmpCN_R, m_Z, c_NNZ_L, c_NNZ_R, P0, c_NNP0_L, c_NNP0_R  &
      & , S0, c_NNS0_L, c_NNS0_R, m32, c_NGP, c_NGZ, c_NGH, 1 )
    Call check_charge(Chi0(i1)%id,Chi0(i1)%id2)

    If (Chi0(i1)%g.Lt. fac3*Abs(Chi0(i1)%m)) Then
     Chi0(i1)%gi2 = 0._dp
     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 4000, .False.)

    Else ! calculate only 3-body via virtual particles

     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .True., 0._dp, 200    &
      & , 4000, .False.)

    End If

   Else ! calculation of 3-body decay modes is enforced

    Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u  &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 4000, .False.)

   End If ! CTBD
   Call check_charge(Chi0(i1)%id,Chi0(i1)%id3)

  End Do

  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
!   If (GenerationMixing) Then
!    Write(ErrCan,*) "Warning 3-body decays of ~t_1 are not in case of"
!    Write(ErrCan,*) "generation mixing"
   Do i1=1,6
    Sup(i1)%gi3 = 0
    Sup(i1)%bi3 = 0
   End Do
!   Else 
!    If ((Sup(5)%g.Lt.fac3*Sup(5)%m).Or.CTBD) Then ! calculation including widths
!     Sup(5)%gi2 = 0._dp
!     !--------------------------------------------
!     ! check if flavour violating decays exist
!     !--------------------------------------------
!     Do i1=4,5 ! the first three are neutrinos
!      Call ScalarToTwoFermions(Sup(5)%m, mf_u(2), Chi0(i1)%m &
!             & , c_UNSu_L(2,i1,5), c_UNSu_R(2,i1,5), gam)
!      Sup(5)%gi2(i1) = gam
!      Sup(5)%id2(i1,1) = Chi0(i1)%id
!      Sup(5)%id2(i1,2) = id_u(2)
!     End Do
!     Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d  &
!        & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0(4:7), gauge(2)             &
!        & , c_UNSu_L(:,4:7,:), c_UNSu_R(:,4:7,:), c_DNSd_L(:,4:7,:)          &
!        & , c_DNSd_R(:,4:7,:), c_SdSuW(:,:,1), Sneut, ChiPm(4:5), c_CLSn_L   &
!        & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)     &
!        & , Slept, c_CNuSl_R, epsI, .False.)
! 
!   Else  ! calculation excluding widths
!     Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d     &
!        & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
!        & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
!        & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
!        & , Slept, c_CNuSl_R, epsI, .True.)
! 
!   End If
! 
!   End If ! generation mixing

  Iname = Iname - 1

 End Subroutine CalculateBR_RPeps

 Subroutine CalculateBR_RPlam(n_nu, id_nu, n_l, id_l, n_d, id_d, n_u, id_u    &
    & , n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_g, n_s0   &
    & , n_p0, n_Spm, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu, ChiPm, U, V &
    & , Chi0, N, Sup, RSup, Sdown, RSdown, uD_L, uD_R, uU_L, uU_R, S0, RS0    &
    & , P0, RP0, Spm, RSpm, epsI, deltaM, CTBD, fac3, Y_d, A_d, Y_l, A_l      &
    & , Y_u, A_u, mu, eps, RP_lam, RP_lamp, vevSM, vL, Fgmsb, m32, grav_fac) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 ! 28.02.14: changing units of m32 to GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu, n_sle, n_Sd &
     & , n_su, n_n, n_c, n_g, n_s0, n_p0, n_Spm, id_grav, id_gl, id_ph 
  Integer, Intent(in), Dimension(1) :: id_Z, id_W
  Integer, Intent(in), Dimension(3) :: id_nu, id_l, id_d, id_u

  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: RP0(5,5), RS0(5,5)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(8,8), U(5,5), V(5,5), N(7,7) &
         & , RSup(6,6), RSdown(6,6), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu, eps(3), RP_lam(3,3,3), RP_lamp(3,3,3)
  Real(dp), Intent(in) :: vevSM(2), Fgmsb, m32, vL(3), grav_fac
  Logical, Intent(in) :: CTBD

  Type(particle2), Intent(inout) :: Sdown(6), Spm(8), P0(5)
  Type(particle23), Intent(inout) :: Chi0(7), Sup(6), ChiPm(5), Glu, S0(5)

  Type(particle2) :: Sup2(6) ! contains only the information on the 2-body
                             ! decays of the scalar up
  Type(particle23) :: Sneut(0), Slept(0) ! dummy arguments
  Type(particle2) :: SMp(8) ! contains the identies of the negative charged ones
  Type(particle23) :: ChiM(5) ! contains the identies of the negative charged ones

  Integer :: i1, k_neut
  Real(dp) :: sinW2

 !----------------------------
 ! dummy couplings
 !----------------------------
 Complex(dp) :: c_P0SlSl(5,6,6)=0._dp, c_CNuSl_R(5,3,6)=0._dp &
      & , c_CNuSl_L(5,3,6)=0._dp, c_CLSn_R(5,3,3)=0._dp       &
      & , c_CLSn_L(5,3,3)=0._dp, c_LLP0_R(3,3,5)=0._dp        &
      & , c_LLP0_L(3,3,5)=0._dp, c_LLS0_R(3,3,5)=0._dp        &
      & , c_LLS0_L(3,3,5)=0._dp, c_NuNSn_R(3,7,3)=0._dp       &
      & , c_NuNSn_L(3,7,3)=0._dp, c_LNSl_R(3,7,6)=0._dp       &
      & , c_LNSl_L(3,7,6)=0._dp, c_SmpLNu_R(8,3,3)=0._dp      &
      & , c_SmpLNu_L(8,3,3)=0._dp, c_SmpSlSn(8,6,3)=0._dp     &
      & , c_S0SlSl(5,6,6)=0._dp, c_S0SnSn(5,3,3)=0._dp        &
      & , c_LNuW(3,3)=0._dp
 Real(dp) :: c_LLZ_L = 0._dp, c_LLZ_R = 0._dp, c_NuNuZ_L = 0._dp &
      & , c_NuNuZ_R = 0._dp
 !----------------------------
 ! couplings
 !----------------------------
  Complex(dp) :: c_SmpSdSu(8,6,6), c_SmpSuSd(8,6,6)       &
      & , c_SmpP03(8,8,5), c_SmpP0W(8,5,1), c_SmpS03(8,8,5)            &
      & , c_SmpS0W(8,5,1), c_SmpDU_L(8,3,3), c_SmpDU_R(8,3,3)          &
      & , c_SmpZ(8,8,1), c_DUW(3,3)
  Real(dp) :: c_DDZ_L = 0._dp, c_DDZ_R = 0._dp, c_UUZ_L = 0._dp      &
      & , c_UUZ_R = 0._dp
  Complex(dp) :: c_CCZ_L(5,5,1), c_CCZ_R(5,5,1)           &
      & , c_NNZ_L(7,7,1), c_NNZ_R(7,7,1), c_NNS0_L(7,7,5)            &
      & , c_NNS0_R(7,7,5), c_NNP0_L(7,7,5), c_NNP0_R(7,7,5) 
  Complex(dp) :: c_GDSd_L(3,6), c_GDSd_R(3,6)          &
      & , c_DNSd_L(3,7,6), c_DNSd_R(3,7,6), c_GUSu_L(3,6)         &
      & , c_GUSu_R(3,6), c_UNSu_L(3,7,6), c_UNSu_R(3,7,6)         &
      & , c_DDP0_L(3,3,5), c_UUP0_L(3,3,5), c_DDP0_R(3,3,5)       &
      & , c_UUP0_R(3,3,5), c_DDS0_L(3,3,5), c_UUS0_L(3,3,5)       &
      & , c_DDS0_R(3,3,5), c_UUS0_R(3,3,5)
  Complex(dp) :: c_CUSd_L(5,3,6), c_CUSd_R(5,3,6)      &
      & , c_CDSu_L(5,3,6), c_CDSu_R(5,3,6)
  Complex(dp) :: c_GlGlP0(5), c_GGP0(5), c_GlGlS0(5), c_GGS0(5)
  Complex(dp) :: c_P0SdSd(5,6,6), c_P0SuSu(5,6,6), c_P0S0Z(5,5,1) 
  Real(dp) :: c_P0S03(5,5,5)
  Complex(dp) :: c_S0SdSd(5,6,6), c_S0SuSu(5,6,6)
  Real(dp) :: c_S03(5,5,5), c_S0WW(5,1), c_S0ZZ(5,1), c_S0WWvirt(5,1) &
      & , c_S0ZZvirt(5,1)
  Complex(dp) :: c_SdSuW(6,6,1), c_SuSdW(6,6,1), c_SdSdZ(6,6,1) &
      & , c_SuSuZ(6,6,1), c_NGP, c_NGZ(2), c_NGH, c_CGW_L(2), c_CGW_R(2)
  Complex(dp) :: c_CCP0_L(5,5,5), c_CCP0_R(5,5,5)    &
      & , c_CCS0_L(5,5,5), c_CCS0_R(5,5,5), c_CNW_L(5,7,1)        &
      & , c_CNW_R(5,7,1), c_SmpCN_L(8,5,7), c_SmpCN_R(8,5,7)

  Complex(dp), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R

  Real(dp) :: tanb, sinW, cosW, vev, m_Z(1), m_W(1), F_eff, sinb, cosb
  Complex(dp) :: bi(4)
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_RPlam'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  bi(1) = mu
  bi(2:4) = eps

  m_W = mW
  m_Z = mZ
  Spm(1)%g = gamW
  P0(1)%g = gamZ

  Sup2%m = Sup%m 
  Sup2%id = Sup%id 
  ChiM%m = ChiPm%m
  ChiM%id = ChiPm%id + 1
  SMp%m = SPm%m
  SMp%id = SPm%id + 1

  tanb = vevSM(2) / vevSM(1)
  cosb = 1._dp / Sqrt(1._dp + tanb**2)
  sinb = cosb * tanb 
   
  If (grav_fac.Ne.0._dp) Then
   F_eff = grav_fac*Fgmsb ! effective SUSY breaking scale in case of GSMB
  Else 
   F_eff = Fgmsb
  End If

  Call AllCouplingsLam(gauge, Y_l, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, A_l, RSup &
    & , A_u, RSdown, A_d, vevSM, vL, RP_lam, RP_lamp, RSpm, RP0, RS0, U, V, N  &
    & , bi, PhaseGlu, GenerationMixing                                         &
    & , c_SmpSdSu, c_SmpSuSd, c_SmpP03, c_SmpP0W(:,:,1), c_SmpS03              &
    & , c_SmpS0W(:,:,1), c_SmpDU_L, c_SmpDU_R, c_SmpZ(:,:,1), c_DUW, c_DDZ_L   &
    & , c_DDZ_R, c_UUZ_L, c_UUZ_R, c_CCZ_L(:,:,1), c_CCZ_R(:,:,1)              &
    & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_NNS0_L, c_NNS0_R, c_NNP0_L, c_NNP0_R &
    & , c_GDSd_L, c_GDSd_R, c_DNSd_L, c_DNSd_R, c_GUSu_L, c_GUSu_R, c_UNSu_L   &
    & , c_UNSu_R, c_DDP0_L, c_UUP0_L, c_DDP0_R, c_UUP0_R, c_DDS0_L, c_UUS0_L   &
    & , c_DDS0_R, c_UUS0_R, c_CUSd_L, c_CUSd_R, c_CDSu_L, c_CDSu_R, c_GlGlS0   &
    & , c_P0SdSd, c_P0SuSu, c_P0S0Z, c_P0S03, c_S0SdSd, c_S0SuSu, c_S03        &
    & , c_S0WW(:,1), c_S0ZZ(:,1), c_SdSuW(:,:,1), c_SuSdW(:,:,1)               &
    & , c_SdSdZ(:,:,1), c_SuSuZ(:,:,1), c_CCP0_L, c_CCP0_R, c_CCS0_L, c_CCS0_R &
    & , c_CNW_L(:,:,1), c_CNW_R(:,:,1), c_SmpCN_L, c_SmpCN_R)

  c_GraDSd_L = 0
  c_GraDSd_R = 0
  c_GraUSu_L = 0
  c_GraUSu_R = 0

  Do i1=1,6
   Sup2(i1)%id2 = 0
  End Do
  Call SfermionTwoBodyDecays(-1, n_Su, n_u, id_u, n_n, n_g, n_c, n_d, id_d &
          & , n_W, id_W, n_Sd, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sup2 &
          & , mf_u, mf_d, Chi0, c_UNSu_L, c_UNSu_R, ChiPm, c_CDSu_L  &
          & , c_CDSu_R, Sdown, m_W, c_SuSdW, m_Z, c_SuSuZ, Spm       &
          & , c_SmpSuSd, P0, c_P0SuSu, S0, c_S0SuSu, m32, F_eff   &
          & , c_GraUSu_L, c_GraUSu_R, 0, Glu, c_GUSu_L, c_GUSu_R)
  Do i1=1,6
   Sup(i1)%g = Sup2(i1)%g
   Sup(i1)%gi2 = Sup2(i1)%gi2
   Sup(i1)%bi2 = Sup2(i1)%bi2
   Sup(i1)%id2 = Sup2(i1)%id2
   Call check_charge(Sup(i1)%id,Sup(i1)%id2)
  End Do

  Call SfermionTwoBodyDecays(-1, n_Sd, n_d, id_d, n_n, n_g, n_c, n_u, id_u    &
          & , n_W, id_W+1, n_Su, n_Spm, n_Z, id_Z, n_P0, n_S0, id_grav, Sdown &
          & , mf_d, mf_u, Chi0, c_DNSd_L, c_DNSd_R, ChiM, c_CUSd_L      &
          & , c_CUSd_R, Sup2, m_W, c_SdSuW, m_Z, c_SdSdZ, Smp           &
          & , c_SmpSdSu, P0, c_P0SdSd, S0, c_S0SdSd, m32, F_eff      &
          & , c_GraDSd_L, c_GraDSd_R, 0, Glu, c_GDSd_L, c_GDSd_R)
  Do i1=1,6
   Call check_charge(Sdown(i1)%id,Sdown(i1)%id2)
  End Do


  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(n_d, id_d, n_Sd, n_u, id_u, n_Su, id_grav, id_gl &
     & , Glu, Sdown, c_GDSd_L, c_GDSd_R, mf_d, Sup, c_GUSu_L, c_GUSu_R, mf_u &
     & , m32, Fgmsb, 0 )
   Call check_charge(Glu%id,Glu%id2)

   If (Glu%g.Lt.fac3*Glu%m) Then
    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .False. )

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
       & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
       & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
       & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
       & , c_GDSd_R, m_W, epsI, deltaM, .True. )

   End If

  Else ! calculation of 3-body decay modes is enforced

   Call GluinoThreeBodyDecays(n_d, id_d, n_u, id_u, n_n, n_c, id_gl, n_su  &
      & , n_sd, n_W, Glu, Chi0, ChiPm, mf_u, g_T, mf_d, Sup, gauge(3)      &
      & , c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R, c_GUSu_L, c_GUSu_R       &
      & , c_SdSuW, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L, c_CUSd_R, c_GDSd_L &
      & , c_GDSd_R, m_W, epsI, deltaM, .True. )

  End If
  Call check_charge(Glu%id,Glu%id3)

  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  c_S0WWvirt = c_S0WW / vev
  c_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * c_S0ZZ / vev
  c_GGS0 = 0._dp
  c_GlGlS0 = 0._dp

  Call ScalarTwoBodyDecays(-1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u    &
      & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c, n_p0 &
      & , n_Spm, id_ph, id_gl, S0, c_S03, c_GlGlS0, c_GGS0, mf_l, c_LLS0_L     &
      & , c_LLS0_R, mf_d, c_DDS0_L, c_DDS0_R, mf_u, c_UUS0_L, c_UUS0_R, Slept  &
      & , c_S0SlSl, Sneut, c_S0SnSn, Sdown, c_S0SdSd, Sup, c_S0SuSu, Chi0      & 
      & , c_NNS0_L, c_NNS0_R, ChiPm, c_CCS0_L, c_CCS0_R, m_W, c_S0WW           &
      & , c_S0WWvirt, m_Z, c_S0ZZ, c_S0ZZvirt, Spm, c_SmpS03, P0, c_P0S03      &
      & , c_P0S0Z, c_SmpS0W, Glu%m)

  Do i1=1,n_S0
   Call check_charge(S0(i1)%id,S0(i1)%id2)
  End Do

  c_GGP0 = 0._dp
  c_GlGlP0 = 0._dp

  Do i1=2,5
   Call PseudoscalarTwoBodyDecays(i1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d  &
       & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_sle, n_Sd, n_su, n_n, n_c, n_p0  &
       & , n_Spm, id_ph, id_gl, P0, mf_l, c_LLP0_L, c_LLP0_R, mf_d, c_DDP0_L   &
       & , c_DDP0_R, mf_u, c_UUP0_L, c_UUP0_R, Slept, c_P0SlSl, Sdown, c_P0SdSd&
       & , Sup, c_P0SuSu, Chi0, c_NNP0_L, c_NNP0_R, ChiPm, c_CCP0_L, c_CCP0_R  &
       & , Spm, c_SmpP03, S0, c_P0S03, m_Z, c_P0S0Z, m_W, c_SmpP0W, c_GlGlP0   &
       & , c_GGP0, Glu%m )
   Call check_charge(P0(i1)%id,P0(i1)%id2)
  End Do

  Do i1=2,8
   Call ChargedscalarTwoBodyDecays(i1, n_s0, n_nu, id_nu, n_l, id_l, n_d, id_d &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c  &
      & , n_p0, n_Spm, Spm, mf_l, c_SmpLNu_L, c_SmpLNu_R, mf_d, mf_u           &
      & , c_SmpDU_L, c_SmpDU_R, Slept, Sneut, c_SmpSlSn, Sdown, Sup, c_SmpSdSu &
      & , Chi0, ChiPm, c_SmpCN_L, c_SmpCN_R, m_W, m_Z, c_SmpZ, P0, c_SmpP03    &
      & , c_SmpP0W, S0, c_SmpS03, c_SmpS0W, 1)
   Call check_charge(Spm(i1)%id,SPm(i1)%id2)

!     if ((abs(Chi0(4)%m).gt.Spm(2)%m).and.(i1.le.4)) then
!      z1 = 181 
!      z2 = 191
!      Do i2=1,3
!       Do i3=i2,3
!        call Slepton_to_Slepton_ll(i1, 2, i2, i3, Spm%m, Chi0(4:7)%m &
!                & , c_SmpCN_L(:,:,4:7) , c_SmpCN_R(:,:,4:7) , 1.e-5_dp, gam)
!        gP_Spm(i1,z1) = gam(1)
!        z1 = z1 + 1
!        if (i2.ne.i3) gP_Spm(i1,z1) = gam(1)
!        if (i2.ne.i3) z1 = z1 + 1
!        gP_Spm(i1,z2) = gam(2)
!        z2 = z2 + 1
!       end do
!      end do
!      gT_Spm(i1) = Sum(gP_Spm(i1,:))
!      BR_Spm(i1,:) = gP_Spm(i1,:) / gT_Spm(i1)
!     end if

  End Do

  Do i1=4,5
   c_CGW_L(1) = U(i1,1) /F_eff 
   c_CGW_R(1) = V(i1,1) /F_eff 
   c_CGW_L(2) = oosqrt2 * U(i1,2) * cosb / F_eff
   c_CGW_R(2) = oosqrt2 * V(i1,2) * sinb / F_eff
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, n_nu, id_nu, n_l, id_l, n_d, id_d, n_u      &
       & , id_u, n_Z, id_Z, n_W, id_W, n_snu, n_sle, n_Sd, n_su, n_n, n_c      &
       & , n_s0, n_p0, n_Spm, id_grav, ChiPm, Slept, c_CNuSl_L, c_CNuSl_R      &
       & , Sneut, c_CLSn_L, c_CLSn_R, mf_l, Sdown, c_CUSd_L, c_CUSd_R, mf_u    &
       & , Sup, c_CDSu_L, c_CDSu_R, mf_d, Chi0, m_W, c_CNW_L, c_CNW_R, Spm     &
       & , c_SmpCN_L, c_SmpCN_R, m_Z, c_CCZ_L, c_CCZ_R, P0, c_CCP0_L, c_CCP0_R &
       & , S0, c_CCS0_L, c_CCS0_R, m32, c_CGW_L, c_CGW_R, 1)

    k_neut =0
    Call check_charge(ChiPm(i1)%id,ChiPm(i1)%id2)

    If (ChiPm(i1)%g.Lt. fac3*Abs(ChiPm(i1)%m)) Then
     ChiPm(i1)%gi2 = 0._dp
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

    Else ! calculate 3-body decays with virtual interemdiate states only

     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .True. )
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u   &
       & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su      &
       & , n_S0, n_P0, n_Spm, ChiPm, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l      &
       & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R    &
       & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), Chi0, mW, gamW, c_LNuW, c_DUW       &
       & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm &
       & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R  &
       & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R      &
       & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R      &
       & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R          &
       & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R     &
       & , Glu, c_GUSu_L, c_GUSu_R, Sdown,  c_DNSd_L, c_DNSd_R, c_CUSd_L       &
       & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L &
       & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R           &
       & , GenerationMixing, k_neut, epsI, deltaM, .False. )

   End If
   Call check_charge(ChiPm(i1)%id,ChiPm(i1)%id3)

  End Do

  OnlySM = .True.

  Do i1=4,7
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    c_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / F_eff
    c_NGZ(1) = (- N(i1,1)*sinW+ N(i1,2)*cosW) / F_eff
    c_NGZ(2) = oosqrt2 * (cosb * N(i1,3) - sinb * N(i1,4) ) / F_eff
    c_NGH =  oosqrt2 * (N(i1,3) * RS0(1,1) + N(i1,4)*RS0(1,2)) / F_eff
    Call NeutralinoTwoBodyDecays(i1, Chi0, n_nu, id_nu, n_l, id_l, n_d, id_d   &
      & , n_u, id_u, n_Z, id_Z, n_W, id_W , n_snu, n_sle, n_Sd, n_su, n_n, n_c &
      & , n_s0, n_p0, n_Spm, id_ph, id_grav, Slept, c_LNSl_L, c_LNSl_R, mf_l   &
      & , Sneut, c_NuNSn_L, c_NuNSn_R, Sdown, c_DNSd_L, c_DNSd_R, mf_d, Sup    &
      & , c_UNSu_L, c_UNSu_R, mf_u, ChiPm, m_W, c_CNW_L, c_CNW_R, Spm          &
      & , c_SmpCN_L, c_SmpCN_R, m_Z, c_NNZ_L, c_NNZ_R, P0, c_NNP0_L, c_NNP0_R  &
      & , S0, c_NNS0_L, c_NNS0_R, m32, c_NGP, c_NGZ, c_NGH, 1 )
    Call check_charge(Chi0(i1)%id,Chi0(i1)%id2)

    If (Chi0(i1)%g.Lt. fac3*Abs(Chi0(i1)%m)) Then
     Chi0(i1)%gi2 = 0._dp
     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 4000, .False.)

    Else ! calculate only 3-body via virtual particles

     Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .True., 0._dp, 200    &
      & , 4000, .False.)

    End If

   Else ! calculation of 3-body decay modes is enforced

    Call NeutralinoThreeBodyDecays(i1, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u  &
      & , id_u, n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sle, n_Snu, n_Sd, n_Su, n_S0 &
      & , n_P0, n_Spm, id_ph, Chi0, mZ, gamZ, c_NuNuZ_L, c_NuNuZ_R, mf_l       &
      & , c_LLZ_L, c_LLZ_R, mf_u, c_UUZ_L, c_UUZ_R, mf_d, c_DDZ_L, c_DDZ_R     &
      & , c_NNZ_L(:,:,1), c_NNZ_R(:,:,1), ChiPm, mW, gamW, c_LNuW, c_DUW       &
      & , c_CCZ_L(:,:,1), c_CCZ_R(:,:,1), c_CNW_L(:,:,1), c_CNW_R(:,:,1), Spm  &
      & , c_SmpCN_L, c_SmpCN_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R   &
      & , S0, c_NNS0_L, c_NNS0_R, c_CCS0_L, c_CCS0_R, c_DDS0_L, c_DDS0_R       &
      & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, P0, c_NNP0_L, c_NNP0_R       &
      & , c_CCP0_L, c_CCP0_R, c_DDP0_L, c_DDP0_R, c_LLP0_L, c_LLP0_R           &
      & , c_UUP0_L, c_UUP0_R, Sup, c_UNSu_L, c_UNSu_R, c_CDSu_L, c_CDSu_R      &
      & , Glu, c_GUSu_L, c_GUSu_R, Sdown, c_DNSd_L, c_DNSd_R, c_CUSd_L         &
      & , c_CUSd_R, c_GDSd_L, c_GDSd_R, Sneut, c_NuNSn_L, c_NuNSn_R, c_CLSn_L  &
      & , c_CLSn_R, Slept, c_LNSl_L, c_LNSl_R, c_CNuSl_L, c_CNuSl_R, gauge(2)  &
      & , sinW2, GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp, 200   &
      & , 4000, .False.)

   End If ! CTBD
   Call check_charge(Chi0(i1)%id,Chi0(i1)%id3)

  End Do

  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
!   If (GenerationMixing) Then
!    Write(ErrCan,*) "Warning 3-body decays of ~t_1 are not in case of"
!    Write(ErrCan,*) "generation mixing"
   Do i1=1,6
    Sup(i1)%gi3 = 0
    Sup(i1)%bi3 = 0
   End Do
!   Else 
!    If ((Sup(5)%g.Lt.fac3*Sup(5)%m).Or.CTBD) Then ! calculation including widths
!     Sup(5)%gi2 = 0._dp
!     !--------------------------------------------
!     ! check if flavour violating decays exist
!     !--------------------------------------------
!     Do i1=4,5 ! the first three are neutrinos
!      Call ScalarToTwoFermions(Sup(5)%m, mf_u(2), Chi0(i1)%m &
!             & , c_UNSu_L(2,i1,5), c_UNSu_R(2,i1,5), gam)
!      Sup(5)%gi2(i1) = gam
!      Sup(5)%id2(i1,1) = Chi0(i1)%id
!      Sup(5)%id2(i1,2) = id_u(2)
!     End Do
!     Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d  &
!        & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0(4:7), gauge(2)             &
!        & , c_UNSu_L(:,4:7,:), c_UNSu_R(:,4:7,:), c_DNSd_L(:,4:7,:)          &
!        & , c_DNSd_R(:,4:7,:), c_SdSuW(:,:,1), Sneut, ChiPm(4:5), c_CLSn_L   &
!        & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)     &
!        & , Slept, c_CNuSl_R, epsI, .False.)
! 
!   Else  ! calculation excluding widths
!     Call StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d     &
!        & , id_d, id_W, n_n, n_c,Sup, Sdown, Chi0, gauge(2), c_UNSu_L, c_UNSu_R &
!        & , c_DNSd_L, c_DNSd_R, c_SdSuW(:,:,1), Sneut, ChiPm, c_CLSn_L          &
!        & , c_CLSn_R, c_CDSu_L, c_CDSu_R, c_CNW_L(:,:,1), c_CNW_R(:,:,1)        &
!        & , Slept, c_CNuSl_R, epsI, .True.)
! 
!   End If
! 
!   End If ! generation mixing

  Iname = Iname - 1

 End Subroutine CalculateBR_RPlam

End Module BranchingRatios

