Module LoopCouplings

Use Control
Use Mathematics
Use LoopFunctions
Use RGEs
Use StandardModel
Use Couplings

Logical, Save :: MZ_input = .False.
 Real(dp), Private :: vevSM(2)

Contains

 Real(dp) Function AlphaEwDR(Q, mHpm, mSqU, mSqD, mSl, mChar)
 !-----------------------------------------------------------------------
 ! This function calculates the electromagnetig alpha(Q) in the DR scheme.
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B 491, 3, (1997)
 ! The input is:
 !  - Alpha at Q=0, something like 1/137
 !  - Q the energy
 !  - mW, the W-boson mass
 !  - mT, the top mass
 !  - mHpm, the mass of the charged Higgs boson
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 !  - mSl(i), i=1,..,6 the masses of the sleptons
 !  - mChar(i), i=1,2 the masses of the charginos
 ! written by Werner Porod, 19.7.1999
 ! 16.11.01: - writting it such that R-parity violation is included
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Q, mHpm(:), mSqU(:), mSqD(:), mSl(:), mChar(:)

  Integer :: i, n_C, n_H
  Real(dp) :: DeltaAlpha,sumI(4)

  n_H = Size( mHpm )
  n_C = Size( mChar )

  sumI = 0._dp
  Do i=1,6
   sumI(1) = sumI(1) + Log( mSqU(i) / Q )
   sumI(2) = sumI(2) + Log( mSqD(i) / Q )
  End Do

  Do i=2,n_H
   sumI(3) = sumI(3) + Log( mHpm(i) / Q )
  End Do
  Do i=1,(5-n_C)*2
   sumI(3) = sumI(3) + Log( mSl(i) / Q )
  End Do
  Do i=n_C-1,n_C
   sumI(4) = sumI(4) + Log( mChar(i) / Q )
  End Do

  If (MZ_input) then
   DeltaAlpha = 1._dp - Alpha / Alpha_MZ_MS  ! MSbar value ^= mW + light fermions
   DeltaAlpha = DeltaAlpha + alpha / (6._dp * Pi) ! conversion to DRbar
   DeltaAlpha = DeltaAlpha - ( 16._dp * Log(mf_u(3) / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp ) * Alpha / (2._dp * pi) ! SUSY contr.
  else
   DeltaAlpha = 7._dp * Log(Q / mW) + 16._dp * Log(mf_u(3) / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp
   DeltaAlpha = Delta_Alpha_Lepton + Delta_Alpha_Hadron &
            & - alpha * DeltaAlpha / ( 2._dp * Pi)
  End If

  AlphaEwDR = Alpha / (1._dp - DeltaAlpha)

 End Function AlphaEwDR


 Real(dp) Function AlphaEwMS(Q, alphaDR, mHpm, mSqU, mSqD, mSl, mChar, mt)
 !-----------------------------------------------------------------------
 ! This function calculates the electromagnetig alpha(Q) in the MS scheme
 ! starting from the DRbar value
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B 491, 3, (1997)
 ! The input is:
 !  - Alpha at Q=0, something like 1/137
 !  - Q the energy
 !  - mW, the W-boson mass
 !  - mT, the top mass
 !  - mHpm, the mass of the charged Higgs boson
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 !  - mSl(i), i=1,..,6 the masses of the sleptons
 !  - mChar(i), i=1,2 the masses of the charginos
 !  - mt, top-quark mass
 ! written by Werner Porod, 05.02.2010
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Q, alphaDR, mHpm(:), mSqU(:), mSqD(:), mSl(:) &
                      & , mChar(:), mt

  Integer :: i, n_C, n_H
  Real(dp) :: DeltaAlpha,sumI(4)

  n_H = Size( mHpm )
  n_C = Size( mChar )

  sumI = 0._dp
  Do i=1,6
   sumI(1) = sumI(1) + Log( mSqU(i) / Q )
   sumI(2) = sumI(2) + Log( mSqD(i) / Q )
  End Do

  Do i=2,n_H
   sumI(3) = sumI(3) + Log( mHpm(i) / Q )
  End Do
  Do i=1,(5-n_C)*2
   sumI(3) = sumI(3) + Log( mSl(i) / Q )
  End Do
  Do i=n_C-1,n_C
   sumI(4) = sumI(4) + Log( mChar(i) / Q )
  End Do

  DeltaAlpha = 1._dp / (6._dp * Pi) ! conversion from DRbar
  DeltaAlpha = DeltaAlpha - ( 16._dp * Log(mt / Q ) / 9._dp &
            & + sumI(3) / 3._dp + (4._dp * sumI(1) + sumI(2) ) / 9._dp    &
            & + 4._dp * sumI(4) / 3._dp )  / (2._dp * pi) ! SUSY contr.

  AlphaEwMS = AlphaDR / (1._dp + AlphaDR * DeltaAlpha)

 End Function AlphaEwMS


 Real(dp) Function Alpha_MSbar(Q, mW, mt)
 Implicit None
  Real(dp), Intent(in) :: Q, mW
  Real(dp), Intent(in), Optional :: mt

  Real(dp) :: DeltaAlpha

  If (MZ_input) Then
   Alpha_MSbar = Alpha_mZ_MS
   If (Present(mt)) Then
    DeltaAlpha = - 8._dp * Log(Q / mt) / (9._dp * Pi)
    Alpha_MSbar = Alpha_MSbar / ( 1._dp + DeltaAlpha * alpha)
   End If 
  Else
   DeltaAlpha = 3.5_dp * Log(Q / mW) / Pi + 0.5_dp * oo3pi
   If (Present(mt)) DeltaAlpha = DeltaAlpha - 8._dp * Log(Q / mt) / (9._dp * Pi)
   Alpha_MSbar = Alpha / ( 1._dp - Delta_Alpha_Lepton - Delta_Alpha_Hadron  &
              &         + DeltaAlpha * alpha)
   Alpha_MZ_MS = Alpha_MSbar
  End If

 End Function Alpha_MSbar


 Real(dp) Function AlphaSDR(Q, mG, mSqU, mSqD, mt_in)
 !-----------------------------------------------------------------------
 ! This function calculates the strong coupling alpha_s(Q) in the DR scheme.
 ! The formula is taken from J. Bagger et al., Nucl.Phys.B
 ! The input is:
 !  - Q the energy
 !  - mG, the gluino mass
 !  - mT_in, the top mass
 !  - mSqU(i), i=1,..,6 the masses of the u-squarks
 !  - mSqD(i), i=1,..,6 the masses of the d-squarks
 ! written by Werner Porod, 19.7.1999
 ! 16.11.01: portation to f90 
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Q,mG,mSqU(6),mSqD(6)
  Real(dp), Intent(in), Optional :: mt_in

  Integer :: i
  Real(dp) :: DeltaAlpha, sumI, mt

  sumI = 0._dp
  Do i=1,6
   sumI = sumI + Log( mSqU(i) / Q ) + Log( mSqD(i) / Q )
  Enddo

  If (Present(mt_in)) Then
   mt = mt_in
  Else
   mt = mf_u(3)
  End If

  DeltaAlpha = 0.5_dp - 2._dp * Log(mG /Q) - 2._dp * Log(mt/Q ) / 3._dp &
           & - sumI / 6._dp
  DeltaAlpha = AlphaS_mZ * DeltaAlpha / ( 2._dp * Pi)
  
  AlphaSDR = AlphaS_mZ / (1._dp - DeltaAlpha)

 End Function AlphaSDR

 Subroutine Coup_DDA_1Leff(i, j, k, e, yuk_f, RP0, vevSM, mSf2, mglu, mN, mSfp2 &
    & , mC, c_CDSu_L, c_CDSu_R, c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R          &
    & , cpl_CCP0_L, cpl_CCP0_R, cpl_NNP0_L, cpl_NNP0_R, cpl_P0SdSd, cpl_P0SuSu  &
    & , coupL, coupR)
 !-----------------------------------------------------------------------
 ! -calculates the coupling between fermions and the pseudoscalars
 !  at 1-loop level, valid for 3-generation MSSM
 ! - based on the formulas of A.Buras et al., NPB 659, 3 (2003)
 ! - here it is assumed that this is an effective coupling, setting the outer
 !   momenta zero
 ! input: 
 !  i,j ........ generation index of the fermions
 !  k .......... index of pseudoscalar boson
 !  e .......... charge of the fermions
 !  Yuk_f ...... fermion yukawa couplings (3*3 matrix)
 !  RP0(i,j) ... mixing matrix of the pseudoscalar bosons
 !  vevSM(i) ... MSSM vevs (v_d, v_u)
 !  mSf2(i) .... sfermion masses squared
 !  RSF(i,j) ... sfermion mixing matrix
 !  mglu ....... mass of the mu parameter
 !  mu ......... mu parameter of the superpotential
 !  mN(i) ...... neutralino masses
 !  mSfp2(i) ... sfermion' masses squared
 !  mC(i) ...... chargino masses
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) P0(k)
 ! written by Werner Porod
 !  22.06.2003: taking  CoupFermionScalar3 as starting point
 !  10.07.03: C0 function in Buras is defined with a relative minus sign
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: e, RP0(2,2), vevSM(2), mglu, mSf2(6), mN(:)   &
      & , mC(:), mSfp2(6)
  Complex(dp), Intent(in) :: yuk_f(3,3), c_CDSu_L(:,:,:), c_CDSu_R(:,:,:) &
      & , c_DGSd_L(:,:), c_DGSd_R(:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)  &
      & , cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:), cpl_NNP0_L(:,:,:)         &
      & , cpl_NNP0_R(:,:,:), cpl_P0SdSd(:,:,:), cpl_P0SuSu(:,:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2, i3, i4, n_fp, n_char, n_neut
  Complex(dp) :: dc_L, dc_R, sigL(3,3), sigR(3,3), sigS(3,3) &
     & , DeltaM(3,3), YukR_L(3,3), YukR_R(3,3), Mass(3,3)
  Real(dp) :: C0m, mglu2, mN2(Size(mN)), C2m, mC2(Size(mC)), B0m2, B1m2
  logical :: l_QCD

  coupL = ZeroC
  coupR = ZeroC

  l_QCD = .False. ! check if gluino part is calculated
  If ((abs(e).lt.1._dp).and.(abs(e).gt.0._dp)) l_QCD = .True.
  n_char = Size(mC)
  n_neut = Size(mN)

  !------------------------------
  ! 1-loop part, 2-point function
  ! p^2=0
  !------------------------------
  sigS = ZeroC 
  sigL = ZeroC 
  sigR = ZeroC 
  !------------------------------
  ! gluino part
  !------------------------------
  if (l_QCD) then
   mglu2 = mglu**2
   Do i1=1,6
    B1m2 = -8._dp * B1(0._dp, mglu2, mSf2(i1)) / 3._dp
    B0m2 = 16._dp * mglu * B0(0._dp, mglu2, mSf2(i1)) / 3._dp
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_DGSd_R(i4,i1) ) * c_DGSd_R(i3,i1) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_DGSd_L(i4,i1) ) * c_DGSd_L(i3,i1) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_DGSd_R(i4,i1) ) * c_DGSd_L(i3,i1) * B0m2
     End Do
    End Do
   End Do
  end if
  !--------------
  ! charginos 
  !--------------
  mC2 = mC**2
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos

  Do i1=1,n_char
   Do i2=1,n_fp
    B1m2 = - 0.5_dp * B1(0._dp,mC2(i1), mSfp2(i2))
    B0m2 = mC(i1) * B0(0._dp,mC2(i1), mSfp2(i2))
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_CDSu_R(i1,i4,i2) ) * c_CDSu_R(i1,i3,i2) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_CDSu_L(i1,i4,i2) ) * c_CDSu_L(i1,i3,i2) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_CDSu_R(i1,i4,i2) ) * c_CDSu_L(i1,i3,i2) * B0m2
     End Do
    End Do
  End Do
  End Do
  !-------------
  ! neutralinos
  !-------------
  mN2 = mN**2

  Do i1=1,n_neut
   Do i2=1,6
    B0m2 = mN(i1) * B0(0._dp,mN2(i1), mSf2(i2))
    B1m2 = - 0.5_dp * B1(0._dp,mN2(i1), mSf2(i2))
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_DNSd_R(i4,i1,i2) ) * c_DNSd_R(i3,i1,i2) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_DNSd_L(i4,i1,i2) ) * c_DNSd_L(i3,i1,i2) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_DNSd_R(i4,i1,i2) ) * c_DNSd_L(i3,i1,i2) * B0m2
     End Do
    End Do
   End Do
  End Do

  !------------------------------
  ! 1-loop part, 3-point function
  !------------------------------
  dc_L = 0._dp
  dc_R = 0._dp
  !---------
  ! gluinos
  !---------
  If (l_QCD) then
   Do i1=1,6
    Do i2=1,6
     C0m = - cpl_P0SdSd(k,i1,i2) * C0_3m(mglu2, mSf2(i1), mSf2(i2)) 
     dc_L = dc_L - Conjg(c_DGSd_R(j,i2)) * c_DGSd_L(i,i1) * C0m
     dc_R = dc_R - Conjg(c_DGSd_L(j,i2)) * c_DGSd_R(i,i1) * C0m
    End Do
   End Do
   dc_L = 16._dp * mglu * dc_L / 3._dp
   dc_R = 16._dp * mglu * dc_R / 3._dp
  end if

  !--------------
  ! charginos 
  !--------------
  Do i1=1,n_char
   Do i2=1,6
    Do i3=1,6
     C0m = -mC(i1) * cpl_P0SuSu(k,i2,i3) * C0_3m(mC2(i1), mSfp2(i2), mSfp2(i3)) 
     dc_L = dc_L - Conjg(c_CDSu_R(i1,j,i3) ) * c_CDSu_L(i1,i,i2) * C0m
     dc_R = dc_R - Conjg(c_CDSu_L(i1,j,i3) ) * c_CDSu_R(i1,i,i2) * C0m
    End Do
    Do i3=1,n_char
     C0m = - mC(i1) *mC(i3) * C0_3m(mSfp2(i2), mC2(i1) , mC2(i3))
     C2m = C_2(mSfp2(i2), mC2(i1) , mC2(i3)) 
     dc_L = dc_L - Conjg(c_CDSu_R(i3,j,i2) ) * c_CDSu_L(i1,i,i2)               &
          &      * ( cpl_CCP0_R(i1, i3, k) * C2m + cpl_CCP0_L(i1, i3, k) * C0m )
     dc_R = dc_R - Conjg(c_CDSu_L(i3,j,i2) ) * c_CDSu_R(i1,i,i2)               &
          &      * ( cpl_CCP0_L(i1, i3, k) * C2m + cpl_CCP0_R(i1, i3, k) * C0m )
    End Do
   End Do
  End Do

  !--------------
  ! neutralinos 
  !--------------
  mN2 = mN**2
  Do i1=1,n_neut
   Do i2=1,6
    Do i3=1,6
     C0m = - mN(i1) * cpl_P0SdSd(k,i2,i3) * C0_3m(mN2(i1), mSf2(i2), mSf2(i3)) 
     dc_L = dc_L - Conjg(c_DNSd_R(j,i1,i3)) * c_DNSd_L(i,i1,i2) * C0m
     dc_R = dc_R - Conjg(c_DNSd_L(j,i1,i3)) * c_DNSd_R(i,i1,i2) * C0m
    End Do
    Do i3=1,n_neut
     C0m = - mN(i1)*mN(i3) * C0_3m(mSf2(i2), mN2(i1) , mN2(i3)) 
     C2m = C_2(mSf2(i2), mN2(i1) , mN2(i3)) 
     dc_L = dc_L - Conjg(c_DNSd_R(j,i3,i2)) * c_DNSd_L(i,i1,i2)               &
          &      * ( cpl_NNP0_R(i1, i3, k) * C2m + cpl_NNP0_L(i1, i3, k) * C0m )
     dc_R = dc_R - Conjg(c_DNSd_L(j,i3,i2)) * c_DNSd_R(i,i1,i2)               &
          &      * ( cpl_NNP0_L(i1, i3, k) * C2m + cpl_NNP0_R(i1, i3, k) * C0m )
    End Do
   End Do
  End Do

  coupL = oo16pi2 * dc_L
  coupR = oo16pi2 * dc_R

  DeltaM = oo16pi2 * (SigS + vevSM(1) * oosqrt2 &
         &                   * (Matmul(SigL,Yuk_f) + Matmul(Yuk_f,SigR) ) )
  Mass = oosqrt2 * yuk_f * vevSM(1) - DeltaM
  YukR_L = Matmul(SigL,Yuk_f) + Matmul(Yuk_f,SigR)
  YukR_R = Matmul(SigL,Conjg(Yuk_f)) + Matmul(Conjg(Yuk_f),SigR)
  YukR_L = oo16pi2 * RP0(k,1) * YukR_L
  YukR_R = oo16pi2 * RP0(k,1) * YukR_R

  coupL = coupL + YukR_L(i,j) &
      & - Cmplx(0._dp,RP0(k,1),dp) * DeltaM(i,j) / (vevSM(1)-sqrt2*DeltaM(i,i)/yuk_f(i,i)) 
  coupR = coupR + YukR_R(j,i) &
      & -  Cmplx(0._dp,RP0(k,1),dp) * DeltaM(j,i) / (vevSM(1)-sqrt2*DeltaM(i,i)/yuk_f(i,i)) 

 End Subroutine Coup_DDA_1Leff


 Subroutine Coup_DDH_1Leff(i, j, k, e, yuk_f, RS0, vevSM, mSf2, mglu, mN       &
    & , mSfp2, mC, c_CDSu_L, c_CDSu_R, c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R  &
    & , cpl_CCS0_L, cpl_CCS0_R, cpl_NNS0_L, cpl_NNS0_R, cpl_S0SdSd, cpl_S0SuSu &
    & , coupL, coupR)
 !-----------------------------------------------------------------------
 ! -calculates the coupling between fermions and the scalars
 !  at 1-loop level, valid for 3-generation MSSM
 ! - based on the formulas of A.Buras et al., NPB 659, 3 (2003)
 ! - here it is assumed that this is an effective coupling at low energies,
 !   setting the outer momenta zero
 ! input: 
 !  i,j ........ generation index of the fermions
 !  k .......... index of scalar boson
 !  e .......... charge of the fermions
 !  vevSM(i) ... MSSM vevs (v_d, v_u)
 !  mSf2(i) .... sfermion masses squared
 !  mglu ....... mass of the mu parameter
 !  mN(i) ...... neutralino masses
 !  N(i,j) ..... neutralino mixing matrix
 !  mSfp2(i) ... sfermion' masses squared
 !  mC(i) ...... chargino masses
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) S0(k)
 ! written by Werner Porod
 !  22.06.2003: taking  CoupFermionScalar3 as starting point
 !  10.07.03: C0 function in Buras is defined with a relative minus sign
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: e, RS0(2,2), vevSM(2), mglu, mSf2(6), mN(:)   &
      & , mC(:), mSfp2(6)
  Complex(dp), Intent(in) :: yuk_f(3,3), c_CDSu_L(:,:,:), c_CDSu_R(:,:,:) &
      & , c_DGSd_L(:,:), c_DGSd_R(:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)  &
      & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:), cpl_NNS0_L(:,:,:)         &
      & , cpl_NNS0_R(:,:,:), cpl_S0SdSd(:,:,:), cpl_S0SuSu(:,:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2, i3, i4, n_fp, n_char, n_neut
  Complex(dp) :: dc_L, dc_R, sigL(3,3), sigR(3,3), sigS(3,3) &
     & , DeltaM(3,3), YukR_L(3,3), YukR_R(3,3), Mass(3,3)
  Real(dp) :: C0m, mglu2, mN2(Size(mN)), C2m, mC2(Size(mC)), B0m2, B1m2
  logical :: l_QCD

  coupL = ZeroC
  coupR = ZeroC

  l_QCD = .False. ! check if gluino part is calculated
  If ((abs(e).lt.1._dp).and.(abs(e).gt.0._dp)) l_QCD = .True.
  n_char = Size(mC)
  n_neut = Size(mN)

  !------------------------------
  ! 1-loop part, 2-point function
  ! p^2=0
  !------------------------------
  sigS = ZeroC 
  sigL = ZeroC 
  sigR = ZeroC 
  !------------------------------
  ! gluino part
  !------------------------------
  if (l_QCD) then
   mglu2 = mglu**2
   Do i1=1,6
    B1m2 = -8._dp * B1(0._dp, mglu2, mSf2(i1)) / 3._dp
    B0m2 = 16._dp * mglu * B0(0._dp, mglu2, mSf2(i1)) / 3._dp
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_DGSd_R(i4,i1) ) * c_DGSd_R(i3,i1) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_DGSd_L(i4,i1) ) * c_DGSd_L(i3,i1) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_DGSd_R(i4,i1) ) * c_DGSd_L(i3,i1) * B0m2
     End Do
    End Do
   End Do
  end if
  !--------------
  ! charginos 
  !--------------
  mC2 = mC**2
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos

  Do i1=1,n_char
   Do i2=1,n_fp
    B1m2 = - 0.5_dp * B1(0._dp,mC2(i1), mSfp2(i2))
    B0m2 = mC(i1) * B0(0._dp,mC2(i1), mSfp2(i2))
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_CDSu_R(i1,i4,i2) ) * c_CDSu_R(i1,i3,i2) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_CDSu_L(i1,i4,i2) ) * c_CDSu_L(i1,i3,i2) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_CDSu_R(i1,i4,i2) ) * c_CDSu_L(i1,i3,i2) * B0m2
     End Do
    End Do
  End Do
  End Do
  !-------------
  ! neutralinos
  !-------------
  mN2 = mN**2

  Do i1=1,n_neut
   Do i2=1,6
    B0m2 = mN(i1) * B0(0._dp,mN2(i1), mSf2(i2))
    B1m2 = - 0.5_dp * B1(0._dp,mN2(i1), mSf2(i2))
    Do i3=1,3
     Do i4=1,3
      SigL(i3,i4) = SigL(i3,i4) &
                & + Conjg( c_DNSd_R(i4,i1,i2) ) * c_DNSd_R(i3,i1,i2) * B1m2
      SigR(i3,i4) = SigR(i3,i4) &
                & + Conjg( c_DNSd_L(i4,i1,i2) ) * c_DNSd_L(i3,i1,i2) * B1m2
      SigS(i3,i4) = SigS(i3,i4) &
                & + Conjg( c_DNSd_R(i4,i1,i2) ) * c_DNSd_L(i3,i1,i2) * B0m2
     End Do
    End Do
   End Do
  End Do

  !------------------------------
  ! 1-loop part, 3-point function
  !------------------------------
  dc_L = 0._dp
  dc_R = 0._dp
  !---------
  ! gluinos
  !---------
  If (l_QCD) then
   Do i1=1,6
    Do i2=1,6
     C0m = - cpl_S0SdSd(k,i1,i2) * C0_3m(mglu2, mSf2(i1), mSf2(i2)) 
     dc_L = dc_L - Conjg(c_DGSd_R(j,i2)) * c_DGSd_L(i,i1) * C0m
     dc_R = dc_R - Conjg(c_DGSd_L(j,i2)) * c_DGSd_R(i,i1) * C0m
    End Do
   End Do
   dc_L = 16._dp * mglu * dc_L / 3._dp
   dc_R = 16._dp * mglu * dc_R / 3._dp
  end if
  !--------------
  ! charginos 
  !--------------
  Do i1=1,n_char
   Do i2=1,6
    Do i3=1,6
     C0m = -mC(i1) * cpl_S0SuSu(k,i2,i3) * C0_3m(mC2(i1), mSfp2(i2), mSfp2(i3))
     dc_L = dc_L - Conjg(c_CDSu_R(i1,j,i3) ) * c_CDSu_L(i1,i,i2) * C0m
     dc_R = dc_R - Conjg(c_CDSu_L(i1,j,i3) ) * c_CDSu_R(i1,i,i2) * C0m
    End Do
    Do i3=1,n_char
     C0m = - mC(i1)*mC(i3) * C0_3m(mSfp2(i2), mC2(i1) , mC2(i3))
     C2m = C_2(mSfp2(i2), mC2(i1) , mC2(i3)) 
     dc_L = dc_L - Conjg(c_CDSu_R(i3,j,i2) ) * c_CDSu_L(i1,i,i2)               &
          &      * ( cpl_CCS0_R(i1, i3, k) * C2m + cpl_CCS0_L(i1, i3, k) * C0m )
     dc_R = dc_R - Conjg(c_CDSu_L(i3,j,i2) ) * c_CDSu_R(i1,i,i2)               &
          &      * ( cpl_CCS0_L(i1, i3, k) * C2m + cpl_CCS0_R(i1, i3, k) * C0m )
    End Do
   End Do
  End Do

  !--------------
  ! neutralinos 
  !--------------
  mN2 = mN**2
  Do i1=1,n_neut
   Do i2=1,6
    Do i3=1,6
     C0m = - mN(i1) * cpl_S0SdSd(k,i2,i3) * C0_3m(mN2(i1), mSf2(i2), mSf2(i3)) 
     dc_L = dc_L - Conjg(c_DNSd_R(j,i1,i3)) * c_DNSd_L(i,i1,i2) * C0m
     dc_R = dc_R - Conjg(c_DNSd_L(j,i1,i3)) * c_DNSd_R(i,i1,i2) * C0m
    End Do
    Do i3=1,n_neut
     C0m = - mN(i1)*mN(i3) * C0_3m(mSf2(i2), mN2(i1) , mN2(i3)) 
     C2m = C_2(mSf2(i2), mN2(i1) , mN2(i3)) 
     dc_L = dc_L - Conjg(c_DNSd_R(j,i3,i2)) * c_DNSd_L(i,i1,i2)               &
          &      * ( cpl_NNS0_R(i1, i3, k) * C2m + cpl_NNS0_L(i1, i3, k) * C0m )
     dc_R = dc_R - Conjg(c_DNSd_L(j,i3,i2)) * c_DNSd_R(i,i1,i2)               &
          &      * ( cpl_NNS0_L(i1, i3, k) * C2m + cpl_NNS0_R(i1, i3, k) * C0m )
    End Do
   End Do
  End Do

  coupL = oo16pi2 * dc_L
  coupR = oo16pi2 * dc_R

  DeltaM = oo16pi2 * (SigS + vevSM(1) * oosqrt2 &
         &                   * (Matmul(SigL,Yuk_f) + Matmul(Yuk_f,SigR) ) )
  Mass = oosqrt2 * yuk_f * vevSM(1) - DeltaM
  YukR_L = Matmul(SigL,Yuk_f) + Matmul(Yuk_f,SigR)
  YukR_R = Matmul(SigL,Conjg(Yuk_f)) + Matmul(Conjg(Yuk_f),SigR)
  YukR_L = oo16pi2 * RS0(k,1) * YukR_L
  YukR_R = oo16pi2 * RS0(k,1) * YukR_R

  coupL = coupL + YukR_L(i,j) &
      & + RS0(k,1) * DeltaM(i,j) / (vevSM(1)-sqrt2*DeltaM(i,i)/yuk_f(i,i)) 
  coupR = coupR + YukR_R(j,i) &
      & + RS0(k,1) * DeltaM(j,i) / (vevSM(1)-sqrt2*DeltaM(i,i)/yuk_f(i,i)) 

 End Subroutine Coup_DDH_1Leff


 Subroutine CoupFermionPseudoScalar31L_eff(i, j, k, T3, e, g, yuk_f, Rfl, Rfr &
    & , RP0, vevSM, mSf2, RSf, A_f, phi_glu, mglu, mu, mN, N, mSfp2, RSfp     &
    & , Yuk_fp, A_fp, mC, U, V, coupL, coupR)
 !-----------------------------------------------------------------------
 ! -calculates the coupling between fermions and the pseudoscalars
 !  at 1-loop level, valid for 3-generation MSSM
 ! - based on the formulas of A.Buras et al., NPB 659, 3 (2003)
 ! - here it is assumed that this is an effective coupling, setting the outer
 !   momenta zero
 ! input: 
 !  i,j ........ generation index of the fermions
 !  k .......... index of pseudoscalar boson
 !  T3 ......... isospin of the fermions
 !  e .......... charge of the fermions
 !  g(i) ....... gauge couplings [ U(1), SU(2), SU(3) ]
 !  Yuk_f ...... fermion yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RP0(i,j) ... mixing matrix of the pseudoscalar bosons
 !  vevSM(i) ... MSSM vevs (v_d, v_u)
 !  mSf2(i) .... sfermion masses squared
 !  RSF(i,j) ... sfermion mixing matrix
 !  A_f(i,j) ... trilinear sfermion Higgs couplings
 !  phi_glu .... phase of M_3
 !  mglu ....... mass of the mu parameter
 !  mu ......... mu parameter of the superpotential
 !  mN(i) ...... neutralino masses
 !  N(i,j) ..... neutralino mixing matrix
 !  mSfp2(i) ... sfermion' masses squared
 !  RSFp(i,j) .. sfermion' mixing matrix
 !  Yuk_fp ..... fermion' yukawa couplings (3*3 matrix)
 !  A_fp(i,j) .. trilinear sfermion' Higgs couplings
 !  mC(i) ...... chargino masses
 !  U(i,j) ..... chargino mixing matrix
 !  V(i,j) ..... chargino mixing matrix
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) P0(k)
 ! written by Werner Porod
 !  22.06.2003: taking  CoupFermionScalar3 as starting point
 !  10.07.03: C0 function in Buras is defined with a relative minus sign
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3, e, RP0(2,2), vevSM(2), g(3), mglu, mSf2(6)   &
      & , mN(4), mC(2), mSfp2(6)
  Complex(dp), Intent(in) :: yuk_f(3,3), RFl(3,3), RFr(3,3), phi_glu, mu   &
      & , RSf(6,6), A_f(3,3), N(4,4), RSfp(6,6), U(2,2), V(2,2), A_fp(3,3) &
      & , Yuk_fp(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2, i3, n_fp, n_char, n_neut
  Complex(dp) :: dc_L, dc_R, coupC, coupRC1, coupLC1, coupRC2, coupLC2     &
      & , coupRC3, coupLC3, sigLR, sigRL, bi(1)
  Real(dp) :: C0m, mglu2, mN2(4), C2m, mC2(2), e_u, B0m2
  logical :: l_QCD

  coupL = ZeroC
  coupR = ZeroC

  l_QCD = .False. ! check if gluino part is calculated
  If ((abs(e).lt.1._dp).and.(abs(e).gt.0._dp)) l_QCD = .True.
  n_char = Size(mC)
  n_neut = Size(mN)
  bi = mu
  !--------------------------------
  ! tree level part
  !--------------------------------
  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk_f(j,i)
     coupR = Conjg(yuk_f(i,j)) 

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   Do i2=1,3
    coupL = coupL + yuk_f(j,i2) * Conjg(RFr(i,i2))
    coupR = coupR + Conjg(yuk_f(i,i2)) * RFr(j,i2)    
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk_f(i1,i)
    coupR = coupR + RFl(i,i1) * Conjg(yuk_f(i1,j))
   End Do

  Else
   Do i1=1,3
    Do i2=1,3
     coupL = coupL + Conjg(RFl(j,i1)) * yuk_f(i1,i2) * Conjg(RFr(i,i2))
     coupR = coupR + RFl(i,i1) * Conjg(yuk_f(i1,i2)) * RFr(j,i2)
    End Do
   End Do
  End If

  If (T3.Gt.0._dp) Then
   coupL = coupL * RP0(k,2) * Cmplx(0._dp,-oosqrt2,dp)
   coupR = coupR * RP0(k,2) * Cmplx(0._dp,oosqrt2,dp)
  Else
   coupL = coupL * RP0(k,1) * Cmplx(0._dp,-oosqrt2,dp)
   coupR = coupR * RP0(k,1) * Cmplx(0._dp,oosqrt2,dp)
  End If

  !------------------------------
  ! 1-loop part, 2-point function
  ! p^2=0
  !------------------------------
  sigLR = 0._dp
  sigRL = 0._dp
  !------------------------------
  ! gluino part
  !------------------------------
  if (l_QCD) then
   mglu2 = mglu**2
   Do i1=1,6
    B0m2 = 16._dp * mglu * B0(0._dp, mglu2, mSf2(i1)) / 3._dp
    Call CoupGluinoSquark(g(3), phi_glu, i, i1, RSf, RfL, RfR, coupLC1, coupRC1)
    Call CoupGluinoSquark(g(3), phi_glu, j, i1, RSf, RfL, RfR, coupLC2, coupRC2)
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  end if
  !--------------
  ! charginos 
  !--------------
  mC2 = mC**2
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos

  Do i1=1,n_char
   Do i2=1,n_fp
    B0m2 = mC(i1) * B0(0._dp,mC2(i1), mSfp2(i2))
    If (T3.Gt.0._dp) Then
     Call CoupCharginoSfermion(i1, i, i2, g(2), T3, RSfp, Yuk_fp, Yuk_f &
                                &, RfL, RfR, U, V, coupLC1, coupRC1)
     Call CoupCharginoSfermion(i1, j, i2, g(2), T3, RSfp, Yuk_fp, Yuk_f &
                                &, RfL, RfR, U, V, coupLC2, coupRC2)
    Else
     Call CoupCharginoSfermion(i1, i, i2, g(2), T3, RSfp, Yuk_f, Yuk_fp &
                                &, RfL, RfR, U, V, coupLC1, coupRC1)
     Call CoupCharginoSfermion(i1, j, i2, g(2), T3, RSfp, Yuk_f, Yuk_fp &
                                &, RfL, RfR, U, V, coupLC2, coupRC2)
    End If
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  End Do
  !-------------
  ! neutralinos
  !-------------
  mN2 = mN**2
  Do i1=1,n_neut
   Do i2=1,6
    B0m2 = mN(i1) * B0(0._dp,mN2(i1), mSf2(i2))
    Call CoupNeutralinoSfermion(i, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC1, coupRC1)
    Call CoupNeutralinoSfermion(j, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  End Do

  SigRL = oosqrt2 * yuk_f(i,i) * oo16pi2 * SigRL
  SigLR = oosqrt2 * Conjg(yuk_f(j,j)) * oo16pi2 * SigLR

  !------------------------------
  ! 1-loop part, 3-point function
  !------------------------------
  dc_L = 0._dp
  dc_R = 0._dp
  !---------
  ! gluinos
  !---------
  If (l_QCD) then
   Do i1=1,6
    Call CoupGluinoSquark3(g(3), phi_glu, i, i1, RSf, RfL, RfR, coupLC1, coupRC1)
    Do i2=1,6
     Call CoupGluinoSquark3(g(3), phi_glu, j, i2, RSf, RfL,RfR, coupLC2, coupRC2)
     Call CoupPseudoScalarSfermion3(k, i1, i2, RP0, T3, yuk_f, RSf, A_f, bi &
                                & , coupC)
     C0m = - C0_3m(mglu2, mSf2(i1), mSf2(i2)) 
     dc_L = dc_L - coupC * Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - coupC * Conjg(coupLC2) * coupRC1 * C0m
    End Do
   End Do
   dc_L = 16._dp * mglu * dc_L / 3._dp
   dc_R = 16._dp * mglu * dc_R / 3._dp
  end if

  !--------------
  ! charginos 
  !--------------
  If (e.Eq.-1._dp) Then
   e_u = 0._dp
  Else If (e.Eq.-1._dp/3._dp) Then
   e_u = 2._dp / 3._dp
  Else If (e.Eq.2._dp/3._dp) Then
   e_u = -1._dp / 3._dp
  End If

  Do i1=1,2
   Do i2=1,6
    If (T3.gt.0._dp) then
     Call CoupCharginoSfermion3(i1, i, i2, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                               &, RfR, U, V, coupLC1, coupRC1)
    else
     Call CoupCharginoSfermion3(i1, i, i2, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                               &, RfR, U, V, coupLC1, coupRC1)
    end if
    Do i3=1,6
     If (T3.gt.0._dp) then
      Call CoupCharginoSfermion3(i1, j, i3, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                               &, RfR, U, V, coupLC2, coupRC2)
     else
      Call CoupCharginoSfermion3(i1, j, i3, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                               &, RfR, U, V, coupLC2, coupRC2)
     end if
     Call CoupPseudoScalarSfermion3(k, i2, i3, RP0, -T3, yuk_fp, RSf, A_fp, bi &
         & , coupC)
     C0m = - mC(i1) * coupC * C0_3m(mC2(i1), mSfp2(i2), mSfp2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * C0m
    End Do
    Do i3=1,2
     If (T3.gt.0._dp) then
      Call CoupCharginoSfermion3(i3, j, i2, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                                &, RfR, U, V, coupLC2, coupRC2)
     else
      Call CoupCharginoSfermion3(i3, j, i2, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                                &, RfR, U, V, coupLC2, coupRC2)
     end if
     Call CoupCharginoPseudoScalar(i1, i3, k, U, V , RP0, g(2), coupLC3, coupRC3)
     C0m = - mC(i1) *mC(i3) * C0_3m(mSfp2(i2), mC2(i1) , mC2(i3))
     C2m = C_2(mSfp2(i2), mC2(i1) , mC2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * ( coupRC3 * C2m + coupLC3 * C0m )
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * ( coupLC3 * C2m + coupRC3 * C0m )
    End Do
   End Do
  End Do

  !--------------
  ! neutralinos 
  !--------------
  mN2 = mN**2
  Do i1=1,4
   Do i2=1,6
    Call CoupNeutralinoSfermion3(i, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC1, coupRC1)
    Do i3=1,6
     Call CoupNeutralinoSfermion3(j, i1, i3, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
     Call CoupPseudoScalarSfermion3(k, i2, i3, RP0, T3, yuk_f, RSf, A_f, bi &
                               & , coupC)
     C0m = - mN(i1) * coupC * C0_3m(mN2(i1), mSf2(i2), mSf2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * C0m
    End Do
    Do i3=1,4
     Call CoupNeutralinoSfermion3(j, i3, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
     Call CoupNeutralinoPseudoScalar(i1, i3, k, N, RP0, g(1), g(2) &
                              & , coupLC3, coupRC3)
     C0m = - mN(i1)*mN(i3) * C0_3m(mSf2(i2), mN2(i1) , mN2(i3)) 
     C2m = C_2(mSf2(i2), mN2(i1) , mN2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * ( coupRC3 * C2m + coupLC3 * C0m )
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * ( coupLC3 * C2m + coupRC3 * C0m )
    End Do
   End Do
  End Do

  If (T3.gt.0._dp) then
   coupL = coupL - SigRL * RP0(k,2) / vevSM(2)
   coupR = coupR - SigLR * RP0(k,2) / vevSM(2)
  else
   coupL = coupL - SigRL * RP0(k,1) / vevSM(1)
   coupR = coupR - SigLR * RP0(k,1) / vevSM(1)
  end if
  coupL = coupL + oo16pi2 * dc_L
  coupR = coupR + oo16pi2 * dc_R

 End Subroutine CoupFermionPseudoScalar31L_eff


 Subroutine CoupFermionScalar31L_eff(i, j, k, T3, e, g, yuk_f, Rfl, Rfr, RS0 &
    & , vevSM, mSf2, RSf, A_f, phi_glu, mglu, mu, mN, N, mSfp2, RSfp, Yuk_fp &
    & , A_fp, mC, U, V, coupL, coupR)
 !-----------------------------------------------------------------------
 ! -calculates the coupling between fermions and the scalars
 !  at 1-loop level, valid for 3-generation MSSM
 ! - based on the formulas of A.Buras et al., NPB 659, 3 (2003)
 ! - here it is assumed that this is an effective coupling, setting the outer
 !   momenta zero
 ! input: 
 !  i,j ........ generation index of the fermions
 !  k .......... index of scalar boson
 !  T3 ......... isospin of the fermions
 !  e .......... charge of the fermions
 !  g(i) ....... gauge couplings [ U(1), SU(2), SU(3) ]
 !  Yuk_f ...... fermion yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  vevSM(i) ... MSSM vevs (v_d, v_u)
 !  mSf2(i) .... sfermion masses squared
 !  RSF(i,j) ... sfermion mixing matrix
 !  A_f(i,j) ... trilinear sfermion Higgs couplings
 !  phi_glu .... phase of M_3
 !  mglu ....... mass of the mu parameter
 !  mu ......... mu parameter of the superpotential
 !  mN(i) ...... neutralino masses
 !  N(i,j) ..... neutralino mixing matrix
 !  mSfp2(i) ... sfermion' masses squared
 !  RSFp(i,j) .. sfermion' mixing matrix
 !  Yuk_fp ..... fermion' yukawa couplings (3*3 matrix)
 !  A_fp(i,j) .. trilinear sfermion' Higgs couplings
 !  mC(i) ...... chargino masses
 !  U(i,j) ..... chargino mixing matrix
 !  V(i,j) ..... chargino mixing matrix
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) S0(k)
 ! written by Werner Porod
 !  22.06.2003: taking  CoupFermionScalar3 as starting point
 !  10.07.03: C0 function in Buras is defined with a relative minus sign
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: T3, e, RS0(2,2), vevSM(2), g(3), mglu, mSf2(6)   &
      & , mN(4), mC(2), mSfp2(6)
  Complex(dp), Intent(in) :: yuk_f(3,3), RFl(3,3), RFr(3,3), phi_glu, mu   &
      & , RSf(6,6), A_f(3,3), N(4,4), RSfp(6,6), U(2,2), V(2,2), A_fp(3,3) &
      & , Yuk_fp(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2, i3, n_fp, n_char, n_neut
  Complex(dp) :: dc_L, dc_R, coupC, coupRC1, coupLC1, coupRC2, coupLC2     &
      & , coupRC3, coupLC3, sigLR, sigRL
  Real(dp) :: C0m, mglu2, mN2(4), C2m, mC2(2), e_u, B0m2
  logical :: l_QCD

  coupL = ZeroC
  coupR = ZeroC

  l_QCD = .False. ! check if gluino part is calculated
  If ((abs(e).lt.1._dp).and.(abs(e).gt.0._dp)) l_QCD = .True.
  n_char = Size(mC)
  n_neut = Size(mN)

  !--------------------------------
  ! tree level part
  !--------------------------------
  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk_f(j,i)
     coupR = Conjg(yuk_f(i,j)) 

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   Do i2=1,3
    coupL = coupL + yuk_f(j,i2) * Conjg(RFr(i,i2))
    coupR = coupR + Conjg(yuk_f(i,i2)) * RFr(j,i2)    
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk_f(i1,i)
    coupR = coupR + RFl(i,i1) * Conjg(yuk_f(i1,j))
   End Do

  Else
   Do i1=1,3
    Do i2=1,3
     coupL = coupL + Conjg(RFl(j,i1)) * yuk_f(i1,i2) * Conjg(RFr(i,i2))
     coupR = coupR + RFl(i,i1) * Conjg(yuk_f(i1,i2)) * RFr(j,i2)
    End Do
   End Do
  End If

  If (T3.Gt.0._dp) Then
   coupL = - coupL * RS0(k,2) * oosqrt2
   coupR = - coupR * RS0(k,2) * oosqrt2
  Else
   coupL = - coupL * RS0(k,1) * oosqrt2
   coupR = - coupR * RS0(k,1) * oosqrt2
  End If

  !------------------------------
  ! 1-loop part, 2-point function
  ! p^2=0
  !------------------------------
  sigLR = 0._dp
  sigRL = 0._dp
  !------------------------------
  ! gluino part
  !------------------------------
  if (l_QCD) then
   mglu2 = mglu**2
   Do i1=1,6
    B0m2 = 16._dp * mglu * B0(0._dp, mglu2, mSf2(i1)) / 3._dp
    Call CoupGluinoSquark(g(3), phi_glu, i, i1, RSf, RfL, RfR, coupLC1, coupRC1)
    Call CoupGluinoSquark(g(3), phi_glu, j, i1, RSf, RfL, RfR, coupLC2, coupRC2)
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  end if
  !--------------
  ! charginos 
  !--------------
  mC2 = mC**2
  n_fp = 6 ! number of sfermions'
  If (e.Eq.-1._dp) n_fp = 3  ! sneutrinos

  Do i1=1,n_char
   Do i2=1,n_fp
    B0m2 = mC(i1) * B0(0._dp,mC2(i1), mSfp2(i2))
    If (T3.Gt.0._dp) Then
     Call CoupCharginoSfermion(i1, i, i2, g(2), T3, RSfp, Yuk_fp, Yuk_f &
                                &, RfL, RfR, U, V, coupLC1, coupRC1)
     Call CoupCharginoSfermion(i1, j, i2, g(2), T3, RSfp, Yuk_fp, Yuk_f &
                                &, RfL, RfR, U, V, coupLC2, coupRC2)
    Else
     Call CoupCharginoSfermion(i1, i, i2, g(2), T3, RSfp, Yuk_f, Yuk_fp &
                                &, RfL, RfR, U, V, coupLC1, coupRC1)
     Call CoupCharginoSfermion(i1, j, i2, g(2), T3, RSfp, Yuk_f, Yuk_fp &
                                &, RfL, RfR, U, V, coupLC2, coupRC2)
    End If
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  End Do
  !-------------
  ! neutralinos
  !-------------
  mN2 = mN**2
  Do i1=1,n_neut
   Do i2=1,6
    B0m2 = mN(i1) * B0(0._dp,mN2(i1), mSf2(i2))
    Call CoupNeutralinoSfermion(i, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC1, coupRC1)
    Call CoupNeutralinoSfermion(j, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
    sigLR = sigLR + Conjg( coupLC2 ) * coupRC1 * B0m2
    sigRL = sigRL + Conjg( coupRC2 ) * coupLC1 * B0m2
   End Do
  End Do

  SigRL = Cmplx(0._dp,oosqrt2,dp) * yuk_f(i,i) * oo16pi2 * SigRL
  SigLR = Cmplx(0._dp,-oosqrt2,dp) * Conjg(yuk_f(j,j)) * oo16pi2 * SigLR

  !------------------------------
  ! 1-loop part, 3-point function
  !------------------------------
  dc_L = 0._dp
  dc_R = 0._dp
  !---------
  ! gluinos
  !---------
  If (l_QCD) then
   Do i1=1,6
    Call CoupGluinoSquark3(g(3), phi_glu, i, i1, RSf, RfL, RfR, coupLC1, coupRC1)
    Do i2=1,6
     Call CoupGluinoSquark3(g(3), phi_glu, j, i2, RSf, RfL,RfR, coupLC2, coupRC2)
     Call CoupScalarSfermion3(k, i1, i2, RS0, T3, e, yuk_f, RSf, A_f, mu, &
                       &          vevSM, g(1), g(2), coupC)
     C0m = - C0_3m(mglu2, mSf2(i1), mSf2(i2)) 
     dc_L = dc_L - coupC * Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - coupC * Conjg(coupLC2) * coupRC1 * C0m
    End Do
   End Do
   dc_L = 16._dp * mglu * dc_L / 3._dp
   dc_R = 16._dp * mglu * dc_R / 3._dp
  end if
  !--------------
  ! charginos 
  !--------------
  If (e.Eq.-1._dp) Then
   e_u = 0._dp
  Else If (e.Eq.-1._dp/3._dp) Then
   e_u = 2._dp / 3._dp
  Else If (e.Eq.2._dp/3._dp) Then
   e_u = -1._dp / 3._dp
  End If

  Do i1=1,2
   Do i2=1,6
    If (T3.gt.0._dp) then
     Call CoupCharginoSfermion3(i1, i, i2, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                               &, RfR, U, V, coupLC1, coupRC1)
    else
     Call CoupCharginoSfermion3(i1, i, i2, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                               &, RfR, U, V, coupLC1, coupRC1)
    end if
    Do i3=1,6
     If (T3.gt.0._dp) then
      Call CoupCharginoSfermion3(i1, j, i3, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                               &, RfR, U, V, coupLC2, coupRC2)
     else
      Call CoupCharginoSfermion3(i1, j, i3, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                               &, RfR, U, V, coupLC2, coupRC2)
     end if
     Call CoupScalarSfermion3(k, i2, i3, RS0, -T3, e_u, yuk_fp, RSf, A_fp, mu, &
                      &          vevSM, g(1), g(2), coupC)
     C0m = - mC(i1) * coupC * C0_3m(mC2(i1), mSfp2(i2), mSfp2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * C0m
    End Do
    Do i3=1,2
     If (T3.gt.0._dp) then
      Call CoupCharginoSfermion3(i3, j, i2, g(2), T3, RSf, yuk_fp, yuk_f, RfL &
                                &, RfR, U, V, coupLC2, coupRC2)
     else
      Call CoupCharginoSfermion3(i3, j, i2, g(2), T3, RSf, yuk_f, yuk_fp, RfL &
                                &, RfR, U, V, coupLC2, coupRC2)
     end if
     Call CoupCharginoScalar(i1, i3, k, U, V , RS0, g(2), coupLC3, coupRC3)
     C0m = - mC(i1)*mC(i3) * C0_3m(mSfp2(i2), mC2(i1) , mC2(i3))
     C2m = C_2(mSfp2(i2), mC2(i1) , mC2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * ( coupRC3 * C2m + coupLC3 * C0m )
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * ( coupLC3 * C2m + coupRC3 * C0m )
    End Do
   End Do
  End Do

  !--------------
  ! neutralinos 
  !--------------
  mN2 = mN**2
  Do i1=1,4
   Do i2=1,6
    Call CoupNeutralinoSfermion3(i, i1, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC1, coupRC1)
    Do i3=1,6
     Call CoupNeutralinoSfermion3(j, i1, i3, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
     Call CoupScalarSfermion3(k, i2, i3, RS0, T3, e, yuk_f, RSf, A_f, mu, &
                      &          vevSM, g(1), g(2), coupC)
     C0m = - mN(i1) * coupC * C0_3m(mN2(i1), mSf2(i2), mSf2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * C0m
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * C0m
    End Do
    Do i3=1,4
     Call CoupNeutralinoSfermion3(j, i3, i2, g(1), g(2), T3, e, RSf, RfL &
                               &, RfR, Yuk_f, N, coupLC2, coupRC2)
     Call CoupNeutralinoScalar(i1, i3, k, N, RS0, g(1), g(2), coupLC3, coupRC3)
     C0m = - mN(i1)*mN(i3) * C0_3m(mSf2(i2), mN2(i1) , mN2(i3)) 
     C2m = C_2(mSf2(i2), mN2(i1) , mN2(i3)) 
     dc_L = dc_L - Conjg(coupRC2) * coupLC1 * ( coupRC3 * C2m + coupLC3 * C0m )
     dc_R = dc_R - Conjg(coupLC2) * coupRC1 * ( coupLC3 * C2m + coupRC3 * C0m )
    End Do
   End Do
  End Do

  If (T3.gt.0._dp) then
   coupL = coupL - SigRL * RS0(k,2) / vevSM(2)
   coupR = coupR - SigLR * RS0(k,2) / vevSM(2)
  else
   coupL = coupL - SigRL * RS0(k,1) / vevSM(1)
   coupR = coupR - SigLR * RS0(k,1) / vevSM(1)
  end if
  coupL = coupL + oo16pi2 * dc_L
  coupR = coupR + oo16pi2 * dc_R

 End Subroutine CoupFermionScalar31L_eff


 Subroutine CoupPseudoScalarGluon(m_H2, mf_u2, g_u, mf_d2, g_d, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a pseudoscalar
 ! to two gluons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! based on  CoupScalarGluon by Werner Porod
 ! written by Florian Staub, 18.06.2010
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mf_u2(3), mf_d2(3)
  Complex(dp), Intent(in) :: g_u(3), g_d(3)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + AP_onehalf(mH2p/ mf_d2(i1))+ AP_onehalf(mH2p/ mf_u2(i1))
   End Do
   coupSM = 0.75_dp * coupSM
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_d(i1) * AP_onehalf(mH2p/ mf_d2(i1))  &
         &      + g_u(i1) * AP_onehalf(mH2p/ mf_u2(i1))
   End Do
  
  coup = 0.75_dp * coup
 
 End Subroutine CoupPseudoScalarGluon


 Subroutine CoupPseudoScalarPhoton(m_H2,mf_u2, g_u, mf_d2, g_d, mf_l2    &
    & , g_l, mC2, g_C, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a pseudoscalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! based on  CoupScalarPhoton by Werner Porod
 ! written by Florian Staub, 18.06.2010
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2,  mf_u2(3), mf_d2(3), mf_l2(3), mC2(2)
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_l(3), g_C(2)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + AP_onehalf(mH2p/ mf_l2(i1)) + (AP_onehalf(mH2p/ mf_d2(i1)) &
                                + 4._dp * AP_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_l(i1) * AP_onehalf(mH2p/ mf_l2(i1))    &
         &      + ( g_d(i1) * AP_onehalf(mH2p/ mf_d2(i1))  &
         &        + 4._dp * g_u(i1) * AP_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  !--------------------------------------
  ! chargino contributions
  !--------------------------------------
  Do i1=1,2
   coup = coup + g_C(i1) * AP_onehalf(mH2p/ mC2(i1))
  End Do
  
  
 End Subroutine CoupPseudoScalarPhoton


 Subroutine CoupScalarGluon(m_H2, mf_u2, g_u, mf_d2, g_d, mSup2, g_su &
                          & , mSdown2, g_sd, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a scalar
 ! to two gluons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 3.1.2009
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mf_u2(3), mf_d2(3), mSup2(6), mSdown2(6)
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_su(6), g_sd(6)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + A_onehalf(mH2p/ mf_d2(i1))+ A_onehalf(mH2p/ mf_u2(i1))
   End Do
   coupSM = 0.75_dp * coupSM
  End If

  coup = 0._dp
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_d(i1) * A_onehalf(mH2p/ mf_d2(i1))  &
         &      + g_u(i1) * A_onehalf(mH2p/ mf_u2(i1))
   End Do
  !--------------------------------------
  ! Squark contributions
  !--------------------------------------
  Do i1=1,6
   coup = coup + g_sd(i1) * A_zero(mH2p/ mSdown2(i1)) &
        &      + g_su(i1) * A_zero(mH2p/ mSup2(i1))
  End Do
  
  coup = 0.75_dp * coup
 
 End Subroutine CoupScalarGluon

 Subroutine CoupScalarPhoton(m_H2, mW2, g_W, mf_u2, g_u, r_T, mf_d2, g_d, mf_l2 &
    & , g_l, mC2, g_C, m_Hp2, g_Hp, mSup2, g_su, mSdown2, g_sd, r_sq            &
    & , mSlept2, g_sl, coup, coupSM)
 !----------------------------------------------------------------------------
 ! calculates the lowest order coupling for the radiative decay of a scalar
 ! to two photons
 ! based on the formulas of M.Spira et al., Nucl. Phys. B 453 (1995) 17
 ! written by Werner Porod, 27.12.2008
 ! 8.6.2012: adding approximate QCD corrections based on above paper,
 !           using form factors which are however only useful for light Higgs
 !           which should be set to 1 for heavier Higgs bosons
 !----------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_H2, mW2, mf_u2(3), mf_d2(3), mf_l2(3), m_Hp2 &
    & , mSup2(6), mSdown2(6), mSlept2(6), mC2(2)
  Real(dp), Intent(in) :: g_W, g_Hp, r_T, r_sq 
  Complex(dp), Intent(in) :: g_u(3), g_d(3), g_l(3), g_C(2), g_su(6), g_sd(6) &
    & , g_sl(6)
  Complex(dp), Intent(out) :: coup
  Complex(dp), Optional, Intent(out) :: coupSM
 
  Integer :: i1
  Real(dp) :: Mh2p
  !--------------------------------------
  ! W-boson contribution
  !--------------------------------------
  mH2p = 0.25_dp * m_H2

  If (Present(coupSM)) Then
   coupSM = A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coupSM = coupSM + A_onehalf(mH2p/ mf_l2(i1)) &
         & + ( A_onehalf(mH2p/ mf_d2(i1))  &
         &   + r_T * 4._dp * A_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  End If

  coup = g_W * A_one(mH2p / mW2)
  !--------------------------------------
  ! SM-Fermion contributions
  !--------------------------------------
   Do i1=1,3
    coup = coup + g_l(i1) * A_onehalf(mH2p/ mf_l2(i1))   &
        &      + ( g_d(i1) * A_onehalf(mH2p/ mf_d2(i1))  &
        &        + r_T * 4._dp * g_u(i1) * A_onehalf(mH2p/ mf_u2(i1)) ) / 3._dp
   End Do
  !--------------------------------------
  ! chargino contributions
  !--------------------------------------
  Do i1=1,2
   coup = coup + g_C(i1) * A_onehalf(mH2p/ mC2(i1))
  End Do
  !--------------------------------------
  ! charged Higgs contributions
  !--------------------------------------
  coup = coup + g_Hp * A_zero(mH2p/ m_Hp2)
  !--------------------------------------
  ! Sfermion contributions
  !--------------------------------------
  Do i1=1,6
   coup = coup + g_sl(i1) * A_zero(mH2p/ mSlept2(i1))          &
        &      + r_sq * ( g_sd(i1) * A_zero(mH2p/ mSdown2(i1)) &
        &               + 4._dp * g_su(i1) * A_zero(mH2p/ mSup2(i1)) ) / 3._dp
  End Do
  
 End Subroutine CoupScalarPhoton

 Subroutine InitializeLoopCouplings(vevs)
 implicit none
  Real(dp), intent(in) :: vevs(2)

  vevSM = vevs
 end Subroutine InitializeLoopCouplings

 Subroutine RunningCouplings(Q,g,yuk)
 !-----------------------------------------------------------------------
 ! calculates running gauge couplings and the Yukawa couplings of the
 ! third generation to two-loop order
 ! written by Werner Porod, 6.10.2000
 ! - correct change to DRbar scheme is missing, using On-shell and MSbar
 !   values as input
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Q
  Real(dp), Intent(out) :: g(3), yuk(3)

  Integer :: kont
  Real(dp) :: g1(6), tz, dt, cosW2, sinW2

  g1 = 0._dp
  cosW2 = mW2 / mZ2
  sinW2 = 1._dp - cosW2
  !----------------------------------
  ! running from Z-scale to scale Q
  !----------------------------------
  g1(1) = Sqrt( 20._dp*pi*alpha_mZ/(3._dp*cosW2) )
  g1(2) = Sqrt( 4._dp*pi*alpha_mZ/sinW2)
  g1(3) = Sqrt( 4._dp * pi * alphaS_mZ )

  g1(4) = sqrt2 * mf_l(3) / vevSM(1)
  g1(5) = sqrt2 * mf_d(3) / vevSM(1)
  g1(6) = sqrt2 * mf_u(3) / vevSM(2)

  tz = Log(mZ/Q)
  dt = - tz / 50._dp

  Call odeint(g1, 6, tz, 0._dp, 1.e-5_dp, dt, 0._dp, rge6, kont)

  g1(1) = Sqrt(3._dp/5._dp) * g1(1)

  g = g1(1:3)
  yuk = g1(4:6)

 End Subroutine RunningCouplings


 Real(Dp) Function SUALFE(QS)
 !----------------------------------------------------------------------
 !     Returns the running EM coupling alpha_em(q**2)
 !     taken from ISASUSY
 !  is used to calculate the coupling at low energies 
 !
 ! SEE BARGER/PHILLIPS, P. 202 
 !----------------------------------------------------------------------
 Implicit None

  Real(Dp), Intent(in) :: qs
  Integer :: i1
  Real(Dp) :: sum, qu2, qd2, qs4

  SUM=0._dp
  QU2=4._dp/3._dp
  QD2=1._dp/3._dp
  qs4 = 0.25_dp * QS 
  Do i1=1,3
   If (qs4.Gt.mf_l2(i1)) SUM=SUM+Log(qs4/mf_l2(i1) )
   If (qs4.Gt.mf_u2(i1)) SUM=SUM+QU2*Log(qs4/mf_u2(i1) )
   If (qs4.Gt.mf_d2(i1)) SUM=SUM+QD2*Log(qs4/mf_d2(i1) )
  End Do

  SUALFE = Alpha/(1._dp- Alpha * sum /(3._dp*PI) )

 End Function SUALFE

End Module LoopCouplings
