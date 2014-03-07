Module LHC_observables
!----------------------------------------------------------
! contains functions for kinematical observables at LHC
!----------------------------------------------------------
Use Control
Use StandardModel
Use Couplings
Use SusyDecays

Contains


 Subroutine Calc_LHC_observables(mN, N, mC, U, V, mSlepton, Rslepton        &
     & , mSdown, RSdown, mSup, RSup, mGlu, PhaseGlu, mS0, RS0, mP0, RP0     &
     & , mSpm, RSpm, gauge, Y_u, Y_d, A_u, A_d, mu, vevSM, GenerationMixing &
     & , observ)
 !--------------------------------------------------------------
 ! calculates various kinematical LHC observables
 ! input:
 !  mN(:) ......... neutralino masses
 !  mSlepton(:) ... slepton masses
 !  mSdown(6) ....... u-squark masses
 !  mSup(6) ....... u-squark masses
 ! written by Werner Porod, 5.3.2010
 !--------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: mN(:), mSlepton(:), mSup(6), mSdown(6), gauge(3) &
      & , vevSM(2), mGlu, mC(:), mS0(:), RS0(:,:), mP0(:), RP0(:,:)    &
      & , mSpm(:)
  Complex(dp), Intent(in) :: N(:,:), PhaseGlu, RSup(6,6), RSdown(6,6), U(:,:) &
      & , V(:,:), RSpm(:,:), mu, Rslepton(:,:)
  Complex(dp), Intent(in), Dimension(3,3) :: Y_u, Y_d, A_u, A_d
  Logical, intent(in) :: GenerationMixing
  Real(dp), intent(out) :: observ(:)

  Real(dp) :: mSq_av_2, mSlepton2(Size(mSlepton)), lq_low, lq_near, lq_far 
  Integer :: i1, i2, i_n, i_count

  i_n = size(mN)
  mSlepton2 = mSlepton**2

  If (GenerationMixing) then
   observ = 0._dp
  Else
   !---------------------------------------
   ! e+ e- edges
   !---------------------------------------
   i_count = 1
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(1),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(2),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
   !---------------------------------------
   ! mu+ mu- edges
   !---------------------------------------
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(3),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(4),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
   !---------------------------------------
   ! tau+ tau- edges
   !---------------------------------------
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(5),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
   Do i2=2,i_n
    Do i1=1,i2-1
     observ(i_count) = Mll_max(Abs(mN(i2)),mSlepton(6),Abs(mN(i1)))
     i_count = i_count + 1
    end do
   end do
  !-------------------------------------------------------------
  ! l+ l- q edge, averaging over d_L, s_L, u_L, c_L = m^max_llq
  !-------------------------------------------------------------
  mSq_av_2 = (mSup(2)+mSup(4)+mSdown(2)+mSdown(4))/4._dp
  observ(i_count) = Mllq_max(mSq_av_2,mSlepton(1),Abs(mN(2)),Abs(mN(1)))
  i_count = i_count + 1
  !----------------------------------------------------------------
  ! l+ l- q thres, averaging over d_L, s_l, u_L, c_L = m^min_llq
  !----------------------------------------------------------------
  observ(i_count) = Mllq_min(mSq_av_2,mSlepton(1),Abs(mN(2)),Abs(mN(1)))
  i_count = i_count + 1
  !------------------------------------------------------
  ! l+-_near q edge, averaging over d_L, s_l, u_L, c_L
  !------------------------------------------------------
  lq_near = Mqlnear_max( mSq_av_2,Abs(mN(2)),mSlepton(1))
  !-------------------------------------------------
  ! l+-_far q edge, averaging over d_L, s_l, u_L, c_L
  !-------------------------------------------------
  lq_far = Mqlfar_max( mSq_av_2,Abs(mN(2)),mSlepton(1),Abs(mN(1)))

  !------------------------------
  ! l q low edge
  !------------------------------
  lq_low  = Mql_low( mSq_av_2,Abs(mN(2)),mSlepton(1),Abs(mN(1)))
  observ(i_count) = Min(lq_near,lq_low)
  i_count = i_count + 1
  !------------------------------
  ! l q high edge
  !------------------------------
  observ(i_count) = Max(lq_near,lq_far) 
  i_count = i_count + 1

  !-------------------------------------------------
  ! l+ l- b thres
  !-------------------------------------------------
  observ(i_count) = Mllq_min(mSdown(5),mSlepton(1),Abs(mN(2)),Abs(mN(1)))
  i_count = i_count + 1
  !-----------------
  ! m_sq_R - mSl
  !----------------
  observ(i_count) = 0.25_dp * (msdown(1)+msdown(3)+msup(1)+msup(3)) - mN(1)
  i_count = i_count + 1
  observ(i_count) = mSlepton(2) - mN(1)
  i_count = i_count + 1

  !---------------------------------------
  ! M_tb by Hisano et al., hep-ph/0304214
  !---------------------------------------
!   observ(i_count) = M_w_tb(gauge, Y_u, Y_d, A_u, A_d, vevSM, mu, mGlu, PhaseGlu &
!     & , mSup, RSup, mSdown, RSdown, mC, U, V, mN, N, mS0, RS0, mP0, RP0    &
!     & , mSpm, RSpm )
  End If


 end Subroutine Calc_LHC_observables

 Subroutine Write_LHC_observables(io, observ, names)
 implicit none
  Integer, Intent(in) :: io
  Real(dp), intent(in) :: observ(:) 
  Character(len=*), intent(in) :: names(:)
  
  Write(io,100) "Block LHCobservables"


100 Format(a)

 end Subroutine Write_LHC_observables

  Real(dp) Function Mll_max(mchi2, msl, mchi1)
  !----------------------------------------------------
  ! function to calculate m^2_ll^max
  ! written by Werner Porod, 5.4.2008
  !----------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mchi2, msl, mchi1
   Real(dp) :: mchi22, msl2, mchi12

   If ((mchi2.Le.msl).Or.(msl.Le.mchi1)) Then
    Mll_max = -1._dp
    
   Else
    mchi22 = mchi2**2
    mchi12 = mchi1**2
    msl2 = msl**2
    Mll_max = mchi2 * Sqrt( (1._dp - msl2/mchi22) * (1._dp - mchi12/msl2) )
   End If
  End Function Mll_max


  Real(dp) Function Mllq_max(msq, msl, mchi2, mchi1)
  !-------------------------------------------------------------
  ! calculates the llq edge
  ! written by Werner Porod, 5.4.2008
  !-------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: msq, mchi2, msl, mchi1
   Real(dp) :: mchi22, msl2, mchi12, msq2, w1, w2, w3

   If ((msq.Le.mchi2).Or.(mchi2.Le.msl).Or.(msl.Le.mchi1)) Then
    Mllq_max = -1._dp
    
   Else
    mchi22 = mchi2**2
    mchi12 = mchi1**2
    msl2 = msl**2
    msq2 = msq**2

    If (       (msl2.Lt.msq*mchi1).And.(msq*mchi1.Lt.mchi22) &
       & .And. (mchi22*mchi1.Lt.msq*msl2) ) Then
     Mllq_max = msq - mchi1

    Else

     w1 = Mll_max(msq,mchi2,mchi1)
     w2 = Mll_max(msq,msl,mchi1)
     w3 = (msq2*msl2-mchi22*mchi12) * (mchi22-msl2) / (mchi22*msl2)
     If (w3.Gt.0._dp) w3 = sqrt(w3)

     Mllq_max = Max(w1,w2,w3)
    End If

   End If

  End Function Mllq_max


  Real(dp) Function Mllq_min(msq, msl, mchi2, mchi1)
  !-------------------------------------------------------------
  ! calculates the llq threshold 
  ! written by Werner Porod, 5.4.2008
  !-------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: msq, mchi2, msl, mchi1
   Real(dp) :: mchi22, msl2, mchi12, msq2, wert

   If ((msq.Le.mchi2).Or.(mchi2.Le.msl).Or.(msl.Le.mchi1)) Then
    Mllq_min = -1._dp
    
   Else
    mchi22 = mchi2**2
    mchi12 = mchi1**2
    msl2 = msl**2
    msq2 = msq**2

    wert = ((mchi22+msl2)*(msl2+mchi12))**2 - 16._dp*mchi22*msl2**2*mchi12
    
    Mllq_min = ( 2._dp*msl2*(msq2-mchi22)*(mchi22-mchi12)     &
             & + (msq2+mchi22)*(mchi22-msl2)*(msl2-mchi12)    &
             & - (msq2-mchi22)*Sqrt(wert)                      &
             & ) / (4._dp*msl2*mchi22)

    Mllq_min = Sqrt(Mllq_min)
   End If

  End Function Mllq_min


  Real(dp) Function MmaxHQ(mh, mq, mchi1, mchi2)
  !--------------------------------------------------------
  ! calculates the invariant mass of q and X appearing
  ! in the chain
  !        ~q -> q chi^0_2 -> q X chi^0_1
  ! where X is either a Z- or a higgs boson
  ! written by Werner Porod, 5.4.2008
  !--------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mh, mq, mchi1, mchi2
   Real(dp) :: wert, mh2, mq2, mchi12, mchi22

   If ((mq.Lt.Abs(mchi2)).Or.(Abs(mchi2).Lt.(mh+abs(mchi1)))) then
    MmaxHQ = -1._dp
    Return
   End If    

   mh2 = mh**2
   mq2 = mq**2
   mchi12 = mchi1**2
   mchi22 = mchi2**2

   wert = (mchi22 - mh2 - mchi12)**2 - 4._dp * mchi12 * mh2

   wert = (mchi22 + mh2 - mchi12 + wert) * 0.5_dp / mchi22

   MmaxHQ = Sqrt(mh2 + (mq2 - mchi22) * wert)

  End Function MmaxHQ


  Real(dp) Function Mqlfar_max(msq, mchi2, msl, mchi1)
  !----------------------------------------------------
  ! function to calculate m^2_ql_far^max
  ! written by Werner Porod, 5.4.2008
  !----------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: msq, mchi2, msl, mchi1
   Real(dp) :: msq2, mchi22, msl2, mchi12

   If ((msq.Le.mchi2).Or.(mchi2.Le.msl).Or.(msl.Le.mchi1)) Then
    Mqlfar_max = -1._dp
   Else
    msq2 = msq**2
    msl2 = msl**2
    mchi12 = mchi1**2
    mchi22 = mchi2**2
    Mqlfar_max = Sqrt( (msq2-mchi22)*(msl2-mchi12) ) / msl
   End If

  End Function Mqlfar_max


  Real(dp) Function Mql_low(msq, mchi2, msl, mchi1)
  !----------------------------------------------------
  ! function to calculate m^2_ql_far^max
  ! written by Werner Porod, 5.4.2008
  !----------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: msq, mchi2, msl, mchi1
   Real(dp) :: msq2, mchi22, msl2, mchi12

   If ((msq.Le.mchi2).Or.(mchi2.Le.msl).Or.(msl.Le.mchi1)) Then
    Mql_low = -1._dp
   Else
    msq2 = msq**2
    msl2 = msl**2
    mchi12 = mchi1**2
    mchi22 = mchi2**2
    Mql_low = Sqrt( (msq2-mchi22)*(msl2-mchi12)  / (2._dp*msl2-mchi12) )
   End If

  End Function Mql_low


  Real(dp) Function Mqlnear_max(msq, mchi, msl)
  !----------------------------------------------------
  ! function to calculate m^2_ql_near^max
  ! written by Werner Porod, 5.4.2008
  !----------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: msq, mchi, msl
   Real(dp) :: msq2, mchi2, msl2

   If ((msq.Le.mchi).Or.(mchi.Le.msl)) Then
    Mqlnear_max = -1._dp
   Else
    msq2 = msq**2
    msl2 = msl**2
    mchi2 = mchi**2
    Mqlnear_max = Sqrt( (msq2-mchi2)*(mchi2-msl2) ) / mchi
   End If

  End Function Mqlnear_max


 Real(dp) Function M_w_tb(g, Y_u, Y_d, A_u, A_d, vevSM, mu, mGlu, PhaseGlu  &
    & , mSup, RSup, mSdown, RSdown, mC, U, V, mN, N, mS0, RS0, mP0, RP0     &
    & , mSpm, RSpm )
 !--------------------------------------------------------------
 ! calculates the weighted edge M^w_tb proposed by J.Hisano et al
 ! hep-ph/0304421
 !--------------------------------------------------------------
  Implicit None
  Real(dp), Intent(in) :: mGlu, msup(6), msdown(6), g(3), mN(4), mC(2)  &
     & , mS0(2), mP0(2), mSpm(2), RS0(2,2), RP0(2,2), vevSM(2)
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6), Y_u(3,3), Y_d(3,3)  &
     & , A_u(3,3), A_d(3,3), mu, U(2,2), V(2,2), N(4,4), RSpm(2,2), PhaseGlu 

  Complex(dp) :: coupC, coupLC, coupRC, cpl_SmpSdSu(2,6,6), RSd(2,2)         &
    & , Rsu(2,2), cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6), cpl_CDSu_L(2,3,6)    &
    & , cpl_CDSu_R(2,3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6)              &
    & , cpl_UNSu_L(3,4,6), cpl_UNSu_R(3,4,6), cpl_GDSd_L(3,6)                &
    & , cpl_GDSd_R(3,6), cpl_GUSu_L(3,6), cpl_GUSu_R(3,6), cpl_SdSuW(6,6)    &
    & , cpl_SuSdW(6,6), cpl_P0SdSd(2,6,6), cpl_P0SuSu(2,6,6)                 &
    & , cpl_S0SdSd(2,6,6), cpl_S0SuSu(2,6,6), cpl_SdSdZ(6,6), cpl_SuSuZ(6,6) &
    & , bi(1), cpl_SmpSuSd(2,6,6)
  Real(dp) :: M_tb_1, M_tb_11, BR_1, BR_11, sinW2, gP_Sd(6,54), gT_Sd(6)   &
    & , BR_Sd(6,54), gP_Su(6,66), gT_Su(6), BR_Su(6,66), gP_G2(36), gT_G   &
    & , BR_G2(36), mSt2, mG2, mC2, mT2, mSb2
  Integer :: i1, i2, i3
  Real(dp), Parameter :: e_d = -1._dp / 3._dp, e_u = 2._dp / 3._dp
  Complex(dp), Parameter :: ZeroC = (0._dp, 0._dp)


  M_W_TB = 0._dp
  !--------------------------------------------
  ! couplings
  !--------------------------------------------
  sinW2 = g(1)**2 / (g(1)**2 + g(2)**2)
  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSdSu = ZeroC
  Rsd = RSdown(5:6,5:6)
  Rsu = RSup(5:6,5:6)
  Do i1=1,2
   Do i2=1,2
    Call CoupChargedScalarSfermion3(2, i1, i2, RSpm, g(2), vevSM, mu, Y_d(3,3)&
                       &,  A_d(3,3), Rsd, Y_u(3,3), A_u(3,3),Rsu,coupC)
     cpl_SmpSdSu(2,4+i1,4+i2) = coupC
     cpl_SmpSuSd(2,4+i2,4+i1) = Conjg(coupC)
    End Do
  End Do
  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  Do i1=1,2
   Do i2=1,2
    Call CoupCharginoSfermion(i1, i2, g(2), 0.5_dp, RSd, Y_d(3,3), Y_u(3,3) &
                             & , U, V, coupLC, coupRC)
    cpl_CUSd_L(i1, 3, 4 + i2) = coupLC
    cpl_CUSd_R(i1, 3, 4 + i2) = coupRC
    Call CoupCharginoSfermion(i1, i2, g(2), -0.5_dp, RSu, Y_d(3,3), Y_u(3,3) &
                            & , U, V, coupLC, coupRC)
    cpl_CDSu_L(i1,3, 4 + i2) = coupLC
    cpl_CDSu_R(i1,3, 4 + i2) = coupRC
   End Do
  End Do
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC
  Do i1=1,2
   Do i2=1,4
    Call CoupNeutralinoSdown(i2, i1, g(1), g(2), RSd, Y_d(3,3), N  &
                            & , coupLC, coupRC)
    cpl_DNSd_L(3, i2, 4 + i1 ) = coupLC
    cpl_DNSd_R(3, i2, 4 + i1 ) = coupRC
    Call CoupNeutralinoSup(i2, i1, g(1), g(2), RSu, Y_u(3,3), N    &
                            & , coupLC, coupRC)
    cpl_UNSu_L(3, i2, 4 + i1 ) = coupLC
    cpl_UNSu_R(3, i2, 4 + i1 ) = coupRC
   End Do
  End Do
  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC
  Do i1 =1,2
   Call CoupGluinoSquark(g(3), phaseglu, i1, Rsd, coupLC, coupRC)
   cpl_GDSd_L(3, 4 + i1) = coupLC
   cpl_GDSd_R(3, 4 + i1) = coupRC
   Call CoupGluinoSquark(g(3), phaseglu, i1, Rsu, coupLC, coupRC)
   cpl_GUSu_L(3, 4 + i1) = coupLC
   cpl_GUSu_R(3, 4 + i1) = coupRC
  End Do  
  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  Do i1 = 1,2
   Do i2 = 1,2
    Call CoupSfermionW3(i1, i2, g(2), RSd, RSu, cpl_SdSuW( 4+i1, 4+i2 ) )
   End Do
  End Do
  Call Adjungate(cpl_SdSuW, cpl_SuSdW)
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp
  bi(1) = mu
  Do i1=1,2
   Do i2=1,2
    Call CoupPseudoScalarSfermion3(2, i1, i2, RP0, -0.5_dp, Y_d(3,3), Rsd   &
                                   &, A_d(3,3) , bi, coupC )
    cpl_P0SdSd(2,4+i1,4+i2) = coupC
    Call CoupPseudoScalarSfermion3(2, i1, i2, RP0, 0.5_dp, Y_u(3,3), Rsu   &
                                   &, A_u(3,3), bi, coupC )
    cpl_P0SuSu(2,4+i1,4+i2) = coupC
   End Do
  End Do
  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp
  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
     Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d(3,3), Rsd   &
                            &, A_d(3,3), mu, vevSM, g(1), g(2), coupC )
     cpl_S0SdSd(i1,4+i2,4+i3) = coupC
     Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u(3,3), Rsu   &
                               &, A_u(3,3), mu, vevSM, g(1), g(2), coupC )
     cpl_S0SuSu(i1,4+i2,4+i3) = coupC
    End Do
   End Do
  End Do
  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SuSuZ = 0._dp
  Do i1=1,2
   Do i2=1,2
    Call CoupSdownZ(i1, i2, g(2), sinW2, Rsd, cpl_SdSdZ( 4+i1, 4+i2 ) )
    Call CoupSupZ(i1, i2, g(2), sinW2, Rsu, cpl_SuSuZ( 4+i1, 4+i2 ) )
   End Do
  End Do


  !--------------------------------------------
  ! branching ratios
  !--------------------------------------------
  gP_Su = 0._dp
  gT_Su = 0._dp
  BR_Su = 0._dp
  Call SfermionTwoBodyDecays_old(5, mSup, mf_u, mf_d                           &
          &, mN, cpl_UNSu_L, cpl_UNSu_R, mC, cpl_CDSu_L, cpl_CDSu_R         &
          &, mSdown, mW, cpl_SuSdW, mZ, cpl_SuSuZ, mSpm, cpl_SmpSuSd        &
          &, mP0, cpl_P0SuSu, ms0, cpl_S0SuSu                               &
          &, 0, GenerationMixing                                            &
          &, gP_Su, gT_Su, BR_Su, mGlu, cpl_GUSu_L, cpl_GUSu_R)

  gP_G2 = 0._dp
  gT_G = 0._dp
  BR_G2 = 0._dp
  Call GluinoTwoBodyDecays_old(mGlu, mSdown, cpl_GDSd_L        &  
          &, cpl_GDSd_R, mf_d, mSup, cpl_GUSu_L, cpl_GUSu_R, mf_u         &  
          &, 0, GenerationMixing, gP_G2, gT_G, BR_G2 )
  !-----------------------------
  ! only 3-body decays
  !-----------------------------
  If (gT_G.Eq.0._dp) Return

  gP_Sd = 0._dp
  gT_Sd = 0._dp
  BR_Sd = 0._dp
  Call SfermionTwoBodyDecays_old(5, mSdown, mf_d, mf_u                         &
          &, mN, cpl_DNSd_L, cpl_DNSd_R, mC, cpl_CUSd_L, cpl_CUSd_R         &
          &, mSup, mW, cpl_SdSuW, mZ, cpl_SdSdZ, mSpm, cpl_SmpSdSu          &
          &, mP0, cpl_P0SdSd, mS0, cpl_S0SdSd                               &
          &, 0, GenerationMixing                                            &
          &, gP_Sd, gT_Sd, BR_Sd, mGlu, cpl_GDSd_L, cpl_GDSd_R)
  Call SfermionTwoBodyDecays_old(6, mSdown, mf_d, mf_u                         &
          &, mN, cpl_DNSd_L, cpl_DNSd_R, mC, cpl_CUSd_L, cpl_CUSd_R         &
          &, mSup, mW, cpl_SdSuW, mZ, cpl_SdSdZ, mSpm, cpl_SmpSdSu          &
          &, mP0, cpl_P0SdSd, mS0, cpl_S0SdSd                               &
          &, 0, GenerationMixing                                            &
          &, gP_Sd, gT_Sd, BR_Sd, mGlu, cpl_GDSd_L, cpl_GDSd_R)


  mSt2 = mSup(5)**2
  mG2 = mGlu**2
  mC2 = mC(1)**2
  mT2 = mf_u2(3)

  M_tb_1 = (mG2 - (mSup(5)-mf_u(3))**2) * (mG2 - (mSup(5)+mf_u(3))**2)
  If ( M_tb_1.Lt.0._dp) Return
  M_tb_1 = mT2 + 0.5_dp * (mSt2-mC2) * (mG2 -mSt2 -mT2 + Sqrt( M_tb_1)) / mSt2
  If ( M_tb_1.Lt.0._dp) Return
  M_tb_1 = Sqrt(M_tb_1)

  mSb2 = mSdown(5)**2
  M_tb_11 = (mSb2 - (mC(1)-mf_u(3))**2) * (mSb2 - (mC(1)+mf_u(3))**2)
  If ( M_tb_11.Lt.0._dp) Return
  M_tb_11 = mT2 + 0.5_dp * (mG2-mSb2) * (mSb2 -mC2 +mT2 + Sqrt( M_tb_11)) /mSb2
  If ( M_tb_11.Lt.0._dp) Return
  M_tb_11 = Sqrt(M_tb_11)

  BR_1 = BR_G2(5) * BR_Su(5,5)  &
     & + BR_G2(11) * BR_Sd(5,8) * BR_Su(5,5) &
     & + BR_G2(12) * BR_Sd(6,8) * BR_Su(5,5)
  BR_11 = BR_G2(11) * BR_Sd(5,5)
  
  M_w_tb = 0._dp
  if ((BR_1 + BR_11).gt.0._dp) &
       & M_w_tb = ( BR_1 * M_tb_1 + BR_11 * M_tb_11) / (BR_1 + BR_11)

 End Function M_w_tb









End Module LHC_observables
