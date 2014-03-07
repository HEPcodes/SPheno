Module Stop3BodyDecays
!---------------------------------------------------------------
! 23.11: adding decay modes as given in Djouadi and Mambrini
!        PRD63, 115005 (2001) hep-ph/0011364
!---------------------------------------------------------------
! load modules
Use Control
Use DecayFunctions, only: ScalarToTwoFermions
Use Mathematics
Use StandardModel, Only: CKM, mW, mW2, mZ2, mf_u, mf_u2, mf_d, mf_d2  &
  & , mf_l, mf_l2, mZ
Use ThreeBodyPhaseSpaceS
! load modules

! private variables
Integer, private :: istosle
Real(dp), Private :: mL2, mb2, mSl2, mst2, gammc(2), gammc2(2), gammsb(2)   &
         & , gammsb2(2), gamtmt, gamtmt2, amc2(2), SumM(3), amN2, msbot2(2) &
         & , mt2
Complex(dp), Private :: CoeffSl(3,5), ca(3,9), cb(2,8), cc(2,2,7), cd(3,6) &
                      & , ce(3,2)

! private routines
Private :: Ftilde
! private variables

Contains



 Real(dp) Function Ftilde(muChi, muPhi, muSq)
 implicit none
  Real(dp), Intent(in) :: muChi, muPhi, muSq
  
  Real(dp) :: parts(2), rp, rm

  parts(1) = 1._dp + muPhi - MuSq
  parts(2) = Sqrt( (1._dp - muPhi - muSq)**2 - 4._dp * muPhi * muSq)
  rp = 0.5_dp * (parts(1) + parts(2) )
  rm = 0.5_dp * (parts(1) - parts(2) )

  FTilde = - f(1._dp) + f(rp) + f(rm) &
       & + Log(muPhi) * Log((muChi-muSq)/(muChi-1._dp) )

 contains
  real(dp) function f(x)
  implicit none
   Real(dp), intent(in) :: x
    f = Log(muChi - 1._dp + x) * Log((1._dp - muChi)/(muSq-muChi))  &
      & + Li2((muChi-muSq)/(muchi - 1._dp + x) )                    &
      & - Li2((muChi-1._dp)/(muchi - 1._dp + x) )
  end function 
 End Function Ftilde

 
 Real(dp) Function Sq2toSq1Chi1(muChi, muSq)
 implicit none
  Real(dp), intent(in) :: muChi(2), muSq

  Real(dp) :: lam1, sq_lam1, lam2, sq_lam2, ProdChi

  if (muChi(1).eq.muChi(2)) then
   lam1 = -1._dp + 2._dp * (muChi(1) + muSq) - (muSq - muChi(1))**2
   sq_lam1 = Sqrt(abs(lam1))

   Sq2toSq1Chi1 = (muSq - 1._dp) * (5._dp - 6._dp * muChi(1)                 &
              &                    + 5._dp * muSq- 2._dp * muSq / muChi(1) ) &
              &  - 2._dp * (muSq/muchi(1))**2 * Log( muSq )                  &
              &  + 2._dp * (muSq-muChi(1)) * (muChi(1)-1._dp)                &
              &    * (muSq + muSq*muChi(1) + muChi(1) - 3._dp*muChi(1)**2 )  &
              &    * Log((muChi(1)-1._dp)/(muChi(1)-muSq)) / muChi(1)**2

  else
   lam1 = -1._dp + 2._dp * (muChi(1) + muSq) - (muSq - muChi(1))**2
   sq_lam1 = Sqrt(abs(lam1))
   lam2 = -1._dp + 2._dp * (muChi(2) + muSq) - (muSq - muChi(2))**2
   sq_lam2 = Sqrt(abs(lam2))
   prodChi = muChi(1) * muChi(2)
   Sq2toSq1Chi1 = (muSq - 1._dp) * (3._dp +3._dp * muSq -2._dp * Sum(muChi))  &
              & - 2._dp * muSq**2/ProdChi * Log( muSq )                       &
              &  + ( (muSq-muChi(1))**2 * (muChi(1)-1._dp)**2                 &
              &            * Log((muChi(1)-muSq)/(muChi(1)-1._dp)) / muchi(1) &
              &    -muChi(1) * (muSq-muChi(2))**2 * (muChi(2)-1._dp)**2       &
              &        * Log((muChi(2)-muSq)/(muChi(2)-1._dp)) / muchi(2)     &
              &    ) / (0.5_dp * (muChi(2)-muChi(1)))
  end if

  Sq2toSq1Chi1 = 0.25_dp * Sq2toSq1Chi1

 end Function Sq2toSq1Chi1

 Real(dp) Function Sq2toSq1Chi2(muChi, muSq)
 implicit none
  Real(dp), intent(in) :: muChi(2), muSq

  Real(dp) :: lam1, sq_lam1, lam2, sq_lam2, prodChi

  if (muChi(1).eq.muChi(2)) then
   lam1 = -1._dp + 2._dp * (muChi(1) + muSq) - (muSq - muChi(1))**2
   sq_lam1 = Sqrt(abs(lam1))
   Sq2toSq1Chi2 = 0.5_dp * (muSq - 1._dp)                                   &
              &   * (muchi(1) - 2._dp * muchi(1)**2 - 2._dp * muSq          &
              &     + muSq * muChi(1) )                                     &
              & + muSq * (1._dp + muSq - muSq/muChi(1) ) * Log( muSq )      &
              & + (muSq- muChi(1)) * (muChi(1)-1._dp) * (muSq-muChi(1)**2)  &
              &   * Log( (muChi(1)-1._dp) / (muChi(1)-muSq) ) / muChi(1)

  else
   lam1 = -1._dp + 2._dp * (muChi(1) + muSq) - (muSq - muChi(1))**2
   sq_lam1 = Sqrt(abs(lam1))
   lam2 = -1._dp + 2._dp * (muChi(2) + muSq) - (muSq - muChi(2))**2
   sq_lam2 = Sqrt(abs(lam2))
   prodChi = muChi(1) * muChi(2)
   Sq2toSq1Chi2 = muSq * (1._dp + muSq - 0.5_dp *muSq * Sum(muChi) / ProdChi) &
             &         * Log(muSq)                                            &
             &  + 0.5_dp * (muSq - 1._dp) * (muSq + ProdChi)                  &
             &  + (muChi(2) * (muSq-muChi(1))**2 * (muChi(1)-1._dp)**2        &
             &    * Log((muChi(1)-muSq)/(muChi(1)-1._dp)) / muchi(1)          &
             &    -muChi(1) * (muSq-muChi(2))**2 * (muChi(2)-1._dp)**2        &
             &    * Log((muChi(2)-muSq)/(muChi(2)-1._dp)) / muchi(2)          &
             &    ) / (2._dp * (muChi(2)-muChi(1)) )

  end if

  Sq2toSq1Chi2 = Sq2toSq1Chi2 / (muChi(1)*muChi(2))

 end Function Sq2toSq1Chi2

 Subroutine Sq2ToSq1ff(mSq, C_Sq2Sq1Z, Lf, Rf, mS0, c_Sq2Sq1S0, c_FFS0 &
           & , mP0, c_Sq2Sq1P0, c_FFP0, mC, C_CFSfp_L, C_CFSfp_R, mN   &
           & , C_FNSf_L, C_FNSf_R, mGlu, C_FGSf_L, C_FGSf_R, Nc        &
           & , SameGeneration, width)
 !--------------------------------------------------------------------------
 ! calculates the decay sq_2 -> sq_1 f \bar{f}
 ! uses the formulas given in A. Djouadi and Y. Mambrini, hep-ph/0011364
 ! 23.11: implementing gauge boson exchange
 !--------------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: mSq(2), Lf, Rf, mS0(2), mP0(2), c_FFS0(2)  &
       &  , c_FFP0(2), mC(2), mN(4), mGlu, Nc
  Complex(dp), Intent(in) :: C_Sq2Sq1Z, c_Sq2Sq1S0(2), c_Sq2Sq1P0(2)     &
       & , C_CFSfp_L(2,2), C_CFSfp_R(2,2), C_FNSf_L(4,2), C_FNSf_R(4,2)  &
       & , C_FGSf_L(2), C_FGSf_R(2)
  logical :: SameGeneration(2)
  Real(dp), intent(out) :: width

  Integer :: i1, i2, i_count
  Real(dp) :: parts(16), muV, muSq, muPhi(2)

  width = 0._dp
  parts = 0._dp
  !------------------
  ! Z boson exchange
  !------------------
  muV = mZ / mSq(2)
  muSq = mSq(1) / mSq(2)
  parts(1) = 2._dp * (Rf**2 + Lf**2 ) * Abs(C_Sq2Sq1Z)**2 &
         &         * Sq2ToSq1VV(muV, muSq)
  !-------------------------------------------------------------------------
  ! h^0 h^0 exchange, note that Higgs-F-F coupling is defined differently
  !-------------------------------------------------------------------------
  muPhi = mS0(1) / mSq(2)
  parts(2) = 4._dp * c_FFS0(1)**2 * Abs(c_Sq2Sq1S0(1))**2    &
           & * Sq2ToSq1PhiPhi(muPhi, muSq)
  !-------------------------------------------------------------------------
  ! h^0 H^0 exchange, note that Higgs-F-F coupling is defined differently
  !-------------------------------------------------------------------------
  muPhi(1) = mS0(1) / mSq(2)
  muPhi(2) = mS0(2) / mSq(2)
  parts(3) = 8._dp * c_FFS0(1) * c_FFS0(2)                    &
           & * Real(c_Sq2Sq1S0(1)*Conjg(c_Sq2Sq1S0(2)),dp)**2    &
           & * Sq2ToSq1PhiPhi(muPhi, muSq)
  !-------------------------------------------------------------------------
  ! H^0 H^0 exchange, note that Higgs-F-F coupling is defined differently
  !-------------------------------------------------------------------------
  muPhi = mS0(2) / mSq(2)
  parts(4) = 4._dp * c_FFS0(2)**2 * Abs(c_Sq2Sq1S0(2))**2    &
           & * Sq2ToSq1PhiPhi(muPhi, muSq)
  !-------------------------------------------------------------------------
  ! A^0 A^0 exchange, note that Higgs-F-F coupling is defined differently
  !-------------------------------------------------------------------------
  muPhi = mP0(2) / mSq(2)
  parts(5) = 4._dp * c_FFP0(2)**2 * Abs(c_Sq2Sq1P0(2))**2    &
           & * Sq2ToSq1PhiPhi(muPhi, muSq)
  !--------------------------------
  ! chargino contribution
  !--------------------------------
  If (SameGeneration(1)) then
   muPhi = mC(1) / mSq(2)
   parts(6) = ( Abs( C_CFSfp_L(1,1) ) * Abs( C_CFSfp_L(1,2) )     &
          &   + Abs( C_CFSfp_R(1,1) ) * Abs( C_CFSfp_R(1,2) )  )  &
          &  * Sq2toSq1Chi1(muPhi, muSq)                          &
          & + Real( C_CFSfp_L(1,2) * C_CFSfp_L(1,1)               &
          &         * Conjg( C_CFSfp_R(1,2) * C_CFSfp_R(1,1) ) ,dp)  &
          &  * Sq2toSq1Chi2(muPhi, muSq)
   muPhi = mC / mSq(2)
   parts(7) = ( Real( C_CFSfp_L(1,1) * C_CFSfp_L(2,2)                &
          &          * Conjg( C_CFSfp_L(1,2) * C_CFSfp_L(2,1) ),dp )    &
          &   + Real( C_CFSfp_R(1,1) * C_CFSfp_R(2,2)                &
          &          * Conjg( C_CFSfp_R(1,2) * C_CFSfp_R(2,1) ),dp )  ) &
          & * Sq2toSq1Chi1(muPhi, muSq)                              &
          & + ( Real( C_CFSfp_L(1,2) * C_CFSfp_L(2,1)                &
          &         * Conjg( C_CFSfp_R(1,1) * C_CFSfp_R(2,2) ),dp )     &
          &   +  Real( C_CFSfp_R(1,2) * C_CFSfp_R(2,1)               &
          &         * Conjg( C_CFSfp_L(1,1) * C_CFSfp_L(2,2) ),dp) )   &
          &  * Sq2toSq1Chi2(muPhi, muSq)
   parts(7) = 2._dp * parts(7)
   muPhi = mC(2) / mSq(2)
   parts(8) = ( Abs( C_CFSfp_L(2,1) ) * Abs( C_CFSfp_L(2,2) )     &
          &   + Abs( C_CFSfp_R(2,1) ) * Abs( C_CFSfp_R(2,2) )  )  &
          & * Sq2toSq1Chi1(muPhi, muSq)                           &
          & + Real( C_CFSfp_L(2,2) * C_CFSfp_L(2,1)               &
          &         * Conjg( C_CFSfp_R(2,2) * C_CFSfp_R(2,1) ),dp )  &
          &  * Sq2toSq1Chi2(muPhi, muSq)
  end if
  !--------------------------------
  ! neutralino contribution
  !--------------------------------
  If (SameGeneration(2)) then
   i_count = 5
   Do i1=1,4
    muPhi = mN(i1) / mSq(2)
    parts(i_count) = ( Abs( C_FNSf_L(i1,1) ) * Abs( C_FNSf_L(i1,2) )     &
                 &   + Abs( C_FNSf_R(i1,1) ) * Abs( C_FNSf_R(i1,2) )  )  &
                 & * Sq2toSq1Chi1(muPhi, muSq)                           &
                 & + Real( C_FNSf_L(i1,2) * C_FNSf_L(i1,1)               &
                 &         * Conjg( C_FNSf_R(i1,2) * C_FNSf_R(i1,1) ),dp )  &
                 &  * Sq2toSq1Chi2(muPhi, muSq)
    i_count = i_count + 1
    Do i2=i1+1,4
     muPhi(1) = mN(i1) / mSq(2)
     muPhi(2) = mN(i2) / mSq(2)
     parts(i_count) = ( Real( C_FNSf_L(i1,1) * C_FNSf_L(i2,2)        &
          &          * Conjg( C_FNSf_L(i1,2) * C_FNSf_L(i2,1) ) )    &
          &   + Real( C_FNSf_R(i1,1) * C_FNSf_R(i2,2)                &
          &          * Conjg( C_FNSf_R(i1,2) * C_FNSf_R(i2,1) ),dp)) &
          & * Sq2toSq1Chi1(muPhi, muSq)                              &
          & + ( Real( C_FNSf_L(i1,2) * C_FNSf_L(i2,1)                &
          &         * Conjg( C_FNSf_R(i1,1) * C_FNSf_R(i2,2) ),dp )     &
          &   +  Real( C_FNSf_R(i1,2) * C_FNSf_R(i2,1)               &
          &         * Conjg( C_FNSf_L(i1,1) * C_FNSf_L(i2,2) ),dp ) )   &
          &  * Sq2toSq1Chi2(muPhi, muSq)
     parts(i_count) = 2._dp * parts(i_count)
     i_count = i_count + 1
    end do
   end do
  end if
  !----------------------
  ! gluino contribution
  !----------------------
  If (SameGeneration(2)) then
   muPhi = mGlu / mSq(2)
   parts(i_count) = ( Abs( C_FGSf_L(1) ) * Abs( C_FGSf_L(2) )     &
                &   + Abs( C_FGSf_R(1) ) * Abs( C_FGSf_R(2) )  )  &
                & * Sq2toSq1Chi1(muPhi, muSq)                           &
                & + Real( C_FGSf_L(2) * C_FGSf_L(1)               &
                &         * Conjg( C_FGSf_R(2) * C_FGSf_R(1) ),dp )  &
                &  * Sq2toSq1Chi2(muPhi, muSq)
  end if

  width = oo128pi3 * mSq(2) * Nc * Sum(parts)

 end subroutine Sq2ToSq1ff


 Real(dp) Function Sq2ToSq1PhiChi(muChi, muPhi, muSq)
 !-------------------------------------------------------------------
 ! calculates the gaugino - scalar interference for sfermion decays
 ! input: muChi .... m_gaugino / m_sfermion_2
 !        muPhi .... m_scalar / m_sfermion_2
 !        muSq ..... m_sfermion_1 /  m_sfermion_2
 ! uses: function Ftilde
 ! written by Werner Porod, 09.12.02
 !-------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: muChi, muPhi, muSq

  Sq2ToSq1PhiChi = muSq - 1._dp + muSq * Log(muSq) / muChi  &
        & + (muSQ-muchi)*(muChi-1._dp)/muChi                &
        &    *  Log( (muChi-1._dp) / (muChi-muSq) )         &
        & + muPhi * Ftilde(muChi, muPhi, muSq) 

 End Function Sq2ToSq1PhiChi

 Real(dp) Function Sq2ToSq1PhiPhi(muPhi, muSq)
 implicit none
  Real(dp), intent(in) :: muPhi(2), muSq

  Real(dp) :: lam1, sq_lam1, lam2, sq_lam2

  if (muPhi(1).eq.muPhi(2)) then
   lam1 = -1._dp + 2._dp * (muPhi(1) + muSq) - (muSq - muPhi(1))**2
   sq_lam1 = Sqrt(abs(lam1))
   Sq2ToSq1PhiPhi = 2._dp * (muSq - 1._dp)                                 & 
                & + 0.5_dp * (1._dp + muSq - Sum(muPhi)) * Log(muSq)       &
                & - ( muPhi(1) * (1._dp+muSq-muPhi(1)) + lam1)             &
                &   * ATan( (1-muSq) * sq_lam1                             &
                &         / (muPhi(1)*(1._dp - muPhi(1) + muSq) - lam1) )  &
                &   / sq_lam1
  else
   lam1 = -1._dp + 2._dp * (muPhi(1) + muSq) - (muSq - muPhi(1))**2
   sq_lam1 = Sqrt(abs(lam1))
   lam2 = -1._dp + 2._dp * (muPhi(2) + muSq) - (muSq - muPhi(2))**2
   sq_lam2 = Sqrt(abs(lam2))
   Sq2ToSq1PhiPhi = muSq - 1._dp                                           &
                & + 0.5_dp * (1._dp + muSq - Sum(muPhi)) * Log(muSq)       &
                & + ( muPhi(1) * sq_lam1                                   &
                &   * ATan( (1-muSq) * sq_lam1                             &
                &         / (muPhi(1)*(1._dp - muPhi(1) + muSq) - lam1) )  &
                &   - muPhi(2) * sq_lam2                                   &
                &   * ATan( (1-muSq) * sq_lam2                             &
                &         / (muPhi(2)*(1._dp - muPhi(2) + muSq) - lam2) )  &
                &   ) / (muPhi(2) - muPhi(1))
   
  end if

 end Function Sq2ToSq1PhiPhi


 Real(dp) Function Sq2ToSq1VChi(muChi, muV, muSq)
 !-------------------------------------------------------------------
 ! calculates the gaugino - vector boson interference for sfermion decays
 ! input: muChi .... m_gaugino / m_sfermion_2
 !        muV ...... m_V / m_sfermion_2
 !        muSq ..... m_sfermion_1 /  m_sfermion_2
 ! uses: function Ftilde
 ! written by Werner Porod, 09.12.02
 !-------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: muChi, muV, muSq

  real(dp) :: lam, sqlam

  lam = -1._dp + 2._dp * (muV + muSq) - (muV - muSq)**2
  sqlam = Sqrt(lam)
  Sq2ToSq1VChi = (muSq - 1._dp) * ( 0.25_dp * (1._dp+muSQ+2._dp*muV)          &
        &                         + (1._dp+muSq-2._dp*(muChi+muV))*Log(muV))  &
        & - (muSQ-muchi)*(muChi-1._dp) *  Log( (muChi-1._dp) / (muChi-muSq) ) &
        & - 0.25_dp * (lam + 2._dp * muchi * (1._dp - muV + muSq)             &
        &             - 6._dp * muSq) * Log(muSq)                             &
        & + 0.5_dp * sqlam * (1._dp + muSq - 2._dp + muChi - muV)             &
        &     * Atan( (1._dp - muSq) * sqlam / (lam + muV*(muV-1._dp-muSq)) ) &
        & + muV * Ftilde(muChi, muV, muSq) 

 End Function Sq2ToSq1VChi

 real(dp) function Sq2ToSq1VV(muV, muSq)
 implicit none
  Real(dp), intent(in) :: muV, muSq
  Real(dp) :: lam, sq_lam

  lam = -1._dp + 2._dp * (muV + muSq) - (muSq - muV)**2
  sq_lam = Sqrt(abs(lam))
  Sq2ToSq1VV = (1-muSq) * ( 5._dp * (1._dp + muSq) - 4._dp * muV  &
           &              + 2._dp * lam / muV) / 3._dp            &
           & + (lam - 2._dp * muSq) * Log( muSq )                 &
           & + 2._dp * (1._dp - muV + muSq) * sq_lam              &
           &         * Atan( (1-muSq) * sq_lam                    &
           &               / (muV*(1._dp - muV + muSq) - lam) )
  Sq2ToSq1VV = 0.25_dp * Sq2ToSq1VV

 end function Sq2ToSq1VV

 Real(dp) Function stbWchi(s)
 !-----------------------------------------------------
 ! Integrand for stop_1 -> b + W+ + neutralino_1
 ! The following functions are called:
 !   -dambda
 !   -I2dsmg1tmg2
 !    I2dtm1g1
 !    I2dtm12g12
 !    I2tdsmg1tmg2
 !    I2tdtm1g1
 !    I2tdtm12g12
 !    I2t2dtm1g1
 !    I2t2dtm12g12
 ! written by Werner Porod
 ! 21.01.02: portation to f90
 !-----------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: s

  Integer j, k
  Real(dp) :: s2, s3, ws, tmin, tmax, vwert, uwert, fakt
  Complex(dp) :: intI(4,3)


  Iname = Iname + 1
  NameOfUnit(Iname) = "StBWChi"

  s2 = s*s
  s3 = s2*s

  ws = kappa(s, mst2, mw2) * kappa(s, amN2, mb2) / (2._dp*s)
 !-------------------------
 ! Range of t-Integration
 !-------------------------
  tmin = 0.5_dp * (sumM(1)-s) - (mst2-mW2)*(amN2-mb2) / (2._dp*s) - ws
  tmax = tmin + 2._dp * ws

  stbWchi = 0._dp
 !--------------------------------------------------
 ! t-Integrals for chargino-exchange,
 ! chargino exchange: j=1 ... chargino_1 chargino_1
 !                    j=2 ... chargino_2 chargino_2
 !                    j=3 ... chargino_1 chargino_2
 !--------------------------------------------------
  Do j=1,2
   vwert = sumM(j+1)-s
   intI(j,1) = I2dtm1g1(tmin,tmax,vwert,gammc2(j))
   intI(j,2) = I2tdtm1g1(tmin,tmax,vwert,gammc2(j))
   intI(j,3) = I2t2dtm1g1(tmin,tmax,vwert,gammc2(j))
  End Do

  vwert = sumM(2)-s
  uwert = sumM(3)-s
  intI(3,1) = I2dtm12g12(tmin,tmax,vwert,gammc(1),uwert,gammc(2))
  intI(3,2) = I2tdtm12g12(tmin,tmax,vwert,gammc(1),uwert,gammc(2))
  intI(3,3) = I2t2dtm12g12(tmin,tmax,vwert,gammc(1),uwert,gammc(2))

!Write(40,*) " "
!Write(40,123) s,ws,tmin,tmax,vwert
!123 Format(5f15.6)
  Do j = 1,3
   stbWchi = stbWchi                                                      &
         & + ( ca(j,1) + ca(j,2)*s + ca(j,3)*s2 + ca(j,4)*s3) * intI(j,1) &
         & + ( ca(j,5) + ca(j,6)*s + ca(j,7)*s2) * intI(j,2)              &
         & + ( ca(j,8) + ca(j,9)*s) * intI(j,3)
  End Do
 !-------------------
 ! chargino_j top interference term
 !-------------------
  Do j=1,2
   uwert = sumM(j+1) - s
   intI(j,1) = I2dtm12g12(tmin,tmax,uwert,-gammc(j),mt2,gamtmt)
   intI(j,2) = I2tdtm12g12(tmin,tmax,uwert,-gammc(j),mt2,gamtmt)
   intI(j,3) = I2t2dtm12g12(tmin,tmax,uwert,-gammc(j),mt2,gamtmt)

 !--------------------------------------------------
 ! neagitive sign because chargino is u-channel
 !--------------------------------------------------
   stbWchi = stbWchi - ( cb(j,1) + cb(j,2)*s + cb(j,3)*s2) * intI(j,1)   &
           &         - ( cb(j,4) + cb(j,5)*s + cb(j,6)*s2) * intI(j,2)   &
           &         - ( cb(j,7) + cb(j,8)*s) * intI(j,3) 
  End Do
 !-------------------
 ! chargino_j sbottom_k interference term
 !-------------------
  Do j=1,2
   uwert = sumM(j+1) - s
   Do k=1,2
    intI(j,1) = I2dsmg1tmg2(s,tmin,tmax,msbot2(k),gammsb(k),uwert,-gammc(j))
    intI(j,2) = I2tdsmg1tmg2(s,tmin,tmax,msbot2(k),gammsb(k),uwert,-gammc(j))

 !--------------------------------------------------
 ! negative sign because chargino is u-channel
 !--------------------------------------------------
    stbWchi = stbWchi                                                       &
          & - (cc(k,j,1)+cc(k,j,2)*s+cc(k,j,3)*s2+cc(k,j,4)*s3) * intI(j,1) &
          & - ( cc(k,j,5) + cc(k,j,6)*s + cc(k,j,7)*s2) * intI(j,2)
   End Do
  End Do
 !--------------------------------------------------
 ! t-Integrals for top-top term,
 !--------------------------------------------------
  intI(1,1) = I2dtm1g1(tmin,tmax,mt2,gamtmt2)
  intI(1,2) = I2tdtm1g1(tmin,tmax,mt2,gamtmt2)
  intI(1,3) = I2t2dtm1g1(tmin,tmax,mt2,gamtmt2)
!j=1
!Write(40,*) j,Real(inti(j,:))
!Write(*,*) "s,tmin,tmax,mt2,gamtmt2,intI(1,:)"
!Write(42,123) s,tmin,tmax,mt2,gamtmt2,real(intI(1,:))
!123 format(4e16.7)
!write(42,*) " "
!stop 99
  stbWchi = stbWchi + (cd(1,1) + cd(1,2)*s ) * intI(1,1)  &
          &         + (cd(1,3) + cd(1,4)*s ) * intI(1,2)  &
          &         + (cd(1,5) + cd(1,6)*s ) * intI(1,3)
 !-------------------
 ! top-sbottom_k interference term
 !-------------------
  Do k=1,2
   intI(k,1) = I2dsmg1tmg2(s,tmin,tmax,msbot2(k),gammsb(k),mt2,gamtmt)
   intI(k,2) = I2tdsmg1tmg2(s,tmin,tmax,msbot2(k),gammsb(k),mt2,gamtmt)

   stbWchi = stbWchi + (cd(k+1,1) + cd(k+1,2)*s + cd(k+1,3)*s2) * intI(k,1) &
           &         + (cd(k+1,4) + cd(k+1,5)*s + cd(k+1,6)*s2) * intI(k,2)
  End Do
 !--------------------------------------------------
 ! t-Integrals for sbottom_k-sbottom_j terms,
 ! difference if decay-width of sbottom_k is known
 !--------------------------------------------------
 ! k=j :
 !-------
  fakt = kappa(s,mst2,mW2)**2 / mW2

  Do k=1,2 
   If (gammsb(k).Eq.0.0_dp) Then
    intI(k,1) = fakt * (ce(k,1) + ce(k,2)*s) / (s-msbot2(k))**2
   Else
    intI(k,1) = fakt * (ce(k,1) + ce(k,2)*s) / ((s-msbot2(k))**2 + gammsb2(k))
   End If
  End Do
 !-----------
 ! k=1, j=2 :
 !-----------
  If ((gammsb(1).Eq.0._dp).Or.(gammsb(2).Eq.0)) Then
   intI(3,1) = fakt * (ce(3,1) + ce(3,2)*s) / ((s-msbot2(1))*(s-msbot2(2)))
  Else
   intI(3,1) = fakt * (ce(3,1) + ce(3,2)*s)                                 &
           &    * ( (s-msbot2(1)) * (s-msbot2(2)) + gammsb(1)*gammsb(2) )   &
           &    / ( ((s-msbot2(1))**2 + gammsb(1)**2)                       &
           &        * ((s-msbot2(2))**2 + gammsb(2)**2) )
  End If

  stbWchi = stbWchi + 2._dp * (intI(1,1) + intI(2,1) + intI(3,1) ) * ws

  Iname = Iname - 1

 End Function stbWchi

 Subroutine StopDecays3(n_l, id_l, n_nu, id_nu, n_su, n_sd, n_sle, n_snu, n_d &
        & , id_d, id_W, n_n, n_c, Sup, Sdown, Chi0, g, c_UNSu_L, c_UNSu_R     &
        & , c_DNSd_L, c_DNSd_R, c_SdSuW, Sneut, ChiPm, c_CLSn_L, c_CLSn_R     &
        & , c_CDSu_L, c_CDSu_R, c_CNW_L, c_CNW_R, Slept, c_CNuSl_R, prec      &
        & , Check_Real_States)
 !-----------------------------------------------------------------------
 ! calculates the 3-body stop decays
 ! written by Werner Porod, 08.01.02
 ! 22.09.03: including the possiblity to check for real intermediates states
 ! 20.09.10: adapting to new variable types
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n_l, id_l(:), n_nu, id_nu(:), n_su, n_sd, n_sle  &
     & , n_snu, n_d, id_d(:), id_W(:), n_n, n_c
  Real(dp), Intent(in) :: prec, g
  Complex(dp), Intent(in) :: c_CLSn_L(:,:,:), c_CLSn_R(:,:,:), c_CDSu_L(:,:,:) &
     & , c_CDSu_R (:,:,:), c_CNuSl_R(:,:,:), c_UNSu_L(:,:,:), c_UNSu_R(:,:,:)  &
     & , c_CNW_L(:,:), c_CNW_R(:,:), c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)          &
     & , c_SdSuW(:,:)
  Logical, Intent(in) :: Check_Real_States
  Type(particle2), Intent(in) :: Sdown(:)
  Type(particle23), Intent(in) :: Chi0(:), ChiPm(:), Slept(:), Sneut(:)
  Type(particle23), Intent(inout) :: Sup(:)

  Integer :: i1, i2, i3, i, j, k, i_c
  Real(dp) ::  mStop(n_su), mStop2(n_su), mSbottom2(n_sd) &
     & , mN(n_n), mN2(n_n), mSneut2(n_snu), mC(n_c), mC2(n_c), g_C(n_c)  &
     & , mSlept2(n_sle), g_sb(n_sd)
  Real(dp) :: kl2kt2, kl2lt2 &
     & , kll2, ktl2, mbmcha, mtaumc, mbmtau, mbmcha1                          &
     & , mbmcha2, smin, smax, Abslt11Sq, Abskt11Sq, AbslL11Sq, mSl4, f_stop   &
     & , qlqr, ql2qr2, kl, k2l2, lqlkqr, kqllqr,  mbchp1, mchp1chi, mbchi     &
     & , mbt, mb4, mb6, mchi2, mchi4, mchit, mchp1chp2, mst4   &
     & , mb2chi2, mb2st2, mchi2st2, oomw, mbomw, mbtomw, mbchiomw, mchitomw   &
     & , mstomw, mchiomw, kt1l1CQliCQri, kt1Cl1QliCQri, mbchp2, mchp2chi      &
     & , mchp1t, mbchp1omw, mchp1chiomw, ft2, ht2, fht, mSb(2)
  Complex(dp) :: lt11, kt11  &
     & , lt12C, kt12C, lL11C, kL11C, lL12, kL12, mcha1mcha2, kt11C, lt11C   &
     & , l1l2Ql2Qr1, k1l2Ql2Qr1, k2l1Ql2Qr1, k1k2Ql2Qr1, l1l2Ql1Ql2         &
     & , k1l2Ql1Ql2, k2l1Ql1Ql2, k1k2Ql1Ql2, fltQl1, fktQl1, fltQr1, fktQr1 &
     & , hltQl1, hktQl1, hltQr1, hktQr1, hbltQl1, fbltQl1, hbktQl1, fbktQl1 &
     & , ffb1, fhb1, hfb1, hhb1
  Logical :: calc

  Iname = Iname + 1
  NameOfUnit(Iname) = "StopDecays3"

  mStop = Sup%m
  mStop2 = Sup%m2
  mSbottom2 = Sdown%m2
  g_sb = Sdown%g
  mSneut2 = Sneut%m2
  mSlept2 = Slept%m2
  mN = Chi0%m
  mN2 = Chi0%m2
  mC = ChiPm%m
  mC2 = ChiPm%m2
  g_C = ChiPm%g

  Do i1=1,6
   Sup(i1)%gi3 = 0._dp
   Sup(i1)%bi3 = 0._dp
  End Do

  amc2 = mC2
  MSbot2 = mSbottom2(5:6)
  If (Check_Real_States) Then
   GamMSB2 = 0._dp
   gammc = 0._dp
   gamtmt2 = 0._dp
  Else
   GamMSB2 = mSbottom2(5:6) * g_sb(5:6)**2
   gammc = mC * g_C
   gamtmt2 = 2737.37816_dp
  End If
  gammc2 = gammc**2
  gamtmt = Sqrt(gamtmt2)
  GamMSB = Sqrt( GamMSB2 )
!  Write(*,*) "hier",mc,g_c,gammc

  i_c = 1
  If (GenerationMixing) Then
  Else
   f_stop = 2._dp * oo512pi3 / mStop(5)**3
  End If
  !---------------------------------------------
  ! decays into d-quarks, sneutrinos and leptons
  !---------------------------------------------
  istosle = 1 ! needed for integration routine

  If (GenerationMixing) Then
  Else
   Sup(5)%gi3 = 0._dp

   lt11 = c_CDSu_R(1,3,5)
   kt11 = c_CDSu_L(1,3,5)
   lt12C = Conjg(c_CDSu_R(2,3,5))
   kt12C = Conjg(c_CDSu_L(2,3,5) )
   mst2 = mStop2(5)
   mcha1mcha2 = mC(1) * mC(2)
   mb2 = mf_d2(3)
   mbmcha1 = mf_d(3) * mC(1)
   mbmcha2 = mf_d(3) * mC(2)

   Do i1=1,3 ! lepton generation
    mSl2 = mSneut2(i1)
    !-----------------------------------------
    ! checking kinematics 
    !-----------------------------------------
    If (mStop(5).Lt.(Sqrt(mSl2)+mf_d(3)+mf_l(i1))) Cycle

    mbmtau = mf_d(3) * mf_l(i1)
    coeffsl = 0._dp
    Do i2=1,2 ! number of charginos
     kl2kt2 = 0.5_dp * (Abs(c_CLSn_L(i2,i1,i1))**2 * Abs(c_CDSu_L(i2,3,5))**2 &
          &            + Abs(c_CLSn_R(i2,i1,i1))**2 * Abs(c_CDSu_R(i2,3,5))**2)
     kl2lt2 = 0.5_dp * (Abs(c_CLSn_L(i2,i1,i1))**2 * Abs(c_CDSu_R(i2,3,5))**2 &
          &            + Abs(c_CLSn_R(i2,i1,i1))**2 * Abs(c_CDSu_L(i2,3,5))**2)
     kll2 = 2._dp * Real( Conjg(c_CDSu_L(i2,3,5)) * c_CDSu_R(i2,3,5),dp )   &
        & * ( Abs(c_CLSn_L(i2,i1,i1))**2 + Abs(c_CLSn_R(i2,i1,i1))**2 )
     ktl2 = 2._dp * Real(Conjg(c_CLSn_L(i2,i1,i1)) * c_CLSn_R(i2,i1,i1)) & 
        &   * (Abs(c_CDSu_L(i2,3,5))**2 + Abs(c_CDSu_R(i2,3,5))**2)

     mbmcha = mf_d(3) * mC(i2)
     mtaumc = mf_l(i1) * mC(i2)

     CoeffSl(i2,1) = - kll2 * mbmcha - kl2lt2 * mC2(i2) - ktl2 * mtaumc     &
        &  + kl2kt2 * (mstop2(5) + mSneut2(i1) - mf_d2(3)-mf_l2(i1))        &
        &  - 4._dp * Real( Conjg(c_CDSu_L(i2,3,5) * c_CLSn_R(i2,i1,i1))     &
        &                * c_CDSu_R(i2,3,5) * c_CLSn_L(i2,i1,i1),dp ) * mbmtau 

     CoeffSl(i2,2) = kl2lt2 * mC2(i2) * (mf_d2(3)-mstop2(5))  &
        &                  * (mSneut2(i1)-mf_l2(i1))

     CoeffSl(i2,3) = kl2kt2 * (mstop2(5)-mf_d2(3))*(mf_l2(i1)-mSneut2(i1))    &
        &    + kl2lt2 * mC2(i2) * (mstop2(5)+mSneut2(i1)-mf_d2(3)-mf_l2(i1))  &
        &    + ktl2 * mtaumc * (mstop2(5)-mf_d2(3))                           &
        &    + kll2 * mbmcha * (mSneut2(i1)-mf_l2(i1))                        &
        &    - 4._dp * Real( Conjg(c_CDSu_L(i2,3,5) * c_CLSn_L(i2,i1,i1))     &
        &         * c_CDSu_R(i2,3,5) * c_CLSn_R(i2,i1,i1),dp ) *mbmtau *mC2(i2)

     CoeffSl(i2,4) = - kl2kt2
    End Do

    ! this results from a Mathematica output
    lL11C = Conjg(c_CLSn_R(1,i1,i1))
    kL11C = Conjg(c_CLSn_L(1,i1,i1))
    lL12 = c_CLSn_R(2,i1,i1)
    kL12 = c_CLSn_L(2,i1,i1)
    mL2 = mf_l2(i1)

    CoeffSl(3,1) =  - kL11C*kL12*kt11*kt12C*mb2 -                 &
     &  lL11C*lL12*lt11*lt12C*mb2 -                               &
     &  2*kL11C*kL12*kt12C*lt11*mbmcha1 -                         &
     &  2*kt11*lL11C*lL12*lt12C*mbmcha1 -                         &
     &  2*kt12C*lL11C*lL12*lt11*mbmcha2 -                         &
     &  2*kL11C*kL12*kt11*lt12C*mbmcha2 -                         &
     &  kt11*kt12C*lL11C*lL12*mcha1mcha2 -                        &
     &  kL11C*kL12*lt11*lt12C*mcha1mcha2 +                        &
     &  kL11C*kL12*kt11*kt12C*mL2 + lL11C*lL12*lt11*lt12C*mL2 +   &
     &  kL11C*kL12*kt11*kt12C*mSl2 +                              &
     &  lL11C*lL12*lt11*lt12C*mSl2 +                              &
     &  kL11C*kL12*kt11*kt12C*mst2 + lL11C*lL12*lt11*lt12C*mst2


    CoeffSl(3,2) = -(kt11*kt12C*lL11C*lL12 + kL11C*kL12*lt11*lt12C) &
      &    * mcha1mcha2 * (mL2 - mSl2) * (mb2 - mst2)

    CoeffSl(3,3) = kt11*(2*lt12C*(lL11C*lL12*mbmcha1 + kL11C*kL12*mbmcha2)*  &
     &      mSl2 + kt12C*(-(lL11C*lL12*mb2*mcha1mcha2) +                     &
     &        kL11C*kL12*mb2*mL2 - lL11C*lL12*mcha1mcha2*mL2 +               &
     &        kL11C*kL12*mb2*mSl2 +                                          &
     &        lL11C*lL12*mcha1mcha2*mSl2 +                                   &
     &        lL11C*lL12*mcha1mcha2*mst2 -                                   &
     &        kL11C*kL12*mL2*mst2 - kL11C*kL12*mSl2*mst2)) +                 &
     &  lt11*(kL11C*kL12*(-(lt12C*mb2*mcha1mcha2) -                          &
     &        lt12C*mcha1mcha2*mL2 + 2*kt12C*mbmcha1*mSl2 +                  &
     &        lt12C*mcha1mcha2*mSl2 + lt12C*mcha1mcha2*mst2) +               &
     &     lL11C*lL12*(lt12C*mb2*mL2 + lt12C*mb2*mSl2 +                      &
     &        2*kt12C*mbmcha2*mSl2 - lt12C*mL2*mst2 -                        &
     &        lt12C*mSl2*mst2))

    CoeffSl(3,4) = -(kL11C*kL12*kt11*kt12C + lL11C*lL12*lt11*lt12C)

    smin = (Sqrt(msneut2(i1)) + mf_l(i1))**2
    smax = (mstop(5) - mf_d(3))**2

    calc = .True.
    If (Check_Real_States) Then
     If (mStop(5).Gt.(mC(2)+mf_d(3))) Then  ! all intermediate states are real
      calc = .False.           
     Else If (mStop(5).Gt.(mC(1)+mf_d(3))) Then ! only chi^+_2 contributes
      CoeffSl(1,:) = 0._dp
      CoeffSl(3,:) = 0._dp
     End If
    End If

    If (calc) Then
      Sup(5)%gi3(i_c) = f_stop * Dgauss(stoslep,smin,smax,prec)
      Sup(5)%id3(i_c,1) = Sneut(i1)%id
      Sup(5)%id3(i_c,2) = id_l(i1) + 1
      Sup(5)%id3(i_c,3) = id_d(3) 
      i_c = i_c +1
    End If
    
   End Do
  End If

  !---------------------------------------------
  ! decays into d-quarks, neutrinos and sleptons
  !---------------------------------------------
  istosle = 2 ! needed for integration routine

  If (GenerationMixing) Then
  Else
   mb2 = mf_d2(3)
   mst2 = mstop2(5)

   mst2 = mstop2(5)
   mcha1mcha2 = mC(1) * mC(2)
   mb2 = mf_d2(3)
   mbmcha1 = mf_d(3) * mC(1)
   mbmcha2 = mf_d(3) * mC(2)

   Do i1=1,3 ! lepton generation
    Do i2=1,2 ! number of sleptons
     mSl2 = mSlept2(2*(i1-1)+i2)
     !-----------------------------------------
     ! checking kinematics 
     !-----------------------------------------
     If (mStop(5).Lt.(Sqrt(mSl2)+mf_d(3))) Cycle
     mSl4 = mSl2**2

     Do i3=1,2 ! number of charginos
      lt11 = c_CDSu_R(i3,3,5)
      kt11 = c_CDSu_L(i3,3,5)
      lt11C = Conjg(c_CDSu_R(i3,3,5))
      kt11C = Conjg(c_CDSu_L(i3,3,5))
      Abslt11Sq = lt11 * lt11C
      Abskt11Sq = kt11 * kt11C
      mbmcha1 = mf_d(3) * mC(i3)
      AbslL11Sq = Abs(c_CNuSl_R(i3,i1,2*(i1-1)+i2))**2

      CoeffSl(i3,1) = 0.5_dp * Abslt11Sq * mC2(i3)                  &
                  &          * (mst2 + 2._dp * mSl2 - mb2)          & 
                  & + Abskt11Sq*mSl2*(mb2 - mst2 - 0.5_dp * mSl2)   &
                  & + 4._dp*Real(kt11C*lt11,dp)*mbmcha1*mSl2

      CoeffSl(i3,2) = 0.5_dp * Abslt11Sq*mC2(i3)*mSl4*(mst2 - mb2)

      CoeffSl(i3,3) = - 0.5_dp * mSl2*(mSl2*                 &
       &        (2*(kt11C*lt11 + kt11*lt11C)*mbmcha1 +                  &
       &          Abskt11Sq*(mb2 - mst2)) +                             &
       &       Abslt11Sq*mC2(i3)*(-2*mb2 + mSl2 + 2*mst2))

      CoeffSl(i3,4) = - 0.5_dp * (4*Real(kt11C*lt11,dp)*mbmcha1        &
       &       + Abslt11Sq*mC2(i3) + Abskt11Sq*(mb2 - 2*mSl2 - mst2))

      CoeffSl(i3,5) = - 0.5_dp * Abskt11Sq
      coeffSl(i3,:) = AbslL11Sq * coeffSl(i3,:)
     End Do

     lt11 = c_CDSu_R(1,3,5)
     kt11 = c_CDSu_L(1,3,5)
     lt11C = Conjg(c_CDSu_R(1,3,5))
     kt11C = Conjg(c_CDSu_L(1,3,5))
     lt12C = Conjg(c_CDSu_R(2,3,5))
     kt12C = Conjg(c_CDSu_L(2,3,5) )
     Abslt11Sq = lt11 * lt11C
     Abskt11Sq = kt11 * kt11C
     lL11C = Conjg(c_CNuSl_R(1,i1,2*(i1-1)+i2))
     lL12 = c_CNuSl_R(2,i1,2*(i1-1)+i2)

     CoeffSl(3,1) = -(lL11C*lL12*lt11*lt12C*mb2*mcha1mcha2) +  &
       &  2*kt11*kt12C*lL11C*lL12*mb2*mSl2 +                   &
       &  4*kt12C*lL11C*lL12*lt11*mbmcha1*mSl2 +               &
       &  4*kt11*lL11C*lL12*lt12C*mbmcha2*mSl2 +               &
       &  2*lL11C*lL12*lt11*lt12C*mcha1mcha2*mSl2 -            &
       &  kt11*kt12C*lL11C*lL12*mSl4 +                         &
       &  lL11C*lL12*lt11*lt12C*mcha1mcha2*mst2 -              &
       &  2*kt11*kt12C*lL11C*lL12*mSl2*mst2

     CoeffSl(3,2) = lL11C*lL12*lt11*lt12C*mcha1mcha2*mSl4*(mst2 -mb2)

     CoeffSl(3,3) = -(lL11C*lL12*mSl2*                           &
       &    (kt11*mSl2*(2*lt12C*mbmcha2 + kt12C*(mb2 - mst2)) +  &
       &      lt11*(2*kt12C*mbmcha1*mSl2 +                       &
       &         lt12C*mcha1mcha2*(-2*mb2 + mSl2 + 2*mst2))))

     CoeffSl(3,4) = - lL11C*lL12*(lt11*                              &
       &       (2*kt12C*mbmcha1 + lt12C*mcha1mcha2) +               &
       &      kt11*(2*lt12C*mbmcha2 + kt12C*(mb2 - 2*mSl2 - mst2)) )

     CoeffSl(3,5) = - kt11*kt12C*lL11C*lL12

     smin = mSl2
     smax = (mstop(5) - mf_d(3))**2

     calc = .True.
     If (Check_Real_States) Then
      If (mStop(5).Gt.(mC(2)+mf_d(3))) Then  ! all intermediate states are real
       calc = .False.           
      Else If (mStop(5).Gt.(mC(1)+mf_d(3))) Then ! only chi^+_2 contributes
       CoeffSl(1,:) = 0._dp
       CoeffSl(3,:) = 0._dp
      End If
     End If

     If (calc) Then
      Sup(5)%gi3(i_c) = f_stop * Dgauss(stoslep,smin,smax,prec)
      Sup(5)%id3(i_c,1) = Slept(2*(i1-1)+i2)%id + 1
      Sup(5)%id3(i_c,2) = id_nu(i1)
      Sup(5)%id3(i_c,3) = id_d(3) 
      i_c = i_c +1
     End If
    End Do
   End Do
  End If

  !---------------------------
  ! stop_1 -> b + W + chi^0_j
  !---------------------------
  If (GenerationMixing) Then
  Else ! .not.GenerationMixing
   j = 1 ! lightest neutralino
   amN2 = mN2(j)

   If (mStop(5).Gt.(mW + mf_d(3) + Abs(mN(j)) ) ) Then ! kinematics 
    !--------------------------------------------------
    ! sums that are needed for chargino-exchange terms
    !--------------------------------------------------
    sumM(1) = mf_d2(3) + mN2(j) + mW2 + mStop2(5)
    sumM(2) = sumM(1) - mC2(1)
    sumM(3) = sumM(1) - mC2(2)

    !--------------------------------------------------
    ! mass-combinations that are needed in the following
    !--------------------------------------------------
     mbchi = mf_d(3) * mN(j)
     mbt = mf_d(3) * mf_u(3)
     mb2 = mf_d2(3)
     mb4 = mf_d2(3)**2
     mb6 = mb4 * mb2
     mchi2 = mN2(j)
     mchi4 = mchi2 * mchi2
     mchit = mN(j) * mf_u(3)
     mchp1chp2 = mC(1) * mC(2)
     mst2 = mStop2(5)
     mst4 = mst2 * mst2
     mt2 = mf_u2(3)

     mb2chi2 = mf_d2(3) * mchi2
     mb2st2 = mf_d2(3) * mst2
     mchi2st2 = mchi2 * mst2

     oomw = 1._dp / mW2
     mbomw = mf_d2(3) * oomw
     mbtomw = mbt * oomw
     mbchiomw = mbchi * oomw
     mchiomw = mchi2 * oomw
     mchitomw = mchit * oomw
     mstomw = mst2 * oomw
     !--------------------------------------------------
     ! coefficients for chargino-exchange terms
     !--------------------------------------------------
     ! 11,22
     !-------
     Do i=1,2
      qlqr = Real( c_CNW_L(i,j) * Conjg( c_CNW_R(i,j) ), dp )
      ql2qr2 = Abs( c_CNW_L(i,j) )**2 + Abs( c_CNW_R(i,j) )**2
      kl = Real(  c_CDSu_L(i,3,5) * Conjg( c_CDSu_R(i,3,5) ), dp )
      k2l2 = Abs( c_CDSu_L(i,3,5) )**2  + Abs( c_CDSu_R(i,3,5) )**2
      lqlkqr = Abs( c_CDSu_R(i,3,5) )**2 * Abs( c_CNW_R(i,j) )**2    &
           & + Abs( c_CDSu_L(i,3,5) )**2 * Abs( c_CNW_L(i,j) )**2
      kqllqr = Abs( c_CDSu_L(i,3,5) )**2 * Abs( c_CNW_R(i,j) )**2    &
           & + Abs( c_CDSu_R(i,3,5) )**2 * Abs( c_CNW_L(i,j) )**2
      kt1l1CQliCQri = Real( c_CNW_L(i,j) *  c_CDSu_L(i,3,5)             &
            &             * Conjg(  c_CNW_R(i,j) * c_CDSu_R(i,3,5)), dp )
      kt1Cl1QliCQri = Real( c_CNW_R(i,j) *  c_CDSu_L(i,3,5)             &
            &             * Conjg(  c_CNW_L(i,j) * c_CDSu_R(i,3,5)), dp )

      mbchp1 = mf_d(3) * mC(i)
      mchp1chi = mC(i) * mN(j)
 
      ca(i,1) = ( 6._dp*qlqr* k2l2 * mchp1chi * (2._dp * mb2 + mchi2 + mW2) &
            &   + 12._dp * kt1l1CQliCQri * mbchi * (mb2+mchi2+mst2+mW2)     &
            &   + 12._dp * kt1Cl1QliCQri * mbchi * mC2(i)                   &
            &   - 2._dp * kl * ql2qr2 * mbchp1                              &
            &           * ( 3._dp*mb2 + 2._dp*mchi2 + 3._dp*mst2            &
            &             + (mst2+mb2)**2 * oomw)                           &
            &  - lqlkqr * ( (4._dp*mb2+mW2) * (mchi2+mb2) + mchi4           &
            &             + mst2 * (2._dp*mW2 + 4._dp*mb2 + mchi2)          &
            &             + (mb2+mst2) * (mb4+mb2st2+mchi2st2) * oomw )     &
            &  - kqllqr * mC2(i) * ( mchi2 + mst2 + 2._dp * mb2             &
            &                      + (mb4+mb2st2)*oomw)  )

      ca(i,2) = - 12._dp * kt1l1CQliCQri * mbchi                            &
            & - 6._dp * k2l2 * qlqr * mchp1chi                              &
            & + 2._dp * kl* ql2qr2 * mbchp1 * (3._dp+2._dp*(mb2+mst2)*oomw) &
            & + lqlkqr * (6._dp * mb2 + 2._dp * mchi2 + 2._dp * mst2 + mW2  &
            &            + (3._dp * mb4 + mb2chi2 + 4.*mb2st2               &
            &            + 2._dp * mchi2st2 + mst4)*oomw )                  &
            & + kqllqr * mC2(i) * (2._dp + mb2 * oomw)    
                
      ca(i,3) = - 2._dp * kl * ql2qr2 * mbchp1 * oomw         &
            & - lqlkqr * (2._dp + (mchi2+3._dp*mb2+2._dp*mst2) * oomw)

      ca(i,4) = lqlkqr * oomw

      ca(i,5) = -12._dp * kl * qlqr * mbchi                                 &
            & - 6._dp * k2l2 * qlqr * mchp1chi                              &
            & + 2._dp * kl * ql2qr2 * mbchp1 *(3._dp+2._dp*(mb2+mst2)*oomw) &
            & + lqlkqr * ( 6._dp*mb2 + 3._dp*mchi2 + 2._dp*mst2 + 2._dp*mW2 &
            &            + (2._dp*mb4 + 2._dp*mb2st2 + mchi2st2) * oomw)    &
            & + kqllqr * mC2(i) * (1._dp + (2._dp*mb2+mst2) * oomw)

      ca(i,6) = - 4.* kl * ql2qr2 * mbchp1 * oomw                          &
            & - lqlkqr * (4._dp + (mchi2 + 4.*mb2 + 2._dp*mst2 ) * oomw )  &
            & - kqllqr * mC2(i) * oomw

      ca(i,7) = 2._dp * lqlkqr * oomw

      ca(i,8) = -2._dp * kl * ql2qr2 * mbchp1 * oomw                    &
            & - lqlkqr * (2._dp + mb2*oomw) - kqllqr * mC2(i) * oomw

      ca(i,9) = lqlkqr * oomw
     End Do
     !-------
     ! 12
     !-------
     l1l2Ql2Qr1 = lt11 * lt12C * Conjg( c_CNW_R(2,j) ) * c_CNW_L(1,j) &
              & + kt11 * kt12C * c_CNW_R(1,j) * Conjg( c_CNW_L(2,j) )
     k1l2Ql2Qr1 = kt11 * lt12C * Conjg( c_CNW_R(2,j) ) *c_CNW_L(1,j)  &
              & + kt12C * lt11 * c_CNW_R(1,j) * Conjg( c_CNW_L(2,j) )
     k2l1Ql2Qr1 = kt12C * lt11 * Conjg( c_CNW_R(2,j) ) * c_CNW_L(1,j) &
              & + kt11 * lt12C * c_CNW_R(1,j) * Conjg( c_CNW_L(2,j) )
     k1k2Ql2Qr1 = kt11 * kt12C * Conjg( c_CNW_R(2,j) ) * c_CNW_L(1,j) &
              & + lt11 * lt12C * c_CNW_R(1,j) * Conjg( c_CNW_L(2,j) )
     l1l2Ql1Ql2 = lt11 * lt12C * c_CNW_R(1,j) * Conjg( c_CNW_R(2,j) ) &
              & + kt11 * kt12C * c_CNW_L(1,j) * Conjg( c_CNW_L(2,j) )
     k1l2Ql1Ql2 = kt11 * lt12C * c_CNW_R(1,j) * Conjg( c_CNW_R(2,j) ) &
              & + kt12C * lt11 * c_CNW_L(1,j) * Conjg( c_CNW_L(2,j) )
     k2l1Ql1Ql2 = kt12C * lt11 * c_CNW_R(1,j) * Conjg( c_CNW_R(2,j) ) &
              & + kt11 * lt12C * c_CNW_L(1,j) * Conjg( c_CNW_L(2,j) )
     k1k2Ql1Ql2 = kt11 * kt12C * c_CNW_R(1,j) * Conjg( c_CNW_R(2,j) ) &
              & + lt11 * lt12C * c_CNW_L(1,j) * Conjg( c_CNW_L(2,j) )

     mbchp1 = mf_d(3) * mC(1)
     mbchp2 = mf_d(3) * mC(2)
     mchp1chi = mC(1) * mN(j)
     mchp2chi = mC(2) * mN(j)

     ca(3,1) = 3._dp * l1l2Ql2Qr1 * mchp1chi * (2._dp*mb2+mchi2+mW2)   &
           & + 6._dp * k1l2Ql2Qr1 * mbchi * (mb2+mchi2+mst2+mW2)       &
           & + 6._dp * k2l1Ql2Qr1 * mbchi * mchp1chp2                  &
           & + 3._dp * k1k2Ql2Qr1 * mchp2chi * (2._dp*mb2+mchi2+mW2)   &
           & - l1l2Ql1Ql2 * ( (mchi2+mb2)*(mchi2+mst2+mW2+mb2)         &
           &                + mb2 * (3._dp*mb2+2._dp*mchi2+3._dp*mst2) &
           &                + 2._dp*mst2*mW2                           &
           &                + (mb2+mst2)*(mb4+mst2*(mchi2+mb2))*oomw ) &
           & - k1l2Ql1Ql2 * mbchp1 * (3._dp*mb2+2._dp*mchi2+3._dp*mst2 &
           &                         + (mb2+mst2)*(mb2+mst2)*oomw )    &
           & - k2l1Ql1Ql2 * mbchp2 * (3._dp*mb2+2._dp*mchi2+3._dp*mst2 &
           &                         + (mb2+mst2)*(mb2+mst2)*oomw )    &
           & - k1k2Ql1Ql2 * mchp1chp2 * (2._dp*mb2+mchi2+mst2          &
           &                            + (mb2+mst2)*mbomw )

     ca(3,2) = -3._dp * l1l2Ql2Qr1 * mchp1chi                            &
           & -6._dp * k1l2Ql2Qr1 * mbchi                                 &
           & -3._dp * k1k2Ql2Qr1 * mchp2chi                              &
           & + l1l2Ql1Ql2 * ( 6._dp*mb2 + 2._dp*mchi2 + 2._dp*mst2 + mW2 &
           &                + ((2._dp*mb2+mchi2)*(mb2+2._dp*mst2)        &
           &                   +mb4+mst4)*oomw )                         &
           & + k1l2Ql1Ql2 * mbchp1 * (3._dp+2._dp*(mb2+mst2)*oomw)       &
           & + k2l1Ql1Ql2 * mbchp2 * (3._dp+2._dp*(mb2+mst2)*oomw)       &
           & + k1k2Ql1Ql2 * mchp1chp2 * (2._dp+mbomw)

     ca(3,3) = -l1l2Ql1Ql2 * (2._dp+(mchi2+3._dp*mb2+2._dp*mst2)*oomw)  &
           & - ( k1l2Ql1Ql2 * mbchp1 + k2l1Ql1Ql2 * mbchp2 ) * oomw

     ca(3,4) = l1l2Ql1Ql2 * oomw

     ca(3,5) = - 3._dp * l1l2Ql2Qr1 * mchp1chi                            &
           & - 6._dp * k1l2Ql2Qr1 * mbchi                                 &
           & - 3._dp * k1k2Ql2Qr1 * mchp2chi                              &
           & + l1l2Ql1Ql2 * ( 6._dp*mb2+3._dp*mchi2+2._dp*mst2+2._dp*mW2  &
           &                + (2._dp*mb4+mst2*(mchi2+2._dp*mb2))*oomw )   &
           & + k1l2Ql1Ql2 * mbchp1 * (3._dp+2._dp*(mb2+mst2)*oomw )       &
           & + k2l1Ql1Ql2 * mbchp2 * (3._dp+2._dp*(mb2+mst2)*oomw )       &
           & + k1k2Ql1Ql2 * mchp1chp2 * (1._dp+(2._dp*mb2+mst2)*oomw )

     ca(3,6) = - l1l2Ql1Ql2 * (4.+(mchi2+4.*mb2+2._dp*mst2)*oomw)        &
           & - ( 2._dp*k1l2Ql1Ql2 * mbchp1 + 2._dp*k2l1Ql1Ql2 * mbchp2   &
           &   + k1k2Ql1Ql2 * mchp1chp2 ) * oomw

     ca(3,7) = 2._dp*l1l2Ql1Ql2 * oomw

     ca(3,8) = - l1l2Ql1Ql2 * (2._dp+mbomw)                    &
           & - ( k1l2Ql1Ql2 * mbchp1 + k2l1Ql1Ql2 * mbchp2     &
           &   + k1k2Ql1Ql2 * mchp1chp2 ) * oomw 

     ca(3,9) = l1l2Ql1Ql2 * oomw

     ca(3,:) = 2._dp * ca(3,:)  ! factor because of interference

     !--------------------------------------------------
     ! coefficients for top-chargino_i terms
     !--------------------------------------------------
     Do i=1,2
      fltQl1 = Conjg( c_UNSu_R(3,j,5) ) * c_CDSu_R(i,3,5) * c_CNW_R(i,j)
      fktQl1 = Conjg( c_UNSu_R(3,j,5) ) * c_CDSu_L(i,3,5) * c_CNW_R(i,j)
      fltQr1 = Conjg( c_UNSu_R(3,j,5) ) * c_CDSu_R(i,3,5) * c_CNW_L(i,j)
      fktQr1 = Conjg( c_UNSu_R(3,j,5) ) * c_CDSu_L(i,3,5) * c_CNW_L(i,j)
      hltQl1 = Conjg( c_UNSu_L(3,j,5) ) * c_CDSu_R(i,3,5) * c_CNW_R(i,j)
      hktQl1 = Conjg( c_UNSu_L(3,j,5) ) * c_CDSu_L(i,3,5) * c_CNW_R(i,j)
      hltQr1 = Conjg( c_UNSu_L(3,j,5) ) * c_CDSu_R(i,3,5) * c_CNW_L(i,j)
      hktQr1 = Conjg( c_UNSu_L(3,j,5) ) * c_CDSu_L(i,3,5) * c_CNW_L(i,j)

      mbchp1 = mf_d(3) * mC(i)
      mchp1chi = mC(i) * mN(j)
      mchp1t = mC(i) * mf_u(3) 
 
      cb(i,1) = fltQl1 * ( 2._dp * (mchi2-mst2) * (mW2+mchi2+mb2)         &
            &            - mb2 * mst2 * (1._dp+(mchi2-mb2-mst2)*oomw) )   &
            & + 3._dp * fktQl1 * mbchp1 * (mchi2-mst2)                    &
            & + 3._dp * hltQl1 * mchit * (2._dp*mb2+mchi2+mW2)            &
            & + 6._dp * hktQl1 * mchit * mbchp1                           &
            & + fltQr1 * mchp1chi * (2._dp*mW2-mb2*(1._dp+mbomw))         &
            & + fktQr1 * mbchi                                            &
            &    * (3._dp*mW2+mst2*(2._dp-mbomw)-mb2*(2._dp+mbomw)-mchi2) &
            & - hltQr1 * mchp1t                                           &
            &    * (mchi2+mst2*(1._dp+mbomw)+mb2*(2._dp+mbomw))           &
            & - hktQr1 * mbt                                              &
            &    * (2._dp*mchi2 + (mb2+mst2)*(3._dp+(mb2+mst2)*oomw))  

      cb(i,2) = fltQl1 * (mst2-mchi2)*(2._dp-mbomw)               &
             & - 3._dp* hltQl1 * mchit                            &
             & + fktQr1 * mbchi * (mbomw-1._dp)                   &
             & + hltQr1 * mchp1t * (2._dp+mbomw)                  &
             & + hktQr1 * mbt * (3._dp+2._dp*(mb2+mst2)*oomw)

      cb(i,3) = - hktQr1 * mbt * oomw 

      cb(i,4) = fltQl1 * (4.*mb2+2._dp*mW2+mchi2*mst2*oomw       &
            &              +mst2*(2._dp-mbomw) )                 &
            & + 3._dp * fktQl1 * mbchp1                          &
            & - 3._dp * hltQl1 * mchit                           &
            & - fltQr1 * mchp1chi * (1._dp-2._dp*mbomw)          &
            & + fktQr1 * mbchi * (1._dp+2._dp*mbomw+mst2*oomw)   &
            & + hltQr1 * mchp1t * (1._dp+2._dp*mbomw+mst2*oomw)  &
            & + hktQr1 * mbt * (3._dp+2._dp*(mb2+mst2)*oomw)

      cb(i,5) = - fltQl1 * (3._dp+mbomw+(mst2+mchi2)*oomw)    &
            & - ( fktQr1 * mbchi  + hltQr1 * mchp1t           &
            &   + 2._dp * hktQr1 * mbt ) * oomw

      cb(i,6) = fltQl1 * oomw

      cb(i,7) = - 2._dp*fltQl1                                &
            & - ( fltQr1 * mchp1chi  + fktQr1 * mbchi         &
            &   + hltQr1 * mchp1t  + hktQr1 * mbt ) * oomw

      cb(i,8) = fltQl1 * oomw
     End Do

     cb = - sqrt2 * g * cb  ! factor because of interference

     !--------------------------------------------------
     ! coefficients for sbottom_k-chargino_i terms
     !--------------------------------------------------
     Do i=1,2
      mbchp1 = mf_d(3) * mC(i)
      mchp1chi = mC(i) * mN(j)
      mbchp1omw = mbchp1 * oomw
      mchp1chiomw = mchp1chi * oomw

      Do k=1,2
       hbltQl1 = Conjg( c_DNSd_L(3,j,k+4) ) * c_CDSu_R(i,3,5) * c_CNW_R(i,j) &
             & + Conjg( c_DNSd_R(3,j,k+4) ) * c_CDSu_L(i,3,5) * c_CNW_L(i,j)
       fbltQl1 = Conjg( c_DNSd_R(3,j,k+4) ) * c_CDSu_R(i,3,5) * c_CNW_R(i,j) &
             & + Conjg( c_DNSd_L(3,j,k+4) ) * c_CDSu_L(i,3,5) * c_CNW_L(i,j)
       hbktQl1 = Conjg( c_DNSd_L(3,j,k+4) ) * c_CDSu_L(i,3,5) * c_CNW_R(i,j) &
             & + Conjg( c_DNSd_R(3,j,k+4) ) * c_CDSu_R(i,3,5) * c_CNW_L(i,j)
       fbktQl1 = Conjg( c_DNSd_R(3,j,k+4) ) * c_CDSu_L(i,3,5) * c_CNW_R(i,j) &
             & + Conjg( c_DNSd_L(3,j,k+4) ) * c_CDSu_R(i,3,5) * c_CNW_L(i,j)

       cc(k,i,1) = hbltQl1 *mbchi *(mb2-mst2-2._dp*mchi2+mst2*(mbomw+mstomw)) &
              & + fbltQl1 * ( mb2*mW2+mchi2*mst2*(mstomw-1._dp)               &
              &             + mb2*mst2*(mbomw-2._dp)+mb4-2._dp*mb2*mchi2      &
              &             + mb2*mst2*mstomw )                               &
              & + hbktQl1 *mchp1chi *(mb2-mW2-2._dp*mchi2+mst2*(1._dp+mbomw)) &
              & + fbktQl1 * mbchp1 *(mb2-2._dp*mchi2+mst2*(mbomw+mstomw-1._dp))

       cc(k,i,2) = - hbltQl1 * mbchi * (1._dp+mbomw+2._dp*mstomw)            &
              & - fbltQl1 * (mchi2*(1._dp+2._dp*mstomw)+mst2*mstomw          &
              &             +mW2+mb2*(3._dp+mbomw+3._dp*mstomw) )            &
              & - hbktQl1 * mchp1chi * (mbomw-1._dp)                         &
              & - fbktQl1 * mbchp1 * (1._dp+mbomw+2._dp*mstomw) 

       cc(k,i,3) = hbltQl1 * mbchiomw + fbktQl1 * mbchp1omw  &
              & + fbltQl1 * (2._dp+2._dp*mbomw+mchiomw+2._dp*mstomw)

       cc(k,i,4) = - fbltQl1 * oomw

       cc(k,i,5) = ( hbktQl1 * mchp1chi + fbktQl1 * mbchp1            &
               &   + hbltQl1 * mbchi + fbltQl1 * mb2 ) * (1._dp-mstomw)

       cc(k,i,6) = hbltQl1 * mbchiomw + fbltQl1 * (1._dp+mbomw+mstomw) &
               & + hbktQl1 * mchp1chiomw + fbktQl1 * mbchp1omw

       cc(k,i,7) = cc(k,i,4)

       cc(k,i,:) = cc(k,i,:) * Conjg( C_SdSuW(k+4,5) )

      End Do
     End Do
     cc = -2._dp * cc

     !--------------------------------------------------
     ! coefficients for top-top term
     !--------------------------------------------------
     ft2 = Abs( c_UNSu_R(3,j,5) )**2
     ht2 = Abs( c_UNSu_L(3,j,5) )**2
     fht = Real( c_UNSu_R(3,j,5) * Conjg( c_UNSu_L(3,j,5) ),dp )

     cd(1,1) = ft2 * (mchi2-mst2) * (2._dp*mW2-mb2*(1._dp+mbomw))    &
           & - ht2 * mt2 * (mchi2+mst2+mb2*(2._dp+mbomw+mst2*oomw))  &
           & + 2._dp * fht * mchit * (2._dp*mW2-mb2*(1._dp+mbomw))

     cd(1,2) = ht2 * mt2 * (2._dp+mbomw)

     cd(1,3) = ft2 * (mb2+2._dp*mst2+2._dp*mW2+(2._dp*mchi2-mst2)*mbomw)  &
           & + ht2 * mt2 *(1._dp+(2._dp*mb2+mst2)*oomw)                   &
           & - 2._dp * fht * mchit * (1._dp-2._dp*mbomw)

     cd(1,4) = - ft2 * (2._dp+mbomw) - ht2 * mt2 * oomw

     cd(1,5) = - ft2 * (2._dp+mchi2*oomw) &
           & - ( 2._dp*fht * mchit  + ht2 * mt2 ) * oomw

     cd(1,6) = ft2 * oomw

     cd(1,:) = 0.5_dp *g**2 * cd(1,:)

     !--------------------------------------------------
     ! coefficients for top-sbottom_k terms
     !--------------------------------------------------
     Do i=1,2
      ffb1 = Conjg( c_DNSd_R(3,j,i+4) ) * c_UNSu_R(3,j,5)
      fhb1 = Conjg( c_DNSd_R(3,j,i+4) ) * c_UNSu_L(3,j,5)
      hfb1 = Conjg( c_DNSd_L(3,j,i+4) ) * c_UNSu_R(3,j,5)
      hhb1 = Conjg( c_DNSd_L(3,j,i+4) ) * c_UNSu_L(3,j,5)

      cd(1+i,1) = ffb1 * (mb2*(mchi2+mst2)+2._dp*mchi2*(mst2-mchi2-mW2)  &
              &          +(mchi2-mst2)*mb2*mstomw )                      &
              & + (fhb1 * mchit + hfb1 * mbchi)                          &
              &              * (mb2+mst2-2._dp*mchi2-mW2+mb2*mstomw)     &
              & + hhb1 * mbt * (mb2-2._dp*mchi2-mst2+(mb2+mst2)*mstomw)

      cd(1+i,2) = ffb1 * (mchi2-mst2)*(2._dp-mbomw)                     &
              & + (fhb1 * mchit + hfb1 * mbchi) * (1._dp-mbomw)         &
              & - hhb1 * mbt * (1._dp+mbomw+2._dp*mstomw)

      cd(1+i,3) = hhb1 * mbt * oomw

      cd(1+i,4) = ( ffb1*mchi2 + fhb1 * mchit + hfb1 * mbchi  &
              & + hhb1 * mbt) * (1._dp-mstomw)

      cd(1+i,5) = ffb1 * (1._dp+mchiomw+mstomw)                     &
              & + fhb1 * mchitomw + hfb1 * mbchiomw + hhb1 * mbtomw

      cd(1+i,6) = -ffb1 * oomw

      cd(1+i,:) =  sqrt2 * g * Conjg( c_SdSuW(i+4,5) ) * cd(1+i,:)

     End Do

     !--------------------------------------------------
     ! coefficients for sbottom_i-sbottom_k terms
     !--------------------------------------------------
     ! i=k :
     !--------
     Do i=1,2
      ft2 = Abs( c_SdSuW(i+4,5) * c_DNSd_R(3,j,i+4) )**2
      ht2 = Abs( c_SdSuW(i+4,5) * c_DNSd_L(3,j,i+4) )**2
      fht = Real( c_DNSd_R(3,j,i+4) * Conjg( c_DNSd_L(3,j,i+4) ),dp ) &
          & * Abs( c_SdSuW(i+4,5) )**2

      ce(i,1) = -( (ft2 +ht2 ) * (mb2+mchi2) + 4._dp * fht * mbchi )

      ce(i,2) = (ft2 +ht2 )
     End Do
     !-----------
     ! i=1, k=2 :
     !-----------
     ce(3,2) = 2._dp * ( c_DNSd_R(3,j,5) * Conjg( c_DNSd_R(3,j,6) )    &
           &     + c_DNSd_L(3,j,5) * Conjg( c_DNSd_L(3,j,6) ) )  &
           &    * c_SdSuW(5,5) * Conjg( c_SdSuW(6,5) )

     ce(3,1) = - ce(3,2) * (mb2+mchi2)                           &
           & - 4._dp * ( c_DNSd_R(3,j,5) * Conjg( c_DNSd_L(3,j,6) )    &
           &           + c_DNSd_L(3,j,5) * Conjg( c_DNSd_R(3,j,6) ) )  &
           &         * c_SdSuW(5,5) * Conjg( c_SdSuW(6,5) ) * mbchi
 
     smin = (mf_d(3) + mN(j))**2
     smax = (mstop(5) - mW)**2

     calc = .True.
     If (Check_Real_States) Then
      mSb = Sqrt(MSbot2)

      If (mStop(5).Gt.Max(mC(2)+mf_d(3),mf_u(3)+mN(j),mW+mSb(2))) Then
       calc = .False.   ! all intermediate states are real 
      Else ! check now for intermediate real states          
       If (mStop(5).Gt.(mC(2)+mf_d(3))) Then ! no chargino contribution
        ca = 0._dp
        cb = 0._dp
        cc = 0._dp
       Else If (mStop(5).Gt.(mC(1)+mf_d(3))) Then ! only chi^+_2 contributes
        ca(1,:) = 0._dp
        ca(3,:) = 0._dp
        cb(1,:) = 0._dp
        cc(1,:,:) = 0._dp
       End If

       If (mStop(5).Gt.(mN(j)+mf_u(3))) Then ! no top contribution
        cb = 0._dp
        cd = 0._dp
       End If

       If (mStop(5).Gt.(mSb(2)+mW)) Then ! no sbottom contribution
        cc = 0._dp
        cd(2:3,:) = 0._dp
        ce = 0._dp
        
       Else If (mStop(5).Gt.(mSb(1)+mW)) Then ! only sbottom_2 contributes
        cc(:,1,:) = 0._dp
        cd(2,:) = 0._dp
        ce(1,:) = 0._dp
        ce(3,:) = 0._dp
       End If

      End If
     End If

    If (calc)  Then
      Sup(5)%gi3(i_c) = f_stop * dgauss(stbWchi, smin, smax, prec)
      Sup(5)%id3(i_c,1) = Chi0(1)%id
      Sup(5)%id3(i_c,2) = id_W(1)
      Sup(5)%id3(i_c,3) = id_d(3) 
      i_c = i_c +1
    End If
   End If ! kinematic check
  End If
  Sup(5)%g = Sum(Sup(5)%gi2) + Sum(Sup(5)%gi3)
  If (Sup(5)%g.Ne.0._dp) Then
   Sup(5)%bi2 = Sup(5)%gi2 / Sup(5)%g
   Sup(5)%bi3 = Sup(5)%gi3 / Sup(5)%g
  End If

  Iname = Iname - 1

 End Subroutine StopDecays3


 Real(dp) Function  stoslep(s)
 !----------------------------------------------------------------------
 ! Integrand for the decay stop_1 -> b + slepton + lepton'.
 ! written by Werner Porod
 ! 24.01.02:portation to f90
 !----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s
  
  Integer :: i
  Real(dp) :: wsr, s2
  Complex(dp) :: inte(3)
 !-----------------------------------------------------
 ! the following variables from the commons are needed:
 !      integer istosle
 !  -coeffSl,istosle
 !   (coffecients for the integrals in the function stoslep) and
 !   istosle controls the kind of decay:
 !     istosle=1 : stop_1 -> b + sneutrino + electron
 !     istosle=2 : stop_1 -> b + tau-sneutrino + tau
 !     istosle=3 : stop_1 -> b + selectron_L + neutrino
 !     istosle=4 : stop_1 -> b + stau_1 + tau-neutrino
 !     istosle=5 : stop_1 -> b + stau_2 + tau-neutrino
 !  -amc2(i,j), i,j=1,2
 !  -mb2
 !  -amtau2
 !  -msne2,mtsneut2
 !  -mStop2(5)
 !-----------------------------------------------------

  Iname = Iname + 1
  NameOfUnit(Iname) = "StoSle"

  s2 = s*s

  If (istosle.Eq.1) Then ! sneutrino + lepton
   wsr = kappa(s,mst2,mb2) * kappa(s, mSl2, mL2)

  Elseif (istosle.Eq.2) Then
   wsr = kappa(s,mst2,mb2)

  Else
   Write (ErrCan,*) 'Error in function stoslep: istolse =',istosle
   Call TerminateProgram()
  End If

  inte(1) = wsr / ( (s-amc2(1))**2 + gammc2(1) )
  inte(2) = wsr / ( (s-amc2(2))**2 + gammc2(2) )
!  inte(3) = wsr * ( (s-amc2(1))*(s-amc2(2)) + gammc(1)*gammc(2) )         &
!        &  / ( ((s-amc2(1))**2 + gammc2(1)) *((s-amc2(2))**2 + gammc2(2)) )
  inte(3) = wsr / ( (s-amc2(1) + (0._dp,1._dp) * gammc(1) )         &
          &       * (s-amc2(2) - (0._dp,1._dp) * gammc(2) ) )

  stoslep = 0.
  If (istosle.Eq.1) Then
   Do i=1,3
    stoslep = stoslep                                                     &
          &  + Real( inte(i) * ( coeffSl(i,1) + coeffSl(i,2)/s2           &
          &                    + coeffSl(i,3)/s + coeffSl(i,4)*s),dp )
   Enddo
  Else
   Do i=1,3
    stoslep = stoslep                                                       &
          &  + Real( inte(i) * ( coeffSl(i,1) + coeffSl(i,2)/s2             &
          &        + coeffSl(i,3)/s + coeffSl(i,4) * s + coeffSl(i,5)*s2),dp )
   Enddo
  End If
  If (stoslep.Lt.0) Then
   Write (ErrCan,*) 'istosle,s,stoslep'
   Write (ErrCan,*) istosle,s,stoslep
   Write (ErrCan,*) 'inte(1),inte(2),inte(3)'
   Write (ErrCan,*) inte(1),inte(2),inte(3)
   Write (ErrCan,*) 'coeffSl(i,j)'
   Write (ErrCan,*) coeffSl(1,:)
   Write (ErrCan,*) coeffSl(2,:)
   Write (ErrCan,*) coeffSl(2,:)
  End If

  Iname = Iname - 1

 End Function  stoslep

End Module Stop3BodyDecays

