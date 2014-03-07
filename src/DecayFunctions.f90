Module DecayFunctions
! comments
!This module is a collection of subroutines with general formulas for the
!calculation of 2-body decays.

! load modules
Use Control
! load modules

! interfaces
 Interface VectorBosonToTwoFermions
  Module Procedure VectorBosonToTwoFermionsC, VectorBosonToTwoFermionsR
 End Interface

Contains


  Subroutine FermionToFermionPi(mF,mF1,mS,f_pi, G_F, kL,kR,width)
  !-----------------------------------------------------------------------
  ! FermionToFermionPi(mF,mF1,kL,kR,width) calculates the two body decay
  ! width of fermion decaying to another fermion and a pion. mF is the mass
  ! of the decaying fermion, mF1 the mass of the fermion in the
  ! final state and mS the mass of the pion. 
  ! kL and kR are the left and right couplings, respectively. All
  ! couplings are complex.
  ! written by Werner Porod, 12.06.2013
  !-----------------------------------------------------------------------
  implicit none
  Real(dp), Intent(in) :: mF, mF1, mS, f_pi, G_F
  real(dp), intent(out) :: width
  complex(dp), intent(in) :: kL,kR

  Real(dp) :: mFsq,mF1sq,mSsq,kappa

  if ( abs(mF).le.( abs(mF1) + mS ) ) then
   width = 0._dp

  elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
   width = 0._dp

  else
   mFsq = mF * mF
   mF1sq = mF1 * mF1
   mSsq = mS * mS
   kappa = Sqrt( (mFsq-mF1sq-mSsq)**2 - 4._dp * mF1sq*mSsq )
   If (kL.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kR)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else If (kR.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kL)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else
    width = oo4Pi * kappa * (f_pi * G_F)**2                      &
      & * ( Abs(kL+Conjg(kR))**2 * ( (mFsq - mF1sq)**2 - mSsq * (mF-mF1)**2 ) &
      &   + Abs(kL-Conjg(kR))**2 * ( (mFsq - mF1sq)**2 - mSsq * (mF+mF1)**2 ) &
      &   )    / Abs(mF)**3 
   End If

  endif

  End Subroutine FermionToFermionPi

  Subroutine FermionToFermionScalar(mF,mF1,mS,kL,kR,width)
  !-----------------------------------------------------------------------
  ! FermionToFermionScalar(mF,mF1,mS,kL,kR,width) calculates the two body decay
  ! width of fermion decaying to another fermion and a scalar. mF is the mass
  ! of the decaying fermion, mF1 (mS) the mass of the fermion (scalar) in the
  ! final state. kL and kR are the left and right couplings respectively. All
  ! couplings are complex.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
  real(dp), intent(in) :: mF,mF1,mS
  real(dp), intent(out) :: width
  complex(dp), intent(in) :: kL,kR

  real(dp) :: mFsq,mF1sq,mSsq,kappa

  if ( abs(mF).le.( abs(mF1) + mS ) ) then
   width = 0._dp

  elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
   width = 0._dp

  elseif ((mF1.eq.0._dp).and.(mS.eq.0._dp)) then
   If (kL.eq.0._dp) Then
    width = Abs(kR)**2 * mF * oo32Pi
   Else If (kR.eq.0._dp) Then
    width = Abs(kL)**2 * mF * oo32Pi
   Else
    width = (Abs(kL)**2 + Abs(kR)**2) * mF * oo32Pi
   End If

  elseif (mF1.eq.0._dp) then
   mFsq = mF * mF
   mSsq = mS * mS
   If (kL.eq.0._dp) Then
    width =  oo32pi * Abs(kR)**2 * (mFsq - mSsq)**2 /Abs(mF)**3 
   Else If (kR.eq.0._dp) Then
    width =  oo32pi * Abs(kL)**2 * (mFsq - mSsq)**2 /Abs(mF)**3 
   Else
    width =  oo32pi * (Abs(kL)**2 + Abs(kR)**2) * (mFsq - mSsq)**2 /Abs(mF)**3 
   End If

  elseif (mS.eq.0._dp) then
   mFsq = mF * mF
   mF1sq = mF1 * mF1
   If (kL.eq.0._dp) Then
    width = oo32Pi * (mFsq - mF1sq) * Abs(kR)**2 * (mFsq + mF1sq) / Abs(mF)**3
   Else If (kR.eq.0._dp) Then
    width = oo32Pi * (mFsq - mF1sq) * Abs(kL)**2 * (mFsq + mF1sq) / Abs(mF)**3
   Else
    width = oo32Pi * (mFsq - mF1sq)                         &
         & * ( (Abs(kL)**2 + Abs(kR)**2) * (mFsq + mF1sq)  &
         &   + 4._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 ) / Abs(mF)**3  
   End If

  else
   mFsq = mF * mF
   mF1sq = mF1 * mF1
   mSsq = mS * mS
   kappa = Sqrt( (mFsq-mF1sq-mSsq)**2 - 4._dp * mF1sq*mSsq )
   If (kL.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kR)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else If (kR.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kL)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else
    width = oo32Pi * kappa                                       &
         & * ( (Abs(kL)**2 + Abs(kR)**2) * (mFsq + mF1sq - mSsq) &
         &   + 4._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 )    / Abs(mF)**3  
   End If

  endif

  End Subroutine FermionToFermionScalar


  Subroutine FermionToFermionVectorBoson(mF,mF1,mV,kL,kR,width)
  !-----------------------------------------------------------------------
  ! FermionToFermionVectorBoson calculates the two body decay 
  ! width of fermion decaying to another fermion and a vectorboson. mF is the
  ! mass of the decaying fermion, mF1 (mV) the mass of the fermion 
  ! (vector boson) in the final state. kL and kR are the left and right 
  ! couplings respectively. All couplings are complex.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
   real(dp), intent(in) :: mF,mF1,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: kL,kR

   real(dp) :: mFsq,mF1sq,mVsq,kappa

   if ( abs(mF).le.( abs(mF1) + mV ) ) then
    width = 0._dp

   elseif (mV.eq.0._dp) then
    write(ErrCan,*) 'Warning from subroutine FermionToFermionVectorBoson.'
    write(ErrCan,*) ' Vectorboson mass = 0 has occured, setting width to 0!!!' 
    width = 0._dp

   elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
    width = 0._dp

   elseif (mF1.eq.0._dp) then
    mFsq = mF * mF
    mVsq = mV * mV
    width =  oo16pi * (Abs(kL)**2 + Abs(kR)**2) * (mFsq - mVsq)         &
          & * ( 0.5_dp * ( mFsq**2 /  mVsq + mFsq ) - mVsq ) / Abs(mF)**3 

   else
    mFsq = mF * mF
    mF1sq = mF1 * mF1
    mVsq = mV * mV
    kappa = Sqrt( (mFsq-mF1sq-mVsq)**2 - 4._dp * mF1sq*mVsq )
    width = oo16Pi * kappa                                                   &
        & * ( (Abs(kL)**2 + Abs(kR)**2)                                      &
        &     * ( 0.5_dp * ( (mFsq-mF1sq)**2 / mVsq + mFsq + mF1sq) - mVsq )  &
        &   - 6._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 )  / Abs(mF)**3  

   endif

  End Subroutine FermionToFermionVectorBoson


  Subroutine ScalarToScalarVectorBoson(mS,mS1,mV,coup,width)
  !-----------------------------------------------------------------------
  ! ScalarToScalarVectorBoson calculates the two body decay 
  ! width of a scalar decaying to another scalar and a vectorboson. mS is the
  ! mass of the decaying scalar, mS1 (mV) the mass of the scalar (vectorboson)
  ! in the final state. coup is a complex coupling.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
   real(dp), intent(in) :: mS,mS1,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: coup

   real(dp) :: mSsq,mS1sq,mVsq,kappa

   if ( abs(mS).le.( abs(mS1) + mV ) ) then
    width = 0._dp

   elseif (mV.eq.0._dp) then
    write(ErrCan,*) 'Warning from subroutine ScalarToScalarVectorBoson.'
    write(ErrCan,*) 'Vectorboson mass = 0 has occured, setting width to 0!!!' 
    width = 0._dp

  elseif (Abs(coup).eq.0._dp) then
   width = 0._dp

   elseif (mS1.eq.0._dp) then
    mSsq = mS * mS
    mVsq = mV * mV
    width = oo16pi * Abs(coup)**2  * (mSsq - mVsq)**3 / ( mS**3 * mVsq ) 

   else
    mSsq = mS * mS
    mS1sq = mS1 * mS1
    mVsq = mV * mV
    kappa = Sqrt( (mSsq-mS1sq-mVsq)**2 - 4._dp * mS1sq*mVsq )
    width = oo16Pi * Abs(coup)**2 * kappa**3 / ( mS**3 * mVsq ) 

   endif

  End Subroutine ScalarToScalarVectorBoson


 Subroutine ScalarToTwoFermions(mS,mF1,mF2,kL,kR,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 fermions. mS is the mass of the decaying scalar and
 ! mF1 (mF2) are the masses of the fermions.
 ! kL and kR are the left and right couplings respectively. All
 ! couplings can be complex.
 ! written by Werner Porod, 3.11.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mF1,mF2
  real(dp), intent(out) ::width
  complex(dp), intent(in) :: kL,kR

  real(dp) :: mSsq,mF1sq,mF2sq,kappa

  if ( abs(mS).le.( abs(mF1) + abs(mF2) ) ) then
   width = 0._dp

  elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
   width = 0._dp

  elseif ((mF1.eq.0._dp).and.(mF2.eq.0._dp)) then
   width = (Abs(kL)**2 + Abs(kR)**2) * mS * oo16Pi

  elseif ((mF1.eq.0._dp).or.(mF2.eq.0._dp)) then
   mSsq = mS * mS
   if (mF1.eq.0._dp) mF1sq = mF2 * mF2
   if (mF2.eq.0._dp) mF1sq = mF1 * mF1
   width = (Abs(kL)**2 + Abs(kR)**2) * (mSsq - mF1sq)**2 * oo16pi / mS**3 

  elseif (mF1.eq.mF2) then
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   kappa = Sqrt( mSsq**2 - 4._dp * mSsq * mF1sq )
   width = oo16Pi * kappa                                        &
     &   *  ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - 2._dp * mF1sq)  &
     &      - 4._dp * Real( Conjg(kL) * kR,dp ) * mF1sq )  / mS**3  

  else
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   mF2sq = mF2 * mF2
   kappa = Sqrt( (mSsq-mF1sq-mF2sq)**2 - 4._dp * mF1sq*mF2sq )
   width = oo16Pi * kappa                                            &
     &   * ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - mF1sq - mF2sq)      &
     &     - 4._dp * Real( Conjg(kL) * kR,dp ) * mF2 * mF1 ) / mS**3  

  endif

 End Subroutine ScalarToTwoFermions

 Subroutine ScalarToFermionGravitino(mS,mF1,m32,F,kL,kR,width)
 !-----------------------------------------------------------------------
 ! using the formula of 
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 fermions. mS is the mass of the decaying scalar and
 ! mF1 (m32) are the masses of the fermion (gravitino).
 ! F is the SUSY breaking F-term.
 ! kL and kR are the left and right couplings respectively. All
 ! couplings can be complex.
 ! written by Werner Porod, 3.11.1999
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mS,mF1,m32,F
  Real(dp), Intent(out) ::width
  Complex(dp), Intent(in) :: kL,kR

  Real(dp) :: mSsq,mF1sq,m32sq,kappa

  width = 0._dp
  If ( Abs(mS).Le.( Abs(mF1) + Abs(m32) ) ) Then
   Return

  Elseif ((Abs(kL).Eq.0._dp).And.(Abs(kR).Eq.0._dp)) Then
   Return

  Elseif (mF1.Eq.0._dp) Then
   mSsq = mS * mS
   m32sq = m32 * m32
   width = (Abs(kL)**2 + Abs(kR)**2) * (mSsq - m32sq)**4 / mSsq**2

  Elseif (mF1.Eq.m32) Then
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   kappa = Sqrt( mSsq**2 - 4._dp * mSsq * mF1sq )
   width = kappa**3                                               &
     &   *  ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - 2._dp * mF1sq)  &
     &      - 4._dp * Real( Conjg(kL) * kR,dp ) * mF1sq )  / mSsq**2

  Else
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   m32sq = m32 * m32
   kappa = Sqrt( (mSsq-mF1sq-m32sq)**2 - 4._dp * mF1sq*m32sq )
   width = kappa**3                                                  &
     &   * ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - mF1sq - m32sq)      &
     &     - 4._dp * Real( Conjg(kL) * kR,dp ) * m32 * mF1 ) / mSsq**2

  Endif

  width = oo16Pi * mS * width / F**2

 End Subroutine ScalarToFermionGravitino


 Subroutine ScalarToTwoScalars(mS,mS1,mS2,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 scalars. mS is the mass of the decaying scalar and
 ! mS1 (mS2) are the masses of the scalars in the final state.
 ! coup is a complex coupling.
 ! written by Werner Porod, 3.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mS1,mS2
  real(dp), intent(out) ::width
  complex(dp), intent(in) :: coup

  real(dp) :: mSsq,mS1sq,mS2sq,kappa

  if ( abs(mS).le.( abs(mS1) + abs(mS2) ) ) then
   width = 0._dp

  elseif (Abs(coup).eq.0._dp) then
   width = 0._dp

  elseif ((mS1.eq.0._dp).and.(mS2.eq.0._dp)) then
   width = Abs(coup)**2 * oo16Pi / mS 

  elseif ((mS1.eq.0._dp).or.(mS2.eq.0._dp)) then
   mSsq = mS * mS
   if (mS1.eq.0._dp) mS1sq = mS2 * mS2
   if (mS2.eq.0._dp) mS1sq = mS1 * mS1
   width = Abs(coup)**2 * (mSsq - mS1sq) * oo16pi / mS**3 

  else
   mSsq = mS * mS
   mS1sq = mS1 * mS1
   mS2sq = mS2 * mS2
   kappa = Sqrt( (mSsq-mS1sq-mS2sq)**2 - 4._dp * mS1sq*mS2sq )
   width = oo16Pi * Abs(coup)**2 * kappa / mS**3 

  endif

  End Subroutine ScalarToTwoScalars
Subroutine ScalarToTwoVectorbosons(mS,mV1,mV2,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 Vectorbosons. mS is the mass of the decaying scalar and
 ! mV are the mass of the vectorboson, coup is a real coupling.
 ! written by Werner Porod, 5.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mS,mV1,mV2
   Real(dp), Intent(out) :: width
   Complex(dp), Intent(in) :: coup

   Real(dp) :: mSsq,mV1sq, mV2sq,kappa, x

   If ( Abs(mS).Le.( Abs(mV1)+Abs(mV2) ) ) Then
    width = 0._dp

   Elseif (mV1.Eq.0._dp) Then
    Write(ErrCan,*) 'Warning, in subroutine ScalarToTwoVectorbosons'
    Write(ErrCan,*) 'm_V1 = 0, setting width to 0'
    width = 0._dp

   Elseif (mV2.Eq.0._dp) Then
    Write(ErrCan,*) 'Warning, in subroutine ScalarToTwoVectorbosons'
    Write(ErrCan,*) 'm_V2 = 0, setting width to 0'
    width = 0._dp


   Elseif (Abs(coup).Eq.0._dp) Then
    width = 0._dp

   Elseif (mV1.Eq.mV2) Then
    mSsq = mS * mS
    mV1sq = mV1**2
    x = mSsq / mV1sq
    kappa = Sqrt( 1._dp - 4._dp / x )
    width = oo64Pi * Abs(coup)**2 * kappa * (x*x - 4._dp*x + 12._dp) / mS 

   Else
    mSsq = mS * mS
    mV1sq = mV1**2
    mV2sq = mV2**2
    kappa = Sqrt( (mSsq-mV1sq-mV2sq)**2 - 4._dp * mV1sq*mV2sq )/(mS**3)
    width = Abs(coup)**2/(mV1sq*mV2sq)*(mV1sq**2 + 10._dp*mV1sq*mV2sq - &
             & 2._dp*(mV1sq+mV2sq)*mSsq + mV2sq**2 +mSsq**2)
    width = oo64Pi*width*kappa

   Endif

  End Subroutine ScalarToTwoVectorbosons


 Subroutine ScalarToVectorbosonsVR(mS,mV,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 Vectorbosons, where on of them is real and the other is
 ! virtuell. mS is the mass of the decaying scalar and
 ! mV are the mass of the vectorboson. coup is an effective real coupling,
 ! which is the scalar-vectorboson coupling times g times a correction
 ! factor for the Z-boson. This correction factor squared reads:
 ! (m_Z/m_W)**2 * (3.5-20.*sinw2/3.+80.*sinw4/9.) * cosW2
 ! written by Werner Porod, 5.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mV,coup
  real(dp), intent(out) :: width

  real(dp) :: RHiggs,x,x2,sum1,sum2,sum3

  width = 0._dp

  if ((mS.lt.mV).or.(Abs(coup).eq.0._dp) ) return
  If (mS.Gt.(2._dp*mV-0.1_dp)) return ! in this case a refined calculation
                                      ! is necessary

  x = (mV / mS )**2
  x2 = x*x
  sum1 = 3._dp * (1._dp-8._dp*x+20._dp*x2) / sqrt(4._dp*x-1._dp)   &
     & * acos( (3._dp*x-1._dp) / (2._dp*sqrt(x**3)) )
  sum2 = -0.5_dp * (1._dp-x) * (2._dp-13._dp*x+47._dp*x2) / x
  sum3 = - 1.5_dp * (1._dp-6._dp*x+4._dp*x2) * log(x)

  Rhiggs = sum1 + sum2 + sum3

  width = 3._dp * oo128pi3 * coup**2 * mS * Rhiggs 

 End Subroutine ScalarToVectorbosonsVR


  Subroutine VectorBosonToTwoFermionsC(mV,mF1,mF2,Nc,cL,cR,width)
  !-----------------------------------------------------------------------
  ! VectorBosonToTwoScalars calculates the two body decay 
  ! width of a vector poson decaying to two fermionss. mV is the
  ! mass of the decaying vector boson, mF2 and mF2 are the masses of the fermions
  ! in the final state. cL and cR are complex couplings and Nc a color factor
  ! written by Werner Porod, 03.09.2010
  !-----------------------------------------------------------------------
  Implicit None
   Integer :: Nc ! color factor
   Real(dp), Intent(in) :: mF2,mF1,mV
   Real(dp), Intent(out) :: width
   Complex(dp), Intent(in) :: cL, cR

   Real(dp) :: mF1sq, mF2sq, mVsq, kappa, r

   If (     ( mV.Le.( Abs(mF1) + Abs(mF2) ) ) &
      & .Or.((Abs(cL)+Abs(cR)).Eq.0._dp) ) Then
    width = 0._dp
    Return

   Else If ((mF1.Eq.0._dp).And.(mF2.Eq.0._dp)) Then
    width = 2._dp * (Abs(cL)**2 + Abs(cR)**2)

   Else If (mF1.Eq.0._dp) Then
    r = (mF2 / mV)**2 
    width = (Abs(cL)**2 + Abs(cR)**2) * (1._dp-r) * (2._dp - r**2 - r)

   Else If (mF2.Eq.0._dp) Then
    r = (mF1 / mV)**2 
    width = (Abs(cL)**2 + Abs(cR)**2) * (1._dp-r) * (2._dp - r**2 - r)

   Else If (mF1.Eq.mF2) Then
    r = (mF1 / mV)**2 
    width = ( (Abs(cL)**2 + Abs(cR)**2) * (1._dp-r)   &
        &   + 6._dp * Real(cL*Conjg(cR),dp) * r )            &
        &   * 2._dp * Sqrt(1._dp - 4._dp * r)

   Else
    mF1sq = mF1 * mF1
    mF2sq = mF2 * mF2
    mVsq = mV * mV
    kappa = (mVsq-mF1sq-mF2sq)**2 - 4._dp * mF1sq*mF2sq 
    width = (Abs(cL)**2 + Abs(cR)**2)                                    &
        &   * (2._dp * mVsq - mF1sq - mF2sq - (mF1sq-mF2sq)**2 / mVsq )  &
        & + 12._dp * Real(cL*Conjg(cR),dp) * mF1 * mF2 
    width = width * Sqrt(kappa) / mVsq**2

   End If

   width = Nc * oo48Pi * width * mV

  End Subroutine VectorBosonToTwoFermionsC


  Subroutine VectorBosonToTwoFermionsR(mV,mF1,mF2,Nc,cL,cR,width)
  !-----------------------------------------------------------------------
  ! VectorBosonToTwoScalars calculates the two body decay 
  ! width of a vector poson decaying to two fermionss. mV is the
  ! mass of the decaying vector boson, mF2 and mF2 are the masses of the fermions
  ! in the final state. cL and cR are complex couplings and Nc a color factor
  ! written by Werner Porod, 03.09.2010
  !-----------------------------------------------------------------------
  Implicit None
   Integer :: Nc ! color factor
   Real(dp), Intent(in) :: mF2,mF1,mV, cL, cR
   Real(dp), Intent(out) :: width

   Real(dp) :: mF1sq, mF2sq, mVsq, kappa, r

   If (( mV.Le.( mF1 + mF2 ) ).Or.((Abs(cL)+Abs(cR)).Eq.0._dp) ) Then
    width = 0._dp
    Return

   Else If ((mF1.Eq.0._dp).And.(mF2.Eq.0._dp)) Then
    width = 2._dp * (cL**2 + cR**2)

   Elseif (mF1.Eq.0._dp) Then
    r = (mF2 / mV)**2 
    width = (cL**2 + cR**2) * (1._dp-r) * (2._dp - r**2 - r)

   Elseif (mF2.Eq.0._dp) Then
    r = (mF1 / mV)**2 
    width = (cL**2 + cR**2) * (1._dp-r) * (2._dp - r**2 - r)

   Elseif (mF1.Eq.mF2) Then
    r = (mF1 / mV)**2 
    width = ( (cL**2 + cR**2) * (1._dp-r) + 6._dp * cL * cR * r )   &
        &   * 2._dp * Sqrt(1._dp - 4._dp * r)

   Else
    mF1sq = mF1 * mF1
    mF2sq = mF2 * mF2
    mVsq = mV * mV
    kappa = (mVsq-mF1sq-mF2sq)**2 - 4._dp * mF1sq*mF2sq 
    width = (cL**2 + cR**2)                                              &
        &   * (2._dp * mVsq - mF1sq - mF2sq - (mF1sq-mF2sq)**2 / mVsq )  &
        & + 12._dp * cL * cR  * mF1 * mF2
    width = width * Sqrt(kappa) / mVsq**2

   Endif

   width = Nc * oo48Pi * width * mV

  End Subroutine VectorBosonToTwoFermionsR

  Subroutine VectorBosonToTwoScalars(mV,mS1,mS2,Nc,coup,width)
  !-----------------------------------------------------------------------
  ! ScalarToScalarVectorBoson calculates the two body decay 
  ! width of a vector poson decaying to two scalars. mV is the
  ! mass of the decaying vector boson, mS and mS1 are the masses of the scalars
  ! in the final state. coup is a complex coupling and Nc a color factor
  ! written by Werner Porod, 03.09.2010
  !-----------------------------------------------------------------------
  Implicit None
   Integer :: Nc ! color factor
   Real(dp), Intent(in) :: mS2,mS1,mV
   Real(dp), Intent(out) :: width
   Complex(dp), Intent(in) :: coup

   Real(dp) :: r, r2

   If (( mV.Le.( mS1 + mS2 ) ).Or.(Abs(coup).Eq.0._dp) ) Then
    width = 0._dp
    Return

   Else If (mS1.Eq.0._dp) Then
    r = (mS2/mV)**2    
    width = (1._dp - r )**3

   Else If (mS2.Eq.0._dp) Then
    r = (mS1/mV)**2    
    width = (1._dp - r )**3

   Else If (mS1.Eq.mS2) Then
    r = (mS1/mV)**2    
    width = (1._dp - 4._dp * r )**1.5_dp

   Else
    r = (mS1/mV)**2    
    r2 = (mS2/mV)**2
    width = ( 1._dp - 2._dp*(r+r2) + (r-r2)**2 )**1.5_dp    

   End If

   width = Nc * oo48Pi * Abs(coup)**2  * width * mV

  End Subroutine VectorBosonToTwoScalars

 Subroutine VectorBosonToTwoVectorBosons(mV,mV1,mV2,c,width)
  !-----------------------------------------------------------------------
  ! Two Body decay of a vector boson into two other vector boson
  ! Florian Staub. 25.08.2011
  !-----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mV2,mV1,mV
   Real(dp), Intent(out) :: width
   Complex(dp) :: c

   Real(dp) :: mV1sq, mV2sq, mVsq, kappa

   If (( mV.Le.( mV1 + mV2 ) ).Or.(c.Eq.0._dp) ) Then
    width = 0._dp
    Return

   Else If ((mV1.Eq.0._dp).Or.(mV2.Eq.0._dp)) Then
    width = 0._dp

   Else
    mV1sq = mV1 * mV1
    mV2sq = mV2 * mV2
    mVsq = mV * mV
    kappa = (mVsq-mV1sq-mV2sq)**2 - 4._dp * mV1sq*mV2sq 
    width = ((mV1-mV2-mV)*(mV1+mV2-mV)*(mV1-mV2+mV)*(mV1+mV2+mV)     &
          &   * (mV1sq*mV1sq + 10._dp*mV1sq*mV2sq + mV2sq*mV2sq      &
          &     + 10._dp*mV1sq*mVsq + 10._dp*mV2sq*mVsq + mVsq*mVsq) &
          & ) / (4._dp*mV1sq*mV2sq*mVsq)
    width = width * Sqrt(kappa) / mV**3

   Endif

   width = Abs(c)**2 *oo48Pi * width

 End Subroutine VectorBosonToTwoVectorBosons

 Subroutine VectorBosonToScalarAndVectorBoson(mV,mS,mV2,c,width)
  !-----------------------------------------------------------------------
  ! Two Body decay of a vector boson into a scalar and a vector boson
  ! Florian Staub. 25.08.2011
  !-----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: mV2,mS,mV
   Real(dp), Intent(out) :: width
   Complex(dp) :: c

   Real(dp) :: mSsq, mV2sq, mVsq, kappa

   If (( mV.Le.( mS + mV2 ) ).Or.(c.Eq.0._dp) ) Then
    width = 0._dp
    Return

   Else If (mV2.Eq.0._dp) Then
    width = 0._dp

   Else
    mSsq = mS * mS
    mV2sq = mV2 * mV2
    mVsq = mV * mV
    kappa = (mVsq-mSsq-mV2sq)**2 - 4._dp * mSsq*mV2sq 
    width = 2._dp + (mSsq-mVsq-mV2sq)**2/(4._dp*mVsq*mV2sq)
    width = width * Sqrt(kappa) / mV**3

   Endif

   width = Abs(c)**2 * oo48Pi * width

  End Subroutine VectorBosonToScalarAndVectorBoson

End Module DecayFunctions

