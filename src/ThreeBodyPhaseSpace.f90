Module ThreeBodyPhaseSpace
! comments
!This module contains general subroutines that are needed for the three-body
!decays of fermions.

! load modules
Use Control
Use Mathematics
! load modules

! private variables
! for integration
 Real(dp), Private :: mS2,mSG2,mf2(4)
 Real(dp), Private :: mS12,mS1gS1,mS22,mS2gS2
 Real(dp), Private :: mSgS,mSgS2,mT2,mTgT,mTgT2
 Real(dp), Private :: mG2,mGgG2
 Integer, Private :: int_v
! private variables

Contains



 Subroutine F3BDgaugeSscalarSint(gauge,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the matrix element squared of
 ! M^2_fi = - Tr[\[Gamma]_\[Mu].(C[1] Pl + C[2] Pr).(P1+mf[1]).
 !           (C[3] Pl + C[4] Pr).(P2+mf[2])]
 !      * Tr[\[Gamma]^\[Mu].(C[5] Pl + C[6] Pr).(P3+mf[3]).
 !           (C[7] Pl + C[8] Pr).(P4+mf[4])] 
 !      / ( (p1-p2)^2 - m^2_G  + I m^2_G \[CapitalGamma]_G )
 !         * (p1-p2)^2 - m^2_S  - I m^2_S \[CapitalGamma]_S ) )
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p2)^2, t=(p1-p3)^2, and u=(p1-p4)^2.
 ! Here {C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8]} 
 ! are the couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses,
 ! gauge = {m_G ,\[CapitalGamma]_G, m_S,\[CapitalGamma]_S} 
 ! are the masses and the total decay widths of the gauge boson and the scalar
 ! boson, respectively.
 ! written by Werner Porod, 9.1.2000
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: gauge(4),mf(4),eps
   Complex(dp), Intent(in) :: coup(8)
   Complex(dp), Intent(inout) :: int1(4)
   Complex(dp), Intent(out) :: erg
   Logical, Intent(inout) :: Integrate

  Integer :: i1,Imin,Imax,len1
  Real(dp) :: smin,smax,sG(4),int1a(8),int1b(4),int1c(2)
  Complex(dp) :: sumI(4)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   erg = (0._dp,0._dp)
   Return
  Endif

  If ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
   int1 = (0._dp,0._dp)
   Integrate = .False.
   erg = (0._dp,0._dp)
   Return
  Endif

  If (      ((coup(1) * coup(3) + coup(2) * coup(4)).Eq.(0._dp,0._dp) ) &
     & .And.((coup(2) * coup(3) + coup(1) * coup(4)).Eq.(0._dp,0._dp) )) Then
   If (Integrate) Integrate = .False.
   erg = (0._dp,0._dp)
   Return

  Elseif (     ((coup(5) * coup(7) + coup(6) * coup(8)).Eq.(0._dp,0._dp) ) &
       & .And.((coup(6) * coup(7) + coup(5) * coup(8)).Eq.(0._dp,0._dp) ) ) Then
   If (Integrate) Integrate = .False.
   erg = (0._dp,0._dp)
   Return
  Endif

  sumI = 0
  mf2 = mf**2

  If (Integrate) Then
   int1 = ZeroC
   smax = (Abs(mf(1)) - Abs(mf(2)))**2
   smin = (Abs(mf(3)) + Abs(mf(4)))**2
   mS2 = gauge(1)**2
   mSgS = gauge(1)*gauge(2)
   mSgS2 = mSgS**2
   mT2 = gauge(3)**2
   mTgT = gauge(3)*gauge(4)
   mTgT2 = mTgT**2
   If (gauge(1).Lt.gauge(3)) Then
    sG(1) = (gauge(1)-2._dp*gauge(2))**2
    sG(2) = (gauge(1)+2._dp*gauge(2))**2
    sG(3) = (gauge(3)-2._dp*gauge(4))**2
    sG(4) = (gauge(3)+2._dp*gauge(4))**2
   Else
    sG(1) = (gauge(3)-2._dp*gauge(4))**2
    sG(2) = (gauge(3)+2._dp*gauge(4))**2
    sG(3) = (gauge(1)-2._dp*gauge(2))**2
    sG(4) = (gauge(1)+2._dp*gauge(2))**2
   Endif
   If (sG(2).Ge.sG(3)) Then
    sG(1) = Min( sG(1),sG(3) )
    sG(2) = Max( sG(2),sG(4) )
    len1 = 3
   Else
    len1 = 5
   Endif
   Imin = len1
   Imax = 0
   Do i1=1,len1-1
    If (smin.Lt.sG(len1-i1)) Imin = len1 - i1
    If (smax.Gt.sG(i1)) Imax = i1
   Enddo
  
   If ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp)) Then
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel6,2,smin,smax,int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel6,2,smin,sG(Imin),int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
     Call DgaussInt(F3BDgaugeSscalarSkernel6,2,sG(Imin),smax,int1c,eps)
     int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel6,2,smin,sG(Imin),int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel6,2,sG(i1),sG(i1+1),int1c,eps)
      int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel6,2,sG(Imax),smax,int1c,eps)
     int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel5,2,smin,smax,int1c,eps)
     int1(2) = int1c(1) + Ic * int1c(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel5,2,smin,sG(Imin),int1c,eps)
     int1(2) = int1c(1) + Ic * int1c(2)
     Call DgaussInt(F3BDgaugeSscalarSkernel5,2,sG(Imin),smax,int1c,eps)
     int1(2) = int1(2) + int1c(1) + Ic * int1c(2)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel5,2,smin,sG(Imin),int1c,eps)
     int1(2) = int1c(1) + Ic * int1c(2)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel5,2,sG(i1),sG(i1+1),int1c,eps)
      int1(2) = int1(2) + int1c(1) + Ic * int1c(2)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel5,2,sG(Imax),smax,int1c,eps)
     int1(2) = int1(2) + int1c(1) + Ic * int1c(2)
    Endif

   Elseif (mf(2).Eq.0._dp) Then
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel4,4,smin,smax,int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel4,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Call DgaussInt(F3BDgaugeSscalarSkernel4,4,sG(Imin),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel4,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel4,4,sG(i1),sG(i1+1),int1b,eps)
      int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
      int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel4,4,sG(Imax),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Endif

   Elseif (mf(3).Eq.0._dp) Then
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel3,4,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel3,4,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Call DgaussInt(F3BDgaugeSscalarSkernel3,4,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel3,4,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel3,4,sG(i1),sG(i1+1),int1b,eps)
      int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
      int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel3,4,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Endif

   Elseif (mf(4).Eq.0._dp) Then
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel2,4,smin,smax,int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel2,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     Call DgaussInt(F3BDgaugeSscalarSkernel2,4,sG(Imin),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel2,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel2,4,sG(i1),sG(i1+1),int1b,eps)
      int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
      int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel2,4,sG(Imax),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
    Endif

   Else
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel1,8,smin,smax,int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarSkernel1,8,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
     Call DgaussInt(F3BDgaugeSscalarSkernel1,8,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)
    Else
     Call DgaussInt(F3BDgaugeSscalarSkernel1,8,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDgaugeSscalarSkernel1,8,sG(i1),sG(i1+1),int1a,eps)
      int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)
     Enddo
     Call DgaussInt(F3BDgaugeSscalarSkernel1,8,sG(Imax),smax,int1a,eps)
     int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)         
    Endif
   Endif
  Endif

  sumI(1) = (coup(2) * coup(3) + coup(1) * coup(4)) &
        & * (coup(6) * coup(7) + coup(5) * coup(8)) * mf(2) * mf(4) * int1(1)
  sumI(2) = (coup(1) * coup(3) + coup(2) * coup(4)) &
        & * (coup(5) * coup(7) + coup(6) * coup(8)) * mf(1) * mf(3) * int1(2)
  sumI(3) = (coup(2) * coup(3) + coup(1) * coup(4)) &
        & * (coup(5) * coup(7) + coup(6) * coup(8)) * mf(2) * mf(3) * int1(3)
  sumI(4) = (coup(1) * coup(3) + coup(2) * coup(4)) &
        & * (coup(6) * coup(7) + coup(5) * coup(8)) * mf(1) * mf(4) * int1(4)

  erg = -2._dp * Sum( sumI )

 End Subroutine F3BDgaugeSscalarSint


 Subroutine F3BDgaugeSscalarSkernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(2)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI - (mf2(1) - mf2(2))*(mf2(3) - mf2(4)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   sum1 = diff2 + (mf2(1) + mf2(3)) * diff
   erg(1) = erg(1) + sum1 * ReProp
   erg(2) = erg(2) + sum1 * ImProp

   sum1 = - diff2 - (mf2(2) + mf2(4)) * diff
   erg(3) = erg(3) + sum1 * ReProp
   erg(4) = erg(4) + sum1 * ImProp

   sum1 = - diff2 + (sbar - mf2(2) - mf2(3)) * diff
   erg(5) = erg(5) + sum1 * ReProp
   erg(6) = erg(6) + sum1 * ImProp

   sum1 = diff2 - (sbar - mf2(1) - mf2(4)) * diff
   erg(7) = erg(7) + sum1 * ReProp
   erg(8) = erg(8) + sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel1


 Subroutine F3BDgaugeSscalarSkernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) - sbar
   m12 = mf2(1)
   m22 = mf2(2)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   sumI = sumI - (mf2(1) - mf2(2)) * mf2(3)  / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   sum1 = - diff2 - mf2(2) * diff
   erg(1) = erg(1) + sum1 * ReProp
   erg(2) = erg(2) + sum1 * ImProp

   sum1 = - diff2 + (sbar - mf2(2) - mf2(3)) * diff
   erg(3) = erg(3) + sum1 * ReProp
   erg(4) = erg(4) + sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel2


 Subroutine F3BDgaugeSscalarSkernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(2)
   diff = kappa(sbar,m12,m22)

   m22 = mf2(4)
   sumI = sumI + (mf2(1) - mf2(2)) * mf2(4) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   sum1 = diff2 + mf2(1) * diff
   erg(1) = erg(1) + sum1 * ReProp
   erg(2) = erg(2) + sum1 * ImProp

   sum1 = diff2 - (sbar - mf2(1) - mf2(4)) * diff
   erg(3) = erg(3) + sum1 * ReProp
   erg(4) = erg(4) + sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel3


 Subroutine F3BDgaugeSscalarSkernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI - (mf2(1) - mf2(2))*(mf2(3) - mf2(4)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   sum1 = - diff2 - mf2(4) * diff
   erg(1) = erg(1) + sum1 * ReProp
   erg(2) = erg(2) + sum1 * ImProp

   sum1 = diff2 - (sbar - mf2(1) - mf2(4)) * diff
   erg(3) = erg(3) + sum1 * ReProp
   erg(4) = erg(4) + sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel4


 Subroutine F3BDgaugeSscalarSkernel5(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   sumI = sumI - mf2(1) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )          &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )          &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   erg(1) = erg(1) - diff2 * ReProp
   erg(2) = erg(2) - diff2 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel5


 Subroutine F3BDgaugeSscalarSkernel6(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,tmin,tmax,sumI,diff,sbar,ReProp,ImProp,diff2

  erg = 0._dp

  Do i1 = 1,2
   sbar = s(i1)
   sumI = mf2(1) +  mf2(4) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m22 = mf2(4)
   sumI = sumI + mf2(1) * mf2(4) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)
   diff2 = 0.5_dp * (tmin**2 - tmax**2)

   ReProp = ( (sbar-mS2) * (sbar-mT2) + mSgS * mTgT )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )
   ImProp = ( (sbar-mS2) * mTgT - mSgS * (sbar-mT2) )           &
        & / ( ((sbar-mS2)**2 + mSgS2) * ((sbar-mT2)**2 + mTgT2) )

   sum1 = diff2 - (sbar - mf2(1) - mf2(4)) * diff
   erg(1) = erg(1) + sum1 * ReProp
   erg(2) = erg(2) + sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarSkernel6


 Subroutine F3BDgaugeSscalarTint(gauge,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the integral of the matrix element squared of
 ! M^2_fi = - Tr[\[Gamma]_\[Mu].(C[1] Pl + C[2] Pr).(P1+mf[1]).
 !           (C[3] Pl + C[4] Pr).(P2+mf[2]).\[Gamma]^\[Mu].
 !           (C[5] Pl + C[6] Pr).(P3+mf[3]).
 !           (C[7] Pl + C[8] Pr).(P4+mf[4])] 
 !           ( (p1-p4)^2 - m^2_G  + I m^2_G \[CapitalGamma]_G )
 !            * (p1-p2)^2 - m^2_S  - I m^2_S \[CapitalGamma]_S ) )
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p4)^2, t=(p1-p2)^2, and u=(p1-p3)^2.
 ! Here {C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8]} 
 ! are the couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses,
 ! gauge = {m_G^2 , m_G \[CapitalGamma]_G, m_S^2, m_S \[CapitalGamma]_S} 
 ! are the masses and the total decay widths of the gauge boson and the scalar
 ! boson, respectively.
 ! written by Werner Porod, 10.1.2000
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(4),mf(4),eps
  Complex(dp), Intent(in) :: coup(8)
  Complex(dp), Intent(inout) :: int1(8)
  Complex(dp), Intent(out) :: erg
  Logical, Intent(inout) :: Integrate

  Integer :: i1,Imin,Imax
  Real(dp) :: smin,smax,sG(2),int1a(16),int1b(8),int1c(4),int1d(2)
  Complex(dp) :: sumI(8)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   Do i1=1,8
    int1(i1) = (0._dp,0._dp)
   Enddo
   erg = (0._dp,0._dp)
   Return
  Endif

  mf2 = mf**2  

  sumI(8) = -8._dp * ( coup(2) * coup(4) * coup(5) * coup(7)  &
          &         + coup(1) * coup(3) * coup(6) * coup(8)) &
          &       * mf(1) * mf(2) * mf(3) * mf(4)
  sumI(7) = 2._dp * ( coup(1) * coup(4) * coup(5) * coup(7)   &
          &        + coup(2) * coup(3) * coup(6) * coup(8)) * mf(2) * mf(3)
  sumI(6) = 2._dp * ( coup(2) * coup(4) * coup(6) * coup(7)   &
          &        + coup(1) * coup(3) * coup(5) * coup(8)) * mf(1) * mf(4)
  sumI(5) = -4._dp * ( coup(2) * coup(3) * coup(5) * coup(7)   &
          &         + coup(1) * coup(4) * coup(6) * coup(8)) *  mf(3) * mf(4)
  sumI(4) = -4._dp * ( coup(1) * coup(3) * coup(6) * coup(7)   &
          &          + coup(2) * coup(4) * coup(5) * coup(8)) * mf(1) * mf(2)
  sumI(3) = -2._dp * ( coup(1) * coup(4) * coup(6) * coup(7)   &
          &         + coup(2) * coup(3) * coup(5) * coup(8))
  sumI(2) = 2._dp * ( coup(2) * coup(3) * coup(6) * coup(7)    &
          &        + coup(1) * coup(4) * coup(5) * coup(8)) * mf(2) * mf(4)
  sumI(1) = 2._dp * ( coup(1) * coup(3) * coup(5) * coup(7)    &
                  + coup(2) * coup(4) * coup(6) * coup(8)) * mf(1) * mf(3)
  
  erg = Sum( Abs(sumI) )
  If (erg.Eq.ZeroC) Then
   If (Integrate) Integrate = .False.
   Return
  Endif


  If (Integrate) Then
   mS2 = gauge(1)**2
   mSgS = gauge(1)*gauge(2)
   mSgS2 = mSgS**2
   mT2 = gauge(3)**2
   mTgT = gauge(3)*gauge(4)
   mTgT2 = mTgT**2

   int1 = ZeroC

   smax = (Abs(mf(1))-Abs(mf(4)))**2
   smin = (Abs(mf(3))+Abs(mf(2)))**2
   sG(1) = (gauge(1)-2._dp*gauge(2))**2
   sG(2) = (gauge(1)+2._dp*gauge(2))**2
   Imin = 3
   Imax = 0
   Do i1=1,2
    If (smin.Lt.sG(3-i1)) Imin = 3 - i1
    If (smax.Gt.sG(i1)) Imax = i1
   Enddo

   If ((mf(2).Ne.0._dp).And.(mf(3).Ne.0._dp).And.(mf(4).Ne.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,smin,smax,int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,sG(Imin),sG(Imax),int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSscalarTkernel1,16,sG(Imax),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,smin,smax,int1d,eps)
     int1(3) = int1d(1) + Ic * int1d(2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(3) = int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,sG(Imin),smax,int1d,eps)
     int1(3) = int1(3) + int1d(1) + Ic * int1d(2)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(3) = int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,sG(Imin),sG(Imax),int1d,eps)
     int1(3) = int1(3) + int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSscalarTkernel8,2,sG(Imax),smax,int1d,eps)
     int1(3) = int1(3) + int1d(1) + Ic * int1d(2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,smin,smax,int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,sG(Imin),smax,int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,sG(Imin),sG(Imax),int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel7,4,sG(Imax),smax,int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,smin,smax,int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(3) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(3) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,sG(Imin),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(3) = int1(3) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(3) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(3) = int1(3) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel6,4,sG(Imax),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(3) = int1(3) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,smin,smax,int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,sG(Imin),smax,int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(3) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,sG(Imin),sG(Imax),int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSscalarTkernel5,4,sG(Imax),smax,int1c,eps)
     int1(3) = int1(3) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif (mf(2).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel4,8,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(3).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,smin,smax,int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,sG(Imin),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,sG(Imin),sG(Imax),int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel3,8,sG(Imax),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(4).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSscalarTkernel2,8,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)
    Endif

   Endif
  Endif

  erg = Sum(sumI * int1)
  
 End Subroutine F3BDgaugeSscalarTint


 Subroutine F3BDgaugeSscalarTkernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(16)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + (mf2(1) - mf2(4))*(mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan +  mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

   sum1 = ( mf2(1) + mf2(3) + mf2(2) + mf2(4)                                 &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)                        &
      & - ( mf2(1) + mf2(3) + mf2(2) + mf2(4)                                 &
      &   - 1.5_dp * mT2- 0.5_dp * tmin) * (tmin - mT2)                         &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(2) - mf2(4) ) *DiffTan &
      & + ( mTgT2 + ( mf2(1) + mf2(2) - mT2)                                  &
                   * ( mT2 - mf2(3) - mf2(4) ) ) * DiffLog
   sum2 = - mTgT * diff                                                       &
      & + ( mTgT2 - ( mT2 - mf2(1) - mf2(2) )                                 &
      &             * ( mT2 - mf2(3) - mf2(4) )  ) * DiffTan                  &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(2) - mf2(4) ) * mTgT * DiffLog

   erg(5) = erg(5) + sum1 * ReProp + sum2 * ImProp
   erg(6) = erg(6) + sum2 * ReProp - sum1 * ImProp

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(3) - mf2(4) ) * DiffLog
   sum2 = ( mT2 - mf2(3) - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(7) = erg(7) + sum1 * ReProp + sum2 * ImProp
   erg(8) = erg(8) + sum2 * ReProp - sum1 * ImProp

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) - mf2(2) ) * DiffLog
   sum2 = - ( mT2 - mf2(1) - mf2(2) ) * DiffTan - mTgT * DiffLog

   erg(9) = erg(9) + sum1 * ReProp + sum2 * ImProp
   erg(10) = erg(10) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * (sbar - mf2(2) - mf2(3))
   sum2 = DiffTan * (sbar - mf2(2) - mf2(3))

   erg(11) = erg(11) + sum1 * ReProp + sum2 * ImProp
   erg(12) = erg(12) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * (-sbar + mf2(1) + mf2(4))
   sum2 = DiffTan * (-sbar + mf2(1) + mf2(4))

   erg(13) = erg(13) + sum1 * ReProp + sum2 * ImProp
   erg(14) = erg(14) + sum2 * ReProp - sum1 * ImProp

   erg(15) = erg(15) + DiffLog * ReProp + DiffTan * ImProp
   erg(16) = erg(16) + DiffTan * ReProp - DiffLog * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel1


 Subroutine F3BDgaugeSscalarTkernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + mf2(1) * (mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = (mf2(1) + mf2(3) + mf2(2) - 1.5_dp* mT2- 0.5_dp* tmax) * (tmax - mT2) &
      & - (mf2(1) + mf2(3) + mf2(2) - 1.5_dp* mT2- 0.5_dp* tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(2) ) * DiffTan         &
      & + ( mTgT2 + ( mf2(1) + mf2(2) - mT2) * ( mT2 - mf2(3) ) ) * DiffLog
   sum2 = - mTgT * diff                                                       &
      & + ( mTgT2 - ( mT2 - mf2(1) - mf2(2) ) * ( mT2 - mf2(3) )  ) * DiffTan &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(2) ) * mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(3) ) * DiffLog
   sum2 = ( mT2 - mf2(3) ) * DiffTan + mTgT * DiffLog

   erg(5) = erg(5) + sum1 * ReProp + sum2 * ImProp
   erg(6) = erg(6) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * (-sbar + mf2(1))
   sum2 = DiffTan * (-sbar + mf2(1))

   erg(7) = erg(7) + sum1 * ReProp + sum2 * ImProp
   erg(8) = erg(8) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel2


 Subroutine F3BDgaugeSscalarTkernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m22 = mf2(2)
   sumI = sumI - (mf2(1) - mf2(4)) * mf2(2) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan +  mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = (mf2(1) + mf2(2) + mf2(4) - 1.5_dp* mT2- 0.5_dp* tmax) * (tmax - mT2) &
      & - (mf2(1) + mf2(2) + mf2(4) - 1.5_dp* mT2- 0.5_dp* tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(2) - mf2(4) ) * DiffTan         &
      & + ( mTgT2 + ( mf2(1) + mf2(2) - mT2) * ( mT2 - mf2(4) ) ) * DiffLog
   sum2 = - mTgT * diff                                                       &
      & + ( mTgT2 - ( mT2 - mf2(1) - mf2(2) ) * ( mT2 - mf2(4) ) ) * DiffTan  &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(2) - mf2(4) ) * mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(4) ) * DiffLog
   sum2 = ( mT2 - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(5) = erg(5) + sum1 * ReProp + sum2 * ImProp
   erg(6) = erg(6) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * (sbar - mf2(2))
   sum2 = DiffTan * (sbar - mf2(2))

   erg(7) = erg(7) + sum1 * ReProp + sum2 * ImProp
   erg(8) = erg(8) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel3


 Subroutine F3BDgaugeSscalarTkernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   sumI = sumI + (mf2(1) - mf2(4)) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = (mf2(1) + mf2(3) + mf2(4) - 1.5_dp* mT2- 0.5_dp* tmax) * (tmax - mT2) &
      & - (mf2(1) + mf2(3) + mf2(4) - 1.5_dp* mT2- 0.5_dp* tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(4) ) * DiffTan         &
      & + ( mTgT2 + ( mf2(1) - mT2) * ( mT2 - mf2(3) - mf2(4) ) ) * DiffLog
   sum2 = - mTgT * diff                                                       &
      & + ( mTgT2 - ( mT2 - mf2(1) ) * ( mT2 - mf2(3) - mf2(4) )  ) * DiffTan &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(3) - mf2(4) ) * mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) ) * DiffLog
   sum2 = - ( mT2 - mf2(1) ) * DiffTan - mTgT * DiffLog

   erg(5) = erg(5) + sum1 * ReProp + sum2 * ImProp
   erg(6) = erg(6) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * (sbar - mf2(3))
   sum2 = DiffTan * (sbar - mf2(3))

   erg(7) = erg(7) + sum1 * ReProp + sum2 * ImProp
   erg(8) = erg(8) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel4


 Subroutine F3BDgaugeSscalarTkernel5(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m22 = mf2(2)
   sumI = sumI - mf2(1) * mf2(2) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = ( mf2(1) + mf2(2) - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2) &
      & - ( mf2(1) + mf2(2) - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(2) ) * DiffTan             &
      & + ( mTgT2 + ( mf2(1) + mf2(2) - mT2) * mT2 ) * DiffLog
   sum2 = - mTgT * diff                                                  &
      & + ( mTgT2 - ( mT2 - mf2(1) - mf2(2) ) * mT2 ) * DiffTan          &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(2) ) * mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = diff - mTgT * DiffTan + mT2  * DiffLog
   sum2 = mT2 * DiffTan + mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel5


 Subroutine F3BDgaugeSscalarTkernel6(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12, sum1, sum2, tmin, tmax, sumI, diff, DiffTan &
         &        ,DiffLog, sbar, ReProp, ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   sumI = sumI + mf2(1) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = ( mf2(1) + mf2(3) - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2) &
      & - ( mf2(1) + mf2(3) - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(3) ) * DiffTan             &
      & + ( mTgT2 + ( mf2(1) - mT2) * ( mT2 - mf2(3) ) ) * DiffLog
   sum2 = - mTgT * diff                                                  &
      & + ( mTgT2 - ( mT2 - mf2(1) ) * ( mT2 - mf2(3) )  ) * DiffTan     &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(3) ) * mTgT * DiffLog

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel6


 Subroutine F3BDgaugeSscalarTkernel7(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,DiffTan &
         &        ,DiffLog,sbar,ReProp,ImProp

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
           &                / ((tmin - mT2)**2 + mTgT2 ) ) )

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = ( mf2(1) + mf2(4) - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2) &
      & - ( mf2(1) + mf2(4) - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) - mf2(4) ) * DiffTan             &
      & + ( mTgT2 + ( mf2(1) - mT2) * ( mT2 - mf2(4) ) ) * DiffLog
   sum2 = - mTgT * diff                                                  &
      & + ( mTgT2 - ( mT2 - mf2(1) ) * ( mT2 - mf2(4) )  ) * DiffTan     &
      & - ( 2._dp * mT2 -  mf2(1) - mf2(4) ) * mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

   sum1 = DiffLog * sbar
   sum2 = DiffTan * sbar

   erg(3) = erg(3) + sum1 * ReProp + sum2 * ImProp
   erg(4) = erg(4) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel7


 Subroutine F3BDgaugeSscalarTkernel8(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSscalarT for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12, sum1, sum2, tmax, sumI, diff, DiffTan &
         &        ,DiffLog, sbar, ReProp, ImProp, wert1

  erg = 0._dp
  if (mTgT.gt.0._dp) wert1 =  Atan( -mT2 / mTgT )
  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   tmax = sumI

   if (mTgT.gt.0._dp) then
    DiffTan = Atan( (tmax-mT2) / mTgT ) - wert1 
   else
    DiffTan = 0._dp
   end if
   DiffLog = 0.5_dp  * Log( ( (tmax - mT2)**2 + mTgT2 ) /( mT2**2 + mTgT2 ))

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = ( mf2(1) - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2) &
      & - ( mf2(1) - 1.5_dp * mT2) * ( - mT2) &
      & + mTgT * ( 2._dp * mT2 -  mf2(1) ) * DiffTan             &
      & + ( mTgT2 + ( mf2(1) - mT2) * mT2 ) * DiffLog
   sum2 = - mTgT * diff                                         &
      & + ( mTgT2 - ( mT2 - mf2(1) ) * mT2 ) * DiffTan          &
      & - ( 2._dp * mT2 -  mf2(1) ) * mTgT * DiffLog

   erg(1) = erg(1) + sum1 * ReProp + sum2 * ImProp
   erg(2) = erg(2) + sum2 * ReProp - sum1 * ImProp

  Enddo

 End Subroutine F3BDgaugeSscalarTkernel8


 Subroutine F3BDgaugeSSint(gauge,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the matrix element squared of
 ! M^2_fi = Tr[\[Gamma]_\[Mu].(C[1] Pl + C[2] Pr).(P1+mf[1]).\[Gamma]_\[Sigma].
 !        (Cc[1] Pl + Cc[2] Pr).(P2+mf[2])] *
 !      Tr[\[Gamma]^\[Mu].(C[3] Pl + C[4] Pr).(P3+mf[3]).\[Gamma]^\[Sigma].
 !        (Cc[3] Pl + Cc[4] Pr).(P4+mf[4])] /
 !      ( (p1-p2)^2 - m^2_G )^2 + m^2_G \[CapitalGamma]^2_G )\n
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p2)^2, t=(p1-p3)^2, and u=(p1-p4)^2.
 ! Here {C[1],C[2],C[3],C[4]} are the couplings, Cc are the complex conjugated
 ! couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses, and
 ! {m_G,\[CapitalGamma]_G} are the mass and the total decay width of the gauge
 ! boson. It is integrated over the kinematical allowed phase space.
 ! written by Werner Porod, 15.11.1999
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(2),mf(4),eps
  Real(dp), Intent(inout) :: int1(4)
  Real(dp), Intent(out) :: erg
  Complex(dp), Intent(in) :: coup(4)
  Logical, Intent(inout) :: integrate

  Integer :: i1,Imax,Imin
  Real(dp) :: coup2(4),smin,smax,sG(2),int1a(4),int1b(2),testR,mR    &
         &       , sminG,smaxG
  Complex(dp) :: coupC(4)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   int1 = 0._dp
   erg = 0._dp
   Return
  Endif

  Do i1=1,4
   coupC(i1) = Conjg( coup(i1) )
   coup2(i1) = Abs( coup(i1) )**2
   mf2(i1) = mf(i1)**2
  Enddo

  If (((coup2(1) + coup2(2)).Eq.0._dp).Or.((coup2(3) + coup2(4)).Eq.0._dp)) Then
   If (Integrate) Integrate = .False.
   erg = 0._dp
   Return
  Endif


  If (Integrate) Then 
   mG2 = gauge(1)**2
   mGgG2 = mG2 * gauge(2)**2
   smax = (Abs(mf(1)) - Abs(mf(2)))**2
   smin = (Abs(mf(3)) + Abs(mf(4)))**2
   sG(1) = (gauge(1)-5._dp*gauge(2))**2
   sG(2) = (gauge(1)+5._dp*gauge(2))**2
   Imin = 3
   Imax = 0
   Do i1=1,2
    If (smin.Lt.sG(3-i1)) Imin = 3 - i1
    If (smax.Gt.sG(i1)) Imax = i1
   Enddo

   int1 = 0._dp

   If ((mf(2).Ne.0._dp).And.(mf(3).Ne.0._dp).And.(mf(4).Ne.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSSkernel1,4,smin,smax,int1,eps)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSSkernel1,4,smin,sG(Imin),int1,eps)
     Call DgaussInt(F3BDgaugeSSkernel1,4,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a
    Else
     Call DgaussInt(F3BDgaugeSSkernel1,4,smin,sG(1),int1,eps)
     Call DgaussInt(F3BDgaugeSSkernel1,4,sG(1),sG(2),int1a,eps)
     int1 = int1 + int1a
     Call DgaussInt(F3BDgaugeSSkernel1,4,sG(2),smax,int1a,eps)
     int1 = int1 + int1a
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    mR = mf(1)
    Call IntGaugeSS1(gauge,mR,testR)
    int1(1)=testR

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     int1(1) = dgauss(F3BDgaugeSSkernel2,smin,smax,eps)
    Elseif (Imin.Eq.Imax) Then
     sminG = sG(Imin)
     int1(1) = dgauss(F3BDgaugeSSkernel2,smin,sminG,eps) &
           & + dgauss(F3BDgaugeSSkernel2,sminG,smax,eps)
    Else
     sminG = sG(1)
     sminG = sG(1)
     smaxG = sG(2)
     int1(1) = dgauss(F3BDgaugeSSkernel2,smin,sminG,eps)  &
           & + dgauss(F3BDgaugeSSkernel2,sminG,smaxG,eps) &
           & + dgauss(F3BDgaugeSSkernel2,smaxG,smax,eps)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     int1(1) = dgauss(F3BDgaugeSSkernel3,smin,smax,eps)
    Elseif (Imin.Eq.Imax) Then
     sminG = sG(Imin)
     int1(1) = dgauss(F3BDgaugeSSkernel3,smin,sminG,eps) &
           & + dgauss(F3BDgaugeSSkernel3,sminG,smax,eps)
    Else
     sminG = sG(1)
     sminG = sG(1)
     smaxG = sG(2)
     int1(1) = dgauss(F3BDgaugeSSkernel3,smin,sminG,eps)  &
           & + dgauss(F3BDgaugeSSkernel3,sminG,smaxG,eps) &
           & + dgauss(F3BDgaugeSSkernel3,smaxG,smax,eps)
    Endif

   Elseif ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSSkernel4,2,smin,smax,int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSSkernel4,2,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel4,2,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Else
     Call DgaussInt(F3BDgaugeSSkernel4,2,smin,sG(1),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel4,2,sG(1),sG(2),int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel4,2,sG(2),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Endif

   Elseif (mf(4).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSSkernel5,2,smin,smax,int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSSkernel5,2,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel5,2,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Else
     Call DgaussInt(F3BDgaugeSSkernel5,2,smin,sG(1),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel5,2,sG(1),sG(2),int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel5,2,sG(2),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Endif

   Elseif (mf(3).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSSkernel6,2,smin,smax,int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSSkernel6,2,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel6,2,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Else
     Call DgaussInt(F3BDgaugeSSkernel6,2,smin,sG(1),int1b,eps)
     int1(1) = int1b(1)
     int1(2) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel6,2,sG(1),sG(2),int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel6,2,sG(2),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(2) = int1(2) + int1b(2)
    Endif

   Elseif (mf(2).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSSkernel7,2,smin,smax,int1b,eps)
     int1(1) = int1b(1)
     int1(3) = int1b(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSSkernel7,2,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1)
     int1(3) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel7,2,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(3) = int1(3) + int1b(2)
    Else
     Call DgaussInt(F3BDgaugeSSkernel7,2,smin,sG(1),int1b,eps)
     int1(1) = int1b(1)
     int1(3) = int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel7,2,sG(1),sG(2),int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(3) = int1(3) + int1b(2)
     Call DgaussInt(F3BDgaugeSSkernel7,2,sG(2),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1)
     int1(3) = int1(3) + int1b(2)
    Endif
   Endif
 
  Endif

  If ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
   erg = 2._dp * ( coup2(1) + coup2(2) ) * ( coup2(3) + coup2(4) ) * int1(1)

  Elseif ((mf(2).Eq.0._dp).And.((mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp))) Then
   erg = 2._dp * ( coup2(1) + coup2(2) ) * ( coup2(3) + coup2(4) ) * int1(1)

  Elseif ((mf(2).Ne.0._dp).And.((mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp))) Then

   erg = 2._dp * ( coup2(1) + coup2(2) ) * ( coup2(3) + coup2(4) ) * int1(1) &
     & - 8._dp * Real( coup(1) * coupC(2),dp )                               &
     &        * ( coup2(3) + coup2(4) ) * mf(1) * mf(2) * int1(2)

  Elseif ((mf(2).Eq.0._dp).And.(mf(3).Ne.0._dp).And.(mf(4).Ne.0._dp)) Then
  
   erg = 2._dp * ( coup2(1) + coup2(2) ) * ( coup2(3) + coup2(4) ) * int1(1) &
     & - 8._dp * ( coup2(1) + coup2(2) )                                     &
     &        * Real( coup(3) * coupC(4),dp ) * mf(3) * mf(4) * int1(3)

  Else

   erg = 64._dp * Real(coup(2) * coupC(1),dp ) * Real(coup(3) * coupC(4),dp ) &
     &         * mf(1) * mf(2) * mf(3) * mf(4) * int1(4)                      &
     & - 8._dp * ( coup2(1) + coup2(2) )                                      &
     &        * Real( coup(3) * coupC(4),dp ) * mf(3) * mf(4) * int1(3)       &
     & - 8._dp * Real( coup(1) * coupC(2),dp )                                &
     &        * ( coup2(3) + coup2(4) ) * mf(1) * mf(2) * int1(2)             &
     & + 2._dp * ( coup2(1) + coup2(2) ) * ( coup2(3) + coup2(4) ) * int1(1)

  Endif

 End Subroutine F3BDgaugeSSint


 Subroutine F3BDgaugeSSKernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,ooProp,sbar

  erg = 0._dp
 
  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sum1 = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(4)
   sum1 = sum1 *  kappa(sbar,m12,m22) / sbar

   sum2 = sum1 * ( -2._dp * sbar**4                                   &
        &        + (mf2(1)+mf2(2)+mf2(3)+mf2(4)) * sbar**3           &
        &        + ( (mf2(1)-mf2(2))**2 + (mf2(3)-mf2(4))**2         &
        &          - 2._dp * (mf2(1)+mf2(2)) * (mf2(3)+mf2(4))        &
        &          ) * sbar**2                                       &
        &        + ( (mf2(1)+mf2(2)) * (mf2(3)-mf2(4))**2            &
        &          + (mf2(3)+mf2(4)) * (mf2(1)-mf2(2))**2 ) * sbar   &
        &        - 2._dp * (mf2(1)-mf2(2))**2 * (mf2(3)-mf2(4))**2    &
        &        ) / (3._dp * sbar**2)

   ooProp = 1._dp / ( (sbar-mG2)**2 + mGgG2 )
   erg(1) = erg(1) + sum2 * ooProp
   erg(2) = erg(2) + sum1 * ( sbar - mf2(3) - mf2(4) ) * ooProp
   erg(3) = erg(3) + sum1 * (-sbar + mf2(1) + mf2(2) ) * ooProp
   erg(4) = erg(4) + sum1 * ooProp
  Enddo

 End Subroutine F3BDgaugeSSKernel1


 Real(dp) Function F3BDgaugeSSkernel2(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  s

  Real(dp) :: sum1, sum2

  sum1 = Abs(s-mf2(1))*Abs(s-mf2(4))

  sum2 = sum1 * ( -2._dp * s**4                                       &
     &           + (mf2(1)+mf2(4)) * s**3                            &
     &           + ( mf2(1) - mf2(4) )**2 * s**2                     &
     &           + ( mf2(1) * mf2(4)**2 + mf2(4) * mf2(1)**2 ) * s   &
     &           - 2._dp * mf2(1)**2 * mf2(4)**2                      &
     &           ) / (3._dp * s**3)

  F3BDgaugeSSkernel2 = sum2 / ( (s-mG2)**2 + mGgG2 )

 End Function F3BDgaugeSSkernel2


 Real(dp) Function F3BDgaugeSSkernel3(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  s

  Real(dp) :: sum1,sum2

  sum1 = Abs(s-mf2(1))*Abs(s-mf2(3))

  sum2 = sum1 * ( -2._dp * s**4                                      &
     &           + (mf2(1)+mf2(3)) * s**3                           &
     &           + (mf2(1) - mf2(3))**2 * s**2                      &
     &           + ( mf2(1) * mf2(3)**2 + mf2(3) * mf2(1)**2 ) * s  &
     &           - 2._dp * mf2(1)**2 * mf2(3)**2                     &
     &           ) / (3._dp * s**3)

  F3BDgaugeSSkernel3 = sum2 / ( (s-mG2)**2 + mGgG2 )

 End Function F3BDgaugeSSkernel3


 Subroutine F3BDgaugeSSKernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,ooProp,sbar

  erg = 0._dp
 
  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sum1 = kappa(sbar,m12,m22)

   sum2 = sum1 * ( -2._dp * sbar**2 + (mf2(1)+mf2(2)) * sbar &
        &        +  (mf2(1)-mf2(2))**2 ) / 3._dp 

   ooProp = 1._dp / ( (sbar-mG2)**2 + mGgG2 )
   erg(1) = erg(1) + sum2 * ooProp
   erg(2) = erg(2) + sum1 *  sbar * ooProp
  Enddo

 End Subroutine F3BDgaugeSSKernel4


 Subroutine F3BDgaugeSSKernel5(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,ooProp,sbar

  erg = 0._dp
 
  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sum1 = kappa(sbar,m12,m22)

   m12 = mf2(3)
   sum1 = sum1 * Abs(sbar-m12) / sbar
   sum2 = sum1 * ( -2._dp * sbar**4                                   &
        &        + (mf2(1)+mf2(2)+mf2(3)) * sbar**3                  &
        &        + ( (mf2(1)-mf2(2))**2 + mf2(3)**2                  &
        &          - 2._dp * (mf2(1)+mf2(2)) * mf2(3) ) * sbar**2     &
        &        + ( (mf2(1)+mf2(2)) * mf2(3)**2                     &
        &          + mf2(3) * (mf2(1)-mf2(2))**2 ) * sbar            &
        &        - 2._dp * (mf2(1)-mf2(2))**2 * mf2(3)**2             &
        &        ) / (3._dp * sbar**2)

   ooProp = 1._dp / ( (sbar-mG2)**2 + mGgG2 )
   erg(1) = erg(1) + sum2 * ooProp
   erg(2) = erg(2) + sum1 * ( sbar - mf2(3) ) * ooProp
  Enddo

 End Subroutine F3BDgaugeSSKernel5


 Subroutine F3BDgaugeSSKernel6(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,ooProp,sbar

  erg = 0._dp
 
  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sum1 = kappa(sbar,m12,m22)

   m12 = mf2(4)
   sum1 = sum1 *  Abs(sbar-m12) / sbar
   sum2 = sum1 * ( -2._dp * sbar**4                            &
      &          + (mf2(1)+mf2(2)+mf2(4)) * sbar**3           &
      &          + ( (mf2(1)-mf2(2))**2 + mf2(4)**2           &
      &            - 2._dp * (mf2(1)+mf2(2)) * mf2(4)          &
      &            ) * sbar**2                                &
      &          + ( (mf2(1)+mf2(2)) * mf2(4)**2              &
      &            + mf2(4) * (mf2(1)-mf2(2))**2 ) * sbar     &
      &          - 2._dp * (mf2(1)-mf2(2))**2 * mf2(4)**2      &
      &          ) / (3._dp * sbar**2)

   ooProp = 1._dp / ( (sbar-mG2)**2 + mGgG2 )
   erg(1) = erg(1) + sum2 * ooProp
   erg(2) = erg(2) + sum1 * ( sbar - mf2(4) ) * ooProp
  Enddo

 End Subroutine F3BDgaugeSSKernel6


 Subroutine F3BDgaugeSSKernel7(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 3.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,ooProp,sbar

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   sum1 = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(4)
   sum1 = sum1 *  kappa(sbar,m12,m22) / sbar
   sum2 = sum1 * ( -2._dp * sbar**4                          &
      &          + (mf2(1)+mf2(3)+mf2(4)) * sbar**3         &
      &          + ( mf2(1)**2 + (mf2(3)-mf2(4))**2         &
      &            - 2._dp * mf2(1) * (mf2(3)+mf2(4))        &
      &            ) * sbar**2                              &
      &          + ( mf2(1) * (mf2(3)-mf2(4))**2            &
      &            + (mf2(3)+mf2(4)) * mf2(1)**2 ) * sbar   &
      &          - 2._dp * mf2(1)**2 * (mf2(3)-mf2(4))**2    &
      &          ) / (3._dp * sbar**2)

   ooProp = 1._dp / ( (sbar-mG2)**2 + mGgG2 )
   erg(1) = erg(1) + sum2 * ooProp
   erg(2) = erg(2) + sum1 * (-sbar + mf2(1) ) * ooProp
  Enddo

 End Subroutine F3BDgaugeSSKernel7


 Subroutine F3BDgaugeSTint(gauge,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the integral of the matrix element squared of
 ! M^2_fi = Tr[\[Gamma]_\[Mu].(C[1] Pl + C[2] Pr).(P1+mf[1]).\[Gamma]_\[Sigma].
 !         (C[3] Pl + C[4] Pr).(P2+mf[2]).\[Gamma]^\[Mu].
 !         (C[5] Pl + C[6] Pr).(P3+mf[3]).
 !         \[Gamma]^\[Sigma].(C[7] Pl + C[8] Pr).(P4+mf[4])] 
 !    / ( (p1-p4)^2 - m^2_S  + I m^2_S \[CapitalGamma]_S )
 !      (p1-p3)^2 - m^2_T  - I m^2_T \[CapitalGamma]_T ) )
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p4)^2, t=(p1-p2)^2, and u=(p1-p3)^2.
 ! Here {C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8]} 
 ! are the couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses,
 ! gauge = {m_S^2, m_S \[CapitalGamma]_S, m_T^2, m_T \[CapitalGamma]_T} 
 ! are the masses squared and the mass times the total decay widths of the
 ! gauge boson in the s- and t-channel, respectively.
 ! written by Werner Porod, 4.1.2000
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(4),mf(4),eps
  Complex(dp), Intent(in) :: coup(8)
  Complex(dp), Intent(inout) :: int1(8)
  Complex(dp), Intent(out) :: erg
  Logical, Intent(inout) :: Integrate

  Integer :: i1,Imin,Imax
  Real(dp) :: smin,smax,sG(2),int1a(16),int1b(8),int1c(4),int1d(2)
  Complex(dp) :: sumI(8)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   erg = (0._dp,0._dp)
   Return
  Endif

   sumI(8) = -16._dp * ( coup(2) * coup(4) * coup(5) * coup(7)   &
           &          + coup(1) * coup(3) * coup(6) * coup(8))  &
           &        * mf(1) * mf(2) * mf(3) * mf(4)
   sumI(7) = 4._dp * ( coup(1) * coup(4) * coup(5) * coup(7)   &
           &        + coup(2) * coup(3) * coup(6) * coup(8)) * mf(2) * mf(3)
   sumI(6) = 4._dp * ( coup(2) * coup(3) * coup(5) * coup(7)   &
           &        + coup(1) * coup(4) * coup(6) * coup(8)) *  mf(1) * mf(4)
   sumI(5) = 4._dp * ( coup(2) * coup(4) * coup(6) * coup(7)   &
           &        + coup(1) * coup(3) * coup(5) * coup(8)) * mf(3) * mf(4)
   sumI(4) = 4._dp * ( coup(1) * coup(3) * coup(6) * coup(7)   &
           &        + coup(2) * coup(4) * coup(5) * coup(8)) * mf(1) * mf(2)
   sumI(3) = 4._dp * ( coup(2) * coup(3) * coup(6) * coup(7)   &
           &        + coup(1) * coup(4) * coup(5) * coup(8)) * mf(2) * mf(4)
   sumI(2) = 4._dp * ( coup(1) * coup(4) * coup(6) * coup(7)   &
           &        + coup(2) * coup(3) * coup(5) * coup(8)) * mf(1) * mf(3)
   sumI(1) = -4._dp * ( coup(1) * coup(3) * coup(5) * coup(7)   &
           &         + coup(2) * coup(4) * coup(6) * coup(8) )

  erg = Sum( Abs(sumI) )

  If ( erg.Eq.ZeroC ) Then
   If (Integrate) Integrate = .False.
   erg = (0._dp,0._dp)
   Return
  Endif

  mf2 = mf**2

  If (Integrate) Then
   mS2 = gauge(1)**2
   mSgS = gauge(1)*gauge(2)
   mSgS2 = mSgS**2
   mT2 = gauge(3)**2
   mTgT = gauge(3)*gauge(4)
   mTgT2 = mTgT**2
   int1 = (0._dp,0._dp)
   smax = (Abs(mf(1))-Abs(mf(4)))**2
   smin = (Abs(mf(3))+Abs(mf(2)))**2
   sG(1) = (gauge(1)-2._dp*gauge(2))**2
   sG(2) = (gauge(1)+2._dp*gauge(2))**2
   Imin = 3
   Imax = 0
   Do i1=1,2
    If (smin.Lt.sG(3-i1)) Imin = 3 - i1
    If (smax.Gt.sG(i1)) Imax = i1
   Enddo

   If ((mf(2).Ne.0._dp).And.(mf(3).Ne.0._dp).And.(mf(4).Ne.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel1,16,smin,smax,int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSTkernel1,16,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)

    Else
     Call DgaussInt(F3BDgaugeSTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSTkernel1,16,sG(Imin),sG(Imax),int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDgaugeSTkernel1,16,sG(Imax),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel8,2,smin,smax,int1d,eps)
     int1(1)= int1d(1) + Ic * int1d(2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(1)= int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSTkernel8,2,sG(Imin),smax,int1d,eps)
     int1(1)= int1(1) + int1d(1) + Ic * int1d(2)

    Else
     Call DgaussInt(F3BDgaugeSTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(1)= int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSTkernel8,2,sG(Imin),sG(Imax),int1d,eps)
     int1(1)= int1(1) + int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDgaugeSTkernel8,2,sG(Imax),smax,int1d,eps)
     int1(1)= int1(1) + int1d(1) + Ic * int1d(2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel5,4,smin,smax,int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(6)= int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(6)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel5,4,sG(Imin),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(6)= int1(6) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(6)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel5,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(6)= int1(6) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel5,4,sG(Imax),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(6)= int1(6) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel6,4,smin,smax,int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(2)= int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(2)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel6,4,sG(Imin),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(2)= int1(2) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(2)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel6,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(2)= int1(2) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel6,4,sG(Imax),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(2)= int1(2) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel7,4,smin,smax,int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(4)= int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(4)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel7,4,sG(Imin),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(4)= int1(4) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDgaugeSTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(1)= int1c(1) + Ic * int1c(2)
     int1(4)= int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel7,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(4)= int1(4) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDgaugeSTkernel7,4,sG(Imax),smax,int1c,eps)
     int1(1)= int1(1) + int1c(1) + Ic * int1c(2)
     int1(4)= int1(4) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif (mf(2).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel2,8,smin,smax,int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(5)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(5)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel2,8,sG(Imin),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(5)= int1(5) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(5)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel2,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(5)= int1(5) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel2,8,sG(Imax),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(5)= int1(5) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(3).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel3,8,smin,smax,int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(3)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(3)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel3,8,sG(Imin),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(3)= int1(3) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(3)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(6)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel3,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(3)= int1(3) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel3,8,sG(Imax),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(3)= int1(3) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(6)= int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(4).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDgaugeSTkernel4,8,smin,smax,int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(7)= int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDgaugeSTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(7)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel4,8,sG(Imin),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(7)= int1(7) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDgaugeSTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1)= int1b(1) + Ic * int1b(2)
     int1(2)= int1b(3) + Ic * int1b(4)
     int1(4)= int1b(5) + Ic * int1b(6)
     int1(7)= int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel4,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(7)= int1(7) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDgaugeSTkernel4,8,sG(Imax),smax,int1b,eps)
     int1(1)= int1(1) + int1b(1) + Ic * int1b(2)
     int1(2)= int1(2) + int1b(3) + Ic * int1b(4)
     int1(4)= int1(4) + int1b(5) + Ic * int1b(6)
     int1(7)= int1(7) + int1b(7) + Ic * int1b(8)
    Endif

   Endif
  Endif

  erg = Sum(sumI * int1)

 End Subroutine F3BDgaugeSTint

 Subroutine F3BDgaugeSTkernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(16)

  Integer :: i1
  Real(dp) :: sbar,DiffLog,DiffTan,Prop,ReProp,ImProp
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + (mf2(1) - mf2(4))*(mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(3) + mf2(2) + mf2(4) - 2._dp * sbar        &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)             &
      & - ( mf2(1) + mf2(3) + mf2(2) + mf2(4) - 2._dp * sbar        &
      &   - 1.5_dp * mT2- 0.5_dp * tmin) * (tmin - mT2)              &
      & + mTgT * ( 2._dp * (mT2 + sbar)                             &
      &          -  mf2(1) - mf2(3) - mf2(2) - mf2(4) ) * DiffTan  &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) - mf2(4) )               &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffLog

   sum2 = - mTgT * diff                                             &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) - mf2(4) )                &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffTan  &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(3)             &
      &   - mf2(2) - mf2(4) ) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog 

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(5) = erg(5) + ReProp * sum1 + ImProp * sum2 
   erg(6) = erg(6) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(3) - mf2(4) ) * DiffLog
   sum2 = ( mT2 - mf2(3) - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(7) = erg(7) + ReProp * sum1 + ImProp * sum2 
   erg(8) = erg(8) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) - mf2(2) ) * DiffLog
   sum2 = - ( mT2 - mf2(1) - mf2(2) ) * DiffTan - mTgT  * DiffLog

   erg(9) = erg(9) + ReProp * sum1 + ImProp * sum2 
   erg(10) = erg(10) + ReProp * sum2 - ImProp * sum1 
 
   sum1 = DiffLog * (sbar - mf2(2) - mf2(3))
   sum2 = DiffTan * (sbar - mf2(2) - mf2(3))

   erg(11) = erg(11) + ReProp * sum1 + ImProp * sum2 
   erg(12) = erg(12) + ReProp * sum2 - ImProp * sum1 

   sum1 = DiffLog * (-sbar + mf2(1) + mf2(4))
   sum2 = DiffTan * (-sbar + mf2(1) + mf2(4))

   erg(13) = erg(13) + ReProp * sum1 + ImProp * sum2 
   erg(14) = erg(14) + ReProp * sum2 - ImProp * sum1 

   erg(15) = erg(15) + ReProp * DiffLog + ImProp * DiffTan
   erg(16) = erg(16) + ReProp * DiffTan - ImProp * DiffLog

  Enddo

 End Subroutine F3BDgaugeSTkernel1


 Subroutine F3BDgaugeSTkernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: sbar,DiffLog,DiffTan,Prop,ReProp,ImProp
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

    m12 = mf2(3)
    sumI = sumI + (mf2(1) - mf2(4)) * mf2(3) / sbar
    diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(3) + mf2(4) - 2._dp * sbar                   &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)               &
      & - ( mf2(1) + mf2(3) + mf2(4) - 2._dp * sbar                   &
      &   - 1.5_dp * mT2- 0.5_dp * tmin) * (tmin - mT2)                &
      & + mTgT * ( 2._dp * (mT2 + sbar)                               &
      &          -  mf2(1) - mf2(3) - mf2(4) ) * DiffTan             &
      & + ( mTgT2 - ( mT2 + sbar - mf2(4) )                          &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffLog

   sum2 = - mTgT * diff                                              &
      & + ( mTgT2 - ( mT2 + sbar - mf2(4) )                          &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffTan   &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(3) - mf2(4) )   &
      &   * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog 

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) ) * DiffLog
   sum2 = - ( mT2 - mf2(1) ) * DiffTan -  mTgT  * DiffLog

   erg(5) = erg(5) + ReProp * sum1 + ImProp * sum2 
   erg(6) = erg(6) + ReProp * sum2 - ImProp * sum1 
 
   sum1 = DiffLog * (sbar - mf2(3))
   sum2 = DiffTan * (sbar - mf2(3))

   erg(7) = erg(7) + ReProp * sum1 + ImProp * sum2 
   erg(8) = erg(8) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel2


 Subroutine F3BDgaugeSTkernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: sbar,DiffLog,DiffTan,Prop,ReProp,ImProp
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m22 = mf2(2)
   sumI = sumI - (mf2(1) - mf2(4)) * mf2(2) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(2) + mf2(4) - 2._dp * sbar         &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)     &
      & - ( mf2(1) + mf2(2) + mf2(4) - 2._dp * sbar         &
      &   - 1.5_dp * mT2- 0.5_dp * tmin) * (tmin - mT2)      &
      & + mTgT * ( 2._dp * (mT2 + sbar)                     &
      &          -  mf2(1) - mf2(2) - mf2(4) ) * DiffTan   &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) - mf2(4) )       &
      &             * ( mT2 + sbar - mf2(1) ) ) * DiffLog

   sum2 = - mTgT * diff                                    &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) - mf2(4) )       &
      &             * ( mT2 + sbar - mf2(1) ) ) * DiffTan  &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1)             &
      &   - mf2(2) - mf2(4) ) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(4) ) * DiffLog
   sum2 = ( mT2 - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(5) = erg(5) + ReProp * sum1 + ImProp * sum2 
   erg(6) = erg(6) + ReProp * sum2 - ImProp * sum1 

   sum1 = DiffLog * (sbar - mf2(2))
   sum2 = DiffTan * (sbar - mf2(2))

   erg(7) = erg(7) + ReProp * sum1 + ImProp * sum2 
   erg(8) = erg(8) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel3


 Subroutine F3BDgaugeSTkernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: sbar,DiffLog,DiffTan,Prop,ReProp,ImProp
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + mf2(1) * (mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(3) + mf2(2) - 2._dp * sbar                   &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)               &
      & - ( mf2(1) + mf2(3) + mf2(2) - 2._dp * sbar                   &
      &   - 1.5_dp * mT2- 0.5_dp * tmin) * (tmin - mT2)                &
      & + mTgT * ( 2._dp * (mT2 + sbar)                               &
      &          -  mf2(1) - mf2(3) - mf2(2) ) * DiffTan             &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) )                          &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffLog

   sum2 = - mTgT * diff                                              &
      & + ( mTgT2 - ( mT2 + sbar - mf2(2) )                          &
      &             * ( mT2 + sbar - mf2(1) - mf2(3) ) ) * DiffTan   &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(3) - mf2(2) )   &
      &   * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog 

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(3)) * DiffLog
   sum2 = ( mT2 - mf2(3) - mf2(4) ) * DiffTan + mTgT * DiffLog

   erg(5) = erg(5) + ReProp * sum1 + ImProp * sum2 
   erg(6) = erg(6) + ReProp * sum2 - ImProp * sum1 

   sum1 = DiffLog * (-sbar + mf2(1))
   sum2 = DiffTan * (-sbar + mf2(1))

   erg(7) = erg(7) + ReProp * sum1 + ImProp * sum2 
   erg(8) = erg(8) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel4


 Subroutine F3BDgaugeSTkernel5(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: sbar, DiffLog, DiffTan, Prop, ReProp, ImProp
  Real(dp) :: m12, m22, sum1, sum2, tmin, tmax, sumI, diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(4) - 2._dp * sbar                                 &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)                    &
      & - ( mf2(1) + mf2(4) - 2._dp * sbar                                 & 
      &   - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2)                    &
      & + mTgT * ( 2._dp * (mT2 + sbar) -  mf2(1) - mf2(4) ) * DiffTan     &
      & + (mTgT2 - ( mT2 + sbar - mf2(4)) * (mT2 + sbar - mf2(1)) ) * DiffLog

   sum2 = - mTgT * diff                                       &
      & + ( mTgT2 - ( mT2 + sbar - mf2(4) )                   &
      &             * ( mT2 + sbar - mf2(1) ) ) * DiffTan     &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(4) ) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = DiffLog * sbar
   sum2 = DiffTan * sbar

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel5


 Subroutine F3BDgaugeSTkernel6(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: sbar, DiffLog, DiffTan, Prop, ReProp, ImProp
  Real(dp) :: m12, sum1, sum2, tmin, tmax, sumI, diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   sumI = sumI + mf2(1) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(3) - 2._dp * sbar                                    &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)                       &
      & - ( mf2(1) + mf2(3) - 2._dp * sbar                                    &
      &   - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2)                       &
      & + mTgT * ( 2._dp * (mT2 + sbar) -  mf2(1) - mf2(3) ) * DiffTan        &
      & + ( mTgT2 - (mT2 + sbar) * (mT2 + sbar - mf2(1) - mf2(3)) ) * DiffLog

   sum2 = - mTgT * diff                                                      &
      & + (mTgT2 - (mT2 + sbar) * (mT2 + sbar - mf2(1) - mf2(3)) ) * DiffTan &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(3) ) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog 

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel6


 Subroutine F3BDgaugeSTkernel7(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: sbar, DiffLog, DiffTan, Prop, ReProp, ImProp
  Real(dp) :: m12, sum1, sum2, tmin, tmax, sumI, diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(2)
   sumI = sumI - mf2(1) * mf2(2) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = ( mf2(1) + mf2(2) - 2._dp * sbar                                    &
      &   - 1.5_dp * mT2 - 0.5_dp * tmax) * (tmax - mT2)                       &
      & - ( mf2(1) + mf2(2) - 2._dp * sbar                                    &
      &   - 1.5_dp * mT2 - 0.5_dp * tmin) * (tmin - mT2)                       &
      & + mTgT * ( 2._dp * (mT2 + sbar) -  mf2(1) - mf2(2) ) * DiffTan        &
      & + (mTgT2 - (mT2 + sbar - mf2(2)) * (mT2 + sbar - mf2(1)) ) * DiffLog

   sum2 = - mTgT * diff                                                      &
      & + (mTgT2 - (mT2 + sbar - mf2(2)) * (mT2 + sbar - mf2(1)) ) * DiffTan &
      & - (2._dp * mT2 + 2._dp * sbar -  mf2(1) - mf2(2) ) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

   sum1 = diff - mTgT * DiffTan + mT2 * DiffLog
   sum2 = mT2 * DiffTan + mTgT * DiffLog

   erg(3) = erg(3) + ReProp * sum1 + ImProp * sum2 
   erg(4) = erg(4) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel7


 Subroutine F3BDgaugeSTkernel8(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeST for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: sbar, DiffLog, DiffTan, Prop, ReProp, ImProp
  Real(dp) :: sum1, sum2, tmax, sumI, diff

  erg = 0._dp

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) - sbar
   diff = Abs(sbar-mf2(1))

   tmax = sumI

   If (mTgT.Eq.0._dp) Then
    DiffLog = 0.5_dp * Log( (tmax - mT2)**2/ mT2**2)
    DiffTan = 0._dp
   Else
    DiffLog = 0.5_dp * Log( ( (tmax - mT2)**2 + mTgT2 ) / (mT2**2 + mTgT2 ) )
    DiffTan = Atan( (tmax-mT2) / mTgT )  - Atan( (-mT2) / mTgT )
   Endif

   Prop = 1._dp / ( (sbar-mS2)**2 + mSgS2 )
   ReProp = (sbar-mS2) * Prop
   ImProp = mSgS * Prop

   sum1 = (mf2(1) - 2._dp * sbar - 1.5_dp * mT2 - 0.5_dp * tmax)*(tmax - mT2) &
      & - (mf2(1) - 2._dp * sbar - 1.5_dp * mT2 ) * ( - mT2) &
      & + mTgT * ( 2._dp * (mT2 + sbar) - mf2(1) ) * DiffTan                  &
      & + ( mTgT2 - (mT2 + sbar) * (mT2 + sbar - mf2(1)) ) * DiffLog

   sum2 = - mTgT * diff                                                      &
      & + ( mTgT2 - (mT2 + sbar) * (mT2 + sbar - mf2(1)) ) * DiffTan         &
      & - ( 2._dp * mT2 + 2._dp * sbar -  mf2(1)) * mTgT * DiffLog

   erg(1) = erg(1) + ReProp * sum1 + ImProp * sum2 
   erg(2) = erg(2) + ReProp * sum2 - ImProp * sum1 

  Enddo

 End Subroutine F3BDgaugeSTkernel8



 Subroutine F3BDscalarS1S2int(scalar,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the integral of the matrix element squared of
 ! M^2_fi = Tr[(C[1] Pl + C[2] Pr).(P1+mf[1]).(C[3] Pl + C[4] Pr).(P2+mf[2])]
 !      Tr[(C[5] Pl + C[6] Pr).(P3+mf[3]).(C[7] Pl + C[8] Pr).(P4+mf[4])] /
 !      ( (p1-p2)^2 - m^2_S1 + I m_S1 \[CapitalGamma]_S1 )
 !        (p1-p2)^2 - m^2_S2 - I m_S2 \[CapitalGamma]_S2 ) )\n
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p2)^2, t=(p1-p3)^2, and u=(p1-p4)^2.
 ! Here {C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8]} are the couplings,
 ! {mf[1],mf[2],mf[3],mf[4]} are the fermion masses, and
 ! {m_S1,\[CapitalGamma]_S1,mS2,\[CapitalGamma]_S2} are the mass and the total 
 ! decay widths of the scalar bosons.
 ! written by Werner Porod, 16.2.2000
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: scalar(4),mf(4),eps
  Complex(dp), Intent(in) :: coup(8)
  Complex(dp), Intent(inout) :: int1(4)
  Complex(dp), Intent(out) :: erg
  Logical, Intent(inout) :: integrate

  Real(dp) :: int1a(8),int1b(4),int1c(2),smin,smax,sG(4)

  Integer :: i1,len1,Imin,Imax
  Complex(dp) :: coupC(4)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   erg = ZeroC
   Return
  Endif
  
  mf2 = mf**2

  erg = ZeroC

  If (Integrate) Then
   If (scalar(1).Lt.scalar(3)) Then
    sG(1) = (scalar(1)-2._dp*scalar(2))**2
    sG(2) = (scalar(1)+2._dp*scalar(2))**2
    sG(3) = (scalar(3)-2._dp*scalar(4))**2
    sG(4) = (scalar(3)+2._dp*scalar(4))**2
   Else
    sG(1) = (scalar(3)-2._dp*scalar(4))**2
    sG(2) = (scalar(3)+2._dp*scalar(4))**2
    sG(3) = (scalar(1)-2._dp*scalar(2))**2
    sG(4) = (scalar(1)+2._dp*scalar(2))**2
   Endif
   If (sG(2).Ge.sG(3)) Then
    sG(1) = Min( sG(1),sG(3) )
    sG(2) = Max( sG(2),sG(4) )
    len1 = 3
   Else
    len1 = 5
   Endif
   Imin = len1
   Imax = 0
  Endif
  
  If ( (mf(2).Eq.0._dp).And.((mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp))) Then
   coupC(1) = coup(1) * coup(4) + coup(2) * coup(3)
   coupC(3) = coup(5) * coup(8) + coup(6) * coup(7)
   If ((coupC(1).Eq.ZeroC).Or.(coupC(3).Eq.ZeroC)) Then
    Integrate = .False.
    Return
   Endif

   If (Integrate) Then
    mS12 = scalar(1)**2
    mS1gS1 = scalar(1) * scalar(2)
    mS22 = scalar(3)**2
    mS2gS2 = scalar(3) * scalar(4)
    smax = mf2(1)
    smin = (Abs(mf(3)) + Abs(mf(4)))**2

    int1(1) = 0._dp
    int1(2) = 0._dp
    int1(3) = 0._dp
    Do i1=1,len1-1
     If (smin.Lt.sG(len1-i1)) Imin = len1 - i1
     If (smax.Gt.sG(i1)) Imax = i1
    Enddo
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarS1S2kernel4,2,smin,smax,int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarS1S2kernel4,2,smin,sG(Imin),int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
     Call DgaussInt(F3BDscalarS1S2kernel4,2,sG(Imin),smax,int1c,eps)
     int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
    Else
     Call DgaussInt(F3BDscalarS1S2kernel4,2,smin,sG(Imin),int1c,eps)
     int1(4) = int1c(1) + Ic * int1c(2)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDscalarS1S2kernel4,2,sG(i1),sG(i1+1),int1c,eps)
      int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
     Enddo
     Call DgaussInt(F3BDscalarS1S2kernel4,2,sG(Imax),smax,int1c,eps)
     int1(4) = int1(4) + int1c(1) + Ic * int1c(2)
    Endif
   Endif
   erg = coupC(1) * coupC(3) * int1(4)

  Elseif ( (mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp)) Then
   coupC(1) = coup(1) * coup(4) + coup(2) * coup(3)
   coupC(2) = coup(2) * coup(4) + coup(1) * coup(3) 
   coupC(3) = coup(5) * coup(8) + coup(6) * coup(7)
   If ( (coupC(3).Eq.ZeroC).Or.                         &
    &   ((coupC(1).Eq.ZeroC).And.(coupC(2).Eq.ZeroC)) ) Then
    Integrate = .False.
    Return
   Endif
   If (Integrate) Then
    mS12 = scalar(1)**2
    mS1gS1 = scalar(1) * scalar(2)
    mS22 = scalar(3)**2
    mS2gS2 = scalar(3) * scalar(4)
    smax = (Abs(mf(1)) - Abs(mf(2)))**2
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    int1(1) = 0._dp
    int1(3) = 0._dp
    Do i1=1,len1-1
     If (smin.Lt.sG(len1-i1)) Imin = len1 - i1
     If (smax.Gt.sG(i1)) Imax = i1
    Enddo
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarS1S2kernel3,4,smin,smax,int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarS1S2kernel3,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Call DgaussInt(F3BDscalarS1S2kernel3,4,sG(Imin),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Else
     Call DgaussInt(F3BDscalarS1S2kernel3,4,smin,sG(Imin),int1b,eps)
     int1(2) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDscalarS1S2kernel3,4,sG(i1),sG(i1+1),int1b,eps)
      int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
      int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
     Enddo
     Call DgaussInt(F3BDscalarS1S2kernel3,4,sG(Imax),smax,int1b,eps)
     int1(2) = int1(2) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Endif
   Endif
   erg = coupC(1) * coupC(3) * int1(4)      &
     & + 2._dp * coupC(2) * coupC(3) * mf(1) * mf(2) * int1(2)

  Elseif (mf(2).Eq.0._dp) Then
   coupC(1) = coup(1) * coup(4) + coup(2) * coup(3)
   coupC(3) = coup(5) * coup(8) + coup(6) * coup(7)
   coupC(4) = coup(6) * coup(8) + coup(5) * coup(7) 
   If ( (coupC(1).Eq.ZeroC).Or.                       &
      & ((coupC(3).Eq.ZeroC).And.(coupC(4).Eq.ZeroC)) ) Then
    Integrate = .False.
    Return
   Endif
   If (Integrate) Then
    mS12 = scalar(1)**2
    mS1gS1 = scalar(1) * scalar(2)
    mS22 = scalar(3)**2
    mS2gS2 = scalar(3) * scalar(4)
    smax = mf2(1)
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    int1(1) = 0._dp
    int1(2) = 0._dp
    Do i1=1,len1-1
     If (smin.Lt.sG(len1-i1)) Imin = len1 - i1
     If (smax.Gt.sG(i1)) Imax = i1
    Enddo
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarS1S2kernel2,4,smin,smax,int1b,eps)
     int1(3) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarS1S2kernel2,4,smin,sG(Imin),int1b,eps)
     int1(3) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Call DgaussInt(F3BDscalarS1S2kernel2,4,sG(Imin),smax,int1b,eps)
     int1(3) = int1(3) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Else
     Call DgaussInt(F3BDscalarS1S2kernel2,4,smin,sG(Imin),int1b,eps)
     int1(3) = int1b(1) + Ic * int1b(2)
     int1(4) = int1b(3) + Ic * int1b(4)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDscalarS1S2kernel2,4,sG(i1),sG(i1+1),int1b,eps)
      int1(3) = int1(3) + int1b(1) + Ic * int1b(2)
      int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
     Enddo
     Call DgaussInt(F3BDscalarS1S2kernel2,4,sG(Imax),smax,int1b,eps)
     int1(3) = int1(3) + int1b(1) + Ic * int1b(2)
     int1(4) = int1(4) + int1b(3) + Ic * int1b(4)
    Endif

   Endif
   erg = coupC(1) * coupC(3) * int1(4)      &
     & + 2._dp * coupC(1) * coupC(4) * mf(3) * mf(4) * int1(3)

  Else
   coupC(1) = coup(1) * coup(4) + coup(2) * coup(3)
   coupC(2) = coup(2) * coup(4) + coup(1) * coup(3) 
   coupC(3) = coup(5) * coup(8) + coup(6) * coup(7)
   coupC(4) = coup(6) * coup(8) + coup(5) * coup(7) 
   If ( ((coupC(1).Eq.ZeroC).And.(coupC(2).Eq.ZeroC)) .Or. &
      & ((coupC(3).Eq.ZeroC).And.(coupC(4).Eq.ZeroC))     ) Then
    Integrate = .False.
    Return
   Endif
   If (Integrate) Then
    mS12 = scalar(1)**2
    mS1gS1 = scalar(1) * scalar(2)
    mS22 = scalar(3)**2
    mS2gS2 = scalar(3) * scalar(4)
    smax = (Abs(mf(1)) - Abs(mf(2)))**2
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    Do i1=1,len1-1
     If (smin.Lt.sG(len1-i1)) Imin = len1 - i1
     If (smax.Gt.sG(i1)) Imax = i1
    Enddo
    If ((Imin.Eq.len1).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarS1S2kernel1,8,smin,smax,int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarS1S2kernel1,8,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
     Call DgaussInt(F3BDscalarS1S2kernel1,8,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)
    Else
     Call DgaussInt(F3BDscalarS1S2kernel1,8,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:7:2) + Ic * int1a(2:8:2)
     Do i1=Imin,Imax-1
      Call DgaussInt(F3BDscalarS1S2kernel1,8,sG(i1),sG(i1+1),int1a,eps)
      int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)
     Enddo
     Call DgaussInt(F3BDscalarS1S2kernel1,8,sG(Imax),smax,int1a,eps)
     int1 = int1 + int1a(1:7:2) + Ic * int1a(2:8:2)
    Endif
   Endif
   erg = coupC(1) * coupC(3) * int1(4)                             &
     & + 2._dp * ( coupC(1) * coupC(4) * mf(3) * mf(4) * int1(3)    &
     &          + coupC(2) * coupC(3) * mf(1) * mf(2) * int1(2) )  &
     & + 4._dp * coupC(2) * coupC(4) * mf(1) * mf(2) * mf(3) * mf(4) * int1(1)

  Endif

 End Subroutine F3BDscalarS1S2int


 Subroutine F3BDscalarS1S2kernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar,multA,multB

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sumI = kappa(sbar,m12,m22)

   If ((mf2(3).Ne.0._dp).Or.(mf2(4).Ne.0._dp)) Then
    m12 = mf2(3)
    m22 = mf2(4)
    sumI = sumI * kappa(sbar,m12,m22) / sbar
   Endif

   sumI = sumI / (((sbar-mS12)**2 + mS1gS1**2) * ((sbar-mS22)**2 + mS2gS2**2))

   multA = sumI * ( (sbar-mS12) * (sbar-mS22) + mS1gS1 * mS2gS2)
   multB = sumI * ( (sbar-mS12) * mS2gS2 - (sbar-mS22) * mS1gS1)
   erg(1) = erg(1) + multA
   erg(2) = erg(2) + multB
   erg(3) = erg(3) + multA * (sbar - mf2(3) - mf2(4))
   erg(4) = erg(4) + multB * (sbar - mf2(3) - mf2(4))
   erg(5) = erg(5) + multA * (mf2(1) + mf2(2) - sbar)
   erg(6) = erg(6) + multB * (mf2(1) + mf2(2) - sbar)
   erg(7) = erg(7) + multA *(mf2(1) + mf2(2) -sbar) *(sbar - mf2(3) - mf2(4))
   erg(8) = erg(8) + multB *(mf2(1) + mf2(2) -sbar) *(sbar - mf2(3) - mf2(4))
  Enddo

 End Subroutine F3BDscalarS1S2kernel1


 Subroutine F3BDscalarS1S2kernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
  Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar,multA,multB

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   sumI = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI * kappa(sbar,m12,m22) / sbar

   sumI = sumI * (mf2(1) - sbar) &
        &  / (((sbar-mS12)**2 + mS1gS1**2) * ((sbar-mS22)**2 + mS2gS2**2))

   multA = sumI * ( (sbar-mS12) * (sbar-mS22) + mS1gS1 * mS2gS2)
   multB = sumI * ( (sbar-mS12) * mS2gS2 - (sbar-mS22) * mS1gS1)
   erg(1) = erg(1) + multA
   erg(2) = erg(2) + multB
   erg(3) = erg(3) + multA * (sbar - mf2(3) - mf2(4))
   erg(4) = erg(4) + multB * (sbar - mf2(3) - mf2(4))
  Enddo

 End Subroutine F3BDscalarS1S2kernel2


 Subroutine F3BDscalarS1S2kernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar,multA,multB

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   m12 = mf2(1)
   m22 = mf2(2)
   sumI = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI * kappa(sbar,m12,m22) / sbar

   sumI = sumI  * (sbar - mf2(3) - mf2(4)) &
        & / (((sbar-mS12)**2 + mS1gS1**2) * ((sbar-mS22)**2 + mS2gS2**2))

   multA = sumI * ( (sbar-mS12) * (sbar-mS22) + mS1gS1 * mS2gS2)
   multB = sumI * ( (sbar-mS12) * mS2gS2 - (sbar-mS22) * mS1gS1)
   erg(1) = erg(1) + multA
   erg(2) = erg(2) + multB
   erg(3) = erg(3) + multA *(mf2(1) + mf2(2) -sbar)
   erg(4) = erg(4) + multB *(mf2(1) + mf2(2) -sbar)
  Enddo

 End Subroutine F3BDscalarS1S2kernel3


 Subroutine F3BDscalarS1S2kernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: sumI, sbar

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = Abs(sbar-mf2(1))

   If ((mf2(3).Ne.0._dp).And.(mf2(4).Ne.0._dp)) Then
    sumI = sumI * kappa(sbar,mf2(3),mf2(4)) / sbar
   Else If (mf2(3).Ne.0._dp) Then
    SumI = SumI * Abs(1._dp - mf2(3) / sbar )
   Else If (mf2(4).Ne.0._dp) Then
    SumI = SumI * Abs(1._dp - mf2(4) / sbar )
   Endif

   sumI = sumI  * (mf2(1) - sbar) * (sbar - mf2(3) - mf2(4)) &
        & / (((sbar-mS12)**2 + mS1gS1**2) * ((sbar-mS22)**2 + mS2gS2**2 ))

   erg(1) = erg(1) + sumI * ( (sbar-mS12) * (sbar-mS22) + mS1gS1 * mS2gS2)
   erg(2) = erg(2) + sumI * ( (sbar-mS12) * mS2gS2 - (sbar-mS22) * mS1gS1)
  Enddo

 End Subroutine F3BDscalarS1S2kernel4

 Complex(dp) Function F3BDscalarS1S2kernel5(mi2, mgS)
 Implicit None
  Real(dp), Intent(in) :: mi2, mgS(4)

  Real(dp) :: ratio, mS12, mS22, eps, Log1, mGS1, mGS2, nen, sum1, sum2, sum3 &
    & , sum4, Atan1, Atan2, Log2

  F3BDscalarS1S2kernel5 = 0._dp

  mS12 = mgS(1)**2
  mS22 = mgS(3)**2
  If ((mgS(2).Eq.0._dp).And.(mgS(4).Eq.0._dp)) Then ! 0 widths
   ratio = mS12 / mS22
   If (Abs(ratio-1._dp).Lt.1.e-6_dp) Then
    eps = mS22 - mS12
    Log1 = Log(1._dp - mi2/mS12) 
    F3BDscalarS1S2kernel5 = mi2 * (3._dp * mS12 - 2.5_dp * mi2)     &
      & + (3._dp * mS12**2 - 4._dp * mS12 * mi2 + mi2**2) * Log1    &
      & + ( mi2 * (6._dp * mS12 - mi2)                              &
      &   + mS12 * (6._dp * mS12 - 4._dp * mi2) * Log1              &
      &   ) * 0.5_dp * eps / mS12                                   &
      & + ( mi2 * (mi2**2 + 3._dp * mi2 * mS12 - 6._dp * mS12**2)   &
      &   + 6._dp * mS12**2 * (mi2 - mS12) * Log1                   &
      &   ) * eps**2 / (6._dp * mS12**2 * (mi2-mS12) )              &
      & - mi2**4 * eps**3 / (12._dp ** mS12**3 * (mi2-mS12)**2 )    &
      & + mi2**4 * (3._dp * mi2 - 5._dp * mS12) * eps**4            &
      &    / (60._dp * mS12**4 * (mi2-mS12)**3 ) 
   Else
    F3BDscalarS1S2kernel5 = mi2 * (mS12 + mS22 - 1.5_dp * mi2)     &
         &  + ( mS12 * (mi2 - mS12)**2 * Log( 1._dp - mi2 / mS12)  &
         &    - mS22 * (mi2 - mS22)**2 * Log( 1._dp - mi2 / mS22)  &
         &    ) / (mS12-mS22)
   End If
 
  Else ! including decay widths
    mGS1 = mgS(1) * mgS(2)
    mGS2 = mgS(3) * mgS(4)

    nen = 1._dp / ( (mS12 - mS22)**2 + (mGS1 + mGS2)**2)
    sum1 = mGS1**3 * (mGS1 + mGS2)                                            &
       & - mGS1 * mGS2 * (3._dp * mS12 - mi2) * (mS12 - mi2)                  &
       & - mGS1*2 * (mi2**2 + 3._dp * mS12 * mS22 - 2._dp*mi2 * (mS12+mS22) ) &
       & - mS12 * (mS12-mS22) * (mi2-mS12)**2
    sum2 = mGS2**3 * (mGS1 + mGS2)                                            &
       & - mGS1 * mGS2 * (3._dp * mS22 - mi2) * (mS22 - mi2)                  &
       & - mGS2*2 * (mi2**2 + 3._dp * mS12 * mS22 - 2._dp*mi2 * (mS12+mS22) ) &
       & - mS22 * (mS22-mS12) * (mi2-mS22)**2
    sum3 = mGS1**2 * mGS2 * (2._dp * mi2 - 3._dp * mS12)                   &
       & + mGS2 * mS12 * (mi2 - mS12)**2                                   &
       & + mGS1**3 * (2._dp * (mi2-mS12) - mS22)                           &
       & + mGS1 * (mi2-mS12) * (2._dp* mS12 *(mS12-mS22) + (mi2-mS12) * mS22 )
    sum4 = mGS2**2 * mGS1 * (2._dp * mi2 - 3._dp * mS22)                   &
       & + mGS1 * mS22 * (mi2 - mS22)**2                                   &
       & + mGS2**3 * (2._dp * (mi2-mS22) - mS12)                           &
       & + mGS2 * (mi2-mS22) * (2._dp* mS22 *(mS22-mS12) + (mi2-mS22) * mS12 )
    Atan1 = Atan( mGS1 * mi2 / (mGS1**2 + mS12 * (mS12 - mi2) ) )
    Atan2 = Atan( mGS2 * mi2 / (mGS2**2 + mS22 * (mS22 - mi2) ) )

    If (mi2.Gt.mS12) Atan1 = Atan1 + Pi
    If (mi2.Gt.mS22) Atan2 = Atan2 + Pi
    Log1 = 0.5_dp * Log( ((mi2-mS12)**2 + mGS1**2) / (mS12**2 + mGS1**2 ) )
    Log2 = 0.5_dp * Log( ((mi2-mS22)**2 + mGS2**2) / (mS22**2 + mGS2**2 ) )

    F3BDscalarS1S2kernel5 = mi2 * (mS12 + mS22 - 1.5_dp * mi2)     &
       &     + nen * (- sum1 * Log1 - sum2 * Log2                  &
       &             + sum3 * Atan1 + sum4 * Atan2)                &
       &   + Ic * ( mi2 * (mGS2 - mGS1)                            &
       &          + nen * ( sum3 * Log1 + sum4 * Log2              &
       &                  - sum1 * Atan1 + sum2 * Atan2) )
  End If

 End Function F3BDscalarS1S2kernel5


 Complex(dp) Function F3BDscalarS1S2kernel6(mi2, mgS)
 implicit none
  Real(dp), Intent(in) :: mi2, mgS(4)

  Real(dp) :: ratio, mS12, mS22, eps, Log1, mGS1, mGS2, Atan1, Atan2, Log2
  Real(dp) :: re, im, nenR, nenI, sum1R, sum1I, sum2R, sum2I 

  F3BDscalarS1S2kernel6 = 0._dp

  mS12 = mgS(1)**2
  mS22 = mgS(3)**2
  if ((mgS(2).eq.0._dp).and.(mgS(4).eq.0._dp)) then ! 0 widths
   ratio = mS12 / mS22
   If (Abs(ratio-1._dp).Lt.1.e-6_dp) then
    eps = mS22 - mS12
    Log1 = Log(1._dp - mi2/mS12) 
    F3BDscalarS1S2kernel6 = mi2 * (3._dp * mS12 - 2.5_dp * mi2)     &
      & + (3._dp * mS12**2 - 4._dp * mS12 * mi2 + mi2**2) * Log1    &
      & + ( mi2 * (6._dp * mS12 - mi2)                              &
      &   + mS12 * (6._dp * mS12 - 4._dp * mi2) * Log1              &
      &   ) * 0.5_dp * eps / mS12                                   &
      & + ( mi2 * (mi2**2 + 3._dp * mi2 * mS12 - 6._dp * mS12**2)   &
      &   + 6._dp * mS12**2 * (mi2 - mS12) * Log1                   &
      &   ) * eps**2 / (6._dp * mS12**2 * (mi2-mS12) )              &
      & - mi2**4 * eps**3 / (12._dp ** mS12**3 * (mi2-mS12)**2 )    &
      & + mi2**4 * (3._dp * mi2 - 5._dp * mS12) * eps**4            &
      &    / (60._dp * mS12**4 * (mi2-mS12)**3 ) 
   else
    F3BDscalarS1S2kernel6 = mi2 * (mS12 + mS22 - 1.5_dp * mi2)     &
         &  + ( mS12 * (mi2 - mS12)**2 * Log( 1._dp - mi2 / mS12)  &
         &    - mS22 * (mi2 - mS22)**2 * Log( 1._dp - mi2 / mS22)  &
         &    ) / (mS12-mS22)
   end if
 
  else ! including decay widths
    mGS1 = mgS(1) * mgS(2)
    mGS2 = mgS(3) * mgS(4)

    re = mi2 * (mS12 + mS22 - 1.5*mi2)
    im =  mi2 *(mGS1-mGS2)

    nenR = 0.5_dp * (mGS1 + mGS2) / ( (mGS1 + mGS2)**2 + (mS12-mS22)**2 )
    nenI = 0.5_dp * (mS12-mS22) /  ( (mGS1 + mGS2)**2 + (mS12-mS22)**2 )

    sum1R = - mGS1 * (mGS1**2 + (mi2-mS12) * (mi2-3._dp*mS12) )
    sum1I = (mi2-mS12) * (mS12 * (mi2-mS12) + 2._dp * mGS1**2) - mS12 * mGS1**2
    sum2R = - mGS2 * (mGS2**2 + (mi2-mS22) * (mi2-3._dp*mS22) )
    sum2I = (mi2-mS22) * (mS22 * (mi2-mS22) + 2._dp * mGS2**2) - mS22 * mGS2**2

    Log1 = Log( ((mi2-mS12)**2 + mGS1**2) / (mS12**2 + mGS1**2 ) )
    Log2 = Log( ((mi2-mS22)**2 + mGS2**2) / (mS22**2 + mGS2**2 ) )
    Atan1 = 2._dp * Atan( -mGS1 * mi2 / (mGS1**2 + mS12 * (mS12 - mi2) ) )
    Atan2 = 2._dp * Atan( -mGS2 * mi2 / (mGS2**2 + mS22 * (mS22 - mi2) ) )
    If (mi2.Gt.mS12) Atan1 = Atan1 + Pi
    If (mi2.Gt.mS22) Atan2 = Atan2 + Pi

    re = re + (nenR * sum1R - nenI * sum1I) * Log1         &
       &    - (nenI * sum1R + nenR * sum1I) * Atan1        &
       &    + (nenR * sum2R - nenI * sum2I) * Log2         &
       &    - (nenI * sum2R + nenR * sum2I) * Atan2
    F3BDscalarS1S2kernel6 = re + ic * im

  end if

 end Function F3BDscalarS1S2kernel6


 Real(dp) Function F3BDscalarSSa4(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s

  Real(dp) :: sumI

  sumI = kappa(s,mf2(1),mf2(2))

  If ((mf2(3).Ne.0._dp).And.(mf2(4).Ne.0._dp)) Then
   sumI = sumI * kappa(s,mf2(3),mf2(4)) / s
  Else If (mf2(3).Ne.0._dp) Then
   SumI = SumI * Abs(1._dp - mf2(3) / s )
  Else If (mf2(4).Ne.0._dp) Then
   SumI = SumI * Abs(1._dp - mf2(4) / s )
  Endif

  F3BDscalarSSa4 = sumI * (mf2(1) + mf2(2) -s)  &
               & * (s - mf2(3) - mf2(4)) / ( (s-mS2)**2 + mSG2 )

 End Function F3BDscalarSSa4


 Subroutine F3BDscalarSSint(scalar,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the integral of the matrix element squared of
 ! M^2_fi = Tr[(C[1] Pl + C[2] Pr).(P1+mf[1]).(C[3] Pl + C[4] Pr).(P2+mf[2])]
 !      Tr[(C[5] Pl + C[6] Pr).(P3+mf[3]).(C[7] Pl + C[8] Pr).(P4+mf[4])] /
 !    / ((p1-p2)^2 - m^2_S1)^2 + m_S1^2 \[CapitalGamma]_S1^2 )
 ! in terms of the Mandelstam variable s which are given by s=(p1-p2)^2.
 ! Here {C[1],C[2],C[3],C[4]} are the couplings, Cc are the complex conjugated
 ! couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses, and
 ! {m_S1,\[CapitalGamma]_S1} are the mass and the total decay widths of the
 ! scalar boson.
 ! written by Werner Porod, 17.2.2000
 ! changing the integraion from dgauss to DgaussInt: 14.6.2000
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: scalar(2),mf(4),eps
  Real(dp), Intent(inout) :: int1(4)
  Real(dp), Intent(out) :: erg
  Complex(dp), Intent(in) :: coup(4)
  Logical, Intent(inout) ::integrate

  Real(dp) :: smin,smax,sminG,smaxG,int2a(2),int2b(2) &
               & , int2c(2),mr,resR, coupC(4)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   int1 = 0._dp
   erg = 0._dp
   Return
  Endif

  mf2 = mf**2

  erg = 0._dp
  coupC(1) = Abs( coup(1) )**2 + Abs( coup(2) )**2
  coupC(3) = Abs( coup(3) )**2 + Abs( coup(4) )**2
  If ( (coupC(1).Eq.0._dp).Or.(coupC(3).Eq.0._dp) ) Then
   If (Integrate) Integrate = .False.
   Return
  Endif

  If ( (mf(2).Eq.0._dp).And.((mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp))) Then

   If (Integrate) Then
    int1 = 0._dp
    If ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
     mR = mf(1)
     Call IntScalarSS1(scalar,mR,resR)
     int1(4) = resR
    Else
     mS2 = scalar(1)**2
     mSG2 = mS2 * scalar(2)**2
     smax = mf2(1)
     smin = (Abs(mf(3)) + Abs(mf(4)))**2
     sminG = (scalar(1)-5._dp*scalar(2))**2
     smaxG = (scalar(1)+5._dp*scalar(2))**2
     If ((smin.Gt.mS2).Or.(smax.Lt.mS2)) Then
      int1(4) = dgauss(F3BDscalarSSa4,smin,smax,eps)
     Elseif ((smin.Lt.sminG).And.(smax.Gt.smaxG)) Then
      int1(4) = dgauss(F3BDscalarSSa4,smin,sminG,eps)  &
            & + dgauss(F3BDscalarSSa4,sminG,smaxG,eps) &
            & + dgauss(F3BDscalarSSa4,smaxG,smax,eps)
     Elseif (smin.Lt.sminG) Then
      int1(4) = dgauss(F3BDscalarSSa4,smin,sminG,eps)  &
            & + dgauss(F3BDscalarSSa4,sminG,smax,eps)
     Else
      int1(4) = dgauss(F3BDscalarSSa4,smin,smaxG,eps)  &
            & + dgauss(F3BDscalarSSa4,smaxG,smax,eps)
     Endif 
    Endif
   Endif
   erg = coupC(1) * coupC(3) * int1(4)

  Elseif ( (mf(3).Eq.0._dp).Or.(mf(4).Eq.0._dp)) Then
   coupC(2) = Real(coup(1) * Conjg( coup(2) ),dp )
   If (Integrate) Then
    mS2 = scalar(1)**2
    mSG2 = mS2 * scalar(2)**2
    smax = (Abs(mf(1)) - Abs(mf(2)))**2
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    sminG = (scalar(1)-5._dp*scalar(2))**2
    smaxG = (scalar(1)+5._dp*scalar(2))**2
    int1(1) = 0._dp
    int1(3) = 0._dp
    If ((smin.Gt.mS2).Or.(smax.Lt.mS2)) Then ! resonance is outside
     Call DgaussInt(F3BDscalarSSkernel3,2,smin,smax,int2a,eps)
     int1(2) = int2a(1)
     int1(4) = int2a(2)
    Elseif ((smin.Lt.sminG).And.(smax.Gt.smaxG)) Then
     Call DgaussInt(F3BDscalarSSkernel3,2,smin,sminG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel3,2,sminG,smaxG,int2b,eps)
     Call DgaussInt(F3BDscalarSSkernel3,2,smaxG,smax,int2c,eps)
     int1(2) = int2a(1)+int2b(1)+int2c(1)
     int1(4) = int2a(2)+int2b(2)+int2c(2)
    Elseif (smin.Lt.sminG) Then
     Call DgaussInt(F3BDscalarSSkernel3,2,smin,sminG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel3,2,sminG,smax,int2b,eps)
     int1(2) = int2a(1)+int2b(1)
     int1(4) = int2a(2)+int2b(2)
    Elseif (smaxG.Lt.smax) Then
     Call DgaussInt(F3BDscalarSSkernel3,2,smin,smaxG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel3,2,smaxG,smax,int2b,eps)
     int1(2) = int2a(1)+int2b(1)
     int1(4) = int2a(2)+int2b(2)
    Else
     Call DgaussInt(F3BDscalarSSkernel3,2,smin,smax,int2a,eps)
     int1(2) = int2a(1)
     int1(4) = int2a(2)
    Endif 
   Endif
   erg = coupC(1) * coupC(3) * int1(4)  &
     & + 4._dp * coupC(2) * coupC(3) * mf(1) * mf(2) * int1(2)

  Elseif (mf(2).Eq.0._dp) Then
   coupC(4) = Real(coup(3) * Conjg( coup(4) ),dp )

   If (Integrate) Then
    mS2 = scalar(1)**2
    mSG2 = mS2 * scalar(2)**2
    smax = mf2(1)
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    sminG = (scalar(1)-5._dp*scalar(2))**2
    smaxG = (scalar(1)+5._dp*scalar(2))**2
    int1(1) = 0._dp
    int1(2) = 0._dp
    If ((smin.Gt.mS2).Or.(smax.Lt.mS2)) Then
     Call DgaussInt(F3BDscalarSSkernel2,2,smin,smax,int2a,eps)
     int1(3) = int2a(1)
     int1(4) = int2a(2)
    Elseif ((smin.Lt.sminG).And.(smax.Gt.smaxG)) Then
     Call DgaussInt(F3BDscalarSSkernel2,2,smin,sminG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel2,2,sminG,smaxG,int2b,eps)
     Call DgaussInt(F3BDscalarSSkernel2,2,smaxG,smax,int2c,eps)
     int1(3) = int2a(1)+int2b(1)+int2c(1)
     int1(4) = int2a(2)+int2b(2)+int2c(2)
    Elseif (smin.Lt.sminG) Then
     Call DgaussInt(F3BDscalarSSkernel2,2,smin,sminG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel2,2,sminG,smax,int2b,eps)
     int1(3) = int2a(1)+int2b(1)
     int1(4) = int2a(2)+int2b(2)
    Elseif (smaxG.Lt.smax) Then
     Call DgaussInt(F3BDscalarSSkernel2,2,smin,smaxG,int2a,eps)
     Call DgaussInt(F3BDscalarSSkernel2,2,smaxG,smax,int2b,eps)
     int1(3) = int2a(1)+int2b(1)
     int1(4) = int2a(2)+int2b(2)
    Else
     Call DgaussInt(F3BDscalarSSkernel2,2,smin,smax,int2a,eps)
     int1(3) = int2a(1)
     int1(4) = int2a(2)
    Endif 
   Endif
   erg = coupC(1) * coupC(3) * int1(4)    &
     & + 4._dp * coupC(1) * coupC(4) * mf(3) * mf(4) * int1(3)

  Else
   coupC(2) = Real(coup(1) * Conjg( coup(2) ),dp )
   coupC(4) = Real(coup(3) * Conjg( coup(4) ),dp )

   If (Integrate) Then
    mS2 = scalar(1)**2
    mSG2 = mS2 * scalar(2)**2
    smax = (Abs(mf(1)) - Abs(mf(2)))**2
    smin = (Abs(mf(3)) + Abs(mf(4)))**2
    sminG = (scalar(1)-5._dp*scalar(2))**2
    smaxG = (scalar(1)+5._dp*scalar(2))**2
    If ((smin.Gt.mS2).Or.(smax.Lt.mS2)) Then
     Do int_v=1,4     
      int1(int_v) = Dgauss(F3BDscalarSSkernel1a,smin,smax,eps)
     End Do
    Elseif ((smin.Lt.sminG).And.(smax.Gt.smaxG)) Then
     Do int_v=1,4     
      int1(int_v) = Dgauss(F3BDscalarSSkernel1a,smin,sminG,eps)  &
                & + Dgauss(F3BDscalarSSkernel1a,sminG,smaxG,eps)  &
                & + Dgauss(F3BDscalarSSkernel1a,smaxG,smax,eps)
     End Do
    Elseif (smin.Lt.sminG) Then
     Do int_v=1,4     
      int1(int_v) = Dgauss(F3BDscalarSSkernel1a,smin,sminG,eps)  &
                & + Dgauss(F3BDscalarSSkernel1a,sminG,smax,eps)
     End Do
    Elseif (smaxG.Lt.smax) Then
     Do int_v=1,4     
      int1(int_v) = Dgauss(F3BDscalarSSkernel1a,smin,smaxG,eps)  &
                & + Dgauss(F3BDscalarSSkernel1a,smaxG,smax,eps)
     End Do
    Else
     Do int_v=1,4     
      int1(int_v) = Dgauss(F3BDscalarSSkernel1a,smin,smax,eps)
     End Do
    Endif 
   Endif
   erg = coupC(1) * coupC(3) * int1(4)                               &
     & + 4._dp * ( coupC(1) * coupC(4) * mf(3) * mf(4) * int1(3)      &
     &          + coupC(2) * coupC(3) * mf(1) * mf(2) * int1(2) )    &
     & + 16._dp * coupC(2) * coupC(4) * mf(1) * mf(2) * mf(3) * mf(4) * int1(1)

  Endif

 End Subroutine F3BDscalarSSint


 Real(dp) Function F3BDscalarSSkernel1a(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s

  F3BDscalarSSkernel1a = 0._dp
  If (int_v.Eq.1) Then
   F3BDscalarSSkernel1a =  kappa(s,mf2(1),mf2(2)) * kappa(s,mf2(3),mf2(4))   &
                        &  / (s * ( (s-mS2)**2 + mSG2 ) )
  Else If (int_v.Eq.2) Then
   F3BDscalarSSkernel1a = kappa(s,mf2(1),mf2(2)) * kappa(s,mf2(3),mf2(4))    &
                  &  * (s - mf2(3) - mf2(4)) / (s * ( (s-mS2)**2 + mSG2 ) )
  Else If (int_v.Eq.3) Then
   F3BDscalarSSkernel1a =  kappa(s,mf2(1),mf2(2)) * kappa(s,mf2(3),mf2(4))   &
                  &  * (mf2(1) + mf2(2) - s) / (s * ( (s-mS2)**2 + mSG2 ) )
  Else If (int_v.Eq.4) Then
   F3BDscalarSSkernel1a =  kappa(s,mf2(1),mf2(2)) * kappa(s,mf2(3),mf2(4))   &
                       &  * (mf2(1) + mf2(2) - s)  * (s - mf2(3) - mf2(4))   &
                       &  / (s * ( (s-mS2)**2 + mSG2 ) )
  End If
 End Function F3BDscalarSSkernel1a


 Subroutine  F3BDscalarSSkernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar

  erg = 0._dp

  Do i1=1,2     
   m12 = mf2(1)
   m22 = mf2(2)
   sbar = s(i1)
   sumI = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI * kappa(sbar,m12,m22) / sbar

   sumI = sumI / ( (sbar-mS2)**2 + mSG2 )
   erg(1) = erg(1) + sumI 
   erg(2) = erg(2) + sumI * (sbar - mf2(3) - mf2(4))
   erg(3) = erg(3) + sumI * (mf2(1) + mf2(2) - sbar)
   erg(4) = erg(4) + sumI * (mf2(1) + mf2(2) -sbar)  * (sbar - mf2(3) - mf2(4))
  Enddo

 End Subroutine  F3BDscalarSSkernel1


 Subroutine  F3BDscalarSSkernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar
 
  erg = 0._dp

  Do i1=1,2     
   m12 = mf2(1)
   sbar = s(i1)
   sumI = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI * kappa(sbar,m12,m22) / sbar

   sumI = sumI* (mf2(1) - sbar) / ( (sbar-mS2)**2 + mSG2 )
   erg(1) = erg(1) + sumI 
   erg(2) = erg(2) + sumI * (sbar - mf2(3) - mf2(4)) 
  Enddo

 End Subroutine  F3BDscalarSSkernel2


 Subroutine  F3BDscalarSSkernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSSint for the integration
 ! written by Werner Porod, 17.2.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,m22,sumI,sbar

  erg = 0._dp

  Do i1=1,2     
   m12 = mf2(1)
   m22 = mf2(2)
   sbar = s(i1)
   sumI = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(4)
   sumI = sumI * kappa(sbar,m12,m22) / sbar

   sumI = sumI * (sbar - mf2(3) - mf2(4)) / ( (sbar-mS2)**2 + mSG2 )
   erg(1) = erg(1) + sumI 
   erg(2) = erg(2) + sumI * (mf2(1) + mf2(2) -sbar)
  Enddo

 End Subroutine  F3BDscalarSSkernel3


 Subroutine F3BDscalarSTint(scalar,mf,coup,eps,Integrate,int1,erg)
 !-----------------------------------------------------------------------
 ! gives the integral of the matrix element squared of
 ! M^2_fi = Tr[(C[1] Pl + C[2] Pr).(P1+mf[1]).(C[3] Pl + C[4] Pr).(P2+mf[2]).
 !         (C[5] Pl + C[6] Pr).(P3+mf[3]).(C[7] Pl + C[8] Pr).(P4+mf[4])] /
 !           ( (p1-p4)^2 - m^2_G  + I m^2_G \[CapitalGamma]_G )
 !            * (p1-p2)^2 - m^2_S  - I m^2_S \[CapitalGamma]_S ) )
 ! in terms of the Mandelstam variables {s,t,u}, which are given by
 ! s=(p1-p4)^2, t=(p1-p2)^2, and u=(p1-p3)^2.
 ! Here {C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8]} 
 ! are the couplings, {mf[1],mf[2],mf[3],mf[4]} are the fermion masses,
 ! scalar = {m_S^2, m_S \[CapitalGamma]_S, m_T^2, m_T \[CapitalGamma]_T} 
 ! are the masses squared and the mass times the total decay widths of the
 ! scalar boson in the s- and t-channel, respectively.
 ! written by Werner Porod, 7.2.2000
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: scalar(4),mf(4),eps
  Complex(dp), Intent(in) :: coup(8)
  Complex(dp), Intent(inout) :: int1(8)
  Complex(dp), Intent(out) :: erg
  Logical, Intent(inout) :: Integrate

  Integer :: i1,Imin,Imax
  Real(dp) :: smin,smax,sG(2),int1a(16),int1b(8),int1c(4),int1d(2)
  Complex(dp) :: sumI(8)

  If (Abs(mf(1)).Lt.( Abs(mf(2))+Abs(mf(3))+Abs(mf(4))  ) ) Then
   erg = (0._dp,0._dp)
   Return
  Endif

  sumI(8) = 2._dp * ( coup(1) * coup(3) * coup(5) * coup(7)    &
          &        + coup(2) * coup(4) * coup(6) * coup(8))   &
          &       * mf(1) * mf(2) * mf(3) * mf(4)
  sumI(7) = ( coup(2) * coup(3) * coup(5) * coup(7)           &
          & + coup(1) * coup(4) * coup(6) * coup(8)) * mf(2) * mf(3)
  sumI(6) = ( coup(1) * coup(3) * coup(6) * coup(7)           &
          & + coup(2) * coup(4) * coup(5) * coup(8)) * mf(1) * mf(4)
  sumI(5) = ( coup(1) * coup(4) * coup(5) * coup(7)           &
          & + coup(2) * coup(3) * coup(6) * coup(8)) * mf(3) * mf(4)
  sumI(4) = ( coup(2) * coup(4) * coup(6) * coup(7)           &
          & + coup(1) * coup(3) * coup(5) * coup(8)) * mf(1) * mf(2)
  sumI(3) = ( coup(1) * coup(4) * coup(6) * coup(7)           &
          & + coup(2) * coup(3) * coup(5) * coup(8)) * mf(2) * mf(4)
  sumI(2) = ( coup(2) * coup(4) * coup(5) * coup(7)           &
          & + coup(1) * coup(3) * coup(6) * coup(8)) * mf(1) * mf(3)
  sumI(1) = 0.5_dp * ( coup(2) * coup(3) * coup(6) * coup(7)   &
          &         + coup(1) * coup(4) * coup(5) * coup(8) )

  erg = Sum( Abs(sumI) )

  If (erg.Eq.ZeroC) Then
   If (Integrate) Integrate = .False.
   Return
  Endif

  mf2 = mf**2

  If (Integrate) Then
   mS2 = scalar(1)**2
   mSgS = scalar(1)*scalar(2)
   mSgS2 = mSgS**2
   mT2 = scalar(3)**2
   mTgT = scalar(3)*scalar(4)
   mTgT2 = mTgT**2
   int1 = ZeroC

   smax = (Abs(mf(1))-Abs(mf(4)))**2
   smin = (Abs(mf(3))+Abs(mf(2)))**2
   sG(1) = (scalar(1)-2._dp*scalar(2))**2
   sG(2) = (scalar(1)+2._dp*scalar(2))**2
   Imin = 3
   Imax = 0
   Do i1=1,2
    If (smin.Lt.sG(3-i1)) Imin = 3 - i1
    If (smax.Gt.sG(i1)) Imax = i1
   Enddo

   If ((mf(2).Ne.0._dp).And.(mf(3).Ne.0._dp).And.(mf(4).Ne.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel1,16,smin,smax,int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDscalarSTkernel1,16,sG(Imin),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)

    Else
     Call DgaussInt(F3BDscalarSTkernel1,16,smin,sG(Imin),int1a,eps)
     int1 = int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDscalarSTkernel1,16,sG(Imin),sG(Imax),int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
     Call DgaussInt(F3BDscalarSTkernel1,16,sG(Imax),smax,int1a,eps)
     int1 = int1 + int1a(1:15:2) + Ic * int1a(2:16:2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel8,2,smin,smax,int1d,eps)
     int1(1) = int1d(1) + Ic * int1d(2)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(1) = int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDscalarSTkernel8,2,sG(Imin),smax,int1d,eps)
     int1(1) = int1(1) + int1d(1) + Ic * int1d(2)

    Else
     Call DgaussInt(F3BDscalarSTkernel8,2,smin,sG(Imin),int1d,eps)
     int1(1) = int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDscalarSTkernel8,2,sG(Imin),sG(Imax),int1d,eps)
     int1(1) = int1(1) + int1d(1) + Ic * int1d(2)
     Call DgaussInt(F3BDscalarSTkernel8,2,sG(Imax),smax,int1d,eps)
     int1(1) = int1(1) + int1d(1) + Ic * int1d(2)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(3).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel7,4,smin,smax,int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel7,4,sG(Imin),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDscalarSTkernel7,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(6) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel7,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel7,4,sG(Imax),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(6) = int1(6) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(2).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel6,4,smin,smax,int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(2) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(2) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel6,4,sG(Imin),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(2) = int1(2) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDscalarSTkernel6,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(2) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel6,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(2) = int1(2) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel6,4,sG(Imax),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(2) = int1(2) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif ((mf(3).Eq.0._dp).And.(mf(4).Eq.0._dp)) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel5,4,smin,smax,int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel5,4,sG(Imin),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)

    Else
     Call DgaussInt(F3BDscalarSTkernel5,4,smin,sG(Imin),int1c,eps)
     int1(1) = int1c(1) + Ic * int1c(2)
     int1(4) = int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel5,4,sG(Imin),sG(Imax),int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)
     Call DgaussInt(F3BDscalarSTkernel5,4,sG(Imax),smax,int1c,eps)
     int1(1) = int1(1) + int1c(1) + Ic * int1c(2)
     int1(4) = int1(4) + int1c(3) + Ic * int1c(4)
    Endif

   Elseif (mf(2).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel4,8,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel4,8,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDscalarSTkernel4,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(5) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel4,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel4,8,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(5) = int1(5) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(3).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel3,8,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel3,8,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDscalarSTkernel3,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(3) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(6) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel3,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel3,8,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(3) = int1(3) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(6) = int1(6) + int1b(7) + Ic * int1b(8)
    Endif

   Elseif (mf(4).Eq.0._dp) Then
    If ((Imin.Eq.3).Or.(Imax.Eq.0).Or.((Imin-1).Eq.Imax)) Then
     Call DgaussInt(F3BDscalarSTkernel2,8,smin,smax,int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)

    Elseif (Imin.Eq.Imax) Then
     Call DgaussInt(F3BDscalarSTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel2,8,sG(Imin),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)

    Else
     Call DgaussInt(F3BDscalarSTkernel2,8,smin,sG(Imin),int1b,eps)
     int1(1) = int1b(1) + Ic * int1b(2)
     int1(2) = int1b(3) + Ic * int1b(4)
     int1(4) = int1b(5) + Ic * int1b(6)
     int1(7) = int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel2,8,sG(Imin),sG(Imax),int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)
     Call DgaussInt(F3BDscalarSTkernel2,8,sG(Imax),smax,int1b,eps)
     int1(1) = int1(1) + int1b(1) + Ic * int1b(2)
     int1(2) = int1(2) + int1b(3) + Ic * int1b(4)
     int1(4) = int1(4) + int1b(5) + Ic * int1b(6)
     int1(7) = int1(7) + int1b(7) + Ic * int1b(8)
    Endif

   Endif
  Endif

  erg = Sum(sumI * int1)

 End Subroutine F3BDscalarSTint


 Subroutine F3BDscalarSTkernel1(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(16)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar, &
       &           factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + (mf2(1) - mf2(4))*(mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = DiffLog
   sum2 = DiffTan

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(16) = erg(16) + factI
   erg(15) = erg(15) + factR

   erg(14) = erg(14) + factI * (-sbar + mf2(1) + mf2(4))
   erg(13) = erg(13) + factR * (-sbar + mf2(1) + mf2(4))

   erg(12) = erg(12) + factI * (sbar - mf2(2) - mf2(3)) 
   erg(11) = erg(11) + factR * (sbar - mf2(2) - mf2(3)) 

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) - mf2(2) ) * DiffLog
   sum2 = - ( mT2 - mf2(1) - mf2(2) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(10) = erg(10) + factI
   erg(9) = erg(9) + factR

   sum1 = diff - mTgT * DiffTan               &
        & + ( mT2 - mf2(3) - mf2(4) ) * DiffLog
   sum2 = ( mT2 - mf2(3) - mf2(4) ) * DiffTan +  mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2
   factI = ReProp * sum2 - ImProp * sum1

   erg(8) = erg(8) + factI
   erg(7) = erg(7) + factR

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan  + mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(6) = erg(6) + factI
   erg(5) = erg(5) + factR

   sum1 = - diff + mTgT * DiffTan                    &
        & - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   sum1 = 2._dp * ( sbar * ( diff - mTgT * DiffTan )                  &
        &        + (mT2 * sbar - mf2(1) * mf2(3) - mf2(2) * mf2(4) ) * DiffLog)

   sum2 = 2._dp * (( mT2 * sbar - mf2(1) * mf2(3) - mf2(2) * mf2(4) ) *DiffTan &
        &        + mTgT * sbar * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel1


 Subroutine F3BDscalarSTkernel2(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar, &
       &           factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   m22 = mf2(2)
   sumI = sumI + mf2(1) * (mf2(3) - mf2(2)) / sbar
   diff = diff * kappa(sbar,m12,m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = DiffLog
   sum2 = DiffTan

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(8) = erg(8) + factI * (-sbar + mf2(1) )
   erg(7) = erg(7) + factR * (-sbar + mf2(1) )

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(3) ) * DiffLog
   sum2 = ( mT2 - mf2(3) ) * DiffTan +  mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2
   factI = ReProp * sum2 - ImProp * sum1

   erg(6) = erg(6) + factI
   erg(5) = erg(5) + factR

   sum1 = - diff + mTgT * DiffTan                    &
        & - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog

   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   sum1 = 2._dp * ( sbar * ( diff - mTgT * DiffTan )          &
        &        + (mT2 * sbar - mf2(1) * mf2(3) ) * DiffLog )

   sum2 = 2._dp * ( ( mT2 * sbar - mf2(1) * mf2(3) ) * DiffTan &
        &        + mTgT * sbar * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel2


 Subroutine F3BDscalarSTkernel3(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar, &
       &           factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m22 = mf2(2)
   sumI = sumI - (mf2(1) - mf2(4))* mf2(2) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = DiffLog
   sum2 = DiffTan

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(8) = erg(8) + factI * (sbar - mf2(2) ) 
   erg(7) = erg(7) + factR * (sbar - mf2(2) ) 

   sum1 = diff - mTgT * DiffTan + ( mT2 - mf2(4) ) * DiffLog

   sum2 = ( mT2 - mf2(4) ) * DiffTan + mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2
   factI = ReProp * sum2 - ImProp * sum1

   erg(6) = erg(6) + factI
   erg(5) = erg(5) + factR

   sum1 = diff - mTgT * DiffTan + ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffLog
   sum2 = ( mT2 + sbar - mf2(2) - mf2(4) ) * DiffTan + mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   sum1 = 2._dp * ( sbar * ( diff - mTgT * DiffTan )          &
        &        + (mT2 * sbar - mf2(2) * mf2(4) ) * DiffLog )

   sum2 = 2._dp * ( ( mT2 * sbar - mf2(2) * mf2(4) ) * DiffTan &
        &        + mTgT * sbar * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel3


 Subroutine F3BDscalarSTkernel4(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar, &
       &           factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   m12 = mf2(3)
   sumI = sumI + (mf2(1) - mf2(4)) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = DiffLog
   sum2 = DiffTan

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(8) = erg(8) + factI * (sbar - mf2(3)) 
   erg(7) = erg(7) + factR * (sbar - mf2(3)) 

   sum1 = - diff + mTgT * DiffTan - ( mT2 - mf2(1) ) * DiffLog

   sum2 = - ( mT2 - mf2(1) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(6) = erg(6) + factI
   erg(5) = erg(5) + factR

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   sum1 = 2._dp * ( sbar * ( diff - mTgT * DiffTan )          &
        &        + (mT2 * sbar - mf2(1) * mf2(3) ) * DiffLog )

   sum2 = 2._dp * (( mT2 * sbar - mf2(1) * mf2(3) ) * DiffTan &
        &        + mTgT * sbar * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel4


 Subroutine F3BDscalarSTkernel5(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar,  &
       &     factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(2) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m22 = mf2(2)
   sumI = sumI - mf2(1) * mf2(2) / sbar
   diff = diff * Abs(sbar-m22) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = diff - mTgT * DiffTan +  mT2 * DiffLog

   sum2 = mT2 * DiffTan +  mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2
   factI = ReProp * sum2 - ImProp * sum1

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   erg(2) = erg(2) + 2._dp * sbar * factI
   erg(1) = erg(1) + 2._dp * sbar * factR

  Enddo

 End Subroutine F3BDscalarSTkernel5


 Subroutine F3BDscalarSTkernel6(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
  Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,sum1,sum2,tmin,tmax,sumI,diff,sbar,   &
      &    factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(3) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   m12 = mf2(3)
   sumI = sumI + mf2(1) * mf2(3) / sbar
   diff = diff * Abs(sbar-m12) / sbar

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = - diff + mTgT * DiffTan - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffLog
   sum2 = - ( mT2 + sbar - mf2(1) - mf2(3) ) * DiffTan - mTgT * DiffLog

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI
   erg(3) = erg(3) + factR

   sum1 = 2._dp * ( sbar * ( diff - mTgT * DiffTan )          &
        &        + (mT2 * sbar - mf2(1) * mf2(3) ) * DiffLog )

   sum2 = 2._dp * (( mT2 * sbar - mf2(1) * mf2(3) ) * DiffTan &
        &        + mTgT * sbar * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel6


 Subroutine F3BDscalarSTkernel7(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1
  Real(dp) :: m12,m22,sum1,sum2,tmin,tmax,sumI,diff,sbar,  &
       &   factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) + mf2(4) - sbar
   m12 = mf2(1)
   m22 = mf2(4)
   diff = kappa(sbar,m12,m22)

   tmin = 0.5_dp * ( sumI - diff)
   tmax = 0.5_dp * ( sumI + diff)

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = Log( (tmax - mT2) / (tmin - mT2) )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (tmin-mT2) / mTgT )
    DiffLog = 0.5_dp  * ( Log( ((tmax - mT2)**2 + mTgT2 ) &
            &                / ((tmin - mT2)**2 + mTgT2 ) ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = DiffLog
   sum2 = DiffTan

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(4) = erg(4) + factI * sbar
   erg(3) = erg(3) + factR * sbar

   sum1 = 2._dp * sbar * ( diff - mTgT * DiffTan + mT2 * DiffLog )

   sum2 = 2._dp * sbar * ( mT2 * DiffTan + mTgT * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel7


 Subroutine F3BDscalarSTkernel8(s,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDscalarSTint for the integration
 ! written by Werner Porod, 21.6.00
 ! 28.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: m12,sum1,sum2,tmax,sumI,diff,sbar, &
         &  factR,factI,DiffTan,DiffLog,ReProp,ImProp

  erg = 0

  Do i1=1,2
   sbar = s(i1)
   sumI = mf2(1) - sbar
   m12 = mf2(1)
   diff = Abs(sbar-m12)

   tmax = sumI

   If (mTgT.Eq.0._dp) Then
    DiffTan = 0._dp
    DiffLog = 0.5_dp * Log( (tmax - mT2)**2 / mT2**2 )
   Else
    DiffTan = Atan( (tmax-mT2) / mTgT ) - Atan( (-mT2) / mTgT )
    DiffLog = 0.5_dp * Log( ((tmax - mT2)**2 + mTgT2 ) / (mT2**2 + mTgT2 ) )
   Endif

   ReProp = (sbar-mS2) / ( (sbar-mS2)**2 + mSgS2 )
   ImProp = mSgS / ( (sbar-mS2)**2 + mSgS2 )

   sum1 = 2._dp * sbar * ( diff - mTgT * DiffTan  + mT2 * DiffLog )
   sum2 = 2._dp * sbar * ( mT2 * DiffTan + mTgT * DiffLog )

   factR = ReProp * sum1 + ImProp * sum2 
   factI = ReProp * sum2 - ImProp * sum1 

   erg(2) = erg(2) + factI
   erg(1) = erg(1) + factR

  Enddo

 End Subroutine F3BDscalarSTkernel8

 Subroutine IntegrateGaugeSscalarS(Scalar,Mass,coup2,deltaM,epsI &
                                & ,IntegralsC4,n_C,resC, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 29.10.2000: porting to f90, the second dimension of IntegralsC4 has to
 !             be 12 in the calling program
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: n_C
  Real(dp), Intent(in) :: Scalar(4),deltaM,epsI
  Real(dp), Intent(inout) :: mass(4)
  Complex(dp), Intent(in) :: coup2(8)
  Complex(dp), Intent(inout) :: IntegralsC4(:,:)
  Complex(dp), Intent(out) :: resC
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM
  Complex(dp) :: intC4(4)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateGaugeSscalarS'

  resC = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If (Scalar(3).eq.0._dp) then
     Iname = Iname - 1
     Return
    End If    
    If ( (Abs(mass(1)).Gt.(Abs(mass(2))+Scalar(1))) .And. &
       & (Scalar(1).Gt.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
    If ( (Abs(mass(1)).Gt.(Abs(mass(2))+Scalar(3))) .And. &
       & (Scalar(3).Gt.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_C
   If ( (Abs(IntegralsC4(i1,2)-mass(2)).Lt.deltaM).And.   &
      & (Abs(IntegralsC4(i1,3)-mass(3)).Lt.deltaM).And.   &
      & (Abs(IntegralsC4(i1,4)-mass(4)).Lt.deltaM).And.   &
      & (IntegralsC4(i1,5).Eq.scalar(1)).And.             &
      & (IntegralsC4(i1,6).Eq.scalar(2)).And.             &
      & (IntegralsC4(i1,7).Eq.scalar(3)).And.             &
      & (IntegralsC4(i1,8).Eq.scalar(4))                  ) Then

    intC4 = IntegralsC4(i1,9:12)
    Integrate = .False.
    Call F3BDGaugeSscalarSint(Scalar,mass,coup2,epsI,Integrate,intC4,resC)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDGaugeSscalarSint(Scalar,mass,coup2,epsI,Integrate,intC4,resC)
  If (Integrate) Then
   n_C = N_C + 1
   IntegralsC4(n_C,1:4) = mass
   IntegralsC4(n_C,5:8) = scalar
   IntegralsC4(n_C,9:12) = intC4
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateGaugeSscalarS

 Subroutine IntegrateGaugeSscalarT(gauge,Mass,coup2,deltaM,epsI &
                                & ,IntegralsC8,n_C,resC, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 29.10.2000: porting to f90, the second dimension of IntegralsC8 has to
 !             be 16 in the calling program
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: n_C
  Real(dp), Intent(in) :: gauge(4),deltaM,epsI
  Complex(dp), Intent(in) :: coup2(8)
  Real(dp), Intent(inout) :: mass(4)
  Complex(dp), Intent(inout) :: IntegralsC8(:,:)
  Complex(dp), Intent(out) :: resC
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM
  Complex(dp) :: intC8(8)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateGaugeSsclarT'

  resC = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If (Gauge(3).eq.0._dp) then
     Iname = Iname - 1
     Return
    End If    
    If ( (Abs(mass(1)).Gt.(Abs(mass(4))+Gauge(1))) .And. &
       & (Gauge(1).Gt.(Abs(mass(3))+Abs(mass(2)))) ) Then
     Iname = Iname - 1
     Return
    End If
    If ( (Abs(mass(1)).Gt.(Abs(mass(2))+Gauge(3))) .And. &
       & (Gauge(3).Gt.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_C
   If ( (Abs(IntegralsC8(i1,2)-mass(2)).Lt.deltaM).And.  &
      & (Abs(IntegralsC8(i1,3)-mass(3)).Lt.deltaM).And.  &
      & (Abs(IntegralsC8(i1,4)-mass(4)).Lt.deltaM).And.  &
      & (Abs(IntegralsC8(i1,5)-gauge(1)).Lt.deltaM).And. &
      & (Abs(IntegralsC8(i1,6)-gauge(2)).Lt.deltaM).And. &
      & (Abs(IntegralsC8(i1,7)-gauge(3)).Lt.deltaM).And. &
      & (Abs(IntegralsC8(i1,8)-gauge(4)).Lt.deltaM) ) Then

    intC8 = IntegralsC8(i1,9:16)
    Integrate = .False.
    Call F3BDgaugeSscalarTint(gauge,mass,coup2,epsI,Integrate,intC8,resC)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDgaugeSscalarTint(gauge,mass,coup2,epsI,Integrate,intC8,resC)
  If (Integrate) Then
   n_C = N_C + 1
   IntegralsC8(n_C,1:4) = mass
   IntegralsC8(n_C,5:8) = gauge
   IntegralsC8(n_C,9:16) = intC8
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateGaugeSscalarT

 Subroutine IntegrateGaugeSS(Gauge,Mass,coup1,deltaM,epsI &
                            & ,IntegralsR4,n_IR4,resR, check)
 !-----------------------------------------------------------------------
 ! input:  gauge(1) ... mass of gauge boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 29.10.2000: porting to f90, the second dimension of  IntegralsR4 has to
 !             8 in the calling subroutine
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: n_IR4
  Real(dp), Intent(in) :: Gauge(2),deltaM,epsI
  Real(dp), Intent(inout) :: IntegralsR4(:,:),mass(4)
  Real(dp), Intent(out) :: resR
  Complex(dp), Intent(in) :: coup1(4)
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM,intR4(4)
  Logical Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateGaugeSS'

  resR = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Gauge(1))) .And. &
       & (Gauge(1).Ge.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_IR4
   If ( (Abs(IntegralsR4(i1,2)-mass(2)).Lt.deltaM).And. &
      & (Abs(IntegralsR4(i1,3)-mass(3)).Lt.deltaM).And. &
      & (Abs(IntegralsR4(i1,4)-mass(4)).Lt.deltaM) ) Then

    intR4 = IntegralsR4(i1,5:8)
    Integrate = .False.
    Call F3BDgaugeSSint(Gauge,mass,coup1,epsI,Integrate,intR4,resR)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDgaugeSSint(Gauge,mass,coup1,epsI,Integrate,intR4,resR)
  If (Integrate) Then
   n_IR4 = N_IR4 + 1
   IntegralsR4(n_IR4,1:4) = mass
   IntegralsR4(n_IR4,5:8) = intR4
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateGaugeSS

 Subroutine IntegrategaugeST(gauge,Mass,coup2,deltaM,epsI,IntegralsC8,n_C &
                            &,resC, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 29.10.2000: porting to f90, the second dimension of IntegralsC8 has to 
 !             be 12
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: n_C
  Real(dp), Intent(in) :: gauge(4),deltaM,epsI
  Real(dp), Intent(inout) :: mass(4)
  Complex(dp), Intent(in) :: coup2(8)
  Complex(dp), Intent(inout) :: IntegralsC8(:,:)
  Complex(dp), Intent(out) :: resC
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM
  Complex(dp) :: intC8(8)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateGaugeST'

  resC = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If ( (Abs(mass(1)).Ge.(Abs(mass(4))+Gauge(1))) .And. &
       & (Gauge(1).Ge.(Abs(mass(3))+Abs(mass(2)))) )  Then
     Iname = Iname - 1
     Return
    End If
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Gauge(3))) .And. &
       & (Gauge(3).Ge.(Abs(mass(3))+Abs(mass(4)))) )  Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_C
   If ( (Abs(IntegralsC8(i1,2)-mass(2)).Lt.deltaM).And.  &
      & (Abs(IntegralsC8(i1,3)-mass(3)).Lt.deltaM).And.  &
      & (Abs(IntegralsC8(i1,4)-mass(4)).Lt.deltaM) ) Then

    intC8 = IntegralsC8(i1,5:12)
    Integrate = .False.
    Call F3BDgaugeSTint(gauge,mass,coup2,epsI,Integrate,intC8,resC)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDgaugeSTint(gauge,mass,coup2,epsI,Integrate,intC8,resC)
  If (Integrate) Then
   n_C = N_C + 1
   IntegralsC8(n_C,1:4) = mass
   IntegralsC8(n_C,5:12) = intC8
  Endif

  Iname = Iname - 1

 End Subroutine IntegrategaugeST

 Subroutine IntegrateScalarS1S2(Scalar,Mass,coup2,deltaM,epsI,IntegralsC4,n_C &
                             & ,resC, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 27.10.2000: porting to f90, second dimension of IntegralsC4 has to be 12
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: n_C
  Real(dp), Intent(in) :: Scalar(4),deltaM,epsI
  Real(dp), Intent(inout) :: mass(4)
  Complex(dp), Intent(in) :: coup2(8)
  Complex(dp), Intent(inout) :: IntegralsC4(:,:)
  Complex(dp), Intent(out) :: resC
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM
  Complex(dp) :: intC4(4)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateScalarS1S2'

  resC = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If ((Scalar(1).Eq.0._dp).Or.(Scalar(3).Eq.0._dp)) Then
     Iname = Iname - 1
     Return
    End If    
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Scalar(1))) .And. &
       & (Scalar(1).Ge.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Scalar(3))) .And. &
       & (Scalar(3).Ge.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_C
   If ( (Abs(IntegralsC4(i1,2)-mass(2)).Lt.deltaM).And.     &
      & (Abs(IntegralsC4(i1,3)-mass(3)).Lt.deltaM).And.     &
      & (Abs(IntegralsC4(i1,4)-mass(4)).Lt.deltaM).And.     &
      & (IntegralsC4(i1,5).Eq.scalar(1)).And.               &
      & (IntegralsC4(i1,6).Eq.scalar(2)).And.               &
      & (IntegralsC4(i1,7).Eq.scalar(3)).And.               &
      & (IntegralsC4(i1,8).Eq.scalar(4))                    ) Then

    intC4 = IntegralsC4(i1,9:12)
    Integrate = .False.
    Call F3BDscalarS1S2int(Scalar,mass,coup2,epsI,Integrate,intC4,resC)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDscalarS1S2int(Scalar,mass,coup2,epsI,Integrate,intC4,resC)
  If (Integrate) Then
   n_C = N_C + 1
   IntegralsC4(n_C,1:4) = mass
   IntegralsC4(n_C,5:8) = scalar
   IntegralsC4(n_C,9:12) = intC4
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateScalarS1S2

 Subroutine IntegrateScalarSS(Scalar,Mass,coup1,deltaM,epsI,IntegralsR4,n_IR4 &
                             & ,resR, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 !         Scalar(2) ... decy width of the scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 27.10.2000: portation to f90, second dimension of  IntegralsR4 has
 !             to be 10 
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: n_IR4
  Real(dp), Intent(in) :: Scalar(2),deltaM,epsI
  Real(dp), Intent(inout) :: IntegralsR4(:,:),mass(4)
  Real(dp), Intent(out) :: resR
  Complex(dp), Intent(in) :: coup1(4)
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM,intR4(4)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateScalarSS'

  resR = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Scalar(1))) .And. &
       & (Scalar(1).Ge.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
    If ( Scalar(1).Eq.0._dp ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_IR4
   If ( (Abs(IntegralsR4(i1,2)-mass(2)).Lt.deltaM).And.   &
    &   (Abs(IntegralsR4(i1,3)-mass(3)).Lt.deltaM).And.   &
    &   (Abs(IntegralsR4(i1,4)-mass(4)).Lt.deltaM).And.   &
    &   (IntegralsR4(i1,5).Eq.scalar(1)).And.             &
    &   (IntegralsR4(i1,6).Eq.scalar(2))                  ) Then

    intR4 = IntegralsR4(i1,7:10)
    Integrate = .False.
    Call F3BDScalarSSint(Scalar,mass,coup1,epsI,Integrate,intR4,resR)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDScalarSSint(Scalar,mass,coup1,epsI,Integrate,intR4,resR)
  If (Integrate) Then
   n_IR4 = N_IR4 + 1
   IntegralsR4(n_IR4,1:4) = mass
   IntegralsR4(n_IR4,5:6) = scalar
   IntegralsR4(n_IR4,7:10) = intR4
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateScalarSS

 Subroutine IntegrateScalarST(scalar,Mass,coup2,deltaM,epsI,IntegralsC8,n_C &
                           & ,resC, check)
 !-----------------------------------------------------------------------
 ! input:  Scalar(1) ... mass of Scalar boson
 ! output:
 ! written by Werner Porod, 10.6.2000
 ! 28.10.2000: porting to f90, the second dimension of IntegralsC8 has to
 !             be 16 in the calling program
 ! 12.09.03: new optional variable check, if .TRUE. then a check will be
 !           done if the intermediate state is real or virtual. In the
 !           former case resR is set to 0 
 !-----------------------------------------------------------------------
  Implicit None

  Integer, Intent(inout) :: n_C
  Real(dp), Intent(in) :: scalar(4),deltaM,epsI
  Real(dp), Intent(inout) :: mass(4)
  Complex(dp), Intent(in) :: coup2(8)
  Complex(dp), Intent(inout) :: IntegralsC8(:,:)
  Complex(dp), Intent(out) :: resC
  Logical, Optional :: check

  Integer :: i1
  Real(dp) :: diffM
  Complex(dp) :: intC8(8)
  Logical :: Integrate

  Iname = Iname + 1
  NameOfUnit(Iname) = 'IntegrateScalarST'

  resC = 0._dp

  DiffM = Abs(mass(1)) - Abs(mass(2)) - Abs(mass(3)) - Abs(mass(4))
  If (DiffM.Le.0._dp) Then
   Iname = Iname - 1
   Return
  Endif

  If ( (Abs(mass(2))/DiffM).Lt.deltaM ) mass(2) = 0._dp
  If ( (Abs(mass(3))/DiffM).Lt.deltaM ) mass(3) = 0._dp
  If ( (Abs(mass(4))/DiffM).Lt.deltaM ) mass(4) = 0._dp

  If (Present(check)) Then
   If (check) Then
    If ((Scalar(1).Eq.0._dp).Or.(Scalar(3).Eq.0._dp)) Then
     Iname = Iname - 1
     Return
    End If    
    If ( (Abs(mass(1)).Ge.(Abs(mass(4))+Scalar(1))) .And. &
       & (Scalar(1).Ge.(Abs(mass(3))+Abs(mass(2)))) ) Then
     Iname = Iname - 1
     Return
    End If
    If ( (Abs(mass(1)).Ge.(Abs(mass(2))+Scalar(3))) .And. &
       & (Scalar(3).Ge.(Abs(mass(3))+Abs(mass(4)))) ) Then
     Iname = Iname - 1
     Return
    End If
   End If
  End If

  Do i1=1,n_C
   If ( (Abs(IntegralsC8(i1,2)-mass(2)).Lt.deltaM).And.    &
      & (Abs(IntegralsC8(i1,3)-mass(3)).Lt.deltaM).And.    &
      & (Abs(IntegralsC8(i1,4)-mass(4)).Lt.deltaM).And.    &
      & (IntegralsC8(i1,5).Eq.scalar(1)).And.              &
      & (IntegralsC8(i1,6).Eq.scalar(2)).And.              &
      & (IntegralsC8(i1,7).Eq.scalar(3)).And.              &
      & (IntegralsC8(i1,8).Eq.scalar(4))                   ) Then

    intC8 = IntegralsC8(i1,9:16)
    Integrate = .False.
    Call F3BDscalarSTint(scalar,mass,coup2,epsI,Integrate,intC8,resC)
    Iname = Iname - 1
    Return
   Endif
  Enddo

  Integrate = .True.
  Call F3BDscalarSTint(scalar,mass,coup2,epsI,Integrate,intC8,resC)
  If (Integrate) Then
   n_C = N_C + 1
   IntegralsC8(n_C,1:4) = mass
   IntegralsC8(n_C,5:8) = scalar
   IntegralsC8(n_C,9:16) = intC8
  Endif

  Iname = Iname - 1

 End Subroutine IntegrateScalarST


 Subroutine IntGaugeSS1(gauge,mf,erg)
 !-----------------------------------------------------------------------
 ! auxiliary function for subroutine F3BDgaugeSS for the integration
 ! written by Werner Porod, 8.1.00
 ! 29.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf,gauge(2)
  Real(dp), Intent(out) :: erg

  Real(dp) :: mS,gS,mS2,mf2a,mgS,mgS2,ratio,mf4
  Complex(dp) :: ergC

  mS = gauge(1)
  gS = gauge(2)
  mf2a = mf**2
  mS2 = mS**2
  If (gS.Eq.0._dp) Then
   ratio = mf2a / mS2
   If (ratio.Gt.0.01_dp) Then
    mf4 = mf2a**2
    erg = ratio * ( mS2 * ( 2._dp * mS2 - mf2a) - mf4 / 3._dp) &
      & + 2._dp * mS2 * (mS2-mf2a) * Log(1._dp - ratio)
   Else
    ratio = 1._dp / mS2**2 
    erg = mf2a**4 * ratio *(1._dp/ 6._dp + 0.1_dp*mf2a*mS2*ratio    &
        &                 + mf2a**2 * ratio / 15._dp             &
        &                 + mf2a**3 * ratio**2 * mS2 / 21._dp )
   Endif

  Else
   mgS = mS * gS
   mgS2 = mgS**2
   If ( Abs(mf/mS).Gt.0.15_dp) Then
    ergC = Ic * ( Log( mgS / (mgS - Ic * (mf2a-mS2)) )  &
         &      - Log(mgS / (mgS + Ic * (mf2a-mS2)))   &
         &      + Log(mgS / (mgS - Ic * mS2 ) )       &
         &      - Log(mgS / (mgS + Ic * mS2 ) ) )
    erg = ( 3._dp * mS2 * (mS2 - mf2a) - mgS2 )                        &
      &   * Log( ((mS2 - mf2a)**2 + mgS2) / (mS2**2+mgS2) ) / 3._dp   &
      & + 2._dp * mf2a * ( 2._dp * mS2 - mf2a) / 3._dp                    &
      & + ( 3._dp * mgS2 * (2._dp * mS2 - mf2a)                         &
      & - (2._dp * mS2 + mf2a) * (mS2 - mf2a)**2 ) * Real(ergC,dp)     &
      &    / (6._dp*mgS)
   Else
    ratio = 1._dp / (mS2**2 + mgS2)
    erg = mf2a**4 * ratio *(1._dp/ 6._dp + 0.1_dp*mf2a*mS2*ratio             &
        &       + mf2a**2 * ratio *(3._dp - 4._dp *mgS2 *ratio)/45._dp       &
        &       + mf2a**3 * ratio**2 * mS2 * (1._dp-2._dp*ratio*mgS2) / 21._dp)
   Endif
  Endif
  
 End Subroutine IntGaugeSS1


 Subroutine IntScalarSS1(scalar,mf,erg)
 !-----------------------------------------------------------------------
 ! written by Werner Porod, 3.1.00
 ! 27.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf,scalar(2)
  Real(dp), Intent(out) :: erg

  Real(dp) :: mS,gS,mS2,mgS,mgS2,ratio,mf2a
  Complex(dp) :: ergC

  mS = scalar(1)
  gS = scalar(2)
  mf2a = mf**2
  mS2 = mS**2
  If (gS.Eq.0._dp) Then
   If (Abs(mf/mS).Gt.0.2_dp) Then
    erg = mf2a * ( 3._dp * mS2 - 2.5_dp * mf2a)               &
      & + (3._dp*mS2-mf2a)*(mS2-mf2a) * Log(1._dp - mf2a/mS2)
   Else
    ratio = 1._dp / mS2**2
    erg = mf2a**4 * ratio * ( 1._dp/12._dp + mS2*ratio*mf2a/15._dp       &
        &                  + 3._dp*mS2**2 *(mf2a*ratio)**2 /6.d1      &
        &                  + 4._dp * mS2**3 * (mf2a*ratio)**3 /1.05d2)
   Endif

  Else
   mgS = mS * gS
   mgS2 = mgS**2
   If (Abs(mf/mS).Gt.0.2_dp) Then
    ergC = Ic * ( Log(mgS / (mgS - Ic * (mf2a-mS2)) )  &
         &    - Log(mgS / (mgS + Ic * (mf2a-mS2)))    &
         &    + Log(mgS / (mgS - Ic * mS2 ) )        &
         &    - Log(mgS / (mgS + Ic * mS2 ) ) )
    erg = ( mS2 * (1.5_dp * mS2 - 2._dp*mf2a) + 0.5_dp*(mf2a**2-mgS2) )    &
      &    * Log( ((mS2 - mf2a)**2 + mgS2) / (mS2**2+mgS2) )          &
      & + mf2a * ( 2._dp * mS2 - 1.5_dp * mf2a)                           &
      & + ( mgS2 * (3._dp * mS2 - 2._dp * mf2a) - mS2 * (mS2 - mf2a)**2 ) &
      &   * Real(ergC,dp)  / (2._dp*mgS)
   Else
    ratio = 1._dp / (mgS2+mS2**2)
    erg = mf2a**4 * ratio * ( 1._dp/12._dp + mS2*ratio*mf2a/15._dp         &
        &                + (3._dp*mS2**2-mgS2) *(mf2a*ratio)**2 /6.d1   &
        &                + 4._dp * mS2 * (mS2**2-mgS2)                 &
        &                       * (mf2a*ratio)**3 /1.05d2 )
   Endif

  Endif
  
 End Subroutine IntScalarSS1

End Module ThreeBodyPhaseSpace

