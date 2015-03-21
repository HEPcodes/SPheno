Module ThreeBodyPhaseSpaceS
! comments
!-----------------------------------------------------------
! contains the phase space function needed for 3-body decays
! involving fermions and bosons in intial and/or final state 
! 16.08.2012: adding contributions by Lukas Mitzka
!-----------------------------------------------------------

! load modules
Use Control
Use Mathematics, Only: Kappa, DGaussInt
! load modules

Real(dp), Private ::  rmt2, rjt2, rat2, rbt2, rkt2, Rsquared

Contains

 
!!!!! 
!Integration Routine for the Slepton charge chaning channel
 Subroutine ChiChiInterfIntegrands(s, erg)
 Implicit None

  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(8)

  Integer :: i1

  Real(dp) :: x, x2, yplus, yminus, fOnexrarb, fTwoxrarb, l1abbr, kabbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2

   yplus = 1._dp / ( 2._dp * ( 1 - x + rkt2 ) )                               &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   
   l1abbr =      Log( Abs(  ( yplus  + rbt2 - 1._dp - rmt2 )         &
            &             / ( yminus + rbt2 - 1._dp - rmt2 )   ) )   

   kabbr = (   yplus - yminus  ) / ( 1._dp - x + rkt2 -rat2 )

   fOnexrarb = - l1abbr / ( 1._dp - x + rkt2 - rat2  )

   fTwoxrarb = ( 1._dp + rmt2 -rbt2 )*fOnexrarb - kabbr

   erg(1) = erg(1) + (-rjt2 + rkt2 + rmt2 - 1._dp)*fOnexrarb

   erg(2) = erg(2) + (-rkt2 + x - 1._dp )*fTwoxrarb                           &
                &  + (-rjt2 + rkt2*(2._dp*rmt2 + 1._dp)                       &
                &  -  (rmt2 + 1._dp)*(x - 1._dp)    )*fOnexrarb
   
   erg(3) = erg(3) + (-rjt2 + rkt2 + rmt2 + 1._dp - x)*fOnexrarb
   
   erg(4) = erg(4) + (2._dp*rkt2 - x )*fOnexrarb
   
   erg(5) = erg(5) + fTwoxrarb + (rjt2 - rkt2 - rmt2 - 1._dp)*fOnexrarb 
   
   erg(6) = erg(6) + fTwoxrarb - 2._dp*rmt2*fOnexrarb
   
   erg(7) = erg(7) + ( rjt2 - rkt2 - rmt2 - 1._dp + x)*fOnexrarb + fTwoxrarb
   
   erg(8) = erg(8) + fOnexrarb
   
  End Do  

 End Subroutine ChiChiInterfIntegrands
  

 Subroutine IntegrateChiChiInterference(mass,m_in,r_out,coup,smin,smax,eps,gam)
 Implicit None  

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_ChiaChibInterf(8)
  
  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)
  
  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2  
   
  Call DGaussInt(ChiChiInterfIntegrands,8,smin,smax,I_ChiaChibInterf(:),eps)
  I_ChiaChibInterf(1) = Sqrt(rkt2*rmt2)*I_ChiaChibInterf(1)
  I_ChiaChibInterf(3) = Sqrt(rkt2*rbt2)*I_ChiaChibInterf(3)
  I_ChiaChibInterf(4) = Sqrt(rmt2*rbt2)*I_ChiaChibInterf(4)
  I_ChiaChibInterf(5) = Sqrt(rmt2*rat2)*I_ChiaChibInterf(5)
  I_ChiaChibInterf(6) = Sqrt(rkt2*rat2)*I_ChiaChibInterf(6)
  I_ChiaChibInterf(7) = Sqrt(rat2*rbt2)*I_ChiaChibInterf(7)
  I_ChiaChibInterf(8) = - 2._dp*Sqrt(rat2*rbt2*rkt2*rmt2)*I_ChiaChibInterf(8)
  
  
  gam = ( coup(1) * coup(4) * Conjg(coup(7) * coup(6))                         &
    &   + coup(2) * coup(3) * Conjg(coup(8) * coup(5)) ) * I_ChiaChibInterf(1) &
    & + ( coup(1) * coup(4) * Conjg(coup(8) * coup(5))                         &
    &   + coup(2) * coup(3) * Conjg(coup(7) * coup(6)) ) * I_ChiaChibInterf(2) &
    & + ( coup(1) * coup(4) * Conjg(coup(8) * coup(6))                         &
    &   + coup(2) * coup(3) * Conjg(coup(7) * coup(5)) ) * I_ChiaChibInterf(3) &
    & + ( coup(1) * coup(4) * Conjg(coup(7) * coup(5))                         &
    &   + coup(2) * coup(3) * Conjg(coup(8) * coup(6)) ) * I_ChiaChibInterf(4) &
    & + ( coup(1) * coup(3) * Conjg(coup(8) * coup(5))                         &
    &   + coup(2) * coup(4) * Conjg(coup(7) * coup(6)) ) * I_ChiaChibInterf(5) &
    & + ( coup(1) * coup(3) * Conjg(coup(7) * coup(6))                         &
    &   + coup(2) * coup(4) * Conjg(coup(8) * coup(5)) ) * I_ChiaChibInterf(6) &
    & + ( coup(1) * coup(3) * Conjg(coup(7) * coup(5))                         &
    &   + coup(2) * coup(4) * Conjg(coup(8) * coup(6)) ) * I_ChiaChibInterf(7) &
    & + ( coup(1) * coup(3) * Conjg(coup(8) * coup(6))                         &
    &   + coup(2) * coup(4) * Conjg(coup(7) * coup(5)) ) * I_ChiaChibInterf(8)

  gam = 2._dp  *  oo512pi3  *  m_in  *  gam  

 End Subroutine IntegrateChiChiInterference


 Subroutine FFIntegrands(s, erg)
 Implicit None

  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(6)

  Integer :: i1
  Real(dp) :: x, x2, fxrarb

  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2

   fxrarb = Sqrt(x2-4._dp*rkt2) * kappa(1._dp-x+rkt2,rmt2,rjt2)               &
     & / ( (1._dp - x + rkt2 - rat2)                                          &
     &   * (1._dp - x + rkt2 - rbt2) * (1._dp - x + rkt2)**2  ) 
     
   erg(2) = erg(2) + (x - 2._dp * rkt2) * (Rsquared - x ) * fxrarb
   fxrarb = (1._dp - x + rkt2) * fxrarb 
   erg(1) = erg(1) + (x - 2._dp * rkt2) * (Rsquared - x ) * fxrarb
   erg(3) = erg(3) + (x - 2._dp * rkt2) * fxrarb
   erg(4) = erg(4) + (Rsquared - x) * fxrarb
   erg(5) = erg(5) + fxrarb
   erg(6) = erg(6) + (1._dp - x + rkt2) * fxrarb 
                          !^^the square gets lost in the redefinition of fxrarb
  End Do

 End Subroutine FFIntegrands
  

 Subroutine IntegrateFF(mass , m_in, r_out, coup, smin, smax, eps, gam)
 Implicit None
  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam
  
  Real(dp) :: I_FaFb(6),I_FaFbm(6)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3) 

  rat2 = (mass(1) / m_in )**2
  rbt2 = (mass(2) / m_in )**2

  Rsquared = 1._dp - rjt2 + rkt2 + rmt2   

  Call DgaussInt(FFIntegrands,6,smin,smax,I_FaFb(:),eps)

  I_FaFbm(1) = I_FaFb(1)
  I_FaFbm(2) = Sqrt(rat2 * rbt2) * I_FaFb(2)
  I_FaFbm(3) = 4._dp * Sqrt(rmt2 * rbt2) * I_FaFb(3)
  I_FaFbm(4) = 4._dp * Sqrt(rkt2 * rbt2) * I_FaFb(4)
  I_FaFbm(5) = 8._dp * Sqrt(rkt2 * rmt2 * rat2 * rbt2) * I_FaFb(5)
  I_FaFbm(6) = 8._dp * Sqrt(rkt2 * rmt2) * I_FaFb(6)


  gam = ( Conjg( coup(1) * coup(2) ) * coup(3) * coup(4)                    &
    &   + Conjg( coup(5) * coup(6) ) * coup(7) * coup(8) ) * I_FaFbm(1)     &
    & + ( Conjg( coup(1) * coup(6) ) * coup(7) * coup(4)                    &
    &   + Conjg( coup(5) * coup(2) ) * coup(3) * coup(8) ) * I_FaFbm(2)     &
    & + Real( Conjg( coup(1) * coup(6) ) * coup(3) * coup(4)                &
    &       + Conjg( coup(5) * coup(2) ) * coup(7) * coup(8) ) * I_FaFbm(3) &
    & - Real( Conjg( coup(1) * coup(2) ) * coup(3) * coup(8)                &
    &       + Conjg( coup(5) * coup(6) ) * coup(7) * coup(4) ) * I_FaFbm(4) &
    & - Real( Conjg( coup(1) * coup(2) ) * coup(7) * coup(8) )* I_FaFbm(5)  &
    & - Real( Conjg( coup(1) * coup(6) ) * coup(3) * coup(8) )* I_FaFbm(6)

  gam = oo512pi3 * m_in * gam

 End Subroutine IntegrateFF



 Subroutine FFIntegrandsLM(s, erg)
 Implicit None

  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(6)

  Integer :: i1
  Real(dp) :: x, x2, fxrarb, gxrarb, yplus, yminus

  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2

   yplus = ( (2-x) * (Rsquared - x)                            &
       &   + Sqrt(x2-4._dp*rkt2) * kappa(1-x+rkt2,rmt2,rjt2) ) &
       & / ( 2._dp * ( 1._dp - x + rkt2 ) )
      
   yminus = ( (2-x) * (Rsquared - x)                            &
        &   - Sqrt(x2-4._dp*rkt2) * kappa(1-x+rkt2,rmt2,rjt2) ) &
        & / ( 2._dp * ( 1._dp - x + rkt2 ) )
     
   fxrarb =      (yplus - yminus)   &
      &     /  ( (1._dp - x + rkt2 - rat2) * (1._dp - x + rkt2 - rbt2)   )
   
   gxrarb =       (yplus**2 - yminus**2)   &
      &      /( 2._dp * (1._dp - x + rkt2 - rat2) * (1._dp - x + rkt2 - rbt2) )
   
   erg(1) = erg(1) + (rkt2 - 1._dp)*(rjt2- rkt2 -rmt2 + x - 1._dp)*fxrarb &
        & + (x-rkt2-1._dp)*gxrarb
   erg(2) = erg(2) + (rjt2 - rkt2 - rmt2 + x - 1._dp) * fxrarb + gxrarb
   erg(3) = erg(3) + (x - 2._dp*rkt2)*fxrarb
   erg(4) = erg(4) + (rjt2-rkt2-rmt2+x-1)*fxrarb
   erg(5) = erg(5) + fxrarb
   erg(6) = erg(6) + (rkt2 - x + 1._dp)*fxrarb
  
  End Do

 End Subroutine FFIntegrandsLM
  
  
 Subroutine IntegrateFFLM( mass , m_in, r_out, coup, smin,smax,eps,gam)
  Implicit None
  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam
  
  Real(dp) :: I_FaFbs(8), I_FaFb(6)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3) 

  rat2 = (mass(1) / m_in )**2
  rbt2 = (mass(2) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2

  Call DgaussInt(FFIntegrandsLM,6,smin,smax,I_FaFb(:),eps)
   I_FaFbs(1) = I_FaFb(1)
   I_FaFbs(2) = Sqrt(rat2 * rbt2) * I_FaFb(2)
   I_FaFbs(3) = Sqrt(rat2 * rmt2) * I_FaFb(3)
   I_FaFbs(4) = Sqrt(rbt2 * rmt2) * I_FaFb(3)
   I_FaFbs(5) = Sqrt(rat2 * rkt2) * I_FaFb(4)
   I_FaFbs(6) = Sqrt(rbt2 * rkt2) * I_FaFb(4)
   I_FaFbs(7) = - 2._dp * Sqrt(rkt2 * rmt2 * rat2 * rbt2) * I_FaFb(5)
   I_FaFbs(8) = - 2._dp * Sqrt(rkt2 * rmt2) * I_FaFb(6)
 
 
  gam = ( coup(1) * coup(4) * Conjg(coup(8) * coup(5))                 &
    &   + coup(2) * coup(3) * Conjg(coup(7) * coup(6)) ) * I_FaFbs(1)  &
    & + ( coup(1) * coup(3) * Conjg(coup(7) * coup(5))                 &
    &   + coup(2) * coup(4) * Conjg(coup(8) * coup(6)) ) * I_FaFbs(2)  &
    & + ( coup(1) * coup(3) * Conjg(coup(8) * coup(5))                 &
    &   + coup(2) * coup(4) * Conjg(coup(7) * coup(6)) ) * I_FaFbs(3)  &
    & + ( coup(1) * coup(4) * Conjg(coup(7) * coup(5))                 &
    &   + coup(2) * coup(3) * Conjg(coup(8) * coup(6)) ) * I_FaFbs(4)  &
    & + ( coup(1) * coup(3) * Conjg(coup(7) * coup(6))                 &
    &   + coup(2) * coup(4) * Conjg(coup(8) * coup(5)) ) * I_FaFbs(5)  &  
    & + ( coup(1) * coup(4) * Conjg(coup(8) * coup(6))                 &
    &   + coup(2) * coup(3) * Conjg(coup(7) * coup(5)) ) * I_FaFbs(6)  &  
    & + ( coup(1) * coup(3) * Conjg(coup(8) * coup(6))                 &
    &   + coup(2) * coup(4) * Conjg(coup(7) * coup(5)) ) * I_FaFbs(7)  &
    & + ( coup(1) * coup(4) * Conjg(coup(7) * coup(6))                 &
    &   + coup(2) * coup(3) * Conjg(coup(8) * coup(5)) ) * I_FaFbs(8)

  gam = 2._dp*oo512pi3 * m_in * gam

 End Subroutine IntegrateFFLM

 
 Subroutine SaSaIntegrands(s,erg)
  Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1 
 
  Real(dp) :: x, x2, yplus, yminus, iOnexra, iTwoxra, l1abbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
   
   yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                           &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp -x+rkt2,rmt2,rjt2))
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp -x+rkt2,rmt2,rjt2))

   l1abbr =  Log( Abs( ( yplus  + x - 1._dp + rjt2 - rat2 )         &
          &            / ( yminus + x - 1._dp + rjt2 - rat2 )   ) )

   iOnexra =  (yplus - yminus)                            &
           & / ( (yplus + x - 1._dp + rjt2 - rat2 )       &
           &   * (yminus + x - 1._dp + rjt2 - rat2 ) )
      
   iTwoxra = l1abbr + (1._dp - x - rjt2 + rat2)* iOnexra

   erg(1) = erg(1) + iTwoxra + (x - 1._dp +rjt2 - rmt2 - rkt2)*iOnexra
   erg(2) = erg(2) + iOnexra

  End Do 

 End Subroutine SaSaIntegrands


 Subroutine IntegrateSaSa(mass,m_in, r_out, coup,smin,smax,eps,gam)
  Implicit None
  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam
  
  Real(dp) :: I_SaSa(2)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)

  rat2 = (mass(1) / m_in)**2
  rbt2 = (mass(2) / m_in)**2
  
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2

  Call DgaussInt(SaSaIntegrands,2,smin,smax,I_SaSa(:),eps)

  gam = (coup(2)*Conjg(coup(5)) + coup(3)*Conjg(coup(6))) * I_SaSa(1)  &
    & - (coup(2)*Conjg(coup(6)) + coup(3)*Conjg(coup(5)))               &
    &   * 2._dp * Sqrt(rmt2*rkt2) * I_SaSa(2)

  gam = 2._dp * oo512pi3 / m_in * coup(1)*Conjg(coup(4))*  gam

 End Subroutine IntegrateSaSa


 Subroutine SaSbIntegrands (s,erg)
  Implicit None  
  
   Real(dp), Intent(in) :: s(2)
   Real(dp), Intent(out) :: erg(2)

   Integer :: i1 
 
   Real(dp) :: x, x2, yplus, yminus, hOnexrarb, hTwoxrarb, l1abbr, l2abbr
 
   erg = 0._dp

   Do i1=1,2
    x = s(i1)
    x2 = x**2
   
    yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                         &
    &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)                &
    &                                 * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                         &
    &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)                &
    &                                 * kappa(1._dp-x+rkt2,rmt2,rjt2) )
      
   l1abbr = Log( Abs( ( yplus  + x - 1._dp + rjt2 - rat2 )         &
          &           / ( yminus + x - 1._dp + rjt2 - rat2 )   ) ) 
 
   l2abbr = Log( Abs( ( yplus  + x - 1._dp + rjt2 - rbt2 )         &
           &          / ( yminus + x - 1._dp + rjt2 - rbt2 )   ) )  
            
   hOnexrarb =  ( l1abbr - l2abbr  ) / (rat2 - rbt2)
   
   hTwoxrarb =   l2abbr + (1._dp - x -rjt2 + rat2 )*hOnexrarb
      
   erg(1) = erg(1) +  hTwoxrarb + (x - 1._dp +rjt2 - rmt2 - rkt2)*hOnexrarb
   erg(2) = erg(2) +  hOnexrarb
   
  End Do 

 End Subroutine SaSbIntegrands


 Subroutine IntegrateSaSb(mass,m_in, r_out, coup,smin,smax,eps,gam)
  Implicit None
  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam
  
  Real(dp) :: I_SaSb(2)
  
  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)

  rat2 = (mass(1) / m_in)**2
  rbt2 = (mass(2) / m_in)**2
  
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2  

  Call DgaussInt(SaSbIntegrands,2,smin,smax,I_SaSb(:),eps)

  gam = (coup(2)*Conjg(coup(5)) + coup(3)*Conjg(coup(6))) * I_SaSb(1)  &
    & - (coup(2)*Conjg(coup(6)) + coup(3)*Conjg(coup(5)))              &
    &   * 2._dp * Sqrt(rmt2*rkt2) * I_SaSb(2)
  
  gam = 2._dp * oo512pi3 / m_in * coup(1) * Conjg(coup(4)) * gam
 
 End Subroutine IntegrateSaSb
 
 
 Subroutine SFIntegrands(s,erg)
 Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1 

  Real(dp) :: x, x2, yplus, yminus, gOnexrarb, gTwoxrarb, l1abbr, kabbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
 
   yplus = 1._dp / ( 2._dp * ( 1 - x + rkt2 ) )                               &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               & 
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          & 
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               & 
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )

   l1abbr =      Log( Abs(  ( yplus  + x - 1._dp + rjt2 - rbt2 )         &
            &             / ( yminus + x - 1._dp + rjt2 - rbt2 )   ) ) 
            
   kabbr = (   yplus - yminus  ) / ( 1._dp - x + rkt2 -rat2 )

   gOnexrarb =  l1abbr / ( 1._dp - x + rkt2 - rat2  )

   gTwoxrarb =  kabbr + ( 1._dp - x -rjt2 + rbt2 )*gOnexrarb

   erg(1) = erg(1) + (-rjt2 + rkt2 + rmt2 - x + 1._dp)*gOnexrarb
   erg(2) = erg(2) + (-x + 2._dp*rkt2)*gOnexrarb
   erg(3) = erg(3) + gOnexrarb
   erg(4) = erg(4) + gTwoxrarb + (rjt2 - rkt2 - rmt2 + x - 1._dp)*gOnexrarb
   
  End Do  

 End Subroutine SFIntegrands


 Subroutine IntegrateSF(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_SF(4)


  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)
  
  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2  
   
  Call DGaussInt(SFIntegrands,4,smin,smax,I_SF(:),eps)

  I_SF(1) = I_SF(1)
  I_SF(2) = I_SF(2)
  I_SF(3) = - I_SF(3)  
  I_SF(4) = I_SF(4)  
  
  gam = ( coup(1) * coup(4) * Conjg(coup(7))                          &
    &   + coup(2) * coup(3) * Conjg(coup(6)) ) * Sqrt(rkt2) * I_SF(1) &
    & + ( coup(1) * coup(4) * Conjg(coup(6))                          &
    &   + coup(2) * coup(3) * Conjg(coup(7)) ) * Sqrt(rmt2) * I_SF(2) &
    & - ( coup(1) * coup(3) * Conjg(coup(7))                          &
    &   + coup(2) * coup(4) * Conjg(coup(6)) )                        &
    &                        * 2._dp * Sqrt(rkt2*rmt2*rat2) * I_SF(3) &
    & + ( coup(1) * coup(3) * Conjg(coup(6))                          &
    &   + coup(2) * coup(4) * Conjg(coup(7)) ) * Sqrt(rat2) * I_SF(4)

  gam = 2._dp * oo512pi3 * Conjg(coup(5)) * gam

 End Subroutine IntegrateSF
 
 Subroutine VFIntegrands(s,erg)
 Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(4)

  Integer :: i1 
 
  Real(dp) :: x, x2, yplus, yminus, gOnexrarb, gTwoxrarb, l1abbr, kabbr
  
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
   
   yplus = 1._dp / ( 2._dp * ( 1 - x + rkt2 ) )                               &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) ) &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               & 
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   
   l1abbr =      Log( Abs(  ( yplus  + x - 1._dp + rjt2 - rbt2 )         &
            &             / ( yminus + x - 1._dp + rjt2 - rbt2 )   ) ) 
            
   kabbr = (   yplus - yminus  ) / ( 1._dp - x + rkt2 -rat2  )

   gOnexrarb =  l1abbr / ( 1._dp - x + rkt2 - rat2  )

   gTwoxrarb =  kabbr + ( 1._dp - x -rjt2 + rbt2 )*gOnexrarb

   erg(1) = erg(1) -   2._dp*(rkt2 - x + 1._dp)*gTwoxrarb                   &
         &         + ( rjt2 *(rkt2 - 2._dp ) -  rkt2 *rkt2                  &
         &         +   rkt2 *(rmt2 + x + 1._dp)                             &
         &         -   rmt2*(x - 2._dp) -2._dp*x  +2._dp   )*gOnexrarb
   erg(2) = erg(2) + gTwoxrarb + (rjt2 + rkt2 - rmt2 - x - 1._dp)*gOnexrarb
   erg(3) = erg(3) + (rjt2 + rkt2 - rmt2 -2._dp*x + 3._dp)*gOnexrarb
   erg(4) = erg(4) + gTwoxrarb + (1._dp - rjt2 + rkt2 - rmt2 - x)*gOnexrarb   

  End Do  

 End Subroutine VFIntegrands
 
 
 Subroutine IntegrateVF(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None  

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_VF(4)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)

  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2   

  Call DGaussInt(VFIntegrands,4,smin,smax,I_VF(:),eps)

  I_VF(2)=I_VF(2)
  I_VF(3)=-I_VF(3)  
  I_VF(4)=I_VF(4)

  gam = ( coup(1) * coup(4) * Conjg(coup(7))                               &
    &   + coup(2) * coup(3) * Conjg(coup(6)) ) * I_VF(1)                   &
    & + ( coup(1) * coup(3) * Conjg(coup(7))                               &
    &   + coup(2) * coup(4) * Conjg(coup(6)) ) * Sqrt(rat2*rmt2) * I_VF(2) &
    & - ( coup(1) * coup(4) * Conjg(coup(6))                               &
    &   + coup(2) * coup(3) * Conjg(coup(7)) ) * Sqrt(rkt2*rmt2) * I_VF(3) &
    & + ( coup(1) * coup(3) * Conjg(coup(6))                               &
        + coup(2) * coup(4) * Conjg(coup(7)) ) * Sqrt(rkt2*rat2) * I_VF(4)      

  gam = - 2._dp * oo512pi3 * Conjg(coup(5)) * m_in * gam  

 End Subroutine IntegrateVF
 
 Subroutine VSIntegrands(s,erg)
  Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1 
 
  Real(dp) :: x, x2, yplus, yminus, hOnexrarb, hTwoxrarb, l1abbr, l2abbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
   
   yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                           &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
 
   l1abbr =      Log( Abs(  ( yplus  + x - 1._dp + rjt2 - rat2 )         &
            &             / ( yminus + x - 1._dp + rjt2 - rat2 )   ) ) 
 
   l2abbr =      Log( Abs(  ( yplus  + x - 1._dp + rjt2 - rbt2 )         &
            &             / ( yminus + x - 1._dp + rjt2 - rbt2 )   ) )  
            
   hOnexrarb = ( l1abbr - l2abbr )  / (rat2 - rbt2)
   
   hTwoxrarb = l2abbr + (1._dp - x -rjt2 + rat2 )*hOnexrarb
 
   erg(1) = erg(1) + hTwoxrarb + (rkt2 - rjt2 - rmt2 - x +1._dp)*hOnexrarb
   erg(2) = erg(2) + hTwoxrarb + (rjt2 - rmt2 - x - 1._dp +rkt2)*hOnexrarb
 
  End Do 

 End Subroutine VSIntegrands
  
 
 Subroutine IntegrateVS(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None  
  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_VS(2)
  
  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)  
  
  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in)**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2   

  Call DGaussInt(VSIntegrands,2,smin,smax,I_VS(:),eps)

  gam = (coup(2)*Conjg(coup(5)) + coup(3)*Conjg(coup(6)))*Sqrt(rkt2)*I_VS(1) &
    & + (coup(2)*Conjg(coup(6)) + coup(3)*Conjg(coup(5)))*Sqrt(rmt2)*I_VS(2)         

  gam = 2._dp * oo512pi3 * coup(1) * Conjg(coup(4)) * gam

 End Subroutine IntegrateVS
 

 Subroutine VSIntegrandsGoldstone(s,erg)
 Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1 

  Real(dp) :: x, x2, yplus, yminus, iOnexra, iTwoxra, l1abbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
   
   yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                           &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
 
   l1abbr =      Log( Abs(  ( yplus  + x - 1._dp + rjt2 - rat2 )         &
            &             / ( yminus + x - 1._dp + rjt2 - rat2 )   ) ) 
          
   iOnexra =                           (yplus - yminus)                       &
    & / ((yplus + x - 1._dp + rjt2 - rat2)*(yminus + x - 1._dp + rjt2 - rat2))
      
   iTwoxra = l1abbr + (1._dp - x - rjt2 + rat2)* iOnexra   
  
   erg(1) = erg(1) + iTwoxra + (rkt2 - rjt2 - rmt2 - x +1._dp)*iOnexra
   erg(2) = erg(2) + iTwoxra + (rjt2 - rmt2 - x - 1._dp +rkt2)*iOnexra

  End Do 

 End Subroutine VSIntegrandsGoldstone
 
  
 Subroutine IntegrateVSGoldstone(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None  

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_VS(2)

  
  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)  
  
  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in)**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2  

  Call DGaussInt(VSIntegrandsGoldstone,2,smin,smax,I_VS(:),eps)

  gam = (coup(2)*Conjg(coup(5)) + coup(3)*Conjg(coup(6)))*Sqrt(rkt2)*I_VS(1) &
    & + (coup(2)*Conjg(coup(6)) + coup(3)*Conjg(coup(5)))*Sqrt(rmt2)*I_VS(2)

  gam = 2._dp * oo512pi3 * coup(1) * Conjg(coup(4)) * gam

 End Subroutine IntegrateVSGoldstone
 
 Subroutine VVIntegrands (s,erg)
  Implicit None  
  
  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: x, x2, yplus, yminus, iOnexra, iTwoxra, l1abbr
 
  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2
   
   yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                           &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) ) 
    
   l1abbr = Log( Abs( ( yplus  + x - 1._dp + rjt2 - rat2 )         &
          &         / ( yminus + x - 1._dp + rjt2 - rat2 )   ) ) 
          
   iOnexra = (yplus - yminus)                                 &
           & / ( (yplus  + x - 1._dp + rjt2 - rat2 )          &
           &   *  (yminus + x - 1._dp + rjt2 - rat2 ) )
      
   iTwoxra = l1abbr + (1._dp - x - rjt2 + rat2)* iOnexra   

   erg(1) = erg(1) +   ( -3._dp * rkt2 + rmt2 + 4._dp * x - 4._dp) *iTwoxra   &
                &  +   (  rjt2 * (rkt2 + rmt2 - 4._dp)                        &
                &  -         rkt2*rkt2                                        &
                &  +      rkt2 * (2._dp * rmt2 + x + 3._dp)                   &
                &  -      rmt2 * (rmt2 + 3._dp * x - 3._dp)                   &
                &  -       4._dp * x + 4._dp                      )*iOnexra
                
   erg(2) = erg(2) + iTwoxra + (x - rjt2 - 3._dp)*iOnexra
   
  End Do
 
 End Subroutine VVIntegrands
 
 
 Subroutine IntegrateVV(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None  

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_VV(2)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)

  rat2 = ( mass(1) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2 

  Call DGaussInt(VVIntegrands,2,smin,smax,I_VV(:),eps)
  
  gam = (coup(2)*Conjg(coup(2)) + coup(3)*Conjg(coup(3)) ) * I_VV(1)  &
    & + Real(coup(2)*coup(3))  * 4._dp*Sqrt( rmt2 * rkt2 )*I_VV(2) 
  
  gam = 2._dp * oo512pi3 * m_in * coup(1) *Conjg(coup(1)) *gam


 End Subroutine IntegrateVV
 
 Complex(dp) Function I2dsmg1tmg2(s, tmin, tmax, m12, g1, m22, g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     1                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: s, tmin, tmax, m12, m22, g1, g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2dsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0) Then
   ReT = Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   ReT = 0.5_dp * Log(( g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    ImT = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    ImT = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    ImT = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
  End If

  I2dsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

  End Function I2dsmg1tmg2


 Complex(dp) Function I2dtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     1                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3),g12,g22, x1, x2, y1, y2

  I2dtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2dtm12g12 = ( Log(Abs( (tmax-m12) / (tmin-m12) ) )          &
              & - Log(Abs( (tmax-m22) / (tmin-m22) ) ) ) / (m12 - m22)
  Else
   g12 = g1 * g1
   g22 = g2 * g2
   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * (m12-m22)
   fa(3) = (g1+g2)
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2dtm12g12 =                                                            &
     &   fa(2) * ( Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )    &
     &          - Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22) )  )  &
     & - fa(3) * ( y1 + y2  )
   I2dtm12g12 = I2dtm12g12 / fa(1)

  End If

 End Function I2dtm12g12


 Real(dp) Function I2dtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (            1
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin, tmax, m12, g12
  Real(dp) :: g1, x1, x2

  If (g12.Le.0._dp) Then
   I2dtm1g1 = 1._dp / (tmin-m12) - 1._dp / (tmax-m12)
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   I2dtm1g1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) I2dtm1g1 = I2dtm1g1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) I2dtm1g1 = I2dtm1g1 - pi
   I2dtm1g1 = I2dtm1g1 / g1
  End If

  End Function I2dtm1g1


 Complex(dp) Function I2t2dsmg1tmg2(s,tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t^2                 )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s,tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2t2dsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0._dp) Then
   ReT = 0.5_dp * ( tmax*(2._dp*m22+tmax) - tmin*(2._dp*m22+tmin)) &
     & + m22**2 * Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   fa(1) = Log((g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   ReT = 0.5_dp * ( tmax*(2._dp*m22+tmax) - tmin*(2._dp*m22+tmin)) &
     & + 0.5_dp * (m22**2 - g2**2) * fa(1)
   ImT = -g2 * (tmax - tmin + m22 * fa(1))
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    fa(2) = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    fa(2) = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    fa(2) = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
   ReT = ReT + 2._dp * g2 * m22 * fa(2)
   ImT = ImT + (m22**2 - g2**2) * fa(2)
  End If

  I2t2dsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

 End Function I2t2dsmg1tmg2


 Complex(dp) Function I2tdtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(4),g12,g22, x1, x2, y1, y2

  I2tdtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2tdtm12g12 = ( m12 * Log(Abs( (tmax-m12) / (tmin-m12) ) )      &
               & - m22 * Log(Abs( (tmax-m22) / (tmin-m22) ) ) ) / (m12 - m22)
  Else
   g12 = g1 * g1
   g22 = g2 * g2
   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * (m12 * (m12-m22) + g1 *(g1+g2) )
   fa(3) = 0.5_dp * (m22 * (m22-m12) + g2 *(g1+g2) )
   fa(4) = (g1*m22+g2*m12)
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2tdtm12g12 = fa(2) * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) ) &
             & + fa(3) * Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22) ) &
             & - fa(4) * ( y1 + y2  )
   I2tdtm12g12 = I2tdtm12g12 / fa(1)
  End If

 End Function I2tdtm12g12


 Complex(dp) Function I2t2dtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                    t^2                  )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: m14,m24,fa(5),g12,g22, x1, x2, y1, y2

  I2t2dtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2t2dtm12g12 = tmax - tmin                                           &
              & + ( m12**2 * Log(Abs( (tmax-m12) / (tmin-m12) ) )       &
              &   - m22**2 * Log(Abs( (tmax-m22) / (tmin-m22) ))) / (m12 - m22)
  Else
   m14 = m12 * m12
   m24 = m22 * m22
   g12 = g1 * g1
   g22 = g2 * g2

   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * ( (m14-g12)*(m12-m22) + 2._dp* m12* g1* (g1+g2) )
   fa(3) = -( (m14-g12)*(g1+g2) - 2._dp* m12* g1* (m12-m22) )
   fa(4) = 0.5_dp * ( (m24-g22)*(m22-m12) + 2._dp* m22* g2* (g1+g2) )
   fa(5) = -( (m24-g22)*(g1+g2) + 2._dp* m22* g2* (m12-m22) )
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2t2dtm12g12 = tmax - tmin                                                 &
              & + fa(2) * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12)) &
              & + fa(3) * y1                                                  &
              & + fa(4) * Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22)) &
              & + fa(5) * y2
   I2t2dtm12g12 = I2t2dtm12g12 / fa(1)
  End If

 End Function I2t2dtm12g12



 Real(dp) Function I2t2dtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (           t^2
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp),Intent(in) :: tmin, tmax, m12, g12
  Real(dp) :: g1, x1, x2, y

  If (g12.Le.0._dp) Then
   I2t2dtm1g1 = m12**2 * ( 1._dp/(tmin-m12) - 1._dp/(tmax-m12) )    &
            & + 2._dp * m12 * Log(Abs( (tmax-m12) / (tmin-m12) ) )  &
            & + tmax - tmin
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   y = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y - pi

   I2t2dtm1g1 = (m12**2-g12) * y / g1 &
            & + m12 * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )   &
            & + tmax - tmin
  End If

 End Function I2t2dtm1g1


 Complex(dp) Function I2tdsmg1tmg2(s,tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s,tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2tdsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0._dp) Then
   ReT = tmax - tmin + m22 * Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   fa(1) = 0.5_dp * Log((g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   ReT = tmax - tmin + m22 * fa(1)
   ImT = -g2 * fa(1)
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    fa(2) = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    fa(2) = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    fa(2) = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
   ReT = ReT + g2 * fa(2)
   ImT = ImT + m22 * fa(2)
  End If

  I2tdsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

 End Function I2tdsmg1tmg2


 Real(dp) Function I2tdtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (            t
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: tmin,tmax,m12,g12
  Real(dp) :: g1, x1, x2, y


  If (g12.Le.0._dp) Then
   I2tdtm1g1 = m12 / (tmin-m12) - m12 / (tmax-m12)   &
           & + Log(Abs( (tmax-m12) / (tmin-m12) ) )
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   y = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y - pi

   I2tdtm1g1 = m12 * y / g1          &
           & + 0.5_dp * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )
  End If

 End Function I2tdtm1g1

!---------------------------------------------------------------------
! integrals for multiple vector bosons in 3-body decays of fermions
! provided by L. Mitzka, 04.03.2015
!---------------------------------------------------------------------
  Subroutine VaVbIntegrands (s,erg)
  Implicit None

  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(2)

  Integer :: i1
  Real(dp) :: x, x2, yplus, yminus, hOnexrarb, hTwoxrarb, l1abbr, l2abbr

  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2

   yplus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                           &
      &    * ( (2._dp-x) * (Rsquared - x) + Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )
   yminus = 1._dp / ( 2._dp * ( 1._dp - x + rkt2 ) )                          &
      &    * ( (2._dp-x) * (Rsquared - x) - Sqrt(x2-4._dp*rkt2)               &
      &                                     * kappa(1._dp-x+rkt2,rmt2,rjt2) )

   l1abbr = Log( Abs( ( yplus  + x - 1._dp + rjt2 - rat2 )         &
          &           / ( yminus + x - 1._dp + rjt2 - rat2 )   ) )

   l2abbr = Log( Abs( ( yplus  + x - 1._dp + rjt2 - rbt2 )         &
           &          / ( yminus + x - 1._dp + rjt2 - rbt2 )   ) )

   hOnexrarb =  ( l1abbr - l2abbr  ) / (rat2 - rbt2)

   hTwoxrarb =   l2abbr + (1._dp - x -rjt2 + rat2 )*hOnexrarb

   erg(1) = erg(1) +   ( -3._dp * rkt2 + rmt2 + 4._dp * x - 4._dp) *hTwoxrarb &
                &  +   (  rjt2 * (rkt2 + rmt2 - 4._dp)                        &
                &  -         rkt2*rkt2                                        &
                &  +      rkt2 * (2._dp * rmt2 + x + 3._dp)                   &
                &  -      rmt2 * (rmt2 + 3._dp * x - 3._dp)                   &
                &  -       4._dp * x + 4._dp                      )*hOnexrarb

   erg(2) = erg(2) + hTwoxrarb + (x - rjt2 - 3._dp)*hOnexrarb

  End Do

 End Subroutine VaVbIntegrands


 Subroutine IntegrateVaVb(mass,m_in, r_out, coup,smin,smax,eps,gam)
 Implicit None

  Real(dp), Intent(in)    :: mass(:)
  Real(dp), Intent(in)    :: m_in
  Real(dp), Intent(in)    :: r_out(:)
  Complex(dp), Intent(in) :: coup(:)
  Real(dp), Intent(in)    :: smin
  Real(dp), Intent(in)    :: smax
  Real(dp), Intent(in)    :: eps
  Complex(dp), Intent(out)   :: gam

  Real(dp) :: I_VaVb(2)

  rjt2 = r_out(1)
  rkt2 = r_out(2)
  rmt2 = r_out(3)

  rat2 = ( mass(1) / m_in )**2
  rbt2 = ( mass(2) / m_in )**2
  Rsquared = 1._dp - rjt2 + rkt2 + rmt2

  Call DGaussInt(VaVbIntegrands,2,smin,smax,I_VaVb(:),eps)

  gam = (coup(2)*Conjg(coup(5)) + coup(3)*Conjg(coup(6)) )  * I_VaVb(1)  &
    & + (coup(2)*Conjg(coup(6)) + coup(3)*Conjg(coup(5)) )  * 2._dp      &
    &    *Sqrt( rmt2 * rkt2 )*I_VaVb(2)

  gam = 2._dp * oo512pi3 * m_in * coup(1) *Conjg(coup(4)) *gam


 End Subroutine IntegrateVaVb

End Module ThreeBodyPhaseSpaceS

