Module Slepton3BodyDecays

Use Control
Use Mathematics
Use ThreeBodyPhaseSpaceS

Contains


 Subroutine Slepton_3body(n_in, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u  &
   & , n_c, n_n, n_W, id_W, n_Z, id_Z, n_Sl, n_Snu, n_S0, n_P0, n_Spm         &
   & , Slept, mZ, mW, mf_l, mf_u, mf_d, Chi0, ChiPm, P0, S0, Spm, Sneut       &
   & , c_P0SlSl, c_S0SlSl, c_SlSlZ, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn    &
   & , c_CNuSl_L, c_CNuSl_R, L_nu, R_nu, L_e, R_e, L_u, R_u, L_d, R_d, c_LNuW &
   & , c_DUW, c_LLP0_L, c_LLP0_R, c_UUP0_L, c_UUP0_R, c_DDP0_L, c_DDP0_R      &
   & , c_LLS0_L, c_LLS0_R, c_UUS0_L, c_UUS0_R, c_DDS0_L, c_DDS0_R             &
   & , c_CLSn_L, c_CLSn_R, c_NuNSn_L, c_NuNSn_R, c_SmpLNu_L, c_SmpLNu_R       &
   & , c_SmpDU_L, c_SmpDU_R, GenerationMixing, epsI)
 !-----------------------------------------------------------------------------
 ! calculates the various slepton 3-body decays using the routines provided
 ! by Lukas Mitzka
 ! written by Werner Porod, 17.08.2012
 !-----------------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_in, n_nu, n_l, n_d, n_u, n_Z, n_W, n_snu  &
       & , n_sl, n_n, n_c, n_s0, n_p0, n_Spm
  Integer, Intent(in) :: id_nu(:), id_l(:), id_d(:), id_u(:), id_Z(:), id_W(:)
  Real(dp), Intent(in) :: mZ, L_nu, R_nu, mf_l(3), L_e, R_e, mf_u(3), L_u, R_u &
     & , mf_d(3), L_d, R_d, mW, epsI
  Complex(dp), Intent(in) :: c_S0SlSl(:,:,:), c_LLS0_L(:,:,:), c_LLS0_R(:,:,:) &
     & , c_P0SlSl(:,:,:), c_LLP0_L(:,:,:), c_LLP0_R(:,:,:), c_LNSl_L(:,:,:)    &
     & , c_LNSl_R(:,:,:), c_SlSlZ(:,:,:), c_UUS0_L(:,:,:), c_UUS0_R(:,:,:)     &
     & , c_UUP0_L(:,:,:), c_UUP0_R(:,:,:), c_DDS0_L(:,:,:), c_DDS0_R(:,:,:)    &
     & , c_DDP0_L(:,:,:), c_DDP0_R(:,:,:), c_CNuSl_L(:,:,:), c_CNuSl_R(:,:,:)  &
     & , c_LNuW(:,:), c_SmpLNu_L(:,:,:), c_SmpLNu_R(:,:,:), c_SlSnW(:,:,:)     &
     & , c_SmpSlSn(:,:,:), c_NuNSn_L(:,:,:), c_NuNSn_R(:,:,:), c_CLSn_L(:,:,:) &
     & , c_CLSn_R(:,:,:), c_DUW(:,:), c_SmpDU_L(:,:,:), c_SmpDU_R(:,:,:)
  Logical, Intent(in) :: GenerationMixing
  Type(particle2), Intent(in) :: Spm(:), P0(:)
  Type(particle23), Intent(in) :: Sneut(:), ChiPm(:), Chi0(:), S0(:)
  Type(particle23), Intent(inout) :: Slept(:)

  Integer :: i_start, i_end, i_run, i2, i3, i4, i_c, i2_max
  Real(dp) :: mN(n_n), mSlept(n_sl), mSneut(n_snu), mC(n_c), mS0(n_S0), gam    &
     & , mSpm(n_Spm), mP0(n_P0), m_nu(3), gam2(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Slepton_3body'

  If (n_in.Lt.0) Then
   i_start = 1
   i_end = n_sl

  Else If ( (n_in.Ge.1).And.(n_in.Le.n_sl) ) Then 
   i_start = n_in 
   i_end = n_in

  Else
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write (ErrCan,*) 'Value of n_in out of range, (n_in,n_sl) = ',n_in,n_sl
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  !------------------------
  ! initialisation
  !------------------------
  m_nu = 0._dp

  mC = ChiPm%m
  mN = Chi0%m
  mSlept = Slept%m
  mSneut = Sneut%m
  mP0 = P0%m
  mS0 = S0%m
  mSpm = Spm%m

  Do i_run = i_start, i_end
   Slept(i_run)%gi3 = 0._dp
   Slept(i_run)%bi3 = 0._dp
   i_c = 1
   !---------------------------------------------------------------------------
   ! decays into another slepton and two leptons, it is easier to run over all
   ! sleptons and check if the decay is allowed than disthingishing between
   ! the cases with and without flavour violation
   !---------------------------------------------------------------------------
   If (GenerationMixing) Then
    i2_max = i_run - 1
   Else
    i2_max = n_sl
   End If

   If (mSlept(i_run).Lt.Abs(mN(1))) Then !
    Do i2=1,i2_max
     Do i3=1,3
      Do i4=1,3
       Call Slepton_to_Slepton_ll(i_run, i2, i3, i4, mSlept, mN, mZ, mS0, mP0, mf_l &
           & , c_LNSl_L, c_LNSl_R, c_SlSlZ, c_S0SlSl, c_P0SlSl, c_LLS0_L       &
           & , c_LLS0_R , c_LLP0_L, c_LLP0_R, L_e, R_e, epsI, gam2)
       Slept(i_run)%gi3(i_c) = gam2(1)
       Slept(i_run)%id3(i_c,1) = Slept(i2)%id
       Slept(i_run)%id3(i_c,2) = id_l(i3)
       Slept(i_run)%id3(i_c,3) = id_l(i4) + 1
       i_c = i_c + 1
       
       If(i3>=i4) Then
        Slept(i_run)%gi3(i_c) = gam2(2)
        Slept(i_run)%id3(i_c,1) = Slept(i2)%id + 1
        Slept(i_run)%id3(i_c,2) = id_l(i3)
        Slept(i_run)%id3(i_c,3) = id_l(i4)
        i_c = i_c + 1
       End If
      End Do

     End Do
    End Do 
   End If ! check if m_N > m_Sl

   !---------------------------------------------------
   ! to slepton + nu nu
   !---------------------------------------------------
   Do i2=1,i2_max
    Do i3=1,3
     Do i4=1,3
      Call Slepton_to_Slepton_nunu(i_run, i2, i3, i4, mSlept, m_nu, mC, mZ    &
          & , c_CNuSl_L, c_CNuSl_R, c_SlSlZ, L_nu, R_nu, epsI, gam)
       Slept(i_run)%gi3(i_c) = gam
       Slept(i_run)%id3(i_c,1) = Slept(i2)%id
       Slept(i_run)%id3(i_c,2) = id_nu(i3)
       Slept(i_run)%id3(i_c,3) = id_nu(i4) + 1
       i_c = i_c + 1
     End Do
    End Do
   End Do  

   !---------------------------------------------------
   ! to slepton + uu, currently only diagonal channels
   !---------------------------------------------------
   Do i2=1,i2_max
    Do i3=1,3
!     Do i4=1,3
     i4 = i3
     Call Slepton_to_Slepton_qq(i_run, i2, i3, i4, mSlept, mf_u, mZ, mS0   &
          & , mP0, c_SlSlZ, c_S0SlSl, c_P0SlSl, c_UUS0_L &
          & , c_UUS0_R , c_UUP0_L, c_UUP0_R, L_u, R_u, epsI, gam)
     Slept(i_run)%gi3(i_c) = gam
     Slept(i_run)%id3(i_c,1) = Slept(i2)%id
     Slept(i_run)%id3(i_c,2) = id_u(i3)
     Slept(i_run)%id3(i_c,3) = id_u(i4) + 1
     i_c = i_c + 1
!      end do
    End Do
   End Do
  
   !---------------------------------------------------
   ! to slepton + dd, currently only diagonal channels
   !---------------------------------------------------
   Do i2=1,i2_max
    Do i3=1,3
!     Do i4=1,3
     i4=i3
     Call Slepton_to_Slepton_qq(i_run, i2, i3, i4, mSlept, mf_d, mZ, mS0   &
          & , mP0, c_SlSlZ, c_S0SlSl, c_P0SlSl, c_DDS0_L &
          & , c_DDS0_R , c_DDP0_L, c_DDP0_R, L_d, R_d, epsI, gam)
     Slept(i_run)%gi3(i_c) = gam
     Slept(i_run)%id3(i_c,1) = Slept(i2)%id
     Slept(i_run)%id3(i_c,2) = id_d(i3)
     Slept(i_run)%id3(i_c,3) = id_d(i4) + 1
     i_c = i_c + 1
!      End Do
    End Do
   End Do

   !---------------------------------------------------
   ! to sneutrino + l nu
   !---------------------------------------------------
   Do i2=1,n_snu
    Do i3=1,3
     Do i4=1,3
      Call Slepton_to_Sneutrino_lnu(i_run, i2, i3, i4, mN, mSlept, mSneut     &
           & , mW, mSPm, mf_l, m_nu, c_LNSl_L, c_LNSl_R, c_NuNSn_L, c_NuNSn_R &
           & , c_SlSnW, c_SmpSlSn, c_LNuW(i3,i4), c_SmpLNu_L, c_SmpLNu_R, epsI, gam)

      Slept(i_run)%gi3(i_c)  = gam
      Slept(i_run)%id3(i_c,1) = Sneut(i2)%id 
      Slept(i_run)%id3(i_c,2) = id_l(i3)
      Slept(i_run)%id3(i_c,3) = id_nu(i4) + 1
      i_c = i_c + 1
     End Do
    End Do
   End Do

   !---------------------------------------------------
   ! to anti-sneutrino + l nu
   !---------------------------------------------------
   Do i2=1,n_snu
    Do i3=1,3
     Do i4=1,3
      Call Slepton_to_Antisneutrino_lnu(i_run, i2, i3, i4, mN, mSlept, mSneut &
             & , mC, mf_l, m_nu, c_LNSl_L, c_LNSl_R, c_NuNSn_L , c_NuNSn_R     &
             & , c_CLSn_L, c_CLSn_R, c_CNuSl_L, c_CNuSl_R, epsI, gam)

      Slept(i_run)%gi3(i_c)  = gam
      Slept(i_run)%id3(i_c,1) = Sneut(i2)%id +1
      Slept(i_run)%id3(i_c,2) = id_l(i3)
      Slept(i_run)%id3(i_c,3) = id_nu(i4)
      i_c = i_c + 1
     End Do
    End Do
   End Do

   !---------------------------------------------------
   ! to sneutrino + q bar(q)
   !---------------------------------------------------
   Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       Call Slepton_to_Sneutrino_qq(i_run, i2, i3, i4, mSlept, mSneut, mSpm, mW   &
         & , mf_d, mf_u, c_SmpDU_L, c_SmpDU_R, c_DUW, c_SmpSlSn, c_SlSnW, epsI &
         & , gam)
       Slept(i_run)%gi3(i_c) = gam
       Slept(i_run)%id3(i_c,1) = Sneut(i2)%id
       Slept(i_run)%id3(i_c,2) = id_d(i3) 
       Slept(i_run)%id3(i_c,3) = id_u(i4) + 1 
       i_c = i_c + 1
      End Do
     End Do
    End Do
  
   Slept(i_run)%g = Sum(Slept(i_run)%gi2) + Sum(Slept(i_run)%gi3)
   If (Slept(i_run)%g.Ne.0._dp) Then
    Slept(i_run)%bi2 = Slept(i_run)%gi2 / Slept(i_run)%g
    Slept(i_run)%bi3 = Slept(i_run)%gi3 / Slept(i_run)%g
   End If
   
  End Do ! i_run

  Iname = Iname - 1

 End Subroutine Slepton_3body

 Subroutine Slepton_to_Slepton_ll(i, j, k, m, mSl, mN, mZ, mH, mA, mf_l, C_L, C_R, &
                              &   C_Z, C_H, C_A, CL_H_L, CL_H_R, CL_A_L,     &
                              &   CL_A_R, CL_Z_L, CL_Z_R , eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of lepton m
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mH(:)    ! higgs masses 
  Real(dp), Intent(in) :: mA(:)    ! pseudo scalar higgs masses 
  Real(dp), Intent(in) :: mZ       ! suprise surprise Z-mass!!
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Complex(dp), Intent(in) :: C_R(:,:,:), C_L(:,:,:) ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_Z(:,:,:) ! couplings of sleptons to the Z 
  Complex(dp), Intent(in) :: C_H(:,:,:) ! couplings of sleptons to the scalar   
  Complex(dp), Intent(in) :: C_A(:,:,:) ! couplings of sleptons to the pseudosc
  Real(dp), Intent(in) :: CL_Z_L ! couplings of leptons to the Z 
  Real(dp), Intent(in) :: CL_Z_R ! couplings of leptons to the Z 
  Complex(dp), Intent(in) :: CL_H_L(:,:,:) ! left couplings of leptons to sc.
  Complex(dp), Intent(in) :: CL_A_L(:,:,:) ! left couplings of leptons to psc.
  Complex(dp), Intent(in) :: CL_H_R(:,:,:) ! right couplings of leptons to sc,     
  Complex(dp), Intent(in) :: CL_A_R(:,:,:) ! right couplings of leptons to psc.
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam(2)  ! partial widths: 1 .. same sign sleptons
                                   !                 2 .. opposite sign sleptons
  
  
  Integer :: a, b, i_n, i_h, i_a
  Real(dp) :: smin, smax, smin2, smax2, rj2, rk2, rm2
  Real(dp) :: r_out(3),r_outcrossed(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  
  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Slepton_ll"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSl(j)+mf_l(k)+mf_l(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n = Size(mN)
  i_h = Size(mH)
  i_a = Size(mA)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSl(i))**2
  rm2 = (mf_l(m) / mSl(i))**2
  rj2 = (mSl(j) / mSl(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  smin2 = 2._dp * Sqrt(rm2)
  smax2 = 1._dp + rm2 -  rj2 - rk2 - 2._dp * Sqrt(rj2*rk2)

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!!!!!!Slepton charge preserving !!!!!!!!!!!!!!!!!!!!!!!!!!!XXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
 r_out(1) = rj2
 r_out(2) = rk2
 r_out(3) = rm2
 
 r_outcrossed(1) = rj2
 r_outcrossed(2) = rm2
 r_outcrossed(3) = rk2


!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   

  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)
      coup(1) = C_L(k,a,i)
      coup(2) = C_R(k,a,i)
      coup(3) = Conjg(C_R(m,a,j))
      coup(4) = Conjg(C_L(m,a,j))

      coup(5) = C_L(k,b,i)
      coup(6) = C_R(k,b,i)
      coup(7) = Conjg(C_R(m,b,j))
      coup(8) = Conjg(C_L(m,b,j))
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam(1) = gam(1) + gamTemp
      Else
       gam(1) = gam(1) + 2._dp * gamTemp
      End If
      
   End Do
  End Do

If(k.Eq.m) then

!-----------------------------------------------------------    
! 2.   HH contribution  ------------------------------------
!-----------------------------------------------------------                                          

  Do a = 1, i_h                                 
   Do b = a, i_h                                         
   mass(1) = mH(a)
   mass(2) = mH(b)
   m_in    = mSl(i)
   coup(1) = C_H(a,i,j)
   coup(2) = CL_H_L(k,m,a)
   coup(3) = CL_H_R(k,m,a)
   coup(4) = C_H(b,i,j)
   coup(5) = CL_H_L(k,m,b)
   coup(6) = CL_H_R(k,m,b)
   
   If(a.Eq.b) then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    gamTemp = 2._dp * gamTemp
   End If
   
   gam(1) = gam(1) + Real(gamTemp)
   End Do
  End Do                                   


!-----------------------------------------------------------    
! 3.   AA contribution  ------------------------------------
!-----------------------------------------------------------                                          

  Do a = 1, i_a                                 
   Do b = a, i_a                                         
   mass(1) = mA(a)
   mass(2) = mA(b)
   m_in    = mSl(i)
   coup(1) = C_A(a,i,j)
   coup(2) = CL_A_L(k,m,a)
   coup(3) = CL_A_R(k,m,a)
   
   coup(4) = C_A(b,i,j)
   coup(5) = CL_A_L(k,m,b)
   coup(6) = CL_A_R(k,m,b)

   If(a.Eq.b) then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
    gamTemp = 2._dp * gamTemp
   End If

   gam(1) = gam(1) + Real(gamTemp)

   End Do
  End Do                                                    
                             
                   
!-----------------------------------------------------------    
! 4.   HA contribution  ------------------------------------
!-----------------------------------------------------------                                          

 Do a = 1, i_h                                 
   Do b = 1, i_a                                         
   mass(1) = mH(a)
   mass(2) = mA(b)
   m_in    = mSl(i)
   coup(1) = C_H(a,i,j)
   coup(2) = CL_H_L(k,m,a)
   coup(3) = CL_H_R(k,m,a)
   
   coup(4) = C_A(b,i,j)
   coup(5) = CL_A_L(k,m,b)
   coup(6) = CL_A_R(k,m,b)
   
   Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

   End Do
  End Do    
  
                            
!-----------------------------------------------------------    
! 5.   ZZ contribution  ------------------------------------
!-----------------------------------------------------------                                          

   mass(1) = mZ
   m_in    = mSl(i)
   coup(1) = C_Z(j,i,1)
   coup(2) = CL_Z_L
   coup(3) = CL_Z_R
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam(1) = gam(1) + Real(gamTemp)

!-----------------------------------------------------------    
! 6.   ZChi0 contribution  ------------------------------------
!-----------------------------------------------------------                                     

  Do a = 1,i_n
   mass(1) = mN(a)
   mass(2) = mZ

   m_in = mSl(i)
   
   coup(1) = C_L(k,a,i)
   coup(2) = C_R(k,a,i)
   coup(3) = Conjg(C_R(m,a,j))
   coup(4) = Conjg(C_L(m,a,j))  
   
   coup(5) = C_Z(j,i,1)
   coup(6) = CL_Z_L
   coup(7) = CL_Z_R


   
   Call IntegrateVF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

  End Do 
  

!-----------------------------------------------------------    
! 7.   ZH contribution  ---------------------------------
!-----------------------------------------------------------                                         
                         
  Do b = 1,i_h
   mass(1) = mZ
   mass(2) = mH(b)
   m_in = mSl(i)
   coup(1) = C_Z(j,i,1)
   coup(2) = CL_Z_L
   coup(3) = CL_Z_R
   
   coup(4) = C_H(b,i,j)
   coup(5) = CL_H_L(k,m,b)
   coup(6) = CL_H_R(k,m,b)

   Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

  End Do                                    

!-----------------------------------------------------------    
! 8.   ZA contribution  ---------------------------------
!-----------------------------------------------------------                                         
                           
  Do b = 1,i_a
   mass(1) = mZ
   mass(2) = mA(b)
   m_in = mSl(i)
   coup(1) = C_Z(i,j,1)
   coup(2) = CL_Z_L
   coup(3) = CL_Z_R
   
   coup(4) = C_A(b,i,j)
   coup(5) = CL_A_L(k,m,b)
   coup(6) = CL_A_R(k,m,b)
  
   If(mass(1).Eq.mass(2)) then
   Call IntegrateVSGoldstone(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)                                
   else
   Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
   End If
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

  End Do                                                                      


!-----------------------------------------------------------    
! 9.  HF  contribution  ---------------------------------
!-----------------------------------------------------------                                       

  Do a = 1, i_n
   Do b = 1 , i_h
   mass(1) = mN(a)
   mass(2) = mH(b)
   
   coup(1) = C_L(k,a,i)
   coup(2) = C_R(k,a,i)
   coup(3) = Conjg(C_R(m,a,j))
   coup(4) = Conjg(C_L(m,a,j))

   coup(5) = C_H(b,i,j)
   coup(6) = CL_H_L(k,m,b)
   coup(7) = CL_H_R(k,m,b)

   Call IntegrateSF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

   End Do  
  End Do

!-----------------------------------------------------------    
!10.  AF  contribution  ---------------------------------
!-----------------------------------------------------------                                       

  Do a = 1, i_n
   Do b = 1 , i_a
   mass(1) = mN(a)
   mass(2) = mA(b)
   
   coup(1) = C_L(k,a,i)
   coup(2) = C_R(k,a,i)
   coup(3) = Conjg(C_R(m,a,j))
   coup(4) = Conjg(C_L(m,a,j))

   coup(5) = C_A(b,i,j)
   coup(6) = CL_A_L(k,m,b)
   coup(7) = CL_A_R(k,m,b)

   Call IntegrateSF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam(1) = gam(1) + 2._dp*Real(gamTemp)

   End Do  
  End Do
End If  
  

 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!!!!!!Slepton charge changing !!!!!!!!!!!!!!!!!!!!!!!!!!!XXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
!-----------------------------------------------------------    
!1.  uncrossed contribution ---------------------------------
!-----------------------------------------------------------    

Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)
      coup(1) = C_L(k,a,i)
      coup(2) = C_R(k,a,i)
      coup(3) = C_L(m,a,j)
      coup(4) = C_R(m,a,j)

      coup(5) = C_L(k,b,i)
      coup(6) = C_R(k,b,i)
      coup(7) = C_L(m,b,j)
      coup(8) = C_R(m,b,j)
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam(2) = gam(2) + gamTemp
      Else
       gam(2) = gam(2) + 2._dp * gamTemp
      End If

   End Do
  End Do

!-----------------------------------------------------------    
!2. flipped contribution squared ---------------------------------
!-----------------------------------------------------------    

  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)
      coup(1) = C_L(m,a,i)
      coup(2) = C_R(m,a,i)
      coup(3) = C_L(k,a,j)
      coup(4) = C_R(k,a,j)

      coup(5) = C_L(m,b,i)
      coup(6) = C_R(m,b,i)
      coup(7) = C_L(k,b,j)
      coup(8) = C_R(k,b,j)
      Call IntegrateFFLM(mass,m_in,r_outcrossed, coup,smin2, smax2, eps, gamTemp)
      If (a.Eq.b) Then
       gam(2) = gam(2) + gamTemp
      Else
       gam(2) = gam(2) + 2._dp * gamTemp
      End If

   End Do
  End Do

!-----------------------------------------------------------    
!3.  interference contribution  ---------------------------------
!-----------------------------------------------------------    
Do a = 1, i_n                                
   Do b = 1, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)
      
      coup(1) = C_L(k,a,i)
      coup(2) = C_R(k,a,i)
      coup(3) = C_L(m,a,j)
      coup(4) = C_R(m,a,j)
      
      coup(5) = C_L(m,b,i)
      coup(6) = C_R(m,b,i)
      coup(7) = C_L(k,b,j)
      coup(8) = C_R(k,b,j)
      Call IntegrateChiChiInterference(mass,m_in,r_out,coup,smin,smax,eps,gamTemp)
   
      gam(2) = gam(2) + 2._dp*Real(gamTemp)

   End Do
  End Do

  
  If(m.Eq.k) gam(2) = 0.5_dp*gam(2)

  Iname = Iname - 1

 End Subroutine Slepton_to_Slepton_ll
 

 Subroutine Slepton_to_Slepton_nunu(i, j, k, m, mSl,m_nu, mC, mZ &
                         &, C_Lnu, C_Rnu, C_Z, Cnu_Z_L, Cnu_Z_R , eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of neutrino k
  Integer, Intent(in) :: m         ! index of neutrino m
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: m_nu(:)  ! neutrxino masses
  Real(dp), Intent(in) :: mC(:)    ! chargino masses
  Real(dp), Intent(in) :: mZ       ! suprise surprise Z-mass!!
  Complex(dp), Intent(in) :: C_Z(:,:,:) ! couplings of sleptons to the Z 
  Complex(dp), Intent(in) :: C_Lnu(:,:,:), C_Rnu(:,:,:) !LR couplings
  Real(dp), Intent(in) :: Cnu_Z_L ! couplings of leptons to the Z 
  Real(dp), Intent(in) :: Cnu_Z_R ! couplings of leptons to the Z 
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  
  
  Integer :: a, b, i_cha
  Real(dp) :: smin, smax, rk2, rm2, rj2
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Slepton_nunu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSl(j)+m_nu(k)+m_nu(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).gt.(abs(mC(1)))) then
   Iname = Iname - 1
   Return
  Else If (mSl(i).Gt.(mSl(j)+mZ)) Then
   Iname = Iname - 1
   Return
  End If

  i_cha = Size(mC)

  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (m_nu(k) / mSl(i))**2
  rm2 = (m_nu(m) / mSl(i))**2
  rj2 = (mSl(j) / mSl(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)

  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
! -----------------------------------------------------------    
! 1.   ChaCha contribution  ------------------------------------
! -----------------------------------------------------------    
  Do a = 1, i_cha                                
   Do b = a, i_cha                              
      mass(1) = mC(a)
      mass(2) = mC(b)
      m_in    = mSl(i)
      coup(1) = C_Lnu(a,k,i)
      coup(2) = C_Rnu(a,k,i)
      coup(3) = Conjg(C_Rnu(a,m,j))
      coup(4) = Conjg(C_Lnu(a,m,j))

      coup(5) = C_Lnu(b,k,i)
      coup(6) = C_Rnu(b,k,i)
      coup(7) = Conjg(C_Rnu(b,m,j))
      coup(8) = Conjg(C_Lnu(b,m,j))
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If
      
   End Do
  End Do

  
  If(k.Eq.m) Then 
!-----------------------------------------------------------    
! 2.   ZZ contribution  ------------------------------------
!-----------------------------------------------------------     

   mass(1) = mZ
   m_in    = mSl(i)
   coup(1) = C_Z(j,i,1)
   coup(2) = Cnu_Z_L
   coup(3) = Cnu_Z_R
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam = gam + Real(gamTemp)

   
!-----------------------------------------------------------    
! 2.   ZCha contribution  ----------------------------------
!-----------------------------------------------------------       

   Do a = 1,i_cha
    mass(1) = mC(a)
    mass(2) = mZ
    m_in = mSl(i)

    coup(1) = C_Lnu(a,k,i)
    coup(2) = C_Rnu(a,k,i)
    coup(3) = Conjg(C_Rnu(a,m,j))
    coup(4) = Conjg(C_Lnu(a,m,j))
      
    coup(5) = C_Z(j,i,1)
    coup(6) = Cnu_Z_L
    coup(7) = Cnu_Z_R
    Call IntegrateVF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
    gam = gam + 2._dp*Real(gamTemp)
   End Do  
 
  End If

   Iname = Iname - 1

 End Subroutine Slepton_to_Slepton_nunu

 
 Subroutine Slepton_to_Slepton_qq(i, j, k, m, mSl, mf_q, mZ, mH, mA,      &
                              &   C_Z, C_H, C_A, CU_H_L, CU_H_R, CU_A_L,  &
                              &   CU_A_R, CU_Z_L, CU_Z_R , eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 ! 17.08.2012: slight modifications by Werner Porod
 !             - take mf_q as input
 !             - generalize it from uu to qq
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of quark k
  Integer, Intent(in) :: m         ! index of quark m
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mH(:)    ! higgs masses 
  Real(dp), Intent(in) :: mf_q(:)  ! quark masses 
  Real(dp), Intent(in) :: mA(:)    ! pseudo scalar higgs masses 
  Real(dp), Intent(in) :: mZ       ! suprise surprise Z-mass!!
  Complex(dp), Intent(in) :: C_Z(:,:,:) ! couplings of sleptons to the Z 
  Complex(dp), Intent(in) :: C_H(:,:,:) ! couplings of sleptons to the scalar   
  Complex(dp), Intent(in) :: C_A(:,:,:) ! couplings of sleptons to the pseudosc
  Real(dp), Intent(in) :: CU_Z_L ! couplings of leptons to the Z 
  Real(dp), Intent(in) :: CU_Z_R ! couplings of leptons to the Z 
  Complex(dp), Intent(in) :: CU_H_L(:,:,:) ! left couplings of leptons to sc.
  Complex(dp), Intent(in) :: CU_A_L(:,:,:) ! left couplings of leptons to psc.
  Complex(dp), Intent(in) :: CU_H_R(:,:,:) ! right couplings of leptons to sc,     
  Complex(dp), Intent(in) :: CU_A_R(:,:,:) ! right couplings of leptons to psc.

  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  

  Integer :: a, b, i_h, i_a
  Real(dp) :: smin, smax
  Real(dp) :: r_out(3), rj2, rk2, rm2
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Slepton_qq"
  
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSl(j)+mf_q(k)+mf_q(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).gt.(mSl(j)+Min(mZ,mH(1))) ) Then
   Iname = Iname - 1
   Return
  Else If (k.Ne.m) Then
   Write(ErrCan,*) &
     & "Warning, the case m.ne.k is not yet included in "//NameOfUnit(Iname)
   Iname = Iname - 1
   Return
  End If

  i_h = Size(mH)
  i_a = Size(mA)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_q(k) / mSl(i))**2
  rm2 = (mf_q(m) / mSl(i))**2
  rj2 = (mSl(j) / mSl(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
 
!-----------------------------------------------------------    
! 1.   ZZ contribution  ------------------------------------
!-----------------------------------------------------------   
   mass(1) = mZ
   m_in    = mSl(i)
   coup(1) = C_Z(j,i,1)
   coup(2) = CU_Z_L
   coup(3) = CU_Z_R
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   gam = gam + Real(gamTemp)

!-----------------------------------------------------------    
! 2.   HH contribution  ------------------------------------
!-----------------------------------------------------------    
  Do a = 1, i_h                                 
   Do b = a, i_h               
   mass(1) = mH(a)
   mass(2) = mH(b)
   m_in    = mSl(i)
   
   coup(1) = C_H(a,i,j)
   coup(2) = CU_H_L(k,m,a)
   coup(3) = CU_H_R(k,m,a)
   
   coup(4) = C_H(b,i,j)
   coup(5) = CU_H_L(k,m,b)
   coup(6) = CU_H_R(k,m,b)
   
   If(a.Eq.b) Then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   Else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    gamTemp = 2._dp * gamTemp
   End If

   gam = gam + Real(gamTemp)

   End Do
  End Do                     

!-----------------------------------------------------------    
! 2.   AA contribution  ------------------------------------
!-----------------------------------------------------------    
  Do a = 1, i_a                                 
   Do b = a, i_a                        
   mass(1) = mA(a)
   mass(2) = mA(b)
   m_in    = mSl(i)
   coup(1) = C_A(a,i,j)
   coup(2) = CU_A_L(k,m,a)
   coup(3) = CU_A_R(k,m,a)
   
   coup(4) = C_A(b,i,j)
   coup(5) = CU_A_L(k,m,b)
   coup(6) = CU_A_R(k,m,b)
   
   If(a.Eq.b) Then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   Else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    gamTemp = 2._dp * gamTemp
   End If
   
   gam = gam + Real(gamTemp)

   End Do
  End Do          

!-----------------------------------------------------------    
! 3.   HA contribution  ------------------------------------
!-----------------------------------------------------------                                          

  Do a = 1, i_h                                 
   Do b = 1, i_a          
   mass(1) = mH(a)
   mass(2) = mA(b)
   m_in    = mSl(i)

   coup(1) = C_H(a,i,j)
   coup(2) = CU_H_L(k,m,a)
   coup(3) = CU_H_R(k,m,a)
   
   coup(4) = C_A(b,i,j)
   coup(5) = CU_A_L(k,m,b)
   coup(6) = CU_A_R(k,m,b)
   
   Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam = gam + 2._dp*Real(gamTemp)

   End Do
  End Do     

!-----------------------------------------------------------    
! 4.   ZH contribution  ---------------------------------
!-----------------------------------------------------------                                             
  Do b = 1,i_h
   mass(1) = mZ
   mass(2) = mH(b)
   m_in = mSl(i)
   
   coup(1) = C_Z(i,j,1)
   coup(2) = CU_Z_L
   coup(3) = CU_Z_R
   
   coup(4) = C_H(b,i,j)
   coup(5) = CU_H_L(k,m,b)
   coup(6) = CU_H_R(k,m,b)
   
   Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)  
   
   gam = gam + 2._dp*Real(gamTemp)

  End Do   

!-----------------------------------------------------------    
! 5.   ZA contribution  ---------------------------------
!-----------------------------------------------------------                                         
  Do b = 1,i_a
   mass(1) = mZ
   mass(2) = mA(b)
   m_in = mSl(i)
   coup(1) = C_Z(i,j,1)
   coup(2) = CU_Z_L
   coup(3) = CU_Z_R
   
   coup(4) = C_A(b,i,j)
   coup(5) = CU_A_L(k,m,b)
   coup(6) = CU_A_R(k,m,b)
   
   If(mass(1).Eq.mass(2)) Then
   Call IntegrateVSGoldstone(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)                                
   Else
   Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
   End If
   
   gam = gam + 2._dp*Real(gamTemp)

  End Do      

  gam = 3._dp * gam ! color factor
 
  Iname = Iname - 1
 
  End Subroutine Slepton_to_Slepton_qq
  

 Subroutine Slepton_to_Antisneutrino_lnu(i, j, k, m, mN, mSl, mSn, mC, mf_l  &
    & , mf_nu, C_NSlL_L, C_NSlL_R , C_NSnNu_L, C_NSnNu_R, C_CSnL_L, C_CSnL_R &
    & , C_CNuSl_L, C_CNuSl_R,  eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
  Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state antisneutrino
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of neutrino m
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mC(:)    ! chargino masses 
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Real(dp), Intent(in) :: mf_nu(:) ! neutrino masses
  Complex(dp), Intent(in) :: C_NSlL_L(:,:,:) ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_NSlL_R(:,:,:) !    to neutralinos 
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Complex(dp), Intent(in) :: C_CSnL_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_CSnL_R(:,:,:) !    to chargninos
  Complex(dp), Intent(in) :: C_CNuSl_L(:,:,:) ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_CNuSl_R(:,:,:) !    to chargninos
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam     ! partial width

  Integer :: a, b, i_n, i_cha
  Real(dp) :: smin, smax,smin2,smax2, rk2, rm2, rj2
  Real(dp) :: r_out(3),r_outcrossed(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp
  
  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Antisneutrino_lnu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSn(j)+ mf_l(k)+ mf_nu(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).Gt.(Abs(mN(1))+ mf_l(k)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n   = Size(mN)
  i_cha = Size(mC)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSl(i))**2
  rm2 = (mf_nu(m) / mSl(i))**2
  rj2 = (mSn(j) / mSl(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
  
  smin2 = 2._dp * Sqrt(rm2)
  smax2 = 1._dp + rm2 -  rj2 - rk2 - 2._dp * Sqrt(rj2*rk2)
 
  r_outcrossed(1) = rj2
  r_outcrossed(2) = rm2
  r_outcrossed(3) = rk2

!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)

      coup(1) = C_NSlL_L(k,a,i)
      coup(2) = C_NSlL_R(k,a,i)
      coup(3) = C_NSnNu_L(m,a,j)
      coup(4) = C_NSnNu_R(m,a,j)
      
      coup(5) = C_NSlL_L(k,b,i)
      coup(6) = C_NSlL_R(k,b,i)
      coup(7) = C_NSnNu_L(m,b,j)
      coup(8) = C_NSnNu_R(m,b,j)
      
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do



!-----------------------------------------------------------    
! 2.  ChaCha contribution  ---------------------------------
!-----------------------------------------------------------     

  Do a = 1, i_cha                                
   Do b = a, i_cha                              
      mass(1) = mC(a)
      mass(2) = mC(b)
      m_in    = mSl(i)
      coup(1) = C_CNuSl_L(a,m,i)
      coup(2) = C_CNuSl_R(a,m,i)
      coup(3) = C_CSnL_L(a,k,j)
      coup(4) = C_CSnL_R(a,k,j)

      coup(5) = C_CNuSl_L(b,m,i)
      coup(6) = C_CNuSl_R(b,m,i)
      coup(7) = C_CSnL_L(b,k,j)
      coup(8) = C_CSnL_R(b,k,j)
      Call IntegrateFFLM(mass,m_in,r_outcrossed, coup,smin2, smax2, eps, gamTemp)
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do


!-----------------------------------------------------------    
! 2.  ChaChi contribution  ---------------------------------
!-----------------------------------------------------------  
  Do a = 1, i_n                                
   Do b = 1, i_cha                              
      mass(1) = mN(a)
      mass(2) = mC(b)
      m_in    = mSl(i)
      
      coup(1) = C_NSlL_L(k,a,i)
      coup(2) = C_NSlL_R(k,a,i)
      coup(3) = C_NSnNu_L(m,a,j)
      coup(4) = C_NSnNu_R(m,a,j)

      coup(5) = C_CNuSl_L(b,m,i)
      coup(6) = C_CNuSl_R(b,m,i)
      coup(7) = C_CSnL_L(b,k,j)
      coup(8) = C_CSnL_R(b,k,j)
      Call IntegrateChiChiInterference(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      
      gam = gam + 2._dp*Real(gamTemp)

   End Do
  End Do 

  Iname = Iname - 1

 End Subroutine Slepton_to_Antisneutrino_lnu  
 

 
 Subroutine Slepton_to_Sneutrino_lnu(i,j,k,m, mN, mSl, mSn, mW, mHpm, mf_l     &
       & , mf_nu, C_NSlL_L, C_NSlL_R, C_NSnNu_L, C_NSnNu_R, C_WSlSn, C_HpmSlSn &
       & , C_WLNu, C_HpmLNu_L, C_HpmLNu_R, eps, gam                     )
  !----------------------------------------------------------------------------
  ! written by Lukas Mitzka
  !----------------------------------------------------------------------------
  Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of lepton m
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mHpm(:)  ! charged higgs masses 
  Real(dp), Intent(in) :: mW       ! W mass
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Real(dp), Intent(in) :: mf_nu(:) ! neutrino masses
  Complex(dp), Intent(in) :: C_NSlL_L(:,:,:)  ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_NSlL_R(:,:,:)  !    to neutralinos 
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Complex(dp), Intent(in) :: C_WSlSn(:,:,:)     ! coupling of sfermions to W
  Complex(dp), Intent(in) :: C_HpmSlSn(:,:,:) ! coupling of sfermions to Hpm
  Complex(dp), Intent(in) :: C_WLNu     ! coupling of leptons to W
  Complex(dp), Intent(in) :: C_HpmLNu_L(:,:,:)  ! left coupling of leptons to Hpm
  Complex(dp), Intent(in) :: C_HpmLNu_R(:,:,:)  ! left coupling of leptons to Hpm
  
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  ! partial widths 
 
  Integer :: a, b, i_n, i_hpm
  Real(dp) :: smin, smax, rk2, rm2, rj2
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Sneutrino_lnu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSn(j)+mf_l(k)+mf_nu(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).Gt.(Abs(mN(1))+mf_l(k)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n = Size(mN)
  i_hpm = Size(mHpm)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSl(i))**2
  rm2 = (mf_nu(m) / mSl(i))**2
  rj2 = (mSn(j) / mSl(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
   
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
                            
!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   

  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSl(i)
      coup(1) = C_NSlL_L(k,a,i)
      coup(2) = C_NSlL_R(k,a,i)
      coup(3) = Conjg(C_NSnNu_R(m,a,j))
      coup(4) = Conjg(C_NSnNu_L(m,a,j))

      coup(5) = C_NSlL_L(k,b,i)
      coup(6) = C_NSlL_R(k,b,i)
      coup(7) = Conjg(C_NSnNu_R(m,b,j))
      coup(8) = Conjg(C_NSnNu_L(m,b,j))
      
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do
  
!-----------------------------------------------------------    
! 2.   HpmHpm contribution  --------------------------------
!-----------------------------------------------------------
  Do a = 1, i_hpm                                 
   Do b = a, i_hpm                                         
   mass(1) = mHpm(a)
   mass(2) = mHpm(b)
   m_in    = mSl(i)
   
   coup(1) = C_HpmSlSn(a,i,j)
   coup(2) = C_HpmLNu_L(a,k,m)
   coup(3) = C_HpmLNu_R(a,k,m)

   coup(4) = C_HpmSlSn(b,i,j)
   coup(5) = C_HpmLNu_L(b,k,m)
   coup(6) = C_HpmLNu_R(b,k,m)

   If(a.Eq.b) Then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   Else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    gamTemp = 2._dp * gamTemp
   End If
   
   gam = gam + Real(gamTemp)

   End Do
  End Do       

!-----------------------------------------------------------    
! 3.    WW    contribution  --------------------------------
!-----------------------------------------------------------            

   mass(1) = mW
   m_in    = mSl(i)
   coup(1) = Conjg(C_WSlSn(i,j,1))
   coup(2) = C_WLNu
   coup(3) = 0._dp
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)

   gam = gam + Real(gamTemp)
  
!-----------------------------------------------------------    
! 4.  WHpm   contribution  ---------------------------------
!-----------------------------------------------------------                                         
                              
  Do b = 1,i_hpm
   mass(1) = mW
   mass(2) = mHpm(b)
   m_in = mSl(i)
   
   coup(1) = Conjg(C_WSlSn(i,j,1))
   coup(2) = C_WLNu
   coup(3) = 0._dp
   
   coup(4) = C_HpmSlSn(b,i,j)
   coup(5) = C_HpmLNu_L(b,k,m)
   coup(6) = C_HpmLNu_R(b,k,m)
   
   If(mass(1).Eq.mass(2)) Then
    Call IntegrateVSGoldstone(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)                                
   Else
    Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)  
   End If

   gam = gam + 2._dp*Real(gamTemp)

  End Do                   

!-----------------------------------------------------------    
! 5.   WChi0 contribution  ---------------------------------
!-----------------------------------------------------------                                    
  Do a = 1,i_n
   mass(1) = mN(a)
   mass(2) = mW

   m_in = mSl(i)

   coup(1) = C_NSlL_L(k,a,i)
   coup(2) = C_NSlL_R(k,a,i)
   coup(3) = Conjg(C_NSnNu_R(m,a,j))
   coup(4) = Conjg(C_NSnNu_L(m,a,j))
   
   coup(5) = Conjg(C_WSlSn(i,j,1))
   coup(6) = C_WLNu
   coup(7) = 0._dp
   
   Call IntegrateVF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam = gam + 2._dp*Real(gamTemp)
  End Do   
  
!-----------------------------------------------------------    
! 6.  HpmChi0  contribution  --------------------------------
!-----------------------------------------------------------                                       
  Do a = 1, i_n
   Do b = 1 , i_hpm
    mass(1) = mN(a)
    mass(2) = mHpm(b)
    m_in    = mSl(i)
   
    coup(1) = C_NSlL_L(k,a,i)
    coup(2) = C_NSlL_R(k,a,i)
    coup(3) = Conjg(C_NSnNu_R(m,a,j))
    coup(4) = Conjg(C_NSnNu_L(m,a,j))
   
    coup(5) = C_HpmSlSn(b,i,j)
    coup(6) = C_HpmLNu_L(b,k,m)
    coup(7) = C_HpmLNu_R(b,k,m)

    Call IntegrateSF(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)

    gam = gam + 2._dp*Real(gamTemp)
   End Do  
  End Do 

  Iname = Iname - 1

 End Subroutine Slepton_to_Sneutrino_lnu
 
 
 Subroutine Slepton_to_Sneutrino_qq(i,j,k,m, mSl, mSn, mHpm, mW, mf_d, mf_u &
      & , C_HpmUD_L, C_HpmUD_R, C_WUD, C_HpmSlSn, C_WSlSn, eps, gam )
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
  Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state sneutrino
  Integer, Intent(in) :: k         ! index of down type quark k
  Integer, Intent(in) :: m         ! index of up type m
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mHpm(:)  ! higgs masses 
  Real(dp), Intent(in) :: mW       ! suprise surprise W-mass!!
  Real(dp), Intent(in) :: mf_d(:)  ! d-quark masses
  Real(dp), Intent(in) :: mf_u(:)  ! u-quark masses
  Complex(dp), Intent(in) :: C_HpmUD_L(:,:,:)  ! left coupling of quarks to Hpm
  Complex(dp), Intent(in) :: C_HpmUD_R(:,:,:)  ! right coupling of quarks to Hpm
  Complex(dp), Intent(in) :: C_WUD(:,:)        ! coupling of quarks to W
  Complex(dp), Intent(in) :: C_HpmSlSn(:,:,:) ! coupling of sfermions to Hpm\
  Complex(dp), Intent(in) :: C_WSlSn(:,:,:)     ! coupling of fermions to W
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam     ! partial width

  Integer :: a, b, i_hpm
  Real(dp) :: smin, smax, rk2, rm2, rj2
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Sneutrino_qq"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSn(j)+mf_d(k)+mf_u(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).gt.(mSn(j)+mW) ) Then
   Iname = Iname - 1
   Return
  End If

  i_hpm = Size(mHpm)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_d(k) / mSl(i))**2
  rm2 = (mf_u(m) / mSl(i))**2
  rj2 = ( mSn(j) / mSl(i))**2
  
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)

  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2

!-----------------------------------------------------------    
! 1.   HpmHpm contribution  --------------------------------
!-----------------------------------------------------------     
  Do a = 1, i_hpm                                 
   Do b = a, i_hpm                                         
   mass(1) = mHpm(a)
   mass(2) = mHpm(b)
   m_in    = mSl(i)
   
   coup(1) = C_HpmSlSn(a,i,j)
   coup(2) = C_HpmUD_L(a,k,m)
   coup(3) = C_HpmUD_R(a,k,m)
   
   coup(4) = C_HpmSlSn(b,i,j)
   coup(5) = C_HpmUD_L(b,k,m)
   coup(6) = C_HpmUD_R(b,k,m)
   
   If(a.Eq.b) Then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   Else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
    gamTemp = 2._dp * gamTemp
   End If
   
   gam = gam +3._dp * Real(gamTemp)

   End Do
  End Do                                
 
!-----------------------------------------------------------    
! 2.   WW contribution  ------------------------------------
!-----------------------------------------------------------   
   mass(1) = mW
   m_in    = mSl(i)
   coup(1) = Conjg(C_WSlSn(i,j,1))
   coup(2) = C_WUD(m,k)
   coup(3) = 0._dp
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam = gam + 3._dp *Real(gamTemp)

!-----------------------------------------------------------    
! 3. HpmW  contribution  -----------------------
!-----------------------------------------------------------       
  Do b = 1,i_hpm
   mass(1) = mW
   mass(2) = mHpm(b)
   m_in = mSl(i)
   
   coup(1) = Conjg(C_WSlSn(i,j,1))
   coup(3) = C_WUD(m,k)
   coup(2) = 0._dp
   
   coup(4) = Conjg(C_HpmSlSn(b,i,j))
   coup(5) = Conjg(C_HpmUD_R(b,m,k))
   coup(6) = Conjg(C_HpmUD_L(b,m,k))
   
   If(mass(1).Eq.mass(2)) Then
    Call IntegrateVSGoldstone(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)                                
   Else
    Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
   End If
   
   gam = gam + 3._dp *2._dp*Real(gamTemp)
   
  End Do               

  Iname = Iname - 1

 End Subroutine Slepton_to_Sneutrino_qq
  

 Subroutine Sneutrino_3body(n_in, n_l, id_l, n_nu, id_nu, n_d, id_d, n_u, id_u &
   & , n_c, n_n, n_W, id_W, n_Sl, n_Snu, n_Spm, Sneut, mW, mf_l, mf_u, mf_d    &
   & , Chi0, ChiPm, Spm, Slept, c_LNSl_L, c_LNSl_R, c_SlSnW, c_SmpSlSn         &
   & , c_CNuSl_L, c_CNuSl_R, c_LNuW, c_DUW, c_CLSn_L, c_CLSn_R                 &
   & , c_NuNSn_L, c_NuNSn_R, c_SmpLNu_L, c_SmpLNu_R, c_SmpDU_L, c_SmpDU_R      &
   & , GenerationMixing, epsI)
 !-----------------------------------------------------------------------------
 ! calculates the various slepton 3-body decays using the routines provided
 ! by Lukas Mitzka
 ! written by Werner Porod, 17.08.2012
 !-----------------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_in, n_nu, n_l, n_d, n_u, n_W, n_snu  &
       & , n_sl, n_n, n_c, n_Spm
  Integer, Intent(in) :: id_nu(:), id_l(:), id_d(:), id_u(:), id_W(:)
  Real(dp), Intent(in) :: mf_l(3), mf_u(3), mf_d(3), mW, epsI
  Complex(dp), Intent(in) :: c_LNSl_L(:,:,:), c_LNSl_R(:,:,:), c_LNuW(:,:)    &
     & , c_CNuSl_L(:,:,:), c_CNuSl_R(:,:,:), c_SlSnW(:,:,:), c_SmpSlSn(:,:,:) &
     & , c_SmpLNu_L(:,:,:), c_SmpLNu_R(:,:,:), c_DUW(:,:), c_NuNSn_L(:,:,:)   &
     & , c_NuNSn_R(:,:,:), c_CLSn_L(:,:,:), c_CLSn_R(:,:,:), c_SmpDU_L(:,:,:) &
     & , c_SmpDU_R(:,:,:)
  Logical, Intent(in) :: GenerationMixing
  Type(particle2), Intent(in) :: Spm(:)
  Type(particle23), Intent(in) :: Slept(:), ChiPm(:), Chi0(:)
  Type(particle23), Intent(inout) :: Sneut(:)

  Integer :: i_start, i_end, i_run, i2, i3, i4, i_c, i2_max
  Real(dp) :: mN(n_n), mSlept(n_sl), mSneut(n_snu), mC(n_c), gam    &
     & , mSpm(n_Spm), m_nu(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Sneutrino_3body'

  If (n_in.Lt.0) Then
   i_start = 1
   i_end = n_snu

  Else If ( (n_in.Ge.1).And.(n_in.Le.n_snu) ) Then 
   i_start = n_in 
   i_end = n_in

  Else
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write (ErrCan,*) 'Value of n_in out of range, (n_in,n_sl) = ',n_in,n_sl
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   Iname = Iname - 1
   Return
  End If

  !------------------------
  ! initialisation
  !------------------------
  m_nu = 0._dp

  mC = ChiPm%m
  mN = Chi0%m
  mSlept = Slept%m
  mSneut = Sneut%m
  mSpm = Spm%m

  Do i_run = i_start, i_end
   Sneut(i_run)%gi3 = 0._dp
   Sneut(i_run)%bi3 = 0._dp
   i_c = 1
   If (GenerationMixing) Then
    i2_max = i_run - 1
   Else
    i2_max = n_snu
   End If
   !---------------------------------------------------------------------------
   ! sneutrino l- l+
   !---------------------------------------------------------------------------
    Do i2=1,i2_max
     Do i3=1,3
      Do i4=1,3
       Call Sneutrino_to_Sneutrino_ll(i_run, i2, i3, i4, mSneut, mC, mf_l   &
                &  , c_CLSn_L, c_CLSn_R, epsI, gam )

       Sneut(i_run)%gi3(i_c)   = gam
       Sneut(i_run)%id3(i_c,1) = Sneut(i2)%id 
       Sneut(i_run)%id3(i_c,2) = id_l(i3)
       Sneut(i_run)%id3(i_c,3) = id_l(i4) + 1
       i_c = i_c + 1
      End Do
     End Do
    End Do
   !---------------------------------------------------------------------------
   ! sneutrino nu bar(nu)
   !---------------------------------------------------------------------------
    Do i2=1,i2_max
     Do i3=1,3
      Do i4=1,3
       Call Sneutrino_to_Sneutrino_nunu(i_run, i2, i3, i4, mN, mSneut, m_nu &
                            &  , c_NuNSn_L, c_NuNSn_R, epsI, gam )

       Sneut(i_run)%gi3(i_c)   = gam
       Sneut(i_run)%id3(i_c,1) = Sneut(i2)%id 
       Sneut(i_run)%id3(i_c,2) = id_nu(i3)
       Sneut(i_run)%id3(i_c,3) = id_nu(i4) + 1
       i_c = i_c + 1
      End Do
     End Do
    End Do
   !---------------------------------------------------------------------------
   ! sneutrino^* nu nu
   !---------------------------------------------------------------------------
    Do i2=1,i2_max
     Do i3=1,3
      Do i4=i3,3
       Call Sneutrino_to_AntiSneutrino_nunu(i_run, i2, i3, i4, mN, mSneut  &
                            & , m_nu, c_NuNSn_L, c_NuNSn_R, epsI, gam )

       Sneut(i_run)%gi3(i_c)   = gam
       Sneut(i_run)%id3(i_c,1) = Sneut(i2)%id +1
       Sneut(i_run)%id3(i_c,2) = id_nu(i3)
       Sneut(i_run)%id3(i_c,3) = id_nu(i4)
       i_c = i_c + 1
      End Do
     End Do
    End Do
   !---------------------------------------------------------------------------
   ! slepton^- l+ nu
   !---------------------------------------------------------------------------
    Do i2=1,n_sl
     Do i3=1,3
      Do i4=1,3
       Call Sneutrino_to_Slepton_lnu(i_run, i2, i3, i4, mN, mSlept, mSneut, mW &
          & , mSpm, mf_l, m_nu, c_LNSl_L, c_LNSl_R, c_NuNSn_L, c_NuNSn_R       &
          & , c_SlSnW, c_SmpSlSn, c_LNuW(i3,i4), c_SmpLNu_L, c_SmpLNu_R, epsI  &
          & , gam )

       Sneut(i_run)%gi3(i_c)   = gam
       Sneut(i_run)%id3(i_c,1) = Slept(i2)%id +1
       Sneut(i_run)%id3(i_c,2) = id_l(i3)
       Sneut(i_run)%id3(i_c,3) = id_nu(i4)
       i_c = i_c + 1
      End Do
     End Do
    End Do

   !---------------------------------------------------------------------------
   ! slepton^+ l- nu
   !---------------------------------------------------------------------------
    Do i2=1,n_sl
     Do i3=1,3
      Do i4=1,3
       Call Sneutrino_to_Antislepton_lnu(i_run, i2, i3, i4, mN, mSlept, mSneut &
                  & , mC, mf_l, m_nu, c_LNSl_L, c_LNSl_R, c_NuNSn_L, c_NuNSn_R &
                  & , c_CLSn_L, c_CLSn_R, c_CNuSl_L, c_CNuSl_R, epsI, gam )

       Sneut(i_run)%gi3(i_c)   = gam
       Sneut(i_run)%id3(i_c,1) = Slept(i2)%id +1
       Sneut(i_run)%id3(i_c,2) = id_l(i3)
       Sneut(i_run)%id3(i_c,3) = id_nu(i4)
       i_c = i_c + 1
      End Do
     End Do
    End Do
   !---------------------------------------------------------------------------
   ! slepton q bar(q)'
   !---------------------------------------------------------------------------
    Do i2=1,n_sl
     Do i3=1,3
      If (GenerationMixing) Then
       Do i4=1,3
        Call Sneutrino_to_Slepton_qq(i_run, i2, i3, i4, mSlept, mSneut, mSpm  &
          & , mW, mf_d, mf_u, c_SmpDU_L, c_SmpDU_R, c_DUW, c_SmpSlSn, c_SlSnW &
          & , epsI, gam)
        Sneut(i_run)%gi3(i_c) = gam
        Sneut(i_run)%id3(i_c,1) = Slept(i2)%id
        Sneut(i_run)%id3(i_c,2) = id_u(i3) + 1
        Sneut(i_run)%id3(i_c,3) = id_d(i4)
        i_c = i_c + 1
       End Do

      Else ! .not.GenerationMixing

       Call Sneutrino_to_Slepton_qq(i_run, i2, i3, i3, mSlept, mSneut, mSpm   &
          & , mW, mf_d, mf_u, c_SmpDU_L, c_SmpDU_R, c_DUW, c_SmpSlSn, c_SlSnW &
          & , epsI, gam)
       Sneut(i_run)%gi3(i_c) = gam
       Sneut(i_run)%id3(i_c,1) = Slept(i2)%id
       Sneut(i_run)%id3(i_c,2) = id_u(i3) + 1
       Sneut(i_run)%id3(i_c,3) = id_d(i3)
        i_c = i_c + 1

      End If ! GenerationMixing
     End Do
    End Do 
  
   Sneut(i_run)%g = Sum(Sneut(i_run)%gi2) + Sum(Sneut(i_run)%gi3)
   If (Sneut(i_run)%g.Ne.0._dp) Then
    Sneut(i_run)%bi2 = Sneut(i_run)%gi2 / Sneut(i_run)%g
    Sneut(i_run)%bi3 = Sneut(i_run)%gi3 / Sneut(i_run)%g
   End If
   
  End Do ! i_run

  Iname = Iname - 1

 End Subroutine Sneutrino_3body

 Subroutine Sneutrino_to_Sneutrino_ll(i, j, k, m, mSn,  mC, mf_l, C_CSnL_L &
                                   & , C_CSnL_R, eps,gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of lepton m
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Real(dp), Intent(in) :: mC(:)    ! chargino masses
  Complex(dp), Intent(in) :: C_CSnL_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_CSnL_R(:,:,:) !    to chargninos
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  ! partial width

  Integer :: a, b, i_cha
  Real(dp) :: smin, smax, rj2, rk2, rm2 
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_Sneutrino_ll"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.(mSn(j)+mf_l(k)+mf_l(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).gt.abs(mC(1))) then
   Iname = Iname - 1
   Return
  End If

  i_cha = Size(mC)

  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSn(i))**2
  rm2 = (mf_l(m) / mSn(i))**2
  rj2 = (mSn(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2

! -----------------------------------------------------------    
! 1.   ChaCha contribution  ---------------------------------
! -----------------------------------------------------------    

  Do a = 1, i_cha                                
   Do b = a, i_cha                              
      mass(1) = mC(a)
      mass(2) = mC(b)
      m_in    = mSn(i)
      coup(1) = C_CSnL_L(a,k,i)
      coup(2) = C_CSnL_R(a,k,i)
      coup(3) = Conjg(C_CSnL_R(a,m,j))
      coup(4) = Conjg(C_CSnL_L(a,m,j))

      coup(5) = C_CSnL_L(b,k,i)
      coup(6) = C_CSnL_R(b,k,i)
      coup(7) = Conjg(C_CSnL_R(b,m,j))
      coup(8) = Conjg(C_CSnL_L(b,m,j))
      
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do

  Iname = Iname - 1

 End Subroutine Sneutrino_to_Sneutrino_ll
 

 Subroutine Sneutrino_to_Sneutrino_nunu(i, j, k, m,mN, mSn, mf_nu    &
                          & , C_NSnNu_L, C_NSnNu_R,  eps, gam )
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka 
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of lepton m
  Real(dp), Intent(in) :: mSn(:)   ! slepton masses
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mf_nu(:) ! neutrino masses
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  ! partial widths

  Integer :: a, b, i_n
  Real(dp) :: smin, smax, rj2, rk2, rm2
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_Sneutrino_nunu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.mSn(j) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n = Size(mN)

  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_nu(k) / mSn(i) )**2
  rm2 = (mf_nu(m) / mSn(i) )**2
  rj2 = (mSn(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
 
!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)
      coup(1) = C_NSnNu_L(k,a,i)
      coup(2) = C_NSnNu_R(k,a,i)
      coup(3) = Conjg(C_NSnNu_R(m,a,j))
      coup(4) = Conjg(C_NSnNu_L(m,a,j))

      coup(5) = C_NSnNu_L(k,b,i)
      coup(6) = C_NSnNu_R(k,b,i)
      coup(7) = Conjg(C_NSnNu_R(m,b,j))
      coup(8) = Conjg(C_NSnNu_L(m,b,j))
      
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)

      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do
 
  Iname = Iname - 1

 End Subroutine Sneutrino_to_Sneutrino_nunu 
 
 Subroutine Sneutrino_to_AntiSneutrino_nunu(i, j, k, m, mN, mSn, m_nu         &
                                      &     , C_NSnNu_L, C_NSnNu_R, eps, gam) 
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of lepton m
  Real(dp), Intent(in) :: mSn(:)   ! slepton masses
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: m_nu(:)  ! neutrino masses
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  ! partial width

  Integer :: a, b, i_n
  Real(dp) :: r_out(3),smin, smax, rj2, rk2, rm2
  Real(dp) :: r_outcrossed(3), smin2, smax2
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_AntiSneutrino_nunu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.(mSn(j)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n = Size(mN)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (m_nu(k) / mSn(i) )**2
  rm2 = (m_nu(m) / mSn(i) )**2
  rj2 = (mSn(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  smin2 = 2._dp * Sqrt(rm2)
  smax2 = 1._dp + rm2 -  rj2 - rk2 - 2._dp * Sqrt(rj2*rk2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2

  r_outcrossed(1) = rj2
  r_outcrossed(2) = rm2
  r_outcrossed(3) = rk2
 
!-----------------------------------------------------------    
!1.  11  contribution  ---------------------------------
!-----------------------------------------------------------  
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)
      coup(1) = C_NSnNu_L(k,a,i)
      coup(2) = C_NSnNu_R(k,a,i)
      coup(3) = C_NSnNu_L(m,a,j)
      coup(4) = C_NSnNu_R(m,a,j)

      coup(5) = C_NSnNu_L(k,b,i)
      coup(6) = C_NSnNu_R(k,b,i)
      coup(7) = C_NSnNu_L(m,b,j)
      coup(8) = C_NSnNu_R(m,b,j)
      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)

      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do

!-----------------------------------------------------------
!2.  22  contribution  ---------------------------------
!-----------------------------------------------------------
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)
      coup(1) = C_NSnNu_L(m,a,i)
      coup(2) = C_NSnNu_R(m,a,i)
      coup(3) = C_NSnNu_L(k,a,j)
      coup(4) = C_NSnNu_R(k,a,j)

      coup(5) = C_NSnNu_L(m,b,i)
      coup(6) = C_NSnNu_R(m,b,i)
      coup(7) = C_NSnNu_L(k,b,j)
      coup(8) = C_NSnNu_R(k,b,j)
      Call IntegrateFFLM(mass,m_in,r_outcrossed,coup,smin2, smax2, eps, gamTemp)

      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do

!-----------------------------------------------------------    
!2.  21  contribution  ---------------------------------
!-----------------------------------------------------------    
  Do a = 1, i_n                                
   Do b = 1, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)
      
      coup(1) = C_NSnNu_L(k,a,i)
      coup(2) = C_NSnNu_R(k,a,i)
      coup(3) = C_NSnNu_L(m,a,j)
      coup(4) = C_NSnNu_R(m,a,j)
      
      coup(5) = C_NSnNu_L(m,b,i)
      coup(6) = C_NSnNu_R(m,b,i)
      coup(7) = C_NSnNu_L(k,b,j)
      coup(8) = C_NSnNu_R(k,b,j)
      
      Call IntegrateChiChiInterference(mass,m_in,r_out,coup,smin,smax,eps &
        & ,gamTemp)
      gam = gam + 2._dp*Real(gamTemp)

   End Do
  End Do
  
  If(m.Eq.k) gam = 0.5_dp*gam

  Iname = Iname - 1
 
 
 End Subroutine Sneutrino_to_AntiSneutrino_nunu
 

 Subroutine Sneutrino_to_Slepton_lnu(i,j,k,m, mN, mSl, mSn, mW, mHpm, mf_l   &
                & , mf_nu, C_NSlL_L, C_NSlL_R, C_NSnNu_L, C_NSnNu_R, C_WSlSn &
                & , C_HpmSlSn, C_WLNu, C_HpmLNu_L, C_HpmLNu_R, eps, gam )
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
  Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of antineutrino m
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mHpm(:)  ! charged higgs masses 
  Real(dp), Intent(in) :: mW       ! W mass
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Real(dp), Intent(in) :: mf_nu(:) ! neutrino masses
  Complex(dp), Intent(in) :: C_NSlL_L(:,:,:)  ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_NSlL_R(:,:,:)  !    to neutralinos 
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Complex(dp), Intent(in) :: C_WSlSn(:,:,:)     ! coupling of sfermions to W
  Complex(dp), Intent(in) :: C_HpmSlSn(:,:,:) ! coupling of sfermions to Hpm
  Complex(dp), Intent(in) :: C_WLNu      ! coupling of leptons to W
  Complex(dp), Intent(in) :: C_HpmLNu_L(:,:,:)  ! left coupling of leptons to Hpm
  Complex(dp), Intent(in) :: C_HpmLNu_R(:,:,:)  ! left coupling of leptons to Hpm
  
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam  ! partial widths 
 
  Integer :: a, b, i_n, i_hpm
  Real(dp) :: smin, smax,smin2,smax2, rj2, rk2, rm2
  Real(dp) :: r_out(3),r_outcrossed(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_Slepton_lnu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.(mSl(j)+mf_l(k)+mf_nu(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n = Size(mN)
  i_hpm = Size(mHpm)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSn(i))**2
  rm2 = 0._dp
  rj2 = (mSl(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  smin2 = 2._dp * Sqrt(rm2)
  smax2 = 1._dp + rm2 -  rj2 - rk2 - 2._dp * Sqrt(rj2*rk2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2

  r_outcrossed(1) = rj2
  r_outcrossed(2) = rm2
  r_outcrossed(3) = rk2

!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)

      coup(1) = C_NSnNu_L(m,a,i)
      coup(2) = C_NSnNu_R(m,a,i)
      coup(3) = Conjg(C_NSlL_R(k,a,j))
      coup(4) = Conjg(C_NSlL_L(k,a,j))
      
      coup(5) = C_NSnNu_L(m,b,i)
      coup(6) = C_NSnNu_R(m,b,i)
      coup(7) = Conjg(C_NSlL_R(k,b,j))
      coup(8) = Conjg(C_NSlL_L(k,b,j))
    
      Call IntegrateFFLM(mass,m_in,r_outcrossed, coup,smin2, smax2, eps, gamTemp)
      
      If (a.Eq.b) Then
       gam = gam + gamTemp
      Else
       gam = gam + 2._dp * gamTemp
      End If

   End Do
  End Do

!----------------------------------------------------------- 
! 2.   HpmHpm contribution  --------------------------------
!-----------------------------------------------------------
  Do a = 1, i_hpm                                 
   Do b = a, i_hpm                                         
    mass(1) = mHpm(a)
    mass(2) = mHpm(b)
    m_in    = mSn(i)
    coup(1) = Conjg(C_HpmSlSn(a,j,i))
    coup(2) = Conjg(C_HpmLNu_R(a,k,m))
    coup(3) = Conjg(C_HpmLNu_L(a,k,m))
   
    coup(4) = Conjg(C_HpmSlSn(b,j,i))
    coup(5) = Conjg(C_HpmLNu_R(b,k,m))
    coup(6) = Conjg(C_HpmLNu_L(b,k,m))
   
    If(a.Eq.b) Then
     Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    Else
     Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
     gamTemp = 2._dp * gamTemp  
    End If
 
    gam = gam + Real(gamTemp)

   End Do
  End Do       
  
!-----------------------------------------------------------    
! 3.    WW    contribution  --------------------------------
!-----------------------------------------------------------            

   mass(1) = mW
   m_in    = mSn(i)
   
   coup(1) = C_WSlSn(j,i,1)
   coup(2) = C_WLNu
   coup(3) = 0._dp
   
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   
   gam = gam + Real(gamTemp)

!-----------------------------------------------------------
! 4.  WHpm   contribution  ---------------------------------
!-----------------------------------------------------------
  Do b = 1,i_hpm
   mass(1) = mW
   mass(2) = mHpm(b)
   m_in = mSn(i)
   
   coup(1) = C_WSlSn(j,i,1)
   coup(2) = C_WLNu
   coup(3) = 0._dp
   
   coup(4) = Conjg(C_HpmSlSn(b,j,i))
   coup(5) = Conjg(C_HpmLNu_R(b,k,m))
   coup(6) = Conjg(C_HpmLNu_L(b,k,m))
   
   If(mass(1).Eq.mass(2)) Then
    Call IntegrateVSGoldstone(mass,m_in,r_outcrossed,coup,smin2,smax2,eps, gamTemp)
   Else
    Call IntegrateVS(mass,m_in,r_outcrossed,coup,smin2,smax2,eps, gamTemp)   
   End If

   gam = gam + 2._dp*Real(gamTemp)

  End Do                   

!-----------------------------------------------------------
! 5.   WChi0 contribution  ---------------------------------
!-----------------------------------------------------------
  Do a = 1,i_n
   mass(1) = mN(a)
   mass(2) = mW
   m_in = mSn(i)
   
   coup(1) = C_NSnNu_L(m,a,i)
   coup(2) = C_NSnNu_R(m,a,i)
   coup(3) = Conjg(C_NSlL_R(k,a,j))
   coup(4) = Conjg(C_NSlL_L(k,a,j))
  
   coup(5) = C_WSlSn(j,i,1)
   coup(6) = C_WLNu
   coup(7) = 0._dp
   
   Call IntegrateVF(mass,m_in,r_outcrossed,coup,smin2,smax2,eps, gamTemp)
   
   gam = gam + 2._dp*Real(gamTemp)

  End Do   
  
!-----------------------------------------------------------    
! 6.  HpmChi0  contribution  --------------------------------
!-----------------------------------------------------------                                       
  Do a = 1, i_n
   Do b = 1 , i_hpm
   mass(1) = mN(a)
   mass(2) = mHpm(b)
   m_in = mSn(i)
      
   coup(1) = C_NSnNu_L(m,a,i)
   coup(2) = C_NSnNu_R(m,a,i)
   coup(3) = Conjg(C_NSlL_R(k,a,j))
   coup(4) = Conjg(C_NSlL_L(k,a,j))
   
   coup(5) = Conjg(C_HpmSlSn(b,j,i))
   coup(6) = Conjg(C_HpmLNu_R(b,k,m))
   coup(7) = Conjg(C_HpmLNu_L(b,k,m))
   
   Call IntegrateSF(mass,m_in,r_outcrossed,coup,smin2,smax2,eps, gamTemp)
   
   gam = gam + 2._dp*Real(gamTemp)

   End Do  
  End Do 

  Iname = Iname - 1

 End Subroutine Sneutrino_to_Slepton_lnu 
 
 
 Subroutine Sneutrino_to_Antislepton_lnu(i,j,k,m, mN, mSl, mSn, mC, mf_l    &
                     & , mf_nu, C_NSlL_L, C_NSlL_R, C_NSnNu_L, C_NSnNu_R    &
                     & , C_CSnL_L, C_CSnL_R , C_CNuSl_L, C_CNuSl_R, eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
  Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state antisneutrino
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: m         ! index of neutrino m
  Real(dp), Intent(in) :: mN(:)    ! neutralino masses
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mC(:)    ! chargino masses 
  Real(dp), Intent(in) :: mf_l(:)  ! lepton masses
  Real(dp), Intent(in) :: mf_nu(:) ! neutrino masses
  Complex(dp), Intent(in) :: C_NSlL_L(:,:,:) ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_NSlL_R(:,:,:) !    to neutralinos 
  Complex(dp), Intent(in) :: C_NSnNu_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_NSnNu_R(:,:,:) !    to neutralinos
  Complex(dp), Intent(in) :: C_CSnL_L(:,:,:) ! LR couplings of sneutrinos
  Complex(dp), Intent(in) :: C_CSnL_R(:,:,:) !    to chargninos
  Complex(dp), Intent(in) :: C_CNuSl_L(:,:,:) ! LR couplings of sleptons
  Complex(dp), Intent(in) :: C_CNuSl_R(:,:,:) !    to chargninos
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam     ! partial width

 
 Integer :: a, b, i_n, i_cha
  Real(dp) :: smin, smax,smin2,smax2, rj2, rk2, rm2
  Real(dp) :: r_out(3),r_outcrossed(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_Antislepton_lnu"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.(mSl(j)+ mf_l(k)+mf_nu(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  i_n   = Size(mN)
  i_cha = Size(mC)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_l(k) / mSn(i))**2
  rm2 = (mf_nu(m) / mSn(i))**2
  rj2 = (mSl(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)
  
  smin2 = 2._dp * Sqrt(rm2)
  smax2 = 1._dp + rm2 -  rj2 - rk2 - 2._dp * Sqrt(rj2*rk2)
  
  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2
  
  r_outcrossed(1) = rj2
  r_outcrossed(2) = rm2
  r_outcrossed(3) = rk2

!-----------------------------------------------------------    
!1.   Chi0Chi0 contribution  ------------------------------
!-----------------------------------------------------------                                   
  Do a = 1, i_n                                
   Do b = a, i_n                              
      mass(1) = mN(a)
      mass(2) = mN(b)
      m_in    = mSn(i)

      coup(1) = C_NSnNu_L(m,a,i)
      coup(2) = C_NSnNu_R(m,a,i)
      coup(3) = C_NSlL_L(k,a,j)
      coup(4) = C_NSlL_R(k,a,j)
      
      coup(5) = C_NSnNu_L(m,b,i)
      coup(6) = C_NSnNu_R(m,b,i)
      coup(7) = C_NSlL_L(k,b,j)
      coup(8) = C_NSlL_R(k,b,j)
      Call IntegrateFFLM(mass,m_in,r_outcrossed, coup,smin2, smax2, eps, gamTemp)
      
     If (a.Eq.b) Then 
      gam = gam + gamTemp
     Else
      gam = gam + 2._dp * gamTemp
     End If

   End Do
  End Do
!-----------------------------------------------------------    
! 2.  ChaCha contribution  ---------------------------------
!-----------------------------------------------------------     
  Do a = 1, i_cha                                
   Do b = a, i_cha                              
      mass(1) = mC(a)
      mass(2) = mC(b)
      m_in    = mSn(i)
      
      coup(1) = C_CSnL_L(a,k,i)
      coup(2) = C_CSnL_R(a,k,i)
      coup(3) = C_CNuSl_L(a,m,j)
      coup(4) = C_CNuSl_R(a,m,j)
      
      coup(5) = C_CSnL_L(b,k,i)
      coup(6) = C_CSnL_R(b,k,i)
      coup(7) = C_CNuSl_L(b,m,j)
      coup(8) = C_CNuSl_R(b,m,j)

      Call IntegrateFFLM(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)

     If (a.Eq.b) Then 
      gam = gam + gamTemp
     Else
      gam = gam + 2._dp * gamTemp
     End If

   End Do
  End Do
!-----------------------------------------------------------    
! 3.  ChaChi contribution  ---------------------------------
!-----------------------------------------------------------     
  Do a = 1, i_n                                
   Do b = 1, i_cha                              
      mass(1) = mN(a)
      mass(2) = mC(b)
      m_in    = mSn(i)
      coup(1) = C_NSnNu_L(m,a,i)
      coup(2) = C_NSnNu_R(m,a,i)
      coup(3) = C_NSlL_L(k,a,j)
      coup(4) = C_NSlL_R(k,a,j)

      coup(5) = C_CSnL_L(b,k,i)
      coup(6) = C_CSnL_R(b,k,i)
      coup(7) = C_CNuSl_L(b,m,j)
      coup(8) = C_CNuSl_R(b,m,j)
      Call IntegrateChiChiInterference(mass,m_in,r_out, coup,smin, smax, eps, gamTemp)
      
      gam = gam + 2._dp*Real(gamTemp)

   End Do
  End Do 
 
  Iname = Iname - 1

 End Subroutine Sneutrino_to_Antislepton_lnu 


 Subroutine Sneutrino_to_Slepton_qq(i,j,k,m, mSl, mSn, mHpm, mW, mf_d, mf_u    &
                   &, C_HpmUD_L, C_HpmUD_R ,C_WUD, C_HpmSlSn, C_WSlSn, eps, gam)
 !-----------------------------------------------------------------------------
 ! written by Lukas Mitzka
 !-----------------------------------------------------------------------------
 Implicit None
 
  Integer, Intent(in) :: i         ! index of decaying snutrino
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of up type quark k
  Integer, Intent(in) :: m         ! index of down type m
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mSn(:)   ! sneutrino masses
  Real(dp), Intent(in) :: mHpm(:)  ! higgs masses 
  Real(dp), Intent(in) :: mW       ! suprise surprise Z-mass!!
  Real(dp), Intent(in) :: mf_d(:)  ! d-quark masses
  Real(dp), Intent(in) :: mf_u(:)  ! u-quark masses
  Complex(dp), Intent(in) :: C_HpmUD_L(:,:,:)  ! left coupling of quarks to Hpm
  Complex(dp), Intent(in) :: C_HpmUD_R(:,:,:)  ! right coupling of quarks to Hpm
  Complex(dp), Intent(in) :: C_WUD(:,:)        ! coupling of quarks to W
  Complex(dp), Intent(in) :: C_HpmSlSn(:,:,:) ! coupling of sfermions to Hpm\
  Complex(dp), Intent(in) :: C_WSlSn(:,:,:)     ! coupling of fermions to W
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam     ! partial width
  
  
  Integer :: a, b, i_hpm
  Real(dp) :: smin, smax, rj2, rk2, rm2
  Real(dp) :: r_out(3)
  Real(dp) :: mass(3), m_in
  Complex(dp) :: coup(10), gamTemp

  mass = 0._dp
  coup = 0._dp
  gamTemp = 0._dp
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "Sneutrino_to_Slepton_qq"
  gam = 0._dp
  
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSn(i).Lt.(mSl(j)+mf_u(k)+mf_d(m)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSn(i).gt.(mSl(j)+mW) ) Then
   Iname = Iname - 1
   Return
  End If

  i_hpm = Size(mHpm)
  
  !--------------------------------------
  ! kinematical functions
  ! a and b are the neutralino indizes
  !--------------------------------------
  rk2 = (mf_u(k) / mSn(i))**2
  rm2 = (mf_d(m) / mSn(i))**2
  rj2 = (mSl(j) / mSn(i) )**2
  smin = 2._dp * Sqrt(rk2)
  smax = 1._dp + rk2 -  rj2 - rm2 - 2._dp * Sqrt(rj2*rm2)

  r_out(1) = rj2
  r_out(2) = rk2
  r_out(3) = rm2

!-----------------------------------------------------------    
! 1.   HpmHpm contribution  --------------------------------
!-----------------------------------------------------------     

  Do a = 1, i_hpm                                 
   Do b = a, i_hpm                                         
   mass(1) = mHpm(a)
   mass(2) = mHpm(b)
   m_in    = mSn(i)
   
   coup(1) = Conjg(C_HpmSlSn(a,j,i))
   coup(2) = Conjg(C_HpmUD_R(a,m,k))
   coup(3) = Conjg(C_HpmUD_L(a,m,k))
   
   coup(4) = Conjg(C_HpmSlSn(b,j,i))
   coup(5) = Conjg(C_HpmUD_R(b,m,k))
   coup(6) = Conjg(C_HpmUD_L(b,m,k))
   
   If(a.Eq.b) Then
    Call IntegrateSaSa(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
   Else
    Call IntegrateSaSb(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)
    gamTemp = 2._dp * gamTemp
   End If
   
   gam = gam + Real(gamTemp)

   End Do
  End Do                                

!-----------------------------------------------------------    
! 2.   WW contribution  ------------------------------------
!-----------------------------------------------------------   

   mass(1) = mW
   m_in    = mSn(i)
   coup(1) = C_WSlSn(j,i,1)
   coup(2) = C_WUD(m,k)
   coup(3) = 0._dp
   Call IntegrateVV(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)

   gam = gam + Real(gamTemp)

!-----------------------------------------------------------    
! 3. HpmW  contribution  -----------------------
!-----------------------------------------------------------       
  Do b = 1,i_hpm
   mass(1) = mW
   mass(2) = mHpm(b)
   m_in = mSn(i)

   coup(1) = C_WSlSn(j,i,1)
   coup(3) = C_WUD(m,k)
   coup(2) = 0._dp

   coup(4) = Conjg(C_HpmSlSn(b,j,i))
   coup(5) = Conjg(C_HpmUD_R(b,m,k))
   coup(6) = Conjg(C_HpmUD_L(b,m,k))

   If(mass(1).Eq.mass(2)) Then
    Call IntegrateVSGoldstone(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)                                
   Else
    Call IntegrateVS(mass,m_in,r_out,coup,smin,smax,eps, gamTemp)   
   End If

   gam = gam + 2._dp*Real(gamTemp)

  End Do
 
  gam = 3._dp * gam

  Iname = Iname - 1

 End Subroutine Sneutrino_to_Slepton_qq 
 

End Module Slepton3bodyDecays
