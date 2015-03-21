Module LowEnergy

! load modules
Use Control
Use LoopFunctions
Use StandardModel
Use Couplings
Use LoopCouplings
Use LoopMasses
! load modules

! interfaces

Interface C_10_10p
 Module Procedure C_QdQdLL_AAp
End Interface

Interface C_11_11p
 Module Procedure C_QdQdNuNu_LR
End Interface

Interface Gminus2
 Module Procedure Gminus2a, Gminus2b
End Interface

Interface LtoLpGamma
 Module Procedure  LtoLpGammaMSSM, LtoLpGammaEPS3, LtoLpGammaRPcons
End Interface
! interfaces 

Real(dp) :: L_edm_contri(5,2)
! private variables

Real(dp), Parameter, Private :: e_u = 2._dp/3._dp, e_d=-1._dp/3._dp
Logical, Private :: WriteLLpGamma=.False.
! private variables
! Wilson coefficients at different scales
 Integer, Parameter, Private :: n_i=2
 Real(dp), Private, Save :: Q_out(n_i)
 Complex(dp), Private, Save :: C7_out(n_i,2) = ZeroC, C7p_out(n_i,2) = ZeroC  &
   & , C8_out(n_i,2) = ZeroC, C8p_out(n_i,2) = ZeroC, C9_out(n_i,3,2) = ZeroC &
   & , C9p_out(n_i,3,2) = ZeroC, C10_out(n_i,3,2) = ZeroC                     &
   & , C10p_out(n_i,3,2) = ZeroC, C11_out(n_i,3,2) = ZeroC                    &
   & , C11p_out(n_i,3,2) = ZeroC
 Complex(dp), Private, Save :: &
   &   WC_c7(n_i,3,2,2) = ZeroC, WC_c7p(n_i,3,2,2) = ZeroC           &
   & , WC_c8(n_i,3,2,2) = ZeroC, WC_c8p(n_i,3,2,2) = ZeroC           &
   & , WC_c9p(n_i,3,2,3,2) = ZeroC, WC_c9(n_i,3,2,3,2) = ZeroC       &
   & , WC_c10(n_i,3,2,3,2) = ZeroC, WC_c10p(n_i,3,2,3,2) = ZeroC     &
   & , WC_c11(n_i,3,2,3,2) = ZeroC, WC_c11p(n_i,3,2,3,2) = ZeroC     &
   & , WC_4d_VLL(n_i,3,3) = ZeroC, WC_4d_VRR(n_i,3,3) = ZeroC        &
   & , WC_4d_LR1(n_i,3,3) = ZeroC, WC_4d_LR2(n_i,3,3) = ZeroC        &
   & , WC_4d_SLL1(n_i,3,3) = ZeroC, WC_4d_SRR1(n_i,3,3) = ZeroC      &
   & , WC_4d_SLL2(n_i,3,3) = ZeroC, WC_4d_SRR2(n_i,3,3) = ZeroC      &
   & , WC_2d2l_CS(n_i,3,3,3) = ZeroC, WC_2d2l_CSp(n_i,3,3,3) = ZeroC &
   & , WC_2d2l_CP(n_i,3,3,3) = ZeroC, WC_2d2l_CPp(n_i,3,3,3) = ZeroC
 Logical, Private, Save :: l_wc(n_i) = .False.
! for internal use
 Logical, Private, Save :: WriteDetails = .False.

Contains 


 Subroutine BR_lj_to_3li(j, i, gp, g, Y_l, uL_L, uL_R, mSlep2, RSl, mN, N &
                      & , mSnu2, RSn, mC, U, V, mS02, RS0, mP02, RP0, A_l &
                      & , mu, vevSM, Br)
 !-------------------------------------------------------------------------
 ! calculates the lepton number violating decay of lepton into 3 lighter
 ! leptons
 ! Formulas by E. Arganda and M.J. Herrero, Phys. Rev. D73, 055003 (2006).
 ! Based on routines by E. Arganda and M.J. Herrero
 ! 20.07.2009: adding lepton mixing matrices, thanks to António Figueiredo
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: j, i ! index of decaying lepton and the index of
                              ! daughters
  Real(dp), Intent(in) :: gp, g             ! U(1) and SU(2) gauge couplings
  Complex(dp), Intent(in) :: Y_L(3,3)       ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: uL_L(3,3), uL_R(3,3) ! the corresponding mixing matrices
  Real(dp), Intent(in) :: mSlep2(6)         ! slepton masses squared
  Complex(dp), Intent(in) :: RSl(6,6)       ! slepton mixing matrix
  Real(dp), Intent(in) :: mN(4)             ! neutralino masses
  Complex(dp), Intent(in) :: N(4,4)         ! neutralino mixing matrix 
  Real(dp), Intent(in) :: mSnu2(3)          ! sneutrino masses squared
  Complex(dp), Intent(in) :: RSn(3,3)       ! sneutrino mixing matrix
  Real(dp), Intent(in) :: mC(2)             ! chargino masses
  Complex(dp), Intent(in) :: U(2,2), V(2,2) ! chargino mixing matrices
  Real(dp), Intent(in) :: mS02(2)           ! neutral Higgs boson masses squared
  Real(dp), Intent(in) :: RS0(2,2)          ! neutral Higgs mixing matrix
  Real(dp), Intent(in) :: mP02(2)           ! pseudoscalar Higgs masses squared
  Real(dp), Intent(in) :: RP0(2,2)          ! pseudoscalar Higgs mixing matrix
  Complex(dp), Intent(in) :: A_l(3,3)     ! trilinear Higgs - slepton couplings
  Complex(dp), Intent(in) :: mu             ! mu parameter of the superpotential
  Real(dp), Intent(in) ::  vevSM(2)         ! Higgs vevs (v_d, v_u)
  Real(dp), Intent(out) :: BR
  
  Integer :: i1, i2, i3, i4
  Real(dp) :: mN2(4), mC2(2), fun1, fun2, X_ax, cosW, cosW2, sinW2 &
      & , c_Zee_L, c_Zee_R, fun3, fun4
  Complex(dp) :: c_LNSl_L(3,4,6), c_LNSl_R(3,4,6), QSlepton(6,6), QSnu(3,3) &
            &  , c_CLSn_L(2,3,3), c_CLSn_R(2,3,3), ELn(4,4), ERn(4,4)       &
            &  , ELc(2,2), ERc(2,2), DR(3,4,4), DL(3,4,4)                   &
            &  , WR(3,2,2), WL(3,2,2), GSlepton(3,6,6), bi(1)               &
            &  , c_HLL_L(3,3), c_HLL_R(3,3), GSnu(3,3,3)
  Complex(dp) :: A1L, A1R, A2L, A2R, FL, FR, FLL, FLR, FRL, FRR &
      & , B1L, B1R, B2L, B2R, B3L, B3R, B4L, B4R                &
      & , H2L, H2R, H3L, H3R, wertL, wertR, fun1c, fun2c

  Real(dp), Parameter :: oo576Pi2 = 1._dp / (576._dp * Pi2)
  Complex(dp), Parameter :: mat0(3,3) = 0._dp

  Iname = Iname + 1
  NameOfUnit(Iname) = "BR_lj_to_3li"

  BR = 0._dp

  !------------
  ! couplings
  !------------
  c_LNSl_L = 0._dp
  c_LNSl_R = 0._dp
  Do i1=1,6
   Do i2=1,4
    Call CoupNeutralinoSlepton3(i, i2, i1, gp, g, RSl, uL_L, uL_R, Y_l, N &
                                  &, c_LNSl_L(i,i2,i1), c_LNSl_R(i,i2,i1) )
    Call CoupNeutralinoSlepton3(j, i2, i1, gp, g, RSl, uL_L, uL_R, Y_l, N &
                                  &, c_LNSl_L(j,i2,i1), c_LNSl_R(j,i2,i1) )
   End Do
  End Do

  c_CLSn_L = 0._dp
  c_CLSn_R = 0._dp
  Do i1=1,3
   Do i2=1,2
    Call CoupCharginoSfermion3(i2, j, i1, g, -0.5_dp, RSn, Y_l, mat0, uL_L &
                  & , uL_R, U, V, c_CLSn_L(i2, j, i1), c_CLSn_R(i2, j, i1) )
    Call CoupCharginoSfermion3(i2, i, i1, g, -0.5_dp, RSn, Y_l, mat0, uL_L &
                  & , uL_R, U, V, c_CLSn_L(i2, i, i1), c_CLSn_R(i2, i, i1) )
   End Do
  End Do

  cosW2 = mW2 / mZ2 !g**2 / (gp**2 + g**2)
  cosW = Sqrt(cosW2)
  Do i1=1,4
   Do i2=1,4
    Call CoupNeutralinoZ(i1, i2, N, g, cosW, ELn(i1,i2), ERn(i1,i2) )
   End Do
  End Do

  Do i1=1,2
   Do i2=1,2
    Call CoupCharginoZ(i1, i2, U, V, g, cosW, ELc(i1,i2), ERc(i1,i2) )
   End Do
  End Do

  sinW2 = 1._dp - cosW2
  Do i1=1,6
   Do i2=1,6
    Call CoupSleptonZ(i1, i2, g, sinW2, RSl, QSlepton(i1,i2) )
   End Do
  End Do
  Do i1=1,3
   Do i2=1,3
    Call CoupSneutrinoZ3(i1, i2, g, sinW2, RSn, QSnu(i1,i2) )
   End Do
  End Do
  Qslepton = - Qslepton ! sign difference in paper by Arganda/Herrero
  Qsnu = - Qsnu ! sign difference in paper by  Arganda/Herrero

  Call CoupFermionZ(-0.5_dp, -1._dp, g, sinW2, c_Zee_L,c_Zee_R)

  DR = 0._dp
  DL = 0._dp
  Do i2=1,4
   Do i3=1,4
    Do i1=1,2
     Call CoupNeutralinoScalar(i2,i3,i1, N, RS0, gp, g &
                             &, DL(i1,i2,i3), DR(i1,i2,i3) )
    End Do
     Call CoupNeutralinoPseudoScalar(i2,i3,2, N, RP0, gp, g &
                             &, DL(3,i2,i3), DR(3,i2,i3) )
   End Do
  End Do

  WR = 0._dp
  WL = 0._dp
  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
     Call CoupCharginoScalar(i1,i2,i3, U, V, RS0, g, WL(i3,i2,i1), WR(i3,i2,i1))
    End Do
    Call CoupCharginoPseudoScalar(i1,i2,2, U, V, RP0, g &
                                & , WL(3,i2,i1), WR(3,i2,i1))
   End Do
  End Do

  GSlepton = 0._dp
  bi = mu
  Do i1=1,6
   Do i2=1,6
    Do i3=1,2
     Call CoupScalarSfermion3MSSM_3(i3, i1, i2, RS0, -0.5_dp, -1._dp, Y_L &
                 & , Rsl, A_l, mu, vevSM, gp, g, GSlepton(i3, i1, i2) )
    End Do
    Call CoupPseudoScalarSfermion3_3(2, i1, i2, RP0, -0.5_dp, Y_l, Rsl &
                 & , A_l, bi, GSlepton(3, i1, i2) )
   End Do
  End Do

  GSnu = 0._dp
  bi = mu
  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Call CoupScalarSfermion3MSSM_3(i3, i1, i2, RS0, 0.5_dp, 0._dp, mat0 &
                 & , Rsn, mat0, mu, vevSM, gp, g, GSnu(i3, i1, i2) )
    End Do
   End Do
  End Do

  c_HLL_L = 0._dp
  c_HLL_R = 0._dp
  Do i1=1,3
   Do i2=1,2
    Call CoupFermionScalar3(i1, i1, i2, -0.5_dp, Y_l, uL_L, uL_R, RS0 &
                          & , c_HLL_L(i2,i1), c_HLL_R(i2,i1) )
   End Do
   Call CoupFermionPseudoScalar3(i1, i1, 2, -0.5_dp, Y_l, uL_L, uL_R, RP0 &
                          & , c_HLL_L(3,i1), c_HLL_R(3,i1) )
  End Do

  !---------------------------------------------
  ! Penguin contributions
  !---------------------------------------------
  ! neutralinos
  !--------------
  A1L = 0._dp
  A1R = 0._dp
  A2L = 0._dp
  A2R = 0._dp
  mN2 = mN**2

  Do i1=1,6
   Do i2=1,4
    X_ax = mN2(i2)/mSlep2(i1)
    fun1 = ( 2._dp - 9._dp*X_ax + 18._dp*X_ax**2 - 11._dp*X_ax**3 & 
        & + 6._dp*X_ax**3 * Log(X_ax) ) /  (1._dp-X_ax)**4
    A1L = A1L &
      & + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1)) * fun1 / mSlep2(i1)
    A1R = A1R &
      & + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1)) * fun1 / mSlep2(i1)

    fun1 = 2._dp * F2(X_ax)
    fun2 = 2._dp * F4(X_ax)
    A2L = A2L + ( ( c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))      & 
        &         + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))      &
        &           *mf_l(i)/mf_l(j) ) * fun1                         &
        &       + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))        &
        &         * mN(i2)/mf_l(j)*fun2 ) / mSlep2(i1)
    A2R = A2R + ( ( c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))      & 
        &         + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))      &
        &           *mf_l(i)/mf_l(j) ) * fun1                         &
        &       + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))        &
        &         * mN(i2)/mf_l(j)*fun2 ) / mSlep2(i1)
   End Do
  End Do

  !--------------
  ! charginos
  !--------------
  mC2 = mC**2
  Do i1=1,3
   Do i2=1,2
    X_ax=mC2(i2)/mSnu2(i1)
    fun1 = ( 16._dp - 45._dp*X_ax + 36._dp * X_ax**2 - 7._dp*X_ax**3 &
         & + 6._dp * (2._dp-3._dp*X_ax) * Log(X_ax) ) / (1._dp-X_ax)**4
    A1L = A1L - c_CLSn_R(i2,i,i1) * Conjg(c_CLSn_R(i2,j,i1)) * fun1 / mSnu2(i1)
    A1R = A1R - c_CLSn_L(i2,i,i1) * Conjg(c_CLSn_L(i2,j,i1)) * fun1 / mSnu2(i1)
    fun1 = 2._dp * F1(X_ax)
    fun2 = 2._dp * F3(X_ax)
    A2L = A2L - ( ( c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1))                 & 
        &         + c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1))     &
        &           *mf_l(i)/mf_l(j) ) *fun1                       &
        &       + c_CLSn_L(i2,i,i1) * Conjg(c_CLSn_R(i2,j,i1))     &
        &          * mC(i2) / mf_l(j) * fun2 ) / mSnu2(i1)
    A2R = A2R - ( ( c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1))                 & 
        &         + c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1))     &
        &           *mf_l(i)/mf_l(j) ) *fun1                       &
        &       + c_CLSn_R(i2,i,i1) * Conjg(c_CLSn_L(i2,j,i1))     &
        &          * mC(i2) / mf_l(j) * fun2 ) / mSnu2(i1)
   End Do
  End Do

  A1L = oo576Pi2 * A1L 
  A1R = oo576Pi2 * A1R 
  A2L = oo32Pi2 * A2L 
  A2R = oo32Pi2 * A2R 

  !-----------
  ! Z-penguin
  !-----------
  FL = 0._dp
  FR = 0._DP

  Do i1=1,6
   Do i2=1,4
    Do i3=1,4
     fun1 = 0.5_dp * vertexC0tilde(mSlep2(i1),mN2(i2),mN2(i3))
     fun2 = mN(i2) * mN(i3) * C0_3m(mSlep2(i1),mN2(i2),mN2(i3))
     FL = FL - c_LNSl_R(i,i3,i1)*Conjg(c_LNSl_R(j,i2,i1))             & 
      &   * ( ERn(i3,i2) * fun1 - ELn(i3,i2) * fun2 )
     FR = FR - c_LNSl_L(i,i3,i1)*Conjg(c_LNSl_L(j,i2,i1))             & 
      &   * ( ELn(i3,i2) * fun1 - ERn(i3,i2) * fun2 )
    End Do
   End Do
  End Do

   Do i1=1,4
    Do i2=1,6
     Do i3=1,6
      fun1 = 0.5_dp * Qslepton(i2,i3)  &
           &        * vertexC0tilde(mN2(i1),mSlep2(i2),mSlep2(i3)) 
      FL = FL - c_LNSl_R(i,i1,i2) * Conjg(c_LNSl_R(j,i1,i3)) * fun1
      FR = FR - c_LNSl_L(i,i1,i2) * Conjg(c_LNSl_L(j,i1,i3)) * fun1
     End Do
    End Do
   End Do

   Do i1=1,6
    Do i2=1,4
     fun1 = B1(0._dp, mN2(i2), mSlep2(i1) )
     FL = FL - c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1)) * c_Zee_L * fun1
     FR = FR - c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1)) * c_Zee_R * fun1
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mSnu2(i1), mC2(i2), mC2(i3) )
      fun2 = C0_3m(mSnu2(i1), mC2(i2), mC2(i3) )
      FL = FL - c_CLSn_R(i3,i,i1)*Conjg(c_CLSn_R(i2,j,i1) )             & 
         &  * ( ERc(i3,i2) * fun1 - ELc(i3,i2) * mC(i2) * mC(i3) * fun2 )
      FR = FR - c_CLSn_L(i3,i,i1)*Conjg(c_CLSn_L(i2,j,i1) )             & 
         &  * ( ELc(i3,i2) * fun1 - ERc(i3,i2) * mC(i2) * mC(i3) * fun2 )
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mC2(i3), mSnu2(i1), mSnu2(i2) )
      FL = FL - c_CLSn_R(i3,i,i1) * Conjg(c_CLSn_R(i3,j,i2)) &
         &       * Qsnu(i1,i2) * fun1
      FR = FR - c_CLSn_L(i3,i,i1) * Conjg(c_CLSn_L(i3,j,i2)) &
         &       * Qsnu(i1,i2) * fun1
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     fun1 = B1(0._dp, mC2(i2), mSnu2(i1) )
     FL = FL - c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1)) * c_Zee_L * fun1
     FR = FR - c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1)) * c_Zee_R * fun1
    End Do
   End Do

  FL = oo16pi2 * FL
  FR = oo16pi2 * FR

  fun1 = 1._dp / (g**2 * sinW2 * mZ2)
  FLL = - fun1 * c_Zee_L * FL 
  FLR = - fun1 * c_Zee_R * FL 
  FRL = - fun1 * c_Zee_L * FR 
  FRR = - fun1 * c_Zee_R * FR

  !-------------------
  ! box contributions
  !-------------------
  B1L = 0._dp
  B1R = 0._dp
  B2L = 0._dp
  B2R = 0._dp
  B3L = 0._dp
  B3R = 0._dp
  B4L = 0._dp
  B4R = 0._dp


  Do i1=1,6
   Do i2=1,6
    Do i3=1,4
     Do i4=1,4
      fun1 = D27_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSlep2(i2) )
      fun2 = mN(i3) * mN(i4)  &
         &   * D0_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSlep2(i2) )

      B1L = B1L + 2._dp * fun1 * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1)) &
          &                    * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2)) & 
          &      + fun2 * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1)               &
          &             * Conjg(c_LNSl_R(j,i3,i1)) * Conjg(c_LNSl_R(i,i3,i2))

      B1R = B1R + 2._dp * fun1 * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1)) &
          &                    * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2)) & 
          &      + fun2 * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1)               &
          &             * Conjg(c_LNSl_L(j,i3,i1)) * Conjg(c_LNSl_L(i,i3,i2))

      B2L = B2L + fun1 * ( c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &                * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))      &
          &              - c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &                * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i3,i2))      &
          &              + c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(i,i3,i2))        &
          &                *Conjg(c_LNSl_R(j,i3,i1)) * c_LNSl_L(i,i4,i1) )     &
          &     - 0.5_dp * fun2 * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1)) &
          &                     * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))

      B2R = B2R + fun1 * ( c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &                * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))      &
          &              - c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &                * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i3,i2))      &
          &              + c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(i,i3,i2))        &
          &                *Conjg(c_LNSl_L(j,i3,i1)) * c_LNSl_R(i,i4,i1) )     &
          &     - 0.5_dp * fun2 * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1)) &
          &                     * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))

      B3L = B3L + fun2 * ( c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &              * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))        &
          &              + 0.5_dp * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1)      &
          &               * Conjg(c_LNSl_R(j,i3,i1))* Conjg(c_LNSl_R(i,i3,i2)) )

      B3R = B3R + fun2 * ( c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &              * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))        &
          &              + 0.5_dp * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1)      &
          &               * Conjg(c_LNSl_L(j,i3,i1))* Conjg(c_LNSl_L(i,i3,i2)) )

      B4L = B4L + fun2 * Conjg(c_LNSl_R(j,i3,i1)) * Conjg(c_LNSl_R(i,i3,i2))  &
         &            * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1) / 8._dp

      B4R = B4R + fun2 * Conjg(c_LNSl_L(j,i3,i1)) * Conjg(c_LNSl_L(i,i3,i2))  &
          &            * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1) / 8._dp

     End Do
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Do i4=1,2
      fun1 = D27_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSnu2(i2) )
      fun2 = mC(i3) * mC(i4) * D0_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSnu2(i2))

      B1L = B1L + 2._dp * fun1 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1)) &
          &                    * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B2L = B2L + fun1 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1))          &
          &            * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))          & 
          &     - 0.5_dp * fun2 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1)) &
          &                     * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

      B3L = B3L + fun2 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1))    &
          &            * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B1R = B1R + 2._dp * fun1 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1)) &
          &                    * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

      B2R = B2R + fun1 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1))          &
          &            * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))          & 
          &     - 0.5_dp * fun2 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1)) &
          &                     * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B3R = B3R + fun2 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1))    &
          &            * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

     End Do
    End Do
   End Do
  End Do

  fun1 = oo16pi2 / (g**2*sinW2)
  B1L = fun1 * B1L
  B1R = fun1 * B1R
  B2L = fun1 * B2L
  B2R = fun1 * B2R
  B3L = fun1 * B3L
  B3R = fun1 * B3R
  B4L = fun1 * B4L
  B4R = fun1 * B4R
  !---------------------------
  ! Higgs contributions
  !---------------------------
  H2L = 0._dp
  H2R = 0._dp
  H3L = 0._dp
  H3R = 0._dp

  Do i1=1,3
   wertL = 0._dp
   wertR = 0._dp
   Do i2=1,4
    Do i3=1,4
     fun1 = B0(0._dp, mN2(i2), mN2(i3))
     Do i4=1,6
      fun2 = C0_3m(mSlep2(i4), mN2(i2), mN2(i3) )
      fun3 = vertexC12(mSlep2(i4), mN2(i2), mN2(i3) )
      fun4 = vertexC11(mSlep2(i4), mN2(i2), mN2(i3) )

      wertL = wertL &
        &   + ( fun1 + mSlep2(i4) * fun2 +mf_l2(j) * fun3                    &
        &       + mf_l2(i)*( fun4 - fun3 ) )                                 &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mf_l(j) *( fun4 + fun2 )                             &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mN(i3) *( fun4 - fun3 + fun2 )                       &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(j) * mN(i3) * fun3                                        &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mN(i2) * ( fun4 - fun3 )                             &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(j)* mN(i2) * ( fun3 + fun2 )                              &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mN(i2) * mN(i3) * fun2                                         &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4))

      wertR = wertR &
        &   + ( fun1 + mSlep2(i4) * fun2 +mf_l2(j) * fun3                    &
        &       + mf_l2(i)*( fun4 - fun3 ) )                                 &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mf_l(j) *( fun4 + fun2 )                             &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mN(i3) *( fun4 - fun3 + fun2 )                       &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(j) * mN(i3) * fun3                                        &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mN(i2) * ( fun4 - fun3 )                             &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(j)* mN(i2) * ( fun3 + fun2 )                              &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mN(i2) * mN(i3) * fun2                                         &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4))
     End Do
    End Do
   End Do

   Do i2=1,6
    Do i3=1,6
     Do i4=1,4
      fun1 = vertexC11(mN2(i4), mSlep2(i2), mSlep2(i3))
      fun2 = vertexC12(mN2(i4), mSlep2(i2), mSlep2(i3))
      fun3 = C0_3m(mN2(i4), mSlep2(i2), mSlep2(i3))

      wertL = wertL                                                            &
          &  - ( mf_l(i) * (fun1-fun2)                                         &
          &              *c_LNSl_R(i,i4,i2) *Conjg(c_LNSl_R(j,i4,i3))          &
          &    + mf_l(j) * fun2 * c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_L(j,i4,i3)) &
          &    - mN(i4) * fun3 * c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(j,i4,i3))  &
          &    ) * GSlepton(i1,i2,i3)
      wertR = wertR                                                            &
          &  - ( mf_l(i) * (fun1-fun2)                                         &
          &              *c_LNSl_L(i,i4,i2) *Conjg(c_LNSl_L(j,i4,i3))          &
          &    + mf_l(j) * fun2 * c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_R(j,i4,i3)) &
          &    - mN(i4) * fun3 * c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(j,i4,i3))  &
          &    ) * GSlepton(i1,i2,i3)
     End Do
    End Do
   End Do

   Do i2=1,6
    Do i3=1,4
     fun1 = B1(0._dp, mN2(i3), mSlep2(i2))
     fun2 = B0(0._dp, mN2(i3), mSlep2(i2))

     wertL = wertL                                                     &
         & + (-mf_l2(i) * fun1                                         &
         &              * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2)) &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  & 
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             *c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))   & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             *c_LNSl_L(i,i3,i2)*Conjg(c_LNSl_R(j,i3,i2))     &
         &   ) * c_HLL_L(i1,j) / (mf_l2(i) - mf_l2(j) )                &
         & + (-mf_l2(j) * fun1                                         &
         &              * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2)) & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   ) * c_HLL_L(i1,i) / (mf_l2(j)-mf_l2(i))
     wertR = wertR                                                     &
         & + (-mf_l2(i) * fun1                                         &
         &              * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2)) &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  & 
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             *c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))   & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             *c_LNSl_R(i,i3,i2)*Conjg(c_LNSl_L(j,i3,i2))     &
         &   ) * c_HLL_R(i1,j) / (mf_l2(i) - mf_l2(j) )                &
         & + (-mf_l2(j) * fun1                                         &
         &              * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2)) & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   ) * c_HLL_R(i1,i) / (mf_l2(j)-mf_l2(i))
      End Do
   End Do

   Do i2=1,2
    Do i3=1,2
     fun1 = B0(0._dp, mC2(i2), mC2(i3))
     Do i4=1,3
      fun2 = C0_3m(mSnu2(i4), mC2(i2), mC2(i3))
      fun3 = vertexC12(mSnu2(i4), mC2(i2), mC2(i3))
      fun4 = vertexC11(mSnu2(i4), mC2(i2), mC2(i3))

      wertL = wertL                                                         &
        & + ( fun1 + mSnu2(i4) * fun2 + mf_l2(j) * fun3                     &
        &   + mf_l2(i) * (fun4 -fun3) )                                     &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(i) * mf_l(j) * (fun4  +fun2)                             &
        &     *c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4))  &
        &   + mf_l(i) * mC(i3) * (fun4 - fun3 + fun2)                       &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(j) * mC(i3) * fun3                                       &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(i) * mC(i2) *(fun4 -fun3)                                &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(j) * mC(i2) * (fun3 + fun2)                              &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mC(i2) * mC(i3) * fun2                                        &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4))

      wertR = wertR                                                         &
        & + ( fun1 + mSnu2(i4) * fun2 + mf_l2(j) * fun3                     &
        &   + mf_l2(i) * (fun4 -fun3) )                                     &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(i) * mf_l(j) * (fun4  +fun2)                             &
        &     *c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4))  &
        &   + mf_l(i) * mC(i3) * (fun4 - fun3 + fun2)                       &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(j) * mC(i3) * fun3                                       &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(i) * mC(i2) *(fun4 -fun3)                                &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(j) * mC(i2) * (fun3 + fun2)                              &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mC(i2) * mC(i3) * fun2                                        &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4))
     End Do
    End Do
   End Do

   If (i1.Lt.3) Then ! pseudoscalar does not contribute
    Do i2=1,3
     Do i3=1,3
      Do i4=1,2
       fun1 = vertexC11(mC2(i4), mSnu2(i2), mSnu2(i3))
       fun2 = vertexC12(mC2(i4), mSnu2(i2), mSnu2(i3))
       fun3 = C0_3m(mC2(i4), mSnu2(i2), mSnu2(i3))

       wertL = wertL                                                        &
        & - ( mf_l(i) * (fun1 -fun2)                                        &
        &             * c_CLSn_R(i4,i,i2) *Conjg(c_CLSn_R(i4,j,i3))         & 
        &   + mf_l(j) * fun2 * c_CLSn_L(i4,i,i2) * Conjg(c_CLSn_L(i4,j,i3)) & 
        &   - mC(i4) * fun3 * c_CLSn_L(i4,i,i2) * Conjg(c_CLSn_R(i4,j,i3))  &
        &   ) * GSnu(i1,i2,i3)

       wertR = wertR                                                        &
        & - ( mf_l(i) * (fun1 -fun2)                                        &
        &             * c_CLSn_L(i4,i,i2) *Conjg(c_CLSn_L(i4,j,i3))         & 
        &   + mf_l(j) * fun2 * c_CLSn_R(i4,i,i2) * Conjg(c_CLSn_R(i4,j,i3)) & 
        &   - mC(i4) * fun3 * c_CLSn_R(i4,i,i2) * Conjg(c_CLSn_L(i4,j,i3))  &
        &   ) * GSnu(i1,i2,i3)
      End Do
     End Do
    End Do 
   End If

   Do i2=1,3
    Do i3=1,2
     fun1 = B1(0._dp, mC2(i3), mSnu2(i2))
     fun2 = B0(0._dp, mC2(i3), mSnu2(i2))

     wertL = wertL                                                   &
     & + ( - mf_l2(i) * fun1                                         &
     &       * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))          & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   + mf_l(j) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   )*c_HLL_L(i1,j) / (mf_l2(i)-mf_l2(j))                       &
     & + ( - mf_l2(j) * fun1                                         &
     &                * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2)) & 
     &   + mf_l(j) * mC(i3) * fun2                                   & 
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  & 
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   ) * c_HLL_L(i1,i) / (mf_l2(j)-mf_l2(i))

     wertR = wertR                                                   &
     & + ( - mf_l2(i) * fun1                                         &
     &       * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))          & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   + mf_l(j) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   )*c_HLL_R(i1,j) / (mf_l2(i)-mf_l2(j))                       &
     & + ( - mf_l2(j) * fun1                                         &
     &                * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2)) & 
     &   + mf_l(j) * mC(i3) * fun2                                   & 
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  & 
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   ) * c_HLL_R(i1,i) / (mf_l2(j)-mf_l2(i))
    End Do
   End Do

   If (i1.Lt.3) Then
    fun1c = c_HLL_R(i1,i) / mS02(i1)
    fun2c = c_HLL_L(i1,i) / mS02(i1)
   Else
    fun1c = c_HLL_R(i1,i) / mP02(2)
    fun2c = c_HLL_L(i1,i) / mP02(2)
   End If
   H2L = H2L + 0.5_dp * fun1c * wertL 
   H2R = H2R + 0.5_dp * fun2c * wertR 
   H3L = H3L - fun2c * wertL 
   H3R = H3R - fun1c * wertR 
  End Do
  !------------------------------------------
  ! adding Bbox and Higgs contributions
  ! as the contribute to the same operator
  !------------------------------------------
  fun1 = oo16pi2 / (g**2*sinW2)

  B2L = B2L + fun1 * H2L
  B2R = B2R + fun1 * H2R
  B3L = B3L + fun1 * H3L
  B3R = B3R + fun1 * H3R


  !-------------------------
  ! now the branching ratio
  !-------------------------
  BR = Abs(A1L)**2 + Abs(A1R)**2                                         &
   & - 4._dp* Real(A1L*Conjg(A2R)+A2L*Conjg(A1R),dp)                     &
   & + (Abs(A2L)**2 + Abs(A2R)**2)                                       &
   &     * (16._dp * Log(mf_l(j)/mf_l(i)) - 22._dp ) /3._dp              &
   & + ( Abs(B1L)**2 + Abs(B1R)**2 ) / 6._dp                             &
   & + ( Abs(B2L)**2 + Abs(B2R)**2 ) / 3._dp                             &
   & + ( Abs(B3L)**2 + Abs(B3R)**2 ) / 24._dp                            &
   & + 6._dp * ( Abs(B4L)**2 + Abs(B4R)**2 )                             &
   & - Real( B3L * Conjg(B4L) + B3R * Conjg(B4R),dp)                     &
   & + 2._dp * Real( A1L * Conjg(B1L) + A1R * Conjg(B1R)                 &
   &               + A1L * Conjg(B2L) + A1R * Conjg(B2R), dp) / 3._dp    &
   & - 4._dp * Real( A2R * Conjg(B1L) + A2L * Conjg(B1R)                 &
   &               + A2L * Conjg(B2R) + A2R * Conjg(B2L), dp) / 3._dp    &
   & + ( 2._dp * (Abs(FLL)**2 + Abs(FRR)**2) + Abs(FLR)**2 + Abs(FRL)**2 &
   &   + 2._dp * Real( B1L * Conjg(FLL) + B1R * Conjg(FRR)               &
   &                 + B2L * Conjg(FLR) + B2R * Conjg(FRL), dp )         &
   &   + 4._dp * Real( A1L * Conjg(FLL) + A1R * Conjg(FRR), dp )         &
   &   + 2._dp * Real( A1L * Conjg(FLR) + A1R * Conjg(FRL), dp )         &
   &   - 8._dp * Real( A2R * Conjg(FLL) + A2L * Conjg(FRR), dp )         &
   &   - 4._dp * Real( A2L * Conjg(FRL) + A2R * Conjg(FLR), dp ) ) / 3._dp

  !----------------------------------------------------------------------
  ! taking alpha(Q=0) instead of alpha(m_Z) as this contains most of the
  ! running of the Wilson coefficients
  !----------------------------------------------------------------------
  BR = oo32pi * Alpha**2 * mf_l(j)**5 * BR

  If (j.Eq.2) Then
   BR = BR / GammaMu
  Else If  (j.Eq.3) Then
   BR = BR / GammaTau
  End If

  Iname = Iname - 1

 End Subroutine BR_lj_to_3li


 Real(dp) Function Bm_to_l_nu(i, j, mH2, tanb, RSpm, Yd, RuL, RdR, Yl, vevSM, ratio)
 !-------------------------------------------------------------------
 ! calculates the branching ratio of the rare decay B- -> tau nu
 ! based on the formula by W.Hou, PRD48 (1993) 2342 
 !-------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i           ! index of the lepton generation
  Integer, Intent(in) :: j           ! j=1 -> B_u, j=2-> B_c
  Real(dp), Intent(in) :: mH2, tanb  ! m^2_H+, tan(beta)
  Complex(dp), Intent(in), Optional :: RSpm(2,2) & ! H+ mixing matrix
                                   & , Yd(3,3)   & ! d-quark Yukawas
                                   & , RuL(3,3)  & ! L u-quark mixing matrix
                                   & , RdR(3,3)  & ! R d-quark mixing matrix
                                   & , Yl(3,3)     ! lepton Yukawa coupling
  Real(dp), Intent(in), Optional :: vevSM(2)       !
  Logical, Intent(in), Optional :: ratio

  Real(dp) :: r_fac, pre_fac
  Complex(dp) :: cL, cR
  Complex(dp), Parameter :: Zero3(3,3) = (0._dp,0._dp)

  If (    Present(RSpm).And.Present(Yd).And.Present(RuL).And.Present(RdR) &
     .And.Present(Yl).And.Present(vevSM)) Then
   Call CoupChargedScalarFermion3(2, 3, j, RSpm, Yd, id3C, RdR, Zero3 &
                                    &, RuL, id3C, cL, cR)
   r_fac = 0.25_dp * Abs(cL * Yl(i,i) / (CKM(j,3) * RSpm(2,1)))**2 &
         &          * vevSM(1)**4 / (mf_l2(i) * mf_d_mZ(3)**2 ) 
  Else
   r_fac = 1
  End If

  If (Present(ratio)) Then
   If (ratio) Then
    pre_fac = 1._dp
   Else
    pre_fac =  oo8pi * G_F**2 * mf_l2(i) * MassBm(j)                &
            &   * (1._dp -  mf_l2(i) / MassBm(j)**2)**2             &
            &   * Abs(CKM(j,3))**2 * FBhatB(j)**2 *  tauBm(j) / hbar
   End If
  Else
   pre_fac =  oo8pi * G_F**2 * mf_l2(i) * MassBm(j)                 &
            &   * (1._dp -  mf_l2(i) / MassBm(j)**2)**2             &
            &   * Abs(CKM(j,3))**2 * FBhatB(j)**2 *  tauBm(j) / hbar
  End If

  Bm_to_l_nu = pre_fac * (1._dp - r_fac * (tanb * MassBm(j))**2 / mH2)**2

 End Function Bm_to_l_nu


 Subroutine B_to_Q_Gamma(I_f, CKM, c7, c7p, c8, c8p, Bratio             &
       & , A_CP, i_scheme, NNLO_SM_in)
 !-----------------------------------------------------------------------------
 ! gives 10^4 Br(b->s gamma)
 ! input:
 !  I_f .................. if I_f==1 -> Q=d, I_f==1 -> Q=s
 !  CKM(j,j) ............. CKM mixing matrix 
 !  Wilson coefficients at m_W (including the various contributions)
 !                              all are optional
 !  c7(i) ................ C_7: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c7p(i) ............... C'_7: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 !  c8(i) ................ C_8: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c8p(i) ............... C'_8: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 !  i_scheme (optional) ..  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 ! output:
 !  Bratio ............... 10^4 BR( b -> Q gamma)
 !           i_scheme=i_t:  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 ! written by Werner Porod, 11.11.02
 ! 13.03.02: -adding the possiblity to switch between d and s by the switch if
 !           -adding gluino and neutralino contributions
 !           -creating interface for neutral Higgs boson, are set to 0 for
 !            the moment being
 ! 02.10.03: changing numbers in front of R_7,8 as given in coll. with 
 !           Enrico and Tobias, adding i_scheme for the various possiblities
 !           of the a^x_y - arrays, default = 3
 !           i_scheme=i_t:  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 ! 23.08.2008: adding  i_t = 0, using the formula of
 !                     E.Lunghi, J.Matias, hep-ph/0612166, eq.18 still missing
 !                      1-loop part including the NNLO SM value of 2.98
 ! 16.02.2010: adding a different NNLO SM value as optional input for option
 !             it=0
 ! 16.07.2010: changing input from mixing matices to couplings with the idea to
 !             enlarge the number of models where this routine can be used
 ! 14.03.2014: taking Wilson coefficients as input
 !-------------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: I_f
  Complex(dp), Intent(in) :: CKM(3,3), c7(2), c7p(2), c8(2), c8p(2)

  Real(dp), Intent(out) :: Bratio
  Real(dp), Intent(out), Optional :: A_CP
  Integer, Intent(in), Optional :: i_scheme
  Real(dp), Intent(in), Optional :: NNLO_SM_in

  Integer :: i_t
  Real(dp) :: NNLO_SM
  Complex(dp) :: r7, r7p, r8, r8p, epsq, epsqC
  Complex(dp) :: delta_C7_0, delta_C7p_0, delta_C8_0, delta_C8p_0, delta_C7_1 &
               & , delta_C7p_1, delta_C8_1, delta_C8p_1
  Real(dp), Parameter :: T3=-0.5_dp, ed = -1._dp / 3._dp
  Real(dp), Dimension(5), Parameter ::                                         &
    &   a = (/  7.8221_dp, 6.9120_dp, 8.1819_dp, 7.1714_dp, 7.9699_dp /)       &
    & , a_77 = (/ 0.8161_dp, 0.8161_dp, 0.8283_dp, 0.8283_dp, 0.9338_dp /)     &
    & , ar_7 = (/ 4.8802_dp, 4.5689_dp, 4.9228_dp, 4.6035_dp, 5.3314_dp /)     &
    & , ai_7 = (/ 0.3546_dp, 0.2167_dp, 0.3322_dp, 0.2029_dp, 0.0_dp    /)     &
    & , a_88 = (/ 0.0197_dp, 0.0197_dp, 0.0986_dp, 0.0986_dp, 0.0066_dp /)     &
    & , ar_8 = (/ 0.5680_dp, 0.5463_dp, 0.7810_dp, 0.7600_dp, 0.4498_dp /)     &
    & , ai_8 = (/-0.0987_dp,-0.1105_dp,-0.0963_dp,-0.1091_dp, 0.0_dp /)        &
    & , a_eps = (/ 0.4384_dp, 0.3787_dp, 0.8598_dp, 0.7097_dp, 0.0_dp /)       &
    & , ar_eps = (/-1.6981_dp,-2.6679_dp,-1.3329_dp,-2.4935_dp, 0.0_dp /)      &
    & , ai_eps = (/ 2.4997_dp, 2.8956_dp, 2.5274_dp, 2.9127_dp, 0.0_dp /)      &
    & , ar_87 = (/ 0.1923_dp, 0.1923_dp, 0.2025_dp, 0.2025_dp, 0.1576_dp /)    &
    & , ar_7eps = (/-0.7827_dp,-1.0940_dp,-0.8092_dp,-1.1285_dp, 0.0_dp /)     &
    & , ar_8eps = (/ -0.0601_dp,-0.0819_dp,-0.0573_dp,-0.0783_dp, 0.0_dp /)    &
    & , ai_87 = (/-0.0487_dp,-0.0487_dp,-0.0487_dp,-0.0487_dp, 0.0_dp /)       &
    & , ai_7eps = (/-0.9067_dp,-1.0447_dp,-0.9291_dp,-1.0585_dp, 0.0_dp /)     &
    & , ai_8eps = (/-0.0661_dp,-0.0779_dp,-0.0637_dp,-0.0765_dp, 0.0_dp /)

  Iname = Iname + 1
  NameOfUnit(Iname) = "B_to_Q_Gamma"

  !---------------
  ! couplings
  !---------------
  r7 = c7(1) / c7(2)
  r7p = c7p(1) / c7(2)
  r8 = c8(1) / c8(2)
  r8p = c8p(1) / c8(2)

  epsq = Conjg(CKM(1,I_f)) * CKM(1,3) / (Conjg(CKM(3,I_f)) * CKM(3,3))
  epsqC = Conjg(epsq)

  If (Present(i_scheme)) Then
   i_t = i_scheme
  Else
   i_t = 3
  End If

  If (i_t.Ne.0) Then
   Bratio = a(i_t) + a_77(i_t) * ( Abs(r7)**2 + Abs(r7p)**2)                    &
       & + ar_7(i_t) * Real(r7,dp) + ai_7(i_t) * Aimag(r7)                     &
       & + a_88(i_t) * ( Abs(r8)**2 + Abs(r8p)**2)                             &
       & + ar_8(i_t) * Real(r8,dp) + ai_8(i_t) * Aimag(r8)                     &
       & + a_eps(i_t) * Abs(epsq)**2 + ar_eps(i_t) * Real(epsq,dp)             &
       & + ai_eps(i_t) * Aimag(epsq)                                           &
       & + ar_87(i_t) * Real(r8*Conjg(r7)+r8p*Conjg(r7p),dp)                   &
       & + ai_87(i_t) * Aimag(r8*Conjg(r7)+r8p*Conjg(r7p))                     &
       & + ar_7eps(i_t) * Real(r7*epsqC,dp) + ar_8eps(i_t) * Real(r8*epsqC,dp) &
       & + ai_7eps(i_t) * Aimag(r7*epsqC) + ai_8eps(i_t) * Aimag(r8*epsqC) 

  !-----------------------
  ! is 10^4 BR
  !-----------------------
   Bratio = 2.567e-1_dp * Abs(Conjg(CKM(3,I_f))*CKM(3,3)/CKM(2,3))**2 * Bratio

  Else ! E.Lunghi, J.Matias, hep-ph/0612166, eq.18 still missing 1-loop part
   delta_C7_0 = C7(1) - C7(2) ! Sum(C7(3:7))
   delta_C7p_0 = C7p(1)
   delta_C8_0 = C8(1) - C8(2) ! Sum(C8(3:7))
   delta_C8p_0 = C8p(1)
   delta_C7_1 = 0._dp
   delta_C7p_1 = 0._dp
   delta_C8_1 = 0._dp
   delta_C8p_1 = 0._dp
   If (Present(NNLO_SM_in)) Then
    NNLO_SM = NNLO_SM_in
   Else
    NNLO_SM = 2.98_dp
   End If

   Bratio = NNLO_SM + 4.743_dp * (Abs(delta_C7_0)**2 + Abs(delta_C7p_0)**2 ) &
        & + 0.789_dp *  (Abs(delta_C8_0)**2 + Abs(delta_C8p_0)**2 )          &
        & + Real( (-7.184_dp,0.612_dp) * delta_C7_0                          &
        &       + (-2.225_dp,-0.557_dp) * delta_C8_0                         &
        &       + (2.454_dp,-0.884_dp) * ( delta_C7_0 * Conjg(delta_C8_0)    &
        &                                +delta_C7p_0 * Conjg(delta_C8p_0) ) &
        &       , dp)  
  End If

  If (Present(A_CP)) Then
   A_CP = Aimag( ai_7(i_t)*(r7) + ai_8(i_t)*(r8)                       &
      &        + ai_87(i_t) * (r8*Conjg(r7)+r8p*Conjg(r7p))            &
      &        + ai_7eps(i_t) * r7*epsqC + ai_8eps(i_t) * r8*epsqC )   &
      & / ( a(i_t) + a_77(i_t) * ( Abs(r7)**2 + Abs(r7p)**2)           &
      &   + ar_7(i_t) * Real(r7,dp) + ar_8(i_t) * Real(r8,dp)          &
      &   + a_88(i_t) * ( Abs(r8)**2 + Abs(r8p)**2)                    &
      &   + ar_87(i_t) * Real(r8*Conjg(r7)+r8p*Conjg(r7p),dp) )
  End If

  Iname = Iname - 1

 End Subroutine B_to_Q_Gamma


 Subroutine BToSLL(c7, c7p, c8, c8p, c9, c9p, c10, c10p, BtoSEE, BtoSMuMu)
 !-----------------------------------------------------------------
 ! gives BR(b->s e+ e-) and BR(b->s e+ e-) using formulas 12/13 of
 ! T. Huber et al., hep-ph/0512066
 ! input: all running quantities are given at Q=160._dp
 !  Wilson coefficnets c7, c7p, c8, c8p, c9, c9p, c10, c10p
 ! output:
 !  BtoSEE ............... BR(b -> s e+ e-) 
 !  BtoSMuMu ............. BR(b -> s mu+ mu-)
 ! written by Werner Porod, 07.12.05
 ! 15.03.2014: taking now Wilson coefficients as input
 !----------------------------------------------------------------
 Implicit None
  !------------
  ! Wilson coefficients
  !------------
  Complex(dp), Intent(in) :: c7(2), c7p(2), c8(2), c8p(2), c9(3,2), c9p(3,2) &
    & , c10(3,2), c10p(3,2)
  !------------
  ! results
  !------------
  Real(dp), Intent(out) :: BtoSEE, BToSMuMu
  Complex(dp) :: r7, r7p, r8, r8p, r9(2), r9p(2), r10(2), r10p(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = "BToSLL"

  r7 = c7(1) / c7(2)
  r7p = c7p(1) / c7(2)
  r8 = c8(1) / c8(2)
  r8p = c8p(1) / c8(2)

  r9 = c9(1:2,1) / c9(1:2,2)
  r9p = c9p(1:2,1) / c9(1:2,2)
  r10 = c10(1:2,1) / c10(1:2,2)
  r10p = c10p(1:2,1) / c10(1:2,2)

  BtoSEE = (2.3148_dp - 1.658e-3_dp * Aimag(R10(1))                            &
       & + 5.e-4_dp * Aimag(r10(1)*Conjg(r8) + r10p(1)*Conjg(r8p) )            &
       & + 5.23e-2_dp * Aimag(r7) + 5.18e-3_dp * Aimag(r8)                     &
       & + 2.266e-2_dp * Aimag(r7 * Conjg(r8) + r7p * Conjg(r8p) )             &
       & + 4.96e-3_dp * Aimag(r7 * Conjg(r9(1)) + r7p * Conjg(r9p(1)) )        &
       & + 2.61e-2_dp * Aimag(r8 * Conjg(r9(1)) + r8p * Conjg(r9p(1)) )        &
       & - 6.21e-3_dp * Aimag(r9(1)) - 0.5420_dp * Real( r10(1), dp)           &
       & - 3.340e-2_dp * Real(r7,dp) - 1.35e-2_dp * Real(r8,dp)                &
       & + 1.53e-2_dp * Real(r7*Conjg(r10(1)) + r7p*Conjg(r10p(1)), dp )       &
       & + 6.73e-2_dp * Real(r7 * Conjg(r8) + r7p * Conjg(r8p), dp )           &
       & - 0.86916_dp * Real(r7*Conjg(r9(1)) + r7p*Conjg(r9p(1)), dp )         &
       & + 1.85e-3_dp * Real(r8*Conjg(r10(1)) + r8p*Conjg(r10p(1)), dp )       &
       & - 9.921e-2_dp * Real(r8*Conjg(r9(1)) + r8p*Conjg(r9p(1)), dp )        &
       & + 2.833_dp * Real(r9(1),dp) + 0.2804_dp * (Abs(r7)**2 + Abs(r7p)**2 ) &
       & - 0.10698_dp * Real( r9(1) * Conjg(r10(1))                            &
       &                    + r9p(1) * Conjg(r10p(1)), dp)                     &
       & + 11.0348_dp * (Abs(r10(1))**2 + Abs(r10p(1))**2 )                    &
       & + 1.527_dp * (Abs(r9(1))**2 + Abs(r9p(1))**2 )                        &
       & + 3.763e-3_dp * (Abs(r8)**2 + Abs(r8p)**2 ) ) * 1.e-7_dp

  BtoSMuMu = (2.1774_dp - 1.658e-3_dp * Aimag(R10(2))                          &
      & + 5.e-4_dp * Aimag(r10(2)*Conjg(r8) + r10p(2)*Conjg(r8p) )             &
      & + 5.34e-2_dp * Aimag(r7) + 5.27e-3_dp * Aimag(r8)                      &
      & + 2.266e-2_dp * Aimag(r7 * Conjg(r8) + r7p * Conjg(r8p) )              &
      & + 4.96e-3_dp * Aimag(r7 * Conjg(r9(2)) + r7p * Conjg(r9p(2)) )         &
      & + 2.61e-2_dp * Aimag(r8 * Conjg(r9(2)) + r8p * Conjg(r9p(2)) )         &
      & - 1.15e-2_dp * Aimag(r9(2)) - 0.5420_dp * Real( r10(2), dp)            &
      & + 2.08e-2_dp * Real(r7,dp) - 9.38e-3_dp * Real(r8,dp)                  &
      & + 1.53e-2_dp * Real(r7*Conjg(r10(2)) + r7p*Conjg(r10p(2)), dp )        &
      & + 6.848e-2_dp * Real(r7 * Conjg(r8) + r7p * Conjg(r8p), dp )           &
      & - 0.8545_dp * Real(r7*Conjg(r9(2)) + r7p*Conjg(r9p(2)), dp )           &
      & + 1.85e-3_dp * Real(r8*Conjg(r10(2)) + r8p*Conjg(r10p(2)), dp )        &
      & - 9.81e-2_dp * Real(r8*Conjg(r9(2)) + r8p*Conjg(r9p(2)), dp )          &
      & + 2.6917_dp * Real(r9(2),dp) + 0.2880_dp * (Abs(r7)**2 + Abs(r7p)**2 ) &
      & - 0.10698_dp * Real( r9(2) * Conjg(r10(2))                             &
      &                    + r9p(2) * Conjg(r10p(2)), dp)                      &
      & + 10.7652_dp * (Abs(r10(2))**2 + Abs(r10p(2))**2 )                     &
      & + 1.4884_dp * (Abs(r9(2))**2 + Abs(r9p(2))**2 )                        &
      & + 3.81e-3_dp * (Abs(r8)**2 + Abs(r8p)**2 ) ) * 1.e-7_dp

  Iname = Iname - 1

 End Subroutine BToSLL


  Subroutine B_To_SNuNu(CKM, C11, C11p, BtoSNuNu)
 !-----------------------------------------------------
 ! gives BR(b->s nu nu) using formulas 5.10 of
 ! Bobeth et al., NPB630 (2002) 87
 ! input: 
 !  CKM ............. CKM Matrix
 !  C11 ............. Wilson coefficients C_11
 !  C11p ............ Wilson coefficients C_11'
 ! output:
 !  BtoSNuNu ............. BR(b -> s nu nu) 
 ! written by Werner Porod, 10.10.07
 ! 16.03.2014: taking now Wilson coefficients as input
 !-----------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: CKM(3,3), c11(3), c11p(3)
  Real(dp), Intent(out) :: BToSNuNu

  Integer :: i1
 
  Iname = Iname + 1
  NameOfUnit(Iname) = "B_To_SNuNu"

  BToSNuNu = 0._dp
  Do i1=1,3
   BToSNuNu = BToSNuNu + Abs(C11(i1))**2 + Abs(C11p(i1))**2 &
            &          - 0.08_dp * Real(C11(i1)*Conjg(C11p(i1)),dp)
  End Do
  BToSNuNu = 5.39e-6_dp * Abs(CKM(3,2)/CKM(2,3))**2 * BToSNuNu

  Iname = Iname - 1

 End Subroutine B_To_SNuNu

 Subroutine Calculate_Wilson_Coeff_MSSM(mf_d, mf_u, mf_l, mW, mSpm2, mC, mSup2  &
       & , mSdown2, mSneut2, mSlepton2, mglu, mN, mS02, mP02, gauge, vevSM, CKM &
       & , c_uWd, c_CSQQp_L, c_CSQQp_R, c_CDSu_L, c_CDSu_R, c_DGSd_L, c_DGSd_R  &
       & , c_DNSd_L, c_DNSd_R, c_DDS0_L, c_DDS0_R, c_DDP0_L, c_DDP0_R, c_LLS0_L &
       & , c_LLS0_R, c_LLP0_L, c_LLP0_R, c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R &
       & , c_CNuSl_R, c_NuNSn_R, ZNN, ZUU, ZVV, ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR  &
       & , c_DDS0_1L_L, c_DDS0_1L_R, c_DDP0_1L_L, c_DDP0_1L_R, Q_in, i_in       &
       & , WC_c7, WC_c7p, WC_c8, WC_c8p, WC_c9, WC_c9p, WC_c10, WC_c10p         &
       & , WC_c11, WC_c11p)
 !------------------------------------------------
 ! calculates all Wilson coefficients if required
 ! the default is to use Q=m_Z, which however can be changed using the optional
 ! arguments Qin which gives the scale and i_in which defines the set to which 
 ! coefficients are saved
 ! input:
 !  - mf_d, mf_u, mW, mSpm2, mC, mSup2      &
 !      , mSdown2, mglu, mN, mS02, mP02, CKM, c_uWd, c_CSQQp_L           &
 !       , c_CSQQp_R, c_CDSu_L &
 !       , c_CDSu_R, c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R         &
 !       , c_DDS0_L, c_DDS0_R, c_DDP0_L, c_DDP0_R
 ! output:
 !  - C7(3,2), C7(3,1), C7'(3,2), C7'(3,1) 
 !  - C8(3,2), C8(3,1), C8'(3,2), C8'(3,1) 
 ! written by Werner Porod, 06.03.2014
 !------------------------------------------------
 Implicit None
  
  Real(dp), Intent(in) :: mf_d(3), mf_u(3), mf_l(3), mW, mSpm2(:), mC(:), mSup2(6)   &
               & , mSdown2(6), mglu, mN(:), mS02(:), mP02(:), mSneut2(:)    &
               & , mSlepton2(:), gauge(3), vevSM(2)
  Complex(dp), Intent(in) :: CKM(3,3), c_uWd(3,3), c_CSQQp_L(:,:,:)           &
    & , c_CSQQp_R(:,:,:), c_CDSu_L(:,:,:), c_CDSu_R(:,:,:), c_DGSd_R(:,:)     &
    & , c_DNSd_L(:,:,:), c_DGSd_L(:,:), c_DNSd_R(:,:,:), c_DDS0_L(:,:,:)      &
    & , c_DDS0_R(:,:,:), c_DDP0_L(:,:,:), c_DDP0_R(:,:,:), c_CLSn_L(:,:,:)    &
    & , c_CLSn_R(:,:,:), c_LNSl_L(:,:,:), c_LNSl_R(:,:,:), ZNN(:,:), ZUU(:,:) &
    & , ZVV(:,:), ZSdSdL(:,:), ZSdSdR(:,:), ZSuSuL(:,:), ZSuSuR(:,:)          &
    & , c_CNuSl_R(:,:,:), c_NuNSn_R(:,:,:), c_DDS0_1L_L(:,:,:)                &
    & , c_DDS0_1L_R(:,:,:), c_DDP0_1L_L(:,:,:), c_DDP0_1L_R(:,:,:)            &
    & , c_LLS0_L(:,:,:), c_LLS0_R(:,:,:), c_LLP0_L(:,:,:), c_LLP0_R(:,:,:)
  Real(dp), Intent(in) :: Q_in
  Integer, Intent(in), Optional :: i_in

  Complex(dp), Intent(inout) :: WC_c7(:,:,:,:), WC_c7p(:,:,:,:)                &
     & , WC_c8(:,:,:,:) , WC_c8p(:,:,:,:), WC_c9(:,:,:,:,:), WC_c9p(:,:,:,:,:) &
     & , WC_c10(:,:,:,:,:), WC_c10p(:,:,:,:,:), WC_c11(:,:,:,:,:)              &
     & , WC_c11p(:,:,:,:,:)

  Real(dp) :: mC2(Size(mC)), mN2(Size(mN)), mG2, a_s, sW2, tanb, e2
  Integer :: is, i1, i2, i3
  Complex(dp) :: c7(7), c7p(6), c8(7), c8p(6), norm, kappa_q
  Complex(dp) :: c9(3), c9p(3), c10(3), c10p(3), c9_c(3,6), c9p_c(3,6)     &
     & , c10_c(3,6), c10p_c(3,6), c11(3), c11p(3), c11_c(3,6), c11p_c(3,6) &
     & , B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2, coupLC &
     & , coupRC, fac

  If (Present(i_in)) Then
   is = i_in
  Else 
   is = 1
  End If

  Q_out(is) = Q_in
  WC_c7(is,:,:,:) = 0._dp
  WC_c7p(is,:,:,:) = 0._dp
  WC_c8(is,:,:,:) = 0._dp
  WC_c8p(is,:,:,:) = 0._dp
  WC_c9(is,:,:,:,:) = 0._dp
  WC_c9p(is,:,:,:,:) = 0._dp
  WC_c10(is,:,:,:,:) = 0._dp
  WC_c10p(is,:,:,:,:) = 0._dp
  WC_c11(is,:,:,:,:) = 0._dp
  WC_c11p(is,:,:,:,:) = 0._dp
  !---------------------------------------
  ! coefficients for b -> s gamma  
  !---------------------------------------

  Do i1=1,2
   norm = 0.5_dp * c_uWd(3,3) * Conjg(c_uWd(3,i1)) / mW2
   Call C_7(3, i1, mW, mf_u, c_uWd, mf_d, mSpm2, c_CSQQp_L, c_CSQQp_R &
        & , mC, mSup2, c_CDSu_L, c_CDSu_R, mGlu, mSdown2, c_DGSd_L    &
        & , c_DGSd_R, mN,  c_DNSd_L, c_DNSd_R, mS02, c_DDS0_L       &
        & , c_DDS0_r, mP02, c_DDP0_L, c_DDP0_R, c7 )
   Call C_7p(3, i1, mf_u, mf_d, mSpm2, c_CSQQp_L, c_CSQQp_R             &
        & , mC, mSup2, c_CDSu_L, c_CDSu_R, mGlu, mSdown2, c_DGSd_L    &
        & , c_DGSd_R, mN,  c_DNSd_L, c_DNSd_R, mS02, c_DDS0_L       &
        & , c_DDS0_r, mP02, c_DDP0_L, c_DDP0_R, c7p  )
   Call C_8(3, i1, mW, mf_u, c_uWd, mf_d, mSpm2, c_CSQQp_L, c_CSQQp_R &
        & , mC, mSup2, c_CDSu_L, c_CDSu_R, mGlu, mSdown2, c_DGSd_L    &
        & , c_DGSd_R, mN,  c_DNSd_L, c_DNSd_R, mS02, c_DDS0_L       &
        & , c_DDS0_r, mP02, c_DDP0_L, c_DDP0_R, c8)
   Call C_8p(3, i1, mf_u, mf_d, mSpm2, c_CSQQp_L, c_CSQQp_R             &
        & , mC, mSup2, c_CDSu_L, c_CDSu_R, mGlu, mSdown2, c_DGSd_L    &
        & , c_DGSd_R, mN,  c_DNSd_L, c_DNSd_R, mS02, c_DDS0_L       &
        & , c_DDS0_r, mP02, c_DDP0_L, c_DDP0_R, c8p )

   WC_c7(is,3,i1,:) = c7(1:2) / norm
   WC_c7p(is,3,i1,1) = c7p(1) / norm
   WC_c7p(is,3,i1,2) = 0._dp
   WC_c8(is,3,i1,:) = c8(1:2) / norm
   WC_c8p(is,3,i1,1) = c8p(1) / norm
   WC_c8p(is,3,i1,2) = 0._dp
  End Do

  mC2 = mC**2
  mN2 = mN**2
  mG2 = mGlu**2

  a_s = gauge(3)**2 * oo4pi
  sW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  tanb = vevSM(2) / vevSM(1)
  e2 = gauge(2)**2 * sW2

  Do i2=1,2
   kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,3) * Conjg( CKM(3,i2) )
   kappa_q = 1._dp / kappa_q
   Do i1=1,3
    Call C_9(3, i2, i1, Q_in, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2, mN2, 4 &
       & , mSneut2, mSlepton2, mSup2, mSdown2, tanb, ZUU, ZVV, ZNN, ZSdSdL    &
       & , ZSdSdR, ZSuSuL, ZSuSuR, c_CDSu_L, c_CDSu_R, c_CLSn_L         &
       & , c_CLSn_R, c_LNSl_L, c_LNSl_R, c_DGSd_L, c_DGSd_R         &
       & , c_DNSd_L, c_DNSd_R, a_s, kappa_q, e2, sW2, .False.             &
       & , c9(i1), c9p(i1), c9_c(i1,:), c9p_c(i1,:) )

    WC_C9(is,3,i2,i1,1) = c9(i1)
    WC_C9(is,3,i2,i1,2) = c9_c(i1,2)
    WC_C9p(is,3,i2,i1,1) = c9p(i1)
    WC_C9p(is,3,i2,i1,2) = c9p_c(i1,2)

    Call C_10_10p(3, i2, i1, Q_in, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2  &
     & , mN2, 4, mSneut2, mSlepton2, mSup2, mSdown2, tanb, ZUU, ZVV, ZNN    &
     & , ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR, c_CDSu_L, c_CDSu_R, c_CLSn_L &
     & , c_CLSn_R, c_LNSl_L, c_LNSl_R, c_DGSd_L, c_DGSd_R         &
     & , c_DNSd_L, c_DNSd_R, a_s, kappa_q, e2, sW2, .False.             &
     & , c10(i1), c10p(i1), c10_c(i1,:), c10p_c(i1,:))

    WC_C10(is,3,i2,i1,1) = c10(i1)
    WC_C10(is,3,i2,i1,2) = c10_c(i1,2)
    WC_C10p(is,3,i2,i1,1) = c10p(i1)
    WC_C10p(is,3,i2,i1,2) = c10p_c(i1,2)
   End Do
  End Do

  Do i1=1,3
   Do i2=2,3
    If (i2.Eq.3) kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,3) * Conjg( CKM(3,2) )
    If (i2.Eq.2) kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,2) * Conjg( CKM(3,1) )
    kappa_q = 1._dp / kappa_q
    Call C_11_11p(i2, i2-1, i1, Q_in, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2      &
     & , mN2, mSneut2, mSlepton2, mSup2, mSdown2, tanb &
     & , ZNN, ZUU, ZVV, ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR &
     & , c_CDSu_L, c_CDSu_R, c_CNuSl_R, c_NuNSn_R, c_DNSd_L         &
     & , c_DNSd_R, c_DGSd_L, c_DGSd_R, a_s, kappa_q, e2, sW2, .False.     &
     & , c11(i1), c11p(i1), c11_c(i1,:), c11p_c(i1,:))
    WC_C11(is,i2,i2-1,i1,1) = c11(i1)
    WC_C11(is,i2,i2-1,i1,2) = c11_c(i1,2)
    WC_C11p(is,i2,i2-1,i1,1) = c11p(i1)
    WC_C11p(is,i2,i2-1,i1,2) = c11p_c(i1,2)
   End Do
  End Do
 
  !-------------------------------------------------------------
  ! penguin for 2d2l operator
  !-------------------------------------------------------------
  WC_2d2l_CS(is,:,:,:) = 0._dp
  WC_2d2l_CSp(is,:,:,:) = 0._dp
  WC_2d2l_CP(is,:,:,:) = 0._dp
  WC_2d2l_CPp(is,:,:,:) = 0._dp

  Do i1=1,2
   fac = -sqrt2 * 4._dp * Pi2 /(G_F * gauge(2)**2 * CKM(3,i1) * Conjg(CKM(3,3)) )
   Do i2=1,3
    Do i3=1,Size(mS02)
     WC_2d2l_CS(is,3,i1,i2) = WC_2d2l_CS(is,3,i1,i2) &
      & + fac * c_DDS0_1L_L(3,i1,i3) * c_LLS0_L(i2,i2,i3) / (mS02(i3)*mf_d(3))
     WC_2d2l_CSp(is,3,i1,i2) = WC_2d2l_CSp(is,3,i1,i2) &
      & + fac * c_DDS0_1L_R(3,i1,i3) * c_LLS0_R(i2,i2,i3) / (mS02(i3)*mf_d(i1))
     WC_2d2l_CP(is,3,i1,i2) = WC_2d2l_CP(is,3,i1,i2) &
      & + fac * c_DDP0_1L_L(3,i1,i3) * c_LLP0_L(i2,i2,i3) / (mP02(i3)*mf_d(3))
     WC_2d2l_CPp(is,3,i1,i2) = WC_2d2l_CPp(is,3,i1,i2) &
      & + fac * c_DDP0_1L_R(3,i1,i3) * c_LLP0_R(i2,i2,i3) / (mP02(i3)*mf_d(i1))
    End Do
   End Do
  End Do
  !--------------------------------------------------------------------
  ! box diagrams for meson mixing
  ! note, that these are the B_i of appendix 4A of Buras et al
  ! NPB 659 (2003) 3 which differ by a factor (G_F M_W V_tb^* V_tq)^2
  ! from the Wilson coefficients
  !--------------------------------------------------------------------
  WC_4d_VLL(is,:,:) = 0._dp
  WC_4d_VRR(is,:,:) = 0._dp
  WC_4d_LR1(is,:,:) = 0._dp
  WC_4d_LR2(is,:,:) = 0._dp
  WC_4d_SLL1(is,:,:) = 0._dp
  WC_4d_SRR1(is,:,:) = 0._dp
  WC_4d_SLL2(is,:,:) = 0._dp
  WC_4d_SRR2(is,:,:) = 0._dp
  Do i1=2,3
   Do i2=1,i1-1
    Call Delta_F2_Boxes(i1, i2, -0.5_dp, mf_u, mf_d, mC, mC2, mN, mN2  &
       & , mGlu, mSpm2, mSup2, mSdown2, c_uWd, c_CSQQp_L, c_CSQQp_R    &
       & , c_CDSu_L, c_CDSu_R, c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R  &
       & , B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2  )
    !---------------------------------------------
    ! adding now the double penguin contributions
    !---------------------------------------------
    Do i3=1,Size(mS02)
     coupLC = c_DDS0_1L_L(i1,i2,i3)
     coupRC = c_DDS0_1L_R(i1,i2,i3)
     B_LR2 = B_LR2 + 16._dp * Pi2 * coupLC * coupRC / mS02(i3) 
     B_SLL1 = B_SLL1 + 8._dp * Pi2 * coupLC**2 / mS02(i3) 
     B_SRR1 = B_SRR1 + 8._dp * Pi2 * coupRC**2 / mS02(i3)
    End Do

    Do i3=2,Size(mP02)
     coupLC = c_DDP0_1L_L(i1,i2,i3)
     coupRC = c_DDP0_1L_R(i1,i2,i3)
     B_LR2 = B_LR2 + 16._dp * Pi2  * coupLC * coupRC / mP02(i3) 
     B_SLL1 = B_SLL1 + 8._dp * Pi2 * coupLC**2 / mP02(i3) 
     B_SRR1 = B_SRR1 + 8._dp * Pi2 * coupRC**2 / mP02(i3)
    End Do
    WC_4d_VLL(is,i1,i2) = B_VLL
    WC_4d_VRR(is,i1,i2) = B_VRR
    WC_4d_LR1(is,i1,i2) = B_LR1
    WC_4d_LR2(is,i1,i2) = B_LR2
    WC_4d_SLL1(is,i1,i2) = B_SLL1
    WC_4d_SRR1(is,i1,i2) = B_SRR1
    WC_4d_SLL2(is,i1,i2) = B_SLL2
    WC_4d_SRR2(is,i1,i2) = B_SRR2
   End Do
  End Do

 End Subroutine Calculate_Wilson_Coeff_MSSM


 Subroutine C_7(i, j, mW, mf_u, cpl_uWd  &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R         &
        & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res   )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  e m_i / (16 Pi^2) \bar{d}_L,j,a \sigma^\mu\nu d_R,i,a F_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mW, mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_uWd(3,3), cpl_CSQQp_L(:,:,:)         &
     & , cpl_CSQQp_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)     &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(7)
  Integer :: i1, i2, n_i
  Real(dp) :: xt
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_7"

  res = 0._dp
  !----------------
  ! standard model
  !----------------
  Do i1=1,3
   xt = (mf_u(i1)/mW)**2
   res(2) = res(2) - Conjg( cpl_uWd(i1, j ) ) * cpl_uWd(i1, i)   &
       &             * xt * (2._dp * F1(xt) + 3._dp * F2(xt) )
  End Do
  res(2) = 0.25_dp * res(2) / mW**2

  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(3) = res(3) &
      & - ( Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &     * (e_u * F1(xt)   + F2(xt) )                          &
      &   + Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &     * mf_u(i2) * (e_u * F3(xt)  + F4(xt) ) / mf_d(i) )   / mSpm2(i1)
   End Do
  End Do
  res(3) = 0.25_dp * res(3) 
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(4) = res(4) &
      & + ( Conjg(cpl_CDSu_R(i1,i,i2)) * cpl_CDSu_R(i1,j,i2)    &
      &    * (e_u * F2(xt)  + F1(xt) )                          &
      &   + Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2)  &
      &      * mC(i1) * (e_u * F4(xt)  + F3(xt) ) / mf_d(i)  )  / mSup2(i2)
   End Do
  End Do
  res(4) = 0.25_dp * res(4) 

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(5) = res(5) + ( Conjg(cpl_DGSd_R(i,i2)) * cpl_DGSd_R(j,i2) * F2(xt) &
      &               + Conjg(cpl_DGSd_L(i,i2)) * cpl_DGSd_R(j,i2)          &
      &                 * mGlu * F4(xt) / mf_d(i) )   / mSdown2(i2)
   End Do
   res(5) = e_d * 4._dp * res(5) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(6) = res(6) &
      & + ( Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2) * F2(xt)   &
      &   + Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2)            &
      &         * mN(i1) * F4(xt) / mf_d(i) )   / mSdown2(i2)
    End Do
   End Do
   res(6) = 0.25_dp * e_d * res(6) 
   !----------
   ! Higgs
   !----------
   n_i = Size( mS02 )
   Do i1=1,3
    Do i2=1,n_i
     xt = mf_d(i1)**2/mS02(i2)
     res(7) = res(7) &
       & + ( Conjg(cpl_DDS0_R(i,i1,i2)) * cpl_DDS0_R(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mS02(i2)
    End Do
    Do i2=2,n_i
     xt = mf_d(i1)**2/mP02(i2)
     res(7) = res(7) &
       & + ( Conjg(cpl_DDP0_R(i1,i,i2)) * cpl_DDP0_R(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDP0_L(i1,i,i2) ) * cpl_DDP0_R(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mP02(i2)
    End Do
   End Do
   res(7) = 0.25_dp * e_d * res(7) 
  End If

  res(1) = Sum( res(2:7) )
  
  Iname = Iname - 1

 End Subroutine C_7

 Subroutine C_7p(i, j, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res  )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  e m_i / (16 Pi^2) \bar{d}_R,j,a \sigma^\mu\nu d_L,i,a F_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_CSQQp_L(:,:,:), cpl_CSQQp_R(:,:,:)   &
      & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)                        &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(6)

  Integer :: i1, i2, n_i
  Real(dp) :: xt

  Iname = Iname + 1
  NameOfUnit(Iname) = "C_7p"

  !----------------
  ! standard model
  !----------------
  res = 0._dp
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(2) = res(2) &
      & - ( Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)    &
      &          * (e_u * F1(xt)  + F2(xt) )                        &
      &   + Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)    &
      &          * mf_u(i2) * (e_u * F3(xt)  + F4(xt) ) / mf_d(i) ) &
      &  / mSpm2(i1)
   End Do
  End Do
  res(2) = 0.25_dp * res(2)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(3) = res(3) &
      & + ( Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)    &
      &                * (e_u * F2(xt)  + F1(xt) )                &
      &   + Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)    &
      &         * mC(i1) * (e_u * F4(xt)  + F3(xt) ) / mf_d(i))   &
      &  / mSup2(i2)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_L(j,i2) * F2(xt)  &
      &   + Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_L(j,i2)           &
      &     * mGlu * F4(xt) / mf_d(i) )    / mSdown2(i2)
   End Do
   res(4) = e_d * 4._dp * res(4) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(5) = res(5) &
      & + ( Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2)           &
      &         * mN(i1) * F4(xt) / mf_d(i)  )     / mSdown2(i2)
    End Do
   End Do
   res(5) = 0.25_dp * e_d * res(5) 
   !----------
   ! Higgs
   !----------
   n_i = Size( mS02 )
   Do i1=1,3
    Do i2=1,n_i
     xt = mf_d(i1)**2/mS02(i2)
     res(6) = res(6) &
       & + ( Conjg(cpl_DDS0_L(i1,i,i2)) * cpl_DDS0_L(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDS0_R(i1,i,i2) ) * cpl_DDS0_L(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mS02(i2)
    End Do
    Do i2=2,n_i
     xt = mf_d(i1)**2/mP02(i2)
     res(6) = res(6) &
       & + ( Conjg(cpl_DDP0_L(i,i1,i2)) * cpl_DDP0_L(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDP0_R(i,i1,i2) ) * cpl_DDP0_L(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mP02(i2)
    End Do
   End Do
   res(6) = 0.25_dp * e_d * res(6) 
  End If

  res(1) = Sum( res(2:6) )
  
  Iname = Iname - 1

 End Subroutine C_7p

 Subroutine C_8(i, j, mW, mf_u, cpl_uWd  &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R   &
        & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res)
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  g_s m_i / (16 Pi^2) \bar{d}_L,j,a \sigma^\mu\nu d_R,i,b t_c^ab G^c_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mW, mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_uWd(3,3), cpl_CSQQp_L(:,:,:)         &
     & , cpl_CSQQp_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)     &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(7)

  Integer :: i1, i2, n_i
  Real(dp) :: xt
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_8"

  res = 0._dp
  !----------------------------
  ! standard model
  !----------------------------
  Do i1=1,3
   xt = (mf_u(i1)/mW)**2
   res(2) = res(2) - Conjg( cpl_uWd(i1, j ) ) * cpl_uWd(i1, i) * xt * F1(xt)
  End Do
  res(2) = 0.75_dp * res(2) / mW**2
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(3) = res(3) &
      & - ( Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2) * F1(xt)  &
      &   + Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &                *  mf_u(i2) * F3(xt) / mf_d(i) )      / mSpm2(i1)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2) *  F2(xt)    &
      &   + Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2)              &
      &                          * mC(i1) * F4(xt) / mf_d(i))    / mSup2(i2)
   End Do
  End Do
  res(4) = 0.25_dp * res(4)

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(5) = res(5) &
      & + ( Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_R(j,i2)       &
      &                * (9._dp * F1(xt)  + F2(xt) )           &
      &      + Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_R(j,i2)    &
      &         * mGlu * (9._dp * F3(xt) + F4(xt) ) / mf_d(i)) &
      &  / mSdown2(i2)
   End Do
   res(5) = 0.5_dp * res(5) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(6) = res(6) &
      & + ( Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2)           &
      &         * mN(i1) * F4(xt) / mf_d(i) )   / mSdown2(i2)
    End Do
   End Do
   res(6) = 0.25_dp * res(6)
   !----------
   ! Higgs
   !----------
   n_i = Size( mS02 )
   Do i1=1,3
    Do i2=1,n_i
     xt = mf_d(i1)**2/mS02(i2)
     res(7) = res(7) &
       & - ( Conjg( cpl_DDS0_R(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mS02(i2)
    End Do
    Do i2=2,n_i
     xt = mf_d(i1)**2/mP02(i2)
     res(7) = res(7) &
       & - ( Conjg( cpl_DDP0_R(i,i1,i2) ) * cpl_DDP0_R(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDP0_L(i,i1,i2) ) * cpl_DDP0_R(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mP02(i2)
    End Do
   End Do
   res(7) = 0.25_dp * res(7)
  End If

  res(1) = Sum(res(2:7))

  Iname = Iname - 1

 End Subroutine C_8

 Subroutine C_8p(i, j, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  g_s m_i / (16 Pi^2) \bar{d}_R,j,a \sigma^\mu\nu d_L,i,b t_c^ab G^c_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_CSQQp_L(:,:,:), cpl_CSQQp_R(:,:,:)   &
      & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)                        &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(6)

  Integer :: i1, i2, n_i
  Real(dp) :: xt

  Iname = Iname + 1
  NameOfUnit(Iname) = "C_8p"

  !----------------------------
  ! standard model, tree level
  !----------------------------
  res = 0._dp
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(2) = res(2) &
      & - ( Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2) * F1(xt)   &
      &   + Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)            &
      &                * mf_u(i2) * F3(xt) / mf_d(i) ) / mSpm2(i1)
   End Do
  End Do
  res(2) = 0.25_dp * res(2)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(3) = res(3) &
      & + ( Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2) * F2(xt)    &
      &   + Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)             &
      &                   * mC(i1) * F4(xt) / mf_d(i)    )     / mSup2(i2)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)

  If (GenerationMixing) Then
   !--------
   ! gluino
   !--------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_L(j,i2)         &
      &                * (9._dp * F1(xt)  + F2(xt) )             &
      &   + Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_L(j,i2)         &
      &         * mGlu * (9._dp * F3(xt) + F4(xt) ) / mf_d(i))   &
      &  / mSdown2(i2)
   End Do
   res(4) = 0.5_dp * res(4) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(5) = res(5) &
      & + ( Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2)           &
      &                      * mN(i1) * F4(xt) / mf_d(i)     )   / mSdown2(i2)
    End Do
   End Do
   res(5) = 0.25_dp * res(5)
   !----------
   ! Higgs
   !----------
   n_i = Size( mS02 )
   Do i1=1,3
    Do i2=1,n_i
     xt = mf_d(i1)**2/mS02(i2)
     res(6) = res(6) &
       & - ( Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_L(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDS0_R(i,i1,i2) ) * cpl_DDS0_L(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mS02(i2)
    End Do
    Do i2=2,n_i
     xt = mf_d(i1)**2/mP02(i2)
     res(6) = res(6) &
       & - ( Conjg( cpl_DDP0_L(i,i1,i2) ) * cpl_DDP0_L(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDP0_R(i,i1,i2) ) * cpl_DDP0_L(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mP02(i2)
    End Do
   End Do
   res(6) = 0.25_dp * res(6)
  End If

  res(1) = Sum( res(2:6) )
  
  Iname = Iname - 1

 End Subroutine C_8p

 Subroutine C_9(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2, mC2, mN2     &
       & , n_n, mSnu2, mSlept2, mSqu2, mSqd2, tanb                            &
       & , ZUU, ZVV, ZNN, ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR, c_CDSu_L, c_CDSu_R  &
       & , c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R, c_DGSd_L, c_DGSd_R         &
       & , c_DNSd_L, c_DNSd_R, a_s, kappa_q, e2, sW2, l_QCD                   &
       & , c9, c9p, c9_c, c9p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_A and c_A' of the effective operator
 !
 !             2 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c_9 O_A + c_9' O_A')
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! with O_A = (\bar(q) \gamma^\mu P_L b) (\bar(l) \gamma_\mu l)
 !      O_A' = (\bar(q) \gamma^\mu P_R b) (\bar(l) \gamma_\mu l)
 ! for the transitions b -> q l+ l- with q=s,d;  s -> d mu+ mu-;
 ! B_q -> l+ l- and K_L -> mu+ mu- are calculated.
 ! The formulas include QCD corrections and are based on
 ! C.bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 14.12.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of the lepton l
  Integer, Intent(in) :: n_n               ! number of neutralinos
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(:)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  ! reduced couplings of Z to a pair squarks, neutralinos and charginos,
  ! respectivly
  Complex(dp), Intent(in) :: ZSuSuL(6,6), ZSuSuR(6,6), ZSdSdL(6,6), ZSdSdR(6,6)
  Complex(dp), Intent(in) :: ZNN(:,:), ZUU(:,:), ZVV(:,:)
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(:,:,:), c_DNSd_R(:,:,:)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - lepton - sneutrino
  Complex(dp), Intent(in) ::  c_CLSn_L(2,3,3), c_CLSn_R(2,3,3)
  ! coupling neutralino - lepton - slepton
  Complex(dp), Intent(in) ::  c_LNSl_L(:,:,:), c_LNSl_R(:,:,:)
  Complex(dp), Intent(in) :: kappa_Q    ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Real(dp), Intent(in) :: e2            ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2           ! sin^2(theta_W)
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included

  Complex(dp), Intent(out) :: C9, C9p ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C9_c(6), C9p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, x_ij(2,2), y_ai(6,2), v_fi(3,2), r_a(6), s_ai(6,n_n) &
      & , n_bi(6,n_n), n_ij(n_n,n_n), Lt, Lua(6), a_sp, x_in
  Complex(dp) :: B_L(4), B_R(4), C_L(5), C_R(5), D_L(5), D_R(5), fact, fac(2)

  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_9"

  c9 = 0._dp
  c9p = 0._dp
  B_L = 0._dp
  B_R = 0._dp
  C_L = 0._dp
  C_R = 0._dp
  D_L = 0._dp
  D_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2

  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
   End Do
   Do i2=1,3
    v_fi(i2,i1) = mSnu2(i2) / mC2(i1)
   End Do
  End Do
  Do i1=1,n_n
   n_ij(i1,i1) = 1._dp
   Do i2=i1+1,n_n
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
   End Do
   Do i2=1,6
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    n_bi(i2,i1) = mSlept2(i2) / mN2(i1)
   End Do
  End Do

  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = 0.25_dp * (f_1_0(x) + a_sp * (f_1_1(x)+ 8._dp* x * Df_1_0(x) * Lt)) 
   B_L(1) = 0.25_dp * (f_2_0(x) + a_sp * (f_10_1(x) + 8._dp* x *Df_2_0(x) * Lt))
  Else
   C_L(1) = 0.25_dp * f_1_0(x)
   B_L(1) = 0.25_dp * f_2_0(x)
  End If 
  D_L(1) = h_0_1(x)

  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = - 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
         &       / (16._dp * mHp2 * mW2)

  D_L(2) = y * (2._dp/3._dp * h_0_5(y) - h_0_6(y)) / tanb**2
  D_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * mf_l(k)**2  / (16._dp * mHp2**2) &
         & * y * (2._dp/3._dp * h_0_5(y) - h_0_6(y))

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     C_L(3) = C_L(3) - c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))     &
         &             * ( fac(1) * ZUU(i1,i2) - fac(2) * ZVV(i2,i1) )
     C_R(3) = C_R(3) - c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))     &
          &             * ( fac(1) * ZVV(i1,i2) - fac(2) * ZUU(i2,i1) )
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      End If
      B_L(3) = B_L(3) + Conjg(c_CDSu_R(i2,i,i3))* c_CDSu_R(i1,j,i3) &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_R(i1,k,i4))          &
          &                    * c_CLSn_R(i2,k,i4)                          &
          &           - fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_L(i2,k,i4)       &
          &                    * Conjg(c_CLSn_L(i1,k,i4))    ) / mC2(i2)
      B_R(3) = B_R(3) - Conjg(c_CDSu_L(i2,i,i3))*c_CDSu_L(i1,j,i3) &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_L(i1,k,i4))           &
          &                    * c_CLSn_L(i2,k,i4)                           &
          &           - fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_R(i2,k,i4)        &
          &                    * Conjg(c_CLSn_R(i1,k,i4))    ) / mC2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    x_in = 1._dp / y_ai(i2,i1)
    fac(1) = (2._dp/3._dp * h_0_6(x_in) - h_0_5(x_in)) / mSqu2(i2)
    D_L(3) = D_L(3) + fac(1) * Conjg(c_CDSu_R(i1,i,i2))* c_CDSu_R(i1,j,i2)
    D_R(3) = D_R(3) + fac(1) * Conjg(c_CDSu_L(i1,i,i2))* c_CDSu_L(i1,j,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) - fac(1) * ZSuSuL(i3,i2)   &
         &               * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) - fac(1) * ZSuSuR(i3,i2)   &
         &               * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
  End Do
  B_L(3) = B_L(3) * Kappa_Q * sW2
  B_R(3) = B_R(3) * Kappa_Q * sW2
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = C_L(3) * fact
  C_R(3) = C_R(3) * fact
  fact = 16._dp * fact * mW2
  D_L(3) = D_L(3) * fact
  D_R(3) = D_R(3) * fact

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,n_n
   Do i2=1,n_n
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If

     C_L(4) = C_L(4) - (fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1) )&
         &            * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) - (fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1) )&
         &            * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3)) 

     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) - Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * (- Conjg(c_LNSl_R(k,i1,i4)) * c_LNSl_R(k,i2,i4)   &
          &                 + Conjg(c_LNSl_L(k,i2,i4)) * c_LNSl_L(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_L(k,i2,i4) * Conjg(c_LNSl_L(k,i1,i4))   &
          &                 - c_LNSl_R(k,i1,i4) * Conjg(c_LNSl_R(k,i2,i4)) ) &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) + Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_L(k,i1,i4)) * c_LNSl_L(k,i2,i4)   &
          &                 + Conjg(c_LNSl_R(k,i2,i4)) * c_LNSl_R(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_R(k,i2,i4) * Conjg(c_LNSl_R(k,i1,i4))   &
          &                 - c_LNSl_L(k,i1,i4) * Conjg(c_LNSl_L(k,i2,i4)) ) &
          &           ) / mN2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    x_in = 1._dp/s_ai(i2,i1)
    fac(1) = h_0_6(x_in) / mSqd2(i2)
    D_L(4) = D_L(4) + fac(1) * Conjg(c_DNSd_R(i,i1,i2))* c_DNSd_R(j,i1,i2)
    D_R(4) = D_R(4) + fac(1) * Conjg(c_DNSd_L(i,i1,i2))* c_DNSd_L(j,i1,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) - fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) - fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
    
  End Do
  B_L(4) = B_L(4) * Kappa_Q * sW2
  B_R(4) = B_R(4) * Kappa_Q * sW2
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = C_L(4) * fact
  C_R(4) = C_R(4) * fact
  fact = - 4._dp * kappa_Q * e2 / 3._dp
  D_L(4) = D_L(4) * fact
  D_R(4) = D_R(4) * fact

  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
   Do i2=1,6
    x_in = 1._dp/r_a(i2)
    fac(1) = h_0_6(x_in) / mSqd2(i2)
    D_L(5) = D_L(5) + fac(1) * Conjg(c_DGSd_R(i,i2))* c_DGSd_R(j,i2)
    D_R(5) = D_R(5) + fac(1) * Conjg(c_DGSd_L(i,i2))* c_DGSd_L(j,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) - fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) - fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)
  C_L(5) = C_L(5) * fact
  C_R(5) = C_R(5) * fact
  fact = - 64._dp * kappa_Q * e2 / 9._dp
  D_L(5) = D_L(5) * fact
  D_R(5) = D_R(5) * fact

  C9 = ((1._dp - 4._dp * sW2) *Sum(C_L) - Sum(B_L)) / sW2 - Sum(D_L)
  C9p = ((1._dp - 4._dp * sW2) * Sum(C_R) - Sum(B_R)) / sW2 - Sum(D_R)
  !---------------------
  ! c-quark contribution
  !---------------------

  C9 = c9 + 38._dp/27._dp - 4._dp * Log(mu**2/mW2)/ 9._dp

  If (Present(c9_c)) Then
   c9_c(1) = c9
   c9_c(2:6) = ((1._dp - 4._dp * sW2) * C_L ) / sW2 - D_L
   c9_c(2) = c9_c(2) + 38._dp/27._dp - 4._dp * Log(mu**2/mW2)/ 9._dp
   c9_c(2:5) =c9_c(2:5) - B_L / sW2
  End If
  If (Present(c9p_c)) Then
   c9p_c(1) = c9p
   c9p_c(2:6) = ((1._dp - 4._dp * sW2) * C_R ) / sW2 - D_R
   c9p_c(2:5) =c9p_c(2:5) - B_R / sW2
  End If

  Iname = Iname - 1

 End Subroutine C_9


 Subroutine C_QdQdLL_AAp(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2     &
       & , mC2, mN2, n_n, mSnu2, mSlept2, mSqu2, mSqd2, tanb                 &
       & , ZUU, ZVV, ZNN, ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR, c_CDSu_L, c_CDSu_R &
       & , c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R, c_DGSd_L, c_DGSd_R        &
       & , c_DNSd_L, c_DNSd_R, a_s, kappa_q, e2, sW2, l_QCD                  &
       & , c_A, c_Ap, c10_c, c10p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_A and c_A' of the effective operator
 !
 !             2 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c_A O_A + c_A' O_A')
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! with O_A = (\bar(q) \gamma^\mu P_L b) (\bar(l) \gamma_\mu gamma_5 l)
 !      O_A' = (\bar(q) \gamma^\mu P_R b) (\bar(l) \gamma_\mu gamma_5 l)
 ! for the transitions b -> q l+ l- with q=s,d;  s -> d mu+ mu-;
 ! B_q -> l+ l- and K_L -> mu+ mu- are calculated.
 ! The formulas include QCD corrections and are based on
 ! C.Bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 14.12.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of the lepton l
  Integer, Intent(in) :: n_n               ! number of neutralinos
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(:)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  ! reduced couplings of Z to a pair squarks, neutralinos and charginos,
  ! respectivly
  Complex(dp), Intent(in) :: ZSuSuL(6,6), ZSuSuR(6,6), ZSdSdL(6,6), ZSdSdR(6,6)
  Complex(dp), Intent(in) :: ZNN(:,:), ZUU(:,:), ZVV(:,:)
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - lepton - sneutrino
  Complex(dp), Intent(in) ::  c_CLSn_L(2,3,3), c_CLSn_R(2,3,3)
  ! coupling neutralino - lepton - slepton
  Complex(dp), Intent(in) ::  c_LNSl_L(3,4,6), c_LNSl_R(3,4,6)
  Complex(dp), Intent(in) :: kappa_Q    ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Real(dp), Intent(in) :: e2            ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2           ! sin^2(theta_W)
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included

  Complex(dp), Intent(out) :: C_A, C_Ap ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C10_c(6),C10p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, x_ij(2,2), y_ai(6,2), v_fi(3,2), r_a(6), s_ai(6,n_n) &
      & , n_bi(6,n_n), n_ij(n_n,n_n), Lt, Lua(6), a_sp
  Complex(dp) :: fact, fac(2), c_L(5), c_R(5), B_L(4), B_R(4)

  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_QdQdLL_AAp"

  C_A = 0._dp
  C_Ap = 0._dp
  B_L = 0._dp
  B_R = 0._dp
  C_L = 0._dp
  C_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2
  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
   End Do
   Do i2=1,3
    v_fi(i2,i1) = mSnu2(i2) / mC2(i1)
   End Do
  End Do
  Do i1=1,n_n
   n_ij(i1,i1) = 1._dp
   Do i2=i1+1,n_n
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
   End Do
   Do i2=1,6
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    n_bi(i2,i1) = mSlept2(i2) / mN2(i1)
   End Do
  End Do

  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = -0.25_dp *(f_1_0(x) + a_sp * (f_1_1(x)+ 8._dp* x * Df_1_0(x) * Lt))
   B_L(1) = 0.25_dp *(f_2_0(x) + a_sp * (f_10_1(x)+ 8._dp* x * Df_2_0(x) * Lt))
  Else
   C_L(1) = - 0.25_dp * f_1_0(x)
   B_L(1) = 0.25_dp * f_2_0(x)
  End If
                         
  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = - 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
      &              / (16._dp * mHp2 * mW2)

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     C_L(3) = C_L(3) + c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))     &
         &             * ( fac(1) * ZUU(i1,i2) - fac(2) * ZVV(i2,i1) )
     C_R(3) = C_R(3) + c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))     &
          &             * ( fac(1) * ZVV(i1,i2) - fac(2) * ZUU(i2,i1) )
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      End If
      B_L(3) = B_L(3) + Conjg(c_CDSu_R(i2,i,i3))* c_CDSu_R(i1,j,i3)      &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_R(i1,k,i4))       &
          &                    * c_CLSn_R(i2,k,i4)                       &
          &           + fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_L(i2,k,i4)    &
          &                    * Conjg(c_CLSn_L(i1,k,i4))    ) / mC2(i2)
      B_R(3) = B_R(3) - Conjg(c_CDSu_L(i2,i,i3))*c_CDSu_L(i1,j,i3)       &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_L(i1,k,i4))       &
          &                    * c_CLSn_L(i2,k,i4)                       &
          &           + fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_R(i2,k,i4)    &
          &                    * Conjg(c_CLSn_R(i1,k,i4))    ) / mC2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) + fac(1) * ZSuSuL(i3,i2)   &
         &               * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) + fac(1) * ZSuSuR(i3,i2)   &
         &               * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
    
  End Do
  B_L(3) = kappa_Q *sW2 * B_L(3)
  B_R(3) = kappa_Q *sW2 * B_R(3)
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = fact * C_L(3)
  C_R(3) = fact * C_R(3)

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,n_n
   Do i2=1,n_n
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If
     fac(1) = fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1)
     C_L(4) = C_L(4) + fac(1) * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) + fac(1) * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3)) 
     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) + Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_R(k,i1,i4)) * c_LNSl_R(k,i2,i4)   &
          &                 + Conjg(c_LNSl_L(k,i2,i4)) * c_LNSl_L(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_L(k,i2,i4) * Conjg(c_LNSl_L(k,i1,i4))   &
          &                 + c_LNSl_R(k,i1,i4) * Conjg(c_LNSl_R(k,i2,i4)) ) &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) - Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_L(k,i1,i4)) * c_LNSl_L(k,i2,i4)   &
          &                 + Conjg(c_LNSl_R(k,i2,i4)) * c_LNSl_R(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_R(k,i2,i4) * Conjg(c_LNSl_R(k,i1,i4))   &
          &                 + c_LNSl_L(k,i1,i4) * Conjg(c_LNSl_L(k,i2,i4)) ) &
          &           ) / mN2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) + fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) + fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
    
  End Do

  B_L(4) = kappa_Q *sW2 * B_L(4)
  B_R(4) = kappa_Q *sW2 * B_R(4)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = fact * C_L(4)
  C_R(4) = fact * C_R(4)
  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)

   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) + fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) + fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
   c_L(5) = c_L(5) * fact
   c_R(5) = c_R(5) * fact

  c_a = (Sum(C_L) + Sum(B_L)) / sW2
  c_ap= (Sum(C_R) + Sum(B_R)) / sW2

  If (Present(c10_c)) Then
   c10_c(1) = c_a
   c10_c(2:6) = C_L / sW2
   c10_c(2:5) =c10_c(2:5) + B_L / sW2
  End If
  If (Present(c10p_c)) Then
   c10p_c(1) = c_ap
   c10p_c(2:6) = C_R / sW2
   c10p_c(2:5) =c10p_c(2:5) + B_R / sW2
  End If

  Iname = Iname - 1

 End Subroutine C_QdQdLL_AAp


 Subroutine C_QdQdNuNu_LR(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2   &
          & , mC2, mN2, mSnu2, mSlept2, mSqu2, mSqd2, tanb, ZNN, ZUU, ZVV   &
          & , ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR, c_CDSu_L, c_CDSu_R, c_CNuSl_R &
          & , c_NuNSn_R, c_DNSd_L, c_DNSd_R, c_DGSd_L, c_DGSd_R, a_s        &
          & , kappa_q, e2, sW2, l_QCD, c11, c11p, C11_c, C11p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_L and c_R of the effective operator
 !
 !             4 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c11 O_L + c11p O_R)
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! for the transitions b -> q nu nu with q=s,d and s -> d nu nu are
 ! calculated. 
 ! The formulas include QCD corrections and are based on
 ! C.Bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 10.11.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of neutrinos
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(4)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  Real(dp), Intent(in) :: e2               ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2              ! sin^2(theta_W)
  ! reduced couplings of Z to a pair squarks, neutralinos and charginos,
  ! respectivly
  Complex(dp), Intent(in) :: ZSuSuL(6,6), ZSuSuR(6,6), ZSdSdL(6,6), ZSdSdR(6,6)
  Complex(dp), Intent(in) :: ZNN(:,:), ZUU(:,:), ZVV(:,:)
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - neutrino - slepton
  Complex(dp), Intent(in) ::  c_CNuSl_R(2,3,6)
  ! coupling neutralino - neutrino - sneutrino
  Complex(dp), Intent(in) ::  c_NuNSn_R(3,4,3)
  Complex(dp), Intent(in) :: kappa_Q     ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included

  Complex(dp), Intent(out) :: c11, c11p ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C11_c(6),C11p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, x_ij(2,2), y_ai(6,2), z_bi(6,2), r_a(6), s_ai(6,4) &
      & , n_bi(3,4), n_ij(4,4), Lt, Lua(6), a_sp
  Complex(dp) :: fact, c_L(5), c_R(5), b_L(5), b_R(5), fac(2)
  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_QdQdNuNu_LR"

  c_L = 0._dp
  c_R = 0._dp
  b_L = 0._dp
  b_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2
  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
    z_bi(i2,i1) = mSlept2(i2)**2 / mC2(i1)
   End Do
  End Do
  Do i1=1,4
   n_ij(i1,i1) = 1._dp
   Do i2=i1+1,4
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
   End Do
   Do i2=1,3
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    s_ai(i2+3,i1) = mSqd2(i2+3) / mN2(i1)
    n_bi(i2,i1) = mSnu2(i2) / mN2(i1)
   End Do
  End Do
  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = 0.25_dp * (f_1_0(x) + a_sp* (f_1_1(x) + 8._dp* x* Df_1_0(x) * Lt) )
   B_L(1) = - f_2_0(x) - a_sp * (f_6_1(x) + 8._dp * x * Df_2_0(x) * Lt)
  Else
   C_L(1) = 0.25_dp * f_1_0(x)
   B_L(1) = - f_2_0(x)
  End If
  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = - 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
         &           / (16._dp * mHp2 * mW2)

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     ! Z contribution
     C_L(3) = C_L(3) - c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))   &
         &             * ( fac(1) * ZUU(i1,i2) - fac(2) * ZVV(i2,i1) )
     C_R(3) = C_R(3) - c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))   &
          &             * ( fac(1) * ZVV(i1,i2) - fac(2) * ZUU(i2,i1) )
     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
      End If
      B_L(3) = B_L(3) + fac(1)                    &
          &       * Conjg(c_CNuSl_R(i1,k,i4) * c_CDSu_R(i2,i,i3) )   &
          &       * c_CNuSl_R(i2,k,i4) * c_CDSu_R(i1,j,i3) / mC2(i2)
      B_R(3) = B_R(3) - fac(2) * Sqrt(x_ij(i1,i2))   &
          &       * Conjg(c_CNuSl_R(i1,k,i4) * c_CDSu_L(i2,i,i3) )     &
          &       * c_CNuSl_R(i2,k,i4) * c_CDSu_L(i1,j,i3) / mC2(i2)
     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) - fac(1) * ZSuSuL(i3,i2)    &
         &                    * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) - fac(1) * ZSuSuR(i3,i2)    &
         &                    * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
       
  End Do

  B_L(3) = 0.5_dp * kappa_Q * sW2 * B_L(3)
  B_R(3) = kappa_Q * sW2 * B_R(3)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = fact * C_L(3)
  C_R(3) = fact * C_R(3)

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If
     fac(1) = - fac(1) * ZNN(i1,i2) - fac(2) * ZNN(i2,i1)
     C_L(4) = C_L(4) + fac(1) * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) + fac(1) * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3))
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = 0.5_dp * f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = 0.5_dp * f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) + Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( fac(1) * Conjg(c_NuNSn_R(k,i1,i4)) * c_NuNSn_R(k,i2,i4) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                              &
          &               * c_NuNSn_R(k,i1,i4) * Conjg(c_NuNSn_R(k,i2,i4))      &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) - Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( fac(1) * Conjg(c_NuNSn_R(k,i2,i4)) * c_NuNSn_R(k,i1,i4) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                              &
          &               * c_NuNSn_R(k,i2,i4) * Conjg(c_NuNSn_R(k,i1,i4))      &
          &           ) / mN2(i2)
     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) - fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) - fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
       
  End Do

  B_L(4) = kappa_Q *sW2 * B_L(4)
  B_R(4) = kappa_Q *sW2 * B_R(4)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = fact * C_L(4)
  C_R(4) = fact * C_R(4)

  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) - fac(1) * ZSdSdR(i3,i2)   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) - fac(1) * ZSdSdL(i3,i2)   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)
  c_L(5) = c_L(5) * fact
  c_R(5) = c_R(5) * fact

  C11 = Sum(C_L) + Sum(B_L)
  C11p = Sum(C_R) + Sum(B_R)

  If (Present(c11_c)) Then
   C11_c(1) = C11
   C11_c(2:6) = C_L + B_L
  End If
  If (Present(c11p_c)) Then
   C11p_c(1) = C11p
   C11p_c(2:6) = C_R + B_R
  End If

  Iname = Iname - 1

 End Subroutine C_QdQdNuNu_LR


 Subroutine Delta_F2_Boxes(i, j, T3, mf_u, mf_d, mC, mC2, mN, mN2, mGlu, mSpm2 &
   & , mSup2, mSdown2, c_Wqqp_L, c_Hqqp_Lin, c_Hqqp_Rin                        &
   & , c_CQSq_Lin, c_CQSq_Rin, c_QGSq_Lin, c_QGSq_Rin, c_QNSq_Lin, c_QNSq_Rin  &
   & , B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2  )
 !---------------------------------------------------------------------------
 ! gives the box contribution to Delta_F =2 box contributions:
 !     B_VLL  -> bar(q)_Li \gamma_mu q_Lj  * \bar(q)_Li \gamma^mu q_Lj
 !     B_VRR  -> bar(q)_Ri \gamma_mu q_Rj  * \bar(q)_Ri \gamma^mu q_Rj
 !     B_LR1  -> bar(q)_Li \gamma_mu q_Lj  * \bar(q)_Ri \gamma^mu q_Rj
 !     B_LR2  -> bar(q)_Ri q_Lj  * \bar(q)_Li q_Rj
 !     B_SLL1 -> bar(q)_Ri q_Lj  * \bar(q)_Ri q_Lj
 !     B_SLL2 -> bar(q)_Ri \sigma^\mu\nu q_Lj  * \bar(q)_Ri \sigma_\mu\nu q_Lj 
 !     B_SRR1 -> bar(q)_Li q_Rj  * \bar(q)_Li q_Rj
 !     B_SRR2 -> bar(q)_Li \sigma^\mu\nu q_Rj  * \bar(q)_Li \sigma_\mu\nu q_Rj 
 ! using formulas of Buras et al., NPB659 (2003) 3, for charged Higgs and
 ! charginos, Appendix A4
 ! for neutralinos/gluinos the formulas of S.Baek et al. are used,
 ! PRD64 (2001) 095001
 ! Input: i, j .... indices of the outside quarks
 !        T3 ...... isospin of the outside left quarks
 !         mC, mN, mGlu, mSpm, mS0, mP0, RS0, RP0, mSup2, mSdown2
 !         U, V, N, phi_g, RSpm, RSup, RSdown
 !        
 ! Output: G^2_F m^2_W * (V^tb* V^tq)**2 *
 !               (C_VLL, C_VRR, C_SLL1, C_SLL2, C_SRR1, C_SRR2, C_LR1, C_LR2)
 !         
 ! written by Werner Porod, 08 Jul 2003
 !   08.07.03: charged Higgs, for down-quarks
 !---------------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: c_Hqqp_Lin(2,3,3), c_Hqqp_Rin(2,3,3)           &
     & , c_Wqqp_L(3,3), c_QGSq_Lin(3,6), c_QGSq_Rin(3,6), c_CQSq_Lin(:,:,:) &
     & , c_CQSq_Rin(:,:,:), c_QNSq_Lin(:,:,:), c_QNSq_Rin(:,:,:)
  Real(dp), Intent(in) ::  mC(:), mN(:), mC2(:), mN2(:), mGlu        &
     & , mSpm2(:), mSup2(6), mSdown2(6), T3, mf_u(3), mf_d(3)
  Integer, Intent(in) :: i, j

  Complex(dp), Intent(out) ::  B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2   &
     & , B_SRR1, B_SRR2

  Real(dp) :: D0m2, D2m2, mq2(3), mq(3), mqij, mSf2(6), mSf2p(6)    &
     & , mf(3), mg2
  complex(dp) :: c_Hqqp_L(2,3,3), c_Hqqp_R(2,3,3), c_CQSq_L(Size(mC),3,6)   &
     & , c_CQSq_R(Size(mC),3,6), c_QGSq_L(3,6), c_QGSq_R(3,6)               &
     & , c_QNSq_L(3,Size(mN),6), c_QNSq_R(3,Size(mN),6)
  Integer :: i1, i2, i3, i4, n_c, n_n, n_sp

  B_VLL = 0._dp
  B_VRR = 0._dp
  B_LR1 = 0._dp
  B_LR2 = 0._dp
  B_SLL1 = 0._dp
  B_SLL2 = 0._dp
  B_SRR1 = 0._dp
  B_SRR2 = 0._dp

  If (T3.Lt.0._dp) Then
   mq = mf_u
   mq2 = mq**2
   mf = mf_d
   mSf2 = mSdown2
   mSf2p = mSup2
   c_Hqqp_L = Conjg(c_Hqqp_Rin)
   c_Hqqp_R = Conjg(c_Hqqp_Lin)
  Else ! T3 >0
   mq = mf_d
   mq2 = mq**2
   mf = mf_u
   mSf2p = mSdown2
   mSf2 = mSup2
   c_Hqqp_L = c_Hqqp_Lin
   c_Hqqp_R = c_Hqqp_Rin
  End If 
  !-------------------------------------------------
  ! I defined left/right in view of SUSY fermion 
  !-------------------------------------------------
  c_CQSq_L = Conjg(c_CQSq_Rin)
  c_CQSq_R = Conjg(c_CQSq_Lin)
  c_QGSq_L = Conjg(c_QGSq_Rin)
  c_QGSq_R = Conjg(c_QGSq_Lin)
  c_QNSq_L = Conjg(c_QNSq_Rin)
  c_QNSq_R = Conjg(c_QNSq_Lin)


  !------------------------------
  ! W+ H+ + H+ H+ contributions
  !------------------------------

  Do i1=1,3
   Do i2=1,3
    mqij = mq(i1) * mq(i2)
    Do i3=1,2
     Do i4=1,2
      D0m2 = - mqij * D0_Bagger(mSpm2(i3), mSpm2(i4), mq2(i1), mq2(i2) )
      D2m2 = 4._dp * D27_Bagger(mSpm2(i3), mSpm2(i4), mq2(i1), mq2(i2) )
      If ((i3.Eq.1).And.(i4.Eq.2)) Then
!       B_VLL = B_VLL - Conjg(c_Wqqp_L(i1,i)) * c_Wqqp_L(i2,j)                &
       B_VLL = B_VLL - Conjg(c_Wqqp_L(i1,i)) * c_Wqqp_L(i2,j) / 2            &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D0m2  &
!            &        + 0.25_dp *  Conjg(c_Hqqp_L(1,i,i1)) * c_Hqqp_L(1,j,i2) &
            &        - 0.25_dp *  Conjg(c_Hqqp_L(1,i,i1)) * c_Hqqp_L(1,j,i2) &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D2m2  
      Else If ((i3.Eq.2).And.(i4.Eq.2)) Then
!       B_VLL = B_VLL + 0.125_dp * Conjg(c_Hqqp_L(2,i,i1)) * c_Hqqp_L(2,j,i2) &
       B_VLL = B_VLL - 0.125_dp * Conjg(c_Hqqp_L(2,i,i1)) * c_Hqqp_L(2,j,i2) &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D2m2  
       End If
!      B_VRR = B_VRR + 0.125_dp * Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
      B_VRR = B_VRR - 0.125_dp * Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &          * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_R(i3,j,i1) * D2m2
      B_SLL1 = B_SLL1 + 0.5_dp * Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_L(i4,j,i2) &
            &      * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D0m2
      B_SRR1 = B_SRR1 + 0.5_dp * Conjg(c_Hqqp_L(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &      * Conjg(c_Hqqp_L(i3,i,i2)) * c_Hqqp_R(i3,j,i1) * D0m2
!      B_LR1 = B_LR1 + 0.25_dp * Conjg(c_Hqqp_L(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
      B_LR1 = B_LR1 - 0.25_dp * Conjg(c_Hqqp_L(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &         * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D2m2
      B_LR2 = B_LR2 + Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &     * Conjg(c_Hqqp_L(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D0m2
      If (i3.Eq.1) Then
       B_LR2 = B_LR2 - Conjg(c_Wqqp_L(i1,i)) * c_Wqqp_L(i2,j) &
            &            * Conjg(c_Hqqp_R(i4,i,i2)) * c_Hqqp_R(i4,j,i1) * D2m2
      End If
     End Do
    End Do
   End Do
  End Do
  !------------------------------
  ! Chargino contributions
  !------------------------------
  n_c = Size(mC)
  Do i1=1,n_c
   Do i2=1,n_c
    mqij = mC(i1) * mC(i2)
    Do i3=1,6
     Do i4=1,6
      D0m2 = - mqij * D0_Bagger(mC2(i1), mC2(i2), mSf2p(i3), mSf2p(i4) )
!      D2m2 = 4._dp * D27_Bagger(mC2(i1), mC2(i2), mSf2p(i3), mSf2p(i4) )
      D2m2 = - 4._dp * D27_Bagger(mC2(i1), mC2(i2), mSf2p(i3), mSf2p(i4) )
      B_VLL = B_VLL + 0.125_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &                  * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D2m2
      B_VRR = B_VRR + 0.125_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &                  * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D2m2
      B_SLL1 = B_SLL1 - 0.25_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D0m2
      B_SLL2 = B_SLL2 + 0.0625_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4))&
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D0m2
      B_SRR1 = B_SRR1 - 0.25_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &           * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_SRR2 = B_SRR2 + 0.0625_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4))&
            &           * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_LR1 = B_LR1 - 0.5_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_LR2 = B_LR2 - 0.5_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_R(i2,j,i4) * D2m2
     End Do
    End Do
   End Do
  End Do

  !-----------------------------------------------------------
  ! gluino contributions
  !-----------------------------------------------------------
  mg2 = mglu**2
  Do i1=1,6
   Do i2=1,6
    D0m2 = mG2 * D0_Bagger(mg2, mg2, mSf2(i1), mSf2(i2) )
    D2m2 = D27_Bagger(mg2, mg2, mSf2(i1), mSf2(i2) )
    B_VLL = B_VLL - 4._dp *(D0m2 + 11*D2m2)                      &
          &               * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &               * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 9._dp
    B_VRR = B_VRR - 4._dp *(D0m2 + 11*D2m2)                      &
          &               * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &               * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 9._dp
    B_SLL1 = B_SLL1 - 37._dp * D0m2 * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &               * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 9._dp
!    B_SLL2 = B_SLL2 - D0m2 * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
    B_SLL2 = B_SLL2 + D0m2 * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &                * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 12._dp
    B_SRR1 = B_SRR1 - 37._dp * D0m2 * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &               * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 9._dp
!    B_SRR2 = B_SRR2 - D0m2 * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
    B_SRR2 = B_SRR2 + D0m2 * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &                * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 12._dp
    B_LR1 = B_LR1 + 2._dp * (20._dp * D2m2                             &
          &                    * Conjg(c_QGSq_L(i,i2)*c_QGSq_R(i,i1))  &
          &                    * c_QGSq_L(j,i2)*c_QGSq_R(j,i1)         &
          &                 +  Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))    &
          &                    *  (-30._dp * D2m2                      &
          &                            *c_QGSq_L(j,i2)*c_QGSq_R(j,i1)  & 
          &                       + D0m2                               &
          &                           *c_QGSq_L(j,i1)*c_QGSq_R(j,i2) ) &
          &                  ) / 9._dp
    B_LR2 = B_LR2 + ( 48._dp * D2m2                                &
          &              * Conjg(c_QGSq_L(i,i2)*c_QGSq_R(i,i1))    &
          &              * c_QGSq_L(j,i2)*c_QGSq_R(j,i1)           &
          &         + 4._dp * Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2)) &
          &              * (22._dp * D2m2                          &
          &                        * c_QGSq_L(j,i2)*c_QGSq_R(j,i1) & 
          &                -21._dp * D0m2                          &
          &                        * c_QGSq_L(j,i1)*c_QGSq_R(j,i2) &
          &                 ) ) / 9._dp
   End Do 
  End Do
! if (l_mbq(12))  return

  !-------------------------------------------------------------------
  ! gluino/neutralino contributions, the sign in front B27_Bagger is 
  ! the relativ sign to S.Baek et al.
  !-------------------------------------------------------------------
  n_n = Size(mN)
  Do i1=1,6
   Do i2=1,6
    Do i3=1,n_n
     D0m2 = mN(i3) * mGlu * D0_Bagger(mN2(i3), mg2, mSf2(i1), mSf2(i2) )
     D2m2 = D27_Bagger(mN2(i3), mg2, mSf2(i1), mSf2(i2) )
     B_VLL = B_VLL - ( D0m2 * Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2)) &
           &                * c_QGSq_R(j,i1)*c_QGSq_R(j,i2)              & 
           &         + (4._dp * D2m2 * c_QGSq_R(j,i1)                    &
           &                  * Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1))  &
           &           + D0m2 * c_QNSq_R(j,i3,i1)                        &
           &                  * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))     &
           &           ) * c_QNSq_R(j,i3,i2)                             &
           &         ) / 3._dp
     B_VRR = B_VRR - ( D0m2 * Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2)) &
           &                * c_QGSq_L(j,i1)*c_QGSq_L(j,i2)              & 
           &         + (4._dp * D2m2 * c_QGSq_L(j,i1)                    &
           &                  * Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1))  &
           &           + D0m2 * c_QNSq_L(j,i3,i1)                        &
           &                  * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))     &
           &           ) * c_QNSq_L(j,i3,i2)                             &
           &         ) / 3._dp
     B_SRR1 = B_SRR1 + D0m2 * ( Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))      &
           &                      *c_QGSq_L(j,i1)*c_QGSq_L(j,i2)                &
           &                  + (-7._dp*Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1)) &
           &                       *c_QGSq_L(j,i1)                              &
           &                    + Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))          &
           &                      *c_QNSq_L(j,i3,i1) ) *c_QNSq_L(j,i3,i2)       &
           &                  ) / 3._dp
!     B_SRR2 = B_SRR2 - D0m2 * ( Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
     B_SRR2 = B_SRR2 + D0m2 * ( Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
           &                    * c_QGSq_L(j,i1)*c_QGSq_L(j,i2)             &
           &                  + ( Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1))   &
           &                     * c_QGSq_L(j,i1)                           & 
           &                    + Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))      &
           &                      *c_QNSq_L(j,i3,i1) ) * c_QNSq_L(j,i3,i2)  &
            &                 ) / 12._dp
     B_SLL1 = B_SLL1 + D0m2 * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))      &
           &                      *c_QGSq_R(j,i1)*c_QGSq_R(j,i2)                &
           &                  + (-7._dp*Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1)) &
           &                       *c_QGSq_R(j,i1)                              &
           &                    + Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))          &
           &                      *c_QNSq_R(j,i3,i1) ) *c_QNSq_R(j,i3,i2)       &
           &                  ) / 3._dp
!     B_SLL2 = B_SLL2 - D0m2 * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))  &
     B_SLL2 = B_SLL2 + D0m2 * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))  &
           &                    * c_QGSq_R(j,i1)*c_QGSq_R(j,i2)             &
           &                  + ( Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1))   &
           &                     * c_QGSq_R(j,i1)                           & 
           &                    + Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))      &
           &                      *c_QNSq_R(j,i3,i1) ) * c_QNSq_R(j,i3,i2)  &
            &                 ) / 12._dp
     B_LR1 = B_LR1 - 2._dp * ( ( D0m2                                          &
           &                    * ( Conjg(c_QGSq_L(i,i2)*c_QNSq_R(i,i3,i1))    &
           &                          * c_QGSq_R(j,i1)*c_QNSq_L(j,i3,i2)       &
           &                      + Conjg(c_QGSq_R(i,i1)*c_QNSq_L(i,i3,i2))    &
           &                          *c_QGSq_L(j,i2)* c_QNSq_R(j,i3,i1)       &
           &                      ) ) / 6._dp                                  &
           &                  + D2m2                                           &
           &                    * (( Conjg(c_QNSq_R(i,i3,i2))*c_QGSq_L(j,i2)   &
           &                       + Conjg(c_QGSq_R(i,i2))*c_QNSq_L(j,i3,i2))  &
           &                      * (Conjg(c_QNSq_L(i,i3,i1))*c_QGSq_R(j,i1)   &
           &                        + Conjg(c_QGSq_L(i,i1))*c_QNSq_R(j,i3,i1)) &
           &                      + (Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_R(i,i3,i2))&
           &                          *c_QGSq_L(j,i1)* c_QGSq_R(j,i2)          &
           &                        + Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))     &
           &                             *c_QNSq_L(j,i3,i1)*c_QNSq_R(j,i3,i2) )&
           &                        / 3._dp )                                   &
           &                      ) 
     B_LR2 = B_LR2 - 4._dp * ( D0m2                                            &
           &                   * ( Conjg(c_QGSq_L(i,i2)*c_QNSq_R(i,i3,i1))     &
           &                       * c_QGSq_R(j,i1)*c_QNSq_L(j,i3,i2)          &
           &                     + Conjg(c_QGSq_R(i,i1)*c_QNSq_L(i,i3,i2))     &
           &                         *c_QGSq_L(j,i2)* c_QNSq_R(j,i3,i1)        &
           &                     )  / 2._dp                                    &
           &                 + D2m2                                            &
           &                   * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
           &                        * c_QGSq_L(j,i1)*c_QGSq_R(j,i2)            &
           &                     + (( Conjg(c_QNSq_R(i,i3,i2))*c_QGSq_L(j,i2)  &
           &                        + Conjg(c_QGSq_R(i,i2))*c_QNSq_L(j,i3,i2)) &
           &                      * (Conjg(c_QNSq_L(i,i3,i1))*c_QGSq_R(j,i1)   &
           &                        + Conjg(c_QGSq_L(i,i1))*c_QNSq_R(j,i3,i1)) &
           &                       ) /3._dp                                    &
           &                     + Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))        &
           &                        *c_QNSq_L(j,i3,i1)* c_QNSq_R(j,i3,i2)      &
           &                       )  )
       
    End Do 
   End Do 
  End Do

  !-------------------------------------------------------------------
  ! neutralino contributions, the sign in front B27_Bagger is 
  ! the relativ sign to S.Baek et al.
  !-------------------------------------------------------------------
  Do i1=1,6
   Do i2=1,6
    Do i3=1,n_n
     Do i4=1,n_n
      D0m2 = mN(i3) * mN(i4) * D0_Bagger(mN2(i3), mN2(i4), mSf2(i1), mSf2(i2) )
      D2m2 = D27_Bagger(mN2(i3), mN2(i4), mSf2(i1), mSf2(i2) )
      B_VLL = B_VLL - Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_R(j,i4,i2)              &
            &    * ( 2._dp * D2m2 *Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_R(j,i3,i1) & 
            &      + D0m2*Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_R(j,i4,i1) ) / 4._dp
      B_VRR = B_VRR - Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_L(j,i4,i2)              &
            &    * ( 2._dp * D2m2 *Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_L(j,i3,i1) & 
            &      + D0m2*Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_L(j,i4,i1) ) / 4._dp
      B_SLL1 = B_SLL1 + 0.25_dp * D0m2 * c_QNSq_R(j,i4,i2)                     &
             &           * (2._dp *Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))  &
             &                    *c_QNSq_R(j,i4,i1)                           & 
             &             - Conjg(c_QNSq_L(i,i3,i2))                          &
             &                * ( Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_R(j,i4,i1)   & 
             &                  - Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_R(j,i3,i1) ) ) 
!      B_SLL2 = B_SLL2 - D0m2 * Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_R(j,i4,i2)    &
      B_SLL2 = B_SLL2 + D0m2 * Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_R(j,i4,i2)    &
             &            * ( Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_R(j,i4,i1)     & 
             &              - Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_R(j,i3,i1) )   &
             &            / 16._dp
      B_SRR1 = B_SRR1 + 0.25_dp * D0m2 * c_QNSq_L(j,i4,i2)                     &
             &           * (2._dp *Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
             &                    *c_QNSq_L(j,i4,i1)                           & 
             &             - Conjg(c_QNSq_R(i,i3,i2))                          &
             &                * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)   & 
             &                  - Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) ) ) 
!      B_SRR2 = B_SRR2 - D0m2 * Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_L(j,i4,i2)    &
      B_SRR2 = B_SRR2 + D0m2 * Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_L(j,i4,i2)    &
             &            * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)     & 
             &              - Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) )   &
             &            / 16._dp
      B_LR1 = B_LR1 + 0.5_dp * Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_R(j,i4,i2)     &
            &    * ( 2._dp * D2m2*Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_L(j,i4,i1) & 
            &      + D0m2 * Conjg(c_QNSq_L(i,i4,i1)) * c_QNSq_L(j,i3,i1) )
      B_LR2 = B_LR2 + 2._dp*D2m2*Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_R(j,i4,i2)  &
            &           * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)       &
            &             + Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) ) 
     End Do 
    End Do 
   End Do 
  End Do

 End Subroutine Delta_F2_Boxes

 Subroutine Delta_MB(i, m_t, mW2, CKM, B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2 &
      & , B_SRR1, B_SRR2, res )
 !---------------------------------------------------------------------------
 ! Input: mf_u, mC, mN, mGlu, mS0, mP0, mSpm
 !        U, V, N, C
 !        
 ! Output: 
 !         res
 !         
 ! written by Werner Porod, 01 Jul 2003
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: m_t, mW2
  Complex(dp), Intent(in) :: CKM(3,3), B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1 &
     & , B_SLL2, B_SRR1, B_SRR2 
  Integer, Intent(in) :: i

  real(dp) :: xt
  Complex(dp), Intent(out) :: res
  Real(dp), Parameter :: oo4r = 0.25_dp / 0.985_dp, P1_LR = -0.71_dp   &
    & , P2_LR = 0.9_dp, P1_SLL = -0.37_dp, P2_SLL = -0.72_dp 

  xt = m_t**2 / mW2

  res = oo6pi2 * MassBq(i) * etaB * FBhatB(i)**2                            &
      &     * ( (G_F*mW)**2 * CKM(3,i)**2 * Conjg(CKM(3,3))**2 * S0low(xt)  &
      &          + oo4r * (B_VLL + B_VRR) + P1_LR * B_LR1  + P2_LR * B_LR2  &
      &          + P1_SLL * (B_SLL1 + B_SRR1) + P2_SLL * (B_SLL2 + B_SRR2)  )

 End Subroutine Delta_MB


 Real(dp) Function DeltaRho(mZ2, mW2, mP02, RP0, mSneutrino2, RSneutrino     &
           & , mSlepton2, RSlepton, mUsquark2, RUsquark, mDSquark2, RDSquark &
           & , mC, U, V, mN, N)
 !---------------------------------------------------------------------------
 ! calculates the rho-parameter taking into account the complete SM + SUSY
 ! contributions at 1-loop, 2-loop QCD corrections due to t-quark
 ! using tree-level Higgs masses because otherwise there is a problem 
 ! with gauge invariance.
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mZ2, mW2, mP02(2), RP0(2,2), mSneutrino2(3)    &
     & , mSlepton2(6), mUsquark2(6), mDSquark2(6), mC(2), mN(4)
  Complex(dp), Intent(in) :: RSneutrino(3,3), RSlepton(6,6)  &
    & , RUsquark(6,6), RDSquark(6,6), U(2,2), V(2,2), N(4,4)

  Real(dp) :: mC2(2), mN2(4), vev2, v2(2), mHp2(2), mh2(2) &
          & , gSU2, sinW2, t2a, alpha, Rh0(2,2), mu_old, Drho_top
  Complex(dp) :: RSpm(2,2), dmZ2, dmW2

  mC2 = mC**2
  mN2 = mN**2

  gSU2 = Sqrt(4._dp * Sqrt2 * G_F * mW2)

  sinW2 = 1._dp - mW2/mZ2 

  vev2 = 2._dp * Sqrt(mW2)/gSU2
  v2(1) = -vev2 * rp0(1,1)
  v2(2) = vev2  * rp0(1,2)

  mHp2(1) = mW2
  mHp2(2) = mP02(2) + mW2
  RSpm = RP0

  mh2(1) = 0.5_dp * (mP02(2)+mZ2 - Sqrt( (mP02(2)+mZ2)**2 &
         &                   - 4._dp *mZ2*mP02(2)*(rp0(1,1)**2-rp0(1,2)**2)**2))
  mh2(2) = 0.5_dp * (mP02(2)+mZ2 + Sqrt( (mP02(2)+mZ2)**2 &
         &                   - 4._dp *mZ2*mP02(2)*(rp0(1,1)**2-rp0(1,2)**2)**2))
  t2a = (mP02(2)+mZ2)/(mP02(2)-mZ2) &
      &  *2._dp*rp0(1,1)*rp0(1,2)/(rp0(1,1)**2-rp0(1,2)**2)
  alpha = Atan(t2a)*0.5_dp

  rh0(1,1) = -Sin(alpha)
  rh0(1,2) = Cos(alpha)
  rh0(2,1) = Rh0(1,2)
  rh0(2,2) = -rh0(1,1)
  !--------------------------------------------------------------
  ! precision observable defined at the Z-pole, performing DRbar
  ! renormalization
  ! 1-loop contributions
  !--------------------------------------------------------------
  mu_old = SetRenormalizationScale( mZ2 )
  Call SetLoopMassModel(2,4,2,2,2,6,3)
  Call PiZZT1(0._dp, gSU2, sinW2, v2, mZ2, mW2, mh2, Rh0, mP02, RP0        &
   & , mHp2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUsquark2 &
   & , RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V   &
   & , mN, mN2, N, dmZ2)

  Call PiWWT1(0._dp, gSU2, sinW2, mh2, Rh0, mHp2, RSpm , v2                &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2, mW2, dmW2)

  !-----------------------------------------------------------------
  ! 1-loop top contributions, I am only interested in the SUSY part
  !-----------------------------------------------------------------
  Drho_top = 3 * G_F * mf_u2(3) * oosqrt2 * oo8pi2
  !--------------------------------------------------------------
  ! resetting renormalization scale
  !--------------------------------------------------------------
  mu_old = SetRenormalizationScale( mu_old )

  DeltaRho = dmZ2 / mZ2 - dmW2 / mW2 - Drho_top

 End Function DeltaRho


 Subroutine epsilon_K(m_d, m_s, mf_u, mW2, CKM, K_VLL, K_VRR, K_LR1, K_LR2 &
      & , K_SLL1, K_SLL2, K_SRR1, K_SRR2, res, res_SM, Delta_MK )
 !---------------------------------------------------------------------------
 ! Input: mf_u, mC, mN, mGlu, mS0, mP0, mSpm
 !        U, V, N, C
 !        
 ! Output: 
 !         res
 !         
 ! written by Werner Porod, 07 Aug. 2010
 ! for QCD corrections, see Buras, 
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_u(3), m_s, m_d, mW2
  Complex(dp), Intent(in) :: CKM(3,3), K_VLL, K_VRR, K_LR1, K_LR2, K_SLL1 &
      & , K_SLL2, K_SRR1, K_SRR2
  Real(dp), Intent(out) :: res
  Real(dp), Intent(out), optional :: Delta_MK(2), res_SM

  Real(dp) :: xt, xc, r_m

  Real(dp), Parameter :: BK = 0.725_dp !  +- 0.008 +- 0.028 
  Real(dp), Parameter :: DMK_b_VLL = 0.61_dp, DMK_b_SLL1 = 0.76_dp  &
    & , DMK_b_SLL2 = 0.51_dp, DMK_b_LR1 = 0.96_dp, DMK_b_LR2 = 1.3_dp
  !--------------------------------------------------------
  ! by S.Herrlich, U.Nierste, NPB476 (1996) 27
  !--------------------------------------------------------
  Real(dp), Parameter :: DMK_eta_tt = 0.57_dp, DMK_eta_ct = 0.47_dp &
    & , DMK_eta_cc = 1.44_dp
  Complex(dp) :: KK_mat

  xt = mf_u(3)**2 / mW**2
  xc = mf_u(2)**2 / mW**2

  !----------------------------
  ! SM contribution
  !----------------------------
  KK_mat = DMK_eta_cc * (Conjg(CKM(2,2))*CKM(2,1))**2 * xc          &
       & + DMK_eta_tt * (Conjg(CKM(3,2))*CKM(3,1))**2 * S0low(xt)   &
       & + Conjg(CKM(2,2)*CKM(3,2))*(CKM(2,1)*CKM(3,1)) &
       &   * 2._dp * DMK_eta_ct * S0_2(xc,xt)

  KK_mat = MassK0 * FK**2 * BK * oo4pi2 * (G_F*mW)**2 * KK_mat / 3._dp
  If (Present(res_SM)) res_SM = - Aimag(KK_mat) / (sqrt2*DeltaMK)
  If (Present(Delta_MK)) Delta_MK(1) = 2._dp * Real(KK_mat,dp)

  !------------------------------
  ! adding SUSY contributions
  !------------------------------
  r_m = (MassK0 / (m_s + m_d))**2
  KK_mat = KK_mat + oo16pi2 * MassK0 * FK**2  / 24._dp               &
           &    * ( 8._dp * DMK_b_VLL * (K_VLL + K_VRR)              &
           &      -  r_m * (  5._dp * DMK_b_SLL1 * (K_SLL1 + K_SRR1) &
           &              + 12._dp * DMK_b_SLL2 * (K_SLL2 + K_SRR2)  &
           &              +  4._dp * DMK_b_LR1 * K_LR1               &
           &              -  6._dp * DMK_b_LR2 * K_LR2             ) )

  res = - Aimag(KK_mat) / (sqrt2*DeltaMK)
  If (Present(Delta_MK)) Delta_MK(2) = 2._dp * Real(KK_mat,dp)
  
 End Subroutine epsilon_K


 Subroutine Gminus2a(gp,g, mf, Y_l, mSn2, mSl2, Rsl, mC, U, V, mN, N, a_mu)
 !-----------------------------------------------------------------------
 ! Calculates g-2 of a lepton using the formulas given in
 ! T.Ibrahi and P.Nath, PRD61, 095008 (2000)
 ! written by Werner Porod, 21.03.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gp, g, mf, mSn2, mSl2(:), mC(:), mN(:)
  Real(dp), Intent(out) :: a_mu
  Complex(dp), Intent(in) :: Y_l, Rsl(:,:), U(:,:), V(:,:), N(:,:)

  Integer :: i1, i2, n_C, n_N
  Real(dp) :: ratio
  Complex(dp) :: coupL, coupR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Gminus2a'
       
  n_C = 2
  n_N = 4

  a_mu = 0._dp
  !-----------------------
  ! chargino contribution
  !-----------------------
  Do i1=1,n_C
   Call CoupCharginoSfermion(i1, 1, g, -0.5_dp, id2C, Y_l, ZeroC, U, V, &
                             coupL, coupR)
   ratio = mSn2 / mC(i1)**2
   a_mu = a_mu                                                        &
     &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
     &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC(i1)**2 
  End Do
  !-------------------------
  ! neutralino contribution
  !-------------------------
  Do i1=1,n_N
   Do i2 = 1,2
    Call CoupNeutralinoSlepton(i1, i2, gp, g, RSl, Y_l, N, coupL, coupR)
    ratio = mSl2(i2) / mN(i1)**2
    a_mu = a_mu                                                      &
      &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
      &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN(i1)**2
   End Do
  End Do

  a_mu = a_mu * mf * oo16pi2

  Iname = Iname - 1

 End Subroutine Gminus2a

 Subroutine Gminus2b(i_in, mSn2, mC, cpl_CLSn_L, cpl_CLSn_R, mSl2, mN &
                   &, cpl_LNSl_L, cpl_LNSl_R, a_mu, GenerationMixing)
 !-----------------------------------------------------------------------
 ! Calculates g-2 of a lepton using the formulas given in
 ! T.Ibrahim and P.Nath, PRD61, 095008 (2000)
 ! input:
 !  i_in ................ index of the lepton: 1 -> electron, 2 -> muon,
 !                                             3 -> tau
 !  mSn2(i) ............. sneutrino masses squared
 !  mC(i) ............... chargino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino - lepton -sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino - lepton -sneutrino coupling
 !  mSl2(i) ............. slepton masses squared
 !  mN(i) ............... neutralino masses squared
 !  cpl_LNSl_L(i,j,k) ... left lepton - neutralino - slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton - neutralino - slepton coupling
 ! output:
 !  a_mu ................ SUSY contribution to g-2
 ! written by Werner Porod, 07.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mSn2(3), mSl2(6), mC(:), mN(:)
  Real(dp), Intent(out) :: a_mu
  Complex(dp), Intent(in) :: cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:) &
                         & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)
  Logical, Intent(in) :: GenerationMixing

  Integer :: i1, i2, n_C, n_N, n_sf
  Real(dp) :: mf, ratio, mC2, mN2
  Complex(dp) :: coupL, coupR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Gminus2b'
       
  n_C = Size(mC)
  n_N = Size(mN)

  a_mu = 0._dp
  mf = mf_l(i_in)
  !-----------------------
  ! chargino contribution
  !-----------------------
  Do i1=1,n_C
   mC2 =  mC(i1)**2
   If (GenerationMixing) Then
    Do i2 =1,3
     coupL = cpl_CLSn_L(i1,i_in,i2)
     coupR = cpl_CLSn_R(i1,i_in,i2)
     ratio = mSn2(i2) / mC2
     a_mu = a_mu                                                        &
       &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
       &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC2 
    End Do
   Else
    coupL = cpl_CLSn_L(i1,i_in,i_in)
    coupR = cpl_CLSn_R(i1,i_in,i_in)
    ratio = mSn2(i_in) / mC2
    a_mu = a_mu                                                        &
      &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
      &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC2 
   End If
  End Do
  !-------------------------
  ! neutralino contribution
  !-------------------------
  Do i1=1,n_N
   mN2 = mN(i1)**2
   If (GenerationMixing) Then
    Do i2 = 1,6
     coupL = cpl_LNSl_L(i_in,i1,i2)
     coupR = cpl_LNSl_R(i_in,i1,i2)
     ratio = mSl2(i2) / mN2
     a_mu = a_mu                                                      &
       &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
       &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN2
    End Do
   Else
    Do i2 = 1,2
     n_sf = (i_in-1)*2 + i2
     coupL = cpl_LNSl_L(i_in,i1,n_sf)
     coupR = cpl_LNSl_R(i_in,i1,n_sf)
     ratio = mSl2(n_sf) / mN2
     a_mu = a_mu                                                      &
       &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
       &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN2
    End Do
   End If
  End Do

  a_mu = a_mu * mf * oo16pi2

  Iname = Iname - 1

 End Subroutine Gminus2b

 Subroutine Lepton_EDM(i_in, mN, mSl2, cpl_LNSl_L, cpl_LNSl_R    &
                     & , mC, mSnu2, cpl_CLSn_L, cpl_CLSn_R, EDM  )
 !----------------------------------------------------------------------
 ! calculates the EDM-contributions to leptons in the MSSM
 ! input:
 ! output: EDM(i) .... i=1 Neutralino contribution
 !                     i=2 Chargino contribution
 !                     i=3 total EDM
 ! written by Werner Porod, 31.07.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mN(4), mSl2(6), mC(2), mSnu2(3)
  Complex(dp), Intent(in) :: cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:) &
       & , cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)
  Real(dp), Intent(out) :: EDM(3)

  Integer :: i1, i2, i_sl
  Real(dp) :: mF2, r

  EDM = 0._dp
  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  Do i1=1,4
   mF2 = mN(i1)**2
   Do i2=1,2
    i_sl = 2*(i_in-1)+i2
    r = mF2 / mSl2(i_sl)  
    EDM(1) = EDM(1)                                                   &
   &     - Aimag(cpl_LNSl_L(i_in,i1,i_sl) * Conjg( cpl_LNSl_R(i_in,i1,i_sl)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i_sl)  
    L_edm_contri(i1,i2) =   - Aimag(cpl_LNSl_L(i_in,i1,i_sl) * Conjg( cpl_LNSl_R(i_in,i1,i_sl)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i_sl)  
   Enddo
  Enddo

  !------------------------------------
  ! chargino contribution
  !------------------------------------
  Do i1=1,2
   mF2 = mC(i1)**2
   r = mF2 / mSnu2(i_in) 
   EDM(2) = EDM(2)                  &
   &     + Aimag(cpl_CLSn_R(i1,i_in,i_in) * Conjg( cpl_CLSn_L(i1,i_in,i_in)) )&
   &     * FeynFunctionA(r) * mC(i1) / mSnu2(i_in)
   l_edm_contri(5,i1) = Aimag(cpl_CLSn_R(i1,i_in,i_in) * Conjg( cpl_CLSn_L(i1,i_in,i_in)) )&
   &     * FeynFunctionA(r) * mC(i1) / mSnu2(i_in)
  Enddo

  EDM(3) = EDM(1) + EDM(2)
  EDM = oo16pi2 * ecmfactor * EDM 
  l_edm_contri =  oo16pi2 * ecmfactor *l_edm_contri 

 End Subroutine Lepton_EDM

 Subroutine Lepton_EDM3(i_in, mN, mSl2, cpl_LNSl_L, cpl_LNSl_R    &
                     & , mC, mSnu2, cpl_CLSn_L, cpl_CLSn_R, EDM  )
 !----------------------------------------------------------------------
 ! calculates the EDM-contributions to leptons in the MSSM
 ! input:
 ! output: EDM(i) .... i=1 Neutralino contribution
 !                     i=2 Chargino contribution
 !                     i=3 total EDM
 ! written by Werner Porod, 31.07.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mN(4), mSl2(6), mC(2), mSnu2(3)
  Complex(dp), Intent(in) :: cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:) &
       & , cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)
  Real(dp), Intent(out) :: EDM(3)

  Integer :: i1, i2
  Real(dp) :: mF2, r

  EDM = 0._dp
  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  Do i1=1,4
   mF2 = mN(i1)**2
   Do i2=1,6
    r = mF2 / mSl2(i2)  
    EDM(1) = EDM(1)                                                   &
   &     - Aimag(cpl_LNSl_L(i_in,i1,i2) * Conjg( cpl_LNSl_R(i_in,i1,i2)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i2)  
   Enddo
  Enddo

  !------------------------------------
  ! chargino contribution
  !------------------------------------
  Do i1=1,2
   mF2 = mC(i1)**2
   Do i2=1,3
    r = mF2 / mSnu2(i2) 
    EDM(2) = EDM(2)                  &
    &     + Aimag(cpl_CLSn_R(i1,i_in,i2) * Conjg( cpl_CLSn_L(i1,i_in,i2)) )&
    &     * FeynFunctionA(r) * mC(i1) / mSnu2(i2)
   End Do
  Enddo

  EDM(3) = EDM(1) + EDM(2)
  EDM = oo16pi2 * ecmfactor * EDM 

 End Subroutine Lepton_EDM3

 Subroutine Neutron_EDM(mN, mC, mglu, mSup2, c_UNSu_L, c_UNSu_R              &
            & , c_GUSu_L, c_GUSu_R, c_CDSu_L, c_CDSu_R, mSdown2              &
            & , c_DNSd_L, c_DNSd_R, c_GDSd_L, c_GDSd_R, c_CUSd_L, c_CUSd_R   &
            & , GenerationMixing, L_cq, L_rqp, EDMs, EDMS_parts, edm_parts2)
 !-------------------------------------------------------------------
 ! calculates the neutron EDM
 ! written by Werner Porod
 ! 06.07.03: adding charged Higgs boson
 !           adding optional argument EDMS_parts(2,i=1,..,4)
 !           for the contributions: i=1...neutralino
 !                                  i=2...chargino
 !                                  i=3...gluino
 !                                  i=4...charged Higgs
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mN(4), mC(2), mglu, mSup2(6), mSdown2(6)
  Complex(dp), Intent(in) :: c_UNSu_L(3,4,6), c_UNSu_R(3,4,6), c_GUSu_L(3,6) &
       & , c_GUSu_R(3,6), c_CDSu_L(2,3,6), c_CDSu_R(2,3,6), c_DNSd_L(3,4,6)  &
       & , c_DNSd_R(3,4,6), c_GDSd_L(3,6), c_GDSd_R(3,6), c_CUSd_L(2,3,6)    &
       & , c_CUSd_R(2,3,6)
  Logical, Intent(in) :: GenerationMixing, L_cq, L_rqp
  Real(dp), Intent(out) :: EDMs(2)
  Real(dp), Optional, Intent(out) :: EDMS_parts(2,3), edm_parts2(2,6)

  Integer :: i1, i2
  Real(dp) :: r, B_r, A_r, mF2, g_s, d_d_n, dh_d_n, d_u_n, dh_u_n &
    & , d_s_n, d_d_c, dh_d_c, d_u_c, dh_u_c, d_s_c, d_d_g, dh_d_g, d_u_g    &
    & , dh_u_g, d_s_g

  !----------------------------
  ! Initialization
  !----------------------------
  EDMs = 0._dp
  d_d_n = 0._dp  ! neutralino contribution to d-quark EDM
  dh_d_n = 0._dp ! neutralino contribution to d-quark chromoelectric EDM
  d_u_n = 0._dp  ! neutralino contribution to u-quark EDM
  dh_u_n = 0._dp ! neutralino contribution to u-quark chromoelectric EDM
  d_s_n = 0._dp  ! neutralino contribution to s-quark EDM

  d_d_c = 0._dp  ! chargino contribution to d-quark EDM
  dh_d_c = 0._dp ! chargino contribution to d-quark chromoelectric EDM
  d_u_c = 0._dp  ! chargino contribution to u-quark EDM
  dh_u_c = 0._dp ! chargino contribution to u-quark chromoelectric EDM
  d_s_c = 0._dp  ! chargino contribution to s-quark EDM

  d_d_g = 0._dp  ! gluino contribution to d-quark EDM
  dh_d_g = 0._dp ! gluino contribution to d-quark chromoelectric EDM
  d_u_g = 0._dp  ! gluino contribution to u-quark EDM
  dh_u_g = 0._dp ! gluino contribution to u-quark chromoelectric EDM
  d_s_g = 0._dp  ! gluino contribution to s-quark EDM

  g_s = Sqrt(4._dp * Pi * AlphaS_mZ)

  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  If (GenerationMixing) Then
   Do i1=1,4
    mF2 = mN(i1)**2
    Do i2=1,6
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2)
     d_d_n = d_d_n - Aimag(c_DNSd_L(1,i1,i2) * Conjg(c_DNSd_R(1,i1,i2)) ) * B_r
     If (L_rqp) d_s_n = d_s_n   &
          & - Aimag(c_DNSd_L(2,i1,i2) * Conjg(c_DNSd_R(2,i1,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSup2(i2)
     d_u_n = d_u_n - Aimag(c_UNSu_L(1,i1,i2) * Conjg(c_UNSu_R(1,i1,i2)) ) * B_r
    Enddo
   Enddo
 
  Else ! .not.GenerationMixing
   Do i1=1,4
    mF2 = mN(i1)**2
    Do i2=1,2
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2)
     d_d_n = d_d_n - Aimag(c_DNSd_L(1,i1,i2) * Conjg(c_DNSd_R(1,i1,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSup2(i2)
     d_u_n = d_u_n - Aimag(c_UNSu_L(1,i1,i2) * Conjg(c_UNSu_R(1,i1,i2)) ) * B_r
     If (L_rqp) Then
      r = mF2 / mSdown2(i2+2)
      B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2+2)
      d_s_n = d_s_n &
          & - Aimag(c_DNSd_L(2,i1,i2+2) * Conjg(c_DNSd_R(2,i1,i2+2)) ) * B_r
     End If
    Enddo
   Enddo
  End If

  If (L_cq) Then
   dh_d_n = g_s * d_d_n
   dh_u_n = g_s * d_u_n
  Endif
  d_d_n = e_d * d_d_n
  d_s_n = e_d * d_s_n
  d_u_n = e_u * d_u_n
  !------------------------------------
  ! chargino contribution
  !------------------------------------
  If (GenerationMixing) Then
   Do i1=1,2
    mF2 = mC(i1)**2
    Do i2=1,6
     r = mF2 / mSup2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2)
     d_d_c = d_d_c - Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) &
           &         * (e_u * B_r + (e_d - e_u) * A_r)
     If (L_cq) dh_d_c = dh_d_c                        &
          & + g_s * Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) * B_r
     If (L_rqp) d_s_c = d_s_c   &
          & - Aimag(c_CDSu_L(i1,2,i2) * Conjg(c_CDSu_R(i1,2,i2)) ) &
          &         * (e_u * B_r + (e_d - e_u) * A_r)

     r = mF2 / mSdown2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSdown2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSdown2(i2)
     d_u_c = d_u_c - Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) &
          &         * (e_d * B_r + (e_u - e_d) * A_r)
     If (L_cq) dh_u_c = dh_u_c                        &
         & + g_s * Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) * B_r
    Enddo
   Enddo
 
  Else ! .not.GenerationMixing
   Do i1=1,2
    mF2 = mC(i1)**2
    Do i2=1,2
     r = mF2 / mSup2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2)
     d_d_c = d_d_c - Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) &
           &         * (e_u * B_r + (e_d - e_u) * A_r)
     If (L_cq) dh_d_c = dh_d_c                        &
          & + g_s * Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) * B_r

     r = mF2 / mSdown2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSdown2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSdown2(i2)
     d_u_c = d_u_c - Aimag(c_CUSd_L(1,i1,i2) * Conjg(c_CUSd_R(1,i1,i2)) ) &
          &         * (e_d * B_r + (e_u - e_d) * A_r)
     If (L_cq) dh_u_c = dh_u_c                        &
         & + g_s * Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) * B_r

     If (L_rqp) Then
      r = mF2 / mSup2(i2+2)
      A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2+2)
      B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2+2)
      d_s_c = d_s_c &
          & - Aimag(c_CDSu_L(i1,2,i2+2) * Conjg(c_CDSu_R(i1,2,i2+2)) ) &
          &         * (e_u * B_r + (e_d - e_u) * A_r)
     End If
    Enddo
   Enddo
  End If
  !------------------------------------
  ! gluino contribution
  !------------------------------------
  mF2 = mGlu**2
  If (GenerationMixing) Then
    Do i2=1,6
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSdown2(i2)
     d_d_g = d_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSdown2(i2) - B_r/3._dp
      dh_d_g = dh_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * A_r
     end if
     If (L_rqp) d_s_g = d_s_g   &
          & - Aimag(c_GDSd_L(2,i2) * Conjg(c_GDSd_R(2,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSup2(i2)
     d_u_g = d_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSup2(i2) - B_r/3._dp
      dh_u_g = dh_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * A_r
     end if
    Enddo
 
  Else ! .not.GenerationMixing
    Do i2=1,2
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSdown2(i2)
     d_d_g = d_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSdown2(i2) - B_r/3._dp
      dh_d_g = dh_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * A_r
     end if

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSup2(i2)
     d_u_g = d_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSup2(i2) - B_r/3._dp
      dh_u_g = dh_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * A_r
     end if
     If (L_rqp) Then
      r = mF2 / mSdown2(i2+2)
      B_r = FeynFunctionB(r) * mGlu / mSdown2(i2+2)
      d_s_g = d_s_g &
          & - Aimag(c_GDSd_L(2,i2+2) * Conjg(c_GDSd_R(2,i2+2)) ) * B_r
     End If
    Enddo
  End If

  d_d_g = 16._dp * e_d * d_d_g / 3._dp
  d_s_g = 16._dp * e_d * d_s_g / 3._dp
  d_u_g = 16._dp * e_u * d_u_g / 3._dp

  !---------------------------------------
  ! neutron EDM in the chiral quark model
  !---------------------------------------
  If (L_cq) Then
   EDMs(1) = qcdCorrection * ( 4._dp * (d_d_n + d_d_c + d_d_g)   &
         &                   - d_u_n - d_u_c -d_u_g) / 3._dp     &
         & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)             &
         &              * ( 4._dp * (dh_d_n + dh_d_c + dh_d_g)   &
         &                - dh_u_n - dh_u_c -dh_u_g) / 3._dp 
if (present(edm_parts2)) then
   edm_parts2(1,1) = d_d_n
   edm_parts2(1,2) = dh_d_n
   edm_parts2(1,3) = d_d_c
   edm_parts2(1,4) = dh_d_c
   edm_parts2(1,5) = d_d_g
   edm_parts2(1,6) = dh_d_g
   edm_parts2(2,1) = d_u_n
   edm_parts2(2,2) = dh_u_n
   edm_parts2(2,3) = d_u_c
   edm_parts2(2,4) = dh_u_c
   edm_parts2(2,5) = d_u_g
   edm_parts2(2,6) = dh_u_g
end if
  End If
  !----------------------------------------------------
  ! neutron EDM in the relativistic quark-parton model
  ! Deltas refer to proton, not to neutron
  !----------------------------------------------------
  If (L_rqp) Then
   EDMs(2) = qcdCorrection * ( DeltaUp      * (d_d_n + d_d_c + d_d_g) &
         &                   + DeltaDown    * (d_u_n + d_u_c + d_u_g) &
         &                   + DeltaStrange * (d_s_n + d_s_c + d_s_g) )
  End If

  EDMs = oo16pi2 * ecmfactor * EDMs

  !-----------------------------------
  ! individual contributions
  !-----------------------------------
  If (present(EDMS_parts)) then
   !  chiral quark model
   EDMs_parts(1,1) = qcdCorrection * ( 4._dp * d_d_n  - d_u_n ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_n - dh_u_n ) / 3._dp
   EDMs_parts(1,2) = qcdCorrection * ( 4._dp * d_d_c  - d_u_c ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_c - dh_u_c ) / 3._dp
   EDMs_parts(1,3) = qcdCorrection * ( 4._dp * d_d_g  - d_u_g ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_g - dh_u_g ) / 3._dp
 
   ! relativistic quark-parton model
   EDMs_parts(2,1) = qcdCorrection * ( DeltaDown * d_d_n + DeltaUp * d_u_n  &
                   &                 + DeltaStrange * d_s_n  )
   EDMs_parts(2,2) = qcdCorrection * ( DeltaDown * d_d_c + DeltaUp * d_u_c  &
                   &                 + DeltaStrange * d_s_c  )
   EDMs_parts(2,3) = qcdCorrection * ( DeltaDown * d_d_g + DeltaUp * d_u_g  &
                   &                 + DeltaStrange * d_s_g  )
   
   EDMs_parts = oo16pi2 * ecmfactor * EDMs_parts
  end if

 End Subroutine Neutron_EDM


 Subroutine Low_Energy_Constraints_MSSM(Qin, gi, Y_l, Y_d, Y_u, T_l, T_d, T_u &
   & , Mi, mu, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, B, tanb_Q, mP02, mS02      &
   & , mSpm2, CKM_in, kont, GenMix, rho_parameter, DeltaMBd, DeltaMBs         &
   & , BRbtosgamma, Bs_ll, Bd_ll, BrBToSLL, BtoSNuNu, BR_Bu_TauNu, R_Bu_TauNu &
   & , epsK, DeltaMK, K0toPi0NuNu, KptoPipNuNu, a_e, a_mu, a_tau, d_e, d_mu   &
   & , d_tau, BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma, BrMu3e, BrTau3e    &
   & , BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau)

 Implicit None
  !---------------------------------
  ! input
  !---------------------------------
  Logical, Intent(in) :: GenMix
  Real(dp), Intent(in) :: Qin, gi(3), M2_H(2), mP02(2), mS02(2), mSpm2(2), tanb_Q
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Y_l, T_d, T_u, T_l  &
                                 & , M2_E, M2_L, M2_D, M2_Q, M2_U, CKM_in
  Complex(dp), Intent(in) :: Mi(3), mu, B
  !---------------------------------
  ! output
  !---------------------------------
  Integer, Intent(out) :: kont
  Real(dp), Intent(out) :: BRbtosgamma, BrBToSLL, BR_Bu_TauNu, a_mu, a_e       &
   & , a_tau, d_e, d_mu, d_tau, BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma    &
   & , BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau, BtoSNuNu &
   & , Bs_ll(3), Bd_ll(3), rho_parameter, epsK, K0toPi0NuNu, KptoPipNuNu       &
   & , R_Bu_TauNu, DeltaMK
  Complex(dp), Intent(out) :: DeltaMBd, DeltaMBs

  !------------------
  ! local variables
  !---------------------------------------------------------------
  ! 160 scale GeV for the calculation of the hadronic observables
  !---------------------------------------------------------------
  Real(dp) :: gi_160(3), M2_H_160(2), tanb_160
  Complex(dp) ::  Y_l_160(3,3), Y_d_160(3,3), Y_u_160(3,3), Mi_160(3)    &
      & , T_l_160(3,3), T_d_160(3,3), T_u_160(3,3), M2_E_160(3,3)        &
      & , M2_L_160(3,3), M2_D_160(3,3), M2_Q_160(3,3), M2_U_160(3,3)     &
      & , mu_160, B_160, T_l_s(3,3), T_d_s(3,3), T_u_s(3,3), M2_E_s(3,3) &
      & , M2_L_s(3,3), M2_D_s(3,3), M2_Q_s(3,3), M2_U_s(3,3), Y_d_s(3,3) &
      & , Y_u_s(3,3), Y_l_s(3,3), CKM_160(3,3)
  !----------------------------------------------------------
  ! at m_Z
  !----------------------------------------------------------
  Real(dp) :: gi_mZ(3), M2_H_mZ(2), mudim_old, tanb_mZ
  Complex(dp) ::  Y_l_mZ(3,3), Y_d_mZ(3,3), Y_u_mZ(3,3), Mi_mZ(3)    &
      & , T_l_mZ(3,3), T_d_mZ(3,3), T_u_mZ(3,3), M2_E_mZ(3,3)        &
      & , M2_L_mZ(3,3), M2_D_mZ(3,3), M2_Q_mZ(3,3), M2_U_mZ(3,3)     &
      & , mu_mZ, B_mZ
  !----------------------------------------------------------
  ! scale independent
  !----------------------------------------------------------
  Complex(dp), Dimension(3,3) :: uD_L, uD_R, uU_L, uU_R, CKM
  Real(dp) :: dt, tz, g2(214), vev2, g, gp, gs,cpl_LLZ_L,cpl_LLZ_R &
     & , vevSM(2), g1(57)

  Complex(dp) :: cpl_uWd(3,3), cpl_CsDU_L(2,3,3), cpl_CsDU_R(2,3,3)       &
     & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)          &
     & , cpl_CLSn_R(2,3,3), cpl_DDS0_L(3,3,2), cpl_DDS0_R(3,3,2)          &
     & , cpl_DDP0_L(3,3,2), cpl_DDP0_R(3,3,2), cpl_DDS0_1L_L(3,3,2)       &
     & , cpl_DDS0_1L_R(3,3,2), cpl_DDP0_1L_L(3,3,2), cpl_DDP0_1L_R(3,3,2) &
     & , cpl_LLS0_L(3,3,2), cpl_LLS0_R(3,3,2), cpl_LLP0_L(3,3,2)          &
     & , cpl_LLP0_R(3,3,2)
  Complex(dp) :: cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6), cpl_NuNSn_R(3,4,3)      &
    & , cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6), cpl_NNZ_L(4,4), cpl_NNZ_R(4,4) &
    & , cpl_DGSd_L(3,6), cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6) &
    & , cpl_CCZ_L(2,2), cpl_CCZ_R(2,2), cpl_SnSnZ(3,3), cpl_SlSlZ(6,6)         &
    & , cpl_CCP0_L(2,2,2), cpl_CCP0_R(2,2,2), cpl_CCS0_L(2,2,2)                &
    & , cpl_CCS0_R(2,2,2), cpl_NNP0_L(4,4,2), cpl_NNP0_R(4,4,2)                &
    & , cpl_NNS0_L(4,4,2), cpl_NNS0_R(4,4,2), cpl_P0SdSd(2,6,6)                &
    & , cpl_P0SuSu(2,6,6), cpl_S0SdSd(2,6,6), cpl_S0SuSu(2,6,6), cpl_dWu(3,3)  &
    & , cpl_UUS0_L(3,3,2), cpl_UUS0_R(3,3,2), cpl_UUP0_L(3,3,2)                &
    & , cpl_UUP0_R(3,3,2), cpl_CsP0W(2,2), cpl_CsS0W(2,2), cpl_CsCsS0(2,2,2)
  Real(dp) :: cpl_S0WW(2)

  Real(dp), Parameter :: T3_d=-0.5_dp, T3_u=0.5_dp

  Real(dp) :: mGlu_T, mC_T(2), mC2_T(2), mN_T(4), mN2_T(4), mSneut_T(3)   &
     & , mSneut2_T(3), mSlept_T(6), mSlept2_T(6), mSdown_T(6), mSdown2_T(6) &
     & , mSup_T(6), mSup2_T(6), mP0_T(2), mP02_T(2), RP0_T(2,2), mS0_T(2) &
     & , mS02_T(2), RS0_T(2,2), mSpm_T(2), mSpm2_T(2),mZ2_run       &
     & , mW2_run, DMK(2), work
  Complex(dp) :: Phase_Glu_T, U_T(2,2), V_T(2,2), N_T(4,4), Rsneut_T(3,3)  &
     & , RSlept_T(6,6), RSdown_T(6,6), RSup_T(6,6), RSpm_T(2,2), bi(1)   &
     & , ZNN(4,4), ZUU(2,2), ZVV(2,2), ZSuSuL(6,6), ZSuSuR(6,6), ZSdSdL(6,6) &
     & , ZSdSdR(6,6)
  Integer :: i1,i2,i3,i4,i5
  Real(dp) :: GMutoEGamma, GTautoEGamma, GTautoMuGamma
  Real(dp) :: BtoSEE, EDM_e(3), EDM_mu(3), EDM_tau(3), gU1, gSU2  &
     & , cosW, sinW2, mf_u_in(3), abs_mu2, mf_d_in(3)

  Real(dp) :: mf_u_Q(3), s12, s13, s23, c12, c23, c13, phase
  Complex(dp) :: epsD(3), epsL(3), epsD_FC
  Logical :: GenMix_save
! should be shifted to input/output
  Real(dp) :: epsK_SM

  Iname = Iname + 1
  NameOfUnit(Iname) = "Low_Energy_Constraints_MSSM"

  !----------------------------------------------------------------
  ! initialisation, in case that somewhere a problem appears
  !----------------------------------------------------------------
  kont = 0
  BRbtosgamma = 0._dp
  BToSNuNu = 0._dp
  BrBToSLL = 0._dp
  DeltaMBd = 0._dp
  DeltaMBs = 0._dp
  Bd_ll = 0._dp
  Bs_ll = 0._dp
  BR_Bu_TauNu = 0._dp
  R_Bu_TauNu = 0._dp

  epsK = 0._dp
  DeltaMK = 0._dp
  K0toPi0NuNu = 0._dp
  KptoPipNuNu = 0._dp

  a_e = 0._dp
  a_mu = 0._dp
  a_tau = 0._dp
  d_e = 0._dp
  d_mu = 0._dp
  d_tau = 0._dp
  BrMutoEGamma = 0._dp
  BrTautoEGamma = 0._dp
  BrTautoMuGamma = 0._dp
  BrMu3e = 0._dp
  BrTau3e = 0._dp
  BrTau3Mu = 0._dp
  BR_Z_e_mu = 0._dp
  BR_Z_e_tau = 0._dp
  BR_Z_mu_tau = 0._dp
  rho_parameter = 0._dp


  !-----------------------------------------------
  ! run to Q=160 GeV if necessary
  !-----------------------------------------------
  Y_l_160 = Transpose(Y_l) 
  Y_d_160 = Transpose(Y_d) 
  Y_u_160 = Transpose(Y_u) 
  T_l_160 = Transpose(T_l) 
  T_d_160 = Transpose(T_d) 
  T_u_160 = Transpose(T_u) 
  M2_E_160 = M2_E
  M2_L_160 = M2_L
  M2_D_160 = M2_D
  M2_U_160 = M2_U
  M2_Q_160 = M2_Q
  Mi_160 = Mi
  mu_160 = mu
  B_160 = B
  M2_H_160 = M2_H

  CKM = CKM_in
  !------------------------------------------------------------
  ! add flavour mixing if necessary
  ! using partially in this first stage variable names for 160 GeV
  ! all of them are reset below 
  ! in this first version only epsilon_D_FC is considered and the
  ! remaining 1-loop corrections to the CKM are neglected;
  ! will be improved later
  !------------------------------------------------------------
  If (.Not.Genmix) Then
   sinW2 = 1._dp - mW2 / mZ2
   vev2 =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev2 / Sqrt(1._dp + tanb_Q**2)
   vevSM(2) = tanb_Q * vevSM(1)

   Call  CouplingsToG(gi, y_l_160, y_d_160, y_u_160, g1)
 
   tz = Log(Qin/mZ)
   dt = - tz / 50._dp

   g1(1) = Sqrt(5._dp/3._dp) * g1(1)
   Call odeint(g1, 57, tz, 0._dp, delta_mass, dt, 0._dp, rge57, kont)
   g1(1) = Sqrt(3._dp/5._dp) * g1(1)
 
   Call  GToCouplings(g1, gi_mZ, y_l_mZ, y_d_mZ, y_u_mZ)
   y_u_160 = Matmul(y_u_mZ, CKM)  ! I am reusing the Yukawas at m_Z below

   Call  CouplingsToG(gi_mZ, Y_l_mZ, Y_d_mZ, y_u_160, g1)

   tz = Log(mZ/Qin)
   dt = - tz / 50._dp

   g1(1) = Sqrt(5._dp/3._dp) * g1(1)
   Call odeint(g1, 57, tz, 0._dp, delta_mass, dt, 0._dp, rge57, kont)
   g1(1) = Sqrt(3._dp/5._dp) * g1(1)

   Call GToCouplings(g1, gi_160, y_l_160, y_d_160, y_u_160)
   Y_l_160 = Transpose(Y_l_160) 
   Y_d_160 = Transpose(Y_d_160) 
   Y_u_160 = Transpose(Y_u_160) 

   Call QuarkMasses_and_PhaseShifts(Y_d_160, Y_u_160, vevSM, mf_d_in, uD_L &
                                & , uD_R, mf_u_in, uU_L, uU_R, CKM_160)

   Call Sigma_SM_chirally_enhanced(gi, vevSM, mu, CKM_160, Y_l, Y_d &
      & , Y_u, Mi, T_l, T_d, T_u, M2_E, M2_L, M2_D, M2_Q, M2_U     &
      & , epsD, epsL, epsD_FC, kont)
   !-------------------------------------------------
   ! recalculate Yukawa, taking epsD_FC into account
   !-------------------------------------------------
   Y_u_mZ = Matmul(Y_u_mZ, CKM)  
   Y_u_mZ(1:2,3) = Y_u_mZ(1:2,3) / (1._dp-epsD_FC)
   Y_u_mZ(3,1:2) = Y_u_mZ(3,1:2) / (1._dp-Conjg(epsD_FC))
   Call  CouplingsToG(gi_mZ, Y_l_mZ, Y_d_mZ, y_u_mZ, g1)

   tz = Log(mZ/Qin)
   dt = - tz / 50._dp

   g1(1) = Sqrt(5._dp/3._dp) * g1(1)
   Call odeint(g1, 57, tz, 0._dp, delta_mass, dt, 0._dp, rge57, kont)
   g1(1) = Sqrt(3._dp/5._dp) * g1(1)

   Call GToCouplings(g1, gi_160, y_l_160, y_d_160, y_u_160)

   Y_d_160 = Transpose(Y_d_160) 
   Y_u_160 = Transpose(Y_u_160) 
   Call QuarkMasses_and_PhaseShifts(Y_d_160, Y_u_160, vevSM, mf_d_in, uD_L &
                                & , uD_R, mf_u_in, uU_L, uU_R, CKM_160)

   Y_d_160 = Transpose(Y_d_160) 
   Y_u_160 = Transpose(Y_u_160) 
   
   T_u_160 = Matmul(Transpose(uU_R),Matmul(T_u_160,uU_L))
   T_d_160 = Matmul(Transpose(uD_R),Matmul(T_d_160,uD_L))
   M2_Q_160 = Matmul(Conjg(Transpose(uD_L)),Matmul(M2_Q_160,uD_L))
   M2_U_160 = Matmul(Transpose(uU_R),Matmul(M2_U_160,Conjg(uU_R)))
   M2_D_160 = Matmul(Transpose(uD_R),Matmul(M2_D_160,Conjg(uD_R)))
  End If

  !-----------------------------------------------------------
  ! bug fix as GenerationMixing is used in various routines
  !-----------------------------------------------------------
  GenMix_save = GenerationMixing
  GenerationMixing = .True.
  Call ParametersToG(gi, y_l_160, y_d_160, y_u_160, Mi_160, T_l_160, T_d_160 &
        & , T_u_160, M2_E_160, M2_L_160, M2_D_160, M2_Q_160, M2_U_160        &
        & , M2_H_160, mu_160, B_160, g2(1:213))

  g2(214) = Log(tanb_Q)

  tz = Log(160._dp/Qin)

  If (tz.Ne.0._dp) Then
   dt = tz / 100._dp
   g2(1) = Sqrt(5._dp / 3._dp ) * g2(1)
   Call odeint(g2, 214, 0._dp, tz, delta_mass, dt, 0._dp, rge214, kont)
   g2(1) = Sqrt(3._dp / 5._dp ) * g2(1)
   tanb_160 = Exp(g2(214))

   Call GToParameters(g2(1:213), gi_160, Y_l_160, Y_d_160, Y_u_160, Mi_160   &
       & , T_l_160, T_d_160, T_u_160, M2_E_160, M2_L_160, M2_D_160, M2_Q_160 &
       & , M2_U_160, M2_H_160, mu_160, B_160)
  Else
   tanb_160 = tanb_Q
   gi_160 = gi
  End If

  tz = Log(mZ / 160._dp)
  dt = tz / 100._dp

  g2(1) = Sqrt(5._dp / 3._dp ) * g2(1)
  Call odeint(g2, 214, 0._dp, tz, delta_mass, dt, 0._dp, rge214, kont)
  tanb_mZ = Exp(g2(214))
  g2(1) = Sqrt(3._dp / 5._dp ) * g2(1)

  Call GToParameters(g2(1:213), gi_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, T_l_mZ &
                  & , T_d_mZ, T_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ     &
                  & , M2_U_mZ, M2_H_mZ, mu_mZ, B_mZ)

   Y_l_160 = Transpose(Y_l_160) 
   Y_d_160 = Transpose(Y_d_160) 
   Y_u_160 = Transpose(Y_u_160) 
   T_l_160 = Transpose(T_l_160) 
   T_d_160 = Transpose(T_d_160) 
   T_u_160 = Transpose(T_u_160) 

   Y_l_mZ = Transpose(Y_l_mZ) 
   Y_d_mZ = Transpose(Y_d_mZ) 
   Y_u_mZ = Transpose(Y_u_mZ) 
   T_l_mZ = Transpose(T_l_mZ) 
   T_d_mZ = Transpose(T_d_mZ) 
   T_u_mZ = Transpose(T_u_mZ) 

  !-------------------------------------
  ! calculate running masses at 160 GeV
  ! has to be improved
  !-------------------------------------
  mudim_old = SetRenormalizationScale( 160._dp**2 )

  sinW2 = 1._dp - mW2 / mZ2
  vev2 =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
  vevSM(1) = vev2 / Sqrt(1._dp + tanb_160**2)
  vevSM(2) = tanb_160 * vevSM(1)

  gp = gi_160(1)
  g = gi_160(2)
  gs = gi_160(3)
  mZ2_run = (gp**2+g**2)*0.25_dp*(vevSM(1)**2+vevSM(2)**2)
  mW2_run = g**2*0.25_dp*(vevSM(1)**2+vevSM(2)**2)

  abs_mu2= (M2_H_160(2) * tanb_160**2 - M2_H_160(1) ) / (1._dp-tanb_160**2) &
         & - 0.5_dp * mZ2_run
  mu_160 = phase_mu * Sqrt(abs_mu2)
  B_160 = (M2_H_160(1) + M2_H_160(2) + 2._dp * Abs_Mu2) * tanb_160 &
      & / (1+tanb_160**2)

  Call TreeMassesMSSM2(gp, g, vevSM, Mi_160(1), Mi_160(2), Mi_160(3)         &
    & , mu_160, B_160, tanb_160, M2_E_160, M2_L_160, T_l_160, Y_l_160        &
    & , M2_D_160, M2_U_160, M2_Q_160, T_d_160, T_u_160, Y_d_160, Y_u_160     &
    & , uU_L, uU_R ,uD_L, uD_R, uL_L, uL_R                                   &
    & , mGlu_T, Phase_Glu_T, mC_T, mC2_T, U_T, V_T, mN_T, mN2_T, N_T         &
    & , mSneut_T, mSneut2_T, Rsneut_T, mSlept_T, mSlept2_T       &
    & , RSlept_T, mSdown_T, mSdown2_T, RSdown_T, mSup_T, mSup2_T, RSup_T   &
    & , mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T, mSpm_T, mSpm2_T, RSpm_T  &
    & , mZ2_run, mW2_run, .True., kont, .False., .False., mf_u_Q)

  If (kont.Ne.0) Then ! there is a problem with the running masses
   Iname = Iname - 1
   kont = -700
   Call AddError(700)
   Return
  End If 

  CKM_160 =  Matmul(uU_L, Transpose(Conjg(uD_L)) )

  If (Aimag(CKM(1,3)).Eq.0._dp) Then ! real CKM, need to 
    s12 = lam_wolf                   ! recalculate with phase
    s23 = s12**2 * A_wolf
    s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2)
    phase = Atan(eta_wolf/rho_wolf)

    c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    CKM(1,1) = c12 * c13
    CKM(1,2) = s12 * c13
    CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
    CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM(2,3) = s23 * c13
    CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM(3,3) = c23 * c13
    !---------------------------------------------
    ! and now the one at 160 GeV
    !---------------------------------------------
    s13 = Abs(CKM_160(1,3))              ! to be recalculated
    c13 = Sqrt(1._dp - s13**2)
    s12 = Abs(CKM_160(1,2)) / c13
    c12 = Sqrt(1._dp - s12**2)
    s23 = Abs(CKM_160(2,3)) / c13
    c23 = Sqrt(1._dp - s23**2)
    phase = Atan(eta_wolf/rho_wolf)

    CKM_160(1,1) = c12 * c13
    CKM_160(1,2) = s12 * c13
    CKM_160(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
    CKM_160(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM_160(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM_160(2,3) = s23 * c13
    CKM_160(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM_160(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    CKM_160(3,3) = c23 * c13

    Y_u_160 = Matmul(Conjg(uU_L),Y_u_160) ! partially rotate
    uU_L = Matmul(CKM_160,uD_L)               ! new mixing with phases
    Y_u_160 = Matmul(Transpose(uU_L),Y_u_160)     ! rotate back
    !---------------------------------------------------------------
    ! recalculate u-squark masses and mixing
    !---------------------------------------------------------------
    Call SfermionMass3mssm(M2_Q_160, M2_U_160, T_u_160, mu_160, vevSM  &
      & , Y_u_160, 0.5_dp,  1._dp / 3._dp, -4._dp / 3._dp, g, gp, kont &
      & , mSup_T, mSup2_T, RSup_T)
  End If

  !----------------------------
  ! couplings 
  !----------------------------

  cpl_uWd = g * oosqrt2 * CKM_160

  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Call CoupChargedScalarFermion3(i3, i1, i2, RSpm_T, Y_d_160, uD_L, uD_R   &
        & , Y_u_160, uU_L, uU_R, cpl_CsDU_L(i3,i1,i2), cpl_CsDU_R(i3,i1,i2) )
     Call CoupFermionScalar3(i1, i2, i3, T3_u, Y_u_160, uU_L, uU_R, RS0_T     &
                & , cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3) )
     Call CoupFermionPseudoScalar3(i1, i2, i3, T3_u, Y_u_160, uU_L, uU_R      &
                & , RP0_T, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3) )
     Call CoupFermionScalar3(i1, i2, i3, T3_d, Y_l_160, uL_L, uL_R, RS0_T &
               & , cpl_LLS0_L(i1,i2,i3), cpl_LLS0_R(i1,i2,i3))
     Call CoupFermionPseudoScalar3(i1, i2, i3, T3_d, Y_l_160, uL_L, uL_R &
               & , RP0_T, cpl_LLP0_L(i1,i2,i3), cpl_LLP0_R(i1,i2,i3))
     Call CoupFermionScalar3(i1, i2, i3, T3_d, Y_d_160, uD_L, uD_R, RS0_T &
               & , cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
     Call CoupFermionPseudoScalar3(i1, i2, i3, T3_d, Y_d_160, uD_L, uD_R &
               & , RP0_T, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
    End Do
   End Do
  End Do

  Do i1=1,2
   Do i2=1,3
    Do i3=1,3     
     Call CoupCharginoSfermion(i1, i2, i3, g, -0.5_dp, RSneut_T               &
             & , Y_l_160, Zero33C, uL_L, uL_R, U_T, V_T, cpl_CLSn_L(i1,i2,i3) &
             & , cpl_CLSn_R(i1,i2,i3) )
    End Do
    Do i3=1,6
     Call CoupCharginoSfermion(i1, i2, i3, g, 0.5_dp, RSlept_T              &
      & , Y_l_160, Zero33C, id3C, id3C, U_T, V_T, cpl_CNuSl_L(i1,i2,i3)     &
      & , cpl_CNuSl_R(i1,i2,i3) )
     Call CoupCharginoSfermion(i1, i2, i3, g, -0.5_dp, RSup_T               &
           & , Y_D_160, Y_U_160, uD_L, uD_R, U_T, V_T, cpl_CDSu_L(i1,i2,i3) &
           & , cpl_CDSu_R(i1,i2,i3))
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,4
    Do i3=1,6  
     Call CoupNeutralinoSlepton(i1, i2, i3, gp, g, RSlept_T &
      & , uL_L, uL_R, Y_l_160, N_T, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
     Call CoupNeutralinoSdown(i1, i2, i3, gp, g, RSdown_T &
      & , uD_L, uD_R, Y_d_160, N_T, cpl_DNSd_L(i1,i2,i3), cpl_DNSd_R(i1,i2,i3) )
    End Do
    Do i3=1,3  
     Call CoupNeutralinoSneutrino(i1, i2, i3, gp, g, N_T &
           & , RSneut_T, id3C, cpl_NuNSn_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i3=1,6  
    Call CoupGluinoSquark3(gs, phase_Glu_T, i1, i3, RSdown_T, uD_L, uD_R &
           & , cpl_DGSd_L(i1,i3), cpl_DGSd_R(i1,i3) )
   End Do
  End Do
  
  !--------------------------------------------------------
  ! Z-neutralino-neutralino couplings without gauge factor
  !--------------------------------------------------------
  Do i1=1,4
   ZNN(i1,i1) = Abs(N_T(i1,4))**2 - Abs(N_T(i1,3))**2
   Do i2=i1+1,4
    ZNN(i1,i2) = N_T(i1,4) * Conjg( N_T(i2,4) ) - N_T(i1,3) * Conjg( N_T(i2,3) )
    ZNN(i2,i1) = Conjg( ZNN(i1,i2) ) 
   End Do
  End Do
  !-------------------------------------------------------
  ! Z-chargino-chargino couplings without gauge factor
  ! taking only relevant parts
  !--------------------------------------------------------
  Do i1=1,2
   ZUU(i1,i1) = Abs(U_T(i1,1))**2
   ZVV(i1,i1) = Abs(V_T(i1,1))**2
   Do i2=i1+1,2
    ZUU(i1,i2) = U_T(i1,1) * Conjg( U_T(i2,1) )
    ZUU(i2,i1) = Conjg( ZUU(i1,i2) ) 
    ZVV(i1,i2) = V_T(i1,1) * Conjg( V_T(i2,1) )
    ZVV(i2,i1) = Conjg( ZVV(i1,i2) ) 
   End Do
  End Do

  !-------------------------------------------------------
  ! Z-squark-squark couplings without gauge factor
  ! taking only relevant parts
  !--------------------------------------------------------
  Do i1=1,6
   ZSuSuL(i1,i1) = RSup_T(i1,1) * Conjg( RSup_T(i1,1) ) &
               & + RSup_T(i1,2) * Conjg( RSup_T(i1,2) ) &
               & + RSup_T(i1,3) * Conjg( RSup_T(i1,3) )
   ZSuSuR(i1,i1) = RSup_T(i1,4) * Conjg( RSup_T(i1,4) ) &
               & + RSup_T(i1,5) * Conjg( RSup_T(i1,5) ) &
               & + RSup_T(i1,6) * Conjg( RSup_T(i1,6) )
   ZSdSdL(i1,i1) = RSdown_T(i1,1) * Conjg( RSdown_T(i1,1) ) &
               & + RSdown_T(i1,2) * Conjg( RSdown_T(i1,2) ) &
               & + RSdown_T(i1,3) * Conjg( RSdown_T(i1,3) )
   ZSdSdR(i1,i1) = RSdown_T(i1,4) * Conjg( RSdown_T(i1,4) ) &
               & + RSdown_T(i1,5) * Conjg( RSdown_T(i1,5) ) &
               & + RSdown_T(i1,6) * Conjg( RSdown_T(i1,6) )
   Do i2=i1+1,6
    ZSuSuL(i1,i2) = RSup_T(i1,1) * Conjg( RSup_T(i2,1) ) &
                & + RSup_T(i1,2) * Conjg( RSup_T(i2,2) ) &
                & + RSup_T(i1,3) * Conjg( RSup_T(i2,3) )
    ZSuSuR(i1,i2) = RSup_T(i1,4) * Conjg( RSup_T(i2,4) ) &
                & + RSup_T(i1,5) * Conjg( RSup_T(i2,5) ) &
                & + RSup_T(i1,6) * Conjg( RSup_T(i2,6) )
    ZSdSdL(i1,i2) = RSdown_T(i1,1) * Conjg( RSdown_T(i2,1) ) &
                & + RSdown_T(i1,2) * Conjg( RSdown_T(i2,2) ) &
                & + RSdown_T(i1,3) * Conjg( RSdown_T(i2,3) )
    ZSdSdR(i1,i2) = RSdown_T(i1,4) * Conjg( RSdown_T(i2,4) ) &
                & + RSdown_T(i1,5) * Conjg( RSdown_T(i2,5) ) &
                & + RSdown_T(i1,6) * Conjg( RSdown_T(i2,6) )
    ZSuSuL(i2,i1) = Conjg( ZSuSuL(i1,i2) ) 
    ZSuSuR(i2,i1) = Conjg( ZSuSuR(i1,i2) ) 
    ZSdSdL(i2,i1) = Conjg( ZSdSdL(i1,i2) ) 
    ZSdSdR(i2,i1) = Conjg( ZSdSdR(i1,i2) ) 
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
    Call CoupCharginoPseudoScalar(i1, i2, i3, U_T, V_T, RP0_T, g  &
                     & , cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
    Call CoupCharginoScalar(i1, i2, i3, U_T, V_T, RS0_T, g  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  Do i1=1,2
   Do i2=1,6
    Do i3=1,6
     Call CoupScalarSfermion3(i1, i2, i3, RS0_T, -0.5_dp, e_d, Y_d_160   &
      & , RSdown_T, T_d_160, mu_160, vevSM, gp, g, cpl_S0SdSd(i1,i2,i3) )
     Call CoupScalarSfermion3(i1, i2, i3, RS0_T, 0.5_dp, e_u, Y_u_160  &
      & , RSup_T, T_u_160, mu_160, vevSM, gp, g, cpl_S0SuSu(i1,i2,i3) )
    End Do
   End Do
  End Do


  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  bi(1) = mu_160

  Do i1=1,2
   Do i2=1,6
    Do i3=1,6
     Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0_T, -0.5_dp, Y_d_160 &
      & , RSdown_T, T_d_160, bi, cpl_P0SdSd(i1,i2,i3) )
     Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0_T, 0.5_dp, Y_u_160  &
      & , RSup_T, T_u_160, bi, cpl_P0SuSu(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,2
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, N_T, RP0_T, gp, g, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,2
    Call CoupNeutralinoScalar(i1, i2, i3, N_T, RS0_T, gp, g, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  Do i1=1,2
   Call CoupScalarW(i1, g, vevSM, RS0_T, cpl_S0WW(i1))
   Do i2=1,2
    Call CoupChargedScalarPseudoscalarW(i1, i2, g, RSpm_T, RP0_T, cpl_CsP0W(i1,i2))
    Call CoupChargedScalarScalarW(i1, i2, g, RSpm_T, RS0_T, cpl_CsS0W(i1,i2))
    Do i3=1,2
     Call CoupChargedScalarScalar3(i1, i2, i3, RSpm_T, RS0_T, vevSM, gp, g &
                                & , cpl_CsCsS0(i1,i2,i3))
    End Do
   End Do
  End Do

   Call Switch_to_superCKM(Y_d_160, Y_u_160, T_d_160, T_u_160, M2_D_160   &
     & , M2_Q_160, M2_U_160, T_d_s, T_u_s, M2_D_s, M2_Q_s, M2_U_s, .False. )

   Call QuarkMasses_and_PhaseShifts(Y_d_160, Y_u_160, vevSM, mf_d_in, uD_L &
                                   & , uD_R, mf_u_in, uU_L, uU_R, CKM_160)
   Y_d_s = Matmul(Conjg(uD_L),Y_d_160)
   Y_d_s = Matmul(Y_d_s,Transpose(Conjg(uD_R)))
   Y_u_s = Matmul(Conjg(uU_L),Y_u_160)
   Y_u_s = Matmul(Y_u_s,Transpose(Conjg(uU_R)))
   Call Chop(Y_d_s)
   Call Chop(Y_u_s)
   Y_l_s = Y_l_160
   T_l_s = T_l_160

   Call Sigma_SM_chirally_enhanced(gi_160, vevSM, mu_160, CKM_160, Y_l_s &
      & , Y_d_s, Y_u_s, Mi_160, T_l_s, T_d_s, T_u_s, M2_E_160            &
      & , M2_L_160, M2_D_s, M2_Q_s, M2_U_s, epsD, epsL, epsD_FC, kont    &
      & , RS0_T, RP0_T, mf_d_mt, cpl_DDS0_1L_R, cpl_DDP0_1L_R)

   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      cpl_DDS0_1L_L(i1,i2,i3) = Conjg(cpl_DDS0_1L_R(i2,i1,i3))
      cpl_DDP0_1L_L(i1,i2,i3)= Conjg(cpl_DDP0_1L_R(i2,i1,i3))
     End Do
    End Do
   End Do
   !---------------------------------
   ! back to the original basis
   !---------------------------------
   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      cpl_DDS0_1L_R(i1,i2,i3) = 0._dp
      cpl_DDP0_1L_R(i1,i2,i3) = 0._dp
      Do i4=1,3
       Do i5=1,3
        cpl_DDS0_1L_R(i1,i2,i3) = cpl_DDS0_1L_R(i1,i2,i3) + cpl_DDS0_1L_L(i4,i5,i3) &
                               & * Conjg(uD_L(i2,i5)) * Conjg(uD_R(i1,i4))
        cpl_DDP0_1L_R(i1,i2,i3) = cpl_DDP0_1L_R(i1,i2,i3) + cpl_DDP0_1L_L(i4,i5,i3) &
                               & * Conjg(uD_L(i2,i5)) * Conjg(uD_R(i1,i4))
       End Do
      End Do
     End Do
    End Do
   End Do
   cpl_DDS0_1L_L = cpl_DDS0_1L_R
   cpl_DDP0_1L_L = cpl_DDP0_1L_R
   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      cpl_DDS0_1L_R(i1,i2,i3) = Conjg(cpl_DDS0_1L_L(i2,i1,i3))
      cpl_DDP0_1L_R(i1,i2,i3)= Conjg(cpl_DDP0_1L_L(i2,i1,i3))
     End Do
    End Do
   End Do

  !---------------------------------------
  ! Wilson coefficients
  !---------------------------------------
  Call Calculate_Wilson_Coeff_MSSM(mf_d_mt, mf_u_mt, mf_l_mt, mW, mSpm2, mC_T  &
   & , mSup2_T, mSdown2_T, mSneut2_T, mSlept2_T, mglu_T, mN_T, mS02_T, mP02_T  &
   & , gi_160, vevSM, CKM_160, cpl_uWd, cpl_CsDU_L, cpl_CsDU_R, cpl_CDSu_L     &
   & , cpl_CDSu_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R, cpl_DDS0_L  &
   & , cpl_DDS0_R, cpl_DDP0_L, cpl_DDP0_R, cpl_LLS0_L, cpl_LLS0_R, cpl_LLP0_L  &
   & , cpl_LLP0_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_R &
   & , cpl_NuNSn_R, ZNN, ZUU, ZVV, ZSdSdL, ZSdSdR, ZSuSuL, ZSuSuR              &
   & , cpl_DDS0_1L_L, cpl_DDS0_1L_R, cpl_DDP0_1L_L, cpl_DDP0_1L_R, 160._dp, 2  &
   & , WC_c7, WC_c7p, WC_c8, WC_c8p, WC_c9, WC_c9p, WC_c10, WC_c10p            &
   & , WC_c11, WC_c11p)
  l_wc(2) = .True. ! print out values of Wilson coefficient
  !---------------------------------------
  ! BR(b-> s gamma)
  !---------------------------------------
  Call B_to_Q_Gamma(2, CKM, WC_c7(2,3,2,:), WC_c7p(2,3,2,:)       &
    & , WC_c8(2,3,2,:), WC_c8p(2,3,2,:), BRbtosgamma, i_scheme=0  &
    & , NNLO_SM_in=3.2_dp)
  !-------------------------
  ! b -> s l+ l-
  !-------------------------
  Call BToSLL(WC_c7(2,3,2,:), WC_c7p(2,3,2,:), WC_c8(2,3,2,:) &
    & , WC_c8p(2,3,2,:), WC_c9(2,3,2,:,:), WC_c9p(2,3,2,:,:)  &
    & , WC_c10(2,3,2,:,:), WC_c10p(2,3,2,:,:) &
    & , BtoSEE, BrBtoSLL)
  !---------------------------------
  ! b -> s nu nu, no QCD corrections
  !---------------------------------
  Call B_To_SNuNu(CKM, WC_c11(2,3,2,:,1), WC_c11p(2,3,2,:,1), BtoSNuNu)

  Do i1=1,3 
   !-------------------------
   ! B_d -> l+ l-
   !-------------------------
   Call Bq_to_ll(1, i1, mf_d, mf_l, CKM_160, WC_c10(2,3,1,i1,1)              &
       & ,  WC_c10p(2,3,1,i1,1), WC_2d2l_CS(2,3,1,i1), WC_2d2l_CSp(2,3,1,i1) &
       & , WC_2d2l_CP(2,3,1,i1), WC_2d2l_CPp(2,3,1,i1), Bd_ll(i1) )
   !-------------------------
   ! B_s -> mu+ mu-
   !-------------------------
   Call Bq_to_ll(2, i1, mf_d, mf_l, CKM_160, WC_c10(2,3,2,i1,1)              &
       & ,  WC_c10p(2,3,2,i1,1), WC_2d2l_CS(2,3,2,i1), WC_2d2l_CSp(2,3,2,i1) &
       & , WC_2d2l_CP(2,3,2,i1), WC_2d2l_CPp(2,3,2,i1), Bs_ll(i1) )
  End Do
  !-------------------
  ! B^-_u -> tau nu
  !-------------------
  BR_Bu_TauNu = Bm_to_l_nu(3,1, mSpm2(2), tanb_160, RSpm_T, Y_d_160, uU_L &
              &           , uD_R , Y_l_160, vevSM)
  R_Bu_TauNu = Bm_to_l_nu(3,1, mSpm2(2), tanb_160, RSpm_T, Y_d_160, uU_L  &
              &           , uD_R , Y_l_160, vevSM, .True.)
  !-------------------
  ! Delta(M_Bd)
  !-------------------
  Call Delta_MB(1, mf_u_mt(3), mW2, CKM, WC_4d_VLL(2,3,1), WC_4d_VRR(2,3,1)  &
      & , WC_4d_LR1(2,3,1), WC_4d_LR2(2,3,1), WC_4d_SLL1(2,3,1) &
      & , WC_4d_SLL2(2,3,1), WC_4d_SRR1(2,3,1), WC_4d_SRR2(2,3,1), DeltaMBd)
  !-------------------
  ! Delta(M_Bs)
  !-------------------
  Call Delta_MB(2, mf_u_mt(3), mW2, CKM, WC_4d_VLL(2,3,2), WC_4d_VRR(2,3,2)  &
      & , WC_4d_LR1(2,3,2), WC_4d_LR2(2,3,2), WC_4d_SLL1(2,3,2) &
      & , WC_4d_SLL2(2,3,2), WC_4d_SRR1(2,3,2), WC_4d_SRR2(2,3,2), DeltaMBs)
  ! conversion to pico-seconds
  DeltaMBd = 1.e-12_dp*DeltaMBd/hbar
  DeltaMBs = 1.e-12_dp*DeltaMBs/hbar

  !------------------------
  ! K -> pi nu nu
  !------------------------
   Call K_To_PiNuNu(CKM, WC_c11(2,2,1,:,1), WC_c11p(2,2,1,:,1) &
    & , K0toPi0NuNu, KptoPipNuNu)
  !------------------------
  ! epsilon_K 
  !------------------------
   Call epsilon_K(mf_d(1), mf_d(2), mf_u, mW2, CKM, WC_4d_VLL(2,2,1)  &
      & , WC_4d_VRR(2,2,1), WC_4d_LR1(2,2,1), WC_4d_LR2(2,2,1) &
      & , WC_4d_SLL1(2,2,1), WC_4d_SLL2(2,2,1), WC_4d_SRR1(2,2,1)     &
      & , WC_4d_SRR2(2,2,1), epsK, epsK_SM, DMK )
   DeltaMK = DMK(2)

  !-------------------------------------
  ! calculate running masses at m_Z
  ! needs to be improved
  !-------------------------------------
  sinW2 = 1._dp - mW2 / mZ2
  vev2 =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
  vevSM(1) = vev2 / Sqrt(1._dp + tanb_mZ**2)
  vevSM(2) = tanb_mZ * vevSM(1)

  gp = gi_mZ(1)
  g = gi_mZ(2)
  gs = gi_mZ(3)
  mZ2_run = (gp**2+g**2)*0.25_dp*(vevSM(1)**2+vevSM(2)**2)
  mW2_run = g**2*0.25_dp*(vevSM(1)**2+vevSM(2)**2)

  abs_mu2 = (M2_H_mZ(2) * tanb_mz**2 - M2_H_mZ(1) ) / (1._dp-tanb_mZ**2)&
          &  - 0.5_dp * mZ2_run
  mu_mZ = phase_mu * Sqrt(abs_mu2)
  B_mZ = (M2_H_mZ(1) + M2_H_mZ(2) + 2._dp * Abs_Mu2) * tanb_mz / (1+tanb_mz**2)

  Call TreeMassesMSSM2(gp, g, vevSM, Mi_mZ(1), Mi_mZ(2)   &
     & , Mi_mZ(3), mu_mZ, B_mZ, tanb_mZ, M2_E_mZ, M2_L_mZ, T_l_mZ, Y_l_mZ    &
     & , M2_D_mZ, M2_U_mZ, M2_Q_mZ, T_d_mZ, T_u_mZ, Y_d_mZ, Y_u_mZ           &
     & , uU_L, uU_R ,uD_L, uD_R, uL_L, uL_R                                  &
     & , mGlu_T, Phase_Glu_T, mC_T, mC2_T, U_T, V_T, mN_T, mN2_T, N_T        &
     & , mSneut_T, mSneut2_T, Rsneut_T, mSlept_T, mSlept2_T      &
     & , RSlept_T, mSdown_T, mSdown2_T, RSdown_T, mSup_T, mSup2_T, RSup_T  &
     & , mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T, mSpm_T, mSpm2_T, RSpm_T &
     & , mZ2_run, mW2_run, .True., kont, .False., .False., mf_u_Q)

  If (kont.Ne.0) Then 
   Iname = Iname - 1
   kont = -701
   Call AddError(701)
   Return
  End If 

  !-------------------------------
  ! setting renormalisation scale
  !-------------------------------
  work = SetRenormalizationScale( mZ**2 )

  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cosW = g / Sqrt(gp**2 + g**2)
  sinW2 = 1._dp - cosW**2

  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,2
   Do i2=1,2
    Call CoupCharginoZ(i1, i2, U_T, V_T, g, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,4
   Do i2=1,4
    Call CoupNeutralinoZ(i1, i2, N_T, g, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  Do i1=1,2
   Do i2=1,3
    Do i3=1,3     
     Call CoupCharginoSfermion(i1, i2, i3, g, -0.5_dp, RSneut_T  &
             & , Y_l_mZ, Zero33C, uL_L, uL_R, U_T, V_T, cpl_CLSn_L(i1,i2,i3) &
             & , cpl_CLSn_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,4
    Do i3=1,6  
     Call CoupNeutralinoSlepton(i1, i2, i3, gp, g, RSlept_T &
         & , uL_L, uL_R, Y_l_mZ, N_T, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------------------------------------------
  ! leptonic electric dipole moments
  !------------------------------------------------------------------
  Call Lepton_EDM3(1, mN_T, mSlept2_T, cpl_LNSl_L, cpl_LNSl_R          &
                & , mC_T, mSneut2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_e   )
  Call Lepton_EDM3(2, mN_T, mSlept2_T, cpl_LNSl_L, cpl_LNSl_R          &
                & , mC_T, mSneut2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_mu  )
  Call Lepton_EDM3(3, mN_T, mSlept2_T, cpl_LNSl_L, cpl_LNSl_R          &
                & , mC_T, mSneut2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_tau )
  d_e = EDM_e(3)
  d_mu = EDM_mu(3)
  d_tau = EDM_tau(3)
  !------------------------------------------------------------------
  ! leptonic anomalous magnetic moments
  !------------------------------------------------------------------
  Call Gminus2(1, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_e, GenMix)
  Call Gminus2(2, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_mu, GenMix)
  Call Gminus2(3, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_tau, GenMix)
  !------------------------------------------------------------------
  ! rare decays of leptons: l -> l' gamma
  !------------------------------------------------------------------
  BrMutoEGamma = 0._dp
  BrTautoEGamma = 0._dp
  BrTautoMuGamma = 0._dp
  If (GenMix) Then
   Call LtoLpGamma(2, 1, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T, mN_T &
                  &, cpl_LNSl_L, cpl_LNSl_R, GMutoEGamma, BrMutoEGamma)
   Call LtoLpGamma(3, 1, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T, mN_T &
                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoEGamma, BrTautoEGamma)
   Call LtoLpGamma(3, 2, mSneut2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlept2_T, mN_T &
                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoMuGamma, BrTautoMuGamma)
  End If
  !------------------------------------------------------------------
  ! rare decays of leptons: l -> 3 l' 
  !------------------------------------------------------------------
  BrMu3E = 0._dp
  BrTau3E = 0._dp
  BrTau3Mu = 0._dp
  gU1 = gi_mZ(1)
  gSU2 = gi_mZ(2)
  If (GenMix) Then
   Call BR_lj_to_3li(2, 1, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlept2_T, RSlept_T   &
          & , mN_T, N_T, mSneut2_T, RSneut_T, mC_T, U_T, V_T, mS02, RS0_T     &
          & , mP02, RP0_T, T_l_mZ, mu_mZ, vevSM, BrMu3e)
   Call BR_lj_to_3li(3, 1, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlept2_T, RSlept_T   &
          & , mN_T, N_T, mSneut2_T, RSneut_T, mC_T, U_T, V_T, mS02, RS0_T     &
          & , mP02, RP0_T, T_l_mZ, mu_mZ, vevSM, BrTau3e)
   Call BR_lj_to_3li(3, 2, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlept2_T, RSlept_T   &
          & , mN_T, N_T, mSneut2_T, RSneut_T, mC_T, U_T, V_T, mS02, RS0_T     &
          & , mP02, RP0_T, T_l_mZ, mu_mZ, vevSM, BrTau3Mu)
  End If
  !-------------
  ! delta(rho)
  !-------------
  rho_parameter = DeltaRho(mZ2, mW2, mP02_T, RP0_T, mSneut2_T, RSneut_T  &
                &         , mSlept2_T, RSlept_T, mSup2_T, RSup_T, mSdown2_T    &
                &         , RSdown_T, mC_T, U_T, V_T, mN_T, N_T)
  !----------------------------------
  ! rare Z-boson decays into leptons
  !----------------------------------
  BR_Z_e_mu = 0._dp
  BR_Z_e_tau = 0._dp
  BR_Z_mu_tau = 0._dp

  If (GenMix) Then

   cpl_SnSnZ = 0._dp
   Call CoupSneutrinoZ(gSU2, sinW2, cpl_SnSnZ(1,1))
   cpl_SnSnZ(2,2) = cpl_SnSnZ(1,1)
   cpl_SnSnZ(3,3) = cpl_SnSnZ(1,1)

   Do i1=1,6
    Do i2=1,6
     Call CoupSleptonZ(i1, i2, gSU2, sinW2, RSlept_T, cpl_SlSlZ(i1,i2))
    End Do
   End Do

   Call CoupFermionZ(-0.5_dp,-1._dp, gSU2,sinW2,cpl_LLZ_L,cpl_LLZ_R)

   Call ZtoLiLj(1, 2, .False., cpl_LLZ_L, cpl_LLZ_R  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneut2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlept2_T, cpl_SlSlZ   &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, BR_Z_e_mu)

   Call ZtoLiLj(1, 3, .True., cpl_LLZ_L, cpl_LLZ_R  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneut2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlept2_T, cpl_SlSlZ   &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, BR_Z_e_tau)

   Call ZtoLiLj(2, 3, .True., cpl_LLZ_L, cpl_LLZ_R  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneut2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlept2_T, cpl_SlSlZ   &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, BR_Z_mu_tau)

  End If ! GenerationMixing

  !-------------------------------------------------------------------------
  ! neutrino masses and mixings, if the dim-5 operator has non-zero entries
  ! or if neutrino mass matrix is given from the outside
  !-------------------------------------------------------------------------
  If ((Maxval(Abs(MnuL5)).Gt.0._dp)) Then
   MnuL5 = MnuL5 * vevSM(2)**2
   Call NeutrinoMass_1L(MnuL5, gi(1), gi(2), Y_l_mZ, mC2_T, U_T, V_T, mN2_T &
         & , N_T, mSlept2_T, RSlept_T, mSneut2_T, Rsneut_T, mf_nu, Unu, kont)
  Elseif ((Maxval(Abs(MatNu)).Gt.0._dp).And.(.Not.fake_m_nu)) Then
   Call NeutrinoMasses(MatNu, mf_nu, Unu, kont)
  End If

  If (WriteDetails) Close(45)

  !---------------------------------
  ! re-setting renormalisation scale
  !---------------------------------
  GenerationMixing = GenMix_save
  mudim_old = SetRenormalizationScale( mudim_old )

  Iname = Iname - 1

 End Subroutine Low_Energy_Constraints_MSSM


! Interface LtoLpGamma
 Subroutine LtoLpGammaMSSM(i, j, mSn2, mC, cpl_CLSn_L, cpl_CLSn_R, mSl2, mN &
                   &, cpl_LNSl_L, cpl_LNSl_R, width, Br, Amplitude)
 !-----------------------------------------------------------------------
 ! Calculates the width l -> l' gamma using the formula of
 ! J.Hisano et al., PRD 53, 2442 (1996)
 ! input:
 !  i ................... index of decaying lepton: 2 -> muon, 3 -> tau
 !  j ................... index of final state lepton: 1 -> electron, 2 -> muon
 !  mSn2(i) ............. sneutrino masses squared
 !  mC(i) ............... chargino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino - lepton -sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino - lepton -sneutrino coupling
 !  mSl2(i) ............. slepton masses squared
 !  mN(i) ............... neutralino masses 
 !  cpl_LNSl_L(i,j,k) ... left lepton - neutralino - slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton - neutralino - slepton coupling
 ! output:
 !  width ............... decay width in GeV
 !  Br .................. branching ratio, optional
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mSn2(3), mSl2(6), mC(:), mN(:)
  Real(dp), Intent(out) :: width
  Real(dp), Intent(out), Optional :: Br
  Complex(dp), Intent(out), Optional :: Amplitude(2,2)
  Complex(dp), Intent(in) :: cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:) &
                         & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)

  Integer :: i1, i2, n_C, n_N
  Real(dp) :: ratio, mN2, mC2, part(2)
  Complex(dp) :: AL, AR, amp(2,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'LtoLpGammaMSSM'
       
  n_C = Size(mC)
  n_N = Size(mN)

  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! neutralino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_N
   mN2 = mN(i1)**2
   Do i2 = 1, 6
    ratio = mSl2(i2) / mN2
    part(1) = 2._dp * mf_l(i) * F1(ratio) / mN2
    part(2) = 2._dp * F4(ratio) / mN(i1)
    AL = AL + cpl_LNSl_L(j,i1,i2) * ( Conjg( cpl_LNSl_L(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_R(i,i1,i2) ) * part(2) )
    AR = AR + cpl_LNSl_R(j,i1,i2) * ( Conjg( cpl_LNSl_R(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_L(i,i1,i2) ) * part(2) )
    
   End Do
  End Do
  amp(1,1) = AL
  amp(1,2) = AR
  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! chargino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_C
   mC2 = mC(i1)**2
   Do i2 = 1, 3
    ratio = mSn2(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CLSn_L(i1,j,i2) * ( Conjg( cpl_CLSn_L(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CLSn_R(i1,j,i2) * ( Conjg( cpl_CLSn_R(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_L(i1,i,i2) ) * part(2) )
   End Do
  End Do
  amp(2,1) = AL
  amp(2,2) = AR

  AL = oo32pi2 * (amp(1,1)+amp(2,1))
  AR = oo32pi2 * (amp(1,2)+amp(2,2))

  width = 0.25_dp * Alpha * mf_l(i)**3 * (Abs(AL)**2 + Abs(AR)**2)

  If (Present(Br)) Then
   If (i.Eq.2) Br = width / (width+GammaMu)
   If (i.Eq.3) Br = width / (width+GammaTau)
  End If

  If (Present(Amplitude))  Amplitude = oo32pi2 * amp

  Iname = Iname - 1

 End Subroutine LtoLpGammaMSSM

 Subroutine LtoLpGammaEps3(i, j, mC, mP02, cpl_CCP0_L, cpl_CCP0_R   &
                         &, mS02, cpl_CCS0_L, cpl_CCS0_R            &
                         &, mN, mSpm2, cpl_SpmCN_L, cpl_SpmCN_R     &
                         &, mUsquark2, cpl_CDSu_L, cpl_CDSu_R       &
                         &, mDsquark2, cpl_CUSd_L, cpl_CUSd_R       &
                         &, width, Br, GenerationMixing, Amplitude)
 !-----------------------------------------------------------------------
 ! Calculates the width l -> l' gamma based on the formula of
 ! J.Hisano et al., PRD 53, 2442 (1996), this is an extension to the
 ! 3-generation R-parity model
 ! input:
 !  i ................... index of decaying lepton: 2 -> muon, 3 -> tau
 !  j ................... index of final state lepton: 1 -> electron, 2 -> muon
 !  mC(i) ............... chargino masses
 !  mP02(i) ............. pseudoscalar masses squared
 !  cpl_CCP0_L(i,j,k) ... left chargino - chargino - pseudoscalar coupling
 !  cpl_CCP0_R(i,j,k) ... right chargino - chargino - pseudoscalar coupling
 !  mS02(i) ............. scalar masses squared
 !  cpl_CCS0_L(i,j,k) ... left chargino - chargino - scalar coupling
 !  cpl_CCS0_R(i,j,k) ... right chargino - chargino - scalar coupling
 !  mN(i) ............... neutralino masses 
 !  mSpm2(i) ............. charged scalar masses squared
 !  cpl_SpmCN_L(i,j,k) .. left charged scalar - chargino - neutralino  coupling
 ! cpl_SpmCN_R(i,j,k) .. right charged scalar - chargino - neutralino  coupling
 !  GenerationMixing ... generation mixing between quarks/squarks is taken into
 !                       account if .TRUE., optional
 !  Amplitude(i,j) ..... contains the neutralino the chargino,, the u-quark
 !                       and the d-quark contribution to the amplitude if
 !                       present, optional, the second index specifies L,R
 ! output:
 !  width ............... decay width in GeV
 !  Br .................. branching ratio, optional
 ! written by Werner Porod, 14.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mP02(:), mS02(:), mSpm2(:), mC(:), mN(:)   &
         &, mUsquark2(6), mDsquark2(6)
  Real(dp), Intent(out) :: width
  Real(dp), Intent(out), Optional :: Br
  Complex(dp), Intent(in) :: cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:)    &
                         & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:)    &
                         & , cpl_SpmCN_L(:,:,:), cpl_SpmCN_R(:,:,:)  &
                         & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)    &
                         & , cpl_CUSd_L(:,:,:), cpl_CUSd_R(:,:,:)
  Complex(dp), Intent(out), Optional :: Amplitude(4,2)
  Logical, Intent(in), Optional :: GenerationMixing

  Integer :: i1, i2, n_C, n_N, n_P0, n_S0, n_Spm, i_gen
  Real(dp) :: ratio, mN2, mC2, part(2), mF2
  Complex(dp) :: AL, AR, amp(4,2)
  Logical :: generation

  Iname = Iname + 1
  NameOfUnit(Iname) = 'LtoLpGammaEps3'
       
  n_C = Size(mC)
  n_N = Size(mN)
  n_P0 = Size(mP02)
  n_S0 = Size(mS02)
  n_Spm = Size(mSpm2)

  AL = ZeroC
  AR = ZeroC
  amp = ZeroC

  If (Present(GenerationMixing) ) Then
   generation = GenerationMixing
  Else
   generation = .False.
  End If

  !-------------------------------------
  ! neutralino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_N
   mN2 = mN(i1)**2
   Do i2 = 1, n_Spm
    ratio = mSpm2(i2) / mN2
    part(1) = 2._dp * mf_l(i) * F1(ratio) / mN2
    part(2) = 2._dp * F4(ratio) / mN(i1)
    AL = AL + cpl_SpmCN_L(i2,j,i1) * ( Conjg(cpl_SpmCN_L(i2,i,i1) ) * part(1) &
       &                            + Conjg( cpl_SpmCN_R(i2,i,i1) ) * part(2) )
    AR = AR + cpl_SpmCN_R(i2,j,i1) * ( Conjg(cpl_SpmCN_R(i2,i,i1) ) * part(1) &
       &                            + Conjg( cpl_SpmCN_L(i2,i,i1) ) * part(2) )
   End Do
  End Do
  amp(1,1) = AL
  amp(1,2) = AR
  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! chargino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_C
   mC2 = mC(i1)**2
   Do i2 = 1, n_P0
    ratio = mP02(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CCP0_L(i1,j,i2) * ( Conjg( cpl_CCP0_L(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCP0_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CCP0_R(i1,j,i2) * ( Conjg( cpl_CCP0_R(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCP0_L(i1,i,i2) ) * part(2) )
   End Do
   Do i2 = 1, n_S0
    ratio = mS02(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CCS0_L(i1,j,i2) * ( Conjg( cpl_CCS0_L(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCS0_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CCS0_R(i1,j,i2) * ( Conjg( cpl_CCS0_R(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCS0_L(i1,i,i2) ) * part(2) )
   End Do
  End Do
  amp(2,1) = AL !- amp(1,1)
  amp(2,2) = AR !- amp(1,2)
  AL = ZeroC
  AR = ZeroC
  !----------------------------------------------------------------------
  ! d-quark contribution to A_L, A_R
  ! note that the colour factor 3 cancels the 1/3 contribution of the
  ! quark charges
  !----------------------------------------------------------------------
  If (generation) Then
   Do i1 = 1, 3
    mF2 = mf_d(i1)**2
    Do i2 = 1,6
     ratio = mUsquark2(i2) / mF2
     part(1) = mf_l(i) * ( -2._dp * F2(ratio) - 4._dp * F1(ratio) ) / mF2
     part(2) = (F3gamma(ratio) + 4._dp * F4(ratio) ) / mf_d(i1)
     AL = AL + cpl_CDSu_L(j,i1,i2) * ( Conjg( cpl_CDSu_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CDSu_R(j,i1,i2) * ( Conjg( cpl_CDSu_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  Else
   Do i1 = 1, 3
    mF2 = mf_d(i1)**2
    i_gen = 2*i1-1
    Do i2 = i_gen,i_gen+1
     ratio = mUsquark2(i2) / mF2
     part(1) = mf_l(i) * (- 2._dp * F2(ratio) - 4._dp * F1(ratio) ) / mF2
     part(2) = (F3gamma(ratio) + 4._dp * F4(ratio) ) / mf_d(i1)
     AL = AL + cpl_CDSu_L(j,i1,i2) * ( Conjg( cpl_CDSu_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CDSu_R(j,i1,i2) * ( Conjg( cpl_CDSu_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  End If
  amp(3,1) = AL !- amp(1,1)
  amp(3,2) = AR !- amp(1,2)
  AL = ZeroC
  AR = ZeroC
  !----------------------------------------------------------------------
  ! u-quark contribution to A_L, A_R
  ! note that the colour factor 3 cancels the 1/3 contribution of the
  ! quark charges
  !----------------------------------------------------------------------
  If (generation) Then
   Do i1 = 1, 3
    mF2 = mf_u(i1)**2
    Do i2 = 1,6
     ratio = mDsquark2(i2) / mF2
     part(1) = mf_l(i) * (4._dp * F2(ratio) + 2._dp * F1(ratio) ) / mF2
     part(2) = (- 2._dp * F3gamma(ratio) - 2._dp * F4(ratio) ) / mf_u(i1)
     AL = AL + cpl_CUSd_L(j,i1,i2) * ( Conjg( cpl_CUSd_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CUSd_R(j,i1,i2) * ( Conjg( cpl_CUSd_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  Else
   Do i1 = 1, 3
    mF2 = mf_u(i1)**2
    i_gen = 2*i1-1
    Do i2 = i_gen,i_gen+1
     ratio = mDsquark2(i2) / mF2
     part(1) = mf_l(i) * (4._dp * F2(ratio) + 2._dp * F4(ratio) ) / mF2
     part(2) = (- 2._dp * F3gamma(ratio) - 2._dp * F4(ratio) ) / mf_u(i1)
     AL = AL + cpl_CUSd_L(j,i1,i2) * ( Conjg( cpl_CUSd_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CUSd_R(j,i1,i2) * ( Conjg( cpl_CUSd_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  End If
  amp(4,1) = AL 
  amp(4,2) = AR 

  AL = oo32pi2 * Sum(amp(:,1))
  AR = oo32pi2 * Sum(amp(:,2))

  width = 0.25_dp * Alpha * mf_l(i)**3 * (Abs(AL)**2 + Abs(AR)**2)

  If (Present(Br)) Then
   If (i.Eq.2) Br = width / (width+GammaMu)
   If (i.Eq.3) Br = width / (width+GammaTau)
  End If

  If (Present(Amplitude))  Amplitude = oo32pi2 * amp

  Iname = Iname - 1
 
 End Subroutine LtoLpGammaEps3

! Interface LtoLpGamma
 Subroutine LtoLpGammaRPcons(i, j, mSn2, mC, cpl_CLSn_L, cpl_CLSn_R, mSl2, mN  &
                  & , cpl_LNSl_L, cpl_LNSl_R, mW2, mNu, cpl_LNuW_L, cpl_LNuW_R &
                  & , width, Br, Amplitude)
 !-----------------------------------------------------------------------------
 ! Calculates the width l -> l' gamma using the formula of
 ! J.Hisano et al., PRD 53, 2442 (1996)
 ! input:
 !  i ................... index of decaying lepton: 2 -> muon, 3 -> tau
 !  j ................... index of final state lepton: 1 -> electron, 2 -> muon
 !  mSn2(i) ............. sneutrino masses squared
 !  mC(i) ............... chargino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino - lepton -sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino - lepton -sneutrino coupling
 !  mSl2(i) ............. slepton masses squared
 !  mN(i) ............... neutralino masses
 !  cpl_LNSl_L(i,j,k) ... left lepton - neutralino - slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton - neutralino - slepton coupling
 !  mW2(i) .............. W-boson masses squared
 !  mNu(i) .............. neutrino masses
 !  cpl_LNuW_L(i,j,k) ... left lepton - neutrino - W coupling
 !  cpl_LNuW_R(i,j,k) ... right lepton - neutrino - W coupling
 ! output:
 !  width ............... decay width in GeV
 !  Br .................. branching ratio, optional
 !  Amplitude ........... neutralino, chargino and neutrino contribtions to the 
 !                        amplitude, optional
 ! written by Werner Porod, 09.05.2001
 ! 21.09.2011: extending LtoLpGammaMSSM by the W-nu loops, using formulas of
 !             Illakova, Pilaftis, NPB 437 (1995) 491
 ! 22.09.2011: implementing mixing part using formulas by L.Lavoura
 !             hep-ph/0302221 
 !-----------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mSn2(:), mSl2(:), mC(:), mN(:), mW2(:), mNu(:)
  Real(dp), Intent(out) :: width
  Real(dp), Intent(out), Optional :: Br
  Complex(dp), Intent(out), Optional :: Amplitude(3,2)
  Complex(dp), Intent(in) :: cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:) &
                         & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:) &
                         & , cpl_LNuW_L(:,:,:), cpl_LNuW_R(:,:,:)

  Integer :: i1, i2, n_C, n_N, n_Sn, n_Sl, n_W, n_Nu
  Real(dp) :: ratio, mN2, mC2, part(2), mNu2
  Complex(dp) :: AL, AR, amp(3,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'LtoLpGammaRPcons'
       
  n_C = Size(mC)
  n_N = Size(mN)
  n_Nu = Size(mNu)
  n_Sn = Size(mSn2)
  n_Sl = Size(mSl2)
  n_W = Size(mW2)

  AL = ZeroC
  AR = ZeroC
  amp(3,2) = ZeroC
  !-------------------------------------
  ! neutralino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_N
   mN2 = mN(i1)**2
   Do i2 = 1, n_Sl
    ratio = mSl2(i2) / mN2
    part(1) = 2._dp * mf_l(i) * F1(ratio) / mN2
    part(2) = 2._dp * F4(ratio) / mN(i1)
    AL = AL + cpl_LNSl_L(j,i1,i2) * ( Conjg( cpl_LNSl_L(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_R(i,i1,i2) ) * part(2) )
    AR = AR + cpl_LNSl_R(j,i1,i2) * ( Conjg( cpl_LNSl_R(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_L(i,i1,i2) ) * part(2) )
    
   End Do
  End Do
  amp(1,1) = AL
  amp(1,2) = AR
  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! chargino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_C
   mC2 = mC(i1)**2
   Do i2 = 1, n_Sn
    ratio = mSn2(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CLSn_L(i1,j,i2) * ( Conjg( cpl_CLSn_L(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CLSn_R(i1,j,i2) * ( Conjg( cpl_CLSn_R(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_L(i1,i,i2) ) * part(2) )
   End Do
  End Do
  amp(2,1) = AL
  amp(2,2) = AR
  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! Neutrino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_Nu
   mNu2 = mNu(i1)**2
   Do i2 = 1, n_W
    ratio = mNu2 / mW2(i2) 
    part(1) = 3._dp * mf_l(i) * ratio * F2(ratio) / mW2(i2)
    part(2) = 3._dp * mNu(i1) * F4(ratio) / mW2(i2)
    AL = AL + cpl_LNuW_L(j,i1,i2) * ( Conjg( cpl_LNuW_L(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNuW_R(i,i1,i2) ) * part(2) )
    AR = AR + cpl_LNuW_R(j,i1,i2) * ( Conjg( cpl_LNuW_R(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNuW_L(i,i1,i2) ) * part(2) )
   End Do
  End Do
  amp(3,1) = AL
  amp(3,2) = AR

  AL = oo32pi2 * Sum(amp(:,1))
  AR = oo32pi2 * Sum(amp(:,2))

  width = 0.25_dp * Alpha * mf_l(i)**3 * (Abs(AL)**2 + Abs(AR)**2)

  If (Present(Br)) Then
   If (i.Eq.2) Br = width / GammaMu
   If (i.Eq.3) Br = width / GammaTau
  End If

  If (Present(Amplitude)) Amplitude = oo32pi2 * amp

  Iname = Iname - 1

 End Subroutine LtoLpGammaRPcons
! Interface LtoLpGamma

 Subroutine ZtoLiLj(i, j, UseSavedLoopFunctions, cpl_LLZ_L, cpl_LLZ_R  &
       & , mC, mC2, cpl_CLSn_L, cpl_CLSn_R, mSnu2, cpl_SnSnZ           &
       & , mN, mN2, cpl_LNSl_L, cpl_LNSl_R, mSle2, cpl_SlSlZ           &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, BR) 
 !-------------------------------------------------------------------------
 ! in this subroutine the decay branching ratio for the decay
 ! Z -> l_i  l_j is calculated.
 ! The formulas of X.Bi et al, PRD63, 096008, (2001) are used.
 ! Input: i,j ....................... indices of the leptons
 !         UseSavedLoopFunctions .... the Loopfunctions calculated in a
 !                                    previous run are used if set .true.
 ! written by Werner Porod, 17.10.01 
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Logical, Intent(in) :: UseSavedLoopFunctions
  Real(dp), Intent(in) :: mC(2), mC2(2), mSnu2(3), mN(4), mN2(4), mSle2(6) &
     & , cpl_LLZ_L, cpl_LLZ_R
  Complex(dp), Intent(in) ::  cpl_SnSnZ(:,:), cpl_CLSn_L(:,:,:)            &
     & , cpl_CLSn_R(:,:,:), cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)           &
     & , cpl_SlSlZ(:,:), cpl_CCZ_L(:,:), cpl_CCZ_R(:,:), cpl_NNZ_L(:,:)    &
     & , cpl_NNZ_R(:,:)
  Real(dp), Intent(out) :: BR

  !--------------------------
  ! to save calculation time
  !--------------------------
  Complex(dp), Save :: C00_C_i_Snu_jk(2,3), C00_N_i_Sle_jk(4,6,6)   &
     & , C00_Snu_i_C_jk(3,2,2), C00_Sle_i_N_jk(6,4,4)               &
     & , C0_C_i_Snu_jk(2,3), C0_N_i_Sle_jk(4,6,6)                 &
     & , C1_C_i_Snu_jk(2,3), C1_N_i_Sle_jk(4,6,6)                 &
     & , C2_C_i_Snu_jk(2,3), C2_N_i_Sle_jk(4,6,6)                 &
     & , C0_Snu_i_C_jk(3,2,2), C0_Sle_i_N_jk(6,4,4)                 &
     & , C1_Snu_i_C_jk(3,2,2), C1_Sle_i_N_jk(6,4,4)                 &
     & , C2_Snu_i_C_jk(3,2,2), C2_Sle_i_N_jk(6,4,4)                 &
     & , B0_Z_CC(2,2), B0_0_CSn(2,3), B0_Z_NN(4,4), B0_0_NSl(4,6)   &
     & , B1_Z_CC(2,2), B1_0_CSn(2,3), B1_Z_NN(4,4), B1_0_Nsl(4,6)!     &

  Integer :: i1, i2, i3
  Real(dp) :: mZ2a, gam
  Complex(dp) :: Ai(4)

  Iname = Iname + 1
  NameOfUnit(Iname) = "ZtoLiLj"

  !------------------------------------------------------
  ! Calculate first the loop functions if necessary
  !------------------------------------------------------
  If (.Not.UseSavedLoopFunctions) Then
   mZ2a = 0._dp ! the leptons are on-shell
   Do i1=1,2
    B0_Z_CC(i1,i1) = B0(mZ2, mC2(i1), mC2(i1) )
    B1_Z_CC(i1,i1) = B1(mZ2, mC2(i1), mC2(i1) )
    Do i2=i1+1,2
     B0_Z_CC(i1,i2) = B0(mZ2, mC2(i1), mC2(i2) )
     B1_Z_CC(i1,i2) = B1(mZ2, mC2(i1), mC2(i2) )
     B0_Z_CC(i2,i1) = B0(mZ2, mC2(i2), mC2(i1) )
     B1_Z_CC(i2,i1) = B1(mZ2, mC2(i2), mC2(i1) )
    End Do
    Do i2=1,3
     B0_0_CSn(i1,i2) = B0(mZ2a, mSnu2(i2), mC2(i1) )
     B1_0_CSn(i1,i2) = B1(mZ2a, mSnu2(i2), mC2(i1) )
     C00_C_i_Snu_jk(i1,i2) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C0_C_i_Snu_jk(i1,i2) = &
          & C0( 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C1_C_i_Snu_jk(i1,i2) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C2_C_i_Snu_jk(i1,i2) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C0_Snu_i_C_jk(i2,i1,i1)  = &
          & C0( 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C1_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C2_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C00_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     If (i1.Eq.1) Then 
      C0_Snu_i_C_jk(i2,i1,i1+1)  = &
            & C0( 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C0_Snu_i_C_jk(i2,i1+1,i1) = C0_Snu_i_C_jk(i2,i1,i1+1)
      C1_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C1  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C1_Snu_i_C_jk(i2,i1+1,i1) = C1_Snu_i_C_jk(i2,i1,i1+1)
      C2_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C2  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C2_Snu_i_C_jk(i2,i1+1,i1) = C2_Snu_i_C_jk(i2,i1,i1+1)
      C00_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C00 ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C00_Snu_i_C_jk(i2,i1+1,i1) = C00_Snu_i_C_jk(i2,i1,i1+1)
     End If
    End Do
   End Do
   Do i1=1,4
    B0_Z_NN(i1,i1) = B0(mZ2, mN2(i1), mN2(i1) )
    B1_Z_NN(i1,i1) = B1(mZ2, mN2(i1), mN2(i1) )
    Do i2=i1+1,4
     B0_Z_NN(i1,i2) = B0(mZ2, mN2(i1), mN2(i2) )
     B1_Z_NN(i1,i2) = B1(mZ2, mN2(i1), mN2(i2) )
     B0_Z_NN(i2,i1) = B0(mZ2, mN2(i2), mN2(i1) )
     B1_Z_NN(i2,i1) = B1(mZ2, mN2(i2), mN2(i1) )
    End Do
    Do i2=1,6
     B0_0_NSl(i1,i2) = B0(mZ2a, mSle2(i2), mN2(i1) )
     B1_0_NSl(i1,i2) = B1(mZ2a, mSle2(i2), mN2(i1) )
     C00_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C0_N_i_Sle_jk(i1,i2,i2) = &
          & C0( 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C1_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C2_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     Do i3=i2+1,6
      C00_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C00_N_i_Sle_jk(i1,i3,i2) = C00_N_i_Sle_jk(i1,i2,i3)
      C0_N_i_Sle_jk(i1,i2,i3) = &
          & C0( 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C0_N_i_Sle_jk(i1,i3,i2) = C0_N_i_Sle_jk(i1,i2,i3)
      C1_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C1_N_i_Sle_jk(i1,i3,i2) = C1_N_i_Sle_jk(i1,i2,i3)
      C2_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C2_N_i_Sle_jk(i1,i3,i2) = C2_N_i_Sle_jk(i1,i2,i3)
     End Do
     C0_Sle_i_N_jk(i2,i1,i1) = &
          & C0( 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C1_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C2_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C00_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     Do i3=i1+1,4
      C0_Sle_i_N_jk(i2,i1,i3) = &
          & C0( 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C0_Sle_i_N_jk(i2,i3,i1) = C0_Sle_i_N_jk(i2,i1,i3)
      C1_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C1_Sle_i_N_jk(i2,i3,i1) = C1_Sle_i_N_jk(i2,i1,i3)
      C2_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C2_Sle_i_N_jk(i2,i3,i1) = C2_Sle_i_N_jk(i2,i1,i3)
      C00_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C00_Sle_i_N_jk(i2,i3,i1) = C00_Sle_i_N_jk(i2,i1,i3)
     End Do
    End Do
   End Do
  End If

  !------------------------------------------------------
  ! now the effective couplings
  !------------------------------------------------------
  Ai = 0._dp
  Do i1=1,2
   Do i2=1,3
    Ai(1) = Ai(1)                                                           &
        & + 2._dp * cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_R(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * C00_C_i_Snu_jk(i1,i2)
    Ai(2) = Ai(2)                                                           &
        & + 2._dp * cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_L(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * C00_C_i_Snu_jk(i1,i2)
    Ai(3) = Ai(3)                                                           &
        & + 2._dp * cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_L(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * mC(i1) * ( C0_C_i_Snu_jk(i1,i2)      &
        &                                       + C1_C_i_Snu_jk(i1,i2)      &
        &                                       + C2_C_i_Snu_jk(i1,i2) )
    Ai(4) = Ai(4)                                                           &
        & + 2._dp * cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_R(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * mC(i1) * ( C0_C_i_Snu_jk(i1,i2)      &
        &                                       + C1_C_i_Snu_jk(i1,i2)      &
        &                                       + C2_C_i_Snu_jk(i1,i2) )
    Do i3=1,2
     Ai(1) = Ai(1)                                                           &
      &  - cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_R(i3, j, i2) )            &
      &    * ( cpl_CCZ_R(i1,i3) * ( mSnu2(i2) * C0_Snu_i_C_jk(i2,i1,i3)      &
      &                           + B0_Z_CC(i1,i3)                           &
      &                           - 2._dp * C00_Snu_i_C_jk(i2,i1,i3) )       & 
      &      - mC(i1) *mC(i3) *cpl_CCZ_L(i1,i3) * C0_Snu_i_C_jk(i2,i1,i3) )
     Ai(2) = Ai(2)                                                           &
      &  - cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_L(i3, j, i2) )            &
      &    * ( cpl_CCZ_L(i1,i3) * ( mSnu2(i2) * C0_Snu_i_C_jk(i2,i1,i3)      &
      &                           + B0_Z_CC(i1,i3)                           &
      &                           -  2._dp * C00_Snu_i_C_jk(i2,i1,i3) )      & 
      &      - mC(i1) *mC(i3) *cpl_CCZ_R(i1,i3) * C0_Snu_i_C_jk(i2,i1,i3) )
     Ai(3) = Ai(3)                                                           &
      &  + cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_L(i3, j, i2) )            &
      &    * ( cpl_CCZ_L(i1,i3) * mC(i3) *  C1_Snu_i_C_jk(i2,i1,i3)          &
      &      + cpl_CCZ_R(i1,i3) * mC(i1) *  C2_Snu_i_C_jk(i2,i1,i3)  )
     Ai(4) = Ai(4)                                                           &
      &  + cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_R(i3, j, i2) )            &
      &    * ( cpl_CCZ_R(i1,i3) * mC(i3) *  C1_Snu_i_C_jk(i2,i1,i3)          &
      &      + cpl_CCZ_L(i1,i3) * mC(i1) *  C2_Snu_i_C_jk(i2,i1,i3)  )
    End Do
   End Do
  End Do
  Do i1=1,4
   Do i2=1,6
    Ai(1) = Ai(1)                                                           &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * C00_N_i_Sle_jk(i1,i2,i2)
    Ai(2) = Ai(2)                                                           &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * C00_N_i_Sle_jk(i1,i2,i2)
    Ai(3) = Ai(3)                                                           &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i2) )
    Ai(4) = Ai(4)                                                           &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i2) )
    Do i3=i2+1,6
     Ai(1) = Ai(1)                                                          &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * C00_N_i_Sle_jk(i1,i2,i3)             &
        & + 2._dp * cpl_LNSl_R(i, i1, i3) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * C00_N_i_Sle_jk(i1,i3,i2)
     Ai(2) = Ai(2)                                                          &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * C00_N_i_Sle_jk(i1,i2,i3)             &
        & + 2._dp * cpl_LNSl_L(i, i1, i3) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * C00_N_i_Sle_jk(i1,i3,i2)
     Ai(3) = Ai(3)                                                          &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i3) ) &
        & + 2._dp * cpl_LNSl_R(i, i1, i3) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i3,i2) )
     Ai(4) = Ai(4)                                                          &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i3) ) &
        & + 2._dp * cpl_LNSl_L(i, i1, i3) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i3,i2) )
    End Do
    Do i3=1,4
     Ai(1) = Ai(1)                                                           &
      &  - 0.5_dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i3, i2) )   &
      &    * ( cpl_NNZ_R(i1,i3) * ( mSle2(i2) * C0_Sle_i_N_jk(i2,i1,i3)      &
      &                           + B0_Z_NN(i1,i3)                           &
      &                           -  2._dp * C00_Sle_i_N_jk(i2,i1,i3) )      & 
      &      - mN(i1) *mN(i3) *cpl_NNZ_L(i1,i3) * C0_Sle_i_N_jk(i2,i1,i3) )
     Ai(2) = Ai(2)                                                           &
      &  - 0.5_dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i3, i2) )   &
      &    * ( cpl_NNZ_L(i1,i3) * ( mSle2(i2) * C0_Sle_i_N_jk(i2,i1,i3)      &
      &                           + B0_Z_NN(i1,i3)                           &
      &                           -  2._dp * C00_Sle_i_N_jk(i2,i1,i3) )      & 
      &      - mN(i1) *mN(i3) *cpl_NNZ_R(i1,i3) * C0_Sle_i_N_jk(i2,i1,i3) )
     Ai(3) = Ai(3)                                                           &
      &  + cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i3, i2) )            &
      &    * ( cpl_NNZ_L(i1,i3) * mN(i3) *  C1_Sle_i_N_jk(i2,i1,i3)          &
      &      + cpl_NNZ_R(i1,i3) * mN(i1) *  C2_Sle_i_N_jk(i2,i1,i3)  )
     Ai(4) = Ai(4)                                                           &
      &  + cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i3, i2) )            &
      &    * ( cpl_NNZ_R(i1,i3) * mN(i3) *  C1_Sle_i_N_jk(i2,i1,i3)          &
      &      + cpl_NNZ_L(i1,i3) * mN(i1) *  C2_Sle_i_N_jk(i2,i1,i3)  )
    End Do
   End Do
  End Do
  !----------------------------------------------
  ! wave functions
  !----------------------------------------------
  Do i1=1,2
   Do i2=1,3
    Ai(1) = Ai(1) - 2._dp * cpl_LLZ_L * (B0_0_CSn(i1,i2) + B1_0_CSn(i1,i2))    &
          &           * cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_R(i1, j, i2) )
   End Do
  end do
  Do i1=1,4
   Do i2=1,6
    Ai(1) = Ai(1) - 2._dp * cpl_LLZ_L * (B0_0_NSl(i1,i2) + B1_0_NSl(i1,i2))    &
          &           * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i2) )
    Ai(2) = Ai(2) - 2._dp * cpl_LLZ_R * (B0_0_NSl(i1,i2) + B1_0_NSl(i1,i2))    &
          &           * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i2) )
   End Do
  end do

  Ai = oo16pi2 * Ai
  gam = oo48Pi * ( 2._dp * (Abs(Ai(1))**2 + Abs(Ai(2))**2 ) * mZ         &
     &          + 0.25_dp * (Abs(Ai(3))**2 + Abs(Ai(4))**2 ) * mZ * mZ2 ) 
  BR = gam / (gam+gamZ)

  Iname = Iname - 1
 
 End Subroutine ZtoLiLj

  Subroutine Bq_to_ll(i, l, mf_d, mf_l, CKM, c10, c10p, cs, csp, cp, cpp, res )
 !---------------------------------------------------------------------------
 ! Calculates the Branching ratio for the decay of a Bd or Bs meson
 ! into an l+l- pair
 ! Input: i=1,2 ..... d,s
 !        l=1,2,3 ... e,mu,tau
 !        mf_d ...... d-type quark masses
 !        CKM ....... CKM matrix
 !        c10, c10p . C_10, C_10'
 !        cs, csp ... C_S,C_S'
 !        cp, cpp ... C_P,C_P'
 !        
 ! Output: res, branching ratio
 !         
 ! written by Werner Porod, 01 April 2014
 ! the eqn numbers ref. to Bobeth et al., NPB640, 87 (2002)
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_d(3), mf_l(3)
  Complex(dp), Intent(in) :: CKM(3,3), c10, c10p, cs, csp, cp, cpp
  ! reduced Z-squark- mixing couplings
  Integer, Intent(in) :: i, l

  Real(dp), Intent(out) :: res

  Real(dp) :: sW2, mhat
  Complex(dp) :: F_A, F_S, F_P

  sW2 = 1._dp - mW2 / mZ2

  !-----------------------------
  ! axial vector contributions
  !-----------------------------
  F_A = sW2 *(c10 - c10p) / 2._dp
  !------------------------------------------
  ! double penguins
  !------------------------------------------
  mhat = mf_d(i) / mf_d(3)

  F_S = 0.5_dp * MassBq(i)**2 * ( cs - csP*mhat) /(1._dp + mhat)
  F_P = 0.5_dp * MassBq(i)**2 * ( cp - cpP*mhat) /(1._dp + mhat)
  !-------------------------------------------------------------
  ! the formulas of NPB 630 are for \bar{B}_q whereas the ones of
  ! NPB 659 are for B_q 
  !-------------------------------------------------------------
  F_A = Conjg(F_A)
  !---------------------------------------------
  ! the life-time is in ps -> division by hbar
  !---------------------------------------------
  res = (G_F * Alpha * Fb(i) * Abs(CKM(3,3)*Conjg(CKM(3,i))))**2 * MassBq(i) &
    &   / (16._dp*Pi**3 * sW2**2)                                            &
    &    * TauBq(i) * Sqrt(1._dp - 4._dp * mf_l2(l) / MassBq(i)**2)          &
    &    * ( (1._dp - 4._dp * mf_l2(l) / MassBq(i)**2) * Abs(F_S)**2         &
    &      + Abs(F_P + 2._dp * mf_l(l) * F_A)**2 ) / hbar

 End Subroutine Bq_to_ll


 Subroutine K_To_PiNuNu(CKM, c11, c11p, BrK0toPi0NuNu, BrKptoPipNuNu)
 !-------------------------------------------
 ! gives BR(K->pi nu nu) using formulas of
 ! Buras et al., NPB 714 (2005) 103
 ! input: all running quantities are given at Q=M_W
 !  CKM ............. CKM Matrix
 !  C11 ............. Wilson coefficients C_11
 !  C11p ............ Wilson coefficients C_11'
 ! output:
 !  BrK0toPi0NuNu ..................... BR(K^0 -> pi^0 nu nu) 
 !  BrKptoPipNuNuBtoSNuNu ............. BR(K^+ -> pi^+ nu nu) 
 ! written by Werner Porod, 03.02.2014
 ! 16.03.2014: take now Wilson coefficients as input 
 !-------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: CKM(3,3), c11(3), c11p(3)

  Real(dp), Intent(out) :: BrK0toPi0NuNu, BrKptoPipNuNu

  Integer :: i1
  Real(dp) :: Re_Xi(3), Im_Xi(3), lam_wolf, p_c
 
  Iname = Iname + 1
  NameOfUnit(Iname) = "K_To_PiNuNu"

  Do i1=1,3
   Re_Xi(i1) = Real( Conjg(CKM(3,1)) * CKM(3,2) &
               &        * ( C11(i1)+C11p(i1) ), dp )
   Im_Xi(i1) = Aimag( Conjg(CKM(3,1)) * CKM(3,2) &
               &        * ( C11(i1)+C11p(i1) ) )
  End Do

  lam_wolf = Abs(CKM(1,2)/Sqrt(1-Abs(CKM(1,3))**2))

  BrK0ToPi0NuNu = Sum(Im_Xi**2)/lam_wolf**10
  BrKpToPipNuNu =  BrK0ToPi0NuNu

  BrK0ToPi0NuNu = 0.710e-10_dp * BrK0ToPi0NuNu

  ! this includes the NLO c-quark contribution
  p_c = 0.39_dp * Real(Conjg(CKM(2,2)) * CKM(2,1),dp)/lam_wolf
  Do i1=1,3
   BrKpToPipNuNu = BrKpToPipNuNu + ( p_c + Re_Xi(i1)/lam_wolf**5)**2
  End Do

  BrKpToPipNuNu = 1.72e-11_dp * BrKpToPipNuNu

  Iname = Iname - 1

 End Subroutine K_To_PiNuNu


 Subroutine Write_Wilson_Coeff(n, i, kont)
 !---------------------------------------------------
 ! write the Wilson coefficients at the energy Q(i)
 ! to the channel n
 !---------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n, i
  Integer, Intent(inout) :: kont

  Integer :: i2
  Logical :: l_im

  Iname = Iname + 1
  NameOfUnit(Iname) = "Write_Wilson_Coeff"
  !------------------------------------------------
  ! test, if set of coefficients exists
  !------------------------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_i)) Then
   Write(ErrCan,*) "Error in routine Write_Wilson_Coeff"
   Write(ErrCan,*) "set with index",i,"does not exist"
   kont = -1000
   Iname = Iname - 1
   Return
  End If

  ! write only in case that at least one coefficient is non-zero
  If (l_wc(i)) Then
   l_im = .False.

   Write(n,106) "Block FWCOEF Q=",Q_out(i),"# Wilson coefficients at scale Q"
   Write(n,100) "#    id        order  M        value         comment"
   l_im = WriteLineR(n, "0305", "4422", "00", 0, WC_C7(i,3,2,2), "C7", l_im)
   l_im = WriteLineR(n, "0305", "4422", "00", 2, WC_C7(i,3,2,1), "C7", l_im)
   l_im = WriteLineR(n, "0305", "4322", "00", 2, WC_C7p(i,3,2,1), "C7'", l_im)
   l_im = WriteLineR(n, "0305", "6421", "00", 0, WC_C8(i,3,2,2), "C8", l_im)
   l_im = WriteLineR(n, "0305", "6421", "00", 2, WC_C8(i,3,2,1), "C8", l_im)
   l_im = WriteLineR(n, "0305", "6321", "00", 2, WC_C8p(i,3,2,1), "C8'", l_im)

   l_im = WriteLineR(n, "03051111", "4133", "00", 0, WC_C9(i,3,2,1,2) &
        &           , "C9 e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4133", "00", 2, WC_C9(i,3,2,1,1) &
        &           , "C9 e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4233", "00", 0, WC_C9p(i,3,2,1,2) &
        &           , "C9' e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4233", "00", 2, WC_C9p(i,3,2,1,1) &
        &           , "C9' e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4137", "00", 0, WC_C10(i,3,2,1,2) &
        &           , "C10 e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4137", "00", 2, WC_C10(i,3,2,1,1) &
        &           , "C10 e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4237", "00", 0, WC_C10p(i,3,2,1,2) &
        &           , "C10' e+e-", l_im)
   l_im = WriteLineR(n, "03051111", "4237", "00", 2, WC_C10p(i,3,2,1,1) &
        &           , "C10' e+e-", l_im)

   l_im = WriteLineR(n, "03051313", "4133", "00", 0, WC_C9(i,3,2,2,2) &
        &           , "C9 mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4133", "00", 2, WC_C9(i,3,2,2,1) &
        &           , "C9 mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4233", "00", 0, WC_C9p(i,3,2,2,2) &
        &           , "C9' mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4233", "00", 2, WC_C9p(i,3,2,2,1) &
        &           , "C9' mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4137", "00", 0, WC_C10(i,3,2,2,2) &
        &           , "C10 mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4137", "00", 2, WC_C10(i,3,2,2,1) &
        &           , "C10 mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4237", "00", 0, WC_C10p(i,3,2,2,2) &
        &           , "C10' mu+mu-", l_im)
   l_im = WriteLineR(n, "03051313", "4237", "00", 2, WC_C10p(i,3,2,2,1) &
        &           , "C10' mu+mu-", l_im)

   Do i2=1,3
    l_im = WriteLineR(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4137", "00" &
         & , 0, WC_C11(i,3,2,i2,2), "C11 nu_"//Bu(i2)//" nu_"//Bu(i2), l_im)
    l_im = WriteLineR(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4137", "00" &
         & , 2, WC_C11(i,3,2,i2,1), "C11 nu_"//Bu(i2)//" nu_"//Bu(i2), l_im)
    l_im = WriteLineR(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4237", "00" &
         & , 0, WC_C11p(i,3,2,i2,2), "C11' nu_"//Bu(i2)//" nu_"//Bu(i2), l_im)
    l_im = WriteLineR(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4237", "00" &
         & , 2, WC_C11p(i,3,2,i2,1), "C11' nu_"//Bu(i2)//" nu_"//Bu(i2), l_im)
   End Do

   If (l_im) Then
    Write(n,106) "Block IMFWCOEF Q=",Q_out(i) &
              & ,"# Im(Wilson coefficients) at scale Q"
    Write(n,100) "#    id        order  M        value         comment"

    l_im = WriteLineI(n, "0305", "4422", "00", 0, WC_C7(i,3,2,2), "C7")
    l_im = WriteLineI(n, "0305", "4422", "00", 2, WC_C7(i,3,2,1), "C7")
    l_im = WriteLineI(n, "0305", "4322", "00", 2, WC_C7p(i,3,2,1), "C7'")
    l_im = WriteLineI(n, "0305", "6421", "00", 0, WC_C8(i,3,2,2), "C8")
    l_im = WriteLineI(n, "0305", "6421", "00", 2, WC_C8(i,3,2,1), "C8")
    l_im = WriteLineI(n, "0305", "6321", "00", 2, WC_C8p(i,3,2,1), "C8'")

    l_im = WriteLineI(n, "03051111", "4133", "00", 0, WC_C9(i,3,2,1,2) &
        &           , "C9 e+e-")
    l_im = WriteLineI(n, "03051111", "4133", "00", 2, WC_C9(i,3,2,1,1) &
        &           , "C9 e+e-")
    l_im = WriteLineI(n, "03051111", "4233", "00", 0, WC_C9p(i,3,2,1,2) &
        &           , "C9' e+e-")
    l_im = WriteLineI(n, "03051111", "4233", "00", 2, WC_C9p(i,3,2,1,1) &
        &           , "C9' e+e-")
    l_im = WriteLineI(n, "03051111", "4137", "00", 0, WC_C10(i,3,2,1,2) &
        &           , "C10 e+e-")
    l_im = WriteLineI(n, "03051111", "4137", "00", 2, WC_C10(i,3,2,1,1) &
        &           , "C10 e+e-")
    l_im = WriteLineI(n, "03051111", "4237", "00", 0, WC_C10p(i,3,2,1,2) &
        &           , "C10' e+e-")
    l_im = WriteLineI(n, "03051111", "4237", "00", 2, WC_C10p(i,3,2,1,1) &
        &           , "C10' e+e-")

    l_im = WriteLineI(n, "03051313", "4133", "00", 0, WC_C9(i,3,2,2,2) &
        &           , "C9 mu+mu-")
    l_im = WriteLineI(n, "03051313", "4133", "00", 2, WC_C9(i,3,2,2,1) &
        &           , "C9 mu+mu-")
    l_im = WriteLineI(n, "03051313", "4233", "00", 0, WC_C9p(i,3,2,2,2) &
        &           , "C9' mu+mu-")
    l_im = WriteLineI(n, "03051313", "4233", "00", 2, WC_C9p(i,3,2,2,1) &
        &           , "C9' mu+mu-")
    l_im = WriteLineI(n, "03051313", "4137", "00", 0, WC_C10(i,3,2,2,2) &
        &           , "C10 mu+mu-")
    l_im = WriteLineI(n, "03051313", "4137", "00", 2, WC_C10(i,3,2,2,1) &
        &           , "C10 mu+mu-")
    l_im = WriteLineI(n, "03051313", "4237", "00", 0, WC_C10p(i,3,2,2,2) &
        &           , "C10' mu+mu-")
    l_im = WriteLineI(n, "03051313", "4237", "00", 2, WC_C10p(i,3,2,2,1) &
        &           , "C10' mu+mu-")

    Do i2=1,3
     l_im = WriteLineI(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4137", "00" &
         & , 0, WC_C11(i,3,2,i2,2), "C11 nu_"//Bu(i2)//" nu_"//Bu(i2))
     l_im = WriteLineI(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4137", "00" &
          & , 2, WC_C11(i,3,2,i2,1), "C11 nu_"//Bu(i2)//" nu_"//Bu(i2))
    l_im = WriteLineI(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4237", "00" &
         & , 0, WC_C11p(i,3,2,i2,2), "C11' nu_"//Bu(i2)//" nu_"//Bu(i2))
     l_im = WriteLineI(n, "03051"//Bu(2*i2)//"1"//Bu(2*i2), "4237", "00" &
         & , 2, WC_C11p(i,3,2,i2,1), "C11' nu_"//Bu(i2)//" nu_"//Bu(i2))
    End Do
   End If
  End If
  Iname = Iname - 1

  100 Format(a)
  106 Format(a,1P,e16.8,2x,a)

 Contains 

  Logical Function WriteLineR(n, id, id2, ord, i_in,val, com, l_im)
   Implicit None
   Integer, Intent(in) :: n, i_in
   Character(len=*), Intent(in) :: id, id2, ord, com
   Complex(dp), Intent(in) :: val
   Logical, Intent(in) :: l_im

   If (Real(val,dp).Ne.0._dp) &
    & Write(n,100) Trim(id),Trim(id2),Trim(ord),i_in,Real(val,dp),Trim(com)
   WriteLineR = .False.
   If ((Aimag(val).Ne.0._dp).Or.l_im) WriteLineR = .True.

   100 Format(1x,a8,1x,a4,3x,a2,3x,i1,3x,1P,e16.8,0P,3x,'#',1x,a)

  End Function WriteLineR

  Logical Function WriteLineI(n, id, id2, ord, i_in,val, com)
   Implicit None
   Integer, Intent(in) :: n, i_in
   Character(len=*), Intent(in) :: id, id2,ord, com
   Complex(dp), Intent(in) :: val

   If (Aimag(val).Ne.0._dp) &
     & Write(n,100) Trim(id),Trim(id2),Trim(ord),i_in,Aimag(val),Trim(com)
   WriteLineI = .False.

   100 Format(1x,a8,1x,a4,3x,a2,3x,i1,3x,1P,e16.8,0P,3x,'#',1x,a)

  End Function WriteLineI

 End Subroutine Write_Wilson_Coeff

End Module LowEnergy
