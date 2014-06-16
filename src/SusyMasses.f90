Module SusyMasses

! load modules
Use Control
Use Mathematics
Use StandardModel, Only: mW, mW2, mZ, mZ2, mf_l, mf_l2, mf_u, mf_u2, mf_d &
    & , mf_d2, AlphaS_mB, alphaS_mZ, G_F, FermionMass, QuarkMasses_and_PhaseShifts
! load modules


! interfaces
 Interface ChargedScalarMass
  Module Procedure ChargedScalarMassEps1nt, ChargedScalarMassEps1  &
               & , ChargedScalarMassEps3nt, ChargedScalarMassEps3  &
               & , ChargedScalarMassLam3nt, ChargedScalarMassLam3
 End Interface

 Interface CharginoMass
  Module Procedure CharginoMass2, CharginoMass3, CharginoMass3a &
               & , CharginoMass5, CharginoMass5Lam
 End Interface

 Interface NeutralinoMass
  Module Procedure NeutralinoMass4, NeutralinoMass5, NeutralinoMass5NMSSM &
               & , NeutralinoMass7
 End Interface

 Interface PseudoScalarMass
  Module Procedure PseudoScalarMassMSSMnT, PseudoScalarMassEps1nT, &
                 & PseudoScalarMassEps3nT,  PseudoScalarMassMSSM,  &
                 & PseudoScalarMassEps1, PseudoScalarMassEps3
 End Interface

 Interface ScalarMass
  Module Procedure ScalarMassMSSMnT, ScalarMassEps1nT, &
                 & ScalarMassEps3nT,  ScalarMassMSSM,  &
                 & ScalarMassEps1, ScalarMassEps3,     &
                 & ScalarMassMSSMeff, ScalarMassNMSSMeff
 End Interface

 Interface SfermionMass
  Module Procedure SfermionMass1mssm, SfermionMass1eps1, SfermionMass1eps3, &
                 & SfermionMass3MSSM 
 End Interface

 Interface TreeMasses
  Module Procedure TreeMassesMSSM, TreeMassesEps3, TreeMassesLam3           &
               &,  TreeMassesMSSM2, TreeMassesNMSSM
 End Interface
! interfaces

Contains



 Subroutine ChargedScalarMassEps1(bi, B, vevSM, vevL, ME2, A &
                                &, Yuk, gp, g, mSpm, mSpm2, RSpm, kont)
!                                &, Yuk, gp, g, tad_i, mSpm, mSpm2, RSpm, kont)
 !-----------------------------------------------------------------------
 ! calculates the charged scalar masses in the 3-generation epsilon model
 ! Input:
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  ML2(i,j) ..... M_Lij squared
 !  ME2(i,j) ..... M_Eij squared
 !  A(i,j) ....... A_e(i,j)
 !  Yuk(i,j) ..... Y_e(i,j)
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mSpm(i) ...... the masses
 !  mSpm2(i) ..... the masses squared
 !  RSpm(i,J) .... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 21.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: ME2, vevSM(2), vevL, g, gp!, tad_i(5)
  Complex(Dp), Intent(in) :: A, bi(2), B(2), yuk
  Real(Dp), Intent(out) :: mSpm(4), mSpm2(4)
  Complex(Dp), Intent(out) :: RSpm(4,4)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm, mL2a
  Complex(Dp) :: sumC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassEps1'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2 + vevL**2 )

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * bi(2) * vevL,dp) ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1) - Real( B(2),dp ) * vevL ) / vevSM(2)
   

  sumC = - vevSM(2) * B(2) + vevSM(1) * Conjg(bi(1)) * bi(2)  &
       & - Conjg( bi(2) ) * bi(2) * vevL
  ML2a = - DTerm + Real( sumC,dp ) / vevL

  Call ChargedScalarMassEps1nt(MH1sq, MH2sq, bi, B,vevSM, vevL, mL2a, mE2,A &
  &                           ,yuk, gp, g, mSpm, mSpm2, RSpm, kont)

  Iname = Iname - 1

 End Subroutine ChargedScalarMassEps1

 Subroutine ChargedScalarMassEps1nt(MH1sq,MH2sq,bi,Bmu,vevSM,vevL,mL2,mR2,A &
  &                                ,yuk,gp,g,mSpm,mSpm2,RSpm,kont           &
  &                                ,NoSymmetryBreaking)
 !------------------------------------------------------------------------
 ! calculates charged scalar masses + mixing matrix R in the 3-generation
 ! epsilon model without refering to the tadpoles.
 ! input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  bi(i) ...... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  Bmu(i) ..... i=1 mu*B
 !              i=2 epsilon_1 * B_1
 !              i=3 epsilon_2 * B_2
 !              i=4 epsilon_3 * B_3
 !  vevSM(i) ... i=1 v_d
 !               i=2 v_u
 !  vevL(i) .... i=1 v_1
 !               i=2 v_2
 !               i=3 v_3
 !  mL(i,j) ... left slepton mass matrix, i,j=1,3
 !  mR(i,j) ... right slepton mass matrix, i,j=1,3
 !  A(i,j) .... trilinear slepton Higgs parameters, i,j=1,3
 !  yuk(i) .... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! output:
 !  mSpm(i) .... charged scalar masses, mSpm(1)=mW, xi=1 gauge
 !  mSpm2(i) ... = mSpm(i)**2
 !  RSpm(i,j) .. mixing matrix of charged scalar
 ! written by Werner Porod, 18.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2), vevL, g, gp, MH1sq, MH2sq, mL2, mR2
  Real(Dp), Intent(out) :: mSpm(4), mSpm2(4)
  Complex(Dp), Intent(in) :: bi(2), Bmu(2), A, yuk
  Complex(Dp), Intent(out) :: RSpm(4,4)
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Integer :: i1,i2,ierr!, nrot, i, j, ii, jj, imais, imenos
  Complex(Dp) :: mat4(4,4),vec, yukA , yuk2p!, matOut(4,4) ,RT(4,4)
  Real(Dp) :: gZ,v12,v22,v32, DTerm,g2,gp2, test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassEps1nt'

  gp2 = gp**2
  g2 = g**2
  gZ = 0.25_dp * g2
  YukA = Conjg( Yuk )
  yuk2p = Abs( Yuk )**2

  mSpm = 0._dp
  mSpm2 = 0._dp
  RSpm = (0._dp,0._dp)

  mat4 = (0._dp,0._dp)

  v12 = vevSM(1)**2
  v22 = vevSM(2)**2
  v32 = vevL**2

  DTerm = 0.125_dp * ( (g2+gp2) * v12 + (g2-gp2) * ( v22 - v32) )

  mat4(1,1) = MH1sq + Abs( bi(1) )**2 + DTerm  + 0.5_dp * yuk2p * v32
  mat4(1,2) = Bmu(1) + gZ * vevSM(1) * vevSM(2)
  vec = 0.5_dp * vevSM(1) * vevL * yuk2p
  mat4(1,3) = gZ * vevSM(1) * vevL - Conjg( bi(1) ) * bi(2) - vec
  mat4(1,4) = - (vevSM(2) * bi(2) * yukA  + vevL*Conjg(A) ) * oosqrt2

  DTerm = 0.125_dp * ( (g2+gp2) * v22 + (g2-gp2) * ( v12 + v32) )
  mat4(2,2) = MH2sq + Dot_product(Conjg(bi), bi) + DTerm
  mat4(2,3) = - Bmu(2) + gZ * vevSM(2) * vevL
  mat4(2,4) = - (vevSM(1) * bi(2) * yukA + Conjg(bi(1)) * vevL * yukA) *oosqrt2

  DTerm = 0.125_dp * (g2-gp2) * ( v12 - v22 + v32 )
  mat4(3,3) = mL2 + 0.5_dp * v12 * yuk2p - DTerm + Abs(bi(2))**2 + gZ * v32

  mat4(3,4) = ( - bi(1) * vevSM(2) * yukA + Conjg(A) * vevSM(1) ) * oosqrt2

  DTerm = 0.25_dp * gp**2 * (v22 - v12 - v32)
  mat4(4,4) = mR2 + 0.5_dp * (v12 + v32) * yuk2p + Dterm 

  Do i1=2,4
   Do i2 = 1,i1-1
    mat4(i1,i2) = Conjg( mat4(i2,i1) )
   End Do
  End Do

  Call EigenSystem(mat4,mSpm2,RSpm,ierr, test)
  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'In subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'subroutine Eigensystem gives ierr = ',ierr
   Write(ErrCan,*) 'MH1sq ',MH1sq
   Write(ErrCan,*) 'MH2sq ',MH2sq
   Write(ErrCan,*) 'bi ',bi
   Write(ErrCan,*) 'Bmu ',Bmu
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) 'g ',g
   Write(ErrCan,*) 'yuk ',yuk
   Write(ErrCan,*) 'ML2 ',ML2
   Write(ErrCan,*) 'MR2 ',MR2
   Write(ErrCan,*) 'A ',A
   Write(ErrCan,*) '  '
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,4
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
      Write(ErrCan,*) 'MH1sq ',MH1sq
      Write(ErrCan,*) 'MH2sq ',MH2sq
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'Bmu ',Bmu
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'g ',g
      Write(ErrCan,*) 'yuk ',yuk
      Write(ErrCan,*) 'ML2 ',ML2
      Write(ErrCan,*) 'MR2 ',MR2
      Write(ErrCan,*) 'A ',A
      Write(ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -201
    Call AddError(201)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mSpm(1) = mW
   mSpm2(1) = mW**2
   Do i1=2,4
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
      Write(ErrCan,*) 'MH1sq ',MH1sq
      Write(ErrCan,*) 'MH2sq ',MH2sq
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'Bmu ',Bmu
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'g ',g
      Write(ErrCan,*) 'yuk ',yuk
      Write(ErrCan,*) 'ML2 ',ML2
      Write(ErrCan,*) 'MR2 ',MR2
      Write(ErrCan,*) 'A ',A
      Write(ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -201
    Call AddError(201)
    End If
   End Do
  End If

  Iname = Iname - 1

 End Subroutine ChargedScalarMassEps1nt

 Subroutine ChargedScalarMassEps3(bi, B, vevSM, vevL, ML2, ME2, A &
                                &, Yuk, gp, g, mSpm, mSpm2, RSpm, kont)
!                                &, Yuk, gp, g, tad_i, mSpm, mSpm2, RSpm, kont)
 !-----------------------------------------------------------------------
 ! calculates the charged scalar masses in the 3-generation epsilon model
 ! Input:
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  ML2(i,j) ..... M_Lij squared
 !  ME2(i,j) ..... M_Eij squared
 !  A(i,j) ....... A_e(i,j)
 !  Yuk(i,j) ..... Y_e(i,j)
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mSpm(i) ...... the masses
 !  mSpm2(i) ..... the masses squared
 !  RSpm(i,J) .... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 21.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp!, tad_i(5)
  Complex(Dp), Intent(in) :: ML2(3,3), ME2(3,3), A(3,3), bi(4), B(4) &
                               &, yuk(3,3)
  Real(Dp), Intent(out) :: mSpm(8), mSpm2(8)
  Complex(Dp), Intent(out) :: RSpm(8,8)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm
  Complex(Dp) :: ML2a(3,3), sumC
  Integer :: i1,i2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassEps3'

  DTerm = 0.125_dp * (g**2+gp**2)  &
      & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * Dot_product(bi(2:4),vevL),dp ) ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1)              &
      &   - Dot_product( Real(B(2:4),dp) , vevL) ) / vevSM(2)
   
  ML2a = ML2
  Do i1=1,3
   sumC = - vevSM(2) * B(i1+1) + vevSM(1) * Conjg(bi(1)) * bi(i1+1)
   Do i2=1,3
    sumC = sumC - Conjg( bi(i1+1) ) * bi(i2+1) * vevL(i2)
    If (i1.Ne.i2) sumC = sumC - vevL(i2) * ML2a(i2,i1)
   End Do
   ML2a(i1,i1) = - DTerm + Real( sumC,dp ) / vevL(i1) !- tad_i(i1+2) 
  End Do

!  Mh1sq = mh1sq + tad_i(1)
!  Mh2sq = mh1sq + tad_i(2)
  Call ChargedScalarMassEps3nt(MH1sq, MH2sq, bi, B,vevSM, vevL, mL2a, mE2,A &
  &                           ,yuk, gp, g, mSpm, mSpm2, RSpm, kont)

  Iname = Iname - 1

 End Subroutine ChargedScalarMassEps3

 Subroutine ChargedScalarMassEps3nt(MH1sq,MH2sq,bi,Bmu,vevSM,vevL,mL2,mR2,A &
  &                                ,yuk,gp,g,mSpm,mSpm2,RSpm,kont           &
  &                                ,NoSymmetryBreaking)
 !------------------------------------------------------------------------
 ! calculates charged scalar masses + mixing matrix R in the 3-generation
 ! epsilon model without refering to the tadpoles.
 ! input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  bi(i) ...... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  Bmu(i) ..... i=1 mu*B
 !              i=2 epsilon_1 * B_1
 !              i=3 epsilon_2 * B_2
 !              i=4 epsilon_3 * B_3
 !  vevSM(i) ... i=1 v_d
 !               i=2 v_u
 !  vevL(i) .... i=1 v_1
 !               i=2 v_2
 !               i=3 v_3
 !  mL(i,j) ... left slepton mass matrix, i,j=1,3
 !  mR(i,j) ... right slepton mass matrix, i,j=1,3
 !  A(i,j) .... trilinear slepton Higgs parameters, i,j=1,3
 !  yuk(i) .... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! output:
 !  mSpm(i) .... charged scalar masses, mSpm(1)=mW, xi=1 gauge
 !  mSpm2(i) ... = mSpm(i)**2
 !  RSpm(i,j) .. mixing matrix of charged scalar
 ! written by Werner Porod, 18.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp, MH1sq, MH2sq
  Real(Dp), Intent(out) :: mSpm(8), mSpm2(8)
  Complex(Dp), Intent(in) :: bi(4), Bmu(4), A(3,3), yuk(3,3), mL2(3,3) &
                               &, mR2(3,3)
  Complex(Dp), Intent(out) :: RSpm(8,8)
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Integer :: i1,i2,ierr!, nrot, i, j, ii, jj, imais, imenos
  Complex(Dp) :: mat8(8,8),vec(3),yukA(3,3),yuk2p(3,3)!,matOut(8,8),RT(8,8)
  Real(Dp) :: gZ,v12,v22,v32,v42,v52,DTerm,g2,gp2, test(2)
!  Real(dp) :: mat8R(8,8), d1(8), v1(8,8),v2(8,8),mass8(8,8), work(8)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassEps3nt'

  gp2 = gp**2
  g2 = g**2
  gZ = 0.25_dp * g2
  Call Adjungate(yuk,yukA)
  yuk2p = Matmul( Conjg(yuk), Transpose(yuk) )

  mSpm = 0._dp
  mSpm2 = 0._dp
  RSpm = (0._dp,0._dp)

  mat8 = (0._dp,0._dp)

  v12 = vevSM(1)**2
  v22 = vevSM(2)**2
  v32 = vevL(1)**2
  v42 = vevL(2)**2
  v52 = vevL(3)**2

  DTerm = 0.125_dp * ( (g2+gp2) * v12                      &
      &              + (g2-gp2) * ( v22 - v32 - v42 - v52) )

  mat8(1,1) = MH1sq + Abs( bi(1) )**2 + DTerm &
          & + 0.5_dp * Dot_product(Matmul(Matmul(vevL,yukA),yuk),vevL)
  mat8(1,2) = Bmu(1) + gZ * vevSM(1) * vevSM(2)
  vec = 0.5_dp * vevSM(1) * Matmul(vevL,yuk2p)
  mat8(1,3:5) = gZ * vevSM(1) * vevL - Conjg( bi(1) ) * bi(2:4) - vec
  mat8(1,6:8) = vevSM(2) * Conjg( Matmul(bi(2:4), yuk ) )  &
            & + Matmul(vevL,Conjg(A))
  mat8(1,6:8) = - mat8(1,6:8) * oosqrt2

  DTerm = 0.125_dp * ( (g2+gp2) * v22                      &
      &              + (g2-gp2) * ( v12 + v32 + v42 + v52) )
  mat8(2,2) = MH2sq + Dot_product(Conjg(bi), bi) + DTerm
  mat8(2,3:5) = - Bmu(2:4) + gZ * vevSM(2) * vevL
  mat8(2,6:8) = vevSM(1) * Matmul(bi(2:4),Conjg(yuk)) &
            & + Conjg(bi(1)) * Matmul(vevL,Conjg(yuk))
  mat8(2,6:8) = - mat8(2,6:8) * oosqrt2

  DTerm = 0.125_dp * (g2-gp2) * ( v12 - v22 + v32 + v42 + v52)
  mat8(3:5,3:5) = mL2 + 0.5_dp * v12 * yuk2p
  Do i1=1,3
   mat8(i1+2,i1+2) =  mat8(i1+2,i1+2) - DTerm
   Do i2=i1,3
    mat8(i1+2,i2+2) = mat8(i1+2,i2+2) + Conjg( bi(i1+1) ) * bi(i2+1) &
                    &                 + gZ * vevL(i1) * vevL(i2)
   End Do
  End Do

  mat8(3:5,6:8) = ( - bi(1) * vevSM(2) * Conjg(yuk) &
                &   + Conjg(A) * vevSM(1) ) * oosqrt2

  DTerm = 0.25_dp * gp**2 * (v22 - v12 - v32 - v42 - v52)
  mat8(6:8,6:8) = mR2 + 0.5_dp * v12 * Matmul( Transpose(yuk), Conjg(yuk) )
  vec = Matmul(vevL,yuk)
  Do i1=1,3
   mat8(i1+5,i1+5) = mat8(i1+5,i1+5) + Dterm
   Do i2=i1,3
    mat8(i1+5,i2+5) = mat8(i1+5,i2+5) + 0.5_dp * vec(i1) * Conjg(vec(i2))
   End Do
  End Do
  Do i1=2,8
   Do i2 = 1,i1-1
    mat8(i1,i2) = Conjg( mat8(i2,i1) )
   End Do
  End Do

  Call EigenSystem(mat8,mSpm2,RSpm,ierr, test)
  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (ErrCan,*) 'In subroutine ',NameOfUnit(Iname)
   Write (ErrCan,*) 'subroutine EigenSystem gives ierr = ',ierr
   Write (ErrCan,*) 'MH1sq ',MH1sq
   Write (ErrCan,*) 'MH2sq ',MH2sq
   Write (ErrCan,*) 'bi ',bi
   Write (ErrCan,*) 'Bmu ',Bmu
   Write (ErrCan,*) 'vevSM ',vevSM
   Write (ErrCan,*) 'vevL ',vevL
   Write (ErrCan,*) 'g ',g
   Write (ErrCan,*) 'yuk ',yuk
   Write (ErrCan,*) 'ML2(1,i1) ',(ML2(1,i1),i1=1,3)
   Write (ErrCan,*) 'ML2(2,i1) ',(ML2(2,i1),i1=1,3)
   Write (ErrCan,*) 'ML2(3,i1) ',(ML2(3,i1),i1=1,3)
   Write (ErrCan,*) 'MR2(1,i1) ',(MR2(1,i1),i1=1,3)
   Write (ErrCan,*) 'MR2(2,i1) ',(MR2(2,i1),i1=1,3)
   Write (ErrCan,*) 'MR2(3,i1) ',(MR2(3,i1),i1=1,3)
   Write (ErrCan,*) 'A(1,i1) ',(A(1,i1),i1=1,3)
   Write (ErrCan,*) 'A(2,i1) ',(A(2,i1),i1=1,3)
   Write (ErrCan,*) 'A(3,i1) ',(A(3,i1),i1=1,3)
   Write (ErrCan,*) '  '
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,8
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write (ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write (ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
      Write (ErrCan,*) 'MH1sq ',MH1sq
      Write (ErrCan,*) 'MH2sq ',MH2sq
      Write (ErrCan,*) 'bi ',bi
      Write (ErrCan,*) 'Bmu ',Bmu
      Write (ErrCan,*) 'vevSM ',vevSM
      Write (ErrCan,*) 'vevL ',vevL
      Write (ErrCan,*) 'g ',g
      Write (ErrCan,*) 'yuk ',yuk
      Write (ErrCan,*) 'ML2(1,i2) ',(ML2(1,i2),i2=1,3)
      Write (ErrCan,*) 'ML2(2,i2) ',(ML2(2,i2),i2=1,3)
      Write (ErrCan,*) 'ML2(3,i2) ',(ML2(3,i2),i2=1,3)
      Write (ErrCan,*) 'MR2(1,i2) ',(MR2(1,i2),i2=1,3)
      Write (ErrCan,*) 'MR2(2,i2) ',(MR2(2,i2),i2=1,3)
      Write (ErrCan,*) 'MR2(3,i2) ',(MR2(3,i2),i2=1,3)
      Write (ErrCan,*) 'A(1,i2) ',(A(1,i2),i2=1,3)
      Write (ErrCan,*) 'A(2,i2) ',(A(2,i2),i2=1,3)
      Write (ErrCan,*) 'A(3,i2) ',(A(3,i2),i2=1,3)
      Write (ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     mSpm2(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -202
    Call AddError(202)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mSpm(1) = mW
   mSpm2(1) = mW**2
   Do i1=2,8
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write (ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write (ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
!       Write (ErrCan,*) 'MH1sq ',MH1sq
!       Write (ErrCan,*) 'MH2sq ',MH2sq
!       Write (ErrCan,*) 'bi ',bi
!       Write (ErrCan,*) 'Bmu ',Bmu
!       Write (ErrCan,*) 'vevSM ',vevSM
!       Write (ErrCan,*) 'vevL ',vevL
!       Write (ErrCan,*) 'g ',g
!       Write (ErrCan,*) 'yuk ',yuk
!       Write (ErrCan,*) 'ML2(1,i2) ',(ML2(1,i2),i2=1,3)
!       Write (ErrCan,*) 'ML2(2,i2) ',(ML2(2,i2),i2=1,3)
!       Write (ErrCan,*) 'ML2(3,i2) ',(ML2(3,i2),i2=1,3)
!       Write (ErrCan,*) 'MR2(1,i2) ',(MR2(1,i2),i2=1,3)
!       Write (ErrCan,*) 'MR2(2,i2) ',(MR2(2,i2),i2=1,3)
!       Write (ErrCan,*) 'MR2(3,i2) ',(MR2(3,i2),i2=1,3)
!       Write (ErrCan,*) 'A(1,i2) ',(A(1,i2),i2=1,3)
!       Write (ErrCan,*) 'A(2,i2) ',(A(2,i2),i2=1,3)
!       Write (ErrCan,*) 'A(3,i2) ',(A(3,i2),i2=1,3)
     If (Maxval(Abs(Aimag(mat8))).Eq.0._dp) Then
       Do i2=1,8
        Write(ErrCan,876) Real(mat8(i2,:),dp)
       End Do
     Else
       Do i2=1,8
        Write(ErrCan,876) mat8(i2,:)
       End Do
     End If 

876 Format(1p,8e12.3)
      Write (ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     mSpm2(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -202
    Call AddError(202)
    End If
   End Do
  End If

  Iname = Iname - 1

 End Subroutine ChargedScalarMassEps3nt


 Subroutine ChargedScalarMassLam3(bi, B, vevSM, vevL, ML2, ME2, A &
                           &, Yuk, lam, Alam, gp, g, mSpm, mSpm2, RSpm, kont)
 !-----------------------------------------------------------------------
 ! calculates the charged scalar masses in the 3-generation epsilon model
 ! including trilinear couplings
 ! Input:
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  ML2(i,j) ..... M_Lij squared
 !  ME2(i,j) ..... M_Eij squared
 !  A(i,j) ....... A_e(i,j)
 !  Yuk(i,j) ..... Y_e(i,j)
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mSpm(i) ...... the masses
 !  mSpm2(i) ..... the masses squared
 !  RSpm(i,J) .... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 25.08.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp
  Complex(Dp), Intent(in) :: ML2(3,3), ME2(3,3), A(3,3), bi(4), B(4) &
                               &, yuk(3,3), lam(3,3,3), Alam(3,3,3)
  Real(Dp), Intent(out) :: mSpm(8), mSpm2(8)
  Complex(Dp), Intent(out) :: RSpm(8,8)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm
  Complex(Dp) :: ML2a(3,3), sumC
  Integer :: i1,i2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassLam3'

  DTerm = 0.125_dp * (g**2+gp**2)  &
      & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * Dot_product(bi(2:4),vevL),dp ) ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1)              &
      &   - Dot_product( Real(B(2:4),dp) , vevL) ) / vevSM(2)
   
  ML2a = ML2
  Do i1=1,3
   sumC = - vevSM(2) * B(i1+1) + vevSM(1) * Conjg(bi(1)) * bi(i1+1)
   Do i2=1,3
    sumC = sumC - Conjg( bi(i1+1) ) * bi(i2+1) * vevL(i2)
    If (i1.Ne.i2) sumC = sumC - vevL(i2) * ML2a(i2,i1)
   End Do
   ML2a(i1,i1) = - DTerm + Real( sumC,dp ) / vevL(i1) 
  End Do

  Call ChargedScalarMassLam3nt(MH1sq, MH2sq, bi, B,vevSM, vevL, mL2a, mE2,A &
  &                           ,yuk, lam, Alam, gp, g, mSpm, mSpm2, RSpm, kont)

  Iname = Iname - 1

 End Subroutine ChargedScalarMassLam3

 Subroutine ChargedScalarMassLam3nt(MH1sq,MH2sq,bi,Bmu,vevSM,vevL,mL2,mR2,A   &
  &                                ,yuk, lam, Alam, gp,g,mSpm,mSpm2,RSpm,kont &
  &                                ,NoSymmetryBreaking)
 !------------------------------------------------------------------------
 ! calculates charged scalar masses + mixing matrix R in the 3-generation
 ! epsilon model without refering to the tadpoles.
 ! input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  bi(i) ...... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  Bmu(i) ..... i=1 mu*B
 !              i=2 epsilon_1 * B_1
 !              i=3 epsilon_2 * B_2
 !              i=4 epsilon_3 * B_3
 !  vevSM(i) ... i=1 v_d
 !               i=2 v_u
 !  vevL(i) .... i=1 v_1
 !               i=2 v_2
 !               i=3 v_3
 !  mL(i,j) ... left slepton mass matrix, i,j=1,3
 !  mR(i,j) ... right slepton mass matrix, i,j=1,3
 !  A(i,j) .... trilinear slepton Higgs parameters, i,j=1,3
 !  yuk(i) .... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! output:
 !  mSpm(i) .... charged scalar masses, mSpm(1)=mW, xi=1 gauge
 !  mSpm2(i) ... = mSpm(i)**2
 !  RSpm(i,j) .. mixing matrix of charged scalar
 ! written by Werner Porod, 18.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp, MH1sq, MH2sq
  Real(Dp), Intent(out) :: mSpm(8), mSpm2(8)
  Complex(Dp), Intent(in) :: bi(4), Bmu(4), A(3,3), yuk(3,3), mL2(3,3) &
                               &, mR2(3,3), lam(3,3,3), Alam(3,3,3)
  Complex(Dp), Intent(out) :: RSpm(8,8)
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Integer :: i1, i2, i3, ierr
  Complex(Dp) :: mat8(8,8), vec(3), yukA(3,3), yukV(3,3), yukVA(3,3)
  Real(Dp) :: gZ,v12,v22,v32,v42,v52,DTerm,g2,gp2, test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedScalarMassLam3nt'

  gp2 = gp**2
  g2 = g**2
  gZ = 0.25_dp * g2
  Call Adjungate(yuk,yukA)
  yukV = vevSM(1) * yuk
  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     yukV(i1,i2) = yukV(i1,i2) + vevL(i3) * ( lam(i3,i1,i2) - lam(i1,i3,i2) )
    End Do
   End Do
  End Do
  Call Adjungate(yukV,yukVA)

  mSpm = 0._dp
  mSpm2 = 0._dp
  RSpm = (0._dp,0._dp)

  mat8 = (0._dp,0._dp)

  v12 = vevSM(1)**2
  v22 = vevSM(2)**2
  v32 = vevL(1)**2
  v42 = vevL(2)**2
  v52 = vevL(3)**2

  DTerm = 0.125_dp * ( (g2+gp2) * v12                      &
      &              + (g2-gp2) * ( v22 - v32 - v42 - v52) )

  mat8(1,1) = MH1sq + Abs( bi(1) )**2 + DTerm &
          & + 0.5_dp * Dot_product(Matmul(Matmul(vevL,yukA),yuk),vevL)
  mat8(1,2) = Bmu(1) + gZ * vevSM(1) * vevSM(2)
  vec = 0.5_dp * Matmul(vevL,Matmul( Conjg(yuk), Transpose(yukV) ))
  mat8(1,3:5) = gZ * vevSM(1) * vevL - Conjg( bi(1) ) * bi(2:4) - vec
  mat8(1,6:8) = vevSM(2) * Conjg( Matmul(bi(2:4), yuk ) )  &
            & + Matmul(vevL,Conjg(A))
  mat8(1,6:8) = - mat8(1,6:8) * oosqrt2

  DTerm = 0.125_dp * ( (g2+gp2) * v22                      &
      &              + (g2-gp2) * ( v12 + v32 + v42 + v52) )
  mat8(2,2) = MH2sq + Dot_product(Conjg(bi), bi) + DTerm
  mat8(2,3:5) = - Bmu(2:4) + gZ * vevSM(2) * vevL
  mat8(2,6:8) = Matmul(bi(2:4),Conjg(yukV)) &
            & + Conjg(bi(1)) * Matmul(vevL,Conjg(yuk))
  mat8(2,6:8) = - mat8(2,6:8) * oosqrt2

  DTerm = 0.125_dp * (g2-gp2) * ( v12 - v22 + v32 + v42 + v52)
  mat8(3:5,3:5) = mL2 + 0.5_dp * Matmul(yukVA,yukV)
  Do i1=1,3
   mat8(i1+2,i1+2) =  mat8(i1+2,i1+2) - DTerm
   Do i2=i1,3
    mat8(i1+2,i2+2) = mat8(i1+2,i2+2) + Conjg( bi(i1+1) ) * bi(i2+1) &
                    &                 + gZ * vevL(i1) * vevL(i2)
   End Do
  End Do

  mat8(3:5,6:8) = ( - bi(1) * vevSM(2) * Conjg(yuk) &
                &   + Conjg(A) * vevSM(1) ) * oosqrt2
  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     mat8(i1+2,i2+5) = mat8(i1+2,i2+5)                                     &
       &  + ( Conjg(Alam(i3,i1,i2) - Alam(i1,i3,i2))*vevL(i3)              &
       &    + vevSM(2) * bi(1+i3) * Conjg(lam(i3,i1,i2) - lam(i1,i3,i2)) ) &
       &    * oosqrt2
    End Do
   End Do
  End Do

  DTerm = 0.25_dp * gp**2 * (v22 - v12 - v32 - v42 - v52)
  mat8(6:8,6:8) = mR2 + 0.5_dp * Matmul(yukV,yukVA)
  vec = Matmul(vevL,yuk)
  Do i1=1,3
   mat8(i1+5,i1+5) = mat8(i1+5,i1+5) + Dterm
   Do i2=i1,3
    mat8(i1+5,i2+5) = mat8(i1+5,i2+5) + 0.5_dp * vec(i1) * Conjg(vec(i2))
   End Do
  End Do

  Do i1=2,8
   Do i2 = 1,i1-1
    mat8(i1,i2) = Conjg( mat8(i2,i1) )
   End Do
  End Do

!  write (*,*) 'mSpm'
!  write (*,800) Real(mat8(1,:))
!  write (*,800) Real(mat8(2,:))
!  write (*,800) Real(mat8(3,:))
!  write (*,800) Real(mat8(4,:))
!  write (*,800) Real(mat8(5,:))
!  write (*,800) Real(mat8(6,:))
!  write (*,800) Real(mat8(7,:))
!  write (*,800) Real(mat8(8,:))
!  write (*,*) ' '
! 800 format(8e13.5)

  Call EigenSystem(mat8,mSpm2,RSpm,ierr, test)
  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If

  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'In subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'subroutine EigenSystem gives ierr = ',ierr
   Write(ErrCan,*) 'MH1sq ',MH1sq
   Write(ErrCan,*) 'MH2sq ',MH2sq
   Write(ErrCan,*) 'bi ',bi
   Write(ErrCan,*) 'Bmu ',Bmu
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) 'g ',g
   Write(ErrCan,*) 'yuk ',yuk
   Write(ErrCan,*) 'ML2(1,i1) ',(ML2(1,i1),i1=1,3)
   Write(ErrCan,*) 'ML2(2,i1) ',(ML2(2,i1),i1=1,3)
   Write(ErrCan,*) 'ML2(3,i1) ',(ML2(3,i1),i1=1,3)
   Write(ErrCan,*) 'MR2(1,i1) ',(MR2(1,i1),i1=1,3)
   Write(ErrCan,*) 'MR2(2,i1) ',(MR2(2,i1),i1=1,3)
   Write(ErrCan,*) 'MR2(3,i1) ',(MR2(3,i1),i1=1,3)
   Write(ErrCan,*) 'A(1,i1) ',(A(1,i1),i1=1,3)
   Write(ErrCan,*) 'A(2,i1) ',(A(2,i1),i1=1,3)
   Write(ErrCan,*) 'A(3,i1) ',(A(3,i1),i1=1,3)
   Write(ErrCan,*) '  '
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,8
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
      Write(ErrCan,*) 'MH1sq ',MH1sq
      Write(ErrCan,*) 'MH2sq ',MH2sq
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'Bmu ',Bmu
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'g ',g
      Write(ErrCan,*) 'yuk ',yuk
      Write(ErrCan,*) 'ML2(1,i2) ',(ML2(1,i2),i2=1,3)
      Write(ErrCan,*) 'ML2(2,i2) ',(ML2(2,i2),i2=1,3)
      Write(ErrCan,*) 'ML2(3,i2) ',(ML2(3,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(1,i2) ',(MR2(1,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(2,i2) ',(MR2(2,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(3,i2) ',(MR2(3,i2),i2=1,3)
      Write(ErrCan,*) 'A(1,i2) ',(A(1,i2),i2=1,3)
      Write(ErrCan,*) 'A(2,i2) ',(A(2,i2),i2=1,3)
      Write(ErrCan,*) 'A(3,i2) ',(A(3,i2),i2=1,3)
      Write(ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -203
    Call AddError(203)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mSpm(1) = mW
   mSpm2(1) = mW**2
   Do i1=2,8
    If (mSpm2(i1).Ge.0._dp) Then
     mSpm(i1) = Sqrt( mSpm2(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mSpm2(i1) 
      Write(ErrCan,*) 'MH1sq ',MH1sq
      Write(ErrCan,*) 'MH2sq ',MH2sq
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'Bmu ',Bmu
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'g ',g
      Write(ErrCan,*) 'yuk ',yuk
      Write(ErrCan,*) 'ML2(1,i2) ',(ML2(1,i2),i2=1,3)
      Write(ErrCan,*) 'ML2(2,i2) ',(ML2(2,i2),i2=1,3)
      Write(ErrCan,*) 'ML2(3,i2) ',(ML2(3,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(1,i2) ',(MR2(1,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(2,i2) ',(MR2(2,i2),i2=1,3)
      Write(ErrCan,*) 'MR2(3,i2) ',(MR2(3,i2),i2=1,3)
      Write(ErrCan,*) 'A(1,i2) ',(A(1,i2),i2=1,3)
      Write(ErrCan,*) 'A(2,i2) ',(A(2,i2),i2=1,3)
      Write(ErrCan,*) 'A(3,i2) ',(A(3,i2),i2=1,3)
      Write(ErrCan,*) '  '
     End If
     mSpm(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -203
    Call AddError(203)
    End If
   End Do
  End If

  Iname = Iname - 1

 End Subroutine ChargedScalarMassLam3nt

 Subroutine CharginoMass2(M2,mu,vevs,g,mC,U,V,kont)
 !-----------------------------------------------------------------
 ! calculates chargino masses + mixing matrices U,V in the MSSM,
 ! input:
 !  M2 ........ SU(2) gaugino mass
 !  mu ........ mu-parameter
 !  vevs(i) ... i=1 v_d
 !              i=2 v_u
 !  g ......... SU(2) gauge coupling
 ! written by Werner Porod, 2.10.2000
 !----------------------------------------------------------------------
 Implicit None

  Complex(Dp), Intent(in) :: M2,mu
  Real(Dp), Intent(in) :: vevs(2), g
  Complex(Dp), Intent(out) :: V(2,2), U(2,2)
  Real(Dp), Intent(out) :: mC(2)
  Integer, Intent(inout) :: kont

  Real(dp) :: u1(2,2), V1(2,2), mC2(2), test(2)
  Complex(Dp) :: mat2(2,2)=0, mat22(2,2)=0, U2(2,2), V2(2,2), phaseM
  Integer :: ierr, i1

  Iname = Iname + 1
  NameOfUnit(Iname) = "CharginoMass2"

  mC = 0
  U = 0
  V = 0

  mat2(1,1) = M2
  mat2(1,2) = g * vevs(2) * oosqrt2
  mat2(2,1) = g * vevs(1) * oosqrt2
  mat2(2,2) = mu

  mat22 = Matmul( Transpose( Conjg( mat2 ) ), mat2 )
  If ( Maxval( Abs( Aimag(mat22) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat22,dp), mC2, v1, ierr, test)
   v2 = v1
  Else
   Call EigenSystem(mat22, mC2, v2, ierr, test)
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat22 = Matmul( mat2, Transpose( Conjg( mat2 ) ) )
  If ( Maxval( Abs( Aimag(mat22) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat22,dp), mC2, u1, ierr, test)
   u2 = u1
  Else
   Call EigenSystem(mat22, mC2, u2, ierr, test)
  End If
  u2 = Conjg(u2)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat22 = Matmul( Matmul( Conjg(u2), mat2), Transpose( Conjg(v2) ) )
  Do i1=1,2
   phaseM =   mat22(i1,i1)   / Abs( mat22(i1,i1) )
   v2(i1,:) = phaseM * v2(i1,:)
  End Do

  If (Abs(v2(1,2)).Gt.0._dp) Then
   phaseM = v2(1,2) / Abs(v2(1,2))
   v2(1,2) = Abs(v2(1,2))
   v2(1,1) = v2(1,1) * Conjg(phaseM)
   u2(1,1) = u2(1,1) * phaseM
   u2(1,2) = u2(1,2) * phaseM
  End If

  If (Abs(v2(2,1)).Gt.0._dp) Then
   phaseM = v2(2,1) / Abs(v2(2,1))
   v2(2,1) = Abs(v2(2,1))
   v2(2,2) = v2(2,2) * Conjg(phaseM)
   u2(2,2) = u2(2,2) * phaseM
   u2(2,1) = u2(2,1) * phaseM
  End If

!  If (WriteOut) then
!   mat22 = Matmul( Matmul( Conjg(u2), mat2), Transpose( Conjg(v2) ) )
!   Write(ErrCan,*) "test",sqrt(mc2)
!   Write(ErrCan,*) "    ",mat22(1,:)
!   Write(ErrCan,*) "    ",mat22(2,:)
!   Write(errcan,*) "ratios",1._dp - mat22(1,1)/Sqrt(mc2(1)) &
!                 &         ,1._dp - mat22(2,2)/Sqrt(mc2(2))
!   Write(ErrCan,*) " "
!  End if

  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'Warning in subroutine CharginoMass2, ierr = ',ierr
   Write(ErrCan,*) 'M2,g ',M2,g
   Write(ErrCan,*) 'mu ',mu
   Write(ErrCan,*) 'vevs ',vevs
   Write(ErrCan,*) ' '
   kont = ierr
   mC(1) = Abs(M2)
   mC(2) = Abs(mu)
   U = id2C
   V = U
   Iname = Iname - 1
   Return
  Endif

  mC = Sqrt( mC2 )
  U = u2
  V = v2

  Iname = Iname - 1

 End Subroutine CharginoMass2


 Subroutine CharginoMass3a(M2,bi,vevSM,vevL,g,mTau,mC,U,V,htau,kont)
 Implicit None

  Complex(Dp), Intent(in) :: M2,bi(2)
  Real(Dp), Intent(in) :: vevSM(2),vevL(:),g,mTau
  Complex(Dp), Intent(out) :: V(3,3),U(3,3)
  Real(Dp), Intent(out) :: mC(3),htau
  Integer, Intent(inout) :: kont

  Real(Dp) vL1
  
  vL1 = vevL(1)
  Call CharginoMass3(M2,bi,vevSM,vL1,g,mTau,mC,U,V,htau,kont)

 End Subroutine CharginoMass3a

 Subroutine CharginoMass3(M2,bi,vevSM,vevL,g,mTau,mC,U,V,htau,kont)
 !-----------------------------------------------------------------
 ! calculates chargino masses + mixing matrices U,V in the MSSM,
 ! input:
 !  M2 ........ SU(2) gaugino mass
 !  bi ........ i=1 mu-parameter
 !              i=2 epsilon_3
 !  vevSM(i) .. i=1 v_d
 !              i=2 v_u
 !  vevL ...... v_3
 !  g ......... SU(2) gauge coupling
 ! written by Werner Porod, 2.10.2000
 !----------------------------------------------------------------------
 Implicit None
  Complex(Dp), Intent(in) :: M2,bi(2)
  Real(Dp), Intent(in) :: vevSM(2),vevL,g,mTau
  Complex(Dp), Intent(out) :: V(3,3),U(3,3)
  Real(Dp), Intent(out) :: mC(3),htau
  Integer, Intent(inout) :: kont

  Real(Dp) :: mu2,eps2,M22,mT2,g2,vD2,vU2,v32,g4,nen1,zaehl1,REM2eps, &
       &    ReM2mu,ReCepsmu,ht2,v3a,vD,vU, v1(3,3), u1(3,3), mC2(3), test(2)
  Complex(Dp) :: mat3(3,3)=0, mat32(3,3)=0, U3(3,3), V3(3,3), phaseM
  Integer :: ierr,i1

  Iname = Iname + 1
  NameOfUnit(Iname) = "CharginoMass3"

  mC = 0
  U = 0
  V = 0

  mT2 = mTau**2
  mu2 = Abs( bi(1) )**2
  eps2 = Abs( bi(2) )**2
  M22 = Abs( M2 )**2
  g2 = g*g
  g4 = g2 * g2
  vD = vevSM(1)
  vU = vevSM(2)
  v3a = vevL
  vD2 = vevSM(1)**2
  vU2 = vevSM(2)**2
  v32 = vevL**2
  ReM2eps = Real( M2*bi(2),dp )
  ReM2mu = Real( M2*bi(1),dp )
  ReCepsmu = Real( Conjg(bi(2)) * bi(1),dp )
  zaehl1 = mT2 * ( 4._dp * mT2* (mT2 - eps2 - mu2 - M22 )                    &
     &     - 2._dp * g2 * mT2 * ( v32 + vU2 + vD2 )                          &
     &     + g4 * vU2 * (v32 + vD2 )                                         &
     &     + 4._dp * g2 * vU * ( v3a * ReM2eps - vD * ReM2mu )                &
     &     + 2._dp * g2 * ( vD2 * eps2 + v32 *  mu2 + 2 * v3a * vD * ReCepsmu)&
     &     + 4._dp * M22 * ( eps2 + mu2 ) )
                      
  nen1 = (v32 + vD2)                                                         &
     & * ( 2._dp * mT2 * (2*mT2 - g2 * (v32 + vU2 + vD2) - 2*M22 )            &
     &   + g4 * vU2 * (v32 + vD2 )                                           &
     &   + 4._dp * g2 * vU * ( v3a * ReM2eps - vD * ReM2mu ) )                &
     &   - 4._dp*(mT2 - M22) * (v32* eps2 + vD2* mu2 - 2._dp* v3a*vD*ReCepsmu) 
      
  ht2 = 2._dp * zaehl1 / nen1
  If (ht2.Gt.0._dp) Then
   htau = Sqrt( ht2 )
  Else
   Write(ErrCan,*) 'Servere Warning from routine CharginoMass3 !'
   Write(ErrCan,*) 'Abs(h_tau)**2 < 0 :',ht2
   Write(ErrCan,*) 'Taking the square root from the negative !'
   Write(ErrCan,*) 'The Parameters are :'
   Write(ErrCan,*) ' g_SU(2) ',g
   Write(ErrCan,*) ' vevSM    ',vevSM
   Write(ErrCan,*) ' vevL    ',vevL
   Write(ErrCan,*) ' mu,eps_3 ',bi
   Write(ErrCan,*) ' M2 ',M2
   kont = -204
   Call AddError(204)
   htau = Sqrt( -ht2 )
  Endif

  mat3(1,1) = M2
  mat3(1,2) = g * vU * oosqrt2
  mat3(1,3) = 0._dp
  mat3(2,1) = g * vD * oosqrt2
  mat3(2,2) = bi(1)
  mat3(2,3) = - htau * v3a * oosqrt2
  mat3(3,1) = g * v3a * oosqrt2
  mat3(3,2) = -bi(2)
  mat3(3,3) = htau * vD * oosqrt2

  mat32 = Matmul( Transpose( Conjg( mat3 ) ), mat3 )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, v1, ierr, test)
   v3 = v1
  Else
   Call EigenSystem(mat32, mC2, v3, ierr, test)
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( mat3, Transpose( Conjg( mat3 ) ) )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, u1, ierr, test)
   u3 = u1
  Else
   Call EigenSystem(mat32, mC2, u3, ierr, test)
  End If
  u3 = Conjg(u3)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( Matmul( Conjg(u3), mat3), Transpose( Conjg(v3) ) )
  Do i1=1,3
   phaseM =   mat32(i1,i1)   / Abs( mat32(i1,i1) )
   v3(i1,:) = phaseM * v3(i1,:)
  End Do

  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'Warning in subroutine CharginoMass3, ierr = ',ierr
   Write(ErrCan,*) 'M2,g ',M2,g
   Write(ErrCan,*) 'mu,eps ',bi
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) ' '
   kont = ierr
   Iname = Iname - 1
   Return
  Endif

  mC = Sqrt( mC2 )
  U = u3
  V = v3

 Iname = Iname - 1

 End Subroutine CharginoMass3

 Subroutine CharginoMass5(M2,bi,vevSM,vevL,g,mL,mC,U,V,yuk_L,kont)
 !-----------------------------------------------------------------
 ! calculates chargino masses + mixing matrices U,V in the  3-generation 
 ! epsilon model. In addition the lepton yukawas are calculated.
 ! input:
 !  M2 ........ SU(2) gaugino mass
 !  bi(i) ..... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  vevSM(i) .. i=1 v_d
 !              i=2 v_u
 !  vevL(i) ... i=1 v_1
 !              i=2 v_2
 !              i=3 v_3
 !  g ......... SU(2) gauge coupling
 !  mL(i) ..... i=1 electron mass
 !              i=2 muon mass
 !              i=3 tau mass
 ! written by Werner Porod, 2.10.2000
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2),vevL(3),g,mL(3)
  Real(Dp), Intent(out) :: mC(5),yuk_L(3)
  Complex(Dp), Intent(in) :: M2,bi(4)
  Complex(Dp), Intent(out) :: U(5,5),V(5,5)

  Integer :: i1,ierr
  Complex(Dp) :: mat5(5,5),mat52(5,5), U5(5,5), V5(5,5), PhaseM
  Real(Dp) :: mu2,eps2,M22,mT2,g2,vD2,vU2,v32,g4,nen1,zaehl1,REM2eps, &
       &  ReM2mu,ReCepsmu,ht2,v3,vD,vU, U1(5,5), V1(5,5), mC2(5), test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = "CharginoMass5"

  mC = 0
  U = 0
  V = 0

  mu2 = Abs( bi(1) )**2
  eps2 = Abs( bi(2) )**2
  M22 = Abs( M2 )**2
  g2 = g*g
  g4 = g2 * g2
  vD = vevSM(1)
  vU = vevSM(2)
  vD2 = vevSM(1)**2
  vU2 = vevSM(2)**2
  ReM2mu = Real( M2*bi(1),dp )
  Do i1=1,3
   v3 = vevL(i1)
   v32 = v3**2
   mT2 = mL(i1)**2
   eps2 = Abs( bi(i1+1) )**2
   ReM2eps = Real( M2*bi(i1+1),dp )
   ReCepsmu = Real( Conjg(bi(i1+1)) * bi(1),dp )
      
   zaehl1 = mT2 * ( 4._dp * mT2* (mT2 - eps2 - mu2 - M22 )                    &
     &     - 2._dp * g2 * mT2 * ( v32 + vU2 + vD2 )                           &
     &     + g4 * vU2 * (v32 + vD2 )                                          &
     &     + 4._dp * g2 * vU * ( v3 * ReM2eps - vD * ReM2mu )                 &
     &     + 2._dp * g2 * ( vD2 * eps2 + v32 *  mu2 + 2 * v3 * vD * ReCepsmu) &
     &     + 4._dp * M22 * ( eps2 + mu2 ) )
                      
   nen1 = (v32 + vD2)                                                        &
     & * ( 2._dp * mT2 * (2*mT2 - g2 * (v32 + vU2 + vD2) - 2*M22 )           &
     &   + g4 * vU2 * (v32 + vD2 )                                           &
     &   + 4._dp * g2 * vU * ( v3 * ReM2eps - vD * ReM2mu ) )                 &
     &   - 4._dp*(mT2 - M22) * (v32* eps2 + vD2* mu2 - 2._dp* v3* vD*ReCepsmu) 
      
   ht2 = 2._dp * zaehl1 / nen1
   If (ht2.Gt.0._dp) Then
    yuk_L(i1) = Sqrt( ht2 )
   Else
    Write(ErrCan,*) 'Servere Warning from routine CharginoMass5!!!'
    Write(ErrCan,*) 'Abs(h_tau)**2 < 0 :',i1,ht2
    Write(ErrCan,*) 'Taking the square root from the negative!!!'
    Write(ErrCan,*) 'The Parameters are :'
    Write(ErrCan,*) ' g_SU(2) ',g
    Write(ErrCan,*) ' vevSM    ',vevSM
    Write(ErrCan,*) ' vevL    ',vevL
    Write(ErrCan,*) ' mu,eps_i ',bi
    Write(ErrCan,*) ' M2 ',M2
    kont = -205
    Call AddError(205)
    yuk_L(i1) = Sqrt( -ht2 )
   Endif
  Enddo

  mat5(1,1) = M2
  mat5(1,2) = g * vevSM(2) * oosqrt2
  mat5(1,3) = 0._dp
  mat5(1,4) = 0._dp
  mat5(1,5) = 0._dp
  mat5(2,1) = g * vevSM(1) * oosqrt2
  mat5(2,2) = bi(1)
  mat5(2,3) = - yuk_L(1) * vevL(1) * oosqrt2
  mat5(2,4) = - yuk_L(2) * vevL(2) * oosqrt2
  mat5(2,5) = - yuk_L(3) * vevL(3) * oosqrt2
  mat5(3,1) = g * vevL(1) * oosqrt2
  mat5(3,2) = -bi(2)
  mat5(3,3) = yuk_L(1) * vevSM(1) * oosqrt2
  mat5(3,4) = 0._dp
  mat5(3,5) = 0._dp
  mat5(4,1) = g * vevL(2) * oosqrt2
  mat5(4,2) = -bi(3)
  mat5(4,3) = 0._dp
  mat5(4,4) = yuk_L(2) * vevSM(1) * oosqrt2
  mat5(4,5) = 0._dp
  mat5(5,1) = g * vevL(3) * oosqrt2
  mat5(5,2) = -bi(4)
  mat5(5,3) = 0._dp
  mat5(5,4) = 0._dp
  mat5(5,5) = yuk_L(3) * vevSM(1) * oosqrt2

  mat52 = Matmul( Transpose( Conjg( mat5 ) ), mat5 )
  If ( Maxval( Abs( Aimag(mat52) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat52,dp), mC2, v1, ierr, test)
   v5 = v1
  Else
   Call EigenSystem(mat52, mC2, v5, ierr, test)
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat52 = Matmul( mat5, Transpose( Conjg( mat5 ) ) )
  If ( Maxval( Abs( Aimag(mat52) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat52,dp), mC2, u1, ierr, test)
   u5 = u1
  Else
   Call EigenSystem(mat52, mC2, u5, ierr, test)
  End If
  u5 = Conjg(u5)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat52 = Matmul( Matmul( Conjg(u5), mat5), Transpose( Conjg(v5) ) )
  Do i1=1,5
   phaseM =   mat52(i1,i1)   / Abs( mat52(i1,i1) )
   v5(i1,:) = phaseM * v5(i1,:)
  End Do

 If (ierr.Ne.0) Then
   Write(ErrCan,*) 'Warning in subroutine CharginoMass5, ierr = ',ierr
   Write(ErrCan,*) 'M2,g ',M2,g
   Write(ErrCan,*) 'mu,eps ',bi
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) ' '
   kont = ierr
   Iname = Iname - 1
   Return
  Endif

  mC = Sqrt( mC2 )
  U = u5
  V = v5

  Iname = Iname - 1

 End Subroutine CharginoMass5


 Subroutine CharginoMass5Lam(M2,bi,vevSM,vevL,g,yuk_L,lam,mL,mC,U,V,kont)
 !-----------------------------------------------------------------
 ! calculates chargino masses + mixing matrices U,V in the  3-generation 
 ! epsilon model.
 ! input:
 !  M2 ........ SU(2) gaugino mass
 !  bi(i) ..... i=1 mu-parameter
 !              i=2 epsilon_1
 !              i=3 epsilon_2
 !              i=4 epsilon_3
 !  vevSM(i) .. i=1 v_d
 !              i=2 v_u
 !  vevL(i) ... i=1 v_1
 !              i=2 v_2
 !              i=3 v_3
 !  g ......... SU(2) gauge coupling
 !  mL(i) ..... i=1 electron mass
 !              i=2 muon mass
 !              i=3 tau mass
 ! written by Werner Porod, 25.08.2001
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, mL(3)
  Real(Dp), Intent(out) :: mC(5)
  Complex(Dp), Intent(in) :: M2, bi(4), lam(3,3,3)
  Complex(Dp), Intent(inout) :: yuk_L(3,3)
  Complex(Dp), Intent(out) :: U(5,5), V(5,5)

  Integer :: i1,i2,i3,ierr, i_count
  Complex(Dp) :: mat5(5,5),mat52(5,5),phaseM,  U5(5,5),V5(5,5)
  Real(Dp) :: U1(5,5), V1(5,5), mC2(5), test(2)
!  Real(Dp) :: mu2,eps2,M22,mT2,g2,vD2,vU2,v32,g4,nen1,zaehl1,REM2eps, &
!       &           ReM2mu,ReCepsmu,ht2,v3,vD,vU

  Iname = Iname + 1
  NameOfUnit(Iname) = "CharginoMass5Lma"

  mC = 0
  U = 0
  V = 0

  i_count = 1
  Do 
   mat5(1,1) = M2
   mat5(1,2) = g * vevSM(2) * oosqrt2
   mat5(1,3) = 0._dp
   mat5(1,4) = 0._dp
   mat5(1,5) = 0._dp
   mat5(2,1) = g * vevSM(1) * oosqrt2
   mat5(2,2) = bi(1)
   mat5(2,3:5) = - Matmul(vevL,yuk_L)  * oosqrt2
   mat5(3:5,1) = g * vevL * oosqrt2
   mat5(3:5,2) = -bi(2:4)
   mat5(3:5,3:5) = yuk_L * vevSM(1) * oosqrt2
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      mat5(2+i1,2+i2) = mat5(2+i1,2+i2)  &
                    & +(lam(i3,i1,i2) - lam(i1,i3,i2)) * vevL(i3) * oosqrt2
     End Do
    End Do
   End Do

   mat52 = Matmul( Transpose( Conjg( mat5 ) ), mat5 )
   If ( Maxval( Abs( Aimag(mat52) ) ).Eq.0._dp) Then ! reel matrix
    Call EigenSystem(Real(mat52,dp), mC2, v1, ierr, test)
    v5 = v1
   Else
    Call EigenSystem(mat52, mC2, v5, ierr, test)
   End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
   mat52 = Matmul( mat5, Transpose( Conjg( mat5 ) ) )
   If ( Maxval( Abs( Aimag(mat52) ) ).Eq.0._dp) Then ! reel matrix
    Call EigenSystem(Real(mat52,dp), mC2, u1, ierr, test)
    u5 = u1
   Else
    Call EigenSystem(mat52, mC2, u5, ierr, test)
   End If
   u5 = Conjg(u5)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
   mat52 = Matmul( Matmul( Conjg(u5), mat5), Transpose( Conjg(v5) ) )
   Do i1=1,5
    phaseM =   mat52(i1,i1)   / Abs( mat52(i1,i1) )
    v5(i1,:) = phaseM * v5(i1,:)
   End Do

   If (ierr.Ne.0) Then
    Write(ErrCan,*) 'Warning in subroutine CharginoMass5Lam, ierr = ',ierr
    Write(ErrCan,*) 'M2,g ',M2,g
    Write(ErrCan,*) 'mu,eps ',bi
    Write(ErrCan,*) 'vevSM ',vevSM
    Write(ErrCan,*) 'vevL ',vevL
    Write(ErrCan,*) ' '
    kont = ierr
    Iname = Iname - 1
    Return
   Endif

   mC = Sqrt( mC2 )
   ! a numerical problem occurs here, because m_e^2 /mC(2)^2 is in the order
   ! of 10^-14
   If (      ( (Abs(mL(2)-Abs(mC(2)))/mL(2)).Lt.1.e-3_dp)          &
      & .And.( (Abs(mL(3)-Abs(mC(3)))/mL(3)).Lt.1.e-3_dp) ) Exit
    yuk_L(1,:) = yuk_L(1,:) * mL(1) / mC(1)
    yuk_L(2,:) = yuk_L(2,:) * mL(2) / mC(2)
    yuk_L(3,:) = yuk_L(3,:) * mL(3) / mC(3)
   i_count = i_count + 1
   If (i_count.Gt.50) Exit
  End Do

  U = u5
  V = v5

  Do i1=1,3
   If(Abs(mL(i1)-Abs(mC(i1)))/mL(i1).Gt.0.01) Then
    Write(ErrCan,*) 'Problem with lepton masses in Subroutine CharginoMass5Lam:'
    Write(ErrCan,*) 'm_l(i1), mC(i1) :',mL(i1), mC(i1)
   End If
  End Do

  Iname = Iname - 1

 End Subroutine CharginoMass5Lam


 Subroutine NeutralinoMass4(M1,M2,mu,vevs,gp,g,mN,N,kont)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the MSSM. 
 ! input:
 !  M1 .......... U(1) gaugino mass
 !  M2 .......... SU(2) gaugino mass
 !  mu .......... mu-parameter
 !  vevs(i) ..... i=1 v_d
 !                i=2 v_u
 !  g ........... SU(2) gauge coupling
 !  gp .......... U(1) gauge coupling
 ! written by Werner Porod, 6.10.2000
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevs(2),g,gp
  Real(Dp), Intent(out) :: mN(4)
  Complex(Dp), Intent(in) :: M1,M2,mu
  Complex(Dp), Intent(out) :: N(4,4)

  Integer :: i1,i2,ierr
  Complex(Dp) :: mat4(4,4), mat42(4,4), E4(4), phaseM
  Real(Dp) :: g1, gp1, N4a(4,4), eig(4), test(2)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass4'

  g1 = 0.5_dp * g 
  gp1 = 0.5_dp * gp 

  mN = 0._dp
  N = (0._dp,0._dp)

  mat4(1,1) = M1
  mat4(1,2) = 0._dp
  mat4(1,3) = - gp1 * vevs(1)
  mat4(1,4) = gp1 * vevs(2)
  mat4(2,2) =  M2
  mat4(2,3) = g1 * vevs(1)
  mat4(2,4) = - g1 * vevs(2)
  mat4(3,3) = 0._dp
  mat4(3,4) = - mu
  mat4(4,4) = 0._dp

  Do i1=2,4
   Do i2=1,i1-1
    mat4(i1,i2) = mat4(i2,i1)
   Enddo
  Enddo

  If (Maxval(Abs(Aimag(mat4))).Eq.0._dp) Then ! matrix is reel
   Call EigenSystem(Real(mat4,dp), Eig, N4a, ierr, test)

   Do i1=1,4
    If (Eig(i1).Lt.0._dp) Then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N4a(i1,:)
    Else
     mN(i1) = Eig(i1)
     N(i1,:) =N4a(i1,:)
    End If
   End Do

   Do i1=1,3
    Do i2=i1+1,4
     If (mN(i1).Gt.mN(i2)) Then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E4 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E4
     End If
    End Do
   End Do

  Else

   mat42 = Matmul( Transpose(Conjg( mat4 ) ), mat4 )
   Call EigenSystem(mat42, Eig, N, ierr, test)
   mat42 = Matmul(Conjg(N), Matmul( mat4, Transpose( Conjg( N ) ) ) )
   Do i1=1,4
    phaseM =   Sqrt( mat42(i1,i1)   / Abs( mat42(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN = Sqrt( Eig )

  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (Errcan,*) 'Warning in subroutine NeutralinoMass4, ierr =',ierr
   Write (Errcan,*) 'M1,M2 ',M1,M2
   Write (Errcan,*) 'gp,g ',gp,g
   Write (Errcan,*) 'mu ',mu
   Write (Errcan,*) 'vevs ',vevs
   Write (Errcan,*) ' '
   kont = ierr
   Iname = Iname - 1
   Return
  Endif

  Iname = Iname - 1

 End Subroutine NeutralinoMass4


 Subroutine NeutralinoMass5a(M1,M2,bi,vevSM,vevL,g,gp,mN,N,kont)
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2),vevL(:),g,gp
  Real(Dp), Intent(out) :: mN(5)
  Complex(Dp), Intent(in) :: M1,M2,bi(2)
  Complex(Dp), Intent(inout) :: N(5,5)

  Real(Dp) vL1
 
  vL1 = vevL(1)
  Call NeutralinoMass5(M1,M2,bi,vevSM,vL1,g,gp,mN,N,kont)

 End Subroutine NeutralinoMass5a

 Subroutine NeutralinoMass5(M1,M2,bi,vevSM,vevL,g,gp,mN,N,kont)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the  1-generation epsilon
 ! input:
 !  M1 .......... U(1) gaugino mass
 !  M2 .......... SU(2) gaugino mass
 !  bi(i) ....... i=1 mu-parameter
 !                i=2 epsilon_3
 !  vevSM(i) ..... i=1 v_d
 !                 i=2 v_u
 !  vevL ......... v_3
 !  g ........... SU(2) gauge coupling
 !  gp .......... U(1) gauge coupling
 ! written by Werner Porod, 2.8.99
 ! 2.1.00: changing diagonalization from subroutine eisch1 to subroutine csvdc
 ! 1.11.2000: portation to f90
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2),vevL,g,gp
  Real(Dp), Intent(out) :: mN(5)
  Complex(Dp), Intent(in) :: M1,M2,bi(2)
  Complex(Dp), Intent(inout) :: N(5,5)

  Integer :: i1,i2,ierr
  Complex(Dp) :: mat5(5,5), mat52(5,5), phaseM, E5(5)
  Real(Dp) :: g1, gp1, N5a(5,5), eig(5), test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass5'

  g1 = g / 2._dp
  gp1 = gp / 2._dp

  mN = 0._dp
  N = (0._dp,0._dp)

  mat5(1,1) = M1
  mat5(1,2) = 0._dp
  mat5(1,3) = - gp1 * vevSM(1)
  mat5(1,4) = gp1 * vevSM(2)
  mat5(1,5) = - gp1 * vevL
  mat5(2,2) =  M2
  mat5(2,3) = g1 * vevSM(1)
  mat5(2,4) = - g1 * vevSM(2)
  mat5(2,5) = g1 * vevL
  mat5(3,3) = 0._dp
  mat5(3,4) = - bi(1)
  mat5(3,5) = 0._dp
  mat5(4,4) = 0._dp
  mat5(4,5) = bi(2)
  mat5(5,5) = 0._dp

  Do i1=2,5
   Do i2=1,i1-1
    mat5(i1,i2) = mat5(i2,i1)
   Enddo
  Enddo

  If (Maxval(Abs(Aimag(mat5))).Eq.0._dp) Then ! matrix is reel

   Call EigenSystem(Real(mat5,dp), Eig, N5a, ierr, test)

   Do i1=1,5
    If (Eig(i1).Lt.0._dp) Then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N5a(i1,:)
    Else
     mN(i1) = Eig(i1)
     N(i1,:) =N5a(i1,:)
    End If
   End Do

   Do i1=1,4
    Do i2=i1+1,5
     If (mN(i1).Gt.mN(i2)) Then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E5 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E5
     End If
    End Do
   End Do

  Else

   mat52 = Matmul( Transpose(Conjg( mat5 ) ), mat5 )
   Call EigenSystem(mat52, Eig, N, ierr, test)
   mat52 = Matmul(Conjg(N), Matmul( mat5, Transpose( Conjg( N ) ) ) )
   Do i1=1,5
    phaseM =   Sqrt( mat52(i1,i1)   / Abs( mat52(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN = Sqrt( Eig )

  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (ErrCan,*) 'Warning in subroutine NeutralinoMass5, ierr =',ierr
   Write (ErrCan,*) 'M1,M2 ',M1,M2
   Write (ErrCan,*) 'gp,g ',gp,g
   Write (ErrCan,*) 'bi ',bi
   Write (ErrCan,*) 'vevSM ',vevSM
   Write (ErrCan,*) 'vevL ',vevL
   Write (ErrCan,*) ' '
   kont = ierr
  Endif

  Iname = Iname - 1

 End Subroutine NeutralinoMass5

 Subroutine NeutralinoMass5NMSSM(M1, M2, mu, gp, g, h0, lam, vevSM, vevN &
                               &, mN, N, kont)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the  1-generation epsilon
 ! input:
 !  M1 .......... U(1) gaugino mass
 !  M2 .......... SU(2) gaugino mass
 !  mu ..........  mu-parameter
 !  gp .......... U(1) gauge coupling
 !  g ........... SU(2) gauge coupling
 !  h0 .......... trilinear H_d H_u phi coupling
 !  lam ......... trilinear self interaction of phi 
 !  vevSM(i) .... i=1 v_d
 !                i=2 v_u
 !  vevN ........ v_phi 
 ! output:
 !  mN(i) ....... vector with neutralino masses m_i<m_j for i<j, all positive
 !  N(i,j) ...... mixing matrix
 !  kont ........ is 0 if everything is fine, otherwise there has been a 
 !                numerical problem
 ! written by Werner Porod, 14.02.05
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2),g,gp,vevN
  Real(Dp), Intent(out) :: mN(5)
  Complex(Dp), Intent(in) :: M1,M2,mu,h0,lam
  Complex(Dp), Intent(inout) :: N(5,5)

  Integer :: i1, i2, ierr
  Complex(Dp) :: mat5(5,5), mat52(5,5), phaseM, E5(5)
  Real(Dp) :: g1, gp1, N5a(5,5), eig(5), test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass5NMSSM'

  kont = 0

  g1 = g / 2._dp
  gp1 = gp / 2._dp

  mN = 0._dp
  N = (0._dp,0._dp)

  mat5(1,1) = M1
  mat5(1,2) = 0._dp
  mat5(1,3) = - gp1 * vevSM(1)
  mat5(1,4) = gp1 * vevSM(2)
  mat5(1,5) = 0._dp
  mat5(2,2) =  M2
  mat5(2,3) = g1 * vevSM(1)
  mat5(2,4) = - g1 * vevSM(2)
  mat5(2,5) = 0._dp
  mat5(3,3) = 0._dp
  mat5(3,4) = - mu - h0 * vevN * oosqrt2
  mat5(3,5) =  - h0 * vevSM(2) * oosqrt2
  mat5(4,4) = 0._dp
  mat5(4,5) =  - h0 * vevSM(1) * oosqrt2
  mat5(5,5) = oosqrt2 * lam * vevN

  Do i1=2,5
   Do i2=1,i1-1
    mat5(i1,i2) = mat5(i2,i1)
   Enddo
  Enddo

  If (MaxVal(Abs(Aimag(mat5))).eq.0._dp) then ! matrix is reel

   Call EigenSystem(Real(mat5,dp), Eig, N5a, ierr, test)

   Do i1=1,5
    If (Eig(i1).lt.0._dp) then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N5a(i1,:)
    else
     mN(i1) = Eig(i1)
     N(i1,:) =N5a(i1,:)
    end if
   end do

   Do i1=1,4
    Do i2=i1+1,5
     if (mN(i1).gt.mN(i2)) then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E5 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E5
     end if
    end do
   end do

  else

   mat52 = Matmul( Transpose(Conjg( mat5 ) ), mat5 )
   Call EigenSystem(mat52, Eig, N, ierr, test)
   mat52 = Matmul(Conjg(N), Matmul( mat5, Transpose( Conjg( N ) ) ) )
   Do i1=1,5
    phaseM =   Sqrt( mat52(i1,i1)   / Abs( mat52(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   end do
   mN = Sqrt( Eig )

  end if

  If (ierr.Ne.0) Then
   If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
   End If
   Write (ErrCan,*) 'Warning in subroutine NeutralinoMass, ierr =',ierr
   Write (ErrCan,*) 'M1,M2 ',M1,M2
   Write (ErrCan,*) 'gp,g ',gp,g
   Write (ErrCan,*) "h0, lam",h0,lam
   Write (ErrCan,*) 'mu ',mu
   Write (ErrCan,*) 'vevSM ',vevSM
   Write (ErrCan,*) 'vevN ',vevN
   Write (ErrCan,*) ' '
   kont = ierr
  Endif

  Iname = Iname - 1

 End Subroutine NeutralinoMass5NMSSM


 Subroutine NeutralinoMass7(M1,M2,bi,vevSM,vevL,g,gp,mN,N,kont)
 !-----------------------------------------------------------------
 ! calculates neutralino masses + mixing matrix N in the 3-generation epsilon
 ! model
 ! input:
 !  M1 .......... U(1) gaugino mass
 !  M2 .......... SU(2) gaugino mass
 !  bi(i) ....... i=1 mu-parameter
 !                i=2 epsilon_1
 !                i=3 epsilon_2
 !                i=4 epsilon_3
 !  vevSM(i) .... i=1 v_d
 !                i=2 v_u
 !  vevL(i) ..... i=1 v_1
 !                i=2 v_2
 !                i=3 v_3
 !  g ........... SU(2) gauge coupling
 !  gp .......... U(1) gauge coupling
 ! written by Werner Porod, 2.8.99
 ! 2.1.00: changing diagonalization from subroutine eisch1 to subroutine csvdc
 ! 23.3.2000: first steps for the implementation of the full model
 ! 1.11.2000: portation to f90
 !----------------------------------------------------------------------
 Implicit None
  Integer, Intent(inout) :: kont
  Real(Dp), Intent(in) :: vevSM(2),vevL(3),g,gp
  Real(Dp), Intent(out) :: mN(7)
  Complex(Dp), Intent(in) :: M1, M2, bi(4)
  Complex(Dp), Intent(inout) :: N(7,7)

  Integer :: i1,i2,ierr
  Complex(Dp) :: mat7(7,7), mat72(7,7), phaseM, E7(7)
  Real(Dp) :: g1, gp1, N7a(7,7), eig(7), test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoMass7'

  g1 = g / 2._dp
  gp1 = gp / 2._dp

  mN = 0._dp
  N = (0._dp,0._dp)

  mat7  = 0._dp
  mat7(1,1) = M1
  mat7(1,3) = - gp1 * vevSM(1)
  mat7(1,4) = gp1 * vevSM(2)
  mat7(1,5) = - gp1 * vevL(1)
  mat7(1,6) = - gp1 * vevL(2)
  mat7(1,7) = - gp1 * vevL(3)
  mat7(2,2) =  M2
  mat7(2,3) = g1 * vevSM(1)
  mat7(2,4) = - g1 * vevSM(2)
  mat7(2,5) = g1 * vevL(1)
  mat7(2,6) = g1 * vevL(2)
  mat7(2,7) = g1 * vevL(3)
  mat7(3,4) = - bi(1)
  mat7(4,5) = bi(2)
  mat7(4,6) = bi(3)
  mat7(4,7) = bi(4)

  Do i1=2,7
   Do i2=1,i1-1
    mat7(i1,i2) = mat7(i2,i1)
   Enddo
  Enddo

  If (Maxval(Abs(Aimag(mat7))).Eq.0._dp) Then ! matrix is reel

   Call EigenSystem(Real(mat7,dp), Eig, N7a, ierr, test)

   Do i1=1,7
    If (Eig(i1).Lt.0._dp) Then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N7a(i1,:)
    Else
     mN(i1) = Eig(i1)
     N(i1,:) =N7a(i1,:)
    End If
   End Do

   Do i1=1,6
    Do i2=i1+1,7
     If (mN(i1).Gt.mN(i2)) Then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E7 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E7
     End If
    End Do
   End Do

  Else

   mat72 = Matmul( Transpose(Conjg( mat7 ) ), mat7 )
   Call EigenSystem(mat72, Eig, N, ierr, test)
   mat72 = Matmul(Conjg(N), Matmul( mat7, Transpose( Conjg( N ) ) ) )
   Do i1=1,7
    phaseM =   Sqrt( mat72(i1,i1)   / Abs( mat72(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   !---------------------------------------------------------------
   ! use abs to avoid problems mit numerical zeros
   !---------------------------------------------------------------
   mN = Sqrt( Abs(Eig) )

  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (ErrCan,*) 'Warning in subroutine NeutralinoMass7, ierr =',ierr
   Write (ErrCan,*) 'M1,M2 ',M1,M2
   Write (ErrCan,*) 'gp,g ',gp,g
   Write (ErrCan,*) 'bi ',bi
   Write (ErrCan,*) 'vevSM ',vevSM
   Write (ErrCan,*) 'vevL ',vevL
   Write (ErrCan,*) ' '
   kont = ierr
  Endif

  Iname = Iname - 1

 End Subroutine NeutralinoMass7

 Subroutine OneLoopTadpolesEps3(gU1, gSU2, vevSM, vevL, bi &
            & , mSup2, Y_u, A_u, RUsquark, mSdown2, Y_d, A_d, RDsquark  &
            & , tadpoles)
 !-----------------------------------------------------
 ! Calculates the 1-loop tadpoles in the epsilon model
 ! written by Werner Porod
 !  - 24.10.01: including squark/quark contributions
 !-----------------------------------------------------
 Use Couplings
 Implicit None

  Complex(dp), Intent(in) :: Y_d(3,3), Y_u(3,3), RUsquark(6,6), RDsquark(6,6) &
      & , A_u(3,3), A_d(3,3), bi(4)
  Real(dp), Intent(in) :: gU1, gSU2, mSup2(6), mSdown2(6), vevSM(2), vevL(3)
  Real(dp), Intent(out) :: tadpoles(5)

  Integer :: i1, i2, i3 
  Real(dp) :: e_d, e_u, id5r(5,5), A0m
  Complex(dp) :: sumI(5), coupC, coupRC, Rsf(2,2)
  !----------------------
  ! Inititalization 
  !----------------------
  tadpoles = 0._dp
  sumI = 0._dp
  e_d = -1._dp / 3._dp
  e_u = 2._dp / 3._dp
  id5R = 0._dp
  Do i1=1,5
   id5R(i1,i1) = 1._dp
  End Do 
  !----------------------
  ! quark contribution
  !-----------------------
  Do i1=1,3
   Call CoupFermionScalar(1, -0.5_dp, Y_d(i1,i1), id5R, coupC, coupRC)
   sumI(1) = 12._dp * mf_d(i1) * Real( CoupC,dp ) * A0(mf_d2(i1))
   Call CoupFermionScalar(2, 0.5_dp, Y_u(i1,i1), id5R, coupC, coupRC)
   sumI(2) = 12._dp * mf_u(i1) * Real( coupC,dp ) * A0(mf_u2(i1))
   tadpoles(1:2) = tadpoles(1:2) + sumI(1:2)
  End Do
  !----------------------
  ! squark contribution
  !-----------------------

  Do i1=1,3
   i2 = 2*i1-1
   Rsf = RUsquark(i2:i2+1,i2:i2+1)
   Do i2=1,2
    Write(ErrCan,*) 2*(i1 -1) +i2,mSup2(2*(i1 -1) +i2),A0( mSup2(2*(i1 -1) +i2) )
    A0m = A0( mSup2(2*(i1 -1) +i2) )
    sumI = 0._dp
    Do i3=1,5
     Call CoupScalarSfermion3(i3, i2, i2, id5R, 0.5_dp, e_u, Y_u(i1,i1), Rsf  &
                            &, A_u(i1,i1), bi, vevSM, vevL, gU1, gSU2, coupC )
     sumI(i3) = - 3._dp * coupC * A0m
    End Do
    tadpoles = tadpoles + sumI
   End Do
  End Do

  Do i1=1,3
   i2 = 2*i1-1
   Rsf = RDsquark(i2:i2+1,i2:i2+1)
   Do i2=1,2
    A0m = A0( mSdown2(2*(i1 -1) +i2) )
    sumI = 0._dp
    Do i3=1,5
     Call CoupScalarSfermion3(i3, i2, i2, id5R, -0.5_dp, e_d, Y_d(i1,i1), Rsf &
                            &, A_d(i1,i1), bi, vevSM, vevL, gU1, gSU2, coupC )
     sumI(i3) = - 3._dp * coupC * A0m
    End Do
    tadpoles = tadpoles + sumI
   End Do
  End Do

  tadpoles = oo16pi2 * tadpoles

 Contains

  Real(dp) Function A0(m2)
   Implicit None
   Real(dp), Intent(in) :: m2
   A0 = - m2 * ( Log(m2/mZ2) - 1)
  End Function a0

 End Subroutine OneLoopTadpolesEps3


 Subroutine PseudoScalarMassEps1(bi,B,vevSM,vevL,gp,g,mP0,mP02,RP0,kont)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the 1-generation epsilon model
 ! Input:
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL ......... v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL, g, gp
  Real(Dp), Intent(out) :: mP0(3), mP02(3), RP0(3,3)
  Complex(Dp), Intent(in) ::  bi(2), B(2)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, ML2a, DTerm
  Complex(Dp) :: sumC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassEps1'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2 + vevL**2) 

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * bi(2),dp) * vevL ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1) + Real(B(2),dp) * vevL) / vevSM(2)
   
  sumC = - vevSM(2) * B(2) + vevSM(1) * Conjg(bi(1)) * bi(2)
  ML2a = - DTerm  - Abs( bi(2) )**2 + Real( sumC,dp ) / vevL

  Call PseudoScalarMassEps1nT(MH1sq,MH2sq,ML2a,bi,B,vevSM,vevL,gp,g, &
                            & mP0,mP02,RP0,kont)

  Iname = Iname - 1

 End Subroutine PseudoScalarMassEps1

 Subroutine PseudoScalarMassEps1nT(MH1sq,MH2sq,ML2,bi,B,vevSM,vevL,gp,g, &
                                 & mP0,mP02,RP0,kont,NoSymmetryBreaking)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the 1-generation epsilon model
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  ML2 .......... M_L33 squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL ......... v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, ML2, vevSM(2), vevL, g, gp
  Real(Dp), Intent(out) :: mP0(3), mP02(3), RP0(3,3)
  Complex(Dp), Intent(in) ::  bi(2), B(2)
  Integer, Intent(inout) :: kont
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Real(Dp) :: test(2), mat3(3,3), DTerm
  Integer :: i1, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassEps1nT'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2 + vevL**2 ) 

  mat3(1,1) = MH1sq + DTerm + Abs(bi(1))**2
  mat3(1,2) = Real( B(1),dp )
  mat3(1,3) = - Real( Conjg(bi(1)) * bi(2),dp )
  mat3(2,1) = mat3(1,2)
  mat3(2,2) = MH2sq - DTerm + Abs(bi(1))**2
  mat3(2,3) = - Real( B(2),dp )
  mat3(3,1) = mat3(1,3)
  mat3(3,2) = mat3(2,3)
  mat3(3,3) = ML2 + DTerm + Abs( bi(2) )**2

  Call EigenSystem(mat3,mP02,RP0,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) 'gp, g', g, gp
   Write(ErrCan,*) 'bi ',bi 
   Write(ErrCan,*) 'B ',B
   Write(ErrCan,*) 'ML2 ',ML2
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,3
    If (mP02(i1).Ge.0._dp) Then
     mP0(i1) = Sqrt( mP02(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mP02(i1) 
      Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq
      Write(ErrCan,*) 'ML2 ',ML2 
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'gp, g', g, gp
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'B ',B
     End If
     mP0(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -206
     Call AddError(206)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mP0(1) = mZ
   mP02(1) = mZ**2
   Do i1=2,3
    If (mP02(i1).Ge.0._dp) Then
     mP0(i1) = Sqrt( mP02(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mP02(i1) 
      Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'gp, g', g, gp
      Write(ErrCan,*) 'bi ',bi 
      Write(ErrCan,*) 'B ',B
      Write(ErrCan,*) 'ML2 ',ML2
     End If
     mP0(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -206
     Call AddError(206)
    End If
   End Do
  End If
  Iname = Iname - 1

 End Subroutine PseudoScalarMassEps1nT


 Subroutine PseudoScalarMassEps3(ML2,bi,B,vevSM,vevL,gp,g, mP0,mP02,RP0,kont)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the 3-generation epsilon model
 ! Input:
 !  ML2(i,j) ..... M_Lij squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp
  Real(Dp), Intent(out) :: mP0(5), mP02(5), RP0(5,5)
  Complex(Dp), Intent(in) :: ML2(3,3), bi(4), B(4)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm
  Complex(Dp) :: ML2a(3,3), sumC
  Integer :: i1,i2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassEps3'

  DTerm = 0.125_dp * (g**2+gp**2)  &
      & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * Dot_product(bi(2:4),vevL),dp ) ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1)              &
      &   - Dot_product( Real(B(2:4),dp) , vevL) ) / vevSM(2)
   
  ML2a = ML2
  Do i1=1,3
   sumC = - vevSM(2) * B(i1+1) + vevSM(1) * Conjg(bi(1)) * bi(i1+1)
   Do i2=1,3
    sumC = sumC - Conjg( bi(i1+1) ) * bi(i2+1) * vevL(i2)
    If (i1.Ne.i2) sumC = sumC - vevL(i2) * ML2a(i2,i1)
   End Do
   ML2a(i1,i1) = - DTerm + Real( sumC,dp ) / vevL(i1) 
  End Do

  Call PseudoScalarMassEps3nT(MH1sq,MH2sq,ML2a,bi,B,vevSM,vevL,gp,g, &
                            & mP0,mP02,RP0,kont)

  Iname = Iname - 1

 End Subroutine PseudoScalarMassEps3

 Subroutine PseudoScalarMassEps3nT(MH1sq,MH2sq,ML2,bi,B,vevSM,vevL,gp,g, &
                                 & mP0,mP02,RP0,kont,NoSymmetryBreaking)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the 3-generation epsilon model
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  ML2(i,j) ..... M_Lij squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, vevSM(2), vevL(3), g, gp
  Real(Dp), Intent(out) :: mP0(5), mP02(5), RP0(5,5)
  Complex(Dp), Intent(in) :: ML2(3,3), bi(4), B(4)
  Integer, Intent(inout) :: kont
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Real(Dp) :: mat5(5,5), DTerm, test(2)
  Integer :: i1,i2, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassEps3nT'

  DTerm = 0.125_dp * (g**2+gp**2)  &
      & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  mat5(1,1) = MH1sq + DTerm + Abs(bi(1))**2
  mat5(1,2) = Real( B(1),dp )
  mat5(1,3:5) = - Real( Conjg(bi(1)) * bi(2:4),dp )
  mat5(2,2) = MH2sq - DTerm + Dot_product(bi,bi)
  mat5(2,3:5) = - Real( B(2:4),dp )
  mat5(3,3) = Real(ML2(1,1),dp) + DTerm + Abs( bi(2) )**2
  mat5(3,4:5) = Real( Conjg(bi(2)) * bi(3:4) + ML2(1,2:3),dp )
  mat5(4,4) = Real(ML2(2,2),dp) + DTerm + Abs( bi(3) )**2
  mat5(4,5) = Real( Conjg(bi(3)) * bi(4) + ML2(2,3),dp )
  mat5(5,5) = Real(ML2(3,3),dp) + DTerm + Abs( bi(4) )**2

  Do i1=2,5
   Do i2=1,i1-1
    mat5(i1,i2)=mat5(i2,i1)
   End Do
  End Do

  Call EigenSystem(mat5,mP02,RP0,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) 'gp, g', g, gp
   Write(ErrCan,*) 'bi ',bi 
   Write(ErrCan,*) 'B ',B
   Write(ErrCan,*) 'ML2 ',ML2
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,5
    If (mP02(i1).Ge.0._dp) Then
     mP0(i1) = Sqrt( mP02(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mP02(i1) 
      Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq
      Write(ErrCan,*) 'ML2 ',ML2 
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'gp, g', g, gp
      Write(ErrCan,*) 'bi ',bi
      Write(ErrCan,*) 'B ',B
     End If
     mP0(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -207
     Call AddError(207)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mP0(1) = mZ
   mP02(1) = mZ**2
   Do i1=2,5
    If (mP02(i1).Ge.0._dp) Then
     mP0(i1) = Sqrt( mP02(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mP02(i1) 
      Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'vevL ',vevL
      Write(ErrCan,*) 'gp, g', g, gp
      Write(ErrCan,*) 'bi ',bi 
      Write(ErrCan,*) 'B ',B
      Write(ErrCan,*) 'ML2 ',ML2
     End If
     mP0(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -207
     Call AddError(207)
    End If
   End Do
  End If

  Iname = Iname - 1

 End Subroutine PseudoScalarMassEps3nT


 Subroutine PseudoScalarMassMSSM(mu,B,vevSM,gp,g,mP0,mP02,RP0,kont)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the MSSM
 ! Input:
 !  mu ........... i=1 the mu parameter
 !  B ............ bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), g, gp
  Real(Dp), Intent(out) :: mP0(2), mP02(2), RP0(2,2)
  Complex(Dp), Intent(in) ::  mu, B
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassMSSM'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2) 

  MH1sq = - Abs( mu )**2 - DTerm + Real( B,dp ) * vevSM(2) / vevSM(1)
  MH2sq = - Abs( mu )**2 + DTerm + Real( B,dp ) * vevSM(1) / vevSM(2)
   
  Call PseudoScalarMassMSSMnT(MH1sq,MH2sq,mu,B,vevSM,gp,g,mP0,mP02,RP0,kont)

  Iname = Iname - 1

 End Subroutine PseudoScalarMassMSSM

 Subroutine PseudoScalarMassMSSMnT(MH1sq,MH2sq,mu,B,vevSM,gp,g, &
                                 & mP0,mP02,RP0,kont,NoSymmetryBreaking)
 !-----------------------------------------------------------------------
 ! calculates the pseudoscalar masses in the MSSM
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  mu ........... the mu parameter
 !  B ............ bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 !  NoSymmetryBreaking ... optional parameter, if present, it is not assumed
 !                         that the lightest state is a Goldstone boson
 ! OutPut:
 !  mP0(i) ....... the masses
 !  mP02(i) ...... the masses squared
 !  RP0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 11.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, vevSM(2), g, gp
  Real(Dp), Intent(out) :: mP0(2), mP02(2), RP0(2,2)
  Complex(Dp), Intent(in) ::  Mu, B 
  Integer, Intent(inout) :: kont
  Integer, Optional, Intent(in) :: NoSymmetryBreaking

  Real(Dp) :: mat2(2,2), DTerm, test(2)
  Integer :: i1, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoScalarMassMSSMnT'

  kont = 0
  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2 )

  mat2(1,1) = MH1sq + DTerm + Abs(mu)**2
  mat2(1,2) = Real( B,dp )
  mat2(2,2) = MH2sq - DTerm + Abs(mu)**2
  mat2(2,1) = mat2(1,2)

  Call EigenSystem(mat2,mP02,RP0,ierr, test)
  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'gp, g', g, gp
   Write(ErrCan,*) 'mu ',mu
   Write(ErrCan,*) 'B ',B
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  If (Present(NoSymmetryBreaking)) Then
   Do i1=1,2
    If (mP02(i1).Ge.0._dp) Then
     mP0(i1) = Sqrt( mP02(i1) )
    Else
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
      Write(ErrCan,*) 'a mass squared is negative :',i1,mP02(i1) 
      Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
      Write(ErrCan,*) 'vevSM ',vevSM
      Write(ErrCan,*) 'gp, g', g, gp
      Write(ErrCan,*) 'mu ',mu
      Write(ErrCan,*) 'B ',B
     End If
     mP0(i1) = 0._dp
     If (ErrorLevel.Eq.2) Call TerminateProgram
     kont = -208
     Call AddError(208)
    End If
   End Do

  Else
   !--------------------
   ! xsi=1 gauge
   !--------------------
   mP0(1) = mZ
   mP02(1) = mZ**2
   If (mP02(2).Ge.0._dp) Then
    mP0(2) = Sqrt( mP02(2) )
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) 'a mass squared is negative :',2,mP02(2) 
     Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
     Write(ErrCan,*) 'vevSM ',vevSM
     Write(ErrCan,*) 'gp, g', g, gp
     Write(ErrCan,*) 'mu ',mu
     Write(ErrCan,*) 'B ',B
    End If
    mP0(2) = 0._dp
    If (ErrorLevel.Eq.2) Call TerminateProgram
    kont = -208
    Call AddError(208)
   End If
  End If

  Iname = Iname - 1

 End Subroutine PseudoScalarMassMSSMnT


 Subroutine ScalarMassEps1(bi,B,vevSM,vevL,gp,g, mStop2,mT2, mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the 1-generation epsilon model
 ! Input:
 !                 i=2 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL ......... v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL, g, gp, mStop2(2), mT2
  Real(Dp), Intent(out) :: mS0(3), mS02(3), RS0(3,3)
  Complex(Dp), Intent(in) :: bi(2), B(2)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, ML2a, DTerm
  Complex(Dp) ::  sumC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassEps1'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2 + vevL**2)

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * bi(2),dp) * vevL ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1) + Real(B(2),dp) * vevL) / vevSM(2)
   
  sumC = - vevSM(2) * B(2) + vevSM(1) * Conjg(bi(1)) * bi(2)
  ML2a = - DTerm  - Abs( bi(2) )**2 + Real( sumC,dp ) / vevL

  Call ScalarMassEps1nT(MH1sq,MH2sq,ML2a,bi,B,vevSM,vevL,gp,g, mStop2, mT2, &
                      & mS0,mS02,RS0,kont)

  Iname = Iname - 1

 End Subroutine ScalarMassEps1

 Subroutine ScalarMassEps1nT(MH1sq,MH2sq,ML2,bi,B,vevSM,vevL,gp,g,mStop2,mT2, &
                            & mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the 1-generation epsilon model
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  ML2 .......... M_L33 squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL ......... v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, ML2, vevSM(2), vevL, g, gp, mStop2(2) &
                         & , mT2
  Real(Dp), Intent(out) :: mS0(3), mS02(3), RS0(3,3)
  Complex(Dp), Intent(in) :: bi(2), B(2)
  Integer, Intent(inout) :: kont

  Real(Dp) :: mat3(3,3), DTerm, g2, c2, test(2)
  Integer :: i1, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassEps1nT'

  g2 = 0.25 * (g**2+gp**2)
  DTerm = 0.5_dp * g2 * ( vevSM(1)**2 - vevSM(2)**2 + vevL**2 ) 

  mat3(1,1) = MH1sq + DTerm + Abs(bi(1))**2 + g2 * vevSM(1)**2
  mat3(1,2) = - Real( B(1),dp ) - g2 * vevSM(1) * vevSM(2)
  mat3(1,3) = - Real( Conjg(bi(1)) * bi(2),dp ) + g2 * vevSM(1) * vevL
  mat3(2,1) = mat3(1,2)
  mat3(2,2) = MH2sq - DTerm + Abs(bi(1))**2 + g2 * vevSM(2)**2
  mat3(2,3) = Real( B(2),dp ) -  g2 * vevSM(2) * vevL
  mat3(3,1) = mat3(1,3)
  mat3(3,2) = mat3(2,3)
  mat3(3,3) = ML2 + DTerm + Abs( bi(2) )**2 + g2 * vevL**2

  c2 = 1._dp - vevL**2 / (vevSM(1)**2 + vevSM(2)**2 + vevL**2)
  mat3(2,2) = mat3(2,2) + 3._dp * mT2**2 *Log(mStop2(1)*mStop2(2)/mT2**2) &
            &             / (4._dp * Pi2 * vevSM(2)**2 * c2)

  Call EigenSystem(mat3,mS02,RS0,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write (ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write (ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write (ErrCan,*) 'vevSM ',vevSM
   Write (ErrCan,*) 'vevL ',vevL
   Write (ErrCan,*) 'gp, g', g, gp
   Write (ErrCan,*) 'bi ',bi 
   Write (ErrCan,*) 'B ',B
   Write (ErrCan,*) 'ML2 ',ML2
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  Do i1=1,3
   If (mS02(i1).Ge.0._dp) Then
    mS0(i1) = Sqrt( mS02(i1) )
   Else
    If (ErrorLevel.Ge.0) Then
     Write (ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) 'a mass squared is negative :',i1,mS02(i1) 
     Write (ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq
     Write (ErrCan,*) 'ML2 ',ML2 
     Write (ErrCan,*) 'vevSM ',vevSM
     Write (ErrCan,*) 'vevL ',vevL
     Write (ErrCan,*) 'gp, g', g, gp
     Write (ErrCan,*) 'bi ',bi
     Write (ErrCan,*) 'B ',B
    End If
    mS0(i1) = 0._dp
    If (ErrorLevel.Eq.2) Call TerminateProgram
    kont = -210
    Call AddError(210)
   End If
  End Do

  Iname = Iname - 1

 End Subroutine ScalarMassEps1nT


 Subroutine ScalarMassEps3(ML2,bi,B,vevSM,vevL,gp,g, mStop2, mTop2 &
                         &, mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the 3-generation epsilon model
 ! Input:
 !  ML2(i,j) ..... M_Lij squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), mStop2(2), mTop2, g, gp
  Real(Dp), Intent(out) :: mS0(5), mS02(5), RS0(5,5)
  Complex(Dp), Intent(in) :: ML2(3,3), bi(4), B(4)
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm
  Complex(Dp) :: ML2a(3,3), sumC
  Integer :: i1,i2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassEps3'

  DTerm = 0.125_dp * (g**2+gp**2)  &
      & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  MH1sq = - Abs( bi(1) )**2 - DTerm  &
      & + ( Real( B(1),dp ) * vevSM(2)  &
      &   + Real( Conjg(bi(1)) * Dot_product(bi(2:4),vevL),dp ) ) / vevSM(1)
  MH2sq = - Dot_product( Conjg(bi), bi) + DTerm  &
      & + ( Real( B(1),dp ) * vevSM(1)              &
      &   - Dot_product( Real(B(2:4),dp) , vevL) ) / vevSM(2)
   
  ML2a = ML2
  Do i1=1,3
   sumC = - vevSM(2) * B(i1+1) + vevSM(1) * Conjg(bi(1)) * bi(i1+1)
   Do i2=1,3
    sumC = sumC - Conjg( bi(i1+1) ) * bi(i2+1) * vevL(i2)
    If (i1.Ne.i2) sumC = sumC - vevL(i2) * ML2a(i2,i1)
   End Do
   ML2a(i1,i1) = - DTerm + Real( sumC,dp ) / vevL(i1) 
  End Do

  Call ScalarMassEps3nT(MH1sq,MH2sq,ML2a,bi,B,vevSM,vevL,gp,g, &
                      & mStop2,mTop2,mS0,mS02,RS0,kont)

  Iname = Iname - 1

 End Subroutine ScalarMassEps3

 Subroutine ScalarMassEps3nT(MH1sq,MH2sq,ML2,bi,B,vevSM,vevL,gp,g, &
                            & mStop2,mTop2,mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the 3-generation epsilon model
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  ML2(i,j) ..... M_Lij squared
 !  bi(i) ........ i=1 the mu parameter
 !                 i=2 epsilon_1
 !                 i=3 epsilon_2
 !                 i=4 epsilon_3
 !  B(i) ......... bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared  i=1 ... B_0
 !                               i=2 ... B_1
 !                               i=3 ... B_2
 !                               i=4 ... B_3
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  vevL(i) ...... i=1 v_1
 !                 i=2 v_2
 !                 i=3 v_3
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, vevSM(2), vevL(3), mStop2(2), mTop2 &
                    & , g, gp
  Real(Dp), Intent(out) :: mS0(5), mS02(5), RS0(5,5)
  Complex(Dp), Intent(in) :: ML2(3,3), bi(4), B(4)
  Integer, Intent(inout) :: kont

  Real(Dp) :: mat5(5,5), DTerm, g2, c2, test(2)
  Integer :: i1,i2, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassEps3nT'

  g2 = 0.25 * (g**2+gp**2)
  DTerm = 0.5_dp * g2 * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL)) 

  mat5(1,1) = MH1sq + DTerm + Abs(bi(1))**2 + g2 * vevSM(1)**2
  mat5(1,2) = - Real( B(1),dp ) - g2 * vevSM(1) * vevSM(2)
  mat5(1,3:5) = - Real( Conjg(bi(1)) * bi(2:4),dp ) + g2 * vevSM(1) * vevL
  mat5(2,2) = MH2sq - DTerm + Abs(bi(1))**2 + g2 * vevSM(2)**2
  mat5(2,3:5) = Real( B(2:4),dp ) -  g2 * vevSM(2) * vevL
  mat5(3,3) = Real(ML2(1,1),dp) + DTerm + Abs( bi(2) )**2 + g2 * vevL(1)**2
  mat5(3,4:5) = Real( Conjg(bi(2)) * bi(3:4) + ML2(1,2:3),dp ) &
            &  + g2 * vevL(1) * vevL(2:3)
  mat5(4,4) = Real(ML2(2,2),dp) + DTerm + Abs( bi(3) )**2 + g2 * vevL(2)**2
  mat5(4,5) = Real( Conjg(bi(3)) * bi(4) + ML2(2,3),dp ) + g2 *vevL(2) *vevL(3)
  mat5(5,5) = Real(ML2(3,3),dp) + DTerm + Abs( bi(4) )**2 + g2 * vevL(3)**2

  c2 = 1._dp - (vevL(1)**2 + vevL(2)**2 + vevL(3)**2)        &
     &     / (vevSM(1)**2 + vevSM(2)**2 + vevL(1)**2 + vevL(2)**2 +vevL(3)**2)
  mat5(2,2) = mat5(2,2) + 3._dp * mTop2**2 *Log(mStop2(1)*mStop2(2)/mTop2**2) &
            &             / (4._dp * Pi2 * vevSM(2)**2 * c2)
 
  Do i1=2,5
   Do i2=1,i1-1
    mat5(i1,i2)=mat5(i2,i1)
   End Do
  End Do

  Call EigenSystem(mat5,mS02,RS0,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'vevL ',vevL
   Write(ErrCan,*) 'gp, g', g, gp
   Write(ErrCan,*) 'bi ',bi 
   Write(ErrCan,*) 'B ',B
   Write(ErrCan,*) 'ML2 ',ML2
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  Do i1=1,5
   If (mS02(i1).Ge.0._dp) Then
    mS0(i1) = Sqrt( mS02(i1) )
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) 'a mass squared is negative :',i1,mS02(i1) 
     Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq
     Write(ErrCan,*) 'ML2 ',ML2 
     Write(ErrCan,*) 'vevSM ',vevSM
     Write(ErrCan,*) 'vevL ',vevL
     Write(ErrCan,*) 'gp, g', g, gp
     Write(ErrCan,*) 'bi ',bi
     Write(ErrCan,*) 'B ',B
    End If
    mS0(i1) = 0._dp
    If (ErrorLevel.Eq.2) Call TerminateProgram
    kont = -211
    Call AddError(211)
   End If
  End Do

  Iname = Iname - 1

 End Subroutine ScalarMassEps3nT

 Subroutine ScalarMassMSSMeff(mA2, mZ2, tanb, vevSM                  &
                            &, mT, mStop2, Atop, mu, mB, mSbot2, Abot   &
                            &, mS0, mS02, RS0, kont, LoopHiggs)
 !-----------------------------------------------------------------------
 ! calculates the scalar  Higgs masses in the MSSM
 ! using the effective potential method
 ! Ellis, Ridolfi, Zwirner, Phys. Lett. B262 (1991) 477
 ! Input:
 !  mu ........... i=1 the mu parameter
 !  B ............ bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !   note that the dimension is not fixed, so that I can latter on 
 !   include CP violation
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: mA2, mZ2, tanb, vevSM(2)
  Real(Dp), Intent(in) :: mT, mStop2(2), mB, mSbot2(2)
  Complex(Dp), Intent(in) ::  mu, Atop, Abot
  Logical, Optional :: LoopHiggs
  Real(Dp), Intent(out) :: mS0(2), mS02(2), RS0(2,2)
  Integer, Intent(inout) :: kont

  Integer :: ierr
  Logical :: LoopHiggs_l
  Real(Dp) :: dmst, dmsb, sin2b, sinb2, cosb2, higdel
  Real(Dp) :: fakt1, fakt2, abmutb, atmucb, delmat(2,2), mat(2,2)
  Real(Dp) :: amtop4, ambot4, E2(2), S2(2,2), test(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassMSSMeff'

  if (Present(LoopHiggs)) then
    LoopHiggs_l = LoopHiggs
  else
    LoopHiggs_l = .True.
  end if

  dmst = mstop2(1) - mstop2(2)
  dmsb = msbot2(1) - msbot2(2)

  cosb2 = 1._dp / (1._dp + tanb**2)
  sin2b = 2._dp * tanb * cosb2
  sinb2 = tanb**2 * cosb2

  higdel = mA2 * sin2b

  amtop4 = mT**4
  ambot4 = mB**4
  fakt1 = amtop4 / sinb2
  fakt2 = ambot4 / cosb2
  abmutb = abot - mu * tanb
  atmucb = atop - mu / tanb

  delmat(1,1) = fakt2 * ( Log(msbot2(1)*msbot2(2)/ambot4)                     &
     &              + 2._dp * abot * abmutb * Log(msbot2(1)/msbot2(2)) / dmsb &
     &              + (abot*abmutb/dmsb)**2 * gg(msbot2(1),msbot2(2)) )       &
     &              + fakt1 * (mu*atmucb/dmst)**2 * gg(mstop2(1),mstop2(2))
  delmat(2,2) = fakt1 * (Log(mstop2(1)*mstop2(2)/amtop4)                      &
     &              + 2._dp * atop * atmucb * Log(mstop2(1)/mstop2(2)) / dmst &
     &              + (atop*atmucb/dmst)**2 * gg(mstop2(1),mstop2(2)) )       &
     &              + fakt2 * (mu*abmutb/dmsb)**2 * gg(msbot2(1),msbot2(2))
  delmat(1,2) = - fakt1 * mu * atmucb / dmst                                  &
     &            * ( Log(mstop2(1)/mstop2(2))                                &
     &              + atop*atmucb/dmst * gg(mstop2(1),mstop2(2)) )            &
     &          - fakt2 * mu * abmutb / dmsb                                  &
     &            * ( Log(msbot2(1)/msbot2(2))                                &
     &              + abot*abmutb/dmsb * gg(msbot2(1),msbot2(2)) )
  delmat(2,1) = delmat(1,2)

  fakt1 = 1.5_dp / (pi2 * (vevSM(1)**2 + vevSM(2)**2) )

  mat(1,1) = mZ2 * 2._dp * cosb2 + tanb * higdel
  mat(1,2) = -mZ2 * sin2b - higdel
  mat(2,1) = mat(1,2)
  mat(2,2) = mZ2 * 2._dp * sinb2 + higdel / tanb

  if (LoopHiggs_l) then
   mat = 0.5_dp * ( mat + fakt1 * delmat)
  else
   mat = 0.5_dp * mat
  end if

  Call RealEigenSystem(Mat,E2,S2,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Eq.0).And.(E2(1).Ge.0._dp)) Then
   mS02 = E2
   mS0 = Sqrt(E2)
   RS0 = S2
  Else If ((ierr.Eq.0).And.(E2(1).Le.0._dp)) Then
   kont = -212
   Call AddError(212)
   If (ErrorLevel.Ge.0) Then
    Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname)
    Write(ErrCan,*) 'm_h^2 = ',E2(1)
    Write(ErrCan,*) 'Setting m_h  to the sqrt(abs(m^2_h))'
    If (ErrorLevel.Eq.2) Call TerminateProgram
   End If
   mS02 = E2
   mS0 = Sqrt(Abs(E2))
   RS0 = S2
  Else 
   kont = ierr
   If (ErrorLevel.Ge.0) Then
    Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname)
    Write(ErrCan,*) 'diagonalization failed'
    If (ErrorLevel.Eq.2) Call TerminateProgram
   End If
  End If

  !----------------------------------
  ! phase convention
  !----------------------------------
  If (RS0(1,2).Lt.0._dp) RS0 = - RS0
  Iname = Iname - 1

 Contains
  Real(dp) Function gg(x,y)
  Implicit None
   Real(dp), Intent(in) :: x, y
    If (x.Eq.y) Then
     gg = 0._dp
    Else If (Abs(x-y)/x.Lt.1.e-8_dp) Then
     gg = 2._dp * (1._dp - y)
    Else
     gg = 2._dp - (x+y)/(x-y)*Log(x/y)
    End If
  End Function gg

 End Subroutine ScalarMassMSSMeff


 Subroutine ScalarMassMSSM(mu,B,vevSM,gp,g,mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the MSSM
 ! Input:
 !  mu ........... i=1 the mu parameter
 !  B ............ bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), g, gp
  Real(Dp), Intent(out) :: mS0(2), mS02(2), RS0(2,2)
  Complex(Dp), Intent(in) ::  mu, B
  Integer, Intent(inout) :: kont

  Real(Dp) :: MH1sq, MH2sq, DTerm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassMSSM'

  DTerm = 0.125_dp * (g**2+gp**2) * ( vevSM(1)**2 - vevSM(2)**2) 

  MH1sq = - Abs( mu )**2 - DTerm + Real( B,dp ) * vevSM(2) / vevSM(1)
  MH2sq = - Abs( mu )**2 + DTerm + Real( B,dp ) * vevSM(1) / vevSM(2)
   
  Call ScalarMassMSSMnT(MH1sq,MH2sq,mu,B,vevSM,gp,g,mS0,mS02,RS0,kont)

  Iname = Iname - 1

 End Subroutine ScalarMassMSSM

 Subroutine ScalarMassMSSMnT(MH1sq,MH2sq,mu,B,vevSM,gp,g,mS0,mS02,RS0,kont)
 !-----------------------------------------------------------------------
 ! calculates the scalar masses in the MSSM 
 ! without referring to the tad-pole equations.
 ! Input:
 !  MH1sq ........ M_H_1 squared
 !  MH2sq ........ M_H_2 squared
 !  mu ........... the mu parameter
 !  B ............ bilinear soft SUSY breaking Parameter with dimension 
 !                 mass squared 
 !  vevSM(i) ..... i=1 v_D
 !                 i=2 v_U
 !  gp ........... U(1)_Y gauge coupling g'
 !  g ............ SU(2)_L gauge coupling g
 ! OutPut:
 !  mS0(i) ....... the masses
 !  mS02(i) ...... the masses squared
 !  RS0(i,J) ..... mixing matrix
 !  kont ......... control parameter, =0 if everything is o.k.
 ! written by Werner Porod: 16.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: MH1sq, MH2sq, vevSM(2), g, gp
  Real(Dp), Intent(out) :: mS0(2), mS02(2), RS0(2,2)
  Complex(Dp), Intent(in) :: mu, B
  Integer, Intent(inout) :: kont

  Real(Dp) :: mat2(2,2), DTerm, g2, test(2)
  Integer :: i1, ierr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarMassMSSMnT'

  g2 = 0.25 * (g**2+gp**2)
  DTerm = 0.5_dp * g2 * ( vevSM(1)**2 - vevSM(2)**2 ) 

  mat2(1,1) = MH1sq + DTerm + Abs(mu)**2 + g2 * vevSM(1)**2
  mat2(1,2) = - Real( B,dp ) - g2 * vevSM(1) * vevSM(2)
  mat2(2,1) = mat2(1,2)
  mat2(2,2) = MH1sq - DTerm + Abs(mu)**2 + g2 * vevSM(2)**2

  Call EigenSystem(mat2,mS02,RS0,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.-1)) Then
   Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Diagonalization failed, ierr :',ierr
   Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq 
   Write(ErrCan,*) 'vevSM ',vevSM
   Write(ErrCan,*) 'gp, g', g, gp
   Write(ErrCan,*) 'mu ',mu
   Write(ErrCan,*) 'B ',B
   kont = ierr
   Iname = Iname - 1
   Return
  End If

  Do i1=1,2
   If (mS02(i1).Ge.0._dp) Then
    mS0(i1) = Sqrt( mS02(i1) )
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from Subroutine '//NameOfUnit(Iname)
     Write(ErrCan,*) 'a mass squared is negative :',i1,mS02(i1) 
     Write(ErrCan,*) 'MH1sq, MH2sq ', MH1sq, MH2sq
     Write(ErrCan,*) 'vevSM ',vevSM
     Write(ErrCan,*) 'gp, g', g, gp
     Write(ErrCan,*) 'mu ',mu
     Write(ErrCan,*) 'B ',B
    End If
    mS0(i1) = 0._dp
    If (ErrorLevel.Eq.2) Call TerminateProgram
    kont = -213
    Call AddError(213)
   End If
  End Do

  Iname = Iname - 1

 End Subroutine ScalarMassMSSMnT

 Subroutine ScalarMassNMSSMeff(h0, lam, mu_in, vevSM, vP &
                & , Al_in, Ak, mB, mT, mStop2, RStop, mSbot2, RSbot        &
                & , M2_L1, M2_E1, M2_Q1, M2_U1, M2_D1, M2_E, M2_L &
                & , M2_Q, M2_U, M2_D, At, Ab, M2  &
                & , mS0, mS02, RS0, mP0, mP02, RP0, mSpm, mSpm2, RSpm, kont)
 !-----------------------------------------------------------------------------
 !
 !	This subroutine computes the Higgs masses and couplings in
 !	the NMSSM. The relevant input parameters are
 !
 !	PAR(1) = lambda
 !	PAR(2) = kappa
 !	PAR(3) = tan(beta)	at the scale M_top
 !	PAR(4) = mu		effective mu term (=lambda*s)
 !	PAR(5) = A_lambda
 !	PAR(6) = A_kappa
 !	PAR(7) = MQ3           Left-handed Squarks
 !	PAR(8) = MU3           Right-handed Stops
 !	PAR(9) = MD3           Right-handed Sbottoms
 !	PAR(12) = A_TOP
 !	PAR(13) = A_BOT         
 !	PAR(20) = M1
 !	PAR(21) = M2  
 !
 !	On output: 
 !
 !	mS0(1-3): CP-even masses (ordered)
 !
 !	RS0(1-3,1-3): Mixing angles: if HB(I) are the bare states,
 !	  HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates, 
 !	  the convention is HB(I) = SUM_(J=1,3) RS0(J,I)*HM(J)
 !	  which is equivalent to HM(I) = SUM_(J=1,3) RS0(I,J)*HB(J)
 !
 !	mP0(1-2): CP-odd masses (ordered)
 !
 !	RP0(1-2,1-2): Mixing angles: if AB(I) are the bare states,
 !	  AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates, 
 !	  the convention is 
 !	  AM(I) = RP0(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
 !			+ RP0(I,2)*AB(3)
 !
 !	mSpm: Charged Higgs mass
 !
 !	kont:  =   0	        OK
 !		=   1,3,5,7	mS0(1)**2 < 0
 !		=   2,3,6,7	mP0(1)**2 < 0
 !		=   4,5,6,7	mSpm**2 < 0
 !		=   8		MST1 or MSB1 <0
 !		=   10   	Singular parameters (l, k, mu or tanb = 0)
 !
 !	The precision in the computation of the lightest Higgs mass is:
 !
 !	Terms ~ ht^4/hb^4 from (s)top/(s)bottom-loops are computed
 !	exactly in the mixing parameters.
 !
 !	Terms ~ g^2*(ht^2/hb^2) (where g is an electro-weak gauge coupling)
 !	are taken into account due to the wave function renormalizations
 !	factors, finite self energies (pole masses) and corrections from
 !       stop/sbottom D terms
 !
 !	Leading logs ~ (g, l, k)^4 are added explicitely.
 !
 !	Leading double logs from two loops ~ ht/hb^6 and ht/hb^4*alpha_s 
 !	are taken into account.
 !
 !	For heavy higgses (with masses mhh > mtop) the leading log
 !	contributions ~ (ht^2/hb^2)*log(mhh/mtop) to the pole masses
 !	are taken into account, but not the corresponding effects on the
 !	mixing angles.
 !
 !	All mixing angles are at the scale m_top.
 !
 !	The dominant errors come from one loop terms ~ (g,l,k)^4 without
 !	large logs, and from two loop terms without large double logs.
 !
 ! things to be done:
 !   - replace the fixed value of alpha_s(m_b) by a calculation
 !-----------------------------------------------------------------------------
 Implicit None
  Integer, intent(inout) :: kont
  Real(dp), Intent(in) :: h0, lam
  Real(dp), Intent(in) :: mB, mT, M2_Q, M2_U, M2_D, mStop2(2), mSbot2(2)
  Real(dp), Intent(in) :: M2_Q1, M2_U1, M2_D1, M2_E1, M2_L1, M2_E, M2_L
  Complex(dp), Intent(in) :: Rsbot(2,2), RStop(2,2)
  Real(dp), intent(in) :: mu_in, vevSM(2), vP
  Real(dp), Intent(in) :: Ak, Al_in, M2, At, Ab 

  Real(dp), Intent(out) :: mS0(3), mS02(3), RS0(3,3), mP0(3), mP02(3), RP0(3,3) &
     & , mSpm(2), mSpm2(2)
  Complex(dp), Intent(out) :: RSpm(2,2)

  
  Integer :: I
  Real(dp) :: k, mT2, L, tanb, tanbq, Zht, Zhb, QSTSB, cb, sb, s2, t, g
  Real(dp) :: nu, B, c2tw, s2tw, h1q, h2q, htq, hbq, Al, h1, h2
  Real(dp) :: rt, rb, ALSMT, MA2, MS2, MP2, M12, ALSQ, mtopq, mbotq, s2t, s2b
  Real(dp) :: emt, fmt, gmt, Xt, ht
  Real(dp) :: emb, fmb, gmb, Xb, hb, asmb
  Real(dp) :: Alshift, gb, gt, Db, Dt, D1, D2,X,BT,BB
  Real(dp) :: mdia(3),moff(3),rdia(3),roff(3), mat3(3,3), test(2)
  Real(dp) :: mp(3),subdet
  Real(dp) :: sferm,bos, g1, g2
  Real(dp) :: LM2,Lmu,Lnu,LM2mu,Lmunu,LA,LP,LS,LPP,P1,P2
  real(dp) :: mu ! later complex

  kont = 0

  k = 0.5_dp * lam ! Cyrill Hugonie et al. use a different convention

  mT2 = mT**2
  tanb = vevSM(2) / vevSM(1)
  mu = mu_in + oosqrt2 * h0 * vP
 !-----------------------------------------------------------
 ! these parts corresponds to the ones of MSFER in NMHDECAY
 !-----------------------------------------------------------
  QSTSB = Max( (2._dp*M2_Q + M2_U + M2_D)/4._dp, mT2)
  t = LOG(QSTSB/MT**2)

 !   Alphas at MT and QSTSB (for 2 loop corrs.)
  ALSMT = Alphas_mZ / (1._dp + oo4pi *23._dp /3._dp * Alphas_mZ*Log(MT2/MZ2))
  ALSQ = ALSMT / (1._dp + 7._dp * oo4pi * ALSMT * T)

  h1 = tanb * Sqrt(1._dp / (2._dp * sqrt2 * (1._dp+tanb**2)*G_F) )
  ht = MT/(1._dp+4._dp*ALSMT/(3._dp*PI)+11._dp*(ALSMT/PI)**2)/h1

  h2 = Sqrt(1._dp / (2._dp * sqrt2 * (1._dp+tanb**2)*G_F) )
  asmb = AlphaS_mb
  hb = mf_d(3) * (ALsmt/AsMB)**(12._dp/23._dp)  &
     &         *(1._dp +7462._dp*(Alsmt-AsMB)/ (4._dp*PI*1587._dp)) / h2

  htq = ht * (1._dp + oo4pi * 7._dp * ALSMT * t)**(-4._dp/7._dp)  &
      &    * (1._dp + oo64pi2 * (9._dp*ht**2+hb**2) * t)
  hbq = hb * (ALSQ/ALsMT)**(4._dp/7._dp)                            &
      &    * (1._dp + 7398._dp * (ALSQ-ALsMT)/ (4._dp*PI*1323._dp)) &
      &    * (1._dp + 0*oo64pi2 * (9._dp*hb**2+ht**2) * t)

 !   Approximate value for the MSSM-like CP odd Higgs mass squared MA^2:
  MA2=(Al_in+k/h0*MU)*MU*(tanb+1._dp/tanb)
  IF (MA2.GT.MT**2) THEN
    htq=htq*(1._dp + oo32pi2 * hb**2 * Log(Min(MA2,QSTSB)/MT2) )
    hbq=hbq*(1._dp + 0*oo32pi2 * ht**2 * Log(Min(MA2,QSTSB)/MT2) )
  END IF

  Zht = 1._dp + 3._dp * oo16pi2 * htq**2 * t
  Zhb = 1._dp + 3._dp * oo16pi2 * hbq**2 * t

  tanbq = tanb * Sqrt(Zhb/Zht)

  h1q = tanb * Sqrt(1._dp / (2._dp * sqrt2 * (1._dp+tanb**2)*G_F*Zht) )
  h2q = Sqrt(1._dp / (2._dp * sqrt2 * (1._dp+tanb**2)*G_F*Zhb) )

  L=h0*SQRT(Zht*Zhb)

  g1 = 4._dp * SQRT2 * G_F *(MZ2-MW2)
  g2 = 4._dp * SQRT2 * G_F * MW2
  g = 2._dp * sqrt2 * G_F * mZ2
  C2TW = mW2 / mZ2
  s2tw = 1._dp - C2TW 
  !----------------------------------------
  ! calculate mu(QSTSB)
  !----------------------------------------
  mu = mu * Sqrt(Zht*Zhb) &
    &    * (1._dp + oo16pi2 * ( L**2+(mf_l2(3)/H2Q**2-g1-3._dp*g2)/2._dp) &
    &                                 * LOG(QSTSB/MZ2) )

 Al = Al_in + 3._dp*oo16pi2*(HTq**2*AT+HBq**2*AB)*T

  If(L*k*tanbq*mu == 0._dp)Then
    kont = -214
    Call AddError(214)
    Return
  End If

 !   Trig. functions of beta
  cb = 1./Sqrt(1.+tanbq**2)
  sb = tanbq*cb
  s2 = 2._dp*sb*cb

  nu = k/l*mu
  B = Al+nu

 !   Weak angle theta_W (S2TW  =  sin(theta_W)):
 !   Approximate value for the MSSM-like CP odd Higgs mass squared MA^2:
  MA2 = mu*B*(tanbq + 1._dp / tanbq)
 !   Approximate value for the Singlet-like CP even Higgs mass squared:
  MS2 = Max(nu*(4._dp*nu+Ak),mZ2)
 !   Approximate value for the Singlet-like CP odd Higgs mass squared:
  MP2 = -3._dp*nu*Ak
 !   Approximate value for the off-diag. CP odd mass matrix element:
  M12 = l*Sqrt(h1q**2+h2q**2)*(Al-2._dp*nu)
 !   Alphas at MT and QSTSB (for 2 loop corrs.)

 !   One loop functions for stop/sbottom loop corrections
  rt = 3._dp * oo32pi2 * htq**2
  mtopq = htq * h1q
  s2t = 2._dp * Rstop(1,1) * Rstop(1,2)
  If(mStop2(1) /= mStop2(2))Then
    fmt = (mStop2(2)*Log(mStop2(2)/QSTSB)-mStop2(1)*Log(mStop2(1)/QSTSB)) &
        & / (mStop2(2)-mStop2(1)) - 1._dp
    gmt = s2t**2*((mStop2(2)+mStop2(1))/(mStop2(2)-mStop2(1)) &
        & *Log(mStop2(2) / mStop2(1))-2._dp)
    emt = -mtopq*s2t*Log(mStop2(2)/mStop2(1))
  Else
    fmt = Log(mStop2(1)/QSTSB)
    gmt = 0._dp
    emt = 0._dp
  End If

  rb = 3._dp * oo32pi2 * hbq**2
  mbotq = oosqrt2 * hbq * h2q
  s2b = 2._dp * Rsbot(1,1) * Rsbot(1,2)
  If(mSbot2(1) /= mSbot2(2))Then
    fmb = (mSbot2(2)*Log(mSbot2(2)/QSTSB)-mSbot2(1)*Log(mSbot2(1)/QSTSB)) &
        & / (mSbot2(2)-mSbot2(1)) - 1._dp
    gmb = s2b**2*((mSbot2(2)+mSbot2(1))/(mSbot2(2)-mSbot2(1))   &
        & * Log(mSbot2(2)/mSbot2(1))-2._dp)
    emb = -mbotq*s2b*Log(mSbot2(2)/mSbot2(1))
  Else
    fmb = Log(mSbot2(1)/QSTSB)
    gmb = 0._dp
    emb = 0._dp
  End If
  
 !  The subsequent shifts simplify the expressions for the 
 !  one loop rad.corrs. below.
 !  The parameter Al is defined at the scale QSTSB
  If (mStop2(1) /= mStop2(2)) Then
    Alshift = Al + 2._dp * rt *At  &
      & * ((mStop2(2)*Log(mStop2(2)/QSTSB)-mStop2(1)*Log(mStop2(1)/QSTSB)) &
      &    / (mStop2(2)-mStop2(1))-1.)
  Else
    Alshift = Al + 2._dp * rt * At * Log(mStop2(1)/QSTSB)
  End If
  If (mSbot2(1) /= mSbot2(2)) Then
    Alshift = Alshift + 2._dp * rb * Ab &
      & * ((mSbot2(2)*Log(mSbot2(2)/QSTSB)-mSbot2(1)*Log(mSbot2(1)/QSTSB)) &
      &   / (mSbot2(2)-mSbot2(1))-1.)
  Else
    Alshift = Alshift + 2._dp * rb * Ab * Log(mSbot2(1)/QSTSB)
  End If
  B = Alshift+nu

 !   Tree level CP-even Higgs mass matrix

  mdia(1) = g*h1q**2 + mu*B/tanbq
  mdia(2) = g*h2q**2 + mu*B*tanbq
  mdia(3) = 4._dp*nu**2 + L**2*Alshift*h1q*h2q/mu + Ak*nu
  moff(1) = (2._dp*L**2-g)*h1q*h2q - mu*B
  moff(2) = L*(2._dp*mu*h1q - (B+nu)*h2q)
  moff(3) = L*(2._dp*mu*h2q - (B+nu)*h1q)

 !   1-loop radiative corrections
  rdia(1) = - At**2*gmt + 4._dp*At*emt &
        &  + 4._dp*mtopq**2*Log(mStop2(1)*mStop2(2)/mtopq**4)
  rdia(2) = - mu**2*gmt 
  rdia(3) = -l**2*h2q**2*gmt
  roff(1) = mu*(At*gmt - 2._dp*emt)
  roff(2) = l*h2q*(At*gmt - 2._dp*emt)
  roff(3) = - l*h2q*mu*gmt + 4._dp*l*h2q*mu*fmt

  mdia = mdia + rt * rdia
  moff = moff + rt * roff

  rdia(1) = - mu**2*gmb
  rdia(2) = - Ab**2*gmb + 4._dp*Ab*emb &
        & + 4._dp*mbotq**2*Log(mSbot2(1)*mSbot2(2)/mbotq**4)
  rdia(3) = -l**2*h1q**2*gmb
  roff(1) = mu*(Ab*gmb - 2._dp*emb)
  roff(2) = - l*h1q*mu*gmb + 4._dp*l*h1q*mu*fmb
  roff(3) = l*h1q*(Ab*gmb - 2._dp*emb)
  
  mdia = mdia + rt * rdia
  moff = moff + rt * roff

 !   Corrections from higgs/stop/sbottom couplings from D-terms
  gt =  g/2._dp - 2._dp * g1 / 3._dp
  Xt = At - mu / tanbq
  Dt = (M2_Q-M2_U+gt*(h2q**2-h1q**2))/2._dp
  D1 = oo16pi2 * 3._dp * (-g/4._dp*emt + gt/2._dp*Dt/Xt*gmt)
  D2 = oo16pi2 * 3._dp * ( -gt*Dt/Xt*emt &
     &                    - g/2._dp*mtopq**2*Log(mStop2(1)*mStop2(2)/QSTSB**2))
  mdia(1) =  mdia(1) + 2._dp*At*D1 + 2.*D2
  mdia(2) =  mdia(2) + 2._dp*mu/tanbq*D1
  moff(1) =  moff(1) - 1._dp/tanbq*D2 - (mu+At/tanbq)*D1
  moff(2) =  moff(2) - l*h2q*D1
  moff(3) =  moff(3) + l*h2q/tanbq*D1

  gb =  g/2._dp - 2._dp * g1/3._dp
  Xb = At - mu * tanbq
  Db = (M2_Q-M2_D+gb*(h1q**2-h2q**2))/2._dp
  D1 = oo16pi2 * 3._dp * (-g/4._dp*emb + gb/2._dp*Db/Xb*gmb)
  D2 = oo16pi2 * 3._dp * (-gb*Db/Xb*emb &
     &                   -g/2._dp*mbotq**2*Log(mSbot2(1)*mSbot2(2)/QSTSB**2))
  mdia(1) =  mdia(1) + 2._dp*mu*tanbq*D1
  mdia(2) =  mdia(2) + 2._dp*Ab*D1 + 2._dp*D2
  moff(1) =  moff(1) - tanbq * D2 - (mu+Ab*tanbq)*D1
  moff(2) =  moff(2) + l*h1q*tanbq*D1
  moff(3) =  moff(3) - l*h1q*D1

 !   2-loop terms
  mdia(1) = mdia(1) + oo8pi2**2 * 3._dp *htq**4 &
          &  *(32._dp*pi * ALSQ - 1.5_dp * htq**2 - hbq**2/2._dp)*h1q**2*t**2
  mdia(2) = mdia(2) + oo8pi2**2 * 3. * hbq**4 &
          & * (32._dp*pi*ALSQ-1.5_dp*hbq**2-htq**2/2._dp)*h2q**2*t**2

 !   Take care of the Z factors
  mdia(1) = mdia(1)/Zht
  mdia(2) = mdia(2)/Zhb
  moff(1) = moff(1)/Sqrt(Zht*Zhb)
  moff(2) = moff(2)/Sqrt(Zht)
  moff(3) = moff(3)/Sqrt(Zhb)

 !   Leading-log electroweak contributions (1 loop):
  
 !   a) Sfermion contributions
  sferm = g * mZ2 /(12._dp *pi2) &
    & * ( (S2TW**2/6._dp+C2TW**2*1.5_dp)*(Log(M2_Q/mZ2)+2._dp *Log(M2_Q1/mZ2)) &
    &   + 4._dp/3._dp * S2TW**2 * (Log(M2_U/mZ2)+2._dp *Log(M2_U1/mZ2))        &
    &   + 1._dp/3._dp * S2TW**2 * (Log(M2_D/mZ2)+2._dp *Log(M2_D1/mZ2))        &
    &   + 0.5_dp * (S2TW**2+C2TW**2) * (Log(M2_L/mZ2)+2._dp *Log(M2_L1/mZ2))   &
    &   + S2TW**2 * (Log(M2_E/mZ2)+2._dp *Log(M2_E1/mZ2))   )

  mdia(1) = mdia(1) + sferm*sb**2
  mdia(2) = mdia(2) + sferm*cb**2
  moff(1) = moff(1) - sferm*sb*cb

 !    b) Chargino/neutralino contributions
  LM2 = Log(Max(M2**2,mZ2)/mZ2)
  Lmu = Log(Max(mu**2,MZ2)/MZ2)
  Lnu = Log(Max(4._dp*nu**2,MZ2)/MZ2)
  LM2mu = Log(Max(M2**2,mu**2,MZ2)/MZ2)
  Lmunu = Log(Max(mu**2,4._dp*nu**2,MZ2)/MZ2)

  rdia(1) = -8._dp/3._dp*g*MZ2*sb**2*C2TW**2*LM2                             &
    &+(2._dp*l*k*mu**2/tanbq-4._dp/3._dp*g*MZ2*sb**2* (C2TW**2+S2TW**2))*Lmu &
    &+ 2._dp*l*k*nu**2/tanbq*Lnu                                             &
    & - (g*MZ2*sb**2*(16._dp*S2TW**2-28._dp*S2TW+14._dp) &
    & + g*mu*nu/tanbq*(2._dp*S2TW-3._dp))*LM2mu- 3._dp*l*k*mu**2/tanbq*Lmunu
  rdia(2) = -8._dp/3._dp*g*MZ2*cb**2*C2TW**2*LM2                             &
    &+(2._dp*l*k*mu**2*tanbq-4._dp/3._dp*g*MZ2*cb**2* (C2TW**2+S2TW**2))*Lmu &
    &+ 2._dp*l*k*nu**2*tanbq*Lnu                                             &
    & - (g*MZ2*cb**2*(16._dp*S2TW**2-28._dp*S2TW+14._dp)                     &
    & + g*mu*nu*tanbq*(2._dp*S2TW-3._dp))*LM2mu- 3._dp*l*k*mu**2*tanbq*Lmunu
  rdia(3) = 8._dp*mu**2*(k**2-l**2)*Lmu - 16._dp*k**2*nu**2*Lnu
  roff(1) = 4._dp/3._dp*g*MZ2*s2*C2TW**2*LM2      &
    & +(2._dp/3._dp*g*MZ2*s2*(C2TW**2+S2TW**2)   &
    & -2._dp*l*k*mu**2)*Lmu- 2._dp*l*k*nu**2*Lnu &
    &+ (g*MZ2*s2*(2._dp*S2TW-5._dp)-l**2*MZ2*s2*(2._dp*S2TW-3._dp) &
    & + g*mu*nu*(2._dp*S2TW-3._dp))*LM2mu+ (3._dp*l*k*mu**2-l**2*MZ2*s2)*Lmunu
  roff(2) = MZ/Sqrt(g)*mu*(4._dp*(l*sb-k*cb)*(l**2*Lmu+k**2*Lnu) &
    & + g*(l*sb+k*cb)*(4._dp*S2TW-6._dp)*LM2mu                   &
    &                  - l*(16._dp*k**2*sb+2._dp*l**2*sb-6._dp*l*k*cb)*Lmunu)
  roff(3) = MZ/Sqrt(g)*mu*(4._dp*(l*cb-k*sb)*(l**2*Lmu+k**2*Lnu) &
    & + g*(l*cb+k*sb)*(4._dp*S2TW-6._dp)*LM2mu                     &
    & - l*(16._dp*k**2*cb+2._dp*l**2*cb-6._dp*l*k*sb)*Lmunu)

  mdia = mdia + oo16pi2 * rdia
  moff = moff + oo16pi2 * roff

 !    c) Higgs loop contributions
 !    (Only if all masses squared are positive, and only to the lighter
 !    CP-even doublet-like state)
  subdet = mdia(3)*(sb**2*mdia(1)+cb**2*mdia(2)+s2*moff(1)) &
       & -(sb*moff(2)+cb*moff(3))**2
  If(subdet>0._dp)Then
    P2 = Max(MA2+MP2,MZ2)
    P1 = Max((MA2*MP2-M12**2)/P2,MZ2)
    LA = Log(Max(MA2,MZ2)/MZ2)
    LS = Log(MS2/MZ2)
    LP = Log(P2/MZ2)
    LPP = Log(P2/P1)
    bos = MZ2/(16._dp*PI**2*g)*((g**2/3._dp*(S2TW**2*(4._dp+2._dp*s2**2)        &
      -S2TW*(4._dp+8._dp*s2**2)-33._dp/4._dp*s2**4+16._dp*s2**2+5._dp/4._dp)    &
      +g*l**2*(2._dp*S2TW*s2**2+5.5_dp*s2**4-15._dp/2._dp*s2**2-1._dp)          &
      +l**4*(-11._dp/4._dp*s2**4+5._dp/2._dp*s2**2+1._dp))*LA+(l**2*(l-k*s2)**2 &
      &+3._dp*l**2/MS2*(g+(l**2-g)*s2**2)*(2._dp*mu-s2*(Al+2._dp*nu))**2        &
      &-l**4/MS2**2*(2._dp*mu-s2*(Al+2._dp*nu))**4)*LS                          &
      &+(g**2/4._dp*(1._dp-s2**4)+g*l**2*(1._dp/2._dp*s2**4+1._dp/2._dp*s2**2   &
      & -1._dp) &
      +l**4*(-1._dp/4._dp*s2**4-1._dp/2._dp*s2**2+1._dp)+l**2*(l+k*s2)**2)*LP&
      &-((g-l**2)/2._dp*MP2/P2*s2**2*(1._dp-s2**2) &
      +(l*MA2*(l+k*s2)-l**2*(Al-2._dp*nu)**2-MP2/2._dp*(g*(1._dp-s2**2) &
      -l**2*(2._dp-s2**2)))**2/P2**2-l*MA2*MP2*(l+k*s2) &
      *(g*(1._dp-s2**2)-l**2*(2._dp-s2**2))/P2**2)*LPP)
    mdia(1) = mdia(1) + bos*sb**2
    mdia(2) = mdia(2) + bos*cb**2
    moff(1) = moff(1) + bos*sb*cb
  End If

 !   Diagonalization,  SPheno uses a different basis coompared to NMHDECAY
  mat3(1,1) = mdia(2)
  mat3(2,2) = mdia(1)
  mat3(3,3) = mdia(3)
  mat3(1,2) = moff(1)
  mat3(2,1) = mat3(1,2)
  mat3(1,3) = moff(3)
  mat3(3,1) = mat3(1,3)
  mat3(2,3) = moff(2)
  mat3(3,2) = mat3(2,3)
  call EigenSystem( mat3, mS02, RS0, kont, test)
 !   CP even pole masses
  Do I = 1,3
    If(mS02(I)>0._dp)Then
      If(mS02(I)>4._dp*MT**2)Then
        X = Sqrt(1._dp-4._dp*MT**2/mS02(I))
        BT = 2._dp-X*Log((1._dp+X)/(1._dp-X))
      Else
        X = Sqrt(mS02(I)/(4._dp*MT**2-mS02(I)))
        BT = 2._dp*(1._dp-ATAN(X)/X)
      End If
      If(mS02(I)>4._dp*MB**2)Then
        X = Sqrt(1._dp-4._dp*MB**2/mS02(I))
        BB = 2._dp-X*Log((1._dp+X)/(1._dp-X))
      Else
        X = Sqrt(mS02(I)/(4._dp*MB**2-mS02(I)))
        BB = 2._dp*(1._dp-ATAN(X)/X)
      End If
      mS02(I) = mS02(I) - 2._dp*rt*RS0(I,1)**2*(mS02(I)-4._dp*MT**2)*BT  &
        & - 2._dp*rb*RS0(I,2)**2*(mS02(I)*Log(MT**2/MB**2)               &
        & +(mS02(I)-4._dp*MB**2)*BB)
    End If
  End Do
  If(mS02(1)<0._dp)Then
    kont = -215
    Call AddError(215)
  Else
    Do I = 1,3
      mS0(I) = Sqrt(mS02(I))
    End Do
  End If

 !   CP-odd Higgs mass matrix including
 !   1-loop top/bottom radiative corrections
  mp(1) = mu*B*(h1q/(Zhb*h2q)+h2q/(Zht*h1q))
  mp(2) = L**2*h1q*h2q/mu*(3._dp*nu+B)-3._dp*Ak*nu
  mp(3) = L*(Alshift-2._dp*nu)*Sqrt(h1q**2/Zhb+h2q**2/Zht)
  
 !   Dominant chargino/neutralino contribution ~<S>**2
  mp(1) = mp(1)+ g*mu*nu*(3._dp-2._dp*S2TW)*LM2mu*(tanbq+1/tanbq) * oo16pi2

 !   Diagonalization
  mP02(1) = mZ2  ! goldstone boson
  if (mp(3).eq.0._dp) then
   mP02(2) = Min(mp(1), mp(2))
   mP02(3) = Max(mp(1), mp(2))
   RP0 = Real(id3c, dp)
   RP0(1,1) = -cb
   RP0(1,2) = sb
   RP0(2,1) = sb
   RP0(2,2) = cb
  else
   mp02(2) = 0.5_dp * (mp(1)+mp(2) - Sqrt( (mp(1)-mp(2))**2 + 4._dp * mp(3)**2))
   mp02(3) = 0.5_dp * (mp(1)+mp(2) + Sqrt( (mp(1)-mp(2))**2 + 4._dp * mp(3)**2))
   RP0(1,1) = -cb
   RP0(1,2) = sb
   RP0(1,3) = 0._dp
   RP0(2,2) = - mp(3) / Sqrt((mp(1) - mP02(2))**2 + mp(3)**2)
   RP0(2,3) = (mp(1) - mP02(2) ) / Sqrt((mp(1) - mP02(2))**2 + mp(3)**2)
   RP0(3,3) = RP0(2,2)
   RP0(3,2) = - RP0(2,3)
   RP0(2,1) = RP0(2,2) * sb
   RP0(2,2) = RP0(2,2) * cb
   RP0(3,1) = RP0(3,2) * sb
   RP0(3,2) = RP0(3,2) * cb
  end if

 !   CP odd pole masses
  Do I = 2,3
    If(mP02(I)>0._dp)Then
      If(mP02(I)>4._dp*MT**2)Then
        X = Sqrt(1._dp-4._dp*MT**2/mP02(I))
        BT = 2._dp-X*Log((1._dp+X)/(1._dp-X))
      Else
        X = Sqrt(mP02(I)/(4._dp*MT**2-mP02(I)))
        BT = 2._dp*(1._dp-ATAN(X)/X)
      End If
      If(mP02(I)>4._dp*MB**2)Then
        X = Sqrt(1._dp-4._dp*MB**2/mP02(I))
        BB = 2._dp-X*Log((1._dp+X)/(1._dp-X))
      Else
        X = Sqrt(mP02(I)/(4._dp*MB**2-mP02(I)))
        BB = 2._dp*(1._dp-ATAN(X)/X)
      End If
      mP02(I) = mP02(I) - 2._dp*rt*RP0(I,1)**2*cb**2*(mP02(I)-4._dp*MT**2)*BT&
        &- 2._dp*rb*RP0(I,1)**2*sb**2*(mP02(I)*Log(MT**2/MB**2) &
        +(mP02(I)-4._dp*MB**2)*BB)
    End If
  End Do

  If(mP02(2)<0._dp)Then
    kont = -216
    Call AddError(216)
  Else
    Do I = 1,3
      mP0(I) = Sqrt(mP02(I))
    End Do
  End If

  !------------------------------------------------------------
  ! rotate RP0 to electroweak basis including Goldstone bosons
  !------------------------------------------------------------
   RSpm(1,1) = - cb
   RSpm(1,2) = sb 
   RSpm(2,1) = sb 
   RSpm(2,2) = cb

 !   Charged Higgs mass including 1-loop radiative corrections
   mSpm2(2) = ((g2/2._dp-L**2)*h1q*h2q+mu*B) &
    *(h1q**2*Zht+h2q**2*Zhb)/(h1q*h2q*Zht*Zhb) &
    +oo16pi2*3._dp*htq**2*hbq**2/(Sqrt2*G_F)*t&
    &+oo16pi2*g2*MW**2/3._dp*(12._dp*t+3._dp*(5*S2TW/C2TW-1._dp)*LM2mu &
    &                         -4._dp*LM2-2._dp*Lmu)
  
 !   Dominant chargino/neutralino contribution ~<S>**2
  mSpm2(2) = mSpm2(2)+ oo16pi2*g*mu*nu*(3._dp-2._dp*S2TW)*LM2mu*(tanbq+1/tanbq)
 !   Charged Higgs pole mass (in the LLA only)
  If(mSpm2(2)>MT**2)Then
    mSpm2(2) = mSpm2(2)/ (1._dp-(2._dp*rt*cb**2+2._dp*rb*sb**2) &
           &  *Log(mSpm2(2)/MT**2))
  End If
  If(mSpm2(2)<0._dp)Then
    kont = -217
    Call AddError(217)
  Else
    mSpm(2) = Sqrt(mSpm2(2))
  End If
  mSpm(1) = mW
  mSpm2(1) = mW2

 End Subroutine ScalarMassNMSSMeff
  

 Subroutine SdownMass3Lam(M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, lamp, &
   &  Alam, g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates squark masses and mixing angles for complex parameters
 ! for 3-generation epsilon model
 ! input:
 !  M_L2(i,j) . left sfermion squared mass matrix
 !  M_R2(i,j) . right sfermion squared mass matrix
 !  Af(i,j) ... trilinear A-parameter 
 !  mu ........ mu-parameter
 !  vevSM(i) ... i=1 v_d
 !              i=2 v_u
 !  Yuk ....... Yukawa coupling
 !  T3 ........ weak isospin of the left sfermion
 !  Yl ........ hypercharge of the left sfermion
 !  Yr ........ hypercharge of the right sfermion
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  mf ........ fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod: 25.08.2001
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: vevSM(2), vevL(3), g, gp
  Real(Dp), Intent(out) :: msf(6), msf2(6)
  Complex(Dp), Intent(in) :: M_L2(3,3), M_R2(3,3), Af(3,3), Yuk(3,3) &
                              & , bi(4), lamp(3,3,3), Alam(3,3,3)
  Complex(Dp), Intent(out) :: Rsf(6, 6)
  Integer, Intent(inout) :: kont
 
  Integer :: i1, i2, i3, ierr
  Real(Dp) :: T3, Yl, Yr, vev2, test(2)
  Complex(Dp) :: mat6(6,6),YukA(3,3),AfA(3,3),off(3,3), yukp(3,3)

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SdownMass3Lam'

  kont = 0

  vev2 = 0.25_dp * ( vevSM(1)**2 - vevSM(2)**2   &
       &           + vevL(1)**2 + vevL(2)**2 + vevL(3)**2 )
  AfA = Af
  Call Adjungate(AfA)
  yukp = vevSM(1) * yuk 
  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     yukp(i1,i2) = yukp(i1,i2) + vevL(i3) * lamp(i3,i1,i2)
    End Do
   End Do
  End Do
  YukA = Yukp
  Call Adjungate(YukA)

  T3 = -0.5_dp
  Yl = 1._dp / 3._dp
  Yr = 2._dp / 3._dp

  mat6(1:3,1:3) = M_L2 + 0.5_dp * Matmul(YukA,Yukp) &
      &         + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2 * id3C
  mat6(4:6,4:6) = M_R2 + 0.5_dp * Matmul(Yukp,YukA) &
      &         - 0.5_dp * Yr * gp**2 * vev2 * id3C
  off = (vevSM(1) * AfA - bi(1) * vevSM(2) * YukA ) * oosqrt2
  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     off(i1,i2) = off(i1,i2) + ( Conjg(Alam(i3,i1,i2))*vevL(i3) &
                &   + vevSM(2) * bi(1+i3) * Conjg(lamp(i3,i1,i2) ) ) * oosqrt2
    End Do
   End Do
  End Do

   mat6(1:3,4:6) = off

   Call Adjungate(off)
   mat6(4:6,1:3) = off


  Call EigenSystem(mat6,msf2,Rsf,ierr, test)

   If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
 If ((ierr.Ne.0).And.(ErrorLevel.Ge.0)) Then
   Write(ErrCan,*) 'Diagonalization did not work in routine SdownMass3Lam!'
    Write(ErrCan, * ) 'msf2 ', msf2
    Write(ErrCan, * ) 'M_L2, ',M_L2
    Write(ErrCan, * ) 'M_R2, ',M_R2
    Write(ErrCan, * ) 'A_f, ',Af
    Write(ErrCan, * ) 'Y_f, ',Yuk
    Write(ErrCan,*) 'A_lam',Alam
    Write(ErrCan,*) 'lamp',lamp
    Write(ErrCan, * ) 'bi, vevSM, vevL, g, gp'
    Write(ErrCan, * ) bi, vevSM, vevL, g, gp
   kont = ierr
    If (ErrorLevel.Eq.2) Call TerminateProgram
   Iname = Iname - 1
   Return
  End If

  Do i1=1,6
   If (mSf2(i1).Gt.0._dp) Then
    mSf(i1) = Sqrt( mSf2(i1) )
   Else If (ErrorLevel.Ge.0) Then
    Write(ErrCan, * ) 'Warning from routine SdownMass3Lam!'
    Write(ErrCan, * ) 'in the calculation of the masses'
    Write(ErrCan, * ) 'occurred a negative mass squared!!!'
    Write(ErrCan, * ) 'msf2 ', msf2
    Write(ErrCan, * ) 'M_L2, ',M_L2
    Write(ErrCan, * ) 'M_R2, ',M_R2
    Write(ErrCan, * ) 'A_f, ',Af
    Write(ErrCan, * ) 'Y_f, ',Yuk
    Write(ErrCan,*) 'A_lam',Alam
    Write(ErrCan,*) 'lamp',lamp
    Write(ErrCan, * ) 'bi, vevSM, vevL, g, gp'
    Write(ErrCan, * ) bi, vevSM, vevL, g, gp

    If (ErrorLevel.Eq.2) Call TerminateProgram

    kont = -219
    Call AddError(219)
    msf(i1) = 0._dp
   End If
  End Do
  Iname = Iname - 1

 End Subroutine SdownMass3Lam


 Subroutine SfermionMass1eps1(M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3,  &
   & Yl, Yr, g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates sfermion masses and mixing angles for complex parameters
 ! for 1 generation
 ! input:
 !  M_L2 ....... left sfermion mass squared
 !  M_R2 ....... left sfermion mass squared
 !  Af ......... trilinear A-parameter
 !  bi(i) ...... i=1 mu-parameter
 !              i=2 epsilon_1
 !  vevSM(i) ... i=1 v_d
 !              i=2 v_u
 !  vevL ...... v_3
 !  Yuk ....... Yukawa coupling
 !  T3 ........ weak isospin of the left sfermion
 !  Yl ........ hypercharge of the left sfermion
 !  Yr ........ hypercharge of the right sfermion
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  mf ........ fermion mass
 !  mf ........ fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod, 9.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: M_L2, M_R2, T3, Yl, Yr, vevSM(2), g, gp, vevL
  Real(Dp), Intent(out) :: msf(2), msf2(2)
  Complex(Dp), Intent(in) :: Af, Yuk, bi(2)
  Complex(Dp), Intent(out) :: Rsf(2, 2)
  Integer, Intent(inout) :: kont

  Real(Dp) :: vev2, diag(2), abs_offdiag, trace, det2, det,     &
    &  nen2, nen, cost, sint
  Complex(Dp) :: offdiag, phase

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SfermionMass1eps1'
  kont = 0

  vev2 = 0.25_dp * (vevSM(1)**2 - vevSM(2)**2 + vevL**2)
  If (T3.Gt.0) Then
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevSM(2)**2        &
     &     + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevSM(2)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = (Conjg(Af) * vevSM(2)                          &
     &       - Yuk * (bi(1) * vevSM(1) - bi(2) * vevL) ) * oosqrt2

  Else
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevSM(1)**2        &
     &     + (T3 * g**2 - 0.5_dp * YL * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevSM(1)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = (Conjg(Af) * vevSM(1) - bi(1) * vevSM(2) * Yuk) * oosqrt2
  End If
  !------------
  ! Masses
  !------------
  abs_offdiag = Abs(offdiag)
  trace = diag(1) + diag(2)
  det2 = (diag(1) - diag(2) )**2 + 4._dp * abs_offdiag**2

  det = Sqrt(det2)
  msf2(1) = 0.5_dp * (trace-det)
  msf2(2) = 0.5_dp * (trace+det)
  If (msf2(1) .Le.0.) Then
   Write(ErrCan, * ) 'Warning from routine SfermionMass1!'
   Write(ErrCan, * ) 'in the calculation of the masses'
   Write(ErrCan, * ) 'occurred a negative mass squared!!!'
   Write(ErrCan, * ) 'msf2 ', msf2
   !       write(10,*) 'setting it to 0.'
   !       msf2(1) = 0._dp
   Write(ErrCan, * ) 'M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3, Yl, Yr, g, gp'
   Write(ErrCan, * ) M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3, Yl, Yr, g, gp

   kont = -220
   Call AddError(220)
   msf(1) = 0._dp
   msf(2) = Sqrt(msf2(2) )
  Else
   msf(1) = Sqrt(msf2(1) )
   msf(2) = Sqrt(msf2(2) )
  End If

  !---------------
  ! mixing matrix
  !---------------
  If (offdiag.Eq. (0._dp, 0._dp) ) Then
   phase = (0._dp, 0._dp)
  Else
   phase = (0._dp, 0.5_dp) * Arg(offdiag)
  End If

  nen2 = (diag(1) - msf2(1) )**2 + abs_offdiag**2

  nen = Sqrt(nen2)
  If (nen.Lt.1.d-7) Then
   cost = 1._dp
   sint = 0._dp
  Else
   cost = - abs_offdiag / nen
   sint = (diag(1) - msf2(1) ) / nen
  End If

  Rsf(1, 1) = cost * Exp(phase)
  Rsf(1, 2) = sint * Exp( - phase)
  Rsf(2, 1) = - sint * Exp(phase)
  Rsf(2, 2) = cost * Exp( - phase)

  Iname = Iname - 1

 End Subroutine SfermionMass1eps1

 Subroutine SfermionMass1eps3(M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3,  &
   & Yl, Yr, g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates sfermion masses and mixing angles for complex parameters
 ! for 1 generation
 ! input:
 !  M_L2 ....... left sfermion mass squared
 !  M_R2 ....... left sfermion mass squared
 !  Af ......... trilinear A-parameter
 !  bi(i) ...... i=1 mu-parameter
 !               i=2 epsilon_1
 !               i=3 epsilon_2
 !               i=4 epsilon_3
 !  vevSM(i) ... i=1 v_d
 !               i=2 v_u
 !  vevL(i) .... i=1 v_1
 !               i=2 v_2
 !               i=3 v_3
 !  Yuk ........ Yukawa coupling
 !  T3 ......... weak isospin of the left sfermion
 !  Yl ......... hypercharge of the left sfermion
 !  Yr ......... hypercharge of the right sfermion
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  mf ......... fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod, 10.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: M_L2, M_R2, T3, Yl, Yr, vevSM(2), g, gp, &
                               &  vevL(3)
  Real(Dp), Intent(out) :: msf(2), msf2(2)
  Complex(Dp), Intent(in) :: Af, Yuk, bi(4)
  Complex(Dp), Intent(out) :: Rsf(2, 2)
  Integer, Intent(inout) :: kont

  Real(Dp) :: vev2, diag(2), abs_offdiag, trace, det2, det,     &
    &  nen2, nen, cost, sint
  Complex(Dp) :: offdiag, phase

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SfermionMass1eps3'
  kont = 0

  vev2 = 0.25_dp * (vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL) )
  If (T3.Gt.0) Then
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevSM(2)**2        &
     &     + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevSM(2)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = - bi(1) * vevSM(1) + bi(2) * vevL(1) + bi(3) * vevL(2) &
           & + bi(4) * vevL(3)
   offdiag = (Conjg(Af) * vevSM(2) + offdiag * Yuk) * oosqrt2

  Else
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevSM(1)**2        &
     &     + (T3 * g**2 - 0.5_dp * YL * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevSM(1)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = (Conjg(Af) * vevSM(1) - bi(1) * vevSM(2) * Yuk) * oosqrt2
  End If

  !------------
  ! Masses
  !------------
  abs_offdiag = Abs(offdiag)
  trace = diag(1) + diag(2)
  det2 = (diag(1) - diag(2) )**2 + 4._dp * abs_offdiag**2

  det = Sqrt(det2)
  msf2(1) = 0.5_dp * (trace-det)
  msf2(2) = 0.5_dp * (trace+det)
  If (msf2(1) .Le.0.) Then
   Write(ErrCan, * ) 'Warning from routine SfermionMass1!'
   Write(ErrCan, * ) 'in the calculation of the masses'
   Write(ErrCan, * ) 'occurred a negative mass squared!!!'
   Write(ErrCan, * ) 'msf2 ', msf2
   !       write(10,*) 'setting it to 0.'
   !       msf2(1) = 0._dp
   Write(ErrCan, * ) 'M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3, Yl, Yr, g, gp'
   Write(ErrCan, * ) M_L2, M_R2, Af, bi, vevSM, vevL, Yuk, T3, Yl, Yr, g, gp

   kont = -221
   Call AddError(221)
   msf(1) = 0._dp
   msf(2) = Sqrt(msf2(2) )
  Else
   msf(1) = Sqrt(msf2(1) )
   msf(2) = Sqrt(msf2(2) )
  End If

  !---------------
  ! mixing matrix
  !---------------
  If (offdiag.Eq. (0._dp, 0._dp) ) Then
   phase = (1._dp, 0._dp)
  Else If (Aimag(offdiag).eq.0._dp) Then
   If (Real(offdiag).gt.0._dp) Then
    phase = (1._dp, 0._dp)
   Else
    phase = (-1._dp, 0._dp)
   End If
  Else
   phase = (0._dp, 1._dp) * Arg(offdiag)
   Phase =  Exp(phase)
  End If

  nen2 = (diag(1) - msf2(1) )**2 + abs_offdiag**2
  nen = Sqrt(nen2)
  If (nen.Lt.1.e-7_dp) Then
   cost = 1._dp
   sint = 0._dp
  Else
   cost = - abs_offdiag / nen
   sint = (diag(1) - msf2(1) ) / nen
  End If

  Rsf(1, 1) = cost * phase
  Rsf(1, 2) = sint 
  Rsf(2, 1) = - sint
  Rsf(2, 2) = cost / phase

  Iname = Iname - 1

 End Subroutine SfermionMass1eps3

 Subroutine SfermionMass1mssm(M_L2, M_R2, Af, mu, vevs, Yuk, T3, Yl, Yr, &
   &  g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates sfermion masses and mixing angles for complex parameters
 ! for 1 generation
 ! input:
 !  M_L2 ...... left sfermion mass squared
 !  M_R2 ...... right sfermion mass squared
 !  Af ........ trilinear A-parameter
 !  mu ........ mu-parameter
 !  vevs(i) ... i=1 v_d
 !              i=2 v_u
 !  Yuk ....... Yukawa coupling
 !  T3 ........ weak isospin of the left sfermion
 !  Yl ........ hypercharge of the left sfermion
 !  Yr ........ hypercharge of the right sfermion
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  mf ........ fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod: 9.11.2000
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: M_L2, M_R2, T3, Yl, Yr, vevs(2), g, gp
  Real(Dp), Intent(out) :: msf(2), msf2(2)
  Complex(Dp), Intent(in) :: Af, Yuk, mu
  Complex(Dp), Intent(out) :: Rsf(2, 2)
  Integer, Intent(inout) :: kont
 
  Real(Dp) :: vev2, diag(2), abs_offdiag, trace, det2, det,     &
    &  nen2, nen, cost, sint
  Complex(Dp) :: offdiag, phase

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SfermionMass1MSSM'

  kont = 0

  vev2 = 0.25_dp * (vevs(1)**2 - vevs(2)**2)
  If (T3.Gt.0) Then
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevs(2)**2        &
     &     + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevs(2)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = (Conjg(Af) * vevs(2) - mu * vevs(1) * conjg(Yuk)) * oosqrt2

  Else
   diag(1) = M_L2 + 0.5_dp * Abs(Yuk)**2 * vevs(1)**2        &
     &     + (T3 * g**2 - 0.5_dp * YL * gp**2) * vev2
   diag(2) = M_R2 + 0.5_dp * Abs(Yuk)**2 * vevs(1)**2        &
     &     - 0.5_dp * Yr * gp**2 * vev2

   offdiag = (Conjg(Af) * vevs(1) - Conjg(Yuk) * mu * vevs(2) ) * oosqrt2
  End If
  !------------
  ! Masses
  !------------
  abs_offdiag = Abs(offdiag)
  trace = diag(1) + diag(2)
  det2 = (diag(1) - diag(2) )**2 + 4._dp * abs_offdiag**2

  det = Sqrt(det2)
  msf2(1) = 0.5_dp * (trace-det)
  msf2(2) = 0.5_dp * (trace+det)
  If (msf2(1) .Le.0.) Then
   Write(ErrCan, * ) 'Warning from routine SfermionMass1mssm!'
   Write(ErrCan, * ) 'in the calculation of the masses'
   Write(ErrCan, * ) 'occurred a negative mass squared!!!'
   Write(ErrCan, * ) 'msf2 ', msf2
   Write(ErrCan, * ) 'M_L2, M_R2, Af, mu, vevs, Yuk, T3, Yl, Yr, g, gp'
   Write(ErrCan, * ) M_L2, M_R2, Af, mu, vevs, Yuk, T3, Yl, Yr, g, gp

   kont = - 222
   Call AddError(222)
   msf = Sqrt( Abs(msf2) )

  Else
   msf = Sqrt(msf2 )
  End If

  !---------------
  ! mixing matrix
  !---------------
  If (offdiag.Eq. (0._dp, 0._dp) ) Then
   phase = (1._dp, 0._dp)
  Else If (Aimag(offdiag).eq.0._dp) Then
   If (Real(offdiag).gt.0._dp) Then
    phase = (1._dp, 0._dp)
   Else
    phase = (-1._dp, 0._dp)
   End If
  Else
   phase = (0._dp, 1._dp) * Arg(offdiag)
   Phase =  Exp(phase)
  End If

  nen2 = (diag(1) - msf2(1) )**2 + abs_offdiag**2

  nen = Sqrt(nen2)
  If (nen.Lt.1.e-9_dp) Then
   cost = 1._dp
   sint = 0._dp
  Else
   cost = - abs_offdiag / nen
   sint = (diag(1) - msf2(1) ) / nen
  End If

  Rsf(1, 1) = cost * phase
  Rsf(1, 2) = sint 
  Rsf(2, 1) = - sint
  Rsf(2, 2) = cost / phase

  Iname = Iname - 1

 End Subroutine SfermionMass1mssm

 Subroutine SfermionMass3mssm(M_L2, M_R2, Af, mu, vevs, Yuk, T3, Yl, Yr, &
   &  g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates sfermion masses and mixing angles for complex parameters
 ! for 1 generation
 ! input:
 !  M_L2(i,j) . left sfermion squared mass matrix
 !  M_R2(i,j) . right sfermion squared mass matrix
 !  Af(i,j) ... trilinear A-parameter 
 !  mu ........ mu-parameter
 !  vevs(i) ... i=1 v_d
 !              i=2 v_u
 !  Yuk ....... Yukawa coupling
 !  T3 ........ weak isospin of the left sfermion
 !  Yl ........ hypercharge of the left sfermion
 !  Yr ........ hypercharge of the right sfermion
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  mf ........ fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod: 9.11.2000
 ! 16.03.03: including a test, if there is really generation mixing
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: T3, Yl, Yr, vevs(2), g, gp
  Real(Dp), Intent(out) :: msf(6), msf2(6)
  Complex(Dp), Intent(in) :: M_L2(3,3), M_R2(3,3), Af(3,3), Yuk(3,3), mu
  Complex(Dp), Intent(out) :: Rsf(6, 6)
  Integer, Intent(inout) :: kont
 
  Integer :: i1, i2, ierr
  Real(Dp) :: vev2, test(2), Rsfa(6,6), m2(2), m22(2), Ml, Mr
  Complex(Dp) :: mat6(6,6), off(3,3), Rsf2(2,2), A, Y &
     & , vec(6), YukT(3,3), YukC(3,3)

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SfermionMass3MSSM'

  kont = 0

  vev2 = 0.25_dp * (vevs(1)**2 - vevs(2)**2)

  !-----------------------------------------
  ! test for generation mixing
  !-----------------------------------------
  test(1) = 0._dp
  Do i1=1,2
   Do i2=i1+1,3
    test(1) = test(1) + Abs(M_L2(i1,i2))+Abs(M_r2(i1,i2))+Abs(af(i1,i2))  &
            &         + Abs(Af(i2,i1)) + Abs(Yuk(i1,i2)) + Abs(yuk(i2,i1))
   End Do
  End Do

  If (test(1).Eq.0._dp) Then
   Rsf = 0._dp
   msf2 = 10000._dp
   Do i1=1,3
    Ml = M_L2(i1,i1)
    Mr = M_R2(i1,i1)
    A = Af(i1,i1)
    Y = Yuk(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevs, Y, T3, Yl, Yr, g, gp, ierr &
                    &, m2, m22, Rsf2 )

    If ((ierr.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SfermionMass3mssm!'
     msf2(2*i1-1:2*i1) = m22
     Write(ErrCan, * ) 'msf2 ', msf2
     Write(ErrCan, * ) 'M_L2, ',M_L2
     Write(ErrCan, * ) 'M_R2, ',M_R2
     Write(ErrCan, * ) 'A_f, ',Af
     Write(ErrCan, * ) 'Y_f, ',Yuk
     Write(ErrCan, * ) 'mu, vevs, T3, Yl, Yr, g, gp'
     Write(ErrCan, * ) mu, vevs, T3, Yl, Yr, g, gp
     kont = ierr
     If (ErrorLevel.Eq.2) Call TerminateProgram
     Iname = Iname - 1
     Return
    End If
    msf(2*i1-1:2*i1) = m2
    Rsf(2*i1-1,i1) = Rsf2(1,1)
    Rsf(2*i1-1,i1+3) = Rsf2(1,2)
    Rsf(2*i1,i1) = Rsf2(2,1)
    Rsf(2*i1,i1+3) = Rsf2(2,2)
   End Do

   !----------------------
   ! reordering
   !----------------------
   Do i1=1,5
    Do i2=i1+1,6
     If (msf(i1).Gt.msf(i2)) Then
      test(1) = msf(i1)
      msf(i1) = msf(i2)
      msf(i2) = test(1)
      vec = Rsf(i1,:)
      Rsf(i1,:) = Rsf(i2,:)
      Rsf(i2,:) = vec
     End If
    End Do
   End Do
   msf2 = msf**2

  Else  ! there is really generation mixing
   test = 0._dp

   YukT = Transpose(Yuk)
   YukC = Conjg(Yuk)

   If (T3.Gt.0) Then
    mat6(1:3,1:3) = M_L2 + 0.5_dp * vevs(2)**2 * Matmul(YukC,YukT) &
        &         + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2 * id3C
    mat6(4:6,4:6) = M_R2 + 0.5_dp * vevs(2)**2 * Matmul(YukT,YukC) &
       &         - 0.5_dp * Yr * gp**2 * vev2 * id3C
    off = (vevs(2) * Af - Conjg(mu) * vevs(1) * Yuk ) * oosqrt2
   Else
    mat6(1:3,1:3) = M_L2 + 0.5_dp * vevs(1)**2 * Matmul(YukC,YukT) &
        &         + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2 * id3C
    mat6(4:6,4:6) = M_R2 + 0.5_dp * vevs(1)**2 * Matmul(YukT,YukC) &
       &         - 0.5_dp * Yr * gp**2 * vev2 * id3C
    off = (vevs(1) * Af - Conjg(mu) * vevs(2) * Yuk ) * oosqrt2
   End If

   off = Conjg(off)
   mat6(1:3,4:6) = off
   Call Adjungate(off)
   mat6(4:6,1:3) = off      

   If (Maxval(Abs(Aimag(mat6))).Eq.0._dp) Then
    Call EigenSystem(Real(mat6,dp),msf2,Rsfa,ierr, test)
    Rsf = Rsfa
   Else
    Call EigenSystem(mat6,msf2,Rsf,ierr, test)
   End If
 
   If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) "T_3, Y_l",T3,yl
     Write(ErrCan,*) "M_L2",M_L2(1,:)
     Write(errcan,*) "    ",M_L2(2,:) 
     Write(errcan,*) "    ",M_L2(3,:) 
     Write(ErrCan,*) "M_R2",M_R2(1,:)
     Write(errcan,*) "    ",M_r2(2,:) 
     Write(errcan,*) "    ",M_r2(3,:)
     off =  Matmul(YukC,YukT)
     Write(Errcan,*) "Y^* Y^T",off(1,:)
     Write(errcan,*) "     ",off(2,:)
     Write(errcan,*) "     ",off(3,:)
     off =  Matmul(YukT,YukC)
     Write(Errcan,*) "Y^T Y^*",off(1,:)
     Write(errcan,*) "     ",off(2,:)
     Write(errcan,*) "     ",off(3,:)
     Write(ErrCan,*) "mat 1",mat6(1,:)
     Write(ErrCan,*) "mat 2",mat6(2,:)
     Write(ErrCan,*) "mat 3",mat6(3,:)
     Write(ErrCan,*) "mat 4",mat6(4,:)
     Write(ErrCan,*) "mat 5",mat6(5,:)
     Write(ErrCan,*) "mat 6",mat6(6,:)
     Write(ErrCan,*) " "
     Write(ErrCan,*) "msf2",msf2
     Write(ErrCan,*) "Rsf 1",rsf(1,:)
     Write(ErrCan,*) "Rsf 2",rsf(2,:)
     Write(ErrCan,*) "Rsf 3",rsf(3,:)
     Write(ErrCan,*) "Rsf 4",rsf(4,:)
     Write(ErrCan,*) "Rsf 5",rsf(5,:)
     Write(ErrCan,*) "Rsf 6",rsf(6,:)
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     ierr = 0
   End If
  

   If ((ierr.Ne.0).And.(ErrorLevel.Ge.0)) Then
     Write(ErrCan,*) 'Diagonalization did not work in routine SfermionMass3mssm!'
     Write(ErrCan, * ) 'msf2 ', msf2
     Write(ErrCan, * ) 'M_L2, ',M_L2
     Write(ErrCan, * ) 'M_R2, ',M_R2
     Write(ErrCan, * ) 'A_f, ',Af
     Write(ErrCan, * ) 'Y_f, ',Yuk
     Write(ErrCan, * ) 'mu, vevs, T3, Yl, Yr, g, gp'
     Write(ErrCan, * ) mu, vevs, T3, Yl, Yr, g, gp
     kont = ierr
     If (ErrorLevel.Eq.2) Call TerminateProgram
     Iname = Iname - 1
     Return
   End If

   Do i1=1,6
    If (mSf2(i1).Gt.0._dp) Then
     mSf(i1) = Sqrt( mSf2(i1) )
    Else
     kont = -223
     Call AddError(223)
     msf(i1) = 0._dp
     If (ErrorLevel.Ge.0) Then
      Write(ErrCan, * ) 'Warning from routine SfermionMass3mssm!'
      Write(ErrCan, * ) 'in the calculation of the masses'
      Write(ErrCan, * ) 'occurred a negative mass squared!!!'
      Write(ErrCan, * ) 'msf2 ', msf2
      Write(ErrCan, * ) 'M_L2, ',M_L2
      Write(ErrCan, * ) 'M_R2, ',M_R2
      Write(ErrCan, * ) 'A_f, ',Af
      Write(ErrCan, * ) 'Y_f, ',Yuk
      Write(ErrCan, * ) 'mu, vevs, T3, Yl, Yr, g, gp'
      Write(ErrCan, * ) mu, vevs, T3, Yl, Yr, g, gp

      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
    End If
   End Do
  End If ! test for generation mixing

  Iname = Iname - 1

 End Subroutine SfermionMass3mssm


 Subroutine SquarkMass3Eps(M_L2, M_R2, Af, bi, vevSM, T3, vevL, Yuk, &
   &  g, gp, kont, msf, msf2, Rsf)
 !----------------------------------------------------------------------
 ! calculates squark masses and mixing angles for complex parameters
 ! for 3-generation epsilon model
 ! input:
 !  M_L2(i,j) . left sfermion squared mass matrix
 !  M_R2(i,j) . right sfermion squared mass matrix
 !  Af(i,j) ... trilinear A-parameter 
 !  mu ........ mu-parameter
 !  vevSM(i) ... i=1 v_d
 !              i=2 v_u
 !  Yuk ....... Yukawa coupling
 !  T3 ........ weak isospin of the left sfermion
 !  Yl ........ hypercharge of the left sfermion
 !  Yr ........ hypercharge of the right sfermion
 !  g ......... SU(2) gauge coupling
 !  gp ........ U(1) gauge coupling
 !  mf ........ fermion mass
 !  msf(i) .... sfermion mass i=1,2 msf(1) < msf(2)
 !  msf2(i) ... sfermion mass squared i=1,2
 !  Rsf(i,j) .. sfermion mixing matrix, i,j =1,2
 ! written by Werner Porod: 25.08.2001
 !----------------------------------------------------------------------
 Implicit None
  Real(Dp), Intent(in) :: T3, vevSM(2), vevL(3), g, gp
  Real(Dp), Intent(out) :: msf(6), msf2(6)
  Complex(Dp), Intent(in) :: M_L2(3,3), M_R2(3,3), Af(3,3), Yuk(3,3) &
                              & , bi(4)
  Complex(Dp), Intent(out) :: Rsf(6, 6)
  Integer, Intent(inout) :: kont
 
  Integer :: i1, ierr
  Real(Dp) :: Yl, Yr, vev2, test(2)
  Complex(Dp) :: mat6(6,6),YukA(3,3),AfA(3,3),off(3,3)

  Iname = Iname + 1

  NameOfUnit(Iname) = 'SquarkMass3Eps'

  kont = 0

  vev2 = 0.25_dp * ( vevSM(1)**2 - vevSM(2)**2   &
       &           + vevL(1)**2 + vevL(2)**2 + vevL(3)**2 )
  AfA = Af
  Call Adjungate(AfA)
  YukA = Yuk
  Call Adjungate(YukA)

  If (T3.Gt.0) Then
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   mat6(1:3,1:3) = M_L2 + 0.5_dp * vevSM(2)**2 * Matmul(YukA,Yuk) &
       &         + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2 * id3C
   mat6(4:6,4:6) = M_R2 + 0.5_dp * vevSM(2)**2 * Matmul(Yuk,YukA) &
       &         - 0.5_dp * Yr * gp**2 * vev2 * id3C
   off = ( vevSM(2) * AfA &
       & - (bi(1) * vevSM(1) - Dot_product(vevL,bi(2:4)) ) * Yuk ) * oosqrt2
   mat6(1:3,4:6) = off
   Call Adjungate(off)
   mat6(4:6,1:3) = off

  Else
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   mat6(1:3,1:3) = M_L2 + 0.5_dp * vevSM(1)**2 * Matmul(YukA,Yuk) &
       &         + (T3 * g**2 - 0.5_dp * Yl * gp**2) * vev2 * id3C
   mat6(4:6,4:6) = M_R2 + 0.5_dp * vevSM(1)**2 * Matmul(Yuk,YukA) &
       &         - 0.5_dp * Yr * gp**2 * vev2 * id3C
   off = (vevSM(1) * AfA - bi(1) * vevSM(2) * Yuk ) * oosqrt2
   mat6(1:3,4:6) = off
   Call Adjungate(off)
   mat6(4:6,1:3) = off
  End If

  Call EigenSystem(mat6,msf2,Rsf,ierr, test)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If ((ierr.Ne.0).And.(ErrorLevel.Ge.0)) Then
   Write(ErrCan,*) 'Diagonalization did not work in routine SquarkMass3Eps!'
    Write(ErrCan, * ) 'msf2 ', msf2
    Write(ErrCan, * ) 'M_L2, ',M_L2
    Write(ErrCan, * ) 'M_R2, ',M_R2
    Write(ErrCan, * ) 'A_f, ',Af
    Write(ErrCan, * ) 'Y_f, ',Yuk
    Write(ErrCan, * ) 'bi, vevSM, vevL, T3, Yl, Yr, g, gp'
    Write(ErrCan, * ) bi, vevSM, vevL, T3, Yl, Yr, g, gp
    kont = ierr
    If (ErrorLevel.Eq.2) Call TerminateProgram
    Iname = Iname - 1
   Return
  End If

  Do i1=1,6
   If (mSf2(i1).Gt.0._dp) Then
    mSf(i1) = Sqrt( mSf2(i1) )
   Else If (ErrorLevel.Ge.0) Then
    Write(ErrCan, * ) 'Warning from routine SquarkMass3Eps!'
    Write(ErrCan, * ) 'in the calculation of the masses'
    Write(ErrCan, * ) 'occurred a negative mass squared!!!'
    Write(ErrCan, * ) 'msf2 ', msf2
    Write(ErrCan, * ) 'M_L2, ',M_L2
    Write(ErrCan, * ) 'M_R2, ',M_R2
    Write(ErrCan, * ) 'A_f, ',Af
    Write(ErrCan, * ) 'Y_f, ',Yuk
    Write(ErrCan, * ) 'bi, vevSM, vevL, T3, Yl, Yr, g, gp'
    Write(ErrCan, * ) bi, vevSM, vevL, T3, Yl, Yr, g, gp

    If (ErrorLevel.Eq.2) Call TerminateProgram

    kont = -224
    Call AddError(224)
    msf(i1) = 0._dp
   End If
  End Do
  Iname = Iname - 1

 End Subroutine SquarkMass3Eps

 Subroutine TreeMassesEps1(gp, g, vevSM, vevL, M1, M2, M3, mu, B              &
          &, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u       &
          &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N                        &
          &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup                       &
          &, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2, RSneut           &
          &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm                &
          &, GenerationMixing, kont)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the 3-generation epsilon model
 ! written by Werner Porod, 21.05.2001
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu(2), B(2)
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), vevL
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(3), mC2(3), mN(5), mN2(5)               &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6), mP0(3) &
                      &  , mP02(3), RP0(3,3), mS0(3), mS02(3), RS0(3,3)     &
                      &  , mSpm(4), mSpm2(4), mSlepton(6), mSlepton2(6)     &
                      &  , mSneut(3), mSneut2(3)
  Complex(dp), Intent(out) :: PhaseGlu, U(3,3), V(3,3), RSdown(6,6), N(5,5) &
                      &  , RSup(6,6), RSpm(4,4), RSlepton(6,6), Rsneut(3,3)

  Integer :: i1
  Real(dp) :: D_sneut, T3, Yl, Yr, Ml, Mr, msf(2), msf2(2)
  Real(dp) :: mStop2(2), mSbottom2(2), mT, Yukl
  Complex(dp) :: A, Y, Rsf(2,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesEps1'

  kont = 0
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, vevL, g, mf_l(3), mC, U, V, Yukl, kont)
  mC2 = mC**2
  call NeutralinoMass(M1,M2,mu,vevSM,vevL,g,gp,mN,N,kont)
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
!    Call SquarkMass3Eps(M2_Q, M2_D, A_d, mu, vevSM, T3, vevL, Y_d &
!                    &  , g, gp, kont, mSdown, mSdown2, Rsdown)
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
!    Call SquarkMass3Eps(M2_Q, M2_U, A_u, mu, vevSM, T3, vevL, Y_u &
!                    &  , g, gp, kont, mSup, mSup2, Rsup)
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mSneut = 0
    mSneut2 = 0
    Do i1=1,2
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -225
      Call AddError(225)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   mSlepton = 0
   mSlepton2 = 0
   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
   Do i1 = 1,2
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = A_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSlepton(2*(i1-1)+1:2*i1) = msf
    mSlepton2(2*(i1-1)+1:2*i1) = msf2
    RSlepton(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
  End If

  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   Call PseudoScalarMass(mu,B,vevSM,vevL,gp,g, mP0,mP02,RP0,kont)
  !---------------
  ! charged
  !---------------
!   write(*,*) A_l
   Call ChargedScalarMass(mu, B, vevSM, vevL, Real(M2_E(3,3),dp), A_l(3,3) &
                        &, Y_l(3,3), gp, g, mSpm, mSpm2, RSpm, kont)
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
   mStop2 = msf2
  Else
   mStop2 = mSup2(5:6) 
  End If
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3)
  Call ScalarMass(mu, B, vevSM, vevL, gp, g, mStop2, mT**2 &
                &, mS0, mS02, RS0, kont)

  Iname = Iname - 1

 End Subroutine TreeMassesEps1

 Subroutine TreeMassesEps3(gp, g, vevSM, vevL, M1, M2, M3, mu, B              &
                         &, M2_E, M2_L, A_l, Y_l                              &
                         &, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u              &
                         &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N         &
                         &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup        &
                         &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm &
                         &, GenerationMixing, kont)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the 3-generation epsilon model
 ! written by Werner Porod, 21.05.2001
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu(4), B(4)
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), vevL(3)
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(5), mC2(5), mN(7), mN2(7)               &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6), mP0(5) &
                      &  , mP02(5), RP0(5,5), mS0(5), mS02(5), RS0(5,5)     &
                      &  , mSpm(8), mSpm2(8)
  Complex(dp), Intent(out) :: PhaseGlu, U(5,5), V(5,5), RSdown(6,6)         &
                      &  , RSup(6,6), RSpm(8,8), N(7,7)

  Integer :: i1
  Real(dp) :: T3, Yl, Yr, Ml, Mr, msf(2), msf2(2)
  Real(dp) :: mStop2(2), mT, Yukl(3)
  Complex(dp) :: A, Y, Rsf(2,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesEps3'

  kont = 0
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, vevL, g, mf_l, mC, U, V, Yukl, kont)
  mC2 = mC**2
  Call NeutralinoMass7(M1, M2, mu, vevSM, vevL, g, gp, mN, N, kont)
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
    Call SquarkMass3Eps(M2_Q, M2_D, A_d, mu, vevSM, T3, vevL, Y_d &
                    &  , g, gp, kont, mSdown, mSdown2, Rsdown)
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Call SquarkMass3Eps(M2_Q, M2_U, A_u, mu, vevSM, T3, vevL, Y_u &
                    &  , g, gp, kont, mSup, mSup2, Rsup)
  Else
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
  End If

  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   Call PseudoScalarMass(M2_L,mu,B,vevSM,vevL,gp,g, mP0,mP02,RP0,kont)
  !---------------
  ! charged
  !---------------
!   write(*,*) A_l
   call ChargedScalarMass(mu, B, vevSM, vevL,M2_L,M2_E, A_l &
                        &, Y_l, gp, g, mSpm, mSpm2, RSpm, kont)
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
   mStop2 = msf2
  Else
   mStop2 = mSup2(5:6) 
  End If
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3)
  Call ScalarMass(M2_L, mu, B, vevSM, vevL, gp, g, mStop2, mT**2 &
                &, mS0, mS02, RS0, kont)

  Iname = Iname - 1

 End Subroutine TreeMassesEps3

 Subroutine TreeMassesLam3(gp, g, vevSM, vevL, M1, M2, M3, mu, B               &
                         &, M2_E, M2_L, A_l, Y_l, lam, A_lam                   &
                         &, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u, lamp, A_lamp &
                         &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N          &
                         &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup         &
                         &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm  &
                         &, GenerationMixing, kont)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the 3-generation epsilon model
 ! including trilinear couplings
 ! written by Werner Porod, 03.09.2001
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu(4), B(4)
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3), A_lam(3,3,3)  &
                         & , A_lamp(3,3,3)
  Complex(dp), Intent(inout) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), lam(3,3,3)    &
                         & , lamp(3,3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)  & 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), vevL(3)
  Logical, Intent(in) :: GenerationMixing
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(5), mC2(5), mN(7), mN2(7)               &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6), mP0(5) &
                      &  , mP02(5), RP0(5,5), mS0(5), mS02(5), RS0(5,5)     &
                      &  , mSpm(8), mSpm2(8)
  Complex(dp), Intent(out) :: PhaseGlu, U(5,5), V(5,5), N(7,7), RSdown(6,6) &
                      &  , RSup(6,6), RSpm(8,8)

  Integer :: i1
  Real(dp) :: T3, Yl, Yr, Ml, Mr, msf(2), msf2(2)
  Real(dp) :: mStop2(2), mT
  Complex(dp) :: A, Y, Rsf(2,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesLam3'

  kont = 0
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, vevL, g, y_L, lam, mf_l, mC, U, V, kont)
  mC2 = mC**2
  Call NeutralinoMass7(M1, M2, mu, vevSM, vevL, g, gp, mN, N, kont)
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !--------------
   ! D-squarks
   !--------------
    Call SdownMass3Lam(M2_Q, M2_D, A_d, mu, vevSM, vevL, Y_d, lamp &
                    & , A_lamp, g, gp, kont, mSdown, mSdown2, Rsdown)
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Call SquarkMass3Eps(M2_Q, M2_U, A_u, mu, vevSM, T3, vevL, Y_u &
                    &  , g, gp, kont, mSup, mSup2, Rsup)
  Else
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
  End If

  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   Call PseudoScalarMass(M2_L,mu,B,vevSM,vevL,gp,g, mP0,mP02,RP0,kont)
  !---------------
  ! charged
  !---------------
!   write(*,*) A_l
   Call ChargedScalarMass(mu, B, vevSM, vevL, M2_L, M2_E, A_l &
                        &, Y_l, lam, A_lam, gp, g, mSpm, mSpm2, RSpm, kont)
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, vevL, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
   mStop2 = msf2
  Else
   mStop2 = mSup2(5:6) 
  End If
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3)
  Call ScalarMass(M2_L, mu, B, vevSM, vevL, gp, g, mStop2, mT**2 &
                 &, mS0, mS02, RS0, kont)

  Iname = Iname - 1

 End Subroutine TreeMassesLam3

 Subroutine TreeMassesMSSM(gp, g, vevSM, M1, M2, M3, mu, B, tanb              &
                         &, M2_E, M2_L, A_l, Y_l                              &
                         &, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u              &
                         &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N         &
                         &, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2      &
                         &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2    &
                         &, RSup, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2 &
                         &, RSpm                                              &
                         &, GenerationMixing, kont, Slopy)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the MSSM
 ! written by Werner Porod, 24.03.2001 
 ! 27.03.02: adding the logical variable slopy, which is useful in the
 !           context of Sugra models, because there the tree-level
 !           results for scalar particles can lead to to tachonic
 !           states whereas the 1-loop calculation leads to correct
 !           results
 !           Slopy = .True. -> take the modules of the scalar masses squared
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu, B
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), tanb
  Logical, Intent(in) :: GenerationMixing, Slopy
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(2), mC2(2), mN(4), mN2(4)     &
                      &  , mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6), mP0(2) &
                      &  , mP02(2), RP0(2,2), mS0(2), mS02(2), RS0(2,2)     &
                      &  , mSpm(2), mSpm2(2)
  Complex(dp), Intent(out) :: PhaseGlu, U(2,2), V(2,2), N(4,4), Rsneut(3,3) &
                         &  , RSlepton(6,6), RSdown(6,6), RSup(6,6), RSpm(2,2)

  Integer :: i1, ierr
  Real(dp) :: D_sneut, T3, Yl, Yr, Ml, Mr, msf(2), msf2(2), cosb, sinb
  Real(dp) :: mStop2(2), mSbottom2(2), mT, mB, test(2), Rsn(3,3)
  Complex(dp) :: A, Y, Rsf(2,2), Ap, mat3(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesMSSM'

  kont = 0
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, g, mC, U, V, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mC2 = mC**2
  Call NeutralinoMass(M1, M2, mu, vevSM, gp, g, mN, N, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 =  M2_L + D_sneut * id3C
   If (Maxval(Abs(Aimag(mat3))).Eq.0._dp) Then
    Call EigenSystem(Real(mat3,dp),mSneut2,Rsn,ierr, test)
    Rsneut = Rsn
   Else
    Call EigenSystem(mat3,mSneut2,Rsneut,ierr, test)
   End If

    If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     ierr = 0
    End If
  
    If (ierr.Ne.0) Then
     Write(ErrCan,*) 'Problems with the diagonalization of sneutrinos'
     Write(ErrCan,*) 'in routine ',NameOfUnit(Iname),'. ierr = ',ierr
     kont = ierr
     Iname = Iname - 1
     Return
    Else
     Do i1=1,3
      If (mSneut2(i1).Ge.0._dp) Then
       mSneut(i1) = Sqrt( mSneut2(i1) )
      Else
       If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname),' mSneut2 ',i1
        Write(ErrCan,*) '< 0, : ',mSneut2(i1),'is set to its modulus '
       End If
       mSneut2(i1) = Abs(msneut2(i1))
       mSneut(i1) = Sqrt( mSneut2(i1) )
       kont = -226
       Call AddError(226)
       If (ErrorLevel.Eq.2) Call TerminateProgram
      End If
     Enddo
    End If
    If ((kont.Ne.0).And.(.Not.slopy)) Then
     Iname = Iname - 1
     Return
    End If
   !--------------
   ! Sleptons
   !--------------
    T3 = -0.5_dp
    Yl = -1._dp 
    Yr = 2._dp
    Call SfermionMass(M2_L, M2_E, A_l, mu, vevSM, Y_l, T3, Yl, Yr, g,gp, kont &
                    &, mSlepton, mSlepton2, Rslepton) 
    If ((kont.Ne.0).And.(Slopy)) Then
     mSlepton2 = Abs(mSlepton2) + 1._dp
     mSlepton = Sqrt( mSlepton2 )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
    Yl = 1._dp / 3._dp
    Yr = 2._dp / 3._dp
    Call SfermionMass(M2_Q, M2_D, A_d, mu, vevSM, Y_d, T3, Yl, Yr, g,gp, kont &
                    &, mSdown, mSdown2, Rsdown) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSdown2 = Abs(mSdown2) + 1._dp
    mSdown = Sqrt( mSdown2 )
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Yl = 1._dp / 3._dp
    Yr = -4._dp / 3._dp
    Call SfermionMass(M2_Q, M2_U, A_u, mu, vevSM, Y_u, T3, Yl, Yr, g,gp, kont &
                    &, mSup, mSup2, Rsup) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSup2 = Abs(mSup2) + 1._dp
    mSup = Sqrt( mSup2 )
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -226
      Call AddError(226)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     If ((kont.Ne.0).And.(.Not.slopy)) Then
      Iname = Iname - 1
      Return
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = A_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSlepton(2*(i1-1)+1:2*i1) = msf
    mSlepton2(2*(i1-1)+1:2*i1) = msf2
    RSlepton(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSlepton2(2*(i1-1)+1:2*i1) = Abs(mSlepton2(2*(i1-1)+1:2*i1)) + 1._dp
     mSlepton(2*(i1-1)+1:2*i1) = Sqrt( mSlepton2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSdown2(2*(i1-1)+1:2*i1) = Abs(mSdown2(2*(i1-1)+1:2*i1)) + 1._dp
     mSdown(2*(i1-1)+1:2*i1) = Sqrt( mSdown2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSup2(2*(i1-1)+1:2*i1) = Abs(mSup2(2*(i1-1)+1:2*i1)) + 1._dp
     mSup(2*(i1-1)+1:2*i1) = Sqrt( mSup2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
  End If
  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   cosb = 1._dp / Sqrt(1._dp + tanb**2)
   sinb = cosb * tanb
   mP0(1) = mZ 
   mP02(1) = mZ**2
   mP02(2) = Real(B,dp) / (sinb * cosb)
   If (mP02(2).Gt.0._dp) Then
    mP0(2) = Sqrt(mP02(2) )
   Else
    If (.Not.Slopy) Then
     kont =-227
     Call AddError(227)
     Iname = Iname - 1
     Return
    End If
    mP0(2) = Sqrt(Abs(mP02(2)) )
   End If    
   RP0(1,1) = - cosb
   RP0(1,2) = sinb 
   RP0(2,1) = sinb 
   RP0(2,2) = cosb
  !---------------
  ! charged
  !---------------
   mSpm(1) = mW
   mSpm2(1) = mW**2
   mSpm2(2) = mSpm2(1) + mP02(2)
   If (mSpm2(2).Gt.0._dp) Then
    mSpm(2) = Sqrt(mSpm2(2) )
   Else
    If (.Not.Slopy) Then
     kont =-228
     Call AddError(228)
     Iname = Iname - 1
     Return
    End If
    mSpm(2) = Sqrt(Abs(mSpm2(2)) )
   End If    
   RSpm = RP0
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If ((kont.Ne.0).And.(slopy)) Then
    msf2 = Abs(msf2) + 1._dp
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   mStop2 = msf2
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_D(3,3)
   A = A_d(3,3)
   Y = Y_d(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If ((kont.Ne.0).And.(slopy)) Then
    msf2 = Abs(msf2) + 1._dp
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   mSbottom2 = msf2
  Else
   mStop2 = mSup2(5:6) 
   mSbottom2 = mSdown2(5:6) 
  End If
  mB = Abs(Y_d(3,3)) * vevSM(1) * oosqrt2
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3) / Y_u(3,3)
  Ap = A_d(3,3) / Y_d(3,3)
  Call ScalarMassMSSMeff(mP02(2), mZ2, tanb, vevSM               &
                       &, mT, mStop2, A, mu, mB, mSbottom2, Ap   &
                       &, mS0, mS02, RS0, kont)
  Iname = Iname - 1

 End Subroutine TreeMassesMSSM

 Subroutine TreeMassesMSSM2(gp, g, vevSM, M1, M2, M3, mu, B, tanb             &
                         &, M2_E, M2_L, A_l, Y_l                              &
                         &, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u              &
                         &, uU_L, uU_R ,uD_L, uD_R, uL_L, uL_R                &
                         &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N         &
                         &, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2      &
                         &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2    &
                         &, RSup, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2 &
                         &, RSpm, mZ2, mW2                                    &
                         &, GenerationMixing, kont, Slopy, loophiggs, mf_u_Q)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the MSSM
 ! written by Werner Porod, 24.03.2001 
 ! 10.10.01: including diagonalization of SM fermions
 ! 27.03.02: adding the logical variable slopy, which is useful in the
 !           context of Sugra models, because there the tree-level
 !           results for scalar particles can lead to to tachonic
 !           states whereas the 1-loop calculation leads to correct
 !           results
 !           Slopy = .True. -> take the modules of the scalar masses squared
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu, B
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), tanb, mZ2, mW2
  Logical, Intent(in) :: GenerationMixing, Slopy
  Logical, Intent(in), Optional ::  loophiggs
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(2), mC2(2), mN(4), mN2(4)     &
                      &  , mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6), mP0(2) &
                      &  , mP02(2), RP0(2,2), mS0(2), mS02(2), RS0(2,2)     &
                      &  , mSpm(2), mSpm2(2)
  Complex(dp), Intent(out) :: PhaseGlu, U(2,2), V(2,2), N(4,4), Rsneut(3,3) &
       &  , RSlepton(6,6), RSdown(6,6), RSup(6,6), RSpm(2,2)                &
       &  , uU_L(3,3), uU_R(3,3) ,uD_L(3,3), uD_R(3,3), uL_L(3,3), uL_R(3,3)
  Real(dp), Intent(out), Optional :: mf_u_Q(3)

  Integer :: i1, ierr
  Real(dp) :: D_sneut, T3, Yl, Yr, Ml, Mr, msf(2), msf2(2), cosb, sinb
  Real(dp) :: mStop2(2), mSbottom2(2), mT, mB, mf_t(3), test(2), mf_ta(3)
  Complex(dp) :: A, Y, Rsf(2,2), Ap, mat3(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesMSSM2'

  kont = 0
  !-------------------------------
  ! SM-fermions
  !-------------------------------
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_t,uL_L,uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_t, uD_L, uD_R &
                                     & , mf_ta, uU_L, uU_R)
   If (Present(mf_u_Q)) mf_u_Q = mf_ta
  Else
   uL_L = id3C
   uL_R = id3C
   uD_L = id3C
   uD_R = id3C
   uU_L = id3C
   uU_R = id3C
   If (Present(mf_u_Q)) Then
    Do i1=1,3
     mf_u_Q(i1) = oosqrt2 * vevSM(2) * Abs(Y_u(i1,i1))
    End Do
   End If
  End If
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, g, mC, U, V, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mC2 = mC**2
  Call NeutralinoMass(M1, M2, mu, vevSM, gp, g, mN, N, kont)
  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return
  End If
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 =  M2_L + D_sneut * id3C
    Call ComplexEigenSystem(mat3,mSneut2,Rsneut,ierr, test)

    If (ierr.Eq.-14) Then
      Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
      Write(ErrCan,*) "test =",test
      Write(ErrCan,*) " "
      If (ErrorLevel.Eq.2) Call TerminateProgram
      ierr = 0
    End If
  
    If (ierr.Ne.0) Then
     Write(ErrCan,*) 'Problems with the diagonalization of sneutrinos'
     Write(ErrCan,*) 'in routine ',NameOfUnit(Iname),'. ierr = ',ierr
     kont = ierr
     Iname = Iname - 1
     Return
    Else
     Do i1=1,3
      If (mSneut2(i1).Ge.0._dp) Then
       mSneut(i1) = Sqrt( mSneut2(i1) )
      Else
       kont = -229
      Call AddError(229)
       If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname),' mSneut2 ',i1
        Write(ErrCan,*) '< 0, : ',mSneut2(i1),'is set to 0 '
       End If
       mSneut2(i1) = 0
       mSneut(i1) = 0
       If (ErrorLevel.Eq.2) Call TerminateProgram
      End If
     Enddo
    End If
    If ((kont.Ne.0).And.(.Not.slopy)) Then
     Iname = Iname - 1
     Return
    End If
   !--------------
   ! Sleptons
   !--------------
    T3 = -0.5_dp
    Yl = -1._dp 
    Yr = 2._dp
    Call SfermionMass(M2_L, M2_E, A_l, mu, vevSM, Y_l, T3, Yl, Yr, g,gp, kont &
                    &, mSlepton, mSlepton2, Rslepton) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSlepton2 = Abs(mSlepton2) + 1._dp
    mSlepton = Sqrt( mSlepton2 )
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
    Yl = 1._dp / 3._dp
    Yr = 2._dp / 3._dp
    Call SfermionMass(M2_Q, M2_D, A_d, mu, vevSM, Y_d, T3, Yl, Yr, g,gp, kont &
                    &, mSdown, mSdown2, Rsdown) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSdown2 = Abs(mSdown2) + 1._dp
    mSdown = Sqrt( mSdown2 )
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Yl = 1._dp / 3._dp
    Yr = -4._dp / 3._dp
    Call SfermionMass(M2_Q, M2_U, A_u, mu, vevSM, Y_u, T3, Yl, Yr, g,gp, kont &
                    &, mSup, mSup2, Rsup) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSup2 = Abs(mSup2) + 1._dp
    mSup = Sqrt( mSup2 )
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -229
      Call AddError(229)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
    If ((kont.Ne.0).And.(.Not.slopy)) Then
     Iname = Iname - 1
     Return
    End If
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = A_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSlepton(2*(i1-1)+1:2*i1) = msf
    mSlepton2(2*(i1-1)+1:2*i1) = msf2
    RSlepton(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSlepton2(2*(i1-1)+1:2*i1) = Abs(mSlepton2(2*(i1-1)+1:2*i1)) + 1._dp
     mSlepton(2*(i1-1)+1:2*i1) = Sqrt( mSlepton2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSdown2(2*(i1-1)+1:2*i1) = Abs(mSdown2(2*(i1-1)+1:2*i1)) + 1._dp
     mSdown(2*(i1-1)+1:2*i1) = Sqrt( mSdown2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSup2(2*(i1-1)+1:2*i1) = Abs(mSup2(2*(i1-1)+1:2*i1)) + 1._dp
     mSup(2*(i1-1)+1:2*i1) = Sqrt( mSup2(2*(i1-1)+1:2*i1) )
    Else If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
   End Do
  End If
  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   cosb = Sqrt(1._dp / (1._dp + tanb**2))
   sinb = cosb * tanb
   mP02(1) = mZ2
   mP0(1) = Sqrt(mZ2)
   mP02(2) = Real(B,dp) / (sinb * cosb)
   If (mP02(2).Gt.0._dp) Then
    mP0(2) = Sqrt(mP02(2) )
   Else
    If (.Not.Slopy) Then
     kont =-230
     Call AddError(230)
     Iname = Iname - 1
     Return
    End If
    mP0(2) = Sqrt(Abs(mP02(2)) )
   End If    
   mP0(2) = Sqrt( Abs(mP02(2)) )
   RP0(1,1) = - cosb
   RP0(1,2) = sinb 
   RP0(2,1) = sinb 
   RP0(2,2) = cosb
  !---------------
  ! charged
  !---------------
   mSpm2(1) = mW2
   mSpm(1) = Sqrt(mW2)
   mSpm2(2) = mSpm2(1) + mP02(2)
   If (mSpm2(2).Gt.0._dp) Then
    mSpm(2) = Sqrt(mSpm2(2) )
   Else
    If (.Not.Slopy) Then
     kont =-231
     Call AddError(231)
     Iname = Iname - 1
     Return
    End If
    mSpm(2) = Sqrt(Abs(mSpm2(2)) )
   End If    
   mSpm(2) = Sqrt(Abs(mSpm2(2)))
   RSpm = RP0
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If ((kont.Ne.0).And.(slopy)) Then
    msf2 = Abs(msf2) + 1._dp
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   mStop2 = msf2
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_D(3,3)
   A = A_d(3,3)
   Y = Y_d(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If ((kont.Ne.0).And.(slopy)) Then
    msf2 = Abs(msf2) + 1._dp
   Else If (kont.Ne.0) Then
    Iname = Iname - 1
    Return
   End If
   mSbottom2 = msf2
  Else
   mStop2 = mSup2(5:6) 
   mSbottom2 = mSdown2(5:6) 
  End If
  mB = Abs(Y_d(3,3)) * vevSM(1) * oosqrt2
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3) / Y_u(3,3)
  Ap = A_d(3,3) / Y_d(3,3)
  If (Present(LoopHiggs)) Then
   Call ScalarMassMSSMeff(mP02(2), mZ2, tanb, vevSM              &
                       &, mT, mStop2, A, mu, mB, mSbottom2, Ap   &
                       &, mS0, mS02, RS0, kont, LoopHiggs)
  Else
   Call ScalarMassMSSMeff(mP02(2), mZ2, tanb, vevSM              &
                       &, mT, mStop2, A, mu, mB, mSbottom2, Ap   &
                       &, mS0, mS02, RS0, kont)
  End If

  Iname = Iname - 1

 End Subroutine TreeMassesMSSM2

 Subroutine TreeMassesMSSM3(gp, g, vevSM, M1, M2, M3, mu, tanb             &
                         &, M2_E, M2_L, A_l, Y_l                              &
                         &, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u              &
                         &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N         &
                         &, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2      &
                         &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2    &
                         &, RSup, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2 &
                         &, RSpm                                              &
                         &, GenerationMixing, kont, Slopy)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the MSSM
 ! written by Werner Porod, 24.03.2001 
 ! 10.10.01: including diagonalization of SM fermions
 ! 27.03.02: adding the logical variable slopy, which is useful in the
 !           context of Sugra models, because there the tree-level
 !           results for scalar particles can lead to to tachonic
 !           states whereas the 1-loop calculation leads to correct
 !           results
 !           Slopy = .True. -> take the modules of the scalar masses squared
 ! 25.10.03: takes m_A as input instead as output
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
  Real(dp), Intent(in) :: g, gp, vevSM(2), tanb, mP0(2), mP02(2)
  Logical, Intent(in) :: GenerationMixing, Slopy
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(2), mC2(2), mN(4), mN2(4)     &
                      &  , mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6) &
                      &  , RP0(2,2), mS0(2), mS02(2), RS0(2,2)     &
                      &  , mSpm(2), mSpm2(2)
  Complex(dp), Intent(out) :: PhaseGlu, U(2,2), V(2,2), N(4,4), Rsneut(3,3) &
       &  , RSlepton(6,6), RSdown(6,6), RSup(6,6), RSpm(2,2)  

  Integer :: i1, ierr
  Real(dp) :: D_sneut, T3, Yl, Yr, Ml, Mr, msf(2), msf2(2), cosb, sinb
  Real(dp) :: mStop2(2), mSbottom2(2), mT, mB, mf_t(3), test(2), mf_ta(3)
  Complex(dp) :: A, Y, Rsf(2,2), Ap, mat3(3,3)
  Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesMSSM3'

  kont = 0
  !-------------------------------
  ! SM-fermions
  !-------------------------------
  If (GenerationMixing) Then
   Call FermionMass(Y_l,vevSM(1),mf_t,uL_L,uL_R,kont)
   Call QuarkMasses_and_PhaseShifts(Y_d, Y_u, vevSM, mf_t, uD_L, uD_R &
                                     & , mf_ta, uU_L, uU_R)
  Else
   uL_L = id3C
   uL_R = id3C
   uD_L = id3C
   uD_R = id3C
   uU_L = id3C
   uU_R = id3C
  End If
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0.d0,1.d0) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  Call CharginoMass(M2, mu, vevSM, g, mC, U, V, kont)
  mC2 = mC**2
  Call NeutralinoMass(M1, M2, mu, vevSM, gp, g, mN, N, kont)
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 = M2_L + D_sneut * id3C
    Call ComplexEigenSystem(mat3,mSneut2,Rsneut,ierr, test)

    If (ierr.Eq.-14) Then
      Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
      Write(ErrCan,*) "test =",test
      Write(ErrCan,*) " "
      If (ErrorLevel.Eq.2) Call TerminateProgram
      ierr = 0
    End If
  
    If (ierr.Ne.0) Then
     Write(ErrCan,*) 'Problems with the diagonalization of sneutrinos'
     Write(ErrCan,*) 'in routine ',NameOfUnit(Iname),'. ierr = ',ierr
     kont = ierr
    Else
     Do i1=1,3
      If (mSneut2(i1).Ge.0._dp) Then
       mSneut(i1) = Sqrt( mSneut2(i1) )
      Else
       kont = -232
      Call AddError(232)
       If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Warning from ',NameOfUnit(Iname),' mSneut2 ',i1
        Write(ErrCan,*) '< 0, : ',mSneut2(i1),'is set to 0 '
       End If
       mSneut2(i1) = 0
       mSneut(i1) = 0
       If (ErrorLevel.Eq.2) Call TerminateProgram
      End If
     Enddo
    End If
   !--------------
   ! Sleptons
   !--------------
    T3 = -0.5_dp
    Yl = -1._dp 
    Yr = 2._dp
    Call SfermionMass(M2_L, M2_E, A_l, mu, vevSM, Y_l, T3, Yl, Yr, g,gp, kont &
                    &, mSlepton, mSlepton2, Rslepton) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSlepton2 = Abs(mSlepton2) + 1._dp
    mSlepton = Sqrt( mSlepton2 )
   End If
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
    Yl = 1._dp / 3._dp
    Yr = 2._dp / 3._dp
    Call SfermionMass(M2_Q, M2_D, A_d, mu, vevSM, Y_d, T3, Yl, Yr, g,gp, kont &
                    &, mSdown, mSdown2, Rsdown) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSdown2 = Abs(mSdown2) + 1._dp
    mSdown = Sqrt( mSdown2 )
   End If
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Yl = 1._dp / 3._dp
    Yr = -4._dp / 3._dp
    Call SfermionMass(M2_Q, M2_U, A_u, mu, vevSM, Y_u, T3, Yl, Yr, g,gp, kont &
                    &, mSup, mSup2, Rsup) 
   If ((kont.Ne.0).And.(Slopy)) Then
    mSup2 = Abs(mSup2) + 1._dp
    mSup = Sqrt( mSup2 )
   End If
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -232
      Call AddError(232)
      If (ErrorLevel.Ge.0) Then
        Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write(ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write(ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = A_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSlepton(2*(i1-1)+1:2*i1) = msf
    mSlepton2(2*(i1-1)+1:2*i1) = msf2
    RSlepton(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSlepton2(2*(i1-1)+1:2*i1) = Abs(mSlepton2(2*(i1-1)+1:2*i1)) + 1._dp
     mSlepton(2*(i1-1)+1:2*i1) = Sqrt( mSlepton2(2*(i1-1)+1:2*i1) )
    End If
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSdown2(2*(i1-1)+1:2*i1) = Abs(mSdown2(2*(i1-1)+1:2*i1)) + 1._dp
     mSdown(2*(i1-1)+1:2*i1) = Sqrt( mSdown2(2*(i1-1)+1:2*i1) )
    End If
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
    If ((kont.Ne.0).And.(Slopy)) Then
     mSup2(2*(i1-1)+1:2*i1) = Abs(mSup2(2*(i1-1)+1:2*i1)) + 1._dp
     mSup(2*(i1-1)+1:2*i1) = Sqrt( mSup2(2*(i1-1)+1:2*i1) )
    End If
   End Do
  End If
  !---------------
  ! Higgs bosons
  !---------------
  ! pseudoscalar
  !---------------
   cosb = Sqrt(1._dp / (1._dp + tanb**2))
   sinb = cosb * tanb
   RP0(1,1) = - cosb
   RP0(1,2) = sinb 
   RP0(2,1) = sinb 
   RP0(2,2) = cosb
  !---------------
  ! charged
  !---------------
   mSpm(1) = mW
   mSpm2(1) = mW**2
   mSpm2(2) = mSpm2(1) + mP02(2)
   mSpm(2) = Sqrt(mSpm2(2))
   RSpm = RP0
  !---------------
  ! scalar
  !---------------
  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If (kont.Ne.0) msf2 = Abs(msf2) + 1._dp
   mStop2 = msf2
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_D(3,3)
   A = A_d(3,3)
   Y = Y_d(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, Rsf)
 ! tree level results are not always reliable in mSugra
   If (kont.Ne.0) msf2 = Abs(msf2) + 1._dp
   mSbottom2 = msf2
  Else
   mStop2 = mSup2(5:6) 
   mSbottom2 = mSdown2(5:6) 
  End If
  mB = Abs(Y_d(3,3)) * vevSM(1) * oosqrt2
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3) / Y_u(3,3)
  Ap = A_d(3,3) / Y_d(3,3)
  Call ScalarMassMSSMeff(mP02(2), mZ2, tanb, vevSM               &
                       &, mT, mStop2, A, mu, mB, mSbottom2, Ap   &
                       &, mS0, mS02, RS0, kont)

  Iname = Iname - 1

 End Subroutine TreeMassesMSSM3


 Subroutine TreeMassesNMSSM(gp, g, vevSM, vP, M1, M2, M3, mu_in, B, h0, lam   &
             & , Ah0, Alam, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u  &
             & , Y_d, Y_u, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N          &
             & , mSneut, mSneut2, Rsneut, mSlepton, mSlepton2, RSlepton       &
             & , mSdown, mSdown2, RSdown, mSup, mSup2, RSup, mP0, mP02, RP0   &
             & , mS0, mS02, RS0, mSpm, mSpm2, RSpm, GenerationMixing, kont    &
             & , higgs_loop)
 !-----------------------------------------------------------------
 ! calculates all SusyMasses in the NMSSM at tree level except for
 ! neutral Higgs bosons where the 1-loop effective potential is used
 ! written by Werner Porod, 04.03.05 
 !-----------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: M1, M2, M3, mu_in, B
  Complex(dp), Intent(in) :: h0, lam, Ah0, Alam
  Complex(dp), Intent(in) :: A_l(3,3), A_d(3,3), A_u(3,3)
  Complex(dp), Intent(in) :: Y_l(3,3), Y_d(3,3), Y_u(3,3)
  Complex(dp), Intent(in) :: M2_E(3,3), M2_L(3,3), M2_D(3,3), M2_U(3,3)& 
                        &  , M2_Q(3,3)
!  real(dp), intent(in) :: M2_H(2)
  Real(dp), Intent(in) :: g, gp, vevSM(2), vP
  Logical, Intent(in) :: GenerationMixing, higgs_loop
  Integer, Intent(inout) :: kont
  Real(dp), Intent(out) :: mGlu, mC(2), mC2(2), mN(5), mN2(5)     &
                      &  , mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
                      &  , mSdown(6), mSdown2(6), mSup(6), mSup2(6)
  Complex(dp), Intent(out) :: PhaseGlu, U(2,2), V(2,2), N(5,5), Rsneut(3,3) &
                         &  , RSlepton(6,6), RSdown(6,6), RSup(6,6)
  !----------------------------------------------------
  ! Higgs sector might be external, therefore inout
  !----------------------------------------------------
  Real(dp), Intent(inout) :: mP0(3), mP02(3), RP0(3,3), mS0(3), mS02(3)     &
                      &  , RS0(3,3), mSpm(2), mSpm2(2)
  Complex(dp), Intent(inout) :: RSpm(2,2)

  Integer :: i1, ierr
  Real(dp) :: D_sneut, T3, Yl, Yr, Ml, Mr, msf(2), msf2(2)
  Real(dp) :: mStop2(2), mSbottom2(2), mT, mB, test(2)
  Complex(dp) :: A, Y, Rsf(2,2), Ap, mat3(3,3), mu, Rstop(2,2), RSbot(2,2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'TreeMassesNMSSM'

  kont = 0
  !-------------------------------
  ! gluino
  !-------------------------------
  mGlu = Abs( M3 )
  PhaseGlu = Exp( (0._dp,1._dp) * Arg(M3) )
  !-------------------------------
  ! charginos + neutralinos
  !-------------------------------
  mu = mu_in + oosqrt2 * vP * h0
  Call CharginoMass(M2, mu, vevSM, g, mC, U, V, kont)
  mC2 = mC**2
  Call NeutralinoMass(M1, M2, mu_in, gp, g, h0, lam, vevSM, vP, mN, N, kont)
  mN2 = mN**2
  !-----------
  ! sfermions
  !-----------
  If (GenerationMixing) Then
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    mat3 =  M2_L + D_sneut * id3C
    Call ComplexEigenSystem(mat3,mSneut2,Rsneut,ierr, test)

    If (ierr.Eq.-14) Then
     Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
     Write(ErrCan,*) "test =",test
     Write(ErrCan,*) " "
     If (ErrorLevel.Eq.2) Call TerminateProgram
     ierr = 0
    End If
  
    If (ierr.Ne.0) Then
     Write (ErrCan,*) 'Problems with the diagonalization of sneutrinos'
     Write (ErrCan,*) 'in routine ',NameOfUnit(Iname),'. ierr = ',ierr
     kont = ierr
    Else
     Do i1=1,3
      If (mSneut2(i1).Ge.0._dp) Then
       mSneut(i1) = Sqrt( mSneut2(i1) )
      Else
       If (ErrorLevel.Ge.0) Then
        Write (ErrCan,*) 'Warning from ',NameOfUnit(Iname),' mSneut2 ',i1
        Write (ErrCan,*) '< 0, : ',mSneut2(i1),'is set to its modulus '
       End If
       mSneut2(i1) = abs(msneut2(i1))
       mSneut(i1) = Sqrt( mSneut2(i1) )
       kont = -233
       Call AddError(233)
       If (ErrorLevel.Eq.2) Call TerminateProgram
      End If
     Enddo
    End If
   !--------------
   ! Sleptons
   !--------------
    T3 = -0.5_dp
    Yl = -1._dp 
    Yr = 2._dp
    Call SfermionMass(M2_L, M2_E, A_l, mu, vevSM, Y_l, T3, Yl, Yr, g,gp, kont &
                    &, mSlepton, mSlepton2, Rslepton) 
   !--------------
   ! D-squarks
   !--------------
    T3 = -0.5_dp
    Yl = 1._dp / 3._dp
    Yr = 2._dp / 3._dp
    Call SfermionMass(M2_Q, M2_D, A_d, mu, vevSM, Y_d, T3, Yl, Yr, g,gp, kont &
                    &, mSdown, mSdown2, Rsdown) 
   !--------------
   ! U-squarks
   !--------------
    T3 = 0.5_dp
    Yl = 1._dp / 3._dp
    Yr = -4._dp / 3._dp
    Call SfermionMass(M2_Q, M2_U, A_u, mu, vevSM, Y_u, T3, Yl, Yr, g,gp, kont &
                    &, mSup, mSup2, Rsup) 
  Else
   !-----------
   ! Sneutrino
   !-----------
    D_sneut = 0.125_dp * (g**2 + gp**2) * (vevSM(1)**2 - vevSM(2)**2)
    Do i1=1,3
     mSneut2(i1) = Real(M2_L(i1,i1),dp)  + D_sneut
     If (mSneut2(i1).Lt.0._dp) Then
      kont = -233
       Call AddError(233)
      If (ErrorLevel.Ge.0) Then
        Write (ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
        Write (ErrCan,*) 'mSneutrino^2 ',i1,' <= 0 :',mSneut2(i1)
        Write (ErrCan,*) 'setting it to 10.'
        mSneut2(i1) = 10._dp
      End If
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If
     mSneut(i1) = Sqrt(mSneut2(i1))
    End Do
    RSneut = id3C
   !----------------
   ! sleptons
   !---------------
   RSlepton = 0
   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
   Do i1 = 1,3
    Ml = M2_L(i1,i1)
    Mr = M2_E(i1,i1)
    A = A_l(i1,i1)
    Y = Y_l(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSlepton(2*(i1-1)+1:2*i1) = msf
    mSlepton2(2*(i1-1)+1:2*i1) = msf2
    RSlepton(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! down-squarks
   !---------------
   RSdown = 0
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_D(i1,i1)
    A = A_d(i1,i1)
    Y = Y_d(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSdown(2*(i1-1)+1:2*i1) = msf
    mSdown2(2*(i1-1)+1:2*i1) = msf2
    RSdown(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
   !---------------
   ! up-squarks
   !---------------
   RSup = 0
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Do i1 = 1,3
    Ml = M2_Q(i1,i1)
    Mr = M2_U(i1,i1)
    A = A_u(i1,i1)
    Y = Y_u(i1,i1)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    mSup(2*(i1-1)+1:2*i1) = msf
    mSup2(2*(i1-1)+1:2*i1) = msf2
    RSup(2*(i1-1)+1:2*i1,2*(i1-1)+1:2*i1) = Rsf
   End Do
  End If

  If (external_higgs) Return

  If (GenerationMixing) Then
   !----------------------------------------------------------------------
   ! I need the stops and sbottoms for the 1-loop calculation for neutral
   ! Higgs bosons
   !----------------------------------------------------------------------
   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_U(3,3)
   A = A_u(3,3)
   Y = Y_u(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, RStop)
 ! tree level results are not always reliable in mSugra
   If (kont.Ne.0) msf2 = Abs(msf2) + 1._dp
   mStop2 = msf2
   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2_Q(3,3)
   Mr = M2_D(3,3)
   A = A_d(3,3)
   Y = Y_d(3,3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                     &, msf, msf2, RSbot)
 ! tree level results are not always reliable in mSugra
   If (kont.Ne.0) msf2 = Abs(msf2) + 1._dp
   mSbottom2 = msf2
  Else
   mStop2 = mSup2(5:6) 
   RStop = RSup(5:6,5:6)
   mSbottom2 = mSdown2(5:6) 
   RSbot = RSdown(5:6,5:6)
  End If
  mB = Abs(Y_d(3,3)) * vevSM(1) * oosqrt2
  mT = Abs(Y_u(3,3)) * vevSM(2) * oosqrt2
  A = A_u(3,3) / Y_u(3,3)
  Ap = A_d(3,3) / Y_d(3,3)

  Call ScalarMassNMSSMeff(Real(h0,dp), Real(lam,dp)  &
          & , Real(mu_in,dp), vevSM, vP, Real(Ah0,dp)    &
          & , Real(Alam,dp), mB, mT, mStop2, RStop, mSbottom2, RSbot          &
          & , Real(M2_L(1,1),dp), Real(M2_E(1,1),dp), Real(M2_Q(1,1),dp)       &
          & , Real(M2_U(1,1),dp), Real(M2_D(1,1),dp), Real(M2_L(3,3),dp)       &
          & , Real(M2_E(3,3),dp), Real(M2_Q(3,3),dp), Real(M2_U(3,3),dp)       &
          & , Real(M2_D(3,3),dp), Real(A,dp), Real(Ap,dp), Real(M2,dp)         &
          & , mS0, mS02, RS0, mP0, mP02, RP0, mSpm, mSpm2, RSpm, kont)

  Iname = Iname - 1

 End Subroutine TreeMassesNMSSM

End Module SusyMasses

