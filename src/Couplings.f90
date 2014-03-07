Module Couplings
! To do:
! 04.02.02: The self interaction of sfermions needs to be checked
!  if all particles are mass eigenstates

! load modules
Use Control
Use Mathematics, Only: CompareMatrices, Adjungate
Use StandardModel, Only: CKM
! load modules

! interfaces
 Interface AllCouplings
  Module Procedure AllCouplingsMSSM1G, AllCouplingsMSSM, AllCouplingsEps1  &
               & , AllCouplingsEps3, AllCouplingsLam, AllCouplingsNMSSM
 End Interface

 Interface CoupCSCharginoNeutralino
  Module Procedure CoupCSCharginoNeutralinoMSSM, CoupCSCharginoNeutralinoEps &
               & , CoupCSCharginoNeutralinoLam, CoupCSCharginoNeutralinoLam2 &
               & , CoupCSCharginoNeutralinoSpon, Coup_Chip_Chi0_Sm_NMSSM
 End Interface

 Interface CoupChargedScalarFermion
  Module Procedure CoupChargedScalarFermion1, CoupChargedScalarFermion3  &
         & , CoupChargedScalarFermion3Lam
 End Interface

 Interface CoupChargedScalarSfermion3
  Module Procedure CoupChargedScalarSfermion3MSSM1                       &
       & , CoupChargedScalarSfermion3Eps1, CoupChargedScalarSfermion3Eps3  &
       & , CoupChargedScalarSfermion3MSSM3
 End Interface
 
 Interface CoupChargedScalarSfermion4
  Module Procedure CoupChargedScalarSfermion4MSSM1                       &
       & , CoupChargedScalarSfermion4MSSM3
 End Interface
 
 Interface CoupCharginoPseudoScalar
  Module Procedure CoupCharginoPseudoScalarMSSM,   &
         &         CoupCharginoPseudoScalarEps1,   &
         &         CoupCharginoPseudoScalarEps3,   &
         &         CoupCharginoPseudoScalarLam3,   &
         &         CoupCharginoPseudoScalarSpon1,  &
         &         Coup_Chip_Chim_P0_NMSSM
 End Interface
 
 Interface CoupCharginoPseudoScalara
  Module Procedure CoupCharginoPseudoScalarMSSMa
 End Interface
 
 Interface CoupCharginoScalar
  Module Procedure CoupCharginoScalarMSSM,  CoupCharginoScalarEps1  &
    &   , CoupCharginoScalarEps3, CoupCharginoScalarLam3            &
    &   , CoupCharginoScalarSpon1, Coup_Chip_Chim_S0_NMSSM
 End Interface
 
 Interface CoupCharginoScalara
  Module Procedure CoupCharginoScalarMSSMa
 End Interface
 
 Interface CoupCharginoSfermion
  Module Procedure CoupCharginoSfermion1, CoupCharginoSfermion3   &
    &             , CoupCharginoSfermion3Lam 
 End Interface
 
 Interface CoupFermionPseudoScalar
  Module Procedure CoupFermionPseudoScalar1, CoupDquarkPseudoScalar1expl, &
                 & CoupDquarkPseudoScalar3expl, CoupFermionPseudoScalar3
 End Interface
 
 Interface CoupFermionScalar
  Module Procedure CoupFermionScalar1, CoupDquarkScalar1expl, &
                 & CoupDquarkScalar3expl, CoupFermionScalar3
 End Interface
 
 Interface CoupGluinoSquark
  Module Procedure CoupGluinoSquark1, CoupGluinoSquark3
 End Interface
 
 Interface CoupGravitinoSfermion
  Module Procedure CoupGravitinoSfermion1, CoupGravitinoSfermion3
 End Interface
 
 Interface CoupNeutralinoSfermion
  Module Procedure CoupNeutralinoSfermion1, CoupNeutralinoSfermion3 
 End Interface
 
 Interface CoupNeutralinoSdown
  Module Procedure CoupNeutralinoSdown1, CoupNeutralinoSdown3  &
               & , CoupNeutralinoSdown3Expl
 End Interface
 
 Interface CoupNeutralinoSlepton
  Module Procedure CoupNeutralinoSlepton1, CoupNeutralinoSlepton3 
 End Interface
 
 Interface CoupNeutralinoSneutrino
  Module Procedure CoupNeutralinoSneutrino1, CoupNeutralinoSneutrino3 
 End Interface
 
 Interface CoupNeutralinoSup
  Module Procedure CoupNeutralinoSup1, CoupNeutralinoSup3 
 End Interface
 
 Interface CoupPseudoScalarSfermion3
  Module Procedure CoupPseudoScalarSfermion3_1, CoupPseudoScalarSfermion3_3  &
                & , Coup_P0_Sf_Sf_NMSSM_1
 End Interface

 Interface CoupPseudoScalarSfermion3a
  Module Procedure CoupPseudoScalarSfermion3a_1, CoupPseudoScalarSfermion3a_3
 End Interface

 Interface CoupPseudoScalarSfermion4
  Module Procedure CoupPseudoScalarSfermion4_1, CoupPseudoScalarSfermion4_3
 End Interface

 Interface CoupScalarSfermion3
  Module Procedure CoupScalarSfermion3MSSM_1, CoupScalarSfermion3MSSM_3 &
                &, CoupScalarSfermion3Eps_1, Coup_S0_Sf_Sf_NMSSM_1      
 End Interface

 Interface CoupScalarSfermion3a
  Module Procedure CoupScalarSfermion3MSSMa_1, CoupScalarSfermion3MSSMa_3 &
                &, CoupScalarSfermion3Epsa_1
 End Interface

 Interface CoupScalarSfermion4
  Module Procedure CoupScalarSfermion4MSSM_1, CoupScalarSfermion4MSSM_3
 End Interface
 
 Interface CoupScalarW
  Module Procedure CoupScalarWMSSM, CoupScalarWRP
 End Interface
 
 Interface CoupScalarZ
  Module Procedure CoupScalarZMSSM, CoupScalarZrp
 End Interface

 Interface CoupSneutrinoZ
  Module Procedure CoupSneutrinoZ1, CoupSneutrinoZ3
 End Interface

 Interface CoupSfermion4Y
  Module Procedure CoupSfermion4Y_1, CoupSfermionSelf4Y_1, CoupSfermion4Y_3 &
                &, CoupSfermionSelf4Y_3
 End Interface

 Interface CoupSfermion4G
  Module Procedure CoupSfermion4G_13, CoupSfermionSelf4G_13
 End Interface
! interfaces
 
Contains


 Subroutine CoupChargedScalarFermion1(i, RSpm, yukD, yukU, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalars and fermions
 ! valid for the MSSM, 1-generation epsilon model, 
 ! 3-generation epsilon model, and the model where R-parity is broken 
 ! spontaneously
 !
 ! the lagrangian: H^-_i \bar(d) [ yuk_d R^*_i1 P_L +  yuk^*_u R_i2 P_L ] u
 !
 ! written by Werner Porod, 01.05.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i                 ! index of charged scalar
 Complex(dp), Intent(in) :: RSpm(:,:)     ! mixing matrix of charged scalars
 Complex(dp), Intent(in) ::  yukD, yukU   ! d-type and u-type Yukawa couplings
 Complex(dp), Intent(out) :: coupL, coupR ! left + right coupling defined above

 Integer :: n_Spm

 n_Spm = Size(RSpm, Dim=1)

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarFermion1'

 If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
  Write(ErrCan,*) 'Problem in Subroutine '//NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
  Call TerminateProgram
 End If

 coupL = yukD * Conjg( RSpm(i,1) )
 coupR = Conjg(yukU) * RSpm(i,2)

 Iname = Iname - 1

 End Subroutine CoupChargedScalarFermion1

 Subroutine CoupChargedScalarFermion3(i, j, k, RSpm, yukD, RdL, RdR, yukU &
                                    &, RuL, RuR, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalars and fermions
 ! valid for the 3-generation MSSM, and the
 ! 3-generation epsilon model, and the model where R-parity is broken 
 ! spontaneously
 !  i .......... index of charged scalar 
 !  j .......... index of down-type fermion
 !  k .......... index of up-type fermion
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  yukD(i,j) .. down type Yukawa coupling
 !  RdL(i,j) ... left mixing matrix for down type fermions
 !  RdR(i,j) ... right mixing matrix for down type fermions
 !  yukU(i,j) .. up type Yukawa coupling
 !  RuL(i,j) ... left mixing matrix for up type fermions
 !  RuR(i,j) ... right mixing matrix for up type fermions
 ! the lagrangian: H^-_i \bar(d_j) [ yuk_d R^*_i1 P_L +  yuk^*_u R_i2 P_R ] u_k
 ! output
 !  coupL ....... the left coupling(i,j,k)
 !  coupR ....... the right coupling(i,j,k)
 ! written by Werner Porod, 01.05.2001
 !-----------------------------------------------------------------------
 Implicit None

 Complex(dp), Intent(in) :: RSpm(:,:), yukD(3,3), RdL(3,3), RdR(3,3) &
            & , yukU(3,3), RuL(3,3), RuR(3,3)
 Complex(dp), Intent(out) :: coupL, coupR
 Integer, Intent(in) :: i, j, k

 Integer :: n_Spm, i1, i2, n_test

 n_Spm = Size(RSpm, Dim=1)

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarFermion3'

 If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
  Write(ErrCan,*) 'Problem in Subroutine '//NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
  Call TerminateProgram
 End If

 coupL = 0._dp
 coupR = 0._dp

 n_test = 0
 If (CompareMatrices(RdL,id3C,NearlyZero) ) n_test = n_test + 1
 If (CompareMatrices(RdR,id3C,NearlyZero) ) n_test = n_test + 10
 If (CompareMatrices(RuL,id3C,NearlyZero) ) n_test = n_test + 100
 If (CompareMatrices(RuR,id3C,NearlyZero) ) n_test = n_test + 1000

 Select Case(n_test)
 Case (1)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1)
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) )
   End Do
  End Do

 Case (10)
  Do i1=1,3
    coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) )
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2)
   End Do
  End Do

 Case (11)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1)
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) )
  End Do

 Case (100)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) )
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2)
   End Do
  End Do

 Case (101)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) )
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1)
  End Do

 Case (110)
  coupL = yukD(k,j) 
  Do i1=1,3
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2)
   End Do
  End Do

 Case (111)
  coupL = yukD(k,j)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1)
  End Do

 Case (1001)
  coupR = Conjg( yukU(j,k) )
  Do i1=1,3
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) )
   End Do
  End Do

 Case (1010)
  Do i1=1,3
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) )
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1)
  End Do

 Case (1011)
  coupR = Conjg( yukU(j,k) )
  Do i1=1,3
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) )
  End Do

 Case (1100)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) )
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1)
  End Do

 Case (1101)
  coupR = Conjg( yukU(j,k) )
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) )
  End Do

 Case (1110)
  coupL = yukD(k,j)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1)
  End Do

 Case (1111)
  coupL = yukD(k,j)
  coupR = Conjg( yukU(j,k) )

 Case Default
  Do i1=1,3
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) )
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2)
   End Do
  End Do

 End Select
 
 coupL = coupL * Conjg( RSpm(i,1) )
 coupR = coupR * RSpm(i,2)

 Iname = Iname - 1

 End Subroutine CoupChargedScalarFermion3

 Subroutine CoupChargedScalarPseudoscalarW(i, j, g, RSpm, RP0, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalar, pseudoscalars and
 ! the W-boson
 ! valid for the MSSM, NMSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of charged scalar 
 !  j .......... index of pseudoscalar 
 !  g .......... SU(2) gauge coupling
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 01.05.2001
 ! 14.12.03: correcting sign for MSSM case
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j
 Real(dp), Intent(in) :: g, RP0(:,:)
 Complex(dp), Intent(in) :: RSpm(:,:)
 Complex(dp), Intent(out) :: coup

 Integer :: n_Spm, n_P0

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarPseudoscalarW'

 n_P0 = Size( RP0, Dim=1)
 n_Spm = Size( RSpm, Dim=1)

 If (.Not.(   ((n_P0.Eq.2).And.(n_Spm.Eq.2))   & ! MSSM
    &     .Or.((n_P0.Eq.3).And.(n_Spm.Eq.2))   & ! NMSSM
    &     .Or.((n_P0.Eq.3).And.(n_Spm.Eq.4))   & ! 1-generation model
    &     .Or.((n_P0.Eq.5).And.(n_Spm.Eq.8))   & ! 3-generation model
    &     .Or.((n_P0.Eq.12).And.(n_Spm.Eq.8))  & ! full model
    &     .Or.((n_P0.Eq.6).And.(n_Spm.Eq.8))   & ! Munyoz Model
    &     .Or.((n_P0.Eq.7).And.(n_Spm.Eq.8))   & ! Munyoz Model 2 nuC
    &     ) ) Then 
  Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'model is not consistent, (n_P0,n_Spm) = ',n_P0,n_Spm
  Call TerminateProgram
 End If
 If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
  Call TerminateProgram
 End If
 If ((j.Lt.1).Or.(j.Gt.n_P0)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index j out of range: (j,n_P0) = ',j,n_P0
  Call TerminateProgram
 End If

 coup = 0._dp

 If ( (n_P0.Eq.2).And.(n_Spm.Eq.2) ) Then
  coup = Conjg( RSpm(i,1) ) * RP0(j,1) + Conjg( RSpm(i,2) ) * RP0(j,2)

 Else If ( (n_P0.Eq.3).And.(n_Spm.Eq.2) ) Then
  coup = Conjg( RSpm(i,1) ) * RP0(j,1) + Conjg( RSpm(i,2) ) * RP0(j,2)

 Else If ( (n_P0.Eq.3).And.(n_Spm.Eq.4) ) Then
  coup = Conjg( RSpm(i,1) ) * RP0(j,1)  &
    &  + Conjg( RSpm(i,2) ) * RP0(j,2)  &
    &  + Conjg( RSpm(i,3) ) * RP0(j,3)  

 Else If ( ( ((n_P0.ge.5).and.(n_P0.le.7)).Or.(n_P0.Eq.12)) .And. (n_Spm.Eq.8) )&
 & Then
  coup = Conjg( RSpm(i,1) ) * RP0(j,1)  &
    &  + Conjg( RSpm(i,2) ) * RP0(j,2)  &
    &  + Conjg( RSpm(i,3) ) * RP0(j,3)  &
    &  + Conjg( RSpm(i,4) ) * RP0(j,4)  &
    &  + Conjg( RSpm(i,5) ) * RP0(j,5)

 Else
  Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)//'  !'
  Write(ErrCan,*) 'The combination of ',n_Spm,' charged scalars'
  Write(ErrCan,*) 'and ',n_P0,' pseudoscalars is not possible.'
  Call TerminateProgram
 End If

 coup =  (0._dp,-0.5_dp) * g * coup

 Iname = Iname - 1

 End Subroutine CoupChargedScalarPseudoscalarW

 Subroutine CoupChargedScalarScalarW(i, j, g, RSpm, RS0, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalar, scalar bosons and
 ! the W-boson
 ! valid for the MSSM, NMSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of charged scalar 
 !  j .......... index of scalar boson
 !  g .......... SU(2) gauge coupling
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j
 Real(dp), Intent(in) :: g, RS0(:,:)
 Complex(dp), Intent(in) :: RSpm(:,:)
 Complex(dp), Intent(out) :: coup

 Integer :: n_Spm, n_S0

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarScalarW'

 n_S0 = Size( RS0, Dim=1)
 n_Spm = Size( RSpm, Dim=1)

 If (.Not.(   ((n_S0.Eq.2).And.(n_Spm.Eq.2))   & ! MSSM
    &     .Or.((n_S0.Eq.3).And.(n_Spm.Eq.2))   & ! NMSSM
    &     .Or.((n_S0.Eq.3).And.(n_Spm.Eq.4))   & ! 1-generation model
    &     .Or.((n_S0.Eq.5).And.(n_Spm.Eq.8))   & ! 3-generation model
    &     .Or.((n_S0.Eq.12).And.(n_Spm.Eq.8))  & ! full model
    &     .Or.((n_S0.Eq.6).And.(n_Spm.Eq.8))   & ! Munyoz model
    &     .Or.((n_S0.Eq.7).And.(n_Spm.Eq.8))   & ! Munyoz model 2 nuC
    &     ) ) Then
  Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'model is not consistent, (n_S0,n_Spm) = ',n_S0,n_Spm
  Call TerminateProgram
 End If
 If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
  Call TerminateProgram
 End If
 If ((j.Lt.1).Or.(j.Gt.n_S0)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index j out of range: (j,n_S0) = ',j,n_S0
  Call TerminateProgram
 End If

 coup = 0._dp

 If ( (n_S0.Eq.2).And.(n_Spm.Eq.2) ) Then
  coup = Conjg( RSpm(i,1) ) * RS0(j,1) - Conjg( RSpm(i,2) ) * RS0(j,2)

 Else If ( (n_S0.Eq.3).And.(n_Spm.Eq.2) ) Then
  coup = Conjg( RSpm(i,1) ) * RS0(j,1) - Conjg( RSpm(i,2) ) * RS0(j,2)

 Else If ( (n_S0.Eq.3).And.(n_Spm.Eq.4) ) Then
  coup = Conjg( RSpm(i,1) ) * RS0(j,1)  &
    &  - Conjg( RSpm(i,2) ) * RS0(j,2)  &
    &  + Conjg( RSpm(i,3) ) * RS0(j,3)  

 Else If ( ( ((n_S0.Ge.5).And.(n_S0.Le.7)).Or.(n_S0.Eq.12)) .And. (n_Spm.Eq.8) )&
 &Then
  coup = Conjg( RSpm(i,1) ) * RS0(j,1)  &
    &  - Conjg( RSpm(i,2) ) * RS0(j,2)  &
    &  + Conjg( RSpm(i,3) ) * RS0(j,3)  &
    &  + Conjg( RSpm(i,4) ) * RS0(j,4)  &
    &  + Conjg( RSpm(i,5) ) * RS0(j,5)

 Else
  Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)//'  !'
  Write(ErrCan,*) 'The combination of ',n_Spm,' charged scalars'
  Write(ErrCan,*) 'and ',n_S0,' scalars is not possible.'
  Call TerminateProgram
 End If

 coup = coup * g * 0.5_dp

 Iname = Iname - 1

 End Subroutine CoupChargedScalarScalarW

  Subroutine CoupChargedScalarZ(i, j, g, sinW2, RSpm, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalar bosons and the Z-boson
 ! valid for the MSSM, 1-generation epsilon model, 
 ! 3-generation epsilon model, and the spontaneously broken model
 !  i .......... index of pseudoscalar boson
 !  j .......... index of scalar boson
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 02.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(in) :: RSpm(:,:)
  Complex(dp), Intent(out) :: coup

  Integer :: i1, n_Spm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarZ'

  n_Spm = Size( RSpm, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.n_Spm)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range: (j,n_Spm) = ',j,n_Spm 
   Call TerminateProgram
  End If

  coup = 0._Dp

  Do i1=1,2
   coup = coup - RSpm(i,i1) * Conjg( RSpm(j,i1) )
  End Do
  If (n_Spm.Eq.4) Then
   coup = coup - RSpm(i,3) * Conjg( RSpm(j,3) )
  Else If (n_Spm.Eq.8) Then
   Do i1=3,5
    coup = coup - RSpm(i,i1) * Conjg( RSpm(j,i1) )
   End Do
  End If

  coup = coup * (0.5_dp - sinW2)

  If (n_Spm.Eq.4) Then
   coup = coup + RSpm(i,4) * Conjg( RSpm(j,4) ) * sinW2
  Else If (n_Spm.Eq.8) Then
   Do i1=6,8
    coup = coup + RSpm(i,i1) * Conjg( RSpm(j,i1) ) * sinW2
   End Do
  End If

  coup = coup * g / Sqrt(1._dp -sinW2)

 Iname = Iname - 1

 End Subroutine CoupChargedScalarZ

 Subroutine CoupCharginoNeutralinoW(i, j, N, U, V, g, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and the W-boson
 ! valid for the MSSM, NMSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model and the full R-parity model
 !  i .......... index of chargino
 !  j .......... index of the neutralino
 !  N.. ........ mixing matrix of the neutralino
 !  U,V ........ mixing matrices of the chargino
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! the lagrangian is given by
 !  W^\mu,- \bar{\chim(i)} \gamma_mu (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 5.8.1999
 ! - including variable boundaries for the matrices N,U,V
 !   simplifies life for different models, before the calculation is
 !   done, the range of the indices is checked: 23.3.2000
 ! 26.04.2001: porting to f90
 ! 09.11.2004: adding spontaneous model with one generation of singlets
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: g
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:)
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: n_char, n_neut

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoNeutralinoW'

  coupL = ZeroC
  coupR = ZeroC

  n_char = Size(U, Dim=1 )
  n_neut = Size(N, Dim=1 )
  !------
  ! MSSM
  !------
  If ((n_char.Eq.2).And.(n_neut.Eq.4)) Then
   If ((i.Lt.1).Or.(i.Gt.2)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.4)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

   coupL = - U(i,1) * Conjg( N(j,2) ) - U(i,2) * Conjg( N(j,3) ) * oosqrt2
 
  !------
  ! NMSSM
  !------
  Else If ((n_char.Eq.2).And.(n_neut.Eq.5)) Then
   If ((i.Lt.1).Or.(i.Gt.2)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.5)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

   coupL = - U(i,1) * Conjg( N(j,2) ) - U(i,2) * Conjg( N(j,3) ) * oosqrt2
 
  !----------------------------
  ! 1-generation epsilon model
  !----------------------------
  Else If ((n_char.Eq.3).And.(n_neut.Eq.5)) Then
   If ((i.Lt.1).Or.(i.Gt.3)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.5)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

   coupL = - U(i,1) * Conjg( N(j,2) )         & 
         & - ( U(i,2) * Conjg( N(j,3) ) + U(i,3) * Conjg( N(j,5) ) ) * oosqrt2

  !----------------------------
  ! 3-generation epsilon model
  !----------------------------
  Else If ((n_char.Eq.5).And.(n_neut.Eq.7)) Then
   If ((i.Lt.1).Or.(i.Gt.5)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.7)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

  coupL = - U(i,1) * Conjg( N(j,2) )                                & 
        & - ( U(i,2) * Conjg( N(j,3) ) + U(i,3) * Conjg( N(j,5) )   &
        &   + U(i,4) * Conjg( N(j,6) ) + U(i,5) * Conjg( N(j,7) ) ) * oosqrt2

  !----------------------------
  ! Munyoz model 1 NuC
  !----------------------------
  Else If ((n_char.Eq.5).And.(n_neut.Eq.8)) Then
   If ((i.Lt.1).Or.(i.Gt.5)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.8)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

  coupL = - U(i,1) * Conjg( N(j,2) )                                & 
        & - ( U(i,2) * Conjg( N(j,3) ) + U(i,3) * Conjg( N(j,5) )   &
        &   + U(i,4) * Conjg( N(j,6) ) + U(i,5) * Conjg( N(j,7) ) ) * oosqrt2

  !----------------------------
  ! Munyoz model 2 nuC
  !----------------------------
  Else If ((n_char.Eq.5).And.(n_neut.Eq.9)) Then
   If ((i.Lt.1).Or.(i.Gt.5)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.9)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

  coupL = - U(i,1) * Conjg( N(j,2) )                                & 
        & - ( U(i,2) * Conjg( N(j,3) ) + U(i,3) * Conjg( N(j,5) )   &
        &   + U(i,4) * Conjg( N(j,6) ) + U(i,5) * Conjg( N(j,7) ) ) * oosqrt2

  !------------
  ! full model
  !------------
  Else If ((n_char.Eq.5).And.((n_neut.Eq.14).or.(n_neut.eq.10))) Then
   If ((i.Lt.1).Or.(i.Gt.5)) Then
    Write(ErrCan,*) 'Chargino index out of range, i, n_char',i,n_char
    Call TerminateProgram
   End If
   If ((j.Lt.1).Or.(j.Gt.n_neut)) Then
    Write(ErrCan,*) 'Neutralino index out of range, j, n_neut',j,n_neut
    Call TerminateProgram
   End If

  coupL = - U(i,1) * Conjg( N(j,2) )                                & 
        & - ( U(i,2) * Conjg( N(j,3) ) + U(i,3) * Conjg( N(j,5) )   &
        &   + U(i,4) * Conjg( N(j,6) ) + U(i,5) * Conjg( N(j,7) ) ) * oosqrt2

  Else
   Write(ErrCan,*) 'Error in Subroutine CoupCharginoNeutralinoW.'
   Write(ErrCan,*) 'Model (n_char,n_neut) (',n_char,n_neut,')'
   Write(ErrCan,*) 'is not defined.'
   Call TerminateProgram
  End If

  coupR = - Conjg( V(i,1) ) * N(j,2) + Conjg( V(i,2) ) * N(j,4)  * oosqrt2 

  coupL = g * coupL
  coupR = g * coupR

  Iname = Iname - 1

 End Subroutine CoupCharginoNeutralinoW

 Subroutine CoupCharginoZa(i, j, U, V, g, cosW, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and the Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of chargino, first index is an electroweak eigenstate
 !  U,V ........ mixing matrices of the chargino
 !  g .......... SU(2) gauge coupling
 !  cosW ....... cos(weak mixing angle) 
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 28.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g, cosW
  Complex(dp), Intent(in) :: U(:,:), V(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j

  Integer :: n_char

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoZa'

  n_char = Size( U, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: (i,n_char) = ',i,n_char
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_char)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range: (j,n_char) = ',j,n_char
   Call TerminateProgram
  End If

  If (i.Eq.1) Then 
   coupL = g * cosW * Conjg( U(j,1) )
   coupR = g * cosW * V(j,1)
  Else
   coupL = g * (cosW - 0.5_dp / cosW)
   coupR = coupL * V(j,2)
   coupL = coupL * Conjg( U(j,2) )
  End If

  Iname = Iname - 1

 End Subroutine CoupCharginoZa


 Subroutine CoupCharginoZ(i, j, U, V, g, cosW, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and the Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of chargino
 !  U,V ........ mixing matrices of the chargino
 !  g .......... SU(2) gauge coupling
 !  cosW ....... cos(weak mixing angle) 
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 28.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g, cosW
  Complex(dp), Intent(in) :: U(:,:), V(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j

  Integer :: n_char, i1
  Real(dp) :: sinW2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoZ'

  n_char = Size( U, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: (i,n_char) = ',i,n_char
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_char)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range: (j,n_char) = ',j,n_char
   Call TerminateProgram
  End If
 
  coupL = U(i,1) * Conjg( U(j,1) )
  Do i1=2,n_char
   coupL = coupL + 0.5_dp * U(i,i1) * Conjg( U(j,i1) )
  End Do

  coupR = V(j,1) * Conjg( V(i,1) ) + 0.5_dp * V(j,2) * Conjg( V(i,2) )

  If (i.Eq.j) Then
   sinW2 = 1.0_dp - cosW**2
   coupL = coupL - sinW2
   coupR = coupR - sinW2
  End If

  coupL = coupL * g / cosW
  coupR = coupR * g / cosW

  Iname = Iname - 1

 End Subroutine CoupCharginoZ


 Subroutine CoupFermionZ(T3,e,g,sinW2,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between SM-fermions and Z
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  T3 ......... weak isospin of the fermion
 !  e .......... electric charge of the fermion
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 ! output
 !  coupL ...... the left coupling
 !  coupR ...... the right coupling
 ! written by Werner Porod, 6.8.1999
 ! 2.10.2000: changing to f90 form
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: T3,e,g,sinW2
  Real(dp), Intent(out) :: coupL,coupR
  Real(dp) :: oocosW
 
  oocosW = 1._dp / Sqrt(1._dp - sinW2)
  coupL = g * ( -T3 + e*sinW2) * oocosW
  coupR = g * e*sinW2 * oocosW

 End Subroutine CoupFermionZ


 Subroutine CoupNeutralinoPseudoScalara(i,j,k,N,RP0,gp,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and pseudocalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model, assuming that the first neutralino index i
 ! corresponds to a electroweak eigenstate
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RP0(i,j) ... mixing matrix of pseudocalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g,RP0(:,:),gp
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_P0, n_N
  Complex(dp) :: sumI, sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoPseudoScalara'

  n_P0 = Size(RP0, Dim=1)
  n_N = Size(N, Dim=1)

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.(    ((n_N.Eq.4).And.(n_P0.Eq.2)) & ! MSSM
     &      .Or.((n_N.Eq.5).And.(n_P0.Eq.3)) & ! 1-generation epsilon model
     &      .Or.((n_N.Eq.7).And.(n_P0.Eq.5)) & ! 3-generation epsilon model
     & ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with (n_N,n_P0) = ',n_N,n_P0
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_P0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If

  If (i.Gt.2) Then
   sumI = (0._dp,0.5_dp) * ( gp * N(j,1) - g * N(j,2) )
   sumI_j = Conjg( sumI )
  
  Else If (n_N.Eq.4) Then ! MSSM 
   sumI_j = Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)

  Else If (n_N.Eq.5) Then ! 1-generation epsilon model
   sumI_j = Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)  &
      &   + Conjg(N(j,5)) * RP0(k,3)

  Else If (n_N.Eq.7) Then ! 3-generation epsilon model
   sumI_j = Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)  &
      &   + Conjg(N(j,5)) * RP0(k,3) + Conjg(N(j,6)) * RP0(k,4)  &
      &   + Conjg(N(j,7)) * RP0(k,5)

  End If 

  Select Case(i)
  Case(1)
   coupL = (0._dp,-0.5_dp) * gp * sumI_j
  Case(2)
   coupL = (0._dp,0.5_dp) * g * sumI_j
  Case Default
   coupL = sumI_j * RP0(k,i-2)
   If (i.Eq.4) coupL = - coupL
  End Select

  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoPseudoScalara


 Subroutine CoupNeutralinoPseudoScalar(i,j,k,N,RP0,gp,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and pseudocalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RP0(i,j) ... mixing matrix of pseudocalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g,RP0(:,:),gp
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_P0, n_N
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoPseudoScalar'

  n_P0 = Size(RP0, Dim=1)
  n_N = Size(N, Dim=1)

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.(    ((n_N.Eq.4).And.(n_P0.Eq.2)) & ! MSSM
     &      .Or.((n_N.Eq.5).And.(n_P0.Eq.3)) & ! 1-generation epsilon model
     &      .Or.((n_N.Eq.7).And.(n_P0.Eq.5)) & ! 3-generation epsilon model
     & ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with (n_N,n_P0) = ',n_N,n_P0
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_P0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = (0._dp,0.5_dp) * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = (0._dp,0.5_dp) * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  If (n_N.Eq.4) Then ! MSSM 
   coupL = sumI_i * ( Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2) ) &
      &  + sumI_j * ( Conjg(N(i,3)) * RP0(k,1) - Conjg(N(i,4)) * RP0(k,2) )

  Else If (n_N.Eq.5) Then ! 1-generation epsilon model
   coupL = sumI_i * ( Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)  &
      &             + Conjg(N(j,5)) * RP0(k,3) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RP0(k,1) - Conjg(N(i,4)) * RP0(k,2)  &
      &             + Conjg( N(i,5) ) * RP0(k,3) )

  Else If (n_N.Eq.7) Then ! 3-generation epsilon model
   coupL = sumI_i * ( Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)  &
      &             + Conjg(N(j,5)) * RP0(k,3) + Conjg(N(j,6)) * RP0(k,4)  &
      &             + Conjg(N(j,7)) * RP0(k,5) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RP0(k,1) - Conjg(N(i,4)) * RP0(k,2)  &
      &             + Conjg(N(i,5)) * RP0(k,3) + Conjg(N(i,6)) * RP0(k,4)  &
      &             + Conjg(N(i,7)) * RP0(k,5) )

  End If 
  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoPseudoScalar


 Subroutine CoupNeutralinoScalara(i,j,k,N,RS0,gp,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model, assuming that the first neutralino index i
 ! corresponds to a electroweak eigenstate
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RS0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 27.12.01
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,RS0(:,:),gp
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_S0, n_N
  Complex(dp) :: sumI, sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoScalara'

  n_S0 = Size(RS0, Dim=1)
  n_N = Size(N, Dim=1)

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.(    ((n_N.Eq.4).And.(n_S0.Eq.2)) & ! MSSM
     &      .Or.((n_N.Eq.5).And.(n_S0.Eq.3)) & ! 1-generation epsilon model
     &      .Or.((n_N.Eq.7).And.(n_S0.Eq.5)) & ! 3-generation epsilon model
     & ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with (n_N,n_S0) = ',n_N,n_S0
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_S0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  If (i.Gt.2) Then
   sumI = 0.5_dp * ( gp * N(j,1) - g * N(j,2) )
   sumI_j = Conjg( sumI )
  
  Else If (n_N.Eq.4) Then ! MSSM 
   sumI_j = Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)

  Else If (n_N.Eq.5) Then ! 1-generation epsilon model
   sumI_j = Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)  &
      &   + Conjg(N(j,5)) * RS0(k,3)

  Else If (n_N.Eq.7) Then ! 3-generation epsilon model
   sumI_j = Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)  &
      &   + Conjg(N(j,5)) * RS0(k,3) + Conjg(N(j,6)) * RS0(k,4)  &
      &   + Conjg(N(j,7)) * RS0(k,5)

  End If 

  Select Case(i)
  Case(1)
   coupL = 0.5_dp * gp * sumI_j
  Case(2)
   coupL = -0.5_dp * g * sumI_j
  Case Default
   coupL = sumI_j * RS0(k,i-2)
   If (i.Eq.4) coupL = - coupL
  End Select

  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoScalara


Subroutine CoupNeutralinoScalar(i,j,k,N,RS0,gp,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RS0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 6.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,RS0(:,:),gp
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_S0, n_N
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoScalar'

  n_S0 = Size(RS0, Dim=1)
  n_N = Size(N, Dim=1)

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.(    ((n_N.Eq.4).And.(n_S0.Eq.2)) & ! MSSM
     &      .Or.((n_N.Eq.5).And.(n_S0.Eq.3)) & ! 1-generation epsilon model
     &      .Or.((n_N.Eq.7).And.(n_S0.Eq.5)) & ! 3-generation epsilon model
     & ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with (n_N,n_S0) = ',n_N,n_S0
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_S0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = 0.5_dp * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = 0.5_dp * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  If (n_N.Eq.4) Then ! MSSM
   coupL = sumI_i * ( Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2) ) &
      &  + sumI_j * ( Conjg(N(i,3)) * RS0(k,1) - Conjg(N(i,4)) * RS0(k,2) )

  Else If (n_N.Eq.5) Then ! 1-generation epsilon model
   coupL = sumI_i * ( Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)  &
      &             + Conjg(N(j,5)) * RS0(k,3) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RS0(k,1) - Conjg(N(i,4)) * RS0(k,2)  &
      &             + Conjg( N(i,5) ) * RS0(k,3) )

  Else If (n_N.Eq.7) Then ! 3-generation epsilon model
   coupL = sumI_i * ( Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)  &
      &             + Conjg(N(j,5)) * RS0(k,3) + Conjg(N(j,6)) * RS0(k,4)  &
      &             + Conjg(N(j,7)) * RS0(k,5) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RS0(k,1) - Conjg(N(i,4)) * RS0(k,2)  &
      &             + Conjg(N(i,5)) * RS0(k,3) + Conjg(N(i,6)) * RS0(k,4)  &
      &             + Conjg(N(i,7)) * RS0(k,5) )

  End If 

  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoScalar


 Subroutine CoupNeutralinoZa(i,j,N,g,cosW,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model, assuming that the first neutralino index i
 ! corresponds to a electroweak eigenstate
 ! input
 !  i,j ........ index of Neutralino
 !  N .......... mixing matrix of the Neutralino
 !  g .......... SU(2) gauge coupling
 !  cosW ....... cos(theta_Weinberg)
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 28.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, cosW
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j

  Integer :: n_N

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoZa'

  n_N = Size(N, Dim=1)

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.( (n_N.Eq.4).Or.(n_N.Eq.5).Or.(n_N.Eq.7).Or.(n_N.Eq.14) ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with n_N = ',n_N
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  End If

  If ((i.Ge.3).And.(i.Le.7)) Then
   coupL =  - 0.5_dp * g * Conjg( N(j,i) ) / cosW
   If (i.Eq.4) coupL = - coupL
  Else
   coupL = 0._dp
  End If

  coupR = - Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoZa

 Subroutine CoupNeutralinoZ(i,j,N,g,cosW,coupL,coupR, NMSSM)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  N .......... mixing matrix of the Neutralino
 !  g .......... SU(2) gauge coupling
 !  cosW ....... cos(theta_Weinberg)
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 28.11.2000
 ! 09.11.2004: adding spontaneous model with one generation of singlets
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, cosW
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j
  Logical, Optional, intent(in) :: NMSSM

  Integer :: n_N
  logical :: l_nmssm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoZ'

  n_N = Size(N, Dim=1)

  l_nmssm = .False.
  if (present(NMSSM)) l_nmssm = NMSSM

  !---------------------------
  ! Checking first for model
  !---------------------------
  If (.Not.( (n_N.Eq.4).Or.(n_N.Eq.5).Or.((n_N.Ge.7).and.(n_N.Le.10)) &
           .Or.(n_N.Eq.14) ) ) Then
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'Model with n_N = ',n_N
    Write(ErrCan,*) 'not yet included. Setting coupling to 0.'
    coupL = 0
    coupR = 0
   End If
   If (ErrorLevel.Ge.1) Then
    Call TerminateProgram
   Else
    Iname = Iname - 1
    Return
   End If
  End If
  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  End If
 
  If (n_N.Eq.4) Then ! MSSM
   coupL = N(i,4) * Conjg( N(j,4) ) - N(i,3) * Conjg( N(j,3) )

  Else If (n_N.Eq.5) Then ! 1-generation epsilon model
   if (l_nmssm) then
    coupL = N(i,4) * Conjg( N(j,4) ) - N(i,3) * Conjg( N(j,3) )
   else
    coupL = N(i,4) * Conjg( N(j,4) ) - N(i,3) * Conjg( N(j,3) ) &
        & - N(i,5) * Conjg( N(j,5) )
   end if
  ! 3-generation epsilon model or full model, or Munyoz model
  Else If ( ((n_N.Ge.7).and.(n_N.Le.10)).Or.(n_N.Eq.14)) Then
   coupL = N(i,4) * Conjg( N(j,4) ) - N(i,3) * Conjg( N(j,3) )  &
       & - N(i,5) * Conjg( N(j,5) ) - N(i,6) * Conjg( N(j,6) )  &
       & - N(i,7) * Conjg( N(j,7) )

  End If 

  coupL = 0.5_dp * g * coupL / cosW
  coupR = - Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoZ

 Subroutine CoupPseudoscalarScalarZ(i,j,g,cosW,RP0,RS0,coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between pseudoscalar boson, scalar boson, and
 ! Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model and Munyoz model
 !  i .......... index of pseudoscalar boson
 !  j .......... index of scalar boson
 !  g .......... SU(2) gauge coupling
 !  cosW ....... cos(weak mixing angle)
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 4.8.1999
 ! 24.10.2000: porting tp f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g, cosW, RP0(:,:), RS0(:,:)
  Complex(dp), Intent(out) :: coup
  Integer, Intent(in) :: i,j

  Integer :: n_S0

  n_S0 = Size(RS0, dim=1)

  If (n_S0.Eq.2) Then
   coup = RP0(i,1)*RS0(j,1) - RP0(i,2)*RS0(j,2) 

  Else If (n_S0.Eq.3) Then
   coup = RP0(i,1)*RS0(j,1) - RP0(i,2)*RS0(j,2) + RP0(i,3)*RS0(j,3)

  Else If (((n_S0.Ge.5).and.(n_S0.Le.7)).Or.(n_S0.Eq.12)) Then
   coup = RP0(i,1)*RS0(j,1) - RP0(i,2)*RS0(j,2)    &
   &    + RP0(i,3)*RS0(j,3) + RP0(i,4)*RS0(j,4)    &
   &    + RP0(i,5)*RS0(j,5)

  Else
   Write(ErrCan,*) 'Problem in Subroutine CoupPseudoscalarScalarZ'
   Write(ErrCan,*) 'Model with n_S0 = ',n_S0,' is not yet included '
   Write(ErrCan,*) 'Setting the coupling to 0 '
  End If

  coup = (0._dp,-0.5_dp) * g *coup  / cosW

 End Subroutine CoupPseudoscalarScalarZ



 Subroutine CoupSdownZ(i, j, g, sinW2, RSf, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sdowns and Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of sdown
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RSf(i,j) ... mixing matrix of sdown
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: Rsf(:,:)
  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(out) ::coup
  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSdownZ'

  T3 = -0.5_dp
  e = -1._dp/3._dp
  Call CoupSfermionZ(i, j, g, sinW2, e, T3, RSf, coup)

  Iname = Iname - 1

 End Subroutine CoupSdownZ


 Subroutine CoupSfermionW3(i, j, g, RSf, RSfp, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sfermions and a W-boson
 ! valid for the 1- and 3- generation MSSM
 !  i .......... index of sfermion with isospin -1/2
 !  j .......... index of sfermion' with isospin 1/2
 !  g .......... SU(2) gauge coupling
 !  RSf(i,j) ... mixing matrix of sfermions
 !  RSfp(i,j) .. mixing matrix of sfermions'
 ! output
 !  coup ....... the coupling(i,j)
 ! the lagrangian is given by
 ! [I \tilde{d}^*_,\mu(i) \tilde{u}(j) - I \tilde{d}^*(i) \tilde{u}_,\mu(j) ]
 ! * W^-,\mu * coup
 ! written by Werner Porod, 23.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: Rsf(:,:), Rsfp(:,:)
  Real(dp), Intent(in) :: g
  Complex(dp), Intent(out) :: coup

  Integer :: i1, n_sfer, n_sferp

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSfermionW3'

  coup = 0._dp

  n_sfer = Size( Rsf, Dim=1)
  n_sferp = Size( Rsfp, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.n_Sfer)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index i = ',i
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_sferp)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index j = ',j
   Call TerminateProgram
  End If

  If (n_Sfer.Eq.2) Then
   coup = g * Rsf(i,1) * Conjg( Rsfp(j,1) ) * oosqrt2

  Else
   coup = 0._dp
   Do i1 = 1,3  
    coup = coup + g * Rsf(i,i1) * Conjg( Rsfp(j,i1) ) * oosqrt2
   End Do

  End If

  Iname = Iname - 1

 End Subroutine CoupSfermionW3


 Subroutine CoupSfermionZa(i, j, g, sinW2, e, T3, RSf, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sfermions and Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model 
 !  i .......... index of sfermion being an electroweak eigenstate
 !  j .-........ index of sfermion being a mass eigenstate
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  e .......... charge of the sfermion
 !  T3 ......... isospin of the left sfermion
 !  RSf(i,j) ... mixing matrix of sfermions
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 ! 23.04.2001: including MSSM, 3-generation case
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Complex(dp), Intent(in) :: Rsf(:,:)
  Real(dp), Intent(in) :: g,sinW2,e,T3
  Complex(dp), Intent(out) ::coup

  Integer :: n_sfer

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSfermionZa'

  coup = 0._dp

  n_sfer = Size( Rsf, Dim=1 )

  If ((i.Lt.1).Or.(i.Gt.n_sfer)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index i = ',i
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_sfer)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index j = ',j
   Call TerminateProgram
  End If

  If (n_sfer.Eq.2) Then
   If (i.Eq.1) Then 
    coup = (T3 - e * sinW2 ) * Rsf(j,1) 
   Else
    coup = - e * sinW2 * Rsf(j,2)
   End If

  Else If (n_sfer.Eq.3) Then ! sneutrinos
   coup = T3 * Rsf(j,i)

  Else If (n_sfer.Eq.6) Then 
   If (i.Le.3) Then
    coup = (T3 - e * sinW2 ) * Rsf(j,i)
   Else
    coup =  - e * sinW2 * Rsf(j,i)
   End If
  End If

  coup = coup * g / Sqrt(1._dp -sinW2)

  Iname = Iname - 1

 End Subroutine CoupSfermionZa

 Subroutine CoupSfermionZ(i, j, g, sinW2, e, T3, RSf, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sfermions and Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of sfermion
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  e .......... charge of the sfermion
 !  T3 ......... isospin of the left sfermion
 !  RSf(i,j) ... mixing matrix of sfermions
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 ! 23.04.2001: including MSSM, 3-generation case
 ! 15.11.04: simplifying the 3-gen. case
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Complex(dp), Intent(in) :: Rsf(:,:)
  Real(dp), Intent(in) :: g,sinW2,e,T3
  Complex(dp), Intent(out) ::coup

  Integer :: n_sfer, i1

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSfermionZ'

  coup = 0._dp

  n_sfer = Size( Rsf, Dim=1 )

  If ((i.Lt.1).Or.(i.Gt.n_sfer)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index i = ',i
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_sfer)) Then
   Write(ErrCan,*) 'Error in Soubroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Sfermion index j = ',j
   Call TerminateProgram
  End If

  If (n_sfer.Eq.2) Then  
   if (i.eq.j) coup = - e * sinW2
   coup = coup + T3 * Conjg( Rsf(i,1) ) * Rsf(j,1)

  Else If (n_sfer.Eq.3) Then ! sneutrinos
   if (i.eq.j) coup = T3

  Else If (n_sfer.Eq.6) Then 
   if (i.eq.j)  coup = - e * sinW2
   Do i1=1,3
    coup = coup + T3 * Conjg( Rsf(i,i1) ) * Rsf(j,i1)
   End Do
  End If

  coup = coup * g / Sqrt(1._dp -sinW2)

  Iname = Iname - 1

 End Subroutine CoupSfermionZ


 Subroutine CoupSleptonZ(i, j, g, sinW2, RSf, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sleptonss and Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of sleptons
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RSf(i,j) ... mixing matrix of sleptons
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: Rsf(:,:)
  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(out) ::coup

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSleptonZ'

  T3 = -0.5_dp
  e = -1._dp
  Call CoupSfermionZ(i, j, g, sinW2, e, T3, RSf, coup)

  Iname = Iname - 1

 End Subroutine CoupSleptonZ


 Subroutine CoupSneutrinoZ1(g, sinW2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sneutrino and Z-boson
 ! valid for the MSSM, 1-generation epsilon model
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 ! output
 !  coup ....... the coupling
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(out) :: coup

  coup = 0.5_dp * g / Sqrt(1._dp -sinW2)

 End Subroutine CoupSneutrinoZ1


 Subroutine CoupSneutrinoZ3(i, j, g, sinW2, Rsneut, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sneutrino and Z-boson
 ! valid for the MSSM, 1-generation epsilon model
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 ! output
 !  coup ....... the coupling
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: Rsneut(:,:)
  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(out) :: coup

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSneutZ3'

  T3 = 0.5_dp
  e = 0._dp
  Call CoupSfermionZ(i, j, g, sinW2, e, T3, RSneut, coup)

  Iname = Iname - 1

 End Subroutine CoupSneutrinoZ3



 Subroutine CoupSupZ(i, j, g, sinW2, RSf, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between sups and Z-boson
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of sup
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RSf(i,j) ... mixing matrix of sups
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 8.8.1999
 ! 2.10.2000: portation fo f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: Rsf(:,:)
  Real(dp), Intent(in) :: g, sinW2
  Complex(dp), Intent(out) :: coup

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupSupZ'

  T3 = 0.5_dp
  e = 2._dp/3._dp
  Call CoupSfermionZ(i, j, g, sinW2, e, T3, RSf, coup)

  Iname = Iname - 1

 End Subroutine CoupSupZ


 Subroutine AllCouplingsMSSM1G(g, Y_l, Y_d, Y_u, vevs                        &
    & , RSpm, RP0, RS0, Umat, Vmat, Nmat, mu, phase_glu                      &
    & , RSlepton, Ae, RSup, Au, RSdown, Ad                                   &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03       &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R       &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R    &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R   &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L           &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R           &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L         &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R           &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R           &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R           &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0         &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03 &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW    &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW        &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L    &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R             &
    & , cpl_SmpCN_L, cpl_SmpCN_R)
 !-----------------------------------------------------------------
 ! Routine for calculating all couplings of the MSSM
 ! output:
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 ! In addition, an _L or _R is added in the case of fermions indicated
 ! if it is the left-handed or right-handed coupling. 
 ! written by Werner Porod, 
 ! 01.03.2001: taking InitializeDecaysMSSM as basis
 !          - charged Higgs -sfermion-sfermion, 1 Gen., MSSM
 ! 12.09.2002: changing interface to avoid global variables
 !-------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Real(dp), Intent(in) :: vevS(2)     ! MSSM Higgs vevs [v_d, v_u]
  Complex(dp), Intent(in) :: RSpm(2,2) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(2,2)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(2,2)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(2,2), Vmat(2,2) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(4,4) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  mu   ! superpotential bilinear mu
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: RSlepton(6,6) ! slepton mixing matrix
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters

  Complex(dp), Intent(out) :: cpl_SmpSlSn(2,6,3), cpl_SmpSdSu(2,6,6)   &
      & , cpl_SmpSnSl(2,3,6), cpl_SmpSuSd(2,6,6), cpl_SmpP03(2,2,2)    &
      & , cpl_SmpP0W(2,2), cpl_SmpS03(2,2,2), cpl_SmpS0W(2,2)          &
      & , cpl_SmpLNu_L(2,3,3), cpl_SmpLNu_R(2,3,3), cpl_SmpDU_L(2,3,3) &
      & , cpl_SmpDU_R(8,3,3), cpl_SmpZ(2,2), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(2,2), cpl_CCZ_R(2,2)           &
      & , cpl_NNZ_L(4,4), cpl_NNZ_R(4,4), cpl_NNS0_L(4,4,2)            &
      & , cpl_NNS0_R(4,4,2), cpl_NNP0_L(4,4,2), cpl_NNP0_R(4,4,2) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,4,6), cpl_UNSu_R(3,4,6)         &
      & , cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6), cpl_NuNSn_L(3,4,3)      & 
      & , cpl_NuNSn_R(3,4,3), cpl_DDP0_L(3,3,2), cpl_LLP0_L(3,3,2)      &
      & , cpl_UUP0_L(3,3,2), cpl_DDP0_R(3,3,2), cpl_LLP0_R(3,3,2)       &
      & , cpl_UUP0_R(3,3,2), cpl_DDS0_L(3,3,2), cpl_LLS0_L(3,3,2)       &
      & , cpl_UUS0_L(3,3,2), cpl_DDS0_R(3,3,2), cpl_LLS0_R(3,3,2)       &
      & , cpl_UUS0_R(3,3,2)
  Complex(dp), Intent(out) :: cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6)      &
      & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)       &
      & , cpl_CLSn_R(2,3,3), cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(2)
  Complex(dp) :: cpl_P0SdSd(2,6,6), cpl_P0SuSu(2,6,6), cpl_P0SlSl(2,6,6) &
      & , cpl_P0SnSn(2,3,3), cpl_P0S0Z(2,2) 
  Real(dp), Intent(out) :: cpl_P0S03(2,2,2)
  Complex(dp), Intent(out) :: cpl_S0SdSd(2,6,6), cpl_S0SuSu(2,6,6) &
      & , cpl_S0SlSl(2,6,6), cpl_S0SnSn(2,3,3), cpl_LNuW(3,3)
  Real(dp), Intent(out) :: cpl_S03(2,2,2), cpl_S0WW(2), cpl_S0ZZ(2), cpl_FFpW
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(2,2,2), cpl_CCP0_R(2,2,2)    &
      & , cpl_CCS0_L(2,2,2), cpl_CCS0_R(2,2,2), cpl_CNW_L(2,4)        &
      & , cpl_CNW_R(2,4), cpl_SmpCN_L(2,2,4), cpl_SmpCN_R(2,2,4)

  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3, i4
  Real(dp) :: gU1, gSU2, gSU3, e_d, e_u, sinW2, cosW2, cosW
  Complex(dp) :: Rsd(2,2), Rsu(2,2), Rsl(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, bi(1)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsMSSM1G'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 2
  n_neut = 4
  n_P0 = 2
  n_S0 = 2
  n_Spm = 2

  !--------------------
  ! some constants
  !--------------------
  e_u = 2._dp / 3._dp
  e_d = - 1._dp / 3._dp

  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm     &
               &, gU1, gSU2, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpLNu_L = 0._dp
  cpl_SmpLNu_R = 0._dp
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp

   Do i1=1,2
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_l(i2,i2), ZeroC               &
                       &, cpl_SmpLNu_L(i1,i2,i2), cpl_SmpLNu_R(i1,i2,i2) )
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  !------------------------------
  ! charged scalar - pseudoscalar - W
  !------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_Spm
    Do i3 = 1,n_P0
     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
                                       &, cpl_SmpP03(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !--------------------------------------
  ! charged scalar - charged scalar - scalar
  !--------------------------------------
   Do i1 = 1, n_Spm
    Do i2 = 1, n_Spm
     Do i3 = 1, n_S0
      Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
                                   &, cpl_SmpS03(i1,i2,i3) )
    End Do
   End Do
  End Do
  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSlSn = ZeroC
  cpl_SmpSdSu = ZeroC

   Do i1 =1,n_Spm
    Do i2=1,3
     !----------
     ! Sleptons
     !----------
     Yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Call CoupChargedScalarSfermion3(i1, i3, 1, RSpm, gSU2, vevs, mu,    &
                       &             Yuk, A, Rsl, ZeroC, ZeroC, Id2C, coupC)
      cpl_SmpSlSn(i1,2*(i2-1)+i3,i2) = coupC
      cpl_SmpSnSl(i1,i2,2*(i2-1)+i3) = Conjg(coupC)
     End Do
     !----------
     ! Squarks
     !----------
     Yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Yukp = Y_u(i2,i2)
     Ap = Au(i2,i2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Do i4=1,2
       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
      End Do
     End Do
    End Do ! i2
   End Do ! i1

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, gSU2  &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, gSU2  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp
  cpl_CLSn_L = 0._dp
  cpl_CLSn_R = 0._dp
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    YukP = 0._dp
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSl, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CNuSl_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2C, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
     cpl_CLSn_L(i2, i1, i1) = coupLC
     cpl_CLSn_R(i2, i1, i1) = coupRC
    End Do

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_LLP0_L = 0._dp
  cpl_LLP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_d(i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_l(i1,i1), RP0   &
                                &, cpl_LLP0_L(i1,i1,i3), cpl_LLP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_LLS0_L = 0._dp
  cpl_LLS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, -0.5_dp, Y_d(i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, -0.5_dp, Y_l(i1,i1), RS0   &
                          &, cpl_LLS0_L(i1,i1,i3), cpl_LLS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  Call CoupFermionZ(-0.5_dp, -1._dp, gSU2, sinW2, cpl_LLZ_L, cpl_LLZ_R)
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)
  Call CoupFermionZ(0.5_dp, 0._dp, gSU2, sinW2, cpl_NuNuZ_L, cpl_NuNuZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  cpl_FFpW = - gSU2 * oosqrt2
  cpl_LNuW = id3C * cpl_FFpW

   cpl_DUW = id3C * cpl_FFpW
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_LNSl_L = 0._dp
  cpl_LNSl_R = 0._dp
  cpl_NuNSn_L = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

   Do i1 = 1,4
    Call CoupNeutralinoSneutrino(i1, gU1, gSU2, Nmat, coupRC)
    Do i2 = 1,3
     cpl_NuNSn_R(i2, i1, i2 ) = coupRC
    End Do
   End Do

   Do i1 = 1,3
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    Do i2 = 1,2
     Do i3 = 1,4
      Call CoupNeutralinoSlepton(i3, i2, gU1, gSU2, RSl, Yuk, Nmat &
                               &, coupLC, coupRC)
      cpl_LNSl_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_LNSl_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, Nmat, RP0, gU1, gSU2, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalar(i1, i2, i3, Nmat, RS0, gU1, gSU2, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  cpl_SlSnW = 0._dp

   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do

    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Call CoupSfermionW3(i2, 1, gSU2, RSl, id2C, coupC )
      cpl_SlSnW( (i1-1)*2+i2, i1) = coupC
    End Do
   End Do

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)
  Call Adjungate(cpl_SlSnW, cpl_SnSlW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
  Do i1=1,n_P0
   Do i2=1,n_P0
    Do i3=1,n_S0
     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
                                &, cpl_P0S03(i1,i2,i3) )
     End Do
    End Do
   End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
    End Do
   End Do
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp
  cpl_P0SlSl = 0._dp
  cpl_P0SnSn = 0._dp

  bi(1) = mu

   Do i1=1,n_P0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsl   &
                                   &, A, bi, coupC )
       cpl_P0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
                                   &, A, bi, coupC )
       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
                                   &, A, bi, coupC )
       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
  Do i1=1,n_S0
   Do i2=1,n_S0
    Do i3=1,n_S0
     Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
    End Do
   End Do
  End Do

  !-------------------
  ! scalar - W+ - W-
  !-------------------
  cpl_S0WW = 0._dp
  Do i1=1,n_S0
   Call CoupScalarW(i1, gSU2, vevs, RS0, cpl_S0WW(i1) )
  End Do

  !----------------
  ! scalar - Z - Z
  !----------------
  cpl_S0ZZ = 0._dp
  Do i1=1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevs, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp
  cpl_S0SlSl = 0._dp
  cpl_S0SnSn = 0._dp

   Do i1=1,n_S0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

     Call CoupScalarSfermion3(i1, 1, 1, RS0, 0.5_dp, 0._dp, ZeroC, id3C  &
                               &, ZeroC, mu, vevs, gU1, gSU2, coupC )
     cpl_S0SnSn(i1,i2,i2) = coupC

     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, -1._dp, Yuk, Rsl  &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  !----------
  ! scalar W
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarW(i1, gSU2, vevs, RS0, cpl_S0WW(i1) )
  End Do

  !----------
  ! scalar Z
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevs, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SlSlZ = 0._dp
  cpl_SnSnZ = 0._dp
  cpl_SuSuZ = 0._dp

   Call CoupSneutrinoZ(gSU2, sinW2, cpl_SnSnZ(1,1) )
   Do i1=1,3
    cpl_SnSnZ(i1,i1) = cpl_SnSnZ(1,1)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSleptonZ(i2, i3, gSU2, sinW2, Rsl, coupC )
      cpl_SlSlZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  Iname = Iname - 1

 End Subroutine AllCouplingsMSSM1G

 Subroutine CoupChargedScalarPseudoscalar3(i, j, k, RSpm, RP0, vevs, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between pseudoscalar and charged scalars
 ! valid for the MSSM
 !  i,j ........ indices of charged scalar bosons: i=-, j=+
 !  k .......... index of the pseuoscalar boson
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vacuum expectation values (v_d,v_u)
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 01.05.2001
 !-----------------------------------------------------------------------
 Implicit None

 Real(dp), Intent(in) :: RP0(2,2), vevs(2), g
 Integer, Intent(in) :: i, j, k
 Complex(dp), Intent(in) :: RSpm(2,2)
 Complex(dp), Intent(out) :: coup

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarPseudoscalar3'

 If ((i.Lt.1).Or.(i.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i
  Call TerminateProgram
 End If
 If ((j.Lt.1).Or.(j.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index j out of range: (j,n_Spm) = ',j
  Call TerminateProgram
 End If
 If ((k.Lt.1).Or.(k.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index k out of range: (k,n_S0) = ',k
  Call TerminateProgram
 End If

 coup = (0._dp, 0.25_dp) * g**2 * (vevs(2) * RP0(k,1) + vevs(1) * RP0(k,2) ) &
      &  * (Conjg( Rspm(i,2) ) * Rspm(j,1) - Conjg( Rspm(i,1) ) * Rspm(j,2) )
 
 Iname = Iname - 1

 End Subroutine CoupChargedScalarPseudoscalar3

 Subroutine CoupChargedScalarPseudoscalar4(i, j, k, l, RSpm, RP0, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between pseudoscalar and charged scalars
 ! valid for the MSSM
 !  i,j ........ indices of charged scalar bosons: i=-, j=+
 !  k,l ........ index of the pseuoscalar boson
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 10.10.01
 !-----------------------------------------------------------------------
 Implicit None

 Real(dp), Intent(in) :: RP0(2,2), g, gp
 Integer, Intent(in) :: i, j, k, l
 Complex(dp), Intent(in) :: RSpm(2,2)
 Complex(dp), Intent(out) :: coup

 Real(dp) :: g2sum, g2diff, g2
 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarPseudoscalar4'

 g2 = g**2
 g2diff = g2 - gp**2
 g2sum = g2 + gp**2

  coup = g2sum * ( Conjg( RSpm(i,1) ) * RSpm(j,1) * RP0(k,1) * RP0(l,1)    &
     &           + Conjg( RSpm(i,2) ) * RSpm(j,2) * RP0(k,2) * RP0(l,2) )  &
     & + g2diff * ( Conjg( RSpm(i,1) ) * RSpm(j,1) * RP0(k,2) * RP0(l,2)   &
     &            + Conjg( RSpm(i,2) ) * RSpm(j,2) * RP0(k,1) * RP0(l,1) )

  If (k.Eq.l) Then
   coup = coup                                                             &
     & - g2 * (Conjg(RSpm(i,2))*RSpm(j,1) + Conjg(RSpm(i,1))* RSpm(j,2) ) &
     &      * 2._dp * RP0(k,1) * RP0(l,2)
  Else
   coup = coup                                                             &
     & - g2 * (Conjg(RSpm(i,2))*RSpm(j,1) + Conjg(RSpm(i,1))* RSpm(j,2) ) &
     &      * (RP0(k,1) * RP0(l,2) + RP0(k,2) * RP0(l,1) )
  End If
 
 coup = 0.25_dp * coup

 Iname = Iname - 1

 End Subroutine CoupChargedScalarPseudoscalar4

 Subroutine CoupChargedScalarScalar3(i, j, k, RSpm, RS0, vevs, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between neutral scalar and charged scalars
 ! valid for the MSSM
 !  i,j ........ indices of charged scalar bosons: i=-, j=+
 !  k .......... index of the scalar boson
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vacuum expectation values
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Real(dp), Intent(in) :: RS0(2,2), vevs(2), g, gp
 Integer, Intent(in) :: i, j, k
 Complex(dp), Intent(in) :: RSpm(2,2)
 Complex(dp), Intent(out) :: coup

 Complex(dp) :: sum1, sum2
 Real(dp) :: g2diff, g2sum, g2, gp2

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarScalar3'

 If ((i.Lt.1).Or.(i.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i
  Call TerminateProgram
 End If
 If ((j.Lt.1).Or.(j.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index j out of range: (j,n_Spm) = ',j
  Call TerminateProgram
 End If
 If ((k.Lt.1).Or.(k.Gt.2)) Then
  Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
  Write(ErrCan,*) 'index k out of range: (k,n_S0) = ',k
  Call TerminateProgram
 End If

 g2 = 0.25_dp * g**2 
 gp2 = 0.25_dp * gp**2
 g2diff = g2 - gp2
 g2sum =  g2 + gp2

 sum1 = - g2 * vevs(2) * ( Rspm(i,2) * Conjg( Rspm(j,1) )           &
   &                     + Rspm(i,1) * Conjg( Rspm(j,2) ) )         &
   &  - vevs(1) * ( g2sum *  Rspm(i,1) * Conjg( Rspm(j,1) )         &
   &              + g2diff * Rspm(i,2) * Conjg( Rspm(j,2) ) )
 sum1 = sum1 * RS0(k,1)

 sum2 = - g2 * vevs(1) * ( Rspm(i,2) * Conjg( Rspm(j,1) )           &
   &                     + Rspm(i,1) * Conjg( Rspm(j,2) ) )         &
   &  - vevs(2) * ( g2diff *  Rspm(i,1) * Conjg( Rspm(j,1) )        &
   &              + g2sum * Rspm(i,2) * Conjg( Rspm(j,2) ) )
 sum2 = sum2 * RS0(k,2)

 coup =  sum1 + sum2 
 
 Iname = Iname - 1

 End Subroutine CoupChargedScalarScalar3

 Subroutine CoupChargedScalarScalar4(i, j, k, l, RSpm, RS0, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between pseudoscalar and charged scalars
 ! valid for the MSSM
 !  i,j ........ indices of charged scalar bosons: i=-, j=+
 !  k,l ........ index of the pseuoscalar boson
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 10.10.01
 !-----------------------------------------------------------------------
 Implicit None

 Real(dp), Intent(in) :: RS0(2,2), g, gp
 Integer, Intent(in) :: i, j, k, l
 Complex(dp), Intent(in) :: RSpm(2,2)
 Complex(dp), Intent(out) :: coup

 Real(dp) :: g2sum, g2diff, g2
 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarScalar4'

 g2 = g**2
 g2diff = g2 - gp**2
 g2sum = g2 + gp**2

  coup = g2sum * ( Conjg( RSpm(i,1) ) * RSpm(j,1) * RS0(k,1) * RS0(l,1)    &
     &           + Conjg( RSpm(i,2) ) * RSpm(j,2) * RS0(k,2) * RS0(l,2) )  &
     & + g2diff * ( Conjg( RSpm(i,1) ) * RSpm(j,1) * RS0(k,2) * RS0(l,2)   &
     &            + Conjg( RSpm(i,2) ) * RSpm(j,2) * RS0(k,1) * RS0(l,1) )

  If (k.Eq.l) Then
   coup = coup                                                             &
     & + g2 * (Conjg(RSpm(i,2))*RSpm(j,1) +  Conjg(RSpm(i,1))* RSpm(j,2) ) &
     &      * 2._dp * RS0(k,1) * RS0(l,2)
  Else
   coup = coup                                                             &
     & + g2 * (Conjg(RSpm(i,2))*RSpm(j,1) +  Conjg(RSpm(i,1))* RSpm(j,2) ) &
     &      * (RS0(k,1) * RS0(l,2) + RS0(k,2) * RS0(l,1) )
  End If
 
 coup = 0.25_dp * coup

 Iname = Iname - 1

 End Subroutine CoupChargedScalarScalar4



 Subroutine CoupCharginoPseudoScalarMSSMa(i,j,k,U,V,RP0,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for the MSSM
 !  i,j ........ index of chargino, first index is a electroweak eigenstate
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RP0(2,2)
  Complex(dp), Intent(in) :: U(2,2), V(2,2)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarMSSMa'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  If (i.Eq.1) Then
   coupL = RP0(k,1) * Conjg(U(j,2) )
   coupR = RP0(k,2) * V(j,2)
  Else
   coupL = RP0(k,2) * Conjg(U(j,1))
   coupR = RP0(k,1) * V(j,1)
  End If
 
  coupL = g * coupL * (0._dp,1._dp) * oosqrt2
  coupR = g * coupR * (0._dp,-1._dp) * oosqrt2 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarMSSMa


 Subroutine CoupCharginoPseudoScalarMSSM(i,j,k,U,V,RP0,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for the MSSM
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RP0(2,2)
  Complex(dp), Intent(in) :: U(2,2),V(2,2)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = g * ( RP0(k,1) * Conjg( V(i,1) * U(j,2) )             &
        &     + RP0(k,2) * Conjg( V(i,2) * U(j,1) ) ) * oosqrt2
  coupR = g * ( RP0(k,1) * V(j,1) * U(i,2)             &
        &     + RP0(k,2) * V(j,2) * U(i,1) ) * oosqrt2
 
  coupL = coupL * (0._dp,1._dp)
  coupR = coupR * (0._dp,-1._dp) 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarMSSM



 Subroutine CoupCharginoScalarMSSMa(i,j,k,U,V,RS0,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for the MSSM
 !  i,j ........ index of chargino, first index is an electroweak eigenstate
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 09.11.01
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RS0(2,2)
  Complex(dp), Intent(in) :: U(2,2), V(2,2)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarMSSMa'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  If (i.Eq.1) Then
   coupL = RS0(k,1) * Conjg( U(j,2) )
   coupR = RS0(k,2) * Conjg( V(j,2) )
  Else
   coupL = RS0(k,2) * Conjg( U(j,1) )
   coupR = RS0(k,1) * Conjg( V(j,1) )
  End If
 
  coupL = - g * coupL * oosqrt2
  coupR = - g * Conjg( coupR ) * oosqrt2

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarMSSMa


 Subroutine CoupCharginoScalarMSSM(i,j,k,U,V,RS0,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for the MSSM
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RS0(2,2)
  Complex(dp), Intent(in) :: U(2,2),V(2,2)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = g * ( RS0(k,1) * Conjg( V(i,1) * U(j,2) )             &
        &     + RS0(k,2) * Conjg( V(i,2) * U(j,1) ) ) * oosqrt2
  coupR = g * ( RS0(k,1) * Conjg( V(j,1) * U(i,2) )             &
        &     + RS0(k,2) * Conjg( V(j,2) * U(i,1) ) ) * oosqrt2
 
  coupL = - coupL
  coupR = - Conjg( coupR )

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarMSSM


 Subroutine CoupCharginoSfermion1(i,j,g,T3,RSf,YukD,YukU,U,V,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between chargino, fermion, and sfermion'
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of chargino
 !  j .......... index of sfermion
 !  g .......... SU(2) gauge coupling
 !  T3 ......... weak isospin of the fermion
 !  RSf(i,j) ... mixing matrix of the sfermion
 !  YukD ....... Yukawa coupling of the T3=-0.5 fermion 
 !  YukU ....... Yukawa coupling of the T3=0.5 fermion 
 !  U,V ........ mixing matrices of the chargino
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! the corresponding part of the Lagrangian is:
 !  \bar{f} (coupL P_L + coupR P_R) chi^\pm_i \tilde{f'}_j
 ! written by Werner Porod, 5.8.1999
 ! 6.4.2000: changing error system
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,T3
  Complex(dp), Intent(in) :: RSf(2,2),U(:,:),V(:,:),YukD,YukU
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j

  Integer :: n_char

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoSfermion1' 

  n_char = Size(U, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: (i,n_char) = ',i,n_char
   Call TerminateProgram
  End If 
  If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range: ',j
   Call TerminateProgram
  End If 

  If ( T3.Gt.0.) Then
   coupL = YukU * Conjg( RSf(j,1) ) * Conjg( V(i,2) )

   coupR = - g * Conjg( RSf(j,1) ) * U(i,1)           &
     &   + Conjg( YukD) * Conjg( RSf(j,2) ) * U(i,2)

  Else
   coupL = YukD * Conjg( RSf(j,1) ) * Conjg( U(i,2) )
   coupR = - g * Conjg( RSf(j,1) ) * V(i,1)           &
     &   + Conjg( YukU) * Conjg( RSf(j,2) ) * V(i,2)

  End If

  Iname = Iname - 1

 End Subroutine CoupCharginoSfermion1

 Subroutine CoupCSCharginoNeutralinoMSSM(k, i, j, N, U, V, RSpm &
                                       &, gp, g, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and
 ! charged scalar bosons 
 ! valid for the MSSM
 !  i .......... index of chargino
 !  j .......... index of the neutralino
 !  k .......... index of the charged scalar
 !  N.. ........ mixing matrix of the neutralino
 !  U,V ........ mixing matrices of the chargino
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coupL ...... the left coupling(k,i,j)
 !  coupR ...... the right coupling(k,i,j)
 ! the lagrangian is given by
 !  S^-(k) \bar{\chim(i)} (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 5.8.1999
 ! 26.04.2001: porting to f90, including also possibilities to extend for
 !             R-parity violation
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) ::  i, j, k
  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:)
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sumI, sumIC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCSCharginoNeutralinoMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Chargino index out of range (i,n_char) = ',i
   Call TerminateProgram
  Elseif ((j.Lt.1).Or.(j.Gt.4)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Neutralino index out of range (i,n_neut) = ',j
   Call TerminateProgram
  Elseif ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Scalar index out of range (i,n_Spm) = ',k
   Call TerminateProgram
  End If

  coupL = ZeroC
  coupR = ZeroC

  sumI = oosqrt2 * ( gp * N(j,1) + g * N(j,2) )
  sumIC = Conjg(sumI)

  coupL = - RSpm(k,2) * ( g * Conjg( V(i,1) * N(j,4) ) &
        &               + Conjg( V(i,2) ) * sumIC)
  coupR = RSpm(k,1) * (- g * U(i,1) * N(j,3) + U(i,2) * sumI )
 
  Iname = Iname - 1

 End Subroutine CoupCSCharginoNeutralinoMSSM



 Subroutine CoupFermionPseudoscalar1(i,T3,yuk,RP0,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between SM-fermions and pseudoscalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of the pseudoscalar boson
 !  T3 ......... weak isospin of the fermion
 !  Yuk ........ lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3,RP0(:,:)
  Complex(dp), Intent(in) :: yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i
 
  If (T3.Gt.0._dp) Then
   coupL = - (0._dp,1._dp) * yuk * RP0(i,2) * oosqrt2
  Else
   coupL = - (0._dp,1._dp) * yuk * RP0(i,1) * oosqrt2
  End If
  coupR = Conjg( coupL )

 End Subroutine CoupFermionPseudoscalar1


 Subroutine CoupFermionScalar1(i,T3,yuk,RS0,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between SM-fermions and scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model, and the spontaneously broken model
 !  i .......... index of scalar boson
 !  T3 ......... weak isospin of the fermion
 !  Yuk ........ yukawa couplings
 !  RS0(i,j) ... mixing matrix of scalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3,RS0(:,:)
  Complex(dp), Intent(in) :: yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i
 
  If (T3.Gt.0._dp) Then
   coupL = - yuk * RS0(i,2) * oosqrt2
  Else
   coupL = - yuk * RS0(i,1) * oosqrt2
  End If
  coupR = Conjg(coupL) 

 End Subroutine CoupFermionScalar1

 Subroutine CoupGravitinoSfermion1(j, RSf, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the overall coupling gravitino-sfermion-fermion 
 ! and without the generator T^a_ij, the Lagrangian is given by:
 !  L = \bar{f_i} (coupR P_R + coupL P_L) gamma_mu gamma^nu G^mu D_\nu 
 !            *  \tilde f_j / (Sqrt[2] m_Planck)
 ! input:
 !  g3 ............ SU(3) coupling
 !  phase_M3 ...... phase of the parameter M3
 !  i ............. index of the fermion
 !  j ............. index of the sfermion
 !  RSf(i,j) ...... sfermion mixing matrix
 !  RL(i,j) ....... left fermion mixing matrix
 !  RR(i,j). ...... right fermion mixing matrix
 ! output
 !  coupL ......... left coupling
 !  coupR ......... right coupling
 ! written by Werner Porod, 18.09.10
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: j
  Complex(dp), Intent(in) :: RSf(2,2)
  Complex(dp), Intent(out) :: coupL, coupR

  coupL = Conjg( RSf(j,2) )
  coupR = Conjg( RSf(j,1) )

 End Subroutine CoupGravitinoSfermion1

 Subroutine CoupGluinoSquark1(g3, phase_M3, i_sq, Rsq, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling gluino-squark-quark without generation mixing
 ! and without the generator T^a_ij, the Lagrangian is given by:
 !  L = g3 T^a_ij \bar{q} (coupR P_R + coupL P_L) g^a \tilde q_(i_sq)
 ! input:
 !  g3 ............ SU(3) coupling
 !  phase_M3 ...... phase of the parameter M3
 !  i_sq .......... index of the squark
 !  Rsq(i,j) ...... squark mixing matrix
 ! output
 !  coupL ......... left coupling
 !  coupR ......... right coupling
 ! the lagrangian is given by
 !  Lambda^a \bar{q} ( coupL P_L + coupR P_R) g^a \tilde{q}_i_sq
 ! where Lambda^a are the Gell-Mann matrices
 ! written by Werner Porod, 21.1.01
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g3
  Integer, Intent(in) :: i_sq
  Complex(dp), Intent(in) :: phase_M3, Rsq(2,2)
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sqrt_phase

  sqrt_phase = Sqrt( phase_M3 )
  coupL = g3 * Conjg( Rsq(i_sq,2) ) * oosqrt2 / sqrt_phase
  coupR = - g3 * oosqrt2 * sqrt_phase * Conjg( Rsq(i_sq,1) )

 End Subroutine CoupGluinoSquark1


 Subroutine CoupNeutralinoSdown1(i,j,gp,g,RSf,Yuk,N,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, down, and sdown
 ! valid for the MSSM, 1-generation epsilon model, 
 ! 3-generation epsilon model, and the full spontaneously broken model
 !  i .......... index of neutralino
 !  j .......... index of sdown
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RSf(i,j) ... mixing matrix of the sdown
 !  Yuk ........ Yukawa coupling of the down 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 !    L = \bar(f) ( coupL P_L + coupR P_R) chi^0_i \tilde(f)_j  
 ! written by Werner Porod, 6.8.1999
 ! 23.3.2000: including the full model with spontaneously broken R-parity
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: gp,g
  Complex(dp), Intent(in) :: RSf(2,2),N(:,:),Yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j

  Real(dp) :: e,T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSdown1'

  T3 = -0.5_dp
  e = -1._dp / 3._dp

  Call CoupNeutralinoSfermion1(i,j,gp,g,T3,e,RSf,Yuk,N,coupL,coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSdown1


 Subroutine CoupNeutralinoSfermion1(i,j,gp,g,T3,e,RSf,Yuk,N,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, fermion, and sfermion
 ! valid for the MSSM, 1-generation epsilon model,
 ! 3-generation epsilon model, and the full spontaneously broken model
 !  i .......... index of neutralino
 !  j .......... index of sfermion
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  T3 ......... weak isospin of the fermion
 !  e .......... charge of the fermion
 !  RSf(i,j) ... mixing matrix of the sfermion
 !  Yuk ........ Yukawa coupling of the fermion 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 !                                      k=2 ... right coupling
 ! written by Werner Porod, 6.8.1999
 ! 23.3.2000: including the full model with spontaneously broken R-parity
 ! 9.10.2000: porting to f90
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f) ( coupL P_L + coupR P_R) chi^0_i \tilde(f)_j  
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,T3,gp,e
  Complex(dp), Intent(in) :: RSf(2,2),N(:,:),Yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j

  Complex(dp) :: fl,fr,hl,hr

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSfermion1'

  If ( T3.Gt.0.) Then
   hl = - Yuk * Conjg( N(i,4) )
  Else
   hl = - Yuk * Conjg( N(i,3) )
  End If
  hr = Conjg(hl)

  fL = - sqrt2 *( (e-T3) * gp * N(i,1) + T3 * g * N(i,2) )
  fR = sqrt2 * e * gp * Conjg( N(i,1) ) 

  coupL = Conjg( Rsf(j,1) ) * hl +  Conjg( Rsf(j,2) ) * fr 
  coupR = Conjg( Rsf(j,1) ) * fl +  Conjg( Rsf(j,2) ) * hr 

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSfermion1


 Subroutine CoupNeutralinoSlepton1(i,j,gp,g,RSf,Yuk,N,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, lepton, and slepton
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of neutralino
 !  j .......... index of slepton
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RSf(i,j) ... mixing matrix of the slepton
 !  Yuk ........ Yukawa coupling of the lepton 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 !    L = \bar(f) ( coupL P_L + coupR P_R) chi^0_i \tilde(f)_j  
 ! written by Werner Porod, 6.8.1999
 ! 23.3.2000: including the full model with spontaneously broken R-parity
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: gp,g
  Complex(dp), Intent(in) :: RSf(2,2),N(:,:),Yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j

  Real(dp) :: e,T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSlepton1'

  T3 = -0.5_dp
  e = -1._dp 

  Call CoupNeutralinoSfermion1(i,j,gp,g,T3,e,RSf,Yuk,N,coupL,coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSlepton1


 Subroutine CoupNeutralinoSneutrino1(i,gp,g,N,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, neutrino, and sneutrino
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of neutralino
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupR ...... the right coupling(i)
 ! the lagrangian is given by:
 !  \bar{\nu} ( coupR P_R ) \chi^0_i \tilde{\nu}
 ! written by Werner Porod, 6.8.1999
 ! 23.3.2000: including the full model with spontaneously broken R-parity
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: gp, g
  Complex(dp), Intent(in) :: N(:,:)
  Complex(dp), Intent(out) :: coupR
  Integer, Intent(in) :: i

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSneutrino1'

  coupR = sqrt2 * 0.5_dp * ( gp * N(i,1) - g * N(i,2) )

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSneutrino1


 Subroutine CoupNeutralinoSup1(i,j,gp,g,RSf,Yuk,N,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, up, and sup
 ! valid for the MSSM, 1-generation epsilon model,
 ! 3-generation epsilon model, and the full spontaneously broken model
 !  i .......... index of neutralino
 !  j .......... index of sup
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RSf(i,j) ... mixing matrix of the sup
 !  Yuk ........ Yukawa coupling of the up 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 !                                      k=2 ... right coupling
 ! written by Werner Porod, 6.8.1999
 ! 23.3.2000: including the full model with spontaneously broken R-parity
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: gp,g
  Complex(dp), Intent(in) :: RSf(2,2),N(:,:),Yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j

  Real(dp) :: e,T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSup1'

  T3 = 0.5_dp
  e = 2._dp / 3._dp

  Call CoupNeutralinoSfermion(i,j,gp,g,T3,e,RSf,Yuk,N,coupL,coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSup1


 Subroutine CoupPseudoScalarScalar3(i, j, k, RP0, RS0, gp, g, vevs, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between neutral scalar and pseudoscalars
 ! valid for the MSSM
 !  i,j ........ indices of pseudo scalar bosons
 !  k .......... index of the scalar boson
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  vevs(i) .... vacuum expectation values
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j, k
 Real(dp), Intent(in) :: RS0(2,2), RP0(2,2), vevs(2), g, gp
 Real(dp), Intent(out) :: coup

 coup = - 0.125_dp * (g**2 + gp**2)                               & 
      &           * (RP0(i,1) * RP0(j,1) - RP0(i,2) * RP0(j,2))  &
      &           * (vevs(1) * RS0(k,1) - vevs(2) * RS0(k,2) )
 
 End Subroutine CoupPseudoScalarScalar3


 Subroutine CoupPseudoScalarScalar4(i, j, k, l, RP0, RS0, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between neutral scalar and pseudoscalars
 ! valid for the MSSM
 !  i,j ........ indices of pseudo scalar bosons
 !  k,l ........ index of the scalar boson
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 10.10.01
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j, k, l
 Real(dp), Intent(in) :: RS0(2,2), RP0(2,2), g, gp
 Real(dp), Intent(out) :: coup

 coup = - 0.125_dp * (g**2 + gp**2)                                & 
      &            * (RP0(i,1) * RP0(j,1) - RP0(i,2) * RP0(j,2) )  &
      &            * (RS0(l,1) * RS0(k,1) - RS0(l,2) * RS0(k,2) )
 
 End Subroutine CoupPseudoScalarScalar4


 Subroutine CoupPseudoScalarSfermion3_1(i, j, k, RP0, T3, yuk, Rsf, A, bi &
                                     & , coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a pseudoscalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 15.9.1999
 ! 15.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2), yuk, A, bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RP0(:,:), T3
  Integer, Intent(in) :: i, j, k

  Integer :: n_P0, i1
  Complex(dp) :: ayuk, aA, abi(Size(bi)), parts(Size(RP0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion3_1'

  coup = ZeroC

  n_P0 =  Size(RP0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  If (T3.Gt.0._dp) Then
   parts(1) = yuk * abi(1) * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          & - ayuk * bi(1) * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(2) =  A * Rsf(k,2) * Conjg( Rsf(j,1) )     &
          & -  aA * Rsf(k,1) * Conjg( Rsf(j,2) )

   Do i1=3,n_P0
    parts(i1) = - yuk * abi(i1-1) * Rsf(k,2) * Conjg( Rsf(j,1) ) &
            &   + ayuk * bi(i1-1) * Rsf(k,1) * Conjg( Rsf(j,2) )
   End Do

  Else

   parts(1) =  A * Rsf(k,2) * Conjg( Rsf(j,1) ) &
          & - aA * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(2) = yuk * abi(1) * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          & - ayuk * bi(1) * Rsf(k,1) * Conjg( Rsf(j,2) )

  End If

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * parts(i1)
  End Do

  coup = (0._dp, -1._dp) * coup *oosqrt2

  Iname = Iname - 1

 End Subroutine CoupPseudoScalarSfermion3_1

 Subroutine CoupPseudoScalarSfermion3a_1(i, j, k, RP0, T3, yuk, Rsf, A, bi &
                                     & , coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a pseudoscalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}, is electroweak eigenstate
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 30.12.01
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2), yuk, A, bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RP0(:,:), T3
  Integer, Intent(in) :: i, j, k

  Integer :: n_P0, i1
  Complex(dp) :: ayuk, aA, abi(Size(bi)), parts(Size(RP0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion3a_1'

  coup = ZeroC

  n_P0 =  Size(RP0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  If (T3.Gt.0._dp) Then
   parts(1) = yuk * abi(1) * Rsf(k,2) * Conjg( id2c(j,1) )  &
          & - ayuk * bi(1) * Rsf(k,1) * Conjg( id2c(j,2) )

   parts(2) =  A * Rsf(k,2) * Conjg( id2c(j,1) )     &
          & -  aA * Rsf(k,1) * Conjg( id2c(j,2) )

   Do i1=3,n_P0
    parts(i1) = - yuk * abi(i1-1) * Rsf(k,2) * Conjg( id2c(j,1) ) &
            &   + ayuk * bi(i1-1) * Rsf(k,1) * Conjg( id2c(j,2) )
   End Do

  Else

   parts(1) =  A * Rsf(k,2) * Conjg( id2c(j,1) ) &
          & - aA * Rsf(k,1) * Conjg( id2c(j,2) )

   parts(2) = yuk * abi(1) * Rsf(k,2) * Conjg( id2c(j,1) )  &
          & - ayuk * bi(1) * Rsf(k,1) * Conjg( id2c(j,2) )

  End If

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * parts(i1)
  End Do

  coup = (0._dp, -1._dp) * coup *oosqrt2

  Iname = Iname - 1

 End Subroutine CoupPseudoScalarSfermion3a_1

 Subroutine CoupPseudoscalarSfermion4_1(i, j, k, l, RP0, T3, e, yuk,  &
   &                                    Rsf, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between pseudoscalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of pseudoscalar boson
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  n_P0 ....... specifies the model: 2 -> MSSM
 !                                    3 -> 1-generation epsilon model
 !                                    5 -> 3-generation epsilon model
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 29.11.1999
 ! 10.10.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2), yuk
  Real(dp), Intent(in) :: RP0(:,:), g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: n_P0, i1
  Real(dp) :: g2, gp2, YL, YR, yuk2
  Complex(dp) :: parts(12), Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion4_1'

  coup = Cmplx( 0._dp,0._dp,dp )

  n_P0 = Size(RP0, Dim=1)

  yuk2 = Abs( yuk )**2

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = ZeroC

  Dterm = 0.25_dp *( (T3* g2 - YL* gp2) * Rsf(l,1) *Conjg( Rsf(k,1) )  &
        &          - YR * gp2 * Rsf(l,2) * Conjg( Rsf(k,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = -  Dterm

   parts(2) = Dterm                  & 
          & - 0.5_dp * yuk2          &
          &   * ( Rsf(l,1) * Conjg( Rsf(k,1) ) + Rsf(l,2) * Conjg( Rsf(k,2) ) )

  Else

   parts(1) = - Dterm              &
         &  - 0.5_dp * yuk2        &
         &    * ( Rsf(l,1) * Conjg( Rsf(k,1) ) + Rsf(l,2) * Conjg( Rsf(k,2) ) )

   parts(2) = Dterm

  Endif

  Do i1=3,n_P0
   parts(i1) = - Dterm
  Enddo

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * RP0(j,i1) * parts(i1)
  Enddo

  Iname = Iname - 1

 End Subroutine CoupPseudoscalarSfermion4_1


 Subroutine CoupScalar3a(i, j, k, RS0, gp, g, vevSM, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between neutral scalars
 ! valid for the MSSM and the epsilon models
 ! 3-generation epsilon model
 !  i,j,k ...... indices of scalar bosons
 !               here I assume that the index i belongs to an electroweak
 !               eigenstate
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  vevs(i) .... vacuum expectation values
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 10.10.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j, k
 Real(dp), Intent(in) :: RS0(2,2), vevSM(2), g, gp
! real(dp), intent(in), optional :: vevL(:)
 Real(dp), Intent(out) :: coup

! integer :: i1, n_S0
 Real(dp) :: sumJK, sumIK, sumIJ

! n_S0 = Size( RS0, Dim = 1)

 If (i.Eq.1) Then
  sumIJ = RS0(j,1) 
  sumIK = RS0(k,1) 
  sumJK = RS0(j,1) * RS0(k,1) - RS0(j,2) * RS0(k,2)

  coup = vevSM(1) * ( sumJK + RS0(j,1)*sumIK + RS0(k,1)*sumIJ ) &
     & - vevSM(2) * ( RS0(j,2)*sumIK +  RS0(k,2)*sumIJ )

 Else If (i.Eq.2) Then
  sumIJ = - RS0(j,2)
  sumIK = - RS0(k,2)
  sumJK = RS0(j,1) * RS0(k,1) - RS0(j,2) * RS0(k,2)

  coup = vevSM(1) * ( RS0(j,1)*sumIK +  RS0(k,1)*sumIJ )       &
     & - vevSM(2) * ( sumJK +  RS0(j,2)*sumIK +  RS0(k,2)*sumIJ )

 End If

 coup = - 0.25_dp * (g**2 + gp**2) * coup
 
 End Subroutine CoupScalar3a


 Subroutine CoupScalar3(i, j, k, RS0, gp, g, vevSM, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between neutral scalars
 ! valid for the MSSM and the epsilon models
 ! 3-generation epsilon model
 !  i,j,k ...... indices of scalar bosons
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  vevs(i) .... vacuum expectation values
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i, j, k
 Real(dp), Intent(in) :: RS0(2,2), vevSM(2), g, gp
! real(dp), intent(in), optional :: vevL(:)
 Real(dp), Intent(out) :: coup

! integer :: i1, n_S0
 Real(dp) :: sumJK, sumIK, sumIJ

! n_S0 = Size( RS0, Dim = 1)

 sumIJ = RS0(i,1) * RS0(j,1) - RS0(i,2) * RS0(j,2)
 sumIK = RS0(i,1) * RS0(k,1) - RS0(i,2) * RS0(k,2)
 sumJK = RS0(j,1) * RS0(k,1) - RS0(j,2) * RS0(k,2)

! do i1 =3, n_S0
!  sumIJ = sumIJ + RS0(i,i1) * RS0(j,i1)
!  sumIK = sumIK + RS0(i,i1) * RS0(k,i1)
!  sumJK = sumJK + RS0(j,i1) * RS0(k,i1)
! end do

 coup = vevSM(1) * ( RS0(i,1)*sumJK +  RS0(j,1)*sumIK +  RS0(k,1)*sumIJ ) &
    & - vevSM(2) * ( RS0(i,2)*sumJK +  RS0(j,2)*sumIK +  RS0(k,2)*sumIJ )

! do i1 =3, n_S0
!  coup = coup + vevL(i1-2)  &
!                * ( RS0(i,i1)*sumJK +  RS0(j,i1)*sumIK +  RS0(k,i1)*sumIJ )
! end do

 coup = - 0.25_dp * (g**2 + gp**2) * coup
 
 End Subroutine CoupScalar3


 Subroutine CoupScalar4a(i, j, k, l, RS0, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between neutral scalars
 ! valid for the MSSM
 !  i,j,k,l .... indices of  scalar bosons
 !               assuming that i and j are electroweak states
 !  RS0(i,j) ... mixing matrix of  scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 2.10.1999
 ! 09.11.01: poratation to f90 
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: RS0(2,2), g, gp
  Integer, Intent(in) :: i, j, k, l
  Real(dp), Intent(out) :: coup

  Real(dp) :: sum1

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalar4a'

  If ((i.Eq.1).And.(j.Eq.1)) Then
   sum1 = 3._dp * RS0(k,1) * RS0(l,1) -  RS0(k,2) * RS0(l,2)
  Elseif ((i.Eq.2).And.(j.Eq.2)) Then
   sum1 = 3._dp * RS0(k,2) * RS0(l,2) -  RS0(k,1) * RS0(l,1)
  Else
   sum1 = - RS0(k,1) * RS0(l,2) -  RS0(k,2) * RS0(l,1)
  Endif

  coup = - (g**2 + gp**2) * sum1 / 4._dp

  Iname = Iname - 1

 End Subroutine CoupScalar4a



 Subroutine CoupScalarSfermion3MSSM_1(i, j, k, RS0, T3, e, yuk, Rsf, A, mu, &
                      &          vevSM, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the 1-generation MSSM
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  mu ......... bilinear parameters of the superpotential
 !  vevSM(i) ... vevs of the Higgs bosons
 !  gU1 ........ U(1) coupling
 !  gSU2 ....... SU(2) coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 25.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2) , yuk , A, mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(2,2), T3, e, vevSM(2), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: i1
  Complex(dp) :: ayuk, aA, muC, parts(2)
  Real(dp) :: g2, gp2, YL, YR, Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3MSSM_1'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  muC = Conjg( mu )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,1) * Conjg( Rsf(j,1) ) &
        &         - YR * gp2 * Rsf(k,2) * Conjg( Rsf(j,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = ( yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & -  vevSM(1) * Dterm

   parts(2) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                   &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2      &
          & + vevSM(2) * Dterm                                        &
          & - vevSM(2) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) ) &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

  Else

   parts(1) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                    &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2       &
          & - vevSM(1) * Dterm                                         &
          & - vevSM(1) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) )  &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

   parts(2) = ( yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & + vevSM(2) * Dterm

  End If

  Do i1=1,2
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3MSSM_1


 Subroutine CoupScalarSfermion3MSSMa_1(i, j, k, RS0, T3, e, yuk, Rsf, A, mu, &
                      &          vevSM, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the 1-generation MSSM
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}, is a electroweak eigenstate
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  mu ......... bilinear parameters of the superpotential
 !  vevSM(i) ... vevs of the Higgs bosons
 !  gU1 ........ U(1) coupling
 !  gSU2 ....... SU(2) coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 25.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2) , yuk , A, mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(2,2), T3, e, vevSM(2), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: i1
  Complex(dp) :: ayuk, aA, muC, parts(2), Dterm
  Real(dp) :: g2, gp2, YL, YR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3MSSMa_1'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  muC = Conjg( mu )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,1) * Conjg( id2c(j,1) ) &
        &         - YR * gp2 * Rsf(k,2) * Conjg( id2c(j,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = ( yuk * muC * Rsf(k,2) * Conjg( id2c(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2  &
          & -  vevSM(1) * Dterm

   parts(2) = - ( A * Rsf(k,2) * Conjg( id2c(j,1) )                   &
          &     + aA * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2      &
          & + vevSM(2) * Dterm                                        &
          & - vevSM(2) * yuk * ayuk * ( Rsf(k,1) * Conjg( id2c(j,1) ) &
          &                           + Rsf(k,2) * Conjg( id2c(j,2) ) )

  Else

   parts(1) = - ( A * Rsf(k,2) * Conjg( id2c(j,1) )                    &
          &     + aA * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2       &
          & - vevSM(1) * Dterm                                         &
          & - vevSM(1) * yuk * ayuk * ( Rsf(k,1) * Conjg( id2c(j,1) )  &
          &                           + Rsf(k,2) * Conjg( id2c(j,2) ) )

   parts(2) = ( yuk * muC * Rsf(k,2) * Conjg( id2c(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2  &
          & + vevSM(2) * Dterm

  End If

  Do i1=1,2
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3MSSMa_1

 Subroutine CoupScalarSfermion4MSSM_1(i, j, k, l, RS0, T3, e, yuk, Rsf &
                                    &, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between neutral scalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of scalar boson
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 29.11.1999
 ! -09.11.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2), yuk
  Real(dp), Intent(in) :: RS0(2,2), g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: i1
  Complex(dp) :: parts(2), Dterm
  Real(dp) :: g2, gp2, YL, YR, yuk2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion4MSSM_1'

  coup = Cmplx( 0._dp,0._dp,dp )

  yuk2 = Abs( yuk )**2

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = Cmplx( 0._dp,0._dp,dp )

  Dterm = 0.25_dp * (( T3 * g2 - YL* gp2)* Rsf(l,1)* Conjg( Rsf(k,1) )  &
        &             - YR * gp2 * Rsf(l,2) * Conjg( Rsf(k,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = - Dterm

   parts(2) = Dterm  &
          & - 0.5_dp * yuk2 * ( Rsf(l,1) * Conjg( Rsf(k,1) )  &
          &                   + Rsf(l,2) * Conjg( Rsf(k,2) ) )

  Else

   parts(1) = - Dterm                                           &
            & - 0.5_dp * yuk2 * ( Rsf(l,1) * Conjg( Rsf(k,1) )  &
            &                   + Rsf(l,2) * Conjg( Rsf(k,2) ) )

   parts(2) = Dterm

  Endif

  Do i1=1,2
   coup = coup + RS0(i,i1) * RS0(j,i1) * parts(i1)
  Enddo
  Iname = Iname - 1

 End Subroutine CoupScalarSfermion4MSSM_1

 Subroutine CoupScalarWMSSM(i, g, vevs, RS0, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between scalar boson and W-boson
 ! valid for the MSSM + NMSSM
 !  i .......... index of scalar boson
 !  g .......... SU(2) gauge coupling
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  vevs ....... vevs of the neutral scalar fields: i=1 -> H_d
 !                                             i=2 -> H_u
 ! output
 !  coup ....... the coupling(i)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i
 Real(dp), Intent(in) :: g, vevs(2), RS0(:,:)
 Real(dp), Intent(out) :: coup
 
 coup = 0.5_dp * g**2 * (vevs(1) * RS0(i,1) + vevs(2) * RS0(i,2) )

 End Subroutine CoupScalarWMSSM

 Subroutine CoupScalarZMSSM(i,g,cosW2,vevs,RS0,coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between pseudoscalar boson, scalar boson, and
 ! Z-boson
 ! valid for the MSSM + NMSSM
 !  i .......... index of pseudoscalar boson
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  vevs ....... vevs of the neutral scalar fields: i=1 -> H_d
 !                                                  i=2 -> H_u
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 4.8.1999
 ! 24.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g, cosW2, vevs(2), RS0(:,:)
  Real(dp), Intent(out) :: coup
  Integer, Intent(in) :: i

  coup = 0.5_dp * g**2 * (vevs(1) * RS0(i,1) + vevs(2) * RS0(i,2)) / cosW2

 End Subroutine CoupScalarZMSSM

 Subroutine CoupSfermion4G_13(i, j, k, l, gp,g, T3_1, e1, Rsf1, T3_2, e2 &
                            &, Rsf2, k_sf2, Wpart, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion self coupling for different flavours of
 ! sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i ........... index of the sfermion 1
 !  j ........... index of the adjungated sfermion 1
 !  k ........... index of the sfermion 2
 !  l ........... index of the adjungated sfermion 2
 !  gp .......... U(1) gauge coupling 
 !  g ........... SU(2) gauge coupling 
 !  T3_1, e1 .... isospin and charge of sfermion 1
 !  Rsf1(i,j) ... mixing matrix of the sfermions 1
 !  T3_2, e2 .... isospin and charge of sfermion 2
 !  Rsf2(i,j) ... mixing matrix of the sfermions 2
 !  k_sf2 ....... coloured particle in the loops require factor of 3
 !  Wpart ....... if .true. then also the contribution of W1 and W2 in the
 !                D-term is included
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 01.01.02
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, k_sf2
  Real(dp), Intent(in) :: T3_1, T3_2, gp, g, e1, e2
  Complex(dp), Intent(in) :: Rsf1(:,:), Rsf2(:,:)
  Logical, Intent(in) :: Wpart
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2, n_sf1, n_sf2
  Real(dp) :: YL12, YR12, YL1R2, YR1L2, g2sq
  
  Iname = Iname + 1 
  NameOfUnit(Iname) = "CoupSfermion4G_13"

  YL12 = gp**2 * (T3_1 - e1) * (T3_2 - e2) * k_sf2 &
     & + T3_1 * T3_2 * k_sf2 * g**2
  YR12 = gp**2 * e1 * e2 * k_sf2
  YL1R2 = gp**2 * (T3_1 - e1) * e2 * k_sf2
  YR1L2 = gp**2 * e1 * (T3_2 - e2) * k_sf2

  n_sf1 = Size( Rsf1, dim=1)
  n_sf2 = Size( Rsf2, dim=1)

  If (Wpart) g2sq = 0.5_dp * g**2 ! otherwise I do not need it

  coup = ZeroC
  If ((n_Sf1.Eq.2).And.(n_Sf2.Eq.2)) Then
   If (Wpart) Then
    coup = - (g2sq + YL12) * Rsf1(j,1) * Rsf2(l,1)                        &
       &                               * Conjg(Rsf1(i,1) * Rsf2(k,1))     &
       & - Yr12 * Rsf1(j,2) * Rsf2(l,2) * Conjg(Rsf1(i,2) * Rsf2(k,2))    &
       & - YL1R2 * Rsf1(j,1) * Rsf2(l,2) * Conjg(Rsf1(i,1) * Rsf2(k,2))   &
       & - YR1L2 * Rsf1(j,2) * Rsf2(l,1) * Conjg(Rsf1(i,2) * Rsf2(k,1))
   Else
    coup = - YL12 * Rsf1(j,1) * Rsf2(l,1) * Conjg(Rsf1(i,1) * Rsf2(k,1))  &
       & - Yr12 * Rsf1(j,2) * Rsf2(l,2) * Conjg(Rsf1(i,2) * Rsf2(k,2))    &
       & - YL1R2 * Rsf1(j,1) * Rsf2(l,2) * Conjg(Rsf1(i,1) * Rsf2(k,2))   &
       & - YR1L2 * Rsf1(j,2) * Rsf2(l,1) * Conjg(Rsf1(i,2) * Rsf2(k,1))
   End If

  Else If ((n_Sf1.Eq.3).And.(n_Sf2.Eq.3)) Then
   Do i1=1,3
    Do i2=1,3
     If (Wpart) coup = coup - g2sq * Rsf1(j,i2) * Rsf2(l,i2)              &
           &                      * Conjg(Rsf1(i,i1) * Rsf2(k,i1))
     coup = coup                                                            &
        & - YL12 * Rsf1(j,i1) * Rsf2(l,i2) * Conjg(Rsf1(i,i1) * Rsf2(k,i2))
    End Do
   End Do

  Else If ((n_Sf1.Eq.3).And.(n_Sf2.Eq.6))  Then
   Do i1=1,3
    Do i2=1,3
     If (Wpart) coup = coup - g2sq * Rsf1(j,i2) * Rsf2(l,i2)              &
           &                      * Conjg(Rsf1(i,i1) * Rsf2(k,i1))
     coup = coup                                                            &
        & - YL12 * Rsf1(j,i1) * Rsf2(l,i2) * Conjg(Rsf1(i,i1) * Rsf2(k,i2)) &
        & - YL1R2 * Rsf1(j,i1) * Rsf2(l,i2+3) * Conjg(Rsf1(i,i1) * Rsf2(k,i2+3))
    End Do
   End Do

  Else If ((n_Sf1.Eq.6).And.(n_Sf2.Eq.3))  Then
   Do i1=1,3
    Do i2=1,3
     If (Wpart) coup = coup - g2sq * Rsf1(j,i2) * Rsf2(l,i2)              &
           &                      * Conjg(Rsf1(i,i1) * Rsf2(k,i1))
      coup = coup                                                            &
         & - YL12 * Rsf1(j,i1) * Rsf2(l,i2) * Conjg(Rsf1(i,i1) * Rsf2(k,i2)) &
         & - YR1L2 * Rsf1(j,i1+3) * Rsf2(l,i2) * Conjg(Rsf1(i,i1+3) * Rsf2(k,i2))
    End Do
   End Do

  Else If ((n_Sf1.Eq.6).And.(n_Sf2.Eq.6)) Then
   Do i1=1,3
    If (Wpart) coup = coup - g2sq * Rsf1(j,i1) * Rsf2(l,i1)              &
           &                      * Conjg(Rsf1(i,i1) * Rsf2(k,i1))
    Do i2=1,3
     coup = coup                                                            &
        & - YL12 * Rsf1(j,i1) * Rsf2(l,i2) * Conjg(Rsf1(i,i1) * Rsf2(k,i2)) &
        & - Yr12 * Rsf1(j,i1+3) * Rsf2(l,i2+3)                              &
        &        * Conjg(Rsf1(i,i1+3) * Rsf2(k,i2+3))                       &
        & - YL1R2 * Rsf1(j,i1) * Rsf2(l,i2+3)                               &
        &        * Conjg(Rsf1(i,i1) * Rsf2(k,i2+3))                         &
        & - YR1L2 * Rsf1(j,i1+3) * Rsf2(l,i2)                               &
        &        * Conjg(Rsf1(i,i1+3) * Rsf2(k,i2))
    End Do
   End Do

  Else  If (ErrorLevel.Gt.-2) Then! there is a problem
   Write(ErrCan,*) "Error in subroutine "//NameOfUnit(Iname)
   Write(ErrCan,*) "number of sfermions n_sf1, n_sf2 =",n_sf1, n_sf2
   Write(ErrCan,*) "does not match any theory"
   If (ErrorLevel.Ge.0) Call TerminateProgram
  End If

  Iname = Iname - 1 

 End Subroutine CoupSfermion4G_13


 Subroutine CoupSfermion4Y_1(i, j, k, l, T3_1, Yuk1, Rsf1, T3_2, yuk2, Rsf2, nc &
                            &, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion self coupling for different flavours of
 ! sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i ........ index of the sfermion 1
 !  j ........ index of the adjungated sfermion 1
 !  k ........ index of the sfermion 2
 !  l ........ index of the adjungated sfermion 2
 !  Yuk1 ........ Yukawa coupling of sfermion1
 !  Rsf1(i,j) ... mixing matrix of the sfermions1
 !  Yuk2 ........ Yukawa coupling of sfermion2
 !  Rsf2(i,j) ... mixing matrix of the sfermions2
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 11.11.1999
 ! last change: 11.11.1999
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, nc
  Real(dp), Intent(in) :: T3_1, T3_2
  Complex(dp), Intent(in) :: Yuk1, Rsf1(2,2), Yuk2, Rsf2(2,2)
  Complex(dp), Intent(out) :: coup

  If ((T3_1*T3_2).Gt.0._dp) Then
   coup = - Yuk1 *Conjg(Yuk2 * Rsf1(i,1) * Rsf2(k,2)) * Rsf1(j,2) * Rsf2(l,1) &
        & - Yuk2 *Conjg(Yuk1 * Rsf2(k,1) * Rsf1(i,2)) * Rsf1(j,1) * Rsf2(l,2)
  Else
   coup = -Abs(Yuk1)**2 *Conjg(Rsf1(i,2) * Rsf2(k,1)) * Rsf1(j,2) * Rsf2(l,1) &
        & -Abs(Yuk2)**2 *Conjg(Rsf1(i,1) * Rsf2(k,2)) * Rsf1(j,1) * Rsf2(l,2)
  End If

!  coup = nc * coup

 End Subroutine CoupSfermion4Y_1


 Subroutine CoupSfermionSelf4G_13(i, j, k, l, gp, g, e, T3, Rsf, k_sf &
                                &, coup, self)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion gauge self coupling 
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,k ........ index of sfermions
 !  j,l ........ index of the adjungated sfermions
 !  Yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  self ....... integer, if present it is assumed that the first two
 !               indices belong to an electroweak eigenstate sfermion
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 01.01.02
 ! 04.02.02: correcting the case where the first two indices belong to
 !           an electroweak eigenstate
 !           The case where all states are mass eigenstates needs to be checked
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, k_sf
  Integer, Optional :: self
  Real(dp), Intent(in) :: gp, g, e, T3
  Complex(dp), Intent(in) :: Rsf(:,:)
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2, n_sf
  Real(dp) :: YL2, YR2, YLR2

  Iname = Iname + 1 
  NameOfUnit(Iname) = "CoupSfermionSelf4G_13"

  YL2 = 0.25_dp * g**2 + gp**2 * (T3 - e)**2
  YR2 = gp**2 * e**2
  YLR2 = gp**2 * e * (T3 - e)
  n_sf = Size( Rsf, dim=1 )

  coup = ZeroC

  If (n_sf.Eq.2) Then !no generation mixing
   If (Present(self)) Then
    If ((i.Eq.1).And.(j.Eq.1)) Then
     coup = - Yl2 * (1._dp + k_sf) * Rsf(l,1) * Conjg(Rsf(k,1))     &
          & - YLR2 * k_sf * Rsf(l,2) * Conjg(Rsf(k,2))
    Else If ((i.Eq.2).And.(j.Eq.2)) Then
     coup = - YR2 * (1._dp + k_sf) * Rsf(l,2) * Conjg(Rsf(k,2))     &
          & - YLR2 * k_sf * Rsf(l,1) * Conjg(Rsf(k,1))
    Else 
     coup = - YLR2 * Rsf(l,i) * Conjg(Rsf(k,j))
    End If
   Else
    coup = - Yl2 * Rsf(j,1) * Rsf(l,1) * Conjg(Rsf(i,1) * Rsf(k,1))    &
         & - YR2 * Rsf(j,2) * Rsf(l,2) * Conjg(Rsf(i,2) * Rsf(k,2))    &
         & - YLR2 * ( Rsf(j,1) * Rsf(l,2) * Conjg(Rsf(i,1) * Rsf(k,2)) &
         &          + Rsf(j,2) * Rsf(l,1) * Conjg(Rsf(i,2) * Rsf(k,1)) )
   End If

  Else If (n_sf.Eq.3) Then ! sneutrinos
   If (Present(self)) Then
    coup = Conjg(Rsf(k,j)) * Rsf(l,i)
    If ((i.Eq.j).And.(k.Eq.l)) coup = coup + 1

   Else
    Do i1=1,3
     Do i2=1,3
      coup = coup + Rsf(j,i1) * Rsf(l,i2) * Conjg(Rsf(i,i1) *Rsf(k,i2))  &
           &      + Rsf(j,i1) * Rsf(l,i2) * Conjg(Rsf(i,i2) *Rsf(k,i1))
     End Do
    End Do
   End If
   coup = - YL2 * coup 

  Else If (n_sf.Eq.6) Then ! sfermion with left right mixing
   If (Present(self)) Then
    If ((i.le.3).and.(j.le.3)) then ! LL
     coup = - Yl2 * Rsf(l,i) * Conjg(Rsf(k,j))
     If (i.eq.j) then
      Do i2=1,3
       coup = coup - k_sf * ( Yl2 * Rsf(l,i2) * Conjg(Rsf(k,i2))      &
           &               + YLR2 * Rsf(l,i2+3) * Conjg(Rsf(k,i2+3)) )
      End Do
     End If

    Else If ((i.ge.4).and.(j.ge.4)) then ! RR
     coup = - YR2 * Rsf(l,i) * Conjg(Rsf(k,j))
     If (i.eq.j) then
      Do i2=1,3
       coup = coup - k_sf * ( YLR2 * Rsf(l,i2) * Conjg(Rsf(k,i2))      &
           &               + YR2 * Rsf(l,i2+3) * Conjg(Rsf(k,i2+3)) )
      End Do
     End If

    Else ! LR or RL
     coup = - YLR2 * Rsf(l,i) * Conjg(Rsf(k,j))
    End If

   Else
    Do i1=1,3
     Do i2=1,3
       coup = coup                                                            &
         & - Yl2 * Rsf(j,i1) * Rsf(l,i2) * Conjg(Rsf(i,i1) *Rsf(k,i2))        &
         & - YR2 * Rsf(j,i1+3) * Rsf(l,i2+3) * Conjg(Rsf(i,i1+3)*Rsf(k,i2+3)) &
         & - YLR2 * ( Rsf(j,i1) * Rsf(l,i2+3) * Conjg(Rsf(i,i1)*Rsf(k,i2+3))  &
         &          + Rsf(j,i+31) * Rsf(l,i2) * Conjg(Rsf(i,i1+3)*Rsf(k,i2)) )
     End Do
    End Do
   End If

  Else If (ErrorLevel.Gt.-2) Then! there is a problem
   Write(ErrCan,*) "Error in subroutine "//NameOfUnit(Iname)
   Write(ErrCan,*) "number of sfermions n_sf =",n_sf
   Write(ErrCan,*) "does not match any theory"
   If (ErrorLevel.Ge.0) Call TerminateProgram
  End If
  
  coup = coup
  
  Iname = Iname - 1 

 End Subroutine CoupSfermionSelf4G_13

 Subroutine CoupSfermionSelf4Y_1(i, j, k, l, Yuk, Rsf, nc, coup, self)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion Yukawa self coupling 
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,k ........ index of sfermions
 !  j,l ........ index of the adjungated sfermions
 !  Yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  self ....... integer, if present it is assumed that the first two
 !               indices belong to an electroweak eigenstate sfermion
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 31.12.01
 ! 17.03.03: bug fixing due to SU(3) structure, assumed that colour
 !           index for i is the samke as j, and k has the same colour
 !           index as l 
 !     if (self==1) then we assume that the sfermions with indices i and j
 !     are electroweak eigenstates
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, nc
  Integer, Intent(in) :: self
  Complex(dp), Intent(in) :: Yuk, Rsf(2,2)
  Complex(dp), Intent(out) :: coup

  If (self.Eq.1) Then
   If ((i.eq.1).and.(j.eq.2)) then
    coup = - Abs(yuk)**2 * Conjg( RSf(k,2) ) * Rsf(l,1) * nc
   else if ((i.eq.2).and.(j.eq.1)) then
    coup = - Abs(yuk)**2 * Conjg( RSf(k,1) ) * Rsf(l,2) * nc
   else if ((i.eq.1).and.(j.eq.1)) then
    coup = - Abs(yuk)**2 * Conjg( RSf(k,2) ) * Rsf(l,2)
   else if ((i.eq.2).and.(j.eq.2)) then
    coup = - Abs(yuk)**2 * Conjg( RSf(k,1) ) * Rsf(l,1)
   end if
  Else
   coup = - Abs(Yuk)**2 * ( Conjg(Rsf(i,1) * Rsf(k,2)) * Rsf(j,2) * Rsf(l,1)  &
        &                 + Conjg(Rsf(k,1) * Rsf(i,2)) * Rsf(l,2) * Rsf(j,1) )
  End If

 End Subroutine CoupSfermionSelf4Y_1

 Subroutine AllCouplingsMSSM(g, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R             &
    & , Y_u, uU_L, uU_R, vevs                                                &
    & , RSpm, RP0, RS0, Umat, Vmat, Nmat, mu, phase_glu                      &
    & , RSlepton, Ae, Rsneut, RSup, Au, RSdown, Ad                           &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03       &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R       &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R    &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R   &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L           &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R           &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L         &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R           &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R           &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R           &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0         &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03 &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW    &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW        &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L    &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R             &
    & , cpl_SmpCN_L, cpl_SmpCN_R, c_GraDSd_L, c_GraDSd_R, c_GraUSu_L, c_GraUSu_R &
    & , c_GraLSl_L, c_GraLSl_R, c_GraNuSn_L, c_GraNuSn_R, GenerationMixing)
 !-----------------------------------------------------------------
 ! Routine for calculating all couplings of the MSSM
 ! output:
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 ! In addition, an _L or _R is added in the case of fermions indicated
 ! if it is the left-handed or right-handed coupling. 
 ! written by Werner Porod, 
 ! 01.03.2001: taking InitializeDecaysMSSM as basis
 !          - charged Higgs -sfermion-sfermion, 1 Gen., MSSM
 ! 12.09.2002: changing interface to avoid global variables
 !-------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: uL_L(3,3) ! mixing matrix of left leptons
  Complex(dp), Intent(in) :: uL_R(3,3) ! mixing matrix of right leptons
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3) ! mixing matrix of left d-quarks
  Complex(dp), Intent(in) :: uD_R(3,3) ! mixing matrix of right d-quarks
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3) ! mixing matrix of left u-quarks
  Complex(dp), Intent(in) :: uU_R(3,3) ! mixing matrix of right u-quarks
  Real(dp), Intent(in) :: vevS(2)     ! MSSM Higgs vevs [v_d, v_u]
  Complex(dp), Intent(in) :: RSpm(2,2) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(2,2)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(2,2)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(2,2), Vmat(2,2) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(4,4) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  mu   ! superpotential bilinear mu
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: RSlepton(6,6) ! slepton mixing matrix
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSneut(3,3)   ! sneutrino mixing matrix
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters
  Logical, Intent(in) :: GenerationMixing ! if .true. generation mixing of
                                          ! (s)fermions is taken into account

  Complex(dp), Intent(out) :: cpl_SmpSlSn(2,6,3), cpl_SmpSdSu(2,6,6)   &
      & , cpl_SmpSnSl(2,3,6), cpl_SmpSuSd(2,6,6), cpl_SmpP03(2,2,2)    &
      & , cpl_SmpP0W(2,2), cpl_SmpS03(2,2,2), cpl_SmpS0W(2,2)          &
      & , cpl_SmpLNu_L(2,3,3), cpl_SmpLNu_R(2,3,3), cpl_SmpDU_L(2,3,3) &
      & , cpl_SmpDU_R(8,3,3), cpl_SmpZ(2,2), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(2,2), cpl_CCZ_R(2,2)           &
      & , cpl_NNZ_L(4,4), cpl_NNZ_R(4,4), cpl_NNS0_L(4,4,2)            &
      & , cpl_NNS0_R(4,4,2), cpl_NNP0_L(4,4,2), cpl_NNP0_R(4,4,2) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,4,6), cpl_UNSu_R(3,4,6)         &
      & , cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6), cpl_NuNSn_L(3,4,3)      & 
      & , cpl_NuNSn_R(3,4,3), cpl_DDP0_L(3,3,2), cpl_LLP0_L(3,3,2)      &
      & , cpl_UUP0_L(3,3,2), cpl_DDP0_R(3,3,2), cpl_LLP0_R(3,3,2)       &
      & , cpl_UUP0_R(3,3,2), cpl_DDS0_L(3,3,2), cpl_LLS0_L(3,3,2)       &
      & , cpl_UUS0_L(3,3,2), cpl_DDS0_R(3,3,2), cpl_LLS0_R(3,3,2)       &
      & , cpl_UUS0_R(3,3,2)
  Complex(dp), Intent(out) :: cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6)      &
      & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)       &
      & , cpl_CLSn_R(2,3,3), cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(2)
  Complex(dp) :: cpl_P0SdSd(2,6,6), cpl_P0SuSu(2,6,6), cpl_P0SlSl(2,6,6) &
      & , cpl_P0SnSn(2,3,3), cpl_P0S0Z(2,2) 
  Real(dp), Intent(out) :: cpl_P0S03(2,2,2)
  Complex(dp), Intent(out) :: cpl_S0SdSd(2,6,6), cpl_S0SuSu(2,6,6) &
      & , cpl_S0SlSl(2,6,6), cpl_S0SnSn(2,3,3), cpl_LNuW(3,3)
  Real(dp), Intent(out) :: cpl_S03(2,2,2), cpl_S0WW(2), cpl_S0ZZ(2), cpl_FFpW
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(2,2,2), cpl_CCP0_R(2,2,2)    &
      & , cpl_CCS0_L(2,2,2), cpl_CCS0_R(2,2,2), cpl_CNW_L(2,4)        &
      & , cpl_CNW_R(2,4), cpl_SmpCN_L(2,2,4), cpl_SmpCN_R(2,2,4)
  Complex(dp), Intent(out), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R, c_GraLSl_L, c_GraLSl_R
  Complex(dp), Intent(out), Dimension(3,3) :: c_GraNuSn_L, c_GraNuSn_R
  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3, i4
  Real(dp) :: gU1, gSU2, gSU3, e_d, e_u, sinW2, cosW2, cosW
  Complex(dp) :: Rsd(2,2), Rsu(2,2), Rsl(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, mat6(6,6), mat3(3,3), bi(1)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsMSSM'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 2
  n_neut = 4
  n_P0 = 2
  n_S0 = 2
  n_Spm = 2

  !--------------------
  ! some constants
  !--------------------
  e_u = 2._dp / 3._dp
  e_d = - 1._dp / 3._dp

  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm     &
               &, gU1, gSU2, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpLNu_L = 0._dp
  cpl_SmpLNu_R = 0._dp
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp
  mat3 = 0._dp

  If (GenerationMixing) Then
   Do i1=1,2
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_l, uL_L, uL_R, mat3 &
              &, id3C, id3C, cpl_SmpLNu_L(i1,i2,i3), cpl_SmpLNu_R(i1,i2,i3) )
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_D, uD_L, uD_R, Y_U &
              &, uU_L, uU_R, cpl_SmpDU_L(i1,i2,i3), cpl_SmpDU_R(i1,i2,i3) )
     End Do
    End Do
   End Do
  Else
   Do i1=1,2
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_l(i2,i2), ZeroC               &
                       &, cpl_SmpLNu_L(i1,i2,i2), cpl_SmpLNu_R(i1,i2,i2) )
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  End If
  !------------------------------
  ! charged scalar - pseudoscalar - W
  !------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_Spm
    Do i3 = 1,n_P0
     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
                                       &, cpl_SmpP03(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !--------------------------------------
  ! charged scalar - charged scalar - scalar
  !--------------------------------------
   Do i1 = 1, n_Spm
    Do i2 = 1, n_Spm
     Do i3 = 1, n_S0
      Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
                                   &, cpl_SmpS03(i1,i2,i3) )
    End Do
   End Do
  End Do
  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSlSn = ZeroC
  cpl_SmpSnSl = ZeroC
  cpl_SmpSuSd = ZeroC
  cpl_SmpSdSu = ZeroC
  mat3 = zeroC

  If (GenerationMixing) Then
   Do i1 = 1, n_Spm
    Do i2 = 1,6
     Do i3 = 1,6
      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
                                  & , Y_d, Ad, Rsdown, Y_u, Au, Rsup, coupC)
      cpl_SmpSdSu(i1, i2, i3) = coupC
      cpl_SmpSuSd(i1, i3, i2) = Conjg(coupC)
     End Do
     Do i3 = 1,3
      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
                           & , Y_l, Ae, Rslepton, mat3, mat3, Rsneut, coupC)
      cpl_SmpSlSn(i1, i2, i3) = coupC
      cpl_SmpSnSl(i1, i3, i2) = Conjg(coupC)
     End Do
    End Do
   End Do

  Else
   Do i1 =1,n_Spm
    Do i2=1,3
     !----------
     ! Sleptons
     !----------
     Yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Call CoupChargedScalarSfermion3(i1, i3, 1, RSpm, gSU2, vevs, mu,    &
                       &             Yuk, A, Rsl, ZeroC, ZeroC, Id2C, coupC)
      cpl_SmpSlSn(i1,2*(i2-1)+i3,i2) = coupC
      cpl_SmpSnSl(i1,i2,2*(i2-1)+i3) = Conjg(coupC)
     End Do
     !----------
     ! Squarks
     !----------
     Yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Yukp = Y_u(i2,i2)
     Ap = Au(i2,i2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Do i4=1,2
       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
      End Do
     End Do
    End Do ! i2
   End Do ! i1
  End If

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, gSU2  &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, gSU2  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp
  cpl_CLSn_L = 0._dp
  cpl_CLSn_R = 0._dp
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  If (GenerationMixing) Then

   mat6 = 0._dp
   mat6(1:3,1:3) = Rsneut
   mat3 = 0._dp
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSlepton, Y_l &
                              &, mat3, id3c, id3c, Umat, Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i1, i2, i3) = coupLC
      cpl_CNuSl_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                              &, Y_u, uU_L, uU_R, Umat, Vmat, coupLC, coupRC)
      cpl_CUSd_L(i1, i2, i3) = coupLC
      cpl_CUSd_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                              &, Y_u, uD_L, uD_R, Umat, Vmat, coupLC, coupRC)
      cpl_CDSu_L(i1, i2, i3) = coupLC
      cpl_CDSu_R(i1, i2, i3) = coupRC
     End Do
     Do i3=1,3
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, mat6, Y_l, mat3 &
                              &, uL_L, uL_R, Umat, Vmat, coupLC, coupRC)
      cpl_CLSn_L(i1, i2, i3) = coupLC
      cpl_CLSn_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
   
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    YukP = 0._dp
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSl, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CNuSl_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2C, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
     cpl_CLSn_L(i2, i1, i1) = coupLC
     cpl_CLSn_R(i2, i1, i1) = coupRC
    End Do

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do
  End If 

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_LLP0_L = 0._dp
  cpl_LLP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_P0
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RP0 &
                                 &, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_l, uL_L, uL_R, RP0 &
                                 &, cpl_LLP0_L(i1,i2,i3), cpl_LLP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_d(i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_l(i1,i1), RP0   &
                                &, cpl_LLP0_L(i1,i1,i3), cpl_LLP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_LLS0_L = 0._dp
  cpl_LLS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_S0
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RS0 &
                           &, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_l, uL_L, uL_R, RS0 &
                           &, cpl_LLS0_L(i1,i2,i3), cpl_LLS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RS0 &
                           &, cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, -0.5_dp, Y_d(i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, -0.5_dp, Y_l(i1,i1), RS0   &
                          &, cpl_LLS0_L(i1,i1,i3), cpl_LLS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  Call CoupFermionZ(-0.5_dp, -1._dp, gSU2, sinW2, cpl_LLZ_L, cpl_LLZ_R)
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)
  Call CoupFermionZ(0.5_dp, 0._dp, gSU2, sinW2, cpl_NuNuZ_L, cpl_NuNuZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  cpl_FFpW = - gSU2 * oosqrt2
  cpl_LNuW = id3C * cpl_FFpW

  If (GenerationMixing) Then
   cpl_DUW = CKM * cpl_FFpW
  Else
   cpl_DUW = id3C * cpl_FFpW
  End If
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_LNSl_L = 0._dp
  cpl_LNSl_R = 0._dp
  cpl_NuNSn_L = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,4
     Do i3 = 1,6
      Call CoupNeutralinoSlepton(i1, i2, i3, gU1, gSU2, RSlepton, uL_L, uL_R &
                               &, Y_l, Nmat, coupLC, coupRC)
      cpl_LNSl_L(i1, i2, i3 ) = coupLC
      cpl_LNSl_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, Nmat, coupLC, coupRC)
      cpl_DNSd_L(i1, i2, i3 ) = coupLC
      cpl_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, Nmat, coupLC, coupRC)
      cpl_UNSu_L(i1, i2, i3 ) = coupLC
      cpl_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
     Do i3 = 1,3
      Call CoupNeutralinoSneutrino(i1, i2, i3, gU1, gSU2, Nmat, RSneut, id3C &
                               &, coupRC)
      cpl_NuNSn_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,4
    Call CoupNeutralinoSneutrino(i1, gU1, gSU2, Nmat, coupRC)
    Do i2 = 1,3
     cpl_NuNSn_R(i2, i1, i2 ) = coupRC
    End Do
   End Do

   Do i1 = 1,3
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    Do i2 = 1,2
     Do i3 = 1,4
      Call CoupNeutralinoSlepton(i3, i2, gU1, gSU2, RSl, Yuk, Nmat &
                               &, coupLC, coupRC)
      cpl_LNSl_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_LNSl_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,4
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If
  !--------------------------------
  ! Gravitino - fermion - sfermion
  !--------------------------------
  c_GraNuSn_L = 0._dp
  c_GraNuSn_R = 0._dp
  c_GraLSl_L = 0._dp
  c_GraLSl_R = 0._dp
  c_GraDSd_L = 0._dp
  c_GraDSd_R = 0._dp
  c_GraUSu_L = 0._dp
  c_GraUSu_R = 0._dp

  If (GenerationMixing) Then

   Do i1=1,3
    Do i2=1,3
     Call CoupGravitinoSfermion(i1, i2, RSneut, RSneut, id3c &
                  & , c_GraNuSn_L(i1,i2), c_GraNuSn_R(i1,i2))
    End Do
    Do i2=1,6
     Call CoupGravitinoSfermion(i1, i2, RSlepton, uL_L, uL_R &
                  & , c_GraLSl_L(i1,i2), c_GraLSl_R(i1,i2))
     Call CoupGravitinoSfermion(i1, i2, RSdown, uD_L, uD_R &
                  & , c_GraDSd_L(i1,i2), c_GraDSd_R(i1,i2))
     Call CoupGravitinoSfermion(i1, i2, RSup, uU_L, uU_R &
                  & , c_GraUSu_L(i1,i2), c_GraUSu_R(i1,i2))
    End Do
   End Do

  Else
   
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    c_GraNuSn_R(i1, i1) = 1._dp
    Do i2=1,2
     Call CoupGravitinoSfermion(i2, RSd, coupLC, coupRC)
     c_GraDSd_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGravitinoSfermion(i2, RSu, coupLC, coupRC)
     c_GraUSu_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraUSu_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGravitinoSfermion(i2, RSl, coupLC, coupRC)
     c_GraLSl_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraLSl_R(i1, (i1-1)*2 + i2) = coupRC
    End Do
   End Do

  End If 

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsdown, uD_L, uD_R, &
                          & coupLC, coupRC)
     cpl_GDSd_L(i1, i2) = coupLC
     cpl_GDSd_R(i1, i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsup, uU_L, uU_R, &
                          & coupLC, coupRC)
     cpl_GUSu_L(i1, i2) = coupLC
     cpl_GUSu_R(i1, i2) = coupRC
    End Do
   End Do
  Else

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  End If 

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, Nmat, RP0, gU1, gSU2, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalar(i1, i2, i3, Nmat, RS0, gU1, gSU2, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  cpl_SlSnW = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,6
    Do i2 = 1,6
     Call CoupSfermionW3(i1, i2, gSU2, RSdown, RSup, cpl_SdSuW(i1,i2) )
    End Do
    Do i2 = 1,3
     Call CoupSfermionW3(i1, i2, gSU2, RSlepton, RSneut, cpl_SlSnW(i1,i2) )
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do

    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Call CoupSfermionW3(i2, 1, gSU2, RSl, id2C, coupC )
      cpl_SlSnW( (i1-1)*2+i2, i1) = coupC
    End Do
   End Do

  End If

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)
  Call Adjungate(cpl_SlSnW, cpl_SnSlW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
  Do i1=1,n_P0
   Do i2=1,n_P0
    Do i3=1,n_S0
     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
                                &, cpl_P0S03(i1,i2,i3) )
     End Do
    End Do
   End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
    End Do
   End Do
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp
  cpl_P0SlSl = 0._dp
  cpl_P0SnSn = 0._dp

  bi(1) = mu

  If (GenerationMixing) Then
   Do i1=1,n_P0
    Do i2=1,6
     Do i3=1,6
      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_d, Rsdown   &
                                   &, Ad, bi, cpl_P0SdSd(i1,i2,i3) )
      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, 0.5_dp, Y_u, Rsup      &
                                   &, Au, bi, cpl_P0SuSu(i1,i2,i3) )
      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_l, Rslepton &
                                   &, Ae, bi, cpl_P0SlSl(i1,i2,i3) )
     End Do
    End Do
   End Do

  Else
   Do i1=1,n_P0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsl   &
                                   &, A, bi, coupC )
       cpl_P0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
                                   &, A, bi, coupC )
       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
                                   &, A, bi, coupC )
       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  End If

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
  Do i1=1,n_S0
   Do i2=1,n_S0
    Do i3=1,n_S0
     Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
    End Do
   End Do
  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp
  cpl_S0SlSl = 0._dp
  cpl_S0SnSn = 0._dp
  mat3 = zeroC

  If (GenerationMixing) Then
   Do i1=1,n_S0
    Do i2=1,6
     Do i3=1,6
      Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d, Rsdown   &
                           &, Ad, mu, vevs, gU1, gSU2, cpl_S0SdSd(i1,i2,i3) )
      Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u, Rsup      &
                           &, Au, mu, vevs, gU1, gSU2, cpl_S0SuSu(i1,i2,i3) )
      Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, -1._dp, Y_l        &
                 &, Rslepton, Ae, mu, vevs, gU1, gSU2, cpl_S0SlSl(i1,i2,i3) )
     End Do
    End Do
    Do i2=1,3
     Do i3=1,3
      Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, 0._dp, mat3, Rsneut &
                         &, mat3, mu, vevs, gU1, gSU2, cpl_S0SnSn(i1,i2,i3) )
     End Do
    End Do
   End Do

  Else
   Do i1=1,n_S0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

     Call CoupScalarSfermion3(i1, 1, 1, RS0, 0.5_dp, 0._dp, ZeroC, id3C  &
                               &, ZeroC, mu, vevs, gU1, gSU2, coupC )
     cpl_S0SnSn(i1,i2,i2) = coupC

     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, -1._dp, Yuk, Rsl  &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
                               &, A, mu, vevs, gU1, gSU2, coupC )
       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  End If

  !----------
  ! scalar W
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarW(i1, gSU2, vevs, RS0, cpl_S0WW(i1) )
  End Do

  !----------
  ! scalar Z
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevs, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SlSlZ = 0._dp
  cpl_SnSnZ = 0._dp
  cpl_SuSuZ = 0._dp

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,3
     Call CoupSneutrinoZ(i1, i2, gSU2, sinW2, Rsneut, cpl_SnSnZ(i1, i2) )
    End Do
   End Do 
   Do i1=1,6
    Do i2=1,6
     Call CoupSleptonZ(i1, i2, gSU2, sinW2, Rslepton, cpl_SlSlZ(i1, i2) )
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, cpl_SdSdZ(i1, i2) )
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, cpl_SuSuZ(i1, i2) )
    End Do
   End Do
 
  Else
   Call CoupSneutrinoZ(gSU2, sinW2, cpl_SnSnZ(1,1) )
   Do i1=1,3
    cpl_SnSnZ(i1,i1) = cpl_SnSnZ(1,1)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSleptonZ(i2, i3, gSU2, sinW2, Rsl, coupC )
      cpl_SlSlZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  End If

  Iname = Iname - 1

 End Subroutine AllCouplingsMSSM

 Subroutine CoupChargedScalarSfermion3MSSM3(i, j, k, RSpm, g, vevs, mu    &
                                  & , yukd, Ad, Rsfd, yuku, Au, Rsfu, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between charged scalar and sfermions
 ! valid for the 3-gneration MSSM
 !  i .......... index of charged scalar boson, positiv
 !  j .......... index of the down-type sfermion \tilde{f}
 !  k .......... index of the up-type sfermion \conjugate{\tilde{f}}
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vevs of the Higgs fields
 !  mu ......... Higgs mass parameter of the superpotential
 !  yukd(i,j) .. Yukawa coupling of the down-type sfermion
 !  Ad(i,j) .... trilinear coupling of the down-type sfermion
 !  Rsfd(i,j) .. mixing matrix of the down-type sfermion
 !  yuku(i,j) .. Yukawa coupling of the up-type sfermion
 !  Au(i,j) .... trilinear coupling of the up-type sfermion
 !  Rsfu(i,j) .. mixing matrix of the down-type sfermion
 ! output
 !  coup ....... the coupling(i,j,k)
 ! the lagrangian is given by
 ! coup \tilde{d}(j) \tilde{u}^*(k) H^+(i)
 ! written by Werner Porod, 23.04.2001
 ! 17.03.03: fixing bug in the right-right part, where the transposed
 !           Yukawa had been used
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(2,2), RSfd(6,6),RSfu(:,:) &
                           &   , yukd(3,3), yuku(3,3), Ad(3,3), Au(3,3), mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: g, vevs(2)
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2, i3, n_sfer
  Complex(dp) :: sum1, RSda(6,6), Yukda(3,3), Yukua(3,3), muC   &
                 & , Rsp1, Rsp2C, sum2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion3MSSM3'

  n_sfer = Size(RSfu, Dim=1)

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: i = ',i
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.6)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'down-type sfermion index j out of range: j = ',j
   Call TerminateProgram
  End If
  If ((k.Lt.1).Or.(k.Gt.n_sfer)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'up-type sfermion index k out of range: k = ',k
   Call TerminateProgram
  End If

  RSda = Conjg( RSfd )
  Yukda = Conjg( Yukd )
  Yukua = Conjg( Yuku )
  muC = Conjg( mu )
  Rsp1 = RSpm(i,1)
  Rsp2C = Conjg( RSpm(i,2) )
  sum1 = oosqrt2 * 0.5_dp * g**2 * (vevs(1) * RSp1 + vevs(2) * RSp2C )
  sum2 = oosqrt2 * (vevs(2) * RSp1 + vevs(1) * RSp2C )

  coup = ZeroC

  Do i1 = 1,3
   !-----------------
   ! left d left u
   !-----------------
   coup = coup - sum1 * RSfu(k,i1) * RSda(j,i1)
   Do i2 = 1,3
    Do i3 = 1,3
     coup = coup &
        & + oosqrt2 * RSfu(k,i3) * RSda(j,i1)                        &
        &           * ( yukD(i1,i2) * yukDa(i3,i2) * vevs(1) * RSp1  &
        &             + yukU(i1,i2) * yukUa(i3,i2) * vevs(2) * RSp2C ) 
    End Do
   !-----------------
   ! right d left u
   !-----------------
    coup = coup + ( mu * YukDa(i1,i2) * Rsp2C + Conjg( Ad(i1,i2) ) * RSp1) &
         &        * RSfu(k,i1) * RSda(j,i2+3)
   End Do
   If (n_sfer.Eq.6) Then
    !-----------------
    ! left d right u
    !-----------------
    Do i2 = 1,3
     coup = coup + ( muC * YukU(i1,i2) * RSp1 + Au(i1,i2) * RSp2C )   &
          &        * RSfu(k,3+i2) * RSda(j,i1)
    End Do
    !-----------------
    ! right d right u
    !-----------------
    Do i2 = 1,3
     Do i3 = 1,3
     coup = coup  &
        & + sum2 * YukU(i2,i1) * YukDa(i2,i3) * RSfu(k,3+i1) * RSda(j,3+i3)
     End Do
    End Do
   End If
  End Do

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion3MSSM3

 Subroutine CoupChargedScalarSfermion4MSSM3(i, j, k, l, RSpm, T3, e, yuk  &
                                     &, yukb, Rsf, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between charged scalar and sfermions
 ! valid for the 3-generation MSSM
 !  i,j ........ index of charged scalar boson (i=+, j=- )
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk(i,j) ... Yukawa coupling of the sfermion
 !  yukb(i,j) .. Yukawa coupling of the other sfermion in the SU(2) doublett 
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod,  09.11.01
 !  - sneutrinos need special treatment
 ! 17.03.03: bug fix: Yukawas in parts(2) had both been complex conjugated
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(2,2), RSf(:,:), yuk(3,3), yukb(3,3)
  Real(dp), Intent(in) :: g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2, i3, i_sf
  Complex(dp) :: parts(2), Dterm
  Real(dp) :: g2, gp2, YL, YR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion4MSSM3'

  coup = Cmplx( 0._dp,0._dp,dp )

  i_sf = Size(Rsf, dim=1)

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = Cmplx( 0._dp,0._dp,dp )

  Dterm = 0._dp

  If (i_sf.Eq.6) Then 
   Do i1=1,3
    Dterm = Dterm                                                             &
        & + 0.5_dp * ( (T3 * g2 + YL * gp2) * Rsf(l,i1) * Conjg( Rsf(k,i1) )  &
        &            + YR * gp2 * Rsf(l,i1+3) * Conjg( Rsf(k,i1+3) ) )
   End Do
  Else ! Sneutrinos
   Dterm = 0._dp
   if (l.eq.k)  Dterm = 0.5_dp * (T3 * g2 + YL * gp2)
  End If

  If (T3.Gt.0._dp) Then   
   parts(1) = Dterm
   parts(2) = - Dterm
   If (i_sf.Eq.6) Then
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       parts(1) = parts(1) - Conjg(Rsf(k,i1)) * yukb(i1,i2)         &
                &            * Conjg(yukb(i3,i2)) * Rsf(l,i3)
       parts(2) = parts(2) - Conjg(Rsf(k,i1+3)) * Conjg(yuk(i2,i1))  &
                &            * yuk(i2,i3) * Rsf(l,i3+3)
      End Do
     End Do
    End Do
   Else ! sneutrinos
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       parts(1) = parts(1) - Conjg(Rsf(k,i1)) * yukb(i1,i2)         &
                &            * Conjg(yukb(i3,i2)) * Rsf(l,i3)
      End Do
     End Do
    End Do
   End If

  Else ! T_3 < 0
   parts(1) = Dterm
   parts(2) = - Dterm
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      parts(1) = parts(1) - Conjg(Rsf(k,i1+3)) * Conjg(yuk(i2,i1))  &
               &            * yuk(i2,i3) * Rsf(l,i3+3)
      parts(2) = parts(2) - Conjg(Rsf(k,i1)) * yukb(i1,i2)         &
               &            * Conjg(yukb(i3,i2)) * Rsf(l,i3)
     End Do
    End Do
   End Do

  Endif

  Do i1=1,2
   coup = coup + RSpm(i,i1) * Conjg( RSpm(j,i1) ) * parts(i1)
  Enddo

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion4MSSM3

 Subroutine CoupChargedScalarSfermion3MSSM1(i, j, k, RSpm, g, vevs, mu    &
                                  & , yukd, Ad, Rsfd, yuku, Au, Rsfu, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between charged scalar and sfermions
 ! valid for the MSSM
 !  i .......... index of charged scalar boson, positiv
 !  j .......... index of the down-type sfermion \tilde{f}
 !  k .......... index of the up-type sfermion \conjugate{\tilde{f}}
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vevs of the Higgs fields
 !  mu ......... Higgs mass parameter of the superpotential
 !  yukd ....... Yukawa coupling of the down-type sfermion
 !  Ad ......... trilinear coupling of the down-type sfermion
 !  Rsfd(i,j) .. mixing matrix of the down-type sfermion
 !  yuku ....... Yukawa coupling of the up-type sfermion
 !  Au ......... trilinear coupling of the up-type sfermion
 !  Rsfu(i,j) .. mixing matrix of the down-type sfermion
 ! output
 !  coup ....... the coupling(i,j,k)
 ! the lagrangian is given by
 ! coup \tilde{d}(j) \tilde{u}^*(k) H^+(i)
 ! written by Werner Porod, 15.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(2,2), RSfd(2,2), RSfu(2,2) &
                           &   ,yukd, yuku, Ad, Au, mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: g, vevs(2)
  Integer :: i, j, k

  Integer :: i1, i2
  Complex(dp) :: parts(2,2)
  Real(dp) :: g2, yukd2, yuku2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion3MSSM1'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: i = ',i
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'down-type sfermion index j out of range: j = ',j
   Call TerminateProgram
  End If
  If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'up-type sfermion index k out of range: k = ',k
   Call TerminateProgram
  End If

  coup = ZeroC
  parts = ZeroC

  g2 = g**2
  yukd2 = Abs( yukd )**2
  yuku2 = Abs( yuku )**2

  parts(1,1) = - ( vevs(1) * Rspm(i,1) * (0.5_dp*g2 - yukd2) &
           &     + vevs(2) * Rspm(i,2) * (0.5_dp*g2 - yuku2) ) * oosqrt2
 
  parts(1,2) = yuku * Conjg( mu ) * Rspm(i,1) + Au * Rspm(i,2) 

  parts(2,1) = Conjg(Ad) * RSpm(i,1) + Conjg( yukd ) * mu * RSpm(i,2)

  parts(2,2) = yuku * Conjg(yukd) * (vevs(1) *Rspm(i,2) + vevs(2) *RSpm(i,1) )
  parts(2,2) = oosqrt2 * parts(2,2)

  Do i1=1,2
   Do i2=1,2
    coup = coup + Conjg( Rsfd(j,i1) ) * parts(i1,i2) * Rsfu(k,i2)
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion3MSSM1

 Subroutine CoupChargedScalarSfermion4MSSM1(i, j, k, l, RSpm, T3, e, yuk  &
                                     &, yukb, Rsf, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between charged scalar and sfermions
 ! valid for the 1-generation MSSM
 !  i,j ........ index of charged scalar boson (i=+, j=- )
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk ........ Yukawa coupling of the sfermion
 !  yukb ....... Yukawa coupling of the other sfermion in the SU(2) doublett 
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 29.11.1999
 ! 09.11.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(2,2), RSf(2,2), yuk, yukb
  Real(dp), Intent(in) :: g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: i1
  Complex(dp) :: parts(2), Dterm
  Real(dp) :: g2, gp2, YL, YR, yuk2, yukb2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion4MSSM1'

  coup = Cmplx( 0._dp,0._dp,dp )

  yuk2 = Abs( yuk )**2
  yukb2 = Abs( yukb )**2

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = Cmplx( 0._dp,0._dp,dp )

  Dterm = 0.5_dp * ( (T3 * g2 + YL * gp2) * Rsf(l,1) * Conjg( Rsf(k,1) )  &
        &          + YR * gp2 * Rsf(l,2) * Conjg( Rsf(k,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = Dterm - yukb2 * Rsf(l,1) * Conjg( Rsf(k,1) )
   parts(2) = - Dterm - yuk2 * Rsf(l,2) * Conjg( Rsf(k,2) )

  Else
   parts(1) = Dterm - yuk2 * Rsf(l,2) * Conjg( Rsf(k,2) )
   parts(2) = - Dterm - yukb2 * Rsf(l,1) * Conjg( Rsf(k,1) )

  Endif

  Do i1=1,2
   coup = coup + RSpm(i,i1) * Conjg( RSpm(j,i1) ) * parts(i1)
  Enddo

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion4MSSM1
 Subroutine CoupCharginoSfermion3(i, j, k, g, T3, RSf, YukD, YukU, RfL, RfR &
                                &, U, V, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between chargino, fermion, and sfermion'
 ! valid for the 3-gernation MSSM
 !  i .......... index of chargino
 !  j .......... index of fermion
 !  k .......... index of sfermion
 !  g .......... SU(2) gauge coupling
 !  T3 ......... weak isospin of the fermion
 !  RSf(i,j) ... mixing matrix of the sfermion
 !  YukD ....... Yukawa coupling of the T3=-0.5 fermion 
 !  YukU ....... Yukawa coupling of the T3=0.5 fermion 
 !  U,V ........ mixing matrices of the chargino
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding part of the Lagrangian is:
 !  \bar{f_j} (coupL P_L + coupR P_R) chi^\pm_i \tilde{f'}_k
 ! written by Werner Porod, 16.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, T3
  Complex(dp), Intent(in) :: RSf(:,:), U(:,:), V(:,:)
  Complex(dp), Intent(in), Dimension(3,3) :: YukD, YukU, RfL, RfR
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Integer :: n_char, i1, i2, n_sf

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoSfermion3' 

  n_char = Size(U, Dim=1)
  n_sf = Size(Rsf, Dim=1)

  If (ErrorLevel.Ge.-1) Then
   If ((i.Lt.1).Or.(i.Gt.n_char)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index i out of range: (i,n_char) = ',i,n_char
    Call TerminateProgram
   End If 
   If ((j.Lt.1).Or.(j.Gt.3)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index j out of range: ',j
    Call TerminateProgram
   End If 
   If ((k.Lt.1).Or.(k.Gt.6)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index k out of range: ',k
    Call TerminateProgram
   End If 
  End If

  coupL = 0._dp
  coupR = 0._dp

  If ( T3.Gt.0._dp) Then

   If ( CompareMatrices(RfL,id3C,NearlyZero) .And. &
      & CompareMatrices(RfR,id3C,NearlyZero)        ) Then
    coupR = - g  * Conjg( RSf(k,j) ) * U(i,1)
    Do i1=1,3
     coupL = coupL + YukU(i1,j) * Conjg( RSf(k,i1) )
     coupR = coupR + Conjg( YukD(j,i1) * RSf(k,i1+3) ) * U(i,2)
    End Do

   Else If (CompareMatrices(RfL,id3C,NearlyZero)) Then
    coupR = - g  * Conjg( RSf(k,j) ) * U(i,1)
    Do i1=1,3
     coupR = coupR + Conjg( YukD(j,i1) * RSf(k,i1+3) ) * U(i,2)
     Do i2=1,3
      coupL = coupL + YukU(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
     End Do
    End Do

   Else If (CompareMatrices(RfR,id3C,NearlyZero)) Then
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * U(i,1)
     coupL = coupL + YukU(i1,j) * Conjg( RSf(k,i1) )
     Do i2=1,3
      coupR = coupR + Conjg( YukD(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * U(i,2)
     End Do
    End Do

   Else
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * U(i,1)
     Do i2=1,3
      coupL = coupL + YukU(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
      coupR = coupR + Conjg( YukD(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * U(i,2)
     End Do
    End Do

   End If
   coupL = coupL * Conjg( V(i,2) )

  Else If (n_sf.Eq.3) Then ! T3 < 0, and sneutrinos 
 
   If ( CompareMatrices(RfL,id3C,NearlyZero) .And. &
      & CompareMatrices(RfR,id3C,NearlyZero)       ) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) )
    End Do

   Else If (CompareMatrices(RfL,id3C,NearlyZero)) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     Do i2=1,3
      coupL = coupL + YukD(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
     End Do
    End Do

   Else If (CompareMatrices(RfR,id3C,NearlyZero)) Then
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) )
    End Do

   Else
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     Do i2=1,3
      coupL = coupL + YukD(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
     End Do
    End Do
   End If
   coupL = coupL * Conjg( U(i,2) )

  Else  ! T3 < 0
 
   If ( CompareMatrices(RfL,id3C,NearlyZero) .And. &
      & CompareMatrices(RfR,id3C,NearlyZero)       ) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) )
     coupR = coupR + Conjg( YukU(j,i1) * RSf(k,i1+3) ) * V(i,2)
    End Do

   Else If (CompareMatrices(RfL,id3C,NearlyZero)) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     coupR = coupR + Conjg( YukU(j,i1) * RSf(k,i1+3) ) * V(i,2)
     Do i2=1,3
      coupL = coupL + YukD(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
     End Do
    End Do

   Else If (CompareMatrices(RfR,id3C,NearlyZero)) Then
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) )
     Do i2=1,3
      coupR = coupR + Conjg( YukU(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * V(i,2)
     End Do
    End Do

   Else
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     Do i2=1,3
      coupL = coupL + YukD(i1,i2) * Conjg( RSf(k,i1)  * RfR(j,i2) )
      coupR = coupR + Conjg( YukU(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * V(i,2)
     End Do
    End Do
   End If
   coupL = coupL * Conjg( U(i,2) )

  End If

  Iname = Iname - 1

 End Subroutine CoupCharginoSfermion3

 Subroutine CoupFermionPseudoScalar3(i, j, k, T3, yuk, Rfl, Rfr, RP0 &
                                   &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between fermions and the pseudoscalars
 ! valid for 3-generation MSSM
 !  i,j ........ generation index of the fermions
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) P0(k)
 ! written by Werner Porod, 27.04.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3, RP0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3), RFl(3,3), RFr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2

  coupL = ZeroC
  coupR = ZeroC

  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk(j,i)
     coupR = conjg(yuk(i,j))

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   Do i2=1,3
    coupL = coupL + yuk(j,i2) * Conjg(RFr(i,i2))
    coupR = coupR + Conjg(yuk(i,i2)) * RFr(j,i2)
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
    coupR = coupR + RFl(i,i1) * Conjg(yuk(i1,j))
   End Do

  Else
   Do i1=1,3
    Do i2=1,3
     coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i2) * Conjg(RFr(i,i2))
     coupR = coupR + RFl(i,i1) * Conjg(yuk(i1,i2)) * RFr(j,i2)
    End Do
   End Do
  End If

  If (T3.Gt.0._dp) Then
   coupL = - (0._dp,1._dp) * coupL * RP0(k,2) * oosqrt2
   coupR = (0._dp,1._dp) * coupR * RP0(k,2) * oosqrt2
  Else
   coupL = - (0._dp,1._dp) * coupL * RP0(k,1) * oosqrt2
   coupR = (0._dp,1._dp) * coupR * RP0(k,1) * oosqrt2
  End If

 End Subroutine CoupFermionPseudoScalar3

 Subroutine CoupFermionPseudoScalar3a(i, j, k, T3, yuk, Rfl, Rfr, RP0 &
                                   &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between fermions and the pseudoscalars
 ! valid for 3-generation MSSM
 !  i .......... generation index of the fermions (electroweak eigenstate)
 !  j .......... generation index of the fermions (mass eigenstate)
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) P0(k)
 ! written by Werner Porod, 15.06.2008
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3, RP0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3), RFl(3,3), RFr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2

  coupL = ZeroC
  coupR = ZeroC

  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk(j,i)
     coupR = conjg(yuk(i,j))

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   coupL = yuk(j,i)
   Do i2=1,3
    coupR = coupR + Conjg(yuk(i,i2)) * RFr(j,i2)
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   coupR = conjg(yuk(i,j))
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
   End Do

  Else
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
    coupR = coupR + Conjg(yuk(i,i1)) * RFr(j,i1)
   End Do
  End If

  If (T3.Gt.0._dp) Then
   coupL = - (0._dp,1._dp) * coupL * RP0(k,2) * oosqrt2
   coupR = (0._dp,1._dp) * coupR * RP0(k,2) * oosqrt2
  Else
   coupL = - (0._dp,1._dp) * coupL * RP0(k,1) * oosqrt2
   coupR = (0._dp,1._dp) * coupR * RP0(k,1) * oosqrt2
  End If

 End Subroutine CoupFermionPseudoScalar3a

 Subroutine CoupFermionScalar3(i, j, k, T3, yuk, Rfl, Rfr, RS0 &
                                   &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between fermions and the scalars
 ! valid for 3-generation MSSM
 !  i,j ........ generation index of the fermions
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RS0(i,j) ... mixing matrix of calar bosons
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) S0(k)
 ! written by Werner Porod, 27.04.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3, RS0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3), RFl(3,3), RFr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2

  coupL = ZeroC
  coupR = ZeroC

  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk(j,i)
     coupR = Conjg(yuk(i,j)) 

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   Do i2=1,3
    coupL = coupL + yuk(j,i2) * Conjg(RFr(i,i2))
    coupR = coupR + Conjg(yuk(i,i2)) * RFr(j,i2)    
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
    coupR = coupR + RFl(i,i1) * conjg(yuk(i1,j))
   End Do

  Else
   Do i1=1,3
    Do i2=1,3
     coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i2) * Conjg(RFr(i,i2))
     coupR = coupR + RFl(i,i1) * conjg(yuk(i1,i2)) * RFr(j,i2)
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

 End Subroutine CoupFermionScalar3

 Subroutine CoupFermionScalar3a(i, j, k, T3, yuk, Rfl, Rfr, RS0 &
                                   &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between fermions and the scalars
 ! valid for 3-generation MSSM
 !  i .......... generation index of the fermions (electroweak eigenstate)
 !  j .......... generation index of the fermions (mass eigenstate)
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  RFl ........ mixing matrix for the left-handed fermions
 !  RFr ........ mixing matrix for the right-handed fermions
 !  RS0(i,j) ... mixing matrix of calar bosons
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by
 !  \bar{f}(i) (coupL P_L + coupR P_R) f(j) S0(k)
 ! written by Werner Porod, 27.04.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: T3, RS0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3), RFl(3,3), RFr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k
 
  Integer :: i1, i2

  coupL = ZeroC
  coupR = ZeroC

  If ( CompareMatrices(id3C,RFl,NearlyZero).And. &
     & CompareMatrices(id3C,RFr,NearlyZero)      ) Then
     coupL = yuk(j,i)
     coupR = Conjg(yuk(i,j)) 

  Else If (CompareMatrices(id3C,RFl,NearlyZero) ) Then
   coupL = yuk(j,i)
   Do i2=1,3
    coupR = coupR + Conjg(yuk(i,i2)) * RFr(j,i2)    
   End Do

  Else If (CompareMatrices(id3C,RFr,NearlyZero) ) Then
   coupR = Conjg(yuk(i,j)) 
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
   End Do

  Else
   Do i1=1,3
    coupL = coupL + Conjg(RFl(j,i1)) * yuk(i1,i)
    coupR = coupR + Conjg(yuk(i,i1)) * RFr(j,i1)
   End Do
  End If

  If (T3.Gt.0._dp) Then
   coupL = - coupL * RS0(k,2) * oosqrt2
   coupR = - coupR * RS0(k,2) * oosqrt2
  Else
   coupL = - coupL * RS0(k,1) * oosqrt2
   coupR = - coupR * RS0(k,1) * oosqrt2
  End If

 End Subroutine CoupFermionScalar3a


 Subroutine CoupGluinoSquark3(g3, phase_M3, i, j, Rsq, RqL, RqR, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling gluino-squark-quark without generation mixing
 ! and without the generator T^a_ij, the Lagrangian is given by:
 !  L = g3 T^a_ij \bar{q} (coupR P_R + coupL P_L) g^a \tilde q_(i_sq)
 ! input:
 !  g3 ............ SU(3) coupling
 !  phase_M3 ...... phase of the parameter M3
 !  i ............. index of the quark
 !  j ............. index of the squark
 !  Rsq(i,j) ...... squark mixing matrix
 !  RqL(i,j) ...... left quark mixing matrix
 !  RqR(i,j) ...... right quark mixing matrix
 ! output
 !  coupL ......... left coupling
 !  coupR ......... right coupling
 ! the lagrangian is given by
 !  Lambda^a \bar{q_i} ( coupL P_L + coupR P_R) g^a \tilde{q}_j
 ! where Lambda^a are the Gell-Mann matrices
 ! written by Werner Porod, 21.1.01
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g3
  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: phase_M3, Rsq(6,6), RqL(3,3), RqR(3,3)
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: i1
  Complex(dp) :: sqrt_phase

  sqrt_phase = Sqrt( phase_M3 )
  coupL = 0._dp
  coupR = 0._dp

  If ( CompareMatrices(RqR,id3C,NearlyZero).And. &
     & CompareMatrices(RqL,id3C,NearlyZero)      ) Then
   coupL = Conjg( Rsq(j,i+3) )
   coupR = Conjg( Rsq(j,i) )

  Else If (CompareMatrices(RqR,id3C,NearlyZero)) Then
   coupL = Conjg( Rsq(j,i+3) )
   Do i1 = 1,3
    coupR = coupR + Conjg( Rsq(j,i1) ) * RqL(i,i1)
   End Do

  Else If (CompareMatrices(RqL,id3C,NearlyZero)) Then
   coupR = Conjg( Rsq(j,i) )
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsq(j,i1+3) * RqR(i,i1) )
   End Do

  Else
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsq(j,i1+3) * RqR(i,i1) )
    coupR = coupR + Conjg( Rsq(j,i1) ) * RqL(i,i1)
   End Do
  End If

  coupL = g3 * coupL * oosqrt2 / sqrt_phase
  coupR = - g3 * oosqrt2 * sqrt_phase * coupR

 End Subroutine CoupGluinoSquark3

 Subroutine CoupGravitinoSfermion3(i, j, RSf, RL, RR, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the overall coupling gravitino-sfermion-fermion 
 ! and without the generator T^a_ij, the Lagrangian is given by:
 !  L = \bar{f_i} (coupR P_R + coupL P_L) gamma_mu gamma^nu G^mu D_\nu 
 !            *  \tilde f_j / (Sqrt[2] m_Planck)
 ! input:
 !  g3 ............ SU(3) coupling
 !  phase_M3 ...... phase of the parameter M3
 !  i ............. index of the fermion
 !  j ............. index of the sfermion
 !  RSf(i,j) ...... sfermion mixing matrix
 !  RL(i,j) ....... left fermion mixing matrix
 !  RR(i,j). ...... right fermion mixing matrix
 ! output
 !  coupL ......... left coupling
 !  coupR ......... right coupling
 ! written by Werner Porod, 18.09.10
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Complex(dp), Intent(in) :: RSf(:,:), RL(:,:), RR(:,:)
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: i1, len

  coupL = 0._dp
  coupR = 0._dp

  len = Size( Rsf, dim=1 )

  If ( CompareMatrices(RR,id3C,NearlyZero).And. &
     & CompareMatrices(RL,id3C,NearlyZero)      ) Then
   If (len.Eq.6) coupL = Conjg( RSf(j,i+3) )
   coupR = Conjg( RSf(j,i) )

  Else If (CompareMatrices(RR,id3C,NearlyZero)) Then
   If (len.Eq.6) coupL = Conjg( RSf(j,i+3) )
   Do i1 = 1,3
    coupR = coupR + Conjg( RSf(j,i1) ) * RL(i,i1)
   End Do

  Else If (CompareMatrices(RL,id3C,NearlyZero)) Then
   coupR = Conjg( RSf(j,i) )
   If (len.Eq.6) Then
    Do i1 = 1,3
     coupL = coupL + Conjg( RSf(j,i1+3) * RR(i,i1) )
    End Do
   End If

  Else
   Do i1 = 1,3
    If (len.Eq.6) coupL = coupL + Conjg( RSf(j,i1+3) * RR(i,i1) )
    coupR = coupR + Conjg( RSf(j,i1) ) * RL(i,i1)
   End Do
  End If

 End Subroutine CoupGravitinoSfermion3

 Subroutine CoupNeutralinoSdown3(i, j, k, gp, g, RSf, RfL, RfR, Yuk, N &
                                  &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, down-quark, and sdowns
 ! valid for the 3-generation MSSM, epsilon model and spontaneous R-parity
 ! breaking
 !  i .......... index of down quark
 !  j .......... index of neutralino
 !  k .......... index of sdown
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RfL(i,j) ... mixing matrix of the left-handed down quarks
 !  RfR(i,j) ... mixing matrix of the right-handed down quarks
 !  RSf(i,j) ... mixing matrix of the sdowns
 !  Yuk ........ down Yukawa coupling 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f_i) ( coupL P_L + coupR P_R) chi^0_j \tilde(f)_k  
 ! written by Werner Porod, 19.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in), Dimension(3,3) :: RfL, RfR, Yuk
  Complex(dp), Intent(in) :: RSf(6,6), N(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSdown3'

  T3 = -0.5_dp
  e = -1._dp / 3._dp

  Call CoupNeutralinoSfermion3(i, j, k, gp, g, T3, e, RSf, RfL, RfR, Yuk, N &
                              &, coupL, coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSdown3

 Subroutine CoupNeutralinoSfermion3(i, j, k, gp, g, T3, e, RSf, RfL, RfR, Yuk &
                                  &, N, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, fermion, and sfermion
 ! valid for the 3-generation MSSM, epsilon model and model of spontaneous
 ! R-parity breaking
 !  i .......... index of fermion
 !  j .......... index of neutralino
 !  k .......... index of sfermion
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  T3 ......... weak isospin of the fermion
 !  e .......... charge of the fermion
 !  RfL(i,j) ... mixing matrix of the left-handed fermions
 !  RfR(i,j) ... mixing matrix of the right-handed fermions
 !  RSf(i,j) ... mixing matrix of the sfermion
 !  Yuk ........ Yukawa coupling of the fermion 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f_i) ( coupL P_L + coupR P_R) chi^0_j \tilde(f)_k  
 ! written by Werner Porod, 11.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, T3, gp, e
  Complex(dp), Intent(in), Dimension(3,3) :: RfL, RfR, Yuk
  Complex(dp), Intent(in) :: RSf(6,6), N(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2
  Complex(dp) :: fl, fr, hl(3,3), hr(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSfermion3'

  If ( T3.Gt.0.) Then
   hl = - Yuk * Conjg( N(j,4) )
  Else
   hl = - Yuk * Conjg( N(j,3) )
  End If
  hr = Conjg(hl)

  fL = - sqrt2 *( (e-T3) * gp * N(j,1) + T3 * g * N(j,2) )
  fR = sqrt2 * e * gp * Conjg( N(j,1) ) 

  coupL = 0._dp
  coupR = 0._dp

  If ( (CompareMatrices(id3C,RfR,NearlyZero)) .And. &
     & (CompareMatrices(id3C,RfL,NearlyZero))       )Then
   coupL = Conjg( Rsf(k,i+3) ) * fr
   coupR = Conjg( Rsf(k,i)   ) * fl
   Do i2 = 1,3
    coupL = coupL + Conjg( Rsf(k,i2) ) * hl(i2,i) 
    coupR = coupR + Conjg( Rsf(k,i2+3)) * hr(i,i2)
   End Do

  Else If (CompareMatrices(id3C,RfL,NearlyZero)) Then
   coupR = Conjg( Rsf(k,i) ) * fl
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1+3) * RfR(i,i1) ) * fr
    coupR = coupR + Conjg( Rsf(k,i1+3)) * hr(i,i1)
    Do i2 = 1,3
     coupL = coupL + Conjg( Rsf(k,i1) * RfR(i,i2) ) * hl(i1,i2) 
    End Do
   End Do

  Else If (CompareMatrices(id3C,RfR,NearlyZero)) Then
   coupL = Conjg( Rsf(k,i+3) ) * fr
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1) ) * hl(i1,i) 
    coupR = coupR + Conjg( Rsf(k,i1) ) * RfL(i,i1) * fl
    Do i2 = 1,3
     coupR = coupR + Conjg( Rsf(k,i2+3) ) * RfL(i,i1) * hr(i1,i2)
    End Do
   End Do

  Else
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1+3) * RfR(i,i1) ) * fr
    coupR = coupR + Conjg( Rsf(k,i1) ) * RfL(i,i1) * fl
    Do i2 = 1,3
     coupL = coupL + Conjg( Rsf(k,i1) * RfR(i,i2) ) * hl(i1,i2) 
     coupR = coupR + Conjg( Rsf(k,i2+3)) * RfL(i,i1) * hr(i1,i2)
    End Do
   End Do

  End If

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSfermion3

 Subroutine CoupNeutralinoSlepton3(i, j, k, gp, g, RSf, RfL, RfR, Yuk, N &
                                  &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, lepton, and slepton
 ! valid for the 3-generation MSSM
 !  i .......... index of lepton
 !  j .......... index of neutralino
 !  k .......... index of slepton
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RfL(i,j) ... mixing matrix of the left-handed leptons
 !  RfR(i,j) ... mixing matrix of the right-handed leptons
 !  RSf(i,j) ... mixing matrix of the slepton
 !  Yuk ........ Yukawa coupling of the lepton 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f_i) ( coupL P_L + coupR P_R) chi^0_j \tilde(f)_k  
 ! written by Werner Porod, 15.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in), Dimension(3,3) :: RfL, RfR, Yuk
  Complex(dp), Intent(in) :: RSf(6,6), N(4,4)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSlepton3'

  T3 = -0.5_dp
  e = -1._dp 
  Call CoupNeutralinoSfermion3(i, j, k, gp, g, T3, e, RSf, RfL, RfR, Yuk, N &
                              &, coupL, coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSlepton3


 Subroutine CoupNeutralinoSneutrino3(i,j,k,gp,g,N,Rsneut,Rnu,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, neutrino, and sneutrino
 ! valid for the 3-generation MSSM
 !  i .......... index of neutrino
 !  j .......... index of neutralino
 !  k .......... index of sneutrino
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupR ...... the right coupling(i,j,k)
 ! the lagrangian is given by:
 !  \bar{\nu}_i ( coupR P_R ) \chi^0_j \tilde{\nu}_k
 ! written by Werner Porod: 18.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: gp,g
  Complex(dp), Intent(in) :: N(4,4), Rsneut(3,3), Rnu(3,3)
  Complex(dp), Intent(out) :: coupR
  Integer, Intent(in) :: i, j, k

  Integer :: i1
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSneutrino3'

  coupR = 0._dp
  If ( CompareMatrices(Rsneut,id3C,NearlyZero) .And. &
     & CompareMatrices(Rnu,id3C,NearlyZero)           ) Then
   If (i.Eq.k) coupR = 1._dp

  Else If ( CompareMatrices(Rsneut,id3C,NearlyZero) ) Then
   coupR = Rnu(i,k)
 
  Else If ( CompareMatrices(Rnu,id3C,NearlyZero) ) Then
   coupR = Conjg( Rsneut(k,i) )
 
  Else
   Do i1= 1,3
    coupR = coupR + Conjg( Rsneut(k,i1) ) * Rnu(i,i1)
   End Do
  End If

  coupR = sqrt2 * 0.5_dp * ( gp * N(j,1) - g * N(j,2) ) * coupR

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSneutrino3

 Subroutine CoupNeutralinoSup3(i, j, k, gp, g, RSf, RfL, RfR, Yuk, N &
                                  &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, up-quark, and sups
 ! valid for the 3-generation MSSM, epsilon model and model of spontaneous
 ! R-parity breaking
 !  i .......... index of up quark
 !  j .......... index of neutralino
 !  k .......... index of sup
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  RfL(i,j) ... mixing matrix of the left-handed up quarks
 !  RfR(i,j) ... mixing matrix of the right-handed up quarks
 !  RSf(i,j) ... mixing matrix of the sups
 !  Yuk ........ up Yukawa coupling 
 !  N .......... mixing matrix of the neutralino
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f_i) ( coupL P_L + coupR P_R) chi^0_j \tilde(f)_k  
 ! written by Werner Porod, 19.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in), Dimension(3,3) :: RfL, RfR, Yuk
  Complex(dp), Intent(in) :: RSf(6,6), N(:,:)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSup3'

  T3 = 0.5_dp
  e = 2._dp / 3._dp

  Call CoupNeutralinoSfermion3(i, j, k, gp, g, T3, e, RSf, RfL, RfR, Yuk, N &
                              &, coupL, coupR)

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSup3

 Subroutine CoupPseudoScalarSfermion3_3(i,j,k,RP0,T3,yuk,Rsf,A,bi,coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a pseudoscalar and sfermions
 ! valid for the 3-generation MSSM and the 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk(i,j) ... Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A(i,j) ..... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 24.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(6,6), yuk(3,3), A(3,3), bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RP0(:,:), T3
  Integer, Intent(in) :: i, j, k

  Integer :: n_P0, i1, i2, i3
  Complex(dp) :: ayuk(3,3), aA(3,3), abi(Size(bi)), parts(Size(RP0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion3_3'

  coup = ZeroC

  n_P0 =  Size(RP0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  If (T3.Gt.0._dp) Then
   Do i1 = 1,3
    Do i2 = 1,3
     parts(1) = parts(1)                                                  &
            & - ayuk(i1,i2) * bi(1) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) )    &
            & + yuk(i1,i2) * abi(1) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) ) 

     parts(2) = parts(2) + A(i1,i2) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )    &
            &            - aA(i1,i2) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) )

     Do i3=3,n_P0
      parts(i3) = parts(i3) &
             &  - yuk(i1,i2) * abi(i3-1) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) ) &
             &  + ayuk(i1,i2) * bi(i3-1) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) )
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,3
    Do i2 = 1,3
     parts(1) = parts(1) + A(i1,i2) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )    &
            &            - aA(i1,i2) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) )

     parts(2) = parts(2)                                                  &
            & - ayuk(i1,i2) * bi(1) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) )    &
            & + yuk(i1,i2) * abi(1) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) ) 

    End Do
   End Do

  End If

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * parts(i1)
  End Do

  coup = (0._dp, -1._dp) * coup * oosqrt2

  Iname = Iname - 1

 End Subroutine CoupPseudoScalarSfermion3_3


 Subroutine CoupPseudoScalarSfermion3a_3(i,j,k,RP0,T3,yuk,Rsf,A,bi,coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a pseudoscalar and sfermions
 ! valid for the 3-generation MSSM and the 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk(i,j) ... Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A(i,j) ..... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 24.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(6,6), yuk(3,3), A(3,3), bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RP0(:,:), T3
  Integer, Intent(in) :: i, j, k

  Integer :: n_P0, i1, i2, i3
  Complex(dp) :: ayuk(3,3), aA(3,3), abi(Size(bi)), parts(Size(RP0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion3a_3'

  coup = ZeroC

  n_P0 =  Size(RP0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  If (T3.Gt.0._dp) Then
   Do i1 = 1,3
    Do i2 = 1,3
     parts(1) = parts(1)                                                  &
            & - ayuk(i1,i2) * bi(1) * Rsf(k,i1) * Conjg( id6c(j,i2+3) )    &
            & + yuk(i1,i2) * abi(1) * Rsf(k,i2+3) * Conjg( id6c(j,i1) ) 

     parts(2) = parts(2) + A(i1,i2) * Rsf(k,i2+3) * Conjg( id6c(j,i1) )    &
            &            - aA(i1,i2) * Rsf(k,i1) * Conjg( id6c(j,i2+3) )

     Do i3=3,n_P0
      parts(i3) = parts(i3) &
             &  - yuk(i1,i2) * abi(i3-1) * Rsf(k,i2+3) * Conjg( id6c(j,i1) ) &
             &  + ayuk(i1,i2) * bi(i3-1) * Rsf(k,i1) * Conjg( id6c(j,i2+3) )
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,3
    Do i2 = 1,3
     parts(1) = parts(1) + A(i1,i2) * Rsf(k,i2+3) * Conjg( id6c(j,i1) )    &
            &            - aA(i1,i2) * Rsf(k,i1) * Conjg( id6c(j,i2+3) )

     parts(2) = parts(2)                                                  &
            & - ayuk(i1,i2) * bi(1) * Rsf(k,i1) * Conjg( id6c(j,i2+3) )    &
            & + yuk(i1,i2) * abi(1) * Rsf(k,i2+3) * Conjg( id6c(j,i1) ) 

    End Do
   End Do

  End If

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * parts(i1)
  End Do

  coup = (0._dp, -1._dp) * coup * oosqrt2

  Iname = Iname - 1

 End Subroutine CoupPseudoScalarSfermion3a_3


 Subroutine CoupPseudoscalarSfermion4_3(i, j, k, l, RP0, T3, e, yuk,  &
   &                                    Rsf, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between pseudoscalar and sfermions
 ! valid for the 3-generation MSSM,  3-generation epsilon model and the 
 ! full RP model
 !  i,j ........ index of pseudoscalar boson
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  n_P0 ....... specifies the model: 2 -> MSSM
 !                                    3 -> 1-generation epsilon model
 !                                    5 -> 3-generation epsilon model
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 29.11.1999
 ! 10.10.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(:,:), yuk(3,3)
  Real(dp), Intent(in) :: RP0(:,:), g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: n_P0, i1, n_Sf, i2 ,i3
  Real(dp) :: g2, gp2, YL, YR, yuk2
  Complex(dp) :: parts(12), Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupPseudoScalarSfermion4_3'

  coup = Cmplx( 0._dp,0._dp,dp )

  n_P0 = Size(RP0, Dim=1)
  n_Sf = Size(RSf, Dim=1)

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = ZeroC
  yuk2 = 0
  Dterm = 0._dp

  If (n_Sf.Eq.6) Then ! left and right sfermions
   Do i1=1,3
    Do i2 = 1,3
     Do i3 = 1,3
      yuk2 = yuk2 &
         & + Conjg( Rsf(k,i1) ) * yuk(i1,i2) * Conjg(yuk(i3,i2)) * Rsf(l,i3) &
         & + Conjg( Rsf(k,i1+3) ) * Conjg(yuk(i2,i1)) * yuk(i2,i3)           &
         &                        * Rsf(l,i3+3)
     End Do
    End Do
    Dterm = Dterm                  &
      &   + 0.25_dp * ( (T3 * g2 - YL * gp2) * Rsf(l,i1) * Conjg( Rsf(k,i1) ) &
      &              - YR * gp2 * Rsf(l,i1+3) * Conjg( Rsf(k,i1+3) ) )
   End Do

  Else ! Sneutrinos
   if (l.eq.k) then
    Dterm = 0.125_dp * (g2 + gp2)
   else 
    Dterm = 0._dp
   end if

  End If


  If (T3.Gt.0._dp) Then
   parts(1) = -  Dterm
   parts(2) = Dterm - 0.5_dp * yuk2

  Else
   parts(1) = - Dterm - 0.5_dp * yuk2
   parts(2) = Dterm

  Endif

  Do i1=3,n_P0
   parts(i1) = - Dterm
  Enddo

  Do i1=1,n_P0
   coup = coup + RP0(i,i1) * RP0(j,i1) * parts(i1)
  Enddo

  Iname = Iname - 1

 End Subroutine CoupPseudoscalarSfermion4_3


 Subroutine CoupScalarSfermion3MSSM_3(i, j, k, RS0, T3, e, yuk, Rsf, A, mu, &
                      &          vevSM, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the 3-generation MSSM
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk(i,j) ... Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A(i,j) ..... A-parameter
 !  mu ......... bilinear parameters of the superpotential
 !  vevSM(i) ... vevs of the Higgs bosons
 !  gU1 ........ U(1) coupling
 !  gSU2 ....... SU(2) coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 25.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(:,:) , yuk(3,3) , A(3,3), mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(2,2), T3, e, vevSM(2), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2, i3, n_sfer
  Complex(dp) :: ayuk(3,3), aA(3,3), muC, parts(2)
  Real(dp) :: g2, gp2, YL, YR, Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3MSSM_3'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  muC = Conjg( mu )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  n_sfer = Size( Rsf, Dim=1 )
  If (T3.Gt.0._dp) Then
   If (n_sfer.Eq.3) Then ! Sneutrino
    if (j.eq.k) then
     Dterm = 0.5_dp * (g2*T3 - YL * gp2) 
     parts(1) = parts(1) - vevSM(1) * Dterm
     parts(2) = parts(2) + vevSM(2) * Dterm
    else
     parts(1:2) = 0._dp
    end if

   Else 
    Do i1=1,3
     Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,i1) * Conjg( Rsf(j,i1) ) &
         &         - YR * gp2 * Rsf(k,i1+3) * Conjg( Rsf(j,i1+3) ) )
     parts(1) = parts(1) - vevSM(1) * Dterm
     parts(2) = parts(2) + vevSM(2) * Dterm
 
     Do i2=1,3 
      parts(1) = parts(1)                                                 &
        & + ( yuk(i1,i2) * muC * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )         &
        &   + ayuk(i1,i2) * mu * Rsf(k,i1) * Conjg( Rsf(j,i2+3) ) ) * oosqrt2
 
      parts(2) = parts(2)                                                 &
        & - ( A(i1,i2) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )                 &
        &   + aA(i1,i2) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) ) ) * oosqrt2
      Do i3=1,3
      parts(2) = parts(2)                                                 &
             & - vevSM(2) * ( yuk(i1,i3) * ayuk(i2,i3)                    &
             &                         * Rsf(k,i2) * Conjg( Rsf(j,i1) )   &
             &              + yuk(i3,i1) * ayuk(i3,i2)                    &
             &                       * Rsf(k,i1+3) * Conjg( Rsf(j,i2+3) ) )
      End Do
     End Do
    End Do
   End If
  Else
   Do i1=1,3
    Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,i1) * Conjg( Rsf(j,i1) ) &
        &         - YR * gp2 * Rsf(k,i1+3) * Conjg( Rsf(j,i1+3) ) )
    parts(1) = parts(1) - vevSM(1) * Dterm
    parts(2) = parts(2) + vevSM(2) * Dterm

    Do i2=1,3 
     parts(1) = parts(1)                                                 &
       & - ( A(i1,i2) * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )                 &
       &   + aA(i1,i2) * Rsf(k,i1) * Conjg( Rsf(j,i2+3) ) ) * oosqrt2

     parts(2) = parts(2)                                                 &
       & + ( yuk(i1,i2) * muC * Rsf(k,i2+3) * Conjg( Rsf(j,i1) )         &
       &   + ayuk(i1,i2) * mu * Rsf(k,i1) * Conjg( Rsf(j,i2+3) ) ) * oosqrt2
     Do i3=1,3
     parts(1) = parts(1)                                                 &
            & - vevSM(1) * ( yuk(i1,i3) * ayuk(i2,i3)                    &
            &                         * Rsf(k,i2) * Conjg( Rsf(j,i1) )   &
            &              + yuk(i3,i1) * ayuk(i3,i2)                    &
            &                       * Rsf(k,i1+3) * Conjg( Rsf(j,i2+3) ) )
     End Do
    End Do
   End Do

  End If

  Do i1=1,2
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3MSSM_3


 Subroutine CoupScalarSfermion3MSSMa_3(i, j, k, RS0, T3, e, yuk, Rsf, A, mu, &
                      &          vevSM, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the 3-generation MSSM
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk(i,j) ... Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A(i,j) ..... A-parameter
 !  mu ......... bilinear parameters of the superpotential
 !  vevSM(i) ... vevs of the Higgs bosons
 !  gU1 ........ U(1) coupling
 !  gSU2 ....... SU(2) coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 25.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(:,:) , yuk(3,3) , A(3,3), mu
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(2,2), T3, e, vevSM(2), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2, i3, n_sfer
  Complex(dp) :: ayuk(3,3), aA(3,3), muC, parts(2)
  Real(dp) :: g2, gp2, YL, YR, Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3MSSMa_3'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  muC = Conjg( mu )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  n_sfer = Size( Rsf, Dim=1 )
  If (T3.Gt.0._dp) Then
   If (n_sfer.Eq.3) Then ! Sneutrino
    Do i1=1,3
     Dterm = 0.5_dp * (g2*T3 - YL * gp2) * Rsf(k,i1) * Conjg( id3c(j,i1) ) 
     parts(1) = parts(1) - vevSM(1) * Dterm
     parts(2) = parts(2) + vevSM(2) * Dterm
   End Do

   Else 
    Do i1=1,3
     Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,i1) * Conjg( id6c(j,i1) ) &
         &         - YR * gp2 * Rsf(k,i1+3) * Conjg( id6c(j,i1+3) ) )
     parts(1) = parts(1) - vevSM(1) * Dterm
     parts(2) = parts(2) + vevSM(2) * Dterm
 
     Do i2=1,3 
      parts(1) = parts(1)                                                 &
        & + ( yuk(i1,i2) * muC * Rsf(k,i2+3) * Conjg( id6c(j,i1) )         &
        &   + ayuk(i1,i2) * mu * Rsf(k,i1) * Conjg( id6c(j,i2+3) ) ) * oosqrt2
 
      parts(2) = parts(2)                                                 &
        & - ( A(i1,i2) * Rsf(k,i2+3) * Conjg( id6c(j,i1) )                 &
        &   + aA(i1,i2) * Rsf(k,i1) * Conjg( id6c(j,i2+3) ) ) * oosqrt2
      Do i3=1,3
      parts(2) = parts(2)                                                 &
             & - vevSM(2) * ( yuk(i1,i3) * ayuk(i2,i3)                    &
             &                         * Rsf(k,i2) * Conjg( id6c(j,i1) )   &
             &              + yuk(i3,i1) * ayuk(i3,i2)                    &
             &                       * Rsf(k,i1+3) * Conjg( id6c(j,i2+3) ) )
      End Do
     End Do
    End Do
   End If
  Else
   Do i1=1,3
    Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,i1) * Conjg( id6c(j,i1) ) &
        &         - YR * gp2 * Rsf(k,i1+3) * Conjg( id6c(j,i1+3) ) )
    parts(1) = parts(1) - vevSM(1) * Dterm
    parts(2) = parts(2) + vevSM(2) * Dterm

    Do i2=1,3 
     parts(1) = parts(1)                                                 &
       & - ( A(i1,i2) * Rsf(k,i2+3) * Conjg( id6c(j,i1) )                 &
       &   + aA(i1,i2) * Rsf(k,i1) * Conjg( id6c(j,i2+3) ) ) * oosqrt2

     parts(2) = parts(2)                                                 &
       & + ( yuk(i1,i2) * muC * Rsf(k,i2+3) * Conjg( id6c(j,i1) )         &
       &   + ayuk(i1,i2) * mu * Rsf(k,i1) * Conjg( id6c(j,i2+3) ) ) * oosqrt2
     Do i3=1,3
     parts(1) = parts(1)                                                 &
            & - vevSM(1) * ( yuk(i1,i3) * ayuk(i2,i3)                    &
            &                         * Rsf(k,i2) * Conjg( id6c(j,i1) )   &
            &              + yuk(i3,i1) * ayuk(i3,i2)                    &
            &                       * Rsf(k,i1+3) * Conjg( id6c(j,i2+3) ) )
     End Do
    End Do
   End Do

  End If

  Do i1=1,2
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3MSSMa_3

 Subroutine CoupScalarSfermion4MSSM_3(i, j, k, l, RS0, T3, e, yuk, Rsf &
                                    &, gp, g, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-coupling between neutral scalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,j ........ index of scalar boson
 !  k .......... index of the sfermion \tilde{f}
 !  l .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  e .......... electric charge of sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 09.11.01
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(:,:), yuk(3,3)
  Real(dp), Intent(in) :: RS0(2,2), g, T3, gp, e
  Integer, Intent(in) :: i, j, k, l
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2, i3, i_sf
  Complex(dp) :: parts(2), Dterm
  Real(dp) :: g2, gp2, YL, YR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion4MSSM_3'

  coup = Cmplx( 0._dp,0._dp,dp )

  g2 = g**2
  gp2 = gp**2
  YL = e - T3
  YR = - e

  parts = Cmplx( 0._dp,0._dp,dp )

  i_sf = Size(Rsf, dim=1)

  Dterm = 0._dp
  If (i_sf.Eq.6) Then
   Do i1=1,3
    Dterm = Dterm                                              &
        & + 0.25_dp * ( (T3 * g2 - YL * gp2)                    &
        &             * Rsf(l,i1) * Conjg( Rsf(k,i1) )  &
        &             - YR * gp2 * Rsf(l,i1+3) * Conjg( Rsf(k,i1+3) ) )
   End Do
  Else ! sneutrinos
   Do i1=1,3
    Dterm = Dterm                                        &
        & + 0.25_dp * ( T3 * g2 - YL * gp2)        &
        &          * Rsf(l,i1) * Conjg( Rsf(k,i1) )
   End Do
  End If 

  If (T3.Gt.0._dp) Then
   parts(1) = - Dterm

   parts(2) = Dterm
   If (i_sf.Eq.6) Then
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       parts(2) = parts(2)  &
              & - 0.5_dp * ( Conjg(Rsf(k,i1)) * yuk(i1,i2)              &
              &              * Conjg(yuk(i3,i2)) * Rsf(l,i3)            &
              &            + Conjg(Rsf(k,i1+3)) * Conjg(yuk(i2,i1))     &
              &              * yuk(i2,i3) * Rsf(l,i3+3)   )
      End Do
     End Do
    End Do
   End If 

  Else

   parts(1) = - Dterm
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      parts(1) = parts(1)  &
             & - 0.5_dp * ( Conjg(Rsf(k,i1)) * yuk(i1,i2)              &
             &              * Conjg(yuk(i3,i2)) * Rsf(l,i3)            &
             &            + Conjg(Rsf(k,i1+3)) * Conjg(yuk(i2,i1))     &
             &              * yuk(i2,i3) * Rsf(l,i3+3)   )
     End Do
    End Do
   End Do

   parts(2) = Dterm

  Endif

  Do i1=1,2
   coup = coup + RS0(i,i1) * RS0(j,i1) * parts(i1)
  Enddo
  Iname = Iname - 1

 End Subroutine CoupScalarSfermion4MSSM_3

 Subroutine CoupSfermion4Y_3(i, j, k, l, T3_1, Yuk1, Rsf1, T3_2, yuk2, Rsf2, nc &
                            &, coup)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion self coupling for different flavours of
 ! sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i ........ index of the adjungated sfermion 1
 !  j ........ index of the sfermion 1
 !  k ........ index of the adjungated sfermion 2
 !  l ........ index of the sfermion 2
 !  Yuk1 ........ Yukawa coupling of sfermion1
 !  Rsf1(i,j) ... mixing matrix of the sfermions1
 !  Yuk2 ........ Yukawa coupling of sfermion2
 !  Rsf2(i,j) ... mixing matrix of the sfermions2
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 11.11.1999
 ! last change: 11.11.1999
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, nc
  Real(dp), Intent(in) :: T3_1, T3_2
  Complex(dp), Intent(in) :: Yuk1(3,3), Rsf1(:,:), Yuk2(3,3), Rsf2(:,:)
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2 , i3, i4, n_sf1, n_sf2

  coup = ZeroC
  n_sf1 = Size(Rsf1, dim=1)
  n_sf2 = Size(Rsf2, dim=1)

  If ((T3_1*T3_2).Gt.0._dp) Then
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       coup = coup &
          & - Yuk1(i1,i2) * Conjg(Yuk2(i3,i4) * Rsf1(i,i1) * Rsf2(k,i4+3)) &
          &               * Rsf1(j,i2+3) * Rsf2(l,i3)                      &
          & - Yuk2(i3,i4) * Conjg(Yuk1(i1,i2) * Rsf1(i,i2+3) * Rsf2(k,i3)) &
          &               * Rsf1(j,i1) * Rsf2(l,i4+3)
      End Do
     End Do
    End Do
   End Do

  Else If (n_sf1.Eq.3) Then ! sneutrinos
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       coup = coup &
          & - Yuk2(i3,i4) * Conjg(Yuk2(i1,i2) * Rsf1(i,i1) * Rsf2(k,i4+3)) &
          &               * Rsf1(j,i3) * Rsf2(l,i2+3)
      End Do
     End Do
    End Do
   End Do

  Else If (n_sf2.Eq.3) Then ! sneutrinos
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       coup = coup &
          & - Yuk1(i1,i2) * Conjg(Yuk1(i3,i4) * Rsf1(i,i4+3) * Rsf2(k,i1)) &
          &               * Rsf1(j,i2+3) * Rsf2(l,i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       coup = coup &
          & - Yuk1(i1,i2) * Conjg(Yuk1(i3,i4) * Rsf1(i,i4+3) * Rsf2(k,i1)) &
          &               * Rsf1(j,i2+3) * Rsf2(l,i3)                      &
          & - Yuk2(i3,i4) * Conjg(Yuk2(i1,i2) * Rsf1(i,i1) * Rsf2(k,i4+3)) &
          &               * Rsf1(j,i3) * Rsf2(l,i2+3)
      End Do
     End Do
    End Do
   End Do
  End If


 End Subroutine CoupSfermion4Y_3


 Subroutine CoupSfermionSelf4Y_3(i, j, k, l, Yuk, Rsf, nc, coup, self)
 !-----------------------------------------------------------------------
 ! calculates the 4-sfermion self coupling 
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i,k ........ index of sfermions
 !  j,l ........ index of the adjungated sfermions
 !  Yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  nc ......... number of colours 
 !  self ....... integer, if present it is assumed that the first two
 !               indices belong to an electroweak eigenstate sfermion
 ! output
 !  coup ....... the coupling(i,j,k,l)
 ! written by Werner Porod, 31.12.01
 ! 17.03.03: bug fixing due to SU(3) structure, assumed that colour
 !           index for i is the samke as j, and k has the same colour
 !           index as l 
 !     if (self==1) then we assume that the sfermions with indices i and j
 !     are electroweak eigenstates
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k, l, nc
  Integer, Intent(in) :: self
  Complex(dp), Intent(in) :: Yuk(3,3), Rsf(6,6)
  Complex(dp), Intent(out) :: coup

  Integer :: i1, i2, i3, i4, ii(9), jj(9)
  Complex(dp) :: parts(3,3), coupP, coupN
  Real(dp) :: maxi, matr33(3,3)

  coup = ZeroC
  If (self.Eq.1) Then
   If ((i.Le.3).And.(j.Ge.4)) Then
     Do i1=1,3
      Do i2=1,3
       parts(i1,i2) =  - Conjg( Yuk(i1,i2) * Rsf(k,i2+3)) * Rsf(l,i1) * nc
      End Do
     End Do
     parts = parts * Yuk(i,j-3)

   Else If  ((j.Le.3).And.(i.Ge.4)) Then
     Do i1=1,3
      Do i2=1,3
       parts(i1,i2) = - Yuk(i1,i2) * Conjg(Rsf(k,i1) ) * Rsf(l,i2+3) * nc
      End Do
     End Do
     parts = parts * Conjg( Yuk(j,i-3) )

   Else If  ((j.Le.3).And.(i.Le.3)) Then
     Do i1=1,3
      Do i2=1,3
       parts(i1,i2) = - Yuk(i,i2)* Conjg( Yuk(j,i1) ) &
            &        * Conjg(Rsf(k,i1+3) ) * Rsf(l,i2+3)
      End Do
     End Do

   Else If  ((j.Ge.4).And.(i.Ge.4)) Then
     Do i1=1,3
      Do i2=1,3
       parts(i1,i2) = - Yuk(i2,j-3)* Conjg( Yuk(i1,i-3) ) &
            &        * Conjg(Rsf(k,i2) ) * Rsf(l,i1)
      End Do
     End Do

   End If
   !-----------------------------------------
   ! to stabilize the numerical calculation
   !-----------------------------------------
   matr33 = Abs(Real(parts,dp))
   Do i1=1,9
    maxi = Maxval(matr33)
    Call FindPosition(3, matr33, Maxi, ii(i1), jj(i1))
    matr33(ii(i1), jj(i1)) = 0._dp
   End Do
   coupN = ZeroC
   coupP = ZeroC
   matr33 = Real(parts,dp)
   Do i1=9,1,-1
    if (matr33(ii(i1), jj(i1)).lt.0._dp) then
     coupN = coupN + parts(ii(i1), jj(i1)) 
    else
     coupP = coupP + parts(ii(i1), jj(i1)) 
    end if
   End Do
   coup = coupN + coupP

  Else

   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       coup = coup &
          & - Yuk(i1,i2) * Conjg( Yuk(i3,i4) )                              &
          &   * ( Conjg(Rsf(i,i1) * Rsf(k,i4+3)) * Rsf(j,i2+3) * Rsf(l,i3)  &
          &     + Conjg(Rsf(k,i1) * Rsf(i,i4+3)) * Rsf(l,i2+3) * Rsf(j,i3) )
      End Do
     End Do
    End Do
   End Do
  End If

 End Subroutine CoupSfermionSelf4Y_3

 Subroutine AllCouplingsEps1(g, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R             &
    & , Y_u, uU_L, uU_R, vevSM, vevL                                         &
    & , RSpm, RP0, RS0, Umat, Vmat, Nmat, bi, phase_glu                      &
    & , RSlepton, Ae, RSup, Au, RSdown, Ad, GenerationMixing                 &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03       &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R       &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R    &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R   &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L           &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R           &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L         &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R           &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R           &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R           &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0         &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03 &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW    &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW        &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L    &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R             &
    & , cpl_SmpCN_L, cpl_SmpCN_R)
 !--------------------------------------------------------------------------
 ! Routine for calculating all couplings of the one-generation espilon model
 ! output:
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 !  in the case that fermions are involved, _L and _R refer to left- and
 !  right couplings
 ! written by Werner Porod, 
 ! 05.09.2001: taking AllCouplingsMSSM as basis
 ! 12.09.2002: changing interface to avoid global variables
 !--------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: uL_L(3,3) ! mixing matrix of left leptons
  Complex(dp), Intent(in) :: uL_R(3,3) ! mixing matrix of right leptons
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3) ! mixing matrix of left d-quarks
  Complex(dp), Intent(in) :: uD_R(3,3) ! mixing matrix of right d-quarks
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3) ! mixing matrix of left u-quarks
  Complex(dp), Intent(in) :: uU_R(3,3) ! mixing matrix of right u-quarks
  Real(dp), Intent(in) :: vevSM(2)     ! MSSM Higgs vevs [v_d, v_u]
  Real(dp), Intent(in) :: vevL         ! tau-sneutrino vev
  Complex(dp), Intent(in) :: RSpm(4,4) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(3,3)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(3,3)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(3,3), Vmat(3,3) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(5,5) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  bi(2)   ! superpotential bilinears [mu, epsilon_3]
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: RSlepton(6,6) ! slepton mixing matrix
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters
  Logical, Intent(in) :: GenerationMixing ! if .true. generation mixing of
                                          ! (s)fermions is taken into account

  Complex(dp), Intent(out) :: cpl_SmpSlSn(4,6,3), cpl_SmpSdSu(4,6,6)   &
      & , cpl_SmpSnSl(4,3,6), cpl_SmpSuSd(4,6,6), cpl_SmpP03(4,4,3)    &
      & , cpl_SmpP0W(4,3), cpl_SmpS03(4,4,3), cpl_SmpS0W(4,3)          &
      & , cpl_SmpLNu_L(4,3,3), cpl_SmpLNu_R(4,3,3), cpl_SmpDU_L(4,3,3) &
      & , cpl_SmpDU_R(4,3,3), cpl_SmpZ(4,4), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(3,3), cpl_CCZ_R(3,3)           &
      & , cpl_NNZ_L(5,5), cpl_NNZ_R(5,5), cpl_NNS0_L(5,5,3)            &
      & , cpl_NNS0_R(5,5,3), cpl_NNP0_L(5,5,3), cpl_NNP0_R(5,5,3) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,5,6), cpl_DNSd_R(3,5,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,5,6), cpl_UNSu_R(3,5,6)         &
      & , cpl_LNSl_L(3,5,6), cpl_LNSl_R(3,5,6), cpl_NuNSn_L(3,5,3)      & 
      & , cpl_NuNSn_R(3,5,3), cpl_DDP0_L(3,3,3), cpl_LLP0_L(3,3,3)      &
      & , cpl_UUP0_L(3,3,3), cpl_DDP0_R(3,3,3), cpl_LLP0_R(3,3,3)       &
      & , cpl_UUP0_R(3,3,3), cpl_DDS0_L(3,3,3), cpl_LLS0_L(3,3,3)       &
      & , cpl_UUS0_L(3,3,3), cpl_DDS0_R(3,3,3), cpl_LLS0_R(3,3,3)       &
      & , cpl_UUS0_R(3,3,3)
  Complex(dp), Intent(out) :: cpl_CUSd_L(3,3,6), cpl_CUSd_R(3,3,6)      &
      & , cpl_CDSu_L(3,3,6), cpl_CDSu_R(3,3,6), cpl_CLSn_L(3,3,3)       &
      & , cpl_CLSn_R(3,3,3), cpl_CNuSl_L(3,3,6), cpl_CNuSl_R(3,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(3)
  Complex(dp) :: cpl_P0SdSd(3,6,6), cpl_P0SuSu(3,6,6), cpl_P0SlSl(3,6,6) &
      & , cpl_P0SnSn(3,3,3), cpl_P0S0Z(3,3) 
  Real(dp), Intent(out) :: cpl_P0S03(3,3,3)
  Complex(dp), Intent(out) :: cpl_S0SdSd(3,6,6), cpl_S0SuSu(3,6,6) &
      & , cpl_S0SlSl(3,6,6), cpl_S0SnSn(3,3,3) 
  Real(dp), Intent(out) :: cpl_S03(3,3,3), cpl_S0WW(3), cpl_S0ZZ(3), cpl_FFpW &
      & , cpl_LNuW(3,3)
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(3,3,3), cpl_CCP0_R(3,3,3)    &
      & , cpl_CCS0_L(3,3,3), cpl_CCS0_R(3,3,3), cpl_CNW_L(3,5)        &
      & , cpl_CNW_R(3,5), cpl_SmpCN_L(4,3,5), cpl_SmpCN_R(4,3,5)

  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3, i4
  Real(dp) :: gU1, gSU2, gSU3, e_d, e_u, sinW2, cosW2, cosW, vL(1)
  Complex(dp) :: Rsd(2,2), Rsu(2,2), Rsl(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, mu, YukL(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsEps1'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 3
  n_neut = 5
  n_P0 = 3
  n_S0 = 3
  n_Spm = 4
  vL(1) = vevL
  !--------------------
  ! some constants
  !--------------------
  e_u = 2._dp / 3._dp
  e_d = - 1._dp / 3._dp

  Do i1=1,3
   yukL(i1)=y_L(i1,i1)
  End Do
  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm, yukL  &
               &, gU1, gSU2, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpLNu_L = 0._dp
  cpl_SmpLNu_R = 0._dp
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_D, uD_L, uD_R, Y_U &
              &, uU_L, uU_R, cpl_SmpDU_L(i1,i2,i3), cpl_SmpDU_R(i1,i2,i3) )
     End Do
    End Do
   End Do
  Else
   Do i1=1,n_Spm
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  End If
  !------------------------------
  ! charged scalar - pseudoscalar - W
  !------------------------------
  cpl_SmpP0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  cpl_SmpP03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_P0
!     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
!                                       &, cpl_SmpP03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  cpl_SmpS0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------
  ! charged scalar - charged scalar - scalar
  !------------------------------------------
  cpl_SmpS03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_S0
!     Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
!                                  &, cpl_SmpS03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSlSn = ZeroC
  cpl_SmpSnSl = ZeroC
  cpl_SmpSdSu = ZeroC
  cpl_SmpSuSd = ZeroC

!  If (GenerationMixing) Then
!   Do i1 = 1, n_Spm
!    Do i2 = 1,6
!     Do i3 = 1,6
!      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
!                                  & , Y_d, Ad, Rsdown, Y_u, Au, Rsup, coupC)
!      cpl_SmpSdSu(i1, i2, i3) = coupC
!      cpl_SmpSuSd(i1, i3, i2) = Conjg(coupC)
!     End Do
!    End Do
!   End Do

!  Else
!   Do i1 =1,n_Spm
!    Do i2=1,3
!     Yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Yukp = Y_u(i2,i2)
!     Ap = Au(i2,i2)
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
!                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
!       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
!       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
!      End Do
!     End Do
!    End Do ! i2
!   End Do ! i1
!  End If

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  cpl_SmpZ = 0._dp
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, YukL(3), gSU2  &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, YukL(3), gSU2  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp
  cpl_CLSn_L = 0._dp
  cpl_CLSn_R = 0._dp
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_char
    Do i2=1,3
     Do i3=1,6
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                              &, Y_u, uU_L, uU_R, Umat, Vmat, coupLC, coupRC)
      cpl_CUSd_L(i1, i2, i3) = coupLC
      cpl_CUSd_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                              &, Y_u, uD_L, uD_R, Umat, Vmat, coupLC, coupRC)
      cpl_CDSu_L(i1, i2, i3) = coupLC
      cpl_CDSu_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,n_char
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do
   Do i1=1,2
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    YukP = 0._dp
    Do i2=1,n_char
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSl, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CNuSl_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2C, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
     cpl_CLSn_L(i2, i1, i1) = coupLC
     cpl_CLSn_R(i2, i1, i1) = coupRC
    End Do
   End Do

  End If 

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_LLP0_L = 0._dp
  cpl_LLP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_P0
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RP0 &
                                 &, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_d(i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_LLS0_L = 0._dp
  cpl_LLS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_S0
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RS0 &
                           &, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RS0 &
                           &, cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, -0.5_dp, Y_d(i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  cpl_LLZ_L = 0._dp
  cpl_LLZ_R = 0._dp
  cpl_NuNuZ_L = 0._dp
  cpl_NuNuZ_R = 0._dp
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  cpl_FFpW = - gSU2 * oosqrt2
  cpl_LNuW = 0._dp

  If (GenerationMixing) Then
   cpl_DUW = CKM * cpl_FFpW
  Else
   cpl_DUW = id3C * cpl_FFpW
  End If
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_LNSl_L = 0._dp
  cpl_LNSl_R = 0._dp
  cpl_NuNSn_L = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,n_neut
     Do i3 = 1,6
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, Nmat, coupLC, coupRC)
      cpl_DNSd_L(i1, i2, i3 ) = coupLC
      cpl_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, Nmat, coupLC, coupRC)
      cpl_UNSu_L(i1, i2, i3 ) = coupLC
      cpl_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_neut
    Call CoupNeutralinoSneutrino(i1, gU1, gSU2, Nmat, coupRC)
    Do i2 = 1,2
     cpl_NuNSn_R(i2, i1, i2 ) = coupRC
    End Do
   End Do

   Do i1 = 1,2
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    Do i2 = 1,2
     Do i3 = 1,n_neut
      Call CoupNeutralinoSlepton(i3, i2, gU1, gSU2, RSl, Yuk, Nmat &
                               &, coupLC, coupRC)
      cpl_LNSl_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_LNSl_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsdown, uD_L, uD_R, &
                          & coupLC, coupRC)
     cpl_GDSd_L(i1, i2) = coupLC
     cpl_GDSd_R(i1, i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsup, uU_L, uU_R, &
                          & coupLC, coupRC)
     cpl_GUSu_L(i1, i2) = coupLC
     cpl_GUSu_R(i1, i2) = coupRC
    End Do
   End Do
  Else

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  End If 

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, Nmat, RP0, gU1, gSU2, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalar(i1, i2, i3, Nmat, RS0, gU1, gSU2, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  cpl_SlSnW = 0._dp
  cpl_SnSlW = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,6
    Do i2 = 1,6
     Call CoupSfermionW3(i1, i2, gSU2, RSdown, RSup, cpl_SdSuW(i1,i2) )
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do
   End Do
  End If

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
!  Do i1=1,n_P0
!   Do i2=1,n_P0
!    Do i3=1,n_S0
!     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
!                                &, cpl_P0S03(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
    End Do
   End Do 
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp
  cpl_P0SlSl = 0._dp
  cpl_P0SnSn = 0._dp

!  If (GenerationMixing) Then
!   Do i1=1,n_P0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_d, Rsdown   &
!                                   &, Ad, bi, cpl_P0SdSd(i1,i2,i3) )
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, 0.5_dp, Y_u, Rsup      &
!                                   &, Au, bi, cpl_P0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
!
!  Else
!   Do i1=1,n_P0
!    Do i2=1,3
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
!                                   &, A, bi, coupC )
!       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!     yuk = Y_u(i2,i2)
!     A = Au(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
!                                   &, A, bi, coupC )
!       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!    End Do
!   End Do
!
!  End If

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
!  Do i1=1,n_S0
!   Do i2=1,n_S0
!    Do i3=1,n_S0
!     Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
!    End Do
!   End Do
!  End Do

  !-------------------
  ! scalar - W+ - W-
  !-------------------
  cpl_S0WW = 0._dp
  Do i1=1,n_S0
   Call CoupScalarW(i1, gSU2, vevSM, vL, RS0, cpl_S0WW(i1) )
  End Do

  !----------------
  ! scalar - Z - Z
  !----------------
  cpl_S0ZZ = 0._dp
  Do i1=1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevSM, vL, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp
  cpl_S0SlSl = 0._dp
  cpl_S0SnSn = 0._dp

  If (GenerationMixing) Then
!   Do i1=1,n_S0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d, Rsdown   &
!                           &, Ad, mu, vevs, gU1, gSU2, cpl_S0SdSd(i1,i2,i3) )
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u, Rsup      &
!                           &, Au, mu, vevs, gU1, gSU2, cpl_S0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!
  Else
   Do i1=1,n_S0
    Do i2=1,3
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
                               &, A, bi, vevSM, vL, gU1, gSU2, coupC )
       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
                               &, A, bi, vevSM, vL, gU1, gSU2, coupC )
       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  End If

  !----------
  ! scalar W
  !----------
  cpl_S0WW = 0._dp
!  Do i1 = 1,n_S0
!   Call CoupScalarW(i1, gSU2, vevs, RS0, cpl_S0WW(i1) )
!  End Do

  !----------
  ! scalar Z
  !----------
  cpl_S0ZZ = 0._dp
!  Do i1 = 1,n_S0
!   Call CoupScalarZ(i1, gSU2, cosW2, vevs, RS0, cpl_S0ZZ(i1) )
!  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SlSlZ = 0._dp
  cpl_SnSnZ = 0._dp
  cpl_SuSuZ = 0._dp

  If (GenerationMixing) Then
   Do i1=1,6
    Do i2=1,6
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, cpl_SdSdZ(i1, i2) )
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, cpl_SuSuZ(i1, i2) )
    End Do
   End Do
 
  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  End If

  Iname = Iname - 1

 End Subroutine AllCouplingsEps1

 Subroutine CoupChargedScalarSfermion3Eps1(i, j, k, RSpm, g, vevs, bi    &
                       &      , yukL, yukd, Ad, Rsfd, yuku, Au, Rsfu, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between charged scalar and sfermions
 ! valid for the 1-generation epsilon model
 !  i .......... index of charged scalar boson, positiv
 !  j .......... index of the down-type sfermion \tilde{f}
 !  k .......... index of the up-type sfermion \conjugate{\tilde{f}}
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vevs (v_d, v_u, v_L3)
 !  bi(i) ...... bilinear superpotential parameters (mu, epsilon_3)
 !  yukL ....... tau-Yukawa coupling
 !  yukd ....... Yukawa coupling of the down-type sfermion
 !  Ad ......... trilinear couplings of down-type sfermion
 !  Rsfd(i,j) .. mixing matrix of the down-type sfermion
 !  yuku ....... Yukawa coupling of the up-type sfermion
 !  Au ......... trilinear couplings of up-type sfermion
 !  Rsfu(i,j) .. mixing matrix of the down-type sfermion
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 15.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(4,4),RSfd(2,2),RSfu(2,2)  &
                           &   ,yukd,yuku,Ad,Au,bi(2),yukL
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: g,vevs(3)
  Integer :: i,j,k

  Integer :: i1,i2
  Complex(dp) :: parts(2,2)
  Real(dp) :: g2,yukd2,yuku2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion3Eps1'

  If ((i.Lt.1).Or.(i.Gt.4)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: i = ',i
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'down-type sfermion index j out of range: j = ',j
   Call TerminateProgram
  End If
  If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'up-type sfermion index k out of range: k = ',k
   Call TerminateProgram
  End If

  coup = ZeroC
  parts = ZeroC

  g2 = g**2
  yukd2 = Abs( yukd )**2
  yuku2 = Abs( yuku )**2

  parts(1,1) = vevs(1) * Rspm(i,1) * (0.5_dp * g2 - yukd2)  &
           & + vevs(2) * Rspm(i,2) * (0.5_dp * g2 - yuku2)  &
           & + vevs(3) * Rspm(i,3) * 0.5_dp * g2
  parts(1,1) = - parts(1,1) * oosqrt2
 
  parts(1,2) = yuku * ( Conjg( bi(1) ) * Rspm(i,1)        &
             &        + Au * Rspm(i,2)                     &
             &        - Conjg( bi(2) ) * Rspm(i,3) )

  parts(2,1) = Conjg( yukd ) * ( Conjg( Ad ) * RSpm(i,1) &
             &                  + bi(1) * RSpm(i,2)        &
             &                  - yukL * vevs(3) * RSpm(i,4) * oosqrt2 )

  parts(2,2) = yuku * Conjg(yukd) * (vevs(1) *Rspm(i,2) + vevs(2) *RSpm(i,1) )

  Do i1=1,2
   Do i2=1,2
    coup = coup + Conjg( Rsfd(j,i1) ) * parts(i1,i2) * Rsfu(k,i2)
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion3Eps1


 Subroutine CoupCharginoPseudoScalarEps1(i,j,k,U,V,RP0,Yuk,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseduoscalar bosons 
 ! valid for 1-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  Yuk ........ tau yukawa couplings: 
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RP0(3,3)
  Complex(dp), Intent(in) :: U(3,3),V(3,3),yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Complex(dp) :: yukC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarEps1'

  If ((i.Lt.1).Or.(i.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  yukC = Conjg( yuk )
  coupL = ( RP0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )       &
        &              - yuk * Conjg( V(i,3) * U(j,3) ) )   &
        & + g * RP0(k,2) * Conjg( V(i,2) * U(j,1) )         &
        & + RP0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )       &
        &              + yuk * Conjg( V(i,3) * U(j,2) ) ) ) * oosqrt2
  coupR = ( RP0(k,1) * ( g *  V(j,1) * U(i,2) - yukC *  V(j,3) * U(i,3) ) &
        & + g * RP0(k,2) *  V(j,2) * U(i,1)                               &
        & + RP0(k,3) * ( g * V(j,1) * U(i,3)                              &
        &              + yukC * V(j,3) * U(i,2)  ) ) * oosqrt2

  coupL = coupL * (0._dp,1._dp)
  coupR = coupR * (0._dp,-1._dp) 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarEps1


 Subroutine CoupCharginoScalarEps1(i,j,k,U,V,RS0,Yuk,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for 1-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  Yuk ........ tau yukawa couplings: 
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RS0(3,3)
  Complex(dp), Intent(in) :: U(3,3),V(3,3),yuk
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarEps1'

  If ((i.Lt.1).Or.(i.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = ( RS0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )        &
     &                 + yuk * Conjg( V(i,3) * U(j,3) ) )    &
     &    + g * RS0(k,2) * Conjg( V(i,2) * U(j,1) )          &
     &    + RS0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )        &
     &                 - yuk * Conjg( V(i,3) * U(j,2) ) )    &
     &    ) * oosqrt2
  coupR = ( RS0(k,1) * ( g * Conjg( V(j,1) * U(i,2) )        &
     &                 + yuk * Conjg( V(j,3) * U(i,3) ) )    &
     &    + g * RS0(k,2) * Conjg( V(j,2) * U(i,1) )          &
     &    + RS0(k,3) * ( g * Conjg( V(j,1) * U(i,3) )        &
     &                 - yuk * Conjg( V(j,3) * U(i,2) ) )    &
     &    ) * oosqrt2

  coupL = - coupL
  coupR = - Conjg( coupR )

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarEps1

 Subroutine CoupCSCharginoNeutralinoEps(k, i, j, N, U, V, RSpm, Yuk, gp, g &
                                       &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and
 ! charged scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model,
 ! 3-generation epsilon model
 !   k .......... index of the charged scalar
 !  i .......... index of chargino
 !  j .......... index of the neutralino
 !  N.. ........ mixing matrix of the neutralino
 !  U,V ........ mixing matrices of the chargino
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! the lagrangian is given by
 !  S^-(k) \bar{\chim(i)} (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 21.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k
  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), yuk(3)
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sumI, sumIC, ql(8), qr(8)
  Integer :: n_Spm,i1,n_neut,n_char

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCSCharginoNeutralinoEps'

  n_char = Size(U, Dim=1)
  n_neut = Size(N, Dim=1)
  n_Spm = Size(Rspm, Dim=1)

  If (.Not.( ((n_char.Eq.2).And.(n_neut.Eq.4).And.(n_Spm.Eq.2))    & ! MSSM
     &     .Or.((n_char.Eq.3).And.(n_neut.Eq.5).And.(n_Spm.Eq.4))  & ! 1-gen.
     &     .Or.((n_char.Eq.5).And.(n_neut.Eq.7).And.(n_Spm.Eq.8))  & ! 3-gen.
     &     ) ) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'model not consistent defined (n_char,n_neut,n_Spm)=' &
              & ,n_char,n_neut,n_Spm
   Call TerminateProgram
  End If

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Chargino index out of range (i,n_char) = ',i,n_char
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_neut)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Neutralino index out of range (i,n_neut) = ',j,n_neut
   Call TerminateProgram
  Else If ((k.Lt.1).Or.(k.Gt.n_Spm)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Scalar index out of range (i,n_Spm) = ',k,n_Spm
   Call TerminateProgram
  End If

  coupL = ZeroC
  coupR = ZeroC

  sumI = gp * N(j,1) + g * N(j,2)
  sumIC = Conjg(sumI)
  !------
  ! MSSM
  !------
  If (n_Spm.Eq.2) Then
   ql(2) = - (g * Conjg( V(i,1) * N(j,4)) + Conjg( V(i,2) ) * sumIC * oosqrt2)
   qr(1) = - g * U(i,1) * N(j,3) + U(i,2) * sumI * oosqrt2
   coupL = RSpm(k,2) * ql(2) 
   coupR = RSpm(k,1) * qr(1)
 
  !----------------------------
  ! 1-generation epsilon model
  !----------------------------
  Else If (n_Spm.Eq.4) Then
   ql(1) = Yuk(3) * Conjg( V(i,3) * N(j,5) )
   ql(2) = - (g * Conjg( V(i,1) * N(j,4)) + Conjg( V(i,2) ) * sumIC * oosqrt2)
   ql(3) = - Yuk(3) * Conjg( V(i,3) * N(j,3) )
   ql(4) = - sqrt2 * gp * Conjg( V(i,3) * N(j,1) )
   qr(1) = - g * U(i,1) * N(j,3) + U(i,2) * sumI * oosqrt2 
   qr(2) = ZeroC
   qr(3) = - g * U(i,1) * N(j,5) + U(i,3) * sumI * oosqrt2 
   qr(4) = Conjg(Yuk(3) )  * ( U(i,2) * N(j,5) - U(i,3) * N(j,3) )

   Do i1=1,4
    coupL = coupL + RSpm(k,i1) * ql(i1)
    coupR = coupR + RSpm(k,i1) * qr(i1)
   End Do
  !-----------------------------------------------
  ! 3-generation epsilon model
  !-----------------------------------------------
  Else If (n_Spm.Eq.8) Then
   ql(1) = Yuk(1) * Conjg( V(i,3) * N(j,5) )   &
       & + Yuk(2) * Conjg( V(i,4) * N(j,6) )   &
       & + Yuk(3) * Conjg( V(i,5) * N(j,7) )
   ql(2) = - (g * Conjg( V(i,1) * N(j,4)) + Conjg( V(i,2) ) * sumIC * oosqrt2)
   ql(3) = - Yuk(1) * Conjg( V(i,3) * N(j,3) )
   ql(4) = - Yuk(2) * Conjg( V(i,4) * N(j,3) )
   ql(5) = - Yuk(3) * Conjg( V(i,5) * N(j,3) )
   ql(6) = - sqrt2 * gp * Conjg( V(i,3) * N(j,1) )
   ql(7) = - sqrt2 * gp * Conjg( V(i,4) * N(j,1) )
   ql(8) = - sqrt2 * gp * Conjg( V(i,5) * N(j,1) )
   qr(1) = - g * U(i,1) * N(j,3) + U(i,2) * sumI * oosqrt2 
   qr(2) = ZeroC
   qr(3) = - g * U(i,1) * N(j,5) + U(i,3) * sumI * oosqrt2 
   qr(4) = - g * U(i,1) * N(j,6) + U(i,4) * sumI * oosqrt2 
   qr(5) = - g * U(i,1) * N(j,7) + U(i,5) * sumI * oosqrt2 
   qr(6) = Conjg(Yuk(1) ) * ( U(i,2) * N(j,5) - U(i,3) * N(j,3) )
   qr(7) = Conjg(Yuk(2) ) * ( U(i,2) * N(j,6) - U(i,4) * N(j,3) )
   qr(8) = Conjg(Yuk(3) ) * ( U(i,2) * N(j,7) - U(i,5) * N(j,3) )

   Do i1=1,8
    coupL = coupL + RSpm(k,i1) * ql(i1)
    coupR = coupR + RSpm(k,i1) * qr(i1)
   End Do
  Else
   Write(ErrCan,*) 'Error in Subroutine '//NameofUnit(Iname)
   Write(ErrCan,*) 'Model ',n_Spm,' is not defined.'
   Call TerminateProgram
  End If

 Iname = Iname - 1

 End Subroutine CoupCSCharginoNeutralinoEps

 Subroutine CoupScalarSfermion3Eps_1(i, j, k, RS0, T3, e, yuk, Rsf, A, bi, &
                      &          vevSM, vevL, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 15.9.1999
 ! 15.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2) , yuk , A, bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(:,:), T3, e, vevSM(2), vevL(:), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: n_S0,i1
  Complex(dp) :: ayuk,aA, abi(Size(bi))
  Real(dp) :: g2,gp2,YL,YR,Dterm, parts(Size(RS0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3Eps_1'

  coup = ZeroC

  n_S0 =  Size(RS0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,1) * Conjg( Rsf(j,1) ) &
        &         - YR * gp2 * Rsf(k,2) * Conjg( Rsf(j,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = ( yuk * abi(1) * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * bi(1) * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & -  vevSM(1) * Dterm

   parts(2) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                   &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2      &
          & + vevSM(2) * Dterm                                        &
          & - vevSM(2) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) ) &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

   Do i1=3,n_S0
    parts(i1) = - ( yuk * abi(i1-1) * Rsf(k,2) * Conjg(Rsf(j,1))             &
            &     + ayuk * bi(i1-1) * Rsf(k,1) * Conjg(Rsf(j,2)) ) * oosqrt2 &
            & -  vevL(i1-2) * Dterm
   End Do

  Else

   parts(1) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                    &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2       &
          & - vevSM(1) * Dterm                                         &
          & - vevSM(1) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) )  &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

   parts(2) = ( yuk * abi(1) * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * bi(1) * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & + vevSM(2) * Dterm

   parts(3:5) = - vevL * Dterm

  End If
  Do i1=1,n_S0
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3Eps_1

 Subroutine CoupScalarSfermion3Epsa_1(i, j, k, RS0, T3, e, yuk, Rsf, A, bi, &
                      &          vevSM, vevL, gU1, gSU2, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the MSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  bi(i) ...... bilinear parameters of the superpotential
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 15.9.1999
 ! 15.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2) , yuk , A, bi(:)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(:,:), T3, e, vevSM(2), vevL(:), gU1, gSU2
  Integer, Intent(in) :: i, j, k

  Integer :: n_S0,i1
  Complex(dp) :: ayuk,aA, abi(Size(bi))
  Real(dp) :: g2,gp2,YL,YR,Dterm, parts(Size(RS0, Dim=1))

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupScalarSfermion3Epsa_1'

  coup = ZeroC

  n_S0 =  Size(RS0, Dim=1) 

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  abi = Conjg( bi )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,1) * Conjg( id2c(j,1) ) &
        &         - YR * gp2 * Rsf(k,2) * Conjg( id2c(j,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = ( yuk * abi(1) * Rsf(k,2) * Conjg( id2c(j,1) )              &
          &   + ayuk * bi(1) * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2  &
          & -  vevSM(1) * Dterm

   parts(2) = - ( A * Rsf(k,2) * Conjg( id2c(j,1) )                   &
          &     + aA * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2      &
          & + vevSM(2) * Dterm                                        &
          & - vevSM(2) * yuk * ayuk * ( Rsf(k,1) * Conjg( id2c(j,1) ) &
          &                           + Rsf(k,2) * Conjg( id2c(j,2) ) )

   Do i1=3,n_S0
    parts(i1) = - ( yuk * abi(i1-1) * Rsf(k,2) * Conjg(id2c(j,1))             &
            &     + ayuk * bi(i1-1) * Rsf(k,1) * Conjg(id2c(j,2)) ) * oosqrt2 &
            & -  vevL(i1-2) * Dterm
   End Do

  Else

   parts(1) = - ( A * Rsf(k,2) * Conjg( id2c(j,1) )                    &
          &     + aA * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2       &
          & - vevSM(1) * Dterm                                         &
          & - vevSM(1) * yuk * ayuk * ( Rsf(k,1) * Conjg( id2c(j,1) )  &
          &                           + Rsf(k,2) * Conjg( id2c(j,2) ) )

   parts(2) = ( yuk * abi(1) * Rsf(k,2) * Conjg( id2c(j,1) )              &
          &   + ayuk * bi(1) * Rsf(k,1) * Conjg( id2c(j,2) ) ) * oosqrt2  &
          & + vevSM(2) * Dterm

   parts(3:5) = - vevL * Dterm

  End If
  Do i1=1,n_S0
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine CoupScalarSfermion3Epsa_1

 Subroutine CoupScalarWRP(i, g, vevs, vevL, RS0, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between scalar boson and W-boson
 ! valid for the MSSM + NMSSM
 !  i .......... index of scalar boson
 !  g .......... SU(2) gauge coupling
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  vevs ....... vevs of the neutral scalar fields: i=1 -> H_d
 !                                             i=2 -> H_u
 ! output
 !  coup ....... the coupling(i)
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

 Integer, Intent(in) :: i
 Real(dp), Intent(in) :: g, vevs(2), vevL(:), RS0(:,:)
 Real(dp), Intent(out) :: coup
 
 Integer :: n, i1

 n = size(vevL)

 coup = vevs(1) * RS0(i,1) + vevs(2) * RS0(i,2) 

 Do i1=1,n
  coup = coup + vevL(i1) * RS0(i,i1+2)
 End do

 coup = 0.5_dp * g**2 * coup

 End Subroutine CoupScalarWRP




 Subroutine CoupScalarZRP(i,g,cosW2,vevSM,vevL,RS0,coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between pseudoscalar boson, scalar boson, and
 ! Z-boson
 ! valid for the 1-generation epsilon model, and 3-generation epsilon model
 !  i .......... index of pseudoscalar boson
 !  g .......... SU(2) gauge coupling
 !  sinW2 ...... sin(weak mixing angle) squared
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  vevSM ...... vevs of the neutral scalar fields: i=1 -> H_d
 !                                                  i=2 -> H_u
 !  in case 1-generation model
 !  vevL ...... vevs of the neutral scalar fields: i=1 -> sneutrino_tau
 !  in case 3-generation model
 !  vevL ...... vevs of the neutral scalar fields: i=1 -> sneutrino_e
 !                                                 i=2 -> sneutrino_mu
 !                                                 i=3 -> sneutrino_tau
 ! output
 !  coup ....... the coupling(i,j)
 ! written by Werner Porod, 24.10.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g,cosW2,vevSM(2),RS0(:,:),vevL(:)
  Real(dp), Intent(out) :: coup
  Integer, Intent(in) :: i

  Integer :: i1,n_vev

  n_vev = Size( vevL )

  coup = vevSM(1) * RS0(i,1) + vevSM(2) * RS0(i,2)
  Do i1=1,n_vev
   coup = coup + RS0(i,2+i1) * vevL(i1)
  End Do

  coup = 0.5_dp * g**2 * coup / cosW2

 End Subroutine CoupScalarZRP

 Subroutine AllCouplingsEps3(g, Y_l, Y_d, uD_L, uD_R , Y_u, uU_L, uU_R        &
    & , vevSM, vevL, RSpm, RP0, RS0, Umat, Vmat, Nmat, bi, phase_glu          &
    & , Ae, RSup, Au, RSdown, Ad, GenerationMixing                            &
    & , cpl_SmpSdSu, cpl_SmpSuSd, cpl_SmpP03, cpl_SmpP0W, cpl_SmpS03          &
    & , cpl_SmpS0W, cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_DDZ_L    &
    & , cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L      &
    & , cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R, cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L &
    & , cpl_GDSd_R, cpl_DNSd_L, cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R            &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_DDP0_L, cpl_UUP0_L, cpl_DDP0_R            &
    & , cpl_UUP0_R, cpl_DDS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_UUS0_R            &
    & , cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R, cpl_GlGlS0            &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0S0Z, cpl_P0S03, cpl_S0SdSd, cpl_S0SuSu  &
    & , cpl_S03, cpl_S0WW, cpl_S0ZZ, cpl_SdSuW, cpl_SuSdW, cpl_SdSdZ          &
    & , cpl_SuSuZ, cpl_CCP0_L, cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L  &
    & , cpl_CNW_R, cpl_SmpCN_L, cpl_SmpCN_R)
 !----------------------------------------------------------------------------
 ! Routine for calculating all couplings of the three-generation espilon model
 ! output:
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 !  in the case that fermions are involved, _L and _R refer to left- and
 !  right couplings
 ! written by Werner Porod, 
 ! 05.09.2001: taking AllCouplingsMSSM as basis
 ! 12.09.2002: changing interface to avoid global variables
 !-----------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3) ! mixing matrix of left d-quarks
  Complex(dp), Intent(in) :: uD_R(3,3) ! mixing matrix of right d-quarks
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3) ! mixing matrix of left u-quarks
  Complex(dp), Intent(in) :: uU_R(3,3) ! mixing matrix of right u-quarks
  Real(dp), Intent(in) :: vevSM(2)     ! MSSM Higgs vevs [v_d, v_u]
  Real(dp), Intent(in) :: vevL(3)      ! sneutrino vevs [e, mu, tau]
  Complex(dp), Intent(in) :: RSpm(8,8) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(5,5)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(5,5)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(5,5), Vmat(5,5) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(7,7) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  bi(4)   ! superpotential bilinears [mu, epsilon_i]
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters
  Logical, Intent(in) :: GenerationMixing ! if .true. generation mixing of
                                          ! (s)fermions is taken into account

  Complex(dp), Intent(out) :: cpl_SmpSdSu(8,6,6), cpl_SmpSuSd(8,6,6)       &
      & , cpl_SmpP03(8,8,5), cpl_SmpP0W(8,5), cpl_SmpS03(8,8,5)            &
      & , cpl_SmpS0W(8,5), cpl_SmpDU_L(8,3,3), cpl_SmpDU_R(8,3,3)          &
      & , cpl_SmpZ(8,8), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(5,5), cpl_CCZ_R(5,5)           &
      & , cpl_NNZ_L(7,7), cpl_NNZ_R(7,7), cpl_NNS0_L(7,7,5)            &
      & , cpl_NNS0_R(7,7,5), cpl_NNP0_L(7,7,5), cpl_NNP0_R(7,7,5) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,7,6), cpl_DNSd_R(3,7,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,7,6), cpl_UNSu_R(3,7,6)         &
      & , cpl_DDP0_L(3,3,5), cpl_UUP0_L(3,3,5), cpl_DDP0_R(3,3,5)       &
      & , cpl_UUP0_R(3,3,5), cpl_DDS0_L(3,3,5), cpl_UUS0_L(3,3,5)       &
      & , cpl_DDS0_R(3,3,5), cpl_UUS0_R(3,3,5)
  Complex(dp), Intent(out) :: cpl_CUSd_L(5,3,6), cpl_CUSd_R(5,3,6)      &
      & , cpl_CDSu_L(5,3,6), cpl_CDSu_R(5,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(5)
  Complex(dp) :: cpl_P0SdSd(5,6,6), cpl_P0SuSu(5,6,6), cpl_P0S0Z(5,5) 
  Real(dp), Intent(out) :: cpl_P0S03(5,5,5)
  Complex(dp), Intent(out) :: cpl_S0SdSd(5,6,6), cpl_S0SuSu(5,6,6)
  Real(dp), Intent(out) :: cpl_S03(5,5,5), cpl_S0WW(5), cpl_S0ZZ(5)
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6) &
      & , cpl_SdSdZ(6,6), cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(5,5,5), cpl_CCP0_R(5,5,5)    &
      & , cpl_CCS0_L(5,5,5), cpl_CCS0_R(5,5,5), cpl_CNW_L(5,7)        &
      & , cpl_CNW_R(5,7), cpl_SmpCN_L(8,5,7), cpl_SmpCN_R(8,5,7)

  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3, i4
  Real(dp) :: gU1, gSU2, gSU3, e_d, e_u, sinW2, cosW2, cosW
  Complex(dp) :: Rsd(2,2), Rsu(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, mu, YukL(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsEps3'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 5
  n_neut = 7
  n_P0 = 5
  n_S0 = 5
  n_Spm = 8

  !--------------------
  ! some constants
  !--------------------
  e_u = 2._dp / 3._dp
  e_d = - 1._dp / 3._dp

  Do i1=1,3
   yukL(i1)=y_L(i1,i1)
  End Do
  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm, yukL  &
               &, gU1, gSU2, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_D, uD_L, uD_R, Y_U &
              &, uU_L, uU_R, cpl_SmpDU_L(i1,i2,i3), cpl_SmpDU_R(i1,i2,i3) )
     End Do
    End Do
   End Do
  Else
   Do i1=1,n_Spm
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  End If
  !------------------------------
  ! charged scalar - pseudoscalar - W
  !------------------------------
  cpl_SmpP0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  cpl_SmpP03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_P0
!     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
!                                       &, cpl_SmpP03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  cpl_SmpS0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------
  ! charged scalar - charged scalar - scalar
  !------------------------------------------
  cpl_SmpS03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_S0
!     Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
!                                  &, cpl_SmpS03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSdSu = ZeroC
  cpl_SmpSuSd = ZeroC

!  If (GenerationMixing) Then
!   Do i1 = 1, n_Spm
!    Do i2 = 1,6
!     Do i3 = 1,6
!      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
!                                  & , Y_d, Ad, Rsdown, Y_u, Au, Rsup, coupC)
!      cpl_SmpSdSu(i1, i2, i3) = coupC
!      cpl_SmpSuSd(i1, i3, i2) = Conjg(coupC)
!     End Do
!    End Do
!   End Do

!  Else
!   Do i1 =1,n_Spm
!    Do i2=1,3
!     Yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Yukp = Y_u(i2,i2)
!     Ap = Au(i2,i2)
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
!                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
!       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
!       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
!      End Do
!     End Do
!    End Do ! i2
!   End Do ! i1
!  End If

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  cpl_SmpZ = 0._dp
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, YukL, gSU2  &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, YukL, gSU2  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_char
    Do i2=1,3
     Do i3=1,6
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                              &, Y_u, uU_L, uU_R, Umat, Vmat, coupLC, coupRC)
      cpl_CUSd_L(i1, i2, i3) = coupLC
      cpl_CUSd_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                              &, Y_u, uD_L, uD_R, Umat, Vmat, coupLC, coupRC)
      cpl_CDSu_L(i1, i2, i3) = coupLC
      cpl_CDSu_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,n_char
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do
  End If 

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_P0
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RP0 &
                                 &, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_d(i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_S0
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RS0 &
                           &, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RS0 &
                           &, cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, -0.5_dp, Y_d(i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  If (GenerationMixing) Then
   cpl_DUW = - CKM * gSU2 * oosqrt2
  Else
   cpl_DUW = - id3C * gSU2 * oosqrt2
  End If
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,n_neut
     Do i3 = 1,6
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, Nmat, coupLC, coupRC)
      cpl_DNSd_L(i1, i2, i3 ) = coupLC
      cpl_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, Nmat, coupLC, coupRC)
      cpl_UNSu_L(i1, i2, i3 ) = coupLC
      cpl_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsdown, uD_L, uD_R, &
                          & coupLC, coupRC)
     cpl_GDSd_L(i1, i2) = coupLC
     cpl_GDSd_R(i1, i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsup, uU_L, uU_R, &
                          & coupLC, coupRC)
     cpl_GUSu_L(i1, i2) = coupLC
     cpl_GUSu_R(i1, i2) = coupRC
    End Do
   End Do
  Else

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  End If 

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, Nmat, RP0, gU1, gSU2, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalar(i1, i2, i3, Nmat, RS0, gU1, gSU2, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  If (GenerationMixing) Then
   Do i1 = 1,6
    Do i2 = 1,6
     Call CoupSfermionW3(i1, i2, gSU2, RSdown, RSup, cpl_SdSuW(i1,i2) )
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do
   End Do
  End If

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
!  Do i1=1,n_P0
!   Do i2=1,n_P0
!    Do i3=1,n_S0
!     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
!                                &, cpl_P0S03(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
   End Do
  End Do
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp

!  If (GenerationMixing) Then
!   Do i1=1,n_P0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_d, Rsdown   &
!                                   &, Ad, bi, cpl_P0SdSd(i1,i2,i3) )
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, 0.5_dp, Y_u, Rsup      &
!                                   &, Au, bi, cpl_P0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
!
!  Else
!   Do i1=1,n_P0
!    Do i2=1,3
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
!                                   &, A, bi, coupC )
!       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!     yuk = Y_u(i2,i2)
!     A = Au(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
!                                   &, A, bi, coupC )
!       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!    End Do
!   End Do
!
!  End If

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
!  Do i1=1,n_S0
!   Do i2=1,n_S0
!    Do i3=1,n_S0
!     Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
!    End Do
!   End Do
!  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp

  If (GenerationMixing) Then
!   Do i1=1,n_S0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d, Rsdown   &
!                           &, Ad, mu, vevs, gU1, gSU2, cpl_S0SdSd(i1,i2,i3) )
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u, Rsup      &
!                           &, Au, mu, vevs, gU1, gSU2, cpl_S0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!
  Else
   Do i1=1,n_S0
    Do i2=1,3
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
                               &, A, bi, vevSM, vevL, gU1, gSU2, coupC )
       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
                               &, A, bi, vevSM, vevL, gU1, gSU2, coupC )
       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

  End If

  !----------
  ! scalar W
  !----------
  cpl_S0WW = 0._dp
  Do i1 = 1,n_S0
   Call CoupScalarW(i1, gSU2, vevSM, vevL, RS0, cpl_S0WW(i1) )
  End Do

  !----------
  ! scalar Z
  !----------
  cpl_S0ZZ = 0._dp
  Do i1=1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevSM, vevL, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SuSuZ = 0._dp

  If (GenerationMixing) Then
   Do i1=1,6
    Do i2=1,6
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, cpl_SdSdZ(i1, i2) )
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, cpl_SuSuZ(i1, i2) )
    End Do
   End Do
 
  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  End If

  Iname = Iname - 1

 End Subroutine AllCouplingsEps3

 Subroutine CoupChargedScalarSfermion3Eps3(i, j, k, RSpm, g, vevs, bi    &
                       &      ,yukL, yukd, Ad, Rsfd, yuku, Au, Rsfu, coup)
 !-----------------------------------------------------------------------
 ! calculates the 3-coupling between charged scalar and sfermions
 ! valid for the 3-generation epsilon model
 !  i .......... index of charged scalar boson, positiv
 !  j .......... index of the down-type sfermion \tilde{f}
 !  k .......... index of the up-type sfermion \conjugate{\tilde{f}}
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  vevs(i) .... vevs (v_d, v_u, v_L1, v_L2, v_L3)
 !  bi(i) ...... bilinear superpotential parameters
 !                             (mu, epsilon_1, epsilon_2, epsilon_3
 !  yukL(i) .... lepton-Yukawa couplings
 !  yukd ....... Yukawa coupling of the down-type sfermion
 !  Ad ......... trilinear couplings of down-type sfermion
 !  Rsfd(i,j) .. mixing matrix of the down-type sfermion
 !  yuku ....... Yukawa coupling of the up-type sfermion
 !  Au ......... trilinear couplings of up-type sfermion
 !  Rsfu(i,j) .. mixing matrix of the down-type sfermion
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 15.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSpm(8,8),RSfd(2,2),RSfu(2,2)  &
                           &   ,yukd,yuku,Ad,Au,bi(4),yukL(3)
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: g,vevs(5)
  Integer :: i,j,k

  Integer :: i1,i2
  Complex(dp) :: parts(2,2)
  Real(dp) :: g2,yukd2,yuku2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChargedScalarSfermion3Eps3'

  If ((i.Lt.1).Or.(i.Gt.8)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range: i = ',i
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'down-type sfermion index j out of range: j = ',j
   Call TerminateProgram
  End If
  If ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'up-type sfermion index k out of range: k = ',k
   Call TerminateProgram
  End If

  coup = ZeroC
  parts = ZeroC

  g2 = g**2
  yukd2 = Abs( yukd )**2
  yuku2 = Abs( yuku )**2

  parts(1,1) = vevs(1) * Rspm(i,1) * (0.5_dp * g2 - yukd2)   &
           & + vevs(2) * Rspm(i,2) * (0.5_dp * g2 - yuku2)   &
           & + 0.5_dp * ( vevs(3) * Rspm(i,3)                &
           &           + vevs(4) * Rspm(i,4)                &
           &           + vevs(5) * Rspm(i,5) )
  parts(1,1) = - parts(1,1) * oosqrt2

  parts(1,2) = yuku * ( Conjg( bi(1) ) * Rspm(i,1)  &
             &        + Au * Rspm(i,2)               &
             &        - Conjg( bi(2) ) * Rspm(i,3)  &
             &        - Conjg( bi(3) ) * Rspm(i,4)  &
             &        - Conjg( bi(4) ) * Rspm(i,5) )

  parts(2,1) = Conjg( yukd ) * ( Conjg( Ad ) * RSpm(i,1)         &
             &                  + bi(1) * RSpm(i,2)                &
             &                  - ( yukL(1) * vevs(3) * RSpm(i,6)  &
             &                    + yukL(2) * vevs(4) * RSpm(i,7)  &
             &                    + yukL(3) * vevs(5) * RSpm(i,8) ) * oosqrt2 )

  parts(2,2) = yuku * Conjg(yukd) * (vevs(1) *Rspm(i,2) + vevs(2) *RSpm(i,1) )

  Do i1=1,2
   Do i2=1,2
    coup = coup + Conjg( Rsfd(j,i1) ) * parts(i1,i2) * Rsfu(k,i2)
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CoupChargedScalarSfermion3Eps3



 Subroutine CoupCharginoPseudoScalarEps3(i,j,k,U,V,RP0,Yuk,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for 3-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RP0(5,5)
  Complex(dp), Intent(in) :: U(5,5),V(5,5),yuk(3)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Complex(dp) :: yukC(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarEps3'

  If ((i.Lt.1).Or.(i.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If

  yukC = Conjg( yuk ) 

  coupL = ( RP0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )         &
        &              - yuk(1) * Conjg( V(i,3) * U(j,3) )    &
        &              - yuk(2) * Conjg( V(i,4) * U(j,4) )    &
        &              - yuk(3) * Conjg( V(i,5) * U(j,5) ) )  &
        & + g * RP0(k,2) * Conjg( V(i,2) * U(j,1) )           &
        & + RP0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )         &
        &              + yuk(1) * Conjg( V(i,3) * U(j,2) ) )  &
        & + RP0(k,4) * ( g * Conjg( V(i,1) * U(j,4) )         &
        &              + yuk(2) * Conjg( V(i,4) * U(j,2) ) )  &
        & + RP0(k,5) * ( g * Conjg( V(i,1) * U(j,5) )         &
        &              + yuk(3) * Conjg( V(i,5) * U(j,2) ) ) ) * oosqrt2
  coupR = ( RP0(k,1) * ( g * V(j,1) * U(i,2)             &
        &              - yukC(1) * V(j,3) * U(i,3)       &
        &              - yukC(2) * V(j,4) * U(i,4)       &
        &              - yukC(3) * V(j,5) * U(i,5) )     &
        & + g * RP0(k,2) * V(j,2) * U(i,1)               &
        & + RP0(k,3) * ( g * V(j,1) * U(i,3)             &
        &              + yukC(1) *  V(j,3) * U(i,2) )    &
        & + RP0(k,4) * ( g *  V(j,1) * U(i,4)            &
        &              + yukC(2) * V(j,4) * U(i,2) )     &
        & + RP0(k,5) * ( g * V(j,1) * U(i,5)             &
        &              + yukC(3) * V(j,5) * U(i,2) ) ) * oosqrt2 

  coupL = coupL * (0._dp,1._dp)
  coupR = coupR * (0._dp,-1._dp) 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarEps3



 Subroutine CoupCharginoScalarEps3(i,j,k,U,V,RS0,Yuk,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for 3-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RS0(5,5)
  Complex(dp), Intent(in) :: U(5,5),V(5,5),yuk(3)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarEps3'

  If ((i.Lt.1).Or.(i.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = ( RS0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )             &
     &                 + yuk(1) * Conjg( V(i,3) * U(j,3) )        &
     &                 + yuk(2) * Conjg( V(i,4) * U(j,4) )        &
     &                 + yuk(3) * Conjg( V(i,5) * U(j,5) ) )      &
     &     + g * RS0(k,2) * Conjg( V(i,2) * U(j,1) )              &
     &     + RS0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )            &
     &                  - yuk(1) * Conjg( V(i,3) * U(j,2) ) )     &
     &     + RS0(k,4) * ( g * Conjg( V(i,1) * U(j,4) )            &
     &                 - yuk(2) * Conjg( V(i,4) * U(j,2) ) )      &
     &     + RS0(k,5) * ( g * Conjg( V(i,1) * U(j,5) )            &
     &                  - yuk(3) * Conjg( V(i,5) * U(j,2) ) )     &
     &     ) * oosqrt2
  coupR = ( RS0(k,1) * ( g * Conjg( V(j,1) * U(i,2) )             &
     &                 + yuk(1) * Conjg( V(j,3) * U(i,3) )        &
     &                 + yuk(2) * Conjg( V(j,4) * U(i,4) )        &
     &                 + yuk(3) * Conjg( V(j,5) * U(i,5) ) )      &
     &    + g * RS0(k,2) * Conjg( V(j,2) * U(i,1) )               &
     &    + RS0(k,3) * ( g * Conjg( V(j,1) * U(i,3) )             &
     &                 - yuk(1) * Conjg( V(j,3) * U(i,2) ) )      &
     &    + RS0(k,4) * ( g * Conjg( V(j,1) * U(i,4) )             &
     &                 - yuk(2) * Conjg( V(j,4) * U(i,2) ) )      &
     &    + RS0(k,5) * ( g * Conjg( V(j,1) * U(i,5) )             &
     &                 - yuk(3) * Conjg( V(j,5) * U(i,2) ) )      &
     &    ) * oosqrt2

  coupL = - coupL
  coupR = - Conjg( coupR )

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarEps3

 Subroutine AllCouplingsLam(g, Y_l, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R          &
    & , Ae, RSup, Au, RSdown, Ad, vevSM, vevL, lam, lamp                      &
    & , RSpm, RP0, RS0, Umat, Vmat, Nmat, bi, phase_glu, GenerationMixing     &
    & , cpl_SmpSdSu, cpl_SmpSuSd, cpl_SmpP03, cpl_SmpP0W, cpl_SmpS03          &
    & , cpl_SmpS0W, cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_DDZ_L    &
    & , cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L      &
    & , cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R, cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L &
    & , cpl_GDSd_R, cpl_DNSd_L, cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R            &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_DDP0_L, cpl_UUP0_L, cpl_DDP0_R            &
    & , cpl_UUP0_R, cpl_DDS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_UUS0_R            &
    & , cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R, cpl_GlGlS0            &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0S0Z, cpl_P0S03, cpl_S0SdSd, cpl_S0SuSu  &
    & , cpl_S03, cpl_S0WW, cpl_S0ZZ, cpl_SdSuW, cpl_SuSdW, cpl_SdSdZ          &
    & , cpl_SuSuZ, cpl_CCP0_L, cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L  &
    & , cpl_CNW_R, cpl_SmpCN_L, cpl_SmpCN_R)
 !-----------------------------------------------------------------
 ! Routine for calculating all couplings of the model with explicitly broken
 ! R-parity
 ! output:
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 ! written by Werner Porod, 
 ! 07.09.2001: taking AllCouplingsEps as basis
 ! 12.09.2002: changing interface to avoid global variables
 !-----------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3) ! mixing matrix of left d-quarks
  Complex(dp), Intent(in) :: uD_R(3,3) ! mixing matrix of right d-quarks
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3) ! mixing matrix of left u-quarks
  Complex(dp), Intent(in) :: uU_R(3,3) ! mixing matrix of right u-quarks
  Real(dp), Intent(in) :: vevSM(2)     ! MSSM Higgs vevs [v_d, v_u]
  Real(dp), Intent(in) :: vevL(3)      ! sneutrino vevs [e, mu, tau]
  Complex(dp), Intent(in) :: lam(3,3,3) ! lambda_ijk
  Complex(dp), Intent(in) :: lamp(3,3,3) ! lambda'_ijk
  Complex(dp), Intent(in) :: RSpm(8,8) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(5,5)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(5,5)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(5,5), Vmat(5,5) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(7,7) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  bi(4)   ! superpotential bilinears [mu, epsilon_i]
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters
  Logical, Intent(in) :: GenerationMixing ! if .true. generation mixing of
                                          ! (s)fermions is taken into account

  Complex(dp), Intent(out) :: cpl_SmpSdSu(8,6,6), cpl_SmpSuSd(8,6,6)       &
      & , cpl_SmpP03(8,8,5), cpl_SmpP0W(8,5), cpl_SmpS03(8,8,5)            &
      & , cpl_SmpS0W(8,5), cpl_SmpDU_L(8,3,3), cpl_SmpDU_R(8,3,3)          &
      & , cpl_SmpZ(8,8), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(5,5), cpl_CCZ_R(5,5)           &
      & , cpl_NNZ_L(7,7), cpl_NNZ_R(7,7), cpl_NNS0_L(7,7,5)            &
      & , cpl_NNS0_R(7,7,5), cpl_NNP0_L(7,7,5), cpl_NNP0_R(7,7,5) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,7,6), cpl_DNSd_R(3,7,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,7,6), cpl_UNSu_R(3,7,6)         &
      & , cpl_DDP0_L(3,3,5), cpl_UUP0_L(3,3,5), cpl_DDP0_R(3,3,5)       &
      & , cpl_UUP0_R(3,3,5), cpl_DDS0_L(3,3,5), cpl_UUS0_L(3,3,5)       &
      & , cpl_DDS0_R(3,3,5), cpl_UUS0_R(3,3,5)
  Complex(dp), Intent(out) :: cpl_CUSd_L(5,3,6), cpl_CUSd_R(5,3,6)      &
      & , cpl_CDSu_L(5,3,6), cpl_CDSu_R(5,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(5)
  Complex(dp) :: cpl_P0SdSd(5,6,6), cpl_P0SuSu(5,6,6), cpl_P0S0Z(5,5) 
  Real(dp), Intent(out) :: cpl_P0S03(5,5,5)
  Complex(dp), Intent(out) :: cpl_S0SdSd(5,6,6), cpl_S0SuSu(5,6,6)
  Real(dp), Intent(out) :: cpl_S03(5,5,5), cpl_S0WW(5), cpl_S0ZZ(5)
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6) &
      & , cpl_SdSdZ(6,6), cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(5,5,5), cpl_CCP0_R(5,5,5)    &
      & , cpl_CCS0_L(5,5,5), cpl_CCS0_R(5,5,5), cpl_CNW_L(5,7)        &
      & , cpl_CNW_R(5,7), cpl_SmpCN_L(8,5,7), cpl_SmpCN_R(8,5,7)

  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3 !, i4
  Real(dp) :: gU1, gSU2, gSU3, e_d, e_u, sinW2, cosW2, cosW
  Complex(dp) :: Rsd(2,2), Rsu(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, mu, YukL(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsLam'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 5
  n_neut = 7
  n_P0 = 5
  n_S0 = 5
  n_Spm = 8

  !--------------------
  ! some constants
  !--------------------
  e_u = 2._dp / 3._dp
  e_d = - 1._dp / 3._dp

  Do i1=1,3
   yukL(i1)=y_L(i1,i1)
  End Do
  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm, Y_L  &
             &, gU1, gSU2, lam, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_Spm
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_D, uD_L, uD_R, Y_U &
              &, uU_L, uU_R, cpl_SmpDU_L(i1,i2,i3), cpl_SmpDU_R(i1,i2,i3) )
     End Do
    End Do
   End Do
  Else
   Do i1=1,n_Spm
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  End If
  !------------------------------
  ! charged scalar - pseudoscalar - W
  !------------------------------
  cpl_SmpP0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  cpl_SmpP03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_P0
!     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
!                                       &, cpl_SmpP03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  cpl_SmpS0W = 0._dp
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------
  ! charged scalar - charged scalar - scalar
  !------------------------------------------
  cpl_SmpS03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_S0
!     Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
!                                  &, cpl_SmpS03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSdSu = ZeroC
  cpl_SmpSuSd = ZeroC

!  If (GenerationMixing) Then
!   Do i1 = 1, n_Spm
!    Do i2 = 1,6
!     Do i3 = 1,6
!      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
!                                  & , Y_d, Ad, Rsdown, Y_u, Au, Rsup, coupC)
!      cpl_SmpSdSu(i1, i2, i3) = coupC
!      cpl_SmpSuSd(i1, i3, i2) = Conjg(coupC)
!     End Do
!    End Do
!   End Do

!  Else
!   Do i1 =1,n_Spm
!    Do i2=1,3
!     Yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Yukp = Y_u(i2,i2)
!     Ap = Au(i2,i2)
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
!                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
!       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
!       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
!      End Do
!     End Do
!    End Do ! i2
!   End Do ! i1
!  End If

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  cpl_SmpZ = 0._dp
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, Y_L, lam, gSU2 &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, Y_L, lam, gSU2  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_char
    Do i2=1,3
     Do i3=1,6
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                    &, Y_u, uU_L, uU_R, lamp, Umat, Vmat, coupLC, coupRC)
      cpl_CUSd_L(i1, i2, i3) = coupLC
      cpl_CUSd_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                     &, Y_u, uD_L, uD_R, lamp, Umat, Vmat, coupLC, coupRC)
      cpl_CDSu_L(i1, i2, i3) = coupLC
      cpl_CDSu_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,n_char
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do
  End If 

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_P0
      Call CoupFermionPseudoScalar(i1, i2, i3, Y_d, lamp, uD_L, uD_R, RP0 &
                                 &, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, Y_d(i1,i1), lamp(i1,i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_S0
      Call CoupFermionScalar(i1, i2, i3, Y_d, lamp, uD_L, uD_R, RS0 &
                           &, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RS0 &
                           &, cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, Y_d(i1,i1), lamp(i1,i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  If (GenerationMixing) Then
   cpl_DUW = - CKM * gSU2 * oosqrt2
  Else
   cpl_DUW = - id3C * gSU2 * oosqrt2
  End If

  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,n_neut
     Do i3 = 1,6
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, Nmat, lamp, coupLC, coupRC)
      cpl_DNSd_L(i1, i2, i3 ) = coupLC
      cpl_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, Nmat, coupLC, coupRC)
      cpl_UNSu_L(i1, i2, i3 ) = coupLC
      cpl_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,n_neut
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsdown, uD_L, uD_R, &
                          & coupLC, coupRC)
     cpl_GDSd_L(i1, i2) = coupLC
     cpl_GDSd_R(i1, i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsup, uU_L, uU_R, &
                          & coupLC, coupRC)
     cpl_GUSu_L(i1, i2) = coupLC
     cpl_GUSu_R(i1, i2) = coupRC
    End Do
   End Do
  Else

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  End If 

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call CoupNeutralinoPseudoscalar(i1, i2, i3, Nmat, RP0, gU1, gSU2, &
                       & cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call CoupNeutralinoScalar(i1, i2, i3, Nmat, RS0, gU1, gSU2, &
                       & cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2) )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,6
    Do i2 = 1,6
     Call CoupSfermionW3(i1, i2, gSU2, RSdown, RSup, cpl_SdSuW(i1,i2) )
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do
   End Do
  End If

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
!  Do i1=1,n_P0
!   Do i2=1,n_P0
!    Do i3=1,n_S0
!     Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
!                                &, cpl_P0S03(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
    End Do
   End Do
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp

!  If (GenerationMixing) Then
!   Do i1=1,n_P0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_d, Rsdown   &
!                                   &, Ad, bi, cpl_P0SdSd(i1,i2,i3) )
!      Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, 0.5_dp, Y_u, Rsup      &
!                                   &, Au, bi, cpl_P0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!   End Do
!
!  Else
!   Do i1=1,n_P0
!    Do i2=1,3
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
!                                   &, A, bi, coupC )
!       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!     yuk = Y_u(i2,i2)
!     A = Au(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
!                                   &, A, bi, coupC )
!       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!    End Do
!   End Do
!
!  End If

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
!  Do i1=1,n_S0
!   Do i2=1,n_S0
!    Do i3=1,n_S0
!     Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
!    End Do
!   End Do
!  End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp

!  If (GenerationMixing) Then
!   Do i1=1,n_S0
!    Do i2=1,6
!     Do i3=1,6
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d, Rsdown   &
!                           &, Ad, mu, vevs, gU1, gSU2, cpl_S0SdSd(i1,i2,i3) )
!      Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u, Rsup      &
!                           &, Au, mu, vevs, gU1, gSU2, cpl_S0SuSu(i1,i2,i3) )
!     End Do
!    End Do
!
!  Else
!   Do i1=1,n_S0
!    Do i2=1,3
!     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
!
!     yuk = Y_d(i2,i2)
!     A = Ad(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
!                               &, A, mu, vevs, gU1, gSU2, coupC )
!       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!
!     yuk = Y_u(i2,i2)
!     A = Au(i2,i2)
!     Do i3=1,2
!      Do i4=1,2
!       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
!                               &, A, mu, vevs, gU1, gSU2, coupC )
!       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
!      End Do
!     End Do
!    End Do
!   End Do
!
!  End If

  !----------
  ! scalar W
  !----------
  cpl_S0WW = 0._dp
  Do i1 = 1,n_S0
   Call CoupScalarW(i1, gSU2, vevSM, vevL, RS0, cpl_S0WW(i1) )
  End Do

  !----------
  ! scalar Z
  !----------
  cpl_S0ZZ = 0._dp
  Do i1 = 1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevSM, vevL, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SuSuZ = 0._dp

  If (GenerationMixing) Then
   Do i1=1,6
    Do i2=1,6
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, cpl_SdSdZ(i1, i2) )
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, cpl_SuSuZ(i1, i2) )
    End Do
   End Do
 
  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  End If

  Iname = Iname - 1

 End Subroutine AllCouplingsLam

 Subroutine CoupChargedScalarFermion3Lam(i, j, k, RSpm, yukD, RdL, RdR, yukU &
                                    &, RuL, RuR, lamp, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charged scalars and fermions
 ! valid for the 3-generation MSSM, and the
 ! 3-generation epsilon model, and the model where R-parity is broken 
 ! spontaneously
 !  i ............. index of charged scalar 
 !  j ............. index of down-type fermion
 !  k ............. index of up-type fermion
 !  RSpm(i,j) ..... mixing matrix of charged scalar bosons
 !  yukD(i,j) ..... down type Yukawa coupling
 !  RdL(i,j) ...... left mixing matrix for down type fermions
 !  RdR(i,j) ...... right mixing matrix for down type fermions
 !  yukU(i,j) ..... up type Yukawa coupling
 !  RuL(i,j) ...... left mixing matrix for up type fermions
 !  RuR(i,j) ...... right mixing matrix for up type fermions
 !  lamp(i,j,k) ... trilinear couplings lambda_ijk L_i Q_j D^c_k
 !
 ! the lagrangian: H^-_i \bar(d_j) [ yuk_d R^*_i1 P_L +  yuk^*_u R_i2 P_L ] u_k
 !
 ! output
 !  coupL ....... the left coupling(i,j,k)
 !  coupR ....... the right coupling(i,j,k)
 ! written by Werner Porod, 01.05.2003
 !-----------------------------------------------------------------------
 Implicit None

 Complex(dp), Intent(in) :: RSpm(:,:), yukD(3,3), RdL(3,3), RdR(3,3) &
            & , yukU(3,3), RuL(3,3), RuR(3,3), lamp(3,3,3)
 Complex(dp), Intent(out) :: coupL, coupR
 Integer, Intent(in) :: i, j, k

 Integer :: n_Spm, i1, i2, n_test, i3

 n_Spm = Size(RSpm, Dim=1)

 Iname = Iname + 1
 NameOfUnit(Iname) = 'CoupChargedScalarFermion3Lam'

 If ((i.Lt.1).Or.(i.Gt.n_Spm)) Then
  Write(ErrCan,*) 'Problem in Subroutine '//NameOfUnit(Iname)
  Write(ErrCan,*) 'index i out of range: (i,n_Spm) = ',i,n_Spm 
  Call TerminateProgram
 End If

 coupL = 0._dp
 coupR = 0._dp

 n_test = 0
 If (CompareMatrices(RdL,id3C,NearlyZero) ) n_test = n_test + 1
 If (CompareMatrices(RdR,id3C,NearlyZero) ) n_test = n_test + 10
 If (CompareMatrices(RuL,id3C,NearlyZero) ) n_test = n_test + 100
 If (CompareMatrices(RuR,id3C,NearlyZero) ) n_test = n_test + 1000

 Select Case(n_test)
 Case (1)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1) * RSpm(i,2)
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) * RSpm(i,1) )
   End Do
  End Do

 Case (10)
  Do i1=1,3
    coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) * RSpm(i,1) )
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2) * RSpm(i,2)
   End Do
  End Do

 Case (11)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1) * RSpm(i,2)
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) * RSpm(i,1) )
  End Do

 Case (100)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) * RSpm(i,1) )
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2) * RSpm(i,2)
   End Do
  End Do

 Case (101)
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1) * RSpm(i,2)
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) * RSpm(i,1) )
  End Do

 Case (110)
  coupL = yukD(k,j) * Conjg( RSpm(i,1) )
  Do i1=1,3
   Do i2 = 1,3
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2) * RSpm(i,2)
   End Do
  End Do

 Case (111)
  coupL = yukD(k,j) * Conjg( RSpm(i,1) )
  Do i1=1,3
   coupR = coupR + Conjg( yukU(j,i1) ) * RuR(k,i1) * RSpm(i,2)
  End Do

 Case (1001)
  coupR = coupR + Conjg( yukU(j,k) ) * RSpm(i,2)
  Do i1=1,3
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) * RSpm(i,1) )
   End Do
  End Do

 Case (1010)
  Do i1=1,3
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) * RSpm(i,1) )
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1) * RSpm(i,2)
  End Do

 Case (1011)
  coupR = coupR + Conjg( yukU(j,k) ) * RSpm(i,2)
  Do i1=1,3
   coupL = coupL + yukD(i1,j) * Conjg( RuL(k,i1) * RSpm(i,1) )
  End Do

 Case (1100)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) * RSpm(i,1) )
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1) * RSpm(i,2)
  End Do

 Case (1101)
  coupR = Conjg( yukU(j,k) ) * RSpm(i,2)
  Do i1=1,3
   coupL = coupL + yukD(k,i1) * Conjg( RdR(j,i1) * RSpm(i,1) )
  End Do

 Case (1110)
   coupL = yukD(k,j) * Conjg( RSpm(i,1) )
  Do i1=1,3
   coupR = coupR + Conjg( yukU(i1,k) ) * RdL(j,i1) * RSpm(i,2)
  End Do

 Case (1111)
  coupL = yukD(k,j) * Conjg( RSpm(i,1) )
  coupR = Conjg( yukU(j,k) ) * RSpm(i,2)

 Case Default
  Do i1=1,3
   Do i2 = 1,3
    coupL = coupL + yukD(i1,i2) * Conjg( RdR(j,i2) * RuL(k,i1) * RSpm(i,1) )
    coupR = coupR + Conjg( yukU(i1,i2) ) * RdL(j,i1) * RuR(k,i2) * RSpm(i,2)
   End Do
  End Do

 End Select

 Do i1=1,3
  Do i2=1,3
   Do i3=1,3
    coupL = coupL &
        & + lamp(i1,i2,i3) * Conjg( RdR(j,i3) * RuL(k,i2) * RSpm(i,2+i1) )
   end do
  end do
 end do

 Iname = Iname - 1

 End Subroutine CoupChargedScalarFermion3Lam

 Subroutine CoupCharginoPseudoScalarLam3(i,j,k,U,V,RP0,Yuk,lam,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for 3-generation epsilon model including trilinear couplings
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  Yuk(i,j) ... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  lam(i,j,k) . trilinear coupling lambda LLE
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 3.11.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RP0(5,5)
  Complex(dp), Intent(in) :: U(5,5), V(5,5), yuk(3,3), lam(3,3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k


  Integer :: i1, i2, i3
  Complex(dp) :: part(2,5)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarLam3'

  If ((i.Lt.1).Or.(i.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If


  coupL = 0
  coupR = 0

  part(1,1) =  g * Conjg( V(i,1) * U(j,2) )
  part(2,1) =  g * V(j,1) * U(i,2)
  part(1,2) =  g * Conjg( V(i,2) * U(j,1) )
  part(2,2) =  g * V(j,2) * U(i,1)

  Do i1=1,3
   part(1,2+i1) = g * Conjg( V(i,1) * U(j,2+i1) )
   part(2,2+i1) = g * V(j,1) * U(i,2+i1)
   Do i2=1,3
    part(1,1) = part(1,1) - yuk(i1,i2) * Conjg( V(i,i2+2) * U(j,i1+2) )
    part(1,2+i1) = part(1,2+i1) + yuk(i1,i2) * Conjg( V(i,i2+2) * U(j,2) )
    part(2,1) = part(2,1) - Conjg( yuk(i1,i2) ) * V(j,i2+2) * U(i,i1+2)
    part(2,2+i1) = part(2,2+i1) + Conjg( yuk(i1,i2) ) * V(j,i2+2) * U(i,2)
    Do i3=1,3
     part(1,2+i1) = part(1,2+i1) + ( lam(i2,i1,i3) - lam(i1,i2,i3) ) &
                  &                * Conjg( V(i,i3+2) * U(j,i2+2) )
     part(2,2+i1) = part(2,2+i1) + Conjg( lam(i2,i1,i3) - lam(i1,i2,i3) ) &
                  &                * V(j,i3+2) * U(i,i2+2)
    End Do
   End Do
  End Do

  Do i1=1,5
   coupL = coupL + RP0(k,i1) * part(1,i1)
   coupR = coupR + RP0(k,i1) * part(2,i1)
  End Do

  coupL = coupL * (0._dp,1._dp) * oosqrt2
  coupR = coupR * (0._dp,-1._dp) * oosqrt2 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarLam3

 Subroutine CoupCharginoScalarLam3(i,j,k,U,V,RS0,Yuk,lam,g,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for 3-generation epsilon model including trilinear couplings
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  Yuk(i,j).... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  lam(i,j,k) . trilinear coupling lambda LLE
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RS0(5,5)
  Complex(dp), Intent(in) :: U(5,5), V(5,5), yuk(3,3), lam(3,3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2, i3
  Complex(dp) :: part(2,5)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarLam3'

  If ((i.Lt.1).Or.(i.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.5)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If

  coupL = 0
  coupR = 0

  part(1,1) =  g * Conjg( V(i,1) * U(j,2) )
  part(2,1) =  g * V(j,1) * U(i,2)
  part(1,2) =  g * Conjg( V(i,2) * U(j,1) )
  part(2,2) =  g * V(j,2) * U(i,1)

  Do i1=1,3
   part(1,2+i1) = g * Conjg( V(i,1) * U(j,2+i1) )
   part(2,2+i1) = g * V(j,1) * U(i,2+i1)
   Do i2=1,3
    part(1,1) = part(1,1) + yuk(i1,i2) * Conjg( V(i,i2+2) * U(j,i1+2) )
    part(1,2+i1) = part(1,2+i1) - yuk(i1,i2) * Conjg( V(i,i2+2) * U(j,2) )
    part(2,1) = part(2,1) + Conjg( yuk(i1,i2) ) * V(j,i2+2) * U(i,i1+2)
    part(2,2+i1) = part(2,2+i1) - Conjg( yuk(i1,i2) ) * V(j,i2+2) * U(i,2)
    Do i3=1,3
     part(1,2+i1) = part(1,2+i1) - ( lam(i2,i1,i3) - lam(i1,i2,i3) ) &
                  &                * Conjg( V(i,i3+2) * U(j,i2+2) )
     part(2,2+i1) = part(2,2+i1) - Conjg( lam(i2,i1,i3) - lam(i1,i2,i3) ) &
                  &                * V(j,i3+2) * U(i,i2+2)
    End Do
   End Do
  End Do

  Do i1=1,5
   coupL = coupL + RS0(k,i1) * part(1,i1)
   coupR = coupR + RS0(k,i1) * part(2,i1)
  End Do

  coupL = - coupL * oosqrt2
  coupR = - coupR * oosqrt2

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarLam3

 Subroutine CoupCharginoSfermion3Lam(i, j, k, g, T3, RSf, YukD, YukU, RfL &
                                &, RfR, lamp, U, V, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between chargino, fermion, and sfermion'
 ! valid for the 3-generation model with explicit R-parity breaking 
 !  i .......... index of chargino
 !  j .......... index of fermion
 !  k .......... index of sfermion
 !  g .......... SU(2) gauge coupling
 !  T3 ......... weak isospin of the fermion
 !  RSf(i,j) ... mixing matrix of the sfermion
 !  YukD ....... Yukawa coupling of the T3=-0.5 fermion 
 !  YukU ....... Yukawa coupling of the T3=0.5 fermion 
 !  U,V ........ mixing matrices of the chargino
 !  lamp(i,j,k) . lambda' coupling LQD
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding part of the Lagrangian is:
 !  \bar{f_j} (coupL P_L + coupR P_R) chi^\pm_i \tilde{f'}_k
 ! written by Werner Porod, 25.08.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, T3
  Complex(dp), Intent(in) :: RSf(6,6), U(:,:), V(:,:), lamp(3,3,3)
  Complex(dp), Intent(in), Dimension(3,3) :: YukD, YukU, RfL, RfR
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Integer :: n_char, i1, i2, i3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoSfermion3Lam' 

  n_char = Size(U, Dim=1)

  If (ErrorLevel.Ge.-1) Then
   If ((i.Lt.1).Or.(i.Gt.n_char)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index i out of range: (i,n_char) = ',i,n_char
    Call TerminateProgram
   End If 
   If ((j.Lt.1).Or.(j.Gt.3)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index j out of range: ',j
    Call TerminateProgram
   End If 
   If ((k.Lt.1).Or.(k.Gt.6)) Then
    Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
    Write(ErrCan,*) 'index k out of range: ',k
    Call TerminateProgram
   End If 
  End If

  coupL = 0._dp
  coupR = 0._dp

  If ( T3.Gt.0._dp) Then

   If ( CompareMatrices(RfL,id3C,NearlyZero) .And. &
      & CompareMatrices(RfR,id3C,NearlyZero)        ) Then
    coupR = - g  * Conjg( RSf(k,j) ) * U(i,1)
    Do i1=1,3
     coupL = coupL + YukU(i1,j) * Conjg( RSf(k,i1) )
     coupR = coupR + Conjg( YukD(j,i1) * RSf(k,i1+3) ) * U(i,2)
     Do i2=1,3
      coupR = coupR + Conjg(lamp(i2,i1,j) * RSf(k,3+i1)) * U(i,2+i2)
     End Do
    End Do

   Else If (CompareMatrices(RfL,id3C,NearlyZero)) Then
    coupR = - g  * Conjg( RSf(k,j) ) * U(i,1)
    Do i1=1,3
     coupR = coupR + Conjg( YukD(j,i1) * RSf(k,i1+3) ) * U(i,2)
     Do i2=1,3
      coupL = coupL + YukU(i1,i2) * Conjg( RSf(k,i1) ) * RfR(j,i2) 
      coupR = coupR + Conjg(lamp(i2,i1,j) * RSf(k,3+i1)) * U(i,2+i2)
     End Do
    End Do

   Else If (CompareMatrices(RfR,id3C,NearlyZero)) Then
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * U(i,1)
     coupL = coupL + YukU(i1,j) * Conjg( RSf(k,i1) )
     Do i2=1,3
      coupR = coupR + Conjg( YukD(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * U(i,2)
      Do i3=1,3
       coupR = coupR   &
           & + Conjg(lamp(i3,i1,i2) * RSf(k,3+i1)) * RfL(j,i2) * U(i,2+i3)
      End Do
     End Do
    End Do

   Else
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * U(i,1)
     Do i2=1,3
      coupL = coupL + YukU(i1,i2) * Conjg( RSf(k,i1) ) * RfR(j,i2) 
      coupR = coupR + Conjg( YukD(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * U(i,2)
      Do i3=1,3
       coupR = coupR   &
           & + Conjg(lamp(i3,i1,i2) * RSf(k,3+i1)) * RfL(j,i2) * U(i,2+i3)
      End Do
     End Do
    End Do

   End If
   coupL = coupL * Conjg( V(i,2) )

  Else  ! T3 < 0
 
   If ( CompareMatrices(RfL,id3C,NearlyZero) .And. &
      & CompareMatrices(RfR,id3C,NearlyZero)       ) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) ) * Conjg( U(i,2) )
     coupR = coupR + Conjg( YukU(j,i1) * RSf(k,i1+3) ) * V(i,2)
     Do i2=1,3
      coupL = coupL   &
          & + lamp(i2,i1,j) * Conjg(RSf(k,i1)) * Conjg(U(i,2+i2))
     End Do
    End Do

   Else If (CompareMatrices(RfL,id3C,NearlyZero)) Then
    coupR = - g * Conjg( RSf(k,j) ) * V(i,1)
    Do i1=1,3
     coupR = coupR + Conjg( YukU(j,i1) * RSf(k,i1+3) ) * V(i,2)
     Do i2=1,3
      coupL = coupL  &
          & + YukD(i1,i2) * Conjg( RSf(k,i1) ) * RfR(j,i2) * Conjg( U(i,2) )
      Do i3=1,3
       coupL = coupL   &
           & + lamp(i3,i1,i2) * Conjg(RSf(k,i1)) * RfR(j,i2) * Conjg(U(i,2+i3))
      End Do
     End Do
    End Do

   Else If (CompareMatrices(RfR,id3C,NearlyZero)) Then
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     coupL = coupL + YukD(i1,j) * Conjg( RSf(k,i1) ) * Conjg( U(i,2) )
     Do i2=1,3
      coupL = coupL   &
          & + lamp(i2,i1,j) * Conjg(RSf(k,i1)) * Conjg(U(i,2+i2))
      coupR = coupR + Conjg( YukU(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * V(i,2)
     End Do
    End Do

   Else
    Do i1=1,3
     coupR = coupR - g * Conjg( RSf(k,i1) ) * RfL(j,i1) * V(i,1)
     Do i2=1,3
      coupL = coupL   &
          & + YukD(i1,i2) * Conjg( RSf(k,i1) ) * RfR(j,i2) * Conjg( U(i,2) )
      Do i3=1,3
       coupL = coupL   &
           & + lamp(i3,i1,i2) * Conjg(RSf(k,i1)) * RfR(j,i2) * Conjg(U(i,2+i3))
      End Do
      coupR = coupR + Conjg( YukU(i1,i2) * RSf(k,i2+3) ) * RfL(j,i1) * V(i,2)
     End Do
    End Do
   End If

  End If

  Iname = Iname - 1

 End Subroutine CoupCharginoSfermion3Lam


 Subroutine CoupCSCharginoNeutralinoLam2(k, i, j, N, U, V, RSpm, Yuk, gp, g &
                                       &, lam, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and
 ! charged scalar bosons 
 ! valid for the 3 generation model with explicit R-parity breaking
 !   k ............ index of the charged scalar
 !  i ............ index of chargino
 !  j ............ index of the neutralino
 !  N.. .......... mixing matrix of the neutralino
 !  U,V .......... mixing matrices of the chargino
 !  RSpm(i,j) .... mixing matrix of charged scalar bosons
 !  Yuk(i,j) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  gp ........... U(1) gauge coupling
 !  g ............ SU(2) gauge coupling
 !  lam(i,j,k) ... lambda coupling LLE
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! the lagrangian is given by
 !  S^-(k) \bar{\chim(i)} (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 25.08.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k
  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), yuk(3,3) &
                          &, lam(3,3,3)   
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sumI, sumIC, ql(8), qr(8)
  Integer :: n_Spm, i1, i2, i3, n_neut, n_char

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCSCharginoNeutralinoLam2'

  n_char = Size(U, Dim=1)
  n_neut = Size(N, Dim=1)
  n_Spm = Size(Rspm, Dim=1)

  If (.Not.((n_char.Eq.5).And.(n_neut.Eq.7).And.(n_Spm.Eq.8)) ) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'model not consistent defined (n_char,n_neut,n_Spm)=' &
              & ,n_char,n_neut,n_Spm
   Call TerminateProgram
  End If

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Chargino index out of range (i,n_char) = ',i,n_char
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_neut)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Neutralino index out of range (i,n_neut) = ',j,n_neut
   Call TerminateProgram
  Else If ((k.Lt.1).Or.(k.Gt.n_Spm)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Scalar index out of range (i,n_Spm) = ',k,n_Spm
   Call TerminateProgram
  End If

  coupL = ZeroC
  coupR = ZeroC

  sumI = gp * N(j,1) + g * N(j,2)
  sumIC = Conjg(sumI)

  ql = zeroC
  qr = zeroC
  ql(2) = - (g * Conjg( V(i,1) * N(j,4)) + Conjg( V(i,2) ) * sumIC * oosqrt2)
  qr(1) = - g * U(i,1) * N(j,3) + U(i,2) * sumI * oosqrt2 
   Do i1=1,3
    ql(5+i1) = - sqrt2 * gp * Conjg( V(i,i1+2) * N(j,1) )
    qr(2+i1) = - g * U(i,1) * N(j,i1+4) + U(i,i1+2) * sumI * oosqrt2 
    Do i2=1,3
     ql(1) = ql(1) + Yuk(i1,i2) * Conjg( V(i,i2+2) * N(j,i1+4) )
     ql(2+i1) = ql(2+i1) - Yuk(i1,i2) * Conjg( V(i,i2+2) * N(j,3) )
     qr(5+i1) = qr(5+i1) + Conjg(Yuk(i2,i1)) * U(i,2) * N(j,i2+4)    &
              &          - Conjg(Yuk(i1,i2)) * U(i,i2+2) * N(j,3) 
     Do i3=1,3
      qr(5+i1) = qr(5+i1) + Conjg(lam(i2,i3,i1) - lam(i3,i2,i1))     &
               &            * U(i,i2+2) * N(j,i3+4) 
     End Do
    End Do
   End Do

   Do i1=1,8
    coupL = coupL + RSpm(k,i1) * ql(i1)
    coupR = coupR + RSpm(k,i1) * qr(i1)
   End Do

 Iname = Iname - 1

 End Subroutine CoupCSCharginoNeutralinoLam2

  Subroutine CoupCSCharginoNeutralinoLam(k1, k2, k3, gp, g, RSpm, N, U, V &
                                       &, lambda, Yl, coupL, coupR)
 !---------------------------------------------------------------------------
 ! Calculates the coupling between charged scalars, charginos and neutralinos
 ! in models with explicitly broken R-parity
 ! input:
 !  k1 ............. charged scalar index
 !  k2 ............. chargino index
 !  k3 ............. neutralino index
 !  gp ............. U(1) gauge coupling
 !  g .............. SU(2) gauge coupling
 !  RSpm(:,:) ...... charged scalar mixing matrix
 !  N(:,:) ......... neutralino mixing matrix
 !  U(:,:) ......... chargino mixing matrix
 !  V(:,:) ......... chargino mixing matrix
 !  lambda(3,3,3) .. trilinear coupling
 !  Yl(3,3) ........ lepton Yukawa coupling
 ! output:
 !  coupL,coupR
 ! the Lagrangian has the form
 !  \bar{chi^-_k2} ( coupL P_L + coupR P_R) chi^0_k3 S^-_k1
 ! written by Werner Porod : 16.01.03
 !--------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: k1, k2, k3
  Real(dp), Intent(in) :: gp, g
  Complex(dp), Intent(in) :: RSpm(:,:), N(:,:), U(:,:), V(:,:),        &
       &  lambda(3,3,3), Yl(3,3)
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: i, j, k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCSCharginoNeutralinoLam'

  coupL = - g * Conjg(N(k3, 4) * V(k2, 1)) * RSpm(k1, 2)                &
      & - oosqrt2 * Conjg( ( gp * N(k3, 1) +  g * N(k3, 2)) * V(k2, 2)) &
      &           * RSpm(k1, 2)

  coupR = - g * N(k3, 3) * U(k2, 1) * RSpm(k1, 1)                           &
      & + oosqrt2 * (gp * N(k3, 1) + g * N(k3, 2)) * U(k2, 2) * RSpm(k1, 1) 

  Do i=1,3
   coupL = coupL - gp * sqrt2 * Conjg(N(k3,1) * V(k2,2 + i)) * RSpm(k1, 5 + i)
   coupR = coupR -  g * N(k3, 4+i) * U(k2, 1) * RSpm(k1, 2+i)                &
      &  + oosqrt2 * (gp * N(k3,1) + g * N(k3,2)) * U(k2, 2+i) * RSpm(k1, 2+i)
   Do j=1,3
    coupL = coupL + Yl(i,j) * ( Conjg(N(k3,4+i) * V(k2,2 + j)) * RSpm(k1, 1) &
        &                     - Conjg( N(k3,3) * V(k2,2 + j)) * RSpm(k1,2+i) )
    coupR = coupR + Conjg(Yl(i,j)) * RSpm(k1, 5 + j)                         &
          &          * ( N(k3,4+i) * U(k2, 2) - N(k3, 3) * U(k2,2 + i)) 

    Do k=1,3
     coupL = coupL - Conjg(N(k3,4+i) * V(k2,2+k)) *  RSpm(k1,2+j) &
         &                              * lambda(i,j,k)  *2._dp
     coupR = coupR + N(k3, 4 + i) * U(k2, 2 + j) * RSpm(k1, 5 + k)        &
         &                        * Conjg(lambda(j, i, k)) * 2._dp 
    end do
   end do
  end do

 Iname = Iname - 1

 End Subroutine CoupCSCharginoNeutralinoLam


 Subroutine CoupDquarkPseudoScalar1expl(i,yuk,yukRP,RP0,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between d-quarks and pseudoscalar bosons 
 ! valid for 1-generation model with explicitly broken R-parity
 !  i .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings
 !  YukRP ...... lambda' yukawa couplings
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 24.10.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: RP0(:,:)
  Complex(dp), Intent(in) :: yuk,yukRP
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i
 
  coupL = - (0._dp,1._dp) * (yuk * RP0(i,1) + yukRP * RP0(i,3) ) * oosqrt2

  coupR = Conjg(coupL) 

 End Subroutine CoupDquarkPseudoScalar1expl



 Subroutine CoupDquarkPseudoScalar3expl(i,j,k,yuk,yukRP,Rdl,Rdr,RP0 &
                                   &   ,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between d-quarks and pseudoscalar bosons 
 ! valid for 3-generation model with explicitly broken R-parity
 !  i,j ........ generation index of the d-quarks
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  YukRP ...... lambda' yukawa couplings (3*3 matrix)
 !  Rdl ........ mixing matrix for the left-handed d-quarks
 !  Rdr ........ mixing matrix for the right-handed d-quarks
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 24.10.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: RP0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3),yukRP(3,3,3),Rdl(3,3),Rdr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i,j,k
 
  Integer :: i1, i2, i3
  Complex(dp) :: wert

  coupL = ZeroC
  Do i1=1,3
   Do i2=1,3
    wert = yuk(i1,i2) * RP0(k,1)
    Do i3=1,3
     wert = wert + yukRP(i3,i1,i2) * RP0(k,2+i3)
    End Do
    coupL = coupL &
        & + Conjg(Rdl(i,i1)) * wert * Conjg(Rdr(j,i2))
   End Do
  End Do
  coupL = - (0._dp,1._dp) * coupL * oosqrt2

  coupR = Conjg(coupL) 

 End Subroutine CoupDquarkPseudoScalar3expl



 Subroutine CoupDquarkScalar1expl(i,yuk,yukRP,RS0,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between d-quarks and scalar bosons 
 ! valid for 1-generation model with explicitly broken R-parity
 !  i .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings
 !  YukRP ...... lambda' yukawa couplings
 !  RS0(i,j) ... mixing matrix of scalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 24.10.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: RS0(:,:)
  Complex(dp), Intent(in) :: yuk,yukRP
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i
 
  coupL = - (yuk * RS0(i,1) + yukRP * RS0(i,3)) * oosqrt2

  coupR = Conjg(coupL) 

 End Subroutine CoupDquarkScalar1expl



 Subroutine CoupDquarkScalar3expl(i,j,k,yuk,yukRP,Rdl,Rdr,RS0,coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between d-quarks and scalar bosons 
 ! valid for 3-generation model with explicitly broken R-parity
 !  i,j ........ generation index of the d-quarks
 !  k .......... index of scalar boson
 !  Yuk ........ d-quark yukawa couplings (3*3 matrix)
 !  YukRP ...... lambda' yukawa couplings (3*3 matrix)
 !  Rdl ........ mixing matrix for the left-handed d-quarks
 !  Rdr ........ mixing matrix for the right-handed d-quarks
 !  RS0(i,j) ... mixing matrix of scalar bosons
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 24.10.2000
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: RS0(:,:)
  Complex(dp), Intent(in) :: yuk(3,3),yukRP(3,3,3),Rdl(3,3),Rdr(3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i,j,k
 
  Integer :: i1, i2, i3
  Complex(dp) :: wert

  coupL = ZeroC
  Do i1=1,3
   Do i2=1,3
    wert = yuk(i1,i2) * RS0(k,1)
    Do i3=1,3
     wert = wert + yukRP(i3,i1,i2) * RS0(k,2+i3)
    End Do
    coupL = coupL &
        & + Conjg(Rdl(i,i1)) * wert * Conjg(Rdr(j,i2))
   End Do
  End Do
  coupL = - coupL * oosqrt2

  coupR = Conjg(coupL) 

 End Subroutine CoupDquarkScalar3expl

 Subroutine CoupNeutralinoSdown3Expl(i, j, k, gp, g, RSf, RfL, RfR, Yuk, N &
                                  &, lamp, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between neutralino, down-quark, and sdowns
 ! valid for the 3-generation explicit R-parity breaking with lambda'
 !  i ............ index of down quark
 !  j ............ index of neutralino
 !  k ............ index of sdown
 !  g ............ SU(2) gauge coupling
 !  gp ........... U(1) gauge coupling
 !  RfL(i,j) ..... mixing matrix of the left-handed down quarks
 !  RfR(i,j) ..... mixing matrix of the right-handed down quarks
 !  RSf(i,j) ..... mixing matrix of the sdowns
 !  Yuk .......... down Yukawa coupling 
 !  N ............ mixing matrix of the neutralino
 !  lamp(a,b,c) .. trilinear R-partiy breaking coupling lambda'
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! the corresponding Lagrangian is given by:
 !    L = \bar(f_i) ( coupL P_L + coupR P_R) chi^0_j \tilde(f)_k  
 ! written by Werner Porod, 19.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in), Dimension(3,3) :: RfL, RfR, Yuk
  Complex(dp), Intent(in) :: RSf(6,6), N(7,7), lamp(3,3,3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Integer :: i1, i2, i3
  Complex(dp) :: fl, fr, hl(3,3), hr(3,3)
  Real(dp) :: e, T3

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupNeutralinoSdown3Expl'

  T3 = -0.5_dp
  e = -1._dp / 3._dp

  hl = - Yuk * Conjg( N(j,3) )


  Do i1=1,3
   Do i2=1,3
    Do i3=1,3
     hl(i1,i2) = hl(i1,i2) - Conjg(N(j,4+i3)) * lamp(i3, i1, i2)
    End Do
   End Do
  End Do
  hr = Conjg(hl)

  fL = - sqrt2 *( (e-T3) * gp * N(j,1) + T3 * g * N(j,2) )
  fR = sqrt2 * e * gp * Conjg( N(j,1) ) 

  coupL = 0._dp
  coupR = 0._dp

  If ( (CompareMatrices(id3C,RfR,NearlyZero)) .And. &
     & (CompareMatrices(id3C,RfL,NearlyZero))       )Then

   coupL = Conjg( Rsf(k,i+3) ) * fr
   coupR = Conjg( Rsf(k,i)   ) * fl
   Do i1=1,3
    coupL = coupL + Conjg( Rsf(k,i1) ) * hl(i1,i)
    coupR = coupR + Conjg( Rsf(k,i1+3) ) * hr(i,i1)
   end do

  Else If (CompareMatrices(id3C,RfL,NearlyZero)) Then
   coupR = coupR + Conjg( Rsf(k,i) ) * fl
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1+3) * RfR(i,i1) ) * fr
    coupR = coupR + Conjg( Rsf(k,i1+3)) * hr(i,i1)
    Do i2 = 1,3
     coupL = coupL + Conjg( Rsf(k,i1) * RfR(i,i2) ) * hl(i1,i2) 
    End Do
   End Do

  Else If (CompareMatrices(id3C,RfR,NearlyZero)) Then
   coupL = coupL + Conjg( Rsf(k,i+3) ) * fr
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1) ) * hl(i1,i) 
    coupR = coupR + Conjg( Rsf(k,i1) ) * RfL(i,i1) * fl
    Do i2 = 1,3
     coupR = coupR + Conjg( Rsf(k,i2+3) ) * RfL(i,i1) * hr(i1,i2)
    End Do
   End Do

  Else
   Do i1 = 1,3
    coupL = coupL + Conjg( Rsf(k,i1+3) * RfR(i,i1) ) * fr
    coupR = coupR + Conjg( Rsf(k,i1) ) * RfL(i,i1) * fl
    Do i2 = 1,3
     coupL = coupL + Conjg( Rsf(k,i1) * RfR(i,i2) ) * hl(i1,i2) 
     coupR = coupR + Conjg( Rsf(k,i2+3)) * RfL(i,i1) * hr(i1,i2)
    End Do
   End Do

  End If

  Iname = Iname - 1

 End Subroutine CoupNeutralinoSdown3Expl

 Subroutine CoupChi0P0(i, j, k, N, RP0, gp, g, h0, h, Y_nu, lam, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RP0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u Phi
 !  h .......... trilinear coupling Phi nu^C S
 !  Y_nu ....... neutrino Yukawa coupling
 !  lam ........ self coupling of the phi field
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 10.09.2004
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,RP0(8,8),gp, h0, h, lam
  Complex(dp), Intent(in) :: N(10,10), Y_nu(3)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_P0, n_N, i1
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChi0P0'

  n_P0 = Size(RP0, Dim=1)
  n_N = Size(N, Dim=1)

  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_P0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = -0.5_dp * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = -0.5_dp * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  !---------------------------------
  ! gauge part
  !---------------------------------
   coupL = sumI_i * ( Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2)  &
      &             + Conjg(N(j,5)) * RP0(k,3) + Conjg(N(j,6)) * RP0(k,4)  &
      &             + Conjg(N(j,7)) * RP0(k,5) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RP0(k,1) - Conjg(N(i,4)) * RP0(k,2)  &
      &             + Conjg(N(i,5)) * RP0(k,3) + Conjg(N(i,6)) * RP0(k,4)  &
      &             + Conjg(N(i,7)) * RP0(k,5) )

  !---------------------------------
  ! Yukawa part
  !---------------------------------
  coupL = coupL + oosqrt2 * h0                                             &
      &   * ( (Conjg(N(i,3)*N(j,4)) + Conjg(N(j,3)*N(i,4))) * RP0(k,6)     &
      &     + (Conjg(N(i,3)*N(j,10)) + Conjg(N(j,3)*N(i,10))) * RP0(k,2)   &
      &     + (Conjg(N(i,4)*N(j,10)) + Conjg(N(j,4)*N(i,10))) * RP0(k,1) ) &
      & - oosqrt2 * h                                                      &
      &   * ( (Conjg(N(i,8)*N(j,9)) + Conjg(N(j,8)*N(i,9))) * RP0(k,6)     &
      &     + (Conjg(N(i,8)*N(j,10)) + Conjg(N(j,8)*N(i,10))) * RP0(k,7)   &
      &     + (Conjg(N(i,9)*N(j,10)) + Conjg(N(j,9)*N(i,10))) * RP0(k,8) )
  Do i1=1,3
   coupL = coupL - oosqrt2 * Y_nu(i1)                                         &
      &   * ( (Conjg(N(i,4)*N(j,i1+4)) + Conjg(N(j,4)*N(i,i1+4))) * RP0(k,8)    &
      &     + (Conjg(N(i,4)*N(j,8)) + Conjg(N(j,4)*N(i,8))) * RP0(k,2+i1)      &
      &     + (Conjg(N(i,i1+4)*N(j,8)) + Conjg(N(j,i1+4)*N(i,8))) * RP0(k,2) )
  End Do
  coupL = coupL - oosqrt2 * lam * Conjg(N(i,10)*N(j,10)) * RP0(k,6) 

  coupL = (0._dp,1._dp) * coupL
  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupChi0P0


 Subroutine CoupChi0S0(i, j, k, N, RS0, gp, g, h0, h, Y_nu, lam, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RS0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u Phi
 !  h .......... trilinear coupling Phi nu^C S
 !  Y_nu ....... neutrino Yukawa coupling
 !  lam ........ self coupling of the phi field
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 10.09.2004
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: g,RS0(8,8),gp, h0, h, lam
  Complex(dp), Intent(in) :: N(10,10), Y_nu(3)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Integer :: n_S0, n_N, i1
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupChi0S0'

  n_S0 = Size(RS0, Dim=1)
  n_N = Size(N, Dim=1)

  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_S0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = 0.5_dp * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = 0.5_dp * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  !---------------------------------
  ! gauge part
  !---------------------------------
   coupL = sumI_i * ( Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2)  &
      &             + Conjg(N(j,5)) * RS0(k,3) + Conjg(N(j,6)) * RS0(k,4)  &
      &             + Conjg(N(j,7)) * RS0(k,5) )                           &
      &  + sumI_j * ( Conjg(N(i,3)) * RS0(k,1) - Conjg(N(i,4)) * RS0(k,2)  &
      &             + Conjg(N(i,5)) * RS0(k,3) + Conjg(N(i,6)) * RS0(k,4)  &
      &             + Conjg(N(i,7)) * RS0(k,5) )
  !---------------------------------
  ! Yukawa part
  !---------------------------------
  coupL = coupL + oosqrt2 * h0                                             &
      &   * ( (Conjg(N(i,3)*N(j,4)) + Conjg(N(j,3)*N(i,4))) * RS0(k,6)     &
      &     + (Conjg(N(i,3)*N(j,10)) + Conjg(N(j,3)*N(i,10))) * RS0(k,2)   &
      &     + (Conjg(N(i,4)*N(j,10)) + Conjg(N(j,4)*N(i,10))) * RS0(k,1) ) &
      & - oosqrt2 * h                                                      &
      &   * ( (Conjg(N(i,8)*N(j,9)) + Conjg(N(j,8)*N(i,9))) * RS0(k,6)     &
      &     + (Conjg(N(i,8)*N(j,10)) + Conjg(N(j,8)*N(i,10))) * RS0(k,7)   &
      &     + (Conjg(N(i,9)*N(j,10)) + Conjg(N(j,9)*N(i,10))) * RS0(k,8) )
  Do i1=1,3
   coupL = coupL - oosqrt2 * Y_nu(i1)                                        &
      &   * ( (Conjg(N(i,4)*N(j,i1+4)) + Conjg(N(j,4)*N(i,i1+4))) * RS0(k,8)   &
      &     + (Conjg(N(i,4)*N(j,8)) + Conjg(N(j,4)*N(i,8))) * RS0(k,2+i1)      &
      &     + (Conjg(N(i,i1+4)*N(j,8)) + Conjg(N(j,i1+4)*N(i,8))) * RS0(k,2) )
  End Do
  coupL = coupL - oosqrt2 * lam * Conjg(N(i,10)*N(j,10)) * RS0(k,6) 

  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine CoupChi0S0


 Subroutine CoupCSCharginoNeutralinoSpon(k, i, j, N, U, V, RSpm, Yuk, gp, g, h0 &
                                       &, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and
 ! charged scalar bosons 
 ! valid for the spontaneous model of R-parity breaking
 !   k .......... index of the charged scalar
 !  i .......... index of chargino
 !  j .......... index of the neutralino
 !  N.. ........ mixing matrix of the neutralino
 !  U,V ........ mixing matrices of the chargino
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  gp ......... U(1) gauge coupling
 !  g .......... SU(2) gauge coupling
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! the lagrangian is given by
 !  S^-(k) \bar{\chim(i)} (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 21.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k
  Real(dp), Intent(in) :: g, gp, h0
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), yuk(3)
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sumI, sumIC, ql(8), qr(8)
  Integer :: n_Spm,i1,n_neut,n_char

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCSCharginoNeutralinoSpon'

  n_char = Size(U, Dim=1)
  n_neut = Size(N, Dim=1)
  n_Spm = Size(Rspm, Dim=1)

  If ( .not.((n_char.Eq.5).And.(n_neut.Eq.10).And.(n_Spm.Eq.8))  & ! 1-gen.,spon
     &     )  Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'model not consistent defined (n_char,n_neut,n_Spm)=' &
              & ,n_char,n_neut,n_Spm
   Call TerminateProgram
  End If

  If ((i.Lt.1).Or.(i.Gt.n_char)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Chargino index out of range (i,n_char) = ',i,n_char
   Call TerminateProgram
  Else If ((j.Lt.1).Or.(j.Gt.n_neut)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Neutralino index out of range (i,n_neut) = ',j,n_neut
   Call TerminateProgram
  Else If ((k.Lt.1).Or.(k.Gt.n_Spm)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Scalar index out of range (i,n_Spm) = ',k,n_Spm
   Call TerminateProgram
  End If

  coupL = ZeroC
  coupR = ZeroC

  sumI = gp * N(j,1) + g * N(j,2)
  sumIC = Conjg(sumI)

   ql(1) = Yuk(1) * Conjg( V(i,3) * N(j,5) )   &
       & + Yuk(2) * Conjg( V(i,4) * N(j,6) )   &
       & + Yuk(3) * Conjg( V(i,5) * N(j,7) )   &
       & - h0 * Conjg( V(i,2) * N(j,10) )
   ql(2) = - (g * Conjg( V(i,1) * N(j,4)) + Conjg( V(i,2) ) * sumIC * oosqrt2)
   ql(3) = - Yuk(1) * Conjg( V(i,3) * N(j,3) )
   ql(4) = - Yuk(2) * Conjg( V(i,4) * N(j,3) )
   ql(5) = - Yuk(3) * Conjg( V(i,5) * N(j,3) )
   ql(6) = - sqrt2 * gp * Conjg( V(i,3) * N(j,1) )
   ql(7) = - sqrt2 * gp * Conjg( V(i,4) * N(j,1) )
   ql(8) = - sqrt2 * gp * Conjg( V(i,5) * N(j,1) ) 
   qr(1) = - g * U(i,1) * N(j,3) + U(i,2) * sumI * oosqrt2 
   qr(2) = - h0 * U(i,2) * N(j,10)
   qr(3) = - g * U(i,1) * N(j,5) + U(i,3) * sumI * oosqrt2 
   qr(4) = - g * U(i,1) * N(j,6) + U(i,4) * sumI * oosqrt2 
   qr(5) = - g * U(i,1) * N(j,7) + U(i,5) * sumI * oosqrt2 
   qr(6) = Conjg(Yuk(1) ) * ( U(i,2) * N(j,5) - U(i,3) * N(j,3) )
   qr(7) = Conjg(Yuk(2) ) * ( U(i,2) * N(j,6) - U(i,4) * N(j,3) )
   qr(8) = Conjg(Yuk(3) ) * ( U(i,2) * N(j,7) - U(i,5) * N(j,3) )

   Do i1=1,8
    coupL = coupL + RSpm(k,i1) * ql(i1)
    coupR = coupR + RSpm(k,i1) * qr(i1)
   End Do

 Iname = Iname - 1

 End Subroutine CoupCSCharginoNeutralinoSpon


 Subroutine CoupCharginoPseudoScalarSpon1(i, j, k, U, V, RP0, Yuk, g, h0, hnu &
                                        & , coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for 3-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g .......... SU(2) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u phi
 !  hnu(i) ..... neutrino Yukawa couplings
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 10.02.2006
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, h0, RP0(8,8)
  Complex(dp), Intent(in) :: U(5,5), V(5,5), yuk(3), hnu(3)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Complex(dp) :: yukC(3), hnuC(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoPseudoScalarSpon1'


  yukC = Conjg( yuk ) 
  hnuC = Conjg( hnu ) 

  coupL = ( RP0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )         &
        &              - yuk(1) * Conjg( V(i,3) * U(j,3) )    &
        &              - yuk(2) * Conjg( V(i,4) * U(j,4) )    &
        &              - yuk(3) * Conjg( V(i,5) * U(j,5) ) )  &
        & + g * RP0(k,2) * Conjg( V(i,2) * U(j,1) )           &
        & + RP0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )         &
        &              + yuk(1) * Conjg( V(i,3) * U(j,2) ) )  &
        & + RP0(k,4) * ( g * Conjg( V(i,1) * U(j,4) )         &
        &              + yuk(2) * Conjg( V(i,4) * U(j,2) ) )  &
        & + RP0(k,5) * ( g * Conjg( V(i,1) * U(j,5) )         &
        &              + yuk(3) * Conjg( V(i,5) * U(j,2) ) )  &
        & - h0 * RP0(k,6) *  Conjg( V(i,2) * U(j,2) )         &
        & + RP0(k,8) * ( hnu(1) * Conjg( V(i,2) * U(j,3) )    &
        &              + hnu(2) * Conjg( V(i,2) * U(j,4) )    &
        &              + hnu(3) * Conjg( V(i,2) * U(j,5) ) )  &
        &     ) * oosqrt2
  coupR = ( RP0(k,1) * ( g * V(j,1) * U(i,2)             &
        &              - yukC(1) * V(j,3) * U(i,3)       &
        &              - yukC(2) * V(j,4) * U(i,4)       &
        &              - yukC(3) * V(j,5) * U(i,5) )     &
        & + g * RP0(k,2) * V(j,2) * U(i,1)               &
        & + RP0(k,3) * ( g * V(j,1) * U(i,3)             &
        &              + yukC(1) *  V(j,3) * U(i,2) )    &
        & + RP0(k,4) * ( g *  V(j,1) * U(i,4)            &
        &              + yukC(2) * V(j,4) * U(i,2) )     &
        & + RP0(k,5) * ( g * V(j,1) * U(i,5)             &
        &              + yukC(3) * V(j,5) * U(i,2) )     &
        & - h0 * RP0(k,6) * V(j,2) * U(i,2)              &
        & + RP0(k,8) * ( hnuC(1) * V(j,2) * U(i,3)       &
        &              + hnuC(2) * V(j,2) * U(i,4)       &
        &              + hnuC(3) * V(j,2) * U(i,5)    )  &
        &     ) * oosqrt2

  coupL = coupL * (0._dp,1._dp)
  coupR = coupR * (0._dp,-1._dp) 

  Iname = Iname - 1

 End Subroutine CoupCharginoPseudoScalarSpon1

 Subroutine CoupCharginoScalarSpon1(i, j, k, U, V, RS0, Yuk, g, h0, hnu &
                                 & ,coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for 3-generation epsilon model
 !  i,j ........ index of chargino
 !  k .......... index of the scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  Yuk(i) ..... lepton yukawa couplings: i=1 -> e, i=2 -> mu, i=3 -> tau
 !  g .......... SU(2) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u phi
 !  hnu(i) ..... neutrino Yukawa couplings
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 10.02.2006
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g,RS0(8,8), h0
  Complex(dp), Intent(in) :: U(5,5), V(5,5), yuk(3), hnu(3)
  Complex(dp), Intent(out) :: coupL,coupR
  Integer, Intent(in) :: i,j,k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CoupCharginoScalarSpon1'

 
  coupL = ( RS0(k,1) * ( g * Conjg( V(i,1) * U(j,2) )             &
     &                 + yuk(1) * Conjg( V(i,3) * U(j,3) )        &
     &                 + yuk(2) * Conjg( V(i,4) * U(j,4) )        &
     &                 + yuk(3) * Conjg( V(i,5) * U(j,5) ) )      &
     &     + g * RS0(k,2) * Conjg( V(i,2) * U(j,1) )              &
     &     + RS0(k,3) * ( g * Conjg( V(i,1) * U(j,3) )            &
     &                  - yuk(1) * Conjg( V(i,3) * U(j,2) ) )     &
     &     + RS0(k,4) * ( g * Conjg( V(i,1) * U(j,4) )            &
     &                 - yuk(2) * Conjg( V(i,4) * U(j,2) ) )      &
     &     + RS0(k,5) * ( g * Conjg( V(i,1) * U(j,5) )            &
     &                  - yuk(3) * Conjg( V(i,5) * U(j,2) ) )     &
     &     - h0 * RS0(k,6) *  Conjg( V(i,2) * U(j,2) )            &
     &     + RS0(k,8) * ( hnu(1) * Conjg( V(i,2) * U(j,3) )       &
     &                  + hnu(2) * Conjg( V(i,2) * U(j,4) )       &
     &                  + hnu(3) * Conjg( V(i,2) * U(j,5) ) )     &
     &     ) * oosqrt2
  coupR = ( RS0(k,1) * ( g * Conjg( V(j,1) * U(i,2) )             &
     &                 + yuk(1) * Conjg( V(j,3) * U(i,3) )        &
     &                 + yuk(2) * Conjg( V(j,4) * U(i,4) )        &
     &                 + yuk(3) * Conjg( V(j,5) * U(i,5) ) )      &
     &    + g * RS0(k,2) * Conjg( V(j,2) * U(i,1) )               &
     &    + RS0(k,3) * ( g * Conjg( V(j,1) * U(i,3) )             &
     &                 - yuk(1) * Conjg( V(j,3) * U(i,2) ) )      &
     &    + RS0(k,4) * ( g * Conjg( V(j,1) * U(i,4) )             &
     &                 - yuk(2) * Conjg( V(j,4) * U(i,2) ) )      &
     &    + RS0(k,5) * ( g * Conjg( V(j,1) * U(i,5) )             &
     &                 - yuk(3) * Conjg( V(j,5) * U(i,2) ) )      &
     &    - h0 * RS0(k,6) *  Conjg( V(j,2) * U(i,2) )             &
     &    + RS0(k,8) * ( hnu(1) * Conjg( V(j,2) * U(i,3) )        &
     &                 + hnu(2) * Conjg( V(j,2) * U(i,4) )        &
     &                 + hnu(3) * Conjg( V(j,2) * U(i,5) ) )      &
     &    ) * oosqrt2

  coupL = - coupL
  coupR = - Conjg( coupR )

  Iname = Iname - 1

 End Subroutine CoupCharginoScalarSpon1


 Subroutine CoupPseudoScalarScalar3Spon1G(i, j, k, gp, g, hn, h, h0, lmbd  &
       &  , vevSM, vL, vR, vS, vP, Ahn, Ah, Ah0, Almbd, muhat, Mphi, MR      &
       &  , RS0, RP0, coup)
  !-----------------------------------------------------------------------
  ! i, j  ........... pseudoscalar indices
  ! k ............... scalar index
  !-----------------------------------------------------------------------
  Implicit None
   Integer, Intent(in) :: i, j, k
   Real(dp), Intent(in) :: gp, g, hn(3), h, h0, lmbd, vevSM(2), vL(3), vR &
         & , vS, vP, Ahn(3), Ah, Ah0, Almbd, Mphi, MR, RS0(8,8), RP0(8,8)
   Complex(dp), Intent(in) :: muhat
   Real(dp), Intent(out) :: coup

   Real(dp) :: g2, gp2, C_H0A0A0(8,8,8), vd, vu, Dot_vL_hn
   Integer :: i1,i2,i3

   C_H0A0A0 = 0._dp

   g2=g*g
   gp2=gp*gp
   vd = vevSM(1)
   vu = vevSM(2)
   Dot_vL_hn = Dot_product(vL, hn)

   C_H0A0A0(1,1,1)=-(g2*vd)/4._dp-(gp2*vd)/4._dp

   C_H0A0A0(1,2,2)=(g2*vd)/4._dp+(gp2*vd)/4._dp-h0**2*vd

   C_H0A0A0(1,2,6)=(h0*lmbd*vp)/2._dp-Ah0/sqrt2+(h0*Mphi) /sqrt2
   C_H0A0A0(1,2,7)=(h*h0*vR)/2._dp
   C_H0A0A0(1,2,8)=(h*h0*vS)/2._dp

   C_H0A0A0(1,3,3)=-(g2*vd)/4._dp-(gp2*vd)/4._dp

   C_H0A0A0(1,3,6)=(h0*vR*hn(1))/2._dp

   C_H0A0A0(1,3,8)=-(h0*vp*hn(1))/2._dp-(muhat*hn(1)) /sqrt2

   C_H0A0A0(1,4,4)=-(g2*vd)/4._dp-(gp2*vd)/4._dp

   C_H0A0A0(1,4,6)=(h0*vR*hn(2))/2._dp

   C_H0A0A0(1,4,8)=-(h0*vp*hn(2))/2._dp-(muhat*hn(2)) /sqrt2

   C_H0A0A0(1,5,5)=-(g2*vd)/4._dp-(gp2*vd)/4._dp
   C_H0A0A0(1,5,6)=(h0*vR*hn(3))/2._dp
   C_H0A0A0(1,5,8)=-(h0*vp*hn(3))/2._dp-(muhat*hn(3)) /sqrt2

   C_H0A0A0(1,6,6)=-(h0**2*vd)-(h0*lmbd*vu)/2._dp
   C_H0A0A0(1,6,8)= 0.5_dp * h0 * Dot_vL_hn 

   C_H0A0A0(1,7,8)=-(h*h0*vu)/2._dp

   C_H0A0A0(2,1,1)=(g2*vu)/4._dp+(gp2*vu)/4._dp-h0**2*vu
   C_H0A0A0(2,1,6)=(h0*lmbd*vp)/2._dp-Ah0/sqrt2+(h0*Mphi) /sqrt2
   C_H0A0A0(2,1,7)=(h*h0*vR)/2._dp
   C_H0A0A0(2,1,8)=(h*h0*vS)/2._dp

   C_H0A0A0(2,2,2)=-(g2*vu)/4._dp-(gp2*vu)/4._dp

   C_H0A0A0(2,3,3)=(g2*vu)/4._dp+(gp2*vu)/4._dp-vu*hn(1)**2
   C_H0A0A0(2,3,4)=-(vu*hn(1)*hn(2))
   C_H0A0A0(2,3,5)=-(vu*hn(1)*hn(3))
   C_H0A0A0(2,3,6)=-(h*vS*hn(1))/2._dp
   C_H0A0A0(2,3,7)=-(h*vp*hn(1))/2._dp-(MR*hn(1))/sqrt2
   C_H0A0A0(2,3,8)=Ahn(1)/sqrt2

   C_H0A0A0(2,4,4)=(g2*vu)/4._dp+(gp2*vu)/4._dp-vu*hn(2)**2
   C_H0A0A0(2,4,5)=-(vu*hn(2)*hn(3))
   C_H0A0A0(2,4,6)=-(h*vS*hn(2))/2._dp
   C_H0A0A0(2,4,7)=-(h*vp*hn(2))/2._dp-(MR*hn(2))/sqrt2
   C_H0A0A0(2,4,8)=Ahn(2)/sqrt2

   C_H0A0A0(2,5,5)=(g2*vu)/4._dp+(gp2*vu)/4._dp-vu*hn(3)**2
   C_H0A0A0(2,5,6)=-(h*vS*hn(3))/2._dp
   C_H0A0A0(2,5,7)=-(h*vp*hn(3))/2._dp-(MR*hn(3))/sqrt2
   C_H0A0A0(2,5,8)=Ahn(3)/sqrt2

   C_H0A0A0(2,6,6)=-(h0*lmbd*vd)/2._dp-h0**2*vu
   C_H0A0A0(2,6,7)= 0.5_dp * h * Dot_vL_hn 

   C_H0A0A0(2,7,8)=-(h*h0*vd)/2._dp

   C_H0A0A0(2,8,8)=-(vu*hn(1)**2)-vu*hn(2)**2-vu*hn(3)**2

   C_H0A0A0(3,1,1)=-((g2+gp2)*vL(1))/4._dp
   C_H0A0A0(3,1,6)=-(h0*vR*hn(1))/2._dp
   C_H0A0A0(3,1,8)=(h0*vp*hn(1))/2._dp+(muhat*hn(1))/sqrt2

   C_H0A0A0(3,2,2)= ((gp2+g2)*vL(1))/4._dp - hn(1)* Dot_vL_hn
   C_H0A0A0(3,2,6)=-(h*vS*hn(1))/2._dp
   C_H0A0A0(3,2,7)=-(h*vp*hn(1))/2._dp-(MR*hn(1))/sqrt2
   C_H0A0A0(3,2,8)=Ahn(1)/sqrt2

   C_H0A0A0(3,3,3)=-((g2+gp2)*vL(1))/4._dp

   C_H0A0A0(3,4,4)=-((g2+gp2)*vL(1))/4._dp

   C_H0A0A0(3,5,5)=-((gp2+g2)*vL(1))/4._dp-(gp2*vL(1))/4._dp

   C_H0A0A0(3,6,7)=(h*vu*hn(1))/2._dp
   C_H0A0A0(3,6,8)=(h0*vd*hn(1))/2._dp

   C_H0A0A0(3,8,8)=-hn(1) * Dot_vL_hn

   C_H0A0A0(4,1,1)=-((gp2+g2)*vL(2))/4._dp
   C_H0A0A0(4,1,6)=-(h0*vR*hn(2))/2._dp
   C_H0A0A0(4,1,8)=(h0*vp*hn(2))/2._dp+(muhat*hn(2))/sqrt2

   C_H0A0A0(4,2,2)= ((g2+gp2)*vL(2))/4._dp - hn(2) *Dot_vL_hn
   C_H0A0A0(4,2,6)=-(h*vS*hn(2))/2._dp
   C_H0A0A0(4,2,7)=-(h*vp*hn(2))/2._dp-(MR*hn(2))/sqrt2
   C_H0A0A0(4,2,8)=Ahn(2)/sqrt2

   C_H0A0A0(4,3,3)=-((gp2+g2)*vL(2))/4._dp

   C_H0A0A0(4,4,4)=-((gp2+g2)*vL(2))/4._dp

   C_H0A0A0(4,5,5)=-((gp2+g2)*vL(2))/4._dp

   C_H0A0A0(4,6,7)=(h*vu*hn(2))/2._dp
   C_H0A0A0(4,6,8)=(h0*vd*hn(2))/2._dp

   C_H0A0A0(4,8,8)=-hn(2)*Dot_vL_hn

   C_H0A0A0(5,1,1)=-((g2+gp2)*vL(3))/4._dp
   C_H0A0A0(5,1,6)=-(h0*vR*hn(3))/2._dp
   C_H0A0A0(5,1,8)=(h0*vp*hn(3))/2._dp+(muhat*hn(3))/sqrt2

   C_H0A0A0(5,2,2)=((g2+gp2)*vL(3))/4._dp-hn(3)*Dot_vL_hn 
   C_H0A0A0(5,2,6)=-(h*vS*hn(3))/2._dp
   C_H0A0A0(5,2,7)=-(h*vp*hn(3))/2._dp-(MR*hn(3))/sqrt2
   C_H0A0A0(5,2,8)=Ahn(3)/sqrt2

   C_H0A0A0(5,3,3)=-(g2*vL(3))/4._dp-(gp2*vL(3))/4._dp

   C_H0A0A0(5,4,4)=-((gp2+g2)*vL(3))/4._dp

   C_H0A0A0(5,5,5)=-((gp2+g2)*vL(3))/4._dp

   C_H0A0A0(5,6,7)=(h*vu*hn(3))/2._dp
   C_H0A0A0(5,6,8)=(h0*vd*hn(3))/2._dp

   C_H0A0A0(5,8,8)=-hn(3)*Dot_vL_hn

   C_H0A0A0(6,1,1)=-(h0**2*vp)-h0*muhat*sqrt2
   C_H0A0A0(6,1,2)=-(h0*lmbd*vp)/2._dp-Ah0/sqrt2-(h0*Mphi) /sqrt2
   C_H0A0A0(6,1,3)=(h0*vR*hn(1))/2._dp
   C_H0A0A0(6,1,4)=(h0*vR*hn(2))/2._dp
   C_H0A0A0(6,1,5)=(h0*vR*hn(3))/2._dp
   C_H0A0A0(6,1,6)=(h0*lmbd*vu)/2._dp
   C_H0A0A0(6,1,8)=0.5_dp * h0 * Dot_vL_hn
 
   C_H0A0A0(6,2,2)=-(h0**2*vp)-h0*muhat*sqrt2
   C_H0A0A0(6,2,3)=(h*vS*hn(1))/2._dp
   C_H0A0A0(6,2,4)=(h*vS*hn(2))/2._dp
   C_H0A0A0(6,2,5)=(h*vS*hn(3))/2._dp
   C_H0A0A0(6,2,6)=(h0*lmbd*vd)/2._dp
   C_H0A0A0(6,2,7)=-0.5_dp * h * Dot_vL_hn 

   C_H0A0A0(6,3,7)=-(h*vu*hn(1))/2._dp
   C_H0A0A0(6,3,8)=-(h0*vd*hn(1))/2._dp

   C_H0A0A0(6,4,7)=-(h*vu*hn(2))/2._dp
   C_H0A0A0(6,4,8)=-(h0*vd*hn(2))/2._dp

   C_H0A0A0(6,5,7)=-(h*vu*hn(3))/2._dp
   C_H0A0A0(6,5,8)=-(h0*vd*hn(3))/2._dp

   C_H0A0A0(6,6,6)=-(lmbd**2*vp)/2._dp+Almbd/sqrt2-(lmbd*Mphi)/sqrt2
   C_H0A0A0(6,6,7)=-(h*lmbd*vR)/2._dp
   C_H0A0A0(6,6,8)=-(h*lmbd*vS)/2._dp

   C_H0A0A0(6,7,7)=-(h**2*vp)-h*MR*sqrt2
   C_H0A0A0(6,7,8)=(h*lmbd*vp)/2._dp+Ah/sqrt2+(h*Mphi)/sqrt2

   C_H0A0A0(6,8,8)=-(h**2*vp)-h*MR*sqrt2

   C_H0A0A0(7,1,2)=-(h*h0*vR)/2._dp
   C_H0A0A0(7,1,8)=(h*h0*vu)/2._dp

   C_H0A0A0(7,2,3)=(h*vp*hn(1))/2._dp+(MR*hn(1))/sqrt2
   C_H0A0A0(7,2,4)=(h*vp*hn(2))/2._dp+(MR*hn(2))/sqrt2
   C_H0A0A0(7,2,5)=(h*vp*hn(3))/2._dp+(MR*hn(3))/sqrt2
   C_H0A0A0(7,2,6)=-0.5_dp * h * Dot_vL_hn
   C_H0A0A0(7,2,8)=(h*h0*vd)/2._dp

   C_H0A0A0(7,3,6)=-(h*vu*hn(1))/2._dp

   C_H0A0A0(7,4,6)=-(h*vu*hn(2))/2._dp

   C_H0A0A0(7,5,6)=-(h*vu*hn(3))/2._dp

   C_H0A0A0(7,6,6)=(h*lmbd*vR)/2._dp-h**2*vS
   C_H0A0A0(7,6,8)=-(h*lmbd*vp)/2._dp+Ah/sqrt2-(h*Mphi) /sqrt2

   C_H0A0A0(7,8,8)=-(h**2*vS)

   C_H0A0A0(8,1,2)=-(h*h0*vS)/2._dp
   C_H0A0A0(8,1,3)=(h0*vp*hn(1))/2._dp+(muhat*hn(1))/sqrt2
   C_H0A0A0(8,1,4)=(h0*vp*hn(2))/2._dp+(muhat*hn(2))/sqrt2
   C_H0A0A0(8,1,5)=(h0*vp*hn(3))/2._dp+(muhat*hn(3))/sqrt2
   C_H0A0A0(8,1,6)=-0.5_dp * h0 * Dot_vL_hn
   C_H0A0A0(8,1,7)=(h*h0*vu)/2._dp

   C_H0A0A0(8,2,2)=-(vR*hn(1)**2)-vR*hn(2)**2-vR*hn(3)**2
   C_H0A0A0(8,2,3)=Ahn(1)/sqrt2
   C_H0A0A0(8,2,4)=Ahn(2)/sqrt2
   C_H0A0A0(8,2,5)=Ahn(3)/sqrt2
   C_H0A0A0(8,2,7)=(h*h0*vd)/2._dp

   C_H0A0A0(8,3,3)=-(vR*hn(1)**2)
   C_H0A0A0(8,3,4)=-(vR*hn(1)*hn(2))
   C_H0A0A0(8,3,5)=-(vR*hn(1)*hn(3))
   C_H0A0A0(8,3,6)=(h0*vd*hn(1))/2._dp

   C_H0A0A0(8,4,4)=-(vR*hn(2)**2)
   C_H0A0A0(8,4,5)=-(vR*hn(2)*hn(3))
   C_H0A0A0(8,4,6)=(h0*vd*hn(2))/2._dp

   C_H0A0A0(8,5,5)=-(vR*hn(3)**2)
   C_H0A0A0(8,5,6)=(h0*vd*hn(3))/2._dp

   C_H0A0A0(8,6,6)=-(h**2*vR)+(h*lmbd*vS)/2._dp
   C_H0A0A0(8,6,7)=-(h*lmbd*vp)/2._dp+Ah/sqrt2-(h*Mphi) /sqrt2

   C_H0A0A0(8,7,7)=-(h**2*vR)

   Do i1=1,8
    Do i2=1,8
     Do i3=i2+1,8
      C_H0A0A0(i1,i3,i2) = C_H0A0A0(i1,i2,i3)
     End Do
    End Do
   End Do

   coup = 0._dp
   Do i1=1,8
    Do i2=1,8
     Do i3=1,8
      coup = coup + C_H0A0A0(i1,i2,i3) * RS0(k,i1) * RP0(j,i2) * RP0(i,i3) 
     End Do
    End Do
   End Do

  End Subroutine CoupPseudoScalarScalar3Spon1G

 Subroutine CoupScalar3Spon1Gen(i, j, k, gp, g, h0, h, hn, lmbd, Ahn, Ah, Ah0 &
                   & , Almbd, Mphi, MR, muhat, vevSM, vL, vR, vS, vP, RS0, coup)
 !-------------------------------------------------------------------------
 ! coupling between three neutral scalars in model of spontaneous R-parity
 ! breaking with one generations of singlets
 ! input:
 !  i, j, k ......... indices of the scalar particles
 !  gp .............. U(1) gauge coupling
 !  g ............... SU(2) gauge coupling
 !  h0 .............. superpotential coupling H_u H_d phi
 !  h ............... superpotential coupling nu^C S phi
 !  hn(i) ........... superpotential coupling H_u L nu^C
 !  lmbd ............superpotential coupling phi^3
 !  vevSM(i) ........ doublet Higgs vevs (v_d, v_u)
 !  vL(i) ........... sneutrino_L vevs
 !  vR .............. sneutrino_R vev
 !  vS .............. vev of the ~S field
 !  vP .............. vev of the phi field
 !  RS0(i,j) ........ mixing matrix of the neutral scalars
 ! output
 !  coup ............ coupling of S0_i  S0_j  S0_k
 !-------------------------------------------------------------------------
 Implicit none

   Integer, Intent(in) :: i, j, k
   Real(dp), Intent(in) :: gp, g, hn(3), h, h0, lmbd, vevSM(2), vL(3), vR &
         & , vS, vP, Ahn(3), Ah, Ah0, Almbd, Mphi, MR, RS0(8,8)
   Complex(dp), Intent(in) :: muhat
   Real(dp), Intent(out) :: coup

   Real(dp) :: g2, gp2, C_H0H0H0(8,8,8), Dot_vL_hn, h02, h2, hn2(3)
   Integer :: i1,i2,i3

   g2 = g**2
   gp2 = gp**2
   h02 = h0**2
   h2 = h**2
   hn2 = hn**2
   Dot_vL_hn = Dot_product(vL, hn)

   c_H0H0H0 = 0._dp

   c_H0H0H0(1,1,1) = (-3*(g2+gp2)*vevSM(1))/4._dp
   c_H0H0H0(1,1,2) = ((g2+gp2-4*h02)*vevSM(2))/4._dp
   c_H0H0H0(1,1,3) = -((g2+gp2)*vL(1))/4._dp
   c_H0H0H0(1,1,4) = -((g2+gp2)*vL(2))/4._dp
   c_H0H0H0(1,1,5) = -((g2+gp2)*vL(3))/4._dp
   c_H0H0H0(1,1,6) = -(h0*muhat*sqrt2)-h02*vp

   c_H0H0H0(1,2,2) = ((g2+gp2-4*h02)*vevSM(1))/4._dp
   c_H0H0H0(1,2,6) = Ah0*ooSqrt2+h0*Mphi*ooSqrt2+(h0*lmbd*vp)/2._dp
   c_H0H0H0(1,2,7) = (h*h0*vR)/2._dp
   c_H0H0H0(1,2,8) = (h*h0*vS)/2._dp

   c_H0H0H0(1,3,3) = -((g2+gp2)*vevSM(1))/4._dp
   c_H0H0H0(1,3,6) = (h0*vR*hn(1))/2._dp
   c_H0H0H0(1,3,8) = ((2*muhat*ooSqrt2+h0*vp)*hn(1))/2._dp

   c_H0H0H0(1,4,4) = -((g2+gp2)*vevSM(1))/4._dp
   c_H0H0H0(1,4,6) = (h0*vR*hn(2))/2._dp
   c_H0H0H0(1,4,8) = ((2*muhat*ooSqrt2+h0*vp)*hn(2))/2._dp

   c_H0H0H0(1,5,5) = -((g2+gp2)*vevSM(1))/4._dp
   c_H0H0H0(1,5,6) = (h0*vR*hn(3))/2._dp
   c_H0H0H0(1,5,8) = ((2*muhat*ooSqrt2+h0*vp)*hn(3))/2._dp

   c_H0H0H0(1,6,6) = -(h02*vevSM(1))+(h0*lmbd*vevSM(2))/2._dp
   c_H0H0H0(1,6,8) = (Dot_vL_hn*h0)/2._dp

   c_H0H0H0(1,7,8) = (h*h0*vevSM(2))/2._dp

   c_H0H0H0(2,2,2) = (-3*(g2+gp2)*vevSM(2))/4._dp
   c_H0H0H0(2,2,3) = (-4*Dot_vL_hn*hn(1)+(g2+gp2)*vL(1))/4._dp
   c_H0H0H0(2,2,4) = (-4*Dot_vL_hn*hn(2)+(g2+gp2)*vL(2))/4._dp
   c_H0H0H0(2,2,5) = (-4*Dot_vL_hn*hn(3)+(g2+gp2)*vL(3))/4._dp
   c_H0H0H0(2,2,6) = -(h0*muhat*sqrt2)-h02*vp
   c_H0H0H0(2,2,8) = -(vR*(hn2(1)+hn2(2)+hn2(3)))

   c_H0H0H0(2,3,3) = (vevSM(2)*(g2+gp2-4*hn2(1)))/4._dp
   c_H0H0H0(2,3,4) = -(vevSM(2)*hn(1)*hn(2))
   c_H0H0H0(2,3,5) = -(vevSM(2)*hn(1)*hn(3))
   c_H0H0H0(2,3,6) = -(h*vS*hn(1))/2._dp
   c_H0H0H0(2,3,7) = -((2*MR*ooSqrt2+h*vp)*hn(1))/2._dp
   c_H0H0H0(2,3,8) = -(ooSqrt2*Ahn(1))

   c_H0H0H0(2,4,4) = (vevSM(2)*(g2+gp2-4*hn2(2)))/4._dp
   c_H0H0H0(2,4,5) = -(vevSM(2)*hn(2)*hn(3))
   c_H0H0H0(2,4,6) = -(h*vS*hn(2))/2._dp
   c_H0H0H0(2,4,7) = -((2*MR*ooSqrt2+h*vp)*hn(2))/2._dp
   c_H0H0H0(2,4,8) = -(ooSqrt2*Ahn(2))

   c_H0H0H0(2,5,5) = (vevSM(2)*(g2+gp2-4*hn2(3)))/4._dp
   c_H0H0H0(2,5,6) = -(h*vS*hn(3))/2._dp
   c_H0H0H0(2,5,7) = -((2*MR*ooSqrt2+h*vp)*hn(3))/2._dp
   c_H0H0H0(2,5,8) = -(ooSqrt2*Ahn(3))

   c_H0H0H0(2,6,6) = (h0*lmbd*vevSM(1))/2._dp-h02*vevSM(2)
   c_H0H0H0(2,6,7) = -(Dot_vL_hn*h)/2._dp

   c_H0H0H0(2,7,8) = (h*h0*vevSM(1))/2._dp

   c_H0H0H0(2,8,8) = -(vevSM(2)*(hn2(1)+hn2(2)+hn2(3)))

   c_H0H0H0(3,3,3) = (-3*(g2+gp2)*vL(1))/4._dp
   c_H0H0H0(3,3,4) = -((g2+gp2)*vL(2))/4._dp
   c_H0H0H0(3,3,5) = -((g2+gp2)*vL(3))/4._dp
   c_H0H0H0(3,3,8) = -(vR*hn2(1))

   c_H0H0H0(3,4,4) = -((g2+gp2)*vL(1))/4._dp
   c_H0H0H0(3,4,8) = -(vR*hn(1)*hn(2))

   c_H0H0H0(3,5,5) = -((g2+gp2)*vL(1))/4._dp
   c_H0H0H0(3,5,8) = -(vR*hn(1)*hn(3))

   c_H0H0H0(3,6,7) = -(h*vevSM(2)*hn(1))/2._dp
   c_H0H0H0(3,6,8) = (h0*vevSM(1)*hn(1))/2._dp

   c_H0H0H0(3,8,8) = -(Dot_vL_hn*hn(1))

   c_H0H0H0(4,4,4) = (-3*(g2+gp2)*vL(2))/4._dp
   c_H0H0H0(4,4,5) = -((g2+gp2)*vL(3))/4._dp
   c_H0H0H0(4,4,8) = -(vR*hn2(2))

   c_H0H0H0(4,5,5) = -((g2+gp2)*vL(2))/4._dp
   c_H0H0H0(4,5,8) = -(vR*hn(2)*hn(3))

   c_H0H0H0(4,6,7) = -(h*vevSM(2)*hn(2))/2._dp
   c_H0H0H0(4,6,8) = (h0*vevSM(1)*hn(2))/2._dp

   c_H0H0H0(4,8,8) = -(Dot_vL_hn*hn(2))

   c_H0H0H0(5,5,5) = (-3*(g2+gp2)*vL(3))/4._dp
   c_H0H0H0(5,5,8) = -(vR*hn2(3))

   c_H0H0H0(5,6,7) = -(h*vevSM(2)*hn(3))/2._dp
   c_H0H0H0(5,6,8) = (h0*vevSM(1)*hn(3))/2._dp

   c_H0H0H0(5,8,8) = -(Dot_vL_hn*hn(3))

   c_H0H0H0(6,6,6) = ((-2*Almbd-6*lmbd*Mphi)*ooSqrt2-3*lmbd**2*vp)/2._dp
   c_H0H0H0(6,6,7) = (-(h*lmbd*vR)-2*h2*vS)/2._dp
   c_H0H0H0(6,6,8) = (-2*h2*vR-h*lmbd*vS)/2._dp

   c_H0H0H0(6,7,7) = -(h*MR*sqrt2)-h2*vp
   c_H0H0H0(6,7,8) = (-2*Ah*ooSqrt2-2*h*Mphi*ooSqrt2-h*lmbd*vp)/2._dp

   c_H0H0H0(6,8,8) = -(h*MR*sqrt2)-h2*vp

   c_H0H0H0(7,7,8) = -(h2*vR)
   c_H0H0H0(7,8,8) = -(h2*vS)

   Do i1=1,8
    Do i2=i1+1,8
     Do i3=i2+1,8
      c_H0H0H0(i1,i3,i2) = c_H0H0H0(i1,i2,i3)
      c_H0H0H0(i2,i1,i3) = c_H0H0H0(i1,i2,i3)
      c_H0H0H0(i2,i3,i1) = c_H0H0H0(i1,i2,i3)
      c_H0H0H0(i3,i1,i2) = c_H0H0H0(i1,i2,i3)
      c_H0H0H0(i3,i2,i1) = c_H0H0H0(i1,i2,i3)
    End Do
    End Do
   End Do

   coup = 0._dp
   Do i1=1,8
    Do i2=1,8
     Do i3=1,8
      coup = coup + C_H0H0H0(i1,i2,i3) * RS0(k,i1) * RS0(j,i2) * RS0(i,i3) 
     End Do
    End Do
   End Do

 End Subroutine CoupScalar3Spon1Gen

 Subroutine AllCouplingsNMSSM(g, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R            &
    & , Y_u, uU_L, uU_R, vevs, h0, A_h0, lam, A_lam, vP                      &
    & , RSpm, RP0, RS0, Umat, Vmat, Nmat, mu_in, phase_glu                   &
    & , RSlepton, Ae, Rsneut, RSup, Au, RSdown, Ad                           &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03       &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R       &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R    &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R   &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L           &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R           &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L         &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R           &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R           &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R           &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0         &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03 &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW    &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW        &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L    &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R             &
    & , cpl_SmpCN_L, cpl_SmpCN_R, c_GraDSd_L, c_GraDSd_R, c_GraUSu_L, c_GraUSu_R &
    & , c_GraLSl_L, c_GraLSl_R, c_GraNuSn_L, c_GraNuSn_R, GenerationMixing)
 !-----------------------------------------------------------------
 ! Routine for calculating all couplings of the MSSM
 ! output: 
 ! Couplings, the start generically with cpl_ and the remaining letters
 ! indicate the particles involved: C....chargino
 !                                  D....d-quark
 !                                  G....Gluino
 !                                  Gl...Gluon
 !                                  L....charged leptons
 !                                  N....neutralino
 !                                  Nu...neutrino
 !                                  P0...P^0
 !                                  Sl...slepton
 !                                  Smp..S^-
 !                                  Sn...sneutrino
 !                                  Sd...d-squark
 !                                  Su...u-squark
 !                                  S0...S^0
 !                                  W....W-boson
 !                                  U....u-quark
 !                                  Z....Z-boson
 ! In addition, an _L or _R is added in the case of fermions indicated
 ! if it is the left-handed or right-handed coupling. 
 ! written by Werner Porod, 
 ! 01.03.2001: taking InitializeDecaysMSSM as basis
 !          - charged Higgs -sfermion-sfermion, 1 Gen., MSSM
 ! 12.09.2002: changing interface to avoid global variables
 !-------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g(3)         ! gauge couplings [U(1), SU(2), SU(3)]
  Complex(dp), Intent(in) :: Y_l(3,3)  ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: uL_L(3,3) ! mixing matrix of left leptons
  Complex(dp), Intent(in) :: uL_R(3,3) ! mixing matrix of right leptons
  Complex(dp), Intent(in) :: Y_d(3,3)  ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3) ! mixing matrix of left d-quarks
  Complex(dp), Intent(in) :: uD_R(3,3) ! mixing matrix of right d-quarks
  Complex(dp), Intent(in) :: Y_u(3,3)  ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3) ! mixing matrix of left u-quarks
  Complex(dp), Intent(in) :: uU_R(3,3) ! mixing matrix of right u-quarks
  Real(dp), Intent(in) :: vevS(2)     ! MSSM Higgs vevs [v_d, v_u]
  Complex(dp), Intent(in) :: h0  ! trilinear superpotential coupling H_u H_d S
  Complex(dp), Intent(in) :: A_h0 !trilinear soft-potentialtrilinear  H_u H_d S
  Complex(dp), Intent(in) :: lam ! trilinear superpotential coupling S^3
  Complex(dp), Intent(in) :: A_lam !trilinear soft-potentialtrilinear
  Real(dp), Intent(in) :: vP  ! vev of the singlet S
  Complex(dp), Intent(in) :: RSpm(2,2) ! mixing matrix of charged scalars
  Real(dp), Intent(in) :: RP0(3,3)     ! mixing matrix of neutral pseudoscalars
  Real(dp), Intent(in) :: RS0(3,3)     ! mixing matrix of neutral scalars
  Complex(dp), Intent(in) :: Umat(2,2), Vmat(2,2) ! chargino mixing matrices
  Complex(dp), Intent(in) :: Nmat(5,5) ! neutralino mixing matrix
  Complex(dp), Intent(in) ::  mu_in   ! superpotential bilinear mu
  Complex(dp), Intent(in) ::  phase_glu ! phase of the gluino parameter M_3
  Complex(dp), Intent(in) :: RSlepton(6,6) ! slepton mixing matrix
  Complex(dp), Intent(in) :: Ae(3,3)     ! trilinear Higgs-slepton parameters
  Complex(dp), Intent(in) :: RSneut(3,3)   ! sneutrino mixing matrix
  Complex(dp), Intent(in) :: RSup(6,6)     ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Au(3,3)     ! trilinear Higgs - u-squark parameters
  Complex(dp), Intent(in) :: RSdown(6,6)   ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Ad(3,3)     ! trilinear Higgs - d-squark parameters
  Logical, Intent(in) :: GenerationMixing ! if .true. generation mixing of
                                          ! (s)fermions is taken into account

  Complex(dp), Intent(out) :: cpl_SmpSlSn(2,6,3), cpl_SmpSdSu(2,6,6)   &
      & , cpl_SmpSnSl(2,3,6), cpl_SmpSuSd(2,6,6), cpl_SmpP03(2,2,3)    &
      & , cpl_SmpP0W(2,3), cpl_SmpS03(2,2,3), cpl_SmpS0W(2,3)          &
      & , cpl_SmpLNu_L(2,3,3), cpl_SmpLNu_R(2,3,3), cpl_SmpDU_L(2,3,3) &
      & , cpl_SmpDU_R(8,3,3), cpl_SmpZ(2,2), cpl_DUW(3,3)
  Real(dp), Intent(out) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp), Intent(out) :: cpl_CCZ_L(2,2), cpl_CCZ_R(2,2)           &
      & , cpl_NNZ_L(5,5), cpl_NNZ_R(5,5), cpl_NNS0_L(5,5,3)            &
      & , cpl_NNS0_R(5,5,3), cpl_NNP0_L(5,5,3), cpl_NNP0_R(5,5,3) 
  Complex(dp), Intent(out) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,5,6), cpl_DNSd_R(3,5,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,5,6), cpl_UNSu_R(3,5,6)         &
      & , cpl_LNSl_L(3,5,6), cpl_LNSl_R(3,5,6), cpl_NuNSn_L(3,5,3)      & 
      & , cpl_NuNSn_R(3,5,3), cpl_DDP0_L(3,3,3), cpl_LLP0_L(3,3,3)      &
      & , cpl_UUP0_L(3,3,3), cpl_DDP0_R(3,3,3), cpl_LLP0_R(3,3,3)       &
      & , cpl_UUP0_R(3,3,3), cpl_DDS0_L(3,3,3), cpl_LLS0_L(3,3,3)       &
      & , cpl_UUS0_L(3,3,3), cpl_DDS0_R(3,3,3), cpl_LLS0_R(3,3,3)       &
      & , cpl_UUS0_R(3,3,3)
  Complex(dp), Intent(out) :: cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6)      &
      & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)       &
      & , cpl_CLSn_R(2,3,3), cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)
  Complex(dp), Intent(out) :: cpl_GlGlS0(3)
  Complex(dp) :: cpl_P0SdSd(3,6,6), cpl_P0SuSu(3,6,6), cpl_P0SlSl(3,6,6) &
      & , cpl_P0SnSn(3,3,3), cpl_P0S0Z(3,3) 
  Real(dp), Intent(out) :: cpl_P0S03(3,3,3)
  Complex(dp), Intent(out) :: cpl_S0SdSd(3,6,6), cpl_S0SuSu(3,6,6) &
      & , cpl_S0SlSl(3,6,6), cpl_S0SnSn(3,3,3), cpl_LNuW(3,3)
  Real(dp), Intent(out) :: cpl_S03(3,3,3), cpl_S0WW(3), cpl_S0ZZ(3), cpl_FFpW
  Complex(dp), Intent(out) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp), Intent(out) :: cpl_CCP0_L(2,2,3), cpl_CCP0_R(2,2,3)    &
      & , cpl_CCS0_L(2,2,3), cpl_CCS0_R(2,2,3), cpl_CNW_L(2,5)        &
      & , cpl_CNW_R(2,5), cpl_SmpCN_L(2,2,5), cpl_SmpCN_R(2,2,5)
  Complex(dp), Intent(out), Dimension(3,6) :: c_GraDSd_L, c_GraDSd_R, c_GraUSu_L &
      & , c_GraUSu_R, c_GraLSl_L, c_GraLSl_R
  Complex(dp), Intent(out), Dimension(3,3) :: c_GraNuSn_L, c_GraNuSn_R

  Integer :: n_char, n_neut, n_S0, n_P0, n_Spm 
  Integer :: i1, i2, i3, i4
  Real(dp) :: gU1, gSU2, gSU3, sinW2, cosW2, cosW
  Complex(dp) :: Rsd(2,2), Rsu(2,2), Rsl(2,2), coupLC, coupRC, Yuk, Yukp, &
               & coupC, A, Ap, mat6(6,6), mat3(3,3), bi(1), mu
  Real(dp), Parameter ::  e_u = 2._dp / 3._dp,  e_d = - 1._dp / 3._dp

  Iname = Iname + 1
  NameOfUnit(Iname) = 'AllCouplingsNMSSM'

  !---------------------------
  ! specifying the couplings
  !---------------------------
  gU1 = g(1)
  gSU2 = g(2)
  gSU3 = g(3)
  sinW2 = gU1**2 / (gU1**2 + gSU2**2)
  cosW2 = 1._dp - sinW2
  cosW = Sqrt(cosW2)

  !------------------------
  ! specifying the model
  !------------------------
  n_char = 2
  n_neut = 5
  n_P0 = 3
  n_S0 = 3
  n_Spm = 2

  !--------------------
  ! abbriviations
  !--------------------
  mu = mu_in + oosqrt2 * h0 * vP

  !----------------------------------------
  ! charged scalar - chargino - neutralino
  !----------------------------------------
  cpl_SmpCN_L = 0._dp
  cpl_SmpCN_R = 0._dp

  Do i1 = 1,n_Spm
   Do i2 = 1,n_char
    Do i3 = 1,n_neut
     Call CoupCSCharginoNeutralino(i1, i2, i3, Nmat, Umat, Vmat, RSpm     &
               &, gU1, gSU2, h0, cpl_SmpCN_L(i1,i2,i3), cpl_SmpCN_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !----------------------------------------
  ! charged scalar - fermion - fermion'
  !----------------------------------------
  cpl_SmpLNu_L = 0._dp
  cpl_SmpLNu_R = 0._dp
  cpl_SmpDU_L = 0._dp
  cpl_SmpDU_R = 0._dp
  mat3 = 0._dp

  If (GenerationMixing) Then
   Do i1=1,2
    Do i2=1,3
     Do i3=1,3
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_l, uL_L, uL_R, mat3 &
              &, id3C, id3C, cpl_SmpLNu_L(i1,i2,i3), cpl_SmpLNu_R(i1,i2,i3) )
      Call CoupChargedScalarFermion(i1, i2, i3, RSpm, Y_D, uD_L, uD_R, Y_U &
              &, uU_L, uU_R, cpl_SmpDU_L(i1,i2,i3), cpl_SmpDU_R(i1,i2,i3) )
     End Do
    End Do
   End Do
  Else
   Do i1=1,2
    Do i2=1,3
     Call CoupChargedScalarFermion(i1, RSpm, Y_l(i2,i2), ZeroC               &
                       &, cpl_SmpLNu_L(i1,i2,i2), cpl_SmpLNu_R(i1,i2,i2) )
     Call CoupChargedScalarFermion(i1, RSpm, Y_D(i2,i2), Y_U(i2,i2)          &
                         &, cpl_SmpDU_L(i1,i2,i2), cpl_SmpDU_R(i1,i2,i2) )
    End Do
   End Do
  End If
  !-------------------------------------
  ! charged scalar - pseudoscalar - W
  !-------------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_P0
     Call CoupChargedScalarPseudoscalarW(i1, i2, gSU2, RSpm, RP0 &
                                      &, cpl_SmpP0W(i1,i2) )
   End Do
  End Do

  !------------------------------------------------
  ! charged scalar - charged scalar - pseudoscalar
  !-----------------------------------------------
  cpl_SmpP03 = 0._dp
!  Do i1 = 1,n_Spm
!   Do i2 = 1,n_Spm
!    Do i3 = 1,n_P0
!     Call CoupChargedScalarPseudoscalar3(i1, i2, i3, RSpm, RP0, vevs, gSU2 &
!                                       &, cpl_SmpP03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do

  !------------------------------
  ! charged scalar - scalar - W
  !------------------------------
  Do i1 = 1,n_Spm
   Do i2 = 1,n_S0
     Call CoupChargedScalarScalarW(i1, i2, gSU2, RSpm, RS0, cpl_SmpS0W(i1,i2) )
   End Do
  End Do

  !--------------------------------------
  ! charged scalar - charged scalar - scalar
  !--------------------------------------
  cpl_SmpS03 = 0._dp
!   Do i1 = 1, n_Spm
!    Do i2 = 1, n_Spm
!     Do i3 = 1, n_S0
!      Call CoupChargedScalarScalar3(i1, i2, i3, RSpm, RS0, vevs, gU1, gSU2 &
!                                   &, cpl_SmpS03(i1,i2,i3) )
!    End Do
!   End Do
!  End Do
  !--------------------------------------
  ! charged scalar - sfermion - sfermion
  !--------------------------------------
  cpl_SmpSlSn = ZeroC
  cpl_SmpSnSl = ZeroC
  cpl_SmpSdSu = ZeroC
  cpl_SmpSuSd = ZeroC
  mat3 = zeroC

  If (GenerationMixing) Then
   Do i1 = 1, n_Spm
    Do i2 = 1,6
     Do i3 = 1,6
      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
                                  & , Y_d, Ad, Rsdown, Y_u, Au, Rsup, coupC)
      cpl_SmpSdSu(i1, i2, i3) = coupC
      cpl_SmpSuSd(i1, i3, i2) = Conjg(coupC)
     End Do
     Do i3 = 1,3
      Call CoupChargedScalarSfermion3(i1, i2, i3, RSpm, gSU2, vevs, mu     &
                           & , Y_l, Ae, Rslepton, mat3, mat3, Rsneut, coupC)
      cpl_SmpSlSn(i1, i2, i3) = coupC
      cpl_SmpSnSl(i1, i3, i2) = Conjg(coupC)
     End Do
    End Do
   End Do

  Else
   Do i1 =1,n_Spm
    Do i2=1,3
     !----------
     ! Sleptons
     !----------
     Yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Call CoupChargedScalarSfermion3(i1, i3, 1, RSpm, gSU2, vevs, mu,    &
                       &             Yuk, A, Rsl, ZeroC, ZeroC, Id2C, coupC)
      cpl_SmpSlSn(i1,2*(i2-1)+i3,i2) = coupC
      cpl_SmpSnSl(i1,i2,2*(i2-1)+i3) = Conjg(coupC)
     End Do
     !----------
     ! Squarks
     !----------
     Yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Yukp = Y_u(i2,i2)
     Ap = Au(i2,i2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Do i3=1,2
      Do i4=1,2
       Call CoupChargedScalarSfermion3(i1,i3,i4,RSpm,gSU2,vevs,mu,    &
                       &               yuk,A,Rsd,yukp,Ap,Rsu,coupC)
       cpl_SmpSdSu(i1,2*(i2-1)+i3,2*(i2-1)+i4) = coupC
       cpl_SmpSuSd(i1,2*(i2-1)+i4,2*(i2-1)+i3) = Conjg(coupC)
      End Do
     End Do
    End Do ! i2
   End Do ! i1
  End If

  !-------------------------------------
  ! charged scalar - Z
  !-------------------------------------
  Do i1=1,n_Spm
   Do i2=1,n_Spm
    Call CoupChargedScalarZ(i1, i2, gSU2, sinW2, RSpm, cpl_SmpZ(i1,i2) )
   End Do
  End Do

  !-------------------------------------
  ! chargino - chargino - pseudoscalar
  !-------------------------------------
  cpl_CCP0_L = 0.0_dp
  cpl_CCP0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_P0
    Call CoupCharginoPseudoScalar(i1, i2, i3, Umat, Vmat, RP0, gSU2, h0 &
                     &, cpl_CCP0_L(i1,i2,i3), cpl_CCP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !-------------------------------------
  ! chargino - chargino - scalar
  !-------------------------------------
  cpl_CCS0_L = 0.0_dp
  cpl_CCS0_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Do i3=1,n_S0
    Call CoupCharginoScalar(i1, i2, i3, Umat, Vmat, RS0, gSU2, h0  &
                     &, cpl_CCS0_L(i1,i2,i3), cpl_CCS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  !---------------------------
  ! chargino - chargino - Z
  !---------------------------
  cpl_CCZ_L = 0.0_dp
  cpl_CCZ_R = 0.0_dp
  Do i1=1,n_char
   Do i2=1,n_char
    Call CoupCharginoZ(i1, i2, Umat, Vmat, gSU2, cosW  &
                     &, cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2) )
   End Do
  End Do

  !--------------------------------
  ! chargino - fermion - sfermion
  !--------------------------------
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp
  cpl_CLSn_L = 0._dp
  cpl_CLSn_R = 0._dp
  cpl_CUSd_L = 0._dp
  cpl_CUSd_R = 0._dp
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp

  If (GenerationMixing) Then

   mat6 = 0._dp
   mat6(1:3,1:3) = Rsneut
   mat3 = 0._dp
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSlepton, Y_l &
                              &, mat3, id3c, id3c, Umat, Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i1, i2, i3) = coupLC
      cpl_CNuSl_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, 0.5_dp, RSdown, Y_d &
                              &, Y_u, uU_L, uU_R, Umat, Vmat, coupLC, coupRC)
      cpl_CUSd_L(i1, i2, i3) = coupLC
      cpl_CUSd_R(i1, i2, i3) = coupRC
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, RSup, Y_d &
                              &, Y_u, uD_L, uD_R, Umat, Vmat, coupLC, coupRC)
      cpl_CDSu_L(i1, i2, i3) = coupLC
      cpl_CDSu_R(i1, i2, i3) = coupRC
     End Do
     Do i3=1,3
      Call  CoupCharginoSfermion(i1, i2, i3, gSU2, -0.5_dp, mat6, Y_l, mat3 &
                              &, uL_L, uL_R, Umat, Vmat, coupLC, coupRC)
      cpl_CLSn_L(i1, i2, i3) = coupLC
      cpl_CLSn_R(i1, i2, i3) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    YukP = 0._dp
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSl, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CNuSl_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CNuSl_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
     Call CoupCharginoSfermion(i2, 1, gSU2, -0.5_dp, id2C, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
     cpl_CLSn_L(i2, i1, i1) = coupLC
     cpl_CLSn_R(i2, i1, i1) = coupRC
    End Do

    Yuk = Y_d(i1,i1)
    YukP = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,2
      Call CoupCharginoSfermion(i2, i3, gSU2, 0.5_dp, RSd, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CUSd_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CUSd_R(i2, i1, (i1-1)*2 + i3) = coupRC
      Call CoupCharginoSfermion(i2, i3, gSU2, -0.5_dp, RSu, Yuk, Yukp, Umat, &
                              & Vmat, coupLC, coupRC)
      cpl_CDSu_L(i2, i1, (i1-1)*2 + i3) = coupLC
      cpl_CDSu_R(i2, i1, (i1-1)*2 + i3) = coupRC
     End Do
    End Do

   End Do
  End If 

  !-------------------------
  ! chargino - neutralino W
  !-------------------------
  cpl_CNW_L = 0._dp
  cpl_CNW_R = 0._dp

  Do i1 = 1,n_char
   Do i2= 1,n_neut
    Call CoupCharginoNeutralinoW(i1, i2, Nmat, Umat, Vmat, gSU2 &
                               &, cpl_CNW_L(i1,i2), cpl_CNW_R(i1,i2) )
   End Do
  End Do

  !----------------------------------
  ! fermion - fermion - pseudoscalar
  !--------------------------------
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_LLP0_L = 0._dp
  cpl_LLP0_R = 0._dp
  cpl_UUP0_L = 0._dp
  cpl_UUP0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_P0
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RP0 &
                                 &, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, -0.5_dp, Y_l, uL_L, uL_R, RP0 &
                                 &, cpl_LLP0_L(i1,i2,i3), cpl_LLP0_R(i1,i2,i3))
      Call CoupFermionPseudoScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RP0 &
                                 &, cpl_UUP0_L(i1,i2,i3), cpl_UUP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_P0
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_d(i1,i1), RP0   &
                                &, cpl_DDP0_L(i1,i1,i3), cpl_DDP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, -0.5_dp, Y_l(i1,i1), RP0   &
                                &, cpl_LLP0_L(i1,i1,i3), cpl_LLP0_R(i1,i1,i3))
     Call CoupFermionPseudoScalar(i3, 0.5_dp, Y_u(i1,i1), RP0   &
                                &, cpl_UUP0_L(i1,i1,i3), cpl_UUP0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !----------------------------------
  ! fermion - fermion - scalar
  !--------------------------------
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_LLS0_L = 0._dp
  cpl_LLS0_R = 0._dp
  cpl_UUS0_L = 0._dp
  cpl_UUS0_R = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,3
     Do i3 = 1,n_S0
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_d, uD_L, uD_R, RS0 &
                           &, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, -0.5_dp, Y_l, uL_L, uL_R, RS0 &
                           &, cpl_LLS0_L(i1,i2,i3), cpl_LLS0_R(i1,i2,i3))
      Call CoupFermionScalar(i1, i2, i3, 0.5_dp, Y_u, uU_L, uU_R, RS0 &
                           &, cpl_UUS0_L(i1,i2,i3), cpl_UUS0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  Else
   Do i1 = 1,3
    Do i3 = 1,n_S0
     Call CoupFermionScalar(i3, -0.5_dp, Y_d(i1,i1), RS0   &
                          &, cpl_DDS0_L(i1,i1,i3), cpl_DDS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, -0.5_dp, Y_l(i1,i1), RS0   &
                          &, cpl_LLS0_L(i1,i1,i3), cpl_LLS0_R(i1,i1,i3))
     Call CoupFermionScalar(i3, 0.5_dp, Y_u(i1,i1), RS0   &
                          &, cpl_UUS0_L(i1,i1,i3), cpl_UUS0_R(i1,i1,i3))
    End Do
   End Do
  End If

  !-----------------------
  ! fermion - fermion - Z
  !-----------------------
  Call CoupFermionZ(-0.5_dp, -1._dp, gSU2, sinW2, cpl_LLZ_L, cpl_LLZ_R)
  Call CoupFermionZ(-0.5_dp, e_d, gSU2, sinW2, cpl_DDZ_L, cpl_DDZ_R)
  Call CoupFermionZ(0.5_dp, e_u, gSU2, sinW2, cpl_UUZ_L, cpl_UUZ_R)
  Call CoupFermionZ(0.5_dp, 0._dp, gSU2, sinW2, cpl_NuNuZ_L, cpl_NuNuZ_R)

  !-----------------------
  ! femion - fermion' - W
  !-----------------------
  cpl_FFpW = - gSU2 * oosqrt2
  cpl_LNuW = id3C * cpl_FFpW

  If (GenerationMixing) Then
   cpl_DUW = CKM * cpl_FFpW
  Else
   cpl_DUW = id3C * cpl_FFpW
  End If
  !--------------------------------------
  ! fermion - neutralino - sfermion
  !--------------------------------------
  cpl_LNSl_L = 0._dp
  cpl_LNSl_R = 0._dp
  cpl_NuNSn_L = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_DNSd_L = ZeroC
  cpl_DNSd_R = ZeroC
  cpl_UNSu_L = ZeroC
  cpl_UNSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1 = 1,3
    Do i2 = 1,5
     Do i3 = 1,6
      Call CoupNeutralinoSlepton(i1, i2, i3, gU1, gSU2, RSlepton, uL_L, uL_R &
                               &, Y_l, Nmat, coupLC, coupRC)
      cpl_LNSl_L(i1, i2, i3 ) = coupLC
      cpl_LNSl_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSdown(i1, i2, i3, gU1, gSU2, RSdown, uD_L, uD_R &
                               &, Y_d, Nmat, coupLC, coupRC)
      cpl_DNSd_L(i1, i2, i3 ) = coupLC
      cpl_DNSd_R(i1, i2, i3 ) = coupRC
      Call CoupNeutralinoSup(i1, i2, i3, gU1, gSU2, RSup, uU_L, uU_R &
                               &, Y_u, Nmat, coupLC, coupRC)
      cpl_UNSu_L(i1, i2, i3 ) = coupLC
      cpl_UNSu_R(i1, i2, i3 ) = coupRC
     End Do
     Do i3 = 1,3
      Call CoupNeutralinoSneutrino(i1, i2, i3, gU1, gSU2, Nmat, RSneut, id3C &
                               &, coupRC)
      cpl_NuNSn_R(i1, i2, i3 ) = coupRC
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,5
    Call CoupNeutralinoSneutrino(i1, gU1, gSU2, Nmat, coupRC)
    Do i2 = 1,3
     cpl_NuNSn_R(i2, i1, i2 ) = coupRC
    End Do
   End Do

   Do i1 = 1,3
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Yuk = Y_l(i1,i1)
    Do i2 = 1,2
     Do i3 = 1,5
      Call CoupNeutralinoSlepton(i3, i2, gU1, gSU2, RSl, Yuk, Nmat &
                               &, coupLC, coupRC)
      cpl_LNSl_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_LNSl_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_d(i1,i1)
    Do i2=1,2
     Do i3=1,5
      Call CoupNeutralinoSdown(i3, i2, gU1, gSU2, RSd, Yuk, Nmat, &
                             & coupLC, coupRC)
      cpl_DNSd_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_DNSd_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

    Yuk = Y_u(i1,i1)
    Do i2=1,2
     Do i3=1,5
      Call CoupNeutralinoSup(i3, i2, gU1, gSU2, RSu, Yuk, Nmat, &
                           & coupLC, coupRC)
      cpl_UNSu_L(i1, i3, (i1-1)*2 + i2 ) = coupLC
      cpl_UNSu_R(i1, i3, (i1-1)*2 + i2 ) = coupRC
     End Do
    End Do

   End Do ! i1

  End If

  !--------------------------------
  ! Gravitino - fermion - sfermion
  !--------------------------------
  c_GraNuSn_L = 0._dp
  c_GraNuSn_R = 0._dp
  c_GraLSl_L = 0._dp
  c_GraLSl_R = 0._dp
  c_GraDSd_L = 0._dp
  c_GraDSd_R = 0._dp
  c_GraUSu_L = 0._dp
  c_GraUSu_R = 0._dp

  If (GenerationMixing) Then

   Do i1=1,3
    Do i2=1,3
     Call CoupGravitinoSfermion(i1, i2, RSneut, RSneut, id3c &
                  & , c_GraNuSn_L(i1,i2), c_GraNuSn_R(i1,i2))
    End Do
    Do i2=1,6
     Call CoupGravitinoSfermion(i1, i2, RSlepton, uL_L, uL_R &
                  & , c_GraLSl_L(i1,i2), c_GraLSl_R(i1,i2))
     Call CoupGravitinoSfermion(i1, i2, RSdown, uD_L, uD_R &
                  & , c_GraDSd_L(i1,i2), c_GraDSd_R(i1,i2))
     Call CoupGravitinoSfermion(i1, i2, RSup, uU_L, uU_R &
                  & , c_GraUSu_L(i1,i2), c_GraUSu_R(i1,i2))
    End Do
   End Do

  Else
   
   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    c_GraNuSn_R(i1, i1) = 1._dp
    Do i2=1,2
     Call CoupGravitinoSfermion(i2, RSd, coupLC, coupRC)
     c_GraDSd_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGravitinoSfermion(i2, RSu, coupLC, coupRC)
     c_GraUSu_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraUSu_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGravitinoSfermion(i2, RSl, coupLC, coupRC)
     c_GraLSl_L(i1, (i1-1)*2 + i2) = coupLC
     c_GraLSl_R(i1, (i1-1)*2 + i2) = coupRC
    End Do
   End Do

  End If 

  !--------------------------
  ! Gluino
  !--------------------------
  cpl_GDSd_L = ZeroC
  cpl_GDSd_R = ZeroC
  cpl_GUSu_L = ZeroC
  cpl_GUSu_R = ZeroC

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,6
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsdown, uD_L, uD_R, &
                          & coupLC, coupRC)
     cpl_GDSd_L(i1, i2) = coupLC
     cpl_GDSd_R(i1, i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i1, i2, Rsup, uU_L, uU_R, &
                          & coupLC, coupRC)
     cpl_GUSu_L(i1, i2) = coupLC
     cpl_GUSu_R(i1, i2) = coupRC
    End Do
   End Do
  Else

   Do i1=1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)

    Do i2 =1,2
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsd, coupLC, coupRC)
     cpl_GDSd_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GDSd_R(i1, (i1-1)*2 + i2) = coupRC
     Call CoupGluinoSquark(gSU3, phase_glu, i2, Rsu, coupLC, coupRC)
     cpl_GUSu_L(i1, (i1-1)*2 + i2) = coupLC
     cpl_GUSu_R(i1, (i1-1)*2 + i2) = coupRC
    End Do  
   End Do

  End If 

  !-------------------
  ! Gluon Gluon scalar
  !-------------------
  cpl_GlGlS0 = 0._dp

  !-----------------------------------------
  ! neutralino - neutralino - pseudoscalar
  !-----------------------------------------
  cpl_NNP0_L = 0.0_dp
  cpl_NNP0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_P0
    Call Coup_Chi0_Chi0_P0_NMSSM(i1, i2, i3, Nmat, RP0, gU1, gSU2, h0, lam &
                       & , cpl_NNP0_L(i1,i2,i3), cpl_NNP0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !-----------------------------------------
  ! neutralino - neutralino - scalar
  !-----------------------------------------
  cpl_NNS0_L = 0.0_dp
  cpl_NNS0_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Do i3=1,n_S0
    Call Coup_Chi0_Chi0_S0_NMSSM(i1, i2, i3, Nmat, RS0, gU1, gSU2, h0, lam &
                       & , cpl_NNS0_L(i1,i2,i3), cpl_NNS0_R(i1,i2,i3) )
    End Do
   End Do
  End Do

  !------------------------------
  ! neutralino - neutralino - Z
  !------------------------------
  cpl_NNZ_L = 0.0_dp
  cpl_NNZ_R = 0.0_dp
  Do i1=1,n_neut
   Do i2=1,n_neut
    Call CoupNeutralinoZ(i1, i2, Nmat, gSU2, cosW, &
                       & cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2), .True. )
   End Do
  End Do

  !------------------------------
  ! sfermion - sfermion - W
  !------------------------------
  cpl_SdSuW = 0._dp
  cpl_SlSnW = 0._dp

  If (GenerationMixing) Then
   Do i1 = 1,6
    Do i2 = 1,6
     Call CoupSfermionW3(i1, i2, gSU2, RSdown, RSup, cpl_SdSuW(i1,i2) )
    End Do
    Do i2 = 1,3
     Call CoupSfermionW3(i1, i2, gSU2, RSlepton, RSneut, cpl_SlSnW(i1,i2) )
    End Do
   End Do

  Else
   Do i1 = 1,3
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Do i3 = 1,2
      Call CoupSfermionW3(i2, i3, gSU2, RSd, RSu, coupC )
      cpl_SdSuW( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do

    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2 = 1,2
     Call CoupSfermionW3(i2, 1, gSU2, RSl, id2C, coupC )
      cpl_SlSnW( (i1-1)*2+i2, i1) = coupC
    End Do
   End Do

  End If

  Call Adjungate(cpl_SdSuW, cpl_SuSdW)
  Call Adjungate(cpl_SlSnW, cpl_SnSlW)

  !---------------------------
  ! Pseudoscalar - scalar 
  !---------------------------
  cpl_P0S03 = 0._dp
!   Do i1=1,n_P0
!    Do i2=1,n_P0
!     Do i3=1,n_S0
!      Call CoupPseudoScalarScalar3(i1, i2, i3, RP0, RS0, gU1, gSU2, vevs &
!                                 &, cpl_P0S03(i1,i2,i3) )
!      End Do
!     End Do
!    End Do
  !---------------------------
  ! Pseudoscalar - scalar - Z
  !---------------------------
  cpl_P0S0Z = 0._dp

  Do i1=1,n_P0
   Do i2=1,n_S0
    Call CoupPseudoscalarScalarZ(i1, i2, gSU2, cosW, RP0, RS0,cpl_P0S0Z(i1,i2))
    End Do
   End Do
  !-------------------------------------
  ! Pseudoscalar - sfermion - sfermion 
  !-------------------------------------
  cpl_P0SdSd = 0._dp
  cpl_P0SuSu = 0._dp
  cpl_P0SlSl = 0._dp
  cpl_P0SnSn = 0._dp

  bi(1) = mu

!   If (GenerationMixing) Then
!    Do i1=1,n_P0
!     Do i2=1,6
!      Do i3=1,6
!       Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_d, Rsdown   &
!                                    &, Ad, bi, cpl_P0SdSd(i1,i2,i3) )
!       Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, 0.5_dp, Y_u, Rsup      &
!                                    &, Au, bi, cpl_P0SuSu(i1,i2,i3) )
!       Call CoupPseudoScalarSfermion3(i1, i2, i3, RP0, -0.5_dp, Y_l, Rslepton &
!                                    &, Ae, bi, cpl_P0SlSl(i1,i2,i3) )
!      End Do
!     End Do
!    End Do
! 
!   Else
   Do i1=1,n_P0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsl   &
                                   &, A, mu_in, h0, vevS, vP, coupC )
       cpl_P0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, -0.5_dp, Yuk, Rsd   &
                                   &, A, mu_in, h0, vevS, vP,  coupC )
       cpl_P0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupPseudoScalarSfermion3(i1, i3, i4, RP0, 0.5_dp, Yuk, Rsu   &
                                   &, A, mu_in, h0, vevS, vP, coupC )
       cpl_P0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

!   End If

  !-------------------------------------
  ! scalar - scalar -scalar
  !-------------------------------------
  cpl_S03 = 0._dp
!   Do i1=1,n_S0
!    Do i2=1,n_S0
!     Do i3=1,n_S0
!      Call CoupScalar3(i1, i2, i3, RS0, gU1, gSU2, vevs, cpl_S03(i1,i2,i3))
!     End Do
!    End Do
!   End Do

  !-------------------------------------
  ! scalar - sfermion - sfermion 
  !-------------------------------------
  cpl_S0SdSd = 0._dp
  cpl_S0SuSu = 0._dp
  cpl_S0SlSl = 0._dp
  cpl_S0SnSn = 0._dp
  mat3 = zeroC

!   If (GenerationMixing) Then
!    Do i1=1,n_S0
!     Do i2=1,6
!      Do i3=1,6
!       Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, e_d, Y_d, Rsdown   &
!                            &, Ad, mu, vevs, gU1, gSU2, cpl_S0SdSd(i1,i2,i3) )
!       Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, e_u, Y_u, Rsup      &
!                            &, Au, mu, vevs, gU1, gSU2, cpl_S0SuSu(i1,i2,i3) )
!       Call CoupScalarSfermion3(i1, i2, i3, RS0, -0.5_dp, -1._dp, Y_l        &
!                  &, Rslepton, Ae, mu, vevs, gU1, gSU2, cpl_S0SlSl(i1,i2,i3) )
!      End Do
!     End Do
!     Do i2=1,3
!      Do i3=1,3
!       Call CoupScalarSfermion3(i1, i2, i3, RS0, 0.5_dp, 0._dp, mat3, Rsneut &
!                          &, mat3, mu, vevs, gU1, gSU2, cpl_S0SnSn(i1,i2,i3) )
!      End Do
!     End Do
!    End Do
! 
!   Else
   Do i1=1,n_S0
    Do i2=1,3
     Rsl = RSlepton(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsd = RSdown(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)
     Rsu = RSup(2*(i2-1)+1:2*(i2-1)+2, 2*(i2-1)+1:2*(i2-1)+2)

     Call CoupScalarSfermion3(i1, 1, 1, RS0, 0.5_dp, 0._dp, ZeroC, id3C  &
                               &, ZeroC, mu_in, vevs, vP, gU1, gSU2, h0, coupC )
     cpl_S0SnSn(i1,i2,i2) = coupC

     yuk = Y_l(i2,i2)
     A = Ae(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, -1._dp, Yuk, Rsl  &
                               &, A, mu_in, vevs, vP, gU1, gSU2, h0, coupC )
       cpl_S0SlSl(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_d(i2,i2)
     A = Ad(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, -0.5_dp, e_d, Yuk, Rsd   &
                               &, A, mu_in, vevs, vP, gU1, gSU2, h0, coupC )
       cpl_S0SdSd(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do

     yuk = Y_u(i2,i2)
     A = Au(i2,i2)
     Do i3=1,2
      Do i4=1,2
       Call CoupScalarSfermion3(i1, i3, i4, RS0, 0.5_dp, e_u, Yuk, Rsu   &
                               &, A, mu_in, vevs, vP, gU1, gSU2, h0, coupC )
       cpl_S0SuSu(i1,(i2-1)*2+i3,(i2-1)*2+i4) = coupC
      End Do
     End Do
    End Do
   End Do

!   End If

  !----------
  ! scalar W
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarW(i1, gSU2, vevs, RS0, cpl_S0WW(i1) )
  End Do

  !----------
  ! scalar Z
  !----------
  Do i1 = 1,n_S0
   Call CoupScalarZ(i1, gSU2, cosW2, vevs, RS0, cpl_S0ZZ(i1) )
  End Do

  !-------------------------
  ! sfermion - sfermion - Z
  !-------------------------
  cpl_SdSdZ = 0._dp
  cpl_SlSlZ = 0._dp
  cpl_SnSnZ = 0._dp
  cpl_SuSuZ = 0._dp

  If (GenerationMixing) Then
   Do i1=1,3
    Do i2=1,3
     Call CoupSneutrinoZ(i1, i2, gSU2, sinW2, Rsneut, cpl_SnSnZ(i1, i2) )
    End Do
   End Do 
   Do i1=1,6
    Do i2=1,6
     Call CoupSleptonZ(i1, i2, gSU2, sinW2, Rslepton, cpl_SlSlZ(i1, i2) )
     Call CoupSdownZ(i1, i2, gSU2, sinW2, Rsdown, cpl_SdSdZ(i1, i2) )
     Call CoupSupZ(i1, i2, gSU2, sinW2, Rsup, cpl_SuSuZ(i1, i2) )
    End Do
   End Do
 
  Else
   Call CoupSneutrinoZ(gSU2, sinW2, cpl_SnSnZ(1,1) )
   Do i1=1,3
    cpl_SnSnZ(i1,i1) = cpl_SnSnZ(1,1)
    Rsd = RSdown(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsu = RSup(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Rsl = RSlepton(2*(i1-1)+1:2*(i1-1)+2, 2*(i1-1)+1:2*(i1-1)+2)
    Do i2=1,2
     Do i3=1,2
      Call CoupSleptonZ(i2, i3, gSU2, sinW2, Rsl, coupC )
      cpl_SlSlZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSdownZ(i2, i3, gSU2, sinW2, Rsd, coupC )
      cpl_SdSdZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
      Call CoupSupZ(i2, i3, gSU2, sinW2, Rsu, coupC )
      cpl_SuSuZ( (i1-1)*2+i2, (i1-1)*2+i3 ) = coupC
     End Do
    End Do    
   End Do

  End If

  Iname = Iname - 1

 End Subroutine AllCouplingsNMSSM

 Subroutine Coup_Chip_Chim_P0_NMSSM(i, j, k, U, V, RP0, g, h0, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and pseudoscalar bosons 
 ! valid for the NMSSM
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RP0(i,j) ... mixing matrix of pseudoscalar bosons
 !  g .......... SU(2) gauge coupling
 !  h0 ......... trilinear coupling h_u H_d phi
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 03.03.2005
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RP0(3,3)
  Complex(dp), Intent(in) :: h0, U(2,2), V(2,2)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_Chip_Chim_P0_NMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = ( g * ( RP0(k,1) * Conjg( V(i,1) * U(j,2) )             &
        &       + RP0(k,2) * Conjg( V(i,2) * U(j,1) ) )           &
        & - h0 * RP0(k,3) * Conjg( V(i,2) * U(j,2) )    ) * oosqrt2
  coupR = ( g * ( RP0(k,1) * V(j,1) * U(i,2)             &
        &       + RP0(k,2) * V(j,2) * U(i,1) )           &
        & - conjg(h0) * RP0(k,3) * V(j,2) * U(i,2)   )  * oosqrt2
 
  coupL = coupL * (0._dp,1._dp)
  coupR = coupR * (0._dp,-1._dp) 

  Iname = Iname - 1

 End Subroutine Coup_Chip_Chim_P0_NMSSM

 Subroutine Coup_Chip_Chim_S0_NMSSM(i, j, k, U, V, RS0, g, h0, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos and scalar bosons 
 ! valid for the MSSM
 !  i,j ........ index of chargino
 !  k .......... index of the pseudo scalar
 !  U,V ........ mixing matrices of the chargino
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  g .......... SU(2) gauge coupling
 !  h0 ......... trilinear coupling h_u H_d phi
 ! output
 !  coupL ...... the left coupling(i,j)
 !  coupR ...... the right coupling(i,j)
 ! written by Werner Porod, 5.8.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) ::  g, RS0(3,3)
  Complex(dp), Intent(in) :: h0, U(2,2), V(2,2)
  Complex(dp), Intent(out) :: coupL, coupR
  Integer, Intent(in) :: i, j, k

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_Chip_Chim_S0_NMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.2)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.3)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  coupL = ( g * ( RS0(k,1) * Conjg( V(i,1) * U(j,2) )             &
        &       + RS0(k,2) * Conjg( V(i,2) * U(j,1) ) )           &
        & - h0 * RS0(k,3) * Conjg( V(i,2) * U(j,2) )    ) * oosqrt2
  coupR = ( g * ( RS0(k,1) * V(j,1) * U(i,2)             &
        &       + RS0(k,2) * V(j,2) * U(i,1) )           &
        & - Conjg(h0) * RS0(k,3) *  V(j,2) * U(i,2)    ) * oosqrt2
 
  coupL = - coupL
  coupR = - coupR 

  Iname = Iname - 1

 End Subroutine Coup_Chip_Chim_S0_NMSSM

 Subroutine Coup_Chip_Chi0_Sm_NMSSM(k, i, j, N, U, V, RSpm &
                                  &, gp, g, h0, coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between charginos, neutralinos, and
 ! charged scalar bosons 
 ! valid for the NMSSM
 !  i .......... index of chargino
 !  j .......... index of the neutralino
 !  k .......... index of the charged scalar
 !  N.. ........ mixing matrix of the neutralino
 !  U,V ........ mixing matrices of the chargino
 !  RSpm(i,j) .. mixing matrix of charged scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  h0 ......... trilinear coupling H_u H_d phi
 ! output
 !  coupL ...... the left coupling(k,i,j)
 !  coupR ...... the right coupling(k,i,j)
 ! the lagrangian is given by
 !  S^-(k) \bar{\chim(i)} (coupL P_L + coupR P_R) \chi0(j)
 ! written by Werner Porod, 25.03.05
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) ::  i, j, k
  Real(dp), Intent(in) :: g, gp
  Complex(dp), Intent(in) :: N(:,:), U(:,:), V(:,:), RSpm(:,:), h0
  Complex(dp), Intent(out) :: coupL, coupR

  Complex(dp) :: sumI, sumIC

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_Chip_Chi0_Sm_NMSSM'

  If ((i.Lt.1).Or.(i.Gt.2)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Chargino index out of range (i,n_char) = ',i
   Call TerminateProgram
  Elseif ((j.Lt.1).Or.(j.Gt.5)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Neutralino index out of range (i,n_neut) = ',j
   Call TerminateProgram
  Elseif ((k.Lt.1).Or.(k.Gt.2)) Then
   Write(ErrCan,*) 'Error in Subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'Scalar index out of range (i,n_Spm) = ',k
   Call TerminateProgram
  End If

  coupL = ZeroC
  coupR = ZeroC

  sumI = oosqrt2 * ( gp * N(j,1) + g * N(j,2) )
  sumIC = Conjg(sumI)

  coupL = - RSpm(k,2) * ( g * Conjg( V(i,1) * N(j,4) ) &
        &               + Conjg( V(i,2) ) * sumIC)     &
        & - h0 * RSpm(k,1) * Conjg( V(i,2) * N(j,5) )
  coupR = RSpm(k,1) * (- g * U(i,1) * N(j,3) + U(i,2) * sumI )     &
      & - Conjg(h0) * RSpm(k,2) * U(i,2) * N(j,5)
 
  Iname = Iname - 1

 End Subroutine Coup_Chip_Chi0_Sm_NMSSM

 Subroutine Coup_Chi0_Chi0_P0_NMSSM(i, j, k, N, RP0, gp, g, h0, lam, coupL,coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the NMSSM
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RP0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u Phi
 !  lam ........ trilinear coupling Phi^3
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 10.03.2005
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k
  Real(dp), Intent(in) :: g, gp, RP0(3,3)
  Complex(dp), Intent(in) :: N(5,5), h0, lam
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: n_P0, n_N
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_Chi0_Chi0_P0_NMSSM'

  n_P0 = Size(RP0, Dim=1)
  n_N = Size(N, Dim=1)

  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_P0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = 0.5_dp * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = 0.5_dp * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  !---------------------------------
  ! gauge part
  !---------------------------------
   coupL = sumI_i * ( Conjg(N(j,3)) * RP0(k,1) - Conjg(N(j,4)) * RP0(k,2) ) &
      &  + sumI_j * ( Conjg(N(i,3)) * RP0(k,1) - Conjg(N(i,4)) * RP0(k,2) )

  !---------------------------------
  ! Yukawa part
  !---------------------------------
  coupL = coupL + oosqrt2 * h0                                           &
      &   * ( (Conjg(N(i,3)*N(j,4)) + Conjg(N(j,3)*N(i,4))) * RP0(k,3)   &
      &     + (Conjg(N(i,3)*N(j,5)) + Conjg(N(j,3)*N(i,5))) * RP0(k,2)   &
      &     + (Conjg(N(i,4)*N(j,5)) + Conjg(N(j,4)*N(i,5))) * RP0(k,1) ) &
      & - oosqrt2 * lam * Conjg(N(j,5)*N(i,5)) * RP0(k,3)


  coupL = (0._dp,1._dp) * coupL
  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine Coup_Chi0_Chi0_P0_NMSSM


 Subroutine Coup_Chi0_Chi0_S0_NMSSM(i, j, k, N, RS0, gp, g, h0, lam &
               & , coupL, coupR)
 !-----------------------------------------------------------------------
 ! calculates the coupling between Neutralinos and Scalar bosons 
 ! valid for the MSSM, 1-generation epsilon model, and 
 ! 3-generation epsilon model
 ! input
 !  i,j ........ index of Neutralino
 !  k .......... index of the scalar
 !  N .......... mixing matrix of the Neutralino
 !  RS0(i,j) ... mixing matrix of Scalar bosons
 !  g .......... SU(2) gauge coupling
 !  gp ......... U(1) gauge coupling
 !  h0 ......... trilinear coupling H_d H_u Phi
 !  lam ........ trilinear coupling Phi^3
 ! output
 !  coupL ...... the left coupling(i,j,k)
 !  coupR ...... the right coupling(i,j,k)
 ! written by Werner Porod, 10.09.2004
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j, k
  Real(dp), Intent(in) :: g, gp, RS0(3,3)
  Complex(dp), Intent(in) :: N(5,5), h0, lam
  Complex(dp), Intent(out) :: coupL, coupR

  Integer :: n_S0, n_N
  Complex(dp) :: sumI,sumI_i,sumI_j

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_Chi0_Chi0_S0_NMSSM'

  n_S0 = Size(RS0, Dim=1)
  n_N = Size(N, Dim=1)

  !----------------------------------
  ! Checking for correct index range
  !----------------------------------
  If ((i.Lt.1).Or.(i.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index i out of range:',i
   Call TerminateProgram

  Else If ((j.Lt.1).Or.(j.Gt.n_N)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index j out of range:',j
   Call TerminateProgram

  Else If ((k.Lt.1).Or.(k.Gt.n_S0)) Then
   Write(ErrCan,*) 'Problem in subroutine ',NameOfUnit(Iname)
   Write(ErrCan,*) 'index k out of range:',k
   Call TerminateProgram
  End If
 
  sumI = 0.5_dp * ( gp * N(i,1) - g * N(i,2) )
  sumI_i = Conjg( sumI )
  sumI = 0.5_dp * ( gp * N(j,1) - g * N(j,2) )
  sumI_j = Conjg( sumI )

  !---------------------------------
  ! gauge part
  !---------------------------------
   coupL = sumI_i * ( Conjg(N(j,3)) * RS0(k,1) - Conjg(N(j,4)) * RS0(k,2) ) &
      &  + sumI_j * ( Conjg(N(i,3)) * RS0(k,1) - Conjg(N(i,4)) * RS0(k,2) )
  !---------------------------------
  ! Yukawa part
  !---------------------------------
  coupL = coupL + oosqrt2 * h0                                           &
      &   * ( (Conjg(N(i,3)*N(j,4)) + Conjg(N(j,3)*N(i,4))) * RS0(k,3)   &
      &     + (Conjg(N(i,3)*N(j,5)) + Conjg(N(j,3)*N(i,5))) * RS0(k,2)   &
      &     + (Conjg(N(i,4)*N(j,5)) + Conjg(N(j,4)*N(i,5))) * RS0(k,1) ) &
      & - oosqrt2 * lam * Conjg(N(j,5)*N(i,5)) * RS0(k,3)

  coupR = Conjg( coupL )

  Iname = Iname - 1

 End Subroutine Coup_Chi0_Chi0_S0_NMSSM


 Subroutine Coup_P0_Sf_Sf_NMSSM_1(i, j, k, RP0, T3, yuk, Rsf, A &
                                & , mu_in, h0, vevSM, vP, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a pseudoscalar and sfermions
 ! valid for the NMSSM, 1-generation epsilon model, and
 ! 3-generation epsilon model
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RP0(i,j) ... mixing matrix of pseudo scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  mu_in ...... bilinear Higgs parameter of the superpotential
 !  h0 ......... trilinear coupling H_u H_d phi
 !  vevSM ...... (v_d, v_u)
 !  vP ......... <S>
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 22.12.04
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: vevSM(2), vP
  Complex(dp), Intent(in) :: RSf(2,2), yuk, A, mu_in, h0
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RP0(3,3), T3
  Integer, Intent(in) :: i, j, k

  Integer :: i1
  Complex(dp) :: ayuk, aA, muC, mu, parts(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_P0_Sf_Sf_NMSSM_1'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  mu = mu_in + vP * h0 * oosqrt2
  muC = Conjg(mu)

  parts = ZeroC

  If (T3.Gt.0._dp) Then
   parts(1) = yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          & - ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(2) =  A * Rsf(k,2) * Conjg( Rsf(j,1) )     &
          & -  aA * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(3) = oosqrt2 * vevSM(1) &
          & * ( yuk * Conjg(h0) * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          &   - ayuk * h0 * Rsf(k,1) * Conjg( Rsf(j,2) ) )

  Else

   parts(1) =  A * Rsf(k,2) * Conjg( Rsf(j,1) ) &
          & - aA * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(2) = yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          & - ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) )

   parts(3) = oosqrt2 * vevSM(2) &
          & * ( yuk * Conjg(h0) * Rsf(k,2) * Conjg( Rsf(j,1) )  &
          &   - ayuk * h0 * Rsf(k,1) * Conjg( Rsf(j,2) ) )

  End If

  Do i1=1,3
   coup = coup + RP0(i,i1) * parts(i1)
  End Do

  coup = (0._dp, -1._dp) * coup *oosqrt2

  Iname = Iname - 1

 End Subroutine Coup_P0_Sf_Sf_NMSSM_1

 Subroutine Coup_S0_Sf_Sf_NMSSM_1(i, j, k, RS0, T3, e, yuk, Rsf, A, mu_in,&
                      &          vevSM, vP, gU1, gSU2, h0, coup)
 !-----------------------------------------------------------------------
 ! calculates the coupling between a neutral scalar and sfermions
 ! valid for the 1-generation MSSM
 !  i .......... index of pseudo scalar boson
 !  j .......... index of the sfermion \tilde{f}
 !  k .......... index of the sfermion \conjugate{\tilde{f}}
 !  RS0(i,j) ... mixing matrix of scalar bosons
 !  T3 ......... isospin of the left-sfermion
 !  yuk ........ Yukawa coupling
 !  Rsf(i,j) ... mixing matrix of the sfermions
 !  A .......... A-parameter
 !  mu_in ...... bilinear parameters of the superpotential
 !  vevSM(i) ... vevs of the Higgs doublets
 !  vP ......... vev of the Higgs singlet
 !  gU1 ........ U(1) coupling
 !  gSU2 ....... SU(2) coupling
 !  h0 ......... trilinear coupling H_u H_d phi
 ! output
 !  coup ....... the coupling(i,j,k)
 ! written by Werner Porod, 25.01.2005
 !-----------------------------------------------------------------------
 Implicit None

  Complex(dp), Intent(in) :: RSf(2,2) , yuk , A, mu_in, h0
  Complex(dp), Intent(out) :: coup
  Real(dp), Intent(in) :: RS0(3,3), T3, e, vevSM(2), gU1, gSU2, vP
  Integer, Intent(in) :: i, j, k

  Integer :: i1
  Complex(dp) :: ayuk, aA, muC, parts(3), mu
  Real(dp) :: g2, gp2, YL, YR, Dterm

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Coup_S0_Sf_Sf_NMSSM_1'

  coup = ZeroC

  ayuk = Conjg( yuk )
  aA = Conjg( A )

  mu = mu_in + h0 * vP * oosqrt2
  muC = Conjg( mu )
  parts = ZeroC

  g2 = gSU2**2
  gp2 = gU1**2
  YL = e - T3
  YR = - e

  Dterm = 0.5_dp * ( (g2*T3 - YL * gp2) * Rsf(k,1) * Conjg( Rsf(j,1) ) &
        &         - YR * gp2 * Rsf(k,2) * Conjg( Rsf(j,2) ) )

  If (T3.Gt.0._dp) Then
   parts(1) = ( yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & -  vevSM(1) * Dterm

   parts(2) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                   &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2      &
          & + vevSM(2) * Dterm                                        &
          & - vevSM(2) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) ) &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

   parts(3) = 0.5_dp * vevSM(1) * ( ayuk *h0 * Rsf(k,2) * Conjg( Rsf(j,1) ) &
            &            + Conjg(h0) * yuk * Rsf(k,1) * Conjg( Rsf(j,2) ) )
  Else

   parts(1) = - ( A * Rsf(k,2) * Conjg( Rsf(j,1) )                    &
          &     + aA * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2       &
          & - vevSM(1) * Dterm                                         &
          & - vevSM(1) * yuk * ayuk * ( Rsf(k,1) * Conjg( Rsf(j,1) )  &
          &                           + Rsf(k,2) * Conjg( Rsf(j,2) ) )

   parts(2) = ( yuk * muC * Rsf(k,2) * Conjg( Rsf(j,1) )              &
          &   + ayuk * mu * Rsf(k,1) * Conjg( Rsf(j,2) ) ) * oosqrt2  &
          & + vevSM(2) * Dterm
   parts(3) = vevSM(2) * ( ayuk *h0 * Rsf(k,2) * Conjg( Rsf(j,1) )  &
            &            + Conjg(h0) * yuk * Rsf(k,1) * Conjg( Rsf(j,2) ) )

  End If

  Do i1=1,3
   coup = coup + RS0(i,i1) * parts(i1)
  End Do

  Iname = Iname - 1

 End Subroutine Coup_S0_Sf_Sf_NMSSM_1


End Module Couplings

