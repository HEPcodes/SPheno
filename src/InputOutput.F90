Module InputOutput
!--------------------------------------------------------------
! this module contains routines for input/output
!--------------------------------------------------------------

Use Control
Use LHC_observables
Use Experiment
Use LowEnergy
Use Model_Data
Use MSSM_Data
Use NMSSM_Data
Use RP_Data
Use SugraRuns

Interface WriteMatrixBlock
 Module Procedure WriteMatrixBlockC, WriteMatrixBlockR
End Interface

Character(len=15) :: HighScaleModel, HighScaleModel2
! input/output according to SUSY Les Houches accord
Logical, Save :: LesHouches_Format
! if R-parity is added at low energies
Logical, Save ::  Add_Rparity = .False.  
! transfer of GMSB info
Real(dp), Save :: grav_fac = 1._dp
! used in combination with Fittino
Character(len=80), Save :: Old_Data=""
! for input/output at various and different scales
Real(dp), Save :: Qout=0._dp, Qin, Q_PDG_out(36)
Integer, Save :: n_Q_out = 0
Logical, Save :: l_Q_out = .False., l_PDG_out(36) = .False., Calc_Mass=.False.
! using 1st SLHA2 output with flavour ordered states
Logical, Save, Private :: Use_Flavour_States = .False.
! branching ratios larger than BrMin are written out
Real(dp), Save, Private :: BrMin=1.e-4_dp
! cross sections larger than SigMin [in fb] are written out
Real(dp), Save, Private :: SigMin=1.e-3_dp
! contains information on possible inconsitencies in the input
Integer, Save, Private :: in_kont(2)
! name of 'input-program'
Character(len=40), Private :: sp_info 
! tempory variables for Higgs mixing in case of NMSSM
Real(dp), Private, Dimension(3,3) :: RS03_save, RP03_save
! for Les Houches input/output
Logical :: Write_SLHA1 = .False. ! write a second SLHA output file
                                 ! using SLHA1 standard only called SPheno_1.spc
Integer, Private :: i_cpv=0
Logical, Private :: l_RP_Pythia = .False. ! Pythia only takes 4x4 matrix 
                                         ! for neutralinos and 2x2 for charginos
Integer, Private :: io_RP=0 ! 0 standard SLHA way of putting RP state numbering
                            ! 1 flavour ordered way of putting RP state numbering
                            !   using the SLHA number
                            ! 2 using the old SPheno way, only for compatabilty
                            !   with old projects
                            ! is set by using entry 94 in SPhenoInput
Logical, Private :: LWrite_LHC_Observables = .False. ! give LHC observables in the output
Logical, Private :: l_pmns_in = .False. ! in case that PMNS is given
Logical, Private :: minpar_set = .False. ! in case that PMNS is given
! write branching ratios h-> V V^* with off-shell vector bosons instead
! of folding it with the branching ratios of vector bosons
Logical, Private :: BR_Higgs_with_offshell_V = .False. 
! for the nat-SUSY project  with SASHA
Logical, Private :: MADGraph_style = .False. 
! for gravitino mass as input
Logical, Private :: l_m32_in = .False.
                                        
Contains



 Subroutine LesHouches_Input(inFile, kont, HighScaleModel &
                           & , Ecms, Pm, Pp, l_ISR, Fgmsb)
 !--------------------------------------------------------------------
 ! reads in data using the Les Houches standard as defined in 
 ! hep-ph/03111231 and arXiv:0801.0045
 ! input:
 !  - inFile: name of input file
 ! output:
 !   - kont .................. is 0 if the input is consistent, non-zero if
 !                             there is a problem
 !   - HighScaleModel ........ string specifiying the model
 !   - Ecms .................. center of mass energy in GeV
 !   - Pm .................... degree of polarisation of incoming electrons
 !   - Pp .................... degree of polarisation of incoming positrons
 !   - l_ISR ................. if .true. then calculate initial state rediation
 !                             stemming from the incoming elctron/positron beams
 !   - Fgmsb ................. the vev of the F-component in the GMSB model
 !--------------------------------------------------------------------
 Implicit None
  Character(len=60) :: inFile
  Integer, Intent(out) :: kont
  Real(dp), Intent(out) :: Fgmsb, Ecms(:), Pm(:), Pp(:)
  Character(len=15), Intent(out) :: HighScaleModel
  Logical, Intent(out) :: l_ISR(:)
  
  Character(len=80) :: read_line
  Integer :: i_mod=-1, set_mod_par(27)=-1 &
    & , i1, p_max, p_act, i_sp, i_model=-1, i_particles=-1, i_rp=0
  Real(dp) :: wert, Abs_Mu2, cosb2, cos2b, sinb2, RG0(3,3) &
    & , mat_D(3,3), R2(2,2), s12, s13, s23, c12, c13, c23
  Logical :: check, calc_ferm, check_alpha(2), test_l
  Complex(dp) :: lam_vS, vec3C(3), wertC, MMnu(3,3), vec1C(1), mat1(1,1)
  Logical, Save :: l_open = .True.
  Logical :: l_pmns(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = "LesHouches_Input"

  check_alpha = .False. ! used to check consistency of alpha(mZ) calculation
  in_kont = 0
  l_pmns = .False. ! to construct the dim 5 neutrino operator from data
  vec1C = ZeroC

  Call InitializeStandardModel
  Call InitializeLoopFunctions
  !-------------------------------------------
  ! this has been shifted in case of a loop
  !-------------------------------------------
  i_mod = -1
  set_mod_par = -1

  ErrorLevel = -1
  GenerationMixing=.False.
  L_BR=.True.
  L_CS=.False.
  L_ISR = .False. 
  If (l_open) Then
   Open(ErrCan,file="Messages.out",status="unknown")
   l_open = .False.
  End If

  !-------------------------------------------------------
  ! set all model parameters to zero
  !-------------------------------------------------------
  Call Set_All_Parameters_0()
  lam_vs = 0._dp
  sp_info = " "
  !-------------------------------------------------------
  ! necessary to exclude right handed neutrinos from RGEs
  ! is set positive in the corresponding model
  !-------------------------------------------------------
  MNuR = - 1.e-9_dp
  !-------------------------------------------------------
  !take highest precision, will be changed at a later stage
  !-------------------------------------------------------
  TwoLoopRGE = .True.
  !-----------------------------------------------------------------------
  ! these variables are only used in GMSB and will be set correctly below
  !-----------------------------------------------------------------------
  Fgmsb = 1.e12_dp
  m32 = 1.e10_dp ! set an abitrary large gravitino mass in GeV

  kont = 0

  Open(99,file=Trim(inFile),status="old",err=200)

  Do ! reading file
   Read(99,"(a80)",End=200,err=200) read_line
! Write(*,*) trim(read_line)
   If (read_line(1:1).Eq."#") Cycle ! ignore comments for the moment
   If (read_line.Eq." ") Cycle ! ignore empty lines for the moment

   Call PutUpperCase(read_line)
! Write(*,*) trim(read_line)
   If (read_line(1:5).Eq."BLOCK") Then ! assigning values for the select case
    If (read_line(7:12).Eq."MODSEL") Then
     Call Read_MODSEL(99, i_particles, i_model, i_cpv, i_rp, kont)
     If (i_cpv.Eq.0) Then ! one has to recalculated the CKM to the real case
                                     ! because InitializeStandardModel assumes a non-zero phase
     s12 = lam_wolf
     s23 = s12**2 * A_wolf
     s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2) 
     c12 = Sqrt(1._dp-s12*s12)
     c23 = Sqrt(1._dp-s23*s23)
     c13 = Sqrt(1._dp-s13*s13)

     CKM(1,1) = c12 * c13
     CKM(1,2) = s12 * c13
     CKM(1,3) = s13      
     CKM(2,1) = -s12*c23 -c12*s23*s13 
     CKM(2,2) = c12*c23 -s12*s23*s13 
     CKM(2,3) = s23 * c13
     CKM(3,1) = s12*s23 -c12*c23*s13 
     CKM(3,2) = -c12*s23 - s12*c23*s13 
     CKM(3,3) = c23 * c13
    End If

    Else If (read_line(7:14).Eq."SMINPUTS") Then
     Call Read_SMinput(99)

    Else If (read_line(7:14).Eq."IMMINPAR") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMINPAR") 
      Cycle
     End If
     Call Read_MINPAR(99, 1, i_model, set_mod_par, kont)

    Else If (read_line(7:12).Eq."MINPAR") Then
     Call Read_MINPAR(99, 0, i_model, set_mod_par, kont)

    Else If (read_line(7:14).Eq."IMEXTPAR") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMEXTPAR") 
      Cycle
     End If
     Call Read_EXTPAR(99, 1, i_model, set_mod_par, kont)

    Else If (read_line(7:12).Eq."EXTPAR") Then
     Call Read_EXTPAR(99, 0, i_model, set_mod_par, kont)

    Else If (read_line(7:17).Eq."SPHENOINPUT") Then
     Call Read_SPhenoInput(99)

    Else If (read_line(7:12).Eq."SPINFO") Then
     Call  Read_SPINFO(99, kont)

    Elseif ((read_line(7:12).Eq."HIGMIX").Or.(read_line(7:12).Eq."NMHMIX")) Then
     Call ReadMatrixR(99, 3, RS03_save, "RS03_save", kont)

    Else If ((read_line(7:10).Eq."AMIX").Or.(read_line(7:12).Eq."NMAMIX")) Then
     Call ReadMatrixR(99, 3, RP03_save, "RP03_save", kont)

    Else If (read_line(7:17).Eq."IMRVKAPPAIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMRVKAPPAIN") 
      Cycle
     End If
     Call ReadVectorC(99, 3, eps, 1, "Im(epsilon)", kont)
     set_mod_par(26) = 1

    Else If (read_line(7:15).Eq."RVKAPPAIN") Then
     Call ReadVectorC(99, 3, eps, 0, "Re(epsilon)", kont)
     
    Else If (read_line(7:15).Eq."RVSNVEVIN") Then
     Call ReadVectorR(99, 3, vevL, "v_L", kont)  
     set_mod_par(27) = 1

    Else If (read_line(7:19).Eq."IMRVLAMPBDAIN") Then
     If (i_cpv.Lt.2) Then
       Call Warn_CPV(i_cpv, "IMRVLAMPBDAIN") 
       Cycle
     End If
     Call ReadTensorC(99, 3, Rp_lam, 1, "Im(lambda_ijk)", kont, 12)
     RP_trilinear = .True.
 
    Else If (read_line(7:17).Eq."RVLAMPBDAIN") Then
     Call ReadTensorC(99, 3, Rp_lam, 0, "Re(lambda_ijk)", kont, 12)
     RP_trilinear = .True.
 
     Else If (read_line(7:18).Eq."RVLAMPBDAPIN") Then
      Call ReadTensorC(99, 3, Rp_lamp, 0, "Re(lambda'_ijk)", kont)
      RP_trilinear = .True.

     Else If (read_line(7:20).Eq."IMRVLAMPBDAPIN") Then
      If (i_cpv.Lt.2) Then
       Call Warn_CPV(i_cpv, "IMRVLAMPBDAPIN") 
       Cycle
      End If
      Call ReadTensorC(99, 3, Rp_lamp, 1, "Im(lambda'_ijk)", kont)
      RP_trilinear = .True.

    Else If (read_line(7:12).Eq."VCKMIN") Then
     Call Read_CKM(99,i_cpv)

    Else If (read_line(7:10).Eq."MSL2") Then
     Call ReadMatrixC(99, 3, M2L_pmns, 0, "Re(M2_L)", kont, 1)
     M2L_0_pmns = M2L_pmns
     l_ML = .True.
     set_mod_par(11:13) = 1

    Else If (read_line(7:12).Eq."IMMSL2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSL2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2L_pmns, 1, "Im(M2_L)", kont, 1)
     M2L_0_pmns = M2L_pmns
     l_ML = .True.
     set_mod_par(11:13) = 1

    Else If (read_line(7:10).Eq."MSE2") Then
     Call ReadMatrixC(99, 3, M2E_pmns, 0, "Re(M2_E)", kont, 1)
     M2E_0_pmns = M2E_pmns
     l_ME = .True.
     set_mod_par(14:16) = 1

    Else If (read_line(7:12).Eq."IMMSE2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSE2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2E_pmns, 1, "Im(M2_E)", kont, 1)
     M2E_0_pmns = M2E_pmns
     l_ME = .True.
     set_mod_par(14:16) = 1

    Else If (read_line(7:10).Eq."MSQ2") Then
     Call ReadMatrixC(99, 3, M2Q_sckm, 0, "Re(M2_Q)", kont, 1)
     M2Q_0_sckm = M2Q_sckm
     l_MQ = .True.
     set_mod_par(17:19) = 1

    Else If (read_line(7:12).Eq."IMMSQ2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSQ2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2Q_sckm, 1, "Im(M2_Q)", kont, 1)
     M2Q_0_sckm = M2Q_sckm
     l_MQ = .True.
     set_mod_par(17:19) = 1

    Else If (read_line(7:10).Eq."MSU2") Then
     Call ReadMatrixC(99, 3, M2U_sckm, 0, "Re(M2_U)", kont, 1)
     M2U_0_sckm = M2U_sckm
     l_MU = .True.
     set_mod_par(20:22) = 1

    Else If (read_line(7:12).Eq."IMMSU2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSU2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2U_sckm, 1, "Im(M2_U)", kont, 1)
     M2U_0_sckm = M2U_sckm
     l_MU = .True.
     set_mod_par(20:22) = 1

    Else If (read_line(7:10).Eq."MSD2") Then
     Call ReadMatrixC(99, 3, M2D_sckm, 0, "Re(M2_D)", kont, 1)
     M2D_0_sckm = M2D_sckm
     l_MD = .True.
     set_mod_par(23:25) = 1

    Else If (read_line(7:12).Eq."IMMSD2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSD2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2D_sckm, 1, "Im(M2_D)", kont, 1)
     M2D_0_sckm = M2D_sckm
     l_MD = .True.
     set_mod_par(23:25) = 1

    Else If (read_line(7:11).Eq."M15IN") Then

     Call ReadVectorC2(99, 1, vec1C, 0, "Re(M_15)", kont)
     MTM0 = vec1C(1)

     m_H3 = Abs(MTM0)
     ! as initalization, will be computed more precisely later
     MS15_mH3 = Abs(MTM0)
     MT15_mH3 = Abs(MTM0)
     MZ15_mH3 = Abs(MTM0)
     MTM_GUT = Abs(MTM0)

    Else If (read_line(7:13).Eq."IMM15IN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMM15IN") 
      Cycle
     End If
     Call ReadVectorC2(99, 1, vec1C, 1, "Im(M_15)", kont)
     MTM0 = vec1C(1)
     m_H3 = Abs(MTM0)
     ! as initalization, will be computed more precisely later
     MS15_mH3 = Abs(MTM0)
     MT15_mH3 = Abs(MTM0)
     MZ15_mH3 = Abs(MTM0)
     MTM_GUT = Abs(MTM0)

    Else If (read_line(7:11).Eq."M15T15TBIN") Then

     Call ReadVectorC2(99, 1, vec1C, 0, "Re(M_T)", kont)
     MTM0 = vec1C(1)

     m_H3 = Abs(MTM0)
     ! as initalization, will be computed more precisely later
     MS15_mH3 = Abs(MTM0)
     MT15_mH3 = Abs(MTM0)
     MZ15_mH3 = Abs(MTM0)
     MTM_GUT = Abs(MTM0)

    Else If (read_line(7:13).Eq."IMM15T15TBIN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMM15T15TBIN") 
      Cycle
     End If
     Call ReadVectorC2(99, 1, vec1C, 1, "Im(M_T)", kont)
     MTM0 = vec1C(1)
     m_H3 = Abs(MTM0)
     ! as initalization, will be computed more precisely later
     MS15_mH3 = Abs(MTM0)
     MT15_mH3 = Abs(MTM0)
     MZ15_mH3 = Abs(MTM0)
     MTM_GUT = Abs(MTM0)

    Else If (read_line(7:10).Eq."TUIN") Then
     Call ReadMatrixC(99, 3, AU_0_sckm, 0, "Re(T_U)", kont)
     AU_sckm = Transpose(AU_0_sckm) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation
     l_Au = .True.
     set_mod_par(4) = 1

    Else If (read_line(7:12).Eq."IMTUIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTU") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, AU_0_sckm, 1, "Im(T_U)", kont)
     AU_sckm = Transpose(AU_0_sckm) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation
     l_Au = .True.
     set_mod_par(4) = 1


    Else If (read_line(7:10).Eq."TDIN") Then
     Call ReadMatrixC(99, 3, AD_0_sckm, 0, "Re(T_D)", kont)
     AD_sckm = Transpose(AD_0_sckm) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation

     l_Ad = .True.
     set_mod_par(5) = 1

    Else If (read_line(7:12).Eq."IMTDIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTD") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, AD_0_sckm, 1, "Im(T_D)", kont)
     AD_sckm = Transpose(AD_0_sckm) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation
     l_Ad = .True.
     set_mod_par(5) = 1

    Else If (read_line(7:10).Eq."TEIN") Then
     Call ReadMatrixC(99, 3, Al_0_pmns, 0, "Re(T_E)", kont)
     Al_pmns = Transpose(Al_0_pmns) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation
     l_Al = .True.
     set_mod_par(6) = 1

    Else If (read_line(7:12).Eq."IMTEIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTE") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, AL_0_pmns, 1, "Im(T_E)", kont)
     Al_pmns = Transpose(Al_0_pmns) ! unfortunatly there is a transpose due
                                    ! to the RGE implemenation
     l_Al = .True.
     set_mod_par(6) = 1

    Else If (read_line(7:15).Eq."YNURLHUIN")Then

     If (.Not.Ynu_at_MR3) Fixed_Nu_Yukawas = .True.
     Call ReadMatrixC2(99, 3, Y_nu_0, 0,3, "Re(Y_nu_0)", kont)
     If (Maxval(Abs(Y_nu_0))**2.Gt.6._dp) Then
      Write(ErrCan,*) "Y_nu is non-perturbative at M_GUT"
      kont = -314
      Call AddError(314)      
      Return
     End If

    Else If (read_line(7:17).Eq."IMYNURLHUIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYNURLHU") 
      Cycle
     End If
     If (.Not.Ynu_at_MR3) Fixed_Nu_Yukawas = .True.
     Call ReadMatrixC2(99, 3, Y_nu_0, 1,3, "Im(Y_nu_0)", kont)
     If (Maxval(Abs(Y_nu_0))**2.Gt.6._dp) Then
      Write(ErrCan,*) "Y_nu is non-perturbative at M_GUT"
      kont = -314
      Call AddError(314)
      Return
     End If

    Else If (read_line(7:15).Eq."MNURNURIN") Then

     Call ReadVectorR2(99, 3, MnuR, "M_nu_R", kont)

    Else If (read_line(7:13).Eq."UPMNSIN") Then
     l_pmns_in = .True.
     Call Read_PMNS(99)

! Florian Staub
    Else If (read_line(7:11).Eq."Y24IN") Then

     Call ReadMatrixC2(99, 3, Yb3_H24_gut,0,1, "Y_24", kont)
     If (Maxval(Abs(Yb3_H24_gut))**2.Gt.6._dp) Then
      Write(ErrCan,*) "Y_24 is non-perturbative at M_GUT"
      kont = -314
      Call AddError(314)
      Return
     End If
      
     Yb30_H24(3,:,:) = Yb3_H24_gut
     Yb30_H24(2,:,:) = Yb30_H24(3,:,:)
     Yb30_H24(2,:,3) = 0._dp
     Yb30_H24(1,:,:) = Yb30_H24(2,:,:)
     Yb30_H24(1,:,2) = 0._dp
     Yw30_H24 = Yb30_H24
     Yx30_H24 = Yb30_H24


    Else If (read_line(7:13).Eq."IMY24IN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMY24IN") 
      Cycle
     End If

     Call ReadMatrixC2(99, 3, Yb3_H24_gut,1,1, "Im(Y_24)", kont)
     If (Maxval(Abs(Yb3_H24_gut))**2.Gt.6._dp) Then
      Write(ErrCan,*) "Y_24 is non-perturbative at M_GUT"
      kont = -314
      Call AddError(314)
      Return
     End If
     Yb30_H24(3,:,:) = Yb3_H24_gut
     Yb30_H24(2,:,:) = Yb30_H24(3,:,:)
     Yb30_H24(2,:,3) = 0._dp
     Yb30_H24(1,:,:) = Yb30_H24(2,:,:)
     Yb30_H24(1,:,2) = 0._dp
     Yw30_H24 = Yb30_H24
     Yx30_H24 = Yb30_H24


    Else If (read_line(7:16).Eq."YHD15THDIN") Then

     mat1(1,1) = Lambda1_gut
     Call ReadMatrixC2(99, 1, mat1,0,2, "lambda_1", kont)
     Lambda1_gut =  mat1(1,1)
     lam12_0(1) = Lambda1_gut

    Else If (read_line(7:18).Eq."IMYHD15THDIN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYHD15THDIN") 
      Cycle
     End If
     mat1(1,1) = Lambda1_gut
     Call ReadMatrixC2(99, 1, mat1,1,2, "lambda_1", kont)
     Lambda1_gut =  mat1(1,1)
     lam12_0(1) = Lambda1_gut

    Else If (read_line(7:17).Eq."YHU15TBHUIN") Then

     mat1(1,1) = Lambda2_gut
     Call ReadMatrixC2(99, 1, mat1,0,2, "lambda_2", kont)
     Lambda2_gut =  mat1(1,1)
     lam12_0(2) = Lambda2_gut

    Else If (read_line(7:19).Eq."IMYHU15TBHUIN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYHU15TBHUIN") 
      Cycle
     End If
     mat1(1,1) = Lambda2_gut
     Call ReadMatrixC2(99, 1, mat1,1,2, "lambda_2", kont)
     Lambda2_gut =  mat1(1,1)
     lam12_0(2) = Lambda2_gut

    Else If (read_line(7:14).Eq."YL15TLIN") Then

     Call ReadMatrixC2(99, 3, YT_H15_gut,0,2, "Y_T", kont)
     YT_H15_GUT(2,1) = YT_H15_GUT(1,2)
     YT_H15_GUT(3,1) = YT_H15_GUT(1,3)
     YT_H15_GUT(3,2) = YT_H15_GUT(2,3)
     YT0_H15 = YT_H15_GUT
     YZ0_H15 = YT_H15_GUT
     YS0_H15 = YT_H15_GUT
     Y_T_0 = YT_H15_gut

    Else If (read_line(7:16).Eq."IMYL15TLIN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYL15TLIN") 
      Cycle
     End If
     Call ReadMatrixC2(99, 3, YT_H15_gut,1,2, "Im(Y_T)", kont)
     YT_H15_GUT(2,1) = YT_H15_GUT(1,2)
     YT_H15_GUT(3,1) = YT_H15_GUT(1,3)
     YT_H15_GUT(3,2) = YT_H15_GUT(2,3)
     YT0_H15 = YT_H15_GUT
     YZ0_H15 = YT_H15_GUT
     YS0_H15 = YT_H15_GUT
     Y_T_0 = YT_H15_gut

    Else If (read_line(7:11).Eq."Y15IN") Then

     Call ReadMatrixC2(99, 3, YT_H15_gut,0,2, "Y_T", kont)
     YT_H15_GUT(2,1) = YT_H15_GUT(1,2)
     YT_H15_GUT(3,1) = YT_H15_GUT(1,3)
     YT_H15_GUT(3,2) = YT_H15_GUT(2,3)
     YT0_H15 = YT_H15_GUT
     YZ0_H15 = YT_H15_GUT
     YS0_H15 = YT_H15_GUT
     Y_T_0 = YT_H15_gut

    Else If (read_line(7:13).Eq."IMY15IN") Then

     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMY15IN") 
      Cycle
     End If
     Call ReadMatrixC2(99, 3, YT_H15_gut,1,2, "Im(Y_T)", kont)
     YT_H15_GUT(2,1) = YT_H15_GUT(1,2)
     YT_H15_GUT(3,1) = YT_H15_GUT(1,3)
     YT_H15_GUT(3,2) = YT_H15_GUT(2,3)
     YT0_H15 = YT_H15_GUT
     YZ0_H15 = YT_H15_GUT
     YS0_H15 = YT_H15_GUT
     Y_T_0 = YT_H15_gut

    Else If (read_line(7:11).Eq."M24IN") Then
     Call ReadMatrixC(99, 3, MWM30,0, "M24IN", kont)
     MWM3_gut = MWM30
     MWM3Running(1,:,:) = MWM30
     MWM3Running(1,3,:) = 0._dp
     MWM3Running(1,2,:) = 0._dp
     MWM3Running(1,:,3) = 0._dp
     MWM3Running(1,:,2) = 0._dp
     MWM3Running(2,:,:) = MWM30
     MWM3Running(2,3,:) = 0._dp
     MWM3Running(2,:,3) = 0._dp
     MWM3Running(3,:,:) = MWM30

     MGM3 = MWM30
     MWM3 = MWM30
     MBM3 = MWM30
     MXM3 = MWM30

     Do i1=1,3
      MassMWM3(i1) = Abs(MWM30(i1,i1))
      MassMGM3(i1) = Abs(MWM30(i1,i1))
      MassMBM3(i1) = Abs(MWM30(i1,i1))
      MassMXM3(i1) = Abs(MWM30(i1,i1))
     End Do

    Else If (read_line(7:13).Eq."IMMWMIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMWMIN") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, MWM30,1, "MWM3", kont)
     MWM3_gut = MWM30
     MWM3Running(1,:,:) = MWM30
     MWM3Running(1,3,:) = 0._dp
     MWM3Running(1,2,:) = 0._dp
     MWM3Running(1,:,3) = 0._dp
     MWM3Running(1,:,2) = 0._dp
     MWM3Running(2,:,:) = MWM30
     MWM3Running(2,3,:) = 0._dp
     MWM3Running(2,:,3) = 0._dp
     MWM3Running(3,:,:) = MWM30

     MGM3 = MWM30
     MWM3 = MWM30
     MBM3 = MWM30
     MXM3 = MWM30

     Do i1=1,3
      MassMWM3(i1) = Abs(MWM30(i1,i1))
      MassMGM3(i1) = Abs(MWM30(i1,i1))
      MassMBM3(i1) = Abs(MWM30(i1,i1))
      MassMXM3(i1) = Abs(MWM30(i1,i1))
     End Do

! Florian Staub

    Else If (read_line(7:14).Eq."HIGGS3IN") Then
     Call Read_Higgs3(99)

    Else If (read_line(7:10).Eq."MASS") Then
     Call Read_MASS(99)

    Else If (read_line(7:11).Eq."FMASS") Then
     Call Read_FMASS(99)

    Else If (read_line(7:11).Eq."FLIFE") Then
     Call Read_FLIFE(99)

    Else If (read_line(7:12).Eq."FCONST") Then
     Call Read_FCONST(99)

    Else If (read_line(7:10).Eq."UMIX") Then
     Call ReadMatrixC(99, 2, U, 0, "Re(U)", kont)

    Else If (read_line(7:12).Eq."IMUMIX") Then
     Call ReadMatrixC(99, 2, U, 1, "Im(U)", kont)

    Else If (read_line(7:10).Eq."VMIX") Then
     Call ReadMatrixC(99, 2, V, 0, "Re(V)", kont)

    Else If (read_line(7:12).Eq."IMVMIX") Then
     Call ReadMatrixC(99, 2, V, 1, "Im(V)", kont)

    Else If ((read_line(7:10).Eq."NMIX").Or. (read_line(7:12).Eq."NMNMIX")) Then
     If (HighScaleModel.Eq."NMSSM") Then
      Call ReadMatrixC(99, 5, N5, 0, "Re(N)", kont)
     Else
      Call ReadMatrixC(99, 4, N, 0, "Re(N)", kont)
     End If

    Else If (read_line(7:12).Eq."IMNMIX") Then
     If (HighScaleModel.Eq."NMSSM") Then
      Call ReadMatrixC(99, 5, N5, 1, "Im(N)", kont)
     Else
      Call ReadMatrixC(99, 4, N, 1, "Im(N)", kont)
     End If

    Else If (read_line(7:13).Eq."STOPMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~t)", kont)
     If (Rsup(1,1).Eq.0._dp) Rsup = Id6C             ! first initialization
     RSup(5:6,5:6) = Cmplx(R2, Aimag(RSup(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSTOPMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~t)", kont)
     If (Rsup(1,1).Eq.0._dp) Rsup = Id6C             ! first initialization
     RSup(5:6,5:6) = Cmplx(Real(RSup(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:13).Eq."SBOTMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~b)", kont)
     If (Rsdown(1,1).Eq.0._dp) Rsdown = Id6C
     RSdown(5:6,5:6) = Cmplx(R2, Aimag(RSdown(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSBOTMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~b)", kont)
     If (Rsdown(1,1).Eq.0._dp) Rsdown = Id6C
     RSdown(5:6,5:6) = Cmplx(Real(RSdown(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:13).Eq."STAUMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~tau)", kont)
     If (Rslepton(1,1).Eq.0._dp) Rslepton = Id6C
     RSlepton(5:6,5:6) = Cmplx(R2, Aimag(RSlepton(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSTAUMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~tau)", kont)
     If (Rslepton(1,1).Eq.0._dp) Rslepton = Id6C
     RSlepton(5:6,5:6) = Cmplx(Real(RSlepton(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:24).Eq."NEUTRINOBOUNDSIN") Then
     Call Read_Neutrino_Bounds(99)

    Else If (read_line(7:19).Eq."STARTDATAFILE") Then
     Read(99,*) Old_Data
     Old_data= Trim(Old_data) ! to avoid trailing blanks

    Else
     If (output_screen) Write(*,*) "Warning, the following block is ignored"
     If (output_screen) Write(*,*) Trim(read_line)
     Write(ErrCan,*) "Warning, the following block is ignored"
     Write(ErrCan,*) Trim(read_line)

    End If
   End If

   If (kont.Ne.0) Return

  End Do 
200 Close(99)

!-----------------------------------------------
! now some checks and additional settings
!-----------------------------------------------
  If ((Phase_mu.Eq.0._dp).And.(Abs(mu).Gt.0._dp)) phase_mu = mu/Abs(mu)

  If (SPA_convention.And.(.Not.tanb_in_at_Q)) Then
   Write(Errcan,*) &
     & "Warning: in case of SPA conventions, tan(beta) should be given at 1 TeV"
   Write(Errcan,*) "and not at m_Z"
   If (Write_warning_to_screen) Then
    Write(*,*) &
     & "Warning: in case of SPA conventions, tan(beta) should be given at 1 TeV"
    Write(*,*) "and not at m_Z"
   End If
  End If

  If (i_particles.Eq.1) Then  ! MSSM particle content
   If (i_model.Eq.0) Then 
    If (i_rp.Ne.0) HighScaleModel2 = "RPexplicit"
    If ((set_mod_par(7).Eq.1).And.(set_mod_par(8).Eq.1)) Then
     If (i_rp.Eq.0) HighScaleModel = "MSSM1"
     If (set_mod_par(9).Eq.1) Then
      Write(ErrCan,*)  "m^2_H1 and m^2_H2 have been specified together with mu"
      Write(ErrCan,*)  "mu will be ignored"
      set_mod_par(9) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(9).Eq.-1) set_mod_par(9) = 0 ! to avoid problems with the sum
     If (set_mod_par(10).Eq.1) Then
      Write(ErrCan,*) &
           &  "m^2_H1 and m^2_H2 have been specified together with mA^2"
      Write(ErrCan,*)  "mA^2 will be ignored"
      in_kont(2) = 1
      set_mod_par(10) = 0
     End If
     If (set_mod_par(10).Eq.-1) set_mod_par(10) = 0 ! to avoid problems with the sum
     !-------------------------------
     ! first guess of mu and B, mA
     !-------------------------------
     If (i_rp.Eq.0) Then
      cosb2 = 1._dp / (1._dp + tanb**2)
      sinb2 = tanb**2 * cosb2
      cos2b = cosb2 - sinb2
      Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2 )/ cos2b - 0.5_dp * mZ2
      If (Abs_mu2.Lt.0._dp) Abs_mu2 = 1.e4_dp
      mu = Sqrt(abs_mu2) * phase_mu
      B = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanb / (1+tanb**2)
      mP02(2) = Abs(B) * (1._dp/tanb + tanb)
      mP0(2) = Sqrt(mp02(2))
     Else
     End If
    Else If ((set_mod_par(9).Eq.1).And.(set_mod_par(10).Eq.1)) Then
!     HighScaleModel = "MSSM"
     If (set_mod_par(7).Eq.1) Then
      Write(ErrCan,*)  "mu and m_A0 have been specified together with M^2_Hd"
      Write(ErrCan,*)  "M^2_Hd will be ignored"
      set_mod_par(7) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(7).Eq.-1) set_mod_par(7) = 0 ! to avoid problems with the sum
     If (set_mod_par(8).Eq.1) Then
      Write(ErrCan,*)  "mu and m_A0 have been specified together with M^2_Hu"
      Write(ErrCan,*)  "M^2_Hu will be ignored"
      set_mod_par(8) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(8).Eq.-1) set_mod_par(8) = 0 ! to avoid problems with the sum
     !-------------------------------
     ! first guess of B, M^2_H_i
     !-------------------------------
     mA2_Q =  mP02(2)
     B = mP02(2) * tanb / (1._dp + tanb**2)

    Else 
     Write(ErrCan,*) "Higgs sector has not been specified consistently. Aborting"
     kont = -307
     Call AddError(307)
     Return
    End If

   End If ! i_mod.eq.0

     kont = -305 ! model has not specified completly

     If ((i_model.Eq.0).And.(Sum(set_mod_par(1:25)).Eq.23)) kont = 0 ! MSSM
     If ((i_model.Eq.1).And.(Sum(set_mod_par(1:5)).Eq.5)) kont = 0 ! mSugra 
     If ((i_model.Eq.2).And.(Sum(set_mod_par(1:5)).Eq.5)) Then ! GMSB 
      kont = 0
      AoY_d_0 = 0._dp
      AoY_l_0 = 0._dp
      AoY_nu_0 = 0._dp
      AoY_u_0 = 0._dp
      Fgmsb = Lambda * MlambdaS ! needed for the calculation of NLSP
      If (grav_fac.Ge.0) Then
       m32 = grav_fac * 2.4e-19_dp * Fgmsb  ! gravitino mass in GeV  
      Else
       m32 = 2.4e-19_dp * Fgmsb  ! gravitino mass in GeV  
      End If

      If (Lambda.Gt.MlambdaS) Then
       Write(ErrCan,*) "Inconsistent GMSB input, Lambda > M_M"
       Write(ErrCan,*) "Lambda: ",Lambda,"M_M: " ,MlambdaS
       kont = -313
       Call AddError(313)
       Return
      End If
     End If
     If ((i_model.Eq.3).And.(Sum(set_mod_par(1:4)).Eq.4)) kont = 0 ! AMSB
     If (kont.Ne.0) Call AddError(Abs(kont))

  Else If (i_particles.Eq.2) Then  ! MSSM + nu_R particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2L_pmns(1,1) = M2L_pmns(3,3)
   If (set_mod_par(12).Ne.1) M2L_pmns(2,2) = M2L_pmns(3,3)
   If (set_mod_par(14).Ne.1) M2E_pmns(1,1) = M2E_pmns(3,3)
   If (set_mod_par(15).Ne.1) M2E_pmns(2,2) = M2E_pmns(3,3)
   If (set_mod_par(17).Ne.1) M2Q_sckm(1,1) = M2Q_sckm(3,3)
   If (set_mod_par(18).Ne.1) M2Q_sckm(2,2) = M2Q_sckm(3,3)
   If (set_mod_par(20).Ne.1) M2U_sckm(1,1) = M2U_sckm(3,3)
   If (set_mod_par(21).Ne.1) M2U_sckm(2,2) = M2U_sckm(3,3)
   If (set_mod_par(23).Ne.1) M2D_sckm(1,1) = M2D_sckm(3,3)
   If (set_mod_par(24).Ne.1) M2D_sckm(2,2) = M2D_sckm(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   A_nu(1,:) = AoY_nu(1,:) * Y_nu(1,:) 

  Else If (i_particles.Eq.5) Then  ! MSSM + 2 nu_R particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2L_pmns(1,1) = M2L_pmns(3,3)
   If (set_mod_par(12).Ne.1) M2L_pmns(2,2) = M2L_pmns(3,3)
   If (set_mod_par(14).Ne.1) M2E_pmns(1,1) = M2E_pmns(3,3)
   If (set_mod_par(15).Ne.1) M2E_pmns(2,2) = M2E_pmns(3,3)
   If (set_mod_par(17).Ne.1) M2Q_sckm(1,1) = M2Q_sckm(3,3)
   If (set_mod_par(18).Ne.1) M2Q_sckm(2,2) = M2Q_sckm(3,3)
   If (set_mod_par(20).Ne.1) M2U_sckm(1,1) = M2U_sckm(3,3)
   If (set_mod_par(21).Ne.1) M2U_sckm(2,2) = M2U_sckm(3,3)
   If (set_mod_par(23).Ne.1) M2D_sckm(1,1) = M2D_sckm(3,3)
   If (set_mod_par(24).Ne.1) M2D_sckm(2,2) = M2D_sckm(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   A_nu(1:2,:) = AoY_nu(1:2,:) * Y_nu(1:2,:) 
   A_h02 = Ao_h02 * h02
   A_lam222 = Ao_lam222 * lam2
   A_lam112 = Ao_lam112 * lam112
   A_lam122 = Ao_lam122 * lam122

  Else If (i_particles.Eq.3) Then  ! NMSSM particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !--------------------------------------
   ! the input is term of an effective mu
   !--------------------------------------
   vP = sqrt2 * lam_vS / h0
   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2L_pmns(1,1) = M2L_pmns(3,3)
   If (set_mod_par(12).Ne.1) M2L_pmns(2,2) = M2L_pmns(3,3)
   If (set_mod_par(14).Ne.1) M2E_pmns(1,1) = M2E_pmns(3,3)
   If (set_mod_par(15).Ne.1) M2E_pmns(2,2) = M2E_pmns(3,3)
   If (set_mod_par(17).Ne.1) M2Q_sckm(1,1) = M2Q_sckm(3,3)
   If (set_mod_par(18).Ne.1) M2Q_sckm(2,2) = M2Q_sckm(3,3)
   If (set_mod_par(20).Ne.1) M2U_sckm(1,1) = M2U_sckm(3,3)
   If (set_mod_par(21).Ne.1) M2U_sckm(2,2) = M2U_sckm(3,3)
   If (set_mod_par(23).Ne.1) M2D_sckm(1,1) = M2D_sckm(3,3)
   If (set_mod_par(24).Ne.1) M2D_sckm(2,2) = M2D_sckm(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   If (      (Maxval(mSpm).Gt.0._dp).And.(Maxval(mP03).Gt.0._dp)            &
      & .And.(Maxval(mS03).Gt.0._dp).And.(Maxval(Abs(RP03_save)).Gt.0._dp) &
      & .And.(Maxval(Abs(RS03_save)).Gt.0._dp) ) Then
!    External_Higgs = .True.
    mP03(1) = mZ
    mP032 = mP03**2
    mS032 = mS03**2
    mSpm(1) = mW
    mSpm2 = mSpm**2

    RS03 = RS03_save
    RG0 = 0._dp
    RG0(1,1) = - 1._dp / Sqrt(1._dp + tanb**2)
    RG0(2,2) = RG0(1,1)
    RG0(3,3) = 1._dp
    RG0(1,2) = tanb * RG0(1,1)
    RG0(2,1) = RG0(1,2)
!    RP03 = Matmul(RP03_save, RG0)
    RP03 = RP03_save
    RSpm = RG0(1:2,1:2)
   End If

  Else If (i_particles.Eq.4) Then  ! NMSSM + nu_R + S  particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam
   A_nu(1,:) = AoY_nu(1,:) * Y_nu(1,:) 
   A_pns = Ao_hpns * h_pns
   If (lam_vs.Ne.0._dp) Then ! calulate either h0 or v_P
    If ((h0.Ne.0._dp).And.(vP.Eq.0._dp)) vp = sqrt2 * lam_vs / h0
    If ((h0.Eq.0._dp).And.(vP.Ne.0._dp)) h0 = sqrt2 * lam_vs / vp
   End If 
   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2L_pmns(1,1) = M2L_pmns(3,3)
   If (set_mod_par(12).Ne.1) M2L_pmns(2,2) = M2L_pmns(3,3)
   If (set_mod_par(14).Ne.1) M2E_pmns(1,1) = M2E_pmns(3,3)
   If (set_mod_par(15).Ne.1) M2E_pmns(2,2) = M2E_pmns(3,3)
   If (set_mod_par(17).Ne.1) M2Q_sckm(1,1) = M2Q_sckm(3,3)
   If (set_mod_par(18).Ne.1) M2Q_sckm(2,2) = M2Q_sckm(3,3)
   If (set_mod_par(20).Ne.1) M2U_sckm(1,1) = M2U_sckm(3,3)
   If (set_mod_par(21).Ne.1) M2U_sckm(2,2) = M2U_sckm(3,3)
   If (set_mod_par(23).Ne.1) M2D_sckm(1,1) = M2D_sckm(3,3)
   If (set_mod_par(24).Ne.1) M2D_sckm(2,2) = M2D_sckm(3,3)

  End If
  !---------------------------------------------
  ! warning if alpha(mZ) and alpha(0) are given
  !---------------------------------------------
  If (check_alpha(1).And.check_alpha(2)) Then
   Write(ErrCan,*) "Warning: alpha(0) and alpha(mZ) have been specified,"
   Write(ErrCan,*) "the consisteny has not been checked!"
   in_kont(1) = 1 
  End If
  !------------------------------------------
  ! recalculate quantities to be sure
  !------------------------------------------
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2

 !---------
 ! W-boson, first rough estimate
 !---------
  mW2 = mZ2 * (0.5_dp + Sqrt(0.25_dp-Alpha_Mz*pi / (sqrt2*G_F*mZ2))) / 0.985_dp

  mW = Sqrt(mW2)      ! mass
  gamW = 2.06_dp     ! width
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2
  !------------------------------------------------------------------
  ! the running fermion masses at m_Z need to be recalculated
  !------------------------------------------------------------------
  Alpha_mZ = Alpha_MSbar(mZ, mW)
  If (calc_ferm) Call CalculateRunningMasses(mf_l, mf_d, mf_u            &
                           &  , Q_light_quarks, alpha_mZ, alphas_mZ, mZ  &
                           &  , mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)

  If (HighScaleModel.Eq."SUGRA_NuR") Then
   Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
   !--------------------------------------------
   ! check if R-neutrino masses are ordered
   !--------------------------------------------
   If (Abs(MnuR(1)).Gt.Abs(MnuR(2))) Then
    wert = MnuR(1)
    vec3C = Y_nu_0(1,:)
    MnuR(1) = MnuR(2)
    Y_nu_0(1,:) = Y_nu_0(2,:)
    MnuR(2) = wert
    Y_nu_0(2,:) = vec3C
   End If
   If (Abs(MnuR(1)).Gt.Abs(MnuR(3))) Then
    wert = MnuR(1)
    vec3C = Y_nu_0(1,:)
    MnuR(1) = MnuR(3)
    Y_nu_0(1,:) = Y_nu_0(3,:)
    MnuR(3) = wert
    Y_nu_0(3,:) = vec3C
   End If
   If (Abs(MnuR(2)).Gt.Abs(MnuR(3))) Then
    wert = MnuR(2)
    vec3C = Y_nu_0(2,:)
    MnuR(2) = MnuR(3)
    Y_nu_0(2,:) = Y_nu_0(3,:)
    MnuR(3) = wert
    Y_nu_0(3,:) = vec3C
   End If
    
   Y_nu_mR(1,:,1) = Y_nu_0(:,1)
   Y_nu_mR(2,:,1:2) = Y_nu_0(:,1:2)
   Y_nu_mR(3,:,:) = Y_nu_0
  End If

  If (HighScaleModel.Eq."SEESAW_III_3G") Then
   !--------------------------------------------
   ! check if 24-plet masses are ordered
   !--------------------------------------------
   If (Abs(MWM30(1,1)).Gt.Abs(MWM30(2,2))) Then
    wertC = MWM30(1,1)
    vec3C = Yb3_H24_gut(1,:)
    MWM30(1,1) = MWM30(2,2)
    Yb3_H24_gut(1,:) = Yb3_H24_gut(2,:)
    MWM30(2,2) = wertC
    Yb3_H24_gut(2,:) = vec3C
   End If
   If (Abs(MWM30(1,1)).Gt.Abs(MWM30(3,3))) Then
    wertC = MWM30(1,1)
    vec3C = Yb3_H24_gut(1,:)
    MWM30(1,1) = MWM30(3,3)
    Yb3_H24_gut(1,:) = Yb3_H24_gut(3,:)
    MWM30(3,3) = wertC
    Yb3_H24_gut(3,:) = vec3C
   End If
   If (Abs(MWM30(2,2)).Gt.Abs(MWM30(3,3))) Then
    wertC = MWM30(2,2)
    vec3C = Yb3_H24_gut(2,:)
    MWM30(2,2) = MWM30(3,3)
    Yb3_H24_gut(2,:) = Yb3_H24_gut(3,:)
    MWM30(3,3) = wertC
    Yb3_H24_gut(3,:) = vec3C
   End If
    
   Yb30_H24(3,:,:) = Yb3_H24_gut
   Yb30_H24(2,:,:) = Yb30_H24(3,:,:)
   Yb30_H24(2,:,3) = 0._dp
   Yb30_H24(1,:,:) = Yb30_H24(2,:,:)
   Yb30_H24(1,:,2) = 0._dp
   Yw30_H24 = Yb30_H24
   Yx30_H24 = Yb30_H24

    MWM3_gut = MWM30
    MWM3Running(1,:,:) = MWM30
    MWM3Running(1,3,:) = 0._dp
    MWM3Running(1,2,:) = 0._dp
    MWM3Running(1,:,3) = 0._dp
    MWM3Running(1,:,2) = 0._dp
    MWM3Running(2,:,:) = MWM30
    MWM3Running(2,3,:) = 0._dp
    MWM3Running(2,:,3) = 0._dp
    MWM3Running(3,:,:) = MWM30

    MGM3 = MWM30
    MWM3 = MWM30
    MBM3 = MWM30
    MXM3 = MWM30

    Do i1=1,3
     MassMWM3(i1) = Abs(MWM30(i1,i1))
     MassMGM3(i1) = Abs(MWM30(i1,i1))
     MassMBM3(i1) = Abs(MWM30(i1,i1))
     MassMXM3(i1) = Abs(MWM30(i1,i1))
    End Do

  End If

  If ((HighScaleModel(1:9).Eq."SUGRA_NuR").And.(D_SO_10.Ne.0._dp)) Then
    mat_D = 0._dp
    mat_D(1,1) = D_SO_10
    mat_D(2,2) = mat_D(1,1)
    mat_D(3,3) = mat_D(1,1)
    M2E_0_pmns = M2E_0_pmns + mat_D
    M2L_0_pmns = M2L_0_pmns - 3._dp * mat_D
    M2_R_0 = M2_R_0 + 5._dp * mat_D
    M2D_0_sckm = M2D_0_sckm - 3._dp * mat_D
    M2Q_0_sckm = M2Q_0_sckm + mat_D
    M2U_0_sckm = M2U_0_sckm + mat_D
    M2_H_0(1) =  M2_H_0(1) + 2._dp * D_SO_10
    M2_H_0(2) =  M2_H_0(2) - 2._dp * D_SO_10
  End If

  If (External_spectrum) Then
   If (GenerationMixing) Then
    Write(*,*) "Sorry, but the use of an external spectrum with generationmixing"
    Write(*,*) "is not implemented"
    Stop 99
   Else
    RSneut = id3C
    phaseGlu = 1._dp
    !----------------------------------------------------
    ! check for the order of sfermions
    !----------------------------------------------------
    If (mSup(1).Gt.mSup(2)) Then
     wert = mSup(2)
     mSup(2) = mSup(1)
     mSup(1) = wert
     RSup(1,1) = 0._dp
     RSup(2,2) = 0._dp
     Rsup(1,2) = 1
     Rsup(2,1) = -1._dp
    End If
    If (mSup(3).Gt.mSup(4)) Then
     wert = mSup(4)
     mSup(4) = mSup(3)
     mSup(3) = wert
     RSup(3,3) = 0._dp
     RSup(4,4) = 0._dp
     Rsup(3,4) = 1
     Rsup(4,3) = -1._dp
    End If
    If (mSdown(1).Gt.mSdown(2)) Then
     wert = mSdown(2)
     mSdown(2) = mSdown(1)
     mSdown(1) = wert
     RSdown(1,1) = 0._dp
     RSdown(2,2) = 0._dp
     Rsdown(1,2) = 1
     Rsdown(2,1) = -1._dp
    End If
    If (mSdown(3).Gt.mSdown(4)) Then
     wert = mSdown(4)
     mSdown(4) = mSdown(3)
     mSdown(3) = wert
     RSdown(3,3) = 0._dp
     RSdown(4,4) = 0._dp
     Rsdown(3,4) = 1
     Rsdown(4,3) = -1._dp
    End If
    If (mSlepton(1).Gt.mSlepton(2)) Then
     wert = mSlepton(2)
     mSlepton(2) = mSlepton(1)
     mSlepton(1) = wert
     RSlepton(1,1) = 0._dp
     RSlepton(2,2) = 0._dp
     Rslepton(1,2) = 1
     Rslepton(2,1) = -1._dp
    End If
    If (mSlepton(3).Gt.mSlepton(4)) Then
     wert = mSlepton(4)
     mSlepton(4) = mSlepton(3)
     mSlepton(3) = wert
     RSlepton(3,3) = 0._dp
     RSlepton(4,4) = 0._dp
     Rslepton(3,4) = 1
     Rslepton(4,3) = -1._dp
    End If
   ! changing neutralino masses to positive values
    If (HighScaleModel.Eq."NMSSM") Then
     Do i1=1,5
      If (mN5(i1).Lt.0._dp) Then
       mN5(i1) = - mN5(i1)
       N5(i1,:) = (0._dp,1._dp) * N5(i1,:)
      End If
     End Do
    Else
     Do i1=1,4
      If (mN(i1).Lt.0._dp) Then
       mN(i1) = - mN(i1)
       N(i1,:) = (0._dp,1._dp) * N(i1,:)
      End If
     End Do
    End If
   End If
  End If

  !----------------------------------------------------------------------------
  ! construction of dim 5 neutrino operator, either with neutrino masses or
  ! with fake masses, as at this stage only the flavour structure is important
  !----------------------------------------------------------------------------
  If (l_pmns(1)) Then
   MMnu = 0._dp
   If (l_pmns(2)) Then
    Forall(i1=1:3) MMnu(i1,i1) = mf_nu(i1)
    fake_m_nu = .False.
   Else
    Forall(i1=1:3) MMnu(i1,i1) = i1
    fake_m_nu = .True.
   End If
   MatNu = 2._dp * Matmul(Matmul(Transpose(Unu),MMnu),Unu)
  Else
   MatNu = 0._dp
  End If

  !----------------------------------------------------------------
  ! check if T_f and A_f given, if yes, then A_f gets overwritten
  !----------------------------------------------------------------
  If (Au_sckm(3,3).Ne.ZeroC) At_save = ZeroC
  If (Ad_sckm(3,3).Ne.ZeroC) Ab_save = ZeroC
  If (Al_pmns(3,3).Ne.ZeroC) Atau_save = ZeroC

  Iname = Iname - 1

 Contains

  Subroutine Read_Neutrino_Bounds(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(1)
      m2_atm_min = wert
     Case(2)
      m2_atm_max = wert
     Case(3)
      tan2_atm_min = wert
     Case(4)
      tan2_atm_max = wert
     Case(5)
      m2_sol_min = wert
     Case(6)
      m2_sol_max = wert
     Case(7)
      tan2_sol_min = wert
     Case(8)
      tan2_sol_max = wert
     Case(9)
      Ue32_min = wert
     Case(10)
      Ue32_max = wert
     Case default
      Write(ErrCan,*) "Reading block NeutrinoBoundsIn"
      Write(ErrCan,*) "Particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned value is",wert
     End Select

    End Do

  200 Return

  End Subroutine Read_Neutrino_Bounds


  Subroutine Warn_CPV(i_cpv, name)
  Implicit None 
   Integer, Intent(in) :: i_cpv
   Character(len=*), Intent(in) :: name
   If (i_cpv.Eq.0) Write(ErrCan,*) "CP violation is switched off"
   If (i_cpv.Eq.1) Write(ErrCan,*) "CP violation beyond CKM is switched off"
   Write(ErrCan,*) "Ignoring block "//Trim(name)
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End Subroutine Warn_CPV


  Subroutine Read_MASS(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(25)
      mS0(1) = wert
      mS02(1) = mS0(1)**2
      mS03(1) = wert
      mS032(1) = mS02(1)
      mS05(1) = wert
      mS052(1) = mS02(1)
     Case(35)
      mS0(2) = wert
      mS02(2) = mS0(2)**2
      mS03(2) = wert
      mS032(2) = mS02(2)
      mS05(2) = wert
      mS052(2) = mS02(2)
     Case(36)
      P0(2)%m = wert
      P0(2)%m2 = wert**2
! to be chekced
      mP0(2) = wert
      mP02(2) = mP0(2)**2
      mP03(2) = wert
      mP032(2) = mP02(2)
      mP05(2) = wert
      mP052(2) = mP02(2)
     Case(37)
      mSpm(2) = wert
      mSpm2(2) = mSpm(2)**2
     Case(45)
      mS03(3) = wert
      mS032(3) = wert**2
     Case(46)
      mP03(3) = wert
      mP032(3) = wert**2
     Case(1000001)
      mSdown(1) = wert
      mSdown2(1) = wert**2
     Case(1000002)
      mSup(1) = wert
      mSup2(1) = wert**2
     Case(1000003)
      mSdown(3) = wert
      mSdown2(3) = wert**2
     Case(1000004)
      mSup(3) = wert
      mSup2(3) = wert**2
     Case(1000005)
      mSdown(5) = wert
      mSdown2(5) = wert**2
     Case(1000006)
      mSup(5) = wert
      mSup2(5) = wert**2
     Case(1000011)
      mSlepton(1) = wert
      mSlepton2(1) = wert**2
     Case(1000012)
      mSneut(1) = wert
      mSneut2(1) = wert**2
     Case(1000013)
      mSlepton(3) = wert
      mSlepton2(3) = wert**2
     Case(1000014)
      mSneut(2) = wert
      mSneut2(2) = wert**2
     Case(1000015)
      mSlepton(5) = wert
      mSlepton2(5) = wert**2
     Case(1000016)
      mSneut(3) = wert
      mSneut2(3) = wert**2
     Case(1000022)
      mN(1) = wert
      mN2(1) = wert**2
      mN5(1) = wert
      mN52(1) = mN2(1)
     Case(1000023)
      mN(2) = wert
      mN2(2) = wert**2
      mN5(2) = wert
      mN52(2) = mN2(2)
     Case(1000024)
      mC(1) = wert
      mC2(1) = wert**2
     Case(1000025)
      mN(3) = wert
      mN2(3) = wert**2
      mN5(3) = wert
      mN52(3) = mN2(3)
     Case(1000035)
      mN(4) = wert
      mN2(4) = wert**2
      mN5(4) = wert
      mN52(4) = mN2(4)
     Case(1000037)
      mC(2) = wert
      mC2(2) = wert**2
     Case(1000045)
      mN5(5) = wert
      mN52(5) = wert**2
     Case(1000021)
      mGlu = wert
     Case(2000001)
      mSdown(2) = wert
      mSdown2(2) = wert**2
     Case(2000002)
      mSup(2) = wert
      mSup2(2) = wert**2
     Case(2000003)
      mSdown(4) = wert
      mSdown2(4) = wert**2
     Case(2000004)
      mSup(4) = wert
      mSup2(4) = wert**2
     Case(2000005)
      mSdown(6) = wert
      mSdown2(6) = wert**2
     Case(2000006)
      mSup(6) = wert
      mSup2(6) = wert**2
     Case(2000011)
      mSlepton(2) = wert
      mSlepton2(2) = wert**2
     Case(2000013)
      mSlepton(4) = wert
      mSlepton2(4) = wert**2
     Case(2000015)
      mSlepton(6) = wert
      mSlepton2(6) = wert**2
     Case default
      Write(ErrCan,*) "Particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned mass is",wert," GeV"
     End Select

    End Do

  200 Return

  End Subroutine Read_MASS


  Subroutine Read_FCONST(io)
  Implicit None
   Integer, Intent(in) :: io

   Integer :: i1, i2

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, i2, wert, read_line
     Select Case(i1)
     Case(311)
      If (i2.Eq.1) FK = wert
     Case(321)
      If (i2.Eq.1) FK = wert
     Case(511)
      If (i2.Eq.1) FB(1) = wert
     Case(521)
      If (i2.Eq.1) FB(1) = wert
     Case(531)
      If (i2.Eq.1) FB(2) = wert
     Case default
      Write(ErrCan,*) "Block FCONST, particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned mass is",wert," GeV"
     End Select

    End Do

  200 Return

  End Subroutine Read_FCONST


  Subroutine Read_FMASS(io)
  Implicit None
   Integer, Intent(in) :: io

   Integer :: i_scheme
   Real(dp) :: scale

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, i_scheme, scale, read_line
     Select Case(i1)
     Case(311)
      MassK0 = wert
     Case(511)
      MassBq(1) = wert
     Case(521)
      MassBm(1) = wert
     Case(531)
      MassBq(2) = wert
     Case(541)
      MassBm(1) = wert
     Case default
      Write(ErrCan,*) "Block FMASS, particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned mass is",wert," GeV"
     End Select

    End Do

  200 Return

  End Subroutine Read_FMASS


  Subroutine Read_FLIFE(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(511)
      TauBq(1) = wert
     Case(521)
      TauBm(1) = wert
     Case(531)
      TauBq(2) = wert
     Case(541)
      TauBm(2) = wert
     Case default
      Write(ErrCan,*) "Block FLIFE, particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned life time is",wert," s"
     End Select

    End Do

  200 Return

  End Subroutine Read_FLIFE


  Subroutine Read_Higgs3(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(1)
      M_H3 = wert
      ! as initalization, will be computed more precisely later
      MS15_mH3 = wert
      MT15_mH3 = wert
      MZ15_mH3 = wert
      MTM0 = wert
      MTM_GUT = MTM0
     Case(2)
      lam12_0(1) = Cmplx(wert, Aimag(lam12_0(1) ), dp)
     Case(3)
      lam12_0(1) = Cmplx(Real(lam12_0(1), dp ), wert, dp) 
     Case(4)
      lam12_0(2) = Cmplx(wert, Aimag(lam12_0(2) ), dp)
     Case(5)
      lam12_0(2) = Cmplx(Real(lam12_0(2), dp ), wert, dp) 
     Case(6)
      If (wert.Eq.1) Then
       Fifteen_plet = .True.
      Else
       Fifteen_plet = .False.
      End If
     End Select

    End Do

   200 Return

  End Subroutine Read_Higgs3

  Subroutine Read_PMNS(io)
  Implicit None
   Integer, Intent(in) :: io

   Real(dp) :: s12, s13, s23, c12, c13, c23, phase, alpha, beta

    l_pmns(1) = .True.

    s12 = 0._dp
    s13 = 0._dp
    s23 = 0._dp
    phase = 0._dp
    alpha = 0._dp
    beta = 0._dp

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)     
     Case(1)
      theta_12 = wert
      s12 = Sin(wert)
     Case(2)
      theta_23 = wert
      s23 = Sin(wert)
     Case(3)
      theta_13 = wert
      s13 = Sin(wert)
     Case(4)
      delta_nu = wert
      phase = wert
     Case(5)
      alpha_nu1 = wert
      alpha = wert
     Case(6)
      alpha_nu2 = wert
      beta = wert
     Case default
     End Select

    End Do

 200 c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    Unu(1,1) = c12 * c13
    Unu(1,2) = s12 * c13
    Unu(2,3) = s23 * c13
    Unu(3,3) = c23 * c13
    If (phase.Ne.0._dp) Then
     Unu(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
     Unu(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    Else
     Unu(1,3) = s13
     Unu(2,1) = -s12*c23 -c12*s23*s13
     Unu(2,2) = c12*c23 -s12*s23*s13
     Unu(3,1) = s12*s23 -c12*c23*s13
     Unu(3,2) = -c12*s23 - s12*c23*s13
    End If
    ! Majorana phases
    If (alpha.Ne.0._dp) Unu(1,:) =  Unu(1,:) * Exp( (0._dp,-0.5_dp) * alpha )
    If (beta.Ne.0._dp) Unu(2,:) =  Unu(2,:) * Exp( (0._dp,-0.5_dp) * beta )

  End Subroutine Read_PMNS


  Subroutine Read_CKM(io, i_cpv)
  Implicit None
   Integer, Intent(in) :: io, i_cpv

   Real(dp) :: s12, s13, s23, c12, c13, c23, phase
    
    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)     
     Case(1)
      lam_wolf = wert
     Case(2)
      A_wolf = wert
     Case(3)
      rho_wolf = wert
     Case(4)
      eta_wolf = wert
     Case default
     End Select

    End Do

 200   s12 = lam_wolf
    s23 = s12**2 * A_wolf
    s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2)

    If (i_cpv.Eq.0) Then
     Write(ErrCan,*) "Warning: CP violation is switched of, ignoring CKM phase."
     phase = 0._dp
    Else
     phase = Atan(eta_wolf/rho_wolf)
    End If


    c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    CKM(1,1) = c12 * c13
    CKM(1,2) = s12 * c13
    CKM(2,3) = s23 * c13
    CKM(3,3) = c23 * c13
    If (phase.Ne.0._dp) Then
     CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
     CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    Else
     CKM(1,3) = s13
     CKM(2,1) = -s12*c23 -c12*s23*s13
     CKM(2,2) = c12*c23 -s12*s23*s13
     CKM(3,1) = s12*s23 -c12*c23*s13
     CKM(3,2) = -c12*s23 - s12*c23*s13
    End If

  End Subroutine Read_CKM

  Subroutine Read_SPINFO(io, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(inout) :: kont

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line

     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sp, read_line

     If (i_sp.Eq.1) Then
      sp_info = Trim(read_line)//" "//Trim(sp_info)
     Else If (i_sp.Eq.2) Then
      sp_info = Trim(sp_info)//" "//Trim(read_line)
     Else If (i_sp.Eq.4) Then ! there is some inconsistency, exit
      kont = -306
      Call AddError(306)
      Iname = Iname - 1
      Return
     Else
      Write(ErrCan,*) "Unknown entry in BLOCK SPINFO, is ignored:"
      Write(ErrCan,*) i_sp, read_line
     End If
    End Do

   200 Return

  End Subroutine Read_SPINFO

  Subroutine Read_MODSEL(io, i_particles, i_model, i_cpv, i_rp, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(out) :: i_particles, i_model, i_cpv, i_rp
   Integer, Intent(inout) :: kont

   Integer :: i_mod, i_test
   Real(dp) :: r_mod
   Character(len=80) :: read_line

   i_cpv = 0
   i_rp = 0

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_test,r_mod ! ,read_line
     If (i_test.Ne.12) Then
      Backspace(io)
      Read(io,*) i_test,i_mod ! ,read_line
     End If

     If (i_test.Eq.1) Then
      i_particles = i_test
      i_model = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.3)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "MSSM, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      Else If (i_mod.Eq.0) Then
       HighScaleModel = "MSSM"
      Else If (i_mod.Eq.1) Then
       HighScaleModel = "mSugra"
       check = SetHighScaleModel("SUGRA")
       minpar_set = .True.

      Else If (i_mod.Eq.2) Then
       HighScaleModel = "GMSB"
       check = SetHighScaleModel("GMSB")
      Else If (i_mod.Eq.3) Then
       HighScaleModel = "AMSB"
       check = SetHighScaleModel("AMSB")
      End If

     Else If (i_test.Eq.3) Then
      External_Higgs = .False.
      mP03 = 0._dp
      mS03 = 0._dp
      mSpm = 0._dp
      RP03_save = 0._dp
      RS03_save = 0._dp
      If (i_mod.Eq.1) Then
       i_particles = i_test
       i_model = i_mod
       HighScaleModel = "NMSSM"
      Else If (i_mod.Eq.2) Then  ! adding SU(5)
       HighScaleModel = "SUGRA_SU5"
       check = SetHighScaleModel("SUGRA_SU5")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.3) Then  ! adding nu_R, one scale only
       HighScaleModel = "SUGRA_NuR1"
       check = SetHighScaleModel("SUGRA_NuR1")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.111) Then  ! adding nu_R, three scales
       HighScaleModel = "SUGRA_NuR"
       check = SetHighScaleModel("SUGRA_NuR")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.5) Then  ! adding 15-plet, using formula by Ana Rossi
       HighScaleModel = "SEESAW_II"
       check = SetHighScaleModel("SEESAW_II")
       i_particles = 1
       i_model = 1
       Fifteen_plet = .True.
      Else If (i_mod.Eq.114) Then  ! adding Higgs triplet 
       HighScaleModel = "SEESAW_II"
       check = SetHighScaleModel("SEESAW_II")
       i_particles = 1
       i_model = 1
       Fifteen_plet = .False.
      Else If (i_mod.Eq.6) Then  ! adding one nu_R
       HighScaleModel = "NURRP1"
       check = SetHighScaleModel("NURRP1")
       i_particles = 2
       i_model = 4
      Else If (i_mod.Eq.7) Then  ! adding one nu_R, S, Phi
       HighScaleModel = "RPspon"
       check = SetHighScaleModel("RPspon")
       i_particles = 4
       i_model = 6
      Else If (i_mod.Eq.8) Then  ! adding two nu_R
       HighScaleModel = "NURRP2"
       check = SetHighScaleModel("NURRP2")
       i_particles = 5
       i_model = 5
! Florian Staub
      Else If (i_mod.Eq.113) Then  ! adding three 24-plet 
       HighScaleModel = "SEESAW_III_3G"
#ifdef SEESAWIII
       check = SetHighScaleModel("SEESAW_III_3G")
       i_particles = 1
       i_model = 1
#else
       Write(ErrCan,*) "You need to compile with the option -DSEESAWIII"
       Write(ErrCan,*) "to run the seesaw III model"
       Write(*,*) "You need to compile with the option -DSEESAWIII"
       Write(*,*) "to run the seesaw III model"
       Call TerminateProgram()
#endif
      Else If (i_mod.Eq.112) Then  ! adding one 15-plet
#ifdef SEESAWIII
       HighScaleModel = "SEESAW_II_SARAH"
       check = SetHighScaleModel("SEESAW_II_SARAH")
       i_particles = 1
       i_model = 1
#else
       Write(ErrCan,*) "You need to compile with the option -DSEESAWIII"
       Write(ErrCan,*) "to run the seesaw II model with 15-plets"
       Write(*,*) "You need to compile with the option -DSEESAWIII"
       Write(*,*) "to run the seesaw II model with 15-plets"
       Call TerminateProgram()
#endif
! Florian Staub
      Else If (i_mod.Ne.0) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "NMSSM, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.4) Then
      If (i_mod.Eq.1) Then
       i_rp = 1

      Else If (i_mod.Ne.0) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.5) Then
      i_cpv = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.2)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.6) Then
      If (i_mod.Eq.0) Then
       GenerationMixing = .False.
      Else If ((i_mod.Ge.1).And.(i_mod.Le.3)) Then
       GenerationMixing = .True.
      Else
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "GenerationMixing, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.11) Then
      If (i_mod.Lt.0) Then
       Write(*,*) &
          & "You must not give a negative step number for entry 11 of Model",i_mod 
       Write(*,*) "It will be ignored"
      Else If (i_mod.Ge.1) Then ! do nothing in case of 0 
       n_Q_out = i_mod
       l_Q_out = .True.
      End If

     Else If (i_test.Eq.12) Then
      Qout = r_mod
      If (.Not.l_Q_out) Then ! set inital value only in case that entry 11
       n_Q_out = 1           ! was not yet set
       l_Q_out = .True.
      End If
       
!      Call SetRGEScale(r_mod**2)  ! set Q_EWSB

     Else If (i_test.Eq.21) Then
      Select Case(i_mod)
      Case(6)
       l_PDG_out(36) = .True. 
      Case(25) ! h0
       l_PDG_out(8) = .True. 
      Case(35) ! H0
       l_PDG_out(9) = .True. 
      Case(36) ! A0
       l_PDG_out(10) = .True. 
      Case(37) ! H+
       l_PDG_out(11) = .True. 
      Case(45) ! H03
       l_PDG_out(36) = .True. 
      Case(46) ! A02
       l_PDG_out(36) = .True. 
      Case(1000001) ! ~d_1
       l_PDG_out(18) = .True. 
      Case(1000002)  ! ~u_1
       l_PDG_out(12) = .True. 
      Case(1000003)  ! ~d_2
       l_PDG_out(20) = .True. 
      Case(1000004)  ! ~u_2
       l_PDG_out(14) = .True. 
      Case(1000005) ! ~d_3
       l_PDG_out(22) = .True. 
      Case(1000006)  ! ~u_3
       l_PDG_out(16) = .True. 
      Case(1000011)  ! ~l_1
       l_PDG_out(24) = .True. 
      Case(1000012)  ! ~nu_1
       l_PDG_out(30) = .True. 
      Case(1000013)  ! ~l_2
       l_PDG_out(26) = .True. 
      Case(1000014)  ! ~nu_2
       l_PDG_out(31) = .True. 
      Case(1000015)  ! ~l_3
       l_PDG_out(28) = .True. 
      Case(1000016)  ! ~nu_3
       l_PDG_out(33) = .True. 
      Case(1000022) ! neutralino_1
       l_PDG_out(4) = .True. 
      Case(1000023) ! neutralino_2
       l_PDG_out(5) = .True. 
      Case(1000024) ! chargino_1
       l_PDG_out(2) = .True. 
      Case(1000025) ! neutralino_3
       l_PDG_out(6) = .True. 
      Case(1000035) ! neutralino_4
       l_PDG_out(7) = .True. 
      Case(1000037) ! chargino_2
       l_PDG_out(3) = .True. 
      Case(1000045) ! neutralino_5
       l_PDG_out(36) = .True. 
      Case(1000021)  ! gluino
       l_PDG_out(1) = .True. 
      Case(2000001) ! ~d_4
       l_PDG_out(19) = .True. 
      Case(2000002) ! ~u_4
       l_PDG_out(13) = .True. 
      Case(2000003) ! ~d_5
       l_PDG_out(21) = .True. 
      Case(2000004) ! ~u_5
       l_PDG_out(15) = .True. 
      Case(2000005) ! ~d_6
       l_PDG_out(23) = .True. 
      Case(2000006) ! ~u_6
       l_PDG_out(17) = .True. 
      Case(2000011) ! ~l_4
       l_PDG_out(25) = .True. 
      Case(2000013)  ! ~l_5
       l_PDG_out(27) = .True. 
      Case(2000015)  ! ~l_6
       l_PDG_out(29) = .True. 
      Case default
       Write(ErrCan,*) "In block MODSEL, entry 21"
       Write(ErrCan,*) "The number",i_mod,"either does not exist or corresponds"
       Write(ErrCan,*) "to a particle with mass below m_Z. It is thus ignored."
      End Select

     End If

    End Do ! i_mod

    If ((i_rp.Eq.1).And.(i_model.Eq.0)) Then
     HighScaleModel = "RPexplicit"
     check = SetHighScaleModel("RPexplicit")
    Else If ((i_rp.Eq.1).And.(i_model.Ge.1).And.(i_model.Le.3)) Then
     Add_Rparity = .True.
    End If

  End Subroutine Read_MODSEL

  Subroutine Read_SPhenoInput(io)
  Implicit None
   Integer, Intent(in) :: io

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    ! This initialization is necessary for the arrar of production infos
    p_max = Size(Ecms)
    p_act = 0
    Ecms = 0._dp
    Pm = 0._dp
    Pp = 0._dp
    l_ISR = .False.
    Do 
     Read(io,*,End=200,err=200) read_line
!     Write(*,*) trim(read_line)
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*,End=200) i_par,wert ! ,read_line
!     write(*,*) i_par,wert,trim(read_line)
     Select Case(i_par)
     Case(1)
      ErrorLevel = Int(wert)

     Case(2)
      If (Int(wert).Ne.0) Then
       SPA_convention = .True.
       Call SetRGEScale(1.e3_dp**2)
      End If

     Case(3)
      If (Int(wert).Ne.0) Then 
       External_Spectrum = .True.
       External_Higgs = .True.
      End If

     Case(4)
      If (Int(wert).Ne.0) Use_Flavour_States = .True.

     Case(5)
      If (Int(wert).Ne.0) FermionMassResummation = .False.

     Case(6)
      If (Int(wert).Ne.0) Then
       Ynu_at_MR3 = .True.
       Fixed_Nu_Yukawas = .False.
      End If

     Case(7)
      If (Int(wert).Ne.0) Then
       Only_1loop_Higgsmass = .True.
      Else
       Only_1loop_Higgsmass = .False.
      End If

     Case(8) ! calculates Masses for extra scales if required
      If (Int(wert).Ne.0) Then
       Calc_Mass = .True.
      Else
       Calc_Mass = .False.
      End If

     Case(9) ! use old version of BoundaryEW
      If (Int(wert).Eq.1) Then
       UseNewBoundaryEW = .False.
      Else
       UseNewBoundaryEW = .True.
      End If

      Case(10) ! use old version to calculate scale
      If (Int(wert).Eq.1) Then
       UseNewScale = .False.
      Else
       UseNewScale = .True.
      End If

     Case(11)  ! whether to calculate  branching ratios or not
      If (Int(wert).Eq.1) L_BR = .True.
      If (Int(wert).Eq.0) L_BR = .False.

     Case(12) ! minimal value such that a branching ratio is written out
      Call SetWriteMinBR(wert)

     Case(13) ! whether the output of h-> V V* should be folded with
              ! branching ratios of the V*
      If (Int(wert).Ne.0) BR_Higgs_with_offshell_V = .True.
      If (Int(wert).Eq.0) BR_Higgs_with_offshell_V = .False.
 
     Case(21)  ! whether to calculate cross sections or not
      If (Int(wert).Eq.1) L_CS = .True.
      If (Int(wert).Eq.0) L_CS = .False.

     Case(22) ! cms energy
      p_act = p_act + 1
      ! this test is necessary to avoid a memory violation
      If (p_act.Le.p_max) Then
       Ecms(p_act) = wert
      Else
       If (output_screen) &
           & Write(*,*) "The number of required points for the calculation"// &
           &  " of cross sections exceeds",p_max
       If (output_screen) &
           & Write(*,*) "Ignoring this information"
       If (output_screen) &
     &  Write(*,*) "Please enlarge the corresponding arrays in the main program."
       Write(ErrCan,*) "The number of required points for the calculation"// &
               &   " of cross sections exceeds",p_max
       Write(ErrCan,*) "Ignoring this information"
       Write(ErrCan,*) &
         &"Please enlarge the corresponding arrays in the main program."
      End If

     Case (23) ! polarisation of incoming e- beam
      If (Abs(wert).Gt.1._dp) Then
       If (output_screen) Write(*,*) &
           & "e- beam polarisation has to between -1 and 1 and not",wert
       If (output_screen) &
           & Write(*,*) "using now unpolarised e- beam"
       Write(ErrCan,*) &
          & "e- beam polarisation has to between -1 and 1 and not",wert
       Write(ErrCan,*) "using now unpolarised e- beam"
       If (p_act.Le.p_max) Pm(p_act) = 0
      Else
       If (p_act.Le.p_max) Pm(p_act) = wert
      End If

     Case (24) ! polarisation of incoming e+ beam
      If (Abs(wert).Gt.1._dp) Then
       If (output_screen) Write(*,*) &
           & "e+ beam polarisation has to between -1 and 1 and not",wert
       If (output_screen) &
           & Write(*,*) "using now unpolarised e+ beam"
       Write(ErrCan,*) &
          & "e+ beam polarisation has to between -1 and 1 and not",wert
       Write(ErrCan,*) "using now unpolarised e+ beam"
       If (p_act.Le.p_max) Pp(p_act) = 0
      Else
       If (p_act.Le.p_max) Pp(p_act) = wert
      End If

     Case(25)
      If ((wert.Eq.1._dp).And.(p_act.Le.p_max)) L_ISR(p_act) = .True.

     Case(26) ! minimal value such that a cross section is written out
      Call SetWriteMinSig(wert)

     Case(31) ! setting a fixed GUT scale if wert > 0
      If (wert.Gt.0._dp) Call SetGUTScale(wert)

     Case(32) ! requires strict unification
      If (Int(wert).Ne.0) check = SetStrictUnification(.True.)

     Case(34) ! precision of mass calculation
      delta_mass = wert

     Case(35) ! maximal number of iterations
      n_run = Int(wert)

     Case(36) ! write out debug information
      If (wert.Eq.0) Then
       WriteOut = .False.
      Else
       WriteOut = .True.
      End If

     Case(37) ! if =1 -> CKM thourgh V_u, if =2 CKM through V_d 
      If ((wert.Eq.1._dp).Or.(wert.Eq.2._dp)) i1 =  SetYukawaScheme(Int(wert))

     Case(38) ! set looplevel of RGEs
      If (wert.Ne.2._dp) Then
       TwoLoopRGE=.False.
      Else
       TwoLoopRGE=.True.
      End If

     Case(39) ! write additional SLHA1 file
      If (wert.Eq.1._dp) Write_SLHA1 = .True.

     Case(40) ! alpha(0)
      check_alpha(2) = .True.
      Alpha = 1._dp / wert

     Case(41) ! Z-boson width
      gamZ = wert

     Case(42) ! W-boson width
      gamW = wert


     Case(80) ! exit for sure with non-zero value if a problem occurs
      If (wert.Eq.1) Non_Zero_Exit = .True.      

     Case(89) ! quick and dirty way to implement model by Suchita Kulkarni
      If (wert.Eq.1) Model_Suchita = .True.      

     Case(90) ! add R-parity at low energies
      If (wert.Eq.1) Add_Rparity = .True.      

     Case(91) ! fit RP parameters such, that neutrino data are o.k.
      If (wert.Eq.1) l_fit_RP_parameters = .True.      

     Case(92) ! for Pythia input
      If (wert.Eq.1) l_RP_Pythia = .True.      

     Case(93) ! calculates cross section in case of RP, only partially implemented
      If (wert.Eq.1) l_CSrp = .True.      

     Case(94) ! calculates cross section in case of RP, only partially implemented
      If ((wert.Eq.1).Or.(wert.Eq.2)) io_RP = wert

     Case(99) ! MADGraph output style, some additional information
      If (wert.Eq.1) MADGraph_style = .True. 

     Case(100) ! use bsstep instead of 
      If (wert.Eq.1) test_l = Set_Use_bsstep_instead_of_rkqs(.True.)

     Case(101) ! use bsstep instead of 
      If (wert.Eq.1) test_l = Set_Use_rzextr_instead_of_pzextr(.True.)

     Case(110) ! write output for LHC observables
      If (wert.Eq.1) Then
       LWrite_LHC_Observables = .True.
      Else
       LWrite_LHC_Observables = .False.
      End If

     Case Default
      If (output_screen) Write(*,*) &
           & "Problem while reading SPhenoInput, ignoring unknown entry" &
           & ,i_par,wert
      Write(Errcan,*) &
          & "Problem while reading  SPhenoInput, ignoring unknown entry" &
               & ,i_par,wert
     End Select ! i_par

    End Do  ! i_par 

   200 Return

  End Subroutine Read_SPhenoInput

  Subroutine Read_EXTPAR(io, i_c, i_model, set_mod_par, kont)
  Implicit None
   Integer, Intent(in) :: io, i_c, i_model
   Integer, Intent(inout) :: kont, set_mod_par(:)

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    If (i_model.Lt.0) Then ! check if model is already defined
     Write(ErrCan,*) &
     "You must first specify the model before the model parameters can be set."
     kont = -303
     Call AddError(-kont)
     Return
    End If
!    Write(ErrCan,*) "Reading EXTPAR"

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) trim(read_line)
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_par,wert ! ,read_line
     !-------------------------
     ! save data
     !-------------------------
     If (i_c.Eq.0) Then
      If (i_par.Ne.1000039) r_extpar(i_par) = wert
      If (i_par.Ne.1000039) in_extpar(i_par,1) = 1
     Else
      If (i_par.Ne.1000039) i_extpar(i_par) = wert
      If (i_par.Ne.1000039) in_extpar(i_par,2) = 1
     End If

     If ((i_par.Eq.0).And.(i_c.Eq.0)) Then
      n_extpar(0) = "scale for input parameters"
      If (i_model.Eq.0) Call SetRGEScale(wert**2)  ! in case of MSSM
      If (i_model.Eq.1) Call SetGUTScale(wert)     ! Sugra
      If (i_model.Eq.3) Call SetGUTScale(wert)     ! AMSB
     Else If (i_par.Eq.1) Then 
      n_extpar(1) = "M_1"
      If (i_c.Eq.0) Mi(1) = Cmplx(wert, Aimag(Mi(1)), dp) 
      If (i_c.Eq.0) Mi_0(1) = Cmplx(wert, Aimag(Mi_0(1)), dp) 
      If (i_c.Eq.1) Mi(1) = Cmplx(Real(Mi(1),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(1) = Cmplx(Real(Mi_0(1),dp), wert, dp)  
      set_mod_par(1) = 1
     Else If (i_par.Eq.2) Then

      n_extpar(2) = "M_2"
      If (i_c.Eq.0) Mi(2) = Cmplx(wert, Aimag(Mi(2)), dp)
      If (i_c.Eq.0) Mi_0(2) = Cmplx(wert, Aimag(Mi_0(2)), dp)
      If (i_c.Eq.1) Mi(2) = Cmplx(Real(Mi(2),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(2) = Cmplx(Real(Mi_0(2),dp), wert, dp)  
      set_mod_par(2) = 1
     Else If (i_par.Eq.3) Then 
      n_extpar(3) = "M_3"
      If (i_c.Eq.0) Mi(3) = Cmplx(wert, Aimag(Mi(3)), dp) 
      If (i_c.Eq.0) Mi_0(3) = Cmplx(wert, Aimag(Mi_0(3)), dp) 
      If (i_c.Eq.1) Mi(3) = Cmplx(Real(Mi(3),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(3) = Cmplx(Real(Mi_0(3),dp), wert, dp)  
      set_mod_par(3) = 1
     Else If (i_par.Eq.11) Then 
      n_extpar(11) = "A_t"
      If (i_c.Eq.0) AoY_u(3,3) = Cmplx(wert, Aimag(AoY_u(3,3)), dp) 
      If (i_c.Eq.1) AoY_u(3,3) = Cmplx(Real(AoY_u(3,3),dp), wert, dp) 
      At_save = AoY_u(3,3)
      AoY_u_0 = AoY_u 
      set_mod_par(4) = 1
     Else If (i_par.Eq.12) Then 
      n_extpar(12) = "A_b"
      If (i_c.Eq.0) AoY_d(3,3) = Cmplx(wert, Aimag(AoY_d(3,3)), dp)
      If (i_c.Eq.1) AoY_d(3,3) = Cmplx(Real(AoY_d(3,3),dp), wert, dp) 
      Ab_save = AoY_d(3,3)
      AoY_d_0 = AoY_d 
      set_mod_par(5) = 1
     Else If (i_par.Eq.13) Then 
      n_extpar(13) = "A_tau"
      If (i_c.Eq.0) AoY_l(3,3) = Cmplx(wert, Aimag(AoY_l(3,3)), dp) 
      If (i_c.Eq.1) AoY_l(3,3) = Cmplx(Real(AoY_l(3,3),dp), wert, dp) 
      Atau_save = AoY_l(3,3)
      AoY_l_0 = AoY_l 
      set_mod_par(6) = 1
     Else If ((i_par.Eq.21).And.(i_c.Eq.0)) Then 
      n_extpar(21) = "M^2_Hd"
      M2_H(1) = wert
      M2_H_0(1) = wert
      set_mod_par(7) = 1

     Else If ((i_par.Eq.22).And.(i_c.Eq.0)) Then 
      n_extpar(22) = "M^2_Hu"
      M2_H(2) = wert
      M2_H_0(2) = wert
      set_mod_par(8) = 1

     Else If (i_par.Eq.23) Then
      n_extpar(23) = "mu"
      If ((i_model.Eq.0).Or.(HighScaleModel.Eq."NMSSM")) Then 
       If (i_c.Eq.0) mu = Cmplx(wert, Aimag(mu), dp)
       If (i_c.Eq.1) mu = Cmplx(Real(mu,dp), wert, dp)
       set_mod_par(9) = 1
      Else If (i_model.Eq.1) Then
       If (i_c.Eq.0) mu = Cmplx(wert, Aimag(mu), dp)
       If (i_c.Eq.1) mu = Cmplx(Real(mu,dp), wert, dp)
       set_mod_par(9) = 1
      Else
       Write(ErrCan,*) &
         &   "mu can only be specified in the general (N)MSSM and mSUGRA"
       Write(ErrCan,*) "and is ignored for the high scale models GMSB and AMSB"
      End If 

     Else If ((i_par.Eq.24).And.(i_c.Eq.0)) Then 
      n_extpar(24) = "M^2_A(Q)"
      If ((i_model.Eq.0).Or.(i_model.Eq.0)) Then 
       mP02(2) = wert
       If (mP02(2).Ge.0._dp) mP0(2) = Sqrt(mP02(2))
       set_mod_par(10) = 1
      Else
       Write(ErrCan,*) "m_A0(Q) can only be specified in the general MSSM and is"
       Write(ErrCan,*) "ignored for high scale models like SUGRA, GMSB or AMSB"
      End If
 
     Else If ((i_par.Eq.25).And.(i_c.Eq.0)) Then 
      n_extpar(25) = "tan(beta)"
      tanb = wert
      tanb_Q = wert
      tanb_in_at_Q = .True.

     Else If ((i_par.Eq.26).And.(i_c.Eq.0)) Then 
      n_extpar(26) = "m_A, pole mass"
      If (i_model.Eq.0) Then 
       HighScaleModel = "pMSSM"
       check = SetHighScaleModel("pMSSM")
       P0(2)%m = wert
       P0(2)%m2 = wert**2
! to be checked
       mP0(2) = wert
       mP02(2) = wert**2
       set_mod_par(10) = 1
      Else If (i_model.Eq.1) Then 
       P0(2)%m = wert
       P0(2)%m2 = wert**2
! to be checked
       mP0(2) = wert
       mP02(2) = wert**2
       set_mod_par(10) = 1
      Else If (HighScaleModel.Eq."RPspon") Then 
       mP0(2) = wert
       mP02(2) = wert**2
      Else
       Write(ErrCan,*) "m_A0 can only be specified in the general MSSM and is"
       Write(ErrCan,*) "ignored for high scale models like SUGRA, GMSB or AMSB"
      End If

     Else If ((i_par.Eq.31).And.(i_c.Eq.0)) Then 
      n_extpar(31) = "M^2_L11"
      M2L_pmns(1,1) = wert**2
      M2L_0_pmns(1,1) = M2L_pmns(1,1)
      set_mod_par(11) = 1
     Else If ((i_par.Eq.32).And.(i_c.Eq.0)) Then 
      n_extpar(32) = "M^2_L22"
      M2L_pmns(2,2) = wert**2
      M2L_0_pmns(2,2) = M2L_pmns(2,2)
      set_mod_par(12) = 1
     Else If ((i_par.Eq.33).And.(i_c.Eq.0)) Then 
      n_extpar(33) = "M^2_L33"
      M2L_pmns(3,3) = wert**2
      M2L_0_pmns(3,3) = M2L_pmns(3,3)
      set_mod_par(13) = 1
     Else If ((i_par.Eq.34).And.(i_c.Eq.0)) Then 
      n_extpar(34) = "M^2_E11"
      M2E_pmns(1,1) = wert**2
      M2E_0_pmns(1,1) = M2E_pmns(1,1)
      set_mod_par(14) = 1
     Else If ((i_par.Eq.35).And.(i_c.Eq.0)) Then 
      n_extpar(35) = "M^2_E22"
      M2E_pmns(2,2) = wert**2
      M2E_0_pmns(2,2) = M2E_pmns(2,2)
      set_mod_par(15) = 1
     Else If ((i_par.Eq.36).And.(i_c.Eq.0)) Then 
      n_extpar(36) = "M^2_E33"
      M2E_pmns(3,3) = wert**2
      M2E_0_pmns(3,3) = M2E_pmns(3,3)
      set_mod_par(16) = 1
     Else If ((i_par.Eq.41).And.(i_c.Eq.0)) Then 
      n_extpar(41) = "M^2_Q11"
      M2Q_sckm(1,1) = wert**2
      M2Q_0_sckm(1,1) = M2Q_sckm(1,1)
      set_mod_par(17) = 1
     Else If ((i_par.Eq.42).And.(i_c.Eq.0)) Then 
      n_extpar(42) = "M^2_Q22"
      M2Q_sckm(2,2) = wert**2
      M2Q_0_sckm(2,2) = M2Q_sckm(2,2)
      set_mod_par(18) = 1
     Else If ((i_par.Eq.43).And.(i_c.Eq.0)) Then 
      n_extpar(43) = "M^2_Q33"
      M2Q_sckm(3,3) = wert**2
      M2Q_0_sckm(3,3) = M2Q_sckm(3,3)
      set_mod_par(19) = 1
     Else If ((i_par.Eq.44).And.(i_c.Eq.0)) Then 
      n_extpar(44) = "M^2_U11"
      M2U_sckm(1,1) = wert**2
      M2U_0_sckm(1,1) = M2U_sckm(1,1)
      set_mod_par(20) = 1
     Else If ((i_par.Eq.45).And.(i_c.Eq.0)) Then 
      n_extpar(45) = "M^2_U22"
      M2U_sckm(2,2) = wert**2
      M2U_0_sckm(2,2) = M2U_sckm(2,2)
      set_mod_par(21) = 1
     Else If ((i_par.Eq.46).And.(i_c.Eq.0)) Then 
      n_extpar(46) = "M^2_U33"
      M2U_sckm(3,3) = wert**2
      M2U_0_sckm(3,3) = M2U_sckm(3,3)
      set_mod_par(22) = 1
     Else If ((i_par.Eq.47).And.(i_c.Eq.0)) Then 
      n_extpar(47) = "M^2_D11"
      M2D_sckm(1,1) = wert**2
      M2D_0_sckm(1,1) = M2D_sckm(1,1)
      set_mod_par(23) = 1
     Else If ((i_par.Eq.48).And.(i_c.Eq.0)) Then 
      n_extpar(48) = "M^2_D22"
      M2D_sckm(2,2) = wert**2
      M2D_0_sckm(2,2) = M2D_sckm(2,2)
      set_mod_par(24) = 1
     Else If ((i_par.Eq.49).And.(i_c.Eq.0)) Then 
      n_extpar(49) = "M^2_D33"
      M2D_sckm(3,3) = wert**2
      M2D_0_sckm(3,3) = M2D_sckm(3,3)
      set_mod_par(25) = 1

     !------------------------------
     ! NMSSM
     !------------------------------
     Else If (i_par.Eq.61) Then 
      n_extpar(61) = "lambda"
      If (i_c.Eq.0) h0 = Cmplx(wert, Aimag(h0),dp) 
      If (i_c.Eq.1) h0 = Cmplx(Real(h0,dp), wert, dp)
     Else If (i_par.Eq.62) Then 
      n_extpar(62) = "kappa"
      ! different convention with respect to Cyril
      If (i_c.Eq.0) lam = Cmplx(2._dp * wert, Aimag(lam),dp)
      If (i_c.Eq.1) lam = Cmplx(Real(lam,dp), 2._dp * wert, dp)
     Else If (i_par.Eq.63) Then 
      n_extpar(63) = "A_lambda"
      If (i_c.Eq.0) Ao_h0 = Cmplx(wert, Aimag(Ao_h0),dp) 
      If (i_c.Eq.1) Ao_lam = Cmplx(Real(Ao_lam,dp), wert, dp)
     Else If (i_par.Eq.64) Then 
      n_extpar(64) = "A_kappa"
      If (i_c.Eq.0) Ao_lam = Cmplx(wert, Aimag(Ao_lam),dp) 
      If (i_c.Eq.1) Ao_lam = Cmplx(Real(Ao_lam,dp), wert, dp)
     Else If (i_par.Eq.65) Then
      n_extpar(65) = "mu_eff"
      If ((HighScaleModel.Eq."NMSSM").Or.(HighScaleModel.Eq."RPspon")) Then 
       If (i_c.Eq.0) lam_vS = Cmplx(wert, Aimag(lam_vS),dp)
       If (i_c.Eq.1) lam_vS = Cmplx(Real(lam_vS,dp), wert, dp)
       set_mod_par(9) = 1
      Else
       Write(ErrCan,*) "Attempt to use i_par == 65"
       Write(ErrCan,*) "mu can only be specified in this way in the NMSSM and"
       Write(ErrCan,*) "is ignored in model "//Trim(HighScaleModel)
      End If
     
     !---------------------------------------
     ! explicit R-parity breaking a la Munoz
     ! and spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.71) Then
        M2_P = wert
      Else If (i_par.Eq.72) Then 
       vP = wert ! has to be replaced by v_R later in case of Munoz
      Else If (i_par.Eq.73) Then 
       Y_nu(1,1) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.74) Then 
       Y_nu(1,2) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.75) Then 
       Y_nu(1,3) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.76) Then 
       AoY_nu(1,1) = wert ! has to be replaced by AoY_nu(1) later
      Else If (i_par.Eq.77) Then 
       AoY_nu(1,2) = wert ! has to be replaced by AoY_nu(1) later
      Else If (i_par.Eq.78) Then 
       AoY_nu(1,3) = wert ! has to be replaced by AoY_nu(1) later
     !---------------------------------------
     ! spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.81) Then
        M2_S = wert
      Else If (i_par.Eq.82) Then
        M2_R(1,1) = wert
      Else If (i_par.Eq.83) Then 
       vS = wert 
      Else If (i_par.Eq.84) Then 
       vR = wert  
      Else If (i_par.Eq.85) Then 
       h_pns = wert  
      Else If (i_par.Eq.86) Then 
       Ao_hpns = wert  
      Else If (i_par.Eq.87) Then 
       M_rs = wert  
      Else If (i_par.Eq.88) Then 
       M_phi = wert  
      Else If (i_par.Eq.89) Then 
       BM_rs = wert  
      Else If (i_par.Eq.90) Then 
       BM_phi = wert  
     !---------------------------------------
     ! explicit R-parity breaking a la Munoz
     ! and spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.91) Then
       If (i_c.Eq.0) h02 = Cmplx(wert, Aimag(h02),dp) 
       If (i_c.Eq.1) h02 = Cmplx(Real(h02,dp), wert, dp)
      Else If (i_par.Eq.92) Then 
        ! different convention  with respect to Cyril
       If (i_c.Eq.0) lam2 = Cmplx(2._dp * wert, Aimag(lam2),dp) 
       If (i_c.Eq.1) lam2 = Cmplx(Real(lam2,dp), 2._dp * wert , dp)
      Else If (i_par.Eq.93) Then 
       If (i_c.Eq.0) lam112 = Cmplx(wert, Aimag(lam112),dp) 
       If (i_c.Eq.1) lam112 = Cmplx(Real(lam112,dp), wert, dp)
      Else If (i_par.Eq.94) Then 
       If (i_c.Eq.0) lam122 = Cmplx(wert, Aimag(lam122),dp) 
       If (i_c.Eq.1) lam122 = Cmplx(Real(lam122,dp), wert, dp)
      Else If (i_par.Eq.95) Then 
       If (i_c.Eq.0) Ao_h02 = Cmplx(wert, Aimag(Ao_h02),dp)
       If (i_c.Eq.1) Ao_h02 = Cmplx(Real(Ao_h02,dp), wert, dp)
      Else If (i_par.Eq.96) Then 
       If (i_c.Eq.0) Ao_lam222 = Cmplx(wert, Aimag(Ao_lam222),dp) 
       If (i_c.Eq.1) Ao_lam222 = Cmplx(Real(Ao_lam222,dp), wert, dp)
      Else If (i_par.Eq.97) Then 
       If (i_c.Eq.0) Ao_lam112 = Cmplx(wert, Aimag(Ao_lam112),dp)
       If (i_c.Eq.1) Ao_lam112 = Cmplx(Real(Ao_lam112,dp), wert, dp) 
      Else If (i_par.Eq.98) Then 
       If (i_c.Eq.0) Ao_lam122 = Cmplx(wert, Aimag(Ao_lam122),dp)
       If (i_c.Eq.1) Ao_lam122 = Cmplx(Real(Ao_lam122,dp), wert, dp)
      Else If (i_par.Eq.99) Then 
       M2_R(2,2) = wert**2
      Else If (i_par.Eq.100) Then 
       vP2 = wert ! has to be replaced by v_R2 later
      Else If (i_par.Eq.101) Then 
       Y_nu(2,1) = wert ! has to be replaced by hnu2(1) later
      Else If (i_par.Eq.102) Then 
       Y_nu(2,2) = wert ! has to be replaced by hnu2(2) later
      Else If (i_par.Eq.103) Then 
       Y_nu(2,3) = wert ! has to be replaced by hnu2(3) later
      Else If (i_par.Eq.104) Then 
       AoY_nu(2,1) = wert ! has to be replaced by AoY_nu2(1) later
      Else If (i_par.Eq.105) Then 
       AoY_nu(2,2) = wert ! has to be replaced by AoY_nu2(2) later
      Else If (i_par.Eq.106) Then 
       AoY_nu(2,3) = wert ! has to be replaced by AoY_nu2(3) later

! Florian Staub, Seesaw II+III
      Else If (i_par.Eq.200) Then 
       If (i_c.Eq.0) MTM0 = Cmplx(wert, Aimag(MTM0),dp)
       If (i_c.Eq.1) MTM0 = Cmplx(Real(MTM0,dp), wert,dp)
       If (i_c.Eq.0) Then
        m_H3 = wert
       ! as initalization, will be computed more precisely later
        MS15_mH3 = wert
        MT15_mH3 = wert
        MZ15_mH3 = wert
       End If
       MTM_GUT = MTM0
      Else If (i_par.Eq.202) Then 
       If (i_c.Eq.0) Lambda1_gut = Cmplx(wert, Aimag(Lambda1_gut),dp)
       If (i_c.Eq.1) Lambda1_gut = Cmplx(Real(Lambda1_gut,dp), wert,dp)
       lam12_0(1) = Lambda1_gut
      Else If (i_par.Eq.203) Then 
       If (i_c.Eq.0) Lambda2_gut = Cmplx(wert, Aimag(Lambda2_gut),dp)
       If (i_c.Eq.1) Lambda2_gut = Cmplx(Real(Lambda2_gut,dp), wert,dp)
       lam12_0(2) = Lambda2_gut

! Florian Staub, Seesaw II+III

     Else If (i_par.Eq.1000039) Then ! gravitino mass
      m32 = wert
      Fgmsb = m32 / 2.4e-19_dp   ! needed for calculation of decay widths
      l_m32_in = .True.

     Else
      If (i_c.Eq.0) Then
       If (output_screen)  Write(*,*) &
        & "Problem while reading EXTPAR, ignoring unknown (unsupported) entry"&
        & ,i_par,wert
       Write(Errcan,*) &
       & "Problem while reading EXTPAR, ignoring unknown  (unsupported) entry"&
       & ,i_par,wert
      Else If (i_c.Eq.1) Then
       If (output_screen)  Write(*,*) &
      & "Problem while reading IMEXTPAR, ignoring unknown (unsupported) entry"&
        & ,i_par,wert
       Write(Errcan,*) &
     & "Problem while reading IMEXTPAR, ignoring unknown  (unsupported) entry"&
       & ,i_par,wert
      End If
     End If
!     Write(errcan,*) i_par,wert
    End Do  ! i_par 

    200 Return

  End Subroutine Read_EXTPAR 


  Subroutine Read_MINPAR(io, i_c, i_model, set_mod_par, kont)
  Implicit None
   Integer, Intent(in) :: io, i_c, i_model
   Integer, Intent(inout) :: kont, set_mod_par(:)

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    If (i_model.Lt.0) Then ! check if model is already defined
     Write(ErrCan,*) &
     "You must first specify the model before the model parameters can be set."
     kont = -303
     Call AddError(-kont)
     Return
    End If

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_par,wert ! ,read_line
! Write(*,*) i_c,i_par,wert ! ,read_line
! Write(*,*) i_model
     If ((i_par.Eq.1).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.1) Then ! mSugra, M_0
       set_mod_par(1) = 1
       M2D_0_sckm = 0._dp
       Do i1=1,3
        M2D_0_sckm(i1,i1) = wert**2
       End Do
       M2E_0_pmns = M2D_0_sckm
       M2L_0_pmns = M2D_0_sckm
       M2_R_0 = M2D_0_sckm
       M2Q_0_sckm = M2D_0_sckm
       M2U_0_sckm = M2D_0_sckm
       M2_H_0 = wert**2
       M2_T_0 = M2_H_0
      Else If (i_model.Eq.2) Then ! GMSB, Lambda
       set_mod_par(1) = 1
       Lambda = wert
      Else If (i_model.Eq.3) Then ! AMSB, common scalar mass
       set_mod_par(1) = 1
       M0_amsb = wert
      End If

     Else If (i_par.Eq.2) Then 
      If (i_model.Eq.1) Then ! mSugra, M_1/2
       set_mod_par(2) = 1
       If (i_c.Eq.0) Mi_0 =  Cmplx(wert, Aimag(Mi_0),dp) 
       If (i_c.Eq.1) Mi_0 =  Cmplx(Real(Mi_0,dp),  wert, dp)
      Else If ((i_model.Eq.2).And.(i_c.Eq.0)) Then ! GMSB, M_M
       set_mod_par(2) = 1
       MlambdaS = wert
      Else If ((i_model.Eq.3).And.(i_c.Eq.0)) Then ! AMSB, Gravitino mass
       set_mod_par(2) = 1
       M_32 = wert
      End If

     Else If ((i_par.Eq.3).And.(i_c.Eq.0)) Then 
      If ((i_model.Ge.0).And.(i_model.Le.3)) Then ! MSSM, mSugra, GMSB, AMSB
       set_mod_par(3) = 1
       tanb = wert
       tanb_mZ = wert
      End If

     Else If (i_par.Eq.4) Then 
      If ((i_model.Ge.0).And.(i_model.Le.3)) Then ! MSSM, mSugra, GMSB, AMSB, sign_mu
       set_mod_par(4) = 1
       If (i_c.Eq.0) phase_mu = Cmplx(wert, Aimag(phase_mu),dp)
       If (i_c.Eq.1) phase_mu = Cmplx(Real(phase_mu, dp),  wert, dp)
      End If

     Else If (i_par.Eq.5) Then 
      If (i_model.Eq.1) Then ! mSugra, A_0
       set_mod_par(5) = 1
       If (i_c.Eq.0) AoY_d_0 = Cmplx(wert, Aimag(AoY_d_0),dp) 
       If (i_c.Eq.1) AoY_d_0 = Cmplx( Real(AoY_d_0,dp) , wert, dp)
       AoY_l_0 = AoY_d_0
       AoY_u_0 = AoY_d_0
       AoY_nu_0 = AoY_d_0
       AoT_0 = AoY_d_0
       Aolam12_0 = AoY_d_0(1,1)
!       If (i_c.Eq.0) Alam12_0 =  Cmplx(wert, Aimag(Alam12_0),dp)
!       If (i_c.Eq.1) Alam12_0 =  Cmplx(Real(Alam12_0,dp) , wert, dp)
      Else If ((i_model.Eq.2).And.(i_c.Eq.0)) Then ! GMSB, n_5
       set_mod_par(5) = 1
       n5plets = wert
       n10plets = 0
      End If

     Else If ((i_par.Eq.6).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.2) Then ! GMSB, Gravitino mass factor
       set_mod_par(6) = 1
       grav_fac = wert
      End If

     Else If ((i_par.Eq.7).And.(i_c.Eq.0)) Then ! SUGRA_SU5, SO(10) scale
      M_SO_10 = wert

     Else If ((i_par.Eq.8).And.(i_c.Eq.0)) Then ! SUGRA_SU5 or SUGRA_NuR, D-term
      D_SO_10 = wert

     Else If (i_par.Eq.9) Then ! SUGRA_SU5, real(lambda(m_GUT))
      If (i_c.Eq.0) lam_0 = Cmplx(wert, Aimag(lam_0),dp) 
      If (i_c.Eq.1) lam_0 = Cmplx( Real(lam_0,dp), wert, dp)
 
     Else If (i_par.Eq.10) Then ! SUGRA_SU5, real(lambda'(m_GUT))
      If (i_c.Eq.0) lamp_0 = Cmplx(wert, Aimag(lamp_0), dp) 
      If (i_c.Eq.1) lamp_0 = Cmplx( Real(lamp_0,dp), wert, dp) 
 
     Else 
      Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
      If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block MINPAR ",i_par 
      If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMMINPAR ",i_par 
      Call AddError(304)
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If 

    End Do ! i_par

  200  Return

  End Subroutine Read_MINPAR

  Subroutine Read_SMinput(io)
  Implicit None
   Integer, Intent(in) :: io
   
   Integer :: i_sm
   Real(dp) :: wert
   Character(len=80) :: read_line

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sm,wert ! ,read_line

     Select Case(i_sm)
     Case(1)
      check_alpha(1) = .True.
      MZ_input = .True.
      Alpha_MZ_MS = 1._dp / wert

     Case(2) ! G_F
      G_F = wert

     Case(3) ! alpha_s(m_Z)
      alphaS_mZ = wert

     Case(4) ! m_Z
      mZ = wert
      mZ2 = mZ**2
      calc_ferm = .True.

     Case(5) ! m_b(m_b)^MSbar
      mf_d(3) = wert
      mf_d2(3) = mf_d(3)**2
      calc_ferm = .True.

     Case(6) ! m_t^pole
      mf_u(3) = wert
      mf_u2(3) = mf_u(3)**2

     Case(7) ! m_tau^pole
      mf_l(3) = wert
      mf_l2(3) = mf_l(3)**2
      calc_ferm = .True.

     Case(8) ! m_nu_3, input is in GeV
      Mf_nu(3) = wert
      l_pmns(2) = .True.

     Case(11) ! electron mass
      mf_l(1) = wert
      mf_l2(1) = wert**2
      calc_ferm = .True.

     Case(12) ! m_nu_1, input is in GeV
      Mf_nu(1) = wert 
      l_pmns(2) = .True.

     Case(13) ! muon mass
      mf_l(2) = wert
      mf_l2(2) = wert**2
      calc_ferm = .True.

     Case(14) ! m_nu_2, input is in eV, transform to GeV
      Mf_nu(2) = wert 
      l_pmns(2) = .True.

     Case(21) ! d-quark mass at 2 GeV
      mf_d(1) = wert
      mf_d2(1) = wert**2
      calc_ferm = .True.

     Case(22) ! u-quark mass at 2 GeV
      mf_u(1) = wert
      mf_u2(1) = wert**2
      calc_ferm = .True.

     Case(23) ! s-quark mass at 2 GeV
      mf_d(2) = wert
      mf_d2(2) = wert**2
      calc_ferm = .True.

     Case(24) ! c-quark mass at Q=m_c
      mf_u(2) = wert
      mf_u2(2) = wert**2
      calc_ferm = .True.

     Case Default
      If (output_screen) &
           & Write(*,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
      Write(ErrCan,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
     End Select

    End Do ! i_sm

  End Subroutine Read_SMinput
 
 End  Subroutine LesHouches_Input


 Subroutine LesHouches_Out(io_L, kont, id_p, names, HighScaleModel, M_GUT &
      & , BRbtosgamma, Bs_ll, Bd_ll, DeltaMBd, DeltaMBs, BrBToSLL         &
      & , BtoSNuNu, BR_Bu_TauNu, R_Bu_TauNu  &
      & , a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma, BrTauToEGamma  &
      & , BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau &
      & , BR_Z_mu_tau, epsK, DeltaMK, K0toPi0NuNu, KptoPipNuNu  &
      & , Rho_parameter, Ecms, Pm, Pp, ISR, SigSup, SigSdown, SigSle &
      & , SigSn, SigChi0, SigC, SigS0, SigSP, SigHp, omega, f_name)
 !--------------------------------------------------------------------
 ! writes out data using the Les Houches standard as defined in 
 ! hep-ph/0311123
 ! input:
 !   - HighScaleModel ........ string specifiying the model 
 !   - M_GUT ................. scale of high scale theory, in general a GUT theory
 !   - BRbtosgamma ........... 10^4 BR(b -> s gamma)
 !   - a_mu .................. SUSY contribution to (g-2)_muon
 !   - Rho_parameter ......... rho parameter
 !   - Ecms .................. center of mass energy in GeV
 !   - Pm .................... degree of polarisation of incoming electrons
 !   - Pp .................... degree of polarisation of incoming positrons
 !   - ISR ................... if .true. then calculate initial state rediation
 !                             stemming from the incoming elctron/positron beams
 !   - Sigsup ................ cross sections of u-type squarks
 !   - SigSdown .............. cross sections of d-type squarks
 !   - SigSle ................ cross sections of sleptons
 !   - SigSn ................. cross sections of sneutrinos
 !   - SigChi0 ............... cross sections of neutralinos
 !   - SigC .................. cross sections of charginos
 !   - SigSP ................. cross sections of neutral Higgs bosons
 !   - SigHp ................. cross section of charged Higgs boson
 !   - omega ................. relic density omega h^2, optional
 !   - f_name ................ alternative name of output file, optional
 !--------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: io_L
  Integer, Intent(inout) :: id_p(:), kont
  Real(dp), Intent(in) :: M_GUT, BRbtosgamma, Bs_ll(3), BrBToSLL, BR_Bu_TauNu &
      & , BtoSNuNu, a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma          &
      & , BrTauToEGamma, BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu, Bd_ll(3)  &
      & , BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau                                  &
      & , epsK, DeltaMK, K0toPi0NuNu, KptoPipNuNu                             &
      & , R_Bu_TauNu, Rho_parameter, Ecms(:), Pm(:)                           &
      & , Pp(:), SigSup(:,:,:), SigSdown(:,:,:), SigSle(:,:,:), SigSn(:,:,:)  &
      & , SigChi0(:,:,:), SigC(:,:,:), SigS0(:,:), SigSP(:,:,:), SigHp(:,:,:)
  Complex(dp), Intent(in) :: DeltaMBd, DeltaMBs
  Character(len=15), Intent(in) :: HighScaleModel
  Character(len=12) :: names(:)
  Logical, Intent(in) :: ISR(:)
  Real(dp), Intent(in), Optional :: Omega
  Character(len=*), Intent(in), Optional :: f_name

  Integer :: i1, i2, i3, id_sle(6), id_su(6) &
      & , id_sd(6), i_zaehl, n_min, ii, jj, p_id, p_id1, p_id2, p_id3
  Complex(dp) :: nr(10,10), mat8C(8,8)
  Real(dp) :: mnr(10), Q, RG0(3,3), mat6R(6,6), mat3R(3,3) &
      & , mat5R(5,5), mat8R(8,8)
  Real(dp) :: MaxCont
  Character(len=10) :: c_snu(3), c_sle(6), c_su(6), c_sd(6), c_slep(6) &
     & , zeit
  Character(len=6) :: c_lm(3), c_lp(3), c_cp(5), c_cm(5), c_c0(8), c_sp(7) &
     & , c_sm(7)
  Character(len=4) :: c_S0(6), c_P0(5)
  Character(len=8) :: Datum
  Integer :: n_n, n_c, n_s0, n_p0, n_spm, n_sl, n_sn, n_sd, n_su 
  Integer :: id_n(8), id_sp(7), id_sm(7), id_S0(6), id_P0(5), id_c(5), id_snu(3)
  Integer, Parameter ::id_A0 = 36, id_Hp = 37                            &
      & , id_W = 24, id_Z = 23, id_ph = 22, id_gl = 21                   &
      & , id_l(3) = (/ 11, 13, 15 /), id_nu(3) = (/ 12, 14, 16 /)        &
      & , id_u(3) = (/ 2, 4, 6 /), id_d(3) = (/ 1, 3, 5 /)               &
      & , id_grav = 1000039, id_glu = 1000021
  Logical :: non_minimal ! checks if there is deviations from the minimal models
  Logical, Save :: l_open = .True. ! in case of a loop I want to open the
                                   ! output only once 
  !-------------------------------------------------------------------
  ! these variables are useful for the case of R-parity violation
  !------------------------------------------------------------------
  Integer :: i1_min
  Logical :: file_exists,  use_new_RP
  !--------------------------------------------------------------------- 
  ! mixing matrices for shifts to super-CKM and super-PMNS basis 
  !--------------------------------------------------------------------- 
  Integer :: ierr, i_errors(1100), id_check(2)
  Real(dp) :: Yu(3), Yd(3), Yl(3)
  Complex(dp), Dimension(3,3) :: CKM_Q, PMNS_Q, RSn_pmns
  Complex(dp), Dimension(6,6) :: RUsq_ckm, RDsq_ckm, RSl_pmns
  !--------------------------------------------------------------------- 
  ! LHC edge variables
  !--------------------------------------------------------------------- 
  Real(dp) :: LHC_observ(50), mSle(6), mSu(6), mSd(6), m4(4), m5(5), m7(7)
  !--------------------------------------------------------------------- 
  ! for output at various scales
  !--------------------------------------------------------------------- 
  Real(dp) :: delta_Q, Q_act, tb_save
  !--------------------------------------------------------------------- 
  ! dummy particles
  !--------------------------------------------------------------------- 
  Type(particle2) :: part2  
  Type(particle23) :: part23

  c_lm(1) = "e-"
  c_lp(1) = "e+"
  c_lm(2) = "mu-"
  c_lp(2) = "mu+"
  c_lm(3) = "tau-"
  c_lp(3) = "tau+"
  c_cp(1) = "chi_1+"
  c_cp(2) = "chi_2+"
  c_cm(1) = "chi_1-"
  c_cm(2) = "chi_2-"
  c_c0(1) = "chi_10"
  c_c0(2) = "chi_20"
  c_c0(3) = "chi_30"
  c_c0(4) = "chi_40"
  id_n(1) = 1000022
  id_n(2) = 1000023
  id_n(3) = 1000025
  id_n(4) = 1000035

  c_sp(1) = "H+"
  c_sm(1) = "H-"
  id_sp = 0
  id_sp(1) = 37
  id_sm = - id_sp  

  n_sd = 6
  n_su = 6
 
  id_c(1) = 1000024
  id_c(2) = 1000037

  If (HighScaleModel.Eq."NMSSM") Then
   c_c0(5) = "chi_50"
   id_n(5) = 1000045
   c_s0(1) = "S0_1"
   c_s0(2) = "S0_2"
   c_s0(3) = "S0_3"
   id_S0(1) = 25
   id_S0(2) = 35
   id_S0(3) = 45
   c_P0(1) = "P0_1"
   c_P0(2) = "P0_2"
   id_P0(1) = 36
   id_P0(2) = 46
   n_n = 5
   n_c = 2
   n_s0 = 3
   n_p0 = 2 
   n_spm = 1 
   n_sl = 6
   n_sn = 3
   n_min = 1   ! first relevant neutralino index

  Else If (HighScaleModel.Eq."RPexplicit") Then
   id_c(4:5) = id_c(1:2)
   id_c(1:3) = - id_l
   c_cp(4:5) = c_cp(1:2)
   c_cp(1:3) = c_lp
   c_cm(4:5) = c_cm(1:2)
   c_cm(1:3) = c_lm

   id_n(4:7) = id_n(1:4)
   id_n(1:3) = id_nu
   c_c0(4:7) = c_c0(1:4)
   c_c0(1) = "nu_1"
   c_c0(2) = "nu_2"
   c_c0(3) = "nu_3"
!   l_CS = .False.
!   L_BR = .False.
   n_n = 7
   n_c = 5
   n_s0 = 5
   n_p0 = 4 
   n_spm = 7 
   n_sl = 0
   n_sn = 0 

   If (io_RP.Eq.1) Then  ! flavour ordered output
    mat5R = Abs(RP05)
    mat5R(1,:) = 0._dp ! this is the Goldstone boson

    Do i1=1,4
     c_P0(i1) = "P0_"//Bu(i1)
     MaxCont = Maxval(mat5R)
     Call FindPosition(5, mat5R, MaxCont, ii, jj)
     Select Case(jj)
     Case(3)
      id_p0(ii-1) = 1000017   ! Im(snu_e)
     Case(4)
      id_p0(ii-1) = 1000018   ! Im(snu_mu)
     Case(5)
      id_p0(ii-1) = 1000019   ! Im(snu_tau)
     Case default
      id_p0(ii-1) = 36        ! A_0
     End Select
     id_p(46+ii-1) =  id_p0(ii-1)
     mat5R(ii,:) = 0._dp
     mat5R(:,jj) = 0._dp
    End Do

    i_zaehl = 1
    mat5R = Abs(RS05)
    Do i1=1,5
     c_s0(i1) = "S0_"//Bu(i1)
     MaxCont = Maxval(mat5R)
     Call FindPosition(5, mat5R, MaxCont, ii, jj)
     Select Case(jj)
     Case(3)
      id_s0(ii) = 1000012   ! Re(snu_e)
     Case(4)
      id_s0(ii) = 1000014   ! Re(snu_mu)
     Case(5)
      id_s0(ii) = 1000016   ! Re(snu_tau)
     Case default
      id_s0(ii) = 15 + i_zaehl*10
      i_zaehl = i_zaehl + 1
     End Select
     id_p(41+ii) = id_s0(ii)
     mat5R(ii,:) = 0._dp
     mat5R(:,jj) = 0._dp
    End Do

    i_zaehl = 1 
    mat8R = Abs(RSpm8)
    mat8R(1,:) = 0._dp ! this is the Goldstone boson
    Do i1=1,7
     c_sp(i1) = "S^+_"//Bu(i1)
     c_sm(i1) = "S^-_"//Bu(i1)
     MaxCont = Maxval(mat8R)
     Call FindPosition(8, mat8R, MaxCont, ii, jj)
     Select Case(jj)
     Case(1,2)
      id_Sp(ii-1) = 37           ! H^+
     Case(3)
      id_sp(ii-1) = -1000011     ! ~e_L
     Case(4)
      id_sp(ii-1) = -1000013     ! ~mu_L
     Case(6)
      id_sp(ii-1) = -2000011     ! ~e_R
     Case(7)
      id_sp(ii-1) = -2000013     ! ~mu_R
     Case default              ! stau_(i_zaehl)
      id_sp(ii-1) = -(1000000 * i_zaehl + 15)
      i_zaehl = i_zaehl + 1
     End Select
     id_p(49+2*(ii-1)) =  id_sp(ii-1)
     id_p(50+2*(ii-1)) = -id_sp(ii-1)
     mat8R(ii,:) = 0._dp
     mat8R(:,jj) = 0._dp
    End Do
    
   Else If (io_RP.Eq.2) Then ! old way of putting 'PDG' numbers

    Do i1=1,4
     c_s0(i1) = "S0_"//Bu(i1)
     c_P0(i1) = "P0_"//Bu(i1)
     id_s0(i1) = 15 + 10 * i1
     id_p0(i1) = 26 + 10 * i1
    End Do
    c_s0(5) = "S0_5"
    id_s0(5) = 65
    Do i1=1,7
     c_sp(i1) = "S^+_"//Bu(i1)
     c_sm(i1) = "S^-_"//Bu(i1)
     id_sp(i1) = 27 + 10 * i1
    End Do

   Else ! standard SLHA output

    m4 = P05(2:5)%m
    ii=4
    mat5R = RP05
    Do i1=1,4
     MaxCont = Maxval(m4)
     Do i2=1,4
      If (m4(i2).Eq.MaxCont) Then
       m4(i2)=0._dp
       If (ii.Eq.4) id_P0(i2) = 1000019
       If (ii.Eq.3) id_P0(i2) = 1000018
       If (ii.Eq.2) id_P0(i2) = 1000017
       If (ii.Eq.1) id_P0(i2) = 36
       RP05(ii+1,:) = mat5R(i2+1,:)
       ii = ii - 1
      End If
     End Do
    End Do
    Do i1=1,4
     c_P0(i1) = "P0_"//Bu(i1)
     id_p(46+i1) =  id_p0(i1)
    End Do
    m5 = S05%m
    ii=5
    mat5R = RS05
    Do i1=1,5
     MaxCont = Maxval(m5)
     Do i2=1,5
      If (m5(i2).Eq.MaxCont) Then
       m5(i2)=0._dp
       If (ii.Eq.5) id_S0(i2) = 1000016
       If (ii.Eq.4) id_S0(i2) = 1000014
       If (ii.Eq.3) id_S0(i2) = 1000012
       If (ii.Eq.2) id_S0(i2) = 35
       If (ii.Eq.1) id_S0(i2) = 25
       RS05(ii,:) = mat5R(i2,:)
       ii = ii - 1
      End If
     End Do
    End Do
    Do i1=1,5
     c_s0(i1) = "S0_"//Bu(i1)
     id_p(41+i1) = id_s0(i1)
    End Do
    m7 = Spm8(2:8)%m
    mat8C = RSpm8
    ii=7
    Do i1=1,7
     MaxCont = Maxval(m7)
     Do i2=1,7
      If (m7(i2).Eq.MaxCont) Then
       m7(i2)=0._dp
       If (ii.Eq.7) id_Sp(i2) = -2000015
       If (ii.Eq.6) id_Sp(i2) = -2000013
       If (ii.Eq.5) id_Sp(i2) = -2000011
       If (ii.Eq.4) id_Sp(i2) = -1000015
       If (ii.Eq.3) id_Sp(i2) = -1000013
       If (ii.Eq.2) id_Sp(i2) = -1000011
       If (ii.Eq.1) id_Sp(i2) = 37
       RSpm8(ii+1,:) = mat8C(i2+1,:)
       ii = ii - 1
      End If
     End Do
    End Do
    Do i1=1,7
     c_sp(i1) = "S^+_"//Bu(i1)
     c_sm(i1) = "S^-_"//Bu(i1)
     id_p(49+2*i1) =  id_sp(i1)
     id_p(50+2*i1) = -id_sp(i1)
    End Do

   End If

   id_sm = - id_sp
   n_min = 4   ! first relevant neutralino index

  Else If (HighScaleModel.Eq."NURRP1") Then
   id_c(4:5) = id_c(1:2)
   id_c(1:3) = - id_l
   c_cp(4:5) = c_cp(1:2)
   c_cp(1:3) = c_lp
   c_cm(4:5) = c_cm(1:2)
   c_cm(1:3) = c_lm

   id_n(4:7) = id_n(1:4)
   id_n(1:3) = id_nu
   c_c0(4:7) = c_c0(1:4)
   c_c0(1) = "nu_1"
   c_c0(2) = "nu_2"
   c_c0(3) = "nu_3"

   n_n = 5
   n_c = 2
   n_s0 = 6
   n_p0 = 5 
   n_spm = 7 
   n_sl = 0
   n_sn = 0 

   i_zaehl = 1 
   mat8R = Abs(RSpm8)
   mat8R(1,:) = 0._dp ! this is the Goldstone boson
   Do i1=1,7
    c_sp(i1) = "S^+_"//Bu(i1)
    c_sm(i1) = "S^-_"//Bu(i1)
    MaxCont = Maxval(mat8R)
    Call FindPosition(8, mat8R, MaxCont, ii, jj)
    Select Case(jj)
    Case(1,2)
     id_Sp(ii) = 37           ! H^+
    Case(3)
     id_sp(i1) = -1000011     ! ~e_L
    Case(4)
     id_sp(i1) = -1000013     ! ~mu_L
    Case(6)
     id_sp(i1) = -2000011     ! ~e_R
    Case(7)
     id_sp(i1) = -2000013     ! ~mu_R
    Case default              ! stau_(i_zaehl)
     id_sp(i1) = -(1000000 * i_zaehl + 15)
     i_zaehl = i_zaehl + 1
    End Select
    mat8R(ii,:) = 0._dp
    mat8R(:,jj) = 0._dp
   End Do

   id_sm = - id_sp
   n_min = 4   ! first relevant neutralino index

  Else
   c_s0(1) = "h0"
   c_s0(2) = "H0"
   id_S0(1) = 25
   id_S0(2) = 35
   c_P0(1) = "A0"
   id_P0(1) = 36
   n_n = 4
   n_c = 2
   n_s0 = 2
   n_p0 = 1 
   n_spm = 1 
   n_sl = 6
   n_sn = 3
   n_min = 1   ! first relevant neutralino index
  End If

  Q = Sqrt( GetRenormalizationScale() )

  Call Date_and_time(datum,zeit)
  If (l_open) Then
   If (Present(f_name)) Then
    Open(io_L,file=Trim(f_name),status="unknown")
   Else
    Open(io_L,file="SPheno.spc",status="unknown")
   End If
   l_open = .False.
  End If
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2 - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno "//version
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) "# in case of problems send email to porod@physik.uni-wuerzburg.de"
   Write(io_L,100) "# Created: "//Datum(7:8)//"."//Datum(5:6)//"."//Datum(1:4) &
     & //",  "//Zeit(1:2)//":"//Zeit(3:4)
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   "//version//"    # version number"
   !-----------------------------------------------
   ! check if somewhere a problem has had happened
   !-----------------------------------------------
   Call GetError(i_errors)
   !--------------------------------------
   ! a numerical problem might have happen
   !--------------------------------------
   If ((i_errors(1)+i_errors(3)+i_errors(5)+i_errors(7)+i_errors(8) &
     & + i_errors(10) + i_errors(12)+ Sum(i_errors(14:19))).Gt.0)   &
     & Write(io_L,100) &
 & "     3               # potential numerical problem, check file Messages.out"
   If (in_kont(1).Eq.1) Write(io_L,99) 3, &
    & "alpha(0) and alpha(mZ) have both been specified without check for"// &
    &  " consistency"
   If (in_kont(2).Eq.1) Write(io_L,99) 3, &
    & "redundant specification in Higgs sector"
   If (kont.Ne.0)   Write(io_L,100)  &
     "     4               # internal problem, see Messages.out for infos"
   Write(io_L,100) "#"
   Write(io_L,100) "Block SPhenoINFO     # SPheno specific information"
   If (TwoLoopRGE) Then
    Write(io_L,100) "    1      2         # using 2-loop RGEs"
   Else 
    Write(io_L,100) "    1      1         # using 1-loop RGEs"
   End If
   If (YukScen.Eq.1) Then
    Write(io_L,100) &
     &"    2      1         # using running masses for boundary conditions at mZ"
   Else
    Write(io_L,100) &
      &"    2      2         # using pole masses for boundary conditions at mZ"
   End If

  If (kont.Ne.0) Then
   Write(*,*) "There has been a problem during the run."
   Write(*,*) "Please check the file Messages.out for further information."
   Return
  End If

   !--------------------------------------
   ! model information
   !--------------------------------------
   If (HighScaleModel.Eq."mSugra") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (i_cpv.Gt.0) Write(io_L,110) 5,i_cpv,"switching on CP violation"
    If (GenerationMixing) Write(io_L,100) &
      &             "    6    3    # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    If (tanb_mZ.Ne.0._dp) Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# cos(phase_mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    Write(io_L,100) "#"
    If ((Aimag(mu).Ne.0._dp).Or.(Aimag(AoY_l_0(1,1)).Ne.0._dp)) Then
     Write(io_L,100) "Block IMMINPAR  # Input parameters, imaginary part"
     Write(io_L,101) 4, Aimag(phase_mu),"# sin(phase_mu)"
     Write(io_L,101) 5, Aimag(AoY_l_0(1,1)),"# Im(A0)"
    End If

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(M_GUT)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(M_GUT)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(M_GUT)^DRbar"
    
   !----------------------------------------------------------------
   ! in the case of GenerationMixing Yukawas
   ! are given in the SuperCKM basis 
   !----------------------------------------------------------------
    If (GenerationMixing) Then
     Call Switch_to_superCKM(Y_d_0, Y_u_0, A_d_0, A_u_0, M2_D_0, M2_Q_0, M2_U_0         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .True. &
              &, CKM_out=CKM_Q, Yd=Yd, Yu=Yu )

     Do i1=1,3
      Yl(i1) = Real(Y_l_0(i1,i1),dp)
     End Do


    Else ! .non.GenerationMixing

     Do i1=1,3
      Yu(i1) = Real(Y_u_0(i1,i1),dp)
      Yd(i1) = Real(Y_d_0(i1,i1),dp)
      Yl(i1) = Real(Y_l_0(i1,i1),dp)
     End Do

    End If

    Write(io_L,106) "Block Yu Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yu(1), "# Y_u(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yu(2), "# Y_c(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yu(3), "# Y_t(M_GUT)^DRbar"

    Write(io_L,106) "Block Yd Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yd(1), "# Y_d(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yd(2), "# Y_s(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yd(3), "# Y_b(M_GUT)^DRbar"

    If (GenerationMixing) Call WriteMatrixBlockC(io_L,3,CKM_Q,M_GUT &
                         & ,"VCKM","V_CKM at the GUT scale","V_(")

    Write(io_L,106) "Block Ye Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yl(1), "# Y_e(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yl(2), "# Y_mu(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yl(3), "# Y_tau(M_GUT)^DRbar"

    If (Model_Suchita) Then
     Call chop(Au_sckm)
     Call chop(Ad_sckm)
     Call WriteMatrixBlockC(io_L,3,Au_sckm,M_GUT,"Tu","M_GUT","T_(u,",tr=.True.)
     Call WriteMatrixBlockC(io_L,3,Ad_sckm,M_GUT,"Td","M_GUT","T_(d,",tr=.True.)
     Call chop(M2Q_sckm)
     Call WriteMatrixBlockC(io_L,3,M2Q_SCKM,M_GUT,"MSQ2" &
           & ,"M^2_Q soft SUSY breaking masses","M^2_(Q,")
     Call chop(M2U_sckm)
     Call WriteMatrixBlockC(io_L,3,M2U_SCKM,M_GUT,"MSU2" &
           & ,"M^2_U soft SUSY breaking masses","M^2_(U,")
     Call chop(M2D_sckm)
     Call WriteMatrixBlockC(io_L,3,M2D_SCKM,M_GUT,"MSD2" &
           & ,"M^2_D soft SUSY breaking masses","M^2_(D,")
    End If ! Model_Suchita

   Else If ((HighScaleModel(1:5).Eq."SUGRA").Or.   &
          & (HighScaleModel(1:5).Eq."SEESA")) Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (i_cpv.Gt.0) Write(io_L,111) 5,i_cpv,"switching on CP violation"
    If (GenerationMixing) Write(io_L,100) &
      &             "    6    3    # switching on flavour violation"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,100) "    3    2    # mSUGRA model + SU(5)"
    Else If (HighScaleModel.Eq."SUGRA_NuR1") Then
     Write(io_L,100) "    3    3    # mSUGRA model + nu_R at a common scale"
     Write(io_L,100) "Block MnuRnuR  # mass scale of the right handed neutrinos"
     Write(io_L,107) 1,1,MnuR(1),"# m_nu_R_1"     
     Write(io_L,107) 2,2,MnuR(1),"# m_nu_R_2"     
     Write(io_L,107) 3,3,MnuR(1),"# m_nu_R_3"     
    Else If (HighScaleModel.Eq."SUGRA_NuR") Then
     Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
     Write(io_L,100) "    3    4    # mSUGRA model + three  nu_R, ~nu_R"
     Write(io_L,100) "Block SEESAWGENERATIONS"
     Write(io_l,98) 1,3,"number of nu^C"

     Call WriteMatrixBlockC3(io_L,3,3,Y_nu_0,m_GUT &
                         & ,"YNURLHU","GUT scale","Y_(nu,")

    Else If (HighScaleModel.Eq."SEESAW_II") Then
     Write(io_L,100) "    3    5    # mSUGRA model + Higgs triplett"
     Write(io_L,100) "Block SEESAWGENERATIONS"
     Write(io_l,98) 15,1,"number of 15-plets"

     Call WriteMatrixBlockC(io_L,3,Y_T_0,m_GUT &                ! symmetric
                         & ,"YT0","GUT scale","Y_(T,", .True.)  ! matrix

     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,M_H3(1),"# m_H3      "
     Write(io_L,101) 2,Real(lam12_0(1),dp),"# Re(lambda_1)"
     If (Aimag(lam12_0(1)).Ne.0._dp) &
         Write(io_L,101) 3,Aimag(lam12_0(1)),"# Im(lambda_1)"
     Write(io_L,101) 4,Real(lam12_0(2),dp),"# Re(lambda_2)"
     If (Aimag(lam12_0(2)).Ne.0._dp) &
         Write(io_L,101) 5,Aimag(lam12_0(2)),"# Im(lambda_2)"
     If (Fifteen_plet) Then
      Write(io_L,100) "    6    1               # using RGEs for 15-plet"
     Else
      Write(io_L,100) "    6    0               # using RGEs for 3-plet"
     End If
! Florian SARAH, Seesaw II+III
    Else If (HighScaleModel.Eq."SEESAW_III_3G") Then
     Write(io_L,100) "    3    113   # mSUGRA model + 3 24-plets"
     Write(io_L,100) "Block SEESAWGENERATIONS"
     Write(io_l,98) 24,3,"number of 24-plets"

     Call WriteMatrixBlockC(io_L,3,MWM3_gut,m_GUT,"M24" &
           & ,"24-plet masses at M_GUT","M_(24,")
     Call WriteMatrixBlockC3(io_L,3,1,Yb3_h24_gut,m_GUT &
                         & ,"Y24","GUT scale","Y_(24,")

!     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
!     Write(io_L,101) 1,Real(MWM3_gut(1,1),dp),"# MWM      "
     Call WriteMatrixBlockC(io_L,3,MBM3running(3,:,:),Abs(MWM30(3,3)),"M24B24B" &
           & ,"singlet mass matrix at M_W3(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,3,MGM3running(3,:,:),Abs(MWM30(3,3)),"M24G24G" &
           & ,"octet mass matrix at M_W3(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,3,MWM3running(3,:,:),Abs(MWM30(3,3)),"M24W24W" &
           & ,"triplet mass matrix at M_W3(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,3,MXM3running(3,:,:),Abs(MWM30(3,3)),"M24X24X" &
           & ,"X mass matrix at M_W3(GUT)","M_(W,")
     Call WriteMatrixBlockC3(io_L,3,1,Yb30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YHU24BL","Triplet scale 3","Y_(b,")
     Call WriteMatrixBlockC3(io_L,3,1,Yw30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YHU24WL","Triplet scale 3","Y_(w,")
     Call WriteMatrixBlockC3(io_L,3,1,Yx30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YHU24XD","Triplet scale 3","Y_(x,")


     Call WriteMatrixBlockC(io_L,2,MBM3running(2,1:2,1:2),Abs(MWM30(2,2)),"M24B24B" &
           & ,"singlet mass matrix at M_W2(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,2,MGM3running(2,1:2,1:2),Abs(MWM30(2,2)),"M24G24G" &
           & ,"octet mass matrix at M_W2(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,2,MWM3running(2,1:2,1:2),Abs(MWM30(2,2)),"M24W24W" &
           & ,"triplet mass matrix at M_W2(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,2,MXM3running(2,1:2,1:2),Abs(MWM30(2,2)),"M24X24X" &
           & ,"X mass matrix at M_W2(GUT)","M_(W,")
     Call WriteMatrixBlockC3(io_L,3,1,Yb30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YHU24BL","Triplet scale 2","Y_(b,")
     Call WriteMatrixBlockC3(io_L,3,1,Yw30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YHU24WL","Triplet scale 2","Y_(w,")
     Call WriteMatrixBlockC3(io_L,3,1,Yx30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YHU24XD","Triplet scale 2","Y_(x,")

     Call WriteMatrixBlockC(io_L,1,MBM3running(1,1,1),Abs(MWM30(1,1)),"M24B24B" &
           & ,"singlet mass matrix at M_W1(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,1,MGM3running(1,1,1),Abs(MWM30(1,1)),"M24G24G" &
           & ,"octet mass matrix at M_W1(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,1,MWM3running(1,1,1),Abs(MWM30(1,1)),"M24W24W" &
           & ,"triplet mass matrix at M_W1(GUT)","M_(W,")
     Call WriteMatrixBlockC(io_L,1,MXM3running(1,1,1),Abs(MWM30(1,1)),"M24X24X" &
           & ,"X mass matrix at M_W1(GUT)","M_(W,")
     Call WriteMatrixBlockC3(io_L,3,1,Yb30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YHU24BL","Triplet scale 1","Y_(b,")
     Call WriteMatrixBlockC3(io_L,3,1,Yw30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YHU24WL","Triplet scale 1","Y_(w,")
     Call WriteMatrixBlockC3(io_L,3,1,Yx30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YHU24XD","Triplet scale 1","Y_(x,")


    Else If (HighScaleModel.Eq."SEESAW_II_SARAH") Then
     Write(io_L,100) "    3    112   # mSUGRA model + 1 15-plet"
     Write(io_L,100) "Block SEESAWGENERATIONS"
     Write(io_l,98) 15,1,"number of 15-plets"

     Call WriteMatrixBlockC3(io_L,3,2,YT_h15_gut,m_GUT &          ! symmetric
                         & ,"Y15","GUT scale","Y_(15,", .True.)  ! matrix
     Call WriteMatrixBlockC3a(io_L,Lambda1_GUT,m_GUT,"YHD15THD","GUT scale" &
                           & ,"lambda_1")
     Call WriteMatrixBlockC3a(io_L,Lambda2_GUT,m_GUT,"YHU15TBHU","GUT scale" &
                           & ,"lambda_2")
     
     Write(io_L,106) "Block M15 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,107) 1,1,Abs(MTM_gut),"# MWM      "
     Write(io_L,106) "Block M15T15TB Q=",Abs(MTM_GUT),"# (Triplet Scale)"
     Write(io_L,107) 1,1,Abs(MTM),"# M_T   "
     Write(io_L,106) "Block M15Z15ZB Q=",Abs(MTM_GUT),"# (Triplet Scale)"
     Write(io_L,107) 1,1,Abs(MZM),"# M_Z   "
     Write(io_L,106) "Block M15S15SB Q=",Abs(MTM_GUT),"# (Triplet Scale)"
     Write(io_L,107) 1,1,Abs(MSM),"# M_Z   "
     
     Call WriteMatrixBlockC3(io_L,3,2,YT0_h15,Abs(MTM_GUT) &            ! symmetric
                         & ,"YL15TL","triplet scale","Y_(T,", .True.)! matrix
     Call WriteMatrixBlockC3(io_L,3,2,YS0_h15,Abs(MTM_GUT) &            ! symmetric
                         & ,"YD15SD","triplet scale","Y_(S,", .True.)! matrix
     Call WriteMatrixBlockC3(io_L,3,2,YZ0_h15,Abs(MTM_GUT) &
                         & ,"YD15ZL","triplet scale","Y_(Z,")

     Call WriteMatrixBlockC3a(io_L,Lambda10,Abs(MTM_GUT),"YHD15THD","Triplet scale" &
                           & ,"lambda_1")
     Call WriteMatrixBlockC3a(io_L,Lambda20,Abs(MTM_GUT),"YHU15TBHU","Triplet scale" &
                           & ,"lambda_2")

! Florian SARAH, Seesaw II+III      
    End If
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,101) 7,Max(M_SO_10,m_GUT),"# SO(10) scale"
     Write(io_L,101) 8,D_SO_10,"# D-terms at SO(10) scale"
    End If
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"

   !----------------------------------------------------------------
   ! in the case of GenerationMixing Yukawas
   ! are given in the SuperCKM basis 
   !----------------------------------------------------------------
    If (GenerationMixing) Then
     Call Switch_to_superCKM(Y_d_0, Y_u_0, A_d_0, A_u_0, M2_D_0, M2_Q_0, M2_U_0         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .True. &
              &, CKM_out=CKM_Q, Yd=Yd, Yu=Yu )

     Do i1=1,3
      Yl(i1) = Real(Y_l_0(i1,i1),dp)
     End Do


    Else ! .non.GenerationMixing

     Do i1=1,3
      Yu(i1) = Real(Y_u_0(i1,i1),dp)
      Yd(i1) = Real(Y_d_0(i1,i1),dp)
      Yl(i1) = Real(Y_l_0(i1,i1),dp)
     End Do

    End If

    Write(io_L,106) "Block Yu Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yu(1), "# Y_u(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yu(2), "# Y_c(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yu(3), "# Y_t(M_GUT)^DRbar"

    Write(io_L,106) "Block Yd Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yd(1), "# Y_d(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yd(2), "# Y_s(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yd(3), "# Y_b(M_GUT)^DRbar"

    If (GenerationMixing) Call WriteMatrixBlockC(io_L,3,CKM_Q,M_GUT &
                         & ,"VCKM","V_CKM at the GUT scale","V_(")

    Write(io_L,106) "Block Ye Q=",M_GUT,"# (GUT scale)"
    Write(io_L,107) 1,1,Yl(1), "# Y_e(M_GUT)^DRbar"
    Write(io_L,107) 2,2,Yl(2), "# Y_mu(M_GUT)^DRbar"
    Write(io_L,107) 3,3,Yl(3), "# Y_tau(M_GUT)^DRbar"

    If (Model_Suchita) Then
     Call chop(Au_sckm)
     Call chop(Ad_sckm)
     Call WriteMatrixBlockC(io_L,3,Au_sckm,M_GUT,"Tu","M_GUT","T_(u,",tr=.True.)
     Call WriteMatrixBlockC(io_L,3,Ad_sckm,M_GUT,"Td","M_GUT","T_(d,",tr=.True.)
     Call chop(M2Q_sckm)
     Call WriteMatrixBlockC(io_L,3,M2Q_SCKM,M_GUT,"MSQ2" &
           & ,"M^2_Q soft SUSY breaking masses","M^2_(Q,")
     Call chop(M2U_sckm)
     Call WriteMatrixBlockC(io_L,3,M2U_SCKM,M_GUT,"MSU2" &
           & ,"M^2_U soft SUSY breaking masses","M^2_(U,")
     Call chop(M2D_sckm)
     Call WriteMatrixBlockC(io_L,3,M2D_SCKM,M_GUT,"MSD2" &
           & ,"M^2_D soft SUSY breaking masses","M^2_(D,")
    End If ! Model_Suchita

   Else If (HighScaleModel.Eq."GMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    2    # mGMSB model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Lambda,"# Lambda, scale of SUSY breaking"
    Write(io_L,101) 2,MlambdaS,"# M_M, messenger mass scale"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(n5plets+3*n10plets,dp) , &
                      & "# N_5 number of effective five-plets"
    Write(io_L,101) 6,grav_fac,"# c_grav, Gravitino mass factor"

    Write(io_L,106) "Block gauge Q=",MlambdaS,"# (GMSB scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If (HighScaleModel.Eq."AMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    3    # mAMSB model"
    If (GenerationMixing) Write(io_L,100) &
      &             "    6    3    # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,m0_amsb,"# M_0" 
    Write(io_L,101) 2,m_32,"# m_3/2, gravitino mass"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

   Else If (HighScaleModel.Eq."NMSSM") Then
    non_minimal = .True.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    3    1    # NMSSM model"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    If (phase_mu.Ne.0._dp) Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

   Else If (HighScaleModel.Eq."RPexplicit") Then
    non_minimal = .True.
    If (minpar_set) Then
     Write(io_L,100) "Block MODSEL  # Model selection"
     Write(io_L,100) "    1    1    # mSUGRA model"
     Write(io_L,100) "    4    1    # MSSM with explicit R-parity violation"
     Write(io_L,100) "Block MINPAR  # Input parameters"
     Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
     Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
     If (tanb_mZ.Ne.0._dp) Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
     Write(io_L,101) 4, Real(phase_mu,dp),"# cos(phase_mu)"
     Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
     Write(io_L,100) "#"
     If ((Aimag(mu).Ne.0._dp).Or.(Aimag(AoY_l_0(1,1)).Ne.0._dp)) Then
      Write(io_L,100) "Block IMMINPAR  # Input parameters, imaginary part"
      Write(io_L,101) 4, Aimag(phase_mu),"# sin(phase_mu)"
      Write(io_L,101) 5, Aimag(AoY_l_0(1,1)),"# Im(A0)"
     End If

     Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
     Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
     Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
     Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
    Else
     Write(io_L,100) "Block MODSEL  # Model selection"
     Write(io_L,100) "    4    1    # MSSM with explicit R-parity violation"
     Write(io_L,100) "Block MINPAR  # Input parameters"
     If (tanb_mZ.Ne.0._dp) Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
     If (phase_mu.Ne.0._dp) Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    End If

   Else
    non_minimal = .True.
    Write(io_L,100) "# Either the general MSSM or a model has been used"
    Write(io_L,100) &
      & "# which has not yet been implemented in the LesHouches standard"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
   End If 

   If (Sum(in_extpar(:,1)).Gt.0) Then
    Write(io_L,100) "Block EXTPAR  # non-universal input parameters"
    Do i1=1,203
     If (in_extpar(i1,1).Eq.1) &
        & Write(io_L,101) i1,r_extpar(i1) ,"# "//Trim(n_extpar(i1))
    End Do
   End If

   If (Sum(in_extpar(:,2)).Gt.0) Then
    Write(io_L,100) "Block IMEXTPAR  # non-universal input parameters"
    Do i1=1,203
     If (in_extpar(i1,2).Eq.1) &
        & Write(io_L,101) i1,i_extpar(i1) ,"# "//Trim(n_extpar(i1))
    End Do
   End If

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"
   Write(io_L,102) 8,mf_nu(3),"# m_nu_3"
   Write(io_L,102) 11,mf_l(1),"# m_e(pole)"
   Write(io_L,102) 12,mf_nu(1),"# m_nu_1"
   Write(io_L,102) 13,mf_l(2),"# m_muon(pole)"
   Write(io_L,102) 14,mf_nu(2),"# m_nu_2"
   Write(io_L,102) 21,mf_d(1),"# m_d(2 GeV), MSbar"
   Write(io_L,102) 22,mf_u(1),"# m_u(2 GeV), MSbar"
   Write(io_L,102) 23,mf_d(2),"# m_s(2 GeV), MSbar"
   Write(io_L,102) 24,mf_u(2),"# m_c(m_c), MSbar"
   !----------------------------------------------------------------
   ! in the case of GenerationMixing all parameters and mixings
   ! are given in the SuperCKM basis for squarks
   !----------------------------------------------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
    Call Switch_to_superCKM(Y_d, Y_u, A_d, A_u, M2_D, M2_Q, M2_U         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .False. &
              &, RSdown, RSup, Rdsq_ckm, RUsq_ckm, CKM_Q, Yd, Yu )

    If (l_pmns_in) Then
     Write(io_L,100) "Block UPMNSIN  # PMNS matrix, PDG parameterization"
     Write(io_L,102) 1,theta_12,"# theta_12, solar"
     Write(io_L,102) 2,theta_23,"# theta_23, atmospheric"
     Write(io_L,102) 3,theta_13,"# theta_13"
     Write(io_L,102) 4,delta_nu,"# delta_nu"
     Write(io_L,102) 5,alpha_nu1,"# alpha_1"
     Write(io_L,102) 6,alpha_nu2,"# alpha_2"
    End If

    If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
     Call Switch_to_superPMNS(Y_l, MnuL5, A_l, M2_E, M2_L, Al_pmns, M2E_pmns  &
        &, M2L_pmns, .False., RSlepton, RSneut, Rsl_pmns, RSn_pmns, PMNS_Q, Yl)

    Else If (Maxval(Abs(MatNu)).Gt.0._dp) Then
     Call Switch_to_superPMNS(Y_l, MatNu, A_l, M2_E, M2_L, Al_pmns, M2E_pmns  &
        &, M2L_pmns, .False., RSlepton, RSneut, Rsl_pmns, RSn_pmns, PMNS_Q, Yl)

    Else
     Call Switch_to_superPMNS(Y_l, id3C, A_l, M2_E, M2_L, Al_pmns, M2E_pmns  &
        &, M2L_pmns, .False., RSlepton, RSneut, Rsl_pmns, RSn_pmns, PMNS_Q, Yl)

    End If

   Else ! .non.GenerationMixing

    Do i1=1,3
     Yu(i1) = Real(Y_u(i1,i1),dp)
     Yd(i1) = Real(Y_d(i1,i1),dp)
     Yl(i1) = Real(Y_l(i1,i1),dp)
    End Do
    Al_pmns = A_l
    Ad_sckm = A_d
    Au_sckm = A_u

    If (HighScaleModel.Ne."RPexplicit") Then
     M2D_SCKM = M2_D
     M2U_SCKM = M2_U
     M2Q_SCKM = M2_Q
     M2E_pmns = M2_E
     M2L_pmns = M2_L
    End If

    RUsq_ckm = RSup
    RDsq_ckm = RSdown

    RSn_pmns = RSneut
    RSl_pmns = RSlepton

   End If

      
! couplings
  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,gauge(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,gauge(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,gauge(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  Write(io_L,106) "Block Ye Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yl(1), "# Y_e(Q)^DRbar"
  Write(io_L,107) 2,2,Yl(2), "# Y_mu(Q)^DRbar"
  Write(io_L,107) 3,3,Yl(3), "# Y_tau(Q)^DRbar"

  If (GenerationMixing) Then 
                             
   Call WriteMatrixBlockC(io_L,3,CKM_Q,Q &
                         & ,"VCKM","V_CKM at the SUSY scale","V_(")

   Call WriteMatrixBlockC(io_L,3,PMNS_Q,Q &
                         & ,"UPMNS","U_PMNS at the SUSY scale","V_(")

  End If ! generationmixing


  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,Au_sckm,Q,"Tu","SUSY scale","T_(u,",tr=.True.)
   Call WriteMatrixBlockC(io_L,3,Ad_sckm,Q,"Td","SUSY scale","T_(d,",tr=.True.)
   Call WriteMatrixBlockC(io_L,3,Al_pmns,Q,"Te","SUSY scale","T_(l,", tr=.True.)

  Else 
   Write(io_L,106) "Block Au Q=",Q,"# (SUSY scale)"
   If (Abs(y_u(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Au_sckm(1,1)/y_u(1,1),dp), "# A_u(Q)^DRbar"
   If (Abs(y_u(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Au_sckm(2,2)/y_u(2,2),dp), "# A_c(Q)^DRbar"
   If (Abs(y_u(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t(Q)^DRbar"
   If (Maxval(Abs(Aimag(Au_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAu Q=",Q,"# (SUSY scale)"
    If (Abs(y_u(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Au_sckm(1,1)/y_u(1,1)), "# Im(A_u)(Q)^DRbar"
    If (Abs(y_u(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Au_sckm(2,2)/y_u(2,2)), "# Im(A_c)(Q)^DRbar"
    If (Abs(y_u(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Au_sckm(3,3)/y_u(3,3)), "# Im(A_t)(Q)^DRbar"
   End If

   Write(io_L,106) "Block Ad Q=",Q,"# (SUSY scale)"
   If (Abs(y_d(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Ad_sckm(1,1)/y_d(1,1),dp), "# A_d(Q)^DRbar"
   If (Abs(y_d(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Ad_sckm(2,2)/y_d(2,2),dp), "# A_s(Q)^DRbar"
   If (Abs(y_d(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b(Q)^DRbar"
   If (Maxval(Abs(Aimag(Ad_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAd Q=",Q,"# (SUSY scale)"
    If (Abs(y_d(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Ad_sckm(1,1)/y_d(1,1)), "# Im(A_d)(Q)^DRbar"
    If (Abs(y_d(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Ad_sckm(2,2)/y_d(2,2)), "# Im(A_s)(Q)^DRbar"
    If (Abs(y_d(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Ad_sckm(3,3)/y_d(3,3)), "# Im(A_b)(Q)^DRbar"
   End If

   Write(io_L,106) "Block Ae Q=",Q,"# (SUSY scale)"
   If (Abs(Yl(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Al_pmns(1,1)/Yl(1),dp), "# A_e(Q)^DRbar"
   If (Abs(Yl(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Al_pmns(2,2)/Yl(2),dp), "# A_mu(Q)^DRbar"
   If (Abs(Yl(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Al_pmns(3,3)/Yl(3),dp), "# A_tau(Q)^DRbar"
   If (Maxval(Abs(Aimag(Al_pmns))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAe Q=",Q,"# (SUSY scale)"
    If (Abs(Yl(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Al_pmns(1,1)/Yl(1)), "# Im(A_e)(Q)^DRbar"
    If (Abs(Yl(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Al_pmns(2,2)/Yl(2)), "# Im(A_mu)(Q)^DRbar"
    If (Abs(Yl(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Al_pmns(3,3)/Yl(3)), "# Im(A_tau)(Q)^DRbar"
   End If
  End If


  Write(io_L,106) "Block MSOFT Q=",Q,"# soft SUSY breaking masses at Q"
  Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
  Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
  Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
  Write(io_L,104) 21,M2_H(1),"# M^2_(H,d)"
  Write(io_L,104) 22,M2_H(2),"# M^2_(H,u)"

  Write(io_L,104) 31,Sqrt(Real(M2L_pmns(1,1),dp)),"# M_(L,11)"
  Write(io_L,104) 32,Sqrt(Real(M2L_pmns(2,2),dp)),"# M_(L,22)"
  Write(io_L,104) 33,Sqrt(Real(M2L_pmns(3,3),dp)),"# M_(L,33)"
  Write(io_L,104) 34,Sqrt(Real(M2E_pmns(1,1),dp)),"# M_(E,11)"
  Write(io_L,104) 35,Sqrt(Real(M2E_pmns(2,2),dp)),"# M_(E,22)"
  Write(io_L,104) 36,Sqrt(Real(M2E_pmns(3,3),dp)),"# M_(E,33)"
  Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
  Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
  Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
  Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
  Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
  Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
  Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
  Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
  Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"
  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,104) 61,Real(h0,dp),"# lambda"
   Write(io_L,104) 62,0.5_dp*Real(lam,dp),"# kappa"
   Write(io_L,104) 63,Real(ao_h0,dp),"# A_lambda"
   Write(io_L,104) 64,Real(Ao_lam,dp),"# A_kappa"
   Write(io_L,104) 65,Real(oosqrt2*vP*h0,dp),"# mu_eff"
  End If

  If (Maxval(Abs(Aimag(Mi))).Ne.0._dp) Then
   Write(io_L,106) "Block IMMSOFT Q=",Q,"# soft SUSY breaking masses at Q, imaginary parts"
   Write(io_L,104) 1,Aimag(Mi(1)),"# M_1"
   Write(io_L,104) 2,Aimag(Mi(2)),"# M_2"
   Write(io_L,104) 3,Aimag(Mi(3)),"# M_3"
  End If

  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,M2L_pmns,Q,"MSL2" &
           & ,"M^2_L soft SUSY breaking masses","M^2_(L,")

   Call WriteMatrixBlockC(io_L,3,M2E_pmns,Q,"MSE2" &
           & ,"M^2_E soft SUSY breaking masses","M^2_(E,")

   Call WriteMatrixBlockC(io_L,3,M2Q_SCKM,Q,"MSQ2" &
           & ,"M^2_Q soft SUSY breaking masses","M^2_(Q,")

   Call WriteMatrixBlockC(io_L,3,M2U_SCKM,Q,"MSU2" &
           & ,"M^2_U soft SUSY breaking masses","M^2_(U,")

   Call WriteMatrixBlockC(io_L,3,M2D_SCKM,Q,"MSD2" &
           & ,"M^2_D soft SUSY breaking masses","M^2_(D,")

  End If


  If (HighScaleModel.Eq."RPexplicit") Then
   Write(io_L,106) "Block RVKAPPA Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(eps(1),dp),"# epsilon_1"
   Write(io_L,102) 2,Real(eps(2),dp),"# epsilon_2"
   Write(io_L,102) 3,Real(eps(3),dp),"# epsilon_3"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVKAPPA Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(eps(1)),"# Im(epsilon_1)"
    Write(io_L,102) 2,Aimag(eps(2)),"# Im(epsilon_2)"
    Write(io_L,102) 3,Aimag(eps(3)),"# Im(epsilon_3)"
   End If

   Write(io_L,106) "Block RVD Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(Beps(1),dp),"# Re( B_1 epsilon_1)"
   Write(io_L,102) 2,Real(Beps(2),dp),"# Re( B_2 epsilon_2)"
   Write(io_L,102) 3,Real(Beps(3),dp),"# Re( B_3 epsilon_3)"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVD Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(Beps(1)),"# Im( B_1 epsilon_1)"
    Write(io_L,102) 2,Aimag(Beps(2)),"# Im( B_2 epsilon_2)"
    Write(io_L,102) 3,Aimag(Beps(3)),"# Im( B_3 epsilon_3)"
   End If

   Write(io_L,106) "Block RVSNVEV Q=",Q,"# sneutrino vevs at Q"
   Write(io_L,102) 1,vevL(1),"# v_L_1"
   Write(io_L,102) 2,vevL(2),"# v_L_2"
   Write(io_L,102) 3,vevL(3),"# v_L_3"

   If (Maxval(Abs(Rp_lam)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDA Q=",Q,"# lambda_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lam(i1,i2,i3),dp) &
                     & ,"# lambda_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDA Q=",Q,"# Im(lambda_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lam(i1,i2,i3)) &
            & ,"# Im(lambda_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   If (Maxval(Abs(Rp_lamp)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDAP Q=",Q,"# lambda'_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lamp(i1,i2,i3),dp) &
                     & ,"# lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDAP Q=",Q,"# Im(lambda'_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lamp(i1,i2,i3)) &
            & ,"# Im(lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   Write(io_L,100) "Block SPhenoRP  # additional RP parameters"
   Write(io_L,102) 4,Lam_Ex(1),"# Lambda_1 = v_d epsilon_1 + mu v_L1"
   Write(io_L,102) 5,Lam_Ex(2),"# Lambda_2 = v_d epsilon_2 + mu v_L2"
   Write(io_L,102) 6,Lam_Ex(3),"# Lambda_3 = v_d epsilon_3 + mu v_L3"
   Write(io_L,102) 7,(Chi07(3)%m**2-Chi07(2)%m**2)*1.e18_dp,"# m^2_atm [eV^2]"
   Write(io_L,102) 8,(Chi07(2)%m**2-Chi07(1)%m**2)*1.e18_dp,"# m^2_sol [eV^2]"
   Write(io_L,102) 9,Abs(N7(3,6)/N7(3,7))**2,"# tan^2 theta_atm"
   Write(io_L,102) 10,Abs(N7(2,5)/N7(1,5))**2,"# tan^2 theta_sol"
   Write(io_L,102) 11,Abs(N7(3,5))**2,"# U_e3^2"
   Write(io_L,102) 15,vevSM(1),"# v_d"
   Write(io_L,102) 16,vevSM(2),"# v_u"
  End If


   Write(io_L,100) "Block MASS  # Mass spectrum"
   Write(io_L,100) "#   PDG code      mass          particle"
!   Write(io_L,102) id_u(2),mf_u(2),"# m_c(m_c), MSbar"
!   Write(io_L,102) id_d(3),mf_d(3),"# m_b(m_b), MSbar"
   If (MADGraph_style) Write(io_L,102) 5,4.89_dp,"# m_b(pole)"
   Write(io_L,102) id_u(3),mf_u(3),"# m_t(pole)"
   Write(io_L,102) id_Z,mZ,"# m_Z(pole)"
   Write(io_L,102) id_W,mW,"# W+"
   If (HighScaleModel.Ne."RPexplicit") &
        & Write(io_L,102) id_l(3),mf_l(3),"# m_tau(pole)"
   If (HighScaleModel.Eq."NMSSM") Then
    Write(io_L,102) 25,S03(1)%m,"# leightest neutral scalar" 
    Write(io_L,102) 35,S03(2)%m,"# second neutral scalar" 
    Write(io_L,102) 45,S03(3)%m,"# third neutral scalar" 
    Write(io_L,102) 36,P03(2)%m,"# leighter pseudoscalar" 
    Write(io_L,102) 46,P03(3)%m,"# heavier pseudoscalar" 
    Write(io_L,102) 37,Spm(2)%m,"# H+"
   Else If (HighScaleModel.Eq."RPexplicit") Then
!    Write(io_L,102) 12,mN7(1),"# lightest neutrino"
!    Write(io_L,102) 14,mN7(2),"# second lightest neutrino"
!    Write(io_L,102) 16,mN7(3),"# heaviest neutrino"
    If (io_rp.Eq.0) Then  ! the standard SLHA output requires some searching 
     Do i1=1,5
      If (id_s0(i1).Eq.25) Exit
     End Do
     Write(io_L,102) id_s0(i1),S05(i1)%m,"# leightest neutral scalar" 
     Do i1=1,5
      If (id_s0(i1).Eq.35) Exit
     End Do
     Write(io_L,102) id_s0(i1),S05(i1)%m,"# 2nd neutral scalar" 
     Do i1=1,5
      If (id_s0(i1).Eq.1000012) Exit
     End Do
     Write(io_L,102) id_s0(i1),S05(i1)%m,"# 3rd neutral scalar" 
     Do i1=1,5
      If (id_s0(i1).Eq.1000014) Exit
     End Do
     Write(io_L,102) id_s0(i1),S05(i1)%m,"# 4th neutral scalar" 
     Do i1=1,5
      If (id_s0(i1).Eq.1000016) Exit
     End Do
     Write(io_L,102) id_s0(i1),S05(i1)%m,"# 5th neutral scalar" 

     Do i1=1,4
      If (id_p0(i1).Eq.36) Exit
     End Do
     Write(io_L,102) id_P0(i1),P05(i1+1)%m,"# leightest pseudoscalar" 
     Do i1=1,4
      If (id_p0(i1).Eq.1000017) Exit
     End Do
     Write(io_L,102) id_P0(i1),P05(i1+1)%m,"# 2nd pseudoscalar" 
     Do i1=1,4
      If (id_p0(i1).Eq.1000018) Exit
     End Do
     Write(io_L,102) id_P0(i1),P05(i1+1)%m,"# 3rd pseudoscalar" 
     Do i1=1,4
      If (id_p0(i1).Eq.1000019) Exit
     End Do
     Write(io_L,102) id_p0(i1),P05(i1+1)%m,"# 4th pseudoscalar"

     Do i1=1,7
      If (id_Sp(i1).Eq.37) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# leightest charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-1000011) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 2nd charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-1000013) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 3rd charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-1000015) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 4th charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-2000011) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 5th charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-2000013) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 6th charged scalar" 
     Do i1=1,7
      If (id_Sp(i1).Eq.-2000015) Exit
     End Do
     Write(io_L,102) Abs(id_sp(i1)),Spm8(i1+1)%m,"# 7th charged scalar" 

    Else 
     Write(io_L,102) id_s0(1),S05(1)%m,"# leightest neutral scalar" 
     Write(io_L,102) id_s0(2),S05(2)%m,"# 2nd neutral scalar" 
     Write(io_L,102) id_s0(3),S05(3)%m,"# 3rd neutral scalar" 
     Write(io_L,102) id_s0(4),S05(4)%m,"# 4th neutral scalar" 
     Write(io_L,102) id_s0(5),S05(5)%m,"# 5th neutral scalar" 
     Write(io_L,102) id_P0(1),P05(2)%m,"# leightest pseudoscalar" 
     Write(io_L,102) id_P0(2),P05(3)%m,"# 2nd pseudoscalar" 
     Write(io_L,102) id_P0(3),P05(4)%m,"# 3rd pseudoscalar" 
     Write(io_L,102) id_p0(4),P05(5)%m,"# 4th pseudoscalar"
     Write(io_L,102) Abs(id_sp(1)),Spm8(2)%m,"# leightest charged scalar" 
     Write(io_L,102) Abs(id_sp(2)),Spm8(3)%m,"# 2nd charged scalar" 
     Write(io_L,102) Abs(id_sp(3)),Spm8(4)%m,"# 3rd charged scalar" 
     Write(io_L,102) Abs(id_sp(4)),Spm8(5)%m,"# 4th charged scalar" 
     Write(io_L,102) Abs(id_sp(5)),Spm8(6)%m,"# 5th charged scalar" 
     Write(io_L,102) Abs(id_sp(6)),Spm8(7)%m,"# 6th charged scalar" 
     Write(io_L,102) Abs(id_sp(7)),Spm8(8)%m,"# 7th charged scalar" 
    End If
    
   Else
    Write(io_L,102) 25,S0(1)%m,"# h0" 
    Write(io_L,102) 35,S0(2)%m,"# H0" 
    Write(io_L,102) 36,P0(2)%m,"# A0" 
    Write(io_L,102) 37,Spm(2)%m,"# H+"
   End If
! squarks

  If (GenerationMixing) Then
   If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
     i_zaehl = 1
     mat6R = Abs(RDsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_sd(ii) = 1000001
       c_sd(ii) = "~d_L"
      Case(2)
       id_sd(ii) = 1000003
       c_sd(ii) = "~s_L-"
      Case(4)
       id_sd(ii) = 2000001
       c_sd(ii) = "~d_R"
      Case(5)
       id_sd(ii) = 2000003
       c_sd(ii) = "~s_R"
      Case default
       id_sd(ii) = 1000000 * i_zaehl + 5
       c_sd(ii) = "~b_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of sbottoms
      If (id_sd(ii).Eq.1000005)  id_check(1) = ii
      If (id_sd(ii).Eq.2000005)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_sd(ii) = 1000005
      c_sd(ii) = "~b_1"
      ii = id_check(1)  ! the heavier one
      id_sd(ii) = 2000005
      c_sd(ii) = "~b_2"
     End If
     Do ii=1,6
      Write(io_L,102) id_sd(ii),Sdown(ii)%m,"# "//Trim(c_sd(ii))
     End Do

     i_zaehl = 1
     mat6R = Abs(RUsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_su(ii) = 1000002
       c_su(ii) = "~u_L"
      Case(2)
       id_su(ii) = 1000004
       c_su(ii) = "~c_L-"
      Case(4)
       id_su(ii) = 2000002
       c_su(ii) = "~u_R"
      Case(5)
       id_su(ii) = 2000004
       c_su(ii) = "~c_R"
      Case default
       id_su(ii) = 1000000 * i_zaehl + 6
       c_su(ii) = "~t_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of stops
      If (id_su(ii).Eq.1000006)  id_check(1) = ii
      If (id_su(ii).Eq.2000006)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_su(ii) = 1000006
      c_su(ii) = "~t_1"
      ii = id_check(1)  ! the heavier one
      id_su(ii) = 2000006
      c_su(ii) = "~t_2"
     End If

     Do ii=1,6
      Write(io_L,102) id_su(ii),Sup(ii)%m,"# "//Trim(c_su(ii))
     End Do
      
    Else ! use mass ordering

     id_sd(1) = 1000001
     id_sd(2) = 1000003
     id_sd(3) = 1000005
     id_sd(4) = 2000001
     id_sd(5) = 2000003
     id_sd(6) = 2000005
     Do i1=1,6
      c_sd(i1) = "~d_"//Bu(i1)
      Write(io_L,102) id_sd(i1),Sdown(i1)%m,"# "//Trim(c_sd(i1))
     End Do

     id_su(1) = 1000002
     id_su(2) = 1000004
     id_su(3) = 1000006
     id_su(4) = 2000002
     id_su(5) = 2000004
     id_su(6) = 2000006
     Do i1=1,6
      c_su(i1) = "~u_"//Bu(i1)
      Write(io_L,102) id_su(i1),Sup(i1)%m,"# "//Trim(c_su(i1))
     End Do
    End If
  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,Sdown(1)%m,"# ~d_L"
    Write(io_L,102) 2000001,Sdown(2)%m,"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,Sdown(2)%m,"# ~d_L"
    Write(io_L,102) 2000001,Sdown(1)%m,"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,Sup(1)%m,"# ~u_L"
    Write(io_L,102) 2000002,Sup(2)%m,"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,Sup(2)%m,"# ~u_L"
    Write(io_L,102) 2000002,Sup(1)%m,"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,Sdown(3)%m,"# ~s_L"
    Write(io_L,102) 2000003,Sdown(4)%m,"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,Sdown(4)%m,"# ~s_L"
    Write(io_L,102) 2000003,Sdown(3)%m,"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,Sup(3)%m,"# ~c_L"
    Write(io_L,102) 2000004,Sup(4)%m,"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,Sup(4)%m,"# ~c_L"
    Write(io_L,102) 2000004,Sup(3)%m,"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,Sdown(5)%m,"# ~b_1"
   Write(io_L,102) 2000005,Sdown(6)%m,"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,Sup(5)%m,"# ~t_1"
   Write(io_L,102) 2000006,Sup(6)%m,"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
   
   If (HighScaleModel.Ne."RPexplicit") Then

! sleptons
    If (GenerationMixing) Then
     If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
 
      mat3R = Abs(RSn_pmns)
      Do i1=1,3
       MaxCont = Maxval(mat3R)
       Call FindPosition(3, mat3R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_snu(ii) = 1000012
        c_snu(ii) = "~nu_eL"
       Case(2)
        id_snu(ii) = 1000014
        c_snu(ii) = "~nu_muL"
       Case(3)
        id_snu(ii) = 1000016
        c_snu(ii) = "~nu_tauL"
       End Select
       mat3R(ii,:) = 0._dp
       mat3R(:,jj) = 0._dp
      End Do
      Do ii=1,3
       Write(io_L,102) id_snu(ii),Sneut(ii)%m,"# "//Trim(c_snu(ii))
      End Do

      i_zaehl = 1
      mat6R = Abs(RSl_pmns)
      Do i1=1,6
       MaxCont = Maxval(mat6R)
       Call FindPosition(6, mat6R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_sle(ii) = 1000011
        c_sle(ii) = "~e_L-"
        c_slep(ii) = "~e_L+"
       Case(2)
        id_sle(ii) = 1000013
        c_sle(ii) = "~mu_L-"
        c_slep(ii) = "~mu_L+"
       Case(4)
        id_sle(ii) = 2000011
        c_sle(ii) = "~e_R-"
        c_slep(ii) = "~e_R+"
       Case(5)
        id_sle(ii) = 2000013
        c_sle(ii) = "~mu_R-"
        c_slep(ii) = "~mu_R+"
       Case default
        id_sle(ii) = 1000000 * i_zaehl + 15
        c_sle(ii) = "~tau_"//bu(i_zaehl)//"-"
        c_slep(ii) = "~tau_"//bu(i_zaehl)//"+"
        i_zaehl = I_zaehl + 1
       End Select
       mat6R(ii,:) = 0._dp
       mat6R(:,jj) = 0._dp
      End Do
      Do ii=1,6 ! check ordering of staus
       If (id_sle(ii).Eq.1000015)  id_check(1) = ii
       If (id_sle(ii).Eq.2000015)  id_check(2) = ii
      End Do
      If (id_check(1).Gt.id_check(2)) Then ! switch ordering
       ii = id_check(2)  ! the lighter one
       id_sle(ii) = 1000015
       c_sle(ii) = "~tau_1-"
       c_slep(ii) = "~tau_1+"
       ii = id_check(1)  ! the heavier one
       id_sle(ii) = 2000015
       c_sle(ii) = "~tau_2-"
       c_slep(ii) = "~tau_2+"
      End If

      Do ii=1,6
       Write(io_L,102) id_sle(ii),Slepton(ii)%m,"# "//Trim(c_sle(ii))
      End Do

     Else ! using mass ordering

      id_snu(1) = 1000012
      id_snu(2) = 1000014
      id_snu(3) = 1000016
      Do i1=1,3
       c_snu(i1) = "~nu_"//Bu(i1)
       Write(io_L,102) id_snu(i1),Sneut(i1)%m,"# "//Trim(c_snu(i1))
      End Do

      id_sle(1) = 1000011
      id_sle(2) = 1000013
      id_sle(3) = 1000015
      id_sle(4) = 2000011
      id_sle(5) = 2000013
      id_sle(6) = 2000015
      Do i1=1,6
       c_sle(i1) = "~l_"//Bu(i1)
       c_slep(i1) = "~l+_"//Bu(i1)
       Write(io_L,102) id_sle(i1),Slepton(i1)%m,"# "//Trim(c_sle(i1))
      End Do
     End If

    Else ! .not.GenerationMixing

     id_snu = (/ 1000012, 1000014, 1000016 /)
     If (Abs(RSl_pmns(1,1)).Gt.0.5_dp) Then
      Write(io_L,102) 1000011,slepton(1)%m,"# ~e_L-"
      Write(io_L,102) 2000011,slepton(2)%m,"# ~e_R-"
      id_sle(1) = 1000011
      id_sle(2) = 2000011
      c_sle(1) = "~e_L-"
      c_sle(2) = "~e_R-"
      c_slep(1) = "~e_L+"
      c_slep(2) = "~e_R+"
     Else
      Write(io_L,102) 1000011,slepton(2)%m,"# ~e_L-"
      Write(io_L,102) 2000011,slepton(1)%m,"# ~e_R-"
      id_sle(2) = 1000011
      id_sle(1) = 2000011
      !------------------------
      ! change ordering
      !------------------------
      c_sle(2) = "~e_L-"
      c_sle(1) = "~e_R-"
      c_slep(2) = "~e_L+"
      c_slep(1) = "~e_R+"
     End If
     Write(io_L,102) 1000012,Sneut(1)%m,"# ~nu_eL"
     c_snu(1) = "~nu_eL"
     If (Abs(RSl_pmns(3,3)).Gt.0.5_dp) Then
      Write(io_L,102) 1000013,slepton(3)%m,"# ~mu_L-"
      Write(io_L,102) 2000013,slepton(4)%m,"# ~mu_R-"
      id_sle(3) = 1000013
      id_sle(4) = 2000013
      c_sle(3) = "~mu_L-"
      c_sle(4) = "~mu_R-"
      c_slep(3) = "~mu_L+"
      c_slep(4) = "~mu_R+"
     Else
      Write(io_L,102) 1000013,slepton(4)%m,"# ~mu_L-"
      Write(io_L,102) 2000013,slepton(3)%m,"# ~mu_R-"
      id_sle(4) = 1000013
      id_sle(3) = 2000013
      !------------------------
      ! change ordering
      !------------------------
      c_sle(4) = "~mu_L-"
      c_sle(3) = "~mu_R-"
      c_slep(4) = "~mu_L+"
      c_slep(3) = "~mu_R+"
     End If
     Write(io_L,102) 1000014,Sneut(2)%m,"# ~nu_muL"
     c_snu(2) = "~nu_muL"
     Write(io_L,102) 1000015,slepton(5)%m,"# ~tau_1-"
     Write(io_L,102) 2000015,slepton(6)%m,"# ~tau_2-"
     id_sle(5) = 1000015
     id_sle(6) = 2000015
     c_sle(5) = "~tau_1-"
     c_sle(6) = "~tau_2-"
     c_slep(5) = "~tau_1+"
     c_slep(6) = "~tau_2+"
     Write(io_L,102) 1000016,Sneut(3)%m,"# ~nu_tauL"
     c_snu(3) = "~nu_tauL"
    End If
   End If ! GenerationMixing

 ! gauginos/higgsinos
   Write(io_L,102) 1000021,Glu%m,"# ~g"
   ! checking for negative sign
   nr = ZeroC
   If (HighScaleModel.Eq."NMSSM") Then
    Do i1=1,5
     If (Sum(Abs(Real(N5(i1,:)))).Eq.0._dp) Then
      mNr(i1) = - Chi05(i1)%m
      nr(i1,1:5) = (0._dp,-1._dp) * n5(i1,:)
     Else   
      mNr(i1) =  Chi05(i1)%m
      nr(i1,1:5) = n5(i1,:)
     End If
    End Do

   Else If (HighScaleModel.Eq."RPexplicit") Then

    If (l_RP_Pythia) Then  ! Pythia only takes 4x4 matrix for neutralinos
                           ! and 2x2 for charginos
     mNr(1:7) = Abs(Chi07(1:7)%m)
     Do i1=1,4
      If (Sum(Abs(Real(N(i1,:)))).Lt.0.1_dp) Then
       nr(i1,1:4) = Aimag(n(i1,:))
      Else   
       nr(i1,1:4) = n(i1,:)
      End If
     End Do

     U = U5(4:5,1:2)
     V = V5(4:5,1:2)

    Else
     Do i1=1,7
      If (Sum(Abs(Real(N7(i1,:)))).Eq.0._dp) Then
       mNr(i1) = - Chi07(i1)%m
       nr(i1,1:7) = Aimag(n7(i1,:))
      Else   
       mNr(i1) =  Chi07(i1)%m
       nr(i1,1:7) = n7(i1,:)
      End If
     End Do
    End If

   Else
    Do i1=1,4
     If (Sum(Abs(Real(N(i1,:)))).Eq.0._dp) Then
      mNr(i1) = - Chi0(i1)%m
      nr(i1,1:4) = Aimag(n(i1,:))
     Else   
      mNr(i1) =  Chi0(i1)%m
      nr(i1,1:4) = n(i1,:)
     End If
    End Do
   End If

   If (HighScaleModel.Eq."RPexplicit") Then
    Do i1=1,7
     Write(io_L,102) id_n(i1),mnr(i1),"# "//Trim(c_c0(i1))
    End Do
    Do i1=1,5
     Write(io_L,102) id_c(i1),ChiPm5(i1)%m,"# "//Trim(c_cp(i1))
    End Do

   Else
    Write(io_L,102) 1000022,mnr(1),"# ~chi_10" 
    Write(io_L,102) 1000023,mnr(2),"# ~chi_20" 
    Write(io_L,102) 1000025,mnr(3),"# ~chi_30" 
    Write(io_L,102) 1000035,mnr(4),"# ~chi_40" 
    If (HighScaleModel.Eq."NMSSM") Write(io_L,102) 1000045,mnr(5),"# ~chi_50"
    Write(io_L,102) 1000024,ChiPm(1)%m,"# ~chi_1+" 
    Write(io_L,102) 1000037,ChiPm(2)%m,"# ~chi_2+"
   End If

   If ((HighScaleModel.Eq."GMSB").or.l_m32_in) &
     & Write(io_L,102) 1000039, m32,"# ~G"

   If (Maxval(MnuR).Gt.0._dp) Then
    Write(io_L,100) "# masses of right handed neutrinos"
    Write(io_L,100) "Block MnuRnuR"
    Write(io_L,107) 1,1,mnur(1),"# m_nu_R_1"
    Write(io_L,107) 2,2,mnur(2),"# m_nu_R_2"
    Write(io_L,107) 3,3,mnur(3),"# m_nu_R_3"
   End If

! Mixing matrices 
  Write(io_L,100) "# Higgs mixing"
  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,103) "Block HIGMIX Q=",Q, "# neutral scalar Higgs mixing"
   If (.Not.External_Higgs) Then
    RS03_save = RS03
    RG0 = 0._dp
    RG0(1,1) = - 1._dp / Sqrt(1._dp + tanb**2)
    RG0(2,2) = - RG0(1,1)
    RG0(3,3) = 1._dp
    RG0(1,2) = tanb * RG0(1,1)
    RG0(2,1) = RG0(1,2)
    RSpm = RG0(1:2,1:2)
!    RP03_save = Matmul(Transpose(RG0),RP03)
   End If

   Do i1=1,n_s0
    Do i2=1,n_s0
     Write(io_L,105) i1,i2,RS03_save(i1,i2),"# R_S0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   Write(io_L,103) "Block AMIX Q=",Q, "# pseudoscalar Higgs mixing"
   Do i1=1,n_p0+1
    Do i2=1,n_p0+1
     Write(io_L,105) i1,i2,RP03_save(i1,i2) &
             & ,"# R_P0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
  Else If (HighScaleModel.Eq."RPexplicit") Then
   If (l_RP_Pythia) Then
    Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
    Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
    Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
    Write(io_L,104) 1,Real(mu,dp),"# mu"
    Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
    Write(io_L,104) 3,vev_Q,"# v(Q)"
    Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"
    If (Aimag(mu).Ne.0._dp) Then
     Write(io_L,103) "Block IMHmix Q=",Q, "# Higgs mixing parameters"
     Write(io_L,104) 1,Aimag(mu),"# Im(mu)"
    End If
    Write(io_L,100) "Block staumix  # stau mixing matrix"
    Write(io_L,105) 1,1,Real(RSl_pmns(5,5),dp),"# R_sta(1,1)"
    Write(io_L,105) 1,2,Real(RSl_pmns(5,6),dp),"# R_sta(1,2)"
    Write(io_L,105) 2,1,Real(RSl_pmns(6,5),dp),"# R_sta(2,1)"
    Write(io_L,105) 2,2,Real(RSl_pmns(6,6),dp),"# R_sta(2,2)"

   Else

   Write(io_L,103) "Block RVHMIX  Q=",Q, "# neutral scalar Higgs mixing"

   Do i1=1,n_s0
    Do i2=1,n_s0
     Write(io_L,105) i1,i2,RS05(i1,i2),"# R_S0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   Write(io_L,103) "Block RVAMIX  Q=",Q, "# pseudoscalar Higgs mixing"
   Do i1=1,n_p0+1
    Do i2=1,n_p0+1
     Write(io_L,105) i1,i2,RP05(i1,i2) &
             & ,"# R_P0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do 
   Write(io_L,103) "Block RVLMIX Q=",Q, "# charged Higgs mixing"
   Do i1=1,8
    Do i2=1,8
     Write(io_L,105) i1,i2,Real(Rspm8(i1,i2),dp) &
             & ,"# R_Spm("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   End If

  Else
   Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
   Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
   Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
   Write(io_L,104) 1,Real(mu,dp),"# mu"
   Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
   Write(io_L,104) 3,vev_Q,"# v(Q)"
   Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"
   If (Aimag(mu).Ne.0._dp) Then
    Write(io_L,103) "Block IMHmix Q=",Q, "# Higgs mixing parameters"
    Write(io_L,104) 1,Aimag(mu),"# Im(mu)"
   End If
  End If
! for MADGraph use
  If (MADGraph_style) Then
   Write(io_L,100) "Block fralpha # Effective Higgs mixing angle"
   Write(io_L,104) 1,-Asin(RS0(1,1)),"# alpha"
  End If

  If (generationmixing) Then
   Call WriteMatrixBlockC2(io_L, 6, RUsq_ckm, "USQmix" &
                         &, "u-sqark mixing matrix", "R_Su(")

   Call WriteMatrixBlockC2(io_L, 6, RDsq_ckm, "DSQmix" &
                         &, "d-sqark mixing matrix", "R_Sd(")

   If (HighScaleModel.Ne."RPexplicit") Then

    Call WriteMatrixBlockC2(io_L, 6, RSl_pmns, "SELmix" &
                         &, "slepton mixing matrix", "R_Sl(")

    Call WriteMatrixBlockC2(io_L, 3, RSn_pmns, "SNUmix" &
                         &, "sneutrino mixing matrix", "R_Sn(")

   End If

  Else ! .not.GenerationMixing

   Call WriteMatrixBlockC2(io_L, 2, RUsq_ckm(5:6,5:6), "stopmix" &
                         &, "stop mixing matrix", "R_st(")

   Call WriteMatrixBlockC2(io_L, 2, RDsq_ckm(5:6,5:6), "sbotmix" &
                         &, "sbottom mixing matrix", "R_sb(")

   If (HighScaleModel.Ne."RPexplicit") & 
    & Call WriteMatrixBlockC2(io_L, 2, RSl_pmns(5:6,5:6), "staumix" &
                            &, "stau mixing matrix", "R_sta(")
  End If

  If ((HighScaleModel.Eq."RPexplicit").And.(.Not.l_RP_Pythia)) Then
   Call WriteMatrixBlockC2(io_L, n_n, nr(1:n_n,1:n_n) &
                         &, "RVNmix", "/neutrino/neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, 5, U5, "RVUmix" &
                         &, "lepton/chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, 5, V5, "RVVmix" &
                         &, "lepton/chargino mixing matrix", "V(")
  Else If (HighScaleModel.Eq."RPexplicit")  Then
   Call WriteMatrixBlockC2(io_L, 4, nr(1:4,1:4), "Nmix" &
                         &, "neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, 2, U, "Umix" &
                         &, "chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, 2, V, "Vmix" &
                         &, "chargino mixing matrix", "V(")
  Else   
   Call WriteMatrixBlockC2(io_L, n_n, nr(1:n_n,1:n_n), "Nmix" &
                         &, "neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, n_c, U, "Umix" &
                         &, "chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, n_c, V, "Vmix" &
                         &, "chargino mixing matrix", "V(")
  End If

 !-------------------------------------------------------------------------
 ! before writing the decays and other information, the above information
 ! is given at different scales, if required
 !-------------------------------------------------------------------------
 If (l_Q_out) then ! need parameters and potentially also the masses at another scale
  If (Qout.eq.0._dp) Qout=Q 
  delta_Q = Exp( Log(Qout/mZ) / n_Q_out )
  Q_act = mZ
  Do i1=0,n_Q_out
   if (i1.gt.0) Q_act = Q_act * delta_Q
   Write(io_L,100) "###########################################################"
   If (Calc_Mass) then
    Write(io_L,301) Q_act
   Else
    Write(io_L,300) Q_act
   End If
   Write(io_L,100) "###########################################################"
   tb_save = tanb_Q  ! unfortunately, this value gets changed in the
                     ! current implementation and needs to be re-adjusted
                     ! afterwards. Requires a major modification of the code
                     ! to avoid this trick here
   Call WriteParametersAtQ(Q, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, tanb_Q, phase_mu      &
          & , GenerationMixing, Q_act, 1.e-4_dp, io_L)
   tanb_Q = tb_save 
   Write(io_L,100) "###########################################################"
   If (Calc_Mass) then
    Write(io_L,311) Q_act
   Else
    Write(io_L,310) Q_act
   End If
   Write(io_L,100) "###########################################################"
  End Do
 End If

 Do i1=1,36 ! at the pol-masses of the differnt particles
  Q_PDG_out = abs(Q_PDG_out) ! in case that a fermion mass is negative
  If (l_PDG_out(i1)) Then
   Write(io_L,100) "###########################################################"
   If (Calc_Mass) then
    Write(io_L,301) Q_PDG_out(i1)
   Else
    Write(io_L,300) Q_PDG_out(i1)
   End If
   Write(io_L,100) "###########################################################"
   tb_save = tanb_Q  ! unfortunately, this value gets changed in the
                     ! current implementation and needs to be re-adjusted
                     ! afterwards. Requires a major modification of the code
                     ! to avoid this trick here
   Call WriteParametersAtQ(Q, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, tanb_Q, phase_mu      &
          & , GenerationMixing, Q_PDG_out(i1), 1.e-4_dp, io_L)
   tanb_Q = tb_save 
   Write(io_L,100) "###########################################################"
   If (Calc_Mass) then
    Write(io_L,311) Q_PDG_out(i1)
   Else
    Write(io_L,310) Q_PDG_out(i1)
   End If
   Write(io_L,100) "###########################################################"
  End If
 End Do
 !------------------
 ! branching ratios
 !------------------
 If (L_BR) Then
! for MADGraph use
!  If (MADGraph_style) Then
   Write(io_L,200) id_Z,2.4952_dp,"Z"
   Write(io_L,200) id_W,2.085_dp,"W"
!  End If
  !---------------------------------
  ! sleptons
  !---------------------------------
  Do i1=1,n_sl
   p_id = Slepton(i1)%id
   Write(io_L,200) id_p(p_id),Slepton(i1)%g,Trim(names(p_id))
   If (Slepton(i1)%g.Gt.0._dp) Then 
    If (Sum(Slepton(i1)%gi2).Gt.0._dp) Then ! 2-body decays
     Write(io_L,100) "#    BR                NDA      ID1      ID2"
     Do i2=1,200
      p_id1 = Slepton(i1)%id2(i2,1) 
      p_id2 = Slepton(i1)%id2(i2,2) 
      If (Slepton(i1)%bi2(i2).Gt.BrMin) Then
       Write(io_L,401) Slepton(i1)%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
      End If
     End Do
    End If

    If (Sum(Slepton(i1)%gi3).Gt.0._dp) Then ! 3-body decays
     Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
     Do i2=1,600
      p_id1 =Slepton(i1)%id3(i2,1) 
      p_id2 =Slepton(i1)%id3(i2,2) 
      p_id3 =Slepton(i1)%id3(i2,3) 
      If (Slepton(i1)%bi3(i2).Gt.BrMin) Then
       Write(io_L,402) Slepton(i1)%bi3(i2),3,id_p(p_id1),id_p(p_id2)      &
                    & ,id_p(p_id3), Trim(names(p_id)), Trim(names(p_id1)) &
                    & , Trim(names(p_id2)), Trim(names(p_id3))
      End If
     End Do
    End If

   End If
  End Do   

  !---------------------------------
  ! sneutrinos
  !---------------------------------
  Do i1=1,n_sn
   p_id = Sneut(i1)%id
   Write(io_L,200) id_p(p_id),Sneut(i1)%g,Trim(names(p_id))
   If (Sneut(i1)%g.Gt.0._dp) Then ! 2-body decays
    If (Sum(Sneut(i1)%gi2).Gt.0._dp) Then ! 2-body decays
     Write(io_L,100) "#    BR                NDA      ID1      ID2"
     Do i2=1,200
      p_id1 = Sneut(i1)%id2(i2,1) 
      p_id2 = Sneut(i1)%id2(i2,2) 
      If (Sneut(i1)%bi2(i2).Gt.BrMin) Then
       Write(io_L,401) Sneut(i1)%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                   & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
      End If
     End Do
    End If

    If (Sum(Sneut(i1)%gi3).Gt.0._dp) Then ! 3-body decays
     Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
     Do i2=1,600
      p_id1 =Sneut(i1)%id3(i2,1) 
      p_id2 =Sneut(i1)%id3(i2,2) 
      p_id3 =Sneut(i1)%id3(i2,3) 
      If (Sneut(i1)%bi3(i2).Gt.BrMin) Then
       Write(io_L,402) Sneut(i1)%bi3(i2),3,id_p(p_id1),id_p(p_id2)      &
                     & ,id_p(p_id3), Trim(names(p_id)), Trim(names(p_id1)) &
                     & , Trim(names(p_id2)), Trim(names(p_id3))
      End If
     End Do
    End If

   End If
  End Do   

  !---------------------------------
  ! d-squarks
  !---------------------------------
  Do i1=1,n_sd
   p_id = Sdown(i1)%id
   Write(io_L,200) id_p(p_id),Sdown(i1)%g,Trim(names(p_id))
   If (Sdown(i1)%g.Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 = Sdown(i1)%id2(i2,1) 
     p_id2 = Sdown(i1)%id2(i2,2) 
     If (Sdown(i1)%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) Sdown(i1)%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
  End Do   

  !---------------------------------
  ! u-squarks
  !---------------------------------
  Do i1=1,n_su
   p_id = Sup(i1)%id
   Write(io_L,200) id_p(p_id),Sup(i1)%g,Trim(names(p_id))
   If (Sup(i1)%g.Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 = Sup(i1)%id2(i2,1) 
     p_id2 = Sup(i1)%id2(i2,2) 
     If (Sup(i1)%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) Sup(i1)%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
    If (Sum(Sup(i1)%gi3).Gt.0._dp) Then ! 3-body decays
     Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
     Do i2=1,600
      p_id1 =Sup(i1)%id3(i2,1) 
      p_id2 =Sup(i1)%id3(i2,2) 
      p_id3 =Sup(i1)%id3(i2,3) 
      If (Sup(i1)%bi3(i2).Gt.BrMin) Then
       Write(io_L,402) Sup(i1)%bi3(i2),3,id_p(p_id1),id_p(p_id2),id_p(p_id3) &
                    & , Trim(names(p_id)), Trim(names(p_id1))                 &
                    & , Trim(names(p_id2)), Trim(names(p_id3))
      End If
     End Do
    End If
   End If

  End Do   ! loop over u-squarks

  !--------------
  ! charginos
  !--------------
  If (HighScaleModel.Eq."RPexplicit") Then
    i1_min = 4
  Else
    i1_min = 1
  End If

  Do i1=i1_min,n_c
   If (HighScaleModel.Eq."RPexplicit") Then
    part23 = ChiPm5(i1)

   Else   ! MSSM output
    part23 = ChiPm(i1)

   End If ! model selection 

   p_id = part23%id
   Write(io_L,200) id_p(p_id),part23%g,Trim(names(p_id))
   If (Sum(part23%gi2).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 =part23%id2(i2,1) 
     p_id2 =part23%id2(i2,2) 
     If (part23%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) part23%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
   If (Sum(part23%gi3).Gt.0._dp) Then ! 3-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,600
     p_id1 =part23%id3(i2,1) 
     p_id2 =part23%id3(i2,2) 
     p_id3 =part23%id3(i2,3) 
     If (part23%bi3(i2).Gt.BrMin) Then
      Write(io_L,402) part23%bi3(i2),3,id_p(p_id1),id_p(p_id2),id_p(p_id3) &
                    & , Trim(names(p_id)), Trim(names(p_id1))                 &
                    & , Trim(names(p_id2)), Trim(names(p_id3))
     End If
    End Do
   End If
  End Do ! i1


  !--------------
  ! neutralinos
  !--------------
  If (HighScaleModel.Eq."RPexplicit") Then
    i1_min = 4
  Else
    i1_min = 1
  End If

  Do i1=i1_min,n_n
   If (HighScaleModel.Eq."NMSSM") Then
    part23 = Chi05(i1)

   Else If (HighScaleModel.Eq."RPexplicit") Then
    part23 = Chi07(i1)

   Else   ! MSSM output
    part23 = Chi0(i1)

   End If ! model selection 

   p_id = part23%id
   Write(io_L,200) id_p(p_id),part23%g,Trim(names(p_id))
   If (Sum(part23%gi2).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 =part23%id2(i2,1) 
     p_id2 =part23%id2(i2,2) 
     If (part23%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) part23%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                   & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
   If (Sum(part23%gi3).Gt.0._dp) Then ! 3-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,600
     p_id1 =part23%id3(i2,1) 
     p_id2 =part23%id3(i2,2) 
     p_id3 =part23%id3(i2,3) 
     If (part23%bi3(i2).Gt.BrMin) Then
!      If ((id_p(p_id2).Eq.12).And.(id_p(p_id2).Eq.-12)) 
      Write(io_L,402) part23%bi3(i2),3,id_p(p_id1),id_p(p_id2),id_p(p_id3) &
                   & , Trim(names(p_id)), Trim(names(p_id1))               &
                   & , Trim(names(p_id2)), Trim(names(p_id3))
     End If
    End Do
   End If
  End Do ! loop over neutralinos

  !--------------
  ! gluino
  !--------------
   p_id = Glu%id
   Write(io_L,200) id_p(p_id),Glu%g,Trim(names(p_id))
   If (Sum(Glu%gi2).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 =Glu%id2(i2,1) 
     p_id2 =Glu%id2(i2,2) 
     If (Glu%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) Glu%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
   If (Sum(Glu%gi3).Gt.0._dp) Then ! 3-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,600
     p_id1 =Glu%id3(i2,1) 
     p_id2 =Glu%id3(i2,2) 
     p_id3 =Glu%id3(i2,3) 
     If (Glu%bi3(i2).Gt.BrMin) Then
      Write(io_L,402) Glu%bi3(i2),3,id_p(p_id1),id_p(p_id2),id_p(p_id3)    &
                 & ,Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2)) &
                 & ,Trim(names(p_id3))
     End If
    End Do
   End If

  !-----------------
  ! neutral scalars
  !-----------------
  Do i1=1,n_s0
   If ( HighScaleModel.Eq."NMSSM") Then
    part23 = S03(i1)
   Else If ( HighScaleModel.Eq."RPexplicit") Then
    part23 = S05(i1)
   Else ! MSSM
    part23 = S0(i1)
   End If

   p_id = part23%id
   Write(io_L,200) id_p(p_id),part23%g,Trim(names(p_id))
   If (part23%g.Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 = part23%id2(i2,1) 
     p_id2 = part23%id2(i2,2) 
     If (part23%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) part23%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                  & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
   If(BR_Higgs_with_offshell_V) Then
    Do i2=1,10
     p_id1 = part23%id3(i2,1) 
     p_id2 = part23%id3(i2,2)
     If (part23%bi3(i2).Gt.BrMin) Then
      Write(io_L,401) part23%bi3(i2),2,id_p(p_id1),id_p(p_id2), &
              & Trim(names(p_id)),Trim(names(p_id1)),"("//Trim(names(p_id2))//")^* "
     End If
    End Do
   Else
  !---------------------------------------------------------------------------
  ! fold with branching ratios of Z and W
  ! the following needs to be changed in case of extended gauge symmetries
  !---------------------------------------------------------------------------
    If (Sum(part23%gi3).Gt.0._dp) Then
     Write(io_L,100) "# writing decays into V V* as 3-body decays"
     Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
     Do i2=1,10
      p_id1 = part23%id3(i2,1) 
      p_id2 = part23%id3(i2,2)
      If (p_id1.Eq.0) Cycle
      If ( id_p(p_id1).Eq.24) Then ! positive charged W-boson is on-shell
       Do i3=1,3
        Write(io_L,402) BrWln(i3) * part23%bi3(i2),3,id_p(p_id1),id_lept(i3) &
                  & ,id_neut(i3), Trim(names(p_id)),Trim(names(p_id1))       &
                  & ,Trim(names(10+2*i3)), Trim(names(4+2*i3))
       End Do
       Do i3=1,2
        Write(io_L,402) BrWqq(i3) * part23%bi3(i2),3,id_p(p_id1),id_d_quark(i3) &
                  & , -id_u_quark(i3), Trim(names(p_id)),Trim(names(p_id1))     &
                  & , Trim(names(22+2*i3)), Trim(names(17+2*i3))
       End Do
      Else If ( id_p(p_id1).Eq.-24) Then ! negative charged W-boson is on-shell
       Do i3=1,3
        Write(io_L,402) BrWln(i3) * part23%bi3(i2),3,id_p(p_id1),-id_lept(i3) &
                  & ,-id_neut(i3), Trim(names(p_id)),Trim(names(p_id1))       &
                  & ,Trim(names(11+2*i3)), Trim(names(5+2*i3))
       End Do
       Do i3=1,2
        Write(io_L,402) BrWqq(i3) * part23%bi3(i2),3,id_p(p_id1),-id_d_quark(i3)&
                   & , id_u_quark(i3), Trim(names(p_id)),Trim(names(p_id1))     &
                   & , Trim(names(23+2*i3)), Trim(names(16+2*i3))
       End Do
      Else ! the Z-boson
       Do i3=1,3
        Write(io_L,402) BrZll(i3)*part23%bi3(i2),3,id_p(p_id1),id_lept(i3)  &
                   & ,-id_lept(i3), Trim(names(p_id)),Trim(names(p_id1))    &
                   & ,Trim(names(10+2*i3)), Trim(names(11+2*i3))
       End Do
       Write(io_L,402) BrZinv*part23%bi3(i2),3,id_p(p_id1),id_neut(1),-id_neut(1)&
                  & , Trim(names(p_id)),Trim(names(p_id1)),Trim(names(6))        &
                  & , Trim(names(7))
       Do i3=1,3
        Write(io_L,402) BrZqq(i3)*part23%bi3(i2),3,id_p(p_id1),id_d_quark(i3) &
                   & , -id_d_quark(i3), Trim(names(p_id)),Trim(names(p_id1))  &
                   & , Trim(names(22+2*i3)), Trim(names(23+2*i3))
       End Do
       Do i3=1,2
        Write(io_L,402) BrZqq(i3+3)*part23%bi3(i2),3,id_p(p_id1),id_u_quark(i3) &
                   & , -id_u_quark(i3), Trim(names(p_id)),Trim(names(p_id1))    &
                   & , Trim(names(16+2*i3)), Trim(names(17+2*i3))
       End Do
      End If

     End Do

    End If ! check of 3-body decays contribute

   End If ! if h->V V* decays is folded with BRs of V*

  End Do ! loop over scalars

  !----------------------
  ! neutral pseudoscalar
  !----------------------
  Do i1=2,n_p0+1
   If ( HighScaleModel.Eq."NMSSM") Then
    part2 = P03(i1)
   Else If ( HighScaleModel.Eq."RPexplicit") Then
    part2 = P05(i1)
   Else
    part2 = P0(i1)
   End If ! model selection 

   p_id = part2%id
   Write(io_L,200) id_p(p_id),part2%g,Trim(names(p_id))
   If (part2%g.Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 = part2%id2(i2,1) 
     p_id2 = part2%id2(i2,2) 
     If (part2%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) part2%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                 & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If

  End Do ! loop over pseudoscalars


  !-----------------
  ! charged scalars
  !-----------------
  Do i1=2,n_spm+1
   If ( HighScaleModel.Eq."RPexplicit") Then
    part2 = Spm8(i1)
   Else
    part2 = Spm(i1)
   End If

   p_id = part2%id
   Write(io_L,200) id_p(p_id),part2%g,Trim(names(p_id))
   If (part2%g.Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,200
     p_id1 = part2%id2(i2,1) 
     p_id2 = part2%id2(i2,2) 
     If (part2%bi2(i2).Gt.BrMin) Then
      Write(io_L,401) part2%bi2(i2),2,id_p(p_id1),id_p(p_id2), &
                   & Trim(names(p_id)),Trim(names(p_id1)),Trim(names(p_id2))
     End If
    End Do
   End If
  End Do ! loop of charged scalars

 End If ! L_BR  

 !---------------------------------------------
 ! cross section information
 !---------------------------------------------
  If ((L_CS).Or.(L_CSrp)) Then
   Write(io_L,100) "Block SPhenoCrossSections  # cross sections"
   Do i3=1,Size(Ecms)
    If (Ecms(i3).Eq.0._dp) Exit
   If (ISR(i3)) Then
    Write(io_L,4712) Ecms(i3),Pm(i3),Pp(i3)," 1  # e+ e- XS, Pe-, Pe+,  including ISR"
   Else
    Write(io_L,4712) Ecms(i3),Pm(i3),Pp(i3)," 0  # e+ e- XS"
   End If  
   Write(io_L,100) "#      Sigma [fb]    NDA        ID1     ID2"
   ! u-squarks
   Do i1=1,6
    Do i2=1,6
     If (SigSup(i3,i2,i1).Gt.SigMin) Write(io_L,4711) SigSup(i3,i2,i1), 2 &
        &      , id_su(i1), -id_su(i2), c_su(i1)//" "//Trim(c_su(i2))//"*"
    End Do
   End Do
   ! d-squarks
   Do i1=1,6
    Do i2=1,6
     If (SigSdown(i3,i2,i1).Gt.SigMin) Write(io_L,4711) SigSdown(i3,i2,i1), 2 &
        &      , id_sd(i1), -id_sd(i2), c_sd(i1)//" "//Trim(c_sd(i2))//"*"
    End Do
   End Do
   If (HighScaleModel.Eq."RPexplicit") Then
   Else
    ! sleptons
    Do i1=1,6
     Do i2=1,6
      If (SigSle(i3,i2,i1).Gt.SigMin) Write(io_L,4711) SigSle(i3,i2,i1), 2 &
        &      , id_sle(i1), -id_sle(i2), c_sle(i1)//" "//c_slep(i2)
     End Do
    End Do
    ! sneutrinos
    Do i1=1,3
     Do i2=1,3
      If (SigSn(i3,i2,i1).Gt.SigMin) Write(io_L,4711) SigSn(i3,i2,i1), 2 &
        &      , id_snu(i1), -id_snu(i2), c_snu(i1)//" "//Trim(c_snu(i2))//"*"
     End Do
    End Do
   End If
   ! neutralinos
   Do i1=1,n_n
    Do i2=i1,n_n
     If (SigChi0(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigChi0(i3,i1,i2), 2 &
        &      , id_n(i1), id_n(i2), c_c0(i1)//" "//c_c0(i2)
    End Do
   End Do
   ! charginos
   Do i1=1,n_c
    Do i2=1,n_c
     If (SigC(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigC(i3,i1,i2), 2 &
        &      , id_c(i1), -id_c(i2), c_cm(i1)//" "//c_cp(i2)
    End Do
   End Do
   
   If (HighScaleModel.Eq."RPexplicit") Then
    Do i1=1,n_s0
     If (SigS0(i3,i1).Gt.SigMin) Write(io_L,4711) SigS0(i3,i1),2, id_S0(i1) &
                                                & , id_Z , Trim(c_s0(i1))//" Z"
    End Do
    Do i1=1,n_s0
     Do i2=1,n_p0
      If (SigSP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSP(i3,i1,i2),2 &
                & , id_S0(i1), id_p0(i2) , Trim(c_s0(i1))//" "//Trim(c_p0(i2))
     End Do
    End Do
    Do i1=1,n_spm
     Do i2=1,n_spm
      If (SigHP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigHP(i3,i1,i2),2 &
                & , id_Sp(i1), id_sm(i2) , Trim(c_sp(i1))//" "//Trim(c_sm(i2))
     End Do
    End Do

   Else If (HighScaleModel.Eq."NMSSM") Then
    Do i1=1,n_s0
     If (SigS0(i3,i1).Gt.SigMin) Write(io_L,4711) SigS0(i3,i1),2, id_S0(i1) &
                                                & , id_Z , Trim(c_s0(i1))//" Z"
    End Do
    Do i1=1,n_s0
     Do i2=1,n_p0
      If (SigSP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSP(i3,i1,i2),2 &
                & , id_S0(i1), id_p0(i2) , Trim(c_s0(i1))//" "//Trim(c_p0(i2))
     End Do
    End Do
    If (SigHP(i3,1,1).Gt.SigMin) &
       & Write(io_L,4711) SigHP(i3,1,1), 2, id_Hp, -id_Hp , "H+ H-"

   Else
    If (SigS0(i3,1).Gt.SigMin) Write(io_L,4711) SigS0(i3,1),2, id_S0(1), id_Z &
                                                & , "h0 Z"
    If (SigS0(i3,2).Gt.SigMin) &
      & Write(io_L,4711) SigS0(i3,2), 2, id_S0(2), id_Z, "H0 Z"
    If (SigSP(i3,1,1).Gt.SigMin) &
      & Write(io_L,4711) SigSP(i3,1,1),2, id_S0(1), id_A0, "h0 A0"
    If (SigSP(i3,2,1).Gt.SigMin) &
       &  Write(io_L,4711) SigSP(i3,2,1),2, id_S0(2), id_A0 , "H0 A0"
    If (SigHP(i3,1,1).Gt.SigMin) &
       & Write(io_L,4711) SigHP(i3,1,1), 2, id_Hp, -id_Hp , "H+ H-"
    End If
   End Do
  End If ! L_CS

  !--------------------------------------------------------
  !information for HiggsBounds, for the moment only MSSM
  !--------------------------------------------------------
  If (HighScaleModel.Eq."NMSSM") Then  

  Else If (HighScaleModel.Eq."RPexplicit") Then
 
  Else
   If (L_BR) Call Write_HiggsBoundsInput(io_L,RS0, RP0)
  End If

  Write(io_L,100) "Block SPhenoLowEnergy  # low energy observables"
  Write(io_L,101) 1,1.e-4*BrBToSGamma," # BR(b -> s gamma)"
  Write(io_L,101) 2,BrBToSLL," # BR(b -> s mu+ mu-)"
  Write(io_L,101) 3,BToSNuNu," # BR(b -> s nu nu)"
  Write(io_L,101) 4,Bd_ll(1)," # BR(Bd -> e+ e-)"
  Write(io_L,101) 5,Bd_ll(2)," # BR(Bd -> mu+ mu-)"
  Write(io_L,101) 6,Bd_ll(3)," # BR(Bd -> tau+ tau-)"
  Write(io_L,101) 7,Bs_ll(1)," # BR(Bs -> e+ e-)"
  Write(io_L,101) 8,Bs_ll(2)," # BR(Bs -> mu+ mu-)"
  Write(io_L,101) 9,Bs_ll(3)," # BR(Bs -> tau+ tau-)"
  Write(io_L,101) 10,BR_Bu_TauNu," # BR(B_u -> tau nu)"
  Write(io_L,101) 11,R_Bu_TauNu," # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM"
  Write(io_L,101) 12,Abs(DeltaMbd)," # |Delta(M_Bd)| [ps^-1] "
  Write(io_L,101) 13,Abs(DeltaMbs)," # |Delta(M_Bs)| [ps^-1] "

  Write(io_L,101) 16, epsK       , " # epsilon_K" 
  Write(io_L,101) 17, DeltaMK    , " # Delta(M_K)" 
  Write(io_L,101) 18, K0toPi0NuNu, " # BR(K^0 -> pi^0 nu nu)" 
  Write(io_L,101) 19, KptoPipNuNu, " # BR(K^+ -> pi^+ nu nu)" 

  Write(io_L,101) 20,A_e," # Delta(g-2)_electron/2"
  Write(io_L,101) 21,A_mu," # Delta(g-2)_muon/2"
  Write(io_L,101) 22,A_tau," # Delta(g-2)_tau/2"
  Write(io_L,101) 23,d_e," # electric dipole moment of the electron"
  Write(io_L,101) 24,d_mu," # electric dipole moment of the muon"
  Write(io_L,101) 25,d_tau," # electric dipole moment of the tau"
  Write(io_L,101) 26,BrMutoEGamma," # Br(mu -> e gamma)"
  Write(io_L,101) 27,BrTautoEGamma," # Br(tau -> e gamma)"
  Write(io_L,101) 28,BrTautoMuGamma," # Br(tau -> mu gamma)"
  Write(io_L,101) 29,BrMu3E," # Br(mu -> 3 e)"
  Write(io_L,101) 30,BrTau3E," # Br(tau -> 3 e)"
  Write(io_L,101) 31,BrTau3Mu," # Br(tau -> 3 mu)"

  Write(io_L,101) 39,rho_parameter," # Delta(rho_parameter)"
  Write(io_L,101) 40,BR_Z_e_mu," # BR(Z -> e mu)"
  Write(io_L,101) 41,BR_Z_e_tau," # BR(Z -> e tau)"
  Write(io_L,101) 42,BR_Z_mu_tau," # BR(Z -> mu tau)"
  If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
   Write(io_L,100) "Block NuMass # neutrino masses"
   Write(io_L,101) 1,mf_nu(1)*1.e9_dp," # m_nu_1 in eV"
   Write(io_L,101) 2,mf_nu(2)*1.e9_dp," # m_nu_2 in eV"
   Write(io_L,101) 3,mf_nu(3)*1.e9_dp," # m_nu_3 in eV"
   Write(io_L,101) 4,(mf_nu(3)**2-mf_nu(2)**2)*1.e18_dp &
                & ," # Delta(m^2_atm) in eV^2"
   Write(io_L,101) 5,(mf_nu(2)**2-mf_nu(1)**2)*1.e18_dp &
                & ," # Delta(m^2_sol) in eV^2"
!   Call WriteMatrixBlockC2(io_L, 3, Unu, "Unu" &
!                         &, "neutrino mixing matrix", "U_(nu,")
  End If
  
  Call Write_Wilson_Coeff(io_L,1,kont)  ! Wilson coefficients at m_Z
  Call Write_Wilson_Coeff(io_L,2,kont)  ! Wilson coefficients at Q=160 ~ m_t(m_t)

  If (Present(omega)) Then
   Write(io_L,100) "Block Omega # omega h^2"
   Write(io_L,101) 1,omega," # omega h^2"   
  End If

  If (LWrite_LHC_Observables) Then
   If (HighScaleModel.Eq."RPexplicit") Then
    Do i1=1,7
     If (Abs(id_sp(i1)).Eq.1000011) Then
      mSle(2) = Spm8(i1)%m
     Else If (Abs(id_sp(i1)).Eq.2000011) Then
      mSle(1) = Spm8(i1)%m
     Else If (Abs(id_sp(i1)).Eq.1000013) Then
      mSle(4) = Spm8(i1)%m
     Else If (Abs(id_sp(i1)).Eq.2000013) Then
      mSle(3) = Spm8(i1)%m
     Else If (Abs(id_sp(i1)).Eq.1000015) Then
      mSle(5) = Spm8(i1)%m
     Else If (Abs(id_sp(i1)).Eq.2000015) Then
      mSle(6) = Spm8(i1)%m
     End If
    End Do
   Else ! conserved R-parity
    Do i1=1,6
     If (id_sle(i1).Eq.1000011) Then
      mSle(2) = Slepton(i1)%m
     Else If (id_sle(i1).Eq.2000011) Then
      mSle(1) = Slepton(i1)%m
     Else If (id_sle(i1).Eq.1000013) Then
      mSle(4) = Slepton(i1)%m
     Else If (id_sle(i1).Eq.2000013) Then
      mSle(3) = Slepton(i1)%m
     Else If (id_sle(i1).Eq.1000015) Then
      mSle(5) = Slepton(i1)%m
     Else If (id_sle(i1).Eq.2000015) Then
      mSle(6) = Slepton(i1)%m
     End If
    End Do
   End If ! R-parity

   Do i1=1,6
    If (id_su(i1).Eq.1000002) Then
     mSu(2) = Sup(i1)%m
    Else If (id_su(i1).Eq.2000002) Then
     mSu(1) = Sup(i1)%m
    Else If (id_su(i1).Eq.1000004) Then
     mSu(4) = Sup(i1)%m
    Else If (id_su(i1).Eq.2000004) Then
     mSu(3) = Sup(i1)%m
    Else If (id_su(i1).Eq.1000006) Then
     mSu(5) = Sup(i1)%m
    Else If (id_su(i1).Eq.2000006) Then
     mSu(6) = Sup(i1)%m
    End If

    If (id_sd(i1).Eq.1000001) Then
     mSd(2) = Sdown(i1)%m
    Else If (id_sd(i1).Eq.2000001) Then
     mSd(1) = Sdown(i1)%m
    Else If (id_sd(i1).Eq.1000003) Then
     mSd(4) = Sdown(i1)%m
    Else If (id_sd(i1).Eq.2000003) Then
     mSd(3) = Sdown(i1)%m
    Else If (id_sd(i1).Eq.1000005) Then
     mSd(5) = Sdown(i1)%m
    Else If (id_sd(i1).Eq.2000005) Then
     mSd(6) = Sdown(i1)%m
    End If
   End Do

   If (HighScaleModel.Eq."RPexplicit") Then
    Call Calc_LHC_observables(Chi07(4:7)%m, N(1:4,1:4), ChiPm%m, U, V, mSle, Rslepton &
      & , mSd, RSdown, mSu, RSup, Glu%m, PhaseGlu, S0%m, RS0, P0%m, RP0         &
      & , Spm%m, RSpm, gauge, Y_u, Y_d, A_u, A_d, mu, vevSM, .False. & ! GenerationMixing &
      & , LHC_observ)
    n_min = 4
   Else ! conserved R-parity
    Call Calc_LHC_observables(Chi0%m, N, ChiPm%m, U, V, mSle, Rslepton        &
      & , mSd, RSdown, mSu, RSup, Glu%m, PhaseGlu, S0%m, RS0, P0%m, RP0     &
      & , Spm%m, RSpm, gauge, Y_u, Y_d, A_u, A_d, mu, vevSM, .False. & ! GenerationMixing &
      & , LHC_observ)
    n_min = 1
   End If ! R-parity

   Write(io_L,100) "Block LHCobservables # edge observables for LHC"
   i_zaehl = 1
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # e+e- edge with right selectron, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # e+e- edge with left selectron, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # mu+mu- edge with right smuon, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # mu+mu- edge with left smuon, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # tau+tau- edge with lighter stau, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=n_min+1,n_n
    Do i2=n_min,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # tau+tau- edge with heavier stau, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- q edge, averaging over d_L, s_L, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- q threshold, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+-_near q edge, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+-_far q edge, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- b threshold"
   i_zaehl = i_zaehl + 1


  End If
!  Close(io_L)

 98 Format(2x,i2,3x,i3,2x,"# ",a)
 99 Format(1x,i5,3x,a)
100 Format(a)
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a)
103 Format(a13,1P,e16.8,2x,a)
104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,3x,a)
108 Format(9x,1P,E16.8,0P,3x,a)
109 Format(1x,3i3,3x,1P,e16.8,3x,a)
110 Format(3x,2i3,3x,"# ",a)
111 Format(2i5,4x,"# ",a)
200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a)
300 Format("# start writing paramaters at Q=",1P,E16.8)
301 Format("# start writing paramaters and masses at Q=",1P,E16.8)
310 Format("# stop writing paramaters at Q=",1P,E16.8)
311 Format("# stop writing paramaters and masses at Q=",1P,E16.8)

401 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x,"# BR(",a," -> ",a," ",a,")")
402 Format(3x,1P,e16.8,0p,3x,I2,3x,3(i9,1x),2x,"# BR(",a," -> ",a," ",a," ",a,")")

4711 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x," # ",A)
4712 Format("XS 11 -11 ",F7.1," ",F5.2," ",F5.2," ",A)

 End Subroutine LesHouches_Out


 Subroutine LesHouches_out_start(io_L, io, kont, use_rge, Q)
 Implicit None
  Integer, Intent(in) :: io_L, io, kont, use_rge
  Real(dp), Intent(out) :: Q

  Character(len=8) :: Datum
  Character(len=10) :: Zeit
  Logical, Save :: l_open = .True. ! in case of a loop I want to open the
                                   ! output only once 
  Integer :: i_errors(1100)

  Q = Sqrt( GetRenormalizationScale() )

  Call Date_and_time(datum,zeit)
  If (l_open) Then
   Open(io_L,file="SPheno.spc",status="unknown")
   l_open = .False.
  End If
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2.beta - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno "//version
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) &
     & "# in case of problems send email to porod@physik.uni-wuerzburg.de"
   Write(io_L,100) "# Created: "//Datum(7:8)//"."//Datum(5:6)//"."//Datum(1:4) &
     & //",  "//Zeit(1:2)//":"//Zeit(3:4)
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   "//version//"    # version number"
   !-----------------------------------------------
   ! check if somewhere a problem has had happened
   !-----------------------------------------------
   Call GetError(i_errors)
   !--------------------------------------
   ! a numerical problem might have happen
   !--------------------------------------
   If ((i_errors(1)+i_errors(3)+i_errors(5)+i_errors(7)+i_errors(8) &
     & + i_errors(10) + i_errors(12)+ Sum(i_errors(14:19))).Gt.0)   &
     & Write(io_L,100) &
 & "     3               # potential numerical problem, check file Messages.out"
   If (in_kont(1).Eq.1) Write(io_L,99) 3, &
    & "alpha(0) and alpha(mZ) have both been specified without check for"// &
    &  " consistency"
   If (in_kont(2).Eq.1) Write(io_L,99) 3, &
    & "redundant specification in Higgs sector"
   If (use_rge.Eq.0) Then
    Write(io_L,100) "#"
    Write(io_L,100) "Block SPhenoINFO     # SPheno specific information"
    If (TwoLoopRGE) Then
     Write(io_L,100) "    1      2         # using 2-loop RGEs"
    Else 
     Write(io_L,100) "    1      1         # using 1-loop RGEs"
    End If
    If (YukScen.Eq.1) Then
     Write(io_L,100) &
    &"    2      1         # using running masses for boundary conditions at mZ"
    Else
     Write(io_L,100) &
       &"    2      2         # using pole masses for boundary conditions at mZ"
    End If
   End If
   !------------------------------
   ! SPheno.out
   !------------------------------
  Write(io,123) "SPheno output file"
  Write(io,123) "Version "//version//" ,  "//"created: "//Datum(7:8)//"."// &
     & Datum(5:6)//"."//Datum(1:4)//",  "//Zeit(1:2)//":"//Zeit(3:4)
  Write(io,*) " "
123 Format(a)

  If (kont.Ne.0) Then
   Write(io,*) "There has been a problem during the run."
   Write(io,*) "Please check the file Messages.out for further information."
   Write(*,*) "There has been a problem during the run."
   Write(*,*) "Please check the file Messages.out for further information."
   Write(io,*) " "
  End If

  Write(io,*) " "
  If (YukScen.Eq.1) Then
   Write(io,*) &
     & "Running masses have been used for the boundary conditions at mZ" 
  Else If (YukScen.Eq.2) Then
   Write(io,*) "Pole masses have been used for the boundary conditions at mZ" 
  End If
  Write(io,*) "  Using Yukawa scheme :",GetYukawaScheme()
  Write(io,*) " "

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"
   !------------------------
   ! input CKM
   !------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
   End If

  99 Format(1x,i5,3x,a)
  100 Format(a)
  102 Format(1x,i5,3x,1P,e16.8,2x,a)

 End Subroutine LesHouches_out_start

 Subroutine LH_Write_GY(io_L, Q, g, Yd, Yu, Yl, Ynu, YT)
 implicit none
  Integer, intent(in) :: io_L
  Real(dp), Intent(in) :: Q, g(3), Yd(3), Yu(3)
  Complex(dp), Optional, Intent(in) :: Yl(3,3), Ynu(3,3), YT(3,3)

  Integer :: i1, i2, ierr

  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,g(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,g(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,g(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  Write(io_L,106) "Block Ye Q=",Q,"# (SUSY scale)"

  ierr = 0
  If (GenerationMixing) Then
   !-------------------------------------------
   ! check if any off-diagonal term is non-zero
   !-------------------------------------------
   Do i1=1,3
    Do i2=1,3
     If ((i1.Ne.i2).And.(Abs(yl(i2,i1)).Ne.0._dp)) ierr = ierr + 1
    End Do
   End Do
  End If

  If (ierr.Ne.0) Then
   Do i1=1,3
    Do i2=1,3
     Write(io_L,105) i2,i1,Real(yl(i2,i1),dp),"# Y_(l,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
  Else 
   Write(io_L,107) 1,1,Real(Yl(1,1),dp), "# Y_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Yl(2,2),dp), "# Y_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Yl(3,3),dp), "# Y_tau(Q)^DRbar"
  End If
  If (Maxval(Abs(Aimag(yl))).Gt.0._dp) Then
   Write(io_L,106) "Block IMYe Q=",Q,"# (SUSY scale)"
   If (GenerationMixing) Then
    Do i1=1,3
      Do i2=1,3
      Write(io_L,105) i2,i1,Aimag(yl(i2,i1)),"# Im(Y_(l,"//bu(i2)//bu(i1)//"))"
     End Do
    End Do
   Else 
    Write(io_L,107) 1,1,Aimag(Yl(1,1)), "# Im(Y_e)(Q)^DRbar"
    Write(io_L,107) 2,2,Aimag(Yl(2,2)), "# Im(Y_mu)(Q)^DRbar"
    Write(io_L,107) 3,3,Aimag(Yl(3,3)), "# Im(Y_tau)(Q)^DRbar"
   End If
  End If

104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,4x,a)

 end Subroutine LH_Write_GY

 Subroutine LH_write_decays(io_L, io, n_c, n_n, n_sl, n_sn, n_sd, n_su, n_S0 &
          & , n_P0, n_Spm, id_l, id_nu, id_d, id_u, id_W, id_Z, id_ph, id_gl &
          & , id_glu, id_grav, id_c, id_n, id_sle, id_snu, c_min, n_min      &
          & , id_sd, id_su, id_S0, id_P0, id_Sm, id_Sp, c_gl, c_lm           &
          & , c_lp, c_nu, c_d, c_u, c_Cm, c_Cp, c_c0, c_sle, c_slep, c_Snu   &
          & , c_Sd, c_Su, c_S0, c_P0, c_Sm, c_Sp, c_grav, c_phot             &
          & , gT_C, gP_C2, gP_C3, BR_C2, BR_C3                               &
          & , gT_N, gP_N2, gP_N3, BR_N2, BR_N3, gT_S0, gP_S0, BR_S0, gT_P0   &
          & , gP_P0, BR_P0, gT_Spm, gP_Spm, BR_Spm )

 implicit none
  Integer, Intent(in) :: io_L, io, n_c, n_n, n_sl, n_sn, n_sd, n_su, n_S0     &
   & , n_P0, n_Spm, id_l(3), id_nu(3), id_d(3), id_u(3), id_c(:), id_n(:)     &
   & , id_sle(:), id_snu(:), id_sd(:), id_su(:), id_S0(:), id_P0(:), id_Sm(:) &
   & , id_Sp(:), id_W, id_Z, id_glu, id_grav, id_gl, id_ph, c_min, n_min 
  Real(dp), Intent(in) :: gT_N(:), gP_N2(:,:), gP_N3(:,:), BR_N2(:,:)       &
   & , BR_N3(:,:), gT_S0(:), gP_S0(:,:), BR_S0(:,:), gT_P0(:), gP_P0(:,:)   &
   & , BR_P0(:,:), gT_Spm(:), gP_Spm(:,:), BR_Spm(:,:), gT_C(:), gP_C2(:,:) &
   & , gP_C3(:,:), BR_C2(:,:), BR_C3(:,:)
  Character(len=10), intent(in) :: c_snu(:), c_sle(:), c_slep(:)
  Character(len=6), intent(in) :: c_lm(3), c_lp(3), c_nu(3), c_sp(:), c_sm(:) &
   & , c_phot
  Character(len=1), intent(in) :: c_u(3), c_d(3)
  Character(len=2), intent(in) :: c_gl, c_grav
  Character(len=4), intent(in) :: c_S0(:), c_P0(:), c_su(6), c_sd(6)
  Character(len=5), intent(in) :: c_cp(:), c_cm(:), c_c0(:)

  Character(len=13) :: c_sfermion 

  Integer :: i_zaehl, id_d_2(200,2), id_d_3(400,3), i1, i2, i3, i4, id_f &
     & , id_fp, i_c2, i_eff
  Real(dp) :: BRmin100, BRtot, gT
  Logical :: l_3body
  Integer, Parameter :: n_max=500
  Real(dp), Dimension(n_max) :: Br, gP
  Character(len=30), Dimension(n_max) :: Fnames, Lnames
  Character(len=10) :: c_m

  Write(io,*) " "
  Write(io,*) " Anti particles are marked with a * in case of"
  Write(io,*) " (s)neutrinos and (s)quarks in the decay section."
  Write(io,*) "                    Decay widths (GeV) and branching ratios"
  Write(io,*) " "

  BRmin100 = 100._dp * BRmin
  !---------------------------------
  ! sleptons
  !---------------------------------
  Do i1=1,n_sl
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Do i3 = 1,3
     Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
     Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_grav
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do

    Call WriteDecays2(io, " slepton_"//Bu(i1) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    i2=id_fp
    If (i2.Eq.1) Then
     c_sfermion = "e-sneutrino"
    Else If (i2.Eq.2) Then
     c_sfermion = "mu-sneutrino"
    Else If (i2.Eq.3) Then
     c_sfermion = "tau-sneutrino"
    End If
    Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
    Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
    id_d_2(i_zaehl,1) = id_snu(i2)
    id_d_2(i_zaehl,2) = -id_W
    i_zaehl = i_zaehl+1
    Do i3=1,n_spm
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     i_zaehl = i_zaehl+1
    End Do
    If (i1.Le.2) Then
     c_sfermion = "selectron"
    Else If (i1.Le.4) Then
     c_sfermion = "smuon"
    Else 
     c_sfermion = "stau"
    End If

    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    i3 = (i1+1)/2
    Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
    Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
    id_d_2(i_zaehl,1) = id_grav
    id_d_2(i_zaehl,2) = id_l(i3)
    i_zaehl = i_zaehl+1

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_sle(i1),gT_sl(i1),Trim( c_sle(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sl(i1,i2).Gt.BrMin) Write(io_L,201) BR_Sl(i1,i2),2,id_d_2(i2,:), &
            &                  Trim( c_sle(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! sneutrinos
  !---------------------------------
  Do i1=1,n_sn
   If (GenerationMixing) Then
    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_nu(i3)
      i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sl
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sl
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" Z"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " Sneutrino_"//Bu(i1) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i1
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_l(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    If (i1.Eq.1) Then
     c_sfermion = "selectron_"
    Else If (i1.Eq.2) Then
     c_sfermion = "smuon_"
    Else If (i1.Eq.3) Then
     c_sfermion = "stau_"
    End If
    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Eq.1) c_sfermion="e-sneutrino"
    If (i1.Eq.2) c_sfermion="mu-sneutrino"
    If (i1.Eq.3) c_sfermion="tau-sneutrino"

    Call WriteDecays2(io, Trim(c_sfermion) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_snu(i1),gT_sn(i1),Trim(c_snu(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sn(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sn(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_snu(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! d-squarks
  !---------------------------------
  Do i1=1,n_sd
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_c(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-down_"//Bu(i1) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_u(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_d(i3)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
     id_fp = 2
    Else 
     c_sfermion = "stop_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2 + id_fp))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2 + id_fp)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2+id_fp))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
    Else 
     c_sfermion = "sbottom_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_sd(i1),gT_sd(i1),Trim(c_sd(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sd(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sd(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_sd(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! u-squarks
  !---------------------------------
  Do i1=1,n_su
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-up_"//Bu(i1) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_d(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_u(id_fp))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(id_fp))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_u(id_fp)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
     id_fp = 2
    Else 
     c_sfermion = "sbottom_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2 + id_fp))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2 + id_fp)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2+id_fp))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
    Else 
     c_sfermion = "stop_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2

    If (i1.Eq.5) Then ! 3-body decays of lighter stop
     id_d_3 = 0
     Fnames(55) = "neutralino_1 c"
     Lnames(55) = Trim(c_c0(1))//" c"
     id_d_2(55,1) = id_n(1)
     id_d_2(55,2) = id_u(2)
     Fnames(56) = "neutralino_2 c"
     Lnames(56) = Trim(c_c0(2))//" c"
     id_d_2(56,1) = id_n(2)
     id_d_2(56,2) = id_u(2)
     i_zaehl = 57
     Fnames(57) = "neutralino_1 W b"
     Lnames(57) = Trim(c_c0(1))//" W+ b"
     id_d_3(1,1) = id_n(1)
     id_d_3(1,2) = id_W
     id_d_3(1,3) = id_d(3)
     Fnames(58) = "e-sneutrino e+ b"
     Lnames(58) = Trim(c_snu(1))//" e+ b"
     id_d_3(2,1) = id_snu(1)
     id_d_3(2,2) = -id_l(1)
     id_d_3(2,3) = id_d(3)
     Fnames(59) = "mu-sneutrino mu+ b"
     Lnames(59) = Trim(c_snu(2))//" mu+ b"
     id_d_3(3,1) = id_snu(2)
     id_d_3(3,2) = -id_l(2)
     id_d_3(3,3) = id_d(3)
     Fnames(60) = "tau-sneutrino tau+ b"
     Lnames(60) = Trim(c_snu(3))//" tau+ b"
     id_d_3(4,1) = id_snu(3)
     id_d_3(4,2) = -id_l(3)
     id_d_3(4,3) = id_d(3)
     Fnames(61) = "selectron_1 nu_e b"
     Lnames(61) = Trim(c_slep(1))//" nu_e b"
     id_d_3(5,1) = -id_sle(1)
     id_d_3(5,2) = id_nu(1)
     id_d_3(5,3) = id_d(3)
     Fnames(62) = "selectron_2 nu_e b"
     Lnames(62) = Trim(c_slep(2))//" nu_e b"
     id_d_3(6,1) = -id_sle(2)
     id_d_3(6,2) = id_nu(1)
     id_d_3(6,3) = id_d(3)
     Fnames(63) = "smuon_1 nu_mu b"
     Lnames(63) = Trim(c_slep(3))//" nu_mu b"
     id_d_3(7,1) = -id_sle(3)
     id_d_3(7,2) = id_nu(2)
     id_d_3(7,3) = id_d(3)
     Fnames(64) = "smuon_2 nu_mu b"
     Lnames(64) = Trim(c_slep(4))//" nu_mu b"
     id_d_3(8,1) = -id_sle(4)
     id_d_3(8,2) = id_nu(2)
     id_d_3(8,3) = id_d(3)
     Fnames(65) = "stau_1 nu_tau b"
     Lnames(65) = Trim(c_slep(5))//" nu_tau b"
     id_d_3(9,1) = -id_sle(5)
     id_d_3(9,2) = id_nu(3)
     id_d_3(9,3) = id_d(3)
     Fnames(66) = "stau_2 nu_tau b"
     Lnames(66) = Trim(c_slep(6))//" nu_tau b"
     id_d_3(10,1) = -id_sle(6)
     id_d_3(10,2) = id_nu(3)
     id_d_3(10,3) = id_d(3)
    End If

    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_su(i1),gT_su(i1),Trim(c_su(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Su(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Su(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_su(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do

   If ((i1.Eq.5).And.(Maxval(BR_Su(i1,57:66)).Gt.BRmin)) Then
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,10
      If  (BR_Su(i1,56+i2).Gt.BrMin) Write(io_L,202) BR_Su(i1,56+i2),3 &
             & ,id_d_3(i2,:), Trim(c_su(i1))//" -> "//Trim(Lnames(56+i2))//")"
    End Do
   End If

  End Do   

  !--------------
  ! charginos
  !--------------

  Do i1=c_min,n_c
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     i3 = (i2+1)/2
     Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Do i3 = 1,3
      Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_lp(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-up_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-down_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_sd(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_lp(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = -id_sd(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do

   End If ! GenerationMixing

   Do i2=1,n_n
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" W+"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" W+"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_W
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_spm
    Do i2=1,n_n
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_sp(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i2=1,i1-1
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" Z"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp

   i_eff = i1 - (n_c-2)

   gP(1:i_c2) = gP_C2(i_eff,1:i_c2)
   BR(1:i_c2) = BR_C2(i_eff,1:i_c2)

   Write(io_L,200) id_c(i1),gT_C(i_eff),Trim(c_cp(i1))
   If (Sum(gp_C2(i_eff,:)).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR_C2(i_eff,i2).Gt.BrMin) Write(io_L,201) BR_C2(i_eff,i2) &
            & ,2,id_d_2(i2,:), Trim(c_cp(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   If (Maxval(BR_C3(i_eff,:)).Gt.BRmin) Then
    If (GenerationMixing) Then
     Do i2=1,n_n
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_n(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3-(n_c-2)
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_glu
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,5-n_c
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
        id_d_3(i_zaehl-i_c2,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      If ((n_c-2).Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_j"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_j"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(3)
       i_zaehl = i_zaehl+1    
      End If
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_n
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      id_d_3(i_zaehl-i_c2,1) = id_glu
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i3)
      i_zaehl = i_zaehl+1    
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      If ((n_c-2).Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_i"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_i"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(1)
       i_zaehl = i_zaehl+1    
      End If
     End Do

    End If ! GenerationMixing

    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,i_zaehl - i_c2 - 1
     If (BR_C3(i_eff,i2).Gt.BrMin) Write(io_L,202) BR_C3(i_eff,i2),3 &
           & ,id_d_3(i2,:) , Trim(c_cp(i1))//" -> "//Trim(Lnames(i2+i_c2))//")"
    End Do

    BR(i_c2+1:i_zaehl-1) = BR_C3(i_eff,1:i_zaehl-i_c2-1)
    gP(i_c2+1:i_zaehl-1) = gP_C3(i_eff,1:i_zaehl-i_c2-1)
   End If ! check for maxval(BR)

   Call WriteDecays2(io, " chargino_"//Bu(i1) , Fnames &
                       &, gP, 100*BR, gT_C(i_eff), BrMin100)

  End Do ! i1


  !--------------
  ! neutralinos
  !--------------
  Do i1=n_min,n_n
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     Do i3 = 1,3
      Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_lp(i3))
      Fnames(i_zaehl+1) = "slepton_"//Bu(i2)//"^+ "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
      Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      id_d_2(i_zaehl+1,1) = -id_sle(i2)
      id_d_2(i_zaehl+1,2) = id_l(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sn
     i3 = i2
     Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = "sneutrino_"//Bu(i2)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      id_d_2(i_zaehl+1,1) = -id_su(i2)
      id_d_2(i_zaehl+1,2) = id_u(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
      Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      id_d_2(i_zaehl+1,1) = -id_sd(i2)
      id_d_2(i_zaehl+1,2) = id_d(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_lp(i3))
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^+ "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
     Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     id_d_2(i_zaehl+1,1) = -id_sle(i2)
     id_d_2(i_zaehl+1,2) = id_l(i3)
     i_zaehl = i_zaehl+2
    End Do

    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
    
   End If ! GenerationMixing

   Do i2=1,n_c
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ W-"
    Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- W+"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" W-"
    Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" W+"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = -id_W
    id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl+2 
   End Do

   Do i3=1,n_spm
    Do i2=1,n_c
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = -id_sp(i3)
     id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl+2   
    End Do
   End Do


   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" Z"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Fnames(i_zaehl) = "Gravitino photon"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_ph
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino Z"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino h0"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_s0(1)
   i_zaehl = i_zaehl + 1

   If (HighScaleModel.Eq."NMSSM") Then
    i_zaehl = 150 
   Else If (HighScaleModel.Eq."RPexplicit") Then
    i_zaehl = 186
   Else If (HighScaleModel.Eq."NURRP1") Then
    i_zaehl = 190
   Else
    i_zaehl = 150 
   End If

   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" photon"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" photon"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_ph
    i_zaehl = i_zaehl+1    
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp
   
   i_eff = i1 + 1 - n_min

   gT = gT_N(i_eff)
   gP(1:i_c2) = gP_N2(i_eff,1:i_c2)
   BR(1:i_c2) = BR_N2(i_eff,1:i_c2)

   Write(io_L,200) id_n(i1),gT,Trim(c_c0(i1))
   If (Sum(gP).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,201) BR(i2),2,id_d_2(i2,:), &
            &                   Trim(c_c0(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   l_3body = (Maxval(BR_N3(i_eff,:)).Gt.BRmin)
   i_zaehl = i_zaehl + (n_n-i1+1)
   i_c2 = i_zaehl - 1

   If (l_3body) Then
    If (GenerationMixing) Then
     Do i2=1,n_c
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_c(i2)
        id_d_3(i_zaehl,2) = id_d(i3)
        id_d_3(i_zaehl,3) = -id_u(i4)
        id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
        i_zaehl = i_zaehl+2
       End Do
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1+(n_n-4)
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_j"
      Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_j"
      id_d_3(i_zaehl,1) = id_n(i2)
      id_d_3(i_zaehl,2) = - id_nu(1)
      id_d_3(i_zaehl,3) = id_nu(3)
      i_zaehl = i_zaehl+1    
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
     
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_c
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_d(i3)
       id_d_3(i_zaehl,3) = -id_u(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

      Do i3=1,3
       i4 = i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i3=1,3
       i4=i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i2=1,i1-1
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
       If ((n_n-4).Lt.3) Then
        Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_i"
        Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_i"
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_nu(1)
        id_d_3(i_zaehl,3) = id_nu(3)
        i_zaehl = i_zaehl+1    
       End If
       Do i3=1,7-n_n
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do    
      End Do
 
     If (     (HighScaleModel.Eq."RPexplicit") &
        & .or.(HighScaleModel(1:5).Eq."NURRP") ) Then
      Do i2=1,n_c
       Do i3=1,Min(i2,3)
        If (i2.Eq.i3) Then
         Fnames(i_zaehl) = &
            & "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_l(i2)
         id_d_3(i_zaehl,3) = id_l(i3)
         i_zaehl = i_zaehl+1
        Else
         Fnames(i_zaehl) = "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Fnames(i_zaehl+1) = &
            & "neutrino_e"//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl+1) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_c(i2)
         id_d_3(i_zaehl,3) = id_c(i3)
         id_d_3(i_zaehl+1,1) = id_n(1)
         id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
         i_zaehl = i_zaehl+2
        End If
       End Do
      End Do
      Do i2=4,i1-1
       Do i3=1,3
        Do i4=1,i3
         If (i3.Eq.i4) Then
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          i_zaehl = i_zaehl+1
         Else
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Fnames(i_zaehl+1) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl+1) = &
            & Trim(c_c0(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          id_d_3(i_zaehl+1,1) = id_n(i2)
          id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
          i_zaehl = i_zaehl+2
         End If
        End Do
       End Do
      End Do
      Do i2=4,i1-1
       Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
       Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
       id_d_3(i_zaehl,1) = id_n(i2)
       id_d_3(i_zaehl,2) = id_n(1)
       id_d_3(i_zaehl,3) = id_n(1)
       i_zaehl = i_zaehl+1
      End Do
      Fnames(i_zaehl) = Trim(c_nu(1))//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
      Lnames(i_zaehl) =  Trim(c_c0(1))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
      id_d_3(i_zaehl,1) = id_n(1)
      id_d_3(i_zaehl,2) = id_n(1)
      id_d_3(i_zaehl,3) = id_n(1)
      i_zaehl = i_zaehl+1
     End If
    End If ! GenerationMixing
   
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"

    BR(i_c2+1:i_zaehl-1) = BR_N3(i_eff,1:i_zaehl-i_c2-1)
    gP(i_c2+1:i_zaehl-1) = gP_N3(i_eff,1:i_zaehl-i_c2-1)

    Do i2=i_c2 + 1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,202) BR(i2),3,id_d_3(i2,:) &
             &                , Trim(c_c0(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do

   End If

  Call WriteDecays2(io, " neutralino_"//Bu(i1), Fnames, gP, 100*BR, gT, BrMin100)
  End Do

  !--------------
  ! gluino
  !--------------
  i_zaehl = 1
  If (GenerationMixing) Then
   Do i2=1,n_su
    Do i3 = 1,3
     Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
    
   Do i2=1,n_sd
    Do i3 = 1,3
     Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
 
    
  Else ! GenerationMixing

   Do i2=1,n_su
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-up_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-charm_"
     id_f = i2 - 2
    Else 
     c_sfermion="stop_"
     id_f = i2 - 2
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
    Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
    id_d_2(i_zaehl,1) = id_su(i2)
    id_d_2(i_zaehl,2) = -id_u(i3)
    id_d_2(i_zaehl+1,1) = -id_su(i2)
    id_d_2(i_zaehl+1,2) = id_u(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   Do i2=1,n_sd
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-down_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-strange_"
     id_f = i2 - 2
    Else 
     c_sfermion="sbottom_"
     id_f = i2 - 4
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_sd(i2)
    id_d_2(i_zaehl,2) = -id_d(i3)
    id_d_2(i_zaehl+1,1) = -id_sd(i2)
    id_d_2(i_zaehl+1,2) = id_d(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   id_f = 1
   c_sfermion="stop_"
   Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(2))//"^*"
   Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(2))
   Lnames(i_zaehl) =  Trim(c_su(5))//" "//Trim(c_u(2))//"^*"
   Lnames(i_zaehl+1) =  Trim(c_su(5))//" "//Trim(c_u(2))
   id_d_2(i_zaehl,1) = id_su(5)
   id_d_2(i_zaehl,2) = -id_u(2)
   id_d_2(i_zaehl+1,1) = -id_su(5)
   id_d_2(i_zaehl+1,2) = id_u(2)
   i_zaehl = i_zaehl+2

  End If ! GenerationMixing

  Do i2=1,n_n
   Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" gluon"
   Lnames(i_zaehl) =  Trim(c_c0(i2))//" gluon"
   id_d_2(i_zaehl,1) = id_n(i2)
   id_d_2(i_zaehl,2) = id_gl
   i_zaehl = i_zaehl+1    
  End Do

  i_c2 = i_zaehl-1
  gP = 0._dp
  BR = 0._dp
  gP(1:i_c2) = gP_Glu2(1:i_c2)
  BR(1:i_c2) = BR_Glu2(1:i_c2)

  Write(io_L,200) id_glu,gT_Glu,Trim(c_gl)
  If (Sum(gp_Glu2).Gt.0._dp) Then ! 2-body decays
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1 
    If (BR_Glu2(i2).Gt.BrMin) Write(io_L,201) BR_Glu2(i2),2,id_d_2(i2,:), &
            &                   Trim(c_gl)//" -> "//Trim(Lnames(i2))//")"
   End Do
  End If

  !-------------------------------------------------
  ! 3-body decays
  !-------------------------------------------------
  If (Maxval(BR_Glu3).Gt.BRmin) Then
   If (GenerationMixing) Then
    Do i2=1,n_n
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
    End Do

    Do i2=1,n_c
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
       id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do
    End Do 

    Do i2=1,6
     Do i3=1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//"^* W+ "//Trim(c_d(i3))
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//" W- "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
      Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
      id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
      id_d_3(i_zaehl-i_c2,2) = id_W
      id_d_3(i_zaehl-i_c2,3) = id_d(i3)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do  

   Else ! GenerationMixing
    Do i2=1,n_n
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i4)
      i_zaehl = i_zaehl+1    
     End Do
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_d(i4)
      i_zaehl = i_zaehl+1    
     End Do
    End Do

    Do i2=1,n_c
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_c(i2)
      id_d_3(i_zaehl-i_c2,2) = id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do 

    Do i2=5,6
     i3=3
     Fnames(i_zaehl) = "stop_"//Bu(i2-4)//"^* W+ "//Trim(c_d(i3))
     Fnames(i_zaehl+1) = "stop_"//Bu(i2-4)//" W- "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
     Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
     id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
     id_d_3(i_zaehl-i_c2,2) = id_W
     id_d_3(i_zaehl-i_c2,3) = id_d(i3)
     id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
     i_zaehl = i_zaehl+2
    End Do  
   End If ! GenerationMixing

   Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
   Do i2=1,i_zaehl - i_c2 - 2
     If (BR_Glu3(i2).Gt.BrMin) Write(io_L,202) BR_glu3(i2),3,id_d_3(i2,:) &
             &                , Trim(c_gl)//" -> "//Trim(Lnames(i2+i_c2))//")"
   End Do

   BR(i_c2+1:i_zaehl-2) = BR_Glu3(1:i_zaehl-i_c2-2)
   gP(i_c2+1:i_zaehl-2) = gP_glu3(1:i_zaehl-i_c2-2)

  End If

  Call WriteDecays2(io, " gluino_" , Fnames, gP, 100*BR, gT_Glu, BrMin100)


  !-----------------
  ! neutral scalars
  !-----------------
  Do i1=1,n_s0
   c_m = c_s0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 5 - n_c
     Do i3 = i2, 5 - n_c
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1, n_sn
     Do i3 = i2, n_sn
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'sneutrino^*_'//Bu(i2)//' sneutrino_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,5 - n_c
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! sneutrinos 
    Do i2=1,n_sn
     Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i2))
     id_d_2(i_zaehl,1) = -id_snu(i2)
     id_d_2(i_zaehl,2) = id_snu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If ((n_c-2).Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If ((n_c-2).Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If ((n_c-2).Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     Fnames(17) = 'e-Sneutrino'
     Fnames(18) = 'mu-Sneutrino'
     i2 = 18
    Else If ((n_c-2).Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     Fnames(22) = 'e-Sneutrino'
     Fnames(23) = 'mu-Sneutrino'
     Fnames(24) = 'tau-Sneutrino'
     i2 = 24
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n
    Do i3=i2,n_n
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c
    Do i3=i2, n_c
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Fnames(i_zaehl) = 'Z0 Z0'
   Lnames(i_zaehl) = 'Z0 Z0'
   id_d_2(i_zaehl,1) = id_Z
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1

   Fnames(i_zaehl) = 'W+ W-'
   Lnames(i_zaehl) = 'W- W+'
   id_d_2(i_zaehl,1) = id_W
   id_d_2(i_zaehl,2) = - id_W
   i_zaehl = i_zaehl + 1

   Do i2=1,n_p0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_P0(i2)
    i_zaehl = i_zaehl + 1
   End Do
   Do i2=1,n_p0
    Do i3=i2,n_p0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_P0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,i1-1
    Do i3=i2,i1-1
     Fnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_S0(i2)
     id_d_2(i_zaehl,2) = id_S0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Fnames(i_zaehl) = "g g"
   Fnames(i_zaehl + 1) = Trim(c_phot)//" "//Trim(c_phot)
   Lnames(i_zaehl) = Fnames(i_zaehl)
   Lnames(i_zaehl+1) = Fnames(i_zaehl+1)
   id_d_2(i_zaehl,:) = id_gl
   id_d_2(i_zaehl+1,:) = id_ph
   i_zaehl = i_zaehl + 2

   gT = gT_S0(i1)
   BR(1:250) = BR_S0(i1,:)
   gP(1:250) = gP_S0(i1,:)

   Write(io_L,200) id_S0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   
   Fnames(i_zaehl) = "W^+ W^-*" 
   Fnames(i_zaehl+1) = "W^- W^+*"
   Fnames(i_zaehl+2) = "Z Z^*" 

   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)
   If ((BR(i_zaehl).Gt.BrMin).Or.(BR(i_zaehl+2).Gt.BrMin)) Then
    Write(io_L,100) "# writing decays into V V* as 3-body decays"
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    If (BR(i_zaehl).Gt.BrMin) Then
     Do i2=1,3
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,id_w,id_l(i2),-id_nu(i2) &
      & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_lm(i2))//" "//Trim(c_nu(i2))//")"
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,-id_w,-id_l(i2),id_nu(i2) &
      & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_lp(i2))//" "//Trim(c_nu(i2))//")"
     End Do
     BRtot = 0.5_dp * Sum(BrWqq)
     Do i2=1,3
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(1),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & ,3,-id_w,id_u(1),-id_d(i2) &
        & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(2),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
       & ,3,-id_w,id_u(2),-id_d(i2) &
       & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
     End Do
    End If
    i_zaehl = i_zaehl + 2
    If (BR(i_zaehl).Gt.BrMin) Then 
     Write(io_L,202) BrZinv*BR(i_zaehl),3,id_Z,id_nu(1),-id_nu(1), &
        & Trim(c_m)//" -> Z0 nu_i nu_i)"
     Do i2=1,3  
      Write(io_L,202) BrZll(i2)*BR(i_zaehl),3,id_Z,id_l(i2),-id_l(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_lm(i2))//" "//Trim(c_lp(i2))//")"
     End Do
     Do i2=1,3  
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_d(i2),-id_d(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_d(i2))//" "//Trim(c_d(i2))//")"
     End Do
     Do i2=1,2 
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_u(i2),-id_u(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_u(i2))//" "//Trim(c_u(i2))//")"
     End Do
    End If
   End If

  End Do


  !----------------------
  ! neutral pseudoscalar
  !----------------------
  Do i1=1,n_p0
   c_m = c_p0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 5 - n_c
     Do i3 = i2, 5 - n_c
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,5 - n_c
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If ((n_c-2).Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If ((n_c-2).Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If ((n_c-2).Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     i2 = 16
    Else If ((n_c-2).Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     i2 = 21
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n
    Do i3=i2,n_n
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c
    Do i3=i2, n_c
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do

   Do i2=1,n_s0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,i1-1
    Do i3=1,n_s0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   gT = gT_P0(i1+1)
   BR(1:250) = BR_P0(i1+1,:)
   gP(1:250) = gP_P0(i1+1,:)

   Write(io_L,200) id_p0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   If (BR(i_zaehl+3).Gt.BrMin)  Write(io_L,201) BR(i_zaehl+3),2,id_gl  &
           &  ,id_gl,Trim(c_m)//" -> g g)"
 
   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)

  End Do


  !-----------------
  ! charged scalars
  !-----------------
  Do i1=1,n_spm
   c_m = c_sm(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    !-----------------------------------------------------------------
    ! leptons, if RP is violated, these BRs are included in the
    ! chargino/neutralino final  states
    !-----------------------------------------------------------------
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do

    !--------
    ! quarks
    !--------
    Do i2 = 1, 3
     Do i3 = 1,3
      Lnames(i_zaehl) = Trim(c_u(i3))//" "//Trim(c_d(i2))
      id_d_2(i_zaehl,1) = id_d(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

    !-----------------------------------------------------------------
    ! sleptons, in case of RP violation -> S0 S-, P0 S- final states
    !-----------------------------------------------------------------
    Do i2 = 1,2*(5-n_c)
     Do i3 = 1,5-n_c
      Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_snu(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
    !------------------
    ! into squarks
    !------------------
    Do i2 = 1,6
     Do i3 = 1,6
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

   Else ! no Generation Mixing
    Do i2=1,5-n_c
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = id_d(i2)
     id_d_2(i_zaehl,2) = -id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1, 2*(5 - n_c)
     Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu((i2+1)/2))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_snu((i2+1)/2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,6
     If (i2.Le.2) i3=1
     If (i2.Le.4) i3=3
     If (i2.Le.6) i3=5
     Do i4=i3,i3+1
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i4))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i4)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
   End If

   Do i2=1,n_c
    Do i3=1,n_n
     Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = - id_c(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   Do i2=1,n_p0
    Lnames(i_zaehl) = " W- "//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_p0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,n_s0
    Lnames(i_zaehl) = " W- "//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   gT = gT_Spm(i1+1)
   BR(1:200) = BR_Spm(i1+1,:)
   gP(1:200) = gP_Spm(i1+1,:)

   Write(io_L,200) id_sm(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    

  End Do

 100 Format(a)
 200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a)
 201 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x,"# BR(",a)
 202 Format(3x,1P,e16.8,0p,3x,I2,3x,3(i9,1x),2x,"# BR(",a)
 
 end subroutine LH_write_decays

 Subroutine LH_write_ext(io_L, Q, bi)
 Implicit None
  Integer, Intent(in) :: io_L
  Logical, Intent(in) :: bi
  Real(dp), intent(in) :: Q

   Write(io_L,100) "Block EXTPAR  # "
   Write(io_L,104) 0,Q , "# scale Q where the parameters below are defined"
   Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
   Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
   Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
   Write(io_L,104) 11,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t"
   Write(io_L,104) 12,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b"
   Write(io_L,104) 13,Real(A_l(3,3)/y_l(3,3),dp), "# A_l"
   If (bi) Write(io_L,104) 23,Real(mu ,dp), "# mu "
   Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
   Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
   Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
   Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
   Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
   Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
   Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
   Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
   Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
   Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
   Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
   Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
   Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
   Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
   Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

 100 Format(a)
 104 Format(i4,2x,1P,e16.8,2x,a)
      
 End Subroutine LH_write_ext

 Subroutine LH_write_msq(io_L, c_sd, id_sd, c_su, id_su)
 Implicit None
  Integer, Intent(in) :: io_L
  Character(len=4), Intent(out) :: c_sd(6), c_su(6)
  Integer, Intent(out) :: id_sd(6), id_su(6)

  Integer :: i1


  If (GenerationMixing) Then

   id_sd(1) = 1000001
   id_sd(2) = 1000003
   id_sd(3) = 1000005
   id_sd(4) = 2000001
   id_sd(5) = 2000003
   id_sd(6) = 2000005
   Do i1=1,6
    c_sd(i1) = "~d_"//Bu(i1)
    Write(io_L,102) id_sd(i1),msdown(i1),"# "//Trim(c_sd(i1))
   End Do

   id_su(1) = 1000002
   id_su(2) = 1000004
   id_su(3) = 1000006
   id_su(4) = 2000002
   id_su(5) = 2000004
   id_su(6) = 2000006
   Do i1=1,6
    c_su(i1) = "~u_"//Bu(i1)
    Write(io_L,102) id_su(i1),msup(i1),"# "//Trim(c_su(i1))
   End Do

  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,msdown(1),"# ~d_L"
    Write(io_L,102) 2000001,msdown(2),"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,msdown(2),"# ~d_L"
    Write(io_L,102) 2000001,msdown(1),"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,msup(1),"# ~u_L"
    Write(io_L,102) 2000002,msup(2),"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,msup(2),"# ~u_L"
    Write(io_L,102) 2000002,msup(1),"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,msdown(3),"# ~s_L"
    Write(io_L,102) 2000003,msdown(4),"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,msdown(4),"# ~s_L"
    Write(io_L,102) 2000003,msdown(3),"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,msup(3),"# ~c_L"
    Write(io_L,102) 2000004,msup(4),"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,msup(4),"# ~c_L"
    Write(io_L,102) 2000004,msup(3),"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,msdown(5),"# ~b_1"
   Write(io_L,102) 2000005,msdown(6),"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,msup(5),"# ~t_1"
   Write(io_L,102) 2000006,msup(6),"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
  
 102 Format(1x,i9,3x,1P,e16.8,2x,a)

 End Subroutine LH_write_msq

 Subroutine LH_write_rsq(io_L, RUsq_ckm, RDsq_ckm)
 Implicit None
  Integer, Intent(in) :: io_L
  Complex(dp), Intent(in) :: RUsq_ckm(6,6), RDsq_ckm(6,6)
  Integer :: i1, i2

  If (generationmixing) Then
   Write(io_L,100) "Block USQmix  # u-sqark mixing matrix"
   Do i1=1,6
    Do i2=1,6
     Write(io_L,105) i1,i2,Real(RUsq_ckm(i1,i2),dp) &
                   & ,"# R_Su("//bu(i1)//","//bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag((RUsq_ckm)))).Ne.0._dp) Then
    Write(io_L,100) "Block IMUSQmix  # imaginiary parts of u-sqark mixing matrix"
    Do i1=1,6
     Do i2=1,6
      Write(io_L,105) i1,i2,Aimag(RUsq_ckm(i1,i2)) &
                    & ,"# Im[R_Su("//bu(i1)//","//bu(i2)//")]"
     End Do
    End Do
   End If
   Write(io_L,100) "Block DSQmix  # d-squark mixing matrix"
   Do i1=1,6
    Do i2=1,6
     Write(io_L,105) i1,i2,Real(RDsq_ckm(i1,i2),dp) &
                   & ,"# R_Sd("//bu(i1)//","//bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag((RDsq_ckm)))).Ne.0._dp) Then
    Write(io_L,100) "Block IMDSQmix  # imaginiary parts of d-sqark mixing matrix"
    Do i1=1,6
     Do i2=1,6
      Write(io_L,105) i1,i2,Aimag(Rdsq_ckm(i1,i2)) &
                    & ,"# Im[R_Sd("//bu(i1)//","//bu(i2)//")]"
     End Do
    End Do
   End If

  Else ! .not.GenerationMixing

   Write(io_L,100) "Block stopmix  # stop mixing matrix"
   Write(io_L,105) 1,1,Real(RUsq_ckm(5,5),dp),"# R_st(1,1)"
   Write(io_L,105) 1,2,Real(RUsq_ckm(5,6),dp),"# R_st(1,2)"
   Write(io_L,105) 2,1,Real(RUsq_ckm(6,5),dp),"# R_st(2,1)"
   Write(io_L,105) 2,2,Real(RUsq_ckm(6,6),dp),"# R_st(2,2)"
   Write(io_L,100) "Block sbotmix  # sbottom mixing matrix"
   Write(io_L,105) 1,1,Real(RDsq_ckm(5,5),dp),"# R_sb(1,1)"
   Write(io_L,105) 1,2,Real(RDsq_ckm(5,6),dp),"# R_sb(1,2)"
   Write(io_L,105) 2,1,Real(RDsq_ckm(6,5),dp),"# R_sb(2,1)"
   Write(io_L,105) 2,2,Real(RDsq_ckm(6,6),dp),"# R_sb(2,2)"

  End If

 100 Format(a)
 105 Format(1x,2i3,3x,1P,e16.8,3x,a)

 End Subroutine LH_write_rsq

 Subroutine PutUpperCase(name)
 Implicit None
  Character(len=80), Intent(inout) :: name
  Integer :: len=80, i1
  Do i1=1,len
   If (name(i1:i1).Eq."a") name(i1:i1) = "A"
   If (name(i1:i1).Eq."b") name(i1:i1) = "B"
   If (name(i1:i1).Eq."c") name(i1:i1) = "C"
   If (name(i1:i1).Eq."d") name(i1:i1) = "D"
   If (name(i1:i1).Eq."e") name(i1:i1) = "E"
   If (name(i1:i1).Eq."f") name(i1:i1) = "F"
   If (name(i1:i1).Eq."g") name(i1:i1) = "G"
   If (name(i1:i1).Eq."h") name(i1:i1) = "H"
   If (name(i1:i1).Eq."i") name(i1:i1) = "I"
   If (name(i1:i1).Eq."j") name(i1:i1) = "J"
   If (name(i1:i1).Eq."k") name(i1:i1) = "K"
   If (name(i1:i1).Eq."l") name(i1:i1) = "L"
   If (name(i1:i1).Eq."m") name(i1:i1) = "M"
   If (name(i1:i1).Eq."n") name(i1:i1) = "N"
   If (name(i1:i1).Eq."o") name(i1:i1) = "O"
   If (name(i1:i1).Eq."p") name(i1:i1) = "P"
   If (name(i1:i1).Eq."q") name(i1:i1) = "Q"
   If (name(i1:i1).Eq."r") name(i1:i1) = "R"
   If (name(i1:i1).Eq."s") name(i1:i1) = "S"
   If (name(i1:i1).Eq."t") name(i1:i1) = "T"
   If (name(i1:i1).Eq."u") name(i1:i1) = "U"
   If (name(i1:i1).Eq."v") name(i1:i1) = "V"
   If (name(i1:i1).Eq."w") name(i1:i1) = "W"
   If (name(i1:i1).Eq."x") name(i1:i1) = "X"
   If (name(i1:i1).Eq."y") name(i1:i1) = "Y"
   If (name(i1:i1).Eq."z") name(i1:i1) = "Z"
  End Do
 End Subroutine PutUpperCase


  Subroutine ReadMatrixC(io, nmax, mat, ic, mat_name, kont, fill)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io, ic
   Integer, Intent(in), Optional :: fill
   Complex(dp), Intent(inout) :: mat(nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert   ! , read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) Then
     mat(i1,i2) = Cmplx(0._dp,Aimag(mat(i1,i2)),dp) + wert
     If (Present(fill).And.(i1.Ne.i2)) &
       &  mat(i2,i1) = Cmplx(0._dp, Aimag(mat(i2,i1)), dp) + wert
    Else If (ic.Eq.1) Then
     mat(i1,i2) = Real(mat(i1,i2),dp) + Cmplx(0._dp, wert, dp)
     !-------------------------------------------------------------
     ! if fill==1 -> matrix is hermitian
     ! if fill==2 -> matrix is complex symmetric
     !-------------------------------------------------------------
     If (Present(fill).And.(i1.Ne.i2)) Then
      If (fill.Eq.1) mat(i2,i1) = Real(mat(i2,i1),dp) - Cmplx(0._dp, wert, dp)
      If (fill.Eq.2) mat(i2,i1) = Real(mat(i2,i1),dp) + Cmplx(0._dp, wert, dp)
     End If
    End If

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadMatrixC

  Subroutine ReadMatrixC2(io, nmax, mat, ic, i_in, mat_name, kont, fill)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io, ic, i_in
   Integer, Intent(in), Optional :: fill
   Complex(dp), Intent(inout) :: mat(nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2, i3
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixC2"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, i3, wert   ! , read_line

    !----------------------------------------------------------------------
    ! reshuffling indices if needed as one is fixed by the model used
    ! and, thus, not needed
    !----------------------------------------------------------------------
    If (i_in.eq.1) then
     i1= i2
     i2 = i3
    Else If (i_in.eq.2) then
     i2 = i3
    End if
    
    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) Then
     mat(i1,i2) = Cmplx(0._dp,Aimag(mat(i1,i2)),dp) + wert
     If (Present(fill).And.(i1.Ne.i2)) &
       &  mat(i2,i1) = Cmplx(0._dp, Aimag(mat(i2,i1)), dp) + wert
    Else If (ic.Eq.1) Then
     mat(i1,i2) = Real(mat(i1,i2),dp) + Cmplx(0._dp, wert, dp)
     !-------------------------------------------------------------
     ! if fill==1 -> matrix is hermitian
     ! if fill==2 -> matrix is complex symmetric
     !-------------------------------------------------------------
     If (Present(fill).And.(i1.Ne.i2)) Then
      If (fill.Eq.1) mat(i2,i1) = Real(mat(i2,i1),dp) - Cmplx(0._dp, wert, dp)
      If (fill.Eq.2) mat(i2,i1) = Real(mat(i2,i1),dp) + Cmplx(0._dp, wert, dp)
     End If
    End If

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadMatrixC2

  Subroutine ReadMatrixR(io, nmax, mat, mat_name, kont)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: mat(nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert   !, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -309
     Call AddError(309) 
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -309
     Call AddError(309) 
     Call TerminateProgram()
    End If

    mat(i1,i2) = wert

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadMatrixR
  
  Subroutine ReadTensorC(io, nmax, mat, ic, mat_name, kont, AS)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: mat(nmax, nmax, nmax)
   Integer, Intent(out) :: kont
   Integer, Intent(in), Optional :: AS ! if present, specifies which entries
                                       ! are antisymmetri: 12 -> i1 i2 
                                       !                   13 -> i1 i3
                                       !                   23 -> i2 i3
   Character(len=80) :: read_line
   Integer :: i1, i2, i3
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadTensorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, i3, wert  ! , read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -312
     Call AddError(312)
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -312 
     Call AddError(312)
     Call TerminateProgram()
    End If
    If ((i3.Lt.1).Or.(i3.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i3=",i3
     Iname = Iname - 2
     kont = -312
      Call AddError(312)
    Call TerminateProgram()
    End If

    If (ic.Eq.0) mat(i1,i2,i3) = Cmplx(0._dp, Aimag(mat(i1,i2,i3)), dp) + wert
    If (ic.Eq.1) mat(i1,i2,i3) = mat(i1,i2,i3) + Cmplx(0._dp, wert, dp)

    If (present(as)) then
     If (as.eq.12) mat(i2,i1,i3) = - mat(i1,i2,i3)
     If (as.eq.13) mat(i3,i2,i1) = - mat(i1,i2,i3)
     If (as.eq.23) mat(i1,i3,i2) = - mat(i1,i2,i3)
    End if

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadTensorC
  
  Subroutine ReadVectorC(io, nmax, vec, ic, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert  !, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -310
     Call AddError(310)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) vec(i1) = Cmplx(0._dp, Aimag(vec(i1)), dp) + wert
    If (ic.Eq.1) vec(i1) = Real(vec(i1),dp) + Cmplx(0._dp, wert, dp)

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorC
  
  Subroutine ReadVectorC2(io, nmax, vec, ic, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorC2"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) Trim(read_line)
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert  ! , read_line

    If (i1.Ne.i2) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))
     Write(ErrCan,*) "indices i1 and i2 are expected to be equal",i1,i2
     Write(ErrCan,*) "ignoring this entry with value",wert
     Call AddError(315)
     Cycle
    End If
     
    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -311
     Call AddError(311)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) vec(i1) = Cmplx(0._dp, Aimag(vec(i1)), dp) + wert
    If (ic.Eq.1) vec(i1) = Real(vec(i1),dp) + Cmplx(0._dp, wert, dp)

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorC2
  
  Subroutine ReadVectorR(io, nmax, vec, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert  ! , read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -311
     Call AddError(311)
     Call TerminateProgram()
    End If

    vec(i1) = wert

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorR
  
  Subroutine ReadVectorR2(io, nmax, vec, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorR2"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) trim(read_line)
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert  ! , read_line

    If (i1.ne.i2) then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))
     Write(ErrCan,*) "indices i1 and i2 are expected to be equal",i1,i2
     Write(ErrCan,*) "ignoring this entry with value",wert
     Call AddError(315)
     Cycle
    End If
     
    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -311
     Call AddError(311)
     Call TerminateProgram()
    End If

    vec(i1) = wert

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorR2

 Subroutine SetWriteMinBR(wert)
 !-------------------------------------------------------------------
 ! sets the minimal branching ratio (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  BrMin = wert
 End Subroutine SetWriteMinBR


 Subroutine SetWriteMinSig(wert)
 !-------------------------------------------------------------------
 ! sets the minimal cross section (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  SigMin = wert
 End Subroutine SetWriteMinSig


 Subroutine WriteSPhenoOutputLHA1(io_L, m_GUT)
 Implicit None

  Integer, Intent(in) :: io_L
  real(dp) :: m_GUT

  Integer :: kont, i1,i2,i3,ierr
  Real(dp) :: g, gp
  !------------------------------
  ! masses and mixing angles
  !------------------------------
  Logical :: non_minimal
  Real(dp) :: Yu(3), Yd(3), Q, MaxCont, mat6r(6,6), nr(4,4), mnr(4)
  Complex(dp), Dimension(3,3) :: CKM_Q
  Complex(dp), Dimension(6,6) :: RUsq_ckm, RDsq_ckm
  Character(len=10) :: name, c_snu(3), c_sle(6)     &
      & , c_su(6), c_sd(6)
  Integer :: n_n, ii, jj
  Integer :: id_sp(7), id_S0(6), id_P0(5), id_snu(3), i_zaehl, id_sle(6) &
      & , id_sd(6), id_check(2), id_su(6)
  Integer, Parameter ::id_A0 = 36, id_Hp = 37                            &
      & , id_W = 24, id_Z = 23, id_ph = 22, id_gl = 21                   &
      & , id_l(3) = (/ 11, 13, 15 /), id_nu(3) = (/ 12, 14, 16 /)        &
      & , id_u(3) = (/ 2, 4, 6 /), id_d(3) = (/ 1, 3, 5 /)               &
      & , id_grav = 1000039, id_glu = 1000021
  Logical :: use_flavour_states=.True.
  Real(dp) :: T3, YL, YR, ML, MR, mSf(2), mSf2(2), mSdown(6), mSup(6)    &
      & , mSlepton(6), mSneut(3)
  Complex(dp) :: Rsf(2,2), A, Y
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "WriteSPhenoOutputLHA1"


  id_sp = 0
  id_sp(1) = 37

   id_S0(1) = 25
   id_S0(2) = 35
   id_P0(1) = 36
   n_n = 4

  Q = Sqrt( GetRenormalizationScale() )
  gp = gauge(1)
  g = gauge(2)
  Open(io_L,file="SPheno_1.spc",status="unknown")
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2.beta - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno 3.xx beta"
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) "# in case of problems send email to porod@ific.uv.es"
   Write(io_L,100) "# Created: today"
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   3.xx beta   # version number"

   !--------------------------------------
   ! model information
   !--------------------------------------
   If (HighScaleModel.Eq."mSugra") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If ((HighScaleModel(1:5).Eq."SUGRA").Or.   &
          & (HighScaleModel(1:5).Eq."SEESA")) Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,100) "    3    2    # mSUGRA model + SU(5)"
    Else If (HighScaleModel.Eq."SUGRA_NuR1") Then
     Write(io_L,100) "    3    3    # mSUGRA model + nu_R at a common scale"
     Write(io_L,100) "Block MnuR  # mass scale of the right handed neutrinos"
     Write(io_L,101) 1,MnuR(1),"# m_nu_R      "     
    Else If (HighScaleModel.Eq."SUGRA_NuR") Then
     Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
     Write(io_L,100) "    3    4    # mSUGRA model + three  nu_R, ~nu_R"
     Write(io_L,106) "Block Ynu0 Q=",m_GUT,"# (GUT scale)"
     Do i1=1,3
      Do i2=1,3
       Write(io_L,105) i2,i1,Real(Y_nu_0(i2,i1),dp), &
             & "# Y_(nu,"//bu(i2)//bu(i1)//")"  
      End Do
     End Do
     If (Maxval(Abs(Aimag(Y_nu_0))).Gt.0._dp) Then
      Write(io_L,106) "Block IMYnu0 Q=",m_GUT,"# (GUT scale)"
      Do i1=1,3
       Do i2=1,3
        Write(io_L,105) i2,i1,Aimag(Y_nu_0(i2,i1)) &
             & ,"# Im(Y_(nu,"//bu(i2)//bu(i1)//"))"  
       End Do
      End Do
     End If
    Else If (HighScaleModel.Eq."SEESAW_II") Then
     Write(io_L,100) "    3    5    # mSUGRA model + Higgs triplett"
     Write(io_L,106) "Block YT0 Q=",m_GUT,"# (GUT scale)"
     Do i1=1,3
      Do i2=i1,3
       Write(io_L,105) i2,i1,Real(Y_T_0(i2,i1),dp),"# Y_(T,"//bu(i2)//bu(i1)//")"  
      End Do
     End Do
     If (Maxval(Abs(Aimag(Y_T_0))).Gt.0._dp) Then
      Write(io_L,106) "Block IMYT0 Q=",m_GUT,"# (GUT scale)"
      Do i1=1,3
       Do i2=i1,3
        Write(io_L,105) i2,i1,Aimag(Y_T_0(i2,i1)) &
             & ,"# Im(Y_(T,"//bu(i2)//bu(i1)//"))"  
       End Do
      End Do
     End If
     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,M_H3(1),"# m_H3      "
     Write(io_L,101) 2,Real(lam12_0(1),dp),"# Re(lambda_1)"
     If (Aimag(lam12_0(1)).Ne.0._dp) &
         Write(io_L,101) 3,Aimag(lam12_0(1)),"# Im(lambda_1)"
     Write(io_L,101) 4,Real(lam12_0(2),dp),"# Re(lambda_2)"
     If (Aimag(lam12_0(2)).Ne.0._dp) &
         Write(io_L,101) 5,Aimag(lam12_0(2)),"# Im(lambda_2)"
     If (Fifteen_plet) Then
      Write(io_L,100) "    6    1               # using RGEs for 15-plet"
     Else
      Write(io_L,100) "    6    0               # using RGEs for 3-plet"
     End If
      
    End If
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,101) 7,Max(M_SO_10,m_GUT),"# SO(10) scale"
     Write(io_L,101) 8,D_SO_10,"# D-terms at SO(10) scale"
    End If
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If (HighScaleModel.Eq."AMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    3    # mAMSB model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,m0_amsb,"# M_0" 
    Write(io_L,101) 2,m_32,"# m_3/2, gravitino mass"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

    
   Else
    non_minimal = .True.
    Write(io_L,100) "# Either the general MSSM or a model has been used"
    Write(io_L,100) &
      & "# which has not yet been implemented in the LesHouches standard"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
   End If 

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"

   !----------------------------------------------------------------
   ! in the case of GenerationMixing all parameters and mixings
   ! are given in the SuperCKM basis for squarks
   !----------------------------------------------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
    Call Switch_to_superCKM(Y_d, Y_u, A_d, A_u, M2_D, M2_Q, M2_U         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .False. &
              &, RSdown, RSup, Rdsq_ckm, RUsq_ckm, CKM_Q, Yd, Yu )

   Else ! .non.GenerationMixing

    Do i1=1,3
     Yu(i1) = Real(Y_u(i1,i1),dp)
     Yd(i1) = Real(Y_d(i1,i1),dp)
    End Do
    Ad_sckm = A_d
    Au_sckm = A_u

    M2D_SCKM = M2_D
    M2U_SCKM = M2_U
    M2Q_SCKM = M2_Q

    RUsq_ckm = RSup
    RDsq_ckm = RSdown

   End If

   If (non_minimal) Then
    Write(io_L,100) "Block EXTPAR  # "
    Write(io_L,104) 0,Q , "# scale Q where the parameters below are defined"
    Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
    Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
    Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
    Write(io_L,104) 11,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t"
    Write(io_L,104) 12,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b"
    Write(io_L,104) 13,Real(A_l(3,3)/y_l(3,3),dp), "# A_l"
    Write(io_L,104) 23,Real(mu ,dp), "# mu "
    Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
    Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
    Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
    Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
    Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
    Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
    Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
    Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
    Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
    Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
    Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
    Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
    Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
    Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
    Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"
   End If
      
! couplings
  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,gauge(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,gauge(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,gauge(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  If (GenerationMixing) Then 
                             
   Write(io_L,106) "Block VCKM Q=",Q,"# Re(CKM) at the SUSY scale"
   Do i1=1,3
    Do i2=1,3
     Write(io_L,107) i1,i2,Real(CKM_Q(i1,i2),dp),"# Re(V_"//Bu(i1)//Bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag(CKM_Q))).Gt.0._dp) Then
    Write(io_L,106) "Block IMVCKM Q=",Q,"# Im(CKM) at the SUSY scale"
    Do i1=1,3
     Do i2=1,3
      Write(io_L,107) i1,i2,Aimag(CKM_Q(i1,i2)),"# Im(V_"//Bu(i1)//Bu(i2)//")"
     End Do
    End Do
   End If

  End If ! generationmixing

  Write(io_L,106) "Block Ye Q=",Q,"# (SUSY scale)"

  ierr = 0
  If (GenerationMixing) Then
   !-------------------------------------------
   ! check if any off-diagonal term is non-zero
   !-------------------------------------------
   Do i1=1,3
    Do i2=1,3
     If ((i1.Ne.i2).And.(Abs(Y_l(i2,i1)).Ne.0._dp)) ierr = ierr + 1
    End Do
   End Do
  End If

  If (ierr.Ne.0) Then
   Do i1=1,3
    Do i2=1,3
     Write(io_L,105) i2,i1,Real(Y_l(i2,i1),dp),"# Y_(l,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
  Else 
   Write(io_L,107) 1,1,Real(y_l(1,1),dp), "# Y_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(y_l(2,2),dp), "# Y_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(y_l(3,3),dp), "# Y_tau(Q)^DRbar"
  End If
  If (Maxval(Abs(Aimag(Y_l))).Gt.0._dp) Then
   Write(io_L,106) "Block IMYe Q=",Q,"# (SUSY scale)"
   If (GenerationMixing) Then
    Do i1=1,3
      Do i2=1,3
      Write(io_L,105) i2,i1,Aimag(Y_l(i2,i1)),"# Im(Y_(l,"//bu(i2)//bu(i1)//"))"
     End Do
    End Do
   Else 
    Write(io_L,107) 1,1,Aimag(y_l(1,1)), "# Im(Y_e)(Q)^DRbar"
    Write(io_L,107) 2,2,Aimag(y_l(2,2)), "# Im(Y_mu)(Q)^DRbar"
    Write(io_L,107) 3,3,Aimag(y_l(3,3)), "# Im(Y_tau)(Q)^DRbar"
   End If
  End If

   Write(io_L,106) "Block Au Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(Au_sckm(1,1)/y_u(1,1),dp), "# A_u(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Au_sckm(2,2)/y_u(2,2),dp), "# A_c(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t(Q)^DRbar"

   Write(io_L,106) "Block Ad Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(Ad_sckm(1,1)/y_d(1,1),dp), "# A_d(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Ad_sckm(2,2)/y_d(2,2),dp), "# A_s(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b(Q)^DRbar"

   Write(io_L,106) "Block Ae Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(A_l(1,1)/y_l(1,1),dp), "# A_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(A_l(2,2)/y_l(2,2),dp), "# A_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(A_l(3,3)/y_l(3,3),dp), "# A_tau(Q)^DRbar"


  Write(io_L,106) "Block MSOFT Q=",Q,"# soft SUSY breaking masses at Q"
  Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
  Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
  Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
  Write(io_L,104) 21,M2_H(1),"# M^2_(H,d)"
  Write(io_L,104) 22,M2_H(2),"# M^2_(H,u)"

  Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
  Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
  Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
  Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
  Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
  Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
  Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
  Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
  Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
  Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
  Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
  Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
  Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
  Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
  Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

  If (GenerationMixing) Then
   Write(io_L,106) "Block MSL2 Q=",Q,"# M^2_L soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2_L(i2,i1),dp)  &
                  & ,"# M^2_(L,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2_L))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSL2 Q=",Q  &
                  & ,"# Im(M^2_L) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2_L(i2,i1)) &
                 &   ,"# Im(M^2_(L,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSE2 Q=",Q,"# M^2_E soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2_E(i2,i1),dp)  &
                  & ,"# M^2_(E,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2_E))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSE2 Q=",Q  &
                  & ,"# Im(M^2_E) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2_E(i2,i1)) &
                 &   ,"# Im(M^2_(E,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSQ2 Q=",Q,"# M^2_Q soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2Q_SCKM(i2,i1),dp) &
                   & ,"# M^2_(Q,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2Q_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSQ2 Q=",Q  &
                  & ,"# Im(M^2_Q) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2Q_SCKM(i2,i1)) &
                 &   ,"# Im(M^2_(Q,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSU2 Q=",Q,"# M^2_U soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2U_SCKM(i1,i2),dp)  &
                  & ,"# M^2_(U,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2U_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSU2 Q=",Q  &
                  & ,"# Im(M^2_U) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2U_SCKM(i1,i2)) &
                 &   ,"# Im(M^2_(U,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSD2 Q=",Q,"# M^2_D soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2D_SCKM(i1,i2),dp)  &
                  & ,"# M^2_(D,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2D_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSD2 Q=",Q  &
                  & ,"# Im(M^2_D) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2D_SCKM(i1,i2)) &
                 &   ,"# Im(M^2_(D,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

  End If

  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,104) 61,Real(h0,dp),"# lambda"
   Write(io_L,104) 62,0.5_dp*Real(lam,dp),"# kappa"
   Write(io_L,104) 63,Real(ao_h0,dp),"# A_lambda"
   Write(io_L,104) 64,Real(Ao_lam,dp),"# A_kappa"
   Write(io_L,104) 65,Real(oosqrt2*vP*h0,dp),"# mu_eff"

  Else If (HighScaleModel.Eq."RPexplicit") Then
   Write(io_L,106) "Block RVKAPPA Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(eps(1),dp),"# epsilon_1"
   Write(io_L,102) 2,Real(eps(2),dp),"# epsilon_2"
   Write(io_L,102) 3,Real(eps(3),dp),"# epsilon_3"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVKAPPA Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(eps(1)),"# Im(epsilon_1)"
    Write(io_L,102) 2,Aimag(eps(2)),"# Im(epsilon_2)"
    Write(io_L,102) 3,Aimag(eps(3)),"# Im(epsilon_3)"
   End If

   Write(io_L,106) "Block RVD Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(Beps(1),dp),"# Re( B_1 epsilon_1)"
   Write(io_L,102) 2,Real(Beps(2),dp),"# Re( B_2 epsilon_2)"
   Write(io_L,102) 3,Real(Beps(3),dp),"# Re( B_3 epsilon_3)"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVD Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(Beps(1)),"# Im( B_1 epsilon_1)"
    Write(io_L,102) 2,Aimag(Beps(2)),"# Im( B_2 epsilon_2)"
    Write(io_L,102) 3,Aimag(Beps(3)),"# Im( B_3 epsilon_3)"
   End If

   Write(io_L,106) "Block RVSNVEV Q=",Q,"# sneutrino vevs at Q"
   Write(io_L,102) 1,vevL(1),"# v_L_1"
   Write(io_L,102) 2,vevL(2),"# v_L_2"
   Write(io_L,102) 3,vevL(3),"# v_L_3"

   If (Maxval(Abs(Rp_lam)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDA Q=",Q,"# lambda_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lam(i1,i2,i3),dp) &
                     & ,"# lambda_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDA Q=",Q,"# Im(lambda_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lam(i1,i2,i3)) &
            & ,"# Im(lambda_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   If (Maxval(Abs(Rp_lamp)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDAP Q=",Q,"# lambda'_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lamp(i1,i2,i3),dp) &
                     & ,"# lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDAP Q=",Q,"# Im(lambda'_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lamp(i1,i2,i3)) &
            & ,"# Im(lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   Write(io_L,100) "Block SPhenoRP  # additional RP parameters"
   Write(io_L,102) 4,Lam_Ex(1),"# Lambda_1 = v_d epsilon_1 + mu v_L1"
   Write(io_L,102) 5,Lam_Ex(2),"# Lambda_2 = v_d epsilon_2 + mu v_L2"
   Write(io_L,102) 6,Lam_Ex(3),"# Lambda_3 = v_d epsilon_3 + mu v_L3"
   Write(io_L,102) 7,(mN7(3)**2-mN7(2)**2)*1.e18_dp,"# m^2_atm [eV^2]"
   Write(io_L,102) 8,(mN7(2)**2-mN7(1)**2)*1.e18_dp,"# m^2_sol [eV^2]"
   Write(io_L,102) 9,Abs(N7(3,6)/N7(3,7))**2,"# tan^2 theta_atm"
   Write(io_L,102) 10,Abs(N7(2,5)/N7(1,5))**2,"# tan^2 theta_sol"
   Write(io_L,102) 11,Abs(N7(3,5))**2,"# U_e3^2"
   Write(io_L,102) 15,vevSM(1),"# v_d"
   Write(io_L,102) 16,vevSM(2),"# v_u"
  End If


   Write(io_L,100) "Block MASS  # Mass spectrum"
   Write(io_L,100) "#   PDG code      mass          particle"
   Write(io_L,102) id_u(2),mf_u(2),"# m_c(m_c), MSbar"
   Write(io_L,102) id_d(3),mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) id_u(3),mf_u(3),"# m_t(pole)"
   Write(io_L,102) id_Z,mZ,"# m_Z(pole)"
   Write(io_L,102) id_W,mW,"# W+"
   If (HighScaleModel.Ne."RPexplicit") &
        & Write(io_L,102) id_l(3),mf_l(3),"# m_tau(pole)"
   If (HighScaleModel.Eq."NMSSM") Then
    Write(io_L,102) 25,S03(1)%m,"# leightest neutral scalar" 
    Write(io_L,102) 35,S03(2)%m,"# second neutral scalar" 
    Write(io_L,102) 45,S03(3)%m,"# third neutral scalar" 
    Write(io_L,102) 36,P03(2)%m,"# leighter pseudoscalar" 
    Write(io_L,102) 46,P03(3)%m,"# heavier pseudoscalar" 
    Write(io_L,102) 37,Spm(2)%m,"# H+"
   Else If (HighScaleModel.Eq."RPexplicit") Then
!    Write(io_L,102) 12,mN7(1),"# lightest neutrino"
!    Write(io_L,102) 14,mN7(2),"# second lightest neutrino"
!    Write(io_L,102) 16,mN7(3),"# heaviest neutrino"
    Write(io_L,102) id_s0(1),S05(1)%m,"# leightest neutral scalar" 
    Write(io_L,102) id_s0(2),S05(2)%m,"# 2nd neutral scalar" 
    Write(io_L,102) id_s0(3),S05(3)%m,"# 3rd neutral scalar" 
    Write(io_L,102) id_s0(4),S05(4)%m,"# 4th neutral scalar" 
    Write(io_L,102) id_s0(5),S05(5)%m,"# 5th neutral scalar" 
    Write(io_L,102) id_P0(1),P05(2)%m,"# leightest pseudoscalar" 
    Write(io_L,102) id_P0(2),P05(3)%m,"# 2nd pseudoscalar" 
    Write(io_L,102) id_P0(3),P05(4)%m,"# 3rd pseudoscalar" 
    Write(io_L,102) id_p0(4),P05(5)%m,"# 4th pseudoscalar"
    Write(io_L,102) Abs(id_sp(1)),Spm8(2)%m,"# leightest charged scalar" 
    Write(io_L,102) Abs(id_sp(2)),Spm8(3)%m,"# 2nd charged scalar" 
    Write(io_L,102) Abs(id_sp(3)),Spm8(4)%m,"# 3rd charged scalar" 
    Write(io_L,102) Abs(id_sp(4)),Spm8(5)%m,"# 4th charged scalar" 
    Write(io_L,102) Abs(id_sp(5)),Spm8(6)%m,"# 5th charged scalar" 
    Write(io_L,102) Abs(id_sp(6)),Spm8(7)%m,"# 6th charged scalar" 
    Write(io_L,102) Abs(id_sp(7)),Spm8(8)%m,"# 7th charged scalar" 
    
   Else
    Write(io_L,102) 25,S0(1)%m,"# h0" 
    Write(io_L,102) 35,S0(2)%m,"# H0" 
    Write(io_L,102) 36,P0(2)%m,"# A0" 
    Write(io_L,102) 37,Spm(2)%m,"# H+"
   End If
! squarks

  mSdown = Sdown%m
  mSup = Sup%m

  If (GenerationMixing) Then
   If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
     i_zaehl = 1
     mat6R = Abs(RDsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_sd(ii) = 1000001
       c_sd(ii) = "~d_L"
      Case(2)
       id_sd(ii) = 1000003
       c_sd(ii) = "~s_L-"
      Case(4)
       id_sd(ii) = 2000001
       c_sd(ii) = "~d_R"
      Case(5)
       id_sd(ii) = 2000003
       c_sd(ii) = "~s_R"
      Case default
       id_sd(ii) = 1000000 * i_zaehl + 5
       c_sd(ii) = "~b_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of sbottoms
      If (id_sd(ii).Eq.1000005)  id_check(1) = ii
      If (id_sd(ii).Eq.2000005)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_sd(ii) = 1000005
      c_sd(ii) = "~b_1"
      ii = id_check(1)  ! the heavier one
      id_sd(ii) = 2000005
      c_sd(ii) = "~b_2"
     End If
     Do ii=1,6
      Write(io_L,102) id_sd(ii),msdown(ii),"# "//Trim(c_sd(ii))
     End Do

     i_zaehl = 1
     mat6R = Abs(RUsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_su(ii) = 1000002
       c_su(ii) = "~u_L"
      Case(2)
       id_su(ii) = 1000004
       c_su(ii) = "~c_L-"
      Case(4)
       id_su(ii) = 2000002
       c_su(ii) = "~u_R"
      Case(5)
       id_su(ii) = 2000004
       c_su(ii) = "~c_R"
      Case default
       id_su(ii) = 1000000 * i_zaehl + 6
       c_su(ii) = "~t_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of stops
      If (id_su(ii).Eq.1000006)  id_check(1) = ii
      If (id_su(ii).Eq.2000006)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_su(ii) = 1000006
      c_su(ii) = "~t_1"
      ii = id_check(1)  ! the heavier one
      id_su(ii) = 2000006
      c_su(ii) = "~t_2"
     End If
     Do ii=1,6
      Write(io_L,102) id_su(ii),msup(ii),"# "//Trim(c_su(ii))
     End Do
   
    Else ! use mass ordering

     id_sd(1) = 1000001
     id_sd(2) = 1000003
     id_sd(3) = 1000005
     id_sd(4) = 2000001
     id_sd(5) = 2000003
     id_sd(6) = 2000005
     Do i1=1,6
      c_sd(i1) = "~d_"//Bu(i1)
      Write(io_L,102) id_sd(i1),msdown(i1),"# "//Trim(c_sd(i1))
     End Do

     id_su(1) = 1000002
     id_su(2) = 1000004
     id_su(3) = 1000006
     id_su(4) = 2000002
     id_su(5) = 2000004
     id_su(6) = 2000006
     Do i1=1,6
      c_su(i1) = "~u_"//Bu(i1)
      Write(io_L,102) id_su(i1),msup(i1),"# "//Trim(c_su(i1))
     End Do
    End If
  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,msdown(1),"# ~d_L"
    Write(io_L,102) 2000001,msdown(2),"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,msdown(2),"# ~d_L"
    Write(io_L,102) 2000001,msdown(1),"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,msup(1),"# ~u_L"
    Write(io_L,102) 2000002,msup(2),"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,msup(2),"# ~u_L"
    Write(io_L,102) 2000002,msup(1),"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,msdown(3),"# ~s_L"
    Write(io_L,102) 2000003,msdown(4),"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,msdown(4),"# ~s_L"
    Write(io_L,102) 2000003,msdown(3),"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,msup(3),"# ~c_L"
    Write(io_L,102) 2000004,msup(4),"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,msup(4),"# ~c_L"
    Write(io_L,102) 2000004,msup(3),"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,msdown(5),"# ~b_1"
   Write(io_L,102) 2000005,msdown(6),"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,msup(5),"# ~t_1"
   Write(io_L,102) 2000006,msup(6),"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
   
   If (HighScaleModel.Ne."RPexplicit") Then

! sleptons
    mSlepton = Slepton%m
    mSneut = Sneut%m
    If (GenerationMixing) Then
     If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
 
      Do i1=1,3
       MaxCont = Maxval(Abs(Rsneut(i1,:))**2)
       If (Abs(Rsneut(i1,1))**2.Eq.maxCont) Then
        id_snu(i1) = 1000012
        c_snu(i1) = "~nu_eL"
       Else If (Abs(Rsneut(i1,2))**2.Eq.maxCont) Then
        id_snu(i1) = 1000014
        c_snu(i1) = "~nu_muL"
       Else
        id_snu(i1) = 1000016
        c_snu(i1) = "~nu_tauL"
       End If
       Write(io_L,102) id_snu(i1),msneut(i1),"# "//Trim(c_snu(i1))
      End Do

      i_zaehl = 1
      mat6R = Abs(RSlepton)
      Do i1=1,6
       MaxCont = Maxval(mat6R)
       Call FindPosition(6, mat6R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_sle(ii) = 1000011
        c_sle(ii) = "~e_L-"
       Case(2)
        id_sle(ii) = 1000013
        c_sle(ii) = "~mu_L-"
       Case(4)
        id_sle(ii) = 2000011
        c_sle(ii) = "~e_R-"
       Case(5)
        id_sle(ii) = 2000013
        c_sle(ii) = "~mu_R-"
       Case default
        id_sle(ii) = 1000000 * i_zaehl + 15
        c_sle(ii) = "~tau_"//bu(i_zaehl)//"-"
        i_zaehl = I_zaehl + 1
       End Select
       mat6R(ii,:) = 0._dp
       mat6R(:,jj) = 0._dp
      End Do
      Do ii=1,6 ! check ordering of staus
       If (id_sle(ii).Eq.1000015)  id_check(1) = ii
       If (id_sle(ii).Eq.2000015)  id_check(2) = ii
      End Do
      If (id_check(1).Gt.id_check(2)) Then ! switch ordering
       ii = id_check(2)  ! the lighter one
       id_sle(ii) = 1000015
       c_sle(ii) = "~tau_1-"
       ii = id_check(1)  ! the heavier one
       id_sle(ii) = 2000015
       c_sle(ii) = "~tau_2-"
      End If
      Do ii=1,6
       Write(io_L,102) id_sle(ii),mslepton(ii),"# "//Trim(c_sle(ii))
      End Do

     Else ! using mass ordering

      id_snu(1) = 1000012
      id_snu(2) = 1000014
      id_snu(3) = 1000016
      Do i1=1,3
       c_snu(i1) = "~nu_"//Bu(i1)
       Write(io_L,102) id_snu(i1),mSneut(i1),"# "//Trim(c_snu(i1))
      End Do

      id_sle(1) = 1000011
      id_sle(2) = 1000013
      id_sle(3) = 1000015
      id_sle(4) = 2000011
      id_sle(5) = 2000013
      id_sle(6) = 2000015
      Do i1=1,6
       c_sle(i1) = "~l_"//Bu(i1)
       Write(io_L,102) id_sle(i1),mSlepton(i1),"# "//Trim(c_sle(i1))
      End Do
     End If

    Else ! .not.GenerationMixing

     id_snu = (/ 1000012, 1000014, 1000016 /)
     If (Abs(rslepton(1,1)).Gt.0.5_dp) Then
      Write(io_L,102) 1000011,mslepton(1),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(2),"# ~e_R-"
      id_sle(1) = 1000011
      id_sle(2) = 2000011
      c_sle(1) = "~e_L-"
      c_sle(2) = "~e_R-"
     Else
      Write(io_L,102) 1000011,mslepton(2),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(1),"# ~e_R-"
      id_sle(2) = 1000011
      id_sle(1) = 2000011
      c_sle(2) = "~e_L-"
      c_sle(1) = "~e_R-"
     End If
     Write(io_L,102) 1000012,msneut(1),"# ~nu_eL"
     c_snu(1) = "~nu_eL"
     If (Abs(rslepton(3,3)).Gt.0.5_dp) Then
      Write(io_L,102) 1000013,mslepton(3),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(4),"# ~mu_R-"
      id_sle(3) = 1000013
      id_sle(4) = 2000013
      c_sle(3) = "~mu_L-"
      c_sle(4) = "~mu_R-"
     Else
      Write(io_L,102) 1000013,mslepton(4),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(3),"# ~mu_R-"
      id_sle(4) = 1000013
      id_sle(3) = 2000013
      c_sle(4) = "~mu_L-"
      c_sle(3) = "~mu_R-"
     End If
     Write(io_L,102) 1000014,msneut(2),"# ~nu_muL"
     c_snu(2) = "~nu_muL"
     Write(io_L,102) 1000015,mslepton(5),"# ~tau_1-"
     Write(io_L,102) 2000015,mslepton(6),"# ~tau_2-"
     id_sle(5) = 1000015
     id_sle(6) = 2000015
     c_sle(5) = "~tau_1-"
     c_sle(6) = "~tau_2-"
     Write(io_L,102) 1000016,msneut(3),"# ~nu_tauL"
     c_snu(3) = "~nu_tauL"
    End If
   End If ! GenerationMixing

 ! gauginos/higgsinos
   Write(io_L,102) 1000021,glu%m,"# ~g"
   ! checking for negative sign

    Do i1=1,4
     If (Sum(Abs(Real(N(i1,:)))).Lt.0.1_dp) Then
      mNr(i1) = - Chi0(i1)%m
      nr(i1,1:4) = Aimag(n(i1,:))
     Else   
      mNr(i1) =  Chi0(i1)%m
      nr(i1,1:4) = Real(n(i1,:),dp)
     End If
    End Do

    Write(io_L,102) 1000022,mnr(1),"# ~chi_10" 
    Write(io_L,102) 1000023,mnr(2),"# ~chi_20" 
    Write(io_L,102) 1000025,mnr(3),"# ~chi_30" 
    Write(io_L,102) 1000035,mnr(4),"# ~chi_40" 
    Write(io_L,102) 1000024,ChiPm(1)%m,"# ~chi_1+" 
    Write(io_L,102) 1000037,ChiPm(2)%m,"# ~chi_2+"

   If (Maxval(MnuR).Gt.0._dp) Then
    Write(io_L,100) "# masses of right handed neutrinos"
    Write(io_L,100) "Block MnuR"
    Write(io_L,102) 1,mnur(1),"# m_nu_R_1"
    Write(io_L,102) 2,mnur(2),"# m_nu_R_2"
    Write(io_L,102) 3,mnur(3),"# m_nu_R_3"
   End If

! Mixing matrices 
  Write(io_L,100) "# Higgs mixing"
   Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
   Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
   Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
   Write(io_L,104) 1,Real(mu,dp),"# mu"
   Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
   Write(io_L,104) 3,vev_Q,"# v(Q)"
   Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"

!--------------------
! add correct order
!--------------------


   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2Q_sckm(3,3)
   Mr = M2U_sckm(3,3)
   A = Au_sckm(3,3)
   Y = Yu(3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
   Write(io_L,100) "Block stopmix  # stop mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_st(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_st(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_st(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_st(2,2)"

   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2Q_sckm(3,3)
   Mr = M2d_sckm(3,3)
   A = Ad_sckm(3,3)
   Y = Yd(3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
   Write(io_L,100) "Block sbotmix  # sbottom mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_sb(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_sb(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_sb(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_sb(2,2)"

   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
    Ml = M2_L(3,3)
    Mr = M2_E(3,3)
    A = A_l(3,3)
    Y = Y_l(3,3)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    Write(io_L,100) "Block staumix  # stau mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_sta(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_sta(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_sta(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_sta(2,2)"

   Write(io_L,100) "Block Nmix  # neutralino mixing matrix"
   Do i1=1,n_n
    Do i2=1,n_n
     name = "# N("//Char(48+i1)//","//Char(48+i2)//")"
     Write(io_L,105) i1,i2,nr(i1,i2),name
    End Do
   End Do

   Write(io_L,100) "Block Umix  # chargino U mixing matrix"
   Write(io_L,105) 1,1,Real(U(1,1),dp),"# U(1,1)"
   Write(io_L,105) 1,2,Real(U(1,2),dp),"# U(1,2)"
   Write(io_L,105) 2,1,Real(U(2,1),dp),"# U(2,1)"
   Write(io_L,105) 2,2,Real(U(2,2),dp),"# U(2,2)"
   Write(io_L,100) "Block Vmix  # chargino V mixing matrix"
   Write(io_L,105) 1,1,Real(V(1,1),dp),"# V(1,1)"
   Write(io_L,105) 1,2,Real(V(1,2),dp),"# V(1,2)"
   Write(io_L,105) 2,1,Real(V(2,1),dp),"# V(2,1)"
   Write(io_L,105) 2,2,Real(V(2,2),dp),"# V(2,2)"

100 Format(a)
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a)
103 Format(a13,1P,e16.8,2x,a)
104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,3x,a)
108 Format(9x,1P,E16.8,0P,3x,a)
109 Format(1x,3i3,3x,1P,e16.8,3x,a)
  Close(io_L)
  

!   if (f_c.lt.10) then
!    f_n = "omg000"//bu(f_c)//".out"
!   else if (f_c.lt.100) then
!    k1 = f_c/10
!    f_n = "omg00"//bu(k1)//bu(f_c-k1*10)//".out"
!   else if (f_c.lt.1000) then
!    k1 = f_c/100
!    k2 = (f_c - k1 * 100)/10
!    f_n = "omg0"//bu(k1)//bu(k2)//bu(f_c-k1*100-k2*10)//".out"
!   else
!    k1 = f_c/1000
!    k2 = (f_c - k1 * 1000)/100
!    k3 = (f_c - k1 * 1000 - k2*100)/10
!    f_n = "omg"//bu(k1)//bu(k2)//bu(k3)//bu(f_c-k1*1000-k2*100-k3*10)//".out"
!   end if   

  Iname = Iname - 1
  
 End Subroutine WriteSPhenoOutputLHA1


 Subroutine Write_HiggsBoundsInput(io, RS0, RP0)
 !---------------------------------------------------------
 ! writes information for HiggsBounds, at the moment only
 ! for the MSSM
 ! written by Werner Porod, 29.04.2011
 !---------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: RS0(:,:), RP0(:,:)

  Real(dp) :: wert, wert1, wert2, wert3

  Iname = Iname + 1
  NameOfUnit(Iname) = "Write_HiggsBoundsInput"

  Write(io,200) 6,2.43_dp,"top"
  Write(io,100) "#    BR                NDA      ID1      ID2"
  Write(io,401) 1.0_dp,2,5,24,"t","b","W"
  
  Write(io,100) "Block HiggsBoundsInputHiggsCouplingsFermions"
  Write(io,100) "# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2"
  wert1 = (RS0(1,1) / RP0(1,1))**2
  Write(io,101) wert1,0.0_dp,3,25,5,5,"h0-b-b"
  wert2 = (RS0(2,1) / RP0(1,1))**2
  Write(io,101) wert2,0.0_dp,3,35,5,5,"H0-b-b"
  wert3 = (RP0(2,1) / RP0(1,1))**2
  Write(io,101) 0.0_dp,wert3,3,36,5,5,"A0-b-b"
  Write(io,100) "#"
  wert = (RS0(1,2) / RP0(1,2))**2
  Write(io,101) wert,0.0_dp,3,25,6,6,"h0-t-t"
  wert = (RS0(2,2) / RP0(1,2))**2
  Write(io,101) wert,0.0_dp,3,35,6,6,"H0-t-t"
  wert = (RP0(2,2) / RP0(1,2))**2
  Write(io,101) 0.0_dp,wert,3,36,6,6,"A0-t-t"
  Write(io,100) "#"
  Write(io,101) wert1,0.0_dp,3,25,15,15,"h0-tau-tau"
  Write(io,101) wert2,0.0_dp,3,35,15,15,"H0-tau-tau"
  Write(io,101) 0.0_dp,wert3,3,36,15,15,"A0-tau-tau"
  Write(io,100) "#"
  Write(io,100) "Block HiggsBoundsInputHiggsCouplingsBosons"
  wert1 = (RP0(2,2) * RS0(1,1) + RP0(1,2) * RS0(1,2) )**2
  Write(io,102) wert1,3,25,24,24,"h0-W-W"
  wert2 = (RP0(2,2) * RS0(2,1) + RP0(1,2) * RS0(2,2) )**2
  Write(io,102) wert2,3,35,24,24,"H0-W-W"
  Write(io,102) 0.0_dp,3,36,24,24,"A0-W-W"
  Write(io,100) "#"
  Write(io,102) wert1,3,25,23,23,"h0-Z-Z"
  Write(io,102) wert2,3,35,23,23,"H0-Z-Z"
  Write(io,102) 0.0_dp,3,36,23,23,"A0-Z-Z"
  Write(io,100) "#"
  Write(io,102) r_GlGlS0(1),3,25,21,21,"h0-g-g"
  Write(io,102) r_GlGlS0(2),3,35,21,21,"H0-g-g"
  Write(io,102) r_GlGlP0,3,36,21,21,"A0-g-g"
  Write(io,100) "#"
  Write(io,102) 0.0_dp,3,25,25,23,"h0-h0-Z"
  Write(io,102) 0.0_dp,3,35,25,23,"H0-h0-Z"
  wert = (RP0(2,1)*RS0(1,1) - RP0(2,2)*RS0(1,2))**2 
  Write(io,102) wert,3,36,25,23,"A0-h0-Z"
  Write(io,102) 0.0_dp,3,35,35,23,"H0-H0-Z"
  wert = (RP0(2,1)*RS0(2,1) - RP0(2,2)*RS0(2,2))**2 
  Write(io,102) wert,3,36,35,23,"A0-H0-Z"
  Write(io,102) 0.0_dp,3,36,36,23,"A0-A0-Z"
  Write(io,100) "#"
  Write(io,103) 0.0_dp,4,25,21,21,23,"h0-g-g-Z"
  Write(io,103) 0.0_dp,4,35,21,21,23,"H0-g-g-Z"
  Write(io,103) 0.0_dp,4,36,21,21,23,"A0-g-g-Z"

  Iname = Iname - 1

100 Format(a)
101 Format(1P,2x,e16.8,2x,e16.8,0P,5x,4i4,"  # ",a," eff. coupling^2, normalised to SM")
102 Format(1P,2x,e16.8,0P,5x,4i4,"  # ",a," eff. coupling^2, normalised to SM")
103 Format(1P,2x,e16.8,0P,5x,5i4,"  # ",a," eff. coupling^2, normalised to SM")

200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a)
401 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x,"# BR(",a," -> ",a," ",a,")")

 End Subroutine Write_HiggsBoundsInput

 Subroutine WriteComplexMatrix(n, matrix, OnlyDiagonal)
 !---------------------------------------------------------------
 ! simplifies the writing of complex matrices in various places
 ! written by Werner Porod
 ! 15.11.01
 !---------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Complex(dp), Intent(in) :: matrix(:,:)
  Logical, Intent(in) :: OnlyDiagonal

  Integer :: i1, i2, dim_matrix
  Complex(dp), Allocatable :: vec(:)

  dim_matrix = Size( matrix, dim=1)

  If (OnlyDiagonal) Then
   Allocate( vec(dim_matrix) )
   Do i1=1,dim_matrix
    vec(i1) = matrix(i1,i1)
   End Do
   If (Maxval( Abs( Aimag(vec) ) ).Eq.0._dp) Then
    Write(n,*) Real(vec)
   Else
    Write(n,*) vec
   End If

   Deallocate( vec )
    
  Else

   If (Maxval( Abs( Aimag(matrix) )) .Lt.1.e-15_dp) Then
    If (dim_matrix.Le.4) Then
     Do i1=1,dim_matrix
      If (dim_matrix.Eq.2) Write(n,902) Real(matrix(i1,:))
      If (dim_matrix.Eq.3) Write(n,903) Real(matrix(i1,:))
      If (dim_matrix.Eq.4) Write(n,904) Real(matrix(i1,:))
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,901) i1,i2,Real(matrix(i1,i2))
      End Do
     End Do
    End If
   Else
    If (dim_matrix.Eq.2) Then
     Do i1=1,dim_matrix
      Write(n,802) matrix(i1,:)
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,801) i1,i2,matrix(i1,i2)
      End Do
     End Do
    End If
   End If
  End If

 801 Format(t2,2i4,"  (",f16.5,",",f16.5,")")
 802 Format(t2,"(",f16.5,",",f16.5,")     (",f16.5,",",f16.5,")")


 901 Format(t2,2i4,f16.5)
 902 Format(t2,2f16.5)
 903 Format(t2,3f16.5)
 904 Format(t2,4f16.5)

 End Subroutine WriteComplexMatrix


 Subroutine WriteCrossSection(io, name, sigma, sigmin)
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: sigma(:,:), sigmin
  Character(len=*), Intent(in) :: name

  Integer :: len1, len2, i1, i2

  len1= Size(sigma, dim=1)
  len2= Size(sigma, dim=2)

  If (Maxval(sigma).Eq.0._dp) Then
    Write(io,*) " "//name//" : kinematically not possible"
   Else If (Maxval(sigma).Lt.SigMin) Then
    Write(io,*) " "//name//" : cross sections below",SigMin,"fb"
   Else 
    Write(io,*) " "//name
    Do i1=1,len1
     Do i2=i1,len2
      If (Sigma(i1,i2).Gt.SigMin) Write(io,900) i1,i2,Sigma(i1,i2) 
     End Do
    End Do
   End If
  Write(io,*) " "

900 Format(2i4,f16.7," fb")

 End Subroutine WriteCrossSection


 Subroutine WriteCrossSection1(io, name, sigma, sigmin)
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: sigma, sigmin
  Character(len=*), Intent(in) :: name

  If (sigma.Eq.0._dp) Then
    Write(io,*) " "//name//" : kinematically not possible"
   Else If (sigma.Lt.SigMin) Then
    Write(io,*) " "//name//" : cross section below",Real(SigMin),"fb"
   Else 
    Write(io,*) " "//name
    Write(io,900) sigma
   End If
  Write(io,*) " "

900 Format(t9,f16.7," fb")

 End Subroutine WriteCrossSection1

 
 Subroutine WriteDecays2(n, titel, names, gP, BR, gT, prez)
 Implicit None
  Integer, Intent(in) :: n
  Character (len=*), Intent(in) :: titel
  Character (len=30), Intent(in) :: names(:)
  Real(dp), Intent(in) :: gP(:), BR(:), gT, prez

  Integer :: i1, dim1

  If (Sum(gp).Eq.0._dp) Then! case of stable particles
   Write (n,*) titel//" : stable"
    Write (n,*) ' '
   Return 
  End If

  dim1 = Size( gP )

  If (Sum(gP).Lt.1.e-7_dp) Then
  ! precision smaller than smallest partial width -> rescaling
   Write (n,*) titel
   Write(n,*) " Attention, the widths are given in eV!"
   Do i1=1,dim1
    If (BR(i1).Ge.(prez+Tiny(1._dp))) Then
     Write (n,100) names(i1), 1.e9_dp*gP(i1), BR(i1)
    End If
   End Do
   If (gT.Gt.0._dp) Then ! to allow the mixing of 2- and 3-body decays
    Write (n,101) 'Total width :                  ',1.e9_dp*gT
    Write (n,*) ' '
   End If
  Else
   Write (n,*) titel
   Do i1=1,dim1
    If (BR(i1).Ge.(prez+Tiny(1._dp))) Then
     Write (n,100) names(i1), gP(i1), BR(i1)
    End If
   End Do
   If (gT.Gt.0._dp) Then ! to allow the mixing of 2- and 3-body decays
    Write (n,101) 'Total width :                  ',gT
    Write (n,*) ' '
   End If
  End If

 100 Format(t3,a31,2f14.8)
 101 Format(t3,a31,f14.8)

 End Subroutine WriteDecays2



 Subroutine WriteMatrixBlockC(io, len, mat_in, scale, name, com1, com2, sym, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Real(dp), Intent(in) :: scale
  Complex(dp), Intent(in) :: mat_in(len,len)
  Character(len=*), Intent(in) :: name, com1, com2
  Logical, Intent(in), Optional :: tr, sym

  Integer :: i1, i2
  Real(dp) :: mat(len,len)
  Logical :: trans, symmetric

  trans = .False.
  If (Present(tr)) trans = tr
  symmetric = .False.
  If (Present(sym)) symmetric = sym

  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  If (trans) mat = Transpose(mat)

  If (symmetric) Then
   Do i2=1,len
    Do i1=i2,len
     Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do

  Else

   Do i2=1,len
    Do i1=1,len
     Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do
  End If

  mat = Aimag(mat_in)
  If (Maxval(Abs(mat)).ne.0._dp) then
   Write(io,106) "Block IM"//Trim(name)//" Q=",scale,"# "//Trim(com1)
   If (trans) mat = Transpose(mat)

   If (symmetric) Then
    Do i2=1,len
     Do i1=i2,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do

   Else

    Do i2=1,len
     Do i1=1,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do
   End If
  End If

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC


 Subroutine WriteMatrixBlockC2(io, len, mat_in, name, com1, com2, sym, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Complex(dp), Intent(in) :: mat_in(len,len)
  Character(len=*), Intent(in) :: name, com1, com2
  Logical, Intent(in), Optional :: tr, sym

  Integer :: i1, i2
  Real(dp) :: mat(len,len)
  Logical :: trans, symmetric

  trans = .False.
  If (Present(tr)) trans = tr
  symmetric = .False.
  If (Present(sym)) symmetric = sym

  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" # "//Trim(com1)
  If (trans) mat = Transpose(mat)

  If (symmetric) Then
   Do i2=1,len
    Do i1=i2,len
     Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do

  Else

   Do i2=1,len
    Do i1=1,len
     Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do
  End If

  mat = Aimag(mat_in)
  If (Maxval(Abs(mat)).ne.0._dp) then
   Write(io,106) "Block IM"//Trim(name)//" # imaginary parts of "//Trim(com1)
   If (trans) mat = Transpose(mat)

   If (symmetric) Then
    Do i2=1,len
     Do i1=i2,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do

   Else

    Do i2=1,len
     Do i1=1,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do
   End If
  End If

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC2


 Subroutine WriteMatrixBlockC3(io, len, i_in, mat_in, scale, name, com1, com2 &
                            &, sym, tr)
 Implicit None
  Integer, Intent(in) :: io, len, i_in
  Real(dp), Intent(in) :: scale
  Complex(dp), Intent(in) :: mat_in(len,len)
  Character(len=*), Intent(in) :: name, com1, com2
  Logical, Intent(in), Optional :: tr, sym

  Integer :: i1, i2
  Real(dp) :: mat(len,len)
  Logical :: trans, symmetric

  trans = .False.
  If (Present(tr)) trans = tr
  symmetric = .False.
  If (Present(sym)) symmetric = sym

  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  If (trans) mat = Transpose(mat)

  If (symmetric) Then
   Do i2=1,len
    Do i1=i2,len
     If (i_in.Eq.1) Then
      Write(io,103) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     Else If (i_in.Eq.2) Then
      Write(io,104) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     Else If (i_in.Eq.3) Then
      Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End If
    End Do
   End Do

  Else

   Do i2=1,len
    Do i1=1,len
     If (i_in.Eq.1) Then
      Write(io,103) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     Else If (i_in.Eq.2) Then
      Write(io,104) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     Else If (i_in.Eq.3) Then
      Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End If
    End Do
   End Do
  End If

  mat = Aimag(mat_in)
  If (Maxval(Abs(mat)).Ne.0._dp) Then
   Write(io,106) "Block IM"//Trim(name)//" Q=",scale,"# "//Trim(com1)
   If (trans) mat = Transpose(mat)

   If (symmetric) Then
    Do i2=1,len
     Do i1=i2,len
      If (i_in.Eq.1) Then
       Write(io,103) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      Else If (i_in.Eq.2) Then
       Write(io,104) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      Else If (i_in.Eq.3) Then
       Write(io,105) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      End If
     End Do
    End Do

   Else

    Do i2=1,len
     Do i1=1,len
      If (i_in.Eq.1) Then
       Write(io,103) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      Else If (i_in.Eq.2) Then
       Write(io,104) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      Else If (i_in.Eq.3) Then
       Write(io,105) i2,i1,mat(i2,i1) &
            & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
      End If
     End Do
    End Do
   End If
  End If

103 Format(1x,"  1",2i3,3x,1P,e16.8,3x,a)
104 Format(1x,i3,"  1",i3,3x,1P,e16.8,3x,a)
105 Format(1x,2i3,"  1",3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC3


 Subroutine WriteMatrixBlockC3a(io, mat_in, scale, name, com1, com2)
 ! simulates a matrix notation with 3 indices
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: scale
  Complex(dp), Intent(in) :: mat_in
  Character(len=*), Intent(in) :: name, com1, com2

  Real(dp) :: mat


  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  Write(io,103) mat,"# Re["//Trim(com2)//"]"

  mat = Aimag(mat_in)
  If (Abs(mat).Ne.0._dp) Then
   Write(io,106) "Block IM"//Trim(name)//" Q=",scale,"# "//Trim(com1)
   Write(io,103) mat,"# In["//Trim(com2)//"]"
  End If

103 Format(1x,"  1  1  1   ",1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC3a


 Subroutine WriteMatrixBlockR(io, len, mat, scale, name, com1, com2, com3, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Real(dp), Intent(in) :: mat(len,len), scale
  Character(len=*), Intent(in) :: name, com1, com2, com3
  Logical, Intent(in), Optional :: tr

  Integer :: i1, i2
  Logical :: trans

  trans = .False.
  If (Present(tr)) trans = tr

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  Do i2=1,len
   Do i1=1,len
    If (trans) Then
     Write(io,105) i2,i1,mat(i1,i2),"# "//Trim(com2)//bu(i2)//","//bu(i1)//Trim(com3)
    Else
     Write(io,105) i2,i1,mat(i2,i1),"# "//Trim(com2)//bu(i2)//","//bu(i1)//Trim(com3)
    End If

   End Do
  End Do

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockR

 Subroutine WriteMSSMParameters(n, CKM, Yd, Yu, WithComments)
 !-----------------------------------------------------------------------
 ! writes out all MSSM parameters
 ! written by Werner Porod
 ! 15.11.01: writting MSSM part,
 !           WriteBilinear and WriteTrilinear should be used later for 
 !           R-parity violation
 !-----------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: CKM(3,3)
  Integer, Intent(in) :: n
  Logical, Intent(in) :: WithComments
  Real(dp), Intent(in) :: Yd(3), Yu(3)

  Logical :: OnlyDiagonal

  OnlyDiagonal = .Not.GenerationMixing

103 Format(1P,3e16.8)

  If (WithComments) Write(n,*) "       g'             g             g_3"
  Write(n,103) gauge
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_e            Y_mu          Y_tau"
  Write(n,903) Real(Y_l(1,1),dp),Real(Y_l(2,2),dp),Real(Y_l(3,3),dp)
!   Call WriteComplexMatrix(n, Y_l, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_u            Y_c            Y_t"
  Write(n,903) Yu
 903 Format(t2,3f16.5)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_d            Y_s            Y_b"
  Write(n,903) Yd
!   Call WriteComplexMatrix(n, Y_d, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) " CKM matrix"
  If (GenerationMixing) Call WriteComplexMatrix(n, CKM, OnlyDiagonal)
  If (WithComments) Write(n,*) " "


  If (WithComments) Write(n,*) "Gaugino mass parameters"
  If (Maxval(Abs(Aimag(Mi))).Eq.0._dp) Then
   Write(n,103) Real( Mi )
  Else
   Write(n,*) Mi
  End If
  If (WithComments) Write(n,*) " "

  
  If (WithComments) Write(n,*) "tan(beta), mu, B"
  If ((Aimag(mu).Eq.0._dp).And.(Aimag(B).Eq.0._dp)) Then
   Write(n,103) tanb, Real(mu), Real(B)
  Else
   Write(n,*) tanb, mu, B
  End If
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Slepton mass parameters"
  If (WithComments) Write(n,*) "A_l"
  Call WriteComplexMatrix(n, A_l, OnlyDiagonal)
  If (WithComments) Write(n,*) "M2_E"
  Call WriteComplexMatrix(n, M2_E, OnlyDiagonal)
  If (WithComments) Write(n,*) "M2_L"
  Call WriteComplexMatrix(n, M2_L, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Squark mass parameters"
  If (GenerationMixing) Then
   If (WithComments) Write(n,*) "A_d"
   Call WriteComplexMatrix(n, Ad_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "A_u"
   Call WriteComplexMatrix(n, Au_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_D"
   Call WriteComplexMatrix(n, M2D_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_U"
   Call WriteComplexMatrix(n, M2U_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_Q"
   Call WriteComplexMatrix(n, M2Q_sckm, OnlyDiagonal)
  Else
   If (WithComments) Write(n,*) "A_d"
   Call WriteComplexMatrix(n, A_d, OnlyDiagonal)
   If (WithComments) Write(n,*) "A_u"
   Call WriteComplexMatrix(n, A_u, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_D"
   Call WriteComplexMatrix(n, M2_D, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_U"
   Call WriteComplexMatrix(n, M2_U, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_Q"
   Call WriteComplexMatrix(n, M2_Q, OnlyDiagonal)
  End If
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Higgs mass parameters"
  Write(n,102) M2_H
  If (WithComments) Write(n,*) " "

102 Format(1P,2e16.8)

 End Subroutine WriteMSSMParameters


 Subroutine WriteParametersAtQ(Qin, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, tanb, phase_mu            &
          & , GenerationMixing, Qout, delta, io_L)
 !-----------------------------------------------------------------------
 ! runs the parameters given at a scale Qin to the scale Qout, calculates
 ! the masses within the MSSM and prints them out 
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: io_L
  Logical, Intent(in) :: GenerationMixing
  Real(dp), Intent(in) :: Qin, Qout, gauge(3), M2_H(2), tanb, delta
  Complex(dp), Intent(in) :: Mi(3), phase_mu
  Complex(dp), Intent(in), Dimension(3,3) :: Y_l, Y_d, Y_u, A_l, A_d, A_u &
    & , M2_E, M2_L, M2_D, M2_Q, M2_U

  ! local variables
  Integer :: kont, i1, id_su(6), id_sd(6), id_sle(6), id_snu(3)
  Real(dp) :: g2(214), gQ(3), M2_H_Q(2), tanb_Q, tz, dt, old_scale, tanb_Ql
  Complex(dp) :: Mi_Q(3), B, mu
  Complex(dp), Dimension(3,3) :: Y_l_Q, Y_d_Q, Y_u_Q, A_l_Q, A_d_Q, A_u_Q  &
    & , M2_E_Q, M2_L_Q, M2_D_Q, M2_Q_Q, M2_U_Q, Ad_sckm, Au_sckm, M2D_sckm &
    & , M2Q_sckm, M2U_sckm, Al_pmns , M2E_pmns, M2L_pmns, PMNS_Q
  Real(dp) :: mC(2), mC2(2), mN(4), mN2(4), mS0(2), mS02(2), RS0(2,2)       &
     & , mP0(2), mP02(2), RP0(2,2), mSpm(2), mSpm2(2), mSdown(6)          &
     & , mSdown2(6), mSup(6), mSup2(6), mSlepton(6), mSlepton2(6) &
     & , mSneut(3), mSneut2(3), mGlu, Yd_Q(3), Yu_Q(3), Yl_Q(3)
  Complex(dp) :: U(2,2), V(2,2), N(4,4), RSpm(2,2), RSneutrino(3,3) &
     & , phase_glu, CKM_Q(3,3), RSn_pmns(3,3), nr(4,4)
  Complex(dp), Dimension(6,6) ::  RDsq_ckm,  RUsq_ckm, RUsquark, RDsquark &
     & , RSlepton, RSl_pmns
  Integer, Parameter :: n_n=4, n_c=2

  Iname = Iname + 1
  NameOfUnit(Iname) = "WriteParametersAtQ"

  tz = Log(Qin/Qout)
  dt = tz / 50._dp

  B = ZeroC  ! will be calculated from tadpole equaitons
  mu = ZeroC

  Y_u_Q = Transpose(Y_u)
  Y_d_Q = Transpose(Y_d)
  Y_l_Q = Transpose(Y_l)
  A_u_Q = Transpose(A_u)
  A_d_Q = Transpose(A_d)
  A_l_Q = Transpose(A_l)
  Call ParametersToG(gauge, Y_l_Q, Y_d_Q, Y_u_Q, Mi, A_l_Q, A_d_Q, A_u_Q &
          & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2(1:213))
  g2(214) = Log( tanb )

  g2(1) = Sqrt(5._dp/3._dp) * g2(1)
  Call odeint(g2, 214, tz, 0._dp, 0.1_dp * delta, dt, 0._dp, rge214, kont)

  g2(1) = Sqrt(3._dp/5._dp) * g2(1)
  Call GToParameters(g2(1:213), gQ, Y_l_Q, Y_d_Q, Y_u_Q, Mi_Q, A_l_Q, A_d_Q &
            & , A_u_Q, M2_E_Q, M2_L_Q, M2_D_Q, M2_Q_Q, M2_U_Q, M2_H_Q, mu, B)
  Y_u_Q = Transpose(Y_u_Q)
  Y_d_Q = Transpose(Y_d_Q)
  Y_l_Q = Transpose(Y_l_Q)
  A_u_Q = Transpose(A_u_Q)
  A_d_Q = Transpose(A_d_Q)
  A_l_Q = Transpose(A_l_Q)
  tanb_Q = Exp( g2(214) )

  If (Calc_Mass) Then ! do this only when required, extends the standard
   kont = 0
   old_scale = SetRenormalizationScale(Qout**2)

   tanb_in_at_Q = .True.
   Call LoopMassesMSSM(0.1_dp*delta, tanb_Ql, tanb_Q, gQ, Y_l_Q, Y_d_Q, Y_u_Q &
     & , Mi_Q, A_l_Q, A_d_Q, A_u_Q, M2_E_Q, M2_L_Q, M2_D_Q, M2_Q_Q, M2_U_Q    &
     & , M2_H_Q, phase_mu, mu, B, 2, uU_L, uU_R, uD_L, uD_R, uL_L, uL_R       &
     & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm      &
     & , mSpm2, RSpm, mSdown, mSdown2, RDsquark, mSup, mSup2, RUsquark        &
     & , mSlepton, mSlepton2, RSlepton, mSneut, mSneut2, RSneutrino           &
     & , mGlu, phase_glu, kont)
    Do i1=1,4
     If (Sum(Abs(Real(N(i1,:)))).Eq.0._dp) Then
      nr(i1,1:4) = Aimag(n(i1,:))
     Else   
      nr(i1,1:4) = n(i1,:)
     End If
    End Do
    mN = Abs(mN)

  End If

  If (GenerationMixing) Then
    Call Switch_to_superCKM(Y_d_Q, Y_u_Q, A_d_Q, A_u_Q, M2_D_Q, M2_Q_Q, M2_U_Q &
              & , Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .False. &
              & , RDSquark, RUSquark, Rdsq_ckm, RUsq_ckm, CKM_Q, Yd_Q, Yu_Q )
    Call Switch_to_superPMNS(Y_l_Q, id3C, A_l_Q, M2_E_Q, M2_L_Q, Al_pmns  &
        & , M2E_pmns, M2L_pmns, .False., RSlepton, RSneutrino, RSl_pmns   &
        & , RSn_pmns, PMNS_Q, Yl_Q)
  Else ! .non.GenerationMixing

    Do i1=1,3
     Yu_Q(i1) = Real(Y_u_Q(i1,i1),dp)
     Yd_Q(i1) = Real(Y_d_Q(i1,i1),dp)
     Yl_Q(i1) = Real(Y_l_Q(i1,i1),dp)
    End Do
    Al_pmns = A_l_Q
    Ad_sckm = A_d_Q
    Au_sckm = A_u_Q
    M2D_SCKM = M2_D_Q
    M2U_SCKM = M2_U_Q
    M2Q_SCKM = M2_Q_Q
    M2E_pmns = M2_E_Q
    M2L_pmns = M2_L_Q
    RDsq_ckm = RDSquark
    RUsq_ckm = RUSquark
    RSl_pmns = RSlepton
    RSn_pmns = RSneutrino

  End If

! couplings
  Write(io_L,106) "Block gauge Q=",Qout,"# (SUSY scale)"
  Write(io_L,104) 1,gQ(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,gQ(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,gQ(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Qout,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu_Q(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu_Q(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu_Q(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Qout,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd_Q(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd_Q(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd_Q(3), "# Y_b(Q)^DRbar"

  Write(io_L,106) "Block Ye Q=",Qout,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yl_Q(1), "# Y_e(Q)^DRbar"
  Write(io_L,107) 2,2,Yl_Q(2), "# Y_mu(Q)^DRbar"
  Write(io_L,107) 3,3,Yl_Q(3), "# Y_tau(Q)^DRbar"

  Write(io_L,100) "Block EXTPAR  # non-universal input parameters"
  Write(io_L,104) 0,Qout,"# Q"
  If (Calc_Mass) Write(io_L,104) 23,Real(mu,dp),"# Re(mu)"
  Write(io_L,104) 25,tanb_Q,"# tan(beta)"
  If (Calc_Mass) Then
   If (Aimag(mu).Ne.0._dp) Then
    Write(io_L,100) "Block IMEXTPAR  # non-universal input parameters"
    Write(io_L,104) 23,Aimag(mu),"# Im(mu)"
   End If
  End If

  If (GenerationMixing) Then 
                             
   Call WriteMatrixBlockC(io_L,3,CKM_Q,Qout &
                         & ,"VCKM","V_CKM at the SUSY scale","V_(")

   Call WriteMatrixBlockC(io_L,3,PMNS_Q,Qout &
                         & ,"UPMNS","U_PMNS at the SUSY scale","V_(")

  End If ! generationmixing


  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,Au_sckm,Qout,"Tu","SUSY scale","T_(u,",tr=.True.)
   Call WriteMatrixBlockC(io_L,3,Ad_sckm,Qout,"Td","SUSY scale","T_(d,",tr=.True.)
   Call WriteMatrixBlockC(io_L,3,Al_pmns,Qout,"Te","SUSY scale","T_(l,", tr=.True.)

  Else 
   Write(io_L,106) "Block Au Q=",Qout,"# (SUSY scale)"
   If (Abs(Yu_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Au_sckm(1,1)/Yu_Q(1),dp), "# A_u(Q)^DRbar"
   If (Abs(Yu_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Au_sckm(2,2)/Yu_Q(2),dp), "# A_c(Q)^DRbar"
   If (Abs(Yu_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Au_sckm(3,3)/Yu_Q(3),dp), "# A_t(Q)^DRbar"
   If (Maxval(Abs(Aimag(Au_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAu Q=",Qout,"# (SUSY scale)"
    If (Abs(Yu_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Au_sckm(1,1)/Yu_Q(1)), "# Im(A_u)(Q)^DRbar"
    If (Abs(Yu_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Au_sckm(2,2)/Yu_Q(2)), "# Im(A_c)(Q)^DRbar"
    If (Abs(Yu_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Au_sckm(3,3)/Yu_Q(3)), "# Im(A_t)(Q)^DRbar"
   End If

   Write(io_L,106) "Block Ad Q=",Qout,"# (SUSY scale)"
   If (Abs(Yd_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Ad_sckm(1,1)/Yd_Q(1),dp), "# A_d(Q)^DRbar"
   If (Abs(Yd_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Ad_sckm(2,2)/Yd_Q(2),dp), "# A_s(Q)^DRbar"
   If (Abs(Yd_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Ad_sckm(3,3)/Yd_Q(3),dp), "# A_b(Q)^DRbar"
   If (Maxval(Abs(Aimag(Ad_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAd Q=",Qout,"# (SUSY scale)"
    If (Abs(Yd_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Ad_sckm(1,1)/Yd_Q(1)), "# Im(A_d)(Q)^DRbar"
    If (Abs(Yd_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Ad_sckm(2,2)/Yd_Q(2)), "# Im(A_s)(Q)^DRbar"
    If (Abs(Yd_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Ad_sckm(3,3)/Yd_Q(3)), "# Im(A_b)(Q)^DRbar"
   End If

   Write(io_L,106) "Block Ae Q=",Qout,"# (SUSY scale)"
   If (Abs(Yl_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Al_pmns(1,1)/Yl_Q(1),dp), "# A_e(Q)^DRbar"
   If (Abs(Yl_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Al_pmns(2,2)/Yl_Q(2),dp), "# A_mu(Q)^DRbar"
   If (Abs(Yl_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Al_pmns(3,3)/Yl_Q(3),dp), "# A_tau(Q)^DRbar"
   If (Maxval(Abs(Aimag(Al_pmns))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAe Q=",Qout,"# (SUSY scale)"
    If (Abs(Yl_Q(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Al_pmns(1,1)/Yl_Q(1)), "# Im(A_e)(Q)^DRbar"
    If (Abs(Yl_Q(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Al_pmns(2,2)/Yl_Q(2)), "# Im(A_mu)(Q)^DRbar"
    If (Abs(Yl_Q(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Al_pmns(3,3)/Yl_Q(3)), "# Im(A_tau)(Q)^DRbar"
   End If
  End If


  Write(io_L,106) "Block MSOFT Q=",Qout,"# soft SUSY breaking masses at Q"
  Write(io_L,104) 1,Real(Mi_Q(1),dp),"# M_1"
  Write(io_L,104) 2,Real(Mi_Q(2),dp),"# M_2"
  Write(io_L,104) 3,Real(Mi_Q(3),dp),"# M_3"
  Write(io_L,104) 21,M2_H_Q(1),"# M^2_(H,d)"
  Write(io_L,104) 22,M2_H_Q(2),"# M^2_(H,u)"
  Write(io_L,104) 31,Sqrt(Real(M2L_pmns(1,1),dp)),"# M_(L,11)"
  Write(io_L,104) 32,Sqrt(Real(M2L_pmns(2,2),dp)),"# M_(L,22)"
  Write(io_L,104) 33,Sqrt(Real(M2L_pmns(3,3),dp)),"# M_(L,33)"
  Write(io_L,104) 34,Sqrt(Real(M2E_pmns(1,1),dp)),"# M_(E,11)"
  Write(io_L,104) 35,Sqrt(Real(M2E_pmns(2,2),dp)),"# M_(E,22)"
  Write(io_L,104) 36,Sqrt(Real(M2E_pmns(3,3),dp)),"# M_(E,33)"
  Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
  Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
  Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
  Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
  Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
  Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
  Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
  Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
  Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

  If (Maxval(Abs(Aimag(Mi_Q))).Ne.0._dp) Then
   Write(io_L,106) "Block IMMSOFT Q=",Qout &
           & ,"# soft SUSY breaking masses at Q, imaginary parts"
   Write(io_L,104) 1,Aimag(Mi_Q(1)),"# M_1"
   Write(io_L,104) 2,Aimag(Mi_Q(2)),"# M_2"
   Write(io_L,104) 3,Aimag(Mi_Q(3)),"# M_3"
  End If

  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,M2L_pmns,Qout,"MSL2" &
           & ,"M^2_L soft SUSY breaking masses","M^2_(L,")

   Call WriteMatrixBlockC(io_L,3,M2E_pmns,Qout,"MSE2" &
           & ,"M^2_E soft SUSY breaking masses","M^2_(E,")

   Call WriteMatrixBlockC(io_L,3,M2Q_SCKM,Qout,"MSQ2" &
           & ,"M^2_Q soft SUSY breaking masses","M^2_(Q,")

   Call WriteMatrixBlockC(io_L,3,M2U_SCKM,Qout,"MSU2" &
           & ,"M^2_U soft SUSY breaking masses","M^2_(U,")

   Call WriteMatrixBlockC(io_L,3,M2D_SCKM,Qout,"MSD2" &
           & ,"M^2_D soft SUSY breaking masses","M^2_(D,")

  End If

  If (.Not.(Calc_Mass.And.(kont.Eq.0))) Then ! return
   If (Calc_Mass) old_scale = SetRenormalizationScale(old_scale)
   Iname = Iname - 1
   Return
  End If
  ! extends the standard

   Write(io_L,103) "Block Hmix Q=",Qout, "# Higgs mixing parameters"
   Write(io_L,104) 1,Real(mu,dp),"# mu"
   Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
   Write(io_L,104) 3,vev_Q,"# v(Q)"
   Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"
   If (Aimag(mu).Ne.0._dp) Then
    Write(io_L,103) "Block IMHmix Q=",Qout, "# Higgs mixing parameters"
    Write(io_L,104) 1,Aimag(mu),"# Im(mu)"
   End If

   Write(io_L,100) "Block MASS  # Mass spectrum"
   Write(io_L,100) "#   PDG code      mass          particle"
   Write(io_L,102) 25,mS0(1),"# h0" 
   Write(io_L,102) 35,mS0(2),"# H0" 
   Write(io_L,102) 36,mP0(2),"# A0" 
   Write(io_L,102) 37,mSpm(2),"# H+"

   If (GenerationMixing) Then
    id_sd(1) = 1000001
    id_sd(2) = 1000003
    id_sd(3) = 1000005
    id_sd(4) = 2000001
    id_sd(5) = 2000003
    id_sd(6) = 2000005
    Do i1=1,6
     Write(io_L,102) id_sd(i1),mSdown(i1),"#~d_"//Bu(i1) 
    End Do

    id_su(1) = 1000002
    id_su(2) = 1000004
    id_su(3) = 1000006
    id_su(4) = 2000002
    id_su(5) = 2000004
    id_su(6) = 2000006
    Do i1=1,6
     Write(io_L,102) id_su(i1),mSup(i1),"# ~u_"//Bu(i1)
    End Do

    id_snu(1) = 1000012
    id_snu(2) = 1000014
    id_snu(3) = 1000016
    Do i1=1,3
     Write(io_L,102) id_snu(i1),mSneut(i1),"# ~nu_"//Bu(i1)
    End Do

    id_sle(1) = 1000011
    id_sle(2) = 1000013
    id_sle(3) = 1000015
    id_sle(4) = 2000011
    id_sle(5) = 2000013
    id_sle(6) = 2000015
    Do i1=1,6
     Write(io_L,102) id_sle(i1),mSlepton(i1),"# ~l_"//Bu(i1)
    End Do

  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,mSdown(1),"# ~d_L"
    Write(io_L,102) 2000001,mSdown(2),"# ~d_R"
   Else
    Write(io_L,102) 1000001,mSdown(2),"# ~d_L"
    Write(io_L,102) 2000001,mSdown(1),"# ~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,mSup(1),"# ~u_L"
    Write(io_L,102) 2000002,mSup(2),"# ~u_R"
   Else
    Write(io_L,102) 1000002,mSup(2),"# ~u_L"
    Write(io_L,102) 2000002,mSup(1),"# ~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,mSdown(3),"# ~s_L"
    Write(io_L,102) 2000003,mSdown(4),"# ~s_R"
   Else
    Write(io_L,102) 1000003,mSdown(4),"# ~s_L"
    Write(io_L,102) 2000003,mSdown(3),"# ~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,mSup(3),"# ~c_L"
    Write(io_L,102) 2000004,mSup(4),"# ~c_R"
   Else
    Write(io_L,102) 1000004,mSup(4),"# ~c_L"
    Write(io_L,102) 2000004,mSup(3),"# ~c_R"
   End If
   Write(io_L,102) 1000005,mSdown(5),"# ~b_1"
   Write(io_L,102) 2000005,mSdown(6),"# ~b_2"
   Write(io_L,102) 1000006,mSup(5),"# ~t_1"
   Write(io_L,102) 2000006,mSup(6),"# ~t_2"

   If (Abs(Rslepton(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000011,mslepton(1),"# ~e_L-"
    Write(io_L,102) 2000011,mslepton(2),"# ~e_R-"
   Else
    Write(io_L,102) 1000011,mslepton(2),"# ~e_L-"
    Write(io_L,102) 2000011,mslepton(1),"# ~e_R-"
   End If
   Write(io_L,102) 1000012,mSneut(1),"# ~nu_eL"
   If (Abs(Rslepton(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000013,mslepton(3),"# ~mu_L-"
    Write(io_L,102) 2000013,mslepton(4),"# ~mu_R-"
   Else
    Write(io_L,102) 1000013,mslepton(4),"# ~mu_L-"
    Write(io_L,102) 2000013,mslepton(3),"# ~mu_R-"
   End If
   Write(io_L,102) 1000014,mSneut(2),"# ~nu_muL"
   Write(io_L,102) 1000015,mslepton(5),"# ~tau_1-"
   Write(io_L,102) 2000015,mslepton(6),"# ~tau_2-"
   Write(io_L,102) 1000016,mSneut(3),"# ~nu_tauL"
  End If ! GenerationMixing

  Write(io_L,102) 1000021,mGlu,"# ~g"
  Write(io_L,102) 1000022,mN(1),"# ~chi_10" 
  Write(io_L,102) 1000023,mN(2),"# ~chi_20" 
  Write(io_L,102) 1000025,mN(3),"# ~chi_30" 
  Write(io_L,102) 1000035,mN(4),"# ~chi_40" 
  Write(io_L,102) 1000024,mC(1),"# ~chi_1+" 
  Write(io_L,102) 1000037,mC(2),"# ~chi_2+"

  If (generationmixing) Then
   Call WriteMatrixBlockC2(io_L, 6, RUsq_ckm, "USQmix" &
                         &, "u-sqark mixing matrix", "R_Su(")

   Call WriteMatrixBlockC2(io_L, 6, RDsq_ckm, "DSQmix" &
                         &, "d-sqark mixing matrix", "R_Sd(")

!   If (HighScaleModel.Ne."RPexplicit") Then

    Call WriteMatrixBlockC2(io_L, 6, RSl_pmns, "SELmix" &
                         &, "slepton mixing matrix", "R_Sl(")

    Call WriteMatrixBlockC2(io_L, 3, RSn_pmns, "SNUmix" &
                         &, "sneutrino mixing matrix", "R_Sn(")

!   End If

  Else ! .not.GenerationMixing

   Call WriteMatrixBlockC2(io_L, 2, RUsq_ckm(5:6,5:6), "stopmix" &
                         &, "stop mixing matrix", "R_st(")

   Call WriteMatrixBlockC2(io_L, 2, RDsq_ckm(5:6,5:6), "sbotmix" &
                         &, "sbottom mixing matrix", "R_sb(")

!   If (HighScaleModel.Ne."RPexplicit") & 
   Call WriteMatrixBlockC2(io_L, 2, RSl_pmns(5:6,5:6), "staumix" &
                            &, "stau mixing matrix", "R_sta(")
  End If

!  If ((HighScaleModel.Eq."RPexplicit").And.(.Not.l_RP_Pythia)) Then
!   Call WriteMatrixBlockC2(io_L, n_n, nr(1:n_n,1:n_n) &
!                         &, "RVNmix", "/neutrino/neutralino mixing matrix", "N(")
!   Call WriteMatrixBlockC2(io_L, 5, U5, "RVUmix" &
!                         &, "lepton/chargino mixing matrix", "U(")
!   Call WriteMatrixBlockC2(io_L, 5, V5, "RVVmix" &
!                         &, "lepton/chargino mixing matrix", "V(")
!  Else If (HighScaleModel.Eq."RPexplicit")  Then
!   Call WriteMatrixBlockC2(io_L, 4, nr(1:4,1:4), "Nmix" &
!                         &, "neutralino mixing matrix", "N(")
!   Call WriteMatrixBlockC2(io_L, 2, U, "Umix" &
!                         &, "chargino mixing matrix", "U(")
!   Call WriteMatrixBlockC2(io_L, 2, V, "Vmix" &
!                         &, "chargino mixing matrix", "V(")
!  Else   
   Call WriteMatrixBlockC2(io_L, n_n, nr(1:n_n,1:n_n), "Nmix" &
                         &, "neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, n_c, U, "Umix" &
                         &, "chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, n_c, V, "Vmix" &
                         &, "chargino mixing matrix", "V(")
!  End If


  If (Calc_Mass) old_scale = SetRenormalizationScale(old_scale)

  Iname = Iname - 1

100 Format(a)
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a)
103 Format(a13,1P,e16.8,2x,a)
104 Format(i4,2x,1P,e16.8,2x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,3x,a)

 End Subroutine WriteParametersAtQ

 Subroutine WriteRealMatrix(n, matrix, OnlyDiagonal)
 !---------------------------------------------------------------
 ! simplifies the writing of complex matrices in various places
 ! written by Werner Porod
 ! 15.11.01
 !---------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Real(dp), Intent(in) :: matrix(:,:)
  Logical, Intent(in) :: OnlyDiagonal

  Integer :: i1, i2, dim_matrix
  Real(dp), Allocatable :: vec(:)

  dim_matrix = Size( matrix, dim=1)

  If (OnlyDiagonal) Then
   Allocate( vec(dim_matrix) )
   Do i1=1,dim_matrix
    vec(i1) = matrix(i1,i1)
   End Do
   Write(n,*) vec
   Deallocate( vec )
    
  Else
    If (dim_matrix.Le.4) Then
     Do i1=1,dim_matrix
      If (dim_matrix.Eq.2) Write(n,902) matrix(i1,:)
      If (dim_matrix.Eq.3) Write(n,903) matrix(i1,:)
      If (dim_matrix.Eq.4) Write(n,904) matrix(i1,:)
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,901) i1,i2,matrix(i1,i2)
      End Do
     End Do
    End If
  End If

 901 Format(t2,2i4,f16.5)
 902 Format(t2,2f16.5)
 903 Format(t2,3f16.5)
 904 Format(t2,4f16.5)

 End Subroutine WriteRealMatrix


End Module InputOutput

