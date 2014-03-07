#define BEAMSTRAHLUNG
Module EplusEminusProduction
! comments
! This module is a collection of subroutines and functions that are
! needed to calculate of cross sections in e+ e- annihilation. 

! load modules
 Use Control
 Use Mathematics, Only: DGauss, Li2, Vegas1, PolInt, sq_kappa => kappa
 Use StandardModel, Only: mZ, mZ2, gmZ, gmZ2, mW2, mW, Alpha, &
                & mf_l2, mf_l, KFactorLee, alpha_mz, mf_u2, mf_u
 Use Couplings
 Use LoopCouplings, Only: RunningCouplings, InitializeLoopCouplings
! load modules

! interfaces
 Interface CalculateCrossSections
  Module Procedure CalculateCrossSectionsMSSM, CalculateCrossSectionsNMSSM &
                 & , CalculateCrossSectionsRPeps
 End Interface

 Interface EpEmToSfermions
#ifdef GENERATIONMIXING
  Module Procedure EpEmToSfermionsZG, EpEmToSelectronsMSSM                    &
   &             , EpEmToESneutrinosMSSM, EpEmToSquarksZG, EpEmToSleptonsMSSM &
   &             , EpEmToEsneutrinosMSSM3
#else
  Module Procedure EpEmToSfermionsZG, EpEmToSelectronsMSSM                    &
   &             , EpEmToESneutrinosMSSM, EpEmToSquarksZG
#endif
 End Interface

 Interface EpEmToScalarZ
  Module Procedure EpEmScalarZMSSM, EpEmScalarZrp
 End Interface

 Interface EpEmtoHpHm
  Module Procedure EpEmtoHpHmMSSM
 End Interface
! interfaces

! public variables
  Logical :: Only_Photon=.False., Only_Z=.False.
! private variables
! internal indices of produced particles + sizes of actual matrices
  Integer, Private :: ind_1, ind_2, n_S0, n_sf
! to handle polarisation information internally
  Real(dp), Private :: P_m, P_p
! internal vevs
  Real(dp), Private :: vevSM_in(2), vL_in(3)
! internal mixing matrices, most general as needed
  Complex(dp), Private :: Rspm_in(8,8), Rsf_in(6,6), V_in(5,5), U_in(5,5) &
    & , N_in(7,7)
  Real(dp), Private :: RP0_in(5,5), RS0_in(5,5)
! quantum numbers
  Real(dp), Private :: T3_in, e_in

! for sfermion production
  Integer, private :: n_n  ! number of neutralinos in case of sleptons 
  Real(dp), Private :: gg, gZ, ZZ, mSfer2(2), smax, smin,  &
                          & fgg, fgz, fzz
! for QCD corrections
  Integer, Private :: icase
  Real(dp), Private :: msq1, msq2, msq12, msq22, mq2, cos2th, sin2th &
     &      ,mglu2, a(2,2), vq, aq, Qq, mglumq!, mq, QCDfact
  Complex(dp), Private :: costh, sinth
! for the process e+ e- -> selectrons, charged scalars 
  Real(dp), Private :: mk(7), mk2(7), fgN(7), fZNr(7), fZNi(7), fNN(2,7,7)
! for the process e+ e- -> sneutrinos
  Real(dp), Private :: fZV(2), fVV(2,2), ml2(2), mSnu2, mSn2(3)
! for the process e+ e- -> S^0_i P^0_j
  Real(dp), Private :: mS2, mP2
! for the process e+ e- -> chi0_i chi0_j, chi^-_i chi^+_j
  Integer, Private :: n1, n2
  Real(dp), Private :: mi, mj, mi2, mj2, PolFactorZ(2), AngularFactorZ(3)  &
     & , AngularFactorZsel(2), mSel2(2), SumMass2(2), ProdMass(2)          &
     & , cT1, cT2, PolFactorSel(2), PolFactorSel2(5), ZZi(2)
  Complex(dp), Private :: Rse(2,2)
  Complex(dp), Private, Allocatable :: Nn(:,:)
  Logical, private :: l_nmssm = .False.
! for e+e- -> q bar(q)
Real(dp), private :: QCDfact
! for e+e- -> W+ W-
Real(dp), Private :: RmZmW2, RmZmW4, FacWW(3)
#ifdef BEAMSTRAHLUNG
! for beamstrahlung
  Real(dp), Private :: Upsilon, KappaBeam, Ngeff, radprob, EtaBeam, &
                         &  ISRfactor, zmin, zmax, beta, ymin
  Real(dp) :: Spline(201,2), DeltaS, DeltaSplineX(5), DeltaSplineY(5)
  Logical, Private, Save :: CalculateSpline
#else
  Real(dp), Private :: zmin, zmax, beta
#endif
! private Routines
#ifdef BEAMSTRAHLUNG
  Private :: Lee, ILee, BeamElectronDistribution, ISRElectronDistribution 
#else
  Private :: Lee, ILee
#endif
! private variables
 
Contains


#ifdef GENERATIONMIXING
 Subroutine CalculateCrossSectionsMSSM(E, Pm, Pp, ISR, Beam, Design           &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mSlepton, RSlepton, Ylp, mSneut, RSneut, SigSle, SigSn         &
           & , mC, U, V, mN, N, SigC, SigmaN                                  &
           & , mS0, RS0, vevSM, mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp  )
#else
 Subroutine CalculateCrossSectionsMSSM(E, Pm, Pp, ISR                         &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mSlepton, RSlepton, mSneut, RSneut, SigSle, SigSn         &
           & , mC, U, V, mN, N, SigC, SigmaN                                  &
           & , mS0, RS0, vevSM, mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp )
#endif
 !------------------------------------------------------------------------
 ! Calculates all SUSY cross sections for a specific energy
 ! 16.09.02: loopcouplings requires now the information of the MSSM vevs
 ! 12.10.02: starting to implement beamstrahlung
 !-----------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: E, Pm, Pp, mSup(6), mf_u(3), mSdown(6), mf_d(3) &
      & , mglu, mSlepton(6), mSneut(3), mC(2), mN(4), mS0(2), RS0(2,2)    &
      & , mP0(2), RP0(2,2), mSpm(2), vevSM(2)
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6), RSlepton(6,6)        &
      & , RSneut(3,3), U(2,2), V(2,2), N(4,4), RSpm(2,2)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: Ylp(3,3)
#endif
  Logical, Intent(in) :: ISR
!  Logical, Intent(in) :: ISR, Beam
!  Character (Len=*), Intent(in) :: Design

#ifdef GENERATIONMIXING
  Character (Len=8), intent(in) :: Design
  Logical, intent(in) :: Beam
#endif

  Real(dp), Intent(out) :: SigSup(6,6), SigSdown(6,6), SigSle(6,6)  &
      & , SigSn(3,3) , SigC(2,2), SigmaN(4,4), SigS0(2), SigSP(2), SigHp

  Integer :: i1, i2, i3
  Real(dp) :: Emax2, mSf(2), mSne, Rsn2
  Complex(dp) :: Rsf(2,2), id3C(3,3)
  character(len=9) :: specie
  ! for later implementatio of beamstrahlung
#ifdef GENERATIONMIXING
#else
  Character (Len=8) :: Design = "TESLA800"
  Logical :: Beam = .False.
#endif


  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateCrossSectionsMSSM"
  Emax2 = E**2

  SigSup = 0._dp
  SigSdown = 0._dp
  SigSn = 0._dp
  SigC = 0._dp
  SigmaN = 0._dp
  SigS0 = 0._dp
  SigSP = 0._dp
  SigHp = 0._dp

  call InitializeLoopCouplings(vevSM)
  id3C = 0._dp
  id3C(1,1) = 1._dp
  id3C(2,2) = 1._dp
  id3C(3,3) = 1._dp

  !-----------
  ! u-squarks
  !-----------
  specie = 'u-squark'
  SigSup = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSup, Rsup, mf_U((i1+1)/2), mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSup(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSup(2*i1-1:2*i1)
    RSf = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_U(i1), mglu, Pm, Pp &
            & , Emax2 , ISR, Beam , SigSup(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !-----------
  ! d-squarks
  !-----------
  specie = 'd-squark'
  SigSdown = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSdown, Rsdown, mf_D((i1+1)/2),mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSdown(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSdown(2*i1-1:2*i1)
    RSf = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_D(i1), mglu, Pm, Pp &
           & , Emax2 , ISR, Beam , SigSdown(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !-----------
  ! sleptons
  !-----------
  specie = 'slepton'
  SigSle = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSleptonsMSSM(i1, i2, mSlepton, Rslepton, Ylp, id3C, id3C  &
          &, mN, N, Pm, Pp, Emax2, ISR, Beam , SigSle(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSlepton(2*i1-1:2*i1)
    RSf = RSlepton(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      if (i1.eq.1) then
       call EpEmToSfermions(i2, i3, mSf, Rsf, mN, N, Pm, Pp, Emax2, ISR, Beam &
                            &, sigSle(i2,i3), design )
      else
       Call EpEmToSfermions(i2, i3, i1, specie, mSf, Rsf, Pm, Pp, Emax2 &
                   &  , ISR, Beam , SigSle(2*(i1-1)+i2, 2*(i1-1)+i3), design )
      end if
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !------------------
  ! sneutrinos
  !------------------
  SigSn = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,3
    Do i2=1,3
     call EpEmToSfermions(i1, i2, mSneut, Rsneut, id3C, mC, V, Pm, Pp &
                        & , Emax2, ISR, Beam, SigSn(i1,i2), Design)
    end do
   end do
  else
#endif
   call EpEmToSfermions( mSneut(1), mC, V, Pm, Pp, Emax2, ISR, Beam &
                      & , SigSn(1,1), Design)
   Rsf = id2C
   mSf = 1.e6_dp
   specie = 'sneutrino'
   Do i2=2,3
    mSf(1) = mSneut(i2)
    Call EpEmToSfermions(1, 1, i2, specie, mSf, Rsf, Pm, Pp, Emax2 &
                   &  , ISR, Beam , SigSn(i2, i2), design )
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !--------------------------
  ! neutralino production
  !--------------------------
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   i2=1
   Do i1=1,6
    If ( (Abs(Rslepton(i1,1))**2 + Abs(Rslepton(i1,4))**2).gt.0.5_dp) then
     mSf(i2) = mSlepton(i1)
     Rsf(i2,1) = Rslepton(i1,1)
     Rsf(i2,2) = Rslepton(i1,4)
     i2 = i2 + 1
     if (i2.eq.3) exit
    end if
   end do
   if (i2.lt.3) then
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Could not determine Selectrons for Neutralino production."
    If (ErrorLevel.Ge.-1) call TerminateProgram
   end if

  else
#endif
   mSf = mSlepton(1:2)
   Rsf = RSlepton(1:2,1:2)
#ifdef GENERATIONMIXING
  end if
#endif
  SigmaN = 0._dp
  Do i1=1,4
   Do i2=i1,4
    Call EpEmToNeutralinos(i1, i2, mN, N, mSf, RSf, Pm, Pp, Emax2 &
            & , ISR, Beam, -1._dp, 1._dp, SigmaN(i1,i2), design )
    If (i1.Ne.i2) SigmaN(i2,i1) = SigmaN(i1,i2)
   End Do
  End Do

  !---------------------
  ! chargino production
  !---------------------
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Rsn2 = Abs(RSneut(1,1))**2
   Do i1=2,3
    If (Abs(RSneut(i1,1))**2.gt.Rsn2) then
     mSne = mSneut(i1)
     Rsn2 = Abs(RSneut(i1,1))**2
    end if
   end do
  else
#endif
   mSne = mSneut(1)
#ifdef GENERATIONMIXING
  end if
#endif

  SigC = 0._dp
  Do i1=1,2
   Do i2=1,2
    Call EpEmToCharginos(i1, i2, mC, U, V, mSne, Pm, Pp, Emax2, ISR  &
            & ,Beam , SigC(i1,i2), design)
   End Do
  End Do

  !--------------------------
  ! (h0, H0) + A0 production
  !--------------------------
  Do i1=1,2
   call EpEmPseudoScalarScalar(2, i1, mP0, RP0, mS0, RS0, Pm, Pp, Emax2, ISR &
                                   &     ,sigSP(i1) )
  end do

  !--------------------------
  ! (h0, H0) + Z production
  !--------------------------
  Do i1=1,2
   Call EpEmToScalarZ(i1, mS0, RS0, vevSM, Pm, Pp, Emax2, ISR, sigS0(i1) )
  end do

  !---------------------
  ! charged Higgs boson
  !---------------------
  Call EpEmtoHpHm(mSpm, RSpm, Pm, Pp, Emax2, ISR, SigHp)

  Iname = Iname - 1

 End Subroutine CalculateCrossSectionsMSSM


#ifdef GENERATIONMIXING
 Subroutine CalculateCrossSectionsNMSSM(E, Pm, Pp, ISR, Beam, Design           &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mSlepton, RSlepton, Ylp, mSneut, RSneut, SigSle, SigSn         &
           & , mC, U, V, mN, N, SigC, SigmaN                                  &
           & , mS0, RS0, vevSM, mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp  )
#else
 Subroutine CalculateCrossSectionsNMSSM(E, Pm, Pp, ISR                         &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mSlepton, RSlepton, mSneut, RSneut, SigSle, SigSn         &
           & , mC, U, V, mN, N, SigC, SigmaN                                  &
           & , mS0, RS0, vevSM, mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp )
#endif
 !------------------------------------------------------------------------
 ! Calculates all SUSY cross sections for a specific energy
 ! 16.09.02: loopcouplings requires now the information of the NMSSM vevs
 ! 12.10.02: starting to implement beamstrahlung
 !-----------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: E, Pm, Pp, mSup(6), mf_u(3), mSdown(6), mf_d(3) &
      & , mglu, mSlepton(6), mSneut(3), mC(2), mN(5), mS0(3), RS0(3,3)    &
      & , mP0(3), RP0(3,3), mSpm(2), vevSM(2)
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6), RSlepton(6,6)        &
      & , RSneut(3,3), U(2,2), V(2,2), N(5,5), RSpm(2,2)
#ifdef GENERATIONMIXING
  Complex(dp), Intent(in) :: Ylp(3,3)
#endif
  Logical, Intent(in) :: ISR
!  Logical, Intent(in) :: ISR, Beam
!  Character (Len=*), Intent(in) :: Design

#ifdef GENERATIONMIXING
  Character (Len=8), intent(in) :: Design
  Logical, intent(in) :: Beam
#endif

  Real(dp), Intent(out) :: SigSup(6,6), SigSdown(6,6), SigSle(6,6)  &
      & , SigSn(3,3) , SigC(2,2), SigmaN(5,5), SigS0(3), SigSP(3,2), SigHp

  Integer :: i1, i2, i3
  Real(dp) :: Emax2, mSf(2), mSne, Rsn2
  Complex(dp) :: Rsf(2,2), id3C(3,3)
  character(len=9) :: specie
  ! for later implementatio of beamstrahlung
#ifdef GENERATIONMIXING
#else
  Character (Len=8) :: Design = "TESLA800"
  Logical :: Beam = .False.
#endif

  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateCrossSectionsNMSSM"
  Emax2 = E**2

  SigSup = 0._dp
  SigSdown = 0._dp
  SigSn = 0._dp
  SigC = 0._dp
  SigmaN = 0._dp
  SigS0 = 0._dp
  SigSP = 0._dp
  SigHp = 0._dp

  call InitializeLoopCouplings(vevSM)
  id3C = 0._dp
  id3C(1,1) = 1._dp
  id3C(2,2) = 1._dp
  id3C(3,3) = 1._dp

  !-----------
  ! u-squarks
  !-----------
  specie = 'u-squark'
  SigSup = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSup, Rsup, mf_U((i1+1)/2), mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSup(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSup(2*i1-1:2*i1)
    RSf = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_U(i1), mglu, Pm, Pp &
            & , Emax2 , ISR, Beam , SigSup(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !-----------
  ! d-squarks
  !-----------
  specie = 'd-squark'
  SigSdown = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSdown, Rsdown, mf_D((i1+1)/2),mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSdown(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSdown(2*i1-1:2*i1)
    RSf = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_D(i1), mglu, Pm, Pp &
           & , Emax2 , ISR, Beam , SigSdown(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !-----------
  ! sleptons
  !-----------
  specie = 'slepton'
  SigSle = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSleptonsMSSM(i1, i2, mSlepton, Rslepton, Ylp, id3C, id3C  &
          &, mN, N, Pm, Pp, Emax2, ISR, Beam , SigSle(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSlepton(2*i1-1:2*i1)
    RSf = RSlepton(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      if (i1.eq.1) then
       call EpEmToSfermions(i2, i3, mSf, Rsf, mN, N, Pm, Pp, Emax2, ISR, Beam &
                            &, sigSle(i2,i3), design )
      else
       Call EpEmToSfermions(i2, i3, i1, specie, mSf, Rsf, Pm, Pp, Emax2 &
                   &  , ISR, Beam , SigSle(2*(i1-1)+i2, 2*(i1-1)+i3), design )
      end if
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !------------------
  ! sneutrinos
  !------------------
  SigSn = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,3
    Do i2=1,3
     call EpEmToSfermions(i1, i2, mSneut, Rsneut, id3C, mC, V, Pm, Pp &
                        & , Emax2, ISR, Beam, SigSn(i1,i2), Design)
    end do
   end do
  else
#endif
   call EpEmToSfermions( mSneut(1), mC, V, Pm, Pp, Emax2, ISR, Beam &
                      & , SigSn(1,1), Design)
   Rsf = id2C
   mSf = 1.e6_dp
   specie = 'sneutrino'
   Do i2=2,3
    mSf(1) = mSneut(i2)
    Call EpEmToSfermions(1, 1, i2, specie, mSf, Rsf, Pm, Pp, Emax2 &
                   &  , ISR, Beam , SigSn(i2, i2), design )
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !--------------------------
  ! neutralino production
  !--------------------------
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   i2=1
   Do i1=1,6
    If ( (Abs(Rslepton(i1,1))**2 + Abs(Rslepton(i1,4))**2).gt.0.5_dp) then
     mSf(i2) = mSlepton(i1)
     Rsf(i2,1) = Rslepton(i1,1)
     Rsf(i2,2) = Rslepton(i1,4)
     i2 = i2 + 1
     if (i2.eq.3) exit
    end if
   end do
   if (i2.lt.3) then
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Could not determine Selectrons for Neutralino production."
    If (ErrorLevel.Ge.-1) call TerminateProgram
   end if

  else
#endif
   mSf = mSlepton(1:2)
   Rsf = RSlepton(1:2,1:2)
#ifdef GENERATIONMIXING
  end if
#endif
  SigmaN = 0._dp
  l_NMSSM = .True. ! is necessary to use the correct chi0-chi0-Z coupling
  Do i1=1,5
   Do i2=i1,5
    Call EpEmToNeutralinos(i1, i2, mN, N, mSf, RSf, Pm, Pp, Emax2 &
            & , ISR, Beam, -1._dp, 1._dp, SigmaN(i1,i2), design )
    If (i1.Ne.i2) SigmaN(i2,i1) = SigmaN(i1,i2)
   End Do
  End Do

  !---------------------
  ! chargino production
  !---------------------
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Rsn2 = Abs(RSneut(1,1))**2
   Do i1=2,3
    If (Abs(RSneut(i1,1))**2.gt.Rsn2) then
     mSne = mSneut(i1)
     Rsn2 = Abs(RSneut(i1,1))**2
    end if
   end do
  else
#endif
   mSne = mSneut(1)
#ifdef GENERATIONMIXING
  end if
#endif

  SigC = 0._dp
  Do i1=1,2
   Do i2=1,2
    Call EpEmToCharginos(i1, i2, mC, U, V, mSne, Pm, Pp, Emax2, ISR  &
            & ,Beam , SigC(i1,i2), design)
   End Do
  End Do

  !--------------------------
  ! (h0, H0) + A0 production
  !--------------------------
  Do i1=1,3
   Do i2=1,2
   Call EpEmPseudoScalarScalar(i2+1, i1, mP0, RP0, mS0, RS0, Pm, Pp, Emax2, ISR &
                             &     ,sigSP(i1,i2) )
   end do
  end do
  !--------------------------
  ! (h0, H0) + Z production
  !--------------------------
  Do i1=1,3
   Call EpEmToScalarZ(i1, mS0, RS0, vevSM, Pm, Pp, Emax2, ISR, sigS0(i1) )
  end do
  !---------------------
  ! charged Higgs boson
  !---------------------
  Call EpEmtoHpHm(mSpm, RSpm, Pm, Pp, Emax2, ISR, SigHp)

  Iname = Iname - 1

 End Subroutine CalculateCrossSectionsNMSSM

#ifdef GENERATIONMIXING
 Subroutine CalculateCrossSectionsRPeps(E, Pm, Pp, ISR, Beam, Design          &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mC, U, V, mN, N, SigC, SigmaN, mS0, RS0, vevSM, vevL           &
           & , mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp  )
#else
 Subroutine CalculateCrossSectionsRPeps(E, Pm, Pp, ISR                        &
           & , mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu, SigSup, SigSdown &
           & , mC, U, V, mN, N, SigC, SigmaN, mS0, RS0, vevSM, vevL           &
           & , mP0, RP0, mSpm, RSpm, SigS0, SigSP, SigHp )
#endif
 !------------------------------------------------------------------------
 ! Calculates all SUSY cross sections for a specific energy
 ! 16.09.02: loopcouplings requires now the information of the MSSM vevs
 ! 12.10.02: starting to implement beamstrahlung
 !-----------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) :: E, Pm, Pp, mSup(6), mf_u(3), mSdown(6), mf_d(3) &
      & , mglu, mC(5), mN(7), mS0(5), RS0(5,5), mP0(5), RP0(5,5), mSpm(8) &
      & , vevSM(2), vevL(3)
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6), U(5,5), V(5,5)       &
      & , N(7,7), RSpm(8,8)
!#ifdef GENERATIONMIXING
!  Complex(dp), Intent(in) :: Ylp(3,3)
!#endif
  Logical, Intent(in) :: ISR
!  Logical, Intent(in) :: ISR, Beam
!  Character (Len=*), Intent(in) :: Design

#ifdef GENERATIONMIXING
  Character (Len=8), intent(in) :: Design
  Logical, intent(in) :: Beam
#endif

  Real(dp), Intent(out) :: SigSup(6,6), SigSdown(6,6), SigC(5,5), SigmaN(7,7) &
   & , SigS0(5), SigSP(5,4), SigHp(7,7)

  Integer :: i1, i2, i3
  Real(dp) :: Emax2, mSf(2), mSne
  Complex(dp) :: Rsf(2,2), id3C(3,3)
  character(len=9) :: specie
  ! for later implementatio of beamstrahlung
#ifdef GENERATIONMIXING
#else
  Character (Len=8) :: Design = "TESLA800"
  Logical :: Beam = .False.
#endif


  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateCrossSectionsRPeps"
  Emax2 = E**2

  SigSup = 0._dp
  SigSdown = 0._dp
  SigC = 0._dp
  SigmaN = 0._dp
  SigS0 = 0._dp
  SigSP = 0._dp
  SigHp = 0._dp

  call InitializeLoopCouplings(vevSM)
  id3C = 0._dp
  id3C(1,1) = 1._dp
  id3C(2,2) = 1._dp
  id3C(3,3) = 1._dp

  !-----------
  ! u-squarks
  !-----------
  specie = 'u-squark'
  SigSup = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSup, Rsup, mf_U((i1+1)/2), mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSup(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSup(2*i1-1:2*i1)
    RSf = RSup(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_U(i1), mglu, Pm, Pp &
            & , Emax2 , ISR, Beam , SigSup(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !-----------
  ! d-squarks
  !-----------
  specie = 'd-squark'
  SigSdown = 0._dp
#ifdef GENERATIONMIXING
  if (GenerationMixing) then
   Do i1=1,6
    Do i2=1,6
     Call EpEmToSquarksZG(i1, i2, specie, mSdown, Rsdown, mf_D((i1+1)/2),mglu &
                   & , Pm, Pp, Emax2 , ISR, Beam , SigSdown(i1,i2), design )
    end do
   end do
  else
#endif
   Do i1=1,3
    mSf = mSdown(2*i1-1:2*i1)
    RSf = RSdown(2*i1-1:2*i1, 2*i1-1:2*i1)
    Do i2=1,2
     Do i3=1,2
      Call EpEmToSquarksZG(i2, i3, specie, mSf, Rsf, mf_D(i1), mglu, Pm, Pp &
           & , Emax2 , ISR, Beam , SigSdown(2*(i1-1)+i2, 2*(i1-1)+i3), design )
     end do
    end do
   end do
#ifdef GENERATIONMIXING
  end if
#endif

  !--------------------------
  ! neutralino production
  !--------------------------
  Rsf = id2C
  Do i1=2,8
   If (Abs(RSpm(i1,3)).gt.0.5) mSf(1) = mSpm(i1)
   If (Abs(RSpm(i1,6)).gt.0.5) mSf(2) = mSpm(i1)
  end do
  
  SigmaN = 0._dp
  Do i1=1,7
   Do i2=1,7
    Call EpEmToNeutralinos(i1, i2, mN, N, mSf, RSf, Pm, Pp, Emax2 &
            & , ISR, Beam, -1._dp, 1._dp, SigmaN(i1,i2), design )
   End Do
  End Do

  !---------------------
  ! chargino production
  !---------------------
  Do i1=1,5
   If (Abs(RS0(i1,3)).gt.0.5) mSne = mS0(i1)
  end do
  SigC = 0._dp
  Do i1=1,5
   Do i2=1,5
    Call EpEmToCharginos(i1, i2, mC, U, V, mSne, Pm, Pp, Emax2, ISR  &
            & ,Beam , SigC(i1,i2), design)
   End Do
  End Do

  !--------------------------
  ! (h0, H0) + A0 production
  !--------------------------
  Do i1=1,5
   Do i2=1,4
    Call EpEmPseudoScalarScalar(i2+1, i1, mP0, RP0, mS0, RS0, Pm, Pp &
                              &  , Emax2, ISR, sigSP(i1,i2) )
   end do
  end do

  !--------------------------
  ! (h0, H0) + Z production
  !--------------------------
  Do i1=1,5
   Call EpEmScalarZRP(i1, mS0, RS0, vevSM, vevL, Pm, Pp, Emax2, ISR, sigS0(i1))
  end do

  !-------------------------------------
  ! charged scalars, currently missing
  !-------------------------------------
  SigHp = 0._dp
  Do i1=2,8
   Do i2=2,8
    Call EpEmToChargedScalarsRPeps(i1, i2, mSpm, Rspm, mN, N, U, V, Pm, Pp &
                               & , Emax2, ISR, Beam, SigHp(i1-1,i2-1) )
   end do
  end do

  Iname = Iname - 1

 End Subroutine CalculateCrossSectionsRPeps


 Subroutine InitializeCrossSections(Ecms, Pm, Pp, ISR)
 !------------------------------------------------------------------------------
 ! reads in the data from the file CrossSections.in
 ! for the initialisation of the e+ e- annihilation processes
 ! output:
 !  Ecms(:) ...... vector containing the c.m.s. energies
 !  Pm(:) ........ vector containing the degree of electron polarisation
 !  Pp(:) ........ vector containing the degree of positron polarisation
 !  ISR(:) ....... vector containing the information if ISR shall be calculated
 !------------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(inout) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(inout) :: ISR(:)

  Integer :: i1
  Ecms = 0._dp
  Pm = 0._dp
  Pp = 0._dp
  ISR = .False.

  Open(92, file="CrossSections.in", status="old", err=310)
  Do i1=1,Size(Ecms)
   Read(92, *, End=300) Ecms(i1)
   Read(92,*, End=300) Pm(i1)
   Read(92,*, End=300) Pp(i1)
   Read(92,*, End=300) ISR(i1)
  End Do
  300 Close(92)
  310 Return

 End Subroutine InitializeCrossSections


 Real(dp) Function EpEmToChargedScalarsTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of selectron production in
 ! e+e- annihilation. Is called by eeChargedScalars and related functions.
 ! written by Werner Porod, 14.10.2005
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Integer :: i1,i2
  Real(dp) :: kappa,kappa3d2,sqrts,gin(3),yukin(3), &
        &   propN(7),kappa1d2,sumGauge2,sumGaugeN,sumN2,invPropZ,   &
        &   propNN(2,7,7),logNN(7),sumMkMsf(7)

  Real(dp) :: Le, Re, ProdPmPp(2), sinW2, cosW2
  Complex(dp) :: coup, yukC(3)
  Complex(dp) :: ai(7), bi(7), aiC(7), biC(7), aj(7), bj(7), ajC(7) &
      &  , bjC(7), coupC

  If (mSfer2(1).Eq.mSfer2(2)) Then
   kappa = s*(s-4._dp*mSfer2(1))
   kappa1d2 = Sqrt(kappa)
   kappa3d2 = kappa**1.5_dp
   Do i1=1,7
    logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-2._dp*mSfer2(1) ) /   &
         &           (s+2._dp*mk2(i1)-kappa1d2-2._dp*mSfer2(1) ) )
    sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))**2
    propN(i1) = (mSfer2(1)-0.5_dp*s-mk2(i1) )*kappa1d2              &
         &    + sumMkMsf(i1) * logNN(i1)
   End Do
   Do i1=1,7
    Do i2=i1,7
! the second condition is due to the fact that neutrino masses
! are signifcantly smaller than SUSY masses implying that one runs
! into numerical problems with the second formula
     If ((mk2(i1).Eq.mk2(i2)).or.(i2.le.3)) Then
      propNN(1,i1,i2) = -2._dp * kappa1d2    &
       &             + (s+2._dp*mk2(i1)-2._dp*mSfer2(1) ) * logNN(i1)
      propNN(2,i1,i2) = s * kappa1d2 / sumMkMsf(i1)
     Else
      propNN(1,i1,i2) = - kappa1d2                                          &
          &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
          &             / (mk2(i1)-mk2(i2))
      propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
     End If
    End Do
   End Do

  Else
   kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
   kappa1d2 = Sqrt(kappa)
   kappa3d2 = kappa**1.5_dp
   Do i1=1,7
    logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-mSfer2(1)-mSfer2(2) ) /      &
         &           (s+2._dp*mk2(i1)-kappa1d2-mSfer2(1)-mSfer2(2) ) )
    sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))*(mk2(i1)-mSfer2(2))
    propN(i1) = (0.5_dp*(mSfer2(1)+mSfer2(2)-s)-mk2(i1) )*kappa1d2          &
         &    + sumMkMsf(i1) * logNN(i1)
   End Do
   Do i1=1,7
    Do i2=i1,7
! the second condition is due to the fact that neutrino masses
! are signifcantly smaller than SUSY masses implying that one runs
! into numerical problems with the second formula
     If ((mk2(i1).Eq.mk2(i2)).or.(i2.le.3)) Then
      propNN(1,i1,i2) = -2._dp * kappa1d2    &
       &             + (s+2._dp*mk2(i1)-mSfer2(1)-mSfer2(2) ) * logNN(i1)
      propNN(2,i1,i2) = s * kappa1d2 /sumMkMsf(i1)
     Else
      propNN(1,i1,i2) = - kappa1d2                                          &
          &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
          &             / (mk2(i1)-mk2(i2))
      propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
     End If
    End Do
   End Do
  End If 
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  yukC = 0._dp
  yukC(3) = yukin(3)

  cosW2 = gin(2)**2/ (gin(1)**2 + gin(2)**2)
  sinW2 = 1._dp - cosW2
  Call CoupChargedScalarZ(ind_1, ind_2, gin(2), sinW2, Rspm_in, coup)
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ProdPmPp(1) = (1-P_m)*(1+P_p)
  ProdPmPp(2) = (1+P_m)*(1-P_p)

  If (ind_1.Eq.ind_2) Then
   fgg = oo16pi * gin(2)**4 * sinW2**2 * (1._dp - P_m*P_p) / 3._dp
   fgZ = -oo16pi * gin(2)**2 * sinW2 * Conjg(coup)   &
       &         * (Le*ProdPmPp(1)+Re*ProdPmPp(2)) /3._dp

!---------------------------------------
! neglecting electron Yukawa coupling
!---------------------------------------
   Do i1=1,7
    Call CoupCSCharginoNeutralino(ind_1, 1, i1, N_in, U_in, V_in, RSpm_in, YukC &
           & , gin(1), gin(2), bi(i1), ai(i1))
    fgN(i1) =  Abs( ai(i1) )**2 * ProdPmPp(1) + Abs( bi(i1) )**2 * ProdPmPp(2)
    fZNr(i1) = Le * Abs( ai(i1) )**2 * ProdPmPp(1)  &
             + Re * Abs( bi(i1) )**2 * ProdPmPp(2)
    fZNi(i1) = 0._dp
   End Do
   fgN = oo16pi * gin(2)**2 * sinW2 * fgN 
   fZNr = - oo16pi * coup * fZNr
   Do i1=1,7
    Do i2=i1,7
     fNN(1,i1,i2) = ( Abs( ai(i1) )**2 *Abs( ai(i2) )**2 * ProdPmPp(1) &
         + Abs( bi(i1) )**2 *Abs( bi(i2) )**2 * ProdPmPp(2) )
     coupC = ai(i1) * bi(i2) * Conjg( ai(i2) * bi(i1) )
     fNN(2,i1,i2) = 2._dp * mk(i1) * mk(i2)                                  &
       * ( Real(coupC,dp) * (1 + P_m*P_p) + Aimag(coupC) * (P_m+P_p) )
    End Do
   End Do

  Else
   fgg = 0._dp
   fgZ = 0._dp
!---------------------------------------
! neglecting electron Yukawa coupling
!---------------------------------------
   Do i1=1,7
    Call CoupCSCharginoNeutralino(ind_1, 1, i1, N_in, U_in, V_in, RSpm_in, YukC &
                             & , gin(1), gin(2), bi(i1), ai(i1))
    Call CoupCSCharginoNeutralino(ind_2, 1, i1, N_in, U_in, V_in, RSpm_in, YukC &
                             & , gin(1), gin(2), bj(i1), aj(i1))
    fgN(i1) = 0._dp
    aiC(i1) = Conjg( ai(i1) )
    biC(i1) = Conjg( bi(i1) )
    ajC(i1) = Conjg( aj(i1) )
    bjC(i1) = Conjg( bj(i1) )
    coupC = Le * aiC(i1) * aj(i1) * ProdPmPp(1) &
          + Re * biC(i1) * bj(i1) * ProdPmPp(2)
    coupC = coup * coupC
    fZNr(i1) = Real( coupC,dp )
    fZNi(i1) = Aimag( coupC ) 
   End Do
   fZNr = - oo16pi * fZNr
   fZNi = - oo16pi * gmZ * fZNi

   Do i1=1,7
    Do i2=i1,7
     fNN(1,i1,i2) = Real( ai(i1) *aiC(i2) *ajC(i1) *aj(i2) * ProdPmPp(1)    &
                  &     + bi(i1) *biC(i2) *bjC(i1) *bj(i2) * ProdPmPp(2),dp )
     fNN(2,i1,i2) = mk(i1) * mk(i2)                                            &
                  & * ( aj(i2) *ajC(i1) *bi(i1) *biC(i2) *(1 - P_m) *(1 - P_p) &
                  &   + ai(i1) *aiC(i2) *bj(i2) *bjC(i1) *(1 + P_m) *(1 + P_p) )
    End Do
   End Do
  End If

  Do i1=1,7
   Do i2=i1,7
    If (i1.Eq.i2) Then
     fNN(1,i1,i2) = oo64pi * fNN(1,i1,i2) 
     fNN(2,i1,i2) = oo64pi * fNN(2,i1,i2) 
    Else
     fNN(1,i1,i2) = oo32pi * fNN(1,i1,i2) 
     fNN(2,i1,i2) = oo32pi * fNN(2,i1,i2) 
    End If
   End Do
  End Do
  fzz = oo16pi * Abs(coup)**2 * (Le**2 *ProdPmPp(1)+Re**2 *ProdPmPp(2)) / 6._dp

  invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )
  sumGauge2 = ( fgg + ( fgZ * (s-mZ2) + fZZ * s ) * s * invPropZ ) * kappa3d2
  sumGauge2 = sumGauge2 / s**4 

  sumGaugeN = 0._dp
  sumN2 = 0._dp
  Do i1=1,7
   sumGaugeN = sumGaugeN &
         &   + (fgN(i1)+ s*((s-mZ2)*fZNr(i1)+fZNi(i1))*invPropZ)*propN(i1)
   Do i2=i1,7
    sumN2 = sumN2 + ( fNN(1,i1,i2)*propNN(1,i1,i2)   &
      &             + fNN(2,i1,i2)*propNN(2,i1,i2) )
   End Do
  End Do
  sumGaugeN =  sumGaugeN / s**3
  sumN2 = sumN2 / s**2

  EpEmToChargedScalarsTree = (sumGauge2+sumGaugeN+sumN2)

 End Function EpEmToChargedScalarsTree


 Subroutine EpEmToChargedScalarsRPeps(i, j, mSpm, Rspm, mN, N, U, V, Pm, Pp, s &
                               & , ISR, Beam, sigma, Design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of charged scalars in the epsilon model
 ! input:
 !  i,j .............. the charged scalar combination
 !  mSpm(i) .......... charged scalar masses
 !  Rspm(i,j) ........ Mixing matrices of charged scalars
 !  mN(i) ............ Neutralino masses
 !  N(i,j) ........... Neutralino mixing matrix
 !  U(i,j), V(i,j) ... chargino mixing matrices
 !  Pm ............... degree of e- polarisation
 !  Pp ............... degree of e+ polarisation
 !  s ................ c.m.s. energy squared
 !  ISR .............. ISR corrections are included if .TRUE.
 !  Beam ............. beam strahlung is included if .TRUE.
 !  Design ........... collider type, necessary in case of beam strahlung
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod, 14.10.2005
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: mSpm(8), Pm, Pp, s, mN(7)
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rspm(8,8), N(7,7), U(5,5), V(5,5)
  Logical, Intent(in) :: ISR, Beam
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToChargedScalarsRPeps'
  sigma = 0._dp

  If ( (i.Lt.1).Or.(i.Gt.8).Or.(j.Lt.1).Or.(j.Gt.8) ) Then
   Write(ErrCan,*) 'Error: in subroutine'//NameOfUnit(Iname)//'the combination'
   Write(ErrCan,*) i,j,' should be calculated. IMPOSSIBLE !!!!'
   Call TerminateProgram
  End If

  !-------------
  ! kinematics
  !-------------
  If ( (mSpm(i)+mSpm(j))**2.Ge.s) Then
   iname = iname - 1
   Return 
  End If

  !--------------------------------------------------------------------
  ! internal variables
  !--------------------------------------------------------------------
  ind_1 = i
  ind_2 = j
  P_m = Pm
  P_p = Pp
  RSpm_in =  RSpm
  N_in = N
  U_in = U
  V_in = V
  mK = mN
  mk2 = mN**2    
  mSfer2(1) = mSpm(i)**2
  mSfer2(2) = mSpm(j)**2

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   ! due to polarization it is possible that all couplings vanish
   ! check therefore first, if tree level is non-zero
   sigma = EpEmToSelectronsTree(s)
   If (sigma.Eq.0._dp) Then
    Iname = Iname - 1
    Return
   End If
   ! now do the real calculation

   smin = (mSpm(i)+mspm(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (mSpm(i)+mSpm(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToSelectronsBeamA,init,ncall,itmx,nprn,1.e-3_dp,erg &
              &, sd, chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = (mSpm(i)+mSpm(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)//' '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToSelectronsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToSelectronsTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToSelectronsISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToChargedScalarsTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * sigma  

  Iname = Iname - 1

 End Subroutine EpEmToChargedScalarsRPeps

 
 Subroutine EpEmToCharginos(i, j, mC, U, V, mSn, Pm, Pp, s, ISR  &
               & ,Beam , sigma, design)
 !--------------------------------------------------------------------------
 ! calculates the cross section for chargino production in e+e- annihilation
 ! the formulas are checked against A.Bartl et al, Z.f.Ph.C 30, 441 (1986)
 ! written by Werner Porod, 17.07.01
 !--------------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j

  Real(dp), Intent(in) :: mC(:), mSn, Pm, Pp, s
  Real(dp), Intent(out) :: sigma

  Complex(dp), Intent(in) :: U(:,:), V(:,:)

  Logical, Intent(in) :: ISR, Beam

  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn, n_c
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToCharginos'

  !-------------------------------
  ! initialisation
  !-------------------------------
  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (Abs(mC(i))+Abs(mC(j)))**2.Ge.s) Then
   Iname = Iname - 1
   Return
  End If

  mi = mC(i)
  mj = mC(j)
  mi2 = mC(i)**2
  mj2 = mC(j)**2
  mSnu2 = mSn**2

  !--------------------------------------------------------------------
  ! internal variables
  !--------------------------------------------------------------------
  ind_1 = i
  ind_2 = j

  P_m = Pm
  P_p = Pp

  n_c = Size(mC)
  V_in = 0
  U_in = 0
  V_in(1:n_c,1:n_c) = V
  U_in(1:n_c,1:n_c) = U
 

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   smin = (mC(i)+mC(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (mC(i)+mC(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToCharginosBeam,init,ncall,itmx,nprn,1.e-3_dp,erg &
             ,sd, chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = (mC(i)+mC(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToCharginosTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToCharginosTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToCharginosISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToCharginosTree(s)
  End If

  !----------------------------
  ! 0.38939e12 gives fb
  !----------------------------
  sigma = 0.38939e12_dp * sigma

  Iname = Iname - 1

 End Subroutine EpEmToCharginos


 Real(dp) Function EpEmToCharginosTree(s)
 Implicit None

  Real(dp), Intent(in) :: s

  Real(dp) :: SqrtS, EiEj, p_char, sigma_ij(6), invPropZ, sp_char, s2           &
      &     , Al, Bl, LogAlBl, hab, gin(3), yukin(3), Al2, Bl2, p2_char, AldBl  &
    & , sinW2, cosW, Le, Re, g2, g4, gg, gZ, gSn, ZZi(2), Zsn(2), SnSn
  Complex(dp) :: OLij, ORij

  If (mi2.Eq.mj2) Then
   p2_char = 0.25_dp * s -  mi2
  Else
   p2_char = 0.25_dp * ( (s-mi2-mj2)**2 - 4._dp * mi2 * mj2 ) / s
  End If
  SqrtS = Sqrt( s )
  sp_char = Sqrt( p2_char ) * SqrtS
  p_char = sp_char / s
  
  sigma_ij = 0._dp

  invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )

  !------------------------
  ! running couplings
  !------------------------
  Call runningCouplings(Sqrts,gin,yukin)
  g2 = gin(2)**2
  g4 = g2**2

  sinW2 = gin(1)**2 / (gin(1)**2 + g2)
  cosW = Sqrt(1._dp - sinW2)

  Call CoupCharginoZ(ind_1, ind_2, U_in, V_in, gin(2), cosW, OLij, ORij)
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)

  If (ind_1.Eq.ind_2) Then
   gg = g4 * (1._dp - P_m * P_p) * sinW2**2 / (2._dp * Pi )
   gZ = g2 * sinW2 * oo4pi * ( Le * (1._dp - P_m) * (1._dp + P_p)            &
      &                      + Re * (1._dp + P_m) * (1._dp - P_p) )          &
      &    * Real(OLij + ORij,dp)
   gSn = oo16pi * g4 * sinW2 * Abs( V_in(ind_1,1) )**2 &
       &        * (1._dp - P_m) * (1._dp + P_p) 
  Else
   gg = 0._dp
   gZ = 0._dp
   gSn = 0._dp
  End If

  ZZi(1) = oo4pi * ( Le**2 * (1._dp - P_m) * (1._dp + P_p)       &
         &         + Re**2 * (1._dp + P_m) * (1._dp - P_p) )
  ZZi(2) = ZZi(1) * Real( OLij * Conjg(ORij),dp ) * mi * mj
  ZZi(1) = 0.5_dp * ZZi(1) * (Abs(OLij)**2 + Abs(ORij)**2 )
  ZSn(1) = oo16pi * g2 * Le * (1._dp - P_m) * (1._dp + P_p)
  ZSn(2) = ZSn(1) * mi * mj * Real(OLij * V_in(ind_1,1)           &
         &                         * Conjg( V_in(ind_2,1) ),dp )
  ZSn(1) = ZSn(1) * Real(ORij * V_in(ind_1,1) * Conjg( V_in(ind_2,1) ),dp )
  SnSn = g4 * oo32pi * Abs( V_in(ind_1,1) )**2 * Abs( V_in(ind_2,1) )**2   &
       &    * (1._dp - P_m) * (1._dp + P_p) / mSnu2**2

  If (mi2.Eq.mj2) Then
   s2 = s**2
   EiEj = p2_char + mi2
   Al = (0.5_dp * s + mSnu2 - mi2) / mSnu2
   Bl = sp_char / mSnu2
   AldBl = (0.5_dp * s + mSnu2 - mi2) / sp_char
   LogAlBl = Log(Abs(Al+Bl) / Abs(Al-Bl) )
   hab = 2._dp * sp_char  - 2._dp * p2_char * AldBl              &
     & + (EiEj + p2_char * AldBl**2 - sp_char * AldBl ) * LogAlBl

   sigma_ij(1) = gg * p_char * (EiEj + p2_char / 3._dp + mi2) / s2
   sigma_ij(2) = gZ * (s-mZ2) * invPropZ * p_char &
             & * (EiEj + p2_char / 3._dp + mi2) / s
   sigma_ij(4) = - gSn * ( hab + mi2 * LogAlBl) / s2
  Else
   EiEj = Sqrt((p2_char + mi2)*(p2_char + mj2))
   Al = (s + 2._dp * mSnu2 - mi2 - mj2) / (2._dp * mSnu2)
   Bl = sp_char / mSnu2
   AldBl = 0.5_dp * (s + 2._dp * mSnu2 - mi2 - mj2) / sp_char
   LogAlBl = Log(Abs(Al+Bl) / Abs(Al-Bl) )
   hab = 2._dp * sp_char - 2._dp * p2_char * AldBl     &
     & + (EiEj + p2_char * AldBl**2 - sp_char * AldBl ) * LogAlBl
  End If
  Al2 = Al**2
  Bl2 = Bl**2

  sigma_ij(3) = invPropZ * p_char  &
           &  * ( ZZi(1) * (EiEj + p2_char / 3._dp) + ZZi(2) ) 
  sigma_ij(5) = - (s-mZ2) * invPropZ * ( ZSn(1) * hab + ZSn(2) * LogAlBl) / s
  sigma_ij(6) = SnSn * p_char                                         &
            & * ( (EiEj + p2_char  - sp_char * AldBl) / (Al2 - Bl2)   &
            &   + 2._dp * p2_char / Bl2                               &
            &   + (sp_char - 2._dp * p2_char * AldBl) * LogAlBl       &
            &     / (2._dp*Bl2)   )

  EpEmToCharginosTree = Sum( sigma_ij )

 End Function EpEmToCharginosTree


  Real(dp) Function EpEmToCharginosISR(x)

  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToCharginosISR'

  s = smax * (1._dp - Exp(x))

  EpEmToCharginosISR = Lee(beta,x) * EpEmToCharginosTree(s)

  Iname = Iname - 1

 End Function EpEmToCharginosISR


#ifdef BEAMSTRAHLUNG
 Real(dp) Function EpEmToCharginosBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeCharginos and related functions.
 ! for the calculation of ISR and Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 

  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToCharginosBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToCharginosTree(s) ! EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToCharginosBeam = BeamFactor * ISRfactor2 * sigmaT

 End Function EpEmToCharginosBeam
#endif


 Subroutine EpEmToNeutralinos(i, j, mN, N, mSel, RSel, Pm, Pp, s, ISR  &
               & ,Beam , cosT1, cosT2, sigma, design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of neutralinos, the case of R-parity 
 ! violaton is included under the assumption that selectrons hardly mix with
 ! the other sleptons and the charged Higgs
 ! input:
 !  i,j ........ indices of neutralinos
 !  mN(i) ...... neutralinos masses 
 !  N(i,j) ..... nixing matrix of neutralinos
 !  mSel(i) .... selectron masses
 !  RSel(i,j) .. mixing matrix for selectrons
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 !  Beam ....... Beam strahlung is included if .TRUE.
 !  cosT1 ...... cos(t_1)
 !  cosT2 ...... cos(t_2) where t_1 and t_2 is the angular range for the
 !               integration
 ! output:
 !  sigma ...... production cross section in fb
 ! written by Werner Porod, 02.07.2001
 !  02.07.2001: taking f77 routines and traansform them to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j

  Real(dp), Intent(in) :: mN(:), mSel(2), Pm, Pp, s, cosT1, cosT2
  Real(dp), Intent(out) :: sigma

  Complex(dp), Intent(in) :: N(:,:), RSel(:,:)

  Logical, Intent(in) :: ISR, Beam

  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: i1, n_neut
  Real(dp) :: z1, z2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToNeutralinos'
  !-----------------------
  ! checking the model
  !-----------------------
  n_neut = Size(mN)
  If ((i.Lt.1).Or.(i.Gt.n_neut)) Then
   Write(ErrCan,*) 'Neutralino index i out of range in subroutine ', &
     &              NameOfUnit(Iname),'i= ',i
   Call TerminateProgram
  End If
  If ((j.Lt.1).Or.(j.Gt.n_neut)) Then
   Write(ErrCan,*) 'Neutralino index j out of range in subroutine ', &
     &              NameOfUnit(Iname),'j= ',j
   Call TerminateProgram
  End If

  !-------------------------------
  ! initialisation
  !-------------------------------
  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (Abs(mN(i))+Abs(mN(j)))**2.Ge.s) Then
   Iname = Iname - 1
   Return
  End If
  !-------------------------------
  ! constant factors for integrals
  !-------------------------------
  n1 = i
  n2 = j
  mi = Abs(mN(i))
  mj = Abs(mN(j))
  mi2 = mN(i)**2
  mj2 = mN(j)**2
  cT1 = cosT1
  cT2 = cosT2
  Allocate( Nn(n_neut,n_neut) )
  Nn = N

  PolFactorZ(1) = (1._dp - Pm) * (1._dp + Pp)
  PolFactorZ(2) = (1._dp + Pm) * (1._dp - Pp)
  AngularFactorZ(3) = (CosT2 - CosT1)
  AngularFactorZ(1) = AngularFactorZ(3) / 24._dp
  AngularFactorZ(2) = 2._dp * ( cosT2**2 + cosT1**2 + cosT2*cosT1 - 1._dp)
  AngularFactorZ(3) = AngularFactorZ(3) * mN(i) * mN(j)

  Do i1=1,2
   mSel2(i1) = mSel(i1)**2
   SumMass2(i1) = mi2 + mj2 - 2._dp *  mSel2(i1)
   ProdMass(i1) = (mi2 - mSel2(i1)) * (mj2 - mSel2(i1) )
  End Do
  AngularFactorZsel(1) = CosT1 - CosT2
  AngularFactorZsel(2) = (CosT1**2 - CosT2**2)
  PolFactorSel(1) = (1._dp - Pm) * (1._dp + Pp) * mN(i) * mN(j)
  PolFactorSel(2) = (1._dp + Pm) * (1._dp - Pp) * mN(i) * mN(j)
  PolFactorSel2(1) = (1._dp - Pm)
  PolFactorSel2(2) = (1._dp + Pm)
  PolFactorSel2(3) = (1._dp - Pp)
  PolFactorSel2(4) = (1._dp + Pp)
  PolFactorSel2(5) = (1._dp + Pm*Pp)
  Rse = Rsel
  !----------------------------------
  ! calculation of the cross section
  !----------------------------------
  If (ISR) Then
   smax = s
   smin = (mi+mj)**2
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   Zmax = Log(1._dp - smin/s)
   If (Zmax.Lt.-8._dp) Then
    sigma = eeNeutralinoTree(s) * ILee(beta,Zmax)
   Else If (smin.Lt.(0.75_dp*mZ2) ) Then
    z1 = Log(1._dp - 0.75_dp*mZ2/s)
    z2 = Log(1._dp - 1.25_dp*mZ2/s)
    Zmin = -10._dp
    sigma = eeNeutralinoTree(s) * ILee(beta,Zmin)        &
     &    + dgauss(eeNeutralinoISR,zmin,z1,1.e-4_dp)     &
     &    + dgauss(eeNeutralinoISR,z1,z2,1.e-4_dp)       &
     &    + dgauss(eeNeutralinoISR,z2,zmax,1.e-4_dp)
   Else If (smin.Lt.(1.25_dp*mZ2) ) Then
    z2 = Log(1._dp - 1.25_dp*mZ2/s)
    Zmin = -10._dp
    sigma = eeNeutralinoTree(s) * ILee(beta,Zmin)        &
     &    + dgauss(eeNeutralinoISR,zmin,z2,1.e-4_dp)     &
     &    + dgauss(eeNeutralinoISR,z2,zmax,1.e-4_dp)
   Else
    Zmin = -10._dp
    sigma = eeNeutralinoTree(s) * ILee(beta,Zmin)        &
     &    + dgauss(eeNeutralinoISR,zmin,zmax,1.e-4_dp)
   Endif
  Else
   sigma = eeNeutralinoTree(s)
  Endif
  !----------------------------
  ! 0.38939e12 gives fb
  !----------------------------
  sigma = oo16pi * 0.38939e12_dp * sigma
  !------------------------------
  ! symmetry factor
  !------------------------------
  If (i.Eq.j) sigma = 0.5_dp * sigma

  Deallocate( Nn )

  Iname = Iname - 1

 End Subroutine EpEmToNeutralinos


  Real(dp) Function eeNeutralinoTree(s)

  Implicit None

  Real(dp), Intent(in) :: s

  Integer i1
  Real(dp) :: Q,g_Q(3),g,gp,sW2,Le,Re,cW,s2,lambda,partZZ,partZsel(2)  &
     &      ,yuk_Q(3),lambda2,ReZprop,LogPcT(2),LogMcT(2),partSel2(2)
  Complex(dp) :: Ol_ij,Or_ij,yuk=0,fLi(2),fRi(2),fLj(2),fRj(2)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'eeNeutralinoTree'
  !-----------------
  ! couplings
  !-----------------
  Q = Sqrt(s)
  Call runningCouplings(Q,g_Q,yuk_Q)
  gp = g_Q(1)
  g = g_Q(2)
  sW2 = gp**2 / (gp**2 + g**2 )

  Call CoupFermionZ(-0.5_dp,-1._dp,g,sW2,Le,Re)
  cW = Sqrt(1._dp - sW2)
  if (l_NMSSM) then
   Call CoupNeutralinoZ(n2,n1,Nn,g,cW,Ol_ij,Or_ij, l_NMSSM)
  else
   Call CoupNeutralinoZ(n2,n1,Nn,g,cW,Ol_ij,Or_ij)
  end if

  yuk = 0._dp 
  Call CoupNeutralinoSlepton(n1,1,gp,g,RSe,Yuk,Nn,fLi(1),fRi(1))
  Call CoupNeutralinoSlepton(n1,2,gp,g,RSe,Yuk,Nn,fLi(2),fRi(2))

  If (n1.Eq.n2) Then
   fLj = fLi
   fRj = fRi
  Else
   Call CoupNeutralinoSlepton(n2,1,gp,g,RSe,Yuk,Nn,fLj(1),fRj(1))
   Call CoupNeutralinoSlepton(n2,2,gp,g,RSe,Yuk,Nn,fLj(2),fRj(2))
  End If


  s2 = s**2

  lambda2 = (s - mi2 - mj2)**2 - 4._dp*mi2*mj2
  lambda = Sqrt( lambda2 )

  partZZ = ( PolFactorZ(1) * Le**2 + PolFactorZ(2) * Re**2 )        &
     &    * lambda * ( Abs(Ol_ij)**2 * AngularFactorZ(1)            &
     &                 * ( 8._dp * s2 - 4._dp * (mi2+mj2) * s       &
     &                   - 4._dp * (mi2-mj2)**2                     &
     &                   + ((mi-mj)**2 - s) * ((mi+mj)**2 - s )     &
     &                     * AngularFactorZ(2) )                    &
     &               - Real( Ol_ij**2,dp ) * AngularFactorZ(3) * s     &
     &               ) / ( (s-mZ2)**2 + gmZ2 )

  ReZprop = (s-mZ2) / ( (s-mZ2)**2 + gmZ2 ) 
  Do i1=1,2
   LogPcT(i1) = Log( (SumMass2(i1) - s + lambda * cT1)     &
     &             / (SumMass2(i1) - s + lambda * cT2) )
   LogMcT(i1) = Log( (SumMass2(i1) - s - lambda * cT1)     &
     &             / (SumMass2(i1) - s - lambda * cT2) )

   partZsel(i1) = - ( PolFactorZ(1) * Le                                  &
     &                     * Real( Conjg(fRi(i1)) * fRj(i1) * Or_ij,dp )  &
     &              + PolFactorZ(2) * Re                                  &
     &                     * Real( Conjg(fLi(i1)) * fLj(i1) * Ol_ij,dp )) &
     &              * (0.5_dp * lambda * (SumMass2(i1) + s)               &
     &                            * AngularFactorZsel(1)                  & 
     &                - ProdMass(i1) * (LogPcT(i1) -  LogMcT(i1))         & 
     &                )                                                   &
     &            + ( PolFactorSel(1) * Le                                &
     &                     * Real( Conjg(fRi(i1)) * fRj(i1) * Ol_ij,dp )  &
     &                   + PolFactorSel(2) * Re                           &
     &                     * Real( Conjg(fLi(i1)) * fLj(i1) * Or_ij,dp )) & 
     &                 * s * ( LogPcT(i1) - LogMcT(i1) )


   partZsel(i1) = ReZprop * partZsel(i1)

   partSel2(i1) = ( Abs(fRi(i1))**2 * PolFactorSel2(4)                     &
     &                 + Abs(fLi(i1))**2 * PolFactorSel2(3) )              &
     &               * ( Abs(fRj(i1))**2 * PolFactorSel2(1)                &
     &                 + Abs(fLj(i1))**2 * PolFactorSel2(2) )              &
     &               * ( - 0.125_dp * AngularFactorZsel(1) * lambda        &
     &                 + 0.25_dp * SumMass2(i1) * LogPcT(i1)               &
     &                 - 0.5_dp * ProdMass(i1)                             & 
     &                   * ( 1._dp / (SumMass2(i1)-s+lambda*cT2)           &
     &                     - 1._dp / (SumMass2(i1)-s+lambda*cT1) )         &
     &                 ) 

  partSel2(i1) = partSel2(i1)                                               &
     &             + ( Abs(fRi(i1))**2 * PolFactorSel2(1)                   &
     &                 + Abs(fLi(i1))**2 * PolFactorSel2(2) )               &
     &               * ( Abs(fRj(i1))**2 * PolFactorSel2(4)                 &
     &                 + Abs(fLj(i1))**2 * PolFactorSel2(3) )               &
     &               * ( - 0.125_dp * AngularFactorZsel(1) * lambda         &
     &                 - 0.25_dp * SumMass2(i1) * LogMcT(i1)                &
     &                 + 0.5_dp * ProdMass(i1)                              &
     &                   * ( 1._dp / (SumMass2(i1)-s-lambda*cT2)            &
     &                     - 1._dp / (SumMass2(i1)-s-lambda*cT1) )          &
     &                 ) 

  partSel2(i1) = partSel2(i1)                                               &
     &             - 0.5_dp * Real( fRi(i1) * fLi(i1)                       & 
     &                          * Conjg( fRj(i1) * fLj(i1) ),dp )           &
     &               * PolFactorSel2(5)                                     &
     &               * ( AngularFactorZsel(1) * lambda                      & 
     &                 + ( 2._dp * ProdMass(i1) / (SumMass2(i1)-s) - s      &
     &                   + s * (s-mi2-mj2) / (SumMass2(i1)-s) )             &
     &                   * (- LogPcT(i1) + LogMcT(i1) )                     &
     &                 )                                                    &
     &              - 0.5_dp * ( PolFactorSel(1)                            &
     &                         * Real( fRi(i1)**2 * Conjg(fRj(i1)**2),dp )  &
     &                       + PolFactorSel(2)                              &
     &                         * Real( fLi(i1)**2 * Conjg(fLj(i1)**2),dp )  & 
     &                       )                                              &
     &                     * s * (- LogPcT(i1) + LogMcT(i1) )               &
     &                     / (SumMass2(i1)-s) 

      End Do 

  eeNeutralinoTree = ( partZZ - partZsel(1) - partZsel(2)   &
     &               + partSel2(1) + partSel2(2) ) / s2

  Iname = Iname - 1

 End Function eeNeutralinoTree


  Real(dp) Function eeNeutralinoISR(x)

  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'eeNeutralinoISR'

  s = smax * (1._dp - Exp(x))

  eeNeutralinoISR = Lee(beta,x) * eeNeutralinoTree(s)

  Iname = Iname - 1

 End Function eeNeutralinoISR
 

 Subroutine EpEmPseudoScalarScalar(i,j,mP0,RP0,mS0,RS0,Pm,Pp,s,ISR &
                                   &     ,sigma)
 !-----------------------------------------------------------------------
 ! calculates the production cross of pseudoscalar and scalar Higgs
 ! input:
 !  i .......... index of the pseudoscalar boson
 !  j .......... index of scalar boson
 !  mP0(i) ..... masses of the pseudoscalar bosons
 !  RP0(i,j) ... Mixing matrices of the pseudoscalar bosons
 !  mS0(i) ..... masses of the scalar bosons
 !  RS0(i,j) ... Mixing matrices of the scalar bosons
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... production cross section in fb
 ! written by Werner Porod, 2.1.00
 ! 24.10.2000: porting the code to f90
 ! 20.02.03: fixing bug: for the computation of smin the indices i and j 
 !                       had to be interchanged
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: mS0(:),RS0(:,:),mP0(:),RP0(:,:),Pm,Pp,s
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR

  Integer :: n_S0,n_P0

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmPseudoScalarScalar'

  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )

  If (n_S0.Ne.n_P0) Then
   Write(ErrCan,*) 'Error in subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'n_S0 =!= n_P0 ',n_S0,n_P0
   Call TerminateProgram
  Else If ( (i.Lt.1).Or.(i.Gt.n_P0).Or.(j.Lt.1).Or.(j.Gt.n_S0) ) Then
   Write(ErrCan,*) 'Error: in subroutine //NameOfUnit(Iname)// the'
   Write(ErrCan,*) 'production of pseudoscalar ',i,' and'
   Write(ErrCan,*) 'scalar ',j,'should be calculated. IMPOSSIBLE !'
   Call TerminateProgram
  End If

  If (n_S0.Gt.3) Then
   Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'model with n_S0 =',n_S0,' is considered.'
   Write(ErrCan,*) 'Chargino exchange is not yet included.'
  End If

  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (mP0(i)+mS0(j))**2.Ge.s) Then
   Iname = Iname - 1
   Return 
  End If
  !--------------------------------------------------------------------
  ! internal information
  !--------------------------------------------------------------------
  ind_1 = i
  ind_2 = j
  RS0_in = 0
  RP0_in = 0
  RS0_in(1:n_s0, 1:n_s0) = RS0
  RP0_in(1:n_p0, 1:n_p0) = RP0
  P_m = Pm
  P_p = Pp

  mP2 = mP0(i)**2
  mS2 = mS0(j)**2

  If (ISR) Then
   smin = (mS0(j)+mP0(i))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = eePseudoScalarScalarTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma = eePseudoScalarScalarTree(s) * ILee(beta,zmin)   &
    &     + dgauss(eePseudoScalarScalarISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = eePseudoScalarScalarTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
   sigma = 0.38939e12_dp * pi * sigma  / 3._dp

  Iname = Iname - 1

 End Subroutine EpEmPseudoScalarScalar


 Real(dp) Function eePseudoScalarScalarTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of CP-odd CP-even Higgs  production
 ! in e+e- annihilation. Is called by EpEmPseudoScalarScalar.
 ! written by Werner Porod, 5.1.00
 ! 24.10.2000: porting the code to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Complex(dp) :: coupC
  Real(dp) :: kappa, sumI, gin(3), yukin(3), Sqrts, sinW2, cosW, ae, ve &
    & , Le, Re, ZZ

!  kappa = (s-mS2-mP2)**2 - 4._dp * mS2 * mP2
  kappa = sq_kappa(s,mS2,mP2)
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)

  sinW2 = gin(1)**2 / (gin(1)**2 + gin(2)**2) 
  cosW = Sqrt(1 - sinW2)

  Call CoupPseudoScalarScalarZ(ind_1, ind_2, gin(2), cosW, RP0_in, RS0_in, coupC)
  Call CoupFermionZ(-0.5_dp,- 1._dp, gin(2), sinW2, Le, Re)
  ae = Re - Le
  ve = - (Re + Le)

  ZZ = Abs(coupC)**2 * 0.25_dp  &
     &  * ( (ve**2 + ae**2) * (1._dp - P_m*P_p) - 2._dp * ae * ve * (P_m - P_p) )

  sumI = oo16pi2 * ZZ * kappa**3 / ( (s-mZ2)**2 + gmZ2 )

  eePseudoScalarScalarTree = sumI / s**2

 End Function eePseudoScalarScalarTree


 Real(dp) Function eePseudoScalarScalarISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of Higgs Z production in
 ! e+e- annihilation. Is called by EpEmPseudoScalarScalar.
 ! for the calculation of ISR corrections
 ! written by Werner Porod, 5.1.00
 ! 24.10.2000: porting the code to f90
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  eePseudoScalarScalarISR = Lee(beta,x) * eePseudoScalarScalarTree(s)

 End Function eePseudoScalarScalarISR


 Subroutine EpEmScalarZMSSM(i,mS0,RS0,vevSM,Pm,Pp,s,ISR,sigma)
 !-----------------------------------------------------------------------
 ! calculates the production cross of scalar Higgs + Z for the MSSM
 ! and the NMSSM
 ! input:
 !  i .......... index of Higgs
 !  mS0(i) ..... Higgs masses
 !  RS0(i,j) ... Mixing matrices of Higgs Bosons
 !  vevs(i) .... vevs of the neutral scalars
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... production cross section in fb
 ! written by Werner Porod, 2.1.00
 ! 24.10.2000: porting the code to f90
 ! 26.08.2006: extension to include NMSSM
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i
  Real(dp), Intent(in) :: mS0(:),RS0(:,:),Pm,Pp,s,vevSM(2)
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmScalarZMSSM'
  n_S0 = Size(mS0)
  If ( (i.Lt.1).Or.(i.Gt.n_S0) ) Then
   Write(ErrCan,*) 'Error: in subroutine '//NameOfUnit(Iname)//' the production'
   Write(ErrCan,*) 'of Higgs ',i,' should be calculated. IMPOSSIBLE !'
   Call TerminateProgram
  End If

  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (ms0(i)+mZ)**2.Ge.s) Then
   Iname = Iname - 1
   Return 
  End If
  !--------------------------------------------------------------------
  ! internal information
  !--------------------------------------------------------------------
  ind_1 = i
  vevSM_in = vevSM
  vL_in = 0._dp
  RS0_in(1:n_S0,1:n_S0) = RS0
  P_m = Pm
  P_p = Pp

  mS2 = mS0(i)**2

  If (ISR) Then
   smin = (mS0(i)+mZ)**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = eeScalarZTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma = eeScalarZTree(s) * ILee(beta,zmin)   &
     &    + dgauss(eeScalarZISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = eeScalarZTree(s)
  End If

  !----------------------------
  ! 0.38939e12 gives fb
  !----------------------------
  sigma = 0.38939e12_dp * pi * sigma  / 3._dp

  Iname = Iname - 1

 End Subroutine EpEmScalarZMSSM


 Subroutine EpEmScalarZrp(i,mS0,RS0,vevSM,vevL,Pm,Pp,s,ISR,sigma)
 !-----------------------------------------------------------------------
 ! calculates the production cross of scalar Higgs + Z
 ! input:
 !  i .......... index of Higgs
 !  mS0(i) ..... Higgs masses
 !  RS0(i,j) ... Mixing matrices of Higgs Bosons
 !  vevs(i) .... vevs of the neutral scalars
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... production cross section in fb
 ! written by Werner Porod, 2.1.00
 ! 24.10.2000: porting the code to f90
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i
  Real(dp), Intent(in) :: mS0(:),RS0(:,:),Pm,Pp,s,vevSM(2),vevL(:)
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR

  integer :: n_vL

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmScalarZrp'

  n_S0 = Size( mS0 )
  If ( (i.Lt.1).Or.(i.Gt.n_S0) ) Then
   Write(ErrCan,*) 'Error: in subroutine '//NameOfUnit(Iname)//' the production'
   Write(ErrCan,*) 'of Higgs ',i,' should be calculated. IMPOSSIBLE !'
   Call TerminateProgram
  End If

  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (ms0(i)+mZ)**2.Ge.s) Then
   Iname = Iname - 1
   Return 
  End If
  !--------------------------------------------------------------------
  ! setting g to 1 because the running coupling is included later
  !--------------------------------------------------------------------
  n_vL = Size(vevL)
  n_S0 = 2 + n_vL 
  ind_1 = i
  vevSM_in = vevSM
  vL_in = 0._dp
  vL_in(1:n_vL) = vevL
  RS0_in(1:n_s0,1:n_s0) = RS0
  P_m = Pm
  P_p = Pp

  mS2 = mS0(i)**2

  If (ISR) Then
   smin = (mS0(i)+mZ)**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = eeScalarZTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma = eeScalarZTree(s) * ILee(beta,zmin)   &
     &    + dgauss(eeScalarZISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = eeScalarZTree(s)
  End If

  !----------------------------
  ! 0.38939e12 gives fb
  !----------------------------
  sigma = 0.38939e12_dp * pi * sigma  / 3._dp

  Iname = Iname - 1

 End Subroutine EpEmScalarZrp


 Real(dp) Function eeScalarZISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of Higgs Z production in
 ! e+e- annihilation. Is called by EpEmScalarZMSSM and 
 ! EpEmScalarZrp for the calculation of ISR corrections
 ! written by Werner Porod, 2.1.00
 ! 24.10.2000: porting the code to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  eeScalarZISR = Lee(beta,x) * eeScalarZTree(s)

 End Function eeScalarZISR


 Real(dp) Function eeScalarZTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of Higgs Z production in
 ! e+e- annihilation. Is called by EpEmScalarZMSSM and EpEmScalarZrp.
 ! written by Werner Porod, 2.1.00
 ! 24.10.2000: porting the code to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s

  Real(dp) :: kappa, sumI, gin(3), yukin(3), Sqrts, cosW2, sinW2, coup &
    & , Le, Re, ae, ve, ZZ

  kappa = (s-mS2-mZ2)**2 - 4._dp * mS2 * mZ2
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)

  cosW2 = gin(2)**2 / (gin(1)**2 + gin(2)**2)
  sinW2 = 1._dp - cosW2

  if ((n_S0.eq.2).or.(n_S0.eq.3)) then  
   Call CoupScalarZ(ind_1, gin(2), cosW2, vevSM_in, RS0_in, coup)
   coup = coup / (2._dp * gin(2) * Sqrt( Dot_product(vevSM_in, vevSM_in) ) )
  else 
   Call CoupScalarZ(ind_1, gin(2), cosW2, vevSM_in, vL_in, RS0_in, coup)
   coup = coup / (2._dp * gin(2) * Sqrt( Dot_product(vevSM_in, vevSM_in) &
        &                              + Dot_product(vL_in, vL_in)       ) )
  end if
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ae = Re - Le
  ve = - (Re + Le)

  ZZ = coup**2 * cosW2 * ( (ve**2 + ae**2) * (1._dp - P_m*P_p)       &
                         - 2._dp * ae * ve * (P_m - P_p) )

  sumI = ZZ * Sqrt(kappa) * (kappa + 12._dp * mZ2 * s) / ( (s-mZ2)**2 + gmZ2 )

  eeScalarZTree = oo16pi2 * sumI / s**2

 End Function eeScalarZTree


 Subroutine EpEmToHpHmMSSM(mSpm,RSpm,Pm,Pp,s,ISR,sigma)
 !-----------------------------------------------------------------------
 ! calculates the production cross of charged Higgs the formula for polarized
 ! beams is taken from the Ph.D. thesis of Sabine Kraml (1999).
 ! input:
 !  i,j ........ the Higgs combination
 !  mSpm(i) .... Higgs masses
 !  RSpm(i,j) .. Mixing matrices of charged Higgs
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! written by Werner Porod, 5.1.00
 ! 14.07.02: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

!  Integer, intent(in) :: i, j
  Real(dp), Intent(in) :: mSpm(2), Pm, Pp, s
  Complex(dp), Intent(in) :: RSpm(2,2)
  Logical, Intent(in) :: ISR

  Real(dp), Intent(out) :: sigma

  Real(dp) :: zmin, zmax, smin

  sigma = 0._dp

 !-------------
 ! kinematics
 !-------------
  If ( (mSpm(2)+mSpm(2))**2.Gt.s) Return 

 !--------------------------------------------------------------------
 ! transfering information to internal variables
 !--------------------------------------------------------------------
  ind_1=2
  ind_2=2
  P_m = Pm
  P_p = Pp
  RSpm_in(1:2,1:2) = RSpm

  mSfer2 = mSpm(2)**2

  If (ISR) Then
   smin = 4._dp * mSfer2(1)
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write (ErrCan,*) 'Warning from subroutine eeHiggs '
    Write (ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write (ErrCan,*) 'The result has to be taken with great care!!!'
   Endif
   beta = 2._dp * Alpha * ( Log(s/mf_l(1)**2) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = eeHiggsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  eeHiggsTree(s) * ILee(beta,zmin)  &
        & + dgauss(eeHiggsISR,zmin,zmax,1.e-3_dp)
   Endif

  Else
   sigma = eeHiggsTree(s)
  Endif

 !----------------------------
 ! 0.38939e12 gives fb
 !----------------------------
   sigma = 0.38939e12_dp * pi * sigma  / 3._dp

 End Subroutine EpEmtoHpHmMSSM


 Real(dp) Function eeHiggsTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of Higgs production in
 ! e+e- annihilation. Is called by eeHiggs and related functions.
 ! written by Werner Porod, 27.8.99
 ! 14.07.02: portation to f90
 !-----------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: s

  Complex(dp) :: coupC
  Real(dp) :: kappa, kappa3d2, sumI, sqrts, gin(3),yukin(3), alpha, sinW2 &
    & ,Le, Re, ae, ve, gg, gZ, ZZ

  kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
  kappa3d2 = kappa**1.5_dp 

 !------------------------
 ! running couplings
 !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  sinW2 = gin(1)**2 / (gin(1)**2 + gin(2)**2)
  alpha = sinW2 * gin(2)**2

  Call CoupChargedScalarZ(ind_1, ind_2, gin(2), sinW2, RSpm_in(1:2,1:2), coupC)
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ae = 2._dp * (Re - Le)
  ve = -2._dp * (Re + Le)

  gg = alpha**2 * (1._dp - P_m*P_p)
  gZ = alpha * coupC * 0.5_dp * ( ve * (1._dp - P_m*P_p) - ae * (P_m - P_p) )
  ZZ = Abs(coupC)**2 * 6.25e-2_dp * ( (ve**2 + ae**2) * (1._dp - P_m*P_p)  &
            - 2._dp * ae * ve * (P_m - P_p) )

  sumI = gg + ( gZ * (s-mZ2) + ZZ * s ) * s / ( (s-mZ2)**2 + gmZ2 )

  eeHiggsTree = oo16pi2 * kappa3d2 * sumI / s**4

 End Function eeHiggsTree


 Real(dp) Function eeHiggsISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeHiggs and related functions.
 ! for the calculation of ISR corrections
 ! written by Werner Porod, 26.9.99
 ! last change: 26.9.99
 !-----------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  eeHiggsISR = Lee(beta,x) *eeHiggsTree(s)

 End Function eeHiggsISR


 Subroutine EpEmToEsneutrinosMSSM( mSn, mC, V, Pm, Pp, s, ISR, Beam &
                               & , sigma, Design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of electron sneutrinos.
 ! input:
 !  mSn ........ sneutrino mass
 !  mC(i) ...... chargino masses
 !  V(i,j) ..... chargino mixing matrix V
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod, 14.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: mSn, Pm, Pp, s, mC(2)
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: V(2,2)
  Logical, Intent(in) :: ISR, Beam
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn, n_c
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToESneutrinosMSSM'
  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  mSnu2 = mSn**2
  mSn2 = 0._dp   ! this is needed for the sneutrino functions to check
                 ! whether this is the 1-generation or 3-generation model
  If ( (4._dp*mSnu2).Ge.s) Then
   iname = iname - 1
   Return 
  End If

  !--------------------------------------------------------------------
  ! transfering information to internal variables
  !--------------------------------------------------------------------
  n_c = size(mC)
  V_in(1:n_c, 1:n_c) = V
  P_m = Pm
  P_p = Pp
  ml2 = mC**2

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   smin = 4._dp *  mSnu2 !(msf(i)+msf(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = 2._dp * mSn / Ebeam - 1._dp ! (msf(i)+msf(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   itmx = 10
   nprn = -1
   ncall = 15000
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToSneutrinosBeam,init,ncall,itmx,nprn,1.e-3_dp,erg &
             &,sd,chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = 4._dp*mSnu2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)//' '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToEsneutrinosTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToEsneutrinosTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToEsneutrinosISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToEsneutrinosTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * sigma
  Iname = Iname - 1

 End Subroutine EpEmToEsneutrinosMSSM

#ifdef GENERATIONMIXING
 Subroutine EpEmToEsneutrinosMSSM3(i, j, mSn, Rsnu,  RfL, mC, V, Pm, Pp &
                               & , s, ISR, Beam, sigma, Design)
 !------------------------------------------------------------------------
 ! calculates the production cross of sneutrinos in the 3 generation case.
 ! input:
 !  i,j ........... indices of sneutrinos
 !  mSn(i) ........ sneutrino mass
 !  Rsnu(i,j) ..... sneutrino mixing matrix
 !  mC(i) ...... chargino masses
 !  V(i,j) ..... chargino mixing matrix V
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod, 17.07.01
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mSn(:), Pm, Pp, s, mC(2)
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rsnu(:,:), RfL(:,:), V(2,2)
  Logical, Intent(in) :: ISR, Beam
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: i1, init, itmx, ncall, nprn, n_c
  Real(dp) :: erg, region(8), Ebeam, fac(2), testC
  Real(dp) :: Le, Re, sinW2, cosW2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToESneutrinosMSSM3'
  sigma = 0._dp

  !-------------
  ! kinematics
  !-------------
  If ( (mSn(i)+mSn(j))**2.Ge.s) Then
   iname = iname - 1
   Return 
  End If

  mSn2 = mSn**2
  ind_1 = i
  ind_2 = j
  n_c = size(mC)
  V_in(1:n_c, 1:n_c) = V
  P_m = Pm
  P_p = Pp
  ml2 = mC**2
  RSf_in = 0._dp
  RSf_in(1:3,1:3) = RSnu
  !--------------------------------------------------------------------
  ! couplings
  ! setting g to 1 because the running coupling is included later
  !--------------------------------------------------------------------
  cosW2 = mW2 / mZ2
  sinW2 = 1._dp - cosW2
  If (i.Eq.j) Then 
   Call CoupFermionZ(-0.5_dp,-1._dp,1._dp,sinW2,Le,Re)
   fzz = oo64pi * (Le**2 *(1-Pm)*(1+Pp) +Re**2 *(1+Pm)*(1-Pp)) / (6._dp *cosW2)
   fac = 0._dp
   Do i1=1,3
    fac(1) = fac(1) + Abs( Conjg( RSnu(i,i1) ) * RfL(1,i1) )**2 
   End Do
   fac(2) = fac(1)
   Do i1=1,2 
    fZV(i1) = Abs( V(i1,1) )**2
   End Do
   fZV = - oo32pi * Le * (1-Pm) * (1+Pp) * fac(1) * fZV / Sqrt(cosW2)
  Else
   fZZ = 0._dp
   fZV = 0._dp
   fac = 0._dp
   Do i1=1,3
    fac(1) = fac(1) + Abs( Conjg( RSnu(i,i1) ) * RfL(1,i1) )**2 
    fac(2) = fac(2) + Abs( Conjg( RSnu(j,i1) ) * RfL(1,i1) )**2 
   End Do
  End If
  fVV(1,1) = oo64pi * Abs( V(1,1) )**4 
  fVV(1,2) = oo32pi * Abs( V(1,1) )**2 * Abs( V(2,1) )**2
  fVV(2,1) = fVV(1,2)
  fVV(2,2) = oo64pi * Abs( V(2,1) )**4
  fVV = fVV * (1-Pm) * (1+Pp) * fac(1) * fac(2)

!  If all couplings are 0, there is no need for further calculation
  testC = Abs(fZZ) + Sum( Abs(fZV) ) + Sum( Abs( fVV ) )
  If (testC.Eq.0._dp) Then
   iname = iname - 1
   Return 
  End If

  ml2 = mC**2
! for testing
!  ml2(1) = 38.1480675_dp**2
!  ml2(2) = 204.997101_dp**2

  If (Beam) Then
   smin = (mSn(i)+mSn(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (mSn(i)+mSn(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToSneutrinosBeam,init,ncall,itmx,nprn,1.e-3_dp,erg)
   sigma = erg * zmax**2

  Else If (ISR) Then
   smin = (mSn(i)+mSn(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)//' '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToEsneutrinosTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToEsneutrinosTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToEsneutrinosISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToEsneutrinosTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * sigma
  Iname = Iname - 1

 End Subroutine EpEmToEsneutrinosMSSM3
#endif


 Real(dp) Function EpEmToEsneutrinosTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sneutrino production in
 ! e+e- annihilation. Is called by eeEsneutrinos and related functions.
 ! written by Werner Porod, 27.8.99
 ! change to f90 standard, 1.10.2000
 ! 03.03.06: adding generation mixing. It is assumed that we are in the
 ! basis where the charged lepton Yukawa coupling is diagonal and real.
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Integer :: i1, i2
  Real(dp) :: kappa, kappa3d2,sqrts,gin(3),yukin(3), &
        &   propN(2),kappa1d2,sumGauge2,sumGaugeN,sumN2,invPropZ,   &
        &   propNN(2,2),logNN(2),sumMkMsf(2), prodSmn(2), mSn2_s
  Real(dp) :: sinW2, cosW2, Le, Re, fzz, fZV(2), fVV(2,2), fac(2)

  If (GenerationMixing) then
   if (ind_1.eq.ind_2) then
    kappa = s*(s-4._dp*mSn2(ind_1))
    prodSmn(1) = (ml2(1)-mSn2(ind_1))**2
    prodSmn(2) = (ml2(2)-mSn2(ind_1))**2
   else
    kappa = (s-mSn2(ind_1)-mSn2(ind_2))**2 - 4._dp * mSn2(ind_1) * mSn2(ind_2)
    prodSmn(1) = (ml2(1)-mSn2(ind_1)) * (ml2(1)-mSn2(ind_2))
    prodSmn(2) = (ml2(2)-mSn2(ind_1)) * (ml2(2)-mSn2(ind_2))
   end if
  else
   mSnu2 = mSn2(1)
   kappa = s*(s-4._dp*mSnu2)
   prodSmn(1) = (ml2(1)-mSnu2)**2
   prodSmn(2) = (ml2(2)-mSnu2)**2
  end if

  kappa1d2 = Sqrt(kappa)
  kappa3d2 = kappa**1.5_dp

  If (GenerationMixing) then
   mSn2_s = mSn2(ind_1) + mSn2(ind_2)
   Do i1=1,2
    logNN(i1) = Log( (s+2._dp*ml2(i1)+kappa1d2-mSn2_s ) /   &
          &          (s+2._dp*ml2(i1)-kappa1d2-mSn2_s ) )
    sumMkMsf(i1) = s*ml2(i1) +  prodSmn(i1)
    propN(i1) = 0.5_dp * (mSn2_s-s-2._dp * ml2(i1) )*kappa1d2 &
            & + sumMkMsf(i1) * logNN(i1)
   End Do

   propNN(1,1) = -2._dp * kappa1d2 + (s+2._dp*ml2(1)-mSn2_s ) * logNN(1)
   If (ml2(1).Eq.ml2(2)) Then
    propNN(1,2) = propNN(1,1)
    propNN(2,2) = propNN(1,1)
   Else
    propNN(2,2) = -2._dp * kappa1d2 + (s+2._dp*ml2(2)-mSn2_s ) * logNN(2)
    propNN(1,2) = - kappa1d2                                          &
          &   + (sumMkMsf(1)* logNN(1)-sumMkMsf(2)* logNN(2)) / (ml2(1)-ml2(2))
   End If

  else ! no generation mixing
   Do i1=1,2
    logNN(i1) = Log( (s+2._dp*ml2(i1)+kappa1d2-2._dp*mSnu2 ) /   &
          &          (s+2._dp*ml2(i1)-kappa1d2-2._dp*mSnu2 ) )
    sumMkMsf(i1) = s*ml2(i1) +  prodSmn(i1)
    propN(i1) = (mSnu2-0.5_dp*s-ml2(i1) )*kappa1d2 + sumMkMsf(i1) * logNN(i1)
   End Do

   propNN(1,1) = -2._dp * kappa1d2 + (s+2._dp*ml2(1)-2._dp*mSnu2 ) * logNN(1)
   If (ml2(1).Eq.ml2(2)) Then
    propNN(1,2) = propNN(1,1)
    propNN(2,2) = propNN(1,1)
   Else
    propNN(2,2) = -2._dp * kappa1d2 + (s+2._dp*ml2(2)-2._dp*mSnu2 ) * logNN(2)
    propNN(1,2) = - kappa1d2                                          &
          &   + (sumMkMsf(1)* logNN(1)-sumMkMsf(2)* logNN(2)) / (ml2(1)-ml2(2))
   End If
  end if

  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)

  cosW2 = gin(2)**2 / (gin(1)**2 + gin(2)**2)
  sinW2 = 1._dp - cosW2

  If ((.not.GenerationMixing).or.(ind_1.eq.ind_2)) then
   Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
   fzz = oo64pi * gin(2)**2 * (Le**2 *(1-P_m)*(1+P_p) +Re**2 *(1+P_m)*(1-P_p))  &
       &        / (6._dp * cosW2)
   if (GenerationMixing) then
    fac(1) = Abs(Rsf_in(ind_1,1))**2
    fac(2) = fac(1)
   else 
    fac = 1._dp
   end if
   Do i1=1,2
    fZV(i1) = Abs( V_in(i1,1) )**2
   End Do
   fZV = - oo32pi * gin(2)**3 * Le * (1-P_m) * (1+P_p) * fac(1) * fZV  &
       &          / Sqrt(cosW2)

  else
   fZZ = 0._dp
   fZV = 0._dp
   fac(1) = Abs(Rsf_in(ind_1,1))**2
   fac(2) = Abs(Rsf_in(ind_2,1))**2
  end if

  fVV(1,1) = oo64pi * Abs( V_in(1,1) )**4
  fVV(1,2) = oo32pi * Abs( V_in(1,1) )**2 * Abs( V_in(2,1) )**2
  fVV(2,1) = fVV(1,2)
  fVV(2,2) = oo64pi * Abs( V_in(2,1) )**4
  fVV = gin(2)**4 * fVV * (1-P_m) * (1+P_p) * fac(1) * fac(2)

  If (fZZ.Ne.0._dp) Then
   invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )
   sumGauge2 = fZZ * kappa3d2 * invPropZ  

   sumGaugeN = (s-mZ2) * invPropZ * ( fZV(1) * propN(1)  &
             &                      + fZV(2) * propN(2)  )
  Else
   sumGaugeN = 0._dp
   sumGauge2 = 0._dp
  End If 

  sumN2 = 0._dp
  Do i1=1,2
   Do i2=i1,2
    sumN2 = sumN2 + fVV(i1,i2)*propNN(i1,i2)
   End Do
  End Do

  EpEmToEsneutrinosTree = (sumGauge2+sumGaugeN+sumN2) / s**2

 End Function EpEmToEsneutrinosTree


 Real(dp) Function EpEmToEsneutrinosISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeEsneutrinos and related functions.
 ! written by Werner Porod, 27.8.99
 ! change to f90 standard, 1.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  EpEmToEsneutrinosISR = Lee(beta,x) * EpEmToEsneutrinosTree(s)

 End Function EpEmToEsneutrinosISR


#ifdef BEAMSTRAHLUNG
 Real(dp) Function EpEmToSneutrinosBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sneutrino production in
 ! e+e- annihilation. Is called by eeSelectrons and related functions.
 ! written by Werner Porod, 16.07.01
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 
  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToSneutrinosBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToSneutrinosBeam = BeamFactor * ISRfactor2 * sigmaT

 End Function EpEmToSneutrinosBeam
#endif


 Subroutine EpEmToSelectronsMSSM(i, j, mSf, Rsf, mN, N, Pm, Pp, s, ISR, Beam &
                               &, sigma, Design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of selectrons
 ! input:
 !  i,j ........ the selectron combination
 !  mSf(i) ..... selectron masses
 !  Rsf(i,j) ... Mixing matrices of selectrons
 !  mN(i) ...... Neutralino masses
 !  N(i,j) ..... Neutralino mixing matrix
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod, 14.10.2000
 ! 15.1.01: adding beamstrahlung
 ! 26.08.2006: extension to include NMSSM, is rather trival as only the
 !             dimensions of N and mN are changed
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: mSf(2), Pm, Pp, s, mN(:)
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rsf(2,2), N(:,:)
  Logical, Intent(in) :: ISR, Beam
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToSelectronsMSSM'
  sigma = 0._dp

  If ( (i.Lt.1).Or.(i.Gt.2).Or.(j.Lt.1).Or.(j.Gt.2) ) Then
   Write(ErrCan,*) 'Error: in subroutine'//NameOfUnit(Iname)//'the combination'
   Write(ErrCan,*) i,j,' should be calculated. IMPOSSIBLE !!!!'
   Call TerminateProgram
  End If

  !-------------
  ! kinematics
  !-------------
  If ( (msf(i)+msf(j))**2.Ge.s) Then
   iname = iname - 1
   Return 
  End If

  !--------------------------------------------------------------------
  ! internal variables
  !--------------------------------------------------------------------
  n_n = Size(mN)
  ind_1 = i
  ind_2 = j
  P_m = Pm
  P_p = Pp
  Rsf_in = 0
  Rsf_in(1:2,1:2) = Rsf
  N_in(1:n_n,1:n_n) = N
  mK(1:n_n) = mN
  mk2(1:n_n) = mN**2    
  mSfer2(1) = msf(i)**2
  mSfer2(2) = msf(j)**2

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   ! due to polarization it is possible that all couplings vanish
   ! check therefore first, if tree level is non-zero
   sigma = EpEmToSelectronsTree(s)
   If (sigma.Eq.0._dp) Then
    Iname = Iname - 1
    Return
   End If
   ! now do the real calculation

   smin = (msf(i)+msf(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

!   zmax = Log(1._dp - smin / s)
!   zmax = Log(1._dp - min(msf(i),msf(j)) / sqrt(s))
!   zmin = - 0.01_dp
!   zmin = min(msf(i),msf(j)) / sqrt(s)
!   zmax = 1
   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (msf(i)+msf(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToSelectronsBeamA,init,ncall,itmx,nprn,1.e-3_dp,erg &
              &, sd, chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = (msf(i)+msf(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)//' '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToSelectronsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToSelectronsTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToSelectronsISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToSelectronsTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * sigma  

  Iname = Iname - 1

 End Subroutine EpEmToSelectronsMSSM


#ifdef GENERATIONMIXING
 Subroutine EpEmToSleptonsMSSM(i, j, mSf, Rsf, Y_l, RlL, RlR, mN, N &
                               &, Pm, Pp, s, ISR, Beam, sigma, Design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of sleptons in the 3-generation MSSM
 ! input:
 !  i,j ........ the selectron combination
 !  mSf(i) ..... slepton masses
 !  Rsf(i,j) ... Mixing matrices of sleptons
 !  Y_l(i,j) ... lepton Yukawa coupling divided by g_SU(2)
 !  RlL(i,j) ... left lepton mixing matrix
 !  RlR(i,j) ... right lepton mixing matrix
 !  mN(i) ...... Neutralino masses
 !  N(i,j) ..... Neutralino mixing matrix
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod, 10.05.2001
 ! 26.08.2006: extension to include NMSSM, is rather trival as only the
 !             dimensions of N and mN are changed
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: mSf(6), Pm, Pp, s, mN(:)
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rsf(6,6), Y_L(3,3), RlL(3,3), RlR(3,3) &
                               &, N(:,:)
  Logical, Intent(in) :: ISR, Beam
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToSleptonsMSSM'
  sigma = 0._dp

  If ( (i.Lt.1).Or.(i.Gt.6).Or.(j.Lt.1).Or.(j.Gt.6) ) Then
   Write(ErrCan,*) 'Error: in subroutine'//NameOfUnit(Iname)//'the combination'
   Write(ErrCan,*) i,j,' should be calculated. IMPOSSIBLE !!!!'
   Call TerminateProgram
  End If

  !-------------
  ! kinematics
  !-------------
  If ( (msf(i)+msf(j))**2.Ge.s) Then
   iname = iname - 1
   Return 
  End If
  
  !--------------------------------------------------------------------
  ! information which are used later on
  !--------------------------------------------------------------------
  n_n = Size(mN)
  ind_1 = i
  ind_2 = j
  N_in(1:n_n,1:n_n) = N
  Rsf_in = Rsf
  P_m = Pm
  P_p = Pp
  mK = 0._dp
  mK(1:n_n) = mN
  mk2 = mK**2
  
  mSfer2(1) = msf(i)**2
  mSfer2(2) = msf(j)**2

  If (Beam) Then
   smin = (msf(i)+msf(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

!   zmax = Log(1._dp - smin / s)
!   zmax = Log(1._dp - min(msf(i),msf(j)) / sqrt(s))
!   zmin = - 0.01_dp
!   zmin = min(msf(i),msf(j)) / sqrt(s)
!   zmax = 1
   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (msf(i)+msf(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   Call Vegas1(region,EpEmToSelectronsBeamA,init,ncall,itmx,nprn,1.e-3_dp,erg &
              &, sd, chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
   smin = (msf(i)+msf(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine '//NameOfUnit(Iname)//' '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToSelectronsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToSelectronsTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToSelectronsISR,zmin,zmax,1.e-3_dp)
   End If

  Else
!   write (52,*) i,j
   sigma = EpEmToSelectronsTree(s)
  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * sigma  

  Iname = Iname - 1

 End Subroutine EpEmToSleptonsMSSM
#endif

 Real(dp) Function EpEmToSelectronsTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of selectron production in
 ! e+e- annihilation. Is called by eeSelectrons and related functions.
 ! written by Werner Porod, 14.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Integer :: i1,i2
  Real(dp) :: kappa,kappa3d2,sqrts,gin(3),yukin(3), &
        &   propN(5),kappa1d2,sumGauge2,sumGaugeN,sumN2,invPropZ,   &
        &   propNN(2,5,5),logNN(5),sumMkMsf(5)

  Real(dp) :: Le, Re, ProdPmPp(2), sinW2, cosW2
  Complex(dp) :: coup
  Complex(dp) :: ai(5), bi(5), aiC(5), biC(5), aj(5), bj(5), ajC(5) &
      &  , bjC(5), coupC, Rsf(2,2), N(n_n,n_n)
  Complex(dp), Parameter :: mat3C(3,3) = 0._dp

  If (mSfer2(1).Eq.mSfer2(2)) Then
   kappa = s*(s-4._dp*mSfer2(1))
   kappa1d2 = Sqrt(kappa)
   kappa3d2 = kappa**1.5_dp
   Do i1=1,n_n
    logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-2._dp*mSfer2(1) ) /   &
         &           (s+2._dp*mk2(i1)-kappa1d2-2._dp*mSfer2(1) ) )
    sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))**2
    propN(i1) = (mSfer2(1)-0.5_dp*s-mk2(i1) )*kappa1d2              &
         &    + sumMkMsf(i1) * logNN(i1)
   End Do
   Do i1=1,n_n
    Do i2=i1,n_n
     If (mk2(i1).Eq.mk2(i2)) Then
      propNN(1,i1,i2) = -2._dp * kappa1d2    &
       &             + (s+2._dp*mk2(i1)-2._dp*mSfer2(1) ) * logNN(i1)
      propNN(2,i1,i2) = s * kappa1d2 / sumMkMsf(i1)
     Else
      propNN(1,i1,i2) = - kappa1d2                                          &
          &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
          &             / (mk2(i1)-mk2(i2))
      propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
     End If
    End Do
   End Do
  Else
   kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
   kappa1d2 = Sqrt(kappa)
   kappa3d2 = kappa**1.5_dp
   Do i1=1,n_n
    logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-mSfer2(1)-mSfer2(2) ) /      &
         &           (s+2._dp*mk2(i1)-kappa1d2-mSfer2(1)-mSfer2(2) ) )
    sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))*(mk2(i1)-mSfer2(2))
    propN(i1) = (0.5_dp*(mSfer2(1)+mSfer2(2)-s)-mk2(i1) )*kappa1d2          &
         &    + sumMkMsf(i1) * logNN(i1)
   End Do
   Do i1=1,n_n
    Do i2=i1,n_n
     If (mk2(i1).Eq.mk2(i2)) Then
      propNN(1,i1,i2) = -2._dp * kappa1d2    &
       &             + (s+2._dp*mk2(i1)-mSfer2(1)-mSfer2(2) ) * logNN(i1)
      propNN(2,i1,i2) = s * kappa1d2 /sumMkMsf(i1)
     Else
      propNN(1,i1,i2) = - kappa1d2                                          &
          &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
          &             / (mk2(i1)-mk2(i2))
      propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
     End If
    End Do
   End Do
  End If 
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  
  cosW2 = gin(2)**2/ (gin(1)**2 + gin(2)**2)
  sinW2 = 1._dp - cosW2
  If (GenerationMixing) Then
   Call CoupSfermionZ(ind_1, ind_2, gin(2), sinW2, -1._dp, -0.5_dp, RSf_in, coup)
  Else
   Rsf = Rsf_in(1:2,1:2)
   Call CoupSfermionZ(ind_1, ind_2, gin(2), sinW2, -1._dp, -0.5_dp, RSf, coup)
  End If

  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ProdPmPp(1) = (1-P_m)*(1+P_p)
  ProdPmPp(2) = (1+P_m)*(1-P_p)

  N = N_in(1:n_n,1:n_n)

  If (ind_1.Eq.ind_2) Then
   fgg = oo16pi * gin(2)**4 * sinW2**2 * (1._dp - P_m*P_p) / 3._dp
   fgZ = -oo16pi * gin(2)**2 * sinW2 * Conjg(coup)   &
       &         * (Le*ProdPmPp(1)+Re*ProdPmPp(2)) /3._dp

!---------------------------------------
! neglecting electron Yukawa coupling
!---------------------------------------
   fZNr = 0._dp
   fZNi = 0._dp
   fgN = 0._dp
   Do i1=1,n_n
    If (GenerationMixing) Then
     Call CoupNeutralinoSlepton(1,i1,ind_1, gin(1), gin(2), RSf_in, id3C, id3C &
                              & ,mat3C, N, bi(i1), ai(i1))
    Else
     Call CoupNeutralinoSlepton(i1,ind_1, gin(1), gin(2), RSf, ZeroC, N &
        &                      , bi(i1), ai(i1))
    End If
    fgN(i1) =  Abs( ai(i1) )**2 * ProdPmPp(1) + Abs( bi(i1) )**2 * ProdPmPp(2)
    fZNr(i1) = Le * Abs( ai(i1) )**2 * ProdPmPp(1)  &
             + Re * Abs( bi(i1) )**2 * ProdPmPp(2)
    fZNi(i1) = 0._dp
   End Do
   fgN = oo16pi * gin(2)**2 * sinW2 * fgN 
   fZNr = - oo16pi * coup * fZNr
   Do i1=1,n_n
    Do i2=i1,n_n
     fNN(1,i1,i2) = ( Abs( ai(i1) )**2 *Abs( ai(i2) )**2 * ProdPmPp(1) &
         + Abs( bi(i1) )**2 *Abs( bi(i2) )**2 * ProdPmPp(2) )
     coupC = ai(i1) * bi(i2) * Conjg( ai(i2) * bi(i1) )
     fNN(2,i1,i2) = 2._dp * mk(i1) * mk(i2)                                  &
       * ( Real(coupC,dp) * (1 + P_m*P_p) + Aimag(coupC) * (P_m+P_p) )
    End Do
   End Do

  Else
   fgg = 0._dp
   fgZ = 0._dp
!---------------------------------------
! neglecting electron Yukawa coupling
!---------------------------------------
   Do i1=1,n_n
    If (GenerationMixing) Then
     Call CoupNeutralinoSlepton(1,i1,ind_1, gin(1), gin(2), RSf_in, id3C, id3C &
                              & ,mat3C, N, bi(i1), ai(i1))
     Call CoupNeutralinoSlepton(1,i1,ind_2, gin(1), gin(2), RSf_in, id3C, id3C &
                              & ,mat3C, N, bj(i1), aj(i1))
    Else
     Call CoupNeutralinoSlepton(i1,ind_1, gin(1), gin(2), RSf, ZeroC &
       &                      , N, bi(i1), ai(i1))
     Call CoupNeutralinoSlepton(i1,ind_2, gin(1), gin(2), RSf, ZeroC &
       &                      , N, bj(i1), aj(i1))
    End If
    fgN(i1) = 0._dp
    aiC(i1) = Conjg( ai(i1) )
    biC(i1) = Conjg( bi(i1) )
    ajC(i1) = Conjg( aj(i1) )
    bjC(i1) = Conjg( bj(i1) )
    coupC = Le * aiC(i1) * aj(i1) * ProdPmPp(1) &
          + Re * biC(i1) * bj(i1) * ProdPmPp(2)
    coupC = coup * coupC
    fZNr(i1) = Real( coupC,dp )
    fZNi(i1) = Aimag( coupC ) 
   End Do
   fZNr = - oo16pi * fZNr
   fZNi = - oo16pi * gmZ * fZNi

   Do i1=1,n_n
    Do i2=i1,n_n
     fNN(1,i1,i2) = Real( ai(i1) *aiC(i2) *ajC(i1) *aj(i2) * ProdPmPp(1)    &
                  &     + bi(i1) *biC(i2) *bjC(i1) *bj(i2) * ProdPmPp(2),dp )
     fNN(2,i1,i2) = mk(i1) * mk(i2)                                            &
                  & * ( aj(i2) *ajC(i1) *bi(i1) *biC(i2) *(1 - P_m) *(1 - P_p) &
                  &   + ai(i1) *aiC(i2) *bj(i2) *bjC(i1) *(1 + P_m) *(1 + P_p) )
    End Do
   End Do
  End If

  Do i1=1,n_n
   Do i2=i1,n_n
    If (i1.Eq.i2) Then
     fNN(1,i1,i2) = oo64pi * fNN(1,i1,i2) 
     fNN(2,i1,i2) = oo64pi * fNN(2,i1,i2) 
    Else
     fNN(1,i1,i2) = oo32pi * fNN(1,i1,i2) 
     fNN(2,i1,i2) = oo32pi * fNN(2,i1,i2) 
    End If
   End Do
  End Do
  fzz = oo16pi * Abs(coup)**2 * (Le**2 *ProdPmPp(1)+Re**2 *ProdPmPp(2)) / 6._dp

  invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )
  sumGauge2 = ( fgg + ( fgZ * (s-mZ2) + fZZ * s ) * s * invPropZ ) * kappa3d2
  sumGauge2 = sumGauge2 / s**4 

  sumGaugeN = 0._dp
  sumN2 = 0._dp
  Do i1=1,n_n
   sumGaugeN = sumGaugeN &
         &   + (fgN(i1)+ s*((s-mZ2)*fZNr(i1)+fZNi(i1))*invPropZ)*propN(i1)
   Do i2=i1,n_n
    sumN2 = sumN2 + ( fNN(1,i1,i2)*propNN(1,i1,i2)   &
      &             + fNN(2,i1,i2)*propNN(2,i1,i2) )
   End Do
  End Do
  sumGaugeN =  sumGaugeN / s**3
  sumN2 = sumN2 / s**2

  EpEmToSelectronsTree = (sumGauge2+sumGaugeN+sumN2)

 End Function EpEmToSelectronsTree


 Real(dp) Function EpEmToSelectronsISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of selectron production in
 ! e+e- annihilation. Is called by eeSelectrons and related functions.
 ! written by Werner Porod, 14.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  EpEmToSelectronsISR = Lee(beta,x) * EpEmToSelectronsTree(s)

 End Function EpEmToSelectronsISR


#ifdef BEAMSTRAHLUNG
 Real(dp) Function EpEmToSelectronsBeamA(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of selectron production in
 ! e+e- annihilation. Is called by eeSelectrons and related functions.
 ! written by Werner Porod, 15.1.01
 ! - 23.05.2001: trying to use a polynomial interpolation to save
 !               computation time for the tree-cross section
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i1,i2, i_run, i_near
  Real(dp) :: kappa,kappa3d2,sqrts,gin(3),yukin(3),s, &
        &   propN(4),kappa1d2,sumGauge2,sumGaugeN,sumN2,invPropZ,   &
        &   propNN(2,4,4),logNN(4),sumMkMsf(4), x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp !/ (1._dp - ymin)
  Else
!   x1 = 1._dp - (1._dp - z(1) * (1._dp - ymin) / radprob)**3  
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If

  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
!   BeamFactor = BeamFactor  / (1._dp - ymin)
  Else If (z(1).Eq.z(3)) Then
   BeamFactor = BeamFactor**2

  Else
!   x3 = 1._dp - (1._dp - z(3) * (1._dp - ymin)  / radprob)**3  
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
!  x1 = 1._dp
!  x3 = 1._dp
!  BeamFactor = 1._dp 
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------

  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  If (z(2).Eq.z(4)) Then
   ISRfactor2 = ISRfactor2**2
  Else
   x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
   ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)
  End If

!  x2 = 1._dp ! z(2)
!  x4 = 1._dp ! z(4)
!  ISRfactor2 = 1._dp
  
  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 
  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToSelectronsBeamA = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    If (mSfer2(1).Eq.mSfer2(2)) Then
     kappa = s*(s-4._dp*mSfer2(1))
     kappa1d2 = Sqrt(kappa)
     kappa3d2 = kappa**1.5_dp
     Do i1=1,4
      logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-2._dp*mSfer2(1) ) /   &
           &           (s+2._dp*mk2(i1)-kappa1d2-2._dp*mSfer2(1) ) )
      sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))**2
      propN(i1) = (mSfer2(1)-0.5_dp*s-mk2(i1) )*kappa1d2              &
           &    + sumMkMsf(i1) * logNN(i1)
     End Do
     Do i1=1,4
      Do i2=i1,4
       If (mk2(i1).Eq.mk2(i2)) Then
        propNN(1,i1,i2) = -2._dp * kappa1d2    &
         &             + (s+2._dp*mk2(i1)-2._dp*mSfer2(1) ) * logNN(i1)
        propNN(2,i1,i2) = s * kappa1d2 / sumMkMsf(i1)
       Else
        propNN(1,i1,i2) = - kappa1d2                                          &
            &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
            &             / (mk2(i1)-mk2(i2))
        propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
       End If
      End Do
     End Do
    Else
     kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
     kappa1d2 = Sqrt(kappa)
     kappa3d2 = kappa**1.5_dp
     Do i1=1,4
      logNN(i1) = Log( (s+2._dp*mk2(i1)+kappa1d2-mSfer2(1)-mSfer2(2) ) /      &
           &           (s+2._dp*mk2(i1)-kappa1d2-mSfer2(1)-mSfer2(2) ) )
      sumMkMsf(i1) = s*mk2(i1)+(mk2(i1)-mSfer2(1))*(mk2(i1)-mSfer2(2))
      propN(i1) = (0.5_dp*(mSfer2(1)+mSfer2(2)-s)-mk2(i1) )*kappa1d2          &
           &    + sumMkMsf(i1) * logNN(i1)
     End Do
     Do i1=1,4
      Do i2=i1,4
       If (mk2(i1).Eq.mk2(i2)) Then
        propNN(1,i1,i2) = -2._dp * kappa1d2    &
         &             + (s+2._dp*mk2(i1)-mSfer2(1)-mSfer2(2) ) * logNN(i1)
        propNN(2,i1,i2) = s * kappa1d2 /sumMkMsf(i1)
       Else
        propNN(1,i1,i2) = - kappa1d2                                          &
            &           + (sumMkMsf(i1)* logNN(i1)-sumMkMsf(i2)* logNN(i2))   &
            &             / (mk2(i1)-mk2(i2))
        propNN(2,i1,i2) = s * (logNN(i2)-logNN(i1)) / (mk2(i1)-mk2(i2))
       End If
      End Do
     End Do
    End If 
    !------------------------
    ! running couplings
    !------------------------
    Sqrts = Sqrt(s)
    Call runningCouplings(Sqrts,gin,yukin)

    invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )
    sumGauge2 = ( fgg + ( fgZ * (s-mZ2) + fZZ * s ) * s * invPropZ ) * kappa3d2
    sumGauge2 = sumGauge2 / s**4

    sumGaugeN = 0._dp
    sumN2 = 0._dp
    Do i1=1,4
     sumGaugeN = sumGaugeN &
           &   + (fgN(i1)+ s*((s-mZ2)*fZNr(i1)+fZNi(i1))*invPropZ)*propN(i1)
     Do i2=i1,4
      sumN2 = sumN2 + ( fNN(1,i1,i2)*propNN(1,i1,i2)   &
        &             + fNN(2,i1,i2)*propNN(2,i1,i2) )
     End Do
    End Do
    sumGaugeN =  sumGaugeN / s**3
    sumN2 = sumN2 / s**2
    Spline(i_run,1) = s
    Spline(i_run,2) = gin(2)**4 * (sumGauge2+sumGaugeN+sumN2)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToSelectronsBeamA = BeamFactor * ISRfactor2 * sigmaT


 End Function EpEmToSelectronsBeamA
#endif


 Subroutine EpEmToSfermionsZG(i, j, gen, specie, mSf, Rsf, Pm, Pp, s, ISR &
                            &, Beam, sigma, design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of sfermions, the formula for polarized
 ! beams is taken from the Ph.D. thesis of Sabine Kraml (1999).
 ! If QCD corrections should be calculated, EpEmToSquarksZG must be used.
 ! input:
 !  i,j ........ the sfermion combination
 !  gen ........ the generation, needed as warning in case of selectrons
 !  specie ..... type of sfermion: 'u-squark', 'd-squark', 'slepton'
 !                                 'sneutrino'
 !  mSf(i) ..... sfermion masses
 !  Rsf(i,j) ... Mixing matrices of sfermions
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 !  Beam ....... Beam and ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod: 24.10.2000
 ! - 05.09.2001: adding beamstrahlung
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j,gen
  Real(dp), Intent(in) :: mSf(2),Pm,Pp,s
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rsf(2,2)
  Logical, Intent(in) :: ISR, Beam
  Character (len=9), Intent(in) :: specie
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToSfermionsZG'
  sigma = 0._dp

  If ( (i.Lt.1).Or.(i.Gt.2).Or.(j.Lt.1).Or.(j.Gt.2) ) Then
   Write(ErrCan,*) 'Error: in subroutine EpEmToSfermionsZG the combination'
   Write(ErrCan,*) i,j,' should be calculated. IMPOSSIBLE !!!!'
   Call TerminateProgram
  End If

  If ((gen.Lt.1).Or.(gen.Gt.3)) Then
   Write(ErrCan,*) 'Error: in subroutine EpEmToSfermionsZG the'
   Write(ErrCan,*) 'Generation ',gen,' does not exist.'
   Call TerminateProgram
  End If

  !-------------
  ! kinematics
  !-------------
  If ( (msf(i)+msf(j))**2.Ge.s) Then
   Iname = Iname - 1
   Return 
  End If
  !---------------------------
  ! using species informaiton
  !---------------------------
  Select Case (specie)
   Case ('u-squark')
     T3_in = 0.5_dp
     e_in = 2._dp / 3._dp
   Case ('d-squark')
     T3_in = -0.5_dp
     e_in = -1._dp / 3._dp
   Case ('slepton')
     T3_in = -0.5_dp
     e_in = -1._dp 
   Case ('sneutrino')
     T3_in = 0.5_dp
     e_in = 0._dp 
   Case Default
    Write(ErrCan,*) 'Error in Subroutine EpEmToSfermionsZG.'
    Write(ErrCan,*) 'Specie '//specie//' not include.'
    Call TerminateProgram 
  End Select
  
  !-------------------------
  ! some warnings
  !-------------------------
  If ( (e_in.Eq.0._dp).And.( (i.Eq.2).Or.(j.Eq.2) ) ) Then
   Write(ErrCan,*) 'Warning from Subroutine EpEmToSfermionsZG.'
   Write(ErrCan,*) 'The production of right handed sneutrinos'
   Write(ErrCan,*) 'has been calculated. i,j = ',i,j
  End If

  If ( (gen.Eq.1).And.((e_in.Eq.0._dp).Or.(e_in.Eq.-1._dp)) ) Then
   Write(ErrCan,*) 'Warning: in subroutine EpEmToSfermionsZG for first '
   Write(ErrCan,*) 'generation sleptons only Z and photon is included!!!'
  End If

  !--------------------------------------------------------------------
  ! internal variables
  !--------------------------------------------------------------------
  ind_1 = i
  ind_2 = j
  P_m = Pm
  P_p = Pp
  Rsf_in = 0
  Rsf_in(1:2,1:2) = Rsf

  mSfer2(1) = msf(i)**2
  mSfer2(2) = msf(j)**2

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   smin = (msf(i)+msf(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (msf(i)+msf(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToSfermionsBeam,init,ncall,itmx,nprn,1.e-3_dp,erg &
             &,sd,chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = (msf(i)+msf(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine EpEmToSfermionsZG '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToSfermionsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToSfermionsTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToSfermionsISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToSfermionsTree(s)
  End If

  !----------------------------
  ! 0.38939e12 gives fb
  !----------------------------
  If ((e_in.Eq.0._dp).Or.(e_in.Eq.-1._dp)) Then
   sigma = 0.38939e12_dp * pi * sigma  / 3._dp
  Else
   sigma = 0.38939e12_dp * pi * sigma  
  End If

  Iname = Iname - 1

 End Subroutine EpEmToSfermionsZG


 Subroutine EpEmToSquarksZG(i, j, specie, mSf, Rsf, mq, mglu, Pm, Pp, s, ISR &
                          &, Beam , sigma, Design)
 !-----------------------------------------------------------------------
 ! calculates the production cross of squarks including QCD corrections, 
 ! the formula for polarized beams is taken from the Ph.D. thesis of 
 ! Sabine Kraml (1999).
 ! input:
 !  i,j ........ the sfermion combination
 !  specie ..... string, specifying the isospin 'u-squark' or 'd-squark'
 !  mSf(i) ..... sfermion masses
 !  Rsf(i,j) ... Mixing matrices of sfermions
 !  mq ......... mass of the corresponding quark
 !  mglu ....... gluino mass
 !  Pm ......... degree of e- polarisation
 !  Pp ......... degree of e+ polarisation
 !  s .......... c.m.s. energy squared
 !  ISR ........ ISR corrections are included if .TRUE.
 !  Beam ....... Beam and ISR corrections are included if .TRUE.
 ! output:
 !  sigma ...... cross section in fb
 ! written by Werner Porod: 24.10.2000
 ! - 05.09.2001: adding beamstrahlung
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i,j
  Real(dp), Intent(in) :: mSf(:), Pm, Pp, s, mq, mglu
  Real(dp), Intent(out) :: sigma
  Complex(dp), Intent(in) :: Rsf(:,:)
  Logical, Intent(in) :: ISR, Beam
  Character (len=9), Intent(in) :: specie
  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, chi2a, sd, region(8), Ebeam

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToSquarksZG'
  sigma = 0._dp

  n_sf = Size(mSf)
  If ( (i.Lt.1).Or.(i.Gt.n_sf).Or.(j.Lt.1).Or.(j.Gt.n_sf) ) Then
   Write(ErrCan,*) 'Error: in subroutine EpEmToSquarksZG the combination'
   Write(ErrCan,*) i,j,' should be calculated. IMPOSSIBLE !!!!'
   Call TerminateProgram
  End If

  !-------------
  ! kinematics
  !-------------
  If ( (msf(i)+msf(j))**2.Ge.s) Then
   Iname = Iname - 1
   Return 
  End If
  !---------------------------
  ! using species informaiton
  !---------------------------
  Select Case (specie)
   Case ('u-squark')
     T3_in = 0.5_dp
     e_in = 2._dp / 3._dp
   Case ('d-squark')
     T3_in = -0.5_dp
     e_in = -1._dp / 3._dp
   Case Default
    Write(ErrCan,*) 'Error in Subroutine EpEmToSquarkssZG.'
    Write(ErrCan,*) 'Specie '//specie//' not include.'
    Call TerminateProgram 
  End Select
  
  !--------------------------------------------------------------------
  ! internal variables
  !--------------------------------------------------------------------
  icase = i + j - 1
  msq1 = mSf(1)
  msq2 = mSf(2)
  msq12 = msq1**2
  msq22 = msq2**2
  costh = Rsf(1,1)
  sinth = Rsf(1,2)
  cos2th = Abs(costh)**2 - Abs(sinth)**2
  sin2th = 2._dp * Real( costh * sinth, dp )
  Qq = e_in

  a(1,2) = - sin2th
  a(2,1) = - sin2th
  mglu2 = mglu**2
  mglumq = mglu * mq
  mq2 = mq**2

  ind_1 = i
  ind_2 = j
  P_m = Pm
  P_p = Pp
  RSf_in = 0._dp
  RSf_in(1:n_sf,1:n_sf) = RSf

  mSfer2(1) = msf(i)**2
  mSfer2(2) = msf(j)**2

#ifdef BEAMSTRAHLUNG
  If (Beam) Then
   smin = (msf(i)+msf(j))**2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = (msf(i)+msf(j)) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1 
   CalculateSpline = .True.
  Call Vegas1(region,EpEmToSfermionsQCDpBeam,init,ncall,itmx,nprn,1.e-3_dp &
             &,erg,sd,chi2a)
   sigma = erg * zmax**2

  Else If (ISR) Then
#else
  If (ISR) Then
#endif
   smin = (msf(i)+msf(j))**2
   smax = s
   If ((smin.Lt.(mZ2+25._dp)).And.(smax.Gt.(mZ2-25._dp))) Then
    Write(ErrCan,*) 'Warning from subroutine EpEmToSfermionsZG '
    Write(ErrCan,*) 'Inclusion of ISR corrections near m_Z '
    Write(ErrCan,*) 'The result has to be taken with great care!!!'
   End If
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToSfermionsQCD(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToSfermionsQCD(s) * ILee(beta,zmin)       &
        & + dgauss(EpEmToSfermionsQCDpISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToSfermionsQCD(s)

  End If

!----------------------------
! 0.38939e12 gives fb
!----------------------------
  sigma = 0.38939e12_dp * pi * sigma  

  Iname = Iname - 1

 End Subroutine EpEmToSquarksZG


 Real(dp) Function EpEmToSfermionsTree(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! written by Werner Porod, 27.8.99
 ! change to f90 standard, 1.10.2000
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Real(dp) :: kappa, kappa3d2, sumI, sqrts, gin(3), yukin(3), sinW2
  Real(dp) :: ae, ve, Le, Re, gg, fgZ, gZ, ZZ
  Complex(dp) :: coup, Rsf(2,2)

  kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
  kappa3d2 = kappa**1.5_dp 
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  sinW2 = gin(1)**2 / (gin(1)**2 + gin(2)**2)

  Rsf = RSf_in(1:2,1:2)
  Call CoupSfermionZ(ind_1, ind_2, gin(2), sinW2, e_in, T3_in, RSf, coup)
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ae = 2._dp * (Re - Le)
  ve = -2._dp * (Re + Le)

  If ((ind_1.Eq.ind_2).And.(e_in.Ne.0._dp) ) Then
   gg = e_in**2 * gin(2)**4 * sinW2**2 * (1._dp - P_m*P_p)
   fgZ = - 0.5_dp * ( ve * (1._dp - P_m*P_p) - ae * (P_m - P_p) )
   gZ = e_in * gin(2)**2 * sinW2 * Conjg( coup ) * fgZ
  Else
   gg = 0._dp
   fgZ = 0._dp
   gZ = 0._dp
  End If

  ZZ = 6.25e-2_dp *Abs(coup)**2 * ( (ve**2 + ae**2) * (1._dp - P_m*P_p) &
     &                            - 2._dp * ae * ve * (P_m - P_p) )

  sumI = gg + ( gZ * (s-mZ2) + ZZ * s ) * s / ( (s-mZ2)**2 + gmZ2 )

  EpEmToSfermionsTree = oo16pi2 * kappa3d2 * sumI / s**4

 End Function EpEmToSfermionsTree


 Real(dp) Function EpEmToSfermionsISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! for the calculation of ISR corrections
 ! written by Werner Porod, 26.9.99
 ! 1.10.2000: porting to f90 standard
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  EpEmToSfermionsISR = Lee(beta,x) * EpEmToSfermionsTree(s)

 End Function EpEmToSfermionsISR


 Real(dp) Function EpEmToSfermionsQCD(s)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! written by Werner Porod, 27.9.99
 ! last change: 27.9.99
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: s

  Real(dp) :: gin(3), yukin(3), kappa, kappa3d2, sumI, erg, sinW2
  Real(dp) :: alpha3, SqrtS, m12, m22, del_eq, del_a, factor, aij
  Real(dp) :: vq, aq, Le, Re, ae, ve 
  Complex(dp) :: coup

  kappa = (s-mSfer2(1)-mSfer2(2))**2 - 4._dp * mSfer2(1) * mSfer2(2)
  kappa3d2 = kappa**1.5_dp 
  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  sinW2 = gin(1)**2 / (gin(1)**2 + gin(2)**2)
  alpha3 = oo4pi * gin(3)**2

  vq = 1._dp - 4._dp * Qq * sinw2
  aq = 1._dp
  a(1,1) = 2._dp * Abs(costh)**2 - 4._dp * Qq * sinw2
  a(2,2) = 2._dp * Abs(sinth)**2 - 4._dp * Qq * sinw2

  Call CoupSfermionZ(ind_1, ind_2, gin(2), sinW2, e_in, T3_in  &
                 & , RSf_in(1:n_sf,1:n_sf), coup)
  Call CoupFermionZ(-0.5_dp, -1._dp, gin(2), sinW2, Le, Re)
  ae = 2._dp * (Re - Le)
  ve = -2._dp * (Re + Le)

  If (ind_1.Eq.ind_2) Then
   fgg = (1._dp - P_m*P_p)
   gg = e_in**2 * gin(2)**4 * sinW2**2 * (1._dp - P_m*P_p)
   fgZ = - 0.5_dp * ( ve * (1._dp - P_m*P_p) - ae * (P_m - P_p) )
   gZ = e_in * gin(2)**2 * sinW2 * Conjg(coup) * fgZ
  Else
   fgg = 0._dp
   gg = 0._dp
   fgZ = 0._dp
   gZ = 0._dp
  End If

  fZZ = 6.25e-2_dp * ( (ve**2 + ae**2) * (1._dp - P_m*P_p) &
      &              - 2._dp * ae * ve * (P_m - P_p) )
  ZZ = Abs(coup)**2 * fZZ

  !------------------
  ! QCD correction      
  !------------------
  m12 = mSfer2(1)
  m22 = mSfer2(2)
  Call eeSfermion_fQCD(s,m12,m22,erg)

!  Call delta_gluino(s, del_eq, del_a, aij)
  del_eq = 0._dp
  del_a = 0._dp
  aij = 0._dp
  factor = fo3pi*alpha3 
  sumI = gg + ( gZ * (s-mZ2) + ZZ * s ) * s / ( (s-mZ2)**2 + gmZ2 )

  EpEmToSfermionsQCD = oo16pi2 * kappa3d2  &
     &      * ( sumI * (1._dp + factor *erg )  &   ! gluon
     &        + factor * ( fgg * Qq * del_eq  &   ! gluino
     &                   + ( fgz * (Qq * del_a + aij * del_eq) * (s-mZ2) &
     &                     + 0.5_dp * fzz * aij * del_a * s )             &
     &          * 0.125_dp / ( (s-mZ2)**2 + gmZ2 )  &
     &                    )                        &
     &        ) / s**4

 End Function EpEmToSfermionsQCD


#ifdef BEAMSTRAHLUNG
 Real(dp) Function EpEmToSfermionsBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! for the calculation of ISR and Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 

  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToSfermionsBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToSfermionsTree(s) ! EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToSfermionsBeam = BeamFactor * ISRfactor2 * sigmaT!            &
!                     &  * kappa3d2 * alpha2 * sumI / s**4

 End Function EpEmToSfermionsBeam
#endif

 Real(dp) Function EpEmToSfermionsQCDpISR(x)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! includes QCD + ISR corrections
 ! written by Werner Porod, 29.12.99
 ! last change: 29.12.99
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  s = smax * (1._dp - Exp(x) )

  EpEmtoSfermionsQCDpISR = Lee(beta,x) * EpEmToSfermionsQCD(s)

 End Function EpEmToSfermionsQCDpISR


#ifdef BEAMSTRAHLUNG
 Real(dp) Function EpEmToSfermionsQCDpBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeSfermions and related functions.
 ! includes QCD + ISR + Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 
  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToSfermionsQCDpBeam = 0._dp
   Return
  End If


  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToSfermionsQCD(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmtoSfermionsQCDpBeam =  BeamFactor * ISRfactor2 * sigmaT
!  EpEmtoSfermionsQCDpBeam =  BeamFactor * ISRfactor2 * EpEmToSfermionsQCDpBeam 

 End Function EpEmToSfermionsQCDpBeam
#endif


 Subroutine eeSfermion_fqcd(skin,m12,m22,erg)
 !-------------------------------------------------------------
 ! Calculates total alpha_s corrections to e+e- -> squark_1 squark_2
 ! including all two and three jets events as given by H. Eberl
 ! written by Porod Werner, 22.10.1995
 ! 1.10.2000: change to f90 standard
 !-------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: skin, m12, m22
  Real(dp), Intent(out) :: erg
  Real(dp) :: mui, muj, mui2, muj2, sq_laij, laij, adsq,       &
   & la0, la1, la2, ln0, ln1, ln2, la12, la22, aux1, auxp, auxm, laij2 

  mui2 = m12/skin
  muj2 = m22/skin
  mui = Sqrt(mui2)
  muj = Sqrt(muj2)
  laij = (1._dp - mui2 -muj2)**2 - 4._dp * mui2 * muj2
  sq_laij = Sqrt(laij)
  laij2 = laij**2
  aux1 = 1._dp - mui2 - muj2 
  adsq = aux1 / sq_laij
  auxp = aux1 + sq_laij
  auxm = aux1 - sq_laij

  la0 = auxp/(2._dp*mui*muj)
  la1 = (1._dp - mui2 + muj2 - sq_laij)/(2._dp*muj)
  la12 = la1**2
  la2 = (1._dp + mui2 - muj2 - sq_laij)/(2._dp*mui)
  la22 = la2**2
  ln0 = Log(la0)
  ln1 = Log(la1)
  ln2 = Log(la2)

  ! \Delta^V_{ij}+\Delta^R_{ij} together:
   erg = Log(mui*muj) + 2._dp                                          &
    &  + adsq * ( Log(auxp/auxm) + 2._dp*pi2/3._dp + Li2(la12)            &
    &           + Li2(la22) - 2._dp*Log(sq_laij)*ln0 + 2._dp*ln0**2      &
    &           + 2._dp*Li2(1._dp-la0**2)-Li2(1._dp-la12)-Li2(1._dp-la22) )   &
    &  + (adsq*ln0 -1._dp)*Log(laij2/(mui2*muj2))                       &
    &  + (1._dp + mui2 + muj2) / laij                                   &
    &  + 4._dp/(sq_laij**3)*(mui2*ln2 + muj2*ln1 +  mui2*muj2*ln0)       &
    &  + (1._dp+2*muj2)/sq_laij*ln2                                      &
    &  + (1._dp+2*mui2)/sq_laij*ln1                                      &
    &  + (2._dp+mui2+muj2)/sq_laij*ln0

 End Subroutine eeSfermion_fqcd



 Real(dp) Function EpEmToTauLeptonsBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeTauLeptons and related functions.
 ! for the calculation of ISR and Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 

  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToTauLeptonsBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToTauLeptonsTree(s) ! EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToTauLeptonsBeam = BeamFactor * ISRfactor2 * sigmaT

 End Function EpEmToTauLeptonsBeam

 
 Subroutine EpEmToTauLeptons(Pm, Pp, s, ISR,Beam , sigma, design)
 !--------------------------------------------------------------------------
 ! calculates the cross section for top-quarj production in e+e- annihilation
 ! written by Werner Porod, 03.10.01
 ! 15.10.01: I take the formulas for top production, and replace the
 !           couplings to Z and photon
 !--------------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Pm, Pp, s
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR, Beam

  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, region(8), Ebeam, sinW2
  Real(dp) :: Le, Re

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToTauLeptons'

  !-------------------------------
  ! initialisation
  !-------------------------------
  sigma = 0._dp
  !-------------
  ! kinematics
  !-------------
  If ( (4._dp * mf_l2(3)) .Ge.s) Then
   Iname = Iname - 1
   Return
  End If

  !--------------------------------------------------------------------
  ! couplings
  ! setting g to 1 because the running coupling is included later
  !--------------------------------------------------------------------
  sinW2 = 1._dp - mW2 / mZ2
  Call CoupFermionZ(-0.5_dp, -1._dp, 1._dp, sinW2, Le, Re)

  gg = (1._dp - Pm * Pp) * sinW2**2 / (2._dp * Pi )
  gZ = - ( Le * (1._dp - Pm) * (1._dp + Pp)                        &
    &    + Re * (1._dp + Pm) * (1._dp - Pp) ) * (Le + Re)   &
    & * sinW2 * oo4pi
  ZZ = oo4pi * ( Le**2 * (1._dp - Pm) * (1._dp + Pp)     &
     &         + Re**2 * (1._dp + Pm) * (1._dp - Pp) )
  ZZi(2) = ZZi(1) * Le * Re * mf_l2(3)
  ZZi(1) = 0.5_dp * ZZi(1) * (Le**2 + Re**2 )

  vQ = 0.5_dp * (Re - Le)
  aQ = - 0.5_dp * (Re + Le)

  If (Beam) Then
   smin = 4._dp * mf_l2(3)
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = 2._dp * mf_l(3) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToTauLeptonsBeam,init,ncall,itmx,nprn,1.e-3_dp,erg)
   sigma = erg * zmax**2

  Else If (ISR) Then
   smin = 4._dp * mf_l2(3)
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToTauLeptonsTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToTauLeptonsTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToTauLeptonsISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToTauLeptonsTree(s)
  End If

  !--------------------------------------------
  ! 0.38939e12 gives fb , times colour factor
  !--------------------------------------------
  sigma = 3._dp * 0.38939e12_dp * sigma

  Iname = Iname - 1

 End Subroutine EpEmToTauLeptons


  Real(dp) Function EpEmToTauLeptonsISR(x)

  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToTauLeptonsISR'

  s = smax * (1._dp - Exp(x))

  EpEmToTauLeptonsISR = Lee(beta,x) * EpEmToTauLeptonsTree(s)

  Iname = Iname - 1

 End Function EpEmToTauLeptonsISR


 Real(dp) Function EpEmToTauLeptonsTree(s)
 Implicit None

  Real(dp), Intent(in) :: s

  Real(dp) :: SqrtS, EiEj, p_char, sigma_ij(3), invPropZ, sp_char, s2    &
      &     , gin(3), yukin(3), p2_char

  p2_char = 0.25_dp * s -  mf_l2(3)
  SqrtS = Sqrt( s )
  sp_char = Sqrt( p2_char ) * SqrtS
  p_char = sp_char / s
  
  sigma_ij = 0._dp

  invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )

   s2 = s**2
   EiEj = p2_char + mf_l2(3)
  !------------------------
  ! running couplings
  !------------------------
  Call runningCouplings(Sqrts,gin,yukin)

  !---------------
  ! cross section
  !---------------
  sigma_ij(1) = gg * p_char * (EiEj + p2_char / 3._dp + mf_l2(3)) / s2
  sigma_ij(2) = gZ * (s-mZ2) * invPropZ * p_char &
            & * (EiEj + p2_char / 3._dp + mf_l2(3)) / s

  sigma_ij(3) = ZZ * invPropZ * p_char                             &
         &  * ( 0.5_dp * (vQ**2 + aQ**2 ) * (EiEj + p2_char/3._dp) &
         &    + (-vQ**2 + aQ**2 ) *  mf_l2(3) )

  EpEmToTauLeptonsTree = gin(2)**4 * Sum( sigma_ij )

 End Function EpEmToTauLeptonsTree


 Real(dp) Function EpEmToTopQuarksBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeTopQuarks and related functions.
 ! for the calculation of ISR and Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 

  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToTopQuarksBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToTopQuarksTree(s) ! EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToTopQuarksBeam = BeamFactor * ISRfactor2 * sigmaT

 End Function EpEmToTopQuarksBeam


 Subroutine EpEmToTopQuarks(Pm, Pp, s, ISR,Beam , sigma, design)
 !--------------------------------------------------------------------------
 ! calculates the cross section for top-quarj production in e+e- annihilation
 ! written by Werner Porod, 03.10.01
 ! 03.10.01: I take the formulas for chargino production, and replace the
 !           couplings to Z and photon, at the moment only tree-level,
 !           will be improved soon 
 !--------------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Pm, Pp, s
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR, Beam

  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, region(8), Ebeam, sinW2
  Real(dp) :: Le, Re, Lt, Rt, eT

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToTopQuarks'

  !-------------------------------
  ! initialisation
  !-------------------------------
  sigma = 0._dp
  !-------------
  ! kinematics
  !-------------
  If ( (4._dp * mf_u2(3)) .Ge.s) Then
   Iname = Iname - 1
   Return
  End If

  eT = 2._dp / 3._dp
  !--------------------------------------------------------------------
  ! couplings
  ! setting g to 1 because the running coupling is included later
  !--------------------------------------------------------------------
  sinW2 = 1._dp - mW2 / mZ2
  Call CoupFermionZ(-0.5_dp, -1._dp, 1._dp, sinW2, Le, Re)
  Call CoupFermionZ(0.5_dp, eT, 1._dp, sinW2, Lt, Rt)

  gg = eT**2 * (1._dp - Pm * Pp) * sinW2**2 / (2._dp * Pi )
  gZ = eT * ( Le * (1._dp - Pm) * (1._dp + Pp)                 &
    &       + Re * (1._dp + Pm) * (1._dp - Pp) ) * (Lt + Rt)   &
    & * sinW2 * oo4pi
  ZZ = oo4pi * ( Le**2 * (1._dp - Pm) * (1._dp + Pp)     &
     &         + Re**2 * (1._dp + Pm) * (1._dp - Pp) )
  ZZi(2) = ZZi(1) * Lt * Rt * mf_u2(3)
  ZZi(1) = 0.5_dp * ZZi(1) * (Lt**2 + Rt**2 )

  vQ = 0.5_dp * (Rt - Lt)
  aQ = - 0.5_dp * (Rt + Lt)
  QCDfact = 0.5_dp * Pi - 0.75_dp / Pi

  If (Beam) Then
   smin = 4._dp * mf_u2(3)
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = 2._dp * mf_u(3) / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToTopQuarksBeam,init,ncall,itmx,nprn,1.e-3_dp,erg)
   sigma = erg * zmax**2

  Else If (ISR) Then
   smin = 4._dp * mf_u2(3)
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToTopQuarksTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToTopQuarksTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToTopQuarksISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToTopQuarksTree(s)
  End If

  !--------------------------------------------
  ! 0.38939e12 gives fb
  !--------------------------------------------
  sigma = 0.38939e12_dp * sigma

  Iname = Iname - 1

 End Subroutine EpEmToTopQuarks


  Real(dp) Function EpEmToTopQuarksISR(x)

  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToTopQuarksISR'

  s = smax * (1._dp - Exp(x))

  EpEmToTopQuarksISR = Lee(beta,x) * EpEmToTopQuarksTree(s)

  Iname = Iname - 1

 End Function EpEmToTopQuarksISR


 Real(dp) Function EpEmToTopQuarksTree(s)
 Implicit None

  Real(dp), Intent(in) :: s

  Real(dp) :: SqrtS, EiEj, p_char, sigma_ij(3), invPropZ, sp_char, s2    &
      &     , gin(3), yukin(3), p2_char, as, KV, KA, betaT, betaT2

  p2_char = 0.25_dp * s -  mf_u2(3)
  SqrtS = Sqrt( s )
  sp_char = Sqrt( p2_char ) * SqrtS
  p_char = sp_char / s
  
  sigma_ij = 0._dp

  invPropZ = 1._dp / ( (s-mZ2)**2 + gmZ2 )

   s2 = s**2
   EiEj = p2_char + mf_u2(3)

  !------------------------
  ! running couplings
  !------------------------
  Call runningCouplings(Sqrts,gin,yukin)
  !-----------------
  ! QCD corrections
  !-----------------
  as = oo4pi * gin(3)**2
  betaT = Sqrt( 1._dp - 4._dp * mf_u2(3)/s)
  betaT2 = betaT**2
  KV = 1._dp + 4._dp * as * (0.5_dp * Pi / betaT                   &
     &                      - 0.25_dp * (3._dp + betaT) * QCDfact  &
     &                      ) / 3._dp
  KA = 1._dp                                                                  &
     & + 4._dp * as * (0.5_dp * Pi / betaT                                    &
     &                - (1.9_dp - 4.4_dp * betaT + 3.5_dp * betaT2) * QCDfact &
     &                ) / 3._dp

  !----------------------------------
  ! cross section including QCD
  !----------------------------------
  sigma_ij(1) = gg * p_char * (EiEj + p2_char / 3._dp + mf_u2(3)) / s2
  sigma_ij(1) = sigma_ij(1) * KV
  sigma_ij(2) = gZ * (s-mZ2) * invPropZ * p_char &
            & * (EiEj + p2_char / 3._dp + mf_u2(3)) / s
  sigma_ij(2) = sigma_ij(2) * KV

  sigma_ij(3) = ZZ * invPropZ * p_char                                       &
         &  * ( 0.5_dp * (vQ**2 * KV + aQ**2 * KA ) * (EiEj + p2_char/3._dp) &
         &    + (-vQ**2 * KV + aQ**2 * KA ) *  mf_u2(3) )

  EpEmToTopQuarksTree = gin(2)**4 * Sum( sigma_ij )

 End Function EpEmToTopQuarksTree


 Real(dp) Function EpEmToWpWmBeam(z,wgt)
 !-----------------------------------------------------------------------
 ! auxiliary function for the calculation of sfermion production in
 ! e+e- annihilation. Is called by eeWpWm and related functions.
 ! for the calculation of ISR and Beam corrections
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: z(:), wgt

  Integer :: i_run, i_near
  Real(dp) :: s, x1, x2, x3, x4
  Real(dp), Save :: BeamFactor, ISRfactor2, sigmaT, dsigmaT, s_save

  !------------------------------------------------------
  ! Parameter for Beamstrahlung
  !------------------------------------------------------
  If (z(1).Ge.radprob) Then
   x1 = 1._dp 
   BeamFactor = 1._dp 
  Else
   x1 = 1._dp - (1._dp - z(1) / radprob)**3  
   BeamFactor = BeamElectronDistribution( z(1) ) 
  End If
  If (z(3).Ge.radprob) Then
   x3 = 1._dp 
  Else
   x3 = 1._dp - (1._dp - z(3)  / radprob)**3  
   BeamFactor = BeamFactor * BeamElectronDistribution( z(3) ) 
  End If
  !------------------------------------------------------
  ! Parameter for ISR
  !------------------------------------------------------
  x2 = 1._dp - (zmax*z(2))**(2._dp / EtaBeam)
  ISRfactor2 = ISRElectronDistribution(x2)
  x4 = 1._dp - (zmax*z(4))**(2._dp / EtaBeam)
  ISRfactor2 = ISRfactor2 * ISRElectronDistribution(x4)

  s = 0.25_dp * smax * (x1 * x2 + x3 * x4)**2 

  If (s.Le.smin) Then ! kinematically forbidden
   EpEmToWpWmBeam = 0._dp
   Return
  End If

  If (CalculateSpline) Then ! in the first call, calculate the 200 points
   DeltaS = (smax - smin) / 2.e2_dp
   Spline(1,1) = smin
   Spline(1,2) = 0._dp  ! kinematically
   s_save = s
   Do i_run = 2,201
    s = smin + (i_run-1) * DeltaS
    Spline(i_run,1) = s
    Spline(i_run,2) = EpEmToWpWmTree(s) ! EpEmToEsneutrinosTree(s)
   End Do ! i_run
   s = s_save
   CalculateSpline = .False.
  End If !CalculateSpline 

  i_near = Int( (s-smin)/DeltaS )
  If (i_near.Le.2) Then
   DeltaSplineX = Spline(1:5,1)
   DeltaSplineY = Spline(1:5,2)
  Else If (i_near.Ge.199) Then
   DeltaSplineX = Spline(197:201,1)
   DeltaSplineY = Spline(197:201,2)
  Else
   DeltaSplineX = Spline(i_near-2:i_near+2,1)
   DeltaSplineY = Spline(i_near-2:i_near+2,2)
  End If

  Call polint(DeltaSplineX,DeltaSplineY,s,sigmaT,dsigmaT) 

  EpEmToWpWmBeam = BeamFactor * ISRfactor2 * sigmaT

 End Function EpEmToWpWmBeam


 Subroutine EpEmToWpWm(Pm, Pp, s, ISR,Beam , sigma, design)
 !--------------------------------------------------------------------------
 ! calculates the cross section for W+W- production in e+e- annihilation
 ! The formulas of A.Aeppli et al.,DESY 92-123A, page 151ff are used
 ! written by Werner Porod, 04.10.01
 !--------------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: Pm, Pp, s
  Real(dp), Intent(out) :: sigma
  Logical, Intent(in) :: ISR, Beam

  Character (Len=*), Optional, Intent(in) :: Design

  Integer :: init, itmx, ncall, nprn
  Real(dp) :: erg, region(8), Ebeam, sinW2

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToWpWm'

  !-------------------------------
  ! initialisation
  !-------------------------------
  sigma = 0._dp
  !-------------
  ! kinematics
  !-------------
  If ( (4._dp * mW2) .Ge.s) Then
   Iname = Iname - 1
   Return
  End If

  RmZmW2 = 0.25_dp * mZ2 / mW2
  RmZmW4 = RmZmW2**2

  sinW2 = 1._dp - mW2 / mZ2
  FacWW(1) = 0.25_dp * (1._dp-Pm) * (1._dp+Pp) 
  FacWW(2) = 0.5_dp * sinW2 * (1._dp-Pm) * (1._dp+Pp)
  FacWW(3) = sinW2**2 * (1._dp-Pm*Pp)
  FacWW = 0.25_dp * oo16Pi * FacWW

  If (Beam) Then
   smin = 4._dp * mW2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   EtaBeam = - 6._dp * Log(1._dp - beta / 6._dp)
   ISRfactor = 0.5_dp * (1._dp + 0.5_dp * EtaBeam )                        &
           & * Exp(- 0.125_dp * (EtaBeam + (pi2/6._dp - 1._dp) *EtaBeam**2) )

   Ebeam = 0.5_dp * Sqrt(s)
   If (Present(Design)) Then
    Call BeamStrahlungInitzialization(Design,Ebeam)
   Else
    If (ErrorLevel.Ge.0) Then
     Write(ErrCan,*) 'Warning from subroutine ',NameOfUnit(Iname)
     Write(ErrCan,*) 'Machine design for Beamstrahlung not defined!'
     Write(ErrCan,*) 'Using therefore the design: TESLA500.'
     If (ErrorLevel.Eq.2) Call TerminateProgram
    End If
    Call BeamStrahlungInitzialization('TESLA500',Ebeam)
   End If
!----------------------------------------------------------------------
! minimal energy one beam needs to produce particles if the other beam
! still has its full energy, improves calculation near threshold
!----------------------------------------------------------------------
   ymin = 2._dp * mW / Ebeam - 1._dp
   If (ymin.Lt.0._dp) Then
    ymin = 0._dp 
    zmax = 1._dp
   Else
    zmax = (1._dp - ymin)**(0.5_dp * EtaBeam)
   End If
   region(1:4) = 0._dp 
   region(5:8) = 1._dp
   init = 0
   ncall = 15000
   itmx = 10
   nprn = -1
   CalculateSpline = .True.
   Call Vegas1(region,EpEmToWpWmBeam,init,ncall,itmx,nprn,1.e-3_dp,erg)
   sigma = erg * zmax**2

  Else If (ISR) Then
   smin = 4._dp * mW2
   smax = s
   beta = 2._dp * Alpha * ( Log(s/mf_l2(1)) - 1._dp ) / Pi
   zmax = Log(1._dp - smin / s)
   If (zmax.Lt.-8._dp) Then
    sigma = EpEmToWpWmTree(s) * ILee(beta,zmax)
   Else
    zmin = -10._dp 
    sigma =  EpEmToWpWmTree(s) * ILee(beta,zmin) &
        & + dgauss(EpEmToWpWmISR,zmin,zmax,1.e-3_dp)
   End If

  Else
   sigma = EpEmToWpWmTree(s)
  End If

  !---------------------
  ! 0.38939e12 gives fb 
  !---------------------
  sigma = 0.38939e12_dp * sigma

  Iname = Iname - 1

 End Subroutine EpEmToWpWm


  Real(dp) Function EpEmToWpWmISR(x)

  Implicit None

  Real(dp), Intent(in) :: x

  Real(dp) :: s

  Iname = Iname + 1
  NameOfUnit(Iname) = 'EpEmToWpWmISR'

  s = smax * (1._dp - Exp(x))

  EpEmToWpWmISR = Lee(beta,x) * EpEmToWpWmTree(s)

  Iname = Iname - 1

 End Function EpEmToWpWmISR


 Real(dp) Function EpEmToWpWmTree(s)
 Implicit None

  Real(dp), Intent(in) :: s
  Real(dp) :: w, w2, betaT, lnwbeta, z, oo3z2, SigmaII, SigmaIQ, SigmaQQ
  Real(dp) :: SqrtS, gin(3), yukin(3)

  w = 4._dp * mW2 / s
  w2 = w**2
  betaT = Sqrt(1._dp - w)
  lnwbeta = Log(w / (1._dp+betaT)**2 )
  z = mZ2 / s
  z = 1._dp - z
  oo3z2 = 1._dp / (3._dp * z**2)

  SigmaII = ( 2._dp * betaT * (-15._dp - RmZmW2 * (4._dp-10._dp*w + 3._dp*w2) &
          &                   + 4._dp * RmZmW4 * (1._dp+5._dp*w - 3._dp*w2) ) &
          & - 3._dp * z * ( 8._dp + 4._dp * w + w2                            &
          &               - 4._dp * (1._dp-z) * (2._dp - w) ) * lnwbeta       &
          & ) * oo3z2
  SigmaIQ = RmZmW2 * oo3z2                                                   &
        & * ( 2._dp * betaT * (4._dp +20._dp*w + 3._dp*w2)                   &
        &                   * (1._dp - RmZmW2 * (2._dp - w) )                &
        &   + 3._dp * z * w2  * (8._dp +  w) * lnwbeta   )
  SigmaQQ = RmZmW4 * oo3z2 * 2._dp * betaT * (1._dp -  w)   &
        &   * (4._dp +20._dp*w + 3._dp*w2)

  !------------------------
  ! running couplings
  !------------------------
  Sqrts = Sqrt(s)
  Call runningCouplings(Sqrts,gin,yukin)
  !-----------------
  ! QCD corrections
  !-----------------
  EpEmToWpWmTree = gin(2)**4 * ( FacWW(1) * SigmaII   &
                 &             + FacWW(2) * SigmaIQ   &
                 &             + FacWW(3) * SigmaQQ   ) / s

 End Function EpEmToWpWmTree

 Real(dp) Function Lee(beta,x)
 !-----------------------------------------------------------------------
 ! calculates the energy distribution for the ISR corrections
 ! for particle production at e+ e- colliders. beta is defined as
 ! 2 alpha_em ( Log[s/m^2_e] - 1 ) / Pi. (1-Exp[x]) is the fraction of energy
 ! in the c.m.s. after the photon(s) has (have) been emitted. s is the total 
 ! c.m.s.-energy squared, and m_e is the electron mass.
 ! Formulas are given by M.~Drees and K.~Hikasa, Phys.~Lett.~B252, 127 (1990)
 ! written by Werner Porod, 27.8.99
 ! 26.9.2000: changing to Fortran 90 form
 !-----------------------------------------------------------------------
  Implicit None

  Real(dp) , Intent(in) :: beta,x

  Lee = KFactorLee * beta &
  & * ( Exp(beta*x) * (1._dp + 0.75_dp * beta) - Exp(x) + 0.5_dp*Exp(2._dp*x))

  End Function Lee


 Real(dp) Function ILee(beta,x)
 !-----------------------------------------------------------------------
 ! calculates the Integral of Lee(beta,x) in the range -Infinity <= x1 <= x.
 ! beta is defined as 2 alpha_em ( Log[s/m^2_e] - 1 ) / Pi. s is the total
 ! c.m.s. energy squared, and m_e is the electron mass.
 ! written by Werner Porod, 27.8.99
 ! 26.9.2000: changing to Fortran 90 form
 !-----------------------------------------------------------------------
 Implicit None
 
  Real(dp) , Intent(in) :: beta,x

  ILee = KFactorLee * ( Exp(beta*x) * (1._dp + 0.75_dp * beta)     &
   &                  - beta * (Exp(x) - 0.25_dp * Exp(2._dp*x) ) )

 End Function ILee


#ifdef BEAMSTRAHLUNG
 Real(dp) Function BeamElectronDistribution(x)
 !-----------------------------------------------------------------
 ! calculates the electron distribution in case of beamstrahlung.
 ! Is based on the program pandora by M.Peskin and his
 ! paper 'Consistent Yokoya-Chen Approximation to Beamstrahlung'
 ! LCC-0010
 ! Input
 !  x1 .... energy fraction after beamstrahlung
 ! written by Werner Porod, 15.01.01
 !-----------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Real(dp) :: Beam, x1, mx1, y, H, y3, y6, y9, z

  !---------------------------
  ! delta peak
  !---------------------------
  If (x.Ge.radprob) Then
   BeamElectronDistribution = 1._dp
   Return
  End If

  x1 = 1._dp - (1._dp - x / radprob)**3
  If (x1.Eq.0._dp) x1 = 1.e-14_dp  ! to avoid division by 0
  !--------------------------
  ! beamstrahlung part
  !--------------------------
  mx1 = (1._dp - x1)**(1._dp / 3._dp) 
  If (mx1.Eq.0._dp) mx1 = 1.e-14_dp  ! to avoid division by 0
  y = KappaBeam * (1._dp - x1) / x1
  Beam = Exp(-Ngeff-y) / (x1 * mx1 )
  y = Ngeff * y**(1._dp / 3._dp)

  ! computes h(x) with a relative error < 0.001 everywhere
  If (y .Lt. 8.89_dp) Then
   y3 = y**3
   y6 = y3**2
   y9 = y3*y6
   H = (y/Gamma_1o3) * (1._dp + 0.125_dp*y3  + y6/2240._dp + y9/3763200._dp) &
   & + (y**2/Gamma_2o3) * (0.5_dp + 0.0125_dp * y3 + y6 / 44800._dp          &
   &                     + y9 / 118272000._dp )                             &
   & + y3 * (1._dp/6._dp + y3 / 720._dp + y6 / 725760._dp                   &
   &        + y9 / 28740096000._dp)
  Else
   z = (y/3._dp)**0.75_dp
   H = Sqrt(3._dp * oo8pi * z) * Exp(4._dp * z)                    &
   & * (1._dp - 35._dp / (288._dp * z) - 1295._dp / (16588._dp * z**2) )
  End If
  BeamElectronDistribution = 3._dp * Beam * H / radprob

 End Function BeamElectronDistribution
#endif


#ifdef BEAMSTRAHLUNG
 Subroutine BeamStrahlungInitzialization(Design,E)
 !------------------------------------------------------------------
 ! initalizes the variables for the different colliders as needed
 ! for the calculation of Beamstrahlung. The data are taken from 
 ! the program pandora 2.1. by M. Peskin
 ! written by Werner Porod, 26.12.2000
 !------------------------------------------------------------------
 Implicit None
  Character (Len=*), Intent(in) :: Design
  Real(dp), Intent(in) :: E
  Real(dp) :: N, sigmax, sigmay, sigmaz, betax, betay, Re, GammaFactor, &
                 & Dx, Dy, Hx, Hy, NuCl, NuGamma, Ncl, Ngamma, sigmaxB,      &
                 & sigmayB 


  Select Case(Design)
   Case ('NLC500A')        ! 500 GeV JLC/NLC A   2/24/00
       N = 0.7e10_dp
       sigmax = 310.0e-9_dp
       sigmay = 4.0e-9_dp
       sigmaz = 90.0e-6_dp
       betax = 12.0e-3_dp
       betay = 0.12e-3_dp
   Case ('NLC500B')        !  500 GeV JLC/NLC B   2/24/00 
       N = 0.82e10_dp
       sigmax = 330.0e-9_dp
       sigmay = 4.6e-9_dp
       sigmaz = 120.0e-6_dp
       betax = 12.0e-3_dp
       betay = 0.12e-3_dp
   Case ('NLC500C')        !  500 GeV JLC/NLC C   2/24/00 
       N = 1.0e10_dp
       sigmax = 365.0e-9_dp
       sigmay = 6.2e-9_dp
       sigmaz = 140.0e-6_dp
       betax = 13.0e-3_dp
       betay = 0.15e-3_dp
   Case ('NLC500H')        !  500 GeV JLC/NLC H   2/24/00 
       N = 0.75e10_dp
       sigmax = 245.0e-9_dp
       sigmay = 2.7e-9_dp
       sigmaz = 110.0e-6_dp
       betax =  8.0e-3_dp
       betay = 0.10e-3_dp
   Case ('NLC1000A')        !  1000 GeV JLC/NLC A  2/24/00
       N = 0.7e10_dp
       sigmax = 220.0e-9_dp
       sigmay = 2.8e-9_dp
       sigmaz = 90.0e-6_dp
       betax = 12.0e-3_dp
       betay = 0.12e-3_dp
   Case ('NLC1000B')        !  1000 GeV JLC/NLC B  2/24/00
       N = 0.82e10_dp
       sigmax = 235.0e-9_dp
       sigmay = 3.2e-9_dp
       sigmaz = 120.0e-6_dp
       betax = 12.0e-3_dp
       betay = 0.15e-3_dp
   Case ('NLC1000C')        !  1000 GeV JLC/NLC C  2/24/00
       N = 1.0e10_dp
       sigmax = 260.0e-9_dp
       sigmay = 4.4e-9_dp
       sigmaz = 140.0e-6_dp
       betax = 13.0e-3_dp
       betay = 0.15e-3_dp
   Case ('NLC1000H')        !  1000 GeV JLC/NLC H  2/24/00
       N = 0.75e10_dp
       sigmax = 200.0e-9_dp
       sigmay = 2.2e-9_dp
       sigmaz = 110.0e-6_dp
       betax = 10.0e-3_dp
       betay = 0.12e-3_dp
   Case ('NLC1500')        !  1500 GeV  modified RRuth Snowmass 96 
       N = 1.10e10_dp
       sigmax = 202.0e-9_dp
       sigmay = 5.1e-9_dp
       sigmaz = 150.0e-6_dp
       betax = 12.0e-3_dp
       betay = 0.20e-3_dp
   Case ('TESLA500')        !  500 GeV TESLA as of 6/28/00 
       N = 2.0e10_dp
       sigmax = 553.2e-9_dp
       sigmay = 5.0e-9_dp
       sigmaz = 30.0e-6_dp
       betax = 15.0e-3_dp
       betay = 0.40e-3_dp
    Case ('TESLA800')      !  800 GeV TESLA as of 6/28/00 
       N = 1.4e10
       sigmax = 391.e-9_dp
       sigmay = 2.0e-9_dp
       sigmaz = 30.0e-6_dp
       betax = 15.0e-3_dp
       betay = 0.30e-3_dp
   Case ('CLIC500')        !  500 GeV CLIC as of  7/11/00  
       N = 4.0e9
       sigmax = 202.0e-9_dp
       sigmay = 2.5e-9_dp
       sigmaz = 30.0e-6_dp
       betax = 10.0e-3_dp
       betay = 0.15e-3_dp
   Case ('CLIC1000')        !  1000 GeV CLIC as of  7/11/00  
       N = 4.0e9_dp
       sigmax = 115.0e-9_dp
       sigmay = 1.75e-9_dp
       sigmaz = 30.0e-6_dp
       betax = 10.0e-3_dp
       betay = 0.15e-3_dp
   Case ('CLIC3000')        !  3000 GeV CLIC as of  7/11/00  
       N = 4.0e9_dp
       sigmax = 43.0e-9_dp
       sigmay = 1.0e-9_dp
       sigmaz = 30.0e-6_dp
       betax = 8.0e-3_dp
       betay = 0.15e-3_dp
   Case ('CLIC5000')        !  5000 GeV CLIC as of  7/11/00  
       N = 4.0e9
       sigmax = 31.0e-9_dp
       sigmay = 0.78e-9_dp
       sigmaz = 25.0e-6_dp
       betax = 6.0e-3_dp
       betay = 0.10e-3_dp
   Case Default
    Write(ErrCan,*) 'The design you have chosen for the initialization of'
    Write(ErrCan,*) ' BeamStrahlung is not defined.'
    Call TerminateProgram
   End Select

   Re = 2.817940e-15_dp ! classical electron radius 
   GammaFactor = E / mf_l(1)

   Dx = 2 * N * Re * sigmaz / (GammaFactor * sigmax * (sigmax+sigmay) )
   Dy = 2 * N * Re * sigmaz / (GammaFactor * sigmay * (sigmax+sigmay) )

   Hx = 1 + Dx**0.25_dp * Dx**3 &
      &     * ( Log(1._dp+Sqrt(Dx))+2*Log(0.8_dp*betax/sigmaz) ) / (1+Dx**3)
   Hy = 1 + Dy**0.25_dp * Dy**3 &
      &     * ( Log(1._dp+Sqrt(Dy))+2*Log(0.8_dp*betay/sigmaz) ) / (1+Dy**3)

   sigmaxB = sigmax / Sqrt(Hx)
   sigmayB = sigmay / Hy**(1._dp / 3._dp)

   Upsilon = 5 * Re**2 * GammaFactor * N &
           & / (6 * Alpha * sigmaz * (sigmaxB +sigmayB) )

   KappaBeam = 2._dp / (3._dp * Upsilon) 

   NuCl = 2.5_dp * Alpha**2 * Upsilon / (Sqrt(3._dp) * Re * GammaFactor )
   NuGamma = NuCl / Sqrt(1._dp + Upsilon**(2._dp/3._dp) )

   Ncl = Sqrt(3._dp) * sigmaz * NuCl
   Ngamma = Sqrt(3._dp) * sigmaz * NuGamma

   Ngeff = NGamma * 0.5_dp ! = NGamma * Ngfraction
   radprob = 1._dp - Exp(-Ngeff)

 End Subroutine BeamStrahlungInitzialization
#endif

#ifdef BEAMSTRAHLUNG
 Real(dp) Function ISRElectronDistribution(x1)
 !-----------------------------------------------------------------
 ! calculates the electron distribution in case of ISR. 
 ! Is based on the program pandora by M.Peskin and his
 ! paper 'Consistent Yokoya-Chen Approximation to Bemastrahlung'
 ! LCC-0010
 ! Input
 !  x1 .... energy fraction after ISR
 ! written by Werner Porod, 15.01.01
 !-----------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x1

  Real(dp) :: ISR, x12
  !--------------------------
  ! ISR part
  !--------------------------
  x12 = x1**2
  ISR = (1._dp + x12)                                                  &
    & - 0.25_dp * EtaBeam * (0.5_dp *(1._dp + 3._dp *x12 ) * Log(x1)   &
    &                       + (1._dp - x1)**2 )
  ISRElectronDistribution = ISRfactor * ISR 

 End Function ISRElectronDistribution
#endif
 
End Module EplusEminusProduction

