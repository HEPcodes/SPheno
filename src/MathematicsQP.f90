Module MathematicsQP

! here we collect mathematical routines using multi precision software

! load modules
Use Control
! load modules

! interfaces
 Interface EigenSystemQP
  Module Procedure ComplexEigenSystem_DP, RealEigenSystem_DP   &
      & , ComplexEigenSystem_QP, RealEigenSystem_QP            &
      & , ComplexEigenSystem_DP1
 End Interface
! interfaces

! private variables
 Real(dp), Parameter, Private :: MinimalPrecision = 1.e-25_dp
 Private :: ComplexEigenSystem_DP, RealEigenSystem_DP   &
      & , ComplexEigenSystem_QP, RealEigenSystem_QP,Tred2a_QP &
      & , ComplexEigenSystem_DP1, Tqli_QP, Pythag_QP    &
      & , HTRIBK_QP, HTRIDI_QP
! private variables

 Real(qp), Private :: Zero=0._qp, One=1._qp, PointTwo=0.2_qp, PointFive=0.5_qp &
    &  , Hundred=1.e2_qp
 Complex(qp), Private :: IOne = (0._qp,1._qp)

Contains

 Real(qp) Function ArgQP(z)
 Implicit None

 Complex(qp), Intent(in) :: z
 Real(qp) :: x,y

  If (Abs(z).Eq.Zero) Then
   ArgQP = Zero
  Else
   x = Real(z)
   y = Aimag(z)

   ArgQP = Atan2(y,x)
  Endif

 End Function ArgQP

 Subroutine ComplexEigenSystem_DP(Matrix, EigenValues, EigenVectors, kont, test)
 !---------------------------------------------------------------------
 ! Subroutine for diagonalization of complex hermitian matrices, based on the
 ! Householder algorithm. Is a portation of  EISCH1 to F90
 ! Input:
 !  Matrix ..... n times n matrix
 ! Output
 !  EigenValues ..... n sorted EigenValues: |m_1| < |m_2| < .. < |m_n|
 !  EigenVectors .... n times n matrix with the eigenvectors
 ! written by Werner Porod, 10.11.2000
 ! 19.07.02: adapting to multi precision
 !---------------------------------------------------------------------
 Implicit None
  !-------
  ! input
  !-------
  Complex(dp), Intent(in) :: Matrix(:,:)
  !--------
  ! output
  !--------
  Integer, Intent(inout) :: kont
  Complex(dp), Intent(out) :: EigenVectors(:,:)
  Real(dp), Intent(out) :: EigenValues(:), test(:)

  !-----------------
  ! local variables
  !-----------------
  Integer :: i1,N1,N2,N3, i2, i3, i4, nrot
  Real(qp) :: AbsAi
  Real(qp), Allocatable :: AR(:,:),AI(:,:), WR(:), ZR(:,:),  WORK(:) &
          & , ZR_in(:,:), testR(:,:), Ar2(:,:), Ai2(:,:)
  Complex(qp), Allocatable ::  ctest(:,:), Rot(:,:)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'ComplexEigenSystem_DP'

  kont = 0

  N1 = Size(Matrix, Dim=1)
  N2 = Size(EigenValues)
  N3 = Size(EigenVectors, Dim=1)
  If ((N1.Ne.N2).Or.(N1.Ne.N3)) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) 'Dimensions to not match: ',N1,N2,N3
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -1001
   Call AddError(1001)
   Return
  End If

  Allocate(AR(N1,N1))
  Allocate(AI(N1,N1))

  AbsAi = Zero
  AR = Real( Matrix, qp)
  Ai = Aimag( Matrix )
  AbsAi = Sum( Abs( Ai ) )

  !--------------------------------------------------------------------------
  ! check first whether the matrix is really complex
  ! if not, I use the only real diagonalization because it is more accurate
  !--------------------------------------------------------------------------
  If (AbsAi .Eq. Zero) Then ! real matrix
   Allocate(WR(N1))
   Allocate(Work(N1))
   Allocate(ZR_in(N1,N1))
   Allocate(testR(N1,N1))
   Do i1=1,N1
    Do i2=1,N1
     ZR_in(i1,i2) = AR(i1,i2)
     AR(i1,i2) = Zero
    End Do
    wr(i1) = Zero
   End Do
   Call JacobiQP(ZR_in, n1, n1, wr, ar, nrot)

   Do n2=1,n1-1
    Do n3=n2+1,n1
     If (Abs(wr(n2)).Gt.Abs(wr(n3))) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      Do i1=1,n1
       work(1) = ar(i1,n2)
       ar(i1,n2) = ar(i1,n3)
       ar(i1,n3) = work(1)
      End Do
     End If
    End Do
   End Do

   Do i1=1,N1
    EigenValues(i1) = WR(i1)
    Do i2=1,n2
     EigenVectors(i1,i2) = AR(i2,i1) 
    End Do
   End Do
   ! now  a test

   ZR_in = Real( Matrix, qp )

   test = 0._dp
   Do i1=1,n1
    Do i2=1,n1
     testR(i1,i2) = Zero
     Do i3=1,n1
      Do i4=1,n1
       testR(i1,i2) = testR(i1,i2) &
          & + Ar(i3,i1)* ZR_in(i3,i4) *  Ar(i4,i2)
      End Do
     End Do
     If (i1.Eq.i2) Then
      If ( test(1).Lt.Abs( testR(i1,i2) ) ) test(1) = Abs( testR(i1,i2) )
     Else
      If ( test(2).Lt.Abs( testR(i1,i2) ) ) test(2) = Abs( testR(i1,i2) )
     End If 
    End Do
   End Do

   If (test(1).Gt.0._dp) Then
    If ( (test(2)/test(1)).Gt.MinimalPrecision) Then
     kont = -1002
     Call AddError(1002)
    End If
   End If

   Deallocate( testR )

  Else ! complex matrix

   Allocate(ZR(2*N1,2*N1))
   Allocate(ZR_in(2*N1,2*N1))
   Allocate(WR(2*N1))
   Allocate(Work(2*N1))
   Allocate(CTest(N1,N1))
   Allocate(Rot(N1,N1))
   Allocate(ar2(N1,N1))
   Allocate(ai2(N1,N1))

   Do i1=1,n1
    Do i2=1,n1
     Ar2(i1,i2) = Zero
     Ai2(i1,i2) = Zero
     Do i3=1,n1
      Ar2(i1,i2) = ar2(i1,i2) + ar(i3,i1)*ar(i3,i2) + ai(i3,i1)*ai(i3,i2)
      Ai2(i1,i2) = ai2(i1,i2) + ar(i3,i1)*ai(i3,i2) - ai(i3,i1)*ar(i3,i2)
     End Do
    End Do
   End Do

   ZR_in(1:N1,1:N1) = AR2
   ZR_in(N1+1:2*N1,N1+1:2*N1) = AR2
   ZR_in(N1+1:2*N1,1:N1) = AI2
   Do i1=1,n1
    Do i2=1,n1
     ZR_in(i1,N1+i2) = - AI2(i1,i2)
    End Do
   End Do

   Call JacobiQP(ZR_in, 2*n1, 2*n1, wr, zr, nrot)

   Do n2=1,2*n1-1
    Do n3=n2+1,2*n1
     If (wr(n2).Gt.wr(n3)) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      work = zr(:,n2)
      zr(:,n2) = zr(:,n3)
      zr(:,n3) = work
     End If
    End Do
   End Do

   Do i1=1,n1
    eigenvalues(i1) = Sqrt( wr(2*i1-1) )
    Do i2=1,n1
     eigenvectors(i1,i2) =  zr(i2,2*i1-1) - Ione * zr(n1+i2,2*i1-1)
     rot(i2,i1) =  zr(i2,2*i1-1) + IOne * zr(n1+i2,2*i1-1)
    End Do
   End Do

   test = 0._dp
   Do i1=1,n1
    Do i2=1,n1
     Ctest(i1,i2) = Zero
     Do i3=1,n1
      Do i4=1,n1
       Ctest(i1,i2) = Ctest(i1,i2) &
          & + Conjg( rot(i3,i1) )* (ar2(i3,i4)+Ione*ai2(i3,i4)) * rot(i4,i2)
      End Do
     End Do
     If (i1.Eq.i2) Then
      If ( test(1).Lt.Abs( Ctest(i1,i2) ) ) test(1) = Abs( Ctest(i1,i2) )
     Else
      If ( test(2).Lt.Abs( Ctest(i1,i2) ) ) test(2) = Abs( Ctest(i1,i2) )
     End If 
    End Do
   End Do

   If (test(1).Gt.0._dp) Then
    If ( (test(2)/test(1)).Gt.MinimalPrecision) Then
     kont = -1002
     Call AddError(1002)
    End If
   End If

   Deallocate(ZR, Ctest, Rot, ar2, ai2)

  End If ! decision whether real or complex matrix

  Deallocate(AR,AI,WR,Work, ZR_in)

  Iname = Iname - 1

 End Subroutine ComplexEigenSystem_DP

 Subroutine ComplexEigenSystem_QP(Matrix, EigenValues, EigenVectors, kont)
 !---------------------------------------------------------------------
 ! Subroutine for diagonalization of complex hermitian matrices, based on the
 ! Householder algorithm. Is a portation of  EISCH1 to F90
 ! Input:
 !  Matrix ..... n times n matrix
 ! Output
 !  EigenValues ..... n sorted EigenValues: |m_1| < |m_2| < .. < |m_n|
 !  EigenVectors .... n times n matrix with the eigenvectors
 ! written by Werner Porod, 10.11.2000
 ! 19.07.02: adapting to multi precision
 !---------------------------------------------------------------------
 Implicit None
  !-------
  ! input
  !-------
  Complex(qp), Intent(in) :: Matrix(:,:)
  !--------
  ! output
  !--------
  Complex(qp), Intent(out) :: EigenVectors(:,:)
  Real(qp), Intent(out) :: EigenValues(:)
  Integer, Intent(inout) :: kont

  !-----------------
  ! local variables
  !-----------------
  Integer :: i1,N1,N2,N3, i2, i3, i4, nrot
  Real(qp) :: AbsAi, AbsTest
  Real(qp), Allocatable :: AR(:,:),AI(:,:), WR(:), ZR(:,:),  WORK(:) &
          & , ZR_in(:,:), testR(:,:)
  Complex(qp), Allocatable ::  test(:,:)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'ComplexEigenSystem_QP'

  kont = 0

  N1 = Size(Matrix, Dim=1)
  N2 = Size(EigenValues)
  N3 = Size(EigenVectors, Dim=1)
  If ((N1.Ne.N2).Or.(N1.Ne.N3)) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) 'Dimensions to not match: ',N1,N2,N3
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -1003
   Call AddError(1003)
   Return
  End If

  Allocate(AR(N1,N1))
  Allocate(AI(N1,N1))
  Allocate(Test(N1,N1))


  AbsAi = Zero
  AR = Real( Matrix, qp)
  Ai = Aimag( Matrix )
  AbsAi = Sum( Abs( Ai ) )

  !--------------------------------------------------------------------------
  ! check first whether the matrix is really complex
  ! if not, I use the only real diagonalization because it is more accurate
  !--------------------------------------------------------------------------
  If (AbsAi .Eq. Zero) Then ! real matrix

   Allocate(WR(N1))
   Allocate(Work(N1))
   Allocate(ZR_in(N1,N1))
   Allocate(testR(N1,N1))
   ZR_in = AR
   Call JacobiQP(ZR_in, n1, n1, wr, ar, nrot)

   Do n2=1,n1-1
    Do n3=n2+1,n1
     If (Abs(wr(n2)).Gt.Abs(wr(n3))) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      work = ar(:,n2)
      ar(:,n2) = ar(:,n3)
      ar(:,n3) = work
     End If
    End Do
   End Do

   EigenValues = WR
   Do i1=1,N1
    Do i2=1,n2
     EigenVectors(i1,i2) = AR(i2,i1) 
    End Do
   End Do

   Do i1=1,n1
    Do i2=1,n1
     testR(i1,i2) = Zero
     Do i3=1,n1
      Do i4=1,n1
       testR(i1,i2) = testR(i1,i2) &
          & + Ar(i3,i1)* ZR_in(i3,i4) *  AR(i4,i2)
      End Do
     End Do
     AbsTest = absTest + Abs( testR(i1,i2) )
    End Do
   End Do
   AbsTest = AbsTest / EigenValues(n1)
   If (Abstest.Gt.1.e-36) Then
     Write(ErrCan,*) "Problem for real diagonlization in"//NameOfUnit(Iname)
     Write(ErrCan,*) "relative precision is ",AbsTest
     Write(ErrCan,*) " "
     kont = -1004
     Call AddError(1004)
   End If

   Deallocate( testR )
  Else ! complex matrix

   Allocate(ZR(2*N1,2*N1))
   Allocate(ZR_in(2*N1,2*N1))
   Allocate(WR(2*N1))
   Allocate(Work(2*N1))

   ZR_in(1:N1,1:N1) = AR
   ZR_in(N1+1:2*N1,N1+1:2*N1) = AR
   ZR_in(N1+1:2*N1,1:N1) = AI
   Do i1=1,n1
    Do i2=1,n1
     ZR_in(i1,N1+i2) = - AI(i1,i2)
    End Do
   End Do

   Call JacobiQP(ZR_in, 2*n1, 2*n1, wr, zr, nrot)

   Do n2=1,2*n1-1
    Do n3=n2+1,2*n1
     If (wr(n2).Gt.wr(n3)) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      work = zr(:,n2)
      zr(:,n2) = zr(:,n3)
      zr(:,n3) = work
     End If
    End Do
   End Do

   Do i1=1,n1
    eigenvalues(i1) = wr(2*i1-1)
    Do i2=1,n1
     eigenvectors(i1,i2) =  zr(i2,2*i1-1) - IOne * zr(n1+i2,2*i1-1)
    End Do
   End Do

   Do i1=2,n1
    Do i2=1,i1-1
     work(1) = eigenvectors(i1,i2)
     eigenvectors(i1,i2) = eigenvectors(i2,i1)
     eigenvectors(i2,i1) = work(1)
    End Do
   End Do

   Do i1=1,n1
    Do i2=1,n1
     test(i1,i2) = Zero
     Do i3=1,n1
      Do i4=1,n1
       test(i1,i2) = test(i1,i2) &
          & + Eigenvectors(i1,i3)* Matrix(i3,i4)* Conjg( EigenVectors(i2,i4))
      End Do
     End Do
    End Do
   End Do

   Deallocate(ZR)

  End If ! decision whether real or complex matrix

  Deallocate(AR,AI,WR,Work, ZR_in)

  Iname = Iname - 1

 End Subroutine ComplexEigenSystem_QP


 Subroutine ComplexEigenSystem_DP1(Matrix, EigenValues, EigenVectors, kont, test, ii)
 !---------------------------------------------------------------------
 ! Subroutine for diagonalization of complex hermitian matrices, based on the
 ! Householder algorithm. Is a portation of  EISCH1 to F90
 ! Input:
 !  Matrix ..... n times n matrix
 ! Output
 !  EigenValues ..... n sorted EigenValues: |m_1| < |m_2| < .. < |m_n|
 !  EigenVectors .... n times n matrix with the eigenvectors
 ! written by Werner Porod, 10.11.2000
 !---------------------------------------------------------------------
 Implicit None
  !-------
  ! input
  !-------
  Complex(Dp), Intent(in) :: Matrix(:,:)
  !--------
  ! output
  !--------
  Complex(Dp), Intent(out) :: EigenVectors(:,:)
  Real(Dp), Intent(out) :: EigenValues(:), test(:)
  Integer, Intent(in) :: ii
  Integer, Intent(inout) :: kont

  !-----------------
  ! local variables
  !-----------------
  Integer :: N1,N2,N3
  Real(qp), Allocatable :: AR(:,:),AI(:,:), WR(:), ZR(:,:),  WORK(:)  &
    & , work2(:,:), ZI(:,:)
  Complex(qp), Allocatable :: Ctest(:,:), Ctest2(:,:), CtestA(:,:)
  Logical :: l_complex = .False.

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ComplexEigenSystem_DP1'

  If (ii.Eq.1) Write(*,*) "ComplexEigenSystem_DP1"
  N1 = Size(Matrix, Dim=1)
  N2 = Size(EigenValues)
  N3 = Size(EigenVectors, Dim=1)
  If ((N1.Ne.N2).Or.(N1.Ne.N3)) Then
   Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'Dimensions to not match: ',N1,N2,N3
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -13
   Call AddError(13)
   Iname = Iname - 1
   Return
  End If

  If (Is_NaN(Real(Matrix,dp)).or.Is_NaN(Aimag(Matrix))) Then !  
   Write(ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write(ErrCan,*) 'matrix contains NaN'
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -31
   Call AddError(31)
   Iname = Iname - 1
   Return 
  End If

  Allocate(AR(N1,N1)) 
  Allocate(AI(N1,N1))
  Allocate(Ctest(N1,N1))
  Allocate(Ctest2(N1,N1))
  Allocate(CtestA(N1,N1))

  AR = Real( Matrix,qp )
  AI = Aimag( Matrix )

  Eigenvectors = ZeroC
  Eigenvalues = 0._dp
  test = 0._dp
  !--------------------------------------------------------------------------
  ! check first whether the matrix is really complex
  ! if not, I use the only real diagonalization because it is more accurate
  !--------------------------------------------------------------------------
  If (Maxval( Abs(AI) ).Eq.0._qp) Then ! real matrix

   Allocate(WR(N1))
   Allocate(Work(N1))

   Call Tred2A_QP(AR, WR, Work)
   Call TQLi_QP(WR,WORK,AR,kont)

   Do n2=1,n1-1
    Do n3=n2+1,n1
     If (wr(n2).Gt.wr(n3)) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      work = ar(:,n2)
      ar(:,n2) = ar(:,n3)
      ar(:,n3) = work
     End If
    End Do
   End Do

   EigenValues = WR
   Ctest = Ar
   EigenVectors = Transpose(AR) 

  Else ! complex matrix
   l_complex = .True.

   Allocate(ZR(N1,N1))
   Allocate(ZI(N1,N1))
   Allocate(WR(N1))
   Allocate(Work(N1))
   Allocate(Work2(2,N1))

   Call HTRIDI_QP(AR, AI, WR, Work, WORK2)
   ZR = 0._qp
   Do n2=1,n1
    ZR(n2,n2) = 1._qp
   End Do
   Call TQLi_QP(WR,WORK,ZR,kont)
   If(KONT/=0) Then
    Iname = Iname - 1
    Deallocate(AR,AI,WR,Work,Ctest)
    Deallocate(ZR, zi, work2)
    Return
   End If
   Call HTRIBK_QP(AR, AI, WORK2, ZR, ZI)
   

   Do n2=1,n1-1
    Do n3=n2+1,n1
     If (wr(n2).Gt.wr(n3)) Then
      work(1) = wr(n2) 
      wr(n2) = wr(n3)
      wr(n3) = work(1)
      work = zr(:,n2)
      zr(:,n2) = zr(:,n3)
      zr(:,n3) = work
      work = zi(:,n2)
      zi(:,n2) = zi(:,n3)
      zi(:,n3) = work
     End If
    End Do
   End Do

   eigenvalues = wr
   eigenvectors = Cmplx(zr,zi, dp)

   Ctest = Cmplx(zr,zi, qp)
   CtestA = Transpose(Ctest)
   CtestA = Conjg(CtestA)
   Eigenvectors = CtestA

   Deallocate(ZR, zi, work2)

  End If ! decision whether real or complex matrix

  !----------------------------------
  ! test of diagonalisation
  !----------------------------------
  Ctest2 = Matmul(Matrix, Ctest )
  ! Ctest^\dagger = Transpose(Cmplx(zr,-zi, qp))
  Ctest = Matmul( CtestA, Ctest2 )

  test = 0._qp
  Do n2=1,n1
   Do n3=1,n1
    If (n2.Eq.n3) Then
     test(1) = Max(test(1), Abs( Ctest(n2,n3) ) )
    Else
     test(2) = Max(test(2), Abs( Ctest(n2,n3) ) )
    End If
   End Do
  End Do
  If (test(1).Gt.0._dp) Then
   If (l_complex) Then
    If ( (test(2)/test(1)).Gt.1.e5_qp*MinimalPrecision) then
     kont = -14
     Call AddError(14)
    End If
   Else 
    If ( (test(2)/test(1)).Gt.MinimalPrecision) then
     kont = -14
     Call AddError(14)
    End If
   End If
  End If

  Deallocate(AR,AI,WR,Work,Ctest,Ctest2,CtestA)

  Iname = Iname - 1

 End Subroutine ComplexEigenSystem_DP1


 Subroutine JacobiQP(a,n,np,d,v,nrot)
 Implicit None
  Integer :: n, np, nrot
  Real(qp) :: a(np,np),d(np),v(np,np)
  Integer, Parameter :: NMAX=500
  Integer :: i, ip, iq, j
  Real(qp) :: c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)

  Do ip=1,n
    Do iq=1,n
      v(ip,iq)=Zero
    End Do
    v(ip,ip)=One
  End Do
  Do ip=1,n
    b(ip)=a(ip,ip)
    d(ip)=b(ip)
    z(ip)=Zero
  End Do
  nrot=0
  Do i=1,50
    sm=Zero
    Do ip=1,n-1
      Do iq=ip+1,n
       sm=sm+Abs(a(ip,iq))
      End Do
    End Do
    If(sm.Eq.Zero)Return
    If(i.Lt.4)Then
      tresh=PointTwo*sm/n**2
    Else
      tresh=Zero
    Endif
    Do ip=1,n-1
     Do iq=ip+1,n
      g=Hundred*Abs(a(ip,iq))
      If(       (i.Gt.4).And.(Abs(d(ip))+g.Eq.Abs(d(ip)))   &
        & .And.(Abs(d(iq))+g.Eq.Abs(d(iq))))Then
        a(ip,iq)=Zero
      Else If(Abs(a(ip,iq)).Gt.tresh)Then
        h=d(iq)-d(ip)
        If(Abs(h)+g.Eq.Abs(h))Then
          t=a(ip,iq)/h
        Else
         theta=PointFive*h/a(ip,iq)
         t=One/(Abs(theta)+Sqrt(One+theta**2))
         If(theta.Lt.Zero)t=-t
        Endif
        c=One/Sqrt(One+t**2)
        s=t*c
        tau=s/(One+c)
        h=t*a(ip,iq)
        z(ip)=z(ip)-h
        z(iq)=z(iq)+h
        d(ip)=d(ip)-h
        d(iq)=d(iq)+h
        a(ip,iq)=Zero
        Do j=1,ip-1
          g=a(j,ip)
          h=a(j,iq)
          a(j,ip)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        End Do
        Do j=ip+1,iq-1
          g=a(ip,j)
          h=a(j,iq)
          a(ip,j)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        End Do
        Do j=iq+1,n
          g=a(ip,j)
          h=a(iq,j)
          a(ip,j)=g-s*(h+g*tau)
          a(iq,j)=h+s*(g-h*tau)
        End Do
        Do j=1,n
          g=v(j,ip)
          h=v(j,iq)
          v(j,ip)=g-s*(h+g*tau)
          v(j,iq)=h+s*(g-h*tau)
        End Do
        nrot=nrot+1
      Endif
     End Do
    End Do
    Do ip=1,n
      b(ip)=b(ip)+z(ip)
      d(ip)=b(ip)
      z(ip)=Zero
    End Do
  End Do
  Write(ErrCan,*) 'too many iterations in JacobiQP'

 End Subroutine JacobiQP

 Function Pythag_QP(a,b)
  Implicit None
  Real(qp), Intent(IN) :: a,b
  Real(qp) :: Pythag_QP
  Real(qp) :: absa, absb
  absa=Abs(a)
  absb=Abs(b)

  If (absa > absb) Then
    Pythag_QP=absa*Sqrt(One+(absb/absa)**2)
  Else
    If (absb == Zero) Then
      Pythag_QP= Zero
    Else
      Pythag_QP=absb*Sqrt(One+(absa/absb)**2)
    End If
  End If
 End Function Pythag_QP

 Function outerprod_QP(a,b)
  Real(qp), Dimension(:), Intent(IN) :: a,b
  Real(qp), Dimension(Size(a),Size(b)) :: outerprod_QP
   outerprod_QP = Spread(a,dim=2,ncopies=Size(b))  &
               & * Spread(b,dim=1,ncopies=Size(a))
 End Function outerprod_QP


 Subroutine RealEigenSystem_DP(Matrix,EigenValues,EigenVectors,kont, test)
 !---------------------------------------------------------------------
 ! Subroutine for diagonalization of real symmetric matrices, based on the
 ! Householder algorithm. Is a portation of  EISRS1 to F90
 ! Input:
 !  Matrix ..... n times n matrix
 ! Output
 !  EigenValues ..... n sorted EigenValues: |m_1| < |m_2| < .. < |m_n|
 !  EigenVectors .... n times n matrix with the eigenvectors
 ! written by Werner Porod, 11.11.2000 
 !---------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: Matrix(:,:)
  Real(dp), Intent(out) :: EigenVectors(:,:), EigenValues(:), test(2)
  Integer, Intent(inout) :: kont

  Integer :: N1,N2,N3, i1, i2, i3, i4, n4
  Real(qp), Allocatable :: AR(:,:), WR(:), WORK(:), testR(:,:), work2(:,:)
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RealEigenSystem_DP'

  kont = 0

  N1 = Size(Matrix, Dim=1)
  N2 = Size(EigenValues)
  N3 = Size(EigenVectors, Dim=1)
  If ((N1.Ne.N2).Or.(N1.Ne.N3)) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) 'Dimensions to not match: ',N1,N2,N3
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -1005
   Call AddError(1005)
   Return 
  End If

  Allocate(AR(N1,N1))
  Allocate(testR(N1,N1))
  Allocate(WR(N1))
  Allocate(Work(N1))
  Allocate(Work2(N1,N1))

  Wr = Zero
  Work = Zero
  AR = Zero
  AR = Real( Matrix, qp )

  Call Tred2A_QP(AR, WR, Work)
  Call TQLi_QP(WR,WORK,AR,kont)

  Do n2=1,n1-1
   Do n3=n2+1,n1
    If (wr(n2).Gt.wr(n3)) Then
     work(1) = wr(n2) 
     wr(n2) = wr(n3)
     wr(n3) = work(1)
     Do n4=1,n1
      work(n4) = ar(n4,n2)
      ar(n4,n2) = ar(n4,n3)
      ar(n4,n3) = work(n4)
     End Do
    End If
   End Do
  End Do

  Do n2=1,n1
   EigenValues(n2) = WR(n2)
   Do n3=1,n1
    EigenVectors(n2,n3) = AR(n3,n2)
   End Do
  End Do

  work2 = Real( Matrix, qp)

   test = 0._dp
   Do i1=1,n1
    Do i2=1,n1
     testR(i1,i2) = Zero
     Do i3=1,n1
      Do i4=1,n1
       testR(i1,i2) = testR(i1,i2) &
          & + ar(i3,i1)* work2(i3,i4) *  ar(i4,i2)
      End Do
     End Do
     If (i1.Eq.i2) Then
      If ( test(1).Lt.Abs( testR(i1,i2) ) ) test(1) = Abs( testR(i1,i2) )
     Else
      If ( test(2).Lt.Abs( testR(i1,i2) ) ) test(2) = Abs( testR(i1,i2) )
     End If 
    End Do
   End Do

  If (test(1).Gt.0._dp) Then
   If ( (test(2)/test(1)).Gt.MinimalPrecision) Then
    kont = -1006
    Call AddError(1006)
   End If
  End If

  Deallocate(AR,WR,Work,testR,work2)

  Iname = Iname - 1

 End Subroutine RealEigenSystem_DP


 Subroutine RealEigenSystem_QP(Matrix,EigenValues,EigenVectors,kont)
 !---------------------------------------------------------------------
 ! Subroutine for diagonalization of real symmetric matrices, based on the
 ! Householder algorithm. Is a portation of  EISRS1 to F90
 ! Input:
 !  Matrix ..... n times n matrix
 ! Output
 !  EigenValues ..... n sorted EigenValues: |m_1| < |m_2| < .. < |m_n|
 !  EigenVectors .... n times n matrix with the eigenvectors
 ! written by Werner Porod, 11.11.2000 
 !---------------------------------------------------------------------
 Implicit None
  Real(qp), Intent(in) :: Matrix(:,:)
  Real(qp), Intent(out) :: EigenVectors(:,:), EigenValues(:)
  Integer, Intent(inout) :: kont

  Integer :: N1,N2,N3, nrot
  Real(qp), Allocatable :: AR(:,:), WR(:), WORK(:,:)
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RealEigenSystem_QP'

  kont = 0

  N1 = Size(Matrix, Dim=1)
  N2 = Size(EigenValues)
  N3 = Size(EigenVectors, Dim=1)
  If ((N1.Ne.N2).Or.(N1.Ne.N3)) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) 'Dimensions to not match: ',N1,N2,N3
   If (ErrorLevel.Ge.-1) Call TerminateProgram
   kont = -1007
   Call AddError(1007)
   Return 
  End If

  Allocate(AR(N1,N1))
  Allocate(WR(N1))
  Allocate(Work(N1,N1))

  Work = Matrix

  Call JacobiQP(Work, n1, n1, wr, ar, nrot)

  Do n2=1,n1-1
   Do n3=n2+1,n1
    If (wr(n2).Gt.wr(n3)) Then
     work(1,1) = wr(n2) 
     wr(n2) = wr(n3)
     wr(n3) = work(1,1)
     work(:,1) = ar(:,n2)
     ar(:,n2) = ar(:,n3)
     ar(:,n3) = work(:,1)
    End If
   End Do
  End Do

  EigenValues = WR
  EigenVectors = Transpose(AR) 

  Deallocate(AR,WR,Work)

  Iname = Iname - 1

 End Subroutine RealEigenSystem_QP

 Subroutine tred2A_QP(a,d,e,novectors)

 Implicit None
  Real(qp), Dimension(:,:), Intent(INOUT) :: a
  Real(qp), Dimension(:), Intent(OUT) :: d, e
  Logical, Optional, Intent(IN) :: novectors

  Integer :: i, j, l, n
  Real(qp) :: f, g, h, hh, scale
  Real(qp), Dimension(Size(a,1)) :: gg
  Logical, Save :: yesvec=.True.

  n = Size(a,1)

  If ((n.Ne.Size(a,2)).Or.(n.Ne.Size(d)).Or.(n.Ne.Size(e)) ) Then
   Write(ErrCan,*) "Error in tred2A_QP",n,Size(a,2),Size(d),Size(e)
   If (ErrorLevel.Gt.-2) Call TerminateProgram
  End If

  If (Present(novectors)) yesvec=.Not. novectors

  Do i=n,2,-1
    l=i-1
    h=0.0_qp
    If (l > 1) Then
      scale=Sum(Abs(a(i,1:l)))
      If (scale == 0._qp) Then
        e(i)=a(i,l)
      Else
        a(i,1:l)=a(i,1:l)/scale
        h=Sum(a(i,1:l)**2)
        f=a(i,l)
        g=-Sign(Sqrt(h),f)
        e(i)=scale*g
        h=h-f*g
        a(i,l)=f-g
        If (yesvec) a(1:l,i)=a(i,1:l)/h
        Do j=1,l
          e(j)=(Dot_product(a(j,1:j),a(i,1:j)) &
          +Dot_product(a(j+1:l,j),a(i,j+1:l)))/h
        End Do
        f=Dot_product(e(1:l),a(i,1:l))
        hh=f/(h+h)
        e(1:l)=e(1:l)-hh*a(i,1:l)
        Do j=1,l
          a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
        End Do
      End If
    Else
      e(i)=a(i,l)
    End If
    d(i)=h
  End Do

  If (yesvec) d(1)=0.0_qp
  e(1)=0.0_qp
  Do i=1,n
    If (yesvec) Then
      l=i-1
      If (d(i) /= 0.0_qp) Then
        gg(1:l)=Matmul(a(i,1:l),a(1:l,1:l))
        a(1:l,1:l)=a(1:l,1:l)-outerprod_QP(a(1:l,i),gg(1:l))
      End If
      d(i)=a(i,i)
      a(i,i)=1.0_qp
      If (i.Gt.1) Then
       a(i,1:l)=0.0_qp
       a(1:l,i)=0.0_qp
      End If
    Else
      d(i)=a(i,i)
    End If
  End Do

 End Subroutine tred2A_QP


 Subroutine tqli_QP(d,e,z,kont)

 Implicit None
  Integer, Intent(inout) :: kont
  Real(qp), Dimension(:), Intent(INOUT) :: d,e
  Real(qp), Dimension(:,:), Optional, Intent(INOUT) :: z
  Integer :: i,iter,l,m,n
  Real(qp) :: b,c,dd,f,g,p,r,s
  Real(qp), Dimension(Size(e)) :: ff

  n = Size(d)
  kont = 0
  If (n.Ne.Size(e)) Then
   Write(ErrCan,*) "Error in tqli_QP",n,Size(e)
   If (ErrorLevel.Gt.-2) Call TerminateProgram
   kont = -17
   Call AddError(17)
  End If

  If (Present(z)) Then
   If ((n.Ne.Size(z,dim=1)).Or.(n.Ne.Size(z,dim=2)) ) Then
    Write(ErrCan,*) "Error in tqli_QP",n,Size(z,dim=1),Size(z,dim=2)
    If (ErrorLevel.Gt.-2) Call TerminateProgram
    kont = -17
    Call AddError(17)
   End If
  End If

  e(:)=Eoshift(e(:),1)

  Do l=1,n
    iter=0
    iterate: Do
      Do m=l,n-1
        dd=Abs(d(m))+Abs(d(m+1))
        If (Abs(e(m))+dd == dd) Exit
      End Do
      If (m == l) Exit iterate
      If (iter == 30) Then
       Write(ErrCan,*) "Problem in tqli_QP, too many iterations"
       kont = -18
       Call AddError(18)
       Return
      End If
      iter=iter+1
      g=(d(l+1)-d(l))/(2.0_qp*e(l))
      r=pythag_qp(g,1.0_qp)
      g=d(m)-d(l)+e(l)/(g+Sign(r,g))
      s=1.0_qp
      c=1.0_qp
      p=0.0_qp
      Do i=m-1,l,-1
        f=s*e(i)
        b=c*e(i)
        r=pythag_QP(f,g)
        e(i+1)=r
        If (r == 0.0_qp) Then
          d(i+1)=d(i+1)-p
          e(m)=0.0_qp
          Cycle iterate
        End If
        s=f/r
        c=g/r
        g=d(i+1)-p
        r=(d(i)-g)*s+2.0_qp*c*b
        p=s*r
        d(i+1)=g+p
        g=c*r-b
        If (Present(z)) Then
          ff(1:n)=z(1:n,i+1)
          z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
          z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
        End If
      End Do
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.0_qp
    End Do iterate
  End Do

 End Subroutine tqli_QP

  
  Subroutine HTRIBK_QP(AR, AI, TAU, ZR, ZI)
  Implicit None
   Real(qp) :: AR(:,:), AI(:,:), TAU(:,:), ZR(:,:), ZI(:,:)

   Integer :: n, m, k, j, i, l
   Real(qp) :: s, si, h

   n = Size(ar,2)
   m = Size(zr,2)

   Do K = 1, N
    Do J = 1, M
      ZI(K,J) = - ZR(K,J) * TAU(2,K)
      ZR(K,J) = ZR(K,J) * TAU(1,K)
    Enddo
   Enddo
   If (N == 1) Return
  Do I = 2, N
    L = I - 1
    H = AI(I,I)
    If (H == 0._qp) Cycle
    Do J = 1, M
      S = 0._qp
      SI = 0._qp
      Do K = 1, L
        S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
        SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
      Enddo
      S = S / H
      SI = SI / H
      Do K = 1, L
        ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
        ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
      Enddo
    Enddo
  End Do

  End Subroutine HTRIBK_QP



  Subroutine HTRIDI_QP(AR, AI, D, E, TAU)
  Implicit None
   Real(qp) :: AR(:,:), AI(:,:), D(:), E(:), TAU(:,:)

   Integer :: n, i, ii, l, k, j, jp1
   Real(qp) :: scalei, h, g, f, si, hh, fi, gi

   n = Size(ar,2)

   TAU(1,N) = 1._qp
   TAU(2,N) = 0._qp
   Do I = 1, N
    D(I) = AR(I,I)
   End Do

   Do II = 1, N
    I = N + 1 - II
    L = I - 1
    H = 0._qp
    SCALEI = 0._qp
    If (L < 1) GO TO 130
    Do K = 1, L
      SCALEI = SCALEI + Abs(AR(I,K)) + Abs(AI(I,K))
    Enddo
    If (SCALEI /= 0._qp) GO TO 140
    TAU(1,L) = 1._qp
    TAU(2,L) = 0._qp
  130 E(I) = 0._qp
    GO TO 290
  140 Do K = 1, L
      AR(I,K) = AR(I,K) / SCALEI
      AI(I,K) = AI(I,K) / SCALEI
      H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
    Enddo
    G = Sqrt(H)
    E(I) = SCALEI * G
    F = Abs(Cmplx(AR(I,L),AI(I,L),qp))
    If (F == 0._qp) GO TO 160
    TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
    SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
    H = H + F * G
    G = 1._qp + G / F
    AR(I,L) = G * AR(I,L)
    AI(I,L) = G * AI(I,L)
    If (L == 1) GO TO 270
    GO TO 170
  160 TAU(1,L) = -TAU(1,I)
    SI = TAU(2,I)
    AR(I,L) = G
  170 F = 0._qp
    Do J = 1, L
      G = 0._qp
      GI = 0._qp
      Do K = 1, J
        G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
        GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
      Enddo
      JP1 = J + 1
      If (L < JP1) GO TO 220
      Do K = JP1, L
        G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
        GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
      Enddo
  220 E(J) = G / H
      TAU(2,J) = GI / H
      F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
    Enddo
    HH = F / (H + H)
    Do J = 1, L
      F = AR(I,J)
      G = E(J) - HH * F
      E(J) = G
      FI = -AI(I,J)
      GI = TAU(2,J) - HH * FI
      TAU(2,J) = -GI
      Do K = 1, J
        AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K) &
          + FI * TAU(2,K) + GI * AI(I,K)
        AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K) &
          - FI * E(K) - GI * AR(I,K)
      Enddo
    Enddo
  270 Do K = 1, L
      AR(I,K) = SCALEI * AR(I,K)
      AI(I,K) = SCALEI * AI(I,K)
    Enddo
    TAU(2,L) = -SI
  290 HH = D(I)
    D(I) = AR(I,I)
    AR(I,I) = HH
    AI(I,I) = SCALEI * SCALEI * H
  Enddo
  Return

  End Subroutine HTRIDI_QP



End Module MathematicsQP
