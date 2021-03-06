! REGE_NE_R.F Ales Ziberna, 2006 - NORMALIZED EQUIVALENCES NORMALIZED MATRICES version of REGE (Douglas R. White, 1985)
!  THIS VERSION ALLOWS USER TO SET THE NUMBER OF ITERATIONS 
      subroutine regene(R,B,N,NR,ITER)
      DOUBLE PRECISION   R, B, DEG, SUM, xxmax, row, col, SUMM, CMIKJM, DM, XMAX, CM
      INTEGER NR, N, ITER, KR, JJ, II, NumIter
      DIMENSION  DEG (N), SUM (N,N), R (N,N, NR), B (N,N), row(N), col(N)

!     COMPUTE DEGREE, SUMS FOR I--&gt;K, INITIAL STRUCTURAL EQUIV.
      DO 100 I=1,N
      DEG(I)=0.0
      DO 99 J=1,N
      SUM(I,J)=0.0
      DO 50 KR=1,NR
      SUM(I,J)=SUM(I,J)+R(I,J,KR)+R(J,I,KR)
   50 END DO
      DEG(I)=DEG(I)+SUM(I,J)
 99   END DO
100   END DO
!     BEGIN ITERATIONS
      DO 700 L=1,ITER
!     INITIALIZE DIFFERENCE IN SUCCESSIVE SE MATRICES
      D = 0.0
!     TAKE POINT I
      DO 520 II = 1, N-1
      I=II
!     IF DEGREE ZERO NEXT I
!     IF(DEG(I).EQ.0.0) GO TO 520
!     TAKE POINT J
      DO 510 JJ= II+1, N
      CM      = 0.0
      J=JJ
!     IF DEGREE ZERO NEXT J
      IF(DEG(J).EQ.0.0) GO TO 506
      I=II
! TAKE EACH OF THE TWO POINTS AS REFERENT
      DO 505 IJ=1,2
      IF (IJ.EQ.1) GOTO 120
      J=II
      I=JJ
!     TAKE POINT K (I--&gt;K, K--&gt;I TO BE EXAMINED)
 120  DO 500 K=1,N
      IF(SUM(I,K).EQ.0.0) GO TO 500
      XMAX=0.0
!     FIND BEST MATCHING POINT M
      DO 400 M=1,N
      IF(SUM(J,M).EQ.0.0) GO TO 400
      SUMM=0.0
      DO 300 KR=1,NR
      SUMM = SUMM +min (R(I,K,KR),r(j,m,kr)) +min (R(K,I,KR),r(m,j,kr))
  300 CONTINUE
      CMIKJM = SUMM * b (max (k,m), min (k,m))
!     IF PERFECT MATCH DESIRED, CORRECT MATCH
!     IF(SUMM.NE.SUM(I,K).AND.NOERRS.EQ.1)  CMIKJM=0.0
      IF(CMIKJM.GT.XMAX) XMAX= CMIKJM
 
      IF(XMAX.EQ.SUM(I,K)) GO TO 450
 
  400 CONTINUE
!     ADD BEST MATCHES TO REGULAR EQUIVALENCE NUMERATOR FOR I,J
  450 CM=CM+XMAX
  500  CONTINUE
  505  CONTINUE
!     COMPUTE REGULAR EQUIVALENCE
  506 DM = DEG(II)+DEG(JJ)
      B (II,JJ)= 1.0
      IF(DM.NE.0.0) B (II,JJ)=CM/DM
!     IF(B (II,JJ).LE.CUT) B (II,JJ)=0.0
  510  CONTINUE
  520  CONTINUE

! symmetrize : to lower half matrix
      DO 600 I = 2, N
      DO 599 J = 1, i-1
      B(i,j) = B(j,i) 
 599  END DO
 600  END DO
! Start normalization
      NumIter=15
      DO K = 1, NumIter
        Xxmax=0.0
  
  
! compute row and col totals of B
      DO I = 1, N
        B(I,I)=0.0
        Row(i)= 0.0
        Col(i)= 0.0  
      ENDDO
      
      DO I = 1, N
      DO J = 1, N
        IF (xxmax.lt.B(I,J)) then
        xxmax=B(I,J)
        ENDIF
        Row(i)= Row(i)+B(I,j)
        Col(j)= Col(j)+B(I,j)
      ENDDO
      ENDDO
      
! normalize the B matrix and symmetrize
      DO I = 2, N
      DO J = 1, i-1
        If (row(i).gt.0.and.col(j).gt.0) then
       B(I,j)=(B(I,j)/Row(i)**.5) /col(j)**.5
         B(J,I)=B(I,J)
        ENDIF
        ENDDO
        ENDDO

        ENDDO ! end of normalization
      
      DO I = 1, N
        B(I,I)=xxmax
        ENDDO
  
  
  700 CONTINUE

      END
