! REGDI.FOR 3/18/85 - DOUG WHITE'S REGULAR DISTANCES PROGRAM
      subroutine regd(R,B,N,NR,ITER)
      DOUBLE PRECISION   R, B, DEG, SUM, CM
      INTEGER NR, N, ITER, KR, JJ, II
      DIMENSION  DEG (N), SUM (N,N), R (N,N, NR), B (N,N)

!     COMPUTE DEGREE, SUMS FOR I-->K, INITIAL STRUCTURAL DISTANCE
      DO 100 I=1,N
      DEG(I)=0.0
      DO 100 J=1,N
      SUM(I,J)=0.0
      DO 50 KR=1,NR
      SM = R(I,J,KR)**2 + R(J,I,KR)**2
   50 SUM(I,J)=SUM(I,J) + sm
  100 DEG(I)=DEG(I)+SUM(I,J)

      IQUIT=0

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
!     IF(DEG(J).EQ.0.0) GO TO 506
      I=II
! TAKE EACH OF THE TWO POINTS AS REFERENT
      DO 505 IJ=1,2
      IF (IJ.EQ.1) GOTO 120
      J=II
      I=JJ
!     TAKE POINT K (I-->K, K-->I TO BE EXAMINED)
!                   J-->K, K-->J IN SECOND ITERATION
 120  DO 500 K=1,N
      IF(SUM(I,K).EQ.0.0) GO TO 500
      XMIN=10000000000.0
!     FIND BEST MATCHING POINT M
      DO 400 M=1,N
! 0 should be allowed as a best fit for small values      IF(SUM(J,M).EQ.0.0) GO TO 400
      SUMM=0.0
      DO 300 KR=1,NR
300   summ = summ + (R(I,K,KR) - R(J,M,KR)) **2 + (R(K,i,KR) - R(M,j,KR)) **2 
      CMIKJM = max (Summ, sum(i,k) * b (max (k,m), min (k,m)))
!     IF PERFECT MATCH DESIRED, CORRECT MATCH
!     IF(SUMM.NE.SUM(I,K).AND.NOERRS.EQ.1)  CMIKJM=DEG(II)+DEG(JJ)
      IF(CMIKJM.LT.XMIN) XMIN= CMIKJM
 
      IF(XMIN.EQ.0) GO TO 450
 
  400 CONTINUE
!     ADD BEST MATCHES TO REGULAR DISTANCE NUMERATOR FOR I,J
  450 CM=CM+XMIN
  500  CONTINUE
  505  CONTINUE
!     COMPUTE REGULAR DISTANCE
  506 DM = DEG(II)+DEG(JJ)
! REMEMBER BOTH POINTS TAKEN AS REFERENCE
      if(cm.gt.dm) cm=DM
      IF(DM.NE.0.0) B (II,JJ)=CM/DM
!     IF(B (II,JJ).LE.CUT) B (II,JJ)=0.0
!      DIFF = B(II,JJ) - B (JJ,II)
!      IF(DIFF.LT.0.0) DIFF = -DIFF
!      D = D + DIFF
 510  CONTINUE
 520  CONTINUE

!     (D.EQ.0.0.AND.L.NE.1).OR.
      IF(L.EQ.ITER) IQUIT=1
! symmetrize : to lower half matrix
      DO 650 I = 2, N
      DO 600 J = 1, i-1
  600 B(i,j) = B(j,i) 
  650 CONTINUE
  
      IF(IQUIT.EQ.1) GO TO 800
 
  700 CONTINUE
  800 CONTINUE
 
      END
