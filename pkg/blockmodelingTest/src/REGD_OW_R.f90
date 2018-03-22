! REGD_OW_R.F Ales Ziberna, 2006 - ONEWAY version of REGD (Douglas R. White, 1985)
      subroutine regdow(R,B,N,NR,ITER)
      DOUBLE PRECISION   R, B, DEG, SUM, SUMM1, SUMM2, XMIN1, XMIN2, CMIKJM1, CMIKJM2, CM
      INTEGER NR, N, ITER, KR, JJ, II
      DIMENSION  DEG (N), SUM (N,N), R (N,N, NR), B (N,N)

!     COMPUTE DEGREE, SUMS FOR I-->K, INITIAL STRUCTURAL DISTANCE
      DO 100 I=1,N
      DO 99 J=1,N
      SUM(I,J)=0.0
      DO 50 KR=1,NR
      SM = R(I,J,KR)**2
      SUM(I,J)=SUM(I,J) + SM
   50 END DO
   99 END DO
  100 END DO
      DO 102 I=1,N
      DEG(I)=0.0
      DO 101 J=1,N
      DEG(I)=DEG(I)+SUM(I,J)+SUM(J,I)
  101 END DO
  102 END DO
!      IQUIT=0

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
      IF((SUM(I,K)+SUM(K,I)).EQ.0.0) GO TO 500
      XMIN1=10000000000.0
      XMIN2=10000000000.0
!     FIND BEST MATCHING POINT M
      DO 400 M=1,N
! 0 should be allowed as a best fit for small values      IF((SUM(J,M)+SUM(M,J)).EQ.0.0) GO TO 400
      SUMM1=0.0
      SUMM2=0.0
      DO 300 KR=1,NR
      IF(R(I,K,KR).NE.0.0) summ1 = summ1 + (R(I,K,KR) - R(J,M,KR)) **2
      IF(R(K,I,KR).NE.0.0) summ2 = summ2 + (R(K,I,KR) - R(M,J,KR)) **2
300   END DO
      CMIKJM1 = max (summ1, sum(i,k) * b (max (k,m), min (k,m)))
      CMIKJM2 = max (summ2, sum(k,i) * b (max (k,m), min (k,m)))
!     IF PERFECT MATCH DESIRED, CORRECT MATCH
!     IF(SUMM.NE.SUM(I,K).AND.NOERRS.EQ.1)  CMIKJM=DEG(II)+DEG(JJ)
      IF(CMIKJM1.LT.XMIN1) XMIN1= CMIKJM1
      IF(CMIKJM2.LT.XMIN2) XMIN2= CMIKJM2
!      call intpr("I",-1,I,1)
!      call intpr("K",-1,K,1)
!      call intpr("J",-1,J,1)
!      call intpr("M",-1,M,1)
!      call dblepr("XMIN1",-1,XMIN1,1)
!      call dblepr("XMIN2",-1,XMIN2,1)

      

      
 
      IF((XMIN1+XMIN2).EQ.0) GO TO 450
 
  400 CONTINUE
!     ADD BEST MATCHES TO REGULAR DISTANCE NUMERATOR FOR I,J
  450 CM=CM+XMIN1+XMIN2
!       call dblepr("CM",-1,CM,1)
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
!	(D.EQ.0.0.AND.L.NE.1).OR.
! symmetrize : to lower half matrix
      DO 650 I = 2, N
      DO 600 J = 1, i-1
      B(i,j) = B(j,i) 
  600 END DO
  650 CONTINUE
      

  700 CONTINUE
  800 CONTINUE
 
      END
