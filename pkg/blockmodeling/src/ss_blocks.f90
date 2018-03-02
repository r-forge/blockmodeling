subroutine critfunsscom(M,n,clu,k,diag,err,E,BM)
	INTEGER n, clu, k, i, j, nA, nAD
	DOUBLE PRECISION M, E, BM, sumA, sumAD, sumA2, sumAD2, err, mean
	LOGICAL diag
	DIMENSION M(n,n), clu(n), E(k,k), BM(k,k), sumA(k,k), sumAD(k),sumA2(k,k), sumAD2(k), nA(k,k), nAD(k)
	
	do i = 1, k
		nAD(i) = 0
		sumAD(i) = 0.0
		sumAD2(i) = 0.0
		do j = 1, k
			 nA(i,j) = 0
			 sumA(i,j) = 0.0
			 sumA2(i,j) = 0.0
		end do
	end do
	
	
	do i = 1, n
		do j = 1, n
			if((.not.diag) .or. (i.ne.j)) then
				nA(clu(i),clu(j)) = nA(clu(i),clu(j)) + 1
				sumA(clu(i),clu(j)) = sumA(clu(i),clu(j))+M(i,j)
				sumA2(clu(i),clu(j)) = sumA2(clu(i),clu(j))+M(i,j)**2
				else 
				nAD(clu(i)) = nAD(clu(i)) + 1
				sumAD(clu(i)) = sumAD(clu(i))+M(i,i)
				sumAD2(clu(i)) = sumAD2(clu(i))+M(i,i)**2
			endif
		end do
	end do
	
	err = 0.0

	do i = 1, k
		do j = 1, k
			if(diag.and.i.eq.j) then
				if(nA(i,j).eq.0) then
					nA(i,j) = 1
				end if
				mean = sumA(i,j)/nA(i,j)
				E(i,j) = sumA2(i,j)-nA(i,j)*mean**2 + sumAD2(i)-sumAD(i)**2/nAD(i)
				BM(i,j) = mean
				err = err + E(i,j)
			else 
				mean = sumA(i,j)/nA(i,j)
				E(i,j) = sumA2(i,j)-nA(i,j)*mean**2
				BM(i,j) = mean
				err = err + E(i,j)
			endif
		end do
	end do
end










 	
subroutine sscom(B,n1,n2,diag,sserr,mean)
	INTEGER n1, n2, ia, iad
	DOUBLE PRECISION B, A, AD, ss, sserr, mean, temp
	LOGICAL diag
	DIMENSION B(n1,n2), A(n1*n2), AD(n1)

	ia = 0
	iad = 0
	do i=1, n1
	do j=1, n2
	if(.not.diag.or.i.ne.j) then
	ia = ia + 1
	A(ia) = B(i,j)
	else 
	iad = iad + 1
	AD(iad) = B(i,j)
	endif
	end do
	end do

	if(diag) then 
	sserr = ss(A,n1*(n1-1),mean) + ss(AD,n1,temp)
	else 
	sserr = ss(A,n1*n2,mean)
	endif
end

DOUBLE PRECISION FUNCTION ss(a,n,m)
      INTEGER n, i
      DOUBLE PRECISION a, s, m
      DIMENSION a(n)  
      s = 0.0
      do i=1, n
      s = s + a(i)
      end do
      m = s / n
      ss = 0.0
      do i=1, n
      ss = ss + (a(i)-m)**2
      end do
end

