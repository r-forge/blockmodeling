subroutine optparsscomtm(M,clu1,clu2,maxiter,n1,n2,k1,k2,err,E,BM,cluM1,cluM2,nbest,iter,printIter)
	INTEGER iter, maxiter, n1, n2, clu1, clu2, nclu1, nclu2, bnclu1, bnclu2, iclu, k1, k2, i,i2, j, ii, jj
	INTEGER nbest, nA, tnA, bnA, bclu1, bclu2, oldiclu, newiclu, cluM1, cluM2
	DOUBLE PRECISION M, E, BM, err, mean, sumA, sumA2, tsumA, tsumA2, tE, tBM, berr, terr, bsumA, bsumA2, bE, bBM
	LOGICAL imp, printIter
	DIMENSION M(n1,n2), clu1(n1), clu2(n2), E(k1,k2), BM(k1,k2), cluM1(50,n1), cluM2(50,n2), sumA(k1,k2), sumA2(k1,k2), nA(k1,k2)
	DIMENSION tE(k1,k2), tBM(k1,k2),  tsumA(k1,k2), tsumA2(k1,k2),tnA(k1,k2), bE(k1,k2), bBM(k1,k2), bclu1(n1), bclu2(n2)
	DIMENSION bsumA(k1,k2), bsumA2(k1,k2), bnA(k1,k2)
	DIMENSION nclu1(k1), nclu2(k2), bnclu1(k1), bnclu2(k2)

	do i = 1, k2
		nclu2(i) = 0
	end do
	
	do i = 1, k1
		nclu1(i) = 0
		do j = 1, k2
			 nA(i,j) = 0
			 sumA(i,j) = 0.0
			 sumA2(i,j) = 0.0
		end do
	end do
	
	do i = 1, n2
		nclu2(clu2(i)) = nclu2(clu2(i)) + 1
	end do
	
	do i = 1, n1
		nclu1(clu1(i)) = nclu1(clu1(i)) + 1
		do j = 1, n2
			nA(clu1(i),clu2(j)) = nA(clu1(i),clu2(j)) + 1
			sumA(clu1(i),clu2(j)) = sumA(clu1(i),clu2(j))+M(i,j)
			sumA2(clu1(i),clu2(j)) = sumA2(clu1(i),clu2(j))+M(i,j)**2			
		end do
	end do	


	err = 0.0

	do i = 1, k1
		do j = 1, k2
			mean = sumA(i,j)/nA(i,j)
			E(i,j) = sumA2(i,j)-nA(i,j)*mean**2
			BM(i,j) = mean
			err = err + E(i,j)
		end do
	end do
	
	nbest = 1
	berr=err
	do i = 1, n1
		bclu1(i)=clu1(i)
		cluM1(nbest,i) = clu1(i)
	end do
	do i = 1, n2
		bclu2(i)=clu2(i)
		cluM2(nbest,i) = clu2(i)
	end do
	
	do i = 1, k2
		bnclu2(i) = nclu2(i)
	end do
	
	do i = 1, k1
		bnclu1(i) = nclu1(i)	
		do j = 1, k2
			bnA(i,j) = nA(i,j)
			bsumA(i,j) = sumA(i,j)
			bsumA2(i,j) = sumA2(i,j)
			bE(i,j)=E(i,j)
			bBM(i,j)=BM(i,j)
		end do
	end do
	
	imp = .TRUE.
	iter = 0

	if (printIter) then	
		call intpr("iter", -1, iter, 1) 
		call intpr("nclu1", -1, nclu1, k1) 
		call intpr("clu1", -1, clu1, n1) 	
		call intpr("nclu2", -1, nclu2, k2) 
		call intpr("clu2", -1, clu2, n2) 			
		call dblepr("err", -1, err, 1) 
		call dblepr("sumA", -1, sumA, k1*k2) 
		call dblepr("sumA2", -1, sumA2, k1*k2) 

	end if
		
	do while (imp .AND. (iter .LT. maxiter))	
		imp = .FALSE.
		iter = iter + 1
		
		do i = 1, n1
			oldiclu=clu1(i)
			if(nclu1(oldiclu).gt.1) then
				do iclu = 1, k1
					if(oldiclu.ne.iclu) then
						newiclu = iclu

						do ii = 1, k1
							do jj = 1, k2
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do
						
						do jj = 1, k2
							tnA(oldiclu,jj) = tnA(oldiclu,jj) - nclu2(jj)
							tnA(newiclu,jj) = tnA(newiclu,jj) + nclu2(jj)
						end do
						
						do j = 1 ,n2
							tsumA(oldiclu,clu2(j)) = tsumA(oldiclu,clu2(j)) - M(i,j)
							tsumA(newiclu,clu2(j)) = tsumA(newiclu,clu2(j)) + M(i,j)
							tsumA2(oldiclu,clu2(j)) = tsumA2(oldiclu,clu2(j)) - M(i,j)**2
							tsumA2(newiclu,clu2(j)) = tsumA2(newiclu,clu2(j)) + M(i,j)**2
						end do


						terr = 0.0

						do ii = 1, k1
							do jj = 1, k2
								if(ii.eq.newiclu .or. ii.eq.oldiclu) then
									mean = tsumA(ii,jj)/tnA(ii,jj)
									tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
									tBM(ii,jj) = mean
								else
									tE(ii,jj) = E(ii,jj)
									tBM(ii,jj) = BM(ii,jj)
								end if

								terr = terr + tE(ii,jj)
							end do
						end do

						if (terr.lt.berr) then
!							call intpr("Move", -1, 1, 1) 
							imp = .TRUE.
							berr=terr
							do ii = 1, n1
								bclu1(ii)=clu1(ii)
							end do

							do ii = 1, k1
								bnclu1(ii) = 0
								do jj = 1, k2
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj)=tE(ii,jj)
									bBM(ii,jj)=tBM(ii,jj)
								end do
							end do

							bclu1(i) = newiclu

							nbest = 1
							do ii = 1, n1
								bnclu1(bclu1(ii)) = bnclu1(bclu1(ii)) + 1
								cluM1(nbest,ii) = bclu1(ii)
							end do

							do jj = 1, n2
								cluM2(nbest,jj) = clu2(jj)
							end do							
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do ii = 1, n1
										cluM1(nbest,ii) = clu1(ii)
									end do
									cluM1(nbest,i) = newiclu
									do jj = 1, n2
										cluM2(nbest,jj) = clu2(jj)
									end do							
								endif
							end if
						end if

					end if
				end do
			end if

			do i2 = 1, i - 1
				if(i.ne.i2) then
					if(clu1(i).ne.clu1(i2)) then
						!tclu = clu
						!tclu(i) = clu(j)
						!tclu(j) = clu(i)
						newiclu = clu1(i2)

						do ii = 1, k1
							do jj = 1, k2
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do

						do jj = 1, n2
							tsumA(oldiclu,clu2(jj)) = tsumA(oldiclu,clu2(jj)) - M(i,jj) + M(i2,jj)
							tsumA(newiclu,clu2(jj)) = tsumA(newiclu,clu2(jj)) + M(i,jj) - M(i2,jj)
							tsumA2(oldiclu,clu2(jj)) = tsumA2(oldiclu,clu2(jj)) - M(i,jj)**2 + M(i2,jj)**2
							tsumA2(newiclu,clu2(jj)) = tsumA2(newiclu,clu2(jj)) + M(i,jj)**2 - M(i2,jj)**2
						end do
						

						terr = 0.0

						do ii = 1, k1
							do jj = 1, k2
								if(ii.eq.newiclu .or. ii.eq.oldiclu) then
									mean = tsumA(ii,jj)/tnA(ii,jj)
									tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
									tBM(ii,jj) = mean
								else
									tE(ii,jj) = E(ii,jj)
									tBM(ii,jj) = BM(ii,jj)
								end if

								terr = terr + tE(ii,jj)
							end do
						end do


						if(terr.lt.berr) then
!							call intpr("Switch", -1, 2, 1) 
							call dblepr("terr", -1, terr, 1)
							call dblepr("berr", -1, berr, 1)
							call dblepr("diff", -1, berr-terr, 1)
							imp = .TRUE.
							berr=terr
							do ii = 1, n1
								bclu1(ii)=clu1(ii)
							end do
							bclu1(i) = newiclu
							bclu1(i2) = oldiclu

							do ii = 1, k1
								bnclu1(ii) = 0
								do jj = 1, k2
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj) = tE(ii,jj)
									bBM(ii,jj) = tBM(ii,jj)
								end do
							end do
							
							nbest = 1
							do ii = 1, n1
								bnclu1(bclu1(ii)) = bnclu1(bclu1(ii)) + 1
								cluM1(nbest,ii) = bclu1(ii)
							end do

							do jj = 1, n2
								cluM2(nbest,jj) = clu2(jj)
							end do							
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do ii = 1, n1
										cluM1(nbest,ii) = clu1(ii)
									end do
									cluM1(nbest,i) = newiclu
									cluM1(nbest,i2) = oldiclu
									do jj = 1, n2
										cluM2(nbest,jj) = clu2(jj)
									end do							
								endif
							end if
						end if
					end if			
				end if
			end do
		end do	


		do j = 1, n2
			oldiclu=clu2(j)
			if(nclu2(oldiclu).gt.1) then
				do iclu = 1, k2
					if(oldiclu.ne.iclu) then
						newiclu = iclu

						do ii = 1, k1
							do jj = 1, k2
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do
						
						do ii = 1, k1
							tnA(ii,oldiclu) = tnA(ii,oldiclu) - nclu1(ii)
							tnA(ii,newiclu) = tnA(ii,newiclu) + nclu1(ii)
						end do
						
						do i = 1 ,n1
							tsumA(clu1(i),oldiclu) = tsumA(clu1(i),oldiclu) - M(i,j)
							tsumA(clu1(i),newiclu) = tsumA(clu1(i),newiclu) + M(i,j)
							tsumA2(clu1(i),oldiclu) = tsumA2(clu1(i),oldiclu) - M(i,j)**2
							tsumA2(clu1(i),newiclu) = tsumA2(clu1(i),newiclu) + M(i,j)**2
						end do


						terr = 0.0

						do ii = 1, k1
							do jj = 1, k2
								if(jj.eq.newiclu .or. jj.eq.oldiclu) then
									mean = tsumA(ii,jj)/tnA(ii,jj)
									tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
									tBM(ii,jj) = mean
								else
									tE(ii,jj) = E(ii,jj)
									tBM(ii,jj) = BM(ii,jj)
								end if

								terr = terr + tE(ii,jj)
							end do
						end do

						if (terr.lt.berr) then
!							call intpr("Move", -1, 1, 1) 
							imp = .TRUE.
							berr=terr
							do jj = 1, n2
								bclu2(jj)=clu2(jj)
							end do

							do jj = 1, k2
								bnclu2(jj) = 0
								do ii = 1, k1
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj)=tE(ii,jj)
									bBM(ii,jj)=tBM(ii,jj)
								end do
							end do

							bclu2(j) = newiclu

							nbest = 1
							do jj = 1, n2
								bnclu2(bclu2(jj)) = bnclu2(bclu2(jj)) + 1
								cluM2(nbest,jj) = bclu2(jj)
							end do

							do ii = 1, n1
								cluM1(nbest,ii) = clu1(ii)
								bclu1(ii)=clu1(ii)
							end do							
							do ii = 1, k1
								bnclu1(ii)=nclu1(ii)
							end do							
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do jj = 1, n2
										cluM2(nbest,jj) = clu2(jj)
									end do
									cluM2(nbest,j) = newiclu
									do ii = 1, n1
										cluM1(nbest,ii) = clu1(ii)
									end do							
								endif
							end if
						end if

					end if
				end do
			end if

			do j2 = 1, j - 1
				if(j.ne.j2) then
					if(clu2(j).ne.clu2(j2)) then
						!tclu = clu
						!tclu(i) = clu(j)
						!tclu(j) = clu(i)
						newiclu = clu2(j2)

						do ii = 1, k1
							do jj = 1, k2
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do

						do ii = 1, n1
							tsumA(clu1(ii),oldiclu) = tsumA(clu1(ii),oldiclu) - M(ii,j) + M(ii,j2)
							tsumA(clu1(ii),newiclu) = tsumA(clu1(ii),newiclu) + M(ii,j) - M(ii,j2)
							tsumA2(clu1(ii),oldiclu) = tsumA2(clu1(ii),oldiclu) - M(ii,j)**2 + M(ii,j2)**2
							tsumA2(clu1(ii),newiclu) = tsumA2(clu1(ii),newiclu) + M(ii,j)**2 - M(ii,j2)**2
						end do
						

						terr = 0.0

						do ii = 1, k1
							do jj = 1, k2
								if(jj.eq.newiclu .or. jj.eq.oldiclu) then
									mean = tsumA(ii,jj)/tnA(ii,jj)
									tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
									tBM(ii,jj) = mean
								else
									tE(ii,jj) = E(ii,jj)
									tBM(ii,jj) = BM(ii,jj)
								end if

								terr = terr + tE(ii,jj)
							end do
						end do

						if(terr.lt.berr) then
!							call intpr("Switch", -1, 2, 1) 
							call dblepr("terr", -1, terr, 1)
							call dblepr("berr", -1, berr, 1)
							call dblepr("diff", -1, berr-terr, 1)
							imp = .TRUE.
							berr=terr
							do jj = 1, n2
								bclu2(jj)=clu2(jj)
							end do
							bclu2(j) = newiclu
							bclu2(j2) = oldiclu

							do jj = 1, k2
								bnclu2(jj) = 0
								do ii = 1, k1
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj) = tE(ii,jj)
									bBM(ii,jj) = tBM(ii,jj)
								end do
							end do
							
							nbest = 1
							do jj = 1, n2
								bnclu2(bclu2(jj)) = bnclu2(bclu2(jj)) + 1
								cluM2(nbest,jj) = bclu2(jj)
							end do

							do ii = 1, n1
								cluM1(nbest,ii) = clu1(ii)
								bclu1(ii)=clu1(ii)
							end do		
							do ii = 1, k1
								bnclu1(ii)=nclu1(ii)
							end do														
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do jj = 1, n2
										cluM2(nbest,jj) = clu2(jj)
									end do
									cluM2(nbest,j) = newiclu
									cluM2(nbest,j2) = oldiclu
									do ii = 1, n1
										cluM1(nbest,ii) = clu1(ii)
									end do							
								endif
							end if
						end if
					end if			
				end if
			end do
		end do	




		err=berr
		do i = 1, n1
			clu1(i)=bclu1(i)
		end do
		do i = 1, n2
			clu2(i)=bclu2(i)
		end do

		do i = 1, k2
			nclu2(i) = bnclu2(i)
		end do

		do i = 1, k1
			nclu1(i) = bnclu1(i)
			do j = 1, k2
				nA(i,j) = bnA(i,j)
				sumA(i,j) = bsumA(i,j)
				sumA2(i,j) = bsumA2(i,j)
				E(i,j) = bE(i,j)
				BM(i,j) = bBM(i,j)
			end do
		end do	
		if (printIter .OR. (iter.GT.(maxiter-10))) then	
			call intpr("iter", -1, iter, 1) 
			call intpr("nclu1", -1, nclu1, k1) 
			call intpr("clu1", -1, clu1, n1) 	
			call intpr("nclu2", -1, nclu2, k2) 
			call intpr("clu2", -1, clu2, n2) 			
			call intpr("nA", -1, nA, k1*k2) 					
			call dblepr("err", -1, err, 1) 
			call dblepr("terr", -1, terr, 1)
			call dblepr("berr", -1, berr, 1)
			call dblepr("sumA", -1, sumA, k1*k2) 
			call dblepr("sumA2", -1, sumA2, k1*k2) 
			call intpr("nA", -1, nA, k1*k2) 
		end if
	enddo  
	
!	call intpr("nAD", -1, nAD, k) 
!	call dblepr("sumAD", -1, sumAD, k) 
!	call dblepr("sumAD2", -1, sumAD2, k) 
!	call intpr("nA", -1, nA, k*k) 
!	call dblepr("sumA", -1, sumA, k*k) 
!	call dblepr("sumA2", -1, sumA2, k*k) 
end


