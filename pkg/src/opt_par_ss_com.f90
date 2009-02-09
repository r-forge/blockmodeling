subroutine optparsscom(M,clu,diag,maxiter,n,k,err,E,BM,cluM,nbest,iter,printIter)
	INTEGER iter, maxiter, n, clu, nclu, bnclu, iclu, k, i, j, ii, jj, cluM, nbest, nA, nAD, tnA, tnAD, bnA, bnAD, bclu, oldiclu
	DOUBLE PRECISION M, E, BM, err, mean, sumA, sumAD, sumA2, sumAD2, tsumA, tsumAD, tsumA2, tsumAD2, tE, tBM, terr, bsumA, bsumAD
	DOUBLE PRECISION bsumA2, bsumAD2, bE, bBM, berr
	LOGICAL diag, imp, printIter
	DIMENSION M(n,n), clu(n), E(k,k), BM(k,k), cluM(50,n), sumA(k,k), sumAD(k), sumA2(k,k), sumAD2(k), nA(k,k), nAD(k)
	DIMENSION tE(k,k), tBM(k,k),  tsumA(k,k), tsumAD(k), tsumA2(k,k), tsumAD2(k), tnA(k,k), tnAD(k), bE(k,k), bBM(k,k), bclu(n)
	DIMENSION bsumA(k,k), bsumAD(k), bsumA2(k,k), bsumAD2(k), bnA(k,k), bnAD(k)		
	DIMENSION nclu(k), bnclu(k)

	do i = 1, k
		nclu(i) = 0
		nAD(i) = 0
		sumAD(i) = 0
		sumAD2(i) = 0
		do j = 1, k
			 nA(i,j) = 0
			 sumA(i,j) = 0
			 sumA2(i,j) = 0
		end do
	end do

	do i = 1, n
		nclu(clu(i)) = nclu(clu(i)) + 1
		do j = 1, n
			if((i.ne.j) .or. (.not.diag)) then
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
					E(i,j) = 0
					mean = sumAD(i)/nAD(i)
				else
					mean = sumA(i,j)/nA(i,j)
					E(i,j) = sumA2(i,j)-nA(i,j)*mean**2 + sumAD2(i)-sumAD(i)**2/nAD(i)
				end if
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
	
	nbest = 1
	berr=err
	do i = 1, n
		bclu(i)=clu(i)
		cluM(nbest,i) = clu(i)
	end do
	
	do i = 1, k
		bnclu(i) = nclu(i)
		bnAD(i) = nAD(i)
		bsumAD(i) = sumAD(i)
		bsumAD2(i) = sumAD2(i)
		
		do j = 1, k
			bnA(i,j) = nA(i,j)
			bsumA(i,j) = sumA(i,j)
			bsumA2(i,j) = sumA2(i,j)
			bE(i,j)=E(i,j)
			bBM(i,j)=BM(i,j)
		end do
	end do
	
	imp = .TRUE.
	iter = 0
	
	do while (imp .AND. (iter .LE. maxiter))	
		imp = .FALSE.
		iter = iter + 1
		
		do i = 1, n
			oldiclu=clu(i)
			if(nclu(oldiclu).gt.1) then
				do iclu = 1, k
					if(oldiclu.ne.iclu) then
						newiclu = iclu

						do ii = 1, k
							tnAD(ii) = nAD(ii)
							tsumAD(ii) = sumAD(ii)
							tsumAD2(ii) = sumAD2(ii)

							do jj = 1, k
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do

						do j = 1 ,n
							if(i.ne.j) then
								tnA(oldiclu,clu(j)) = tnA(oldiclu,clu(j)) - 1
								tnA(newiclu,clu(j)) = tnA(newiclu,clu(j)) + 1
								tsumA(oldiclu,clu(j)) = tsumA(oldiclu,clu(j)) - M(i,j)
								tsumA(newiclu,clu(j)) = tsumA(newiclu,clu(j)) + M(i,j)
								tsumA2(oldiclu,clu(j)) = tsumA2(oldiclu,clu(j)) - M(i,j)**2
								tsumA2(newiclu,clu(j)) = tsumA2(newiclu,clu(j)) + M(i,j)**2

								tnA(clu(j),oldiclu) = tnA(clu(j),oldiclu) - 1
								tnA(clu(j),newiclu) = tnA(clu(j),newiclu) + 1
								tsumA(clu(j),oldiclu) = tsumA(clu(j),oldiclu) - M(j,i)
								tsumA(clu(j),newiclu) = tsumA(clu(j),newiclu) + M(j,i)
								tsumA2(clu(j),oldiclu) = tsumA2(clu(j),oldiclu) - M(j,i)**2
								tsumA2(clu(j),newiclu) = tsumA2(clu(j),newiclu) + M(j,i)**2

							endif
						end do

						tnAD(oldiclu) = tnAD(oldiclu) - 1
						tnAD(newiclu) = tnAD(newiclu) + 1
						tsumAD(oldiclu) = tsumAD(oldiclu) - M(i,i)
						tsumAD(newiclu) = tsumAD(newiclu) + M(i,i)
						tsumAD2(oldiclu) = tsumAD2(oldiclu) - M(i,i)**2
						tsumAD2(newiclu) = tsumAD2(newiclu) + M(i,i)**2


						terr = 0.0

						do ii = 1, k
							do jj = 1, k
								if(diag.and.ii.eq.jj) then
									if(ii.eq.newiclu .or. ii.eq.oldiclu) then
										if(tnA(ii,jj).eq.0)then
											mean = tsumAD(ii)/tnAD(ii)
											tE(ii,jj) = 0
										else
											mean = tsumA(ii,jj)/tnA(ii,jj)
											tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2 + tsumAD2(ii)-tsumAD(ii)**2/tnAD(ii)
										end if
										tBM(ii,jj) = mean
									else
										tE(ii,jj) = E(ii,jj)
										tBM(ii,jj) = BM(ii,jj)
									end if
									terr = terr + tE(ii,jj)
								else 
									if(ii.eq.newiclu .or. ii.eq.oldiclu .or. jj.eq.newiclu .or. jj.eq.oldiclu) then
										mean = tsumA(ii,jj)/tnA(ii,jj)
										tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
										tBM(ii,jj) = mean
									else
										tE(ii,jj) = E(ii,jj)
										tBM(ii,jj) = BM(ii,jj)
									end if

									terr = terr + tE(ii,jj)
								endif
							end do
						end do

						if (terr.lt.berr) then
!							call intpr("Move", -1, 1, 1) 
							imp = .TRUE.
							berr=terr
							do ii = 1, n
								bclu(ii)=clu(ii)
							end do

							do ii = 1, k
								bnclu(ii) = 0
								bnAD(ii) = tnAD(ii)
								bsumAD(ii) = tsumAD(ii)
								bsumAD2(ii) = tsumAD2(ii)

								do jj = 1, k
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj)=tE(ii,jj)
									bBM(ii,jj)=tBM(ii,jj)
								end do
							end do

							bclu(i) = newiclu

							nbest = 1
							do ii = 1, n
								bnclu(bclu(ii)) = bnclu(bclu(ii)) + 1
								cluM(nbest,ii) = bclu(ii)
							end do
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do ii = 1, n
										cluM(nbest,ii) = clu(ii)
									end do
									cluM(nbest,i) = newiclu
								endif
							end if
						end if

					end if
				end do
			end if

			do j = 1, i - 1
				if(i.ne.j) then
					if(clu(i).ne.clu(j)) then
						!tclu = clu
						!tclu(i) = clu(j)
						!tclu(j) = clu(i)
						newiclu = clu(j)

						do ii = 1, k
							tnAD(ii) = nAD(ii)
							tsumAD(ii) = sumAD(ii)
							tsumAD2(ii) = sumAD2(ii)
							do jj = 1, k
								tnA(ii,jj) = nA(ii,jj)
								tsumA(ii,jj) = sumA(ii,jj)
								tsumA2(ii,jj) = sumA2(ii,jj)
							end do
						end do

						do jj = 1, n
							if(jj.ne.i .and. jj.ne.j) then
								tsumA(oldiclu,clu(jj)) = tsumA(oldiclu,clu(jj)) - M(i,jj) + M(j,jj)
								tsumA(newiclu,clu(jj)) = tsumA(newiclu,clu(jj)) + M(i,jj) - M(j,jj)
								tsumA2(oldiclu,clu(jj)) = tsumA2(oldiclu,clu(jj)) - M(i,jj)**2 + M(j,jj)**2
								tsumA2(newiclu,clu(jj)) = tsumA2(newiclu,clu(jj)) + M(i,jj)**2 - M(j,jj)**2

								tsumA(clu(jj),oldiclu) = tsumA(clu(jj),oldiclu) - M(jj,i) + M(jj,j)
								tsumA(clu(jj),newiclu) = tsumA(clu(jj),newiclu) + M(jj,i) - M(jj,j)
								tsumA2(clu(jj),oldiclu) = tsumA2(clu(jj),oldiclu) - M(jj,i)**2 + M(jj,j)**2
								tsumA2(clu(jj),newiclu) = tsumA2(clu(jj),newiclu) + M(jj,i)**2 - M(jj,j)**2
							endif
						end do

						tsumA(oldiclu,newiclu) = tsumA(oldiclu,newiclu) - M(i,j) + M(j,i)
						tsumA(newiclu,oldiclu) = tsumA(newiclu,oldiclu) + M(i,j) - M(j,i)
						tsumA2(oldiclu,newiclu) = tsumA2(oldiclu,newiclu) - M(i,j)**2 + M(j,i)**2
						tsumA2(newiclu,oldiclu) = tsumA2(newiclu,oldiclu) + M(i,j)**2 - M(j,i)**2

						tsumAD(oldiclu) = tsumAD(oldiclu) - M(i,i) + M(j,j)
						tsumAD(newiclu) = tsumAD(newiclu) + M(i,i) - M(j,j)
						tsumAD2(oldiclu) = tsumAD2(oldiclu) - M(i,i)**2 + M(j,j)**2
						tsumAD2(newiclu) = tsumAD2(newiclu) + M(i,i)**2 - M(j,j)**2

						terr = 0.0

						do ii = 1, k
							do jj = 1, k
								if(diag .and. (ii.eq.jj)) then
									if((ii.eq.newiclu) .or. (ii.eq.oldiclu)) then
										if(tnA(ii,jj).eq.0)then
											mean = tsumAD(ii)/tnAD(ii)
											tE(ii,jj) = 0
										else
											mean = tsumA(ii,jj)/tnA(ii,jj)
											tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2 + tsumAD2(ii)-tsumAD(ii)**2/tnAD(ii)
										end if
										tBM(ii,jj) = mean
									else
										tE(ii,jj) = E(ii,jj)
										tBM(ii,jj) = BM(ii,jj)
									end if
									terr = terr + tE(ii,jj)
								else 
									if(ii.eq.newiclu .or. ii.eq.oldiclu .or. jj.eq.newiclu .or. jj.eq.oldiclu) then
										mean = tsumA(ii,jj)/tnA(ii,jj)
										tE(ii,jj) = tsumA2(ii,jj)-tnA(ii,jj)*mean**2
										tBM(ii,jj) = mean
									else
										tE(ii,jj) = E(ii,jj)
										tBM(ii,jj) = BM(ii,jj)
									end if

									terr = terr + tE(ii,jj)
								endif
							end do
						end do

						if(terr.lt.berr)then
!							call intpr("Switch", -1, 2, 1) 
							imp = .TRUE.
							berr=terr
							do ii = 1, n
								bclu(ii)=clu(ii)
							end do
							bclu(i) = newiclu
							bclu(j) = oldiclu

							do ii = 1, k
								bnclu(ii) = 0
								bnAD(ii) = tnAD(ii)
								bsumAD(ii) = tsumAD(ii)
								bsumAD2(ii) = tsumAD2(ii)

								do jj = 1, k
									bnA(ii,jj) = tnA(ii,jj)
									bsumA(ii,jj) = tsumA(ii,jj)
									bsumA2(ii,jj) = tsumA2(ii,jj)
									bE(ii,jj) = tE(ii,jj)
									bBM(ii,jj) = tBM(ii,jj)
								end do
							end do


							nbest = 1
							do ii = 1, n
								bnclu(bclu(ii)) = bnclu(bclu(ii)) + 1
								cluM(nbest,ii) = bclu(ii)
							end do
						else
							if(terr.eq.berr)then
								nbest = nbest + 1
								if(nbest.le.50) then
									do ii = 1, n
										cluM(nbest,ii) = clu(ii)
									end do
									cluM(nbest,i) = newiclu
									cluM(nbest,j) = oldiclu
								end if
							end if
						end if
					end if			
				end if
			end do
		end do	


		err=berr
		do i = 1, n
			clu(i)=bclu(i)
		end do

		do i = 1, k
			nclu(i) = bnclu(i)
			nAD(i) = bnAD(i)
			sumAD(i) = bsumAD(i)
			sumAD2(i) = bsumAD2(i)

			do j = 1, k
				nA(i,j) = bnA(i,j)
				sumA(i,j) = bsumA(i,j)
				sumA2(i,j) = bsumA2(i,j)
				E(i,j) = bE(i,j)
				BM(i,j) = bBM(i,j)
			end do
		end do	
	if (printIter) then	
		call intpr("iter", -1, iter, 1) 
		call intpr("nclu", -1, nclu, k) 
		call intpr("clu", -1, clu, n) 	
		call dblepr("err", -1, err, 1) 
	end if
	enddo  
	
!	call intpr("nAD", -1, nAD, k) 
!	call dblepr("sumAD", -1, sumAD, k) 
!	call dblepr("sumAD2", -1, sumAD2, k) 
!	call intpr("nA", -1, nA, k*k) 
!	call dblepr("sumA", -1, sumA, k*k) 
!	call dblepr("sumA2", -1, sumA2, k*k) 
end


