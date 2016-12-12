c
c
c     Main program
c
c

	subroutine bootMedian(n,m,times,sfailed,ntimes,cen,tau,
     *	                      nSimBoot,plan,percentil,median)


c       Percentil que buscamos. Para no cambiar los programas llamo median al
c                               resultado aunque le paso el parametro percentil


	 implicit none
         integer n,ntimes,maxtimes
	 double precision nu,theta
	 
	 parameter (maxtimes=800)	 
	 
	 integer m(n),i,j,count
	 double precision tau(n),control,Tacum,timesRan,rexp
	 double precision times(ntimes),cen(n)

	 double precision timesPSH(maxtimes),sfailed(ntimes) 

         integer mBoot(n),ntimesBoot,numdistinct,numdistinctFn
	 double precision timesBoot(maxtimes), sfailedBoot(maxtimes)
	 double precision cenBoot(n),distinct(maxtimes)
	 double precision distinctFn(maxtimes),medianPSH,medianWC
	 double precision medianMLE,percentil
	 integer numdeaths(maxtimes),AtRisk(n,maxtimes)
	 
	 double precision FnPSH(maxtimes),FnWC(maxtimes)
	 
	 integer ss,nSimBoot

       
         double precision median(nSimBoot),alphaMLE(nSimBoot)
         integer status,plan	 
     	 
	  
	 double precision rgam,alphaTheor,alpha,lambda0(maxtimes),
     *   lambda(maxtimes),alphaIni,FnMLE(maxtimes),survMLE(maxtimes),
     *   survF0(maxtimes),FnF0(maxtimes),varZ, 
     *   lambdaIni(maxtimes),lambda0Ini(maxtimes)
     
         integer aux,sample2    
     
         	 
         call DistinctFailed(n,m,times,sfailed,ntimes,
     *              cen,numdistinctFn,distinctFn,numdeaths,AtRisk)
         

         call pshFn(n,numdistinctFn,AtRisk,distinctFn,numdeaths,
     *       FnPSH) 
     
     
        call wcFn(n,m,times,ntimes,cen,numdistinctFn,
     *       distinctFn,AtRisk,FnWC)

        
        if (FnWC(numdistinctFn).lt.1.d0) then
          numdistinctFn=numdistinctFn+1
          FnWC(numdistinctFn)=1.d0
          FnPSH(numdistinctFn)=1.d0
          distinctFn(numdistinctFn)=1.d30
        end if
        
        
        
        if ((plan.eq.6).or.(plan.eq.7)) then	   
c
c  MLE Frailty procedure  
c	   
                              
         call mleALL(n,m,numdistinctFn,distinctFn,numdeaths,1,
     *	   AtRisk,lambda0Ini,lambdaIni,alphaIni,survMLE,
     *     percentil,medianMLE,status)  
  
        end if
           

c
c  Inicio de las muestras bootstrap   
c

           
      do ss=1,nSimBoot
       
	 
       if (plan.eq.1) then 
c
c       Boostrap plan I
c         
           call boot1(n,m,tau,ntimes,times,
     *        	  mBoot,ntimesBoot,timesBoot,cenBoot)
	  
	  
	   do i=1,ntimesBoot
            sfailedBoot(i)=timesBoot(i)
	   end do

	   call sorter(sfailedBoot,ntimesBoot)
 
  
         call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)
            
     
         call psh(n,numdistinct,AtRisk,distinct,numdeaths,
     *               percentil,medianPSH) 
         
	    
	   median(ss)=medianPSH
	    
   	 

       elseif (plan.eq.2) then 
c
c       Boostrap plan II nuevo
c         
 
	   call boot6(n,m,tau,numdistinctFn,distinctFn,FnPSH,
     *        	  mBoot,ntimesBoot,timesBoot,cenBoot)
          
           
	   do i=1,ntimesBoot
            sfailedBoot(i)=timesBoot(i)
	   end do

	   call sorter(sfailedBoot,ntimesBoot)

   	   call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)

         call psh(n,numdistinct,AtRisk,distinct,numdeaths,
     *            percentil,medianPSH) 
       
	   median(ss)=medianPSH

        	   
       elseif (plan.eq.3) then
c	 
c       Boostrap plan III nuevo
c         

	   call boot7(n,m,tau,numdistinctFn,distinctFn,FnPSH,
     *        	  mBoot,ntimesBoot,timesBoot,cenBoot)
  

	   do i=1,ntimesBoot
           sfailedBoot(i)=timesBoot(i)
	   end do
       
	   call sorter(sfailedBoot,ntimesBoot)
      
         call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)

         call psh(n,numdistinct,AtRisk,distinct,numdeaths,
     *            percentil,medianPSH) 
	   
    	   median(ss)=medianPSH

       elseif (plan.eq.4) then
c
c       Boostrap plan IV nuevo
c
	 
	   
	   call boot6(n,m,tau,numdistinctFn,distinctFn,FnWC,
     *        	  mBoot,ntimesBoot,timesBoot,cenBoot)
       
	   
	   do i=1,ntimesBoot
           sfailedBoot(i)=timesBoot(i)
	   end do

	   call sorter(sfailedBoot,ntimesBoot)

  	   call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)
     
         
         call wc(n,mBoot,timesBoot,ntimesBoot,cenBoot,numdistinct,
     *	       distinct,AtRisk,percentil,medianWC) 
	   
	   median(ss)=medianWC
        
       	   
       elseif (plan.eq.5) then 
c
c       Boostrap plan V nuevo
c         

	   call boot7(n,m,tau,numdistinctFn,distinctFn,FnWC,
     *        	  mBoot,ntimesBoot,timesBoot,cenBoot)
	  
	   do i=1,ntimesBoot
            sfailedBoot(i)=timesBoot(i)
	   end do
       
	   call sorter(sfailedBoot,ntimesBoot)
      
         call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)

        call wc(n,mBoot,timesBoot,ntimesBoot,cenBoot,numdistinct,
     *	       distinct,AtRisk,percentil,medianWC) 
	   	   
        median(ss)=medianWC       

       elseif (plan.eq.6) then
c
c       Boostrap plan VI nuevo (Semiparametric)
c
	      
	   call boot8(n,m,tau,real(alphaIni),lambdaIni,numdistinctFn,
     *	   distinctFn,mBoot,ntimesBoot,timesBoot,cenBoot,varZ)
         
         	            	   
	   do i=1,ntimesBoot
           sfailedBoot(i)=timesBoot(i)
	   end do

	   call sorter(sfailedBoot,ntimesBoot)

   	   call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)
     
                       
c       En caso que no quiera buscar la semilla	   
c	   alpha=1/varZ   y  0 en mleALL
	   
	   
	   call mleALL(n,mBoot,numdistinct,distinct,numdeaths,1,
     *	   AtRisk,lambda0,lambda,alpha,survMLE,percentil,medianMLE,
     *    status)  
           
         
           median(ss)=medianMLE
	   alphaMLE(ss)=alpha
        
       else
c
c       Boostrap plan VII nuevo (Semiparametric)
c
	   
         call boot9(n,m,tau,real(alphaIni),lambdaIni,numdistinctFn,
     *	   distinctFn,mBoot,ntimesBoot,timesBoot,cenBoot,varZ)
             
         	   
	   do i=1,ntimesBoot
           sfailedBoot(i)=timesBoot(i)
	   end do

	   call sorter(sfailedBoot,ntimesBoot)

   	   call DistinctFailed(n,mBoot,timesBoot,sfailedBoot,ntimesBoot,
     *              cenBoot,numdistinct,distinct,numdeaths,AtRisk)
     
        
c       En caso que no quiera buscar la semilla	   
c	   alpha=1/varZ     y 0 en mleALL
	   
	   
        call mleALL(n,mBoot,numdistinct,distinct,numdeaths,1,
     *	   AtRisk,lambda0,lambda,alpha,survMLE,percentil,medianMLE,
     *    status)  

         
         median(ss)=medianMLE
	 alphaMLE(ss)=alpha
     
       end if 
	
	end do 
      
	
      end subroutine bootMedian

      

c
c
c     Auxiliar functions 
c
c


      function sample(m)
	 
	 implicit none
	 integer sample,m,taux
         double precision aux, rand
c        real aux,rand
	        
c        call DRNUN(1,aux)	  
c        call random_number(aux)
         
c        call system_clock(taux)       
c        aux=rand(taux)
        
         aux=rand(0)
         
        sample=int(aux*m)+1
        return

      end function


       function sample2(n,prob)
	  
	  implicit none
	  integer n,sample2,taux
	  double precision prob(n)
          double precision aux, rand

        
c	call DRNUN(1,aux)
c       call random_number(aux)	  

c       call system_clock(taux)       
c       aux=rand(taux)
        
        aux=rand(0)
	  
        sample2=1
	do while ((prob(sample2).lt.aux).and.(sample2.lt.n))
          sample2=sample2+1
	end do
	 
	return
	 
       end function



c
c   Bootstrap plans
c

 
      subroutine boot1(n,m,tau,ntimes,times,
     *               	  mBoot,ntimesBoot,timesBoot,cenBoot)
       implicit none
       integer n,ntimes
	 integer m(n),i,j,k,l,r,sample,ntimesBoot,pos
	 integer mBoot(n)
	 double precision tau(n),times(ntimes),control
	 double precision timesBoot(800),cenBoot(n)
	 
	 
	 ntimesBoot=0
	
	  do i=1,n
	   
	   if (i.eq.1) then
	    r=1
		do while ((r.le.n).and.(m(r).eq.0))
		   r=r+1
	    end do
	    mBoot(i)=m(r)
	    j=r
	    if (mBoot(i).eq.0) then
           write(*,*) "cagada todos son cero"
	    end if
	   else
	     j=sample(n)
	     mBoot(i)=m(j)
	   end if

	   pos=0
	   do k=1,j-1
	    pos=pos+m(k)
         end do
        
	   control=0.0d0
	   do l=1,mBoot(i)
	    timesBoot(l+ntimesBoot)=times(pos+l)
	    control=control+times(pos+l)
	   end do    

         ntimesBoot=ntimesBoot+mBoot(i) 
	   cenBoot(i)=tau(j)-(control)
	          
	  end do
      
	end subroutine

       subroutine boot6(n,m,tau,ntimes,times,prob,
     *	              mBoot,ntimesBoot,timesBoot,cenBoot)
       implicit none
       integer n,ntimes
	 integer m(n),i,j,count,sample,sample2
	 integer mBoot(n),ntimesBoot
	 double precision tau(n),times(ntimes),prob(ntimes)
	 double precision control,controlOld
	 double precision timesBoot(800),cenBoot(n)
	 
         ntimesBoot=0
	 do i=1,n
        
	  count=0
	  j=sample2(ntimes,prob)
	  control=times(j)

	  do while (control.le.tau(i))
  	   timesBoot(1+count+ntimesBoot)=times(j)
	   count=count+1
	   j=sample2(ntimes,prob)
	   control=control+times(j)
	  end do
	  
	  mBoot(i)=count
	  ntimesBoot=ntimesBoot+count
	  
	  cenBoot(i)=tau(i)-(control-times(j))
	 end do
	
	end subroutine

      subroutine boot7(n,m,tau,ntimes,times,prob,
     *	              mBoot,ntimesBoot,timesBoot,cenBoot)
       implicit none
       integer n,ntimes
	 integer m(n),i,j,l,r,count,sample,sample2
	 integer mBoot(n),ntimesBoot,bucle
	 double precision tau(n),times(ntimes),prob(ntimes),control
	 double precision timesBoot(800),cenBoot(n),tau0
	 
	 
	 ntimesBoot=0
	 bucle=0
	 do while ((ntimesBoot.eq.0).and.(bucle.lt.100))

	  do i=1,n
	   	   
	   l=sample(n)
         tau0=tau(l)

	   count=0
	   j=sample2(ntimes,prob)
	   control=times(j)
	   do while (control.le.tau0)
  	    timesBoot(1+count+ntimesBoot)=times(j)
	    count=count+1
	    j=sample2(ntimes,prob)
	    control=control+times(j)
	   end do

	   mBoot(i)=count
	   ntimesBoot=ntimesBoot+count
         cenBoot(i)=tau0-(control-times(j))
	  end do
	  bucle=bucle+1
	 end do
	
	end subroutine


      subroutine boot8(n,m,tau,alpha,lambda0,numdistinct,distinct,
     *	              mBoot,ntimesBoot,timesBoot,cenBoot,varZ)
       implicit none
       integer n,numdistinct
	 integer m(n),i,j,k,count,sample,sample2
	 integer mBoot(n),ntimesBoot
	 double precision tau(n),control,controlOld
	 double precision timesBoot(800),cenBoot(n)

	 double precision lambda0(numdistinct),distinct(numdistinct)
	 double precision survF0(numdistinct),FnF0(numdistinct)
	 double precision mm,varZ,Z(n),dsum  

	 real random_gamma,rgam,alpha
  
	 ntimesBoot=0
	 do i=1,n
        

c	   call drngam(1,alpha,rgam)
            
         rgam=random_gamma(alpha)

         Z(i)=((1/alpha)*(rgam))
        
	   do k=1,numdistinct
   	      survF0(k)=exp(-Z(i)*dsum(lambda0,numdistinct,k))
	   end do
	   
         
	   do k=1,numdistinct
           FnF0(k)=1-survF0(k)
	   end do
	  
	   	  
	  j=sample2(numdistinct,FnF0)
	  control=distinct(j)
	     
	  count=0
	  

	  do while (control.le.tau(i))
  	   timesBoot(1+count+ntimesBoot)=distinct(j)
	   count=count+1
	   j=sample2(numdistinct,FnF0)
	   control=control+distinct(j)
	  end do
	  
	  mBoot(i)=count
	  ntimesBoot=ntimesBoot+count
	  
	  cenBoot(i)=tau(i)-(control-distinct(j))
	 end do
	 
	   mm=dsum(Z,n,n)/n
         varZ=0
	   do i=1,n
	    varZ=varZ+(Z(i)-mm)**2
	   end do

	end subroutine

      subroutine boot9(n,m,tau,alpha,lambda0,numdistinct,distinct,
     *	              mBoot,ntimesBoot,timesBoot,cenBoot,varZ)
       implicit none
       integer n,numdistinct,bucle,l
	 integer m(n),i,j,k,count,sample,sample2
	 integer mBoot(n),ntimesBoot
	 double precision tau(n),control,controlOld,tau0
	 double precision timesBoot(800),cenBoot(n)

	 double precision lambda0(numdistinct),distinct(numdistinct)
	 double precision survF0(numdistinct),FnF0(numdistinct)
	 double precision mm,varZ,Z(n),dsum

	      
	 real random_gamma,rgam,alpha

       ntimesBoot=0
	 bucle=0
	 do while ((ntimesBoot.eq.0).and.(bucle.lt.100))

	  do i=1,n
	   	   
	   l=sample(n)
         tau0=tau(l)
 

c	   call drngam(1,alpha,rgam)
	   rgam=random_gamma(alpha)
	   
         Z(i)=((1/alpha)*(rgam))
        
	   do k=1,numdistinct
   	      survF0(k)=exp(-Z(i)*dsum(lambda0,numdistinct,k))
	   end do
	   
         
	   do k=1,numdistinct
           FnF0(k)=1-survF0(k)
	   end do
	  
	   	  
	  j=sample2(numdistinct,FnF0)
	  control=distinct(j)
	     
	  count=0
	  

	  do while (control.le.tau0)
  	   timesBoot(1+count+ntimesBoot)=distinct(j)
	   count=count+1
	   j=sample2(numdistinct,FnF0)
	   control=control+distinct(j)
	  end do
	  
	  mBoot(i)=count
	  ntimesBoot=ntimesBoot+count
	  
	  cenBoot(i)=tau0-(control-distinct(j))
	  end do
	  bucle=bucle+1
	 end do

	 
	   mm=dsum(Z,n,n)/n
         varZ=0
	   do i=1,n
	    varZ=varZ+(Z(i)-mm)**2
	   end do

	end subroutine



      subroutine pshFn(n,numdistinct,AtRisk,distinct,numdeaths,Fn)


	 implicit none

	 integer n,numdistinct,i,j,pos
	 integer AtRisk(n,numdistinct),numdeaths(numdistinct)
	 double precision Fn(numdistinct),distinct(numdistinct)
	 double precision AtRiskOk(numdistinct),surv(numdistinct)


	 do i=1,numdistinct
        AtRiskOK(i)=0.d0
	 end do
	 
	 do i=1,numdistinct
	  do j=1,n
	   AtRiskOK(i)=AtRiskOk(i)+AtRisk(j,i)
        end do
       end do

              
       surv(1)=1-(numdeaths(1)/AtRiskOk(1))

	 do i=2,numdistinct
         surv(i)=surv(i-1)*(1-(numdeaths(i)/AtRiskOk(i)))
	 end do
       
       	
       
       
       do i=1,numdistinct
        Fn(i)=1-surv(i)         
       end do
       
       
       

      end subroutine	

      

      subroutine wcFn(n,m,failed,nfailed,censored,numdis,distinct,
     .               AtRisk,Fn)

      implicit none

      integer n,nfailed
      integer m(n),numdis

      double precision failed(nfailed),censored(n)
      double precision distinct(numdis)
      integer vAtRisk(n*numdis)
      integer AtRisk(n,numdis)

      integer i,j,l,jj,mcumold,mcum

      double precision mstar(n)
      double precision dstar(n,numdis),rstar(n,numdis),ple(numdis)

      double precision Fn(numdis),dstarcum(800),rstarcum(800)

      
      do i=1,n
	do j=1,numdis
       dstar(i,j)=0
	end do
	end do

      mcum = 0

      do i=1,n
        mcumold = mcum
        mcum = mcum + m(i)

        if(m(i) .eq. 0) then
           mstar(i) = 1.d0
        else
           mstar(i) = dfloat(m(i))
        endif

        do l=1,numdis
         if(m(i).gt.0) then
          do jj=mcumold+1,mcum
           if(failed(jj) .eq.  distinct(l)) then
              dstar(i, l) = dstar(i, l) + 1
           endif
          end do

          if(censored(i) .ge.  distinct(l)) then
              rstar(i, l) = dfloat(AtRisk(i, l)) - 1.d0
          else
              rstar(i,l) = dfloat(AtRisk(i,l))
          endif
         else
           rstar(i, l) = dfloat(AtRisk(i, l))
         endif
        end do
      end do

      do 40 j=1,numdis
		dstarcum(j)=0.d0
		rstarcum(j)=0.d0
	do 45 i=1,n
		dstarcum(j)=dstarcum(j) + dstar(i,j)/(mstar(i))
		rstarcum(j)=rstarcum(j) + rstar(i,j)/(mstar(i))
45	continue
40	continue
c
		
	ple(1) = 1.d0 - dstarcum(1)/rstarcum(1)
	do 50 l=2,numdis
	   ple(l)=ple(l-1)*(1.d0-dstarcum(l)/rstarcum(l))
50	continue

       do i=1,numdis
        Fn(i)=1.d0-ple(i)
	 end do
     
      end subroutine




       subroutine psh(n,numdistinct,AtRisk,distinct,numdeaths,
     *            percentil,medianPSH) 


	 implicit none

	 integer n,numdistinct,i,j,pos
	 integer AtRisk(n,numdistinct),numdeaths(numdistinct)
	 double precision survfuncPSH(numdistinct),distinct(numdistinct)
	 double precision AtRiskOk(numdistinct),medianPSH,percentil

       do i=1,numdistinct
	  AtRiskOK(i)=0.d0
       end do

	 do i=1,numdistinct
	  do j=1,n
	   AtRiskOK(i)=AtRiskOK(i)+AtRisk(j,i)
        end do 
       end do

       survfuncPSH(1)=1-(numdeaths(1)/AtRiskOk(1))

	 do i=2,numdistinct
         survfuncPSH(i)=survfuncPSH(i-1)*
     .	              (1-(numdeaths(i)/AtRiskOk(i)))
	 end do
	
         	 	 
	 if (survfuncPSH(numdistinct).le.percentil) then
	  pos=1
	  do while (survfuncPSH(pos).gt.percentil)
           pos=pos+1
	  end do
	  medianPSH=distinct(pos)
	 else
	  medianPSH=-1.0d0
	 end if 
      
	end subroutine	
         


      subroutine wc(n,m,failed,nfailed,censored,numdis,distinct,
     .               AtRisk,percentil,medianWC)

      implicit none

      integer n,nfailed
      integer m(n),numdis

      double precision failed(nfailed),censored(n)
      double precision distinct(numdis)
      integer AtRisk(n,numdis)

      integer i,j,l,jj,mcumold,mcum,pos

      double precision mstar(n)
      double precision dstar(n,numdis),rstar(n,numdis),ple(numdis)

      double precision dstarcum(800),rstarcum(800),medianWC,percentil

      do i=1,n
	do j=1,numdis
       dstar(i,j)=0.0d0
	end do
	end do

      mcum = 0

      do i=1,n
        mcumold = mcum
        mcum = mcum + m(i)

        if(m(i) .eq. 0) then
           mstar(i) = 1.d0
        else
           mstar(i) = dfloat(m(i))
        endif

        do l=1,numdis
         if(m(i).gt.0) then
          do jj=mcumold+1,mcum
           if(failed(jj) .eq.  distinct(l)) then
              dstar(i, l) = dstar(i, l) + 1
           endif
          end do

          if(censored(i) .ge.  distinct(l)) then
              rstar(i, l) = dfloat(AtRisk(i, l)) - 1.d0
          else
              rstar(i,l) = dfloat(AtRisk(i,l))
          endif
         else
	       rstar(i, l) = dfloat(AtRisk(i, l))
         endif
        end do
      end do

      do 40 j=1,numdis
		dstarcum(j)=0.d0
		rstarcum(j)=0.d0
	do 45 i=1,n
		dstarcum(j)=dstarcum(j) + dstar(i,j)/mstar(i)
		rstarcum(j)=rstarcum(j) + rstar(i,j)/mstar(i)
	
      
45    continue
40	continue
c

	ple(1) = 1.0d0 - (dstarcum(1)/rstarcum(1))
	do 50 l=2,numdis
	   ple(l)=ple(l-1)*(1.0d0-(dstarcum(l)/rstarcum(l)))
50	continue

	 
	 if (ple(numdis).le.percentil) then
	  pos=1
	  do while (ple(pos).gt.percentil)
           pos=pos+1
	  end do
      	  medianWC=distinct(pos)
      	 else
      	  medianWC=-1.0d0 
      	 end if 

       end subroutine


      subroutine mleALL(n,m,numdistinct,distinct,numdeaths,searchProc,
     *	   AtRisk,lambda0,lambda,alpha,survMLE,percentil,medianMLE,
     *    status)  

       implicit none
	 
	 integer n,i,j,pos,searchProc
	 integer m(n),numdistinct
	 integer numdeaths(numdistinct),AtRisk(n,numdistinct)
	 double precision distinct(numdistinct),medianMLE

	 double precision alpha_min,alpha_max,alpha,lambda0(numdistinct),
     *   lambda(numdistinct),survMLE(numdistinct),
     *   alphadel,alphaseeds(7),AtRiskOK(numdistinct),percentil
	 integer IER,status,ind
      


         do i=1,numdistinct
	   AtRiskOK(i)=0.d0
         end do

	 do i=1,numdistinct
	  do j=1,n
	   AtRiskOK(i)=AtRiskOK(i)+AtRisk(j,i)
          end do 
         end do


         do i=1,numdistinct
          lambda(i)=float(numdeaths(i))/AtRiskOK(i)
	    if (lambda(i).eq.0.0d0) then
             lambda(i)=0.000001d0
	     end if
	 end do
	   
	 do i=1,numdistinct
	   lambda0(i)=lambda(i)
         end do
         
	 alpha_min=0.5d0
	 alpha_max=distinct(numdistinct)

       
	 if (searchProc.eq.1) then
	    call SearchForSeed(n,m,numdistinct,distinct,numdeaths,
     *	      AtRisk,lambda,alpha_min,alpha_max,0.000001,alpha,IER)
	 end if 
	   
	   
	 alphadel=alpha/4

         alphaseeds(1)=alpha
	 alphaseeds(2)=alpha-alphadel
	 alphaseeds(3)=alpha-2*alphadel
	 alphaseeds(4)=alpha-3*alphadel
	 alphaseeds(5)=alpha+alphadel
	 alphaseeds(6)=alpha+2*alphadel
	 alphaseeds(7)=alpha+3*alphadel
	   
	 status=0
	 ind=0
        
	 do while ((status.eq.0).and.(ind.lt.7))
	   ind=ind + 1
           alpha=alphaseeds(ind)
         
           call emalgo(n,m,numdistinct,distinct,numdeaths,
     *		AtRisk,lambda,alpha,0.00001d0,500,status)	   
	
	   
	 end do 
         

	 call mlevalue(numdistinct,alpha,lambda,survMLE)

	 
	 if (survMLE(numdistinct).le.percentil) then
	  pos=1
	  do while (survMLE(pos).gt.percentil)
           pos=pos+1
	  end do
	  medianMLE=distinct(pos)
	 else
	  medianMLE=-1.0d0
	 end if 
	 
      
      return 	
      end subroutine 





      double precision function dsum(vec,n,k)
        implicit none 
	  integer n,k,i
	  double precision vec(n)
	  dsum=0.d0
	  do i=1,k 
         dsum=dsum+vec(i)
	  end do
        return
	end function
      
      
	
	INTEGER FUNCTION atpos(v, n, x)

        !     Determines the position of the scalar x, within the vector v

                integer n
                DOUBLE PRECISION v(n), x

                INTEGER i

                IF (x .LT. v(1)) THEN
                   atpos = 0
                   RETURN
                END IF

                IF (x .GT. v(n)) THEN
                   atpos = n
                   RETURN
                END IF

                DO 100, i=1, n
                   IF (x-v(i) .GE. 0.0d0) THEN
                          atpos = i
                END IF

 100      CONTINUE

                RETURN
      END


      SUBROUTINE sorter(v, n)

C     Sort the vector v, increasingly

C     length of v
      INTEGER n

C     vector to be sorted
      DOUBLE PRECISION v(n)

      LOGICAL qdone
      INTEGER i
      DOUBLE PRECISION temp

      IF (n .EQ. 1) RETURN

 100  qdone = .TRUE.

      DO 200, i=1, n-1

         IF (v(i) .GT. v(i+1)) THEN
            qdone = .FALSE.
            temp = v(i)
            v(i) = v(i+1)
            v(i+1) = temp
         END IF

 200  CONTINUE

      IF (.NOT. qdone) GO TO 100

      RETURN
      END

!    *************************************

!          Random number functions  

!    *************************************


	FUNCTION random_normal() RESULT(fn_val)

c	! Adapted from the following Fortran 77 code
c	!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
c	!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
c	!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

c	!  The function random_normal() returns a normally distributed pseudo-random
c	!  number with zero mean and unit variance.

c	!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
c	!  and J.F. Monahan augmented with quadratic bounding curves.

	IMPLICIT NONE

        REAL   :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0   
c     *                          vsmall = TINY(1.0), vlarge = HUGE(1.0)
c        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

      integer taux 	
	REAL :: fn_val,rand

c	!     Local variables
	REAL :: s=0.449871, t=-0.386595, a = 0.1960, b = 0.25472
       REAL :: r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

c	!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

	DO

c	  CALL RANDOM_NUMBER(u)
c	  CALL RANDOM_NUMBER(v)

c        call system_clock(taux)       
c        u=rand(taux)
         u=rand(0)
         
c        call system_clock(taux)
c        v=rand(taux)
         v=rand(0)

	  v = 1.7156 * (v - half)

c	!     Evaluate the quadratic form
	  x = u - s
	  y = ABS(v) - t
	  q = x**2 + y*(a*y - b*x)

c	!     Accept P if inside inner ellipse
	  IF (q < r1) EXIT
c	!     Reject P if outside outer ellipse
	  IF (q > r2) CYCLE
c	!     Reject P if outside acceptance region
	  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
	END DO

c	!     Return ratio of P's coordinates as the normal deviate
	fn_val = v/u
	RETURN

	END FUNCTION random_normal



c	!FUNCTION random_gamma(s,first) RESULT(fn_val)
	FUNCTION random_gamma(s) RESULT(fn_val)

c	! Adapted from Fortran 77 code from the book:
c	!     Dagpunar, J. 'Principles of random variate generation'
c	!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

c	!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
c	!     CALLS EITHER random_gamma1 (S > 1.0)
c	!     OR random_exponential (S = 1.0)
c	!     OR random_gamma2 (S < 1.0).

c	!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

      IMPLICIT NONE
        REAL:: zero = 0.0, half = 0.5, one = 1.0, two = 2.0   
c     *                          vsmall = TINY(1.0), vlarge = HUGE(1.0)
c        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
	REAL random_gamma1,random_gamma2,random_exponential,rand

c                REAL, INTENT(IN)    :: s
                 REAL::s
c		!LOGICAL, INTENT(IN) :: first
		REAL                :: fn_val

		IF (s <= zero) THEN
		  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
		  STOP
		END IF

		IF (s > one) THEN
c		!  fn_val = random_gamma1(s, first)
		  fn_val = random_gamma1(s)
		ELSE IF (s < one) THEN
c		!  fn_val = random_gamma2(s, first)
		  fn_val = random_gamma2(s)
		ELSE
		  fn_val = random_exponential()
		END IF

		RETURN
		END FUNCTION random_gamma



c	!FUNCTION random_gamma1(s, first) RESULT(fn_val)
	FUNCTION random_gamma1(s) RESULT(fn_val)
c	! Uses the algorithm in
c	! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
c	! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

c	! Generates a random gamma deviate for shape parameter s >= 1.

	IMPLICIT NONE
        REAL  :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0   
c     *                          vsmall = TINY(1.0), vlarge = HUGE(1.0)
c        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
         REAL random_normal 
      
c        REAL, INTENT(IN)    :: s
         REAL:: s
c	!LOGICAL, INTENT(IN) :: first
	REAL                :: fn_val

c	! Local variables
c        REAL, SAVE  :: c, d
        REAL  :: c, d
	REAL        :: u, v, x,rand
      integer taux

c	!IF (first) THEN
	  d = s - one/3.
	  c = one/SQRT(9.0*d)
c	!END IF

c	! Start of main loop
	DO

c	! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

	  DO
	    x = random_normal()
	    v = (one + c*x)**3
	    IF (v > zero) EXIT
	  END DO

c	! Generate uniform variable U

c       CALL RANDOM_NUMBER(u)
      
c        call system_clock(taux)       
c        u=rand(taux)
         u=rand(0)
	  
        IF (u < one - 0.0331*x**4) THEN
	    fn_val = d*v
	    EXIT
	  ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
	    fn_val = d*v
	    EXIT
	  END IF
	END DO

	RETURN
	END FUNCTION random_gamma1



	!FUNCTION random_gamma2(s, first) RESULT(fn_val)
	FUNCTION random_gamma2(s) RESULT(fn_val)

c	! Adapted from Fortran 77 code from the book:
c	!     Dagpunar, J. 'Principles of random variate generation'
c	!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

c	! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
c	! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
c	! GAMMA2**(S-1) * EXP(-GAMMA2),
c	! USING A SWITCHING METHOD.

c	!    S = SHAPE PARAMETER OF DISTRIBUTION
c	!          (REAL < 1.0)

	IMPLICIT NONE
        REAL :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0   
c     *                          vsmall = TINY(1.0), vlarge = HUGE(1.0)
        REAL:: vsmall=2E-34

c        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

	
c        REAL, INTENT(IN)    :: s
         REAL::s

c	!LOGICAL, INTENT(IN) :: first
	REAL                :: fn_val

c	!     Local variables
	REAL       :: r, x, w
c        REAL, SAVE :: a, p, c, uf, vr, d
        REAL :: a, p, c, uf, vr, d,rand
      integer taux  

	IF (s <= zero .OR. s >= one) THEN
	  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
	  STOP
	END IF

c	!IF (first) THEN                        ! Initialization, if necessary
	  a = one - s
	  p = a/(a + s*EXP(-a))
	  IF (s < vsmall) THEN
	    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
	    STOP
	  END IF
	  c = one/s
	  uf = p*(vsmall/a)**s
	  vr = one - vsmall
	  d = a*LOG(a)
c	!END IF

	DO

c        CALL RANDOM_NUMBER(r)
c        call system_clock(taux)       
c        r=rand(taux)
         r=rand(0)
	  
        IF (r >= vr) THEN
	    CYCLE
	  ELSE IF (r > p) THEN
	    x = a - LOG((one - r)/(one - p))
	    w = a*LOG(x)-d
	  ELSE IF (r > uf) THEN
	    x = a*(r/p)**c
	    w = x
	  ELSE
	    fn_val = zero
	    RETURN
	  END IF

c	  CALL RANDOM_NUMBER(r)
c         call system_clock(taux)
c         r=rand(taux)
          r=rand(0) 

	  IF (one-r <= w .AND. r > zero) THEN
	    IF (r*(w + one) >= one) CYCLE
	    IF (-LOG(r) <= w) CYCLE
	  END IF
	  EXIT
	END DO

	fn_val = x
	RETURN

	END FUNCTION random_gamma2




	FUNCTION random_exponential() RESULT(fn_val)

c	! Adapted from Fortran 77 code from the book:
c	!     Dagpunar, J. 'Principles of random variate generation'
c	!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

c	! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
c	! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
c	! TO EXP(-random_exponential), USING INVERSION.

	
        REAL :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0   
c     *              vsmall = TINY(1.0), vlarge = HUGE(1.0)
	REAL  :: fn_val

c	!     Local variable
	REAL  :: r
      integer taux       

 	DO

c        CALL RANDOM_NUMBER(r)
c        call system_clock(taux)       
c        r=rand(taux)
         r=rand(0)

	  IF (r > zero) EXIT
	END DO

	fn_val = -LOG(r)
	RETURN

	END FUNCTION random_exponential


