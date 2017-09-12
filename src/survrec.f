      integer function countge(vec,nvec,val)

        implicit none

        integer nvec,count,i
        double precision vec(nvec),val

         count=0

         do i=1,nvec
           if (vec(i).ge.val) then
            count=count+1
           endif
         end do

         countge=count

      end function countge



      subroutine distinctfailed(n,m,failed,sfailed,nfailed,
     .              censored,numdistinct,distinct,numdeaths,AtRisk)

c
c      failed times for each subject
c      sfailed failed times ordered
c
      implicit none

      integer n,nfailed
      integer m(n)

      double precision censored(n)
      double precision failed(nfailed),sfailed(nfailed)
      double precision tempfailed(nfailed)
      integer ntempfailed

      integer numdeaths(nfailed)
      double precision distinct(nfailed),current
      integer AtRisk(n,nfailed)
      integer numdistinct
      integer countge

      integer i,j,k,contdeaths,mcum,mcumold


      numdistinct = 1
      contdeaths = 1
      current= sfailed(1)


        do i=2,nfailed
          if (sfailed(i).ne.current) then
               distinct(numdistinct) = current
               numdeaths(numdistinct) = contdeaths
               current = sfailed(i)
               contdeaths = 1
               numdistinct = numdistinct + 1
          else
               contdeaths = contdeaths + 1
          end if
               distinct(numdistinct)=current
               numdeaths(numdistinct)=contdeaths
        end do

        mcum = 0

        do i=1,n
          mcumold = mcum + 1
          mcum = mcum + m(i)
          ntempfailed = mcum-mcumold+1

          do k=1,ntempfailed
           tempfailed(k) = failed(mcumold+k-1)
          end do

          do j=1,numdistinct
           if (m(i).ge.0) then
            AtRisk(i, j) = countge(tempfailed,ntempfailed,distinct(j))
           end if
           if (censored(i).ge.distinct(j)) then
            AtRisk(i, j) = AtRisk(i, j) + 1
           end if
          end do
        end do

      end subroutine distinctfailed





      subroutine wcPLE(n,m,failed,nfailed,censored,numdis,distinct,
     .               vAtRisk,dstar,rstar,mstar)

      implicit none

      integer n,nfailed
      integer m(n),numdis

      double precision failed(nfailed),censored(n)
      double precision distinct(numdis)
      integer vAtRisk(n*numdis)
      integer AtRisk(n,numdis)

      integer i,j,l,jj,mcumold,mcum

      double precision mstar(n)
      double precision dstar(n,numdis),rstar(n,numdis)



      do i=1,numdis
        do j=1,n
          AtRisk(j,i)=vAtRisk(j+(n*(i-1)))
         end do
      end do


      mcum = 0

      do i=1,n
        mcumold = mcum
        mcum = mcum + m(i)

        if(m(i).eq.0) then
           mstar(i) = 1
        else
           mstar(i) = m(i)
        endif

        do l=1,numdis
         if(m(i).gt.0) then
          do jj=mcumold+1,mcum
           if(failed(jj) .eq.  distinct(l)) then
              dstar(i, l) = dstar(i, l) + 1
           endif
          end do

          if(censored(i) .ge.  distinct(l)) then
              rstar(i, l) = AtRisk(i, l) - 1
          else
              rstar(i,l) = AtRisk(i,l)
          endif
         else
          rstar(i, l) = AtRisk(i, l)
         endif
        end do
      end do
      
      end subroutine wcPLE
      
      
      subroutine wc2(n,gtime,ctime,count,mc,m,
     *             cen,ucen,nd,udt,tot,gap,event,
     *             r,d,sest,var)
      integer n,mc,nd,cumni,curj,tot,m(n)
      double precision gtime(n,mc),cen(n,mc)
      double precision gap(tot),event(tot)
      double precision ctime(n)
      double precision count(n),ucen(n)
      double precision udt(nd),r(nd),d(nd)
      double precision sest(nd),var(nd)
      double precision hai,fai,phi,w

c     cumni is a counter : m(i-1) where m(0) = 0
c     curj is a counter within a counter

      cumni=0 
      do 11 i=1,n
         do 10 j=1,m(i)
            curj = cumni + j
            gtime(i,j) = gap(curj)
            cen(i,j) = event(curj)   
 10      continue
         cumni=cumni+m(i)
 11   continue 

      do 30 k=1,nd
         r(k)=0.
         d(k)=0.
         do 25 i=1,n
            if (count(i).gt.1.) then
               do 20 j=1,idint(ucen(i))
                  if(gtime(i,j).ge.udt(k)) then
                     r(k)=r(k)+ctime(i)/ucen(i)
                  endif
                  if(gtime(i,j).eq.udt(k)) then
                     d(k)=d(k)+ctime(i)/ucen(i)
                  endif
 20            continue
            else
               if(gtime(i,1).ge.udt(k)) r(k)=r(k)+ctime(i)
               if((gtime(i,1).eq.udt(k)).and.(cen(i,j).gt.0))
     *              d(k)=d(k)+ctime(i)
            endif
 25      continue
 30   continue

c Calculate survivor estimate at time T (i.e. time nd)
      sest(1)=1.-d(1)/r(1)
      do 50 i=2,nd
         sest(i)=sest(i-1)*(1.-d(i)/r(i))
 50   continue

c Calculate standard error
      do 180 k=1,nd
         phi=0.
         do 170 i=1,n
            fai=0.
            w=0.
            do 140 k2=1,k
               hai=0.
               if (count(i).gt.1.) then
                  do 130 j=1,idint(ucen(i))
                     if(gtime(i,j).ge.udt(k2)) then
                        hai=hai+ctime(i)/ucen(i)
                     endif
 130              continue
               elseif (gtime(i,1).ge. udt(k2)) then
                  hai=hai+ctime(i)
               endif
               w=w+hai*d(k2)/(r(k2)*r(k2))
 140        continue
            do 160 j2=1,idint(ucen(i))
               if (gtime(i,j2).lt. udt(k)) then
                  do 150 k3=1,nd
                     if (gtime(i,j2).eq. udt(k3)) then
                        fai=fai+ctime(i)/(ucen(i)*r(k3))
                     endif
 150              continue
               endif
 160        continue
            phi=phi+(w-fai)*(w-fai)
 170     continue
         var(k)=dsqrt(phi)*sest(k)
 180  continue

      return
      end subroutine wc2

     


c
c implements the em algorithm
c


        subroutine emalgo(n,m,numdis,distinct,numdead,AtRisk,lambda
     .                     ,alpha,tol,maxiter,istatus)


        dimension m(n),distinct(numdis),numdead(numdis)
        dimension AtRisk(n,numdis)
        dimension lambda(numdis),Zhat(n),S(n),Zhatold(n)
        dimension lambold(numdis)

        integer n,m,numdis,numdead,AtRisk,maxiter
        integer iter,istatus,istatus2

        double precision distinct,lambda,alpha,tol,Zhat
        double precision S,Zhatold,lambold
        double precision alphaold,distall,dZhat,dlambda,dalpha,stemp

        
        iter = 0
c
        do 5 i=1,n
5       Zhat(i)=0
        istopck=0
c
        do while (istopck.eq.0)
                iter = iter + 1
                do 10 i=1,n
                S(i)=0
                do 15 j=1,numdis
                S(i)=S(i)+dfloat(AtRisk(i,j))*lambda(j)
15              continue
10              continue
c
                do 20 i=1,n
                Zhatold(i)=Zhat(i)
                Zhat(i) = (1.+dfloat(m(i))/alpha)/(1.+S(i)/alpha)
20              continue
c
                do 25 j=1,numdis
                lambold(j) = lambda(j)
25              continue
c
                do 30 j=1,numdis
                stemp=0
                do 35 i=1,n
35              stemp=stemp+dfloat(AtRisk(i,j))*Zhat(i)
                lambda(j)=dfloat(numdead(j))/stemp
30              continue
c
                alphaold = alpha
       
                call estalpha(n,m,numdis,distinct,numdead,AtRisk,lambda,
     1          alpha,tol,maxiter,istatus2)
c
        
                if((istatus2.eq.0)) then
                        istopck=1
                        istatus=0
                endif
                if(iter.gt.maxiter) then
                        istopck=1
                        istatus=-1
                endif
c
c computing the distances of old and new estimates
                dZhat = 0
                do 80 i=1,n
80              dZhat = dZhat + (Zhat(i)-Zhatold(i))**2
                dZhat=dsqrt(dZhat)
c
                dlambda = 0
                do 81 j=1,numdis
81              dlambda = dlambda + (lambda(j)-lambold(j))**2
                dlambda=dsqrt(dlambda)
c
                dalpha1 = dabs(alphaold-alpha)
                dalpha2 = dabs(1/alphaold -1/alpha)
                dalpha = min(dalpha1,dalpha2)
c               if(dalpha.gt.dalpha2) dalpha=dalpha2
c
                distall = max(dZhat,dlambda,dalpha)
c               if(dlambda.gt.distall) distall=dlambda
c               if(dalpha.gt.distall) distall=dalpha
c
                if((distall.le.tol)) then
                        istopck=1
                        istatus=1
                endif
c
        end do
c
       
        return
        end subroutine emalgo




        subroutine estalpha(n,m,numdis,distinct,numdead,AtRisk,
     .                lambda,alpha,tol,maxiter,status)

        dimension m(n),distinct(numdis),numdead(numdis)
        dimension AtRisk(n,numdis),lambda(numdis)
c       dimension S(100),Soalpha(100),moalpha(100)

        integer n,m,numdis,numdead,AtRisk,maxiter,status

c       double precision moalpha,Soalpha,S
        double precision distinct,lambda,tol
        double precision dist,l1, l11,alpha,oldalpha,C1,C2,C3
        double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w
        double precision  fu,fv,fw,fx,x,tol3
        double precision  dabs,dsqrt,d1mach

                
        status=1
c  c is the squared inverse of the golden ratio
        c=0.5d0*(3.0d0-dsqrt(5.0d0))

c  compute starting interval
        ax=max(alpha-50.d0,0.d0)
        bx=alpha+50.d0
        tol=0.0001
c
c  eps is approximately the square root of the relative machine
c  precision.
c

c  10     eps=2.2204460492503131/(10.d0**(16))

 10     eps=d1mach(4)

        tol1=eps+1.0d0
        eps=dsqrt(eps)
c
        a=ax
        b=bx
        v=a+c*(b-a)
        w=v
        x=v
        e=0.0d0

        call loglik(n,m,numdis,distinct,numdead,AtRisk,lambda,x,l1)
                
        fx=-1.d0*l1

        fv=fx
        fw=fx
        tol3=tol/3.0d0
c
c  main loop starts here
c
 20     xm=0.5d0*(a+b)
        tol1=eps*dabs(x)+tol3
        t2=2.0d0*tol1
c
c  check stopping criterion
c
        if (dabs(x-xm).le.(t2-0.5d0*(b-a))) go to 190
        p=0.0d0
        q=0.0d0
        r=0.0d0
        if (dabs(e).le.tol1) go to 50

c       
c  fit parabola
c
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.0d0*(q-r)
        if (q.le.0.0d0) go to 30
        p=-p
        go to 40
 30     q=-q
 40     r=e
        e=d
 50     if ((dabs(p).ge.dabs(0.5d0*q*r)).or.(p.le.q*(a-x))
     .      .or.(p.ge.q*(b-x))) go to 60
c
c  a parabolic-interpolation step
c
        d=p/q
        u=x+d
c
c  f must not be evaluated too close to ax or bx
c
        if (((u-a).ge.t2).and.((b-u).ge.t2)) go to 90
        d=tol1
        if (x.ge.xm) d=-d
        go to 90
c
c  a golden-section step
c
 60     if (x.ge.xm) go to 70
        e=b-x
        go to 80
 70     e=a-x
 80     d=c*e
c
c  f must not be evaluated too close to x
c
 90     if (dabs(d).lt.tol1) go to 100
        u=x+d
        go to 120
 100    if (d.le.0.0d0) go to 110
        u=x+tol1
        go to 120
 110    u=x-tol1

 120    call loglik(n,m,numdis,distinct,numdead,AtRisk,lambda,u,l1)

        fu=-1.d0*l1
c
c  update  a, b, v, w, and x
c
        if (fx.gt.fu) go to 140
        if (u.ge.x) go to 130
        a=u
        go to 140
 130    b=u
 140    if (fu.gt.fx) go to 170
        if (u.ge.x) go to 150
        b=x
        go to 160
 150    a=x
 160    v=w
        fv=fw
        w=x
        fw=fx
        x=u
        fx=fu
        go to 20
 170    if ((fu.gt.fw).and.(w.ne.x)) go to 180
        v=w
        fv=fw
        w=u
        fw=fu
        go to 20
 180    if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) go to 20
        v=u
        fv=fu
        go to 20
c
c  end of main loop
c
 190    fmin=x

c       call dblepr("lower",5,a,1)
c       call dblepr("upper",5,b,1)

        alpha=fmin

        return
        end subroutine estalpha




        subroutine loglik(n,m,numdis,distinct,numdead,AtRisk,
     .                   lambda,alpha,l1)

        dimension m(n),distinct(numdis),numdead(numdis)
        dimension AtRisk(n,numdis),lambda(numdis)
        dimension S(n),Soalpha(n)

c       dimension moalpha(100)
     
        integer n,m,numdis,numdead,AtRisk

c       double precision moalpha
        double precision distinct,lambda,tol,Soalpha,S
        double precision l1, alpha,C1,C2,C3,q1,q2
c
                
        C3 = 0.d0
c
        do 10 i=1,n
        if(m(i) .gt. 0) then
        do 20 j=1,m(i)
           C3 = C3 + dlog(alpha + dfloat(j) - 1.d0)
20      continue
        endif
10      continue

        q2 = 0.d0
        do 25 i=1,n
        S(i) = 0.d0
        do 26 j=1,numdis
        S(i) = S(i) + dfloat(AtRisk(i,j))*lambda(j)
26      continue
        Soalpha(i) = dlog(1.d0+S(i)/alpha)*(alpha+dfloat(m(i)))
        q2=q2+Soalpha(i)
25      continue

        q1 = 0.d0

        do 30 j=1,numdis
           q1=q1+numdead(j)*dlog(lambda(j)/alpha)
30      continue


        l1=q1 - q2 + C3
c        call dblepr("l1",2,l1,1)
                
        end subroutine loglik

          

        subroutine likterms(n,m,alpha,c1,c2,c3)
     
        dimension m(n)
     
        integer n,m

        real alpha,c1,c2,c3
c
        C1 = 0
        C2 = 0
        C3 = 0
c
        do 10 i=1,n
        if(m(i) .gt. 0) then
        do 20 j=1,m(i)
           C1 = C1 + 1/(alpha + j - 1)
           C2 = C2 + 1/(alpha + j - 1)**2
           C3 = C3 + log(alpha + j - 1)
20      continue
        endif
10      continue
c
        return
        end subroutine likterms


        subroutine mlevalue(numdis,alpha,lambda,mle)

        dimension lambda(numdis),mle(numdis)

        integer numdis
        double precision  alpha,lambda,mle
c
        ss=0
        do 10 i=1,numdis
        ss=ss+lambda(i)
        mle(i)=(alpha/(alpha + ss))**alpha
10      continue
c
        return
        end subroutine mlevalue





c
c  SearchForSeed  procedure
c


        subroutine SearchForSeed(n,m,numdis,distinct,numdead,AtRisk
     .               ,lambda,alpha_min,alpha_max,tol,alpha,IER)

c  ojo que falta hacer el buclecillo con alpha.min y alpha.max

        implicit none

        integer n,numdis
        integer m(n),numdead(numdis)
        integer AtRisk(n,numdis)

        double precision distinct(numdis),lambda(numdis)
        double precision alpha_min,alpha_max
        
        integer IP4,IP5
        integer IP3(n+numdis)
        double precision P2(2*numdis),P1(n*numdis)
        
          integer i,j
        
        double precision alpha, mloglik, tol
        integer IER

        external mloglik

        IP4=n
        IP5=numdis

        do i=1,n
          IP3(i)=m(i)
        end do

        do i=1,numdis
          IP3(n+i)=numdead(i)
          P2(i)=distinct(i)
          P2(numdis+i)=lambda(i)
        end do


        do i=1,numdis
         do j=1,n
          P1(j+(n*(i-1)))=AtRisk(j,i)
         end do
        end do


        call ZXGSP(mloglik,P1,P2,IP3,IP4,IP5,
     .                 alpha_min,alpha_max,tol,alpha,IER)

        end subroutine SearchForSeed





        double precision function mloglik(alpha,P1,P2,IP3,IP4,IP5)
c
c  this function is the same that subroutine loglik
c  is adapted for ZXGSP routine
c
        implicit none

        integer IP4,IP5
        integer n,numdis

        integer m(IP4),numdead(IP5)
        integer AtRisk(IP4,IP5)

        double precision distinct(IP5),lambda(IP5)
        double precision alpha,ans


        integer IP3(IP4+IP5)
        double precision P2(2*IP5),P1(IP4*IP5)

        integer i,j

        n=IP4
        numdis=IP5

        do i=1,n
          m(i)=IP3(i)
        end do

        do i=1,numdis
          distinct(i)=P2(i)
          lambda(i)=P2(numdis+i)
          numdead(i)=IP3(n+i)
        end do

        do i=1,numdis
         do j=1,n
          AtRisk(j,i)=P1(j+(n*(i-1)))
         end do
        end do


        call loglik(n,m,numdis,distinct,numdead,AtRisk,
     .              lambda,alpha,ans)


        mloglik=-ans

        return

        end function








CZXGSP
C   IMSL ROUTINE NAME   - ZXGSP                                         ZXGP0010
C                                                                       ZXGP0020
C-----------------------------------------------------------------------ZXGP0030
C                                                                       ZXGP0040
C   COMPUTER            - CDCFT5/SINGLE                                 ZXGP0050
C                                                                       ZXGP0060
C   LATEST REVISION     - JANUARY 1, 1978                               ZXGP0070
C                                                                       ZXGP0080
C   PURPOSE             - ONE-DIMENSIONAL UNIMODAL FUNCTION             ZXGP0090
C                           MINIMIZATION USING THE GOLDEN SECTION       ZXGP0100
C                           SEARCH METHOD - DATA PARAMETERS SPECIFIED   ZXGP0110
C                                                                       ZXGP0120
C   USAGE               - CALL ZXGSP (F,P1,P2,IP3,IP4,IP5,A,B,TOL,XMIN, ZXGP0130
C                           IER)                                        ZXGP0140
C                                                                       ZXGP0150
C   ARGUMENTS    F      - A REAL FUNCTION SUBPROGRAM SUPPLIED BY        ZXGP0160
C                           THE USER. (INPUT)                           ZXGP0170
C                           F MUST BE DECLARED EXTERNAL IN THE CALLING  ZXGP0180
C                           PROGRAM. F DEFINES THE FUNCTION TO BE       ZXGP0190
C                           MINIMIZED AND SHOULD BE OF THE FOLLOWING    ZXGP0200
C                           FORM                                        ZXGP0210
C                             F(X,P1,P2,IP3,IP4,IP5)                    ZXGP0220
C                           WHERE X IS THE INDEPENDENT VARIABLE.        ZXGP0230
C                           P1,P2,IP3,IP4, AND IP5 ARE DESCRIBED BELOW. ZXGP0240
C                           F MUST NOT ALTER X.                         ZXGP0250
C                           F IS ASSUMED TO DEFINE A UNIMODAL FUNCTION. ZXGP0260
C                           THAT IS, A FUNCTION WITH A UNIQUE MINIMUM   ZXGP0270
C                           VALUE, XMIN, IN THE INTERVAL DEFINED BY A   ZXGP0280
C                           AND B. A FUNCTION, F, IS UNIMODAL IF IT     ZXGP0290
C                           SATISFIES THE FOLLOWING CONDITIONS          ZXGP0300
C                             FOR ALL X0, X1, AND X2 IN THE INTERVAL    ZXGP0310
C                             (A,B) INCLUSIVELY, IF X0 IS LESS THAN     ZXGP0320
C                             X1 AND X1 IS LESS THAN X2, THEN           ZXGP0330
C                             1. IF F(X0) IS LESS THAN OR EQUAL TO      ZXGP0340
C                                F(X1), THEN F(X1) IS LESS THAN F(X2).  ZXGP0350
C                             2. IF F(X1) IS GREATER THAN OR EQUAL      ZXGP0360
C                                TO F(X2), THEN F(X0) IS GREATER THAN   ZXGP0370
C                                F(X1).                                 ZXGP0380
C                P1     - PARAMETERS THAT MAY BE USED TO SEND DATA TO   ZXGP0390
C                P2         OR RETURN DATA FROM THE FUNCTION, F.        ZXGP0400
C                IP3        (INPUT/OUTPUT)                              ZXGP0410
C                IP4        P1, P2, AND IP3 ARE VECTORS. IP4 AND IP5    ZXGP0420
C                IP5        ARE SCALERS.                                ZXGP0430
C                A      - (INPUT/OUTPUT)                                ZXGP0440
C                           ON INPUT, A IS THE LOWER ENDPOINT OF THE    ZXGP0450
C                           INTERVAL IN WHICH THE MINIMUM OF F IS TO BE ZXGP0460
C                           LOCATED.                                    ZXGP0470
C                           ON OUTPUT, A IS THE LOWER ENDPOINT OF THE   ZXGP0480
C                           INTERVAL IN WHICH THE MINIMUM OF F IS       ZXGP0490
C                           LOCATED.                                    ZXGP0500
C                B      - (INPUT/OUTPUT)                                ZXGP0510
C                           ON INPUT, B IS THE UPPER ENDPOINT OF THE    ZXGP0520
C                           INTERVAL IN WHICH THE MINIMUM OF F IS TO BE ZXGP0530
C                           LOCATED.                                    ZXGP0540
C                           ON OUTPUT, B IS THE UPPER ENDPOINT OF THE   ZXGP0550
C                           INTERVAL IN WHICH THE MINIMUM OF F IS       ZXGP0560
C                           LOCATED.                                    ZXGP0570
C                TOL    - THE LENGTH OF THE FINAL SUBINTERVAL           ZXGP0580
C                           CONTAINING THE MINIMUM. (INPUT)             ZXGP0590
C                XMIN   - THE APPROXIMATE MINIMUM OF THE FUNCTION F     ZXGP0600
C                           ON THE ORIGINAL INTERVAL (A,B). (OUTPUT)    ZXGP0610
C                           ON OUTPUT, WHEN IER=0, THE FOLLOWING        ZXGP0620
C                           CONDITIONS HOLD                             ZXGP0630
C                           1. (B-A) IS LESS THAN OR EQUAL TO TOL.      ZXGP0640
C                           2. A IS LESS THAN OR EQUAL TO XMIN AND      ZXGP0650
C                              XMIN IS LESS THAN OR EQUAL TO B.         ZXGP0660
C                           3. F(XMIN) IS LESS THAN OR EQUAL TO F(A)    ZXGP0670
C                              AND F(XMIN) IS LESS THAN OR EQUAL TO     ZXGP0680
C                              F(B).                                    ZXGP0690
C                IER    - ERROR PARAMETER. (OUTPUT)                     ZXGP0700
C                         TERMINAL ERROR                                ZXGP0710
C                           IER=129 IMPLIES THAT A IS GREATER THAN      ZXGP0720
C                             OR EQUAL TO B. XMIN IS SET TO A.          ZXGP0730
C                           IER=130 IMPLIES THAT TOL IS GREATER THAN    ZXGP0740
C                             OR EQUAL TO THE LENGTH OF THE INPUT       ZXGP0750
C                             INTERVAL (A,B). XMIN IS SET TO A.         ZXGP0760
C                           IER=131 IMPLIES THAT THE FUNCTION F IS NOT  ZXGP0770
C                             UNIMODAL OR APPEARS TO BE NOT UNIMODAL TO ZXGP0780
C                             THE ROUTINE DUE TO ROUNDING ERRORS IN THE ZXGP0790
C                             EXTERNAL EVALUATION FUNCTION.             ZXGP0800
C                             WHEN THIS ERROR OCCURS THE FOLLOWING      ZXGP0810
C                             CONDITIONS HOLD                           ZXGP0820
C                             1. A IS LESS THAN OR EQUAL TO XMIN AND    ZXGP0830
C                                XMIN IS LESS THAN OR EQUAL TO B.       ZXGP0840
C                             2. F(XMIN) IS GREATER THAN OR EQUAL TO    ZXGP0850
C                                F(A) AND F(XMIN) IS GREATER THAN OR    ZXGP0860
C                                EQUAL TO F(B) (ONLY ONE EQUALITY CAN   ZXGP0870
C                                HOLD).                                 ZXGP0880
C                             FURTHER ANALYSIS OF THE FUNCTION F IS     ZXGP0890
C                             NECESSARY IN ORDER TO DETERMINE WHETHER   ZXGP0900
C                             IT IS NOT UNIMODAL IN THE MATHEMATICAL    ZXGP0910
C                             SENSE OR WHETHER IT APPEARS TO BE NOT     ZXGP0920
C                             UNIMODAL TO THE ROUTINE DUE TO ROUNDING   ZXGP0930
C                             ERRORS IN WHICH CASE THE A,B,AND XMIN     ZXGP0940
C                             RETURNED MAY BE ACCEPTABLE.               ZXGP0950
C                           IER=132 IMPLIES THAT THE INTERVAL HAS BEEN  ZXGP0960
C                             REDUCED AS FAR AS NUMERICALLY POSSIBLE.   ZXGP0970
C                             THIS IS DUE TO TOL BEING TOO SMALL.       ZXGP0980
C                             WHEN THIS ERROR OCCURS, XMIN GIVES        ZXGP0990
C                             THE LOCATION OF THE MINIMUM AS            ZXGP1000
C                             ACCURATELY AS POSSIBLE.                   ZXGP1010
C                                                                       ZXGP1020
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZXGP1030
C                       - SINGLE/H36,H48,H60                            ZXGP1040
C                                                                       ZXGP1050
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 ZXGP1060
C                                                                       ZXGP1070
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZXGP1080
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZXGP1090
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZXGP1100
C                                                                       ZXGP1110
C   REMARKS  1.  THE ARGUMENTS P1,P2,IP3,IP4, AND IP5 ARE NOT           ZXGP1120
C                USED BY THE ROUTINE. THESE ARE PARAMETERS WHICH THE    ZXGP1130
C                USER MAY REQUIRE TO RELAY DATA FROM THE MAIN PROGRAM   ZXGP1140
C                TO THE FUNCTION SUBPROGRAM F. THE CHOICE OF FIVE       ZXGP1150
C                COMMUNICATION PARAMETERS WAS NOT ENTIRELY ARBITRARY.   ZXGP1160
C                FOR EXAMPLE, FIVE PARAMETERS ALLOW THE USE OF          ZXGP1170
C                  A. A REAL 3-DIMENSIONAL ARRAY                        ZXGP1180
C                  B. A REAL VECTOR                                     ZXGP1190
C                  C. AN INTEGER VECTOR                                 ZXGP1200
C                  D. AN INTEGER SCALAR SPECIFYING THE FIRST DIMENSION  ZXGP1210
C                     OF THE ARRAY EXACTLY AS IT APPEARS IN THE MAIN    ZXGP1220
C                     PROGRAM                                           ZXGP1230
C                  E. AN INTEGER SCALAR SPECIFYING THE SECOND DIMENSION ZXGP1240
C                     OF THE ARRAY EXACTLY AS IT APPEARS IN THE MAIN    ZXGP1250
C                     PROGRAM                                           ZXGP1260
C                THIS FEATURE ALLOWS THE USER TO DEFINE A FUNCTION OF   ZXGP1270
C                SOMEWHAT GENERAL APPLICABILITY WITHOUT THE USE OF      ZXGP1280
C                COMMON.                                                ZXGP1290
C            2.  IMSL SUBROUTINE ZXGSN IS IDENTICAL IN                  ZXGP1300
C                PURPOSE TO ZXGSP BUT DOES NOT REQUIRE INPUT            ZXGP1310
C                PARAMETERS P1,P2,IP3,IP4, AND IP5.                     ZXGP1320
C                                                                       ZXGP1330
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       ZXGP1340
C                                                                       ZXGP1350
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZXGP1360
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZXGP1370
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZXGP1380
C                                                                       ZXGP1390
C-----------------------------------------------------------------------ZXGP1400
C                                                                       ZXGP1410
      SUBROUTINE ZXGSP  (F,P1,P2,IP3,IP4,IP5,A,B,TOL,XMIN,IER)          ZXGP1420
C                                  SPECIFICATIONS FOR ARGUMENTS         ZXGP1430
      INTEGER            IP3(1),IP4,IP5,IER                             ZXGP1440
      DOUBLE PRECISION     F,P1(1),P2(1),A,B,TOL,XMIN                   ZXGP1450
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZXGP1460
      DOUBLE PRECISION     C,FA,FB,H,V1,V2,FV1,FV2                      ZXGP1470
C                                  FIRST EXECUTABLE STATEMENT           ZXGP1480
      XMIN = A                                                          ZXGP1490
      IER = 129                                                         ZXGP1500
C                                  B MUST BE GREATER THAN A             ZXGP1510
      IF (B .LE. A) GO TO 9000                                          ZXGP1520
      IER = 130                                                         ZXGP1530
C                                  TOL MUST BE SMALLER THAN THE         ZXGP1540
C                                    INITIAL INTERVAL                   ZXGP1550
      IF (TOL .GE. (B-A)) GO TO 9000                                    ZXGP1560
      IER = 0                                                           ZXGP1570
C                                  COMPUTE THE FIBONACCI CONSTANT       ZXGP1580
      C = (3.0-SQRT(5.0))/2.0                                           ZXGP1590
C                                  COMPUTE THE INITIAL STEP             ZXGP1600
      H = C*(B-A)                                                       ZXGP1610
C                                  COMPUTE THE NEW POINTS               ZXGP1620
      V1 = A+H                                                          ZXGP1630
      V2 = B-H                                                          ZXGP1640
C                                  MAKE THE INITIAL FUNCTION EVALUATIONSZXGP1650
      FA = F(A,P1,P2,IP3,IP4,IP5)                                       ZXGP1660
      FB = F(B,P1,P2,IP3,IP4,IP5)                                       ZXGP1670
      FV1 = F(V1,P1,P2,IP3,IP4,IP5)                                     ZXGP1680
      FV2 = F(V2,P1,P2,IP3,IP4,IP5)                                     ZXGP1690
C                                  EACH ITERATION BEGINS HERE           ZXGP1700
    5 CONTINUE                                                          ZXGP1710
C                                  HAS THE INTERVAL BECOME TOO SMALL    ZXGP1720
      IF (A .GE. V1 .OR. V1 .GE. V2 .OR. V2 .GE. B) GO TO 40            ZXGP1730
C                                  FIND THE CURRENT MINIMUM             ZXGP1740
      IF (FV1 .GE. FV2) GO TO 10                                        ZXGP1750
C                                  V1 IS THE MINIMUM                    ZXGP1760
C                                  CHECK TO SEE IF THE FUNCTION IS      ZXGP1770
C                                    NOT UNIMODAL                       ZXGP1780
      IF (FV2 .GT. FB) GO TO 25                                         ZXGP1790
C                                  UPDATE THE INTERVAL                  ZXGP1800
C                                    V2 BECOMES THE NEW B               ZXGP1810
      B = V2                                                            ZXGP1820
C                                  IS THE INTERVAL SUFFICIENTLY SMALL   ZXGP1830
      IF (TOL .GE. (B-A)) GO TO 15                                      ZXGP1840
C                                  REDUCE THE INTERVAL FURTHER          ZXGP1850
      FB = FV2                                                          ZXGP1860
      V2 = V1                                                           ZXGP1870
      FV2 = FV1                                                         ZXGP1880
      H = C*(B-A)                                                       ZXGP1890
      V1 = A+H                                                          ZXGP1900
      FV1 = F(V1,P1,P2,IP3,IP4,IP5)                                     ZXGP1910
      GO TO 5                                                           ZXGP1920
C                                  V2 IS THE MINIMUM                    ZXGP1930
C                                  CHECK TO SEE IF THE FUNCTION IS      ZXGP1940
C                                    NOT UNIMODAL                       ZXGP1950
   10 IF (FV1 .GT. FA) GO TO 30                                         ZXGP1960
C                                  UPDATE THE INTERVAL                  ZXGP1970
C                                    V1 BECOMES THE NEW A               ZXGP1980
      A = V1                                                            ZXGP1990
C                                  IS THE INTERVAL SUFFICIENTLY SMALL   ZXGP2000
      IF (TOL .GE. (B-A)) GO TO 20                                      ZXGP2010
C                                  REDUCE THE INTERVAL FURTHER          ZXGP2020
      FA = FV1                                                          ZXGP2030
      V1 = V2                                                           ZXGP2040
      FV1 = FV2                                                         ZXGP2050
      H = C*(B-A)                                                       ZXGP2060
      V2 = B-H                                                          ZXGP2070
      FV2 = F(V2,P1,P2,IP3,IP4,IP5)                                     ZXGP2080
      GO TO 5                                                           ZXGP2090
C                                  CONVERGENCE OBTAINED. V1 OR A        ZXGP2100
C                                    IS THE MINIMUM                     ZXGP2110
   15 XMIN = V1                                                         ZXGP2120
      IF (FA .LT. FV1) XMIN = A                                         ZXGP2130
      GO TO 9005                                                        ZXGP2140
C                                  CONVERGENCE OBTAINED. V2 OR B        ZXGP2150
C                                    IS THE MINIMUM                     ZXGP2160
   20 XMIN = V2                                                         ZXGP2170
      IF (FB .LT. FV2) XMIN = B                                         ZXGP2180
      GO TO 9005                                                        ZXGP2190
C                                  FUNCTION IS NOT UNIMODAL. RETURN     ZXGP2200
C                                    THE NECESSARY PARAMETERS           ZXGP2210
   25 XMIN = V2                                                         ZXGP2220
      A = V1                                                            ZXGP2230
      GO TO 35                                                          ZXGP2240
   30 XMIN = V1                                                         ZXGP2250
      B = V2                                                            ZXGP2260
   35 IER = 131                                                         ZXGP2270
      GO TO 9000                                                        ZXGP2280
C                                  THE INTERVAL HAS BECOME TOO SMALL    ZXGP2290
   40 IER = 132                                                         ZXGP2300
      XMIN = A                                                          ZXGP2310
      IF (FB .LT. FA) XMIN = B                                          ZXGP2320
C                                                                       ZXGP2330
 9000 CONTINUE                                                          ZXGP2340
CCC      CALL UERTST (IER,'ZXGSP ')                                        ZXGP2350
 9005 RETURN                                                            ZXGP2360
      END                                                               ZXGP2370



      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J, TEMP
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
*                  WRITE(*,9000)
*                  STOP 779
                   TEMP = 0
                  END IF
            ELSE
*               WRITE(*,9000)
*               STOP 779
                TEMP = 0
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*       switch(*i){
*         case 1: return DBL_MIN;
*         case 2: return DBL_MAX;
*         case 3: return DBL_EPSILON/FLT_RADIX;
*         case 4: return DBL_EPSILON;
*         case 5: return log10((double)FLT_RADIX);
*         }
*       fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*       exit(1); return 0; /* some compilers demand return values */
*}
      END


      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
