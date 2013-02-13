C 
C  bsint.f
C  glbsint
C  
C  Created by Alexander Rudy on 2013-02-12.
C  Copyright 2013 Alexander Rudy. All rights reserved.
C 

       subroutine bsintegrate (derivs,nes,y,t0,t1,tacc,h0,
     &         mxstep,yout,tout)
           implicit real*8(a-h,o-z), integer*4(i-n)
     
           real*8 t0,t1,tacc,h0
           integer*4 mxstep,nes
           real*8 y(nes),yout(mxstep,nes),tout(mxstep)

           external derivs
           
           real*8 dydx(nes),yscal(nes)
           real*8 x
c..    Three body integrator from cartesian ic's
c..    currently set up to solve pythagorean problem.

c..    Uses the Bulirsch-Stoer method.

       
c..    Initial timestep
       htry=h0

c..    number of phase space dimensions
       nvar=nes
c..    number of scattering trials
       nloop=1
       Acc=tacc


c..      number of integration steps in current trial
         iflag=0
         
         x=t0
         iescape=0
         isscape=0
         ibelong1=1
         ibelong2=1
         istop=0

c..      Overall integration loop
2105     continue
c..        Conuter for integration loops
           iflag=iflag+1

c..        Use current timestep to take a Bulirsch-Stoer Integration Step
           call derivs(nes,x,y,dydx)

           do iscale=1,nvar
             yscal(iscale)=abs(y(iscale))+abs(htry*dydx(iscale))+1.e-30
           end do


           call bsstep(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  
           htry=hnext/2
           
           do i = 1, nes
                  yout(iflag,i) = y(i)
           end do
           tout(iflag) = x
           
           
c..        check for printout, or an unbound planet
           if(mod(iflag,1).eq.0) then
           
           
c..        


c..          check if number of integration steps exceeded:
             if(iflag.gt.mxstep) then

                    goto 2106

             end if

c..          check if physical time exceeded:
             if(x.gt.t1) then

                    goto 2106

             end if

c..          Outcome checklist finished
           end if

         goto 2105
c..      End of overall integration loop

2106     continue

c..      Successful termination
       
       
       return


       end

c....................................................................
c
      subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs) 
c
c....................................................................


c..   take a bulirsch-stoer timestep (numerical recipes)

      implicit real*8(a-h,o-z), integer*4(i-n)
      integer*4 nv,nmax,kmaxx,imax
      real*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),safe1,
     +       safe2,redmax,redmin,tiny,scalmx
      parameter (nmax=50,kmaxx=8,imax=kmaxx+1,safe1=0.25,safe2=0.7,
     +       redmax=1.e-5,redmin=0.7,tiny=1.e-30,scalmx=0.1)

    
      integer*4 i,iq,k,kk,km,kmax,kopt,nseq(imax)
      real*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,
     +       xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx),yerr(nmax),
     +       ysav(nmax),yseq(nmax)
      logical first,reduct
      external derivs
      data first/.true./,epsold/-1./
      data nseq /2,4,6,8,10,12,14,16,18/
      save a,alf,epsold, first,kmax,kopt,nseq,xnew

      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=safe1*eps
        a(1)=nseq(1)+1
        do k=1,kmaxx
          a(k+1)=a(k)+nseq(k+1)
        end do
        do iq=2,kmaxx
          do k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/
     +      ((a(iq+1)-a(1)+1.)*(2*k+1)))
          end do
        end do
        epsold=eps
        do kopt=2,kmaxx-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
        end do
1       kmax=kopt
      end if
      h=htry
      do i=1,nv
        ysav(i)=y(i)
      end do
      if(h.ne.hnext.or.x.ne.xnew) then
        first=.true.
        kopt=kmax
      end if
      reduct=.false.
2     do k=1,kmax
        xnew=x+h
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1) then
          errmax=TINY
          do i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
          end do
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/safe1)**(1./(2*km+1))
        end if
        if(k.ne.1.and.(k.ge.kopt-1.or.first)) then
          if(errmax.lt.1.) goto 4
          if(k.eq.kmax.or.k.eq.kopt+1) then
            red=safe2/err(km)
            goto 3
          else if(k.eq.kopt) then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            end if
          else if(kopt.eq.kmax) then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*
     +          safe2/err(km)
               goto 3  
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          end if
        end if
      end do
3     red=min(red,redmin)
      red=max(red,redmax)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do kk=1,km
        fact=max(err(kk),scalmx)
        work=fact*a(kk+1)
        if(work.lt.wrkmin) then
          scale=fact
          wrkmin=work
          kopt=kk+1
        end if
      end do
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),scalmx)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      end if
      return
      end

c....................................................................
c
      subroutine pzextr(iest,xest,yest,yz,dy,nv)
c
c....................................................................

c..   polynoical extrapolation (numerical recipes)
c..   (called by bsstep bulirsch-stoer integrator)

      implicit real*8(a-h,o-z), integer*4(i-n)
      integer*4 iest,nv,imax,nmax
      real*8 xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer*4 j,k1
      real*8 delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)
      save qcol,x
      x(iest)=xest
      do j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
      end do
      if(iest.eq.1) then
        do j=1,nv
          qcol(j,1)=yest(j)
        end do
      else
        do j=1,nv
          d(j)=yest(j)
        end do
        do k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
          end do
        end do
        do j=1,nv
          qcol(j,iest)=dy(j)
        end do
      end if
      return
      end

c.......................................................................
c
      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
c
c.......................................................................

c..   modified midpoint method to take one integration step (numerical recipes)
c..   (called by bsstep bulirsch-stoer integrator)

      integer*4 nstep,nvar,nmax
      real*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      external derivs
      parameter (nmax=50)
      integer*4 i,n
      real*8 h,h2,swap,x,ym(nmax),yn(nmax)
      h=htot/nstep
      do i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
      end do
      x=xs+h
      call derivs(nvar,x,yn,yout)
      h2=2.*h
      do n=2,nstep
        do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
        end do
        x=x+h
        call derivs(nvar,x,yn,yout)
      end do
      do i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
      end do
      return
      end

