     
C
      subroutine advave(curv2d)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion.      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real curv2d(im,jm)
      integer i,j
C
C     u-advection and diffusion:
C
C     Advective fluxes:
C
      do j=1,jm
        do i=1,im
          advua(i,j)=0.e0
        end do
      end do
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125e0*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          tps(i,j)=.25e0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25e0
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do
C
C     u-advection and diffusion:
C
      do j=1,jm
        do i=1,im
          advva(i,j)=0.e0
        end do
      end do
C
C     Advective fluxes:
C
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125e0*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25e0
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do
C
      if(mode.eq.2) then
C
        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5e0*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.25e0*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5e0*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.25e0*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25e0
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do
C
        do j=2,jmm1
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
        end do
C
        do j=3,jmm1
          do i=2,imm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
        end do
C
      endif
C
      return
C
      end
C
      subroutine advct(xflux,yflux,curv)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the horizontal portions of momentum      *
C *                advection well in advance of their use in advu and  *
C *                advv so that their vertical integrals (created in   *
C *                the main program) may be used in the external (2-D) *
C *                mode calculation.                                   *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real curv(im,jm,kb)
      real dtaam
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            curv(i,j,k)=0.e0
            advx(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25e0*((v(i,j+1,k)+v(i,j,k))
     $                         *(dy(i+1,j)-dy(i-1,j))
     $                         -(u(i+1,j,k)+u(i,j,k))
     $                         *(dx(i,j+1)-dx(i,j-1)))
     $                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do
C
C     Calculate x-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)
     $                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))
     $                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)
     $                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))
     $                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
C
            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)
     $                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            advy(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
C
C     Calculate y-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))
     $                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)
     $                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))
     $                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)
C
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)
     $                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=3,jmm1
          do i=2,imm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advq(qb,q,qf,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion, and  *
C *                vertical advection for turbulent quantities.        *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Do horizontal advection:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125e0*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Do horizontal diffusion:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection, add flux terms, then step forward in time:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt1(fb,f,fclim,ff,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is centred scheme, as originally provide in    *
C *                POM (previously called advt).                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
      do j=1,jm
        do i=1,im
           f(i,j,kb)=f(i,j,kbm1)
           fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do
C
C     Do advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do
C
C     Add diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection:
C
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.e0
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and then step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
C
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)
     $                 -dti2*ff(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt2(fb,f,fclim,ff,xflux,yflux,nitera,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is a first-order upstream scheme, which        *
C *                reduces implicit diffusion using the Smolarkiewicz  *
C *                iterative upstream scheme with an antidiffusive     *
C *                velocity.                                           *
C *                                                                    *
C *                It is based on the subroutines of Gianmaria Sannino *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option. It should be noted that    *
C *                this implementation does not include cross-terms    *
C *                which are in the original formulation.              *
C *                                                                    *
C *                fb,f,fclim,ff . as used in subroutine advt1         *
C *                xflux,yflux ... working arrays used to save memory  *
C *                nitera ........ number of iterations. This should   *
C *                                be in the range 1 - 4. 1 is         *
C *                                standard upstream differencing;     *
C *                                3 adds 50% CPU time to POM.         *
C *                sw ............ smoothing parameter. This should    *
C *                                preferably be 1, but 0 < sw < 1     *
C *                                gives smoother solutions with less  *
C *                                overshoot when nitera > 1.          *
C *                                                                    *
C *                Reference:                                          *
C *                                                                    *
C *                Smolarkiewicz, P.K.; A fully multidimensional       *
C *                  positive definite advection transport algorithm   *
C *                  with small implicit diffusion, Journal of         *
C *                  Computational Physics, 54, 325-362, 1984.         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real sw
      integer nitera
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      integer i,j,k,itera
C
C     Calculate horizontal mass fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            xmassflux(i,j,k)=0.e0
            ymassflux(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do
C
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do
C
C     Start Smolarkiewicz scheme:
C
      do itera=1,nitera
C
C     Upwind advection scheme:
C
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5e0
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))
C
              yflux(i,j,k)=0.5e0
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.e0
            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.e0
          end do
        end do
C
        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5e0
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do
C
C     Add net advective fluxes and step forward in time:
C
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
C
C     Calculate antidiffusion velocity:
C
        call smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C
        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do
C
C     End of Smolarkiewicz scheme
C
      end do
C
C     Add horizontal diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advu
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  advu                                                *
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                u-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            uf(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.e0*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advv
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                v-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            vf(i,j,k)=0.e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.e0*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine areas_masks
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates areas and masks.                         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j
C
C     Calculate areas of "t" and "s" cells:
C
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do
C
C     Calculate areas of "u" and "v" cells:
C
      do j=2,jm
        do i=2,im
          aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
C
      do j=1,jm
        aru(1,j)=aru(2,j)
        arv(1,j)=arv(2,j)
      end do
C
      do i=1,im
        aru(i,1)=aru(i,2)
        arv(i,1)=arv(i,2)
      end do
C
C     Initialise and set up free surface mask:
C
      do j=1,jm
        do i=1,im
          fsm(i,j)=0.e0
          dum(i,j)=0.e0
          dvm(i,j)=0.e0
          if(h(i,j).gt.1.e0) fsm(i,j)=1.e0
        end do
      end do
C
C     Set up velocity masks:
C
      do j=2,jm
        do i=2,im
          dum(i,j)=fsm(i,j)*fsm(i-1,j)
          dvm(i,j)=fsm(i,j)*fsm(i,j-1)
        end do
      end do
C
      return
C
      end
C
      subroutine baropg
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
          end do
        end do
      end do
C
C     Calculate x-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do
C
C     Calculate y-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine bcond(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Applies open boundary conditions.                   *
C *                                                                    *
C *                Closed boundary conditions are automatically        *
C *                enabled through specification of the masks, dum,    *
C *                dvm and fsm, in which case the open boundary        *
C *                conditions, included below, will be overwritten.    *
C *                                                                    *
C *                                The C-Grid                          *
C *                                **********                          *
C *                                                                    *
C *                The diagram is for the case where u and v are the   *
C *                primary boundary conditions together with t and     *
C *                s (co-located with el)                              * 
C *                                                                    *
C *                All interpolations are centered in space except     *
C *                those at lateral open boundary where an upstream    *
C *                Horizontal locations of e(el), t and s (etc.) are   *
C *                coincident.                                         *
C *                                                                    *
C *                People not acquainted with sigma coordinates have   *
C *                often asked what kind of boundary condition is      *
C *                applied along closed horizontal boundaries.         *
C *                Although the issue is not as important as it might  *
C *                be  for z-level grids, a direct answer is "half-    *
C *                slip" which, of course, is between free slip and    *
C *                non-slip.                                           *
C
C
C East and West end points for the C-grid in POM.
C
C                      west
C
C           v(1,j+1)=0           v(2,j+1) . . . . 
C          
C     ----<---<----<-----
C     |                 |       
C  u(1,j)   el(1,j)   u(2,j)=BC  el(2,j)   u(3,j) . . . . 
C             |                   |                 
C             -----<----<----<-----            
C
C           v(1,j)=0              v(2,j) . . . . 
C
C                                                    east
C
C                              . . . .  v(im-1,j+1)           v(im,j+1)=0
C                       
C                                           
C                 . . .  .  u(im-1,j)   el(im-1,j)  u(im,j)=BC  el(im,j)
C                                            |                   | 
C                                            ----->----->---->----
C
C                              . . . .   v(im-1,j)             v(im,j)=0
C
C  Notes:
C    1. The suffixes, f  or af, have been deleted.
C    2. All variables NOT designated as boundary condition (=BC) or set to 
C zero or obtained from an interior point are calculated points.
C    3. u(1,j) is never used but is obtained from the interior point for
C cosmetic output. Its counterpart, u(im+1,j), does not exist.
C    4. v=0 at i=1 and i=im are used as open inflow BC's unless specified
C otherwise.
C    5. The south and north extremal points are obtained from the above by 
C permuting u to v, v to u, i to j and j to i.


C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example, the governing boundary conditions are a radiation
C     condition on uaf in the east and in the west, and vaf in the north
C     and south. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do i=1,im
          elf(i,1)=elf(i,2)
          elf(i,jm)=elf(i,jmm1)
        end do
C
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
C
        do j=2,jmm1
C
C     East:
C
          uaf(im,j)=uabe(j)
     $               +rfe*sqrt(grav/h(imm1,j))
     $                         *(el(imm1,j)-ele(j))
          uaf(im,j)=ramp*uaf(im,j)
          vaf(im,j)=0.e0
C
C     West:
C
          uaf(2,j)=uabw(j)
     $              -rfw*sqrt(grav/h(2,j))
     $                        *(el(2,j)-elw(j))
          uaf(2,j)=ramp*uaf(2,j)
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
C
        end do
C
        do i=2,imm1
C
C     North:
C
          vaf(i,jm)=vabn(i)
     $               +rfn*sqrt(grav/h(i,jmm1))
     $                         *(el(i,jmm1)-eln(i))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.e0
C
C     South:
C
          vaf(i,2)=vabs(i)
     $              -rfs*sqrt(grav/h(i,2))
     $                        *(el(i,2)-els(i))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=0.e0
C
        end do
C
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Velocity (radiation conditions; smoothing is used in the direction
C     tangential to the boundaries):
C
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
            ga=sqrt(h(im,j)/hmax)
            uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)
     $                     +.25e0*u(imm1,j+1,k))
     $                  +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)
     $                    +.25e0*u(im,j+1,k))
            vf(im,j,k)=0.e0
C
C     West:
C
            ga=sqrt(h(1,j)/hmax)
            uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k)
     $                    +.25e0*u(3,j+1,k))
     $                 +(1.e0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k)
     $                   +.25e0*u(2,j+1,k))
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            ga=sqrt(h(i,jm)/hmax)
            vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k)
     $                     +.25e0*v(i+1,jmm1,k))
     $                  +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k)
     $                    +.25e0*v(i+1,jm,k))
            uf(i,jm,k)=0.e0
C
C     South:
C
            ga=sqrt(h(i,1)/hmax)
            vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)
     $                    +.25e0*v(i+1,3,k))
     $                 +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)
     $                   +.25e0*v(i+1,2,k))
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=0.e0
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
        do k=1,kbm1
          do j=1,jm
C
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(imm1,j))
                uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
              endif
            endif
C
C     West:
C
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
            else
              uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(2,j))
                uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do i=1,im
C
C     North:
C
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
            else
              uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
              endif
            endif
C
C     South:
C
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
            else
              uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,2))
                uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
        do k=1,kb
          do j=1,jm
C
C     East:
C
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
              uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
            endif
C
C     West:
C
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
            else
              uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
            endif
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
            else
              uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
            endif
C
C     South:
C
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
            else
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
            endif
          end do
        end do
C
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcondorl(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  This is an optional subroutine replacing  bcond and *
C *                using Orlanski's scheme (J. Comp. Phys. 21, 251-269,*
C *                1976), specialized for the seamount problem. To     *
C *                make it work for the seamount problem, I (G.M.)     *
C *                have had to add an extra condition on an "if"       *
C *                statement in the t and s open boundary conditions,  *
C *                which involves the sign of the normal velocity.     *
C *                Thus:                                               *
C *                                                                    *
C *            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k), *
C *                                                                    *
C *                plus 3 others of the same kind.                     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real cl,denom
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example the governing boundary conditions are a radiation
C     condition on uaf(im,j) in the east and an inflow uaf(2,j) in the
C     west. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do  j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
C
        do j=2,jmm1
C
C     West:
C
          uaf(2,j)=ramp*uabw(j)-sqrt(grav/h(2,j))*(el(2,j)-elw(j))
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
C
C     East:
C
          uaf(im,j)=ramp*uabe(j)
     $               +sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
          vaf(im,j)=0.e0
C
        end do
C
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Eastern and western radiation boundary conditions according to
C     Orlanski's explicit scheme:
C
        do k=1,kbm1
          do j=2,jmm1
C
C     West:
C
            denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(3,j,k)-uf(3,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k))
     $                 /(1.e0+cl)
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
C
C     East:
C
            denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k))
     $                  /(1.e0+cl)
            vf(im,j,k)=0.e0
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
        do k=1,kbm1
          do j=1,jm
C
C     West:
C
            ubw(j,k)=ub(2,j,k)
            denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(2,j,k)-uf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k)
C
            denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(2,j,k)-vf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) vf(1,j,k)=sbw(j,k)
C
C     East:
C
            ube(j,k)=ub(im,j,k)
            denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k))
     $                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) uf(im,j,k)=tbe(j,k)
C
            denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k))
     $                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) vf(im,j,k)=sbe(j,k)
C
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
        do k=1,kb
C
          do j=1,jm
            uf(im,j,k)=1.e-10
            vf(im,j,k)=1.e-10
            uf(1,j,k)=1.e-10
            vf(1,j,k)=1.e-10
          end do
C
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
c       subroutine box
c C **********************************************************************
c C *                                                                    *
c C * FUNCTION    :  Sets up conservation box problem.                   *
c C *                                                                    *
c C *                This basin uses the same grid as the seamount       *
c C *                problem, but it has a flat bottom, is surrounded by *
c C *                walls and is initialised with uniform salinity and  *
c C *                temperature. It is forced by a surface input of     *
c C *                water of the same temperature and salinity as the   *
c C *                water in the basin. Therefore, the temperature and  *
c C *                salinity in the basin should not change, and the    *
c C *                free surface should fall at a rate vflux. It is also*
c C *                forced by a steady atmospheric pressure field which *
c C *                depresses the southwestern half of the model by 1 m *
c C *                and elevates the northeastern half of the model by  *
c C *                1 m.                                                *
c C *                                                                    *
c C *                Since this problem defines its own fixed e_atmos,   *
c C *                tatm, satm and e_atmos, comment out corresponding   *
c C *                declarations after the do 9000 statement in main    *
c C *                program.                                            *
c C **********************************************************************
c C
c       implicit none
c C
c       include 'pom2k.c'
c C
c       real depth,delx,tatm,satm
c       integer i,j,k
c C
c C     Water depth:
c C
c       depth=4500.e0
c C
c C     Grid size:
c C
c       delx=8000.e0
c C
c C     Set up grid dimensions, areas of free surface cells, and
c C     Coriolis parameter:
c C
c       do j=1,jm
c         do i=1,im
c C
c C     For constant grid size:
c C
c C         dx(i,j)=delx
c C         dy(i,j)=delx
c C
c C     For variable grid size:
c C
c           dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
c           dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
c C
c           cor(i,j)=1.e-4
c C
c         end do
c       end do
c C
c C     Calculate horizontal coordinates of grid points and rotation
c C     angle.
c C
c C     NOTE that this is introduced solely for the benefit of any post-
c C     processing software, and in order to conform with the requirements
c C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
c C
c C     There are four horizontal coordinate systems, denoted by the
c C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
c C     "e" is an elevation point and "c" is a cell corner), as shown
c C     below. In addition, "east_*" is an easting and "north_*" is a
c C     northing. Hence the coordinates of the "u" points are given by
c C     (east_u,north_u).
c C
c C     Also, if the centre point of the cell shown below is at
c C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
c C     the coordinates of the western of the two "u" points,
c C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
c C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
c C     coordinates of the southwestern corner point of the cell. The
c C     southwest corner of the entire grid is at
c C     (east_c(1,1),north_c(1,1)).
c C
c C                      |              |
c C                    --c------v-------c--
c C                      |              |
c C                      |              |
c C                      |              |
c C                      |              |
c C                      u      e       u
c C                      |              |
c C                      |              |
c C                      |              |
c C                      |              |
c C                    --c------v-------c--
c C                      |              |
c C
c C
c C     NOTE that the following calculation of east_c and north_c only
c C     works properly for a rectangular grid with east and north aligned
c C     with i and j, respectively:
c C
c       do j=1,jm
c         east_c(1,j)=0.e0
c         do i=2,im
c           east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
c         end do
c       end do
c C
c       do i=1,im
c         north_c(i,1)=0.e0
c         do j=2,jm
c           north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
c         end do
c       end do
c C
c C     The following works properly for any grid:
c C
c C     Elevation points:
c C
c       do j=1,jm-1
c         do i=1,im-1
c           east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
c      $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
c           north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
c      $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
c         end do
c       end do
c C
c C     Extrapolate ends:
c C
c       do i=1,im-1
c         east_e(i,jm)
c      $    =((east_c(i,jm)+east_c(i+1,jm))*3.e0
c      $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
c         north_e(i,jm)
c      $    =((north_c(i,jm)+north_c(i+1,jm))*3.e0
c      $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
c       end do
c C
c       do j=1,jm-1
c         east_e(im,j)
c      $    =((east_c(im,j)+east_c(im,j+1))*3.e0
c      $       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
c         north_e(im,j)
c      $    =((north_c(im,j)+north_c(im,j+1))*3.e0
c      $       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
c       end do
c C
c       east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
c      $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
c       north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
c      $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
c C
c C     u-points:
c C
c       do j=1,jm-1
c         do i=1,im
c           east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
c           north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
c         end do
c       end do
c C
c C     Extrapolate ends:
c C
c       do i=1,im
c         east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
c         north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
c       end do
c C
c C     v-points:
c C
c       do j=1,jm
c         do i=1,im-1
c           east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
c           north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
c         end do
c       end do
c C
c C     Extrapolate ends:
c C
c       do j=1,jm
c         east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
c         north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
c       end do
c C
c C     rot is the angle (radians, anticlockwise) of the i-axis relative
c C     to east, averaged to a cell centre:
c C
c C     (NOTE that the following calculation of rot only works properly
c C     for this particular rectangular grid)
c C
c       do j=1,jm
c         do i=1,im
c           rot(i,j)=0.e0
c         end do
c       end do
c C
c C     Define depth:
c C
c       do i=1,im
c         do j=1,jm
c           h(i,j)=depth
c         end do
c       end do
c C
c C     Close the north and south boundaries:
c C
c       do i=1,im
c         h(i,1)=1.e0
c         h(i,jm)=1.e0
c       end do
c C
c C     Close the east and west boundaries:
c C
c       do j=1,jm
c         h(1,j)=1.e0
c         h(im,j)=1.e0
c       end do
c C
c C     Calculate areas and masks:
c C
c       call areas_masks
c C
c C     Adjust bottom topography so that cell to cell variations
c C     in h do not exceed parameter slmax:
c C
c       if(slmax.lt.1.e0) call slpmax
c C
c C     Set tbias and sbias here for test (tbias and sbias would
c C     normally only be set in the main program):
c C
c       tbias=10.e0
c       sbias=20.e0
c       write(6,1) tbias,sbias
c     1 format(/' tbias and sbias changed in subroutine box to:'/
c      $         2f10.3//)
c C
c C     Set initial conditions:
c C
c       do k=1,kbm1
c         do j=1,jm
c           do i=1,im
c             tb(i,j,k)=20.e0-tbias
c             sb(i,j,k)=35.e0-sbias
c             tclim(i,j,k)=tb(i,j,k)
c             sclim(i,j,k)=sb(i,j,k)
c           end do
c         end do
c       end do
c C
c C     Initialise uab and vab as necessary
c C     (NOTE that these have already been initialised to zero in the
c C     main program):
c C
c       do j=1,jm
c         do i=1,im
c C     No conditions necessary for this problem
c         end do
c       end do
c C
c C     Set surface boundary conditions, e_atmos, vflux, wusurf,
c C     wvsurf, wtsurf, wssurf and swrad, as necessary
c C     (NOTE:
c C      1. These have all been initialised to zero in the main program.
c C      2. The temperature and salinity of inflowing water must be
c C         defined relative to tbias and sbias.):
c C
c       do j=1,jm
c         do i=1,im
c           if(i+j-57.le.0) then
c             e_atmos(i,j)=1.e0
c           else
c             e_atmos(i,j)=-1.e0
c           endif
c C
c C     Ensure atmospheric pressure cannot make water depth go negative:
c C
c           e_atmos(i,j)=min(e_atmos(i,j),h(i,j))
c C
c           vfluxf(i,j)=-0.0001e0
c C
c C     See main program, just after "Begin numerical integration", for
c C     an explanation of these terms:
c C 
c           tatm=20.e0
c           satm=35.e0
c C
c         end do
c       end do
c C
c C     Initialise elb, etb, dt and aam2d:
c C
c       do j=1,jm
c         do i=1,im
c           elb(i,j)=-e_atmos(i,j)
c           etb(i,j)=-e_atmos(i,j)
c           dt(i,j)=h(i,j)-e_atmos(i,j)
c           aam2d(i,j)=aam(i,j,1)
c         end do
c       end do
c C
c       call dens(sb,tb,rho)
c C
c C     Generated horizontally averaged density field (in this
c C     application, the initial condition for density is a function
c C     of z (the vertical cartesian coordinate) -- when this is not
c C     so, make sure that rmean has been area averaged BEFORE transfer
c C     to sigma coordinates):
c C
c       do k=1,kbm1
c         do j=1,jm
c           do i=1,im
c             rmean(i,j,k)=rho(i,j,k)
c           end do
c         end do
c       end do
c C
c C     Set lateral boundary conditions, for use in subroutine bcond
c C     (in this problem, all lateral boundaries are closed through
c C     the specification of the masks fsm, dum and dvm):
c C
c       rfe=1.e0
c       rfw=1.e0
c       rfn=1.e0
c       rfs=1.e0
c C
c C     Set thermodynamic boundary conditions (for the seamount
c C     problem, and other possible applications, lateral thermodynamic
c C     boundary conditions are set equal to the initial conditions and
c C     are held constant thereafter - users may, of course, create
c C     variable boundary conditions):
c C
c       do k=1,kbm1
c C
c         do j=1,jm
c           tbe(j,k)=tb(im,j,k)
c           tbw(j,k)=tb(1,j,k)
c           sbe(j,k)=sb(im,j,k)
c           sbw(j,k)=sb(1,j,k)
c         end do
c C
c         do i=1,im
c           tbn(i,k)=tb(i,jm,k)
c           tbs(i,k)=tb(i,1,k)
c           sbn(i,k)=sb(i,jm,k)
c           sbs(i,k)=sb(i,1,k)
c         end do
c C
c       end do
c C
c       return
c C
c       end
c C
      subroutine dens(si,ti,rhoo)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates (density-1000.)/rhoref.                  *
C *                                                                    *
C *                (see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech.,  *
C *                609-611.)                                           *
C *                                                                    *
C *                ti is potential temperature                         *
C *                                                                    *
C *                If using 32 bit precision, it is recommended that   *
C *                cr,p,rhor,sr,tr,tr2,tr3 and tr4 be made double      *
C *                precision, and the "e"s in the constants be changed *
C *                to "d"s.                                            *
C *                                                                    *
C * NOTE: if pressure is not used in dens, buoyancy term (boygr)       *
C *       in profq must be changed (see note in profq)                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      real cr,p,rhor,sr,tr,tr2,tr3,tr4
      integer i,j,k
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
C
            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr
C
C     Approximate pressure in units of bars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5
C
            rhor=-0.157406e0+6.793952e-2*tr
     $            -9.095290e-3*tr2+1.001685e-4*tr3
     $            -1.120083e-6*tr4+6.536332e-9*tr4*tr
C
            rhor=rhor+(0.824493e0-4.0899e-3*tr
     $            +7.6438e-5*tr2-8.2467e-7*tr3
     $            +5.3875e-9*tr4)*sr
     $            +(-5.72466e-3+1.0227e-4*tr
     $            -1.6546e-6*tr2)*abs(sr)**1.5
     $            +4.8314e-4*sr*sr
C
            cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2
     $          +1.34e0*(sr-35.e0)
            rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr))
C
            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)
C
          end do
        end do
      end do
c      write(6,*)'10,10,1,rho',si(10,10,1),ti(10,10,1),rhoo(10,10,1)
c      pause 
C
      return
C
      end
C
      subroutine depth
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Establishes the vertical sigma grid with log        *
C *                distributions at the top and bottom and a linear    *
C *                distribution in between. The number of layers of    *
C *                reduced thickness are kl1-2 at the surface and      *
C *                kb-kl2-1 at the bottom. kl1 and kl2 are defined in  *
C *                the main program. For no log portions, set kl1=2    *
C *                and kl2=kb-1.                                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delz
      integer kdz(12)
      integer k
      integer mend,mdir
      character dir*80
C
      data kdz/1,1,2,4,8,16,32,64,128,256,512,1024/
C
      z(1)=0.e0
C
      do k=2,kl1
        z(k)=z(k-1)+kdz(k-1)
      end do
C
      delz=z(kl1)-z(kl1-1)
C
      do k=kl1+1,kl2
        z(k)=z(k-1)+delz
      end do
C
      do k=kl2+1,kb
        dz(k)=float(kdz(kb-k+1))*delz/float(kdz(kb-kl2))
        z(k)=z(k-1)+dz(k)
      end do
C
      do k=1,kb
        z(k)=-z(k)/z(kb)
      end do
C
      do k=1,kb-1
        zz(k)=0.5e0*(z(k)+z(k+1))
      end do
C
      zz(kb)=2.e0*zz(kb-1)-zz(kb-2)
C
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
C
      dz(kb)=0.e0
      dzz(kb)=0.e0
C
      write(6,1)
    1 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
C
      do k=1,kb
        write(6,2) k,z(k),zz(k),dz(k),dzz(k)
    2   format((' ',i5,4f10.3))
      end do
C
      write(6,3)
    3 format(//)
C
      mend=20190414

      open(70,file=dir(1:mdir)//'sigma_z',
     1        form='unformatted')
            write(70)z(1:kb),zz(1:kb),dz(1:kb),dzz(1:kb),mend
         close(70)

      return
C
      end
C
      subroutine findpsi
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  findpsi                                             *
C *                                                                    *
C * FUNCTION    :  Calculates the stream function, first assuming      *
C *                zero on the southern boundary and then, using the   *
C *                values on the western boundary, the stream function *
C *                is calculated again. If the elevation field is near *
C *                steady state, the two calculations should agree;    *
C *                otherwise not.                                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j
C
      do j=1,jm
        do i=1,im
          psi(i,j)=0.e0
        end do
      end do
C
C     Sweep northward:
C
      do j=2,jmm1
        do i=2,im
          psi(i,j+1)=psi(i,j)
     $                +.25e0*uab(i,j)*(d(i,j)+d(i-1,j))
     $                  *(dy(i,j)+dy(i-1,j))
        end do
      end do
C
      call prxy('Streamfunction, psi from u              ',
     $          time,psi,im,iskp,jm,jskp,0.e0)
C
C    Sweep eastward:
C
      do j=2,jm
        do i=2,imm1
          psi(i+1,j)=psi(i,j)
     $                -.25e0*vab(i,j)*(d(i,j)+d(i,j-1))
     $                  *(dx(i,j)+dx(i,j-1))
        end do
      end do
C
      call prxy('Streamfunction, psi from v              ',
     $          time,psi,im,iskp,jm,jskp,0.e0)
C
      return
C
      end
C
c       subroutine file2ic
c C **********************************************************************
c C *                                                                    *
c C * FUNCTION    :  Sets up my own problem.                             *
c C *                                                                    *
c C * This example read IC from IC.dat file, generated by GRID.f in      *
c C * GRID-DATA directory. Only minimal number of fields are read,       *
c C * while others are calculated here.                                  *
c C *                                                                    *
c C **********************************************************************
c C
c       implicit none
c C
c       include 'pom2k.c'
c C
c       real rad,re,dlat,dlon,cff
c       integer i,j,k,m
c       character*5 field
c       rad=0.01745329
c       re=6371.E3
c C
c       write(6,'(/,'' Read grid and initial conditions '',/)')
c C
c C--- 1D ---
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') z
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') zz
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') dz
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') dzz
c C--- 2D ---
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') east_e
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') north_e
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') h
c C--- 3D ---
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') t
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') s
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') rmean
c C--- Constant wind stress read here
c C (for time dep. read in loop 9000 & interpolate in time)
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') wusurf
c       read(40,'(a5)') field
c       write(6,'(a5)') field
c        read(40,'(8E12.5)') wvsurf
c C
c C --- print vertical grid distribution
c C
c       write(6,2)
c     2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
c       write(6,'(''  '',/)')
c       do k=1,kb
c         write(6,3) k,z(k),zz(k),dz(k),dzz(k)
c     3   format((' ',i5,4f10.3))
c       end do
c       write(6,'(''  '',//)')
c C
c C --- calc. surface & lateral BC from climatology
c C
c         do j=1,jm
c           do i=1,im
c              tsurf(i,j)=t(i,j,1)
c              ssurf(i,j)=s(i,j,1)
c             do k=1,kb
c               tclim(i,j,k)=t(i,j,k)
c               sclim(i,j,k)=s(i,j,k)
c             end do
c           end do
c         end do
c C
c C                    --- EAST & WEST BCs ---
c         do j=1,jm
c               ele(j)=0.
c               elw(j)=0.
c C --- other vel. BCs (fixed in time) can be specified here
c               uabe(j)=0.
c               uabw(j)=0.
c             do k=1,kb
c               ubw(j,k)=0.
c               ube(j,k)=0.
c               tbw(j,k)=tclim(1,j,k)
c               sbw(j,k)=sclim(1,j,k)
c               tbe(j,k)=tclim(im,j,k)
c               sbe(j,k)=sclim(im,j,k)
c             end do
c         end do
c C                    --- NORTH & SOUTH BCs ---
c         do i=1,im
c               els(i)=0.
c               eln(i)=0.
c               vabs(i)=0.
c               vabn(i)=0.
c             do k=1,kb
c               vbs(i,k)=0.
c               vbn(i,k)=0.
c               tbs(i,k)=tclim(i,1,k)
c               sbs(i,k)=sclim(i,1,k)
c               tbn(i,k)=tclim(i,jm,k)
c               sbn(i,k)=sclim(i,jm,k)
c             end do
c         end do
c C
c C     Set initial conditions:
c C
c       do k=1,kb
c         do j=1,jm
c           do i=1,im
c             tb(i,j,k)=t(i,j,k)
c             sb(i,j,k)=s(i,j,k)
c             ub(i,j,k)=0.
c             vb(i,j,k)=0.
c           end do
c         end do
c       end do
c C
c       call dens(sb,tb,rho)
c C
c C --- calc. Curiolis Parameter
c C
c         do j=1,jm
c           do i=1,im
c             cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
c             aam2d(i,j)=aam(i,j,1)
c             elb(i,j)=0.
c             etb(i,j)=0.
c             dt(i,j)=h(i,j)
c           end do
c         end do
c C
c         do j=1,jm
c           do i=2,im-1
c             dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
c      1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
c           end do
c             dx(1,j)=dx(2,j)
c             dx(im,j)=dx(im-1,j)
c         end do
c C
c         do i=1,im
c           do j=2,jm-1
c             dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
c      1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
c           end do
c             dy(i,1)=dy(i,2)
c             dy(i,jm)=dy(i,jm-1)
c         end do
c C
c C     Calculate areas and masks:
c C
c       call areas_masks
c C
c C
c C --- the following grids are needed only for netcdf plotting
c C
c C     Corner of cell points:
c C
c       do j=2,jm
c         do i=2,im
c           east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
c      $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
c           north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
c      $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
c         end do
c       end do
c C
c C
c C     Extrapolate ends (approx.):
c C
c       do i=2,im
c         east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
c         north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
c       end do
c         east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
c C
c       do j=2,jm
c         east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
c         north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
c       end do
c         north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
c C
c C     u-points:
c C
c       do j=1,jm-1
c         do i=1,im
c           east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
c           north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
c         end do
c       end do
c C
c C     Extrapolate ends:
c C
c       do i=1,im
c         east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
c         north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
c       end do
c C
c C     v-points:
c C
c       do j=1,jm
c         do i=1,im-1
c           east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
c           north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
c         end do
c       end do
c C
c C     Extrapolate ends:
c C
c       do j=1,jm
c         east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
c         north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
c       end do
c C
c C     rot is the angle (radians, anticlockwise) of the i-axis relative
c C     to east, averaged to a cell centre: (only needed for CDF plotting)
c C
c       do j=1,jm
c         do i=1,im-1
c           rot(i,j)=0.
c           dlat=north_e(i+1,j)-north_e(i,j)
c           dlon= east_e(i+1,j)- east_e(i,j)
c            if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
c         end do
c        rot(im,j)=rot(im-1,j)
c       end do
c C
c C     Set lateral boundary conditions, for use in subroutine bcond
c C     set all=0 for closed BCs.
c C     Values=0 for vel BC only, =1 is combination of vel+elev.
c       rfe=0.e0
c       rfw=0.e0
c       rfn=0.e0
c       rfs=0.e0
c C
c       return
c       end
c C
      subroutine printall
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  printall                                            *
C *                                                                    *
C * FUNCTION    :  Prints a set of outputs to device 6                 *
C *                                                                    *
C *                Edit as approriate.                                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer io(100),jo(100),ko(100)
C
      include 'pom2k.c'
C
C     2-D horizontal fields:
C
          call prxy('Depth-averaged u, uab                   ',
     $              time,uab,im,iskp,jm,jskp,0.e0)
C
          call prxy('Depth-averaged v, vab                   ',
     $              time,vab,im,iskp,jm,jskp,0.e0)
C
          call prxy('Surface elevation, elb                  ',
     $              time,elb,im,iskp,jm,jskp,0.e0)
C
c         call prxy(' egf ',time,egf,im,iskp,jm,jskp,0.e0)
c         call prxy(' utf ',time,utf,im,iskp,jm,jskp,0.e0)
c         call prxy(' vtf ',time,vtf,im,iskp,jm,jskp,0.e0)
c
C     Calculate and print streamfunction:
C
          call findpsi
C
          if(mode.ne.2) then
C
C     2-D horizontal sections of 3-D fields:
C
C     Set levels for output:
C
            ko(1)=1
            ko(2)=kb/2
            ko(3)=kb-1
C
            call prxyz('x-velocity, u                           ',
     $                 time,u    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
            call prxyz('y-velocity, v                           ',
     $                 time,v    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
            ko(1)=2
            call prxyz('z-velocity, w                           ',
     $                 time,w    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
            ko(1)=1
C
            call prxyz('Potential temperature, t                ',
     $                 time,t    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
C
            call prxyz('Salinity, s                              ',
     $                 time,s    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
C
            call prxyz('(density-1000)/rhoref, rho              ',
     $                 time,rho  ,im,iskp,jm,jskp,kb,ko,3,1.e-5)
C
c           call prxyz('Turbulent kinetic energy x 2, q2        ',
c    $                 time,q2   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
c           call prxyz('Turbulent length scale, l               ',
c    $                 time,l    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
            call prxyz('Horizontal kinematic viscosity, aam     ',
     $                 time,aam  ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
            call prxyz('Vertical kinematic viscosity, km        ',
     $                 time,km   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
c           call prxyz('Vertical kinematic diffusivity, kh      ',
c    $                 time,kh   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
C
C     Vertical sections of 3-D fields, normal to j-axis:
C
C     Set sections for output:
C
            jo(1)=1
            jo(2)=jm/2
            jo(3)=jm-1
C
            call prxz('x-velocity, u                           ',
     $                time,u    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
C
            call prxz('y-velocity, v                           ',
     $                time,v    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
C
            call prxz('z-velocity, w                           ',
     $                time,w    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
C
            call prxz('Potential temperature, t                ',
     $                time,t    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
C
            call prxz('Salinity, s                             ',
     $                time,s    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
C
            call prxz('(density-1000)/rhoref, rho              ',
     $                time,rho  ,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
C
c           call prxz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
C
c           call prxz('Turbulent length scale, l               ',
c    $                time,l    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
C
c           call prxz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
C
c           call prxz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
C
c           call prxz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
C
C     Vertical sections of 3-D fields, normal to i-axis:
C
C     Set sections for output:
C
            io(1)=1
            io(2)=im/2
            io(3)=im-1
C
            call pryz('x-velocity, u                           ',
     $                time,u    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
C
            call pryz('y-velocity, v                           ',
     $                time,v    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
C
            call pryz('z-velocity, w                           ',
     $                time,w    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
C
            call pryz('Potential temperature, t                ',
     $                time,t    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
C
c           call pryz('Salinity x rho / rhoref, s              ',
c    $                time,s    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
C
c           call pryz('(density-1000)/rhoref, rho              ',
c    $                time,rho  ,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
C
c           call pryz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
C
c           call pryz('Turbulent length scale, l               ',
c    $                time,l    ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
C
c           call pryz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
C
c           call pryz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
C
c           call pryz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
C
          endif
C
      return
C
      end
C
      subroutine profq(sm,sh,dh,cc)
C **********************************************************************
C *                                        Updated: Sep. 24, 2003      *
C * FUNCTION    :  Solves for q2 (twice the turbulent kinetic energy), *
C *                q2l (q2 x turbulent length scale), km (vertical     *
C *                kinematic viscosity) and kh (vertical kinematic     *
C *                diffusivity), using a simplified version of the     *
C *                level 2 1/2 model of Mellor and Yamada (1982).      *
C * In this version, the Craig-Banner sub-model whereby breaking wave  * 
C * tke is injected into the surface is included. However, we use an   *
C * analytical solution to the near surface tke equation to solve for  *
C * q2 at the surface giving the same result as C-B diffusion. The new *
C * scheme is simpler and more robust than the latter scheme.          *     
C *                                                                    *
C * References                                                         *
C *   Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence *
C *     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, *
C *     1994.                                                          *
C *   Ezer, T., On the seasonal mixed-layer simulated by a basin-scale *
C *     ocean model and the Mellor-Yamada turbulence scheme,           *
C *     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.                *
C *   Mellor, G.L. and T. Yamada, Development of a turbulence          *
C *     closure model for geophysical fluid fluid problems,            *
C *     Rev. Geophys. Space Phys., 20, 851-875, 1982.                  *
C *   Mellor, G. L., One-dimensional, ocean surface layer modeling,    *
C *     a problem and a solution. J. Phys. Oceanogr., 31(3), 790-809,  *
C *     2001.                                                          *
C *   Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface    *
C *     thermal response, J. Phys. Oceanogr., 2003.                    *
C *   Stacey, M. W., Simulations of the wind-forced near-surface       *
C *     circulation in Knight Inlet: a parameterization of the         *
C *     roughness length. J. Phys. Oceanogr., 29, 1363-1367, 1999.     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real sm(im,jm,kb),sh(im,jm,kb),cc(im,jm,kb)
      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
      real prod(im,jm,kb),kn(im,jm,kb)
      real a1,a2,b1,b2,c1
      real coef1,coef2,coef3,coef4,coef5
      real const1,e1,e2,ghc
      real p,sef,sp,tp
      real l0(im,jm)
      real cbcnst,surfl,shiw
      real utau2, df0,df1,df2 
C
      integer i,j,k,ki
C
      equivalence (prod,kn)
C
      data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
      data e1/1.8e0/,e2/1.33e0/
      data sef/1.e0/
      data cbcnst/100./surfl/2.e5/shiw/0.0/
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
C
C     Surface and bottom boundary conditions:
C
      const1=(16.6e0**(2.e0/3.e0))*sef
C
C initialize fields that are not calculated on all boundaries
C but are later used there
      do i=1,im
        ee(i,jm,1)=0.
        gg(i,jm,1)=0.
        l0(i,jm)=0.
      end do
      do j=1,jm
        ee(im,j,1)=0.
        gg(im,j,1)=0.
        l0(im,j)=0.
      end do
      do i=1,im
      do j=1,jm
       do k=2,kbm1
        prod(i,j,k)=0.
       end do
      end do
      end do
C
      do j=1,jmm1
        do i=1,imm1
          utau2=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                  +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
C Wave breaking energy- a variant of Craig & Banner (1994)
C see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.e0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2 
C Surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2/grav
C
          uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2
     $                   +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
C
C    Calculate speed of sound squared:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
C
C     Calculate pressure in units of decibars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
            cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp -.045e0*tp**2
     $                 +1.34e0*(sp-35.0e0)
            cc(i,j,k)=cc(i,j,k)
     $                 /sqrt((1.e0-.01642e0*p/cc(i,j,k))
     $                   *(1.e0-0.40e0*p/cc(i,j,k)**2))
          end do
        end do
      end do
C
C     Calculate buoyancy gradient:
C
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                    /(dzz(k-1)* h(i,j))
C *** NOTE: comment out next line if dens does not include pressure
     $      +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.028e0)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.e0
          gh(i,j,1)=0.e0
          gh(i,j,kb)=0.e0
        end do
      end do
C
C    Calculate production of turbulent kinetic energy:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25e0*sef
     $                   *((u(i,j,k)-u(i,j,k-1)
     $                      +u(i+1,j,k)-u(i+1,j,k-1))**2
     $                     +(v(i,j,k)-v(i,j,k-1)
     $                      +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                   /(dzz(k-1)*dh(i,j))**2
C   Add shear due to internal wave field
     $             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
C
C  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000),
C  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
      ghc=-6.0e0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.e0
C It is unclear yet if diss. corr. is needed when surf. waves are included.
c           if(gh(i,j,k).lt.0.e0)
c    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
c           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(2.e0*dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
C
      do j=1,jm
        do i=1,im
          ee(i,j,2)=0.e0
          gg(i,j,2)=0.e0
          vf(i,j,kb)=0.e0
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.e0+e2*((1.e0/abs(z(k)-z(1))
     $                               +1.e0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
              uf(i,j,k)=small
              vf(i,j,k)=0.1*dt(i,j)*small
            endif
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves for km and kh:
C
      coef4=18.e0*a1*a1+9.e0*a1*a2
      coef5=9.e0*a1*a2
C
C     Note that sm and sh limit to infinity when gh approaches 0.0288:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
            coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
            coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            kn(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(kn(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
            km(i,j,k)=(kn(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
            kh(i,j,k)=(kn(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
          end do
        end do
      end do
C cosmetics: make boundr. values as interior
C (even if not used, printout otherwise may show strange values)
      do k=1,kb
        do i=1,im
           km(i,jm,k)=km(i,jmm1,k)*fsm(i,jm)
           kh(i,jm,k)=kh(i,jmm1,k)*fsm(i,jm)
           km(i,1,k)=km(i,2,k)*fsm(i,1)
           kh(i,1,k)=kh(i,2,k)*fsm(i,1)
        end do
        do j=1,jm
           km(im,j,k)=km(imm1,j,k)*fsm(im,j)
           kh(im,j,k)=kh(imm1,j,k)*fsm(im,j)
           km(1,j,k)=km(2,j,k)*fsm(1,j)
           kh(1,j,k)=kh(2,j,k)*fsm(1,j)
        end do
      end do
C
      return
C
      end
C
c ---------------------------------------------------------------------
C
      subroutine proft(f,wfsurf,fsurf,nbc,dh)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of temperature and    *
C *                salinity using method described by Richmeyer and    *
C *                Morton.                                             *
C *                                                                    *
C *                Irradiance parameters are from Paulson and Simpson. *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                Paulson, C. A., and J. Simpson, 1977: Irradiance    *
C *                  measurements in the upper ocean, J. Phys.         *
C *                  Oceanogr., 7, 952-956.                            *
C *                                                                    *
C *                NOTES:                                              *
C *                                                                    *
C *                (1) wfsurf and swrad are negative values when water *
C *                    column is warming or salt is being added.       *
C *                                                                    *
C *                (2) nbc may only be 1 and 3 for salinity.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)
      integer nbc
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki
C
C-----------------------------------------------------------------------
C
C     Irradiance parameters after Paulson and Simpson:
C
C       ntp               1      2       3       4       5
C   Jerlov type           i      ia      ib      ii     iii
C
      data r   /       .58e0,  .62e0,  .67e0,  .77e0,  .78e0 /
      data ad1 /       .35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0 /
      data ad2 /       23.e0,  20.e0,  17.e0,  14.e0,  7.9e0 /
C
C-----------------------------------------------------------------------
C
C     Surface boundary condition:
C
C       nbc   prescribed    prescribed   short wave
C             temperature      flux      penetration
C             or salinity               (temperature
C                                           only)
C
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
C     NOTE that only 1 and 3 are allowed for salinity.
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kh*f')'-f=-fb
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C     Calculate penetrative radiation. At the bottom any unattenuated
C     radiation is deposited in the bottom layer:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rad(i,j,k)=0.e0
          end do
        end do
      end do
C
      if(nbc.eq.2.or.nbc.eq.4) then
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)
     $                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))
     $                      +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
C
      endif
C
      if(nbc.eq.1) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
C
      else if(nbc.eq.2) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
C
      else if(nbc.eq.3.or.nbc.eq.4) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.e0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do
C
      endif
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do
C
C     Bottom adiabatic boundary condition:
C
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine profu
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of x-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wusurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do
C

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine profv
C **********************************************************************
C                                                                      *
C * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wvsurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine prxy(label,time,a,im,iskp,jm,jskp,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes a horizontal 2-D field.                      *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm
      real a(im,jm)
      real time,scala
      integer iskp,jskp
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do j=1,jm,jskp
          do i=1,im,iskp
            amx=max(abs(a(i,j)),amx)
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe9.2)
C
      do ib=1,im,cols*iskp
C
        ie=ib+(cols-1)*iskp
        if(ie.gt.im) ie=im
C
        if(scala.ge.0.e0) then
          write(6,3) (i,i=ib,ie,iskp)
    3     format(/,2x,24i5,/)
        else
          write(6,4) (i,i=ib,ie,iskp)
    4     format(/,12i10,/)
        endif
C
        do j=1,jm,jskp
          jwr=jm+1-j
          if(scala.ge.0.e0) then
            write(6,5) jwr,(nint(a(i,jwr)/scale),i=ib,ie,iskp)
    5       format(1x,i3,24i5)
          else
            write(6,6) jwr,(a(i,jwr),i=ib,ie,iskp)
    6       format(1x,i2,12(e10.2))
          endif
        end do
C
        write(6,7)
    7   format(//)
C
      end do
C
      return
C
      end
C
      subroutine prxyz(label,time,a,im,iskp,jm,jskp,kb,ko,nko,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes horizontal layers of a 3-D field with        *
C *                integers or floating point numbers.                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                ko .......... 1-D array of k-indices for output     *
C *                nko ......... number of elements in ko              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C *                (NOTE that this combines functions of old prxyz and *
C *                 eprxyz)                                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      real a(im,jm,kb)
      real time,scala
      integer ko(*)
      integer iskp,jskp,nko
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,k,iko,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do iko=1,nko
          k=ko(iko)
          do j=1,jm,jskp
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe9.2)
C
      do iko=1,nko
C
        k=ko(iko)
C
        write(6,3) k
    3   format(3x,/' Layer k = ',i2)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,2x,24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,12i10,/)
          endif
C
          do j=1,jm,jskp
            jwr=jm+1-j
            if(scala.ge.0.e0) then
              write(6,6) jwr,(nint(a(i,jwr,k)/scale),i=ib,ie,iskp)
    6         format(1x,i3,24i5)
            else
              write(6,7) jwr,(a(i,jwr,k),i=ib,ie,iskp)
    7         format(1x,i2,12(e10.2))
            endif
          end do
C
          write(6,8)
    8     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine prxz(label,time,a,im,iskp,jm,kb,jo,njo,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                x- or i-direction .                                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jo .......... 1-D array of j-indices for output     *
C *                njo ......... number of elements in jo              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer jo(*)
      integer iskp,njo
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,k,ijo,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do ijo=1,njo
            j=jo(ijo)
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe9.2)
C
      do ijo=1,njo
C
        j=jo(ijo)
C
        write(6,3) j
    3   format(3x,/' Section j =',i3)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,'    i =  ',24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,'    i =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),i=ib,ie,iskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),i=ib,ie,iskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),i=ib,ie,iskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine pryz(label,time,a,im,jm,jskp,kb,io,nio,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                y- or j-direction.                                  *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                jskp ........ skipping interval for j               *
C *                io .......... 1-D array of i-indices for output     *
C *                nio ......... number of elements in io              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer io(*)
      integer jskp,nio
      character label*(*)
      real amx,scale
      integer i,j,jb,je,k,iio,cols
C
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do j=1,jm,jskp
            do iio=1,nio
              i=io(iio)
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe9.2)
C
      do iio=1,nio
C
        i=io(iio)
C
        write(6,3) i
    3   format(3x,/' Section i =',i3)
C
        do jb=1,jm,cols*jskp
C
          je=jb+(cols-1)*jskp
          if(je.gt.jm) je=jm
C
          if(scala.ge.0.e0) then
            write(6,4) (j,j=jb,je,jskp)
    4       format(/,'    j =  ',24i5,/)
          else
            write(6,5) (j,j=jb,je,jskp)
    5       format(/,'    j =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),j=jb,je,jskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),j=jb,je,jskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),j=jb,je,jskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
 
C
      subroutine slpmax
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Limits the maximum of:                              *
C *                                                                    *
C *                  <difference of depths>/<sum of depths>            *
C *                                                                    *
C *                for two adjacent cells. The maximum possible value  *
C *                is unity.                                           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real mean,del
      integer i,j,loop
C
      do loop=1,10
C
C     Sweep right:
C
        do j=2,jm-1
C
          do i=2,im-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
C    Sweep left:
C
          do i=im-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
        end do
C
C   Sweep up:
C
        do i=2,im-1
C
          do j=2,jm-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
C   Sweep down:
C
          do j=jm-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the antidiffusive velocity used to       *
C *                reduce the numerical diffusion associated with the  *
C *                upstream differencing scheme.                       *
C *                                                                    *
C *                This is based on a subroutine of Gianmaria Sannino  *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option.                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real ff(im,jm,kb)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      real sw
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
C
      parameter (value_min=1.e-9,epsilon=1.0e-14)
C
C     Apply temperature and salinity mask:
C
      do k=1,kb
        do i=1,im
          do j=1,jm
            ff(i,j,k)=ff(i,j,k)*fsm(i,j)
          end do
        end do
      end do
C
C     Recalculate mass fluxes with antidiffusion velocity:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.e0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.e0
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.e0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.e0
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.e0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine vertvl(xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates vertical velocity.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Reestablish boundary conditions:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do
C
C     NOTE that, if one wishes to include freshwater flux, the
C     surface velocity should be set to vflux(i,j). See also
C     change made to 2-D volume conservation equation which
C     calculates elf.
C
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do
C
      return
C
      end
C
c     include 'pom2k.n'                                       ! *netCDF*
C
C     End of source code
C
C-----------------------------------------------------------------------
C
