      subroutine seamount
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up for seamount problem.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delh,delx,ra,vel,H0
      real elejmid,elwjmid,elsimid,elnimid
      integer i,j,k
      integer mend,mdir
      character dir*80

      H0=2000.e0
C
C     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
C
      delh=0.e0
C
C     Grid size:
C
      delx=8000.e0
C
C     Radius island or seamount:
C
      ra=25000.e0
C
C     Current velocity:
C
      vel=0.2e0
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im
C
C     For constant grid size:
C
         dx(i,j)=delx
         dy(i,j)=delx
C
C     For variable grid size:
C
c          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
c          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
C
          cor(i,j)=1.e-4
C         we'd better not change cor to 0.e0 20190918 gaoj
        end do
      end do
C
C     Calculate horizontal coordinates of grid points and rotation
C     angle.
C
C     NOTE that this is introduced solely for the benefit of any post-
C     processing software, and in order to conform with the requirements
C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
C
C     There are four horizontal coordinate systems, denoted by the
C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
C     "e" is an elevation point and "c" is a cell corner), as shown
C     below. In addition, "east_*" is an easting and "north_*" is a
C     northing. Hence the coordinates of the "u" points are given by
C     (east_u,north_u).
C
C     Also, if the centre point of the cell shown below is at
C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
C     the coordinates of the western of the two "u" points,
C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
C     coordinates of the southwestern corner point of the cell. The
C     southwest corner of the entire grid is at
C     (east_c(1,1),north_c(1,1)).
C
C                      |              |
C                    --c------v-------c--
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                      u      e       u
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                    --c------v-------c--
C                      |              |
C
C
C     NOTE that the following calculation of east_c and north_c only
C     works properly for a rectangular grid with east and north aligned
C     with i and j, respectively:
C
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
C
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
C
C     The following works properly for any grid:
C
C     Elevation points:
C
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
     $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
     $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im-1
        east_e(i,jm)
     $    =((east_c(i,jm)+east_c(i+1,jm))*3.e0
     $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)
     $    =((north_c(i,jm)+north_c(i+1,jm))*3.e0
     $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
C
      do j=1,jm-1
        east_e(im,j)
     $    =((east_c(im,j)+east_c(im,j+1))*3.e0
     $       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)
     $    =((north_c(im,j)+north_c(im,j+1))*3.e0
     $       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
C
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
     $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
     $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre:
C
C     (NOTE that the following calculation of rot only works properly
C     for this particular rectangular grid)
C
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
C
          h(i,j)=H0*(1.e0-delh)
c     $                          *exp(-((east_c(i,j)
c     $                                   -east_c((im+1)/2,j))**2
c     $                                +(north_c(i,j)
c     $                                   -north_c(i,(jm+1)/2))**2)
c     $                                /ra**2))
          if(h(i,j).lt.1.e0) h(i,j)=1.e0
C
        end do
      end do
      write(6,*)h(1,1)
      pause
C
C     Close the north and south boundaries to form a channel:
C
c      do i=1,im
c        h(i,1)=1.e0
c        h(i,jm)=1.e0
c      end do
C
C     Calculate areas and masks:
C
      mend=20190414

c         open(70,file=dir(1:mdir)//'basic.dat',
c     1        form='unformatted')
         open(70,file='basic.dat',form='unformatted')

            write(70)h,mend
         close(70)
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.e0) call slpmax
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
            vb(i,j,k)=0.e0*dvm(i,j)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
          vab(i,j)=0.e0*dvm(i,j)
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
c      do j=1,jm
c        do i=1,im
C     No conditions necessary for this problem
c        end do
c      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in the seamount problem, the east and west boundaries are open,
C     while the south and north boundaries are closed through the
C     specification of the masks fsm, dum and dvm):
C
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
C
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)
C     Set geostrophically conditioned elevations at the boundaries:
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
C
C     Adjust boundary elevations so that they are zero in the middle
C     of the channel:
C
       elejmid=ele(jmm1/2)
       elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do
c
c     the 2-d boundary condition of v  
cc     in South and North boundary        
       do i=2,imm1
         vabs(i)=vab(i,2)
         vabn(i)=vab(i,jmm1)
         eln(i)=eln(i-1)-cor(i,jmm1)*vab(i,jmm1)/grav*dx(i-1,jmm1)
         els(i)=els(i-1)-cor(i,   2)*vab(i,   2)/grav*dx(i-1,   2)
      end do
C
C     Adjust boundary elevations so that they are zero in the middle
C     of the channel:
C
      elsimid=els(imm1/2)
      elnimid=eln(imm1/2)
      do i=2,imm1
        eln(i)=(eln(i)-elnimid)*fsm(i,jm)
        els(j)=(els(j)-elsimid)*fsm(i, 2)
      end do
c


C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
