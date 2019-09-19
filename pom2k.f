       program pom2k
C
C **********************************************************************
C *                                                                    *
C *   The last code change as rcorded in pom2k.change was on           *
C *                                                                    *
C *                     2006-05-03                                     *
C *                                  (adding IC from file)             *
C *                                                                    *
C * FUNCTION    :  This is a version of the three dimensional, time    *
C *                dependent, primitive equation, ocean model          *
C *                developed by Alan Blumberg and George Mellor with   *
C *                subsequent contributions by Leo Oey, Steve Brenner  *
C *                and others. It is now called the Princeton Ocean    *
C *                Model. Two references are:                          *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; Diagnostic and      *
C *                  prognostic numerical circulation studies of the   *
C *                  South Atlantic Bight, J. Geophys. Res. 88,        *
C *                  4579-4592, 1983.                                  *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; A description of a  *
C *                  three-dimensional coastal ocean circulation model,*
C *                  Three-Dimensional Coastal Ocean Models, Coastal   *
C *                  and Estuarine Sciences, 4, N.S. Heaps, ed.,       *
C *                  American Geophysical Union, 1-16, 1987.           *
C *                                                                    *
C *                In subroutine profq the model makes use of the      *
C *                turbulence closure sub-model described in:          *
C *                                                                    *
C *                Mellor, G.L. and T. Yamada; Development of a        *
C *                  turbulence closure model for geophysical fluid    *
C *                  problems, Rev. Geophys. Space Phys., 20, No. 4,   *
C *                  851-875, 1982.                                    *
C *            (note recent profq that includes breaking waves)        *
C *                                                                    *
C *                A user's guide is available:                        *
C *                                                                    *
C *                Mellor, G.L.; User's guide for a three-dimensional, *
C *                  primitive equation, numerical ocean model.        *
C *                  Princeton University Report, 1998.                *
C *                                                                    *
C *                In October 2001, the source code underwent          *
C *                revision by John Hunter of the University of        *
C *                Tasmania. Major aspects of the revision were:       *
C *                                                                    *
C *                (1) The revision was based on pom98 updated to      *
C *                    12/9/2001.                                      *
C *                (2) Declaration of all variables.                   *
C *                (3) Rationalisation of the input of all constants.  *
C *                (4) Modifications to the "printer" output.          *
C *                (5) Output to a netCDF file.                        *
C *                (6) Inclusion of surface freshwater flux.           *
C *                (7) Inclusion of atmospheric pressure.              *
C *                (8) Inclusion of an additional problem to check (6) *
C *                    and (7), above.                                 *
C *                (9) Inclusion of option for Smolarkiewicz           *
C *                    advection scheme.                               *
C *                                                                    *
C *                This revised version is functionally almost         *
C *                equivalent to pom98. The output to device 6 from    *
C *                the "seamount" problem should be almost the same,   *
C *                any differences being due to minor format changes   *
C *                and improvements in rounding.                       *
C *                                                                    *
C *                This revision was helped by the following people:   *
C *                Tal Ezer, Peter Holloway, George Mellor, Rich       *
C *                Signell, Ian Webster, Brian Williams and Emma Young.*
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                  GENERAL NOTES                     *
C *                                                                    *
C *                1. All units are S.I. (M.K.S.) unless otherwise     *
C *                   stated. NOTE that time is in days from the start *
C *                   of the run.                                      *
C *                                                                    *
C *                2. "b", <nothing> and "f" refers to backward,       *
C *                   central and forward time levels.                 *
C *                                                                    *
C *                3. NetCDF output may be used. In order to omit/use  *
C *                   netCDF, comment/uncomment all statements         *
C *                   carrying the comment "*netCDF*" at the end of    *
C *                   the line (or set netcdf_file='nonetcdf')         *
C *                                                                    *
C *                4. NetCDF is version 3. An attempt has been made to *
C *                   conform to the NetCDF Climate and Forecast (CF)  *
C *                   Metadata Conventions, but this may not yet be    *
C *                   complete (see:                                   *
C *                                                                    *
C *          http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html) *
C *                                                                    *
C *                5. In order to use netCDF, the program should be    *
C *                   compiled with the appropriate library. For       *
C *                   example, if using g77, you may need to type:     *
C *                                                                    *
C *                     g77 -o pom2k pom2k.f /usr/lib/libnetcdf.a      *
C *                                                                    *
C *                   You should also have the "include" file of       *
C *                   netCDF subroutines (pom2k.n).                    * 
C *                                                                    *
C *                6. In order to use netCDF, you may need to change   *
C *                   the name of the "include" file in the statement: *
C *                                                                    *
C *                     include '/usr/include/netcdf.inc'              *
C *                                                                    *
C *                   in subroutine write_netcdf                       *
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                SOFTWARE LICENSING                  *
C *                                                                    *
C *                This program is free software; you can redistribute *
C *                it and/or modify it under the terms of the GNU      *
C *                General Public License as published by the Free     *
C *                Software Foundation, either Version 2 of the        *
C *                license, or (at your option) any later version.     *
C *                                                                    *
C *                This program is distributed in the hope that it     *
C *                will be useful, but without any warranty; without   *
C *                even the implied warranty of merchantability or     *
C *                fitness for a particular purpose. See the GNU       *
C *                General Public License for more details.            *
C *                                                                    *
C *                A copy of the GNU General Public License is         *
C *                available at http://www.gnu.org/copyleft/gpl.html   *
C *                or by writing to the Free Software Foundation, Inc.,*
C *                59 Temple Place - Suite 330, Boston, MA 02111, USA. *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
C     New declarations plus ispi,isp2i:
C
      real aam_init,atot
      real cbcmax,cbcmin,darea
      real days,dte2,dvol
      real eaver
      real horcon
      real ispi,isp2i
      real period,prtd1,prtd2
      real saver,smoth,sw,swtch
      real taver,time0
      real vamax,vtot,tsalt
      real z0b
      real tatm,satm
      integer io(100),jo(100),ko(100)
      integer i,iend,iext,imax,ispadv,isplit,iswtch
      integer j,jmax
      integer k
      integer nadv,nbct,nbcs,nitera,nread
      integer iproblem
      integer mdir
      logical lramp
      character*120 netcdf_file
      character*80 dir
C
C***********************************************************************
C
C     source should agree with source_c in pom2k.c and source_n in pom2k.n.
C
      source='pom2k  2006-05-03'
C
      if(source.ne.source_c) then
        write(6,7)
    7   format(/'Incompatible versions of program and include files ',
     $          '..... program terminated'/)
        stop
      endif
C
C***********************************************************************
C
      small=1.e-9           ! Small value
      pi=atan(1.e0)*4.e0    ! PI
C
C***********************************************************************
C     Input of filenames and constants:
C     NOTE that the array sizes im, jm and kb should be set in pom2k.c
C-----------------------------------------------------------------------
      title='Seamount                                   ' ! run's title
C-----------------------------------------------------------------------
C
c      netcdf_file='pom2k.nc'  ! netCDF output file
       netcdf_file='nonetcdf'  ! disable netCDF output
C
C-----------------------------------------------------------------------
C     Problem number:
C     iproblem      problem      initialisation
C                    type          subroutine
C         1        seamount       seamount
C         2        conservation   box
C                  box
C         3        IC from file   file2ic
c      iproblem=1   
c
C       mode                     description
C        2        2-D calculation (bottom stress calculated in advave)
C        3        3-D calculation (bottom stress calculated in profu,v)
C        4        3-D calculation with t and s held fixed
c      mode=3      
c      !describle in params 20190414 gaojie
c
C-----------------------------------------------------------------------
C     Advection scheme:
C      nadv     Advection scheme
C        1       Centred scheme, as originally provide in POM
C        2       Smolarkiewicz iterative upstream scheme, based on
C                subroutines provided by Gianmaria Sannino and Vincenzo
C                Artale
      nadv=1
C
C-----------------------------------------------------------------------
C     Constants for Smolarkiewicz iterative upstream scheme.
C     Number of iterations. This should be in the range 1 - 4. 1 is
C     standard upstream differencing; 3 adds 50% CPU time to POM:
      nitera=2
C
C     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
C     gives smoother solutions with less overshoot when nitera > 1:
      sw=0.5e0
C
C-----------------------------------------------------------------------
C     Index to indicate whether run to start from restart file
C     (nread=0: no restart input file; nread=1: restart input file):
      nread=0
C-----------------------------------------------------------------------
C     Date and time of start of initial run of model in format (i.e.
C     UDUNITS convention)
C
C       YYYY-MM-DD HH:MM:SS <+/->HH:MM
C
C     where "<+/->HH:MM" is the time zone (positive eastwards from
C     Coordinated Universal Time). NOTE that the climatological time
C     axis (i.e. beginning of year zero, which does not exist in the
C     real-world calendar) has been used here. Insert your own date
C     and time as required:
C
      time_start='2000-01-01 00:00:00 +00:00'
C
C-----------------------------------------------------------------------
c
C     External (2-D) time step (secs.) according to CFL:
c      dte    = 6.e0    
c      isplit = 30  
c
c      days   = 10.0
c      prtd1  = 1./48.      ! Initial print interval (days)
c      prtd2  = 1.e0        ! Final print interval (days)
c      swtch  = 1.e0      ! Time to switch from prtd1 to prtd2(days) 
c      iskp=4             ! Printout skip interval in i 
c      jskp=3             ! Printout skip interval in j
c
c      !describle in params 20190414 gaojie
C-----------------------------------------------------------------------
C
c      Reference density (recommended values: 1025 for seawater,
C      1000 for freswater; S.I. units):
c      rhoref=1025.e0
c      tbias=0.e0         ! Temperature bias (deg. C)
c      sbias=0.e0         ! Salinity bias
c      grav=9.806e0       ! gravity constant (S.I. units)
c      kappa=0.4e0        ! von Karman's constant
c      z0b=.01e0          ! Bottom roughness (metres)
c      cbcmin=.0025e0     ! Minimum bottom friction coeff.
c      cbcmax=1.e0        ! Maximum bottom friction coeff.
c      horcon=0.2e0       ! Smagorinsky diffusivity coeff.
c
c      !describle in params 20190414 gaojie
C-----------------------------------------------------------------------
c
C      Logical for inertial ramp (.true. if inertial ramp to be applied
C      to wind stress and baroclinic forcing, otherwise .false.)
       lramp=.false.
c
C-----------------------------------------------------------------------
C
C     Inverse horizontal turbulent Prandtl number
C     (ah/am; dimensionless):
C     NOTE that tprni=0.e0 yields zero horizontal diffusivity!
C
      tprni=.2e0 !changed by gaojie for 2004seamount .2e0
C
C-----------------------------------------------------------------------
C
C     Background viscosity used in subroutines profq, proft, profu and
C     profv (S.I. units):
      umol=2.e-5
C
C-----------------------------------------------------------------------
C
C     Maximum depth used in radiation boundary condition in subroutine
C     bcond (metres):
      hmax=4500.e0
C
C-----------------------------------------------------------------------
C
C     Maximum magnitude of vaf (used in check that essentially tests
C     for CFL violation):
      vmaxl=100.e0
C
C-----------------------------------------------------------------------
C
C     Maximum allowable value of:
C       <difference of depths>/<sum of depths>
C     for two adjacent cells (dimensionless). This is used in subroutine
C     slpmax. If >= 1, then slpmax is not applied:
      slmax=2.e0
C
C-----------------------------------------------------------------------
C
C     Integers defining the number of logarithmic layers at the   
C     surface and bottom (used by subroutine depth). The number of
C     logarithmic layers are kl1-2 at the surface and kb-kl2-1
C     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
      kl1=6
      kl2=kb-2
C
C-----------------------------------------------------------------------
C
C     Water type, used in subroutine proft.
C       ntp    Jerlov water type
C
C        1            i
C        2            ia
C        3            ib
C        4            ii
C        5            iii
C
      ntp=2
C
C-----------------------------------------------------------------------
C
C     Surface temperature boundary condition, used in subroutine proft:
C
C       nbct   prescribed    prescribed   short wave
C              temperature      flux      penetration
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
      nbct=1
C
C-----------------------------------------------------------------------
C
C     Surface salinity boundary condition, used in subroutine proft:
C
C       nbcs   prescribed    prescribed
C               salinity      flux
C
C        1        no           yes
C        3        yes          no
C
C     NOTE that only 1 and 3 are allowed for salinity.
      nbcs=1
C
C-----------------------------------------------------------------------
C
C     Step interval during which external (2-D) mode advective terms are
C     not updated (dimensionless):
      ispadv=5
C
C-----------------------------------------------------------------------
C
C     Constant in temporal filter used to prevent solution splitting
C     (dimensionless):
      smoth=0.10e0
C
C-----------------------------------------------------------------------
C
C     Weight used for surface slope term in external (2-D) dynamic
C     equation (a value of alpha = 0.e0 is perfectly acceptable, but the
C     value, alpha=.225e0 permits a longer time step):
      alpha=0.225e0
C
C-----------------------------------------------------------------------
C     Initial value of aam:
c      aam_init=500.e0
c
c      !describle in params 20190414 gaojie
C-----------------------------------------------------------------------
C
C     End of input of constants
C***********************************************************************
C
C --- Above are the default parameters, alternatively one can 
C --- use parameters from a file created by runscript runpom2k
C
      include 'params'
c     include case,txt 'dir' 
      open(39,file='case.txt') ! the dir is writted by run_pom2k.sh and read here
      read(39,'(a60)')dir
c     write(6,*)dir
      call chop(dir,mdir)
      write(6,*) dir,mdir
c      pause
C
C***********************************************************************
C
C     Calculate some constants:
C
      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2
C
      iend=max0(nint(days*24.e0*3600.e0/dti),2)
      iprint=nint(prtd1*24.e0*3600.e0/dti)
      iswtch=nint(swtch*24.e0*3600.e0/dti)
C
      ispi=1.e0/float(isplit)
      isp2i=1.e0/(2.e0*float(isplit))
C
C-----------------------------------------------------------------------
C
C     Print initial summary:
C
      write(6,'(/,'' source   = '',a40)') source
      write(6,'('' title      = '',a40/)') title
      write(6,'('' iproblem   = '',i10)') iproblem 
      write(6,'('' mode       = '',i10)') mode
      write(6,'('' nadv       = '',i10)') nadv
      write(6,'('' nitera     = '',i10)') nitera
      write(6,'('' sw         = '',f10.4)') sw
      write(6,'('' nread      = '',i10)') nread   
      write(6,'('' dte        = '',f10.2)') dte
      write(6,'('' dti        = '',f10.1)') dti
      write(6,'('' isplit     = '',i10)') isplit 
      write(6,'('' time_start = '',a26)') time_start
      write(6,'('' days       = '',f10.4)') days
      write(6,'('' iend       = '',i10)') iend
      write(6,'('' prtd1      = '',f10.4)') prtd1
      write(6,'('' iprint     = '',i10)') iprint 
      write(6,'('' prtd2      = '',f10.4)') prtd2
      write(6,'('' swtch      = '',f10.2)') swtch 
      write(6,'('' iswtch     = '',i10)') iswtch 
      write(6,'('' iskp, jskp = '',i5'','',i5)') iskp,jskp   
      write(6,'('' lramp      = '',l10)') lramp  
      write(6,'('' rhoref     = '',f10.3)') rhoref
      write(6,'('' tbias      = '',f10.3)') tbias 
      write(6,'('' sbias      = '',f10.3)') sbias 
      write(6,'('' grav       = '',f10.4)') grav  
      write(6,'('' kappa      = '',f10.4)') kappa  
      write(6,'('' z0b        = '',f10.6)') z0b         
      write(6,'('' cbcmin     = '',f10.6)') cbcmin      
      write(6,'('' cbcmax     = '',f10.6)') cbcmax      
      write(6,'('' horcon     = '',f10.3)') horcon      
      write(6,'('' tprni      = '',f10.4)') tprni      
      write(6,'('' umol       = '',f10.4)') umol       
      write(6,'('' hmax       = '',f10.2)') hmax       
      write(6,'('' vmaxl      = '',f10.4)') vmaxl      
      write(6,'('' slmax      = '',f10.4)') slmax      
      write(6,'('' kl1, kl2   = '',i5,'','',i5)') kl1,kl2
      write(6,'('' ntp        = '',i10)') ntp  
      write(6,'('' nbct       = '',i10)') nbct  
      write(6,'('' nbcs       = '',i10)') nbcs  
      write(6,'('' ispadv     = '',i10)') ispadv
      write(6,'('' smoth      = '',f10.4)') smoth     
      write(6,'('' alpha      = '',f10.4)') alpha     
C
C-----------------------------------------------------------------------
C
C     Initialise boundary arrays:
C
      do i=1,im
        vabn(i)=0.e0
        vabs(i)=0.e0
        eln(i)=0.e0
        els(i)=0.e0
        do k=1,kb
          vbn(i,k)=0.e0
          vbs(i,k)=0.e0
          tbn(i,k)=0.e0
          tbs(i,k)=0.e0
          sbn(i,k)=0.e0
          sbs(i,k)=0.e0
        end do
      end do
C
      do j=1,jm
        uabe(j)=0.e0
        uabw(j)=0.e0
        ele(j)=0.e0
        elw(j)=0.e0
        do k=1,kb
          ube(j,k)=0.e0
          ubw(j,k)=0.e0
          tbe(j,k)=0.e0
          tbw(j,k)=0.e0
          sbe(j,k)=0.e0
          sbw(j,k)=0.e0
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Initialise 2-D and 3-D arrays for safety (this may be overwritten
C     later):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=0.e0
          vab(i,j)=0.e0
          elb(i,j)=0.e0
          etb(i,j)=0.e0
          e_atmos(i,j)=0.e0
          vfluxb(i,j)=0.e0
          vfluxf(i,j)=0.e0
          wusurf(i,j)=0.e0
          wvsurf(i,j)=0.e0
          wtsurf(i,j)=0.e0
          wssurf(i,j)=0.e0
          swrad(i,j)=0.e0
          drx2d(i,j)=0.e0
          dry2d(i,j)=0.e0
        end do
      end do
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.e0
            vb(i,j,k)=0.e0
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Set up sigma layers:
C
      if(iproblem.ne.3) call depth
C
C-----------------------------------------------------------------------
C
C     Read in grid data, and initial and lateral boundary conditions:
C
c      if(iproblem.eq.1) then
c
        call seamount
c
c      else if(iproblem.eq.2) then
c        call box
c      else if(iproblem.eq.3) then
c        call file2ic
c      else
c        write(6,8)
c    8   format(/' Invalid value of iproblem ..... program terminated'/)
c        stop
c      endif
C
C     Inertial period for temporal filter:
C
      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0
C
C     Initialise time:
C
      time0=0.e0
      time=0.e0
C
C     Initial conditions:
C
C     NOTE that lateral thermodynamic boundary conditions are often set
C     equal to the initial conditions and are held constant thereafter.
C     Users can of course create variable boundary conditions.
C
      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)
          et(i,j)=etb(i,j)
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            l(i,j,k)=0.1*dt(i,j)
            q2b(i,j,k)=small
            q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
            kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
            km(i,j,k)=kh(i,j,k)
            kq(i,j,k)=kh(i,j,k)
            aam(i,j,k)=aam_init
          end do
        end do
      end do
C
      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do
C
      call dens(s,t,rho)
C
      call baropg
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do
C
C     Calculate bottom friction coefficient:
C
      do j=1,jm
        do i=1,im
          cbc(i,j)=(kappa/log((1.e0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
C
C     If the following is invoked, then it is probable that the wrong
C     choice of z0b or vertical spacing has been made:
C
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
C
C     Calculate external (2-D) CFL time step:
C
      do j=1,jm
        do i=1,im
          tps(i,j)=0.5e0/sqrt(1.e0/dx(i,j)**2+1.e0/dy(i,j)**2)
     $               /sqrt(grav*(h(i,j)+small))*fsm(i,j)
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following data are needed for a seamless restart. if nread=1,
C     data had been created by a previous run (see write(71) at end of
C     this program). nread=0 denotes a first time run.
C
      if(nread.eq.1)
     $  read(70) time0,
     $           wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $           utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,
     $           adx2d,ady2d,advua,advva,
     $           km,kh,kq,l,q2,q2b,aam,q2l,q2lb
C
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do
C
      time=time0
C
C-----------------------------------------------------------------------
C
C     Print geometry and other initial fields (select statements as
C     desired):
C
c      call prxy('grid increment in x, dx                 ',
c     $          time,dx ,im,iskp,jm,jskp,0.e0)
C
c      call prxy('grid increment in y, dy                 ',
c     $          time,dy ,im,iskp,jm,jskp,0.e0)
C
c      call prxy('Easting of elevation points, east_e     ',
c     $          time,east_e ,im,iskp,jm,jskp,0.e0)
C
c      call prxy('Northing of elevation points, north_e   ',
c     $          time,north_e,im,iskp,jm,jskp,0.e0)
C
c      call prxy('Easting of cell corners, east_c         ',
c     $          time,east_c ,im,iskp,jm,jskp,0.e0)
C
c      call prxy('Northing of cell corners, north_c       ',
c     $          time,north_c,im,iskp,jm,jskp,0.e0)
C
c      call prxy('Rotation angle of x-axis wrt. east, rot ',
c     $          time,rot,im,iskp,jm,jskp,0.e0)
C
      call prxy('Undisturbed water depth, h              ',
     $          time,h  ,im,iskp,jm,jskp,1.e1)
C
c      call prxy('Free surface mask, fsm                  ',
c     $          time,fsm,im,iskp,jm,jskp,1.e0)
C
c      call prxy('u-velocity mask, dum                    ',
c     $          time,dum,im,iskp,jm,jskp,1.e0)
C
c      call prxy('v-velocity mask, dvm                    ',
c     $          time,dvm,im,iskp,jm,jskp,1.e0)
C
c      call prxy('External (2-D) CFL time step, tps       ',
c     $          time,tps,im,iskp,jm,jskp,1.e0)
C
C     Set sections for output:
C
      ko(1)=1
      ko(2)=kb/2
      ko(3)=kb-1
C
c      call prxyz('Horizontally-averaged rho, rmean        ',
c     $           time,rmean,im,iskp,jm,jskp,kb,ko,3,1.e-5)
C
C     Set sections for output:
C
      jo(1)=1
      jo(2)=jm/2
      jo(3)=jm-1
C
c      call prxz('Horizontally-averaged rho, rmean        ',
c     $          time,rmean,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
C
C     Set sections for output:
C
      io(1)=1
      io(2)=im/2
      io(3)=im-1
C
c      call pryz('Horizontally-averaged rho, rmean        ',
c     $          time,rmean,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
C
C-----------------------------------------------------------------------
C
C     Initial conditions:
C
C     Select print statements in printall as desired:
C
cc      call printall
C       changed by gaojie 20190415
C-----------------------------------------------------------------------
C
C     Initialise netCDF output and output initial set of data:
C
        if(netcdf_file.ne.'nonetcdf') then
c      call write_netcdf(netcdf_file,1)                        ! *netCDF*
c      call write_netcdf(netcdf_file,2)                        ! *netCDF*
        endif

        open(21,file='boundary_we.dat')
        open(22,file='boundary_sn.dat')

C
C-----------------------------------------------------------------------
C
      do 9000 iint=1,iend      !  Begin internal (3-D) mode
C
        time=dti*float(iint)/86400.e0+time0
C
        if(lramp) then
          ramp=time/period
          if(ramp.gt.1.e0) ramp=1.e0
        else
          ramp=1.e0
        endif
C
C       write(6,2) mode,iint,time
C   2   format(' mode,iint,time =',2i5,f9.2)
C
C-----------------------------------------------------------------------
C
C     Set time dependent, surface and lateral boundary conditions.
C     The latter will be used in subroutine bcond. Users may
C     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
C     wssurf, swrad and vflux.
C
C     Introduce simple wind stress. Value is negative for westerly or
C     southerly winds. The following wind stress has been tapered
C     along the boundary to suppress numerically induced oscilations
C     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
C     To make a healthy surface Ekman layer, it would be well to set
C     kl1=9.
C
        do j=2,jmm1
          do i=2,imm1
c
      if(iproblem.ne.3) then     ! constant wind read in file2ic
c
c           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1))
            wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))
     $                    *.25e0*(dvm(i,j+1)+dvm(i-1,j+1)
     $                          +dvm(i-1,j)+dvm(i,j))
C --- no wind ----
c           wusurf(i,j)=0.e0
            wvsurf(i,j)=0.e0
       endif
            e_atmos(i,j)=0.e0
            vfluxf(i,j)=0.e0
C
C     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
C     the sea surface. See calculation of elf(i,j) below and subroutines
C     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
C     is no net flow across lateral boundaries, the basin volume will be
C     constant; if also vflux(i,j).ne.0, then, for example, the average
C     salinity will change and, unrealistically, so will total salt. 
C
            w(i,j,1)=vfluxf(i,j)
C
C     Set wtsurf to the sensible heat, the latent heat (which involves
C     only the evaporative component of vflux) and the long wave
C     radiation:
C
            wtsurf(i,j)=0.e0
C
C     Set swrad to the short wave radiation:
C
            swrad(i,j)=0.e0
C
C     To account for change in temperature of flow crossing the sea
C     surface (generally quite small compared to latent heat effect)
C
            tatm=t(i,j,1)+tbias    ! an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
C
C     Set the salinity of water vapor/precipitation which enters/leaves
C     the atmosphere (or e.g., an ice cover)
C    
            satm=0.e0              
            wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)  
C
          end do
        end do
C
C-----------------------------------------------------------------------
C
C     Set lateral viscosity:
C
C     If mode=2 then initial values of aam2d are used. If one wishes
C     to use Smagorinsky lateral viscosity and diffusion for an
C     external (2-D) mode calculation, then appropiate code can be
C     adapted from that below and installed just before the end of the
C     "if(mode.eq.2)" loop in subroutine advave.
C
C     Calculate Smagorinsky lateral viscosity:
C
C       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
C                                     +.5*(du/dy+dv/dx)**2) )
C
        if(mode.ne.2) then
          call advct(a,c,ee)
          call baropg
C
          do k=1,kbm1
            do j=2,jmm1
              do i=2,imm1
                aam(i,j,k)=horcon*dx(i,j)*dy(i,j)
     $                      *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                            +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                      +.5e0*(.25e0*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                   -u(i,j-1,k)-u(i+1,j-1,k))
     $                      /dy(i,j)
     $                      +.25e0*(v(i+1,j,k)+v(i+1,j+1,k)
     $                             -v(i-1,j,k)-v(i-1,j+1,k))
     $                      /dx(i,j)) **2)
              end do
            end do
          end do
C
C     Form vertical averages of 3-D fields for use in external (2-D)
C     mode:
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=0.e0
              ady2d(i,j)=0.e0
              drx2d(i,j)=0.e0
              dry2d(i,j)=0.e0
              aam2d(i,j)=0.e0
            end do
          end do
C
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k)
                ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k)
                drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
                dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
                aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k)
              end do
            end do
          end do
C
          call advave(tps)
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)-advua(i,j)
              ady2d(i,j)=ady2d(i,j)-advva(i,j)
            end do
          end do
C
        endif
C
        do j=1,jm
          do i=1,im
            egf(i,j)=el(i,j)*ispi
          end do
        end do
C
        do j=1,jm
          do i=2,im
            utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do
C
C-----------------------------------------------------------------------
C
        do 8000 iext=1,isplit    ! Begin external (2-D) mode 
C
C         write(6,3) iext,time
C   3     format(' iext,time =',i5,f9.2)
C
          do j=2,jm
            do i=2,im
              fluxua(i,j)=.25e0*(d(i,j)+d(i-1,j))
     $                     *(dy(i,j)+dy(i-1,j))*ua(i,j)
              fluxva(i,j)=.25e0*(d(i,j)+d(i,j-1))
     $                     *(dx(i,j)+dx(i,j-1))*va(i,j)
            end do
          end do
C
C     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
C     with pom98.f. See also modifications to subroutine vertvl.
C
          do j=2,jmm1
            do i=2,imm1
              elf(i,j)=elb(i,j)
     $                  +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                          +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                          -vfluxf(i,j))
            end do
          end do
C
          call bcond(1)

          if(mod(iext,ispadv).eq.0) call advave(tps)
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=adx2d(i,j)+advua(i,j)
     $                  -aru(i,j)*.25e0
     $                    *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                     +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $                  +.25e0*grav*(dy(i,j)+dy(i-1,j))
     $                    *(d(i,j)+d(i-1,j))
     $                    *((1.e0-2.e0*alpha)
     $                       *(el(i,j)-el(i-1,j))
     $                      +alpha*(elb(i,j)-elb(i-1,j)
     $                             +elf(i,j)-elf(i-1,j))
     $                      +e_atmos(i,j)-e_atmos(i-1,j))
     $                  +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
            end do
          end do
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                    *aru(i,j)*uab(i,j)
     $                  -4.e0*dte*uaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                     *aru(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=ady2d(i,j)+advva(i,j)
     $                  +arv(i,j)*.25e0
     $                    *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $                     +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $                  +.25e0*grav*(dx(i,j)+dx(i,j-1))
     $                    *(d(i,j)+d(i,j-1))
     $                    *((1.e0-2.e0*alpha)*(el(i,j)-el(i,j-1))
     $                      +alpha*(elb(i,j)-elb(i,j-1)
     $                             +elf(i,j)-elf(i,j-1))
     $                      +e_atmos(i,j)-e_atmos(i,j-1))
     $                  +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                    *vab(i,j)*arv(i,j)
     $                  -4.e0*dte*vaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                     *arv(i,j))
            end do
          end do
C
          call bcond(2)
C
          if(iext.eq.(isplit-2))then
            do j=1,jm
              do i=1,im
                etf(i,j)=.25e0*smoth*elf(i,j)
              end do
            end do
C
          else if(iext.eq.(isplit-1)) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=etf(i,j)+.5e0*(1.-.5e0*smoth)*elf(i,j)
              end do
            end do
C
          else if(iext.eq.isplit) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=(etf(i,j)+.5e0*elf(i,j))*fsm(i,j)
              end do
            end do
C
          endif
C
C     Stop if velocity condition violated (generally due to CFL
C     criterion not being satisfied):
C
          vamax=0.e0
C
          do j=1,jm
            do i=1,im
              if(abs(vaf(i,j)).ge.vamax) then
                vamax=abs(vaf(i,j))
	        imax=i
	        jmax=j
              endif
            end do
          end do
C
          if(vamax.le.vmaxl) then
C
C     Apply filter to remove time split and reset time sequence:
C
            do j=1,jm
              do i=1,im
                ua(i,j)=ua(i,j)
     $                   +.5e0*smoth*(uab(i,j)-2.e0*ua(i,j)+uaf(i,j))
                va(i,j)=va(i,j)
     $                   +.5e0*smoth*(vab(i,j)-2.e0*va(i,j)+vaf(i,j))
                el(i,j)=el(i,j)
     $                   +.5e0*smoth*(elb(i,j)-2.e0*el(i,j)+elf(i,j))
                elb(i,j)=el(i,j)
                el(i,j)=elf(i,j)
                d(i,j)=h(i,j)+el(i,j)
                uab(i,j)=ua(i,j)
                ua(i,j)=uaf(i,j)
                vab(i,j)=va(i,j)
                va(i,j)=vaf(i,j)
              end do
            end do
C
            if(iext.ne.isplit) then
              do j=1,jm
                do i=1,im
                  egf(i,j)=egf(i,j)+el(i,j)*ispi
                end do
              end do
              do j=1,jm
                do i=2,im
                  utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
                end do
              end do
              do j=2,jm
                do i=1,im
                  vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
                end do
              end do
            endif
C
          endif
C
 8000 continue        ! End of external (2-D) mode
C
C-----------------------------------------------------------------------
C
        if(vamax.le.vmaxl) then
C
C     Continue with internal (3-D) mode calculation:
C
          if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then
C
C     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=2,im
                  u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                     (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=2,jm
                do i=1,im
                  v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                     (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
                end do
              end do
            end do
C
C     vertvl calculates w from u, v, dt (h+et), etf and etb:
C
            call vertvl(a,c)
            call bcond(5)
C
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  uf(i,j,k)=0.e0
                  vf(i,j,k)=0.e0
                end do
              end do
            end do
C
C     Calculate q2f and q2lf using uf, vf, a and c as temporary
C     variables:
C
            call advq(q2b,q2,uf,a,c)
            call advq(q2lb,q2l,vf,a,c)
            call profq(a,c,tps,dtef)
            call bcond(6)
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  q2(i,j,k)=q2(i,j,k)
     $                       +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                    -2.e0*q2(i,j,k))
                  q2l(i,j,k)=q2l(i,j,k)
     $                       +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                    -2.e0*q2l(i,j,k))
                  q2b(i,j,k)=q2(i,j,k)
                  q2(i,j,k)=uf(i,j,k)
                  q2lb(i,j,k)=q2l(i,j,k)
                  q2l(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
C
C     Calculate tf and sf using uf, vf, a and c as temporary variables:
C
            if(mode.ne.4) then
C
              if(nadv.eq.1) then
C
                call advt1(tb,t,tclim,uf,a,c)
                call advt1(sb,s,sclim,vf,a,c)
C
              else if(nadv.eq.2) then
C
                call advt2(tb,t,tclim,uf,a,c,nitera,sw)
                call advt2(sb,s,sclim,vf,a,c,nitera,sw)
C
              else
C
                write(6,9)
    9           format(/'Invalid value for nadv ..... ',
     $                 'program terminated'/)
                stop
C
              endif
C
              call proft(uf,wtsurf,tsurf,nbct,tps)
              call proft(vf,wssurf,ssurf,nbcs,tps)
              call bcond(4)
C
              do k=1,kb
                do j=1,jm
                  do i=1,im
                    t(i,j,k)=t(i,j,k)
     $                        +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                     -2.e0*t(i,j,k))
                    s(i,j,k)=s(i,j,k)
     $                        +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                     -2.e0*s(i,j,k))
                    tb(i,j,k)=t(i,j,k)
                    t(i,j,k)=uf(i,j,k)
                    sb(i,j,k)=s(i,j,k)
                    s(i,j,k)=vf(i,j,k)
                  end do
                end do
              end do
C
              call dens(s,t,rho)
C
            endif
C
C     Calculate uf and vf:
C
            call advu
            call advv
            call profu
            call profv
            call bcond(3)
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  u(i,j,k)=u(i,j,k)
     $                      +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                                   -2.e0*u(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  v(i,j,k)=v(i,j,k)
     $                      +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                                   -2.e0*v(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  ub(i,j,k)=u(i,j,k)
                  u(i,j,k)=uf(i,j,k)
                  vb(i,j,k)=v(i,j,k)
                  v(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
C
          endif
C
          do j=1,jm
            do i=1,im
              egb(i,j)=egf(i,j)
              etb(i,j)=et(i,j)
              et(i,j)=etf(i,j)
              dt(i,j)=h(i,j)+et(i,j)
              utb(i,j)=utf(i,j)
              vtb(i,j)=vtf(i,j)
              vfluxb(i,j)=vfluxf(i,j)
            end do
          end do
C
        endif
C
c            
c       check the boundary condition 
       
      do j=1,jm
c         do k=1,kb
c            write(21,"(4(1x,e11.4))")tbw(j,k),tbe(j,k),ubw(j,k),ube(j,k)
c         end do 
         write(21,"(2(1x,e11.4))")uabw(j),uabe(j)
      end do 
     
      do i=1,im
c         do k=1,kb
c            write(22,"(4(1x,e11.4))")tbs(i,k),tbn(i,k),vbs(i,k),vbn(i,k)
c         end do 
        write(22,"(2(1x,e11.4))")vabs(i),vabn(i)
      end do 
c        write(6,*)'t', iint,t(1,1,1),t(10,10,1),t(10,10,3)
C-----------------------------------------------------------------------
C
C     Beginning of print section:
C
        if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
C
        if(mod(iint,iprint).eq.0.or.vamax.gt.vmaxl) then
C
          write(6,4) time,iint,iext,iprint
    4     format(/
     $    '**************************************************',
     $    '**************************************************',
     $    '*************************'//
     $    ' time =',f9.4,', iint =',i8,', iext =',i8,', iprint =',i8,//)
C
C     Select print statements in printall as desired:
C
c          call printall
c
          call snapout(dir,mdir) ! print all data to dir
c         changed by gaojie 2019/04/15
C
          vtot=0.e0
          atot=0.e0
          taver=0.e0
          saver=0.e0
          eaver=0.e0
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                darea=dx(i,j)*dy(i,j)*fsm(i,j)
                dvol=darea*dt(i,j)*dz(k)
                vtot=vtot+dvol
                taver=taver+tb(i,j,k)*dvol
                saver=saver+sb(i,j,k)*dvol
              end do
            end do
          end do
C
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              atot=atot+darea
              eaver=eaver+et(i,j)*darea
            end do
          end do
C
          taver=taver/vtot
          saver=saver/vtot
          eaver=eaver/atot
          tsalt=(saver+sbias)*vtot
C
          write(6,5) vtot,atot,eaver,taver,saver,tsalt
    5     format('vtot = ',e16.7,'   atot = ',e16.7,
     $           '  eaver =',e16.7/'taver =',e16.7,
     $           '   saver =',e16.7,'  tsalt =',e16.7)  
C
C     Write netCDF output:
C
            if(netcdf_file.ne.'nonetcdf') then
c         call write_netcdf(netcdf_file,2)                    ! *netCDF*
            endif
C
          if(vamax.gt.vmaxl) then
C
            write(6,4) time,iint,iext,iprint
C
c            call printall
C
            write(6,6) vamax,imax,jmax
    6       format(///////////////////
     $             '************************************************'/
     $             '************ abnormal job end ******************'/
     $             '************* user terminated ******************'/
     $             '************************************************'/
     $             ' vamax =',e12.3,'   imax,jmax =',2i5)
C
C     Close netCDF file:
C
              if(netcdf_file.ne.'nonetcdf') then
c           call write_netcdf(netcdf_file,3)                  ! *netCDF*
              endif
C
            stop
C
          endif
C
        endif
C
C     End of print section
C
C-----------------------------------------------------------------------
C
 9000 continue       !  End of internal (3-D) mode
          close(21)
          close(22)
 
C
C-----------------------------------------------------------------------
C
      write(6,4) time,iint,iext,iprint
C
C     Set levels for output:
C
      ko(1)=1
      ko(2)=2
      ko(3)=kb/2
      ko(4)=kb-1
      ko(5)=kb
C
C     call prxyz('Vertical velocity, w                    ',
C    $           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
C
C     call prxyz('Turbulent kinetic energy x 2, q2        ',
C    $           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
C
C     Save this data for a seamless restart:
C
      write(71) time,
     $  wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $  utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,
     $  km,kh,kq,l,q2,q2b,aam,q2l,q2lb
C
C     Close netCDF file:
C
        if(netcdf_file.ne.'nonetcdf') then
c       call write_netcdf(netcdf_file,3)                        ! *netCDF*
        endif
C
      stop
C
      end
C
C     End of main program
C
C-----------------------------------------------------------------------
