!Copyright (c) 2017 Lingping Zeng (lpzeng2@gmail.com) , Samuel Huberman (schuberm@mit.edu)

!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:

!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.

!This program solves the quasi-ballistic phonon transport to mimic the transport
!in the transient grating experiments. Spectral properties are used as input into the simulation.
!Boundary conditions are specified to be diffuse adiabactic and periodic boundaries.
!Key parameter: phonon bundle energy

!Deviational energy in simulation domain is first calculated and then divided by
!the number of phonon bundles to be simulated, yielding the phonon bundle energy.

module Dispersion
    implicit none
    integer,parameter:: Nb = 16611
    double precision:: Si_omega(Nb), Si_Vg(Nb), Si_tau(Nb), Si_DOS(Nb), Si_C_omega(Nb), Si_Dom(Nb) !branch-dependent properties
    double precision:: CDF_C(Nb), CDF_tau(Nb), CDF_CV(Nb)
    double precision:: de_dT(Nb), C_om(Nb), C_tau(Nb), Cv(Nb)
    double precision:: SikTC, SiCv !bulk properties
    integer:: mode
end module Dispersion

module Const
    implicit none
    double precision,parameter:: pi = 4.d0*datan(1.d0)
    double precision,parameter:: kB = 1.38d-23
    double precision,parameter:: e = 1.6d-19
    double precision,parameter:: hbar = 1.05d-34
end module Const

module geometry
    implicit none
    double precision:: lambda  ! grating period
end module geometry

program MC
use Dispersion
use Const
use geometry
implicit none

!Parameter definition
integer:: j                          !branch: polarization indicator
double precision:: Vg, Vg0, Sitau, Freq, PhnSign   !PhnSign: particle sign
double precision:: x0, x, kx0, kx  !position and traveling direction
double precision:: Teq
double precision:: Engy, EngyDev, SiEngyDen !Engy: phonon bundle energy
double precision:: factor, factor1, N_Bose
double precision:: Rx, Rt, Rtheta, Rphi
double precision:: dt, t0, tnew, TotalTime
logical:: Flag_A                         !Boundary or interface scattering indicator
integer:: NParticle, ParNum              !number of phonon particles to be simulated
integer,parameter:: NumData=400
double precision:: PeakTemp(NumData), PeakPseudo(NumData), watch(NumData)
double precision:: XIntp, SurfThickness    !XIntp: Interpolated x position
double precision:: kth

!for infinity
real:: log0

!TEMPERATURE
Teq = 300.d0
SurfThickness = 5.d0 !surface layer thickness at the grating peak and valley, in nm

!Load in dispersion and lifetime data
open(16,file='dataSiSheng.txt',action='read')
CDF_C = 0.d0
CDF_CV = 0.d0
CDF_tau = 0.d0
kth = 0.d0
do j = 1, Nb
    read(16,*) Si_Omega(j),Si_DOS(j),Si_Vg(j),Si_C_omega(j), &
              Si_tau(j)
    de_dT(j) = hbar*Si_omega(j)/Teq*(hbar*Si_omega(j)/Teq)/kB &
                    *dexp(hbar*Si_omega(j)/kB/Teq)/(dexp(hbar*Si_omega(j)/kB/Teq)-1)&
                    /(dexp(hbar*Si_omega(j)/kB/Teq)-1)
    if (j==1) then
        Si_Dom(j) = Si_Omega(j)
    else
        Si_Dom(j) = Si_Omega(j) - Si_Omega(j-1)
    endif
    C_om(j) = Si_DOS(j)*Si_Dom(j)*de_dT(j)
    C_tau(j) = Si_DOS(j)*Si_Dom(j)*de_dT(j)/Si_tau(j)
    CV(j) = Si_DOS(j)*Si_Dom(j)*de_dT(j)*Si_Vg(j)
    if (j==1) then
        CDF_C(j) = C_om(j)
        CDF_CV(j) = CV(j)
        CDF_tau(j) = C_tau(j)
        kth = CV(j)*Si_Vg(j)*Si_tau(j)
    else
        CDF_C(j) = CDF_C(j-1) + C_om(j)
        CDF_CV(j) = CDF_CV(j-1) + CV(j)
        CDF_tau(j) = CDF_tau(j-1) + C_tau(j)
        kth = kth + CV(j)*Si_Vg(j)*Si_tau(j)
    endif
enddo
close(16)

!Normalize
CDF_C = CDF_C/maxval(CDF_C)
CDF_CV = CDF_CV/maxval(CDF_CV)
CDF_tau = CDF_tau/maxval(CDF_tau)

Si_Vg = Si_Vg*1.d9
SikTC = kth
SiCv = sum(C_om)
SiEngyDen = SiCv*Teq
NParticle = 20000000   !Maximum bundles: Consistent with phonon bundle Engy
lambda = 100          !grating period: nm

!compute the effective energy per particle
EngyDev = SiCv*lambda/2.d0/pi*1.d-9
Engy = EngyDev/(NParticle/2.d0)           !the factor of 2 accounts for half domain

open(16,file='PeakTemp.txt')
PeakTemp = 0.d0
XIntp = 0.d0

totaltime = 2.0d-10            !total simulation time
dt=totaltime/dble(NumData)
do j=1, NumData
   watch(j) = dt*j
enddo

call RandGen(10,Rx)

!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS!
!$OMP PARALLEL DO PRIVATE(mode,Rtheta,Rphi,Rx,kx,x,Freq,Vg,Sitau, tnew, t0, Flag_A, x0, kx0, Vg0, j, XIntp, PhnSign)&
!$OMP& REDUCTION(+:PeakTemp)
do ParNum = 1, NParticle   !Beginning of particle loop
    !Initialize a particle property
    !____travelling direction____!
    call random_number(Rtheta)
    call random_number(Rphi)
    Rtheta = 2.d0*Rtheta-1.d0
    kx = sqrt(1.d0-(Rtheta**2.d0))*cos(2.d0*pi*Rphi)

    !____initialize position____!
    call random_number(Rx)
    if (ParNum<=NParticle/2) then
        x = lambda/(2.d0*pi)*dasin(Rx) !unit in nm
        PhnSign = 1.d0
    else
        x = lambda/(2.d0*pi)*(dacos(1.d0-Rx)+pi/2.d0)
        PhnSign = -1.d0
    endif

    !Initialize mode
    mode = 1
    call choose(CDF_C,mode,Nb)
    Freq = Si_omega(mode)
    Vg = Si_Vg(mode)
    Sitau = Si_tau(mode)

    !Initialize t0
    t0 = 0.d0
    tnew = 0.d0
    Flag_A = .FALSE.

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
    do while (tnew<totaltime)
    !MOVE Phonons & boundary and interface scattering
    call random_number(Rt)
    dt = -dlog(dble(Rt))*Sitau
    x0 = x
    kx0 = kx
    Vg0 = Vg
    call PhnMove(dt,x,kx,Freq,Vg,Flag_A,tnew)

    !sample the contribution to surface energy
    do j = 1, NumData
        if (t0<watch(j) .AND. tnew>=watch(j)) then
            XIntp=x0+kx0*Vg0*(watch(j)-t0)
            if (XIntp<=SurfThickness) then
                PeakTemp(j) = PeakTemp(j)+Engy*PhnSign
                PeakPseudo(j) = PeakPseudo(j)+Engy*PhnSign/Sitau
            endif
        endif
    enddo

    !Internal scattering if no interface or boundary scattering occurs
    if (.NOT. Flag_A) then
        call choose(CDF_tau,mode,Nb)
        Freq = Si_omega(mode)
        Vg = Si_Vg(mode)
        Sitau = Si_tau(mode)
        call random_number(Rtheta)
        call random_number(Rphi)
        Rtheta = 2.d0*Rtheta-1.d0
        kx = sqrt(1.d0-Rtheta**2)*cos(2.d0*pi*Rphi)
    endif

    t0 = tnew   !update travelling time

    enddo
    !!!!!!!!!!!!!! End of time loop !!!!!!!!!!!!!!!!!!
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!

enddo ! End of the particle loop
!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE!
!$OMP END PARALLEL DO

!Output the peak-valley temperature difference.
PeakTemp=PeakTemp/(SurfThickness)/SiCv*1.d9

do j = 1, NumData
   write(16,*) watch(j)*1.d12, PeakTemp(j)
enddo

close(16)
stop
end

!!!!!!!!!!!! Subroutines Used !!!!!!!!!!!!!!!!!!!
subroutine RandGen(N,R)
    implicit none
    integer N, j
    double precision R

    do j = 1, N
       call random_number(R)
    enddo

    return
end subroutine

!!!!!!!!!!!!! Move Si Phonon !!!!!!!!!!!!!!!!!!!
recursive subroutine PhnMove(time,x,kx,freq,Vg,Flag_A,tnew) !!!! Change interface and boundary scattering here
    use Geometry
    use Const
    use Dispersion
    implicit none

    integer xyz
    logical Flag_A, Flag_XColl
    double precision time, x, kx, freq, Vg
    double precision XDis, TrueDis
    double precision fx, tmp2, tmp3
    double precision XColl
    double precision:: t0, tnew
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if (time<=0.d0) then
                return ! no need to move if time is zero
    endif
    Flag_A = .false.
    Flag_XColl = .false.
    XDis = 9.9d9
    fx = Vg*time*kx+x
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Check boundary scatterings
    if (kx<0.d0) then ! check if the phonon collides with x boundaries
        call BottomScatter(0.d0,x,fx,Flag_XColl,XDis,XColl)
        if(Flag_XColl) then
            Flag_A =.true.
            Truedis = XDis
            time = dabs((x-XColl)/(Vg*kx))
            kx = -kx
            x = XColl
            tnew = tnew +time
            return
        endif
    elseif (kx>0.d0) then
        call TopScatter(lambda/2.d0,x,fx,Flag_XColl,XDis,XColl)
        if(Flag_XColl) then
            Flag_A =.true.
            Truedis = XDis
            time = dabs((x-XColl)/(Vg*kx))
            kx = -kx
            x = XColl
            tnew = tnew +time
            return
        endif
    endif

    tnew=tnew+time
    x=fx
    return
end subroutine

!!!!!!!!!!!!!!! Interface Scattering !!!!!!!!!!!!!!!!!!!!!!!!
subroutine TopScatter(Intface,sz,fz,Flag_Coll,Dis,ZColl)
    use Geometry
    use Const
    implicit none

    double precision sx,sy,sz,fx,fy,fz,Intface
    double precision Length1,Length2
    double precision Dis,XColl,YColl,ZColl
    logical Flag_Coll
    double precision tmp1,tmp11,tmp2,tmp3,tmpx,tmpy,tmpdis

    Flag_Coll = .false.
    dis = 9.9d9
    if(sz<=Intface .and. fz>=Intface) then
       tmpdis = Intface-sz
       Dis = tmpdis
       Flag_Coll = .true.
       ZColl = Intface
    endif

    return
end subroutine

subroutine BottomScatter(Intface,sz,fz,Flag_Coll,Dis,ZColl)
    use Geometry
    use Const

    implicit none
    double precision sz,fz,Intface
    double precision Dis,ZColl
    logical Flag_Coll
    double precision tmp1,tmp11,tmp2,tmp3,tmpx,tmpy,tmpdis

    Flag_Coll = .false.
    dis = 9.9d9
    if(sz>=Intface .and. fz<=Intface) then 
       tmpdis = dabs(Intface-sz)
       Dis = tmpdis
       Flag_Coll = .true.
       ZColl = Intface
    endif
    return
end subroutine

subroutine choose(ArrayIn,idxOut,n)
    integer:: idxOut
    double precision:: ArrayIn(*)
    double precision:: r
    integer:: i1,i2,i3

    call random_number(r)
    i1 = 1
    i2 = n/2
    i3 = n
    do while(i2-i1>0)
        if(r < ArrayIn(i2)) then
            i3 = i2
            i2 = (i2+i1)/2
        else
            i1 = i2
            i2 = (i2+i3)/2
        endif
    enddo
    idxOut = i2
    return
end subroutine
