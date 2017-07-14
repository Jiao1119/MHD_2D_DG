    module shared_data
    implicit none
    save

    integer, parameter::nphi=6,ngua=3,nequ=8
    integer::iter,istop,istart,nx,nx1,ny,ny1,nn,nn1,ktime,ktmax
    real, parameter::pi=atan(1.0)*4.,cflc=0.18
    real::t,dt,tn,amax,bmax,alphamax,betamax,gamma,gp,gamma1,dx,dy,errh1old,errh1,errh2old,errh2,errhinfold,errhinf
    real::rhoerrh1old,rhoerrh1,rhoerrh2old,rhoerrh2,rhoerrhinfold,rhoerrhinf
    real::uxerrh1old,uxerrh1,uxerrh2old,uxerrh2,uxerrhinfold,uxerrhinf
    real::uyerrh1old,uyerrh1,uyerrh2old,uyerrh2,uyerrhinfold,uyerrhinf
    real::uzerrh1old,uzerrh1,uzerrh2old,uzerrh2,uzerrhinfold,uzerrhinf
    real::bxerrh1old,bxerrh1,bxerrh2old,bxerrh2,bxerrhinfold,bxerrhinf
    real::byerrh1old,byerrh1,byerrh2old,byerrh2,byerrhinfold,byerrhinf
    real::bzerrh1old,bzerrh1,bzerrh2old,bzerrh2,bzerrhinfold,bzerrhinf
    real::perrh1old,perrh1,perrh2old,perrh2,perrhinfold,perrhinf
    real::diverrh1old,diverrh1,diverrh2old,diverrh2,diverrhinfold,diverrhinf
    real::alp(ngua),zg(ngua)
    real::zb(2)
    real::phip(nphi)
    real::phixb(nphi,2,ngua)
    real::phiyb(nphi,ngua,2)
    real::phig(nphi,ngua,ngua)
    real::d1phig(nphi,ngua,ngua),d2phig(nphi,ngua,ngua)
    real, allocatable::x(:),y(:)
    real, allocatable::res(:,:,:,:),resp1(:,:,:,:),resp2(:,:,:,:),resq1(:,:,:,:),resq2(:,:,:,:)
    real, allocatable::u(:,:,:,:,:),p1(:,:,:,:,:),p2(:,:,:,:,:),q1(:,:,:,:,:),q2(:,:,:,:,:)
    real, allocatable::f1(:,:,:),f2(:,:,:),g1(:,:,:),g2(:,:,:)

    end module shared_data

    program main
    use shared_data
    !include 'com.txt'

    errh1old=1.
    errh2old=1.
    errhinfold=1.
    rhoerrh1old=1.
    rhoerrh2old=1.
    rhoerrhinfold=1.
    uxerrh1old=1.
    uxerrh2old=1.
    uxerrhinfold=1.
    uyerrh1old=1.
    uyerrh2old=1.
    uyerrhinfold=1.
    uzerrh1old=1.
    uzerrh2old=1.
    uzerrhinfold=1.
    bxerrh1old=1.
    bxerrh2old=1.
    bxerrhinfold=1.
    byerrh1old=1.
    byerrh2old=1.
    byerrhinfold=1.
    bzerrh1old=1.
    bzerrh2old=1.
    bzerrhinfold=1.
    perrh1old=1.
    perrh2old=1.
    perrhinfold=1.
    diverrh1old=1.
    diverrh2old=1.
    diverrhinfold=1.

    call initialdata

    tn=0.0

    !write(*,*) "input the terminal time "
    !
    !read(*,*) t


    do ktime=1,47
        !    do k=5,7
        nx=64
        nx1=nx+1
        ny=nx
        ny1=ny+1

        nn=nx+10
        nn1=nn+1

        allocate(x(0:nn1+1))
        allocate(y(0:nn1+1))
        allocate(res(nphi,nequ,0:nn1,0:nn1))
        allocate(u(nphi,0:nn1,0:nn1,nequ,0:2))
        allocate(p1(nphi,0:nn1,0:nn1,nequ,0:2))
        allocate(p2(nphi,0:nn1,0:nn1,nequ,0:2))
        allocate(q1(nphi,0:nn1,0:nn1,nequ,0:2))
        allocate(q2(nphi,0:nn1,0:nn1,nequ,0:2))

        tn=0.0
        ktmax=1000000

        t=real(ktime)
        !t=5.

        call init

        io=0
        istop=0
        iter=0


1000    continue

        call setbc(io)

        call setdt

        if(tn+dt>t .and. tn<=t) then
            dt=t-tn
            istop=1
        end if

        if (iter>ktmax) istop=1
        if (tn<t .or. istop .ne. 1) then

101         continue

            call setbc(io)

            call fx(io)

            call rk(io)

            io=mod(io+1,3)

            if (io .eq. 1) tn=tn+dt
            if (io .eq. 2) tn=tn-.5*dt
            if (io .eq. 0) tn=tn+.5*dt
            if (io .ne. 0) go to 101

        end if

        iter=iter+1
        write(06,*) nx,ny,iter,tn,dt

        if(io.eq.0 .and. istop.eq.1) goto 1001
        goto 1000
1001    continue

        write(*,*) "rkdg is completed!"

        call outp

        deallocate(x)
        deallocate(y)
        deallocate(res)
        deallocate(u)
        deallocate(p1)
        deallocate(p2)
        deallocate(q1)
        deallocate(q2)

        !    end do
    end do

    end program main

    subroutine initialdata
    use shared_data
    !include 'com.txt'

    !pi=atan(1.0)*4.
    gamma=5./3.
    gamma1=gamma-1.
    !cflc=0.1

    !--------------------------------------------------------the quadrature rule
    gp=sqrt(0.6)*0.5
    zg(1)=-gp
    zg(2)=0.
    zg(3)=-zg(1)
    alp(1)=5./18.
    alp(2)=8./18.
    alp(3)=5./18.

    zb(1)=-1./2
    zb(2)=1./2
    !----------------------------------------------------- orthogonal basis function
    call basisf

    !---------------------------------------------------set error

    return

    end subroutine initialdata


    subroutine basisf
    use shared_data
    !include 'com.txt'

    !----------------------------------------------basis function
    phi1(zx,zy)=1.
    phi2(zx,zy)=sqrt(12.)*zx
    phi3(zx,zy)=sqrt(12.)*zy
    phi4(zx,zy)=12.*zy*zx
    phi5(zx,zy)=sqrt(5.)*(12.*zx*zx-1.)/2.
    phi6(zx,zy)=sqrt(5.)*(12.*zy*zy-1.)/2.

    d1phi1(zx,zy)=0.
    d1phi2(zx,zy)=sqrt(12.)
    d1phi3(zx,zy)=0.
    d1phi4(zx,zy)=12.*zy
    d1phi5(zx,zy)=sqrt(5.)*(12.*2.*zx)/2.
    d1phi6(zx,zy)=0.

    d2phi1(zx,zy)=0.
    d2phi2(zx,zy)=0.
    d2phi3(zx,zy)=sqrt(12.)
    d2phi4(zx,zy)=12.*zx
    d2phi5(zx,zy)=0.
    d2phi6(zx,zy)=sqrt(5.)*(12.*2.*zy)/2.
    !----------------------------------------------phi at the boundary points of the zx=-1,1,zy=yg
    do 11 ixb=1,2
        do 11 jyg=1,ngua
            phixb(1,ixb,jyg)=phi1(zb(ixb),zg(jyg))
            phixb(2,ixb,jyg)=phi2(zb(ixb),zg(jyg))
            phixb(3,ixb,jyg)=phi3(zb(ixb),zg(jyg))
            phixb(4,ixb,jyg)=phi4(zb(ixb),zg(jyg))
            phixb(5,ixb,jyg)=phi5(zb(ixb),zg(jyg))
            phixb(6,ixb,jyg)=phi6(zb(ixb),zg(jyg))
11  continue
    !----------------------------------------------phi at the boundary points of the zy=-1,1,zx=xg
    do 12 ixg=1,ngua
        do 12 jyb=1,2
            phiyb(1,ixg,jyb)=phi1(zg(ixg),zb(jyb))
            phiyb(2,ixg,jyb)=phi2(zg(ixg),zb(jyb))
            phiyb(3,ixg,jyb)=phi3(zg(ixg),zb(jyb))
            phiyb(4,ixg,jyb)=phi4(zg(ixg),zb(jyb))
            phiyb(5,ixg,jyb)=phi5(zg(ixg),zb(jyb))
            phiyb(6,ixg,jyb)=phi6(zg(ixg),zb(jyb))
12  continue
    !----------------------------------------------phi at the Gauss points xg,yg
    do 13 jyg=1,ngua
        do 13 ixg=1,ngua
            phig(1,ixg,jyg)=phi1(zg(ixg),zg(jyg))
            phig(2,ixg,jyg)=phi2(zg(ixg),zg(jyg))
            phig(3,ixg,jyg)=phi3(zg(ixg),zg(jyg))
            phig(4,ixg,jyg)=phi4(zg(ixg),zg(jyg))
            phig(5,ixg,jyg)=phi5(zg(ixg),zg(jyg))
            phig(6,ixg,jyg)=phi6(zg(ixg),zg(jyg))
13  continue
    !----------------------------------------------phi at the grid points x_i,y_i
    zxi=0.
    zyj=0.
    phip(1)=phi1(zxi,zyj)
    phip(2)=phi2(zxi,zyj)
    phip(3)=phi3(zxi,zyj)
    phip(4)=phi4(zxi,zyj)
    phip(5)=phi5(zxi,zyj)
    phip(6)=phi6(zxi,zyj)
    !----------------------------------------------derivative of phi for x at the Gauss points xg,yg
    do 17 jyg=1,ngua
        do 17 ixg=1,ngua
            d1phig(1,ixg,jyg)=d1phi1(zg(ixg),zg(jyg))
            d1phig(2,ixg,jyg)=d1phi2(zg(ixg),zg(jyg))
            d1phig(3,ixg,jyg)=d1phi3(zg(ixg),zg(jyg))
            d1phig(4,ixg,jyg)=d1phi4(zg(ixg),zg(jyg))
            d1phig(5,ixg,jyg)=d1phi5(zg(ixg),zg(jyg))
            d1phig(6,ixg,jyg)=d1phi6(zg(ixg),zg(jyg))
17  continue
    !----------------------------------------------derivative of phi for y at the Gauss points xg,yg
    do 18 jyg=1,ngua
        do 18 ixg=1,ngua
            d2phig(1,ixg,jyg)=d2phi1(zg(ixg),zg(jyg))
            d2phig(2,ixg,jyg)=d2phi2(zg(ixg),zg(jyg))
            d2phig(3,ixg,jyg)=d2phi3(zg(ixg),zg(jyg))
            d2phig(4,ixg,jyg)=d2phi4(zg(ixg),zg(jyg))
            d2phig(5,ixg,jyg)=d2phi5(zg(ixg),zg(jyg))
            d2phig(6,ixg,jyg)=d2phi6(zg(ixg),zg(jyg))
18  continue

    return

    end subroutine basisf

    subroutine init
    use shared_data
    !include 'com.txt'

    ene(rh,rmx,rmy,rmz,bx,by,bz,p)=p/gamma1+0.5*(rmx*rmx+rmy*rmy+rmz*rmz)/rh+0.5*(bx*bx+by*by+bz*bz)
    ux0(zx,zy)=1.-(zy/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
    uy0(zx,zy)=1.+(zx/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
    bx0(zx,zy)=-(zy/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
    by0(zx,zy)=(zx/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
    p0(zx,zy)=1.+((-zx*zx-zy*zy)*exp(1.-zx*zx-zy*zy))/(8.*pi*pi)

    dx=20./nx
    dy=20./ny

    do i=0,nx1+1
        x(i)=(i-0.5)*dx-10.
    end do

    do j=0,ny1+1
        y(j)=(j-0.5)*dy-10.
    end do

    rh0=1.
    uz0=0.
    bz0=0.


    !!$OMP PARALLEL DO
    do i=1,nx
        do j=1,ny
            do m=1,nphi
                u(m,i,j,1,0)=0.
                u(m,i,j,2,0)=0.
                u(m,i,j,3,0)=0.
                u(m,i,j,4,0)=0.
                u(m,i,j,5,0)=0.
                u(m,i,j,6,0)=0.
                u(m,i,j,7,0)=0.
                u(m,i,j,8,0)=0.
                do ixg=1,ngua
                    do jyg=1,ngua
                        zx=x(i)+dx*zg(ixg)
                        zy=y(j)+dy*zg(jyg)
                        rmx0=rh0*ux0(zx,zy)
                        rmy0=rh0*uy0(zx,zy)
                        rmz0=rh0*uz0
                        ene0=ene(rh0,rmx0,rmy0,rmz0,bx0(zx,zy),by0(zx,zy),bz0,p0(zx,zy))
                        u(m,i,j,1,0)=u(m,i,j,1,0)+alp(ixg)*alp(jyg)*rh0*phig(m,ixg,jyg)
                        u(m,i,j,2,0)=u(m,i,j,2,0)+alp(ixg)*alp(jyg)*rmx0*phig(m,ixg,jyg)
                        u(m,i,j,3,0)=u(m,i,j,3,0)+alp(ixg)*alp(jyg)*rmy0*phig(m,ixg,jyg)
                        u(m,i,j,4,0)=u(m,i,j,4,0)+alp(ixg)*alp(jyg)*rmz0*phig(m,ixg,jyg)
                        u(m,i,j,5,0)=u(m,i,j,5,0)+alp(ixg)*alp(jyg)*bx0(zx,zy)*phig(m,ixg,jyg)
                        u(m,i,j,6,0)=u(m,i,j,6,0)+alp(ixg)*alp(jyg)*by0(zx,zy)*phig(m,ixg,jyg)
                        u(m,i,j,7,0)=u(m,i,j,7,0)+alp(ixg)*alp(jyg)*bz0*phig(m,ixg,jyg)
                        u(m,i,j,8,0)=u(m,i,j,8,0)+alp(ixg)*alp(jyg)*ene0*phig(m,ixg,jyg)
                    end do
                end do
            end do
        end do
    end do

    return

    end subroutine init

    subroutine setdt
    use shared_data
    !include 'com.txt'

    call setamax
    dt=cflc/amax

    return

    end subroutine setdt

    subroutine setamax
    use shared_data
    !include 'com.txt'

    fxpress(rho,rmx,rmy,rmz,b1,b2,b3,ene)=gamma1*(ene-0.5*(rmx*rmx+rmy*rmy+rmz*rmz)/rho-0.5*(b1*b1+b2*b2+b3*b3))

    amax=0.
    !!$OMP PARALLEL DO
    do i=1,nx
        do j=1,ny
            rho=u(1,i,j,1,0)
            velx=abs(u(1,i,j,2,0)/rho)
            vely=abs(u(1,i,j,3,0)/rho)
            press=fxpress(rho,u(1,i,j,2,0),u(1,i,j,3,0),u(1,i,j,4,0),u(1,i,j,5,0),u(1,i,j,6,0),u(1,i,j,7,0),u(1,i,j,8,0))
            ca2=gamma*press/rho
            cb2=(u(1,i,j,5,0)*u(1,i,j,5,0)+u(1,i,j,6,0)*u(1,i,j,6,0)+u(1,i,j,7,0)*u(1,i,j,7,0))/rho
            cx1=ca2+cb2
            cx2=sqrt(abs(cx1*cx1-4.*ca2*u(1,i,j,5,0)*u(1,i,j,5,0)/rho))
            cfx=sqrt(abs(0.5*(cx1+cx2)))
            cy2=sqrt(abs(cx1*cx1-4.*ca2*u(1,i,j,6,0)*u(1,i,j,6,0)/rho))
            cfy=sqrt(abs(0.5*(cx1+cy2)))
            amax=max(amax,(cfx+velx)/dx+(cfy+vely)/dy)
        end do
    end do

    return

    end subroutine setamax

    subroutine setbc(io)
    use shared_data
    !include 'com.txt'

    !!$OMP PARALLEL DO
    do j=0,ny1
        do k=1,nequ
            do m=1,nphi
                u(m,0,j,k,io)=u(m,1,j,k,io)
                u(m,nx1,j,k,io)=u(m,nx,j,k,io)
            end do
        end do
    end do

    !!$OMP PARALLEL DO
    do i=0,nx1
        do k=1,nequ
            do m=1,nphi
                u(m,i,0,k,io)=u(m,i,1,k,io)
                u(m,i,ny1,k,io)=u(m,i,ny,k,io)
            end do
        end do
    end do

    return

    end subroutine setbc

    subroutine fx(io)
    use shared_data
    !include 'com.txt'

    real::rhop(ngua),rmxp(ngua),rmyp(ngua),rmzp(ngua),bxp(ngua),byp(ngua),bzp(ngua),enep(ngua)

    real::rhon(ngua),rmxn(ngua),rmyn(ngua),rmzn(ngua),bxn(ngua),byn(ngua),bzn(ngua),enen(ngua)

    fxpress(rho,rmx,rmy,rmz,bx,by,bz,ene)=gamma1*(ene-0.5*(rmx*rmx+rmy*rmy+rmz*rmz)/rho-0.5*(bx*bx+by*by+bz*bz))

    fxrho(a1,a2,b1,b2,aa)=0.5*(a2+b2)-0.5*aa*(a1-b1)
    fxrmx(w1,w2,w3,w4,w5,w6,w7,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,aa,gamma)=0.5*(a1+b1)*(0.5*(gamma-1.)&
        *((w2/w1)*(w2/w1)+(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1))-(w2/w1)*(w2/w1))&
        +0.5*(a2+b2)*(3.-gamma)*(w2/w1)+0.5*(a3+b3)*(1.-gamma)*(w3/w1)&
        +0.5*(a4+b4)*(1.-gamma)*(w4/w1)-0.5*(a5+b5)*(gamma*w5)&
        +0.5*(a6+b6)*(2.-gamma)*w6+0.5*(a7+b7)*(2.-gamma)*w7&
        +0.5*(a8+b8)*(gamma-1.)-0.5*aa*(a2-b2)
    fxrmy(w1,w2,w3,w5,w6,a1,a2,a3,a5,a6,b1,b2,b3,b5,b6,aa)=-0.5*(a1+b1)*(w2*w3)/(w1*w1)+0.5*(a2+b2)*(w3/w1)+0.5*(a3+b3)*(w2/w1)-0.5*(a5+b5)*w6-0.5*(a6+b6)*w5-0.5*aa*(a3-b3)
    fxrmz(w1,w2,w4,w5,w7,a1,a2,a4,a5,a7,b1,b2,b4,b5,b7,aa)=-0.5*(a1+b1)*(w2*w4)/(w1*w1)+0.5*(a2+b2)*(w4/w1)+0.5*(a4+b4)*(w2/w1)-0.5*(a5+b5)*w7-0.5*(a7+b7)*w5-0.5*aa*(a4-b4)
    fxbx(a5,b5,aa)=-0.5*aa*(a5-b5)
    fxby(w1,w2,w3,w5,w6,a1,a2,a3,a5,a6,b1,b2,b3,b5,b6,aa)=0.5*(a1+b1)*(w3*w5-w2*w6)/(w1*w1)&
        +0.5*(a2+b2)*(w6/w1)-0.5*(a3+b3)*(w5/w1)&
        -0.5*(a5+b5)*(w3/w1)+0.5*(a6+b6)*(w2/w1)&
        -0.5*aa*(a6-b6)
    fxbz(w1,w2,w4,w5,w7,a1,a2,a4,a5,a7,b1,b2,b4,b5,b7,aa)=0.5*(a1+b1)*(w4*w5-w2*w7)/(w1*w1)&
        +0.5*(a2+b2)*(w7/w1)-0.5*(a4+b4)*(w5/w1)&
        -0.5*(a5+b5)*(w4/w1)+0.5*(a7+b7)*(w2/w1)&
        -0.5*aa*(a7-b7)
    fxene(w1,w2,w3,w4,w5,w6,w7,a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,aa,gamma,press)=0.5*(a1+b1)&
        *((-w2*gamma*press)/((gamma-1.)*(w1*w1))&
        +0.5*(w2/w1)*(gamma-2.)*((w2/w1)*(w2/w1)+(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1))&
        +(w5*(w3*w6/w1+w4*w7/w1)-(w2/w1)*(w6*w6+w7*w7))/w1)&
        +0.5*(a2+b2)*(0.5*((3.-2.*gamma)*(w2/w1)*(w2/w1)+(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1)&
        +2.*press*gamma/((gamma-1.)*w1))+(w6*w6+w7*w7)/w1)&
        +0.5*(a3+b3)*((1.-gamma)*(w2*w3)/(w1*w1)-(w5*w6)/w1)&
        +0.5*(a4+b4)*((1.-gamma)*(w2*w4)/(w1*w1)-(w5*w7)/w1)&
        +0.5*(a5+b5)*(-w2*w5*gamma/w1-w3*w6/w1-w4*w7/w1)&
        +0.5*(a6+b6)*(-w3*w5/w1-(gamma-2.)*w2*w6/w1)&
        +0.5*(a7+b7)*(-w4*w5/w1-(gamma-2.)*w2*w7/w1)&
        +0.5*(a8+b8)*(w2*gamma/w1)-0.5*aa*(a8-b8)
    gyrho(c1,c3,d1,d3,bb)=0.5*(c3+d3)-0.5*bb*(c1-d1)
    gyrmx(w1,w2,w3,w5,w6,c1,c2,c3,c5,c6,d1,d2,d3,d5,d6,bb)=-0.5*(c1+d1)*(w2*w3)/(w1*w1)+0.5*(c2+d2)*(w3/w1)+0.5*(c3+d3)*(w2/w1)-0.5*(c5+d5)*w6-0.5*(c6+d6)*w5-0.5*bb*(c2-d2)
    gyrmy(w1,w2,w3,w4,w5,w6,w7,c1,c2,c3,c4,c5,c6,c7,c8,d1,d2,d3,d4,d5,d6,d7,d8,bb,gamma)=0.5*(c1+d1)*(0.5*(gamma-1.)&
        *((w2/w1)*(w2/w1)+(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1))-(w3/w1)*(w3/w1))&
        +0.5*(c2+d2)*(1.-gamma)*(w2/w1)+0.5*(c3+d3)*(3.-gamma)*(w3/w1)&
        +0.5*(c4+d4)*(1.-gamma)*(w4/w1)+0.5*(c5+d5)*(2.-gamma)*w5&
        -0.5*(c6+d6)*gamma*w6+0.5*(c7+d7)*(2.-gamma)*w7&
        +0.5*(c8+d8)*(gamma-1.)-0.5*bb*(c3-d3)
    gyrmz(w1,w3,w4,w6,w7,c1,c3,c4,c6,c7,d1,d3,d4,d6,d7,bb)=-0.5*(c1+d1)*(w3*w4)/(w1*w1)+0.5*(c3+d3)*(w4/w1)+0.5*(c4+d4)*(w3/w1)-0.5*(c6+d6)*w7-0.5*(c7+d7)*w6-0.5*bb*(c4-d4)
    gybx(w1,w2,w3,w5,w6,c1,c2,c3,c5,c6,d1,d2,d3,d5,d6,bb)=0.5*(c1+d1)*(-w3*w5+w2*w6)/(w1*w1)&
        -0.5*(c2+d2)*(w6/w1)+0.5*(c3+d3)*(w5/w1)&
        +0.5*(c5+d5)*(w3/w1)-0.5*(c6+d6)*(w2/w1)&
        -0.5*bb*(c5-d5)
    gyby(c6,d6,bb)=-0.5*bb*(c6-d6)
    gybz(w1,w3,w4,w6,w7,c1,c3,c4,c6,c7,d1,d3,d4,d6,d7,bb)=0.5*(c1+d1)*(w4*w6-w3*w7)/(w1*w1)&
        +0.5*(c3+d3)*(w7/w1)-0.5*(c4+d4)*(w6/w1)&
        -0.5*(c6+d6)*(w4/w1)+0.5*(c7+d7)*(w3/w1)&
        -0.5*bb*(c7-d7)
    gyene(w1,w2,w3,w4,w5,w6,w7,c1,c2,c3,c4,c5,c6,c7,c8,d1,d2,d3,d4,d5,d6,d7,d8,bb,gamma,press)=0.5*(c1+d1)&
        *((-w3*gamma*press)/((gamma-1.)*(w1*w1))&
        +0.5*(w3/w1)*(gamma-2.)*((w2/w1)*(w2/w1)+(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1))&
        +((w2/w1)*w5*w6+(w4/w1)*w6*w7-(w3/w1)*(w5*w5+w7*w7))/w1)&
        +0.5*(c2+d2)*((1.-gamma)*(w2*w3)/(w1*w1)-(w5*w6)/w1)&
        +0.5*(c3+d3)*(0.5*((w2/w1)*(w2/w1)+(3.-2.*gamma)*(w3/w1)*(w3/w1)+(w4/w1)*(w4/w1)&
        +2.*press*gamma/((gamma-1.)*w1))+(w5*w5+w7*w7)/w1)&
        +0.5*(c4+d4)*((1.-gamma)*(w3*w4)/(w1*w1)-(w6*w7)/w1)&
        +0.5*(c5+d5)*(-w2*w6/w1-(gamma-2.)*w3*w5/w1)&
        +0.5*(c6+d6)*(-w3*w6*gamma/w1-w2*w5/w1-w4*w7/w1)&
        +0.5*(c7+d7)*(-w4*w6/w1-(gamma-2.)*w3*w7/w1)&
        +0.5*(c8+d8)*(w3*gamma/w1)-0.5*bb*(c8-d8)

    allocate(f1(0:nx1,nequ,ngua))
    allocate(f2(0:nx1,nequ,ngua))
    allocate(g1(0:ny1,nequ,ngua))
    allocate(g2(0:ny1,nequ,ngua))
    allocate(resp1(nphi,nequ,0:nn1,0:nn1))
    allocate(resp2(nphi,nequ,0:nn1,0:nn1))
    allocate(resq1(nphi,nequ,0:nn1,0:nn1))
    allocate(resq2(nphi,nequ,0:nn1,0:nn1))

    !---------------init res
    do i=1,nx
        do j=1,ny
            do k=1,nequ
                do m=1,nphi
                    res(m,k,i,j)=0.
                    resp1(m,k,i,j)=0.
                    resp2(m,k,i,j)=0.
                    resq1(m,k,i,j)=0.
                    resq2(m,k,i,j)=0.
                end do
            end do
        end do
    end do

    !---------------------numerical flux

    !!$OMP PARALLEL DO
    do j=1,ny
        do i=0,nx
            do ixg=1,ngua
                rhop(ixg)=0.
                rmxp(ixg)=0.
                rmyp(ixg)=0.
                rmzp(ixg)=0.
                bxp(ixg)=0.
                byp(ixg)=0.
                bzp(ixg)=0.
                enep(ixg)=0.
                rhon(ixg)=0.
                rmxn(ixg)=0.
                rmyn(ixg)=0.
                rmzn(ixg)=0.
                bxn(ixg)=0.
                byn(ixg)=0.
                bzn(ixg)=0.
                enen(ixg)=0.
                do m=1,nphi
                    rhon(ixg)=rhon(ixg)+u(m,i,j,1,io)*phixb(m,2,ixg)
                    rmxn(ixg)=rmxn(ixg)+u(m,i,j,2,io)*phixb(m,2,ixg)
                    rmyn(ixg)=rmyn(ixg)+u(m,i,j,3,io)*phixb(m,2,ixg)
                    rmzn(ixg)=rmzn(ixg)+u(m,i,j,4,io)*phixb(m,2,ixg)
                    bxn(ixg)=bxn(ixg)+u(m,i,j,5,io)*phixb(m,2,ixg)
                    byn(ixg)=byn(ixg)+u(m,i,j,6,io)*phixb(m,2,ixg)
                    bzn(ixg)=bzn(ixg)+u(m,i,j,7,io)*phixb(m,2,ixg)
                    enen(ixg)=enen(ixg)+u(m,i,j,8,io)*phixb(m,2,ixg)
                    rhop(ixg)=rhop(ixg)+u(m,i+1,j,1,io)*phixb(m,1,ixg)
                    rmxp(ixg)=rmxp(ixg)+u(m,i+1,j,2,io)*phixb(m,1,ixg)
                    rmyp(ixg)=rmyp(ixg)+u(m,i+1,j,3,io)*phixb(m,1,ixg)
                    rmzp(ixg)=rmzp(ixg)+u(m,i+1,j,4,io)*phixb(m,1,ixg)
                    bxp(ixg)=bxp(ixg)+u(m,i+1,j,5,io)*phixb(m,1,ixg)
                    byp(ixg)=byp(ixg)+u(m,i+1,j,6,io)*phixb(m,1,ixg)
                    bzp(ixg)=bzp(ixg)+u(m,i+1,j,7,io)*phixb(m,1,ixg)
                    enep(ixg)=enep(ixg)+u(m,i+1,j,8,io)*phixb(m,1,ixg)
                end do
                f1(i,1,ixg)=rhop(ixg)
                f2(i,1,ixg)=rhon(ixg)
                f1(i,2,ixg)=rmxp(ixg)
                f2(i,2,ixg)=rmxn(ixg)
                f1(i,3,ixg)=rmyp(ixg)
                f2(i,3,ixg)=rmyn(ixg)
                f1(i,4,ixg)=rmzp(ixg)
                f2(i,4,ixg)=rmzn(ixg)
                f1(i,5,ixg)=bxp(ixg)
                f2(i,5,ixg)=bxn(ixg)
                f1(i,6,ixg)=byp(ixg)
                f2(i,6,ixg)=byn(ixg)
                f1(i,7,ixg)=bzp(ixg)
                f2(i,7,ixg)=bzn(ixg)
                f1(i,8,ixg)=enep(ixg)
                f2(i,8,ixg)=enen(ixg)
            end do
        end do
        do i=1,nx
            do k=1,nequ
                do m=1,nphi
                    do ixg=1,ngua
                        resp1(m,k,i,j)=resp1(m,k,i,j)+alp(ixg)*(f1(i,k,ixg)*phixb(m,2,ixg)-f1(i-1,k,ixg)*phixb(m,1,ixg))/dx
                        resp2(m,k,i,j)=resp2(m,k,i,j)+alp(ixg)*(f2(i,k,ixg)*phixb(m,2,ixg)-f2(i-1,k,ixg)*phixb(m,1,ixg))/dx
                    end do
                end do
            end do
        end do
    end do

    deallocate(f1)
    deallocate(f2)

    !!$OMP PARALLEL DO
    do i=1,nx
        do j=0,ny
            do ixg=1,ngua
                rhop(ixg)=0.
                rmxp(ixg)=0.
                rmyp(ixg)=0.
                rmzp(ixg)=0.
                bxp(ixg)=0.
                byp(ixg)=0.
                bzp(ixg)=0.
                enep(ixg)=0.
                rhon(ixg)=0.
                rmxn(ixg)=0.
                rmyn(ixg)=0.
                rmzn(ixg)=0.
                bxn(ixg)=0.
                byn(ixg)=0.
                bzn(ixg)=0.
                enen(ixg)=0.
                do m=1,nphi
                    rhon(ixg)=rhon(ixg)+u(m,i,j,1,io)*phiyb(m,ixg,2)
                    rmxn(ixg)=rmxn(ixg)+u(m,i,j,2,io)*phiyb(m,ixg,2)
                    rmyn(ixg)=rmyn(ixg)+u(m,i,j,3,io)*phiyb(m,ixg,2)
                    rmzn(ixg)=rmzn(ixg)+u(m,i,j,4,io)*phiyb(m,ixg,2)
                    bxn(ixg)=bxn(ixg)+u(m,i,j,5,io)*phiyb(m,ixg,2)
                    byn(ixg)=byn(ixg)+u(m,i,j,6,io)*phiyb(m,ixg,2)
                    bzn(ixg)=bzn(ixg)+u(m,i,j,7,io)*phiyb(m,ixg,2)
                    enen(ixg)=enen(ixg)+u(m,i,j,8,io)*phiyb(m,ixg,2)
                    rhop(ixg)=rhop(ixg)+u(m,i,j+1,1,io)*phiyb(m,ixg,1)
                    rmxp(ixg)=rmxp(ixg)+u(m,i,j+1,2,io)*phiyb(m,ixg,1)
                    rmyp(ixg)=rmyp(ixg)+u(m,i,j+1,3,io)*phiyb(m,ixg,1)
                    rmzp(ixg)=rmzp(ixg)+u(m,i,j+1,4,io)*phiyb(m,ixg,1)
                    bxp(ixg)=bxp(ixg)+u(m,i,j+1,5,io)*phiyb(m,ixg,1)
                    byp(ixg)=byp(ixg)+u(m,i,j+1,6,io)*phiyb(m,ixg,1)
                    bzp(ixg)=bzp(ixg)+u(m,i,j+1,7,io)*phiyb(m,ixg,1)
                    enep(ixg)=enep(ixg)+u(m,i,j+1,8,io)*phiyb(m,ixg,1)
                end do
                g1(j,1,ixg)=rhop(ixg)
                g2(j,1,ixg)=rhon(ixg)
                g1(j,2,ixg)=rmxp(ixg)
                g2(j,2,ixg)=rmxn(ixg)
                g1(j,3,ixg)=rmyp(ixg)
                g2(j,3,ixg)=rmyn(ixg)
                g1(j,4,ixg)=rmzp(ixg)
                g2(j,4,ixg)=rmzn(ixg)
                g1(j,5,ixg)=bxp(ixg)
                g2(j,5,ixg)=bxn(ixg)
                g1(j,6,ixg)=byp(ixg)
                g2(j,6,ixg)=byn(ixg)
                g1(j,7,ixg)=bzp(ixg)
                g2(j,7,ixg)=bzn(ixg)
                g1(j,8,ixg)=enep(ixg)
                g2(j,8,ixg)=enen(ixg)
            end do
        end do
        do j=1,ny
            do k=1,nequ
                do m=1,nphi
                    do ixg=1,ngua
                        resq1(m,k,i,j)=resq1(m,k,i,j)+alp(ixg)*(g1(j,k,ixg)*phiyb(m,ixg,2)-g1(j-1,k,ixg)*phiyb(m,ixg,1))/dy
                        resq2(m,k,i,j)=resq2(m,k,i,j)+alp(ixg)*(g2(j,k,ixg)*phiyb(m,ixg,2)-g2(j-1,k,ixg)*phiyb(m,ixg,1))/dy
                    end do
                end do
            end do
        end do
    end do

    deallocate(g1)
    deallocate(g2)

    !------------------------------surface integral

    !!$OMP PARALLEL DO
    do i=1,nx
        do j=1,ny
            do ixg=1,ngua
                do iyg=1,ngua
                    rhog=0.
                    rmxg=0.
                    rmyg=0.
                    rmzg=0.
                    bxg=0.
                    byg=0.
                    bzg=0.
                    eneg=0.
                    do m=1,nphi
                        rhog=rhog+u(m,i,j,1,io)*phig(m,ixg,iyg)
                        rmxg=rmxg+u(m,i,j,2,io)*phig(m,ixg,iyg)
                        rmyg=rmyg+u(m,i,j,3,io)*phig(m,ixg,iyg)
                        rmzg=rmzg+u(m,i,j,4,io)*phig(m,ixg,iyg)
                        bxg=bxg+u(m,i,j,5,io)*phig(m,ixg,iyg)
                        byg=byg+u(m,i,j,6,io)*phig(m,ixg,iyg)
                        bzg=bzg+u(m,i,j,7,io)*phig(m,ixg,iyg)
                        eneg=eneg+u(m,i,j,8,io)*phig(m,ixg,iyg)
                    end do
                    do m=1,nphi
                        resp1(m,1,i,j)=resp1(m,1,i,j)-alp(ixg)*alp(iyg)*rhog*d1phig(m,ixg,iyg)/dx
                        resp2(m,1,i,j)=resp2(m,1,i,j)-alp(ixg)*alp(iyg)*rhog*d1phig(m,ixg,iyg)/dx
                        resp1(m,2,i,j)=resp1(m,2,i,j)-alp(ixg)*alp(iyg)*rmxg*d1phig(m,ixg,iyg)/dx
                        resp2(m,2,i,j)=resp2(m,2,i,j)-alp(ixg)*alp(iyg)*rmxg*d1phig(m,ixg,iyg)/dx
                        resp1(m,3,i,j)=resp1(m,3,i,j)-alp(ixg)*alp(iyg)*rmyg*d1phig(m,ixg,iyg)/dx
                        resp2(m,3,i,j)=resp2(m,3,i,j)-alp(ixg)*alp(iyg)*rmyg*d1phig(m,ixg,iyg)/dx
                        resp1(m,4,i,j)=resp1(m,4,i,j)-alp(ixg)*alp(iyg)*rmzg*d1phig(m,ixg,iyg)/dx
                        resp2(m,4,i,j)=resp2(m,4,i,j)-alp(ixg)*alp(iyg)*rmzg*d1phig(m,ixg,iyg)/dx
                        resp1(m,5,i,j)=resp1(m,5,i,j)-alp(ixg)*alp(iyg)*bxg*d1phig(m,ixg,iyg)/dx
                        resp2(m,5,i,j)=resp2(m,5,i,j)-alp(ixg)*alp(iyg)*bxg*d1phig(m,ixg,iyg)/dx
                        resp1(m,6,i,j)=resp1(m,6,i,j)-alp(ixg)*alp(iyg)*byg*d1phig(m,ixg,iyg)/dx
                        resp2(m,6,i,j)=resp2(m,6,i,j)-alp(ixg)*alp(iyg)*byg*d1phig(m,ixg,iyg)/dx
                        resp1(m,7,i,j)=resp1(m,7,i,j)-alp(ixg)*alp(iyg)*bzg*d1phig(m,ixg,iyg)/dx
                        resp2(m,7,i,j)=resp2(m,7,i,j)-alp(ixg)*alp(iyg)*bzg*d1phig(m,ixg,iyg)/dx
                        resp1(m,8,i,j)=resp1(m,8,i,j)-alp(ixg)*alp(iyg)*eneg*d1phig(m,ixg,iyg)/dx
                        resp2(m,8,i,j)=resp2(m,8,i,j)-alp(ixg)*alp(iyg)*eneg*d1phig(m,ixg,iyg)/dx
                        resq1(m,1,i,j)=resq1(m,1,i,j)-alp(ixg)*alp(iyg)*rhog*d2phig(m,ixg,iyg)/dy
                        resq2(m,1,i,j)=resq2(m,1,i,j)-alp(ixg)*alp(iyg)*rhog*d2phig(m,ixg,iyg)/dy
                        resq1(m,2,i,j)=resq1(m,2,i,j)-alp(ixg)*alp(iyg)*rmxg*d2phig(m,ixg,iyg)/dy
                        resq2(m,2,i,j)=resq2(m,2,i,j)-alp(ixg)*alp(iyg)*rmxg*d2phig(m,ixg,iyg)/dy
                        resq1(m,3,i,j)=resq1(m,3,i,j)-alp(ixg)*alp(iyg)*rmyg*d2phig(m,ixg,iyg)/dy
                        resq2(m,3,i,j)=resq2(m,3,i,j)-alp(ixg)*alp(iyg)*rmyg*d2phig(m,ixg,iyg)/dy
                        resq1(m,4,i,j)=resq1(m,4,i,j)-alp(ixg)*alp(iyg)*rmzg*d2phig(m,ixg,iyg)/dy
                        resq2(m,4,i,j)=resq2(m,4,i,j)-alp(ixg)*alp(iyg)*rmzg*d2phig(m,ixg,iyg)/dy
                        resq1(m,5,i,j)=resq1(m,5,i,j)-alp(ixg)*alp(iyg)*bxg*d2phig(m,ixg,iyg)/dy
                        resq2(m,5,i,j)=resq2(m,5,i,j)-alp(ixg)*alp(iyg)*bxg*d2phig(m,ixg,iyg)/dy
                        resq1(m,6,i,j)=resq1(m,6,i,j)-alp(ixg)*alp(iyg)*byg*d2phig(m,ixg,iyg)/dy
                        resq2(m,6,i,j)=resq2(m,6,i,j)-alp(ixg)*alp(iyg)*byg*d2phig(m,ixg,iyg)/dy
                        resq1(m,7,i,j)=resq1(m,7,i,j)-alp(ixg)*alp(iyg)*bzg*d2phig(m,ixg,iyg)/dy
                        resq2(m,7,i,j)=resq2(m,7,i,j)-alp(ixg)*alp(iyg)*bzg*d2phig(m,ixg,iyg)/dy
                        resq1(m,8,i,j)=resq1(m,8,i,j)-alp(ixg)*alp(iyg)*eneg*d2phig(m,ixg,iyg)/dy
                        resq2(m,8,i,j)=resq2(m,8,i,j)-alp(ixg)*alp(iyg)*eneg*d2phig(m,ixg,iyg)/dy
                    end do
                end do
            end do
        end do
    end do

    !------------------------source term p1,p2
    do i=1,nx
        do j=1,nphi
            do k=1,nequ
                do m=1,nphi
                    do n=1,nphi
                        p1(m,i,j,k,io)=resp1(m,k,i,j)
                        p2(m,i,j,k,io)=resp2(m,k,i,j)
                        q1(m,i,j,k,io)=resq1(m,k,i,j)
                        q2(m,i,j,k,io)=resq2(m,k,i,j)
                    end do
                end do
            end do
        end do
    end do

    alphamax=0.
    betamax=0.
    !    !!$OMP PARALLEL DO
    do i=1,nx
        do j=1,ny
            rho=u(1,i,j,1,io)
            velx=abs(u(1,i,j,2,io)/rho)
            vely=abs(u(1,i,j,3,io)/rho)
            pres=fxpress(rho,u(1,i,j,2,io),u(1,i,j,3,io),u(1,i,j,4,io),u(1,i,j,5,io),u(1,i,j,6,io),u(1,i,j,7,io),u(1,i,j,8,io))
            ca2=gamma*pres/rho
            cb2=(u(1,i,j,5,io)*u(1,i,j,5,io)+u(1,i,j,6,io)*u(1,i,j,6,io)+u(1,i,j,7,io)*u(1,i,j,7,io))/rho
            cx1=ca2+cb2
            cx2=sqrt(abs(cx1*cx1-4.*ca2*u(1,i,j,5,io)*u(1,i,j,5,io)/rho))
            cy2=sqrt(abs(cx1*cx1-4.*ca2*u(1,i,j,6,io)*u(1,i,j,6,io)/rho))
            cfx=sqrt(abs(0.5*(cx1+cx2)))
            cfy=sqrt(abs(0.5*(cx1+cy2)))
            alphamax=max(alphamax,cfx+velx)
            betamax=max(betamax,cfy+vely)
        end do
    end do

    !!$OMP PARALLEL DO
    do i=1,nx
        do j=1,ny
            do ixg=1,ngua
                do iyg=1,ngua
                    rhog=0.
                    rmxg=0.
                    rmyg=0.
                    rmzg=0.
                    bxg=0.
                    byg=0.
                    bzg=0.
                    eneg=0.
                    p1rhog=0.
                    p1rmxg=0.
                    p1rmyg=0.
                    p1rmzg=0.
                    p1bxg=0.
                    p1byg=0.
                    p1bzg=0.
                    p1eneg=0.
                    p2rhog=0.
                    p2rmxg=0.
                    p2rmyg=0.
                    p2rmzg=0.
                    p2bxg=0.
                    p2byg=0.
                    p2bzg=0.
                    p2eneg=0.
                    q1rhog=0.
                    q1rmxg=0.
                    q1rmyg=0.
                    q1rmzg=0.
                    q1bxg=0.
                    q1byg=0.
                    q1bzg=0.
                    q1eneg=0.
                    q2rhog=0.
                    q2rmxg=0.
                    q2rmyg=0.
                    q2rmzg=0.
                    q2bxg=0.
                    q2byg=0.
                    q2bzg=0.
                    q2eneg=0.
                    do m=1,nphi
                        rhog=rhog+u(m,i,j,1,io)*phig(m,ixg,iyg)
                        rmxg=rmxg+u(m,i,j,2,io)*phig(m,ixg,iyg)
                        rmyg=rmyg+u(m,i,j,3,io)*phig(m,ixg,iyg)
                        rmzg=rmzg+u(m,i,j,4,io)*phig(m,ixg,iyg)
                        bxg=bxg+u(m,i,j,5,io)*phig(m,ixg,iyg)
                        byg=byg+u(m,i,j,6,io)*phig(m,ixg,iyg)
                        bzg=bzg+u(m,i,j,7,io)*phig(m,ixg,iyg)
                        eneg=eneg+u(m,i,j,8,io)*phig(m,ixg,iyg)
                        p1rhog=p1rhog+resp1(m,1,i,j)*phig(m,ixg,iyg)
                        p1rmxg=p1rmxg+resp1(m,2,i,j)*phig(m,ixg,iyg)
                        p1rmyg=p1rmyg+resp1(m,3,i,j)*phig(m,ixg,iyg)
                        p1rmzg=p1rmzg+resp1(m,4,i,j)*phig(m,ixg,iyg)
                        p1bxg=p1bxg+resp1(m,5,i,j)*phig(m,ixg,iyg)
                        p1byg=p1byg+resp1(m,6,i,j)*phig(m,ixg,iyg)
                        p1bzg=p1bzg+resp1(m,7,i,j)*phig(m,ixg,iyg)
                        p1eneg=p1eneg+resp1(m,8,i,j)*phig(m,ixg,iyg)
                        p2rhog=p2rhog+resp2(m,1,i,j)*phig(m,ixg,iyg)
                        p2rmxg=p2rmxg+resp2(m,2,i,j)*phig(m,ixg,iyg)
                        p2rmyg=p2rmyg+resp2(m,3,i,j)*phig(m,ixg,iyg)
                        p2rmzg=p2rmzg+resp2(m,4,i,j)*phig(m,ixg,iyg)
                        p2bxg=p2bxg+resp2(m,5,i,j)*phig(m,ixg,iyg)
                        p2byg=p2byg+resp2(m,6,i,j)*phig(m,ixg,iyg)
                        p2bzg=p2bzg+resp2(m,7,i,j)*phig(m,ixg,iyg)
                        p2eneg=p2eneg+resp2(m,8,i,j)*phig(m,ixg,iyg)
                        q1rhog=q1rhog+resq1(m,1,i,j)*phig(m,ixg,iyg)
                        q1rmxg=q1rmxg+resq1(m,2,i,j)*phig(m,ixg,iyg)
                        q1rmyg=q1rmyg+resq1(m,3,i,j)*phig(m,ixg,iyg)
                        q1rmzg=q1rmzg+resq1(m,4,i,j)*phig(m,ixg,iyg)
                        q1bxg=q1bxg+resq1(m,5,i,j)*phig(m,ixg,iyg)
                        q1byg=q1byg+resq1(m,6,i,j)*phig(m,ixg,iyg)
                        q1bzg=q1bzg+resq1(m,7,i,j)*phig(m,ixg,iyg)
                        q1eneg=q1eneg+resq1(m,8,i,j)*phig(m,ixg,iyg)
                        q2rhog=q2rhog+resq2(m,1,i,j)*phig(m,ixg,iyg)
                        q2rmxg=q2rmxg+resq2(m,2,i,j)*phig(m,ixg,iyg)
                        q2rmyg=q2rmyg+resq2(m,3,i,j)*phig(m,ixg,iyg)
                        q2rmzg=q2rmzg+resq2(m,4,i,j)*phig(m,ixg,iyg)
                        q2bxg=q2bxg+resq2(m,5,i,j)*phig(m,ixg,iyg)
                        q2byg=q2byg+resq2(m,6,i,j)*phig(m,ixg,iyg)
                        q2bzg=q2bzg+resq2(m,7,i,j)*phig(m,ixg,iyg)
                        q2eneg=q2eneg+resq2(m,8,i,j)*phig(m,ixg,iyg)
                    end do
                    pressg=fxpress(rhog,rmxg,rmyg,rmzg,bxg,byg,bzg,eneg)

                    frhog=fxrho(p1rhog,p1rmxg,p2rhog,p2rmxg,alphamax)
                    frmxg=fxrmx(rhog,rmxg,rmyg,rmzg,bxg,byg,bzg,p1rhog,p1rmxg,p1rmyg,p1rmzg,p1bxg,p1byg,p1bzg,p1eneg,p2rhog,p2rmxg,p2rmyg,p2rmzg,p2bxg,p2byg,p2bzg,p2eneg,alphamax,gamma)
                    frmyg=fxrmy(rhog,rmxg,rmyg,bxg,byg,p1rhog,p1rmxg,p1rmyg,p1bxg,p1byg,p2rhog,p2rmxg,p2rmyg,p2bxg,p2byg,alphamax)
                    frmzg=fxrmz(rhog,rmxg,rmzg,bxg,bzg,p1rhog,p1rmxg,p1rmzg,p1bxg,p1bzg,p2rhog,p2rmxg,p2rmzg,p2bxg,p2bzg,alphamax)
                    fbxg=fxbx(p1bxg,p2bxg,alphamax)
                    fbyg=fxby(rhog,rmxg,rmyg,bxg,byg,p1rhog,p1rmxg,p1rmyg,p1bxg,p1byg,p2rhog,p2rmxg,p2rmyg,p2bxg,p2byg,alphamax)
                    fbzg=fxbz(rhog,rmxg,rmzg,bxg,bzg,p1rhog,p1rmxg,p1rmzg,p1bxg,p1bzg,p2rhog,p2rmxg,p2rmzg,p2bxg,p2bzg,alphamax)
                    feneg=fxene(rhog,rmxg,rmyg,rmzg,bxg,byg,bzg,p1rhog,p1rmxg,p1rmyg,p1rmzg,p1bxg,p1byg,p1bzg,p1eneg,p2rhog,p2rmxg,p2rmyg,p2rmzg,p2bxg,p2byg,p2bzg,p2eneg,alphamax,gamma,pressg)

                    grhog=gyrho(q1rhog,q1rmyg,q2rhog,q2rmyg,betamax)
                    grmxg=gyrmx(rhog,rmxg,rmyg,bxg,byg,q1rhog,q1rmxg,q1rmyg,q1bxg,q1byg,q2rhog,q2rmxg,q2rmyg,q2bxg,q2byg,betamax)
                    grmyg=gyrmy(rhog,rmxg,rmyg,rmzg,bxg,byg,bzg,q1rhog,q1rmxg,q1rmyg,q1rmzg,q1bxg,q1byg,q1bzg,q1eneg,q2rhog,q2rmxg,q2rmyg,q2rmzg,q2bxg,q2byg,q2bzg,q2eneg,betamax,gamma)
                    grmzg=gyrmz(rhog,rmyg,rmzg,byg,bzg,q1rhog,q1rmyg,q1rmzg,q1byg,q1bzg,q2rhog,q2rmyg,q2rmzg,q2byg,q2bzg,betamax)
                    gbxg=gybx(rhog,rmxg,rmyg,bxg,byg,q1rhog,q1rmxg,q1rmyg,q1bxg,q1byg,q2rhog,q2rmxg,q2rmyg,q2bxg,q2byg,betamax)
                    gbyg=gyby(q1byg,q2byg,betamax)
                    gbzg=gybz(rhog,rmyg,rmzg,byg,bzg,q1rhog,q1rmyg,q1rmzg,q1byg,q1bzg,q2rhog,q2rmyg,q2rmzg,q2byg,q2bzg,betamax)
                    geneg=gyene(rhog,rmxg,rmyg,rmzg,bxg,byg,bzg,q1rhog,q1rmxg,q1rmyg,q1rmzg,q1bxg,q1byg,q1bzg,q1eneg,q2rhog,q2rmxg,q2rmyg,q2rmzg,q2bxg,q2byg,q2bzg,q2eneg,betamax,gamma,pressg)

                    zx=x(i)+dx*zg(ixg)
                    zy=y(j)+dy*zg(iyg)

                    !s1=0.
                    !s2=-0.5*bxg*(p1bxg+p2bxg+q1byg+q2byg)
                    !s3=-0.5*byg*(p1bxg+p2bxg+q1byg+q2byg)
                    !s4=-0.5*bzg*(p1bxg+p2bxg+q1byg+q2byg)
                    !s5=-0.5*(rmxg/rhog)*(p1bxg+p2bxg+q1byg+q2byg)
                    !s6=-0.5*(rmyg/rhog)*(p1bxg+p2bxg+q1byg+q2byg)
                    !s7=-0.5*(rmzg/rhog)*(p1bxg+p2bxg+q1byg+q2byg)
                    !s8=-0.5*(rmxg*bxg/rhog+rmyg*byg/rhog+rmzg*bzg/rhog)*(p1bxg+p2bxg+q1byg+q2byg)

                    s1=0.
                    s2=0.
                    s3=0.
                    s4=0.
                    s5=0.
                    s6=0.
                    s7=0.
                    s8=0.

                    do m=1,nphi
                        res(m,1,i,j)=res(m,1,i,j)-alp(ixg)*alp(iyg)*(frhog+grhog)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s1*phig(m,ixg,iyg)
                        res(m,2,i,j)=res(m,2,i,j)-alp(ixg)*alp(iyg)*(frmxg+grmxg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s2*phig(m,ixg,iyg)
                        res(m,3,i,j)=res(m,3,i,j)-alp(ixg)*alp(iyg)*(frmyg+grmyg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s3*phig(m,ixg,iyg)
                        res(m,4,i,j)=res(m,4,i,j)-alp(ixg)*alp(iyg)*(frmzg+grmzg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s4*phig(m,ixg,iyg)
                        res(m,5,i,j)=res(m,5,i,j)-alp(ixg)*alp(iyg)*(fbxg+gbxg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s5*phig(m,ixg,iyg)
                        res(m,6,i,j)=res(m,6,i,j)-alp(ixg)*alp(iyg)*(fbyg+gbyg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s6*phig(m,ixg,iyg)
                        res(m,7,i,j)=res(m,7,i,j)-alp(ixg)*alp(iyg)*(fbzg+gbzg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s7*phig(m,ixg,iyg)
                        res(m,8,i,j)=res(m,8,i,j)-alp(ixg)*alp(iyg)*(feneg+geneg)*phig(m,ixg,iyg)+alp(ixg)*alp(iyg)*s8*phig(m,ixg,iyg)
                    end do
                end do
            end do
        end do
    end do

    deallocate(resp1)
    deallocate(resp2)
    deallocate(resq1)
    deallocate(resq2)

    return

    end subroutine fx


    subroutine rk(io)
    use shared_data
    !include 'com.txt'

    !-------------------------------------------THE RKP2 DISCONTINUOUS GALERKIN METHOD

    if(io.eq.0) then
        !!$OMP PARALLEL DO
        do 13 i=1,nx
            do 13 j=1,ny
                do 13 k=1,nequ
                    do 13 m=1,nphi
                        u(m,i,j,k,1)=u(m,i,j,k,0)+dt*res(m,k,i,j)
13      continue


    else if(io.eq.1) then
        !!$OMP PARALLEL DO
        do 14 i=1,nx
            do 14 j=1,ny
                do 14 k=1,nequ
                    do 14 m=1,nphi
                        u(m,i,j,k,2)=0.75*u(m,i,j,k,0)+0.25*(u(m,i,j,k,io)+dt*res(m,k,i,j))
14      continue

    else
        !!$OMP PARALLEL DO
        do 15 i=1,nx
            do 15 j=1,ny
                do 15 k=1,nequ
                    do 15 m=1,nphi
                        u(m,i,j,k,0)=(u(m,i,j,k,0)+2.*(u(m,i,j,k,io)+dt*res(m,k,i,j)))/3.
15      continue
    endif


    return


    end subroutine rk


    subroutine outp
    use shared_data
    !include 'com.txt'

    character(len=512)::cFile
    real::globaldiv(nphi)

    errh1=0.0
    errh2=0.0
    errhinf=0.0
    rhoerrh1=0.0
    rhoerrh2=0.0
    rhoerrhinf=0.0
    uxerrh1=0.0
    uxerrh2=0.0
    uxerrhinf=0.0
    uyerrh1=0.0
    uyerrh2=0.0
    uyerrhinf=0.0
    uzerrh1=0.0
    uzerrh2=0.0
    uzerrhinf=0.0
    bxerrh1=0.0
    bxerrh2=0.0
    bxerrhinf=0.0
    byerrh1=0.0
    byerrh2=0.0
    byerrhinf=0.0
    bzerrh1=0.0
    bzerrh2=0.0
    bzerrhinf=0.0
    perrh1=0.0
    perrh2=0.0
    perrhinf=0.0
    diverrh1=0.0
    diverrh2=0.0
    diverrhinf=0.0

    write(cFile, 101) nx
101 FORMAT('rho',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=11, file=cFile)
    write(11, *) 'VARIABLES="X","Y","<greek>r<\greek>"'
    write(11, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 102) nx
102 FORMAT('ux',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=22, file=cFile)
    write(22, *) 'VARIABLES="X","Y","ux"'
    write(22, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 103) nx
103 FORMAT('uy',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=33, file=cFile)
    write(33, *) 'VARIABLES="X","Y","uy"'
    write(33, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 104) nx
104 FORMAT('uz',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=44, file=cFile)
    write(44, *) 'VARIABLES="X","Y","uz"'
    write(44, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 105) nx
105 FORMAT('bx',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=55, file=cFile)
    write(55, *) 'VARIABLES="X","Y","bx"'
    write(55, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 106) nx
106 FORMAT('by',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=66, file=cFile)
    write(66, *) 'VARIABLES="X","Y","by"'
    write(66, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 107) nx
107 FORMAT('bz',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=77, file=cFile)
    write(77, *) 'VARIABLES="X","Y","bz"'
    write(77, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 108) nx
108 FORMAT('pressure',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=88, file=cFile)
    write(88, *) 'VARIABLES="X","Y","p"'
    write(88, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 109) nx
109 FORMAT('energy',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=99, file=cFile)
    write(99, *) 'VARIABLES="X","Y","E"'
    write(99, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 110) nx
110 FORMAT('magnetic pressure',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=00, file=cFile)
    write(00, *) 'VARIABLES="X","Y","Magnetic pressure"'
    write(00, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 111) nx
111 FORMAT('divergence',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=100, file=cFile)
    write(100, *) 'VARIABLES="X","Y","Divergence"'
    write(100, *) 'ZONE F=POINT, I=', nx, 'J=', ny

    write(cFile, 112) nx
112 FORMAT('global divergence',I3)

    cFile=Trim(AdjustL(cFile))//'.dat'

    open(unit=200, file=cFile)
    write(200, *) 'VARIABLES="T","Global Divergence"'
    write(200, *) 'I=75, J=75'

    !!$OMP PARALLEL DO

    do m=1,nphi
        globaldiv(m)=0.
    end do

    do i=1,nx
        do j=1,ny
            zx=x(i)
            zy=y(j)

            rho=0.
            rmx=0.
            rmy=0.
            rmz=0.
            bx=0.
            divx1=0.
            divx2=0.
            divy1=0.
            divy2=0.
            by=0.
            bz=0.
            ene=0.
            do m=1,nphi
                rho=rho+u(m,i,j,1,0)*phip(m)
                rmx=rmx+u(m,i,j,2,0)*phip(m)
                rmy=rmy+u(m,i,j,3,0)*phip(m)
                rmz=rmz+u(m,i,j,4,0)*phip(m)
                bx=bx+u(m,i,j,5,0)*phip(m)
                divx1=divx1+p1(m,i,j,5,0)*phip(m)
                divx2=divx2+p2(m,i,j,5,0)*phip(m)
                divy1=divy1+q1(m,i,j,6,0)*phip(m)
                divy2=divy2+q2(m,i,j,6,0)*phip(m)
                by=by+u(m,i,j,6,0)*phip(m)
                bz=bz+u(m,i,j,7,0)*phip(m)
                ene=ene+u(m,i,j,8,0)*phip(m)
            end do

            !do ixg=1,ngua
            !    divxn1=0.
            !    divxn2=0.
            !    divxp1=0.
            !    divxp2=0.
            !    divyn1=0.
            !    divyn2=0.
            !    divyp1=0.
            !    divyp2=0.
            !    do m=1,nphi
            !        if(i-1>=1) then
            !            divxn1=divxn1+0.5*(p1(m,i-1,j,5,0)+p2(m,i-1,j,5,0)+q1(m,i-1,j,6,0)+q2(m,i-1,j,6,0))*phixb(m,2,ixg)
            !        end if
            !        divxn2=divxn2+0.5*(p1(m,i,j,5,0)+p2(m,i,j,5,0)+q1(m,i,j,6,0)+q2(m,i,j,6,0))*phixb(m,2,ixg)
            !        divxp1=divxp1+0.5*(p1(m,i,j,5,0)+p2(m,i,j,5,0)+q1(m,i,j,6,0)+q2(m,i,j,6,0))*phixb(m,1,ixg)
            !        if(i+1<=nx) then
            !            divxp2=divxp2+0.5*(p1(m,i+1,j,5,0)+p2(m,i+1,j,5,0)+q1(m,i+1,j,6,0)+q2(m,i+1,j,6,0))*phixb(m,1,ixg)
            !        end if
            !        if(j-1>=1) then
            !            divyn1=divyn1+0.5*(p1(m,i,j-1,5,0)+p2(m,i,j-1,5,0)+q1(m,i,j-1,6,0)+q2(m,i,j-1,6,0))*phiyb(m,ixg,2)
            !        end if
            !        divyn2=divyn2+0.5*(p1(m,i,j,5,0)+p2(m,i,j,5,0)+q1(m,i,j,6,0)+q2(m,i,j,6,0))*phiyb(m,ixg,2)
            !        divyp1=divyp1+0.5*(p1(m,i,j,5,0)+p2(m,i,j,5,0)+q1(m,i,j,6,0)+q2(m,i,j,6,0))*phiyb(m,ixg,1)
            !        if(j+1<=ny) then
            !            divyp2=divyp2+0.5*(p1(m,i,j+1,5,0)+p2(m,i,j+1,5,0)+q1(m,i,j+1,6,0)+q2(m,i,j+1,6,0))*phiyb(m,ixg,1)
            !        end if
            !    end do
            !
            !    do m=1,nphi
            !        if(i-1>=1) then
            !            globaldiv(m)=globaldiv(m)+alp(ixg)*abs(divxp1-divxn1)*dy
            !        end if
            !        if(i+1<=nx) then
            !            globaldiv(m)=globaldiv(m)+alp(ixg)*abs(divxp2-divxn2)*dy
            !        end if
            !        if(j-1>=1) then
            !            globaldiv(m)=globaldiv(m)+alp(ixg)*abs(divyp1-divyn1)*dy
            !        end if
            !        if(j+1<=ny) then
            !            globaldiv(m)=globaldiv(m)+alp(ixg)*abs(divyp2-divyn2)*dy
            !        end if
            !    end do
            !end do
            
            do ixg=1,ngua
                bxn1=0.
                bxn2=0.
                bxp1=0.
                bxp2=0.
                byn1=0.
                byn2=0.
                byp1=0.
                byp2=0.
                do m=1,nphi
                    bxn1=bxn1+u(m,i-1,j,5,0)*phixb(m,2,ixg)
                    bxn2=bxn2+u(m,i,j,5,0)*phixb(m,2,ixg)
                    bxp1=bxp1+u(m,i,j,5,0)*phixb(m,1,ixg)
                    bxp2=bxp2+u(m,i+1,j,5,0)*phixb(m,1,ixg)
                    byn1=byn1+u(m,i,j-1,6,0)*phiyb(m,ixg,2)
                    byn2=byn2+u(m,i,j,6,0)*phiyb(m,ixg,2)
                    byp1=byp1+u(m,i,j,6,0)*phiyb(m,ixg,1)
                    byp2=byp2+u(m,i,j+1,6,0)*phiyb(m,ixg,1)
                end do

                do m=1,nphi
                    globaldiv(m)=globaldiv(m)+alp(ixg)*abs(bxp1-bxn1)*dy

                    globaldiv(m)=globaldiv(m)+alp(ixg)*abs(bxp2-bxn2)*dy

                    globaldiv(m)=globaldiv(m)+alp(ixg)*abs(byp1-byn1)*dx

                    globaldiv(m)=globaldiv(m)+alp(ixg)*abs(byp2-byn2)*dx
                end do
            end do
            do ixg=1,ngua
                do jyg=1,ngua
                    p1bxg=0.
                    p2bxg=0.
                    q1byg=0.
                    q2byg=0.
                    do m=1,nphi
                        p1bxg=p1bxg+p1(m,i,j,5,0)*phig(m,ixg,jyg)
                        p2bxg=p2bxg+p2(m,i,j,5,0)*phig(m,ixg,jyg)
                        q1byg=q1byg+q1(m,i,j,6,0)*phig(m,ixg,jyg)
                        q2byg=q2byg+q2(m,i,j,6,0)*phig(m,ixg,jyg)
                    end do
                    do m=1,nphi
                        globaldiv(m)=globaldiv(m)+alp(ixg)*alp(jyg)*abs(0.5*(p1bxg+p2bxg+q1byg+q2byg))*dy*dx
                    end do
                end do
            end do

            !write(11, *) zx,zy,rho
            !
            !
            !write(22, *) zx,zy,rmx/rho
            !
            !
            !write(33, *) zx,zy,rmy/rho
            !
            !
            !write(44, *) zx,zy,rmz/rho
            !
            !
            !write(55, *) zx,zy,bx
            !
            !
            !write(100, *) zx,zy,0.5*(divx1+divx2+divy1+divy2)
            !
            !
            !write(66, *) zx,zy,by
            !
            !
            !write(77, *) zx,zy,bz
            !
            press=gamma1*(ene-0.5*(rmx*rmx+rmy*rmy+rmz*rmz)/rho-0.5*(bx*bx+by*by+bz*bz))
            !
            !write(88, *) zx,zy,press
            !
            !write(99, *) zx,zy,ene
            !
            !write(00, *) zx,zy,0.5*(bx*bx+by*by+bz*bz)

            if(zx>=-5..and.zx<=5..and.zy>=-5..and.zy<=5.) then
                zx=zx-tn
                zy=zy-tn

                exrho=1.
                exux=1.-(zy/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
                exuy=1.+(zx/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
                exuz=0.
                exbx=-(zy/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
                exby=(zx/(2.*pi))*exp(0.5*(1.-zx*zx-zy*zy))
                exbz=0.
                express=1.+((-zx*zx-zy*zy)*exp(1.-zx*zx-zy*zy))/(8.*pi*pi)
                exdiv=0.

                rhoerrh=abs(exrho-rho)
                rhoerrh1=rhoerrh1+rhoerrh
                rhoerrh2=rhoerrh2+rhoerrh*rhoerrh
                rhoerrhinf=max(rhoerrh,rhoerrhinf)

                uxerrh=abs(exux-rmx/rho)
                uxerrh1=uxerrh1+uxerrh
                uxerrh2=uxerrh2+uxerrh*uxerrh
                uxerrhinf=max(uxerrh,uxerrhinf)

                uyerrh=abs(exuy-rmy/rho)
                uyerrh1=uyerrh1+uyerrh
                uyerrh2=uyerrh2+uyerrh*uyerrh
                uyerrhinf=max(uyerrh,uyerrhinf)

                uzerrh=abs(exuz-rmz/rho)
                uzerrh1=uzerrh1+uzerrh
                uzerrh2=uzerrh2+uzerrh*uzerrh
                uzerrhinf=max(uzerrh,uzerrhinf)

                bxerrh=abs(exbx-bx)
                bxerrh1=bxerrh1+bxerrh
                bxerrh2=bxerrh2+bxerrh*bxerrh
                bxerrhinf=max(bxerrh,bxerrhinf)

                byerrh=abs(exby-by)
                byerrh1=byerrh1+byerrh
                byerrh2=byerrh2+byerrh*byerrh
                byerrhinf=max(byerrh,byerrhinf)

                bzerrh=abs(exbz-bz)
                bzerrh1=bzerrh1+bzerrh
                bzerrh2=bzerrh2+bzerrh*bzerrh
                bzerrhinf=max(bzerrh,bzerrhinf)

                perrh=abs(express-press)
                perrh1=perrh1+perrh
                perrh2=perrh2+perrh*perrh
                perrhinf=max(perrh,perrhinf)

                diverrh=abs(div-0.5*(divx1+divx2+divy1+divy2))
                diverrh1=diverrh1+diverrh
                diverrh2=diverrh2+diverrh*diverrh
                diverrhinf=max(diverrh,diverrhinf)
            end if
        end do
    end do

    globaldivergence=0.
    do m=1,nphi
        globaldivergence=globaldivergence+globaldiv(m)*phip(m)
    end do

    globaldivergence=0.5*globaldivergence

    rhoerrh1=4.*rhoerrh1/(nx*ny)
    rhoerrh2=sqrt(4.*(rhoerrh2/nx)/ny)
    rhoorderrh1=log(rhoerrh1old/rhoerrh1)/log(2.)
    rhoorderrh2=log(rhoerrh2old/rhoerrh2)/log(2.)
    rhoorderrhinf=log(rhoerrhinfold/rhoerrhinf)/log(2.)
    rhoerrh1old=rhoerrh1
    rhoerrh2old=rhoerrh2
    rhoerrhinfold=rhoerrhinf

    write(11,*) 'tn=',tn
    write(11,*) 'nx=',nx,'ny=',ny
    write(11,*) 'L1-rh=',rhoerrh1,'  order=',rhoorderrh1
    write(11,*) 'Linf-rh=',rhoerrhinf,'  order=',rhoorderrhinf
    write(11,*) 'L2-rh=',rhoerrh2,'  order=',rhoorderrh2

    uxerrh1=4.*(uxerrh1/nx)/ny
    uxerrh2=sqrt(4.*(uxerrh2/nx)/ny)
    uxorderrh1=log(uxerrh1old/uxerrh1)/log(2.)
    uxorderrh2=log(uxerrh2old/uxerrh2)/log(2.)
    uxorderrhinf=log(uxerrhinfold/uxerrhinf)/log(2.)
    uxerrh1old=uxerrh1
    uxerrh2old=uxerrh2
    uxerrhinfold=uxerrhinf

    write(22,*) 'tn=',tn
    write(22,*) 'nx=',nx,'ny=',ny
    write(22,*) 'L1-rh=',uxerrh1,'  order=',uxorderrh1
    write(22,*) 'Linf-rh=',uxerrhinf,'  order=',uxorderrhinf
    write(22,*) 'L2-rh=',uxerrh2,'  order=',uxorderrh2


    uyerrh1=4.*(uyerrh1/nx)/ny
    uyerrh2=sqrt(4.*(uyerrh2/nx)/ny)
    uyorderrh1=log(uyerrh1old/uyerrh1)/log(2.)
    uyorderrh2=log(uyerrh2old/uyerrh2)/log(2.)
    uyorderrhinf=log(uyerrhinfold/uyerrhinf)/log(2.)
    uyerrh1old=uyerrh1
    uyerrh2old=uyerrh2
    uyerrhinfold=uyerrhinf

    write(33,*) 'tn=',tn
    write(33,*) 'nx=',nx,'ny=',ny
    write(33,*) 'L1-rh=',uyerrh1,'  order=',uyorderrh1
    write(33,*) 'Linf-rh=',uyerrhinf,'  order=',uyorderrhinf
    write(33,*) 'L2-rh=',uyerrh2,'  order=',uyorderrh2


    uzerrh1=4.*(uzerrh1/nx)/ny
    uzerrh2=sqrt(4.*(uzerrh2/nx)/ny)
    uzorderrh1=log(uzerrh1old/uzerrh1)/log(2.)
    uzorderrh2=log(uzerrh2old/uzerrh2)/log(2.)
    uzorderrhinf=log(uzerrhinfold/uzerrhinf)/log(2.)
    uzerrh1old=uzerrh1
    uzerrh2old=uzerrh2
    uzerrhinfold=uzerrhinf

    write(44,*) 'tn=',tn
    write(44,*) 'nx=',nx,'ny=',ny
    write(44,*) 'L1-rh=',uzerrh1,'  order=',uzorderrh1
    write(44,*) 'Linf-rh=',uzerrhinf,'  order=',uzorderrhinf
    write(44,*) 'L2-rh=',uzerrh2,'  order=',uzorderrh2


    bxerrh1=4.*(bxerrh1/nx)/ny
    bxerrh2=sqrt(4.*(bxerrh2/nx)/ny)
    bxorderrh1=log(bxerrh1old/bxerrh1)/log(2.)
    bxorderrh2=log(bxerrh2old/bxerrh2)/log(2.)
    bxorderrhinf=log(bxerrhinfold/bxerrhinf)/log(2.)
    bxerrh1old=bxerrh1
    bxerrh2old=bxerrh2
    bxerrhinfold=bxerrhinf

    write(55,*) 'tn=',tn
    write(55,*) 'nx=',nx,'ny=',ny
    write(55,*) 'L1-rh=',bxerrh1,'  order=',bxorderrh1
    write(55,*) 'Linf-rh=',bxerrhinf,'  order=',bxorderrhinf
    write(55,*) 'L2-rh=',bxerrh2,'  order=',bxorderrh2


    byerrh1=4.*(byerrh1/nx)/ny
    byerrh2=sqrt(4.*(byerrh2/nx)/ny)
    byorderrh1=log(byerrh1old/byerrh1)/log(2.)
    byorderrh2=log(byerrh2old/byerrh2)/log(2.)
    byorderrhinf=log(byerrhinfold/byerrhinf)/log(2.)
    byerrh1old=byerrh1
    byerrh2old=byerrh2
    byerrhinfold=byerrhinf

    write(66,*) 'tn=',tn
    write(66,*) 'nx=',nx,'ny=',ny
    write(66,*) 'L1-rh=',byerrh1,'  order=',byorderrh1
    write(66,*) 'Linf-rh=',byerrhinf,'  order=',byorderrhinf
    write(66,*) 'L2-rh=',byerrh2,'  order=',byorderrh2


    bzerrh1=4.*(bzerrh1/nx)/ny
    bzerrh2=sqrt(4.*(bzerrh2/nx)/ny)
    bzorderrh1=log(bzerrh1old/bzerrh1)/log(2.)
    bzorderrh2=log(bzerrh2old/bzerrh2)/log(2.)
    bzorderrhinf=log(bzerrhinfold/bzerrhinf)/log(2.)
    bzerrh1old=bzerrh1
    bzerrh2old=bzerrh2
    bzerrhinfold=bzerrhinf

    write(77,*) 'tn=',tn
    write(77,*) 'nx=',nx,'ny=',ny
    write(77,*) 'L1-rh=',bzerrh1,'  order=',bzorderrh1
    write(77,*) 'Linf-rh=',bzerrhinf,'  order=',bzorderrhinf
    write(77,*) 'L2-rh=',bzerrh2,'  order=',bzorderrh2


    perrh1=4.*(perrh1/nx)/ny
    perrh2=sqrt(4.*(perrh2/nx)/ny)
    porderrh1=log(perrh1old/perrh1)/log(2.)
    porderrh2=log(perrh2old/perrh2)/log(2.)
    porderrhinf=log(perrhinfold/perrhinf)/log(2.)
    perrh1old=perrh1
    perrh2old=perrh2
    perrhinfold=perrhinf

    write(88,*) 'tn=',tn
    write(88,*) 'nx=',nx,'ny=',ny
    write(88,*) 'L1-rh=',perrh1,'  order=',porderrh1
    write(88,*) 'Linf-rh=',perrhinf,'  order=',porderrhinf
    write(88,*) 'L2-rh=',perrh2,'  order=',porderrh2

    diverrh1=4.*(diverrh1/nx)/ny
    diverrh2=sqrt(4.*(diverrh2/nx)/ny)
    divorderrh1=log(diverrh1old/diverrh1)/log(2.)
    divorderrh2=log(diverrh2old/diverrh2)/log(2.)
    divorderrhinf=log(diverrhinfold/diverrhinf)/log(2.)
    diverrh1old=diverrh1
    diverrh2old=diverrh2
    diverrhinfold=diverrhinf

    write(100,*) 'tn=',tn
    write(100,*) 'nx=',nx,'ny=',ny
    write(100,*) 'L1-rh=',diverrh1,'  order=',divorderrh1
    write(100,*) 'Linf-rh=',diverrhinf,'  order=',divorderrhinf
    write(100,*) 'L2-rh=',diverrh2,'  order=',divorderrh2

    write(200,*) ktime, globaldivergence

    write(*,*)
    write(*,*) 'nx=',nx,'ny=',ny
    write(*,*) 'L1-rh=',rhoerrh1,'  order=',rhoorderrh1
    write(*,*) 'Linf-rh=',rhoerrhinf,'  order=',rhoorderrhinf
    write(*,*) 'L2-rh=',rhoerrh2,'  order=',rhoorderrh2

    !close(11)
    !close(22)
    !close(33)
    !close(44)
    !close(55)
    !close(66)
    !close(77)
    !close(88)
    !close(99)
    !close(00)
    !close(100)

    return

    end subroutine outp