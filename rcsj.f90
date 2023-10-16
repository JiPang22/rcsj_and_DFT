program aa
implicit none 
integer j,jmax,n,nmax
real*8 t,dt,phi,dphi,i,sumr,sumi,w,beta,phi0,ic,r,c,phidot,dphidot
!i = normalize current / phi = cooper phase / w = frequancy / jmax = maximum time index / nmax = max. fre. 
!dt = time interval / phi0 = fluxon / r = resistance / ic = critical current / c = capacitance
parameter(dt=1e-3,phi0=2.07e-15,r=1.e3,ic=7.5e-6,c=10.e-9,jmax=1000)
!parameter(beta=6.28**2*ic*r**2*c/phi0)
parameter(beta=4.nmax=jmax)
real*8 xt(jmax+1)

!create file
open(1,file='aa')
open(2,file='bb')

!initial condition
t=0.
phi=1/12.
phidot=1/12.
i=1.5

do j=1,jmax
t=j*dt

!eq. and update
dphidot=i/beta -6.28*phidot/beta
dphi=phidot
phidot=phidot+dphidot*dt
phi=phi+dphi*dt

!recode history
write(1,*) t,phi
write(*,*) t,phi
xt(j)=phi
enddo

!DFT
do n=0,nmax
w=6.28*n/(jmax*dt)

!parameter reset
t=0.
sumr=0.
sumi=0.
!time integral
do j=1,jmax
t=j*dt
sumr=sumr +xt(j)*cos(w*t)*dt
sumi=sumr -xt(j)*sin(w*t)*dt
enddo
write(2,*) w,sqrt(sumi**2+sumr**2)
enddo
end
