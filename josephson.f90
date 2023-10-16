program aa
implicit none
integer*8 j,jmax,n,nmax,m
real*8 alpha,phi,phidot,dphi,dphidot,t,i,dt,di,vsum

open(1,file='aa')
!alpha=1.66e-4
alpha=1.e4
phi=1./12.
phidot=1./12.
dt=1.e-2
di=1.e-2
jmax=500
nmax=500
m=1/2

do j=-jmax,jmax
i=j*di

do n=-nmax,nmax
t=n*dt
vsum=vsum+phidot
dphidot=-alpha*phidot-sin(6.28*phi)*cos(3.14*m)+i
dphi=phidot
phidot=phidot+dphidot*dt
phi=phi+dphi*dt
enddo
write(1,*) i, vsum/(2*jmax)
vsum=0.
enddo
end
