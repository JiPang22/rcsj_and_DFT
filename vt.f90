program aa
implicit none
integer j
real*8 i,t,dt,x,v,dx,dv,t0,x0,phi,phidot,dphi,dphidot,g,q0,om0,omc,kapa,r,i0,phi0,l,b,mr,voltage,kr,c
parameter(r=8.,l=90.e-6,mr=2.3e-12,b=80.e-3,q0=7000.,i0=1.e-6,x0=0.1e-15,t0=1.e-6,kr=1.5,phi0=2.07e-15,c=10.e-9,dt=1.e-8)
parameter(omc=2.*3.14*r*i0/phi0,kapa=l*b/phi0,g=l*b/mr,om0=sqrt(kr/mr))
open(1,file='aa')

x=1.;v=0.;phi=1/12.;phidot=0.;i=1

do j=0,1e6
t=j*dt

dphidot=omc*t0**2*i/(r*c*2*3.14)-t0*omc*phidot/(r*c*om0)-omc*t0**2*sin(2*3.14*phi)*cos(kapa*x0*x)
dphi=phidot

phidot=phidot+dphidot*dt
phi=phi+dphi*dt

dv=g*(t0**2)*i/(2*x0)-(om0**2)*(t0**2)*x-om0*t0*v/q0
dx=v

v=v+dv*dt
x=x+dx*dt

voltage=phi0*phidot/(t0*i0*r)

write(*,*) '1st term', g*(t0**2)*i/(2*x0)
write(*,*) '2nd term', (om0**2)*(t0**2)*x
write(*,*) '3th term', om0*t0*v/q0
enddo
end
