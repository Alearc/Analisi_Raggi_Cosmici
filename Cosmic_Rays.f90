program es5
implicit none
real,dimension(:),allocatable:: t,histo,ts
real:: tm,errtm
real,dimension(:),allocatable::q2,th1,th2
real,dimension(5)::u,ni
real::tmin,tmax,dx
integer::n,m,k,l
integer::ind,nbin,i,fuori
real:: s11,s20,s10,s01,s00
real::d,coeffm,q
real::em,eq,cov
real,dimension(5)::x,y,ei
real,dimension(15000)::like
real,dimension(10)::p,calc

!___________________________________________________________MEDIA ARITMETICA 
n=171
allocate (t(n))

open(unit=1,action='read')
 do i=1,n
read(unit=1,fmt=*)t(i)
 end do
 close(unit=1) 
 
tm=sum(t)/(n*1.0)        !tutte le stime vanno moltiplicate per 10^-4
errtm=tm/sqrt(n*1.0)
print*, 'stima con media aritmetica',tm,'con err', errtm  

!___________________________________________________________ISTOGRAMMA 
tmin=minval(t)
tmax=maxval(t)
nbin=10
dx=(tmax-tmin)/nbin
allocate(ts(nbin))
allocate(histo(nbin))
histo=0
fuori=0

do i = 1,n
ind = int((t(i)-tmin)/dx) + 1
if ( ind>nbin .or. ind <1)then
fuori = fuori + 1
else
histo(ind) = histo(ind) + 1
end if
end do


do i=1,nbin  
ts(i)=tmin+(i-1)*dx+dx/2            !t_star: valore medio del singolo bin
write(2,*) i, ts(i), (histo(i)), sqrt(histo(i))
end do

!plot [:][:] 'fort.2'u  2:3 w boxes
!replot 'fort.2'u  2:3:4 w err
!t1=8330.62500
!f1(x)=exp(-x/t1)*200*1666.12500/t1
!replot f1(x)
!t2=6175
!f2(x)=exp(-x/t2)*200*1666.12500/t2
!replot f2(x)
!t3=6263.75342
!f3(x)=exp(-x/t3)*200*1666.12500/t3
!replot f3(x)
!t4=5012
!f4(x)=exp(-x/t4)*200*1666.12500/t4
!replot f4(x)


print*, "dx:", dx


!_________________________________________________MMQ

ni=histo
l=180000
allocate(q2(l))
allocate(th1(l))
allocate(th2(l))

do k=1,l   
u=0
  do i=1,5
  u(i)=n*exp(-(ts(i)/k))*dx/k  !valore di aspettazione della binomiale per l'isto
  end do

q2(k)=sum(((ni-u)**2)/ni)
write(3,*) k, q2(k)
end do

print*, 'Q2 MIN', minval(q2), 'Q2 MIN+1', minval(q2)+1

!plot [5000:7000][:] 'fort.3'u  1:2 w l, 7.23371840,8.23371887   anzi 1.322298, 2.3229747

print*, 'stima con MMQ', minloc(q2) !'con err',634 (io ho trovato 619.5)
!__________________________________________________MMQL
do i=1,5
y(i)=log(ni(i))
ei(i)=1./sqrt(ni(i)*1.)     !errore su yi (dev std) con prop. var.
x(i)=ts(i)
write(4,*) x(i),y(i),ei(i)
end do
 !uso il metodo delle sommatorie perchè ho una relaz. lin.
s11=sum(x*y/ei**2)
s20=sum(x**2/ei**2)
s10=sum(x/ei**2)
s01=sum(y/ei**2)
s00=sum(1./ei**2)
d=s00*s20-s10**2

coeffm=(s00*s11-s01*s10)/d             !m = -1/tau
q=(s01*s20-s11*s10)/d                  !q = ln(tau) + c 

em=sqrt(s00/d)   !migliore stima di sigma_m 
eq=sqrt(s20/d)   !migliore stima di sigma_q
cov=s10/d
print*, 'stima con MMQL', -1./coeffm,'con err', em*(1./coeffm)/coeffm !di nuovo con prop. var.

!_____________________________________________MML
m=10
calc=0
do k=1,15000
  do i=1,m
  p(i)=exp(-ts(i)/k)*dx/(k*(1-exp(-maxval(ts)/k))) !normalizzazione perchè l'integrale non va a infinito
  !p(i)=exp(-ts(i)/k)*dx/k
   calc(i)=log(p(i))*ni(i)
end do

like(k)=sum(calc)
write(11,*) k, like(k)
end do

print*, 'maximum Likelihood', maxval(like),'maximum Likelihood-0.5',maxval(like)-0.5

!plot [:][:] 'fort.11'u  1:2 w l


!vengono ora scritti i valori corrispondenti alle x cercando sul file (tasto ctrl + f) oppure per via grafica,faccio una media sull'errore
print*,'stima con MML', maxloc(like) !, 'con err',387





end program es5
