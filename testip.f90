program es5_testip
implicit none
real,dimension(:),allocatable:: t,histo,ts
real:: tm,errtm
real,dimension(:),allocatable::q2,th1,th2
real,dimension(5)::u,ni
real::tmin,tmax,dx
integer::n,m,k,l
integer::ind,nbin,i,fuori
real::d,coeffm,q,t1,t2,t3,t4
real::em,eq,cov
real,dimension(5) :: chi1, chi2, chi3, chi4, u1, u2, u3, u4
real :: test1, test2, test3, test4, e1, e2, e3, e4
real,dimension(5)::x,y,ei

n=171
allocate (t(n))
t1 = 5935.58496
t2 = 5865
t3 = 5919.60303
t4 = 4976

open(unit=1,action='read')
 do i=1,n
read(unit=1,fmt=*)t(i)
 end do
 close(unit=1) 

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

do i = 1, 5
 u1(i)=(162/t1)*exp(-ts(i)/t1)*dx         !valori di aspettazione dei singoli bin
 u2(i)=(162/t2)*exp(-ts(i)/t2)*dx
 u3(i)=(162/t3)*exp(-ts(i)/t3)*dx
 u4(i)=(162/t4)*exp(-ts(i)/t4)*dx
 
write(7,*)u1(i),u2(i),u3(i),u4(i) 
end do


!________________________________________ TEST D'IPOTESI DI CHI QUADRO
do i = 1, 5               !confrontare i valori con il grado di libertà 3 (5-2)
 chi1(i) = (((histo(i)-u1(i))**2)/(u1(i)))  
 chi2(i) = (((histo(i)-u2(i))**2)/(u2(i)))  
 chi3(i) = (((histo(i)-u3(i))**2)/(u3(i)))  
 chi4(i) = (((histo(i)-u4(i))**2)/(u4(i)))  
end do
test1 = sum(chi1)
test2 = sum(chi2)
test3 = sum(chi3)
test4 = sum(chi4)

print*, "test1:", test1
print*, "test2:", test2
print*, "test3:", test3
print*, "test4:", test4

write(4,*)test1,test2,test3,test4


end program es5_testip
