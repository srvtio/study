program rhoprogram
  implicit none
  integer i,j
  real(8) x1,x2
  real(8), parameter :: l1 = 1.0d0
  real(8), parameter :: l2 = 0.5d0
  integer, parameter :: n1 = 100
  integer, parameter :: n2 = 100
  real(8), parameter :: dx1 = 2.0d0*l1 / dble(n1)
  real(8), parameter :: dx2 = 2.0d0*l2 / dble(n2)
  
  real(8) rho1(n1,n2),rho2(n1,n2)
  real(8), parameter :: t0  = 1.0d0
  real(8), parameter :: ts0 = 1.0d0
  real(8), parameter :: ts1 = 2.0d0
  real(8), parameter :: rho0 = 1.0d0 
  real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
  real(8), parameter :: c = rho0*sqrt(t0)/(pi**1.5d0)

!--- ファイルの出力 ---  
  open(10,file="rho.dat")

  do j=1,n2
     do i=1,n1

        x1 = -l1 + dble(i-1) * dx1
        x2 = -l2 + dble(j-1) * dx2
        
        rho1(i,j) = c * sqrt(pi) / sqrt(t0) &
             * ( atan((l2-x2)/(l1-x1)) + 2.0d0*pi - atan(-(l2-x2)/(l1+x1)) ) &
             + c * sqrt(pi) / sqrt(ts0) &
             * ( atan((l2-x2)/(0.5d0*l1-x1)) - atan((l2-x2)/(l1-x1)) + atan(-(l2-x2)/(l1+x1)) - atan((l2-x2)/(-0.5d0*l1-x1))) &
             + c * sqrt(pi) / sqrt(ts1) &
             * ( atan((l2-x2)/(-0.5d0*l1-x1)) - atan((l2-x2)/(0.5d0*l1-x1)) )

        x2 = -x2

        rho2(i,j) = c * sqrt(pi) / sqrt(t0) &
             * ( atan((l2-x2)/(l1-x1)) + 2.0d0*pi - atan(-(l2-x2)/(l1+x1)) ) &
             + c * sqrt(pi) / sqrt(ts0) &
             * ( atan((l2-x2)/(0.5d0*l1-x1)) - atan((l2-x2)/(l1-x1)) + atan(-(l2-x2)/(l1+x1)) - atan((l2-x2)/(-0.5d0*l1-x1))) &
             + c * sqrt(pi) / sqrt(ts1) &
             * ( atan((l2-x2)/(-0.5d0*l1-x1)) - atan((l2-x2)/(0.5d0*l1-x1)) )

        x2 = -x2

        write(10,'(3(f12.8))') x1,x2,rho1(i,j),rho2(i,j),rho1(i,j)+rho2(i,j)

     end do
     write(10,*)
  end do

  close(10)

end program rhoprogram
