subroutine readMatrixFromFile(matrix)
   ! opening the file for reading
    use, intrinsic :: iso_fortran_env
    real(kind=real64), dimension(36) :: reada
    real(kind=real64), dimension(6, 6) :: matrix
    open (2, file = '../data/tasm.txt', status = 'old')

    read(2,*) reada
    matrix=RESHAPE(reada, shape(matrix)  )
   close(2)
   
end subroutine readMatrixFromFile

program demodist
      use, intrinsic :: iso_fortran_env
      use mod_dist 
      implicit none


      integer i, npart, strlen, a0, z0;
      real(kind=real64), dimension(6) :: coordinates 
      real(kind=real64), dimension(36) :: tas 
      real(kind=real64), dimension(6) :: closed
      real(kind=real64), dimension(2500,6) :: distribution1
      real(kind=real64), dimension(1000000) :: x,xp,y, yp, sigma, deltap
      real(kind=real64), dimension(1000000) :: xn, pxn,yn, pyn, zn, zpn
      real(kind=real64) momentum, mass, one, e1,e2, e3, dp, betx1, zero, pia2, six
      real(kind=real64), dimension(6, 6) :: identity, results, testm
      real(kind=real64) energy0, mass0
      
      character(len=256) filename, fileout
      
      open (2, file = '../data/tasm.txt', status = 'old')

      read(2,*) tas

      e1 = 1.0d0
      e2 = 2.0d0
      e3 = 0.0d0
      dp = 0.000001d0
      pia2 = 2.00d0*3.1415
      zero = 0.0d0
      energy0 =200
      mass = 0.93800
      one =1d0
      six = 6.000d0
      closed(1)=10
      npart = 1000000
      do i=1, npart
      xn(i) = dp*i
      pxn(i) = dp*i
      yn(i) = dp*i*dp
      pyn(i) = dp*i
      zn(i) = dp*i*dp
      zpn(i) = dp
      enddo
      

      npart =1000

      ! Initialize 3 distributions with dimenstion 6

      call dist_initializedistribution(1)
    ! Set the tas matrix 
         
      call dist_setemitt12(e1,e2)
      call dist_setemitt3(e3)
      call dist_settasmatrix(tas)

      call dist_sete0andmass0(energy0, mass0)
      
      !call dist_settotalsteps(npart)
      !call dist_setscan_para_diagonal(0,0,6,zero,one);
      !call dist_setscan_para_diagonal(1,0,4,zero,pia2);
      !call dist_setscan_para_diagonal(2,0,6,zero,one);
      !call dist_setscan_para_diagonal(3,0,4,zero,pia2);
      !call dist_setscan_para_diagonal(4,0,6,zero,one);
      !call dist_setscan_para_diagonal(5,0,4,zero,pia2);

     ! call dist_setscan_para_grid(0,0,6,zero,one,2);
     ! call dist_setscan_para_grid(1,0,4,zero,pia2,2);
     ! call dist_setscan_para_grid(2,0,6,zero,one,3);
     ! call dist_setscan_para_grid(3,0,4,zero,pia2,3);
     ! call dist_setscan_para_grid(4,0,6,zero,one,3);
     ! call dist_setscan_para_grid(5,0,4,zero,pia2,30);
      call dist_setcooords(xn,pxn,yn,pyn,zn,zpn, npart, 1)
      call dist_get6trackcoord(x,xp,y,yp,sigma,deltap, npart)
   
      do i=1, npart
      print *, x(i),xp(i),y(i),yp(i),sigma(i),deltap(i) 
      enddo

    !Distribution 2: a matched distribution

    ! Change the distribution to 1
!    call setdistribution(1)
!    call settasmatrix(tas)
!    call setemittance12(e1,e2)
!    call setemittance3(e3)
!    call setmassmom(mass, momentum)!

    !call setparameter(1,zero,one,100,6);
    !call setparameter(2,zero,pia2,100,4);
    !call setparameter(3,zero,one,100,6);
    !call setparameter(4,zero,pia2,100,4);
    !call setparameter(5,zero,zero,1,0);
    !call setparameter(6,zero,zero,1,0);



!e1 from fort.10 1.999999999648927940E+00
!e2 from fort.10 9.999999998426114534E-01
 
      end program demodist