!---------------------------------------------------------------------!
    module fitter
    implicit none
    double precision :: pi=3.141592653589793238462643383279502884197d0
    contains 
!---------------------------------------------------------------------!
    function cosine(n,y) result(c) 
    implicit none
    integer,                        intent(in) :: n    ! input array dimension
    double precision, dimension(n), intent(in) :: y    ! input array
    double precision, dimension(0:n-1)         :: c    ! output coefficients
    integer                                    :: i,j  ! indexes
    double precision                           :: cnst
    cnst = pi / dble(2 * n)
    c = 0.d0
    c(0) = sum(y) / sqrt(2.d0)
    do i=1,n-1
     do j=1,n
      c(i)=c(i)+y(j)*cos((2 * j - 1) * i * cnst)
     enddo
    enddo
    c = c * 2.d0 / n
    endfunction cosine 
!---------------------------------------------------------------------!
    endmodule fitter
!---------------------------------------------------------------------!
     program TDCT
     use fitter
     implicit none
     double precision, allocatable, dimension(:) :: x,y,c,dx_work
     double precision                            :: x0,dx,cnst
     integer                                     :: i, j, n, m
     character(8)                                :: rubbish
     call getarg(1,rubbish)
     if(rubbish.eq."--help") then
      write(*,100)
      write(*,200)
      write(*,300) "Truncated Discrete Cosine Tranform Generator:"
      write(*,200)
      write(*,300) "Usage: - Forward Mode"
      write(*,300) "./TDCT.x <n> <m> < <file>"
      write(*,300) " -> where <n> is the number of points in the input vector"
      write(*,300) " -> where <m> is the number of points in the requested TDCT"
      write(*,300) " -> where <file> is the file containing the input vector"
      write(*,300) "       NOTE: this must have two comment lines at the top of the file"
      write(*,200)
      write(*,300) "Usage: - Reverse Mode"
      write(*,300) "./TDCT.x -r < <file>"
      write(*,300) " -> where <file> is the file containing the TDCT"
      write(*,300) "       NOTE: this must have two lines containing info on the"
      write(*,300) "             construction of the TDCT (auto-generated in forward mode"
      write(*,200)
      write(*,300) "Author:"
      write(*,300) "Thomas Pope"
      write(*,300) "thomas.pope2@newcastle.ac.uk"
      write(*,300) "Chemistry, School of Natural and Environmental Sciences,"
      write(*,300) "Newcastle University, NE1 7RU, Newcastle Upon Tyne, UK"
      write(*,200)
      write(*,100) 
     elseif(rubbish.eq."-r") then
      read(*,*) rubbish
      read(*,*) rubbish, n, m, x0, dx
      allocate(x(n), y(n), c(0:n-1))
      do j=0,m
       read(*,*) rubbish, c(j)
      enddo
      forall(j=1:n) x(j) = (j - 1) * dx
      cnst = pi / dble(2 * n)
      y = c(0) / sqrt(2.d0)
      do j=1,n
       do i=1,m
        y(j)=y(j) + c(i) * cos((2 * j - 1) * i * cnst)
       enddo
      enddo
      write(*,101)
      write(*,*) "# "
      do j=1,n
       write(*,301) x(j)+x0, y(j)
      enddo
     else
      read(rubbish,'(i5)') n
      call getarg(2,rubbish)
      read(rubbish,'(i5)') m
      if(m.gt.n-1) stop 'm is too large'
      if(m.lt.1) stop 'm is too small'
      allocate(x(n), y(n), c(0:n-1))
      read(*,*) rubbish
      read(*,*) rubbish
      do j=1,n
       read(*,*) x(j), y(j)
      enddo
      x0 = x(1)
      dx = (x(n) - x(1)) / (n - 1)
      x  = (x - x0)
      c  = cosine(n,y)       
      write(*,102) 100.d0*sum(abs(c(m+1:)))
      write(*,202) n, m, x0, dx
      do j=0,m
       write(*,302) j, c(j)
      enddo
     endif
100  format("+----------------------------------------------------------------------+")
200  format(t1,"|",t72,"|")
300  format(t1,"|",1x,a,t72,"|")
101  format(t1,"# Reverse Truncated Discrete Cosine Tranform")
301  format(t1,2(1x,es20.10))
102  format(t1,"# Truncated Discrete Cosine Tranform Generated with a Predicted Median Percentage Error Below",1x,f15.10)
202  format(t1,"#",2(1x,i0),2(1x,es20.10))
302  format(t1,i0,1x,es20.10)
     endprogram TDCT
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!

