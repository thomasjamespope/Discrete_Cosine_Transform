!---------------------------------------------------------------------!
    module fitter
    implicit none
    double precision :: pi=3.141592653589793238462643383279502884197d0
    contains 
!---------------------------------------------------------------------!
    function forward_cosine(n,y) result(c) 
    implicit none
    integer,                        intent(in) :: n    ! input array dimension
    double precision, dimension(n), intent(in) :: y    ! input array
    double precision, dimension(0:n-1)         :: c    ! output coefficients
    integer                                    :: i,j  ! indexes
    double precision                           :: cnst
    cnst = pi / dble(2 * n)
    c(0) = sum(y) / sqrt(2.d0)
    c(1:)=[(sum([(y(j)*cos((2 * j - 1) * i * cnst),j=1,n)]),i=1,n-1)]
    c = c * 2.d0 / n
    endfunction forward_cosine 
!---------------------------------------------------------------------!
    function reverse_cosine(n,y) result(c) 
    implicit none
    integer,                            intent(in) :: n    ! input array dimension
    double precision, dimension(0:n-1), intent(in) :: y    ! input array
    double precision, dimension(n)                 :: c    ! output coefficients
    integer                                        :: i,j  ! indexes
    double precision                               :: cnst
    cnst = pi / dble(2 * n)
    c = y(0) / sqrt(2.d0) + [(sum([(y(i) * cos((2 * j - 1) * i * cnst),i=1,n-1)]),j=1,n)]
    endfunction reverse_cosine 
!---------------------------------------------------------------------!
    endmodule fitter
!---------------------------------------------------------------------!
     program DCT
     implicit none
     integer                     :: i, istat
     character(200), dimension(2) :: arg(2)=""
     call get_command_argument(1,arg(1),status=istat)
     if (istat.ne.0) stop "ERROR::Please select a file"
     call get_command_argument(2,arg(2),status=istat)
     if(arg(1).eq."--help".or.arg(2).eq."--help") then
      write(*,100)
      write(*,200)
      write(*,300) "Discrete Cosine Tranform Generator:"
      write(*,200)
      write(*,300) "Usage: - Forward Mode"
      write(*,300) "./DCT.x <file>"
      write(*,300) " -> where <file> is the file containing the input vector"
      write(*,200)
      write(*,300) "Usage: - Reverse Mode"
      write(*,300) "./DCT.x --r <file>"
      write(*,300) " -> where <file> is the file containing the input vector"
      write(*,200)
      write(*,300) "Author:"
      write(*,300) "Thomas Pope"
      write(*,300) "thomas.pope2@newcastle.ac.uk"
      write(*,300) "Chemistry, School of Natural and Environmental Sciences,"
      write(*,300) "Newcastle University, NE1 7RU, Newcastle Upon Tyne, UK"
      write(*,100) 
     elseif (istat.eq.0) then
      if(arg(1).eq."--r".or.arg(2).eq."--r") then
       if(arg(1).eq."--r") arg(1)=arg(2)
       call cosine_transform(arg(1),.true.)
      else
       stop "ERROR::Problem with the arguments"
      endif
     else
      call cosine_transform(arg(1),.false.)
     endif
100  format("+----------------------------------------------------------------------+")
200  format(t1,"|",t72,"|")
300  format(t1,"|",1x,a,t72,"|")     
     endprogram DCT
!---------------------------------------------------------------------!
     subroutine cosine_transform(fnm,reverse)
     use fitter
     character(200), intent(in)                  :: fnm
     logical, intent(in)                         :: reverse
     integer, parameter                          :: u=60
     integer                                     :: n=-1,ios=0,i
     character(2)                                :: dum
     double precision, allocatable, dimension(:) :: v
     open(u,file=fnm)
     do while (ios.eq.0)
      n = n+1
      read(u,*,iostat=ios) dum
     enddo
     rewind(u)
     allocate(v(n))
     read(u,*) v
     close(u)
     if(reverse) then
      write(*,100) reverse_cosine(n,v)       
     else
      write(*,100) forward_cosine(n,v)       
     endif
100  format(t1,es20.10)
     endsubroutine
!---------------------------------------------------------------------!

