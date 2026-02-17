subroutine mexFunction(nlhs, plhs, nrhs, prhs)
use complex_prolate_swf
use param
use iso_c_binding
implicit none
integer(c_int), parameter :: mxREAL = 0
integer(c_int), parameter :: mxCOMPLEX = 1

!   Inputs (nrhs = 8):
!       1) cc      : complex size parameter (double complex scalar)
!       2) m       : order (integer scalar)
!       3) lnum    : number of degrees (integer scalar)
!       4) ioprad  : radial flag (0/1/2)
!       5) x1      : x - 1 (double scalar)
!       6) iopang  : angular flag (0/1/2)
!       7) iopnorm : normalization flag (0/1)
!       8) arg     : eta values (double vector, length narg)
!
!   Outputs (nlhs up to 15):
!       1) r1c    (complex)  2) ir1e  (double, integer-like)
!       3) r1dc   (complex)  4) ir1de (double, integer-like)
!       5) r2c    (complex)  6) ir2e  (double, integer-like)
!       7) r2dc   (complex)  8) ir2de (double, integer-like)
!       9) naccr  (double, integer-like)
!       10) s1c    (complex) 11) is1e  (double, integer-like)
!       12) s1dc   (complex) 13) is1de (double, integer-like)
!       14) naccs  (double, integer-like)
!       15) naccds (double, integer-like)
!
!   Each function value is returned as characteristic + exponent
!       value ≈ characteristic * 10.^exponent

      integer nlhs, nrhs
      type(c_ptr) plhs(*), prhs(*)

    interface
        function mxGetPr(pa) bind(C, name="mxGetPr") result(p)
            import :: c_ptr
            type(c_ptr), value :: pa
            type(c_ptr) :: p
        end function mxGetPr
        function mxGetPi(pa) bind(C, name="mxGetPi") result(p)
            import :: c_ptr
            type(c_ptr), value :: pa
            type(c_ptr) :: p
        end function mxGetPi
        function mxGetNumberOfElements(pa) bind(C, name="mxGetNumberOfElements") result(n)
            import :: c_ptr, c_size_t
            type(c_ptr), value :: pa
            integer(c_size_t) :: n
        end function mxGetNumberOfElements
        function mxIsDouble(pa) bind(C, name="mxIsDouble") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: pa
            integer(c_int) :: r
        end function mxIsDouble
        function mxIsComplex(pa) bind(C, name="mxIsComplex") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: pa
            integer(c_int) :: r
        end function mxIsComplex
        function mxCreateDoubleMatrix(m, n, flag) bind(C, name="mxCreateDoubleMatrix") result(pa)
            import :: c_ptr, c_size_t, c_int
            integer(c_size_t), value :: m, n
            integer(c_int), value :: flag
            type(c_ptr) :: pa
        end function mxCreateDoubleMatrix
        subroutine mexErrMsgTxt(msg) bind(C, name="mexErrMsgTxt")
            import :: c_char
            character(kind=c_char), dimension(*) :: msg
        end subroutine mexErrMsgTxt
    end interface

    integer :: m, lnum, ioprad, iopang, iopnorm, narg
    integer(c_size_t) :: lnum_mw, narg_mw
    real(knd) :: x1, cc_r, cc_i
    complex(knd) :: cc
    real(knd), allocatable :: arg(:)

    complex(knd), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    integer, allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    complex(knd), allocatable :: s1c(:,:), s1dc(:,:)
    integer, allocatable :: is1e(:,:), is1de(:,:), naccs(:,:), naccds(:,:)

    integer :: nout

    if (nrhs /= 8) call mexErrMsgTxtF('Expected 8 inputs: cc, m, lnum, ioprad, x1, iopang, iopnorm, arg.')
    if (nlhs >= 15) call mexErrMsgTxtF('Too many outputs. Max 15.')

    if (knd /= kind(1.0d0)) call mexErrMsgTxtF('param.knd must be double for MEX.')

    call read_scalar_real(prhs(1), cc_r, 'cc')
    if (mxIsComplex(prhs(1)) .ne. 0) then
        call read_scalar_imag(prhs(1), cc_i)
    else
        cc_i = 0.0_knd
    endif
    cc = cmplx(cc_r, cc_i, kind=knd)

    call read_scalar_int(prhs(2), m, 'm')
    call read_scalar_int(prhs(3), lnum, 'lnum')
    call read_scalar_int(prhs(4), ioprad, 'ioprad')
    call read_scalar_real(prhs(5), x1, 'x1')
    call read_scalar_int(prhs(6), iopang, 'iopang')
    call read_scalar_int(prhs(7), iopnorm, 'iopnorm')

    if (lnum .lt. 1) call mexErrMsgTxtF('lnum must be >= 1.')

    if (mxIsDouble(prhs(8)) .eq. 0) call mexErrMsgTxtF('arg must be double.')
    if (mxIsComplex(prhs(8)) .ne. 0) call mexErrMsgTxtF('arg must be real.')
    narg_mw = mxGetNumberOfElements(prhs(8))
    if (narg_mw .lt. 1) call mexErrMsgTxtF('arg must be non-empty.')
    narg = int(narg_mw)

    allocate(arg(narg))
    call copy_in_real_vector(prhs(8), arg, narg)

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg), naccds(lnum, narg))

    call cprofcn(cc, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, naccds)

    lnum_mw = int(lnum, c_size_t)
    narg_mw = int(narg, c_size_t)
    nout = nlhs

    if (nout >= 1) then
        plhs(1) = mxCreateDoubleMatrix(lnum_mw, 1, mxCOMPLEX)
        call copy_complex_vector(plhs(1), r1c, lnum)
    endif
    if (nout >= 2) then
        plhs(2) = mxCreateDoubleMatrix(lnum_mw, 1, mxREAL)
        call copy_int_vector(plhs(2), ir1e, lnum)
    endif
    if (nout >= 3) then
        plhs(3) = mxCreateDoubleMatrix(lnum_mw, 1, mxCOMPLEX)
        call copy_complex_vector(plhs(3), r1dc, lnum)
    endif
    if (nout >= 4) then
        plhs(4) = mxCreateDoubleMatrix(lnum_mw, 1, mxREAL)
        call copy_int_vector(plhs(4), ir1de, lnum)
    endif
    if (nout >= 5) then
        plhs(5) = mxCreateDoubleMatrix(lnum_mw, 1, mxCOMPLEX)
        call copy_complex_vector(plhs(5), r2c, lnum)
    endif
    if (nout >= 6) then
        plhs(6) = mxCreateDoubleMatrix(lnum_mw, 1, mxREAL)
        call copy_int_vector(plhs(6), ir2e, lnum)
    endif
    if (nout >= 7) then
        plhs(7) = mxCreateDoubleMatrix(lnum_mw, 1, mxCOMPLEX)
        call copy_complex_vector(plhs(7), r2dc, lnum)
    endif
    if (nout >= 8) then
        plhs(8) = mxCreateDoubleMatrix(lnum_mw, 1, mxREAL)
        call copy_int_vector(plhs(8), ir2de, lnum)
    endif
    if (nout >= 9) then
        plhs(9) = mxCreateDoubleMatrix(lnum_mw, 1, mxREAL)
        call copy_int_vector(plhs(9), naccr, lnum)
    endif
    if (nout >= 10) then
        plhs(10) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxCOMPLEX)
        call copy_complex_matrix(plhs(10), s1c, lnum, narg)
    endif
    if (nout >= 11) then
        plhs(11) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxREAL)
        call copy_int_matrix(plhs(11), is1e, lnum, narg)
    endif
    if (nout >= 12) then
        plhs(12) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxCOMPLEX)
        call copy_complex_matrix(plhs(12), s1dc, lnum, narg)
    endif
    if (nout >= 13) then
        plhs(13) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxREAL)
        call copy_int_matrix(plhs(13), is1de, lnum, narg)
    endif
    if (nout >= 14) then
        plhs(14) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxREAL)
        call copy_int_matrix(plhs(14), naccs, lnum, narg)
    endif
    if (nout >= 15) then
        plhs(15) = mxCreateDoubleMatrix(lnum_mw, narg_mw, mxREAL)
        call copy_int_matrix(plhs(15), naccds, lnum, narg)
    endif

    deallocate(arg)
    deallocate(r1c, r1dc, r2c, r2dc)
    deallocate(ir1e, ir1de, ir2e, ir2de, naccr)
    deallocate(s1c, s1dc)
    deallocate(is1e, is1de, naccs, naccds)

    return

    contains

    subroutine read_scalar_real(arr, val, name)
        type(c_ptr), value :: arr
        real(knd) val
        character*(*) name
        integer(c_size_t) n
        type(c_ptr) p
        real(c_double), pointer :: rptr(:)

        if (mxIsDouble(arr) .eq. 0) call mexErrMsgTxtF(trim(name)//' must be double.')
        n = mxGetNumberOfElements(arr)
        if (n .ne. 1) call mexErrMsgTxtF(trim(name)//' must be scalar.')
        p = mxGetPr(arr)
        call c_f_pointer(p, rptr, (/1/))
        val = real(rptr(1), knd)
    end subroutine read_scalar_real

    subroutine read_scalar_imag(arr, val)
        type(c_ptr), value :: arr
        real(knd) val
        type(c_ptr) p
        real(c_double), pointer :: rptr(:)

        p = mxGetPi(arr)
        if (c_associated(p)) then
            call c_f_pointer(p, rptr, (/1/))
            val = real(rptr(1), knd)
        else
            val = 0.0_knd
        endif
    end subroutine read_scalar_imag

    subroutine read_scalar_int(arr, ival, name)
        type(c_ptr), value :: arr
        integer ival
        character*(*) name
        real(knd) tmp

        call read_scalar_real(arr, tmp, name)
        if (abs(tmp - nint(tmp)) .gt. 0.0_knd) then
            call mexErrMsgTxtF(trim(name)//' must be integer.')
        endif
        ival = int(tmp)
    end subroutine read_scalar_int

    subroutine copy_complex_vector(mxarr, z, n)
        type(c_ptr), value :: mxarr
        complex(knd) z(n)
        integer n
        type(c_ptr) pr, pi
        real(c_double), pointer :: rptr(:), iptr(:)

        pr = mxGetPr(mxarr)
        pi = mxGetPi(mxarr)
        call c_f_pointer(pr, rptr, (/n/))
        if (c_associated(pi)) then
            call c_f_pointer(pi, iptr, (/n/))
        endif
        rptr = real(z, c_double)
        if (c_associated(pi)) then
            iptr = aimag(z)
        endif
    end subroutine copy_complex_vector

    subroutine copy_complex_matrix(mxarr, z, n, m)
        type(c_ptr), value :: mxarr
        complex(knd) z(n, m)
        integer n, m
        integer nelem
        type(c_ptr) pr, pi
        real(c_double), pointer :: rptr(:), iptr(:)

        nelem = n * m
        pr = mxGetPr(mxarr)
        pi = mxGetPi(mxarr)
        call c_f_pointer(pr, rptr, (/nelem/))
        if (c_associated(pi)) then
            call c_f_pointer(pi, iptr, (/nelem/))
        endif
        rptr = reshape(real(z, c_double), (/nelem/))
        if (c_associated(pi)) then
            iptr = reshape(aimag(z), (/nelem/))
        endif
    end subroutine copy_complex_matrix

    subroutine copy_int_vector(mxarr, v, n)
        type(c_ptr), value :: mxarr
        integer v(n)
        integer n
        type(c_ptr) pr
        real(c_double), pointer :: rptr(:)

        pr = mxGetPr(mxarr)
        call c_f_pointer(pr, rptr, (/n/))
        rptr = real(v, c_double)
    end subroutine copy_int_vector

    subroutine copy_int_matrix(mxarr, v, n, m)
        type(c_ptr), value :: mxarr
        integer v(n, m)
        integer n, m
        integer nelem
        type(c_ptr) pr
        real(c_double), pointer :: rptr(:)

        nelem = n * m
        pr = mxGetPr(mxarr)
        call c_f_pointer(pr, rptr, (/nelem/))
        rptr = real(reshape(v, (/nelem/)), c_double)
    end subroutine copy_int_matrix

    subroutine copy_in_real_vector(arr, out, n)
        type(c_ptr), value :: arr
        real(knd) out(n)
        integer n
        type(c_ptr) p
        real(c_double), pointer :: rptr(:)

        p = mxGetPr(arr)
        call c_f_pointer(p, rptr, (/n/))
        out = real(rptr, knd)
    end subroutine copy_in_real_vector

    subroutine mexErrMsgTxtF(msg)
        character(len=*), intent(in) :: msg
        character(kind=c_char), allocatable :: cmsg(:)
        integer :: i, n

        n = len_trim(msg)
        allocate(cmsg(n+1))
        do i = 1, n
            cmsg(i) = msg(i:i)
        enddo
        cmsg(n+1) = c_null_char
        call mexErrMsgTxt(cmsg)
    end subroutine mexErrMsgTxtF

end subroutine mexFunction

