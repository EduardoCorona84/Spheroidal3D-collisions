                               complex_prolate_swf

  Cprofcn is available as both a subroutine version provided as the
  module complex_prolate_swf and a stand alone version cprofcn. It was
  originally developed by arnie lee van buren about 2005 with technical
  support from jeffrey boisvert. Version 1.04 added the capability
  to calculate the radial functions of the first kind and their first
  derivatives for x = 1.0. 

  Table of Contents
  1. Purpose
  2. Introduction
  3. Input and Output
  4. Accuracy of results
  5. Obtaining the expansion d coefficients
  6. Obtaining the eigenvalues

 1. Purpose

  To calculate the first and second kind prolate radial functions r1
  and r2 and their first derivatives r1d and r2d for a given order m,
  a range of degrees l beginning at l = m and for a specific complex size
  parameter c and shape parameter x. To calculate the first kind prolate
  angular functions and their first derivatives with respect to eta for
  the same values of m, l and c and for a set of values of the angular
  coordinate eta. The subroutine version complex_prolate_swf calculates
  values for a single input value of m. The stand alone version calculates
  values for a range of values of m.

 2. Introduction

  Cprofcn is written in free format fortran. It is designed around the
  maximum number of decimal digits ndec and the maximum exponent nex
  available in real arithmetic. Procedures used in profcn allow for
  exponents much larger than nex since the resulting floating point
  radial function values are given as a characteristic and an integer
  exponent.

  Cprofcn can be run in double precision, quadruple precision or a hybrid
  where the Bouwkamp procedure to refine the eigenvalues is run in quadruple
  precision while the remainder of the calculations are performed in double
  precision. In the latter case, cprofcn switches to quadruple precision for
  the Bouwkamp procedure whenever double precision fails to provide highly
  accurate eigenvalues. See the discussion below about when to choose which
  arithmetic. The choice is set in the module param provided in the github
  repository. If this is not available, then create param as follows:

    module param
    integer, parameter :: knd = selected_real_kind(8)
    integer, parameter :: knd1 = selected_real_kind(8) 
    logical, parameter :: debug = .true.
    logical, parameter :: warn = .true.
    logical, parameter :: output = .true.
    logical, parameter :: suffix = .true.
    end module param

  Set the value of knd in the parenthesis to either 8 for double precision
  or 16 for quadruple precision arithmetic. Set the value of knd1 to be used
  for the Bouwkamp procedure. Note that if knd = 16, knd1 should also be 16.
  The logicals in param are described below in the discussion of the output
  files.

  Some computers may have more than 8 bytes for double precision
  data and more than 16 bytes for quadruple precision data or may use
  kind values that do not correspond to the number of bytes. In this
  case just use the appropriate integers for the kind parameters in
  module param. Also change the values of kindd and kindq set in
  statement 5 below the comments section to the kind values for
  double precision data and quadruple precision data, respectively.

  This program for complex c is a major advance in capability
  from cprofcn version 1.01. It is written in fortran 90. It was
  developed around profcn version 1.08 which obtains spheroidal
  function values for real values of the size parameter c = kd/2,
  where k is the wavenumber and d is the interfoocal distance of
  the elliptical cross-section of the proate spheroid. A description
  of the methods used in profcn is provided in two articles: (1)
  A. L. Van Buren and J. E. Boisvert, "Accurate calculation of prolate
  spheroidal radial functions of the first kind and their first
  derivatives," Quart. Appl. Math. 60 (2002), 589-599 and (2) A. L.
  Van Buren and J. E. Boisvert, "Improved calculation of prolate
  spheroidal radial functions of the second kind and their first
  derivatives,' Quart. Appl. Math. 62 (2004), 493-507.
  
  I plan to write a manuscript describing the methods used in cprofcn
  and submit it to arXiv.org. In the meantime the user may want
  to look at the article 'Calculation of oblate spheroidal wave
  functions with complex argument,' available at arXiv.org, identifier
  2009.01618, August 2020. Some of the methods used in cprofcn were
  first developed for use in the program coblfcn that is described in
  this article.
  
  Expressions to calculate the radial functions of the first kind and its
  first derivatives for x = 1.0 when m = 0 are obtained from the limits of
  equations given in the first reference (with c now complex). The radial
  functions of the first kind and their first derivatives are equal to 0.0
  when m is unequal to 0. [The corresponding radial functions of the second
  kind and their first derivatives are infinite for all m.]

  Cprofcn provides function values for c complex = real(c) + i aimag(c)
  = cr + i ci, where the imaginary part ci often accounts for losses in
  wave propagation. Ci is assumed positive in cprofcn. If the user has
  a negative value for ci, just run cprofcn with ci positive instead
  and take the complex conjugate of the results, including the function
  values, eigenvalues, expansion coefficients, and normalization
  factors.

  3. Input and Output

  Following is a description of the input and output parameters in the
  call statement for the subroutine version. After that will be a
  description of the the input and output files associated with the
  stand alone version. Note that these output files, if desired, can
  also be obtained when running the subroutine version. See comments about
  this below.

  A sample input and resulting output from profcn is provided by the
  files cprofcndat (text version of the input file cprofcn.dat for the
  stand alone version), cprofort20 (text version of the output file
  fort.20 giving the resulting radial functions) and cprofort30 (text
  version of the output file fort.30 giving the resulting angular
  functions).

  Subroutine Version of cprofcn

    subroutine profcn(c,m,lnum,ioprad,x1,iopang,iopnorm,narg,arg, &
                      r1c,ir1e,r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr, &
                      s1c,is1e,s1dc,is1de,naccs,naccds)

        complex(knd), intent(in)   ::  cc
        real(knd), intent (in)     ::  x1, arg(narg)
        integer, intent (in)       ::  m, lnum, ioprad, iopang, iopnorm, narg
        complex(knd), intent (out) ::  r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                                       s1c(lnum, narg), s1dc(lnum, narg)
        integer, intent (out)      ::  ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                                       is1e(lnum, narg), is1de(lnum, narg), & 
                                       naccr(lnum), naccs(lnum, narg), naccds(lnum, narg)

      Input and output parameters appearing in the subroutine call
      statement are defined below:

          cc     : desired complex value of the size parameter c (= kd/2,
                   where k =  the complex wavenumber and d = interfocal
                   length) [real(knd)]
          m      : desired value for the order m (integer)
          lnum   : number of values desired for the degree l equal
                   to m, m + 1, m + 2, ..., m + lnum - 1 (integer)
                   if lnum is less than 2*[real(c)+aimag(c)]/pi it
                   should be an even integer.
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed
          x1     : x - 1, where x is the radial coordinate [real(knd)]
                   [When x1 = 0.0, i.e., x = 1.0, only radial functions
                   of the first kind and their first derivatives
                   are calculated. Radial functions of the second
                   kind and their first derivtives are infinite.
                   Thus ioprad must be set equal to 1. Note
                   also that the radial functions of the first kind
                   and their first derivatives are equal to zero
                   unless m = 0.] 
          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed
          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner and
                      Schafke normalization scheme.] This norm
                      becomes very large as m becomes large. The
                      angular functions are computed below as
                      a characteristic and an exponent to avoid
                      overflow.
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.
                      This is very useful since it removes the
                      need to calculate a normalization factor
                      when using the angular function values given
                      here. It also eliminates any chance for
                      overflow when the characteristics and exponents
                      are combined to form the angular functions.
          narg   : number of values of the angular coordinate eta for
                   which angular functions are calculated (integer)
          arg    : vector containing the narg values of eta for which
                   angular functions are desired [real(knd)]
          r1c/   : vectors of length lnum containing the complex characteristics
          r1dc     for the radial functions of the first kind r1 and their
                   first derivatives [complex(knd)]
          ir1e/  : integer vectors of length lnum containing the
          ir1de    exponents corresponding to r1c and r1dc
          r2c/   : vectors of length lnum containing the complex characteristics
          r2dc     for the radial functions of the second kind r2 and their
                   first derivatives [complex(knd)]
          ir2e/  : integer vectors of length lnum containing the
          ir2de    exponents corresponding to r2c and r2dc
          naccr  : integer vector of length lnum containing the estimated
                   accuracy of the radial functions
          s1c,   : two-dimensional arrays s1c(lnum,narg) and
          s1dc     s1dc(lnum,narg) that contain narg calculated
                   complex characteristics for the angular functions and
                   their first derivatives for each of the lnum
                   values of l [complex(knd)]
                   For example, s1(10,1) is the characteristic
                   of the angular function for l = m +10 -1 and
                   for the first value of eta given by arg(1)
          is1e,    integer arrays is1e(lnum,narg) and is1de(lnum,narg)
          is1de    containing the exponents corresponding to s1c and
                   s1dc
          naccs  : two-dimensional array naccs(lnum,narg) that contains
                   narg estimated accuracy values for the angular functions
                   for each of the lnum values of l
          naccds : two-dimensional array naccds(lnum,narg) that contains
                   narg estimated accuracy values for first derivatives
                   of the angular functions for each of the lnum values of l
                   
  Stand alone Version of profcn

     Input Data

     Input parameters are read from unit 1 in the file cprofcn.dat
     assumed to be in the directory of cprofcn.f90. cprofcn.dat
     contains the following lines of data:

       line 1:
          mmin   : minimum value for m. (integer)
          minc   : increment for m. (integer)
          mnum   : number of values of m. (integer)
          lnum   : number of values of l [l=m, l=m+1,
                   ..., l=m+lnum-1]. (integer)
                   if lnum is less than 2*[real(c)+aimag(c)]/pi
                   it should be an even integer. If lnum is chosen
                   to be odd, profcn will increase lnum by one.

       line 2:
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed

          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed

          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner-
                      Schafke normalization scheme.]
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.

       line 3:
          cc     : value of the complex size parameter c (= kd/2, where k =
                   the complex wavenumber and d = interfocal length) [complex(knd)]
          x1     : value of the radial coordinate x minus one [real(knd)]
                   (a nominal value of 10.0e0_knd can be entered for x1
                   if ioprad = 0). If x1 = 0, i.e., x = 1.0, only radial
                   functions of the first kind and their first derivatives
                   can be calculated. Radial functions of the second kind
                   and their first derivatives are infinite. Thus ioprad
                   must be set equal to 1. Note also that the radial functions
                   of the first kind and their first derivatives are equal
                   to 0 unless m is equal to 0.]

       line 4:
          ioparg : (integer)
                 : =0 if both arg1 and darg are angles in degrees
                 : =1 if arg1 and darg are dimensionless values of eta

          arg1   : first value for the angle coordinate (in degrees
                   or dimensionless if eta) for which angular
                   functions are to be computed. [real(knd)]

          darg   : increment used to calculate additional desired
                   arguments for angular functions. [real(knd)]

          narg   : number of desired angle arguments. (integer)
                   (line 4 is not read when iopang = 0)

     Output files

     These output files are also available using the subroutine version
     complex_prolate_swf. Generation of each of the files is controlled by a
     logical specified in the module param. False suppresses the output file
     and true enables it. The logical debug controls fort.30 and fort.40,
     the logical output controls fort.20 and fort.30 and warn controls fort.60.
     The logical suffix controls whether the accuracy estimates given in fort.20
     are followed by a letter designating how the accuracy was determined. 'w'
     indicates it is based on the Wronskian and 'e' indicates it is based on
     subtraction errors involved in the calculations. Setting suffix = false
     suppresses the letter. 

   fort.20

     This file contains values for all radial functions that have
     been calculated.
     The first line in the file contains the values for x, c, and
     m, formatted as follows (see statements 120 and 130 in
     subroutine main):

                x      : e23.14 in 64 bit arithmetic; e39.30
                       : in 128 bit arithmetic
                c      : (e23.14,e23.14) in 64 bit arithmetic;
                       : (e39.30,e39,30) in 128 bit arithmetic
                m      : i5
 
     Each subsequent line in fort.20 contains radial functions
     for given values of l. The first line contains values for l = m,
     the next for l=m+1 and continuing to l=m+lnum-1. The radial
     functions are preceeded by the value of l and followed by the
     accuracy, equal to the estimated number of accurate decimal digits
     in the radial functions as measured using the Wronskian or
     estimated using subtraction errors in the calculations. [see
     comments below regarding naccr]. When ioprad = 1 and r2 and r2d are
     not calculated, an estimate of the accuracy of r1 and r1d is given
     in decimal digits at the end of the line containing their values. 

       The output and corresponding format for each line is as follows
       (see statements 690 and 700 in main).

         l            : value for l (i5)
         r1c(l-m+1)   : complex characteristic of the prolate radial
                        function of the first kind r1 for
                        the given value of l (f17.14,f17.14)
         ir1e(l-m+1)  : exponent of r1 (i6)
         r1dc(l-m+1)  : complex characteristic of the first derivative
                        of the prolate radial function of the
                        first kind r1d (f17.14,f17.14)
         ir1de(l-m+1) : exponent of r1d (i6)
         r2c(l-m+1)   : complex characteristic of the prolate radial
                        function of the second kind r2 for
                        the given value of l (f17.14,f14.7)
         ir2e(l-m+1)  : exponent of the prolate radial function of
                        second kind (i6). If the exponent for any
                        of the functions is greater than 9999,
                        the format can be increased to i7 or
                        higher. Note that the procedures used in
                        this program allow for exponents much
                        larger than those allowed in the arithmetic
                        used for the calculations since the
                        floating point function values provided
                        are given as a characteristic and an
                        integer exponent. Use of ratios in
                        calculating the functions eliminates
                        overflow and underflow during the
                        calculations.
         r2dc(l-m+1)  : complex characteristic of the first derivative
                        of the prolate radial function of second kind
                        r2d (f17.14,f17.4)
         ir2de(l-m+1) : exponent of r2d (i6). [See comment above
                        for ir2e.]
         naccr(l-m+1) : estimated accuracy: usually equal to the
                        number of decimal digits of agreement
                        between the theoretical Wronskian and the
                        calculated Wronskian (i2). This is
                        indicated by the letter w following the
                        integer naccr. r1 and r1d are usually very
                        accurate so that naccr relects the accuracy
                        of r2 and r2d.

                        Several situations are described below where
                        the Wronskian can not be used to estimate the
                        accuracy. Other factors are used here to
                        estimate naccr. This is indicated by using
                        the letter e instead of w following the value
                        for naccr in fort.20.

                        Sometimes near the breakpoint value for l,
                        the Legendre function expansions for r2 and r2d
                        converge with acceptable subtraction error but
                        the leading coefficient for the series (joining
                        factor) is highly inaccurate. Here the Wronskian
                        can sometimes be used for small ci to obtain
                        improved accuracy for this coefficient. The
                        accuracy is estimated using subtraction errors
                        involved in the calculations and the estimated
                        accuracy of the eigenvalue in those cases where
                        the Bouwkamp eigenvalue routine does not converge
                        fully.

                        When ci becomes large, there are often values
                        of l where both r1 and r2 and their first
                        derivatives are large in magnitude. Here values
                        for r2 and r2d are given by -i times values of
                        r1 and r1d. The accuracy is estimated using
                        the magnitudes of r1 and r1d, the magnitude
                        of the theoretical value of the wronskian, and
                        the estimated accuracy of r1 and r1d.

                        Sometimes the calculation of r2 and r2d using the
                        variable eta method involves a sufficiently accurate
                        numerator but a relatively inaccurate denominator.
                        Here a sufficiently accurate denominator is obtained
                        using the Wronskian. The accuracy is estimated using
                        subtraction error estimates in the calculations.                               

                        There is one other situations where the value
                        for naccr can be an estimated value. This is
                        when the accuracy using the Legendre is estimated
                        using subtraction errors and the accuracy estimate
                        provided by cprofcn is taken to be the smaller of
                        this and the Wronskian estimate. Here this designated
                        by w rather than e since the Wronskian accuracy sets
                        an upper bound.

   fort.30

     This file fort.30 contains values for all angular functions
     that have been calculated. Its first line contains the values
     for c and m, formatted as follows (see statements 65 and 70 in
     subroutine main).

                c      : (e23.14,e23.14) in 64 bit arithmetic;
                       : (e39.30,e39,30) in 128 bit arithmetic
                m      : i5

     The second line in fort.30 contains the value for the first l
     (=m), formatted as follows (see statement 140 in subroutine
     main):

                l      : i6

     This is followed by a series of narg lines. Each line contains
     a desired value of angle (ioparg = 0) or angular coordinate eta
     (ioparg =1) followed by the corresponding angular functions and
     accuracy. Specific output and format for each line is as follows.

        for iopang = 1:

               arg    : for ioparg = 0, angle in degrees (f19.14; see
                        statement 750 in subroutine main)
            or barg   ; for ioparg = 1, angular coordinate eta
                        (f17.14)
               s1c    : complex characteristic of the prolate angular
                        function of first kind (f17.14,f17.14)
               is1e   : exponent of the prolate angular function of
                        first kind (i5))

        for iopang = 2, each line also includes:
               s1dc   : complex characteristic of the first derivative
                        of the prolate angular function of first kind
                        (f17.14,f17.14); see statement 1470 in
                        subroutine main)
               is1de  : exponent of the first derivative of the
                        prolate angular function of first kind (i5;
                        see statement 1470 in subroutine main)

        for iopang = 1:
               naccs  : accuracy: estimate of the number of decimal
                        digits of accuracy in the angular function
                        (i2). It is a conservative estimate based on
                        the calculated subtraction error in the
                        Legendre function series for the angular
                        functions and the estimated accuracy of the
                        normaliation factor. When the accuracy estimate
                        is equal to 0, the corresponding angular
                        functions are set equal to zero. (i2; see
                        statements 1460 and 1470 in subroutine main).
                        The calculated angular functions tend to be
                        less accurate the larger the value of cr, the
                        smaller the value of l - m (for values less
                        than approximately 2c divided by pi), and the
                        closer eta is to zero (i.e., the closer theta
                        is to 90 degrees).

        for iopang = 2:
              naccds  : accuracy: includes an estimate of the number
                        of decimal digits of accuracy in the first
                        derivative of the angular functions (i2)   

   fort.40 and fort.50

     These files are diagnostic files that contain information
     about specific techniques used and numbers of terms required
     for the radial function and angular function calculations,
     respectively. They are annotated and should be self explanatory.

   fort.60

     This file may be of interest to the user, especially when using
     this program outside the range for which cprofcn is expected to
     provide useful results. Whenever the estimated accuracy falls below
     a designated integer value during the running of this program,
     the associated values of x, c, m, and l are written to fort.60.
     The integer is currently set equal to 6 in the write statements
     for this file found in subroutine main after the statement numbered 186
     for ioprad = 1 and after the statement numbered 710 for ioprad = 2.
     Whenever cprofcn determines that two or more of the eigenvalues of
     the same parity are duplicates, indicating a failure of the Bouwkamp
     routine to converge to the correct eigenvalue, the values of m,l and
     c where this occurs are written to fort.60.

  4. Accuracy of Results

  The following discussion is provided to help the user choose which
  arithmetic option to use. If the compiler does not support quadruple
  arithmetic, then the only option is to use double precision. If
  the compiler does support quadruple precision arithmetic, then it is
  recommended that the hybrid version be used instead of the entirely
  double precision version whenever double precision is expected to
  provide adequate accuracy. Use of quadruple precision arithmetic
  for the Bouwkamp procedure increases the run time somewhat but the
  program is still reasonably fast and the results are more accurate.
  The improvement in accuracy increases as ci increases. If the input
  parameters are outside those for which the hybrid version is expected
  to provide useful results, then the only choice is to use quadruple
  precision if available. Note that this increases the run time
  considerably.

  An integer called minacc is used in cprofcn. This designates the
  minimum number of accurate digits desired for the radial spheroidal
  functions of the second kind. The value of minacc controls which
  methods are used to calculate the radial functions. Minacc is set
  equal to 8 for double precision. It is recommended that this not
  be changed. For quadruple precision, minacc is set equal to 15 digits
  for values of ci up to 20. This should provide 15 or more digits of
  accuracy here. For ci > 20, minacc is set equal to 8 digits. Minacc
  can be increased in this case but higher accuracy might not always
  be achieved. Also the computation time will likely go up with larger
  values of minacc. Minacc is set below following these introductory
  comment statements.

  cprofcn was tested extensively using a laptop pc and a Fortran
  compiler that provides approximately 15 decimal digits in double
  precision arithmetic and approximately 31 digits in quadruple
  precision arithmetic. The estimated accuracy of the resulting
  function values is given below in terms of decimal digits. It is
  often the Wronskian comparison result but can be an estimate
  based on subtraction errors in the calculation. See the discussion
  above for the accuracy integer naccr in the section for fort.20.
  If the user's computer provides a different number of digits, the
  following estimates should be adjusted up or down depending on
  whether more or fewer digits are provided. Testing included values
  of x ranging from 1.000001 to 5.0, values for cr up to 5000, and
  values of ci up to 80. Results for x greater than 5.0 are expected
  to be similar to those for x = 5.0. Testing for the double precision
  version included all values of the order m from 0 to 200 and from
  210 to 1000 in steps of 10. The same was the case for the quadruple
  precision version when cr was less than about 500. When cr was larger
  than 500, testing for quadruple precision included values of m from
  0 to 500 in steps of 10 and from 500 to 1000 in steps of 50. It also
  included values of m from 1 to 9 for larger values of x, ci and cr,
  since the lowest accuracy often occurred here. The values of the
  degree l ranged from m to m + lnum -1, where lnum was chosen large
  enough for magnitudes of r1 and r1d to be less than 10**(-300).

  In the following discussion the term useful results means that the
  estimated accuracy for the radial functions observed during testing
  never fell below 5 decimal digits unless otherwise stated. I expect
  there are many applications where occasional 5 digit results are
  acceptable. Possibly even an isolated 4 digit result is acceptable.
  Note that there is no guarantee that the estimated accuracy for the
  values for parameter values other than those that I tested will be
  as high as I report below. The following tables provide the maximum
  value of cr = real(c) for which useful results are obtained for a
  range of values of x and ci = aimag(c). The first table applies to
  double precision without the use of quadruple precision for the
  Bouwkamp procedure. The maximum values for cr are not much different
  when using the hybrid version unless ci is greater than about 15
  and x is greater than about 1.11. However, the overall accuracy
  with the hybrid version can be considerably higher, especially near
  the breakpoint. When the maximum value for cr is somewhat larger
  for the hybrid version, it is shown in parenthesis next to the
  value for the entirely double precision version. When ci becomes much
  larger than 25 for double precision, the ranges shrink so markedly
  that quadruple precision is required to obtain useful results unless
  cr is small. If the user is willing to accept a rare 4 digit result,
  cprofcn can be used with cr somehwat larger than the values given
  below. For example, with x = 1.005 and ci = 20, using cr = 250 gives
  only 1 result with an estimated accuracy of 4 digits, none less than
  4 digits. Note that estimated accuracies of less than 5 digits
  observed for cr larger than the values given below tend to occur at
  values of m somewhat less than 200, often less than 100.

     Double precision

     aimag(c)      0        10        15          20         25

        x           approximate upper limit to real(c)

     1.000001    5000      5000      5000        1500        300
     1.00001     5000      5000      5000         800        190
     1.00005     3500      3500      3500         570        130
     1.0001      2500      2500      2500         510        120
     1.0005      1200      1200      1150         300         90
     1.001       1000      1000       990         250         60
     1.002        800       750       730         280         50
     1.003        650       600       600         240         80
     1.004        600       550       530         220         70
     1.005        550       500       490         190         70
     1.006        500       450       450         240         60
     1.007        480       450       420         270        100
     1.008        450       400       380         320         90
     1.009        450       400       370         300         80
     1.01         400       350       350         300         70
     1.02         350       300       250         220         90
     1.03         300       250       240         200         70
     1.04         290       250       210         190        100
     1.05         290       250       220         190        100
     1.06         290       240       200         190         90
     1.07         290       240       230         190         70
     1.08         300       270       210         200         90
     1.09         350       300       220         190         80
     1.10        4000**     290       210         210        100
     1.105       5000      5000      5000        3500        120
     1.11        5000      5000      5000        3500        130
     1.50        5000      5000      5000        1000(3500)  230(450)
     5.00        5000      5000      5000         300(3500)  100(300)

  ** There is a single 4 digit result that is not near a root for c
     = 5000 when x = 1.10 and when x = 1.102 and for c = 3000 and 5000
     when x = 1.103.

     Quadruple precision

     aimag(c)      0        40        50        60        70       80

        x            Approximate upper limit to real(c)

     1.000001    5000       5000      950       200       75       60 
     1.00001     5000       5000      900       170       60       60
     1.00005     5000       4500      700       150       60       60
     1.0001      5000       4000      650       135       60       60
     1.0005      3100       2500      600       145       60       60
     1.001       2900       2000      650       145       60       60
     1.002       2100       1700      550       150       70       60 
     1.003       1700       1400      650       170       65       60
     1.004       1600       1300      650       170       65      110
     1.005       1450       1100      650       165      145      115
     1.006       1350       1000      650       220      145      115
     1.007       1300       1000      650       220      145      115 
     1.008       1200        950      650       220      145      115
     1.009       1150        900      700       220      145      115
     1.01        1000        850      700       220      145      115
     1.02         900*       650*     600       220      145      115
     1.03         800*       600*     550       220      145      115
     1.04         800*       600*     550       220      145      115
     1.05         900*       600*     550       220      145      115
     1.06        1000*       700*     550       220      145      115
     1.07        2000*       800*     700       220      145      115
     1.08        5000       1200*     950       220      145      115 
     1.09        5000       1500*     950       220      145      115
     1.10        5000       1500*     950       220      145      115
     1.105       5000       1500*     950       220      145      115
     1.11        5000       5000      950       220      145      115
     1.50        5000       5000      950       220      145      115  
     5.0         5000       3500      450       210      140      110

  *  For quadruple precision arithmetic and for x larger than about
     1.02, much larger values of cr than those given in the table above
     can provide at least 5 digits of accuracy for values of m up to a
     maximum value less than 1000 but possibly large enough for many
     applications. I summarize this by listing the largest value of m
     for specific values of x and cr that provide 5 or more digits of
     accuracy up to the specified value of m. Folllowing the value of
     x is a series of values of cr followed in parenthesis by the
     maximum value of m for that value of cr and ci. This behavior is
     not seen for cr >= 50.

        ci = 0
     x = 1.02: 1000(90); 1500(130); 2000(170); 3000(230); 4000(290);
               5000(350)
     x = 1.03: 1000(120); 1500(180); 2000(240); 3000(330); 4000(420);
               5000(490)
     x = 1.04: 1000(180); 1500(250); 2000(310); 3000(440); 4000(550);
               5000(650)
     x = 1.05: 1000(190); 1500(310); 2000(420); 3000(550); 4000(700);
               5000(850)
     x = 1.06: 1100(300); 1500(390); 2000(500); 3000(750); 4000(900);
               5000(1000)
     x = 1.07: 3000(950); 4000(1000); 5000(1000)

        ci = 40
     x = 1.02: 800(60); 1000(80); 1500(110); 2000(140); 3000(210);
               4000(260);5000(310)
     x = 1.03: 800(80); 1000(110); 1500(170); 2000(220); 3000(320);
               4000(400); 5000(490)
     x = 1.04: 800(120); 1000(160); 1500(230); 2000(300); 3000(440);
               4000(550); 5000(650)
     x = 1.05: 800(140); 1000(180); 1500(280); 2000(390); 3000(550);
               4000(650); 5000(850)
     x = 1.06: 1000(240); 1500(350); 2000(470); 3000(700); 4000(850);
               5000(950)
     x = 1.07: 1000(280); 1500(450); 2000(550); 3000(850); 4000(1000);
               5000(1000)
     x = 1.08: 1500(550); 2000(800); 3000(850); 4000(1000); 5000(1000)
     x = 1.09: 2000(800); 3000(1000); 4000(1000); 5000(1000)
     x = 1.10: 2000(800); 3000(1000); 4000(1000); 5000(1000)
     x = 1.105: 2000(850); 3000(1000); 4000(1000); 5000(1000)
         
   Estimated r1 and r1d accuracy for x = 1.0  

  When x = 1.0 the resulting function values provided by profcn for r1 and r1d
  for m = 0 have an accuracy based on subtraction errors in the Flammer
  normalization constant. These subtraction errors are zero for real c and
  increase as the imaginary value of c increases. For a given c the subtraction
  errors are non-zero only when l is near the so-called breakpoint value
  lb = 2*cr/pi. The subtraction errors near l = lb are similar in magnitude
  to those that occur when computing the Flammer normalization constant for
  oblate functions when c is real with a value equal to ci. The size of these
  errors prevents the use of double precision for obtaining sufficiently
  accurate values for r1 and r1d when ci is larger than about 20. It is
  recommended that the hybrid version be used, i.e., knd = 8 and knd1 = 16,
  when ci is less than 20 and double precision is desired. Use of quadruple
  precision extends the range of ci to about 60.

   Estimated angular function accuracy

  For both choices of arithmetic, the angular functions and their first
  derivatives can lose accuracy for low l, high c, and eta near unity
  due to subtraction errors in their series calculation. However,
  their magnitude in this case is corresponding smaller than angular
  functions for higher values of l and/or eta not near zero. The loss
  in accuracy due to these subtraction errors should not adversely
  effect numerical results for physical problems using these functions.

  A second source of inaccuracy in the angular functions arises from
  subtraction errors that can occur in the calculation of their
  Meixner-Schafke normalization factor dmsnorm. Note that the loss in
  accuracy here is not in addition to other losses in accuracy for the
  angular functions but rather sets an upper limit to their accuracy.
  These errors are by far the largest for m = 0. They occur for values
  of l - m somewhat less to somewhat more than the breakpoint value
  and grow with increasing ci and to some extent with increasing cr.
  They are zero for small ci and can become as large as 6 digits for ci
  = 20, 8 digits for ci = 25, and 11 digits for ci = 30 as cr increases
  to 5000. This loss in accuracy is not likely a problem when using
  double precision arithmetic with 15 decimal digits since as ci
  becomes larger than 20, the values of cr for which the radial
  functions are accurate to 5 or more digits are progressively smaller,
  being 300 for ci = 25 and less than this for higher ci. For ci = 25
  and cr = 300 the maximum subtraction error in calculating dmsnorm
  is 7 digits. When using double precision, the Meixner and Schafke
  normalization should be accurate to at least 5 digits wherever all
  of the radial functions for a given value of m are also accurate to
  at least 5 digits. The file fort.60 mentioned above can alert the user
  whenever the estimated accuracy for the Meixner and Schafke normalization
  is less than the same integer number of decimal digits selected for alerts
  about the accuracy of the radial functions.

  For higher values of ci when using quadruple precision, the loss of
  accuracy in the normalization factor is even greater. For ci = 60,
  the loss of accuracy can be as large as 25 digits for both cr = 2000
  and cr = 5000. Note that the loss of accuracy for m = 1 is only about
  4 digits here and even less for higher values of m. For ci = 70 it is
  26 digits for cr = 1000. And for ci = 80 it is 25 digits for cr =
  400. This should not be a problem using 33 decimal digits since it
  still allows for accuracies of at least 5 digits for the angular
  functions everywhere the radial functions also have an accuracy of 5
  or more digits.

  A third source of inaccuracy in the angular functions arises from the
  potential loss of accuracy in the eigenvalues at values of l near and
  somewhat below the breakpoint when ci is not very small and ci is
  moderate to large. This is most likely to occur when using double
  precision for the calculations including for the Bouwkamp procedure.
  Use of quadruple precision, if available, for the Bouwkamp procedure
  will help considerably here. Even when using double precision for the
  Bouwkamp procedure, the eigenvalue accuracy never fell below 8 digits
  when real(c) was equal to or less than the table values shown below.
  
  5. Obtaining the d expansion coefficients

  The user may desire values for the d coefficients that appear in
  the expression for the angular functions as well as in many of the
  expressions used to calculate the radial functions. Ratios of
  successive d coefficients are stored in the vector enr where enr(k)
  = d(subscript 2k+ix) divided by d(subscript 2k-2+ix). The vector enr
  is calculated in the subroutine dnorm before statement 20 and passed
  to subroutine main. The number lim2 of d coefficients calculated for
  a given l is chosen to be sufficient to compute radial and angulalar
  functions for that l. The number of d coefficients needed to compute
  r1 and r1d and s1 and s1d can range from less than 100 to somewhat
  more than l/2 for large values of l. The number of d coefficients
  required to compute r2 and r2d are comparable to this unless they are
  computed using one of the Neumann function expansions. Here one can
  require up to 200000 or so coefficients when x is near 0.01. Note
  that the vector enr returned by the subroutine conver contains scaled
  ratios where the scaling has been chosen to produce a symmetric
  matrix for computing eigenvalues. The scaling factors are removed in
  subroutine dnorm to obtain the desired d coefficient ratios.

  The d coefficients themselves can be obtained starting with the value
  for d with the subscript l - m. If iopnorm is set = 0, cprofcn uses
  the Meixner and Schafke normalization scheme for the angular
  functions. Here the angular functions have the same norm as the
  corresponding associated Legendre functions. When c is real the
  calculation of the normalizing factor is accurate with no associated
  subraction errors. As ci increases and to some extent as cr
  increases, the subtraction errors increase for some values of l. The
  subtraction errors are greatest for m = 0 and decrease rapidly with
  increasing m. For m = 0, they are near zero for ci up to 10, then
  increase to about 6 digits at ci = 20 and about 15 digits at ci = 50.
  See the discussion of this in the section 'Estimated angular function
  accuracy' given above. The subroutine s1 computes d(subscript l-m)
  for this normalization and returns it to subroutine main as a
  characteristic dmlms and an exponent idmlmse. Use of an exponent
  avoids possible overflow of d(subscript l-m) for extremely large c
  and m. When the user sets iopnorm = 1 so that the angular functions
  have unit norm, the corresponding characteristic and exponent for
  d(subscript l-m) are calculated in subroutine s1 and returned to
  subroutine main as dmlms1 and idmlms1e. Corresponding values for
  the characteristic and exponent of d(subscript l-m) for the
  Morse and Feshbach and for the Flammer normalizations are computed
  in dnorm and returned to main as dmlmf, idmlmfe and dmlf, idmlfe. Note
  that for c complex, all three of these normalization calculations
  suffer subtraction errors for lower values of l-m and non-small c. The
  values for d(subscript l-m) will have reduced accuracy in this
  case. When c is real only the Morse and Feshbach normalization suffers
  subtraction errors in its calculation.

  6. Obtaining the eigenvalues

  The eigenvalues for the prolate functions are computed in subroutine
  conver and returned to main where they are stored in the vector eig(l+1).
  There is such a vector created for each value of m.
