! *****************************************************************
! MATLAB® mex-interface for the IL2017 model.                     *
!                                                                 *
! By Juho Iipponen (Finnish Meteorological Institute)             *
!                                                                 *
! Building:                                                       *
! mex FCFLAGS="\$FCFLAGS -O3" -output IL2017_mex IL2017_coeff_mod.f90 IL2017_mod.f90  IL2017_mex.F90                                                                                                        *
!                                                                 *
! Usage (See README for further instructions):                    *
! [rho, temperature, components] = IL2017_mex(day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average, ae_integrals);  
! *****************************************************************

#include "fintrf.h"
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
    use IL2017_mod
    implicit none

    mwPointer :: plhs(*), prhs(*) ! Matlab input and output arguments
    integer :: nlhs, nrhs

    ! Matlab utility functions for converting between Fortran and Matlab formats
    mwPointer :: mxCreateDoubleScalar
    mwPointer :: mxGetPr, mxIsDouble, mxGetM, mxGetN
    mwPointer :: mxCreateDoubleMatrix
    mwPointer :: output_ptr, ae_int_ptr
    mwSize    :: output_size
    real(kind = 8) :: mxGetScalar
    integer(kind = 4) :: mexPrintf
    character(len = 80) :: mess
    integer :: k

    ! Inputs
    real(kind = 8) :: day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average
    real(kind = 8) :: ae_integrals(4)

    ! Outputs
    real(kind = 8) :: rho, temperature, components(8)
    
    ! Check input and output argument counts
    if(nrhs /= 8) then
       call mexErrMsgTxt('IL2017: Exactly 8 inputs required!')
    endif

    if(nlhs > 3) then
       call mexErrMsgTxt('IL2017: Too many outputs!')
    endif

    ! Read input values
    if(mxIsDouble(prhs(1)) .eq. 0) call mexErrMsgTxt('IL2017: day_of_year has to be a double precision scalar &
                                                        (MEX will cast it internally to correct int!)')
    day_of_year = mxGetScalar(prhs(1))
    if(mxIsDouble(prhs(2)) .eq. 0) call mexErrMsgTxt('IL2017: altitude has to be a double precision scalar!')
    altitude = mxGetScalar(prhs(2))
    if(mxIsDouble(prhs(3)) .eq. 0) call mexErrMsgTxt('IL2017: latitude has to be a double precision scalar!')
    latitude = mxGetScalar(prhs(3))
    if(mxIsDouble(prhs(4)) .eq. 0) call mexErrMsgTxt('IL2017: longitude has to be a double precision scalar!')
    longitude = mxGetScalar(prhs(4))
    if(mxIsDouble(prhs(5)) .eq. 0) call mexErrMsgTxt('IL2017: UTC_hour has to be a double precision scalar!')
    UTC_hour = mxGetScalar(prhs(5))
    if(mxIsDouble(prhs(6)) .eq. 0) call mexErrMsgTxt('IL2017: F30 has to be a double precision scalar!')
    F30 = mxGetScalar(prhs(6))
    if(mxIsDouble(prhs(7)) .eq. 0) call mexErrMsgTxt('IL2017: F30_average has to be a double precision scalar!')
    F30_average = mxGetScalar(prhs(7))
    if(mxIsDouble(prhs(8)) .eq. 0) call mexErrMsgTxt('IL2017: ae_integral has to be a double precision vector!')
    ae_int_ptr = mxGetPr(prhs(8))

    if (mxGetM(prhs(8))*mxGetN(prhs(8)) /= 4) then
        call mexErrMsgTxt('IL2017: ae_integrals vector must contain exactly 4 elements!')
    endif

    call mxCopyPtrToReal8(ae_int_ptr, ae_integrals, 4)

    ! Now call IL2017 *****************************************************
    call IL2017(int(day_of_year), altitude, latitude, longitude, UTC_hour, F30, F30_average,&
                         ae_integrals, rho, temperature, components)                     
    ! **********************************************************************

    plhs(1) = mxCreateDoubleScalar(dble(rho))
    plhs(2) = mxCreateDoubleScalar(dble(temperature))
    
    ! The third output is a matrix of number concentrations and temperature variables
    plhs(3) = mxCreateDoubleMatrix(8,1,0) ! 8 by 1 real matrix
    output_ptr = mxGetPr(plhs(3))
    output_size = 8
    call mxCopyReal8ToPtr(dble(components), output_ptr, output_size) 

    return
end subroutine mexfunction