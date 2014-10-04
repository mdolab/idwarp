!
!      ******************************************************************
!      *                                                                *
!      * File:          constants.F90                                   *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 05-23-2010                                      *
!      * Last modified: 05-23-2010                                      *
!      *                                                                *
!      ******************************************************************
!
module constants
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Definition of the constants used in the code.                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use precision
  implicit none
  save

  ! Maximum number of families
  integer(kind=intType), parameter :: maxFamilies = 32

  ! Maximum numbers of characters in a string and in a cgns name

  integer(kind=intType), parameter :: maxStringLen   = 256
  integer(kind=intType), parameter :: maxCGNSNameLen =  32

  ! Numerical constants

  real(kind=realType), parameter :: pi    = 3.1415926535897931_realType
  real(kind=realType), parameter :: phi = 1+sqrt(5.0_realType)/2.0_realType
  ! Floating point parameters.
  real(kind=realType), parameter :: zero  = 0.0_realType
  real(kind=realType), parameter :: one   = 1.0_realType
  real(kind=realType), parameter :: two   = 2.0_realType
  real(kind=realType), parameter :: three = 3.0_realType
  real(kind=realType), parameter :: four  = 4.0_realType
  real(kind=realType), parameter :: five  = 5.0_realType
  real(kind=realType), parameter :: six   = 6.0_realType
  real(kind=realType), parameter :: eight = 8.0_realType
  reaL(kind=realType), parameter :: ten   = 10.0_realType

  real(kind=realType), parameter :: half   = 0.5_realType
  real(kind=realType), parameter :: third  = one/three
  real(kind=realType), parameter :: fourth = 0.25_realType
  real(kind=realType), parameter :: quarter= 0.25_realType
  real(kind=realType), parameter :: sixth  = one/six
  real(kind=realType), parameter :: eighth = 0.125_realType


end module constants
