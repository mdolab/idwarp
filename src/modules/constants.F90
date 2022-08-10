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
    integer(kind=intType), parameter :: maxFamilies = 128

    ! Maximum numbers of characters in a string and in a cgns name
    integer(kind=intType), parameter :: maxStringLen = 512
    integer(kind=intType), parameter :: maxCGNSNameLen = 32

    ! Numerical constants
    real(kind=realType), parameter :: pi = 3.1415926535897931_realType
    real(kind=realType), parameter :: phi = 2.1180339887498949
    real(kind=realType), parameter :: eps = 1e-15

    ! Floating point parameters.
    real(kind=realType), parameter :: zero = 0.0_realType
    real(kind=realType), parameter :: one = 1.0_realType
    real(kind=realType), parameter :: two = 2.0_realType
    real(kind=realType), parameter :: three = 3.0_realType
    real(kind=realType), parameter :: four = 4.0_realType
    real(kind=realType), parameter :: five = 5.0_realType
    real(kind=realType), parameter :: six = 6.0_realType
    real(kind=realType), parameter :: sevenp5 = 7.5_realType
    real(kind=realType), parameter :: eight = 8.0_realType
    real(kind=realType), parameter :: ten = 10.0_realType
    real(kind=realType), parameter :: twenty = 20.0_realType
    real(kind=realType), parameter :: forty = 40.0_realType
    real(kind=realType), parameter :: eighty = 80.0_realType
    real(kind=realType), parameter :: onesixty = 160.0_realType
    real(kind=realType), parameter :: threetwenty = 320.0_realType
    real(kind=realType), parameter :: sixforty = 640.0_realType

    real(kind=realType), parameter :: half = 0.5_realType
    real(kind=realType), parameter :: third = one / three
    real(kind=realType), parameter :: fourth = 0.25_realType
    real(kind=realType), parameter :: quarter = 0.25_realType
    real(kind=realType), parameter :: sixth = one / six
    real(kind=realType), parameter :: eighth = 0.125_realType

end module constants
