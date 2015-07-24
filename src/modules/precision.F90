module precision
  ! Define precision variables for pyWarpUStrct. We ensure at least 12
  ! digits for the real and standand integers. 

  implicit none
  save 

  integer, parameter :: realType = selected_real_kind(12)
  integer(kind=4), private :: dummyInt
  integer, parameter :: intType = kind(dummyInt)

end module precision
