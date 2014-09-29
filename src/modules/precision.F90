module precision
  implicit none
  save 

  integer, parameter :: realType = selected_real_kind(12)
  integer(kind=4), private :: dummyInt
  integer, parameter :: intType = kind(dummyInt)

end module precision
