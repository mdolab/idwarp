module cgnsData
  !
  ! A module to hold the various data from the CGNS file needed for reading 
  ! and writing
  !
  use precision
  implicit none
  include 'cgnslib_f.h'
  save

  character*32 :: baseName
  integer(kind=intType), dimension(3)::sizes
end module cgnsData
