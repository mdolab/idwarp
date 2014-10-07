module gridInput
  ! 
  ! A module to hold the input data for the mesh movement
  !

  use constants
  implicit none
  save
  real(kind=realType) :: aExp, bExp
  real(kind=realType):: LdefFact, alpha
  integer(kind=intType)::iSymm
  real(kind=realType) :: symmTol
  real(Kind=realType) :: errTol
  logical :: useRotations
  integer(kind=intType) :: evalMode
  integer(kind=intTYpe), parameter :: EVAL_EXACT = 0
  integer(kind=intTYpe), parameter :: EVAL_FAST = 1
end module gridInput
