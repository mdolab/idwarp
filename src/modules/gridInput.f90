module gridInput
    !
    ! A module to hold the input data for the mesh movement
    !

    use constants
    implicit none
    save
    real(kind=realType) :: aExp, bExp
    real(kind=realType) :: LdefFact, alpha
    real(kind=realType) :: symmTol
    real(kind=realType) :: errTol
    logical :: useRotations
    integer(kind=intType) :: evalMode
    real(kind=realType) :: cornerAngle
    logical :: zeroCornerRotations
    integer(kind=intTYpe), parameter :: EVAL_EXACT = 0
    integer(kind=intTYpe), parameter :: EVAL_FAST = 1
end module gridInput
