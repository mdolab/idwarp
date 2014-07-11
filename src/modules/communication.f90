module communication

  use precision
  implicit none
  save

  ! warp_comm_world: The communicator of this processor group.
  ! warp_comm_self : The single processor communicator 
  ! myID:            My processor number in warp_comm_world.
  ! nProc:           The number of processors in warp_comm_world.

  integer(kind=intType) :: warp_comm_world, warp_comm_self, myID, nProc

end module communication
