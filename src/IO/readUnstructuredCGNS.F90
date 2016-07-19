subroutine readUnstructuredCGNS(cg)

  use precision
  use communication
  use gridData
  use gridInput  
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  integer(kind=intType), intent(in) :: cg

  ! CGNS Variables
  integer(kind=intType) :: i, j, k, m, l, istart, iend, localsize, iProc
  integer(kind=intType):: ierr, base, dims(3), iZone
  integer(kind=intType):: nNodes, nCells
  integer(kind=intType) :: tmpSym, nSymm
  character*32 :: zoneName, bocoName, famName
  character*32 :: secName

  integer(kind=intType) :: nbocos, boco
  integer(kind=intType) :: nVertices,nElements, nzones
  integer(kind=intType) :: zoneType,dataType,sec,type
  integer(kind=intType) :: nSections,nElem,nConn
  integer(kind=intType) :: eBeg,eEnd,nBdry, parentFlag
  integer(kind=intType) :: bocoType,ptsettype,nbcelem
  integer(kind=intType) :: normalIndex,normalListFlag
  integer(kind=intType) :: nDataSet, tmpInt(2)
  integer(kind=intType) :: elementDataSize, curElement, eCount, nPnts
  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
  integer(kind=intType), dimension(:), allocatable :: tmpConn
#ifdef USE_COMPLEX
  complex(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
  complex(kind=realType), dimension(:,:), pointer :: elemNodes
#else
  real(kind=realType), dimension(:, :), allocatable :: allNodes, localNodes
  real(kind=realType), dimension(:,:), pointer :: elemNodes
#endif
  integer(kind=intType), dimension(:), allocatable :: surfaceNodes, localSurfaceNodes
  integer(kind=intType), dimension(:, :), allocatable :: sizes

  integer(kind=intType) :: status(MPI_STATUS_SIZE)
  type(sectionDataType), pointer :: secPtr
  integer(kind=intType), dimension(:), pointer :: elemConn, elemPtr

  integer(kind=intType) :: nElemNotFound, curElem, curElemSize, localIndex
  real(kind=realType) :: symmSum(3)
  logical :: isSurface 

  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then

     base = 1_intType

     ! read the number of zones in the file
     call cg_nzones_f(cg, base, nzones, ierr)
     if (ierr .eq. ERROR) call cg_error_exit_f

     ! allocate the zonal storage
     allocate(zones(nZones))

     ! do a check that we have all unstructured zones
     do iZone = 1,nZones
        ! Check the zone type. This should be Unstructured. 
        call cg_zone_type_f(cg,base, iZone, zoneType, ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f

        if (zoneType .ne. Unstructured) then
           print *,"Unstructured zone expected...stopping"
           stop
        end if
     end do

     ! Now count up the total number of nodes
     nNodes = 0
     nCells = 0
     allocate(sizes(3, nZones))

     do iZone=1, nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        nVertices = dims(1)
        nElements = dims(2)
        zones(iZone)%nVertices = nVertices
        zones(iZone)%nElements = nElements
        zones(iZone)%name = zoneName
        sizes(:, iZone) = dims(1:3)
        nNodes = nNodes + dims(1)
        nCells = nCells + dims(2)
     end do

     ! Now we know the total number of nodes we can allocate the final
     ! required space and read them in. 
     allocate(allNodes(3, nNodes))
     allocate(surfaceNodes(nNodes))
     surfaceNodes = 0

     ! Loop over the zones and read the nodes
     zoneLoop: do iZone = 1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        allocate(coorX(dims(1)), coorY(dims(1)), coorZ(dims(1)))

        ! Read the x,y,z-coordinates. Assume double precision
         call cg_coord_read_f(cg,base,iZone,"CoordinateX",&
             & realDouble, 1, dims(1), coorX, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        call cg_coord_read_f(cg,base,iZone,"CoordinateY",&
             & realDouble, 1, dims(1), coorY, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
           
        call cg_coord_read_f(cg,base,iZone,"CoordinateZ",&
             & realDouble, 1, dims(1), coorZ, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Now stack all of the zones in one array
        do i=1, dims(1)
#ifdef USE_COMPLEX
           allNodes(1, i) = cmplx(coorX(i), 0.0)
           allNodes(2, i) = cmplx(coorY(i), 0.0)
           allNodes(3, i) = cmplx(coorZ(i), 0.0)
#else
           allNodes(1, i) = coorX(i)
           allNodes(2, i) = coorY(i)
           allNodes(3, i) = coorZ(i)
#endif
        end do

        ! Now loop over the sections in the unstructured file and
        ! figure out the boundary nodes.

        ! Determine the total number of sections in this zone. This
        ! will involve surface sections and volume sections. We only
        ! care about the surface sections
        call cg_nsections_f(cg, base, iZone, nSections, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Allocate space for the sections
        allocate(zones(iZone)%sections(nSections))

        do sec=1, nSections
           ! Read the section and determine the element type and store
           ! some generic information.
           call cg_section_read_f(cg, base, iZone, sec, secName, &
                type, eBeg, eEnd, nBdry, parentFlag, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Nullify the elemPtr and elemConn, since it may not be allocated
           nullify(zones(iZone)%sections(sec)%elemConn, &
                   zones(iZone)%sections(sec)%elemPtr)

           ! Number of elements on this section
           nElem = eEnd - eBeg + 1
           zones(iZone)%sections(sec)%nElem = nElem
           zones(iZone)%sections(sec)%name = secName
           zones(iZone)%sections(sec)%elemStart = eBeg
           zones(iZone)%sections(sec)%elemEnd = eEnd   

           ! Currently we can only deal with three and four sided
           ! surfaces. These can show up type == TRI_3 (all tris),
           ! type == QUAD_4 or type == MIXED. In the case of mixed we
           ! need to check each of the types individually. 
           select case (type)

              ! all tris or quads
           case (TRI_3, QUAD_4)
              ! Firstly flag this section as a surface
              zones(iZone)%sections(sec)%isSurface = .True.
             
              ! Get the nodes per element "npe" (constant)
              call cg_npe_f(type, nConn, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Allocate space for the connectivities and the pointer
              allocate(zones(iZone)%sections(sec)%elemConn(nConn*nElem))
              allocate(zones(iZone)%sections(sec)%elemPtr(nElem+1))

              ! This is the actual connectivity real call.
              call cg_elements_read_f(cg, base, iZone, sec, &
                   zones(iZone)%sections(sec)%elemConn, NULL, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Set up the pointer which is simple in this case...it
              ! is just an arithematic list. 
              zones(iZone)%sections(sec)%elemPtr(1) = 1
              do i=2,nElem+1
                 zones(iZone)%sections(sec)%elemPtr(i) = &
                      zones(iZone)%sections(sec)%elemPtr(i-1) + nConn
              end do

           case (MIXED) 

              ! This section *could* be a surface but we don't
              ! actually know yet. 

              ! First get the number of values we need to read 
              call cg_ElementDataSize_f(cg, base, iZone, sec, &
                   ElementDataSize , ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Allocate a temporary for this
              allocate(tmpConn(ElementDataSize))

              ! Now read the 'connectivity'--- not really, it has the
              ! connectivity *AND* the element types.
              call cg_elements_read_f(cg, base, iZone, sec, &
                   tmpConn, NULL, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! NOW we can figure out if this is actually a surface
              ! zone. Assume it is is, and flag and stop if not.
              zones(iZone)%sections(sec)%isSurface = .True.
              i = 1
              elementCheckLoop: do while (i <= ElementDataSize)
                 curElement = tmpConn(i)
                 select case (curElement)
                 case (TRI_3, QUAD_4)
                    ! We're OK...now determine how much we need to
                    ! increment the counter

                    call cg_npe_f(curElement, nConn, ierr)
                    if (ierr .eq. CG_ERROR) call cg_error_exit_f
                    
                    i = i + nConn + 1 ! The plus 1 is for the element Type

                 case default
                    ! Not a surface or we can't read it.
                    zones(iZone)%sections(sec)%isSurface = .False.
                    exit elementCheckLoop
                 end select
              end do elementCheckLoop
              
              ! Only continue if we are now sure the section is only
              ! surfaces
              if (zones(iZone)%sections(sec)%isSurface) then 
                               
                 ! Allocate space for the connectivities and the
                 ! pointer. It should be save to take this as
                 ! ElementDataSize - nElem. 
                 allocate(zones(iZone)%sections(sec)%elemConn(ElementDataSize-nElem))
                 allocate(zones(iZone)%sections(sec)%elemPtr(nElem+1))
                 
                 zones(iZone)%sections(sec)%elemConn = -999999999

                 ! Loop back over the section elements
                 ! again. Essentially we are taking the information in
                 ! tmpConn and splitting up into elemConn and elemPtr
                 i = 1
                 eCount = 1
                 zones(iZone)%sections(sec)%elemPtr(eCount) = 1
                 elementLoop: do while (i <= ElementDataSize)
                    curElement = tmpConn(i)
                    select case (curElement)
                    case (TRI_3, QUAD_4)
                       ! We're OK...now determine how much we need to
                       ! increment the counter
                       
                       call cg_npe_f(curElement, nConn, ierr)
                       if (ierr .eq. CG_ERROR) call cg_error_exit_f
                       
                       ! iStart is the start of the element in elemPtr
                       ! iEnd is one more than the end of the element
                       ! in elemPtr (or the start of the next one)

                       iStart = zones(iZone)%sections(sec)%elemPtr(eCount)
                       iEnd   = iStart + nConn
                       zones(iZone)%sections(sec)%elemPtr(eCount + 1) = iEnd

                       ! The connectivity is easy:
                       zones(iZone)%sections(sec)%elemConn(iStart:iEnd-1) = &
                            tmpConn(i+1: i+nConn)
                       
                       ! Move to the next element
                       i = i + nConn + 1
                       eCount = eCount + 1
                    case default
                       ! Not a surface or we can't read it.
                       zones(iZone)%sections(sec)%isSurface = .False.
                       exit elementLoop
                    end select
                 end do elementLoop
              end if

              ! Delete the temporary CGNS connectivity infor
              deallocate(tmpConn)

           case default
              ! Flag as not a surface
              zones(iZone)%sections(sec)%isSurface = .False.
           end select
        end do

        ! Now below we will be assuming that the elements are ordered
        ! sequentially between the muliple sections and that there are
        ! NO gaps. This is required for efficient searching. If this
        ! *isn't* the case, we will just stop. For future reference,
        ! we could make a pointer array to the sections that *is*
        ! sorted to get around this issue and make a compact pointer
        ! for the elements as well.
        do sec=1, nSections-1
           if ( zones(iZone)%sections(sec+1)%elemStart /= &
                zones(iZone)%sections(sec  )%elemEnd + 1) then 
              print *,'The section element ranges are not in order. This cannot be handled.'
              stop
           end if
        end do

        ! Determine the number of boundary conditions:
        call cg_nbocos_f(cg, base, iZone, nBocos, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Allocate space for bocos
        allocate(zones(iZone)%bocos(nBocos))

        ! On the first pass we just want to figure out the number of
        ! *surface* nodes and the size we need for connectivity. This
        ! is somehwat "dumb" since we are just going to take the
        ! number of elements * number of nodes per element, which will
        ! obviously overpredict by a large fraction. That's not a huge
        ! deal.

        bocoLoop: do boco=1, nBocos

           ! Read the info for this boundary condition. 
           call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
                ptsettype, nPnts, NormalIndex, NormalListFlag, datatype, &
                ndataset, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! We only allow ElementRange and ElementList...its just too
           ! hard to figure out what is going on with the point types.
           if (ptSetType == PointRange .or. ptSetType == pointList) then 
              print *, 'pyWarpUStruct cannot handle boundary conditions & 
                   &defined by "PointRange" or "PointList". Please use &
                   &boundary conditions defined by "ElementRange" or "ElementList" &
                   &instead'

              stop
           end if

           ! Load in the 
           if (ptSetType == ElementRange) then 
              ! We have to read off the start and end elements, this
              ! will always require an array of lenght 2 (tmpInt)

              call cg_boco_read_f(cg, base, iZone, boco, tmpInt, NULL, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f
              nBCElem = tmpInt(2) - tmpInt(1) + 1
              allocate(zones(iZone)%bocos(boco)%BCElements(nBCElem))

              ! We now need to fill it up trivially
              do i=1,nBCElem
                 zones(iZone)%bocos(boco)%BCElements(i) = tmpInt(1) + i - 1
              end do
           else if (ptSetType == ElementList) then

              ! nPnts is the actual number of elements we can allocate
              ! and read directly

              allocate(zones(iZone)%bocos(boco)%BCElements(nPnts))
              call cg_boco_read_f(cg, base, iZone, boco, &
                   zones(iZone)%bocos(boco)%BCElements, NULL, ierr)
              nBCElem = nPnts
           end if

           ! Now we need to deduce the actual connectivity for this
           ! boundary condition. This will require querying the
           ! section info. Since the same BC could reference two
           ! sections (i think in theory)

           ! On the first pass we find the sectionID and index of the
           ! element ON THAT SECTION. This require

           ! We will assume at this point that we have at most quads,
           ! so we can determine the maximum size of elemPr and elemConn
           allocate(zones(iZone)%bocos(boco)%elemConn(4*nBCElem))
           allocate(zones(iZone)%bocos(boco)%elemPtr(nBCElem + 1))
           allocate(zones(iZOne)%bocos(boco)%elemNodes(3, 4*nBCElem))

           ! Set pointers for easier reading
           elemConn => zones(iZone)%bocos(boco)%elemConn
           elemPtr  => zones(iZone)%bocos(boco)%elemPtr
           elemNodes=> zones(iZone)%bocos(boco)%elemNodes
           elemConn = -1
           elemPtr = 0
           elemPtr(1) = 1
           nElemNotFound = 0
           k = 1 ! Running index of elemConn, elemNodes
           l = 1 ! Running index of elemPtr

           do i=1, size(zones(iZone)%bocos(boco)%bcElements)
              
              curElem = zones(iZone)%bocos(boco)%bcElements(i)
              sec = -1

              ! For now do a linear search
              setionLoop: do j=1,size(zones(iZone)%sections)
                 if (curElem >= zones(iZone)%sections(j)%elemStart .and. &
                      curElem <= zones(iZone)%sections(j)%elemEnd) then
                    sec = j
                    exit setionLOop
                 end if
              end do setionLoop

              if (sec == -1)  then
                 nElemNotFound = nElemNotFound + 1
                 ! No element...do not increment the k-counter
              else
                 ! Set pointer for easier code reading
                 secPtr => zones(iZone)%sections(sec)

                 ! localIndex is the index on the section
                 localIndex = curElem - secPtr%elemStart + 1

                 ! iStart/iEnd are the indices on the section
                 iStart = secPtr%elemPtr(localIndex)
                 iEnd   = secPtr%elemPtr(localIndex+1) - 1
                 
                 ! Number of nodes on this section. 
                 curElemSize =  iEnd - iStart + 1

                 ! Since we just blew up the global connectivity, the
                 ! local connectivity is now easy...it is just
                 ! 1,2,3,4...etc. elemPtr will determine the type of
                 ! element (3 or 4 nodes). Now copy the actual nodes
                 do m=iStart,iEnd
                    elemNodes(:, k) = allNodes(:, secPtr%elemConn(m))
                    elemConn(k) = k
                    k = k + 1
                 end do

                 ! Update the elemPtr too
                 elemPtr(l+1) = elemPtr(l) + curElemSize
                 l = l + 1

                 ! Flag surface nodes
                 do m=iStart, iEnd
                    surfaceNodes(secPtr%elemConn(m)) = 1
                 end do
              end if
           end do
        
           ! Now we set the *actual* number of nBCElem and nBCNodes
           zones(iZone)%bocos(boco)%nBCElem = nBCElem - nElemNotFound 
           zones(iZOne)%bocos(boco)%nBCNodes = k-1

           ! Save the boco name
           zones(iZone)%bocos(boco)%name =  bocoName
           zones(iZone)%bocos(boco)%type = bocoType
           ! By Default no family name
           zones(iZone)%bocos(boco)%family = ""

           ! Read the family information
           call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t",&
                boco, "end")
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           if (ierr == 0) then ! Node exists
              call cg_famname_read_f(famName, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              zones(iZone)%bocos(boco)%family = famName
           end if
        end do bocoLoop

        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)
     end do zoneLoop

     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     deallocate(sizes)
  end if

  ! Communicate the total number of nodes to everyone
  call MPI_bcast(nNodes, 1, MPI_INTEGER, 0, warp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (myid == 0) then

     ! calculate the start and end indices for the root process
     istart = 0 + 1
     iend   = int(nNodes*(one/nProc))
     localSize = iend - istart + 1

     ! Allocate the local storage for the nodes
     allocate(localNodes(3, localSize))
     allocate(localSurfaceNodes(localSize))

     ! Just copy for the root proc:
     localNodes(:, :) = allNodes(:, 1:localSize)
     localSurfaceNodes(:) = surfaceNodes(1:localSize)

     ! Loop over the remainer of the procs and send the nodes they need
     do iProc=1, nProc-1
        ! compute the node range for the target processor
        istart = int(nNodes*(dble(iProc)/nProc)) +1 
        iend   = int(nNodes*(dble(iProc+1)/nProc))
        localSize = iend - istart + 1

        ! Post the send message for the nodes
#ifdef USE_COMPLEX
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_COMPLEX16, iProc, &
             11, warp_comm_world, ierr)
#else
        call MPI_Send(allNodes(:, iStart), localSize*3, MPI_REAL8, iProc, &
             11, warp_comm_world, ierr)
#endif
        call EChk(ierr, __FILE__, __LINE__)

        ! post the send message for the surface node map.
        call MPI_Send(surfaceNodes(iStart), localSize, MPI_INTEGER4, iProc, &
             12, warp_comm_world, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
     deallocate(allNodes)
     deallocate(surfaceNodes)
  else

     ! compute the incoming node range for this processor
     istart = int(nNodes*(dble(myid)/nProc)) +1 
     iend   = int(nNodes*(dble(myid+1)/nProc))
     localSize = iend - istart + 1

     ! allocate the storage for the incoming data
     allocate(localNodes(3, localSize))
     allocate(localSurfaceNodes(localSize))

     ! Post the receive on each processor
#ifdef USE_COMPLEX
     call MPI_recv(localNodes, 3*localSize, MPI_COMPLEX16, 0, 11, &
          warp_comm_world, status, ierr)
#else
     call MPI_recv(localNodes, 3*localSize, MPI_REAL8, 0, 11, &
          warp_comm_world, status, ierr)
#endif
     call EChk(ierr, __FILE__, __LINE__)
     call MPI_recv(localSurfaceNodes, localSize, MPI_INTEGER4, 0, 12, &
          warp_comm_world, status, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! All we are doing to do here is to create the Xv vector --- which
  ! is done via a call to createCommonGrid

  call createCommonGrid(localNodes, localSurfaceNodes, size(localNodes, 2))
  deallocate(localNodes, localSurfaceNodes)

end subroutine readUnstructuredCGNS


