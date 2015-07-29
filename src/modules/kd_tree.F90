Module kd_tree
  use constants
  use pointReduce
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel. Only the Euclidean metric works so far.
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,

  Integer(kind=IntType) :: bucket_size = 18
  integer(kind=intType), Parameter :: NERR = 12

  ! ..  .. Derived Type Declarations ..  

  type :: tnp
     type(tree_master_record), pointer :: tp
  end type tnp

  ! This is the master list of trees
  type(tnp), dimension(:), allocatable :: mytrees

  ! Global information about the tree pointer to the actual data array
  ! dimensionality and total # of points permuted index into the data,
  ! so that indexes[l..u] of some bucket represent the indexes of the
  ! actual points in that bucket.  root pointer of the tree an
  ! internal tree node the dimension to cut where to cut the dimension
  ! indices of points included in this node, referring back to indices
  ! child pointers 

  real(kind=realType), parameter, dimension(NERR) :: rstar = &
       (/two, three, four, five, sevenp5, ten, twenty, forty, eighty, onesixty, threetwenty, sixforty/)
  real(kind=realType),  dimension(NERR-1) :: orstar

  Type :: tree_master_record
     real(kind=realType), Dimension (:, :), allocatable :: the_data
     Integer(kind=IntType) :: dim, n
     Integer(kind=IntType), Dimension (:), Pointer :: indexes
     Type (tree_node), Pointer :: root

     ! Additional information for IDW 
     integer(kind=intType) :: maxDepth
     real(kind=realType) :: Ldef, Ldef0
     real(kind=realType) :: alphaToBExp
     real(kind=realType) , dimension(:), pointer :: Ai
     real(kind=realType) , dimension(:, :), pointer :: Bi, Bib 
     real(kind=realType) , dimension(:, :), pointer :: normals, normalsb
     real(kind=realType) , dimension(:, :), pointer :: normals0
     real(kind=realType) , dimension(:, :, :), pointer :: Mi, Mib
     real(kind=realType) , dimension(:, :), pointer :: Xu0, Xu, Xub
     integer(kind=intType), dimension(:), pointer :: XuInd
     integer(kind=intType), dimension(:), allocatable :: facePtr, faceConn
     integer(kind=intType), dimension(:, :), allocatable :: nodeToElem
     integer(kind=intType), dimension(:), allocatable :: invLink

     real(Kind=realType) :: errTol
     integer(kind=intType) :: nFace
  End Type tree_master_record

  Type :: tree_node
     Integer(kind=IntType) :: dnum
     real(kind=realType) :: val
     Integer(kind=IntType) :: l, u, n
     Type (tree_node), Pointer :: left, right

     ! Additional information for IDW
     real(kind=realType) :: Ai
     real(kind=realType) :: Bi(3), Bib(3)
     real(kind=realType) :: Mi(3, 3), Mib(3, 3)
     real(kind=realType) :: X(3)
     real(kind=realType) :: oradius
     integer(kind=intType) :: lvl
     real(kind=realType) :: err(NERR)
  End Type tree_node

  Type :: tree_search_record
     Private
     real(kind=realType) :: bsd
     Integer(kind=IntType) :: bestind
     real(kind=realType), Dimension (:), Pointer :: qv
     Integer(kind=IntType), Dimension (:), Pointer :: il
     real(kind=realType), Dimension (:), Pointer :: dsl
     Integer(kind=IntType) :: nbst, nfound
     Integer(kind=IntType) :: centeridx, correltime
  End Type tree_search_record
Contains

  Subroutine destroy_tree(tp)
    ! Deallocates all memory for the tree
    implicit none
    Type (tree_master_record), Pointer :: tp

    call destroy_node(tp%root)
    deallocate(tp%indexes, tp%the_data)
    deallocate(tp%Ai, tp%Bi, tp%Mi)
    deallocate(tp%normals, tp%normals0)
    deallocate(tp%Xu0, tp%Xu, tp%XuInd)
    deallocate(tp%facePtr, tp%faceConn, tp%nodeToElem)
    deallocate(tp%invLink)
  Contains

    Recursive Subroutine destroy_node(np)
      implicit none
      Type (tree_node), Pointer :: np
      If (ASSOCIATED(np%left)) Then
         Call destroy_node(np%left)
         Deallocate (np%left)
      End If
      If (ASSOCIATED(np%right)) Then
         Call destroy_node(np%right)
         Deallocate (np%right)
      End If
      Return
    End Subroutine destroy_node
  End Subroutine destroy_tree

  subroutine build_tree(tp)
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Local Scalars ..
    Integer(kind=IntType) :: j
    ! ..
    Allocate (tp%indexes(tp%n))
    Do j = 1, tp%n
       tp%indexes(j) = j
    End Do
    tp%root => build_tree_for_range(tp, 1, tp%n)

  contains

    Recursive Function build_tree_for_range(tp, l, u) Result (res)
      ! .. Function Return Value ..
      Type (tree_node), Pointer :: res      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      Integer(kind=IntType) :: c, m
      ! ..
      if (.false.) then 
         if ((l .lt. 1) .or. (l .gt. tp%n)) then
            stop 'illegal L value in build_tree_for_range'
         end if
         if ((u .lt. 1) .or. (u .gt. tp%n)) then
            stop 'illegal u value in build_tree_for_range'
         end if
         if (u .lt. l) then
            stop 'U is less than L, thats illegal.'
         end if
      endif
      If ((u-l)<=bucket_size) Then
         Allocate (res)
         res%dnum = 0
         res%val = 0.0
         res%l = l
         res%u = u
         res%n = u-l+1
         Nullify (res%left, res%right)
      Else
         Allocate (res)
         c = most_spread_coordinate(tp, l, u)
         m = (l+u)/2
         Call select_on_coordinate(tp%the_data, tp%indexes, c, m, l, u)
         ! moves indexes around
         res%dnum = c
         res%val = tp%the_data(c, tp%indexes(m))
         res%l = l
         res%u = u
         res%n = u-l+1
         res%left => build_tree_for_range(tp, l, m)
         res%right => build_tree_for_range(tp, m+1, u)
      End If
    End Function build_tree_for_range

    Subroutine select_on_coordinate(v, ind, c, k, li, ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Structure Arguments ..
      !      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: c, k, li, ui
      ! ..
      ! .. Local Scalars ..
      Integer(kind=IntType) :: i, l, m, s, t, u
      ! ..
      ! .. Local Arrays ..
      real(kind=realType) :: v(:, :)
      Integer(kind=IntType) :: ind(:)
      ! ..
      !v => tp%the_data
      !ind => tp%indexes
      l = li
      u = ui
      Do While (l<u)
         t = ind(l)
         m = l
         Do i = l + 1, u
            If (v(c, ind(i))<v(c, t)) Then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            End If
         End Do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         If (m<=k) l = m + 1
         If (m>=k) u = m - 1
      End Do
    End Subroutine select_on_coordinate

    Function most_spread_coordinate(tp, l, u) Result (res)
      ! Of indices in l..u find the axis which has the largest spread, 
      ! and
      ! return its index
      ! .. Function Return Value ..
      Integer(kind=IntType) :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      real(kind=realType) :: bsp, sp
      Integer(kind=IntType) :: i = 0
      ! ..
      res = 0
      bsp = -1.0
      Do i = 1, tp%dim
         sp = spread_in_coordinate(tp, i, l, u)
         If (sp>bsp) Then
            res = i
            bsp = sp
         End If
      End Do
    End Function most_spread_coordinate

    Function spread_in_coordinate(tp, c, l, u) Result (res)
      ! the spread in coordinate 'c', between l and u
      ! for easier local access
      ! ibid
      ! .. Function Return Value ..
      real(kind=realType) :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(kind=realType) :: last, lmax, lmin, smax, smin, t
      Integer(kind=IntType) :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(kind=realType), Pointer :: v(:, :)
      Integer(kind=IntType), Pointer :: ind(:)
      ! ..
      v => tp%the_data
      ind => tp%indexes
      smin = v(c, ind(l))
      smax = smin

      ! truncate u here or not?????
      ulocal = min(u, l+200) !!!?????

      Do i = l + 2, ulocal, 2
         lmin = v(c, ind(i-1))
         lmax = v(c, ind(i))
         If (lmin>lmax) Then
            t = lmin
            lmin = lmax
            lmax = t
         End If
         If (smin>lmin) smin = lmin
         If (smax<lmax) smax = lmax
      End Do
      If (i==ulocal+1) Then
         last = v(c, ind(ulocal))
         If (smin>last) smin = last
         If (smax<last) smax = last
      End If
      res = smax - smin
    End Function spread_in_coordinate
  end subroutine build_tree
  function pt_in_tree(tp, qv)
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    real(kind=realType), Target, Intent (In) :: qv(:)
    real(kind=realType), parameter :: tol=1e-6

    ! .. Scalar Arguments ..
    Integer(kind=IntType) :: n
    real(kind=realType), target :: distances(1)
    Integer(kind=IntType), target :: indexes(1)
    integer(kind=intType) :: pt_in_tree
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic HUGE
    n = 1
    ! ..
    ! the largest real number
    sr%bsd = HUGE(1.0)
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(1.0)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr

    Call n_search(tp, psr, tp%root)

    if (sqrt(distances(1)) < tol) then
       pt_in_tree = indexes(1)
    else
       pt_in_tree = 0_intType
    end if

  end function pt_in_tree

  Subroutine n_nearest_to(tp, qv, n, indexes, distances)
    ! find the 'n' nearest neighbors to 'qv', returning
    ! their indexes in indexes and squared Euclidean distances in
    ! 'distances'
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer(kind=IntType), Intent (In) :: n
    ! ..
    ! .. Array Arguments ..
    real(kind=realType), Target :: distances(n)
    real(kind=realType), Target, Intent (In) :: qv(:)
    Integer(kind=IntType), Target :: indexes(n)
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic HUGE
    ! ..
    ! the largest real number
    sr%bsd = HUGE(1.0)
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(1.0)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr

    Call n_search(tp, psr, tp%root)
  End Subroutine n_nearest_to

  Recursive Subroutine n_search(tp, sr, node)
    ! This is the innermost core routine of the kd-tree search.
    ! it is thus not programmed in quite the clearest style, but is
    ! designed for speed.  -mbk
    ! .. Structure Arguments ..
    Type (tree_node), Pointer :: node
    Type (tree_search_record), Pointer :: sr
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Local Scalars ..
    real(kind=realType) :: dp, sd, sdp, tmp
    Integer(kind=IntType) :: centeridx, i, ii, j, jmax, k, d, correltime, nbst
    Logical :: not_fully_sized
    ! ..
    ! .. Local Arrays ..
    real(kind=realType), Pointer :: qv(:)
    Integer(kind=IntType), Pointer :: ind(:)
    real(kind=realType), dimension(:, :), pointer :: data
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic ABS, ASSOCIATED
    ! ..
    If ( .Not. (ASSOCIATED(node%left)) .And. ( .Not. ASSOCIATED( &
         node%right))) Then
       ! we are on a terminal node
       ind => tp%indexes     ! save for easy access
       qv => sr%qv
       data => tp%the_data
       centeridx = sr%centeridx
       d = tp%dim
       correltime = sr%correltime
       ! search through terminal bucket.
       nbst = sr%nbst
       mainloop: Do i = node%l, node%u
          ii = ind(i)
          If ( (centeridx<0) .OR. (ABS(ii-centeridx)>=correltime)) Then
             ! 
             ! replace with call to square distance with inline
             ! code, an
             ! a goto.  Tested to be significantly faster.   SPECIFIC
             ! FOR
             ! the EUCLIDEAN METRIC ONLY! BEWARE!

             sd = 0.0
             Do k = 1, d
                tmp = data(k, ii) - qv(k)
                sd = sd + tmp*tmp
                If (sd>sr%bsd) Then
                   cycle mainloop
                End If
             End Do

             ! Note test moved out of loop to improve vectorization
             ! should be semantically identical eitehr way as sr%bsd is
             ! a bound afterwhich it doesn't matter. 

             ! we only consider it if it is better than the 'best' on
             ! the list so far.
             ! if we get here
             ! we know sd is < bsd, and bsd is on the end of the list

             ! special case for nbst = 1 (single nearest neighbor)
             if (nbst .eq. 1) then
                sr%il(1) = ii
                sr%dsl(1) = sd
                sr%bsd = sd
             else
                not_fully_sized = (sr%nfound<nbst)
                If (not_fully_sized) Then
                   jmax = sr%nfound
                   sr%nfound = sr%nfound + 1
                Else
                   jmax = nbst - 1
                End If
                ! add it to the list

                ! find the location j where sd >= sr%dsl(j) and sd <
                ! sr%dsl(j+1...jmax+1)
                Do j = jmax, 1, -1
                   If (sd>=sr%dsl(j)) Exit ! we hit insertion location
                   sr%il(j+1) = sr%il(j)
                   sr%dsl(j+1) = sr%dsl(j)
                End Do
                ! if loop falls through j=0 here.
                sr%il(j+1) = ii
                sr%dsl(j+1) = sd

                If ( .Not. not_fully_sized) Then
                   sr%bsd = sr%dsl(nbst)
                End If
             end if
          End If
       End Do mainloop
    Else
       ! we are not on a terminal node

       ! Alrighty, this section is essentially the content of the
       ! the Sproul method for searching a Kd-tree, in other words
       ! the second nested "if" statements in the two halves below
       ! and when they get activated.
       dp = sr%qv(node%dnum) - node%val
       sdp = dp*dp        ! Euclidean
       If (dp<0.0) Then
          Call n_search(tp, sr, node%left)
          If (sdp<sr%bsd) Call n_search(tp, sr, node%right)
          ! if the distance projected to the 'wall boundary' is less
          ! than the radius of the ball we must consider, then perform
          ! search on the 'other side' as well.
       Else
          Call n_search(tp, sr, node%right)
          If (sdp<sr%bsd) Call n_search(tp, sr, node%left)
       End If
    End If
  End Subroutine n_search

  ! --------------------------------------------------------------------
  !                 Additional routines for IDW Warping
  ! --------------------------------------------------------------------

  subroutine labelNodes(tp)
    ! Recurse through the tree and give a unique identifier to each
    ! node and identify which level each node is. Additionally, we
    ! store tp%maxDepth which is the maximum number of layers in the
    ! tree.
    implicit none
    ! Subroutine arguments
    Type (tree_master_record), Pointer :: tp
    integer(kind=intType) :: level
    integer(kind=intType) :: id
    level = 1
    tp%maxDepth = 1
    id = 1
    call labelNodes_node(tp, tp%root, level, id)
  contains
    recursive subroutine labelNodes_node(tp, np, level, id)
      ! Give each node a unique identifier as well as what level in
      ! the tree it is.
      use constants
      implicit none
      Type (tree_master_record), Pointer :: tp
      Type (tree_node), Pointer :: np 
      integer(kind=intType) :: level
      integer(kind=intType) :: id
      np%lvl = level
      id = id + 1
      tp%maxDepth = max(tp%maxDepth, level)
      if (np%dnum /= 0) then 
         call labelNodes_node(tp, np%left, level+1, id)
         call labelNodes_node(tp, np%right, level+1, id)
      end if
    end subroutine labelNodes_node
  end subroutine labelNodes

  recursive subroutine setAi(tp, np)
    ! Propagate the node weights to all nodes. 
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    integer(kind=intType) :: i
    if (np%dnum == 0) then 
       ! On the leaf node...determine the summed Ai
       np%Ai = zero
       do i=np%l, np%u
          np%Ai = np%Ai + tp%Ai(i)
       end do
    else
       ! Call each of the children
       call setAi(tp, np%left)
       call setAi(tp, np%right)

       ! Now add the sums from the two children 
       np%Ai = (np%left%Ai + np%right%Ai)
    end if
  end subroutine setAi

  recursive subroutine setCenters(tp, np)
    ! Compute the center of each node. This is just the average X
    ! value. We also compute the minimum radius of the
    ! centroid-centered sphere that will contain all nodes. This
    ! routine computes np%X and np%radius

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    integer(kind=intType) :: i
    real(kind=realType) :: Xcen(3), Xmax, dist, dx(3), radius
    
    ! Sum up all the children
    np%X = zero
    do i=np%l, np%u
       np%X = np%X + tp%the_data(:, tp%indexes(i))
    end do
    np%X = np%X / np%n

    ! Now do the radius
    radius = zero
    do i=np%l, np%u
       dx = tp%the_data(:, tp%indexes(i)) - np%X
       dist = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
       if (dist > radius) then 
          radius = dist
       end if
       np%oradius = one/radius
    end do

    ! Now call the children
    if (np%dnum /= 0) then 
       call setCenters(tp, np%left)
       call setCenters(tp, np%right)
    end if
  end subroutine setCenters

  Recursive Subroutine setData(tp, np)
    ! Update the nodal displacements and rotation matrices based on a
    ! new set of displcements.

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    integer(kind=intType) :: i
    if (np%dnum /= 0) then 
       np%Bi = zero
       np%Mi = zero
       do i=np%l, np%u
          np%Bi = np%Bi + tp%Bi(:, i)
          np%Mi = np%Mi + tp%Mi(:, :, i)
       end do
       np%Bi = np%Bi / np%n
       np%Mi = np%Mi / np%n

       ! Call each of the children
       call setData(tp, np%left)
       call setData(tp, np%right)
    end if
  end subroutine setData

  recursive subroutine evalNode(tp, np, r, num, den, approxDen)
    ! This the main fast evaluation routine for the displacements. It
    ! is therefore recrusive. Note that num and den must be zeroed
    ! externally before this call. approxDen is used for checking
    ! error tolerances. 

    implicit none
    ! Subroutine arguments
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout), dimension(3) :: num
    real(kind=realType), intent(inout) :: den
    real(kind=realType), intent(in) :: approxDen

    ! Working variables
    real(kind=realType), dimension(3) :: Si
    real(kind=realType) :: dist, err

    if (np%dnum == 0) then 
       ! Leaf node -> Must do exact calc:
       call EvalNodeExact(tp, np, r, num, den)
    else
       ! Tree Node: Compute the distance from 'r' to the
       ! center of the node, np%X
       dist = sqrt((r(1) - np%X(1))**2 + &
                   (r(2) - np%X(2))**2 + & 
                   (r(3) - np%X(3))**2)

       if (dist*np%oradius < two) then ! Too close...call children
          call evalNode(tp, np%left, r, num, den, approxDen)
          call evalNode(tp, np%right, r, num, den, approxDen)
       else 
          ! Use the first error check:
          call getError(tp, np, dist, err)
          if (err < tp%errTol * approxDen) then
             call evalNodeApprox(tp, np, r, num, den, dist)
          else
             ! Not good enough error so call children
             call evalNode(tp, np%left, r, num, den, approxDen)
             call evalNode(tp, np%right, r, num, den, approxDen)
          end if
       end if
    end if
  end subroutine evalNode

  subroutine evalNodeExact(tp, np, r, num, den)
    ! Loop over the owned nodes and evaluate the exact numerator and
    ! denomenator at the requesed r. Note that this rouine *may* be
    ! called on *any* node including the root node. When called on the
    ! root node this is actually the exact scheme. 
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout), dimension(3) :: num
    real(kind=realType), intent(inout) :: den
    real(kind=realType), dimension(3) :: rr, Si
    real(kind=realType) :: LdefoDist, LdefoDist3
    real(kind=realType) :: Wi
    integer(kind=intType) :: i

    do i=np%l, np%u
       rr(1) = r(1) - tp%Xu0(1, i)
       rr(2) = r(2) - tp%Xu0(2, i)
       rr(3) = r(3) - tp%Xu0(3, i)
       LdefoDist = tp%Ldef/sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+eps)
       Ldefodist3 = LdefoDist**3
       Wi = tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       Si = tp%Mi(:, 1, i)*r(1) + tp%Mi(:, 2, i)*r(2) + tp%Mi(:, 3, i)*r(3) + tp%bi(:, i) - r(:)
       num = num + Wi*Si
       den = den + Wi
    end do
  end subroutine evalNodeExact

  subroutine evalNodeApprox(tp, np, r, num, den, dist)
    ! We are far enough away! Whoo! Just use the data on this node
    ! (np). Note that dist must be passed in. It is assumed that the
    ! dist is already available from a calculation to determine
    ! whether to use this approximation or not.
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in) :: dist, r(3)
    real(kind=realType), intent(out), dimension(3) :: num
    real(kind=realType), intent(out) :: den
    real(kind=realType) :: LdefoDist, LdefoDist3
    real(kind=realType) :: Wi, Si(3)

    LdefoDist = tp%Ldef/dist
    Ldefodist3 = LdefoDist**3
    Wi = np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
    Si(1) = np%Mi(1, 1)*r(1) + np%Mi(1, 2)*r(2) + np%Mi(1, 3)*r(3) + np%bi(1) - r(1)
    Si(2) = np%Mi(2, 1)*r(1) + np%Mi(2, 2)*r(2) + np%Mi(2, 3)*r(3) + np%bi(2) - r(2)
    Si(3) = np%Mi(3, 1)*r(1) + np%Mi(3, 2)*r(2) + np%Mi(3, 3)*r(3) + np%bi(3) - r(3)
    num = num + Wi*Si
    den = den + Wi
  end subroutine evalNodeApprox

  subroutine writeTreeTecplot(tp, fileName)
    ! This is a debuging routine that writes the centers of all of the
    ! tree nodes to a tecplot file. **IT DOES NOT WRITE THE NODES IN
    ! THE LEAVES**.
    implicit none
    Type (tree_master_record), Pointer :: tp    
    character*(*), intent(in) :: fileName
    integer(kind=intType) :: lvl

    open(unit=7, file=fileName, status='replace')
101 format('ZONE T="Level ',I2, '"')
    do lvl=1, tp%maxDepth-1
       write(7, 101), lvl
       call writeLevel(tp, tp%root, lvl)
    end do
    close(7)
  contains
    recursive subroutine writeLevel(tp, np, lvlToWrite)
      use constants
      implicit none
      Type (tree_master_record), Pointer :: tp    
      Type (tree_node), Pointer :: np
      integer(kind=intType) :: lvlToWrite
102   format(g20.12, g20.12, g20.12, g20.12)

      if (np%lvl == lvlToWrite) then
         write(7, 102) np%X(1), np%X(2), np%X(3)
      end if

      if (np%dnum /= 0) then
         call writeLevel(tp, np%left, lvlToWrite)
         call writeLevel(tp, np%right, lvlToWrite)
      end if
    end subroutine writeLevel
  end subroutine writeTreeTecplot

  recursive subroutine getWiEstimate(tp, np, r, den)
    ! This routine needs to be only called once. Essentialy what we are
    ! dong is computing just the denomenator (Wi) computation. This is
    ! used as a reference value for the 'real' mesh warp, to be used for
    ! checking against the errors. A hard coded R-limit of 5.0 is used
    ! here. This should be sufficient to get an good estimate of Wi.
    implicit none
    ! Subroutine arguments
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout) :: den
    real(kind=realType), dimension(3) :: rr
    real(kind=realType) :: dist, LdefoDist, LdefoDist3, err
    integer(kind=intType) :: i
    if (np%dnum == 0) then 
       ! Leaf node. Do regular calc:
       do i=np%l, np%u
          rr = r - tp%Xu0(:, i)
          LdefoDist = tp%Ldef/sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)
          Ldefodist3 = LdefoDist**3
          den = den + tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       end do
    else
       rr = r - np%X
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
       
       if (dist*np%oradius > 5.0) then
          ! Far field calc is ok:
          LdefoDist = tp%Ldef/dist
          Ldefodist3 = LdefoDist**3
          den = den + np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       else
          ! To far away, recursively call the children
          call getWiEstimate(tp, np%left, r, den)
          call getWiEstimate(tp, np%right, r, den)
       end if
    end if
  end subroutine getWiEstimate

  recursive subroutine dryRun(tp, np, r, approxDen, c)
    ! This routine estimates the number of operations that *would* be
    ! needed to compute the displacement for a given node. No
    ! computations are actually done. It is used to determine the
    ! approximate costs so that a load balancing algorithm can used
    ! to rebalance the nodes. 

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(in) :: approxDen
    real(kind=realType), intent(out) :: c
    real(kind=realType), dimension(3) :: rr
    real(kind=realType) :: dist, err

    if (np%dnum == 0) then 
       ! The regular loop gets assigned a cost of 1.0.
       c = c + dble(np%n)
    else
       rr = r - np%X
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
       ! This 0.5 cost accounts for the subroutine calling
       ! overhead and the above dist calc.
       c = c + half
       if (dist*np%oradius < two) then ! Too close...call children
          call dryRun(tp, np%left, r, approxDen, c)
          call dryRun(tp, np%right, r, approxDen, c)
       else 
          ! Use the first error check:
          call getError(tp, np, dist, err)
          if (err < tp%errTol * approxDen) then
             ! This is includes the getError calc
             c = c + 1.1
          else
             ! Not good enough error so call children
             call dryRun(tp, np%left, r, approxDen, c)
             call dryRun(tp, np%right, r, approxDen, c)
          end if
       end if
    end if
  end subroutine dryRun

  recursive subroutine computeErrors(tp, np)
    ! This routine computes an estimate of the error to be induced in
    ! the denomenator at each tree node. Then, when we are evaluating
    ! the displacment we can determine if the errors induced by the
    ! approximation are acceptable or not.

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType) , dimension(3, 20) :: vpts
    real(kind=realType) , dimension(3, 12) :: fpts
    real(kind=realType) :: dExact(20), dApprox(20), rr(3)
    real(kind=realType) :: LdefoDist, LdefODist3, dist
    integer(kind=intType) :: ii, i, j

    if (np%dnum /= 0) then 
       do j=1, NERR
          dExact= zero
          dApprox = zero
          call getSpherePts(np%X, one/np%oradius*rstar(j), vpts, fpts)

          np%err(j) = zero
          do ii=1, 20
             ! Compute exact value:
             do i=np%l, np%u
                rr = vpts(:, ii) - tp%Xu0(:, i)
                dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
                LdefoDist = tp%Ldef/dist
                Ldefodist3 = LdefoDist**3
                dExact(ii) = dExact(ii) + tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
             end do

             ! And the approx value:
             rr = vpts(:, ii) - np%X
             dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
             LdefoDist = tp%Ldef/dist
             Ldefodist3 = LdefoDist**3
             dApprox(ii) = dApprox(ii) + np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)

             ! Now set the nodal error:
             np%err(j) = np%err(j) + (dExact(ii) - dApprox(ii))**2 
          end do

          ! This the RMS Error:
          np%err(j) = sqrt(np%err(j)/20)
       end do
       ! Now call the children:
       call computeErrors(tp, np%left)
       call computeErrors(tp, np%right)
    end if
  end subroutine computeErrors

  subroutine getError(tp, np, dist, err)
    ! Simple routine to use the stored errors to estimate what the
    ! error will be for dist 'dist'.
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in) :: dist
    real(Kind=realType), intent(out) :: err
    real(kind=realType) :: dOvrR
    real(kind=realType) :: fact
    integer(kind=intType) :: i

    dOvrR = dist * np%oradius
    
    ! This is a linear search. However, trust me, it is faster than a
    ! binary search since N is small and fixed
    err = np%err(NERR)
    do i=1,11
       if (dOvrR < rstar(i+1)) then
          fact = (dOvrR - rstar(i))*orstar(i)
          err = (one-fact)*np%err(i) + fact*np%err(i+1)
          exit
       end if
    end do
 
  end subroutine getError

  !------------------------------------------------------------
  !               DERIVATIVE ROUTINES
  !------------------------------------------------------------
  subroutine allocDerivValues(tp)
    ! Allocate the required derivative values
    implicit none
    Type (tree_master_record), Pointer :: tp
    allocate(tp%Bib(3, tp%n), tp%Mib(3, 3, tp%n), tp%normalsb(3, tp%n))
    allocate(tp%Xub(3, tp%n))
    tp%Xub = zero
    tp%Bib = zero
    tp%Mib = zero
    tp%normalsb = zero
  end subroutine allocDerivValues

  subroutine deallocDerivValues(tp)
    ! Allocate the required derivative values
    implicit none
    Type (tree_master_record), Pointer :: tp
    deallocate(tp%Bib, tp%Mib, tp%normalsb, tp%Xub)
  end subroutine deallocDerivValues

  recursive subroutine zeroDeriv(np)
    ! Zero the derivatives in the nodes
    implicit none
    Type (tree_node), Pointer :: np 
    np%Bib = zero
    np%Mib = zero
    if (np%dnum /= 0) then
       call zeroDeriv(np%left)
       call zeroDeriv(np%right)
    end if
  end subroutine zeroDeriv
#ifndef USE_COMPLEX
  recursive subroutine evalNode_b(tp, np, r, numb, approxDen, Bib, Mib)
    implicit none
    ! Subroutine arguments
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), dimension(3) :: numb
    real(kind=realType), intent(in) :: approxDen
    real(kind=realType), dimension(:, :) :: Bib
    real(kind=realType), dimension(:, :, :) :: Mib
    ! Working variables
    real(kind=realType), dimension(3) :: rr
    real(kind=realType) :: dist, err

    if (np%dnum == 0) then 
       ! Leaf node -> Must do exact calc:
       call evalNodeExact_b(tp, np, r, numb, Bib, Mib)
    else
       ! Tree Node: Compute the distance from 'r' to the
       ! center of the node, np%X
       rr(1) = r(1) - np%X(1)
       rr(2) = r(2) - np%X(2)
       rr(3) = r(3) - np%X(3)
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)

       if (dist*np%oradius < two) then ! Too close...call children
          call evalNode_b(tp, np%left, r, numb, approxDen, Bib, Mib)
          call evalNode_b(tp, np%right, r, numb, approxDen, Bib, Mib)
       else 
          ! Use the first error check:
          call getError(tp, np, dist, err)
          if (err < tp%errTol * approxDen) then
             call evalNodeApprox_b(tp, np, numb, r, dist)
          else
             ! Not good enough error so call children
             call evalNode_b(tp, np%left, r, numb, approxDen, Bib, Mib)
             call evalNode_b(tp, np%right, r, numb, approxDen, Bib, Mib)
          end if
       end if
    end if
  end subroutine evalNode_b

  subroutine evalNodeExact_b(tp, np, r, numb, Bib, Mib)
    ! Reverse derviative of evalNodeExact

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    integer(kind=intType) :: i
    real(kind=realType), intent(in), dimension(3) :: r, numb
    real(kind=realType), dimension(3) :: rr
    real(kind=realType), dimension(:, :) :: Bib
    real(kind=realType), dimension(:, :, :) :: Mib
    real(kind=realType) :: dist, LdefoDist, LdefoDist3, sib(3)
    real(kind=realType) :: Wi

    do i=np%l,np%u
       rr(1) = r(1) - tp%xu0(1, i)
       rr(2) = r(2) - tp%xu0(2, i)
       rr(3) = r(3) - tp%xu0(3, i)

       dist = SQRT(rr(1)**2 + rr(2)**2 + rr(3)**2 + eps)
       ldefodist = tp%ldef/dist
       ldefodist3 = ldefodist**3
       wi = tp%Ai(i)*(ldefodist3+tp%alphatobexp*ldefodist3*ldefodist*ldefodist)
       Sib = wi*numb
       mib(:, 1, i) = mib(:, 1, i) + r(1)*sib
       mib(:, 2, i) = mib(:, 2, i) + r(2)*sib
       mib(:, 3, i) = mib(:, 3, i) + r(3)*sib
       bib(:, i) = bib(:, i) + sib
    end do

  end subroutine evalNodeExact_b

  subroutine evalNodeApprox_b(tp, np, numb, r, dist)
    ! Reverse Sensitivity of evalNodeApprox
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in) :: dist, r(3)
    real(kind=realType), dimension(3) :: numb
    real(kind=realType) :: LdefoDist, LdefoDist3, Wi, sib(3)

    ldefodist = tp%ldef/dist
    ldefodist3 = ldefodist**3
    wi = np%ai*(ldefodist3+tp%alphatobexp*ldefodist3*ldefodist*ldefodist)
    Sib = wi*numb
    np%mib(:, 1) = np%mib(:, 1) + r(1)*sib
    np%mib(:, 2) = np%mib(:, 2) + r(2)*sib
    np%mib(:, 3) = np%mib(:, 3) + r(3)*sib
    np%bib = np%bib + sib
  end subroutine evalNodeApprox_b

  recursive subroutine setData_b(tp, np)
    ! This performs the reverse accumulation of the nodeal Bib (and
    ! Mib) into the full Bib array.
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType) :: ovrN
    integer(kind=intType) :: i 
    if (np%dnum /= 0) then 
       ovrN = one/ np%n
       do i=np%l, np%u
          tp%Bib(:, i)    = tp%Bib(:, i)    + np%Bib * ovrN
          tp%Mib(:, :, i) = tp%Mib(:, :, i) + np%Mib * ovrN
       end do
       ! Call each of the children
       call setData_b(tp, np%left)
       call setData_b(tp, np%right)
    end if
  end subroutine setData_b

  ! This originally came from tapenade, but has to extensively
  ! manually modified for it work iwthout having an enitrely new
  ! kd_tree_b. 
#ifndef USE_TAPENADE
  SUBROUTINE COMPUTENODALPROPERTIES_B(tp, initialpoint)
    USE GRIDINPUT
    USE COMMUNICATION
    USE GRIDDATA
    IMPLICIT NONE
    ! Subroutine Variables
    TYPE(TREE_MASTER_RECORD), POINTER :: tp
    LOGICAL :: initialpoint
    ! Working variables
    REAL(kind=realtype), DIMENSION(3, 20) :: points
    REAL(kind=realtype), DIMENSION(3, 20) :: pointsb
    INTEGER(kind=inttype) :: i, j, jj, kk, npts, nelem, ind
    INTEGER(kind=inttype) :: nptsmax
    REAL(kind=realtype) :: facearea, facenormal(3)
    REAL(kind=realtype) :: faceareab, facenormalb(3)
    REAL(kind=realtype) :: sumarea, sumnormal(3), si(3), ds(3), smean(3)&
         &    , da, eta, r(3), dx(3)
    REAL(kind=realtype) :: sumareab, sumnormalb(3), dab
    INTEGER :: ad_to, branch
    ! This performs the copy and mirroring as required
    DO i=1,tp%n
       j = tp%xuind(i)
       tp%xu(:, i) = xsptr(3*j-2:3*j)
    END DO
    ! Now we loop over nodes and fill up Ai, Bi, and Mi
    tp%normals = zero
    tp%normalsb = 0.0_8
    tp%xub = 0.0_8
    pointsb = 0.0_8
    DO i=1,tp%n
       sumarea = zero
       sumnormal = zero
       nelem = tp%nodetoelem(1, i)

       DO jj=1,nelem
          ind = tp%nodetoelem(1+jj, i)
          ! Extract points for this face
          npts = tp%faceptr(ind) - tp%faceptr(ind-1)
          DO kk=1,npts
             CALL PUSHREAL8ARRAY(points(:, kk), realtype*3/8)
             points(:, kk) = tp%xu(:, tp%faceconn(tp%faceptr(ind-1)+kk))
          END DO
          CALL PUSHINTEGER4(kk-1)
          CALL PUSHREAL8ARRAY(facenormal, realtype*3/8)
          CALL GETELEMENTPROPS(points, npts, facearea, facenormal)
          CALL PUSHREAL8ARRAY(da, realtype/8)
          ! For face 'ind' how many nodes are on the element?
          da = facearea/(tp%faceptr(ind)-tp%faceptr(ind-1))
          sumarea = sumarea + da
          sumnormal = sumnormal + da*facenormal
       END DO
       tp%normals(:, i) = sumnormal/sumarea
       IF (.NOT.initialpoint) THEN
          ! Now get the rotation Matrix
          ! Now get the rotation Matrix
          IF (userotations) THEN
             CALL PUSHCONTROL1B(0)
          ELSE
             CALL PUSHCONTROL1B(1)
          END IF
          tp%xub(:, i) = tp%xub(:, i) + tp%bib(:, i)
          tp%mib(:, 1, i) = tp%mib(:, 1, i) - tp%xu0(1, i)*tp%bib(:, i)
          tp%mib(:, 2, i) = tp%mib(:, 2, i) - tp%xu0(2, i)*tp%bib(:, i)
          tp%mib(:, 3, i) = tp%mib(:, 3, i) - tp%xu0(3, i)*tp%bib(:, i)
          tp%bib(:, i) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
             CALL GETROTATIONMATRIX3D_B(tp%normals0(:, i), tp%normals(:, i)&
                  , tp%normalsb(:, i), tp%mi(:, :, i), tp%mib(:, :, i))
             tp%mib(:, :, i) = 0.0_8
          ELSE
             tp%mib(:, :, i) = 0.0_8
          END IF
       END IF

       sumnormalb = 0.0_8
       sumnormalb = tp%normalsb(:, i)/sumarea
       sumareab = SUM(-(sumnormal*tp%normalsb(:, i)/sumarea))/sumarea
       tp%normalsb(:, i) = 0.0_8
       DO jj=nelem,1,-1
          facenormalb = 0.0_8
          dab = sumareab + SUM(facenormal*sumnormalb)
          facenormalb = da*sumnormalb
          ind = tp%nodetoelem(1+jj, i)
          CALL POPREAL8ARRAY(da, realtype/8)
          faceareab = dab/(tp%faceptr(ind)-tp%faceptr(ind-1))
          npts = tp%faceptr(ind) - tp%faceptr(ind-1)
          CALL POPREAL8ARRAY(facenormal, realtype*3/8)
          CALL GETELEMENTPROPS_B(points, pointsb, npts, facearea, &
               &                         faceareab, facenormal, facenormalb)
          CALL POPINTEGER4(ad_to)
          DO kk=ad_to,1,-1
             CALL POPREAL8ARRAY(points(:, kk), realtype*3/8)
             tp%xub(:, tp%faceconn(tp%faceptr(ind-1)+kk)) = tp%xub(:, tp%&
                  &            faceconn(tp%faceptr(ind-1)+kk)) + pointsb(:, kk)
             pointsb(:, kk) = 0.0_8
          END DO
       END DO
    END DO
    !xsptrb = 0.0_8 ! CAN'T ZERO HERE! If we have multiple outside trees
    DO i=1, tp%n
       j = tp%xuind(i)
       xsptrb(3*j-2:3*j) = xsptrb(3*j-2:3*j) + tp%xub(:, i)
       !tp%xub(:, i) = 0.0_8
    END DO
    ! tp%bib = 0.0_8
    ! tp%mib = 0.0_8
  END SUBROUTINE COMPUTENODALPROPERTIES_B

  recursive subroutine evalNode_d(tp, np, r, numd, approxDen)
    implicit none
    ! Subroutine arguments
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np 
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout), dimension(3) :: numd
    real(kind=realType), intent(in) :: approxDen

    ! Working variables
    real(kind=realType), dimension(3) :: rr
    real(kind=realType) :: dist, err

    if (np%dnum == 0) then 
       ! Leaf node -> Must do exact calc:
       call evalNodeExact_d(tp, np, r, numd)
    else
       ! Tree Node: Compute the distance from 'r' to the
       ! center of the node, np%X
       rr(1) = r(1) - np%X(1)
       rr(2) = r(2) - np%X(2)
       rr(3) = r(3) - np%X(3)
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)

       if (dist*np%oradius < two) then ! Too close...call children
          call evalNode_d(tp, np%left, r, numd, approxDen)
          call evalNode_d(tp, np%right, r, numd, approxDen)
       else 
          ! Use the first error check:
          call getError(tp, np, dist, err)
          if (err < tp%errTol * approxDen) then
             call evalNodeApprox_d(tp, np, r, numd, dist)
          else
             ! Not good enough error so call children
             call evalNode_d(tp, np%left, r, numd, approxDen)
             call evalNode_d(tp, np%right, r, numd, approxDen)
          end if
       end if
    end if
  end subroutine evalNode_d

  subroutine evalNodeExact_d(tp, np, r, numd)
    ! Forward derviative of evalNodeExact
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout), dimension(3) :: numd
    real(kind=realType), dimension(3) :: rr, Sib
    real(kind=realType) :: dist, LdefoDist, LdefoDist3
    real(kind=realType) :: Wi
    integer(kind=intType) :: i

    do i=np%l, np%u
       rr(1) = r(1) - tp%Xu0(1, i)
       rr(2) = r(2) - tp%Xu0(2, i)
       rr(3) = r(3) - tp%Xu0(3, i)
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+eps)
       LdefoDist = tp%Ldef/dist
       Ldefodist3 = LdefoDist**3
       Wi = tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       Sib = tp%Mib(:, 1, i)*r(1) + tp%Mib(:, 2, i)*r(2) + tp%Mib(:, 3, i)*r(3) + tp%bib(:, i) 
       numd = numd + Wi*Sib
    end do

  end subroutine evalNodeExact_d

  subroutine evalNodeApprox_d(tp, np, r, numd, dist)
    ! Forward Sensitivity of evalNodeApprox
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in) :: dist
    real(kind=realType), dimension(3) :: numd, r
    real(kind=realType) :: LdefoDist, LdefoDist3, Wi, sib(3)

    LdefoDist = tp%Ldef/dist
    Ldefodist3 = LdefoDist**3
    Wi = np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
    Sib(1) = np%Mib(1, 1)*r(1) + np%Mib(1, 2)*r(2) + np%Mib(1, 3)*r(3) + np%bib(1) 
    Sib(2) = np%Mib(2, 1)*r(1) + np%Mib(2, 2)*r(2) + np%Mib(2, 3)*r(3) + np%bib(2) 
    Sib(3) = np%Mib(3, 1)*r(1) + np%Mib(3, 2)*r(2) + np%Mib(3, 3)*r(3) + np%bib(3) 
    numd = numd + Wi*Sib
  end subroutine evalNodeApprox_d

  recursive subroutine setData_d(tp, np)
    ! This performs the forward mode derivative of setData
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType) :: ovrN
    integer(kind=intType) :: i 
    if (np%dnum /= 0) then 
       np%Bib = zero
       np%Mib = zero
       do i=np%l, np%u
          np%Bib = np%Bib + tp%Bib(:, i)
          np%Mib = np%Mib + tp%Mib(:, :, i)
       end do
       np%Bib = np%Bib / np%n
       np%Mib = np%Mib / np%n

       ! Call each of the children
       call setData_d(tp, np%left)
       call setData_d(tp, np%right)
    end if
  end subroutine setData_d

   !  Differentiation of computenodalproperties in forward (tangent) mode (with options noISIZE i4 dr8 r8):
   !   variations   of useful results: *(*tp.bi) *(*tp.mi)
   !   with respect to varying inputs: *xsptr
   !   RW status of diff variables: *xsptr:in *(*tp.bi):out *(*tp.mi):out
   !   Plus diff mem management of: xsptr:in tp:in *tp.bi:in *tp.normals:in
   !                *tp.mi:in *tp.xu:in
  SUBROUTINE COMPUTENODALPROPERTIES_D(tp, initialpoint)
    USE GRIDDATA
    USE GRIDINPUT
    USE COMMUNICATION
    IMPLICIT NONE
    ! Subroutine Variables
    TYPE(TREE_MASTER_RECORD), POINTER :: tp

    LOGICAL :: initialpoint
    ! Working variables
    REAL(kind=realtype), DIMENSION(3, 20) :: points
    REAL(kind=realtype), DIMENSION(3, 20) :: pointsd
    INTEGER(kind=inttype) :: i, j, jj, kk, npts, nelem, ind
    INTEGER(kind=inttype), SAVE :: nptsmax=10
    REAL(kind=realtype) :: facearea, facenormal(3)
    REAL(kind=realtype) :: facearead, facenormald(3)
    REAL(kind=realtype) :: sumarea, sumnormal(3), si(3), ds(3), smean(3)&
         &   , da, eta, r(3), dx(3)
    REAL(kind=realtype) :: sumaread, sumnormald(3), dad
    tp%xub = 0.0_8
    ! This performs the copy and mirroring as required
    DO i=1,tp%n
       j = tp%xuind(i)
       tp%xub(:, i) = xsptrd(3*j-2:3*j)
       tp%xu(:, i) = xsptr(3*j-2:3*j)
    END DO
    ! Now we loop over nodes and fill up Ai, Bi, and Mi
    tp%normalsb = 0.0_8
    tp%normals = zero
    tp%mib = 0.0_8
    tp%mi = zero
    tp%bib = 0.0_8
    tp%bi = zero
    tp%bib = 0.0_8
    tp%normalsb = 0.0_8
    tp%mib = 0.0_8
    pointsd = 0.0_8
    DO i=1,tp%n
       sumarea = zero
       sumnormal = zero
       nelem = tp%nodetoelem(1, i)
       sumnormald = 0.0_8
       sumaread = 0.0_8
       DO jj=1,nelem
          ind = tp%nodetoelem(1+jj, i)
          ! Extract points for this face
          npts = tp%faceptr(ind) - tp%faceptr(ind-1)
          DO kk=1,npts
             pointsd(:, kk) = tp%xub(:, tp%faceconn(tp%faceptr(ind-1)+kk))
             points(:, kk) = tp%xu(:, tp%faceconn(tp%faceptr(ind-1)+kk))
          END DO
          CALL GETELEMENTPROPS_D(points, pointsd, npts, facearea, &
               &                        facearead, facenormal, facenormald)
          ! For face 'ind' how many nodes are on the element?
          dad = facearead/(tp%faceptr(ind)-tp%faceptr(ind-1))
          da = facearea/(tp%faceptr(ind)-tp%faceptr(ind-1))
          sumaread = sumaread + dad
          sumarea = sumarea + da
          sumnormald = sumnormald + dad*facenormal + da*facenormald
          sumnormal = sumnormal + da*facenormal
       END DO
       tp%normalsb(:, i) = (sumnormald*sumarea-sumnormal*sumaread)/&
            &       sumarea**2
       tp%normals(:, i) = sumnormal/sumarea
       IF (.NOT.initialpoint) THEN
          ! Now get the rotation Matrix
          IF (userotations) THEN
             CALL GETROTATIONMATRIX3D_D(tp%normals0(:, i), tp%normals(:, i)&
                  &                              , tp%normalsb(:, i), tp%mi(:, :, i), tp%&
                  &                              mib(:, :, i))
          ELSE
             ! Set identity
             tp%mib(:, :, i) = 0.0_8
             tp%mi(:, :, i) = zero
             tp%mib(1, 1, i) = 0.0_8
             tp%mi(1, 1, i) = one
             tp%mib(2, 2, i) = 0.0_8
             tp%mi(2, 2, i) = one
             tp%mib(3, 3, i) = 0.0_8
             tp%mi(3, 3, i) = one
          END IF
          tp%bib(:, i) = tp%xub(:, i) - tp%xu0(1, i)*tp%mib(:, 1, i) - tp%&
               &         xu0(2, i)*tp%mib(:, 2, i) - tp%xu0(3, i)*tp%mib(:, 3, i)
          tp%bi(:, i) = tp%xu(:, i) - (tp%mi(:, 1, i)*tp%xu0(1, i)+tp%mi(:&
               &         , 2, i)*tp%xu0(2, i)+tp%mi(:, 3, i)*tp%xu0(3, i))
       END IF
    END DO
  END SUBROUTINE COMPUTENODALPROPERTIES_D

#endif
#endif
  function myCreateTree(nodes, cumNodesProc, faceSizesLocal, faceConnLocal) Result(tp)
    use gridInput
    use communication
    implicit none
#ifndef USE_TAPENADE

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

#endif
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), dimension(:, :) :: nodes
    integer(kind=intType), dimension(:) :: cumNodesProc, faceSizesLocal, faceConnLocal

    ! Working
    integer(kind=intType), dimension(:), allocatable :: surfSizesProc, surfSizesDisp
    integer(kind=intType), dimension(:), allocatable :: surfConnProc, surfConnDisp
    real(kind=realType), dimension(:, :), allocatable ::  uniquePts
    integer(kind=intType), dimension(:), allocatable :: link, tempInt, invIndices
    real(kind=realType), dimension(3) :: Xcen, dx
    integer(kind=intType) :: i, j, ind, curInd, ierr, iProc, maxConnectedFace, nPts
    ! Allocate the actual tree: tp means "tree pointer"
    allocate(tp)

    ! Gather the displacements
    allocate(surfSizesProc(nProc), surfSizesDisp(0:nProc), &
         surfConnProc(nProc), surfConnDisp(0:nProc))

    call MPI_allgather(size(faceSizesLocal), 1, MPI_INTEGER, surfSizesProc, 1, MPI_INTEGER, &
         warp_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call MPI_allgather(size(faceConnLocal), 1, MPI_INTEGER, surfConnProc, 1, MPI_INTEGER, &
         warp_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Finish the displacment calc:
    surfSizesDisp(:) = 0
    surfConnDisp(:) = 0
    do i=1, nProc
       surfSizesDisp(i) = surfSizesDisp(i-1) + surfSizesProc(i)
       surfConnDisp(i) = surfConnDisp(i-1) + surfConnProc(i)
    end do

    tp%nFace = surfSizesDisp(nProc)

    ! Allocate space for the full facePtr and faceConn
    allocate(tp%facePtr(0:tp%nFace), tp%faceConn(surfConnDisp(nProc)))

    ! And now for the super magical all-gather-v's!
    call MPI_Allgatherv(&
         faceSizesLocal, size(faceSizesLocal), MPI_INTEGER, &
         tp%facePtr(1:), surfSizesProc, surfSizesDisp, MPI_INTEGER, warp_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call MPI_Allgatherv(&
         faceConnLocal, size(faceConnLocal), MPI_INTEGER, &
         tp%faceConn   , surfConnProc, surfConnDisp, MPI_INTEGER, warp_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Note that we put the result of the "faceSizes" allgather into
    ! 'facePtr' which implies it should be size pointer...which we will
    ! make it into now. Note that facePtr is 0-based which makes things
    ! easier. 
    tp%facePtr(0) = 0
    do i=2, tp%nFace
       tp%facePtr(i) = tp%facePtr(i) + tp%facePtr(i-1)
    end do

    ! We now have to update faceConn since it is currently uses the
    ! indexing from the mirror coordinates so we need up update by the
    ! offset in the Xs nodes
    do iProc=0, nProc-1
       do i=surfConnDisp(iProc)+1, surfConnDisp(iProc+1)
          ! Note that cumNOdesProc is now ONE-BASED since we passed to a
          ! subroutine
          tp%faceConn(i) = tp%faceConn(i) + cumNodesProc(iProc+1)
       end do
    end do

    ! Now we can point-reduce the nodes:
    allocate(uniquePts(3, size(nodes, 2)), link(size(nodes, 2)))
    call pointReduceFast(nodes, size(nodes, 2), symmTol, uniquePts, link, tp%n)

    ! Compute the inverse of the link arrary: this goes from the tree
    ! back to the original "nodes".
    allocate(tp%invLink(tp%n))
    tp%invLink = 0
    do i=1, size(nodes, 2)
       if (tp%invLink(link(i)) == 0) then
          tp%invLink(link(i)) = i
       end if
    end do

    ! Let the user know the actual number of Surface nodes used for
    ! the calc. You will get one of these for each tree.
    if (myid == 0) then
       write(*,"(a)", advance="no") '#--------------------------------#'
       print "(1x)"  
       write(*,"(a)", advance="no") " Unique Surface Nodes : "
       write(*,"(I9,1x)",advance="no") tp%n
       print "(1x)"  
       write(*,"(a)", advance="no") '#--------------------------------#'
       print "(1x)"   
    end if

    ! We will store a copy of the unique data we just generated. This is
    ! required for the original build tree code. 
    tp%dim = 3
    allocate(tp%the_data(tp%dim, tp%n))
    tp%the_data = uniquePts(:, 1:tp%n)

    ! Now we can call the original build_tree function
    call build_tree(tp)
    call labelNodes(tp)
    call setCenters(tp, tp%root)

    ! Copy uniquePoints into Xu0 which will be fixed for all time. Note
    ! that Xu0 is just a permutation of the_data. 
    allocate(tp%Xu0(3, tp%n), tp%Xu(3, tp%n))
    do i=1, tp%n
       tp%Xu0(:, i) = uniquePts(:, tp%indexes(i))
       tp%Xu(:, i) = tp%Xu(:, i)
    end do
    deallocate(uniquePts)

    ! We use the indices used to sort to the tree to compute an inverse
    ! mapping for Bi and Mi.  By doing this first, this removes the
    ! additional layer of indirection during the tree traversals of the
    ! nodal properties to be interpolated. Since this will be done
    ! BILLIONS and BILLIONS of times, it is worth making these accesses
    ! fast and essentially in order in memory.
    allocate(invIndices(tp%n))
    do i=1, tp%n
       invIndices(tp%indexes(i)) = i
    end do

    ! Now we can go through and update the conn. Note that conn is
    ! currently 0-based so we need the plus one when indexing into link
    do i=1, size(tp%faceConn)
       tp%faceConn(i) = invIndices(link(tp%faceConn(i) + 1))
    end do
    ! We also need to compute the nodeToElem pointer....This data
    ! structure contains the number of elements each node is connected
    ! to as well as the indices of the faces. It is stored as 2D array
    ! for simplicitity...we first have to determine the maximum number
    ! of cells connected to each node.
    allocate(tempInt(tp%n))
    tempInt = 0
    do i=1, tp%nFace
       nPts = tp%facePtr(i) - tp%facePtr(i-1)
       do j=1, nPts
          ind = tp%faceConn(tp%facePtr(i-1) + j)
          tempInt(ind) = tempInt(ind) + 1
       end do
    end do

    maxConnectedFace = maxval(tempInt)
    allocate(tp%nodeToElem(maxConnectedFace + 1, tp%n))
    deallocate(tempInt)

    ! Now we loop back again and actually assign the faces to the elements
    tp%nodeToElem = 0
    do i=1, tp%nFace
       nPts = tp%facePtr(i) - tp%facePtr(i-1)
       do j=1, nPts
          ind = tp%faceConn(tp%facePtr(i-1) + j)
          tp%nodeToElem(1, ind) = tp%nodeToElem(1, ind) + 1
          curInd = tp%nodeToElem(1, ind)
          tp%nodeToElem(curInd+1, ind) = i
       end do
    end do

    ! Now we need to compute the one-time mapping that goes between the
    ! unique set of surface nodes and Xs. This will be required for the
    ! sensitivity calc. 

    ! Use zero to indicate an un-assigned node
    allocate(tp%XuInd(tp%n))
    tp%XuInd(:) = 0

    do i=1, size(nodes, 2)
       if (tp%XuInd(invindices(link(i))) == 0) then ! Un-assigned
          tp%XuInd(invindices(link(i))) = i
       end if
    end do

    ! Compute Ldef based on the size of the mesh
    Xcen = zero
    do i=1, tp%n
       Xcen = Xcen + tp%Xu0(:, i)
    end do
    Xcen = Xcen / tp%n

    ! Now get the max dist from Xcen to any node:
    tp%Ldef0 = zero
    do i=1, tp%n
       dx = tp%the_data(:, i) - Xcen
       tp%Ldef0 = max(tp%Ldef0, sqrt(dx(1)**2 + dx(2)**2 + dx(3)))
    end do

    ! Deallocate memory from this routine
    deallocate(link, surfSizesProc, surfSizesDisp, invIndices)
    deallocate(surfConnProc, surfConnDisp)

    ! Allocate Space for the nodal properties:
    allocate(tp%Mi(3, 3, tp%n), tp%Bi(3, tp%n), tp%Ai(tp%n), &
         tp%normals(3, tp%n), tp%normals0(3, tp%n))

    ! Run the compute nodal properties to initialize the normal vectors.
    call initNodalProperties(tp)
    call setAi(tp, tp%root)

    ! Set the input parameters to the module
    tp%Ldef = tp%Ldef0 * LdefFact
    tp%alphaToBexp = alpha**bExp
    tp%errTol = errTol

    ! Define where we will check errors. Note that there is no
    ! underlying reason for this particular sequence. I just made it
    ! up and it seemed to work so I never touch it after. A good thing
    ! would be to plot the actual error profiles and see how the
    ! linear approximation at the following points matches. 

 

    do i=1,NERR-1
       orstar(i) = one/(rstar(i+1) - rstar(i))
    end do

    ! Perform some initialization on the tree
    call computeErrors(tp, tp%root)
  end function myCreateTree

  subroutine initNodalProperties(tp)

    ! This wrapper routine is to be called from python after the surface
    ! mesh is set to compute the initial nodal properties. 

    use gridInput
    use gridData
    use communication
    implicit none

    Type (tree_master_record), Pointer :: tp
    integer(kind=intType) :: ierr

    ! Scatter Xs into our local vector  
    call VecScatterBegin(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecScatterEnd(XsToXsLocal, Xs, XsLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Extract a pointer from Xs to use in the main routine
    call VecGetArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call computeNodalProperties(tp, .True.) 

    ! Restore all the arrays
    call VecRestoreArrayF90(XsLocal, XsPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine initNodalProperties

  subroutine computeNodalProperties(tp, initialPoint)

    use gridData
    use gridInput
    use communication
    implicit none

    ! Subroutine Variables
    Type (tree_master_record), Pointer :: tp
    logical :: initialPoint

    ! Working variables
    real(kind=realType), dimension(3, 20) :: points
    integer(kind=intType) :: i, j, jj, kk, nPts, nElem, ind
    integer(kind=intType) :: nPtsMax = 10
    real(kind=realType) :: faceArea, faceNormal(3)
    real(kind=realType) :: sumArea, sumNormal(3), Si(3), ds(3), sMean(3), da, eta, r(3), dx(3)

    ! This performs the copy and mirroring as required
    do i=1, tp%n
       j = tp%XuInd(i)
       tp%Xu(:, i) = XsPtr(3*j-2:3*j)
    end do

    ! Now we loop over nodes and fill up Ai, Bi, and Mi
    tp%normals = zero
    tp%mi = zero
    tp%bi = zero

    !$AD II-LOOP
    do i=1, tp%n
       sumArea = zero
       sumNormal = zero
       nElem = tp%nodeToElem(1, i)

       do jj=1, nElem 
          ind = tp%nodeToElem(1+jj, i)

          ! Extract points for this face
          nPts = tp%facePtr(ind) - tp%facePtr(ind-1)
          do kk=1, nPts
             points(:, kk) = tp%Xu(:, tp%faceConn(tp%facePtr(ind-1) + kk))
          end do

          call getElementProps(points, nPts, faceArea, faceNormal)

          ! For face 'ind' how many nodes are on the element?
          dA = faceArea / (tp%facePtr(ind) - tp%facePtr(ind-1))
          sumArea = sumArea + dA
          sumNormal = sumNormal + dA*facenormal
       end do
       tp%normals(:, i) = sumNormal / sumArea

       if (initialPoint) then
#ifndef USE_TAPENADE
          tp%normals0(:, i) = tp%normals(:, i)
          tp%Ai(i) = sumArea
#endif
       else
          ! Now get the rotation Matrix
          if (useRotations) then 
             call getRotationMatrix3d(tp%normals0(:, i), tp%normals(:, i), tp%Mi(:, :, i))
          else
             ! Set identity
             tp%Mi(:, :, i) = zero
             tp%Mi(1, 1, i) = one
             tp%Mi(2, 2, i) = one
             tp%Mi(3, 3, i) = one
          end if
          tp%Bi(:, i) = tp%Xu(:, i) - (tp%Mi(:, 1, i)*tp%Xu0(1, i) + tp%Mi(:, 2, i)*tp%Xu0(2, i) + tp%Mi(:, 3, i)*tp%Xu0(3, i))
       end if
    end do
  end subroutine computeNodalProperties
End Module kd_tree


! recursive subroutine computeErrorsNum(tp, np)
!    ! This routine computes an estimate of the error to be induced in
!    ! the denomenator at each tree node. Then, when we are evaluating
!    ! the displacment we can determine if the errors induced by the
!    ! approximation are acceptable or not.

!    implicit none
!    Type (tree_master_record), Pointer :: tp
!    Type (tree_node), Pointer :: np 
!    real(kind=realType) , dimension(3, 20) :: vpts
!    real(kind=realType) , dimension(3, 12) :: fpts
!    real(kind=realType) :: numExact(3, 20), numApprox(3, 20), rr(3), Si(3)
!    real(kind=realType) :: LdefoDist, LdefODist3, dist, W, r(3), Wi
!    integer(kind=intType) :: ii, i, j

!    if (np%dnum /= 0) then 
!       do j=1, NERR
!          numExact= zero
!          numApprox = zero
!          call getSpherePts(np%X, np%radius*tp%rstar(j), vpts, fpts)

!          np%err(j) = zero
!          do ii=1, 20
!             ! Compute exact value:
!             do i=np%l, np%u
!                rr = vpts(:, ii) - tp%Xu0(:, i)
!                dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
!                LdefoDist = tp%Ldef/dist
!                Ldefodist3 = LdefoDist**3
!                Wi = tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
!                Si(1) = tp%Mi(1, 1, i)*r(1) + tp%Mi(1, 2, i)*r(2) + tp%Mi(1, 3, i)*r(3) + tp%bi(1, i) - r(1)
!                Si(2) = tp%Mi(2, 1, i)*r(1) + tp%Mi(2, 2, i)*r(2) + tp%Mi(2, 3, i)*r(3) + tp%bi(2, i) - r(2)
!                Si(3) = tp%Mi(3, 1, i)*r(1) + tp%Mi(3, 2, i)*r(2) + tp%Mi(3, 3, i)*r(3) + tp%bi(3, i) - r(3)

!                numExact(:, ii) = numExact(:, ii) + Wi*Si
!             end do

!             ! And the approx value:
!             rr = vpts(:, ii) - np%X
!             dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
!             LdefoDist = tp%Ldef/dist
!             Ldefodist3 = LdefoDist**3
!             Wi = np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
!             Si(1) = np%Mi(1, 1)*r(1) + np%Mi(1, 2)*r(2) + np%Mi(1, 3)*r(3) + np%bi(1) - r(1)
!             Si(2) = np%Mi(2, 1)*r(1) + np%Mi(2, 2)*r(2) + np%Mi(2, 3)*r(3) + np%bi(2) - r(2)
!             Si(3) = np%Mi(3, 1)*r(1) + np%Mi(3, 2)*r(2) + np%Mi(3, 3)*r(3) + np%bi(3) - r(3)

!             numApprox(:, ii) = numApprox(:, ii) + Wi*Si

!             ! Now set the nodal error:
!             np%NumErr(j) = np%numErr(j) + &
!                  (numExact(1, ii) - numApprox(1, ii))**2 + &
!                  (numExact(2, ii) - numApprox(2, ii))**2 + &
!                  (numExact(3, ii) - numApprox(3, ii))**2  
!          end do

!          ! This the RMS Error:
!          np%numErr(j) = sqrt(np%numErr(j)/20)
!       end do
!       ! Now call the children:
!       call computeErrorsNum(tp, np%left)
!       call computeErrorsNum(tp, np%right)
!    end if
!  end subroutine computeErrorsNum
