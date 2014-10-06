Module kd_tree
  use constants
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric works so far.

  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 
  !

  ! .. Parameters ..
  ! you choose this.
  Integer(kind=IntType), Parameter :: bucket_size = 32
  integer(kind=intType), Parameter :: NERR = 12
  ! ..  .. Derived Type Declarations ..  

  ! Global information about the tree pointer to the actual data array
  ! dimensionality and total # of points permuted index into the data,
  ! so that indexes[l..u] of some bucket represent the indexes of the
  ! actual points in that bucket.  root pointer of the tree an
  ! internal tree node the dimension to cut where to cut the dimension
  ! indices of points included in this node, referring back to indices
  ! child pointers 

  Type :: tree_master_record
     real(kind=realType), Dimension (:, :), allocatable :: the_data
     Integer(kind=IntType) :: dim, n
     Integer(kind=IntType), Dimension (:), Pointer :: indexes
     Type (tree_node), Pointer :: root

     ! Additional information for IDW 
     integer(kind=intType) :: maxDepth
     real(kind=realType) :: Ldef
     real(kind=realType) :: alphaToBExp
     real(kind=realType) :: farField 
     real(kind=realType) , dimension(:), pointer :: Ai
     real(kind=realType) , dimension(:, :), pointer :: Bi, Bib
     real(kind=realType) , dimension(:, :, :), pointer :: Mi, Mib
     real(kind=realType) , dimension(:, :), pointer :: Xu0
     real(Kind=realType) :: errTol
     Type(tnp), dimension(:), pointer :: flatNodes
     real(kind=realType), dimension(NERR) :: rstar
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
     real(kind=realType) :: radius
     integer(kind=intType) :: lvl
     integer(kind=intType) :: id
     real(kind=realType) :: err(NERR)
  End Type tree_node

  type :: tnp
     type(tree_node), pointer :: tn
  end type tnp

Contains

  Subroutine destroy_tree(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    Call destroy_node(tp%root)
    Deallocate (tp%indexes)
    Deallocate (tp%the_data)
    Nullify (tp%indexes)
    if (associated(tp%flatNodes)) then 
       deallocate(tp%flatNodes)
    end if
    Return

  Contains

    Recursive Subroutine destroy_node(np)
      ! .. Structure Arguments ..
      Type (tree_node), Pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      Intrinsic ASSOCIATED
      ! ..
      If (ASSOCIATED(np%left)) Then
         Call destroy_node(np%left)
         Deallocate (np%left)
         Nullify (np%left)
      End If
      If (ASSOCIATED(np%right)) Then
         Call destroy_node(np%right)
         Deallocate (np%right)
         Nullify (np%right)
      End If
      Return
    End Subroutine destroy_node
  End Subroutine destroy_tree

  Function create_tree(input_data) Result (master_record)
    ! create the actual tree structure, given an input array of data.
    ! Arguments
    ! .. Function Return Value ..
    Type (tree_master_record), Pointer :: master_record
    ! ..
    ! .. Array Arguments ..
    real(kind=realType), Target :: input_data(:, :)
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic SIZE
    ! ..
    Allocate (master_record)

    master_record%dim = SIZE(input_data, 1)
    master_record%n = SIZE(input_data, 2)

    ! We will store a copy of the data:
    allocate(master_record%the_data(master_record%dim, master_record%n))
    master_record%the_data = input_data
    nullify(master_record%flatNodes)
    Call build_tree(master_record)
    master_record%rstar = (/two, three, four, five, 7.5_realType, 10_realType, &
         20_realType, 40_realType, 80_realType, 160_realType, 320_realType, &
         640_realType/)
  Contains

    Subroutine build_tree(tp)
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

    End Subroutine build_tree

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
  End Function create_tree

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
      np%id = id
      id = id + 1
      tp%maxDepth = max(tp%maxDepth, level)
      if (np%dnum /= 0) then 
         call labelNodes_node(tp, np%left, level+1, id)
         call labelNodes_node(tp, np%right, level+1, id)
      end if
    end subroutine labelNodes_node
  end subroutine labelNodes

  subroutine setXu0(tp, Xu0)
    ! Just set the pointer for Xu0
    implicit none
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in) , pointer, dimension(:, :) :: Xu0

    tp%Xu0 => Xu0
  end subroutine setXu0

  subroutine setAi(tp, Ai)
    implicit none
    ! This routine needs to be called only once . It sets the Ai
    ! information for all the nodes as well as computing the
    ! radius for each of the nodes. 
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in) , pointer, dimension(:) :: Ai

    ! Set the internal pointer to Ai. This shouldn't change
    tp%Ai => Ai
    call setAi_node(tp, tp%root)
  contains 
    recursive subroutine setAi_node(tp, np)
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
         call setAi_node(tp, np%left)
         call setAi_node(tp, np%right)

         ! Now add the sums from the two children 
         np%Ai = (np%left%Ai + np%right%Ai)
      end if
    end subroutine setAi_node
  end subroutine setAi

  subroutine setCenters(tp)
    ! Compute the center of each node. This is just the average X
    ! value. We also compute the minimum radius of the
    ! centroid-centered sphere that will contain all nodes. This
    ! routine computes np%X and np%radius

    implicit none
    Type (tree_master_record), Pointer :: tp
    call setCenters_node(tp, tp%root)

  contains 
    recursive subroutine setCenters_node(tp, np)
      implicit none
      Type (tree_master_record), Pointer :: tp
      Type (tree_node), Pointer :: np
      integer(kind=intType) :: i
      real(kind=realType) :: Xcen(3), Xmax, dist, dx(3)

      ! Sum up all the children
      np%X = zero
      do i=np%l, np%u
         np%X = np%X + tp%the_data(:, tp%indexes(i))
      end do
      np%X = np%X / np%n

      ! Now do the radius
      np%radius = zero
      do i=np%l, np%u
         dx = tp%the_data(:, tp%indexes(i)) - np%X
         dist = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
         if (dist > np%radius) then 
            np%radius = dist
         end if
      end do

      ! Now call the children
      if (np%dnum /= 0) then 
         call setCenters_node(tp, np%left)
         call setCenters_node(tp, np%right)
      end if
    end subroutine setCenters_node
  end subroutine setCenters

  subroutine setData(tp, Bi, Mi)
    ! Update the nodal displacements and rotation matrices based on a
    ! new set of displcements. 
    implicit none
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in), pointer, dimension(:, :) :: Bi
    real(kind=realType), intent(in), pointer, dimension(:, :, :) :: Mi

    ! Set pointers to the nominal data
    tp%Bi => Bi
    tp%Mi => Mi
    Call setData_node(tp, tp%root)

  contains 
    Recursive Subroutine setData_node(tp, np)
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
         call setData_node(tp, np%left)
         call setData_node(tp, np%right)
      end if
    end Subroutine setData_node
  end subroutine setData

  subroutine evalDisp(tp, r, num, den, approxDen)
    implicit none
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(out), dimension(3) :: num
    real(kind=realType), intent(in) :: approxDen
    real(kind=realType), intent(out) :: den

    call evalDisp_node(tp, tp%root, r, num, den, approxDen)
  contains
    recursive subroutine evalDisp_node(tp, np, r, num, den, approxDen)
      implicit none
      ! Subroutine arguments
      Type (tree_master_record), Pointer :: tp
      Type (tree_node), Pointer :: np 
      real(kind=realType), intent(in), dimension(3) :: r
      real(kind=realType), intent(inout), dimension(3) :: num
      real(kind=realType), intent(inout) :: den
      real(kind=realType), intent(in) :: approxDen

      ! Working variables
      real(kind=realType), dimension(3) :: rr
      real(kind=realType) :: dist, err

      if (np%dnum == 0) then 
         ! Leaf node -> Must do exact calc:
         call evalNodeExact(tp, np, r, num, den)
      else
         ! Tree Node: Compute the distance from 'r' to the
         ! center of the node, np%X
         rr(1) = r(1) - np%X(1)
         rr(2) = r(2) - np%X(2)
         rr(3) = r(3) - np%X(3)
         dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)

         if (dist/np%radius < two) then ! Too close...call children
            call evalDisp_node(tp, np%left, r, num, den, approxDen)
            call evalDisp_node(tp, np%right, r, num, den, approxDen)
         else 
            ! Use the first error check:
            call getError(tp, np, dist, err)
            if (err < tp%errTol * approxDen) then
               call evalNodeApprox(tp, np, num, den, dist)
            else
               ! Not good enough error so call children
               call evalDisp_node(tp, np%left, r, num, den, approxDen)
               call evalDisp_node(tp, np%right, r, num, den, approxDen)
            end if
         end if
      end if
    end subroutine evalDisp_node
  end subroutine evalDisp

  subroutine evalNodeExact(tp, np, r, num, den)
    ! Loop over the owned nodes and evaluate the exact numerator and
    ! denomenator at the requesed r. Note that this rouine *may* be
    ! called on *any* node including the root node. When called on the
    ! root node this replicates *precisely* the brute force algorithm
    ! (although not as efficently)
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(out), dimension(3) :: num
    real(kind=realType), intent(out) :: den
    real(kind=realType), dimension(3) :: rr
    real(kind=realType) :: dist, LdefoDist, LdefoDist3
    real(kind=realType) :: Wi
    integer(kind=intType) :: i

    do i=np%l, np%u
       rr(1) = r(1) - tp%Xu0(1, i)
       rr(2) = r(2) - tp%Xu0(2, i)
       rr(3) = r(3) - tp%Xu0(3, i)
       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)
       LdefoDist = tp%Ldef/dist
       Ldefodist3 = LdefoDist**3
       Wi = tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       num(1) = num(1) + Wi*tp%Bi(1, i)
       num(2) = num(2) + Wi*tp%Bi(2, i)
       num(3) = num(3) + Wi*tp%Bi(3, i)
       den = den + Wi
    end do
  end subroutine evalNodeExact

  subroutine evalNodeApprox(tp, np, num, den, dist)
    ! We are far enough away! Whoo! Just use the data on this node
    ! (np). Note that dist must be passed in. It is assumed that the
    ! dist is already available from a calculation to determine
    ! whether to use this approximation or not.
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in) :: dist
    real(kind=realType), intent(out), dimension(3) :: num
    real(kind=realType), intent(out) :: den
    real(kind=realType) :: LdefoDist, LdefoDist3
    real(kind=realType) :: Wi

    LdefoDist = tp%Ldef/dist
    Ldefodist3 = LdefoDist**3
    Wi = np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
    num(1) = num(1) + Wi*np%Bi(1)
    num(2) = num(2) + Wi*np%Bi(2)
    num(3) = num(3) + Wi*np%Bi(3)
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

  subroutine getWiEstimate(tp, r, den)
    ! This routine needs to be only called once. Essentialy what we
    ! are dong is computing just the denomenator (Wi) computation with
    ! a (VERY!) low tolerance. . This is used as a reference value for
    ! the 'real' mesh warp, to be used for checking against the
    ! errors. A hard coded R-limit of 20 is used here. This should be
    ! sufficient to get an good estimate of Wi.
    implicit none

    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(inout) :: den
    den = zero
    call getWiEstimate_node(tp, tp%root, r, den)

  contains
    recursive subroutine getWiEstimate_node(tp, np, r, den)
      implicit none
      ! Subroutine arguments
      Type (tree_master_record), Pointer :: tp
      Type (tree_node), Pointer :: np 
      real(kind=realType), intent(in), dimension(3) :: r
      real(kind=realType), intent(inout) :: den
      real(kind=realType), dimension(3) :: rr
      real(kind=realType) :: dist, LdefoDist, LdefoDist3, err
      integer(kind=intType) :: ii, evals, i

      if (np%dnum == 0) then 
         ! Leaf node. Do regular calc:
         do i=np%l, np%u
            rr = r - tp%Xu0(:, i)
            dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)
            LdefoDist = tp%Ldef/dist
            Ldefodist3 = LdefoDist**3
            den = den + tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
         end do
      else
         rr = r - np%X
         dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)
         call getError(tp, np, dist, err)
         if (dist/np%radius > 5.0) then
            ! Far field calc is ok:
            LdefoDist = tp%Ldef/dist
            Ldefodist3 = LdefoDist**3
            den = den + np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
         else
            ! To far away, recursively call the children
            call getWiEstimate_node(tp, np%left, r, den)
            call getWiEstimate_node(tp, np%right, r, den)
         end if
      end if
    end subroutine getWiEstimate_node
  end subroutine getWiEstimate

  subroutine dryRun(tp, r, approxDen, c)
    ! This routine estimates the number of operations that *would* be
    ! needed to compute the displacement for a given node. No
    ! computations are actually done. It is used to determine the
    ! approximate costs so that a load balancing algorithm can used
    ! to rebalance the nodes. 
    implicit none
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), intent(in) :: approxDen
    real(kind=realType), intent(out) :: c
    c = zero
    call dryRun_node(tp, tp%root, r, approxDen, c)

  contains
    recursive subroutine dryRun_node(tp, np, r, approxDen, c)
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
         ! This three cost accounts for the subroutine calling
         ! overhead and the above dist calc.
         c = c + three
         if (dist/np%radius < two) then ! Too close...call children
            call dryRun_node(tp, np%left, r, approxDen, c)
            call dryRun_node(tp, np%right, r, approxDen, c)
         else 
            ! Use the first error check:
            call getError(tp, np, dist, err)
            if (err < tp%errTol * approxDen) then
               ! This is includes the getError calc
               c = c + 1.5
            else
               ! Not good enough error so call children
               call dryRun_node(tp, np%left, r, approxDen, c)
               call dryRun_node(tp, np%right, r, approxDen, c)
            end if
         end if
      end if
    end subroutine dryRun_node
  end subroutine dryRun

  subroutine computeErrors(tp)
    ! This routine computes an estimate of the error to be induced in
    ! the denomenator at each tree node. Then, when we are evaluating
    ! the displacment we can determine if the errors induced by the
    ! approximation are acceptable or not.

    implicit none
    Type (tree_master_record), Pointer :: tp
    call computeErrors_node(tp, tp%root)
  contains

    recursive subroutine computeErrors_node(tp, np)

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
            call getSpherePts(np%X, np%radius*tp%rstar(j), vpts, fpts)

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
         call computeErrors_node(tp, np%left)
         call computeErrors_node(tp, np%right)
      end if
    end subroutine computeErrors_node
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
    integer(kind=intType) :: bin, i

    dOvrR = dist / np%radius
    bin = 11
    do i=2,11
       if (dOvrR < tp%rstar(i)) then
          bin = i-1
          exit
       end if
    end do
    fact = (dOvrR - tp%rstar(bin))/(tp%rstar(bin+1)-tp%rstar(bin))
    err = (one-fact)*np%err(bin) + fact*np%err(bin+1)
  end subroutine getError

!------------------------------------------------------------
!               DERIVATIVE ROUTINES
!------------------------------------------------------------

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

  subroutine evalDisp_b(tp, r, numb, approxDen, Bib)
    implicit none
    Type (tree_master_record), Pointer :: tp
    real(kind=realType), intent(in), dimension(3) :: r
    real(kind=realType), dimension(3) :: numb
    real(kind=realType), intent(in) :: approxDen
    real(kind=realType), dimension(:, :) :: Bib
    call evalDisp_node_b(tp, tp%root, r, numb, approxDen, Bib)
  contains
    recursive subroutine evalDisp_node_b(tp, np, r, numb, approxDen, Bib)
      implicit none
      ! Subroutine arguments
      Type (tree_master_record), Pointer :: tp
      Type (tree_node), Pointer :: np 
      real(kind=realType), intent(in), dimension(3) :: r
      real(kind=realType), dimension(3) :: numb
      real(kind=realType), intent(in) :: approxDen
      real(kind=realType), dimension(:, :) :: Bib
      ! Working variables
      real(kind=realType), dimension(3) :: rr
      real(kind=realType) :: dist, err

      if (np%dnum == 0) then 
         ! Leaf node -> Must do exact calc:
         call evalNodeExact_b(tp, np, r, numb, Bib)
      else
         ! Tree Node: Compute the distance from 'r' to the
         ! center of the node, np%X
         rr(1) = r(1) - np%X(1)
         rr(2) = r(2) - np%X(2)
         rr(3) = r(3) - np%X(3)
         dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)

         if (dist/np%radius < two) then ! Too close...call children
            call evalDisp_node_b(tp, np%left, r, numb, approxDen, Bib)
            call evalDisp_node_b(tp, np%right, r, numb, approxDen, Bib)
         else 
            ! Use the first error check:
            call getError(tp, np, dist, err)
            if (err < tp%errTol * approxDen) then
               call evalNodeApprox_b(tp, np, numb, dist)
            else
               ! Not good enough error so call children
               call evalDisp_node_b(tp, np%left, r, numb, approxDen, Bib)
               call evalDisp_node_b(tp, np%right, r, numb, approxDen, Bib)
            end if
         end if
      end if
    end subroutine evalDisp_node_b
  end subroutine evalDisp_b

  subroutine evalNodeExact_b(tp, np, r, numb, Bib)
    ! Reverse derviative of evalNodeExact

    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in), dimension(3) :: r, numb
    real(kind=realType), dimension(3) :: rr
    real(kind=realType), dimension(:, :) :: Bib
    real(kind=realType) :: dist, LdefoDist, LdefoDist3
    real(kind=realType) :: Wi
    integer(kind=intType) :: i

    do i=np%l, np%u
       rr(1) = r(1) - tp%Xu0(1, i)
       rr(2) = r(2) - tp%Xu0(2, i)
       rr(3) = r(3) - tp%Xu0(3, i)

       dist = sqrt(rr(1)**2 + rr(2)**2 + rr(3)**2+1e-16)
       LdefoDist = tp%Ldef/dist
       Ldefodist3 = LdefoDist**3
       Wi = tp%Ai(i)*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
       ! tp%Bib(1, i) = tp%Bib(1, i) + wi*numb(1)
       ! tp%Bib(2, i) = tp%Bib(2, i) + wi*numb(2)
       ! tp%Bib(3, i) = tp%Bib(3, i) + wi*numb(3)
       Bib(1, i) = Bib(1, i) + wi*numb(1)
       Bib(2, i) = Bib(2, i) + wi*numb(2)
       Bib(3, i) = Bib(3, i) + wi*numb(3)

    end do
  end subroutine evalNodeExact_b

  subroutine evalNodeApprox_b(tp, np, numb, dist)
    ! Reverse Sensitivity of evalNodeApprox
    implicit none
    Type (tree_master_record), Pointer :: tp
    Type (tree_node), Pointer :: np
    real(kind=realType), intent(in) :: dist
    real(kind=realType), dimension(3) :: numb
    real(kind=realType) :: LdefoDist, LdefoDist3
    real(kind=realType) :: Wi

    LdefoDist = tp%Ldef/dist
    Ldefodist3 = LdefoDist**3
    Wi = np%Ai*(Ldefodist3 + tp%alphaToBexp*Ldefodist3*LdefoDist*LdefoDist)
    np%Bib(1) = np%Bib(1) + wi*numb(1)
    np%Bib(2) = np%Bib(2) + wi*numb(2)
    np%Bib(3) = np%Bib(3) + wi*numb(3)
  end subroutine evalNodeApprox_b

  subroutine setData_b(tp)
    ! This performs the reverse accumulation of the nodeal Bib (and
    ! Mib) into the full Bib array.
    implicit none
    Type (tree_master_record), Pointer :: tp
    
    call setData_node_b(tp, tp%root)
    contains
      recursive subroutine setData_node_b(tp, np)
        implicit none
        Type (tree_master_record), Pointer :: tp
        Type (tree_node), Pointer :: np
        real(kind=realType) :: ovrN
        integer(kind=intType) :: i
        if (np%dnum /= 0) then 
           ovrN = one/ np%n
           do i=np%l, np%u
              tp%Bib(1, i) = tp%Bib(1, i) + np%Bib(1)*ovrN
              tp%Bib(2, i) = tp%Bib(2, i) + np%Bib(2)*ovrN
              tp%Bib(3, i) = tp%Bib(3, i) + np%Bib(3)*ovrN
              !tp%Mib(:, i) = tp%Mib(:, i) + np%Mib*ovrN
           end do

           ! Call each of the children
           call setData_node_b(tp, np%left)
           call setData_node_b(tp, np%right)
        end if
      end subroutine setData_node_b
    end subroutine setData_b
  End Module kd_tree

