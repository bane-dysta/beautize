! SPDX-License-Identifier: MIT
! This adapter calls the LGPL-3.0-or-later GFN-FF library.
module beautize_bridge
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use gfnff_interface, only: gfnff_data
#ifdef _OPENMP
  use omp_lib, only: omp_set_num_threads, omp_set_dynamic
#endif
  implicit none
contains
  subroutine bt_gfnff_threads(n,status) bind(C)
    integer(c_int),value :: n
    integer(c_int),intent(out) :: status
    status=0
    if(n<1) then
      status=-1
      return
    end if
#ifdef _OPENMP
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(n)
#else
    if(n/=1) status=-1
#endif
  end subroutine

  function bt_gfnff_create(n,at,xyz,charge,bonds,level,status) result(handle) bind(C)
    integer(c_int),value :: n,charge,level
    integer(c_int),intent(in) :: at(n)
    real(c_double),intent(in) :: xyz(3,n),bonds(n,n)
    integer(c_int),intent(out) :: status
    type(c_ptr) :: handle
    type(gfnff_data),pointer :: calc
    integer :: io,i,j,k
    handle = c_null_ptr
    status = -1
    if (n < 1) return
    if (any(at < 1) .or. any(at > 103)) return
    if (.not.all(ieee_is_finite(xyz))) return
    allocate(calc)
    calc%write_topo = .false.
    calc%restart = .false.
    allocate(calc%userinput)
    calc%userinput%bond_orders = bonds
    call calc%init(n,at,xyz,ichrg=charge,printlevel=level,iostat=io)
    if (io /= 0) then
      status = io
      call calc%deallocate()
      deallocate(calc)
      return
    end if
    ! Audit the actual final neighbor list, NOT a copy of the supplied input.
    do i=1,n
      if (calc%neigh%nb(calc%neigh%numnb,i,1) /= count(bonds(:,i)>0.0_c_double)) io=-21
      do k=1,calc%neigh%nb(calc%neigh%numnb,i,1)
        j=calc%neigh%nb(k,i,1)
        if (j<1 .or. j>n) then
          io=-21
        else if (bonds(j,i)<=0.0_c_double) then
          io=-21
        end if
      end do
    end do
    if (io /= 0) then
      status=io
      call calc%deallocate()
      deallocate(calc)
      return
    end if
    status = 0
    handle = c_loc(calc)
  end function

  subroutine bt_gfnff_destroy(handle) bind(C)
    type(c_ptr),value :: handle
    type(gfnff_data),pointer :: calc
    if (.not.c_associated(handle)) return
    call c_f_pointer(handle,calc)
    call calc%deallocate()
    deallocate(calc)
  end subroutine

  subroutine bt_gfnff_evaluate(handle,n,at,xyz,energy,gradient,status) bind(C)
    type(c_ptr),value :: handle
    integer(c_int),value :: n
    integer(c_int),intent(in) :: at(n)
    real(c_double),intent(in) :: xyz(3,n)
    real(c_double),intent(out) :: energy,gradient(3,n)
    integer(c_int),intent(out) :: status
    type(gfnff_data),pointer :: calc
    integer :: io
    status=-1
    energy=0.0_c_double
    gradient=0.0_c_double
    if (.not.c_associated(handle)) return
    if (.not.all(ieee_is_finite(xyz))) return
    call c_f_pointer(handle,calc)
    call calc%singlepoint(n,at,xyz,energy,gradient,iostat=io)
    status=io
    if (.not.ieee_is_finite(energy) .or. .not.all(ieee_is_finite(gradient))) status=-22
  end subroutine

  subroutine bt_gfnff_graph(handle,n,adj) bind(C)
    type(c_ptr),value :: handle
    integer(c_int),value :: n
    integer(c_int),intent(out) :: adj(n,n)
    type(gfnff_data),pointer :: calc
    integer :: i,j,k
    adj=0
    if (.not.c_associated(handle)) return
    call c_f_pointer(handle,calc)
    do i=1,n
      do k=1,calc%neigh%nb(calc%neigh%numnb,i,1)
        j=calc%neigh%nb(k,i,1)
        adj(j,i)=1
      end do
    end do
  end subroutine

  ! Explicit C ABI avoids relying on compiler-specific LAPACK character ABIs.
  subroutine bt_symmetric_eigen(n,a,w,status) bind(C)
    integer(c_int),value :: n
    real(c_double),intent(inout) :: a(n,n)
    real(c_double),intent(out) :: w(n)
    integer(c_int),intent(out) :: status
    real(c_double),allocatable :: work(:)
    integer :: info,lwork
    external :: dsyev
    if(n==0) then
      status=0
      return
    end if
    lwork=max(1,3*n)
    allocate(work(lwork))
    call dsyev('V','U',n,a,n,w,work,lwork,info)
    status=info
  end subroutine
end module
