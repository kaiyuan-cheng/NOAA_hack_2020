program Main

  use netcdf
  use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type, fv_flags_type
  use sw_core_mod, only : c_sw
  #ifdef _OPENACC
  use openacc
  #endif
  implicit none
include 'mpif.h'


  !namelist
  integer :: nx
  integer :: ny
  integer :: nz
  
  real    :: dt = 30.
  real    :: dt2= 15.

  integer :: i, j, k
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: ng = 3 
  integer :: npx, npy, npz

  real, allocatable,   dimension(:,:,:) :: delpc, delp, ptc, pt, u, v, w
  real, allocatable,   dimension(:,:,:) :: uc, vc, ua, va, ut, vt, omga, divgd
  
  logical :: hydrostatic = .False.

  type(fv_grid_type)        :: gridstruct
  type(fv_grid_bounds_type) :: bd
  type(fv_flags_type)       :: flagstruct

  integer :: info, myrank
  integer :: start_time, end_time, cr, cm
  real    :: rate, total_time
  integer :: ncid, vid

! call mpi_init(info)
! call mpi_comm_rank(mpi_comm_world, myrank, info)

  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  rate = REAL(cr)


  print*, "NOTE from Main program: will read data from file "

!  call handle_error(nf90_open(trim("regional_pt.tile7.nc"),nf90_nowrite,ncid))
!  call get_data_i0d(ncid, 'grid_xt', npx)
!  call get_data_i0d(ncid, 'grid_yt', npy)
!  call get_data_i0d(ncid, 'lev', nz)
  npx = 109
  npy = 109
  nz  = 91
!  npx = npx - 6
!  npy = npy - 6
  nx = npx-1
  ny = npy-1
  npz = nz
  print*, "nx,ny,nz=", nx, ny, nz

  isc = 1; iec = nx
  jsc = 1; jec = ny
  isd = isc-ng; ied = iec+ng
  jsd = jsc-ng; jed = jec+ng

  allocate(delp(isd:ied,jsd:jed,nz))
  allocate(delpc(isd:ied,jsd:jed,nz))
  allocate(ptc(isd:ied,jsd:jed,nz))
  allocate(pt(isd:ied,jsd:jed,nz))
  allocate(ua(isd:ied,jsd:jed,nz))
  allocate(va(isd:ied,jsd:jed,nz))
  allocate(ut(isd:ied,jsd:jed,nz))
  allocate(vt(isd:ied,jsd:jed,nz))
  allocate(w(isd:ied,jsd:jed,nz))
  allocate(omga(isd:ied,jsd:jed,nz))
  
  allocate(divgd(isd:ied+1,jsd:jed+1,nz))
  
  allocate(u(isd:ied,jsd:jed+1,nz))
  allocate(v(isd:ied+1,jsd:jed,nz))

  allocate(vc(isd:ied,jsd:jed+1,nz))
  allocate(uc(isd:ied+1,jsd:jed,nz))

  allocate(gridstruct%rarea(isd:ied,jsd:jed))
  allocate(gridstruct%rsin2(isd:ied,jsd:jed))
  allocate(gridstruct%cosa_s(isd:ied,jsd:jed))
  allocate(gridstruct%dxa(isd:ied,jsd:jed))
  allocate(gridstruct%dya(isd:ied,jsd:jed))

  allocate(gridstruct%rarea_c(isd:ied+1,jsd:jed+1))
  allocate(     gridstruct%fC(isd:ied+1,jsd:jed+1))

  allocate(gridstruct%sin_sg(isd:ied,jsd:jed,9))
  allocate(gridstruct%cos_sg(isd:ied,jsd:jed,9))
  
  allocate(gridstruct%cosa_u(isd:ied+1,jsd:jed  ))
  allocate(gridstruct%cosa_v(isd:ied  ,jsd:jed+1))  
  allocate(gridstruct%sina_u(isd:ied+1,jsd:jed  ))
  allocate(gridstruct%sina_v(isd:ied  ,jsd:jed+1))  
  allocate(    gridstruct%dy(isd:ied+1,jsd:jed  ))
  allocate(    gridstruct%dx(isd:ied  ,jsd:jed+1)) 
  allocate(   gridstruct%dxc(isd:ied+1,jsd:jed  ))
  allocate(   gridstruct%dyc(isd:ied  ,jsd:jed+1)) 
  allocate(  gridstruct%rdxc(isd:ied+1,jsd:jed  ))
  allocate(  gridstruct%rdyc(isd:ied  ,jsd:jed+1)) 
  allocate(gridstruct%rsin_u(isd:ied+1,jsd:jed  ))
  allocate(gridstruct%rsin_v(isd:ied  ,jsd:jed+1))   
  
  
  flagstruct%npx = npx
  flagstruct%npy = npy
  flagstruct%grid_type = 0
  flagstruct%nord = 2

  bd%is = isc
  bd%ie = iec
  bd%js = jsc
  bd%je = jec 
  
  bd%isc = isc
  bd%iec = iec
  bd%jsc = jsc
  bd%jec = jec  
  
  bd%isd = isd
  bd%ied = ied
  bd%jsd = jsd
  bd%jed = jed  

  call handle_error(nf90_open(trim("regional_u.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'u', u)
  call handle_error(nf90_open(trim("regional_vc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'vc', vc)
  call handle_error(nf90_open(trim("regional_v.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'v', v)
  call handle_error(nf90_open(trim("regional_uc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'uc', uc)
  call handle_error(nf90_open(trim("regional_w.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'w', w)
  call handle_error(nf90_open(trim("regional_delp.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'delp', delp)
  call handle_error(nf90_open(trim("regional_pt.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'pt', pt)
  call handle_error(nf90_open(trim("regional_ua.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'ua', ua)
  call handle_error(nf90_open(trim("regional_va.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'va', va)

  call handle_error(nf90_open(trim("regional_sin_sg.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'sin_sg', gridstruct%sin_sg)
  call handle_error(nf90_open(trim("regional_cos_sg.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r3d(ncid, 'cos_sg', gridstruct%cos_sg)

  call handle_error(nf90_open(trim("regional_cosa_u.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'cosa_u', gridstruct%cosa_u)
  call handle_error(nf90_open(trim("regional_cosa_v.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'cosa_v', gridstruct%cosa_v)
  call handle_error(nf90_open(trim("regional_sina_u.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'sina_u', gridstruct%sina_u)
  call handle_error(nf90_open(trim("regional_sina_v.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'sina_v', gridstruct%sina_v)
  call handle_error(nf90_open(trim("regional_dx.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dx', gridstruct%dx)
  call handle_error(nf90_open(trim("regional_dy.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dy', gridstruct%dy)
  call handle_error(nf90_open(trim("regional_dxc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dxc', gridstruct%dxc)
  call handle_error(nf90_open(trim("regional_rdxc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rdxc', gridstruct%rdxc)
  call handle_error(nf90_open(trim("regional_dyc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dyc', gridstruct%dyc)
  call handle_error(nf90_open(trim("regional_rdyc.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rdyc', gridstruct%rdyc)
  call handle_error(nf90_open(trim("regional_dya.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dya', gridstruct%dya)
  call handle_error(nf90_open(trim("regional_dxa.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'dxa', gridstruct%dxa)

  call handle_error(nf90_open(trim("regional_rarea.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rarea', gridstruct%rarea)
  call handle_error(nf90_open(trim("regional_rarea_c.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rarea_c', gridstruct%rarea_c)
  call handle_error(nf90_open(trim("regional_fC.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'fC', gridstruct%fC)
  call handle_error(nf90_open(trim("regional_rsin_u.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rsin_u', gridstruct%rsin_u)
  call handle_error(nf90_open(trim("regional_rsin_v.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rsin_v', gridstruct%rsin_v)
  call handle_error(nf90_open(trim("regional_rsin2.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'rsin2', gridstruct%rsin2)
  call handle_error(nf90_open(trim("regional_cosa_s.tile7.nc"),nf90_nowrite,ncid))
  call get_data_r2d(ncid, 'cosa_s', gridstruct%cosa_s)


call system_clock(start_time)
!$ACC DATA copyin(flagstruct, dt2, hydrostatic, bd, gridstruct) copy(delp, delp, ptc, pt, u, v, w, uc , vc, ua, va, omga, ut, vt, divgd)
  do i=1,100
      call c_sw(npz, delpc, delp,  ptc,    &
                pt,    u,    v,    &
                w,   uc,  vc,    &
                ua,   va, omga,    &
                ut,   vt, divgd,   &
                flagstruct%nord,   dt2,  hydrostatic,  .true., bd,  &
                gridstruct, flagstruct)
  enddo
!$ACC END DATA
call system_clock(end_time)
total_time = (end_time-start_time)/rate

  ! call compare_data("pkz", "pkz_out", pkz, pkz_out)
  ! call compare_data("pt", "pt_out", pt, pt_out)
  ! call print_chksum(pkz, "pkz")
  ! call print_chksum(pt, "pt")

  write(*,*) ' elapsed time (secs) = ', total_time

print*, vc(isc+3,jsc+3,:)

contains

  subroutine get_data_r0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r0d

  subroutine get_data_r2d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, dimension(:,:), intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r2d

  subroutine get_data_r3d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, dimension(:,:,:), intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r3d

  subroutine get_data_i0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     integer, intent(out) :: data

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_i0d

  subroutine get_data_i1d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     integer, intent(out) :: data(:)

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_i1d

  subroutine get_data_r1d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     real, intent(out) :: data(:)

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,data))

  end subroutine get_data_r1d



  subroutine get_data_l0d(ncid, field, data)
     integer, intent(in) :: ncid
     character(len=*), intent(in) :: field
     logical, intent(out) :: data
     integer :: idata

     call handle_error(nf90_inq_varid(ncid,trim(field),vid))
     call handle_error(nf90_get_var(ncid,vid,idata))

     if(idata == 0) then
        data = .false.
     else if(idata==1) then
        data = .true.
     else
        call error_handler("field "//trim(field)//" value is neither 0 nor 1")
     endif

  end subroutine get_data_l0d


  subroutine compare_data(str1, str2, data1, data2)
     real, dimension(:,:,:), intent(in) :: data1, data2
     character(len=*), intent(in)       :: str1, str2

     write(*,'(A,ES18.10)' ) "Sum of "//trim(str1)//" = ", sum(data1)
     write(*,'(A,ES18.10)' ) "Sum of "//trim(str2)//" = ", sum(data2)
     write(*,'(A,ES18.10)' ) "max diff of "//trim(str1)//" = ", maxval(data1-data2)


  end subroutine compare_data


  subroutine error_handler(errmsg)
    character(len=*), intent(in) :: errmsg


    write(*,*) "ERROR: ", trim(errmsg)
    call ABORT()

  end subroutine error_handler

  subroutine handle_error(error_flag,err_string)
    ! Simple error handle for NetCDF
    integer,intent(in) :: error_flag
    character(*), intent(in),optional :: err_string
    if ( error_flag  /= nf90_noerr ) then
       write(*,*) 'FATAL ERROR:',nf90_strerror(error_flag)
       if (present(err_string)) write(*,*) trim(err_string)
       call ABORT()
    endif
  end subroutine handle_error

  subroutine print_chksum(data, str)
     real,             intent(in) :: data(:,:,:)
     character(len=*), intent(in) :: str
     integer(kind=8) :: idata(size(data))
     integer(kind=8) :: mold(1), chksum_data, gchksum
     real :: sum_data, gsum_data
     integer :: isc, iec, jsc, jec, error

     idata = transfer(data, mold)
     chksum_data = sum(idata)

     print*, trim(str)//" chksum is ", chksum_data

  end subroutine print_chksum


end program Main
