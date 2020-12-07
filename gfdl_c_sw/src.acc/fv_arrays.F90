module fv_arrays_mod

  public

  type fv_grid_type
   real, allocatable ::  sin_sg(:,:,:)
	 real, allocatable ::  cos_sg(:,:,:)
   real, allocatable ::  cosa_u(:,:)
	 real, allocatable ::  cosa_v(:,:)
	 real, allocatable ::  cosa_s(:,:)
	 real, allocatable ::  sina_u(:,:)
	 real, allocatable ::  sina_v(:,:)
	 real, allocatable ::  rsin_u(:,:)
	 real, allocatable ::  rsin_v(:,:)
	 real, allocatable ::  rsin2(:,:)
   real, allocatable ::  dx(:,:)
	 real, allocatable ::  dy(:,:)
	 real, allocatable ::  dxa(:,:)
	 real, allocatable ::  dya(:,:)
	 real, allocatable ::  dxc(:,:)
	 real, allocatable ::  dyc(:,:)
	 real, allocatable ::  rdxc(:,:)
	 real, allocatable ::  rdyc(:,:)
	 
	 real, allocatable ::  rarea(:,:)
	 real, allocatable ::  rarea_c(:,:)
	 real, allocatable ::  fC(:,:)
	 logical           :: bounded_domain
	 logical           :: sw_corner
	 logical           :: se_corner
	 logical           :: nw_corner
	 logical           :: ne_corner
  end type fv_grid_type

  type fv_flags_type
   integer ::  npx
	 integer ::  npy
	 integer ::  nord
	 integer ::  grid_type
  end type fv_flags_type
  
  type fv_grid_bounds_type

     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed
     integer :: isc, iec, jsc, jec

     integer :: ng = 3 !default

  end type fv_grid_bounds_type
  
end module fv_arrays_mod
