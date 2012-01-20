program for_elliot

implicit none

include 'netcdf.inc'

integer, parameter :: nlat_all=480,nlon_all=960,nday=366,nyr_max=31
!integer, parameter :: max_points_section = 35000, max_points_all = 130000
integer, parameter :: max_points_section = 18000, max_points_all = 130000
real*4, allocatable :: data(:,:,:)
integer start_yr,end_yr, nday_yr
character fileroot*100,year*4,var_name*20,arg*10,outdir*100,file*200
character str_lat*20, str_lon*20, fname*200
integer :: fid,status,n,iyr,ilat,jlon,varid,nyr, v
real*4 all_data(4,max_points_section,nday*nyr_max)
character(len=20) :: var_list(1:4)
character(len=80) :: tmax_dir, tmin_dir, precip_dir, solar_dir
integer :: day_start
integer :: all_times(nyr_max*nday)
real*4 :: lat(nlat_all), lon(nlon_all)
logical :: lat_lon_known, n_land_points_known
integer :: counter, n_print_points, print_point_lat(max_points_all), print_point_lon(max_points_all)
integer :: n_procs, proc_num, min_land_point_proc, max_land_point_proc, n_land_points_all, n_land_points_proc
! program for reading in netcdf data for tmax, tmin, solar and precip and outputting to ascii files
!
! written by David McInerney, based on code from Gavin Schmidt
! 
! usage: 
! ./elliot $temp_source $precip_source $solar_source $year_start $year_stop $in_dir $out_dir $n_procs $proc_num

! allocate size of array data (will deallocate later)
allocate(data(nlon_all,nlat_all,nday))

call getarg(1,tmax_dir)
print*,"tmax_dir:", tmax_dir

call getarg(2,tmin_dir)
print*,"tmin_dir:", tmin_dir

call getarg(3,precip_dir)
print*,"precip_dir:", precip_dir

call getarg(4,solar_dir)
print*,"solar_dir:", solar_dir

call getarg(5,arg)
read(arg,'(I4)') start_yr
print*,"start yr:", start_yr

call getarg(6,arg)
read(arg,'(I4)') end_yr
print*,"end yr:", end_yr

call getarg(7,fileroot)
print*,"fileroot:", fileroot

call getarg(8,outdir)
print*,"outdir:", outdir

call getarg(9,arg)
read(arg,'(I4)') n_procs
print*,"n_procs:", n_procs

call getarg(10,arg)
read(arg,'(I4)') proc_num
print*,"proc_num:", proc_num

! list of variables we will read in
var_list = (/ 'precip', 'solar', 'tmax', 'tmin' /)

nyr = end_yr-start_yr + 1

! flag for whether lat and lon have been read in yet
lat_lon_known = .false.

! flag for whether number of land points has been calculated
n_land_points_known = .false.

! index used for time in data array
day_start = 1

! loop over years
do iyr = 1, nyr

! calculate number of days in year
    if ( mod(iyr+start_yr-1,4) .eq. 0 ) then
      nday_yr = 366
    else
      nday_yr = 365
    end if

! loop over variables
  do v = 1, 4
    var_name = var_list(v)

    write(year,'(I4)') iyr+start_yr-1
! calculate file names for netcdf files
    if ( (var_name .eq. 'tmax')) then
      file=trim(fileroot)//trim(tmax_dir)//"/tmax_"//year//".nc"
    else if ( (var_name .eq. 'tmin')) then
      file=trim(fileroot)//trim(tmin_dir)//"/tmin_"//year//".nc"
    else if ( var_name .eq. 'precip' ) then
      file=trim(fileroot)//trim(precip_dir)//"/precip_"//year//".nc"
    else if ( var_name .eq. 'solar' ) then
      file=trim(fileroot)//trim(solar_dir)//"/solar_"//year//".nc"
    end if

! open netcdf file
    status = nf_open(file,nf_nowrite,fid) 
    if (status .ne. 0) then
      print*,"cannot open ", file
      stop 1
    end if

! read in lat and lon
    if (.not.(lat_lon_known)) then

      status = nf_inq_varid(fid,'lat',varid)
      if (status .ne. 0) then
        stop 12
      end if
      status = nf_get_var_real(fid,varid,lat)
      if (status .ne. 0) then
        stop 13
      end if
      status = nf_inq_varid(fid,'lon',varid)
      if (status .ne. 0) then
        stop 22
      end if
      status = nf_get_var_real(fid,varid,lon)
      if (status .ne. 0) then
        stop 23
      end if
      lat_lon_known = .true.
    end if  

! read in "daily_data" variable
    status = nf_inq_varid(fid,'daily_data',varid)
    if (status .ne. 0) then
      print*,"nf_inq_varid fail"
      stop 2
    end if
    status = nf_get_var_real(fid,varid,data)
    if (status .ne. 0) then
      print*,"nf_get_var_real fail"
      stop 3
    end if
    status = nf_close(fid)  ! close
    if (status .ne. 0) then
      print*,"nf_close fail"
      stop 4
    end if
! calculate number of land points, and which ones will be used by this processor
! this will be moved to start of program when i have land/ocean mask file
    if (.not.(n_land_points_known)) then
      counter = 0
      do ilat = 1,nlat_all
        do jlon = 1,nlon_all
! currently determining land/ocean mask from precip data
          if (data(jlon,ilat,1)>=0.) then
            counter = counter + 1
            print_point_lat(counter) = ilat
            print_point_lon(counter) = jlon
          end if
        end do
      end do
      n_land_points_all = counter
      print*,"n_land_points_all",n_land_points_all
      min_land_point_proc = int(1+dble(proc_num-1)*dble(n_land_points_all)/dble(n_procs))
      print*,"min_land_point_proc",min_land_point_proc
      max_land_point_proc = int(dble(proc_num)*dble(n_land_points_all)/dble(n_procs))
      print*,"max_land_point_proc",max_land_point_proc
n_land_points_proc = max_land_point_proc - min_land_point_proc + 1
print*,"n_land_points_proc",n_land_points_proc
      n_land_points_known = .true.
    end if

! copy data from netcdf file to big array of data
    do counter = 1, n_land_points_proc
      ilat = print_point_lat(counter+min_land_point_proc-1)
      jlon = print_point_lon(counter+min_land_point_proc-1)
      all_data(v,counter,day_start:day_start+nday_yr-1) = data(jlon,ilat,1:nday_yr)
    end do
  end do
  print*,"loaded data for year",iyr

! produce day value of form YYDDD
  do n=1,nday_yr
    all_times(day_start+n-1) = (iyr-1)*1000 + n
  end do

  day_start = day_start + nday_yr

end do
print*,"loaded all data"

deallocate(data)

do counter = 1,n_land_points_proc
  ilat = print_point_lat(counter+min_land_point_proc-1)
  jlon = print_point_lon(counter+min_land_point_proc-1) 
  write(str_lat,'(i3.3)') ilat 
  write(str_lon,'(i3.3)') jlon 
  fname = trim(outdir)//"/data_"//trim(str_lat)//"_"//trim(str_lon)//".dat"
  open(unit=1,file=fname)
!  write(1,'(A)')"      LAT     LONG"
!  write(1,20)lat(ilat), lon(jlon)
  write(1,'(A)')"*WEATHER DATA : cell XXXXXXX years 1980 -- 2009"
  write(1,'(A)')"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT"
  write(1,'(A,F9.4,F9.4,A)') "    CI", lat(ilat), lon(jlon), "   -99   -99   -99   -99   -99"
  write(1,'(A)')"@DATE  SRAD  TMAX  TMIN  RAIN"
  do n=1,day_start-1
    write(1,10) all_times(n), all_data(2,counter,n)*0.0864, all_data(3,counter,n)-273.16, & 
                all_data(4,counter,n)-273.16, all_data(1,counter,n)
  end do
  close(1)
end do

10 format(I5.5,F6.1,F6.1,F6.1,F6.1)
!20 format(F9.4,F9.4)

end program

