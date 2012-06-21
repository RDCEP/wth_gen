! -------------------------------

module data_module

  real*4, allocatable :: data(:,:,:)

end module data_module

! -------------------------------

module dirnames

  character(len=80) :: temp_source, precip_source, solar_source, bcsd_type, fileroot

end module dirnames

! -------------------------------

program for_nc_wth_gen

use data_module
use dirnames

implicit none

include 'netcdf.inc'

!integer, parameter :: nlat_all=480,nlon_all=960,nday=366!,nyr_max=31
integer, parameter :: nlat_all=254,nlon_all=720,nday=366!,nyr_max=31
integer, parameter :: max_points_section = 18000, max_points_all = 130000
integer start_yr,end_yr, nday_yr
character var_name*20,arg*10,outdir*100,file*300
!character str_lat*20, str_lon*20, fname*20, full_fname*200, dirname*200, str_start_yr*4, str_end_yr*4, str_proc_num*1
character fname*20, full_fname*200, dirname*200, str_start_yr*4, str_end_yr*4, str_proc_num*1

character :: fname_soil_old*200,fname_soil_new*200
integer :: fid,status,n,iyr,ilat,jlon,varid,nyr, v
!real*4 all_data(4,max_points_section,nday*nyr_max)
real*4, allocatable :: all_data(:,:,:)
character(len=20) :: var_list(1:4)
integer :: day_start
!integer :: all_times(nyr_max*nday)
integer, allocatable :: all_times(:)
real*4 :: lat(nlat_all), lon(nlon_all)
logical :: lat_lon_known, n_land_points_known, ex
integer :: counter, n_print_points, print_point_lat(max_points_all), print_point_lon(max_points_all)
!integer :: grid_ID_arr(max_points_all)
integer :: n_procs, proc_num, min_land_point_proc, max_land_point_proc, n_land_points_all, n_land_points_proc
integer :: grid_ID
real*4 data_time_0(4,nlon_all,nlat_all)
character :: str_grid_ID*7, str_lat*7, str_lon*7
real*4, allocatable :: solar(:), precip(:), tmax(:), tmin(:), tave(:)
integer, allocatable :: time(:)
integer :: n_chunks, nyr_chunk,chunk_start,chunk_end, chunk
integer :: day_count
real*4 :: tave_sum, tave_mean, amp
integer :: month_counter(12)
real*4 :: tave_month_sum(12), tave_month_mean(12)
integer :: month, m, day_of_year

!common /dirnames/ temp_source, precip_source, solar_source, bcsd_type, fileroot
!common /data_array/ data

! program for reading in netcdf data for tmax, tmin, solar and precip and outputting to ascii files
!
! written by David McInerney, based on code from Gavin Schmidt
! 
! usage: 
! ./nc_wth_gen $temp_source $precip_source $solar_source $year_start $year_stop $in_dir $out_dir $n_procs $proc_num

! allocate size of array data (will deallocate later)
!allocate(data(nlon_all,nlat_all,nday))

!@!call getarg(1,temp_source)
!@!print*,"temp_source:", temp_source
!@!
!@!call getarg(2,precip_source)
!@!print*,"precip_source:", precip_source
!@!
!@!call getarg(3,solar_source)
!@!print*,"solar_source:", solar_source
!@!
!@!call getarg(4,bcsd_type)
!@!print*,"bcsd_type:", bcsd_type
!@!
!@!call getarg(5,arg)
!@!read(arg,'(I4)') start_yr
!@!print*,"start yr:", start_yr
!@!
!@!call getarg(6,arg)
!@!read(arg,'(I4)') end_yr
!@!print*,"end yr:", end_yr
!@!
!@!call getarg(7,fileroot)
!@!print*,"fileroot:", fileroot
!@!
!@!call getarg(8,outdir)
!@!print*,"outdir:", outdir
!@!
!@!call getarg(9,fname)
!@!print*,"fname:", fname
!@!
!@!call getarg(10,arg)
!@!read(arg,'(I4)') n_procs
!@!print*,"n_procs:", n_procs
!@!
!@!call getarg(11,arg)
!@!read(arg,'(I4)') proc_num
!@!print*,"proc_num:", proc_num

call getarg(1,arg)
read(arg,'(I4)') start_yr
print*,"start yr:", start_yr

call getarg(2,arg)
read(arg,'(I4)') end_yr
print*,"end yr:", end_yr

call getarg(3,fileroot)
print*,"fileroot:", fileroot

call getarg(4,outdir)
print*,"outdir:", outdir

call getarg(5,fname)
print*,"fname:", fname

call getarg(6,arg)
read(arg,'(I4)') n_procs
print*,"n_procs:", n_procs

call getarg(7,arg)
read(arg,'(I4)') proc_num
print*,"proc_num:", proc_num

print*

!allocate(all_times(nyr*nday))

! open file for outputting tmax, tmin diagnostics
write(str_start_yr,'(i4.4)') start_yr
write(str_end_yr,'(i4.4)') end_yr
write(str_proc_num,'(i1.1)') proc_num
open(unit=2,file=trim(outdir)//'/temp_diag.'//str_start_yr//'_'//str_end_yr//'.p'// &
     str_proc_num//'.txt')
open(unit=3,file=trim(outdir)//'/precip_diag.'//str_start_yr//'_'//str_end_yr//'.p'// &
     str_proc_num//'.txt')

! list of variables we will read in
var_list = (/ 'precip', 'solar', 'tmax', 'tmin' /)

! number of years
nyr = end_yr-start_yr + 1

! read in lat and lon from file
call calc_lat_lon(var_list(1),start_yr,lat,lon)

allocate(data(nlon_all,nlat_all,nday))
do v = 1,4
  var_name = var_list(v)
  call read_data(var_name,start_yr)
  data_time_0(v,:,:) = data(:,:,1)
end do
deallocate(data)

n_chunks = ceiling(float(nyr)/10.)
!n_chunks = ceiling(float(nyr)/1.)
!print*,n_chunks

! calculate number of land points, and which ones will be used by this processor
! this will be moved to start of program when i have land/ocean mask file
! this also copies soil data to new directory
call calc_land_points(data_time_0,proc_num,n_procs, lat, lon, outdir, &
              print_point_lat, print_point_lon,n_land_points_all, &
              min_land_point_proc,max_land_point_proc, n_land_points_proc)

!allocate(all_data(4,n_land_points_proc,nyr*nday))
!allocate(all_times(nyr*nday))
!allocate(all_data(4,n_land_points_proc,nday))
!allocate(all_times(nday))

! index used for time in data array
!day_start = 1

!allocate(all_data(4,chunk_end-chunk_start+1,nyr*nday))
allocate(all_times(nyr*nday))
allocate(time(nyr*nday))
allocate(tmax(nyr*nday))
allocate(tmin(nyr*nday))
allocate(precip(nyr*nday))
allocate(solar(nyr*nday))
allocate(tave(nyr*nday))

print*,"reading data"

! break number of points into smaller chunks to reduce size of all_data array
do chunk = 1, n_chunks

! first and last indexes for this chunk
  chunk_start = int(1+dble(chunk-1)*dble(n_land_points_proc)/dble(n_chunks))
  chunk_end = int(dble(chunk)*dble(n_land_points_proc)/dble(n_chunks))
  print*,"chunk start=",chunk_start,"chunk end=",chunk_end

allocate(all_data(4,chunk_end-chunk_start+1,nyr*nday))

  day_start = 1

! loop over years
  do iyr = 1, nyr

    print*,"  year",iyr

! calculate number of days in year
    if ( mod(iyr+start_yr-1,4) .eq. 0 ) then
      nday_yr = 366
    else
      nday_yr = 365
    end if

    allocate(data(nlon_all,nlat_all,nday))

! loop over variables
    do v = 1, 4
      var_name = var_list(v)

! read data from netcdf file
!      print*,"reading ",var_name," data for year ",iyr
      call read_data(var_name,iyr+start_yr-1)

!      print*,"done"

! copy data from netcdf file to big array of data
!    do counter = 1, n_land_points_proc
      do counter = chunk_start,chunk_end
        ilat = print_point_lat(counter+min_land_point_proc-1)
        jlon = print_point_lon(counter+min_land_point_proc-1)
!      all_data(v,counter,day_start:day_start+nday_yr-1) = data(jlon,ilat,1:nday_yr)
        all_data(v,counter-chunk_start+1,day_start:day_start+nday_yr-1) = data(jlon,ilat,1:nday_yr)
      end do
    end do

    deallocate(data)

! produce day value of form YYDDD
    do n=1,nday_yr
       if ( mod(start_yr/10,2) .eq. 0 ) then
          all_times(day_start+n-1) = (iyr-1)*1000 + n
       else
          all_times(day_start+n-1) = (iyr-1+10)*1000 + n
       end if
!       all_times(day_start+n-1) = mod(start_yr+iyr-1,100)*1000 + n
    end do

    day_start = day_start + nday_yr

  end do

  print*,"writing data"
!do counter = 1,n_land_points_proc
  do counter = chunk_start,chunk_end
!  if (mod(counter,int(0.1*float(n_land_points_proc))).eq.0) print*,counter," (",n_land_points_proc,")"
    ilat = print_point_lat(counter+min_land_point_proc-1)
    jlon = print_point_lon(counter+min_land_point_proc-1) 

!    write( str_lat , "(i3.3)" ) ilat - 1./24.
!    write( str_lon , "(i3.3)" ) jlon - 1./24.
!    str_grid_ID = trim(str_lat)//"_"//trim(str_lon)

!    grid_ID = floor( 12*lon(jlon) - 51840*lat(ilat) + 4661280.5 )
!    grid_ID = floor( 12*(lon(jlon)+180.-1./24.) - 51840*(lat(ilat)-1./24.) + 4661280.5 )
    grid_ID = nint(( 90. - lat(ilat)) *51840 +( 180 +lon(jlon)) *12)
    write(str_grid_ID,'(i7.7)') grid_ID
    dirname = trim(outdir)//"/"//str_grid_ID(1:3)//"/"//str_grid_ID
    inquire (file=dirname,exist=ex)
    if (.not.(ex)) then
      call system('mkdir -p '//trim(dirname))
    end if
    full_fname = trim(dirname)//"/"//trim(fname)
    open(unit=1,file=full_fname)
!    write(1,'(A)')"*WEATHER DATA : cell "//str_grid_ID//" years "//str_start_yr//" -- "//str_end_yr
!    write(1,'(A)')"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT"
!    write(1,'(A,F9.4,F9.4,A)') "    CI", lat(ilat), lon(jlon), "   -99   -99   -99   -99   -99"
!    write(1,'(A)')"@DATE  SRAD  TMAX  TMIN  RAIN"

    tave_sum = 0d0
    day_count = 0
    do n=1,day_start-1

      if (all_data(2,counter-chunk_start+1,n) .ne. 0d0) then
        day_count = day_count + 1
        time(day_count) = all_times(n)
        solar(day_count) = all_data(2,counter-chunk_start+1,n)*0.0864
        tmax(day_count) = all_data(3,counter-chunk_start+1,n)-273.16
        tmin(day_count) = all_data(4,counter-chunk_start+1,n)-273.16
        precip(day_count) = all_data(1,counter-chunk_start+1,n) * 24. * 3600.

        if (tmax(day_count) < tmin(day_count)+0.1) then
          write(2,*) "lat=",lat(ilat),"lon=",lon(jlon),"time=",time(day_count), & 
                     "tmax_orig=",tmax(day_count),"tmin_orig=",tmin(day_count)
          tmax(day_count) = tmin(day_count) + 0.1
        end if

        if (precip(day_count) < 0.) then
          write(3,*) "lat=",lat(ilat),"lon=",lon(jlon),"time=",time(day_count),"precip=",precip(day_count)
          precip(day_count) = 0.
        end if

!    write(1,10) all_times(n), all_data(2,counter,n)*0.0864, all_data(3,counter,n)-273.16, & 
!                all_data(4,counter,n)-273.16, all_data(1,counter,n)

        tave(day_count) = 0.5d0*(tmin(day_count)+tmax(day_count))

      end if

    end do

!print*,tave(1:day_count)
!print*,day_count
!print*

    tave_mean = sum(tave(1:day_count))/dble(day_count)

    tave_month_sum(:) = 0d0
    month_counter(:) = 0
    do n = 1, day_count
      day_of_year = mod(time(n),1000)
      month = ceiling(dble(day_of_year)*12./366.)
      tave_month_sum(month) = tave_month_sum(month) + tave(n)
      month_counter(month) = month_counter(month) + 1
!      print*,"day of year = ", day_of_year,"   month =",month, "   tave(n) =", tave(n),  &
!             "   tave_month_sum(month) =",tave_month_sum(month), "   month_counter(month) =",month_counter(month)
    end do  

!print*,month_counter
!print*,tave_month_sum
    
    do m = 1,12
      tave_month_mean(m) = tave_month_sum(m)/dble(month_counter(m))
    end do

!print*,tave_month_mean
!stop 222

!print*, tave_month_mean

    amp = maxval(tave_month_mean) - minval(tave_month_mean)

!print*,amp

    write(1,'(A)')"*WEATHER DATA : cell "//str_grid_ID//" years "//str_start_yr//" -- "//str_end_yr
    write(1,'(A)')"@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT"
    write(1,'(A,F9.4,F9.4,A,F6.1,F6.1,A)') "    CI", lat(ilat), lon(jlon), "   -99", tave_mean, amp,"   -99   -99"
    write(1,'(A)')"@DATE  SRAD  TMAX  TMIN  RAIN"
    
    do n = 1, day_count
      write(1,10) time(n), solar(n), tmax(n), tmin(n), precip(n)
    end do

    close(1)

!    fname_soil_old = '/gpfs/pads/projects/see/data/dssat/grid_hwsd/'//str_grid_ID//'/SOIL.SOL' 
!    fname_soil_new = trim(dirname)//'/SOIL.SOL'
!    inquire (file=fname_soil_new,exist=ex)
!    if (.not.(ex)) then
!      call system('cp '//fname_soil_old//' '//fname_soil_new)
!    end if

  end do

  deallocate(all_data)
!  deallocate(all_times)

end do
10 format(I5.5,F6.1,F6.1,F6.1,F6.1)

end program

! --------------------------------------------
 
subroutine calc_filename(var_name,year,file)

  use dirnames

  implicit none
 
  character :: str_year*4,var_name*20,file*200
  integer :: year

  write(str_year,'(I4)') year

!@!  if ( (var_name .eq. 'tmax') .or. (var_name .eq. 'tmin') ) then
!@!    file=trim(fileroot)//'/'//trim(var_name)//'/'//trim(temp_source)//"/"//trim(bcsd_type)//"/"// & 
!@!          trim(var_name)//"_"//str_year//".nc"
!@!  else if ( var_name .eq. 'precip' ) then
!@!    file=trim(fileroot)//'/'//trim(var_name)//'/'//trim(precip_source)//"/"//trim(bcsd_type)//"/"// & 
!@!          trim(var_name)//"_"//str_year//".nc"
!@!  else if ( var_name .eq. 'solar' ) then
!@!    file=trim(fileroot)//'/'//trim(var_name)//'/'//trim(solar_source)//"/"//trim(bcsd_type)//"/"// & 
!@!          trim(var_name)//"_"//str_year//".nc"
!@!  end if

  file=trim(fileroot)//'/'//trim(var_name)//'/'//trim(var_name)//"_"//str_year//".nc"

  return
 
end subroutine calc_filename

! --------------------------------------------
 
subroutine calc_lat_lon(var_name,year,lat,lon)
 
  include 'netcdf.inc'

!  implicit none

!  integer, parameter :: nlat_all=480,nlon_all=960,nday=366,nyr_max=31
integer, parameter :: nlat_all=254,nlon_all=720,nday=366,nyr_max=31
  character :: var_name*20, file*200
  real*4 :: lat(nlat_all), lon(nlon_all)
  integer :: fid, status, year


  call calc_filename(var_name,year,file)

  status = nf_open(file,nf_nowrite,fid)
  if (status .ne. 0) then
    print*,"cannot open ", file
    stop 101
  end if
 
! read in lat and lon
  status = nf_inq_varid(fid,'lat',varid)
  if (status .ne. 0) then
    stop 112
  end if 
  status = nf_get_var_real(fid,varid,lat)
  if (status .ne. 0) then
    stop 113
  end if 
  status = nf_inq_varid(fid,'lon',varid)
  if (status .ne. 0) then
    stop 122
  end if 
  status = nf_get_var_real(fid,varid,lon)
  if (status .ne. 0) then
    stop 123
  end if
  
  return
 
end subroutine calc_lat_lon

! --------------------------------------------
 
subroutine read_data(var_name,year)
 
!  implicit none

  use data_module, only : data

  include 'netcdf.inc'

  integer, parameter :: nlat_all=480,nlon_all=960,nday=366,nyr_max=31
  character var_name*20
  integer :: fid, status, varid, year
  character file*200, nc_var_name*20
 
  call calc_filename(var_name,year,file)

  if ( var_name .eq. 'tmax' ) then
    nc_var_name = 'tasmax'
  else if ( var_name .eq. 'tmin' ) then
    nc_var_name = 'tasmin'
  else if ( var_name .eq. 'solar' ) then
    nc_var_name ='rsds'
  else if ( var_name .eq. 'precip' ) then
    nc_var_name = 'pr'
  end if

  print*,file

  status = nf_open(file,nf_nowrite,fid)
  if (status .ne. 0) then
    print*,"cannot open ", file
    stop 1
  end if

  print*,nc_var_name

!  status = nf_inq_varid(fid,'daily_data',varid)
  status = nf_inq_varid(fid,trim(nc_var_name),varid)
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
 
  return
 
end subroutine read_data

! --------------------------------------------

subroutine calc_land_points(data_time_0,proc_num,n_procs, lat, lon, outdir, &
              print_point_lat, print_point_lon,n_land_points_all, &
              min_land_point_proc,max_land_point_proc, n_land_points_proc)

  implicit none

!  integer, parameter :: nlat_all=480,nlon_all=960,nday=366,nyr_max=31
  integer, parameter :: nlat_all=254,nlon_all=720,nday=366,nyr_max=31
  integer, parameter :: max_points_section = 18000, max_points_all = 130000
  integer :: min_land_point_proc,max_land_point_proc, n_land_points_proc
  integer :: n_land_points_all
  integer :: counter, n_print_points, print_point_lat(max_points_all), print_point_lon(max_points_all)
  integer :: ilat, jlon, m
  integer :: proc_num, n_procs, grid_ID
  real*4 data_time_0(4,nlon_all,nlat_all)
  character fname_old*200, dirname*200,outdir*100,str_grid_ID*7
  logical :: ex, in_soil_grid
  real*4 :: lat(nlat_all), lon(nlon_all)
!  integer :: grid_ID_arr(max_points_all)
integer, parameter :: n_soil_grid = 1971846
!integer, parameter :: n_soil_grid = 118545
integer :: soil_grid(n_soil_grid)
integer :: data_temp(3,max_points_all)
real*4 :: lon_val

!@!!  open(unit=1,file='/gpfs/pads/projects/see/data/dssat/grid_hwsd.txt')
!@!  open(unit=1,file='/autonfs/home/dmcinern/new_wth_gen/grid_id_globe')
!@!  read(1,*) soil_grid 
!@!  close(1)
!@!
!@!inquire (file='grid_id_count.txt',exist=ex)
!@!if (ex) then
!@!  print*, "reading land points"
!@!  open(unit=1,file='grid_id_count.txt')
!@!  read(1,*) counter
!@!  close(1)
!@!  open(unit=1,file='grid_id_list.txt')
!@!  read(1,*) data_temp(:,:counter)
!@!  close(1)
!@!  print_point_lat(:counter) = data_temp(2,:counter)
!@!  print_point_lon(:counter) = data_temp(3,:counter)
!@!
!@!else

  print*,"calculating land points"

!@!  open(unit=1,file='grid_id_list.txt')

  counter = 0
  do ilat = 1,nlat_all
    do jlon = 1,nlon_all


!lon_val = lon(jlon) + 180.
 
!      grid_ID = floor( 12*lon_val - 51840*lat(ilat) + 4661280.5 )
      grid_ID = floor( 12*(lon(jlon)+180.-1./24.) - 51840*(lat(ilat)-1./24.) + 4661280.5 )
      write(str_grid_ID,'(i7.7)') grid_ID

!        if (data_time_0(1,jlon,ilat)>=-999.) then
!          if (data_time_0(2,jlon,ilat)>=-999.) then
!            if (data_time_0(3,jlon,ilat)>=-999.) then
!              if (data_time_0(4,jlon,ilat)>=-999.) then

        if (data_time_0(1,jlon,ilat)<=999999.) then
          if (data_time_0(2,jlon,ilat)<=999999.) then
            if (data_time_0(3,jlon,ilat)<=999999.) then
              if (data_time_0(4,jlon,ilat)<=999999.) then

!@!                in_soil_grid = .false.
!@!                do m = 1, n_soil_grid
!@!                  if (soil_grid(m) .eq. grid_ID) then
!@!                    in_soil_grid = .true.
!@!                    exit
!@!                  endif
!@!                end do

in_soil_grid = .true.

                if (in_soil_grid) then

                  counter = counter + 1

                  print_point_lat(counter) = ilat
                  print_point_lon(counter) = jlon
 !                 grid_ID_arr(counter) = grid_ID

write(1,'(I7,I7,I7)') counter, print_point_lat(counter), print_point_lon(counter)
!print*, counter, grid_ID, lat(ilat), lon(jlon), print_point_lat(counter), print_point_lon(counter)
!stop

                end if
              end if
            end if
          end if
        end if  

    end do
  end do

!@!  close(1)
!@!  open(unit=1,file='grid_id_count.txt')
!@!  write(1,'(I7)') counter

!@!STOP 2222

!@!end if 

  n_land_points_all = counter
  print*,"  n_land_points_all",n_land_points_all
  min_land_point_proc = int(1+dble(proc_num-1)*dble(n_land_points_all)/dble(n_procs))
  print*,"  min_land_point_proc",min_land_point_proc
  max_land_point_proc = int(dble(proc_num)*dble(n_land_points_all)/dble(n_procs))
  print*,"  max_land_point_proc",max_land_point_proc
  n_land_points_proc = max_land_point_proc - min_land_point_proc + 1
  print*,"  n_land_points_proc",n_land_points_proc
  
  return

end subroutine calc_land_points

