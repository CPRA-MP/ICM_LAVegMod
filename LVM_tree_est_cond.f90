

subroutine tree_establishment_conditions
    ! This subroutine generates the tree establishment conditions array for use in the Vegetation model.
    ! This subroutine determines if there is any period of time during the year in which
    ! tree establishment conditions are met at each 500-m grid cell.
    ! A tree establishment condition is set to 1 if at any point from March 1 throuh Aug 16, a grid cell has two weeks
    ! of dry land (depth <= 0) followed by 2 weeks in which the water depth is no deeper than 10 cm.
    !
    ! This subroutine was originally part of ICM-Hydro (for MP17 and MP23) and ported over here for MP29

      use params      

      implicit none

    ! global arrays updated by subroutine:
    !       tree_establishment

    ! local variables:
    real(sp) :: leapcheck                                                           ! decimal year
    integer :: j                                                                    ! iterator over all 365/366 simulation days
    integer :: g                                                                    ! iterator over number of veg grid cells
    integer :: jj                                                                   ! iterator over days during tree establishment window
    integer :: jjj                                                                  ! day of establishment window converted from simulation day
    integer :: dd                                                                   ! number of days in moving window
    integer :: comp                                                                 ! compartment number for veg grid cell (from lookup table)
    integer :: firstday                                                             ! first day of tree establishment window
    integer :: lastday                                                              ! last day of tree establishment window
    integer :: simdays                                                              ! number of days in year currently simulated (either 365 or 366)
    integer :: thresholdlength                                                      ! number of days to analyze for tree establishment conditions
    real(sp), dimension(:,:), allocatable :: grid_dep_daily                         ! array with daily water depth for each veg grid cell
    integer, dimension(:), allocatable :: drypast_flag                              ! flag (1 or 0) to determine if last 14 days were dry
    integer, dimension(:), allocatable :: dryfuture_flag                            ! flag (1 or 0) to determine if future 14 days have less than 10 cm of ponding
    integer, dimension(:), allocatable :: tree_est_flag                             ! combined flags to see if both conditions are met (dry past AND dry future)
    integer, dimension(:), allocatable :: month_DOY                                 ! array holding the starting index for each month in a daily timeseries
    
    
    tree_establishment = 0                                                          ! initialize entire tree_establishment array to 0
    
   allocate(month_DOY(12))
    
    

    month_DOY(1) = 1
    month_DOY(2) = 32
    month_DOY(3) = 60
    month_DOY(4) = 91
    month_DOY(5) = 121
    month_DOY(6) = 152
    month_DOY(7) = 182
    month_DOY(8) = 213
    month_DOY(9) = 244
    month_DOY(10) = 274
    month_DOY(11) = 305
    month_DOY(12) = 335
    
    ! Check if current year is a leap year
    if ( (start_year + elapsed_year)/4.0 > floor((start_year + elapsed_year)/4.0) )then
        simdays = 365
    else:
        simdays = 366
        do j = 3,12
            month_DOY(j) = month_DOY(j)+1                                           ! Update first day of month for March through December during a leap year
        enddo
    endif


    ! Set first and last day of tree establishment window
    firstday = month_DOY(3)
    lastday = month_DOY(8)+16
    thresholdlength = lastday - firstday + 1
    
    ! Allocate temporary array to be of length equal to number of grid cells - these are deallocated at end of this subroutine
    allocate(grid_dep_daily(ngrid,simdays))
    allocate(drypast_flag(thresholdlength))
    allocate(dryfuture_flag(thresholdlength))
    allocate(tree_est_flag(thresholdlength))


    do g=1,ngrid 
        comp = grid_comp(g) 
        if (comp > 0) then                                                          ! check that grid cell has an allowable ICM-Hydro compartment ID
            if (comp > 0) then                                                      ! check that grid cell has an allowable ICM-Hydro compartment ID
                do j = 1,simdays
                    grid_dep_daily(g,j) = stage_daily(j,comp) - grid_elev(g)        ! map compartment stage values to grid cells for each day and convert to depth
                enddo
            
            ! Loop through days at each grid cell and determine tree establishment criteria is met
            do jj = firstday,lastday          
                jjj=jj-firstday+1                                                   ! convert day of year to day of tree establishment window array
                drypast_flag(jjj) = 1                                               ! initialize  day's drypast flag to values of 1
                dryfuture_flag(jjj) = 1                                             ! initialize  day's dryfuture flag to values of 1
                
                ! Loop over past two weeks and determine if any past day is wet
                do dd=0,13
                    if (grid_dep_daily(g,jj-dd) <= -0.30) then
                        drypast_flag(jjj) = drypast_flag(jjj)*1
                    else
                        drypast_flag(jjj) = drypast_flag(jjj)*0                     ! If any day of the past 2 weeks is wet, drypast_flag is set to 0
                    endif
                
                    if (grid_dep_daily(g,jj+dd) <= -0.20) then                      ! loop over next two weeks and determine if any future day is wet                  
                        dryfuture_flag(jjj) = dryfuture_flag(jjj)*1
                    else
                        dryfuture_flag(jjj) = dryfuture_flag(jjj)*0                 ! if any day of the next 2 weeks is flooded by more than 10 cm, dryfuture_flag is set to 0
                    endif
                    
                enddo
                
                tree_est_flag(jjj) = dryfuture_flag(jjj)*drypast_flag(jjj)
            
            enddo
            
            ! Loop over cell's timeseries of flags and set equal to 1 if any daily flags equal 1
            do jj = firstday,lastday                                            
                jjj = jj-firstday+1
                if (tree_est_flag(jjj) > 0) then
                    tree_establishment(g) = 1
                endif
            enddo
        endif
    enddo
    
    ! write csv of tree establishment conditions if intermediate files are being written
    if (write_intermediate_files == 1) then
        open(unit=903, file='veg/'//trim(adjustL(fnc_tag))//'_'//'N'//'_'//year//'_V_tree_est.csv')
        write(903,'(A)') 'GridCellID,TreeEstablishmentCondition'
        do g=1,ngrid
            write(903,3458) g,tree_establishment(g)
        enddo
        close(903)
    endif
    
    
    write(903,'(A)') trim(adjustL(veg_coverage_file_header))
    do g = 1,ngrid
        write(903,'(I0,','I0)') g, coverages(g,:)
    end do
    close(903)
    
    
    deallocate(grid_dep_daily)
    deallocate(drypast_flag)
    deallocate(dryfuture_flag)
    deallocate(tree_est_flag)
    
3458    format(I0,',',I0)
    
    return
end