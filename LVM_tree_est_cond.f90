subroutine tree_establishment_conditions
    ! This subroutine generates the tree establishment inundation conditions array for use in the establishment subroutine.
    ! This subroutine determines if there is any period of time during the year in which
    ! tree establishment conditions are met at each grid cell.
    ! A tree establishment condition is set to 1 if at any point from March 1 throuh Aug 16, a grid cell has two weeks
    ! of dry land (depth <= 0) followed by 2 weeks in which the water depth is no deeper than 14 cm.
    !
    ! This subroutine was originally part of ICM-Hydro (for MP17 and MP23) and ported over here for MP29

    ! Per Scott Duke-Sylvester's original Python code comments in MP17 code:
    !    # This function determines the probability that a bottomland hardwood forest species will become
    !    # established at the current location. The probability of a bottomland hardwood forest species becoming
    !    # established is depended on three major factors. First, the salinity must be below 1.0 ppt.
    !    # Second, there must be a period of 14 days with no flooding followed by a period of 14 with
    !    # water depths below 14 cm. Finally, the height of the habitat above mean water level determines
    !    # the final probability.
    !    
    !    
    !    # I'm not entirely happy with the current state of this function because the computation
    !    # of the basic establishment conditions (14 day no flood, 14 days water depth < 14 cm) is
    !    # handled external to the model. The problem is that the current approach divides
    !    # the responsibilities for representing the ecology of upland forest species (now bottomland hardwood forest).
    !    # An external program determines the establishment conditions while the model proper handles
    !    # the final computation of the probability of establishment. I would like to fix this.

    ! We have now moved this over to ICM-LAVegMod and we think Scott would approve.
    ! The water level controls are calculated here, the salinity control and height above MWL are applied in the *mort_est_prob* subroutine.


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
    real(sp), dimension(:,:), allocatable :: grid_eff_dep_daily                     ! array with daily water effective depth for each veg grid cell
    integer, dimension(:), allocatable :: drypast_flag                              ! flag (1 or 0) to determine if last 14 days were dry
    integer, dimension(:), allocatable :: shallowfuture_flag                        ! flag (1 or 0) to determine if future 14 days are dry or only have shallow inundation
    integer, dimension(:), allocatable :: tree_est_flag                             ! combined flags to see if both conditions are met (dry past AND dry future)
    integer, dimension(:), allocatable :: month_DOY                                 ! array holding the starting index for each month in a daily timeseries
    character*4 :: year                                                             ! calendar year of model run - used in output file name
    real(sp) :: grid_dry_depth                                                      ! dry depth offset to account for variability of elevation within veg grid cell 
                                                                                    !   - if set to 0.3 meter, than the entire grid will not be considered inundated until
                                                                                    !       the depth is at least 0.3 meter
                                                                                    !   - if set to 0.0, then the entire veg grid cell will be considered inundated when the
                                                                                    !       water surface elevation is just greater than the average elevation of the grid cell
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
    else
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


    grid_dry_depth = 0.3

    do g=1,ngrid 
        comp = grid_comp(g) 
        if (comp > 0) then                                                                      ! check that grid cell has an allowable ICM-Hydro compartment ID
            do j = 1,simdays
                grid_eff_dep_daily(g,j) = stage_daily(j,comp) - grid_elev(g) - grid_dry_depth   ! map compartment stage values to grid cells for each day and convert to depth
            enddo
            
            ! Loop through days at each grid cell and determine tree establishment criteria is met
            do jj = firstday,lastday          
                jjj=jj-firstday+1                                                               ! convert day of year to day of tree establishment window array
                drypast_flag(jjj) = 1                                                           ! initialize  day's drypast flag to values of 1
                shallowfuture_flag(jjj) = 1                                                         ! initialize  day's dryfuture flag to values of 1
                
                ! Loop over past two weeks and determine if any past day is wet
                do dd=0,13
                    if (grid_eff_dep_daily(g,jj-dd) <= -0.0) then                               ! effective depth is negative, so grid is dry
                        drypast_flag(jjj) = drypast_flag(jjj)*1                                 ! loop over past two weeks and determine if past day was dry 
                    else
                        drypast_flag(jjj) = drypast_flag(jjj)*0                                 ! if any day of the past 2 weeks was wet set drypast_flag to 0
                    endif
                
                    if (grid_eff_dep_daily(g,jj+dd) <= -0.14) then                              ! loop over next two weeks and determine if any future day is shallower than 14 cm, or dry
                        shallowfuture_flag(jjj) = shallowfuture_flag(jjj)*1
                    else
                        shallowfuture_flag(jjj) = shallowfuture_flag(jjj)*0                     ! if any day of the next 2 weeks is flooded by more than 14 cm, shallowfuture_flag is set to 0
                    endif
                enddo
                tree_est_flag(jjj) = shallowfuture_flag(jjj)*drypast_flag(jjj)
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
        write(year,'(I0)') start_year + elapsed_year - 1
        open(unit=903, file='veg/'//trim(adjustL(fnc_tag))//'_'//'N'//'_'//year//'_V_tree_est.csv')
        write(903,'(A)') 'GridCellID,TreeEstablishmentCondition'
        do g=1,ngrid
            write(903,3458) g,tree_establishment(g)
        enddo
        close(903)
    endif
    
   
    
    deallocate(grid_dep_daily)
    deallocate(drypast_flag)
    deallocate(dryfuture_flag)
    deallocate(tree_est_flag)
    
3458    format(I0,',',I0)
    
    return
end