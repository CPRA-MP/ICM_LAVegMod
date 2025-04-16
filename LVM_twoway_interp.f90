subroutine twoway_interp(variable1, variable2, table, variable1bins, var1bin_n, variable2bins, var2bin_n, yint)

    ! inputs to the subroutine:
    !   variable1 = the value of the input variable specific to the grid cell. For 2023 LAVegMod, variable1 is either mean annual salinity or water level variability. 
    !   variable2 = the value of the input variable specific to the grid cell. For 2023 LAVegMod, variable1 is either mean annual salinity or water level variability. 
    !   table = establishment or mortality table for one species (2D)
    !   variable1bins = the slice of est_X_bins (or mort_X_bins, est_Y_bins, mort_Y_bins) for one species (1D)
    !   var1bin_n = the number of 'bins' used to discretizethe the interpolation table in the variable1 dimension
    !   var2bin_n = dummy variable to hold the number of 'bins' used to discretizethe the interpolation table in the variable2 dimension
    !   variable2bins = the slice of est_Y_bins (or mort_X_bins, est_Y_bins, mort_Y_bins) for one species (1D)
    
    ! output of the subroutine:
    !   yint = the interpolated value of establishment or mortality probability for the value of variable1 and variable2

    ! global arrays updated by subroutine:
    !   establish_P or mortality_P -- but not directly updated here, but rather the output updates one of those arrays

    ! global arrays used by subroutine:
    !   n_X_bins
    !   n_Y_bins

    ! This subroutine interpolates the establishment or mortality probability based on two inputs. It acts on one species for one grid cell.
    ! For 2023 LAVegMod, this subroutine applies to swamp forest, thick and thin floating marsh, and emergent wetland (fresh, intermediate, brackish, and saline). 
    ! The est and mort probabilities of those species are a function of mean annual salinity (sal_av_yr) and water level variabilty (wlv_smr). 


    use params
    implicit none

    ! dummy local variables populated with arrays passed into subroutine
    integer,intent(in) :: var1bin_n                             ! dummy variable to hold the number of 'bins' used to discretizethe the interpolation table in the variable1 dimension
    integer,intent(in) :: var2bin_n                             ! dummy variable to hold the number of 'bins' used to discretizethe the interpolation table in the variable2 dimension
    real(sp),intent(in) :: variable1                            ! dummy variable to hold value of first variable used in 2d interpolation. For the 2023 Master Plan, variable1 was salinity on the establishment/mortality tables.
    real(sp),intent(in) :: variable2                            ! dummy variable to hold value of second variable used in 2d interpolation. For the 2023 Master Plan, variable2 was water level variability on the establishment/mortality tables.
    real(sp),dimension(var1bin_n,var2bin_n),intent(in) :: table ! dummy variable to hold two-dimensional table interpolation is performed on. variable1 is the X-axis of this table, variable2 is the Y-axis of this table.
    real(sp),dimension(var1bin_n),intent(in) :: variable1bins   ! dummy variable to hold the values (in same units of variable1) defining the lower bound of the discretization 'bins' of the interpolation table in the variable1 dimension
    real(sp),dimension(var2bin_n),intent(in) :: variable2bins   ! dummy variable to hold the values (in same units of variable2) defining the lower bound of the discretization 'bins' of the interpolation table in the variable2 dimension    
    real(sp),intent(inout) :: yint                              ! dummy variable to hold the final interpolated value, returned to parent subroutine
    
    ! local variables 
    integer   :: ig                                             ! iterator over Veg grid cells
    integer   :: ib                                             ! iterator over n_X_bins or n_Y_bins
    real(sp)  :: above                                          ! value of establisment table above the input value
    real(sp)  :: below                                          ! value of establisment table below the input value
    real(sp)  :: left                                           ! value of establisment table to the left of the input value
    real(sp)  :: right                                          ! value of establisment table to the right of the input value
    real(sp)  :: min_dif                                        ! the smallest difference of the differences between variable1 or variable2 and each bin values
    real(sp)  :: dif                                            ! difference between variable1 or variable2 and each bin value
    integer   :: closest_index                                  ! index within variable1bins or variable2bins for either the above/below value or left/right
    real(sp)  :: y1                                             ! variable used in the linear interpolation formula 
    real(sp)  :: x1                                             ! variable used in the linear interpolation formula 
    real(sp)  :: y2                                             ! variable used in the linear interpolation formula 
    real(sp)  :: x2                                             ! variable used in the linear interpolation formula 
    real(sp)  :: xint                                           ! variable used in the linear interpolation formula 
    real(sp)  :: yint_varY1                                     ! variable used in the linear interpolation formula 
    real(sp)  :: yint_varY2                                     ! variable used in the linear interpolation formula 

    ! Find the variable1 bin value closest to variable1    
    min_dif = 3000                                              ! arbitary value, just must be larger than any expected differences
    dif = 0
    closest_index = -9999
    do ib = 1, var1bin_n                                        ! loop through the bin values
        dif = abs(variable1bins(ib) - variable1)                ! calculate the absolte value of the difference between each bin and the given value
        if (dif < min_dif) then
            closest_index = ib                                  ! index for the value closest to the given value
            min_dif = dif                                       ! absolute value of the smallest difference between the bin values and the given value
        end if 
    end do

    ! Figure out if the min_dif bin value is above or below the given value and assign above and below 
    if ( (variable1bins(closest_index)-variable1) < 0) then
        left = variable1bins(closest_index)
        right = variable1bins(closest_index+1)
    elseif ( (variable1bins(closest_index)-variable1) == 0) then    ! same value no interpolation needed 
        left = variable1bins(closest_index)
        right = variable1bins(closest_index)
    else ! (variable1bins(closest_index)-variable1) < 0 then 
        left = variable1bins(closest_index)
        right = variable1bins(closest_index-1)
    end if 

    ! Find the variable2 bin value closest to variable2   
    min_dif = 3000                                              ! arbitary value, just must be larger than any expected differences
    dif = 0 
    closest_index = -9999                                         
    do ib = 1, var2bin_n                                        ! loop through the bin values
        dif = abs(variable2bins(ib) - variable2)                ! calculate the absolte value of the difference between each bin and the given value
        if (dif < min_dif) then
            closest_index = ib                                  ! index for the value closest to the given value
            min_dif = dif                                       ! absolute value of the smallest difference between the bin values and the given value
        end if 
    end do

    ! Figure out if the min_dif bin value is above or below the given value and assign above and below 
    if ( (variable2bins(closest_index)-variable2) < 0 ) then
        below = variable2bins(closest_index)
        above = variable2bins(closest_index+1)
    else if ( (variable2bins(closest_index)-variable2) == 0 ) then    ! same value no interpolation needed 
        below = variable2bins(closest_index)
        above = variable2bins(closest_index)
    else ! (variable2bins(closest_index)-variable2) > 0 then 
        above = variable2bins(closest_index)
        below = variable2bins(closest_index-1)
    end if 
    
    ! Interpolate
    if  above == below .and. left == right then
        yint = table(below, left)
    elseif above == below .and. left /= right then
        yint_varY1 = table(below,left)
        yint_varY2 = table(below,right)
        y1 = yint_varY1
        x1 = variable2bins(left)
        y2 = yint_varY2
        x2 = variable2bins(right)
        xint = variable2       
        yint = y1-(((y1-y2)/(x1-x2))*(x1-xint))
    elseif above /= below .and. left == right then
        y1 = table(below,left)
        x1 = variable1bins(below)
        y2 = table(above,left)
        x2 = variable1bins(above)
        xint = variable1
        yint = y1-(((y1-y2)/(x1-x2))*(x1-xint))
    else: 
        y1 = table(below,left)
        x1 = variable1bins(below)
        y2 = table(above,left)
        x2 = variable1bins(above)
        xint = variable2
        yint_varY1 = y1-(((y1-y2)/(x1-x2))*(x1-xint))
        
        y1 = table(below,right)
        x1 = variable1bins(below)
        y2 = table(above,right)
        x2 = variable1bins(above)
        xint = variable2
        yint_varY2 = y1-(((y1-y2)/(x1-x2))*(x1-xint))
        
        y1 = yint_varY1
        x1 = variable2bins(left)
        y2 = yint_varY2
        x2 = variable2bins(right)
        xint = variable1
        yint = y1-(((y1-y2)/(x1-x2))*(x1-xint))
    end if

end