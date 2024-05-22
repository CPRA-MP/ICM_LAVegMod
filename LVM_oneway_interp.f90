subroutine oneway_interp(variable1,table,variable1bins, yint)
    
    ! inputs to the subroutine:
        ! variable1 = the value of the input variable specific to the grid cell. For 2023 LAVegMod, the variable is elevation of the land portion of the grid cell. 
        ! table = establishment or mortality table for one species (1D)
        ! variable1bins = the slice of est_X_bins (or mort_X_bins, etc) for one species (1D)
    
    !output of the subroutine:
        ! yint = the interpolated value of establishment or mortality probability for the value of variable1 


    ! global arrays updated by subroutine:
    !   establish_P or mortality_P -- but not directly updated here, but rather the output updates one of those arrays

    ! global arrays used by subroutine:
    ! n_X_bins

    ! This subroutine interpolates the establishment or mortality probability based on one input. It acts on one species for one grid cell.
    ! For 2023 LAVegMod, this subroutine applies to barrier island and bottomland hardwood species only. The est and mort probabilities
    ! of those species are a function of elevation only. 

    use params
    implicit none

    ! local variables 
    integer   :: ib                               ! iterator over n_X_bins 
    real(sp)  :: above                            ! value of establisment/mortality table above the input value (variable1)
    real(sp)  :: below                            ! value of establisment/mortality table above the input value (variable1)
    real(sp)  :: min_dif                          ! the smallest difference of the differences between variable1 and each variable1bins values
    real(sp)  :: dif                              ! difference between variable1 and each variable1bins values
    integer   :: closest_index                    ! index within variable1bins for either the above or below value
    real(sp)  :: y1                               ! variable used in the linear interpolation formula 
    real(sp)  :: x1                               ! variable used in the linear interpolation formula 
    real(sp)  :: y2                               ! variable used in the linear interpolation formula 
    real(sp)  :: x2                               ! variable used in the linear interpolation formula 
    real(sp)  :: xint                             ! variable used in the linear interpolation formula 
    real(sp)  :: yint                             ! variable used in the linear interpolation formula and output for the subroutine

    ! Find the variable1 bin value closest to variable1 
    min_dif = 3000                                          ! arbitary value, just must be larger than any expected differences
    dif = 0                                                 ! initialize as 0 
    do ib = 1, n_X_bins                                     ! loop through the bin values
        dif = abs(variable1bins(ib) - variable1)            ! calculate the absolte value of the difference between each bin and the given value
        if dif < min_dif then
            closest_index = ib                              ! index for the value closest to the given value
            min_dif = dif                                   ! absolute value of the smallest difference between the bin values and the given value
        end if 
    end do

    ! Figure out if the min_dif bin value is above or below the given value and assign above and below 
    if (variable1bins(closest_index)-variable1) < 0 then
        below = variable1bins(closest_index)
        above = variable1bins(closest_index+1)
    elseif (variable1bins(closest_index)-variable1) = 0 then ! same value no interpolation needed 
        below = variable1bins(closest_index)
        above = variable1bins(closest_index)
    else ! (variable1bins(closest_index)-variable1) > 0
        above = variable1bins(closest_index)
        below = variable1bins(closest_index-1)
    end if 

    ! Interpolate 
    y1 = table(below)
    x1 = variable1bins(below)
    y2 = table(above)
    x2 = variable1bins(above)
    xint = variable1
    yint = y1-(((y1-y2)/(x1-x2))*(x-xint))
   
end 


