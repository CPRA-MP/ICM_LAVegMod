subroutine round_coverages
    ! global arrays updated by subroutine:
    !   coverages
    ! global arrays used by subroutine:
    !   grid_a
    ! This subroutine rounds every coverage to the nearest 1 m2
   
    use params
    implicit none

    ! local variables
    area_threshold
    coverage_area_m2 
    coverage_area_flr 
    dif 
    ig 
    ic

    ! Round coverages less than 1 m^2 down to 0.0

    area_threshold = 1/grid_a   ! Decimal percent of each grid cell equal to 1 m^2; Assumes grid_a is in m^2
    do ig=1,ngrid
        do ic=1,ncov
            if (coverages(ig,ic,2) < area_threshold(ig)) then
                coverages(ig,ic,2) = 0.0
            end if
        end do
    end do

end 


! OR round everything to the nearest 1 m^2

do ig=1,ngrid
    do ic=1,ncov
        if (coverages(ig,ic,2) > 0.0) then
            coverage_area_m2 = coverages(ig,ic,2)*grid_a(ig)
            coverage_area_flr = floor(coverage_area_m2)
            dif = coverage_area_m2 - coverage_area_flr
            if (dif > 0.5) then   
                coverages(ig,ic,2) = coverage_area_flr + 1      ! round up
            else                                                
                coverages(ig,ic,2) = coverage_area_flr          ! round down
            end if 
        else
            ! print an error message - why do we have a negative coverage value? 
        end if
    end do
end do