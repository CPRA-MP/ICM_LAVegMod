subroutine round_coverages
    ! global arrays updated by subroutine:
    !   coverages
    ! global arrays used by subroutine:
    !   grid_a
    ! This subroutine rounds every coverage to the nearest 1 m2
   
    use params
    implicit none

    ! local variables
    integer :: ig                           ! Iterator over the veg grid cells
    integer :: ic                           ! Iterator over the veg coverages (columns)
    real(sp)  :: coverage_area_m2           ! Value in m2 of one coverage (ic) in one grid cell(ig)
    real(sp)  :: coverage_area_flr          ! Coverage_area_m2 rounded down to the nearest whole number (i.e., a whole m2)
    real(sp)  :: dif                        ! Portion of the m2 being rounded; Difference between overage_area_m2 and coverage_area_flr


    ! Round everything to the nearest 1 m^2 and give an error if anything is negative -> note, it does not fix the negatives (with the python did do)

    do ig=1,ngrid
        do ic=1,ncov
            if (coverages(ig,ic) > 0.0) then                      ! If that cell contains that coverage, then round it. 
                coverage_area_m2 = coverages(ig,ic)*grid_a(ig)    
                
                coverages(ig,ic) = aint(coverage_area_m2)/grid_a(ig)
                
!                coverage_area_flr = floor(coverage_area_m2)
!
!                dif = coverage_area_m2 - coverage_area_flr
!                if (dif > 0.5) then   
!                    coverages(ig,ic,2) = coverage_area_flr + 1      ! round up
!                else                                                
!                    coverages(ig,ic,2) = coverage_area_flr          ! round down
!                end if 
            elseif (coverages(ig,ic) < 0.0) then
                write(*,*) '*****************WARNING************************************'
                write(*,*) 'Grid Cell ID ',ig,' has negative coverage value for',cov_symbol(ic)
                write(*,*) '************************************************************'
            end if
        end do
    end do

end 