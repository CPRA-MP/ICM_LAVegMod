subroutine weedy_establishment
    ! global arrays updated by subroutine:
    !   coverages
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   coverages
    !   bni
    !   cov_disp_class
    !   establish_P

    ! This subroutine calculates the establishment and mortality probability for each species/coverage in each grid cell
   
    use params
    implicit none

    ! local variables
    integer :: ig                           ! iterator over veg grid cells
    integer :: ic                           ! iterator over veg grid coverages (columns)
    real(sp) :: total_est_P                 ! sum of the establishment probabilites of all coverages in one veg grid cell; recalculated for each grid cell
    integer :: dispersal_class              ! dispersal class of a coverage; weedy species are all Class 3

    do ig=1,ngrid
        if (coverages(ig,bni,2) > 0) then       ! if there is new bareground, allow class 3 (weedy) to establish
            total_est_P = 0.0                   ! reset the total establishment probability for each cell
            do ic=1,ncov
                dispersal_class = cov_disp_class(ic)    ! only applies to class 3 (weedy) species
                if (dispersal_class == 3) then
                    total_est_P = total_est_P + establish_P(ig,ic)          ! sum the estalishment Ps of the weedy species 
                endif
            end do
            if (total_est_P > 0.0) then
                do ic=1,ncov
                    dispersal_class = cov_disp_class(ic)
                    if (dispersal_class == 3) then          
                        coverages(ig,ic,2) = coverages(ig,ic,1) +  ((establish_P(ig,ic)/total_est_P)*coverages(ig,bni,2))         ! portion out the available new bareground based on realtive establishment P
                    endif
                end do          
                coverages(ig,bni,2) = 0.0               ! reset the new bareground 
                end do 
            end if
        else 
        endif
    end do

end 