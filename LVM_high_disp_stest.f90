subroutine high_disp_stest
    ! global arrays updated by subroutine:
    !   coverages
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   coverages
    !   bni
    !   boi
    !   cov_disp_class
    !   establish_P
    !   exp_lkd_total
    !   total_unoccupied_lnd

    ! This subroutine allows high dispersal (Class 3, weedy) vegetation to estbalish in cells with bareground after standard establishment
   
    use params
    implicit none

    ! local variables
    integer :: ig                           ! iterator over veg grid cells
    integer :: ic                           ! iterator over veg grid coverages (columns)
    real(sp) :: total_est_P                 ! sum of the establishment probabilites of all coverages in one veg grid cell; recalculated for each grid cell


    ! If exp_lkd_total is 0, then allow dispersal class 3 (weedy) species to establish 
    do ig=1,ngrid
        if (exp_lkd_total(ig) == 0.0) then
            total_est_P = 0.0
            do ic=1,ncov
                if (cov_disp_class(ic) == 3) then
                    total_est_P = total_est_P + establish_P(ig,ic)          ! sum the estalishment Ps of the weedy species 
                end if
            end do
            if (total_est_P>0) then
                do ic=1,ncov
                    if (cov_disp_class(ic) == 3) then
                        coverages(ig,ic) = coverages(ig,ic) +  ((establish_P(ig,ic)/total_est_P)*total_unoccupied_lnd(ig))         ! sum the estalishment Ps of the weedy species 
                    end if
                end do          
                coverages(ig,bni) = 0.0                                   ! reset the new bareground 
                coverages(ig,boi) = 0.0                                   ! reset the old bareground                 
            end if
        end if
    end do