subroutine update_flotant
    ! global arrays updated by subroutine:
    !   coverages
    !   exp_lkd_total_flt
    !   total_flt

    
    ! global arrays used by subroutine:
    !   flt_thn_cnt
    !   flt_thk_cnt
    !   flt_thn_indices
    !   flt_thk_indices
    !   ngrid
    !   ncov
    !   dfi
    !   bfi
    !   coverages

    ! This subroutine updates the coverages of the flotant types based on previous subroutines and calclations here
   
    use params
    implicit none

    ! local variables
    integer :: il                                                                                               ! iterator over flotant species within the thin and thick mat categories
    integer :: ig                                                                                               ! iterator over the veg grid cell 
    real(sp) :: flt_thn_expansion                                                                               ! grid-level variable that will store the proportional increase in thin mat flotant
    real(sp) :: flt_thk_expansion                                                                               ! grid-level variable that will store the proportional increase in thick mat flotant

    exp_lkd_total_flt = 0.0                                                                                     ! initialize array before first use
    total_flt = 0.0                                                                                             ! initialize array before first use
    

    ! Sum expansion liklihood across flotant species
    do il=1,flt_thn_cnt
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thn_indices(il))
    end do
    do il=1,flt_thk_cnt
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thk_indices(il))
    end do
    
    ! Sum the flotant in the cell (both vegetated and bare); helpful to have in the next step so we can skip over cells with no flotant
    total_flt = total_unoccupied_flt                                                                            ! At this point, total_unoccupied_flt is the same as barground flotant coverage (from end of sum_onoccupied_flt subroutine)
                                                                                                                !   - total_flt = old bare float from previous year + dead thin float from current year + dead thick float from current year                
    do il=1,flt_thn_cnt
        total_flt = total_flt + coverages(:,flt_thn_indices(il))
    end do
    do il=1,flt_thk_cnt 
        total_flt = total_flt + coverages(:,flt_thk_indices(il))
    end do

    do ig=1,ngrid
        
        flt_thn_expansion = 0.0                                                                                 ! initialize variable for current grid
        flt_thk_expansion = 0.0                                                                                 ! initialize variable for current grid
        
        if (total_flt(ig) > 0.0) then                                                                           ! There is flotant in the cell 
            old_bare_flt = coverages(:,bfi) - (newly_unoccupied_thn_flt + newly_unoccupied_thk_flt)             ! 1 - partition total bare float into portions that were:
                                                                                                                !       - bare after previous year: candidate for establishment OR for conversion to dead float
                                                                                                                !       - thin flotant that is newly bare this year: will be converted to dead float
                                                                                                                !       - thick flotant that is newly bare this year: will be converted to bare float
            if (exp_lkd_total_flt(ig) == 0.0) then                                                              ! 2 - flotant cannot establish in current conditions
                coverages(ig,dfi) = coverages(ig,dfi) + newly_unoccupied_thn_flt(ig)                            !       - bare thin mat from current year is added to dead flotant
                coverages(ig,dfi) = coverages(ig,dfi) + old_bare_flt                                            !       - old bare flotant is added to dead flotant (since nothing has established)
                coverages(ig,bfi) = newly_unoccupied_thk_flt(ig)                                                !       - bare thick mat becomes bare flotant for potential establishment in the following year
            else                                                                                                ! 3 - Flotant can establish in current conditions on any bare flotant present (old, new bare thin, new bare thick)
                do il=1,flt_thn_cnt                                                                             !       - determine how much thin flotant will establish
                    flt_thn_expansion = (exp_lkd(ig,flt_thn_indices(il))/exp_lkd_total_flt(ig))*coverages(:,bfi)!           - determine proportional area of bare_flt that will be converted to thn_flt
                    coverages(ig,flt_thn_indices(il)) = coverages(ig,flt_thn_indices(il))+ flt_thn_expansion    !           - update thin flotant with the portion of bare_flt that is newly established as thn_flt
                end do
                do il=1,flt_thk_cnt                                                                             !       - determine proportional area of bare_flt that will be converted to thk_flt
                    flt_thk_expansion = (exp_lkd(ig,flt_thk_indices(il))/exp_lkd_total_flt(ig))*coverages(:,bfi)!           - determine proportional area of bare_flt that will be converted to thk_flt
                    coverages(ig,flt_thk_indices(il)) = coverages(ig,flt_thk_indices(il))+ flt_thk_expansion    !           - update thick flotant with the portion of old_bare_flt that is newly established as thk_flt
                end do
                coverages(ig,bfi) = coverages(ig,bfi) - flt_thn_expansion - flt_thk_expansion                   ! 4 - remove the areas of think and thick flotant expansion from the bare flotant coverage - THIS SHOULD ALWAYS BE ZERO
            end if 
        end if
    end do


end 