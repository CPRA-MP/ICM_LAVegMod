subroutine update_flotant(exp_lkd,total_unoccupied_flt,newly_unoccupied_thn_flt,newly_unoccupied_thk_flt)
    ! global arrays updated by subroutine:
    !   coverages
    
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

    ! dummy local variables populated with arrays passed into subroutine
    integer,dimension(ngrid,ncov),intent(in) :: exp_lkd           ! dummy variable to pass in the expansion liklihood 2D array
    integer,dimension(ngrid),intent(in) :: total_unoccupied_flt   ! dummy variable to pass in the total unoccupied flotant area of each grid cell
    integer,dimension(ngrid),intent(in) ::newly_unoccupied_thn_flt
    integer,dimension(ngrid),intent(in) ::newly_unoccupied_thk_flt

    ! local variables
    integer :: il                               ! iterator over flotant species within the thin and thick mat categories
    real(sp) :: exp_lkd_total_flt(ngrid)        ! total expansion liklihood for each grid cell 
    real(sp) :: total_flt(ngrid)                ! total area of flotant for each grid cell 
    integer :: ig                               ! iterator over the veg grid cell 

    ! Sum expansion liklihood across flotant species
    do il=1,flt_thn_cnt
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thn_indices(il))
    end do
    do il=1,flt_thk_cnt
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thk_indices(il))
    end do
    
    ! Sum the flotant in the cell; helpful to have in the next step so we can skip over cells with no flotant
    total_flt = total_unoccupied_flt 
    do il=1,flt_thn_cnt
        total_flt = total_flt + coverages(:,flt_thn_indices(il),2)
    end do
    do il=1,flt_thk_cnt
        total_flt = total_flt + coverages(:,flt_thn_indices(il),2)
    end do

    do ig=1,ngrid
        if (total_flt(ig) > 0.0) then                           ! There is flotant in the cell 
            if (exp_lkd_total_flt(ig) == 0.0) then              ! Flotant cannot establish in current conditions
                coverages(ig,dfi) = coverages(ig,dfi) + coverages(ig,bfi)             ! Dead thin mat is added to dead flotant
                coverages(ig,bfi) = 0.0                                                   ! Reset bareground flotant
                coverages(ig,dfi) = coverages(ig,dfi) + newly_unoccupied_thn_flt(ig)              ! Dead thin mat is added to dead flotant
                coverages(ig,bfi) = newly_unoccupied_thk_flt(ig)                                    ! Dead thick mat becomes bareground flotant
            else                                                ! Flotant can establish in current conditions 
                do il=1,flt_thn_cnt
                    coverages(ig,flt_thn_indices(il)) = coverages(ig,flt_thn_indices(il))+ ((exp_lkd(ig,flt_thn_indices(il))/exp_lkd_total_flt(ig))*total_unoccupied_flt)
                end do
                do il=1,flt_thk_cnt
                    coverages(ig,flt_thk_indices(il)) = coverages(ig,flt_thk_indices(il))+ ((exp_lkd(ig,flt_thk_indices(il))/exp_lkd_total_flt(ig))*total_unoccupied_flt)
                end do
            end if 
        end if
    end do


end 