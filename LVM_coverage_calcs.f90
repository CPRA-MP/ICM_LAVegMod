subroutine coverage_calcs
    ! global arrays updated by subroutine:
        !FFIBS_score              ! weighted FFIBS score of ICM-LAVegMod grid cell - used for accretion
        !pct_vglnd_BLHF           ! percent of vegetated land that is bottomland hardwood forest
        !pct_vglnd_SWF            ! percent of vegetated land that is swamp forest
        !pct_vglnd_FM             ! percent of vegetated land that is fresh (attached) marsh
        !pct_vglnd_IM             ! percent of vegetated land that is intermediate marsh
        !pct_vglnd_BM             ! percent of vegetated land that is brackish marsh
        !pct_vglnd_SM             ! percent of vegetated land that is saline marsh

    ! global arrays used by subroutine:
        !FFIBS
        !cov_grp
        !ncov
        !ngrid
        !coverages
        
    ! This subroutine calculates the above summary values for the grid cell. These values are used by other models in the ICM.
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    integer :: ic                       ! iterator over veg grid coverages (columns) 
    real(sp) ::    FFIBS_times_Area     ! Running total per cell of the product of FFIBS value and coverage area
    real(sp) ::    total_vegetated_lnd  ! Running total per cell of the total vegetated coverage (land)
    real(sp) ::    total_BLHF           ! Running total per cell of the bottomland hardwood forest
    real(sp) ::    total_SWF            ! Running total per cell of the swamp forest
    real(sp) ::    total_FM             ! Running total per cell of the fresh (attached) marsh
    real(sp) ::    total_IM             ! Running total per cell of the intermediate marsh
    real(sp) ::    total_BM             ! Running total per cell of the brackish marsh
    real(sp) ::    total_SM             ! Running total per cell of the saline marsh

    do ig=1,ngrid
        FFIBS_times_Area    = 0.0
        total_vegetated_lnd = 0.0
        total_BLHF          = 0.0
        total_SWF           = 0.0
        total_FM            = 0.0
        total_IM            = 0.0
        total_BM            = 0.0
        total_SM            = 0.0
        do ic=1,ncov
            if (cov_grp(ic) == 8) then                                              !       cov_grp =  8; bottomland hardwood forest
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_BLHF = total_BLHF + coverages(ig,ic)
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic))
            elseif (cov_grp(ic) == 9) then                                          !       cov_grp =  9; swamp forest
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_SWF = total_SWF + coverages(ig,ic)
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic))
            elseif (cov_grp(ic) == 10) then                                         !       cov_grp = 10; fresh emergent wetland vegetation
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_FM = total_FM + coverages(ig,ic)                                 
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic))              
            elseif (cov_grp(ic) == 11) then                                         !       cov_grp = 11; intermediate emergent wetland vegetation
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_IM = total_IM + coverages(ig,ic) 
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic)) 
            elseif (cov_grp(ic) == 12) then                                         !       cov_grp = 12; brackish emergent wetland vegetation
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_BM = total_BM + coverages(ig,ic)  
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic))
            elseif (cov_grp(ic) == 13) then ! COME BACK ONCE WE KNOW IF BARRIER ISLAND SPECIES SHOULD BE INCLUDED HERE     !       cov_grp = 13; saline emergent wetland vegetation
                total_vegetated_lnd = total_vegetated_lnd + coverages(ig,ic)
                total_SM = total_SM + coverages(ig,ic)  
                FFIBS_times_Area = FFIBS_times_Area + (FFIBS(ic)*coverages(ig,ic))
            end if
        end do
        if (total_vegetated_lnd > 0) then
            FFIBS_score(ig)    =  FFIBS_times_Area/total_vegetated_lnd  ! weighted FFIBS score of ICM-LAVegMod grid cell - used for accretion
            pct_vglnd_BLHF(ig) =  total_BLHF/total_vegetated_lnd        ! percent of vegetated land that is bottomland hardwood forest
            pct_vglnd_SWF(ig)  =  total_SWF /total_vegetated_lnd        ! percent of vegetated land that is swamp forest
            pct_vglnd_FM(ig)   =  total_FM  /total_vegetated_lnd        ! percent of vegetated land that is fresh (attached) marsh
            pct_vglnd_IM(ig)   =  total_IM  /total_vegetated_lnd        ! percent of vegetated land that is intermediate marsh
            pct_vglnd_BM(ig)   =  total_BM  /total_vegetated_lnd        ! percent of vegetated land that is brackish marsh
            pct_vglnd_SM(ig)   =  total_SM  /total_vegetated_lnd        ! percent of vegetated land that is saline marsh
        else 
            FFIBS_score(ig)    = -9999      ! If there is no vegetated land, there should be no weighted FFIBS score
            pct_vglnd_BLHF(ig) = 0.0        ! percent of vegetated land that is bottomland hardwood forest
            pct_vglnd_SWF(ig)  = 0.0        ! percent of vegetated land that is swamp forest
            pct_vglnd_FM(ig)   = 0.0        ! percent of vegetated land that is fresh (attached) marsh
            pct_vglnd_IM(ig)   = 0.0        ! percent of vegetated land that is intermediate marsh
            pct_vglnd_BM(ig)   = 0.0        ! percent of vegetated land that is brackish marsh
            pct_vglnd_SM(ig)   = 0.0        ! percent of vegetated land that is saline marsh
        end if

    end do

end                                                                  
                                     
!       cov_grp = 14; barrier island vegetation

         
