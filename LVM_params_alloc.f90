! This file contains two subroutines that are used to allocate memory
! for variables that are defined in PARAMS.
!
! The vast majority of variables are allocated in the PARAMS_ALLOC subroutine.
! However, some variables must be allocated for SET_IO, which is called before PARAMS_ALLOC.

subroutine params_alloc_io
    ! separate subroutine to allocate memory for variables defined in PARAMS but read in subroutine: SET_IO
    
    use params
    implicit none

    return
end

subroutine params_alloc
    ! primary parameter memory allocation subroutine used for allocated global arrays defined in PARAMS
    
    use params
    implicit none

    ! allocate memory for variables read in or clalculated in subroutine: PREPROCESSING
    allocate(grid_comp(ngrid))
    allocate(grid_x(ngrid))
    allocate(grid_y(ngrid))
    allocate(grid_a(ngrid))
    allocate(dem_pixel_proportion(ngrid))
    
    ! allocate memory for coverage attribute variables read in from input attribute table in subroutine: PREPROCESSING
    allocate(cov_symbol(ncov))
    allocate(cov_symbol_check(ncov))
    allocate(cov_grp(ncov))
    allocate(cov_disp_class(ncov))
    allocate(FFIBS(ncov))
    allocate(flt_thn_indices(ncov))
    allocate(flt_thk_indices(ncov))
    
    ! allocate memory for variables read and set in subroutine: PREPROCESSING
    ! these variables are 2D arrays [i,j] where the ith dimension represents the grid cell ID and the jth dimension represents the coverage type column
    allocate(coverages(ngrid,ncov))
    
    ! allocate memory for variables read in from compartment_out ICM-Hydro summary file in subroutine: PREPROCESSING
    allocate(stg_mx_yr(ncomp))
    allocate(stg_av_yr(ncomp))
    allocate(stg_av_smr(ncomp))
    allocate(wlv_smr(ncomp))
    allocate(sal_av_yr(ncomp))
    allocate(sal_av_smr(ncomp))
    allocate(sal_mx_14d_yr(ncomp))
    allocate(tmp_av_yr(ncomp))
    allocate(tmp_av_smr(ncomp))

    ! allocate memory for variables read in from ICM-Morph output files in subroutine: PREPROCESSING
    allocate(grid_elev(ngrid))
    allocate(water_from_morph(ngrid))
    
    ! allocate memory for variables read in in subroutine: PREPROCESSING
    allocate(barrier_island(ngrid)) 
    allocate(tree_establishment(ngrid))
    
    ! allocate memory for variables read in for establishment and mortability tables in subroutine: PREPROCESSING
    allocate(est_X_bins(n_X_bins,ncov))
    allocate(est_Y_bins(n_Y_bins,ncov))
    allocate(establish_tables(n_X_bins,n_Y_bins,ncov))
    allocate(mort_X_bins(n_X_bins,ncov))
    allocate(mort_Y_bins(n_Y_bins,ncov))
    allocate(mortality_tables(n_X_bins,n_Y_bins,ncov))
    
    ! these variables are 2D arrays [i,j] where the ith dimension represents the grid cell ID and the jth dimension represents the species coverage for the previous year [j=1] and for the current model year [j=2]
    allocate(FFIBS_score(ngrid,2))
    allocate(pct_vglnd_BLHF(ngrid,2))
    allocate(pct_vglnd_SWF(ngrid,2))
    allocate(pct_vglnd_FM(ngrid,2))
    allocate(pct_vglnd_IM(ngrid,2))
    allocate(pct_vglnd_BM(ngrid,2))
    allocate(pct_vglnd_SM(ngrid,2))

    
    ! allocate memory for variables read and set in subroutine: NEIGHBORS
    allocate(nearest_neighbors(ngrid,max_neighbors))
    allocate(near_neighbors(ngrid,max_neighbors))
    
    ! allocate memory for variables read and set in subroutine: mort_est_prob
    allocate(establish_P(ngrid,ncov))
    allocate(mortality_P(ngrid,ncov))

    return
end
