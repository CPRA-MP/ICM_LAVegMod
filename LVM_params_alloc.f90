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
    allocate(water_from_morph(ngrid))

    
    ! allocate memory for variables read and set in subroutine: PREPROCESSING
    ! these variables are 2D arrays [i,j] where the ith dimension represents the grid cell ID and the jth dimension represents the species coverage for the previous year [j=1] and for the current model year [j=2]
    allocate(water(ngrid,2))
    allocate(upland(ngrid,2))
    allocate(bare_old(ngrid,2))
    allocate(bare_new(ngrid,2))
    allocate(QULA3(ngrid,2))
    allocate(QULE(ngrid,2))
    allocate(QUNI(ngrid,2))
    allocate(QUTE(ngrid,2))
    allocate(QUVI(ngrid,2))
    allocate(ULAM(ngrid,2))
    allocate(NYAQ2(ngrid,2))
    allocate(SANI(ngrid,2))
    allocate(TADI2(ngrid,2))
    allocate(ELBA2_Flt(ngrid,2))
    allocate(PAHE2_Flt(ngrid,2))
    allocate(bare_Flt(ngrid,2))
    allocate(dead_Flt(ngrid,2))
    allocate(COES(ngrid,2))
    allocate(MOCE2(ngrid,2))
    allocate(PAHE2(ngrid,2))
    allocate(SALA2(ngrid,2))
    allocate(ZIMI(ngrid,2))
    allocate(CLMA10(ngrid,2))
    allocate(ELCE(ngrid,2))
    allocate(IVFR(ngrid,2))
    allocate(PAVA(ngrid,2))
    allocate(PHAU7(ngrid,2))
    allocate(POPU5(ngrid,2))
    allocate(SALA(ngrid,2))
    allocate(SCCA11(ngrid,2))
    allocate(TYDO(ngrid,2))
    allocate(SCAM6(ngrid,2))
    allocate(SCRO5(ngrid,2))
    allocate(SPCY(ngrid,2))
    allocate(SPPA(ngrid,2))
    allocate(AVGE(ngrid,2))
    allocate(DISP(ngrid,2))
    allocate(JURO(ngrid,2))
    allocate(SPAL(ngrid,2))
    allocate(BAHABI(ngrid,2))
    allocate(DISPBI(ngrid,2))
    allocate(PAAM2(ngrid,2))
    allocate(SOSE(ngrid,2))
    allocate(SPPABI(ngrid,2))
    allocate(SPVI3(ngrid,2))
    allocate(STHE9(ngrid,2))
    allocate(UNPA(ngrid,2))
    allocate(FFIBS_score(ngrid,2))
    allocate(pct_vglnd_BLHF(ngrid,2))
    allocate(pct_vglnd_SWF(ngrid,2))
    allocate(pct_vglnd_FM(ngrid,2))
    allocate(pct_vglnd_IM(ngrid,2))
    allocate(pct_vglnd_BM(ngrid,2))
    allocate(pct_vglnd_SM(ngrid,2))

    return
end
