subroutine preprocessing

    use params
    implicit none

    ! local variables
    integer :: i                        ! iterator
    integer :: g                        ! local ICM-LAVegMod grid ID variable

    ! read ICM-LAVegMod grid output file into arrays
    write(  *,*) ' - reading in LAVegMod grid-level output from previous year'
    write(000,*) ' - reading in LAVegMod grid-level output from previous year'

    ! initialize grid data arrays to zero before reading in
    water = 0.0
    upland = 0.0
    bare_old = 0.0
    bare_new = 0.0
    QULA3 = 0.0
    QULE = 0.0
    QUNI = 0.0
    QUTE = 0.0
    QUVI = 0.0
    ULAM = 0.0
    NYAQ2 = 0.0
    SANI = 0.0
    TADI2 = 0.0
    ELBA2_Flt = 0.0
    PAHE2_Flt = 0.0
    bare_Flt = 0.0
    dead_Flt = 0.0
    COES = 0.0
    MOCE2 = 0.0
    PAHE2 = 0.0
    SALA2 = 0.0
    ZIMI = 0.0
    CLMA10 = 0.0
    ELCE = 0.0
    IVFR = 0.0
    PAVA = 0.0
    PHAU7 = 0.0
    POPU5 = 0.0
    SALA = 0.0
    SCCA11 = 0.0
    TYDO = 0.0
    SCAM6 = 0.0
    SCRO5 = 0.0
    SPCY = 0.0
    SPPA = 0.0
    AVGE = 0.0
    DISP = 0.0
    JURO = 0.0
    SPAL = 0.0
    BAHABI = 0.0
    DISPBI = 0.0
    PAAM2 = 0.0
    SOSE = 0.0
    SPPABI = 0.0
    SPVI3 = 0.0
    STHE9 = 0.0
    UNPA = 0.0
    FFIBS_score = 0.0
    pct_vglnd_BLHF = 0.0
    pct_vglnd_SWF = 0.0
    pct_vglnd_FM = 0.0
    pct_vglnd_IM = 0.0
    pct_vglnd_BM = 0.0
    pct_vglnd_SM = 0.0

    open(unit=100, file=trim(adjustL(veg_in_file)))
    read(100,1234) veg_coverage_file_header                 ! dump column header row ! format 1234 must match structure of veg_out_file column headers

    do i = 1,ngrid
        read(120,*) g,                               &      ! CELLID
   &        water(g,1),                              &      ! WATER
   &        upland(g,1),                             &      ! NOTMOD
   &        bare_old(g,1),                           &      ! BAREGRND_OLD
   &        bare_new(g,1),                           &      ! BAREGRND_NEW
   &        QULA3(g,1),                              &      ! QULA3
   &        QULE(g,1),                               &      ! QULE
   &        QUNI(g,1),                               &      ! QUNI
   &        QUTE(g,1),                               &      ! QUTE
   &        QUVI(g,1),                               &      ! QUVI
   &        ULAM(g,1),                               &      ! ULAM
   &        NYAQ2(g,1),                              &      ! NYAQ2
   &        SANI(g,1),                               &      ! SANI
   &        TADI2(g,1),                              &      ! TADI2
   &        ELBA2_Flt(g,1),                          &      ! ELBA2_Flt
   &        PAHE2_Flt(g,1),                          &      ! PAHE2_Flt
   &        bare_Flt(g,1),                           &      ! BAREGRND_Flt
   &        dead_Flt(g,1),                           &      ! DEAD_Flt
   &        COES(g,1),                               &      ! COES
   &        MOCE2(g,1),                              &      ! MOCE2
   &        PAHE2(g,1),                              &      ! PAHE2
   &        SALA2(g,1),                              &      ! SALA2
   &        ZIMI(g,1),                               &      ! ZIMI
   &        CLMA10(g,1),                             &      ! CLMA10
   &        ELCE(g,1),                               &      ! ELCE
   &        IVFR(g,1),                               &      ! IVFR
   &        PAVA(g,1),                               &      ! PAVA
   &        PHAU7(g,1),                              &      ! PHAU7
   &        POPU5(g,1),                              &      ! POPU5
   &        SALA(g,1),                               &      ! SALA
   &        SCCA11(g,1),                             &      ! SCCA11
   &        TYDO(g,1),                               &      ! TYDO
   &        SCAM6(g,1),                              &      ! SCAM6
   &        SCRO5(g,1),                              &      ! SCRO5
   &        SPCY(g,1),                               &      ! SPCY
   &        SPPA(g,1),                               &      ! SPPA
   &        AVGE(g,1),                               &      ! AVGE
   &        DISP(g,1),                               &      ! DISP
   &        JURO(g,1),                               &      ! JURO
   &        SPAL(g,1),                               &      ! SPAL
   &        BAHABI(g,1),                             &      ! BAHABI
   &        DISPBI(g,1),                             &      ! DISPBI
   &        PAAM2(g,1),                              &      ! PAAM2
   &        SOSE(g,1),                               &      ! SOSE
   &        SPPABI(g,1),                             &      ! SPPABI
   &        SPVI3(g,1),                              &      ! SPVI3
   &        STHE9(g,1),                              &      ! STHE9
   &        UNPA(g,1),                               &      ! UNPA
   &        FFIBS_score(g,1),                        &      ! FFIBS
   &        pct_vglnd_BLHF(g,1),                     &      ! pL_BF
   &        pct_vglnd_SWF(g,1),                      &      ! pL_SF
   &        pct_vglnd_FM(g,1),                       &      ! pL_FM
   &        pct_vglnd_IM(g,1),                       &      ! pL_IM
   &        pct_vglnd_BM(g,1),                       &      ! pL_BM
   &        pct_vglnd_SM(g,1)                               ! pL_SM
    end do
    close(100)

    
   ! read ICM-Hydro compartment hydro output data in from file
    write(  *,*) ' - reading in annual ICM-Hydro compartment-level output'
    write(000,*) ' - reading in annual ICM-Hydro compartment-level output'
    
    open(unit=101, file=trim(adjustL(hydro_comp_out_file)))
    
    read(101,*) dump_txt        ! dump header
    do i = 1,ncomp
        read(112,*) dump_txt,               &
   &         stg_mx_yr(i),                  &
   &         stg_av_yr(i),                  &
   &         stg_av_smr(i),                 &
   &         wlv_smr(i),                    &
   &         sal_av_yr(i),                  &
   &         sal_av_smr(i),                 &
   &         sal_mx_14d_yr(i),              &
   &         tmp_av_yr(i),                  &
   &         tmp_av_smr(i),                 &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt
    end do
    close(101)
    

1234    format(A,53(',',A))

    return
end
