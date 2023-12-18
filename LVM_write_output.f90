subroutine write_output

    use params
    implicit none

    ! local variables
    integer :: g                        ! iterator

    ! write ICM-LAVegMod grid output file to file
    write(  *,*) ' - writing out LAVegMod grid-level output for end of current model year'
    write(000,*) ' - writing out LAVegMod grid-level output for end of current model year'


    open(unit=121, file=trim(adjustL(veg_out_file)))
    write(121,'(A)') trim(adjustL(veg_coverage_file_header))
    do g = 1,ngrid
        write(121,2345) g,                           &      ! CELLID
   &        water(g,2),                              &      ! WATER
   &        upland(g,2),                             &      ! NOTMOD
   &        bare_old(g,2),                           &      ! BAREGRND_OLD
   &        bare_new(g,2),                           &      ! BAREGRND_NEW
   &        QULA3(g,2),                              &      ! QULA3
   &        QULE(g,2),                               &      ! QULE
   &        QUNI(g,2),                               &      ! QUNI
   &        QUTE(g,2),                               &      ! QUTE
   &        QUVI(g,2),                               &      ! QUVI
   &        ULAM(g,2),                               &      ! ULAM
   &        NYAQ2(g,2),                              &      ! NYAQ2
   &        SANI(g,2),                               &      ! SANI
   &        TADI2(g,2),                              &      ! TADI2
   &        ELBA2_Flt(g,2),                          &      ! ELBA2_Flt
   &        PAHE2_Flt(g,2),                          &      ! PAHE2_Flt
   &        bare_Flt(g,2),                           &      ! BAREGRND_Flt
   &        dead_Flt(g,2),                           &      ! DEAD_Flt
   &        COES(g,2),                               &      ! COES
   &        MOCE2(g,2),                              &      ! MOCE2
   &        PAHE2(g,2),                              &      ! PAHE2
   &        SALA2(g,2),                              &      ! SALA2
   &        ZIMI(g,2),                               &      ! ZIMI
   &        CLMA10(g,2),                             &      ! CLMA10
   &        ELCE(g,2),                               &      ! ELCE
   &        IVFR(g,2),                               &      ! IVFR
   &        PAVA(g,2),                               &      ! PAVA
   &        PHAU7(g,2),                              &      ! PHAU7
   &        POPU5(g,2),                              &      ! POPU5
   &        SALA(g,2),                               &      ! SALA
   &        SCCA11(g,2),                             &      ! SCCA11
   &        TYDO(g,2),                               &      ! TYDO
   &        SCAM6(g,2),                              &      ! SCAM6
   &        SCRO5(g,2),                              &      ! SCRO5
   &        SPCY(g,2),                               &      ! SPCY
   &        SPPA(g,2),                               &      ! SPPA
   &        AVGE(g,2),                               &      ! AVGE
   &        DISP(g,2),                               &      ! DISP
   &        JURO(g,2),                               &      ! JURO
   &        SPAL(g,2),                               &      ! SPAL
   &        BAHABI(g,2),                             &      ! BAHABI
   &        DISPBI(g,2),                             &      ! DISPBI
   &        PAAM2(g,2),                              &      ! PAAM2
   &        SOSE(g,2),                               &      ! SOSE
   &        SPPABI(g,2),                             &      ! SPPABI
   &        SPVI3(g,2),                              &      ! SPVI3
   &        STHE9(g,2),                              &      ! STHE9
   &        UNPA(g,2),                               &      ! UNPA
   &        FFIBS_score(g,2),                        &      ! FFIBS
   &        pct_vglnd_BLHF(g,2),                     &      ! pL_BF
   &        pct_vglnd_SWF(g,2),                      &      ! pL_SF
   &        pct_vglnd_FM(g,2),                       &      ! pL_FM
   &        pct_vglnd_IM(g,2),                       &      ! pL_IM
   &        pct_vglnd_BM(g,2),                       &      ! pL_BM
   &        pct_vglnd_SM(g,2)                               ! pL_SM

    end do
    close(121)


2345    format(I0,54(',',F0.4))

    return
end
