subroutine write_output

    use params
    implicit none

    ! local variables
    integer :: g                        ! iterator

    ! write ICM-LAVegMod grid output file to file
    write(  *,*) ' - writing out LAVegMod grid-level outputs for end of current model year'
    write(000,*) ' - writing out LAVegMod grid-level outputs for end of current model year'


    open(unit=901, file=trim(adjustL(veg_out_file)))
    write(901,'(A)') trim(adjustL(veg_coverage_file_header))
    do g = 1,ngrid
        write(901,3456) g, coverages(g,:,2)
    end do
    close(901)
   
    open(unit=902, file=trim(adjustL(veg_summary_file)))
    write(902,'(A)') 'GridCellID,WeigtedFFIBS,pct_vglnd_BLHF,pct_vglnd_SWF,pct_vglnd_FM,pct_vglnd_IM,pct_vglnd_BM,pct_vglnd_SM'
    do g = 1,ngrid
        write(902,3457) g,                           &      ! grid cell ID
   &        FFIBS_score(g,2),                        &      ! FFIBS
   &        pct_vglnd_BLHF(g,2),                     &      ! pL_BF
   &        pct_vglnd_SWF(g,2),                      &      ! pL_SF
   &        pct_vglnd_FM(g,2),                       &      ! pL_FM
   &        pct_vglnd_IM(g,2),                       &      ! pL_IM
   &        pct_vglnd_BM(g,2),                       &      ! pL_BM
   &        pct_vglnd_SM(g,2)                               ! pL_SM
    end do
    close(902)


3456    format(I0,<ncov>(',',F0.4))
3457    format(I0,7(',',F0.4))
    return
end
