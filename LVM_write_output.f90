subroutine write_output(fileflag,filetag,write_summary)

    use params
    implicit none

    ! local variables
    integer :: g                        ! iterator
    integer :: write_summary            ! flag to write summary output file that includes FFIBS score and percent habitat type coverages; 1=write summary file; 0=do not write
    integer :: year                     ! calendar year of model run
    character*1 :: fileflag             ! flag to indicate whether file is an output file (fileflag = "O") or an intermediate file (fileflag = "N")
    character*10 :: filetag             ! text string to append to veg_out_file to indicate what point in the landscape update logcial structure the output file was written

    ! write ICM-LAVegMod grid output file to file
    write(  *,*) ' - writing out LAVegMod grid-level outputs for end of current model year'
    write(000,*) ' - writing out LAVegMod grid-level outputs for end of current model year'

    write(year,'I0') start_year + elapsed_year - 1


    open(unit=901, file=trim(adjustL(fnc_tag))//trim(adjustL(fileflag))//'_'//year//'_V_vegty'//trim(adjustL(filetag))//'.csv' )
    write(901,'(A)') trim(adjustL(veg_coverage_file_header))
    do g = 1,ngrid
        write(901,3456) g, coverages(g,:)
    end do
    close(901)
    
    if (write_summary == 1) then
        open(unit=902, file=trim(adjustL(fnc_tag))//trim(adjustL(fileflag))//'_'//year//'_V_vegsm'//trim(adjustL(filetag))//'.csv')
        write(902,'(A)') 'GridCellID,WeigtedFFIBS,pct_vglnd_BLHF,pct_vglnd_SWF,pct_vglnd_FM,pct_vglnd_IM,pct_vglnd_BM,pct_vglnd_SM'
        do g = 1,ngrid
            write(902,3457) g,                         &      ! grid cell ID
       &        FFIBS_score(g),                        &      ! FFIBS
       &        pct_vglnd_BLHF(g),                     &      ! pL_BF
       &        pct_vglnd_SWF(g),                      &      ! pL_SF
       &        pct_vglnd_FM(g),                       &      ! pL_FM
       &        pct_vglnd_IM(g),                       &      ! pL_IM
       &        pct_vglnd_BM(g),                       &      ! pL_BM
       &        pct_vglnd_SM(g)                               ! pL_SM
        end do
    end if
    
    close(902)


3456    format(I0,<ncov>(',',F0.4))
3457    format(I0,7(',',F0.4))
    return
end
