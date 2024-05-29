!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                  
!       ICM-LAVegMod - Vegetative Coverage Model   
!                                                  
!                                                  
!   Fortran version of ICM-LAVegMod developed         
!   for 2029 Coastal Master Plan - LA CPRA         
!                                                  
!   original model (2012 MP): Visser and Duke-Sylvester (2012)
!   revised model (2017 MP): Visser and Duke-Sylvester (2017)
!   revised model (2023 MP): Foster-Martinez et al. (2023)
!   current model (2029 MP): Foster-Martinez and White (in prep)
!                                              
!   Questions: eric.white@la.gov, mrfoster@uno.edu
!   last update: 12/15/2023
!                                                     
!   project site: https://github.com/CPRA-MP      
!   documentation: http://coastal.la.gov/our-plan  
!                                                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
program main
    
    use params
    use sort
    
    implicit none
    
    ! local variables
    integer,dimension(8) :: dtvalues                ! variable to store date time values
    
    character*17 :: dtstrf                          ! string to hold formatted datetime
    character*19 :: dtstr                           ! string to hold formatted datetime

    
    call date_and_time(VALUES=dtvalues)             ! grab simulation start time
    write(dtstrf,8888) dtvalues(1),dtvalues(2),dtvalues(3),'_',dtvalues(5),'.',dtvalues(6),'.',dtvalues(7)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    open(unit=000, file=trim(adjustL('log.log')))!'_ICM-LAVegMod_runlog_')//dtstr//trim('.log')))
    
    write(  *,*)
    write(  *,*) '*************************************************************'
    write(  *,*) '****                                                     ****'
    write(  *,*) '****    ****   STARTING ICM-LAVegMod SIMULATION   ***    ****'
    write(  *,*) '****                                                     ****'
    write(  *,*) '*************************************************************'
    write(  *,*)
    write(  *,*) 'Started ICM-LAVegMod simulation at: ',dtstr
    write(  *,*)

    write(000,*)
    write(000,*) '*************************************************************'
    write(000,*) '****                                                     ****'
    write(000,*) '****    ****   STARTING ICM-LAVegMod SIMULATION   ***    ****'
    write(000,*) '****                                                     ****'
    write(000,*) '*************************************************************'    
    write(000,*)
    write(000,*) 'Started ICM-LAVegMod simulation at: ',dtstr
    write(000,*)

    call set_io                                     ! input/output settings - must be run BEFORE parameter allocation   
    call params_alloc
   
    call preprocessing
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),'_',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Preprocessing subroutine ended at: ',dtstr
    write(000,*) 'Preprocessing subroutine ended at: ',dtstr

    call neighbors
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),'_',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Nearest neighbors subroutine ended at: ',dtstr
    write(000,*) 'Nearest neighbors subroutine ended at: ',dtstr

    
    !call reset_coverages
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Reset Coverages subroutine ended at: ',dtstr
    write(000,*) 'Reset Coverages subroutine ended at: ',dtstr

    !call land_change
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Land Change subroutine ended at: ',dtstr
    write(000,*) 'Land Change subroutine ended at: ',dtstr

    !call high_dispersal_est
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'High Dispersal Establishment subroutine ended at: ',dtstr
    write(000,*) 'High Dispersal Establishment ended at: ',dtstr

    !call veg_mortality
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Vegetation Mortality subroutine ended at: ',dtstr
    write(000,*) 'Vegetation Mortality ended at: ',dtstr

    !call veg_establishment
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Vegetation Establishment subroutine ended at: ',dtstr
    write(000,*) 'Vegetation Establishment subroutine ended at: ',dtstr

    !call flotant_change
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Flotant Change subroutine ended at: ',dtstr
    write(000,*) 'Flotant Change subroutine ended at: ',dtstr

    !call acute_salinity
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Acute Salinity subroutine ended at: ',dtstr
    write(000,*) 'Acute Salinity subroutine ended at: ',dtstr

    !call check_sums
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Check Sums subroutine ended at: ',dtstr
    write(000,*) 'Check Sums subroutine ended at: ',dtstr

    !call coverage_calcs
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Coverage calculations subroutine ended at: ',dtstr
    write(000,*) 'Coverage calculations subroutine ended at: ',dtstr
    
    ! write final coverage file for End of Year landscape - also write summary output file
    call write_output('eoy       ',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Output subroutine ended at: ',dtstr
    write(000,*) 'Write Output subroutine ended at: ',dtstr    


    write(  *,*)
    write(  *,*) 'Ended ICM-LAVegMod simulation at: ',dtstr
    write(  *,*)
    
    write(000,*)
    write(000,*) 'Ended ICM-LAVegMod simulation at: ',dtstr
    write(000,*)
    close(000)


8888    format(I4.4,I2.2,I2.2,a,I2.2,a,I2.2,a,I2.2)
8889    format(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)
    
end program
