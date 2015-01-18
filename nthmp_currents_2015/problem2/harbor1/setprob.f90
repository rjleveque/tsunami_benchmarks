subroutine setprob()

    use geoclaw_module
    use topo_module
    use qinit_module
    use fixedgrids_module
    use refinement_module

    implicit none
    integer :: i
    real(kind=8) :: incident_t(780) , incident_eta(780)
    real(kind=8) :: dh_star, dmu_star 
    common /cominc/ incident_t, incident_eta, dh_star, dmu_star

    call set_geo()                    !# sets basic parameters g and coord system
    call set_refinement()             !# sets refinement control parameters
    call read_dtopo_settings()        !# specifies file with dtopo from earthquake
    call read_topo_settings()         !# specifies topography (bathymetry) files
    call set_qinit()                  !# specifies file with dh if this used instead
    call set_fixed_grids()            !# Fixed grid settings

    open(unit=111, &
         file='/Users/rjl/git/tsunami_benchmarks/nthmp_currents_2015/problem2/se.dat',&
         status='old', form='formatted')

    ! read incident wave and shift to start at t=1, change minutes to seconds
    ! also shift eta to start at sea_level
    do i=1,780
        read(111,*) incident_t(i), incident_eta(i)
        if (i>1) then
             incident_t(i) = 60.d0*(incident_t(i) - incident_t(1))
             incident_eta(i) = incident_eta(i) - incident_eta(1) + sea_level
                             
             endif
        enddo
    incident_t(1) = 0.d0  
    incident_eta(1) = sea_level
    write(6,*) 'Adjusting incident_eta to start at sea_level = ',sea_level

end subroutine setprob
