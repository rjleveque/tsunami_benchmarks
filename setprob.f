c=========================================================================
      subroutine setprob()
c=========================================================================

      use geoclaw_module
      use topo_module
      use qinit_module
      use fixedgrids_module
      use refinement_module

      implicit double precision (a-h,o-z)
      common /wave/ profile(451,2)


      call set_geo()          !# sets basic parameters g and coord system
      call set_refinement()             !# sets refinement control parameters
      call read_dtopo_settings()        !# specifies file with dtopo from earthquake
      call read_topo_settings()         !# specifies topography (bathymetry) files
      call set_qinit()         !# specifies file with dh if this used instead
      call set_fixed_grids()    !# specifies output on arbitrary uniform fixed grids
    

      
      open(unit=76,file='../wave.txt',status='old',form='formatted')

      do it=1,451
         read(76,*) profile(it,1),profile(it,2)
      enddo
      close(unit=76)

      return
      end
