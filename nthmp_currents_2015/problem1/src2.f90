subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: x, q2, q2inflow, xforce

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                x = xlower + (i-0.5d0)*dx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    ! MODIFIED for benchmark 1 problem, drive to hu at inflow
                    if (q(1,i,j) <= friction_depth) then
                        do nman = num_manning, 1, -1
                            if (aux(1,i,j) .lt. manning_break(nman)) then
                                coeff = manning_coefficient(nman)
                            endif
                        enddo
                        ! Calculate source term
                        if (x < xforce) then
                            gamma = 0.0055d0  ! based in inflow values, n=0.01
                            dgamma = 1.d0 + dt * gamma
                            !q2inflow = 0.054d0 * 0.115d0
                            q2inflow = q(1,i,j) * 0.115d0  ! use h here
                            q2 = (1.d0 - x/xforce) * q2inflow 
                            q(2, i, j) = (q(2, i, j) - q2) / dgamma
                            q(3, i, j) = q(3, i, j) / dgamma
                          else
                            ! usual friction formulas
                            gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                                  * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                            dgamma = 1.d0 + dt * gamma
                            q(2, i, j) = q(2, i, j) / dgamma
                            q(3, i, j) = q(3, i, j) / dgamma
                          endif
                    endif
                endif
            enddo
        enddo
    endif
    ! End of friction source term

    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

end subroutine src2
