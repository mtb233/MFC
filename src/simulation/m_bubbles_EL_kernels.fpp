!>
!! @file m_bubbles_EL_kernels.f90
!! @brief Contains module m_bubbles_EL_kernels
#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian bubbles in the Eulerian framework.
module m_bubbles_EL_kernels

    use m_mpi_proxy

    implicit none

contains
    pure elemental subroutine s_get_char_vol(cellx, celly, cellz, Charvol)

        integer, intent(in)   :: cellx, celly, cellz
        real(wp), intent(out) :: Charvol

        if (p > 0) then
            Charvol = dx(cellx)*dy(celly)*dz(cellz)
        else
            if (cyl_coord) then
                Charvol = dx(cellx)*dy(celly)*y_cc(celly)*2._wp*pi
            else
                Charvol = dx(cellx)*dy(celly)*lag_params%charwidth
            end if
        end if

    end subroutine s_get_char_vol

    pure subroutine s_getBubbleCell(bubble_coords, cell_coord)

        real(wp), dimension(3), intent(in) :: bubble_coords
        integer, dimension(3), intent(out) :: cell_coord
        integer                            :: i

        cell_coord(:) = int(bubble_coords(:))
        do i = 1, num_dims
            if (bubble_coords(i) < 0.0_wp) then
                cell_coord(i) = cell_coord(i) - 1
            end if
        end do

    end subroutine s_getBubbleCell

    pure subroutine s_computeKernelDeviation(cell_coord, bubble_volume, kernel_deviation)

        integer, dimension(3), intent(in) :: cell_coord
        real(wp), intent(in)              :: bubble_volume
        real(wp), intent(out)             :: kernel_deviation
        real(wp)                          :: characteristic_distance
        real(wp)                          :: bubble_radius

        if (p == 0) then
            characteristic_distance = sqrt(dx(cell_coord(1))*dy(cell_coord(2)))
        else
            characteristic_distance = (dx(cell_coord(1))*dy(cell_coord(2))*dz(cell_coord(3)))**(1.0_wp/3.0_wp)
        end if

        bubble_radius = (3.0_wp*bubble_volume/(4.0_wp*pi))**(1.0_wp/3.0_wp)

        select case (lag_params%kernel_deviation_mode)
        case (0)
            kernel_deviation = 0
        case (1)
            kernel_deviation = 1.0_wp*lag_params%epsilonb*max(characteristic_distance, bubble_radius)
        case (2)
            kernel_deviation = 1.0_wp*lag_params%epsilonb*bubble_radius
        case (3)
            kernel_deviation = 1.0_wp*lag_params%epsilonb*characteristic_distance
        end select

    end subroutine s_computeKernelDeviation

    pure subroutine s_computeKernelExtent(cell_coord, kernel_deviation, kernel_extent)

        integer, dimension(3), intent(in)  :: cell_coord
        real(wp), intent(in)               :: kernel_deviation
        integer, dimension(3), intent(out) :: kernel_extent
        real(wp)                           :: characteristic_distance

        if (lag_params%kernel_extent == 0) then
            ! Characteristic distance
            if (p == 0) then
                characteristic_distance = sqrt(dx(cell_coord(1))*dy(cell_coord(2)))
            else
                characteristic_distance = (dx(cell_coord(1))*dy(cell_coord(2))*dz(cell_coord(3)))**(1.0_wp/3.0_wp)
            end if

            ! Dynamic kernel extent
            kernel_extent(1:2) = ceiling(mapCells*kernel_deviation/characteristic_distance)
            if (p == 0) then
                kernel_extent(3) = 0
            else
                kernel_extent(3) = ceiling(mapCells*kernel_deviation/characteristic_distance)
            end if
        else
            ! Fixed kernel extent
            kernel_extent(1:2) = lag_params%kernel_extent
            if (p == 0) then
                kernel_extent(3) = 0
            else
                kernel_extent(3) = lag_params%kernel_extent
            end if
        end if

    end subroutine s_computeKernelExtent

    pure subroutine s_isCellOutside(cell_coord, cell_outside)

        integer, dimension(3), intent(inout) :: cell_coord
        logical, intent(out)                 :: cell_outside

        cell_outside = .false.
        if (num_dims == 2) then
            if ((cell_coord(1) < -buff_size) .or. (cell_coord(2) < -buff_size)) then
                cell_outside = .true.
            end if
            if (cyl_coord .and. y_cc(cell_coord(2)) < 0._wp) then
                cell_outside = .true.
            end if
            if ((cell_coord(2) > n + buff_size) .or. (cell_coord(1) > m + buff_size)) then
                cell_outside = .true.
            end if
        else
            if ((cell_coord(3) < -buff_size) .or. (cell_coord(1) < -buff_size) .or. (cell_coord(2) < -buff_size)) then
                cell_outside = .true.
            end if
            if ((cell_coord(3) > p + buff_size) .or. (cell_coord(2) > n + buff_size) .or. (cell_coord(1) > m + buff_size)) then
                cell_outside = .true.
            end if
        end if

    end subroutine s_isCellOutside

    pure subroutine s_applyGaussianShape(bubble_position, cell_position, cell_coord, kernel_deviation, kernel_value)

        real(wp), dimension(3), intent(in) :: bubble_position
        real(wp), dimension(3), intent(in) :: cell_position
        integer, dimension(3), intent(in)  :: cell_coord
        real(wp), intent(in)               :: kernel_deviation
        real(wp), intent(out)              :: kernel_value
        real(wp)                           :: cell_volume
        real(wp)                           :: distance

        if (kernel_deviation == 0.0) then
            if (p == 0) then
                cell_volume = dx(cell_coord(1))*dy(cell_coord(2))*lag_params%charwidth
            else
                cell_volume = dx(cell_coord(1))*dy(cell_coord(2))*dz(cell_coord(3))
            end if
            kernel_value = 1.0/cell_volume
        else
            distance = sqrt((bubble_position(1) - cell_position(1))**2.0_wp + (bubble_position(2) - cell_position(2))**2.0_wp &
                            & + (bubble_position(3) - cell_position(3))**2.0_wp)
            kernel_value = exp(-0.5_wp*(distance/kernel_deviation)**2.0_wp)/(sqrt(2.0_wp*pi)*kernel_deviation)**num_dims
        end if

    end subroutine s_applyGaussianShape

    pure subroutine s_shiftCellSymmetric(target_cell_coord, cell_coord, kernel_extent)

        integer, dimension(3), intent(inout) :: target_cell_coord
        integer, dimension(3), intent(in)    :: cell_coord
        integer, dimension(3), intent(in)    :: kernel_extent

        ! x-dir

        if (bc_x%beg == BC_REFLECTIVE .and. (cell_coord(1) <= kernel_extent(1) - 1)) then
            target_cell_coord(1) = abs(target_cell_coord(1)) - 1
        end if
        if (bc_x%end == BC_REFLECTIVE .and. (cell_coord(1) >= m + 1 - kernel_extent(1))) then
            target_cell_coord(1) = target_cell_coord(1) - (2*(target_cell_coord(1) - m) - 1)
        end if

        ! y-dir
        if (bc_y%beg == BC_REFLECTIVE .and. (cell_coord(2) <= kernel_extent(2) - 1)) then
            target_cell_coord(2) = abs(target_cell_coord(2)) - 1
        end if
        if (bc_y%end == BC_REFLECTIVE .and. (cell_coord(2) >= n + 1 - kernel_extent(2))) then
            target_cell_coord(2) = target_cell_coord(2) - (2*(target_cell_coord(2) - n) - 1)
        end if

        if (p > 0) then
            ! z-dir
            if (bc_z%beg == BC_REFLECTIVE .and. (cell_coord(3) <= kernel_extent(3) - 1)) then
                target_cell_coord(3) = abs(target_cell_coord(3)) - 1
            end if
            if (bc_z%end == BC_REFLECTIVE .and. (cell_coord(3) >= p + 1 - kernel_extent(3))) then
                target_cell_coord(3) = target_cell_coord(3) - (2*(target_cell_coord(3) - p) - 1)
            end if
        end if

    end subroutine s_shiftCellSymmetric

    pure subroutine s_updateVoidFractionVariables(variables, cell_coord, volume, volume_derivative, factor)

        type(scalar_field), dimension(:), intent(inout) :: variables
        integer, dimension(3), intent(in)               :: cell_coord
        real(wp), intent(in)                            :: volume
        real(wp), intent(in)                            :: volume_derivative
        real(wp), intent(in)                            :: factor

        ! Update void fraction field

        variables(1)%sf(cell_coord(1), cell_coord(2), cell_coord(3)) = variables(1)%sf(cell_coord(1), cell_coord(2), &
                  & cell_coord(3)) + real(volume*factor, kind=stp)

        ! Update time derivative of void fraction
        variables(2)%sf(cell_coord(1), cell_coord(2), cell_coord(3)) = variables(2)%sf(cell_coord(1), cell_coord(2), &
                  & cell_coord(3)) + real(volume_derivative*factor, kind=stp)

    end subroutine s_updatevoidfractionvariables

    pure subroutine s_applyKernel(bubble_count, bubble_radii, bubble_rdots, bubble_coords, bubble_positions, variables)

        integer, intent(in)                                             :: bubble_count
        real(wp), dimension(1:lag_params%nBubs_glb,1:2), intent(in)     :: bubble_radii
        real(wp), dimension(1:lag_params%nBubs_glb,1:2), intent(in)     :: bubble_rdots
        real(wp), dimension(1:lag_params%nBubs_glb,1:3,1:2), intent(in) :: bubble_coords
        real(wp), dimension(1:lag_params%nBubs_glb,1:3,1:2), intent(in) :: bubble_positions
        type(scalar_field), dimension(:), intent(inout)                 :: variables
        real(wp), dimension(3)                                          :: bubble_position
        real(wp), dimension(3)                                          :: bubble_coord
        real(wp)                                                        :: bubble_volume
        real(wp)                                                        :: bubble_vdot
        real(wp), dimension(3)                                          :: cell_position
        integer, dimension(3)                                           :: cell_coord
        real(wp)                                                        :: kernel_deviation
        integer, dimension(3)                                           :: kernel_extent
        real(wp)                                                        :: kernel_value
        real(wp), dimension(3)                                          :: target_cell_position
        integer, dimension(3)                                           :: target_cell_coord
        logical                                                         :: target_cell_outside
        integer                                                         :: l
        integer                                                         :: i
        integer                                                         :: j
        integer                                                         :: k

        do l = 1, bubble_count
            bubble_position(1:2) = bubble_positions(l,1:2,2)
            if (p > 0) then
                bubble_position(3) = bubble_positions(l, 3, 2)
            end if
            bubble_coord = bubble_coords(l,1:3,2)
            bubble_volume = 4.0_wp/3.0_wp*pi*bubble_radii(l, 2)**3.0_wp
            bubble_vdot = 4.0_wp*pi*bubble_radii(l, 2)**2.0_wp*bubble_rdots(l, 2)

            call s_getBubbleCell(bubble_coord, cell_coord)
            call s_computeKernelDeviation(cell_coord, bubble_volume, kernel_deviation)
            call s_computeKernelExtent(cell_coord, kernel_deviation, kernel_extent)

            do i = 1, (1 + 2*kernel_extent(1))
                do j = 1, (1 + 2*kernel_extent(2))
                    do k = 1, (1 + 2*kernel_extent(3))
                        target_cell_coord(1) = cell_coord(1) + i - (kernel_extent(1) + 1)
                        target_cell_coord(2) = cell_coord(2) + j - (kernel_extent(2) + 1)
                        if (p == 0) then
                            target_cell_coord(3) = 0
                        else
                            target_cell_coord(3) = cell_coord(3) + k - (kernel_extent(3) + 1)
                        end if

                        call s_isCellOutside(target_cell_coord, target_cell_outside)
                        if (.not. target_cell_outside) then
                            target_cell_position(1) = x_cc(target_cell_coord(1))
                            target_cell_position(2) = y_cc(target_cell_coord(2))
                            if (p == 0) then
                                target_cell_position(3) = 0
                            else
                                target_cell_position(3) = z_cc(target_cell_coord(3))
                            end if

                            ! Shape select
                            select case (lag_params%kernel_shape)
                            case (1)
                                call s_applyGaussianShape(bubble_position, target_cell_position, target_cell_coord, &
                                                          & kernel_deviation, kernel_value)
                            end select

                            ! Relocate cells for bubbles intersecting symmetric boundaries
                            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                                call s_shiftCellSymmetric(target_cell_coord, cell_coord, kernel_extent)
                            end if
                        else
                            kernel_value = 0.0_wp
                            target_cell_coord(1) = cell_coord(1)
                            target_cell_coord(2) = cell_coord(2)
                            if (p == 0) then
                                target_cell_coord(3) = 0
                            else
                                target_cell_coord(3) = cell_coord(3)
                            end if
                        end if
                        call s_updateVoidFractionVariables(variables, target_cell_coord, bubble_volume, bubble_vdot, kernel_value)
                    end do
                end do
            end do
        end do

    end subroutine s_applyKernel

    pure subroutine s_smoothfunction(bubble_count, bubble_radii, bubble_rdots, bubble_coords, bubble_positions, variables)

        integer, intent(in)                                             :: bubble_count
        real(wp), dimension(1:lag_params%nBubs_glb,1:2), intent(in)     :: bubble_radii
        real(wp), dimension(1:lag_params%nBubs_glb,1:2), intent(in)     :: bubble_rdots
        real(wp), dimension(1:lag_params%nBubs_glb,1:3,1:2), intent(in) :: bubble_coords
        real(wp), dimension(1:lag_params%nBubs_glb,1:3,1:2), intent(in) :: bubble_positions
        type(scalar_field), dimension(:), intent(inout)                 :: variables

        call s_applyKernel(bubble_count, bubble_radii, bubble_rdots, bubble_coords, bubble_positions, variables)

    end subroutine s_smoothfunction

end module m_bubbles_EL_kernels
