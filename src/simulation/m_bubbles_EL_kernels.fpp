!>
!! @file m_bubbles_EL_kernels.f90
!! @brief Contains module m_bubbles_EL_kernels
#:include 'macros.fpp'

!> @brief This module contains kernel functions used to map the effect of the lagrangian bubbles in the Eulerian framework.
module m_bubbles_EL_kernels

    use m_mpi_proxy

    implicit none

contains

    subroutine s_smoothfunction(bubble_count, bubble_radii, bubble_rdots, bubble_coords, bubble_positions, variables)

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
        real(wp)                                                        :: correction_factor
        real(wp)                                                        :: blend_weight
        real(wp)                                                        :: flat_value

        do l = 1, bubble_count
            bubble_position(1:2) = bubble_positions(l,1:2,2)
            if (p > 0) then
                bubble_position(3) = bubble_positions(l, 3, 2)
            else
                bubble_position(3) = 0._wp
            end if
            bubble_coord = bubble_coords(l,1:3,2)
            bubble_volume = 4.0_wp/3.0_wp*pi*bubble_radii(l, 2)**3.0_wp
            bubble_vdot = 4.0_wp*pi*bubble_radii(l, 2)**2.0_wp*bubble_rdots(l, 2)

            call s_getBubbleCell(bubble_coord, cell_coord)
            call s_computeKernelDeviation(cell_coord, bubble_volume, kernel_deviation)
            call s_computeKernelExtent(cell_coord, bubble_volume, kernel_deviation, kernel_extent)

            ! Calculate correction factor and Gaussian/flat blend weight
            correction_factor = 1.0_wp
            blend_weight = 0.0_wp
            flat_value = 0.0_wp
            if (lag_params%weight_normalization .or. lag_params%kernel_blending) then
                call s_computeKernelWeights(cell_coord, bubble_position, bubble_volume, kernel_deviation, kernel_extent, &
                                            & correction_factor, blend_weight, flat_value)
            end if

            ! Smear bubble
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
                                call s_applyGaussianShape(bubble_position, target_cell_position, target_cell_coord, cell_coord, &
                                                          & kernel_deviation, kernel_value)
                            end select

                            ! flat_value is pre-divided by correction_factor, which the deposit below multiplies back in
                            if (lag_params%kernel_blending) then
                                kernel_value = (1.0_wp - blend_weight)*kernel_value + blend_weight*flat_value
                            end if

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

                        ! Update void fraction field
                        variables(1)%sf(target_cell_coord(1), target_cell_coord(2), &
                                  & target_cell_coord(3)) = variables(1)%sf(target_cell_coord(1), target_cell_coord(2), &
                                  & target_cell_coord(3)) + real(bubble_volume*kernel_value*correction_factor, kind=stp)

                        ! Update time derivative of void fraction
                        variables(2)%sf(target_cell_coord(1), target_cell_coord(2), &
                                  & target_cell_coord(3)) = variables(2)%sf(target_cell_coord(1), target_cell_coord(2), &
                                  & target_cell_coord(3)) + real(bubble_vdot*kernel_value*correction_factor, kind=stp)
                    end do
                end do
            end do
        end do

    end subroutine s_smoothfunction

    elemental subroutine s_get_char_vol(cellx, celly, cellz, Charvol)

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

    subroutine s_getBubbleCell(bubble_coords, cell_coord)

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

    subroutine s_computeKernelDeviation(cell_coord, bubble_volume, kernel_deviation)

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

    subroutine s_computeKernelExtent(cell_coord, bubble_volume, kernel_deviation, kernel_extent)

        integer, dimension(3), intent(in)  :: cell_coord
        real(wp), intent(in)               :: bubble_volume
        real(wp), intent(in)               :: kernel_deviation
        integer, dimension(3), intent(out) :: kernel_extent
        real(wp)                           :: characteristic_distance
        real(wp)                           :: cell_volume
        integer                            :: flat_extent

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

        ! Widen the stencil so the flat part of a blended kernel can stay below kernel_void_cap
        if (lag_params%kernel_blending) then
            call s_get_char_vol(cell_coord(1), cell_coord(2), cell_coord(3), cell_volume)
            flat_extent = max(0, ceiling(((bubble_volume/(lag_params%kernel_void_cap*cell_volume))**(1.0_wp/real(num_dims, &
                              & wp)) - 1.0_wp)/2.0_wp))
            kernel_extent(1:num_dims) = max(kernel_extent(1:num_dims), flat_extent)
        end if

    end subroutine s_computeKernelExtent

    subroutine s_isCellOutside(cell_coord, cell_outside)

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

    subroutine s_applyGaussianShape(bubble_position, cell_position, cell_coord, bubble_cell_coord, kernel_deviation, kernel_value)

        real(wp), dimension(3), intent(in) :: bubble_position
        real(wp), dimension(3), intent(in) :: cell_position
        integer, dimension(3), intent(in)  :: cell_coord
        integer, dimension(3), intent(in)  :: bubble_cell_coord
        real(wp), intent(in)               :: kernel_deviation
        real(wp), intent(out)              :: kernel_value
        real(wp)                           :: cell_volume
        real(wp)                           :: distance

        if (kernel_deviation == 0.0_wp) then
            ! Delta kernel: all weight in the bubble's own cell
            kernel_value = 0.0_wp
            if (all(cell_coord == bubble_cell_coord)) then
                call s_get_char_vol(cell_coord(1), cell_coord(2), cell_coord(3), cell_volume)
                kernel_value = 1.0_wp/cell_volume
            end if
        else
            distance = sqrt((bubble_position(1) - cell_position(1))**2.0_wp + (bubble_position(2) - cell_position(2))**2.0_wp &
                            & + (bubble_position(3) - cell_position(3))**2.0_wp)
            kernel_value = exp(-0.5_wp*(distance/kernel_deviation)**2.0_wp)/((sqrt(2.0_wp*pi)*kernel_deviation)**num_dims)
        end if

    end subroutine s_applyGaussianShape

    subroutine s_shiftCellSymmetric(target_cell_coord, cell_coord, kernel_extent)

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

    !> Pre-pass over the kernel stencil. Returns the weight normalization factor and, for blended kernels, the weight of the flat
    !! part needed to keep the bubble's peak void fraction at or below kernel_void_cap.
    subroutine s_computeKernelWeights(cell_coord, bubble_position, bubble_volume, kernel_deviation, kernel_extent, &
                                      & correction_factor, blend_weight, flat_value)
        integer, dimension(3), intent(in)  :: cell_coord
        real(wp), dimension(3), intent(in) :: bubble_position
        real(wp), intent(in)               :: bubble_volume
        real(wp), intent(in)               :: kernel_deviation
        integer, dimension(3), intent(in)  :: kernel_extent
        real(wp), intent(out)              :: correction_factor
        real(wp), intent(out)              :: blend_weight
        real(wp), intent(out)              :: flat_value  !< flat kernel value, pre-divided by correction_factor
        real(wp)                           :: smeared_volume
        real(wp)                           :: stencil_volume
        real(wp)                           :: kernel_peak
        real(wp)                           :: void_peak
        real(wp)                           :: void_flat
        integer, dimension(3)              :: target_cell_coord
        logical                            :: target_cell_outside
        real(wp), dimension(3)             :: target_cell_position
        real(wp)                           :: target_cell_volume
        real(wp)                           :: kernel_value
        integer                            :: i
        integer                            :: j
        integer                            :: k

        smeared_volume = 0.0_wp
        stencil_volume = 0.0_wp
        kernel_peak = 0.0_wp
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
                            call s_applyGaussianShape(bubble_position, target_cell_position, target_cell_coord, cell_coord, &
                                                      & kernel_deviation, kernel_value)
                        end select

                        ! Relocate cells for bubbles intersecting symmetric boundaries
                        if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_REFLECTIVE)) then
                            call s_shiftCellSymmetric(target_cell_coord, cell_coord, kernel_extent)
                        end if

                        call s_get_char_vol(target_cell_coord(1), target_cell_coord(2), target_cell_coord(3), target_cell_volume)
                        smeared_volume = smeared_volume + bubble_volume*kernel_value*target_cell_volume
                        stencil_volume = stencil_volume + target_cell_volume
                        kernel_peak = max(kernel_peak, kernel_value)
                    end if
                end do
            end do
        end do

        correction_factor = 1.0_wp
        if (lag_params%weight_normalization) correction_factor = bubble_volume/smeared_volume

        blend_weight = 0.0_wp
        flat_value = 0.0_wp
        if (lag_params%kernel_blending) then
            void_peak = bubble_volume*correction_factor*kernel_peak
            void_flat = bubble_volume/stencil_volume
            if (void_peak > lag_params%kernel_void_cap) then
                if (void_flat >= lag_params%kernel_void_cap) then
                    blend_weight = 1.0_wp
                else
                    blend_weight = (void_peak - lag_params%kernel_void_cap)/(void_peak - void_flat)
                end if
            end if
            flat_value = 1.0_wp/(correction_factor*stencil_volume)
        end if

    end subroutine s_computeKernelWeights

end module m_bubbles_EL_kernels
