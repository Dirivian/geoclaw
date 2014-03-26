
subroutine fgmax_values(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
                   xlower,ylower,mask_patch,values)

    ! Given a grid q (and aux if needed), set the elements of 
    !   values(mv,i,j)  for mv=1:FG_NUM_VAL
    ! to the desired values that will be output and/or monitored on
    ! the fixed grid(s).
    !
    ! Only the elements for which mask_patch(i,j) == .true. need be set.

    ! This library routine expects FG_NUM_VAL to be 1, 2, or 5 and sets:
    !   values(1,i,j) = depth          (if FG_NUM_VAL >= 1)
    !   values(2,i,j) = speed          (if FG_NUM_VAL >= 2)
    !   values(1,i,j) = momentum       (if FG_NUM_VAL == 5)
    !   values(1,i,j) = momentum flux  (if FG_NUM_VAL == 5)
    !   values(1,i,j) = -depth         (if FG_NUM_VAL == 5)
    ! The max of -depth can be used to determin the minimum depth of water
    ! at a point over the computation, useful in harbors where ships may be
    ! grounded if the depth goes too low.


    use fgmax_module
    use geoclaw_module, only: sea_level, dry_tolerance

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower
    logical, intent(in) :: mask_patch(1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: values(FG_NUM_VAL, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    real(kind=8) :: s,hs,hss,s_dry_tol
    real(kind=8), dimension(:,:), allocatable :: u,v,h
    integer :: i,j

    allocate(h(1-mbc:mx+mbc, 1-mbc:my+mbc))
    if (FG_NUM_VAL > 1) then
        allocate(u(1-mbc:mx+mbc, 1-mbc:my+mbc), v(1-mbc:mx+mbc, 1-mbc:my+mbc))
        endif

    ! This version sets only 1 value to monitor, the value eta_tilde
    ! defined to be h+B where wet and sea_level where dry

    !if (FG_NUM_VAL .ne. 1) then
    !    write(6,*) '*** Error FG_NUM_VAL in fgmax_module is ',FG_NUM_VAL
    !    write(6,*) '*** Does not agree with number expected in fgmax_values: 1'
    !    stop
    !    endif 

    h = q(1,:,:)


    where (mask_patch .and. (h >= dry_tolerance))
        values(1,:,:) = h + aux(1,:,:)
    endwhere

    if (maxval(values(1,:,:)) > 1e10) then
        write(26,*) '+++ mx,my: ',mx,my
        do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
                if (values(1,i,j) > 1e10) then
                    write(26,*) 'i,j,values,h,aux: ',i,j, &
                            values(1,i,j),q(1,i,j),aux(1,i,j)
                    endif
                enddo
            enddo
        write(6,*) '+++ aborting, see fort.26'
        stop
        endif

    where (mask_patch .and. (h < dry_tolerance))
        values(1,:,:) = -1.d90   !# to indicate dry region
    endwhere

    if (FG_NUM_VAL == 1) then
        return
        endif

    ! compute velocities:
    s_dry_tol = dry_tolerance
    where (mask_patch .and. (h > s_dry_tol))
        u = q(2,:,:) / h
        v = q(3,:,:) / h
    elsewhere
        u = 0.
        v = 0.
    endwhere

    where (mask_patch)
        values(2,:,:) = sqrt(u**2 + v**2)             ! s,    speed
    endwhere

    if (FG_NUM_VAL == 2) then
        return
        endif

    if (FG_NUM_VAL .ne. 5) then
        write(6,*) '*** Error -- expecting FG_NUM_VAL = 1, 2, or 5'
        write(6,*) '***   in fgmax_values, found FG_NUM_VAL = ',FG_NUM_VAL
        stop
        endif

    where (mask_patch)
        values(3,:,:) = h*values(2,:,:)               ! hs,   momentum
        values(4,:,:) = values(3,:,:)*values(2,:,:)   ! hs^2, momentum flux
        values(5,:,:) = -h                            ! to find min h
    endwhere

end subroutine fgmax_values
