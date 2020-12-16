subroutine computeTimeStep(phi, ph_lo, ph_hi, nVar, &
	&  dx, gam, dt) bind(C, name="computeTimeStep")

	use general_functions module, only : Conserv2Prim

	 implicit none
	 double precision, intent(in) :: dx(2)
     integer, intent(in) :: ph_lo(2), ph_hi(2), nVar
     double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2), 0:nVar-1)
	 double precision, intent (inout) :: dt


	 double precision :: Smax_x = 0.0d0
	 double precision :: Smax_y = 0.0d0
	 double precision :: phiCell(0:nVar-1)
	 double precision :: WCell(0:nVar-1)
	 double precision :: cS
	 double precision :: CFL = 0.9d0

	 do i = ph_lo(1), ph_hi(1)
		do j = ph_lo(2), ph_hi(2)
			do n = 0, nVar-1
				phiCell = phi(i, j, n)
			end do
				call Conserv2Prim(phiCell, gam, WCell, nVar)
				cS = DSQRT(gam*WCell(3)/WCell(0))
				Smax_x = DMAX1(Smax_x, cS + DABS(WCell(1)))
				Smax_y = DMAX1(Smax_y, cS + DABS(WCell(2)))
		end do
	end do
		dt = CFL*DMIN1(dx(1)/Smax_x, dx(2)/Smax_y)

	 end subroutine computeTimeStep
