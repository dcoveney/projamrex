module general_functions

implicit none

private 

public :: Prim2Conserv, Conserv2Prim

contains

	subroutine Prim2Conserv(W, gam, phi, nVar)
	
	integer, intent(in) :: nVar
	double precision, intent(in	) :: phi(0:nVar-1)
	double precision, intent(  out) :: W(0:nVar-1)
	double precision, intent(in) :: gam
	
	phi(0) = W(0)
	phi(1) = W(0)*W(1)
	phi(2) = W(0)*W(2)
	phi(3) = W(3)/(gam - 1.0) + 0.5*W(0)*(W(1)**2.0 + W(2)**2.0)
	
	end subroutine Prim2Conserv
	
	subroutine Conserv2Prim(phi, gam, W, nVar)
	
	integer, intent(in) :: nVar
	double precision, intent(in	) :: phi(0:nVar-1)
	double precision, intent(  out) :: W(0:nVar-1)
	double precision, intent(in) :: gam
	
	W(0) = phi(0)
	W(1) = phi(1)/phi(0)
	W(2) = phi(2)/phi(0)
	W(3) = (gam - 1.0)*(phi(3) - 0.5*(phi(1)**2.0 + phi(2)**2.0)/phi(0)
	
	end subroutine Conserv2Prim
end module general_functions
