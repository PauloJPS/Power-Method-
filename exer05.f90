program exer05
	implicit none

	integer :: i, j
	integer :: N, k
	
	real(8), allocatable :: Matrix(:,:)
	real(8), allocatable :: MatrixRes(:,:), MatrixAux(:,:)
	real(8), allocatable :: x(:)
	real(8), allocatable :: xK(:), xK1(:)
	real(8) :: lambdaK, lambdaK1, num, dem, eps, norm

	read(*,*) eps
	read(*,*) N
	allocate(Matrix(N,N))
	allocate(MatrixRes(N,N))
	allocate(MatrixAux(N,N))
	allocate(x(N))
	allocate(xK(N))
	allocate(xK1(N))


	do i=1, N
		read(*,*) Matrix(i,:)	
		x(i) = 1d0
	end do

	do i=1, N
		do j=1, N
			MatrixAux(i,j) = Matrix(i,j)	
		end do
	end do

	call matrixMultiply(MatrixAux, Matrix, N, MatrixRes)	
	call matrixVectorMultiply(MatrixRes, x, N, xK)
	call matrixVectorMultiply(Matrix, xk, N, xK1)

	call dotProduct(xk, xk1, N, num)
	call dotProduct(xk, xk, N, dem)

	lambdaK1 = num/dem

	lambdaK = lambdaK1
	do i=1, N
		x(i) = xk(i)
		xk(i) = xk1(i)
		do j=1, N
			MatrixAux(i,j) = MatrixRes(i,j)	
		end do
	end do


	do while( 1 > 0)
		call matrixMultiply(MatrixAux, Matrix, N, MatrixRes)	
		call matrixVectorMultiply(MatrixRes, x, N, xK)
		call matrixVectorMultiply(Matrix, xk, N, xK1)

		call dotProduct(xk, xk1, N, num)
		call dotProduct(xk, xk, N, dem)

		lambdaK1 = num/dem

		if( abs(lambdaK1 - lambdaK) <= eps) then
			exit
		endif

		lambdaK = lambdaK1
		do i=1, N
			x(i) = xk(i)
			xk(i) = xk1(i)
			do j=1, N
				MatrixAux(i,j) = MatrixRes(i,j)	
			end do
		end do
	end do

	write(*,*) lambdaK
	do i=1, N
		write(*,*) xk1(i)
	end do

	deallocate(Matrix)
	deallocate(MatrixRes)
	deallocate(MatrixAux)
	deallocate(x)
	deallocate(xK)
	deallocate(xK1)

	contains

	subroutine matrixVectorMultiply(mat, vec, N, vecRes)
		implicit none

		integer, intent(in) :: N
		real(8), dimension(N,N), intent(in) :: mat
		real(8), dimension(N), intent(in) :: vec
		real(8), dimension(N), intent(out) :: vecRes
		integer :: i, j, k
		real(8) :: acumulator

		do i=1, N
			acumulator = 0d0
			do j=1, N
				acumulator = acumulator + mat(i,j) * vec(j)
			end do
			vecRes(i) = acumulator
		end do

	end subroutine matrixVectorMultiply

	subroutine matrixMultiply(mat1, mat2, N, matRes)
		implicit none

		integer, intent(in) :: N
		real(8), dimension(N,N), intent(in) :: mat1, mat2
		real(8), dimension(N,N), intent(out) :: matRes
		integer :: i, j, k
		real(8) :: acumulator

		do i=1, N
			do j=1, N
				acumulator = 0d0
				do k=1, N
					acumulator = acumulator + mat1(i, k)*mat2(k, j)
				end do
			matRes(i, j) = acumulator	
			end do
		end do
	end subroutine matrixMultiply

	subroutine dotProduct(vec1, vec2, N, res)
		implicit none

		integer, intent(in) :: N
		real(8), dimension(N), intent(in) :: vec1
		real(8), dimension(N), intent(in) :: vec2
		real(8), intent(out) :: res
		integer :: i

		res = 0d0
		do i=1, N
			res = res + vec1(i)*vec2(i)
		end do
	end subroutine dotProduct

end program 
