program read_python
    implicit none

    integer(kind=8), allocatable :: a(:,:,:)
    integer :: itot=2, jtot=3, ktot=4
    integer :: i, j, k

    allocate(a(itot, jtot, ktot))

    open(1, file="bla.bin", form="unformatted", status="unknown", action="read", access="stream")
    read(1) a

    do i=1, itot
        do j=1, jtot
            do k=1, ktot
                print*, i, j, k, a(i,j,k)
            end do
        end do
    end do

    close(1)
    deallocate(a)

end program read_python
