module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !COS -> CO + S
    krate(:,45) = integrateXsec(1, tau(:,:))

    !O3 -> O2 + O
    krate(:,46) = integrateXsec(2, tau(:,:))

    !O2 -> O + O
    krate(:,47) = integrateXsec(3, tau(:,:))

    !CS2 -> CS + S
    krate(:,48) = integrateXsec(4, tau(:,:))

    !CS2 -> CS2E
    krate(:,49) = integrateXsec(5, tau(:,:))

    !SO2 -> SO + O
    krate(:,50) = integrateXsec(6, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,51) = integrateXsec(7, tau(:,:))

    !H2S -> SH + H
    krate(:,52) = integrateXsec(8, tau(:,:))

    !SO -> S + O
    krate(:,53) = integrateXsec(9, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber) ,dE
    integer::j

    !loop on cells (stride photobins)
!    do j=1,cellsNumber
!      integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
!          /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
!    end do
    	dE=0.04956071 !nm

  !loop on cells (stride photobins)
    do j=1,cellsNumber
	      integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:)*exp(-tau(:,j))*dE) 
    end do

  end function integrateXsec

end module patmo_photoRates
