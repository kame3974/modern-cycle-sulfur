module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0

    dn(:,patmo_idx_CS2E) = &
        - krate(:,41)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        - krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,49)*n(:,patmo_idx_CS2) &
        + krate(:,94)*n(:,patmo_idx_CS2) &
        + krate(:,95)*n(:,patmo_idx_CS2) &
        + krate(:,96)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_COS) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        + krate(:,8)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,45)*n(:,patmo_idx_COS) &
        + krate(:,54)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,55)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,57)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,60)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        - krate(:,61)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_S2) = &
        + krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,58)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,97)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_HO2) = 0d0
 !       - krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
 !       - krate(:,25)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
 !       + krate(:,29)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
 !       + krate(:,30)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
 !       + krate(:,67)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
 !       + krate(:,78)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
 !       - krate(:,82)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
 !       - krate(:,83)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3)

    dn(:,patmo_idx_NO) = 0d0
!        + krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
!        - krate(:,74)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_N) = 0d0
!        + krate(:,36)*n(:,patmo_idx_N2) &
!        + krate(:,36)*n(:,patmo_idx_N2) &
!        - krate(:,89)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
!        - krate(:,89)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_HSO) = &
        + krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,27)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,28)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,67)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,70)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,81)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CO2) = 0d0
!        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
!        - krate(:,54)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO3) = &
        + krate(:,25)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,30)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,31)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,33)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,51)*n(:,patmo_idx_SO3) &
        - krate(:,78)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,79)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,83)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,84)*n(:,patmo_idx_SO3) &
        + krate(:,86)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_SCSOH) = &
        + krate(:,6)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,59)*n(:,patmo_idx_SCSOH) &
        + krate(:,60)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2)

    dn(:,patmo_idx_H2O) = 0d0
!        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
!        + krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
!        - krate(:,33)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
!        - krate(:,64)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
!        - krate(:,67)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
!        + krate(:,86)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_HSO2) = &
        + krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,29)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,60)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        + krate(:,82)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CO) = 0d0
!        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
!        + krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
!        + krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
!        + krate(:,45)*n(:,patmo_idx_COS) &
!        - krate(:,55)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
!        - krate(:,58)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
!        - krate(:,63)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_O2) = 0d0
!        - krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
!        - krate(:,8)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
!        + krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
!        - krate(:,16)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
!        + krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
!        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
!        - krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
!        - krate(:,22)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
!        + krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
!        + krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
!        - krate(:,27)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
!        + krate(:,28)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
!        + krate(:,28)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
!        - krate(:,29)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
!        - krate(:,30)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
!        - krate(:,35)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
!        - krate(:,41)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
!        - krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
!        + krate(:,46)*n(:,patmo_idx_O3) &
!        - krate(:,47)*n(:,patmo_idx_O2) &
!        + krate(:,60)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
!        + krate(:,61)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
!        - krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
!        + krate(:,69)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
!        - krate(:,70)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
!        - krate(:,71)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
!        + krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
!        + krate(:,75)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
!        - krate(:,76)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
!        - krate(:,79)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
!        + krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
!        - krate(:,81)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
!        - krate(:,81)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
!        + krate(:,82)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
!        + krate(:,83)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
!        + krate(:,88)*n(:,patmo_idx_O3) &
!        + krate(:,94)*n(:,patmo_idx_CS2) &
!        + krate(:,96)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_N2) = 0d0
!        - krate(:,36)*n(:,patmo_idx_N2) &
!        - krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
!        + krate(:,89)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
!        + krate(:,95)*n(:,patmo_idx_CS2)

    dn(:,patmo_idx_CS2) = &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,6)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,41)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        - krate(:,48)*n(:,patmo_idx_CS2) &
        - krate(:,49)*n(:,patmo_idx_CS2) &
        + krate(:,56)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,57)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,58)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,59)*n(:,patmo_idx_SCSOH) &
        - krate(:,94)*n(:,patmo_idx_CS2) &
        - krate(:,95)*n(:,patmo_idx_CS2)

    dn(:,patmo_idx_SO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,15)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,16)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,22)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        + krate(:,50)*n(:,patmo_idx_SO2) &
        - krate(:,53)*n(:,patmo_idx_SO) &
        - krate(:,55)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,56)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,68)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,69)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,71)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,74)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,75)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,76)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,77)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,97)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_OH) = 0d0
!        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
!        - krate(:,6)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
!        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
!        + krate(:,12)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
!        + krate(:,16)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
!        - krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
!        - krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
!        + krate(:,25)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
!        + krate(:,27)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
!        - krate(:,32)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
!        + krate(:,34)*n(:,patmo_idx_H2SO4) &
!        + krate(:,34)*n(:,patmo_idx_H2SO4) &
!        - krate(:,39)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
!        - krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
!        + krate(:,54)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
!        + krate(:,59)*n(:,patmo_idx_SCSOH) &
!        + krate(:,64)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
!        - krate(:,65)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
!        - krate(:,69)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
!        + krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
!        + krate(:,77)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
!        - krate(:,78)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
!        - krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
!        + krate(:,85)*n(:,patmo_idx_HSO3) &
!        - krate(:,87)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
!        - krate(:,87)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
!        + krate(:,92)*n(:,patmo_idx_SO2) &
!        + krate(:,93)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_O) = 0d0
!        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
!        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
!        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
!        - krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
!        + krate(:,8)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
!        - krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
!        - krate(:,12)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
!        - krate(:,15)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
!        + krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
!        + krate(:,22)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
!        - krate(:,31)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
!        - krate(:,35)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
!        - krate(:,38)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
!        - krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
!        + krate(:,46)*n(:,patmo_idx_O3) &
!        + krate(:,47)*n(:,patmo_idx_O2) &
!        + krate(:,47)*n(:,patmo_idx_O2) &
!        + krate(:,50)*n(:,patmo_idx_SO2) &
!        + krate(:,51)*n(:,patmo_idx_SO3) &
!        + krate(:,53)*n(:,patmo_idx_SO) &
!        + krate(:,55)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
!        + krate(:,56)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
!        + krate(:,57)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
!        + krate(:,58)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
!        - krate(:,61)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
!        + krate(:,63)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
!        + krate(:,65)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
!        + krate(:,68)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
!        - krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
!       - krate(:,75)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
!        + krate(:,84)*n(:,patmo_idx_SO3) &
!        + krate(:,88)*n(:,patmo_idx_O3) &
!        + krate(:,91)*n(:,patmo_idx_SO2) &
!        + krate(:,97)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,33)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,34)*n(:,patmo_idx_H2SO4) &
        - krate(:,86)*n(:,patmo_idx_H2SO4) &
        + krate(:,87)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_NO2) = 0d0
!        - krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
!        + krate(:,74)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_SO4) = &
        + krate(:,37)*n(:,patmo_idx_SO2) &
        - krate(:,90)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_S) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,22)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        + krate(:,45)*n(:,patmo_idx_COS) &
        + krate(:,48)*n(:,patmo_idx_CS2) &
        + krate(:,53)*n(:,patmo_idx_SO) &
        - krate(:,57)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,63)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,75)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,76)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,77)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,97)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,38)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,39)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,91)*n(:,patmo_idx_SO2) &
        + krate(:,92)*n(:,patmo_idx_SO2) &
        + krate(:,93)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_SO2) = &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,25)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,27)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,31)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,32)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,34)*n(:,patmo_idx_H2SO4) &
        - krate(:,37)*n(:,patmo_idx_SO2) &
        + krate(:,38)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,39)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,50)*n(:,patmo_idx_SO2) &
        + krate(:,51)*n(:,patmo_idx_SO3) &
        - krate(:,71)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,74)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,78)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,79)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,82)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,84)*n(:,patmo_idx_SO3) &
        + krate(:,85)*n(:,patmo_idx_HSO3) &
        - krate(:,87)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,90)*n(:,patmo_idx_SO4) &
        - krate(:,91)*n(:,patmo_idx_SO2) &
        - krate(:,92)*n(:,patmo_idx_SO2) &
        - krate(:,93)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,96)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CH4O3S) = 0d0
!        + krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
!        - krate(:,93)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,30)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,32)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,83)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,85)*n(:,patmo_idx_HSO3)

    dn(:,patmo_idx_H2) = 0d0
  !      + krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
  !      - krate(:,66)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_H2S) = &
        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,12)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,52)*n(:,patmo_idx_H2S) &
        + krate(:,64)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,65)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,66)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,67)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)

    dn(:,patmo_idx_SH) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,12)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,15)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,16)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,28)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,52)*n(:,patmo_idx_H2S) &
        - krate(:,54)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,64)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,65)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,66)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,68)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,69)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,70)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,81)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CS) = &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,8)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        + krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,48)*n(:,patmo_idx_CS2) &
        - krate(:,56)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,61)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,63)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        - krate(:,96)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_H) = 0d0
!        - krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
!        + krate(:,15)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
!        + krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
!        + krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
!        + krate(:,52)*n(:,patmo_idx_H2S) &
!        + krate(:,66)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
!        - krate(:,68)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
!        - krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
!        - krate(:,77)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_O3) = 0d0
!        - krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
!        - krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
!        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
!        - krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
!        - krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
!        - krate(:,28)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
!        + krate(:,35)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
!        - krate(:,46)*n(:,patmo_idx_O3) &
!        + krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
!        + krate(:,70)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
!        + krate(:,71)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
!        + krate(:,76)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
!        + krate(:,79)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
!        + krate(:,81)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
!        - krate(:,88)*n(:,patmo_idx_O3)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

  !  do i=1,chemSpeciesNumber
  !    dn(:,i) = dn(:,i) &
  !        + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
  !        - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
  !        + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
  !        + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
  !  end do

!Notice that transport of sulfur compounds only 
     !S2
     do j=1,60
      dn(j,2) = dn(j,2) &
          + (k_hp(j)-d_hp(j,2)) * ngas_hpp(j) * n_p(j,2) &
          - ((k_hp(j)+d_hp(j,2)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,2)) * ngas_hmz(j)) * n(j,2) &
          + (k_hm(j)+d_hm(j,2)) * ngas_hmm(j) * n_m(j,2)
     end do
	
    do i=9,12 !HSO2, HSO3, CS2, CH4O3S
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
  	do i=17,20 !SO3, SO2, SO4, CS
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
      !SCSOH
     do j=1,60
      dn(j,i) = dn(j,23) &
          + (k_hp(j)-d_hp(j,23)) * ngas_hpp(j) * n_p(j,23) &
          - ((k_hp(j)+d_hp(j,23)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,23)) * ngas_hmz(j)) * n(j,23) &
          + (k_hm(j)+d_hm(j,23)) * ngas_hmm(j) * n_m(j,23)
     end do
	

 	do i=26,28 !H2SO4, SO3, H2S
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
  	do i=30,33 !SH, SO, CS2E, CH3SCH3
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
		 
    do i=1,60 !COS
     dn(i,1) = dn(i,1) &
          + (k_hp(i)-d_hp(i,1)) * ngas_hpp(i) * n_p(i,1) &
          - ((k_hp(i)+d_hp(i,1)) * ngas_hpz(i) &
          + (k_hm(i)-d_hm(i,1)) * ngas_hmz(i)) * n(i,1) &
         + (k_hm(i)+d_hm(i,1)) * ngas_hmm(i) * n_m(i,1)
  	end do	
	
    do i=1,60 !HSO
     dn(i,7) = dn(i,7) &
          + (k_hp(i)-d_hp(i,7)) * ngas_hpp(i) * n_p(i,7) &
          - ((k_hp(i)+d_hp(i,7)) * ngas_hpz(i) &
          + (k_hm(i)-d_hm(i,7)) * ngas_hmz(i)) * n(i,7) &
         + (k_hm(i)+d_hm(i,7)) * ngas_hmm(i) * n_m(i,7)
  	end do
	
!Calculation of troposphere-stratosphere calculation (not sure if really works!	
	!to stratosphere
		! write(43,*)  (k_hm(13)+d_hm(13,1)) * ngas_hmm(13) * n_m(13,1) 
		       
	!to troposphere 
	    ! write(44,*)  (k_hm(13)-d_hm(13,1)) * ngas_hmz(13)* n(13,1) 	
!!!!!End of transport section

  !emission
	dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 1.03d7/1d5           !OCS 1.3 Tg/y
 	dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 2.20d7/1d5           !CS2 124-250 Khalil
	dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 1.40d8/1d5           !H2S 903.5-1252.6
	dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 7.50d8/1d5           !SO2 7196-7715
	dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 2.50d8/1d5   !DMS

		
  !dry deposition
	dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - 1.585d-8*n(1,1)             !OCS  lifetime 2.0y
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - 1.355d-7*n(1,18)            !SO2  lifetime 5.0d
    dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) - 7.940d-6*n(1,33)    !DMS  lifetime 0.5d
	 
  !aerosol formation
   do i=13,34
      if (va(i) <= n(i,26) .and. pa(i) >= n(i,26)) then
	   dn(i,patmo_idx_H2SO4) = dn(i,patmo_idx_H2SO4)- (n(i,26)-va(i))
       dn(i,patmo_idx_SO4) = dn(i,patmo_idx_SO4)+ (n(i,26)-va(i))
      end if
   end do	
		
  !gravity settling SO4 Aerosol (JAM-Kasten-1968,r=0.3)
	 do j=60,2,-1
        dn(j,patmo_idx_SO4) = dn(j,patmo_idx_SO4)-gd(j)*n(j,patmo_idx_SO4)
	    dn(j-1,patmo_idx_SO4) = dn(j-1,patmo_idx_SO4)+gd(j)*n(j,patmo_idx_SO4)
      end do  
	   dn(1,19) = dn(1,19)-gd(1)*n(1,19)

  !wet deposition
	do j=12,2,-1
     do i=1,chemSpeciesNumber
         dn(j,i) = dn(j,i)-wetdep(j,i)*n(j,i)
	     dn(j-1,i) = dn(j-1,i)+wetdep(j,i)*n(j,i)
     end do   
    end do
    do i=1,chemSpeciesNumber
      dn(1,i) = dn(1,i)-wetdep(1,i)*n(1,i)
    end do  

 !DMS Å® SO2 (96%)
  ! do i=1,60
  !     dn(i,patmo_idx_CH3SCH3) = dn(i,patmo_idx_CH3SCH3)- (n(i,30)*0.96d0)
  !     dn(i,patmo_idx_SO2) = dn(i,patmo_idx_SO2)+ (n(i,30)*0.96d0)
  ! end do
! end if

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
