module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,ntot(cellsNumber)
    integer::icell

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)
    open(59,file="H2SO4_photorate.txt",status="old")
    open(60,file="SO2_O_rate.txt",status="old")
    open(61,file="SO2_OH_rate.txt",status="old")
    open(62,file="CS2_OH_rate.txt",status="old")
    open(63,file="DMS_OH_rate.txt",status="old")

    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !COS + OH -> CO2 + SH
      krate(icell,1) = 1.1d-13*exp(-1200/T)

      !COS + O -> CO + SO
      krate(icell,2) = 2.1d-11*exp(-2200/T)

      !CS2 + O -> CS + SO
      krate(icell,3) = 3.30d-11*exp(-650/T)

      !CS2 + O -> COS + S
      krate(icell,4) = 2.90d-12*exp(-650/T)

      !CS2 + O -> CO + S2
      krate(icell,5) = 1.20d-12*exp(-650/T)

      !CS2 + OH -> SCSOH
      krate(icell,6) = 3.80d-13
      !read(62,*) krate(icell,6)
  
      !SCSOH + O2 -> COS + HSO2
      krate(icell,7) = 2.70d-14

      !CS + O2 -> COS + O
      krate(icell,8) = 2.9d-19

      !CS + O3 -> COS + O2
      krate(icell,9) = 3.0d-16

      !CS + O -> CO + S
      krate(icell,10) = 2.70d-10*exp(-761/T)

      !H2S + OH -> H2O + SH
      krate(icell,11) = 6.10d-12*exp(-75/T)

      !H2S + O -> OH + SH
      krate(icell,12) = 9.22d-12*exp(-1803/T)

      !H2S + H -> H2 + SH
      krate(icell,13) = 8.00d-13

      !H2S + HO2 -> H2O + HSO
      krate(icell,14) = 3.00d-15

      !SH + O -> H + SO
      krate(icell,15) = 1.60d-10
 

      !SH + O2 -> OH + SO
      krate(icell,16) = 4.00d-19

      !SH + O3 -> HSO + O2
      krate(icell,17) = 9.00d-12*exp(-280/T)

      !SO + O3 -> SO2 + O2
      krate(icell,18) = 4.50d-12*exp(-1170/T)

      !SO + O2 -> SO2 + O
      krate(icell,19) = 1.60d-13*exp(-2282/T)

      !SO + OH -> SO2 + H
      krate(icell,20) = 2.70d-11*exp(335/T)

      !SO + NO2 -> SO2 + NO
      krate(icell,21) = 0d0

      !S + O2 -> SO + O
      krate(icell,22) = 2.31d-12

      !S + O3 -> O2 + SO
      krate(icell,23) = 1.20d-11

      !S + OH -> H + SO
      krate(icell,24) = 6.59d-11

      !SO2 + HO2 -> OH + SO3
      krate(icell,25) = 1.00d-18

      !SO2 + O3 -> SO3 + O2
      krate(icell,26) = 3.00d-12*exp(-7000/T)

      !HSO + O2 -> SO2 + OH
      krate(icell,27) = 1.69d-15

      !HSO + O3 -> O2 + O2 + SH
      krate(icell,28) = 2.54d-13*exp(-392.4/T)

      !HSO2 + O2 -> HO2 + SO2
      krate(icell,29) = 3.01d-13

      !HSO3 + O2 -> HO2 + SO3
      krate(icell,30) = 1.30d-12*exp(-330/T)

      !SO2 + O -> SO3
      !krate(icell,31) = 1.80d-33*(T/300)**(2.00)
      read(60,*) krate(icell,31)

      !SO2 + OH -> HSO3
      !krate(icell,32) = 3.30d-31*(T/300)**(-4.30)
      read(61,*) krate(icell,32)	  

      !SO3 + H2O -> H2SO4
      krate(icell,33) = 1.20d-15

      !H2SO4 -> SO2 + OH + OH
      !krate(icell,34) = 1.20d-15
      read(59,*) krate(icell,34)

      !O + O2 -> O3
      krate(icell,35) = 0d0

      !N2 -> N + N
      krate(icell,36) = 0d0

      !SO2 -> SO4
      !krate(icell,37) = 1.2d-5
      	krate(icell,37) =0d0

      !CH3SCH3 + O -> SO2
      krate(icell,38) = 1.0d-11*exp(410/T)

      !CH3SCH3 + OH -> SO2
      krate(icell,39) = 1.2d-11*exp(-260/T)

      !CH3SCH3 + OH -> SO2 + CH4O3S
      !krate(icell,40) = 3.04d-12*exp(350/T)
      read(63,*) krate(icell,40)

      !CS2E + O2 -> CS2
      krate(icell,41) = 2.5d-11

      !CS2E + N2 -> CS2
      krate(icell,42) = 2.5d-11

      !CS2E + O2 -> CS + SO2
      krate(icell,43) = 1.25d-12

      !S2 + O -> S + SO
      krate(icell,44) = 1.66d-11

    end do

    close(59)
    close(60)
    close(61)
    close(62)
    close(63)

  end subroutine computeRates

end module patmo_rates
