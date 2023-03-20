      program thermo_PHDOS
      
      character*200 :: fin, fin2, fout
      integer :: nfreq, T, Tmin, Tmax, Tdelta
      real*8 :: freqnorm(5000), freqscal(5000), dos(5000)
      real*8 :: freqaux(5000), freqdos(5000), vibaux(5000), dosaux(5000)
      real*8 :: fscal, norm, de, svib, ZVE, Fvib, Utherm, cv, Ftot
      real*8 :: svib_end, ZVE_end, Fvib_end, Ftot_end
      real*8 :: Utherm_end, cv_end, norm_end
      real*8 :: c1, a1, a2, a3, zunit

      hbar=6.626068E-34
      bol=1.3806503E-23
      avo=6.0221415E+23
      c1=1.9865E-23

100   format (A200)
101   format (80('*'))
102   format (F8.2, 4(F30.8))
103   format (E18.8, E18.8, E18.8)

      write(*,*) 'Name of the input'
      read(*,100) fin
      write(*,*) 'Name of the output'
      read(*,100) fout
      
      open(5,file=fin,status='old')
      open(6,file=fout,status='new')
      
      read(5,100) fin2
      read(5,*) nfreq, fscal, de, zunit
      read(5,*) Tmin, Tmax, Tdelta

      ! Negative values filter
      open(unit=4,file=fin2,status='old')
      read(4,100)

      j=1
      do i=1, nfreq
      read(4,*) freqaux(i), dos(i)
      if (freqaux(i).gt.0) then
      freqnorm(j)=freqaux(i)
      j=j+1
      endif
      enddo
      
!      k=1
      do i=1, j-1
      freqscal(i)=freqnorm(i)*fscal
!      freqdos(i)=freqscal(i)*dos(i)
!         if (freqdos(i).ne.0) then
!         vibaux(k)=freqdos(i)
!         k=k+1
!         endif
      enddo

      do T=Tmin,Tmax,Tdelta

      Fvib=0
      norm=0
      ZVE=0
      svib=0
      Utherm=0
      cv=0
      
      do i=1,j-4,3
      a1=freqscal(i)*c1
      a2=freqscal(i+1)*c1
      a3=freqscal(i+2)*c1

      ! Defining normalization parameter
      ! Density of states for 3*NAT
      norm = norm + 3*dos(i)+3*dos(i+1)+2*dos(i+2)

      !ZVE
      ZVE = ZVE + 3*dos(i)*(0.5*a1)+    &
                3*dos(i+1)*(0.5*a2)+     &
                2*dos(i+2)*(0.5*a3)

      !Entropy
      svib=svib+3*dos(i)*(bol*((a1/(bol*T))*((exp(a1/(bol*T))-1)**(-1))-      &
                        log(1-exp((-1*a1)/(bol*T))))) +                       &
             3*dos(i+1)*(bol*((a2/(bol*T))*((exp(a2/(bol*T))-1)**(-1))-      &
                             log(1-exp((-1*a2)/(bol*T))))) +                 &
             2*dos(i+2)*(bol*((a3/(bol*T))*((exp(a3/(bol*T))-1)**(-1))-      &
                             log(1-exp((-1*a3)/(bol*T)))))

      !Fvib
      Fvib = Fvib + 3*dos(i)*bol*T*(log(1-exp((-1*a1/(bol*T)))))+           &
                   3*dos(i+1)*bol*T*(log(1-exp((-1*a2/(bol*T)))))+         &
                   2*dos(i+2)*bol*T*(log(1-exp((-1*a3/(bol*T)))))

      !Uthermal
      Utherm = Utherm + 3*dos(i)*a1*(0.5+1/((exp(a1/(bol*T)))-1))+           &
                       3*dos(i+1)*a2*(0.5+1/((exp(a2/(bol*T)))-1))+          &
                       2*dos(i+2)*a3*(0.5+1/((exp(a3/(bol*T)))-1))
     
      !Cv
      cv = Cv + 3*dos(i)*bol*(((a1/bol/T)**2)*((exp(a1/bol/T)))/              &
                       ((exp(a1/bol/T)-1)**2)) +                              &
               3*dos(i+1)*bol*(((a2/bol/T)**2)*((exp(a2/bol/T)))/            &
                       ((exp(a2/bol/T)-1)**2)) +                             &
               2*dos(i+2)*bol*(((a3/bol/T)**2)*((exp(a3/bol/T)))/            &
                       ((exp(a3/bol/T)-1)**2))
      
      enddo
      
      ZVE_end=ZVE*de*avo*3/8/zunit
      svib_end=svib*de*avo*3/8/zunit
      Fvib_end=Fvib*de*avo*3/8/zunit
      Ftot_end=Fvib_end+ZVE_end
      Utherm_end=Utherm*de*avo*3/8/zunit
      cv_end=cv*de*avo*3/8/zunit
      norm_end=norm*de*3/8
      
      if(T.eq.Tmin) then
      
      write(6,*)
      write(6,101)
      write(6,*) 'Thermodynamic parametersa from PHDOS calculations'
      write(6,101)
      write(6,*)

      write(6,'("Zero-point energy:",F14.4, " J/mol")') ZVE_end
      write(6,'("Normalization parameter (3*Nat):",F14.4)') norm_end
      write(6,*)
      write(6,101)
      write(6,'(3X,"T(K)", 4X,"E_thermal (J/mol)",8X, "F_vib (J/mol)",4X,"C_v (J/mol)",2X,"S_vib (J/mol)")')
      write(6,101)

      endif
      
      write(6,'(I10,2F21.6,2F15.4)') T, Utherm_end, Ftot_end, cv_end, svib_end

      enddo
      
      pause
      end program thermo_PHDOS
