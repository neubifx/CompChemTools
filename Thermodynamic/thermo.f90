      program Thermo_castep
      
      character*200 :: fin, fin2, fout
      integer :: nfreq, T, Tmin, Tmax, Tdelta
      real*8 :: freqnorm(5000), freqscal(5000), dos(5000)
      real*8 :: freqaux(5000), freqdos(5000), vibaux(5000), dosaux(5000)
      real*8 :: fscal, svib, ZVE, Fvib, Utherm, cv, Ftot
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

      write(*,*) 'Nome do arquivo de entrada'
      read(*,100) fin
      write(*,*) 'Nome do arquivo de saida'
      read(*,100) fout
      
      open(5,file=fin,status='old')
      open(6,file=fout,status='new')
      
      read(5,100) fin2
      read(5,*) nfreq, fscal, zunit
      read(5,*) Tmin, Tmax, Tdelta

      ! Parte que desconsidera erro numerico e escalona
      open(unit=4,file=fin2,status='old')
      read(4,100)

      j=1
      do i=1, nfreq
      read(4,*) freqaux(i)
      if (freqaux(i).gt.0) then
      freqnorm(j)=freqaux(i)
      j=j+1
      endif
      enddo
      
C      k=1
      do i=1, j-1
      freqscal(i)=freqnorm(i)*fscal
C      freqdos(i)=freqscal(i)*dos(i)
C         if (freqdos(i).ne.0) then
C         vibaux(k)=freqdos(i)
C         k=k+1
C         endif
      enddo

      do T=Tmin,Tmax,Tdelta

      Fvib=0
      ZVE=0
      svib=0
      Utherm=0
      cv=0
      
      do i=1,j-1
      a1=freqscal(i)*c1

      !ZVE
      ZVE = ZVE + (0.5*a1)

      !Entropy
      svib = svib + (bol*((a1/(bol*T))*((exp(a1/(bol*T))-1)**(-1))-
     *                        log(1-exp((-1*a1)/(bol*T)))))

      !Fvib
      Fvib = Fvib + bol*T*(log(1-exp((-1*a1/(bol*T)))))

      !Uthermal
      Utherm = Utherm + a1*(0.5+1/((exp(a1/(bol*T)))-1))
     
      !Cv
      cv = Cv + bol*(((a1/bol/T)**2)*((exp(a1/bol/T)))/
     *                  ((exp(a1/bol/T)-1)**2))
      
      enddo
      
      ZVE_end=ZVE*avo*0.0000103643/zunit
      svib_end=svib*avo/zunit
      Fvib_end=Fvib*avo*0.0000103643/zunit
      Ftot_end=Fvib_end+ZVE_end
      Utherm_end=Utherm*avo*0.0000103643/zunit
      cv_end=cv*avo/zunit
      
      if(T.eq.Tmin) then
      
      write(6,*)
      write(6,101)
      write(6,*) 'Thermodynamic properties from CASTEP phonons ',
     *'calculation'
      write(6,101)
      write(6,*)

      write(6,'("Zero-point energy:",F14.4, " eV")') ZVE_end
      write(6,*)
      write(6,101)
      write(6,'(3X,"T(K)", 4X,"E_thermal (eV)",8X, "F_vib (eV)",
     *4X,"C_v (J/mol/K)",2X,"S_vib (J/mol/K)")')
      write(6,101)

      endif
      
      write(6,'(I6,F18.6,F18.6,F17.4,F17.4)') T, Utherm_end, Ftot_end,
     *cv_end, svib_end

C      write(*,*) T, Utherm_end, Ftot_end, cv_end, svib_end

      enddo
      
      pause
      end program
