      program NEB_extractor

!      implicit none

      character(len=256) :: fin, fin2, line, item1, test1, test2, item2
      character(len=100) :: file_name
      real*8 :: item3, item4, zero
      real*8 :: test2_array(100000), item2_array(100000)
      real*8 :: item2_array_ev(100000)
      integer :: ios, i
      logical :: first_iter = .true. ! initialize logical variable
      !a logical variable is a variable that can hold one of two values: .true. or .false.
      ! Here, the first_iter variable is a logical variable that is initially set to .true.. This variable is used to keep track of whether the first iteration of the loop has occurred
      ! Therefore, only the first iteration of write(3,*) 0, item3 inside the do loop will be made

      101   format (80('#'))

      conv=27.211399

      write(*,*) 'What is the name of the file with the NEB data?'
      read(*,*) file_name

      write(*,*) 'In which iteration the convergence was achieved?'
      read(*,*) fin
      
      write(*,*) 'TST or LST?' !Maybe it is useful for extracting data from LST calculations in CASTEP, or if you change the string so Jmol van visualise
      read(*,*) fin2

      open(2,file=file_name,status='old')
      open(3,file="output.txt",status='new')

      read(2,*)

      write(3,*)
      write(3,101)
      write(3,*) '# NEB data extractor for CASTEP v22.11 onwards'
      write(3,101)
      write(3,*) '#'
      write(3,'(1X,"#",1X,"Reaction coordinates", 11X,"Energy (eV)")')
      write(3,101)
      
      i = 0 ! initialize counter variable

      do

      read(2,*, iostat=ios) line, test1, test2
        if (ios == 0) then
      if (trim(line) /= "" .and. line(1:3) == 'REA' .and.  &
           trim(test1) == "1") then
                    read(2,*) item3
                    
      else if (trim(line) /= "" .and. line(1:3) == 'PRO' .and.  &
           trim(test1) == "1") then
                    read(2,*) item4

      else if (trim(line) /= "" .and. line(1:3) == fin2 .and.  &
           trim(test1) == fin) then
            read(2,*) item1, item2

      i = i + 1 ! increment counter variable

        read(test2,*) test2_array(i) ! add test2 value to array
        read(item2,*) item2_array(i) ! add item2 value to array
        if (first_iter) then
        write(3,'(F21.17, F26.15)') 0.0, item3*conv - item3*conv
        first_iter = .false.
        endif
        item2_array_ev(i) = item2_array(i)*conv - item3*conv
        write(3,'(F21.17, F26.15)') test2_array(i), item2_array_eV(i)

      endif
      else
      exit ! exit the loop when end-of-file is reached
      endif
      end do

       write(3,'(F21.17, F26.15)') 1.0, item4*conv - item3*conv

      ! Close the files
      close(2)
      close(3)

      ! Add the following code after closing the output file
      ! Plot the data using gnuplot
!      write(*,*) 'Plotting data using gnuplot...'
!      write(*,*) 'Press enter to continue'
!      read(*,*)

      ! Call gnuplot and pass it the appropriate commands
      call system("gnuplot -p -e     &
                ""plot 'output.txt' u 1:2 w p ps 2 pt 7 lc 7 """)


      
      ! Delete the temporary file
 !     call system('rm plot_data.tmp')

      pause
      end program NEB_extractor
