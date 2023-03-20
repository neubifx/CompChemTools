program friction_GK

  implicit none

  integer ::  i, ios, j, k, atom_number, t, n, Temp
  character(len=1000) :: line, file_name
  character(len=2) :: atom1, atom2
  character(len=1000) :: misc1(400000), misc2(400000)
  real*8 ::  fx(400000), fy(400000), fz(400000), fr(400000), theta(400000), tsim(400000)
  real*8 ::  fx2d(400000), fy2d(400000), fxy_x_tsim(400000), fxy_y_tsim(400000)
  real*8 ::  c_f(400000), lambda(400000), lambda_SI(400000)
  real*8 ::  fxy_x, fxy_y, au_to_si, dt, kb, area


  write(*,*) 'What is the name of the .md file?'
  read(*,*) file_name

  write(*,*) 'Which atoms do you want? (maximum 2, input "Z" for one of them if you only want one)'
  read(*,*) atom1, atom2
  
  write(*,*) 'How many atoms are you considering in the molecule?'
  read(*,*) atom_number
  
  write(*,*) 'Input the area for the GK equation in m^2 (e.g. 2.563748196535856E-19)'
  read(*,*) area

  open(2,file = file_name,status='old')

  do i=1,4
   read(2,*)
  end do

  open(3,file="data_extract.txt",status='new')

  do
  read(2,'(A)', iostat=ios) line
  if (ios /= 0) exit ! end-of-file condition
     if (index(line, atom1) > 0 .and. index(line, 'F') > 0 ) then
        write(3,*) trim(line)
     else if (index(line, atom2) > 0 .and. index(line, 'F') > 0 ) then
        write(3,*) trim(line)
     end if
  end do

  ! close the input file
  close(2)
  close(3)
  
    open(3,file = "data_extract.txt",status='old')

  i = 0
  do
  i = i + 1
  read(3,*, iostat=ios) misc1(i), misc2(i), fx(i), fy(i), fz(i)
     if (ios /= 0) exit ! end-of-file condition
     fr(i) = sqrt(fx(i)**2 + fy(i)**2 + fz(i)**2)
     theta(i) = atan2(fx(i), fy(i))
     fx2d(i) = fr(i)*cos(theta(i))
     fy2d(i) = fr(i)*sin(theta(i))

!  write(*,*) fr(i), theta(i), fx2d(i), fy2d(i)
  end do

  open(4,file = "force_vector_2D.txt",status='new')
  write(4,'("#", 3 X, "Time (fs)", 10X, "F(x,y)_x (Eh/au)",10X, "F(x,y)_y (Eh/au)")')
  !write(4,80('#'))
  
  fxy_x = 0
  fxy_y = 0
  t = 0
  do j = 1, i-1, atom_number
     do k = j, j + atom_number - 1
        fxy_x = fxy_x + fx2d(k)
        fxy_y = fxy_y + fy2d(k)
  end do

  t = t + 1
  write(4,*)  t, fxy_x, fxy_y
  fxy_x = 0
  fxy_y = 0
  end do

  close(3)
  close(4)

  open(4,file = "force_vector_2D.txt",status='old')
  open(5,file = "friction_data.txt",status='new')
  write(5,'("#", 3 X, "Time (fs)", 7X, "AC function",10X, "<int[<F(t).F(0)>dt]",5X, "Friction coefficient (N s/m^3)")')
  
  read(4,*)

  ! calculate auto-correlation function
c_f(i) = 0
  do i=1, t-1
     read(4,*, iostat=ios) tsim(i), fxy_x_tsim(i), fxy_y_tsim(i)
        do j = 1, t-2
        c_f(i) = c_f(i) + fxy_x_tsim(j)*fxy_x_tsim(j+1) + fxy_y_tsim(j)*fxy_y_tsim(j+1)
        end do
        c_f(i) = c_f(i)/i

  end do

 ! integrate auto-correlation function to get friction coefficient

  dt = 1E-15      ! conversion to seconds (1 fs)
  kb =  1.380649E-23 ! m2 kg s-2 K-1
  Temp = 300  ! K
  au_to_si = (8.2387234983E-8)**2 !(Eh/au)^2 to N^2
!  area = 2.563748196535856E-19  !Lateral area
  lambda(i) = 0.0
  do i = 1, t-1
  lambda(i) = lambda(i) + c_f(i)*au_to_si*dt
  lambda_SI(i) = lambda(i) / (2.0 * area * kb*Temp)
  write(5,*)  i, c_f(i), lambda(i), lambda_SI(i)
  end do

  close(4)
  close(5)

  pause

  end program friction_GK
