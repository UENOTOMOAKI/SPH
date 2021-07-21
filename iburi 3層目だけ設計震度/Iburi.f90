program gen_model
   implicit none
   integer, parameter :: maxnp=10000000 ! max. num. of pacticles (limit)
   real, parameter :: dp = 5.0                                     !m
   real, parameter :: thick =15.0                                  !m
   real, parameter :: rho=1.2296e3  !rho=rho=1.6126e3              !kg/m^3
   real, parameter :: yng=53.9e3   !yng=3.865993e6                 !MN/m^2 → kN/m^2
   real, parameter :: poisn = 0.499 !0.3                           !無単位
   real, parameter :: c=4.6e3, phi=11.0  !c=1.0e4, phi=0.0         !c=kN/m^2 → N/m^2    fhi=deg
   real, parameter :: visc=0.001                                   !粘性係数  Pa・s
   real, parameter :: beta_for_stress = 100.0                      !？？
   real, parameter :: dt=1.0e-3                                   !積分時間間隔
   real, parameter :: dt_output = 0.5                              !出力時間間隔
   real, parameter :: alpha=0.01,beta=0.01                         !？？
   integer, parameter :: nstep=1.5e5                               !継続時間＝dt*nstep*2         !step
   integer :: iunit, ounit
   integer :: np,i,j,k
   integer :: nlayer 
   real, allocatable :: pos(:,:)
   integer, allocatable :: fix(:,:), mat(:)
   real :: x, y, z
   real :: xmin, xmax, ymin, ymax, zmin, zmax

   allocate(pos(3,maxnp))
   allocate(fix(3,maxnp), mat(maxnp))

   open(newunit=ounit, file='material.d', action='write', status='replace')
   write(ounit,'(i2)') 5 
   write(ounit,'(i2,7e12.4)') 1, 1.83*rho, yng*185, 0.6*poisn, 37.7*c, 0.52*phi, visc, beta_for_stress      !基盤
   write(ounit,'(i2,7e12.4)') 2, 1.83*rho, yng*185, 0.6*poisn, 37.7*c, 0.52*phi, visc, beta_for_stress      !壁
   write(ounit,'(i2,7e12.4)') 3, rho, yng, poisn, c, phi, visc, beta_for_stress          !流動粒子
   write(ounit,'(i2,7e12.4)') 4, rho, yng, poisn, c, phi, visc, beta_for_stress          !試料採取地点
   write(ounit,'(i2,7e12.4)') 5, rho, yng, poisn, 4/4.6*c, phi/11.0, visc, beta_for_stress
   close(ounit)



! 地盤モデル作成
   open(newunit=iunit,file='Iburi_Model.csv',action='read',status='old')
   np = 0
   do
      if ( np > maxnp )  then
         write(*,*) 'ERROR: np is greater than maxnp.'
      endif
      read(iunit,*,end=900) x, y, z

 ! 試料採取地点(x=32006.3 y=-138314.6)
       if(x>-32000. .and. x<-31985 .and. y>-138320. .and. y<-138310.+1.*dp ) then
         np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 4
      end if
      
      if(x>-32000. .and. x<-31985 .and. y>-138320. .and. y<-138310.+1.*dp ) then
         np=np+1
         j=2
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z -dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 4
      end if

      if(x>-32000. .and. x<-31985 .and. y>-138320. .and. y<-138310.+1.*dp ) then
         np=np+1
         j=3
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z -dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 5
      end if

  !  ①-1流動地盤粒子の読み込み_左岸側 ok

     ! 1層目(地表)
       if(x>-32015.+1.5*dp .and. x<-31720.-1.5*dp .and. y>-138305. .and. y<-138270.-2.0*dp ) then
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 2層目
        if(x>-32015.+0.5*dp .and. x<-31720.-0.5*dp .and. y>-138305. .and. y<-138270.-1.0*dp ) then
          np=np+1
          j=2
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 3層目
        if(x>-32015. .and. x<-31720. .and. y>-138305. .and. y<-138270. ) then
          np=np+1
          j=3
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 5
      end if

!  ①-2流動地盤粒子の読み込み_上流側 ok

     ! 1層目(地表)
       if(x>-32015.+1.5*dp .and. x<-31995.-1.5*dp .and. y>-138320. .and. y<-138305.-0.5*dp ) then
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 2層目
        if(x>-32015.+0.5*dp .and. x<-31995.-0.5*dp .and. y>-138320. .and. y<-138305.-0.5*dp ) then
          np=np+1
          j=2
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 3層目
        if(x>-32015. .and. x<-31995. .and. y>-138320. .and. y<-138305. ) then
          np=np+1
          j=3
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 5
      end if

!  ①-3流動地盤粒子の読み込み_下流側 ok

     ! 1層目(地表)
       if(x>-31985. .and. x<-31720.-1.5*dp .and. y>-138320. .and. y<-138305. ) then
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if
     ! 2層目
        if(x>-31985. .and. x<-31720.-0.5*dp .and. y>-138320 .and. y<-138305. ) then
          np=np+1
          j=2
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
     end if

     ! 3層目
        if(x>-31985. .and. x<-31720. .and. y>-138320. .and. y<-138305. ) then
          np=np+1
         j=3
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 5
      end if

  !  ①-4流動地盤粒子の読み込み_右岸側 ok

     ! 1層目(地表)
       if(x>-32015.+1.5*dp .and. x<-31720.-1.5*dp .and. y>-138400. .and. y<-138320. ) then
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 2層目
        if(x>-32015.+0.5*dp .and. x<-31720.-0.5*dp .and. y>-138400.-1.0*dp .and. y<-138320.) then
          np=np+1
          j=2
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 3
      end if

     ! 3層目
        if(x>-32015. .and. x<-31720. .and. y>-138400.-2.0*dp .and. y<-138320. ) then
          np=np+1
          j=3
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - dp*float(j-1)
          fix(1,np) = 0
          fix(2,np) = 0
          fix(3,np) = 0
          mat(np) = 5
      end if



  !  ②地形（固定）粒子の読み込みと生成_基盤

    ! 4～6層目
      if(x>-32040. .and. x<-31240. .and. y>-138500. .and. y<-138100. ) then
        do j=1,3
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z  -dp*2.0-1.0*dp- 1.0*dp*float(j-1)
          fix(1,np) = 1
          fix(2,np) = 1
          fix(3,np) = 1
          mat(np) = 1
        enddo
      endif
     
  !  ③囲い壁（②の範囲を参照する_dpに係数が必要な場合あり）

   ! 上流（２つ目のxを変える）
     if(x>-32040.-3.5*dp .and. x<-32040. .and. y>-138500.-1.5*dp .and. y<-138100.+0.5*dp) then
       ! 1～3層目
        do j=1,6
          np=np+1
          pos(1,np) = x
         pos(2,np) = y
         pos(3,np) = z - 1.0*dp*float(j-1)
         fix(1,np) = 1
          fix(2,np) = 1
         fix(3,np) = 1
          mat(np) = 2
        enddo
      endif

   ! 左岸（3つ目のyを変える）
     if(x>-32040.-3.5*dp .and. x<-31240. .and. y>-138100. .and. y<-138100.+3.*dp ) then
       ! 1～3層目
         do j=1,6
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - 1.0*dp*float(j-1)
          fix(1,np) = 1
          fix(2,np) = 1
          fix(3,np) = 1
          mat(np) = 2
        enddo
      endif
      
  
  ! 右岸（4つ目のyを変える）
     if(x>-32040.-3.5*dp .and. x<-31240. .and. y>-138500.-3.0*dp .and. y<-138500. ) then
       ! 1～3層目
         do j=1,6
          np=np+1
          pos(1,np) = x
          pos(2,np) = y
          pos(3,np) = z - 1.0*dp*float(j-1)
          fix(1,np) = 1
          fix(2,np) = 1
          fix(3,np) = 1
          mat(np) = 2
        enddo
      endif
  
  
   end do


   900 close(iunit)
   np = np - 1

   nlayer = int(thick/dp)
   if ( nlayer == 0 ) then
      write(*,*) 'ERROR: thick is less than dp.'
      stop
   endif



   xmin = dp*maxnp
   xmax =-xmin
   ymin = dp*maxnp
   ymax = -ymin
   zmin = dp*maxnp
   zmax = -zmin
   do i=1,np
      if ( pos(1,i) > xmax ) xmax = pos(1,i)
      if ( pos(1,i) < xmin ) xmin = pos(1,i) 
      if ( pos(2,i) > ymax ) ymax = pos(2,i)
      if ( pos(2,i) < ymin ) ymin = pos(2,i) 
      if ( pos(3,i) > zmax ) zmax = pos(3,i)
      if ( pos(3,i) < zmin ) zmin = pos(3,i) 
   end do

   open(newunit=ounit, file='particle.d', action='write', status='replace')
   write(ounit,'(i8)') np
   write(ounit,'(e16.8)') dp
   write(ounit,'(e16.8)') dp*1.3
   do i=1,np
      write(ounit,'(i8, 3e16.8,3i3,i3)') i, pos(1,i), pos(2,i), pos(3,i), &
      &  fix(1,i), fix(2,i), fix(3,i), mat(i)
   end do
   close(ounit)

   open(newunit=ounit, file='control.d', action='write', status='replace')
   write(ounit,'(a)') 'case00'
   write(ounit,'(2i3)') 0,0
   write(ounit,'(i10)') nstep
   write(ounit,'(e12.4)') dt
   write(ounit,'(e12.4)') dt_output
   write(ounit,'(i3,2e12.4)') 1,alpha,beta
   close(ounit)


   open(newunit=ounit,file='space.d', action='write', status='replace')
   write(ounit,'(2e16.8)') xmin-dp, xmax+dp
   write(ounit,'(2e16.8)') ymin-dp*50.0, ymax+dp
   write(ounit,'(2e16.8)') zmin-dp, zmax*1.5
   close(ounit)

end program gen_model
