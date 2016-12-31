      module def
      implicit none
      integer nb
   ! define zone structure
	type zone_t
		double precision, pointer :: 
     $		x(:,:,:),y(:,:,:),z(:,:,:),
     $		xx(:,:,:),yy(:,:,:),zz(:,:,:),
     $      x1(:,:,:),x2(:,:,:),
     $      y1(:,:,:),y2(:,:,:),
     $      z1(:,:,:),z2(:,:,:)
      character(len=1), pointer::flag1(:,:),flag2(:,:),flag3(:,:),
     $                   flag4(:,:),flag5(:,:),flag6(:,:)
      end type zone_t
      type(zone_t), dimension(:), allocatable :: b

      contains

	subroutine allo_array(z_size)	
		implicit none
		integer, dimension(3,nb), intent(in) :: z_size
		integer i
		
		allocate(b(nb))
		
		do i=1,nb
			allocate(b(i).xx(z_size(1,i),z_size(2,i),z_size(3,i)))
			allocate(b(i).yy(z_size(1,i),z_size(2,i),z_size(3,i)))
			allocate(b(i).zz(z_size(1,i),z_size(2,i),z_size(3,i)))
			
			allocate(b(i).x(0:z_size(1,i),0:z_size(2,i),0:z_size(3,i)))
			allocate(b(i).y(0:z_size(1,i),0:z_size(2,i),0:z_size(3,i)))
			allocate(b(i).z(0:z_size(1,i),0:z_size(2,i),0:z_size(3,i)))	
			
			allocate(b(i).flag1(z_size(2,i)-1,z_size(3,i)-1))							
			allocate(b(i).flag2(z_size(2,i)-1,z_size(3,i)-1))							
			allocate(b(i).flag3(z_size(1,i)-1,z_size(3,i)-1))							
			allocate(b(i).flag4(z_size(1,i)-1,z_size(3,i)-1))							
			allocate(b(i).flag5(z_size(1,i)-1,z_size(2,i)-1))							
			allocate(b(i).flag6(z_size(1,i)-1,z_size(2,i)-1))	

			allocate(b(i).x1(3,0:z_size(2,i),0:z_size(3,i)))							
			allocate(b(i).x2(3,0:z_size(2,i),0:z_size(3,i)))							
			allocate(b(i).y1(3,0:z_size(1,i),0:z_size(3,i)))							
			allocate(b(i).y2(3,0:z_size(1,i),0:z_size(3,i)))							
			allocate(b(i).z1(3,0:z_size(1,i),0:z_size(2,i)))							
			allocate(b(i).z2(3,0:z_size(1,i),0:z_size(2,i)))				
		end do		
		return
      end subroutine allo_array
      end module def


      program connect_info 
      use def     
      implicit none
	integer, allocatable :: z_size(:,:)
      double precision, allocatable, dimension(:,:)::
     $                   xmin,xmax,ymin,ymax,zmin,zmax
      integer, parameter:: fdatain = 1
      integer, dimension(:), allocatable:: iL,jL,kL
      character(len=20):: filename
      integer::iltp,jltp,kltp,
     $          n, nf,i,j,k,
     $          mb,nfb,i1,j1,k1
     $     ,idimen
      double precision::
     $          temp1,temp2,temp3,eps
      integer::
     $          start(3),endp(3),istart(3),iendp(3),order(3),width
c------in main process     
      integer::
     $          case_11,case_22,
     $                m2,m22,ic,
     $                i3,j3,k3,
     $                i4,j4,k4
      character*20 bcdir
     

      width=8
calex
      open(2,file='pre_bc7.err')
      write(*,*)'Enter the smallest mesh size in bin files'
      write(*,*)'Enter 1d-8 or lower if you do not know'
      read(*,*) eps
      eps=eps/10.d0
      open(12,file='pre_bc7_datain')

c
      write(*,*)'please input the dimension (2/3):'
      read(*,*) idimen
     
      open(1, fiLe='init.input')
      read(1, *)
      read(1, *) nb
      allocate(z_size(3,nb))      
      allocate (iL(nb),jL(nb),kL(nb))
      allocate (xmin(nb,6),xmax(nb,6),ymin(nb,6),ymax(nb,6),
     $           zmin(nb,6),zmax(nb,6))      
      do i=1,nb
         read(1, *) iL(i),jL(i),kL(i)
         if(idimen .eq. 2) kL(i) = 2
         z_size(1,i)=il(i)+1
         z_size(2,i)=jl(i)+1
         z_size(3,i)=kl(i)+1
      end do
      close(1)
      
      call allo_array(z_size)      

      do n=1,nb
           b(n).x1(:,:,:)=1.0e1*n   !for the point at extended boundary
           b(n).x2(:,:,:)=1.0e2*n   !for the point at extended boundary
           b(n).y1(:,:,:)=1.0e3*n   !for the point at extended boundary
           b(n).y2(:,:,:)=1.0e4*n   !for the point at extended boundary
           b(n).z1(:,:,:)=1.0e5*n   !for the point at extended boundary
           b(n).z2(:,:,:)=1.0e6*n   !for the point at extended boundary
      end do

c----------------------------read mesh----------------------------
      do n=1,nb
c *** read mesh information of current block fn
         if (idimen .eq.2) then
            if(n.le.999) write(filename, '("bin_2d_",i3.3,".dat")') n
            if(n.gt.999) write(filename, '("bin_2d_",i4.4,".dat")') n
         else if(idimen.eq.3) then
            if(n.le.999) write(filename, '("bin_3d_",i3.3,".dat")') n
            if(n.gt.999) write(filename, '("bin_3d_",i4.4,".dat")') n
         endif
c        write(filename, '("bin_3d_",i3.3,".dat")') n
        open(9, file=filename, form='unformatted')

        if(idimen .eq. 2) then !2d 
           read(9) iltp, jltp
           kltp  = 2
           kl(n) = 2
        endif

        if(idimen .eq. 3) read(9) iltp, jltp, kltp

         if(iltp.ne.il(n).or.jltp.ne.jl(n).or.kltp.ne.kl(n)) then
          print *, 'dimension data in datain does not mach mesh'
          print *, 'il,jl,kl should be ',iltp,jltp,kltp
          stop
         end if
         
         print *, 'il,jl,kl is ',iltp,jltp,kltp

        if(idimen .eq. 3) then ! 3d read mesh 
        do k = 1, kl(n)+1
          do j = 1, jl(n)+1
            do i = 1, il(n)+1
              read(9) b(n).xx(i,j,k), 
     $                 b(n).yy(i,j,k), b(n).zz(i,j,k)
            end do
          end do
        end do

        else ! 2d read mesh
           k     = 1
          do j = 1, jl(n)+1
           do i = 1, il(n)+1
              read(9) b(n).xx(i,j,k), 
     $                 b(n).yy(i,j,k)
              b(n).zz(i,j,k)   = 0.0 !layer 1 z =0.0

              b(n).xx(i,j,k+1) = b(n).xx(i,j,k) !layer 2 , z = 0.5
              b(n).yy(i,j,k+1) = b(n).yy(i,j,k)
              b(n).zz(i,j,k+1) = 0.5

              b(n).xx(i,j,k+2) = b(n).xx(i,j,k) !layer 3  z = 1.0
              b(n).yy(i,j,k+2) = b(n).yy(i,j,k)
              b(n).zz(i,j,k+2) = 1.0
           end do
          end do
         
       endif
       close(9)

c----------center points on surface------------------------------------

        i=1
        do k=1,kl(n)
          do j=1,jl(n)
             b(n).x1(1,j,k)=(b(n).xx(i,j,k)+b(n).xx(i,j,k+1)
     $                    +b(n).xx(i,j+1,k)+b(n).xx(i,j+1,k+1))/4.d0
             b(n).x1(2,j,k)=(b(n).yy(i,j,k)+b(n).yy(i,j,k+1)
     $                    +b(n).yy(i,j+1,k)+b(n).yy(i,j+1,k+1))/4.d0
             b(n).x1(3,j,k)=(b(n).zz(i,j,k)+b(n).zz(i,j,k+1)
     $                    +b(n).zz(i,j+1,k)+b(n).zz(i,j+1,k+1))/4.d0    
          end do
        end do
        i=il(n)+1
        do k=1,kl(n)
          do j=1,jl(n)
             b(n).x2(1,j,k)=(b(n).xx(i,j,k)+b(n).xx(i,j,k+1)
     $                    +b(n).xx(i,j+1,k)+b(n).xx(i,j+1,k+1))/4.d0
             b(n).x2(2,j,k)=(b(n).yy(i,j,k)+b(n).yy(i,j,k+1)
     $                    +b(n).yy(i,j+1,k)+b(n).yy(i,j+1,k+1))/4.d0
             b(n).x2(3,j,k)=(b(n).zz(i,j,k)+b(n).zz(i,j,k+1)
     $                    +b(n).zz(i,j+1,k)+b(n).zz(i,j+1,k+1))/4.d0    
          end do
        end do
c--------------------------------------------                   
        j=1
        do k=1,kl(n)
          do i=1,il(n)
             b(n).y1(1,i,k)=(b(n).xx(i,j,k)+b(n).xx(i,j,k+1)
     $                    +b(n).xx(i+1,j,k)+b(n).xx(i+1,j,k+1))/4.d0
             b(n).y1(2,i,k)=(b(n).yy(i,j,k)+b(n).yy(i,j,k+1)
     $                    +b(n).yy(i+1,j,k)+b(n).yy(i+1,j,k+1))/4.d0
             b(n).y1(3,i,k)=(b(n).zz(i,j,k)+b(n).zz(i,j,k+1)
     $                    +b(n).zz(i+1,j,k)+b(n).zz(i+1,j,k+1))/4.d0   
          end do
        end do   
        j=jl(n)+1
        do k=1,kl(n)
          do i=1,il(n)
             b(n).y2(1,i,k)=(b(n).xx(i,j,k)+b(n).xx(i,j,k+1)
     $                    +b(n).xx(i+1,j,k)+b(n).xx(i+1,j,k+1))/4.d0
             b(n).y2(2,i,k)=(b(n).yy(i,j,k)+b(n).yy(i,j,k+1)
     $                    +b(n).yy(i+1,j,k)+b(n).yy(i+1,j,k+1))/4.d0
             b(n).y2(3,i,k)=(b(n).zz(i,j,k)+b(n).zz(i,j,k+1)
     $                    +b(n).zz(i+1,j,k)+b(n).zz(i+1,j,k+1))/4.d0  
          end do
        end do              
c--------------------------------------------                   
        k=1
        do j=1,jl(n)
          do i=1,il(n)
             b(n).z1(1,i,j)=(b(n).xx(i,j,k)+b(n).xx(i,j+1,k)
     $                    +b(n).xx(i+1,j,k)+b(n).xx(i+1,j+1,k))/4.d0
             b(n).z1(2,i,j)=(b(n).yy(i,j,k)+b(n).yy(i,j+1,k)
     $                    +b(n).yy(i+1,j,k)+b(n).yy(i+1,j+1,k))/4.d0
             b(n).z1(3,i,j)=(b(n).zz(i,j,k)+b(n).zz(i,j+1,k)
     $                    +b(n).zz(i+1,j,k)+b(n).zz(i+1,j+1,k))/4.d0     
          end do
        end do                   
        k=kl(n)+1
        do j=1,jl(n)
          do i=1,il(n)
             b(n).z2(1,i,j)=(b(n).xx(i,j,k)+b(n).xx(i,j+1,k)
     $                    +b(n).xx(i+1,j,k)+b(n).xx(i+1,j+1,k))/4.d0
             b(n).z2(2,i,j)=(b(n).yy(i,j,k)+b(n).yy(i,j+1,k)
     $                    +b(n).yy(i+1,j,k)+b(n).yy(i+1,j+1,k))/4.d0
             b(n).z2(3,i,j)=(b(n).zz(i,j,k)+b(n).zz(i,j+1,k)
     $                    +b(n).zz(i+1,j,k)+b(n).zz(i+1,j+1,k))/4.d0  
     
          end do
        end do       
      end do    
c---------------------------END OF read mesh----------------------     
c      print*, "finsh generate surface center"
c---------------------------make a box for each face--------------
      do n=1,nb
         do nf=1,6
            xmin(n,nf)=1e10
            xmax(n,nf)=-1e10
            ymin(n,nf)=1e10
            ymax(n,nf)=-1e10
            zmin(n,nf)=1e10
            zmax(n,nf)=-1e10         
            if (nf.eq.1) then
               i=1
               do k=1,kl(n)
                  do j=1,jl(n)
                     xmin(n,nf)=min(xmin(n,nf),b(n).x1(1,j,k))
                     xmax(n,nf)=max(xmax(n,nf),b(n).x1(1,j,k))
                     ymin(n,nf)=min(ymin(n,nf),b(n).x1(2,j,k))
                     ymax(n,nf)=max(ymax(n,nf),b(n).x1(2,j,k))
                     zmin(n,nf)=min(zmin(n,nf),b(n).x1(3,j,k))
                     zmax(n,nf)=max(zmax(n,nf),b(n).x1(3,j,k))          
                  end do
               end do
            elseif (nf.eq.2) then
               i=il(n)
               do k=1,kl(n)
                  do j=1,jl(n)
                     xmin(n,nf)=min(xmin(n,nf),b(n).x2(1,j,k))
                     xmax(n,nf)=max(xmax(n,nf),b(n).x2(1,j,k))
                     ymin(n,nf)=min(ymin(n,nf),b(n).x2(2,j,k))
                     ymax(n,nf)=max(ymax(n,nf),b(n).x2(2,j,k))
                     zmin(n,nf)=min(zmin(n,nf),b(n).x2(3,j,k))
                     zmax(n,nf)=max(zmax(n,nf),b(n).x2(3,j,k))          
                  end do
               end do
            elseif (nf.eq.3) then
               j=1
               do k=1,kl(n)
                  do i=1,il(n)                     
                     xmin(n,nf)=min(xmin(n,nf),b(n).y1(1,i,k))
                     xmax(n,nf)=max(xmax(n,nf),b(n).y1(1,i,k))
                     ymin(n,nf)=min(ymin(n,nf),b(n).y1(2,i,k))
                     ymax(n,nf)=max(ymax(n,nf),b(n).y1(2,i,k))
                     zmin(n,nf)=min(zmin(n,nf),b(n).y1(3,i,k))
                     zmax(n,nf)=max(zmax(n,nf),b(n).y1(3,i,k))     
                  end do
               end do
            elseif (nf.eq.4) then
               j=jl(n)
               do k=1,kl(n)
                  do i=1,il(n)                     
                     xmin(n,nf)=min(xmin(n,nf),b(n).y2(1,i,k))
                     xmax(n,nf)=max(xmax(n,nf),b(n).y2(1,i,k))
                     ymin(n,nf)=min(ymin(n,nf),b(n).y2(2,i,k))
                     ymax(n,nf)=max(ymax(n,nf),b(n).y2(2,i,k))
                     zmin(n,nf)=min(zmin(n,nf),b(n).y2(3,i,k))
                     zmax(n,nf)=max(zmax(n,nf),b(n).y2(3,i,k))         
                  end do
               end do
            elseif (nf.eq.5) then
               k=1
               do j=1,jl(n)
                  do i=1,il(n)                     
                     xmin(n,nf)=min(xmin(n,nf),b(n).z1(1,i,j))
                     xmax(n,nf)=max(xmax(n,nf),b(n).z1(1,i,j))
                     ymin(n,nf)=min(ymin(n,nf),b(n).z1(2,i,j))
                     ymax(n,nf)=max(ymax(n,nf),b(n).z1(2,i,j))
                     zmin(n,nf)=min(zmin(n,nf),b(n).z1(3,i,j))
                     zmax(n,nf)=max(zmax(n,nf),b(n).z1(3,i,j))          
                  end do
               end do
            else
               k=kl(n)
               do j=1,jl(n)
                  do i=1,il(n)                     
                     xmin(n,nf)=min(xmin(n,nf),b(n).z2(1,i,j))
                     xmax(n,nf)=max(xmax(n,nf),b(n).z2(1,i,j))
                     ymin(n,nf)=min(ymin(n,nf),b(n).z2(2,i,j))
                     ymax(n,nf)=max(ymax(n,nf),b(n).z2(2,i,j))
                     zmin(n,nf)=min(zmin(n,nf),b(n).z2(3,i,j))
                     zmax(n,nf)=max(zmax(n,nf),b(n).z2(3,i,j))          
                  end do
               end do
            end if
         end do
      end do

c      print*, "finsh making box for each surface"
c--------------------END OF-make a box for each face--------------   

c--------------------BEGIN OF-matching blocks---------------------   
      do n=1,NB
        b(n).flag1(:,:)='n'
        b(n).flag2(:,:)='n'
        b(n).flag3(:,:)='n'
        b(n).flag4(:,:)='n'
        b(n).flag5(:,:)='n'
        b(n).flag6(:,:)='n'        
      end do   
         
      do n=1,NB
        write(*,*)'Matching Block=',n
        do nf=1,2       !nf=1,6
        bcdir='xie'
        if (nf.eq.1) then
           i=1
           b(n).x(i,:,:)=b(n).x1(1,:,:)
           b(n).y(i,:,:)=b(n).x1(2,:,:)
           b(n).z(i,:,:)=b(n).x1(3,:,:)           
        elseif (nf.eq.2) then
           i=il(n)
           b(n).x(i,:,:)=b(n).x2(1,:,:)
           b(n).y(i,:,:)=b(n).x2(2,:,:)
           b(n).z(i,:,:)=b(n).x2(3,:,:)                      
        end if
        do j=1,jl(n)
          do k=1,kl(n)
             if ((nf.eq.1 .and.b(n).flag1(j,k).eq.'n').or.
     $            (nf.eq.2 .and.b(n).flag2(j,k).eq.'n') ) then
               do mb=n+1,NB
                  do nfb=1,6
                     if (b(n).x(i,j,k).ge.xmin(mb,nfb)-eps .and. !!!!BOX
     $                   b(n).x(i,j,k).le.xmax(mb,nfb)+eps .and.
     $                   b(n).y(i,j,k).ge.ymin(mb,nfb)-eps .and. 
     $                   b(n).y(i,j,k).le.ymax(mb,nfb)+eps .and.
     $                   b(n).z(i,j,k).ge.zmin(mb,nfb)-eps .and. 
     $                   b(n).z(i,j,k).le.zmax(mb,nfb)+eps) then     
                         if (nfb.eq.1) then
                           i1=1
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x1(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x1(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x1(3,:,:)
                              end if                                  
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag1(j1,k1)='y'    ! 1---nfb=1
                                  goto 10000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do 
                           
                         elseif (nfb.eq.2) then
                           i1=il(mb)
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x2(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x2(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag2(j1,k1)='y'    ! 2---nfb=2
                                  goto 10000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do           
                                          
                         elseif (nfb.eq.3) then
                           j1=1
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y1(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y1(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag3(i1,k1)='y'    ! 3---nfb=3
                                  goto 30000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                                       
                         elseif (nfb.eq.4) then
                           j1=jl(mb)
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y2(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y2(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y2(3,:,:)
                              end if                                  
                               
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then       
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag4(i1,k1)='y'    ! 4---nfb=4
                                  goto 30000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.5) then
                           k1=1
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z1(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z1(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag5(i1,j1)='y'    ! 5---nfb=5
                                  goto 50000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.6) then
                           k1=kl(mb)
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z2(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z2(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.1) then
                                      b(n).flag1(j,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.2) then
                                      b(n).flag2(j,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag6(i1,j1)='y'    ! 6---nfb=6
                                  goto 50000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         end if ! end if of    (nfb.eq.1)               
                     end if !end if !!!!BOX
                  end do ! end of (nfb=1,6)
               end do
c---------------------xi-->xi-------------------------------------------                              
10000          if ((nf.eq.1 .and.b(n).flag1(j,k).eq.'y') .or.
     $             (nf.eq.2 .and.b(n).flag2(j,k).eq.'y')) then !get the position of sub-block
                 order(1)=1
                 order(2)=2
                 order(3)=3
                 
                 if (nf.eq.1) then                 
                    start(1)=1
                    endp(1)=1
                 elseif (nf.eq.2) then 
                    start(1)=il(n)
                    endp(1)=il(n)
                 end if                  
                 start(2)=j
                 endp(2)=j
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.1) then
                     istart(1)=1
                     iendp(1)=1
                 elseif (nfb.eq.2) then
                     istart(1)=il(mb)
                     iendp(1)=il(mb)
                 end if                 
                 istart(2)=j1
                 iendp(2)=j1
                 istart(3)=k1
                 iendp(3)=k1    

                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                 
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 1101
                       order(2)=2
                       order(3)=3
                       endp(3)=endp(3)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 1101
                          order(2)=2
                          order(3)=3
                          endp(3)=endp(3)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          k4=iendp(3)
                          j4=iendp(2)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 1101
                             order(2)=3
                             order(3)=2
                             endp(3)=endp(3)+1
                             iendp(2)=iendp(2)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            j4=iendp(2)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 1101
                               order(2)=3
                               order(3)=2
                               endp(3)=endp(3)+1
                               iendp(2)=iendp(2)-1
                               case_11=4
                            else
                               goto 1101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
1101             case_22=0
                 do m2=1,jl(n)-j ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(3)-start(3)  
                       j3=j+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           j4=j1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 1111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 1111
                                   case_22=2
                               else
                                   goto 1111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           j4=j1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 1111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 1111
                                   case_22=2
                               else
                                   goto 1111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 1111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 1111
                                   case_22=2
                               else
                                   goto 1111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 1111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 1111
                                   case_22=2
                               else
                                   goto 1111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 1101       
1111             if (case_11.eq.1) then
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if   
                    
                    if (nf.eq.1) then                                  
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag1(j3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.2) then 
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag2(j3,k3)='y'
                          end do
                       end do
                    end if                                              
                 
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                             b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.2) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if 
                    if (nf.eq.1) then                   
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag1(j3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.2) then                   
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag2(j3,k3)='y'
                          end do
                       end do                       
                    end if
                                    
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.3) then   
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if   
                    if (nf.eq.1) then                 
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag1(j3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.2) then 
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag2(j3,k3)='y'
                          end do
                       end do
                    end if
                                  
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.4) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if   
                    if (nf.eq.1) then                     
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag1(j3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.2) then                     
                       do j3=start(2),endp(2)
                          do k3=start(3),endp(3)
                             b(n).flag2(j3,k3)='y'
                          end do
                       end do
                    end if

                    if (nfb.eq.1) then               
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic                
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)  
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 12345                 
               end if ! end if of 10000
c------------------------------xi-->eta---------------------------
30000          if ((nf.eq.1 .and.b(n).flag1(j,k).eq.'y') .or.
     $             (nf.eq.2 .and.b(n).flag2(j,k).eq.'y')) then !get the position of sub-block
                 order(1)=2           !xi--->eta
                 order(2)=3
                 order(3)=1
                                  
                 if (nf.eq.1) then                 
                    start(1)=1
                    endp(1)=1
                 elseif (nf.eq.2) then 
                    start(1)=il(n)
                    endp(1)=il(n)
                 end if                  

                 start(2)=j
                 endp(2)=j
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.3) then
                     istart(2)=1
                     iendp(2)=1
                 elseif (nfb.eq.4) then
                     istart(2)=jl(mb)
                     iendp(2)=jl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(3)=k1
                 iendp(3)=k1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 3101
                       order(2)=1
                       order(3)=3
                       endp(3)=endp(3)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 3101   
                          order(2)=1
                          order(3)=3
                          endp(3)=endp(3)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          k4=iendp(3)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 3101
                             order(2)=3
                             order(3)=1
                             endp(3)=endp(3)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 3101
                               order(2)=3
                               order(3)=1
                               endp(3)=endp(3)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 3101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
3101             case_22=0
                 do m2=1,jl(n)-j ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(3)-start(3)  
                       j3=j+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 3111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 3111
                                   case_22=2
                               else
                                   goto 3111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 3111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 3111
                                   case_22=2
                               else
                                   goto 3111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 3111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 3111
                                   case_22=2
                               else
                                   goto 3111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 3111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 3111
                                   case_22=2
                               else
                                   goto 3111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 3101
3111             if (case_11.eq.1) then
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3)
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3),-1
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if                    
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if                        
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)   
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 12345
                 
               end if ! end if of 30000
c               if(idimen .eq. 2) goto 55555
c------------------------------xi-->zeta---------------------------
50000          if ((nf.eq.1 .and.b(n).flag1(j,k).eq.'y') .or.
     $             (nf.eq.2 .and.b(n).flag2(j,k).eq.'y')) then !get the position of sub-block
                 order(1)=3           !xi--->zeta
                 order(2)=1
                 order(3)=2
                                  
                 if (nf.eq.1) then                 
                    start(1)=1
                    endp(1)=1
                 elseif (nf.eq.2) then 
                    start(1)=il(n)
                    endp(1)=il(n)
                 end if                  

                 start(2)=j
                 endp(2)=j
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.5) then
                     istart(3)=1
                     iendp(3)=1
                 elseif (nfb.eq.6) then
                     istart(3)=kl(mb)
                     iendp(3)=kl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(2)=j1
                 iendp(2)=j1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    j4=iendp(2)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 5101
                       order(2)=1
                       order(3)=2
                       endp(3)=endp(3)+1
                       iendp(2)=iendp(2)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       j4=iendp(2)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 5101   
                          order(2)=1
                          order(3)=2
                          endp(3)=endp(3)+1
                          iendp(2)=iendp(2)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          j4=iendp(2)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0 .and.case_11.ne.3)goto 5101
                             order(2)=2
                             order(3)=1
                             endp(3)=endp(3)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 5101
                               order(2)=2
                               order(3)=1
                               endp(3)=endp(3)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 5101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
5101             case_22=0
                 do m2=1,jl(n)-j ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(3)-start(3)  
                       j3=j+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 5111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 5111
                                   case_22=2
                               else
                                   goto 5111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 5111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 5111
                                   case_22=2
                               else
                                   goto 5111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           j4=j1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 5111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 5111
                                   case_22=2
                               else
                                   goto 5111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           j4=j1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 5111
                               case_22=1
                           else
                               j4=j1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 5111
                                   case_22=2
                               else
                                   goto 5111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 5101
5111             if (case_11.eq.1) then
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2)
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2),-1
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                       iendp(2)=j1+m2-1
                       ic=1
                    else
                       iendp(2)=j1-m2+1
                       ic=-1
                    end if                    
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(2)=j+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if                        
                    do j3=start(2),endp(2)
                       do k3=start(3),endp(3)
                          if (nf.eq.1) then
                             b(n).flag1(j3,k3)='y'
                          elseif (nf.eq.2) then
                             b(n).flag2(j3,k3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)  
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 12345

               end if ! end if of 50000

               
             end if   ! end if of (b(n).flag1(j,k).eq.'n' )      
12345   i3=1        
          end do
        end do
        end do ! END of do (nf=1,2)


c--------------343434343434343434----------------------------------------  
55555   do nf=3,4               !nf=1,6
        bcdir='eta'
        if (nf.eq.3) then
           j=1
           b(n).x(:,j,:)=b(n).y1(1,:,:)
           b(n).y(:,j,:)=b(n).y1(2,:,:)
           b(n).z(:,j,:)=b(n).y1(3,:,:)           
        elseif (nf.eq.4) then
           j=jl(n)
           b(n).x(:,j,:)=b(n).y2(1,:,:)
           b(n).y(:,j,:)=b(n).y2(2,:,:)
           b(n).z(:,j,:)=b(n).y2(3,:,:)                      
        end if
        do i=1,il(n)
          do k=1,kl(n)
             if ((nf.eq.3 .and.b(n).flag3(i,k).eq.'n').or.
     $            (nf.eq.4 .and.b(n).flag4(i,k).eq.'n') ) then
               do mb=n+1,NB
                  do nfb=1,6
                     if (b(n).x(i,j,k).ge.xmin(mb,nfb)-eps .and. !!!!BOX
     $                   b(n).x(i,j,k).le.xmax(mb,nfb)+eps .and.
     $                   b(n).y(i,j,k).ge.ymin(mb,nfb)-eps .and. 
     $                   b(n).y(i,j,k).le.ymax(mb,nfb)+eps .and.
     $                   b(n).z(i,j,k).ge.zmin(mb,nfb)-eps .and. 
     $                   b(n).z(i,j,k).le.zmax(mb,nfb)+eps) then
                         if (nfb.eq.1) then
                           i1=1
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x1(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x1(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x1(3,:,:)
                              end if                                  
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag1(j1,k1)='y'    ! 1---nfb=1
                                  goto 20000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do 
                           
                         elseif (nfb.eq.2) then
                           i1=il(mb)
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x2(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x2(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag2(j1,k1)='y'    ! 2---nfb=2
                                  goto 20000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do           
                                          
                         elseif (nfb.eq.3) then
                           j1=1
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y1(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y1(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag3(i1,k1)='y'    ! 3---nfb=3
                                  goto 40000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                                       
                         elseif (nfb.eq.4) then
                           j1=jl(mb)
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y2(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y2(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y2(3,:,:)
                              end if                                  
                               
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag4(i1,k1)='y'    ! 4---nfb=4
                                  goto 40000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.5) then
                           k1=1
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z1(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z1(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag5(i1,j1)='y'    ! 5---nfb=5
                                  goto 60000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.6) then
                           k1=kl(mb)
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z2(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z2(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.3) then
                                      b(n).flag3(i,k)='y'   ! 1---nf=1
                                  elseif (nf.eq.4) then
                                      b(n).flag4(i,k)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag6(i1,j1)='y'    ! 6---nfb=6
                                  goto 60000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         end if ! end if of    (nfb.eq.1)               
                     end if !end if !!!!BOX
                  end do ! end of (nfb=1,6)
               end do
c---------------------xi-->xi-------------------------------------------                              
20000          if ((nf.eq.3 .and.b(n).flag3(i,k).eq.'y') .or.
     $             (nf.eq.4 .and.b(n).flag4(i,k).eq.'y')) then !get the position of sub-block
                 order(2)=1
                 order(1)=3
                 order(3)=2
                 
                 if (nf.eq.3) then                 
                    start(2)=1
                    endp(2)=1
                 elseif (nf.eq.4) then 
                    start(2)=jl(n)
                    endp(2)=jl(n)
                 end if                  
                 start(1)=i
                 endp(1)=i
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.1) then
                     istart(1)=1
                     iendp(1)=1
                 elseif (nfb.eq.2) then
                     istart(1)=il(mb)
                     iendp(1)=il(mb)
                 end if                 
                 istart(2)=j1
                 iendp(2)=j1
                 istart(3)=k1
                 iendp(3)=k1    

                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                 
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 2101
                       order(1)=2
                       order(3)=3
                       endp(3)=endp(3)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 2101
                          order(1)=2
                          order(3)=3
                          endp(3)=endp(3)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          k4=iendp(3)
                          j4=iendp(2)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 2101
                             order(1)=3
                             order(3)=2
                             endp(3)=endp(3)+1
                             iendp(2)=iendp(2)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            j4=iendp(2)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 2101
                               order(1)=3
                               order(3)=2
                               endp(3)=endp(3)+1
                               iendp(2)=iendp(2)-1
                               case_11=4
                            else
                               goto 2101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
2101             case_22=0
                 do m2=1,il(n)-i
                    do m22=1,endp(3)-start(3)  
                       i3=i+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           j4=j1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 2111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 2111
                                   case_22=2
                               else
                                   goto 2111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           j4=j1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 2111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 2111
                                   case_22=2
                               else
                                   goto 2111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 2111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 2111
                                   case_22=2
                               else
                                   goto 2111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 2111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 2111
                                   case_22=2
                               else
                                   goto 2111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 2101       
2111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if   
                    
                    if (nf.eq.3) then                                  
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag3(i3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.4) then 
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag4(i3,k3)='y'
                          end do
                       end do
                    end if                                              
                 
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                             b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if 
                    if (nf.eq.3) then                   
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag3(i3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.4) then                   
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag4(i3,k3)='y'
                          end do
                       end do                       
                    end if
                                    
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if   
                    if (nf.eq.3) then                 
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag3(i3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.4) then 
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag4(i3,k3)='y'
                          end do
                       end do
                    end if
                                  
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if   
                    if (nf.eq.3) then                     
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag3(i3,k3)='y'
                          end do
                       end do
                    elseif (nf.eq.4) then                     
                       do i3=start(1),endp(1)
                          do k3=start(3),endp(3)
                             b(n).flag4(i3,k3)='y'
                          end do
                       end do
                    end if

                    if (nfb.eq.1) then               
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic                
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 22345                 
               end if ! end if of 10000
c------------------------------xi-->eta---------------------------
40000          if ((nf.eq.3 .and.b(n).flag3(i,k).eq.'y') .or.
     $             (nf.eq.4 .and.b(n).flag4(i,k).eq.'y')) then !get the position of sub-block

                 order(2)=2           !eta--->eta
                 order(1)=3           
                 order(3)=1
                                                                                                     
                 if (nf.eq.3) then                 
                    start(2)=1
                    endp(2)=1
                 elseif (nf.eq.4) then 
                    start(2)=jl(n)
                    endp(2)=jl(n)
                 end if                  

                 start(1)=i
                 endp(1)=i
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.3) then
                     istart(2)=1
                     iendp(2)=1
                 elseif (nfb.eq.4) then
                     istart(2)=jl(mb)
                     iendp(2)=jl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(3)=k1
                 iendp(3)=k1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 4101
                       order(1)=1
                       order(3)=3
                       endp(3)=endp(3)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 4101   
                          order(1)=1
                          order(3)=3
                          endp(3)=endp(3)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          k4=iendp(3)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 4101
                             order(1)=3
                             order(3)=1
                             endp(3)=endp(3)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 4101
                               order(1)=3
                               order(3)=1
                               endp(3)=endp(3)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 4101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
4101             case_22=0
                 do m2=1,il(n)-i ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(3)-start(3)  
                       i3=i+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 4111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 4111
                                   case_22=2
                               else
                                   goto 4111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 4111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 4111
                                   case_22=2
                               else
                                   goto 4111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 4111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 4111
                                   case_22=2
                               else
                                   goto 4111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 4111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 4111
                                   case_22=2
                               else
                                   goto 4111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 4101
4111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3)
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3),-1
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if                        
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 22345
                 
               end if ! end if of 40000
c               if(idimen.eq.2) goto 66666
c------------------------------xi-->zeta---------------------------
60000          if ((nf.eq.3 .and.b(n).flag3(i,k).eq.'y') .or.
     $             (nf.eq.4 .and.b(n).flag4(i,k).eq.'y')) then !get the position of sub-block
                 order(2)=3           !zeta--->zeta
                 order(1)=2
                 order(3)=1
                                  
                 if (nf.eq.3) then                 
                    start(2)=1
                    endp(2)=1
                 elseif (nf.eq.4) then 
                    start(2)=jl(n)
                    endp(2)=jl(n)
                 end if                  

                 start(1)=i
                 endp(1)=i
                 start(3)=k
                 endp(3)=k
                 if (nfb.eq.5) then
                     istart(3)=1
                     iendp(3)=1
                 elseif (nfb.eq.6) then
                     istart(3)=kl(mb)
                     iendp(3)=kl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(2)=j1
                 iendp(2)=j1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,kl(n)-k
                    k3=endp(3)+1
                    j4=iendp(2)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 6101
                       order(1)=1
                       order(3)=2
                       endp(3)=endp(3)+1
                       iendp(2)=iendp(2)+1
                       case_11=1
                    else
                       k3=endp(3)+1
                       j4=iendp(2)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 6101   
                          order(1)=1
                          order(3)=2
                          endp(3)=endp(3)+1
                          iendp(2)=iendp(2)-1
                          case_11=2                          
                        else
                          k3=endp(3)+1
                          j4=iendp(2)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0 .and.case_11.ne.3)goto 6101
                             order(1)=2
                             order(3)=1
                             endp(3)=endp(3)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            k3=endp(3)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 6101
                               order(1)=2
                               order(3)=1
                               endp(3)=endp(3)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 6101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
6101             case_22=0
                 do m2=1,il(n)-i ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(3)-start(3)  
                       i3=i+m2
                       k3=k+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 6111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 6111
                                   case_22=2
                               else
                                   goto 6111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 6111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 6111
                                   case_22=2
                               else
                                   goto 6111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           j4=j1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 6111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 6111
                                   case_22=2
                               else
                                   goto 6111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           j4=j1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 6111
                               case_22=1
                           else
                               j4=j1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 6111
                                   case_22=2
                               else
                                   goto 6111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 6101
6111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2)
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2),-1
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(2)=j1+m2-1
                       ic=1
                    else
                       iendp(2)=j1-m2+1
                       ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if                        
                    do i3=start(1),endp(1)
                       do k3=start(3),endp(3)
                          if (nf.eq.3) then
                             b(n).flag3(i3,k3)='y'
                          elseif (nf.eq.4) then
                             b(n).flag4(i3,k3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 end if

                 if(idimen.eq.2) then
                    start(3) = 1
                    endp(3) = 1
                    istart(3) = 1
                    iendp(3) = 1
                 endif

                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif 
               goto 22345

               end if ! end if of 60000

               
             end if   ! end if of (b(n).flag3(i,k).eq.'n' )    
22345   i3=1                                    
          end do
        end do                
        end do         !nf=3,4


c--------------5656565656565656565656556-------------------------------------
66666   if(idimen .eq. 2) goto 54321
        do nf=5,6       !nf=5,6        
        bcdir='zta'
        if (nf.eq.5) then
           k=1
           b(n).x(:,:,k)=b(n).z1(1,:,:)
           b(n).y(:,:,k)=b(n).z1(2,:,:)
           b(n).z(:,:,k)=b(n).z1(3,:,:)           
        elseif (nf.eq.6) then
           k=kl(n)
           b(n).x(:,:,k)=b(n).z2(1,:,:)
           b(n).y(:,:,k)=b(n).z2(2,:,:)
           b(n).z(:,:,k)=b(n).z2(3,:,:)              
        end if
        do i=1,il(n)
          do j=1,jl(n)
             if ((nf.eq.5 .and.b(n).flag5(i,j).eq.'n').or.
     $            (nf.eq.6 .and.b(n).flag6(i,j).eq.'n') ) then
               do mb=n+1,NB
                  do nfb=1,6
                     if (b(n).x(i,j,k).ge.xmin(mb,nfb)-eps .and. !!!!BOX
     $                   b(n).x(i,j,k).le.xmax(mb,nfb)+eps .and.
     $                   b(n).y(i,j,k).ge.ymin(mb,nfb)-eps .and. 
     $                   b(n).y(i,j,k).le.ymax(mb,nfb)+eps .and.
     $                   b(n).z(i,j,k).ge.zmin(mb,nfb)-eps .and. 
     $                   b(n).z(i,j,k).le.zmax(mb,nfb)+eps) then
                         if (nfb.eq.1) then
                           i1=1
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x1(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x1(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x1(3,:,:)
                              end if                                  
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag1(j1,k1)='y'    ! 1---nfb=1
                                  goto 70000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do 
                           
                         elseif (nfb.eq.2) then
                           i1=il(mb)
                           do j1=1,jl(mb)
                             do k1=1,kl(mb)
                              if (j1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(i1,:,:)=b(mb).x2(1,:,:)
                                  b(mb).y(i1,:,:)=b(mb).x2(2,:,:)
                                  b(mb).z(i1,:,:)=b(mb).x2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag2(j1,k1)='y'    ! 2---nfb=2
                                  goto 70000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do           
                                          
                         elseif (nfb.eq.3) then
                           j1=1
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y1(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y1(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag3(i1,k1)='y'    ! 3---nfb=3
                                  goto 80000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                                       
                         elseif (nfb.eq.4) then
                           j1=jl(mb)
                           do i1=1,il(mb)
                             do k1=1,kl(mb)
                              if (i1.eq.1 .and. k1.eq.1) then
                                  b(mb).x(:,j1,:)=b(mb).y2(1,:,:)
                                  b(mb).y(:,j1,:)=b(mb).y2(2,:,:)
                                  b(mb).z(:,j1,:)=b(mb).y2(3,:,:)
                              end if                                  
                               
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag4(i1,k1)='y'    ! 4---nfb=4
                                  goto 80000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.5) then
                           k1=1
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z1(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z1(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z1(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag5(i1,j1)='y'    ! 5---nfb=5
                                  goto 90000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         elseif (nfb.eq.6) then
                           k1=kl(mb)
                           do i1=1,il(mb)
                             do j1=1,jl(mb)
                              if (i1.eq.1 .and. j1.eq.1) then
                                  b(mb).x(:,:,k1)=b(mb).z2(1,:,:)
                                  b(mb).y(:,:,k1)=b(mb).z2(2,:,:)
                                  b(mb).z(:,:,k1)=b(mb).z2(3,:,:)
                              end if                                  
                             
                              temp1=abs(b(n).x(i,j,k)-b(mb).x(i1,j1,k1))
                              temp2=abs(b(n).y(i,j,k)-b(mb).y(i1,j1,k1))
                              temp3=abs(b(n).z(i,j,k)-b(mb).z(i1,j1,k1))
                              if (temp1+temp2+temp3.le.eps) then
                                  if (nf.eq.5) then
                                      b(n).flag5(i,j)='y'   ! 1---nf=1
                                  elseif (nf.eq.6) then
                                      b(n).flag6(i,j)='y'   ! 2---nf=2
                                  end if                              
                                  b(mb).flag6(i1,j1)='y'    ! 6---nfb=6
                                  goto 90000                     
                              end if ! end if (temp1+...+temp3)
                             end do
                           end do                          
                         
                         end if ! end if of    (nfb.eq.1)               
                     end if !end if !!!!BOX
                  end do ! end of (nfb=1,6)
               end do
c---------------------xi-->xi-------------------------------------------                              
70000          if ((nf.eq.5 .and.b(n).flag5(i,j).eq.'y') .or.
     $             (nf.eq.6 .and.b(n).flag6(i,j).eq.'y')) then !get the position of sub-block
                 order(3)=1  ! Zeta---->xi
                 order(2)=3
                 order(1)=2
                 
                 if (nf.eq.5) then                 
                    start(3)=1
                    endp(3)=1
                 elseif (nf.eq.6) then 
                    start(3)=kl(n)
                    endp(3)=kl(n)
                 end if                  
                 start(1)=i
                 endp(1)=i
                 start(2)=j
                 endp(2)=j
                 if (nfb.eq.1) then
                     istart(1)=1
                     iendp(1)=1
                 elseif (nfb.eq.2) then
                     istart(1)=il(mb)
                     iendp(1)=il(mb)
                 end if                 
                 istart(2)=j1
                 iendp(2)=j1
                 istart(3)=k1
                 iendp(3)=k1    

                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                 
                 case_11=0
                 do m2=1,jl(n)-j
                    j3=endp(2)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 7101
                       order(1)=2
                       order(2)=3
                       endp(2)=endp(2)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       j3=endp(2)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 7101
                          order(1)=2
                          order(2)=3
                          endp(2)=endp(2)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          j3=endp(2)+1
                          k4=iendp(3)
                          j4=iendp(2)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 7101
                             order(1)=3
                             order(2)=2
                             endp(2)=endp(2)+1
                             iendp(2)=iendp(2)+1
                             case_11=3                                  
                          else
                            j3=endp(2)+1
                            j4=iendp(2)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 7101
                               order(1)=3
                               order(2)=2
                               endp(2)=endp(2)+1
                               iendp(2)=iendp(2)-1
                               case_11=4
                            else
                               goto 7101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
7101             case_22=0
                 do m2=1,il(n)-i
                    do m22=1,endp(2)-start(2)  
                       i3=i+m2
                       j3=j+m22                                      
                       if (case_11.eq.1) then
                           j4=j1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 7111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 7111
                                   case_22=2
                               else
                                   goto 7111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           j4=j1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 7111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 7111
                                   case_22=2
                               else
                                   goto 7111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 7111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 7111
                                   case_22=2
                               else
                                   goto 7111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 7111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 7111
                                   case_22=2
                               else
                                   goto 7111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 7101       
7111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if   
                    
                    if (nf.eq.5) then                                  
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag5(i3,j3)='y'
                          end do
                       end do
                    elseif (nf.eq.6) then 
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag6(i3,j3)='y'
                          end do
                       end do
                    end if                                              
                 
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                             b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3)
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if 
                    if (nf.eq.5) then                   
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag5(i3,j3)='y'
                          end do
                       end do
                    elseif (nf.eq.6) then                   
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag6(i3,j3)='y'
                          end do
                       end do                       
                    end if
                                    
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag1(j4,k4)='y'
                          end do
                       end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),ic
                          do k4=istart(3),iendp(3),-1
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if   
                    if (nf.eq.5) then                 
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag5(i3,j3)='y'
                          end do
                       end do
                    elseif (nf.eq.6) then 
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag6(i3,j3)='y'
                          end do
                       end do
                    end if
                                  
                    if (nfb.eq.1) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                     elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2)
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if   
                    if (nf.eq.5) then                     
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag5(i3,j3)='y'
                          end do
                       end do
                    elseif (nf.eq.6) then                     
                       do i3=start(1),endp(1)
                          do j3=start(2),endp(2)
                             b(n).flag6(i3,j3)='y'
                          end do
                       end do
                    end if

                    if (nfb.eq.1) then               
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic
                              b(mb).flag1(j4,k4)='y'
                           end do
                        end do
                    elseif (nfb.eq.2) then
                       do j4=istart(2),iendp(2),-1
                          do k4=istart(3),iendp(3),ic                
                              b(mb).flag2(j4,k4)='y'
                          end do                 
                       end do
                    end if
                 end if
                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3) 
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif 
               goto 32345                 
               end if ! end if of 10000
c------------------------------xi-->eta---------------------------
80000          if ((nf.eq.5 .and.b(n).flag5(i,j).eq.'y') .or.
     $             (nf.eq.6 .and.b(n).flag6(i,j).eq.'y')) then !get the position of sub-block

                 order(3)=2           !zeta--->eta
                 order(1)=3           
                 order(2)=1
                                                                                                     
                 if (nf.eq.5) then                 
                    start(3)=1
                    endp(3)=1
                 elseif (nf.eq.6) then 
                    start(3)=kl(n)
                    endp(3)=kl(n)
                 end if                  

                 start(1)=i
                 endp(1)=i
                 start(2)=j
                 endp(2)=j
                 if (nfb.eq.3) then
                     istart(2)=1
                     iendp(2)=1
                 elseif (nfb.eq.4) then
                     istart(2)=jl(mb)
                     iendp(2)=jl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(3)=k1
                 iendp(3)=k1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,jl(n)-j
                    j3=endp(2)+1
                    k4=iendp(3)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 8101
                       order(1)=1
                       order(2)=3
                       endp(2)=endp(2)+1
                       iendp(3)=iendp(3)+1
                       case_11=1
                    else
                       j3=endp(2)+1
                       k4=iendp(3)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 8101   
                          order(1)=1
                          order(2)=3
                          endp(2)=endp(2)+1
                          iendp(3)=iendp(3)-1
                          case_11=2                          
                        else
                          j3=endp(2)+1
                          k4=iendp(3)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0.and.case_11.ne.3) goto 8101
                             order(1)=3
                             order(2)=1
                             endp(2)=endp(2)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            j3=endp(2)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 8101
                               order(1)=3
                               order(2)=1
                               endp(2)=endp(2)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 8101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
8101             case_22=0
                 do m2=1,il(n)-i ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(2)-start(2)  
                       i3=i+m2
                       j3=j+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           k4=k1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 8111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 8111
                                   case_22=2
                               else
                                   goto 8111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           k4=k1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 8111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 8111
                                   case_22=2
                               else
                                   goto 8111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           k4=k1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 8111
                               case_22=1
                           else
                               k4=k1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 8111
                                   case_22=2
                               else
                                   goto 8111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           k4=k1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 8111
                               case_22=1
                           else
                               k4=k1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 8111
                                   case_22=2
                               else
                                   goto 8111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 8101
8111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3)
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do k4=istart(3),iendp(3),-1
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(3)=k1+m2-1
                       ic=1
                    else
                       iendp(3)=k1-m2+1
                       ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(3)=k1+m2-1
                        ic=1
                    else
                        iendp(3)=k1-m2+1
                        ic=-1
                    end if                        
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do k4=istart(3),iendp(3),ic
                          if (nfb.eq.3) then
                              b(mb).flag3(i4,k4)='y'
                          elseif (nfb.eq.4) then
                              b(mb).flag4(i4,k4)='y'
                          end if                 
                       end do
                    end do
                 end if
                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3)
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif
               goto 32345
                 
               end if ! end if of 80000

c------------------------------xi-->zeta---------------------------
90000          if ((nf.eq.5 .and.b(n).flag5(i,j).eq.'y') .or.
     $             (nf.eq.6 .and.b(n).flag6(i,j).eq.'y')) then !get the position of sub-block
                 order(3)=3           !zeta--->zeta
                 order(1)=1
                 order(2)=2
                                  
                                  
                 if (nf.eq.5) then                 
                    start(3)=1
                    endp(3)=1
                 elseif (nf.eq.6) then 
                    start(3)=kl(n)
                    endp(3)=kl(n)
                 end if                  

                 start(1)=i
                 endp(1)=i
                 start(2)=j
                 endp(2)=j
                 if (nfb.eq.5) then
                     istart(3)=1
                     iendp(3)=1
                 elseif (nfb.eq.6) then
                     istart(3)=kl(mb)
                     iendp(3)=kl(mb)
                 end if                 
                 istart(1)=i1
                 iendp(1)=i1
                 istart(2)=j1
                 iendp(2)=j1
                 
                 i3=start(1)         
                 j3=start(2)
                 k3=start(3)
                 i4=istart(1)         
                 j4=istart(2)
                 k4=istart(3)
                     
                 case_11=0
                 do m2=1,jl(n)-j
                    j3=endp(2)+1
                    j4=iendp(2)+1
                    temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                    temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                    temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                    if (temp1+temp2+temp3.le.eps) then
                       if(case_11.ne.0 .and.case_11.ne.1) goto 9101
                       order(1)=1
                       order(2)=2
                       endp(2)=endp(2)+1
                       iendp(2)=iendp(2)+1
                       case_11=1
                    else
                       j3=endp(2)+1
                       j4=iendp(2)-1
                       temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                       temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                       temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                       if (temp1+temp2+temp3.le.eps) then
                          if(case_11.ne.0 .and.case_11.ne.2) goto 9101   
                          order(1)=1
                          order(2)=2
                          endp(2)=endp(2)+1
                          iendp(2)=iendp(2)-1
                          case_11=2                          
                        else
                          j3=endp(2)+1
                          j4=iendp(2)
                          i4=iendp(1)+1
                          temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                          temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                          temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                          if (temp1+temp2+temp3.le.eps) then
                             if(case_11.ne.0 .and.case_11.ne.3)goto 9101
                             order(1)=2
                             order(2)=1
                             endp(2)=endp(2)+1
                             iendp(1)=iendp(1)+1
                             case_11=3                                  
                          else
                            j3=endp(2)+1
                            i4=iendp(1)-1
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                            if (temp1+temp2+temp3.le.eps) then
                              if(case_11.ne.0.and.case_11.ne.4)goto 9101
                               order(1)=2
                               order(2)=1
                               endp(2)=endp(2)+1
                               iendp(1)=iendp(1)-1
                               case_11=4
                            else
                               goto 9101
                            end if
                          end if
                       end if    
                    end if        
                 end do  ! end do of (m2=1,kl(n)-k)
9101             case_22=0
                 do m2=1,il(n)-i ! end do of (m2=1,kl(n)-k)
                    do m22=1,endp(2)-start(2)  
                       i3=i+m2
                       j3=j+m22                                      
                       if (case_11.eq.1) then
                           i4=i1+m2
                           j4=j1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 9111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 9111
                                   case_22=2
                               else
                                   goto 9111
                               end if
                           end if
                       elseif (case_11.eq.2) then 
                           i4=i1+m2
                           j4=j1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 9111
                               case_22=1
                           else
                               i4=i1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 9111
                                   case_22=2
                               else
                                   goto 9111
                               end if
                           end if
                       elseif (case_11.eq.3) then 
                           j4=j1+m2
                           i4=i1+m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 9111
                               case_22=1
                           else
                               j4=j1-m2
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 9111
                                   case_22=2
                               else
                                   goto 9111
                               end if
                           end if                           
                       elseif (case_11.eq.4) then 
                           j4=j1+m2
                           i4=i1-m22
                           temp1=abs(b(n).x(i3,j3,k3)-b(mb).x(i4,j4,k4))
                           temp2=abs(b(n).y(i3,j3,k3)-b(mb).y(i4,j4,k4))
                           temp3=abs(b(n).z(i3,j3,k3)-b(mb).z(i4,j4,k4))
                           if (temp1+temp2+temp3.le.eps) then
                               if (case_22.ne.0.and.case_22.ne.1)
     $                               goto 9111
                               case_22=1
                           else
                               j4=j1-m2                           
                                   temp1=abs(b(n).x(i3,j3,k3)-
     $                                       b(mb).x(i4,j4,k4))
                                   temp2=abs(b(n).y(i3,j3,k3)-
     $                                       b(mb).y(i4,j4,k4))
                                   temp3=abs(b(n).z(i3,j3,k3)-
     $                                       b(mb).z(i4,j4,k4))
                               if (temp1+temp2+temp3.le.eps) then
                                   if (case_22.ne.0.and.case_22.ne.2)
     $                                 goto 9111
                                   case_22=2
                               else
                                   goto 9111
                               end if
                           end if                           
                      end if
                    end do
                 end do     ! end of 9101
9111             if (case_11.eq.1) then
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                                        
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if                             
                       end do
                    end do
                 
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2)
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.2) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(1)=i1+m2-1
                        ic=1
                    else
                        iendp(1)=i1-m2+1
                        ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if                          
                       end do
                    end do
                                    
                    do i4=istart(1),iendp(1),ic
                       do j4=istart(2),iendp(2),-1
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.3) then   
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                       iendp(2)=j1+m2-1
                       ic=1
                    else
                       iendp(2)=j1-m2+1
                       ic=-1
                    end if                    
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if
                       end do
                    end do
                                  
                    do i4=istart(1),iendp(1)
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 elseif (case_11.eq.4) then  
                    endp(1)=i+m2-1
                    if (case_22.eq.1) then
                        iendp(2)=j1+m2-1
                        ic=1
                    else
                        iendp(2)=j1-m2+1
                        ic=-1
                    end if                        
                    do i3=start(1),endp(1)
                       do j3=start(2),endp(2)
                          if (nf.eq.5) then
                             b(n).flag5(i3,j3)='y'
                          elseif (nf.eq.6) then
                             b(n).flag6(i3,j3)='y'
                          end if
                       end do
                    end do
                                   
                    do i4=istart(1),iendp(1),-1
                       do j4=istart(2),iendp(2),ic
                          if (nfb.eq.5) then
                              b(mb).flag5(i4,j4)='y'
                          elseif (nfb.eq.6) then
                              b(mb).flag6(i4,j4)='y'
                          end if                 
                       end do
                    end do
                 end if
                 write(12,998) bcdir,n,start(1),start(2),start(3),
     $                             endp(1),endp(2),endp(3)
                 write(12,999)   mb,istart(1),istart(2),istart(3),
     $                             iendp(1),iendp(2),iendp(3),
     $                             order(1),order(2),order(3) 
calex
c detection of interface error and output warning message
                 if (((endp(1)-start(1)+endp(2)-start(2)).le.width)
     >           .or.((endp(2)-start(2)+endp(3)-start(3)).le.width)
     >           .or.((endp(3)-start(3)+endp(1)-start(1)).le.width))then
                 write(2,*) 'WARNING !!! Possible error between block'
     > ,n,'and',mb,'please check you mesh'
                 endif  
               goto 32345

               end if ! end if of 90000

               
             end if   ! end if of (b(n).flag3(i,k).eq.'n' )   
32345     i3=1
          end do
        end do                
        end do         !nf=5,6
        
        write(12,*)
        write(12,*)                                     
        
c---------------------------------------------------------------------------        
54321 end  do                   ! end of (n=1,NB)




      
c---------------block boundary needs to modify-----------------------------
      do n=1,NB
         nf=1
         bcdir='xie'
         start(1)=1
         endp(1)=1
10       case_11=0
         do j=1,jl(n)
            do k=1,kl(n)
               if (case_11.eq.0 .and. b(n).flag1(j,k).eq.'n') then
                   start(2)=j
                   endp(2)=j
                   start(3)=k
                   endp(3)=k
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag1(j,k).eq.'n') then
                   endp(3)=endp(3)+1
                elseif (case_11.eq.1 .and. b(n).flag1(j,k).eq.'y') then
                    goto 100                
                end if
             end do
             if (case_11.eq.1) goto 100       
          end do
          if (case_11.eq.0) goto 201
100       do m2=1,jl(n)-j
             j3=j+m2
             do k3=start(3),endp(3)
                if (b(n).flag1(j3,k3).eq.'y') goto 1000
             end do
          end do
1000      endp(2)=endp(2)+m2-1
          do j=start(2),endp(2)
             do k=start(3),endp(3)
                b(n).flag1(j,k)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

          write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 10
c------------------------
201      nf=2      

         bcdir='xie'
         start(1)=il(n)
         endp(1)=il(n)
20       case_11=0
         do j=1,jl(n)
            do k=1,kl(n)
               if (case_11.eq.0 .and. b(n).flag2(j,k).eq.'n') then
                   start(2)=j
                   endp(2)=j
                   start(3)=k
                   endp(3)=k
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag2(j,k).eq.'n') then
                   endp(3)=endp(3)+1
                elseif (case_11.eq.1 .and. b(n).flag2(j,k).eq.'y') then
                    goto 200
                end if
             end do
             if (case_11.eq.1) goto 200       
          end do
          if (case_11.eq.0) goto 301
200       do m2=1,jl(n)-j
             j3=j+m2
             do k3=start(3),endp(3)
                if (b(n).flag2(j3,k3).eq.'y') goto 2000
             end do
          end do
2000      endp(2)=endp(2)+m2-1
          do j=start(2),endp(2)
             do k=start(3),endp(3)
                b(n).flag2(j,k)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

           write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 20
c------------------------

301      nf=3

         bcdir='eta'
         start(2)=1
         endp(2)=1
30       case_11=0
         do i=1,il(n)
            do k=1,kl(n)
               if (case_11.eq.0 .and. b(n).flag3(i,k).eq.'n') then
                   start(1)=i
                   endp(1)=i
                   start(3)=k
                   endp(3)=k
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag3(i,k).eq.'n') then
                   endp(3)=endp(3)+1
                elseif (case_11.eq.1 .and. b(n).flag3(i,k).eq.'y') then
                   goto 300
                end if
             end do
             if (case_11.eq.1) goto 300       
          end do
          if (case_11.eq.0) goto 401
300       do m2=1,il(n)-i
             i3=i+m2
             do k3=start(3),endp(3)
                if (b(n).flag3(i3,k3).eq.'y') goto 3000
             end do
          end do
3000      endp(1)=endp(1)+m2-1
          do i=start(1),endp(1)
             do k=start(3),endp(3)
                b(n).flag3(i,k)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

           write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 30
c------------------------
401      nf=4

         bcdir='eta'
         start(2)=jl(n)
         endp(2)=jl(n)
40       case_11=0
         do i=1,il(n)
            do k=1,kl(n)
               if (case_11.eq.0 .and. b(n).flag4(i,k).eq.'n') then
                   start(1)=i
                   endp(1)=i
                   start(3)=k
                   endp(3)=k
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag4(i,k).eq.'n') then
                   endp(3)=endp(3)+1
                elseif (case_11.eq.1 .and. b(n).flag4(i,k).eq.'y') then
                   goto 400
                end if
             end do
             if (case_11.eq.1) goto 400       
          end do
          if (case_11.eq.0) goto 501
400       do m2=1,il(n)-i
             i3=i+m2
             do k3=start(3),endp(3)
                if (b(n).flag4(i3,k3).eq.'y') goto 4000
             end do
          end do
4000      endp(1)=endp(1)+m2-1
          do i=start(1),endp(1)
             do k=start(3),endp(3)
                b(n).flag4(i,k)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

           write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 40
c--------------------------------
501      nf=5

          if(idimen.eq.2) goto 701

         bcdir='zta'
         start(3)=1
         endp(3)=1
50       case_11=0
         do i=1,il(n)
            do j=1,jl(n)
               if (case_11.eq.0 .and. b(n).flag5(i,j).eq.'n') then
                   start(1)=i
                   endp(1)=i
                   start(2)=j
                   endp(2)=j
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag5(i,j).eq.'n') then
                   endp(2)=endp(2)+1
                elseif (case_11.eq.1 .and. b(n).flag5(i,j).eq.'y') then
                   goto 500                   
                end if
             end do
             if (case_11.eq.1) goto 500       
          end do
          if (case_11.eq.0) goto 601
500       do m2=1,il(n)-i
             i3=i+m2
             do j3=start(2),endp(2)
                if (b(n).flag5(i3,j3).eq.'y') goto 5000
             end do
          end do
5000      endp(1)=endp(1)+m2-1
          do i=start(1),endp(1)
             do j=start(2),endp(2)
                b(n).flag5(i,j)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

           write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 50
c--------------------------------
601      nf=6

         bcdir='zta'
         start(3)=kl(n)
         endp(3)=kl(n)
60       case_11=0
         do i=1,il(n)
            do j=1,jl(n)
               if (case_11.eq.0 .and. b(n).flag6(i,j).eq.'n') then
                   start(1)=i
                   endp(1)=i
                   start(2)=j
                   endp(2)=j
                   case_11=1  
                elseif (case_11.eq.1 .and. b(n).flag6(i,j).eq.'n') then
                   endp(2)=endp(2)+1
                elseif (case_11.eq.1 .and. b(n).flag6(i,j).eq.'y') then
                   goto 600                   
                end if
             end do
             if (case_11.eq.1) goto 600       
          end do
          if (case_11.eq.0) goto 701
600       do m2=1,il(n)-i
             i3=i+m2
             do j3=start(2),endp(2)
                if (b(n).flag6(i3,j3).eq.'y') goto 6000
             end do
          end do
6000      endp(1)=endp(1)+m2-1
          do i=start(1),endp(1)
             do j=start(2),endp(2)
                b(n).flag6(i,j)='y'
             end do
          end do

          if(idimen.eq.2) then
             start(3) = 1
             endp(3) = 1
             istart(3) = 1
             iendp(3) = 1
          endif

           write(12,995) bcdir,n,start(1),start(2),start(3),
     $                    endp(1),endp(2),endp(3)
          goto 60
c--------------------------------
701       write(12,*)
          write(12,*)                                     
       end do  ! END of (n=1,NB)       
      write(12,996)

995   format("&bcdef bcdir=","'",a3,"'",", block=",
     > i3,", bctype=3, start=",
     > i3,",",i3,",",i3,", end=",i3,",",i3,",",i3,"/")

      
      
996   format("&bcdef bcdir='end'/")

998   format("&bcdef bcdir=","'",a3,"'",", block=",
     > i3,", bctype=7, start=",
     > i3,",",i3,",",i3,", end=",i3,",",i3,",",i3)
999   format(17x," ,iblock=",i3,",          istart=",
     >      i3,",",i3,",",i3,",iend=",i3,",",i3,",",i3,
     >        ", order=",i3,",",i3,",",i3,"/")

      
      end
      
            
                     
                 
 

