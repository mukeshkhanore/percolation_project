!khanore mueksh
!program to 2D site percolation in square lattice
!-------------------------------------------------------------------------------------    
    !      A: is 2d square arry ;					   	  n: is dimention of array A; 
    !     Rn:all row entries which are 1;				    	 Cn:all col entries which are 1;
    !	   p: probability at witch we finding cluster distribution  	   m:no's of entries are sweeped to 1 from 0 for given 'p'.
    !	   l: total sites						cluster(i):nos of lattice side which labeled by lable 'i'.
    ! cno(i):nos of "i" size clusters : contains all cluster sizes which are percolating
    !      s: average cluster size											
    module global
    implicit none
	integer, allocatable:: A(:,:),Cn(:),Rn(:),cluster(:),cno(:),pr(:),newseed(:)
    real*8,allocatable::sum1(:),avg_cluster(:)
	integer::n,k,i,j,chk,m,l,counter,yes,max11,perc,top,right,bottom,left,t1,r1,b,seed_size,y1,o,I1,J1
	real*8::r,p,start,final,avg_c,s,s1,q
    character(128) seeds
    endmodule global
!    
!------------main program-------------------------------------------------------------------------------------	
!
	program perc_path
    	use global 
!        
	write(*,*)"give the dimensions of array(square) n, and probability,p"
	read(*,*) n,p
    !call cpu_time(start)
    allocate(avg_cluster(10))
    l=n*n
    m=nint(p*l)    
       	allocate(A(n,n),pr(10))
CALL RANDOM_SEED(SIZE=seed_size) 
allocate(newseed(seed_size)) 
   write(*,*)"p=",p,"m=",m,"n=",n  
!        
        	 open(unit=100,file="avg_c.txt")
!                open(unit=110,file="time.txt")
                 open(unit=30,file="percolation_path.txt")
                 open(unit=18,file="avg.txt")
!do o=1,10        
	o=1
  write(seeds,*)o,"Seed.txt"
          open(unit=101,file=trim(Seeds))      
	  open(unit=10,file="row_column_no.txt")
          open(unit=12,file="matrix_mathematica.txt")
!         open(unit=120,file="matrix.txt")
	  open(unit=13,file="clusters and size.txt")
	  open(unit=101,file="seed.txt")
      counter=2
      chk=0
      b=0
      A=0
      k=0        
           CALL RANDOM_SEED ()  
           call random_number(r)
           r=r*50000
           y1=nint(r)
           write(101,*)y1
           close(unit=101)
	   open(unit=106,file="seed.txt",status="old",action="read")
 	   read(106,*)y1
	   newseed=y1
           close(unit=106)
	   CALL RANDOM_SEED(Put=newseed)  
 !-----------------------------------------------        
    	do i=1,n
      	do j=1,n
	CALL RANDOM_NUMBER(r) 
        if(r<p)then
          A(i,j)=1
          k=k+1
          write(10,*)i,j
        endif
            !  write(*,*)"inner",k
       	enddo 
       enddo
     close(unit=10)
     m=k
write(*,*)"for given p value sites are accupied now relabling them"
	allocate(sum1(m),Cn(m),Rn(m),cluster(m),Cno(m))
       cluster=0
            cno=0
             pr=0
           sum1=0
           s1=m
   open(unit=40,file="row_column_no.txt",status="old",action="read")        
    do i=1,m
      read(40,*)rn(i),cn(i)
      enddo 
      close(unit=40)      
!---------------lebeling started-------------------------------------------
	do
        	yes=1
            b=b+1
            write(*,*)"lebel b=",b
        		do i=1,m-1
    	    do j=i+1,m         
    		   if((cn(i)==cn(j)+1).and.(rn(i)==rn(j)))then!checking that downword neighbour is exist in this array
					call label()
                elseif((cn(i)==cn(j)-1).and.(rn(i)==rn(j)))then!checking that upword neighbour is exist in this array
                call label()
    		   elseif((cn(i)==cn(j)).and.(rn(i)==rn(j)+1))then!checking that forword neighbour is exist in this array
               	call label()
                elseif((cn(i)==cn(j)).and.(rn(i)==rn(j)-1))then!checking that backword neighbour is exist in this array
                call label()
              endif
           	enddo
         enddo
         if(yes==1)exit
	enddo
!------------------------sorting of clusters as per label---------------------------------------
write(*,*)"relabling is done now doing calculation"
do i=1,m
  k=0
  do j=1,m
    if(A(Rn(j),Cn(j))==i)then
		k=k+1
        cluster(i)=k
    endif    
  enddo
enddo    
!---------------------------calculating cnumber------------------------------------------
do i=2,m	!startin form 2 since all isolated sites are occupied by label 1 
  k=0
  if(cluster(i)>1)then !avoiding all isolated sites which are occupied
  do j=2,m
    if(cluster(i)==cluster(j))then
      k=k+1
      chk=cluster(i)
      cno(chk)=k
    endif        
  enddo
  endif
enddo    
!---------printing output----------------------------------------------------------------      
write(*,*)"calculations are done now printing output"
        write(12,*)"p=",p,"m=",m,"n=",n        
      	do i=1,n
      write(12,*)"{"
       write(12,*)(A(i,j),",",j=1,n-1)
       write(12,*)A(i,n)
          if(i<n)write(12,*)"},"
          if(i==n)write(12,*)"}"
        enddo
       close(unit=12)
!	do i=1,n
!	do j=1,n
!	write(120,"(100I5.1)")(A(i,j),j=1,n)
!	write(120,*)i,j,A(i,j)
!	enddo
!	enddo
!	close(unit=120)
        write(13,*)"	label	size of cluster	"        
  		do i=1,m
          if(cluster(i)>0)then
            write(13,*)i,cluster(i)
          endif      
         enddo         
!
        write(13,*)"	size  	nos of cluster	"        
	write(13,*)"m=",m
	write(13,*)1,cluster(1)
!
       		do i=2,m               
               if(cno(i)>0)then
         	write(13,*)i,cno(i)   
		if(i>=n)then
		max11=i
		call chk_path()
		endif
               endif      
         enddo		
      close(unit=13)
 !  write(*,*)"	size  	nos of cluster	"        
!do i=1,m
!write(*,*)i,cno(i)
!enddo
!-------------calculating avg cluster size-------------------------------         
			k=0
	do i=1,m
               if(i==1)then
                 sum1(i)=i*i*cluster(i)
               	 endif
                 if(cno(i)>0.and.i>1)then
                   sum1(i)=i*i*cno(i)
                 endif
            enddo
            s=sum(sum1)
         do i=1,m
           if(cno(i)>0)then
             k=k+1
           do j=1,10                       
             if(i==pr(j))then
		s=abs(s-i*i*cno(i))
                s1=abs(s1-i*cno(i))
                k=k-1
             endif
             enddo
             endif
             enddo       
           if(s1>0)then        
			q=s/s1     
            else
              q=0
              endif  
             avg_cluster(o)=q     
           write(100,*)"n,p,q,s,s1,y1",n,p,q,s,s1,y1
		deallocate(sum1,Cn,Rn,cluster,Cno)
!            enddo
            q=sum(avg_cluster)
            q=q!/10
	    i=1	
            !do i=1,10
             write(18,*)avg_cluster(i),q
            !enddo
             !   call cpu_time(final)                
            !write(110,*)"time",start,final
         end program	perc_path
!
!-----------------subroutine label-------------------------------------------------
!
	subroutine label()
    use global    
               	 if(A(rn(i),cn(i))==1.and.A(rn(j),cn(j))==1)then             
						A(rn(i),cn(i))=counter
			            A(rn(j),cn(j))=counter
           				counter=counter+1            	
				 endif
						if((A(Rn(i),Cn(i))==1).and.(A(Rn(j),Cn(j)))>1)then
                      		A(Rn(i),Cn(i))=A(Rn(j),Cn(j))
                      	endif
                      if((A(Rn(j),Cn(j))==1).and.(A(Rn(i),Cn(i)))>1)then
                      A(Rn(j),Cn(j))=A(Rn(i),Cn(i))
                      endif
                     if((A(Rn(j),Cn(j))>1).and.(A(Rn(i),Cn(i))>1))then                       
						if(A(Rn(i),Cn(i))<A(Rn(j),Cn(j)))then
                         A(Rn(j),Cn(j))=A(Rn(i),Cn(i))
                         yes=0
                         elseif(A(Rn(i),Cn(i))>A(Rn(j),Cn(j)))then
                           A(Rn(i),Cn(i))=A(Rn(j),Cn(j))
                           yes=0  
                           endif
                     endif     
                           end subroutine label
!
!--------------------------subroutine for checking percolation for cluster larger than "n"-------------------------------------------------
!
subroutine chk_path()
use global  
write(*,*)"max11",max11
		top=0
    		bottom=0
		perc=0
		right=0
	      left=0
    		k=0
		 t1=0
		  r1=0
		   chk=0
 	do i1=1,m		
            	if(cluster(i1)==MAX11)then
                chk=i1
    !            max11=cluster(i)
                endif
	enddo	
    	do j1=1,m
            if(A(rn(j1),cn(j1))==chk)then
            	if(rn(j1)==1)top=1
                if(rn(j1)==n)bottom=1
                  if(cn(j1)==1)left=1
                    if(cn(j1)==n)right=1
                      t1=top+bottom
                      r1=left+right
            endif
             	if((t1>1).or.(r1>1))exit
         enddo            
	      if((t1>1).or.(r1>1))then
        	perc=1
		pr(k+1)=max11
        	  else
        	perc=0
      		endif
	!if(perc>0)exit   
!         enddo			             
    	if(perc>0)then
		 k=k+1   
    	  write(30,*)o,"we found percolation path for n=",n,"p=",p,"cluster size=",max11,"cluster label=",chk
    	else
    	  write(30,*)o,"percolation path not found for n=",n,"p=",p,"cluster size=",max11,"cluster label=",chk
      endif  
    end subroutine        
