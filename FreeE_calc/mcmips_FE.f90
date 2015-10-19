! MC MIPs code

! started by tc387@cam.ac.uk in Autumn 2012
! more or less final form in Spring 2014
! specialised for Free Energy calculations in Summer 2015 


! CELL LISTS:  celllist(:,:,:,:) 1st dim tells how many (max mnpic) and which colloids are in the cell icol*icolspec, 2nd , 3rd and 4th dim are ii,jj,kk indeces of the cell 
! ipc(:,:,:) ith colloid cell - tells in which cell the given colloid is and on which position in a celllist chain in that cell. ipc(4,maxncol,ncolspec)
! rxcolgeo = rx specie on colloid geometry, last matrix in input file
!posolrx = position of colloid and rx that blong to it 
! bondrx(3,nrxpercol,maxncol,ncolspec) tells with who (colspecie,col,position) is particulat rx bound to , position=0 means anchors then (rxspecie,ianc,0)
program mcmips
implicit none
  
! TRY NOT USING GLOBAL VARIABLES BUT DEFINE EVERYTHING IN IT'S OWN SCOPE !
! USE INTENT IN/OUT IN ALL SUBROUTINES & FUNCTIONS (for easier debugging)!
! REMEMBER COLUMN MAJOR FORM !
integer :: ncellsc(3),ncellsa(3),mnpic,maxncol,nrxpercol,nfreecol,&
     ipcnew(3),icol,ii,jj,kk,jjcycle,excacc,moveacc,bondacc,&
     ncolspec,nrxspec,maxnanc,nimp_anc,icolspec
integer*8 :: icycle,ncycles,nout,totidm
real*8 :: lbox(3),tfinal,tout,rcol,socc(3),soca(3),rcut_bond,rnd,bond_ene,&
     time1,time2,tot_ene,fracgcmc,tot_ene_old,kspring,max_hop_col,fracbm,&
     dspring,max_rot_col,dimp_anc,avnbcoldble,KA_analyte
integer,allocatable :: celllistc(:,:,:,:),ipcc(:,:,:),celllista(:,:,:,:),ipca(:,:,:),&
     bondrx(:,:,:,:),bondanc(:,:,:),ncol(:),nanc(:),rxspec(:,:)
integer*8, allocatable :: avencol(:),avnbcol(:)
real*8, allocatable :: posanc(:,:,:),poscolrx(:,:,:,:),RXBONDENE(:,:),&
     activity(:),rxcolgeo(:,:,:)

logical :: random_seed=.true., read_init_conf=.false., mobile_anchors, &
     read_only_anchor_pos, imprint_anchors
real*8, parameter :: pi=3.14159265d0
character*50 :: outfilename, init_conf_filename

open(unit=101, file='input-par.dat',status='old')
read(101,*)
read(101,*)
read(101,*) (lbox(ii),ii=1,3)
read(101,*) ncolspec,nrxspec
allocate(ncol(ncolspec),activity(ncolspec),nanc(nrxspec))
read(101,*) (ncol(ii),ii=1,ncolspec)
read(101,*) maxncol
read(101,*) (nanc(ii),ii=1,nrxspec)
read(101,*) nrxpercol
read(101,*) kspring, dspring
read(101,*) (activity(ii),ii=1,ncolspec)
read(101,*) ncycles
read(101,*) nout
read(101,*) !---------------------------------------------
read(101,*) rcol
read(101,*) max_hop_col
read(101,*) max_rot_col
read(101,*) fracgcmc, fracbm
read(101,*) rcut_bond
read(101,*) mnpic
read(101,*) !---------------------------------------------
read(101,*) random_seed
read(101,*) mobile_anchors
read(101,*) imprint_anchors
read(101,*) read_init_conf, read_only_anchor_pos
read(101,*) init_conf_filename
read(101,*) outfilename
read(101,*) !--------------------------------------------
read(101,*) ! RX INTERACTION MATRIX
maxnanc=maxval(nanc(:))
allocate(RXBONDENE(nrxspec,nrxspec),rxcolgeo(4,nrxpercol,ncolspec),posanc(3,maxnanc,nrxspec))
do ii=1,nrxspec
   read(101,*) (RXBONDENE(jj,ii),jj=1,nrxspec)
enddo

read(101,*) !--------------------------------------------
read(101,*) ! RX BONDS ON COL GEOMETRY
do ii=1,ncolspec
   do jj=1,nrxpercol
      read(101,*) (rxcolgeo(kk,jj,ii),kk=1,4)
   enddo
enddo
read(101,*) !--------------------------------------------
read(101,*) ! IMPRINTED ANCHOr POSITIONS
if (imprint_anchors) then
   read(101,*) (((posanc(kk,jj,ii),kk=1,3),jj=1,nanc(ii)), ii=1,nrxspec)
endif
! define usefull parameters
activity(1:ncolspec)=exp(activity(1:ncolspec))
socc(:)=2*rcol
ncellsc(:)=floor(lbox(:)/socc(:))
ncellsc(:)=max(ncellsc,3)
socc(:)=lbox(:)/ncellsc(:) ! renormalize socc
ncellsa(:)=floor(lbox(:)/(2*rcol))
ncellsa(:)=max(ncellsa,3)
soca(:)=lbox(:)/ncellsa(:) ! renormalize soca
tot_ene=0

! rescale rx bond energy matrix with the unbound ligand (harmonic) partition function
RXBONDENE(:,:)=RXBONDENE(:,:) -1.5*log(kspring/2/pi)

! do boundary conditions on posanc
do ii=1,nrxspec
do jj=1,nanc(ii)
posanc(:,jj,ii)=posanc(:,jj,ii)-lbox(:)*floor(posanc(:,jj,ii)/lbox(:))
enddo
enddo

! write all parameters to screen
write(*,*)
write(*,*) '================ MC MIPS SIM ================'
write(*,*) '============================================='
write(*,*) '================= INPUT PAR ================='
write(*,*) 'boxsize  ' ,(lbox(ii),ii=1,3)
write(*,*) 'init # of colloids  ',ncol
write(*,*) 'max # of colloids  ',maxncol
write(*,*) '# of anchors(ligands)  ',nanc
write(*,*) 'fraction of insert/delete, bond create/destroy moves  ',fracgcmc,fracbm
write(*,*) 'activity  ', activity
write(*,*) 'tot n cycles  ',ncycles
write(*,*) 'nout   ',nout
write(*,*) '---------------------------------------------'
write(*,*) 'rcol  ', rcol
write(*,*) 'max hop/rot col  ',max_hop_col,max_rot_col
write(*,*) '# of rx sites per col ',nrxpercol
write(*,*) 'k spring ',kspring
write(*,*) 'single bond energy rxobondene(2,1)', RXBONDENE(2,1)
write(*,*) 'make/break bond distance cut off ',rcut_bond
write(*,*) 'max # of particles in each cell  ',mnpic
write(*,*) 'size of cell -colloids:', socc
write(*,*) 'size of cell -anchors:', soca
write(*,*) '---------------------------------------------'
write(*,*) 'random seed  ', random_seed
write(*,*) 'mobile anchors ',mobile_anchors
write(*,*) 'read initial conf', read_init_conf, read_only_anchor_pos
write(*,*) 'init conf filename  ', init_conf_filename
write(*,*) 'outfilename  ', outfilename
write(*,*) '============== END INPUT PAR ================'
write(*,*) '============================================='

! allocate big arrays and initialize to 0
allocate(celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol,ncolspec))
allocate(celllista(mnpic,ncellsa(1),ncellsa(2),ncellsa(3)),ipca(4,maxnanc,nrxspec))
allocate(bondrx(3,nrxpercol,&
     maxncol,ncolspec),bondanc(3,maxnanc,nrxspec),avencol(ncolspec),avnbcol(ncolspec))
allocate(poscolrx(3,1+nrxpercol,maxncol,ncolspec)) 
celllistc=0
celllista=0
bondrx=0
bondanc=0
avencol=0
avnbcol=0
poscolrx=0

! initial configuration set up
call initial_conf(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,nrxpercol,&
     kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,bondanc,poscolrx,&
     posanc,read_init_conf,random_seed,init_conf_filename,read_only_anchor_pos,&
     imprint_anchors,dimp_anc,nimp_anc)

call make_cell_list(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,nrxpercol,&
     mnpic,ncellsc,socc,celllistc,ipcc,ncellsa,soca,celllista,ipca,poscolrx,posanc)
                    
call cpu_time(time1)
      excacc=0
      moveacc=0
      bondacc=0
      totidm=0
      avencol=0
      avnbcol=0
write(*,*) '================= START SIM ================='
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXX  MAIN CYCLE  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
do icycle=1,ncycles

  ! do jjcycle=1, 2*sum(ncol)+1  ! so that there is approx 1 hop per colloid per cycle
      
      call random_number(rnd)
      if (rnd .gt. (fracgcmc+fracbm)) then  !HOP COL
         call mc_move(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,&
              nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,&
              bondanc,poscolrx,posanc,max_hop_col,max_rot_col,mobile_anchors,&
              moveacc,tot_ene)       
        
      elseif (rnd .lt. fracgcmc) then  ! INSERT/DELETE COL
         call insert_delete(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,&
              nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,&
              bondrx,bondanc,poscolrx,rxcolgeo,activity,excacc,tot_ene)       
         
      else   ! BOND CREATE/DESTROY
         call  bond_create_destroy(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,&
              rcol,nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,&
              bondanc,poscolrx,posanc,RXCOLGEO,RXBONDENE,bondacc,rcut_bond,tot_ene)  
      endif
      totidm=totidm+1
      avencol=avencol+sum(ncol)
      
      ! also get the average number of bound colloids
      do icolspec=1,ncolspec
         do icol=1,ncol(icolspec)
            if (any(bondrx(1,:,icol,icolspec) .gt. 0) ) then
               avnbcol=avnbcol+1
            endif
         enddo
      enddo
  ! enddo ! jjcycle

      ! bonds array chec for consistency: bondrx and bondanc
      ! only for debugging
  !    call check_bonds(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,nrxpercol,&
  !         mnpic,bondrx,bondanc,rcut_bond)
  !    call check_cell_list(ncolspec,ncol,maxncol,nrxpercol,lbox,poscolrx, &
  !         mnpic,ncellsc,socc,celllistc,ipcc)

   
   ! OUTPUT CONFIGURATION
   if (icycle/nout .eq. dble(icycle)/nout) then  
      call cpu_time(time2)
      tot_ene_old=tot_ene
      call tot_ene_calc(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,&
           rcol,nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,&
           bondanc,poscolrx,posanc,RXCOLGEO,RXBONDENE,tot_ene)
      avnbcoldble=dble(avnbcol(1))/dble(totidm) ! average number of bound colloids(analytes)
      write(*,*) 
      write(*,"(A,I10,A,F10.4)") 'icycle',icycle,', exe time (s) ',time2-time1
      write(*,"(A,F14.9,A,F14.9)") 'average avncol:',dble(avencol)/dble(totidm),&
           '  avnbcol:',avnbcoldble  
      write(*,"(A,F10.3,A,F10.3)") 'tot ene old:', tot_ene_old,'   tot ene:' , tot_ene
      write(*,"(A,F7.5,A,F7.5,A,F7.5)") 'moveacc:',real(moveacc)/nout/(real(2*&
           sum(avencol))/nout+1)/(1-fracgcmc-fracbm),'   excacc:',real(excacc)/nout/&
           (real(2*sum(avencol))/nout+1)/fracgcmc ,  '   bondacc:',real(bondacc)/&
           nout/(real(2*sum(avencol))/nout+1)/fracbm
      KA_analyte=(1.0d0/activity(1)+product(lbox(:)))*avnbcoldble/(1.0d0-avnbcoldble) 
      write(*,'(A,F10.4)') 'Analyte binding associaton constant K_A = ', KA_analyte
      write(*,'(A,F10.4)') 'Analyte binding Free Energy F_cav = ',-log(KA_analyte) 
      if (avnbcoldble .lt. 1d-5 ) write(*,*) 'WARNING, cavity occupancy very low, f=',&
           avnbcoldble ,'consider increasing the chemical potential'
     if (avnbcoldble .gt. 0.99 ) write(*,*) 'WARNING, cavity occupancy very high, f=',&
          avnbcoldble ,'consider decreasing the chemical potential'

      excacc=0
      moveacc=0
      bondacc=0
      totidm=0
      avencol=0
      avnbcol=0
      call output_conf(poscolrx,posanc,maxncol,maxnanc,ncol,ncolspec,nrxspec,&
           nrxpercol,nanc,activity,nout,icycle,bondanc,outfilename)
      time1=time2
   endif
  
enddo !ncycles

write(*,*) '==================== END ====================='
close(101)
close(909)

end program mcmips
! ========================================================================== !
! ========================================================================== !
! ========================================================================== !
subroutine initial_conf(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,nrxpercol,&
     kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,bondanc,poscolrx,&
     posanc,read_init_conf,random_seed,init_conf_filename,read_only_anchor_pos,& 
     imprint_anchors,dimp_anc,nimp_anc)
  implicit none
  integer, intent(inout) :: ncolspec,nrxspec,ncol(ncolspec),maxncol,nanc(nrxspec),&
       maxnanc,nrxpercol,mnpic,ncellsc(3),celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol,ncolspec),bondrx(3,nrxpercol,maxncol,ncolspec),bondanc(2,&
       maxval(nanc),nrxspec),nimp_anc
  real*8, intent(inout) :: lbox(3),rcol,kspring,dspring,socc(3),poscolrx(3,nrxpercol+1,&
       maxncol,ncolspec),posanc(3,maxnanc,nrxspec),dimp_anc
  logical, intent(in) :: read_init_conf,random_seed,read_only_anchor_pos,imprint_anchors
  character*50, intent(in) :: init_conf_filename 
  real*8, parameter :: pi=3.14159265d0
  integer :: ii,jj,kk,seed,now(3),iicol,icol,inputnanc
  real*8 :: rnd3(3),distij(3)
  logical :: overlap
  character*2 :: atom
 
  if (random_seed) then
     call itime(now)
     seed=10000*now(1)+100*now(2)+now(3)
     call seed_random_number
  else
     seed=10101
  endif
  !call ranz_set(seed)
  poscolrx(:,:,:,:)=0
  
  bondrx(:,:,:,:)=0
  bondanc(:,:,:)=0
  

  ! READ INITIAL CONF   ! NOT WORKING YET !!!!!!!!!!!!!!!!!!!
  if (read_init_conf) then
     if (read_only_anchor_pos) then ! only read anchors....
        open(unit=202, file=trim(init_conf_filename),status='old')
        read(202,*)
        read(202,*) ncol(1),inputnanc
        if(inputnanc.ne.nanc(1)) write(*,*) 'WARNING: INPUT NNANC .NE. NANC',inputnanc,nanc
!!$     if(ncol.gt.maxncol) write(*,*) 'WARNING: NCOL .GT. MAXNCOL',ncol,maxncol
        do ii=1,ncol(1)
           do jj=1,nrxpercol+1
              read(202,*) atom, (poscolrx(kk,jj,ii,1),kk=1,3)
           enddo
        enddo
        poscolrx=0 ! only read anchors this time
        ncol=0
        do ii=1,nanc(1)
           read(202,*) atom, (posanc(kk,ii,1),kk=1,3)
        enddo
        
     else ! read everything       
        write(*,*) 'NOt YET IMPLEMENTED!!!   exiting...... '
        call exit()
!!$     open(unit=202, file=trim(init_conf_filename),status='old')
!!$     read(202,*)
!!$     read(202,*) ncol,inputnanc
!!$     if(inputnanc.ne.nanc) write(*,*) 'WARNING: INPUT NNANC .NE. NANC',inputnanc,nanc
!!$     if(ncol.gt.maxncol) write(*,*) 'WARNING: NCOL .GT. MAXNCOL',ncol,maxncol
!!$     do ii=1,ncol
!!$        read(202,*) atom, (poscolrx(kk,ii),kk=1,3)
!!$     enddo
!!$     do ii=1,nanc
!!$        read(202,*) atom, (posancrx(kk,ii),kk=1,3)
!!$     enddo 
!!$     if(read_only_anchor_pos) then
!!$        poscol=0
!!$        ncol=0
!!$     endif
     endif ! read_only anchor_pos
  else ! INSERT STUFF
     if (imprint_anchors) then
!!$        if ((nanc(1) .ne.1).and.(nanc(2).ne.1)) then
!!$           write(*,*) 'WARNING!! IMPRINTING: r_jl=',dimp_anc,'   nanc(1:2) .ne. 1 '
!!$        endif
!!$        if (nimp_anc .eq. 2) then
!!$           posanc(:,1,1)=lbox(:)/2.0
!!$           posanc(:,1,2)=lbox(:)/2.0
!!$           posanc(1,1,2)=posanc(1,1,2)+dimp_anc
!!$           !  write(*,*) posanc   
!!$        elseif (nimp_anc .eq. 4) then
!!$           posanc(:,1,1)=lbox(:)/2.0
!!$           posanc(1,1,1)=posanc(1,1,1)-dimp_anc/2.0
!!$           posanc(:,1,2)=lbox(:)/2.0
!!$           posanc(2,1,2)=posanc(2,1,2)-dimp_anc/2.0
!!$           posanc(:,1,3)=lbox(:)/2.0
!!$           posanc(1,1,3)=posanc(1,1,3)+dimp_anc/2.0
!!$           posanc(:,1,4)=lbox(:)/2.0
!!$           posanc(2,1,4)=posanc(2,1,4)+dimp_anc/2.0
!!$        else
!!$           write(*,*) 'WARNING: imprinted anchors nimp_anc not recognized....',nimp_anc
!!$        endif
     else 
        posanc=0
        do ii=1,nrxspec
           do jj=1,nanc(ii)
              call random_number(rnd3)
              
              posanc(:,jj,ii)=rnd3(:)*lbox(:)
           enddo
        enddo
     endif
!!$     icol=1
!!$     do while (icol .le. ncol)
!!$        call random_number(rnd3)
!!$        poscol(:,icol)=rnd3(:)*lbox(:)
!!$        overlap=.false.
!!$        do iicol=1,icol-1
!!$           distij=poscol(:,icol)-poscol(:,iicol)
!!$           distij=distij-lbox*nint(distij/lbox)
!!$           if (sum(distij*distij) .lt. 4*rcol*rcol) then
!!$              overlap=.true.
!!$              exit
!!$           endif
!!$        enddo
!!$        if (.not. overlap) icol=icol+1
!!$     enddo
  endif
end subroutine initial_conf

! ========================================================================== !
subroutine insert_delete(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,&
     nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,bondanc,&
     poscolrx,rxcolgeo,activity,excacc,tot_ene)
  
  implicit none
  integer, intent(in) :: mnpic,ncellsc(3),maxncol,maxnanc,nrxpercol,ncolspec,&
       nrxspec,nanc(nrxspec)
  integer, intent(inout) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol,ncolspec),excacc,ncol(ncolspec),bondrx(3,nrxpercol,&
       maxncol,ncolspec),bondanc(3,maxnanc,nrxspec)
  real*8, intent(in) :: lbox(3),rcol,socc(3),activity(ncolspec),kspring,dspring,&
       rxcolgeo(4,nrxpercol,ncolspec)
  real*8, intent(inout) :: poscolrx(3,1+nrxpercol,maxncol,ncolspec),tot_ene
  real*8 :: Eold,Enew,rnd,rnd3(3),rinew(3),riold(3),arg,volume,ctb(3,nrxpercol),&
       S1rot,S2rot,V1rot(2),V2rot(2),ROTMATRIX(3,3),sqrtSrot,q(4),rnd2(2)
  integer :: kk,parti, inewcell(3),jj,ii,icol,ifreecol,iifreecol,ibond,ipcold(3),&
       icolspec,nfreecol
  logical :: overlap
  real*8, parameter :: pi=3.14159265d0, O=0.0d0 


  volume=product(lbox(:))  
  !O=0.0d0 ! define zero, for shorter matrix element writing
  icolspec=1
  if (ncolspec .gt. 1) then  ! choose the colloid specie
     call random_number(rnd)
     icolspec=floor(rnd*ncolspec)+1
  endif
  ! Get number of free colloids
  nfreecol=0
  do ii=1,ncol(icolspec)
     if (all(bondrx(1,:,ii,icolspec) .eq. 0)) nfreecol=nfreecol+1
  enddo
  
  call random_number(rnd)
  if (rnd .lt. 0.5d0) then    ! add particle 
     if (ncol(icolspec) .lt. maxncol) then
        call random_number(rnd3)
        rinew=rnd3(:)*lbox(:) 
        inewcell(:) = floor(rinew(:)/socc(:))+1 
        
        !CHECK FOR OVERLAP
        call col_overlap(overlap,icolspec,ncol(icolspec)+1,rinew,inewcell,mnpic,&
             celllistc,ipcc,ncellsc,lbox,poscolrx,rcol,maxncol,ncolspec,nrxpercol)
        if (.not. overlap) then
           arg=activity(icolspec)*volume/(nfreecol+1)
           call random_number(rnd)
           if (rnd .lt. arg) then
              excacc=excacc+1  ! acceptance rate
              ncol(icolspec)=ncol(icolspec)+1
              poscolrx(:,1,ncol(icolspec),icolspec)=rinew(:)           
              
              ipcc(1:3,ncol(icolspec),icolspec)=inewcell(:)
              celllistc(1,inewcell(1),inewcell(2),inewcell(3)) = &
                   celllistc(1,inewcell(1),inewcell(2),inewcell(3))+1
              ipcc(4,ncol(icolspec),icolspec)=celllistc(1,inewcell(1),&
                   inewcell(2),inewcell(3))
              celllistc(ipcc(4,ncol(icolspec),icolspec)+1,inewcell(1),&
                   inewcell(2),inewcell(3))=&
                   ncol(icolspec)+(icolspec-1)*maxncol

              bondrx(:,:,ncol(icolspec),icolspec)=0
              
              ! add rx binding sites
!!$              if (nrxpercol .eq. 6) then
!!$                 ctb(:,1)=(/-rcol,O,O/)
!!$                 ctb(:,2)=(/O,-rcol,O/)
!!$                 ctb(:,3)=(/rcol,O,O/)
!!$                 ctb(:,4)=(/O,rcol,O/)
!!$                 ctb(:,5)=(/O,O,rcol/)
!!$                 ctb(:,6)=(/O,O,-rcol/)
!!$              elseif (nrxpercol .eq. 4) then
!!$                 ctb(:,1)=(/-rcol,O,O/)
!!$                 ctb(:,2)=(/O,-rcol,O/)
!!$                 ctb(:,3)=(/rcol,O,O/)
!!$                 ctb(:,4)=(/O,rcol,O/)
!!$              elseif  (nrxpercol .eq. 2) then
!!$                 ctb(:,1)=(/-rcol,O,O/)
!!$                 ctb(:,2)=(/rcol,O,O/)
!!$              else
!!$                 write(*,*) 'WARNING: n of rx bind sites per colloid not supported',nrxpercol
!!$                 call exit(1802)
!!$              endif              
              ctb(:,:)=rcol*rxcolgeo(2:4,:,icolspec)
              
              ! GET RANDOM QUATERNION
              S1rot=2
              do while (S1rot .gt. 1.0)
                 call random_number(rnd2)
                 V1rot=2*rnd2-1.0d0
                 S1rot=sum(V1rot*V1rot)
              enddo
              S2rot=2
              do while (S2rot .gt. 1.0)
                 call random_number(rnd2)
                 V2rot=2*rnd2-1.0d0
                 S2rot=sum(V2rot*V2rot)
              enddo
              sqrtSrot=sqrt((1-S1rot)/S2rot)
              ! random rotation quaternion and matrix
              q=(/V1rot(1),V1rot(2),V2rot(1)*sqrtSrot,V2rot(2)*sqrtSrot/)
              ROTMATRIX(1,1:3)=(/q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4),&
                   2*q(2)*q(3)-2*q(1)*q(4),2*q(2)*q(4)+2*q(1)*q(3)/)
              ROTMATRIX(2,1:3)=(/2*q(2)*q(3)+2*q(1)*q(4),q(1)*q(1)-q(2)*q(2)+&
                   q(3)*q(3)-q(4)*q(4),2*q(3)*q(4)-2*q(1)*q(2)/)
              ROTMATRIX(3,1:3)=(/2*q(2)*q(4)-2*q(1)*q(3),2*q(3)*q(4)+&
                   2*q(1)*q(2),q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)/)
              ! rotate rx bind sites
              do ii=1,nrxpercol
                 ctb(:,ii)=matmul(ROTMATRIX,ctb(:,ii))
                 poscolrx(:,ii+1,ncol(icolspec),icolspec)=poscolrx(:,1,&
                      ncol(icolspec),icolspec)+ctb(:,ii)
              enddo              
           endif    
        endif ! NOT OVERLAP
     endif
    
     
     ! DELETE PARTICLE
  elseif (nfreecol .gt. 0)   then   
     call random_number(rnd)
    
     ifreecol=floor(nfreecol*rnd)+1
     ! find real icol
     icol=1
     iifreecol=0
     do while (.true.)
        if (all(bondrx(1,:,icol,icolspec) .eq. 0)) iifreecol=iifreecol+1
        if (iifreecol .eq. ifreecol) exit 
        icol=icol+1
     enddo
     if (icol .gt. ncol(icolspec)) write(*,*) 'WARNING: icol > ncol while trying to delete a particle ..', icol,ncol(icolspec) 
     if (any(bondrx(1,:,icol,icolspec) .ne. 0)) write(*,*) 'WARNING: TRY TO DELETE COL WITH RX' 
     arg=nfreecol/(activity(icolspec))/volume
     
     call random_number(rnd)
     if (rnd .lt. arg) then     
        excacc=excacc+1   ! update acceptnce ratio
        
        
        if (icol .ne. ncol(icolspec)) then ! swap colloid indexes
           poscolrx(:,:,icol,icolspec)=poscolrx(:,:,ncol(icolspec),icolspec)
           bondrx(:,:,icol,icolspec)=bondrx(:,:,ncol(icolspec),icolspec)
           ! update bond bookkeeping 
           do ii=1,nrxpercol
              if (bondrx(1,ii,ncol(icolspec),icolspec) .gt. 0) then ! it is bound    
                 if (bondrx(3,ii,ncol(icolspec),icolspec) .eq. 0) then  ! it is bound to anchor                   
                    bondanc(2,bondrx(2,ii,ncol(icolspec),icolspec),bondrx(1,ii,ncol(icolspec),icolspec))&
                         =icol        
                 else
                    bondrx(2,bondrx(3,ii,icol,icolspec),bondrx(2,ii,icol,icolspec),&
                         bondrx(1,ii,icol,icolspec))=icol
                    write(*,*) 'WARNING: colloids should only bind to anchors, not to each other... (while deleting)'
                 endif
              endif
           enddo
        endif ! icol .ne. ncol
        
        bondrx(:,:,ncol(icolspec),icolspec)=0
        
        !update celllists
        jj=ipcc(4,icol,icolspec)
        ipcold(:)=ipcc(1:3,icol,icolspec)   
        
        celllistc(jj+1, ipcold(1), ipcold(2), ipcold(3)) = &
             celllistc(celllistc(1,ipcold(1),ipcold(2),ipcold(3))+1, ipcold(1),ipcold(2),ipcold(3))
        ipcc(4,celllistc(jj+1,ipcold(1),ipcold(2),ipcold(3)),icolspec) = jj   
        
        celllistc(1,ipcold(1),ipcold(2),ipcold(3)) = celllistc(1,ipcold(1),&
             ipcold(2),ipcold(3))-1 
        celllistc(celllistc(1,ipcold(1),ipcold(2),ipcold(3))+2,ipcold(1),&
             ipcold(2),ipcold(3))=0
        
        ipcc(:,icol,icolspec)=ipcc(:,ncol(icolspec),icolspec)
        if (icol .ne. ncol(icolspec)) then
           celllistc(1+ipcc(4,icol,icolspec),ipcc(1,icol,icolspec),ipcc(2,icol,&
                icolspec),ipcc(3,icol,icolspec))=icol+maxncol*(icolspec-1)
        endif
        ipcc(:,ncol(icolspec),icolspec)=0

        poscolrx(:,:,ncol(icolspec),icolspec)=0
        ncol(icolspec)=ncol(icolspec)-1
     endif
  endif
end subroutine insert_delete
!=============================================================!

subroutine mc_move(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,nrxpercol,kspring,&
     dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,bondanc,poscolrx,posanc,&
     max_hop_col,max_rot_col,mobile_anchors,moveacc,tot_ene)
 
  implicit none
  integer, intent(in) :: ncolspec,nrxspec,ncol(ncolspec),maxncol,nanc(nrxspec),&
       maxnanc,nrxpercol,mnpic,ncellsc(3),bondrx(3,nrxpercol,maxncol,ncolspec),&
       bondanc(3,maxnanc,nrxspec)
  real*8, intent(in) :: lbox(3),rcol,kspring,dspring,socc(3),max_hop_col,max_rot_col
  integer, intent(inout) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol,ncolspec),moveacc
 
  real*8, intent(inout) :: poscolrx(3,1+nrxpercol,maxncol,ncolspec),&
       posanc(3,maxnanc,nrxspec),tot_ene
  logical, intent(in) :: mobile_anchors
  real*8 :: Eold,Enew,rnd1,rnd3(3),rinew(3,1+nrxpercol),riold(3,1+nrxpercol),arg,&
       distij(3),globalarg,hopdist(3),rnd2(2),rnd4(4),qrot(4),Vrot(2),Srot,urot(3),&
       sqrtSrot,sinfi,cosfi,omcosfi,ROTMATRIX(3,3),globalrnd,ctb(3,nrxpercol),q(4)
  integer :: parti, ipcnew(3) ,icol,ibond,ianc,jj,icolspec,ii,iancspec
  logical :: overlap

  globalrnd=0
  globalarg=0.6  ! also move anchors
  if (mobile_anchors) call random_number(globalrnd)
  if (globalrnd .lt. globalarg) then ! MOVE COLLOIDS
  call random_number(rnd2)
  icolspec=floor(ncolspec*rnd2(1))+1 ! random specie
  if (ncol(icolspec) .lt. 1) RETURN   ! return if no colloids
  icol=floor(rnd2(2)*ncol(icolspec))+1  ! random particle
  riold(:,:)=poscolrx(:,:,icol,icolspec)

  call random_number(rnd3) 
  do ii=1,nrxpercol+1
     rinew(:,ii)=riold(:,ii)+(2*rnd3(:)-1.0d0)*max_hop_col
  enddo
  rinew(:,1)=rinew(:,1)-lbox(:)*floor(rinew(:,1)/lbox(:))
  ipcnew=floor(rinew(:,1)/socc)+1  
  call col_overlap(overlap,icolspec,icol,rinew(:,1),ipcnew,mnpic,celllistc,ipcc,&
       ncellsc,lbox,poscolrx,rcol,maxncol,ncolspec,nrxpercol)
  if (.not. overlap) then
     ! ROTATE WITHOUT QUATERNIONS
     ! MARSAGLIA METHOD TO GET RANDOM AXIS, THEN ROTATE BY SMALL RANDOM ANGLE
     Srot=2
     do while (Srot .ge. 1.0)
        call random_number(rnd3)
        Vrot=2*rnd3(1:2)-1.0d0
        Srot=sum(Vrot*Vrot)
     enddo
     sqrtSrot=sqrt(1-Srot)
     urot=(/2*Vrot(1)*sqrtSrot,2*Vrot(2)*sqrtSrot,1-2*Srot/) ! random rotation axis
     sinfi=max_rot_col*(2*rnd3(3)-1.0)
     cosfi=dsqrt(1.0d0-sinfi*sinfi)
     omcosfi=1.0d0-cosfi

     ! ROTATION MATRIX
     ROTMATRIX(1,1:3)=(/cosfi+urot(1)*urot(1)*omcosfi,urot(1)*urot(2)*omcosfi-&
          urot(3)*sinfi,urot(1)*urot(3)*omcosfi+urot(2)*sinfi/)
     ROTMATRIX(2,1:3)=(/urot(2)*urot(1)*omcosfi+urot(3)*sinfi,cosfi+urot(2)*&
          urot(2)*omcosfi,urot(2)*urot(3)*omcosfi-urot(1)*sinfi/)
     ROTMATRIX(3,1:3)=(/urot(3)*urot(1)*omcosfi-urot(2)*sinfi,urot(3)*urot(2)*&
          omcosfi+urot(1)*sinfi,cosfi+urot(3)*urot(3)*omcosfi/)
 !    q=(/cosfi,sinfi*(/urot(1),urot(2),urot(3)/)/)
 !    ROTMATRIX(1,1:3)=(/q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4),&
 !         2*q(2)*q(3)-2*q(1)*q(4),2*q(2)*q(4)+2*q(1)*q(3)/)
 !    ROTMATRIX(2,1:3)=(/2*q(2)*q(3)+2*q(1)*q(4),q(1)*q(1)-q(2)*q(2)+&
 !         q(3)*q(3)-q(4)*q(4),2*q(3)*q(4)-2*q(1)*q(2)/)
 !    ROTMATRIX(3,1:3)=(/2*q(2)*q(4)-2*q(1)*q(3),2*q(3)*q(4)+&
 !         2*q(1)*q(2),q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)/)

     ! GET RELATIVE POSITION VECTORS, BIND SITE TO PARTICLE CENTRE     
     do ii=1,nrxpercol
        ctb(:,ii)=rinew(:,ii+1)-rinew(:,1)
        ctb(:,ii)=ctb(:,ii)-lbox*nint(ctb(:,ii)/lbox)
        ctb(:,ii)=matmul(ROTMATRIX(:,:),ctb(:,ii))
        rinew(:,ii+1)=rinew(:,1)+ctb(:,ii) ! NEW POSITION
     enddo
  !   write(*,*) '================================================='
  !   write(*,*) urot,ROTMATRIX,ctb
     ! CALCULATE OLD & NEW ENERGY
     Eold=0
     Enew=0
     if (any(bondrx(1,:,icol,icolspec) .gt. 0)) then ! CALLCULATE RX BONDING ENERGY
        do ibond=1,nrxpercol
           if (bondrx(1,ibond,icol,icolspec) .gt. 0) then
              if (bondrx(3,ibond,icol,icolspec) .gt. 0) then !attached to other particle
                 distij=riold(:,ibond+1)-poscolrx(:,bondrx(3,ibond,icol,icolspec),&
                      bondrx(2,ibond,icol,icolspec),bondrx(1,ibond,icol,icolspec))
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 Eold=Eold+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2
                 
                 distij=rinew(:,ibond+1)-poscolrx(:,bondrx(2,ibond,icol,icolspec),&
                      bondrx(2,ibond,icol,icolspec),bondrx(1,ibond,icol,icolspec))
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 Enew=Enew+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2              
              else ! attached to anchor
                 distij=riold(:,ibond+1)-posanc(:,bondrx(2,ibond,icol,icolspec),&
                      bondrx(1,ibond,icol,icolspec))
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 Eold=Eold+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2
                 
                 distij=rinew(:,ibond+1)-posanc(:,bondrx(2,ibond,icol,icolspec),&
                      bondrx(1,ibond,icol,icolspec))
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 Enew=Enew+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2
              endif
           endif
        enddo
     endif
     
     if (Enew .le. Eold) then 
        moveacc=moveacc+1      ! update acceptance ratio
        tot_ene=tot_ene+Enew-Eold   
        poscolrx(:,:,icol,icolspec)=rinew(:,:)
      
        if (any(ipcnew .ne. ipcc(1:3,icol,icolspec))) then
           call update_cell_list(icol,maxncol,ipcnew,ipcc,celllistc,&
                ncellsc,mnpic,icolspec,ncolspec)        
        endif
     else     
        call random_number(rnd1)
        if (rnd1 .lt. exp(Eold-Enew)) then
           moveacc=moveacc+1      ! update acceptance ratio
           tot_ene=tot_ene+Enew-Eold   
           poscolrx(:,:,icol,icolspec)=rinew(:,:)
        
           if (any(ipcnew .ne. ipcc(1:3,icol,icolspec))) then
              call update_cell_list(icol,maxncol,ipcnew,ipcc,celllistc,&
                   ncellsc,mnpic,icolspec,ncolspec)     
           endif
        endif
     endif
  endif ! NOT OVERLAP 
  
  ! MOVE ANCHORS
else 
   call random_number(rnd2)
   iancspec=floor(nrxspec*rnd2(1))+1 ! random specie
   if (nanc(iancspec) .lt. 1) RETURN   ! return if no anchors
   ianc=floor(rnd2(2)*nanc(iancspec))+1  ! random particle
   riold(:,1)=posanc(:,ianc,iancspec)
   
   call random_number(rnd3)
   rinew(:,1)=riold(:,1)+2*(2*rnd3-1.0d0)*max_hop_col
   rinew(:,1)=rinew(:,1)-lbox(:)*floor(rinew(:,1)/lbox(:))
   
   Eold=0
   Enew=0
   if (bondanc(1,ianc,iancspec) .gt. 0) then ! ene calc
      distij=riold(:,1)-poscolrx(:,bondanc(3,ianc,iancspec)+1,&
           bondanc(2,ianc,iancspec),bondanc(1,ianc,iancspec))
      distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
      Eold=Eold+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2
      
      distij=rinew(:,1)-poscolrx(:,bondanc(3,ianc,iancspec)+1,&
           bondanc(2,ianc,iancspec),bondanc(1,ianc,iancspec))
      distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
      Enew=Enew+kspring/2.0*(sqrt(sum(distij*distij))-dspring)**2 
   endif
   if (Enew .lt. Eold) then
      moveacc=moveacc+1      ! update acceptance ratio
      tot_ene=tot_ene+Enew-Eold
      posanc(:,ianc,iancspec)=rinew(:,1)
   else
      call random_number(rnd1)
      if (rnd1 .lt. exp(Eold-Enew)) then
         moveacc=moveacc+1      ! update acceptance ratio
         tot_ene=tot_ene+Enew-Eold
         posanc(:,ianc,iancspec)=rinew(:,1)
      endif
   endif
endif

end subroutine mc_move
!===========================================================================!
subroutine col_overlap(overlap,icolspec,icol,rinew,icell,mnpic,celllistc,ipcc,&
     ncellsc,lbox,poscolrx,rcol,maxncol,ncolspec,nrxpercol)
  implicit none
  
  integer, intent(in) :: icolspec,ncolspec,icol,mnpic,icell(3),ncellsc(3),&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),maxncol,&
       ipcc(4,maxncol,ncolspec),nrxpercol
  real*8, intent(in) :: rinew(3), poscolrx(3,nrxpercol+1,maxncol,ncolspec),rcol,lbox(3)
  logical, intent(out) :: overlap
  real*8 ::  distij(3), distij2, distij1,sigma2
  integer :: jj, kk, xc, yc, zc, cxc, cyc, czc,jcol,jspecie
 
  sigma2=4*rcol*rcol
  overlap=.false.
  do zc=icell(3)-1, icell(3)+1
     czc=zc-ncellsc(3)*floor(real(zc-1)/ncellsc(3))
     do yc=icell(2)-1, icell(2)+1
        cyc=yc-ncellsc(2)*floor(real(yc-1)/ncellsc(2))
        do xc=icell(1)-1, icell(1)+1
           cxc=xc-ncellsc(1)*floor(real(xc-1)/ncellsc(1))
          
           do jj=1, celllistc(1,cxc, cyc, czc)  
              if (celllistc(jj+1,cxc,cyc,czc) .ne. (icol+maxncol*(icolspec-1))) then
                 
                 jspecie= ceiling((celllistc(jj+1,cxc,cyc,czc)-1e-10)/dble(maxncol))
                 jcol=(celllistc(jj+1,cxc,cyc,czc)-(jspecie-1)*maxncol)
                                  
                 distij=rinew(:)-poscolrx(:,1,jcol,jspecie)
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 distij2=sum(distij*distij)
                 
                 if (distij2 .lt. sigma2) then                 
                    overlap= .true.
                    return            ! OVERLAP
                 endif
                 
              endif            
           enddo
          
        enddo
     enddo
  enddo 
end subroutine col_overlap
!==========================================================================!
subroutine bond_create_destroy(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,&
     rcol,nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,&
     bondanc,poscolrx,posanc,RXCOLGEO,RXBONDENE,bondacc,rcut_bond,tot_ene)
  implicit none
  integer, intent(in) :: ncolspec,nrxspec,ncol(ncolspec),maxncol,nanc(nrxspec),&
       maxnanc,nrxpercol,mnpic,ncellsc(3),celllistc(mnpic,ncellsc(1),ncellsc(2),&
       ncellsc(3)),ipcc(3,maxncol,ncolspec)
  real*8, intent(in) :: lbox(3),rcol,kspring,dspring,socc(3),rcut_bond,&
       poscolrx(3,1+nrxpercol,maxncol,ncolspec), posanc(3,maxnanc,nrxspec),&
       RXBONDENE(nrxspec,nrxspec),RXCOLGEO(4,nrxpercol,ncolspec)
  integer, intent(inout) ::  bondanc(3,maxnanc,nrxspec),bondrx(3,nrxpercol,&
       maxncol,ncolspec),bondacc
  real*8, intent(inout) :: tot_ene
  integer :: ii,jj,kk,icol,ibond,ianc,xc,yc,zc,cxc,cyc,czc,cnearbond(3,1000),inearb,&
       nnearb,jcellcol(3),nncl,ioldcol,icell(3),ispecie,jspecie,jcol,jrx,&
       joldspecie,joldcol,joldrx
  real*8 :: rnd1,distij(3),ij_ene,lambda,arg,qq,pnearbond(1000),enearbond(1000),pibond
  logical :: breakbond
  
  nncl=ceiling(rcut_bond/minval(socc))  ! number of near cell layers: depth bond moves in cells
  nncl=min(nncl,floor(minval(lbox)/2))
  pnearbond=0 ! probability to form a bond with particular bond on particular colloid
  enearbond=0 ! energy of that bond
  cnearbond=0 ! gives global serial number of colloid
  call random_number(rnd1)
  ispecie=floor(rnd1*nrxspec)+1
  if (nanc(ispecie) .lt. 1) RETURN
  call random_number(rnd1)
  ianc=floor(rnd1*nanc(ispecie))+1
  icell=floor(posanc(:,ianc,ispecie)/socc(:))
  
  ! CHECK If BONDED WITH OTHER THAN NEAR NEIGHBOUR
  if (bondanc(1,ianc,ispecie) .gt. 0) then
     jcellcol=floor(poscolrx(:,1,bondanc(2,ianc,ispecie),bondanc(1,ianc,ispecie))&
          /socc(:))+1
     if (any(abs((icell(:)-jcellcol(:))-ncellsc(:)*nint(real(icell(:)-jcellcol(:))/&
          ncellsc(:))).gt.nncl))  RETURN
  endif   

  inearb=0 ! number of free nearby bonding sites
  do zc=icell(3)-nncl, icell(3)+nncl
     czc=zc-ncellsc(3)*floor(real(zc-1)/ncellsc(3))
     do yc=icell(2)-nncl, icell(2)+nncl
        cyc=yc-ncellsc(2)*floor(real(yc-1)/ncellsc(2))
        do xc=icell(1)-nncl, icell(1)+nncl
           cxc=xc-ncellsc(1)*floor(real(xc-1)/ncellsc(1))          
           do jj=1, celllistc(1,cxc, cyc, czc)  
              jspecie= ceiling((celllistc(jj+1,cxc,cyc,czc)-1e-10)/dble(maxncol))
              jcol=celllistc(1+jj,cxc, cyc, czc)-(jspecie-1)*maxncol 
              do jrx=1,nrxpercol
                 ! if not bound and binding energy less than 50kT
                 if ((bondrx(1,jrx,jcol,jspecie) .eq. 0) .and. &
                      (RXBONDENE(ispecie,nint(RXCOLGEO(1,jrx,jspecie))) .lt. 50)) then
                    inearb=inearb+1
                    cnearbond(:,inearb)=(/jspecie,jcol,jrx/)
                    distij=posanc(:,ianc,ispecie)-poscolrx(:,1+jrx,jcol,jspecie)
                    distij=distij-lbox*nint(distij/lbox) ! periodic boundary
                    ij_ene=RXBONDENE(ispecie,nint(RXCOLGEO(1,jrx,jspecie)))+kspring/2*&
                         (sqrt(sum(distij*distij))-dspring)**2
                    enearbond(inearb)=ij_ene
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  jspecie=0
  jcol=0
  jrx=0
  if (bondanc(1,ianc,ispecie).gt.0) then  ! ADD BOND
     inearb=inearb+1   
     jspecie=bondanc(1,ianc,ispecie)
     jcol=bondanc(2,ianc,ispecie)
     jrx=bondanc(3,ianc,ispecie)
     cnearbond(:,inearb)=(/jspecie,jcol,jrx/)
     distij=posanc(:,ianc,ispecie)-poscolrx(:,1+jrx,jcol,jspecie)
     distij=distij-lbox*nint(distij/lbox) ! periodic boundary
     ij_ene=RXBONDENE(ispecie,nint(RXCOLGEO(1,jrx,jspecie)))+kspring/2*(sqrt(sum(distij*&
          distij))-dspring)**2
     enearbond(inearb)=ij_ene
  endif

  nnearb=inearb
  if (nnearb .eq. 0) RETURN   ! IF NO NEIGHBOUR BONDS THEN RETURN
  ! normalize probability
  pnearbond(1:nnearb)=exp(-enearbond(1:nnearb))
  pnearbond(1:nnearb)=pnearbond(1:nnearb)/(1+sum(pnearbond(1:nnearb)))  
  inearb=1
  breakbond=.false.
  call random_number(rnd1)
  pibond=pnearbond(1)
  do while (pibond .lt. rnd1)
     if (inearb .ge. nnearb) then
        breakbond=.true.
        exit
     endif
     inearb=inearb+1
     pibond=pibond + pnearbond(inearb)
  enddo

  if (breakbond) then ! breakbond (if it existed)
     if (bondanc(1,ianc,ispecie) .gt. 0) then
     !   jspecie=bondnac(1,ianc,ispecie) ! STILL ON MEMORY FROM BEFORE
     !   jcol=bondnac(2,ianc,ispecie)
     !   jrx=bondnac(3,ianc,ispecie)
      
        bondrx(:,jrx,jcol,jspecie)=0
        bondanc(:,ianc,ispecie)=0
        tot_ene=tot_ene-enearbond(nnearb)  ! last colloid is the one bound to anchor
        bondacc=bondacc+1
     endif
  else  ! form/swap bond
   
     jspecie=cnearbond(1,inearb) !new specie to bind to
     jcol=cnearbond(2,inearb) ! new colloid
     jrx=cnearbond(3,inearb) ! new rx bond
   
     if (any((/jspecie,jcol,jrx/) .ne. bondanc(:,ianc,ispecie))) then !swap/create bond
        if (bondanc(1,ianc,ispecie) .gt. 0) then ! delete old bond         
           joldspecie=bondanc(1,ianc,ispecie)
           joldcol=bondanc(2,ianc,ispecie)
           joldrx=bondanc(3,ianc,ispecie)
           bondanc(:,ianc,ispecie)=0
           bondrx(:,joldrx,joldcol,joldspecie)=0
           tot_ene=tot_ene-enearbond(nnearb)  ! last colloid is the one bound to anchor
        endif
        ! form new bond
        bondanc(:,ianc,ispecie)=(/jspecie,jcol,jrx/)
        bondrx(:,jrx,jcol,jspecie)=(/ispecie,ianc,0/)
        tot_ene=tot_ene+enearbond(inearb)
        bondacc=bondacc+1
     endif 
  endif

end subroutine bond_create_destroy
!==========================================================================!
subroutine tot_ene_calc(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,&
     rcol,nrxpercol,kspring,dspring,mnpic,ncellsc,socc,celllistc,ipcc,bondrx,&
     bondanc,poscolrx,posanc,RXCOLGEO,RXBONDENE,tot_ene)
  implicit none
  integer, intent(in) :: ncolspec,nrxspec,ncol(ncolspec),maxncol,nanc(nrxspec),&
       maxnanc,nrxpercol,mnpic,ncellsc(3),celllistc(mnpic,ncellsc(1),ncellsc(2),&
       ncellsc(3)),ipcc(3,maxncol,ncolspec),&
       bondanc(3,maxnanc,nrxspec),bondrx(3,nrxpercol,maxncol,ncolspec)
  real*8, intent(in) :: lbox(3),rcol,kspring,dspring,socc(3),&
       poscolrx(3,1+nrxpercol,maxncol,ncolspec), posanc(3,maxnanc,nrxspec),&
       RXBONDENE(nrxspec,nrxspec),RXCOLGEO(4,nrxpercol,ncolspec)
  real*8, intent(out) :: tot_ene
  real*8 ::  distij(3), distij2
  integer :: ii,jj,kk,ianc,irx,icol,ispecie
  
  tot_ene=0
  do ispecie=1,ncolspec
     do icol=1,ncol(ispecie)
        do irx=1,nrxpercol
           if (bondrx(1,irx,icol,ispecie) .gt. 0) then
              if (bondrx(3,irx,icol,ispecie) .eq. 0) then ! bound to anchor..
                 distij(:)=posanc(:,bondrx(2,irx,icol,ispecie),bondrx(1,irx,icol,&
                      ispecie))-poscolrx(:,irx+1,icol,ispecie)
                 distij=distij-lbox*nint(distij/lbox)
                 distij2=sum(distij*distij)
                 tot_ene=tot_ene+RXBONDENE(bondrx(1,irx,icol,ispecie),nint(&
                      RXCOLGEO(1,irx,ispecie)))+kspring/2*(sqrt(distij2)-dspring)**2
              else
                 write(*,*) 'WARNING: BUG WITH BONDRX/BONDANC '
              endif
           endif
        enddo
     enddo
  enddo
  
end subroutine tot_ene_calc
! ========================================================================== !
SUBROUTINE seed_random_number()
  implicit none
  ! Local variables
  INTEGER              :: ii,k, now(3)
  INTEGER, ALLOCATABLE :: seed(:)
  CALL RANDOM_SEED(SIZE=k)
  ALLOCATE( seed(k) )
  call itime(now)
  do ii=1, k
     seed(ii)=now(mod(ii+1,3)+1)
  end do
  CALL RANDOM_SEED(PUT=seed)
  DEALLOCATE( seed )
  RETURN
END SUBROUTINE seed_random_number
! ========================================================================== !
subroutine make_cell_list(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,rcol,&
     nrxpercol,mnpic,ncellsc,socc,celllistc,ipcc,ncellsa,soca,celllista,ipca,&
     poscolrx,posanc)
  implicit none

  integer, intent(in) :: ncellsc(3),ncellsa(3),ncolspec,nrxspec,ncol(ncolspec),&
       maxncol,mnpic,nanc(nrxspec),nrxpercol,maxnanc
  real*8, intent(in) :: socc(3),soca(3),poscolrx(3,nrxpercol+1,maxncol,ncolspec),&
       posanc(3,maxval(nanc),nrxspec),lbox(3),rcol
  integer, intent(out) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)), ipcc(4,maxncol,ncolspec)
  integer, intent(out) :: celllista(mnpic,ncellsa(1),ncellsa(2),ncellsa(3)), ipca(4,maxval(nanc),nrxspec)
  integer ::  ii,jj,kk,nblob
  
  celllistc(:,:,:,:)=0
  ipcc(:,:,:)=0
  celllista(:,:,:,:)=0
  ipca(:,:,:)=0
  
  do ii=1,ncolspec
     do jj=1, ncol(ii)
        ipcc(1:3,jj,ii)=floor(poscolrx(:,1,jj,ii)/socc(:))+1
     enddo
  enddo
  do ii=1,ncolspec
     do jj=1,ncol(ii)
        celllistc(1, ipcc(1,jj,ii), ipcc(2,jj,ii), ipcc(3,jj,ii)) = & 
             celllistc(1, ipcc(1,jj,ii), ipcc(2,jj,ii), ipcc(3,jj,ii)) + 1
        celllistc((celllistc(1,ipcc(1,jj,ii),ipcc(2,jj,ii),ipcc(3,jj,ii))+1),&
             ipcc(1,jj,ii),ipcc(2,jj,ii), ipcc(3,jj,ii))=jj+(ii-1)*maxncol !NOTE jj+ii*maxncol
        
        ipcc(4,jj,ii)=celllistc(1, ipcc(1,jj,ii), ipcc(2,jj,ii), ipcc(3,jj,ii))
     enddo
  enddo
end subroutine make_cell_list
! ========================================================================== ! 
subroutine update_cell_list(parti,maxnpart,ipcnew,ipc,celllist,ncells,mnpic,&
     ispecie,nspecies)
  implicit none
  integer, intent(in) ::  parti,maxnpart,mnpic,ncells(3),ipcnew(3),ispecie,nspecies
  integer, intent(inout) :: ipc(4,maxnpart,nspecies), celllist(mnpic,ncells(1),&
       ncells(2),ncells(3))
  integer :: jj, ipcold(3)
  jj=ipc(4, parti,ispecie)
  ipcold(:)=ipc(1:3,parti,ispecie)   
  
  celllist(jj+1, ipcold(1), ipcold(2), ipcold(3)) = &
       celllist(celllist(1,ipcold(1),ipcold(2),ipcold(3))+1, ipcold(1),ipcold(2),ipcold(3))
  ipc(4,celllist(jj+1,ipcold(1),ipcold(2),ipcold(3))-maxnpart*(ispecie-1),ispecie) = jj   
  
  celllist(1,ipcold(1),ipcold(2),ipcold(3)) = celllist(1,ipcold(1),ipcold(2),ipcold(3))-1 
  celllist(celllist(1,ipcold(1),ipcold(2),ipcold(3))+2,ipcold(1),ipcold(2),ipcold(3))=0   

  celllist(celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))+2, ipcnew(1),ipcnew(2),ipcnew(3)) = parti+(ispecie-1)*maxnpart
  celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))=celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))+1
  
  ipc(1:3,parti,ispecie)=ipcnew(:)
  ipc(4, parti,ispecie)=celllist(1,ipcnew(1),  ipcnew(2), ipcnew(3))
end subroutine update_cell_list
! ========================================================================== !
subroutine check_bonds(lbox,ncolspec,nrxspec,ncol,maxncol,nanc,maxnanc,nrxpercol,&
     mnpic,bondrx,bondanc,rcut_bond)
  implicit none
  integer, intent(in) :: ncolspec,nrxspec,ncol(ncolspec),maxncol,nanc(nrxspec),&
       maxnanc,nrxpercol,mnpic,bondanc(3,maxnanc,nrxspec),bondrx(3,nrxpercol,maxncol,ncolspec)
  real*8, intent(in) :: lbox(3),rcut_bond
  integer :: ii,jj,kk,ianc,icol,ispecie,irx,jcolspec,jancspec,&
       jcol,jrx,janc,kancspec,kcolspec,kanc,kcol,krx
  real*8 :: rcut_bond2

  rcut_bond2=rcut_bond**2

  ! do anchor_colrxr bond check on anchor side
  do ispecie=1,nrxspec
     do ianc=1,nanc(ispecie)
        if (bondanc(1,ianc,ispecie) .gt. 0) then
           jcolspec=bondanc(1,ianc,ispecie)
           jcol=bondanc(2,ianc,ispecie)
           jrx=bondanc(3,ianc,ispecie)
           
           kancspec=bondrx(1,jrx,jcol,jcolspec)
           kanc=bondrx(2,jrx,jcol,jcolspec)
           
           if ( any((/ispecie,ianc/) .ne. (/kancspec,kanc/))) then
              write(*,*)
              write(*,*) 'anchor side bonding error!!  anc i,speci  points to rx j in colloid j,specj :', &
                   ianc,ispecie, jrx,jcol,jcolspec
              write(*,*) ' but rx j in  col j,specj  points to anc k, speck : ', kanc,kancspec


              write(*,*) 'bondanc: '
              do janc=1,nanc(1)
                 write(*,*) bondanc(:,janc,1)
              enddo
              write(*,*) 'bondrx: '
              do jcol=1,ncol(1)+1
                 write(*,*) bondrx(:,:,jcol,1)
              enddo
              call exit(0)
           endif

           if (jcol .gt. ncol(jcolspec)) then
              write(*,*)  'anchor side bonding error!!  anc i,speci  points to rx j in colloid j,specj :', &
                   ianc,ispecie, jrx,jcol,jcolspec
              write(*,*) ' but ncol= ',ncol(jcolspec)
           endif

        endif
     enddo
  enddo
  ! do lig-receptor bond check on ligand side
  do ispecie=1,ncolspec
     do icol=1,maxncol
        do irx=1,nrxpercol
           jancspec=bondrx(1,irx,icol,ispecie)
           janc=bondrx(2,irx,icol,ispecie)
           if (janc .gt. 0) then
              kcolspec=bondanc(1,janc,jancspec)
              kcol=bondanc(2,janc,jancspec)
              krx=bondanc(3,janc,jancspec)
              
              if (any((/kcolspec,kcol,krx/) .ne. (/ispecie,icol,irx/) )) then
                 write(*,*) 'rx side bonding error!!  rx i on col i points to anc j :', irx,icol,janc
                 write(*,*) ' but anc j points to rx k on col k : ',krx,kcol
              endif
           endif
        enddo
     enddo
  enddo
end subroutine check_bonds
! ========================================================================== !
subroutine check_cell_list(ncolspec,ncol,maxncol,nrxpercol,lbox,poscolrx, &
     mnpic,ncellsc,socc,celllistc,ipcc)
  implicit none  
  integer, intent(in) :: ncolspec,ncol(ncolspec),maxncol,nrxpercol,mnpic,&
       ncellsc(3),celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol,ncolspec)
  real*8, intent(in) :: lbox(3),socc(3),poscolrx(3,nrxpercol+1,maxncol,ncolspec)
  integer :: ii,jj,kk,icol,icolspec,jcolspec,jcol,ipctrial(3),ncolincell,&
       icell,jcell,kcell
  
  do icolspec=1,ncolspec
     do icol=1,ncol(icolspec)
        ipctrial(1:3)=floor(poscolrx(:,1,icol,icolspec)/socc(:))+1
        
        if (any(ipctrial(1:3) .ne. ipcc(1:3,icol,icolspec))) then
           write(*,*) 'WARNING: check_cell_list failed..... ipcc not good with poscolrx'
           write(*,*) 'ipcc= ',ipcc(1:3,icol,icolspec), 'ipctrial = ', ipctrial(1:3)
           write(*,*) 'poscolrx(:,1,icol,icolspec) ', poscolrx(:,1,icol,icolspec)
        endif
        
        ! check celllist array
        if ( celllistc(1+ipcc(4,icol,icolspec), ipcc(1,icol,icolspec), ipcc(2,icol,icolspec), &
             ipcc(3,icol,icolspec)) .ne. icol+(icolspec-1)*maxncol) then
           write(*,*)'WARNING: check_cell_list failed..... '
           write(*,*)  'celllist in cell gives particle', celllistc(1+ipcc(4,icol,icolspec), ipcc(1,icol,icolspec), &
                ipcc(2,icol,icolspec), ipcc(3,icol,icolspec)) ,' but particle is ',  icol+(icolspec-1)*maxncol
        endif
     enddo
  enddo
  
  do icell=1,ncellsc(1)
     do jcell=1,ncellsc(2)
        do kcell=1,ncellsc(3)
           ncolincell=celllistc(1, icell,jcell,kcell)
           do ii=2,ncolincell+1            
              jcol=celllistc(ii, icell,jcell,kcell) ! col serial number, convert to ispecie and icol
              icolspec=ceiling(dble(jcol)/maxncol)
              icol=jcol-(icolspec-1)*maxncol
              
              if (any(ipcc(1:3,icol,icolspec).ne.(/icell,jcell,kcell/))) then
                 write(*,*) 'WARNING: check_cell_list failed.....  '
                 write(*,*)  'ipcc gives cell ',ipcc(1:3,icol,icolspec)  ,' but cell is ',&
                      (/icell,jcell,kcell/)
           
              endif
              if (ipcc(4,icol,icolspec) .ne. ii-1)  then
                 write(*,*)'WARNING: check_cell_list failed..... ',&
                   'particle position in given cell array is wrong'
               
              endif
           enddo
           do ii=ncolincell+2,mnpic
              if (celllistc(ii, icell,jcell,kcell) .ne. 0)  write(*,*)'WARNING: check_cell_list failed..... ',&
                   'cell list has non-zero entries beyond the number of blobs in the cell'
           enddo
        enddo
     enddo
  enddo
end subroutine check_cell_list
! ========================================================================== !
 subroutine output_conf(poscolrx,posanc,maxncol,maxnanc,ncol,ncolspec,nrxspec,&
      nrxpercol,nanc,activity,nout,icycle,bondanc,outfilename)
  implicit none

  integer, intent(in) ::  maxncol,maxnanc,ncolspec,nrxspec, &
       ncol(ncolspec),nanc(nrxspec),nrxpercol,bondanc(3,maxnanc,nrxspec)
  integer*8, intent(in) :: nout, icycle
  real*8, intent(in) :: poscolrx(3,nrxpercol+1,maxncol,ncolspec), &
       posanc(3,maxnanc,nrxspec),activity(ncolspec)
  character*50, intent(in) :: outfilename
  character*50 :: filename
  character*4 :: char1,char2,char3,char4
  integer :: ianc,icol,inputnanc,ispecie,ii
  ! VMD OUTPUT FILE
  write(char1,'(I4)') nanc(1)+1000
  write(char2,'(I4)') 2000+nint(log(activity(1)))
  write(char3,'(I4)') ncol(1)+1000
  write(char4,'(I4)') int(icycle/nout)+1000
  !write(*,*) 'OK'
  filename=trim(outfilename)//'_mu'//char2//'_ncol'//char3//'_out'//char4//'.xyz'
  
  open(unit=222,file=filename)
  write(222,*) sum(ncol)*(nrxpercol+1)+sum(nanc)
  write(222,*) ncol, sum(nanc), 'sandwich brushes MD'
  
  
  do icol=1,ncol(1)         
     write(222,"(A,F8.3,F8.3,F8.3)") 'Xe  ', poscolrx(:,1,icol,1)
     do ii=2,nrxpercol+1
        write(222,"(A,F8.3,F8.3,F8.3)") 'H  ', poscolrx(:,ii,icol,1)
     enddo
  enddo
  if (ncolspec .gt. 1) then
     do icol=1,ncol(2)         
        write(222,"(A,F8.3,F8.3,F8.3)") 'I  ', poscolrx(:,1,icol,2) 
        do ii=2,nrxpercol+1
           write(222,"(A,F8.3,F8.3,F8.3)") 'H  ', poscolrx(:,ii,icol,2)
        enddo
     enddo
  endif
  do ispecie=1,nrxspec
     do ianc=1,nanc(ispecie)
        if (bondanc(1,ianc,ispecie) .gt. 0) then ! ligand is bound
           write(222,"(A,F8.3,F8.3,F8.3)") 'N   ', posanc(:,ianc,ispecie)
        else 
           write(222,"(A,F8.3,F8.3,F8.3)") 'O   ', posanc(:,ianc,ispecie) 
           
        endif
     enddo
  enddo
  close(222)
  ! END VMD OUTPUT FILE
  
end subroutine output_conf
