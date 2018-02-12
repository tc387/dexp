! MC Hairy Nanoparticle code, aka multi-multivalency + mobile ligands + 

! by tc387@cam in Spring 2015
! updated October 2016, Free energy as a function of receptor composition
 
! CELL LISTS:  celllist(:,:,:,:) 1st dim tells how many (max mnpic) and which col are in the cell icol*icolspec, 2nd , 3rd and 4th dim are ii,jj,kk indeces of the cell 
! ipc(:,:,:) ith colloid cell - tells in which cell the given colloid is and on which position in a celllist chain in that cell. ipc(4,maxncol,ncolspec)
! rxcolgeo = rx specie on colloid geometry, last matrix in input file
!posolrx = position of colloid and rx that blong to it 
! bondrx(3,nrxpercol,maxncol,ncolspec) tells with who (specie,col,position) is particulat rx bound to , position=0 means anchors then(rxspecie,iancm,0)
program mchnp
implicit none
  
! TRY NOT USING GLOBAL VARIABLES BUT DEFINE EVERYTHING IN IT'S OWN SCOPE !
! USE INTENT IN/OUT IN ALL SUBROUTINES & FUNCTIONS (for easier debugging)!
! REMEMBER COLUMN MAJOR FORM !
integer :: ncellsc(3),ncellsb(3),ncellsr(2),mnpic,maxncol,nchainspercol,nfreecol,&
     ipcnew(3),icol,ii,jj,kk,jjcycle,excacc,moveacc,bondacc,ncol,nrec,&
     nimp_anc,nblobsperchain,nblob,iiblob,jjblob,kfactor,maxnblob,nrxspec,&
     nWL_points,ibin,n1rWL_max,n1rWL_min,WLacc,nbonds
integer*8 :: icycle,ncycles,nout,totidm,avencol,avnbpc,avnbpci,avnbcol,avnbcolii,avnbpcind
real*8 :: lbox(3),rcol,socc(3),socb(3),socr(2),rcut_ligrec,rcut_blobwall,&
     rcut_blobblob,rcut_global,rnd,rnd2(2),rnd3(3),ligrec_ene,rrec,&
     time1,time2,tot_ene,fracgcmc,tot_ene_old,max_hop_col,fracbm,&
     max_rot_col,max_hop_blob,dimp_anc,mu,activity,cWL,cWL_max,cWL_min,fWL_start,&
     fWL_stop,fWL_reduce,fWL,hWL_conv,max_colzWL
integer, allocatable :: celllistc(:,:,:,:),ipcc(:,:),celllistb(:,:,:,:),ipcb(:,:),&
     celllistr(:,:,:),ipcr(:,:),ligboundto(:,:),recboundto(:,:),nrecs(:),&
     irecspec(:),ichainspec(:),nchainsperspec(:),spectorec(:,:)
integer*8, allocatable :: hWL(:)
real*8, allocatable :: posrec(:,:),posblob(:,:,:,:),poscol(:,:),rxcolgeo(:,:),&
     RXBONDENE(:,:),psiWL(:),histWL(:)
logical, allocatable :: tmparr(:)  ! logical tells if a given rec is bound
logical :: random_seed=.true., read_init_conf=.false., mobile_receptors, & 
     mobile_ligands,read_only_anchor_pos, random_anconcolpos, doWL_flag
real*8, parameter :: pi=3.14159265d0
character*50 :: outfilename, init_conf_filename

open(unit=101, file='input-par.dat',status='old')
read(101,*)
read(101,*)
read(101,*) (lbox(ii),ii=1,3)
read(101,*) nrxspec
allocate(nrecs(nrxspec),nchainsperspec(nrxspec))
read(101,*) (nrecs(ii),ii=1,nrxspec)
read(101,*) maxncol
read(101,*) nchainspercol
allocate(ichainspec(nchainspercol))
read(101,*) (nchainsperspec(ii),ii=1,nrxspec)
kk=0 ! assign species to chains
do ii=1,nrxspec
   do jj=1,nchainsperspec(ii)
      kk=kk+1
      ichainspec(kk)=ii
   enddo
enddo
if (sum(nchainsperspec(:)) .gt. nchainspercol) then !check
   write(*,*) 'WARNING: sum(nchainsperspec) .gt.  nchainspercol,   exiting... '
   call exit()
endif
read(101,*) nblobsperchain
read(101,*) mu
read(101,*) ncycles
read(101,*) nout
read(101,*) !---------------------------------------------
read(101,*) rcol
read(101,*) rrec
read(101,*) rcut_ligrec
read(101,*) !---------------------------------------------
read(101,*) fracgcmc, fracbm
read(101,*) max_hop_blob,max_hop_col, max_rot_col
read(101,*) kfactor
read(101,*) mnpic
read(101,*) random_seed
read(101,*) mobile_receptors,mobile_ligands
read(101,*) random_anconcolpos
read(101,*) read_init_conf
read(101,*) init_conf_filename
read(101,*) outfilename
read(101,*) !--------------------------------------------
read(101,*) max_colzWL
read(101,*) fWL_start,fWL_stop,fWL_reduce
read(101,*) hWL_conv
read(101,*) !--------------------------------------------
read(101,*) ! RX INTERACTION MATRIX
allocate(RXBONDENE(nrxspec,nrxspec))
do ii=1,nrxspec
   read(101,*) (RXBONDENE(jj,ii),jj=1,nrxspec)
enddo
allocate(rxcolgeo(3,nchainspercol))
if (.not. random_anconcolpos) then ! read anc on col geometry
   read(101,*) !--------------------------------------------
   read(101,*) ! RX BONDS ON COL GEOMETRY
   do jj=1,nchainspercol
      read(101,*) ichainspec(jj), (rxcolgeo(kk,jj),kk=1,3)
   enddo
endif

! define usefull parameters
nrec=sum(nrecs(:)) ! total number of receptors
activity=exp(mu)
maxnblob=maxncol*nchainspercol*nblobsperchain ! maximum possible number of blobs in the system
rcut_blobblob=3.0d0
rcut_blobwall=2.0d0
rcut_global=max(rcut_blobblob,rcut_ligrec)

socc(:)=2*rcol
ncellsc(:)=floor(lbox(:)/socc(:))
ncellsc=max(ncellsc,3)
socc(:)=lbox(:)/ncellsc(:) ! renormalize socc

ncellsb=floor(lbox(:)/rcut_blobblob)
ncellsb=max(ncellsb,3)
socb(:)=lbox(:)/ncellsb(:) ! renormalize socb

ncellsr(:)=floor(lbox(1:2)/(rcut_ligrec))
ncellsr=max(ncellsr,3)
socr(:)=lbox(1:2)/ncellsr(:) ! renormalize soca

tot_ene=0
excacc=0
moveacc=0
bondacc=0
totidm=0
avencol=0
avnbpcind=0
avnbpc=0
avnbpci=0
avnbcol=0
avnbcolii=0

ncol=0
nblob=0

  
! write all parameters to screen
write(*,*)
write(*,*) '================ MC MIPS SIM ================'
write(*,*) '============================================='
write(*,*) '================= INPUT PAR ================='
write(*,*) 'boxsize  ' ,(lbox(ii),ii=1,3)
write(*,*) '# of rx species  ', nrxspec
write(*,*) '# of receptos per rx specie  ', nrecs
write(*,*) 'max # of colloids  ', maxncol
write(*,*) '# of chains per collod ',nchainspercol
write(*,*) 'rx specie per chain ', ichainspec
write(*,*) '# of blobs per chain ', nblobsperchain
write(*,*) 'mu  ', mu
write(*,*) 'tot n cycles  ',ncycles
write(*,*) 'nout   ',nout
write(*,*) '---------------------------------------------'
write(*,*) 'col radius rcol  ', rcol
write(*,*) 'receptor radius rrec  ', rrec
write(*,*) 'bond distance cut off ',rcut_ligrec
write(*,*) '---------------------------------------------'
write(*,*) 'fraction of insert/delete, bond create/destroy moves  ',fracgcmc,fracbm
write(*,*) 'max hop blob, col, rot col  ',max_hop_blob,max_hop_col,max_rot_col
write(*,*) '# of Rosenbluth trials -- kfactor ', kfactor
write(*,*) 'max # of particles in each cell+1  ',mnpic
write(*,*) 'size of cell -colloids:', socc
write(*,*) 'size of cell -receptors:', socr
write(*,*) 'size of cell -blobs:', socb
write(*,*) '---------------------------------------------'
write(*,*) 'random seed  ', random_seed
write(*,*) 'mobile receptors, mobile ligands ',mobile_receptors, mobile_ligands 
write(*,*) 'random anchors on collod ',random_anconcolpos
write(*,*) 'read initial conf', read_init_conf
write(*,*) 'init conf filename  ', init_conf_filename
write(*,*) 'outfilename  ', outfilename
write(*,*) '---------------------------------------------'
write(*,*) '  max colloid height for WL  ',max_colzWL
write(*,*) 'f parameters for WL psi update ', fWL_start, fWL_stop, fWL_reduce  
write(*,*) 'WL Histogram convergance criterion ',  hWL_conv
write(*,*) '---------------------------------------------'
! ! RX INTERACTION MATRIX
write(*,*) 'RX interaction matrix'
do ii=1,nrxspec
   write(*,*) (RXBONDENE(jj,ii),jj=1,nrxspec)
enddo
if (.not. random_anconcolpos) then ! write anc on col geometry
   write(*,*) 'RX on colloid geometry: specie, xyz pos'
   do jj=1,nchainspercol
      write(*,*) ichainspec(jj), (rxcolgeo(kk,jj),kk=1,3)
   enddo
endif
write(*,*) '============== END INPUT PAR ================'
write(*,*) '============================================='

! allocate big arrays
allocate(celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol))
allocate(celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob))
allocate(celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec))
allocate(poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol),&
     posrec(2,nrec)) 
allocate(ligboundto(nchainspercol,maxncol),recboundto(2,nrec))
allocate(tmparr(nrec),irecspec(nrec))
allocate(psiWL(nchainspercol+1),hWL(nchainspercol+1),spectorec(nrec,nrxspec))

ligboundto(:,:)=0
recboundto(:,:)=0
!------  WANG LANDAU INIT  ------!
fwl=fwl_start
psiWL(:)=0
hWL(:)=0
doWL_flag=.true.
WLacc=0

call initial_conf(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,nrec,&
     rcol,rrec,poscol,posblob,posrec,lbox,mnpic,nrxspec,nrecs,irecspec,ichainspec,&
     spectorec,random_seed,read_init_conf,init_conf_filename)

call make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,nrec,&
     poscol,posblob,posrec,lbox,&
     celllistc,ipcc,celllistb,ipcb,celllistr,ipcr,ncellsc,ncellsb,ncellsr,socc,socb,socr,mnpic)

call insert_col_WL(lbox,ncol,maxncol,nchainspercol,nblobsperchain,nrec,nrecs,maxnblob,rcol,&
     poscol,posblob,posrec,mnpic,ncellsc,ncellsb,ncellsr,socc,socb,socr,celllistc,celllistb,&
     celllistr,ipcc,ipcb,ipcr,kfactor,ligboundto,recboundto,rcut_blobblob,rcut_blobwall,rcut_ligrec,&
     rxcolgeo,random_anconcolpos,nrxspec,irecspec,ichainspec,RXBONDENE,tot_ene)
nbonds=1  ! the number of formed bonds is 1 after WL col insert

!call make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,nrec,&
!     poscol,posblob,posrec,lbox,&
!     celllistc,ipcc,celllistb,ipcb,celllistr,ipcr,ncellsc,ncellsb,ncellsr,socc,socb,socr,mnpic)
!write(*,*) '============== MADE CELL LIST ==============='

call cpu_time(time1)

write(*,*) '================= START SIM ================='
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXX  MAIN CYCLE  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
do icycle=1,ncycles
   do jjcycle=1, nchainspercol*nblobsperchain+1  ! so that there is approx 1 hop per colloid per cycle
      
      call random_number(rnd)
      if (rnd .gt. (fracgcmc+fracbm)) then  !HOP COL/BLOB
         call mc_move_wl(lbox,ncol,maxncol,nchainspercol,nblobsperchain,nrec,maxnblob,rcol,&
              poscol,posblob,posrec,mnpic,ncellsc,ncellsb,ncellsr,socc,socb,socr,celllistc,&
              celllistb,celllistr,ipcc,ipcb,ipcr,ligboundto,recboundto,rcut_blobblob,&
              rcut_blobwall,rcut_ligrec,ligrec_ene,max_hop_blob,max_hop_col,max_rot_col,&
              moveacc,tot_ene,RXBONDENE,nrxspec,irecspec,ichainspec,mobile_ligands,psiWL,&
              max_colzWL,nbonds)
         
         if (mobile_receptors) then
            ! move receptors, 1 move for each blob move..
            call hop_rec(posrec,poscol,lbox,nrec,ncol,rrec,socr,mnpic,recboundto,&
                 ligboundto,posblob,nchainspercol,maxncol,nblobsperchain,maxnblob,celllistr,&
                 ipcr,ncellsr,rcut_ligrec,rcut_blobblob,rcut_blobwall,max_hop_blob,tot_ene)          
         endif

      elseif (rnd .lt. fracgcmc) then  ! INSERT/DELETE COL
      !   call insert_delete(lbox,ncol,maxncol,nchainspercol,nblobsperchain,nrec,maxnblob,rcol,&
      !        poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,celllistc,celllistb,ipcc,ipcb,kfactor,ligboundto,&
      !        recboundto,rcut_blobblob,rcut_blobwall,rxcolgeo,activity,excacc,random_anconcolpos,tot_ene)
         
      else   ! BOND CREATE/DESTROY
         write(*,*) 'bond create/destroy not yet implemented for WL... exiting...'
         call exit()
         call bond_create_destroy(posrec,posblob,lbox,nrec,ncol,maxncol,nchainspercol,& 
              nblobsperchain,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,bondacc,&
              ligboundto,recboundto,nrxspec,RXBONDENE,irecspec,ichainspec)

      endif
      totidm=totidm+1 ! total number of moves
      avencol=avencol+ncol ! average number of colloids
      

      ! check bonds -- only for debugging
      call check_bonds(nrec,ncol,maxncol,maxnblob,nchainspercol,& 
           nblobsperchain,ligboundto,recboundto,posblob,posrec,rcut_ligrec,lbox)   
      call check_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,lbox,poscol,&
           posblob,mnpic,ncellsc,socc,celllistc,ipcc,ncellsb,socb,celllistb,ipcb)  
   enddo ! jjcycle

   ! do Wang Landau
   call WangLandau_bonds(nchainspercol,nbonds,psiWL,&
     hWL,fWL,fWL_stop,fWL_reduce,hWL_conv,icycle,WLacc,doWL_flag)
   
   ! get the average number of bound ligands per chain
   ! only every 10 cycles, because it's expensive
   if (icycle/10 .eq. dble(icycle)/10) then
      
      ! get the average number of bound colloids
      avnbcolii=0
      do ii=1,ncol
         if (any(ligboundto(:,ii).ne.0)) then
            avnbcolii=avnbcolii+1
         endif
      enddo
      avnbcol=avnbcol+avnbcolii
      
      ! get average number of bonds per bound colloid
      tmparr=(recboundto(2,:).ne.0)      
      avnbpci=0
      do ii=1,nrec
         if (tmparr(ii)) avnbpci=avnbpci+1
      enddo
      if (avnbcolii .gt. 0) then
         avnbpc=avnbpc+avnbpci/dble(avnbcolii)
         avnbpcind=avnbpcind+1
      endif
   endif


!!$   if ((any(ligboundto .gt. 0)).and.(poscol(3,1).gt.  max_colzWL)) then
!!$      write(*,*) ' bonded and col pos .gt. h_0 '
!!$      write(*,*) poscol(3,1), '   ', sum(ligboundto)
!!$   elseif  ((sum(ligboundto) .eq. 0).and.(poscol(3,1).gt.  max_colzWL)) then
!!$      write(*,*) 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
!!$      write(*,*) ' unbonded and col pos .gt. h_0 '
!!$      write(*,*) poscol(3,1), '   ', sum(ligboundto)
!!$   endif
   
   
   ! OUTPUT CONFIGURATION
   if (icycle/nout .eq. dble(icycle)/nout) then  
      call cpu_time(time2)
      tot_ene_old=tot_ene

      call tot_ene_calc(posblob,poscol,lbox,nrec,ncol,maxncol,nchainspercol,& 
           nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
           socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,tot_ene)
      write(*,*)
      write(*,"(A,I8,A,F10.4)") 'icycle',icycle,', exe time (s) ',time2-time1
      write(*,"(A,F11.6)")  'av # bonds_per_col:',dble(avnbpc)/dble(avnbpcind)
      write(*,"(A,F9.3,A,F9.3)")  '   tot ene old:', tot_ene_old,'   tot ene new:' , tot_ene
      write(*,"(A,F7.5,A,F7.5,A,F7.5)") 'moveacc:',dble(moveacc)/totidm/&
           (1-fracgcmc-fracbm),'   excacc:',real(excacc)/dble(totidm)/&
           fracgcmc ,  '   bondacc:',dble(bondacc)/dble(totidm)/fracbm
      write(*,"(A,F9.7,A,F6.2,A,F5.3,A,F6.4)") 'fWL: ',fWL,'   hWLrat: ', &
           dble(maxval(hWL))/minval(hWL),'   comp: ',dble(nrecs(1))/nrec, '  WLacc: ',dble(WLacc)/nout
     ! write(*,*) psiWL
     ! write(*,*) hWL
      excacc=0
      moveacc=0
      bondacc=0
      totidm=0
      avencol=0
      avnbcol=0
      avnbpc=0
      avnbpcind=0
      WLacc=0
      call output_conf_WL(poscol,posblob,posrec,maxncol,nchainspercol,nblobsperchain,ncol,&
           nrec,recboundto,activity,nout,icycle,outfilename,nrxspec,irecspec,ichainspec,&
           cWL_min,cWL_max,psiWL,hWL,doWL_flag)
      time1=time2
   endif
  
enddo !ncycles

write(*,*) '==================== END ====================='
close(101)
!close(909)

end program mchnp
! ========================================================================== !
! ========================================================================== !
! ========================================================================== !
subroutine initial_conf(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,nrec,&
     rcol,rrec,poscol,posblob,posrec,lbox,mnpic,nrxspec,nrecs,irecspec,ichainspec,&
     spectorec,random_seed,read_init_conf,init_conf_filename)

  implicit none
  integer, intent(inout) :: ncol,maxncol,maxnblob,nrec,nchainspercol,&
       nblobsperchain,mnpic,nrxspec,nrecs(nrxspec),ichainspec(nchainspercol),&
       irecspec(nrec),spectorec(nrec,nrxspec)
  real*8, intent(inout) :: lbox(3),rcol,rrec,poscol(3,maxncol),posrec(2,nrec),&
       posblob(3,nblobsperchain,nchainspercol,maxncol)
  logical, intent(in) :: read_init_conf,random_seed
  character*50, intent(in) :: init_conf_filename 
  real*8, parameter :: pi=3.14159265d0
  integer :: ii,jj,kk,seed,now(3),irec,jrec,nrectmp,ispec,nrecstmp(nrxspec)
  real*8 :: rnd3(3),distij(3),distij2,rnd2(2)
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
  poscol(:,:)=0
  posblob(:,:,:,:)=0
  posrec(:,:)=0

  ! READ INITIAL CONF   ! NOT WORKING YET !!!!!!!!!!!!!!!!!!!
  if (read_init_conf) then
     write(*,*) 'NOt YET IMPLEMENTED!!!   exiting...... '
     call exit()
     
  else ! INSERT STUFF
     if (ncol .gt. 0) then
        write(*,*) 'WARNING: particles can only be inserted via gcmc, exiting ...'
        call exit()
     endif
     

     spectorec(:,:)=0
     irecspec(:)=0
     ! insert receptors
     nrecstmp(:)=nrecs(:)
     irec=1
     do ispec=1,nrxspec
        ii=1
        do while (ii .le. nrecs(ispec))
           irecspec(irec)=ispec
           spectorec(ii,ispec)=irec
    !       write(*,*) ii,irec,ispec,irecspec(irec)
           call random_number(rnd2)
           posrec(1:2,irec)=rnd2(:)*lbox(1:2)
           
           nrectmp=irec
           if (irec .gt. 1) then
              do jrec=1,nrectmp-1
                 distij(:)=0
                 distij(1:2)=posrec(:,nrectmp)-posrec(:,jrec)
                 distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                 distij2=sum(distij*distij)
                 if (distij2 .lt. 4*rrec*rrec) then
                    irec=irec-1
                    ii=ii-1
                    exit
                    
                 endif
              enddo
           endif
                   
           irec=irec+1
           ii=ii+1
         
        enddo
     enddo
  endif

end subroutine initial_conf

! ========================================================================== !
! ========================================================================== !
subroutine insert_col_WL(lbox,ncol,maxncol,nchainspercol,nblobsperchain,nrec,nrecs,maxnblob,rcol,&
     poscol,posblob,posrec,mnpic,ncellsc,ncellsb,ncellsr,socc,socb,socr,celllistc,celllistb,&
     celllistr,ipcc,ipcb,ipcr,kfactor,ligboundto,recboundto,rcut_blobblob,rcut_blobwall,rcut_ligrec,&
     rxcolgeo,random_anconcolpos,nrxspec,irecspec,ichainspec,RXBONDENE,tot_ene)
  implicit none
  integer, intent(in) :: mnpic,ncellsc(3),ncellsb(3),ncellsr(2),maxncol,nchainspercol,nrecs(nrxspec),&
       nblobsperchain,kfactor,maxnblob,nrec,nrxspec,irecspec(nrec),ichainspec(nchainspercol)
  integer, intent(inout) :: ncol,ligboundto(nchainspercol,maxncol),recboundto(2,nrec),&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),&
       celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob),&
       celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec)
  real*8, intent(in) :: lbox(3),rcol,socc(3),socb(3),socr(2),rxcolgeo(3,nchainspercol),&
       rcut_blobblob,rcut_blobwall,rcut_ligrec,RXBONDENE(nrxspec,nrxspec)
  real*8, intent(inout) :: poscol(3,maxncol),posrec(2,nrec),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),tot_ene
  logical, intent(in) :: random_anconcolpos
  real*8 :: Eold,Enew,rnd,rnd3(3),rinew(3),riold(3),arg,volume,ctb(3,nchainspercol),&
       S1rot,S2rot,V1rot(2),V2rot(2),ROTMATRIX(3,3),sqrtSrot,q(4),rnd2(2),&
       trialpos(3,kfactor),trialdist,xm,ym,opnemlre,effdist,Vmars(2),Smars,sqrtSmars,&
       dirvector(3),expmUk(kfactor),norexpmUk(kfactor),Wtot,commnorexpmUk,distij(3),&
       Ukext(kfactor),Uext(nblobsperchain,nchainspercol),Uk(kfactor),&
       W(nblobsperchain,nchainspercol),ext_col_ene,lowligz,goodene
  integer :: kk,parti, inewcell(3),jj,ii,ifreecol,iifreecol,ibond,&
       ipcnew(3),nfreecol,iblob,iiblob,jblob,jjblob,iichain,jjchain,itrial,iicol,&
       goodiirec,goodrecs,lowlig,iirec,iirs        
  logical :: overlap
  real*8, parameter :: pi=3.14159265d0, O=0.0d0

  Uext=0
  
  if (ncol .ne. 0) then
     write(*,*) 'something wrong in insert_col_WL routine... colloid should not be here!'
  endif
  
  iicol=ncol+1
  rinew=(/lbox(1)/2,lbox(2)/2,rcol+2 /) ! position the colloid in the centre
  inewcell(:) = floor(rinew(:)/socc(:))+1
  poscol(:,iicol)=rinew(:)
  
  ! ANCHOr CHAINS TO THE PARTICLE
  if (random_anconcolpos) then ! random positions of chain anchors on colloid
     do iichain=1,nchainspercol
        ! get a random point on a sphree surface using Marsaglia method
        Smars=2
        do while (Smars .ge. 1.0)
           call random_number(rnd2)
           Vmars(:)=2*rnd2(1:2)-1.0d0
           Smars=sum(Vmars*Vmars)
        enddo
        sqrtSmars=sqrt(1-Smars)
        dirvector(:)=(/2*Vmars(1)*sqrtSmars,2*Vmars(2)*sqrtSmars,1-2*Smars/) ! random point on sphere surface
        ctb(:,iichain)=rcol*dirvector(:)
        posblob(:,1,iichain,iicol)=poscol(:,iicol)+ctb(:,iichain)
     enddo
     
  else
     ctb(:,:)=rcol*rxcolgeo(1:3,:)
     
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
     do iichain=1,nchainspercol
        ctb(:,iichain)=matmul(ROTMATRIX,ctb(:,iichain))
        posblob(:,1,iichain,iicol)=poscol(:,iicol)+ctb(:,iichain)
     enddo
     
  endif ! random_anconcolpos
  
  ext_col_ene=0 ! no colloid or blobs in the system yet WL
  
  ! DO ROSENBLUTH SAMPLING FOR ALL CHAINS
  do iichain=1,nchainspercol
     ! should randomize the rx chains when calculating Rosenbluth factor
     ! if positions of anchors are random then it doesn't matter
     do iiblob=2,nblobsperchain
        do itrial=1,kfactor
           
           if (iiblob .eq. 2) then
              ! get a distance distribution using a rejection method
              do ! infinite do loop -- bad coding practice but it's so simple that it shouldn't make problems
                 call random_number(rnd2)
                 xm=rnd2(1)*5  ! 6 should be enough as then  P(6)~10-5
                 ym=rnd2(2)*0.5 ! so that we cover the peak of the distribution
                 if (xm*xm*exp(-0.75*xm**2) .ge. ym) exit
              enddo
           else
              ! get a distance distribution using a rejection method
              do ! infinite do loop -- bad coding practice but it's so simple that it shouldn't make problems
                 call random_number(rnd2)
                 xm=rnd2(1)*6  ! 6 should be enough as then  P(6)~10-5
                 ym=rnd2(2)*1.8 ! so that we cover the peak of the distribution
                 if (xm*xm*exp(-0.534*(xm-0.73)**2) .ge. ym) exit
              enddo
           endif ! get distance of a new blob
           
           trialdist=xm
           ! get a point on a sphere surface using marsaglia method
           ! MARSAGLIA METHOD
           Smars=2
           do while (Smars .ge. 1.0)
              call random_number(rnd2)
              Vmars(:)=2*rnd2(1:2)-1.0d0
              Smars=sum(Vmars*Vmars)
           enddo
           sqrtSmars=sqrt(1-Smars)
           dirvector(:)=(/2*Vmars(1)*sqrtSmars,2*Vmars(2)*sqrtSmars,1-2*Smars/) ! random point on sphere surface
           trialpos(:,itrial)=posblob(:,iiblob-1,iichain,iicol)+trialdist*dirvector(:) ! get trial position
           call ext_blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
                nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
                socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
                trialpos(:,itrial),iicol,iichain,iiblob,Ukext(itrial))
           
           Uk(itrial)=Ukext(itrial) 
        enddo ! trials
        expmUk(:)=exp(-Uk(:))
        W(iiblob,iichain)=sum(expmUk(:))
        ! select one trial at random weighted by boltzmann
        norexpmUk=expmUk(:)/W(iiblob,iichain) ! normalised boltzman weight
        itrial=1
        commnorexpmUk=norexpmUk(1)  ! commulative normalised exponent of minus Uk
        call random_number(rnd)
        do while (commnorexpmUk .lt. rnd) 
           itrial=itrial+1
           commnorexpmUk=commnorexpmUk+norexpmUk(itrial)
        enddo
        if (itrial .gt. kfactor) write(*,*) 'ERROR WITH ROSENBLUTH: itrial>kfactor'
        ! end select random trial
        
        ! commit to posblob(:,:,:)
        posblob(:,iiblob,iichain,iicol)=trialpos(:,itrial)
        
        ! get new blob-blob distance
        distij(:)=posblob(:,iiblob,iichain,iicol)-posblob(:,iiblob-1,iichain,iicol)
        
        ! only commit external energy, not the rec-lig bond energy
        if (iiblob .eq. 2 ) then
           Uext(iiblob,iichain)=Ukext(itrial)+0.75*(sum(distij*distij))
        else
           Uext(iiblob,iichain)=Ukext(itrial)+0.534*(sqrt(sum(distij*distij))-0.73)**2
        endif
        ! do periodic boundary
        posblob(1:2,iiblob,iichain,iicol)=posblob(1:2,iiblob,iichain,iicol)-lbox(1:2)*&
             floor(posblob(1:2,iiblob,iichain,iicol)/lbox(1:2))
        
     enddo ! blobs in 1 chain
  enddo ! nchainspercol       
  
  Wtot=product(W(2:nblobsperchain,:)/kfactor)           
  ! END ROSENBLUTH
  
   !ALWAYS ACCEPT MOVE AND INSERT _WL
  ncol=ncol+1
  if (ncol .ne. 1) then
     write(*,*) ' ERROR: Colloid number should be 1 in WL insert, exiting... '
     call exit()
  endif
  poscol(:,ncol)=rinew(:)           
  tot_ene=tot_ene+sum(Uext(2:nblobsperchain,:))+ext_col_ene
  
  ! update celllists - col
  ipcc(1:3,ncol)=inewcell(:)
  celllistc(1,inewcell(1),inewcell(2),inewcell(3)) = &
       celllistc(1,inewcell(1),inewcell(2),inewcell(3))+1
  ipcc(4,ncol)=celllistc(1,inewcell(1),inewcell(2),inewcell(3))
  celllistc(ipcc(4,ncol)+1,inewcell(1),inewcell(2),inewcell(3))=ncol
  
  ! update cell list - blobs
  do iichain=1,nchainspercol
     do iiblob=2,nblobsperchain
        iblob=iiblob+(nblobsperchain)*(iichain-1)+nblobsperchain*nchainspercol*(ncol-1) ! transform to a serial number of a blob
        ipcnew=floor(posblob(:,iiblob,iichain,ncol)/socb(:))+1
        
        if (ipcnew(3).lt.1) then ! needed because soft blobs can slightly penetrate the wall
           ipcnew(3)=1
        elseif (ipcnew(3).gt.ncellsb(3)) then
           ipcnew(3)=ncellsb(3)
        endif
        celllistb(celllistb(1,ipcnew(1),ipcnew(2),ipcnew(3))+2, ipcnew(1),ipcnew(2),&
             ipcnew(3))=iblob
        celllistb(1,ipcnew(1),ipcnew(2),ipcnew(3))=celllistb(1,ipcnew(1),ipcnew(2),ipcnew(3))+1
        
        ipcb(1:3,iblob)=ipcnew
        ipcb(4,iblob)=celllistb(1,ipcnew(1),  ipcnew(2), ipcnew(3))   
        ! end update cell list
     enddo !iiblob
  enddo !iichain

  ! BIND WITH THE SURFACE, ATTACH A LIGAND CLOSEST TO THE SURFACE TO A RECEPTOR
  ! find a ligand closest to the surface
  lowlig=-1
  lowligz=lbox(3)
  
  do iichain=1,nchainspercol
     if (posblob(3,nblobsperchain,iichain,ncol) .lt. lowligz) then
        lowligz=posblob(3,nblobsperchain,iichain,ncol)
        lowlig=iichain
     endif
  enddo
  if (lowlig .lt. 0) then
     write(*,*) 'WARNING: Error while trying to bind the first ligand.....  exiting '
     call exit()
  endif
  if (lowligz .gt. rcut_ligrec) then
     write(*,*) 'WARNING: low ligand height .gt. rcut_ligrec when trying to bind the first ligand.... might crash later'
    ! write(*,*) ' mooving the ligand closer to the surface... '
     
  endif

  ! find a cognate receptor with strongest binding  and bind to it.
  goodene=100 
  goodrecs=-1
  do iirs=1,nrxspec
     if (nrecs(iirs) .gt. 0) then ! consider it
        if (RXBONDENE(iirs,ichainspec(lowlig)) .lt. goodene) then
        ! choese this one
           goodene=RXBONDENE(iirs,ichainspec(lowlig))
           goodrecs=iirs
        endif
     endif
  enddo
  !write(*,*) nrxspec,nrecs
  !write(*,*) lowlig
  if (goodrecs .lt. 0) then
     write(*,*) 'WARNING: could not find a good receptor candidate for the first bond... exiting'
     call exit()
  endif
  
  write(*,*) ' '
  write(*,*) 'Found a good first bond candidate receptor type: recs=',goodrecs,', bond energy=',&
       goodene
  ! find a good receptor 
  do iirec=1,nrec
     if (irecspec(iirec) .eq. goodrecs) then
        goodiirec=iirec
        exit
     endif
  enddo
  ! get ligand-receptor  distance
  distij(1:2)=posblob(1:2,nblobsperchain,lowlig,1)-posrec(1:2,goodiirec)
  distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))  ! do periodic boundary
  
  ! move all receptors uniformly such that the ii receptor is beneth the chosen low ligand
  do iirec=1,nrec
     posrec(1:2,iirec)=posrec(1:2,iirec)+distij(1:2)
     posrec(1:2,iirec)=posrec(1:2,iirec)-lbox(1:2)*floor(posrec(1:2,iirec)/lbox(1:2)) ! ...should work.
  enddo
  write(*,*) ' Found good receptor', goodiirec ,' binding to ligand ', lowlig ,&
       ' at height ', posblob(3,nblobsperchain,lowlig,1)
  write(*,*) 'Lateraly mooved all receptors to the optimal binding position,   must redo receptor cell lists!! '

  do kk=1, nrec
     ipcr(1:2,kk)=floor(posrec(:,kk)/socr(:))+1
     !write(*,*) poscol(:,kk),'   ', ipcc(1:3,kk),'  ',socc
     !write(*,*)
  enddo
  do ii=1,nrec
     celllistr(1, ipcr(1,ii), ipcr(2,ii)) = & 
          celllistr(1, ipcr(1,ii), ipcr(2,ii)) + 1
     celllistr((celllistr(1,ipcr(1,ii),ipcr(2,ii))+1),ipcr(1,ii),ipcr(2,ii))=ii
     
     ipcr(3,ii)=celllistr(1, ipcr(1,ii), ipcr(2,ii))
  enddo
  write(*,*) 'Done new receptor cell lists.'
  
  ! binding together
  ligboundto(lowlig,ncol)=goodiirec
  recboundto(1:2,goodiirec)=(/lowlig,ncol/)
  
end subroutine insert_col_WL

! ========================================================================== !
! ========================================================================== !
subroutine  WangLandau_bonds(nchainspercol,nbonds,psiWL,&
     hWL,fWL,fWL_stop,fWL_reduce,hWL_conv,icycle,WLacc,doWL_flag)
  
  
  implicit none
  integer, intent(in) :: nchainspercol,nbonds
  integer, intent(inout) :: WLacc
  integer*8, intent(in) :: icycle
  integer*8, intent(inout) :: hWL(nchainspercol+1)
  real*8, intent(in) :: fWL_stop,fWL_reduce,hWL_conv
  real*8, intent(inout) :: psiWL(nchainspercol+1),fWL
  logical, intent(inout) :: doWL_flag
  
 
  ! ---- UPDATE HISTOGRAM ------- !
 
  hWL(nbonds+1)=hWL(nbonds+1)+1 ! update histogram

  ! ---- UPDATE PSI ------ !
  if (doWL_flag) then
     psiWL(nbonds+1)=psiWL(nbonds+1)-fWL
     if ((minval(hWL) .gt. 100) .and. (dble(maxval(hWL))/minval(hWL) .lt. hWL_conv )) then ! update fWL
        fWL=fWL/fWL_reduce
        
        write(*,*)
        write(*,*) 'new WL step: icycle=',icycle,',  f=',fWL,',  hist_rat=',dble(maxval(hWL))/minval(hWL)
        hWL(:)=0
        psiWL(:)=psiWL(:)-psiWL(1)
        if (fWL .lt. fWL_stop) then ! stop Wang-Landau
           fWL=0
           doWL_flag=.false.
           
           write (*,*) 'WANG-LANDAU finished at cycle ',icycle,',   running equilibrium simulation....'
        endif
     endif
  endif
end subroutine WangLandau_bonds
!============================================================================!

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
! XXXXXXXXXXXXXXXXXX   MC   MOVE      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
subroutine mc_move_wl(lbox,ncol,maxncol,nchainspercol,nblobsperchain,nrec,maxnblob,rcol,&
     poscol,posblob,posrec,mnpic,ncellsc,ncellsb,ncellsr,socc,socb,socr,celllistc,&
     celllistb,celllistr,ipcc,ipcb,ipcr,ligboundto,recboundto,rcut_blobblob,&
     rcut_blobwall,rcut_ligrec,ligrec_ene,max_hop_blob,max_hop_col,max_rot_col,&
     moveacc,tot_ene,RXBONDENE,nrxspec,irecspec,ichainspec,mobile_ligands,psiWL,max_colzWL,nbonds)
  
  implicit none
  integer, intent(in) :: mnpic,ncellsc(3),ncol,maxncol,nchainspercol,nblobsperchain,&
       ncellsb(3),ncellsr(2),maxnblob,nrec,celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec),&
       nrxspec,irecspec(nrec),ichainspec(nchainspercol)
  integer, intent(inout) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol),celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),&
       ipcb(4,maxnblob),moveacc,ligboundto(nchainspercol,maxncol),&
       recboundto(2,nrec),nbonds
  real*8, intent(in) :: lbox(3),posrec(2,nrec),rcol,socc(3),socb(3),socr(2),rcut_blobblob,&
       rcut_blobwall,rcut_ligrec,ligrec_ene,max_hop_col,max_rot_col,max_hop_blob,&
       RXBONDENE(nrxspec,nrxspec),max_colzWL
  logical, intent(in) :: mobile_ligands
  real*8, intent(inout) :: poscol(3,maxncol),posblob(3,nblobsperchain,&
       nchainspercol,maxncol),tot_ene,psiWL(nchainspercol+1)
  real*8 :: Eold,Enew,rnd,rnd3(3),rinew(3),riold(3),arg,volume,ctb(3,nchainspercol),&
       Srot,Vrot(2),ROTMATRIX(3,3),sqrtSrot,urot(3),sinfi,cosfi,omcosfi,rnd2(2),&
       distij(3),Qold,Qnew,rancnew(3,nchainspercol),garg,q1bi_old(mnpic),q1bi_new(mnpic),&
       rcut_ligrec2
  integer :: kk,parti, inewcell(3),jj,ii,ifreecol,iifreecol,ibond,ipcold(3),&
       ipcnew(3),nfreecol,iblob,iiblob,jblob,jjblob,iichain,jjchain,itrial,iicol,&
       jrecind,jjrec,nrechomies_new,nrechomies_old,whichhomies_old(mnpic),&
       whichhomies_new(mnpic),ihomie
  logical :: overlap,lastbond=.false.
  real*8, parameter :: pi=3.14159265d0, O=0.0d0 

  if (ncol .ne. 1) return
  rcut_ligrec2=rcut_ligrec*rcut_ligrec

  garg=1.d0/(nblobsperchain*nchainspercol)
  call random_number(rnd)  
  if (rnd .lt. garg) then ! move colloid
    
     iicol=1     
     call random_number(rnd3)  ! get new position  
     rinew=poscol(:,iicol)+max_hop_col*2*(rnd3(:)-0.5d0)
     ! cannot go out of height max_colzWL for Wang-Landau and unbound colloid
     if ((rinew(3) .gt. max_colzWL ) .and. (rinew(3) .gt. poscol(3,iicol)) .and. &
          (sum(ligboundto(:,1)) .eq.0))  return 
     rinew(1:2)=rinew(1:2)-lbox(1:2)*floor(rinew(1:2)/lbox(1:2))
     inewcell(:) = floor(rinew(:)/socc(:))+1
     !CHECK FOR OVERLAP WITH OTHER COLLOIDS
     call col_overlap(overlap,iicol,rinew,inewcell,mnpic,&
          celllistc,ipcc,ncellsc,lbox,poscol,rcol,maxncol)
     ! overlap with wall
     if ((rinew(3).lt. rcol).or.(rinew(3).gt.lbox(3)-rcol)) overlap =.true.     
     
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
        
        ! GET RELATIVE POSITION VECTORS, BIND SITE TO PARTICLE CENTRE     
        do ii=1,nchainspercol
           ctb(:,ii)=posblob(:,1,ii,iicol)-poscol(:,iicol) ! centre of colloid to anchor pos vector
           ctb(1:2,ii)=ctb(1:2,ii)-lbox(1:2)*nint(ctb(1:2,ii)/lbox(1:2))
           ctb(:,ii)=matmul(ROTMATRIX(:,:),ctb(:,ii))
           rancnew(:,ii)=rinew(:)+ctb(:,ii) ! NEW POSITION OF ANCHORS
        enddo
    
        ! CALCULATE OLD & NEW ENERGY
        Eold=0
        Enew=0       
        call col_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
             poscol(:,iicol),iicol,Eold)
        call col_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
             rinew(:),iicol,Enew)

        ! ADD HARMONIC ENERGY
        do iichain=1,nchainspercol
           distij(:)=posblob(:,2,iichain,iicol)-posblob(:,1,iichain,iicol)     
           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))       
           Eold=Eold+0.75*sum(distij*distij) ! energy of the first blob-anchor

           distij(:)=posblob(:,2,iichain,iicol)-rancnew(:,iichain)     
           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))       
           Enew=Enew+0.75*sum(distij*distij)! energy of the first blob-anchor
        enddo ! iichain
        
        arg=exp(Eold-Enew)
        call random_number(rnd)
        if (rnd .lt. arg) then ! ACCEPT MOVE
           moveacc=moveacc+1
           posblob(:,1,:,iicol)=rancnew(:,:)
           poscol(:,iicol)=rinew
           tot_ene=tot_ene+Enew-Eold
           if (any(inewcell .ne. ipcc(1:3,iicol))) then

              call update_cell_list(iicol,maxncol,inewcell,ipcc,celllistc,&
                   ncellsc,mnpic)
           endif
        else ! REJECT   
        endif
     endif ! OVERLAP
    
!xxxxxxxxxxxxxxxxxxxxxxxxxx  MOVE BLOBS  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx! 
  else   ! MOVE BLOBS
     call random_number(rnd3)
     iicol=floor(rnd3(1)*ncol)+1
     iichain=floor(rnd3(2)*nchainspercol)+1
     iiblob=floor(rnd3(3)*(nblobsperchain))+1   ! can also move anchors

     if (mobile_ligands .and. (iiblob .eq. 1)) then ! move ligand anchor
        call random_number(rnd3)
        rinew(:)=posblob(:,iiblob,iichain,iicol)+max_hop_blob*2*(rnd3(:)-0.5d0)
        ctb(:,1)=rinew(:)-poscol(:,iicol) ! colloid to blob distance
        ctb(:,1)=ctb(:,1)/sqrt(sum(ctb(:,1)*ctb(:,1)))*rcol ! normalize to colloid surface
        rinew(:)=poscol(:,iicol)+ctb(:,1) ! new position on colloid surface
        
        distij(:)=posblob(:,2,iichain,iicol)-posblob(:,1,iichain,iicol)
        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
        Eold=0.75*sum(distij*distij) 
        distij(:)=posblob(:,2,iichain,iicol)-rinew(:)
        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
        Enew=0.75*sum(distij*distij)

   !     write(*,*)
   !     write(*,*) posblob(:,iiblob,iichain,iicol)-rinew(:) 
   !     write(*,*) rcol, rinew
   !     write(*,*) Enew, Eold
        ! DO METROPOLIS
        arg=exp(Eold-Enew)
        call random_number(rnd)
        if (rnd .lt. arg) then! ACCEPT
           moveacc=moveacc+1
           posblob(:,1,iichain,iicol)=rinew(:)
           tot_ene=tot_ene+Enew-Eold 
        endif

     elseif (iiblob .gt. 1) then ! move blobs   
        ! get new blob position
        call random_number(rnd3)
        rinew(:)=posblob(:,iiblob,iichain,iicol)+max_hop_blob*2*(rnd3(:)-0.5d0)
        rinew(1:2)=rinew(1:2)-lbox(1:2)*floor(rinew(1:2)/lbox(1:2))
        
        Eold=0
        Enew=0
        call blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
             posblob(:,iiblob,iichain,iicol),iicol,iichain,iiblob,Eold)
        call blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
             rinew,iicol,iichain,iiblob,Enew)
        
        ! DO BINDING WITH RECEPTORS
        lastbond=.false. ! flag we are moving the last bound ligand
        Qold=1
        Qnew=1
        nrechomies_old=0
        nrechomies_new=0
        if ((iiblob .eq. nblobsperchain).and.(posblob(3,iiblob,iichain,iicol) .lt. rcut_ligrec)) then
           ! check if colloid is outside of h_0
           if ( poscol(3,iicol) .gt.  max_colzWL) then
              ! check if  this blob has the last bond
              if((nbonds .eq. 1).and.(ligboundto(iichain,iicol).gt.0 )) then ! this ligand has the last bond
                 lastbond=.true.
                 Qold=0 ! unbound state is not possible
                 Qnew=0
              elseif(nbonds .eq. 0) then
                 write(*,*) 'WARNING: nbonds = 0, should not be the case for Wl bonds.. '
                 write(*,*) 'only ok if it heapens at the begining and the initial colloid height is > h_0'
                 return
              endif
             
           endif
           if (posblob(3,nblobsperchain,iichain,iicol) .lt. rcut_ligrec) then
              
              ! can do binding
              ! first unbind
              if (ligboundto(iichain,iicol).gt.0) then
                 distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,ligboundto(iichain,iicol))
                 distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                 distij(3)=posblob(3,nblobsperchain,iichain,iicol)
                 if (sum(distij*distij) .gt. rcut_ligrec2 ) then
                    write(*,*) 'WARNING: trying to delete a bond but distij > rcut_ligrec, exiting...'
                    ! cannot unbind if distij2 > rcut2, no reverse move in configurational bias !!!!
                    !call exit()
                 endif

                 recboundto(1:2,ligboundto(iichain,iicol))=0
                 ligboundto(iichain,iicol)=0
                 nbonds=nbonds-1  ! deleted a bond
              endif
              
              call find_rec_homie(posblob(:,nblobsperchain,iichain,iicol),nchainspercol,nblobsperchain,lbox,posrec,&
                   nrec,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,recboundto,&
                   nrechomies_old,whichhomies_old)
              
              ! get the old binding partition function
                           
              if (nrechomies_old .gt. 0) then
                 q1bi_old(:)=0
                 ! get partition funcition with psi bias
                 do ihomie=1,nrechomies_old
                    distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,whichhomies_old(ihomie))
                    distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                    distij(3)=posblob(3,nblobsperchain,iichain,iicol)
                    q1bi_old(ihomie)=exp(-RXBONDENE(irecspec(whichhomies_old(ihomie)),&
                         ichainspec(iichain))-0.75*sum(distij*distij))
                 enddo
                 q1bi_old(:)=q1bi_old(:)*exp(psiWL(nbonds+2)-psiWL(nbonds+1))
                 Qold=Qold+sum(q1bi_old(1:nrechomies_old))
              endif ! nrechomies > 0
          
              ! endif  posblob < rcut_ligrec
           elseif (ligboundto(iichain,iicol) .gt. 0) then
              write(*,*)' WARNING: ligand bound to receptor but out of range: height > rcut_ligrec '
            !  write(*,*) 'Unbinding... , only ok after initialisation if low ligand height .gt. rcut_ligrec. ' 
            !  recboundto(1:2,ligboundto(iichain,iicol))=0
            !  ligboundto(iichain,iicol)=0
            !  nbonds=nbonds-1  ! deleted a bond
           endif !posblob < rcut_ligrec
           
           if (rinew(3) .lt. rcut_ligrec) then
              ! can do binding
              
              call find_rec_homie(rinew,nchainspercol,nblobsperchain,lbox,posrec,&
                   nrec,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,recboundto,&
                   nrechomies_new,whichhomies_new)
              
              ! get the new binding partition function
              
              if (nrechomies_new .gt. 0) then
                 q1bi_new(:)=0
                 ! get partition funcition
                 do ihomie=1,nrechomies_new
                    distij(1:2)=rinew(1:2)-posrec(:,whichhomies_new(ihomie))
                    distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                    distij(3)=rinew(3)
                    q1bi_new(ihomie)=exp(-RXBONDENE(irecspec(whichhomies_new(ihomie)),&
                         ichainspec(iichain))-0.75*sum(distij*distij))
                 enddo
                 q1bi_new(:)=q1bi_new(:)*exp(psiWL(nbonds+2)-psiWL(nbonds+1))
                 Qnew=Qnew+sum(q1bi_new(1:nrechomies_new))
              endif ! nrechomies > 0

           endif
         !  write(*,*)
         !  write(*,*) nbonds, lastbond
         !  write(*,*) Qold,Qnew, nrechomies_old,nrechomies_new
         !  write(*,*) rinew(3)
           
        endif ! iblob .eq. nblobsperchain
        ! END BiNDING WITH RECEPTORS
   
        ! DO METROPOLIS
        arg=exp(Eold-Enew)*Qnew/Qold
        call random_number(rnd)
        if (rnd .lt. arg) then! ACCEPT
           
           moveacc=moveacc+1
           posblob(:,iiblob,iichain,iicol)=rinew(:)
           tot_ene=tot_ene+Enew-Eold        
           
           inewcell(:) = floor(rinew(:)/socb(:))+1
           inewcell(3)=min(inewcell(3),ncellsb(3)) ! because soft blobs can in principle penetrate the wall slightly
           inewcell(3)=max(inewcell(3),1)
           iblob=(iicol-1)*nblobsperchain*nchainspercol+(iichain-1)*nblobsperchain+iiblob
           if (any(inewcell .ne. ipcb(1:3,iblob))) then
              call update_cell_list(iblob,maxnblob,inewcell,ipcb,celllistb,&
                   ncellsb,mnpic)
           endif
           
           ! bind to some receptor
           if (lastbond) then
              arg=0
           else
              arg=1.0d0/Qnew
           endif
           q1bi_new(:)=q1bi_new(:)/Qnew ! normalise partition functions
           call random_number(rnd)

           ! bind to some receptor
           if (rnd .gt. arg) then !  bind
              nbonds=nbonds+1  ! add bond to global bonds count
              do ihomie=1,nrechomies_new
                 arg=arg+q1bi_new(ihomie)
                 if (rnd .le. arg) then
                    jjrec=whichhomies_new(ihomie)
                    ligboundto(iichain,iicol)=jjrec
                    recboundto(1:2,jjrec)=(/iichain,iicol/)
                    exit
                 endif
              enddo
           endif
           
        else ! REJECT
           ! rebind to some receptor
           if (lastbond) then
              arg=0
           else
              arg=1.0d0/Qold
           endif
           q1bi_old(:)=q1bi_old(:)/Qold ! normalise partition functions
           call random_number(rnd)
           ! bind to some receptor
           if (rnd .gt. arg) then !  bind
              nbonds=nbonds+1  ! add bond to global bonds count
              do ihomie=1,nrechomies_old
                 arg=arg+q1bi_old(ihomie)
                 if (rnd .le. arg) then
                    jjrec=whichhomies_old(ihomie)
                    ligboundto(iichain,iicol)=jjrec
                    recboundto(1:2,jjrec)=(/iichain,iicol/)
                    exit
                 endif
              enddo
           endif ! bind
        endif ! accept/reject
     endif ! move anchor/other blob
  endif ! move col/blob
end subroutine mc_move_wl
!=============================================================!
subroutine hop_rec(posrec,poscol,lbox,nrec,ncol,rrec,socr,mnpic,recboundto,&
     ligboundto,posblob,nchainspercol,maxncol,nblobsperchain,maxnblob,celllistr,&
     ipcr,ncellsr,rcut_ligrec,rcut_blobblob,rcut_blobwall,max_hop_blob,tot_ene)   
  implicit none
  integer, intent(in) :: nrec,nchainspercol,ncellsr(2),mnpic,maxncol,maxnblob,&
       ncol,nblobsperchain,recboundto(2,nrec),ligboundto(nchainspercol,maxncol)
  integer, intent(inout) :: ipcr(3,nrec),celllistr(mnpic,ncellsr(1),ncellsr(2))
  real*8, intent(in) :: lbox(3),socr(2),rrec,rcut_blobwall,rcut_blobblob,&
       rcut_ligrec,max_hop_blob,poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol)
  real*8, intent(inout) :: posrec(2,nrec),tot_ene
  real*8 :: rnd,rnd2(2),rnd3(3),arg,oldpos(3),newpos(3),roldpos(2),rnewpos(2),&
       distij(3),distij2,hopmove(2),new_ene,old_ene,rcut_ligrec2,rrec_cut2
  integer :: ii,jj,kk,irec,jrec,ipcnew(2),ipcold(2),icell,jcell,icelltmp,jcelltmp,&
       iiblob,iichain,iblob,krec,iicol
  logical :: overlap=.false., even
  
  rcut_ligrec2=rcut_ligrec*rcut_ligrec
  rrec_cut2=4*rrec*rrec
  
  if (nrec .lt. 1) return 
  !pick a receptor 
  call random_number(rnd)
  irec=floor(nrec*rnd)+1
  if (recboundto(1,irec) .eq. 0) then ! move receptor to a random pos on the surface 
     
     call random_number(rnd2) ! get new random postion
     rnewpos=rnd2(:)*lbox(1:2)
     ipcnew(1:2)=floor(rnewpos(:)/socr(:))+1
     ipcold(1:2)=ipcr(1:2,irec)
     
     overlap=.false.
     ! check for overlap
     do jcelltmp=ipcnew(2)-1,ipcnew(2)+1
        do icelltmp=ipcnew(1)-1,ipcnew(1)+1     
           icell=icelltmp-ncellsr(1)*floor(dble(icelltmp-1)/ncellsr(1)+1.0d-8)
           jcell=jcelltmp-ncellsr(2)*floor(dble(jcelltmp-1)/ncellsr(2)+1.0d-8)
           do jj=1,celllistr(1,icell,jcell)
              jrec=celllistr(jj+1,icell,jcell)
              if (all((/icell,jcell/).eq.ipcnew(1:2)).and. & 
                   (irec .eq. jrec)) cycle
              distij(1:2)=rnewpos(:)-posrec(:,jrec)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij(1:2)*distij(1:2))
              if (distij2 .lt. rrec_cut2) then
                 overlap=.true.
                 return ! do not move if overlap
              endif
           enddo
        enddo
     enddo
     ! write(*,*) 'AAAAA'
     if (.not. overlap) then ! update position
        distij(1:2)=rnewpos(1:2)-posrec(1:2,irec)
        posrec(1:2,irec)=rnewpos(1:2)
        !write(*,*) newpos
        if (any(ipcold .ne. ipcnew)) then
           ! update cell lists
           jj=ipcr(3,irec) ! position in the cell
           celllistr(jj+1, ipcold(1), ipcold(2)) = &
                celllistr(celllistr(1,ipcold(1),ipcold(2))+1, ipcold(1),ipcold(2))
           ipcr(3,celllistr(jj+1,ipcold(1),ipcold(2))) = jj   
           
           celllistr(1,ipcold(1),ipcold(2)) = celllistr(1,ipcold(1),ipcold(2))-1 
           celllistr(celllistr(1,ipcnew(1),ipcnew(2))+2, ipcnew(1),ipcnew(2)) = irec
           celllistr(1,ipcnew(1),ipcnew(2))=celllistr(1,ipcnew(1),ipcnew(2))+1
           
           ipcr(1:2,irec)=ipcnew
           ipcr(3, irec)=celllistr(1,ipcnew(1),  ipcnew(2))  
           celllistr(celllistr(1,ipcold(1),ipcold(2))+2,ipcold(1),ipcold(2))=0 
        endif
       ! write(*,*) posrec(:,irec),posrec(:,jrec)
     endif ! not overlap
  else  !try to move a bonded receptor 
        ! return
     iicol=recboundto(2,irec)
     iichain=recboundto(1,irec)
    
     call random_number(rnd2)
     hopmove(:)=(2.0d0*rnd2(:)-1.0d0)*max_hop_blob
     rnewpos(1:2)=posrec(1:2,irec)+hopmove(:)
     rnewpos(1:2)=rnewpos(1:2)-lbox(1:2)*floor(rnewpos(1:2)/lbox(1:2)) ! do periodic boundary
     ipcnew(1:2)=floor(rnewpos(:)/socr(:))+1
     ipcold(1:2)=ipcr(1:2,irec)

     ! check for overlap
     overlap=.false.
     do jcelltmp=ipcnew(2)-1,ipcnew(2)+1
        do icelltmp=ipcnew(1)-1,ipcnew(1)+1     
           icell=icelltmp-ncellsr(1)*floor(dble(icelltmp-1)/ncellsr(1)+1.0d-8)
           jcell=jcelltmp-ncellsr(2)*floor(dble(jcelltmp-1)/ncellsr(2)+1.0d-8)
           do jj=1,celllistr(1,icell,jcell)
              jrec=celllistr(jj+1,icell,jcell)
              if (all((/icell,jcell/).eq.ipcnew(1:2)).and. & 
                   (irec .eq. jrec)) cycle
              distij(1:2)=rnewpos(:)-posrec(:,jrec)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij(1:2)*distij(1:2))
              if (distij2 .lt. rrec_cut2) then
                 overlap=.true.
                 return ! do not move if overlap
              endif
           enddo
        enddo
     enddo

     ! old energy    
     distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,irec)
     distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
     distij(3)=posblob(3,nblobsperchain,iichain,iicol)
     old_ene=0.75*sum(distij*distij)  
     ! new energy
     distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-rnewpos(:)
     distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
     distij2=sum(distij*distij)
     if (distij2 .gt. rcut_ligrec2) return ! return if new receptor position is outside of the bond rcut_ligrec
     new_ene=0.75*distij2 


     arg=exp(old_ene-new_ene)
     call random_number(rnd)
     if (arg .gt. rnd) then ! ACCEPT
       ! tot_ene=tot_ene+new_ene-old_ene 
    !    write(*,*) 'XXXXXX',new_ene,old_ene
        posrec(:,irec)=rnewpos(:)
        ! update cell list
        if (any(ipcold .ne. ipcnew)) then
           ! update cell lists
           jj=ipcr(3,irec) ! position in the cell
           celllistr(jj+1, ipcold(1), ipcold(2)) = &
                celllistr(celllistr(1,ipcold(1),ipcold(2))+1, ipcold(1),ipcold(2))
           ipcr(3,celllistr(jj+1,ipcold(1),ipcold(2))) = jj   
           
           celllistr(1,ipcold(1),ipcold(2)) = celllistr(1,ipcold(1),ipcold(2))-1 
           celllistr(celllistr(1,ipcnew(1),ipcnew(2))+2, ipcnew(1),ipcnew(2)) = irec
           celllistr(1,ipcnew(1),ipcnew(2))=celllistr(1,ipcnew(1),ipcnew(2))+1
           
           ipcr(1:2,irec)=ipcnew
           ipcr(3, irec)=celllistr(1,ipcnew(1),  ipcnew(2))  
           celllistr(celllistr(1,ipcold(1),ipcold(2))+2,ipcold(1),ipcold(2))=0 
        endif
     endif
  endif
end subroutine hop_rec
   
!===========================================================================!
subroutine col_overlap(overlap,icol,rinew,icell,mnpic,celllistc,ipcc,&
     ncellsc,lbox,poscol,rcol,maxncol)
  implicit none
  
  integer, intent(in) :: icol,mnpic,icell(3),ncellsc(3),&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),maxncol,&
       ipcc(4,maxncol)
  real*8, intent(in) :: rinew(3), poscol(3,maxncol),rcol,lbox(3)
  logical, intent(out) :: overlap
  real*8 ::  distij(3), distij2, distij1,sigma2
  integer :: jj, kk, xc, yc, zc, cxc, cyc, czc,jcol,jspecie
 
  sigma2=4*rcol*rcol
  overlap=.false.
  do czc=icell(3)-1, icell(3)+1
     if ((czc .lt. 1).or.(czc .gt.ncellsc(3))) cycle
     do yc=icell(2)-1, icell(2)+1
        cyc=yc-ncellsc(2)*floor(real(yc-1)/ncellsc(2))
        do xc=icell(1)-1, icell(1)+1
           cxc=xc-ncellsc(1)*floor(real(xc-1)/ncellsc(1))
          
           do jj=1, celllistc(1,cxc, cyc, czc)  
              if (celllistc(jj+1,cxc,cyc,czc) .ne. icol) then
                 
                 jcol=celllistc(jj+1,cxc,cyc,czc)                                 
                 distij=rinew(:)-poscol(:,jcol)
                 distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
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
subroutine bond_create_destroy(posrec,posblob,lbox,nrec,ncol,maxncol,nchainspercol,& 
     nblobsperchain,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,bindacc,&
     ligboundto,recboundto,nrxspec,RXBONDENE,irecspec,ichainspec)
  implicit none
  integer, intent(in) :: maxncol,ncol,nrec,nblobsperchain,nchainspercol,ncellsr(2),mnpic,&
       celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec),irecspec(nrec),&
       ichainspec(nchainspercol),nrxspec
  integer, intent(inout) :: bindacc,ligboundto(nchainspercol,maxncol),recboundto(2,nrec)
  real*8, intent(in) :: RXBONDENE(nrxspec,nrxspec),posrec(2,nrec),lbox(3),socr(2),rcut_ligrec,posblob(3,&
       nblobsperchain,nchainspercol,maxncol)
  real*8 :: rnd,rnd2(2),rnd3(3),arg,q1b,distij(3),distij2,rcut_ligrec2,q1bi(mnpic)
  integer :: kk,itlr,jj,ii,iicol,iichain,jjrec,nrechomies,whichhomies(mnpic),ihomie,jrecind
  real*8, parameter :: pi=3.14159265d0 
  
  rcut_ligrec2=rcut_ligrec*rcut_ligrec
  if (ncol .lt. 1) return

  ! pick a random ligand
  call random_number(rnd2)
  iicol=floor(ncol*rnd2(1))+1 
  iichain=floor(nchainspercol*rnd2(2))+1 

  if (posblob(3,nblobsperchain,iichain,iicol) .lt. rcut_ligrec) then
  !   write(*,*)
  !   write(*,*) 'binding iicol iichain ',iicol,iichain , ncol
  !   write(*,*) posblob(3,nblobsperchain,iichain,iicol)
     ! can do binding
  
     ! first unbind
     if (ligboundto(iichain,iicol).gt.0) then
        distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,ligboundto(iichain,iicol))
        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
        distij(3)=posblob(3,nblobsperchain,iichain,iicol)
        if (sum(distij*distij) .gt. rcut_ligrec2 ) then
           write(*,*) 'WARNING: bond above rcut_ligrec, in bond_create_destr'
           return ! cannot unbind if distij2 > rcut2, no reverse move in configurational bias !!!!
        endif
        recboundto(1:2,ligboundto(iichain,iicol))=0
        ligboundto(iichain,iicol)=0
     endif
     
     call find_rec_homie(posblob(:,nblobsperchain,iichain,iicol),nchainspercol,nblobsperchain,lbox,posrec,&
          nrec,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,recboundto,&
          nrechomies,whichhomies)
     
     
     if (nrechomies .gt. 0) then
        bindacc=bindacc+1
        q1bi(:)=0
        ! get partition funcition
        do ihomie=1,nrechomies
           distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,whichhomies(ihomie))
           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
           distij(3)=posblob(3,nblobsperchain,iichain,iicol)
           q1bi(ihomie)=exp(-RXBONDENE(irecspec(whichhomies(ihomie)),&
                         ichainspec(iichain))-0.75*sum(distij*distij))
        enddo
       
        q1b=1+sum(q1bi(1:nrechomies))
        arg=1.0d0/q1b
        q1bi(:)=q1bi(:)/q1b ! normalise partition functions
        call random_number(rnd)
        ! bind to some receptor
        if (rnd .gt. arg) then !  bind
           do ihomie=1,nrechomies
              arg=arg+q1bi(ihomie)
              if (rnd .le. arg) then
                 jjrec=whichhomies(ihomie)
                 ligboundto(iichain,iicol)=jjrec
                 recboundto(1:2,jjrec)=(/iichain,iicol/)
                 exit
              endif
           enddo
        endif
     endif ! nrechomies > 0
  endif  
end subroutine bond_create_destroy
!==========================================================================!
subroutine ext_blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
              nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
              trialpos,iicol,iichain,iiblob,pot_ene)
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),iichain,iiblob,&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),iicol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),trialpos(3),&
       rcut_blobblob,rcut_blobwall,rcol,poscol(3,maxncol)
  real*8, intent(out) ::pot_ene
  real*8 :: Eold,Enew,rnd,rnd2(2),rnd3(3),arg,volume,effdist,delene,&
       rcut_blobcol,rcut_blobcol2,rbc2
  integer :: jjcol,icelltmp,jcelltmp,kcelltmp,icell,jcell,kcell,& 
       iblob,jblob,jjblob,jjchain,cellrat(3),tmpcell(3),ii,jj,ipcnew(3),&
       ipcbc(3)
  real*8 :: distij(3),distij2,rbb,rbb2,rbw,&
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob
  real*8, parameter :: pi=3.14159265d0
  
  rcut_blobblob2=rcut_blobblob**2
  rcut_blobcol=rcut_blobwall+rcol
  rcut_blobcol2=rcut_blobcol**2
  eCS_blobwall=3.2*exp(-4.17*(rcut_blobwall-0.5d0))
  eCS_blobblob=1.75*exp(-0.8*rcut_blobblob2)
  
  pot_ene=0  

  ! blob blob energy
  ipcnew=floor(trialpos(:)/socb(:))+1
  if (ipcnew(3) .lt. 1) ipcnew(3)=1
  if (ipcnew(3) .gt. ncellsb(3)) ipcnew(3)=ncellsb(3)
  iblob=(iicol-1)*nchainspercol*nblobsperchain+(iichain-1)*nblobsperchain+iiblob
  do kcell=ipcnew(3)-1,ipcnew(3)+1
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsb(3))) cycle 
     do jcelltmp=ipcnew(2)-1,ipcnew(2)+1
        do icelltmp=ipcnew(1)-1,ipcnew(1)+1
           
           icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
           jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
           
           do jj=1,celllistb(1,icell,jcell,kcell)
              jblob=celllistb(jj+1,icell,jcell,kcell)
              if (all((/icell,jcell,kcell/).eq.ipcnew(1:3)).and. & 
                   (iblob .eq. jblob)) cycle
             ! jjchain=ceiling(jblob/dble(nblobsperchain)-1.0d-8)
              jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)              
              if (iicol .eq.jjcol) cycle ! only do other chains because the trial chain is not in the cell list when inserting particles
              jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
              jjblob=mod(jblob,nblobsperchain)
              if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
              if (jjblob .eq. 1) write(*,*) 'WARNING: anchor should not be in the cell list '
         !     write(*,*) jjcol,jjchain,jjblob
              distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbb2=distij2
              if (rbb2 .lt. rcut_blobblob2) then
                 pot_ene=pot_ene + 1.75*exp(-0.8*rbb2)-eCS_blobblob
              endif
           enddo
        enddo
     enddo
  enddo


  
  ! add intra chains on given colloid interactions
  ! only with blobs that are already present.. 
  ! this works assumming that chains are inserted in order iichain=1,2,3,4. which is only ok for random position of anchors on col
  do jjchain=1,iichain-1
     do jjblob=2,nblobsperchain
        distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,iicol)
        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
        distij2=sum(distij*distij)
        rbb2=distij2
        if (rbb2 .lt. rcut_blobblob2) then
           pot_ene=pot_ene + 1.75*exp(-0.8*rbb2)-eCS_blobblob
        endif
     enddo
  enddo
  jjchain=iichain
  do jjblob=2,iiblob-1
     distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,iicol)
     distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
     distij2=sum(distij*distij)
     rbb2=distij2
     if (rbb2 .lt. rcut_blobblob2) then
        pot_ene=pot_ene + 1.75*exp(-0.8*rbb2)-eCS_blobblob
     endif
  enddo
 
  
  ! blob col energy
  ! with self colloid 
  distij(:)=trialpos(:)-poscol(:,iicol)
  distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
  distij2=sum(distij*distij)
  rbc2=distij2
  if (rbc2 .lt. rcut_blobcol2) then
    pot_ene=pot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall
!     write(*,*) ' ext_ene: ', trialpos, poscol(:,iicol)
!     write(*,*) 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0))     
  endif
  
  ! with other colloids
  ipcbc(:)=floor(trialpos(:)/socc(:))+1  ! tells in which col cell the blob is
  if (ipcbc(3).lt.1) ipcbc(3)=1  ! because blobs can penetrate the wall
  if (ipcbc(3).gt.ncellsc(3)) ipcbc(3)=ncellsc(3) ! because blobs can penetrate the wall
  
  do kcell=ipcbc(3)-1,ipcbc(3)+1
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsc(3))) cycle 
     do jcelltmp=ipcbc(2)-1,ipcbc(2)+1
        do icelltmp=ipcbc(1)-1,ipcbc(1)+1
           
           icell=icelltmp-ncellsc(1)*floor(dble(icelltmp-1)/ncellsc(1)+1.0d-8)
           jcell=jcelltmp-ncellsc(2)*floor(dble(jcelltmp-1)/ncellsc(2)+1.0d-8)
           
           do jj=1,celllistc(1,icell,jcell,kcell)
              jjcol=celllistc(jj+1,icell,jcell,kcell)                         
              if (iicol .eq. jjcol) cycle ! only do other colloids  because the trial col is not yet in the cell list          
              distij(:)=trialpos(:)-poscol(:,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbc2=distij2
              if (rbc2 .lt. rcut_blobcol2) then
                 pot_ene=pot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall

              endif              
           enddo
        enddo
     enddo
  enddo
  
  ! blob-wall energy
  if (ipcnew(3) .eq. 1) then
     rbw=trialpos(3)
     if (rbw .lt. rcut_blobwall) then
        pot_ene=pot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
     endif
  elseif (ipcnew(3) .eq. ncellsb(3)) then
     rbw=lbox(3)-trialpos(3)
     if (rbw .lt. rcut_blobwall) then
        pot_ene=pot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
     endif
  endif

end subroutine ext_blob_energy_calc
!==========================================================================!
subroutine blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
              nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
              trialpos,iicol,iichain,iiblob,pot_ene)
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),iichain,iiblob,&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),iicol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),trialpos(3),&
       rcut_blobblob,rcut_blobwall,rcol,poscol(3,maxncol)
  real*8, intent(out) ::pot_ene
  real*8 :: Eold,Enew,rnd,rnd2(2),rnd3(3),arg,volume,effdist,delene,&
       rcut_blobcol,rcut_blobcol2,rbc2
  integer :: jjcol,icelltmp,jcelltmp,kcelltmp,icell,jcell,kcell,& 
       iblob,jblob,jjblob,jjchain,cellrat(3),tmpcell(3),ii,jj,ipcnew(3),&
       ipcbc(3)
  real*8 :: distij(3),distij2,rbb,rbb2,rbw,&
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob
  real*8, parameter :: pi=3.14159265d0
  
  rcut_blobblob2=rcut_blobblob**2
  rcut_blobcol=rcut_blobwall+rcol
  rcut_blobcol2=rcut_blobcol**2
  eCS_blobwall=3.2*exp(-4.17*(rcut_blobwall-0.5d0))
  eCS_blobblob=1.75*exp(-0.8*rcut_blobblob2)
  
  pot_ene=0  

  ! blob blob energy
  ipcnew=floor(trialpos(:)/socb(:))+1
  if (ipcnew(3) .lt. 1) ipcnew(3)=1
  if (ipcnew(3) .gt. ncellsb(3)) ipcnew(3)=ncellsb(3)
  iblob=(iicol-1)*nchainspercol*nblobsperchain+(iichain-1)*nblobsperchain+iiblob
  do kcell=ipcnew(3)-1,ipcnew(3)+1
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsb(3))) cycle 
     do jcelltmp=ipcnew(2)-1,ipcnew(2)+1
        do icelltmp=ipcnew(1)-1,ipcnew(1)+1
           
           icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
           jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
           
           do jj=1,celllistb(1,icell,jcell,kcell)
              jblob=celllistb(jj+1,icell,jcell,kcell)
              if (iblob .eq. jblob) cycle
             ! jjchain=ceiling(jblob/dble(nblobsperchain)-1.0d-8)
              jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)              
     !         if (iicol .eq.jjcol) cycle ! only do other chains because the trial chain is not in the cell list when inserting particles
              jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
              jjblob=mod(jblob,nblobsperchain)
              if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
              if (jjblob .eq. 1) write(*,*) 'WARNING: anchor should not be in the cell list '
         !    write(*,*) jjcol,jjchain,jjblob
              distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbb2=distij2
              if (rbb2 .lt. rcut_blobblob2) then
                 pot_ene=pot_ene + 1.75*exp(-0.8*rbb2)-eCS_blobblob
              endif
           enddo
        enddo
     enddo
  enddo
 
  ! harmonic energy
  distij(:)=trialpos(:)-posblob(:,iiblob-1,iichain,iicol)     
  distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))       
  if ((iiblob-1) .eq. 1) then! bound to anchor
     pot_ene=pot_ene+0.75*(sum(distij*distij)) ! energy of the first blob to anchor
  else
     pot_ene=pot_ene+0.534*(sqrt(sum(distij*distij))-0.73)**2 ! energy of the first blob
  endif
  if (iiblob .ne. nblobsperchain) then ! IF NOt THE LAST BLOB IN A CHAIN
     distij(:)=trialpos(:)-posblob(:,iiblob+1,iichain,iicol)     
     distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))       
     pot_ene=pot_ene+0.534*(sqrt(sum(distij*distij))-0.73)**2 ! energy of the first blob
  endif
  
  
  ! blob col energy
  ipcbc(:)=floor(trialpos(:)/socc(:))+1  ! tells in which col cell the blob is
  if (ipcbc(3).lt.1) ipcbc(3)=1  ! because blobs can penetrate the wall
  if (ipcbc(3).gt.ncellsc(3)) ipcbc(3)=ncellsc(3) ! because blobs can penetrate the wall
  
  do kcell=ipcbc(3)-1,ipcbc(3)+1
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsc(3))) cycle 
     do jcelltmp=ipcbc(2)-1,ipcbc(2)+1
        do icelltmp=ipcbc(1)-1,ipcbc(1)+1
           
           icell=icelltmp-ncellsc(1)*floor(dble(icelltmp-1)/ncellsc(1)+1.0d-8)
           jcell=jcelltmp-ncellsc(2)*floor(dble(jcelltmp-1)/ncellsc(2)+1.0d-8)
           
           do jj=1,celllistc(1,icell,jcell,kcell)
              jjcol=celllistc(jj+1,icell,jcell,kcell)                         
        !      if (iicol .eq. jjcol) cycle ! only do other colloids  because the trial col is not yet in the cell list          
              distij(:)=trialpos(:)-poscol(:,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbc2=distij2
              if (rbc2 .lt. rcut_blobcol2) then
                 pot_ene=pot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall

              endif              
           enddo
        enddo
     enddo
  enddo
  
  ! blob-wall energy
  if (ipcnew(3) .eq. 1) then
     rbw=trialpos(3)
     if (rbw .lt. rcut_blobwall) then
        pot_ene=pot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
     endif
  elseif (ipcnew(3) .eq. ncellsb(3)) then
     rbw=lbox(3)-trialpos(3)
     if (rbw .lt. rcut_blobwall) then
        pot_ene=pot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
     endif
  endif

end subroutine blob_energy_calc
!==========================================================================!
subroutine ext_col_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
              nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
              trialpos,iicol,pot_ene)
! THIS ROUTINE CALCULATES ONLY EXTERNAL ENERGY OF A GIVEN COLLOID, 
! i.e. ENERGY WITH BLOBS NOT BELONGING TO THE COLLOID
  
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),iicol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),trialpos(3),&
       rcut_blobblob,rcut_blobwall,rcol,poscol(3,maxncol)
  real*8, intent(out) ::pot_ene
  real*8 :: Eold,Enew,rnd,rnd2(2),rnd3(3),arg,volume,effdist,delene,&
       rcut_blobcol,rcut_blobcol2,rbc2
  integer :: jjcol,icelltmp,jcelltmp,kcelltmp,icell,jcell,kcell,& 
       iblob,jblob,jjblob,jjchain,cellrat(3),tmpcell(3),ii,jj,ipcnew(3),&
       ipcbc(3)
  real*8 :: distij(3),distij2,rbb,rbb2,rbw,&
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob
  real*8, parameter :: pi=3.14159265d0
  
  rcut_blobblob2=rcut_blobblob**2
  rcut_blobcol=rcut_blobwall+rcol
  rcut_blobcol2=rcut_blobcol**2
  eCS_blobwall=3.2*exp(-4.17*(rcut_blobwall-0.5d0))
  eCS_blobblob=1.75*exp(-0.8*rcut_blobblob2)
  
  pot_ene=0  

  ! blob col energy
  cellrat(:)=ceiling(rcut_blobcol/socb(:)) ! tells over how many blob cells one needs to go when checking for inserted colloid
  ipcnew=floor(trialpos(:)/socb(:))+1
  if (ipcnew(3) .lt. 1) ipcnew(3)=1
  if (ipcnew(3) .gt. ncellsb(3)) ipcnew(3)=ncellsb(3)
  do kcell=ipcnew(3)-cellrat(3),ipcnew(3)+cellrat(3)
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsb(3))) cycle 
     do jcelltmp=ipcnew(2)-cellrat(2),ipcnew(2)+cellrat(2)
        do icelltmp=ipcnew(1)-cellrat(1),ipcnew(1)+cellrat(1)
           
           icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
           jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
           
           do jj=1,celllistb(1,icell,jcell,kcell)
              jblob=celllistb(jj+1,icell,jcell,kcell)        
              jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)   
              
              if (iicol .eq. jjcol) cycle ! only get energy with other external blobs
              jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
              jjblob=mod(jblob,nblobsperchain)
              if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
              if (jjblob .eq. 1) write(*,*) 'WARNING: anchor should not be in the cell list '
      
              distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbc2=distij2
              if (rbc2 .lt. rcut_blobcol2) then
                 pot_ene=pot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall
              endif
          !    write(*,*) jjcol,jjchain,jjblob,distij2,rcut_blobcol2
           enddo
        enddo
     enddo
  enddo
 ! write(*,*) pot_ene, cellrat,iicol
end subroutine ext_col_energy_calc
!==========================================================================!
subroutine col_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
              nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
              trialpos,iicol,pot_ene)
! THIS ROUTINE CALCULATES TOTAL ENERGY OF A GIVEN COLLOID, 
  
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),iicol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),trialpos(3),&
       rcut_blobblob,rcut_blobwall,rcol,poscol(3,maxncol)
  real*8, intent(out) ::pot_ene
  real*8 :: Eold,Enew,rnd,rnd2(2),rnd3(3),arg,volume,effdist,delene,&
       rcut_blobcol,rcut_blobcol2,rbc2
  integer :: jjcol,icelltmp,jcelltmp,kcelltmp,icell,jcell,kcell,& 
       iblob,jblob,jjblob,jjchain,cellrat(3),tmpcell(3),ii,jj,ipcnew(3),&
       ipcbc(3)
  real*8 :: distij(3),distij2,rbb,rbb2,rbw,&
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob
  real*8, parameter :: pi=3.14159265d0
  
  rcut_blobblob2=rcut_blobblob**2
  rcut_blobcol=rcut_blobwall+rcol
  rcut_blobcol2=rcut_blobcol**2
  eCS_blobwall=3.2*exp(-4.17*(rcut_blobwall-0.5d0))
  eCS_blobblob=1.75*exp(-0.8*rcut_blobblob2)
  
  pot_ene=0  

  ! blob col energy
  cellrat(:)=ceiling(rcut_blobcol/socb(:)) ! tells over how many blob cells one needs to go when checking for inserted colloid
  ipcnew=floor(trialpos(:)/socb(:))+1
  if (ipcnew(3) .lt. 1) ipcnew(3)=1
  if (ipcnew(3) .gt. ncellsb(3)) ipcnew(3)=ncellsb(3)
  do kcell=ipcnew(3)-cellrat(3),ipcnew(3)+cellrat(3)
     if ((kcell .lt. 1) .or. (kcell .gt. ncellsb(3))) cycle 
     do jcelltmp=ipcnew(2)-cellrat(2),ipcnew(2)+cellrat(2)
        do icelltmp=ipcnew(1)-cellrat(1),ipcnew(1)+cellrat(1)
           
           icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
           jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
           
           do jj=1,celllistb(1,icell,jcell,kcell)
              jblob=celllistb(jj+1,icell,jcell,kcell)        
              jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)   
              
  !            if (iicol .eq. jjcol) cycle ! only get energy with other external blobs
              jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
              jjblob=mod(jblob,nblobsperchain)
              if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
              if (jjblob .eq. 1) write(*,*) 'WARNING: anchor should not be in the cell list '
      
              distij(:)=trialpos(:)-posblob(:,jjblob,jjchain,jjcol)
              distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
              distij2=sum(distij*distij)
              rbc2=distij2
              if (rbc2 .lt. rcut_blobcol2) then
                 pot_ene=pot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall
              endif
          !    write(*,*) jjcol,jjchain,jjblob,distij2,rcut_blobcol2
           enddo
        enddo
     enddo
  enddo
 ! write(*,*) pot_ene, cellrat,iicol
end subroutine col_energy_calc
!==========================================================================!
subroutine tot_ene_calc(posblob,poscol,lbox,nrec,ncol,maxncol,nchainspercol,& 
              nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,tot_ene)
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),ncol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),&
       rcut_blobblob,rcut_blobwall,rcol,poscol(3,maxncol)
  real*8, intent(out) ::tot_ene
  real*8 :: Eold,Enew,rnd,rnd2(2),rnd3(3),arg,volume,effdist,delene,&
       rcut_blobcol,rcut_blobcol2,rbc2
  integer :: jjcol,iicol,icelltmp,jcelltmp,kcelltmp,icell,jcell,kcell,& 
       iblob,jblob,jjblob,iiblob,jjchain,iichain,cellrat(3),tmpcell(3),ii,jj,ipcnew(3),&
       ipcbc(3)
  real*8 :: distij(3),distij2,rbb,rbb2,rbw,&
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob
  real*8, parameter :: pi=3.14159265d0
  
  rcut_blobblob2=rcut_blobblob**2
  rcut_blobcol=rcut_blobwall+rcol
  rcut_blobcol2=rcut_blobcol**2
  eCS_blobwall=3.2*exp(-4.17*(rcut_blobwall-0.5d0))
  eCS_blobblob=1.75*exp(-0.8*rcut_blobblob2)
  
  tot_ene=0  

 
  do iicol=1,ncol
     do iichain=1,nchainspercol
        do iiblob=2,nblobsperchain ! first blob is just an anchor
           ! blob blob energy
           
           iblob=(iicol-1)*nchainspercol*nblobsperchain+(iichain-1)*nblobsperchain+iiblob
        
           do kcell=ipcb(3,iblob)-1,ipcb(3,iblob)
              if ((kcell .lt. 1) .or. (kcell .gt. ncellsb(3))) cycle 
              do jcelltmp=ipcb(2,iblob)-1,ipcb(2,iblob)+1
                 do icelltmp=ipcb(1,iblob)-1,ipcb(1,iblob)+1
                    if ((kcell.eq.ipcb(3,iblob)) .and. ((jcelltmp.eq.(ipcb(2,iblob)+1)).or. &
                   ((jcelltmp.eq.ipcb(2,iblob)).and.(icelltmp.eq.(ipcb(1,iblob)+1))))) cycle
                    icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
                    jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
                    
                    do jj=1,celllistb(1,icell,jcell,kcell)
                       jblob=celllistb(jj+1,icell,jcell,kcell)
                       if (all((/icell,jcell,kcell/).eq.ipcb(1:3,iblob)).and. & 
                            (iblob .le. jblob)) cycle
 
                       jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)              
                     
                       jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
                       jjblob=mod(jblob,nblobsperchain)
                       if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
                       !     write(*,*) jjcol,jjchain,jjblob
                       distij(:)=posblob(:,iiblob,iichain,iicol)-posblob(:,jjblob,jjchain,jjcol)
                       distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                       distij2=sum(distij*distij)
                       rbb2=distij2
                       if (rbb2 .lt. rcut_blobblob2) then
                          tot_ene=tot_ene + 1.75*exp(-0.8*rbb2)-eCS_blobblob
                       endif
                    enddo
                 enddo ! icell
              enddo ! jcell
           enddo ! kcell

           ! harmonic energy          
           distij(:)=posblob(:,iiblob,iichain,iicol)-posblob(:,iiblob-1,iichain,iicol)
           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
           distij2=sum(distij*distij)
           if (iiblob .eq. 2) then ! blob-anchor energy
              tot_ene=tot_ene+0.75*distij2
           else
              rbb=sqrt(distij2)  
              tot_ene=tot_ene + 0.534*(rbb-0.73)**2                 
           endif

           ! blob col energy
           ipcbc(:)=floor(posblob(:,iiblob,iichain,iicol)/socc(:))+1  ! tells in which col cell the blob is           
     
           do kcell=ipcbc(3)-1,ipcbc(3)+1
              if ((kcell .lt. 1) .or. (kcell .gt. ncellsc(3))) cycle 
              do jcelltmp=ipcbc(2)-1,ipcbc(2)+1
                 do icelltmp=ipcbc(1)-1,ipcbc(1)+1
                
                    icell=icelltmp-ncellsc(1)*floor(dble(icelltmp-1)/ncellsc(1)+1.0d-8)
                    jcell=jcelltmp-ncellsc(2)*floor(dble(jcelltmp-1)/ncellsc(2)+1.0d-8)
                    
                    do jj=1,celllistc(1,icell,jcell,kcell)
                       jjcol=celllistc(jj+1,icell,jcell,kcell)
           !            if (jjcol .ne. iicol) cycle
                       distij(:)=posblob(:,iiblob,iichain,iicol)-poscol(:,jjcol)
                       distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
                       distij2=sum(distij*distij)
                       rbc2=distij2
                       if (rbc2 .lt. rcut_blobcol2) then
                          tot_ene=tot_ene + 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0)) - eCS_blobwall
!!$                          write(*,*)
!!$                          write(*,*) ' tot_ene: ', posblob(:,iiblob,iichain,iicol),poscol(:,jjcol)
!!$                          write(*,*) 3.2*exp(-4.17*(sqrt(rbc2)-rcol-0.5d0))
!!$                          write(*,*) iicol,iichain,iiblob,jjcol, jj, icell,jcell,kcell,tot_ene
!!$                          write(*,*) socc, lbox,ncellsc
                       endif
                       
                    enddo
                 enddo
              enddo
           enddo ! kcell
           
           ! blob-wall energy
           if (ipcb(3,iblob) .eq. 1) then
              rbw=posblob(3,iiblob,iichain,iicol)
              if (rbw .lt. rcut_blobwall) then
                 tot_ene=tot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
              endif
           elseif (ipcb(3,iblob) .eq. ncellsb(3)) then
              rbw=lbox(3)-posblob(3,iiblob,iichain,iicol)
              if (rbw .lt. rcut_blobwall) then
                 tot_ene=tot_ene + 3.2*exp(-4.17*(rbw-0.5d0)) - eCS_blobwall
              endif
           endif

        enddo ! nblobsperchain
     enddo ! nchainspercol
  enddo ! ncol
  
end subroutine tot_ene_calc
! ========================================================================== !
subroutine find_rec_homie(trialpos,nchainspercol,nblobsperchain,lbox,&
     posrec,nrec,celllistr,ipcr,ncellsr,socr,mnpic,rcut_ligrec,&
     recboundto,nrechomies,whichhomies)

  implicit none
  integer, intent(in) :: nrec,nblobsperchain,ncellsr(2),mnpic,nchainspercol,&
       celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec),recboundto(2,nrec)
  real*8, intent(in) :: posrec(2,nrec),lbox(3),socr(2),rcut_ligrec,trialpos(3)
  integer, intent(out) :: nrechomies,whichhomies(mnpic)
  real*8 :: rnd,rnd2(2),rnd3(3),arg,volume,effdist,pot_ene
  integer :: ii,jj,kk,icell,jcell,kcell,icelltmp,jcelltmp,kcelltmp,jrec,ipcnew(2)
  real*8 :: distij(3),distij2,rbr2,rcut_ligrec2
  
  real*8, parameter :: pi=3.14159265d0
  rcut_ligrec2=rcut_ligrec**2
  nrechomies=0  
  whichhomies(:)=0
  
  ipcnew=floor(trialpos(1:2)/socr(1:2))+1
  do icelltmp=ipcnew(1)-1,ipcnew(1)+1
     do jcelltmp=ipcnew(2)-1,ipcnew(2)+1         
        icell=icelltmp-ncellsr(1)*floor(dble(icelltmp-1)/ncellsr(1)+1.0d-8)
        jcell=jcelltmp-ncellsr(2)*floor(dble(jcelltmp-1)/ncellsr(2)+1.0d-8)
        
        do jj=1,celllistr(1,icell,jcell)
           jrec=celllistr(jj+1,icell,jcell)
           if (recboundto(1,jrec) .ne. 0) cycle  ! receptor is not free
           distij(1:2)=trialpos(1:2)-posrec(1:2,jrec)
           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))              
           distij(3)=trialpos(3)
           distij2=sum(distij*distij)
           rbr2=distij2
           if (rbr2 .lt. rcut_ligrec2) then
              nrechomies=nrechomies+1
              whichhomies(nrechomies)=jrec
           endif
        enddo
     enddo
  enddo

end subroutine find_rec_homie
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
subroutine make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,nrec,&
     poscol,posblob,posrec,lbox,&
     celllistc,ipcc,celllistb,ipcb,celllistr,ipcr,ncellsc,ncellsb,ncellsr,socc,socb,socr,mnpic)
 
  implicit none
  integer, intent(in) :: ncellsc(3),ncellsb(3),ncellsr(2),nrec,nchainspercol,&
       nblobsperchain,maxnblob,mnpic,ncol,maxncol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),socr(2),posrec(2,nrec),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),poscol(3,maxncol)
  integer, intent(out) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol)
  integer, intent(out) :: celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob)
  integer, intent(out) :: celllistr(mnpic,ncellsr(1),ncellsr(2)), ipcr(3,nrec)
  integer ::  ii,jj,kk,nblob,iblob,icol

  celllistc(:,:,:,:)=0
  ipcc(:,:)=0
  celllistb(:,:,:,:)=0
  ipcb(:,:)=0
  celllistr(:,:,:)=0
  ipcr(:,:)=0
  
  nblob=ncol*nchainspercol*nblobsperchain
  do kk=1, nrec
     ipcr(1:2,kk)=floor(posrec(:,kk)/socr(:))+1
     !write(*,*) poscol(:,kk),'   ', ipcc(1:3,kk),'  ',socc
     !write(*,*)
  enddo
  do ii=1,nrec
     celllistr(1, ipcr(1,ii), ipcr(2,ii)) = & 
          celllistr(1, ipcr(1,ii), ipcr(2,ii)) + 1
     celllistr((celllistr(1,ipcr(1,ii),ipcr(2,ii))+1),ipcr(1,ii),ipcr(2,ii))=ii
     
     ipcr(3,ii)=celllistr(1, ipcr(1,ii), ipcr(2,ii))
  enddo
!!$  
  do ii=1,ncol
     ipcc(1:3,ii)=floor(poscol(:,ii)/socc(:))+1        
     do jj=1,nchainspercol
        do kk=2,nblobsperchain ! don't add anchors into the cell list
           iblob=(ii-1)*nblobsperchain*nchainspercol+(nblobsperchain)*(jj-1)+kk
           ipcb(1:3,iblob)=floor(posblob(:,kk,jj,ii)/socb(:))+1
           
           if (ipcb(3,iblob) .lt. 1) ipcb(3,iblob)=1
           if (ipcb(3,iblob) .gt. ncellsb(3)) ipcb(3,iblob)=ncellsb(3)
           
        enddo
     enddo
  enddo
 
   do ii=1,ncol
     celllistc(1, ipcc(1,ii), ipcc(2,ii), ipcc(3,ii)) = & 
     celllistc(1, ipcc(1,ii), ipcc(2,ii), ipcc(3,ii)) + 1
     celllistc((celllistc(1,ipcc(1,ii),ipcc(2,ii),ipcc(3,ii))+1),ipcc(1,ii),ipcc(2,ii),ipcc(3,ii))=ii    
     ipcc(4,ii)=celllistc(1,ipcc(1,ii), ipcc(2,ii), ipcc(3,ii))
  enddo

  
  do icol=1,ncol
     do jj=1,nchainspercol
        do kk=2,nblobsperchain  ! don't add the anchor (kk=1)
           ii=(icol-1)*nblobsperchain*nchainspercol+(nblobsperchain)*(jj-1)+kk
     
           celllistb(1, ipcb(1,ii), ipcb(2,ii), ipcb(3,ii)) = & 
                celllistb(1, ipcb(1,ii), ipcb(2,ii), ipcb(3,ii)) + 1
           celllistb((celllistb(1,ipcb(1,ii),ipcb(2,ii),ipcb(3,ii))+1),ipcb(1,ii),ipcb(2,ii),ipcb(3,ii))=ii    
           ipcb(4,ii)=celllistb(1,ipcb(1,ii), ipcb(2,ii), ipcb(3,ii))
        enddo
     enddo
  enddo
  
end subroutine make_cell_list
! ========================================================================== ! 
subroutine update_cell_list(parti,maxnpart,ipcnew,ipc,celllist,ncells,mnpic)
  implicit none
  integer, intent(in) ::  parti,maxnpart,mnpic,ncells(3),ipcnew(3)
  integer, intent(inout) :: ipc(4,maxnpart), celllist(mnpic,ncells(1),&
       ncells(2),ncells(3))
  integer :: jj, ipcold(3)
  jj=ipc(4, parti)
  ipcold(:)=ipc(1:3,parti)   
  
  celllist(jj+1, ipcold(1), ipcold(2), ipcold(3)) = &
       celllist(celllist(1,ipcold(1),ipcold(2),ipcold(3))+1, ipcold(1),ipcold(2),ipcold(3))
  ipc(4,celllist(jj+1,ipcold(1),ipcold(2),ipcold(3))) = jj   
  
  celllist(1,ipcold(1),ipcold(2),ipcold(3)) = celllist(1,ipcold(1),ipcold(2),ipcold(3))-1 
  celllist(celllist(1,ipcold(1),ipcold(2),ipcold(3))+2,ipcold(1),ipcold(2),ipcold(3))=0   

  celllist(celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))+2, ipcnew(1),ipcnew(2),ipcnew(3)) =&
       parti
  celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))=celllist(1,ipcnew(1),ipcnew(2),ipcnew(3))+1
  
  ipc(1:3,parti)=ipcnew(:)
  ipc(4, parti)=celllist(1,ipcnew(1),  ipcnew(2), ipcnew(3))
end subroutine update_cell_list
! ========================================================================== !
subroutine check_bonds(nrec,ncol,maxncol,maxnblob,nchainspercol,& 
     nblobsperchain,ligboundto,recboundto,posblob,posrec,rcut_ligrec,lbox)
  implicit none
  
  integer, intent(in) :: nrec,maxncol,nchainspercol,nblobsperchain,maxnblob,&
       ligboundto(nchainspercol,maxncol),recboundto(2,nrec),ncol 
  real*8, intent(in) :: posblob(3,nblobsperchain,nchainspercol,maxncol),posrec(2,nrec),&
       rcut_ligrec, lbox(3)
  integer :: ii, jj, kk, iiblob, iichain,jjblob,kkblob,&
       llblob,jjrec,kkrec,irec,jjchain,kkchain,jjcol,iicol,kkcol
  real*8 :: distij(3),distij2,rcut_ligrec2
  
  rcut_ligrec2=rcut_ligrec**2

  ! do lig-receptor bond check on receptor side
  do irec=1,nrec
     if (recboundto(1,irec) .gt. 0) then
        jjchain=recboundto(1,irec)
        jjcol=recboundto(2,irec)
        kkrec=ligboundto(jjchain,jjcol)
        if (irec .ne. kkrec) then
           write(*,*) 'recside bonding error!!  rec i  points to ligand j in colloid j :', &
                irec,jjchain,jjcol
           write(*,*) ' but lig j in  j points to rec k : ',jjchain,jjcol,kkrec
        endif   
        
        distij(1:2)=posrec(:,irec)-posblob(1:2,nblobsperchain,jjchain,jjcol)
        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
        distij(3)=posblob(3,nblobsperchain,jjchain,jjcol)
        distij2=sum(distij*distij)
        if (distij2 .gt. rcut_ligrec2) then
       !    write(*,*) 'WARNING: lig-rec bond is present but distance is too great'
       !    write(*,*) ' iicol ',jjcol ,'jjchain ',jjchain ,' kkrec ',irec
       !    write(*,*) 'poslig ', posblob(1:3,nblobsperchain,jjchain,jjcol), 'posrec ', posrec(:,irec)
       !    write(*,*) ' distij2 =', distij2   
        endif
      
     endif
  enddo
  
    ! do lig-receptor bond check on ligand side
  do iicol=1,ncol
     do iichain=1,nchainspercol
        jjrec=ligboundto(iichain,iicol)
        if (jjrec .gt. 0) then
           kkchain=recboundto(1,jjrec)
           kkcol=recboundto(2,jjrec)
           if (any((/iichain,iicol/) .ne. (/kkchain,kkcol/) )) then
              write(*,*) 'ligside bonding error!!  lig i on chain j  points to rec k :', iichain,iicol,jjrec
              write(*,*) ' but rec k points to lig i on chain j : ',jjrec,kkchain,kkcol
           endif
        endif
     enddo
  enddo

end subroutine check_bonds
! ========================================================================== !
subroutine check_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,lbox,poscol,posblob, &
     mnpic,ncellsc,socc,celllistc,ipcc,ncellsb,socb,celllistb,ipcb)
  implicit none  
  integer, intent(in) :: ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,mnpic,&
       ncellsc(3),celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol),&
       ncellsb(3),celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),&
       ipcb(4,maxnblob)
  real*8, intent(in) :: lbox(3),socc(3),socb(3),poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol)
  integer :: ii,jj,kk,icol,icolspec,jcolspec,jcol,ipctrial(3),ncolincell,&
       icell,jcell,kcell,iserialblob,iblob,ichain,nblobincell
    
  do icol=1,ncol
     ipctrial(1:3)=floor(poscol(:,icol)/socc(:))+1
     
     if (any(ipctrial(1:3) .ne. ipcc(1:3,icol))) then
        write(*,*) 'WARNING: check_cell_list failed..... ipcc not good with poscol'
        write(*,*) 'ipcc= ',ipcc(1:3,icol), 'ipctrial = ', ipctrial(1:3)
        write(*,*) 'poscol(:,icol) ', poscol(:,icol)
     endif
     
     ! check celllist array
     if ( celllistc(1+ipcc(4,icol), ipcc(1,icol), ipcc(2,icol), &
          ipcc(3,icol)) .ne. icol) then
        write(*,*)'WARNING: check_cell_list failed..... '
        write(*,*)  'celllist in cell gives particle', celllistc(1+ipcc(4,icol), ipcc(1,icol), &
             ipcc(2,icol), ipcc(3,icol)) ,' but particle is ',  icol
     endif
     
     ! check blobs
     do ichain=1,nchainspercol
        do iblob=2,nblobsperchain ! blob #1 is anchor and is not in the cell lists
           iserialblob=(icol-1)*nblobsperchain*nchainspercol+(ichain-1)*nblobsperchain+iblob
           ipctrial(1:3)=floor(posblob(:,iblob,ichain,icol)/socb(:))+1
           ipctrial(3)=max(ipctrial(3),1)
           if (any(ipctrial(1:3) .ne. ipcb(1:3,iserialblob))) then
              write(*,*) 'WARNING: check_cell_list failed..... ipcb not good with posblob'
              write(*,*) 'icol ',icol,', ichain ',ichain,', iblob',iblob, ', iserialblob',iserialblob
              write(*,*) 'ipcb= ',ipcb(1:3,iserialblob), 'ipctrial = ', ipctrial(1:3)
              write(*,*) 'posblob(:,iblob,ichain,icol) ', posblob(:,iblob,ichain,icol)

           !   write(*,*) 'ipcb array:', ipcb
                ! write(*,*) celllistb

           endif
           
           ! check celllist array
           if ( celllistb(1+ipcb(4,iblob), ipcb(1,iblob), ipcb(2,iblob), &
                ipcb(3,iblob)) .ne. iblob) then
              write(*,*)'WARNING: check_cell_list failed..... '
              write(*,*)  'celllist in cell gives blob', celllistb(1+ipcb(4,iblob), ipcb(1,iblob), &
                   ipcb(2,iblob), ipcb(3,iblob)) ,' but blob is ',  iblob
           endif
        enddo ! iblob
     enddo ! ichain
  enddo ! ncol
  
  ! check col cells
  do icell=1,ncellsc(1)
     do jcell=1,ncellsc(2)
        do kcell=1,ncellsc(3)
           ncolincell=celllistc(1, icell,jcell,kcell)
           do ii=2,ncolincell+1            
              icol=celllistc(ii, icell,jcell,kcell) ! col serial number, convert to ispecie and icol
              
              if (any(ipcc(1:3,icol).ne.(/icell,jcell,kcell/))) then
                 write(*,*) 'WARNING: check_cell_list failed.....  '
                 write(*,*)  'ipcc gives cell ',ipcc(1:3,icol)  ,' but cell is ',&
                      (/icell,jcell,kcell/)
                 
              endif
              if (ipcc(4,icol) .ne. ii-1)  then
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
  
 ! check blob  cells
  do icell=1,ncellsb(1)
     do jcell=1,ncellsb(2)
        do kcell=1,ncellsb(3)
           nblobincell=celllistb(1, icell,jcell,kcell)
           do ii=2,nblobincell+1            
              iserialblob=celllistb(ii, icell,jcell,kcell) ! blob serial number, convert to iblob,ichain,icol
              icol=ceiling(dble(iserialblob)/(nblobsperchain*nchainspercol))
              ichain=ceiling(dble(iserialblob-(icol-1)*nblobsperchain*nchainspercol)/nblobsperchain)
              iblob=iserialblob-(icol-1)*nblobsperchain*nchainspercol-(ichain-1)*nblobsperchain
              
              if (any(ipcb(1:3,iserialblob).ne.(/icell,jcell,kcell/))) then
                 write(*,*) 'WARNING: check_cell_list failed.....  '
                 write(*,*)  'ipcc gives cell ',ipcb(1:3,iserialblob)  ,' but cell is ',&
                      (/icell,jcell,kcell/)
                 
              endif
              if (ipcb(4,iserialblob) .ne. ii-1)  then
                 write(*,*)'WARNING: check_cell_list failed..... ',&
                      'particle position in given cell array is wrong'
                 
              endif
           enddo
           do ii=nblobincell+2,mnpic
              if (celllistb(ii, icell,jcell,kcell) .ne. 0)  write(*,*)'WARNING: check_cell_list failed..... ',&
                   'cell list has non-zero entries beyond the number of blobs in the cell'
           enddo
        enddo
     enddo
  enddo

end subroutine check_cell_list
! ========================================================================== !
subroutine output_conf_wl(poscol,posblob,posrec,maxncol,nchainspercol,nblobsperchain,ncol,&
     nrec,recboundto,mu,nout,icycle,outfilename,nrxspec,irecspec,ichainspec,&
     cWL_min,cWL_max,psiWL,hWL,doWL_flag)  
  implicit none
  
  integer, intent(in) ::  maxncol,ncol,nrec,nchainspercol,nblobsperchain,recboundto(2,nrec),&
       nrxspec,irecspec(nrec),ichainspec(nchainspercol)
  integer*8, intent(in) :: nout, icycle,hWL(nchainspercol+1)
  real*8, intent(in) :: poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol),&
       posrec(2,nrec), mu, psiWL(nchainspercol+1),cWL_min,cWL_max
  character*50, intent(in) :: outfilename
  logical, intent(in) :: doWL_flag
  character*50 :: filename
  character*4 :: char1,char2
  character*6 :: char3,char4
  integer :: irec,icol,ichain,iblob,ii,ibin
  real*8 :: binsize,freeE(nchainspercol+1)
  ! VMD OUTPUT FILE
  write(char1,'(I4)') int(icycle/nout)+1000 
  write(char2,'(I4)') nint(mu)+2000
  write(char3,'(I6)') ncol+100000
  write(char4,'(I6)') nrec+100000
  !write(*,*) 'OK'
  filename=trim(outfilename)//'WLbonds_nrec'//char4//'_ncol'//char3//'_out'//char1//'.xyz'


  freeE(:)=psiWL(:)-log(dble(hWL(:)))  ! get free energy
  freeE(:)=freeE(:)-freeE(1)  ! normalise by 1st bin

  open(unit=222,file=filename)

  if (.not. doWL_flag) then ! only write the free energy datafiles if WL has finished.
     write(*,*) 'bFreeE ',-log(sum(exp(-freeE(2:nchainspercol+1))))   ! ibin=ibond
     do ibin=1,nchainspercol+1
        write(222,*) ibin-1,'  ',freeE(ibin)   ! ibin=ibond
     enddo
  endif
  
  close(222) ! end WL writing
  filename=trim(outfilename)//'conf_nrec'//char4//'_ncol'//char3//'_out'//char1//'.xyz'
  open(unit=222,file=filename)
  
  write(222,*) ncol*nblobsperchain*nchainspercol +nrec
  write(222,*) ncol, nrec, 'Hairy nanoparticle adsorptopn'
  
  do icol=1,ncol         
     write(222,"(A,F8.3,F8.3,F8.3)") 'Xe  ', poscol(:,icol)
     do ichain=1,nchainspercol
        write(222,"(A,F8.3,F8.3,F8.3)") 'H ', posblob(:,1,ichain,icol) ! chain anchor on particle
        do iblob=2,nblobsperchain-1
           write(222,"(A,F8.3,F8.3,F8.3)") 'O  ', posblob(:,iblob,ichain,icol)  ! inert blobs         
        enddo
        if (ichainspec(ichain) .eq. 3) then
           write(222,"(A,F8.3,F8.3,F8.3)") 'N', posblob(:,nblobsperchain,ichain,icol)   ! ligand   
        else   
           write(222,"(A,F8.3,F8.3,F8.3)") 'F', posblob(:,nblobsperchain,ichain,icol)   ! ligand   
        endif
     enddo
  enddo

  ! fill up so that you can run simulations
 ! do icol=ncol+1,maxncol
 !    do ii=1,nblobsperchain*nchainspercol+1
 !       write(222,"(A,F8.3,F8.3,F8.3)") 'H  ', 0, 0, 0
 !    enddo
 ! enddo

  
  
  do irec=1,nrec
!     if (recboundto(1,irec) .gt. 0) then !this receptor is bound to someone
!        write(222,"(A,F8.3,F8.3,F8.3)") 'F ', posrec(:,irec) , 0.0
!     else
!        write(222,"(A,F8.3,F8.3,F8.3)") 'H ', posrec(:,irec) , 0.0
!     endif
     if (irecspec(irec) .eq. 1) then !this receptor is specie 1
        write(222,"(A,F8.3,F8.3,F8.3)") 'N ', posrec(:,irec) , 0.0
     else
        write(222,"(A,F8.3,F8.3,F8.3)") 'F ', posrec(:,irec) , 0.0
     endif


  enddo
  
  close(222)
  ! END VMD OUTPUT FILE
  
end subroutine output_conf_wl
