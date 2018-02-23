! MC Hairy Nanoparticle code, aka multi-multivalency + mobile ligands + 
! updated from adsorbtion MMV and adapted to colloid colloid interactions, only 2 colloids
! by tc387@cam in Dec 2017
! Wang Landau

!UPDATE
! by tc387@cam in Feb 2018
! added bonded anchors; anchors can form bonds


! CELL LISTS:  celllist(:,:,:,:) 1st dim tells how many (max mnpic) and which col are in the cell icol*icolspec, 2nd , 3rd and 4th dim are ii,jj,kk indeces of the cell 
! ipc(:,:,:) ith colloid cell - tells in which cell the g`iven colloid is and on which position in a celllist chain in that cell. ipc(4,maxncol,ncolspec)
! rxcolgeo = rx specie on colloid geometry, last matrix in input file
!posolrx = position of colloid and rx that blong to it 
! bondrx(3,nrxpercol,maxncol,ncolspec) tells with who (specie,col,position) is particulat rx bound to , position=0 means anchors then(rxspecie,iancm,0)
program mc2cbWL
implicit none
  
! TRY NOT USING GLOBAL VARIABLES BUT DEFINE EVERYTHING IN IT'S OWN SCOPE !
! USE INTENT IN/OUT IN ALL SUBROUTINES & FUNCTIONS (for easier debugging)!
! REMEMBER COLUMN MAJOR FORM !
integer :: ncellsc(3),ncellsb(3),ncellsr(2),mnpic,maxncol,nchainspercol,nfreecol,&
     ipcnew(3),icol,ii,jj,kk,jjcycle,excacc,movebacc,movecacc,bondacc,ncol,nrec,&
     nimp_anc,nblobsperchain,nblob,iiblob,jjblob,kfactor,maxnblob,nrxspec,&
     jjchain,jjcol,nWL_points,WLibin,WLmaxnbonds,nbonds
integer*8 :: icycle,ncycles,nout,totidm,avencol,avnbpci,avnbcol,avnbcolii,&
     avnbpcind, nequil,Xseed,WLacc
real*8 :: lbox(3),rcol,socc(3),socb(3),socr(2),rcut_liglig,rcut_blobwall,&
     rcut_blobblob,rcut_global,rnd,rnd2(2),rnd3(3),ligrec_ene,rrec,&
     time1,time2,tot_ene,fracgcmc,tot_ene_old,max_hop_col,fracbm,&
     max_rot_col,max_hop_blob,dimp_anc,mu,activity, rcc_cut,  maxdelRcol2,avnbpc,&
     rWL,rWL_max,rWL_min,fWL_start,fWL_stop,fWL_reduce,fWL,hWL_conv,WLbinsize,&
     max_colzWL
integer,allocatable :: celllistc(:,:,:,:),ipcc(:,:),celllistb(:,:,:,:),ipcb(:,:),&
     celllistr(:,:,:),ipcr(:,:),ligboundto(:,:,:),recboundto(:,:),nrecs(:),&
     irecspec(:),ichainspec(:,:),nchainsperspec(:,:)
integer*8, allocatable :: hWL(:)
real*8, allocatable :: posrec(:,:),posblob(:,:,:,:),poscol(:,:),RXCOLGEO(:,:,:),&
     RXBONDENE(:,:),psiWL(:),histWL(:)
logical ::  read_init_conf=.false., & 
     mobile_ligands,read_only_anchor_pos, random_anconcolpos, doWL_flag
real*8, parameter :: pi=3.14159265d0
character*50 :: outfilename, init_conf_filename1, init_conf_filename2

open(unit=101, file='input-parwlb.dat',status='old')
read(101,*)
read(101,*)
read(101,*) (lbox(ii),ii=1,3)
rcc_cut = lbox(1) / 2.0 ! colloid-colloid cutoff distance
read(101,*) nrxspec
read(101,*) maxncol
read(101,*) nchainspercol
allocate(ichainspec(nchainspercol, 2),nchainsperspec(nrxspec,2)) ! for both colloids independently
read(101,*) (nchainsperspec(ii, 1),ii=1,nrxspec)
read(101,*) (nchainsperspec(ii, 2),ii=1,nrxspec)

do icol=1,2
   kk=0 ! assign species to chains
   do ii=1,nrxspec
      do jj=1,nchainsperspec(ii,icol)
         kk=kk+1
         ichainspec(kk,icol)=ii
      enddo
   enddo
enddo
read(101,*) nblobsperchain
read(101,*) ncycles
read(101,*) nout
read(101,*) !---------------------------------------------
read(101,*) rcol
read(101,*) rcut_liglig
read(101,*) fracbm
read(101,*) max_hop_blob, max_hop_col, max_rot_col
read(101,*) mnpic
read(101,*) Xseed
read(101,*) !---------------------------------------------
read(101,*) mobile_ligands
read(101,*) nequil 
read(101,*) read_init_conf
read(101,*) init_conf_filename1
read(101,*) init_conf_filename2
read(101,*) outfilename
read(101,*) !--------------------------------------------
read(101,*) max_colzWL
read(101,*) WLmaxnbonds  ! = max number of bonds - 1
read(101,*) fWL_start,fWL_stop,fWL_reduce
read(101,*) hWL_conv
read(101,*) !--------------------------------------------
read(101,*) ! RX INTERACTION MATRIX
allocate(RXBONDENE(nrxspec,nrxspec))
do ii=1,nrxspec
   read(101,*) (RXBONDENE(jj,ii),jj=1,nrxspec)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!        CHECKS         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (maxncol .ne. 2) then
   write(*,*) 'ERROR: current code works for only 2 collods, exiting.. '
endif
if ((lbox(2) .lt. 2*(rcol+nblobsperchain+10)) .or. (lbox(3) .lt. 2*(rcol+nblobsperchain+10))) then
   write(*,*) 'WARNING: lateral box dimensions too small, setting to: 2*(rcol+nblobsperchain+2) '
   
   lbox(2)=2*(rcol+nblobsperchain+10)
   lbox(3)=2*(rcol+nblobsperchain+10)
endif
if (lbox(1) .lt. 8*(rcol + nblobsperchain + 2)) then
   write(*,*) 'WARNING: box is probably too small: ', lbox(1), ' setting to 8*(rcol + nblobsperchain + 2)'
   lbox(1)=8*(rcol + nblobsperchain + 2)
   
endif

!!!!!!!!!!!!!!!!!!!!!!!!       END CHECKS         !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define usefull parameters
maxnblob=maxncol*nchainspercol*nblobsperchain ! maximum possible number of blobs in the system
rcut_blobblob=3.0d0
rcut_blobwall=2.0d0
rcut_global=max(rcut_blobblob,rcut_liglig)

socc(:)=2*rcol
ncellsc(:)=floor(lbox(:)/socc(:))
ncellsc=max(ncellsc,3)
socc(:)=lbox(:)/ncellsc(:) ! renormalize socc

ncellsb=floor(lbox(:)/rcut_blobblob)
ncellsb=floor(lbox(:)/rcut_global) ! added 12 Feb 2018
ncellsb=max(ncellsb,3)
socb(:)=lbox(:)/ncellsb(:) ! renormalize socb

kfactor=8
tot_ene=0
excacc=0
movebacc=0
movecacc=0
bondacc=0
totidm=0
avencol=0
avnbpcind=0
avnbpc=0
avnbpci=0
avnbcol=0
avnbcolii=0

nrec=0
ncol=0
nblob=0

maxdelRcol2 = (lbox(1)/4)**2
! WL CHECKS
rWL_min=max(rWL_min,2*rcol)
max_colzWL = max_colzWL+2*rcol
rWL_max=max_colzWL ! min(rWL_max, lbox(1)/4)
!WLbinsize = (rWL_max-rWL_min)/nWL_points
nWL_points = WLmaxnbonds + 1

! write all parameters to screen
write(*,*)
write(*,*) '================ MC MIPS SIM ================'
write(*,*) '============================================='
write(*,*) '================= INPUT PAR ================='
write(*,*) 'boxsize  ' ,(lbox(ii),ii=1,3)
write(*,*) '# of rx species  ', nrxspec
write(*,*) 'max # of colloids  ', maxncol
write(*,*) '# of chains per collod ',nchainspercol
write(*,*) 'colloid 1 # of chains per rx specie   ', nchainsperspec(:, 1)
write(*,*) 'colloid 1 # of chains per rx specie   ', nchainsperspec(:, 2)

write(*,*) 'colloid 1 rx specie per chain ', ichainspec(:,1)
write(*,*) 'colloid 2 rx specie per chain ', ichainspec(:,2)
write(*,*) '# of blobs per chain ', nblobsperchain
write(*,*) 'tot n cycles  ',ncycles
write(*,*) 'nout   ',nout
write(*,*) '---------------------------------------------'
write(*,*) 'col radius rcol  ', rcol
write(*,*) 'bond distance cut off ',rcut_liglig
write(*,*) 'fraction of bond create/destroy moves  ', fracbm
write(*,*) 'max hop blob, col, rot col  ',max_hop_blob, max_hop_col, max_rot_col
write(*,*) 'max # of particles in each cell+1  ',mnpic
write(*,*) 'size of cell -blobs:', socb
write(*,*) 'seed for rng  ', Xseed
write(*,*) '---------------------------------------------'
write(*,*) 'mobile ligands ', mobile_ligands 
write(*,*) 'annealing equilibration steps ',nequil
write(*,*) 'read initial conf', read_init_conf
write(*,*) 'init conf filename  1', init_conf_filename1
write(*,*) 'init conf filename  2', init_conf_filename2
write(*,*) 'outfilename  ', outfilename
write(*,*) '---------------------------------------------'
write(*,*) '  # of Wang Landau points  ',nWL_points 
write(*,*) 'maximum distance for nonbonded col  ', rWL_max
write(*,*) 'f parameters for WL psi update ', fWL_start, fWL_stop, fWL_reduce  
write(*,*) 'WL Histogram convergance criterion ',  hWL_conv
write(*,*) '---------------------------------------------'
! ! RX INTERACTION MATRIX
write(*,*) 'RX interaction matrix'
do ii=1,nrxspec
   write(*,*) (RXBONDENE(jj,ii),jj=1,nrxspec)
enddo

write(*,*) '============== END INPUT PAR ================'
write(*,*) '============================================='



! allocate big arrays
allocate(celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob))
allocate(celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol))
allocate(poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol)) 
allocate(ligboundto(2,nchainspercol,maxncol), posrec(2,nrec)) ! allocate dummy posrec
allocate(RXCOLGEO(3, nchainspercol, 2))
allocate(psiWL(nWL_points),hWL(nWL_points))
posrec=0

! RANDOM SEED
call seed_random_number(Xseed)

!------  WANG LANDAU INIT  ------!
fwl=0
psiWL(:)=0
hWL(:)=0
doWL_flag=.true.


call initial_conf(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
     rcol,poscol,posblob,lbox,mnpic,nrxspec,ichainspec,RXCOLGEO,random_anconcolpos,&
     read_init_conf,init_conf_filename1,init_conf_filename2)



call make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
    poscol,posblob,lbox,&
    celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,socc,socb,mnpic)

! insert the two colloids, 1st at the origin [0,0,0], 2nd at [lbox/2,0,0]
call insert_2col(lbox,ncol,maxncol,nchainspercol,nblobsperchain,maxnblob,rcol,&
     poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,celllistc,celllistb,ipcc,ipcb,&
     kfactor,ligboundto,rcut_blobblob,rcut_blobwall,RXCOLGEO,random_anconcolpos,tot_ene)

nbonds=0
! EQUILIBRATE
if (nequil .gt. 0) then
   write(*,*) ' =========  STARTING EQUILIBRATION  ========='
   do icycle=1, nequil
      do jjcycle=1, nchainspercol*nblobsperchain+1  ! so that there is approx 1 hop per colloid per cycle
         call mc_move2WL(lbox,ncol,maxncol,nchainspercol,nblobsperchain,maxnblob,rcol,&
              poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,socr,celllistc,&
              celllistb,ipcc,ipcb,ligboundto,rcut_blobblob,&
              rcut_blobwall,rcut_liglig,max_hop_blob,max_hop_col,max_rot_col,&
              movebacc,movecacc,tot_ene,RXBONDENE,nrxspec,ichainspec,.true., maxdelRcol2,&
              nWL_points,max_colzWL,psiWL,.true.,nbonds)
      enddo
   enddo
   write(*,*) ' ============   EQUIL FINISHED   ============'
   write(*,*) ' ============================================'      
endif
! MOVE Colloids closer together if necessary for WL

if (poscol(1,2) - poscol(1,1) .ge. rWL_max) then    
   posblob(1,:,:,2) = posblob(1,:,:,2) - (poscol(1,2)-poscol(1,1) - rWL_max+1.0d-3)
   poscol(1,2) = poscol(1,1) + rWL_max - 1.0d-3
   call make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
        poscol,posblob,lbox,&
        celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,socc,socb,mnpic)
endif

! save topology
icycle=0
call output_conf_lammps(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,&
     psiWL,hWL,lbox,ligboundto)

!------  WANG LANDAU INIT  ------!
fwl=fWL_start
psiWL(:)=0
hWL(:)=0
doWL_flag=.true.

call cpu_time(time1)
    
write(*,*) '================= START SIM ================='
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXX  MAIN CYCLE  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
do icycle=1,ncycles
   do jjcycle=1, nchainspercol*nblobsperchain+1  ! so that there is approx 1 hop per colloid per cycle
      
      call random_number(rnd)
      if (rnd .gt. (fracgcmc+fracbm)) then  !HOP COL/BLOB
         call mc_move2WL(lbox,ncol,maxncol,nchainspercol,nblobsperchain,maxnblob,rcol,&
              poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,socr,celllistc,&
              celllistb,ipcc,ipcb,ligboundto,rcut_blobblob,&
              rcut_blobwall,rcut_liglig,max_hop_blob,max_hop_col,max_rot_col,&
              movebacc,movecacc,tot_ene,RXBONDENE,nrxspec,ichainspec,mobile_ligands, &
              maxdelRcol2,nWL_points,max_colzWL,psiWL,.false.,nbonds)

         
      else   ! BOND CREATE/DESTROY
    !     call bond_create_destroy(posrec,posblob,lbox,nrec,ncol,maxncol,nchainspercol,& 
    !          nblobsperchain,celllistr,ipcr,ncellsr,socr,mnpic,rcut_liglig,bondacc,&
    !          ligboundto,recboundto,nrxspec,RXBONDENE,irecspec,ichainspec)

      endif
      totidm=totidm+1
      avencol=avencol+ncol

   ! do Wang Landau
      call WangLandau_bonds(nchainspercol,nbonds,psiWL,nWL_points,&
           hWL,fWL,fWL_stop,fWL_reduce,hWL_conv,icycle,WLacc,doWL_flag)
      
   enddo
   ! check bonds -- only for debugging
!!         call check_bonds_lig(ncol,maxncol,maxnblob,nchainspercol,& 
!!              nblobsperchain,ligboundto,posblob,rcut_liglig,lbox)   
!!   call check_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,lbox,&
!!        poscol,posblob,mnpic,ncellsc,socc,celllistc,ipcc,ncellsb,socb,celllistb,ipcb) 
      
   ! get the average number of bound ligands per chain
   ! only every 10 cycles, because it's expensive
   if (icycle/10 .eq. dble(icycle)/10) then
      ! get average number of bonds per bound colloid 
      avnbpci=0
      do jjcol=1,ncol
         do jjchain=1,nchainspercol
            if (ligboundto(1,jjchain,jjcol) .gt. 0) avnbpci=avnbpci+1
         enddo
      enddo
      avnbpc=avnbpc + dble(avnbpci) / ( nchainspercol * ncol )
      ! simply divide by 2 colloids
      avnbpcind = avnbpcind + 1
   endif  
   
   ! OUTPUT CONFIGURATION
   if (icycle/nout .eq. dble(icycle)/nout) then  
      call cpu_time(time2)
      tot_ene_old=tot_ene

      call tot_ene_calc(posblob,poscol,lbox,nrec,ncol,maxncol,nchainspercol,& 
           nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
           socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,ligboundto,tot_ene)
      write(*,*)
      write(*,"(A,I10,A,F10.4)") 'icycle',icycle,', exe time (s) ',time2-time1
      write(*,"(A,F11.6,A,F11.6,A,F11.6)") 'avncol:',dble(avencol)/dble(totidm),&
           '   av # boundcol:',dble(avnbcol*10)/dble(nout),&
           '   av # bonds_per_col:',avnbpc/dble(avnbpcind)
      write(*,"(A,F9.3,A,F9.3)")  '   tot ene old:', tot_ene_old,'   tot ene new:' , tot_ene
      write(*,"(A,F7.5,A,F7.5,A,F7.5)") 'movebacc:',dble(movebacc)/dble(totidm)/&
           (1-fracbm)*(nchainspercol*nblobsperchain+1)/dble(nchainspercol*nblobsperchain),&
           '   movecacc:',real(movecacc)/dble(totidm)/&
           (1-fracbm)*(nchainspercol*nblobsperchain+1)/5.0 ,  &
           '   bondacc:',dble(bondacc)/dble(totidm)/fracbm
      write(*,"(A,F10.8,A,F8.2,A,I5)") 'fWL: ',fWL,'   hWLrat: ', &
           dble(maxval(hWL))/minval(hWL),'   nbonds: ', nbonds
      excacc=0
      movebacc=0
      movecacc=0
      bondacc=0
      totidm=0
      avencol=0
      avnbcol=0
      avnbpc=0
      avnbpcind=0
   
      call output_conf_trj(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
           nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,&
           psiWL,hWL,lbox,ligboundto)
      call output_conf_wlbonds(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
           activity,nout,icycle,outfilename,nrxspec,ichainspec,&
           psiWL,hWL,doWL_flag,nWL_points)
      call output_conf_lammps(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,&
     psiWL,hWL,lbox,ligboundto)
      time1=time2
   endif
  
enddo !ncycles

write(*,*) '==================== END ====================='
close(101)
!close(909)

end program mc2cbWL
  
subroutine fbondene(ene,distij2)
  implicit none
  real*8, intent(in) :: distij2
  real*8, intent(out) :: ene
  real*8 :: kspringbond
  kspringbond=1.0
  ! simple harmonic function
  ene = kspringbond/2.0 * distij2
end subroutine fbondene


! ========================================================================== !
! ========================================================================== !
! ========================================================================== !
! ========================================================================== !
subroutine  WangLandau_bonds(nchainspercol,nbonds,psiWL,nWL_points,&
     hWL,fWL,fWL_stop,fWL_reduce,hWL_conv,icycle,WLacc,doWL_flag)
  
  
  implicit none
  integer, intent(in) :: nchainspercol,nbonds,nWL_points
  integer*8, intent(inout) :: WLacc
  integer*8, intent(in) :: icycle
  integer*8, intent(inout) :: hWL(nWL_points)
  real*8, intent(in) :: fWL_stop,fWL_reduce,hWL_conv
  real*8, intent(inout) :: psiWL(nWL_points),fWL
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

subroutine initial_conf(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
     rcol,poscol,posblob,lbox,mnpic,nrxspec,ichainspec,RXCOLGEO,random_anconcolpos,&
     read_init_conf,init_conf_filename1,init_conf_filename2)

  implicit none
  integer, intent(inout) :: ncol,maxncol,maxnblob,nchainspercol,&
       nblobsperchain,mnpic,nrxspec,ichainspec(nchainspercol)
  real*8, intent(inout) :: lbox(3),rcol,poscol(3,maxncol),&
       posblob(3,nblobsperchain,nchainspercol,maxncol), RXCOLGEO(3,nchainspercol,maxncol)
  logical, intent(in) :: read_init_conf
  logical, intent(out) :: random_anconcolpos
  character*50, intent(in) :: init_conf_filename1,init_conf_filename2
  real*8, parameter :: pi=3.14159265d0
  integer :: ii,jj,kk,seed,now(3),irec,jrec,nrectmp,ispec,nrecstmp(nrxspec),&
       nchainspercolread1, nchainspercolread2
  real*8 :: rnd3(3),distij(3),distij2,rnd2(2)
  logical :: overlap, twofiles
  character*2 :: atom

  poscol(:,:)=0
  posblob(:,:,:,:)=0
  !posrec(:,:)=0
  
  ! READ INITIAL CONF   ! NOT WORKING YET !!!!!!!!!!!!!!!!!!!
  if (read_init_conf) then
     random_anconcolpos=.false.
     twofiles=.true.
     if (init_conf_filename1 .eq. init_conf_filename2) twofiles=.false.
     open(unit=301, file=init_conf_filename1, status='old')
     
     if (twofiles) open(unit=302, file=init_conf_filename2, status='old')
     
     read(301,*) nchainspercolread1
     if (twofiles) then ;read(unit=302) nchainspercolread2
     else ; nchainspercolread2 = nchainspercolread1
     endif
     
     if (nchainspercolread1 .ne. nchainspercol) then
        write(*,*) 'nchainspercol must be equal to the number of entries in the input file 1!!'
        write(*,*) nchainspercolread1, '  exiting...'
        call exit()
     elseif (nchainspercolread2 .ne. nchainspercol) then
        write(*,*) 'nchainspercol must be equal to the number of entries in the input file 2!!'
        write(*,*) nchainspercolread2, '  exiting...'
        call exit()
     endif
     read(301,*) 
     if (twofiles) read(unit=302) ! empty line
     
     do ii=1,nchainspercol
        read(301,*) atom, (RXCOLGEO(kk,ii,1), kk=1,3)
        if (twofiles) then
           read(302,*) atom, (RXCOLGEO(kk,ii,2), kk=1,3)
        else
           RXCOLGEO(:,ii,2)=RXCOLGEO(:,ii,1)
        endif
     enddo
     close(301)
     if (twofiles) close(302)
  else
     random_anconcolpos=.true.
     RXCOLGEO(:,:,:)=0
  endif
  
end subroutine initial_conf

! ========================================================================== !
! ========================================================================== !
subroutine insert_2col(lbox,ncol,maxncol,nchainspercol,nblobsperchain,maxnblob,rcol,&
     poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,celllistc,celllistb,ipcc,ipcb,&
     kfactor,ligboundto,rcut_blobblob,rcut_blobwall,rxcolgeo,random_anconcolpos,tot_ene)
  
  implicit none
  integer, intent(in) :: mnpic,ncellsc(3),maxncol,nchainspercol,nblobsperchain,&
       ncellsb(3),kfactor,maxnblob
  integer, intent(inout) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol),celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),&
       ipcb(4,maxnblob),ncol,ligboundto(2,nchainspercol,maxncol)
  real*8, intent(in) :: lbox(3),rcol,socc(3),socb(3),&
       rxcolgeo(3,nchainspercol,maxncol),rcut_blobblob,rcut_blobwall
  real*8, intent(inout) :: poscol(3,maxncol),posblob(3,nblobsperchain,&
       nchainspercol,maxncol),tot_ene
  logical, intent(in) :: random_anconcolpos
  real*8 :: Eold,Enew,rnd,rnd3(3),rinew(3),riold(3),arg,volume,ctb(3,nchainspercol),&
       S1rot,S2rot,V1rot(2),V2rot(2),ROTMATRIX(3,3),sqrtSrot,q(4),rnd2(2),&
       trialpos(3,kfactor),trialdist,xm,ym,opnemlre,effdist,Vmars(2),Smars,sqrtSmars,&
       dirvector(3),expmUk(kfactor),norexpmUk(kfactor),Wtot,commnorexpmUk,distij(3),&
       Ukext(kfactor),Uext(nblobsperchain,nchainspercol),Uk(kfactor),&
       W(nblobsperchain,nchainspercol),ext_col_ene
  integer :: kk,parti, inewcell(3),jj,ii,ifreecol,iifreecol,ibond,ipcold(3),&
       ipcnew(3),nfreecol,iblob,iiblob,jblob,jjblob,iichain,jjchain,itrial,iicol,nrec
  logical :: overlap
  real*8, parameter :: pi=3.14159265d0, O=0.0d0 

  nrec=0
  
  poscol(:,1) = (/ lbox(1)/2, lbox(2)/2, lbox(3)/2 /)
  poscol(:,2) = (/ 3.0*lbox(1)/4, lbox(2)/2, lbox(3)/2 /) 

  ! ANCHOr CHAINS TO THE PARTICLE
  if (random_anconcolpos) then ! random positions of chain anchors on colloid
     iicol=1
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
     iicol=2
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
     
  else  ! ANCHOR POSITIONS DETERMINED BY RXCOLGEO
     do iicol=1,2
        do iichain=1,nchainspercol
           posblob(:,1,iichain,iicol)= poscol(:,iicol) + rcol*rxcolgeo(:,iichain,iicol)
        enddo
     enddo
  endif ! random_anconcolpos

iicol=1           
! DO ROSENBLUTH SAMPLING FOR ALL CHAINS 1
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
! END ROSENBLUTH 1

! ACCEPT INSERTION         
tot_ene=tot_ene+sum(Uext(2:nblobsperchain,:))+ext_col_ene


iicol=2           
! DO ROSENBLUTH SAMPLING FOR ALL CHAINS 2
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
! END ROSENBLUTH 2


! ACCEPT INSERTION         
tot_ene=tot_ene+sum(Uext(2:nblobsperchain,:))+ext_col_ene



! UPDATE COL 1
iicol=1
ncol=1
inewcell(:) = floor(poscol(:,iicol)/socc(:)) + 1
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
   enddo
enddo

! UPDATE COL 2
iicol=2
ncol=2
inewcell(:) = floor(poscol(:,iicol)/socc(:)) + 1
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
   enddo
enddo

end subroutine insert_2col
!============================================================================!

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
! XXXXXXXXXXXXXXXXXX   MC   MOVE      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
subroutine mc_move2WL(lbox,ncol,maxncol,nchainspercol,nblobsperchain,maxnblob,rcol,&
     poscol,posblob,mnpic,ncellsc,ncellsb,socc,socb,socr,celllistc,&
     celllistb,ipcc,ipcb,ligboundto,rcut_blobblob,&
     rcut_blobwall,rcut_liglig,max_hop_blob,max_hop_col,max_rot_col,&
     movebacc,movecacc,tot_ene,RXBONDENE,nrxspec,ichainspec,mobile_ligands, &
     maxdelRcol2,nWL_points,max_colzWL,psiWL,fixcol_flag,nbonds)
  implicit none
  integer, intent(in) :: mnpic,ncellsc(3),ncol,maxncol,nchainspercol,nblobsperchain,&
       ncellsb(3),maxnblob,nrxspec,ichainspec(nchainspercol,maxncol),nWL_points
  integer, intent(inout) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),&
       ipcc(4,maxncol),celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),&
       ipcb(4,maxnblob),movebacc,movecacc,ligboundto(2,nchainspercol,maxncol),&
       nbonds
  real*8, intent(in) :: lbox(3),rcol,socc(3),socb(3),socr(2),rcut_blobblob,&
       rcut_blobwall,rcut_liglig,max_hop_col,max_rot_col,max_hop_blob,&
       RXBONDENE(nrxspec,nrxspec), maxdelRcol2,&
       psiWL(nWL_points),max_colzWL
  logical, intent(in) :: mobile_ligands,fixcol_flag
  real*8, intent(inout) :: poscol(3,maxncol),posblob(3,nblobsperchain,&
       nchainspercol,maxncol),tot_ene
  real*8 :: Eold,Enew,rnd,rnd3(3),rinew(3),riold(3),arg,volume,ctb(3,nchainspercol),&
       Srot,Vrot(2),ROTMATRIX(3,3),sqrtSrot,urot(3),sinfi,cosfi,omcosfi,rnd2(2),&
       distij(3),Qold,Qnew,rancnew(3,nchainspercol),garg,q1bi_old(mnpic),q1bi_new(mnpic),&
       rcut_liglig2,distij2,distijold(3),distij2old,tmpreal8,oldbondelene,newbondelene
  integer :: kk,parti, inewcell(3),jj,ii,ifreecol,iifreecol,ibond,ipcold(3),&
       ipcnew(3),nfreecol,iblob,iiblob,jblob,jjblob,iichain,jjchain,itrial,iicol,&
       jrecind,jjrec,nlighomies_new,nlighomies_old,whichhomies_old(2,mnpic*10),&
       whichhomies_new(2,mnpic*10),ihomie,nrec,jjcol,WLoldbin,WLnewbin
  logical :: overlap,lastbond
  real*8, parameter :: pi=3.14159265d0, O=0.0d0 

  nrec=0
  if (ncol .lt. 1) return
  rcut_liglig2=rcut_liglig*rcut_liglig
  WLoldbin=1
  WLnewbin=1
  
  garg=5.d0/(nblobsperchain*nchainspercol + 1)
  call random_number(rnd)  
  if (rnd .lt. garg) then ! move colloid
    
     call random_number(rnd) ! get random colloid
     iicol=floor(rnd*ncol)+1

     if (iicol .eq. 1) then
        rinew = poscol(:,iicol)
     else
        if (iicol .ne. 2) WRITE(*,*) 'ERROR: iicol must be 2, something wrong in WLmc_move'
        rinew(:)=poscol(:,iicol)

        if (.not. fixcol_flag) then !  move the colloid
           ! get old distance
           distij(1:3) = rinew(1:3) - poscol(1:3,1)
           distij2old=sum(distij(1:3)*distij(1:3))
           
           call random_number(rnd3)  ! get new position  
           rinew(:)=rinew(:)+max_hop_col*2*(rnd3-0.5d0)
           distij(1:3) = rinew(1:3) - poscol(1:3,1)   ! no periodic boundary!! anyway box is large and col 1 is always at the centre. 
           distij2=sum(distij(1:3)*distij(1:3))
           !   rinew(1:3)=rinew(1:3)-lbox(1:3)*floor(rinew(1:3)/lbox(1:3))
           ! check if col-col distance is not too large
           !  distij2 = sum((rinew(:) - poscol(:,1))**2)
          
           if ((distij2 .gt. max_colzWL*max_colzWL ) .and. (distij2 .gt. distij2old) .and. &
                (sum(ligboundto(:,:,iicol)) .eq.0))  rinew(:) = poscol(:,iicol) ! move back to the starting position 
           
        endif
        ! get WL bins
     !   WLnewbin=floor((rinew(1)-poscol(1,1)-rWL_min )/WLbinsize) + 1
     !   WLoldbin=floor((poscol(1,2)- poscol(1,1)-rWL_min )/WLbinsize) + 1
     endif
     inewcell(:) = floor(rinew(:)/socc(:))+1
     !CHECK FOR OVERLAP WITH OTHER COLLOIDS
     call col_overlap(overlap,iicol,rinew,inewcell,mnpic,&
          celllistc,ipcc,ncellsc,lbox,poscol,rcol,maxncol)
     ! overlap with wall
    ! if ((rinew(3).lt. rcol).or.(rinew(3).gt.lbox(3)-rcol)) overlap =.true.     
     
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
       
        if (sum(ligboundto(:,:,iicol)) .eq.0) then ! no bonds present, rotate by a large angle
           sinfi=1*(2*rnd3(3)-1.0) ! max rotate by  pi/2
        else ! rotate by a small angle
           sinfi=max_rot_col*(2*rnd3(3)-1.0)
        endif
        
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
        if (nblobsperchain .gt. 1) then
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
        else !take into account bond stretching energy when moving colloid

           ! AAAAAAAA
           ! loop over all formed bonds
       
           do iichain=1,nchainspercol ! loop over the ii colloid
              jjcol=ligboundto(1,iichain,iicol)
              jjchain=ligboundto(2,iichain,iicol)
              if (jjcol .gt. 0) then  ! bond is present
                 distij=posblob(:,1,iichain,iicol) - posblob(:,1,jjchain,jjcol)
                 distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))
                 distij2=sum(distij*distij)
                 call fbondene(tmpreal8,distij2)
                 Eold=Eold + tmpreal8 

                 ! new distance
                 distij=rancnew(:,iichain) - posblob(:,1,jjchain,jjcol)
                 distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))
                 distij2=sum(distij*distij)
                 call fbondene(tmpreal8,distij2)
                 Enew=Enew + tmpreal8
              endif
           enddo
                      
        endif


        
        arg=exp(Eold-Enew) ! *exp(psiWL(WLnewbin)-psiWL(WLoldbin))
        call random_number(rnd)
        if (rnd .lt. arg) then ! ACCEPT MOVE
           movecacc=movecacc+1
           posblob(:,1,:,iicol)=rancnew(:,:)
           poscol(:,iicol)=rinew
           tot_ene=tot_ene+Enew-Eold
           
           call make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
                poscol,posblob,lbox,&
                celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,socc,socb,mnpic)
                     
!!$           if (any(inewcell .ne. ipcc(1:3,iicol))) then
!!$
!!$              call update_cell_list(iicol,maxncol,inewcell,ipcc,celllistc,&
!!$                   ncellsc,mnpic)
!!$           endif
        else ! REJECT   
        endif
     endif ! OVERLAP
     
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!     
!xxxxxxxxxxxxxxxxxxxxxxxxxx  MOVE BLOBS  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx! 
  else   ! MOVE BLOBS
     call random_number(rnd3)
     iicol=floor(rnd3(1)*ncol)+1
     iichain=floor(rnd3(2)*nchainspercol)+1
     iiblob=floor(rnd3(3)*(nblobsperchain))+1   ! can also move anchors

     if (nblobsperchain .eq. 1) then ! ligand is the anchor  12 feb 2018
        rinew(:)= posblob(:,iiblob,iichain,iicol) ! new pos is the old pos
        if (mobile_ligands) then
           call random_number(rnd3)
           rinew(:)=posblob(:,iiblob,iichain,iicol)+max_hop_blob*2*(rnd3(:)-0.5d0)
           ctb(:,1)=rinew(:)-poscol(:,iicol) ! colloid to blob distance
           ctb(:,1)=ctb(:,1)/sqrt(sum(ctb(:,1)*ctb(:,1)))*rcol ! normalize to colloid surface
           rinew(:)=poscol(:,iicol)+ctb(:,1) ! new position on colloid surface
        endif
        
        Eold=0
        Enew=0
        
        ! DO BINDING WITH OTHER LIGANDS
        lastbond=.false. ! flag we are moving the last bound ligand
        Qold=1
        Qnew=1
        nlighomies_old=0
        nlighomies_new=0
        
        
        ! check if colloid is outside of distance
        distij(1:3) = poscol(1:3,2) - poscol(1:3,1)
        distij2=sum(distij(:)*distij(:))
        if ( distij2 .gt.  max_colzWL * max_colzWL) then
           ! check if  this blob has the last bond
           if((nbonds .eq. 1).and.(ligboundto(1,iichain,iicol).gt.0 )) then ! this ligand has the last bond
              lastbond=.true.
              Qold=0 ! unbound state is not possible
              Qnew=0
           elseif(nbonds .eq. 0) then
              write(*,*) 'WARNING: nbonds = 0, should not be the case for Wl bonds.. '
              write(*,*) 'only ok if it heapens at the begining and the initial colloid height is > h_0'
              return  
           endif
           
        endif
        
        ! can do binding        
        ! first unbind
        
        if (ligboundto(1,iichain,iicol).gt.0) then
           distij(:)=posblob(:,nblobsperchain,iichain,iicol) - posblob(:, &
                nblobsperchain,ligboundto(2,iichain,iicol),ligboundto(1,iichain,iicol))
           distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
           distij2=sum(distij*distij)
           if (distij2 .gt. rcut_liglig2 ) then
              write(*,*) 'WARNING: trying to delete a bond but distij > rcut_liglig' ! cannot unbind if distij2 > rcut2, no reverse move in configurational bias !!!!
              return
           endif
           call fbondene(tmpreal8,distij2)
           tot_ene=tot_ene-tmpreal8  !  delete the bond  
           
           ligboundto(:,ligboundto(2,iichain,iicol),ligboundto(1,iichain,iicol))=0
           ligboundto(:,iichain,iicol)=0
           nbonds=nbonds-1
        endif
        
        call find_lig_homie(posblob(:,nblobsperchain,iichain,iicol),nchainspercol,nblobsperchain,lbox,posblob,&
             celllistb,ipcb,ncellsb,socb,mnpic,rcut_liglig,ligboundto,maxncol,maxnblob,&
             nlighomies_old,whichhomies_old)
        
        
        ! get the old binding partition function
        q1bi_old(:)=0  
        if ((nlighomies_old .gt. 0) .and. (nbonds < nWL_points - 1 )) then
           
           ! get partition funcition with psi bias
           do ihomie=1,nlighomies_old
              if (ichainspec(iichain,iicol) .eq. ichainspec(whichhomies_old(2,ihomie), &
                   whichhomies_old(1,ihomie))) then
                 ! homie is of equal type, must exclude otherwise ligands can bind to themselves!!
                 q1bi_old(ihomie)=0
              else
                 distij(:)=posblob(:,nblobsperchain,iichain,iicol) - posblob(:,nblobsperchain,&
                      whichhomies_old(2,ihomie),whichhomies_old(1,ihomie))
                 distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                 
                 call fbondene(oldbondelene,sum(distij*distij))  ! old bond elastic energy
                 q1bi_old(ihomie)=dexp(-RXBONDENE(ichainspec(whichhomies_old(2,ihomie),whichhomies_old(1,ihomie)),&
                      ichainspec(iichain,iicol)) - oldbondelene &
                      + psiWL(nbonds+2)-psiWL(nbonds+1))
                 !  write(*,*)q1bi_old(ihomie)
              endif
           enddo
           
           ! q1bi_old(:)=q1bi_old(:)*dexp(psiWL(nbonds+2)-psiWL(nbonds+1))  !WL
           !     write(*,*) 'AAAA ', q1bi_old(:nlighomies_old)
           !     write(*,*) 'BBBB ', dexp(700.0d0),dexp(-700.0d0),dexp(RXBONDENE(1,1)),dexp(-RXBONDENE(1,1))
           Qold=Qold+sum(q1bi_old(1:nlighomies_old))
           !  write(*,*) 'AAAQo ', nbonds, Qold, q1bi_old(1:nlighomies_old)
           !  write(*,*) psiWL
        endif ! nlighomies > 0
        

        if (mobile_ligands) then  ! new partition function is different
           ! new position homies
           call find_lig_homie(rinew,nchainspercol,nblobsperchain,lbox,posblob,&
                celllistb,ipcb,ncellsb,socb,mnpic,rcut_liglig,ligboundto,maxncol,maxnblob,&
                nlighomies_new,whichhomies_new)

           ! get the new binding partition function
           q1bi_new(:)=0
           if ((nlighomies_new .gt. 0) .and. (nbonds < nWL_points - 1)) then

              ! get partition funcition
              do ihomie=1,nlighomies_new
                 if (ichainspec(iichain,iicol) .eq. ichainspec(whichhomies_new(2,ihomie), &
                      whichhomies_new(1,ihomie))) then ! homie is of equal type
                    q1bi_new(ihomie)=0
                 else

                    distij(:)=rinew(:)-posblob(:,nblobsperchain,whichhomies_new(2,ihomie),whichhomies_new(1,ihomie))
                    distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                    call fbondene(newbondelene,sum(distij*distij))

                    q1bi_new(ihomie)=dexp(-RXBONDENE(ichainspec(whichhomies_new(2,ihomie),whichhomies_new(1,ihomie)),&
                         ichainspec(iichain,iicol)) -newbondelene &
                         +psiWL(nbonds+2)-psiWL(nbonds+1))
                 endif
              enddo
              !      q1bi_new(:)=q1bi_new(:)*exp(psiWL(nbonds+2)-psiWL(nbonds+1)) !WL
              Qnew=Qnew+sum(q1bi_new(1:nlighomies_new))
              !   write(*,*) 'BBBQn ', nbonds, Qnew, q1bi_new(1:nlighomies_new)
           endif ! nrechomies > 0

        else ! new partition function is the same as the old one
           Qnew=Qold
           q1bi_new=q1bi_old
           nlighomies_new=nlighomies_old
           whichhomies_new=whichhomies_old         
        endif ! mobile_ligands
        
     endif ! iiblob .eq. nblobsperchain
     ! END BiNDING WITH LIGANDS
     
     ! DO METROPOLIS
     arg=exp(Eold-Enew)*Qnew/Qold
     !    if (isnan(arg)) then
     !       write(*,*) Qnew, Qold,(0.5 < arg)
     !    endif
     call random_number(rnd)
     if (rnd .lt. arg) then! ACCEPT
        
        movebacc=movebacc+1
        posblob(:,iiblob,iichain,iicol)=rinew(:)
        tot_ene=tot_ene+Enew-Eold        
        
        inewcell(:) = floor(rinew(:)/socb(:))+1
        !!         inewcell(3)=min(inewcell(3),ncellsb(3)) ! because soft blobs can in principle penetrate the wall slightly
        !!         inewcell(3)=max(inewcell(3),1)
        iblob=(iicol-1)*nblobsperchain*nchainspercol+(iichain-1)*nblobsperchain+iiblob
        if (mobile_ligands) then
           if (any(inewcell .ne. ipcb(1:3,iblob))) then
              call update_cell_list(iblob,maxnblob,inewcell,ipcb,celllistb,&
                   ncellsb,mnpic)
           endif
        endif
        
        ! bind to some ligand
        if (lastbond) then ! unbound state is not possible, last bond and colloid is outside of range.
           arg=0
        else
           arg=1.0d0/Qnew
        endif
        
        !  write(*,*) 'CCCC ', arg
        !  write(*,*)
        if (nbonds < nWL_points - 1) then ! less than maximum number of bonds, can bind more
           q1bi_new(:)=q1bi_new(:)/Qnew ! normalise partition functions
           call random_number(rnd)
           ! bind to some receptor
           if (rnd .gt. arg) then !  bind
              
              do ihomie=1,nlighomies_new
                 arg=arg+q1bi_new(ihomie)
                 if ((rnd .le. arg) .or. (isnan(arg))) then
                    jjcol=whichhomies_new(1,ihomie)
                    jjchain=whichhomies_new(2,ihomie)
                    ligboundto(1:2,iichain,iicol)=(/jjcol,jjchain/)
                    ligboundto(1:2,jjchain,jjcol)=(/iicol,iichain/)
                    nbonds=nbonds+1  ! add bond to global bonds count
                    
                    ! update global energy
                    distij(:)=rinew(:)-posblob(:,nblobsperchain,jjchain,jjcol)
                    distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                    call fbondene(newbondelene,sum(distij*distij))
                    tot_ene=tot_ene+newbondelene
                    
                    !        write(*,*) iicol, iichain
                    !        write(*,*) jjcol, jjchain
                    exit
                 endif
              enddo
           endif
        endif
        
     else ! REJECT
        ! rebind to some ligand
        if (lastbond) then
           arg=0
        else
           arg=1.0d0/Qold
        endif
        !  write(*,*) 'DDDD ', arg
        !  write(*,*)
        if (nbonds < nWL_points - 1) then ! less than maximum number of bonds, can bind more
           q1bi_old(:)=q1bi_old(:)/Qold ! normalise partition functions
           call random_number(rnd)
           ! bind to some receptor
           if (rnd .gt. arg) then !  bind
              
              do ihomie=1,nlighomies_old
                 arg=arg+q1bi_old(ihomie)
                 if ((rnd .le. arg) .or. (isnan(arg))) then
                    jjcol=whichhomies_old(1,ihomie)
                    jjchain=whichhomies_old(2,ihomie)
                    ligboundto(1:2,iichain,iicol)=(/jjcol,jjchain/)
                    ligboundto(1:2,jjchain,jjcol)=(/iicol,iichain/)
                    nbonds=nbonds+1

                    ! update global energy
                    distij(:)=posblob(:,nblobsperchain,iichain,iicol ) - posblob(:,nblobsperchain,jjchain,jjcol)
                    distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
                    call fbondene(newbondelene,sum(distij*distij))
                    tot_ene=tot_ene+newbondelene
                    
                    exit
                 endif
              enddo
           endif ! bind
        endif
     endif ! accept/reject
  

!!$     if (mobile_ligands .and. (iiblob .eq. 1)) then ! move ligand anchor
!!$        call random_number(rnd3)
!!$        rinew(:)=posblob(:,iiblob,iichain,iicol)+max_hop_blob*2*(rnd3(:)-0.5d0)
!!$        ctb(:,1)=rinew(:)-poscol(:,iicol) ! colloid to blob distance
!!$        ctb(:,1)=ctb(:,1)/sqrt(sum(ctb(:,1)*ctb(:,1)))*rcol ! normalize to colloid surface
!!$        rinew(:)=poscol(:,iicol)+ctb(:,1) ! new position on colloid surface
!!$        
!!$        distij(:)=posblob(:,2,iichain,iicol)-posblob(:,1,iichain,iicol)
!!$        distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))
!!$        Eold=0.75*sum(distij*distij) 
!!$        distij(:)=posblob(:,2,iichain,iicol)-rinew(:)
!!$        distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))
!!$        Enew=0.75*sum(distij*distij)
!!$
!!$   !     write(*,*)
!!$   !     write(*,*) posblob(:,iiblob,iichain,iicol)-rinew(:) 
!!$   !     write(*,*) rcol, rinew
!!$   !     write(*,*) Enew, Eold
!!$        ! DO METROPOLIS
!!$        arg=exp(Eold-Enew)
!!$        call random_number(rnd)
!!$        if (rnd .lt. arg) then! ACCEPT
!!$           movebacc=movebacc+1
!!$           posblob(:,1,iichain,iicol)=rinew(:)
!!$           tot_ene=tot_ene+Enew-Eold 
!!$        endif
!!$
!!$     elseif (iiblob .gt. 1) then ! move blobs   
!!$        ! get new blob position
!!$        call random_number(rnd3)
!!$        rinew(:)=posblob(:,iiblob,iichain,iicol)+max_hop_blob*2*(rnd3(:)-0.5d0)
!!$        rinew(1:3)=rinew(1:3)-lbox(1:3)*floor(rinew(1:3)/lbox(1:3))
!!$        
!!$        Eold=0
!!$        Enew=0
!!$        call blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
!!$             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
!!$             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
!!$             posblob(:,iiblob,iichain,iicol),iicol,iichain,iiblob,Eold)
!!$        call blob_energy_calc(posblob,poscol,lbox,nrec,maxncol,nchainspercol,& 
!!$             nblobsperchain,maxnblob,celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,&
!!$             socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,&
!!$             rinew,iicol,iichain,iiblob,Enew)
!!$        
!!$        ! DO BINDING WITH OTHER LIGANDS
!!$        lastbond=.false. ! flag we are moving the last bound ligand
!!$        Qold=1
!!$        Qnew=1
!!$        nlighomies_old=0
!!$        nlighomies_new=0
!!$        if (iiblob .eq. nblobsperchain) then
!!$
!!$           ! check if colloid is outside of distance
!!$           distij(1:3) = poscol(1:3,2) - poscol(1:3,1)
!!$           distij2=sum(distij(:)*distij(:))
!!$           if ( distij2 .gt.  max_colzWL * max_colzWL) then
!!$              ! check if  this blob has the last bond
!!$              if((nbonds .eq. 1).and.(ligboundto(1,iichain,iicol).gt.0 )) then ! this ligand has the last bond
!!$                 lastbond=.true.
!!$                 Qold=0 ! unbound state is not possible
!!$                 Qnew=0
!!$              elseif(nbonds .eq. 0) then
!!$                 write(*,*) 'WARNING: nbonds = 0, should not be the case for Wl bonds.. '
!!$                 write(*,*) 'only ok if it heapens at the begining and the initial colloid height is > h_0'
!!$                 return  
!!$              endif
!!$             
!!$           endif
!!$   
!!$           ! can do binding        
!!$           ! first unbind
!!$           
!!$           if (ligboundto(1,iichain,iicol).gt.0) then
!!$              distij(:)=posblob(:,nblobsperchain,iichain,iicol) - posblob(:, &
!!$                   nblobsperchain,ligboundto(2,iichain,iicol),ligboundto(1,iichain,iicol))
!!$              distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
!!$              if (sum(distij*distij) .gt. rcut_liglig2 ) then
!!$                 write(*,*) 'WARNING: trying to delete a bond but distij > rcut_liglig' ! cannot unbind if distij2 > rcut2, no reverse move in configurational bias !!!!
!!$              endif
!!$              
!!$              ligboundto(:,ligboundto(2,iichain,iicol),ligboundto(1,iichain,iicol))=0
!!$              ligboundto(:,iichain,iicol)=0
!!$              nbonds=nbonds-1
!!$           endif
!!$           
!!$           call find_lig_homie(posblob(:,nblobsperchain,iichain,iicol),nchainspercol,nblobsperchain,lbox,posblob,&
!!$                celllistb,ipcb,ncellsb,socb,mnpic,rcut_liglig,ligboundto,maxncol,maxnblob,&
!!$                nlighomies_old,whichhomies_old)
!!$
!!$           
!!$           ! get the old binding partition function
!!$           q1bi_old(:)=0  
!!$           if ((nlighomies_old .gt. 0) .and. (nbonds < nWL_points - 1 )) then
!!$             
!!$              ! get partition funcition with psi bias
!!$              do ihomie=1,nlighomies_old
!!$                 if (ichainspec(iichain,iicol) .eq. ichainspec(whichhomies_old(2,ihomie), &
!!$                      whichhomies_old(1,ihomie))) then
!!$                    ! homie is of equal type, must exclude otherwise ligands can bind to themselves!!
!!$                    q1bi_old(ihomie)=0
!!$                 else
!!$                    distij(:)=posblob(:,nblobsperchain,iichain,iicol) - posblob(:,nblobsperchain,&
!!$                         whichhomies_old(2,ihomie),whichhomies_old(1,ihomie))
!!$                    distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
!!$                    q1bi_old(ihomie)=dexp(-RXBONDENE(ichainspec(whichhomies_old(2,ihomie),whichhomies_old(1,ihomie)),&
!!$                         ichainspec(iichain,iicol)) - 0.534*(sqrt(sum(distij*distij))-0.73)**2 & 
!!$                         + psiWL(nbonds+2)-psiWL(nbonds+1))
!!$                    !  write(*,*)q1bi_old(ihomie)
!!$                 endif
!!$              enddo
!!$             
!!$             ! q1bi_old(:)=q1bi_old(:)*dexp(psiWL(nbonds+2)-psiWL(nbonds+1))  !WL
!!$         !     write(*,*) 'AAAA ', q1bi_old(:nlighomies_old)
!!$         !     write(*,*) 'BBBB ', dexp(700.0d0),dexp(-700.0d0),dexp(RXBONDENE(1,1)),dexp(-RXBONDENE(1,1))
!!$              Qold=Qold+sum(q1bi_old(1:nlighomies_old))
!!$            !  write(*,*) 'AAAQo ', nbonds, Qold, q1bi_old(1:nlighomies_old)
!!$            !  write(*,*) psiWL
!!$           endif ! nlighomies > 0
!!$           
!!$              
!!$           ! new position homies
!!$           call find_lig_homie(rinew,nchainspercol,nblobsperchain,lbox,posblob,&
!!$                celllistb,ipcb,ncellsb,socb,mnpic,rcut_liglig,ligboundto,maxncol,maxnblob,&
!!$                nlighomies_new,whichhomies_new)
!!$           
!!$           ! get the new binding partition function
!!$           q1bi_new(:)=0
!!$           if ((nlighomies_new .gt. 0) .and. (nbonds < nWL_points - 1)) then
!!$              
!!$              ! get partition funcition
!!$              do ihomie=1,nlighomies_new
!!$                 if (ichainspec(iichain,iicol) .eq. ichainspec(whichhomies_new(2,ihomie), &
!!$                      whichhomies_new(1,ihomie))) then ! homie is of equal type
!!$                    q1bi_new(ihomie)=0
!!$                 else
!!$                    
!!$                    distij(:)=rinew(:)-posblob(:,nblobsperchain,whichhomies_new(2,ihomie),whichhomies_new(1,ihomie))
!!$                    distij(:)=distij(:)-lbox(:)*nint(distij(:)/lbox(:))
!!$                    q1bi_new(ihomie)=dexp(-RXBONDENE(ichainspec(whichhomies_new(2,ihomie),whichhomies_new(1,ihomie)),&
!!$                         ichainspec(iichain,iicol)) - 0.534*(sqrt(sum(distij*distij))-0.73)**2 &
!!$                         +psiWL(nbonds+2)-psiWL(nbonds+1))
!!$                 endif
!!$              enddo
!!$        !      q1bi_new(:)=q1bi_new(:)*exp(psiWL(nbonds+2)-psiWL(nbonds+1)) !WL
!!$              Qnew=Qnew+sum(q1bi_new(1:nlighomies_new))
!!$           !   write(*,*) 'BBBQn ', nbonds, Qnew, q1bi_new(1:nlighomies_new)
!!$           endif ! nrechomies > 0
!!$           
!!$        endif ! iiblob .eq. nblobsperchain
!!$        ! END BiNDING WITH LIGANDS
!!$        
!!$        ! DO METROPOLIS
!!$        arg=exp(Eold-Enew)*Qnew/Qold
!!$    !    if (isnan(arg)) then
!!$    !       write(*,*) Qnew, Qold,(0.5 < arg)
!!$    !    endif
!!$        call random_number(rnd)
!!$        if (rnd .lt. arg) then! ACCEPT
!!$           
!!$           movebacc=movebacc+1
!!$           posblob(:,iiblob,iichain,iicol)=rinew(:)
!!$           tot_ene=tot_ene+Enew-Eold        
!!$           
!!$           inewcell(:) = floor(rinew(:)/socb(:))+1
!!$  !!         inewcell(3)=min(inewcell(3),ncellsb(3)) ! because soft blobs can in principle penetrate the wall slightly
!!$  !!         inewcell(3)=max(inewcell(3),1)
!!$           iblob=(iicol-1)*nblobsperchain*nchainspercol+(iichain-1)*nblobsperchain+iiblob
!!$           if (any(inewcell .ne. ipcb(1:3,iblob))) then
!!$              call update_cell_list(iblob,maxnblob,inewcell,ipcb,celllistb,&
!!$                   ncellsb,mnpic)
!!$           endif
!!$           
!!$           ! bind to some ligand
!!$            if (lastbond) then ! unbound state is not possible, last bond and colloid is outside of range.
!!$              arg=0
!!$           else
!!$              arg=1.0d0/Qnew
!!$           endif
!!$           
!!$         !  write(*,*) 'CCCC ', arg
!!$         !  write(*,*)
!!$           if (nbonds < nWL_points - 1) then ! less than maximum number of bonds, can bind more
!!$              q1bi_new(:)=q1bi_new(:)/Qnew ! normalise partition functions
!!$              call random_number(rnd)
!!$              ! bind to some receptor
!!$              if (rnd .gt. arg) then !  bind
!!$                 
!!$                 do ihomie=1,nlighomies_new
!!$                    arg=arg+q1bi_new(ihomie)
!!$                    if ((rnd .le. arg) .or. (isnan(arg))) then
!!$                       jjcol=whichhomies_new(1,ihomie)
!!$                       jjchain=whichhomies_new(2,ihomie)
!!$                       ligboundto(1:2,iichain,iicol)=(/jjcol,jjchain/)
!!$                       ligboundto(1:2,jjchain,jjcol)=(/iicol,iichain/)
!!$                       nbonds=nbonds+1  ! add bond to global bonds count
!!$                       !        write(*,*) iicol, iichain
!!$                       !        write(*,*) jjcol, jjchain
!!$                       exit
!!$                    endif
!!$                 enddo
!!$              endif
!!$           endif
!!$           
!!$        else ! REJECT
!!$           ! rebind to some ligand
!!$           if (lastbond) then
!!$              arg=0
!!$           else
!!$              arg=1.0d0/Qold
!!$           endif
!!$         !  write(*,*) 'DDDD ', arg
!!$         !  write(*,*)
!!$           if (nbonds < nWL_points - 1) then ! less than maximum number of bonds, can bind more
!!$              q1bi_old(:)=q1bi_old(:)/Qold ! normalise partition functions
!!$              call random_number(rnd)
!!$              ! bind to some receptor
!!$              if (rnd .gt. arg) then !  bind
!!$                 
!!$                 do ihomie=1,nlighomies_old
!!$                    arg=arg+q1bi_old(ihomie)
!!$                    if ((rnd .le. arg) .or. (isnan(arg))) then
!!$                       jjcol=whichhomies_old(1,ihomie)
!!$                       jjchain=whichhomies_old(2,ihomie)
!!$                       ligboundto(1:2,iichain,iicol)=(/jjcol,jjchain/)
!!$                       ligboundto(1:2,jjchain,jjcol)=(/iicol,iichain/)
!!$                       nbonds=nbonds+1
!!$                       
!!$                       exit
!!$                    endif
!!$                 enddo
!!$              endif ! bind
!!$           endif
!!$        endif ! accept/reject
!!$     endif ! move anchor/other blob




     
  endif ! move col/blob
!  write(*,*) 'FFFFF ', nbonds
!  write(*,*) 'ligbonds1 ', ligboundto(1,:,1)
!  write(*,*) 'ligbonds2 ', ligboundto(1,:,2)
end subroutine mc_move2WL
!=============================================================!

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
!!$subroutine bond_create_destroy(posrec,posblob,lbox,nrec,ncol,maxncol,nchainspercol,& 
!!$     nblobsperchain,celllistr,ipcr,ncellsr,socr,mnpic,rcut_liglig,bindacc,&
!!$     ligboundto,recboundto,nrxspec,RXBONDENE,irecspec,ichainspec)
!!$  implicit none
!!$  integer, intent(in) :: maxncol,ncol,nrec,nblobsperchain,nchainspercol,ncellsr(2),mnpic,&
!!$       celllistr(mnpic,ncellsr(1),ncellsr(2)),ipcr(3,nrec),irecspec(nrec),&
!!$       ichainspec(nchainspercol),nrxspec
!!$  integer, intent(inout) :: bindacc,ligboundto(nchainspercol,maxncol),recboundto(2,nrec)
!!$  real*8, intent(in) :: RXBONDENE(nrxspec,nrxspec),posrec(2,nrec),lbox(3),socr(2),rcut_liglig,posblob(3,&
!!$       nblobsperchain,nchainspercol,maxncol)
!!$  real*8 :: rnd,rnd2(2),rnd3(3),arg,q1b,distij(3),distij2,rcut_liglig2,q1bi(mnpic)
!!$  integer :: kk,itlr,jj,ii,iicol,iichain,jjrec,nrechomies,whichhomies(mnpic),ihomie,jrecind
!!$  real*8, parameter :: pi=3.14159265d0 
!!$  
!!$  rcut_liglig2=rcut_liglig*rcut_liglig
!!$  if (ncol .lt. 1) return
!!$
!!$  ! pick a random ligand
!!$  call random_number(rnd2)
!!$  iicol=floor(ncol*rnd2(1))+1 
!!$  iichain=floor(nchainspercol*rnd2(2))+1 
!!$
!!$  if (posblob(3,nblobsperchain,iichain,iicol) .lt. rcut_liglig) then
!!$  !   write(*,*)
!!$  !   write(*,*) 'binding iicol iichain ',iicol,iichain , ncol
!!$  !   write(*,*) posblob(3,nblobsperchain,iichain,iicol)
!!$     ! can do binding
!!$  
!!$     ! first unbind
!!$     if (ligboundto(iichain,iicol).gt.0) then
!!$        distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,ligboundto(iichain,iicol))
!!$        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
!!$        distij(3)=posblob(3,nblobsperchain,iichain,iicol)
!!$        if (sum(distij*distij) .gt. rcut_liglig2 ) then
!!$           write(*,*) 'WARNING: bond above rcut_liglig, in bond_create_destr'
!!$           return ! cannot unbind if distij2 > rcut2, no reverse move in configurational bias !!!!
!!$        endif
!!$        recboundto(1:2,ligboundto(iichain,iicol))=0
!!$        ligboundto(iichain,iicol)=0
!!$     endif
!!$     
!!$     call find_rec_homie(posblob(:,nblobsperchain,iichain,iicol),nchainspercol,nblobsperchain,lbox,posrec,&
!!$          nrec,celllistr,ipcr,ncellsr,socr,mnpic,rcut_liglig,recboundto,&
!!$          nrechomies,whichhomies)
!!$     
!!$     
!!$     if (nrechomies .gt. 0) then
!!$        bindacc=bindacc+1
!!$        q1bi(:)=0
!!$        ! get partition funcition
!!$        do ihomie=1,nrechomies
!!$           distij(1:2)=posblob(1:2,nblobsperchain,iichain,iicol)-posrec(:,whichhomies(ihomie))
!!$           distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
!!$           distij(3)=posblob(3,nblobsperchain,iichain,iicol)
!!$           q1bi(ihomie)=exp(-RXBONDENE(irecspec(whichhomies(ihomie)),&
!!$                         ichainspec(iichain))-0.75*sum(distij*distij))
!!$        enddo
!!$       
!!$        q1b=1+sum(q1bi(1:nrechomies))
!!$        arg=1.0d0/q1b
!!$        q1bi(:)=q1bi(:)/q1b ! normalise partition functions
!!$        call random_number(rnd)
!!$        ! bind to some receptor
!!$        if (rnd .gt. arg) then !  bind
!!$           do ihomie=1,nrechomies
!!$              arg=arg+q1bi(ihomie)
!!$              if (rnd .le. arg) then
!!$                 jjrec=whichhomies(ihomie)
!!$                 ligboundto(iichain,iicol)=jjrec
!!$                 recboundto(1:2,jjrec)=(/iichain,iicol/)
!!$                 exit
!!$              endif
!!$           enddo
!!$        endif
!!$     endif ! nrechomies > 0
!!$  endif  
!!$end subroutine bond_create_destroy
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
              socc,socb,mnpic,rcut_blobblob,rcut_blobwall,rcol,ligboundto,tot_ene)
  implicit none
  integer, intent(in) :: maxncol,nrec,nblobsperchain,ncellsc(3),ncellsb(3),mnpic,&
       celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol),&
       nchainspercol,celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),maxnblob,&
       ipcb(4,maxnblob),ncol,ligboundto(2,nchainspercol,maxncol)
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
       rcut_blobblob2,time1,time2, iforcexyz(3),eCS_blobwall,eCS_blobblob,tmpreal8
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


  !add bond energy
  iicol=1
  do iichain=1,nchainspercol ! loop over the ii colloid
     jjcol=ligboundto(1,iichain,iicol)
     jjchain=ligboundto(2,iichain,iicol)
     if (jjcol .gt. 0) then  ! bond is present        
        ! new distance
        distij=posblob(:,1,iichain,iicol) - posblob(:,1,jjchain,jjcol)
        distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))
        distij2=sum(distij*distij)
        call fbondene(tmpreal8,distij2)
        Enew=Enew + tmpreal8
     endif
  enddo
  tot_ene=tot_ene+Enew
  

  
  
end subroutine tot_ene_calc
! ========================================================================== !
subroutine find_lig_homie(trialpos,nchainspercol,nblobsperchain,lbox,posblob,&
                celllistb,ipcb,ncellsb,socb,mnpic,rcut_liglig,ligboundto,maxncol,maxnblob,&
                nlighomies,whichhomies)
  
  implicit none
  integer, intent(in) :: nblobsperchain,ncellsb(3),mnpic,nchainspercol,&
       celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob),&
       ligboundto(2,nchainspercol,maxncol),maxncol,maxnblob
  real*8, intent(in) :: posblob(3,nblobsperchain,nchainspercol,maxncol),lbox(3),socb(3),rcut_liglig,trialpos(3)
  integer, intent(out) :: nlighomies,whichhomies(2,mnpic*10)
  real*8 :: rnd,rnd2(2),rnd3(3),arg,volume,effdist,pot_ene
  integer :: ii,jj,kk,icell,jcell,kcell,icelltmp,jcelltmp,kcelltmp,jlig,ipcnew(3),&
       jblob,jjblob,jjcol,jjchain,ncellstocheck(3)
  real*8 :: distij(3),distij2,rbr2,rcut_liglig2
  
  real*8, parameter :: pi=3.14159265d0
  rcut_liglig2=rcut_liglig**2
  nlighomies=0  
  whichhomies(:,:)=0
  ncellstocheck=ceiling(rcut_liglig / socb(1:3))
  
  ipcnew=floor(trialpos(1:3)/socb(1:3))+1
  do icelltmp=ipcnew(1)-ncellstocheck(1),ipcnew(1)+ncellstocheck(1)
     do jcelltmp=ipcnew(2)-ncellstocheck(2),ipcnew(2)+ncellstocheck(2)
        do kcelltmp=ipcnew(3)-ncellstocheck(3),ipcnew(3)+ncellstocheck(3)
           
           icell=icelltmp-ncellsb(1)*floor(dble(icelltmp-1)/ncellsb(1)+1.0d-8)
           jcell=jcelltmp-ncellsb(2)*floor(dble(jcelltmp-1)/ncellsb(2)+1.0d-8)
           kcell=kcelltmp-ncellsb(3)*floor(dble(jcelltmp-1)/ncellsb(3)+1.0d-8)
           
           do jj=1,celllistb(1,icell,jcell,kcell)
              jblob=celllistb(jj+1,icell,jcell,kcell)

              jjcol=ceiling(jblob/dble(nblobsperchain*nchainspercol)-1.0d-8)              
              jjchain=ceiling((jblob-(jjcol-1)*nblobsperchain*nchainspercol)/dble(nblobsperchain)-1.0d-8)    
              jjblob=mod(jblob,nblobsperchain)
              if (jjblob .eq. 0) jjblob=nblobsperchain ! must set if  mod gives 0
              !     write(*,*) jjcol,jjchain,jjblob
              if (jjblob .eq. nblobsperchain) then ! blobs is a ligand !!!
                 
                 
                 if (ligboundto(1,jjchain,jjcol) .ne. 0) cycle  ! ligand is not free
                 distij(1:3)=trialpos(1:3)-posblob(1:3,nblobsperchain,jjchain,jjcol)
                 distij(1:3)=distij(1:3)-lbox(1:3)*nint(distij(1:3)/lbox(1:3))             
                 distij2=sum(distij*distij)
                 if (distij2 .lt. rcut_liglig2) then
                    nlighomies=nlighomies+1
                    whichhomies(:,nlighomies) = (/jjcol,jjchain/)
                 endif
              endif ! jjblob == nblobsperchain
           enddo
        enddo
     enddo
  enddo
end subroutine find_lig_homie
!!$ ! ========================================================================== !
SUBROUTINE seed_random_number(Xseed)
  implicit none
  ! Local variables
  integer*8, intent(in) :: Xseed
  INTEGER              :: ii,k, now(3)
  INTEGER, ALLOCATABLE :: seed(:)
  CALL RANDOM_SEED(SIZE=k)
  ALLOCATE( seed(k) )
  call itime(now)
  do ii=1, k
     seed(ii)=Xseed+4*ii
  end do
  CALL RANDOM_SEED(PUT=seed)
  DEALLOCATE( seed )
  RETURN
END SUBROUTINE seed_random_number

 ! ========================================================================== !
subroutine make_cell_list(ncol,maxncol,maxnblob,nchainspercol,nblobsperchain,&
     poscol,posblob,lbox,&
     celllistc,ipcc,celllistb,ipcb,ncellsc,ncellsb,socc,socb,mnpic)
 
  implicit none
  integer, intent(in) :: ncellsc(3),ncellsb(3),nchainspercol,&
       nblobsperchain,maxnblob,mnpic,ncol,maxncol
  real*8, intent(in) :: lbox(3),socc(3),socb(3),&
       posblob(3,nblobsperchain,nchainspercol,maxncol),poscol(3,maxncol)
  integer, intent(out) :: celllistc(mnpic,ncellsc(1),ncellsc(2),ncellsc(3)),ipcc(4,maxncol)
  integer, intent(out) :: celllistb(mnpic,ncellsb(1),ncellsb(2),ncellsb(3)),ipcb(4,maxnblob)
  integer ::  ii,jj,kk,nblob,iblob,icol

  celllistc(:,:,:,:)=0
  ipcc(:,:)=0
  celllistb(:,:,:,:)=0
  ipcb(:,:)=0
!write(*,*) 'adasdasdasdasd'
  
  nblob=ncol*nchainspercol*nblobsperchain
 
  do ii=1,ncol
     ipcc(1:3,ii)=floor(poscol(:,ii)/socc(:))+1        
     do jj=1,nchainspercol
        do kk=min(2,nblobsperchain),nblobsperchain ! ANCHOR IS IN THE CELL LIST IF NBLOBSPERCHAIN==1
           iblob=(ii-1)*nblobsperchain*nchainspercol+(nblobsperchain)*(jj-1)+kk
           ipcb(1:3,iblob)=floor(posblob(:,kk,jj,ii)/socb(:))+1
           
        !   if (ipcb(3,iblob) .lt. 1) ipcb(3,iblob)=1
           !   if (ipcb(3,iblob) .gt. ncellsb(3)) ipcb(3,iblob)=ncellsb(3)
           
           celllistb(1, ipcb(1,iblob), ipcb(2,iblob), ipcb(3,iblob)) = & 
                celllistb(1, ipcb(1,iblob), ipcb(2,iblob), ipcb(3,iblob)) + 1
           celllistb((celllistb(1,ipcb(1,iblob),ipcb(2,iblob),ipcb(3,iblob))+1),ipcb(1,iblob),ipcb(2,iblob),ipcb(3,iblob))=iblob    
           ipcb(4,iblob)=celllistb(1,ipcb(1,iblob), ipcb(2,iblob), ipcb(3,iblob))
           
        enddo
     enddo
  enddo
 
   do ii=1,ncol
     celllistc(1, ipcc(1,ii), ipcc(2,ii), ipcc(3,ii)) = & 
     celllistc(1, ipcc(1,ii), ipcc(2,ii), ipcc(3,ii)) + 1
     celllistc((celllistc(1,ipcc(1,ii),ipcc(2,ii),ipcc(3,ii))+1),ipcc(1,ii),ipcc(2,ii),ipcc(3,ii))=ii    
     ipcc(4,ii)=celllistc(1,ipcc(1,ii), ipcc(2,ii), ipcc(3,ii))
  enddo
!!$
!!$  do icol=1,ncol 
!!$     do jj=1,nchainspercol
!!$        do kk=2,nblobsperchain ! ANCHOR IS NOT IN THE CELL LIST
!!$           ii = (icol-1)*(nblobsperchain*nchainspercol) + (jj-1) * &
!!$                nblobsperchain + kk 
!!$           
!!$           celllistb(1, ipcb(1,ii), ipcb(2,ii), ipcb(3,ii)) = & 
!!$                celllistb(1, ipcb(1,ii), ipcb(2,ii), ipcb(3,ii)) + 1
!!$           celllistb((celllistb(1,ipcb(1,ii),ipcb(2,ii),ipcb(3,ii))+1),ipcb(1,ii),ipcb(2,ii),ipcb(3,ii))=ii    
!!$           ipcb(4,ii)=celllistb(1,ipcb(1,ii), ipcb(2,ii), ipcb(3,ii))
!!$        enddo
!!$     enddo
!!$  enddo

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
subroutine check_bonds_lig(ncol,maxncol,maxnblob,nchainspercol,& 
     nblobsperchain,ligboundto,posblob,rcut_liglig,lbox)
  implicit none
  
  integer, intent(in) :: maxncol,nchainspercol,nblobsperchain,maxnblob,&
       ligboundto(2,nchainspercol,maxncol),ncol 
  real*8, intent(in) :: posblob(3,nblobsperchain,nchainspercol,maxncol),&
       rcut_liglig, lbox(3)
  integer :: ii, jj, kk, iiblob, iichain,jjblob,kkblob,&
       llblob,jjrec,kkrec,irec,jjchain,kkchain,jjcol,iicol,kkcol
  real*8 :: distij(3),distij2,rcut_liglig2
  
  rcut_liglig2=rcut_liglig**2

!!$  ! do lig-receptor bond check on receptor side
!!$  do irec=1,nrec
!!$     if (recboundto(1,irec) .gt. 0) then
!!$        jjchain=recboundto(1,irec)
!!$        jjcol=recboundto(2,irec)
!!$        kkrec=ligboundto(jjchain,jjcol)
!!$        if (irec .ne. kkrec) then
!!$           write(*,*) 'recside bonding error!!  rec i  points to ligand j in colloid j :', &
!!$                irec,jjchain,jjcol
!!$           write(*,*) ' but lig j in  j points to rec k : ',jjchain,jjcol,kkrec
!!$        endif   
!!$        
!!$        distij(1:2)=posrec(:,irec)-posblob(1:2,nblobsperchain,jjchain,jjcol)
!!$        distij(1:2)=distij(1:2)-lbox(1:2)*nint(distij(1:2)/lbox(1:2))
!!$        distij(3)=posblob(3,nblobsperchain,jjchain,jjcol)
!!$        distij2=sum(distij*distij)
!!$        if (distij2 .gt. rcut_liglig2) then
!!$       !    write(*,*) 'WARNING: lig-rec bond is present but distance is too great'
!!$       !    write(*,*) ' iicol ',jjcol ,'jjchain ',jjchain ,' kkrec ',irec
!!$       !    write(*,*) 'poslig ', posblob(1:3,nblobsperchain,jjchain,jjcol), 'posrec ', posrec(:,irec)
!!$       !    write(*,*) ' distij2 =', distij2   
!!$        endif
!!$      
!!$     endif
!!$  enddo
  
    ! do lig-receptor bond check on ligand side
  do iicol=1,ncol
     do iichain=1,nchainspercol
        jjcol=ligboundto(1,iichain,iicol)
        jjchain=ligboundto(2,iichain,iicol)
        if (jjcol .gt. 0) then
           kkcol=ligboundto(1,jjchain,jjcol)
           kkchain=ligboundto(2,jjchain,jjcol)
           if (any((/iichain,iicol/) .ne. (/kkchain,kkcol/) )) then
              write(*,*) 'ligside bonding error!!  lig i on col i  points to lig j on col j :', iichain,iicol,jjchain, jjcol
              write(*,*) ' but lig j on col j points to lig k on col k : ',jjchain,jjcol,kkchain,kkcol
           endif
        endif
     enddo
  enddo

end subroutine check_bonds_lig
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
subroutine output_conf(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,psiWL,hWL)  
  implicit none
  
  integer, intent(in) ::  maxncol,ncol,nchainspercol,nblobsperchain,&
       nrxspec,ichainspec(nchainspercol,maxncol),nWL_points
  integer*8, intent(in) :: nout, icycle,hWL(nWL_points)
  real*8, intent(in) :: poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,&
       maxncol), rWL_min,WLbinsize,psiWL(nWL_points)
  real*8 :: freeE(nWL_points)
  character*50, intent(in) :: outfilename
  character*50 :: filename
  character*4 :: char1,char2
  character*6 :: char3,char4
  integer :: irec,icol,ichain,iblob,ii,ibin
  ! VMD OUTPUT FILE
  write(char1,'(I4)') int(icycle/nout)+1000 
!  write(char2,'(I4)') nint(mu)+2000
!  write(char3,'(I6)') ncol+100000
!  write(char4,'(I6)') nrec+100000
  !write(*,*) 'OK'
  filename=trim(outfilename)//'_cnf_out'//char1//'.xyz'
  
  open(unit=222,file=filename)
  write(222,*) ncol*(nblobsperchain*nchainspercol + 1)
  write(222,*) ncol,  '2 Hairy nanoparticles interaction'
  
  do icol=1,ncol         
     write(222,"(A,F8.3,F8.3,F8.3)") 'Xe  ', poscol(:,icol)
     do ichain=1,nchainspercol
        write(222,"(A,F8.3,F8.3,F8.3)") 'H ', posblob(:,1,ichain,icol) ! chain anchor on particle
        do iblob=2,nblobsperchain-1
           write(222,"(A,F8.3,F8.3,F8.3)") 'O  ', posblob(:,iblob,ichain,icol)  ! inert blobs         
        enddo
        if (ichainspec(ichain,icol) .eq. 1) then
           write(222,"(A,F8.3,F8.3,F8.3)") 'N', posblob(:,nblobsperchain,ichain,icol)   ! ligand   
        elseif (ichainspec(ichain,icol) .eq. 2) then   
           write(222,"(A,F8.3,F8.3,F8.3)") 'F', posblob(:,nblobsperchain,ichain,icol)   ! ligand
        else
           write(222,"(A,F8.3,F8.3,F8.3)") 'Se', posblob(:,nblobsperchain,ichain,icol)   ! ligand
        endif
     enddo
  enddo
  
  close(222)
  ! END VMD OUTPUT FILE



  filename=trim(outfilename)//'_WLFE_out'//char1//'.dat'

  freeE(:)=psiWL(:)-log(dble(hWL(:)))  ! get free energy
  freeE(:)=freeE(:)-freeE(1)  ! normalise by 1st bin

  open(unit=222,file=filename)

  do ibin=1,nWL_points
     write(222,*) ibin-1,'  ',freeE(ibin)
  enddo

  close(222) ! end WL writing
  
  
end subroutine output_conf
! ==================================================================== !

subroutine output_conf_lammps(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,psiWL,hWL,&
     lbox,ligboundto)  
  implicit none
  
  integer, intent(in) ::  maxncol,ncol,nchainspercol,nblobsperchain,&
       nrxspec,ichainspec(nchainspercol,maxncol),nWL_points,&
       ligboundto(2,nchainspercol,maxncol)
  integer*8, intent(in) :: nout, icycle,hWL(nWL_points)
  real*8, intent(in) :: poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,&
       maxncol), rWL_min,WLbinsize,psiWL(nWL_points),lbox(3)
  real*8 :: freeE(nWL_points)
  character*50, intent(in) :: outfilename
  character*50 :: filename, FOA, FOB
  character*4 :: char1,char2
  character*6 :: char3,char4
  integer :: irec,icol,ichain,iblob,ii,ibin,ind, bondsave(3,maxncol*nchainspercol*nblobsperchain), &
       bind, jcol,jchain,jind

!!$  if(nblobsperchain .lt. 2) then
!!$     write(*,*) 'ERROR: nblobsperchain must be at least 2 for lammps output to work.. exiting'
!!$     call exit()
!!$  endif
  
  FOA = "(I6,I6,I6,F9.3,F9.3,F9.3)" ! format output atoms
  FOB = "(I6,I6,I6,I6)" ! format output bonds
  ! VMD OUTPUT FILE
  write(char1,'(I4)') int(icycle/nout)+1000 
!  write(char2,'(I4)') nint(mu)+2000
!  write(char3,'(I6)') ncol+100000
!  write(char4,'(I6)') nrec+100000
  !write(*,*) 'OK'
  filename=trim(outfilename)//'_cnf_out'//char1//'.lammps'

  bind=0
  ! get the total number of bonds
  do icol=1,ncol
     do ichain=1,nchainspercol
        if (ligboundto(1,ichain,icol) .gt. 0 ) then  ! ADD BOND
           ! don't double count
           if (icol .gt. ligboundto(1,ichain,icol)) then
              bind=bind+1
           elseif ((icol .eq. ligboundto(1,ichain,icol)).and.(ichain .gt. ligboundto(2,ichain,icol))) then
              bind=bind+1
           endif
        endif
     enddo
  enddo

  
  open(unit=222,file=filename)
  write(222,"(A)") 'Lammps out for doublexp FreeE calcs  '
  write(222,"(A)") ' '
  write(222,*) ncol*(nblobsperchain*nchainspercol + 1), '    atoms'
  write(222,*) ncol*nchainspercol*(nblobsperchain - 1) + bind, '    bonds'
  write(222,"(A)")' '
  write(222,*) 3+nrxspec, '    atom types'
  write(222,*) 3, '    bond types'
  write(222,"(A)") ' '
  write(222,"(F9.3,F9.3,A)") 0.0, lbox(1), ' xlo xhi'
  write(222,"(F9.3,F9.3,A)") 0.0, lbox(2), ' ylo yhi'
  write(222,"(F9.3,F9.3,A)") 0.0, lbox(3), ' zlo zhi'
  write(222,"(A)")' '
  write(222,"(A)")' '
  write(222,"(A)")'Atoms '
  write(222,"(A)")' '
  ind=0 ! atom index
  bind=0 ! bond index
  do icol=1,ncol
     ind=ind+1
     write(222,FOA) ind, icol, 1,  poscol(:,icol)
     do ichain=1,nchainspercol
        if (nblobsperchain .gt. 1) then
           ind=ind+1
           write(222,FOA) ind, icol, 2, posblob(:,1,ichain,icol) ! chain anchor on particle
           do iblob=2,nblobsperchain-1
              ind=ind+1
              write(222,FOA) ind, icol, 3,  posblob(:,iblob,ichain,icol)  ! inert blobs
              ! BONDS SAVE
              if (iblob .eq. 2) then
                 bind=bind+1
                 bondsave(1:3,bind)=(/1,ind,ind-1/) ! anchor to blob bond type 1
              else
                 bind=bind+1
                 bondsave(1:3,bind)=(/2,ind,ind-1/) ! blob to blob bond type 2
              endif
              
           enddo
        endif !nblobsperchain > 1

        
        ind=ind+1
        write(222,FOA) ind, icol, ichainspec(ichain,icol)+ 3, posblob(:,nblobsperchain,ichain,icol)   ! ligand

        
        ! add ligand bond if present
        if (ligboundto(1,ichain,icol) .gt. 0 ) then  ! ADD BOND
           ! don't double count
           if (icol .eq. 2 ) then  ! only works for two colloids
              bind=bind+1

          !    write(*,*)
              jind= 1 + (ligboundto(2,ichain,icol)-0)*nblobsperchain
          !    write(*,*) ligboundto(2,ichain,icol), ligboundto(1,ichain,icol)
          !    write(*,*) ind, ichain, icol
              
              bondsave(1:3,bind) = (/2,ind,jind/)
           endif
        endif
     enddo
  enddo

  write(222,"(A)") ' ' ; write(222,"(A)") ' ' ;  write(222,"(A)") 'Bonds'
  write(222,"(A)") ' '  
  ! WRITE BONDS
  do ii=1,bind
     write(222,FOB) ii, bondsave(1:3,ii)  
  enddo
  write(222,*) ' ' ! add a blank line
  close(222)
!!!!!!!!!!!!!!!!   END VMD OUTPUT FILE

  filename=trim(outfilename)//'_WLFE_out'//char1//'.dat'

  freeE(:)=psiWL(:)-log(dble(hWL(:)))  ! get free energy
  freeE(:)=freeE(:)-freeE(1)  ! normalise by 1st bin

  open(unit=222,file=filename)

  do ibin=1,nWL_points
     write(222,*) ibin-1,'  ',freeE(ibin)
  enddo

  close(222) ! end WL writing
  
  
end subroutine output_conf_lammps

! ==================================================================== !

subroutine output_conf_trj(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     nout,icycle,outfilename,nrxspec,ichainspec,nWL_points,rWL_min,WLbinsize,psiWL,hWL,&
     lbox,ligboundto)  
  implicit none
  
  integer, intent(in) ::  maxncol,ncol,nchainspercol,nblobsperchain,&
       nrxspec,ichainspec(nchainspercol,maxncol),nWL_points,&
       ligboundto(2,nchainspercol,maxncol)
  integer*8, intent(in) :: nout, icycle,hWL(nWL_points)
  real*8, intent(in) :: poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,&
       maxncol), rWL_min,WLbinsize,psiWL(nWL_points),lbox(3)
  real*8 :: freeE(nWL_points)
  character*50, intent(in) :: outfilename
  character*50 :: filename, FOAL
  character*4 :: char1,char2
  character*6 :: char3,char4
  integer :: irec,icol,ichain,iblob,ii,ibin,ind, bondsave(3,maxncol*nchainspercol*nblobsperchain), &
       bind
  logical :: exists

  write(char1,'(I4)') int(icycle/nout)+1000
  
  filename=trim(outfilename)//'_dump.lammpstrj'
  inquire(file=filename,exist=exists) ! check if file exists
  if (exists) then
     open(555,file=filename,status="old",position="append")
  else
     open(555,file=filename,status="new",action="write")
  endif
 
  
  FOAL = "(I6,I6,F9.3,F9.3,F9.3)" ! format output atoms
 
  write(555,"(A)") 'ITEM: TIMESTEP  '
  write(555,*) icycle
  write(555,"(A)") 'ITEM: NUMBER OF ATOMS'
  write(555,*) ncol*nchainspercol*(nblobsperchain) + ncol
  write(555,"(A)") 'ITEM: BOX BOUNDS pp pp pp'
  write(555,*) 0, lbox(1)
  write(555,*) 0, lbox(2)
  write(555,*) 0, lbox(3)
  write(555,"(A)") 'ITEM: ATOMS id type xs ys zs '
  ind=0 ! atom index
  do icol=1,ncol
     ind=ind+1
     write(555,FOAL) ind, 1,  poscol(:,icol)/lbox(:)
     do ichain=1,nchainspercol
        ind=ind+1
        if (nblobsperchain .gt. 1) then
           write(555,FOAL) ind, 2, posblob(:,1,ichain,icol)/lbox(:) ! chain anchor on particle
        
           do iblob=2,nblobsperchain-1
              ind=ind+1
              write(555,FOAL) ind, 3,  posblob(:,iblob,ichain,icol)/lbox(:)  ! inert blobs
              
           enddo
        endif
        ind=ind+1
        write(555,FOAL) ind, ichainspec(ichain,icol)+ 3, &
             posblob(:,nblobsperchain,ichain,icol)/lbox(:)   ! ligand
        
     enddo
  enddo
  close(555)
!!!!!!!!!!!!!!!!   END VMD OUTPUT FILE

  filename=trim(outfilename)//'_WLFE_out'//char1//'.dat'

  freeE(:)=psiWL(:)-log(dble(hWL(:)))  ! get free energy
  freeE(:)=freeE(:)-freeE(1)  ! normalise by 1st bin

  open(unit=222,file=filename)

  do ibin=1,nWL_points
     write(222,*) ibin-1,'  ',freeE(ibin)
  enddo

  close(222) ! end WL writing

end subroutine output_conf_trj

! ========================================================================== !
subroutine output_conf_wlbonds(poscol,posblob,maxncol,nchainspercol,nblobsperchain,ncol,&
     mu,nout,icycle,outfilename,nrxspec,ichainspec,&
     psiWL,hWL,doWL_flag,nWLpoints)  
  implicit none
  
  integer, intent(in) ::  maxncol,ncol,nchainspercol,nblobsperchain,&
       nrxspec,ichainspec(nchainspercol),nWLpoints
  integer*8, intent(in) :: nout, icycle,hWL(nWLpoints)
  real*8, intent(in) :: poscol(3,maxncol),posblob(3,nblobsperchain,nchainspercol,maxncol),&
        mu, psiWL(nWLpoints)
  character*50, intent(in) :: outfilename
  logical, intent(in) :: doWL_flag
  character*50 :: filename
  character*4 :: char1,char2
  character*4 :: char3,char4
  integer :: icol,ichain,iblob,ii,ibin
  real*8 :: binsize,freeE(nWLpoints)
  ! VMD OUTPUT FILE
  write(char1,'(I4)') int(icycle/nout)+1000 
  write(char2,'(I4)') nint(mu)+2000
  write(char3,'(I4)') ncol+1000
  write(char4,'(I4)') nWLpoints+1000
  !write(*,*) 'OK'
  filename=trim(outfilename)//'WLbonds_nWLp'//char4//'_ncol'//char3//'_out'//char1//'.xyz'


  freeE(:)=psiWL(:)-log(dble(hWL(:)))  ! get free energy
  freeE(:)=freeE(:)-freeE(1)  ! normalise by 1st bin

  open(unit=222,file=filename)

  if (.not. doWL_flag) then ! only write the free energy datafiles if WL has finished.
     write(*,*) 'bFreeE ',-log(sum(exp(-freeE(2:nWLpoints))))   ! ibin=ibond
     do ibin=1,nWLpoints
        write(222,*) ibin-1,'  ',freeE(ibin)   ! ibin=ibond
     enddo
  endif
  
  close(222) ! end WL writing
  filename=trim(outfilename)//'conf_nWLp'//char4//'_ncol'//char3//'_out'//char1//'.xyz'
  open(unit=222,file=filename)
  
  write(222,*) ncol*nblobsperchain*nchainspercol
  write(222,*) ncol, ' Hairy nanoparticle adsorptopn'
  
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
  
!!$  do irec=1,nrec
!!$!     if (recboundto(1,irec) .gt. 0) then !this receptor is bound to someone
!!$!        write(222,"(A,F8.3,F8.3,F8.3)") 'F ', posrec(:,irec) , 0.0
!!$!     else
!!$!        write(222,"(A,F8.3,F8.3,F8.3)") 'H ', posrec(:,irec) , 0.0
!!$!     endif
!!$     if (irecspec(irec) .eq. 1) then !this receptor is specie 1
!!$        write(222,"(A,F8.3,F8.3,F8.3)") 'N ', posrec(:,irec) , 0.0
!!$     else
!!$        write(222,"(A,F8.3,F8.3,F8.3)") 'F ', posrec(:,irec) , 0.0
!!$     endif
!!$  enddo
  
  close(222)
  ! END VMD OUTPUT FILE
  
end subroutine output_conf_wlbonds

