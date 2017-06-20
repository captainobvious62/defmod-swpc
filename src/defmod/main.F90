! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

PROGRAM main

#include <petscversion.h>

  USE global
  USE FD2FE
  IMPLICIT NONE
#include "petsc.h"

  CHARACTER(256) :: input_file,viz,fd,exofile,debug
  LOGICAL :: l,v,w,db,em
  INTEGER,POINTER :: null_i=>null()
  REAL(8),POINTER :: null_r=>null()
  REAL(8) :: fdt
  INTEGER :: i,j,j1,j2,j3,j4,j5,j6,n,n_dyn,nodal_bw,ef_eldof,ng,vout,fdout,dbout

  ! Petsc Meshing/DMDA Stuff
  DM        :: dm,cdm
  Vec       :: Vec_Coordinates
  INTEGER   :: Int_elem, coord_len,exo_dof
  REAL(kind=8), ALLOCATABLE :: PETSc_coords(:,:), PETSc_coords_array(:)
  PetscBool :: interpolate
  REAL(8),POINTER :: PETSc_pntr(:)  
  ! Global Vectors
  Vec :: GVec_displacement,GVec_rhob

  ! Local Vectors
  Vec :: LVec_rhob

  ! Viewer Related
  PetscViewer :: ASCIICoordViewer,ASCIIMeshViewer,ASCIIFSLabelViewer,           &
       ASCIIVSLabelViewer,ASCIICSLabelViewer,ASCIIISViewer,ASCIIFEViewer

  ! Mesh interval
  !PetscInt :: pStart,pEnd

  ! Mesh Sets
  DMLabel :: FSLabel,VSLabel,CSLabel

  ! Index Sets
  IS :: IS_LabelIDs,IS_FaceSets,IS_VertexSets,IS_CellSets

  ! Nodesets, sidesets, etc.
  PetscInt :: num_ns,num_ss,num_bl,labelsize

  ! DM Test
  REAL(8), ALLOCATABLE :: DMCoords(:)
  PetscInt :: numFields,numBC
  PetscInt, TARGET, DIMENSION(3) :: numComp

  ! PETScFE Test
  PetscDS :: DS_Prob
  PetscFE :: FE_fem(2)
  PetscQuadrature :: PQ_q
  PetscInt :: order

  

  ! Debug QS
!  Vec :: Vec_U_db,Vec_F_db
!  Mat :: Mat_K_db
  !-----------------------------------------------------------------------------80  
  CALL PetscInitialize(Petsc_Null_Character,ierr)
  CALL MPI_Comm_Rank(MPI_Comm_World,rank,ierr)
  CALL MPI_Comm_Size(MPI_Comm_World,nprcs,ierr)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=6)
  CALL PetscOptionsGetString(Petsc_Null_Character,'-f',input_file,l,ierr)
  CALL PetscOptionsGetString(Petsc_Null_Character,'-ss',viz,v,ierr)
  CALL PetscOptionsGetString(Petsc_Null_Character,'-fd',fd,w,ierr)
  ! Flag to load .exo directly
  CALL PetscOptionsGetString(Petsc_Null_Character,'-m',exofile,em,ierr)
  ! Debug (no fault) flag
  CALL PetscOptionsGetString(Petsc_Null_Character,'-db',debug,db,ierr)

#else
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-f',      &
       input_file,l,ierr)
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-ss',     &
       viz,v,ierr)
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-fd',     &
       fd,w,ierr)
  ! Flag to load .exo directly     
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-m' ,     &
       exofile,em,ierr)     
  ! Debug (no fault) flag
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-db',     &
       debug,db,ierr)
#endif
  
  ! Read Command Line Flags====================================================80
  ! Respond if no input
  IF (.NOT. l) THEN
     CALL PrintMsg("Usage: [mpiexec -n <np>] defmod -f <input_filename>")
     go to 9
  END IF
  ! Produce VTK output  
  vout=0
  IF (v) THEN
     READ (viz,'(I1.0)')vout 
     IF (vout==1) THEN 
        CALL PrintMsg("Snapshot output will slow the run!")
     ELSE
        CALL PrintMsg("Chosen NOT to output snapshot.")
     END IF
  ELSE
     CALL PrintMsg("Use -ss 1 to turn on the snapshot output.")
  END IF
  ! Mixed FE/FD mode
  fdout=0
  IF (w) THEN
     READ (fd,'(I1.0)')fdout
     IF (fdout==1) CALL PrintMsg("Runing in FE-FD mixed mode.")
  ELSE
     CALL PrintMsg("Use -fd 1 to turn on FE-FD mixed mode.")
  END IF
  ! Debug (no fault) mode
  dbout=0
  IF (db) THEN
     READ (debug,'(I1.0)')dbout
     IF (dbout==1) CALL PrintMsg("Running in debug mode.")
  ELSE
     CALL PrintMsg("Use -db 1 to turn on debug (no fault) mode.")
  END IF
  !---------------------------------------------------------------------------80

  ! Read input file parameters
  OPEN(10,file=input_file,status='old')

  ! Output file name is the same as the input
  IF (INDEX(input_file,".inp")>0) THEN
     output_file=input_file(1:INDEX(input_file,".inp")-1)
  ELSE
     output_file=input_file
  END IF

  CALL PrintMsg("Reading input ...")
  CALL ReadParameters

  ! Set element specific constants
  CALL InitializeElement
  p=0; ef_eldof=eldof
  IF (poro) THEN
     p=1; ef_eldof=eldof+eldofp
     IF (eltype=="tri") nip=3; IF (eltype=="tet") nip=4
  ELSEIF (fault) THEN
     IF (eltype=="tri") nip=3; IF (eltype=="tet") nip=4
  END IF


  !===========================================================================80
  ! Read/Test Direct .exo reading
  !---------------------------------------------------------------------------80
  IF (em) THEN
     interpolate=.TRUE.  
     CALL DMPlexCreateFromFile(Petsc_Comm_World, exofile, Petsc_True, dm, ierr)
     !CALL DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-orig_dm_view");CHKERRQ(ierr)

     PRINT *, exo_dof
     ! Read in mesh particulars, generate reference vectors
     CALL PrintMsg('Reading Nodal Coordinates from PETSc')
     
     ! Global Coordinates
     CALL DMGetCoordinates(dm,Vec_Coordinates,ierr);CHKERRQ(ierr)

     ! Local to Global DM
     CALL DMGetCoordinateDM(dm,cdm,ierr)
     
     CALL PrintMsg('Reading # of Coordinates')
     CALL VecGetSize(Vec_Coordinates,Int_elem,ierr)
     ! Deal with coords in FORTRAN Variables
     CALL PrintMsg('Converting to FORTRAN Variables')
     CALL DMGetDimension(dm,dmn,ierr)
     PRINT *, Int_elem, dmn
     coord_len = Int_elem/dmn
     ALLOCATE( PETSc_coords(coord_len,dmn),PETSc_coords_array(Int_elem),        &
          )
     CALL PrintMsg('Trying to read from Array')
     CALL VecGetArrayReadF90(Vec_Coordinates,PETSc_pntr,ierr)
     PETSc_coords_array=PETSc_pntr
     ! Reshape to (nnds,dmn)
     PETSc_coords=RESHAPE(PETSc_pntr,(/dmn,coord_len/))
     CALL PrintMsg('Restoring Data')
     CALL VecRestoreArrayReadF90(Vec_Coordinates,PETSc_pntr,ierr)
     ! Write to test output file
     ! Requires PETSc Reference Section
     ! CALL DMCreateGlobalVector(dm,GVec_rhob,ierr)

     ! Read in nodesets, sidesets, and such (i hope) from mesh
     CALL DMGetLabel(dm,'Face Sets',FSLabel,ierr)
     CALL DMGetLabel(dm,'Vertex Sets',VSLabel,ierr)
     CALL DMGetLabel(dm,'Cell Sets',CSLabel,ierr)

     ! Again with the DMLabels, this time to index sets
     CALL DMGetLabelIdIS(dm,'Vertex Sets',IS_LabelIDs,ierr)

     ! Next step is to read in node set and side set values directly. This
     ! will take more work/drawing out on paper.
     ! Read in sidesets
     CALL DMGetLabelSize(dm,'Face Sets',labelsize,ierr)

     ! Try to read sets directly/Index Set holds node numbers
     CALL DMGetStratumIS(dm,'Vertex Sets',1,IS_VertexSets,ierr)


     ! Debug PETSc output
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'coordPETSc.out',ASCIICoordViewer,ierr)
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'mesh.out',ASCIIMeshViewer,ierr)
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'FSlabel.out',ASCIIFSLabelViewer,ierr)    
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'VSlabel.out',ASCIIVSLabelViewer,ierr)
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'CSlabel.out',ASCIICSLabelViewer,ierr)
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'ISlabel.out',ASCIIISViewer,ierr)
     CALL PetscViewerASCIIOpen(MPI_Comm_World,'PetscFE.out',ASCIIFEViewer,ierr)
     CALL VecView(Vec_Coordinates,ASCIICoordViewer,ierr)
     CALL DMView(dm, ASCIIMeshViewer, ierr)
     CALL DMLabelView(FSLabel, ASCIIFSLabelViewer,ierr)
     CALL DMLabelView(VSLabel, ASCIIVSLabelViewer,ierr)
     CALL DMLabelView(CSLabel, ASCIICSLabelViewer,ierr)
     CALL ISView(IS_VertexSets, ASCIIISViewer, ierr)
     OPEN(14,file='coordFORTRAN.out',status='replace')
     DO i=1,coord_len
        WRITE(14,*)PETSc_coords(:,i)
     END DO
     CLOSE(14)

     !CALL PetscIntView(1,Int_elem,PETSC_VIEWER_STDOUT_WORLD,ierr)

     ! Attempt to generate PetscFE Framework
     CALL PetscFECreateDefault(dm,dmn,dmn,PETSC_FALSE,"u_",-1,FE_fem(1),ierr)
     CALL PetscObjectSetName(FE_fem(1), 'displacement',ierr)
     CALL PetscFEGetQuadrature(FE_fem(1),PQ_q,ierr)
     CALL PetscQuadratureGetOrder(PQ_q,order,ierr)
     CALL PetscFECreateDefault(dm,dmn,1,PETSC_FALSE,"p_",order,FE_fem(2),ierr)
     CALL PetscObjectSetName(FE_fem(2), 'pressure',ierr)

!     DO WHILE (cdm)
!        CALL DMGetDS(cdm,DS_Prob,ierr)
        

     ! Clean up
     CALL PetscViewerDestroy(ASCIICoordViewer,ierr)
     CALL PetscViewerDestroy(ASCIIMeshViewer,ierr)     
     CALL PetscViewerDestroy(ASCIIFSLabelViewer,ierr)
     CALL PetscViewerDestroy(ASCIIVSLabelViewer,ierr)          
     CALL PetscViewerDestroy(ASCIICSLabelViewer,ierr)     
     CALL PetscViewerDestroy(ASCIIISViewer,ierr)     
     CALL PetscViewerDestroy(ASCIIFEViewer,ierr)     
  END IF
  !-----------------------------------------------------------------------------80



  ! Partition mesh using METIS, create mappings, and read on-rank mesh data
  ALLOCATE(npart(nnds),epart(nels)); epart=0; npart=0
  IF (nprcs>1) THEN
     CALL PrintMsg("Partitioning mesh ...")
     IF (rank==0) THEN
        ALLOCATE(nodes(1,npel*nels),work(nels+1)); work(1)=0
        DO i=1,nels
           j=npel*(i-1)+1; n=npel*i; READ(10,*)nodes(1,j:n); work(i+1)=n
        END DO
        nodes=nodes-1
        CALL METIS_PartMeshNodal(nels,nnds,work,nodes,null_i,null_i,nprcs,     &
             null_r,null_i,n,epart,npart)
        DEALLOCATE(nodes,work)
        ! Return to beginning of input file; read parameters in again
        REWIND(10); CALL ReadParameters
     END IF
     CALL MPI_Bcast(npart,nnds,MPI_Integer,0,MPI_Comm_World,ierr)
     CALL MPI_Bcast(epart,nels,MPI_Integer,0,MPI_Comm_World,ierr)
  END IF
  CALL PrintMsg("Reading mesh data ...")
  ALLOCATE(emap(nels),nmap(nnds)); emap=0; nmap=0
  ! Create original to local element mappings and read on-rank element data
  j=1
  DO i=1,nels
     IF (epart(i)==rank) THEN
        epart(i)=1; emap(i)=j; j=j+1
     ELSE
        epart(i)=0
     END IF
  END DO
  n=SUM(epart); ALLOCATE(nodes(n,npel),id(n)) ! id -> mtrl flag
  j=1
  ! Read nodes associated with element (sequentially), plus elemental
  ! material flag
  DO i=1,nels
     IF (epart(i)==1) THEN
        READ(10,*)nodes(j,:),id(j); j=j+1
     ELSE
        READ(10,*)val
     END IF
  END DO
  nels=n
  ! Create original to global nodal mappings and read on-rank + ghost nodes
  ALLOCATE(work(0:nprcs-1))
  j=0
  DO i=1,nnds
     IF (npart(i)==rank) j=j+1
  END DO
  CALL MPI_AllGather(j,1,MPI_Integer,work,1,MPI_Integer,MPI_Comm_World,ierr)
  IF (rank==0) n=1
  IF (rank/=0) n=SUM(work(0:rank-1))+1
  DO i=1,nnds
     IF (npart(i)==rank) THEN
        nmap(i)=n; n=n+1
     END IF
  END DO
  DEALLOCATE(work)
  ALLOCATE(work(nnds))
  CALL MPI_AllReduce(nmap,work,nnds,MPI_Integer,MPI_Sum,MPI_Comm_World,ierr)
  nmap=work
  npart=0
  DO i=1,nels
     DO j=1,npel
        npart(nodes(i,j))=1
     END DO
  END DO
  j=1; work=0
  DO i=1,nnds
     IF (npart(i)==1) THEN
        work(i)=j; j=j+1
     END IF
  END DO
  n=SUM(npart)
  ! Flag to enable Specified Dirichlet Pressure (to be moved to input file)
  DPFlag = 1 
  !ALLOCATE(coords(n,dmn),bc(n,dmn+p))
  ! BC row now includes specified pore pressure (hence the + 1)/rx_press 
  ! variable
  ALLOCATE(coords(n,dmn),bc(n,dmn+p),rx_press(n))
  j=1
  ! Read [sequentially] nodal coordinates and BC flags: x,y,z and/or p
  DO i=1,nnds
     IF (npart(i)==1) THEN
        IF (DPFlag==1) THEN
           READ(10,*)coords(j,:),bc(j,:),rx_press(j); j=j+1
        ELSE
           READ(10,*)coords(j,:),bc(j,:); j=j+1
        END IF
     ELSE
        READ(10,*)val
     END IF
  END DO
  coords=km2m*coords
  ! Re-number on-rank nodes and create local to global node and dof mappings
  DO i=1,nels
     DO j=1,npel
        nodes(i,j)=work(nodes(i,j))
     END DO
  END DO
  n=SUM(npart); ALLOCATE(nl2g(n,2),indxmap((dmn+p)*n,2))
  IF (fault) ALLOCATE(indxmap_u(dmn*n,2))
  j=1
  DO i=1,nnds
     IF (work(i)==j) THEN
        nl2g(j,1)=j; nl2g(j,2)=nmap(i); j=j+1
     END IF
  END DO
  DO i=1,n
     DO j=1,dmn+p
        indxmap((dmn+p)*i-j+1,:)=(dmn+p)*nl2g(i,:)-j ! 0 based index
     END DO
  END DO
  IF (fault) THEN
     DO i=1,n
        DO j=1,dmn
           indxmap_u(dmn*i-j+1,:)=dmn*nl2g(i,:)-j
        END DO
     END DO
  END IF
  DEALLOCATE(work)
  ! Read material data assuming nels >> nmts, i.e., all ranks store all data
  IF (fault) THEN
     ALLOCATE(mat(nmts,5+4*p+init+2))
  ELSE
     ALLOCATE(mat(nmts,5+4*p))
  END IF
  DO i=1,nmts
     READ(10,*)mat(i,:)
  END DO
  DEALLOCATE(epart,npart)
  !---------------------------------------------------------------------------80

  ! Build Initial Formulation
  !===========================================================================80

  ! Initialize local element variables and global U
  ALLOCATE(ipoint(nip,dmn),weight(nip),k(ef_eldof,ef_eldof),m(eldof,eldof),    &
       f(ef_eldof),indx(ef_eldof),enodes(npel),ecoords(npel,dmn),vvec(dmn+p))
  IF (fault) ALLOCATE(k_dyn(eldof,eldof),indx_dyn(eldof))
  CALL SamPts(ipoint,weight)
  n=(dmn+p)*nnds; IF (stype/="explicit") n=n+nceqs
  CALL VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
  ! Global U for dynamic run
  IF (fault) THEN
     n_dyn=dmn*nnds
     CALL VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n_dyn,Vec_U_dyn,ierr)
  END IF
  IF (visco .OR. (fault .AND. lm_str==1)) ALLOCATE(stress(nels,nip,cdmn))

  ! Set scaling constants
  wt=f2*EXP((LOG(MAXVAL(mat(:,1)))+LOG(MINVAL(mat(:,1))))/f2)
  IF (dmn==3) wt=wt*km2m
  IF (poro) THEN
     scale=SQRT(EXP((LOG(MAXVAL(mat(:,1)))+LOG(MINVAL(mat(:,1))))/f2)/         &
          EXP((LOG(MAXVAL(mat(:,6)))+LOG(MINVAL(mat(:,6))))/f2)/         &
          dt/km2m)
  END IF

  ! Form stiffness matrix (K)
  CALL PrintMsg("Forming [K] ...")
  nodal_bw=(dmn+p)*(nodal_bw+1)
  CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,   &
       Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)
  CALL MatSetOption(Mat_K,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)

  IF (fault) THEN ! Dynamic K
     nodal_bw=nodal_bw/(dmn+p)*dmn
     CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n_dyn,n_dyn, &
          nodal_bw,Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K_dyn,ierr)
     CALL MatSetOption(Mat_K_dyn,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
          ierr)
  END IF

  ! Create vector for RHS - Used for case of Dirichlet Pressure
  CALL VecDuplicate(Vec_U,Vec_F,ierr)
  f = f0

  DO i=1,nels
     ! Form elemental K and M matricies
     IF (poro) THEN
        IF (DPFlag == 1) THEN
           ! Use Dirichlet Pressure form
           CALL FormLocalK_DP(i,k,f,indx,"Kp")
        ELSE
           CALL FormLocalK(i,k,indx,"Kp")
        END IF
     ELSE
        CALL FormLocalK(i,k,indx,"Ke")
     END IF
     indx=indxmap(indx,2)
     CALL MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)

     ! Adjust RHS Values to account for dirichlet pressure
     CALL VecSetValues(Vec_F,ef_eldof,indx,f,Add_Values,ierr)

     ! Separate dynamic K matrix for IMEX (elastodynamic only)
     IF (fault) THEN
        dyn=.TRUE.
        CALL FormLocalK(i,k_dyn,indx_dyn,"Ke")
        indx_dyn=indxmap_u(indx_dyn,2)
        CALL MatSetValues(Mat_K_dyn,eldof,indx_dyn,eldof,indx_dyn,k_dyn,       &
             Add_Values,ierr)
        dyn=.FALSE.
     END IF
  END DO

  ! Initialize and form mass matrix and its inverse
  IF (stype=="explicit" .OR. fault) THEN
     CALL PrintMsg("Forming [M] & [M]^-1 ...")
     IF (fault) THEN
        CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n_dyn,    &
             n_dyn,3,Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_M,ierr)
     ELSEIF (gf) THEN
        CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,3,    &
             Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_M,ierr)
     ELSE
        CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,1,    &
             Petsc_Null_Integer,0,Petsc_Null_Integer,Mat_M,ierr)
     END IF
     CALL MatSetOption(Mat_M,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
     DO i=1,nels
        IF (fault) THEN
           CALL FormLocalM(i,m,indx_dyn)
           indx_dyn=indxmap_u(indx_dyn,2)
           DO j=1,eldof
              val=m(j,j)
              CALL MatSetValue(Mat_M,indx_dyn(j),indx_dyn(j),val,Add_Values,   &
                   ierr)
           END DO
        ELSE
           CALL FormLocalM(i,m,indx)
           indx=indxmap(indx,2)
           DO j=1,eldof
              val=m(j,j)
              CALL MatSetValue(Mat_M,indx(j),indx(j),val,Add_Values,ierr)
           END DO
        END IF
     END DO
     CALL MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     CALL MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
     CALL MatDuplicate(Mat_M,Mat_Do_Not_Copy_Values,Mat_Minv,ierr)
     IF (fault) THEN
        CALL MatGetDiagonal(Mat_M,Vec_U_dyn,ierr) ! Vec_U_dyn -> work vector
        CALL VecReciprocal(Vec_U_dyn,ierr)
        CALL MatDiagonalSet(Mat_Minv,Vec_U_dyn,Insert_Values,ierr)
        CALL VecZeroEntries(Vec_U_dyn,ierr)
     ELSE
        CALL MatGetDiagonal(Mat_M,Vec_U,ierr) ! Vec_U is used as a work vector
        CALL VecReciprocal(Vec_U,ierr)
        CALL MatDiagonalSet(Mat_Minv,Vec_U,Insert_Values,ierr)
        CALL VecZeroEntries(Vec_U,ierr)
     END IF
     CALL MatScale(Mat_M,alpha,ierr)
  END IF

  ! Allocate arrays to store loading history
  ALLOCATE(cval(nceqs,3),fnode(nfrcs),fval(nfrcs,dmn+p+2),telsd(ntrcs,2),      &
       tval(ntrcs,dmn+p+2),rxval(ntrcs,dmn+p+2),rxnode(nrxs))
  cval=f0; fval=f0; tval=f0; rxval=f0

  ! Account for constraint eqn's
  IF (nceqs>0) THEN
     CALL PrintMsg("Applying constraints ...")
     IF (stype=="explicit") THEN
        CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,        &
             nceqs,3,Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_Gt,ierr)
        CALL MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
             ierr)
     ELSEIF(fault .AND. nceqs-nceqs_ncf>0 .AND. hyb>0) THEN
        CALL VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nceqs_ncf/(dmn+p)+nfnd,&
             Vec_lmnd,ierr)
        CALL VecGetLocalSize(Vec_lmnd,n_lmnd,ierr)
        CALL VecGetOwnershipRange(Vec_lmnd,lmnd0,j,ierr)
        CALL VecDestroy(Vec_lmnd,ierr)
        ! Dofs of one fault node are not split by different ranks
        CALL VecCreateMPI(Petsc_Comm_World,n_lmnd*dmn,(nceqs_ncf/(dmn+p)+nfnd)*&
             dmn,Vec_lambda_sta,ierr)
        CALL VecDuplicate(Vec_lambda_sta,Vec_lambda_sta0,ierr)
        CALL VecZeroEntries(Vec_lambda_sta,ierr)
        CALL VecZeroEntries(Vec_lambda_sta0,ierr)
        CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,n_lmnd*dmn,dmn*nnds,   &
             (nceqs_ncf/(dmn+p)+nfnd)*dmn,5,Petsc_Null_Integer,5,                &
             Petsc_Null_Integer,Mat_Gt,ierr)
        CALL MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
             ierr)
        CALL MatZeroEntries(Mat_Gt,ierr)
        ALLOCATE(flt_slip(n_lmnd*dmn),tot_flt_slip(n_lmnd*dmn),                &
             qs_flt_slip(n_lmnd*dmn))
        qs_flt_slip=f0
     END IF
     ! Open cnstrns.tmp to store constraint equation data
     IF (rank==0) OPEN(15,file="cnstrns.tmp",status="replace")
     DO i=1,nceqs
        ! Read in number of constraint equation terms
        READ(10,*)n
        IF (rank==0) WRITE(15,*)n
        DO j=1,n
           READ(10,*)vvec,node; node=nmap(node)
           IF (rank==0) WRITE(15,*)REAL(vvec),node
        END DO
        READ(10,*)cval(i,:)
        IF (poro) THEN
           ! Scale constraint value if nonzero p coefficient
           IF (vvec(dmn+1)/=f0) cval(i,1)=cval(i,1)/scale
        END IF
        IF (rank==0) WRITE(15,*)REAL(cval(i,:))
     END DO
     IF (rank==0) CLOSE(15)

     ! Read fault orientation vectors
     IF (fault .OR. gf) THEN
        ALLOCATE(node_pos(nfnd),node_neg(nfnd),vecf(nfnd,dmn*dmn),fc(nfnd),    &
             fcd(nfnd),dc(nfnd),perm(nfnd),vvec_all(2*nfnd*dmn,dmn),             &
             node_all(2*nfnd*dmn),st_init(nfnd,dmn),xfnd(nfnd,dmn),frc(nfnd),    &
             coh(nfnd),dcoh(nfnd))    
        IF (rsf==1) ALLOCATE(rsfb0(nfnd),rsfV0(nfnd),rsfdtau0(nfnd),rsfa(nfnd),&
             rsfb(nfnd),rsfL(nfnd),rsftheta(nfnd))
        DO i=1,nfnd
           IF (poro) THEN
              IF (rsf==1) THEN
                 READ(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                      rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                      perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),              &
                      dcoh(i)
              ELSE
                 READ(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),    &
                      dc(i),perm(i),st_init(i,:),xfnd(i,:),frc(i),               &
                      coh(i),dcoh(i)
              END IF
           ELSE
              IF (rsf==1) THEN
                 READ(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                      rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                      st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              ELSE
                 READ(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),     &
                      dc(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              END IF
           END IF
           node_pos(i)=nmap(node_pos(i)); node_neg(i)=nmap(node_neg(i))
        END DO
     END IF

     ! Create variables for rotated constraint matrix
     IF (fault .OR. gf) THEN
        DO i=1,nfnd*dmn
           DO j=1,2
              READ(10,*)vvec(1:dmn),node; node=nmap(node)
              vvec_all(2*(i-1)+j,:)=vvec(1:dmn); node_all(2*(i-1)+j)=node
           END DO
        END DO
     END IF
     IF (rank==0) CALL ApplyConstraints
     IF (stype=="explicit" .AND. .NOT. gf) THEN
        CALL MatAssemblyBegin(Mat_Gt,Mat_Final_Assembly,ierr)
        CALL MatAssemblyEnd(Mat_Gt,Mat_Final_Assembly,ierr)
     END IF

     ! For dummy fault without split nodes
  ELSEIF (fault) THEN
     CALL PrintMsg("Dummy fault placed ...")
     ALLOCATE(node_pos(nfnd),vecf(nfnd,dmn*dmn),fc(nfnd),perm(nfnd))
     DO i=1,nfnd
        IF (poro) THEN
           READ(10,*) node_pos(i),vecf(i,:),fc(i),perm(i)
        ELSE
           READ(10,*) node_pos(i),vecf(i,:),fc(i)
        END IF
        node_pos(i)=nmap(node_pos(i))
     END DO
  END IF

  CALL MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
  CALL MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
  IF (fault .AND. nceqs>0) THEN
     CALL MatAssemblyBegin(Mat_K_dyn,Mat_Final_Assembly,ierr)
     CALL MatAssemblyEnd(Mat_K_dyn,Mat_Final_Assembly,ierr)
  END IF

  ! Form RHS
  CALL PrintMsg("Forming RHS ...")
  !CALL VecDuplicate(Vec_U,Vec_F,ierr)
  tstep=0
  ! Read in nodal force / fluid source values and time intervals
  DO i=1,nfrcs
     READ(10,*)fnode(i),fval(i,:)
  END DO
  ! Read in traction BC data; element ID, facet number, traction/flux value,
  ! and applied time interval
  DO i=1,ntrcs
     READ(10,*)telsd(i,:),tval(i,:)
  END DO
  ! Read in Perscribed Node BCs - These may have been a bad idea
  DO i=1,nrxs
     READ(10,*)rxnode(i),rxval(i,:)
  END DO

  ! Remap nodes/els
  fnode=nmap(fnode); telsd(:,1)=emap(telsd(:,1));rxnode=nmap(rxnode)

  CALL FormRHS
  CALL VecAssemblyBegin(Vec_F,ierr)
  CALL VecAssemblyEnd(Vec_F,ierr)

  ! Observation solution space
  IF (nobs>0) THEN
     ALLOCATE(ocoord(nobs,dmn))
     DO i=1,nobs
        READ(10,*) ocoord(i,:)
     END DO
     CALL GetObsNd("ob")
     IF (nobs_loc>0) THEN
        ALLOCATE(uu_obs(nobs_loc,dmn+p),tot_uu_obs(nobs_loc,dmn+p),            &
             uu_dyn_obs(nobs_loc,dmn),tot_uu_dyn_obs(nobs_loc,dmn))
        tot_uu_obs=f0
     END IF
     n_log_dyn=0
  END IF

  ! FD domain grid bocks containing fault nodes, xgp, idgp(_loc), gpl2g, gpnlst, gpshape
  IF (nceqs-nceqs_ncf>0 .AND. fdout==1) THEN
     CALL PrintMsg("Locating FD grid points ...")
     CALL FDInit
     CALL GetFDFnd 
     CALL GetObsNd("fd")
     DEALLOCATE(xgp,idgp)  
     IF (ngp_loc>0) ALLOCATE(uu_fd(ngp_loc,dmn))
     CALL NndFE2FD
     CALL MatFE2FD
  END IF

  !===========================================================================80
  ! Account for absorbing boundaries
  CALL PrintMsg("Absorbing boundary ...")
  IF (stype=="explicit" .OR. (fault .AND. nceqs>0)) THEN
     DO i=1,nabcs
        ! For axis aligned absorbing boundaries
        !read(10,*)el,side,j; el=emap(el)
        ! j (dir ID) has no use for arbitrarily facing AbsC 
        READ(10,*)el,side; el=emap(el)
        IF (el/=0) THEN
           IF (fault) THEN
              !call FormLocalAbsC(el,side,abs(j),m,indx_dyn)
              CALL FormLocalAbsC1(el,side,m,indx_dyn)
              indx_dyn=indxmap_u(indx_dyn,2)
              !do j=1,eldof
              !   val=m(j,j)
              !   call MatSetValue(Mat_M,indx_dyn(j),indx_dyn(j),val,           &
              !      Add_Values,ierr)
              !end do
              !call MatSetValues(Mat_M,eldof,indx_dyn,eldof,indx_dyn,m,         &
              !   Add_Values,ierr)
              DO j1=1,eldof
                 DO j2=1,eldof
                    val=m(j1,j2)
                    IF (ABS(val)>f0) CALL MatSetValue(Mat_M,indx_dyn(j1),      &
                         indx_dyn(j2),val,Add_Values,ierr)
                 END DO
              END DO
           ELSE
              !call FormLocalAbsC(el,side,abs(j),m,indx)
              CALL FormLocalAbsC1(el,side,m,indx)
              indx=indxmap(indx,2)
              !do j=1,eldof
              !   val=m(j,j)
              !   call MatSetValue(Mat_M,indx(j),indx(j),val,Add_Values,ierr)
              !end do
              !call MatSetValues(Mat_M,eldof,indx,eldof,indx,m,Add_Values,ierr)
              DO j1=1,eldof
                 DO j2=1,eldof
                    val=m(j1,j2)
                    IF (ABS(val)>f0) CALL MatSetValue(Mat_M,indx(j1),indx(j2), &
                         val,Add_Values,ierr)
                 END DO
              END DO
           END IF
        END IF
     END DO
     CALL MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     CALL MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
  END IF
  CLOSE(10)
  DEALLOCATE(nmap,emap) ! End of input reading

  ! Initialize arrays to communicate ghost node values
  CALL PrintMsg("Setting up solver ...")
  j=SIZE(indxmap,1)
  CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_U,ierr)
  CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,2),Petsc_Copy_Values,From,  &
       ierr)
  CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,1),Petsc_Copy_Values,To,    &
       ierr)
  CALL VecScatterCreate(Vec_U,From,Seq_U,To,Scatter,ierr)
  ALLOCATE(uu(j),tot_uu(j)); tot_uu=f0
  IF (fault) THEN
     IF (lm_str==1) THEN
        j=SIZE(indxmap_u,1)
        CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_SS,ierr)
        CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_SH,ierr)
        IF (nceqs>0) CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_f2s,ierr)
        CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_dip,ierr)
        CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_nrm,ierr)
     END IF
     j=SIZE(indxmap_u,1)
     CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_U_dyn,ierr)
     CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),Petsc_Copy_Values,  &
          From,ierr)
     CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),Petsc_Copy_Values,  &
          To,ierr)
     CALL VecScatterCreate(Vec_U_dyn,From,Seq_U_dyn,To,Scatter_dyn,ierr)
     ALLOCATE(uu_dyn(j),tot_uu_dyn(j)) 
     IF (lm_str==1) THEN
        ALLOCATE(ss(j),sh(j),f2s(j),dip(j),nrm(j))
     END IF
  END IF
  !===========================================================================80
  ! Debug QS
  IF (db) THEN
     CALL PrintMsg(" Nope, we ain't doin' that. ")
     CALL VecDestroy(Vec_U,ierr)
     CALL MatDestroy(Mat_K,ierr)

     ! Start debug run from scratch
     Call SamPts(ipoint,weight)
     n = (dmn+p)*nnds

     CALL VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
     CALL VecDuplicate(Vec_U,Vec_F,ierr)

     CALL MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,    &
          Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)

     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~80
     ! Generate Index Sets for U and P
     CALL VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
     CALL VecCopy(Vec_U,Vec_Um,ierr)
     CALL VecGetLocalSize(Vec_U,j,ierr)
     IF (poro) THEN
        ! Create Index Set Scheme for displacement, pressure, and LM locations from Vec_U
        j2=0; j3=0; j4=0; j5=0; j6=0
        DO i=1,j
           IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
              j2=j2+1
           END IF
           IF (nceqs>0) THEN
              IF (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) THEN
                 j3=j3+1
              END IF
           END IF
           IF (MOD(j1+i,dmn+1)>0 .AND. j1+i-1<(dmn+1)*nnds) THEN
              j4=j4+1
           END IF
           IF (j1+i-1<(dmn+1)*nnds) THEN
              j5=j5+1
           END IF
        END DO
        ALLOCATE(work(j2),workl(j3),worku(j4))
        j2=0; j3=0; j4=0; j5=0
        DO i=1,j
           IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
              j2=j2+1
              work(j2)=j1+i-1
           END IF
           IF (MOD(j1+i,dmn+1)>0 .AND. j1+i-1<(dmn+1)*nnds) THEN
              j4=j4+1
              worku(j4)=j1+i-1
           END IF
           IF (nceqs>0) THEN
              IF (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) THEN
                 j3=j3+1
                 workl(j3)=j1+i-1
              END IF
           END IF
        END DO
        j=SIZE(work)
        CALL ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,RI,  &
             ierr)
        j=SIZE(worku)
        CALL ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,    &
             RIu,ierr)
        CALL MatGetSubMatrix(Mat_K,RIu,RI,Mat_Initial_Matrix,Mat_H,ierr)
        IF (nceqs > 0) THEN
           j=SIZE(workl)
           CALL ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values, &
                RIl,ierr)
        END IF
     END IF
     ! END ABORTION OF INDEX SET GENERATION--------------------------------------80

     ! Generate Stiffness Matrix K-----------------------------------------------80
     call PrintMsg(" ReForming Stiffness Matrix K" )
     DO i=1,nels
        ! Form elemental K and M matricies
        IF (poro) THEN
           CALL FormLocalK(i,k,indx,"Kp")
        ELSE
           CALL FormLocalK(i,k,indx,"Ke")
        END IF
        indx=indxmap(indx,2)
        CALL MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
     END DO

     ! Generate initial RHS------------------------------------------------------80

     CALL VecZeroEntries(Vec_F,ierr)
     CALL ApplyGravity
     CALL ApplySource
     CALL FormRHS_Debug
     CALL MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     CALL MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)

     ! Set up linear solver      
     CALL KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     CALL SetupKSPSolver
     CALL PrintMsg("Solving ...")
     CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     CALL GetVec_U; tot_uu=tot_uu+uu
     IF (nobs_loc>0) THEN
        CALL GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     END IF
  END IF




















  !===========================================================================80
  ! Implicit Solver
  IF (stype/="explicit" .AND. (.NOT. fault)) THEN
     CALL KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     CALL SetupKSPSolver
     CALL PrintMsg("Solving ...")
     CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     CALL GetVec_U; tot_uu=tot_uu+uu
     IF (nobs_loc>0) THEN
        CALL GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     END IF
     IF (visco) THEN
        ! Recover stress
        CALL PrintMsg("Recovering stress ...")
        DO i=1,nels
           CALL RecoverStress(i,stress)
        END DO
     END IF
     ! Write output
     CALL WriteOutput
     ! Prepare for time stepping
     IF (t>f0 .AND. dt>f0 .AND. t>=dt) THEN
        IF (poro) THEN
           ! Form Kc and Up
           CALL VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
           CALL VecCopy(Vec_U,Vec_Um,ierr)
           CALL VecGetLocalSize(Vec_U,j,ierr)
           CALL VecGetOwnershipRange(Vec_U,j1,j2,ierr)
           j2=0
           DO i=1,j
              IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j2=j2+1
              END IF
           END DO
           ALLOCATE(work(j2))
           j2=0
           DO i=1,j
              IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j2=j2+1
                 work(j2)=j1+i-1
              END IF
           END DO
           j=SIZE(work)
           CALL ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,     &
                RI,ierr)
           CALL MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
           CALL MatZeroEntries(Mat_Kc,ierr)
           ALLOCATE(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
           DO i=1,nels
              CALL FormLocalK(i,k,indx,"Kc")
              kc=k(eldof+1:,eldof+1:)
              indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
              indxp=((indxp+1)/(dmn+1))-1
              CALL MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,           &
                   Add_Values,ierr)
           END DO
           CALL MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
           CALL MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
           CALL VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
           CALL VecDuplicate(Vec_Up,Vec_I,ierr) ! I->KcUp
           CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
           ALLOCATE(uup(SIZE(work)))
        END IF
        steps=INT(CEILING(t/dt))
        ! Start time stepping
        DO tstep=1,steps
           IF (rank==0) PRINT'(A11,I0)'," Time Step ",tstep
           ! Reform stiffness matrix, if needed
           IF (visco .AND. (tstep==1 .OR. MAXVAL(mat(:,4))>f1)) THEN
              CALL PrintMsg(" Reforming [K] ...")
              CALL MatZeroEntries(Mat_K,ierr)
              DO i=1,nels
                 CALL FormLocalK(i,k,indx,"Kv")
                 indx=indxmap(indx,2)
                 CALL MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,        &
                      Add_Values,ierr)
              END DO
              ! Account for constraint eqn's
              IF (rank==0 .AND. nceqs>0) CALL ApplyConstraints
              CALL MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
              CALL MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
           END IF
           ! Reform RHS
           CALL VecZeroEntries(Vec_F,ierr)
           CALL PrintMsg(" Reforming RHS ...")
           IF (poro) THEN
              CALL VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
              CALL MatMult(Mat_Kc,Vec_Up,Vec_I,ierr)
              CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              CALL VecScale(Vec_I,-dt,ierr)
              CALL VecGetArrayF90(Vec_I,pntr,ierr)
              uup=pntr
              CALL VecRestoreArrayF90(Vec_I,pntr,ierr)
              j=SIZE(uup)
              CALL VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              ! Stabilize RHS
              DO i=1,nels
                 CALL FormLocalHs(i,Hs,indxp)
                 indxp=indxmap(indxp,2)
                 CALL VecSetValues(Vec_F,npel,indxp,Hs,Add_Values,ierr)
              END DO
           END IF
           IF (visco) THEN
              DO i=1,nels
                 CALL ReformLocalRHS(i,f,indx)
                 indx=indxmap(indx,2)
                 CALL VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
              END DO
           END IF
           CALL FormRHS
           CALL VecAssemblyBegin(Vec_F,ierr)
           CALL VecAssemblyEnd(Vec_F,ierr)
           ! Solve
           CALL PrintMsg(" Solving ...")
           CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
           CALL GetVec_U; tot_uu=tot_uu+uu
           IF (poro) CALL VecAXPY(Vec_Um,f1,Vec_U,ierr)
           IF (visco) THEN
              ! Recover stress
              CALL PrintMsg(" Recovering stress ...")
              DO i=1,nels
                 CALL RecoverVStress(i,stress)
              END DO
           END IF
           ! Write output
           IF (MOD(tstep,frq)==0) CALL WriteOutput
           IF (nobs_loc>0) THEN
              CALL GetVec_obs
              tot_uu_obs=tot_uu_obs+uu_obs
              CALL WriteOutput_obs
           END IF
           !end if
        END DO
        IF (poro) DEALLOCATE(uup,kc,indxp,Hs,work)
     END IF
     IF (poro) THEN
        CALL VecDestroy(Vec_I,ierr)
        CALL VecDestroy(Vec_Up,ierr)
        CALL MatDestroy(Mat_Kc,ierr)
        CALL ISDestroy(RI,ierr)
        CALL VecDestroy(Vec_Um,ierr)
     END IF
     CALL KSPDestroy(Krylov,ierr)
  END IF
  !---------------------------------------------------------------------------80

  !===========================================================================80
  ! Fault/hybrid solver
  IF (fault) THEN
     ! Local to global fault node map
     IF (nceqs-nceqs_ncf>0) CALL GetFltMap
     IF (bod_frc==1) THEN
        ! Apply gravitational body force
        CALL PrintMsg("Applying gravity ...")
        CALL ApplyGravity
     END IF
     CALL KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     CALL KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     CALL SetupKSPSolver
     CALL PrintMsg("Solving ...")
     CALL VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     ! Debug Nodal Ownership Printout
     print'(A,I0,A,I0,A,I0)',"  Rank ",rank," has dofs ",j1+1," to ",j2
     IF (rank==nprcs-1) PRINT'(I0,A,I0,A)',j2+nceqs," dofs on ", nprcs,        &
          " processors."
     ! Solve system to generate initial U values           
     CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     CALL GetVec_U; tot_uu=tot_uu+uu
     ! Get observation
     IF (nobs_loc>0) THEN
        CALL GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     END IF
     IF (visco .OR. lm_str==1) THEN
        ! Recover stress
        CALL PrintMsg("Recovering stress ...")
        DO i=1,nels
           CALL RecoverStress(i,stress)
        END DO
     END IF
     IF (t>f0 .AND. dt>f0 .AND. t>=dt) THEN

        CALL VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
        CALL VecCopy(Vec_U,Vec_Um,ierr)
        CALL VecGetLocalSize(Vec_U,j,ierr)
        IF (poro) THEN
           ! Create Index Set Scheme for displacement, pressure, and LM locations from Vec_U
           j2=0; j3=0; j4=0; j5=0; j6=0
           DO i=1,j
              IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j2=j2+1
              END IF
              IF (nceqs>0) THEN
                 IF (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) THEN
                    j3=j3+1
                 END IF
              END IF
              IF (MOD(j1+i,dmn+1)>0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j4=j4+1
              END IF
              IF (j1+i-1<(dmn+1)*nnds) THEN
                 j5=j5+1
              END IF
           END DO
           ALLOCATE(work(j2),workl(j3),worku(j4))
           j2=0; j3=0; j4=0; j5=0
           DO i=1,j
              IF (MOD(j1+i,dmn+1)==0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j2=j2+1
                 work(j2)=j1+i-1
              END IF
              IF (MOD(j1+i,dmn+1)>0 .AND. j1+i-1<(dmn+1)*nnds) THEN
                 j4=j4+1
                 worku(j4)=j1+i-1
              END IF
              IF (nceqs>0) THEN
                 IF (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) THEN
                    j3=j3+1
                    workl(j3)=j1+i-1
                 END IF
              END IF
           END DO
           j=SIZE(work)
           CALL ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,RI,  &
                ierr)
           j=SIZE(worku)
           CALL ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,    &
                RIu,ierr)
           CALL MatGetSubMatrix(Mat_K,RIu,RI,Mat_Initial_Matrix,Mat_H,ierr)
           IF (nceqs > 0) THEN
              j=SIZE(workl)
              CALL ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values, &
                   RIl,ierr)
           END IF
           ! END ABORTION OF INDEX SET GENERATION-----------------------------80

           ! Generate Compressibility Submatrix
           CALL MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
           CALL MatZeroEntries(Mat_Kc,ierr)
           ALLOCATE(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
           DO i=1,nels
              CALL FormLocalK(i,k,indx,"Kc")
              kc=k(eldof+1:,eldof+1:)
              indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
              indxp=((indxp+1)/(dmn+1))-1
              CALL MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,           &
                   Add_Values,ierr)
           END DO
           CALL MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
           CALL MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)

           ! Create Initial Pressure Vector
           CALL VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
           IF (init==1) CALL VecDuplicate(Vec_Up,Vec_Up0,ierr) ! Initial p
           CALL VecDuplicate(Vec_Up,Vec_I,ierr) ! I->KcUp
           CALL VecCopy(Vec_Up,Vec_I,ierr) ! Store Up

           ! Create deformation induced fluid source vector
           CALL VecDuplicate(Vec_Up,Vec_qu,ierr) ! qu->Htu
           ! Create constraint induced fluid source vector
           CALL VecDuplicate(Vec_Up,Vec_ql,ierr)
           CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
           ALLOCATE(uup(SIZE(work)))

           ! Initialize space for lambda, p related nodal force
           CALL VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           CALL VecDuplicate(Vec_Uu,Vec_fp,ierr) ! fp->Hp
           CALL VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
           CALL VecCopy(Vec_Uu,Vec_fl,ierr) ! Hold Uu
           CALL VecDuplicate(Vec_Uu,Vec_flc,ierr)
           IF (lm_str==1) THEN
              CALL VecDuplicate(Vec_Uu,Vec_SS,ierr)
              CALL VecDuplicate(Vec_Uu,Vec_SH,ierr)
              CALL VecZeroEntries(Vec_SS,ierr)
              CALL VecZeroEntries(Vec_SH,ierr)
              IF (nceqs>0) THEN
                 CALL VecDuplicate(Vec_Uu,Vec_f2s,ierr)
                 CALL VecZeroEntries(Vec_f2s,ierr)
              END IF
              CALL VecDuplicate(Vec_Uu,Vec_dip,ierr)
              CALL VecZeroEntries(Vec_dip,ierr)
              CALL VecDuplicate(Vec_Uu,Vec_nrm,ierr)
              CALL VecZeroEntries(Vec_nrm,ierr)
           END IF
           CALL VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           j=SIZE(indxmap_u,1)
           IF (nceqs>0) THEN
              CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
              CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
           END IF
           CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_fp,ierr)
           CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),              &
                Petsc_Copy_Values,From_u,ierr)
           CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),              &
                Petsc_Copy_Values,To_u,ierr)
           CALL VecScatterCreate(Vec_fp,From_u,Seq_fp,To_u,Scatter_u,ierr)
           j=SIZE(nl2g,1)
           IF (nceqs>0) CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_ql,ierr)
           CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_qu,ierr)
           CALL ISCreateGeneral(Petsc_Comm_Self,j,nl2g(:,2)-1,                 &
                Petsc_Copy_Values,From_p,ierr)
           CALL ISCreateGeneral(Petsc_Comm_Self,j,nl2g(:,1)-1,                 &
                Petsc_Copy_Values,To_p,ierr)
           CALL VecScatterCreate(Vec_qu,From_p,Seq_qu,To_p,Scatter_q,ierr)
           ALLOCATE(fp(SIZE(indxmap_u,1)))
           ALLOCATE(qu(SIZE(nl2g,1)))
           IF (nceqs>0) THEN
              ALLOCATE(fl(SIZE(indxmap_u,1)))
              ALLOCATE(ql(SIZE(nl2g,1)))
              ALLOCATE(flc(SIZE(indxmap_u,1)))
           END IF
           IF (vout==1) THEN ! Extract nodal force by p, and fluid source by u
              CALL MatMult(Mat_H,Vec_I,Vec_fp,ierr)
              CALL VecScale(Vec_fp,-f1,ierr)
              CALL GetVec_fp
              CALL MatCreateTranspose(Mat_H,Mat_Ht,ierr)
              CALL MatMult(Mat_Ht,Vec_fl,Vec_qu,ierr)
              CALL GetVec_qu
           END IF
           ! Extract lambda
           IF (nceqs>0) THEN
              ! Vector to communicate with dynamic LM
              CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              !call LM_s2d 
              ! Extract lambda induced nodal force
              IF (vout==1) THEN
                 CALL VecZeroEntries(Vec_fl,ierr)
                 CALL VecZeroEntries(Vec_flc,ierr)
                 CALL VecZeroEntries(Vec_ql,ierr)
                 CALL GetVec_ql
                 CALL GetVec_flambda
                 CALL GetVec_fl
                 CALL GetVec_fcoulomb
                 CALL GetVec_flc
              ELSE
                 ! Extract Coulomb force (needed by force-stress translation).  
                 CALL GetVec_fcoulomb
              END IF
              CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
           END IF
           IF (lm_str==1) THEN ! Scatter stress to nodes
              CALL GetVec_Stress
              CALL GetVec_S
              IF (nceqs>0) THEN
                 CALL GetVec_f2s
                 CALL GetVec_f2s_seq
              ELSE
                 CALL GetVec_ft
              END IF
              CALL GetVec_dip_nrm
              IF (vout==1) CALL WriteOutput_f2s
              CALL VecDestroy(Vec_dip,ierr)
              CALL VecDestroy(Seq_dip,ierr)
              CALL VecDestroy(Vec_nrm,ierr)
              CALL VecDestroy(Seq_nrm,ierr)
              CALL VecDestroy(Seq_f2s,ierr)
              DEALLOCATE(f2s,dip,nrm)
           END IF
           ! Write output
           IF (vout==1) CALL WriteOutput_x
           IF (nobs_loc>0) CALL WriteOutput_obs
           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++80
           ! Case of initial pressure specification
           IF (init==1) THEN
              CALL PrintMsg("Pore fluid initialization...")
              CALL VecGetSubVector(Vec_Um,RI,Vec_Up0,ierr)
              ! Zero initial pressure 
              CALL VecZeroEntries(Vec_Up0,ierr)
              CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up0,ierr)
              CALL VecDestroy(Vec_Up0,ierr)
              tot_uu=f0
              CALL VecZeroEntries(Vec_F,ierr)
              !              CALL VecZeroEntries(Vec_U,ierr)
              !              CALL VecZeroEntries(Vec_Um,ierr)

              ! Apply initial force and pressure constraints
              CALL ApplySource
              !              CALL ApplyGravity
              !              CALL FormRHS
              CALL VecAssemblyBegin(Vec_F,ierr)
              CALL VecAssemblyEnd(Vec_F,ierr)           
              CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
              CALL GetVec_U; tot_uu=tot_uu+uu
              CALL VecAXPY(Vec_Um,f1,Vec_U,ierr)
              IF (vout==1) THEN
                 IF (nceqs>0) THEN
                    CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    CALL VecZeroEntries(Vec_flc,ierr)
                    CALL GetVec_fcoulomb
                    CALL GetVec_flc
                    CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 END IF
                 ! Initial state should be analyzed to see if any initial slip 
                 CALL WriteOutput_init
              END IF
           END IF
           IF (nceqs-nceqs_ncf>0) THEN
              ALLOCATE(flt_ss(nfnd,dmn),flt_p(nfnd))
              CALL GetVec_flt_qs 
              IF (rank==0) CALL WriteOutput_flt_qs
           END IF
           CALL Rscdt(fdt)
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~80           
        ELSE ! Not poro
           CALL VecGetLocalSize(Vec_U,j,ierr)
           CALL VecGetOwnershipRange(Vec_U,j1,j2,ierr)
           j3=0; j4=0
           DO i=1,j
              IF (j1+i-1>=dmn*nnds+nceqs_ncf) THEN
                 j3=j3+1
              ELSEIF (j1+i-1<dmn*nnds) THEN
                 j4=j4+1
              END IF
           END DO
           ALLOCATE(workl(j3)); ALLOCATE(worku(j4))
           j3=0; j4=0
           DO i=1,j
              IF (j1+i-1>=dmn*nnds+nceqs_ncf) THEN
                 j3=j3+1
                 workl(j3)=j1+i-1
              ELSEIF (j1+i-1<dmn*nnds) THEN
                 j4=j4+1
                 worku(j4)=j1+i-1
              END IF
           END DO
           j=SIZE(worku)
           CALL ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,    &
                RIu,ierr)
           j=SIZE(workl)
           CALL ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values,    &
                RIl,ierr)
           CALL VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           CALL VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
           IF (nceqs>0) CALL VecDuplicate(Vec_Uu,Vec_flc,ierr)
           IF (lm_str==1) THEN
              CALL VecDuplicate(Vec_Uu,Vec_SS,ierr)
              CALL VecDuplicate(Vec_Uu,Vec_SH,ierr)
              CALL VecZeroEntries(Vec_SS,ierr)
              CALL VecZeroEntries(Vec_SH,ierr)
              IF (nceqs>0) THEN
                 CALL VecDuplicate(Vec_Uu,Vec_f2s,ierr)
                 CALL VecZeroEntries(Vec_f2s,ierr)
              END IF
              CALL VecDuplicate(Vec_Uu,Vec_dip,ierr)
              CALL VecZeroEntries(Vec_dip,ierr)
              CALL VecDuplicate(Vec_Uu,Vec_nrm,ierr)
              CALL VecZeroEntries(Vec_nrm,ierr)
           END IF
           CALL VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           j=SIZE(indxmap_u,1)
           IF (nceqs>0) THEN
              CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
              CALL VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
              CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),           &
                   Petsc_Copy_Values,From_u,ierr)
              CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),           &
                   Petsc_Copy_Values,To_u,ierr)
              CALL VecScatterCreate(Vec_fl,From_u,Seq_fl,To_u,Scatter_u,ierr)
           ELSEIF (lm_str==1) THEN
              CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),           &
                   Petsc_Copy_Values,From_u,ierr)
              CALL ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),           &
                   Petsc_Copy_Values,To_u,ierr)
              CALL VecScatterCreate(Vec_U,From_u,Seq_U,To_u,Scatter_u,ierr)
           END IF
           IF (nceqs>0) THEN
              ALLOCATE(fl(SIZE(indxmap_u,1)))
              ALLOCATE(flc(SIZE(indxmap_u,1)))
              ! Vector to communicate with dynamic LM
              CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              !call LM_s2d
              ! Extract lambda induced nodal force
              IF (vout==1) THEN
                 CALL VecZeroEntries(Vec_fl,ierr)
                 CALL VecZeroEntries(Vec_flc,ierr)
                 CALL GetVec_flambda
                 CALL GetVec_fl
                 CALL GetVec_fcoulomb
                 CALL GetVec_flc
              ELSE
                 ! Extract Coulomb force (needed by force-stress translation).
                 CALL GetVec_fcoulomb
              END IF
              CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
           END IF
           IF (lm_str==1) THEN ! Scatter stress to nodes
              CALL GetVec_Stress
              CALL GetVec_S
              IF (nceqs>0) THEN
                 CALL GetVec_f2s
                 CALL GetVec_f2s_seq
              ELSE
                 CALL GetVec_ft
              END IF
              CALL GetVec_dip_nrm
              IF (vout==1) CALL WriteOutput_f2s
              ! Cleanup
              CALL VecDestroy(Vec_dip,ierr)
              CALL VecDestroy(Seq_dip,ierr)
              CALL VecDestroy(Vec_nrm,ierr)
              CALL VecDestroy(Seq_nrm,ierr)
              CALL VecDestroy(Seq_f2s,ierr)
              DEALLOCATE(f2s,dip,nrm)
              IF (nceqs-nceqs_ncf>0) THEN
                 ALLOCATE(flt_ss(nfnd,dmn))
                 CALL GetVec_flt_qs
                 IF (rank==0) CALL WriteOutput_flt_qs
              END IF
           END IF
           ! Write output
           IF (vout==1) CALL WriteOutput_x
           !if (rank==0 .and. nobs>0) call WriteOutput_obs
           IF (nobs_loc>0) CALL WriteOutput_obs
        END IF ! Poro or not
        !=====================================================================80
        ! Solution space is allocated differently for static and dynamic runs, 
        ! so we keep mat_K and Mat_K_dyn separate instead of
        ! call MatGetSubMatrix(Mat_K,RIu,RIu,Mat_Initial_Matrix,Mat_K_dyn,ierr)

        ! Initialize slip indicator and slip history
        IF (nceqs>0 .AND. hyb>0) THEN
           ALLOCATE(slip(nfnd),slip0(nfnd),slip_sum(nfnd)) 
           slip(:)=0;slip_sum=0
           n_log_wave=0
           n_log_slip=0
           IF (rsf==1) ALLOCATE(mu_cap(nfnd),rsfv(nfnd),rsf_sta(nfnd_loc))
           trunc=f0
           ALLOCATE(mu_hyb(nfnd))
        END IF

        ! Start implicit time step
        steps=INT(CEILING(t/dt)); t_abs=f0
        dyn=.FALSE.; fail=.FALSE.; n_log=0
        DO tstep=1,steps
           t_abs=t_abs+dt
           IF (rank==0) PRINT'(A11,I0)'," Time Step ",tstep
           ! Reform stiffness matrix, if needed
           IF (visco .AND. (tstep==1 .OR. MAXVAL(mat(:,4))>f1)) THEN
              CALL PrintMsg(" Reforming [K] ...")
              CALL MatZeroEntries(Mat_K,ierr)
              DO i=1,nels
                 CALL FormLocalK(i,k,indx,"Kv")
                 indx=indxmap(indx,2)
                 CALL MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,        &
                      Add_Values,ierr)
              END DO
              ! Account for constraint eqn's
              IF (rank==0 .AND. nceqs>0) CALL ApplyConstraints
              CALL MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
              CALL MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
           END IF
           ! Reform RHS
           CALL VecZeroEntries(Vec_F,ierr)
           CALL PrintMsg(" Reforming RHS ...")
           IF (poro) THEN
              CALL VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
              CALL PrintMsg(" Updating f_p ...")
              CALL MatMult(Mat_Kc,Vec_Up,Vec_I,ierr)
              CALL VecScale(Vec_I,-dt,ierr)
              CALL VecGetArrayF90(Vec_I,pntr,ierr)
              uup=pntr
              CALL VecRestoreArrayF90(Vec_I,pntr,ierr)
              j=SIZE(uup)
              CALL VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              ! Stabilize RHS
              DO i=1,nels
                 CALL FormLocalHs(i,Hs,indxp)
                 indxp=indxmap(indxp,2)
                 CALL VecSetValues(Vec_F,npel,indxp,Hs,Add_Values,ierr)
              END DO
           END IF
           IF (visco) THEN
              DO i=1,nels
                 CALL ReformLocalRHS(i,f,indx)
                 indx=indxmap(indx,2)
                 CALL VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
              END DO
           END IF
           CALL FormRHS
           IF (fail) THEN 
              CALL FaultSlip
              ! Backup QS slip 
              qs_flt_slip=qs_flt_slip+tot_flt_slip
              fail=.FALSE.
           END IF
           ! Solve
           CALL VecAssemblyBegin(Vec_F,ierr)
           CALL VecAssemblyEnd(Vec_F,ierr)
           CALL PrintMsg(" Solving ...")
           CALL KSPSolve(Krylov,Vec_F,Vec_U,ierr)
           ! Reset dynamic (slip) solutions 
           IF (nceqs>0 .AND. hyb>0) THEN
              tot_uu_dyn=f0
              IF (nobs_loc>0) tot_uu_dyn_obs=f0
              flt_slip=f0; tot_flt_slip=f0
           END IF
           n_log=n_log+1
           CALL GetVec_U; tot_uu=tot_uu+uu
           IF (nobs_loc>0) THEN 
              CALL GetVec_obs
              tot_uu_obs=tot_uu_obs+uu_obs
           END IF
           ! Record absolute solution
           IF (poro .OR. nceqs>0) CALL VecAXPY(Vec_Um,f1,Vec_U,ierr)
           IF (visco .OR. lm_str==1) THEN
              ! Recover stress
              CALL PrintMsg(" Recovering stress ...")
              DO i=1,nels
                 CALL RecoverVStress(i,stress)
              END DO
              CALL GetVec_Stress
              CALL GetVec_S
              IF (fault) THEN
                 CALL GetVec_flt_qs
                 IF (rank==0) CALL WriteOutput_flt_qs
              END IF
           END IF
           ! Extract nodal force by p
           IF (poro) THEN
              IF (vout==1) THEN
                 ! Extract force by p
                 CALL VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
                 CALL MatMult(Mat_H,Vec_Up,Vec_fp,ierr)
                 CALL VecScale(Vec_fp,-f1,ierr)
                 CALL GetVec_fp
                 CALL VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
                 ! Extract fluid source by u
                 CALL VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
                 CALL MatMult(Mat_Ht,Vec_Uu,Vec_qu,ierr)
                 CALL GetVec_qu
                 CALL VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
                 ! Extract lambda and induced nodal f and q
                 IF (nceqs>0) THEN
                    CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    CALL VecZeroEntries(Vec_fl,ierr)
                    CALL VecZeroEntries(Vec_flc,ierr)
                    CALL VecZeroEntries(Vec_ql,ierr)
                    CALL GetVec_flambda
                    CALL GetVec_fl
                    CALL GetVec_ql
                    CALL GetVec_fcoulomb
                    CALL GetVec_flc
                    CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 END IF
                 ! Write output
                 IF (MOD(tstep,frq)==0) CALL WriteOutput_x
              END IF
              !if (rank==0 .and. nobs>0) call WriteOutput_obs
              IF (nobs_loc>0) CALL WriteOutput_obs
           ELSE
              IF (vout==1) THEN
                 ! Extract lambda and induced nodal f
                 IF (nceqs>0) THEN
                    CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    CALL VecZeroEntries(Vec_fl,ierr)
                    CALL VecZeroEntries(Vec_flc,ierr)
                    CALL GetVec_flambda
                    CALL GetVec_fl
                    CALL GetVec_fcoulomb
                    CALL GetVec_flc
                    CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 END IF
                 ! Write output
                 IF (MOD(tstep,frq)==0) CALL WriteOutput_x
              END IF
              !if (rank==0 .and. nobs>0) call WriteOutput_obs
              IF (nobs_loc>0) CALL WriteOutput_obs
           END IF

           ! Determine if the fault will fail by LM 
           IF (nceqs-nceqs_ncf>0 .AND. hyb>0) THEN
              CALL VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              CALL LM_s2d  
              CALL VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              ! Determine if the fault will fail
              CALL GetSlip_sta
              rslip=REAL(SUM(slip))/REAL(SIZE(slip))
              IF (rank==0) PRINT'(F0.2,A)',rslip*100.0,"% fault critical."
              IF (rslip>f0) THEN ! Failure threshold
                 dyn=.TRUE. 
              END IF
           END IF

           ! Prepare variables for hybrid run
           IF (dyn) THEN
              ! Prepare the working space for the dynamic model
              CALL VecZeroEntries(Vec_U_dyn,ierr)
              CALL VecDuplicate(Vec_U_dyn,Vec_Um_dyn,ierr)
              CALL VecDuplicate(Vec_U_dyn,Vec_Up_dyn,ierr)
              CALL VecDuplicate(Vec_U_dyn,Vec_U_dyn_tot,ierr)
              CALL VecZeroEntries(Vec_Um_dyn,ierr)
              CALL VecZeroEntries(Vec_Up_dyn,ierr)
              CALL VecZeroEntries(Vec_U_dyn_tot,ierr)
              CALL VecDuplicateVecsF90(Vec_U_dyn,6,Vec_W,ierr)
              IF (n_log_dyn==0) THEN! Create full Gt
                 CALL GetMat_Gt
                 CALL MatAssemblyBegin(Mat_Gt,Mat_Final_Assembly,ierr)
                 CALL MatAssemblyEnd(Mat_Gt,Mat_Final_Assembly,ierr)
              END IF
              CALL MatCreateTranspose(Mat_Gt,Mat_G,ierr)
              CALL VecDuplicate(Vec_lambda_sta,Vec_lambda,ierr)
              ! GMinvGt is replaced by its inverse, assuming GMinvGt is diagonal
              CALL MatPtAP(Mat_Minv,Mat_Gt,Mat_Initial_Matrix,f1,              &
                   Mat_GMinvGt,ierr)
              CALL VecDuplicateVecsF90(Vec_lambda,2,Vec_Wlm,ierr)
              CALL MatGetDiagonal(Mat_GMinvGt,Vec_Wlm(1),ierr)
              CALL MatSetOption(Mat_GMinvGt,Mat_New_Nonzero_Allocation_Err,    &
                   PETSC_FALSE,ierr)
              CALL VecReciprocal(Vec_Wlm(1),ierr)
              CALL MatDiagonalSet(Mat_GMinvGt,Vec_Wlm(1),Insert_Values,ierr)
              CALL VecZeroEntries(Vec_Wlm(1),ierr)
              ! Dynamic constraint I=0 
              CALL VecDuplicate(Vec_lambda,Vec_I_dyn,ierr)
              CALL VecZeroEntries(Vec_I_dyn,ierr)
              ! Pass pseudo velocity to U_dyn
              IF (rsf==1 .AND. tstep>1) THEN
                 CALL Rsfv2Dyn 
                 rsf_sta=0
              END IF
              CALL VecDuplicate(Vec_lambda,Vec_lambda_tot,ierr)
              CALL VecZeroEntries(Vec_lambda_tot,ierr)
              ! Form 1/dt^2GMinvGt assuming it doesnt change with time
              CALL MatScale(Mat_GMinvGt,f1/dt_dyn**2,ierr)
              CALL PrintMsg("Hybrid Solving ...")
              steps_dyn=INT(CEILING(t_dyn/dt_dyn)); t_hyb=t_abs
              ih=0; t_sta=f0; fail=.TRUE.;crp=.FALSE.
           END IF

           ! Explicit/implicit hybrid step for rupture propagation
           DO WHILE (dyn .OR. (t_sta>f0 .AND. (t_sta<t_lim) .AND. .NOT. crp))
              ! Explicit time step
              DO tstep_dyn=0,steps_dyn
                 t_hyb=t_hyb+dt_dyn
                 ! Up=Minv(dt^2(F-KU)-dt(C(U-Um)))+2U-Um
                 CALL MatMult(Mat_K_dyn,Vec_U_dyn,Vec_W(1),ierr)
                 CALL VecScale(Vec_W(1),-f1*dt_dyn**2,ierr)
                 CALL VecWAXPY(Vec_W(2),-f1,Vec_Um_dyn,Vec_U_dyn,ierr)
                 CALL MatMult(Mat_K_dyn,Vec_W(2),Vec_W(3),ierr)
                 CALL MatMult(Mat_M,Vec_W(2),Vec_W(4),ierr)
                 CALL VecAXPY(Vec_W(4),beta,Vec_W(3),ierr)
                 CALL VecScale(Vec_W(4),dt_dyn,ierr)
                 CALL VecWAXPY(Vec_W(5),-f1,Vec_W(4),Vec_W(1),ierr)
                 CALL MatMult(Mat_Minv,Vec_W(5),Vec_W(6),ierr)
                 CALL VecWAXPY(Vec_Up_dyn,f2,Vec_U_dyn,Vec_W(6),ierr)
                 CALL VecAXPY(Vec_Up_dyn,-f1,Vec_Um_dyn,ierr)
                 ! Form lambda=(GUp-Flm)/(dt^2GMinvGt)
                 CALL MatMult(Mat_G,Vec_Up_dyn,Vec_Wlm(1),ierr)
                 CALL VecWAXPY(Vec_Wlm(2),-f1,Vec_I_dyn,Vec_Wlm(1),ierr)
                 CALL MatMult(Mat_GMinvGt,Vec_Wlm(2),Vec_lambda,ierr)
                 IF (rsf==1 .AND. tstep>1 .AND. tstep_dyn+ih==0) THEN 
                    ! Skip friction law, and reset constraint 
                    CALL VecZeroEntries(Vec_I_dyn,ierr)
                 ELSE 
                    ! Cap the nodal LM not to exceed max friction
                    CALL CapLM_dyn
                 END IF
                 CALL VecAXPY(Vec_lambda_tot,f1,Vec_lambda,ierr)
                 ! Form Up=Up-dt^2(Minv(Gtlambda))
                 CALL MatMult(Mat_Gt,Vec_lambda,Vec_W(1),ierr)
                 CALL MatMult(Mat_Minv,Vec_W(1),Vec_W(2),ierr)
                 CALL VecAXPY(Vec_Up_dyn,-dt_dyn**2,Vec_W(2),ierr)
                 ! Update solution
                 CALL VecCopy(Vec_U_dyn,Vec_Um_dyn,ierr)
                 CALL VecCopy(Vec_Up_dyn,Vec_U_dyn,ierr)
                 CALL VecAXPY(Vec_U_dyn_tot,f1,Vec_U_dyn,ierr)
                 ! Extract dynamic solution
                 IF (.NOT. dyn) THEN
                    dyn=.TRUE.
                    CALL GetVec_U
                    tot_uu_dyn=tot_uu_dyn+uu_dyn
                    IF (nobs_loc>0) THEN
                       CALL GetVec_obs
                       tot_uu_dyn_obs=tot_uu_dyn_obs+uu_dyn_obs
                       uu_dyn_obs=uu_dyn_obs/dt_dyn
                       IF (MOD(n_log_dyn,frq_wave)==0) CALL WriteOutput_obs 
                    END IF
                    IF (MOD(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
                    dyn=.FALSE.
                 ELSE
                    CALL GetVec_U
                    tot_uu_dyn=tot_uu_dyn+uu_dyn
                    IF (nobs_loc>0) THEN
                       CALL GetVec_obs
                       tot_uu_dyn_obs=tot_uu_dyn_obs+uu_dyn_obs
                       uu_dyn_obs=uu_dyn_obs/dt_dyn
                       IF (MOD(n_log_dyn,frq_wave)==0) CALL WriteOutput_obs
                    END IF
                    IF (MOD(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
                 END IF
                 ! Evaluated FD grid movement
                 IF (ngp_loc>0 .AND. MOD(n_log_dyn,frq_wave)==0 .AND. fdout==1) THEN  
                    CALL GetVec_fd
                    CALL WriteOutput_fd
                 END IF
                 ! Extract and output temporal fault slip
                 CALL MatMult(Mat_G,Vec_U_dyn,Vec_Wlm(1),ierr)
                 CALL VecGetArrayF90(Vec_Wlm(1),pntr,ierr)
                 flt_slip=pntr
                 tot_flt_slip=tot_flt_slip+flt_slip
                 CALL VecRestoreArrayF90(Vec_Wlm(1),pntr,ierr)
                 IF (MOD(n_log_dyn,frq_slip)==0) THEN 
                    CALL WriteOutput_slip
                    n_log_slip=n_log_slip+1
                 END IF
                 ! Export dynamic snapshot
                 dsp_dyn=.TRUE.
                 IF (vout==1) THEN
                    uu_dyn=uu_dyn/dt_dyn
                    IF (MOD(n_log_dyn,frq_dyn)==0) CALL WriteOutput
                 END IF
                 dsp_dyn=.FALSE.
                 n_log_dyn=n_log_dyn+1
              END DO ! Explicit loop 
              ! Assess the fault status 
              IF (dyn) CALL GetSlip_dyn
              rslip=REAL(SUM(slip))/REAL(SIZE(slip))
              IF ((.NOT. SUM(slip)>0) .OR. (.NOT. dyn)) THEN
                 dyn=.FALSE. ! Failure threshold
                 t_sta=t_sta+t_dyn
              END IF
              IF (rank==0 .AND. dyn) PRINT'(F0.2,A)',rslip*100.0,              &
                   "% fault slipping."
              IF (rank==0 .AND. crp) PRINT'(A)',"Aseismic"
              IF (rank==0 .AND. .NOT. (dyn .OR. crp)) PRINT'(A,F0.4,A,F0.4)',  &
                   "Wave traveling ",t_sta,"/",t_lim
              ! Turn off "dyn" for debug if the fault cannot stabilize
              IF (DBLE(ih+1)*t_dyn>=t_lim .AND. dyn) THEN
                 dyn=.FALSE.
                 slip=0
                 t_sta=t_sta+t_dyn
              END IF
              ! Hybrid iteration count
              ih=ih+1
           END DO ! Hybrid run
           IF (fail) THEN
              IF (rank==0 .AND. nobs>0) THEN 
                 CALL WriteOutput_log
                 CALL WriteOutput_log_wave
                 IF (nceqs>0) CALL WriteOutput_log_slip
              END IF
              ! Latest fault stress (Vec_lambda_sta0)
              CALL GetVec_lambda_hyb(trunc) ! Last time truncation
              ! Cleanup dynamics
              CALL VecDestroy(Vec_Um_dyn,ierr)
              CALL VecDestroy(Vec_Up_dyn,ierr)
              CALL VecDestroy(Vec_U_dyn_tot,ierr)
              CALL VecDestroy(Vec_I_dyn,ierr)
              CALL VecDestroyVecsF90(6,Vec_W,ierr)
              CALL VecDestroy(Vec_lambda,ierr)
              CALL VecDestroy(Vec_lambda_tot,ierr)
              CALL VecDestroyVecsF90(2,Vec_Wlm,ierr)
              CALL MatDestroy(Mat_G,ierr)
              CALL MatDestroy(Mat_GMinvGt,ierr)
           ELSE
              ! Backup fault stress
              CALL VecCopy(Vec_lambda_sta,Vec_lambda_sta0,ierr)
           END IF
        END DO ! Implicit time steps
     ELSEIF (.NOT. poro) THEN ! One time static
        CALL VecDuplicate(Vec_U,Vec_SS,ierr)
        CALL VecDuplicate(Vec_U,Vec_SH,ierr)
        CALL VecZeroEntries(Vec_SS,ierr)
        CALL VecZeroEntries(Vec_SH,ierr)
        CALL GetVec_Stress
        Scatter_u=Scatter
        CALL GetVec_S
        IF (vout==1) CALL WriteOutput_x
        IF (nobs_loc>0) CALL WriteOutput_obs
     END IF ! Assert implicit time
     crp=.TRUE.
     IF (rank==0) CALL WriteOutput_log 
     ! Cleanup hybrid
     IF (poro) DEALLOCATE(uup,kc,indxp,Hs,work,fp,qu)
     IF (nceqs>0) THEN 
        DEALLOCATE(fl,flc,workl,worku)
        IF (lm_str==1) DEALLOCATE(ss,sh)
        IF (poro) DEALLOCATE(ql)
     END IF
     CALL MatDestroy(Mat_Kc,ierr)
     CALL MatDestroy(Mat_H,ierr)
     CALL MatDestroy(Mat_Ht,ierr)
     CALL MatDestroy(Mat_K_dyn,ierr)
     CALL MatDestroy(Mat_Gt,ierr)
     CALL MatDestroy(Mat_M,ierr)
     CALL MatDestroy(Mat_Minv,ierr)
     CALL VecDestroy(Vec_I,ierr)
     CALL VecDestroy(Vec_fp,ierr)
     CALL VecDestroy(Vec_Up,ierr)
     CALL VecDestroy(Vec_Uu,ierr)
     CALL VecDestroy(Vec_Um,ierr)
     CALL VecDestroy(Vec_qu,ierr)
     CALL VecDestroy(Vec_Ul,ierr)
     CALL VecDestroy(Vec_lambda_sta,ierr)
     CALL VecDestroy(Vec_lambda_sta0,ierr)
     CALL VecDestroy(Vec_fl,ierr)
     CALL VecDestroy(Vec_flc,ierr)
     CALL VecDestroy(Vec_f2s,ierr)
     CALL VecDestroy(Vec_ql,ierr)
     CALL VecDestroy(Vec_U_dyn,ierr)
     CALL VecDestroy(Vec_SS,ierr)
     CALL VecDestroy(Vec_SH,ierr)
     CALL VecDestroy(Seq_fp,ierr)
     CALL VecDestroy(Seq_qu,ierr)
     CALL VecDestroy(Seq_fl,ierr)
     CALL VecDestroy(Seq_flc,ierr)
     CALL VecDestroy(Seq_ql,ierr)
     CALL VecDestroy(Seq_U_dyn,ierr)
     CALL VecDestroy(Seq_SS,ierr)
     CALL VecDestroy(Seq_SH,ierr)
     CALL ISDestroy(RI,ierr)
     CALL ISDestroy(RIu,ierr)
     CALL ISDestroy(RIl,ierr)
     CALL ISDestroy(From_u,ierr)
     CALL ISDestroy(To_u,ierr)   
     CALL ISDestroy(From_p,ierr)
     CALL ISDestroy(To_p,ierr)
     CALL VecScatterDestroy(Scatter_u,ierr)
     CALL VecScatterDestroy(Scatter_q,ierr)
     CALL VecScatterDestroy(Scatter_s2d,ierr)
     CALL VecScatterDestroy(Scatter_dyn,ierr)
     CALL KSPDestroy(Krylov,ierr)
  END IF ! Fault solver
  !---------------------------------------------------------------------------80

  !===========================================================================80
  ! Explicit (Green's function) Solver
  IF (stype=="explicit" .AND. t>f0 .AND. dt>f0 .AND. t>=dt) THEN
     steps=INT(CEILING(t/dt))
     ! Initialize work vectors
     CALL VecDuplicate(Vec_U,Vec_Um,ierr)
     CALL VecDuplicate(Vec_U,Vec_Up,ierr)
     CALL VecDuplicateVecsF90(Vec_U,6,Vec_W,ierr)
     IF (nceqs>0) THEN
        IF (gf) CALL GetMat_Gt
        CALL MatCreateTranspose(Mat_Gt,Mat_G,ierr)
        CALL VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nceqs,Vec_lambda,ierr)
        CALL VecDuplicate(Vec_lambda,Vec_I,ierr)
        CALL VecDuplicateVecsF90(Vec_lambda,2,Vec_Wlm,ierr)
     END IF
     ! Start time stepping
     CALL PrintMsg("Solving ...") ! Up=Minv(dt^2(F-KU)-dt(C(U-Um)))+2U-Um
     CALL VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     IF (rank==nprcs-1) PRINT'(I0,A,I0,A)',j2+nceqs," dofs on ", nprcs,        &
          " processors."
     ng=1
     IF (nobs_loc>0) tot_uu_dyn_obs=f0
     DO igf=1,nceqs
        IF (ng>10) EXIT ! Max mumber of GFs
        j=(igf-1)/dmn+1
        DO tstep=0,steps
           IF (gf .AND. frc(j)<1) THEN
              EXIT
           ELSEIF (gf .AND. tstep==0) THEN
              IF (MOD(igf,dmn)==0) THEN
                 EXIT
              ELSE
                 ng=ng+1
                 IF (rank==0) PRINT'(A,I0,A,I0)',"Source ",j," component ",     &
                      MOD(igf,dmn)
              END IF
           END IF
           IF (rank==0 .AND. .NOT. gf) PRINT'(A12,I0)'," Time Step ",tstep
           ! Solve
           CALL MatMult(Mat_K,Vec_U,Vec_W(1),ierr)
           CALL VecAYPX(Vec_W(1),-f1,Vec_F,ierr)
           CALL VecScale(Vec_W(1),dt**2,ierr)
           CALL VecWAXPY(Vec_W(2),-f1,Vec_Um,Vec_U,ierr)
           CALL MatMult(Mat_K,Vec_W(2),Vec_W(3),ierr)
           CALL MatMult(Mat_M,Vec_W(2),Vec_W(4),ierr)
           CALL VecAXPY(Vec_W(4),beta,Vec_W(3),ierr)
           CALL VecScale(Vec_W(4),dt,ierr)
           CALL VecWAXPY(Vec_W(5),-f1,Vec_W(4),Vec_W(1),ierr)
           CALL MatMult(Mat_Minv,Vec_W(5),Vec_W(6),ierr)
           CALL VecWAXPY(Vec_Up,f2,Vec_U,Vec_W(6),ierr)
           CALL VecAXPY(Vec_Up,-f1,Vec_Um,ierr)
           IF (nceqs>0) THEN
              IF (tstep==0) THEN
                 CALL MatPtAP(Mat_Minv,Mat_Gt,Mat_Initial_Matrix,f1,           &
                      Mat_GMinvGt,ierr)
                 ! GMinvGt is replaced by its inverse, assuming GMinvGt 
                 ! is diagonal
                 CALL MatGetDiagonal(Mat_GMinvGt,Vec_Wlm(1),ierr)
                 CALL VecReciprocal(Vec_Wlm(1),ierr)
                 CALL MatDiagonalSet(Mat_GMinvGt,Vec_Wlm(1),Insert_Values,ierr)
                 CALL VecZeroEntries(Vec_Wlm(1),ierr)
                 ! Form 1/dt^2GMinvGt assuming it doesnt change with time
                 CALL MatScale(Mat_GMinvGt,f1/dt**2,ierr)
              END IF
              ! Apply constraint values using Flm
              IF (gf) THEN
                 IF (tstep==0) THEN
                    CALL VecZeroEntries(Vec_I,ierr)
                    IF (rank==0) CALL VecSetValue(Vec_I,igf-1,f1,Insert_Values,&
                         ierr) 
                    CALL VecAssemblyBegin(Vec_I,ierr)
                    CALL VecAssemblyEnd(Vec_I,ierr)
                 ELSEIF (tstep==5) THEN
                    CALL VecZeroEntries(Vec_I,ierr)  
                 END IF
              ELSE
                 CALL VecZeroEntries(Vec_I,ierr)
                 IF (rank==0) THEN
                    DO i=1,nceqs
                       j1=NINT(cval(i,2)/dt); j2=NINT(cval(i,3)/dt)
                       IF (tstep>=j1 .AND. tstep<=j2) THEN
                          val=cval(i,1)
                          ! Use a decaying function instead of boxcar ...
                          !if (j2>j1) val=val*(dble(j2-tstep)/dble(j2-j1))**2.5
                          CALL VecSetValue(Vec_I,i-1,val,Insert_Values,ierr)
                       END IF
                    END DO
                 END IF
                 CALL VecAssemblyBegin(Vec_I,ierr)
                 CALL VecAssemblyEnd(Vec_I,ierr)
              END IF
              ! Form lambda=(GUp-Flm)/(dt^2GMinvGt)
              CALL MatMult(Mat_G,Vec_Up,Vec_Wlm(1),ierr)
              CALL VecWAXPY(Vec_Wlm(2),-f1,Vec_I,Vec_Wlm(1),ierr)
              CALL MatMult(Mat_GMinvGt,Vec_Wlm(2),Vec_lambda,ierr)
              ! Form Up=Up-dt^2(Minv(Gtlambda))
              CALL MatMult(Mat_Gt,Vec_lambda,Vec_W(1),ierr)
              CALL MatMult(Mat_Minv,Vec_W(1),Vec_W(2),ierr)
              CALL VecAXPY(Vec_Up,-dt**2,Vec_W(2),ierr)
           END IF
           CALL VecCopy(Vec_U,Vec_Um,ierr)
           CALL VecCopy(Vec_Up,Vec_U,ierr)
           IF (gf) THEN
              IF (nobs_loc>0) THEN
                 CALL GetVec_obs
                 tot_uu_dyn_obs=tot_uu_dyn_obs+uu_obs
                 uu_dyn_obs=uu_obs/dt
                 IF (MOD(n_log_dyn,frq_wave)==0) THEN
                    dyn=.TRUE.
                    CALL WriteOutput_obs;
                    dyn=.FALSE.
                 END IF
              END IF
              IF (MOD(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
              n_log_dyn=n_log_dyn+1
              ! Write output
              CALL GetVec_U; tot_uu=tot_uu+uu
              IF (MOD(tstep,frq)==0) CALL WriteOutput
           ELSE
              CALL VecZeroEntries(Vec_F,ierr)
              CALL FormRHS
              CALL VecAssemblyBegin(Vec_F,ierr)
              CALL VecAssemblyEnd(Vec_F,ierr)
              ! Write output
              CALL GetVec_U; tot_uu=tot_uu+uu
              IF (MOD(tstep,frq)==0) CALL WriteOutput
              IF (nceqs>0) CALL VecZeroEntries(Vec_I,ierr)
           END IF
        END DO ! Explicit run
        IF (.NOT. gf) EXIT 
        IF (rank==0 .AND. frc(j)==1 .AND. MOD(igf,dmn)/=0)                     &
             CALL WriteOutput_log_wave
        ! Reset solutions
        CALL VecZeroEntries(Vec_U,ierr) 
        CALL VecZeroEntries(Vec_Um,ierr)
        CALL VecZeroEntries(Vec_Up,ierr)
        IF (nobs_loc>0) tot_uu_dyn_obs=f0
     END DO ! Unit rupture loop
     CALL PrintMsg("Done")
     IF (nceqs>0) THEN
        CALL MatDestroy(Mat_Gt,ierr)
        CALL MatDestroy(Mat_G,ierr)
        CALL MatDestroy(Mat_GMinvGt,ierr)
        CALL VecDestroy(Vec_lambda,ierr)
        CALL VecDestroy(Vec_I,ierr)
        CALL VecDestroyVecsF90(2,Vec_Wlm,ierr)
     END IF
     CALL MatDestroy(Mat_M,ierr)
     CALL MatDestroy(Mat_Minv,ierr)
     CALL VecDestroy(Vec_Um,ierr)
     CALL VecDestroy(Vec_Up,ierr)
     CALL VecDestroyVecsF90(6,Vec_W,ierr)
  END IF ! Explicit solver
  !---------------------------------------------------------------------------80

  ! Cleanup
  CALL PrintMsg("Cleaning up ...")
  CALL VecScatterDestroy(Scatter,ierr)
  CALL ISDestroy(To,ierr)
  CALL ISDestroy(From,ierr)
  CALL VecDestroy(Seq_U,ierr)
  CALL VecDestroy(Vec_F,ierr)
  CALL VecDestroy(Vec_U,ierr)
  CALL MatDestroy(Mat_K,ierr)
  IF (visco) DEALLOCATE(stress)
  DEALLOCATE(coords,nodes,bc,mat,id,k,m,f,indx,ipoint,weight,enodes,ecoords,   &
       vvec,indxmap,tot_uu,uu,cval,fnode,fval,telsd,tval,nl2g)
  ! Delete cnstr.tmp 
  IF (rank==0 .AND. nceqs>0) THEN
     OPEN(15, file="cnstrns.tmp",status='old')
     CLOSE(15, status='delete')
  END IF
  CALL PrintMsg("Finished")

9 CALL PetscFinalize(ierr)
  !---------------------------------------------------------------------------80

  !===========================================================================80
CONTAINS

  ! Read simulation parameters
  SUBROUTINE ReadParameters
    IMPLICIT NONE
    READ(10,*)stype,eltype,nodal_bw
    fault=.FALSE.
    IF (stype=="fault" .OR. stype=="fault-p" .OR. stype=="fault-pv" .OR.       &
         stype=="fault-v") fault=.TRUE.
    gf=.FALSE.
    IF (stype=="explicit-gf") THEN
       stype="explicit"; gf=.TRUE.
    END IF
    IF (fault .OR. gf) THEN
       READ(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nfnd,nobs,nceqs_ncf,nrxs
    ELSE
       READ(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nobs
    END IF
    READ(10,*)t,dt,frq,dsp
    ! Dynamic run time before jumping back to static model
    IF (fault) THEN 
       READ(10,*)t_dyn,dt_dyn,frq_dyn,t_lim,dsp_hyb,lm_str,bod_frc,            &
            hyb,rsf,init
       IF (init==1) THEN
          fdt=dt/DBLE(3600*24)
          dt=DBLE(3600*24)
       END IF
    END IF
    IF (hyb==1 .AND. rsf==0) READ(10,*)frq_wave,frq_slip 
    IF (hyb==1 .AND. rsf==1) READ(10,*)frq_wave,frq_slip,v_bg
    poro=.FALSE.; visco=.FALSE.
    IF (stype=="implicit-p" .OR. stype=="implicit-pv") poro=.TRUE.
    IF (stype=="implicit-v" .OR. stype=="implicit-pv") visco=.TRUE.
    IF (dt==f0) dt=f1
    IF (stype=="explicit") READ(10,*)alpha,beta
    IF (stype=="fault-p" .OR. stype=="fault-pv") poro=.TRUE.
    IF (stype=="fault-v" .OR. stype=="fault-pv") visco=.TRUE.
    IF (fault) READ(10,*)alpha,beta
  END SUBROUTINE ReadParameters

  ! Setup implicit solver
  SUBROUTINE SetupKSPSolver
    IMPLICIT NONE
    CALL KSPSetType(Krylov,"gmres",ierr)
    CALL KSPGetPC(Krylov,PreCon,ierr)
    CALL PCSetType(PreCon,"asm",ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
    CALL KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Double_Precision,        &
         Petsc_Default_Double_Precision,Petsc_Default_Integer,ierr)
#else
    CALL KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Real,Petsc_Default_Real, &
         Petsc_Default_Integer,ierr)
#endif
    CALL KSPSetFromOptions(Krylov,ierr)
  END SUBROUTINE SetupKSPSolver

  ! Scatter U and get all local values
  SUBROUTINE GetVec_U
    IMPLICIT NONE
    IF (dyn) THEN
       CALL GetVec(Scatter_dyn,Vec_U_dyn,Seq_U_dyn,uu_dyn)
    ELSE
       CALL GetVec(Scatter,Vec_U,Seq_U,uu)
    END IF
  END SUBROUTINE GetVec_U

  SUBROUTINE GetVec_S
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_SS,Seq_SS,ss)
    CALL GetVec(Scatter_u,Vec_SH,Seq_SH,sh)
  END SUBROUTINE GetVec_S

  SUBROUTINE GetVec_dip_nrm
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_dip,Seq_dip,dip)
    CALL GetVec(Scatter_u,Vec_nrm,Seq_nrm,nrm)
  END SUBROUTINE GetVec_dip_nrm

  SUBROUTINE GetVec_f2s_seq
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_f2s,Seq_f2s,f2s)
  END SUBROUTINE GetVec_f2s_seq

  SUBROUTINE GetVec_fp
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_fp,Seq_fp,fp)
  END SUBROUTINE GetVec_fp

  SUBROUTINE GetVec_fl
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_fl,Seq_fl,fl)
  END SUBROUTINE GetVec_fl

  SUBROUTINE GetVec_flc
    IMPLICIT NONE
    CALL GetVec(Scatter_u,Vec_flc,Seq_flc,flc)
  END SUBROUTINE GetVec_flc

  SUBROUTINE GetVec_qu
    IMPLICIT NONE
    CALL GetVec(Scatter_q,Vec_qu,Seq_qu,qu)
  END SUBROUTINE GetVec_qu

  SUBROUTINE GetVec_ql
    IMPLICIT NONE
    CALL GetVec(Scatter_q,Vec_ql,Seq_ql,ql)
  END SUBROUTINE GetVec_ql

  SUBROUTINE GetVec(Scatter,Vec_X,Vec_Y,array)
    IMPLICIT NONE
    Vec :: Vec_X, Vec_Y
    VecScatter :: Scatter
    REAL(8) :: array(:)
    REAL(8),POINTER :: pntr(:)
    CALL VecScatterBegin(Scatter,Vec_X,Vec_Y,Insert_Values,Scatter_Forward,ierr)
    CALL VecScatterEnd(Scatter,Vec_X,Vec_Y,Insert_Values,Scatter_Forward,ierr)
    CALL VecGetArrayF90(Vec_Y,pntr,ierr)
    array=pntr
    CALL VecRestoreArrayF90(Vec_Y,pntr,ierr)
  END SUBROUTINE GetVec

END PROGRAM main
