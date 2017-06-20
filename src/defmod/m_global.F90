! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

MODULE global

#include <petscversion.h>

  USE local
  IMPLICIT NONE
#include "petscdef.h"
  ! Global variables
  INTEGER :: nnds,nels,nmts,nceqs,nfrcs,ntrcs,nabcs,p,frq,dsp,dsp_hyb,lm_str,  &
     bod_frc,steps,tstep,steps_dyn,tstep_dyn,nfnd,hyb,frq_dyn,nobs,n_log,      &
     n_log_dyn,n_log_wave,ih,igf,rsf,n_log_slip,nceqs_ncf,init,frq_wave,       &
     frq_slip,n_lmnd,lmnd0,nobs_loc,nfnd_loc,ngp,ngp_loc,DPFlag,nrxs
  REAL(8) :: alpha,beta,t,dt,t_dyn,dt_dyn,t_lim,val,wt,scale,rslip,t_sta,t_abs,&
     t_hyb,v_bg,vtol,trunc
  INTEGER,ALLOCATABLE :: nodes(:,:),bc(:,:),id(:),work(:),fnode(:),telsd(:,:), &
     worku(:),workl(:),node_pos(:),node_neg(:),node_all(:),slip(:),perm(:),    &
     onlst(:,:),frc(:),slip_sum(:),slip0(:),idgp(:,:),idgp_loc(:,:),           &
     gpnlst(:,:),nnd_fe2fd(:),rsf_sta(:),rxnode(:)
  REAL(8),ALLOCATABLE :: coords(:,:),mat(:,:),stress(:,:,:),vvec(:),cval(:,:), &
     fval(:,:),tval(:,:),vvec_all(:,:),vecf(:,:),fc(:),matf(:,:),st_init(:,:), &
     xfnd(:,:),ocoord(:,:),oshape(:,:),fcd(:),dc(:),rsfb0(:),rsfV0(:),         &
     rsfdtau0(:),rsfa(:),rsfb(:),rsfL(:),rsftheta(:),coh(:),dcoh(:),mu_hyb(:), &
     mu_cap(:),rsfv(:),ocoord_loc(:,:),xgp(:,:),gpshape(:,:),rx_press(:),      &
     row_replace(:,:),rxval(:,:)
  REAL(8),ALLOCATABLE,TARGET :: uu(:),tot_uu(:),uup(:),uu_dyn(:),tot_uu_dyn(:),&
     fl(:),ql(:),flc(:),fp(:),qu(:),ss(:),sh(:),f2s(:),dip(:),nrm(:),          &
     flt_slip(:),tot_flt_slip(:),qs_flt_slip(:)
  CHARACTER(12) :: stype
  CHARACTER(256) :: output_file
  LOGICAL :: poro,visco,fault,dyn,fail,dsp_dyn,crp,gf
  Vec :: Vec_F,Vec_U,Vec_Um,Vec_Up,Vec_lambda,Vec_I,Vec_lambda_tot,            &
     Vec_U_dyn,Vec_Um_dyn,Vec_U_dyn_tot,Vec_Up0,Vec_Up_dyn,Vec_I_dyn,Vec_fp,   &
     Vec_qu,Vec_Uu,Vec_Ul,Vec_fl,Vec_flc,Vec_ql,Vec_SS,Vec_SH,Vec_f2s,Vec_dip, &
     Vec_nrm,Vec_lmnd,Vec_lambda_sta,Vec_lambda_sta0
  Vec,POINTER :: Vec_W(:),Vec_Wlm(:)
  Mat :: Mat_K,Mat_M,Mat_Minv,Mat_Gt,Mat_G,Mat_GMinvGt,Mat_Kc,Mat_K_dyn,Mat_H, &
     Mat_Ht
  KSP :: Krylov
  PC :: PreCon
  ! Local element/side/node variables
  INTEGER :: el,side,node
  REAL(8) :: E,nu,dns,visc,expn,H,B,phi,Kf
  INTEGER,ALLOCATABLE :: indx(:),indxp(:),enodes(:),indx_dyn(:)
  REAL(8),ALLOCATABLE :: k(:,:),m(:,:),f(:),ecoords(:,:),kc(:,:),Hs(:),        &
     k_dyn(:,:),uu_obs(:,:),tot_uu_obs(:,:),uu_dyn_obs(:,:),                   &
     tot_uu_dyn_obs(:,:),flt_ss(:,:),flt_p(:),uu_fd(:,:)
  ! Variables for parallel code
  INTEGER :: nprcs,rank,ierr
  INTEGER,ALLOCATABLE :: epart(:),npart(:) ! Partitioning
  INTEGER,ALLOCATABLE :: nmap(:),emap(:),nl2g(:,:),indxmap(:,:),               &
     indxmap_u(:,:),FltMap(:,:),ol2g(:),gpl2g(:) ! L-G Mapping
  Vec :: Seq_U,Seq_U_dyn,Seq_fp,Seq_fl,Seq_flc,Seq_ql,Seq_qu,Seq_SS,Seq_SH,    &
     Seq_f2s,Seq_dip,Seq_nrm
  IS :: From,To,RI,From_u,To_u,RIu,From_p,To_p,RIl
  VecScatter :: Scatter,Scatter_dyn,Scatter_u,Scatter_q,Scatter_s2d
  REAL(8),POINTER :: pntr(:)

CONTAINS

  ! Form local [K] - Includes Dirichlet Pressure
  SUBROUTINE FormLocalK_DP(el,k,f,indx,strng)
    IMPLICIT NONE
    INTEGER :: el,indx(:)
    REAL(8) :: k(:,:),estress(nip,cdmn),f(:)
    CHARACTER(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    IF (visco) estress=stress(el,:,:)
    IF (dyn) THEN
       E=mat(id(el),5+4*p+init+1); nu=mat(id(el),5+4*p+init+2)
    ELSE
       E=mat(id(el),1); nu=mat(id(el),2)
    END IF
    visc=mat(id(el),3); expn=mat(id(el),4)
    IF ((.NOT. poro) .OR. strng=="Ke") THEN
       CALL FormElK(ecoords,estress,E,nu,visc,expn,dt,k,strng)
    ELSE
       H=mat(id(el),6)
       B=mat(id(el),7); phi=mat(id(el),8); Kf=mat(id(el),9)
       CALL FormElKp(ecoords,estress,E,nu,visc,expn,H,B,phi,Kf,1.0d0,scale,dt, &
          k,strng)
    END IF
    CALL AddWinklerFdn(el,k)
    IF (.NOT. dyn) CALL FixBCinLocalK_DP(el,k,f)
    IF (dyn) THEN
       CALL FormLocalIndx_dyn(enodes,indx)
    ELSE
       CALL FormLocalIndx(enodes,indx)
    END IF
  END SUBROUTINE FormLocalK_DP

  ! Form local [K]
  SUBROUTINE FormLocalK(el,k,indx,strng)
    IMPLICIT NONE
    INTEGER :: el,indx(:)
    REAL(8) :: k(:,:),estress(nip,cdmn)
    CHARACTER(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    IF (visco) estress=stress(el,:,:)
    IF (dyn) THEN
       E=mat(id(el),5+4*p+init+1); nu=mat(id(el),5+4*p+init+2)
    ELSE
       E=mat(id(el),1); nu=mat(id(el),2)
    END IF
    visc=mat(id(el),3); expn=mat(id(el),4)
    IF ((.NOT. poro) .OR. strng=="Ke") THEN
       CALL FormElK(ecoords,estress,E,nu,visc,expn,dt,k,strng)
    ELSE
       H=mat(id(el),6)
       B=mat(id(el),7); phi=mat(id(el),8); Kf=mat(id(el),9)
       CALL FormElKp(ecoords,estress,E,nu,visc,expn,H,B,phi,Kf,1.0d0,scale,dt, &
          k,strng)
    END IF
    CALL AddWinklerFdn(el,k)
    IF (.NOT. dyn) CALL FixBCinLocalK(el,k)
    IF (dyn) THEN
       CALL FormLocalIndx_dyn(enodes,indx)
    ELSE
       CALL FormLocalIndx(enodes,indx)
    END IF
  END SUBROUTINE FormLocalK

 
! Rescale local [Kv] for dt
  SUBROUTINE RscKv(el,k,indx,dtmp)
    IMPLICIT NONE
    INTEGER :: el,indx(:)
    REAL(8) :: k(:,:),estress(nip,cdmn),dtmp
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    estress=stress(el,:,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    H=mat(id(el),6)
    CALL RscElKv(ecoords,estress,E,nu,visc,expn,dt,k,dtmp)
    CALL FormLocalIndx(enodes,indx) 
  END SUBROUTINE RscKv

  ! Rescale local [Kp] for dt
  SUBROUTINE RscKp(el,k,indx,dtmp)
    IMPLICIT NONE
    INTEGER :: el,indx(:)
    REAL(8) :: k(:,:),estress(nip,cdmn),dtmp
    CHARACTER(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    IF (visco) estress=stress(el,:,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    H=mat(id(el),6)
    strng="Kp"
    IF (visco) strng="Kv"
    CALL RscElKp(ecoords,estress,E,nu,visc,expn,H,1.0d0,scale,dt,k,dtmp,strng)
    CALL FormLocalIndx(enodes,indx)
  END SUBROUTINE RscKp 

  ! Rescale [K] for new dt (poro and/or linear visco) 
  SUBROUTINE Rscdt(fdt)
    IMPLICIT NONE
#include "petsc.h" 
    INTEGER :: i,ndof
    REAL(8) :: fdt,dtmp
    dtmp=(fdt-f1)*dt
    ndof=eldof+eldofp
    DO i=1,nels
       IF (visco .AND. .NOT. poro) THEN
          CALL RscKv(i,k,indx,dtmp)
          !call RscRHSv(i,f,indx,dtmp)
       ELSE 
          CALL RscKp(i,k,indx,dtmp)
          !call RscRHSp(i,f,indx,dtmp)
       END IF
       indx=indxmap(indx,2)
       CALL MatSetValues(Mat_K,ndof,indx,ndof,indx,k,Add_Values,ierr)
       !call VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
    END DO
    CALL MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    CALL MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    !call VecAssemblyBegin(Vec_F,ierr)
    !call VecAssemblyEnd(Vec_F,ierr)
    dt=fdt*dt
  END SUBROUTINE Rscdt

  ! Form local [M]
  SUBROUTINE FormLocalM(el,m,indx)
    IMPLICIT NONE
    INTEGER :: el,indx(:)
    REAL(8) :: m(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    dns=mat(id(el),5)
    CALL FormElM(ecoords,dns,m)
    CALL FormElIndx(enodes,indx)
  END SUBROUTINE FormLocalM

  ! Account for constraint eqns
  SUBROUTINE ApplyConstraints
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: i,j,n,j1,j2,j3,j4
    OPEN(15,file="cnstrns.tmp",status='old')
    j4=0
    DO i=1,nceqs
       IF (poro .AND. MOD(i,dmn+1)>0) j4=j4+1
       READ(15,*)n
       DO j=1,n
          READ(15,*)vvec,node
          IF (stype/="explicit") THEN
             DO j1=1,dmn+p
                j2=(dmn+p)*node-(dmn+p)+j1-1; j3=(dmn+p)*nnds+i-1
                CALL MatSetValue(Mat_K,j2,j3,wt*vvec(j1),Add_Values,ierr)
                CALL MatSetValue(Mat_K,j3,j2,wt*vvec(j1),Add_Values,ierr)
             END DO
          END IF
          IF ((stype=="explicit" .AND. .NOT. gf).OR. (fault .AND.              &
             i<=nceqs_ncf)) THEN
             DO j1=1,dmn
                j2=dmn*node-dmn+j1-1; j3=i-1
                IF (poro .AND. MOD(i,dmn+1)>0) THEN
                   CALL MatSetValue(Mat_Gt,j2,j4-1,vvec(j1),Add_Values,ierr)
                ELSEIF (.NOT. poro) THEN
                   CALL MatSetValue(Mat_Gt,j2,j3,vvec(j1),Add_Values,ierr)
                END IF
             END DO
          END IF
       END DO
       IF (stype/="explicit") THEN ! Constraint block diagonals have to be
          j1=(dmn+p)*nnds+i-1      ! explicitly set to zero (PETSc req)
          CALL MatSetValue(Mat_K,j1,j1,f0,Add_Values,ierr)
       END IF
       READ(15,*)cval(i,:)
    END DO
    CLOSE(15)
  END SUBROUTINE ApplyConstraints

  ! Create full Mat_Gt for dynamic constraints
  SUBROUTINE GetMat_Gt
    IMPLICIT NONE
    INTEGER :: j,j1,j2,j3,j4,j5
#include "petsc.h"
    IF (rank==0) THEN
       DO j=1,nfnd
          DO j1=1,dmn
             j3=(j-1)*dmn+j1-1
             DO j2=1,2
                vvec=vvec_all(2*j3+j2,:)
                node=node_all(2*j3+j2)
                DO j4=1,dmn
                   j5=dmn*node-dmn+j4-1
                   CALL MatSetValue(Mat_Gt,j5,j3+nceqs_ncf*dmn/(dmn+p),vvec(j4)&
                      ,Add_Values,ierr)
                END DO
             END DO
          END DO
       END DO
    END IF 
  END SUBROUTINE GetMat_Gt

  ! Scatter LMs to solution space
  SUBROUTINE GetVec_flambda
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j3,j4,u1,u2,workpos(dmn),workneg(dmn),row
    REAL(8) :: lm(dmn),vecfl(dmn),lmq,q
    CALL VecGetOwnershipRange(Vec_Ul,u1,u2,ierr)
    DO j=1,nfnd
       lm=f0; vecfl=f0; lmq=f0; q=f0
       DO j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          IF (poro) THEN
             row=dmn*(j-1)+SUM(perm(1:j-1))+j1-1
          ELSE
             row=dmn*(j-1)+j1-1
          END IF
          IF (row>=u1 .AND. row<u2) THEN
             CALL VecGetValues(Vec_Ul,1,row,lm(j1),ierr)
          END IF
       END DO
       CALL MPI_Reduce(lm,vecfl,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
          ierr)
       IF (poro .AND. perm(j)==1) THEN
           row=row+1
           IF (row>=u1 .AND. row<u2) THEN
              CALL VecGetValues(Vec_Ul,1,row,lmq,ierr)
           END IF
           CALL MPI_Reduce(lmq,q,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,   &
              ierr)
       END IF
       IF (rank==nprcs-1) THEN 
          CALL VecSetValues(Vec_fl,dmn,workpos,wt*vecfl,Insert_Values,ierr)
          CALL VecSetValues(Vec_fl,dmn,workneg,-wt*vecfl,Insert_Values,ierr)
          IF (poro .AND. perm(j)==1) THEN
             j3=node_pos(j)-1
             j4=node_neg(j)-1
             CALL VecSetValue(Vec_ql,j3,wt*q,Insert_Values,ierr)
             CALL VecSetValue(Vec_ql,j4,-wt*q,Insert_Values,ierr)
          END IF
       END IF
    END DO
    CALL VecAssemblyBegin(Vec_fl,ierr)
    CALL VecAssemblyEnd(Vec_fl,ierr)
    IF (poro) THEN
       CALL VecAssemblyBegin(Vec_ql,ierr)
       CALL VecAssemblyEnd(Vec_ql,ierr)
    END IF
  END SUBROUTINE GetVec_flambda

  ! Extract the LM in fault's strike, dip and normal directions
  SUBROUTINE GetVec_fcoulomb
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,u1,u2,workpos(dmn),workneg(dmn),row
    REAL(8) :: lm(dmn),vecfl(dmn),mattmp(dmn,dmn),vectmp(dmn,1)
    CALL VecGetOwnershipRange(Vec_Ul,u1,u2,ierr)
    DO j=1,nfnd
       lm=f0; vecfl=f0
       DO j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          IF (poro) THEN
             row=dmn*(j-1)+SUM(perm(1:j-1))+j1-1
          ELSE
             row=dmn*(j-1)+j1-1
          END IF
          IF (row>=u1 .AND. row<u2) THEN
             CALL VecGetValues(Vec_Ul,1,row,lm(j1),ierr)
          END IF
       END DO
       CALL MPI_Reduce(lm,vecfl,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
          ierr)
       IF (rank==nprcs-1) THEN
          vectmp=RESHAPE(vecfl,(/dmn,1/))
          mattmp=TRANSPOSE(RESHAPE(vecf(j,:),(/dmn,dmn/)))
          vectmp=MATMUL(mattmp,vectmp)
          vecfl(:)=vectmp(:,1)
          CALL VecSetValues(Vec_flc,dmn,workpos,wt*vecfl,Insert_Values,ierr)
          CALL VecSetValues(Vec_flc,dmn,workneg,-wt*vecfl,Insert_Values,ierr)
       END IF
    END DO
    CALL VecAssemblyBegin(Vec_flc,ierr)
    CALL VecAssemblyEnd(Vec_flc,ierr)
  END SUBROUTINE GetVec_fcoulomb

  ! Pass pseudo velocity to dynamic model
  SUBROUTINE Rsfv2dyn
    IMPLICIT NONE
#include "petsc.h"
    INTEGER:: j,j1,j2,j3,rw_loc(dmn)
    REAL(8) :: vec(dmn),flt_qs(dmn)
    REAL(8),TARGET :: flt_ndf(n_lmnd*dmn)
    CALL VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_ndf=pntr
    CALL VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    DO j=1,nfnd_loc
       j1=FltMap(j,1); j3=FltMap(j,2)
       rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
       SELECT CASE(dmn)
       CASE(2)
          vec(1)=rsfv(j3)
          vec(2)=f0 ! Zero normal velocity
       CASE (3)
          flt_qs=flt_ndf(rw_loc)+st_init(j3,:)
          vec(1)=rsfv(j3)*flt_qs(1)/SQRT(flt_qs(1)**2+flt_qs(2)**2)
          vec(2)=rsfv(j3)*flt_qs(2)/SQRT(flt_qs(1)**2+flt_qs(2)**2)
          vec(3)=f0 ! Zero normal velocity
       END SELECT
       rw_loc=lmnd0*dmn+rw_loc-1
       CALL VecSetValues(Vec_I_dyn,dmn,rw_loc,vec*dt_dyn,Insert_Values,ierr)
    END DO
    CALL VecAssemblyBegin(Vec_I_dyn,ierr) 
    CALL VecAssemblyEnd(Vec_I_dyn,ierr)
  END SUBROUTINE Rsfv2Dyn

  ! Force to stress ratio and normal dip vector of the fault nodes
  SUBROUTINE GetVec_f2s
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,workpos(dmn),workneg(dmn)
    REAL(8) :: vecfl(dmn),vecss(dmn),vecsh(dmn),r(dmn),dip(dmn),nrm(dmn),      &
       matrot(dmn,dmn),matst(dmn,dmn),st(dmn,dmn),vec(dmn)
    CALL VecGetOwnershipRange(Vec_flc,j2,j3,ierr)
    DO j=1,nfnd
       matrot=RESHAPE(vecf(j,:),(/dmn,dmn/))
       vecfl=f0; vecss=f0; vecsh=f0
       DO j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          IF (workpos(j1)>=j2 .AND. workpos(j1)<j3) THEN
             CALL VecGetValues(Vec_flc,1,workpos(j1),vecfl(j1),ierr)
             CALL VecGetValues(Vec_SS,1,workpos(j1),vecss(j1),ierr)
             CALL VecGetValues(Vec_SH,1,workpos(j1),vecsh(j1),ierr)
          END IF
       END DO
       CALL MPI_Reduce(vecfl,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecfl=vec
       CALL MPI_Reduce(vecss,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecss=vec
       CALL MPI_Reduce(vecsh,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecsh=vec
       IF (rank==nprcs-1) THEN
          SELECT CASE(dmn)
          CASE(2)
             st(1,1)=vecss(1); st(2,2)=vecss(2)
             st(1,2)=vecsh(1); st(2,1)=vecsh(2)
             dip(:)=matrot(:,1)
             nrm(:)=matrot(:,2)
          CASE(3)
             st(1,1)=vecss(1); st(2,2)=vecss(2); st(3,3)=vecss(3)
             st(1,2)=vecsh(1); st(2,3)=vecsh(2); st(1,3)=vecsh(3)
             st(2,1)=vecsh(1); st(3,2)=vecsh(2); st(3,1)=vecsh(3)
             dip(:)=matrot(:,2)
             nrm(:)=matrot(:,3)
          END SELECT
          matst=MATMUL(MATMUL(TRANSPOSE(matrot),st),matrot)
          vecss(dmn)=matst(dmn,dmn)
          IF (ABS(vecfl(dmn))>f0) r=ABS(vecss(dmn)/vecfl(dmn))
          !if (rank==nprcs-1) then
          CALL VecSetValues(Vec_f2s,dmn,workpos,r,Insert_Values,ierr)
          CALL VecSetValues(Vec_dip,dmn,workpos,dip,Insert_Values,ierr)
          CALL VecSetValues(Vec_nrm,dmn,workpos,nrm,Insert_Values,ierr)
          !end if
       END IF
       matrot=RESHAPE(vecf(j,:),(/dmn,dmn/))
       vecfl=f0; vecss=f0; vecsh=f0
       DO j1=1,dmn
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          IF (workneg(j1)>=j2 .AND. workneg(j1)<j3) THEN
             CALL VecGetValues(Vec_flc,1,workneg(j1),vecfl(j1),ierr)
             CALL VecGetValues(Vec_SS,1,workneg(j1),vecss(j1),ierr)
             CALL VecGetValues(Vec_SH,1,workneg(j1),vecsh(j1),ierr)
          END IF
       END DO
       CALL MPI_Reduce(vecfl,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecfl=vec
       CALL MPI_Reduce(vecss,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecss=vec
       CALL MPI_Reduce(vecsh,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecsh=vec
       IF (rank==nprcs-1) THEN
          SELECT CASE(dmn)
          CASE(2)
             st(1,1)=vecss(1); st(2,2)=vecss(2)
             st(1,2)=vecsh(1); st(2,1)=vecsh(2)
             dip(:)=matrot(:,1)
             nrm(:)=matrot(:,2)
          CASE(3)
             st(1,1)=vecss(1); st(2,2)=vecss(2); st(3,3)=vecss(3)
             st(1,2)=vecsh(1); st(2,3)=vecsh(2); st(1,3)=vecsh(3)
             st(2,1)=vecsh(1); st(3,2)=vecsh(2); st(3,1)=vecsh(3)
             dip(:)=matrot(:,2)
             nrm(:)=matrot(:,3)
          END SELECT
          matst=MATMUL(MATMUL(TRANSPOSE(matrot),st),matrot)
          vecss(dmn)=matst(dmn,dmn)
          IF (ABS(vecfl(dmn))>f0) r=(r+ABS(vecss(dmn)/vecfl(dmn)))/f2
          ! Convert prestress to nodal force
          IF (r(dmn)>f0) THEN
             st_init(j,:)=st_init(j,:)/r
             coh(j)=coh(j)/r(1)
             IF (rsf==1) rsfdtau0(j)=rsfdtau0(j)/r(1)
          END IF
          CALL VecSetValues(Vec_f2s,dmn,workneg,-r,Insert_Values,ierr)
          CALL VecSetValues(Vec_dip,dmn,workneg,dip,Insert_Values,ierr)
          CALL VecSetValues(Vec_nrm,dmn,workneg,nrm,Insert_Values,ierr)
       END IF
    END DO
    CALL VecAssemblyBegin(Vec_f2s,ierr)
    CALL VecAssemblyEnd(Vec_f2s,ierr)
    CALL VecAssemblyBegin(Vec_dip,ierr)
    CALL VecAssemblyEnd(Vec_dip,ierr)
    CALL VecAssemblyBegin(Vec_nrm,ierr)
    CALL VecAssemblyEnd(Vec_nrm,ierr)
    CALL MPI_Bcast(st_init,nfnd*dmn,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
    CALL MPI_Bcast(coh,nfnd,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
    IF (rsf==1) CALL MPI_Bcast(rsfdtau0,nfnd,MPI_Real8,nprcs-1,                &
       MPI_Comm_World,ierr)
  END SUBROUTINE GetVec_f2s

  ! Fault dip and normal vectors
  SUBROUTINE GetVec_ft
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,workpos(dmn)
    REAL(8):: dip(dmn),nrm(dmn),matrot(dmn,dmn)
    DO j=1,nfnd
       DO j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
       END DO
       matrot=RESHAPE(vecf(j,:),(/dmn,dmn/))
       SELECT CASE(dmn)
       CASE(2)
          dip(:)=matrot(:,1)
          nrm(:)=matrot(:,2)
       CASE(3)
          dip(:)=matrot(:,2)
          nrm(:)=matrot(:,3)
       END SELECT
       IF (rank==nprcs-1) THEN
          CALL VecSetValues(Vec_dip,dmn,workpos,dip,Insert_Values,ierr)
          CALL VecSetValues(Vec_nrm,dmn,workpos,nrm,Insert_Values,ierr)
       END IF
    END DO
    CALL VecAssemblyBegin(Vec_dip,ierr)
    CALL VecAssemblyEnd(Vec_dip,ierr)
    CALL VecAssemblyBegin(Vec_nrm,ierr)
    CALL VecAssemblyEnd(Vec_nrm,ierr)
  END SUBROUTINE GetVec_ft

  ! RSF pseudo time update 
  SUBROUTINE RSF_QS_update(flt_ndf0,flt_ndf1,slip_loc,trunc)
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,j4,rw_loc(dmn),nqs,cut_glb,cut,slip_loc(nfnd),ntol,nslip
    REAL(8) :: theta,mu,a,b0,b,V0,L,dd,v_qs,dtpsd,flt_qs0(dmn),flt_qs1(dmn),   &
       flt_qs(dmn),trunc !,vmax
    REAL(8),TARGET :: flt_ndf0(n_lmnd*dmn),flt_ndf1(n_lmnd*dmn)
    REAL(8),ALLOCATABLE :: rsfstate(:,:,:)
    INTEGER,SAVE :: k=0,rsflog=0
    dtpsd=f1
    nqs=INT((dt+trunc)/dtpsd)
    ALLOCATE(rsfstate(nqs,nfnd_loc,3))
    rsfstate=f0; slip_loc=0
    cut=nqs
    vtol=1.D-3 ! Velocity threshold 1D-5
    !vmax=1.D-2 ! Maximum velocity 1D-4
    ntol=10 ! Slip node threshold 10 for 5 m spacing
    DO j4=1,nqs
       DO j=1,nfnd_loc
          j1=FltMap(j,1); j3=FltMap(j,2)
          ! RSF parameters
          a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
          IF (j4==1) THEN
             theta=rsftheta(j3)
          ELSE
             theta=rsfstate(j4-1,j,3)
          END IF
          rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
          flt_qs0=flt_ndf0(rw_loc)+st_init(j3,:)
          flt_qs1=flt_ndf1(rw_loc)+st_init(j3,:)
          ! Match shear with friction by updating v_qs
          flt_qs=flt_qs0+(flt_qs1-flt_qs0)*DBLE(j4)/DBLE(nqs)
          mu=SQRT(SUM(flt_qs(:dmn-1)*flt_qs(:dmn-1)))/ABS(flt_qs(dmn))
          !v_qs=sinh(mu/a)*V0*f2/exp((b0+b*log(V0*theta/L))/a)
          v_qs=MIN(vtol,SINH(mu/a)*V0*f2/EXP((b0+b*LOG(V0*theta/L))/a))
          dd=v_qs*dtpsd
          theta=dtpsd/(f1+dd/f2/L)+theta*(f1-dd/L/f2)/(f1+dd/L/f2)
          rsfstate(j4,j,:)=(/v_qs,mu,theta/)
       END DO
       CALL MPI_AllReduce(SIZE(PACK(rsfstate(j4,:,1),rsfstate(j4,:,1)>=vtol)),nslip,1,MPI_Integer,MPI_Sum,MPI_Comm_World,ierr)
       !if (maxval(rsfstate(j4,:,1))>vtol) then
       IF (nslip>=ntol) THEN ! At least ntol fault nodes nucleate
          cut=j4
          EXIT
       END IF
    END DO
    CALL MPI_AllReduce(cut,cut_glb,1,MPI_Integer,MPI_Min,MPI_Comm_World,ierr) 
    DO j=1,nfnd_loc
       j3=FltMap(j,2)
       rsftheta(j3)=rsfstate(cut_glb,j,3)
       mu_hyb(j3)=rsfstate(cut_glb,j,2)
       rsfv(j3)=rsfstate(cut_glb,j,1)
       IF (rsfv(j3)>=vtol .AND. nslip>=ntol) slip_loc(j3)=1 ! .and. nslip>ntol
    END DO 
    IF (cut_glb<nqs) THEN 
       CALL VecScale(Vec_lambda_sta,DBLE(cut_glb)/DBLE(nqs),ierr)
       CALL VecAXPY(Vec_lambda_sta,f1-DBLE(cut_glb)/DBLE(nqs),Vec_lambda_sta0, &
          ierr)
       trunc=dt*(f1-DBLE(cut_glb)/DBLE(nqs))
    ELSE
       trunc=f0
    END IF
    DO j4=1,cut_glb
       IF (MOD(k,200)==0) THEN
          CALL WriteOutput_rsf(rsfstate(j4,:,:2))
          rsflog=rsflog+1
       END IF
       k=k+1
    END DO
    IF (rank==0) CALL WriteOutput_log_rsf(rsflog,dtpsd*DBLE(200))
  END SUBROUTINE RSF_QS_update

  ! Write rate-state pseudo velocity and friction
  SUBROUTINE WriteOutput_rsf(state)
    IMPLICIT NONE
    INTEGER :: j,j3
    INTEGER,SAVE :: k=0
    REAL(8) :: state(nfnd_loc,2)
    CHARACTER(256) :: name,name0,name1
    CHARACTER(64) :: fmt
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    IF (nfnd_loc>0) THEN 
       WRITE(name,'(A,A,A,I0.6,A)')TRIM(name0),TRIM(name1),"_",rank,"_rsf.txt"
       fmt="(2(ES11.2E3,X))"
       IF (k==0) THEN
          OPEN(10,file=ADJUSTL(name),status='replace')
          WRITE(10,'(I0)')nfnd_loc
          DO j=1,nfnd_loc
             j3=FltMap(j,2)
             SELECT CASE(dmn)
                CASE(2); WRITE(10,'(2(F0.3,X),I0)')xfnd(j3,:),j3
                CASE(3); WRITE(10,'(3(F0.3,X),I0)')xfnd(j3,:),j3
             END SELECT
          END DO
       ELSE
          OPEN(10,file=ADJUSTL(name),status='old',position='append',action=    &
             'write')
       END IF
       DO j=1,nfnd_loc
          WRITE(10,fmt)state(j,:)
       END DO
       CLOSE(10); k=k+1
    END IF
  END SUBROUTINE WriteOutput_rsf

  SUBROUTINE WriteOutput_log_rsf(n,dt)
    IMPLICIT NONE
    INTEGER :: n
    INTEGER,SAVE :: k=0
    REAL(8) :: dt
    CHARACTER(256) :: name,name0,name1
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,A,A,I0.6,A)')TRIM(name0),TRIM(name1),"_rsf.log"
    IF (k==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       WRITE(10,'(F0.3)')dt
    ELSE 
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    WRITE(10,'(I0)')n
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_log_rsf

  ! Update slip from the static model (from Vec_lambda_sta)
  SUBROUTINE GetSlip_sta 
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,slip_loc(nfnd),rw_loc(dmn)
    REAL(8) :: mu,theta102,flt_qs(dmn),fsh,fnrm,rsftau,mattmp(dmn,dmn),d,fcoh, &
       a,b0,b,V0,L
    REAL(8),TARGET :: flt_ndf(n_lmnd*dmn),flt_ndf0(n_lmnd*dmn)
    INTEGER,SAVE :: k=0
    CALL VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_ndf=pntr
    CALL VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    IF (rsf==1 .AND. k>0) THEN
        CALL VecGetArrayF90(Vec_lambda_sta0,pntr,ierr)
        flt_ndf0=pntr
        CALL VecRestoreArrayF90(Vec_lambda_sta0,pntr,ierr)
        CALL RSF_QS_update(flt_ndf0,flt_ndf,slip_loc,trunc)
        go to 250
    END IF
    IF (rsf==1 .AND. k==0) rsfv=v_bg
    slip_loc=0
    rsftau=f0 
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       flt_qs=flt_ndf(rw_loc)+st_init(j3,:)
       IF (rsf==1) THEN  
          IF (k==0) THEN
             ! RSF parameters
             a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
             mu_cap(j3)=a*asinh(v_bg/V0/f2*EXP((b0+b*LOG(V0*rsftheta(j3)/L))/a))
             mu=MIN(mu_cap(j3),SQRT(SUM(flt_qs(:dmn-1)*flt_qs(:dmn-1)))        &
                /ABS(flt_qs(dmn)))
             mu_hyb(j3)=mu
             theta102=rsftheta(j3)
             rsftheta(j3)=L/V0*EXP((a*LOG(f2*SINH(mu/a))-b0-a*LOG(v_bg         &
                /V0))/b)
             ! Only for SCEC102
             mu=mu_cap(j3)
             rsftheta(j3)=theta102
             CALL GetExSt(xfnd(j3,:),t_abs,rsfdtau0(j3),rsftau)
          ELSE
             mu=mu_hyb(j3)
          END IF
       ELSE ! Slip weakening
          mu=fc(j3)
       END IF
       flt_qs(1)=flt_qs(1)+rsftau
       SELECT CASE(dmn)
       CASE(2)
          fsh=ABS(flt_qs(1))
          fnrm=flt_qs(2)
       CASE(3)
          fsh=SQRT(flt_qs(1)**2+flt_qs(2)**2)
          fnrm=flt_qs(3)
       END SELECT
       mattmp=TRANSPOSE(RESHAPE(vecf(j3,:),(/dmn,dmn/)))
       ! Cohesive stress if any
       IF (coh(j3)>f0) THEN
          d=SQRT(SUM(qs_flt_slip(rw_loc(:dmn-1))*qs_flt_slip(rw_loc(:dmn-1))))
          IF (d<dcoh(j3)) THEN 
             fcoh=coh(j3)*(f1-d/dcoh(j3))
          ELSE
             fcoh=f0
          END IF
       ELSE
           fcoh=f0
       END IF
       IF ((fnrm<f0 .AND. ABS((fsh-fcoh)/fnrm)>mu) .OR. fnrm>f0) THEN
          slip_loc(j3)=1
       ELSE
          slip_loc(j3)=0
       END IF
    END DO
250 CALL MPI_AllReduce(slip_loc,slip,nfnd,MPI_Integer,MPI_Sum,MPI_Comm_World,  &
       ierr)
    slip0=slip
    slip_sum=slip
    k=k+1
  END SUBROUTINE GetSlip_sta 

  ! Scatter from Vec_Ul to Vec_lambda_sta 
  SUBROUTINE LM_s2d
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,rw_dyn(dmn),rw_loc(dmn),rw_sta(dmn),                 &
       idxmp(nfnd_loc*dmn,2)
    INTEGER,SAVE :: k=0
    REAL(8) :: vec(dmn),mattmp(dmn,dmn),vectmp(dmn,1)
    REAL(8),TARGET :: flt_sta(n_lmnd*dmn)
    IF (k==0) THEN
       DO j1=1,nfnd_loc
          j=FltMap(j1,1); j3=FltMap(j1,2)
          rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
          rw_dyn=lmnd0*dmn+(/((j-1)*dmn+j2,j2=1,dmn)/)-1
          idxmp(rw_loc,1)=rw_dyn
          IF (poro) THEN
             rw_sta=(/((j3-1)*dmn+SUM(perm(1:j3-1))+j2-1,j2=1,dmn)/)
          ELSE
             rw_sta=(/((j3-1)*dmn+j2-1,j2=1,dmn)/)
          END IF
          idxmp(rw_loc,2)=rw_sta
       END DO
       CALL ISCreateGeneral(Petsc_Comm_World,nfnd_loc*dmn,idxmp(:,2),          &
          Petsc_Copy_Values,From,ierr)
       CALL ISCreateGeneral(Petsc_Comm_World,nfnd_loc*dmn,idxmp(:,1),          &
          Petsc_Copy_Values,To,ierr)
       CALL VecScatterCreate(Vec_Ul,From,Vec_lambda_sta,To,Scatter_s2d,ierr)
    END IF
    CALL VecScatterBegin(Scatter_s2d,Vec_Ul,Vec_lambda_sta,Insert_Values,      &
       Scatter_Forward,ierr)
    CALL VecScatterEnd(Scatter_s2d,Vec_Ul,Vec_lambda_sta,Insert_Values,        &
       Scatter_Forward,ierr)
    ! Rotate and scale
    CALL VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    CALL VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       rw_dyn=lmnd0*dmn+rw_loc-1
       vec=flt_sta(rw_loc)*wt
       vectmp=RESHAPE(vec,(/dmn,1/))
       mattmp=TRANSPOSE(RESHAPE(vecf(j3,:),(/dmn,dmn/)))
       vectmp=MATMUL(mattmp,vectmp)
       vec=vectmp(:,1)
       CALL VecSetValues(Vec_lambda_sta,dmn,rw_dyn,vec,Insert_Values,ierr)
    END DO
    CALL VecAssemblyBegin(Vec_lambda_sta,ierr)
    CALL VecAssemblyEnd(Vec_lambda_sta,ierr)
    k=k+1
  END SUBROUTINE LM_s2d 

  ! Cap dynamic LM by frictional laws 
  SUBROUTINE CapLM_dyn
    IMPLICIT NONE
#include "petsc.h"
    INTEGER ::j,j1,j2,j3,rw_loc(dmn)
    REAL(8) :: d,dd,fr,frs,frd,fsh,mu,fcoh,fnrm,rsftau,vec_init(dmn),vec(dmn), &
       lm_sta(dmn),lm_dyn(dmn),lm_dyn0(dmn),mattmp(dmn,dmn)
    REAL(8),TARGET :: flt_sta(n_lmnd*dmn),flt_dyn(n_lmnd*dmn),                 &
       flt_dyn0(n_lmnd*dmn)
    CALL VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    CALL VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    CALL VecGetArrayF90(Vec_lambda,pntr,ierr)
    flt_dyn=pntr
    CALL VecRestoreArrayF90(Vec_lambda,pntr,ierr)
    CALL VecGetArrayF90(Vec_lambda_tot,pntr,ierr)
    flt_dyn0=pntr
    CALL VecRestoreArrayF90(Vec_lambda_tot,pntr,ierr)
    rsftau=f0
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       vec_init=st_init(j3,:)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       lm_sta=flt_sta(rw_loc)
       lm_dyn0=flt_dyn0(rw_loc)
       lm_dyn=flt_dyn(rw_loc)
       d=SQRT(SUM(tot_flt_slip(rw_loc(:dmn-1))*tot_flt_slip(rw_loc(:dmn-1))))
       vec=vec_init+lm_sta+lm_dyn0+lm_dyn 
       IF (rsf==1) THEN ! Rate state friction
          ! RSF options
          !dd=max(v_bg*dt_dyn,sqrt(sum(flt_slip(rw_loc(:dmn-1))                 &
          !   *flt_slip(rw_loc(:dmn-1)))))
          !dd=v_bg*dt_dyn+sqrt(sum(flt_slip(rw_loc(:dmn-1))                     &
          !   *flt_slip(rw_loc(:dmn-1))))
          !dd=max(rsfv(j3)*dt_dyn,sqrt(sum(flt_slip(rw_loc(:dmn-1))             &
          !   *flt_slip(rw_loc(:dmn-1)))))
          !dd=rsfv(j3)*dt_dyn+sqrt(sum(flt_slip(rw_loc(:dmn-1))                 &
          !   *flt_slip(rw_loc(:dmn-1))))
          dd=MAX(v_bg*dt_dyn,flt_slip(rw_loc(1))) ! Only for SCEC102

          ! Stabilize RSF velocity
          !if (tstep>1) then
          !   if (t_sta>f0 .and. dd<vtol*dt_dyn) rsf_sta(j1)=1 
          !   if (rsf_sta(j1)==1) dd=min(dd,vtol*dt_dyn)
          !end if

          rsftheta(j3)=MAX(dt_dyn,dt_dyn/(f1+dd/f2/rsfL(j3))                   &
             +rsftheta(j3)*(f1-dd/rsfL(j3)/f2)/(f1+dd/rsfL(j3)/f2))
          mu=rsfa(j3)*asinh(dd/dt_dyn/rsfV0(j3)/f2*EXP((rsfb0(j3)              &
             +rsfb(j3)*LOG(rsfV0(j3)*rsftheta(j3)/rsfL(j3)))/rsfa(j3)))
          ! Only for SCEC102
          CALL GetExSt(xfnd(j3,:),t_hyb,rsfdtau0(j3),rsftau)
       ELSE ! Slip weakening
          IF (d<dc(j3)) THEN
             mu=(f1-d/dc(j3))*(fc(j3)-fcd(j3))+fcd(j3)
          ELSE
             mu=fcd(j3)
          END IF
       END IF
       mu_hyb(j3)=mu
       vec(1)=vec(1)+rsftau
       SELECT CASE(dmn)
       CASE(2)
          fsh=ABS(vec(1))
          fnrm=vec(2)
       CASE(3)
          fsh=SQRT(vec(1)**2+vec(2)**2)
          fnrm=vec(3)
       END SELECT
       mattmp=TRANSPOSE(RESHAPE(vecf(j3,:),(/dmn,dmn/)))
       ! Slip weakening cohesion
       IF (coh(j3)>f0) THEN
          ! Cumulative slip 
          d=d+SQRT(SUM(qs_flt_slip(rw_loc(:dmn-1))*qs_flt_slip(rw_loc(:dmn-1))))
          IF (d<dcoh(j3)) THEN
             fcoh=coh(j3)*(f1-d/dcoh(j3))
          ELSE
             fcoh=f0
          END IF
       ELSE
          fcoh=f0
       END IF
       ! If LM exceeds maximum friction (mu*fn)
       IF (fnrm<f0 .AND. ABS(fsh)>mu*ABS(fnrm)+fcoh) THEN
          fr=mu*ABS(fnrm)+fcoh 
          frs=fr*ABS(vec(1)/fsh)
          frs=SIGN(frs,vec(1))
          vec(1)=frs
          IF (dmn==3) THEN
             frd=fr*ABS(vec(2)/fsh)
             frd=SIGN(frd,vec(2))
             vec(2)=frd
          END IF
       ! Zero friction when fault faces detach under cohesion
       ELSEIF (fnrm>=f0) THEN 
          vec=vec/SQRT(SUM(vec*vec))*fcoh
       END IF
       ! Subtract the static LM 
       vec(1)=vec(1)-rsftau
       lm_dyn=vec-vec_init-lm_sta-lm_dyn0
       rw_loc=lmnd0*dmn+rw_loc-1
       ! From the new dynamic LM
       CALL VecSetValues(Vec_lambda,dmn,rw_loc,lm_dyn,Insert_Values,ierr)
    END DO
    CALL VecAssemblyBegin(Vec_lambda,ierr)
    CALL VecAssemblyEnd(Vec_lambda,ierr)
  END SUBROUTINE CapLM_dyn

  ! Record latest moment fault stress from the hybrid run
  SUBROUTINE GetVec_lambda_hyb(trunc)
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: i,j,j1,j2,j3,rw_loc(dmn),nqs
    REAL(8) :: vec(dmn),lm_sta(dmn),lm_dyn(dmn),trunc,dtpsd,a,b0,b,V0,L,mu,dd
    REAL(8),TARGET :: flt_sta(n_lmnd*dmn),flt_dyn(n_lmnd*dmn)
    CALL VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    CALL VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    CALL VecGetArrayF90(Vec_lambda_tot,pntr,ierr)
    flt_dyn=pntr
    CALL VecRestoreArrayF90(Vec_lambda_tot,pntr,ierr)
    dtpsd=f1
    nqs=INT((dt-trunc)/dtpsd)
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       lm_sta=flt_sta(rw_loc)
       lm_dyn=flt_dyn(rw_loc)
       IF (rsf==1) THEN ! Pseudo time healing for RSF
          ! RSF parameters
          a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
          DO i=1,nqs 
             vec=lm_sta+lm_dyn*DBLE(i)/DBLE(nqs)
             mu=SQRT(SUM(vec(:dmn-1)*vec(:dmn-1)))/ABS(vec(dmn))
             rsfv(j3)=MIN(vtol,SINH(mu/a)*V0*f2/EXP((b0+b*LOG(V0*rsftheta(j3)/L))/a))
             dd=dtpsd*rsfv(j3)
             rsftheta(j3)=dtpsd/(f1+dd/f2/L)+rsftheta(j3)*(f1-dd/L/f2)/(f1+dd/L/f2)
          END DO
          !print*, mu,rsfv(j3),rsftheta(j3)
       END IF
       vec=lm_sta+lm_dyn
       rw_loc=lmnd0*dmn+rw_loc-1
       CALL VecSetValues(Vec_lambda_sta0,dmn,rw_loc,vec,Insert_Values,ierr)
    END DO
    CALL VecAssemblyBegin(Vec_lambda_sta0,ierr)
    CALL VecAssemblyEnd(Vec_lambda_sta0,ierr)
  END SUBROUTINE GetVec_lambda_hyb

  ! Time dependent initial stress (SCEC102)
  SUBROUTINE GetExSt(x,t,dst,st)
    IMPLICIT NONE
    REAL(8) :: t,x(dmn),rx2,dst,rsfF,rsfG,R2,rsfT,st
    rsfF=f0;rsfG=f1;R2=9.0;rsfT=f1;st=0
    IF (dmn==3) THEN
       rx2=x(2)**2+(x(3)+7.5)**2
    ELSE
       rx2=(x(2)+7.5)**2
    END IF
    IF (rx2<R2) rsfF=EXP(rx2/(rx2-R2))
    IF (t>f0 .AND. t<rsfT) rsfG=EXP((t-rsfT)**2/t/(t-2*rsfT))
    st=dst*rsfF*rsfG
  END SUBROUTINE GetExSt

  ! Determine if the dynamic slip is stabilized
  SUBROUTINE GetSlip_dyn
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,i,nc,nr,rw_loc(dmn),slip_loc(nfnd),slip_sum_loc(nfnd)
    REAL(8) :: d0,d
    slip_loc=slip
    slip_sum_loc=slip_sum
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       d=SQRT(SUM(flt_slip(rw_loc(:dmn-1))*flt_slip(rw_loc(:dmn-1))))
       d0=SQRT(SUM(tot_flt_slip(rw_loc(:dmn-1))*tot_flt_slip(rw_loc(:dmn-1))))
       ! Stabilization tolerance 
       IF (rsf==1) THEN ! Rate state friction
          IF (d0>5e-2*rsfL(j3)) THEN
             IF ((d/d0)<1E-4) THEN
                slip_loc(j3)=0
             ELSE
                slip_loc(j3)=1
                slip_sum_loc(j3)=1
             END IF
          END IF
          IF (d/dt_dyn>=vtol) THEN
            slip_loc(j3)=1
            slip_sum_loc(j3)=1
          END IF
       ELSE ! Slip weakening
          IF (d0>5E-2*dc(j3)) THEN 
             IF ((d/d0)<1E-4) THEN
                slip_loc(j3)=0
             ELSE
                slip_loc(j3)=1
                slip_sum_loc(j3)=1
             END IF
          END IF
       END IF
    END DO
    ! Zero off-rank slip
    IF (nfnd_loc>0) THEN
       j1=FltMap(1,2)-1; j2=FltMap(nfnd_loc,2)+1
       IF (j1>0) THEN
          slip_loc(:j1)=0
          slip_sum_loc(:j1)=0
       END IF
       IF (j2<=nfnd) THEN
          slip_loc(j2:)=0
          slip_sum_loc(j2:)=0
       END IF
    ELSE
       slip_loc=0
       slip_sum_loc=0
    END IF
    CALL MPI_AllReduce(slip_loc,slip,nfnd,MPI_Integer,MPI_Sum,                 &
       MPI_Comm_World,ierr)
    CALL MPI_AllReduce(slip_sum_loc,slip_sum,nfnd,MPI_Integer,MPI_Sum,         &
       MPI_Comm_World,ierr)
    ! Identify aseismic slip, nc=10 for SCEC10/14 (slow weakening)
    nc=10; nr=15
    IF (ih>nc+rsf*nr .AND. SUM((slip0-slip_sum)*(slip0-slip_sum))==0) THEN 
       slip=0
       crp=.TRUE.
    ELSEIF (ih>nc+rsf*nr) THEN
       i=0
       DO j=1,nfnd
          IF (slip0(j)==0 .AND. slip(j)==1) i=i+1
       END DO
       IF (i==0) slip=0 ! No slip except nucleation patch
    ELSEIF (ih<=nc+rsf*nr .AND. SUM(slip)==0) THEN 
       crp=.TRUE. ! Zero slip within time nc+rsf*nr 
    END IF
  END SUBROUTINE GetSlip_dyn 

  ! Get fault QS state 
  SUBROUTINE GetVec_flt_qs
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: l1,l2,j,j1,r1,r2,row_l,row_f2s,row_p 
    REAL(8) :: lm(dmn),vec(dmn),rvec(dmn),rf2s(dmn),mattmp(dmn,dmn),           &
       vectmp(dmn,1),pval,fltp
    CALL VecGetOwnershipRange(Vec_Um,l1,l2,ierr)
    CALL VecGetOwnershipRange(Vec_f2s,r1,r2,ierr)
    DO j=1,nfnd
       lm=f0; vec=f0; rf2s=0; rvec=f0; pval=f0; fltp=f0
       DO j1=1,dmn 
          IF (poro) THEN
             row_l=(dmn+1)*nnds+nceqs_ncf+(j-1)*dmn+SUM(perm(1:j-1))+j1-1
          ELSE
             row_l=dmn*nnds+nceqs_ncf+(j-1)*dmn+j1-1
          END IF
          row_f2s=dmn*node_pos(j)-dmn+j1-1
          IF (row_l>=l1 .AND. row_l<l2) THEN
             CALL VecGetValues(Vec_Um,1,row_l,lm(j1),ierr)
          END IF
          IF (row_f2s>=r1 .AND. row_f2s<r2) THEN
             CALL VecGetValues(Vec_f2s,1,row_f2s,rf2s(j1),ierr)
          END IF
       END DO
       IF (poro) THEN
          row_p=(dmn+1)*node_pos(j)-1
          IF (row_p>=l1 .AND. row_p<l2) THEN
             CALL VecGetValues(Vec_Um,1,row_p,pval,ierr)
          END IF
       END IF
       CALL MPI_Reduce(lm,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,ierr)
       CALL MPI_Reduce(rf2s,rvec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       CALL MPI_Reduce(pval,fltp,1,MPI_Real8,MPI_Sum,nprcs-1,                  &
          MPI_Comm_World,ierr)
       ! Rotate vec to fault coordinate
       IF (rank==nprcs-1) THEN
          vectmp=RESHAPE(vec,(/dmn,1/))
          mattmp=TRANSPOSE(RESHAPE(vecf(j,:),(/dmn,dmn/)))
          vectmp=MATMUL(mattmp,vectmp)
          vec(:)=vectmp(:,1)
          ! Nodal force to stress
          flt_ss(j,:)=vec*wt*rvec(dmn)
          IF (poro) flt_p(j)=fltp*scale 
       END IF
    END DO
    CALL MPI_Bcast(flt_ss,nfnd*dmn,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
    IF (poro) CALL MPI_Bcast(flt_p,nfnd,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
  END SUBROUTINE GetVec_flt_qs

  ! Add fault slip from dynamic model to static model as constraint functions
  SUBROUTINE FaultSlip
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: j,j1,j2,j3,rw_loc(dmn),rw_sta(dmn)
    REAL(8) :: vec(dmn),vectmp(dmn,1),mattmp(dmn,dmn)
    DO j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2) 
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       IF (poro) THEN
          rw_sta=(dmn+1)*nnds+nceqs_ncf+(/((j3-1)*dmn+SUM(perm(1:j3-1))        &
             +j2-1,j2=1,dmn)/)
       ELSE
          rw_sta=dmn*nnds+nceqs_ncf+(/((j3-1)*dmn+j2-1,j2=1,dmn)/)
       END IF
       vec=tot_flt_slip(rw_loc)
       ! Zero non-slip components, and rotate to Cartesian 
       vec(dmn)=f0 
       vectmp=RESHAPE(vec,(/dmn,1/))
       mattmp=RESHAPE(vecf(j3,:),(/dmn,dmn/))
       vectmp=MATMUL(mattmp,vectmp)
       vec=vectmp(:,1)
       CALL VecSetValues(Vec_F,dmn,rw_sta,vec*wt,Add_Values,ierr)
     END DO 
  END SUBROUTINE FaultSlip

  ! Extract solution at observation locations
  SUBROUTINE GetVec_obs
    IMPLICIT NONE
    INTEGER :: ob,i,j,ind(npel)
    INTEGER,ALLOCATABLE :: row(:)
    REAL(8) :: vecshp(npel,1)
    REAL(8),ALLOCATABLE :: vectmp(:,:),mattmp(:,:)
    IF (dyn .AND. .NOT. gf) THEN
       ALLOCATE(row(dmn*npel),vectmp(dmn,1),mattmp(dmn,npel))
    ELSE
       ALLOCATE(row((dmn+p)*npel),vectmp(dmn+p,1),mattmp(dmn+p,npel))
    END IF
    DO ob=1,nobs_loc
       ind=onlst(ob,:)
       IF (dyn .AND. .NOT. gf) THEN
          DO i=1,npel
             row((/((i-1)*dmn+j,j=1,dmn)/))=(/((ind(i)-1)*dmn+j,j=1,dmn)/)
          END DO
          mattmp=RESHAPE(uu_dyn(row),(/dmn,npel/))
          vecshp=RESHAPE(oshape(ob,:),(/npel,1/))
          vectmp=MATMUL(mattmp,vecshp)
          uu_dyn_obs(ob,:)=vectmp(:,1)
       ELSE
          DO i=1,npel
             row((/((i-1)*(dmn+p)+j,j=1,dmn+p)/))=(/((ind(i)-1)*(dmn+p)+j,j=1, &
               dmn+p)/)
          END DO
          mattmp=RESHAPE(uu(row),(/dmn+p,npel/))
          vecshp=RESHAPE(oshape(ob,:),(/npel,1/))
          vectmp=MATMUL(mattmp,vecshp) 
          uu_obs(ob,:)=vectmp(:,1)
          IF (poro) uu_obs(ob,dmn+1)=vectmp(dmn+1,1)*scale
       END IF
    END DO
  END SUBROUTINE GetVec_obs

  ! Apply nodal force
  SUBROUTINE ApplyNodalForce(node,vvec)
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: node,i,j
    REAL(8) :: vvec(:)
    DO i=1,dmn+p
       j=(dmn+p)*node-(dmn+p)+i-1
       val=vvec(i); IF (i==dmn+1) val=scale*val
       CALL VecSetValue(Vec_F,j,val,Add_Values,ierr)
    END DO
  END SUBROUTINE ApplyNodalForce

  ! Apply traction (EbEAve)
  SUBROUTINE ApplyTraction(el,side,vvec)
    IMPLICIT NONE
    INTEGER :: el,side,i,snodes(nps)
    REAL(8) :: vvec(:),area
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    CALL EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    vvec=vvec*area/DBLE(nps)
    snodes=nl2g(snodes,2)
    DO i=1,nps
       CALL ApplyNodalForce(snodes(i),vvec)
    END DO
  END SUBROUTINE ApplyTraction

  ! Form RHS
  SUBROUTINE FormRHS
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: i,j
    REAL(8) :: t1,t2
    ! Account for LM constraints
    DO i=1,nceqs
       t1=cval(i,2)/dt; t2=cval(i,3)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          j=(dmn+p)*nnds+i-1
          IF (stype/="explicit" .AND. rank==0) THEN
             val=wt*cval(i,1)
             CALL VecSetValue(Vec_F,j,val,Add_Values,ierr)
          END IF
       END IF
    END DO
    ! Account for nodal force / fluid source values
    DO i=1,nfrcs
       t1=fval(i,dmn+p+1)/dt; t2=fval(i,dmn+p+2)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          node=fnode(i); vvec=fval(i,1:dmn+p)
          IF (rank==0) CALL ApplyNodalForce(node,vvec)
       END IF
    END DO
    ! Account for traction BCs
    DO i=1,ntrcs
       t1=tval(i,dmn+p+1)/dt; t2=tval(i,dmn+p+2)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          el=telsd(i,1); side=telsd(i,2); vvec=tval(i,1:dmn+p)
          IF (el/=0) CALL ApplyTraction(el,side,vvec)
       END IF
    END DO
    
  END SUBROUTINE FormRHS

    ! Form RHS
  SUBROUTINE FormRHS_Debug
    IMPLICIT NONE
#include "petsc.h"
    INTEGER :: i,j,rowid,n
    REAL(8) :: t1,t2
    ! Account for LM constraints
    DO i=1,nceqs
       t1=cval(i,2)/dt; t2=cval(i,3)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          j=(dmn+p)*nnds+i-1
          IF (stype/="explicit" .AND. rank==0) THEN
             val=wt*cval(i,1)
             CALL VecSetValue(Vec_F,j,val,Add_Values,ierr)
          END IF
       END IF
    END DO
    ! Account for nodal force / fluid source values
    DO i=1,nfrcs
       t1=fval(i,dmn+p+1)/dt; t2=fval(i,dmn+p+2)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          node=fnode(i); vvec=fval(i,1:dmn+p)
          IF (rank==0) CALL ApplyNodalForce(node,vvec)
       END IF
    END DO
    ! Account for traction BCs
    DO i=1,ntrcs
       t1=tval(i,dmn+p+1)/dt; t2=tval(i,dmn+p+2)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          el=telsd(i,1); side=telsd(i,2); vvec=tval(i,1:dmn+p)
          IF (el/=0) CALL ApplyTraction(el,side,vvec)
       END IF
    END DO
    ! Account for prescribed/dirichlet values
    ! In practice only to be used for pressure
    ALLOCATE(row_replace(1,n))
    DO i=1,nrxs
       t1=rxval(i,dmn+p+1)/dt; t2=rxval(i,dmn+p+2)/dt
       IF (tstep>=NINT(t1) .AND. tstep<=NINT(t2)) THEN
          node=rxnode(i); vvec=rxval(i,1:dmn+p)
          IF (rank==0) THEN
             row_replace=f0
             rowid=(dmn+p)*node-(dmn+p)+(dmn+p)-1
             CALL MatZeroRows(Mat_K,1,rowid-1,0.0,ierr)
             CALL MatSetValue(Mat_K,rowid-1,rowid-1,f1,INSERT_VALUE,ierr)
             CALL ApplyNodalForce(node,vvec)
          END IF
       END IF
    END DO

  END SUBROUTINE FormRHS_Debug

  ! Form local damping matrix for elements with viscous dampers
  SUBROUTINE FormLocalAbsC(el,side,dir,m,indx)
    IMPLICIT NONE
    INTEGER :: el,side,dir,indx(:)
    REAL(8) :: m(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    dns=mat(id(el),5)
    CALL FormElAbsC(enodes,ecoords,side,dir,E,nu,dns,m)
    CALL FormElIndx(enodes,indx)
  END SUBROUTINE FormLocalAbsC

  ! Form local damping matrix for elements with viscous dampers
  SUBROUTINE FormLocalAbsC1(el,side,m,indx)
    IMPLICIT NONE
    INTEGER :: el,side,indx(:)
    REAL(8) :: m(:,:),matabs(dmn,dmn),vec1(dmn),vec2(dmn),vec3(dmn)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    dns=mat(id(el),5)
    SELECT CASE (eltype) 
    CASE("tri")
       SELECT CASE(side) 
       CASE(1); vec1=ecoords(2,:)-ecoords(1,:)
       CASE(2); vec1=ecoords(3,:)-ecoords(2,:)
       CASE(3); vec1=ecoords(1,:)-ecoords(3,:)
       END SELECT
       vec1=vec1/SQRT(SUM(vec1*vec1))
       vec2(1)=vec1(2); vec2(2)=-vec1(1)
       matabs(:,1)=vec1; matabs(:,2)=vec2
    CASE("qua")
       SELECT CASE(side) 
       CASE(1); vec1=ecoords(2,:)-ecoords(1,:)
       CASE(2); vec1=ecoords(3,:)-ecoords(2,:)
       CASE(3); vec1=ecoords(4,:)-ecoords(3,:)
       CASE(4); vec1=ecoords(1,:)-ecoords(4,:)
       END SELECT
       vec1=vec1/SQRT(SUM(vec1*vec1))
       vec2(1)=vec1(2); vec2(2)=-vec1(1)
       matabs(:,1)=vec1; matabs(:,2)=vec2
    CASE("tet")
       SELECT CASE(side) 
       CASE(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       CASE(2)
          vec1=ecoords(4,:)-ecoords(3,:)
          CALL Cross(vec1,ecoords(2,:)-ecoords(3,:),vec3)
       CASE(3)
          vec1=ecoords(4,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(3,:)-ecoords(1,:),vec3)
       CASE(4)
          vec1=ecoords(3,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       END SELECT
       vec1=vec1/SQRT(SUM(vec1*vec1))
       vec3=vec3/SQRT(SUM(vec3*vec3))
       CALL Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    CASE("hex")
       SELECT CASE(side) 
       CASE(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(5,:)-ecoords(1,:),vec3)
       CASE(2)
          vec1=ecoords(3,:)-ecoords(2,:)
          CALL Cross(vec1,ecoords(6,:)-ecoords(2,:),vec3)
       CASE(3)
          vec1=ecoords(4,:)-ecoords(3,:)
          CALL Cross(vec1,ecoords(7,:)-ecoords(3,:),vec3)
       CASE(4)
          vec1=ecoords(5,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       CASE(5)
          vec1=ecoords(4,:)-ecoords(1,:)
          CALL Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       CASE(6)
          vec1=ecoords(6,:)-ecoords(5,:)
          CALL Cross(vec1,ecoords(8,:)-ecoords(5,:),vec3)
       END SELECT
       vec1=vec1/SQRT(SUM(vec1*vec1))
       vec3=vec3/SQRT(SUM(vec3*vec3))
       CALL Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    END SELECT
    CALL FormElAbsC1(enodes,ecoords,side,matabs,E,nu,dns,m)
    CALL FormElIndx(enodes,indx)
  END SUBROUTINE FormLocalAbsC1

  SUBROUTINE Cross(a,b,r)
    IMPLICIT NONE
    REAL(8) :: a(3),b(3),r(3)
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
  END SUBROUTINE Cross

  ! Form local index
  SUBROUTINE FormLocalIndx(enodes,indx)
    IMPLICIT NONE
    INTEGER :: enodes(:),indx(:)
    IF (.NOT. poro) THEN
       CALL FormElIndx(enodes,indx)
    ELSE
       CALL FormElIndxp(enodes,indx)
    END IF
  END SUBROUTINE FormLocalIndx

  SUBROUTINE FormLocalIndx_dyn(enodes,indx)
    IMPLICIT NONE
    INTEGER :: enodes(:),indx(:)
    CALL FormElIndx(enodes,indx)
  END SUBROUTINE FormLocalIndx_dyn

  ! Recover stress
  SUBROUTINE RecoverStress(el,stress)
    IMPLICIT NONE
    INTEGER :: el,indxl(eldof)
    REAL(8) :: stress(:,:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    CALL FormLocalIndx(enodes,indx)
    indxl=indx(1:eldof)
    CALL CalcElStress(ecoords,uu(indxl),E,nu,stress(el,:,:))
  END SUBROUTINE RecoverStress

  ! Scatter stress to vertices (normal: SS, shear: SH)
  SUBROUTINE GetVec_Stress
    IMPLICIT NONE
#include "petsc.h"
    INTEGER:: i,j,indx(eldof),row(dmn)
    REAL(8):: sigma(cdmn)
    DO i=1,nels
       enodes=nodes(i,:)
       CALL FormElIndx(enodes,indx)
       indx=indxmap_u(indx,2)
       sigma=SUM(stress(i,:,:),dim=1)
       DO j=1,nip
          SELECT CASE(dmn)
          CASE(2)
             row=(/indx((j-1)*dmn+1),indx((j-1)*dmn+2)/)
             sigma(:)=stress(i,j,:)
             CALL VecSetValues(Vec_SS,dmn,row,sigma(1:dmn),Insert_Values,ierr)
             CALL VecSetValues(Vec_SH,dmn,row,(/sigma(dmn+1:cdmn),             &
                sigma(dmn+1:cdmn)/),Insert_Values,ierr)
          CASE(3)
             row=(/indx((j-1)*dmn+1),indx((j-1)*dmn+2),indx((j-1)*dmn+3)/)
             sigma(:)=stress(i,j,:)
             CALL VecSetValues(Vec_SS,dmn,row,sigma(1:dmn),Insert_Values,ierr)
             CALL VecSetValues(Vec_SH,dmn,row,sigma(dmn+1:cdmn),Insert_Values, &
                ierr)
          END SELECT
       END DO
    END DO
    CALL VecAssemblyBegin(Vec_SS,ierr)
    CALL VecAssemblyEnd(Vec_SS,ierr)
    CALL VecAssemblyBegin(Vec_SH,ierr)
    CALL VecAssemblyEnd(Vec_SH,ierr)
  END SUBROUTINE GetVec_Stress

  ! Form local Hs
  SUBROUTINE FormLocalHs(el,Hs,indxp)
    IMPLICIT NONE
    REAL(8) :: Hs(:)
    INTEGER :: el,j,indxp(:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    CALL FormElIndxp(enodes,indx)
    indxp=indx(eldof+1:)
    E=mat(id(el),1); nu=mat(id(el),2)
    CALL FormElHs(ecoords,uu(indxp),E,nu,scale,Hs)
    DO j=1,npel
       IF (bc(nodes(el,j),dmn+1)==0) Hs(j)=f0
    END DO
  END SUBROUTINE FormLocalHs

  ! Signed distance point to plane/line
  SUBROUTINE Mix3D(a,b,c,m)
    IMPLICIT NONE
    REAL(8) :: a(3),b(3),c(3),r(3),m 
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
    m=(r(1)*c(1)+r(2)*c(2)+r(3)*c(3))/(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
  END SUBROUTINE Mix3D 
  SUBROUTINE Mix2D(a,b,m)
    IMPLICIT NONE
    REAL(8) :: a(2),b(2),a3(3),b3(3),c3(3),m 
    a3=(/f0,f0,f1/); b3=(/a(:),f0/); c3=(/b(:),f0/)
    CALL Mix3D(a3,b3,c3,m) 
  END SUBROUTINE Mix2D

  ! Observation/FD nodal base 
  SUBROUTINE GetObsNd(strng)
    IMPLICIT NONE
    CHARACTER(2) :: strng
    INTEGER :: typ,nloop,ob,el
    INTEGER,ALLOCATABLE :: nd_full(:,:),pick(:)
    REAL(8) :: xmin,xmax,ymin,ymax,zmin,zmax,xmind,xmaxd,ymind,ymaxd,zmind,    &
       zmaxd,dd,du,dl,dr,df,db,d,eta,nu,psi,xob(dmn),N(npel),c,vec12(dmn),     &
       vec13(dmn),vec14(dmn),vec23(dmn),vec24(dmn),vec34(dmn),vec1o(dmn),      &
       vec2o(dmn),vec3o(dmn),vec15(dmn),vec73(dmn),vec76(dmn),vec78(dmn),      &
       vec7o(dmn)
    REAL(8),ALLOCATABLE :: N_full(:,:)
    LOGICAL :: p_in_dom, p_in_el
    c=0.125d0

    ! Type of the evaluation Obs or FD
    IF (strng=="ob") THEN
       typ=0
       nloop=nobs
       ALLOCATE(pick(nobs))
    ELSEIF (strng=="fd") THEN 
       typ=1
       nloop=ngp
       ALLOCATE(pick(nloop))
    END IF
    pick=0

    SELECT CASE(eltype) 
    CASE("tri"); ALLOCATE(nd_full(nloop,3),N_full(nloop,3))
    CASE("qua"); ALLOCATE(nd_full(nloop,4),N_full(nloop,4))
    CASE("tet"); ALLOCATE(nd_full(nloop,4),N_full(nloop,4))
    CASE("hex"); ALLOCATE(nd_full(nloop,8),N_full(nloop,8))
    END SELECT
    xmind=MINVAL(coords(:,1)); xmaxd=MAXVAL(coords(:,1))
    ymind=MINVAL(coords(:,2)); ymaxd=MAXVAL(coords(:,2))
    IF (dmn>2) THEN
       zmind=MINVAL(coords(:,3)); zmaxd=MAXVAL(coords(:,3))
    END IF

    DO ob=1,nloop ! Observation loop
       p_in_dom=.FALSE.
       IF (typ==0) THEN
          xob=ocoord(ob,:)*km2m
       ELSEIF (typ==1) THEN
          xob=xgp(ob,:)
       END IF
       p_in_dom=(xob(1)>=xmind .AND. xob(1)<=xmaxd .AND. xob(2)>=ymind .AND.   &
          xob(2)<=ymaxd)
       IF (dmn>2) p_in_dom=(p_in_dom .AND. xob(3)>=zmind .AND. xob(3)<=zmaxd)
       IF (p_in_dom) THEN ! Point probably in domain
          DO el=1,nels
             p_in_el=.FALSE.
             enodes=nodes(el,:)
             ecoords=coords(enodes,:)
             SELECT CASE(dmn)
             CASE(2)
                xmin=MINVAL(ecoords(:,1)); xmax=MAXVAL(ecoords(:,1))
                ymin=MINVAL(ecoords(:,2)); ymax=MAXVAL(ecoords(:,2))
                ! Point probably in 2D cell 
                IF (xob(1)>=xmin .AND. xob(1)<=xmax .AND. xob(2)>=ymin .AND.   &
                   xob(2)<=ymax) THEN
                   SELECT CASE(eltype)
                   CASE("tri")
                      vec12=ecoords(3,:)-ecoords(1,:)
                      vec13=ecoords(3,:)-ecoords(1,:)
                      vec23=ecoords(3,:)-ecoords(2,:)
                      vec1o=xob-ecoords(1,:)
                      vec2o=xob-ecoords(2,:)
                      vec3o=xob-ecoords(3,:)
                      CALL Mix2D(vec12,vec1o,dd)
                      CALL Mix2D(-vec13,vec3o,dl)
                      CALL Mix2D(vec23,vec2o,d)
                      ! Point in tri
                      IF (dd>=f0 .AND. dl>=f0 .AND. d>=f0) THEN
                         CALL Mix2D(-vec13,vec12,dr)
                         CALL Mix2D(vec12,vec13,du)
                         eta=dl/dr; nu=dd/du
                         N(1)=f1-eta-nu; N(2)=eta; N(3)=nu
                         p_in_el=.TRUE.
                      END IF
                   CASE("qua")
                     vec12=ecoords(2,:)-ecoords(1,:)
                     vec14=ecoords(4,:)-ecoords(1,:)
                     vec23=ecoords(3,:)-ecoords(2,:)
                     vec34=ecoords(4,:)-ecoords(3,:)
                     vec1o=xob-ecoords(1,:)
                     vec3o=xob-ecoords(3,:)
                     CALL Mix2D(vec12,vec1o,dd)
                     CALL Mix2D(-vec14,vec1o,dl)
                     CALL Mix2D(vec23,vec3o,dr)
                     CALL Mix2D(vec34,vec3o,du)
                     ! Point in quad
                     IF (dd>=f0 .AND. dl>=f0 .AND. dr>=f0 .AND. du>=f0) THEN
                        eta=(dl-dr)/(dl+dr); nu=(dd-du)/(du+dd)
                        N(1)=0.25d0*(f1-eta)*(f1-nu)
                        N(2)=0.25d0*(f1+eta)*(f1-nu)
                        N(3)=0.25d0*(f1+eta)*(f1+nu)
                        N(4)=0.25d0*(f1-eta)*(f1+nu)
                        p_in_el=.TRUE.
                     END IF
                   END SELECT ! Tri/qua
                END IF ! Point probably in 2D cell
             CASE(3) 
                xmin=MINVAL(ecoords(:,1)); xmax=MAXVAL(ecoords(:,1))
                ymin=MINVAL(ecoords(:,2)); ymax=MAXVAL(ecoords(:,2))
                zmin=MINVAL(ecoords(:,3)); zmax=MAXVAL(ecoords(:,3))
                ! Point probably in 3D cell
                IF (xob(1)>=xmin .AND. xob(1)<=xmax .AND. xob(2)>=ymin .AND.   &
                   xob(2)<=ymax .AND. xob(3)>=zmin .AND. xob(3)<=zmax) THEN
                   SELECT CASE(eltype)
                   CASE("tet")
                      vec12=ecoords(2,:)-ecoords(1,:)
                      vec13=ecoords(3,:)-ecoords(1,:)
                      vec14=ecoords(4,:)-ecoords(1,:)
                      vec23=ecoords(3,:)-ecoords(2,:)
                      vec24=ecoords(4,:)-ecoords(2,:)
                      vec1o=xob-ecoords(1,:)
                      vec2o=xob-ecoords(2,:)
                      CALL Mix3D(vec12,vec13,vec1o,dd)
                      CALL Mix3D(vec13,vec14,vec1o,dl)
                      CALL Mix3D(vec14,vec12,vec1o,df)
                      CALL Mix3D(vec24,vec23,vec2o,d)
                      ! Point in tet
                      IF (dd>=f0 .AND. dl>=f0 .AND. df>=f0 .AND. d>=f0) THEN 
                         CALL Mix3D(vec12,vec13,vec14,du)
                         CALL Mix3D(vec13,vec14,vec12,dr)
                         CALL Mix3D(vec14,vec12,vec13,db)
                         eta=dl/dr; nu=df/db; psi=dd/du
                         N(1)=f1-eta-nu-psi; N(2)=eta; N(3)=nu; N(4)=psi
                         p_in_el=.TRUE.
                      END IF
                   CASE("hex")
                      vec12=ecoords(2,:)-ecoords(1,:)
                      vec14=ecoords(4,:)-ecoords(1,:)
                      vec15=ecoords(5,:)-ecoords(1,:)
                      vec73=ecoords(3,:)-ecoords(7,:)
                      vec76=ecoords(6,:)-ecoords(7,:)
                      vec78=ecoords(8,:)-ecoords(7,:)
                      vec1o=xob-ecoords(1,:)
                      vec7o=xob-ecoords(7,:)
                      CALL Mix3D(vec12,vec14,vec1o,dd)
                      CALL Mix3D(vec15,vec12,vec1o,df)
                      CALL Mix3D(vec14,vec15,vec1o,dl)
                      CALL Mix3D(vec76,vec78,vec7o,du)
                      CALL Mix3D(vec78,vec73,vec7o,db)
                      CALL Mix3D(vec73,vec76,vec7o,dr)
                      ! Point in hex
                      IF (dd>=f0 .AND. dl>=f0 .AND. df>=f0 .AND. du>=f0 .AND.  &
                         dr>=f0 .AND. db>=f0) THEN
                         eta=(dl-dr)/(dr+dl); nu=(df-db)/(df+db)
                         psi=(dd-du)/(dd+du)
                         N(1)=c*(f1-eta)*(f1-nu)*(f1-psi)
                         N(2)=c*(f1+eta)*(f1-nu)*(f1-psi)
                         N(3)=c*(f1+eta)*(f1+nu)*(f1-psi)
                         N(4)=c*(f1-eta)*(f1+nu)*(f1-psi)
                         N(5)=c*(f1-eta)*(f1-nu)*(f1+psi)
                         N(6)=c*(f1+eta)*(f1-nu)*(f1+psi)
                         N(7)=c*(f1+eta)*(f1+nu)*(f1+psi)
                         N(8)=c*(f1-eta)*(f1+nu)*(f1+psi)
                         p_in_el=.TRUE.
                      END IF
                   END SELECT
                END IF ! Point probably in 3D cell
             END SELECT ! Dimension
             ! Update local pick list
             IF (p_in_el) THEN
                pick(ob)=ob
                N_full(ob,:)=N
                nd_full(ob,:)=enodes
                EXIT
             END IF
          END DO ! Element loop
       END IF ! Point probably in domain
    END DO ! Observation loop
    IF (typ==0) THEN
       nobs_loc=SIZE(PACK(pick,pick/=0))
       ALLOCATE(ol2g(nobs_loc),ocoord_loc(nobs_loc,dmn))
       SELECT CASE(eltype) 
       CASE("tri"); ALLOCATE(onlst(nobs_loc,3),oshape(nobs_loc,3))
       CASE("qua"); ALLOCATE(onlst(nobs_loc,4),oshape(nobs_loc,4))
       CASE("tet"); ALLOCATE(onlst(nobs_loc,4),oshape(nobs_loc,4))
       CASE("hex"); ALLOCATE(onlst(nobs_loc,8),oshape(nobs_loc,8))
       END SELECT
       ol2g=PACK(pick,pick/=0)
       ocoord_loc=ocoord(ol2g,:)
       onlst=nd_full(ol2g,:)
       oshape=N_full(ol2g,:)
    ELSEIF (typ==1) THEN
       ngp_loc=SIZE(PACK(pick,pick/=0))
       ALLOCATE(gpl2g(ngp_loc),idgp_loc(ngp_loc,dmn))
       SELECT CASE(eltype) 
       CASE("tri"); ALLOCATE(gpnlst(ngp_loc,3),gpshape(ngp_loc,3))
       CASE("qua"); ALLOCATE(gpnlst(ngp_loc,4),gpshape(ngp_loc,4))
       CASE("tet"); ALLOCATE(gpnlst(ngp_loc,4),gpshape(ngp_loc,4))
       CASE("hex"); ALLOCATE(gpnlst(ngp_loc,8),gpshape(ngp_loc,8))
       END SELECT
       gpl2g=PACK(pick,pick/=0)
       gpnlst=nd_full(gpl2g,:)
       gpshape=N_full(gpl2g,:)
       idgp_loc=idgp(gpl2g,:)
    END IF
  END SUBROUTINE GetObsNd

  ! Active fault node map loc->glb 
  SUBROUTINE GetFltMap
    IMPLICIT NONE
    INTEGER :: j,j3,map(n_lmnd,2)
    INTEGER, ALLOCATABLE :: rw(:)
    map=0
    DO j=1,n_lmnd
       IF (lmnd0+j>nceqs_ncf/(dmn+p)) THEN ! Has on rank fault nodes
          j3=lmnd0+j-nceqs_ncf/(dmn+p)
          IF (frc(j3)>0) THEN
             map(j,:)=(/j,j3/)
          END IF
       END IF
    END DO
    nfnd_loc=SIZE(PACK(map(:,1),map(:,1)/=0))
    ALLOCATE(rw(nfnd_loc),FltMap(nfnd_loc,2))
    rw=PACK(map(:,1),map(:,1)/=0)
    FltMap=map(rw,:)
  END SUBROUTINE GetFltMap

  ! Apply body force
  SUBROUTINE ApplyGravity
    IMPLICIT NONE
#include "petsc.h"
    REAL(8) :: dns,gip,gsca,gvec(npel*dmn),ecoords(npel,dmn),detj
    INTEGER :: el,i,j,indx((dmn+p)*npel),row(dmn*npel)
    DO el=1,nels
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       CALL FormElIndx(enodes,indx)
       row=indx(1:dmn*npel)
       row=indxmap(row,2)
       dns=mat(id(el),5)
       gsca=dns*gravity/DBLE(npel)
       gip=f0; gvec=f0
       DO i=1,nip
          CALL FormdetJ(ipoint(i,:),ecoords,detj)
          gip=gip+gsca*weight(i)*detj
       END DO
       ! Assume last dim aligned with gravity
       DO i=1,npel
          j=i*dmn
          IF (bc(nodes(el,i),dmn)/=0) gvec(j)=-gip
       END DO
       CALL VecSetValues(Vec_F,dmn*npel,row,gvec,Add_Values,ierr)
    END DO
    CALL VecAssemblyBegin(Vec_F,ierr)
    CALL VecAssemblyEnd(Vec_F,ierr)
  END SUBROUTINE ApplyGravity


!!$  ! Apply fluid body source
  SUBROUTINE ApplySource
    implicit none
#include "petsc.h"
    REAL(8) :: sdns,sip(npel,1),ssca,svec((dmn+1)*npel),ecoords(npel,dmn),detj,H,s
    REAL(8) :: dN(dmn,npel),gvec(dmn,1),Q(dmn,dmn)
    INTEGER :: el,i,j,indx((dmn+1)*npel),row((dmn+1)*npel)
    s=scale
    do el=1,nels
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       call FormElIndxp(enodes,indx)
       row=indx(1:(dmn+1)*npel)
       row=indxmap(row,2)
       H=mat(id(el),6)
       ! Mobility Tensor
       do i=1,dmn
        Q(i,i) = H
       end do        
       sdns=mat(id(el),10)
       ssca=sdns/DBLE(npel)
       sip=f0; svec=f0; gvec=f0
       gvec(dmn,1) = gravity
       do i=1,nip
          call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
          sip=sip+matmul(transpose(dN),matmul(Q,gvec))*sdns*weight(i)*detj*s
          !         sip=sip+ssca*weight(i)*detj
       end do
       do i=1,npel
          j=i*(dmn+1)
          IF (bc(nodes(el,i),dmn+1)/=0) svec(j)=sip(i,1)
       end do
       call VecSetValues(Vec_F,(dmn+1)*npel,row,svec,Add_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine ApplySource

  ! Reform RHS
  SUBROUTINE ReformLocalRHS(el,f,indx)
    IMPLICIT NONE
    INTEGER :: el,indx(:),j,j1
    REAL(8) :: f(:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    f=f0
    CALL ReformElRHS(ecoords,stress(el,:,:),E,nu,visc,expn,dt,f(1:eldof))
    CALL FormLocalIndx(enodes,indx)
    ! Fix BCs (i.e., zero out entries)
    DO j=1,npel
       DO j1=1,dmn
          IF (bc(nodes(el,j),j1)==0) f(dmn*j-dmn+j1)=f0
       END DO
    END DO
  END SUBROUTINE ReformLocalRHS

  ! Recover Vstress during implicit time stepping
  SUBROUTINE RecoverVStress(el,stress)
    IMPLICIT NONE
    INTEGER :: el,indxl(eldof)
    REAL(8) :: stress(:,:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    CALL FormLocalIndx(enodes,indx)
    indxl=indx(1:eldof)
    CALL CalcElVStress(ecoords,uu(indxl),stress(el,:,:),E,nu,visc,expn,dt)
  END SUBROUTINE RecoverVStress

  ! Account for Winkler foundation(s)
  SUBROUTINE AddWinklerFdn(el,k)
    IMPLICIT NONE
    INTEGER :: el,j,j1,j2,elbc(npel,dmn),snodes(nps)
    REAL(8) :: k(:,:),area
    dns=mat(id(el),5)
    elbc=bc(enodes,1:dmn)
    DO j1=1,dmn
       CALL GetWinklerEdge(elbc,j1,side)
       IF (side==0) CYCLE
       CALL EdgeAreaNodes(enodes,ecoords,side,area,snodes)
       val=dns*gravity*area/DBLE(nps)
       DO j=1,npel
          IF (elbc(j,j1)==-1) THEN
             j2=dmn*j-dmn+j1; k(j2,j2)=k(j2,j2)+val
          END IF
       END DO
    END DO
  END SUBROUTINE AddWinklerFdn

  ! Fix BCs (i.e., zero rows/columns) in local [K/Kp]
  SUBROUTINE FixBCinLocalK(el,k)
    IMPLICIT NONE
    INTEGER :: el,j,j1,j2
    REAL(8) :: k(:,:)
    ! Account for fixed displacement BC
    DO j=1,npel
       DO j1=1,dmn
          IF (bc(nodes(el,j),j1)==0) THEN
             j2=dmn*j-dmn+j1
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          END IF
       END DO
       IF (poro) THEN
          IF (bc(nodes(el,j),dmn+1)==0) THEN
             j2=dmn*npel+j
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          END IF
       END IF
    END DO
  END SUBROUTINE FixBCinLocalK

  ! Fix BCs (i.e., zero rows/columns) in local [K/Kp]
  ! Includes Dirichlet Pressure
  SUBROUTINE FixBCinLocalK_DP(el,k,f)
    IMPLICIT NONE
    INTEGER :: el,j,j1,j2
    REAL(8) :: k(:,:),f(:),rxp
    ! Account for fixed displacement BC
    DO j=1,npel
       DO j1=1,dmn
          IF (bc(nodes(el,j),j1)==0) THEN
             j2=dmn*j-dmn+j1
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          END IF
       END DO
       IF (poro) THEN
          IF (bc(nodes(el,j),dmn+1)==0) THEN
             rxp = rx_press(nodes(el,j))*scale
             j2=dmn*npel+j
             f(:) = f(:) - k(:,j2)*rxp
             f(j2) = rxp
             !val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             !k(j2,j2)=val
             k(j2,j2) = f1
          END IF
       END IF
    END DO
  END SUBROUTINE FixBCinLocalK_DP

  ! DP Print
  SUBROUTINE PrintMsg(msg)
    IMPLICIT NONE
    CHARACTER(*) :: msg
    IF (rank==0) PRINT*,msg
  END SUBROUTINE PrintMsg

  ! Write results in ASCII VTK (legacy) format
  SUBROUTINE WriteOutput
    IMPLICIT NONE
    CHARACTER(64) :: fmt
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0
    INTEGER :: i,j,j1,lnnds,lnels
    REAL(8),POINTER :: field_val(:)
    IF (dsp==0) field_val=>uu
    IF (dsp==1) field_val=>tot_uu
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,I0,3A,I0.6,A)')TRIM(name0),rank,"_",TRIM(name1),"_",k,".vtk"
    IF (dsp_dyn) THEN
       IF (dsp_hyb==1) THEN
          field_val=>tot_uu_dyn
       ELSE
          field_val=>uu_dyn
       END IF
       WRITE(name,'(A,I0,3A,I0.6,A)')TRIM(name0),rank,"_",TRIM(name1),"_dyn_", &
          k,".vtk"
    END IF
    OPEN(10,file=ADJUSTL(name),status='replace')
    lnnds=SIZE(coords,1)
    lnels=SIZE(nodes,1)
    WRITE(10,'(A)')"# vtk DataFile Version 2.0"
    WRITE(10,'(A)')"File written by Defmod-dev"
    WRITE(10,'(A)')"ASCII"
    WRITE(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    WRITE(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m),f0/)
       END DO
    CASE(3)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m)/)
       END DO
    END SELECT
    WRITE(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    SELECT CASE(npel)
    CASE(3); fmt="(I0,3(1X,I0))"
    CASE(4); fmt="(I0,4(1X,I0))"
    CASE(8); fmt="(I0,8(1X,I0))"
    END SELECT
    DO i=1,lnels
       WRITE(10,fmt)npel,nodes(i,:)-1
    END DO
    WRITE(10,'(A,I0)')"CELL_TYPES ",lnels
    DO i=1,lnels
       WRITE(10,'(I0)')vtkid
    END DO
    WRITE(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+p
    IF (dsp_dyn) j=dmn
    IF (poro .AND. .NOT. dsp_dyn) THEN
       WRITE(10,'(A)')"SCALARS pressure double"
       WRITE(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(scale*field_val(j1)) ! Pr
       END DO
    END IF
    WRITE(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          j1=i*j-p
          IF (dsp_dyn) j1=i*j
          WRITE(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       END DO
    CASE(3)
       DO i=1,lnnds
          j1=i*j-p
          IF (dsp_dyn) j1=i*j
          WRITE(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       END DO
    END SELECT
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput

  SUBROUTINE WriteOutput_x
    IMPLICIT NONE
    CHARACTER(64) :: fmt
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0
    INTEGER :: i,j,j1,j2,j3,lnnds,lnels
    REAL(8),POINTER :: field_val(:),field_val_fp(:),field_val_fl(:),           &
       field_val_qu(:),field_val_ql(:),field_val_flc(:),field_val_ss(:),       &
       field_val_sh(:)
    IF (dsp==0) field_val=>uu
    IF (dsp==1) field_val=>tot_uu
    IF (nceqs>0) field_val_fl=>fl; field_val_ql=>ql; field_val_flc=>flc
    IF (poro) field_val_fp=>fp; field_val_qu=>qu
    IF (fault .AND. lm_str==1) field_val_ss=>ss; field_val_sh=>sh
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,I0,3A,I0.6,A)')TRIM(name0),rank,"_",TRIM(name1),"_",k,".vtk"
    OPEN(10,file=ADJUSTL(name),status='replace')
    lnnds=SIZE(coords,1)
    lnels=SIZE(nodes,1)
    WRITE(10,'(A)')"# vtk DataFile Version 2.0"
    WRITE(10,'(A)')"File written by Defmod-dev"
    WRITE(10,'(A)')"ASCII"
    WRITE(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    WRITE(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m),f0/)
       END DO
    CASE(3)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m)/)
       END DO
    END SELECT
    WRITE(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    SELECT CASE(npel)
    CASE(3); fmt="(I0,3(1X,I0))"
    CASE(4); fmt="(I0,4(1X,I0))"
    CASE(8); fmt="(I0,8(1X,I0))"
    END SELECT
    DO i=1,lnels
       WRITE(10,fmt)npel,nodes(i,:)-1
    END DO
    WRITE(10,'(A,I0)')"CELL_TYPES ",lnels
    DO i=1,lnels
       WRITE(10,'(I0)')vtkid
    END DO
    WRITE(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+p; j2=dmn
    IF (poro) THEN
       WRITE(10,'(A)')"SCALARS pressure double"
       WRITE(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(scale*field_val(j1)) ! Pr
       END DO
       WRITE(10,'(A)')"SCALARS qu double"
       WRITE(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       DO i=1,lnnds
          WRITE(10,fmt)(field_val_qu(i)) ! qu
       END DO
       WRITE(10,'(A)')"VECTORS fp double"
       fmt="(3(F0.6,1X))"
       SELECT CASE(dmn)
       CASE(2)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_fp(j1-1),field_val_fp(j1),f0/) ! 2D U
          END DO
       CASE(3)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_fp(j1-2),field_val_fp(j1-1),             &
                field_val_fp(j1)/) ! 3D U
          END DO
       END SELECT
       IF (nceqs > 0) THEN
          WRITE(10,'(A)')"SCALARS ql double"
          WRITE(10,'(A)')"LOOKUP_TABLE default"
          fmt="(1(F0.3,1X))"
          DO i=1,lnnds
             WRITE(10,fmt)(field_val_ql(i)) ! ql
          END DO
       END IF
    END IF
    IF (fault .AND. lm_str==1) THEN
       WRITE(10,'(A)')"VECTORS ss double"
       fmt="(3(E10.4,1X))"
       SELECT CASE(dmn)
       CASE(2)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_ss(j1-1),field_val_ss(j1),f0/) ! 2D U
          END DO
       CASE(3)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_ss(j1-2),field_val_ss(j1-1),             &
                field_val_ss(j1)/) ! 3D U
          END DO
       END SELECT
       WRITE(10,'(A)')"VECTORS sh double"
       fmt="(3(E10.4,1X))"
       SELECT CASE(dmn)
       CASE(2)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_sh(j1-1),field_val_sh(j1),f0/) ! 2D U
          END DO
       CASE(3)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_sh(j1-2),field_val_sh(j1-1),             &
                field_val_sh(j1)/) ! 3D U
          END DO
       END SELECT
    END IF
    IF (nceqs>0) THEN
       DO j3=1,dmn
          IF (j3==1) THEN
             WRITE(10,'(A)')"SCALARS fls double"
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j3
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          ELSEIF(j3==2) THEN
             IF (dmn==2) THEN
                WRITE(10,'(A)')"SCALARS fln double"
             ELSE
                WRITE(10,'(A)')"SCALARS fld double"
             END IF
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j3
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          ELSEIF(j3==3) THEN
             WRITE(10,'(A)')"SCALARS fln double"
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j3
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          END IF
       END DO
       WRITE(10,'(A)')"VECTORS fl double"
       fmt="(3(F0.6,1X))"
       SELECT CASE(dmn)
       CASE(2)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_fl(j1-1),field_val_fl(j1),f0/) ! 2D U
          END DO
       CASE(3)
          DO i=1,lnnds
             j1=i*j2
             WRITE(10,fmt)(/field_val_fl(j1-2),field_val_fl(j1-1),             &
                field_val_fl(j1)/) ! 3D U
          END DO
       END SELECT
    END IF
    WRITE(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          j1=i*j-p
          WRITE(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       END DO
    CASE(3)
       DO i=1,lnnds
          j1=i*j-p
          WRITE(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       END DO
    END SELECT
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_x

  ! Write out deformation and (absolute) Coulomb force due 
  ! to pore pressure initialization
  SUBROUTINE WriteOutput_init
    IMPLICIT NONE
    CHARACTER(64) :: fmt
    CHARACTER(256) :: name,name0,name1
    INTEGER :: i,j,j1,j2,lnnds,lnels
    REAL(8),POINTER :: field_val(:),field_val_flc(:)
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,I0,A,A,A)')TRIM(name0),rank,"_",TRIM(name1),"_init.vtk"
    OPEN(10,file=ADJUSTL(name),status='replace')
    lnnds=SIZE(coords,1)
    lnels=SIZE(nodes,1)
    WRITE(10,'(A)')"# vtk DataFile Version 2.0"
    WRITE(10,'(A)')"File written by Defmod-dev"
    WRITE(10,'(A)')"ASCII"
    WRITE(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    WRITE(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m),f0/)
       END DO
    CASE(3)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m)/)
       END DO
    END SELECT
    WRITE(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    SELECT CASE(npel)
    CASE(3); fmt="(I0,3(1X,I0))"
    CASE(4); fmt="(I0,4(1X,I0))"
    CASE(8); fmt="(I0,8(1X,I0))"
    END SELECT
    DO i=1,lnels
       WRITE(10,fmt)npel,nodes(i,:)-1
    END DO
    WRITE(10,'(A,I0)')"CELL_TYPES ",lnels
    DO i=1,lnels
       WRITE(10,'(I0)')vtkid
    END DO
    WRITE(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+1
    field_val=>uu; field_val_flc=>flc
    WRITE(10,'(A)')"SCALARS pressure double"
    WRITE(10,'(A)')"LOOKUP_TABLE default"
    fmt="(1(F0.3,1X))"
    DO i=1,lnnds
       j1=i*j
       WRITE(10,fmt)(scale*field_val(j1)) ! Pr
    END DO
    IF (fault .AND. nfnd>0) THEN    
       DO j2=1,dmn
          IF (j2==1) THEN
             WRITE(10,'(A)')"SCALARS fls double"
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j2
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          ELSEIF(j2==2) THEN
             IF (dmn==2) THEN
                WRITE(10,'(A)')"SCALARS fln double"
             ELSE
                WRITE(10,'(A)')"SCALARS fld double"
             END IF
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j2
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          ELSEIF(j2==3) THEN
             WRITE(10,'(A)')"SCALARS fln double"
             WRITE(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             DO i=1,lnnds
                j1=i*dmn-dmn+j2
                WRITE(10,fmt)(field_val_flc(j1))
             END DO
          END IF
       END DO
    END IF
    WRITE(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          j1=i*j-p
          WRITE(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       END DO
    CASE(3)
       DO i=1,lnnds
          j1=i*j-p
          WRITE(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       END DO
    END SELECT
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_init

  ! Write out force to stress ratio on fault f2s
  SUBROUTINE WriteOutput_f2s
    IMPLICIT NONE
    CHARACTER(64) :: fmt
    CHARACTER(256) :: name,name0,name1
    INTEGER :: i,j,j1,lnnds,lnels
    REAL(8),POINTER :: field_val(:),field_val_dip(:),field_val_nrm(:)
    field_val=>f2s; field_val_dip=>dip; field_val_nrm=>nrm
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,I0,A,A,A)')TRIM(name0),rank,"_",TRIM(name1),"_f2s.vtk"
    OPEN(10,file=ADJUSTL(name),status='replace')
    lnnds=SIZE(coords,1)
    lnels=SIZE(nodes,1)
    WRITE(10,'(A)')"# vtk DataFile Version 2.0"
    WRITE(10,'(A)')"File written by Defmod-dev"
    WRITE(10,'(A)')"ASCII"
    WRITE(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    WRITE(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m),f0/)
       END DO
    CASE(3)
       DO i=1,lnnds
          WRITE(10,fmt)(/(coords(i,:)/km2m)/)
       END DO
    END SELECT
    WRITE(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    SELECT CASE(npel)
    CASE(3); fmt="(I0,3(1X,I0))"
    CASE(4); fmt="(I0,4(1X,I0))"
    CASE(8); fmt="(I0,8(1X,I0))"
    END SELECT
    DO i=1,lnels
       WRITE(10,fmt)npel,nodes(i,:)-1
    END DO
    WRITE(10,'(A,I0)')"CELL_TYPES ",lnels
    DO i=1,lnels
       WRITE(10,'(I0)')vtkid
    END DO
    WRITE(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn
    ! Write force to stress ratio
    IF (nceqs>0) THEN
       WRITE(10,'(A)')"VECTORS f2s double"
       fmt="(3(F0.6,1X))"
       SELECT CASE(dmn)
       CASE(2)
          DO i=1,lnnds
             j1=i*j
             WRITE(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
          END DO
       CASE(3)
          DO i=1,lnnds
             j1=i*j
             WRITE(10,fmt)(/field_val(j1-2),field_val(j1-1),                   &
                field_val(j1)/) ! 3D U
          END DO
       END SELECT
    END IF
    ! Write fault's dip (strike for 2D) and normal vectors
    SELECT CASE(dmn)
    CASE(2)
       WRITE(10,'(A)')"VECTORS strk double"
    CASE(3)
       WRITE(10,'(A)')"VECTORS dip double"
    END SELECT
    fmt="(3(F0.6,1X))"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(/field_val_dip(j1-1),field_val_dip(j1),f0/) ! 2D U
       END DO
    CASE(3)
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(/field_val_dip(j1-2),field_val_dip(j1-1),              &
             field_val_dip(j1)/) ! 3D U
       END DO
    END SELECT
    WRITE(10,'(A)')"VECTORS nrm double"
    SELECT CASE(dmn)
    CASE(2)
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(/field_val_nrm(j1-1),field_val_nrm(j1),f0/) ! 2D U
       END DO
    CASE(3)
       DO i=1,lnnds
          j1=i*j
          WRITE(10,fmt)(/field_val_nrm(j1-2),field_val_nrm(j1-1),              &
             field_val_nrm(j1)/) ! 3D U
       END DO
    END SELECT
    CLOSE(10)
  END SUBROUTINE WriteOutput_f2s

  ! Output observations
  SUBROUTINE WriteOutput_obs
    IMPLICIT NONE
    INTEGER :: i,j
    CHARACTER(64) :: fmt
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0,k_dyn=0
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    IF (dyn) THEN
       WRITE(name,'(A,A,A,I0.6,A)')TRIM(name0),TRIM(name1),"_",rank,           &
          "_dyn_obs.txt"
       j=k_dyn
    ELSE
       WRITE(name,'(A,A,A,I0.6,A)')TRIM(name0),TRIM(name1),"_",rank,"_obs.txt"
       j=k
    END IF
    IF (j==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       SELECT CASE(dmn)
       CASE(2); fmt="(2(F0.3,1X),I0)"
       CASE(3); fmt="(3(F0.3,1X),I0)"
       END SELECT
       WRITE(10,"(I0)")nobs_loc
       DO i=1,nobs_loc
          WRITE(10,fmt)ocoord_loc(i,:),ol2g(i)
       END DO
    ELSE
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    IF (dyn) THEN
       SELECT CASE(dmn)
          CASE(2); fmt="(2(F0.6,1X))"
          CASE(3); fmt="(3(F0.6,1X))"
       END SELECT
    ELSE
       SELECT CASE(dmn)
       CASE(2)
          IF (poro) THEN
             fmt="(3(F0.6,1X))"
          ELSE
             fmt="(2(F0.6,1X))"
          END IF
       CASE(3)
          IF (poro) THEN
             fmt="(4(F0.6,1X))"
          ELSE
             fmt="(3(F0.6,1X))"
          END IF
       END SELECT
    END IF
    DO i=1,nobs_loc
       IF (dyn) THEN
          IF (dsp_hyb==1) THEN
             WRITE(10,fmt)tot_uu_dyn_obs(i,:)
          ELSE
            WRITE(10,fmt)uu_dyn_obs(i,:)
          END IF
       ELSE
          IF (dsp==1) THEN
             WRITE(10,fmt)tot_uu_obs(i,:)
          ELSE
             WRITE(10,fmt)uu_obs(i,:)
          END IF
       END IF
    END DO
    CLOSE(10) 
    IF (dyn) THEN
       k_dyn=k_dyn+1
    ELSE
       k=k+1
    END IF
  END SUBROUTINE WriteOutput_obs

  ! Write event log file for quasi-static data
  SUBROUTINE WriteOutput_log
    IMPLICIT NONE
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0
    INTEGER :: seis
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,A,A)')TRIM(name0),TRIM(name1),".log"
    IF (k==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       WRITE(10,"(F0.6)")dt
    ELSE
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    IF (crp) THEN
       seis=0
    ELSE
       seis=1
    END IF
    WRITE(10,"(2(I0X))")n_log,seis
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_log

  ! Write event log file for seismic data
  SUBROUTINE WriteOutput_log_wave
    IMPLICIT NONE
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,A,A)')TRIM(name0),TRIM(name1),"_dyn.log"
    IF (k==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       IF (gf) THEN
          WRITE(10,"(F0.6)")dt*frq_wave
       ELSE
          WRITE(10,"(F0.6)")dt_dyn*frq_wave
       END IF
    ELSE
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    WRITE(10,"(I0)")n_log_wave
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_log_wave

  ! Write event log file for seismic data
  SUBROUTINE WriteOutput_log_slip
    IMPLICIT NONE
    CHARACTER(256) :: name,name0,name1
    INTEGER,SAVE :: k=0
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,A,A)')TRIM(name0),TRIM(name1),"_slip.log"
    IF (k==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       WRITE(10,"(F0.6)")dt_dyn*frq_slip
    ELSE
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    WRITE(10,"(I0)")n_log_slip
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_log_slip

  ! Write temporal fault slip (and theta for rate state friction) 
  SUBROUTINE WriteOutput_slip
    IMPLICIT NONE
#include "petsc.h"
    CHARACTER(256) :: name,name0,name1
    CHARACTER(64) :: fmt
    INTEGER :: j,j1,j2,j3,rw_loc(dmn) 
    INTEGER,SAVE :: k=0
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    IF (nfnd_loc>0) THEN ! Has on rank fault nodes
       WRITE(name,'(A,A,A,I0.6,A)')TRIM(name0),TRIM(name1),"_",rank,"_slip.txt"
       IF (rsf<1) THEN
          SELECT CASE(dmn)
             CASE(2); fmt="(2(ES11.2E3,X))"
             CASE(3); fmt="(3(ES11.2E3,X))"
          END SELECT
       ELSE
          SELECT CASE(dmn)
             CASE(2); fmt="(2(ES11.2E3,X),E12.4)"
             CASE(3); fmt="(3(ES11.2E3,X),E12.4)"
          END SELECT
       END IF
       IF (k==0) THEN
          OPEN(10,file=ADJUSTL(name),status='replace')
          WRITE(10,'(I0)')nfnd_loc!n_lmnd-max(0,(nceqs_ncf/(dmn+p)-lmnd0))
          DO j1=1,nfnd_loc
             j3=FltMap(j1,2) 
             SELECT CASE(dmn)
                CASE(2); WRITE(10,'(2(F0.3,X),I0)')xfnd(j3,:),j3
                CASE(3); WRITE(10,'(3(F0.3,X),I0)')xfnd(j3,:),j3
             END SELECT
          END DO
       ELSE
          OPEN(10,file=ADJUSTL(name),status='old',position='append',action=    &
             'write')
       END IF
       IF (rsf<1) THEN ! Slip weakening
          DO j1=1,nfnd_loc
             j=FltMap(j1,1); j3=FltMap(j1,2)
             rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
             IF (dsp_hyb==1) THEN
                 WRITE(10,fmt)tot_flt_slip(rw_loc(:dmn-1)),mu_hyb(j3)
             ELSE
                 WRITE(10,fmt)flt_slip(rw_loc(:dmn-1))/dt_dyn,mu_hyb(j3)
             END IF
          END DO
       ELSE ! Rate and state friction
          DO j1=1,nfnd_loc
             j=FltMap(j1,1); j3=FltMap(j1,2)
             rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
             IF (dsp_hyb==1) THEN 
                WRITE(10,fmt)tot_flt_slip(rw_loc(:dmn-1)),mu_hyb(j3),          &
                   rsftheta(j3)
             ELSE
                WRITE(10,fmt)flt_slip(rw_loc(:dmn-1))/dt_dyn,mu_hyb(j3),    &
                   rsftheta(j3)
             END IF
          END DO
       END IF ! Constitutive models
       CLOSE(10); k=k+1
    END IF ! Has on rank fault nodes
  END SUBROUTINE WriteOutput_slip

  ! Write qs fault data 
  SUBROUTINE WriteOutput_flt_qs
    IMPLICIT NONE
    CHARACTER(256) :: name,name0,name1
    CHARACTER(64) :: fmt
    INTEGER :: i
    INTEGER,SAVE :: k=0
    name0=output_file(:INDEX(output_file,"/",BACK=.TRUE.))
    name1=output_file(INDEX(output_file,"/",BACK=.TRUE.)+1:)
    WRITE(name,'(A,A,A)')TRIM(name0),TRIM(name1),"_fqs.txt"
    SELECT CASE(dmn+p)
       CASE(2); fmt="(2(F0.6,1X))"
       CASE(3); fmt="(3(F0.6,1X))"
       CASE(4); fmt="(4(F0.6,1X))"
    END SELECT
    IF (k==0) THEN
       OPEN(10,file=ADJUSTL(name),status='replace')
       WRITE(10,'(I0)')SUM(frc)
       DO i=1,nfnd
          IF (frc(i)>0) THEN
             SELECT CASE(dmn)
                CASE(2); WRITE(10,'(2(F0.3,1X))')xfnd(i,:)
                CASE(3); WRITE(10,'(3(F0.3,1X))')xfnd(i,:)
             END SELECT
          END IF
       END DO
    ELSE
       OPEN(10,file=ADJUSTL(name),status='old',position='append',action='write')
    END IF
    DO i=1,nfnd
       IF (frc(i)>0) THEN
          IF (poro) THEN
             WRITE(10,fmt)flt_ss(i,:),flt_p(i)
          ELSE
             WRITE(10,fmt)flt_ss(i,:)
          END IF
       END IF
    END DO
    CLOSE(10); k=k+1
  END SUBROUTINE WriteOutput_flt_qs

END MODULE global
