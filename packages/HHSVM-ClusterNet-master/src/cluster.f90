SUBROUTINE standard(nobs,nvars,x,ju,isd,xmean,xnorm,maj)     
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER::  nobs
    INTEGER::nvars
    INTEGER::isd
    INTEGER::ju(nvars)
    DOUBLE PRECISION::  x(nobs,nvars)
    DOUBLE PRECISION::xmean(nvars)
    DOUBLE PRECISION::xnorm(nvars)
    DOUBLE PRECISION::maj(nvars)
    ! - - - local declarations - - -
    INTEGER:: j
! - - - begin - - -                                
    DO j=1,nvars                                  
        IF(ju(j)==1) THEN                         
            xmean(j)=sum(x(:,j))/nobs     !mean                        
            x(:,j)=x(:,j)-xmean(j)    
            maj(j)=dot_product(x(:,j),x(:,j))/nobs                                              
              IF(isd==1) THEN
                xnorm(j)=sqrt(maj(j))    !standard deviation               
                x(:,j)=x(:,j)/xnorm(j)
                maj(j)=1.0D0
            ENDIF                                                        
        ENDIF                                     
    END DO                             
END SUBROUTINE standard
!----------------Derive of loss function

 
      SUBROUTINE DWDdrv (nobs, nvars, x, y, r, vl,delta)
         IMPLICIT NONE
         INTEGER :: nobs, nvars, i
         DOUBLE PRECISION :: x (nobs, nvars), y (nobs)
         DOUBLE PRECISION :: r (nobs), vl (nvars), dl (nobs)
            DOUBLE PRECISION :: delta
         vl = 0.0
         DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)- 1.0D0) / delta
                  END IF
               END DO
         vl = Matmul(dl * y, x) / nobs
      END SUBROUTINE DWDdrv
!----calculate nw et ws
 
! --------------------------------------------------
SUBROUTINE chkvars (nobs, nvars, x, ju)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: ju (nvars)
      DOUBLE PRECISION :: x (nobs, nvars)
    ! - - - local declarations - - -
      INTEGER :: i
      INTEGER :: j
      DOUBLE PRECISION :: t
! - - - begin - - -
      DO j = 1, nvars
         ju (j) = 0
         t = x (1, j)
         DO i = 2, nobs
            IF (x(i, j) /= t) THEN
               ju (j) = 1
               EXIT
            END IF
         END DO
      END DO
END SUBROUTINE chkvars
! --------------------------------------------------
SUBROUTINE hsvmlassoNETc (delta, lam2, nobs, nvars, x, y, KK, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr,cluster,strong )
! --------------------------------------------------
    
      IMPLICIT NONE
    ! - - - arg types - - -
      LOGICAL :: strong
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: KK
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
       DOUBLE PRECISION :: cluster (nvars, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
     !-------kmeans variables
     
    
      INTEGER :: IC2(nvars), NC(KK), NCP(KK), ITRAN(KK), LIVE(KK)
      DOUBLE PRECISION ::  ND(nvars), AN1(KK), AN2(KK), WSS(KK),CENTER(KK,nobs)
      
      
     
      
! - - - begin - - -
! - - - allocate variables - - -
  
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
     
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0
      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
       
      CALL hsvmlassoNETpathc (delta, lam2, maj, nobs, nvars, x, y,KK , ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr, cluster,IC2, NC, NCP, ITRAN, LIVE,CENTER,ND, AN1, AN2, WSS,strong )
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
        IF (isd == 1) THEN
           DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
         b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
        & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE hsvmlassoNETc
! --------------------------------------------------
SUBROUTINE hsvmlassoNETpathc (delta, lam2, maj, nobs, nvars, x, y, KK , ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr, cluster,IC2, NC, NCP, ITRAN, LIVE,CENTER,ND, AN1, AN2, WSS,strong)
! --------------------------------------------------
    
    


      IMPLICIT NONE
        ! - - - arg types - - -
    
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: KK
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: clust (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: A (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: vl (nvars), ga(nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: cluster (nvars, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars),tlam
      
    
!!!!-------parameters for kmeans
      INTEGER ::ITER 
      INTEGER :: IFAULT 
      INTEGER :: IC2(nvars), NC(KK), NCP(KK), ITRAN(KK), LIVE(KK)
      DOUBLE PRECISION ::  ND(nvars), AN1(KK), AN2(KK), WSS(KK),CENTER(KK,nobs)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: ws
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al,al0
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
     
   INTEGER :: jx, jxx (nvars)
      
     
      INTEGER :: i
      INTEGER :: cI
      INTEGER :: k
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: l
      INTEGER :: h
      INTEGER :: nw
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER :: COUNTER 
      LOGICAL :: strong
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
    
     
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
     
     
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      
      r = 0.0D0
       
      b = 0.0D0
      oldbeta = rand(0)
      m = 0
      mm = 0
    
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam - 1.0D0))
       
      END IF
! strong rule
!   jxx = 0 using the strong rule
!   jxx = 1 not using the strong rule
         jxx = 0
         IF (strong .EQV. .TRUE.) THEN
            jxx = 0
         ELSE
            jxx = 1
         END IF

! --------- lambda loop ----------------------------
     
      DO l = 1, nlam
      
     IF (me < KK) THEN

                   CALL  kmeans(TRANSPOSE(x), nvars, nobs, CENTER, KK, clust, IC2, NC, AN1, AN2, NCP, ND,&
                     & ITRAN, LIVE, ITER, WSS, IFAULT)
       
              ELSE 
                   DO j = 1,nvars
                       A(:,j) = x(:,j) * b(j)
                   END DO  
                          CALL kmeans(TRANSPOSE(A), nvars, nobs, CENTER, KK, clust, IC2, NC, AN1, AN2, NCP, ND,&
                          & ITRAN, LIVE, ITER, WSS, IFAULT)
       
           END IF

          
         
         
               
           al0 = al
        
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
           
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al0 = 0.0D0
          CALL DWDdrv (nobs, nvars, x, y, r, vl,delta)
                ga = Abs(vl)
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN                      
                        al0 = Max (al0, ga(j) /pf(j))
                       
                     END IF
                  END IF
               END DO
               al = al0 * alf 
                
            END IF
         END IF
          ctr = 0
!---------- check strong rule (in lamdba loop) ------------------
            tlam = 2.0 * al - al0
            loop_strong_rule: DO j = 1, nvars
               IF (jxx(j) == 1) CYCLE
               IF (ga(j) > pf(j) * tlam) jxx(j) = 1
            END DO loop_strong_rule

         


        ! ---------kmeans + GCD----------------------------
        
     DO  
       
       
          oldbeta  = b 
             
            
    ! --------- outer loop ---------------------------- 
          
        !  Nsort = 0 
         DO
         !  Nsort = Nsort + 1
          !  IF (Nsort > 10) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
           
             Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 10) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                        u = u + dl (i) * y (i) * x (i, k)
                   
                END DO
                  
                  
             ! --somme rho jl beta l-------------------------------------                     
                
            
                  
                  nw=1
                  ws = 0.0D0
                  DO h = 1, nvars
                        
                     IF ((clust(k)==clust(h)).AND.(h/=k)) THEN
                        
                        w = 0.0D0
                        DO cI = 1, nobs
                           w = w +x (cI, k)*x (cI, h)
                       END DO
                        ws = ws + w*b (h)/nobs
                         nw = nw + 1       
                     END IF
                     
                  END DO
                  
  !----------------------------------------------------------- 
                    
                     u = maj (k) * b (k) - u / nobs + pf2(k)*lam2*ws/nw
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                     	b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2*(nw-1.0)/nw)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + y * x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
            END DO
               IF (ni > pmax) EXIT
               d = 0.0D0
               DO i = 1, nobs
               IF (r(i) > 1.0D0) THEN
                  dl (i) = 0.0D0
               ELSE IF (r(i) <= (1-delta)) THEN
                  dl (i) = - 1.0D0
               ELSE
                  dl (i) = (r(i)-1.0D0) / delta
               END IF
                  d = d + dl (i) * y (i)
               END DO
               d = - 0.5D0 * delta * d / nobs
               IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r + y * d
                  dif = Max (dif, 2.0*d**2/delta)
              
               END IF
               IF (dif < eps) EXIT
        ! --inner loop----------------------
               ! NiSort = 0
               DO
                !  NiSort = NiSort + 1
                 ! IF (NiSort > 10) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                        u = u + dl (i) * y (i) * x (i, k)
                  
                     END DO
                   
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
                     u = maj (k) * b (k) - u / nobs + pf2(k) * lam2 * ws / nw
                   
                     v = al * pf (k)
                     v = Abs (u) - v
                    
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2 * (nw-1.0)/nw)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + y * x (:, k) * d
                     END IF
                  END DO
                  d = 0.0D0
                  DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                     d = d + dl (i) * y (i)
                  END DO
                  d = - 0.5D0 * delta * d / nobs
                  IF (d /= 0.0D0) THEN
                     b (0) = b (0) + d
                     r = r + y * d
                     dif = Max (dif, 2.0*d**2/delta)

                  END IF
                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF (2.0*(b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF (2.0*(b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
           if (vrg /= 1) cycle
           CALL DWDdrv (nobs, nvars, x, y, r, vl,delta)
               ga = Abs (vl)
               DO j = 1, nvars
                  IF (jxx(j) == 1) CYCLE
                  If (ga(j) > al * pf(j)) THEN
                     jxx(j) = 1
                     vrg = 0
                  END IF
               END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
            
         END DO



!!!!!the second final chek for our algorithm
           vrg = 1
            IF (2.0*(b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF (2.0*(b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
           if (vrg /= 1) cycle
           CALL DWDdrv (nobs, nvars, x, y, r, vl,delta)
               ga = Abs (vl)
               DO j = 1, nvars
                  IF (jxx(j) == 1) CYCLE
                  If (ga(j) > al * pf(j)) THEN
                     jxx(j) = 1
                     vrg = 0
                  END IF
               END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
            
   END DO 

! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
            
         IF (ni > 0) beta (1:ni, l) = b(m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         cluster(1:nvars, l) = clust
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
      
         IF (me > dfmax) EXIT
 
               
         
    END DO
   
      DEALLOCATE (b, oldbeta, r, mm  )
  
      RETURN
    

END SUBROUTINE hsvmlassoNETpathc
 
!--------------------------------------------------------------------------------------
       
     !--------------------------------------------------------------------------------------
       
     SUBROUTINE sample(m,k,b)

               integer::m
               integer::k
               integer::i,bo,j,h         
               integer::b(k)
               DOUBLE PRECISION::a
             

                i = 1
                do while   (i<=k)
                    call random_number(a)
                    if (i == 1) then 
                           b(i) = floor(a*m) + 1 
                           i = i+1
                    else                         
                       bo = floor(a*m) + 1
                       h=1
                       do  j=1,(i-1)
                           if (bo==b(j)) then 
                                  h=0
                           end if
                       end do
                       if (h==1) then 
                             b(i) = bo
                             i = i+1 
                       end if
                    end if
                     
                end do
             
              RETURN
        END   SUBROUTINE sample
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kmeans(A, M, N, CENTER ,  K, IC1, IC2, NC, AN1, AN2, NCP, D,&
       & ITRAN, LIVE, ITER, WSS, IFAULT)

            INTEGER :: M,N,K,ITER ,IFAULT,jer,l,ll
             INTEGER,PARAMETER::nstart1=50
              INTEGER,PARAMETER::nstart2=50
            INTEGER :: IC(M), IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
            DOUBLE PRECISION :: A(M,N), D(M),  AN1(K), AN2(K), WSS(K),sumwss(nstart1),CENTER(K,N),W(K)
            DOUBLE PRECISION :: best,z    
            INTEGER ,DIMENSION(:),allocatable ::b
            
            DOUBLE PRECISION ,DIMENSION(:, :),allocatable ::C
            ALLOCATE(b(K),stat=jer)
          
            ALLOCATE (C(K,N), STAT=jer)
       DO ll = 1,nstart2
         DO l = 1,nstart1 
            Call sample(M,K,b)
              DO  i=1,K
               DO  j=1,N
                  C(i,j)=A(b(i),j)
              END DO
              END DO
            
            CALL KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
                     & ITRAN, LIVE, ITER, WSS, IFAULT)
              W = WSS
            sumwss(l) = sum(W)
         END DO 
           best = MINVAL(sumwss) 
         
           
                Call sample(M,K,b)
              DO  i=1,K
               DO  j=1,N
                  C(i,j)=A(b(i),j)
              END DO
              END DO
                 CALL KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
                     & ITRAN, LIVE, ITER, WSS, IFAULT)
                  IC = IC1
                  W = WSS
                  z = sum(W)
                  
           
         
             if(z .LE. best) exit
                IC1 = IC
             
       
                      
          END DO    
        
      DEALLOCATE (b , C ) 
       return
  END subroutine kmeans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
        & ITRAN, LIVE, ITER, WSS, IFAULT)
      INTEGER ::ITER 
      INTEGER M,N,K ,IFAULT 
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      DOUBLE PRECISION A(M,N), D(M), C(K,N), AN1(K), AN2(K), WSS(K)

      DOUBLE PRECISION DT(2), ZERO, ONE
      INTEGER I,IL,J,L,INDX,IJ,II
      DOUBLE PRECISION BIG, DA, TEMP, DB, DC,AA

      DATA BIG /1.E30/, ZERO /0.0/, ONE /1.0/

      IFAULT = 3

      IF (K .LE. 1 .OR. K .GE. M) RETURN
      IFAULT = 0

      DO 50 I = 1, M
        IC1(I) = 1
        IC2(I) = 2
        DO 10 IL = 1, 2
          DT(IL) = ZERO
          DO 10 J = 1, N
            DA = A(I,J) - C(IL,J)
            DT(IL) = DT(IL) + DA*DA
   10   CONTINUE
        IF (DT(1) .GT. DT(2)) THEN
          IC1(I) = 2
          IC2(I) = 1
          TEMP = DT(1)
          DT(1) = DT(2)
          DT(2) = TEMP
        END IF
        DO 50 L = 3, K
          DB = ZERO
          DO 30 J = 1, N
            DC = A(I,J) - C(L,J)
            DB = DB + DC*DC
            IF (DB .GE. DT(2)) GO TO 50
   30     CONTINUE
          IF (DB .LT. DT(1)) GO TO 40
          DT(2) = DB
          IC2(I) = L
          GO TO 50
   40     DT(2) = DT(1)
          IC2(I) = IC1(I)
          DT(1) = DB
          IC1(I) = L
   50 CONTINUE

      DO 70 L = 1, K
        NC(L) = 0
        DO 60 J = 1, N
   60   C(L,J) = ZERO
   70 CONTINUE
      DO 90 I = 1, M
        L = IC1(I)
        NC(L) = NC(L) + 1
        DO 80 J = 1, N
   80   C(L,J) = C(L,J) + A(I,J)
   90 CONTINUE

      DO 120 L = 1, K
        IF (NC(L) .EQ. 0) THEN
          IFAULT = 1
          RETURN
        END IF
        AA = NC(L)
        DO 110 J = 1, N
  110   C(L,J) = C(L,J) / AA

        AN2(L) = AA / (AA + ONE)
        AN1(L) = BIG
        IF (AA .GT. ONE) AN1(L) = AA / (AA - ONE)
        ITRAN(L) = 1
        NCP(L) = -1
  120 CONTINUE
      INDX = 0
      DO 140 IJ = 1, ITER

        CALL OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &        ITRAN, LIVE, INDX)

        IF (INDX .EQ. M) GO TO 150

        CALL QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &       ITRAN, INDX)

        IF (K .EQ. 2) GO TO 150

        DO 130 L = 1, K
  130   NCP(L) = 0
  140 CONTINUE

      IFAULT = 2

  150 DO 160 L = 1, K
        WSS(L) = ZERO
        DO 160 J = 1, N
          C(L,J) = ZERO
  160 CONTINUE
      DO 170 I = 1, M
        II = IC1(I)
        DO 170 J = 1, N
          C(II,J) = C(II,J) + A(I,J)
  170 CONTINUE
      DO 190 J = 1, N
        DO 180 L = 1, K
  180   C(L,J) = C(L,J) / FLOAT(NC(L))
        DO 190 I = 1, M
          II = IC1(I)
          DA = A(I,J) - C(II,J)
          WSS(II) = WSS(II) + DA*DA
  190 CONTINUE

      RETURN
      END
!---------------------------------------------------------------------------------------
      SUBROUTINE OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &      ITRAN, LIVE, INDX)

      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      DOUBLE PRECISION    A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      INTEGER L,I,L1,L2,LL,J
      DOUBLE PRECISION ZERO, ONE
      DOUBLE PRECISION BIG,DE,DF,DA,DB,R2,RR,DC,DD,AL1,ALW,AL2,ALT

      DATA BIG /1.0E30/, ZERO /0.0/, ONE/1.0/

      DO 10 L = 1, K
        IF (ITRAN(L) .EQ. 1) LIVE(L) = M + 1
   10 CONTINUE
      DO 100 I = 1, M
        INDX = INDX + 1
        L1 = IC1(I)
        L2 = IC2(I)
        LL = L2

        IF (NC(L1) .EQ. 1) GO TO 90

        IF (NCP(L1) .EQ. 0) GO TO 30
        DE = ZERO
        DO 20 J = 1, N
          DF = A(I,J) - C(L1,J)
          DE = DE + DF*DF
   20   CONTINUE
        D(I) = DE * AN1(L1)

   30   DA = ZERO
        DO 40 J = 1, N
          DB = A(I,J) - C(L2,J)
          DA = DA + DB*DB
   40   CONTINUE
        R2 = DA * AN2(L2)
        DO 60 L = 1, K

          IF (I .GE. LIVE(L1) .AND. I .GE. LIVE(L) .OR. L .EQ. L1 .OR.&
          &   L .EQ. LL) GO TO 60
          RR = R2 / AN2(L)
          DC = ZERO
          DO 50 J = 1, N
            DD = A(I,J) - C(L,J)
            DC = DC + DD*DD
            IF (DC .GE. RR) GO TO 60
   50     CONTINUE
          R2 = DC * AN2(L)
          L2 = L
   60     CONTINUE
          IF (R2 .LT. D(I)) GO TO 70

          IC2(I) = L2
          GO TO 90

   70     INDX = 0
          LIVE(L1) = M + I
          LIVE(L2) = M + I
          NCP(L1) = I
          NCP(L2) = I
          AL1 = NC(L1)
          ALW = AL1 - ONE
          AL2 = NC(L2)
          ALT = AL2 + ONE
          DO 80 J = 1, N
            C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
            C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   80     CONTINUE
          NC(L1) = NC(L1) - 1
          NC(L2) = NC(L2) + 1
          AN2(L1) = ALW / AL1
          AN1(L1) = BIG
          IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
          AN1(L2) = ALT / AL2
          AN2(L2) = ALT / (ALT + ONE)
          IC1(I) = L2
          IC2(I) = L1
   90   CONTINUE
        IF (INDX .EQ. M) RETURN
  100 CONTINUE
      DO 110 L = 1, K

        ITRAN(L) = 0
        LIVE(L) = LIVE(L) - M
  110 CONTINUE

      RETURN
      END
!--------------------------------------------------------------------------------------------
      SUBROUTINE QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &    ITRAN, INDX)

      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K)
      DOUBLE PRECISION A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      DOUBLE PRECISION ZERO, ONE
      INTEGER ICOUN,ISTEP,I,L1,L2,J
      DOUBLE PRECISION BIG,DA,DB,DD,AL1,ALW,AL2,ALT,R2,DE

      DATA BIG /1.0E30/, ZERO /0.0/, ONE /1.0/

      ICOUN = 0
      ISTEP = 0
   10 DO 70 I = 1, M
        ICOUN = ICOUN + 1
        ISTEP = ISTEP + 1
        L1 = IC1(I)
        L2 = IC2(I)

        IF (NC(L1) .EQ. 1) GO TO 60

        IF (ISTEP .GT. NCP(L1)) GO TO 30
        DA = ZERO
        DO 20 J = 1, N
          DB = A(I,J) - C(L1,J)
          DA = DA + DB*DB
   20   CONTINUE
        D(I) = DA * AN1(L1)

   30   IF (ISTEP .GE. NCP(L1) .AND. ISTEP .GE. NCP(L2)) GO TO 60
        R2 = D(I) / AN2(L2)
        DD = ZERO
        DO 40 J = 1, N
          DE = A(I,J) - C(L2,J)
          DD = DD + DE*DE
          IF (DD .GE. R2) GO TO 60
   40   CONTINUE

        ICOUN = 0
        INDX = 0
        ITRAN(L1) = 1
        ITRAN(L2) = 1
        NCP(L1) = ISTEP + M
        NCP(L2) = ISTEP + M
        AL1 = NC(L1)
        ALW = AL1 - ONE
        AL2 = NC(L2)
        ALT = AL2 + ONE
        DO 50 J = 1, N
          C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
          C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   50   CONTINUE
        NC(L1) = NC(L1) - 1
        NC(L2) = NC(L2) + 1
        AN2(L1) = ALW / AL1
        AN1(L1) = BIG
        IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
        AN1(L2) = ALT / AL2
        AN2(L2) = ALT / (ALT + ONE)
        IC1(I) = L2
        IC2(I) = L1

   60   IF (ICOUN .EQ. M) RETURN
   70 CONTINUE
      GO TO 10
      END

!----------------------------------------------------------------------------------
 

