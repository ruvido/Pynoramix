MODULE aux_funcs

CONTAINS

SUBROUTINE dist (pbc_opt,coors1,box1,coors2,n1,n2,matrix)

IMPLICIT NONE

INTEGER,INTENT(IN)::pbc_opt
integer,intent(in)::n1,n2
real,dimension(n1,3),intent(in)::coors1
REAL,DIMENSION(3,3),INTENT(IN)::box1
real,dimension(n2,3),intent(in)::coors2
real,dimension(n1,n2),intent(out)::matrix
integer::i,j
real,dimension(3)::vect


IF (pbc_opt==1) THEN

   vect=0.0d0
   do i=1,n1
      do j=1,n2
         
         vect=(coors1(i,:)-coors2(j,:))
         CALL PBC (vect,box1)
         matrix(i,j)=sqrt(dot_product(vect,vect))
         
      end do
   end do
   
ELSE

   vect=0.0d0
   do i=1,n1
      do j=1,n2
         vect=(coors1(i,:)-coors2(j,:))
         matrix(i,j)=sqrt(dot_product(vect,vect))
      end do
   end do
   
END IF

END SUBROUTINE dist


SUBROUTINE neighbs_limit(pbc_opt,ident,limit,coors1,box1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  INTEGER,INTENT(IN)::limit
  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1
  INTEGER,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::i,j,g
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.


  ALLOCATE(list(limit),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.true.

  DO i=1,n_atoms1
     aux=coors1(i,:)
     DO j=1,n_atoms2
        aux2=coors2(j,:)-aux
        IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
        list_dists(j)=sqrt(dot_product(aux2,aux2))
        list_vects(j,:)=aux2
     END DO

     IF (lident.eqv..true.) filter(i)=.false.

     DO j=1,limit
        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
        list(j)=g
        norm=list_dists(g)
        neighb_dist(i,j)=norm
        neighb_uvect(i,j,:)=list_vects(g,:)/norm
        neighb_list(i,j)=g
        filter(g)=.false.
     END DO

     DO j=1,limit
        g=list(j)
        filter(g)=.true.
     END DO
     filter(i)=.true.
  END DO

  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1

END SUBROUTINE NEIGHBS_LIMIT

SUBROUTINE neighbs_dist2(pbc_opt,ident,ii,dist,coors1,box1,coors2,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  REAL,INTENT(IN)::dist
  INTEGER,INTENT(IN)::n_atoms2,ii
  REAL,DIMENSION(3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1

  INTEGER,DIMENSION(n_atoms2),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms2),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms2,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::j,g,limit
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.

  

  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.false.


  aux=coors1(:)
  limit=0

  DO j=1,n_atoms2
     aux2=coors2(j,:)-aux
     IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
     norm=sqrt(dot_product(aux2,aux2))
     IF (norm<=dist) THEN
        limit=limit+1
        filter(j)=.true.
        list_dists(j)=norm
        list_vects(j,:)=aux2
     END IF
  END DO
  
  IF (lident.eqv..true.) THEN 
     filter(ii)=.false.
     limit=limit-1
  END IF
  
  
  DO j=1,limit
     g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
     list(j)=g
     norm=list_dists(g)
     neighb_dist(j)=norm
     neighb_uvect(j,:)=list_vects(g,:)/norm
     neighb_list(j)=g
     filter(g)=.false.
  END DO
  
  DO j=1,limit
     g=list(j)
     filter(g)=.false.
  END DO
  filter(ii)=.false.


  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1


END SUBROUTINE NEIGHBS_DIST2


SUBROUTINE neighbs_dist(pbc_opt,ident,dist,coors1,box1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  REAL,INTENT(IN)::dist
  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1
  INTEGER,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms1,n_atoms2,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::i,j,g,limit
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.


  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.false.

  DO i=1,n_atoms1
     aux=coors1(i,:)
     limit=0
     DO j=1,n_atoms2
        aux2=coors2(j,:)-aux
        IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
        norm=sqrt(dot_product(aux2,aux2))
        IF (norm<=dist) THEN
           limit=limit+1
           filter(j)=.true.
           list_dists(j)=norm
           list_vects(j,:)=aux2
        END IF
     END DO

     IF (lident.eqv..true.) THEN 
        filter(i)=.false.
        limit=limit-1
     END IF

     DO j=1,limit
        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
        list(j)=g
        norm=list_dists(g)
        neighb_dist(i,j)=norm
        neighb_uvect(i,j,:)=list_vects(g,:)/norm
        neighb_list(i,j)=g
        filter(g)=.false.
     END DO

     DO j=1,limit
        g=list(j)
        filter(g)=.false.
     END DO
     filter(i)=.false.
  END DO

  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1

END SUBROUTINE NEIGHBS_DIST


SUBROUTINE PBC(vector,box)

  IMPLICIT NONE

  REAL,DIMENSION(3),INTENT(INOUT)::vector
  REAL,DIMENSION(3,3),INTENT(IN)::box
  INTEGER::i
  REAL::x,L,Lhalf

  DO i=1,3
     L=box(i,i)
     Lhalf=0.50d0*L
     x=vector(i)
     IF (abs(x)>Lhalf) THEN
        IF (x>Lhalf) THEN
           x=x-L
        ELSE
           x=x+L
        END IF
        vector(i)=x
     END IF
  END DO
  
END SUBROUTINE PBC

SUBROUTINE SYMMETRIZE_NET(newKmax,TT_tau,TT_ind,TT_start,Pe,newKtot,T_ind,T_tau,T_start,N_nodes,Ktot)

  IMPLICIT NONE

  TYPE array_pointer
     INTEGER,DIMENSION(:),POINTER::p1
  END TYPE array_pointer

  INTEGER,INTENT(IN)::N_nodes,Ktot,newKtot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind,T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::Pe
  INTEGER,DIMENSION(N_nodes+1),INTENT(OUT)::TT_start
  INTEGER,DIMENSION(newKtot),INTENT(OUT)::TT_ind,TT_tau
  INTEGER,INTENT(OUT)::newKmax

  INTEGER::i,j,l,h,g,gg,destino
  integer,dimension(:),allocatable::salidas,auxx1,auxx2,SL
  TYPE(array_pointer),DIMENSION(:),POINTER::F_ind,flux
  LOGICAL::interruptor
  


  Pe=0
  TT_start=0
  TT_ind=0
  TT_tau=0
  newKmax=0


  ALLOCATE(F_ind(N_nodes),flux(N_nodes),salidas(N_nodes),SL(N_nodes))
  salidas=0
  SL=0

  DO i=1,N_nodes
     DO j=T_start(i)+1,T_start(i+1)

        IF (T_ind(j)==i) THEN
           SL(i)=T_tau(j)*2
        ELSE
           
           destino=T_ind(j)

           gg=salidas(i)
           IF (gg==0) THEN
              ALLOCATE(F_ind(i)%p1(1),flux(i)%p1(1))
              F_ind(i)%p1(1)=destino
              flux(i)%p1(1)=T_tau(j)
              salidas(i)=1
           ELSE
              interruptor=.false.
              DO h=1,gg
                 IF (F_ind(i)%p1(h)==destino) THEN
                    flux(i)%p1(h)=flux(i)%p1(h)+T_tau(j)
                    interruptor=.true.
                    exit
                 END IF
              END DO
              IF (interruptor.eqv..false.) THEN
                 ALLOCATE(auxx1(gg+1),auxx2(gg+1))
                 auxx1(1:gg)=F_ind(i)%p1(:)
                 auxx2(1:gg)=flux(i)%p1(:)
                 auxx1(gg+1)=destino
                 auxx2(gg+1)=T_tau(j)
                 salidas(i)=gg+1
                 DEALLOCATE(F_ind(i)%p1,flux(i)%p1)
                 ALLOCATE(F_ind(i)%p1(gg+1),flux(i)%p1(gg+1))
                 F_ind(i)%p1(:)=auxx1(:)
                 flux(i)%p1(:)=auxx2(:)
                 DEALLOCATE(auxx1,auxx2)
              END IF
           END IF

           gg=salidas(destino)   
           IF (gg==0) THEN
              ALLOCATE(F_ind(destino)%p1(1),flux(destino)%p1(1))
              F_ind(destino)%p1(1)=i
              flux(destino)%p1(1)=T_tau(j)
              salidas(destino)=1
           ELSE
              interruptor=.false.
              DO h=1,gg
                 IF (F_ind(destino)%p1(h)==i) THEN
                    flux(destino)%p1(h)=flux(destino)%p1(h)+T_tau(j)
                    interruptor=.true.
                    exit
                 END IF
              END DO
              IF (interruptor.eqv..false.) THEN
                 ALLOCATE(auxx1(gg+1),auxx2(gg+1))
                 auxx1(1:gg)=F_ind(destino)%p1(:)
                 auxx2(1:gg)=flux(destino)%p1(:)
                 auxx1(gg+1)=i
                 auxx2(gg+1)=T_tau(j)
                 salidas(destino)=gg+1
                 DEALLOCATE(F_ind(destino)%p1,flux(destino)%p1)
                 ALLOCATE(F_ind(destino)%p1(gg+1),flux(destino)%p1(gg+1))
                 F_ind(destino)%p1(:)=auxx1(:)
                 flux(destino)%p1(:)=auxx2(:)
                 DEALLOCATE(auxx1,auxx2)
              END IF
           END IF
           
        END IF

     END DO
  END DO

  g=0
  DO i=1,N_nodes
     IF (SL(i)>0) THEN
        salidas(i)=salidas(i)+1
     END IF
  END DO
  g=SUM(salidas(:),DIM=1)



  TT_ind=0
  TT_tau=0
  TT_start=0
  Pe=0

  g=0
  DO i=1,N_nodes
     TT_start(i)=g
     gg=g
     g=g+salidas(i)
     salidas(i)=0
     l=SL(i)
     IF (l>0) THEN
        TT_ind(gg+1)=i
        TT_tau(gg+1)=l
        Pe(i)=l
        salidas(i)=1
     END IF
  END DO
  TT_start(N_nodes+1)=g

  DO i=1,N_nodes
     gg=size(F_ind(i)%p1(:),DIM=1)
     h=0
     DO j=1,gg
        g=salidas(i)+1
        salidas(i)=g
        g=g+TT_start(i)
        l=flux(i)%p1(j)
        TT_ind(g)=F_ind(i)%p1(j)
        TT_tau(g)=l
        h=h+l
     END DO
     Pe(i)=Pe(i)+h
  END DO

  newKmax=MAXVAL(salidas(:),DIM=1)


  DEALLOCATE(F_ind,flux,salidas,SL)





END SUBROUTINE SYMMETRIZE_NET



SUBROUTINE GRAD (N_sets,comunidades,Pe,T_ind,T_tau,T_start,N_nodes,Ktot)


  IMPLICIT NONE
  

  INTEGER,INTENT(IN)::N_nodes,Ktot
  INTEGER,DIMENSION(Ktot),INTENT(IN)::T_ind,T_tau
  INTEGER,DIMENSION(N_nodes+1),INTENT(IN)::T_start
  INTEGER,DIMENSION(N_nodes),INTENT(IN)::Pe

  INTEGER,INTENT(OUT)::N_sets
  INTEGER,DIMENSION(N_nodes),INTENT(OUT)::comunidades


  INTEGER:: i,j,jj,g,h,dim

  !! Para minimos:
  INTEGER::dim2,inicio,fin,peso,candidato
  INTEGER,DIMENSION(:),ALLOCATABLE::tope_lista
  INTEGER,DIMENSION(:,:),ALLOCATABLE::lista
  LOGICAL,DIMENSION(:),ALLOCATABLE::label_glob_min,label_loc_min,filtro
  logical::inter

  !! Para las basins:
  INTEGER::hacia,desde
  logical,dimension(:),allocatable::label
  integer,dimension(:),allocatable::vect_aux1


  hacia=0
  desde=0
  inter=.false.

  dim=1
  dim2=dim

  !! Ordeno primero los contactos

  ALLOCATE(lista(N_nodes,dim2),tope_lista(N_nodes))

  lista=0
  tope_lista=0

  DO i=1,N_nodes

     h=0
     g=T_start(i+1)-T_start(i)

     ALLOCATE(filtro(g))
     filtro=.true.

     inicio=T_start(i)+1
     fin=T_start(i+1)
     
     DO j=inicio,fin
        IF (T_ind(j)==i) THEN
           filtro(j-inicio+1)=.false.
           EXIT
        END IF
     END DO

     DO WHILE (h<dim2)
        
        g=inicio-1+MAXLOC(T_tau(inicio:fin),DIM=1,MASK=filtro)
        candidato=g
        peso=Pe(T_ind(candidato))
        DO j=inicio,fin
           
           IF ((T_tau(j)==T_tau(candidato)).and.(filtro(j-inicio+1).eqv..true.)) THEN
              IF (Pe(T_ind(j))>peso) THEN
                 peso=Pe(T_ind(j))
                 candidato=j
              END IF
           END IF
        END DO
        
        h=h+1
        lista(i,h)=T_ind(candidato)
        filtro(candidato-inicio+1)=.false.
        
        IF (COUNT(filtro)==0) THEN
           exit
        END IF
        
     END DO

     DEALLOCATE(filtro)
     tope_lista(i)=h

  END DO


  ALLOCATE(label_glob_min(N_nodes),label_loc_min(N_nodes))
  label_glob_min=.false.
  label_loc_min=.false.


  !! Minimos globales:

  label_glob_min=.false.

  DO i=1,N_nodes

     peso=Pe(i)
     inter=.true.
     DO j=T_start(i)+1,T_start(i+1)
        g=T_ind(j)
        IF (peso<Pe(g)) THEN
           inter=.false.
           exit
        END IF
     END DO
     
     IF (inter.eqv..true.) label_glob_min(i)=.true.
     
  END DO
  

  !! Minimos locales:

  label_loc_min=.false.
  
  DO i=1,N_nodes
     
     peso=Pe(i)
     
     inter=.true.
     DO j=1,tope_lista(i)
        
        IF (peso<=Pe(lista(i,j))) THEN
           inter=.false.
           exit
        END IF
        
     END DO
     
     IF (inter.eqv..true.) label_loc_min(i)=.true.
     
  END DO


  !! Defino las basins

  ALLOCATE(label(N_nodes))
  label=.false.
  label=label_loc_min
  
  allocate(vect_aux1(N_nodes))
  vect_aux1=0
    
  DO i=1,N_nodes
     IF (label_loc_min(i).eqv..true.) THEN
        comunidades(i)=i
     END IF
  END DO


  DO i=1,N_nodes

     desde=i
     h=0
     
     DO WHILE (label(desde).eqv..false.)

        h=h+1
        vect_aux1(h)=desde
        label(desde)=.true.
        
        DO j=1,tope_lista(desde)
           IF (Pe(lista(desde,j))>=Pe(desde)) THEN
              hacia=lista(desde,j)
              inter=.true.
              exit
           END IF
        END DO

        IF (inter.eqv..false.) THEN
           
           exit
        END IF
        
        desde=hacia
        
     END DO
     

     
     IF (comunidades(desde)==0) THEN
        j=desde
     ELSE
        j=comunidades(desde)
     END IF
     
     DO jj=1,h
        comunidades(vect_aux1(jj))=j
     END DO
     
  END DO


  label=.false.
  DO i=1,N_nodes
     label(comunidades(i))=.true.
  END DO
  h=0

  DO i=1,N_nodes
     IF (label(i).eqv..true.) THEN
        h=h+1
     END IF
  END DO
  

  comunidades=comunidades-1

  N_sets=h



END SUBROUTINE GRAD






END MODULE aux_funcs
