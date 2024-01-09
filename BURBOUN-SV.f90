program BURBOUN       
!
!*** BURBOUN (Block-with-fractURes flow analysis by BOUNdary  element/
!    finite  element  methods): for the simulation of the steady-state
!    2D fluid flow in porous blocks with  discrete fractures
!
!*** version accompanying the submission of the manuscript entitled "A 
!    code for the simulation of the 2D steady-state fluid flow in  po-
!    rous blocks containing transmissive fractures", authors Fidelibus
!    C, Xu C, Wang Z, Dowd P, to Computers & Geosciences, January 2024 
!
!*** MASTER routine
!
!*** GENERAL DATA
!    NBLOCS: number of blocs to be analyzed
!    NELSID: # of boundary elements per side
!    NELEMF: # of boundary elements per fracture
!
!*** BLOCK DATA
!    CONDUC: block hydraulic conductivity
!    NBSIDS: # of sides per block 
!    NCORNR: number of fracture corners
!    NTEFRA: # of elements along fractures
!    NTESID: # of elements along block sides
!    NTSIDS: total # of sides (block sides + fractures)
!    TRANSM: fracture transmissivities
!    XCORNR,YCORNR: block corner coordinates 
!    XFREXT,YFREXT: coordinates of fracture extremities 
!
! rectangular blocks. # of sides NBSIDS equal to 4 (equal to NCORNR); ! no fractures along block edges
!
!*** dimension statements and definitions
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
parameter (MaxNumberBlocks=100)
dimension INDEXV(MaxNumberUnknowns),nfract(MaxNumberBlocks),qtgrou(MaxNumberBlocks,4,4),rhsvec(MaxNumberUnknowns,4),&
sysval(MaxNumberUnknowns,MaxNumberUnknowns),sysva2(MaxNumberUnknowns,MaxNumberUnknowns),&
TEMPSV(MaxNumberUnknowns),transm(MaxNumberFractures),xcornr(MaxNumberCorners+1),&
xfrext(MaxNumberFractures,2),ycornr(MaxNumberCorners+1),yfrext(MaxNumberFractures,2),&
Xbe(2),Ybe(2),dInters(MaxNumberFractures,2),IntersFlag(MaxNumberFractures,2),&
elengt(MaxNumberFractures),ninter(MaxNumberFractures),LabelsIF(MaxNumberFractures,4),&
IntExt(MaxNumberFractures,4)
character*80 TitleR
!
pigree=2.*dasin(1.d0)
!
!*** opening I/O files
!
call arkivi
!
print *,'read general data'
write (32,1001) 
1001  format ('BURBOUN flow analyses of Block with fractURes via BOUNdary and finite elements'/)
read (31,'(a80)') TitleR
write (32,*) TitleR                  
read (31,*) nblocs,nelsid,nelemf
write (32,2001) nblocs,nelsid,nelemf
2001  format (/'# blocks to be analyzed:',i6/&
               '# elements per block side:',i3/&
               '# elements per (inner) fracture:',i3)
!
!*** the code performs an automatic discretization of block boundaries and frac-
!    tures, making use of constant elements  
!
!*** reading input data
!    Xbl,Ybl coordinates bottom-left corner, Xtr,Ytr c. top-right corner
!
if (nelemf.lt.4) then 
    print *,'Warning, the number of fracture elements is less than 4! forced to 4'
    nelemf=4
endif
!
print *,'read data and build equation systems'
!
do 310 jblocs=1,nblocs
LNInters=0
LabelsIF=0
IntExt=0
print *,'block #',jblocs
print *,'read block data'
!
!*** number of inner intersections less/equal than the number of fractures
!
read (31,*) Xbl,Ybl,Xtr,Ytr,conduc,nfract(jblocs),ninter(jblocs)
!
if (Xtr.lt.Xbl.or.Ytr.lt.Ybl) then
        print *,'error: wrong block coordinates!'
        stop
endif
!
write (32,1002) jblocs,Xbl,Ybl,Xtr,Ytr,conduc,nfract(jblocs),ninter(jblocs)
1002  format (/'BLOCK #:',1x,i4/&
              'COORDINATES: Xbl:',1X,g10.4,'Ybl:',1X,g10.4,' Xtr:',1X,g10.4,'Ytr:',1X,g10.4,1X/&
              'CONDUCTIVITY:',1x,g10.4,',# of fractures:',1x,i3,',# of intersections:',1x,i3/)
!
! rectangular blocks, # of corners NCORNR equal to 4
!
ncornr=4
xcornr(1)=Xbl
ycornr(1)=Ybl
xcornr(2)=Xtr
ycornr(2)=Ybl
xcornr(3)=Xtr
ycornr(3)=Ytr
xcornr(4)=Xbl
ycornr(4)=Ytr
!
!*** corner 5 fictitious coincident with corner 1
!
xcornr(5)=Xbl
ycornr(5)=Ybl
!
!*** block sizes
!
Deltax=Xtr-Xbl
Deltay=Ytr-Ybl
dDeltax=Deltax*0.000001
dDeltay=Deltay*0.000001
Xbe(1)=Xbl
Xbe(2)=Xtr
Ybe(1)=Ybl
Ybe(2)=Ytr
dR=sqrt(dDeltax*dDeltax+dDeltay*dDeltay)
!
if (nfract(jblocs).eq.0) then
    write (32,1003) 
1003  format ('NO FRACTURES IN BLOCK')
    go to 111
endif
!
write (32,1004)
1004  format ('EXTREMITIES COORDINATES, FRACTURE TRANSMISSIVITY  '/&
              3x,'FRAC.#',5x,'X1',10x,'Y1',10x,'X2',10x,'Y2',11x,'T')
print *,'read fracture extremities coordinates and other data'
!
do 106 ifract=1,nfract(jblocs)
      read (31,*) xfrext(ifract,1),yfrext(ifract,1),xfrext(ifract,2),yfrext(ifract,2),&
      transm(ifract)
!
do j=1,2
if (xfrext(ifract,j).lt.(Xbe(1)-dDeltax).or.xfrext(ifract,j).gt.(Xbe(2)+dDeltax)) then
        print *,'check X coord. extr. ',j,' fracture ',ifract
        print *,xfrext(ifract,j),'; X<Xmin or X>Xmax'
        stop
endif
if (yfrext(ifract,j).lt.(Ybe(1)-dDeltay).or.yfrext(ifract,j).gt.(Ybe(2)+dDeltay)) then
        print *,'check Y coord. extr. ',j,' fracture ',ifract
        print *,yfrext(ifract,j),'; Y<Ymin or Y>Ymax'
        stop
endif
enddo
!
      fractl=(xfrext(ifract,1)-xfrext(ifract,2))**2+(yfrext(ifract,1)-yfrext(ifract,2))**2
      fractl=sqrt(fractl)
      elengt(ifract)=fractl/nelemf
      write (32,1006) ifract,xfrext(ifract,1),yfrext(ifract,1),xfrext(ifract,2),yfrext(ifract,2),&
transm(ifract)
1006  format (6x,i3,5(1x,g10.4,1x))
106   continue   
!
if (ninter(jblocs).eq.0) then
    write (32,1005) 
1005  format ('NO INNER INTERSECTIONS IN BLOCK')
    go to 111
endif
!
print *,'read fracture incidences per intersection'
write (32,1014)
1014  format ('FRACTURE INCIDENCES PER INTERSECTION'/&
              'INTE. #',4x,'F1',3x,'Ex.',4x,'F2',3x,'Ex.',4x,'F3',3x,'Ex.',4x,&
              'F4',3x,'Ex.')
!
do 490 iinter=1,ninter(jblocs)
read (31,*) (LabelsIF(iinter,i),IntExt(iinter,i),i=1,4)
write (32,1016) iinter,(LabelsIF(iinter,i),IntExt(iinter,i),i=1,4) 
490 continue
!
1016  format (4x,i3,8(1x,i5))
!
111 jfract=nfract(jblocs)
nbsids=ncornr
ntsids=nbsids+jfract
ntesid=nbsids*nelsid
ntefra=nelemf*jfract
nunknw=ntesid+2*ntefra
!
!*** block edges (in counter-clockwise rotation): //2Xb (parallel to X 
!    bottom), //2Yr (parallel to Y right), //2Xt (// to X top), //2Yl
!
!*** boundary conditions: BC1, hydraulic head 1 on //2Xb, 0 on //2Xt,
!    linear from 0 to 1 in //2Yr, //2Yl; BC2, h=1 //2Yr, h=0, //2Yl, 
!    linear from 0 to 1 in //2Xb, //2Xt; BC3, h=1 //2Xt, h=0, //2Xb, 
!    linear from 0 to 1 in //2Yr, //2Yl, BC4, h=1 //2Yl, h=0, //2Yr, 
!    linear from 0 to 1 in //2Xb, //2Xt
!
print *,'build BEM equations (block)'
call blkeqs (jfract,nbsids,nelemf,nelsid,&
             ntefra,ntsids,rhsvec,sysval,xcornr,xfrext,&
             ycornr,yfrext)
if (jfract.eq.0) go to 135
!
!*** fracture equations
!
print *,'build FEM equations (fractures)'
!
!*** loop for all the elements of the fracture IBRANS
!
IntersFlag=0
!
       do 130 ibrans=1,jfract
!
!*** ELENGT fracture element length
!
       zlengt=elengt(ibrans)
       tranlv=transm(ibrans)/conduc
!
!*** Node1Frac: label first element/node of IBRANS fracture (it hosts a discharge) 
!    Node1FracAdd: label element node associated to Node1Frac (a hydr.head)
!
       Node1Frac=ntesid+(ibrans-1)*nelemf+1
       Node1FracAdd=Node1Frac+ntefra      
!
!*** find intersections with block edges
!
do 82 j=1,2
!
do 80 k=1,2
if (xfrext(ibrans,j).ge.(Xbe(k)-dDeltax).and.xfrext(ibrans,j).le.(Xbe(k)+dDeltax)) then
        xfrext(ibrans,j)=Xbe(k)
        IntersFlag(ibrans,j)=4-(k-1)*2
        dInters(ibrans,j)=abs((yfrext(ibrans,j)-Ybe(3-k))/(Ybe(2)-Ybe(1)))
endif
if (yfrext(ibrans,j).ge.(Ybe(k)-dDeltay).and.yfrext(ibrans,j).le.(Ybe(k)+dDeltay)) then
        yfrext(ibrans,j)=Ybe(k)
        IntersFlag(ibrans,j)=1+(k-1)*2
        dInters(ibrans,j)=abs((xfrext(ibrans,j)-Xbe(k))/(Xbe(2)-Xbe(1)))
endif
80 continue
!
82 continue
!
LNInters=ninter(jblocs)
!
206 print *,'analysis fracture #',ibrans
       call matbrs (nelemf,Node1Frac,Node1FracAdd,rhsvec,sysval,tranlv,zlengt,&
                    IntersFlag,dInters,ibrans,LabelsIF,jfract,transm,conduc,elengt,&
                    IntExt,LNInters,ntesid,ntefra) 
130    continue
!
!*** lower-upper decomposition of SYSVAL matrix of coefficients
!
135 print *,'lower-upper decomposition of the matrix of coefficients'
!
do 200 loaded=1,4
sysva2=sysval 
print *,'decompose equation system (sub. routine LUDCMP)'
!
call LUDCMP (sysva2,nunknw,INDEXV,DVALUE)
!
do 190 inodes=1,nunknw
TEMPSV(inodes)=rhsvec(inodes,loaded)
190 continue
!
!*** forward substitution and backsubstitution
!
print *,'solve equation system through substitution (sub. LUBKSB)'
call LUBKSB (sysva2,nunknw,INDEXV,TEMPSV)
!
!*** find the boundary discharges
!
write (32,1116) loaded 
1116  format (/'-----'/'LOADING #',i1/'BLOCK-EDGE DISCHARGES',&
              28x,'BUDGET')
call domino (conduc,jblocs,loaded,nelsid,ntesid,nunknw,qtgrou,TEMPSV,xcornr,ycornr,nfract(jblocs),&
transm,IntersFlag,nelemf,xfrext,yfrext,dInters,ntefra,elengt)
200 continue
!
310  continue
!
print *,'End of run'
end  
!
!*** SUBROUTINE BLKEQS: write BEM equations
!
subroutine blkeqs (jfract,nbsids,nelemf,nelsid,&
                   ntefra,ntsids,rhsvec,sysval,xcornr,xfrext,&
                   ycornr,yfrext)
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
dimension cgausp(4),rhsvec(MaxNumberUnknowns,4),sysval(MaxNumberUnknowns,MaxNumberUnknowns)
dimension xcornr(MaxNumberCorners+1),xfrext(MaxNumberFractures,2)
dimension ycornr(MaxNumberCorners+1),yfrext(MaxNumberFractures,2),wfactr(4)
!
data cgausp /-0.861136311594053,-0.339981043584856,+0.339981043584856,+0.861136311594053/
data wfactr /0.347854845137454,0.652145154862546,0.652145154862546,0.347854845137454/
!
pigree=2.*dasin(1.d0)
ngausp=4
!
!*** loop on the boundary, i.e. along the sides of  the  block and
!    along the fractures inside; NCollN is the collocation node
!
NCollN=0
!
do 199 LatoCl=1,ntsids
!
!*** DXLENG,DYLENG x-length and y-length of the element of the coll. node
!    XFIRST,YFIRST coordinates first node on the segment (side or frac.)
!
if (LatoCl.le.nbsids) then
!
!---segment hosting the collocation node is a block side
!
xfirst=xcornr(LatoCl)
xsecon=xcornr(LatoCl+1)
yfirst=ycornr(LatoCl)
ysecon=ycornr(LatoCl+1)
dxleng=xsecon-xfirst
dyleng=ysecon-yfirst
dxleng=dxleng/nelsid
dyleng=dyleng/nelsid
NumElm=nelsid
!
else
!
!---segment hosting the collocation node is a fracture
!
xfirst=xfrext(LatoCl-nbsids,1)
xsecon=xfrext(LatoCl-nbsids,2)
yfirst=yfrext(LatoCl-nbsids,1)
ysecon=yfrext(LatoCl-nbsids,2)
dxleng=xsecon-xfirst
dyleng=ysecon-yfirst
dxleng=dxleng/nelemf
dyleng=dyleng/nelemf
NumElm=nelemf
endif
!
!*** XCOLLC,YCOLLC coordinates of the collocation node 
!
do 189 LCollN=1,NumElm
NCollN=NCollN+1
xcollc=xfirst-dxleng/2+LCollN*dxleng
ycollc=yfirst-dyleng/2+LCollN*dyleng
!
!*** calculation of the free term for the current collocation node
!
freetr=0.5
if (LatoCl.le.nbsids) then
!
    if (LatoCl.eq.1) then
    rhsvec(NCollN,1)=freetr*1.
    rhsvec(NCollN,2)=freetr*(xcollc-xfirst)/(xsecon-xfirst)
    rhsvec(NCollN,3)=freetr*0.
    rhsvec(NCollN,4)=freetr*(xsecon-xcollc)/(xsecon-xfirst)
    endif
!
    if (LatoCl.eq.2) then
    rhsvec(NCollN,1)=freetr*(ysecon-ycollc)/(ysecon-yfirst)
    rhsvec(NCollN,2)=freetr*1.
    rhsvec(NCollN,3)=freetr*(ycollc-yfirst)/(ysecon-yfirst)
    rhsvec(NCollN,4)=freetr*0.
    endif
!
    if (LatoCl.eq.3) then
    rhsvec(NCollN,1)=freetr*0.
    rhsvec(NCollN,2)=freetr*(xcollc-xsecon)/(xfirst-xsecon)
    rhsvec(NCollN,3)=freetr*1.
    rhsvec(NCollN,4)=freetr*(xfirst-xcollc)/(xfirst-xsecon)
    endif
!
    if (LatoCl.eq.4) then
    rhsvec(NCollN,1)=freetr*(yfirst-ycollc)/(yfirst-ysecon)
    rhsvec(NCollN,2)=freetr*0. 
    rhsvec(NCollN,3)=freetr*(ycollc-ysecon)/(yfirst-ysecon)
    rhsvec(NCollN,4)=freetr*1.
    endif
!
else
! 
    LabelUnknown=NCollN+jfract*nelemf
    sysval(NcollN,LabelUnknown)=-freetr*2.
endif
!
!*** NINTEG integration node
!
ninteg=0
!
do 127 jsides=1,ntsids
!
!*** NELINT number of elements to be subjected to the integration  for
!    the current JSIDES integration side
!
if (jsides.le.nbsids) go to 90
nelint=nelemf
xlatp1=xfrext(jsides-nbsids,1)
xlatp2=xfrext(jsides-nbsids,2)
ylatp1=yfrext(jsides-nbsids,1)
ylatp2=yfrext(jsides-nbsids,2)
go to 91
90  nelint=nelsid
xlatp1=xcornr(jsides)
xlatp2=xcornr(jsides+1)
ylatp1=ycornr(jsides)
ylatp2=ycornr(jsides+1)
91 deltx2=(xlatp2-xlatp1)/nelint
delty2=(ylatp2-ylatp1)/nelint
!
!*** calculation of the distance between the side and the current col-
!    location node
!
92 slengt=sqrt((xlatp2-xlatp1)**2+(ylatp2-ylatp1)**2)
!
if (slengt.le.1.e-07) then
        write (*,*) 'ERROR! TOO SMALL BLOCK EDGE OR FRACTURE ',jsides 
        STOP
endif
! 
elcalc=slengt/nelint
xnorm1=(ylatp2-ylatp1)/slengt
ynorm1=-(xlatp2-xlatp1)/slengt
!
if (jsides.eq.LatoCl.or.jsides.gt.nbsids) then
    distan=0.
else
    distan=abs((xlatp1-xcollc)*xnorm1+(ylatp1-ycollc)*ynorm1)
endif
!
elengt=elcalc
!
!*** loop for all the elements on the  side
!
do 125 jelems=1,nelint
ninteg=ninteg+1
!
if (ninteg.eq.NCollN) then
!
!*** quick  logarithmic  integration for constant element                
!
    sysval(NCollN,NCollN)=sysval(NCollN,NCollN)+elcalc*(1+alog(1/(0.5*elcalc)))/2./pigree
    go to 125
endif
!
xcentr=xlatp1+deltx2*0.5+deltx2*(jelems-1)
ycentr=ylatp1+delty2*0.5+delty2*(jelems-1)
!
quantH=0.
quantG=0.
!
101   do 110 jgausp=1,ngausp
csival=cgausp(jgausp)
xgausp=deltx2*0.5*csival+xcentr    
ygausp=delty2*0.5*csival+ycentr   
radius=sqrt((xgausp-xcollc)**2+(ygausp-ycollc)**2)
if (radius.lt.1.e-09) STOP
!
quantH=quantH-distan/radius/radius*elcalc*0.5*wfactr(jgausp)/(2*pigree)
quantG=quantG+alog(1/radius)*elcalc*0.5*wfactr(jgausp)/(2*pigree) 
110 continue
!
if (jsides.le.nbsids) then
sysval(NCollN,ninteg)=sysval(NCollN,ninteg)+quantG
!
    if (jsides.eq.1) then
    rhsvec(NCollN,1)=rhsvec(NCollN,1)+quantH*1.
    rhsvec(NCollN,2)=rhsvec(NCollN,2)+quantH*(xcentr-xlatp1)/(xlatp2-xlatp1)
    rhsvec(NCollN,3)=rhsvec(NCollN,3)+quantH*0.
    rhsvec(NCollN,4)=rhsvec(NCollN,4)+quantH*(xlatp2-xcentr)/(xlatp2-xlatp1)
    endif
!
    if (jsides.eq.2) then
    rhsvec(NCollN,1)=rhsvec(NCollN,1)+quantH*(ylatp2-ycentr)/(ylatp2-ylatp1)
    rhsvec(NCollN,2)=rhsvec(NCollN,2)+quantH*1.
    rhsvec(NCollN,3)=rhsvec(NCollN,3)+quantH*(ycentr-ylatp1)/(ylatp2-ylatp1)
    rhsvec(NCollN,4)=rhsvec(NCollN,4)+quantH*0.
    endif
!
    if (jsides.eq.3) then
    rhsvec(NCollN,1)=rhsvec(NCollN,1)+quantH*0.
    rhsvec(NCollN,2)=rhsvec(NCollN,2)+quantH*(xcentr-xlatp2)/(xlatp1-xlatp2)
    rhsvec(NCollN,3)=rhsvec(NCollN,3)+quantH*1.
    rhsvec(NCollN,4)=rhsvec(NCollN,4)+quantH*(xlatp1-xcentr)/(xlatp1-xlatp2)
    endif
!
    if (jsides.eq.4) then
    rhsvec(NCollN,1)=rhsvec(NCollN,1)+quantH*(ylatp1-ycentr)/(ylatp1-ylatp2)
    rhsvec(NCollN,2)=rhsvec(NCollN,2)+quantH*0.
    rhsvec(NCollN,3)=rhsvec(NCollN,3)+quantH*(ycentr-ylatp2)/(ylatp1-ylatp2)
    rhsvec(NCollN,4)=rhsvec(NCollN,4)+quantH*1.
    endif
!
else
sysval(NCollN,ninteg)=sysval(NCollN,ninteg)+quantG
sysval(NCollN,ninteg+ntefra)=sysval(NCollN,ninteg+ntefra)+quantH
endif
!
125   continue
!
127   continue
!
189   continue
!
199   continue
!
return
end
!
!*** SUBROUTINE MATBRS: write FEM equations for the fractures (traces)
!
subroutine matbrs (nelemf,Node1Frac,Node1FracAdd,rhsvec,sysval,tranlv,zlengt,&
                   IntersFlag,dInters,ibrans,LabelsIF,NFracL,transm,conduc,elengt,&
                   IntExt,LNInters,ntesid,ntefra)
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
parameter (MaxNumberBlocks=100)
dimension rhsvec(MaxNumberUnknowns,4),sysval(MaxNumberUnknowns,MaxNumberUnknowns)
dimension IntersFlag(MaxNumberFractures,2),dInters(MaxNumberFractures,2),LabelsIF(MaxNumberFractures,4),&
IntExt(MaxNumberFractures,4),elengt(MaxNumberFractures),transm(MaxNumberFractures),&
LFAtIn(4),LEAtIn(4)
!
!*** rshvec=0 (except for the first or last node when a fracture extremity lies on a block edge)
!
      do 10 jnelem=1,nelemf
      IPivotFEMEq=Node1FracAdd+jnelem-1
      rhsvec(IPivotFEMEq,1)=0.
      rhsvec(IPivotFEMEq,2)=0.
      rhsvec(IPivotFEMEq,3)=0.
      rhsvec(IPivotFEMEq,4)=0.
10    continue
!
!*** first element/node of the fracture; Node1Frac is the first element/node (hosts Q)
!    Node1FracAdd is the unknown corresponding to the first hydraulic head
!
      IPivotFEMEq=Node1FracAdd
!
!*** check if the extremity 1 coincides with an intersection node
!
      IGotIt=0
      if (LNInters.eq.0) go to 28 
      jI=1
!
15    if (jI.eq.(LNInters+1).and.IGotIt.eq.0) goto 28 
!
      do jE=1,4
      if (LabelsIF(jI,jE).eq.ibrans.and.IntExt(jI,jE).eq.1) then
!
          do jF=1,4
          LFAtIn(jF)=LabelsIF(jI,jF)
          LEAtIn(jF)=IntExt(jI,jF)
          IGotIt=1
          enddo
!
      endif
!
      enddo
!
      if (IGotIt.eq.0) then 
          jI=jI+1 
          goto 15
      endif 
!
!*** 1st extremity coincides with a inner intersection node (i.e. inters. among fractures)
!
    SumNtDivLe=0.
!
    do 18 kOIntF=1,4
    if (LEAtIn(kOIntF).eq.0) go to 18 
    dlengt=elengt(LFAtIN(kOIntF))
    SumNtDivLe=SumNtDivLe+2.0*transm(LFAtIN(kOIntF))/conduc/dlengt
18  continue
!
!*** coefficients equation for the node when coincides with an inner inters.
!
    sysval(IPivotFEMEq,IPivotFEMEq)=-3.*tranlv/zlengt
    sysval(IPivotFEMEq,IPivotFEMEq+1)=1.*tranlv/zlengt
    sysval(IPivotFEMEq,Node1Frac)=-5.*zlengt/8.
    sysval(IPivotFEMEq,Node1Frac+1)=-1.*zlengt/8.
!
    do 26 kOIntF=1,4
    if (LFAtIn(kOIntF).eq.0) go to 26
     CoeffS=1.0
!
    if (LEAtIn(kOIntF).eq.1) then
    INIntersF=ntesid+(LFAtIn(kOIntF)-1)*nelemf+1+ntefra
    INIntersq=ntesid+(LFAtIn(kOIntF)-1)*nelemf+1
    else
    INIntersF=ntesid+LFAtIn(kOIntF)*nelemf+ntefra
    INIntersq=ntesid+LFAtIn(kOIntF)*nelemf
    endif
!
    dlengt=elengt(LFAtIn(kOIntF))
    sysval(IPivotFEMEq,INIntersF)=sysval(IPivotFEMEq,INIntersF)+&
CoeffS*2.0*(tranlv/zlengt)*2.0*transm(LFAtIn(kOIntF))/&
           conduc/dlengt/SumNtDivLe 
    sysval(IPivotFEMEq,INIntersq)=sysval(IPivotFEMEq,INIntersq)-&
CoeffS*2.0*(tranlv/zlengt)*0.25*dlengt/SumNtDivLe 
26  continue
!
    go to 35
!
!*** check if the extremity 1 intersects a block edge              
!
28    if (IntersFlag(ibrans,1).eq.0) go to 30
!
!*** fracture extremity on a block edge or coincides with an intersection
!
      sysval(IPivotFEMEq,IPivotFEMEq)=-3.*tranlv/zlengt
      sysval(IPivotFEMEq,IPivotFEMEq+1)=1.*tranlv/zlengt
      sysval(IPivotFEMEq,Node1Frac)=-5.*zlengt/8.
      sysval(IPivotFEMEq,Node1Frac+1)=-1.*zlengt/8.
      if (IGotIt.ne.0) go to 35 
      c=-2.*tranlv/zlengt
      Lato=IntersFlag(ibrans,1)
      rhsvec(IPivotFEMEq,Lato)=1.*c
      k=Lato+1-4*INT((Lato)/4)
      rhsvec(IPivotFEMEq,k)=dInters(ibrans,1)*c
      k=Lato+2-4*INT((Lato+1)/4)
      rhsvec(IPivotFEMEq,k)=0.
      k=Lato+3-4*INT((Lato+2)/4)
      rhsvec(IPivotFEMEq,k)=(1.-dInters(ibrans,1))*c
      go to 35
!
!*** inner extremity
!
30    sysval(IPivotFEMEq,IPivotFEMEq)=-1.*tranlv/zlengt
      sysval(IPivotFEMEq,IPivotFEMEq+1)=1.*tranlv/zlengt
      sysval(IPivotFEMEq,Node1Frac)=-7.*zlengt/8.
      sysval(IPivotFEMEq,Node1Frac+1)=-1.*zlengt/8.
!
!*** all the nodes but the first and the last
!
35    do 40 inelem=2,nelemf-1
      IPivotFEMEq=Node1FracAdd-1+inelem
      labHnd=Node1FracAdd-1+inelem
      labQnd=Node1Frac-1+inelem
      sysval(IPivotFEMEq,labHnd)=-2.*tranlv/zlengt
      sysval(IPivotFEMEq,labHnd+1)=1.*tranlv/zlengt
      sysval(IPivotFEMEq,labHnd-1)=1.*tranlv/zlengt
      sysval(IPivotFEMEq,labQnd)=-3.*zlengt/4.
      sysval(IPivotFEMEq,labQnd+1)=-1.*zlengt/8.
      sysval(IPivotFEMEq,labQnd-1)=-1.*zlengt/8.
40    continue
!
!*** last element/node of the fracture
!
      IPivotFEMEq=Node1FracAdd+nelemf-1
!
!*** check if the extremity 2 coincides with an intersection node
!
      IGotIt=0
      if (LNInters.eq.0) go to 48 
      jI=1
!
45    if (jI.eq.(LNInters+1).and.IGotIt.eq.0) goto 48 
!
      do jE=1,4
      if (LabelsIF(jI,jE).eq.ibrans.and.IntExt(jI,jE).eq.2) then
!
          do jF=1,4
          LFAtIn(jF)=LabelsIF(jI,jF)
          LEAtIn(jF)=IntExt(jI,jF)
          IGotIt=1
          enddo
!
      endif
!
      enddo
!
      if (IGotIt.eq.0) then 
          jI=jI+1 
          goto 45
      endif 
!
!*** last extremity coincides with a inner intersection node (i.e. inters. among fractures)
!
    SumNtDivLe=0.
!
    do 46 kOIntF=1,4
    if (LEAtIn(kOIntF).eq.0) go to 46 
    dlengt=elengt(LFAtIN(kOIntF))
    SumNtDivLe=SumNtDivLe+2.0*transm(LFAtIN(kOIntF))/conduc/dlengt
46  continue
!
!*** coefficients equation for the node when coincides with an inner inters.
!
    sysval(IPivotFEMEq,IPivotFEMEq)=-3.*tranlv/zlengt
    sysval(IPivotFEMEq,IPivotFEMEq-1)=1.*tranlv/zlengt
    sysval(IPivotFEMEq,Node1Frac+nelemf-1)=-5.*zlengt/8. 
    sysval(IPivotFEMEq,Node1Frac+nelemf-2)=-1.*zlengt/8. 
!
    do 47 kOIntF=1,4
    if (LFAtIn(kOIntF).eq.0) go to 47
    CoeffS=(-1)**(LEAtIn(kOIntF))
!
    if (LEAtIn(kOIntF).eq.1) then
    INIntersF=ntesid+(LFAtIn(kOIntF)-1)*nelemf+1+ntefra
    INIntersq=ntesid+(LFAtIn(kOIntF)-1)*nelemf+1
    else
    INIntersF=ntesid+LFAtIn(kOIntF)*nelemf+ntefra
    INIntersq=ntesid+LFAtIn(kOIntF)*nelemf
    endif
!
    dlengt=elengt(LFAtIn(kOIntF))
    sysval(IPivotFEMEq,INIntersF)=sysval(IPivotFEMEq,INIntersF)+&
    2.0*(tranlv/zlengt)*2.0*transm(LFAtIn(kOIntF))/&
           conduc/dlengt/SumNtDivLe 
    sysval(IPivotFEMEq,INIntersq)=sysval(IPivotFEMEq,INIntersq)-&
       2.0*(tranlv/zlengt)*0.25*dlengt/SumNtDivLe 
47  continue
!
    go to 60
!
48    if (IntersFlag(ibrans,2).gt.0) go to 50
      sysval(IPivotFEMEq,IPivotFEMEq-1)=tranlv/zlengt
      sysval(IPivotFEMEq,IPivotFEMEq)=-1.*tranlv/zlengt
      sysval(IPivotFEMEq,Node1Frac+nelemf-2)=-zlengt/8.
      sysval(IPivotFEMEq,Node1Frac+nelemf-1)=-7.*zlengt/8.
      go to 60
!
50    sysval(IPivotFEMEq,IPivotFEMEq)=-3.*tranlv/zlengt
      sysval(IPivotFEMEq,IPivotFEMEq-1)=1.*tranlv/zlengt
      sysval(IPivotFEMEq,Node1Frac+nelemf-1)=-5.*zlengt/8.
      sysval(IPivotFEMEq,Node1Frac+nelemf-2)=-1.*zlengt/8.
      c=-2.*tranlv/zlengt
      Lato=IntersFlag(ibrans,2)       !correzione 26/10/22
      rhsvec(IPivotFEMEq,Lato)=1.*c
      k=Lato+1-4*INT((Lato)/4)
      rhsvec(IPivotFEMEq,k)=dInters(ibrans,2)*c
      k=Lato+2-4*INT((Lato+1)/4)
      rhsvec(IPivotFEMEq,k)=0.
      k=Lato+3-4*INT((Lato+2)/4)
      rhsvec(IPivotFEMEq,k)=(1.-dInters(ibrans,2))*c
60    return
      end
!
!
subroutine domino (conduc,jblocs,loaded,nelsid,ntesid,nunknw,qtgrou,TEMPSV,xcornr,ycornr,LNumberFractures,&
transm,IntersFlag,nelemf,xfrext,yfrext,dInters,ntefra,elengt)
!
!*** this subroutine gathers all the solutions for output with reference to block JBLOCS and load cond. LOADED
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
parameter (MaxNumberBlocks=100)
dimension qtgrou(MaxNumberBlocks,4,4),TEMPSV(MaxNumberUnknowns),xcornr(MaxNumberCorners+1),&
ycornr(MaxNumberCorners+1),transm(MaxNumberFractures),IntersFlag(MaxNumberFractures,2),&
xfrext(MaxNumberFractures,2),yfrext(MaxNumberFractures,2),dInters(MaxNumberFractures,2),&
elengt(MaxNumberFractures),cmnorm(4,2)
!
!*** CMNORM components of the vectors normal to the block sides
!
cmnorm(1,1)=0.
cmnorm(1,2)=-1.
cmnorm(2,1)=1.
cmnorm(2,2)=0.
cmnorm(3,1)=0.
cmnorm(3,2)=1.
cmnorm(4,1)=-1.
cmnorm(4,2)=0.
!
knodes=0
!
do 60 jsides=1,4
!
!*** loop for each JSIDES side of the block
!
qgloba=0.
xnode1=xcornr(jsides)
ynode1=ycornr(jsides)
xcnode=xnode1
ycnode=ynode1
xlastn=xcornr(jsides+1)
ylastn=ycornr(jsides+1)
deltax=(xlastn-xnode1)/nelsid
deltay=(ylastn-ynode1)/nelsid
zlengt=sqrt(deltax**2+deltay**2)
!
do 40 jelems=1,nelsid
knodes=knodes+1
qgloba=qgloba+zlengt*TEMPSV(knodes)
40    continue
!
qtgrou(jblocs,loaded,jsides)=-qgloba*conduc
!
!*** block of instructions for collecting fracture discharges at the edges
!
QFractures=0.
QSingleF=0.
!
!*** find intersections among JSIDES edge and block fractures
!
do 45 kfract=1,LNumberFractures
if (IntersFlag(kfract,1).ne.jsides.and.IntersFlag(kfract,2).ne.jsides) go to 45
!
do 43 iextre=1,2
if (IntersFlag(kfract,iextre).ne.jsides) go to 43
!
       if (loaded.eq.jsides) VHEdge=1.
     jplus1=jsides+1
     if (jplus1.gt.4) jplus1=jplus1-4
       if (loaded.eq.jplus1) VHEdge=dInters(kfract,iextre)
     jplus2=jsides+2
     if (jplus2.gt.4) jplus2=jplus2-4
       if (loaded.eq.jplus2) VHEdge=0.
     jplus3=jsides+3
     if (jplus3.gt.4) jplus3=jplus3-4
       if (loaded.eq.jplus3) VHEdge=1.-dInters(kfract,iextre)
!
     if (iextre.eq.1) then
     VHExtremeNode=TEMPSV(ntesid+ntefra+(kfract-1)*nelemf+1)
     VqExtremeNode=TEMPSV(ntesid+(kfract-1)*nelemf+1)
     DeltaH=VHExtremeNode-VHEdge
     else
     VHExtremeNode=TEMPSV(ntesid+ntefra+kfract*nelemf)
     VqExtremeNode=TEMPSV(ntesid+kfract*nelemf)
     DeltaH=-1.*(VHExtremeNode-VHEdge)
     endif
43 continue
!
ScalarProduct=cmnorm(jsides,1)*(xfrext(kfract,2)-xfrext(kfract,1))+&
              cmnorm(jsides,2)*(yfrext(kfract,2)-yfrext(kfract,1))
QSingleF=-2.0*transm(kfract)*DeltaH/elengt(kfract)/conduc+0.25*VqExtremeNode*elengt(kfract)
QSingleF=QSingleF*conduc
if (ScalarProduct.lt.0) QSingleF=-QSingleF
QFractures=QFractures+QSingleF
!
!*** end loop for fractures intersecting edge JSIDES
!
45 continue
!
qtgrou(jblocs,loaded,jsides)=qtgrou(jblocs,loaded,jsides)+QFractures
!
!*** end loop for block sides
!
60    continue
!
budget=0.
!
do l=1,4
budget=budget+qtgrou(jblocs,loaded,l)
enddo
!
write (32,1117) (qtgrou(jblocs,loaded,jsides),jsides=1,4),budget
!
1117  format (5(1x,g10.4)/)
!
write (32,8001)
8001  format('NODAL VALUES')
write (32,8002)
8002  format('FLUXES AT THE BLOCK EDGES')
write (32,8003)
8003  format('BLOCK 1')
write (32,8004) (iunknw,TEMPSV(iunknw)*conduc,iunknw=1,nelsid)
8004  format(5(1x,i4,1x,g10.4))
write (32,8005)
8005  format('BLOCK 2')
write (32,8004) (iunknw,TEMPSV(iunknw)*conduc,iunknw=nelsid+1,2*nelsid)
write (32,8006)
8006  format('BLOCK 3')
write (32,8004) (iunknw,TEMPSV(iunknw)*conduc,iunknw=2*nelsid+1,3*nelsid)
write (32,8007)
8007  format('BLOCK 4')
write (32,8004) (iunknw,TEMPSV(iunknw)*conduc,iunknw=3*nelsid+1,4*nelsid)
!
write (32,8008)
8008  format('FLUXES NORMAL TO THE FRACTURES')
!
do ilfrac=1,LNumberFractures
write (32,8009) ilfrac
8009  format('FRACTURE #',i3)
write (32,8004) (iunknw,TEMPSV(iunknw)*conduc,iunknw=ntesid+(ilfrac-1)&
                 *nelemf+1,ntesid+ilfrac*nelemf)
enddo
!
write (32,8010)
8010  format('HYDRAULIC HEADS AT THE FRACTURES')
!
do ilfrac=1,LNumberFractures
write (32,8009) ilfrac
write (32,8004) (iunknw,TEMPSV(iunknw),iunknw=ntesid+ntefra+&
                (ilfrac-1)*nelemf+1,ntesid+ntefra+ilfrac*nelemf)
enddo
!
1115  format ('BLOCK',1x,i3,', SUB-TRACE',1x,i3,', DISCHARGE:',1x,g10.4,', MEAN HEAD:',1x,g10.4)
!
return
end
!
subroutine arkivi
      CHARACTER EXTEN1*4,EXTEN2*4,FORMTI*5,IFNAM2*34,IFNAME*30,&
                OFNAME*34,SINGLC*1,TEXFIL*36 
      LOGICAL*4 FILEXI                                                 
!                                                                     
      MAXLEN=30                                                     
      FORMTI='(A30)'                                                
      EXTEN1='.DAT'
      EXTEN2='.RES'
      SINGLC=' '
!                                                                     
!*** read input file name and open it
!                                                                     
100   WRITE (*,2100) MAXLEN                                           
      READ  (*,FORMTI) IFNAME                                           
!
      IF (IFNAME.EQ.' ')  GO TO 100
!
      LENNIF=0
      LENNIF=INDEX(IFNAME,SINGLC)
      LENNIF=LENNIF-1
      IF (LENNIF.LE.0) LENNIF=MAXLEN
      IFNAM2=IFNAME(1:LENNIF)//EXTEN1
      WRITE (*,2120) IFNAM2                                           
      INQUIRE (FILE=IFNAM2, EXIST=FILEXI)                              
!
      IF (FILEXI) THEN                                                  
         OPEN (31,FILE=IFNAM2,STATUS='OLD')                     
      ELSE                                                            
         WRITE (*,2160)                                               
         GO TO 100                                                    
      ENDIF                                                           
!
      OFNAME=IFNAME(1:LENNIF)//EXTEN2
      OPEN (32,FILE=OFNAME)                     
!
!*** format statements 
!                                                                     
2100  format (/'input filename (without extension)? default = INPUT',&      
                     ' - ',i2,' max. char.' )                   
2120  format (/'input filename: ', a30 )                 
2160  format (/'error _ filename does not exist!')     
!
      RETURN
      END
!
subroutine LUDCMP(sysva2,nunknw,INDEXV,DVALUE)
!
!*** given the matrix SYSVA2, this subroutine replaces it by  the lo-
!    wer-upper decomposition of a rowwise  permutation of itself.  At
!    the end SYSVAL is modified, ready for the  forward  substitution
!    and backsubstitution to be performed by LUBKSB. INDEXV is a out-
!    put vector containing the row permutations effected by  the par-
!    tial pivoting; DVALUE is a flag equal to plus-minus 1 indicating
!    even or odd number of row interchanges respectively.
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
dimension sysva2(MaxNumberUnknowns,MaxNumberUnknowns),INDEXV(MaxNumberUnknowns),vv(MaxNumberUnknowns)
! 
!*** VV stores the implicit scaling of each row
!
tiny=1.E-20
DVALUE=1.
!
!*** loop over the rows to get the implicit scaling information
!
do 12 i=1,nunknw
aamax=0.
!
do 11 j=1,nunknw
if (abs(sysva2(i,j)).gt.aamax) aamax=abs(sysva2(i,j))
11    continue
!
if (aamax.eq.0.) then
WRITE (*,*) '!ERROR: singular matrix'
READ (*,'()')
endif
vv(i)=1./aamax
!
!*** saving the scaling information
!
12    continue
!
!*** loop over the columns as requested by Crout's method
!
do 19 j=1,nunknw
!
do 14 i=1,j-1
sum=sysva2(i,j)
!
do 13 k=1,i-1
sum=sum-sysva2(i,k)*sysva2(k,j)
13    continue
!
sysva2(i,j)=sum
14    continue
!
!*** initializing for the search of the largest pivot
!
aamax=0.
!*** initializing for the search of the largest pivot
!
aamax=0.
!
do 16 i=j,nunknw
sum=sysva2(i,j)
!
do 15 k=1,j-1
sum=sum-sysva2(i,k)*sysva2(k,j)
15    continue
!
sysva2(i,j)=sum
dum=vv(i)*abs(sum)
!
if (dum.ge.aamax) then
imax=i
aamax=dum
endif
!
16    continue
!
!*** interchanging rows, if it is needed
!
if (j.ne.imax) then
!
do 17 k=1,nunknw
dum=sysva2(imax,k)
sysva2(imax,k)=sysva2(j,k)
sysva2(j,k)=dum
17        continue
!
DVALUE=-DVALUE
vv(imax)=vv(j)
endif
!
INDEXV(j)=imax
if (sysva2(j,j).eq.0) sysva2(j,j)=tiny
!
!*** dividing by the pivot
!
if (j.ne.nunknw) then
dum=1./sysva2(j,j)
!
do 18 i=j+1,nunknw
sysva2(i,j)=sysva2(i,j)*dum
18        continue
!
endif
!
19    continue
!
return
end
!
!
subroutine LUBKSB(sysva2,nunknw,INDEXV,TEMPSV)
!
!*** the A matrix is the lower-upper decomposition of the  old  SYSVA2
!    as output of LUDCMP. INDEXV is the permutation vector as produced
!    by LUDCMP. TEMPSV is the right hand side vector. At  the  end  it
!    contains the solutions
!
parameter (MaxNumberCorners=4,MaxNumberFractures=70,MaxNElementsXSide=10)
parameter (MaxNumberUnknowns=MaxNElementsXSide*MaxNumberCorners+2*MaxNElementsXSide*MaxNumberFractures)
dimension sysva2(MaxNumberUnknowns,MaxNumberUnknowns),INDEXV(MaxNumberUnknowns),TEMPSV(MaxNumberUnknowns)
!
ii=0
!
do 12 i=1,nunknw
LL=INDEXV(i)
sum=TEMPSV(LL)
TEMPSV(LL)=TEMPSV(i)
!
if (ii.ne.0) then
!
do 11 j=ii,i-1
sum=sum-sysva2(i,j)*TEMPSV(j)
11        continue
!
else if (sum.ne.0.) then
ii=i
endif
!
TEMPSV(i)=sum
12    continue
!
do 14 i=nunknw,1,-1
sum=TEMPSV(i)
!
if (i.lt.nunknw) then
!
do 13 j=i+1,nunknw
sum=sum-sysva2(i,j)*TEMPSV(j)
13        continue
!
endif
!
TEMPSV(i)=sum/sysva2(i,i)
14    continue
!
return
end
