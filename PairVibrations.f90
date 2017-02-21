       module Stuff
       implicit none
!----- Definition of vector sizes -----!
       integer     , parameter :: lmax=15,mphon=250,mfrag=10,NRMAX=100000
       integer     , parameter :: Nbox=150
!----- Definition of mesh and limit of dispersion relation -----!
       real(kind=8), parameter :: dr=0.1, dW = 0.01d0, EcutFunction = 30.d0
!---------------------------------------!

       real(kind=8) ::  ek,w,betal,BE_Lambda,hcoup,emin,emax,              &
     &                  Energy,eqp1,eqpcorr1,zeta,zeta1,efn,Uini,Vini,     &
     &                  PairAddition,PairRemoval,V,Add,Rem,park,           &
     &                  NormAdd,NormRem,SumAdd,SumRem,Wintegrale,       &
     &                  DiffAddition,DiffRemoval,GAdd,GRem,add1,rem1,      &
     &                  Wintegrale_add,Wintegrale_rem,cadd,crem,jfac,e_sp

       real(kind=8) :: diff_add, diff_rem, dispadd, disprem, ek1, emin_add,&
     &                 emin_rem, xadd, xrem, yadd, yrem, Wminimo


       integer      ::  llk, jjk,Nlivelli,ichoose,nfon,lambda2,       &
     &                  ilph,n,frag,jjk1,llk1,frag1,              &
     &                  igapUp,igapDn, ir_tot

       dimension    ::  ek(mphon),e_sp(mphon),llk(mphon),jjk(mphon),w(lmax,mphon),     &
     &                  betal(lmax,mphon),BE_Lambda(lmax,mphon),           &
     &                  nfon(lmax),frag(mphon),zeta(mphon),V(mphon),       &
     &                  add(mphon),Rem(mphon),                             &
     &                  DispAdd(NRmax),DispRem(NRmax),                     &
     &                  Xadd(lmax),Yadd(lmax),Xrem(lmax),Yrem(lmax)

       real(kind=8), dimension(:), allocatable :: DispersionFunction

       character (len=1) :: line

       contains

!---------------------------!
       subroutine read_define(icalc, flag_efn)

       implicit none

       integer :: i,ii, icalc

       logical :: flag_efn

       real(kind=4) :: park1, park4, park5, park6, park7, park8


       if(icalc.eq.0)then
        open(unit=01,file='Spectroscopic_factors.dat')

        n=0
        do i=1,mphon
 
         read(01,*,end=100)ii,ichoose,llk1,jjk1,frag1,ek1,eqp1,eqpcorr1,   &
     &             park1,Vini,Uini,park4,park5,park6,park7,park8,zeta1

         if(frag1.eq.1)then
          n=n+1                        !Conto il numero dei livelli
          llk(n)=llk1
          jjk(n)=jjk1
           V(n)=Vini
          if(Vini.gt.0.1)ek(n)=abs((ek1+eqp1-eqpcorr1)-efn)
          if(Vini.le.0.1)ek(n)=abs((ek1-eqp1+eqpcorr1)-efn)
          Zeta(n)=Zeta1
         endif
        enddo

100     continue

        NLivelli=N

       elseif(icalc.eq.1)then
        open(unit=01,file='Livelli_Sperimentali.dat')
        read(01,*,end=101)line
        if(line(1:1) /= "#")then
           write(*,*)'Old input! Put Addition and Removal Correlation Energies on top'
           stop
        endif
        read(01,*)DiffAddition, DiffRemoval
        read(01,*)
        Zeta=1.d0
        Nlivelli=0
        do i=1,mphon
         read(01,*,end=101)llk(i),jjk(i),e_sp(i),V(i),Zeta(i)

          if(V(i).gt.0.5)then
             iGapDn=i;iGapUp=i+1
          endif

!         ek(i)= (ek(i)) - efn
         NLivelli=NLivelli+1
        enddo
101     continue

        write(*,*)'Last  Occupied Level ek=',e_sp(iGapDn)
        write(*,*)'First Empty    Level ek=',e_sp(iGapUp)
        if(.not.flag_efn)efn =  (e_sp(iGapDn)+e_sp(iGapUp))/2.d0
        write(*,*)'Fermi Energy           =',efn
        do i=1,NLivelli
         ek(i)=abs(e_sp(i)-efn)
        enddo
 
       else
        write(*,*)'errore:0 oppure 1!'
        stop
       endif
       end subroutine read_define

!----------------------------------------------!
       subroutine disp_relation

       implicit none
       integer :: i,ir

       DispersionFunction = 0.d0

       write(*,*)'Calcolare Dispersioni Pair Addition and Removal Mode'

       PairAddition=2*ek(iGapUp)-DiffAddition
       PairRemoval=2*ek(iGapDn)-DiffRemoval
       write(*,*)DiffAddition,PairAddition,DiffRemoval,PairRemoval
       ir = 0; WIntegrale = 0.0025d0

       emin_add = 999.; emin_rem = 999.

       do while (Wintegrale.lt.EcutFunction)
        SumAdd=0.
        SumRem=0.
        DispAdd(ir) =0.
        DispRem(ir)= 0.


        do i=1,NLivelli
          jfac= 0.5d0*(jjk(i)+1.d0)*Zeta(i)
        if(v(i).lt.0.1)then

         add1= jfac/(2.*ek(i)-WIntegrale)
         rem1= jfac/(2.*ek(i)+WIntegrale)
         DispAdd(ir)=DispAdd(ir)+jfac/(2.*ek(i)-WIntegrale)
         DispRem(ir)=DispRem(ir)+jfac/(2.*ek(i)+WIntegrale)

        else
         add1= jfac/(2.*ek(i)+WIntegrale)
         rem1= jfac/(2.*ek(i)-WIntegrale)
         DispAdd(ir)=DispAdd(ir)+ jfac/(2.*ek(i)+WIntegrale)
         DispRem(ir)=DispRem(ir)+ jfac/(2.*ek(i)-WIntegrale)
        endif

         if(abs(Wintegrale-0.).lt.0.01)write(99,120)Wintegrale,ek(i),jjk(i),add1,rem1,DispAdd(ir),DispRem(ir)
        enddo

 120      format(2(f8.4),2x,i5,6(f10.4,2x))
 130      format(6(e12.5,2x))

 
         diff_add= abs(Wintegrale-PairAddition)
         diff_rem= abs(Wintegrale-PairRemoval)

         if(diff_add.lt.emin_add)then
      
          emin_add= diff_add
          GAdd=1/DispAdd(ir)
!          write(*,*)'GAdd=',Gadd, 'PairAddition =', PairAddition
         endif

         if(diff_rem.lt.emin_rem)then
                emin_rem=diff_rem
          GRem=1/DispRem(ir)
!          write(*,*)'GRem=',GRem, 'PairRemoval =', PairRemoval
         endif

         DispersionFunction (ir) =  DispAdd(ir)
         DispersionFunction (-ir) =  DispRem(ir)
!         if(abs(DispAdd(ir)).lt.5..and.abs(DispRem(ir)).lt.5.)then
!          write(10,145)WIntegrale,DispAdd(ir),-Wintegrale,DispRem(ir)

!         endif

!         if(abs(WIntegrale-2*ek(5)).lt.0.01)
!     1      write(*,*)'polo',2*ek(5),DispRem(ir)

!         if(abs(WIntegrale-2*ek(6)).lt.0.01)
!     1      write(*,*)'polo',2*ek(6),DispAdd(ir)

         Wintegrale=Wintegrale+dW
         ir=ir+1
       enddo

       end subroutine disp_relation
!----------------------------------------------!

!----------------------------------------------!
       subroutine XY1_Calculation

       implicit none
       integer :: i

       do i=1,NLivelli

        SumAdd=0.
        SumRem=0.

        jfac= sqrt(0.5*(jjk(i)+1.)*Zeta(i))
!         jfac=1.
        if(V(i).gt.0.1)then
         Xadd(i)=0.
         Yrem(i)=0.
         Xrem(i)= jfac/(2*ek(i)-PairRemoval)
         Yadd(i)= jfac/(2*ek(i)+PairAddition)
         !write(*,*)Xrem(i),ek(i),PairAddition
         SumAdd= SumAdd- Yadd(i)*Yadd(i)
         SumRem= SumRem + Xrem(i)*Xrem(i)
        else
         Xrem(i)=0.
         Yadd(i)=0.
         Xadd(i)= jfac/(2*ek(i)-PairAddition)
         Yrem(i)= jfac/(2*ek(i)+PairRemoval)
!         write(*,*)'why?',jfac,Yrem(i),ek(i),PairRemoval
         SumAdd= SumAdd+ Xadd(i)*Xadd(i)
         SumRem= SumRem- Yrem(i)*Yrem(i)
        endif
       enddo

       SumAdd=0.
       SumRem=0.
       do i=1,NLivelli
!        jfac= (0.5*(jjk(i)+1.))
        jfac=1.d0
         SumAdd=SumAdd+jfac*Xadd(i)**2-jfac*Yadd(i)**2.
         SumRem=SumRem+jfac*Xrem(i)**2-jfac*Yrem(i)**2
       enddo

        if(SumAdd.lt.0.d0.or.SumRem.lt.0.d0)then
           write(*,*)
           write(*,*)'============================================='
           write(*,*)'==== WARNING : Sum X^2-Y^2 not positive! ===='
           write(*,*)'============================================='
           write(*,*)'SumAdd,SumRem=',SumAdd, SumRem
           write(*,*)
        endif

        cadd= 1./sqrt(abs(SumAdd))
        crem= 1./sqrt(abs(SumRem))


!       SumAdd=sqrt(abs(SumAdd))*SumAdd/abs(SumAdd)
!       SumRem=sqrt(abs(SumRem))*SumRem/abs(SumRem)
       write(*,*)'costanti di Accoppiamento'
       write(*,*)'Lambda, Addition Mode=', cadd
       write(*,*)'Lambda, Removal  Mode=', crem

       write(*,*)
       write(*,*)'Coefficienti X_n e Y_n'

       write(11,*)'    E                L    2J      Z           X         Y'
       do i=1,Nlivelli
        Xadd(i)=Xadd(i)*cadd
        Yadd(i)=Yadd(i)*cadd
        Xrem(i)=Xrem(i)*crem
        Yrem(i)=Yrem(i)*crem

        if(V(i).gt.0.1)then
          write(*,*)'==================',i,'       ================='
          write(*,*)'X_rem = ',Xrem(i),' Y_add = ',Yadd(i)
          write(11,105)e_sp(i),ek(i),llk(i),jjk(i),Zeta(i),XRem(i),YAdd(i)
        else
          write(*,*)'==================',i,'     ================='
          write(*,*)'X_add = ',Xadd(i),' Y_rem = ',Yrem(i)
          write(11,105)e_sp(i),ek(i),llk(i),jjk(i),Zeta(i),XAdd(i),YRem(i)
        endif
       enddo
 105      format(2(f8.4),2x,2(i5),6(f10.4,2x))

       end subroutine XY1_Calculation
!--------------------------------------------!

       subroutine XY(Epho)

       implicit none
       real*8 :: Epho
       integer :: i
       
       write(*,*)'Calculating XY for Phonon', Epho

       do i=1,NLivelli

        jfac= sqrt(0.5*(jjk(i)+1.)*Zeta(i))
!         jfac=1.
        if(Epho.lt.0.)then
          if(V(i).gt.0.1)then
            Yrem(i)=0.
            Xrem(i)= jfac/(2*ek(i)-Abs(Epho))
          else
            Xrem(i)=0.
            Yrem(i)= jfac/(2*ek(i)+abs(Epho))
          endif
        else
         if(V(i).gt.0.1)then
         Xadd(i)=0.
         Yadd(i)= jfac/(2*ek(i)+Epho)
        else
         Yadd(i)=0.
         Xadd(i)= jfac/(2*ek(i)-Epho)
        endif
        endif
       enddo

         SumAdd=0.
         SumRem=0.
       if(Epho.lt.0.)then
         do i=1,NLivelli
!          jfac= (0.5*(jjk(i)+1.))
           jfac=1.d0
           SumRem=SumRem+jfac*Xrem(i)**2-jfac*Yrem(i)**2
         enddo
       else
         do i=1,NLivelli
!          jfac= (0.5*(jjk(i)+1.))
           jfac=1.d0
           SumAdd=SumAdd+jfac*Xadd(i)**2-jfac*Yadd(i)**2.
         enddo
       endif
       
        if(SumAdd.lt.0.d0.or.SumRem.lt.0.d0)then
           write(*,*)
           write(*,*)'============================================='
           write(*,*)'==== WARNING : Sum X^2-Y^2 not positive! ===='
           write(*,*)'==== WARNING : for phonon',    Epho,   ' ===='
           write(*,*)'============================================='
           write(*,*)'SumAdd,SumRem=',SumAdd, SumRem
           write(*,*)
        endif

        cadd= 1./sqrt(abs(SumAdd))
        crem= 1./sqrt(abs(SumRem))

       do i=1,Nlivelli
        Xadd(i)=Xadd(i)*cadd
        Yadd(i)=Yadd(i)*cadd
        Xrem(i)=Xrem(i)*crem
        Yrem(i)=Yrem(i)*crem

        if(Epho.lt.0.)then
          write(12,105)e_sp(i),ek(i),llk(i),jjk(i),Zeta(i),Xrem(i),Yrem(i)
        else
          write(12,105)e_sp(i),ek(i),llk(i),jjk(i),Zeta(i),Xadd(i),Yadd(i)
        endif
       enddo
       
 105      format(2(f8.4),2x,2(i5),6(f10.4,2x))
 
       end subroutine XY

       end module Stuff

!------------------- MAIN -------------------!
       program PairVibrations
       use Stuff

       implicit none
       integer             :: i, icalc, Znucl,inuc
       real (kind=8)       :: Sep_hole,Sep_part,Amass, PairingGap
       logical             :: flag_efn = .false.

       
       ir_tot = ceiling(EcutFunction/dW)
       allocate (DispersionFunction(-ir_tot:ir_tot))

       open(unit=10,file='Dispersione.dat')
       open(unit=11,file='Factors.dat')
       open(unit=12,file='RPA.dat')
       open(unit=15,file='Input_WS.in',status='old')
       open(unit=20,file='mass.mas12', status='old',action = 'read')
       open(unit=21,file='rct2.mas12', status='old',action = 'read')
       open(unit=51,file='outlevels.dat')

       
       read(15,*)Amass,Znucl,    inuc

       call read_mass(int(Amass),Znucl,inuc,DiffAddition,DiffRemoval,Sep_hole,Sep_part,PairingGap)
       
       write(*,*)'-- Separation Energy of hole and particle states =', Sep_hole,Sep_part
       write(*,*)'- Correlation Energy of Addition and Removal =', DiffAddition, DiffRemoval
       write(*,*)'------------- Pairing Gap ------------------ =',PairingGap
       
       write(*,*)'New = 2, Livelli Sperimentali = 1, Livelli Self Energy = 0'
       read(*,*) icalc

       if(icalc.le.1)then
         call read_define(icalc,flag_efn)
       elseif(icalc.eq.2)then
         call make_levels(Amass,Znucl,inuc,Sep_hole,Sep_part)
       else
         write(*,*)'wrong starting option icalc'
         stop
       endif


       call disp_relation

       Wintegrale=-EcutFunction+dW
       do i=-ir_tot+2,ir_tot-2
          if(DispersionFunction(i).lt.DispersionFunction(i-1))then
             if(DispersionFunction(i).lt.DispersionFunction(i+1))then
                if(abs(DispersionFunction(i-1)-DispersionFunction(i+1)).lt.10*dW)then
                   Wminimo=Wintegrale
                   write(*,*)'Wminimo Trovato!!',Wminimo
                endif
             endif
          endif

          Wintegrale=Wintegrale+dW
       enddo

       efn = efn + Wminimo/2.d0
       flag_efn = .true.

       write(*,*)'Calculating with adjusted Fermi Energy'

       call disp_relation
       
       call define_minima

       call XY1_Calculation

       END

!!! Pb208 !!!
!        DiffAddition=1.237
!        DiffRemoval= 0.64

!!! Dummy !!!
!       DiffAddition=1.5
!       DiffRemoval=0.5

!!! Sn132 !!!
!       DiffAddition=1.17
!       DiffRemoval=2.14

!!! Sn100 !!!
!       DiffAddition=2.92
!       DiffRemoval=5.15

!!! Li9  !!!
!       DiffAddition=0.49
!       DiffRemoval=2.02

!!! Be10 !!!
!       DiffAddition=2.11
!       DiffAddition=2.66123
!       DiffRemoval=5.14233
!       DiffRemoval=4.44

!!! O16  !!!
!       DiffAddition=3.9
!       DiffRemoval=2.44

!!! Ca40 !!!
!       DiffAddition=3.12
!       DiffRemoval=2.339

!!! Ca48 !!!
!       DiffAddition=1.214
!       DiffRemoval =2.6734

     subroutine define_minima
     use Stuff
     implicit none
     
     integer :: i
     
       write(12,*)'Pairing strength of Pair Addition and Removal phonons'
       write(12,*)'G(a=+2)=',GAdd,'G(a=-2)=',GRem
       write(12,*) 
       
       Wintegrale=-EcutFunction+dW
       do i=-ir_tot+2,ir_tot-2
         if(i.lt.0)then   !E<0
           if(abs(DispersionFunction(i)- 1.d0/Grem).lt.abs(DispersionFunction(i-1)- 1.d0/Grem).and. &
              abs(DispersionFunction(i)- 1.d0/Grem).lt.abs(DispersionFunction(i+1)- 1.d0/Grem).and. &
          abs(abs(DispersionFunction(i-1)-1.d0/Grem) - abs(DispersionFunction(i+1)- 1.d0/Grem)).lt.10.*dW                        &
              )then
               write(12,*)'Ephonon-Rem',Abs(Wintegrale)
               write(12,*)'    E                L    2J      Z           X         Y'   
               call XY(Wintegrale)  
            endif
         else             !E>0
            if(abs(DispersionFunction(i)- 1.d0/Gadd).lt.abs(DispersionFunction(i-1)- 1.d0/Gadd).and. &
               abs(DispersionFunction(i)- 1.d0/Gadd).lt.abs(DispersionFunction(i+1)- 1.d0/Gadd).and. &
           abs(abs(DispersionFunction(i-1)-1.d0/Gadd) - abs(DispersionFunction(i+1)- 1.d0/Gadd)).lt.10.*dW                        &
              )then
               write(12,*)'Ephonon-Add',Abs(Wintegrale)
               write(12,*)'    E                L    2J      Z           X         Y'    
               call XY(Wintegrale)           
            endif
         endif
         
         
! Fermi Energy Check
          if(DispersionFunction(i).lt.DispersionFunction(i-1))then
             if(DispersionFunction(i).lt.DispersionFunction(i+1))then
                if(abs(DispersionFunction(i-1)-DispersionFunction(i+1)).lt.10.*dW)then
                   Wminimo=Wintegrale
                   write(*,*)'Wminimo Trovato!!',Wminimo
                   if(abs(Wminimo).gt.dW)write(*,*)'!!!WARNING!!!! Fermi Energy Not Correct!!!'
                endif
             endif
          endif

          write(10,*)Wintegrale,DispersionFunction(i)
          Wintegrale=Wintegrale+dW
       enddo

       write(*,*)'Caratteristiche dei Livelli'
       do i=1,Nlivelli
        write(*,*)llk(i),jjk(i),ek(i),Zeta(i)
       enddo

       write(*,*)
       write(*,*)
       write(*,*)'---- Fermi Energy=', efn,' ----'
       write(*,*)'Energy of Pair Addition and Removal phonons'
       write(*,*)'W(a=+2)=',PairAddition,'W(a=-2)=',PairRemoval
       write(*,*)'Pairing strength of Pair Addition and Removal phonons'
       write(*,*)'G(a=+2)=',GAdd,'G(a=-2)=',GRem

       write(11,*)'---- Fermi Energy=', efn,' ----'
       write(11,*)'Energy of Pair Addition and Removal phonons'
       write(11,*)'W(a=+2)=',PairAddition,'W(a=-2)=',PairRemoval
       write(11,*)'Pairing strength of Pair Addition and Removal phonons'
       write(11,*)'G(a=+2)=',GAdd,'G(a=-2)=',GRem
       write(11,*)
end subroutine define_minima

!---------------------------------!
!Make the levels launching a Wood Saxon and readjusting to Separation Energies
     subroutine make_levels(Amass,Znucl,inuc, Sep_hole,Sep_part)
       use Stuff
       implicit none
       integer :: i,j,ir, Pnum
       real*8  :: y,Sep_part,Sep_hole,Diff_Sep
       !------ Input Wood Saxon -------!
       integer :: Znucl,inuc
       real*8  :: Amass,r0,an0,rs0,Vson,EFP !A, Z, r0, a0, rs0 di Wood Saxon
       real*8  :: Vrn_n,Vrn_p
       !output
       DIMENSION PS(mphon,Nbox)
       DIMENSION DVN(mphon)

       real*8  :: PS, DVN
       !-------------------------------!
       Zeta = 1.d0
    !Wood Saxon Definition B&M (2-182)

       r0 = 1.27;rs0 = 1.27; an0 = 0.67

       if(inuc.eq.1)then
          write(*,*)'Protons Pairing Vibration'
          Vrn_p =  -51.d0 - 33.d0*(1.d0*Amass - 2.d0*Znucl)/Amass    
          Vson  = 0.44d0*Vrn_p   

       call boundary(mphon, mphon, dr,                                &
                     Amass,Znucl,inuc,r0,rs0,An0,Vson,Vrn_p,           &
                        0,lmax,EcutFunction,Nbox,                     &
              PS,          Nlivelli,jjk,llk,e_sp,DVN)

       else
          write(*,*)'Neutrons Pairing Vibration'
   
          Vrn_n =  -51.d0 + 33.d0*(1.d0*Amass - 2.d0*Znucl)/Amass
          Vson  = 0.44d0*Vrn_n 
       call boundary(mphon, mphon, dr,                                &
                     Amass,Znucl,inuc,r0,rs0,An0,Vson,Vrn_n,           &
                        0,lmax,EcutFunction,Nbox,                     &
              PS,          Nlivelli,jjk,llk,e_sp,DVN)

       endif
                
             !WaveFunction, N, J, L , E  , V, dim vectors1, dim vectors 2 
       do i=1,50
          write(51,*)e_sp(i),float(llk(i)),float(jjk(i))
       enddo
       
       !Sorting Secondo Energia
         do i=1,NLivelli-1
           do j=i+1,NLivelli
             if(E_sp(i).gt.E_sp(j))then
   	           y=e_sp(i)
	           e_sp(i)=e_sp(j)
	           e_sp(j)=y
               
               y=llk(i)
               llk(i)=llk(j)
               llk(j)=y
               
               y=jjk(i)
               jjk(i)=jjk(j)
               jjk(j)=y
               
               DO IR=1,Nbox
                 y=PS(i,IR)
                 PS(i,IR)=PS(j,IR)
                 PS(j,IR)=y
               END DO
  	         endif
          enddo       
         enddo
       
       Pnum = 0
       V = 0.d0
        do i=1,NLivelli
          V(i) = 1.d0
          Pnum = Pnum + jjk(i)+1    
          write(*,*)llk(i),jjk(i),e_sp(i),Pnum
          
          if    (inuc.eq.0 .and. Pnum.ge.Amass-Znucl)then !neutrons
             iGapDn=i;iGapUp=i+1
             exit
          elseif(inuc.eq.1 .and. Pnum.ge.      Znucl)then !protons
             iGapDn=i;iGapUp=i+1          
             exit
          endif

        enddo
        
        write(*,*)'Last  Occupied Level ek=',e_sp(iGapDn),'Sep_energy=',Sep_hole
        write(*,*)'First Empty    Level ek=',e_sp(iGapUp),'Sep_energy=',Sep_part
        
        Diff_Sep = Sep_hole - e_sp(iGapDn)
        do i=1,iGapDn
           e_sp(i) = e_sp(i)+Diff_Sep
        enddo
        
        Diff_Sep = Sep_part - e_sp(iGapUp)

        do i=iGapUp,NLivelli
           e_sp(i) = e_sp(i)+Diff_Sep
        enddo
        
        efn =  (e_sp(iGapDn)+e_sp(iGapUp))/2.d0
        
        write(*,*)'Fermi Energy           =',efn
        
        do i=1,NLivelli
         ek(i)=abs(e_sp(i)-efn)
         if(e_sp(i).gt.EcutFunction)exit
        enddo
        NLivelli = i
      
     end subroutine
     
!--------------------------------------------------------------------!

!Read Masses from AME 2003 and following, remember to remove # (substitute with .)
!Then calculates DiffAddition and DiffRemoval

!Input: Amass, Znucl = A, Z of the nucleus, 
!       inucl = 0, for neutron calculation, = 1 for proton calculation

!Output: DiffAddition,DiffRemoval = Correlation Energy, Positive
!        Sep_hole,Sep_part = Energy of Hole and Particle state (thus Separation of nucleons from A and A+1 Nucleus), negative.
!        PairingGap = from 3-point formula, Positive.

     subroutine read_mass(Amass,Znucl,inuc,DiffAddition,DiffRemoval,Sep_hole,Sep_part,PairingGap)
       implicit none
       integer                          :: Atable,Ztable,Amass,Znucl,inuc
       real (kind = 8) ,dimension (5)   :: BindingE
       real (kind = 8)                  :: DiffAddition,DiffRemoval,Sep_hole,Sep_part,PairingGap
       
       
       character(LEN=9999) :: parkChar
       integer             :: i, parkInt        
       real (kind = 8)     :: parkReal, ParkReal1,parkReal2


       do while(.true.)
          read(20,*,iostat=parkInt)parkChar
          if(index(trim(parkChar),'N-Z').gt.0)then
            !write(*,*)'found'
            exit
          endif
          if(parkInt < 0)then
             write(*,*)'End of File. Cannot Find Masses, Check Input mass.mas12'
             stop
          endif
       end do
       read(20,*)

       do i=1,5
         do while (.true.)
           read(20,654)parkChar,ParkInt,ParkInt,Ztable,Atable, parkChar, parkChar,ParkReal, ParkReal, BindingE(i)
           if(inuc.eq.0 .and. Atable.eq.Amass-3+i .and. Ztable.eq.Znucl    )then !neutron pairing Vibration, search for A+-2 and Z fixed
             exit
           endif        
           if(inuc.eq.1 .and. Atable.eq.Amass-3+i .and. Ztable.eq.Znucl-3+i)then !proton pairing Vibration, search for Z+-2 and N fixed
             exit  !if has the right mass and charge exit the inner reading cycle and associate BindingE(i) and goes to next i in the outside cycle
           endif
         enddo
       enddo
       
654 format (a1,i3,i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5) 

       do i=1,5
         BindingE(i) = BindingE(i)*(Amass-3+i)/1000.d0  !binding energy of nucleus in MeV
!         write(*,*)(Amass-3+i),BindingE(i)
       enddo
       
           !addition: B(A+2)-B(A)-2(B(A+1)-B(A))
       DiffAddition = BindingE(5) - BindingE(3) - 2.d0*(BindingE(4)-BindingE(3))
           !removal : B(A-2)-B(A)-2(B(A-1)-B(A))
       DiffRemoval  = BindingE(1) - BindingE(3) - 2.d0*(BindingE(2)-BindingE(3))
           !Delta   = -1/2 (B(A-1)+B(A+1)-2B(A))
       PairingGap = - ( BindingE(2) + BindingE(4) - 2.d0*BindingE(3) )/2.d0

       do while(.true.)
          read(21,'(a10)',iostat=i)parkChar
          if(index(trim(parkChar),'1 A').gt.0)then
            !write(*,*)'found'
            exit
          endif
          if(i < 0)then
             write(*,*)'Cannot Find Separation Energies, Check Input rct2.mas12'
             stop
          endif
       end do
       
       do while (.true.) !Seach for separation energies Sep_hole,Sep_part
         read(21,'(a1,i3,a3,i5,1x,a40)')parkChar,Atable,parkChar,Ztable,parkChar!,parkReal
         if(inuc .eq. 0) then
           if(Atable.eq.Amass .and. Ztable.eq.Znucl)then 
            if(index(trim(parkChar(1:20) ),'*').le.0)then
              read(parkChar(1:20),*)Sep_hole
            else
              write(*,*)'Neutron Separation Energy non defined, cannot calculate on the dripline'
              stop
            endif           
           endif
           if(Atable.eq.Amass+1 .and. Ztable.eq.Znucl)then 
            if(index(trim(parkChar(1:20) ),'*').le.0)then
              read(parkChar(1:20),*)Sep_part
              exit
            else
              write(*,*)'Neutron Separation Energy non defined, cannot calculate on the dripline'
              stop
            endif
           endif
                       
         elseif(inuc .eq. 1) then
           if(Atable.eq.Amass .and. Ztable.eq.Znucl)then
            if(index(trim(parkChar(20:40) ),'*').le.0)then
              read(parkChar(20:40),*)Sep_hole
            else
              write(*,*)'Proton Separation Energy non defined, cannot calculate on the dripline'
              stop
            endif
           endif
           if(Atable.eq.Amass+1 .and. Ztable.eq.Znucl+1)then 
            if(index(trim(parkChar(20:40) ),'*').le.0)then
              read(parkChar(20:40),*)Sep_part
              exit
            else
              write(*,*)'Proton Separation Energy non defined, cannot calculate on the dripline'
              stop
            endif
           endif
           
         endif        
       enddo

       Sep_hole = -Sep_hole/1000.d0; Sep_part = -Sep_part/1000.d0
       !write(*,*)Atable,Ztable,Sep_hole,Sep_part

      end subroutine read_mass
      

!*****************************************************************************************
!*****************************************************************************************
      subroutine  boundary(nnp, imaxlev, h,              &
                     Amass,Ipro,inuc,r0,rs0,An0,Vso,Vrn,            &
                        Lmin,Lmax,Ecut,Nmaxt,                   &
              PS,          NNST,JJP,LP,EP,VN)
             !WaveFunction, N, J, L , E  , V, dim vectors1, dim vectors 2 

      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: NLMAX=1000

      integer, PARAMETER :: Nmax=1000

      INTEGER :: S,Nmaxt,Ipro,Ishrod
      DIMENSION U(NNP),HMEN(NNP),VN(NNP),VSON(NNP)
      REAL*8  :: DHMEN(NNP),D2HMEN(NNP),GI(NNP),GIprimo(NNP),MESMN(NNP)
      real*8  :: Pi

      DIMENSION LPS(IMAXLEV) 
      DIMENSION RHO(NNP),VC(NNP) 
      DIMENSION PS(imaxlev,Nnp),EP(imaxlev)
      DIMENSION LP(imaxlev),JJP(imaxlev),EE(NMax),DV(NNP)

      DIMENSION Nmin(0:NLmax+1,2),NmaxV(0:NLmax+1,2)  ! da buttare

      common /Interi/iterasaki
      common /Reali/xmk,efn

!       write(*,*)Amass,Ipro,r0,an0,rs0 !A, Z, r0, a0, rs0 di Wood Saxon
!       write(*,*)Vrn,Vso
!       write(*,*)EFA
!       write(*,*)Nmaxt,dr            !Np = dimensione box
!       write(*,*)Ecut,Lmin,Lmax        ! Cut-off single particle energy spectrum,Momento angolare L minimo e massimo


      emax0=2.*ecut

      N2max=Nmax
      Pi=2.*ASIN(1.)
      
      if(inuc.eq.0)open(unit=9,file='Density_neutroni.dat')  ! file di scrittura
      if(inuc.eq.1)open(unit=9,file='Density_protoni.dat')
      

!-- Controllo su Nmax, basato sui numeri di nodi nel caso di buca qudrata,
!-- Landau pg 88 "Meccanica Quantistica"

       Ntrial=int(abs(sqrt(Ecut*(Nmaxt*h)**2/(20.75*3.14**2))))

       If(N2max.lt.Ntrial) then
         write(*,*)'Aumentare il numero di nodi nella subroutine'
         stop
       endif
       if(Nmaxt.gt.nnp) then
          write(*,*)'Errore dimensione vettori subroutine sp'
          stop
       endif

      RN0=r0*(AMass)**(1./3.)
      RNs0=rs0*(AMass)**(1./3.)

      I_letto=0
       do j=1,nmaxt  
           R=h*j   
           TEMP=(R-RN0)/AN0 
           EX=EXP(-TEMP) 
           VN(j)=VRN*EX/(1+EX)                            ! Wood Saxon 

           TEMP=(R-RNS0)/AN0 
           EX=EXP(-TEMP) 
           VSON(j)=VSO*Rs0**2*EX/(1+EX)/(R*AN0*(1.+EX)) 
!
!           write(100,*)VSON(j),VSO,RNS0,EX/(1+EX)/(R*AN0*(1.+EX)) 
!C    calcola il potenziale di una sfera omogenea carica di raggio
          if(inuc.eq.1) then
          if((R-RN0).le.0) then
            VC(j)=((0.71995*(IPRO-1))/RN0)*(3.-(R/RN0)**2)
          else
             VC(j)=1.4399*(IPRO-1)/R
          endif

            VN(j)=VN(j)+VC(j) 
          endif

           DV(j)= -VRDER*EX/(1+EX)/(1+EX)/AN0     !potenziale per la induced 
           HMEN(J)=20.75d0
           DHMEN(J)=0.d0
           D2HMEN(J)=0.d0
        RHO(j)=0.d0
        write(100,*)R,VN(j),VRN
       end do
       
       nh=1
      write(51,*)'Scrivo i livelli di particella singola'
      write(51,*)'protoni-> isospin 1 o neutroni -> isospin 0'
      write(51,*)
      write(51,*)'Energia , N , L , J , isospin'
      write(51,*)

       do L=Lmin,Lmax
        
         do S=1,0,-1
           if(s.eq.1)ij=1
           if(s.eq.0)ij=2
          
           J=2*L+2*S-1
           if(L.eq.0.and.S.eq.0)GO TO 300
        
           do N=1,Nmax,1        
               
               do kx=1,nmaxt
                 GI(kx)=DHMEN(kx)/HMEN(kx)
                 GIprimo(kx)=HMEN(kx)*D2HMEN(kx)-(DHMEN(kx))**2
                 GIprimo(kx)= GIprimo(kx)/(HMEN(kx))**2                
               end do

               CALL NUMEROV(h,Nmaxt,Emax0,N,L,J,E,U,VN,VSON, &
                            HMEN,MESMN,GI,GIprimo,DHMEN,D2HMEN,no)
            
             if(no.ne.-1.and.e.lt.ecut)then
               EE(n)=e
               if(n.gt.1)then
                 if(EE(n).eq.EE(n-1))then
                 write(*,*)'Errore grossolano'
                 endif
               endif

!             write(*,*)E,N,L,J,inuc
 !            write(12,111)E,N,L,J,inuc
 111         format(1x,F15.8,1x,I3,1x,I3,1x,I3,1x,I3)

             LP(NH)=L
             EP(Nh)=E
             JJP(Nh)=j !2*l-2*s+1
             LPS(NH)=2*l+IJ-2
             IF(L.EQ.0)LPS(NH)=0 
             do ii=1,nmaxt 
              PS(NH,ii)=U(ii)
             enddo
             nh=nh+1
             else
!             write(*,*)'Non giusto',e,n,l,j,no
             go to 300
             endif

           end do !N           

300      CONTINUE
         end do

       end do

       NNST=NH-1

       IF(NNST.GT.IMAXLEV) THEN 
          WRITE(6,*) 'IMAXLEV IS TOO SMALL: THE SUBROUTINE Boundary' 
       STOP 
       END IF 
       if(inuc.eq.0)write(*,*)' N. OF ACTIVE NEUTRON LEVELS = ', NNST 
       if(inuc.eq.0)write(51,*)' N. OF ACTIVE NEUTRON LEVELS = ', NNST 
       if(inuc.eq.1)write(*,*)' N. OF ACTIVE PROTON LEVELS = ', NNST 
       if(inuc.eq.1)write(51,*)' N. OF ACTIVE PROTON LEVELS = ', NNST 

!      do i=1,nnst
!       if(inuc.eq.1) then
!         do ii=1,nmaxt   
!           VN(ii)=VN(ii)-VC(ii) ! eventually removes coulomb
!         enddo
!       endif
!      enddo

      return

      END subroutine


!*****************************************************************************************
!*****************************************************************************************


       SUBROUTINE NUMEROV(h,Nmaxt,Emax0,N,L,J,EFINAL,U,VN,VSON,    &
                   hMEN,mesmn,GI,GIprimo,DHMEN,D2HMEN,no)

       IMPLICIT INTEGER (I-N)
       IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 MA,MESMP(Nmaxt),MESMN(Nmaxt),POT 
       DIMENSION U(Nmaxt),V(Nmaxt),HMEN(Nmaxt),HME(Nmaxt),         &
       HMEP(Nmaxt),VN(Nmaxt),VP(Nmaxt),VSON(Nmaxt),VSOP(Nmaxt),    &
       VSO(Nmaxt),VC(Nmaxt),VCOUL(Nmaxt)
       REAL*8 GI(Nmaxt),GIprimo(Nmaxt),EFFE,KAPPAQ(Nmaxt),W(Nmaxt),WAUX(Nmaxt)
       REAL*8 DHMEN(Nmaxt),D2HMEN(Nmaxt)     
       REAL*8 DE
       INTEGER kstar,flagstar
       REAL*8 Fstar,Fstarprec

       no=0
       !-------- Description of this subroutine --------

       !1 - setting Emin equal to the minimum of the potential
       !2 - first search for the eigenfunction. It stops when the number of
       !    node is equal to N-1. (iprogress=1 for this step)
       !3 - second search for the eigenfunction. It stops when the number of
       !    node is equal to N. (iprogress=2 for this step). ETRIAL starts from
       !    the last value in the "iprogress=1 phase".
       !4 - checking for errors on node numbers (ND2-ND1 should be equal to 1)
       !5 - finer search for the eigenstate. The shooting method from here after is governed
       !    looking at the sign of the wave function at the edge of the box.
       !    The search stops when the difference between EMIN and EMAX is less than
       !    a fixed parameter DE, assigned at the beginning of the subroutine.
       !6 - cleaning the tail of the solutions, both for bound and unbound states.
       !    For bound states the exponentially decaying behaviour for r-->infinity
       !    is built by propagation of the solution from the right, and matching with the
       !    left part of the wave function. For unbound states, an addition of a linear function
       !    is performed in order to shift the last node to the edge of the box.
       !7 - going back to the physical wave function, by dividing by a factor
       !    sqrt((hbar**2/2m)*mbare/mstar) and normalizing this wave function.
       !    The outcoming wave function is called U(i). Moreover, U(i) is multiplied
       !    by the sign of U(kstar+1)

       !---- defining the accuracy for the energy (in MeV) ----
       DE=0.5E-12
       niter=500
       NODE=N-1  !this is the search number of nodes         
       EMAX=EMAX0!ecut+emax0

       WLIM=1E-07

       istamp=0
       


       iprogress=1  !this index indicates the progress of the calculation 

       do i=1,nmaxt
!         if(nt.eq.1) then
!           V(i)=VP(i)
!           VSO(i)=VSOP(i)
!           VC(i)=VCOUL(i)
!           HME(i)=HMEP(i)
!         elseif(nt.eq.0) then
           V(i)=VN(i)
           VSO(i)=VSON(i)
           HME(i)=HMEN(i)
!         end if
       end do


       ESO=(J*(J+2)-4*L*(L+1)-3.)/8.
       if(L.eq.0)ESO=0.d0

  56   CONTINUE
  
       
       !---- setting Emin equal to the minimum of the potential ----
       POTmin=0.d0
       do i=1,nmaxt
         R=H*I
         POT=(V(i)+ESO*VSO(i))+L*(L+1)*HMEN(i)/R**2.d0
         if(POT.lt.POTmin)POTmin=POT
       end do
       EMIN=POTmin

 59    CONTINUE
       
       DO KT=1,NITER

         ETRIAL = (EMIN+EMAX)/2.0

         if(istamp.eq.1)write(*,*)
         if(istamp.eq.1)write(*,*)'kt,etrial = ',kt,etrial
         if(istamp.eq.1)write(*,306)emin,emax
 301     FORMAT('kt,etrial = ',I5,1(X,F12.10))
 306     FORMAT(' emin,emax = ',16X,2(X,F12.6))
         W(1)=1.E-10
         ND=0
         IOUT=0
          
         !---- determining the point (kstar) before which the wave function
         !     cannot have an oscillating behaviour.
         flagstar=0
         do i=1,nmaxt
         R=H*I
           POT=(V(i)+ESO*VSO(i))+L*(L+1)*HMEN(i)/R**2.d0
           Fstar=POT-Etrial
!           if(L.gt.10)then
             if(flagstar.eq.0.and.Fstar.eq.0.d0)then
               kstar=i
               flagstar=1             
             elseif(i.ge.2)then
               if(flagstar.eq.0.and.Fstarprec*Fstar.lt.0.d0)then
                 kstar=i
                 flagstar=1
               end if
             end if
!           end if
           Fstarprec=Fstar

           EFFE=(ETRIAL-POT)/HMEN(i)
           KAPPAQ(i)=-0.5d0*GIprimo(i)-0.25d0*(GI(i))**2.d0+EFFE
         end do

         if(L.le.1)kstar=1
         if(istamp.eq.1)write(*,302)kstar
         if(istamp.eq.1)write(*,304)iprogress
 302     FORMAT('     (kstar= ',I5,')')
 304     FORMAT('   (iprogress= ',I5,')')
         !---- propagation of the solution till U(nmaxt)       
         DO I = 1,NMAXT-1!NMXT-2
         R=H*I

           !---- to compute W(2) it uses W(1) and W(0)=0.0
           if(i.eq.1)then
             W(i+1)=W(i)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(i))
             W(i+1)=W(i+1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(i+1)) 
           else
             W(i+1)=W(i)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(i))-   &
                  W(i-1)*(1.d0+(H**2.d0/12.d0)*KAPPAQ(i-1))
             W(i+1)=W(i+1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(i+1))
           end if

           !---- counting the nodes beyond kstar, to avoid counting
           !     non-physical nodes
           if(i.ge.kstar)then
             IF(W(I+1)*W(I).LT.0)then
               ND=ND+2
             elseif(W(I+1)*W(I).eq.0)then
               ND=ND+1
             end if
           end if  !it needs to divide ND by 2 (made later)
           

         END DO       

  46     CONTINUE
         ND=ND/2

         if(istamp.eq.1)write(*,305)ND
         if(istamp.eq.1)write(*,308)W(nmaxt)
 305     FORMAT('   (nodes= ',I5,')')     
 308     FORMAT('   (W(nmaxt)= ',E14.6,')')

         !if(EMAX-EMIN.lt.DE) GO TO 61


         !---- shooting method (the node at the edge of the box does not count)
         !     It varies the energy (using bisection) until the number of nodes
         !     is equal to:
         !                  N-1 for iprogress=1
         !                  N   for iprogress=2
         if(iprogress.eq.1)then
           if(ND.eq.NODE)then
             GO TO 61    
           elseif(ND.lt.NODE)then
             EMIN=ETRIAL
           elseif(ND.gt.NODE)then
             EMAX=ETRIAL
           end if
         elseif(iprogress.eq.2)then
           if(ND.eq.NODE+1)then
             GO TO 61    
           elseif(ND.lt.NODE+1)then
             EMIN=ETRIAL
           elseif(ND.gt.NODE+1)then
             EMAX=ETRIAL
           end if
         elseif(iprogress.eq.3)then
           E3=ETRIAL
           WE3=W(nmaxt)*W(kstar)
           if(w(kstar).eq.0.d0)then
             write(*,*)'ERROR for N,L,J=',N,L,J
             write(*,*)'w(kstar)=0.0'
             NO=-1 
             GO TO 100
           end if
           ND3=ND
           if(WE1*WE3.gt.0.d0)then     !solution 3 has W(nmaxt) of the same sign as for
                                       !solution 1. That is, their energies are both HIGHER
                                       !than the eigenstate.
             ND1=ND3
             WE1=WE3
             E1=E3
             GO TO 36  !when iprogress=3 the shooting method decisions are taken
                       !outside the cycle on KT. 
             
           elseif(WE1*WE3.lt.0.d0)then !solution 3 has W(nmaxt) of the same sign as for
                                       !solution 2. That is, their energies are both LOWER
                                       !than the eigenstate.

             ND2=ND3
             WE2=WE3
             E2=E3
             do i=1,nmaxt
               Waux(i)=W(i)  !memorizing solution 2 (ND=N+1)
             end do
             GO TO 36  !when iprogress=3 the shooting method decisions are taken
                       !outside the cycle on KT.
           elseif(WE1*WE3.eq.0.d0)then !since WE1 is NOT equal to zero, this means that
                                       !WE3=0.d0. In other words, solution 3 is the
                                       !eigenstate we are searching for.
             EFINAL=E3
             NDfinal=ND3
             do i=1,nmaxt
               Waux(i)=W(i)  !memorizing solution 2 (ND=N+1)
             end do
             GO TO 101
           end if
         end if


  60     CONTINUE          
       end do  !KT

  61   CONTINUE

       if(KT.eq.NITER+1)then  !KT ran 1 through NITER and a solution
                              !has not been found. The energy of the
                              !searched eigenstate is outside the energy range
                              ![EMIN,EMAX]
!         write(*,*)'Error for N,L,J=',N,L,J
!         write(*,*)'SOLUTION NOT FOUND AT iprogress=',iprogress
         NO=-1 
         GO TO 100
       end if

       if(iprogress.eq.1)then
         ND1=ND            !memorizing number of nodes ND,wave function at the edge, energy
         WE1=W(nmaxt)*W(kstar)
         if(w(kstar).eq.0.d0)then
           write(*,*)'ERROR for N,L,J=',N,L,J
           write(*,*)'w(kstar)=0.0'
           NO=-1 
           GO TO 100
         end if
         E1=ETRIAL
         if(istamp.eq.1)write(*,*)
         if(istamp.eq.1)write(*,*)'  FOUND:N,L,J,KT=',N,L,J,KT
         if(istamp.eq.1)write(*,303)
 303     FORMAT(50('-'))
         iprogress=iprogress+1

         GO TO 59
       elseif(iprogress.eq.2)then
         ND2=ND            !memorizing number of nodes ND,wave function at the edge, energy
         WE2=W(nmaxt)*W(kstar)
        if(w(kstar).eq.0.d0)then
           write(*,*)'ERROR for N,L,J=',N,L,J
           write(*,*)'w(kstar)=0.0'
           NO=-1 
           GO TO 100
         end if
         E2=ETRIAL
         if(istamp.eq.1)write(*,*)
         if(istamp.eq.1)write(*,*)'  FOUND:N,L,J,KT=',N,L,J,KT
         if(istamp.eq.1)write(*,303)
         do i=1,nmaxt
           Waux(i)=W(i)  !memorizing solution 2
         end do
       end if

       !---- checking for errors
       !     ND2 should be ND1+1
 36    if(ND2-ND1.gt.1.or.ND2-ND1.lt.1)then
!         write(*,*)'Error for N,L,J=',N,L,J
!         write(*,*)'ND2-ND1=',ND2-ND1
         NO=-1 
         GO TO 100
       end if

       !---- stepping to iprogress=3
       !     The shooting method from here after is governed
       !     looking at the sign of the wave function at the edge of the box.
       iprogress=3
       if(istamp.eq.1)write(*,*)' iprogress = 3'
       if(istamp.eq.1)write(*,*)' ND1,E1= ',ND1,E1
       if(istamp.eq.1)write(*,*)' ND2,E2= ',ND2,E2
       if(WE1*WE2.gt.0.d0)then
!         write(*,*)'Error for N,L,J=',N,L,J
!         write(*,*)'WE1*WE2<0'
         NO=-1 
         GO TO 100         
       elseif(WE1*WE2.eq.0.d0.and.WE1.eq.0.d0)then
!         write(*,*)'Error for N,L,J=',N,L,J
!         write(*,*)'WE1=0'
         NO=-1 
         GO TO 100
       elseif(WE1*WE2.eq.0.d0.and.WE1.ne.0.d0)then  !solution found in iprogress=2
                                                    !is the right solution
         EFINAL=E2 
         NDfinal=ND2
         GO TO 101
       elseif(WE1*WE2.lt.0.d0)then  !we have to keep on playing the shooting method game
                                    !until the requested accuracy on the energy is required.         
         if(EMAX-EMIN.le.DE)then
           EFINAL=E2                                 
           NDfinal=ND2
           GO TO 101
         else 
           EMIN=E1                  !the solution is in between solution 1 (with ND=N) 
           EMAX=E2                  !and solution 2 (with ND=N+1)
           GO TO 59  !let's propagate the solution in the same way used for first and second phases.
                     !A difference occurs: the shooting method is driven by the sign of the
                     !wave function at the edge of the box.

         end if
       end if

       !---- last part of the game: the cleaning of the tail of the wave function
 101   CONTINUE

       do i=1,nmaxt
         W(i)=Waux(i)  !solution 2 (=the one with ND=ND+1), which is memorized
                       !in Waux, is now recorded in W(i)
       end do

       !---- cleaning the tail of wave functions
       ! A   finding the last node of the wave function (coordinate=ii+1 or ii+2)
       do i=nmaxt,kstar,-1
         if(W(i-1)*W(i).lt.0.d0)then
           ii=i-1
           GO TO 70
         elseif(W(i-1)*W(i).eq.0.d0)then
           ii=i-2
           GO TO 70
         end if
       end do
 70    CONTINUE

       if(istamp.eq.1)write(*,*)'  last node at ii=',ii
       ! B   for bound states: propagation of the solution from the edge of the box for i>=ii
       !                       and subsequent matching to the left part (i<ii)
       !     for unbound states: addition of a linear function for i>=ii, in order
       !                         to shift the last node to the edge of the box.
       if(EFINAL.lt.0.d0)then

         !searching for the position (=imax) of the first maximum/minimum 
         !and computing the average point (=KL) between this max/min and
         !the last node (=ii). This point (=KL) is found in the tail of
         !last oscillation of the wave function, before of the (expected)
         !exponentially decaying behaviour.
         do imax=ii,kstar,-1   
           if(dabs(W(imax-1))-dabs(W(imax)).le.0.d0)then
             KL=(imax+ii)/2
             GO TO 71
           end if
         end do
 71      CONTINUE
         if(istamp.eq.1)write(*,*)'  first max/min at KL=',imax
         if(istamp.eq.1)write(*,*)'  first ramp at KL=',KL
         Wtemp=W(KL)           !memorizes the actual value at KL
         W(nmaxt)=0.d0
         W(nmaxt-1)=1.E-10
         k4=nmaxt-KL           !this is the number of points between Kl and the edge of the box
         do k3=2,k4
           k=nmaxt-k3+1        !k runs from nmaxt-1 to nmaxt-(nmaxt-KL)+1=KL+1
           R=H*I
           !---- propagation of the solution from the right to k-1=KL
           W(k-1)=W(k)*2.d0*(1.d0-(5.d0*H**2.d0/12.d0)*KAPPAQ(k))-     &
                  W(k+1)*(1.d0+(H**2.d0/12.d0)*KAPPAQ(k+1))
           W(k-1)=W(k-1)/(1.d0+(H**2.d0/12.d0)*KAPPAQ(k-1))
         end do
         !---- rescaling the right part in order to match the left part at KL
         FAC=Wtemp/W(KL)
         do k=KL,nmaxt
           W(k)=W(k)*FAC
         end do
       elseif(EFINAL.ge.0.d0)then
         if(istamp.eq.1)write(*,*)'  again ii=',ii
         do k=ii,nmaxt
           W(k)=W(k)-dfloat(k-ii)*W(nmaxt)/dfloat(nmaxt-ii)
         end do
       end if

       !---- computing U(i)=W(i)/sqrt(HMEN(i))
       do i=1,nmaxt
         U(i)=W(i)/(dsqrt(HMEN(i)))
       end do                   

       !---- normalization of the radial wave function
       !     (employing the extended Simpson's rule (eq.4.1.13 Fortran Numerical Recipes))
       SUM=0.0d0
       do i=1,nmaxt!-1
         if(i.eq.nmaxt)then!! -1
           SUM=SUM+(1.d0/3.d0)*U(i)**2.d0
           GO TO 81
         elseif(INT(i/2)*2.eq.i)then
           SUM=SUM+(2.d0/3.d0)*U(i)**2.d0
         elseif(INT(i/2)*2.ne.i)then
           SUM=SUM+(4.d0/3.d0)*U(i)**2.d0
         end if

 81      CONTINUE
       end do
       SUM=DSQRT(H*SUM)

       !---- giving a sign to the wave function, multiplying
       !     by the sign of U(kstar+1)
       SIG=U(kstar+1)/dabs(U(kstar+1))
       do i=1,nmaxt
         U(i)=U(i)*SIG/SUM
       end do

       !---- checking for divergent behaviour of the wave function
       !     at the edge of the box
       ABSWF=0.d0
       do i = 1,nmaxt
         if(EFINAL.lt.0.d0.and.dabs(U(i)).gt.ABSWF)then
           IABSWF=i
           ABSWF=dabs(U(i))
         end if
       end do

       IF (IABSWF.eq.nmaxt.and.EFINAL.lt.0.d0)then
!         write(*,*)'ERROR IN WAVEFUNCTIONS FOR BOUND STATE'
!         write(*,*)'N,L,J,E',N,L,J,EFINAL
         GO TO 100
       end if

       RETURN  !exiting without errors

       !---- exiting with error: use NO outside the subroutine
       !     to check exiting status
 100   NO=-1
       RETURN

       END   
