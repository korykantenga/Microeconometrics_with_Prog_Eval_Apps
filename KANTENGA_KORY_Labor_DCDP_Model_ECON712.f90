! ----------------------------------------------------------------------
! The routines for random number generation were obtained from
! the book NUMERICAL RECIPES


SUBROUTINE ran1sub(idum,x)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum
    REAL(8) :: x
    INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
    REAL(8), SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/IM
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    end if
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    if (iy<0) iy=iy+IM
    x=am*ior(iand(IM,ieor(ix,iy)),1)
END SUBROUTINE ran1sub


SUBROUTINE ran1(idumtemp,vv,n)

    implicit none
    external ran1sub

    REAL(8) vv(n),x
    INTEGER idumtemp,kk,n

    vvloop: do kk = 1,n
        call ran1sub(idumtemp,x)
        vv(kk) = x
    end do vvloop

END SUBROUTINE ran1



!------------------------


SUBROUTINE GASDEV(idum,vv,n)

    ! This function returns a normally distributed deviate with
    ! zero mean and unit variance, using ran1(IDUM) as the
    ! source of uniform deviates

    integer iset,n,k,gotoc
    real(8) v1,v2,r,fac,gset,gasdev1
    real(8) temp1(2),vv(n)

    external ran1

    write(*,*) 'want',n,'random numbers'

    vvloop: do k = 1,n

        iset = 0
        v1 = 0.
        v2 = 0.
        r = 0.
        fac = 0.
        gset = 0.
        gasdev1 = 0.
        ix1=0.
        ix2=0.
        gotoc = 0

        if (iset.eq.0) then


1           call ran1(idum,temp1,2)

            v1 = 2.*temp1(1)-1.
            v2 = 2.*temp1(2)-1.

            !          write(*,*) 'v1',v1
            !          write(*,*) 'v2',v2
            !          write(*,*) 'idum',idum

            r = v1**2.+v2**2.
            if (r.ge.1.) then
                gotoc = gotoc + 1
                if (gotoc.gt.n) then
                    write(*,*) 'error in gasdev'
                end if
                go to 1
            end if

            fac = dsqrt(-2.*dlog(r)/r)
            gset = v1*fac
            gasdev1 = v2*fac
            iset = 1
        else
            gasdev1 = gset
            iset = 0
        end if

        vv(k)=gasdev1
    !       write(*,*) 'random normal number',k,vv(k)

    end do vvloop

    return
end SUBROUTINE GASDEV


!-------------------------------------------------------------------------------
!       This program implements a dynamic female labor force participation
!       choice model. There are twenty periods (each representing two years) and
!       women choose in each period whether to work or not and experience
!       accumlates and affects their future wage offers.



program laborchoice

    implicit none

    external ran1,gasdev

    !      Model Parameters

    integer numper,numvar,i,s,k,st1,st2,idnumeps,idnumeta,j,t,r,o,q
    integer calcsim


    double precision bbeta,vareps,vareta,ppi,ggamma0,ggamma1,ggamma2,ddelta, aalpha
    double precision husbinc, edu, numchild, kkappa, h

    !       Store shocks for emax calculation at each state point
    !       in each decision period (shocks are used for numerical
    !       integration)   (need 20*100*20 = 40000 shocks)

    double precision alleps(40000),epsvec(100), alleta(40000), etavec(100)

    double precision utilwork(100),utilleisure(100)
    double precision maxu(100),emaxu(0:19,1:20),sum1, sum2(20),sum3

    !       Simulation Parameters
    !       store shocks for simulation of behavior of
    !       1000 women in 20 time periods (need 1000*20 periods=20000)

    double precision simshockeps(20000), simshocketa(20000)
    double precision sutwork,sutleisure,simeps(20),simeta(20)
    integer simchoice(1000,20),nn(1000,20)
    double precision laborsupply(20), laborexp
    double precision subsidy(7)

    ! Set subsidy
    subsidy(1) = 0.0
    subsidy(2) = 50.0
    subsidy(3) = 100.0
    subsidy(4) = 150.0
    subsidy(5) = 200.0
    subsidy(6) = 250.0
    subsidy(7) = 300.0

    !       Draw set of shocks for each time period for
    !       crude numerical integration (common random numbers)
    idnumeps = 23408   ! random number seed
    idnumeta = 12872   ! random number seed
    call gasdev(idnumeps,alleps,40000)
    call gasdev(idnumeta,alleta,40000)
    write(*,*) ''
    
    !      Draw simulation shocks
    call gasdev(idnumeps,simshockeps,20000)    ! draw leisure shocks
    call gasdev(idnumeta,simshocketa,20000)    ! draw wage shocks
    

    subsidyloop: do o=1,7

        write(*,*) ''
        write(*,*) ''

        ! Initialize parameters

        husbinc  = 30000.0
        numchild = 2.0
        ppi      = 300.0-subsidy(o)
        kkappa   = 22500.0
        bbeta   = 1.0
        aalpha  = 0.07
        ggamma0 = 25000.0
        ggamma1 = 0.002
        ggamma2 = -0.001
        ddelta  = 0.99
        edu     = 14.0

        alleps   = 0.0
        alleta   = 0.0
        vareps = (5000.0)**2
        vareta = (5000.0)**2

        utilwork = 0.0
        utilleisure = 0.0






        !       Solve the model backwards
        !       Will evaluate at all possible points
        !       in the observable state space

        oloop:  do i = 20,1,-1

            !          loop over all possible state points (years of experience)

            sloop:  do s = 0,i

                st1 = 2000*(20-i)+(s)*100+1
                st2 = st1+99

                epsvec = sqrt(vareps)*alleps(st1:st2)
                etavec = sqrt(vareta)*alleta(st1:st2)

                !       convert experience to real number

                h = s*1.0

                !       utility is linear
                !       income if work (exponentiate wage equation components to use standard coefficients)
                utilwork = husbinc + ggamma0 + exp(aalpha*edu + ggamma1*h + ggamma2*h*h) - ppi*numchild + etavec
                !       income if no work
                utilleisure = husbinc + kkappa*bbeta + epsvec




                !       if not the terminal period, need to add in the emax
                !       values. If terminal period, only get current period
                !       utility. In either case, need to store the emax values

                if (i.lt.20) then
                    utilleisure = utilleisure + emaxu(s,(i+1))
                    if (s.le.18) then
                        utilwork = utilwork + emaxu((s+1),(i+1))
                    else
                        utilwork = utilleisure    ! if already have maximum experience
                    end if                      ! cannot have any more

                    maxu = max(utilleisure,utilwork)  ! calculate which choice gives
                                                       ! the maximum utility
                end if

                if (i.eq.20) then          ! if in terminal period, do not
                    if (s.eq.19) then        ! need to add in future emax values
                        maxu = utilleisure
                    else
                        maxu = max(utilleisure,utilwork)
                    end if
                end if



                !            calculate and store the emax values associated with each
                !            state point  - the emax is determined by the mean of
                !            the utilities, taken over the realized shocks (which
                !            is numerically integrating with respect to the density of
                !            the shocks)

                sum1 = 0.0
                kloop: do k = 1,100    ! loop over the shocks to sum
                    sum1 = sum1+maxu(k)
                end do kloop
                emaxu(s,i) = sum1/100.0     ! calculate emax value taking avg


            end do sloop

        end do oloop


        write(*,*) 'Modelled is Solved. Moving onto simulation...'
        write(*,*) ''
    
        calcsim = 1
        if (calcsim.eq.1) then

            simchoice = 0.0

            simloop: do i =1,1000
                st1 = (i-1)*20+1
                st2 = st1+20-1
                simeps=sqrt(vareps)*simshockeps(st1:st2)
                simeta=sqrt(vareta)*simshocketa(st1:st2)

                yrloop: do j = 1,20

                    !       calculate yrs of experience based on choices

                    if (j.eq.1) then
                        h = 0.0
                    else
                        h = sum(simchoice(i,1:(j-1)))
                    end if

                    !       calculate the utility for the current period shock

                    sutwork = husbinc + ggamma0 + exp(aalpha*edu + ggamma1*h + ggamma2*h*h) - ppi*numchild + simeta(j)
                    sutleisure= husbinc + kkappa*bbeta + simeps(j)

                    if (j.lt.20) then
                        sutleisure = sutleisure + ddelta*emaxu(s,j+1)
                    end if

                    simchoice(i,j) = 0

                    if (j.lt.20) then   ! add in emax if not in terminal period
                        sutwork = sutwork + ddelta*emaxu((s+1),j+1)
                    end if
                    if (sutwork.gt.sutleisure) then
                        simchoice(i,j)=  1
                    end if

                end do yrloop
            end do simloop
        end if ! calcsim


        !       write out the simulated choices for 5 of the women, periods 10-20

        !    rloop: do i=1,5
        !        write(*,*) ''
        !        write(*,*) 'Agent No.',i
        !        write(*,*) 'Number of Periods Worked:', sum(simchoice(i,:))
        !        write(*,*) 'Choices for Periods 10-20'
        !        write(*,*) (simchoice(i,j),j=10,20)
        !        write(*,*) ''

        !    end do rloop

        write(*,*) ''
        write(*,'(A,f8.2)') 'The subsidy to childcare is ', subsidy(o)
        write(*,'(A,f8.2)') 'The net cost of childcare is ', ppi
        write(*,*) ''

        !       calculate summary statistics - percentage of women working


        sum2 = 0.0
        ssloop: do r = 1,1000    ! loop over to sum
            sum2 = sum2+simchoice(r,:)
        end do ssloop

        write(*,*) 'Labor Participation Ratio'
        ggloop: do t = 1,20
            laborsupply(t) = sum2(t)/1000.0
            write(*,15) t, laborsupply(t)
        end do ggloop

        write(*,*) ''
        sum3 = 0.0
        qqloop: do q = 1,20    ! loop over to sum
            sum3 = sum3+laborsupply(q)
        end do qqloop
    
        sum3 = sum3/20.0
    
        write(*,'(A,f8.2)') 'The avg female participation ratio over 20 years is ', sum3

        write(*,*) ''

        rrloop: do i = 1,1000
            zzloop: do j = 1,20
                nn(i,j) = sum(simchoice(i,1:j))
            end do zzloop
        end do rrloop

        write(*,*) 'Average Level of Female Experience'
        hhloop: do i = 1,20
            laborexp = sum(nn(:,i))/1000.00
            write(*,15) i, laborexp
        end do hhloop


15      format(i5,f8.2)

    end do subsidyloop

end program laborchoice
