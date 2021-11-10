c======================================================================c

      subroutine start_fam( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE fam_energies;
      USE f02f20matrix;
      USE h02h20matrix;
      USE xyfam;
      USE fam_iter;
      USE dh;
      USE dDelta;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      ! Chebyshev RPA
      integer Nit;
      integer tape_mu;
      double precision Omegab;
      double complex ztmp;
      double precision, dimension(:),     allocatable :: mu;
      double complex,   dimension(:,:,:), allocatable :: alpha_new_x;
      double complex,   dimension(:,:,:), allocatable :: alpha_new_y;
      double complex,   dimension(:,:,:), allocatable :: alpha_old_x;
      double complex,   dimension(:,:,:), allocatable :: alpha_old_y;
      double complex,   dimension(:,:,:), allocatable :: alpha_tmp_x;
      double complex,   dimension(:,:,:), allocatable :: alpha_tmp_y;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN start_fam() ***************************';
      write(6,*) '';
      endif






#ifdef DEBUG
      tol = 1.D-9;  ! Self-consistency tolerance
#else
      tol = 1.D-7;  ! Self-consistency tolerance
#endif






c-----Printing strength.out header
      open( tape_strength , file   = './output/QFAM_output/strength.out'
     &                    , status = 'unknown' );

      write( tape_strength , '(a)' , advance = 'no' )
     &      'S(f,omega) = -1/pi * Im[Tr[hermconj(f)*drho(omega)]] ';

      ! Unit of measurement of S(f,omega)
      select case( J_multipole )
          case( 0 )
              write(tape_strength,*) '[fm^4/MeV]';
          case( 1 )
              select case( ISO )
                  case( 0 )
                      write(tape_strength,*) '[fm^2/MeV]';
                  case( 1 )
                      write(tape_strength,*) '[e^2fm^2/MeV]';
                  case default
                      stop 'Error: Wrong ISO!';
              end select
          case( 2 )
              write(tape_strength,*) '[fm^4/MeV]';
          case( 3 )
              write(tape_strength,*) '[fm^6/MeV]';
          case default
              stop 'Error: J_multipole > 3 not implemented!';
      end select

      call print_header( tape_strength );

      write(tape_strength,'(/,a,21x,a,/)')'omega[MeV/hbar]',
     &                                    'S(f,omega)';

      if( i_calculation_type .eq. 0 ) then
          write(6,'(/,17x,a,/)') 'Free response';
          write(6,'(a,21x,a,/)') 'omega[MeV/hbar]' , 'S(f,omega)';
      endif

      call flush(tape_strength);
      call flush(6);






c-----Free response, energy sweep
      if( i_calculation_type .eq. 0 ) then

          ! Free response is defined as a response
          ! which does not take into account induced
          ! self-consistent Hamiltonian (H02 = H20 = 0).

          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call fam_xy      ( .false. );
              call fam_spurious( .false. );
              call fam_drho    ( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(            6,'(1f15.5,1f31.10)') omega , Sp+Sn;
              write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(            6);
              call flush(tape_strength);

              omega = omega + delta_omega;
          enddo

      endif






c-----Fully self-consistent response, energy sweep
      if( i_calculation_type .eq. 1 ) then

          ! Initial guess for self-consistent solution
          ! vector for initial sweep energy is zero.
          ! Otherwise, we use self-consistent solution
          ! of previous energy as initial guess
          ! for the following energy.
          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          ! One can easily parallelize the following energy
          ! loop for systematic, large scale calculations
          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call iter_fam( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(tape_strength);

              omega = omega + delta_omega;
          enddo

      endif






c-----Calculation for given energy, printing vector density
      if( i_calculation_type .eq. 2 ) then

          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          omega = omega_print;
          call iter_fam( .false. );

          Sn = fam_strength( .false. , 1 );
          Sp = fam_strength( .false. , 2 );

          write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
          call flush(tape_strength);

          call print_dens( .false. );

          ! Calculates and prints the localization function
          !call       nuclocfunc( .false. );
          !call print_nuclocfunc( .false. );

      endif






c-----Chebyshev RPA method
      if( i_calculation_type .eq. 3 ) then

          Nit    = 100000;
          Omegab = +4500.d0;






          open( newunit = tape_mu ,
     &          file    = './output/QFAM_output/mu.out' ,
     &          status  = 'unknown'                       );

          call print_header( tape_mu );

          write(tape_mu,'(a,f9.3,a)')'Omegab = ',Omegab,' [MeV]';
          write(      6,'(a,f9.3,a)')'Omegab = ',Omegab,' [MeV]';

          write( tape_mu , '(/,11x,a,10x,a,/)') 'n' , 'mu(n)';
          write(       6 , '(/,11x,a,10x,a,/)') 'n' , 'mu(n)';

          call flush( tape_mu );
          call flush(       6 );






          allocate( mu( 0 : 2*Nit ) );
          allocate( alpha_new_x( N_total , N_total , 2 ) );
          allocate( alpha_new_y( N_total , N_total , 2 ) );
          allocate( alpha_old_x( N_total , N_total , 2 ) );
          allocate( alpha_old_y( N_total , N_total , 2 ) );
          allocate( alpha_tmp_x( N_total , N_total , 2 ) );
          allocate( alpha_tmp_y( N_total , N_total , 2 ) );






          ! Initializing alpha_old
          do it = 1 , 2
              do j = 1 , N_total
                  do i = 1 , N_total
                      alpha_old_x(i,j,it) = + (-f20(i,j,it));
                      alpha_old_y(i,j,it) = - (-f02(j,i,it));
                  enddo
              enddo
          enddo

          ! Initializing alpha_new
          call applyMapping_XS( N_total     , Omegab      ,
     &                          alpha_old_x , alpha_old_y ,
     &                          alpha_new_x , alpha_new_y   );






          ! Calculating mu(0)
          ztmp = complex(0.d0,0.d0);
          do it = 1 , 2
              do j = 1 , N_total
                  do i = 1 , N_total

                      ztmp = ztmp + conjg(-f20(i,j,it))
     &                            * alpha_old_x(i,j,it);

                      ztmp = ztmp + conjg(-f02(j,i,it))
     &                            * alpha_old_y(i,j,it);

                  enddo
              enddo
          enddo
          mu(0) = DREAL( 0.5d0 * ztmp );






          ! Calculating mu(1)
          ztmp = complex(0.d0,0.d0);
          do it = 1 , 2
              do j = 1 , N_total
                  do i = 1 , N_total

                      ztmp = ztmp + conjg(-f20(i,j,it))
     &                            * alpha_new_x(i,j,it);

                      ztmp = ztmp + conjg(-f02(j,i,it))
     &                            * alpha_new_y(i,j,it);

                  enddo
              enddo
          enddo
          mu(1) = DREAL( ztmp );








          write( tape_mu , * ) 0 , mu(0);
          write(       6 , * ) 0 , mu(0);

          do n = 1 , Nit



              ! Calculating mu(2*n-1)
              ztmp = complex(0.d0,0.d0);
              do it = 1 , 2
                  do j = 1 , N_total
                      do i = 1 , N_total

                          ztmp = ztmp + conjg(alpha_old_x(i,j,it))
     &                                *       alpha_new_x(i,j,it);

                          ztmp = ztmp - conjg(alpha_old_y(i,j,it))
     &                                *       alpha_new_y(i,j,it);

                      enddo
                  enddo
              enddo
              mu(2*n-1) = 2.d0*dreal(ztmp) - mu(1);

              ! Calculating mu(2*n)
              ztmp = complex(0.d0,0.d0);
              do it = 1 , 2
                  do j = 1 , N_total
                      do i = 1 , N_total

                          ztmp = ztmp + conjg(alpha_new_x(i,j,it))
     &                                *       alpha_new_x(i,j,it);

                          ztmp = ztmp - conjg(alpha_new_y(i,j,it))
     &                                *       alpha_new_y(i,j,it);

                      enddo
                  enddo
              enddo
              mu(2*n) = 2.d0*dreal(ztmp) - 2.d0*mu(0);



              write( tape_mu , * ) 2*n-1 , mu(2*n-1);
              write(       6 , * ) 2*n-1 , mu(2*n-1);

              write( tape_mu , * ) 2*n   , mu(2*n);
              write(       6 , * ) 2*n   , mu(2*n);



              ! Updating alpha_old and alpha_new
              if( n .eq. Nit ) cycle;

              ! alpha_tmp = alpha_new
              do it = 1 , 2
                  do j = 1 , N_total
                      do i = 1 , N_total
                          alpha_tmp_x(i,j,it) = alpha_new_x(i,j,it);
                          alpha_tmp_y(i,j,it) = alpha_new_y(i,j,it);
                      enddo
                  enddo
              enddo

              ! alpha_new = 2/Omegab * X * S * alpha_new - alpha_old
              call applyMapping_XS( N_total     , Omegab      ,
     &                              alpha_tmp_x , alpha_tmp_y ,
     &                              alpha_new_x , alpha_new_y   );
              do it = 1 , 2
                  do j = 1 , N_total
                      do i = 1 , N_total

                      alpha_new_x(i,j,it) = 2.d0 * alpha_new_x(i,j,it)
     &                                           - alpha_old_x(i,j,it);

                      alpha_new_y(i,j,it) = 2.d0 * alpha_new_y(i,j,it)
     &                                           - alpha_old_y(i,j,it);

                      enddo
                  enddo
              enddo


              ! alpha_old = alpha_tmp
              do it = 1 , 2
                  do j = 1 , N_total
                      do i = 1 , N_total
                          alpha_old_x(i,j,it) = alpha_tmp_x(i,j,it);
                          alpha_old_y(i,j,it) = alpha_tmp_y(i,j,it);
                      enddo
                  enddo
              enddo



          enddo











          deallocate( mu );
          deallocate( alpha_new_x );
          deallocate( alpha_new_y );
          deallocate( alpha_old_x );
          deallocate( alpha_old_y );
          deallocate( alpha_tmp_x );
          deallocate( alpha_tmp_y );



          close(tape_mu);

      endif






      close(tape_strength);






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END start_fam() *****************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine applyMapping_XS( Nxy , Omegab , xin,yin , xout,yout )

c======================================================================c
      ! Calculates [ xout , yout ] = 1/Omegab * X * S * [ xin ; yin ],
      ! where X = [ I , 0 ; 0 , -I ] and S = [ A , B ; B^* , A^* ].

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE fam_energies;
      USE f02f20matrix;
      USE h02h20matrix;
      USE xyfam;
      USE dh;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)


      character fg1, fg2;
      integer,          intent(in)  :: Nxy;
      double precision, intent(in)  :: Omegab;
      double complex,   intent(in)  ::  xin( Nxy , Nxy , 2 );
      double complex,   intent(in)  ::  yin( Nxy , Nxy , 2 );
      double complex,   intent(out) :: xout( Nxy , Nxy , 2 );
      double complex,   intent(out) :: yout( Nxy , Nxy , 2 );



      x_fam = complex(0.d0,0.d0);
      y_fam = complex(0.d0,0.d0);
      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total

                  x_fam(i,j,it) = - xin(i,j,it);
                  y_fam(i,j,it) = - yin(j,i,it);

              enddo
          enddo
      enddo

      call fam_drho       ( .false. );

      call fam_dkappa     ( .false. );

      call fam_ddensdcurr ( .false. );

      call fam_dpotentials( .false. );

      call fam_dh1        ( .false. );

c-----Construction of matrix dh_2
      dh_2 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
         do ib2 = 1 , N_blocks
           do ib1 = 1 , N_blocks

              j0 = ia_spx(ib2);
              j1 = j0+id_spx(ib2)-1;
              do j = j0 , j1
                 fg2 = fg_spx(j-j0+1,ib2);
                 ml2 = ml_spx(j-j0+1,ib2);

                 i0 = ia_spx(ib1);
                 i1 = i0+id_spx(ib1)-1;
                 do i = i0 , i1
                    fg1 = fg_spx(i-i0+1,ib1);
                    ml1 = ml_spx(i-i0+1,ib1);

                    if( abs(ml1-ml2) .eq. K_multipole ) then
                        if( fg1 .eq. fg2 ) then
                            dh_2(i,j,it) = + dh_1(i,j,it);
                        else
                            dh_2(i,j,it) = - dh_1(j,i,it);
                        endif
                    elseif( abs(ml1+ml2+1) .eq. K_multipole ) then
                        if( fg1 .ne. fg2 ) then
                            dh_2(i,j,it) = - dh_1(j,i,it);
                        endif
                    endif

              enddo

            enddo

          enddo
        enddo
      enddo

      call fam_ddelta     ( .false. );

      call fam_h20h02     ( .false. );

      do it = 1 , 2
          do j = 1 , N_total
              Ej = E_fam(j,it);
              do i = 1 , N_total
                  Ei = E_fam(i,it);

                  xout(i,j,it) = - h20(i,j,it) + (Ei+Ej)*xin(i,j,it);
                  yout(i,j,it) = - h02(j,i,it) + (Ei+Ej)*yin(i,j,it);

              enddo
          enddo
      enddo

      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total

                  xout(i,j,it) = + 1.d0/Omegab * xout(i,j,it);
                  yout(i,j,it) = - 1.d0/Omegab * yout(i,j,it);

              enddo
          enddo
      enddo



      return;
      end;
