function [ kVec_x,kVec_y,eigenvectors_of_kBloch_panel,eigenvalues_of_kBloch_panel ] = Brillouin_paths( non_Bravais_lattice, triangular_lattice, no_of_atoms_per_cell, linear_path, PATH_to_evaluate_along_kx,n1Lim,n2Lim,k0_list,kgrid,a_sp,a_ho,Gamma0_list,A_sc,zeeman_X,zeeman_Y,zeeman_Z,Free_shifts_source_I,Hermitian )

%Brillouin_paths This function calculates the complex determinant dips
%along a specified BZ path.

    k0_A = k0_list(1);
    k0_B = k0_list(2);
    k0_C = k0_list(3);
    Gamma0_A =Gamma0_list(1);
    Gamma0_B =Gamma0_list(2);
    Gamma0_C =Gamma0_list(3);

    Omega_grid = k0_A.*ones(kgrid,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% K_BLOCH DEPENDENT QUANTITIES  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if triangular_lattice == 1
        
        [ kVec_x, kVec_y ] = three_fold_symmetric_BZ( linear_path, a_sp, PATH_to_evaluate_along_kx, kgrid );
    
    else
       
        [ kVec_x, kVec_y ] = two_fold_symmetric_BZ( linear_path, a_sp, PATH_to_evaluate_along_kx, kgrid );
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %initially get diagonal term:
    non_Bravais_phase=0;

    [ summed_free_XX_diag,summed_free_YY_diag,~,summed_free_XY_and_YX_diag ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,mat2cell([NaN,NaN],1,2));
    
    
    
    [ matrix_implicit_xx_AA,matrix_implicit_yy_AA,matrix_implicit_xy_AA,matrix_implicit_yx_AA ] = diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_diag,summed_free_YY_diag,summed_free_XY_and_YX_diag,k0_A,a_ho,Gamma0_A,zeeman_Z );

    
    if non_Bravais_lattice == 0
        
        [ eigenvectors_of_kBloch_panel, eigenvalues_of_kBloch_panel ] = matrix_B_2x2( matrix_implicit_xx_AA, matrix_implicit_xy_AA, matrix_implicit_yx_AA, matrix_implicit_yy_AA, kgrid, no_of_atoms_per_cell, Hermitian );
         
    else
        
        %diagonal Matrix_BB term
        [ matrix_implicit_xx_BB,matrix_implicit_yy_BB,matrix_implicit_xy_BB,matrix_implicit_yx_BB ] = diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_diag,summed_free_YY_diag,summed_free_XY_and_YX_diag,k0_B,a_ho,Gamma0_B,zeeman_Z );

        %info about inter-cell basis vector:
        b_vector = mat2cell([1,2],1,2);
        
        %off diagonal Matrix_AB term:
        non_Bravais_phase=1;

        [ summed_free_XX_AB,summed_free_YY_AB,~,summed_free_XY_and_YX_AB ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
        [ matrix_implicit_xx_AB,matrix_implicit_yy_AB,matrix_implicit_xy_AB,matrix_implicit_yx_AB ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_AB,summed_free_YY_AB,summed_free_XY_and_YX_AB,k0_A,k0_B,a_ho,Gamma0_A,Gamma0_B,zeeman_Z );


        %off diagonal Matrix_BA term:
        non_Bravais_phase=-1;

        [ summed_free_XX_BA,summed_free_YY_BA,~,summed_free_XY_and_YX_BA ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
        [ matrix_implicit_xx_BA,matrix_implicit_yy_BA,matrix_implicit_xy_BA,matrix_implicit_yx_BA ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_BA,summed_free_YY_BA,summed_free_XY_and_YX_BA,k0_A,k0_B,a_ho,Gamma0_A,Gamma0_B,zeeman_Z );

        
        if non_Bravais_lattice == 1
            
            [ eigenvectors_of_kBloch_panel, eigenvalues_of_kBloch_panel ] = matrix_NB_4x4(  matrix_implicit_xx_AA , matrix_implicit_xy_AA ,      matrix_implicit_xx_AB , matrix_implicit_xy_AB , ...
                                                                                            matrix_implicit_yx_AA , matrix_implicit_yy_AA ,      matrix_implicit_yx_AB , matrix_implicit_yy_AB , ...
                                                                                                                                                                                 ...
                                                                                            matrix_implicit_xx_BA , matrix_implicit_xy_BA ,      matrix_implicit_xx_BB , matrix_implicit_xy_BB , ...
                                                                                            matrix_implicit_yx_BA , matrix_implicit_yy_BA ,      matrix_implicit_yx_BB , matrix_implicit_yy_BB, kgrid, no_of_atoms_per_cell, Hermitian );

        elseif non_Bravais_lattice == 2
            

            %diagonal Matrix_CC term
            [ matrix_implicit_xx_CC,matrix_implicit_yy_CC,matrix_implicit_xy_CC,matrix_implicit_yx_CC ] = diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_diag,summed_free_YY_diag,summed_free_XY_and_YX_diag,k0_C,a_ho,Gamma0_C,zeeman_Z );

            
            %info about inter-cell basis vector:
            b_vector = mat2cell([1,3],1,2);
            
            %off diagonal Matrix_AC term:
            non_Bravais_phase=1;

            [ summed_free_XX_AC,summed_free_YY_AC,~,summed_free_XY_and_YX_AC ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
            [ matrix_implicit_xx_AC,matrix_implicit_yy_AC,matrix_implicit_xy_AC,matrix_implicit_yx_AC ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_AC,summed_free_YY_AC,summed_free_XY_and_YX_AC,k0_A,k0_C,a_ho,Gamma0_A,Gamma0_C,zeeman_Z );


            %off diagonal Matrix_CA term:
            non_Bravais_phase=-1;

            [ summed_free_XX_CA,summed_free_YY_CA,~,summed_free_XY_and_YX_CA ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
            [ matrix_implicit_xx_CA,matrix_implicit_yy_CA,matrix_implicit_xy_CA,matrix_implicit_yx_CA ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_CA,summed_free_YY_CA,summed_free_XY_and_YX_CA,k0_A,k0_C,a_ho,Gamma0_A,Gamma0_C,zeeman_Z );

            
            b_vector = mat2cell([2,3],1,2);
            
            %off diagonal Matrix_BC term:
            non_Bravais_phase=1;

            [ summed_free_XX_BC,summed_free_YY_BC,~,summed_free_XY_and_YX_BC ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
            [ matrix_implicit_xx_BC,matrix_implicit_yy_BC,matrix_implicit_xy_BC,matrix_implicit_yx_BC ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_BC,summed_free_YY_BC,summed_free_XY_and_YX_BC,k0_B,k0_C,a_ho,Gamma0_B,Gamma0_C,zeeman_Z );


            %off diagonal Matrix_CA term:
            non_Bravais_phase=-1;

            [ summed_free_XX_CB,summed_free_YY_CB,~,summed_free_XY_and_YX_CB ] = Summed_free_terms(non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector);
            [ matrix_implicit_xx_CB,matrix_implicit_yy_CB,matrix_implicit_xy_CB,matrix_implicit_yx_CB ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_CB,summed_free_YY_CB,summed_free_XY_and_YX_CB,k0_B,k0_C,a_ho,Gamma0_B,Gamma0_C,zeeman_Z );

            
            [ eigenvectors_of_kBloch_panel, eigenvalues_of_kBloch_panel ] = matrix_NB_6x6(  matrix_implicit_xx_AA , matrix_implicit_xy_AA ,      matrix_implicit_xx_AB , matrix_implicit_xy_AB , matrix_implicit_xx_AC , matrix_implicit_xy_AC, ...
                                                                                            matrix_implicit_yx_AA , matrix_implicit_yy_AA ,      matrix_implicit_yx_AB , matrix_implicit_yy_AB , matrix_implicit_yx_AC , matrix_implicit_yy_AC, ...
                                                                                                                                                                                 ...
                                                                                            matrix_implicit_xx_BA , matrix_implicit_xy_BA ,      matrix_implicit_xx_BB , matrix_implicit_xy_BB , matrix_implicit_xx_BC , matrix_implicit_xy_BC, ...
                                                                                            matrix_implicit_yx_BA , matrix_implicit_yy_BA ,      matrix_implicit_yx_BB , matrix_implicit_yy_BB , matrix_implicit_yx_BC , matrix_implicit_yy_BC, ... 
                                                                                                                                                                                 ...
                                                                                            matrix_implicit_xx_CA , matrix_implicit_xy_CA ,      matrix_implicit_xx_CB , matrix_implicit_xy_CB , matrix_implicit_xx_CC , matrix_implicit_xy_CC, ...
                                                                                            matrix_implicit_yx_CA , matrix_implicit_yy_CA ,      matrix_implicit_yx_CB , matrix_implicit_yy_CB , matrix_implicit_yx_CC , matrix_implicit_yy_CC, kgrid, no_of_atoms_per_cell, Hermitian );

            
            
        end
        
    end
    
 
    
end
    
    
