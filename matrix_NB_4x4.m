function [ eigenvectors_of_kBloch_panel, eigenvalue_of_kBloch_panel ] = matrix_NB_4x4(  matrix_implicit_xx_AA , matrix_implicit_xy_AA ,      matrix_implicit_xx_AB , matrix_implicit_xy_AB , ...
                                                                                        matrix_implicit_yx_AA , matrix_implicit_yy_AA ,      matrix_implicit_yx_AB , matrix_implicit_yy_AB , ...
                                                                                                                                                                             ...
                                                                                        matrix_implicit_xx_BA , matrix_implicit_xy_BA ,      matrix_implicit_xx_BB , matrix_implicit_xy_BB , ...
                                                                                        matrix_implicit_yx_BA , matrix_implicit_yy_BA ,      matrix_implicit_yx_BB , matrix_implicit_yy_BB, kgrid, no_of_atoms_per_cell,Hermitian )

    eigenvectors_of_kBloch_panel = zeros(kgrid,2*no_of_atoms_per_cell,2*no_of_atoms_per_cell);
    eigenvalue_of_kBloch_panel = zeros(kgrid,2*no_of_atoms_per_cell);
                                        
    for k_Bloch_index = 1:kgrid

            matrix_element_xx_AA = matrix_implicit_xx_AA(k_Bloch_index);
            matrix_element_xy_AA = matrix_implicit_xy_AA(k_Bloch_index);
            matrix_element_yx_AA = matrix_implicit_yx_AA(k_Bloch_index);
            matrix_element_yy_AA = matrix_implicit_yy_AA(k_Bloch_index);

            matrix_element_xx_BB = matrix_implicit_xx_BB(k_Bloch_index);
            matrix_element_xy_BB = matrix_implicit_xy_BB(k_Bloch_index);
            matrix_element_yx_BB = matrix_implicit_yx_BB(k_Bloch_index);
            matrix_element_yy_BB = matrix_implicit_yy_BB(k_Bloch_index);

            matrix_element_xx_AB = matrix_implicit_xx_AB(k_Bloch_index);
            matrix_element_xy_AB = matrix_implicit_xy_AB(k_Bloch_index);
            matrix_element_yx_AB = matrix_implicit_yx_AB(k_Bloch_index);
            matrix_element_yy_AB = matrix_implicit_yy_AB(k_Bloch_index);

            matrix_element_xx_BA = matrix_implicit_xx_BA(k_Bloch_index);
            matrix_element_xy_BA = matrix_implicit_xy_BA(k_Bloch_index);
            matrix_element_yx_BA = matrix_implicit_yx_BA(k_Bloch_index);
            matrix_element_yy_BA = matrix_implicit_yy_BA(k_Bloch_index);            

            matrix_for_eigenvectors = [ matrix_element_xx_AA , matrix_element_xy_AA ,      matrix_element_xx_AB , matrix_element_xy_AB ; ...
                                        matrix_element_yx_AA , matrix_element_yy_AA ,      matrix_element_yx_AB , matrix_element_yy_AB ; ...
                                                                                                                             ...
                                        matrix_element_xx_BA , matrix_element_xy_BA ,      matrix_element_xx_BB , matrix_element_xy_BB ; ...
                                        matrix_element_yx_BA , matrix_element_yy_BA ,      matrix_element_yx_BB , matrix_element_yy_BB ];


            % %             test if matrix contains NaN:
            find_all_NaNs = isnan(matrix_for_eigenvectors);
            sum_of_no_of_NaNs = sum(find_all_NaNs(:));    
            if sum_of_no_of_NaNs>0
                continue
            end
            
            if Hermitian
                matrix_for_eigenvectors = 1/2.*(matrix_for_eigenvectors+matrix_for_eigenvectors');
            end

            [eigenvector_matrix,eigenvalue_matrix] = eig(matrix_for_eigenvectors);

            eigenvalue_column = diag(eigenvalue_matrix);
            
            %sort by real value!
            [~,IDorig]=sort(real(eigenvalue_column));
            eigenvalue_column=eigenvalue_column(IDorig);
            eigenvector_matrix=eigenvector_matrix(:,IDorig);

            eigenvectors_of_kBloch_panel(k_Bloch_index,:,:)  = eigenvector_matrix;
            eigenvalue_of_kBloch_panel(k_Bloch_index,:)  =  eigenvalue_column;

    end
    
end

