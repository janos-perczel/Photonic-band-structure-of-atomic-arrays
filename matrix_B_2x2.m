function [ eigenvectors_of_kBloch_panel, eigenvalue_of_kBloch_panel ] = matrix_B_2x2( matrix_implicit_xx, matrix_implicit_xy, matrix_implicit_yx, matrix_implicit_yy, kgrid, no_of_atoms_per_cell,Hermitian )

    

    eigenvectors_of_kBloch_panel = zeros(kgrid,2*no_of_atoms_per_cell,2*no_of_atoms_per_cell);
    eigenvalue_of_kBloch_panel = zeros(kgrid,2*no_of_atoms_per_cell);
    
    for k_Bloch_index = 1:kgrid

        matrix_element_xx = matrix_implicit_xx(k_Bloch_index);
        matrix_element_xy = matrix_implicit_xy(k_Bloch_index);
        matrix_element_yx = matrix_implicit_yx(k_Bloch_index);
        matrix_element_yy = matrix_implicit_yy(k_Bloch_index);

        matrix_for_eigenvectors = [ matrix_element_xx , matrix_element_xy ;
                                    matrix_element_yx , matrix_element_yy ];

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

