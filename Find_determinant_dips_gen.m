function [ kVec_x_ALL,kVec_y_ALL,eigenvectors_ALL,eigenvalues_ALL] = Find_determinant_dips_gen(non_Bravais_lattice, triangular_lattice, no_of_atoms_per_cell, linear_path,kgrid,a_sp,a_ho,n1Lim,n2Lim,k0_list,Gamma0_list,zeeman_X,zeeman_Y,zeeman_Z,Hermitian)
%Find_determinant_dips --- This function returns the 1000x100 matrix of
%determinantal values as a function of input parameters, minimized over the
%given Gamma range.

    k0_A=k0_list(1);
    
    %area of unit cell
    if triangular_lattice == 1
        A_sc=sqrt(3)/2*a_sp.^2;
    else
        A_sc = a_sp.^2;
    end

%     k0_A
    [Free_shifts_source_I] = k0_A.^3./(6*pi).*(  ( (-1/2)+ 2 .*(k0_A.*a_ho).^2  )./( 2.*pi^(1/2).*(k0_A.*a_ho).^3 )-  (Faddeeva_erfi(k0_A.*a_ho))./(exp( (k0_A.*a_ho).^2 ))        );
%     Free_shifts_source_I
    
    disp('Greens function at source evaluated.')
    
    
    if linear_path
        
        BZ_path_grid = 3;
        
    else %2D BZ sampling.
        
        BZ_path_grid = kgrid;
        
    end

    kVec_x_ALL =            zeros(BZ_path_grid,kgrid);
    kVec_y_ALL =            zeros(BZ_path_grid,kgrid);

    eigenvectors_ALL =      zeros(BZ_path_grid,kgrid,2*no_of_atoms_per_cell,2*no_of_atoms_per_cell);
    eigenvalues_ALL  =      zeros(BZ_path_grid,kgrid,2*no_of_atoms_per_cell);

    for path_index_along_kx=1:BZ_path_grid %will be 1 if linear problem 
            
            [ kVec_x,kVec_y,eigenvectors_of_ky_panel,eigenvalues_of_ky_panel] = Brillouin_paths (non_Bravais_lattice, triangular_lattice, no_of_atoms_per_cell, linear_path, path_index_along_kx,n1Lim,n2Lim,k0_list,kgrid,a_sp,a_ho,Gamma0_list,A_sc,zeeman_X,zeeman_Y,zeeman_Z,Free_shifts_source_I,Hermitian );

            kVec_x_ALL(path_index_along_kx,:) = kVec_x;
            kVec_y_ALL(path_index_along_kx,:) = kVec_y;
            

            eigenvectors_ALL(path_index_along_kx,:,:,:) = eigenvectors_of_ky_panel;
            eigenvalues_ALL(path_index_along_kx,:,:) = eigenvalues_of_ky_panel;

            disp(['Path ',num2str(path_index_along_kx),' evaluated.'])
            
    end


end

