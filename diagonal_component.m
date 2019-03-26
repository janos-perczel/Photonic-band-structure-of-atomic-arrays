function [ matrix_implicit_xx,matrix_implicit_yy,matrix_implicit_xy,matrix_implicit_yx ] = diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_diag,summed_free_YY_diag,summed_free_XY_and_YX_diag,k0,a_ho,Gamma0,zeeman_Z )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% DIAGONAL MATRICES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    free_shifts_a_ho_xx_diag= Free_shifts_source_I + (1/A_sc).*Omega_grid.^2.*summed_free_XX_diag;
    free_shifts_a_ho_yy_diag= Free_shifts_source_I + (1/A_sc).*Omega_grid.^2.*summed_free_YY_diag;
    free_shifts_a_ho_xy_diag=                      + (1/A_sc).*Omega_grid.^2.*summed_free_XY_and_YX_diag;
    free_shifts_a_ho_yx_diag=                      + (1/A_sc).*Omega_grid.^2.*summed_free_XY_and_YX_diag;

    
    free_shifts_renormalized_xx_diag=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_xx_diag;
    free_shifts_renormalized_yy_diag=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_yy_diag;
    free_shifts_renormalized_xy_diag=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_xy_diag;
    free_shifts_renormalized_yx_diag=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_yx_diag;

    %create implicit matrix
    matrix_implicit_xx=(k0)+Gamma0*(3*pi)/(k0^3).*(free_shifts_renormalized_xx_diag);
    matrix_implicit_yy=(k0)+Gamma0*(3*pi)/(k0^3).*(free_shifts_renormalized_yy_diag);
    matrix_implicit_xy=                 +Gamma0*(3*pi)/(k0^3).*(free_shifts_renormalized_xy_diag)+1i*zeeman_Z;
    matrix_implicit_yx=                 +Gamma0*(3*pi)/(k0^3).*(free_shifts_renormalized_yx_diag)-1i*zeeman_Z;
  
    

end

