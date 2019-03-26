function [ matrix_implicit_xx_off,matrix_implicit_yy_off,matrix_implicit_xy_off,matrix_implicit_yx_off ] = off_diagonal_component( Free_shifts_source_I,A_sc,Omega_grid,summed_free_XX_off,summed_free_YY_off,summed_free_XY_and_YX_off,k0_1,k0_2,a_ho,Gamma0_1,Gamma0_2,zeeman_Z ) 


    free_shifts_a_ho_xx_off= + (1/A_sc).*Omega_grid.^2.*summed_free_XX_off;
    free_shifts_a_ho_yy_off= + (1/A_sc).*Omega_grid.^2.*summed_free_YY_off;
    free_shifts_a_ho_xy_off= + (1/A_sc).*Omega_grid.^2.*summed_free_XY_and_YX_off;
    free_shifts_a_ho_yx_off= + (1/A_sc).*Omega_grid.^2.*summed_free_XY_and_YX_off;

    
    free_shifts_renormalized_xx_off=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_xx_off;
    free_shifts_renormalized_yy_off=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_yy_off;
    free_shifts_renormalized_xy_off=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_xy_off;
    free_shifts_renormalized_yx_off=exp( (Omega_grid.*a_ho).^2  ).*free_shifts_a_ho_yx_off;


    %prefactor for off-diagonal matrix elements proportional to d_A & d_B:
    d_1_d_2_per_eps0_prefactor = sqrt(  (Gamma0_1*(3*pi)/(k0_1^3)) .* (Gamma0_2*(3*pi)/(k0_2^3))  );
    
    matrix_implicit_xx_off= d_1_d_2_per_eps0_prefactor .* (free_shifts_renormalized_xx_off    );
    matrix_implicit_yy_off= d_1_d_2_per_eps0_prefactor .* (free_shifts_renormalized_yy_off    );
    matrix_implicit_xy_off= d_1_d_2_per_eps0_prefactor .* (free_shifts_renormalized_xy_off    );
    matrix_implicit_yx_off= d_1_d_2_per_eps0_prefactor .* (free_shifts_renormalized_yx_off    );
    

end

