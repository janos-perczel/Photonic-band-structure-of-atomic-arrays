%%%NB: counts bands from bottom!
plot_area_within_cone = 1;

%%%% separate area within lightcone and outside %%%%%%
lightline_Gamma_centered = @(x,y) sqrt((x).^2 + (y).^2);
within_light_cone = lightline_Gamma_centered( kVec_x_ALL(:,:),kVec_y_ALL(:,:) ) < k0_A;
outside_light_cone = lightline_Gamma_centered( kVec_x_ALL(:,:),kVec_y_ALL(:,:) ) > k0_A;

%%%% process eigenvector -- independent of plotting %%%%
eigenvector_2D_map = zeros(size(kVec_x_ALL,1),size(kVec_x_ALL,2),2*no_of_atoms_per_cell,2*no_of_atoms_per_cell,2*no_of_atoms_per_cell);

for band_no=1:2*no_of_atoms_per_cell
    for eigvector_entry = 1:2*no_of_atoms_per_cell
        
        eigvector_sliver_within = squeeze(eigenvectors_ALL(:,:,eigvector_entry,band_no));
        eigvector_sliver_outside = squeeze(eigenvectors_ALL(:,:,eigvector_entry,band_no));
        
        eigvector_sliver_within(within_light_cone == 0) = NaN;
        eigvector_sliver_outside(outside_light_cone == 0) = NaN;

        eigenvector_2D_map(:,:,band_no, 1 ,eigvector_entry) = eigvector_sliver_within(:,:);
        eigenvector_2D_map(:,:,band_no, 2 ,eigvector_entry) = eigvector_sliver_outside(:,:);

    end
end




custom_z_axis = 1;

%%%%%%%%%% 3D scatter plot of bands %%%%%%%%
if custom_z_axis
    z_min = -5;
    z_max = 5;
else
    z_max = 1.1*max((squeeze(real(eigenvalues_ALL(:)))-k0_A)./Gamma0_A);
    z_min = 1.1*min((squeeze(real(eigenvalues_ALL(:)))-k0_A)./Gamma0_A);
end

figure
hold on
for band_no=1:2*no_of_atoms_per_cell
    
    eig_values_squeezed = squeeze((real(eigenvalues_ALL(:,:,band_no))-k0_A)./Gamma0_A);
    
%     if plot_area_within_cone
%         eig_values_squeezed(within_light_cone == 0) = NaN;
%         
%     else
%         eig_values_squeezed(outside_light_cone == 0) = NaN;
%     end
    
    scatter3(kVec_x_ALL(:),kVec_y_ALL(:),eig_values_squeezed(:))
    
end
% zlim([z_min z_max])

