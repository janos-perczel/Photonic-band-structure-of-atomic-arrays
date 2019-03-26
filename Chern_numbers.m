%%%%This script calculates the Chern number using the Fukui, Hatsugai and Suzuki method.
%%%NB: counts bands from bottom!

light_cone_separated = 1;
% close all
band_to_calculate = 2;
%whether above or below light line
above_lightline = 1;

BZsizeY=size(kVec_y_ALL,1);
BZsizeX=size(kVec_x_ALL,2);


%%%% create DATA_Matrix containing all eVals and eVecs of discretized
%%%% BZ lattice (where q=3 is the number of bands in Harper):



if above_lightline
    lightline_band = 1;
else
    lightline_band = 2;
end

%  Lattice_DATA_eVals=omega_2D_map(:,:,band_to_calculate);
if light_cone_separated
    Lattice_DATA_eVecs = eigenvector_2D_map(:,:,band_to_calculate,lightline_band,:);
else
    Lattice_DATA_eVecs = eigenvectors_ALL(:,:,:,band_to_calculate);
end
%  
%  size(Lattice_DATA_eVecs)
%  
% squeeze(Lattice_DATA_eVecs(1,2,:))
 

 
 



% figure
% hold on
% %NB: order is special for mesh function - first entry is columns and second is rows
% mesh(kVec_x_ALL,kVec_y_ALL,omega_2D_map) 
% hold off


%%%Create 'lattice field strength' mattrix on lattice:

lattice_Field_Strength=zeros(BZsizeY-1,BZsizeX-1);



%%%Calculate 'lattice field strength' 
%%%in the middle band (Lattice_DATA_eVecs(kx,ky,:,2)):


for j=1:(BZsizeX-1)
    for l=1:(BZsizeY-1)
% for j=10:20
%     for l=10:20
        
        m=j+1; % 2<=m<=BZsizeY in plaquette counting
        n=l;   % 1<=n<=(BZsizeX-1) in plaquette counting
        
        bottomLeftState =squeeze(Lattice_DATA_eVecs(m,n,:));
        bottomRightState=squeeze(Lattice_DATA_eVecs(m,n+1,:));
        
        topLeftState =squeeze(Lattice_DATA_eVecs(m-1,n,:));
        topRightState=squeeze(Lattice_DATA_eVecs(m-1,n+1,:));
        
        %%% Calculate field strength of plaquette (m,n):
        
        bottomLink=bottomRightState'*bottomLeftState;
        Nbottom=abs(bottomLink);
        bottomLink=bottomLink/Nbottom;
        
        rightLink=topRightState'*bottomRightState;
        Nright=abs(rightLink);
        rightLink=rightLink/Nright;
        
        topLink=topLeftState'*topRightState;
        Ntop=abs(topLink);
        topLink=topLink/Ntop;
        
        leftLink=bottomLeftState'*topLeftState;
        Nleft=abs(leftLink);
        leftLink=leftLink/Nleft;
        
%         if Nbottom<10^(-10) || Nright<10^(-10) || Ntop <10^(-10) || Nleft <10^(-10)
%            disp('Zero overlap of link states! Need to redefine mesh! NaNs replaced by 0.')
% %            break
%         end
        
        lattice_Field_Strength(j,l)=-1i*log(bottomLink*rightLink*topLink*leftLink);
                
% % % % % %                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %                  %%%%%%%%%%%% Turn values outside the 1st BZ into ZEROSs %%%%%%%%%%%%%%%
% % % % % %                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % 
% % % % % %                 %index of current path:
% % % % % %                 kx_index_current_path = j;
% % % % % % 
% % % % % %                 % index of left-right top-bottom edge of c/6 (M-K length) and b/2 (G-M length):
% % % % % %                 kx_axis_SQUARE_LEFT_edge_index = round(kgrid*1/4);
% % % % % %                 kx_axis_SQUARE_RIGHT_edge_index = round(kgrid*3/4);
% % % % % % 
% % % % % %          
% % % % % %                 %%top and bottom edge of hexagon
% % % % % %                 %left
% % % % % %                 ky_axis_TOP_LEFT_edge_index =     round( kgrid/2  + 2*( kx_index_current_path ) );
% % % % % %                 ky_axis_BOTTOM_LEFT_edge_index =  round( kgrid/2  - 2*( kx_index_current_path ) );
% % % % % %                 %right
% % % % % %                 ky_axis_TOP_RIGHT_edge_index =    round( kgrid/2  + 2*( kgrid - kx_index_current_path ) );
% % % % % %                 ky_axis_BOTTOM_RIGHT_edge_index = round( kgrid/2  - 2*( kgrid - kx_index_current_path ) );
% % % % % % 
% % % % % %                
% % % % % % 
% % % % % %                     if kx_index_current_path < kx_axis_SQUARE_LEFT_edge_index
% % % % % % 
% % % % % %                         if (ky_axis_TOP_LEFT_edge_index < l) || (ky_axis_BOTTOM_LEFT_edge_index > l) 
% % % % % %                             
% % % % % %                             lattice_Field_Strength(j,l)=0;
% % % % % % 
% % % % % %                         end
% % % % % % 
% % % % % % 
% % % % % %                     elseif   kx_axis_SQUARE_RIGHT_edge_index < kx_index_current_path 
% % % % % % 
% % % % % %                         if (ky_axis_TOP_RIGHT_edge_index < l) || (ky_axis_BOTTOM_RIGHT_edge_index > l) 
% % % % % % 
% % % % % %                             lattice_Field_Strength(j,l)=0;
% % % % % % 
% % % % % %                         end
% % % % % % 
% % % % % % 
% % % % % %                     end
% % % % % %               
% % % % % %                 
% % % % % %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    end
end

lattice_Field_Strength(isnan(lattice_Field_Strength)) = 0 ;
% lattice_Field_Strength(lattice_Field_Strength<0) = 0;

%Calculate summed imaginary part of lattice field strengths on all lattice
%points:
sum_Of_Imag_Field_Strength= sum(sum(imag(lattice_Field_Strength(:,:))));

if sum_Of_Imag_Field_Strength>10^(-10)
    
    disp('Field strength has large complex part. Consider precision.')
    
else
    
    %take real part of lattice field strength:
    lattice_Field_Strength=real(lattice_Field_Strength);
    
    %calculate Chern number:
    
    Chern_Number=sum(sum(lattice_Field_Strength))/(2*pi);
    disp(['The Chern number in band is: ', num2str(Chern_Number)]);
    
    [plaquetteY,plaquetteX]=ndgrid(linspace(1,BZsizeY-1,BZsizeY-1),linspace(1,BZsizeX-1,BZsizeX-1));
    figure
    mesh(plaquetteY,plaquetteX,lattice_Field_Strength)



end
             