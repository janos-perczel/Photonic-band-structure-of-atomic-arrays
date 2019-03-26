%%%% In this script we evaluate the band structure near a flat metal surface
%%%%Note that this code does not readily deal with eps_d=/=1.
%%%%NB: at each initial set-up, the following command is needed:
path('/Applications/MATLAB_R2017a.app/toolbox/matlab/Faddeeva-MATLAB',path)
clear all
tic

%lattice types:
triangular = 1; square =     2;
honeycomb =  3; NB_square =  4;
kagome =     5; lieb =       6;

%run on following lattice and atom spacing:
lattice_type = 1; %pick lattice type
Hermitian = 0;    %run Hermitian simulation
linear_plot = 0; %linear or 2D plot
kgrid = 300; %no of points aling k_B axis
B_Z = 1; %switch on B field?
lambda0=7.9; %corresponds to 790nm
k0_A = 2*pi/lambda0; k0_B = k0_A; k0_C = k0_A; %wavevectors
Gamma0_A_per_k0_A = 1.56e-8; %linewidth
Gamma0_A = Gamma0_A_per_k0_A*k0_A;

%%%% atom spacing per lambda and B-field:
if      lattice_type == 1
                    %Tringular:
    a_sp_per_lambda0 = 1/2;
    zeeman_magnitude_per_Gamma0=.3;
    
elseif  lattice_type == 2
                    %Square
    a_sp_per_lambda0 = 0.45;
    zeeman_magnitude_per_Gamma0=0.5;
    
elseif  lattice_type == 3
                    %Honeycomb
    a_sp_per_lambda0 = 0.05*sqrt(3);
    zeeman_magnitude_per_Gamma0=12;
    k0_B = k0_A + 0*Gamma0_A;
    
elseif  lattice_type == 4
                    %NB Square
    a_sp_per_lambda0 = 1/13;
    zeeman_magnitude_per_Gamma0=20;
    k0_B = k0_A + 30*Gamma0_A;
    
elseif  lattice_type == 5
                    %Kagome
    a_sp_per_lambda0 = 1/10;
    zeeman_magnitude_per_Gamma0=20;
    k0_B = k0_A + 0*Gamma0_A;
    k0_C = k0_A + 0*Gamma0_A;
    
elseif  lattice_type == 6
                    %Lieb
    a_sp_per_lambda0 = 1/10;
    zeeman_magnitude_per_Gamma0=12;
    k0_B = k0_A + 0*Gamma0_A;
    k0_C = k0_A + 0*Gamma0_A;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(lattice_type,2) == 1
    triangular_lattice = 1;
else
    triangular_lattice = 0;
end

if lattice_type <= 2
    non_Bravais_lattice = 0;
    no_of_atoms_per_cell = 1;
elseif lattice_type == 3 || lattice_type == 4
    non_Bravais_lattice = 1;
    no_of_atoms_per_cell = 2;
else
    non_Bravais_lattice = 2;
    no_of_atoms_per_cell = 3; 
end

 
%lattice spacing
a_sp=a_sp_per_lambda0*lambda0; %%i.e. super-cube has lattice side of length a_sp=5/k0.
%ground state oscillation cut-off
a_ho=a_sp/10;

%summation limits (heuristically 0.7*a_sp/a_ho )
AllLim=7;
nxLim=AllLim;
nyLim=AllLim;


%%%%%

k0_B_per_k0_A = k0_B/k0_A;
Gamma0_B = (k0_B_per_k0_A)^3*Gamma0_A;

%%%%%

k0_C_per_k0_A = k0_C/k0_A;
Gamma0_C = (k0_C_per_k0_A)^3*Gamma0_A;

k0_list = [k0_A,k0_B,k0_C];
Gamma0_list = [Gamma0_A,Gamma0_B,Gamma0_C];


B_X=0;
B_Y=0;

if      B_X
            zeeman_X=zeeman_magnitude_per_Gamma0*Gamma0_A;
            zeeman_Y=0;
            zeeman_Z=0;
    if  B_Y
            zeeman_Y=zeeman_magnitude_per_Gamma0*Gamma0_A;
        if B_Z
            zeeman_Z=zeeman_magnitude_per_Gamma0*Gamma0_A;
        end
    end
elseif B_Y
            zeeman_X=0;
            zeeman_Y=zeeman_magnitude_per_Gamma0*Gamma0_A;
            zeeman_Z=0;
        if B_Z
            zeeman_Z=zeeman_magnitude_per_Gamma0*Gamma0_A;
        end
elseif B_Z
            zeeman_X=0;
            zeeman_Y=0;
            zeeman_Z=zeeman_magnitude_per_Gamma0*Gamma0_A;
                
else
            zeeman_X=0;
            zeeman_Y=0;
            zeeman_Z=0;
end
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Find the determinant dips and Gamma values %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kVec_x_ALL,kVec_y_ALL,eigenvectors_ALL,eigenvalues_ALL ] = Find_determinant_dips_gen(non_Bravais_lattice, triangular_lattice, no_of_atoms_per_cell, linear_plot, kgrid,a_sp,a_ho,nxLim,nyLim,k0_list,Gamma0_list,zeeman_X,zeeman_Y,zeeman_Z,Hermitian);


toc



