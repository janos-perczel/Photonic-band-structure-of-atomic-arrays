function [ summed_XX,summed_YY,summed_ZZ,summed_XY_and_YX ] = Summed_free_terms( non_Bravais_lattice,triangular_lattice, n1Lim,n2Lim,kVec_x,kVec_y,k0_A,a_sp,a_ho,non_Bravais_phase,b_vector)
%Summed_terms: the momementum sum comes here
    
    [k_X_GRID,axis_1_grid,axis_2_grid]=ndgrid(kVec_x,linspace(-n1Lim,n1Lim,2*n1Lim+1),linspace(-n2Lim,n2Lim,2*n2Lim+1));
    [k_Y_GRID,~,~]                             =ndgrid(kVec_y,linspace(-n1Lim,n1Lim,2*n1Lim+1),linspace(-n2Lim,n2Lim,2*n2Lim+1));

    
    b_vector = cell2mat(b_vector);
    
    if triangular_lattice
        %reciprocal lattice vectors of TRIANGULAR LATTICE!!!
        b1=2*pi/a_sp.*[1,1/sqrt(3)];
        b2=2*pi/a_sp.*[1,-1/sqrt(3)];
        
        if non_Bravais_lattice == 1
            
            delta3 = -a_sp.*( 1/sqrt(3) ).*[0,1];
            
        elseif non_Bravais_lattice == 2
            
            d1 = a_sp/2.*[1/2,sqrt(3)/2-1/sqrt(3)];
            d2 = a_sp/2.*[-1/2,sqrt(3)/2-1/sqrt(3)];
            d3 = a_sp/2.*[0,-1/sqrt(3)];
            
            if sum(b_vector == [1,2]) == 2 %test for list equality
                delta3 = d2-d1;
            elseif sum(b_vector == [1,3]) == 2 
                delta3 = d3-d1;
            elseif sum(b_vector == [2,3]) == 2 
                delta3 = d3-d2;
            end
            
        end
        
    else %if it's a square lattice
        
        b1=2*pi/a_sp.*[1,0];
        b2=2*pi/a_sp.*[0,1];
        
        if non_Bravais_lattice == 1
            
            delta3 = a_sp.*[1,1]./2;
            
        elseif non_Bravais_lattice == 2
            
            d1 = a_sp*[-1/4, 1/4];
            d2 = a_sp*[-1/4,-1/4];
            d3 = a_sp*[ 1/4,-1/4];
            
            if sum(b_vector == [1,2]) == 2 %test for list equality
                delta3 = d2-d1;
            elseif sum(b_vector == [1,3]) == 2 
                delta3 = d3-d1;
            elseif sum(b_vector == [2,3]) == 2 
                delta3 = d3-d2;
            end
            
        end
    end
    

    
    %projection to x-axis
    G_x=b1(1).*axis_1_grid + b2(1).*axis_2_grid;
    %projection to y-axis
    G_y=b1(2).*axis_1_grid + b2(2).*axis_2_grid;
    
    %subtract Bloch vector:
    G_K_x=G_x-k_X_GRID; %need to write kIn_x
    G_K_y=G_y-k_Y_GRID; %need to write kIn_y
    G_K_2=G_K_x.^2+G_K_y.^2;

  
    %%%%integral values
    
    %the square root
    kzd = PosSqrt(k0_A.^2-G_K_2);
  
    %!!!
    
    %the error function
    erfi_value = Faddeeva_erfi( a_ho.*kzd );
    
    %!!!
    
    %integrals
    Int_0=exp(  -a_ho^2.*kzd.^2  ) .* pi./kzd .*( -1i + erfi_value );
    Int_2= - pi^(1/2)/a_ho + exp(  -a_ho^2.*kzd.^2  ) .* pi .* kzd .*( -1i + erfi_value );
     
     %!!!
   
    %%%%% Consider a non-Bravais lattice
    
    if      non_Bravais_phase== 0
        
                non_Bravais_phase_factor = 1;
        
    elseif  non_Bravais_phase==-1;
        
                non_Bravais_phase_factor = exp( -1i.*( delta3(1).*G_K_x + delta3(2).*G_K_y ) );
        
    elseif  non_Bravais_phase== 1;
        
                non_Bravais_phase_factor = exp( 1i.*( delta3(1).*G_K_x + delta3(2).*G_K_y ) );
        
    else
                disp('Incorrect non-Bravais phase factor entered. It should be from {-1,0,1}.')
                quit
    end
        
    %prefactor for free-space Weyl-decomposed Greens (after k^2 multiplication):
    Weyl_prefactor=1/(2*pi.*k0_A.^2).*exp(  -a_ho^2.*G_K_2  ).*non_Bravais_phase_factor;
        
    
   
    
    %calculate the tensor components separately
    ndgrid_final_values_XX = Weyl_prefactor.*( k0_A.^2 - G_K_x.^2 ).*Int_0;
    summed_XX=sum(sum(ndgrid_final_values_XX,2),3);
    
    ndgrid_final_values_XY_and_YX = Weyl_prefactor.*( -G_K_x.*G_K_y ).*Int_0;
    summed_XY_and_YX=sum(sum(ndgrid_final_values_XY_and_YX,2),3);
    
    ndgrid_final_values_YY = Weyl_prefactor.*( k0_A.^2 - G_K_y.^2 ).*Int_0;
    summed_YY=sum(sum(ndgrid_final_values_YY,2),3);
    
    ndgrid_final_values_ZZ = Weyl_prefactor.*( k0_A.^2.*Int_0 - Int_2 );
    summed_ZZ=sum(sum(ndgrid_final_values_ZZ,2),3);
    
    

end
