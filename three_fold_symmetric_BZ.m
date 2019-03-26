function [ kVec_x, kVec_y ] = three_fold_symmetric_BZ( linear_path, a_sp, PATH_to_evaluate_along_kx, kgrid )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% K_BLOCH DEPENDENT QUANTITIES  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %reciprocal lattice vectors of TRIANGULAR LATTICE!!!
    b1=2*pi/a_sp.*[1,1/sqrt(3)];
    b2=2*pi/a_sp.*[1,-1/sqrt(3)];
    
    %vectors to high symmetry points:
    Gamma_M_vector = 2*pi/a_sp*1/sqrt(3).*[0,1];
    M_K_vector     = 2*pi/a_sp*1/3      .*[1,0]; 

    %points/path in Brillouin zone
    BZ_point_M = Gamma_M_vector;
    BZ_point_K = Gamma_M_vector + M_K_vector;
    
    
    if linear_path
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %vector pointing from X to M

        uVar=linspace(0.001,1,kgrid);
        %path Gamma ---> M
        Gamma_to_M_path=KroneckerProduct(uVar',BZ_point_M);
        %path M ---> K 
        M_to_K_path=KroneckerProduct(ones(1,kgrid)',BZ_point_M)+KroneckerProduct(uVar',M_K_vector);
        %path K ---> Gamma
        K_to_Gamma_path=flipud(KroneckerProduct(uVar',BZ_point_K));

        %BZ points mesh along line:
        if PATH_to_evaluate_along_kx == 1
            kVec=Gamma_to_M_path;
        end
        if PATH_to_evaluate_along_kx == 2
            kVec = M_to_K_path;
        end
        if PATH_to_evaluate_along_kx == 3
            kVec = K_to_Gamma_path;
        end

        kVec_x=kVec(:,1);
        kVec_y=kVec(:,2);
        
    else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %sweep rectangle containing the honeycomb.
        uVar=linspace(-2,2,kgrid);
        kx_axis_line = KroneckerProduct(uVar',M_K_vector);
        uVar=linspace(-1,1,kgrid);
        ky_axis_line = KroneckerProduct(uVar',Gamma_M_vector);

        kVec_x = kx_axis_line(PATH_to_evaluate_along_kx,1)+ky_axis_line(:,1);
        kVec_y = ky_axis_line(:,2);



         %%%%%%%%%%%% Turn values outside the 1st BZ into NaNs %%%%%%%%%%%%%%%


        %index of current path:
        kx_index_current_path = PATH_to_evaluate_along_kx;

        % index of left-right top-bottom edge of c/6 (M-K length) and b/2 (G-M length):
        kx_axis_SQUARE_LEFT_edge_index = round(kgrid*1/4);
        kx_axis_SQUARE_RIGHT_edge_index = round(kgrid*3/4);


        %%top and bottom edge of hexagon
        %left
        ky_axis_TOP_LEFT_edge_index =     round( kgrid/2  + 2*( kx_index_current_path ) );
        ky_axis_BOTTOM_LEFT_edge_index =  round( kgrid/2  - 2*( kx_index_current_path ) );
        %right
        ky_axis_TOP_RIGHT_edge_index =    round( kgrid/2  + 2*( kgrid - kx_index_current_path ) );
        ky_axis_BOTTOM_RIGHT_edge_index = round( kgrid/2  - 2*( kgrid - kx_index_current_path ) );

        for y_index = 1:kgrid

            if kx_index_current_path < kx_axis_SQUARE_LEFT_edge_index

                if (ky_axis_TOP_LEFT_edge_index < y_index) || (ky_axis_BOTTOM_LEFT_edge_index > y_index) 

                    kVec_x(y_index)=NaN;

                end


            elseif   kx_axis_SQUARE_RIGHT_edge_index < kx_index_current_path 

                if (ky_axis_TOP_RIGHT_edge_index < y_index) || (ky_axis_BOTTOM_RIGHT_edge_index > y_index) 

                    kVec_x(y_index)=NaN;

                end


            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end


end

