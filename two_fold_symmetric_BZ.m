function [ kVec_x, kVec_y ] = two_fold_symmetric_BZ( linear_path, a_sp, PATH_to_evaluate_along_kx, kgrid )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% K_BLOCH DEPENDENT QUANTITIES  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %reciprocal lattice vectors
    b1=2*pi/a_sp.*[1,0];
    b2=2*pi/a_sp.*[0,1];

    %points/path in Brillouin zone
    BZ_point_X=1/2.*b1;
    BZ_point_M=1/2.*b1+1/2.*b2;

    if linear_path
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        uVar=linspace(0,1,kgrid);
        %path Gamma ---> X
        Gamma_to_X_path=KroneckerProduct(uVar',BZ_point_X);
        %path X ---> M given by [X + u*(M-X)]
        X_to_M_path=KroneckerProduct(ones(1,kgrid)',BZ_point_X)+KroneckerProduct(uVar',BZ_point_M-BZ_point_X);
        %path M ---> Gamma
        M_to_Gamma_path=flipud(KroneckerProduct(uVar',BZ_point_M));

        %BZ points mesh along line:
        if PATH_to_evaluate_along_kx == 1
            kVec=Gamma_to_X_path;
        end
        if PATH_to_evaluate_along_kx == 2
            kVec = X_to_M_path;
        end
        if PATH_to_evaluate_along_kx == 3
            kVec = M_to_Gamma_path;
        end

        kVec_x=kVec(:,1);
        kVec_y=kVec(:,2);

    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        uVar=linspace(-1,1,kgrid);
        kx_axis_line = KroneckerProduct(uVar',BZ_point_X);
        ky_axis_line = KroneckerProduct(uVar',BZ_point_M-BZ_point_X);

        kVec_x = kx_axis_line(PATH_to_evaluate_along_kx,1)+ky_axis_line(:,1);
        kVec_y = ky_axis_line(:,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end


end

