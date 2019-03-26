figure
tic
color_plot = 1;

no_of_paths = size(eigenvalues_ALL,1);

y_max = 1.4*max((squeeze(real(eigenvalues_ALL(:)))-k0_A)./Gamma0_A);
y_min = 1.1*min((squeeze(real(eigenvalues_ALL(:)))-k0_A)./Gamma0_A);

% y_max = 150;
% y_min = -150;
% 
% y_max = -25.0;
% y_min = -35.50;

y_max = 5.0;
y_min = -5.0;

% y_max = 500;
% y_min = -500;

no_of_color_pts = 200;
col_map=colormap(parula(no_of_color_pts));
max_gamma = 3;

x_axis = linspace(0,1,kgrid);


for path_index=1:no_of_paths
    ax(path_index)=subplot(1,no_of_paths,path_index);
    if path_index==1
        label_y_axis=0;
        add_colorbar=0;
        ticklabels_Y=1;
        reverse_k_plot=1;
    else
        label_y_axis=0;
        add_colorbar=0;
        ticklabels_Y=0;
        if path_index==2
           reverse_k_plot = 1; 
        end
    end    
    
    if color_plot == 0
            if reverse_k_plot
                x_ax = fliplr(x_axis);
            else
                x_ax = x_axis(k_ind);
            end
            
            plot(x_ax,(squeeze(real(eigenvalues_ALL(path_index,:,:)))-k0_A)./Gamma0_A,'LineWidth',3 )
    else
        tic
        for k_ind = 1:kgrid
            for band_ind = 1:2*no_of_atoms_per_cell

                index_of_color = round(no_of_color_pts.*((-imag(eigenvalues_ALL(path_index,k_ind,band_ind))./Gamma0_A)/max_gamma))+1;
                index_of_color(index_of_color>200) = 200;
                edge_color = col_map(index_of_color,:);
                fill_color = 'k';
                marker_size = 15; 

                if reverse_k_plot
                    x_ax = x_axis(kgrid+1-k_ind);
                else
                    x_ax = x_axis(k_ind);
                end
                plot(x_ax,(squeeze(real(eigenvalues_ALL(path_index,k_ind,band_ind)))-k0_A)./Gamma0_A,'.','MarkerEdgeColor',edge_color,'MarkerFaceColor',fill_color,'MarkerSize',marker_size )
                hold on

            end
        end
        toc
    end
        
    if label_y_axis
        ylabel('\delta/\Gamma_{0}','FontSize', 15)
    end

    pbaspect([1 2 1])

    hold off
    ylim([y_min y_max])

    if ~ticklabels_Y
        set(gca,'YTickLabel','')
    end
end




%create colorbar
cbr=colorbar;
%calculate and use max gamma of range
caxis([0 max_gamma]);
%set ticks
set(cbr,'YTick',0:0.25:max_gamma);

%reposition panels
pos_n_siz=zeros(no_of_paths,4);
left_position=0.2;
for sub_pl=[1,3,2]
      pos_n_siz(sub_pl,:)=get(ax(sub_pl), 'Position');
      set(ax(sub_pl),'Position',[left_position pos_n_siz(sub_pl,2) 0.8*pos_n_siz(sub_pl,3) pos_n_siz(sub_pl,4)]);
      %iterative keep track of right edge of previous plot
      left_position=left_position+0.83*pos_n_siz(sub_pl,3);
end                    

%reposition colorbar and label with max gamma.
set(cbr, 'Position', [1.05*left_position .3 .03 .4])
% title(cbr,['max \Gamma: ',num2str(max_gamma_plotted/Gamma0_A,2)],'Fontsize',12)

toc