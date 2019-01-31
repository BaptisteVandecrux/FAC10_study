function [] = PlottingAllRegions_3D(metadata,gridx, gridy, FAC_fit_all_post2010_mid,...
    FAC_fit_all_pre2010_mid,  x1, y1, z1, x2, y2, z2, ind_in_firn, vis,orig)
time_bin = [1949.50       1959.50       1969.50       1979.50 ...
    1989.50       1996.50      2009.5  2019.50];
col = lines(length(time_bin)-1);

    f = figure('Visible', vis, 'outerposition',[1 0 30 30]);
    ha = tight_subplot(1,2,0.18,0.20,0.15);
for i = 1:2

    set(f,'CurrentAxes',ha(i))
    hold on

        % plotting the firn line-dervied FAC
        switch i
            case 1
            % plotting the core derived FAC
%             scatter3(metadata.T_avg(metadata.Year<2010), ...
%             metadata.c_avg(metadata.Year<2010), ...
%             metadata.FAC10(metadata.Year<2010),...
%             100,col(discretize(metadata.Year(metadata.Year<2010),time_bin),:),'fill') 
            scatter3(x1,  y1, z1,...
            100,col(7,:),'fill') 
            temp = FAC_fit_all_pre2010_mid;
            temp(~ind_in_firn) = NaN;
            h_s2 = surf(gridx, gridy,temp,0*gridx+1);
            [~, h_c2 ]= contour3(gridx, gridy,temp ,15);
            title('pre 2010')

            case 2
            % plotting the core derived FAC
%             scatter3(metadata.T_avg(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),...
%             metadata.c_avg(or(metadata.Year>=2010,metadata.T_avg<T_thresh)), ...
%             metadata.FAC10(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),...
%             100,col(discretize(metadata.Year(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),time_bin),:),'fill') 
            scatter3(x2,  y2, z2,...
            100,col(7,:),'fill')
            temp = FAC_fit_all_post2010_mid;
             temp(~ind_in_firn) = NaN;
            h_s2 = surf(gridx, gridy,temp,0*gridx+1);
            [~, h_c2 ]= contour3(gridx, gridy,temp ,15);
            title('post 2010')
            
        end
        alpha(h_s2, 0.2)

        h_c2.Color = 'r';
        h_c2.LineWidth = 1.5;
        shading interp
        colormap hsv
%     ylim([-2 2500])
    box on
    ylabel('$\mathrm{\overline{\dot{c}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    zlabel('FAC_{10} (m)','Interpreter','tex')
    view(gca,[5 20]);
ylim([0,3100])
xlim([-32, -4])
end
    print(f,sprintf('./Output/all_areas_%s_3d',orig),'-dtiff')

end