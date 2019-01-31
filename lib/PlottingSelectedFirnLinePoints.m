function [] = PlottingSelectedFirnLinePoints(metadata, T_thresh, accum_thresh, ind1,...
    poly_100, T_FL_obs, c_FL_obs,XX,YY,DSA,LAPA,HAPA,X_FL_obs,Y_FL_obs, ...
    c_max, source_temp_accum,GL,Ice,Firn, orig, vis)

    col=lines(length(ind1));
    f = figure('Visible',vis,'outerposition',[1 1 20 15]);
    ha = tight_subplot(1,2, 0.01,0.07,[0.07 0]);
    set(gcf, 'CurrentAxes',ha(1))

        a = patch([min(metadata.T_avg)-1 min(metadata.T_avg)-1 T_thresh T_thresh ],...
            [-50 c_max c_max -50],...
            RGB('jaune'));
       b=  patch([T_thresh T_thresh 0 0],...
            [-50 accum_thresh accum_thresh -50],...
            RGB('rouge'));
        c=patch([T_thresh T_thresh 0 0],...
            [accum_thresh c_max c_max accum_thresh],...
            RGB('vert'));
        g = patch(poly_100(:,1),poly_100(:,2), RGB('light light gray'));
        g.FaceColor = 'w';
        g.LineWidth = 1;
        alpha(g,0.35)
        hold on
        d=plot(T_FL_obs,c_FL_obs,'.k');
        dummy=plot(1,NaN,'w');

        e=scatter(T_FL_obs(ind1), c_FL_obs(ind1),100,col,...
            'fill','LineWidth',1,'MarkerEdgeColor','k'); 
    %     dx = 1; dy = 1; % displacement so the text does not overlay the data points
    %     text(T_FL_obs(ind1), C_FL_obs(ind1) , num2str(ind1));
        h = plot( [1 1], [NaN NaN],'xk','MarkerSize',7,'MarkerFaceColor','w'); 

        ha(1).Units = 'normalized';
        ha(1).Position = [0.17 0.13 0.4 0.5];

         ylabel('$\mathrm{\overline{\dot{c}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
        xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
        axis tight square
        ylim([0 c_max])
        if source_temp_accum == 2
            x = [-3 +3]-9;
            y= 2300*[1 1];
        else
            x = [-3 +3]-6;
            y= 4200*[1 1];
        end
        pos = ha(1).Position;
        ta1 = annotation('doublearrow',...
            (x + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1),...
            (y - min(ylim))/diff(ylim) * pos(4) + pos(2));
        ta2 = annotation('doublearrow',[0.16 0.21]-0.07, 0.685*[1 1]);
            ta1.LineWidth = 1.5;
        ta2.LineWidth = 1.5;

        h_leg = legend([a b c d h e dummy],'DSA', 'LAPA', 'HAPA',...
        'Locations where the firn line was detected',...
        'FAC_{10} observations',...
        'Selected firn line locations used as FAC_{10} observations',...
        'Perturbation applied in the sensitivity analysis','Interpreter','tex');
        h_leg.FontSize = 12;
        legend boxoff
        h_leg.Units = 'normalized';
        h_leg.Position(1:2) = [0.08 0.68];


        box on
        set(gca,'Ticklength',[0.08 0.16]/2,'layer','top')

    set(gcf, 'CurrentAxes',ha(2))
        hold on
        PlotBackground(GL,Ice,Firn);
        DSA_temp = DSA;
        LAPA_temp = LAPA;
        HAPA_temp = HAPA;
                DSA_temp(DSA == 0)=NaN;
                HAPA_temp(HAPA == 0)=NaN;
                LAPA_temp(LAPA == 0)=NaN;
        p1 = surf(XX,YY,DSA_temp*0, DSA_temp);
        p2 = surf(XX,YY,HAPA_temp*0, HAPA_temp+1);
        p3 = surf(XX,YY,LAPA_temp*0, LAPA_temp+2);

        p1.LineStyle = 'none';
        p2.LineStyle = 'none';
        p3.LineStyle = 'none';
        p1.FaceColor = RGB('jaune');
        p2.FaceColor = RGB('vert');
        p3.FaceColor = RGB('rouge');
        [~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',1.5);
        [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',1.5);
        [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',1.5); 
        plot(metadata.X, metadata.Y ,'xk','MarkerSize',7,'MarkerFaceColor','w'); 

        scatter(X_FL_obs(ind1), Y_FL_obs(ind1),100,col,...
        'fill','LineWidth',1,'MarkerEdgeColor','k'); 
    % dx = 0.2; dy = 0.2; % displacement so the text does not overlay the data points
    % text(X_FL_obs(ind1), Y_FL_obs(ind1) , num2str(ind1));
    daspect( [1 1 1])
    if source_temp_accum == 2
        xlim(1.0e+05 *[-1.5    3.5])
    else
            xlim(1.0e+05 *[-2.6    2.7])
    end
        ylim(1.0e+06 *[-3.35   -2.4])
        box on
        set(gca,'XTickLabel','','YTickLabel','','layer','top')
    ha(2).Position(1:2) = [0.56 0.08];

        annotation(f,'textbox',...
            [0.22 0.58 0.04 0.05],...
            'String','a)',...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',20,...
            'FitBoxToText','off');
            annotation(f,'textbox',...
            [0.87 0.13 0.04 0.05],...
            'String','b)',...
            'LineStyle','none',...
            'FontWeight','bold',...
            'FontSize',20,...
            'FitBoxToText','off');
    print(f,sprintf('./Output/FirnLine_%s',orig),'-dtiff')
    print(f,sprintf('./Output/FirnLine_%s',orig),'-dpdf')
end