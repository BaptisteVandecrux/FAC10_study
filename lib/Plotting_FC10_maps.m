function [FC10_change] = Plotting_FC10_maps (metadata, accum_thresh, T_thresh, ...
    XX,YY, FC10_map, DSA,LAPA,HAPA,GL,Ice,Firn, orig, vis)
      
    time_bin = [1949.50       1959.50       1969.50       1979.50 ...
        1989.50       1996.50      2009.5  2019.50];

    f=figure('Units','Normalized','outerposition',[0 0 1 0.8],'Visible',vis);
    ha =tight_subplot(1,2,-0.15,[0.03 0.17],[0.2 0.2]);

    set(f,'CurrentAxes',ha(1))
    hold on

    PlotBackground(GL,Ice,Firn);

    h_map = pcolor(XX,YY,FC10_map{2});
    h_map.LineStyle = 'none';
    max(max(FC10_map{2}))
    colbar =contourcmap('blues',0:250:4250,'colorbar','on');
    freezeColors
    freezeColors(colbar)
    ylabel(colbar,'Firn Capacity, top 10 m (kg m^{-2})','Interpreter','tex');
    colbar.Position(1) = 0.85;
    colbar.FontSize = 20;
    
    [~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',1);
    [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',1);
    [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',1); 
    plot(metadata.X, metadata.Y ,'xk','MarkerSize',7,'MarkerFaceColor','w'); 

    xlim(1.0e+05 *[-6.2164    8.4827])
    ylim(1.0e+06 *[-3.3439   -0.6732])

    set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
    daspect( [1 1 1])
             
    % Mapping evolution of pore space
set(f,'CurrentAxes',ha(2))
        FC10_change = FC10_map{2}-FC10_map{1};
    %     ps_change_10m(FAC_10_map{1,1} == 0) = NaN;
    hold on
    PlotBackground(GL,Ice,Firn);
    
    disp('Spatial average of FC10 change')   
    temp2 = FC10_change;
    temp2(LAPA == 0) = NaN;
    disp(nanmean(nanmean(temp2)))

    disp('Max FC10 change')
    disp(min(min(temp2)))
    h_map = pcolor(XX,YY,FC10_change);
    h_map.LineStyle = 'none';
    daspect( [1 1 1])

    colbar2 =contourcmap('tej2',-3500:250:0,'colorbar','on');
    freezeColors
    freezeColors(colbar2)
    [~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',2);
    [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',2);
    [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',2);

    met = metadata(and(ismember(metadata.Year,time_bin(6)+0.5:time_bin(7)+0.5),metadata.c_avg<accum_thresh),:);
    met = met(met.T_avg>T_thresh,:);
    plot(met.X, met.Y ,'ok','MarkerSize',0.00001,'MarkerFaceColor','w'); 
    h1 = plot(met.X, met.Y ,'ok','MarkerSize',5,'MarkerFaceColor','w'); 
    met = metadata(and(ismember(metadata.Year,time_bin(7)+0.5:time_bin(8)+0.5),metadata.c_avg<accum_thresh),:);
    met = met(met.T_avg>T_thresh,:);
    h2 = plot(met.X, met.Y ,'^k','MarkerSize',5,'MarkerFaceColor','w'); 

    h_leg = legendflex([h1 h2],{'1997-2008',...
        '2011-2017'},...
        'title','FAC_{10} observations');

    ylabel(colbar2,'\Delta FC_{10} (kg m^{-2})','Interpreter','tex');
    colbar2.FontSize = 20;
    xlim(1.0e+05 *[-4    1.5])
    ylim(1.0e+06 *[-3   -1.4])
    h_title = title('                LAPA \newline2011-2017 minus 1997-2008','Interpreter','tex');
    h_title.Units = 'Normalized';
    h_title.Position = [0.45          1.03             0];
    box on
    set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'layer','top')
    ha(2).Position(1) = 0.45;
    ha(2).Position(4) = 0.67;
    ha(2).Position(2) = 0.16;
    ha(1).Position(1) = 0.183;

    h_leg.Units = 'normalized';
    h_leg.Position(1) = 0.565;
    h_leg.Position(2) = 0.008;

    colbar.Position(1)=0.455;
    colbar2.Position(1)=0.7;

    for i=2:2:size(colbar.YTickLabel,1)
        colbar.YTickLabel(i,:) = '    ';
    end
    for i=2:2:size(colbar2.YTickLabel,1)
        colbar2.YTickLabel(i,:) = '     ';
    end
            
    annotation(f,'textbox',...
        [0.27 0.78 0.04 0.05],...
        'String','a)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');
    annotation(f,'textbox',...
        [0.5 0.88 0.04 0.05],...
        'String','b)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');

    x = 0.35*[1 1];
    y = [0.88 0.67];
    temp1 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp1.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String',' DSA 1953-2017',...
        'FontSize',14,'FontWeight','bold','LineWidth',1.5)

    x = [0.29 0.335];
    y = [0.35 0.3];
    temp2 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp2.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String','LAPA\newline 2011-2017',...
        'Interpreter','tex','FontSize',14,'FontWeight','bold','LineWidth',1.5)

    x = [0.38 0.35];
    y = [0.15 0.18];
    temp3 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp3.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String','HAPA\newline 2010-2017',...
        'Interpreter','tex','FontSize',14,'FontWeight','bold','LineWidth',1.5)

    print(f,sprintf('Output/FC_map_%s',orig),'-dtiff')
end