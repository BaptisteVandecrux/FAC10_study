function [] = PlottingDifferenceBoxMAR(XX_box,YY_box,diff_MAR_Box_1,diff_MAR_Box_2,...
    XX, YY, DSA, HAPA, LAPA, GL, Ice, Firn, vis)

f=figure('Visible', vis, 'outerposition',[0 0 30 20],'Visible',vis);
    ha =tight_subplot(1,2,0.0001,[0.02 0.14],[0.07 0.3]);
    
set(f,'CurrentAxes',ha(2))
    hold on
    PlotBackground(GL,Ice,Firn); 
       
    h_map = pcolor(XX_box,YY_box,diff_MAR_Box_1);
    h_map.LineStyle = 'none';
    [~, ~] = contour(XX,YY,DSA+HAPA+LAPA, [1 1],'Color','k', 'LineWidth',2);
%     [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',3);
%     [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',3);
    title('LAPA\newline1998-2008','Interpreter','tex')

    daspect( [1 1 1])
    xlim(1.0e+05 *[-6.2164    1.1])
    ylim(1.0e+06 *[-3   -1.3])
    box on
    colbar =contourcmap('hsv',0:0.3:2.5,'colorbar','on');
    ylabel(colbar,'Absolute difference between Box13-derived  \newline           and MAR-derived FAC_{10}(m)','Interpreter','tex');
    set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'layer','top')

set(f,'CurrentAxes',ha(1))
    ha(2).Position(1) = 0.4;

    hold on
    PlotBackground(GL,Ice,Firn);            
    h_map = pcolor(XX_box,YY_box,diff_MAR_Box_2);
    h_map.LineStyle = 'none';
    [~, ~] = contour(XX,YY,DSA+LAPA+HAPA, [1 1],'Color','k', 'LineWidth',2);
%     [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',3);
%     [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',3);

    daspect( [1 1 1])
    xlim(1.0e+05 *[-6.2164    8.4827])
    ylim(1.0e+06 *[-3.3439   -0.6732])

    colbar2 =contourcmap('jet2',0:0.3:2.5,'colorbar','on');
    colbar2.Position(2) = 2;
    set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')

    uistack(ha(1),'top')
    uistack(ha(2),'top')

    ha(1).Units = 'normalized';
    ha(2).Units = 'normalized';
    ha(1).Position = [0.21 0.03 0.46 0.80];
    ha(2).Position = [0.525 0.07 0.2 0.4];
    colbar.Position(1:2) = [0.69 0.05];
    
    x = 0.435*[1 1];
    y = [0.88 0.67];
    temp1 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp1.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String',' DSA 1953-2017',...
        'FontSize',14,'FontWeight','bold','LineWidth',1.5)

    x = [0.34 0.41];
    y = [0.35 0.3];
    temp2 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp2.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String','LAPA\newline 2010-2017',...
        'Interpreter','tex','FontSize',14,'FontWeight','bold','LineWidth',1.5)

    x = [0.47 0.435];
    y = [0.15 0.18];
    temp3 = annotation('textarrow',x,y,'String','','LineWidth',3,'Color','w');
    temp3.HeadStyle = 'plain';
    annotation('textarrow',x,y,'String','HAPA\newline 2010-2017',...
        'Interpreter','tex','FontSize',14,'FontWeight','bold','LineWidth',1.5)

    print(f,('Output/uncertainty_FAC_map_MAR_vs_Box13'),'-dtiff')
end