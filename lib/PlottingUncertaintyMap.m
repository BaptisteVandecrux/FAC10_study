function [] = PlottingUncertaintyMap(XX,YY, delta_FAC10, FAC10_map, ...
    DSA, HAPA, LAPA, GL, Ice, Firn, orig, vis)
%%
c_map = 'parula';
% scale = 5:1:20;
scale = 0:5:100;

      f=figure('Units','Normalized','outerposition',[0 0 0.7 0.8],'Visible',vis);
        ha =tight_subplot(1,2,-0.15,[0.03 0.17],[0.01 0.11]);

    set(f,'CurrentAxes',ha(1))    
        hold on
        PlotBackground(GL,Ice,Firn);

        h_map = pcolor(XX,YY,delta_FAC10{2}./FAC10_map{2}*100);
        h_map.LineStyle = 'none';
        set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')

        [~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',2);
        [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',2);
        [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',2);

        daspect( [1 1 1])
        xlim(1.0e+05 *[-6.2164    8.4827])
        ylim(1.0e+06 *[-3.3439   -0.6732])
        colbar =contourcmap(c_map,scale,'colorbar','on');
    freezeColors
    freezeColors(colbar)
        colbar.Position(1) = 0.81;
        for i=2:2:size(colbar.YTickLabel,1)
            colbar.YTickLabel(i,:) = '   ';
        end
        colbar.FontSize = 20;
        ylabel(colbar,'        Relative uncertainty \newline  on estimated FAC_{10} (%)','Interpreter','tex');
        set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')

        set(f,'CurrentAxes',ha(2))
        hold on
        PlotBackground(GL,Ice,Firn);

        h_map = pcolor(XX,YY,delta_FAC10{1}./FAC10_map{1}*100);
        h_map.LineStyle = 'none';
        [~, h_CA] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',2);
        [~, h_WHA] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',2);
        [~, h_WLA] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',2);
        daspect( [1 1 1])
        xlim(1.0e+05 *[-6.2164    8.4827])
        ylim(1.0e+06 *[-3.3439   -0.6732])
        title(' LAPA  \newline1998-2008','Interpreter','tex')
        colbar2 =contourcmap(c_map,scale,'colorbar','on');
    freezeColors
    freezeColors(colbar)
        colbar2.Position(2) = 1.5;
                            set(gca,'XTickLabel','','YTickLabel','')
        xlim(1.0e+05 *[-4    1.5])
        ylim(1.0e+06 *[-3   -1.4])
        uistack(ha(2),'top')

    ha(1).Units = 'normalized';
    ha(2).Units = 'normalized';

    ha(1).Position = [ 0.21 0.03 0.46 0.80];
    ha(2).Position = [0.51 0.07 0.2 0.4];
    colbar.Position(1:2) = [0.67 0.05];

    ha(2).Box = 'on';
    set(ha(2) ,'Layer','top')

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

    annotation(f,'textbox',...
        [0.32 0.88 0.04 0.05],...
        'String','(c)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');

    print(f,sprintf('Output/uncertainty_FAC_map_overall_%s',orig),'-dtiff','-r400')
    %%
end