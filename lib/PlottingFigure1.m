function [] = PlottingFigure1(metadata, T_thresh, accum_thresh, ...
    GL,Ice,Firn, DSA, HAPA, LAPA, poly_100, I, orig, vis)

    [x_map,y_map]=pixcenters(I);
    XX = repmat(x_map,length(y_map),1);
    YY = repmat(y_map',1,length(x_map));

    f = figure('Visible',vis,'outerposition',[1 0 15 20]);
    ha = tight_subplot(4,2,[0.01 0.07],[0.14 0.5],[0.17 0.17]);
    for i = 3:4
        set(ha(i),'Visible','off')
    end

    set(f,'CurrentAxes',ha(1))
    hold on

    [~, ind_sorted] = sort(metadata.FAC10);
    ind_sorted=flipud(ind_sorted);

    vert = [ 63, 182, 83 ]/255;
    rouge = [ 194, 78, 67 ]/255;
    jaune = [ 223, 199, 79 ]/255;

    patch([T_thresh T_thresh max(metadata.T_avg)+1 max(metadata.T_avg)+1],...
        [accum_thresh 1800 1800 accum_thresh],...
        vert)
    patch([min(metadata.T_avg)-1 min(metadata.T_avg)-1 T_thresh T_thresh ],...
        [-50 1800 1800 -50],...
        jaune)
    patch([T_thresh T_thresh max(metadata.T_avg)+1 max(metadata.T_avg)+1],...
        [-50 accum_thresh accum_thresh -50],...
        rouge)

    g = patch(poly_100(:,1),poly_100(:,2), RGB('light light gray'));
    g.FaceColor = 'w';
    g.LineWidth = 1;
    alpha(g,0.35)

    xlim([min(metadata.T_avg)-1 max(metadata.T_avg)+1])
    ylim([0 max(metadata.c_avg)+1])
    
    scatter(metadata.T_avg(ind_sorted),...
        metadata.c_avg(ind_sorted),...
        140,...
        0*metadata.FAC10(ind_sorted)...
        ,'o','fill','MarkerEdgeColor','w');

    h_s1 = scatter(metadata.T_avg(ind_sorted),...
        metadata.c_avg(ind_sorted),...
        140,...
        metadata.FAC10(ind_sorted)...ones(size(metadata.c_avg(ind_sorted))) * 0.8*[1 1 1], ...
        ...+ (metadata.c_avg(ind_sorted)>accum_thresh) * 0.8*[1 0 0],...0.*col(discretize(metadata.Year,time_bin),:),...
        ,'o','fill','MarkerEdgeColor','w');
    cc = colorbar('NorthOutside');
    cc.Label.String = 'FAC_{10} (m)';
    colormap([0 0.2 0.6;0.0158 0.212 0.606;0.031 0.225 0.612;...
        0.0476 0.238 0.619;0.0634920671582222 0.250793665647507 0.625396847724915;0.0793650820851326 0.263492077589035 0.631746053695679;0.095238097012043 0.276190489530563 0.638095259666443;0.111111111938953 0.288888901472092 0.644444465637207;0.126984134316444 0.30158731341362 0.650793671607971;0.142857149243355 0.314285725355148 0.657142877578735;0.158730164170265 0.326984137296677 0.6634920835495;0.174603179097176 0.339682549238205 0.669841289520264;0.190476194024086 0.352380961179733 0.676190495491028;0.206349208950996 0.365079373121262 0.682539701461792;0.222222223877907 0.37777778506279 0.688888907432556;0.238095238804817 0.390476197004318 0.69523811340332;0.253968268632889 0.403174608945847 0.701587319374084;0.269841283559799 0.415873020887375 0.707936525344849;0.28571429848671 0.428571432828903 0.714285731315613;0.30158731341362 0.441269844770432 0.720634937286377;0.31746032834053 0.45396825671196 0.726984143257141;0.333333343267441 0.466666668653488 0.733333349227905;0.349206358194351 0.479365080595016 0.739682555198669;0.365079373121262 0.492063492536545 0.746031761169434;0.380952388048172 0.504761934280396 0.752380967140198;0.396825402975082 0.517460346221924 0.758730173110962;0.412698417901993 0.530158758163452 0.765079379081726;0.428571432828903 0.54285717010498 0.77142858505249;0.444444447755814 0.555555582046509 0.777777791023254;0.460317462682724 0.568253993988037 0.784126996994019;0.476190477609634 0.580952405929565 0.790476202964783;0.492063492536545 0.593650817871094 0.796825408935547;0.507936537265778 0.606349229812622 0.803174614906311;0.523809552192688 0.61904764175415 0.809523820877075;0.539682567119598 0.631746053695679 0.815873026847839;0.555555582046509 0.644444465637207 0.822222232818604;0.571428596973419 0.657142877578735 0.828571438789368;0.58730161190033 0.669841289520264 0.834920644760132;0.60317462682724 0.682539701461792 0.841269850730896;0.61904764175415 0.69523811340332 0.84761905670166;0.634920656681061 0.707936525344849 0.853968262672424;0.650793671607971 0.720634937286377 0.860317468643188;0.666666686534882 0.733333349227905 0.866666674613953;0.682539701461792 0.746031761169434 0.873015880584717;0.698412716388702 0.758730173110962 0.879365086555481;0.714285731315613 0.77142858505249 0.885714292526245;0.730158746242523 0.784126996994019 0.892063498497009;0.746031761169434 0.796825408935547 0.898412704467773;0.761904776096344 0.809523820877075 0.904761910438538;0.777777791023254 0.822222232818604 0.911111116409302;0.793650805950165 0.834920644760132 0.917460322380066;0.809523820877075 0.84761905670166 0.92380952835083;0.825396835803986 0.860317468643188 0.930158734321594;0.841269850730896 0.873015880584717 0.936507940292358;0.857142865657806 0.885714292526245 0.942857146263123;0.873015880584717 0.898412704467773 0.949206352233887;0.888888895511627 0.911111116409302 0.955555558204651;0.904761910438538 0.92380952835083 0.961904764175415;0.920634925365448 0.936507940292358 0.968253970146179;0.936507940292358 0.949206352233887 0.974603176116943;0.952380955219269 0.961904764175415 0.980952382087708;0.968253970146179 0.974603176116943 0.987301588058472;0.98412698507309 0.987301588058472 0.993650794029236;1 1 1])

    h_leg = legend([g h_s1],'Firn area',...
        'Firn cores',...
        'Location','NorthWest');

    set(gca,'Ticklength',[0.08 0.16]/2,'layer','top','yaxislocation','right')
    box on
    ylabel('$\mathrm{\overline{\dot{c}} \:(mm w.eq. yr^{-1})}$','Interpreter','latex')
    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex') 

    % Explaining dataset division
%     ind_HighAcc_HighTemp = and(metadata.c_avg>accum_thresh,metadata.T_avg>T_thresh);
%     ind_Rest = or(metadata.c_avg<=accum_thresh,metadata.T_avg<=T_thresh);
    ha(1).Position = [0.45 0.42 0.4 0.3];
    ylim([-50       1800])

    set(f,'CurrentAxes',ha(6))
    ha(6).Visible = 'off';

    set(f,'CurrentAxes',ha(8))
    ha(8).Visible = 'off';

    set(f,'CurrentAxes',ha(5))
         ha(5).Position = [ 0.17   0.2  0.68   0.125] ;
        hold off
        ind_jaune = metadata.T_avg<T_thresh;
        scatter(metadata.T_avg(ind_jaune), metadata.FAC10(ind_jaune),30,jaune,'fill')
        hold on

        ind_rouge = and(metadata.T_avg>T_thresh, metadata.c_avg<accum_thresh);
        scatter(metadata.T_avg(ind_rouge), metadata.FAC10(ind_rouge),30,rouge,'fill')
        ind_vert = and(metadata.T_avg>T_thresh, metadata.c_avg>accum_thresh);
        scatter(metadata.T_avg(ind_vert),metadata.FAC10(ind_vert),30,vert,'fill')
        % plot([T_thresh T_thresh],[0 7],'--k','LineWidth',2)
        xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
        box on
        ax=gca;
        set(gca,'Ticklength',[0.08 0.16]/2.5,'XMinorTick','on','YMinorTick','on',...
            'XAxisLocation','bottom','YAxisLocation','left','layer','top')
             ax.XAxis.MinorTickValues = ax.XTick(1):1:ax.XTick(end);

        axis tight
        xlimits = get(gca,'XLim');

        ylim([0 7])
        h_ylab = ylabel('FAC_{10} (m)','Interpreter','tex');

    set(f,'CurrentAxes',ha(7))
    ha(7).Visible = 'off';

    for i=1:8
        set(f,'CurrentAxes',ha(i))
        set(gca,'FontSize',14,'FontName','Times New Roman')
    end

    set(f,'CurrentAxes',ha(2))
    ha(2).Position = [-0.005 0.33 0.53 0.53];
    hold on
    [p4, p5] = PlotBackground(GL,Ice,Firn);

    DSA_temp = DSA;
    HAPA_temp = HAPA;
    LAPA_temp = LAPA;
            DSA_temp(DSA == 0)=NaN;
            HAPA_temp(HAPA == 0)=NaN;
            LAPA_temp(LAPA == 0)=NaN;
    p1 = surf(XX,YY,DSA_temp*0, DSA_temp);
    p2 = surf(XX,YY,HAPA_temp*0, HAPA_temp+1);
    p3 = surf(XX,YY,LAPA_temp*0, LAPA_temp+2);

    p1.LineStyle = 'none';
    p2.LineStyle = 'none';
    p3.LineStyle = 'none';
    p1.FaceColor = jaune;
    p2.FaceColor = vert;
    p3.FaceColor = rouge;
    [~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',1);
    [~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',1);
    [~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',1); 
    plot(metadata.X, metadata.Y,'xk',...
        'MarkerFaceColor','w','MarkerSize',7)
                xlim(1.0e+05 *[-6.2164    8.4827])
                ylim(1.0e+06 *[-3.3439   -0.6732])
                set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
             h_leg_map =  legendflex([p1 p2 p3 p4 p5], {'Dry snow area','High-accumulation percolation area',...
                 'Low-accumulation percolation area', 'Ablation area', 'Land'},...
                            'ref', gcf, ...
                           'anchor', {'n','n'}, ...
                           'buffer',[0.3 -20], ...
                           'nrow',3,...
                           'box','off',...
                           'fontsize',13);  
    daspect( [1 1 1])

    % Create textbox
    annotation(f,'textbox',...
        [0.07 0.78 0.04 0.07],...
        'String','a)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');

    % Create textbox
    annotation(f,'textbox',...
        [0.43 0.82 0.04 0.05],...
        'String','b)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',19,...
        'FitBoxToText','off');

    % Create textbox
    annotation(f,'textbox',...
        [0.17 0.22 0.04 0.05],...
        'String','c)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');

    print(f,sprintf('Output/data_presentation_%s',orig),'-dtiff','-r600')
    print(f,sprintf('Output/data_presentation_%s',orig),'-dpdf','-r0')
end