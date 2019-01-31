function [RMSE, ME] = PlottingDSA(metadata, FAC_fit_DSA, Coef, T_thresh, vis, orig)
xx = linspace(min(metadata.c_avg),max(metadata.c_avg));
yy = linspace(min(metadata.T_avg),max(metadata.T_avg));
time_bin = [1949.50       1959.50       1969.50       1979.50 ...
    1989.50       1996.50      2009.5  2019.50];
col = lines(length(time_bin)-1);

ME = mean(metadata.FAC10-FAC_fit_DSA(metadata.c_avg, metadata.T_avg));
RMSE = sqrt(mean((metadata.FAC10 - FAC_fit_DSA(metadata.c_avg, metadata.T_avg)).^2));


f = figure('Visible',vis,'outerposition',[1 0 20 30]);
ha = tight_subplot(2,1,0.1,[0.1 0.3],0.2);

set(f,'CurrentAxes',ha(1))
    hold on
    h_NH = plot(yy,FAC_fit_DSA(xx,yy),'k','LineWidth',2);

    % plotting the core derived FAC
    time_bin2 = 1949.50 :10:2019.50 ;
    h=[];
    for i = 1: length(time_bin2)-1
        met = metadata(and(metadata.Year<time_bin2(i+1),metadata.Year>=time_bin2(i)),:);

        res = abs(FAC_fit_DSA(met.c_avg, met.T_avg) ...
            - met.FAC10);
        [~, ind_sorted] = sort(res, 1,'descend');
        if length(met.c_avg(ind_sorted))>0
            h(length(h)+1) =  scatter(met.T_avg(ind_sorted),...
                met.FAC10(ind_sorted), 30,... 100*res(ind_sorted),...
                    'MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)) ;
        else
            sprintf('Nothing in %0.1f - %0.1f\n',time_bin2(i),time_bin2(i+1))
            continue
        end
    end

    % plotting uncertainty
    h_unc = patch([yy fliplr(yy)], ...
        [2*RMSE + FAC_fit_DSA(xx,yy), fliplr(-2*RMSE + FAC_fit_DSA(xx,yy))],...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
uistack(h_unc,'bottom')
met = metadata;
ind = and(met.T_avg>=T_thresh,met.c_avg<=600);
met(ind,:) = [];

    T_bins = linspace(min(met.T_avg),max(met.T_avg),5);
    for i=1:length(T_bins)
        h_bin = plot(T_bins(i)*[1 1],[3 6.5],'Color',RGB('light gray'));
        uistack(h_bin,'bottom')
    end

    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    ylabel('FAC_{10} (m^{3} m^{-2})','Interpreter','tex')

    leg_text_1 = {'Linear regression','$\pm 2 \times RMSD$','Temperature bins used for regression'};
    leg_text_2 = {'1950 - 1959','1960 - 1969', '1980 - 1989','1990 - 1999','2000 - 2009','2010 - 2017'};
        
%     ha(1).Position =  [0.30  0.30 0.40   0.40];

    h_leg_1 = legendflex([h_NH h_unc h_bin], leg_text_1,...
        'ncol',1,'box','off','interpreter','latex');
    h_leg_2 = legendflex(h, leg_text_2,...
        'ncol',3,'box','off',...
        'title','Survey year:');

    h_leg_1.Box = 'off';
    h_leg_2.Box = 'off';
    h_leg_1.Units = 'normalized';
    h_leg_2.Units = 'normalized';
    h_leg_1.Position(1:2) = [0.25 0.83];
    h_leg_2.Position(1:2) = [0.18 0.73];
        set(gca,'Ticklength',[0.08 0.16]/4,'layer','top')
 
    axis tight
    ylim([3 6.5])
     box on
     set(gca,'layer','top')
     
set(f,'CurrentAxes',ha(2))
    met = metadata(metadata.T_avg<T_thresh,:);
    res = (FAC_fit_DSA(met.c_avg, met.T_avg) ...
                - met.FAC10); 
     temp = datevec(datenum(met.Date));

    scatter(met.Year+temp(:,2)./12,res,20,[0 0 0],'fill')
    hold on
plot([min(met.Year) max(met.Year)], 0*RMSE*[1 1],'--k')
h_p = patch([min(met.Year) max(met.Year)+1 max(met.Year)+1 min(met.Year)], ...
    [2*RMSE 2*RMSE -2*RMSE -2*RMSE],RGB('light light gray'));
uistack(h_p,'bottom')
h_p.LineStyle = 'none';
axis tight
    box on
    xlabel('Year')
    ylabel('Estimated - observed\newline     FAC_{10} (m)','interpreter','tex')
    ylim([-1 1])
        set(gca,'Ticklength',[0.08 0.16]/4,'layer','top')

annotation(f,'textbox',...
    [0.22 0.5 0.04 0.05],...
    'String','a)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(f,'textbox',...
    [0.22 0.29 0.04 0.05],...
    'String','b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
print(f,sprintf('./Output/DSA_%s',orig),'-dtiff')
print(f,sprintf('./Output/DSA_%s',orig),'-dpdf')
end
