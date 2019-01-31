function [out] = Plotting_FAC_RCM_Comp(FAC_RCM, FAC_obs, ...
    T_avg, c_avg, unc_FAC)

    T_thresh = -19;
    accum_thresh = 600;
    ind_DSA = T_avg<T_thresh;
    ind_LAWSA = and(T_avg>=T_thresh, c_avg<accum_thresh);
    ind_HAWSA = and(T_avg>=T_thresh, c_avg>=accum_thresh);

    FAC_RCM(isnan(FAC_obs)) = NaN;
    hold  on
    plot(FAC_RCM(ind_DSA),FAC_obs(ind_DSA),...
        'o', 'Color','k', 'MarkerFaceColor', RGB('jaune'),'LineWidth',0.5)
    plot(FAC_RCM(ind_LAWSA),FAC_obs(ind_LAWSA),...
        'o', 'Color', 'k', 'MarkerFaceColor', RGB('rouge'),'LineWidth',0.5)
    plot(FAC_RCM(ind_HAWSA),FAC_obs(ind_HAWSA),...
        'o', 'Color', 'k', 'MarkerFaceColor', RGB('vert'),'LineWidth',0.5)
    [x, ind_sorted] = sort(FAC_RCM);
    y = FAC_obs(ind_sorted);
%     [lm, ah] = Plotlm(x, y,'Annotation','off','LineStyle','-','LineWidth',2,...
%         'Color','r');
lm = fitlm(x, y);
    axis tight
    temp1 =0;
    temp2 = 35;
    h_1 = plot([temp1 temp2], [temp1 temp2],'k','LineWidth',2);
    h_p = patch([temp1  temp1 temp2 temp2 ], ...
        [temp1-unc_FAC temp1+unc_FAC temp2+unc_FAC temp2-unc_FAC],...
        RGB('light light gray'));
    h_p.LineStyle = 'none';
    uistack(h_1,'bottom')
    uistack(h_p,'bottom')

    box on
    ME_DSA = nanmean(FAC_RCM(ind_DSA)-FAC_obs(ind_DSA));
    RMSE_DSA = sqrt(nanmean((FAC_RCM(ind_DSA)-FAC_obs(ind_DSA)).^2));
    ME_LAWSA = nanmean(FAC_RCM(ind_LAWSA)-FAC_obs(ind_LAWSA));
    RMSE_LAWSA = sqrt(nanmean((FAC_RCM(ind_LAWSA)-FAC_obs(ind_LAWSA)).^2));
    ME_HAWSA = nanmean(FAC_RCM(ind_HAWSA)-FAC_obs(ind_HAWSA));
    RMSE_HAWSA = sqrt(nanmean((FAC_RCM(ind_HAWSA)-FAC_obs(ind_HAWSA)).^2));
    ME_all = nanmean(FAC_RCM-FAC_obs);
    RMSE_all = sqrt(nanmean((FAC_RCM-FAC_obs).^2));

    out = [ME_DSA RMSE_DSA ME_LAWSA RMSE_LAWSA ME_HAWSA RMSE_HAWSA ...
        ME_all RMSE_all lm.Coefficients.Estimate(1:2)'];
    
    set(gca,'TickLength',[0.02 0.05]*2,'Layer','top')

end