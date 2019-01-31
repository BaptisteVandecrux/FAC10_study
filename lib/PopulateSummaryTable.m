function [summary] = PopulateSummaryTable(map, uncert_map,...
    DSA, LAPA, HAPA, unit, R)

summary = array2table(NaN(4,5));
summary.Properties.VariableNames = {'DSA', 'HAPA', 'LAPA_pre_2010', 'LAPA_post_2010','Gr'};

temp = map{1};
temp(DSA ~= 1) = NaN;
temp2 = uncert_map{2};
temp2(DSA ~= 1) = NaN;
summary.DSA(1) = nanmean(temp(:));
summary.DSA(2) = nanmean(temp2(:));

temp = map{1};
temp(LAPA ~= 1) = NaN;
temp2 = uncert_map{1};
temp2(LAPA ~= 1) = NaN;
summary.LAPA_pre_2010(1) = nanmean(temp(:));
summary.LAPA_pre_2010(2) = nanmean(temp2(:));

temp = map{2};
temp(LAPA ~= 1) = NaN;
temp2 = uncert_map{2};
temp2(LAPA ~= 1) = NaN;
summary.LAPA_post_2010(1) = nanmean(temp(:));
summary.LAPA_post_2010(2) = nanmean(temp2(:));

temp = map{2};
temp(HAPA ~= 1) = NaN;
temp2 = uncert_map{2};
temp2(HAPA ~= 1) = NaN;
summary.HAPA(1) = nanmean(temp(:));
summary.HAPA(2) = nanmean(temp2(:));

fprintf('mean of the difference in LAPA : %0.2f\n\n', ...
nanmean(map{1}(LAPA == 1) - map{2}(LAPA == 1)));
fprintf('check (should be minus the above) : %0.2f\n\n',...
summary.LAPA_post_2010(1) - summary.LAPA_pre_2010(1))

temp = map{2};
temp2 = uncert_map{2};
summary.Gr(1) = nanmean(temp(:));
summary.Gr(2) = nanmean(temp2(:));

% Summing FAC10 each area
SpatialIntegration = @(M,R) ...
    (nansum(M* R.CellExtentInWorldX *R.CellExtentInWorldY))/10^9;

summary.Properties.RowNames{3} = 'Summed FAC10';
summary.Properties.RowNames{4} = 'Uncertainty on summed FAC10';
summary.DSA(3) = SpatialIntegration(map{2}(DSA == 1), R);
summary.DSA(4) = SpatialIntegration(uncert_map{2}(DSA == 1), R);
summary.LAPA_pre_2010(3) = SpatialIntegration(map{1}(LAPA == 1), R);
summary.LAPA_pre_2010(4) = SpatialIntegration(uncert_map{1}(LAPA == 1), R);
summary.LAPA_post_2010(3) = SpatialIntegration(map{2}(LAPA == 1), R);
summary.LAPA_post_2010(4) = SpatialIntegration(uncert_map{2}(LAPA == 1), R);
summary.HAPA(3) = SpatialIntegration(map{2}(HAPA == 1), R);
summary.HAPA(4) = SpatialIntegration(uncert_map{2}(HAPA == 1), R);

summary.Gr(3) = SpatialIntegration(map{2}(:), R);
summary.Gr(4) = SpatialIntegration(uncert_map{2}(~isnan(map{2})), R);

% TotalFAC = (SpatialIntegration(FAC10_map{2}(HAPA == 1), R) + SpatialIntegration(FAC10_map{2}(LAPA == 1), R) + SpatialIntegration(FAC10_map{1}(DSA == 1), R));
% TotalFAC_check = SpatialIntegration(FAC10_map{2}(~isnan(FAC10_map{2})), R);

loss_LAPA = SpatialIntegration(map{1}(LAPA == 1), R) - SpatialIntegration(map{2}(LAPA == 1), R);
uncert_loss_LAPA = SpatialIntegration(uncert_map{1}(LAPA == 1), R) + SpatialIntegration(uncert_map{2}(LAPA == 1), R);
fprintf('before %0.2f\nafter %0.2f\n\n', ...
    summary.LAPA_post_2010(3), summary.LAPA_pre_2010(3));
fprintf('Loss LAPA (%s): %0.2f +/- %0.2f \n\n', unit, loss_LAPA, uncert_loss_LAPA);

end