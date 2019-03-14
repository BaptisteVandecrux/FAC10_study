%% ============== Mapping of the firn air content ========================
% This script calculates the spatial distribution of firn air content down
% to 10 m depth on the Greenland ice sheet. More details are available in:
%
% Vandecrux, B., M. MacFerrin, H. Machguth, W.T. Colgan, D. van As, A. Heilig,
% C. M. Stevens, C. Charalampidis, R. Fausto, E. M. Morris, E. Mosley-Thompson,
% L. Koenig, L.N. Montgomery, C. Miège, S. Simonsen, T. Ingeman-Nielsen,
% J.E. Box, Brief communication: New firn dataset shows the evolution of
% the firn air content on the Greenland ice sheet. Submitted to the Cryosphere.
%
% It requires:
% - all functions located in the lib folder
% - FAC10 dataset available for download on www.promice.dk
% - Firn area delineation available for download on www.promice.dk
% - Temperature and accumulation maps from Box 2013 and Box et al. 2013
%         Citation:
%         1.	Box, J. E. 2013. Greenland ice sheet mass balance reconstruction. Part II: Surface mass balance (1840-2010), Journal of Climate,Vol. 26, No. 18. 6974-6989.  doi:10.1175/JCLI-D-12-00518.1
%         2.	Box, J.E., N. Cressie, D.H. Bromwich, J. Jung, M. van den Broeke, J.H. van Angelen, R.R. Forster, C. Miège, E. Mosley-Thompson, B. Vinther, J.R. McConnell. 2013. Greenland ice sheet mass balance reconstruction. Part I: net snow accumulation (1600-2009). Journal of Climate, 26, 3919–3934. doi:10.1175/JCLI-D-12-00373.1
% - Temperature and accumulation maps Fettweis et al. 2017
%         Citation:
%         Fettweis, X., Box, J., Agosta, C., Amory, C., C., K., Lang, C., . . . Gallée, H.: Reconstructions of the 1900–2015 Greenland ice sheet surface mass balance using the regional climate MAR model, Cryosphere, 11, 2, 1015-1033, doi:doi:10.5194/tc-11-1015-2017, 2017.
%
% Baptiste Vandecrux
% Technical University of Denmark
% Geological Survey of Greenland and Denmark
% b.vandecrux@gmail.com
% 23-08-2018
% =========================================================================

% Graph options and path setting
clear all 
close all
clc
set(0,'defaultfigurepaperunits','centimeters');
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxesTitleFontSize',1)
set(0,'DefaultLegendFontSize',15)
set(0,'defaultfigurecolor','w');
set(0,'defaultfigureinverthardcopy','off');
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'defaultfigurepaperpositionmode','auto');
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultAxesFontName', 'Times New Roman')
   
warning('off','curvefit:fit:noStartPoint');
warning('off','MATLAB:table:RowsAddedExistingVars');
warning('off','MATLAB:table:RowsAddedNewVars');

vis = 'on'; % make plot visible. If 'off' plots are not displayed but still printed
mkdir('.\Output')
addpath(genpath('.\lib'))
addpath(genpath('.\Input'))

%% ========= Here you choose the T and b map source ====================
% 1 is Box13 2 is MAR
 source_temp_accum = 2;
switch source_temp_accum 
    case 1
         accum_thresh = 600;
         T_thresh = -16;
         orig = 'Box13';
     case 2 
         accum_thresh = 600;
         T_thresh = -19;
         orig = 'MAR';
end

filename= sprintf('Output/result_%s.txt', orig);
diary(filename)

% Loading temperature and accumulation data
disp('Loading long-term temperature and accumulation maps')
switch source_temp_accum
    case 1
            c_map = table;
            T_map = table;
        if exist('mean_temperature_box13.csv') == 2
            temp = dlmread('mean_temperature_box13.csv',';');
            T_map.lon = temp(:,1);
            T_map.lat = temp(:,2);
            T_map.T_avg = temp(:,3);
            
            temp = dlmread('mean_accumulation_box13.csv',';');
           	c_map.lon = temp(:,1);
            c_map.lat = temp(:,2);
            c_map.c_avg = temp(:,3);
        else
            % Load accumulation from Box 13
            namefile = 'Box_Greenland_Accumulation_annual_1840-1999_ver20140214.nc';
            finfo = ncinfo(namefile);
            names={finfo.Variables.Name};
            for i= 1:size(finfo.Variables,2)
                eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
            end
            fprintf('\nData extracted from nc files.\n');

            % Load temprature from Box 13
            namefile = 'Box_Greenland_Temperature_monthly_1840-2014_5km_cal_ver20141007.nc';
            finfo = ncinfo(namefile);
            names={finfo.Variables.Name};
            for i= 1:size(finfo.Variables,2)
                eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
            end
            fprintf('\nData extracted from nc files.\n');

            c_map.lon = lon(:);
            T_map.lat = lat(:);
            c_map.lon = lon(:);
            T_map.lat = lat(:);

            c_map.c_avg = mean (acc(:,:,131:end),3);
            T_map.T_avg = mean (Temperature(:,:,131:end),3);
            
            dlmwrite('./Output/Temp accum maps/mean_temperature_box13.csv', [lon(:) lat(:) T_map.T_avg(:)],'Delimiter',';');
            dlmwrite('./Output/Temp accum maps/mean_accumulation_box13.csv', [lon(:) lat(:) c_map.c_avg(:)],'Delimiter',';');
        end
        
    case 2  % MAR
        % temperature
        filename = 'T_avg_1979-2014_MAR.csv';
        delimiter = ',';
        formatSpec = '%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
        fclose(fileID);
        T_map = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','T_avg'});
        clearvars filename delimiter formatSpec fileID dataArray ans;
        
        % snowfall
        filename = 'Net_Snowfall_avg_1979-2014_MAR.csv';
        delimiter = ',';
        formatSpec = '%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
        fclose(fileID);
        c_map = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','c_avg'});
        clearvars filename delimiter formatSpec fileID dataArray ans;
end
        c_map.c_avg(c_map.c_avg==-999) = NaN;
        T_map.T_avg(T_map.T_avg==-999) = NaN;
        lon =  reshape(T_map.lon,561,301);
        lat =  reshape(T_map.lat,561,301);
        
%% Loading FAC10 dataset
disp('Loading FAC10 dataset')

if exist('./Input/Core_all2.mat') == 2
    disp('Calculating FAC10 from firn density')
    [metadata] = Create_FAC10_dataset(T_map,c_map,vis);
else
    disp('Reading FAC10 dataset in file')
    filename = '.\Input\FAC_dataset.csv';
    [metadata] = LoadFAC10Dataset(filename);
end

% The metadata is distributed with average temperature and accumulation
% from MAR.
% If working with  Box (2013) and Box et al. (2013), 
% then they need to be overwritten.
if source_temp_accum == 1
    for i = 1:height(metadata)        
        % accumulation
         [~, ind] = min(distance( metadata.Latitude(i),metadata.Longitude(i),...
             c_map.lat,c_map.lon));
        metadata.c_avg(i) = c_map.c_avg(ind);
        metadata.T_avg(i) = T_map.T_avg(ind);
    end
end

% here we get the coordinates in the Polar Stereographic projection
I = geotiffinfo('mean_temperature_box13_3413_firn.tif');
[metadata.X, metadata.Y] = projfwd(I, metadata.Latitude,metadata.Longitude);
metadata_save = metadata;

%% Defining region to be considered
disp('Defining region to be considered')

% Loading firn grid
[FirnArea, R_firn]= geotiffread('FirnLayer_2000_2017_final_4326.tif');
I = geotiffinfo('FirnLayer_2000_2017_final_4326.tif');
FirnArea(isnan(FirnArea)) = 0;
FirnArea(FirnArea>1) = 0;
[x_firn,y_firn]=pixcenters(I);
[XX_firn,YY_firn] = meshgrid(x_firn,y_firn);

is_firn = interp2(XX_firn,YY_firn,FirnArea,lon(:), lat(:),'nearest');
is_firn(isnan(is_firn)) = 0;

T_map.T_avg(find(~is_firn)) = NaN;
c_map.c_avg(find(~is_firn)) = NaN;

% Loading T and b raster files
Ice = shaperead('IcePolygon_3413.shp');
GL = shaperead('GL_land.shp');
[FirnArea_3413, R_fa3413] = geotiffread('mean_temperature_box13_3413_firn.tif');
Firn = shaperead('FirnLayer2000-2017_final_3413.shp');

switch source_temp_accum
    case 1
        [T_raster, ~] = geotiffread('mean_temperature_box13_3413_firn.tif');
        [c_raster, R] = geotiffread('mean_accumulation_box13_3413_firn.tif');
        I = geotiffinfo('mean_temperature_box13_3413_firn.tif'); 

    case 2
        [T_raster, ~] = geotiffread('T_avg_1979-2014_MAR_3413_firn.tif');
        [c_raster, R] = geotiffread('Net_Snowfall_avg_1979-2014_MAR_3413_firn.tif');
        I = geotiffinfo('T_avg_1979-2014_MAR_3413_firn.tif'); 
end
        [x_map,y_map]=pixcenters(I);
        XX = repmat(x_map,length(y_map),1);
        YY = repmat(y_map',1,length(x_map));
        
        T_raster(T_raster==0)=NaN;
        T_raster(T_raster==-999)=NaN;
        T_raster(T_raster==-9999)=NaN;
        T_raster(T_raster<-9999)=NaN;
        c_raster(c_raster==0)=NaN;
        c_raster(c_raster==-999)=NaN;
        c_raster(c_raster==-9999)=NaN;
        c_raster(c_raster<-9999)=NaN;
      
        % Cold area
        DSA = (T_raster<=T_thresh)+0;
        % Warm High Accumulation area
        HAPA = and(T_raster>T_thresh,c_raster>accum_thresh)+0;
        % Warm Low Accumulation area
        LAPA = and(T_raster>T_thresh,c_raster<=accum_thresh)+0;
        
        format BANK 
disp('Firn area (km2):')
fprintf('%0.2f\n\n',sum(sum(~isnan(T_raster))) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);

disp('DSA (km2):')
fprintf('%0.2f\n',sum(sum(DSA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %% \n',sum(sum(DSA)) /sum(sum(~isnan(T_raster)))*100 );
fprintf('%i observations \n\n',length(metadata.FAC10(metadata.T_avg<T_thresh)) );

disp('LAPA (km2):')
fprintf('%0.2f\n',sum(sum(LAPA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %%\n',sum(sum(LAPA)) /sum(sum(~isnan(T_raster)))*100 );
fprintf('%i observations after 2010\n',length(metadata.FAC10(...
    and(and(metadata.T_avg>=T_thresh, metadata.c_avg<accum_thresh), ...
    metadata.Year>=2010) ) ) );
fprintf('%i observations between 1997 and 2010\n\n',length(metadata.FAC10(...
    and(and(metadata.T_avg>=T_thresh, metadata.c_avg<accum_thresh), ...
    and(metadata.Year<=2009, metadata.Year>=1997) ) ) ) );

metadata(and(and(metadata.T_avg>=T_thresh, metadata.c_avg<accum_thresh), ...
    and(metadata.Year<=2009, metadata.Year>=1997) ),:);

disp('HAPA (km2):')
fprintf('%0.2f\n',sum(sum(HAPA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %%\n',sum(sum(HAPA)) /sum(sum(~isnan(T_raster)))*100 );
fprintf('%i observations after 2010\n',length(metadata.FAC10(...
    and(and(metadata.T_avg>=T_thresh, metadata.c_avg>=accum_thresh), ...
    metadata.Year>=2010) ) ) );
fprintf('%i observations from 1998\n\n',length(metadata.FAC10(...
    and(and(metadata.T_avg>=T_thresh, metadata.c_avg>=accum_thresh), ...
    and(metadata.Year<=1998, metadata.Year>=1998) ) ) ) );

    x = T_map.T_avg(and(~isnan(T_map.T_avg),~isnan(c_map.c_avg)) );
    y = c_map.c_avg(and(~isnan(T_map.T_avg),~isnan(c_map.c_avg)) );
    if source_temp_accum==2
        k = boundary(x, y, 0.5);
    else
        k = boundary(x, y, 0.2);
    end
    poly_100 = [x(k), y(k)];

PlottingFigure1(metadata, T_thresh, accum_thresh, ...
    GL,Ice,Firn, DSA, HAPA, LAPA, poly_100, I, orig, vis);

%% FACtot
metadata = metadata_save;

Harper_latlon = [69.87650	47.01020; ...
69.84802	47.27358; ...
69.81998	47.45050; ...
69.78360	47.67018; ...
69.75693	47.88028; ...
69.73802	48.06097; ...
69.72505	48.19020; ...
69.73908	48.24030; ...
69.71978	48.26740; ...
69.70617	48.34497; ...
69.68743	48.49967; ...
69.67393	48.59112; ...
69.66018	48.68945];

FC_Harper = [3717		9798; ...
3642		10891; ...
3337		8580; ...
2931		10676; ...
3133		10280; ...
2758		12071; ...
2529		7340; ...
2376		6145; ...
2167		5125; ...
1148		2915; ...
1902		2619; ...
1816		3287; ...
517		540];

FAC_Harper = FC_Harper*NaN;
FAC_Harper(:,1) = (10*917 - 10 * 843 + FC_Harper(:,1)) ./ 917;
FAC_Harper(:,2) = (100*917 - 100 * 843 + FC_Harper(:,2)) ./ 917;

f= figure('Visible',vis,'Outerposition',[1 1 13 13]);
ha = tight_subplot(2,1,0.1,[0.3 0],[0.1 0.3]);
    set(gcf,'CurrentAxes', ha(1))
    
    ind_long = find(metadata.DepthMax>=100);
    x = metadata.FAC10(ind_long);
    y = metadata.FACtot(ind_long);
    col = brewermap(length(ind_long)+1,'set1')*0+1;
   
    symbol = {'o'}; %'o','s','d','^','v','>','<'};
    hold on
    for i = 1:length(ind_long)
        plot(metadata.FAC10(ind_long(i)),...
            metadata.FACtot(ind_long(i)),...
            symbol{mod(i,length(symbol))+1},'MarkerFaceColor',col(i,:),...
            'MarkerEdgeColor',[0 0 0],'LineWidth',2)
    end
    i=i+1;
    plot(FAC_Harper(:,1), FAC_Harper(:,2),symbol{mod(i,length(symbol))+1},...
        'MArkerFaceColor',col(i,:),'MArkerEdgeColor',[0 0 0],'LineWidth',2)

        x1 = [x; FAC_Harper(:,1)];
        y1 = [y; FAC_Harper(:,2)];
        [x1,ind] = sort(x1);
        y1=y1(ind);
        lm = fitlm(x1,y1,'Intercept',false);
        h1 = plot(x1, lm.Coefficients.Estimate(1)*x1,'k');
        h2 = patch(x1([1 end end 1]), ...
            [lm.Coefficients.Estimate(1)*x1([1 end])+2*lm.RMSE ;...
            (lm.Coefficients.Estimate(1)*x1([end 1])-2*lm.RMSE)],...
            RGB('light light gray'));
        h2.LineStyle = 'none';
        uistack(h2,'bottom')
unc_FACtot = 2*lm.RMSE;
    box on
    axis tight square
    xlabel('Observed FAC_{10} (m)','Interpreter','tex')
    ylabel('Observed FAC_{tot} (m)','Interpreter','tex')
legend([h1 h2],'Linear regression','2 \times RMSD',...
    'Location','NorthWest','Interpreter','tex')
    set(gca,'Ticklength',[0.08 0.16]/2,'layer','top')

set(gcf,'CurrentAxes', ha(2))
hold on
    [p4, p5] = PlotBackground(GL,Ice,Firn);

for i = 1:length(ind_long)

    plot(metadata.X(ind_long(i)),...
        metadata.Y(ind_long(i)),...
        symbol{mod(i,length(symbol))+1},...
    'MArkerFaceColor',col(i,:),'MArkerEdgeColor',[0 0 0],'LineWidth',2)
%     leg_text{i} = Core{metadata.CoreNumber(ind_long(i))}.Info.Name;
end
i=i+1;

[X_Harper, Y_Harper] = polarstereo_fwd( Harper_latlon(:,1),-Harper_latlon(:,2),...
    6378137,0.08181919, 70,-45);
plot(X_Harper, Y_Harper,...
    symbol{mod(i,length(symbol))+1},...
    'MArkerFaceColor',col(i,:),'MArkerEdgeColor',[0 0 0],'LineWidth',2)
    daspect( [1 1 1])
ha(2).XTickLabel = '';
ha(2).YTickLabel = '';
ha(2).Box = 'on';
ha(1).Units = 'normalized';
ha(1).Position = [0.1 0.17 0.8 0.75];

ha(2).Units = 'normalized';
ha(2).Position = [0.59 0.21 0.2 0.35];
    axis fill
    set(gca,'Ticklength',[0.08 0.16]*0,'layer','top')

print(f,'Output/FAC10_vs_FACtot','-dtiff','-r350') 
print(f,'Output/FAC10_vs_FACtot','-dpdf') 

% Assigning FACtot
metadata.FACtot = lm.Coefficients.Estimate(1) * metadata.FAC10;
metadata.FACtot(ind_long) = y;
a_FAC = lm.Coefficients.Estimate(1);
% b_FAC = lm.Coefficients.Estimate(1);
% b_FAC = 0;
metadata_save = metadata;

%% ===== Temperature regression in the Dry Snow Area ===============
disp('Selecting points in the cold area')
vis_save = vis;
metadata = metadata_save;
met = metadata;
ind = and(met.T_avg>=T_thresh,met.c_avg<=accum_thresh);
% ind = and(met.T_avg>T_thresh,met.c_avg>=accum_thresh);
met(ind,:) = [];
fprintf('\nUsing %i cores to build the DSA\n',size(met,1))
% 3D exploration of accumulation and temperature control

% Fitting DSA surface
    T_bins = linspace(min(met.T_avg),max(met.T_avg),5);
    ind_bins = discretize(met.T_avg,T_bins);
    W = NaN(size(met.T_avg));
    for i =1:length(met.T_avg)
        W(i) = 1/sum(ind_bins==ind_bins(i)) *1/(length(T_bins)-1);
    end

    FAC_fit_temp = fit(met.T_avg,met.FAC10,'poly1',...
        'Weight', W);
    Coef = [FAC_fit_temp.p1 FAC_fit_temp.p2];
    FAC_fit_DSA = @(x,y) reshape(Coef(1).*y(:)+Coef(2),size(y));

[RMSD, bias] = PlottingDSA(met, FAC_fit_DSA, Coef, T_thresh, vis, orig);
disp(RMSD)
uncert_DSA = RMSD*2;

%% climatic control on firn line
%loading firn line
Firn_4326 = shaperead('FirnLayer2000-2017_final_4326.shp');

ind_nonan = find(~isnan(T_map.T_avg));
DT = delaunayn([T_map.lon(ind_nonan), T_map.lat(ind_nonan)]);
ind_firn_line = [];

for i = 1:length(Firn_4326)
    A = polyarea(Firn_4326(i).X(~isnan(Firn_4326(i).X))', ...
        Firn_4326(i).Y(~isnan(Firn_4326(i).X))');

        ind = dsearchn([T_map.lon(ind_nonan), T_map.lat(ind_nonan)], DT, ...
        [Firn_4326(i).X(~isnan(Firn_4326(i).X))' Firn_4326(i).Y(~isnan(Firn_4326(i).X))']);
        ind_firn_line = [ind_firn_line; unique(ind_nonan(ind))];
end

clearvars DT
ind_no_nan=find(~isnan(T_map.T_avg));

%% Plotting
% fitting a polynomial to the firn line

lat_FL_obs = T_map.lat(ind_firn_line);
lon_FL_obs = c_map.lon(ind_firn_line);
[X_FL_obs, Y_FL_obs] = polarstereo_fwd(lat_FL_obs, lon_FL_obs,...
    6378137,0.08181919, 70,-45);

T_FL_obs = T_map.T_avg(ind_firn_line);
c_FL_obs = c_map.c_avg(ind_firn_line);
% ind1 = [1305; 1335; 1271; 1165; 1265; 1425];
ind1 = 1000:5:1400;
if source_temp_accum == 2
    ind1 = [1110; 1240;1330; 1050;1035; 1195];
    c_max = 3000;
else
    ind1 = [1320; 1225; 1365; 1390;1180; 1235];
    c_max = 5000;
end
ind1=ind1';

PlottingSelectedFirnLinePoints(metadata, T_thresh, accum_thresh, ind1,...
    poly_100, T_FL_obs, c_FL_obs,XX,YY,DSA,LAPA,HAPA,X_FL_obs,Y_FL_obs, ...
    c_max, source_temp_accum, GL,Ice,Firn, orig, vis);

T_FL_selec = T_FL_obs(ind1);
c_FL_selec = c_FL_obs(ind1);

%% ============= Empirical function for all regions ==================
disp('Outlook for all regions')
metadata = metadata_save;

FAC_nofirn = 0;

% 3D exploration of accumulation and temperature control
    % Defining the points over which FAC will be extrapolated
    c_extrap = linspace(0,5000,500);
    T_extrap = linspace(-50, -2, 500);
   
    y_all = [metadata.c_avg; c_FL_selec];
    x_all = [metadata.T_avg; T_FL_selec];
    z_all = [metadata.FAC10; FAC_nofirn*ones(size(T_FL_selec))];
    year_all = [metadata.Year; 2017*ones(size(T_FL_selec))];
    
    smoothness = 1;

    % Fitting function for the 1998-2008 period
    % here we exclude all the measurements after 2010 and before 1998
    cond1 = ((x_all > T_thresh) + (z_all~=FAC_nofirn) + or(year_all>=2010,year_all<1998)) ~= 3;
    x1 = x_all(cond1);
    y1 = y_all(cond1);
    z1 = z_all(cond1);

    %     year1 = year_all(cond1);
    [FAC_fit_all_pre2010, ~, ~] = gridfit(x1, y1, z1, T_extrap, c_extrap,'smoothness',smoothness);
        
    % here we exclude all the measurements in the LAPA before 2010
    cond2 = ((y_all < accum_thresh) +...
        (x_all>T_thresh) + ...
        (z_all~=FAC_nofirn) + (year_all<2010)) ~= 4;
    cond2(and(year_all == 1998, y_all > accum_thresh)) = 0;
    x2 = x_all(cond2);
    y2 = y_all(cond2);
    z2 = z_all(cond2);
%     year2 = year_all(cond1);
    [FAC_fit_all_post2010, gridx, gridy] = gridfit(x2, y2,z2, T_extrap, c_extrap,'smoothness',smoothness);

    % Values of FAC derived from our empirical function for DSA
    FAC_DSA_grid = FAC_fit_DSA(gridx,gridx);
    % FAC_BV = FAC_BV_func(bb,TT,280,1,-5);

ind_in_firn = reshape(inpoly([gridx(:) gridy(:)],(poly_100)),size(gridx));
FAC_fit_all_pre2010 = min(FAC_fit_all_pre2010, FAC_DSA_grid);
FAC_fit_all_pre2010(gridx<T_thresh) = FAC_DSA_grid(gridx<T_thresh);
FAC_fit_all_pre2010(FAC_fit_all_pre2010<0) = 0;

FAC_fit_all_post2010 = min(FAC_fit_all_post2010,FAC_DSA_grid);
FAC_fit_all_post2010(gridx<T_thresh) = FAC_DSA_grid(gridx<T_thresh);
FAC_fit_all_post2010(FAC_fit_all_post2010<0) = 0;

% creating the function
FAC10_98_08_func = @(x,y) interp2(gridx,gridy,FAC_fit_all_pre2010,x,y);
FAC10_10_17_func = @(x,y) interp2(gridx,gridy,FAC_fit_all_post2010,x,y);

% 3D plotting
PlottingAllRegions_3D(metadata,gridx, gridy, FAC_fit_all_post2010,...
    FAC_fit_all_pre2010, x1, y1, z1, x2, y2, z2, ind_in_firn, vis,orig)

% 2d plots
PlottingAllRegions_2D(metadata,gridx, gridy, FAC10_98_08_func,...
    FAC10_10_17_func, poly_100, T_thresh,accum_thresh, x1, y1, z1, x2, y2, z2, vis,orig)

%%  Applying on the T_avg map 

c_raster(isnan(T_raster)) = NaN;
c_raster(T_raster<-60) = NaN;
T_raster(T_raster<-60) = NaN;
   
FAC10_map = {};
h_CA_bis= {};
h_WHA_bis= {};
h_WLA_bis= {};

    % pre 2010
    FAC10_map{1} = max(zeros(size(T_raster)), FAC10_98_08_func(T_raster, c_raster));
    %post 2010
    FAC10_map{2} = max(zeros(size(T_raster)), FAC10_10_17_func(T_raster, c_raster));
%     FAC_10_map{2,2} = max(zeros(size(T_raster)), FAC_bav_spread(c_raster,T_raster));
    
    FAC10_map{1}(isnan(T_raster))=NaN;
    FAC10_map{2}(isnan(T_raster))=NaN;
    FAC10_map{1}(and(T_raster>T_thresh, c_raster>accum_thresh))=NaN;   
    
    ind_met_DSA = metadata.T_avg<T_thresh;
    ind_met_LAPA_post_2010 = and( and(metadata.T_avg>T_thresh, ...
        metadata.c_avg<accum_thresh),...
        metadata.Year>=2010);
    ind_met_LAPA_pre_2010 = and(and( and(metadata.T_avg>T_thresh, ...
        metadata.c_avg<accum_thresh),...
        metadata.Year<2010),...
        metadata.Year>=1998);
    ind_met_HAPA_post_2010 = and( and(metadata.T_avg>T_thresh, ...
        metadata.c_avg>=accum_thresh),...
        metadata.Year>=2010);
    
    metadata.FAC10_predict = NaN * metadata.FAC10;
    metadata.FAC10_predict(ind_met_DSA) = ...
        FAC10_98_08_func(metadata.T_avg(ind_met_DSA),...
        metadata.c_avg(ind_met_DSA));
    metadata.FAC10_predict(ind_met_LAPA_post_2010) = ...
        FAC10_10_17_func(metadata.T_avg(ind_met_LAPA_post_2010),...
        metadata.c_avg(ind_met_LAPA_post_2010));
    metadata.FAC10_predict(ind_met_HAPA_post_2010) = ...
        FAC10_10_17_func(metadata.T_avg(ind_met_HAPA_post_2010), ...
        metadata.c_avg(ind_met_HAPA_post_2010));
    metadata.FAC10_predict(ind_met_LAPA_pre_2010) =  ...
        FAC10_98_08_func(metadata.T_avg(ind_met_LAPA_pre_2010), ...
        metadata.c_avg(ind_met_LAPA_pre_2010));
    
    metadata.residual = metadata.FAC10_predict-metadata.FAC10;
    bias(1) = mean(metadata.residual(ind_met_DSA));
    RMSD(1) = sqrt(mean(metadata.residual(ind_met_DSA).^2));
    bias(2) = mean(metadata.residual(ind_met_LAPA_pre_2010));
    RMSD(2) = sqrt(mean(metadata.residual(ind_met_LAPA_pre_2010).^2));
    bias(3) = mean(metadata.residual(ind_met_LAPA_post_2010));
    RMSD(3) = sqrt(mean(metadata.residual(ind_met_LAPA_post_2010).^2));
    bias(4) = mean(metadata.residual(ind_met_HAPA_post_2010));
    RMSD(4) = sqrt(mean(metadata.residual(ind_met_HAPA_post_2010).^2));
    N(1) = sum(ind_met_DSA);
    N(2) = sum(ind_met_LAPA_pre_2010);
    N(3) = sum(ind_met_LAPA_post_2010);
    N(4) = sum(ind_met_HAPA_post_2010);

    summary_residuals = array2table([bias; RMSD; N]);
    summary_residuals.Properties.VariableNames = ...
        {'DSA', 'LAPA_pre_2010','LAPA_post_2010', 'HAPA'};
    summary_residuals.Properties.RowNames = {'bias', 'RMSD','N'};
      writetable(summary_residuals,...
    sprintf('./Output/Summary_residual_%s.csv',orig),'Delimiter',';','WriteRowNames',true)

        metadata_save = metadata;
        I_box = geotiffinfo('mean_temperature_box13_3413_firn.tif'); 
        [x_map,y_map]=pixcenters(I_box);
        XX_box = repmat(x_map,length(y_map),1);
        YY_box = repmat(y_map',1,length(x_map));
        
        I_MAR = geotiffinfo('T_avg_1979-2014_MAR_3413_firn.tif'); 
        [x_map,y_map]=pixcenters(I_MAR);
        XX_mar = repmat(x_map,length(y_map),1);
        YY_mar = repmat(y_map',1,length(x_map));
                
if source_temp_accum==1
    XX = XX_box;
    YY = YY_box;
else
    XX = XX_mar;
    YY = YY_mar;
end
format BANK

% Plotting FAC10 maps
 [FAC10_change] = Plotting_FAC10_maps (metadata, accum_thresh, T_thresh, ...
    XX,YY, FAC10_map, DSA,LAPA,HAPA,GL,Ice,Firn, orig, vis);

%% Uncertainty analysis 

disp('Sensitivity analysis')
close all

plotting = 0;
F = {};
SA_post2010_func = {};
SA_pre2010_func = {};

c_extrap = linspace(0,5000,100);
T_extrap = linspace(-50, -2, 100);

ind_LAPA_pre2010 =  find(((metadata.c_avg < accum_thresh) + ...
    (metadata.T_avg > T_thresh) + (metadata.Year<2010)) == 3);

ind_LAPA_post2010 =  find(((metadata.c_avg < accum_thresh) + ...
    (metadata.T_avg > T_thresh) + (metadata.Year>=2010)) == 3);
    
for ii = 1:1000
%         fprintf('%i %%\n',ii)

    if mod(ii,100) == 0
        fprintf('%i %%\n',ii/10)
    end
    
    met = metadata;

    % Factors to be modified
    FAC_nofirn = rand(1);
    smoothness = rand(1)+0.5;
    perturb_FL = normrnd(0,1);  % Simulate heights;    
    
    % We remove 4 from the LAPA in each time period
    ind_rand1 = randi([1 length(ind_LAPA_pre2010)],1,4)';
    ind_rand2 = randi([1 length(ind_LAPA_post2010)],1,4)';
    met([ind_LAPA_pre2010(ind_rand1); ind_LAPA_post2010(ind_rand2)],:) = [];

    % calculating points to be fitted  
    y_all = [met.c_avg; c_FL_selec];
    x_all = [met.T_avg; T_FL_selec + perturb_FL];
    z_all = [met.FAC10; FAC_nofirn*ones(size(T_FL_selec))];
    year_all = [met.Year; 2017*ones(size(T_FL_selec))];

    cond1 = ((y_all < accum_thresh) + (x_all > T_thresh) + (z_all~=FAC_nofirn) + (year_all>=2010)) ~= 4;
    x1 = x_all(cond1);
    y1 = y_all(cond1);
    z1 = z_all(cond1);
    cond2 = ((y_all < accum_thresh) + (x_all>T_thresh) + (z_all~=FAC_nofirn) + (year_all<2010)) ~= 4;
    x2 = x_all(cond2);
    y2 = y_all(cond2);
    z2 = z_all(cond2);
    
    % Defining the points over which FAC will be extrapolated

    %     year1 = year_all(cond1);
    [SA_pre2010, ~, ~] = gridfit(x1, y1, z1, T_extrap, c_extrap,'smoothness',smoothness);

    %     year2 = year_all(cond1);
    [SA_post2010, gridx, gridy] = gridfit(x2, y2,z2, T_extrap, c_extrap,'smoothness',smoothness);
    FAC_DSA_grid = FAC_fit_DSA(gridx,gridx);
    ind_in_firn = reshape(inpoly([gridx(:) gridy(:)],(poly_100)),size(gridx));

    SA_pre2010 = min(SA_pre2010,FAC_DSA_grid);
    SA_pre2010(gridx<T_thresh) = FAC_DSA_grid(gridx<T_thresh);
    SA_pre2010(SA_pre2010<0) = 0;

    SA_post2010 = min(SA_post2010,FAC_DSA_grid);
    SA_post2010(gridx<T_thresh) = FAC_DSA_grid(gridx<T_thresh);
    SA_post2010(SA_post2010<0) = 0;

   if plotting == 1
       PlottingAllRegions_3D(met,gridx, gridy, SA_post2010,...
            SA_pre2010,  x1, y1, z1, x2, y2, z2, ind_in_firn, vis,orig);
 
      F{ii} = getframe(gcf) ;
      drawnow
   end
   
    SA_pre2010_func{ii} = @(x,y) interp2(gridx,gridy,SA_pre2010,x,y);
    SA_post2010_func{ii} = @(x,y) interp2(gridx,gridy,SA_post2010,x,y);
end

    % create the video writer with 1 fps
      writerObj = VideoWriter('Output/SensitivityAnalysis.avi');
      writerObj.FrameRate = 10;
      % set the seconds per image

    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % conRGB('vert') the image to a frame
        frame = F{i} ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
XX_LAPA = T_raster;
YY_LAPA = c_raster;

M = NaN([size(c_raster),length(SA_post2010_func)]);
for ii = 1:length(SA_post2010_func)
    M(:,:,ii) =SA_post2010_func{ii}(XX_LAPA, YY_LAPA);
end

M_mean1 = mean(M,3);
Uncertainty_post2010 = std(M,0,3);

M = NaN([size(c_raster),length(SA_pre2010_func )]);
for ii = 1:length(SA_pre2010_func)
    M(:,:,ii) =SA_pre2010_func{ii}(XX_LAPA, YY_LAPA);
end
M_mean = mean(M,3);
Uncertainty_pre2010 = std(M,0,3);

clearvars writerObj F met gridx gridy SA_pre2010 SA_post2010 x1 y1 z1 
clearvars x2 y2 z2 x_all y_all z_all

%% Plotting uncertainty
metadata = metadata_save;

figure('Visible',vis)
subplot(1,2,1)
surf(XX,YY,Uncertainty_pre2010)
shading interp
colorbar
title('Uncertainty in the LAPA (1997-2008)')

subplot(1,2,2)
surf(XX,YY,Uncertainty_post2010)
shading interp
colorbar
title('Uncertainty in the LAPA (2011-2017)')

% creating uncertainty map
delta_FAC10 = FAC10_map;
delta_FAC10{1}= NaN*FAC10_map{1};
delta_FAC10{2}= NaN*FAC10_map{1};

delta_FAC10{2} = max(0.3,Uncertainty_post2010*2);
delta_FAC10{1} = max(0.3, Uncertainty_pre2010*2);
    
delta_FAC100 = delta_FAC10;
% delta_FAC100{2} = max(unc_FACtot*ones(size(temp1{1})), a_FAC*delta_FAC10{2});
% delta_FAC100{1} = max(unc_FACtot*ones(size(temp1{1})), a_FAC*delta_FAC10{1});
delta_FAC100{2} = a_FAC*delta_FAC10{2};
delta_FAC100{1} = a_FAC*delta_FAC10{1};

disp('Mapping Uncertainty of pore space')

disp('Mean residual:')
fprintf('%0.2f m3 m-3\n',nanmean(metadata.residual))
    
PlottingUncertaintyMap(XX,YY, delta_FAC10, FAC10_map, ...
    DSA, HAPA, LAPA, GL, Ice, Firn, orig, vis);
% PlottingUncertaintyMap2(XX,YY, delta_FAC10, FAC10_map, ...
%     DSA, HAPA, LAPA, GL, Ice, Firn, orig, vis);

%% Firn capacity
% Caluclating the Firn Capacity (FC) as the amount of infiltration ice is
% needed to bring the upper 10 m of firn to a density of 843 kg/m3
% according to Harper et al. (2012)
FC10{1} = max(0, 10*843 - 917*(10-FAC10_map{1}));
FC10{2} = max(0, 10*843 - 917*(10-FAC10_map{2}));
FC10{1}(isnan(FAC10_map{1})) = NaN;
FC10{2}(isnan(FAC10_map{2})) = NaN;

FC100{1} = max(0, 100*843 - 917*(100-a_FAC*FAC10_map{1}));
FC100{2} = max(0, 100*843 - 917*(100-a_FAC*FAC10_map{2}));
FC100{1}(isnan(FAC10_map{1})) = NaN;
FC100{2}(isnan(FAC10_map{2})) = NaN;

[FC10_change] = Plotting_FC10_maps (metadata, accum_thresh, T_thresh, ...
    XX,YY, FC10, DSA,LAPA,HAPA,GL,Ice,Firn, orig, vis);

% FACpco_map{1} = FAC10_map{1} * a_FAC + b_FAC;
% FACpco_map{2} = FAC10_map{2} * a_FAC + b_FAC;
% FCpco{1} = max(0, 10*843 - 917*(10-FACpco_map{1}));
% FCpco{2} = max(0, 10*843 - 917*(10-FACpco_map{2}));
% FCpco{1}(isnan(FACpco_map{1})) = NaN;
% FCpco{2}(isnan(FACpco_map{2})) = NaN;

%% Spatial averages
% average values
clc
disp(' ============= Spatial FAC10 ============= ')
[summary_FAC10] = PopulateSummaryTable(FAC10_map, delta_FAC10,...
    DSA, LAPA, HAPA, 'km3',R);
summary_FAC10.Properties.RowNames{1} = 'Mean FAC10';
summary_FAC10.Properties.RowNames{2} = 'Uncertainty on mean FAC10';
summary_FAC10.Properties.RowNames{3} = 'Summed FAC10';
summary_FAC10.Properties.RowNames{4} = 'Uncertainty on summed FAC10';
disp(summary_FAC10)

% FACtot
disp(' ============= Spatial FAC100 ============= ')
summary_FAC100 = summary_FAC10;
summary_FAC100{1,:} = summary_FAC10{1,:} *a_FAC;
summary_FAC100{2,:} = summary_FAC10{2,:} *a_FAC;
summary_FAC100{3,:} = summary_FAC10{3,:} *a_FAC;
summary_FAC100{4,:} = summary_FAC10{4,:} *a_FAC;
summary_FAC100.Properties.RowNames{1} = 'Mean FAC100';
summary_FAC100.Properties.RowNames{2} = 'Uncertainty on mean FAC100';
summary_FAC100.Properties.RowNames{3} = 'Summed FAC100';
summary_FAC100.Properties.RowNames{4} = 'Uncertainty on summed FAC100';
disp(summary_FAC100)

disp(' ============= Spatial FC10 ============= ')
temp1 = FC10;
temp1{1} = temp1{1}/1000;
temp1{2} = temp1{2}/1000;
temp2 = delta_FAC10;
temp2{1} = 917*temp2{1}/1000;
temp2{2} = 917*temp2{2}/1000;

[summary_FC10] = PopulateSummaryTable(temp1, temp2,...
    DSA, LAPA, HAPA, 'GT',R);
summary_FC10.Properties.RowNames{1} = 'Mean FC10 (GT)';
summary_FC10.Properties.RowNames{2} = 'Uncertainty on mean FC10 (GT)';
summary_FC10.Properties.RowNames{3} = 'Summed FC10 (GT)';
summary_FC10.Properties.RowNames{4} = 'Uncertainty on summed FC10 (GT)';
temp = summary_FC10(1:2,:);
temp.Properties.RowNames{1} = 'Summed FC10 (mm sle)';
temp.Properties.RowNames{2} = 'Uncertainty FC10 (mm sle)';
summary_FC10 = [summary_FC10; temp];
summary_FC10{5,:} = summary_FC10{3,:}/ 361;
summary_FC10{6,:} = summary_FC10{4,:}/ 361;
disp(summary_FC10)

disp(' ============= Spatial FC100 ============= ')
temp1 = FC100;
temp1{1} = temp1{1}/1000;
temp1{2} = temp1{2}/1000;
temp2 = delta_FAC100;
temp2{1} = 917*temp2{1}/1000;
temp2{2} = 917*temp2{2}/1000;

[summary_FC100] = PopulateSummaryTable(temp1, temp2,...
    DSA, LAPA, HAPA, 'GT',R);
summary_FC100.Properties.RowNames{1} = 'Mean FC100';
summary_FC100.Properties.RowNames{2} = 'Uncertainty on mean FC100';
summary_FC100.Properties.RowNames{3} = 'Summed FC100';
summary_FC100.Properties.RowNames{4} = 'Uncertainty on summed FC100';
temp = summary_FC100(1:2,:);
temp.Properties.RowNames{1} = 'Summed FC100 (mm sle)';
temp.Properties.RowNames{2} = 'Uncertainty FC100 (mm sle)';
summary_FC100 = [summary_FC100; temp];
summary_FC100{5,:} = summary_FC100{3,:} / 361;
summary_FC100{6,:} = summary_FC100{4,:}/ 361;
disp(summary_FC100)

summary_all = [summary_FAC10; ...
    summary_FAC100; ...
    summary_FC10; ...
    summary_FC100];

writetable(summary_all,...
    sprintf('./Output/Summary_output_%s.csv',orig),'Delimiter',';','WriteRowNames',true)

%% Comparison with Harper et al. 2012
[T_sorted, ind_sorted] = sort(T_raster(:),1,'descend');
i_first_nonan = find(~isnan(T_sorted),1,'first');

area_percolation = 1.5 * 10^5; %km2

num_cell = round(area_percolation/(R.CellExtentInWorldX * R.CellExtentInWorldY*10^(-6)));

is_perc_area = zeros(size(T_sorted));
is_perc_area_unsorted = zeros(size(T_sorted));
is_perc_area(i_first_nonan:i_first_nonan+num_cell) = 1;
is_perc_area_unsorted(ind_sorted) = is_perc_area;

is_perc_area_mat = reshape(is_perc_area_unsorted,size(T_raster));
is_perc_area_mat(is_perc_area_mat==0) = 0;
PA_LAPA = and(is_perc_area_mat==1,LAPA == 1);
PA_HAPA = and(is_perc_area_mat==1,HAPA == 1);

disp (' =========== FC10 in Harper''s percolation area =====')
SpatialIntegration = @(M,R) ...
    (nansum(M* R.CellExtentInWorldX *R.CellExtentInWorldY))/10^9;
FC_perc_area = 10^(-3)*SpatialIntegration(FC10{2}(is_perc_area_mat==1),R);
uncert_FC_perc_area = 10^(-3)*SpatialIntegration(917*delta_FAC10{2}(is_perc_area_mat==1), R);
fprintf('FC10 in percolation area: \n%0.2f +/- %0.2f\n', FC_perc_area, uncert_FC_perc_area)

FC_perc_area = 10^(-3)*SpatialIntegration(FC100{2}(is_perc_area_mat==1),R);
uncert_FC_perc_area = 10^(-3)*SpatialIntegration(917*delta_FAC100{2}(is_perc_area_mat==1), R);
fprintf('FC100 in percolation area: \n%0.2f +/- %0.2f\n', FC_perc_area, uncert_FC_perc_area)

FC_perc_area_LAPA = 10^(-3)*SpatialIntegration(FC10{2}(PA_LAPA==1),R);
uncert_FC_perc_area_LAPA = 10^(-3)*SpatialIntegration(917*delta_FAC10{2}(PA_LAPA==1), R);
fprintf('FC10 in LAPA: \n%0.2f +/- %0.2f\n', FC_perc_area_LAPA, uncert_FC_perc_area_LAPA)

FC_perc_area_HAPA = 10^(-3)*SpatialIntegration(FC10{2}(PA_HAPA==1),R);
uncert_FC_perc_area_HAPA = 10^(-3)*SpatialIntegration(917*delta_FAC10{2}(PA_HAPA==1), R);
fprintf('FC10 in HAPA: \n%0.2f +/- %0.2f\n', FC_perc_area_HAPA, uncert_FC_perc_area_HAPA)

fprintf('\nHarper''s perc zone area: %0.2f km2\n',SpatialIntegration(is_perc_area_mat(is_perc_area_mat==1),R)*10^3)
temp = or (LAPA == 1, HAPA ==1);
fprintf('our LAPA + HAPA: %0.2f km2\n',SpatialIntegration(temp(:),R)*10^3)

fprintf('\nHarper''s LAPA area: %0.2f %\n',SpatialIntegration(is_perc_area_mat(PA_LAPA==1),R)/SpatialIntegration(is_perc_area_mat(is_perc_area_mat==1),R)*100)
fprintf('\nHarper''s HAPA area: %0.2f %\n',SpatialIntegration(is_perc_area_mat(PA_HAPA==1),R)/SpatialIntegration(is_perc_area_mat(is_perc_area_mat==1),R)*100)

figure
hold on 
PlotBackground(GL,Ice,Firn);
% h = pcolor(XX,YY,T_raster);
% h.LineStyle = 'none';
temp = PA_LAPA+0;
temp(temp==0) = NaN;
h2 = surf(XX,YY,temp);
h2.LineStyle = 'none';
h2.FaceColor = RGB('rouge');
temp = PA_HAPA+0;
temp(temp==0) = NaN;
h2 = surf(XX,YY,temp+1);
h2.LineStyle = 'none';
h2.FaceColor = RGB('vert');
set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')

[~, ~] = contour(XX,YY,DSA, [1 1],'Color','k', 'LineWidth',1);
[~, ~] = contour(XX,YY,HAPA, [0.5 0.5],'Color','k', 'LineWidth',1);
[~, ~] = contour(XX,YY,LAPA, [0.5 0.5],'Color','k', 'LineWidth',1);
axis tight
daspect( [1 1 1])
xlim(1.0e+05 *[-6.2164    8.4827])
ylim(1.0e+06 *[-3.3439   -0.6732])
title('Harper et al.''s percolation area') 

%% Printing to files
disp('')

disp('Printing to raster files')
    mkdir('Output/FAC10 maps')
    
% metadata
temp = datevec(datenum(metadata.Date));
metadata.Month = temp(:,2);
metadata.Day = temp(:,3);

writetable(metadata,sprintf('Output/metadata_%s.csv',orig),'Delimiter',';');

% FAC10

    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_DSA_1953_2017_%s.tif', orig),...
        FAC10_map{2}, DSA, R)
    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_DSA_1953_2017_u_%s.tif', orig),...
        delta_FAC10{2}, DSA, R)

    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_LAPA_1998_2008_%s.tif', orig),...
        FAC10_map{1}, LAPA, R)
    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_LAPA_1998_2008_u_%s.tif', orig),...
        delta_FAC10{1}, LAPA, R)

    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_PA_2010_2017_%s.tif', orig),...
        FAC10_map{2}, LAPA+HAPA, R)
    PrintFACmap(sprintf('Output/FAC10 maps/FAC10_PA_2010_2017_u_%s.tif', orig),...
        delta_FAC10{2}, LAPA+HAPA, R)



% FC10
    PrintFACmap(sprintf('Output/FAC10 maps/FC10_post2010_%s.tif', orig),...
        FC10{2}, ~isnan(FC10{2}), R)
    PrintFACmap(sprintf('Output/FAC10 maps/FC10_post2010_u_%s.tif', orig),...
        917*delta_FAC10{2}, ~isnan(FAC10_map{2}), R)

    PrintFACmap(sprintf('Output/FAC10 maps/FC10_pre2010_%s.tif', orig),...
        FC10{1}, ~isnan(FC10{1}), R)
    PrintFACmap(sprintf('Output/FAC10 maps/FC10_pre2010_u_%s.tif', orig),...
        917*delta_FAC10{1}, ~isnan(FC10{1}), R)
%FC100
    PrintFACmap(sprintf('Output/FAC10 maps/FC100_post2010_%s.tif', orig),...
        FC100{2}*a_FAC, ~isnan(FC100{2}), R)
    PrintFACmap(sprintf('Output/FAC10 maps/FC100_post2010_u_%s.tif', orig),...
        delta_FAC10{2}*a_FAC, ~isnan(FC100{2}), R)

    PrintFACmap(sprintf('Output/FAC10 maps/FC100_pre2010_%s.tif', orig),...
        FC100{1}*a_FAC, ~isnan(FC100{1}), R)
    PrintFACmap(sprintf('Output/FAC10 maps/FC100_pre2010_u_%s.tif', orig),...
        917*delta_FAC10{1}*a_FAC, ~isnan(FC100{1}), R)

%% difference Box MAR
    if source_temp_accum == 1
        FAC_10_map_Box13 = FAC10_map;
        FAC_10_spread_Box13 = delta_FAC10;
        metadata_Box13 = metadata;
        R_Box13 = R;
        save ('./Output/Box13_result.mat','FAC_10_map_Box13','metadata_Box13','R_Box13','FAC_10_spread_Box13') ;
        load('./Output/MAR_result.mat')
    else
        FAC_10_map_MAR = FAC10_map;
        metadata_MAR = metadata;
        FAC_10_spread_MAR = delta_FAC10;

        R_MAR = R;
        save('./Output/MAR_result.mat','FAC_10_map_MAR','metadata_MAR','R_MAR','FAC_10_spread_MAR') 
        load('./Output/Box13_result.mat') 
    end
    
for i = 1:2
    FAC_10_map_MAR_resampled{i} = ...
        griddata(XX_mar,YY_mar,FAC_10_map_MAR{i},XX_box,YY_box);
    FAC_10_spread_MAR_resampled = ...
        griddata(XX_mar,YY_mar,FAC_10_spread_MAR{i},XX_box,YY_box);
end

diff_MAR_box_1 = abs(FAC_10_map_MAR_resampled{1}-FAC_10_map_Box13{1});
diff_MAR_box_1(FAC_10_map_Box13{1,1} == 0) = NaN;

diff_MAR_box_2 = abs(FAC_10_map_MAR_resampled{2}-FAC_10_map_Box13{2});
diff_MAR_box_2(FAC_10_map_Box13{1,1} == 0) = NaN;
    
PlottingDifferenceBoxMAR(XX_box,YY_box,diff_MAR_box_1,diff_MAR_box_2,...
    XX, YY, DSA, HAPA, LAPA, GL, Ice, Firn, vis);

%% Uncertainty on pore space
disp('Uncertainty on pore space')

years = unique(metadata.Year);
  count_plot = 0 ;  
out_table = table;
for i = 1:length(years)
    ind_year = find(years(i) == metadata.Year);

    for j = 1:length(ind_year)
        dist_cores = distance(metadata.Latitude(ind_year(j)), ...
            metadata.Longitude(ind_year(j)),...
            metadata.Latitude(ind_year), ...
            metadata.Longitude(ind_year));
        if j>1
            if ismember(j,nearby_cores)
                continue
            end
        end
        nearby_cores = find(dist_cores < 3/111);
        if length(nearby_cores)>1
            M = table(metadata.Name(ind_year(nearby_cores)), ...
                metadata.Year(ind_year(nearby_cores)), ...
                metadata.FAC10(ind_year(nearby_cores)),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*mean(metadata.FAC10(ind_year(nearby_cores))),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*std(metadata.FAC10(ind_year(nearby_cores))),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*numel(ind_year(nearby_cores)),...
                metadata.Citation(ind_year(nearby_cores)),'VariableNames',...
                {'Name','Year','FAC_10m','mean_FAC','std_FAC','NumCores','Citation'});
            
            out_table = [out_table; M];
            
%              f =PlotCore(Core,'CoreNumber',metadata.CoreNumber(ind_year(nearby_cores)),...
%                  'PlotStrat','simple',...
%                 'Ylim',15); 
%             pause(0.5)
%             count_plot = count_plot+1;
%             if count_plot ==30
%                 close all
%                 count_plot=0;
%             end

        end
        
 
    end
end
                close all

out_table(out_table.std_FAC==0,:) = [];
out_table(isnan(out_table.std_FAC),:) = [];

out_table([8:29 35:36],:) = [];

[~, ind_sorted] = sort(out_table.Citation(:));
out_table= out_table(ind_sorted,:);

[C,ind_same_site,IC] = unique(out_table.std_FAC,'rows','stable');
[~,~, out_table.Citation_id] = unique(out_table.Citation(:));


% symbol = {'x','o','d','^','v','s','p','h','<','>','+','*','.'};
f = figure;
ha = tight_subplot(1,1,0.01, [.15 .3], 0.1);
set(f,'CurrentAxes',ha(1))

hold on
leg_text = {'Calculated from:','2 cores','4 cores',...
    '10 cores','','Source:','Spencer et al. (2000)'...
    'Forster et al. (2014)', 'Koenig et al. (2014)', ...
    'Harper et al. (2012)', 'Machguth et al. (2016)',...
    'Morris and Wingham (2014)','this study'};
% leg_text = {};
col = lines(max(out_table.Citation_id));

h(1) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));
h(2) = plot(NaN,NaN,'ok','MarkerSize',sqrt(50*2));
h(3) = plot(NaN,NaN,'ok','MarkerSize',sqrt(50*4));
h(4) = plot(NaN,NaN,'ok','MarkerSize',sqrt(50*10));
h(5) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));
h(6) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));
% h(8) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));
% h(9) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));
% h(10) = plot(NaN,NaN,'ow','MarkerSize',sqrt(50*2));

for i = unique(out_table.Citation_id)' %for each study
    ind = find(out_table.Citation_id(ind_same_site)==i);
    scatter(out_table.Year(ind_same_site(ind)),...
        out_table.std_FAC(ind_same_site(ind)),...
        50*out_table.NumCores(ind_same_site(ind)),col(i,:),'MarkerFaceColor',col(i,:))
%     leg_text = {leg_text{:}, out_table.Citation{ind_same_site(ind(1))}};
end

xlimit= get(gca,'xlim');
plot(xlimit,0.3*[1 1],'--k')
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',18,...
            'LineWidth',1.2,'Ticklength',[0.04 0.08])
        legendflex(leg_text, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'ncol',3, ...
                       'fontsize',18,...
                       'box','off',...
                       'Interpreter','tex');
% dx = 0; dy = 0; % displacement so the text does not overlay the data points
% label_x = 2020 * ones(size(out_table.Year(IA)));
% label_y = flipud(linspace(0.02,0.98,size(out_table.Year(IA),1))');
% dx = out_table.Year(IA)- label_x;
% dy = out_table.std_FAC(IA) - label_y;
% axis tight square
% 
% xlimits = get(gca,'XLim');
% ylimits = get(gca,'YLim');
% 
% quiver(out_table.Year(IA),out_table.std_FAC(IA),...
%     -dx,...
%     -dy,0,...
%     'Color','r',...
%     'ShowArrowHead','off')

% text(label_x, label_y, out_table.Var2(IA));
% ylim([0 1])
% xlim([1980 2033])
box on
axis tight square
xlabel('Year of survey')
ylabel('Standard deviation \newlineof FAC_{10} (m)','Interpreter','tex')

print(f,'Output/Standard_Deviation.tif','-dtiff')

writetable(out_table, ...
sprintf('./Output/Uncertainty_table_%s.csv',orig),'Delimiter',';');
