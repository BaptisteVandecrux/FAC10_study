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

%% Graph options and path setting
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

vis = 'off'; % make plot visible. If 'off' plots are not displayed but still printed
mkdir('.\Output')
addpath(genpath('.\lib'))
addpath(genpath('.\Input'))

%% ========= Here you choose the T and b map source ====================
% 1 is Box13 2 is MAR 3 is HIRHAM
 source_temp_accum = 2;
switch source_temp_accum 
    case 1
         accum_thresh = 600;
         T_thresh = -16;
         orig = 'Box13';
     case 2 
         accum_thresh = 600;
         T_thresh = -20;
         orig = 'MAR';
end

filename= sprintf('Output/result_%s.txt', orig);
diary(filename)

%% Loading temperature and accumulation data
disp('Loading long-term temperature and accumulation maps')
switch source_temp_accum
    case 1
            b_map = table;
            T_map = table;
        if exist('mean_temperature_box13.csv') == 2
            temp = dlmread('mean_temperature_box13.csv',';');
            T_map.lon = temp(:,1);
            T_map.lat = temp(:,2);
            T_map.T_avg = temp(:,3);
            
            temp = dlmread('mean_accumulation_box13.csv',';');
           	b_map.lon = temp(:,1);
            b_map.lat = temp(:,2);
            b_map.b_avg = temp(:,3);
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

            b_map.lon = lon(:);
            T_map.lat = lat(:);
            b_map.lon = lon(:);
            T_map.lat = lat(:);

            b_map.b_avg = mean (acc(:,:,131:end),3);
            T_map.T_avg = mean (Temperature(:,:,131:end),3);
            
            dlmwrite('./Output/Temp accum maps/mean_temperature_box13.csv', [lon(:) lat(:) T_map.T_avg(:)],'Delimiter',';');
            dlmwrite('./Output/Temp accum maps/mean_accumulation_box13.csv', [lon(:) lat(:) b_map.b_avg(:)],'Delimiter',';');
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
        b_map = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','b_avg'});
        clearvars filename delimiter formatSpec fileID dataArray ans;
end
        b_map.b_avg(b_map.b_avg==-999) = NaN;
        T_map.T_avg(T_map.T_avg==-999) = NaN;
        lon =  reshape(T_map.lon,561,301);
        lat =  reshape(T_map.lat,561,301);

%% Loading FAC10 dataset
disp('Loading FAC10 dataset')
if exist('FAC10 dataset.xlsx') == 2
    [~, ~, raw] = xlsread('.\Input\FAC10 dataset.xlsx','Sheet1','A2:I345');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[2,9]);
    raw = raw(:,[1,3,4,5,6,7,8]);
    data = reshape([raw{:}],size(raw));
    metadata = table;
    metadata.Year = data(:,1);
    metadata.Name = cellVectors(:,1);
    metadata.Latitude = data(:,2);
    metadata.Longitude = data(:,3);
    metadata.T_avg = data(:,4);
    metadata.b_avg = data(:,5);
    metadata.DepthMax = data(:,6);
    metadata.FAC10 = data(:,7);
    metadata.Citation = cellVectors(:,2);
    clearvars data raw cellVectors;
else
    metadata = Create_FAC10_dataset(T_map,b_map,vis);
end

% The metadata is distributed with average temperature and accumulation
% from Box (2013) and Box et al. (2013).
% If working with MAR, then they need to be overwritten.
if source_temp_accum == 2
    for i = 1:height(metadata)        
        % accumulation
         [~, ind] = min(distance( metadata.Latitude(i),metadata.Longitude(i),...
             b_map.lat,b_map.lon));
        metadata.b_avg(i) = b_map.b_avg(ind);
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
b_map.b_avg(find(~is_firn)) = NaN;

% Loading raster files
disp('Applying on the T_avg map')

Ice = shaperead('IcePolygon_3413.shp');
GL = shaperead('GL_land.shp');
[FirnArea_3413, R_fa3413] = geotiffread('mean_temperature_box13_3413_firn.tif');
Firn = shaperead('FirnLayer2000-2017_final_3413.shp');

switch source_temp_accum
    case 1
        [T_raster, ~] = geotiffread('mean_temperature_box13_3413_firn.tif');
        [b_raster, R] = geotiffread('mean_accumulation_box13_3413_firn.tif');
        I = geotiffinfo('mean_temperature_box13_3413_firn.tif'); 

    case 2
        [T_raster, ~] = geotiffread('T_avg_1979-2014_MAR_3413_firn.tif');
        [b_raster, R] = geotiffread('Net_Snowfall_avg_1979-2014_MAR_3413_firn.tif');
        I = geotiffinfo('T_avg_1979-2014_MAR_3413_firn.tif'); 
end

        T_raster(T_raster==0)=NaN;
        T_raster(T_raster==-999)=NaN;
        T_raster(T_raster==-9999)=NaN;
        T_raster(T_raster<-9999)=NaN;
        b_raster(b_raster==0)=NaN;
        b_raster(b_raster==-999)=NaN;
        b_raster(b_raster==-9999)=NaN;
        b_raster(b_raster<-9999)=NaN;
        
        [x_map,y_map]=pixcenters(I);
        XX = repmat(x_map,length(y_map),1);
        YY = repmat(y_map',1,length(x_map));
      
        % Cold area
        DSA = (T_raster<=T_thresh)+0;
        % Warm High Accumulation area
        HAWSA = and(T_raster>T_thresh,b_raster>accum_thresh)+0;
        % Warm Low Accumulation area
        LAWSA = and(T_raster>T_thresh,b_raster<=accum_thresh)+0;
        
        format BANK 
disp('Firn area (km2):')
fprintf('%0.2f\n\n',sum(sum(~isnan(T_raster))) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);

disp('DSA (km2):')
fprintf('%0.2f\n',sum(sum(DSA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %% \n\n',sum(sum(DSA)) /sum(sum(~isnan(T_raster)))*100 );

disp('LAWSA (km2):')
fprintf('%0.2f\n',sum(sum(LAWSA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %%\n\n',sum(sum(LAWSA)) /sum(sum(~isnan(T_raster)))*100 );

disp('HAWSA (km2):')
fprintf('%0.2f\n',sum(sum(HAWSA)) * R.CellExtentInWorldX *R.CellExtentInWorldY /1000000);
fprintf('%0.1f %%\n\n',sum(sum(HAWSA)) /sum(sum(~isnan(T_raster)))*100 );

    x = T_map.T_avg(~isnan(T_map.T_avg));
    y = b_map.b_avg(~isnan(T_map.T_avg));
    k = boundary(x, y, 0.8);
    poly_100 = [x(k), y(k)];

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
    [accum_thresh max(metadata.b_avg)+1 max(metadata.b_avg)+1 accum_thresh],...
    vert)
patch([min(metadata.T_avg)-1 min(metadata.T_avg)-1 T_thresh T_thresh ],...
    [0 max(metadata.b_avg)+1 max(metadata.b_avg)+1 0],...
    jaune)
patch([T_thresh T_thresh max(metadata.T_avg)+1 max(metadata.T_avg)+1],...
    [0 accum_thresh accum_thresh 0],...
    rouge)

g = patch(poly_100(:,1),poly_100(:,2), RGB('light light gray'));
g.FaceColor = 'w';
g.LineWidth = 1;
alpha(g,0.35)

xlim([min(metadata.T_avg)-1 max(metadata.T_avg)+1])
ylim([0 max(metadata.b_avg)+1])
% plot(x,(x+28)*20+320,'k','LineWidth',2)
% plot(get(gca,'Xlim'),accum_thresh*ones(1,2),'--k','LineWidth',2)
% plot([T_thresh T_thresh],[0 1]*accum_thresh,'--k','LineWidth',2)

h_s1 = scatter(metadata.T_avg(ind_sorted),...
    metadata.b_avg(ind_sorted),...
    140,...
    metadata.FAC10(ind_sorted)...ones(size(metadata.b_avg(ind_sorted))) * 0.8*[1 1 1], ...
    ...+ (metadata.b_avg(ind_sorted)>accum_thresh) * 0.8*[1 0 0],...0.*col(discretize(metadata.Year,time_bin),:),...
    ,'o','fill');

h_s2 = plot(metadata.T_avg(ind_sorted),...
    metadata.b_avg(ind_sorted),'ok','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor','w');

cc = colorbar('NorthOutside');
cc.Label.String = 'Firn air content (m^3 m^{-2})';
colormap('toh')

h_leg = legend([g h_s2],'Firn area',...
    'Firn cores',...
    'Location','NorthWest');

set(gca,'Ticklength',[0.08 0.16]/2,'layer','top','yaxislocation','right')
box on
ylabel('$\mathrm{\overline{\dot{b}} \:(mm w.eq. yr^{-1})}$','Interpreter','latex')
xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex') 

% Explaining dataset division
ind_HighAcc_HighTemp = and(metadata.b_avg>accum_thresh,metadata.T_avg>T_thresh);
ind_Rest = or(metadata.b_avg<=accum_thresh,metadata.T_avg<=T_thresh);
ha(1).Position = [0.45 0.42 0.4 0.3];

set(f,'CurrentAxes',ha(6))
ind_jaunevert = or(metadata.T_avg<T_thresh,metadata.b_avg>accum_thresh);
    % ha(6).Position(2) = 0.28;
    scatter(metadata.T_avg(ind_jaunevert),metadata.FAC10(ind_jaunevert),20,'k','fill')
    hold on
    % plot([T_thresh T_thresh],[0 7],'--k','LineWidth',2)
    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    box on

    axis tight
    xlimits = get(gca,'XLim');
    pp1 = patch([xlimits(1) xlimits(1) T_thresh-0.2 T_thresh-0.2 ],...
        [0 7 7 0],...
        jaune);
    pp2 = patch([T_thresh-0.2 T_thresh-0.2 xlimits(2)   xlimits(2)],...
        [0 7 7 0],...
        vert);
    uistack(pp1,'bottom')
    uistack(pp2,'bottom')

    h_ylab = ylabel('FAC_{10} (m^3 m^{-2})','Interpreter','tex');
set(gca,'Ticklength',[0.08 0.16]/2.5,'Ticklength',[0.08 0.16]/2,...
    'XMinorTick','on','YMinorTick','on','XAxisLocation','bottom','YAxisLocation','right','layer','top')
axis tight

ylim([0 7])

set(f,'CurrentAxes',ha(8))
ha(8).Visible = 'off';

set(f,'CurrentAxes',ha(5))
    % ha(6).Position(2) = 0.28;
    scatter(metadata.T_avg(ind_Rest),metadata.FAC10(ind_Rest),20,'k','fill')
    hold on
    % plot([T_thresh T_thresh],[0 7],'--k','LineWidth',2)
    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    box on
    set(gca,'Ticklength',[0.08 0.16]/2.5,'XMinorTick','on','YMinorTick','on',...
        'XAxisLocation','bottom','YAxisLocation','left','layer','top')
    axis tight
    xlimits = get(gca,'XLim');
    pp1 = patch([xlimits(1) xlimits(1) T_thresh T_thresh ],...
        [0 7 7 0],...
        jaune);
    pp2 = patch([T_thresh T_thresh xlimits(2)   xlimits(2)],...
        [0 7 7 0],...
        rouge);
    uistack(pp1,'bottom')
    uistack(pp2,'bottom')
    ylim([0 7])
    h_ylab = ylabel('FAC_{10} (m^3 m^{-2})','Interpreter','tex');

set(f,'CurrentAxes',ha(7))
ha(7).Visible = 'off';

for i=1:8
    set(f,'CurrentAxes',ha(i))
    set(gca,'FontSize',14,'FontName','Times New Roman')
end
 
set(f,'CurrentAxes',ha(2))
ha(2).Position = [-0.005 0.33 0.53 0.53];
hold on
for i = 1:length(GL)
    ind = isnan(GL(i).X+GL(i).Y);
    GL(i).Y(ind)=[];
    GL(i).X(ind)=[];
    patch(GL(i).X,GL(i).Y,RGB('gray'))
end
for i = 1:length(Ice)
    ind = isnan(Ice(i).X+Ice(i).Y);
    Ice(i).Y(ind)=[];
    Ice(i).X(ind)=[];
    patch(Ice(i).X,Ice(i).Y,RGB('light light gray'))
end

for kk = 1: length(Firn)
    fill(Firn(kk).X,Firn(kk).Y,'w')
end
hold on
% plot(XX(ind0),YY(ind0),'.')
CA_save = DSA;
WHA_save = HAWSA;
WLA_save = LAWSA;
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',3);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',3);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',3);
        DSA(DSA == 0)=NaN;
        HAWSA(HAWSA == 0)=NaN;
        LAWSA(LAWSA == 0)=NaN;
p1 = surf(XX,YY,DSA*0, DSA);
p2 = surf(XX,YY,HAWSA*0, HAWSA+1);
p3 = surf(XX,YY,LAWSA*0, LAWSA+2);
p1.LineStyle = 'none';
p2.LineStyle = 'none';
p3.LineStyle = 'none';
p1.FaceColor = jaune;
p2.FaceColor = vert;
p3.FaceColor = rouge;

plot(metadata.X, metadata.Y,'ok',...
    'MarkerFaceColor','w','MarkerSize',3)
            xlim(1.0e+05 *[-6.2164    8.4827])
            ylim(1.0e+06 *[-3.3439   -0.6732])
            set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
         h_leg_map =  legendflex([p1 p2 p3], {'Dry Snow Area','High Accumulation Wet Snow Area',...
             'Low Accumulation Wet Snow Area'},...
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
    [0.17 0.24 0.04 0.05],...
    'String','c)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');


% Create textbox
annotation(f,'textbox',...
    [0.53 0.24 0.04 0.05],...
    'String','d)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');

print(f,sprintf('Output/data_presentation_%s',orig),'-dtiff')
print(f,sprintf('Output/data_presentation_%s',orig),'-dpdf','-r0')

%% climatic control on firn line
figure('Visible',vis)
hold on
ind_no_nan=find(~isnan(T_map.T_avg));

scatter(T_map.lon(ind_no_nan),T_map.lat(ind_no_nan))
k = boundary(T_map.lon(ind_no_nan),T_map.lat(ind_no_nan),1);
scatter(T_map.lon(ind_no_nan(k)),T_map.lat(ind_no_nan(k)))
plot(T_map.lon(ind_no_nan(k)),T_map.lat(ind_no_nan(k)),'LineWidth',2)

ind_firn_line = ind_no_nan(k);

figure('Visible',vis)
scatter(T_map.T_avg(ind_firn_line),...
    b_map.b_avg(ind_firn_line))
[x, ind_sorted] = sort(T_map.T_avg(ind_firn_line));
temp = b_map.b_avg(ind_firn_line);
y=temp(ind_sorted);
firn_line = fit(x,y,'exp1');
hold on 
plot(x, firn_line(x),'r','LineWidth',2)
box on
ylabel('Accumulation')
xlabel('Temperature')

accum_zerofirn = b_map.b_avg(ind_firn_line);
temp_zerofirn = T_map.T_avg(ind_firn_line);

%% ===== Region 1: Cold elevated regions ===============
disp('Selecting points in the cold area')
metadata = metadata_save;
time_bin = [1949.50       1959.50       1969.50       1979.50 ...
    1989.50       1996.50      2009.5  2019.50];
col = lines(length(time_bin)-1);

ind = metadata.T_avg>=T_thresh;
metadata(ind,:) = [];
fprintf('\nUsing %i cores to build the DSA\n',size(metadata,1))
% 3D exploration of accumulation and temperature control
f = figure('Visible',vis,'outerposition',[1 0 20 30]);
ha = tight_subplot(1,1,0.1,0.3,0.3);
set(f,'CurrentAxes',ha(1))
hold on

% plotting the core derived FAC
scatter3(metadata.b_avg, ...
    metadata.T_avg,metadata.FAC10,...
    100,col(discretize(metadata.Year,time_bin),:),'fill') 

hold on

% Plotting NH surface
    Eqn = 'FAC_NH_func(x,y,a)';
    FAC_fit_CA = fit([metadata.b_avg, metadata.T_avg],metadata.FAC10,Eqn);
      
    disp(FAC_fit_CA.a)
xx = linspace(min(metadata.b_avg),max(metadata.b_avg));
yy = linspace(min(metadata.T_avg),max(metadata.T_avg));
[gridx, gridy] = meshgrid(xx,yy);

h_HL = surf(gridx,gridy,FAC_fit_CA(gridx,gridy),3*ones(size(gridx)));
alpha(h_HL, 0.5);
h_HL.LineStyle = 'none';

ME = mean(metadata.FAC10-FAC_fit_CA(metadata.b_avg, metadata.T_avg));
RMSE = sqrt(mean((metadata.FAC10-FAC_fit_CA(metadata.b_avg, metadata.T_avg)).^2));
axis tight 
    zlim([0 6])   
 box on
 
view(ha,[-76.4 27.6]);
 uncert_DSA = RMSE*2;
h_title = title('a) Dry Snow Area','Interpreter','tex');
h_title.Units = 'Normalized';
h_title.Position(2) = 1.4;
h_xlab = xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
h_xlab.Units = 'Normalized';
h_xlab.Position = [1.2   -0.07 0];
ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
zlabel('FAC_{10} (m^{3} m^{-2})','Interpreter','tex')
time_bin2 = 1949.50 :10:2019.50 ;
colbar =contourcmap('lines',time_bin2,'colorbar','on','YLabelString','Year');
colbar.Position(1) = 0.75;
colbar.YTickLabel = time_bin2+0.5;
colbar.YTickLabel(end,:) = '2017';
colbar.YTickLabel(7,:) = '2010';
colbar.YTickLabel(6,:) = '1998';

annotation(f,'textarrow',0.6*[1 1],[0.75 0.6])

h_leg = annotation('textbox','String',sprintf(['Densification law from Arthern et al. (2010) \\newline', ...
    '\\rho_{fresh snow} = %0.2f kg m^{-3}   RMSD = %0.3f m^{3} m^{-2}'],FAC_fit_CA.a,RMSE),'Interpreter','tex');
h_leg.Units = 'Normalized';
h_leg.Position(2) = 0.73;
h_leg.Position(1) = 0.14;
h_leg.FontSize = 15;
h_leg.FitBoxToText = 'on';
h_leg.BackgroundColor = 'w';


print(f,sprintf('./Output/ColdAreas_%s',orig),'-dtiff')

% optional plot to fit the NH dry compaction function to each year
% individually:
% FAC_fit_CA_yr = {};
% figure('Visible',vis)
% 
% hold on
%  ind_binned = discretize(metadata.Year,time_bin);
%     for i = 1: length(time_bin)
%         x = metadata.b_avg(ind_binned == i);
%         y = metadata.T_avg(ind_binned == i);
%         z = metadata.FAC10(ind_binned == i);
%         if length(x)<5
%             continue
%         end
%              % plotting the core derived FAC
%         scatter3(x,y,z,...
%             70,col(i,:),'fill') 
%         format shortg
%         Eqn = 'FAC_NH_func(x,y,a)';
%         FAC_fit_CA_yr{i} = fit([x, y],z,Eqn);
% 
%         xx = linspace(min(metadata.b_avg),max(metadata.b_avg));
%         yy = linspace(min(metadata.T_avg),max(metadata.T_avg));
%         [gridx, gridy] = meshgrid(xx,yy);
% 
%         h_HL = surf(gridx,gridy,FAC_fit_CA_yr{i}(gridx,gridy),3*ones(size(gridx)));
%         alpha(h_HL, 0.5);
%         h_HL.LineStyle = 'none';
%         h_HL.FaceColor = col(i,:);
%         view([-76.4 27.6]);
%         RMSD_period  = (mean((z - FAC_fit_CA(x,y))));
%         fprintf('%i   %i   %0.2f   %0.2f   %0.2f\n', round(time_bin(i)), round(time_bin(i+1)),...
%             nanmean(nanmean(FAC_fit_CA_yr{i}(gridx,gridy))) ,FAC_fit_CA_yr{i}.a, RMSD_period)
% 
% %         pause
%     end

%% ===== Region 2: High accumulation areas =============
disp('Selecting points in the HAWSA')
metadata = metadata_save;

ind = or(metadata.T_avg<T_thresh-1,metadata.b_avg<accum_thresh);
% ind = metadata.T_avg<T_thresh;
metadata(ind,:) = [];
ind = metadata.Year<2010;
metadata(ind,:) = [];

fprintf('\nUsing %i cores to build the HAWSA\n',size(metadata,1))

x1 = metadata.T_avg;
x2 = metadata.b_avg;
z = metadata.FAC10;

f = figure('Visible', vis, 'outerposition',[0 0 20 20]);
subplot(2,1,1)
scatter(x1,z)
hold on
[lm, h] = Plotlm(x1,z,'Annotation','off','LineStyle','-','LineWidth',1);
title(sprintf('FAC_10 = %0.2f * T_a + %0.2f \n RMSD = %0.2f p-value = %0.2f',...
    lm.Coefficients.Estimate(2), lm.Coefficients.Estimate(1), ...
    lm.RMSE ,max(lm.Coefficients.pValue)))
xlabel('Ta')
ylabel('FAC_10')
box on

subplot(2,1,2)
scatter(x2,lm.Residuals.Raw)
hold on
[lm2, h] = Plotlm(x2,lm.Residuals.Raw,'Annotation','off','LineStyle','-','LineWidth',1);
title(sprintf('FAC_10 = %0.2g * b_a + %0.2f \n RMSD = %0.2f p-value = %0.2f',...
    lm2.Coefficients.Estimate(2), lm2.Coefficients.Estimate(1), ...
    lm2.RMSE ,max(lm2.Coefficients.pValue)))

% title(sprintf('Mean = %0.2f   RMSD = %0.2f',mean(z),std(z)))
xlabel('Accumulation')
ylabel('FAC_10')
box on
print(f,sprintf('./Output/LinReg_HighAccumulationArea_%s',orig),'-dtiff')

% preparing the surfaces
% metadata = metadata_save;
% ind = or(metadata.T_avg<T_thresh,metadata.b_avg<accum_thresh);
% % ind = metadata.T_avg<T_thresh;
% metadata(ind,:) = [];
    T_extrap = min(metadata.T_avg):0.5:-2;
FAC_nofirn = 0;

x_all = [metadata.b_avg; firn_line(T_extrap')];
y_all = [metadata.T_avg; T_extrap'];
z_all = [metadata.FAC10; FAC_nofirn*ones(size(T_extrap'))];

xx = linspace(accum_thresh-100,3000,1000);
yy = linspace(T_thresh-5,-2,1000);

[sf_mid, gridx, gridy] = gridfit(x_all, y_all,z_all, xx,yy,'smoothness',4);
sf_mid = min(sf_mid,lm.Coefficients.Estimate(2)*gridy + lm.Coefficients.Estimate(1));
sf_mid = max(0,sf_mid);
FAC_fit_HA_mid = @(x,y) interp2(xx,yy,sf_mid,x,y);

sf_high = 0.*gridx;
sf_high(gridx>reshape(firn_line(gridy),size(gridy))) = ...
    lm.Coefficients.Estimate(2)*gridy(gridx>reshape(firn_line(gridy),size(gridy))) + lm.Coefficients.Estimate(1);
FAC_fit_HA_high = @(x,y) interp2(xx,yy,sf_high,x,y);

sf_low = 0.*gridx;
sf_low(gridx>reshape(firn_line(gridy),size(gridy))) = ...
    lm.Coefficients.Estimate(2)*gridy(gridx>reshape(firn_line(gridy),size(gridy))) ...
    + lm.Coefficients.Estimate(1);
line_low = fit(metadata.b_avg,metadata.T_avg,'poly1');
sf_low(gridy>reshape(line_low(gridx)+2,size(gridy))) = 0;
sf_low = min(sf_low,sf_mid);

FAC_fit_HA_low = @(x,y) interp2(xx,yy,sf_low,x,y);

xx_plot = linspace(accum_thresh,3000,200);
yy_plot = linspace(T_thresh,-1,1000);
[gridx,gridy] = meshgrid(xx_plot, yy_plot);

%%
f = figure('Visible', vis, 'outerposition',[1 0 30 30]);
ha = tight_subplot(1,2,0.1,0.25,[0.2 0.2]);

    set(f,'CurrentAxes',ha(1))
hold on
plot(b_map.b_avg(ind_firn_line),T_map.T_avg(ind_firn_line),'x','Color',col(end,:),'LineWidth',1.5)
[x, ind_sorted] = sort(T_map.T_avg(ind_firn_line));
temp = b_map.b_avg(ind_firn_line);
y=temp(ind_sorted);
firn_line = fit(x,y,'exp1');
plot(firn_line(x),x,'k','LineWidth',3)
box on
h_xlab = xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
axis tight
ha(1).Position =   [0.07         0.1         0.2          0.25];
title('Location of the firn layer boundary','FontSize',11)
h_l = legend('Remotely-sensed','Idealized','Location','SouthEast');
h_l.FontSize = 11;
% legend boxoff

set(f,'CurrentAxes',ha(2))
ha(2).Position = [0.34         0.23         0.3          0.45];

hold on

% plotting the core derived FAC
h_FAC = scatter3(metadata.b_avg, ...
    metadata.T_avg,metadata.FAC10,...
    100,col(end,:),'fill'); 

hold on
%     h_FL = plot3(accum_zerofirn, temp_zerofirn, ...
%         FAC_nofirn*ones(size(temp_zerofirn)),...
%         'x','Color',col(end,:),'LineWidth',2);
%     alpha(h_FL,0.2);
    T_extrap = [T_thresh-2:1:-2]';
    h_IFL = plot3(firn_line(T_extrap),T_extrap+0.1, ...
        FAC_nofirn*ones(size(T_extrap)),'Color','k','LineWidth',4);
    
% sf(abs(gridx-reshape(firn_line(gridy),size(gridy))) <50) = NaN;
h_s2 = surf(gridx, gridy,FAC_fit_HA_low(gridx,gridy),0*gridy+3);
alpha(h_s2, 0.02)
h_s1 = surf(gridx, gridy,FAC_fit_HA_mid(gridx,gridy),0*gridy+4);
alpha(h_s1, 0.05)

h_s = surf(gridx, gridy,FAC_fit_HA_high(gridx,gridy)+0.01,0*gridy+2);
alpha(h_s, 0.1)


shading interp

[~, h_c ]= contour3(gridx, gridy,FAC_fit_HA_high(gridx,gridy),10);
[~, h_c1 ] = contour3(gridx, gridy, FAC_fit_HA_mid(gridx,gridy),11);
[~, h_c2 ] = contour3(gridx, gridy, FAC_fit_HA_low(gridx,gridy),10);
h_c.Color = 'b';
h_c1.Color = 'r';
h_c2.Color = 'g';
h_c.LineWidth = 1.5;
h_c1.LineWidth = 1.5;
h_c2.LineWidth = 1.5;
colormap('hsv')

RMSE = std(metadata.FAC10);
RMSE1 = std(metadata.FAC10-FAC_fit_HA_high(metadata.b_avg,metadata.T_avg));
RMSE2 = std(metadata.FAC10-FAC_fit_HA_mid(metadata.b_avg,metadata.T_avg));
RMSE3 = std(metadata.FAC10-FAC_fit_HA_low(metadata.b_avg,metadata.T_avg));

title('c) High Accumulation Wet Snow Area');

axis tight 
    zlim([0 6])   
ylim([min(x1)-2 max(x1)+2])
xlim([min(x2)-10 max(x2)+10])

if source_temp_accum == 2
    ylim([min(x1)-2 max(x1)+7])
    xlim([min(x2)-50 max(x2)+10])
end
 box on
view(ha(2),[-66 11]);

h_leg = legend([h_FAC h_IFL h_c h_c1 h_c2],'FAC_{10} observations',...
    'Idealized firn layer boundary',...
    sprintf('\nUpper-range FAC_{10} estimate \n RMSD = %0.2f m^{3} m^{-2}\n',RMSE1),...
    sprintf('\nMid-range FAC_{10} estimate \n RMSD = %0.2f m^{3} m^{-2}\n',RMSE2),...
    sprintf('\nLower-range FAC_{10} estimate \n RMSD = %0.2f m^{3} m^{-2}\n',RMSE3),...
    'Location','East','Interpreter','tex');
h_leg.Position(1) = 0.65;
h_leg.Position(2) = 0.3;
set(gca,'XTickLabel',[' 800'; '    '; '1200'; '    '; '1600';]);
%  title(sprintf('ME = %0.2f   RMSE = %0.2f',ME,RMSE))
h_xlab = xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
h_xlab.Units = 'Normalized';
h_xlab.Position = [1   -0.1         0];
h_ylab = ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex');
h_ylab.Position = [1200      -9      -1.6743];
% Create arrow

if source_temp_accum == 2
    annotation(f,'arrow',[0.4 0.28],...
        [0.28 0.25]);
else
    annotation(f,'arrow',[0.4 0.28],...
        [0.28 0.24]);
end
zlabel('FAC_{10} (m)','Interpreter','tex')
print(f,sprintf('./Output/Extrapolation_HighAccumulationArea_%s',orig),'-dtiff','-r0')

%% ===== Region 3: warm and low accumulation ===========
disp('Selecting points in the warm & low accumulation area')
metadata = metadata_save;
ind = or(metadata.b_avg>accum_thresh,metadata.T_avg<T_thresh-1);
metadata(ind,:) = [];

% 3D exploration of accumulation and temperature control  
f = figure('Visible', vis, 'outerposition',[1 0 20 30]);
ha = tight_subplot(1,1,0.1,0.25,0.25);
set(f,'CurrentAxes',ha(1))
hold on

    % plotting the core derived FAC
%     scatter3(metadata.b_avg, ...
%         metadata.T_avg,metadata.FAC10,...
%         100,col(discretize(metadata.Year,time_bin),:),'fill') 

xlim([min(metadata.b_avg)-2 max(metadata.b_avg)+10])
ylim([min(metadata.T_avg)-1 max(metadata.T_avg)+4])

     box on
h_xlab = xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;

ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
zlabel('FAC_{10} (m)','Interpreter','tex')
    ind_binned = discretize(metadata.Year,time_bin);
    
    i =  6;
    x = metadata.b_avg(ind_binned==i);
    disp('Num points 2000-2009')
    disp(length(x))
    y = metadata.T_avg(ind_binned==i);
    z = metadata.FAC10(ind_binned==i);
    scatter3(x,y, z,...
        100,col(6,:),'fill') 
 fprintf('\nUsing %i cores to build the LAWSA pre2010\n',size(z,1))
   
    xx = linspace(min(metadata.b_avg)-10,accum_thresh+10);
    yy = linspace(T_thresh-5,max(metadata.T_avg)+5);

%     [sf1, gridx, gridy] = gridfit(x,y,z, xx,yy,'smoothness',10);
%     FAC_fit_WHA_pre2010 = @(x,y) interp2(xx,yy,sf1,x,y);


    Eqn = sprintf(['max(0, a*x + ' ...
        'b*min(y+%i,0) + c*max(y+%i,0) + ' ...
        'd*min(y+%i,0) + e*max(y+%i,0) + f)'],...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*2/3),...
        -(min(y)+ (max(y) - min(y))*2/3));
  FAC_fit_WHA_pre2010 = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[Inf 0    0   0   0     Inf],...
        'StartPoint',[0.45414        0.5029       0.72089      0.54005      0.15411      0.12192]);

% mean(mean(SA_WHA_pre2010{ii}(gridx,gridy)))
count = 0;
while sqrt(mean((FAC_fit_WHA_pre2010(x,y)-z).^2))>1
  FAC_fit_WHA_pre2010 = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[Inf 0    0   0   0     Inf]);
    count = count+1;
    disp(count)
end
%     Eqn = 'max(0, a*x + b*y +c*y^2 + d)';
%   FAC_fit_WHA_pre2010 = fit([x, y],z,Eqn,...
%         'Lower',[0     -Inf -Inf  -Inf],...
%         'Upper',[Inf Inf  Inf Inf]);
    xx = linspace(min(metadata.b_avg),accum_thresh+10);
    yy = linspace(T_thresh-1,max(metadata.T_avg)+1);
    [gridx, gridy] = meshgrid(xx,yy);
    
    h_s = surf(gridx, gridy,FAC_fit_WHA_pre2010(gridx, gridy) ,0*gridy+0.5);
    alpha(h_s, 0.2)
    
    [~, h_c ]= contour3(gridx, gridy, FAC_fit_WHA_pre2010(gridx, gridy),15);
    h_c.Color = 'b';
    h_c.LineWidth = 1.5;
RMSE1 = std(z- FAC_fit_WHA_pre2010(x,y));

    i =  7;
    x = metadata.b_avg(ind_binned==i);
    disp('Num points 2010-2017')
    disp(length(x))
    y = metadata.T_avg(ind_binned==i);
    z = metadata.FAC10(ind_binned==i);
        scatter3(x,y, z,...
        100,col(7,:),'fill') 
 fprintf('\nUsing %i cores to build the LAWSA post2010\n',size(z,1))

%     [sf2, gridx, gridy] = gridfit(x,y,z, xx,yy,'smoothness',10);
%     FAC_fit_WHA_post2010 = @(x,y) interp2(xx,yy,sf2,x,y);
    Eqn = sprintf(['max(0, a*x + ' ...
        'b*min(y+%i,0) + c*max(y+%i,0) + ' ...
        'd*min(y+%i,0) + e*max(y+%i,0) + f)'],...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*2/3),...
        -(min(y)+ (max(y) - min(y))*2/3));
    FAC_fit_WHA_post2010 = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[Inf 0    0   0   0     Inf],...
        'StartPoint',[0.45414        0.5029       0.72089      0.54005      0.15411      0.12192]);
    [gridx, gridy] = meshgrid(xx,yy);
%         Eqn = 'max(0, a*x + b*y + c)';
%   FAC_fit_WHA_post2010 = fit([x, y],z,Eqn,...
%         'Lower',[0     -Inf  -Inf],...
%         'Upper',[Inf 0   0]);
    h_s1 = surf(gridx, gridy, FAC_fit_WHA_post2010(gridx, gridy), 0*gridy+3);
    shading interp
    alpha(h_s1, 0.2)
    [~, h_c1 ] = contour3(gridx, gridy, FAC_fit_WHA_post2010(gridx, gridy),15);
    h_c1.Color = 'r';
    h_c1.LineWidth = 1.5;
    zlim([0 6])   
        colormap('hsv')

RMSE2 = std(z- FAC_fit_WHA_post2010(x,y));

title('b) Low Accumulation Wet Snow Area');


annotation(f,'textarrow', [0.5, 0.57], [0.67, 0.58])

temp = {'$$     \widehat{FAC_{10}}$$ $$_{1997-2008}$$' ,['$$RMSD = ' sprintf('%0.2f',RMSE1) 'm^{3} m^{-2}$$'] };
h_leg = text(1,1,temp,...
    'Interpreter','Latex');
h_leg.Units = 'Normalized';
h_leg.Position(2) = 0.89;
h_leg.Position(1) = 0.07;
h_leg.FontSize = 13;
h_leg.BackgroundColor = 'w';
% h_leg.FitBoxToText = 'on';

annotation(f,'textarrow', [0.47, 0.55], [0.58, 0.47])

temp = {'$$     \widehat{FAC_{10}}$$ $$_{2011-2017}$$' ,['$$RMSD = ' sprintf('%0.2f',RMSE2) 'm^{3} m^{-2}$$'] };
h_leg2 = text(1,1,temp,...
    'Interpreter','Latex');
h_leg2.Units = 'Normalized';
h_leg2.Position(2) = 0.73;
h_leg2.Position(1) = 0.07;
h_leg2.FontSize = 13;
h_leg2.BackgroundColor = 'w';
% h_leg.FitBoxToText = 'on';
view(gca,[-95 5]);

print(f,sprintf('./Output/Warm_LowAccumulationArea_%s',orig),'-dtiff')

%% ============= All regions ==================
disp('Outlook for all regions')
metadata = metadata_save;

% 3D exploration of accumulation and temperature control
    % Defining the points over which FAC will be extrapolated
    b_extrap = linspace(0,5000,100);
    T_extrap = linspace(-40, -2, 100);
    gridx  =repmat(b_extrap,length(T_extrap),1);
    gridy = repmat(T_extrap',1,length(b_extrap));

    % Values of FAC derived Naborro Herring equation at these points
    FAC_NH = FAC_NH_func(gridx,gridy,280);
    % FAC_BV = FAC_BV_func(bb,TT,280,1,-5);

    % Selection of the point used for curve fitting 
FAC_fit_all_pre2010_mid = FAC_fit_CA(gridx, gridy);
FAC_fit_all_post2010_mid = FAC_fit_CA(gridx, gridy);
FAC_fit_all_spread = 0* gridy; 

FAC_fit_all_spread(and(gridy>T_thresh,gridx>accum_thresh)) =  FAC_fit_HA_high(gridx(and(gridy>T_thresh,gridx>accum_thresh)), gridy(and(gridy>T_thresh,gridx>accum_thresh)))...
    - FAC_fit_HA_low(gridx(and(gridy>T_thresh,gridx>accum_thresh)), gridy(and(gridy>T_thresh,gridx>accum_thresh)));


FAC_fit_all_post2010_mid(and(gridy>T_thresh,gridx>accum_thresh)) =  FAC_fit_HA_mid(gridx(and(gridy>T_thresh,gridx>accum_thresh)), gridy(and(gridy>T_thresh,gridx>accum_thresh)));

FAC_fit_all_pre2010_mid(and(gridy>=T_thresh,gridx<=accum_thresh)) =  FAC_fit_WHA_pre2010(gridx(and(gridy>=T_thresh,gridx<=accum_thresh)), gridy(and(gridy>=T_thresh,gridx<=accum_thresh)));
FAC_fit_all_post2010_mid(and(gridy>=T_thresh,gridx<=accum_thresh)) =  FAC_fit_WHA_post2010(gridx(and(gridy>=T_thresh,gridx<=accum_thresh)), gridy(and(gridy>=T_thresh,gridx<=accum_thresh)));


% creating the function
FAC_bav_mid_past = @(x,y) interp2(gridx,gridy,FAC_fit_all_pre2010_mid,x,y);
FAC_bav_spread = @(x,y) interp2(gridx,gridy,FAC_fit_all_spread,x,y);
FAC_bav_mid_pres = @(x,y) interp2(gridx,gridy,FAC_fit_all_post2010_mid,x,y);

ind_in_firn = reshape(inpoly([gridx(:) gridy(:)],fliplr(poly_100)),size(gridx));
FAC_fit_all_spread(~ind_in_firn) = NaN;
FAC_fit_all_pre2010_mid(~ind_in_firn) = NaN;
FAC_fit_all_spread(~ind_in_firn) = NaN;
FAC_fit_all_post2010_mid(~ind_in_firn) = NaN;
FAC_fit_all_pre2010_mid(and(gridy>T_thresh,gridx>accum_thresh)) =  NaN;

FAC_fit_all_post2010_high = FAC_fit_all_post2010_mid;
FAC_fit_all_post2010_low = FAC_fit_all_post2010_mid;
FAC_fit_all_post2010_high(and(gridy>T_thresh,gridx>accum_thresh)) =  FAC_fit_HA_high(gridx(and(gridy>T_thresh,gridx>accum_thresh)), gridy(and(gridy>T_thresh,gridx>accum_thresh)));
FAC_fit_all_post2010_low(and(gridy>T_thresh,gridx>accum_thresh)) =  FAC_fit_HA_low(gridx(and(gridy>T_thresh,gridx>accum_thresh)), gridy(and(gridy>T_thresh,gridx>accum_thresh)));
FAC_bav_high_pres = @(x,y) interp2(gridx,gridy,FAC_fit_all_post2010_high,x,y);
FAC_bav_low_pres = @(x,y) interp2(gridx,gridy,FAC_fit_all_post2010_low,x,y);

% plotting
for i = 1:2
    f = figure('Visible', vis, 'outerposition',[1 0 20 30]);
    ha = tight_subplot(1,1,0.1,0.25,0.25);
    set(f,'CurrentAxes',ha(1))
    hold on

        %col(end,:),'fill')
        % plotting the firn line-dervied FAC
        scatter3(accum_zerofirn, temp_zerofirn, ...
            FAC_nofirn*ones(size(temp_zerofirn)),...
            10,col(end,:),'fill')
        ind_binned = discretize(metadata.Year,time_bin);

        if i == 1
            % plotting the core derived FAC
            scatter3(metadata.b_avg(metadata.Year<2010), ...
            metadata.T_avg(metadata.Year<2010),metadata.FAC10(metadata.Year<2010),...
            100,col(discretize(metadata.Year(metadata.Year<2010),time_bin),:),'fill') 

%             h_s = surf(gridx, gridy,FAC_fit_all_spread,0*gridx);
%             [~, h_c ]= contour3(gridx, gridy,FAC_fit_all_spread,15);

            h_s2 = surf(gridx, gridy,FAC_fit_all_pre2010_mid,0*gridx+1);
            [~, h_c2 ]= contour3(gridx, gridy,FAC_fit_all_pre2010_mid ,15);
        else
            % plotting the core derived FAC
            scatter3(metadata.b_avg(or(metadata.Year>=2010,metadata.T_avg<T_thresh)), ...
            metadata.T_avg(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),...
            metadata.FAC10(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),...
            100,col(discretize(metadata.Year(or(metadata.Year>=2010,metadata.T_avg<T_thresh)),time_bin),:),'fill') 

%             h_s = surf(gridx, gridy,FAC_fit_all_spread,0*gridx);
%             [~, h_c ]= contour3(gridx, gridy,FAC_fit_all_spread,15);

            h_s2 = surf(gridx, gridy,FAC_fit_all_post2010_mid,0*gridx+1);
            [~, h_c2 ]= contour3(gridx, gridy,FAC_fit_all_post2010_mid ,15);
        end
        alpha(h_s2, 0.2)
        alpha(h_s, 0.2)
        h_c.Color = 'b';
        h_c.LineWidth = 1.5;
        h_c2.Color = 'r';
        h_c2.LineWidth = 1.5;
        shading interp
        colormap hsv
%         zlim([0 7])
    xlim([-2 2500])
%     ylim([min(metadata.T_avg)-1 max(metadata.T_avg)+4])
    box on
    xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;
    ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    zlabel('FAC_{10} (m)','Interpreter','tex')
    view(gca,[-103 15.4]);
end

%%  Applying on the T_avg map 

% ind2 =  and(b_raster<accum_thresh,~isnan(T_raster));
% T_raster(~ind2) = NaN;
b_raster(isnan(T_raster)) = NaN;
b_raster(T_raster<-60) = NaN;
T_raster(T_raster<-60) = NaN;
   
FAC_10_map = {};
            h_CA_bis= {};
            h_WHA_bis= {};
            h_WLA_bis= {};
    % pre 2010
    FAC_10_map{1} = max(zeros(size(T_raster)), FAC_bav_mid_past(b_raster,T_raster));
    %post 2010
    FAC_10_map{2} = max(zeros(size(T_raster)), FAC_bav_mid_pres(b_raster,T_raster));
    
%     FAC_10_map{2,2} = max(zeros(size(T_raster)), FAC_bav_spread(b_raster,T_raster));
    
    FAC_10_map{1}(isnan(T_raster))=NaN;
    FAC_10_map{1}(and(T_raster>T_thresh, b_raster>accum_thresh))=NaN;
    FAC_10_map{2}(isnan(T_raster))=NaN;
    
    FAC_10_map{1}(FAC_10_map{1,1}==0)=NaN;
    FAC_10_map{2}(isnan(T_raster))=NaN;
    
    for i =1 : size(metadata,1)
        if or(metadata.T_avg(i) <= T_thresh, metadata.b_avg(i) > accum_thresh)
            metadata.FAC10_predict(i) = FAC_bav_mid_pres(metadata.b_avg(i),metadata.T_avg(i));
        elseif metadata.b_avg(i) <= accum_thresh
            if metadata.Year(i) <= 2010
                metadata.FAC10_predict(i) =  FAC_bav_mid_past(metadata.b_avg(i),metadata.T_avg(i));
            else
                metadata.FAC10_predict(i) = FAC_bav_mid_pres(metadata.b_avg(i),metadata.T_avg(i));
            end
        end
    end
    metadata.residual = metadata.FAC10_predict-metadata.FAC10;

    if source_temp_accum == 1
        FAC_10_map_Box13 = FAC_10_map;
        FAC_10_spread_Box13 = FAC_bav_spread(b_raster,T_raster);
        metadata_Box13 = metadata;
        R_Box13 = R;
        save ('./Output/Box13_result.mat','FAC_10_map_Box13','metadata_Box13','R_Box13','FAC_10_spread_Box13') ;
        load('./Output/MAR_result.mat')
    else
        FAC_10_map_MAR = FAC_10_map;
        metadata_MAR = metadata;
        FAC_10_spread_MAR = FAC_bav_spread(b_raster,T_raster);

        R_MAR = R;
        save('./Output/MAR_result.mat','FAC_10_map_MAR','metadata_MAR','R_MAR','FAC_10_spread_MAR') 
        load('./Output/Box13_result.mat') 
    end
        metadata_save = metadata;
        I_box = geotiffinfo('mean_temperature_box13_3413_firn.tif'); 
        [x_map,y_map]=pixcenters(I_box);
        XX_box = repmat(x_map,length(y_map),1);
        YY_box = repmat(y_map',1,length(x_map));
        
        I_MAR = geotiffinfo('T_avg_1979-2014_MAR_3413_firn.tif'); 
        [x_map,y_map]=pixcenters(I_MAR);
        XX_mar = repmat(x_map,length(y_map),1);
        YY_mar = repmat(y_map',1,length(x_map));
        
        for i = 1:2
                FAC_10_map_MAR_resampled{i} = ...
                    griddata(XX_mar,YY_mar,FAC_10_map_MAR{i},XX_box,YY_box);
                FAC_10_spread_MAR_resampled = ...
                    griddata(XX_mar,YY_mar,FAC_10_spread_MAR,XX_box,YY_box);
        end
        
if source_temp_accum==1
    XX = XX_box;
    YY = YY_box;
else
    XX = XX_mar;
    YY = YY_mar;
end
format BANK

%% Plotting
        f=figure('Units','Normalized','outerposition',[0 0 1 0.8]);
        ha =tight_subplot(1,4,-0.15,[0.03 0.17],[0.01 0.11]);

                set(f,'CurrentAxes',ha(1))
                hold on

                PlotBackground(GL,Ice,Firn);

                temp = FAC_10_map{1};
                temp(isnan(DSA)) = NaN;
                h_map = pcolor(XX,YY,temp);
                h_map.LineStyle = 'none';

                colbar =contourcmap('toh',0:0.5:6,'colorbar','on');
                ylabel(colbar,'FAC_{10} (m)','Interpreter','tex');
                   colbar.Position(2) = 2;

                [~, ~] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',1);

                met = metadata(metadata.T_avg<=T_thresh,:);
                plot(met.X, met.Y ,'ok','MarkerSize',5,'MarkerFaceColor','w'); 

                xlim(1.0e+05 *[-6.2164    8.4827])
                ylim(1.0e+06 *[-3.3439   -0.6732])

                h_title = title('         DSA\newline    1953-2017','Interpreter','tex');
                h_title.FontSize = 14;
                h_title.Units = 'Normalized';
                h_title.Position = [0.45          1.01             0];
                set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
                count = count+1;
                daspect( [1 1 1])
                

            set(f,'CurrentAxes',ha(2))
                hold on

                PlotBackground(GL,Ice,Firn);
                temp = FAC_10_map{1};
                temp(isnan(LAWSA)) = NaN;
                h_map = pcolor(XX,YY,temp);
                h_map.LineStyle = 'none';

                colbar =contourcmap('toh',0:0.5:6,'colorbar','on');
                ylabel(colbar,'FAC_{10} (m)','Interpreter','tex');
                   colbar.Position(2) = 2;

                [~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',1);
                [~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',1);
                [~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',1); 

                met = metadata(and(ismember(metadata.Year,time_bin(6)+0.5:time_bin(7)+0.5),metadata.b_avg<accum_thresh),:);
                met = met(met.T_avg>T_thresh,:);
                plot(met.X, met.Y ,'ok','MarkerSize',5,'MarkerFaceColor','w'); 

                xlim(1.0e+05 *[-6.2164    8.4827])
                ylim(1.0e+06 *[-3.3439   -0.6732])

                h_title = title('          LAWSA\newline         1997-2008','Interpreter','tex');
                h_title.FontSize = 14;
                                h_title.Units = 'Normalized';
                h_title.Position = [0.4          1.01             0];
                set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
                count = count+1;
                daspect( [1 1 1])

     set(f,'CurrentAxes',ha(3))
                hold on
              PlotBackground(GL,Ice,Firn);
                temp = FAC_10_map{2};
                temp(~or(~isnan(LAWSA),~isnan(HAWSA))) = NaN;
                h_map = pcolor(XX,YY,temp);
                h_map.LineStyle = 'none';
                
                colbar =contourcmap('toh',0:0.5:6,'colorbar','on');
                ylabel(colbar,'FAC_{10} (m^3 m^{-2})','Interpreter','tex');
                colbar.Position(1) = 0.81;
                colbar.FontSize = 20;

                [~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',1);
                [~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',1);
                [~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',1); 

                met = metadata(or(and(ismember(metadata.Year,2012:2017),metadata.T_avg>T_thresh),...
                    metadata.b_avg>accum_thresh),:);
                plot(met.X, met.Y ,'ok','MarkerSize',5,'MarkerFaceColor','w'); 

                xlim(1.0e+05 *[-6.2164    8.4827])
                ylim(1.0e+06 *[-3.3439   -0.6732])

                h_title = title(['LAWSA 2011-2017',...
                    '\newline             & \newlineHAWSA 2010-2017'],...
                    'Interpreter','tex');
                h_title.FontSize = 14;
                set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
                daspect( [1 1 1])
               
annotation(f,'textbox',...
    [0.08 0.88 0.04 0.05],...
    'String','a)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(f,'textbox',...
    [0.27 0.88 0.04 0.05],...
    'String','b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(f,'textbox',...
    [0.43 0.91 0.04 0.05],...
    'String','c)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');

% Mapping evolution of pore space

set(f,'CurrentAxes',ha(4))
    ps_change_10m = FAC_10_map{2}-FAC_10_map{1};
%     ps_change_10m(FAC_10_map{1,1} == 0) = NaN;
        hold on
            PlotBackground(GL,Ice,Firn);
disp('Spatial average of FAC10 change')   
temp2 = ps_change_10m;
temp2(isnan(LAWSA)) = NaN;
disp(nanmean(nanmean(temp2)))

disp('Max FAC10 change')
disp(min(min(temp2)))
    h_map = pcolor(XX,YY,ps_change_10m);
    h_map.LineStyle = 'none';
    daspect( [1 1 1])

    colbar2 =contourcmap('tej2',-3.6:0.2:0,'colorbar','on');
    ylabel(colbar2,'\Delta FAC_{10} (m^3 m^{-2})','Interpreter','tex');
colbar2.FontSize = 20;
    xlim(1.0e+05 *[-4    1.5])
    ylim(1.0e+06 *[-3.   -1.4])
    h_title = title('                LAWSA \newline2011-2017 minus 1997-2008','Interpreter','tex');
                                h_title.Units = 'Normalized';
                h_title.Position = [0.45          1.03             0];
    box on
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',2);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',2);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',2);

            set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'layer','top')
annotation(f,'textbox',...
    [0.66 0.92 0.04 0.05],...
    'String','d)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
ha(4).Position(1) = 0.62;
colbar.Position(1)=0.62;
colbar2.Position(1)=0.83;
colormap(ha(1),'toh')
colormap(ha(2),'toh')
colormap(ha(3),'toh')
colormap(colbar,'toh')
colormap(ha(4),'tej2')
            for i=2:2:size(colbar2.YTickLabel,1)
            colbar2.YTickLabel(i,:) = '    ';
            end
    print(f,sprintf('Output/FAC_map_%s',orig),'-dtiff')

%% Uncertainty analysis LAWSA

disp('Sensitivity analysis in the LAWSA')
close all
plotting =0;
SA_LAWSA_post2010 = {};
SA_LAWSA_pre2010 = {};

for ii = 1:1000
    if mod(ii,100) == 0
        fprintf('%i %%\n',ii/10)
    end
metadata = metadata_save;
ind = or(metadata.b_avg>accum_thresh,metadata.T_avg<T_thresh-1);
metadata(ind,:) = [];
    ind_binned = discretize(metadata.Year,time_bin);

    ind = find(ind_binned==6);
        x = metadata.b_avg(ind);
        
    ind_remove  = unidrnd(length(x),4,1);
    metadata(ind(ind_remove),:) = [];
    
    ind_binned = discretize(metadata.Year,time_bin);
    ind = find(ind_binned==7);
        x = metadata.b_avg(ind);

    ind_remove  = unidrnd(length(x),4,1);
    metadata(ind(ind_remove),:) = [];
    
    ind_binned = discretize(metadata.Year,time_bin);

    % 3D exploration of accumulation and temperature control
    
    i =  6;
    x = metadata.b_avg(ind_binned==i);
    y = metadata.T_avg(ind_binned==i);
    z = metadata.FAC10(ind_binned==i);

    
    xx = linspace(min(metadata.b_avg),accum_thresh);
    yy = linspace(T_thresh,max(metadata.T_avg)+10);

    Eqn = sprintf(['max(0, a*x + ' ...
        'b*min(y+%i,0) + c*max(y+%i,0) + ' ...
        'd*min(y+%i,0) + e*max(y+%i,0) + f)'],...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*2/3),...
        -(min(y)+ (max(y) - min(y))*2/3));
    SA_LAWSA_pre2010{ii} = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[Inf 0    0   0   0     Inf],...
        'StartPoint',[0.45414        0.5029       0.72089      0.54005      0.15411      0.12192]);

    count = 0;

while sqrt(mean((SA_LAWSA_pre2010{ii}(x,y)-z).^2))>1
     SA_LAWSA_pre2010{ii} = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[Inf 0    0   0   0     Inf]);
    count = count+1;
end

[gridx, gridy] = meshgrid(xx,yy);
    
    if plotting
        f = figure('Visible', vis, 'outerposition',[1 0 20 20]);
        ha = tight_subplot(1,1,0.1,0.25,0.25);
        set(f,'CurrentAxes',ha(1))
        hold on

        xlim([min(metadata.b_avg)-2 max(metadata.b_avg)+10])
        ylim([min(metadata.T_avg)-1 max(metadata.T_avg)+4])

             box on
        h_xlab = xlabel('$\mathrm{\overline{\dot{b}}  \: (mm w.eq. yr^{-1})}$','Interpreter','latex') ;

        ylabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
        zlabel('FAC_{10} (m)','Interpreter','tex')
        scatter3(x,y, z,...
            100,col(6,:),'fill') 
        h_s = surf(gridx, gridy,SA_LAWSA_pre2010{ii}(gridx, gridy) ,0*gridy+0.5);
        alpha(h_s, 0.2)

        [~, h_c ]= contour3(gridx, gridy, SA_LAWSA_pre2010{ii}(gridx, gridy),15);
        h_c.Color = 'b';
        h_c.LineWidth = 1.5;
        view(gca,[-103 15.4]);

end
    i =  7;
    x = metadata.b_avg(ind_binned==i);

    y = metadata.T_avg(ind_binned==i);
    z = metadata.FAC10(ind_binned==i);

    Eqn = sprintf(['max(0, a*x + ' ...
        'b*min(y+%0.2f,0) + c*max(y+%0.2f,0) + ' ...
        'd*min(y+%0.2f,0) + e*max(y+%0.2f,0) + f)'],...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*1/3),...
        -(min(y)+ (max(y) - min(y))*2/3),...
        -(min(y)+ (max(y) - min(y))*2/3));
    SA_LAWSA_post2010{ii} = fit([x, y],z,Eqn,...
        'Lower',[0     -Inf  -Inf   -Inf  -Inf -Inf],...
        'Upper',[0.003 0    0   0   0     Inf],...
        'StartPoint',[0.45414        0.5029       0.72089      0.54005      0.15411      0.12192]);
    [gridx, gridy] = meshgrid(xx,yy);
    
    if plotting
        scatter3(x,y, z,...
        100,col(7,:),'fill') 
            h_s1 = surf(gridx, gridy, SA_LAWSA_post2010{ii}(gridx, gridy), 0*gridy+3);
            shading interp
            alpha(h_s1, 0.2)
            [~, h_c1 ] = contour3(gridx, gridy, SA_LAWSA_post2010{ii}(gridx, gridy),15);
            h_c1.Color = 'r';
            h_c1.LineWidth = 1.5;
            zlim([0 7])   
               colormap('hsv')

        view(gca,[-103 15.4]);
        pause
    end

end

XX_LAWSA =b_raster;
XX_LAWSA(LAWSA~=1) = NaN;
YY_LAWSA =T_raster;
YY_LAWSA(LAWSA~=1) = NaN;

M = NaN([size(b_raster),length( FAC_fit_WHA_post2010)]);
for ii = 1:length( SA_LAWSA_post2010)
    M(:,:,ii) =SA_LAWSA_post2010{ii}(XX_LAWSA, YY_LAWSA);
end

M_mean1 = mean(M,3);
Uncertainty_post2010_LAWSA = std(M,0,3);

M = NaN([size(b_raster),length( FAC_fit_WHA_pre2010)]);
for ii = 1:length(SA_LAWSA_pre2010)
    M(:,:,ii) =SA_LAWSA_pre2010{ii}(XX_LAWSA, YY_LAWSA);
end

M_mean = mean(M,3);
Uncertainty_pre2010_LAWSA = std(M,0,3);

figure('Visible',vis)
subplot(1,2,1)
surf(XX,YY,Uncertainty_pre2010_LAWSA)
shading interp
colorbar
title('Uncertainty in the LAWSA (1997-2008)')

subplot(1,2,2)
surf(XX,YY,Uncertainty_post2010_LAWSA)
shading interp
colorbar
title('Uncertainty in the LAWSA (2011-2017)')

%% Plotting uncertainty
metadata = metadata_save;
disp('Mapping Uncertainty of pore space')

    disp('Mean residual:')
    fprintf('%0.2f m3 m-3\n',nanmean(metadata.residual))
    
    f=figure('Units','Normalized','outerposition',[0 0 0.7 0.8]);
    ha =tight_subplot(1,3,-0.15,[0.03 0.17],[0.01 0.11]);

set(f,'CurrentAxes',ha(1))
    uncertainty_DSA = CA_save*uncert_DSA ;
    uncertainty_DSA(FAC_10_map{1} == 0) = NaN;
    uncertainty_DSA(HAWSA==1) = NaN; 
    uncertainty_DSA(LAWSA==1) = NaN; 
%     delta_FAC_1(DSA==1) = NaN; 
    uncertainty_DSA(isnan(FAC_10_map{1})) = NaN;

        hold on
            PlotBackground(GL,Ice,Firn);

            h_map = pcolor(XX,YY,uncertainty_DSA./FAC_10_map{1}*100);
            h_map.LineStyle = 'none';
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',2);

            daspect( [1 1 1])
            xlim(1.0e+05 *[-6.2164    8.4827])
            ylim(1.0e+06 *[-3.3439   -0.6732])
            title('     DSA\newline1953-2017','Interpreter','tex')
            colbar =contourcmap('jet',0:5:100,'colorbar','on');
                ylabel(colbar,'Uncertainty on FAC_{10} (m)','Interpreter','tex');
            set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
            for i=1:2:size(colbar.YTickLabel,1)
            colbar.YTickLabel(i,:) = '   ';
            end
            colbar.Position(2) = 1.5;
set(f,'CurrentAxes',ha(2))
    uncertainty_LAWSA_pre2010 = max(0.3, Uncertainty_pre2010_LAWSA*2) ;
%     delta_FAC_2(delta_FAC_2<=  0.3) = 0.299;
    uncertainty_LAWSA_pre2010(HAWSA==1) = NaN; 
    uncertainty_LAWSA_pre2010(DSA==1) = NaN; 
%     uncertainty_LAWSA_pre2010(FAC_10_map{1,1} == 0) = NaN;
    uncertainty_LAWSA_pre2010(isnan(FAC_10_map{1,1})) = NaN;

        hold on
            PlotBackground(GL,Ice,Firn);

            h_map = pcolor(XX,YY,uncertainty_LAWSA_pre2010./FAC_10_map{1}*100);
            h_map.LineStyle = 'none';
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',2);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',2);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',2);
            daspect( [1 1 1])
            xlim(1.0e+05 *[-6.2164    8.4827])
            ylim(1.0e+06 *[-3.3439   -0.6732])
            title(' LAWSA  \newline1997-2008','Interpreter','tex')
            colbar =contourcmap('jet',0:5:100,'colorbar','on');
                ylabel(colbar,'Uncertainty on FAC_{10} (m)','Interpreter','tex');
            colbar.Position(2) = 1.5;
             for i=1:2:size(colbar.YTickLabel,1)
            colbar.YTickLabel(i,:) = '   ';
            end
                        set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')

    set(f,'CurrentAxes',ha(3))
           uncertainty_LAWSA_HAWSA_post2010 = abs(FAC_bav_spread(b_raster,T_raster))/2 + Uncertainty_post2010_LAWSA*2;
           ind_nan = isnan(uncertainty_LAWSA_HAWSA_post2010);
           uncertainty_LAWSA_HAWSA_post2010=  max(0.3,uncertainty_LAWSA_HAWSA_post2010);
           uncertainty_LAWSA_HAWSA_post2010(ind_nan)=NaN;
            uncertainty_LAWSA_HAWSA_post2010(DSA==1) = NaN; 

        hold on
            PlotBackground(GL,Ice,Firn);

            h_map = pcolor(XX,YY,uncertainty_LAWSA_HAWSA_post2010./FAC_10_map{2}*100);
            h_map.LineStyle = 'none';
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',2);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',2);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',2);

            daspect( [1 1 1])
            xlim(1.0e+05 *[-6.2164    8.4827])
            ylim(1.0e+06 *[-3.3439   -0.6732])
            title('LAWSA 2011-2017\newline              &   \newlineHAWSA 2010-2017','Interpreter','tex')
            colbar =contourcmap('jet',0:5:100,'colorbar','on');
            colbar.Position(1) = 0.81;
            for i=2:2:size(colbar.YTickLabel,1)
            colbar.YTickLabel(i,:) = '   ';
            end
            colbar.FontSize = 20;
            ylabel(colbar,'        Relative uncertainty \newline  on predicted FAC_{10} (%)','Interpreter','tex');
            set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
uistack(ha(1),'top')
annotation(f,'textbox',...
    [0.077 0.88 0.04 0.05],...
    'String','a)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(f,'textbox',...
    [0.32 0.88 0.04 0.05],...
    'String','b)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(f,'textbox',...
    [0.53 0.9 0.04 0.05],...
    'String','c)',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off');
colormap(ha(1),'parula')
colormap(ha(2),'parula')
colormap(ha(3),'parula')
colormap(colbar,'parula')

print(f,sprintf('Output/uncertainty_FAC_map_overall_%s',orig),'-dtiff')

%% difference Box MAR            
                f=figure('Visible', vis, 'outerposition',[0 0 30 20]);
    ha =tight_subplot(1,2,0.0001,[0.02 0.14],[0.07 0.3]);
            set(f,'CurrentAxes',ha(1))
    delta_FAC_4 = abs(FAC_10_map_MAR_resampled{1}-FAC_10_map_Box13{1});
    delta_FAC_4(FAC_10_map_Box13{1,1} == 0) = NaN;

    hold on
                       PlotBackground(GL,Ice,Firn);            
            h_map = pcolor(XX_box,YY_box,delta_FAC_4);
            h_map.LineStyle = 'none';
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',3);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',3);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',3);

%             scatter(metadata.X, metadata.Y ,50,...
%                 abs(metadata_MAR.FAC10_predict-metadata_Box13.FAC10_predict),'fill'); 
%             scatter(metadata.X, metadata.Y ,50,RGB('light light gray'));

            daspect( [1 1 1])
ha(1).Position(4) = 0.5;
ha(1).Position(1) = 0.22;
ha(1).Position(2) = 0.05;
            xlim(1.0e+05 *[-6.2164    1.1])
            ylim(1.0e+06 *[-3   -1.3])
            box on
            title('LAWSA 1997-2008','Interpreter','tex')
            colbar =contourcmap('hsv',0:0.3:2.5,'colorbar','on');
                ylabel(colbar,' ','Interpreter','tex');
            set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'layer','top')
            colbar.Position(2) = 1.5;

set(f,'CurrentAxes',ha(2))
ha(2).Position(1) = 0.4;
    delta_FAC_5 = abs(FAC_10_map_MAR_resampled{2}-FAC_10_map_Box13{2});
    delta_FAC_5(FAC_10_map_Box13{1,1} == 0) = NaN;

    hold on
                        PlotBackground(GL,Ice,Firn);            
            h_map = pcolor(XX_box,YY_box,delta_FAC_5);
            h_map.LineStyle = 'none';
[~, h_CA] = contour(XX,YY,CA_save, [1 1],'Color','k', 'LineWidth',3);
[~, h_WHA] = contour(XX,YY,WHA_save, [0.5 0.5],'Color','k', 'LineWidth',3);
[~, h_WLA] = contour(XX,YY,WLA_save, [0.5 0.5],'Color','k', 'LineWidth',3);

daspect( [1 1 1])
            xlim(1.0e+05 *[-6.2164    8.4827])
            ylim(1.0e+06 *[-3.3439   -0.6732])
            title('DSA (1953-2017)\newlineLAWSA (2011-2017)\newline HAWSA (2010-2017)','Interpreter','tex')

            colbar =contourcmap('jet2',0:0.3:2.5,'colorbar','on');
            ylabel(colbar,'Absolute difference between Box13-derived  \newline           and MAR-derived FAC_{10}(m^3 m^{^-2})','Interpreter','tex');
            set(gca,'XTickLabel','','YTickLabel','','XColor','w','YColor','w')
uistack(ha,'top')
            colbar.FontSize = 20;

print(f,('Output/uncertainty_FAC_map_MAR_vs_Box13'),'-dtiff')
 
%% Screen printing
% average values

temp = FAC_10_map{1};
temp(isnan(DSA)) = NaN;
temp2 = uncertainty_DSA;
temp2(isnan(DSA)) = NaN;

disp('Spatial average of FAC10 in DSA')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp)))
disp('average uncertainty')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp2)))
fprintf('%0.2f %%\n\n',nanmean(nanmean(temp2))/nanmean(nanmean(temp))*100)

temp = FAC_10_map{1};
temp(isnan(LAWSA)) = NaN;
temp2 = uncertainty_LAWSA_pre2010;
temp2(isnan(LAWSA)) = NaN;

disp('Spatial average of FAC10 in LAWSA')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp)))
disp('average uncertainty')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp2)))
fprintf('%0.2f %%\n\n',nanmean(nanmean(temp2))/nanmean(nanmean(temp))*100)

temp = FAC_10_map{2};
temp(isnan(LAWSA)) = NaN;
temp2 = Uncertainty_post2010_LAWSA;
temp2(isnan(LAWSA)) = NaN;

disp('Spatial average of FAC10 in LAWSA post2010')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp)))
disp('average uncertainty')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp2)))
fprintf('%0.2f %%\n\n',nanmean(nanmean(temp2))/nanmean(nanmean(temp))*100)

temp = FAC_10_map{2};
temp(isnan(HAWSA)) = NaN;
temp2 = uncertainty_LAWSA_HAWSA_post2010;
temp2(isnan(HAWSA)) = NaN;

disp('Spatial average of FAC10 in HAWSA')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp)))
disp('average uncertainty')
fprintf('%0.2f m3 m-3\n',nanmean(nanmean(temp2)))
fprintf('%0.2f %%\n\n',nanmean(nanmean(temp2))/nanmean(nanmean(temp))*100)

%% Summing for each area
if source_temp_accum == 1
    R = R_Box13;
end
SpatialIntegration = @(M,R) (nansum(M* R.CellExtentInWorldX *R.CellExtentInWorldY))/1000000000;

% disp('Mid estimate all pre2010 (km3):')
% disp([SpatialIntegration(FAC_10_map{1}, R),...
%     SpatialIntegration(FAC_10_map{1}, R)-SpatialIntegration(FAC_10_map{1}+0.3, R),...
%     SpatialIntegration(FAC_10_map{1}, R)-SpatialIntegration(FAC_10_map{1}-0.3, R)]);

disp('Mid estimate DSA  (km3):')
disp([SpatialIntegration(FAC_10_map{1}(DSA==1), R),SpatialIntegration(uncertainty_DSA(DSA==1), R)]);
disp('Mid estimate LAWSA pre2010 (km3):')
disp([SpatialIntegration(FAC_10_map{1}(LAWSA==1), R), SpatialIntegration(uncertainty_LAWSA_pre2010(LAWSA==1), R)]);
disp('Mid estimate LAWSA post2010 (km3):')
disp([SpatialIntegration(FAC_10_map{2}(LAWSA==1), R),SpatialIntegration(uncertainty_LAWSA_HAWSA_post2010(LAWSA==1), R)]);
disp('Mid estimate HAWSA post2010 (km3):')
disp([SpatialIntegration(FAC_10_map{2}(HAWSA==1), R),SpatialIntegration(uncertainty_LAWSA_HAWSA_post2010(HAWSA==1), R)]);
nanmedian(FAC_10_map{1}(LAWSA==1) - FAC_10_map{2}(LAWSA==1))

disp('Total FAC post2010 (km3):')
TotalFAC = (SpatialIntegration(FAC_10_map{2}(HAWSA==1), R) + SpatialIntegration(FAC_10_map{2}(LAWSA==1), R) + SpatialIntegration(FAC_10_map{1}(DSA==1), R));
uncert_TotalFAC = SpatialIntegration(uncertainty_DSA(DSA==1), R) + SpatialIntegration(uncertainty_LAWSA_HAWSA_post2010(LAWSA==1), R) +SpatialIntegration(uncertainty_LAWSA_HAWSA_post2010(HAWSA==1), R);
disp([TotalFAC uncert_TotalFAC]);
disp('assuming that it was filled with ice (GT):')
disp([TotalFAC uncert_TotalFAC]* 873 * 10^9 / 10^12) ;  
disp('assuming that it was filled with ice (mm SLE):')
disp([TotalFAC uncert_TotalFAC]* 873 * 10^9 / 10^12/ 361) ;
  
disp('Loss LAWSA (km3):')
loss_LAWSA = SpatialIntegration(FAC_10_map{1}(LAWSA==1), R)-SpatialIntegration(FAC_10_map{2}(LAWSA==1), R);
uncert_loss_LAWSA = SpatialIntegration(uncertainty_LAWSA_pre2010(LAWSA==1), R) + SpatialIntegration(uncertainty_LAWSA_HAWSA_post2010(LAWSA==1), R);
disp([loss_LAWSA uncert_loss_LAWSA]);

disp('assuming that it was filled with ice (GT):')
disp([loss_LAWSA uncert_loss_LAWSA]* 873 * 10^9 / 10^12) ;
    
disp('assuming that it was filled with ice (mm SLE):')
disp([loss_LAWSA uncert_loss_LAWSA]* 873 * 10^9 / 10^12/ 361) ;

%% Printing to files
disp('')

disp('Printing FAC10 map')

% DSA
    temp = FAC_10_map{1};
    temp(isnan(DSA)) = NaN;
    mkdir('Output/FAC10 maps')
    filename= sprintf('Output/FAC10 maps/FAC10_DSA_%s.tif', name_source);
    geotiffwrite(filename, temp,R,'CoordRefSysCode',3413)  

        filename= sprintf('Output/FAC10 maps/FAC10_uncertainty_DSA_%s.tif', name_source);
        temp = uncertainty_DSA;
        temp(isnan(temp))=NaN;
        geotiffwrite(filename, uncertainty_DSA,R,'CoordRefSysCode',3413)  
    
% LAWSA
    %pre
    temp1 = FAC_10_map{1};
    temp1(isnan(LAWSA)) = NaN;
    
    filename= sprintf('Output/FAC10 maps/FAC10_LAWSA_pre2010_%s.tif', name_source);
    geotiffwrite(filename, temp1,R,'CoordRefSysCode',3413)

        filename= sprintf('Output/FAC10 maps/FAC10_uncertainty_LAWSA_pre2010_%s.tif', name_source);
        
        temp = uncertainty_LAWSA_pre2010;
        temp(isnan(temp))=NaN;
    geotiffwrite(filename, temp,R,'CoordRefSysCode',3413)

    %post
    temp2 = FAC_10_map{2};
    temp2(isnan(LAWSA)) = NaN;
    
    filename= sprintf('Output/FAC10 maps/FAC10_LAWSA_post2011_%s.tif', name_source);
    geotiffwrite(filename, temp2,R,'CoordRefSysCode',3413)  

        temp = uncertainty_LAWSA_HAWSA_post2010;
        temp(isnan(LAWSA)) = NaN;

        filename= sprintf('Output/FAC10 maps/FAC10_uncertainty_LAWSA_post2011_%s.tif', name_source);
        geotiffwrite(filename, temp,R,'CoordRefSysCode',3413)
        
    %change
        filename= sprintf('Output/FAC10 maps/FAC10_change_LAWSA_%s.tif', name_source);
        temp = temp1-temp2;
        temp(isnan(temp))=NaN;
        geotiffwrite(filename, temp1-temp2,R,'CoordRefSysCode',3413)  

%HAWSA
    temp = FAC_10_map{2};
    temp(isnan(HAWSA)) = NaN;
    
    filename= sprintf('Output/FAC10 maps/FAC10_HAWSA_%s.tif', name_source);
    geotiffwrite(filename, temp,R,'CoordRefSysCode',3413) 

        temp = uncertainty_LAWSA_HAWSA_post2010;
        temp(isnan(HAWSA)) = NaN;

        filename= sprintf('Output/FAC10 maps/FAC10_uncertainty_HAWSA_%s.tif', name_source);
        geotiffwrite(filename, temp,R,'CoordRefSysCode',3413) 
       
%% Analysis of the temperature thresholds
figure

subplot(2,1,1)
% for Box13
scatter(metadata_Box13.T_avg, metadata_Box13.FAC10)
bins = floor(min(metadata_Box13.T_avg)): 2 : max(metadata_Box13.T_avg);
std_bin = 0*bins(1:end-1);

for i = 1:length(bins)-1
    ind = and( metadata_Box13.T_avg>=bins(i), metadata_Box13.T_avg<bins(i+1));
    std_bin(i) = std(metadata_Box13.FAC10(ind));
end

hold on
stairs(bins(1:end-1), std_bin,'LineWidth',2)
plot(bins([1 end]), [1 1]*0.3,'--k')
h_leg1 = legend('FAC10 measurements','Standard deviation of FAC10\newline measurements in 1K-wide bins',...
    'Criteria defining the temperature \newline threshold: 0.3 m3 m-2','Location','NorthOutside','Interpreter','tex') ;
h_leg1.Units = 'normalized';
xlabel('Long-term average temperature (degC)')
ylabel('FAC_10 (m3 m-2)')
box on
title('Temperature threshold between DSA and LAWSA \newline                              for Box13 product','Interpreter','tex')

% for MAR
bins = floor(min(metadata_MAR.T_avg))-1: 1 : max(metadata_MAR.T_avg);
std_bin = 0*bins(1:end-1);

for i = 1:length(bins)-1
    ind = and( metadata_MAR.T_avg>=bins(i), metadata_MAR.T_avg<bins(i+1));
    std_bin(i) = std(metadata_MAR.FAC10(ind));
end

subplot(2,1,2)
scatter(metadata_MAR.T_avg, metadata_MAR.FAC10)
hold on
stairs(bins(1:end-1), std_bin,'LineWidth',2)
plot(bins([1 end]), [1 1]*0.3,'--k')
h_leg = legend('FAC10 measurements','Standard deviation of FAC10\newline measurements in 1K-wide bins',...
    'Criteria defining the temperature \newline threshold: 0.3 m3 m-2','Location','NorthOutside','Interpreter','tex') ;
h_leg.Units = 'normalized';
h_leg.Parent = gcf;
h_leg1.Parent = gcf;
h_leg.Position = [0.72 0.8 0.26 0.19];
h_leg1.Position = h_leg.Position;
xlabel('Long-term average temperature (degC)')
ylabel('FAC_10 (m3 m-2)')
box on
title('  \newline for MAR products','Interpreter','tex')

%% Slope-break between DSA and HAWSA
figure
ind_jaune = metadata.T_avg<T_thresh;
ind_vert = and(metadata.T_avg>T_thresh,metadata.b_avg>accum_thresh);
ind_rouge = and(metadata.T_avg>T_thresh,metadata.b_avg<=accum_thresh);
x1 = metadata.T_avg(ind_jaune);
y1=metadata.FAC10(ind_jaune);
x2 = metadata.T_avg(ind_vert);
y2=metadata.FAC10(ind_vert);
hold on
    % ha(6).Position(2) = 0.28;
    scatter(x1,y1,20,'k','fill')
    scatter(x2,y2,20,'k','fill')
    [lm, ah] = Plotlm(x1,y1,'Unit','m^3 m^{-2} K^{-1}');
    ah.LineStyle = 'none';
    [lm, ah] = Plotlm(x2,y2,'Unit','m^3 m^{-2} K^{-1}');
        ah.LineStyle = 'none';

    hold on
    % plot([T_thresh T_thresh],[0 7],'--k','LineWidth',2)
    xlabel('$\mathrm{\overline{T_a}}  \:(^oC)$','Interpreter','latex')
    box on

    axis tight
    xlimits = get(gca,'XLim');
    pp1 = patch([xlimits(1) xlimits(1) T_thresh T_thresh ],...
        [0 7 7 0],...
        jaune);
    pp2 = patch([T_thresh T_thresh xlimits(2)   xlimits(2)],...
        [0 7 7 0],...
        vert);
    uistack(pp1,'bottom')
    uistack(pp2,'bottom')

    h_ylab = ylabel('FAC_{10} (m)','Interpreter','tex');
set(gca,'Ticklength',[0.08 0.16]/2.5,'Ticklength',[0.08 0.16]/2,...
    'XMinorTick','on','YMinorTick','on','XAxisLocation','bottom','YAxisLocation','left','layer','top')
axis tight

ylim([0 7])
title('Slope break between DSA and HAWSA')

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
            M = table(metadata.CoreNumber(ind_year(nearby_cores)), ...
                metadata.Name(ind_year(nearby_cores)), ...
                metadata.Year(ind_year(nearby_cores)), ...
                metadata.FAC10(ind_year(nearby_cores)),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*mean(metadata.FAC10(ind_year(nearby_cores))),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*std(metadata.FAC10(ind_year(nearby_cores))),...
                ones(size(metadata.FAC10(ind_year(nearby_cores))))*numel(ind_year(nearby_cores)),...
                metadata.Citation(ind_year(nearby_cores)),'VariableNames',...
                {'CoreNumber','Name','Year','FAC_10m','mean_FAC','std_FAC','NumCores','Citation'});
            
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
ylabel('Standard deviation \newlineof FAC\_10 (m^3 / m^2)','Interpreter','tex')

print(f,'Output/Standard_Deviation.tif','-dtiff')

writetable(out_table,'./Output/Uncertainty_table.csv','Delimiter',';');

