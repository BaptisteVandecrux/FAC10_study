% FAC_RCM_comparison.m
% This script compares the Firn Air Content (FAC) observed in firn cores
% and calculated by three regional climate models (RCM): HIRHAM5, RACMO2.3p2 and
% MARv3.9. The RCM FAC data is not publicly available so to run this code,
% additional files are needed. Contact me or RCM managers to obtain the
% appropriate datasets.
% Baptiste Vandecrux
% Technical University of Denmark
% Geological Survey of Denmark and Greenland
% b.vandecrux@gmail.com
% ========================================================================

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

%% Loading FAC10 metadata

filename = '.\Output\metadata_MAR.csv';
delimiter = ';';
startRow = 2;

formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,5,6,7,8,9,11,12,13,14,16]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
try
    dates{2} = datetime(dataArray{2}, 'Format', 'dd-MMM-yyyy', 'InputFormat', 'dd-MMM-yyyy');
catch
    try
        dataArray{2} = cellfun(@(x) x(2:end-1), dataArray{2}, 'UniformOutput', false);
        dates{2} = datetime(dataArray{2}, 'Format', 'dd-MMM-yyyy', 'InputFormat', 'dd-MMM-yyyy');
    catch
        dates{2} = repmat(datetime([NaN NaN NaN]), size(dataArray{2}));
    end
end

anyBlankDates = cellfun(@isempty, dataArray{2});
anyInvalidDates = isnan(dates{2}.Hour) - anyBlankDates;
dates = dates(:,2);

rawNumericColumns = raw(:, [1,3,5,6,7,8,9,11,12,13,14,16]);
rawCellColumns = raw(:, [4,10,15,17]);
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

metadata = table;
metadata.Year = cell2mat(rawNumericColumns(:, 1));
metadata.Date = dates{:, 1};
metadata.CoreNumber = cell2mat(rawNumericColumns(:, 2));
metadata.Name = rawCellColumns(:, 1);
metadata.Latitude = cell2mat(rawNumericColumns(:, 3));
metadata.Longitude = cell2mat(rawNumericColumns(:, 4));
metadata.Elevation = cell2mat(rawNumericColumns(:, 5));
metadata.T_avg = cell2mat(rawNumericColumns(:, 6));
metadata.c_avg = cell2mat(rawNumericColumns(:, 7));
metadata.Densities = rawCellColumns(:, 2);
metadata.DepthMax = cell2mat(rawNumericColumns(:, 8));
metadata.FAC10 = cell2mat(rawNumericColumns(:, 9));
metadata.FAC10_HL = cell2mat(rawNumericColumns(:, 10));
metadata.SnowDepth = cell2mat(rawNumericColumns(:, 11));
metadata.Citation = rawCellColumns(:, 3);
metadata.FACpco = cell2mat(rawNumericColumns(:, 12));
metadata.VarName17 = rawCellColumns(:, 4);

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns rawCellColumns R;

metadata.FACtot = metadata.FACpco;
metadata.FACpco(metadata.DepthMax<100) = NaN;

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

            c_map.lon = lon_org(:);
            T_map.lat = lat_org(:);
            c_map.lon = lon_org(:);
            T_map.lat = lat_org(:);

            c_map.c_avg = mean (acc(:,:,131:end),3);
            T_map.T_avg = mean (Temperature(:,:,131:end),3);
            
            dlmwrite('./Output/Temp accum maps/mean_temperature_box13.csv', [lon_org(:) lat_org(:) T_map.T_avg(:)],'Delimiter',';');
            dlmwrite('./Output/Temp accum maps/mean_accumulation_box13.csv', [lon_org(:) lat_org(:) c_map.c_avg(:)],'Delimiter',';');
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
        lon_org =  reshape(T_map.lon,561,301);
        lat_org =  reshape(T_map.lat,561,301);
% Loading firn grid
[FirnArea, R_firn]= geotiffread('FirnLayer_2000_2017_final_4326.tif');
I = geotiffinfo('FirnLayer_2000_2017_final_4326.tif');
FirnArea(isnan(FirnArea)) = 0;
FirnArea(FirnArea>1) = 0;
[x_firn,y_firn]=pixcenters(I);
[XX_firn,YY_firn] = meshgrid(x_firn,y_firn);

is_firn = interp2(XX_firn,YY_firn,FirnArea,lon_org(:), lat_org(:),'nearest');
is_firn(isnan(is_firn)) = 0;

T_map.T_avg(find(~is_firn)) = NaN;
c_map.c_avg(find(~is_firn)) = NaN;

        ind_firn = ~isnan(T_map.T_avg);
 ind_DSA_org = and(ind_firn, T_map.T_avg<T_thresh);
 ind_LAPA_org = and(ind_firn, ...
     and(T_map.T_avg>=T_thresh,c_map.c_avg < accum_thresh));
 ind_HAPA_org = and(ind_firn, ...
     and(T_map.T_avg>=T_thresh, c_map.c_avg >= accum_thresh));

%% Extracting from HH 
%Loading
disp('Loading HH')
tic
filename = '.\Input\RCM FAC\HIRHAM5\HH_FAC.csv';
if exist(filename) ~= 2 
    error(['The file %s is not publicly available.\n' ...
        'Contact Baptiste Vandecrux (b.vandecrux@gmail.com) or RCM managers to obtain RCM FACs.'],...
        filename);
end
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

fclose(fileID);
HHFAC = table(dataArray{1:end-1}, 'VariableNames', ...
    {'Date','CoreNumber','Lat','Lon',...
    'rho_10_LIN', 'rho_100_LIN', 'rho_10_MOD', 'rho_100_MOD',...
    'FAC10_LIN',  'FAC100_LIN',  'FAC10_MOD',  'FAC100_MOD'});
clearvars filename delimiter formatSpec fileID dataArray ans;
metadata.FAC10_HH_LIN = HHFAC.FAC10_LIN; 
metadata.FAC100_HH_LIN = HHFAC.FAC100_LIN; 
metadata.FAC10_HH_MOD = HHFAC.FAC10_MOD; 
metadata.FAC100_HH_MOD = HHFAC.FAC100_MOD; 
toc

%% Extracting from MAR 
disp('loading MAR')
tic
path = '.\Input\RCM FAC\MARv3.9\Monthly subsurface grids';
clearvars dir
file_list = dir(path);
file_list(1:2) = [];
rho10_MAR = NaN(96,179,12*size(file_list,1));
time_MAR = NaN(1,12*size(file_list,1));
ind_DSA = metadata.T_avg<T_thresh;
ind_LAPA = and(metadata.T_avg>T_thresh, metadata.c_avg<accum_thresh);
ind_HAPA = and(metadata.T_avg>T_thresh, metadata.c_avg>accum_thresh);
    if isempty(file_list) 
        error(['The content of folder %s is not publicly available.\n' ...
            'Contact Baptiste Vandecrux (b.vandecrux@gmail.com) or RCM managers to obtain RCM FACs.'],...
            path);
    end
    
%     figure
for k =1 : size(file_list,1)
    namefile = sprintf('%s\\%s',path,file_list(k).name);

    finfo = ncinfo(namefile);
%     names={finfo.Variables.Name};
    if k == 1
            names={'LON','LAT','RO_10M_AVE'};
    else
            names={'RO_10M_AVE'};
    end
    for i= 1:length(names) 
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end
    rho10_MAR(:,:,((k-1)*12+1):(k*12)) = RO_10M_AVE;
%     pcolor(LON,LAT,RO_10M_AVE(:,:,6));
%     shading interp
    time_MAR(((k-1)*12+1):(k*12)) = datenum((1957+k),1:12,1);
end
lat_MAR = double(LAT);
lon_MAR = double(LON);
ps = @(m,rho) max(0,m.*(1./rho -1/917));
FAC10_MAR = ps(rho10_MAR*10,rho10_MAR);

ind_time = dsearchn(time_MAR', datenum(metadata.Date));
[i_lat, i_lon] = ind2sub(size(LAT), dsearchn([LAT(:) LON(:)], [metadata.Latitude metadata.Longitude]));
metadata.FAC10_MAR = metadata.FAC10*NaN;
for k = 1: length(metadata.FAC10)
    metadata.FAC10_MAR(k) = FAC10_MAR(i_lat(k),i_lon(k),ind_time(k));
end
toc

%% Extracting from RACMO5
 disp('RACMO5')
tic
dir = '.\Input\RCM FAC\RACMO2.3p2\';
flnm1 = [dir 'FDM_FirnAir_10m_FGRN055_rerun_1960-2017.nc'];
    if exist(flnm1) ~= 2 
        error(['The file %s is not publicly available.\n' ...
            'Contact Baptiste Vandecrux (b.vandecrux@gmail.com) or RCM managers to obtain RCM FACs.'],...
            flnm1);
    end
FAC10_RACMO5 = ncread(flnm1, 'FirnAir_10m'); 
time_RACMO5 = ncread(flnm1, 'time');
lat_RACMO5 = ncread(flnm1, 'lat'); 
lon_RACMO5 = ncread(flnm1, 'lon');

flnm1 = [dir 'FDM_FirnAir_FGRN055_rerun_1960-2017.nc'];
    if exist(flnm1) ~= 2 
        error(['The file %s is not publicly available.\n' ...
            'Contact Baptiste Vandecrux (b.vandecrux@gmail.com) or RCM managers to obtain RCM FACs.'],...
            flnm1);
    end
FAC100_RACMO5 = ncread(flnm1, 'FirnAir'); 
 
[FirnArea, R_firn]= geotiffread('FirnLayer_2000_2017_final_4326.tif');
I = geotiffinfo('FirnLayer_2000_2017_final_4326.tif');
FirnArea(isnan(FirnArea)) = 0;
FirnArea(FirnArea>1) = 0;
[x_firn,y_firn]=pixcenters(I);
[XX_firn,YY_firn] = meshgrid(x_firn,y_firn);

is_firn = interp2(XX_firn,YY_firn,FirnArea,lon_RACMO5(:), lat_RACMO5(:),'nearest');
is_firn(isnan(is_firn)) = 0;

for i = 1:length(FAC10_RACMO5(100,100,:))
    if mod(i,100) == 0
        fprintf('%i  out of %i\n',i, length(FAC10_RACMO5(100,100,:)))
    end
  
    temp = FAC10_RACMO5(:,:,i);
    ind_interp = and(isnan(temp(:)),is_firn);
    F = scatteredInterpolant(lon_RACMO5(~isnan(temp(:))), ...
        lat_RACMO5(~isnan(temp(:))),  temp(~isnan(temp(:))));
    temp(ind_interp) = F(lon_RACMO5(ind_interp),lat_RACMO5(ind_interp));
    FAC10_RACMO5(:,:,i)= reshape(temp,size(lat_RACMO5));

    temp = FAC100_RACMO5(:,:,i);
    ind_interp = and(isnan(temp(:)),is_firn);
    F = scatteredInterpolant(lon_RACMO5(~isnan(temp(:))), ...
        lat_RACMO5(~isnan(temp(:))),  temp(~isnan(temp(:))));
    temp(ind_interp) = F(lon_RACMO5(ind_interp),lat_RACMO5(ind_interp));
    FAC100_RACMO5(:,:,i)= reshape(temp,size(lat_RACMO5));
end
 
DV = datevec(metadata.Date);
DaysInYear = yeardays(DV(:,1));
DV2 = DV;
DV2(:,2:end) = 0;
doy = datenum(DV)-datenum(DV2);
metadata.dec_year = DV(:,1) + doy./DaysInYear;
 
ind_time = dsearchn(time_RACMO5, metadata.dec_year);
[i_lat, i_lon] = ind2sub(size(lat_RACMO5), dsearchn([lat_RACMO5(:) lon_RACMO5(:)], [metadata.Latitude metadata.Longitude]));
metadata.FAC10_RACMO = metadata.FAC10*NaN;
metadata.FAC100_RACMO = metadata.FAC10*NaN;
for k = 1: length(metadata.FAC10)
    metadata.FAC10_RACMO(k) = FAC10_RACMO5(i_lat(k),i_lon(k),ind_time(k));
    metadata.FAC100_RACMO(k) = FAC100_RACMO5(i_lat(k),i_lon(k),ind_time(k));
end

metadata.FAC10_RACMO(metadata.dec_year<1960) = NaN;
metadata.FAC100_RACMO(metadata.dec_year<1960) = NaN;

% metadata.Date = datenum(metadata.Date);
% M = table2array(metadata(:,[2 3 5 6 17:24]));
% csvwrite('RACMO_FAC.csv',M)
toc

%% Loading gridded data HIRHAM
disp('Loading gridded data HIRHAM')
tic
load('.\Input\RCM FAC\HIRHAM5\FAC_HH')

filename = ['.\Input\RCM FAC\HIRHAM5\GL2LIN_Darcy_60m_liqCL_wh1_RHOFIRN_1980_MM.nc'];
    if exist(filename) ~= 2 
        error(['The file %s is not publicly available.\n' ...
            'Contact Baptiste Vandecrux (b.vandecrux@gmail.com) or RCM managers to obtain RCM FACs.'],...
            filename);
    end
lat_HH = double(ncread(filename, 'lat')); 
lon_HH =double( ncread(filename, 'lon'));

FAC10_HH = FAC_10_HH;
FAC100_HH = FAC_100_HH;
clearvars FAC_10_HH FAC_100_HH

time_HH = [];
for i = 1980:2016 
    time_HH = [time_HH datenum(i,1:12,1)];
end
time_HH = [time_HH datenum(2017,1:8,1)];
toc

%% Plotting all RCM
pos_title = [0.03 0.75];
h_tit=[];
f = figure('Visible',vis,'outerposition',[1 1 20 25]);
ha = tight_subplot(2,4,[0.0 0.01],0.17,0.1);
set(f,'CurrentAxes',ha(1))
out_HH10_LIN = Plotting_FAC_RCM_Comp(metadata.FAC10_HH_LIN, metadata.FAC10,...
    metadata.T_avg, metadata.c_avg, 0.3);
title('a) HH_LIN','Interpreter','none');
axis square tight
    ylabel('Observed FAC_{10} (m)','interpreter','tex')

set(f,'CurrentAxes',ha(5))
out_HH100_LIN = Plotting_FAC_RCM_Comp(metadata.FAC100_HH_LIN, metadata.FACpco,...
    metadata.T_avg, metadata.c_avg, 1);
title('e) HH_LIN','Interpreter','none');
axis square tight
    ylabel('Observed FAC_{tot} (m)','interpreter','tex')

set(f,'CurrentAxes',ha(2))
out_HH10_MOD = Plotting_FAC_RCM_Comp(metadata.FAC10_HH_MOD, metadata.FAC10,...
    metadata.T_avg, metadata.c_avg, 0.3);
title('b) HH_MOD','Interpreter','none');
axis square tight

set(f,'CurrentAxes',ha(6))
out_HH100_MOD = Plotting_FAC_RCM_Comp(metadata.FAC100_HH_MOD, metadata.FACpco,...
    metadata.T_avg, metadata.c_avg, 1);
title('f) HH_MOD','Interpreter','none');
axis square tight
h_x1 = xlabel('Modelled FAC_{tot}  (m)','interpreter','tex');


set(f,'CurrentAxes',ha(3))
out_RACMO10 = Plotting_FAC_RCM_Comp(metadata.FAC10_RACMO, metadata.FAC10,...
    metadata.T_avg, metadata.c_avg, 0.3);
title('c) RACMO2.3p2','Interpreter','tex');
axis square tight
    
set(f,'CurrentAxes',ha(7))
out_RACMO100 = Plotting_FAC_RCM_Comp(metadata.FAC100_RACMO, metadata.FACpco,...
    metadata.T_avg, metadata.c_avg, 1);
title('g) RACMO2.3p2','Interpreter','tex');
axis square tight

set(f,'CurrentAxes',ha(4))
out_MAR10 = Plotting_FAC_RCM_Comp(metadata.FAC10_MAR, metadata.FAC10,...
    metadata.T_avg, metadata.c_avg, 0.3);
 title('d) MARv3.9','Interpreter','tex');
axis square tight
h_x2 = xlabel('Modelled FAC_{10}  (m)','interpreter','tex');
h_x2.Units = 'normalized';
h_x2.Position(1) = -1;

h=[];
h(1) = plot(NaN,NaN,'o','MarkerFaceColor',RGB('jaune'),'Color','k','LineWidth',0.5);
h(2) = plot(NaN,NaN,'o','MarkerFaceColor',RGB('rouge'),'Color','k','LineWidth',0.5);
h(3) = plot(NaN,NaN,'o','MarkerFaceColor',RGB('vert'),'Color','k','LineWidth',0.5);
h(4) = plot([NaN NaN], [NaN NaN],'k','LineWidth',2);
h(5) = patch([NaN  NaN NaN NaN ], ...
    [NaN  NaN NaN NaN ],...
    RGB('light light gray'));
% h(4) = plot([NaN NaN], [NaN NaN],'r','LineWidth',2);

legendflex(h,{'DSA','LAPA','HAPA','1:1 line','Uncertainty '},'ref',gcf,'anchor',{'n' 'n'},...
    'ncol',1,'box','on','anchor',{'se' 'se'},'buffer',[-50 180]);

ha(8).Visible = 'off';
ha(4).YAxisLocation = 'right';
ha(7).YAxisLocation = 'right';


for i = 1:4
    if ismember(i,2:3)
    ha(i).YTickLabel = '';
    end
    ha(i).YLim = [0.5 6.5];
    ha(i).XLim = [0.5 6.5];
    ha(i).Title.FontSize = 13;
end
for i = 5:8
    if i==6
    ha(i).YTickLabel = '';
    end
    ha(i).YLim = [15 31];
    ha(i).XLim = [15 31];
    ha(i).Title.FontSize = 13;
end

print(f,'./Output/RCM_comp','-dtiff')
print(f,'./Output/RCM_comp','-dpdf')

M = [out_HH10_LIN;  
    out_HH10_MOD; 
    out_RACMO10;
    out_MAR10;
    out_HH100_LIN;
    out_HH100_MOD;
    out_RACMO100;];
summary_RCMfit = array2table(M);
summary_RCMfit.Properties.VariableNames = ...
    {'bias_DSA' 'RMSD_DSA' 'bias_LAPA' 'RMSD_LAPA' 'bias_HAPA' 'RMSD_HAPA' 'bias_all' 'RMSD_all','Slope','Intercept'};
summary_RCMfit.Properties.RowNames = ...
    {'HIRHAM_10_LIN','HIRHAM_10_MOD','RACMO_10', 'MAR_10', 'HIRHAM_100_LIN','HIRHAM_100_MOD', 'RACMO_100'};

writetable(summary_RCMfit, './Output/summary_RCMfit.csv','Delimiter',';','WriteRowNames',true);

fprintf('%i  %i  %i  %i  \n',sum(metadata.T_avg<T_thresh),...
    sum(and(metadata.T_avg>=T_thresh,metadata.c_avg>=600)),...
    sum(and(metadata.T_avg>=T_thresh,metadata.c_avg<600)),...
    length(metadata.T_avg));
met = metadata;
met(isnan(met.FACpco),:) = [];
fprintf('%i  %i  %i  %i  \n',sum(met.T_avg<T_thresh),...
    sum(and(met.T_avg>=T_thresh,met.c_avg>=600)),...
    sum(and(met.T_avg>=T_thresh,met.c_avg<600)),...
    length(met.T_avg));

%% Calculating spatially integrated FAC in MAR
disp('======= Spatial integration ========')
disp('MAR')
DT = delaunayn([T_map.lat T_map.lon]);
ind_MAR_in_org = dsearchn([T_map.lat T_map.lon], DT, [lat_MAR(:) lon_MAR(:)]); 

ind_DSA_MAR = reshape(ind_DSA_org(ind_MAR_in_org),size(lat_MAR));
ind_LAPA_MAR = reshape(ind_LAPA_org(ind_MAR_in_org),size(lat_MAR));
ind_HAPA_MAR = reshape(ind_HAPA_org(ind_MAR_in_org),size(lat_MAR));
is_firn_org = ~isnan(T_map.T_avg);
ind_GrFirn = reshape(is_firn_org(ind_MAR_in_org),size(lat_MAR));

figure
scatter(lon_MAR(ind_DSA_MAR), lat_MAR(ind_DSA_MAR))
hold on
scatter(lon_MAR(ind_LAPA_MAR), lat_MAR(ind_LAPA_MAR))
scatter(lon_MAR(ind_HAPA_MAR), lat_MAR(ind_HAPA_MAR))
scatter(T_map.lon(~isnan(T_map.T_avg)), T_map.lat(~isnan(T_map.T_avg)),'.')

disp('GrIS')
tic
[FAC10_GrIS_MAR, a_MAR] = CalcFACevo(FAC10_MAR, ind_GrFirn,...
    time_MAR, lat_MAR, lon_MAR,15);

disp('DSA')
[FAC10_DSA_MAR, ~] = CalcFACevo(FAC10_MAR, ind_DSA_MAR,...
    time_MAR, lat_MAR, lon_MAR,15);

disp('LAPA')
[FAC10_LAPA_MAR, ~] = CalcFACevo(FAC10_MAR, ind_LAPA_MAR,...
    time_MAR, lat_MAR, lon_MAR,15);

disp('HAPA')
[FAC10_HAPA_MAR, ~] = CalcFACevo(FAC10_MAR, ind_HAPA_MAR,...
    time_MAR, lat_MAR, lon_MAR,15);
toc

%% Calculating spatially integrated FAC in RACMO5
disp('RACMO5 10m')
DT = delaunayn([T_map.lat T_map.lon]);
ind_RACMO5_in_org = dsearchn([T_map.lat T_map.lon], DT, [lat_RACMO5(:) lon_RACMO5(:)]); 

ind_DSA_RACMO5 = reshape(ind_DSA_org(ind_RACMO5_in_org),size(lat_RACMO5));
ind_LAPA_RACMO5 = reshape(ind_LAPA_org(ind_RACMO5_in_org),size(lat_RACMO5));
ind_HAPA_RACMO5 = reshape(ind_HAPA_org(ind_RACMO5_in_org),size(lat_RACMO5));
is_firn_org = ~isnan(T_map.T_avg);
ind_GrFirn = reshape(is_firn_org(ind_RACMO5_in_org),size(lat_RACMO5));
disp('GrIS')
tic
mask = ind_GrFirn;
[FAC10_GrIS_RACMO5, a_RACMO5] = CalcFACevo(FAC10_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('DSA')
mask = ind_DSA_RACMO5;
[FAC10_DSA_RACMO5, ~] = CalcFACevo(FAC10_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('LAPA')
mask = ind_LAPA_RACMO5;
[FAC10_LAPA_RACMO5, ~] = CalcFACevo(FAC10_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('HAPA')
mask = ind_HAPA_RACMO5;
[FAC10_HAPA_RACMO5, ~] = CalcFACevo(FAC10_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);
toc

%% Calculating spatially integrated FAC100 in RACMO5
disp('RACMO5 100m')
tic
mask = ind_GrFirn;
[FAC100_GrIS_RACMO5, a_RACMO5] = CalcFACevo(FAC100_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('DSA')
mask = ind_DSA_RACMO5;
[FAC100_DSA_RACMO5, ~] = CalcFACevo(FAC100_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('LAPA')
mask = ind_LAPA_RACMO5;
[FAC100_LAPA_RACMO5, ~] = CalcFACevo(FAC100_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);

disp('HAPA')
mask = ind_HAPA_RACMO5;
[FAC100_HAPA_RACMO5, ~] = CalcFACevo(FAC100_RACMO5, mask,...
    time_RACMO5, lat_RACMO5, lon_RACMO5,5.5);
toc

%% Calculating spatially integrated FAC10 in HH
disp('FAC 10 HH')
DT = delaunayn([T_map.lat T_map.lon]);
ind_HH_in_org = dsearchn([T_map.lat T_map.lon], DT, [lat_HH(:) lon_HH(:)]); 

ind_DSA_HH = reshape(ind_DSA_org(ind_HH_in_org),size(lat_HH));
ind_LAPA_HH = reshape(ind_LAPA_org(ind_HH_in_org),size(lat_HH));
ind_HAPA_HH = reshape(ind_HAPA_org(ind_HH_in_org),size(lat_HH));
is_firn_org = ~isnan(T_map.T_avg);
ind_GrFirn = reshape(is_firn_org(ind_HH_in_org),size(lat_HH));
disp('GrIS')
tic
mask = ind_GrFirn;
[FAC10_GrIS_HH, a_HH] = CalcFACevo(FAC10_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('DSA')
mask = ind_DSA_HH;
[FAC10_DSA_HH, ~] = CalcFACevo(FAC10_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('LAPA')
mask = ind_LAPA_HH;
[FAC10_LAPA_HH, ~] = CalcFACevo(FAC10_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('HAPA')
mask = ind_HAPA_HH;
[FAC10_HAPA_HH, ~] = CalcFACevo(FAC10_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);
toc

%% Calculating spatially integrated FAC100 in HH
disp('FAC 100 HH')
% DT = delaunayn([T_map.lat T_map.lon]);
% ind_HH_in_org = dsearchn([T_map.lat T_map.lon], DT, [lat_HH(:) lon_HH(:)]); 

ind_DSA_HH = reshape(ind_DSA_org(ind_HH_in_org),size(lat_HH));
ind_LAPA_HH = reshape(ind_LAPA_org(ind_HH_in_org),size(lat_HH));
ind_HAPA_HH = reshape(ind_HAPA_org(ind_HH_in_org),size(lat_HH));
is_firn_org = ~isnan(T_map.T_avg);
ind_GrFirn = reshape(is_firn_org(ind_HH_in_org),size(lat_HH));
disp('GrIS')
tic
mask = ind_GrFirn;
[FAC100_GrIS_HH, a_HH] = CalcFACevo(FAC100_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('DSA')
mask = ind_DSA_HH;
[FAC100_DSA_HH, ~] = CalcFACevo(FAC100_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('LAPA')
mask = ind_LAPA_HH;
[FAC100_LAPA_HH, ~] = CalcFACevo(FAC100_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);

disp('HAPA')
mask = ind_HAPA_HH;
[FAC100_HAPA_HH, ~] = CalcFACevo(FAC100_HH, mask,...
    time_HH, lat_HH, lon_HH,5.5);
toc

%% FAC Evolution plot
filename = '.\Output\Summary_output_MAR.csv';
delimiter = ';';
startRow = 2;
formatSpec = '%s%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
summary = table(dataArray{1:end-1}, 'VariableNames', {'Row','DSA','HAPA','LAPA_pre_2010','LAPA_post_2010','Gr'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
summary.Properties.RowNames = summary.Row;
summary(:,1) = [];
set(groot,'defaultAxesColorOrder',brewermap(6,'set1'))

f = figure('outerposition',[1 1 20 15]);
ha = tight_subplot(4,2,[0.03 0.01],[0.07 0.1],0.1);
set(gcf,'CurrentAxes',ha(1))
    hold on
    h1 = plot(datenum(time_RACMO5,1,1),0.001*FAC10_GrIS_RACMO5,'LineWidth',1.5);
    h2 = plot(time_HH, 0.001*FAC10_GrIS_HH, 'LineWidth',1.5);
    h = plot(time_MAR,0.001*FAC10_GrIS_MAR,'LineWidth',1.5);
    h3 = plot(time_HH, NaN*FAC10_GrIS_HH, 'w','LineWidth',1.5);
    h4 = plot(datenum([2010 2017],1,1), ...
        0.001*[1 1]*summary.Gr(3),'k','LineWidth',1.5);
    h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
        0.001*(summary.Gr(3)+ ...
        [summary.Gr([4 4]); - summary.Gr([4 4])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';

    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    legendflex([h h2 h1 h3], {'MARv3.9', 'HIRHAM5_MOD', 'RACMO2.3p2'},...
        'nrow',1,'box', 'off','ref',gcf,'anchor',{'n' 'n'},'buffer',[20 0],...
        'Interpreter','none')
    legendflex([h4 h_unc], {'Derived from observations','Uncertainty'},...
        'nrow',1,'box', 'off','ref',gcf,'anchor',{'n' 'n'},'buffer',[30 -20])
    uistack(h,'bottom')
    uistack(h_unc,'bottom')

    h_tit = title('a) All firn area','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position(1:2) = [0.25 0.1  ];
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)
    ylim([5.5 8])

set(gcf,'CurrentAxes',ha(3))
    hold on
    h=plot(datenum(time_RACMO5,1,1),0.001*FAC10_DSA_RACMO5,'LineWidth',1.5);
    plot(time_HH, 0.001*FAC10_DSA_HH, 'LineWidth',1.5);
    plot(time_MAR,0.001*FAC10_DSA_MAR,'LineWidth',1.5)
        uistack(h,'top')
    plot(datenum([1953 2017],1,1), 0.001*[1 1]*summary.DSA(3),'k','LineWidth',1.5)
   h_unc = patch(datenum([1953 2017 2017 1953],1,1), ...
        0.001*(summary.DSA(3)+ ...
        [summary.DSA([4 4]); - summary.DSA([4 4])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
    set_monthly_tick(datenum(time_RACMO5,1,1));
    ylim([4.5 6.5])
    xlim(datenum([1958 2016],1,1))
    h_tit = title('b) DSA  ','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized';
    h_tit.Position(1:2) = [0.17 0.68  ];
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)

set(gcf,'CurrentAxes',ha(5))
    hold on
    h= plot(datenum(time_RACMO5,1,1),0.001*FAC10_LAPA_RACMO5,'LineWidth',1.5);
    plot(time_HH, 0.001*FAC10_LAPA_HH, 'LineWidth',1.5);
    plot(time_MAR,0.001*FAC10_LAPA_MAR,'LineWidth',1.5)
        uistack(h,'top')

    plot(datenum([1997 2010],1,1), ...
        0.001*[1 1]*summary.LAPA_pre_2010(3),'k','LineWidth',1.5)
  h_unc = patch(datenum([1997 2010 2010 1997],1,1), ...
        0.001*(summary.LAPA_pre_2010(3)+ ...
        [summary.LAPA_pre_2010([4 4]); - summary.LAPA_pre_2010([4 4])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')

    plot(datenum([2010 2017],1,1), 0.001*[1 1]*summary.LAPA_post_2010(3),'k','LineWidth',1.5)
   h_unc = patch(datenum([ 2010 2017 2017 2010 ],1,1), ...
        0.001*(summary.LAPA_post_2010(3)+ ...
        [summary.LAPA_post_2010([4 4]); - summary.LAPA_post_2010([4 4])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
    
    ylim([0.3 0.9])
    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)
    h_tit = title('c) LAPA','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position(1:2) = [0.17 0.17  ];
    h_lab=ylabel('Spatially integrated FAC_{10} (10^{3} km^3)','Interpreter','tex');
    h_lab.Units='normalized';
    h_lab.Position(1:2) = [-0.11 1];

set(gcf,'CurrentAxes',ha(7))
    hold on
    h=plot(datenum(time_RACMO5,1,1),0.001*FAC10_HAPA_RACMO5,'LineWidth',1.5);
    h2 =plot(time_HH, 0.001*FAC10_HAPA_HH, 'LineWidth',1.5);
    plot(time_MAR,0.001*FAC10_HAPA_MAR,'LineWidth',1.5)
        uistack(h,'top')
        uistack(h2,'top')
    plot(datenum([2010 2017],1,1), 0.001*[1 1]*summary.HAPA(3),'k','LineWidth',1.5)
   h_unc = patch(datenum([ 2010 2017 2017 2010 ],1,1), ...
        0.001*(summary.HAPA(3)+ ...
        [summary.HAPA([4 4]); - summary.HAPA([4 4])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
        ylim([0.3 0.9])
    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    h_tit = title('d)HAPA','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position(1:2) = [0.17 0.17 ];
    set(gca, 'layer','top','TickLength', get(gca,'TickLength')*1.5)

set(gcf,'CurrentAxes',ha(2))
    hold on
    plot(datenum(time_RACMO5,1,1),0.001*FAC100_GrIS_RACMO5,'LineWidth',1.5)
    plot(time_HH, 0.001*FAC100_GrIS_HH, 'LineWidth',1.5);

    plot(datenum([2010 2017],1,1), ...
        0.001*[1 1]*summary.Gr(7),'k','LineWidth',1.5)
       h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
        0.001*(summary.Gr(7)+ ...
        [summary.Gr([8 8]); - summary.Gr([8 8])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    h_tit = title('e) All firn area','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position(1:2) = [0.25 0.1];
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,...
        'YAxisLocation','right')
    ylim([23 35])

    set(gcf,'CurrentAxes',ha(4))
    hold on
    plot(datenum(time_RACMO5,1,1),0.001*FAC100_DSA_RACMO5,'LineWidth',1.5)
     plot(time_HH, 0.001*FAC100_DSA_HH, 'LineWidth',1.5);
   plot(datenum([1953 2017],1,1), ...
       0.001*[1 1]*summary.DSA(7),'k','LineWidth',1.5)
      h_unc = patch(datenum([1953 2017 2017 1953],1,1), ...
        0.001*(summary.DSA(7)+ ...
        [summary.DSA([8 8]); - summary.DSA([8 8])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
    set_monthly_tick(datenum(time_RACMO5,1,1));
    ylim([20 28])
    xlim(datenum([1958 2016],1,1))
    h_tit = title('f) DSA  ','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position(1:2) = [0.17 0.68  ];
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,...
        'YAxisLocation','right')

set(gcf,'CurrentAxes',ha(6))
    hold on
    plot(datenum(time_RACMO5,1,1),0.001*FAC100_LAPA_RACMO5,'LineWidth',1.5)
    plot(time_HH, 0.001*FAC100_LAPA_HH, 'LineWidth',1.5);

    plot(datenum([1997 2010],1,1), ...
        0.001*[1 1]*summary.LAPA_pre_2010(7),'k','LineWidth',1.5)
       h_unc = patch(datenum([1997 2010 2010 1997],1,1), ...
        0.001*(summary.LAPA_pre_2010(7)+ ...
        [summary.LAPA_pre_2010([8 8]); - summary.LAPA_pre_2010([8 8])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
    
    plot(datenum([2010 2017],1,1), ...
        0.001*[1 1]*summary.LAPA_post_2010(7),'k','LineWidth',1.5)
       h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
        0.001*(summary.LAPA_post_2010(7)+ ...
        [summary.LAPA_post_2010([8 8]); - summary.LAPA_post_2010([8 8])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')   
    ylim([2 3.5])
    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,'YAxisLocation','right')
    h_tit = title('g) LAPA','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    h_tit.Position = [0.18 0.28 0.1];
    h_lab=ylabel('Spatially integrated FAC_{tot} (10^{3} km^3)','Interpreter','tex');
    h_lab.Units='normalized';
    h_lab.Position(1:2) = [1.1 1];

set(gcf,'CurrentAxes',ha(8))
    hold on
    plot(datenum(time_RACMO5,1,1),0.001*FAC100_HAPA_RACMO5,'LineWidth',1.5)
    plot(time_HH, 0.001*FAC100_HAPA_HH, 'LineWidth',1.5);
    plot(datenum([2010 2017],1,1), ...
        0.001*[1 1]*summary.HAPA(7),'k','LineWidth',1.5)
    h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
        0.001*(summary.HAPA(7)+ ...
        [summary.HAPA([8 8]); - summary.HAPA([8 8])]),...
        RGB('light light gray'));
    h_unc.LineStyle = 'none';
    uistack(h_unc,'bottom')
     ylim([1.5 4.3])
    set_monthly_tick(datenum(time_RACMO5,1,1));
    xlim(datenum([1958 2016],1,1))
    h_tit = title('h) HAPA','Interpreter','tex');
    h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
    set(gca, 'layer','top','TickLength', get(gca,'TickLength')*1.5,'YAxisLocation','right')
    h_tit.Position = [0.21 0.17 0.1];

print(f,'./Output/FAC_evolution','-dtiff')
print(f,'./Output/FAC_evolution','-dpdf')

%% Calculating FC
% FC10_DSA_HH         =	max(0, 10*843 - 917*(10-	FAC10_DSA_HH	));
% FC10_DSA_RACMO5     =	max(0, 10*843 - 917*(10-	FAC10_DSA_RACMO5	));
% FC10_DSA_MAR        =	max(0, 10*843 - 917*(10-	FAC10_DSA_MAR	));
% FC10_GrIS_HH        =	max(0, 10*843 - 917*(10-	FAC10_GrIS_HH	));
% FC10_GrIS_RACMO5	=	max(0, 10*843 - 917*(10-	FAC10_GrIS_RACMO5	));
% FC10_GrIS_MAR       =	max(0, 10*843 - 917*(10-	FAC10_GrIS_MAR	));
% FC10_HAPA_HH       =	max(0, 10*843 - 917*(10-	FAC10_HAPA_HH	));
% FC10_HAPA_RACMO5	=	max(0, 10*843 - 917*(10-	FAC10_HAPA_RACMO5	));
% FC10_HAPA_MAR      =	max(0, 10*843 - 917*(10-	FAC10_HAPA_MAR	));
% FC10_LAPA_HH       =	max(0, 10*843 - 917*(10-	FAC10_LAPA_HH	));
% FC10_LAPA_RACMO5	=	max(0, 10*843 - 917*(10-	FAC10_LAPA_RACMO5	));
% FC10_LAPA_MAR      =	max(0, 10*843 - 917*(10-	FAC10_LAPA_MAR	));
% 				
% FC100_DSA_HH        =	max(0, 100*843 - 917*(100-	FAC100_DSA_HH	));
% FC100_DSA_RACMO5	=	max(0, 100*843 - 917*(100-	FAC100_DSA_RACMO5	));
% FC100_GrIS_HH       =	max(0, 100*843 - 917*(100-	FAC100_GrIS_HH	));
% FC100_GrIS_RACMO5	=	max(0, 100*843 - 917*(100-	FAC100_GrIS_RACMO5	));
% FC100_HAPA_HH      =	max(0, 100*843 - 917*(100-	FAC100_HAPA_HH	));
% FC100_HAPA_RACMO5	=	max(0, 100*843 - 917*(100-	FAC100_HAPA_RACMO5	));
% FC100_LAPA_HH      =	max(0, 100*843 - 917*(100-	FAC100_LAPA_HH	));
% FC100_LAPA_RACMO5	=	max(0, 100*843 - 917*(100-	FAC100_LAPA_RACMO5	));
% 
% %% FC evolution
% 
% f = figure('outerposition',[1 1 20 15]);
% ha = tight_subplot(4,2,[0.03 0.01],[0.07 0.1],0.1);
% set(gcf,'CurrentAxes',ha(1))
%     hold on
%     h1 = plot(datenum(time_RACMO5,1,1),10^(-6)*FC10_GrIS_RACMO5,'LineWidth',1.5);
%     h2 = plot(time_HH, 10^(-6)*FC10_GrIS_HH, 'LineWidth',1.5);
%     h = plot(time_MAR,10^(-6)*FC10_GrIS_MAR,'LineWidth',1.5);
%     h3 = plot(time_HH, NaN*FC10_GrIS_HH, 'w','LineWidth',1.5);
%     h4 = plot(datenum([2010 2017],1,1), ...
%         0.001*[1 1]*summary.Gr(11),'k','LineWidth',1.5);
%     h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
%         0.001*(summary.Gr(11)+ ...
%         [summary.Gr([12 12]); - summary.Gr([12 12])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
% 
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     legendflex([h h2 h1 h3], {'MARv3.9', 'HIRHAM5_MOD', 'RACMO2.3p2'},...
%         'nrow',1,'box', 'off','ref',gcf,'anchor',{'n' 'n'},'buffer',[20 0],...
%         'Interpreter','none')
%     legendflex([h4 h_unc], {'Derived from observations','Uncertainty'},...
%         'nrow',1,'box', 'off','ref',gcf,'anchor',{'n' 'n'},'buffer',[30 -20])
%     uistack(h,'bottom')
%     uistack(h_unc,'bottom')
% 
%     h_tit = title('a) All firn area ','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position(1:2) = [0.17 0.68  ];
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)
% %     ylim([5.5 8])
% 
% set(gcf,'CurrentAxes',ha(3))
%     hold on
%     h=plot(datenum(time_RACMO5,1,1),10^(-6)*FC10_DSA_RACMO5,'LineWidth',1.5);
%     plot(time_HH, 10^(-6)*FC10_DSA_HH, 'LineWidth',1.5);
%     plot(time_MAR,10^(-6)*FC10_DSA_MAR,'LineWidth',1.5)
%         uistack(h,'top')
%     plot(datenum([1953 2017],1,1), 0.001*[1 1]*summary.DSA(11),'k','LineWidth',1.5)
%    h_unc = patch(datenum([1953 2017 2017 1953],1,1), ...
%         0.001*(summary.DSA(11)+ ...
%         [summary.DSA([12 12]); - summary.DSA([12 12])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%     set_monthly_tick(datenum(time_RACMO5,1,1));
% %     ylim([4.5 6.5])
%     xlim(datenum([1958 2016],1,1))
%     h_tit = title('b) DSA  ','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized';
%     h_tit.Position(1:2) = [0.17 0.68  ];
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)
% 
% set(gcf,'CurrentAxes',ha(5))
%     hold on
%     h= plot(datenum(time_RACMO5,1,1),10^(-6)*FC10_LAPA_RACMO5,'LineWidth',1.5);
%     plot(time_HH, 10^(-6)*FC10_LAPA_HH, 'LineWidth',1.5);
%     plot(time_MAR,10^(-6)*FC10_LAPA_MAR,'LineWidth',1.5)
%         uistack(h,'top')
% 
%     plot(datenum([1997 2010],1,1), ...
%         0.001*[1 1]*summary.LAPA_pre_2010(11),'k','LineWidth',1.5)
%   h_unc = patch(datenum([1997 2010 2010 1997],1,1), ...
%         0.001*(summary.LAPA_pre_2010(11)+ ...
%         [summary.LAPA_pre_2010([12 12]); - summary.LAPA_pre_2010([12 12])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
% 
%     plot(datenum([2010 2017],1,1), 0.001*[1 1]*summary.LAPA_post_2010(11),'k','LineWidth',1.5)
%    h_unc = patch(datenum([ 2010 2017 2017 2010 ],1,1), ...
%         0.001*(summary.LAPA_post_2010(11)+ ...
%         [summary.LAPA_post_2010([12 12]); - summary.LAPA_post_2010([12 12])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%     
% %     ylim([0.3 0.9])
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5)
%     h_tit = title('c) LAPA','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position(1:2) = [0.17 0.17  ];
%     h_lab=ylabel('Retention capacity in the top 10 m of firn (10^{3} Gt)','Interpreter','tex');
%     h_lab.Units='normalized';
%     h_lab.Position(1:2) = [-0.11 1];
% 
% set(gcf,'CurrentAxes',ha(7))
%     hold on
%     h=plot(datenum(time_RACMO5,1,1),10^(-6)*FC10_HAPA_RACMO5,'LineWidth',1.5);
%     h2 =plot(time_HH, 10^(-6)*FC10_HAPA_HH, 'LineWidth',1.5);
%     plot(time_MAR,10^(-6)*FC10_HAPA_MAR,'LineWidth',1.5)
%         uistack(h,'top')
%         uistack(h2,'top')
%     plot(datenum([2010 2017],1,1), 0.001*[1 1]*summary.HAPA(11),'k','LineWidth',1.5)
%    h_unc = patch(datenum([ 2010 2017 2017 2010 ],1,1), ...
%         0.001*(summary.HAPA(11)+ ...
%         [summary.HAPA([12 12]); - summary.HAPA([12 12])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%         ylim([0.3 0.9])
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     h_tit = title('d)HAPA','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position(1:2) = [0.17 0.17 ];
%     set(gca, 'layer','top','TickLength', get(gca,'TickLength')*1.5)
% 
% set(gcf,'CurrentAxes',ha(2))
%     hold on
%     plot(datenum(time_RACMO5,1,1),10^(-6)*FC100_GrIS_RACMO5,'LineWidth',1.5)
%     plot(time_HH, 10^(-6)*FC100_GrIS_HH, 'LineWidth',1.5);
% 
%     plot(datenum([2010 2017],1,1), ...
%         0.001*[1 1]*summary.Gr(17),'k','LineWidth',1.5)
%        h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
%         0.001*(summary.Gr(17)+ ...
%         [summary.Gr([18 18]); - summary.Gr([18 18])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     h_tit = title('e) All firn area','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position(1:2) = [0.17 0.68  ];
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,...
%         'YAxisLocation','right')
% %     ylim([21 38])
% 
%     set(gcf,'CurrentAxes',ha(4))
%     hold on
%     plot(datenum(time_RACMO5,1,1),10^(-6)*FC100_DSA_RACMO5,'LineWidth',1.5)
%      plot(time_HH, 10^(-6)*FC100_DSA_HH, 'LineWidth',1.5);
%    plot(datenum([1953 2017],1,1), ...
%        0.001*[1 1]*summary.DSA(17),'k','LineWidth',1.5)
%       h_unc = patch(datenum([1953 2017 2017 1953],1,1), ...
%         0.001*(summary.DSA(17)+ ...
%         [summary.DSA([18 18]); - summary.DSA([18 18])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%     set_monthly_tick(datenum(time_RACMO5,1,1));
% %     ylim([15 30])
%     xlim(datenum([1958 2016],1,1))
%     h_tit = title('f) DSA  ','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position(1:2) = [0.17 0.68  ];
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,...
%         'YAxisLocation','right')
% 
% set(gcf,'CurrentAxes',ha(6))
%     hold on
%     plot(datenum(time_RACMO5,1,1),10^(-6)*FC100_LAPA_RACMO5,'LineWidth',1.5)
%     plot(time_HH, 10^(-6)*FC100_LAPA_HH, 'LineWidth',1.5);
% 
%     plot(datenum([1997 2010],1,1), ...
%         0.001*[1 1]*summary.LAPA_pre_2010(17),'k','LineWidth',1.5)
%        h_unc = patch(datenum([1997 2010 2010 1997],1,1), ...
%         0.001*(summary.LAPA_pre_2010(17)+ ...
%         [summary.LAPA_pre_2010([18 18]); - summary.LAPA_pre_2010([18 18])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
%     
%     plot(datenum([2010 2017],1,1), ...
%         0.001*[1 1]*summary.LAPA_post_2010(17),'k','LineWidth',1.5)
%        h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
%         0.001*(summary.LAPA_post_2010(17)+ ...
%         [summary.LAPA_post_2010([18 18]); - summary.LAPA_post_2010([18 18])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')   
% %     ylim([1.5 4])
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     set(gca,'XTickLabel','', 'layer','top','TickLength', get(gca,'TickLength')*1.5,'YAxisLocation','right')
%     h_tit = title('g) LAPA','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     h_tit.Position = [0.21 0.68 0.1];
%     h_lab=ylabel('Retention capacity in whole firn column (10^{3} Gt)','Interpreter','tex');
%     h_lab.Units='normalized';
%     h_lab.Position(1:2) = [1.1 1];
% 
% set(gcf,'CurrentAxes',ha(8))
%     hold on
%     plot(datenum(time_RACMO5,1,1),10^(-6)*FC100_HAPA_RACMO5,'LineWidth',1.5)
%     plot(time_HH, 10^(-6)*FC100_HAPA_HH, 'LineWidth',1.5);
%     plot(datenum([2010 2017],1,1), ...
%         0.001*[1 1]*summary.HAPA(17),'k','LineWidth',1.5)
%     h_unc = patch(datenum([2010 2017 2017 2010],1,1), ...
%         0.001*(summary.HAPA(17)+ ...
%         [summary.HAPA([18 18]); - summary.HAPA([18 18])]),...
%         RGB('light light gray'));
%     h_unc.LineStyle = 'none';
%     uistack(h_unc,'bottom')
% %      ylim([1.5 4.3])
%     set_monthly_tick(datenum(time_RACMO5,1,1));
%     xlim(datenum([1958 2016],1,1))
%     h_tit = title('h) HAPA','Interpreter','tex');
%     h_tit.FontSize = 13; h_tit.Units = 'normalized'; 
%     set(gca, 'layer','top','TickLength', get(gca,'TickLength')*1.5,'YAxisLocation','right')
%     h_tit.Position = [0.21 0.17 0.1];
% 
% print(f,'./Output/FC_evolution','-dtiff')
% print(f,'./Output/FC_evolution','-dpdf')


