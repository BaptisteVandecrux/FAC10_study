function [metadata] = Create_FAC10_dataset(T_map,b_map,vis)
    addpath('../matlab_functions/toolbox_cores')
    load ../matlab_functions/Core_all
    %defining the pore space function
ps = @(m,rho) max(0,m.*(1./rho -1/873));
            %% Assigning long-term average temperature and accumulation
% Core_org = Core;
disp('Assigning long-term average temperature and accumulation')
Core{FindCore(Core,'Name','Dome GRIP')}.Data.Density(1:453)=NaN;

tic
for i = 1:length(Core)        
    % accumulation
     [~, ind] = min(distance(Core{i}.Info.Coord1(1),-abs(Core{i}.Info.Coord1(2)),b_map.lat,b_map.lon));
    Core{i}.Info.b_lt_mean = b_map.b_avg(ind);
    %temperature
     [~, ind] = min(distance(Core{i}.Info.Coord1(1),-abs(Core{i}.Info.Coord1(2)),T_map.lat,T_map.lon));
    Core{i}.Info.T_lt_mean =  T_map.T_avg(ind);
end
toc


%% data treatment for each core
disp('Data treatment for each core')
% figure
if ~isempty(FindCore(Core,'Name','Dome GRIP'))
    Core{FindCore(Core,'Name','Dome GRIP')}.Data.Density(1:453) = ...
        Core{FindCore(Core,'Name','Mayewski_1990')}.Data.Density(1:453);
    Core{FindCore(Core,'Name','Dome GRIP')}.Info.Name = 'GRIP & Mayewski';
end

Core{FindCore(Core,'Name','T41e_Spring_2004')}.Data.Density = flipud(Core{FindCore(Core,'Name','T41e_Spring_2004')}.Data.Density);
Core{FindCore(Core,'Name','T41e_Spring_2004')}.Data.Depth = flipud(Core{FindCore(Core,'Name','T41e_Spring_2004')}.Data.Depth);

count = 1;

for i = 1:length(Core)
    % putting everything in column
    text_cit= Core{i}.Info.Citation;
    ind = strfind(text,';');
    text_cit(ind) = ','; 
    Core{i}.Info.Citation = text_cit;
    if ~iscolumn(Core{i}.Data.Density)
        Core{i}.Data.Density=Core{i}.Data.Density';
    end
    if ~iscolumn(Core{i}.Data.Depth)
        Core{i}.Data.Depth=Core{i}.Data.Depth';
    end
    if strcmp('', Core{i}.Info.Densities)
        Core{i}.Info.Densities = 'y';
    end
    if and(strcmp(Core{i}.Info.Densities,'y'),...
            length(Core{i}.Data.Density)>500)
        % adding nan values to density so that density depth and type have
        % the same length

        depth = Core{i}.Data.Depth;
        density = Core{i}.Data.Density;
        density(density==0) = NaN;
%         ice_frac = Core{i}.Data.Type_perc;
            ice_frac = zeros(size(density)); % ### here we decide to ignore ice fraction data

        % Cropping depth and density arrays at the last non-nan density
        ind_last = find(~isnan(density),1,'last');
        if length(density) - ind_last < 100
            depth = depth(1:ind_last);
            density = density(1:ind_last);
            ice_frac = ice_frac(1:min(length(ice_frac),ind_last));

            Core{i}.Data.Depth = depth;
            Core{i}.Data.Density = density;
            Core{i}.Data.Type_perc = ice_frac;
        end

        % creating artificial ice_frac column when absent
        if isempty(ice_frac) || ...
                strcmp(Core{i}.Info.Name,'core_2015T_A5')|| ...
                strcmp(Core{i}.Info.Name,'WindSled_2016')
            ice_frac = zeros(size(density));
        end
        
       
        % here we're going to partition each section between the ice
        % fraction and the firn fraction, each one of them having their own
        % density
        
        % spotting the cells where ice is present
        is_ice = ice_frac > 1;
        is_ice (min(find(is_ice) +1,length(is_ice))) = 1;
        is_ice (max(find(is_ice)-1,1)) = 1;
        
        % removing their density from the firn-only density
        density_firn = density;
        density_firn(find(is_ice)) = NaN;
        
        % if rho>750 at shallow depth then it must contain ice
%         is_ice = and(density_firn>750 , depth<15);
%         density_firn(is_ice) = NaN;

        % assigning surface density if not present
        rho_surf = density(1);
        if isnan(rho_surf)
            rho_surf = 315;
        end
       
        % first attempt to fit a exponential density-depth function
        density_func = @(z_rho, x) ...
                        873 -(873 - rho_surf)*exp(-(x-1)/z_rho);
        ind_nan = isnan(density_firn);
        if sum(ind_nan) >0
            fitted_density = fit(depth(~ind_nan), ...
                density_firn(~ind_nan), density_func);

            % if not successful (mainly because of reduced depth range of a
            % core) fitting a linear function
            if fitted_density.z_rho<5
                density_func = @(a,  x) ...
                            min(873, a * (x-1) + 315);
    %         density_func = @(z_rho,rho_ice ,x) ...
    %                         rho_ice -(rho_ice - rho_surf)*exp(-x/z_rho);
                fitted_density = fit(depth(~ind_nan), ...
                    density_firn(~ind_nan), density_func);
            end

            % exception for Dome GRIP
            if strcmp(Core{i}.Info.Name,'Dome GRIP')
                fitted_density = fit(depth(~ind_nan), ...
                    density_firn(~ind_nan), 'poly3');
            end
% plot(depth, density_firn)
% hold on
% plot(depth, fitted_density(depth))

            % assigning interpolated densities
            density_firn(ind_nan) = fitted_density(depth(ind_nan));
% plot(depth, density_firn)
% hold off
% pause
            ind_nan = find(isnan(density));
            density_org = zeros(size(depth))==1;

            for ii = 1:length(ind_nan)
                ice = ice_frac(ind_nan(ii))/100;
                density(ind_nan(ii)) = density_firn(ind_nan(ii)) ...
                    + (873 - density_firn(ind_nan(ii))) * ice;

                if isnan(density(ind_nan(ii)))
                    density(ind_nan(ii)) = density_firn(ind_nan(ii));
                end

                density_org(ind_nan(ii)) = 1;
            end
            Core{i}.Data.Density = density;
            Core{i}.Data.Density_origin = density_org;
        end
        % Plotting firn air content versus depth and extrapolating to 10 m
            if max(Core{i}.Data.Depth)>1000
                [Core{i}.Data.DensityHL,~]=densitymodel(Core{i}.Info.T_lt_mean, ...
                    Core{i}.Info.b_lt_mean/1000,315,Core{i}.Data.Depth/100, 'NabarroHerring');
            else
                [Core{i}.Data.DensityHL,~]=densitymodel(Core{i}.Info.T_lt_mean, ...
                    Core{i}.Info.b_lt_mean/1000,315,0.01:0.01:10, 'NabarroHerring');
            end
                

            Core{i}.Data.PoreVolume = cumsum(ps(Core{i}.Data.Density.*0.01,Core{i}.Data.Density));

    end
end

%% Gap-filling pore space
disp('Extrapolating pore space')
f= figure('Visible',vis);
ha = tight_subplot(3,5,0.03, [.1 0.02], 0.08);
plotting = 1;
count = 1;
Core_save = Core;

for i = 1:length(Core) 
    if and(strcmp(Core{i}.Info.Densities,'y'),...
            length(Core{i}.Data.Density)>500)
        
            PV = (Core{i}.Data.PoreVolume);
            depth = Core{i}.Data.Depth;
             
             if max(Core{i}.Data.Depth)<1000
%                 PV_fit_func = fit(depth, PV, 'poly2');
                 Core{i}.Data.Depth = [1:1000]';
                 depth = Core{i}.Data.Depth;

                 if strcmp(Core{i}.Info.Name,'H4-2') || strcmp(Core{i}.Info.Name,'H5-1')
                     PV_fit_func = fit(depth(200:length(PV)), PV(200:length(PV)), 'poly1');
                     PV_2 = PV_fit_func(Core{i}.Data.Depth);
                     Core{i}.Data.PoreVolume = [PV; ...
                         PV_2([(length(PV)+1):1000]')];
                 else
                     % now we go through the available cores to look for
                     % similare PV profiles
                     RMSE = 1000000*ones(length(Core),1);
                     for j = 1:length(Core_save)
                         if max(Core_save{j}.Data.Depth)<1000 || ...
                                 ~ and(strcmp(Core_save{j}.Info.Densities,'y'),...
                                                length(Core_save{j}.Data.Density)>500) || ...
                                                strcmp(Core_save{j}.Info.Name,'WindSled_2016') || ...
                                                strcmp(Core_save{j}.Info.Name,'core_2015T_A5')                                            
                             continue
                         else
                             PV_2 = Core_save{j}.Data.PoreVolume;
                             ind = 1:length(PV);
                             RMSE(j) = sqrt(mean((PV(ind)-PV_2(ind)).^2));
                         end
                     end
                     [RMSE_best, ind_filling_core] = min(RMSE);
                     fprintf('%s - %s - RMSE: %0.2f\n',...
                         Core{i}.Info.Name,...
                         Core{ind_filling_core}.Info.Name,...
                         RMSE_best);
                    PV_2 = Core_save{ind_filling_core}.Data.PoreVolume;
                    % PV_2 is cumulative so we need to extract the vertically
                    % resolved PV and take the part that we need for gap
                    % filling

                    PV_2_z = PV_2;
                    PV_2_z (2:end)= PV_2(2:end) - PV_2(1:end-1);

                    PV_z = PV;
                    PV_z (2:end)= PV(2:end) - PV(1:end-1);

                    PV_z_filled = [PV_z; PV_2_z([(length(PV_z)+1):1000]')];
                     Core{i}.Data.PoreVolume = cumsum(PV_z_filled);
                 end
                 
%                  Core{i}.Data.Density = [Core{i}.Data.Density; ...
%                      ps2rho(PV_2([(length(PV)+1):1000]'))];
             
                 Core{i}.Data.Density_origin((length(PV)+1):1000) = 1;


             
                if length(Core{i}.Data.PoreVolume)<1000
                    error('PV trop petit')
                end
if plotting
                set(f,'CurrentAxes',ha(count))
                hold on
                plot(PV_2(1:1000),-depth(1:1000)/100,'r','LineWidth',1.5)
                plot(PV,-depth(1:length(PV))/100,'b','LineWidth',3)
                text_h = annotation(f,'textbox',...
                    [0.01 0.68 0.1 0.03],...
                    'String',{Core{i}.Info.Name},...
                    'HorizontalAlignment','left',...
                    'FitBoxToText','on',...
                    'Interpreter','none',...
                    'LineStyle','none',...
                    'Color','b');
                set(text_h,'Parent',gca,'Units','Normalized')
                text_h.Position = [1.5    -1    0.5    0.4];
                 if strcmp(Core{i}.Info.Name,'H4-2') || strcmp(Core{i}.Info.Name,'H5-1')
                    text_h = annotation(f,'textbox',...
                        [0.01 0.68 0.1 0.03],...
                        'String','Linear fit',...
                        'HorizontalAlignment','left',...
                        'FitBoxToText','on',...
                        'Interpreter','none',...
                        'LineStyle','none',...
                        'Color','r');
                    set(text_h,'Parent',gca,'Units','Normalized')
                    text_h.Position = [1.7    -2    0.5    0.4];
                 else
                    text_h = annotation(f,'textbox',...
                        [0.01 0.68 0.1 0.03],...
                        'String',{Core{ind_filling_core}.Info.Name},...
                        'HorizontalAlignment','left',...
                        'FitBoxToText','on',...
                        'Interpreter','none',...
                        'LineStyle','none',...
                        'Color','r');
                    set(text_h,'Parent',gca,'Units','Normalized')
                    text_h.Position = [1.7    -2    0.5    0.4];
                 end
                
                if ismember(count,[1 6])
                    set(gca,'XTickLabel','')
                elseif count == 11
                elseif count>11
                    set(gca,'YTickLabel','')
                else
                    set(gca,'XTickLabel','','YTickLabel','')
                end
                if count == 6
                    ylabel('Depth (m)')
                elseif count == 13
                    xlabel('Cumulated firn air content (m^3 m^{-2})','Interpreter','tex')
                end
                xlim([0 6]);
                ylim([-10 0])
                box on
                if count == 15
                    i_file = 1;
                        NameFile = sprintf('Output/Extrapolation pore space/comp_dens_fit_%i.png',i_file)  ;
                    while exist(NameFile, 'file') == 2
                        i_file = i_file + 1;
                        NameFile = sprintf('Output/Extrapolation pore space/comp_dens_fit_%i.png',i_file)  ;
                    end
                    print(f,NameFile,'-dpng');

                    if strcmp(vis,'off')
                        close(f)
                    end
                    f= figure('Visible',vis);
                    ha = tight_subplot(3,5,0.03, [.1 0.02], 0.08);
                    count =1;
                else
                    count = count+1;
                end
end
                if length(Core{i}.Data.PoreVolume)~=length(Core{i}.Data.Depth)
                    error('Depth and density have different length')
                end
             end
    end
end
if plotting

        if count-1<6
        set(f,'CurrentAxes',ha(3))
        xlabel('Cumulated firn air content (m^3/m^2)','Interpreter','tex')
        elseif count-1<13
        set(f,'CurrentAxes',ha(6))
        xlabel('Cumulated firn air content (m^3/m^2)','Interpreter','tex')
        end


        for ii = count:length(ha)
            set(ha(ii),'Visible','off')
        end

                    i_file = 1;
                        NameFile = sprintf('Output/Extrapolation pore space/comp_dens_fit_%i.png',i_file)  ;
                    while exist(NameFile, 'file') == 2
                        i_file = i_file + 1;
                        NameFile = sprintf('Output/Extrapolation pore space/comp_dens_fit_%i.png',i_file)  ;
                    end
                    print(f,NameFile,'-dpng');    
        close all
end

%% Create metadata
disp('Create metadata')
% figure
metadata = table;
metadata.Year = [1:10]';
for i = 1:length(Core)
    metadata.CoreNumber(i) = i;
    metadata.Name{i} = Core{i}.Info.Name;
    metadata.Year(i) = Core{i}.Info.DateCored.Year;
    metadata.Latitude(i) = Core{i}.Info.Coord1(1);
    metadata.Longitude(i) = -abs(Core{i}.Info.Coord1(2));
    metadata.Elevation(i) = Core{i}.Info.Coord1(3);
    metadata.T_avg(i) = Core{i}.Info.T_lt_mean;
    metadata.b_avg(i) = Core{i}.Info.b_lt_mean;
    metadata.Densities{i} = Core{i}.Info.Densities;
    metadata.DepthMax(i) = length(Core{i}.Data.Density)/100;


    
    if strcmp('y', Core{i}.Info.Densities) ...
            && metadata.DepthMax(i)>5
        % pore space at 10 m
        metadata.FAC10(i) = Core{i}.Data.PoreVolume(1000);
       
        % HL porespace at 10 m
        Core{i}.Data.PoreVolume_HL  = cumsum(ps(Core{i}.Data.DensityHL.*0.01,Core{i}.Data.DensityHL));
        metadata.FAC10_HL(i) = Core{i}.Data.PoreVolume_HL(1000);
        
        % finding snow thickness 
        density = Core{i}.Data.Density;
        thickness_weq_mm = 0.01 .* density ;
        depth_weq_mm = cumsum(thickness_weq_mm);
        [~,ind] = min(abs(depth_weq_mm-Core{i}.Info.b_lt_mean));
         metadata.SnowDepth(i) = Core{i}.Data.Depth(ind);
         metadata.FAC10_snow(i) = Core{i}.Data.PoreVolume(ind);
         metadata.FAC10_firn(i) = metadata.FAC10(i) - metadata.FAC10_snow(i);

%          plotting
%          plot(Core{i}.Data.PoreVolume(1:1000),Core{i}.Data.Depth(1:1000))
% hold on
%          plot([metadata.FAC10(i) metadata.FAC10(i)],...
%              Core{i}.Data.Depth([1 1000]))
%          plot([metadata.FAC10_snow(i) metadata.FAC10_snow(i)],...
%              Core{i}.Data.Depth([1 ind]))
%          plot([metadata.FAC10_firn(i) metadata.FAC10_firn(i)],...
%              Core{i}.Data.Depth([ind 1000]))
%          set(gca,'YDir','reverse')
%          ylim([0 1000])
%          xlim([0 8])
%          hold off
%          pause(0.01)

    else
        metadata.FAC10_HL(i) = NaN;
    end
    metadata.Citation{i} = Core{i}.Info.Citation;

end
ind_remove = [];
for i = 1:length(metadata.CoreNumber)
    if isfield(Core{metadata.CoreNumber(i)}.Data,'Density_origin')
        if sum(Core{metadata.CoreNumber(i)}.Data.Density_origin(1:1000))>500
            ind_remove = [ind_remove i];
        end
    end
end
metadata(ind_remove,:) = [];

% metadata(metadata.DepthMax<3) = [];
metadata(strcmp(metadata.Densities,'n'),:) = [];
metadata(isnan(metadata.FAC10),:) = [];
metadata(metadata.DepthMax<5,:) = [];
metadata(strcmp(metadata.Name,'WindSled_2016'),:) = [];
metadata(strcmp(metadata.Name,'core_2015T_A5'),:) = [];
metadata(strcmp(metadata.Name,'62'),:) = [];
metadata(strcmp(metadata.Name,'Inge Lehmann'),:) = [];
metadata(strcmp(metadata.Name,'T19_Spring_2004'),:) = [];
metadata(strcmp(metadata.Name,'62'),:) = [];
metadata(strcmp(metadata.Name,'DYE2 1998 combined'),:) = [];
metadata.Citation{strcmp(metadata.Name,'122')} = 'Steen-Larsen H. C., Masson-Delmotte V., Sjolte J., Johnsen S. J., Vinther B. M., Bre´on F.-M., Clausen H. B., Dahl-Jensen D., Falourd S., Galle´e H., Jouzel J., Kageyama M., Lerche H., Minster B., Picard G., Punge H. J., Risi C., Salas D., Schwander J., Steffen K., Sveinbjo¨rnsdo´ ttir A. E., Svensson A. and White J. (2011) Understanding the climatic signal in the water stable isotope records from the NEEM cores. J. Geophys. Res. 116, D06108. doi:10.1029/2010JD014311.';
metadata.Name{strcmp(metadata.Name,'122')} = 'NEEM07S3';
% metadata(strcmp(metadata.Name,'FA13A'),:) = [];
% metadata(strcmp(metadata.Name,'FA13B'),:) = [];
% metadata(strcmp(metadata.Name,'FA14'),:) = [];
% metadata(strcmp(metadata.Name,'ACT11A2'),:) = [];
% metadata(strcmp(metadata.Name,'ACT11A'),:) = [];

writetable(metadata,'Output/metadata.csv','Delimiter',';');

% cores potentially to exclude
%     SearchCore(Core,'Name','Site'),...
%     FindCore(Core,'Name','H3.5-1'),...
%     SearchCore(Core,'Name','T19_Spring')...
%     SearchCore(Core,'Name','H3-1')...
%     SearchCore(Core,'Name','Benson_1954_10')...

figure('Visible',vis)
subplot(2,3,1)
hist(metadata.Year)
xlabel('Year')
axis tight square
box on
subplot(2,3,2)
hist(metadata.Elevation)
xlabel('Elevation')
axis tight square
box on

subplot(2,3,3)
hist(metadata.DepthMax,60)
xlabel('DepthMax')
axis tight square
box on
xlim([0 25])

subplot(2,3,4)
hist(metadata.T_avg)
xlabel('Average temperature')
axis tight square
box on

subplot(2,3,5)
hist(metadata.b_avg)
xlabel('Average accumulation')
axis tight square
box on

subplot(2,3,6)
hist(metadata.FAC10)
xlabel('Firn air content to 10m')
axis tight square
box on

end