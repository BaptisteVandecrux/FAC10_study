function [FAC_tt, a] = CalcFACevo(fac_3d, mask, time, lat, lon,res)
FAC_tt = NaN(size(time));

        a = nan(size(fac_3d(:,:,1)));
%         Ee = referenceSphere('earth','km');
        
        % calculating the area of each tile
% %up left
%         lat_up_left = lat;
%         lat_up_left(1,:) = NaN;
%         lat_up_left(:,1) = NaN;
%         lat_up_left(2:end,2:end) = lat(1:end-1,1:end-1);
%         
%         lon_up_left = lon;
%         lon_up_left(1,:) = NaN;
%         lon_up_left(:,1) = NaN;
%         lon_up_left(2:end,2:end) = lon(1:end-1,1:end-1);
% %up right
%         lat_up_right = lat;
%         lat_up_right(1,:) = NaN;
%         lat_up_right(:,end) = NaN;
%         lat_up_right(2:end,1:end-1) = lat(1:end-1,2:end);
%         
%         lon_up_right = lon;
%         lon_up_right(1,:) = NaN;
%         lon_up_right(:,end) = NaN;
%         lon_up_right(2:end,1:end-1) = lon(1:end-1,2:end);
% % down right
%         lat_down_right = lat;
%         lat_down_right(end,:) = NaN;
%         lat_down_right(:,end) = NaN;
%         lat_down_right(1:end-1,1:end-1) = lat(2:end,2:end);
%         
%         lon_down_right = lon;
%         lon_down_right(end,:) = NaN;
%         lon_down_right(:,end) = NaN;
%         lon_down_right(1:end-1,1:end-1) = lon(2:end,2:end);
% %down left
%         lat_down_left = lat;
%         lat_down_left(end,:) = NaN;
%         lat_down_left(:,1) = NaN;
%         lat_down_left(1:end-1,2:end) = lat(2:end,1:end-1);
%         
%         lon_down_left = lon;
%         lon_down_left(end,:) = NaN;
%         lon_down_left(:,1) = NaN;
%         lon_down_left(1:end-1,2:end) = lon(2:end,1:end-1);

%         a1 = areaquad( lat(:), lon(:), lat_up_left(:), lon_up_left(:), Ee );
%         a2 = areaquad( lat(:), lon(:), lat_up_right(:), lon_up_right(:), Ee );
%         a3 = areaquad( lat(:), lon(:), lat_down_left(:), lon_down_left(:), Ee );
%         a4 = areaquad( lat(:), lon(:), lat_down_right(:), lon_down_right(:), Ee );
% 
%         a = reshape(mean([a1, a2, a3, a4], 2), size(lat));
% figure
% hold on
    for i = 1:1:length(time)
        factmp = fac_3d(:,:,i);
        m = and(isfinite(factmp),mask==1);
        factmp(m~=1) = 0;
%         facm = NaN(size(factmp));

        % calculating the FAC
%         facm = a.*1000.*1000.*factmp;
        facm = res^2.*1000.*1000.*factmp;
        FAC_tt(i) = nansum(facm(:)) / 10^9;
%         plot(time(i), FAC_tt(i),'o')
%         datetick('x','dd-mm-yyyy')
%         pause(0.01)

    end
end

%         %% serguei's solution
%         tic        
%         for r = 2:size(factmp,1)-1
%             for c = 2:size(factmp,2)-1
%                 if m(r,c) == 0
%                     continue
%                 else
%                     a1 = areaquad( lat(r,c), lon(r,c), lat(r+1,c+1), lon(r+1,c+1), Ee );
%                     a2 = areaquad( lat(r,c), lon(r,c), lat(r-1,c+1), lon(r-1,c+1), Ee );
%                     a3 = areaquad( lat(r,c), lon(r,c), lat(r+1,c-1), lon(r+1,c-1), Ee );
%                     a4 = areaquad( lat(r,c), lon(r,c), lat(r-1,c-1), lon(r-1,c-1), Ee );
%                     a(r,c) = mean([a1 a2 a3 a4]); clear a1 a2 a3 a4
%                 end
%             end
%         end
%         toc
        