function [metadata] = LoadFAC10Dataset(filename)
    delimiter = ';';
    formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[1,4,5,6,7,8,9,10,11,13,14]
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
            % Handle dates surrounded by quotes
            dataArray{2} = cellfun(@(x) x(2:end-1), dataArray{2}, 'UniformOutput', false);
            dates{2} = datetime(dataArray{2}, 'Format', 'dd-MMM-yyyy', 'InputFormat', 'dd-MMM-yyyy');
        catch
            dates{2} = repmat(datetime([NaN NaN NaN]), size(dataArray{2}));
        end
    end

    anyBlankDates = cellfun(@isempty, dataArray{2});
    anyInvalidDates = isnan(dates{2}.Hour) - anyBlankDates;
    dates = dates(:,2);

    rawNumericColumns = raw(:, [1,4,5,6,7,8,9,10,11,13,14]);
    rawCellColumns = raw(:, [3,12,15]);


    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

    metadata = table;
    metadata.Year = cell2mat(rawNumericColumns(:, 1));
    metadata.Date = dates{:, 1};
    metadata.Name = rawCellColumns(:, 1);
    metadata.Latitude = cell2mat(rawNumericColumns(:, 2));
    metadata.Longitude = cell2mat(rawNumericColumns(:, 3));
    metadata.Elevation = cell2mat(rawNumericColumns(:, 4));
    metadata.T_avg = cell2mat(rawNumericColumns(:, 5));
    metadata.c_avg = cell2mat(rawNumericColumns(:, 6));
    metadata.DepthMax = cell2mat(rawNumericColumns(:, 7));
    metadata.FAC10 = cell2mat(rawNumericColumns(:, 8));
    metadata.FACtot = cell2mat(rawNumericColumns(:, 9));
    metadata.Citation = rawCellColumns(:, 2);
    metadata.FAC10_predict = cell2mat(rawNumericColumns(:, 10));
    metadata.residual = cell2mat(rawNumericColumns(:, 11));

    metadata(1,:) =  [];
    metadata = standardizeMissing(metadata,-9999);
    clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns rawCellColumns R;
end