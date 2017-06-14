classdef CyTOFData
    properties
        data
        dataTransformed
        name2Channel
        channel2Name
        markerNames
        channelNames
        masses
        analysisDate
    end
    methods
        function obj = CyTOFData(filename)
            [obj.data, header, ~] = fca_readfcs(filename);
            obj.analysisDate = header.date;
            obj.markerNames = {header.par.name2};
            obj.channelNames = {header.par.name};
            obj = obj.computeNameChannelMaps(header);
            obj = obj.populateMasses();
            obj.dataTransformed = [];
        end
        function obj = computeNameChannelMaps(obj, header)
            header_size = length(header.par);
            obj.name2Channel = containers.Map();
            obj.channel2Name = cell(header_size,1);

            for i=1:header_size
                hname = header.par(i).name2;
                
                if(length(strfind(hname, '('))==0)
                    %stuff like the DNA channels do not have the isotope
                    %names
                    if(length(hname)==0)
                        hname = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end
                    hname = lower(hname);
                    obj.name2Channel(hname) = i;
                    obj.channel2Name{i} = hname;
                else
                    n=strfind(hname, '(');
                    
                    if(n==1)
                        %using isotope name
                        n_end = strfind(hname, ')');
                        real_name = hname(n+1:n_end-1);
                    else
                        real_name = hname(1:n-1);
                    end
                    
                    if(length(real_name)==0)
                        real_name = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end

                    real_name = lower(real_name);
                    obj.name2Channel(real_name) = i;
                    obj.channel2Name{i} = real_name;
                end
            end
        end
        function obj = populateMasses(obj)
            obj.masses = NaN(size(obj.data, 2), 1);
            for i=1:length(obj.masses)
                channelName = obj.channelNames{i};
                startIndex = regexp(channelName, '[\d]+Di');
                if (~isempty(startIndex))
                    obj.masses(i) = str2num(channelName(startIndex:end-2));
                end
            end
        end
        function obj = addClusterAssigments(obj, index, assigment)
           obj.channelNames{length(obj.channelNames) + 1} = 'Cluster';
           obj.markerNames{length(obj.markerNames) + 1} = 'Cluster';
           obj.data(index, end+1) = assigment;
        end
        function writeData(obj, filename)
            fca_writefcs(filename, obj.data, obj.markerNames, obj.channelNames);
        end
        function obj = filterLive(obj, thresholdDead, emitPlot, prefixFilenamePlot, varargin)
            disp('Filtering live cells ...');
            channelDeadLive = 'Pt195Di';

            % parsing of variable argument list
            for i=1:length(varargin)-1
                if (strcmp(varargin{i}, 'maxCisplatin'))
                    maxCisplatin = varargin{i+1};
                end
                if (strcmp(varargin{i}, 'maxFrequency'))
                    maxFrequency = varargin{i+1};
                end
                if (strcmp(varargin{i}, 'channelDeadLive'))
                    channelDeadLive = varargin{i+1};
                end
            end
            columnDeadLive = find(ismember(obj.channelNames, channelDeadLive));

            if (~exist('maxCisplatin'))
                maxCisplatin = max(obj.dataTransformed(:, columnDeadLive));
            end

            if (~isempty(obj.dataTransformed))
                cisplatinValues = obj.dataTransformed(:, [columnDeadLive]);
            else
                cisplatinValues = CyTOFData.transform(obj.data(:, [columnDeadLive], 1));
            end

            indicesLive = cisplatinValues < thresholdDead;

            numEvents = size(obj.data, 1);
            numEventsAfter = sum(indicesLive);
            percentageLive = round(100*numEventsAfter/numEvents);

            if (emitPlot)
                figure
                thresholdHistogram(cisplatinValues, ...
                                   thresholdDead, ...
                                   'colorAssignment', [1 0 0; 0 0 0])
                xlim([0 maxCisplatin]);
                if (exist('maxFrequency'))
                    ylim([0, maxFrequency]);
                end
                set(gcf, 'Color', 'w');
                set(gcf, 'Position', [2560 1 400 300]);
                xlabel('Cisplatin (195pt)');
                ylabel('# Cells');
                set(gca, 'fontsize', 18);
                legend({['dead (' num2str(100-percentageLive) '%)'], ...
                        ['live (' num2str(percentageLive) '%)']});
                if (exist('export_fig', 'file'))
                    export_fig(strcat(prefixFilenamePlot, '_histo_deadlive.png'));
                end
            end

            obj = obj.filterCells(indicesLive);

            disp(['    removed ' num2str(numEvents-numEventsAfter) ' from ' num2str(numEvents) ...
                       ' (' num2str(100-percentageLive) ...
                       '%).']);
        end
        function obj = filterSingleNucleated(obj, intervalDNA1, intervalDNA2, emitPlots, prefixFilenamePlots, varargin)
            disp('Filtering single nucleated cells ...');
            columnDNA1 = find(ismember(obj.channelNames, 'Ir191Di'));
            columnDNA2 = find(ismember(obj.channelNames, 'Ir193Di'));

            % parsing of variable argument list
            for i=1:length(varargin)-1
                if (strcmp(varargin{i}, 'maxDNA1'))
                    maxDNA1 = varargin{i+1};
                end
                if (strcmp(varargin{i}, 'maxDNA2'))
                    maxDNA2 = varargin{i+1};
                end
            end

            if (~exist('maxDNA1'))
                maxDNA1 = max(obj.data(:, columnDNA1));
            end
            if (~exist('maxDNA2'))
                maxDNA2 = max(obj.data(:, columnDNA2));
            end

            if (~isempty(obj.dataTransformed))
                intercalationValues = obj.dataTransformed(:, [columnDNA1 columnDNA2]);
            else
                intercalationValues = CyTOFData.transform(obj.data(:, [columnDNA1 columnDNA2], 1));
            end
            indicesSingleNucleated =   (intercalationValues(:, 1) > intervalDNA1(1)) ...
                                     & (intercalationValues(:, 2) > intervalDNA2(1)) ...
                                     & (intercalationValues(:, 1) < intervalDNA1(2)) ...
                                     & (intercalationValues(:, 2) < intervalDNA2(2));
            numEvents = size(obj.data, 1);
            numEventsAfter = sum(indicesSingleNucleated);
            percentageSingleNucleated = round(100*numEventsAfter/numEvents);
            if (emitPlots)
                %overall
                transparency = 1;
                if (numEvents > 150000)
                    transparency = 0.01;
                elseif (numEvents > 50000)
                    transparency = 0.05;
                end

                figure
                hold on
                scatterX([NaN, NaN], 'colorAssignment', [0 0 0], 'transparency', 1);
                scatterX([NaN, NaN], 'colorAssignment', [1 0 0], 'transparency', 1);
                legend({strcat('single nucleated (', num2str(percentageSingleNucleated), '%)'), ...
                        strcat('rest (', num2str(100-percentageSingleNucleated), '%)')}, ...
                       'location', 'northwest');
                scatterX(intercalationValues(indicesSingleNucleated, :), ...
                         'colorAssignment', repmat([0 0 0], sum(indicesSingleNucleated), 1), ...
                         'sizeAssignment', 0.1*ones(sum(indicesSingleNucleated), 1), ...
                         'transparency', transparency, ...
                         'addJitter', false);
                scatterX(intercalationValues(~indicesSingleNucleated, :), ...
                         'colorAssignment', repmat([1 0 0], sum(~indicesSingleNucleated), 1), ...
                         'sizeAssignment', 0.1*ones(sum(~indicesSingleNucleated), 1), ...
                         'transparency', transparency, ...
                         'addJitter', false);
                contourX(intercalationValues)
                set(gcf, 'Color', 'w');
                set(gcf, 'Position', [2068 1 700 700]);
                set(gca, 'fontsize', 18);
                xlabel('DNA_1 (191ir)');
                ylabel('DNA_2 (193ir)');
                xlim([0 maxDNA1]);
                ylim([0 maxDNA2]);
                if (exist('export_fig', 'file'))
                    export_fig(strcat(prefixFilenamePlots, '_intercalation.png'));
                end
                hold off

                %zoomed in
                figure
                hold on
                scatterX(intercalationValues(indicesSingleNucleated, :), ...
                         'colorAssignment', repmat([0 0 0], sum(indicesSingleNucleated), 1), ...
                         'sizeAssignment', 1*ones(sum(indicesSingleNucleated), 1), ...
                         'transparency', transparency, ...
                         'addJitter', false);
                scatterX(intercalationValues(~indicesSingleNucleated, :), ...
                         'colorAssignment', repmat([1 0 0], sum(~indicesSingleNucleated), 1), ...
                         'sizeAssignment', 1*ones(sum(~indicesSingleNucleated), 1), ...
                         'transparency', transparency, ...
                         'addJitter', false);
                contourX(intercalationValues)
                set(gcf, 'Color', 'w');
                set(gcf, 'Position', [2068 1 700 700]);
                set(gca, 'fontsize', 18);
                xlabel('DNA_1 (191ir)');
                ylabel('DNA_2 (193ir)');
                xslack = (intervalDNA1(2)-intervalDNA1(1))/2;
                yslack = (intervalDNA2(2)-intervalDNA1(2))/2;
                xlim([intervalDNA1(1)-xslack, intervalDNA1(2)+xslack]);
                ylim([intervalDNA2(1)-yslack, intervalDNA2(2)+yslack]);
                if (exist('export_fig', 'file'))
                    export_fig(strcat(prefixFilenamePlots, '_intercalation_detail.png'));
                end
                hold off
            end

            obj = obj.filterCells(indicesSingleNucleated);

            disp(['    removed ' num2str(numEvents-numEventsAfter) ' from ' num2str(numEvents) ...
                       ' (' num2str(100-percentageSingleNucleated) ...
                       '%).']);
        end
        function [obj, filteredClusterAssignment] = filterSmallClusters(obj, clusterAssignment, cutoff)
            % remove clusters smaller than cutoff data points
            clusters = unique(clusterAssignment);
            for i=1:length(clusters)
                if (sum(clusterAssignment == clusters(i)) < cutoff)
                    obj = obj.filterCells(clusterAssignment ~= clusters(i));
                    clusterAssignment = clusterAssignment(clusterAssignment ~= clusters(i));
                end
            end

            % regularlize (cluster assignments are from 1 to #clusters)
            remainingClusters = unique(clusterAssignment);
            for i=1:length(remainingClusters)
                clusterAssignment(clusterAssignment == remainingClusters(i)) = i;
            end

            filteredClusterAssignment = clusterAssignment;
        end
        function obj = transformAllIsotopeChannels(obj, c)
            obj.dataTransformed = zeros(size(obj.data));
            isotopeChannels = [];
            obj.masses = NaN(size(obj.data, 2), 1);
            for i=1:size(obj.channelNames, 2)
                channelName = obj.channelNames{i};
                isotopeChannelName = regexp(channelName, '[\d]+Di');
                if (~isempty(isotopeChannelName))
                    isotopeChannels = [isotopeChannels i];
                end
            end
            obj.dataTransformed(:, isotopeChannels) = CyTOFData.transform(obj.data(:, isotopeChannels), c);
        end
        function result = numCells(obj)
            assert(isempty(obj.dataTransformed) || size(obj.dataTransformed, 1) == size(obj.data, 1));
            result = size(obj.data, 1);
        end
        function obj = filterCells(obj, indices)
            obj.data = obj.data(indices, :);
            if (~isempty(obj.dataTransformed))
                obj.dataTransformed = obj.dataTransformed(indices, :);
            end
        end
        function obj = concat(obj, obj2)
            obj.data = [obj.data ; obj2.data];
            if (~isempty(obj.dataTransformed))
                assert(~isempty(obj2.dataTransformed))
                obj.dataTransformed = [obj.dataTransformed ; obj2.dataTransformed];
            end
        end
    end
    methods(Static)
        function transformedData = transform(data, c)
            transformedData = zeros(size(data));
            for i=1:size(data, 2)
                transformedData(:, i) = asinh(data(:, i)/c);
            end
        end
    end
end
