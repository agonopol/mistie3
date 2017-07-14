function [contractor] = cytd(file, varargin)
   
    options = Options();
    options.epsilonClusterIdentificationMethod = 'constantEpsilon';
    options.frequencyMergingEpsilonClusters = 'always'; %always,uponMetastability%
    options.controlSigmaMethod = 'nuclearNormStabilization'; %nuclearNormStabilization,movementStabilization
    options.numDiffusionSteps = 3;
    options.phateEmbedding = true;
    
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'phateEmbedding'))
            options.phateEmbedding = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'controlSigmaMethod'))
            options.controlSigmaMethod = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'frequencyMergingEpsilonClusters'))
            options.frequencyMergingEpsilonClusters = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'epsilonClusterIdentificationMethod'))
            options.epsilonClusterIdentificationMethod = varargin{i+1};
        end
    end

    path = fullfile(file.folder, file.name);
    obj = CyTOFData(path);
    obj.dataTransformed = CyTOFData.transform(obj.data, 1);
    fields = channels(obj);
    data = obj.dataTransformed(:, cell2mat(fields(:,1))');
    [data, index] = datasample(data, min(length(data), 2000),'Replace', false);
   
    [~, name, ~] = fileparts(path);
    options.destination = fullfile(pwd(), 'results', name, '//');
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
    
     contractor = ContractionClustering(data, fields(:,2), options);
end