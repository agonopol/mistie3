close all;
clear;
clc;

addpath('./lib');
loaddeps();

options = Options();
options.clusterAssignmentMethod = 'none';
options.frequencyMergingEpsilonClusters = 'uponMetastability'; %always,uponMetastability%
options.controlSigmaMethod = 'nuclearNormStabilization';
options.numDiffusionSteps = 3;
options.fastStop = true;
options.phateEmbedding = true;

files = dir('data/*.fcs');

for file = files'
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
    contractor = contractor.contract();
    contractor.heatmap();
    
    obj = obj.addClusterAssigments(index, contractor.clusterAssignments(end, :));
    obj.writeData(strrep(path, 'data', 'clustered'));

    clc;
    close all force;
    close all hidden;
end

