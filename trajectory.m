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
options.phateEmbedding = true;


d1 = 'data/6300_HF_D1_splorm_0_normalized_Clean.fcs';
d3 = 'data/6300_HF_D3_splorm_0_normalized_Clean.fcs';


d1_cytof = CyTOFData(d1);
d1_cytof.dataTransformed = CyTOFData.transform(d1_cytof.data, 1);
d1_fields = channels(d1_cytof);
d1_data = d1_cytof.dataTransformed(:, cell2mat(d1_fields(:,1))');
[d1_data, ~] = datasample(d1_data, min(length(d1_data), 2000),'Replace', false);

[~, name, ~] = fileparts(d1);
options.destination = fullfile(pwd(), 'results', name, '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
    
d1_cluster = ContractionClustering(d1_data, d1_fields(:,2), options);
d1_cluster = d1_cluster.contract();

d3_cytof = CyTOFData(d3);
d3_cytof.dataTransformed = CyTOFData.transform(d3_cytof.data, 1);
d3_fields = channels(d3_cytof);
d3_data = d3_cytof.dataTransformed(:, cell2mat(d3_fields(:,1))');
[d3_data, ~] = datasample(d3_data, min(length(d3_data), 2000),'Replace', false);

[~, name, ~] = fileparts(d3);
options.destination = fullfile(pwd(), 'results', name, '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
    
d3_cluster = ContractionClustering(d3_data, d3_fields(:,2), options);
d3_cluster = d3_cluster.contract();
   
