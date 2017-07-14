close all;
clear;
clc;

addpath('./lib');
loaddeps();

options = Options();
options.epsilonClusterIdentificationMethod = 'constantEpsilon';
options.frequencyMergingEpsilonClusters = 'always'; %always,uponMetastability%
options.controlSigmaMethod = 'nuclearNormStabilization'; %nuclearNormStabilization,movementStabilization
options.numDiffusionSteps = 3;
options.phateEmbedding = true;


from = 'data/6300_Blood_D1_splorm_0_normalized_Clean.fcs';
to = 'data/6300_Blood_D3_splorm_0_normalized_Clean.fcs';


from_cytof = CyTOFData(from);
from_cytof.dataTransformed = CyTOFData.transform(from_cytof.data, 1);
from_fields = channels(from_cytof);
from_data = from_cytof.dataTransformed(:, cell2mat(from_fields(:,1))');
[from_data, ~] = datasample(from_data, min(length(from_data), 2000),'Replace', false);

[~, name, ~] = fileparts(from);
options.destination = fullfile(pwd(), 'results', name, '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
    
from_cluster = ContractionClustering(from_data, from_fields(:,2), options);
from_cluster = from_cluster.contract();

to_cytof = CyTOFData(to);
to_cytof.dataTransformed = CyTOFData.transform(to_cytof.data, 1);
to_fields = channels(to_cytof);
to_data = to_cytof.dataTransformed(:, cell2mat(to_fields(:,1))');
[to_data, ~] = datasample(to_data, min(length(to_data), 2000),'Replace', false);

[~, name, ~] = fileparts(to);
options.destination = fullfile(pwd(), 'results', name, '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
    
to_cluster = ContractionClustering(to_data, to_fields(:,2), options);
to_cluster = to_cluster.contract();
   
[I,J] = kneepoint(from_cluster, to_cluster);
