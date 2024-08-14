function map=load_atlas_mapping(n_nodes,category)

% syntax: map=load_atlas_mapping(n_nodes,category)
% n_nodes are number of nodes: 268 or 278
% category: 'subnetwork' or 'lobe'

load(sprintf('map%0.0f_%s.mat',n_nodes,category));