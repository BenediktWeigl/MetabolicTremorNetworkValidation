% Compute spatial agreement and disagreement maps between a tremor network
% and a PET-derived statistical map.

clear; clc;

% --- File paths ---
folder = '/Volumes/Weigl_lab/02_LabProjects/Goede_Tremor/04_Network_Spatial_similarity/Voxelwise_Similarity';
network_file = fullfile(folder, 'alligned_mulitmodal_tremormap_zscored_AveragedPair.nii');
pet_file     = fullfile(folder, 'simpleTDiff_swbTremor_4D_vox_tstat_c1.nii');
output_info  = niftiinfo(pet_file);

out_agree = fullfile(folder, 'Z-Agreement_Map_Sign.nii');
out_disagree = fullfile(folder, 'Z-Disagreement_Map_PET-Sign.nii');
out_petonly = fullfile(folder, 'Z-PET_Only_Map.nii');

% --- Load images ---
network = double(niftiread(network_file));
pet     = double(niftiread(pet_file));

if ~isequal(size(network), size(pet))
    error('Image dimensions do not match.');
end

% Thresholding to account for numbers near zero after reslicing
threshold = 1e-5; 
if threshold > 0
    network(abs(network) < threshold) = 0;
    pet(abs(pet) < threshold) = 0;
end

% --- Compute agreement/disagreement ---
prod = network .* pet;
sign_net = sign(network);

agree = zeros(size(pet));
agree(prod > 0) = prod(prod > 0) .* sign_net(prod > 0);

disagree = zeros(size(pet));
disagree(prod < 0) = prod(prod < 0) .* sign_net(prod < 0);

petonly = zeros(size(pet));
petonly((pet ~= 0) & (network == 0)) = pet((pet ~= 0) & (network == 0));

% --- Z-score non-zero voxels ---
zscore_nonzero = @(x) (x ~= 0) .* ((x - mean(x(x ~= 0))) / std(x(x ~= 0)));
agree = zscore_nonzero(agree);
disagree = zscore_nonzero(disagree);
petonly = zscore_nonzero(petonly);

% --- Save output maps ---
niftiwrite(single(agree),    out_agree, output_info);
niftiwrite(single(disagree), out_disagree, output_info);
niftiwrite(single(petonly),  out_petonly, output_info);
