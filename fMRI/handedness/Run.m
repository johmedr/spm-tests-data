%% Analysis of handedness data using Bayesian Model Reduction
%
% This dataset includes 36 subjects, 18 left and 18 right handed. Subjects 
% clenched their left or right fists when cued at variable rate. The 
% left-handed subjects are flipped in the x-axis, so all subjects have the 
% following 6 regions:
%
% 1. Dominant SMA  
% 2. Dominant ventral pre-motor (PM)  
% 3. Dominant primary motor (M1)
% 4. Non-dominant SMA
% 5. Non-dominant PM
% 6. Non-dominant M1
%
% Each model has one driving input (Task) and four modulatory inputs:
%
% 1. Dominant fist clench
% 2. Dominant fist clench speed (parametric)
% 3. Non-dominant fist clench
% 4. Non-dominant fist clench (parametric)
%
% Subject metadata including handedness is stored in subjects.mat
%
% Hypotheses: 
% 1. Clenching the dominant hand should increase connectivity to dominant M1.
% 2. Inhibitory inter-hemispheric connections should be stronger for
%    right handed participants
%
% Data from Pool et al., NeuroImage 2014.

%% Settings
base_dir = pwd;

% Subject metadata
load('subjects.mat');
subject_ids = subjects.id;

% Dimensions
n = length(subject_ids);
m = 1; 

% Regions: d = dominant, nd = nondominant
region_names = {'dSMA','ndSMA','dPM','ndPM','dM1','ndrM1'};

% Design matrix for group analyses. Right handed = -1, left handed = 1
HANDEDNESS=2; GENDER=3;
X = ones(n,3); 
X(subjects.handedness == 1, HANDEDNESS) = -1;
X(subjects.gender == 1, GENDER) = -1;

%% Build subjects x models array of DCM filenames
models_dir = fullfile(base_dir,'models_spm12_6regions_1input_leftflipped');

P = cell(n,m);
for s = 1:n
    id = subject_ids(s);
    for model = 1:m
        filename = sprintf('DCM_s%d_model_%d.mat',id,model);
        P{s,model} = fullfile(models_dir, filename);
    end
end

RCM = spm_dcm_load(P);
%% Average model across all subjects
BMA = spm_dcm_bma(RCM);
save('BMA_flipped.mat','BMA');
%% Display average model using heat maps
figure;
subplot(2,2,1); imagesc(BMA.Ep.B(:,:,1));
set(gca,'XTickLabel',region_names); set(gca,'YTickLabel',region_names);
colorbar; colormap(bipolar); title('Dominant fist clench');

subplot(2,2,2); imagesc(BMA.Ep.B(:,:,3));
set(gca,'XTickLabel',region_names); set(gca,'YTickLabel',region_names);
colorbar; colormap(bipolar); title('Non-dominant fist clench');

subplot(2,2,3); imagesc(BMA.Ep.B(:,:,2));
set(gca,'XTickLabel',region_names); set(gca,'YTickLabel',region_names);
colorbar; colormap(bipolar); title('Dominant hand clench speed');

subplot(2,2,4);imagesc(BMA.Ep.B(:,:,4));
set(gca,'XTickLabel',region_names); set(gca,'YTickLabel',region_names);
colorbar; colormap(bipolar); title('Non-dominant hand clench speed');
%% Perform CVA (no empirical Bayes yet)
% Build subjects x params matrix
Q = [];
for i = 1:n
    Q(i,:)=spm_vec(BMA.SUB(i).Ep);
end

% CVA
CVA = spm_cva(Q,X,[],[0 1 0]');
figure,bar([CVA.v X(:,HANDEDNESS)]); legend({'Prediction','Data'});
title(sprintf('CVA (no empirical Bayes) p = %2.2f',CVA.p));
%% Run empirical Bayes: use group membership to optimise 1st level
M = struct('X',X);
[peb,PCM] = spm_dcm_peb(RCM,M,'all');
BMA_PEB = spm_dcm_bma(PCM);

% Build subjects x params matrix
Q = [];
for i = 1:n
    Q(i,:)=spm_vec(BMA_PEB.SUB(i).Ep);
end

% CVA
CVA = spm_cva(Q,X,[],[0 1 0]');
figure,bar([CVA.v X(:,HANDEDNESS)]); legend({'Prediction','Data'});
title(sprintf('CVA with empirical Bayes p = %2.2f',CVA.p));

assert(CVA.p < 0.05, 'Handedness was non-significant');

% Convert CVA weights to DCM format
CVA_params = spm_unvec(CVA.V, RCM{1,1}.Ep);
%% Peform model reduction at the second level - should retain handedness

% Build design matrix. Right handed = -1, left handed = 1
M = struct('X',X);
[BMC,PEB] = spm_dcm_bmc_peb(RCM,M);

% Average
[BMA] = spm_dcm_peb_bmc(PEB);
save('BMA_PEB.mat', 'BMA', 'PEB', 'BMC');

disp('Second level models attempted:');
disp(full(BMC.K));
