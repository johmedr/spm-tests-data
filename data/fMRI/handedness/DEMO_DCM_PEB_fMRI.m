
% Analysis of handedness data using Bayesian Model Reduction
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
 
% Settings and metadata
%==========================================================================
base_dir    = pwd;
load('subjects.mat');
subject_ids = subjects.id;
n  = length(subject_ids);
m  = 1;
 
% Regions: d = dominant, nd = non-dominant
%--------------------------------------------------------------------------
region_names = {'dSMA','ndSMA','dPM','ndPM','dM1','ndrM1'};
 
% Build array of DCM filenames
%--------------------------------------------------------------------------
models_dir = fullfile(base_dir,'models_spm12_6regions_1input_leftflipped');
 
P     = cell(n,m);
for s = 1:n
    id = subject_ids(s);
    for model = 1:m
        filename = sprintf('DCM_s%d_model_%d.mat',id,model);
        P{s,model} = fullfile(models_dir, filename);
    end
end
 
 
% second level model specification with between subject effects;
%==========================================================================
% here, handedness and gender:
% Design matrix for group analyses. Right handed = -1, left handed = 1
%--------------------------------------------------------------------------
HANDEDNESS = 2; 
GENDER     = 3;
X          = ones(n,3); 
X(subjects.handedness == 1, HANDEDNESS) = -1;
X(subjects.gender     == 1, GENDER)     = -1;
 
 
%  invert DCMs:
%==========================================================================
%  This inversion scheme is a group inversion scheme that suppresses local
%  minima by applying empirical shrinkage priors recursively to subject
%  specific fits. The resulting posteriors are then adjusted to give
%  what one would have obtained under the original priors.
%--------------------------------------------------------------------------
 
% first (within subject) level inversion
%--------------------------------------------------------------------------
M.X  = X(:,1);          % preclude any bias during group inversion
M.Q  = 'single';        % assume all parameters share a variance ratio
M.hE = 0;               % prior expectation of between subject precision
M.hC = 1/32;            % prior variance of between subject precision
 
[RCM,REB,M,HCM] = spm_dcm_peb_fit(P,M);
 
save workspace
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
load workspace
 
% show results: free energy over iterations of empirical Bayes
%--------------------------------------------------------------------------
spm_figure('GetWin','Inversion'); clf
 
subplot(3,3,1), bar(RCM{1}.FEB(:) - RCM{1}.FEB(1))
xlabel('iteration'), ylabel('(second level) free energy')
title('Free energy','FontSize',16)
axis square
 
subplot(3,3,2), bar(RCM{1}.EEB(:),'b')
xlabel('iteration'), ylabel('correlation')
title('log precision','FontSize',16)
axis square
 
subplot(3,3,3), bar(RCM{1}.HEB(:) - RCM{1}.HEB(1),'c')
xlabel('iteration'), ylabel('correlation')
title('Posterior uncertainty','FontSize',16)
axis square
 
% Parameter estimates
%--------------------------------------------------------------------------
clear Q
for i = 1:n
    for j = 1:size(HCM,2)
        Q(:,i,j) = full(spm_vec(HCM{i,j}.Ep));
    end
end
 
i     = spm_find_pC(RCM{1},[],{'A','B'});             % parameters to plot                 
[q,j] = sort(mean(Q(i,:,end),2),'descend'); i = i(j); % rank: descending mean
 
subplot(3,1,2), plot(Q(i,:,1),'.','MarkerSize',8)
xlabel('parameter'), ylabel('estimate')
title('First iteration','FontSize',16), set(gca,'YLim',[-2 2])
 
subplot(3,1,3), plot(Q(i,:,end),'.','MarkerSize',8)
xlabel('parameter'), ylabel('estimate')
title('Last iteration','FontSize',16), set(gca,'YLim',[-2 2])
 
% sort connections into intrinsic, inter-hemispheric and intra-hemispheric
%--------------------------------------------------------------------------
i      = spm_find_pC(RCM{1},[],{'A'});
Np     = numel(i);
Afield = spm_fieldindices(RCM{1}.Ep,i);
A      = [
    0 1 1 2 2 2
    1 0 1 2 2 2
    1 1 0 2 2 2
    2 2 2 0 1 1
    2 2 2 1 0 1
    2 2 2 1 1 0];
 
for i = 1:length(Afield), Ai(i,1) = eval(Afield{i}); end
field      = Afield(find(Ai == 1));              % Intra-hemispheric
field      = Afield(find(Ai == 2));              % Inter-hemispheric
field      = Afield(find(Ai == 0));              % Intrinsic
model      = sparse(Ai + 2,1:Np,1);              % simple model space
model(1,:) = 1;                                  % add full model
 
% analysis
%==========================================================================
 
% restore group effects in design matrix
%--------------------------------------------------------------------------
M.X  = X;
M.hE = 0;
M.hC = 1/16;
 
% Omnibus analysis searching for any effects of handedness at the group
% level on fixed connectivity (A). This searches over second
% level models, defined here by the indicator variables in 'model'
%--------------------------------------------------------------------------
BMA  = spm_dcm_peb_bmc(spm_dcm_peb(RCM,M,'A'),model);
 
% One can also perform an exhaustive search over (intrinsic) parameters at 
% boths levels. Here, there is only one (full) first level model.
%--------------------------------------------------------------------------
[BMC,PEB] = spm_dcm_bmc_peb(RCM,M,field);
BMA       = spm_dcm_peb_bmc(PEB);
 
 
 
% auxiliary analyses
%==========================================================================
 
% Confirm the effects of (intrinsic) parameters using cross-validation and
% a leave-one-out scheme
%--------------------------------------------------------------------------
spm_dcm_loo(RCM,M,field);
 
% randomisation analysis
%--------------------------------------------------------------------------
spm_dcm_peb_rnd(RCM,M,field);
 
% randomisation analysis - line search over hyperpriors
%--------------------------------------------------------------------------
% spm_dcm_peb_rnd_search(RCM,M,field)




