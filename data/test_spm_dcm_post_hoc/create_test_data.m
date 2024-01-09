% Script to create the test data for test_spm_dcm_post_hoc. Adapted
% from the example script available on the SPM website. Creates a
% generative model with simulated data based on empirical parameter
% estimates (DCM_fwd_simulated.mat) and a model fitted to the same data
% with a redundant parameter (DCM_fwd_bwd_simulated.mat).
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Peter Zeidman
% $Id: dcm_spm12_batch.m 12 2014-09-29 19:58:09Z guillaume $


% Directory containing the attention data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
fprintf('%-40s:', 'Downloading Attention dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/attention/attention.zip','attention.zip');
unzip(fullfile(data_path,'attention.zip'));
data_path = fullfile(data_path,'attention');
fprintf(' %30s\n', '...done');

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION & INFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factors = load(fullfile(data_path,'factors.mat'));

f = spm_select('FPList', fullfile(data_path,'functional'), '^snf.*\.img$');

clear matlabbatch

% OUTPUT DIRECTORY
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% MODEL SPECIFICATION
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT    = 3.22;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans            = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name     = 'Photic';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset    = [factors.att factors.natt factors.stat];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name     = 'Motion';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset    = [factors.att factors.natt];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name     = 'Attention';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset    = [factors.att];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 10;

% MODEL ESTIMATION
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

% INFERENCE
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(3);
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Photic';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [1 0 0];
matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'Motion';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 1 0];
matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'Attention';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 0 1];

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLUMES OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

% EXTRACTING TIME SERIES: V5
%--------------------------------------------------------------------------
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{1}.spm.util.voi.session = 1; % session 1
matlabbatch{1}.spm.util.voi.name = 'V5';
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 3;  % "Motion" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.contrast = 4; % "Attention" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.mtype = 0; % inclusive
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [-36 -87 -3];
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: V1
%--------------------------------------------------------------------------
matlabbatch{2}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{2}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{2}.spm.util.voi.session = 1; % session 1
matlabbatch{2}.spm.util.voi.name = 'V1';
matlabbatch{2}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{2}.spm.util.voi.roi{1}.spm.contrast = 2;  % "Photic" T-contrast
matlabbatch{2}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{2}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{2}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.centre = [0 -93 18];
matlabbatch{2}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{2}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: SPC
%--------------------------------------------------------------------------
matlabbatch{3}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{3}.spm.util.voi.session = 1; % session 1
matlabbatch{3}.spm.util.voi.name = 'SPC';
matlabbatch{3}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{3}.spm.util.voi.roi{1}.spm.contrast = 4;  % "Attention" T-contrast
matlabbatch{3}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
matlabbatch{3}.spm.util.voi.roi{1}.spm.thresh = 0.001;
matlabbatch{3}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.centre = [-27 -84 36];
matlabbatch{3}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

spm_jobman('run',matlabbatch);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC CAUSAL MODELLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear DCM

% SPECIFICATION DCMs "attentional modulation of backward/forward connection"
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

load(fullfile(data_path,'GLM','SPM.mat'));

% Load regions of interest
%--------------------------------------------------------------------------
load(fullfile(data_path,'GLM','VOI_V1_1.mat'),'xY');
DCM.xY(1) = xY;
load(fullfile(data_path,'GLM','VOI_V5_1.mat'),'xY');
DCM.xY(2) = xY;
load(fullfile(data_path,'GLM','VOI_SPC_1.mat'),'xY');
DCM.xY(3) = xY;

DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

% Time series
%--------------------------------------------------------------------------
DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

% Experimental inputs
%--------------------------------------------------------------------------
DCM.U.dt   =  SPM.Sess.U(1).dt;
DCM.U.name = [SPM.Sess.U.name];
DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...
              SPM.Sess.U(2).u(33:end,1) ...
              SPM.Sess.U(3).u(33:end,1)];

% DCM parameters and options
%--------------------------------------------------------------------------
DCM.delays = repmat(SPM.xY.RT/2,DCM.n,1);
DCM.TE     = 0.04;

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.nograph    = 1;

%--------------------------------------------------------------------------
% Modified for testing post_hoc DCM from here
%--------------------------------------------------------------------------

% Connectivity matrices for model with forward (fwd) modulation
%--------------------------------------------------------------------------
DCM.a = [1 1 0; 1 1 1; 0 1 1];
DCM.b = zeros(3,3,3);  DCM.b(2,1,2) = 1;  DCM.b(2,1,3) = 1; % Fwd
DCM.c = [1 0 0; 0 0 0; 0 0 0];
DCM.d = zeros(3,3,0);
save(fullfile(data_path,'GLM','DCM_mod_fwd.mat'),'DCM');

% DCM Estimation
%--------------------------------------------------------------------------
clear matlabbatch

matlabbatch{1}.spm.dcm.fmri.estimate.dcmmat = ...
    {fullfile(data_path,'GLM','DCM_mod_fwd.mat')};

spm_jobman('run',matlabbatch);

% Simulate data from forward (generative) model
%--------------------------------------------------------------------------
hE              = 2; % Log precision of the observation noise
noise_variance  = 1 / exp(hE);
noise_sd        = repmat(noise_variance ^ 2,1,3);
spm_dcm_create('generative','DCM_mod_fwd.mat',noise_sd);

% Estimate models with simulated data
%--------------------------------------------------------------------------

% Forward (generative) model
load('DCM_generative.mat');
DCM = rmfield(DCM,'M');
DCM = spm_dcm_estimate(DCM);
save('DCM_fwd_simulated.mat','DCM');

% Create model with extra (redundant) parameter
DCM.b(2,3,3) = 1; % Bwd
DCM = rmfield(DCM,'M');
DCM = spm_dcm_estimate(DCM);
save('DCM_fwd_bwd_simulated.mat','DCM');

spm_dcm_compare({'DCM_fwd_simulated.mat','DCM_fwd_bwd_simulated.mat'});

% Copy models to test data directory
output_path = fileparts(mfilename('fullpath'));
copyfile('DCM_fwd_simulated.mat',output_path);
copyfile('DCM_fwd_bwd_simulated.mat',output_path);