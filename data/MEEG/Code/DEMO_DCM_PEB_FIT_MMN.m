% function DEMO_DCM_PEB_FIT
% Test routine to check group DCM for electrophysiology
%--------------------------------------------------------------------------
% This routine illustrates the use of Bayesian model reduction when
% inverting hierarchical (dynamical) models; for example, multisubject DCM
% models. In this demonstration empirical Bayesian model reduction is
% applied recursively in an attempt to finesse the local  minimum problem
% with a nonlinear DCMs. The estimates are compared against standard
% Bayesian model reduction, in terms of the subject specific estimates and
% Bayesian model comparison  (and averages) at the between subject level.
% % 
% This demo considers a single goup (e.g., of normal subjects) and the
% differences between the group average using emprical Bayesian reduction  
% and the Bayesian reduction of the (grand) average response.
%
% See also: DEMO_DCM_PEB_REC.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_DCM_PEB_FIT_MMN.m 16 2015-07-22 11:31:36Z karl $


% change to directory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------

DCM_DIR{1} = 'C:\home\spm\MMN\EEG_MMN_GROUP\odd\ERP_5s';
DCM_DIR{2} = 'C:\home\spm\MMN\EEG_MMN_GROUP\even\ERP_5s';
DCM_DIR{3} = 'C:\home\spm\MMN\EEG_MMN_GROUP\null\ERP_5s';
DCM_DIR{4} = 'C:\home\spm\MMN\EEG_MMN_GROUP\odd\CMC_5s';
DCM_DIR{5} = 'C:\home\spm\MMN\EEG_MMN_GROUP\even\CMC_5s';
DCM_DIR{6} = 'C:\home\spm\MMN\EEG_MMN_GROUP\null\CMC_5s';
cd(DCM_DIR{1})

close all, clear all
clc
spm eeg

rng('default')
corr = @(x,y) subsref(corrcoef(x,y),substruct('()',{1,2})); % Stats tbx

MODEL = 'ERP';

% Set up
%==========================================================================
DCM = [];

DCM.options.trials   = [1 2];
DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'ERP';
DCM.options.model    = MODEL;
DCM.options.Nmax     = 64;
DCM.options.Tdcm     = [0 300];
DCM.options.onset    = 80;
DCM.options.dur      = 16;
DCM.options.Nmodes   = 8;
DCM.options.D        = 1;
DCM.options.h        = 1;
DCM.options.DATA     = 1;
DCM.name             = 'DCM_GROUP';

DCM.xU.X             = [0; 1];
DCM.xU.name          = {'odd'};


DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};

Nareas    = size(DCM.Lpos,2);

DCM.A{1} = zeros(Nareas,Nareas);
DCM.A{1} = zeros(Nareas, Nareas);
DCM.A{1}(3,1) = 1;
DCM.A{1}(4,2) = 1;
DCM.A{1}(5,4) = 1;

DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{2}(1,3) = 1;
DCM.A{2}(2,4) = 1;
DCM.A{2}(4,5) = 1;

DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(4,3) = 1;
DCM.A{3}(3,4) = 1;

DCM.C = [1 1 0 0 0]';

% model space - within subject effects
%--------------------------------------------------------------------------
b{1}  = [1 0 0 0 0;
         0 1 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0]; % intrinsic lower
b{2}  = [0 0 0 0 0;
         0 0 0 0 0;
         0 0 1 0 0;
         0 0 0 1 0;
         0 0 0 0 0]; % intrinsic higher
b{3}  = [0 0 0 0 0;
         0 0 0 0 0;
         1 0 0 0 0;
         0 1 0 0 0;
         0 0 0 1 0]; % forward
b{4}  = [0 0 1 0 0;
         0 0 0 1 0;
         0 0 0 0 0;
         0 0 0 0 1;
         0 0 0 0 0]; % and backward
     
k     = spm_perm_mtx(length(b));
for i = 1:size(k,1)
    B{i}  = zeros(size(b{1}));
    for j = 1:size(k,2)
        B{i} = B{i} + k(i,j)*b{j};
    end
end

% model space
%--------------------------------------------------------------------------
switch MODEL
    case{'ERP'}        
        DCM.options.onset = 80;
    case{'CMC'} 
        DCM.options.onset = 100;
        
    otherwise
        
end


ff = cellstr(spm_select('List', '.', '^fm.*.dat'));

% model space
%--------------------------------------------------------------------------
Nm  = length(B);                      % number of models
Ns  = length(ff);                     % number of subjects

% set up DCM structure array
%--------------------------------------------------------------------------
for i = 1:Ns
    
    % set up data and spatial models
    %----------------------------------------------------------------------
    DCMs = DCM;
    DCMs.xY.Dfile = spm_file(ff{i}, 'ext', 'mat');
    DCMs = spm_dcm_erp_data(DCMs);
    DCMs = spm_dcm_erp_dipfit(DCMs);
    
    
    % parametric models
    %----------------------------------------------------------------------
    for j = 1:Nm
        GCM{i,j}          = rmfield(DCMs,'M');
        GCM{i,j}.M.dipfit = DCMs.M.dipfit;
        GCM{i,j}.B        = B(j);
        GCM{i,j}.name     = ['DCM_S' num2str(i) '_m' num2str(j)];
    end
end


DCMs = DCM;
DCMs.xY.Dfile = spm_select('List', '.', '^MMN_GA.*.mat');
DCMs = spm_dcm_erp_data(DCMs);
DCMs = spm_dcm_erp_dipfit(DCMs);
    
for j = 1:Nm
    gcm{1,j}          = rmfield(DCMs,'M');
    gcm{1,j}.M.dipfit = DCMs.M.dipfit;
    gcm{1,j}.B        = B(j);
    gcm{1,j}.name     = 'DCM_GA';
end



%% inversion of grand average and emprical Bayesian inversion
%==========================================================================
gcm(1)          = spm_dcm_fit(gcm(1));
[RCM,REB,M,HCM] = spm_dcm_peb_fit(GCM);


% Comparison
%==========================================================================

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

save workspace_group
%--------------------------------------------------------------------------

%% load inverted DCMs
%--------------------------------------------------------------------------
DCM_DIR{1} = 'C:\home\spm\MMN\EEG_MMN_GROUP\odd\ERP_5s';
DCM_DIR{2} = 'C:\home\spm\MMN\EEG_MMN_GROUP\even\ERP_5s';
DCM_DIR{3} = 'C:\home\spm\MMN\EEG_MMN_GROUP\null\ERP_5s';
DCM_DIR{4} = 'C:\home\spm\MMN\EEG_MMN_GROUP\odd\CMC_5s';
DCM_DIR{5} = 'C:\home\spm\MMN\EEG_MMN_GROUP\even\CMC_5s';
DCM_DIR{6} = 'C:\home\spm\MMN\EEG_MMN_GROUP\null\CMC_5s';

cd(DCM_DIR{5})
clear all
load workspace_group
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


% BMC/BMA - first level
%--------------------------------------------------------------------------
[rcm,bmc,bma] = spm_dcm_bmr(gcm,{'A','B'});

% BMC/BMA - second level
%==========================================================================

% Empirical Bayes and model comparison (specified model space)
%--------------------------------------------------------------------------
M.hE = -2; REB.Eh;
M.hC = 1/8;
PEB  = spm_dcm_peb(RCM(:,1),M,{'A','B'});
BMA  = spm_dcm_peb_bmc(PEB,gcm(1,:));


% Bayesian family comparison (where model parameters are defined by k)
%--------------------------------------------------------------------------
spm_figure('GetWin','Bayesian family comparison');clf


str = {'sensory','higher','forward','backward'};
K   = k(:,1:4);
nk  = size(K,2);
for i = 1:nk
    pk(i) = sum(bmc.P(K(:,i)));
    PK(i) = sum(BMA.P(K(:,i)));
end

subplot(2,2,1), bar(pk),          hold on
plot([0 (nk + 1)],[.9 .9],'--r'), hold off
text(1:nk,ones(1,nk)/8,str,'color','w','rotation',90)
xlabel('Family'), ylabel('probability')
title('Grand mean','FontSize',16)
set(gca,'YLim',[0 1]), axis square

subplot(2,2,2), bar(PK),          hold on
plot([0 (nk + 1)],[.9 .9],'--r'), hold off
text(1:nk,ones(1,nk)/8,str,'color','w','rotation',90)
xlabel('Family'), ylabel('probability')
title('Emprical Bayes','FontSize',16)
set(gca,'YLim',[0 1]), axis square



%% indices for plotting
%--------------------------------------------------------------------------
iA   = spm_find_pC(GCM{1},{'A'});
iB   = spm_find_pC(GCM{1},{'B'});
jA   = 1:length(iA);
jB   = (1:length(iB)) + jA(end);


% extract and plot results
%==========================================================================
clear Q
for i = 1:Ns
    
    % data - over subjects
    %----------------------------------------------------------------------
    Y(:,i,1) = GCM{i,1}.xY.y{1}*RCM{i,1}.M.U(:,1);
    Y(:,i,2) = GCM{i,1}.xY.y{2}*RCM{i,1}.M.U(:,1);
    
    % Parameter estimates
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Ep);
    Q(:,i,2) = spm_vec(RCM{i,1}.Ep);
    
end


%% first level parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','data and estimates');clf, ALim = 2;

subplot(2,2,1)
plot(Q(iA,:,1),Q(iA,:,2),'.c','MarkerSize',12), hold on
plot(Q(iB,:,1),Q(iB,:,2),'.b','MarkerSize',12), hold on
plot([-1 1],[-1 1],'-.'), hold off
xlabel('standard'), ylabel('empirical Bayes')
title('empirical shrinkage','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

% group means
%--------------------------------------------------------------------------
Qp   = spm_vec(bma.Ep);
r    = corr(Qp(iB),BMA.Ep(jB));
tstr = sprintf('Grand mean vs. Empirical Bayes: cor = %-0.2f',r);

subplot(2,2,2)
plot(Qp(iA),BMA.Ep(jA),'.c','MarkerSize',16), hold on
plot(Qp(iB),BMA.Ep(jB),'.b','MarkerSize',16), hold on
plot([-1 1],[-1 1],'-.'), hold off
xlabel('Grand average estimate'), ylabel('empirical Bayesian estimate')
title(tstr,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

% plot data
%==========================================================================
pst    = gcm{1}.xY.pst;
g(:,1) = gcm{1}.xY.y{1}*gcm{1}.M.U(:,1);
g(:,2) = gcm{1}.xY.y{2}*gcm{1}.M.U(:,1);

subplot(2,2,3)
plot(pst,Y(:,:,1), 'k'), hold on
plot(pst,Y(:,:,2),':k'), hold on
plot(pst,g(:,1),'r','LineWidth',2), hold on
plot(pst,g(:,2),'r','LineWidth',2), hold on
xlabel('pst (ms)'), ylabel('response'), title('Group data','FontSize',16)
axis square, spm_axis tight, a = axis;

subplot(2,2,4)
plot(pst,Y(:,:,2) - Y(:,:,1),'k'), hold on
plot(pst,g(:,2) - g(:,1),'r','LineWidth',2), hold on
xlabel('pst (ms)'), ylabel('differential response'), title('Difference waveforms','FontSize',16)
axis square, spm_axis tight, axis(a)



