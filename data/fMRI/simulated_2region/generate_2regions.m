%% Simulates a simple 2 region DCM with a block design
%
% There are 2 groups of subjects. The self-connection (A) of each model is 
% is sampled from a distribution for the relevant group. Each DCM has 
% 2 inputs: stimulus (driving input) and modulator (modulatory input).
%
% The 1st level model space is as follows:
%
% 1. Modulation of fwd and bwd connections
% 2. Fwd only
% 3. Bwd only
% 4. No modulation
%
% Simulated data are created from Model 2 for two groups of subjects, where
% each group differs in the strength of the modulatory connection from 
% R1->R2. The second level model space includes three covariates:
%
% 1. Mean
% 2. Group membership (half -1s, half 1s)
% 3. Random
% _________________________________________________________________________

% Settings

trial_duration  = 2;   % secs
isi             = 4;   % secs
TR              = 1.5; % secs
num_trials      = 50;
num_null_trials = 20;   
num_extra_scans = 2;

t  = 16;  % Microtime resolution
t0 = 8;   % Microtime onset

SNR_std = 2;  % First level SNR based on the standard deviation

a_mean         = 0.2;       % Exogenous connection strength (Hz)
group_means    = [0.3 0.7]; % Modulatory input stength (Hz)
beta           = 16;        % Within:between subjects variance
group_size     = 15;        % Subjects per group

ns = length(group_means) * group_size; % Number of subjects

%% Build stimulus timeseries
rng(1);

% Onsets and durations (secs)
trial_and_isi   = trial_duration + isi;
onsets_secs     = 0:trial_and_isi:((trial_and_isi*num_trials)-1);
total_scan_time = trial_and_isi * num_trials + (num_extra_scans*TR);
nscans          = ceil(total_scan_time / TR);

% Set some trials as null
null_trials = randperm(num_trials);
null_trials = null_trials(1:num_null_trials);
onsets_secs(null_trials) = [];

% Set half of trials as modulated
ntrials    = length(onsets_secs);
idx        = randperm(ntrials);
idx        = idx(1:ceil(ntrials/2));
mod_onsets = onsets_secs(idx);

% Build dummy SPM
SPM = struct();
SPM.nscan     = nscans;
SPM.xBF.T     = t;
SPM.xBF.dt    = TR / t;
SPM.xBF.units = 'secs';
SPM.xY.RT     = TR;
SPM.Sess.U(1).name = {'Stimulus'};
SPM.Sess.U(1).ons  = onsets_secs';
SPM.Sess.U(1).dur  = repmat(trial_duration, ntrials, 1);
SPM.Sess.U(1).P    = struct('name','none','h',0,'i',1);
SPM.Sess.U(2).name = {'Modulator'};
SPM.Sess.U(2).ons  = mod_onsets';
SPM.Sess.U(2).dur  = repmat(trial_duration, length(mod_onsets), 1);
SPM.Sess.U(2).P    = struct('name','none','h',0,'i',1);

% Secs -> microtime
SPM.Sess.U = spm_get_ons(SPM,1);

%% Build DCM output structure
Y      = struct();
Y.name = {'R1','R2'};
Y.dt   = TR;
X0     = ones(nscans,1);
Y.X0   = X0;
%% Specify DCMs
rng(1);

n  = 2; % Number of regions
nm = 4; % Number of models
nu = 2; % Number of inputs

% b-matrices (one element per model)
b = cell(1,nm);
b{1} = [0 1
        1 0];
b{2} = [0 0
        1 0];
b{3} = [0 1
        0 0];
b{4} = [0 0
        0 0];
    
% Generative B-matrix
b_gen = zeros(n,n,nu);
b_gen(:,:,2) = [0.00 0.00;
                0.50 0.00];

% Build DCM structure
DCM0   = struct();
DCM0.a = ones(n,n); 
DCM0.b = zeros(1,1,nu); 
DCM0.c = [1 0;
          0 0]; 
DCM0.d = zeros(n,n,0);
DCM0.Y = Y;
DCM0.n = n;
DCM0.v = nscans;
DCM0.delays  = repmat(TR * (t0/t), 1, DCM0.n);
DCM0.options = struct();

% Prepare inputs (convert SPM -> DCM .U format)
DCM0 = spm_dcm_U(DCM0,SPM,1,{1 2});

Ep_template = struct('A',[0 1; 1 0], ...
                     'B',b_gen, ...
                     'C',[1 0; 0 0], ...
                     'D',zeros(n,n,0), ...
                     'transit',zeros(n,1), 'decay',zeros(n,1), 'epsilon',0);

%% Create models
rng(1);

DCMs    = cell(length(group_means)*group_size, nm);
subject = 1;
gen_model = 2;

for group = 1:length(group_means)    
    fprintf('Group %d\n',group);
    
    for s = 1:group_size
        fprintf('Subject %d\n',s);
        
        DCM = DCM0;
        
        for model = 1:nm    
            DCM.b = zeros(n,n,nu); 
            DCM.b(:,:,2) = b{model};
            
            if model == gen_model
                [pE,pC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,DCM.options);
                pC      = spm_unvec(diag(pC), pE);
                
                % Update priors on between-region A-matrix connections (temporarily)   
                idx = ~eye(DCM.n);
                pE.A(idx) = a_mean;
                pC.A      = pC.A ./ beta;
                
                % Update priors on between-group connection (temporarily)
                pE.B(2,1,2) = group_means(group);
                pC.B(2,1,2) = pC.B(2,1,2) / beta;
                                
                % Fix C-matrix
                pE.C(1) = 1;
                pC.C(1) = 0;
                
                % Draw from priors
                pE = full(spm_vec(pE));
                pC = full(spm_vec(pC));                                
                Ep = sqrt(pC) .* randn(length(pE),1) + pE;

                DCM.Ep = spm_unvec(Ep,Ep_template);
            end
            
            % High SNR matching DCM prior
            DCM.Ce = repmat(1/exp(6),DCM.n,1);
            
            DCMs{subject,model} = DCM;                        
        end
        
        subject = subject + 1;
    end    
end
%% Generate simulated timeseries (from model 2) under little noise
rng(1);
idx = 2;
[DCMs,gen] = spm_dcm_simulate(DCMs, [], [], idx);
%% Estimate each subject's generative model to get more stable parameters
DCMs = spm_dcm_fit(DCMs);
%% Generate simulated timeseries under realistic noise
rng(1);
idx = 2;
[DCMs,gen] = spm_dcm_simulate(DCMs, 'SNR_std', SNR_std, idx);

GCM = DCMs;
save('models/GCM_simulated.mat','GCM','gen');
%% Estimate under realistic noise

% Clear out old priors
for i = 1:numel(GCM)
    GCM{i} = rmfield(GCM{i},'M');
end

M = struct();
M.X  = ones(ns,1);                  % preclude any bias during group inversion
M.Q  = 'single';                    % assume all parameters share a variance ratio
M.hE = 0;                           % prior expectation of between subject log precision
M.hC = 1/16;                        % prior variance of between subject log precision
 
[GCM,PEB,M,HCM] = spm_dcm_peb_fit(GCM,M);
save workspace;

% Save group model
save('models/GCM_simulated.mat','GCM','gen','M','HCM');

% Save individual DCMs
for s = 1:ns
    for m = 1:nm
        DCM      = GCM{s,m};
        filename = sprintf('DCM_s%d_m%d.mat', s, m);        
        save(fullfile('models',filename),'DCM');
    end
end
%% Create second level design matrix
rng(1);

% Design matrix (constant, group difference, random)
X = ones(ns,3);
X(1:group_size, 2) = -1;
%X(1:group_size, 2) = 0;
X(:,3) = rand(ns,1);
X(:,3) = X(:,3) - mean(X(:,3));

M = struct();
M.X  = X;

PEB = spm_dcm_peb(GCM,M);

save('PEB_test.mat','PEB');
%% Plot estimated parameters against generative ones
clear gEp Ep;
for s = 1:ns
    gEp(s,:) = full(spm_vec(gen{s}.Ep))';
    Ep(s,:)  = full(spm_vec(GCM{s,2}.Ep))';
end
gEp = mean(gEp);
Ep  = mean(Ep);
figure;scatter(gEp,Ep);
lsline;
[R,P]=corrcoef(gEp,Ep);
title(sprintf('Correlation: %2.2f',R(1,2)));

%% Plot timeseries for one subject from each group
spm_figure('GetWin','Example timeseries'); clf
rows = 2; cols = 2;

for group = 1:2

    % Subplot indices
    idx = group:2:4;
    
    % Get example generative model
    if group == 1
        DCM_gen = gen{1};
        DCM = GCM{1,1};
    else
        DCM_gen = gen{end};
        DCM = GCM{end,1};
    end       
    
    % Neuronal response
    subplot(rows,cols,idx(1));
        
    offset = 0.3;
    
    % Task
    t = DCM.U.dt*[0:(length(DCM.U.u)-1)];
    data = DCM.U.u(:,1) .* (offset+0.05);
    data(data == 0) = offset;
    b1 = bar(t,data,'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.75 0.75 0.75]);      
    set(get(b1,'BaseLine'),'BaseValue',offset);
    hold on;
    
    % Task + modulation
    data = DCM.U.u(:,2) .* (offset+0.05);
    data(data == 0) = offset;
    b2 = bar(t,data,'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);        
    set(get(b2,'BaseLine'),'BaseValue',offset);
    
    line([0 310],[offset offset],'Color','k');
    
    % Neuronal response
    t = DCM.Y.dt*[0:DCM.v-1];
    h = plot(t,DCM_gen.x,'LineWidth',2);
    xlabel('Time (s)');
    ylabel('Connection (Hz)');
    title('Neuronal Response','FontSize',16);  
    xlim([0 310]); ylim([0 0.35]);
    legend([b1 b2 h(1) h(2)],{'Task','Task+Modulation','Region 1','Region 2'});

    % Haemodynamic response
    subplot(rows,cols,idx(2));
    t = DCM.Y.dt*[0:DCM.v-1];
    plot(t,DCM.Y.y);
    xlabel('Time (s)');
    ylabel('BOLD response');
    title('Haemodynamics','FontSize',16);  
    xlim([0 310]);% ylim([-1 4]);
    hold on;
    legend({'Region 1','Region 2'});
end
%% Plot the two groups' evidence and parameters

load('PEB_test.mat');
PEB = PEB(2);

load('models/GCM_simulated.mat');

% Design matrix (group difference)
ns = size(GCM,1);
X  = ones(ns,1);
X(1:group_size, 1) = -1;

% Generative and estimated parameters (B only)
clear Ep gEp;
for s = 1:ns
    gEp(s,:) = full(spm_vec(gen{s}.Ep))';
    Ep(s,:)  = full(spm_vec(GCM{s,2}.Ep))';
end
mgEp1 = mean(gEp(X==-1,:)); vgEp1 = var(gEp(X==-1,:));
mgEp2 = mean(gEp(X==1,:)); vgEp2 = var(gEp(X==1,:));

is_b = spm_find_pC(GCM{s,1},'B');
Ep   = Ep(:,is_b);
gEp   = gEp(:,is_b);

% Estimated parameters BPA
eBPA1 = spm_dcm_bpa(GCM(X == -1,2),true);
eBPA2 = spm_dcm_bpa(GCM(X == 1,2),true);

% Generative parameters BPA
gBPA1 = spm_dcm_bpa(gen(X == -1),true);
gBPA2 = spm_dcm_bpa(gen(X == 1),true);

% PEB parameters
pEp = full(PEB.Ep);
pCp = diag(PEB.Cp);
pCp = reshape(pCp,size(pEp));

% Report B parameters
fprintf('Estimated BPA Group 1: N(%2.2f,%2.2f), BPA Group 2: N(%2.2f,%2.2f)\n', ...
            eBPA1.Ep.B(2,1,2), eBPA1.Vp.B(2,1,2),...
            eBPA2.Ep.B(2,1,2), eBPA2.Vp.B(2,1,2));
fprintf('PEB group mean: N(%2.2f,%2.2f), PEB group difference: N(%2.2f,%2.2f)\n', ...
            pEp(5,1), pCp(5,1),...
            pEp(5,2), pCp(5,2));

rows = 5;
cols = 2;

spm_figure('GetWin','Group comparison');
spm_clf;

subplot(rows,cols,1);
post = spm_dcm_bmc(GCM(X == -1,:));
bar(post);ylabel('P(model)');xlabel('Model');
title('Group -1');

subplot(rows,cols,2);
post = spm_dcm_bmc(GCM(X == 1,:));
bar(post);ylabel('P(model)');xlabel('Model');
title('Group 1');

subplot(rows,cols,3);
spm_plot_ci(mgEp1',vgEp1); ylim([-0.2 1.2]);
hold on;
plot(is_b,gEp(X == -1,:),'.','Color','b');
xlabel('Parameter'); title('Generative model mean');

subplot(rows,cols,4);
spm_plot_ci(mgEp2',vgEp2); ylim([-0.2 1.2]);
hold on;
plot(is_b,gEp(X == 1,:),'.','Color','b');
xlabel('Parameter');  title('Generative model mean');

subplot(rows,cols,5);
spm_plot_ci(gBPA1.Ep,diag(gBPA1.Cp)); ylim([-0.2 1.2]);
hold on;
plot(is_b,gEp(X == -1,:),'.','Color','b');
xlabel('Parameter'); title('Generative model BPA');

subplot(rows,cols,6);
spm_plot_ci(gBPA2.Ep,diag(gBPA2.Cp)); ylim([-0.2 1.2]);
hold on;
plot(is_b,gEp(X == 1,:),'.','Color','b');
xlabel('Parameter');  title('Generative model BPA');

subplot(rows,cols,7);
spm_plot_ci(eBPA1.Ep,diag(eBPA1.Cp)); ylim([-0.2 1.2]);
hold on;
plot(is_b,Ep(X == -1,:),'.','Color','b');
xlabel('Parameter'); title('Estimated Model 2 BPA');

subplot(rows,cols,8);
spm_plot_ci(eBPA2.Ep,diag(eBPA2.Cp)); ylim([-0.2 1.2]);
hold on;
plot(is_b,Ep(X == 1,:),'.','Color','b');
xlabel('Parameter');  title('Estimated Model 2 BPA');

subplot(rows,cols,9);
spm_plot_ci(pEp(:,1),pCp(:,1));
set(gca,'XTickLabel',PEB.Pnames);
title('PEB common');

subplot(rows,cols,10);
spm_plot_ci(pEp(:,2),pCp(:,2));
set(gca,'XTickLabel',PEB.Pnames);
title('PEB difference');
%% Generate dummy SPMs & VOI files
load('models/GCM_simulated.mat');
mkdir('GLM');
for i = 1:ns
    subject_dir = fullfile('GLM',sprintf('S%03d', i));

    mkdir(subject_dir);
    
    % Create dummy SPM
    save(fullfile(subject_dir,'SPM.mat'),'SPM');
    
    % Create region 1 VOI
    xY   = GCM{i,1}.xY(1);
    Y    = GCM{i,1}.Y.y(:,1);
    xY.u = Y;
    save(fullfile(subject_dir,'VOI_R1.mat'), 'Y', 'xY');
    
    % Create region 2 VOI
    xY   = GCM{i,1}.xY(2);
    Y    = GCM{i,1}.Y.y(:,2);
    xY.u = Y;
    save(fullfile(subject_dir,'VOI_R2.mat'), 'Y', 'xY');    
end

%% Generate BMA

load('PEB_test.mat');

% Automatic search
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1));
save('BMA_search.mat','BMA');

% Specific models
GCM_templates = load(fullfile('models','GCM_simulated.mat'));
GCM_templates = GCM_templates.GCM;
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),GCM_templates(1,:));
save('BMA_specific_models.mat','BMA');
