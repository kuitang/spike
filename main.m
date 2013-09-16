%% Load data

PREDICTION_VALUE_COL = 2;
LOWER_95_COL         = 3;
UPPER_95_COL         = 4;
LOWER_68_COL         = 5;
UPPER_68_COL         = 6;
VALUE_COL            = 9;

% Subset the data?
subIdxs = 3000:4000;

%load('model2.mat');

% model.HMMmodel.params are just Gamma(a,b) parameters, and c, d are
% parameters for a Beta prior. These are pretty flat priors, so you
% probably don't need to change anything.

% The code ALWAYS uses DP mixture emissions!

impt     = importdata('345_Inner_Join.csv');
resids   = impt.data(:,VALUE_COL)    - impt.data(:,PREDICTION_VALUE_COL);
conf68Size = impt.data(:,UPPER_68_COL) - impt.data(:,LOWER_68_COL);
conf95Size = impt.data(:,UPPER_95_COL) - impt.data(:,LOWER_95_COL);

sigma    = conf68Size / 2;
% Note: Somehow the ratio is 2.0061, but
% norminv(0.95 + 0.025) / norminv(0.68 + 0.16) = 1.9709
figure;

scaled   = (resids - mean(resids)) ./ sigma;

% Rest of code wants a ROW vector.
scaledSmall = scaled(subIdxs)';

%% Hypers
% Observation dimension
d = 1;

obsModelType = 'Gaussian';
priorType = 'IW-N';  % prior on Gaussian N(mu_{k,j},Sigma_{k,j}) emissions (non-conjugate)

% Uninformative prior
sig0 = 1e6;  % covariance of the N(mu0,sig0) prior on mu_{k,j}
meanSigma = eye(d);  % expected mean of IW(nu,nu_delta) prior on Sigma_{k,j}
Kz = 20;  % truncation level of the DP prior on HMM transition distributions pi_k
Ks = 1;  % truncation level of the DPMM on emission distributions pi_s

% Always using DP mixtures emissions, with single Gaussian forced by
% Ks=1...Need to fix.
model.obsModel.mixtureType = 'infinite';

% Sticky HDP-HMM parameter settings:
model.HMMmodel.params.a_alpha=1;  % affects \pi_z
model.HMMmodel.params.b_alpha=0.01;
model.HMMmodel.params.a_gamma=1;  % global expected # of HMM states (affects \beta)
model.HMMmodel.params.b_gamma=0.01;
if Ks>1
    % Set per-state DP concentrations.
    model.HMMmodel.params.a_sigma = 1;
    model.HMMmodel.params.b_sigma = 0.01;
end
% if isfield(settings,'Kr')
%     if settings.Kr > 1
%         model.HMMmodel.params.a_eta = 1;
%         model.HMMmodel.params.b_eta = 0.01;
%     end
% end
model.HMMmodel.params.c=100;  % self trans
model.HMMmodel.params.d=1;
model.HMMmodel.type = 'HDP';

% We could import other priors but this works well enough.

%%
% Setting for inference:

saveDir = 'tmp';

switch obsModelType
    case {'AR','SLDS'}
        settings.Kz = Kz;   % truncation level for mode transition distributions
        settings.Ks = 1;  % truncation level for mode transition distributions
        if strcmp(obsModelType,'SLDS') && ~strcmp(y_priorType,'IW')
            settings.Kr = Kr;  % truncation level for MoG measurement noise
        end
    case 'Gaussian'
        settings.Kz = Kz;   % truncation level for mode transition distributions
        settings.Ks = Ks;  % truncation level for mode transition distributions
end
settings.Niter = 2000;  % Number of iterations of the Gibbs sampler
settings.resample_kappa = 1;  % Whether or not to use sticky model
settings.seqSampleEvery = 100; % How often to run sequential z sampling
settings.saveEvery = 1000;  % How often to save Gibbs sample stats
settings.storeEvery = 1;
settings.storeStateSeqEvery = 100;
settings.ploton = 1;  % Whether or not to plot the mode sequence while running sampler
settings.plotEvery = 20;
settings.plotpause = 0;  % Length of time to pause on the plot
settings.saveDir = saveDir;  % Directory to which to save files

settings.trial = 1;

%%
% Set Hyperparameters

% Type of dynamical system:
model.obsModel.type = obsModelType;

if strcmp(obsModelType,'AR')
    % Order of AR process:
    model.obsModel.r = r;
    m = d*r;
else
    m = d;
end

% Type of prior on dynamic parameters. Choices include matrix normal
% inverse Wishart on (A,Sigma) and normal on mu ('MNIW-N'), matrix normal
% inverse Wishart on (A,Sigma) with mean forced to 0 ('MNIW'), normal on A,
% inverse Wishart on Sigma, and normal on mu ('N-IW-N'), and fixed A,
% inverse Wishart on Sigma, and normal on mu ('Afixed-IW-N').  NOTE: right
% now, the 'N-IW-N' option is only coded for shared A!!!
model.obsModel.priorType = priorType;

switch model.obsModel.priorType
    case 'NIW'

        model.obsModel.params.M  = zeros([d 1]);
        model.obsModel.params.K =  kappa;
        
    case 'IW-N'
        % Mean and covariance for Gaussian prior on mean:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
    
    case 'MNIW'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);

        % Inverse covariance along rows of A (sampled Sigma acts as
        % covariance along columns):
        model.obsModel.params.K =  K(1:m,1:m);
        
    case 'MNIW-N'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);

        % Inverse covariance along rows of A (sampled Sigma acts as
        % covariance along columns):
        model.obsModel.params.K =  K(1:m,1:m);

        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));

    case 'N-IW-N'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);
        model.obsModel.params.Lambda0_A = inv(kron(inv(K),meanSigma));

        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
        
    case 'Afixed-IW-N'
        % Set fixed A matrix:
        model.obsModel.params.A = A_shared;
        
        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
        
    case 'ARD'
        % Gamma hyperprior parameters for prior on precision parameter:
        model.obsModel.params.a_ARD = 10;
        model.obsModel.params.b_ARD = 0.01;
        
        % Placeholder for initializeStructs. Can I get rid of this?
        model.obsModel.params.M  = zeros([d m]);

        % Mean and covariance for mean of process noise:
        model.obsModel.params.zeroMean = 1;
end
        
% Degrees of freedom and scale matrix for covariance of process noise:
model.obsModel.params.nu = 1000; %d + 2;
model.obsModel.params.nu_delta = (model.obsModel.params.nu-d-1)*meanSigma;

if strcmp(obsModelType,'SLDS')
    % Degrees of freedom and scale matrix for covariance of measurement noise:
    model.obsModel.y_params.nu = 1000; %dy + 2;
    model.obsModel.y_params.nu_delta = (model.obsModel.y_params.nu-dy-1)*y_var*eye(dy);
    
    model.obsModel.y_priorType = y_priorType;
    
    switch model.obsModel.y_priorType
        case 'NIW'
            
            model.obsModel.y_params.M  = zeros([dy 1]);
            model.obsModel.y_params.K =  kappa_y;
            
        case 'IW-N'
            % Mean and covariance for Gaussian prior on mean:
            model.obsModel.y_params.mu0 = zeros(dy,1);
            model.obsModel.y_params.cholSigma0 = chol(sig0_y*eye(dy));
    end
    
    % Fixed measurement matrix:
    model.obsModel.params.C = [eye(dy) zeros(dy,d-dy)];
    
    % Initial state covariance:
    model.obsModel.params.P0 = P0*eye(d);
end

%% Inference
data_struct.obs = scaledSmall;
HDPHMMDPinference(data_struct, model, settings);


%% Generate data from the prior:
% data_struct = generateData(model,settings,T);
% while length(unique(data_struct.true_labels))==1
%     data_struct = generateData(model,settings,T);
% end
% 
% close all;
% figure; plot(data_struct.obs'); hold on; plot(data_struct.true_labels,'r','LineWidth',2); hold off;
% title('True State and Mode Sequences')


%%
% for seq=1
%     data_struct(1).test_cases = seq;
%     for t=trial_vec
% 
%         settings.trial = t;  % Defines trial number, which is part of the filename used when saving stats
% 
%         HDPHMMDPinference(data_struct,model,settings)
%     end
% end

%% Run
%data_struct.obs = residsSmall;
