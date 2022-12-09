function [p_RS, tDOF_RS, Yhat_RS, DriftTerms_RS, CNR_RS, nHat_RS, Gamma_RS, Noise_RS, PercChange_RS, ConstReg_RS, CapniaTrigSeries_delay] = RegressionAnalysis(DataIn_temp, CapniaTrigSeries, TimeSec, Tau, Tau2, DelayTime, plot_flag, flag_no_SPION_decay, flag_no_hemodynamic_tau, flag_FirstPerTransient)

sx = size(DataIn_temp,2);
sy = size(DataIn_temp,1);
n = size(DataIn_temp,3); % number of time points

TimeSec = TimeSec - TimeSec(1); % ensure time vector starts at zero
%% Set up regressors:
[ActivationRegSelected,CapniaTrigSeries_delay] = BlockRegressor(TimeSec,CapniaTrigSeries,Tau,Tau2,DelayTime,plot_flag,flag_no_SPION_decay,flag_no_hemodynamic_tau);
ActivationRegSelected = ActivationRegSelected/(max(ActivationRegSelected)-min(ActivationRegSelected));

Const = ones(size(ActivationRegSelected));
Linear = TimeSec/TimeSec(end);
Quad = Linear.^2;

X = [ActivationRegSelected(:),Const(:),Linear(:)];
X = [X, Quad(:)]; 

if flag_FirstPerTransient
    len = find(abs(diff(CapniaTrigSeries))>0,1,'first');   %%%%%%%%%%%% shoudl this be CapniaTrigSeries_delay??? 
    Tau3 = len/3;
    FirstPerTrans = exp(-(1:length(ActivationRegSelected))/Tau3);
    X = [X,FirstPerTrans(:)];
end

num_regressors = size(X,2);

%% Solve:
DataIn_Reshape = reshape(DataIn_temp(:,:,1:end),sy*sx,n)';
Y = DataIn_Reshape;
Betas = inv(X'*X)*X'*Y;
Betas_RS = reshape(Betas',sy,sx,num_regressors);

Yhat = X*Betas;
Yhat_RS = reshape(Yhat',sy,sx,n);

%% Drift terms (all regressors except for activation & constant offset)
DriftTerms = X(:,3:num_regressors)*Betas(3:num_regressors,:); 
DriftTerms_RS = reshape(DriftTerms',sy,sx,n);

nHat = Y-Yhat; % nHat is the error
nHat_RS = reshape(nHat',sy,sx,n);
DOF = length(TimeSec)-size(X,2);
ResidualVariance = reshape(diag(nHat'*nHat)/DOF,sy,sx); % Mean Square Error

%% T-stat & p-value
C = [1,zeros(1,num_regressors-1)];
Gamma = C*Betas; % Gamma = Betas for Activation regressor only
Gamma_RS = reshape(Gamma,sy,sx);
ActivationReg_Energy = C*inv(X'*X)*C';
ContrastVariance_RS = ActivationReg_Energy*ResidualVariance; % error^2/Power
ContrastVariance = reshape(ContrastVariance_RS,sy*sx,1)';
tDOF = Gamma./sqrt(ContrastVariance); % coefficients/sqrt(noise^2/Power) --> how far from mean relative to noise
tDOF_RS = reshape(tDOF,sy,sx);
p = 2*(1-tcdf((abs(tDOF)),DOF))*sx*sy; % p is scaled by number of tests done (Bonferroni correction)
p_RS=reshape(p,sy,sx);

%% CNR:
Noise = abs(std(nHat,0,1)); 
Noise_RS = reshape(Noise,sy,sx);

CNR = Gamma./(Noise);
CNR_RS = reshape(CNR,sy,sx);

%% Constant regressor:
C1 = zeros(size(C)); C1(2) = 1;
ConstReg = C1*Betas; % Gamma = Betas for Constant regressor only
ConstReg_RS = reshape(ConstReg,sy,sx);

%% Percent signal change:
PercChange = Gamma./ConstReg*100;
PercChange_RS = reshape(PercChange,sy,sx);

%% Optional plotting 
if plot_flag
    cmap = 'hot';
    figure,
    ylim([min(ActivationRegSelected) max(Const)]*1.1); xlim([min(TimeSec) max(TimeSec)]),
    [legvals,~] = plot_boxes_BlockActivation(CapniaTrigSeries_delay,TimeSec,0.9);
    hold on, plot(TimeSec,X(:,1),'k','LineWidth',2)
    hold on, plot(TimeSec,X(:,2),'Color',[52, 94, 235]/255,'LineWidth',2),
    hold on, plot(TimeSec,X(:,3),'Color',[245, 78, 78]/255,'LineWidth',2),
    hold on, plot(TimeSec,X(:,4),'Color',[148, 68, 201]/255,'LineWidth',2),
    LegStr = {'Activation','Constant','Linear','Quadratic'};
    if flag_FirstPerTransient
        hold on, plot(TimeSec,X(:,5),'Color',[.2, .7, .2],'LineWidth',2),
        LegStr{length(LegStr)+1} = 'Initial Transient'; 
    end
    xlim([TimeSec(1) TimeSec(end)])
    f1=(get(gca,'Children'));
    legend([f1(5:-1:1);f1(end-1:end)],[LegStr,legvals])
    set(gca,'FontSize',14);
    title('Regressors')
    xlabel('Time [sec]')
    
    figure,
    subplot(num_regressors,1,1), imagesc(Betas_RS(:,:,1)); axis image; title('Activation'); colormap(cmap); colorbar;
    subplot(num_regressors,1,2), imagesc(Betas_RS(:,:,2)); axis image; title('Constant'); colormap(cmap); colorbar;
    subplot(num_regressors,1,3) ,imagesc(Betas_RS(:,:,3)); axis image; title('Linear'); colormap(cmap); colorbar;
    subplot(num_regressors,1,4) ,imagesc(Betas_RS(:,:,4)); axis image; title('Quadratic'); colormap(cmap); colorbar;
    if flag_FirstPerTransient
        subplot(num_regressors,1,5) ,imagesc(Betas_RS(:,:,5)); axis image; title('Initial Transient'); colormap(cmap); colorbar;
    end
    sgtitle('Beta for each regressor');
end

end