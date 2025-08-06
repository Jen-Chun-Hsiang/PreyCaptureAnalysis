% assume your data is in the column vector `x`, with 0 ≤ x ≤ 1
half_gap = mean(diff(edges))/2;
x = [];
num_sample = 10000;
hist_sample= round(hist*num_sample/sum(hist));
for i = 1:numel(hist_sample)
    x = [x; (rand(hist_sample(i), 1)-0.5)*half_gap + edges(i) + half_gap];
end
figure; 
subplot(1, 2, 1);
plot(0.5*(edges(2:end)+edges(1:end-1)), hist);
subplot(1, 2, 2);
histogram(x);
%%

% x is your column vector of data in (0,1)
% 1) Fit Beta with fitdist (no bounds argument)
pd = fitdist(double(x),'Beta');
%    pd is a prob.BetaDistribution object with properties pd.a and pd.b

% 2) Pull out the estimated shape parameters
alphaHat = pd.a;
betaHat  = pd.b;

% 3) Generate N random draws from the fitted Beta
N = 10000;
r = random(pd, N, 1);

% 4) (Optional) Visual check
figure; 
histogram(x,50,'Normalization','pdf','EdgeAlpha',0.2); 
hold on
xx = linspace(0,1,200);
plot(xx, pdf(pd,xx), 'r-', 'LineWidth',2);
xlabel('x'); ylabel('density');
title('Data histogram and fitted Beta PDF');

%% ——————————————————————————————————————————————
% 1) Load or define your data vector x (values in [0,1]).
%    e.g. x = yourData(:);
x = double(x);
% 2) Define the mixture‐PDF handle
mixPdf = @(data, w, a, b, mu, sigma) ...
    w.*betapdf(data,a,b) + (1-w).*normpdf(data,mu,sigma);

% 3) Choose sensible starting values
w0     = 0.5;               % mixing weight
a0     = 2;   b0 = 5;       % Beta shape guesses
mu0    = mean(x);           % Normal mean guess
sigma0 = std(x);            % Normal std‐dev guess

% 4) Run MLE with bounds
phat = mle(x, ...
    'pdf',       mixPdf, ...
    'start',     [w0 a0 b0 mu0 sigma0], ...
    'LowerBound',[0  eps  eps -Inf  eps], ...  % w∈[0,1], a>0, b>0, σ>0
    'UpperBound',[1  Inf  Inf  Inf  Inf]...       % no upper bounds on a,b,μ
    );

wHat     = phat(1);
aHat     = phat(2);
bHat     = phat(3);
muHat    = phat(4);
sigmaHat = phat(5);

% 5) Package into distrib objects (optional)
pdBeta = makedist('Beta',   'a',aHat,'b',bHat);
pdNorm = makedist('Normal', 'mu',muHat,'sigma',sigmaHat);

% 6) Draw new samples from your fitted mixture
N = 10000;
u = rand(N,1) < wHat;         % decide which component
y = zeros(N,1);
y(u)    = random(pdBeta, sum(u),1);
y(~u)   = random(pdNorm, sum(~u),1);

% 7) Visual check
figure; hold on
histogram(x,50,'Normalization','pdf','EdgeColor','none','FaceAlpha',0.3);
xx = linspace(0,1,200);
mixPDF_fitted = wHat*pdf(pdBeta,xx) + (1-wHat)*pdf(pdNorm,xx);
plot(xx,mixPDF_fitted,'r-','LineWidth',2);
legend('Data','Fitted mixture'); xlabel x, ylabel pdf


%%
figure; histogram(r);

save('betaParams_52808.mat','alphaHat','betaHat')

