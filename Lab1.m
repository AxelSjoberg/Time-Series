%-------------------------------- Lab 1 ---------------------------------
%       ....................... 2 .........................
%% 2.4

% Import Data
y_prep = svedala

% Differentation
s = 5;
w_prep = filter([1 zeros(1, s-1) -1], 1, y_prep);
w_prep = w_prep(24:end);

% Parameter estimation
modelinit_prep = idpoly([1 0 0 0], [],[1 zeros(1,12)]);
modelinit_prep.Structure.a.Free = [ 1 1 0 1];
modelinit_prep.Structure.c.Free = [ 1  1 zeros(1,10) 1];
model_armax_prep = pem(w_prep, modelinit_prep)



%        ...................... 3 ........................
%%  ----------------------SECTION 3.1--------------------------------
A1 = [1 -1.79 0.84];
C1 = [1 -0.18 -0.11];

A2 = [1 -1.79];
C2 = [1 -0.18 -0.11];

ARMA_poly1 = idpoly(A1, [], C1);
ARMA_poly2 = idpoly(A2, [], C2);
%%
pzmap(ARMA_poly2)

%%
N = 200;
sigma2 = 1.5;
e1 = sqrt(sigma2)*randn(N, 1);
e2 = sqrt(sigma2)*randn(N, 1);



y1 = filter(ARMA_poly1.c, ARMA_poly1.a, e1);
y2 = filter(ARMA_poly2.c, ARMA_poly2.a, e1);

%%
subplot(211)
plot(y1)
subplot(212)
plot(y2)
% ----Q1----- Converges because the A(z) polynomial has roots oustide of the unit
% circle

%%
m = 200;

rtheo = kovarians(ARMA_poly1.c, ARMA_poly1.a, m);
stem(0:m,rtheo * sigma2);
hold on
rest = covf(y1, m+1)
stem(0:m,rest, 'r')

% -----Q2-----  1. We have limited data. If our N -> inf, the estimated covariance function
%                   coverges to the theoretcial. 
%               2. About 24. The theoretical N/4 = 50. 

%%
subplot(221);
acf(y1, 20,0.05 ,true);
title('acf');
subplot(222);
pacf(y1, 20,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(y1);
% AR(2) seems reasonable, data seems to be Gaussian dist

%%
data = iddata(y1);
na = 3;
nc = 1;
armodel = arx(y1,[na]);
arma_model = armax(y1,[na nc]);
e_hat = filter(arma_model.a, arma_model.c, y1);

% -----------armodel---------
%            Fit          FPE           MSE
% na = 1:   71.85 %      3.226         3.162
% na = 2:   81.29 %      1.454         1.397
% na = 3:   81.52 %      1.447         1.362          
% na = 4:   81.74 %      1.441         1.330
% na = 5:   81.85 %      1.453         1.314
% na = 2 is probably good enough to not overfit


% ----------arma_model--------- (p = 1)
%           Fit           FPE           MSE
% nc = 1:  77.94 %       1.981         1.942    
% nc = 2:  79.46 %       1.735         1.684
% nc = 3:  80.35 %       1.603         1.541


% ----------arma_model--------- (p = 2)
%           Fit           FPE           MSE
% nc = 1:  81.30 %       1.438         1.396    
% nc = 2:  81.36 %       1.442         1.386
% nc = 3:  81.37 %       1.456         1.385

% ----------arma_model--------- (p = 3)
%           Fit           FPE           MSE
% nc = 1:  81.44 %       1.430         1.374    
% nc = 2:  81.45 %       1.444         1.373
% nc = 3:  81.45 %       1.458         1.373

%acf(e_hat, 50, 0.05, true);

% AR(2) makes sense if we look at the resdiual plots.
% The best according to FPE is ARMA(3,1) however I would suggest the AR(2)
% after studying the residual plots. 

%% Resdidual Analysis
subplot(221);
acf(e_hat, 50,0.05 ,true);
title('acf');
subplot(222);
pacf(e_hat, 20,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(e_hat);

%% ------------------------ SECTION 3.2 ----------------------------
n = 500;
A = [1 -1.35 0.43];
sigma2 = 4;
noise = sqrt(sigma2) * randn(n+100,1);
y = filter(1,A, noise);
y = y(101:end);
subplot(211);
plot(y);
subplot(212);
plot(noise);

% Takes time for the process to set. y1 should depend depend on y0 and
% y(-1). These do not exist and we therefor set our y1 and analogously y2
% to 0. 

%%
nest = floor(2/3 * n);
y_est = iddata(y(1:nest));       % Training data
y_val = iddata(y(nest+1:end));   % Validation data

%%
NN = [1:10]';

V = arxstruc(y_est, y_val, NN);
norder = selstruc(V,0);
naic = selstruc(V, 'aic');
%%
n_aic = 0;
n_order = 0;
for i= 1:100
    noise = sqrt(sigma2) * randn(n+100,1);
    y = filter(1,A, noise);
    y = y(101:end);     
    y_est = iddata(y(1:nest)); 
    y_val = iddata(y(nest+1:end));
    V = arxstruc(y_est, y_val, NN);
    n_order(i) = selstruc(V,0);
    n_aic(i) = selstruc(V, 'aic');
end

%%
subplot(211);
histogram(n_order);
title('Least Squares');
subplot(212);
histogram(n_aic);
title('AIC');

% 1. AIC chooses much smaller model order than LS. 
% 2. ?

%%
armodel = arx(y, n_order(end));
armodel.NoiseVariance
armodel.CovarianceMatrix
present(armodel)


%% ------------------ SECTION 3.3 ------------------------
% AR
data_1 = table2array(data1);
data_1=iddata(data_1);

%%
ar1_model = arx(data_1,[1]);
ar2_model = arx(data_1,[2]);
ar3_model = arx(data_1,[3]);
ar4_model = arx(data_1,[4]);
ar5_model = arx(data_1,[5]);

rar1 = resid(ar1_model, data_1);
rar2 = resid(ar2_model, data_1);
rar3 = resid(ar3_model, data_1);
rar4 = resid(ar4_model, data_1);
rar5 = resid(ar5_model, data_1);

% 5.1 the estimation for a1 is fairly good for the AR(1) but poor for the
% AR(2) 
% 5.2 AR(4) is best if we go by FPE, theory?

%%
am11_model = armax(data_1, [1 1]);
am22_model = armax(data_1, [2 2]);


fpe(am11_model); % 1.1818 higher
fpe(am22_model); % 1.1931


%% ------------- SECTION 3.4 ------------------
A   = [1 -1.5 0.7];
C   = [1 zeros(1,11) -0.5];
A12 = [1 zeros(1,11) -1];
A_star = conv(A, A12);
e = randn(600, 1);
y = filter(C, A_star, e);
y = y(100:end);
plot(y);

% Varaince not stationary, transform data
%%
subplot(221);
acf(y, 125,0.05 ,true);
title('acf');
subplot(222);
pacf(y, 125,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(y);

%% Transform data
bcNormPlot(y_1);
%%
[N,temp] = size(y);
for i = 1:N
    if(y(i) >= 0)
        y_stat(i) = sqrt(y(i));
    end
    
    if(y(i) < 0)
        y_stat(i) = -sqrt(-y(i));
    end
 
end
plot(y_stat)
% THis is probably not a good idea

%%
y_s =filter([1 zeros(1,11) -1], 1, y);
y_s = y_s(12:end);
data = iddata(y_s);

%%
subplot(221);
acf(y_s, 50,0.05 ,true);
title('acf');
subplot(222);
pacf(y_s, 50,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(y_s);

% Looks more clean however the seasonality seem to increase (s=12)

%%
modelinit = idpoly(A,[],C)

%%
modelinit = idpoly([1 0 0], [],[]);
model_armax = pem(data, modelinit)
%%
res = filter(model_armax.A, model_armax.C, y_s);
subplot(221);
acf(res, 120,0.05 ,true);
title('acf');
subplot(222);
pacf(res, 120,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(res);

% We have a strong dependecy at s = 12

%%
modelinit = idpoly([1 0 0], [],[1 zeros(1,12)]);
modelinit.Structure.c.Free = [zeros(1,12) 1];
model_armax = pem(data, modelinit)
% Looks good!
% Param est A = [1 -1.53 0.7552] and C = [1 0...0 -04139] although i think
% C = 1 was equally good. Statistically significant? (ASK)

%% ----------------------- SECTION 3.5 ----------------------------

y = svedala;
s = 24;
w = filter([1 zeros(1, s-1) -1], 1, y);
w = w(24:end);
%% PEM modelling
modelinit = idpoly([1 0 0], [],[1 zeros(1,24)]);
modelinit.Structure.c.Free = [zeros(1,24) 1];
model_armax = pem(w, modelinit)

%%
res = filter(model_armax.A, model_armax.C, w);
subplot(221);
acf(res, 30,0.05 ,true);
title('acf');
subplot(222);
pacf(res, 30,0.05 ,true, true);
title('pacf');
subplot(212);
normplot(res);
% Looks about white
% A(z) = 1 - 1.395 z^-1 + 0.4556 z^-2 
% C(z) = 1 - 0.7811 z^-24     

