% Lab2
% Generate Some data
n = 500;
A1 = [1 -0.6]; A2 = [1 .90 .78];
C = 1; B = [0 0 0 0 0.4];
e = sqrt(1.5) * randn(n + 100,1);
w = sqrt(2) * randn(n + 200,1);

A3 = [1 0.5]; C3 = [1 -0.3 0.2];
u = filter(C3, A3,w); u = u(101:end);
y = filter(C,A1,e) + filter(B,A2,u);
u = u(101:end); y = y(101:end);
clear A1 A2 C B e w A3 C3;

%% Look at the acf, pacf and normplot
acfpacf(u);
%% Try some different model orders. 
modelinit = idpoly([1 0], [], [1 0 0]);
model_arma = pem(u, modelinit);
%% Check out the residuals
res = filter(model_arma.A, model_arma.C, u);
acfpacf(res);
% QUESTION 1: Lowest fpe from an ARMA(1,2) model. This is the same model we used to
% generate the input. 

%%
u_pw = res;
y_pw = filter(model_arma.A, model_arma.C, y);
[w temp]= size(u_pw);

M = 40; stem(-M:M, crosscorr(u_pw, y_pw, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;

% Using table 4.7 the impulse weights for the transfer functions seems to
% be something along the lines (4,2,2) as we have a sinosoidal behavior;
% implying that the roots of the A2(Z) are complex valued and r = 2, with a
% decay starting from 4. Assuming d = 4 this gives us s = 0.

%%
d = 4;
r = 2;
s = 0;

A2 = [1 zeros(1,r)];
B = [zeros(1,d) 1 zeros(1,s)];
Mi = idpoly([1], [B], [], [], [A2]);
z_pw = iddata(y_pw, u_pw);
Mba2 = pem(z_pw, Mi); present(Mba2)
vhat = resid(Mba2,z_pw);

%%
a = crosscorr(u_pw, vhat.OutputData, M);
M = 40; stem(-M:M, crosscorr(u_pw, vhat.OutputData, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;

sum(abs(a) > 2/sqrt(n)) / length(a) % Number of observations outside of our conf interval
%%
acfpacf(vhat.OutputData);

% QUESTION 2: delay = 4, B(z); order 0, A(z); order 2. V(t) seems to be uncorrelated
% with u_pw. v(t) does not seem like a white noise process. 

%%
C1 = [1 0];
A1 = [1 0 0];
B = Mba2.B;
A2 = Mba2.F;
Mi_2 = idpoly([1], [B], [C1], [A1], [A2]);
Mba2_2 = pem(z_pw, Mi_2); 
fpe(Mba2_2)
vhat_2 = resid(Mba2_2,z_pw);
acfpacf(vhat_2.OutputData);
% residual looks very white
%%
a_2 = crosscorr(u_pw, vhat_2.OutputData, M);
M = 40; stem(-M:M, crosscorr(u_pw, vhat_2.OutputData, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;
sum(abs(a_2) > 2/sqrt(n)) / length(a_2)
% Still some correlation here. Could try to reduce it further.
% QUESTION 3: ord(A1) = 2, ord(C1) = 1. No all the dependence was not
% removed.

%%
A1 = [1 0 0];
A2 = [1 0 0];
B = [1 0 0 0 0];
C = [1 0];
Mi = idpoly(1,B,C,A1,A2);
z = iddata(y,u);
MboxJ = pem(z,Mi);
present(MboxJ)
ehat = resid(MboxJ,z);

%%
acfpacf(ehat.OutputData)
%%
whitenessTest(ehat.OutputData)
% QUESTION 4. Yes it seems very white if we study acfpacf and the
% whitenesstest.
%% ---------------------- SECTION 3.2 ---------------------------
tork = tork - repmat(mean(tork), length(tork), 1);
y = tork(:,1); u = tork(:,2);
z = iddata(y,u);
n = length(y);
plot(z(1:300))
%%
acfpacf(u)
%% Try some different model orders. (AR(1) is a main suspect)
A = [1 0];
C = [];
modelinit = idpoly(A, [], C);
model_ar1 = pem(u, modelinit);
fpe(model_ar1)

%% Check out the residuals
res = filter(model_ar1.A, model_ar1.C, u);
acfpacf(res);
%Looks good
%%
u_pw = res;
y_pw = filter(model_ar1.A, model_ar1.C, y);
[w temp]= size(u_pw);

M = 40; stem(-M:M, crosscorr(u_pw, y_pw, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;
%%
d = 3;
r = 2;
s = 2;

A2 = [1 zeros(1,r)];
B = [zeros(1,d) 1 zeros(1,s)];
Mi = idpoly([1], [B], [], [], [A2]);
z_pw = iddata(y_pw, u_pw);
Mba2 = pem(z_pw, Mi); present(Mba2)
vhat = resid(Mba2,z_pw);
%%
a = crosscorr(u_pw, vhat.OutputData, M);
M = 40; stem(-M:M, crosscorr(u_pw, vhat.OutputData, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;

sum(abs(a) > 2/sqrt(n)) / length(a) % Number of observations outside of our conf interval

%%
acfpacf(vhat.OutputData)
%%
C1 = [];
A1 = [1 0 0 0];
B = Mba2.B;
A2 = Mba2.F;
Mi_2 = idpoly([1], [B], [C1], [A1], [A2]);
Mba2_2 = pem(z_pw, Mi_2); 
fpe(Mba2_2)
vhat_2 = resid(Mba2_2,z_pw);
acfpacf(vhat_2.OutputData);

%%
a_2 = crosscorr(u_pw, vhat_2.OutputData, M);
M = 40; stem(-M:M, crosscorr(u_pw, vhat_2.OutputData, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;
sum(abs(a_2) > 2/sqrt(n)) / length(a_2)
%%
x = y - filter(B, A2, u);
a_2 = crosscorr(u, x, M);
M = 40; stem(-M:M, crosscorr(u, x, M));
title('Cross Correlation Function'), xlabel('Lag')
hold on;
plot(-M:M, 2/sqrt(n)*ones(1, 2*M+1), '--');
plot(-M:M, -2/sqrt(n)*ones(1, 2*M+1), '--');
hold off;
sum(abs(a_2) > 2/sqrt(n)) / length(a_2)
% Looks uncorrelated
%%
analyzets(x)
% Choose AR(1) -> AR1 = [1 0], C = [] (DID NOT DO IN FIRST PART!)
%%
A1 = [1 0];
A2 = Mba2.F;
B = Mba2.B;
C = [];
Mi = idpoly(1,B,C,A1,A2);
z = iddata(y,u);
MboxJ = pem(z,Mi);
present(MboxJ)
ehat = resid(MboxJ,z);
%%
acfpacf(ehat.OutputData);
%%
whitenessTest(ehat.OutputData);
% Passess all the tests!!
% QUESTION 5: · delay is 3·0.08 = 0.24 seconds
%             · residal is white, All paramters are significant!
%            

%% ------------------- SECTION 3.3 ---------------------------
y = svedala;
A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];
k = 26


%%
[CS,AS] = equalLength(C,A);
[Fk,Gk] = deconv(conv([1, zeros(1,k-1)], CS), AS);
yhat_k = filter(Gk, C, y);

%%
e_noise = var(filter(A,C,y));
% QUESTION 6: the estimated noice variance is 0.3895
mean(yhat_k-y);
% QUESTION 7: 0.0029 and 0
epsilon_variance_theoretical =  sum(Fk.^2)*e_noise;
epsilon_variance_estimated = var(yhat_k - y)
% Question 8: · 4.5795
%             · 4.1435
conf = 1.96*epsilon_variance_theoretical;
% QUESTION 9: conf_interval = [-1.96·V(prediction error),1.96·V(prediction error)]
%           : conf_interval = [-8.9757, 8.9757]
res = yhat_k - y;
sum(abs(res) > epsilon_variance_estimated)/length(y);
% Question 10: 4.92 %
%%
subplot(211)
plot(y)
hold on
plot(yhat_k, 'yellow')
title('The process (blue) and the predictions (yellow)')
subplot(212)
plot(res)
title('residuals')

%%
covf(res,4)
acfpacf(res);
%No the resdiauls does not behave like an MA(k-1) process. 
%% -------------------- SECTION 3.4 --------------------------
A = [1 -1.49 0.57];
B = [0 0 0 0.28 -0.26];
C = [1];
% QUESTION 12: the delay is 3. (indicated by the number of zeros in the
% begining of the B vector)
%%
u = sturup;
% Old estimation
[CS,AS] = equalLength(C,A);
[Fk_y,Gk_y] = deconv(conv([1, zeros(1,k-1)], CS), AS);
yhat_k = filter(Gk_y, C, y);

% New estimation
[CS_u,AS_u] = equalLength(conv(B, Fk_y),C);
[Fk_u,Gk_u] = deconv(conv([1, zeros(1,k-1)], CS_u), AS_u);

uhat_k = filter(Gk_u, C, u);
yhat_k_armax = uhat_k + yhat_k
%%
plot(y)
hold on
plot(yhat_k_armax, 'yellow')
hold on
plot(yhat_k, 'red')
%%
res_armax = yhat_k_armax - y;
mean(res_armax)
var(res_armax)
% The new mean is alot closer to 0 and the varance was reduced by a factor
% of 2. 

%%
acfpacf(res)

%% ----------------------- SECTION 3.5 --------------------------------
acfpacf(y)
% QUESTION 14: The period is 24.
%%
s = 24;
AS = [1 zeros(1, s-1) -1];
A_new = conv(A, AS);
%%
% w = filter([1 zeros(1, s-1) -1], 1, y);
% w = w(24:end);
% %% PEM modelling
% modelinit = idpoly([1 0 0], [],[1 zeros(1,24)]);
% modelinit.Structure.c.Free = [zeros(1,24) 1];
% model_armax = pem(w, modelinit)
% %%
% res_s = filter(model_armax.A, model_armax.C, w);

%%
A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];
k = 26
%%
[CS_s,AS_s] = equalLength(C,conv(A,AS));
[Fk_s,Gk_s] = deconv(conv([1, zeros(1,k-1)], CS_s), AS_s);
yhat_k_s = filter(Gk_s, C, y)
%%
var(yhat_k_s - y)
plot(y)
hold on
plot(yhat_k_armax, 'yellow')
hold on
plot(yhat_k_s, 'red')