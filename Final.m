% clear
% clc
%%
plotOn = 0;
%%

% 
% % %%
% load climate66.dat
% load climate67.dat
% load climate68.dat
% load climate69.dat
% load climate70.dat
% load climate71.dat
% load climate72.dat
% load climate73.dat
%%
%mdata = [climate66; climate67; climate68; climate69; climate70];
%vdata = [climate71; climate72]; 
%tdata =  climate73;
% 
%  data = [climate66; climate67; climate68; climate69; climate70; climate71; climate72; climate73];
%  data = [datenum([data(:,1:4) zeros(length(data),2)]) data];

% clear mdata vdata tdata climate66 climate67 climate68 climate69 climate70 climate71 climate72 climate73

% ------ 

start1 = 12097;
start2 = 27505;

tm  = (start1:start1+24*7*10-1)';
tv  = (tm(end)+1:tm(end)+24*7*2)';
tt1 = (tv(end)+1:tv(end)+24*7)';
tp2 = (start2-7*24*2+1:start2-1)'; % Preparation 
tt2 = (start2:start2+7*24-1)';

% % ------
% 
% clearvars -except data tm tv tt1 tt2 tp2 plotOn

%% Defining data sequences

% y

ym  = data(tm, 9);
yv  = data(tv, 9);
yt1 = data(tt1,9);
yp2 = data(tp2,9);
yt2 = data(tt2,9);

% x

xm  = data(tm, 7);
xv  = data(tv, 7);
xt1 = data(tt1,7);
xp2 = data(tp2,7);
xt2 = data(tt2,7);

if plotOn
    figure
    sgtitle("Net radiation")
    subplot(221)
    plot(data([tm; tv; tt1], 1), [xm; xv; xt1])
    ylabel("Net radiation [W/m^2]");
    datetick("x", "dd mmm")
    title("Data before transformation")
    axis tight
    
    subplot(222)
    bcNormPlot([xm; xv; xt1])
    title("BCN before transformation")
end

% Removing y mean

yMean1  = mean([ym; yv; yt1]);
yMean2  = mean([yp2; yt2]);

ym  = ym  - yMean1;
yv  = yv  - yMean1;
yt1 = yt1 - yMean1;
yp2 = yp2 - yMean2;
yt2 = yt2 - yMean2;


% Logging x

xMin1  = min([xm; xv; xt1]);
xMin2  = min([xp2; xt2]);

c = 11;

xm  = log(xm -xMin1+c);
xv  = log(xv -xMin1+c);
xt1 = log(xt1-xMin1+c);
xp2 = log(xp2-xMin2+c);
xt2 = log(xt2-xMin2+c);

% c = 10;
% 
% pwr = 1/5;
% 
% xm  = (xm -xMin1+c).^pwr;
% xv  = (xv -xMin1+c).^pwr;
% xt1 = (xt1-xMin1+c).^pwr;
% xp2 = (xp2-xMin2+c).^pwr;
% xt2 = (xt2-xMin2+c).^pwr;

clear c xMin1 xMin2

% Removing x mean

xMean1 = mean([xm; xv; xt1]);
xMean2 = mean([xp2; xt2]);

xm  = xm  - xMean1;
xv  = xv  - xMean1;
xt1 = xt1 - xMean1;
xp2 = xp2 - xMean2;
xt2 = xt2 - xMean2;

if plotOn
    subplot(223)
    plot(data([tm; tv; tt1], 1), [xm; xv; xt1])
    ylabel("Net radiation [ln(W/m^2)]");
    datetick("x", "dd mmm")
    title("Data after transformation")
    axis tight
    
    subplot(224)
    bcNormPlot([xm; xv; xt1])
    title("BCN after transformation")
end

clear xMean1 xMean2

%% Analyze xm to model it

%analyzets(xm)
% There is a clear periodicity of 24
xdsmodel = armamodel(24, 0, xm, [zeros(1,24) 1]);
%xdsmodel.A = [1 zeros(1,23) -0.8];
xm_ds = filter(xdsmodel.A,1,xm);
analyzets(xm_ds)

%clear xdsmodel 

%%
[xmodel, xm_pw] = armamodel(12, 24, xm_ds, 1, [zeros(1,24) 1], 1);
%[xmodel, xm_pw] = armamodel(24, 0, xm_ds, 1, 1, 1);
%%
xmodel = removeUnsignificant(xmodel, xm_ds, 1, 1);

%%

A3 = xmodel.A;
A3 = conv(xdsmodel.A,A3);
C3 = xmodel.C;

et = filter(A3, C3, ym);
wt = filter(A3, C3, xm);
if length(xm_pw) < length(wt)
    et = et(25:end);
    wt = wt(25:end);
end

%%
crosscorrelation(wt, et);

%% (iii)
d = 0;
r = 3;
s = 0;

B   = [zeros(1,d) 1 zeros(1,r)];
A2  = [1 zeros(1,r)];


Zyx = iddata(ym, xm);

Mi1 = idpoly(1, B, [], [], A2);
M1  = pem(Zyx, Mi1);

B  = M1.B;
A2 = M1.F;

Hx = filter(B,A2,xm);
e_tilde = ym - Hx;

if plotOn
    crosscorrelation(wt, e_tilde);
end

%%
Hw = filter(B, A2, wt);
vt  = et - Hw;

%if plotOn
    crosscorrelation(wt, vt);
%end

%%
if plotOn
    analyzets(e_tilde);
end

%%
[emodel, eres] = armamodel(6, 2, e_tilde, 1, 1,1);
% (6, 2) works

emodel = removeUnsignificant(emodel, e_tilde, 1, 1, 1);

%% (iv)

A1 = emodel.A;
C1 = emodel.C;

z = iddata(ym, xm);

Mi2 = idpoly(1, B, C1, A1, A2);
M2  = pem(z, Mi2);

A  = conv(M2.D, M2.F);
B2 = conv(M2.D, M2.B);
C  = conv(M2.F, M2.C);
k = 7;

% y part estimation
[Fy, Gy] = polydiv(C, A, k);
yhat_y = filter(Gy, C, [ym; yv]);

% x part estimation
yhat_x = filter(conv(B2,Fy),C, [xm;xv]);

yhat = yhat_y + yhat_x; 
yhat = yhat(length(ym)+1:end);

%%

hold on
plot(yv)
plot(yhat)
hold off

%%
var(yv-yhat)

%%

yhat_t1 = filter(Gy, C, [yv; yt1]) + filter(conv(B2,Fy),C, [xv;xt1]);
yhat_t1 = yhat_t1(length(yv)+1:end);
yhat_t2 = filter(Gy, C, [yp2; yt2]) + filter(conv(B2,Fy),C, [xp2;xt2]);
yhat_t2 = yhat_t2(length(tp2)+1:end);

%%
var(yt1-yhat_t1)
var(yt2-yhat_t2)

if plotOn
    figure
    sgtitle("Test data prediction, ARMAX")
    subplot(211)
    hold on
    plot(data(tt1,1),yt1+yMean1, "DisplayName", "Actual data")
    plot(data(tt1,1),yhat_t1+yMean1, "DisplayName", "7 step prediction")
    hold off
    ylabel("Temperature [°C]")
    axis([min(tt1) max(tt1) min(yt1+yMean1-1) (max(yt1+yMean1)+4)])
    datetick("x", "dd mmm")
    legend

    subplot(212)
    hold on
    plot(data(tt2,1),yt2+yMean2, "DisplayName", "Actual data")
    plot(data(tt2,1),yhat_t2+yMean2, "DisplayName", "7 step prediction")
    hold off
    ylabel("Temperature [°C]")
    axis([min(tt2) max(tt2) min(yt2+yMean2-3) (max(yt2+yMean2)+3)])
    datetick("x", "dd mmm")
    legend
end
    
%%

%A  = [1,-2.949697142913895,3.279110351802082,-1.714304877498572,0.382523391934560,0.125791796958986,-0.183380636410895,0.060671206920013];
%B2 = [1.122765295440385,-2.526948725423966,1.915208293984102,-0.585936041464229,0.019885514746477,0.155135629244599,-0.097445883341481];
%C  = [1,-1.603748228296266,0.632429255997286];

p = length(A) -1;
r = length(B)-1 ;
q = length(C) -1;

AK = eye(p+r+q+1);

sigmaA = 0.3;
sigmaB = 4;
sigmaC = 0.7;

%Re = eye(length(AK))*1e0; % Hidden state noise variance matrix
Re = diag([ones(1,p)*sigmaA ones(1,r+1)*sigmaB ones(1,q)*sigmaC])*1e-1;
Rw = 1e3;                  % Guess for observation variance

% Initial values
Rxx0 = eye(length(AK))*1e-4; % Initial variance
a = armax(z,[9 [r+1] 5 [0]])
%x0 = [a.A(2:end) a.B a.C(2:end)]'; % Initial state
x0 = [-1.35307161509510,0.394182697457719,-6.00532472712315e-07,-1.25048001034894e-06,-1.18867738276424e-06,-9.71372847232800e-07,-1.44915291466681e-06,-1.32665764303290e-06,-7.52720493180325e-07,0.947964908869584,-0.373128691429899,-0.349731364263012,5.69619993212059e-07,-0.344283301814754,1.71854718639431e-08,1.51588810044718e-07,-5.69778836018427e-07,6.05673117255761e-08]'
k = 7; 


[y_hat_kalman, xsave] = superKalman2([ym; yv],[xm; xv],AK,Re,Rw,Rxx0,x0,p,r,q,k);

subplot(311)
hold on
plot(yv)
plot(y_hat_kalman(length(ym)+1:end))
subplot(312)
axis tight
hold on
plot([ym; yv])
plot(y_hat_kalman)
subplot(313)
axis tight
plot(xsave)
axis tight

display("Variance of " + k + " steps: " + var(yv-y_hat_kalman(length(ym)+1:end)))




plot(xsave)

%%
x0 = mean(xsave)';
for i=1:length(x0)
    if abs(x0(i)) < 0.2
        x0(i) = 0;
        Re(i,i) = 0;
    end
end

%%
[y_hat_kalman, xsave] = superKalman2([ym; yv],[xm; xv],AK,Re,Rw,Rxx0,x0,p,r,q,26);

subplot(311)
hold on
plot(yv)
plot(y_hat_kalman(length(ym)+1:end))
subplot(312)
axis tight
hold on
plot([ym; yv])
plot(y_hat_kalman)
subplot(313)
axis tight
plot(xsave)
axis tight

plot(xsave)
display("Variance of " + k + " steps: " + var(yv-y_hat_kalman(length(ym)+1:end)))
yhat_k26 = y_hat_kalman(length(ym)+1:end);
%%

[y_hat_kalman, xsave] = superKalman2([yv; yt1],[xv; xt1],AK,Re,Rw,Rxx0,x0,p,r,q,k);

subplot(211)
hold on
ylabel("Temperature [°C]")
plot(yt1, "DisplayName", "Actual data")
plot(y_hat_kalman(length(yv)+1:end), "DisplayName", "7 step prediction")
title('Temperature Predictions')
datetick("x", "dd mmm")
axis tight
legend
%subplot(312)
% hold on
% plot([yv; yt1])
% plot(y_hat_kalman)
% axis tight
subplot(212)
i = xsave((length(yv)+1:end),:)
plot(i)
title('Parameter Estimations')
axis tight
sgtitle("Test Data 1")

display("Variance of " + k + " steps: " + var(yt1-y_hat_kalman(length(yv)+1:end)))
%acf(yt1-y_hat_kalman(length(yv)+1:end),40,0.05,1)
%etilde1_t1 = yt1-y_hat_kalman(length(yv)+1:end);
%%

[y_hat_kalman, xsave] = superKalman2([yp2; yt2],[xp2; xt2],AK,Re,Rw,Rxx0,x0,p,r,q,k);
subplot(211)
hold on
ylabel("Temperature [°C]")
plot(yt2, "DisplayName", "Actual data")
plot(y_hat_kalman(length(yp2)+1:end), "DisplayName", "7 step prediction")
title('Temperature Predictions')
datetick("x", "dd mmm")
axis tight
legend
%subplot(312)
% hold on
% plot([yp2; yt2])
% plot(y_hat_kalman)
% axis tight
subplot(212)
i = xsave((length(yp2)+1:end),:);
plot(i)
plot(i)
title('Parameter Estimations')
axis tight
sgtitle("Test Data 2")

display("Variance of " + k + " steps: " + var(yt2-y_hat_kalman(length(yp2)+1:end)))
%yhat2_k1 = y_hat_kalman(length(yp2)+1:end);
%%
%% Whole data prediction

Re = eye(length(AK))*1e-8; % Hidden state noise variance matrix
Rw = 1e-2;                  % Guess for observation variance

% Initial values
Rxx0 = eye(length(AK))*1e-2; % Initial variance

k = 7;

[y_hat_kalman, xsave] = superKalman2(data(:,9),data(:,7),AK,Re,Rw,Rxx0,x0,p,r,q,k);

subplot(311)
hold on
plot(data(10000:end,9))
plot(y_hat_kalman(10000:end))
axis tight
subplot(312)
hold on
plot(data(t,9))
plot(y_hat_kalman)
axis tight
subplot(313)
plot(xsave)
axis tight

%%
sgtitle("Validation Data Prediction Using the Kalman Filter")
subplot(311)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k1, "DisplayName", "1 step prediction")
datetick("x", "dd mmm")
legend
axis tight
subplot(312)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k7, "DisplayName", "7 step prediction")
datetick("x", "dd mmm")
legend
axis tight
subplot(313)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k26, "DisplayName", "26 step prediction")
datetick("x", "dd mmm")
legend
axis tight
%%
sgtitle("Test Data 1 Prediction Uing the Kalman Filter")
subplot(211)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k1, "DisplayName", "1 step prediction")
datetick("x", "dd mmm")
legend
axis tight
subplot(312)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k7, "DisplayName", "7 step prediction")
datetick("x", "dd mmm")
legend
axis tight
subplot(313)
hold on
ylabel("Temperature [°C]")
plot(yv, "DisplayName", "Actual data")
plot(yhat_k26, "DisplayName", "7 step prediction")
datetick("x", "dd mmm")
legend
axis tight

%%
% subplot(311)
% 
% acf(etilde1_t1,40,0.05,1)
% title('acf ê_{1}')
% subplot(312)
% acf(etilde1_t7,40,0.05,1)
% title('acf ê_{7}')
% 
% subplot(313)
% acf(etilde1_t26,40,0.05,1)
% title('acf ê_{26}')
% sgtitle("ê for data test 2")