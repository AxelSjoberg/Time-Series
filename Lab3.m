% Lab 3


%% 1.1
P11 = 7/8
P22 = P11
P = [P11 (1-P11); (1-P22) P22]
mc = dtmc(P)
numSteps = 999
u = simulate(mc, numSteps)


%% 1.2

% Example of Kalman filter
% Simulate process

% Example of Kalmanfilter
% Simulate process

N = 200;
A = [1 .5 0.8];
sigma2 = 1.5;
e = sqrt(sigma2)*randn(N,1);
y = filter(1,A,e);



% Data length
N = length(y);

% Define the state space equations

Sigma2_e = 1.5;
sigma2_w = 1.5;

A = eye(2); %???
Re = zeros(2); % Hidden state noise covariance matrix %IS THIS RIGHT???
Rw = sigma2_w; % Observation variance % Initital guess here is 1.

% usually C should be set here to,
% but in this case C is a function of time.

% Set some initial values
Rxx_1 = (1) * eye(2); % Initial variance
xtt_1 = A*[0 0]'; % Initial state .     ????

% Vector to store values in
xsave = zeros(2, N);
% Kalman filter. Start from k=3 ,
% since we need old values of y.
for k=3:N,
    % C is, in our case, a function of time.
    C = [-y(k-1) -y(k-2)];
    % Update
    Ryy = C*Rxx_1*C' + Rw;
    Kt = Rxx_1*C'*inv(Ryy);
    xtt = xtt_1 + Kt* (y(k) - C*xtt_1);
    Rxx = (eye(2) - Kt * C)*Rxx_1;
    % Save
    xsave(:, k) = xtt;
    % Predict
    Rxx_1 = A*Rxx*A' + Re;
    xtt_1 = A*xtt;
end;

xtt

%% -------------------------- SECTION 3 --------------------------------
subplot(211)
plot(tar2)
title('tar2')
subplot(212)
plot(thx)
title('thx')
%%
na = 2;
nb = 0;
nk = 0;
model = [na]; 
lambda = 0.95;
[Aest, yhat, covAest, yprev] = rarx(tar2, model, 'ff', lambda);

%%
plot(Aest)
hold on
plot(thx)
% QUESTION 1, the lower the lambda, the faster and noiser the parameter
% estimates become
%%
n = 100;
lambda_line = linspace(0.85, 1, n);
ls2 =  zeros(n, 1);
for i =1:length(lambda_line)
    [Aest, yhat, CovAest, trash] = ...
    rarx(tar2, [2], 'ff', lambda_line(i));
    ls2(i) = sum((tar2 - yhat).^2);
end
plot(lambda_line, ls2)
%find(ls2 == min(ls2(:)))
% QUESTION 2, the square sum of the prediction errors.  
% QUESTION 3, Re = [sigma2_e 0; 0 0] as e = [sigma_e 0]
%%
sigma2_e = 10^(-3);
Re = sigma2_e*[1 0; 0 0];
Rw = 1;
V0 = 10*eye(2);
m0 = [0 0]';
param = kalman_param(tar2, 2, Re, Rw, V0, m0);
%zsave1 = param(1,:)'
%param2 = param(2,:)'
%param = [zsave1 param2]
plot(param')
hold on 
plot(thx)
% QUESTION 4.1, Increasing Re and decreasing Rw increases the speed and varaince. Also it
% is the relationship between Re and Rw that is important not the absolute
% terms. 
% QUESTION 4.2, Yes we got alot better results than RLS. Kalman might be less 
% computionally efficient and we have to check that our poles are inside of our unit circle 

%% ----------------------- SECTION 3.3 -----------------------
clear y e N
sigma2_e = 1;
sigma2_v = 4;
b = 20;
N = 1000;

RW = [1 -1];
e = sqrt(sigma2_e) * randn(N,1);
v = sqrt(sigma2_v) * randn(N,1);
x = filter(1, RW, e);

x = x(100:end);
e = e(100:end);
v = v(100:end);
u = u(100:end);


for k = 1:length(u)
    y(k) = x(k) + b * u(k) + v(k);
end
y = y';
%%

A = [1 0; 0 1];

sigma2_e = 1;
sigma2_v = 4;

Re = sigma2_e * [1 0; 0 0];
Rw = sigma2_v;

ztt_1 = [0;0];
Rxx_1 = 100*eye(2);

clear zsave
for k = 1:length(y)
    ztt = ztt_1;
    C = [1 u(k)];
    Kt = (Rxx_1*C')*((C*Rxx_1*C' + Rw)^(-1));
    ztt = ztt_1 + Kt * (y(k)-C*ztt_1);
    Rxx = Rxx_1 - Kt*C*Rxx_1;
    zsave(:,k) = ztt;
    ztt_1 = ztt;
    Rxx_1 = Rxx + Re; % ignored A here since it is the identity matrix.
end 
zsave = [zsave(1,:)' zsave(2,:)'];
plot(zsave)
hold on
plot(x, '--')
hold on
plot(repmat(20, [length(x) 1]))


% QUESTION 5.2, the state is the param b
% 5.3 Re = [sigma2_e 0; 0 0] and Rw = sigma2_v


%% --------------------- SECTION 3.4 -------------------
plot(svedala94)
s = 6
y_diff = filter([1 zeros(1, s-1) -1], 1, svedala94)

%%
subplot(211)
plot(y_diff)
subplot(212)
plot(svedala94)
%%
subplot(211)
T = linspace(datenum(1994,1,1), datenum(1994, 12, 31),...
    length(svedala94));
plot(T,svedala94);
datetick('x');
subplot(212)
T = linspace(datenum(1994,1,1), datenum(1994, 12, 31),...
    length(svedala94));
plot(T,y_diff);
datetick('x');

%% 
th = armax(y_diff, [2 2]);
th_winter = armax(y_diff(1:540), [2 2]);
th_summer = armax(y_diff(907:1458), [2 2]);
%% i: ALL the data
present(th)
% FPE = 2.603, MSE = 2.593, All params are significant
th.A % [1 -1.6235 0.6551]
th.c % [1 -0.8313 -0.1361]
%% ii: Jan - Mar
present(th_winter)
% FPE = 1.6, MSE = 1.576, All params are significant
th_winter.A % [1 -1.6662 0.7064]
th_winter.C % [1 -0.8331 -0.1078]

%% iii: Jun - Aug
present(th_summer)
% FPE = 3.515, MSE = 3.464, All params are not significant
th_summer.A % [1 0.2545 -0.4598]
th_summer.C % [1 1.0326 0.0427]
% QUESTION 6, Yes it is likley that th_winter and th_summer are different
% processes given our param estimates

%%
th0 = [th_winter.A(2:end) th_winter.C(2:end)]; 
[thr, yhat] = rarmax(y_diff, [2 2], 'ff', 0.99, th0);
subplot(311);
plot(T, svedala94);
datetick('x');

subplot(312);
plot(thr(:,1:2));
hold on
plot(repmat(th_winter.A(2:end), [length(thr) 1]), 'b');
plot(repmat(th_summer.A(2:end), [length(thr) 1]), 'r');
axis tight
hold off

subplot(313)
plot(thr(:,3:end))
hold on
plot(repmat(th_winter.C(2:end), [length(thr) 1]), 'b');
plot(repmat(th_summer.C(2:end), [length(thr) 1]), 'r');
axis tight
hold off
% QUESTION 7, ??

%% ------------------ SECTION 3.5 -------------------------------

y = svedala94(850:1100)
y_mean = mean(y)

%%
s = 6;
t = (1:length(y))';
U = [sin(2*pi*t/s) cos(2*pi*t/s)];
Z = iddata(y ,U);
model = [3 [1 1] 4 [0 0]];
    %[na [nb_1 nb_2] nc [nk_1 nk_2]];


thx = armax(Z,model);
thx.b;

subplot(211)
plot(U * cell2mat(thx.b)')
subplot(212)

plot(y)

% QUESTION 7.1, increasing the length of the season decreases the frequency
% of U*thx.b ??
% QUESTION 7.2, we might need to take a trend into acount

%%
U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z = iddata(y,U);
m0 = [thx.A(2:end) cell2mat(thx.B) 0 thx.C(2:end)]; %% thx.B
Re = diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re ,m0)

plot(yhat)
hold on
plot(y, 'black')
%%

m = thr(:,6);
a = thr(end, 4);
b = thr(end, 5);
y_mean = m + a*U(:, 1)+ b*U(:, 2);
y_mean = [0;y_mean(1:end - 1)];

%%
plot(y)
hold on
plot(y_mean, '--')

% QUESTION 9, They seem rahter similar except that the true y values have
% an added constant. Moreover we also have some peaks in our y_mean that
% are not present in the true values. WHY?

%%

y = svedala94;
y = y-y(1);
%%
s = 6;
t = (1:length(y))';
U = [sin(2*pi*t/s) cos(2*pi*t/s)];
Z = iddata(y ,U);
model = [3 [1 1] 4 [0 0]];
    %[na [nb_1 nb_2] nc [nk_1 nk_2]];


thx = armax(Z,model);
thx.b;


U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z = iddata(y,U);
m0 = [thx.A(2:end) cell2mat(thx.B) 0 thx.C(2:end)]; %% thx.B
Re = diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re ,m0);

%%
y_tilde = (y - yhat)
%acfpacf(y_tilde)
whitenessTest(y_tilde)
% 
%% Due to that the processes for the different season vary quite alot we could redo the analysis above for the different seasons
% SUMMER
y = svedala94(850:1100);
y = y-y(1);
s = 6;
t = (1:length(y))';
U = [sin(2*pi*t/s) cos(2*pi*t/s)];
Z = iddata(y ,U);
model = [3 [1 1] 4 [0 0]];
    %[na [nb_1 nb_2] nc [nk_1 nk_2]];


thx = armax(Z,model);
thx.b;

U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z = iddata(y,U);
m0 = [thx.A(2:end) cell2mat(thx.B) 0 thx.C(2:end)]; %% thx.B
Re = diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re ,m0);

y_tilde = (y - yhat)
whitenessTest(y_tilde)
% Remove the outlier in observation 8 perhaps by interpolation or something
% var(y_tilde) = 9.7612
%%
% WINTER
y = svedala94(1:450);
y = y-y(1);
s = 6;
t = (1:length(y))';
U = [sin(2*pi*t/s) cos(2*pi*t/s)];
Z = iddata(y ,U);
model = [3 [1 1] 4 [0 0]];
    %[na [nb_1 nb_2] nc [nk_1 nk_2]];


thx = armax(Z,model);
thx.b;

U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))];
Z = iddata(y,U);
m0 = [thx.A(2:end) cell2mat(thx.B) 0 thx.C(2:end)]; %% thx.B
Re = diag([0 0 0 0 0 1 0 0 0 0]);
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]];
[thr, yhat] = rpem(Z, model, 'kf', Re ,m0);

y_tilde = (y - yhat)
whitenessTest(y_tilde)

% Same as above!

