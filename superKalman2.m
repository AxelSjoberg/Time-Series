% y    : The output vector that we're predicting
% u    : Input vector
% A    : Kalman A matrix, here we use the identity matrix
% Re   : The state noise variance matrix
% Rw   : Guess for measurment variance
% Rxx0 : Initial guess for state variance matrix
% x0   : Initial guess for state vector
% p    : Order of Box Jeniks A-polynomial
% r    : Order of Box Jeniks B-polynomial
% q    : Order of Box Jeniks C-polynomial
% k    : Number of time steps we predict ahead

function [y_hat, xsave] = superKalman2(y, u, A, Re, Rw, Rxx0, x0, p, r, q, k)

if nargin < 11
    k = 1;
end

xtt_1 = x0;   % Initial state guess
Rxx_1 = Rxx0; % Initial state variance

len   = length(A); 
N     = length(y); 

y_hat  = zeros(N,1);          % Prediction vector
e      = zeros(N,1);          % Stores prediction errors
y_temp = zeros(N+k,1);        % Temporarily stores predictions
y_temp(1:len) = y(1:len);     % Sets the first elements = y
xsave  = zeros(N-k,length(x0)); % Matrix for storing our states

for t=max([p,r,q])+1:N-k
    % C is a function of time, containing the following
    %    The previous p outputs 
    %    The previous r inputs 
    %    The previous q prediction errors
    C = [-y(t-(1:p)); u(t-(0:r)); e(t+1-(1:q))]';
    
    y_temp(t) = y(t);           % Now we know the current output
    e(t)      = y(t) - C*xtt_1; % Our best guess for the prediction error
    
    % Update
    Ryy = C*Rxx_1*C' + Rw; 
    Kt  = Rxx_1*C'/Ryy;
    xtt = xtt_1 + Kt*e(t);
    Rxx = Rxx_1 - Kt*C*Rxx_1;
    
    % Predict
    Rxx_1 = A*Rxx*A' + Re;
    xtt_1 = A*xtt;
    
    % For loop since we're predicting k steps forward
    for tp = t + (1:k) % Time of prediction
        % Ck consists of:
        %    The p previous best predictions of y
        %    The current and r-1 previous inputs
        %    The q previous prediction errors, or zeros if unknown
        Ck = [-y_temp(tp-(1:p)); u(tp-(0:r)); e(tp-(1:q))]';
        y_temp(tp) = Ck*xtt; % Current prediction k steps forward
    end
    y_hat(t+k) = y_temp(t+k); % Our prediction k steps ahead
    
    xsave(t,:) = xtt; % Save states
end

end