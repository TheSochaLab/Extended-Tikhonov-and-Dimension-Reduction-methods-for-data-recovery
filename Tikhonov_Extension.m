%  Extension of the Tikhonov method for large datasets
%
%  Input:            h  -  The impulse response of the system.
%                    y  -  Recorded data
%                Gamma  -  Regularization parameter
%                    n  -  Window size
%                    N  -  Accepting number from each window (N<n)
%         SamplingRate  -  The sampling rate of the recorded data 
%
%  Optional Intputs:
%        BaseThreshold  -  all the estimation below this value will be significantly suppressed
%              u_exact  -  The exact input signal. if the exact solution is known, we can compute the estimation error
%
%  Output:
%          U_estimated  -  Estimation of the input. The recovered signal will be saved in this variable
%                 Enrm  -  norm of the error (if x_exact is known)
%
%  Ref: Hodjat Pendar, John J. Socha, & Julianne Chung, New 
%       methods for recovering signals in physiological systems with large 
%       datasets, submitted manuscript, 2016.
%
% By: Hodjat Pendar, John Jake Socha & Julianne Chung 
%     Virginia Tech
%
%%

function [U_estimated, Enrm] = Tikhonov_Extension(h, y,Gamma, n, N, SamplingRate, BaseThreshold , u_exact)

%% Read the impulse response
I = h / sum(h) * SamplingRate;   % Normalizing the impulse response 

if length(I)<= n     % Unifying the length of the data and the impulse response.
    I = [I; zeros(n-length(I),1)];
else
    I = I(1:n);
end


%% Extended Tikonov method

% Set up the matrices:
U = zeros(size(y));   % Initializing the estimated input vector

q1 = zeros(n,1); q1(1:3)=[1 -2 1]';
q2 = zeros(n,1); q2(1)=1;
Q = toeplitz(q1,q2);      % Q is the regularization matrix



a = I; 
b = a*0; b(1) = a(1);

H = toeplitz(a,b)/SamplingRate;
M = inv(H'*H + Gamma* Q'*Q) * H';   % Equation (2.3)

% Calculating the solution
Y=y;
time = 0 : 1/SamplingRate : (length(y)-1)/SamplingRate;   % time vector
DataPoints = length(time);


for i = 1 : fix((length(y))/N)
    
    if length(Y) >= n
    U_ = M*Y(1:n);  % estimation of the input. Equation (2.2)
    Y1 = conv(I, U_(1:N))/SamplingRate;
    
    if length(Y) >= length(Y1)
        Y(1:length(Y1)) = Y(1:length(Y1)) - Y1;
    end
    
    Y=Y(N+1:end);
    
    U_acceptable = smooth(U_(1:N), 5,'moving');
    U_acceptable(U_acceptable<BaseThreshold)=smooth(U_acceptable(U_acceptable<BaseThreshold),50,'moving');

    U((i-1)*N+1:i*N) = U_acceptable;
    
    end
end


%% Calculate the the outputs
U_estimated = U;

if ~isempty(u_exact)
    Enrm = norm(u_exact- U_estimated);
end

