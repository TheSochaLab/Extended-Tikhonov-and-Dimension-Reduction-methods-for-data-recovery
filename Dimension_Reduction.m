%  Dimension Reduction method for large datasets
%
%  Input:            h  -  The impulse response of the system.
%                    y  -  Recorded data
%                    m  -  Regularization parameter
%                    n  -  Window size
%                    N  -  Accepting number from each window (N<n)
%         SamplingRate  -  initial guess 
%
%  Optional Intputs:
%        BaseThreshold  -  all the estimation below this value will be significantly suppressed
%              u_exact  -  True input. If the exact solution is known, we can compute the estimation error
%
%  Output:
%          U_estimated  -  Estimation of the input (the solution)
%                 Enrm  -  norm of the error (if x_exact is known)
%
%  Ref: Hodjat Pendar, John J. Socha, & Julianne Chung, New 
%       methods for recovering signals in physiological systems with large 
%       datasets, submitted manuscript, 2016.
%
% By: Hodjat Pendar, John J. Socha & Julianne Chung 
%     Virginia Tech
%
%%

function [U_estimated, Enrm] = Dimension_Reduction(h, y,m, n, N, SamplingRate, BaseThreshold , u_exact)


%% Dimension correction
n = fix(n/m)*m;

%% Read the impulse response
I = h / sum(h) * SamplingRate;   % Normalizing the impulse response 

if length(I)<= n     % Unifying the length of the data and the impulse response.
    I = [I; zeros(n-length(I),1)];
else
    I = I(1:n);
end


%% Dimension Reduction method

% Set up the matrices:
U=zeros(size(y));    % Initializing the input estimation vector

% optional regularization matrix
Gamma = 0.000000001;
q1 = zeros(n/m,1); q1(1:3)=[1 -2 1]';
q2 = zeros(n/m,1); q2(1)=1;
Q = toeplitz(q1,q2);   % regularization matrix


% H matrix (Equation (1.4))
a = I; 
b = a*0; b(1) = a(1);

H = toeplitz(a,b)/SamplingRate;

% L matrix
L = zeros(n,fix(n/m));
for i = 1:fix(n/m)
    L(m*(i-1)+1:i*m,i)=1;     % Equation (3.3)
end

% M matrix
M = L*inv(L'*H'*H*L + Gamma *Q'*Q)*L'*H'; % Equation (3.5)
M_ = M;

for i = 2:m;
    M_(1:end-i+1, 1:end-i+1) = M_(1:end-i+1, 1:end-i+1) + M(i:end, i:end);
end

M = M_/m;  % Equation (3.7)

%% Calculate the solution
Y=y;
time = 0 : 1/SamplingRate : (length(y)-1)/SamplingRate;


for i = 1 : fix((length(y))/N)
    
    if length(Y) >= n
    U_ = M*Y(1:n);  % Estimation of the input. Equation (3.6)
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


%% Output of the function:

U_estimated = U;

if ~isempty(u_exact)
    Enrm = norm(u_exact- U_estimated);
end





