%  Recovering true signals from indirectly recorded signals
%
%  In this code, we use the ‘Extended Tikhonov' method and the 'Dimension 
%  Reduction' method to recover true signals from indirectly recorded 
%  signals. 
%
%  Ref: Hodjat Pendar, John J. Socha, & Julianne Chung, New 
%       methods for recovering signals in physiological systems with large 
%       datasets, submitted manuscript, 2016.
%
% By: Hodjat Pendar, John J. Socha & Julianne Chung 
%     Virginia Tech

%% initializing
clc
clear
close all

%% Parameters
SamplingRate = 10;  % The sampling rate of data collection (Hz)

NoiseLevel = 0.0001;  % Standard deviation of the signal noise: To determine this number, run your system without any animal in it and record the noise of the system. Then, determine the SD (standard deviation) of the noise.

% regularization parameters: the paper explains how to determine the regularization parameter.
Gamma = 0.0000001;   % Tikhonov regularization parameter. 
m = 4;               % Dimension Reduction parameter

BaseThreshold = 20;  % Threshold to cancel the noise: Any value below this threshold will be significantly filtered. To find this threshold, apply the methods on recorded noise. The recovery methods amplify the noise. This threshold could be the SD of the amplified noise or its maximum (it's optional).

n = 1500;  % Window size: the methods will be applied on 'n' data points at each iteration. 'n' should be more than the length of the impulse response of the system. 
N = 900;   % Accepted part: after applying the methods on 'n' data points, only the first 'N' points will be accepted. 'n-N' should be less than the length of the impulse response.


%% impulse response
h = textread('ImpulseResponse.txt');   % Put the name of the impulse response file here. It should be in txt format.
h = h(:,2);

%% Exact input and output
Data = textread('Data.txt');  % Put the text file name that includes your data here.
Time = Data(:,1);
y = Data(:,2);

u_exact = y*0; % Replace this with the exact input if it is known

%% Extended Tikhonov Method
[U_estimated_Tikh, Enrm_Tikh] = Tikhonov_Extension(h, y, Gamma, n, N, SamplingRate, BaseThreshold , u_exact);  % Here we call the 'Tikhonov_Extension' function and the estimated inputs will be saved in the 'U_estimated_Tikh' vector. 

RecoveredData = [Time, U_estimated_Tikh];

save Recovered_Tikhonov.txt RecoveredData -ASCII   % Here we save the estimated input in a text file: 'Recovered_Tikhonov.txt'. The first column of the text file is time and the second column is the recovered data.

%% Dimension Reduction Method
Delay = find(h > 0.0001*max(h),1,'first');  % To avoid singularity, we shift the data to get rid of the delay in the data. Here we find the delay in the system.
h_NoDelay = h(Delay:end);                   % Removing the delay from the impulse response. 
u_NoDelay_exact = u_exact(1:end - Delay+1); % Removing the delay from the exact solution (if it is available).
y_NoDelay = y(Delay:end);                   % Removing the delay from the data.

[U_estimated_DR, Enrm_DR] = Dimension_Reduction(h_NoDelay, y_NoDelay, m, n, N, SamplingRate, BaseThreshold , u_exact);   % Recovering the true inputs, using the dimension reduction method
U_estimated_DR = [U_estimated_DR; zeros(Delay-1,1)];   % Adding the removed delay after finding the true inputs.


RecoveredData = [Time, U_estimated_DR];

save Recovered_DimensionReduction.txt RecoveredData -ASCII   % Here we save the estimated input in a text file: 'Recovered_DimensionReduction.txt'. The first column of the text file is time and the second column is the recovered data.

%% Plots:

close all

%%%%%%%%%%%%%%%%%%
% Tikhonov Method:
u_esmt = U_estimated_Tikh;
u = u_exact;

% time
t = Time(1:length(y));

FigHandle = figure;
set(FigHandle, 'Position', [1, 500, 1280, 400]);

plot(t,u,'r');
hold on
plot(t,y, 'LineWidth',2)
plot(t,u_esmt,'m', 'LineWidth',2)


xlim([t(1) t(end)]); 
ylim([-.2*max(y) 2*max(y)])
xlabel('Time (s)', 'FontSize', 25, 'FontName', 'Times New Roman')
ylabel('CO_2 concentration (ppm)', 'FontSize', 25, 'FontName', 'Times New Roman')

set(gca, 'FontSize',20)

if max(abs(u_exact))>0
    legend('True input', 'Recorded data', 'Recovered - Tikhonov Method')
else
    legend('True input (not available)', 'Recorded data', 'Recovered - Tikhonov Method')
end

%%%%%%%%%%%%%%%%%%
% Tikhonov Method:
u_esmt = U_estimated_DR;

FigHandle = figure;
set(FigHandle, 'Position', [1, 1, 1280, 400]);

plot(t,u,'r');
hold on
plot(t,y, 'LineWidth',2)
plot(t,u_esmt,'m', 'LineWidth',2)

xlim([t(1) t(end)]); 
ylim([-.2*max(y) 2*max(y)])
xlabel('Time (s)', 'FontSize', 25, 'FontName', 'Times New Roman')
ylabel('CO_2 concentration (ppm)', 'FontSize', 25, 'FontName', 'Times New Roman')



set(gca, 'FontSize',20)

if max(abs(u_exact))>0
    legend('True input', 'Recorded data', 'Recovered - Tikhonov Method')
else
    legend('True input (not available)', 'Recorded data', 'Recovered - Tikhonov Method')
end








