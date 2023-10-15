clear; clc; clf;
close all

%% LAB PART ONE %%
load('Balance_Part1_Data.mat')

[vel_elijah, t_total, COP_Elijah] = Balance(Elijah_Balance_Data, 100, 0.0025);

figure
plot(COP_Elijah(1,:)*100 - mean(COP_Elijah(1,:)*100), COP_Elijah(2,:)*100 - mean(COP_Elijah(2,:)*100), 'LineWidth', 1)

title('Stabilogram for Personal Data')
subtitle(['Average Velocity: ', num2str(round(vel_elijah*100,1)), 'cm/s, Total Time: ', num2str(round(t_total,2)), 's'])
xlabel('x-position of COP [cm]')
ylabel('y-position of COP [cm]')

% view the data in 3-space for fun
%figure
%plot3(COP_Elijah(1,:)*100 - mean(COP_Elijah(1,:)*100), COP_Elijah(2,:)*100 - mean(COP_Elijah(2,:)*100), linspace(0,30,3000),'LineWidth', 1)
%% LAB PART TWO %%
load('Balance_Part2_Data.mat')
h = [0.0025, 0.0025, 0.0025, 0.0025+0.06, 0.0025+0.06];
% Loop through subjects and trials
for i = 1:3
    for j = 1:5
        variable_name = sprintf('subject%d_trial%d', i, j);
        
        % Call the 'Balance' function and store the velocity in the v_avg matrix
        v_avg(i, j) = Balance(eval(variable_name), 100, h(j));
    end
end

% Calculate velocity differences for each subject
for i = 1:3
    velocity_diff(i, 1) = v_avg(i, 2) - v_avg(i, 1); % Trial 2 - Trial 1
    velocity_diff(i, 2) = v_avg(i, 3) - v_avg(i, 2); % Trial 3 - Trial 2
    velocity_diff(i, 3) = v_avg(i, 4) - v_avg(i, 1); % Trial 4 - Trial 1
    velocity_diff(i, 4) = v_avg(i, 5) - v_avg(i, 1); % Trial 5 - Trial 1
end

% Calculate the average of velocity differences across the 3 subjects
avg_diff = mean(velocity_diff*100);

figure;
trials = 1:4;

% Calculate the minimum and maximum velocity differences for each trial
min_diff = min(velocity_diff*100);
max_diff = max(velocity_diff*100);

bar(trials, avg_diff, 0.5, 'FaceColor', [0.5 0.5 0.8]);
hold on;

% Adjust the error bars to extend from min to max values
errorbar(trials, avg_diff, avg_diff - min_diff, max_diff - avg_diff, 'k', 'linestyle', 'none', 'LineWidth', 2);

ylabel('Change in COP Mean Velocity [cm/s]');
title('Change in COP Mean Velocity with Impairment of Different Feedback Systems');
subtitle('Data collected for (n=3) subjects')
labels = {'Impaired Visual System', 'Impaired Vestibular System', 'Impaired Proprioceptive System', 'All Systems Impaired'};

for i = 1:length(trials)
    label = sprintf('%.2f', avg_diff(i));
    text(trials(i)+0.05, avg_diff(i) + 0.1, label, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    text(trials(i), -0.5, labels{i}, 'HorizontalAlignment', 'center'); % Add text labels under each bar
end

xlim([0.5, 4.5]);
xticks([]);  % Hide the numbers on the x-axis
hold off;


%% ANALYSIS FUNCTION %%
function [v_avg, t_total, COP] = Balance(data, sample_freq, h)
% Function calculates 2xn array with x and y COP coordinates, calculates
% mean velocity based on the change in COP and returns total sample time
% for use in further analysis. Inputs: 6xn matrix containing data in
% columns, the force plate sample recording frequency, and h, the distance 
% between the force sensor and the top of the force plate. 

    % Import Data
        % all in N
    Fx = data(:,1);
    Fy = data(:,2);
    Fz = data(:,3);
        % all in N-m
    Mx = data(:,4);
    My = data(:,5);
    Mz = data(:,6);
    
    n = size(data,1);
    
    t_total = n/sample_freq; % seconds
    
    % filter data
    [a,b] = butter(2,10/50);

    Fx = filtfilt(a,b,Fx);
    Fy = filtfilt(a,b,Fy);
    Fz = filtfilt(a,b,Fz);

    Mx = filtfilt(a,b,Mx);
    My = filtfilt(a,b,My);
    Mz = filtfilt(a,b,Mz);

    % assign imported data to arrays
    F(1,:) = Fx;
    F(2,:) = Fy;
    F(3,:) = Fz;

    M(1,:) = Mx;
    M(2,:) = My;
    M(3,:) = Mz;

    % use F and M to calculate COP=(x,y)
    COP = zeros(2,n);
    for i = 1:n
        COP(1,i) = (-h*F(1,i) - M(2,i))/F(3,i); % COP x
        COP(2,i) = (-h*F(2,i) - M(1,i))/F(3,i); % COP y
    end
    
    % use (x,y) points to calculate average velocity
    % lengths in meters
    dl = zeros(1,n);
    L = 0;
    for i = 1:n-1
        dl(i) = sqrt( (COP(1,i+1) - COP(1,i))^2 + (COP(2,i+1)-COP(2,i))^2 );
        L = L + dl(i);
    end
    
    v_avg = L/t_total; % m/s
    
end