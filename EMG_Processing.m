clear; clc; clf; close all

load('forcedelaydata.mat')

% Create vectors of data to use in processing
t = table2array(Forcesub6(:,1));
F = table2array(Forcesub6(:,2));
emg = abs(table2array(Forcesub6(:,3))); % rectify the emg data

% low pass filter
[a,b] = butter(2,15/50);
emg = filtfilt(a,b,emg);

% rms envelope
T = 50; % envelope moving mean time frame
emg = sqrt(movmean(emg.^2, T));

a = 450;
% plot results
plot(t,[a*emg,F], 'LineWidth', 1)

% average over range {t = 5.74 - 8.78} range for at rest picked manually 
emg_rest = emg(t > 5.74 & t < 8.78);
F_rest = F(t > 5.74 & t < 8.78);

emg_rest_avg = mean(emg_rest);
F_rest_avg = mean(F_rest);

s_emg = std(emg_rest);
s_F = std(F_rest);

% calculate thresholds for force and emg signal activation
emg_lim = emg_rest_avg + 4*s_emg;
F_lim = F_rest_avg + 4.5*s_F;

emg_ind = [];
F_ind = [];

for i = 2:length(t)
    if emg(i) > emg_lim & emg(i-1) < emg_lim
        if emg(i+100) > emg_lim
            emg_ind = [emg_ind, i];
        end

    else
        emg_ind = [emg_ind];
        
    end
end




F_ind = [];

for i = 2:length(t)
    if F(i) > F_lim & F(i-1) < F_lim
        if F(i+100) > F_lim
            F_ind = [F_ind, i];
        end
        
    else
        F_ind = [F_ind];
        
    end
end


figure
hold on
xlim([-15,70])
plot(t,[a*emg,F], 'LineWidth', 1)
emg_line = yline(a*emg_lim, '-.k', ['EMG threshold: ' num2str(a*emg_lim) 'mV'], 'LineWidth', 1);
F_line = yline(F_lim, '-.k', ['Force threshold: ' num2str(F_lim) 'N'], 'LineWidth', 1);
plot(t(emg_ind), a*emg(emg_ind),'xk','MarkerSize', 8, 'LineWidth', 1.5, 'Color', '#77AC30')
plot(t(F_ind), F(F_ind), 'xk', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color','#7E2F8E')
emg_line.LabelHorizontalAlignment = 'Left';
emg_line.Color = '#77AC30';
F_line.LabelHorizontalAlignment = 'Left';
F_line.LabelVerticalAlignment = 'Bottom';
F_line.Color = "#7E2F8E";
xline(t(emg_ind),'-.','LineWidth', 1, 'Color', '#77AC30')
xline(t(F_ind), '-.','LineWidth', 1, 'Color', '#7E2F8E')

for i = 1:length(emg_ind)
    delays(i) = t(F_ind(i)) - t(emg_ind(i));
end

disp(num2str(1000*delays'))

title('Post-Processed EMG and Grip Force Over Time with Onset Thresholds Marked')
subtitle('Close up view of EMG-Force delay for first 3 contractions')
xlabel('Time (seconds)')
ylabel('Force (N) and 450*EMG (mV)')
legend('450xEMG-RMS Signal','Force')

hold off

