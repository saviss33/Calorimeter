% Created by Sarjot Singh %

% Goal is to analyze calorimeter data from sample file
% then graph, extrapolate and find errors

% Environment
clear;
clc;

%% READING DATA FROM SAMPLE CALORIMETER DATA %
data = readtable('Sample');
% TCM1 and TCM2 are calorimetry
time = data.Var2;
TCM1 = data.Var3;
TCM2 = data.Var6;
% TCW is the boiling water temperature
TW = data.Var4;
% TA is the air temperature
TA = data.Var5;
% the average calorimetry temperature of both calorimeters
TC_avg_data = (TCM1+TCM2)./2;
TC_avg_data_error = std(TC_avg_data/sqrt(length(TC_avg_data)));



%% REGRESSION LINES %
% for regression line at presample
pretime = time(1:425,1);
TC_pre = TC_avg_data(1:425,1);
N_pre = length(pretime);
H_pre = [ones(N_pre,1),pretime];
y_pre = [TC_pre];
W_pre = eye(N_pre);
P_pre = inv(H_pre'*W_pre*H_pre);
x_hat_pre = P_pre*H_pre'*W_pre*y_pre;
A_pre = x_hat_pre(1);
B_pre = x_hat_pre(2);

% calculate sigma y for presample
sum = 0;
    for i=1:N_pre
     sum = sum+((y_pre(i)-A_pre-(B_pre*pretime(i)))^2);
    end
sigma_y_pre = sqrt((1/(N_pre-2))*sum);
    for i=1:N_pre
     W_pre(i,i) = 1/(sigma_y_pre^2);
    end
P_pre = inv(H_pre'*W_pre*H_pre);

% uncertainty in A and uncertainty in B from P matrix in presample
A_pre_error = sqrt(P_pre(1));
B_pre_error = sqrt(P_pre(4));
sigma_q_pre = sqrt((A_pre_error)^2 + (593.663)^2 * (B_pre_error)^2);

% for regression line at near equilbrium of the sample data
endtime = time(495:714,1);
TC_end = TC_avg_data(495:714,1);
N_end = length(endtime);
H_end = [ones(N_end,1),endtime];
y_end = [TC_end];
W_end = eye(N_end);
P_end = inv(H_end'*W_end*H_end);
x_hat_end = P_end*H_end'*W_end*y_end;
A_end = x_hat_end(1);
B_end = x_hat_end(2);

% calculate sigma y for end of sample data
sum_2 = 0;
    for i=1:N_end
     sum_2 = sum_2+((y_end(i)-A_end-(B_end*endtime(i)))^2);
    end
sigma_y_end = sqrt((1/(N_end-2))*sum_2);
    for i=1:N_end
     W_end(i,i) = 1/(sigma_y_end^2);
    end
P_end = inv(H_end'*W_end*H_end);

% uncertainty in A and uncertainty in B from P matrix in presample
A_end_error = sqrt(P_end(1));
B_end_error = sqrt(P_end(4));
sigma_q_end = sqrt((A_end_error)^2 + (728.795)^2 * (B_end_error)^2);



%% ERROR FOR TEMPERATURES %
% For T0, choosing temperature 28.3725 at time 600.76
% For T_ext, extrapolated line from pre-sample 32.6885 at time 600.76
% For T2, choosing temperature 32.6746 at time 619.258 based off avg temperature of two calorimeters
% For T_avg, the temperature is estimated around 30.6255 at time 619.258 on the data
T0 = 28.3725;
T0_error = sqrt((A_pre_error)^2 + (600.76)^2 * (B_pre_error)^2);
T_ext = 32.6885;
T_ext_error = sqrt((A_end_error)^2 + (600.76)^2 * (B_end_error)^2);
T1 = mean(TW(1:444));
T1_error = std(TW(1:444));
T2 = 32.6746;
T2_error = sqrt((A_end_error)^2 + (619.258)^2 * (B_end_error)^2);
T_avg = (T0+T_ext)/2;
T_avg = 30.6255;
T_avg_t = 619.258;



%% SPECIFIC HEAT %
Cc_av = 0.895; % calorimeter specific heat [J/gC]
m_c = 510; % calorimeter mass [g]
m_s = 88.897; % sample D mass [g]
% equation of specific heat of sample
Cs_av = (m_c * Cc_av * (T2-T0))/(m_s * (T1-T2));
% uncertainty for specific heat with error propagation general formula
%delta_m_c = (Cc_av * (T2-T0))/(m_c * (T1-T2));
m_c_error = 5;
%delta_m_s = ((-m_s) * Cc_av * (T2-T0))/((m_c)^2 * (T1-T2));
m_s_error = 0.005;
partial_T0 = ((-m_c) * Cc_av)/(m_s * (T1-T2));
partial_T1 = ((-m_c) * Cc_av * (T2-T0))/(m_s * (T1-T2)^2);
partial_T2 = (m_c * Cc_av * (T1-T0))/(m_s * (T1-T2)^2);
% error of specific heat
Cs_av_error = sqrt((partial_T0 * T0_error)^2 + (partial_T1 * T1_error)^2 + (partial_T2 * T2_error)^2);



% PRINTING VALUES AREA %
fprintf('The temperature when the sample is added: %0.3f ± %0.3f [°C]\n',T0, T0_error);
fprintf('The temperature extrapolated from the sample: %0.3f ± %0.3f [°C]\n',T_ext, T_ext_error);
fprintf(['The average temperature of the initial temperature and extrapolated temperature ' ...
    'is: %0.3f [°C] at time: %0.3f\n'],T_avg, T_avg_t);
fprintf('The final temperature is: %0.3f ± %0.3f [°C]\n',T2,T2_error);
fprintf('The specific heat of Sample: %0.3f ± %0.3f [J/g°C]\n',Cs_av,Cs_av_error);
fprintf('From the data gathered, Sample is made of zinc with a specific heat of 0.402 [J/g°C]\n');



%% GRAPHING PORTION %
% extrapolate domains
pretime_extrapolate = time(1:495,1);
endtime_extrapolate = time(375:714,1);

% linear fit line bottom
pre_line = A_pre + B_pre*pretime_extrapolate;
% linear fit line top
end_line = A_end + B_end*endtime_extrapolate;

hold on
title('Time vs Average Calorimetry Temperature');
xlabel('Time (s)');
ylabel(['Temperature (' char(176) 'C)']);
plot(time, TC_avg_data, '.-.k');
plot(pretime_extrapolate, pre_line,'m','LineWidth',2);
plot(endtime_extrapolate, end_line,'b','LineWidth',2);
xline(600.76,'r')
plot(619.258,30.6255,'g.','MarkerSize',20)
xline(619.258,'g');
yline(30.6255,'g');
legend('Avg Calorimeter Data','Lower Regression Line', ...
    'Upper Regression Line','Sample Added Line','Point of Avg Temp', ...
    'Point of Avg Temp Lines','Location','best');
hold off
