clear
close all
clc
%===simulation time===
simulationTime = 6; %in milliseconds
deltaT=.001;
Fs = 1/deltaT;
t=0:deltaT:simulationTime;
vc = 0;

% Vrest = -60;
% drugs application
pctTTX = 0;
pctTEA = 0;
TTX = 1-pctTTX;
TEA = 1 - pctTEA;

%===constant parameters===%
%All of these can be found in Table 3
C_out_Na = 460;
C_in_Na = 50;
C_out_K = 22;
C_in_K = 400;
% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 6.3; % in kelvin

E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
rat = log(C_out_Na/C_in_Na)/log(C_out_K/C_in_K);
E_L=10.6;
gbar_K=36*TEA; 
gbar_Na=120*TTX; 
g_L=.3;
C=1;

% Vrest = K*((gbar_K*C_out_K+gbar_Na*C_out_Na)/(gbar_K*C_in_K+gbar_Na*C_in_Na))
%===set the initial states===%
V=0; %Baseline voltage
p_K = 1;
p_Na = 0.05;
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

% Vrest = (115+12*rat)/(rat-1);

% V = V + Vrest;

alpha_n = (.01 * ( (10-V) / (exp((10-V)/10)-1) ));
beta_n = .125*exp(-V/80);
alpha_m = (.1*( (25-V) / (exp((25-V)/10)-1) ));
beta_m = 4*exp(-V/18);
alpha_h = .07*exp(-V/20);
beta_h = 1/(exp((30-V)/10)+1);

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18
I_Na = zeros(1,length(t)); %Equations 3 and 14
I_K = zeros(1,length(t)); %Equations 4 and 6
I_L = zeros(1,length(t)); %Equation 5
I_ion = zeros(1,length(t));

I = zeros(1,numel(t));

% Define Vrest
% g_Na = (gbar_Na*(m(1).^3).*h(1)); % .01
% g_K = (gbar_K*(n(1).^4)); % .039
% g_L = .3;


% Vrest = I_ion/g_ion;
% Vrest = ((g_Na*E_Na_orig)+(g_K*E_K_orig))/(g_Na+g_K);
E_K = E_K_orig - Vrest; 
E_Na = E_Na_orig - Vrest;
% E_K = -12;
% E_Na = 115;



% V = zeros(1, numel(t));
%% ======== define stimuli =====================
if vc
    Vrest = -60;
    voltageLevels= 60;
    
    stimS = 0.5; % in ms
    stimE = 5.5; % in ms
    stimDur = round((stimE-stimS))*Fs;
    V(round(stimS*Fs):round(stimE*Fs)) = voltageLevels; %I(2001:numel(t)) = currentLevels;
else
    %Set externally applied current across time
    %Here, first 500 timesteps are at current of 50, next 1500 timesteps at
    %current of zero (resets resting potential of neuron), and the rest of
    %timesteps are at constant current
    currentLevels1 = [100]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
    
    stimS1 = 0.6; % in ms
    stimE1 = 0.8; % in ms
    stimDur1 = round((stimE1-stimS1)*Fs);
    I(round(stimS1*Fs):round(stimE1*Fs)) = currentLevels1; %I(2001:numel(t)) = currentLevels;
    
    currentLevels2 = [100]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
    
    stimS2 = 5; % in ms
    stimE2 = 5.5; % in ms
    stimDur2 = round((stimE2-stimS2)*Fs);
%     I(round(stimS2*Fs):round(stimE2*Fs)) = currentLevels2; %I(2001:numel(t)) = currentLevels;
    
    %plot the stimulus
    figure
    plot(t,I,'r')
    %Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
    %I(1:numel(t)) = currentLevels;
end

for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
    
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = (.01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) ));
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = (.1*( (25-V(i)) / (exp((25-V(i))/10)-1) ));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    
    
    %---calculate the currents---%
    I_Na(i) = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K(i) = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L(i) = g_L *(V(i)-E_L); %Equation 5
    I_ion(i) = I(i) - I_K(i) - I_Na(i) - I_L(i);
    
    
    
    %---calculate the derivatives using Euler first order approximation---%
    if ~vc
        V(i+1) = V(i) + deltaT*I_ion(i)/C; 
    end
    %     I_ion(i) = (C*(-V(i)+V(i+1)))/deltaT;
    n(i+1) = (n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i))); %Equation 7
    m(i+1) = (m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i))); %Equation 15
    h(i+1) = (h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i))); %Equation 16
    
end


V = V+Vrest; %Set resting potential to -70mv

%===plot Voltage===%
figure
plot(t,V,'LineWidth',3)
hold on
legend({'voltage'})
ylabel('Voltage (mv)')
xlabel('time (ms)')
title('Voltage over Time in Simulated Neuron')


%===plot Conductance===%
figure
p1 = plot(t,gbar_K*n.^4,'LineWidth',2);
hold on
p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
ylabel('Conductance')
xlabel('time (ms)')
title('Conductance for Potassium and Sodium Ions in Simulated Neuron')

% ======== plot n, m and h ================
figure
p1 = plot(t,n,'LineWidth',2);
hold on
p2 = plot(t,m,'r','LineWidth',2);
p3 = plot(t,h,'k','LineWidth',2);
title('m h and n')
legend('n','m','h')

% ========== plot individual currents ================
figure
p1 = plot(t,I_Na,'r','LineWidth',2);
hold on
p2 = plot(t,I_K,'b','LineWidth',2);
Itot = I_Na+I_K+I_L-I;
p3 = plot(t,Itot,'k','LineWidth',2);
title('Currents')
legend('Na', 'K', 'Total')

