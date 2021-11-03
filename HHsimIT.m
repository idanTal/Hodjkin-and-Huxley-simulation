%===simulation time===
simulationTime = 30; %in milliseconds
deltaT=.01;
Fs = 1/deltaT;
t=0:deltaT:simulationTime;


%===specify the external current I===
changeTimes = [0]; %in milliseconds
currentLevels = [0]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)

%Set externally applied current across time
%Here, first 500 timesteps are at current of 50, next 1500 timesteps at
%current of zero (resets resting potential of neuron), and the rest of
%timesteps are at constant current
I = zeros(1,numel(t));
stimS = 2; % in ms
stimE = 6; % in ms
stimDur = round((stimE-stimS))*Fs;
I(round(stimS*Fs):round(stimE*Fs)) = currentLevels; %I(2001:numel(t)) = currentLevels;
%plot the stimulus
figure
plot(t,I,'r')
%Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
%I(1:numel(t)) = currentLevels;


%===constant parameters===%
%All of these can be found in Table 3
gbar_K=36; gbar_Na=120; g_L=.3;
E_K = -12; E_Na=115; E_L=10.6;
C=1;


%===set the initial states===%
V=0; %Baseline voltage

alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) ); %Equation 12
beta_n = .125*exp(-V/80); %Equation 13
alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) ); %Equation 20
beta_m = 4*exp(-V/18); %Equation 21
alpha_h = .07*exp(-V/20); %Equation 23
beta_h = 1/(exp((30-V)/10)+1); %Equation 24

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18

I_Na = zeros(1,length(t)); %Equations 3 and 14
I_K = zeros(1,length(t)); %Equations 4 and 6
I_L = zeros(1,length(t)); %Equation 5
I_ion = zeros(1,length(t));

% drugs application
pctTTX = 0;
pctTEA = 0;
TTX = 1-pctTTX;
TEA = 1 - pctTEA;

for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
    
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = TEA*(.01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) ));
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = TTX*(.1*( (25-V(i)) / (exp((25-V(i))/10)-1) ));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    
    
    %---calculate the currents---%
    I_Na(i) = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K(i) = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L(i) = g_L *(V(i)-E_L); %Equation 5
    I_ion(i) = I(i) - I_K(i) - I_Na(i) - I_L(i); 
    
    
    
    %---calculate the derivatives using Euler first order approximation---%
    V(i+1) = V(i) + deltaT*I_ion(i)/C;
    n(i+1) = (n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i))); %Equation 7
    m(i+1) = (m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i))); %Equation 15
    h(i+1) = (h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i))); %Equation 16

end


V = V-70; %Set resting potential to -70mv

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
% p3 = plot(t,I_ion,'k','LineWidth',2);
title('Currents')
legend('Na', 'K', 'Total')

