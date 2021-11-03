
alpha_n = (.01 * ( (10-V) / (exp((10-V)/10)-1) ));
beta_n = .125*exp(-V/80);
alpha_m = (.1*( (25-V) / (exp((25-V)/10)-1) ));
beta_m = 4*exp(-V/18);
alpha_h = .07*exp(-V/20);
beta_h = 1/(exp((30-V)/10)+1);

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18
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
        V(i+1) = V(i) + deltaT*I_ion(i)/C;

    %     I_ion(i) = (C*(-V(i)+V(i+1)))/deltaT;
    n(i+1) = (n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i))); %Equation 7
    m(i+1) = (m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i))); %Equation 15
    h(i+1) = (h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i))); %Equation 16
    
end
