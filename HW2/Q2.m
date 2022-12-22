T = 500; % 500 ms
delT = 0.1; % 0.1ms
lambda = 1e-3; % poisson process rate (per milli-second)
p = lambda*delT; % expected number of spikes in delT length interval
t = 0:delT:T;
Ns = 100;
N = size(t, 2);
pd = makedist('Binomial',2, p);
spikes = random(pd, Ns, N);

I0 = 1; % 1pA
w0 = 50;
sigmaw = 5;
tau = 15; % 15 ms
tau_s = tau/4;

pd = makedist('Normal',w0, sigmaw);
w = random(pd, 1, Ns);
I_app = zeros(1,N);
for i = 1:N
    for j = 1:Ns
        tm = t(spikes(j,1:i) == 1);
        if size(I0*w(j)*sum(exp(-1*(t(i) - tm)/tau) - exp(-1*(t(i) - tm)/tau_s), 2), 1) ~= 0
            I_app(i) = I_app(i) + I0*w(j)*sum(exp(-1*(t(i) - tm)/tau) - exp(-1*(t(i) - tm)/tau_s), 2);
        end
    end
end

% AEF RS Neuron
% equilibrium values from previous HW
Vdc = -69.9999;
Udc = 2.0000e-04;
% neuronal parameters
C = 200;
Gl = 10;
El = -70;
Vt = -50;
delt = 2;
a = 2;
tau = 30;
b = 0;
Vr = -58;
neuronType = "RS";

V = zeros(1,N);
U = zeros(1,N);

V(1) = Vdc;
U(1) = Udc; 
spikecount = 0;
for ind = 1:(N-1)
   dVdt = (-Gl*(V(ind)-El) + Gl*delt*exp((V(ind)-Vt)/delt) - U(ind)+I_app(ind))/C  
   dUdt = (a*(V(ind)-El) - U(ind))/tau;

   V(ind+1) = V(ind) + dVdt*delT;
   U(ind+1) = U(ind) + dUdt*delT; 
   
   if V(ind+1) >= 0
       V(ind+1) = Vr;
       U(ind+1) = U(ind) + b;
       spikecount = spikecount + 1;
   end
end


figure;
plot(t, I_app); % t in ms, I in pA
title(join(['$I_{app}$ vs Time, Number of spikes = ', spikecount, ', $w_0 = $', w0, ", $\sigma_w = $", sigmaw]), 'interpreter', 'latex')
ylabel('Current in pA', 'interpreter', 'latex')
xlabel('Time in ms', 'interpreter', 'latex')

figure
plot(t, V)
title(join(['Voltage vs Time for neuron type ', neuronType, ', Number of spikes = ', spikecount, ', $w_0 = $', w0, ", $\sigma_w = $", sigmaw]), 'interpreter', 'latex')
ylabel('Voltage in mV', 'interpreter', 'latex')
xlabel('Time in ms', 'interpreter', 'latex')

I0 = 1; % 1pA
w0 = 250;
sigmaw = 25;
tau = 15; % 15 ms
tau_s = tau/4;

pd = makedist('Normal',w0, sigmaw);
w = random(pd, 1, Ns);
I_app = zeros(1,N);
for i = 1:N
    for j = 1:Ns
        tm = t(spikes(j,1:i) == 1);
        if size(I0*w(j)*sum(exp(-1*(t(i) - tm)/tau) - exp(-1*(t(i) - tm)/tau_s), 2), 1) ~= 0
            I_app(i) = I_app(i) + I0*w(j)*sum(exp(-1*(t(i) - tm)/tau) - exp(-1*(t(i) - tm)/tau_s), 2);
        end
    end
end

% AEF RS Neuron
% equilibrium values from previous HW
Vdc = -69.9999;
Udc = 2.0000e-04;
% neuronal parameters
C = 200;
Gl = 10;
El = -70;
Vt = -50;
delt = 2;
a = 2;
tau = 30;
b = 0;
Vr = -58;
neuronType = "RS";

V = zeros(1,N);
U = zeros(1,N);

V(1) = Vdc;
U(1) = Udc; 
spikecount = 0;
for ind = 1:(N-1)
   dVdt = (-Gl*(V(ind)-El) + Gl*delt*exp((V(ind)-Vt)/delt) - U(ind)+I_app(ind))/C;   
   dUdt = (a*(V(ind)-El) - U(ind))/tau;

   V(ind+1) = V(ind) + dVdt*delT;
   U(ind+1) = U(ind) + dUdt*delT; 
   
   if V(ind+1) >= 0
       V(ind+1) = Vr;
       U(ind+1) = U(ind) + b;
       spikecount = spikecount + 1;
   end
end


figure;
plot(t, I_app); % t in ms, I in pA
title(join(['$I_{app}$ vs Time, Number of spikes = ', spikecount, ', $w_0 = $', w0, ", $\sigma_w = $", sigmaw]), 'interpreter', 'latex')
ylabel('Current in pA', 'interpreter', 'latex')
xlabel('Time in ms', 'interpreter', 'latex')

figure
plot(t, V)
title(join(['Voltage vs Time for neuron type ', neuronType, ', Number of spikes = ', spikecount, ', $w_0 = $', w0, ", $\sigma_w = $", sigmaw]), 'interpreter', 'latex')
ylabel('Voltage in mV', 'interpreter', 'latex')
xlabel('Time in ms', 'interpreter', 'latex')



