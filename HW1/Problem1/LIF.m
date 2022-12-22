% Problem 1
% LIF Model

p = getModelConstants();
C = p.C;
gL = p.gL;
VT = p.VT;
EL = p.EL;
I_c = p.I_c;

s = getSimConstants();
alpha = s.alpha;
delT = s.delT;
T = s.T;
M = s.M;
t_arr = s.t_arr;
N = s.N;
k_arr = s.k_arr;

I = zeros(N, M);
for k = 1:N
   I(k,:) = (1+k*alpha)*I_c; 
end

[V, avg_t_spike] = LIF_RK2(I);
 
for k = k_arr
    figure
    plot(t_arr*1e3, 1e3*V(k, :)')
    title(['Voltage of neuron ' num2str(k) ' vs Time'])
    xlabel('time in ms', 'interpreter', 'latex')
    ylabel('voltage in mV', 'interpreter', 'latex')    
end

figure
plot(1e9*I(:,1), 1e6*avg_t_spike)
title(['Average Spiking Time vs Current'])
xlabel('Current in nA', 'interpreter', 'latex')
ylabel('Time in $\mu s$', 'interpreter', 'latex')

function consts = getModelConstants()
    % set constants for the model
    % returns a struct
    C = 300e-12; % 300 pF
    gL = 30e-9; % 30 nS
    VT = 20e-3; % 20mV
    EL = -70e-3; % -70mV
    I_c = gL*(VT-EL); % minimum value of steady state current for spiking
    
    consts.C = C;
    consts.gL = gL;
    consts.VT = VT;
    consts.EL = EL;
    consts.I_c = I_c;
end 

function consts = getSimConstants()
    % set constants for simulation
    % returns a struct
    alpha = 0.1;
    delT = 0.1e-3; % 0.1ms
    T = 500e-3; % 500ms
    M = T/delT;
    t_arr = 0:delT:(T-delT);
    N = 10; % number of neurons
    k_arr = [2,4,6,8]; % indices of neurons to plot
    
    consts.alpha = alpha;
    consts.delT = delT;
    consts.T = T;
    consts.M = M;
    consts.t_arr = t_arr;
    consts.N = N;
    consts.k_arr = k_arr;
end

function der = dvdt(v, i)
    % returns dv/dt = f(t, v)
    % i(t) is the only part dependent on t, thus we pass it directly
    p = getModelConstants();
    C = p.C;
    gL = p.gL;
    VT = p.VT;
    EL = p.EL;
    I_c = p.I_c;

    der = (-gL/C)*v +((gL*EL)/C) + i/C;
end

function [V, avg_t_spike] = LIF_RK2(I)
    % input: NxM matrix of current values
    % output: NxM matrix of voltage values
    % using 2nd order Runge Kutta to find voltage
    p = getModelConstants();
    C = p.C;
    gL = p.gL;
    VT = p.VT;
    EL = p.EL;
    I_c = p.I_c;
    
    s = getSimConstants();
    alpha = s.alpha;
    delT = s.delT;
    T = s.T;
    M = s.M;
    t_arr = s.t_arr;
    N = s.N;
    k_arr = s.k_arr;
    
    [N, M] = size(I);
    V = zeros(N,M);
    % calculate resting potential V_resting
    V_resting = EL;
    V(:,1) = V_resting;
    
    % store the time instants of spiking
    t_spikes = cell(N, 1);
    for n = 1:N
       t_spikes{n} = []; 
    end
    
    for m = 2:M
        % simultaneously update all voltage values using Runge-Kutta
        V(:, m) = V(:, m-1) + delT*dvdt(V(:, m-1) + 0.5*delT*dvdt(V(:, m-1), I(:, m-1)), 0.5*(I(:, m-1) + I(:, m))); % assuming current at t + delT/2 can be approximated by the average of i(t) and i(t + delT)
        % check if spiking condition is met for each neuron
        for i = 1:N
           if (V(i, m) > VT) % fire if threshold crossed
               V(i, m) = EL;
               t_spikes{i} = [t_spikes{i}, m]; % store indices of spike times
           end
        end
    end
    
    % compute average time between spikes for each neuron
    avg_t_spike = zeros(N, 1);
    for n = 1:N
        L = length(t_spikes{n});
        avg_t_spike(n) = (1/(L-1))*sum(t_spikes{n}(2:L)) - (1/(L-1))*sum(t_spikes{n}(1:L-1));
        avg_t_spike(n) = delT*avg_t_spike(n);
    end
end


