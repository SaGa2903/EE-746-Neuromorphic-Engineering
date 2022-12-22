close all
clear all
rng(4)
T = 500; % 500 ms
delT = 0.1; % 0.1ms
lambda = 1e-3; % poisson process rate (per milli-second)
p = lambda*delT; % expected number of spikes in delT length interval
t = 0:delT:(T-delT);
Ns = 100;
N = size(t, 2);

I0 = 1; % 1pA
w0 = 250;
sigmaw = 25;
tau = 15; % 15 ms
tau_s = tau/4;

pd = makedist('Normal',w0, sigmaw);
w = random(pd, 1, Ns);
gamma = 1;

spikecount = 1;
iters = 0;

spikecounts_list = [];
iters_list = [];
delwks = {};
deltks = {};
figure;
while(spikecount > 0) 
    pd = makedist('Binomial',2, p);
    spikes = random(pd, Ns, N);
    iters = iters + 1;
    I_app = zeros(1,N);
    tau = 15; % 15 ms
    tau_s = tau/4;
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
    tsp = zeros(1,N);
    
    for ind = 1:(N-1)
       dVdt = (-Gl*(V(ind)-El) + Gl*delt*exp((V(ind)-Vt)/delt) - U(ind)+I_app(ind))/C;   
       dUdt = (a*(V(ind)-El) - U(ind))/tau;

       V(ind+1) = V(ind) + dVdt*delT;
       U(ind+1) = U(ind) + dUdt*delT; 
       
       if V(ind+1) >= 0
           V(ind+1) = Vr;
           U(ind+1) = U(ind) + b;
           tsp(ind+1) = t(ind+1);
           spikecount = spikecount + 1;
       end
    end
    delwk = zeros(Ns, 1);
    deltk = zeros(Ns, 1);
    tspm = tsp(tsp>0);
    tau = 15; % 15 ms
    tau_s = tau/4;
    deltks{iters} = [];
    delwks{iters} = [];
    for i = 1: Ns
        deltkm = inf;
        delwktot = 0;
        for j = 1:size(tspm,2)
           tspikes = t(spikes(i,:)==1);
           prev_spikes= tspikes(tspikes<tspm(j));
           if ~isempty(prev_spikes)
               tktemp = max(prev_spikes);
               if(tspm(j)- tktemp<deltkm)
                   deltkm = tspm(j)- tktemp;
               end
           end
           if deltkm ~= inf
               deltks{iters} = [deltks{iters}, deltkm];
               delwk = -w(i)*gamma*(exp(-1*(deltkm)/tau) - exp(-1*(deltkm)/tau_s));
               delwktot = delwktot + delwk;
               delwks{iters} = [delwks{iters}, delwk];
           end  
        end
        w(i) = max(w(i) + delwktot, 10);
    end
   
    plot(t, V);
    spikecounts_list = [spikecounts_list, spikecount];   
    hold on;
end

title(join(['post synaptic potential vs time for multiple iterations']), 'interpreter', 'latex')
ylabel('$V$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
lgd = {};
for i = 1:iters
    lgd = [lgd, join(['iteration ',num2str(i),', ',num2str(spikecounts_list(i)),' spikes'])];
end
legend(lgd)
hold off;
% part a
disp('weights after training')
w

% part b
figure;
for iter = 1:iters
    hold on;
    scatter(deltks{iter}, delwks{iter});
end
lgd = {};
for i = 1:iters
    lgd = [lgd, num2str(i)];
end
legend(lgd)
title(join(['$\Delta w_k$ vs $\Delta t_k$ for multiple iterations']), 'interpreter', 'latex')
ylabel('$\Delta w_k$', 'interpreter', 'latex')
xlabel('$\Delta t_k$', 'interpreter', 'latex')




