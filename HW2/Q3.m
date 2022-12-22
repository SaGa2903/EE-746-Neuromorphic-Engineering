rng(22)
T = 500; % 500 ms
delT = 0.1; % 0.1ms
lambda = 1e-3; % poisson process rate (per milli-second)
p = lambda*delT; % expected number of spikes in delT length interval
t = 0:delT:(T-delT);
Ns = 100;
N = size(t, 2);

I0 = 1; % 1pA
w0 = 50;
sigmaw = 5;

pd = makedist('Normal',w0, sigmaw);
w = random(pd, 1, Ns);
gamma = 1;

spikecount = 0;
iters = 0;

delwks = [];
deltks = [];

while(spikecount < 1)
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
    
    % store state at peak of membrane potential
    Vmax = Vdc;
    tmax = 0;
    imax = 0;
    
    for ind = 1:(N-1)
       dVdt = (-Gl*(V(ind)-El) + Gl*delt*exp((V(ind)-Vt)/delt) - U(ind)+I_app(ind))/C;   
       dUdt = (a*(V(ind)-El) - U(ind))/tau;

       V(ind+1) = V(ind) + dVdt*delT;
       U(ind+1) = U(ind) + dUdt*delT; 
       
       if V(ind+1) > Vmax
          Vmax = V(ind + 1);
          tmax = t(ind + 1);
          imax = ind + 1;
       end
       
       if V(ind+1) >= 0
           V(ind+1) = Vr;
           U(ind+1) = U(ind) + b;
           spikecount = spikecount + 1;
       end
    end
    
    spikesprev = spikes(:, 1:(imax-1));
    delwk = zeros(Ns, 1);
    deltk = zeros(Ns, 1);
    tau = 15; % 15 ms
    tau_s = tau/4;
    for i = 1:Ns
       tspikes = t(spikesprev(i,:)==1);
       if isempty(tspikes)
           delwk(i) = 0;
           deltk(i) = 0;
       else
           tk = max(tspikes);
           deltk1 = tmax - tk;
           deltk(i) = deltk1;
           delwk(i) = w(i)*gamma*(exp(-1*(deltk1)/tau) - exp(-1*(deltk1)/tau_s));
       end  
       w(i) = min(w(i) + delwk(i), 500);
    end

    deltks = horzcat(deltks, deltk);
    delwks = horzcat(delwks, delwk);
end

% part a
disp('weights after training')
w

% part b
figure;
for iter = 1:iters
    hold on;
    scatter(deltks(:,iter), delwks(:, iter));
end
lgd = {};
for i = 1:iters
    lgd = [lgd, num2str(i)];
end
legend(lgd)
title(join(['$\Delta w_k$ vs $\Delta t_k$ for multiple iterations']), 'interpreter', 'latex')
ylabel('$\Delta w_k$', 'interpreter', 'latex')
xlabel('$\Delta t_k$', 'interpreter', 'latex')




