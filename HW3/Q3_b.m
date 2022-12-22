clear all;
close all;
N = 200;

%Thorughout, the neurons are being addressed as indices!!

% Also, note that in the V graph, as soon as we reach Vt, I coming back to
% El and not showing the spike to peak

%params
We = 1000;
Ws = 3000;
inh_delay = 1; %ms
Wi = -We;
Io = 1; %pA

weight = cell(1,N);   % each cell is a single number
delay = cell(1,N);    % each cell is a single number
fanout = cell(1,N);   % each cell is an array of 50 numbers sorted  in ascending order

exc_delay = randi([1, 20], 0.8*N, 1); % delays in ms
delay(1:0.8*N) = num2cell(exc_delay);  % putting delays in first 80% exc neurons

for i = 1:0.8*N
    fanout{i} = sort(randperm(N,N/10));    % every neuron has 50 connections and decided by random and unique
    weight{i} = We;   % exc 
end

for i = (0.8*N+1):N
    fanout{i} = sort(randperm(0.8*N,N/10));    % every neuron has 50 connections and decided by random and unique
    weight{i} = Wi;  %inhibitory
    delay{i} = inh_delay;
end

spike_time=cell(1,N);
arrival_time = cell(1,N);
pre_neuron=cell(1,N);

% filling preneruon cell array
for i=1:N %each neuron
    for j=1:length(fanout{i})
        pre_neuron{fanout{i}(j)} = [pre_neuron{fanout{i}(j)},i];
    end
end

V = cell(1,N);
Iapp = cell(1,N);
Isyn = cell(1,N);

C=300; %pF
g_l = 30; %nS
Vt = 20; %mV
El = -70; %mV
Rp = 2; %ms
tau = 15;    %ms
tau_s = 3.75;   %ms

end_time = 300; %ms
dt = 0.1;  %ms

%poisson stimulus for first 25 neurons
lambda = 0.1; % poisson process rate (per milli-second)
p = lambda*dt; % expected number of spikes in delT length interval
time_for_poisson = 0:dt:(end_time-dt);
N_p = size(time_for_poisson, 2);
pd = makedist('Binomial',2, p);
stim_times = cell(1,25);

for i=1:25
    spikes_temp = random(pd, 1, N_p);
    stim_times{i} = find(spikes_temp==1);
end

total_iterations = end_time/dt;
timer = zeros(1,total_iterations);

for i=1:N  % for each neuron
    V{i} = -70*ones(1,total_iterations); 
    Isyn{i} = zeros(1,total_iterations);
    
    if i<=25
        Iapp{i} = zeros(1,total_iterations);
    end

    spike_time{i}=[];
    arrival_time{i}=[];
end





ind = 1;  % current time index
for t = 0:dt:(end_time-dt)
  
    for neuron = 1:N
        
        %Iapp calculation
        if neuron<=25
            % apply poisson stimulus here
            stim_done = stim_times{neuron}(stim_times{neuron} < ind);
            Iapp{neuron}(ind) = Io*Ws*sum(  exp(-(t - stim_done*dt)/tau) - exp(-(t - stim_done*dt)/tau_s));
           
            stim_done = [];


        else
            Iapp{neuron}(ind) = 0;
        end
        %Iapp calc ends


        %Isyn coming in calc begins
        if ~isempty(pre_neuron{neuron})
            for pre = pre_neuron{neuron}
                spike_times_temp = spike_time{pre}(spike_time{pre}<=(t-delay{pre})) + delay{pre};
                Isyn{neuron}(ind) = Isyn{neuron}(ind) + Io*weight{pre}*sum(  exp(-(t-spike_times_temp)/tau) - exp(-(t- spike_times_temp)/tau_s)   );
            end 
        end
        %Isyn coming in ends
        
            

        %V calc begins
        if isempty(spike_time{neuron}) == 1

            if V{neuron}(ind) >= Vt
                spike_time{neuron} = [spike_time{neuron},t];
                V{neuron}(ind+1) = El;
            else
                dv = ((-g_l*( V{neuron}(ind) - El)) + Iapp{neuron}(ind) + Isyn{neuron}(ind))*dt/C;
                V{neuron}(ind+1) = V{neuron}(ind) + dv;
            end

        else

            if V{neuron}(ind) >= Vt
                spike_time{neuron} = [spike_time{neuron},t];
                V{neuron}(ind+1) = El;
            elseif (t - spike_time{neuron}(end)) < Rp
                V{neuron}(ind+1) = El;
            else
                dv = ((-g_l*( V{neuron}(ind) - El)) + Iapp{neuron}(ind) + Isyn{neuron}(ind))*dt/C;
                V{neuron}(ind+1) = V{neuron}(ind) + dv;
            end
        end
        %V calc ends

    end  % each neuron end

    timer(ind) = t;
    ind = ind+1;
    
end  % each time instant end

Re = zeros(total_iterations - 10/dt, 1);
Ri = zeros(total_iterations - 10/dt, 1);
time_minus_10 = 0:dt:((end_time - dt) - 10);
iter = 1;
for tk = time_minus_10
    Re(iter) = 0;
    for n = 1:0.8*N
       Re(iter) = Re(iter) + nnz((spike_time{n} >= tk) & (spike_time{n} <= tk+10)); 
    end
    
    Ri(iter) = 0;
    for n = (0.8*N+1):N
       Ri(iter) = Ri(iter) + nnz((spike_time{n} >= tk) & (spike_time{n} <= tk+10)); 
    end
    
    iter = iter + 1;
end

figure
hold on;
plot(time_minus_10, Ri)
plot(time_minus_10, Re)
hold off;
legend({'$R_e(t)$', '$R_i(t)$'}, 'Interpreter',"latex");
xlabel('Time');

spike_times_concat = [];
for i = 1:N
    spike_times_concat = [spike_times_concat, (i-1)*end_time + spike_time{i}];
end
rasterplot(spike_times_concat/dt, N, total_iterations, dt*1e-3)

% figure 
% plot(timer(1:total_iterations),V{1}(1:total_iterations));
% title("kjnx")
% xlabel('Time') 
% ylabel('V{1}') 
% 
% figure 
% plot(timer(1:total_iterations),Isyn{1}(1:total_iterations));
% title("kjnx")
% xlabel('Time') 
% ylabel('Isyn{1}')
% 
% figure 
% plot(timer(1:total_iterations),Iapp{1}(1:total_iterations));
% title("kjnx")
% xlabel('Time') 
% ylabel('Iapp{1}') 
