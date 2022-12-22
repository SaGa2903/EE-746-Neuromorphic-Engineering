close all
Fanout = {[0,0,0,0,0],[1,0,0,0,5],[1,0,0,0,5],[1,0,0,0,5],[0,0,0,0,0]};
Weight = {[0,0,0,0,0],[3000,0,0,0,3000],[3000,0,0,0,3000],[3000,0,0,0,3000],[0,0,0,0,0]};
Delay = {[0,0,0,0,0],[1e-3,0,0,0,8e-3],[5e-3,0,0,0,5e-3],[9e-3,0,0,0,1e-3],[0,0,0,0,0]};

N = 5;

spike_time=cell(1,N);
arrival_time = cell(1,N);
strength=cell(1,N);
pre_neuron=cell(1,N);

V = cell(1,N);
Iapp = cell(1,N);
Isyn = cell(1,N);

Iapp_prev_positive_transition = zeros(1,N);

C=300e-12; %pF
g_l = 30e-9; %nS
Vt = 20e-3; %mV
El = -70e-3; %mV
Rp = 2e-3; %ms

end_time = 50e-3; %ms
dt = 0.01e-3;  %ms

Iapp_start = {end_time+1,7e-3,3e-3,0,end_time+1};  %ms
I_injected = 50e-9;  %nA
Io = 1e-12;   %pA
pulse_duration = 1e-3;  %ms
tau = 15e-3;    %ms
tau_s = tau/4;   %ms

total_iterations = end_time/dt;

timer = zeros(1,total_iterations);



for i=1:N
    V{i} = El*ones(1,total_iterations);
    Iapp{i} = zeros(1,total_iterations);
    Isyn{i} = zeros(1,total_iterations);

    spike_time{i}=[];
    arrival_time{i}=[];
    strength{i}=[];
    pre_neuron{i}=[];

    Iapp_prev_positive_transition(i) = Iapp_start{i};
end

pre_neuron = {[2,3,4],[],[],[],[2,3,4]};
ind = 1;  % current time index

for t = 0:dt:(end_time-dt)
    for neuron = 1:1:5
        
        %Iapp calculation  % or make a function
%         if t>= Iapp_start{neuron}
% 
%             if (t-Iapp_prev_positive_transition(neuron)) < pulse_duration
%                 Iapp{neuron}(ind) = I_injected;
%             elseif (t-Iapp_prev_positive_transition(neuron)) >= pulse_duration && (t-Iapp_prev_positive_transition(neuron)) < 2*pulse_duration
%                 Iapp{neuron}(ind) = 0;
%             elseif (t-Iapp_prev_positive_transition(neuron)) == 2*pulse_duration
%                 Iapp_prev_positive_transition(neuron) = t;
%                 Iapp{neuron}(ind) = I_injected;
%             end
%         else
%             Iapp{neuron}(ind) = 0;
%         end
        
        if ((t >= Iapp_start{neuron})&&(t < (Iapp_start{neuron}+pulse_duration)))
           Iapp{neuron}(ind) = I_injected; 
        end

        %Iapp calc ends

        %Isyn coming in calc begins
        if isempty(pre_neuron{neuron}) == 0
            for pre = 2:4
                for each_spike = 1:length(spike_time{pre})
                    tk_tau = spike_time{pre}(each_spike) + Delay{pre}(neuron);
                    if t >= tk_tau
                        Isyn{neuron}(ind) = Isyn{neuron}(ind) + Io*Weight{pre}(neuron)*(exp(-(t-tk_tau)/tau) - exp(-(t-tk_tau)/tau_s));
                    end 

                end
            end
        end
        %Isyn coming in ends
        
            

        %V calc begins
        if isempty(spike_time{neuron}) == 1
                dv = ((-g_l*( V{neuron}(ind) - El)) + Iapp{neuron}(ind) + Isyn{neuron}(ind))*dt/C;
                V{neuron}(ind+1) = V{neuron}(ind) + dv;
                if V{neuron}(ind + 1) >= Vt
                    spike_time{neuron} = [spike_time{neuron},t];
                    V{neuron}(ind+1) = El;
                end               
        else
           if (t - spike_time{neuron}(end)) > Rp
                dv = ((-g_l*( V{neuron}(ind) - El)) + Iapp{neuron}(ind) + Isyn{neuron}(ind))*dt/C;
                V{neuron}(ind+1) = V{neuron}(ind) + dv;
                if V{neuron}(ind + 1) >= Vt
                    spike_time{neuron} = [spike_time{neuron},t];
                    V{neuron}(ind+1) = El;
                end
            end
        end
        %V calc ends

    end  % each neuron end

    timer(ind) = t;
    ind = ind+1;
    
  
end  % each time instant end


figure
for i = 1:N
    plot(1e3*timer(1:total_iterations),1e3*V{i}(1:total_iterations));
    title("case 2: Voltage vs Time for neuron")
    xlabel('Time (ms)') 
    ylabel('Voltage (mV)')
    hold on
end
legend({'1','2','3','4','5'})

hold off
figure
for i = 1:N
    plot(1e3*timer(1:total_iterations),1e12*(Iapp{i}(1:total_iterations)+Isyn{i}(1:total_iterations)));
    title("case 2: Current vs Time")
    xlabel('Time (ms)') 
    ylabel('Current (pA)')  
    hold on
end
legend({'1','2','3','4','5'})













