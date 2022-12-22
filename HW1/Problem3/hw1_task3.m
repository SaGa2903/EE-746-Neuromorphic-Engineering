% Finding steady state values
syms V

%% RS
a = 2;%ns
Gl = 10;%ns
El = -70;%mV
Vt = -50;%mV
dt = 2; %mV
Vdc = -69.9999;


%% IB
% a = 4;%ns
% Gl = 18;%ns
% El = -58;%mV
% Vt = -50;%mV
% dt = 2; %mV
% Vdc = -57.9696;


%% CH
% a = 2;%ns
% Gl = 10;%ns
% El = -58;%mV
% Vt = -50;%mV
% dt = 2; %mV
% Vdc = -57.9690;



% Udc = a*(Vdc - El)
% Cdvdt = -Gl*(Vdc-El) + Gl*dt*exp((Vdc-Vt)/dt) -  a*(Vdc-El) % should equal 0
eqn1 = -Gl*(V-El) + Gl*dt*exp((V-Vt)/dt) -  a*(V-El) == 0;

S = solve(eqn1);

S