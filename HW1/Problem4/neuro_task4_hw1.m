% in units of mS/cm^2
Gna = 120;
Gk = 36;
Gl = 0.3;

% in mV
Ena = 50;
Ek = -77;
El = -55;

% in uF/cm^2
C = 1;

% in ms
dt = 0.01;

%in uA/cm^2
Io = 15; 
Iext = 0; 

%initialize
m =0.0520;
n =0.3153;
h =0.6016;

Vdc = -65.1560308; % found using another matlab program for dc resting state
V = Vdc;  % initializing V

time = zeros(1,9000);
Na_curr = zeros(1,9000);
K_curr = zeros(1,9000);
L_curr = zeros(1,9000);
Voltage = zeros(1,9000);
Ext_curr = zeros(1,9000);

Pna = zeros(1,9000);
Pk = zeros(1,9000);
Pl = zeros(1,9000);
P_cap = zeros(1,9000);
diss_ener = zeros(1,9000); 

ind = 1;

for t = 30:dt:120
    if t == 60; Iext = 2*Io;end
    if t == 90; Iext = 0; end
    
    am = 0.1*(V+40)/(1 - exp(-(V+40)/10) ) ;
    bm = 4 * exp(-0.0556*(V + 65));
    
    an = 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
    bn = 0.125*exp(-(V+ 65)/80);
    
    ah = 0.07* exp(-0.05*(V + 65));
    bh = 1/(1 + exp(-0.1*(V+ 35)));

    dm = (am*(1-m) - bm*m)*dt;
    m = m+dm;
    dn = (an*(1-n) - bn*n)*dt;
    n = n+dn;
    dh = (ah*(1-h) - bh*h)*dt;
    h = h+dh;

    Ina = Gna * (m^3) *h*(V - Ena);
    Ik = Gk*(n^4) * (V-Ek);
    Il = Gl * (V-El);


    dV = (Iext - Ina-Ik-Il)*dt/C;
    V = V+dV;
    
    time(ind) = t;
    Na_curr(ind) = Ina;
    K_curr(ind) = Ik;
    L_curr(ind) = Il;
    Voltage(ind) = V;
    Ext_curr(ind) = Iext;
    
    Pna(ind) = Ina*(V-Ena);
    Pk(ind) = Ik*(V-Ek);
    Pl(ind) = Il*(V-El);
    P_cap(ind) = V*(Iext-Ina-Ik-Il);



    if ind>1
        diss_ener(ind) = diss_ener(ind-1)+(P_cap(ind)+Pna(ind)+Pk(ind)+Pl(ind))*dt;
    end
    

    ind = ind+1;
end

figure 
plot(time,Ext_curr);
title("External Current")
xlabel('Time (in ms)') 
ylabel('I applied (in nA)') 


figure
plot(time,Voltage);
title("Spikes in Membrane Potential (Hodgkin-Huxley model)")
xlabel('Time (in ms)') 
ylabel('Membrane Potential (in mV)') 

figure
plot(time,Na_curr,time,-K_curr,time,L_curr);
title("Ion currents vs Time")
legend({'Na','K','l'})
xlabel('Time (in ms)') 
ylabel('Ion Current')

figure
plot(time,P_cap,time,Pna,time,Pk,time,Pl)
title("Instantaneous Power dissipated")
legend({'Cap','Na','K','l'})
xlabel('Time (in ms)') 
ylabel('Instantaneous Power')

figure
plot(time,diss_ener)
title("Total Energy Dissipated vs Time")
xlabel('Time (in ms)') 
ylabel('Total Energy dissipated')


