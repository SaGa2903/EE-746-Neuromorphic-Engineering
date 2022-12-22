syms V

% in units of mS/cm^2
Gna = 120;
Gk = 36;
Gl = 0.3;

% in mV
Ena = 50;
Ek = -77;
El = -55;

V = -65.1560308;

am = 0.1*(V+40)/(1 - exp(-(V+40)/10) ) ;
bm = 4 * exp(-0.0556*(V + 65));

an = 0.01*(V + 55)/(1 - exp(-(V + 55)/10));
bn = 0.125*exp(-(V+ 65)/80);

ah = 0.07* exp(-0.05*(V + 65));
bh = 1/(1 + exp(-0.1*(V+ 35)));

m = am/(am+bm);
n = an/(an+bn);
h = ah/(ah+bh);

Ina = Gna * (m^3) *h*(V - Ena);
Ik = Gk*(n^4) * (V-Ek);
Il = Gl * (V-El);
Iext = 0;

%eqn1 = Iext - (Ina + Ik + Il) == 0;

%S = solve(eqn1);

%S

m
n
h
