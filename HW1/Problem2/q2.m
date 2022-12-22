Izhikevich_simulate(3,"RS",[400;500;600],0.1,500);
Izhikevich_simulate(3,"IB",[400;500;600],0.1,500);
Izhikevich_simulate(3,"CH",[400;500;600],0.1,500);
function [U,V] = Izhikevich_simulate(N,neuronType, Iapp, dt, T)
    if neuronType=="RS"
%         C=100e-12;
%         kz=0.7e-6;
%         Er=-60e-3;
%         Et=-40e-3;
%         a=0.03e3;
%         b=-2e-9;
%         c=-50e-3;
%         d=100e-12;
%         vpeak=35e-3;
        C=100;
        kz=0.7;
        Er=-60;
        Et=-40;
        a=0.03;
        b=-2;
        c=-50;
        d=100;
        vpeak=35;
    elseif neuronType=="IB"
%         C=150e-12;
%         kz=1.2e-6;
%         Er=-75e-3;
%         Et=-45e-3;
%         a=0.01e3;
%         b=5e-9;
%         c=-56e-3;
%         d=130e-12;
%         vpeak=50e-3;
        C=150;
        kz=1.2;
        Er=-75;
        Et=-45;
        a=0.01;
        b=5;
        c=-56;
        d=130;
        vpeak=50;       
    elseif neuronType=="CH"
%         C=50e-12;
%         kz=1.5e-6;
%         Er=-60e-3;
%         Et=-40e-3;
%         a=0.03e3;
%         b=1e-9;
%         c=-40e-3;
%         d=150e-12;
%         vpeak=25e-3;  
        C=50;
        kz=1.5;
        Er=-60;
        Et=-40;
        a=0.03;
        b=1;
        c=-40;
        d=150;
        vpeak=25;
    else
        disp("INVALID TYPE");
        return;
    end
    
    % unit conversions
%     C=C*1e-12;
%     kz=kz*1e-6;
%     Er=Er*1e-3;
%     Et=Et*1e-3;
%     a=a*1e3;
%     b=b*1e-9;
%     c=c*1e-3;
%     d=d*1e-12;
%     vpeak=vpeak*1e-3;
%     Iapp = Iapp*1e-12;
%     T = T*1e-3;
%     dt = dt*1e-3;
    %
    
    M = ceil(T/dt);
    t = 0:dt:dt*(M-1);
    %Taking first solution for steady state
    U=((b/kz) + Et)*ones(N,M);
    V=(b*((b/kz) + Et - Er))*ones(N,M);
    
    function der = dvdt(v, u)
        der = (kz*(v-Er).*(v-Et) - u + Iapp)/C;
    end

    function der = dudt(u, v)
        der = a*(b*(v-Er) - u);
    end
    
    for i=1:M-1
%         arr = V(:,i);
%         arr(arr>=vpeak)= c;
%         arr1= U(:,i);
%         arr1(arr>=vpeak) = arr1(arr>=vpeak)+d;
%         V(:,i) = arr;
%         U(:,i) = arr1;

%         V(:,i) = minimum(V(:,i), vpeak);
        %4th order Runge Kutta Method
        
        k1 = dvdt(V(:, i), U(:, i));
        l1 = dudt(U(:, i), V(:, i));
        
        k2 = dvdt(V(:,i) + dt*k1/2, U(:,i)+dt*l1/2);
        l2 = dudt(U(:,i) + dt*l1/2, V(:,i)+dt*k1/2);
        
        k3 = dvdt(V(:,i) + dt*k2/2, U(:,i)+dt*l2/2);
        l3 = dudt(U(:,i) + dt*l2/2, V(:,i)+dt*k2/2);        
        
        k4 = dvdt(V(:,i) + dt*k3, U(:,i)+dt*l3);
        l4 = dudt(U(:,i) + dt*l3, V(:,i)+dt*k3);
        
        V(:,i+1) = V(:,i) + dt*(k1+2*k2+2*k3+k4)/6;
        U(:,i+1) = U(:,i) + dt*(l1+2*l2+2*l3+l4)/6;
        
        for j = 1:N
            if V(j,i+1)>=vpeak
                V(j,i+1) = c;
                U(j,i+1) = U(j,i+1) + d;
            end
        end
    end
    
%     arr = V(:,M);
%     arr(arr>vpeak)= c;
%     arr1= U(:,M);
%     arr1(arr>vpeak) = arr1(arr>vpeak)+d;
%     V(:,M) = arr;
%     U(:,M) = arr1;
    for i= 1:N
        figure
        plot(t,V(i,:));
        title(join(['Membrane Potential vs Time for ', neuronType, ', I_{app} = ',Iapp(i),' pA']))
        xlabel('Time in ms', 'interpreter', 'latex')
        ylabel('Voltage in mV', 'interpreter', 'latex')
%         figure
%         plot(t*1e3,U(i,:)*1e12);
    end
end

    
    
        
        
    
    
        