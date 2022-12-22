for neuron = 1:3
    for Iapp = [250, 350, 450]
        if neuron == 1
            Vdc = -69.9999;
            Udc = 2.0000e-04;
        elseif neuron == 2
            Vdc = -57.9696;
            Udc = 0.1216;               
        elseif neuron == 3
            Vdc = -57.9690;
            Udc = 0.0620;
        end
        plotter(neuron, Iapp,Vdc,Udc)
    end
end

function plotter(type,Iapp,Vdc,Udc)
    
    if type == 1
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

    elseif type == 2
        C = 130;
        Gl = 18;
        El = -58;
        Vt = -50;
        delt = 2;
        a = 4;
        tau = 150;
        b = 120;
        Vr = -50;
        neuronType = "IB";
    else 
        C = 200;
        Gl = 10;
        El = -58;
        Vt = -50;
        delt = 2;
        a = 2;
        tau = 120;
        b = 100;
        Vr = -46;
        neuronType = "CH";
    end

    dt = 0.1;
    h =dt;
    
    V = zeros(1,5000);
    U = zeros(1,5000);
    time = zeros(1,5000);

    V(1) = Vdc;   %%%%%%%%%%put correct value
    U(1) = Udc;    %%%%%%%%%%%%%put correct value
    time(1) = 0;

    for ind = 1:5000

        diffV = (-Gl*(V(ind)-El) + Gl*delt*exp((V(ind)-Vt)/delt) - U(ind)+Iapp)/C;
        
        diffU = (a*(V(ind)-El) - U(ind))/tau;

        V(ind+1) = V(ind) + diffV*h;
        U(ind+1) = U(ind) + diffU*h;
        if V(ind+1) >= 0
            V(ind+1) = Vr;
            U(ind+1) = U(ind) + b;
        end
        
        time(ind+1) = time(ind) + dt;
        
    end
    
    
    figure
    plot(time,V)
    title(join(['Voltage vs Time for neuron type ', neuronType, ", $I_{app} = $", Iapp, " pA"]), 'interpreter', 'latex')
    xlabel('Voltage in mV', 'interpreter', 'latex')
    ylabel('Time in ms', 'interpreter', 'latex')
end


    
    







        




     

     
