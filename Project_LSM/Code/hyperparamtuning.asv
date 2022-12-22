% hyperparameter tuning
% find the best alpha_Gin and alpha_Gres

% Input spike trains
%load preprocessing.mat;

% Network
load('connections.mat');
Nin = size(Gin,2);
Nres = size(Gin,1);

% Paramters
dt = 1e-3;

% Neuron model
Vth = 20;
tau = 64e-3;
Vrst = 0;
RefPer = 2e-3;
RPS = ceil(RefPer/dt);

% Synapse model 2nd order
tau1 = 8e-3;
tau2 = tau1/2;
I0 = 1/(tau1-tau2);
ds = 1e-3;
DS = ceil(ds/dt);
ts = 0:dt:ds+6*(tau1+tau2);
Iwave = I0*(exp(-(ts-ds)/tau1)-exp(-(ts-ds)/tau2));
Iwave(1:ceil(ds/dt)) = 0;

%% Reservoir

alpha_Gin_vals = [0.5, 1, 2, 4, 8, 16];
alpha_Gres_vals = [0.0625, 0.125, 0.25, 0.5, 1, 2];
error_matrix_train = zeros(size(alpha_Gin_vals, 2), size(alpha_Gres_vals, 2));
error_matrix_test = zeros(size(alpha_Gin_vals, 2), size(alpha_Gres_vals, 2));
RR_matrix = cell(size(alpha_Gin_vals, 2), size(alpha_Gres_vals, 2));

i_Gin = 0;
for alpha_Gin = alpha_Gin_vals
    i_Gin = i_Gin + 1;
    i_Gres = 0;
    for alpha_Gres = alpha_Gres_vals
        i_Gres = i_Gres + 1;
        Gnet = [alpha_Gin*Gin, alpha_Gres*Gres'];
        % Gnet = alpha_G*[Gin, Gres'];

        DATA(1).RES = [];

        for insample = 1:numel(DATA)
            n_samples = length(DATA(insample).S);
            in_spikes = DATA(insample).S;

            V = zeros(Nres,1);
            Ibuffer = zeros(Nin+Nres,length(Iwave));
            res_spikes = zeros(Nres,1);
            inres_spikes = zeros(Nin+Nres,1);
            Itotal = zeros(Nin+Nres,1);
            RP = zeros(Nres,1);
            Iin = zeros(Nres,1);

            Iin_all = zeros(Nres,n_samples);
            res_spikes_all = zeros(Nres,n_samples);
            Itotal_all = zeros(Nin+Nres,n_samples);

            for i = 1:n_samples
                inres_spikes = [in_spikes(:,i);res_spikes];
                Iin = Gnet*Itotal;

                res_spikes_all(:,i) = res_spikes;
                Itotal_all(:,i) = Itotal;
                Iin_all(:,i) = Iin;
                Ibuffer = circshift(Ibuffer,-1,2);
                Ibuffer(:,end) = 0;

                Ibuffer = Ibuffer + inres_spikes*Iwave;
                Itotal = Ibuffer(:,1);

                RP = RP-1; RP(RP<0) = 0;
                V = V*(1-dt/tau) + dt*Iin;
                V(RP>0) = 0;
                res_spikes = V>Vth; V(V<0) = Vrst; V(res_spikes) = Vrst; RP(res_spikes) = RPS;
            end  
            DATA(insample).RES = res_spikes_all;
        end

        % Perform TT on specific subset

        sset = 0:9; % mention digits for which to perform TT
        subset_indices = ismember([DATA.type], sset);
        DATA = DATA(subset_indices);

        Nres = 125;
        Nout = length(sset);

        % Paramters
        dt = 1e-3;

        % Neuron model
        Vth = 20;
        tau = 64e-3;
        Vrst = 0;
        RefPer = 2e-3;
        RPS = ceil(RefPer/dt);

        % Synapse model 2nd order
        tau1 = 8e-3;
        tau2 = tau1/2;
        I0 = 1/(tau1-tau2);
        ds = 1e-3;
        DS = ceil(ds/dt);
        ts = 0:dt:ds+6*(tau1+tau2);
        Iwave = I0*(exp(-(ts-ds)/tau1)-exp(-(ts-ds)/tau2));
        Iwave(1:ceil(ds/dt)) = 0;

        rZ=[];Ro=[];
        for sample_i=1:numel(DATA)
            rZ(:,sample_i)=sum(DATA(sample_i).RES,2);
            index_here = find(DATA(sample_i).type == sset);
            Ro(index_here,sample_i)=1;
        end

        avg_spikes_digit = mean(sum(rZ,1));

        Nfold = 5;
        Nstride = floor(numel(DATA)/Nfold);
        Wlim = 8;

        RESULT = struct([]);

        for kfold = 1:Nfold

            testOn = (kfold-1)*Nstride+(1:Nstride);
            trainOn  = setdiff(1:numel(DATA),testOn);

            ZT=rZ(:,trainOn)';RoT=Ro(:,trainOn)';
            W0=[];  
            for i =1:Nout
                options = optimset('Display', 'notify');
                W0(:,i)=lsqlin(ZT,1000*RoT(:,i),[],[],[],[],-Wlim*ones(Nres,1),Wlim*ones(Nres,1),[], options);
            end
            W=W0';

        %     Mn = sparse(1:length(trainOn),1+[DATA(trainOn).type],1,length(trainOn),Nout); % confusion matrix finder;
            Mn = RoT;
        %     MTest = sparse(1:length(testOn),1+[DATA(testOn).type],1,length(testOn),Nout); % confusion matrix finder;
            MTest = Ro(:,testOn)';
            %---------------

            spikeSampleCountTest=zeros(Nout,numel(DATA));

            for insample = 1:numel(DATA)

                n_samples = length(DATA(insample).RES);
                in_spikes = DATA(insample).RES;

                V = zeros(Nout,1);
                Ibuffer = zeros(Nres,length(Iwave));
                out_spikes = zeros(Nout,1);
                Itotal = zeros(Nres,1);
                RP = zeros(Nout,1);
                Iin = zeros(Nout,1);

                Iin_all = zeros(Nout,n_samples);
                Itotal_all = zeros(Nres,n_samples);
                out_spikes_all = zeros(Nout,n_samples);

                for i = 1:n_samples
                    Iin = W*Itotal;

                    Iin_all(:,i) = W*Itotal;
                    Itotal_all(:,i) = Itotal;
                    out_spikes_all(:,i) = out_spikes;

                    Ibuffer = circshift(Ibuffer,-1,2);
                    Ibuffer(:,end) = 0;
                    Ibuffer = Ibuffer + in_spikes(:,i)*Iwave;
                    Itotal = Ibuffer(:,1);

                    RP = RP-1; RP(RP<0) = 0;
                    V = V*(1-dt/tau) + dt*Iin;
                    V(RP>0) = 0;
                    out_spikes = V>Vth; V(V<0) = Vrst; V(out_spikes) = Vrst; RP(out_spikes) = RPS;

                end  
                spikeSampleCountTest(:,insample) = sum(out_spikes_all,2);
            end

            traincounts = spikeSampleCountTest(:,trainOn);
            testcounts = spikeSampleCountTest(:,testOn);

            spikeSampleCountTest = [traincounts, testcounts];

            [M,recognized] = max(spikeSampleCountTest);
            Y = spikeSampleCountTest./repmat(M,Nout,1);

            M  = blkdiag(Mn,MTest);
            Y = Y./repmat(max(Y),Nout,1); Y(Y~=1)=0;

            misClassifiedSamples = find(sum(Y,1)~=1);
            Y(:,misClassifiedSamples) = 0;
            CM = (Y*M); % confusion matrix

            accuracyTrain = 100*trace(CM(:,1:Nout))/length(trainOn);
            numCorrectTrain = trace(CM(:,1:Nout));
            accuracyTest = 100*trace(CM(:,Nout+(1:Nout)))/length(testOn);
            numCorrectTest = trace(CM(:,Nout+(1:Nout)));

            RESULT(kfold).accTest = accuracyTest;
            RESULT(kfold).accTrain = accuracyTrain;
            RESULT(kfold).numCorrectTrain = numCorrectTrain;
            RESULT(kfold).numCorrectTest = numCorrectTest;

            s=sprintf('\r\n Accuracy : kFold(%i) Test %2.2f (%i/%i) Train:%2.2f (%i/%i) \n',kfold,RESULT(kfold).accTest,RESULT(kfold).numCorrectTest,length(testOn),RESULT(kfold).accTrain,RESULT(kfold).numCorrectTrain,length(trainOn));
            fprintf(s);
        end

        disp([RESULT.accTest])
        error_test = 100 - mean([RESULT.accTest]);
        error_matrix_test(i_Gin, i_Gres) = error_test;

        error_train = 100 - mean([RESULT.accTrain]);
        error_matrix_train(i_Gin, i_Gres) = error_train;
        RR_matrix{i_Gin, i_Gres} = DATA;
        
%         save('error_matrix_test.mat', 'error_matrix_test')
%         save('error_matrix_train.mat', 'error_matrix_train')
        %save('RR_matrix.mat', 'RR_matrix')
    end
end