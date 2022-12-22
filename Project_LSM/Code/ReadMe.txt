1. preprocessing.m generates 77 channel input spikes data using TI46-IF.mat, which contains the raw speakers data, and stores it into preprocessing.mat.
2. connections.m generates the input and reservoir connectivity matrices and stores it in connections.mat.
3. RR.m generates the temporal reservoir response using preprocessing.mat and connections.mat and LIF neuron models.
4. TT.m performs 5-fold training and testing using the reservoir response and generates classification accuracy.