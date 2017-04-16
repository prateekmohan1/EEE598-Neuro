function [] = neuron_new_test()

    Vmem_L2 = 0;
    myPoissonSpikeTrain = rand(1, 784) < 100*0.001;
    synapses_L1toL2 = ones(L1_Ne,1);
    weights_L1toL2 = ones(L1_Ne,1);
    
    ref_cnt_L2 = 0;
    lambda_L2 = 0;

    [Vmem_L2] = neuron_new(myPoissonSpikeTrain, synapses_L1toL2', ...
                                               weights_L1toL2',Vth, Vmem_L2, ...
                                               ref_cnt_L2, lambda_L2);     



end