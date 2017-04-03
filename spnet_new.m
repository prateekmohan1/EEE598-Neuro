function [] = spnet_new()
  %50 Excitatory in both Layer 1 and Layer 2. Layer 2 also has 5 Inhibitory
  L1_Ne = 784;
  L2_Ne = 500;
  L3_Ne = 10;
  L2_Ni = 5;
  N_t = L1_Ne + L2_Ne + L3_Ne;

  %Number of Synapses in input to L1
  M_inp = L1_Ne;
  %Number of Synapses in L1 going out to L2 per neuron
  M_L1 = L2_Ne;
  %Number of Synapses in L2 going out to L3 per neuron
  M_L2 = L3_Ne;

  %Total Num synapses
  %syn_cnt = (L1_Ne + L2_Ne) * M;

  %Synapse connections
  synapses_inptoL1 = [ones(M_inp,1)];
  synapses_L1toL2 = [ones(L1_Ne,M_L1)];
  synapses_L2toL3 = [ones(L2_Ne,M_L2)];

  %Weights 
  weights_inptoL1 = [ones(L1_Ne,1)];
  weights_L1toL2 = [ones(L1_Ne,M_L1)];
  weights_L2toL3 = [ones(L2_Ne,M_L2)];

  %Membrane Potentials
  Vmem_L1 = [zeros(L1_Ne,1)];
  Vmem_L2 = [zeros(L2_Ne,1)];
  Vmem_L3 = [zeros(L3_Ne,1)];

  %Counters for tracking STDP
  cnt_max = 100; %max counter value for each neuron
  cnt_L1 = [zeros(L1_Ne,1)*cnt_max];
  cnt_L2 = [zeros(L1_Ne,1)*cnt_max];
  cnt_L3 = [zeros(L1_Ne,1)*cnt_max];
  %Tau used for exponential increase for LTP, the lower 
  %this is the larger the rate of increase for LTP
  exp_tau = 105; 

  %Decay Rates
  lambda_L1 = [ones(L1_Ne,1) * 0.02];
  lambda_L2 = [ones(L2_Ne,1) * 0.02];
  lambda_L3 = [ones(L3_Ne,1) * 0.02];

  %Neuron Config Details
  I = 1; %nA
  C = 1; %nF
  R = 40; %M Ohms
  Vth = 5;  %Spike Threshold

  %Refractory Periods for each neuron
  ref_period = 5;
  ref_cnt_L1 = [zeros(L1_Ne,1)];
  ref_cnt_L2 = [zeros(L2_Ne,1)];
  ref_cnt_L3 = [zeros(L3_Ne,1)];

  %Spikes
  spike_exists_L1 = [zeros(L1_Ne,1)];
  nxt_spike_exists_L2 = [zeros(L2_Ne,1)];
  spike_exists_L2 = [zeros(L2_Ne,1)];
  spike_exists_L3 = [zeros(L3_Ne,1)];
  nxt_spike_exists_L3 = [zeros(L3_Ne,1)];

  %Test input
  inp = [ones(L1_Ne,1)];

  for i=1:L1_Ne
    inp(i) = round(rand);
  end

  %Populate Connections from inp to L1
  for j = 1:M_inp
    synapses_inptoL1(j) = 1;
  end

  %Populate Connections from L1 to L2
  for i=1:L1_Ne
    for j=1:M_L1
      synapses_L1toL2(i,j) = 1;
    end;
  end;
    
  %Populate Connections from L2 to L1
  for i=1:L2_Ne
    for j=1:M_L2
      synapses_L2toL1(i,j) = 1;
    end;
  end;

  %Total Iterations
  for j=1:200

    %If this is the first iteration
    if (j == 1)
      %For all the input neurons, introduce spikes
      for i=1:L1_Ne
        [spike_exists_L1(i),Vmem_L1(i), ...
        ref_cnt_L1(i), cnt_L1(i)] = neuron(zeros(1,L1_Ne), zeros(1,M_L1), zeros(1,M_L1), Vth,
                                           Vmem_L1(i),ref_cnt_L1(i),lambda_L1(i), inp(i), 
                                           synapses_inptoL1(i), weights_inptoL1(i), cnt_L1(i),
                                           cnt_max);
      end
    else
      %At this point, you are not in your first iteration. So, spikes have been created in L1 due
      %to inputs inserted

      %%%%%%%%%%%
      %Calculate Membrane Potentials for L2 and L3
      %%%%%%%%%%%
      %Iterate over all neurons in layer 2
      for i=1:L2_Ne
        [nxt_spike_exists_L2(i), Vmem_L2(i), ...
        ref_cnt_L2(i), cnt_L2(i)] = neuron(spike_exists_L1, synapses_L1toL2(:,i)',
                                           weights_L1toL2(:,i)',Vth, Vmem_L2(i), 
                                           ref_cnt_L2(i), lambda_L2(i), inp(i), 0, 0,
                                           cnt_L2(i), cnt_max);
      end

      %Iterate over all neurons in layer 3
      for i=1:L3_Ne
        [nxt_spike_exists_L3(i),Vmem_L3(i), ...
        ref_cnt_L3(i), cnt_L3(i)] = neuron(spike_exists_L2, synapses_L2toL3(:,i)', 
                                           weights_L2toL3(:,i)', Vth, Vmem_L3(i), ref_cnt_L3(i),
                                           lambda_L3(i), inp(i), 0, 0, cnt_L3(i), cnt_max);
      end

      %%%%%%%%%%%%%%
      %Calculate the synapse weights
      %%%%%%%%%%%%%%
      %Iterate over all L2 neurons
      for i=1:L2_Ne
        %If a spike exists in L2
        if (nxt_spike_exists_L2(i))
          %Set the refractory rate 
          ref_cnt_L2(i) = ref_period;
          %Iterate over all rows in L1 which are the pre-syn (LTP)
          for k=1:L1_Ne
            %If a synapse connection exists between L1 and L2
            if (synapses_L1toL2(k,i) == 1)
              [weights_L1toL2(k,i)] = synapse_calc_pre(cnt_max, cnt_L1(k), weights_L1toL2(k,i), 
                                                      exp_tau);
            end
          end
          %Iterate over all rows in L3 which are the post-syn (LTD)
          for k=1:L3_Ne
            if (synapses_L2toL3(i,k) == 1)
              [weights_L2toL3(i,k)] = synapse_calc_post(cnt_max, cnt_L3(k), weights_L1toL2(i,k), 
                                                        exp_tau);
            end              
          end            
        end
      end

      %Iterate over all L1 neurons
      for i=1:L1_Ne
        %If a spike exists in L1
        if (spike_exists_L1(i))
          %Set the refractory rate
          ref_cnt_L1(i) = ref_period;
          %Iterate over all rows in L2 which are post-syn (LTD)
          for k=1:L2_Ne
            %If a synapse connection exists between L1 and L2
            if (synapses_L1toL2(i,k) == 1)
              [weights_L1toL2(i,k)] = synapse_calc_post(cnt_max, cnt_L2(k), weights_L1toL2(i,k), 
                                                        exp_tau);
            end
          end
        end
      end

      %Iterate over all L3 neurons
      for i=1:L3_Ne
        if (nxt_spike_exists_L3(i))
          %Set the refractory rate
          ref_cnt_L3(i) = ref_period;
          %Iterate over all rows in L2 which are pre-syn (LTP)
          for k=1:L2_Ne
            %If a synapse connection exists between L2 and L3
            if (synapses_L2toL3(k,i) == 1)
              [weights_L2toL3(k,i)] = synapse_calc_pre(cnt_max, cnt_L2(k), weights_L2toL3(k,i), 
                                                        exp_tau);
            end
          end              
        end
      end

      %%%%%%%%%%%
      %Introduce inputs again into the first layer
      %Calculate L1's membrane potentials
      %%%%%%%%%%%

      for i=1:L1_Ne
        [spike_exists_L1(i),Vmem_L1(i), ...
         ref_cnt_L1(i), cnt_L1(i)] = neuron(zeros(L1_Ne,1), zeros(1,M_L1), zeros(1,M_L1), Vth,
                                            Vmem_L1(i),ref_cnt_L1(i),lambda_L1(i), inp(i), 
                                            synapses_inptoL1(i), weights_inptoL1(i), cnt_L1(i), 
                                            cnt_max);
      end

      %%%%%%%%%%
      %Assign spikes for next iteration
      %%%%%%%%%%

      spike_exists_L2 = nxt_spike_exists_L2;
      spike_exists_L3 = nxt_spike_exists_L3;

    end

  end

  %plot(test2)



end

function [weights_post_out] = synapse_calc_post(counter_max, counter_post, weights_post, exp_tau)

  weights_post_out = weights_post - weights_post * exp(-(counter_max - counter_post)/exp_tau);

end

function [weights_pre_out] = synapse_calc_pre(counter_max, counter_pre, weights_pre, exp_tau)

  weights_pre_out = weights_pre + weights_pre * exp(-(counter_max - counter_pre)/exp_tau);

end


function [spike, Vout, ref, counter_out] = neuron (spike_exists, synapses, syn_weights, V_th, Vin, 
                                      ref_in, lambda, spike_exists_ext, synapses_ext, 
                                      syn_weights_ext, counter, counter_max)

  %Scaling factor for increasing membrane potentials
  sc_factor = 0.8;
                                        
  if ~ref_in
    %Vout = Vin - (Vin/(R*C)) + (I/C);
    temp = 0;
    temp2 = 0;
    for i = 1:size(synapses,2)
      %The spike_exists, synapses, and syn_weights are only for the layer that you are interested in
      %For example, if you are looking at a neuron in L2, the spike_exists contains if spike exists in neurons
      %in L1
      temp = temp + spike_exists(i)*synapses(i)*syn_weights(i)*sc_factor;
    end

    %This temp2 calculation is for the external stimulation (mainly for inputs)
    %The spike_exists, synapses, and syn_weights are only for the layer that you are interested in
    %For example, if you are looking at a neuron in L2, the spike_exists contains if spike exists in neurons
    %in L1
    temp2 = temp2 + spike_exists_ext*synapses_ext*syn_weights_ext*sc_factor;

    Vout = Vin + temp - lambda + temp2;
    ref = ref_in;
  else
    ref = ref_in - 1;
    Vout = 0.2*V_th; % reset voltage
  end
     
  if (Vout > V_th && Vout ~= 0.2*V_th)
    Vout = 50;  % emit spike
    spike = 1;
    counter_out = counter_max;
    % ref = abs_ref; % set refractory counter
  else
    spike = 0;
  end

  if (counter > 0)
    counter_out = counter - counter*exp(-1);
  else
    counter_out = counter;
  end

end