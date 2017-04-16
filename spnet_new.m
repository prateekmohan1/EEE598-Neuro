function [weights_L1toL2, weights_L2toL3, freqout, sum_L3_out] = spnet_new(mode, wL1_L2, wL2_L3)

  %Mode = 0 means classification, Mode = 1 means training

  %50 Excitatory in both Layer 1 and Layer 2. Layer 2 also has 5 Inhibitory
  L1_Ne = 784;
  L2_Ne = 500;
  L3_Ne = 10;
  L2_Ni = 5;
  N_t = L1_Ne + L2_Ne + L3_Ne;
  
  %Load the MNIST database
  load mnist_uint8;
  test_x  = double(test_x)  / 255;
  test_y  = double(test_y);

  %Number of Synapses in input to L1
  M_inp = L1_Ne;
  %Number of Synapses in L1 going out to L2 per neuron
  M_L1 = L2_Ne;
  %Number of Synapses in L2 going out to L3 per neuron
  M_L2 = L3_Ne;
  
  %Frequencies for test
  freqout_row = [];
  freqout = [];
  
  %Total Num synapses
  %syn_cnt = (L1_Ne + L2_Ne) * M;
  
  %Sum of spikes
  sum_L1 = [];
  sum_L2 = [];
  sum_L3 = [];

  sum_L1_out = [];
  sum_L2_out = [];
  sum_L3_out = [];
  
  %Synapse connections
  synapses_inptoL1 = [ones(M_inp,1)];
  synapses_L1toL2 = [ones(L1_Ne,M_L1)];
  %synapses_L1toL2 = randi([0 1], L1_Ne, M_L1);
  synapses_L2toL3 = [ones(L2_Ne,M_L2)];

  %Weights 
  weights_inptoL1 = [ones(L1_Ne,1)];
  %weights_L1toL2 = 0.2 + (0.6-0.2)*rand(L1_Ne, M_L1);
  weights_L1toL2 = [ones(L1_Ne,M_L1)];
  weights_L2toL3 = 0.2 + (0.6-0.2)*rand(L2_Ne, M_L2);
  
  %Unnormalized weights
  unnweights_L1toL2 = [ones(L1_Ne,M_L1)];
  unnweights_L2toL3 = [ones(L2_Ne,M_L2)];
  
  %Membrane Potentials
  Vmem_L1 = [zeros(L1_Ne,1)*(-65)];
  Vmem_L2 = [zeros(L2_Ne,1)*(-65)];
  Vmem_L3 = [zeros(L3_Ne,1)*(-65)];

  %Counters for tracking STDP
  cnt_max = 100; %max counter value for each neuron
  cnt_L1 = [zeros(L1_Ne,1)*cnt_max];
  cnt_L2 = [zeros(L1_Ne,1)*cnt_max];
  cnt_L3 = [zeros(L1_Ne,1)*cnt_max];
  
  %Tau used for exponential increase for LTP, the lower 
  %this is the larger the rate of increase for LTP
  exp_tau_LTP = 150; 
  
  %Tau used for exponential decrease for LTD, the lower 
  %this is the larger the rate of increase for LTD
  exp_tau_LTD = 150;

  %Decay Rates
%   lambda_L1 = [ones(L1_Ne,1) * 0.02];
%   lambda_L2 = [ones(L2_Ne,1) * 0.02];
%   lambda_L3 = [ones(L3_Ne,1) * 0.02];
  lambda_L1 = [zeros(L1_Ne,1) * 0.02];
  lambda_L2 = [zeros(L2_Ne,1) * 0.02];
  lambda_L3 = [zeros(L3_Ne,1) * 0.02];  
  
  %Decremenet amount for inhib
  decAmount = 0.05;

  %Neuron Config Details
  Vth = 1;  %Spike Threshold

  %Refractory Periods for each neuron
  ref_period = 2;
  ref_cnt_L1 = [zeros(L1_Ne,1)];
  ref_cnt_L2 = [zeros(L2_Ne,1)];
  ref_cnt_L3 = [zeros(L3_Ne,1)];

  %Spikes
  spike_exists_L1 = [zeros(L1_Ne,1)];
  nxt_spike_exists_L1 = [zeros(L1_Ne,1)];
  nxt_spike_exists_L2 = [zeros(L2_Ne,1)];
  spike_exists_L2 = [zeros(L2_Ne,1)];
  spike_exists_L3 = [zeros(L3_Ne,1)];
  nxt_spike_exists_L3 = [zeros(L3_Ne,1)];
  
  %Populate Connections from inp to L1
  for j = 1:M_inp
    synapses_inptoL1(j) = 1;
  end

  %Populate Connections from L1 to L2
%   for i=1:L1_Ne
%     for j=1:M_L1
%       synapses_L1toL2(i,j) = 1;
%     end
%   end
    
  %Populate Connections from L2 to L3
  for i=1:L2_Ne
    for j=1:M_L2
      synapses_L2toL3(i,j) = 1;
    end
  end

  %If you are classifying
  if (mode == 0)
      weights_L1toL2 = wL1_L2;
      weights_L2toL3 = wL2_L3;
      iteration_count = 35;          %In ms
  else
      iteration_count = 350;        %In ms
  end

  %RasterPlot Data
  raster_L1 = [zeros(iteration_count,L1_Ne)];
  raster_L2 = [zeros(iteration_count,L2_Ne)];
  raster_L3 = [zeros(iteration_count,L3_Ne)];

  %VmemPlot Data
  Vmem_pl_L1 = [zeros(iteration_count,L1_Ne)];
  Vmem_pl_L2 = [zeros(iteration_count,L2_Ne)];
  Vmem_pl_L3 = [zeros(iteration_count,L3_Ne)];
   
  %Factor for counter decrementing in L1. The smaller this 
  %is the slower the decrement.
  c_factor_L1 = 0.8;
  
  %spikeMat = zeros(L1_Ne, iteration_count);
 
  spikeMat = poissonSpikeGen(100,iteration_count/1000,1);

  myPoissonSpikeTrain = rand(1, 784) < 100*0.001;
  
  totalSpikes_L2 = zeros(size(test_x,1),L2_Ne);
  totalSpikes_L3 = zeros(size(test_x,1),L3_Ne);
  
  %Total Iterations
  for j=1:iteration_count

      spike_snapshot = rand(size(test_x));
      inp_image = spike_snapshot <= test_x;
      
      %%%%%%%%%%%
      %Calculate Membrane Potentials for L2
      %%%%%%%%%%%
      %Iterate over all neurons in layer 2
      for ii=1:size(inp_image,1)
          for i=1:L2_Ne
            [Vmem_L2(i), ref_cnt_L2(i)] = neuron_new (inp_image(ii,:), synapses_L1toL2(:,i)', ... 
                                            weights_L1toL2(:,i)', Vmem_L2(i),  ...
                                            ref_cnt_L2(i), lambda_L2(i));
          end
      
          %Check for spikes
          for i=1:L2_Ne
              if (Vmem_L2(i) > Vth)
                  spike_exists_L2(i) = 1;
                  Vmem_L2(i) = 0;
                  ref_cnt_L2(i) = 5;
              else
                  spike_exists_L2(i) = 0;
              end
          end
          
          totalSpikes_L2(ii,:) = totalSpikes_L2(ii,:) + spike_exists_L2';

          %%%%%%%%%%%
          %Calculate Membrane Potentials for L3
          %%%%%%%%%%%

          %Iterate over all neurons in layer 3
          for i=1:L3_Ne
            [Vmem_L3(i), ref_cnt_L3(i)] = neuron_new (spike_exists_L2, synapses_L2toL3(:,i)', ... 
                                            weights_L2toL3(:,i)', Vmem_L3(i),  ...
                                            ref_cnt_L3(i), lambda_L3(i));              
          end

          %Check for spikes
          for i=1:L3_Ne
              if (Vmem_L3(i) > Vth)
                  spike_exists_L3(i) = 1;
                  Vmem_L3(i) = 0;
                  ref_cnt_L3(i) = 5;
              else
                  spike_exists_L3(i) = 0;
              end
          end          
          
          totalSpikes_L3(ii,:) = totalSpikes_L3(ii,:) + spike_exists_L3';
          
        %Set all the Vmems to 0
        Vmem_L2 = [zeroes(L2_Ne,1) * (-65)];
        Vmem_L3 = [zeroes(L3_Ne,1) * (-65)];
        
      end

  end

end

function [] = plotRaster(spikeMat, tVec)
    hold all;
    for trialCount = 1:size(spikeMat,1)
        spikePos = tVec(spikeMat(trialCount, :));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
                [trialCount-0.4 trialCount+0.4], 'k');
        end
    end
    ylim([0 size(spikeMat, 1)+1]);
end

function [weights] = normalize_weights(unnormalized_weights, old_weights)
    sum_oldweights = 0;
    sum_unnormalized_weights = 0;
    
    sum_oldweightsMat = [];
    sum_unnormalized_weightsMat = [];
    
    for j=1:size(old_weights,2)
        for k=1:size(old_weights,1)
            sum_oldweights = sum_oldweights + old_weights(k,j);
            sum_unnormalized_weights = sum_unnormalized_weights + unnormalized_weights(k,j);
        end
        sum_oldweightsMat = [sum_oldweightsMat sum_oldweights];
        sum_unnormalized_weightsMat = [sum_unnormalized_weightsMat sum_unnormalized_weights];
        sum_oldweights = 0;
        sum_unnormalized_weights = 0;
    end
    
    for i=1:size(old_weights,1)
       for j=1:size(old_weights,2)
%            for k=1:size(old_weights,1)
%                sum_oldweights = sum_oldweights + old_weights(k,j);
%            end
%            for k=1:size(old_weights,1)
%                sum_unnormalized_weights = sum_unnormalized_weights + unnormalized_weights(k,j);
%            end
           %weights(i,j) = unnormalized_weights(i,j) * (sum_unnormalized_weightsMat(j) / sum_oldweightsMat(j));
           weights(i,j) = unnormalized_weights(i,j) * (sum_oldweightsMat(j) / sum_unnormalized_weightsMat(j));
       end
    end
end

function [weights_post_out] = synapse_calc_post(counter_max, counter_post, weights_post, exp_tau)

  sc_fac = 0.6;

  weights_post_out = weights_post - weights_post * exp((counter_post+5)/exp_tau) * sc_fac;

end

function [weights_pre_out] = synapse_calc_pre(counter_max, counter_pre, weights_pre, exp_tau)

  sc_fac = 0.2;

  if (weights_pre >= 50)
      weights_pre_out = 50;
  else
      weights_pre_out = weights_pre + weights_pre * exp(-(counter_max + 5 - counter_pre)/exp_tau) * sc_fac;
  end

end

function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
    dt = 1/1000; % s
    nBins = floor(tSim/dt);
    spikeMat = rand(nTrials, nBins) < fr*dt;
    tVec = 0:dt:tSim-dt;
end

function [Vout] = decPot (Vin,decAmt)

  %decrement the potential
  Vout = Vin - decAmt;

end

function [Vout, ref] = neuron_new (spike_exists, synapses, syn_weights, Vin,  ...
                                      ref_in, lambda)

  %Scaling factor for increasing membrane potentials, the larger this the
  %larger the membrane increase
  sc_factor = 1.0;
                                        
  if ~ref_in
    %Vout = Vin - (Vin/(R*C)) + (I/C);
    temp = 0;
    for i = 1:size(synapses,2)
      %The spike_exists, synapses, and syn_weights are only for the layer that you are interested in
      %For example, if you are looking at a neuron in L2, the spike_exists contains if spike exists in neurons
      %in L1
      if (spike_exists(i) && synapses(i))
          temp = temp + syn_weights(i)*sc_factor;
      end 
    end
    Vout = Vin + temp - lambda;
    ref = ref_in;
  else
    ref = ref_in - 1;
    Vout = 0; % reset voltage
  end
    
end

function [spike, Vout, ref, counter_out] = neuron (spike_exists, synapses, syn_weights, V_th, Vin,  ...
                                      ref_in, lambda, spike_exists_ext, synapses_ext,  ...
                                      syn_weights_ext, counter, counter_max)

  %Scaling factor for increasing membrane potentials, the larger this the
  %larger the membrane increase
  sc_factor = 1.0;
  sc_factor_ext = 1.0;
  
  %Factor for reducing counter, the smaller this is the slower the
  %decrement
  c_factor = 1.0;
                                        
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
    temp2 = temp2 + spike_exists_ext*synapses_ext*syn_weights_ext*sc_factor_ext;

    %Basically, the Vmem doesn't go lower than 0
    if ((Vin + temp - lambda + temp2) < 0)
        Vout = 0;
    else
        Vout = Vin + temp - lambda + temp2;
    end
    
    ref = ref_in;
  else
    ref = ref_in - 1;
    Vout = 0; % reset voltage
  end
  
  spike = 0;
  counter_out = 0;
  
%   if (Vout > V_th && Vout ~= 0.2*V_th)
%     Vout = 50;  % emit spike
%     spike = 1;
%     counter_out = counter_max;
%     % ref = abs_ref; % set refractory counter
%   else
%     spike = 0;
%   end
% 
%   if (counter > 0)
%     counter_out = counter - counter*exp(-1/c_factor);
%   elseif (counter == 0 && spike ~= 1)
%     counter_out = counter;
%   end

end