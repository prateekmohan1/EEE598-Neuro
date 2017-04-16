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