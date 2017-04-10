clc;
clear all;
close all;
fr = 100; % firing rate 100Hz
dt = 1/1000; % spikes in 1ms 
spike = zeros(1,10000); % max array size of spikes generated
tsc = zeros(1,10000); % max size of time stamp counter
i = 1; % integer defined
c = 0;
rng('shuffle');
while i<=10000 % current clock compared with previous value of clock (to generate spike for 10s)
    
   % pause(0.001); % 1ms time
    spike(1,i) = rand() < fr*dt; % checking if spiked or not (logical value)
    tsc(1,i) =(i-1)*10; % creating a timestamp for each spike value at interval of 10ms
    i = i+ 1;
end
array = [];
for a = 1:i-1 % run loop to get values
    if spike(1,a) == 1
        c = c+1;
        array(1,c) = tsc(1,a);% store the tsc of 1's in array
    end
end
isi = []; % inter spike intervals array

for i=2:c
    isi(1,i-1)=array(1,i)-array(1,i-1);
end

uniq=unique(isi);
[m n]=size(uniq);

h=histogram(isi,n,'Normalization','probability');

title('Interspike Interval Histogram')
xlabel('Interspike Interval (ms)')
ylabel('Probability of corresponding ISI bins')
% sd=std(isi)
% mn = mean(isi)
% cv=sd/mn
% spike_times=array; % array has time stamp of spikes over 10 seconds%
% fan_f=[]; % array for intervals from 1ms to 100ms
% for k=1:100 % for loop for 1ms,2ms,......,100ms
%     spike_count=[]; % spike count for each number of interval
%     for i=1:10000/k % number of intervals
%         s_count=0; 
%         for j=1:c % c is the size of spike_times ie number of spikes in 10 sec
%             if spike_times(1,j)<=i*k && spike_times(1,j)>(i-1)*k % condition for a specified interval
%                 s_count=s_count+1; %count the spikes in a given interval
%             end
%         end
%         spike_count(i)=s_count;  %store the spikes for each interval     
%     end
%     mean_f =mean(spike_count); % calc mean
%     var_f = var(spike_count); % calc variance 
%     fan_f(k)=var_f/mean_f; % calc fano factor
% end
% f_100=fan_f(100);
% f_50=fan_f(50);
% f_20=fan_f(20);
% f_10=fan_f(10);
% f_5=fan_f(5);
% f_2=fan_f(2);
% f_1=fan_f(1);
% T= table(f_1,f_2,f_5,f_10,f_20,f_50,f_100)
% 
% spike_raster = zeros(100,1000); % max array size of spikes generated
% tsc_spike = zeros(1,1000);
% for z = 1:100
%     i=1;
%     rng(z);
%     while i<=1000 % current clock compared with previous value of clock (to generate spike for 10s)
%         %pause(0.001); % constant rate is 100Hz (0.01s)
%         spike_raster(z,i) = rand() < fr*dt; % checking if spiked or not (logical value)
%         i = i+ 1;
%     end
% end
% i=1;
% while i<=1000
%     tsc_spike(1,i) =(i-1)/1000;
%     i=i+1;
% end
% 
% figure(2)
% 
% for trials=1:100
%     plot(tsc_spike,trials*spike_raster(trials,:),'+')
%     hold on
% end
% title('Spike Raster Plot')
% xlabel('Time (s)')
% ylabel('Poisson spike trains')