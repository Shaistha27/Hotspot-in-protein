clc;
close all;
clear all;
% Define the EIIP values for each amino acid
EIIP_VALUES = containers.Map( ...
    {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}, ...
    [0.0373, 0.0959, 0.0036, 0.1263, 0.0829, 0.0761, 0.0058, 0.0080, 0.0242, 0.0000, ...
     0.0000, 0.0823, 0.0829, 0.0946, 0.0198, 0.0829, 0.0941, 0.0548, 0.0516, 0.0057]);

% Input two protein sequences
sequence1 = 'MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR'; % mouse heamoglobin
sequence2 = 'MVLSADDKTNIKNCWGKIGGHGGEYGEEALQRMFAAFPTTKTYFSHIDVSPGSAQVKAHGKKVADALAKAADHVEDLPGALSTLSDLHAHKLRVDPVNFKFLSHCLLVTLACHHPGDFTPAMHASLDKFLASVSTVLTSKYR'; % rat heamoglobin

% Convert both sequences to EIIP values
eiip_sequence1 = zeros(1, length(sequence1));  % Preallocate array for sequence1
eiip_sequence2 = zeros(1, length(sequence2));  % Preallocate array for sequence2

% Loop over the first sequence
for i = 1:length(sequence1)
    aa = sequence1(i);  % Get the amino acid at position i in sequence1
    eiip_sequence1(i) = EIIP_VALUES(aa);  % Get the corresponding EIIP value
end

% Loop over the second sequence
for i = 1:length(sequence2)
    aa = sequence2(i);  % Get the amino acid at position i in sequence2
    eiip_sequence2(i) = EIIP_VALUES(aa);  % Get the corresponding EIIP value
end

% Pad the sequences to the same length
max_len = max(length(eiip_sequence1), length(eiip_sequence2));  % Get the maximum length
eiip_sequence1 = [eiip_sequence1, zeros(1, max_len - length(eiip_sequence1))];  % Pad sequence1
eiip_sequence2 = [eiip_sequence2, zeros(1, max_len - length(eiip_sequence2))];  % Pad sequence2

% Set maximum period for RFT computation
max_period = 80;

%ramanujan sum
function cq =  ramanujan_sum(q,n);
cq=0;
for k=1:q
   if gcd(k,q)==1;
   cq = cq + exp ( i * 2 * pi * k * n / q );
      end
    end
  end

% ramanujan fourier transform coefficient
function rft_result = rft(signal, max_period)
    rft_result = zeros(1, max_period);
for q = 1:max_period
  phi_q = numel(find(gcd(1:q, q) == 1)); % Compute phi(q)
   xq = 0;
for n = 1:length(signal)
   xq = xq + signal(n) * ramanujan_sum(q, n);  % Accumulate sum
    end
    xq = (1 / phi_q) * xq;
    rft_result(q) = abs(xq);
 end
end

rft_result1 = rft(eiip_sequence1,80);
rft_result2 = rft(eiip_sequence2,80);

Pq = abs(rft_result1 .* rft_result2);

%mean of Pq
mean = sum(Pq) / 80 ;

%SNR
snr = Pq / mean ;

% Step 3: Sort Pq values and find the top 20
sorted_snr = sort(snr, 'descend'); % Sort Pq in descending order
top_20_snr = sorted_snr(1:min(20, length(sorted_snr))); % Take the top 20 values
threshold = min(top_20_snr)+ 0.9; % Threshold is the smallest value among the top 20


figure;
plot(1:80, rft_result1, 'k', 'LineWidth', 1.5);
title('rft result1');
xlabel('Period');
ylabel('Magnitude');
grid on;

figure;
plot(1:80, rft_result2, 'k', 'LineWidth', 1.5);
title('rft result1');
xlabel('Period');
ylabel('Magnitude');
grid on;

figure;
plot(1:80, Pq, 'k', 'LineWidth', 1.5);
hold on;
plot([0, 80], [threshold, threshold], 'r'); % Plot horizontal line at y = threshold
title('cross correlation');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Zoom in (adjust the y-axis limits)
axis([0 80 0 1]);  % X-axis from 0 to 10, Y-axis from -1 to 1

yticks(0:0.1:1);
