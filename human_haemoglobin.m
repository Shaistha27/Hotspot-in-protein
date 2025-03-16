clc;
close all;
clear all;

% Define the EIIP values for each amino acid
EIIP_VALUES = containers.Map( ...
    {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}, ...
    [0.0373, 0.0959, 0.0036, 0.1263, 0.0829, 0.0761, 0.0058, 0.0050, 0.0242, 0.0000, ...
     0.0000, 0.0823, 0.0829, 0.0946, 0.0198, 0.0829, 0.0941, 0.0548, 0.0516, 0.0057]);

% Input two protein sequences
sequence1 = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR';
sequence2 = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH';

% Convert sequences to EIIP values
eiip_sequence1 = zeros(1, length(sequence1));
eiip_sequence2 = zeros(1, length(sequence2));

for i = 1:length(sequence1)
    aa = sequence1(i);
    eiip_sequence1(i) = EIIP_VALUES(aa);
end

for i = 1:length(sequence2)
    aa = sequence2(i);
    eiip_sequence2(i) = EIIP_VALUES(aa);
end

% Pad sequences to the same length
max_len = max(length(eiip_sequence1), length(eiip_sequence2));
eiip_sequence1 = [eiip_sequence1, zeros(1, max_len - length(eiip_sequence1))];
eiip_sequence2 = [eiip_sequence2, zeros(1, max_len - length(eiip_sequence2))];

% Set maximum period for RFT computation
max_period = 50;

% Ramanujan sum function
function cq = ramanujan_sum(q, n)
    cq = 0;
    for k = 1:q
        if gcd(k, q) == 1
            cq = cq + exp(1i * 2 * pi * k * n / q);
        end
    end
end

% Ramanujan Fourier Transform function
function rft_result = rft(signal, max_period)
    rft_result = zeros(1, max_period);
    for q = 1:max_period
        phi_q = numel(find(gcd(1:q, q) == 1)); % Compute phi(q)
        xq = 0;
        for n = 1:length(signal)
            xq = xq + signal(n) * ramanujan_sum(q, n);
        end
        xq = (1 / phi_q) * xq;
        rft_result(q) = abs(xq);
    end
end

% Compute RFT for both sequences
rft_result1 = rft(eiip_sequence1, max_period);
rft_result2 = rft(eiip_sequence2, max_period);

% Cross-correlation
Pq = abs(rft_result1 .* rft_result2);

% Mean of Pq
mean_Pq = sum(Pq) / max_period;

% SNR calculation
snr = Pq / mean_Pq;

% Find top 20 SNR values
sorted_snr = sort(snr, 'descend');
top_20_snr = sorted_snr(1:min(20, length(sorted_snr)));
threshold = min(top_20_snr); % Threshold as the smallest of the top 20

% Find the period with the highest Pq value
[~, characteristic_period] = max(Pq);

% Check if period 1 is an outlier
if characteristic_period == 1
    [~, characteristic_period] = max(Pq(2:end));  % Find the next highest peak
    characteristic_period = characteristic_period + 1;  % Adjust index
end

% Display characteristic period
disp('Characteristic Period:');
disp(characteristic_period);


% Plot results
figure;
plot(1:max_period, rft_result1, 'k', 'LineWidth', 1.5);
title('RFT Result for Sequence 1');
xlabel('Period');
ylabel('Magnitude');
grid on;

figure;
plot(1:max_period, rft_result2, 'k', 'LineWidth', 1.5);
title('RFT Result for Sequence 2');
xlabel('Period');
ylabel('Magnitude');
grid on;

figure;
plot(1:max_period, Pq, 'k', 'LineWidth', 1.5);
hold on;
plot([0, max_period], [threshold, threshold], 'r'); % Threshold line
title('Cross Correlation');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Zoom in and adjust y-axis
axis([0 max_period 0 1]);
yticks(0:0.1:1);

