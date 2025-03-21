clc;
close all;
clear all;

% Define the EIIP values for each amino acid
EIIP_VALUES = containers.Map( ...
    {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}, ...
    [0.0373, 0.0959, 0.0036, 0.1263, 0.0829, 0.0761, 0.0058, 0.0050, 0.0242, 0.0000, ...
     0.0000, 0.0823, 0.0829, 0.0946, 0.0198, 0.0829, 0.0941, 0.0548, 0.0516, 0.0057]);

% Input two protein sequences
sequence1 = 'MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR';
sequence2 = 'MVLSADDKTNIKNCWGKIGGHGGEYGEEALQRMFAAFPTTKTYFSHIDVSPGSAQVKAHGKKVADALAKAADHVEDLPGALSTLSDLHAHKLRVDPVNFKFLSHCLLVTLACHHPGDFTPAMHASLDKFLASVSTVLTSKYR';

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

% Set maximum period for RFT computation
max_period = 50;

% Ramanujan sum
function cq = ramanujan_sum(q, n)
    cq = 0;
    for k = 1:q
        if gcd(k, q) == 1
            cq = cq + exp(i * 2 * pi * k * n / q);
        end
    end
end

% Ramanujan Fourier Transform coefficient
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

% Compute RFT results for both sequences
rft_result1 = rft(eiip_sequence1, 50);
rft_result2 = rft(eiip_sequence2, 50);

% Plot the RFT result for the first sequence
figure;
plot(1:50, rft_result1, 'r', 'LineWidth', 1.5);
hold on;

% Find the characteristic periods by detecting peaks
[~, locs1] = findpeaks(rft_result1); % Get the indices of the peaks

% Highlight the peaks on the plot
plot(locs1, rft_result1(locs1), 'bo', 'MarkerFaceColor', 'b'); % Plot peaks as blue dots

% Add titles, labels, and grid
title('RFT Result 1 with Characteristic Periods');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Optional: Add vertical lines at the characteristic periods
for i = 1:length(locs1)
    xline(locs1(i), '--k'); % Dashed black lines at peak locations
end

% Set yticks
yticks(1:3);

hold off;

% Plot the RFT result for the second sequence
figure;
plot(1:50, rft_result2, 'g', 'LineWidth', 1.5);
hold on;

% Find the characteristic periods by detecting peaks
[~, locs2] = findpeaks(rft_result2); % Get the indices of the peaks

% Highlight the peaks on the plot
plot(locs2, rft_result2(locs2), 'bo', 'MarkerFaceColor', 'b'); % Plot peaks as blue dots

% Add titles, labels, and grid
title('RFT Result 2 with Characteristic Periods');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Optional: Add vertical lines at the characteristic periods
for i = 1:length(locs2)
    xline(locs2(i), '--k'); % Dashed black lines at peak locations
end

% Set yticks
yticks(1:3);

hold off;

