clc;
close all;
clear all;
% Define the EIIP values for each amino acid
EIIP_VALUES = containers.Map( ...
    {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}, ...
    [0.0373, 0.0959, 0.0036, 0.1263, 0.0829, 0.0761, 0.0058, 0.0050, 0.0242, 0.0000, ...
     0.0000, 0.0823, 0.0829, 0.0946, 0.0198, 0.0829, 0.0941, 0.0548, 0.0516, 0.0057]);

% Input two protein sequences
sequence1 = 'MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVEGSCGF'; % human growth hemoglobin
sequence2 = 'MDLWQLLLTLALAGSSDAFSGSEATAAILSRAPWSLQSVNPGLKTNSSKEPKFTKCRSPERETFSCHWTDEVHHGTKNLGPIQLFYTRRNTQEWTQEWKECPDYVSAGENSCYFNSSFTSIWIPYCIKLTSNGGTVDEKCFSVDEIVQPDPPIALNWTLLNVSLTGIHADIQVRWEAPRNADIQKGWMVLEYELQYKEVNETKWKMMDPILTTSVPVYSLKVDKEYEVRVRSKQRNSGNYGEFSEVLYVTLPQMSQFTCEEDFYFPWLLIIIFGIFGLTVMLFVFLFSKQQRIKMLILPPVPVPKIKGIDPDLLKEGKLEEVNTILAIHDSYKPEFHSDDSWVEFIELDIDEPDEKTEESDTDRLLSSDHEKSHSNLGVKDGDSGRTSCCEPDILETDFNANDIHEGTSEVAQPQRLKGEADLLCLDQKNQNNSPYHDACPATQQPSVIQAEKNKPQPLPTEGAESTHQAAHIQLSNPSSLSNIDFYAQVSDITPAGSVVLSPGQKNKAGMSQCDMHPEMVSLCQENFLMDNAYFCEADAKKCIPVAPHIKVESHIQPSLNQEDIYITTESLTTAAGRPGTGEHVPGSEMPVPDYTSIHIVQSPQGLILNATALPLPDKEFLSSCGYVSTDQLNKIMP'; % human growth hemoglobin

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
max_period = 150;

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

rft_result1 = rft(eiip_sequence1,max_period);
rft_result2 = rft(eiip_sequence2,max_period);

Pq = abs(rft_result1 .* rft_result2);

%mean of Pq
mean = sum(Pq) / max_period ;

%SNR
snr = Pq / mean ;

% Step 3: Sort Pq values and find the top 20
sorted_snr = sort(snr, 'descend'); % Sort Pq in descending order
top_20_snr = sorted_snr(1:min(20, length(sorted_snr))); % Take the top 20 values
threshold = min(top_20_snr); % Threshold is the smallest value among the top 20

% Fixed ST-RFT Calculation with Proper Window Application
function ST_RFT = compute_ST_RFT(signal, max_period)
    N = length(signal);                % Full signal length
    window_size = 8;                   % Gaussian window length
    ST_RFT = zeros(1, max_period);     % Preallocate ST-RFT array

    alpha = 2; % Adjust alpha for Gaussian window
    w = exp(-0.5 * (alpha * ((1:window_size) / ((window_size-1)/2))).^2);

    for q = 1:max_period
        phi_q = numel(find(gcd(1:q, q) == 1));  % Compute φ(q)
        xq = 0;

        for n = 1:N
            win_idx = mod(n - 1, window_size) + 1; % Circular window logic
            xq = xq + signal(n) * conj(w(win_idx)) * ramanujan_sum(q, n);
        end

        ST_RFT(q) = (1 / phi_q) * (1 / window_size) * xq;
    end
end

ST_RFT = compute_ST_RFT(Pq, max_period);
ST_RFT_mag = abs(ST_RFT);



% Load the signal package required for findpeaks
pkg load signal

% Detect peaks in ST_RFT magnitude
[peak_vals, peak_locs] = findpeaks(ST_RFT_mag, 'MinPeakHeight', 0.1); % Adjust threshold if needed


% Compute Concentration Measure with Window Size of 8
function CM = compute_CM(ST_RFT, window_size)
    N = length(ST_RFT);                % Total sequence length
    CM = zeros(1, N);                  % Preallocate CM array

    for k = 1:N
        % Define window range ensuring no out-of-bound errors
        start_idx = max(1, k - floor(window_size / 2));
        end_idx = min(N, k + floor(window_size / 2));

        % Compute Concentration Measure within the window
        CM(k) = sum(abs(ST_RFT(start_idx:end_idx)).^1.1); % r = 1.1 as per your formula
    end

    % Normalize for better visualization
    CM = CM / max(CM);
end

% Compute ST-RFT with window size
ST_RFT = compute_ST_RFT(Pq, max_period);

% Compute and Plot CM
CM = compute_CM(ST_RFT, 8);


figure;
plot(1:max_period, rft_result1, 'k', 'LineWidth', 1.5);
title('rft result1');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Zoom in (adjust the y-axis limits)
axis([0 150 0 5]);  % X-axis from 0 to 10, Y-axis from -1 to 1

yticks(0:0.1:5);
xticks(0:1:150);

figure;
plot(1:max_period, rft_result2, 'k', 'LineWidth', 1.5);
title('rft result2');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Zoom in (adjust the y-axis limits)
axis([0 150 0 5]);  % X-axis from 0 to 10, Y-axis from -1 to 1

yticks(0:0.1:5);
xticks(0:1:150);

figure;
plot(1:max_period, Pq, 'k', 'LineWidth', 1.5);
hold on;
plot([0, max_period], [threshold, threshold], 'r'); % Plot horizontal line at y = threshold
title('cross correlation');
xlabel('Period');
ylabel('Magnitude');
grid on;

% Zoom in (adjust the y-axis limits)
axis([0 150 0 5]);  % X-axis from 0 to 10, Y-axis from -1 to 1

yticks(0:0.1:5);
xticks(0:1:150);

ST_RFT = compute_ST_RFT(Pq, max_period);

figure;
plot(1:max_period, abs(ST_RFT), 'k', 'LineWidth', 1.5);
title('ST-RFT (imaginary/magnitude Part)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
drawnow;

% Zoom in (adjust the y-axis limits)
axis([0 150 0 5]);  % X-axis from 0 to 10, Y-axis from -1 to 1
yticks(0:0.1:5);
xticks(0:1:150);

figure;
plot(1:max_period, CM, 'k', 'LineWidth', 1.5);
title('Concentration Measure with Window Size 8');
xlabel('Sample Index');
ylabel('Concentration Measure');
grid on;

% Adjust axis for better visualization
axis([0 150 0 5]);
yticks(0:0.1:5);
xticks(0:1:150);

% Plot ST-RFT magnitude
figure;
plot(1:max_period, ST_RFT_mag, 'k', 'LineWidth', 1.5);
hold on;

% Circle detected hotspots on ST-RFT plot
plot(peak_locs, peak_vals, '*', 'MarkerSize', 3, 'LineWidth', 2);

% Label the hotspot positions
for i = 1:length(peak_locs)
    text(peak_locs(i), peak_vals(i), sprintf('  %d', peak_locs(i)), ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'black');
end

title('ST-RFT Magnitude with Detected Hotspots');
xlabel('Period');
ylabel('Amplitude');
grid on;


axis([0 150 0 5]);
yticks(0:0.1:5);
xticks(0:1:150);
hold off;

