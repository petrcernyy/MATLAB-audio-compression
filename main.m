clc, clear, close all

% Load a table from ISO 11172-3
load('Table_D1c.mat');
load('Map.mat');            % This vector remaps the length of table to 256 values that 
                            % are obtained from FFT

[x, Fs] = audioread('rickroll_smaller.wav');

if (length(x(:,2)) > 1)
    x = x(:,1)';
end

xr = zeros(1,length(x));

FFTLength = 512;
start = 1;

% Creates levels for quantization 4, 6, 8, 12 bits
bins_4 = linspace(min(x), max(x), 2^4);
bins_6 = linspace(min(x), max(x), 2^6);
bins_8 = linspace(min(x), max(x), 2^8);
bins_12 = linspace(min(x), max(x), 2^12);
bins_16 = linspace(min(x), max(x), 2^16);
bins = {bins_4; bins_6; bins_8; bins_12};

c_4bit = 0;         % Counter, how many quantization for 4 bits was kept
c_6bit = 0;         % Counter, how many quantization for 6 bits was kept
c_8bit = 0;         % Counter, how many quantization for 8 bits was kept
c_12bit = 0;        % Counter, how many quantization for 12 bits was kept
c_all = 0;          % How many times was the psychoacoustic model performed

% Segment by 512 samples
while (start <= (length(x) - FFTLength + 1))
   
    s = x(start:start + FFTLength - 1);
    
    % Use psychoacoustic only if the input signal is nonzero
    if nnz(s) > 0

        h = (sqrt(8/3) * hanning(512, 'periodic'))';    % Hanning window
        

        %% FFT a Power density
        x_F = abs(fft(s.*h, FFTLength));
        x_F = x_F(1:(FFTLength/2));
        PSD_dB = 20 * log10(x_F / FFTLength);
    
        PSD_Normed = PSD_dB - max(PSD_dB) + 96;         % Norm the max value to 96 dB
        
        %% Tones searching
        tones = [];
        counter = 1;
        
        % Firt the local maxima are found. Then the way they affect the neighbouring
        % frequencies bands. The bands are chosen based on the tone
        % frequency
        for k = 1:FFTLength/2
            if (k > 2 && k < 250)
                if ((PSD_Normed(k) > PSD_Normed(k-1)) && (PSD_Normed(k) >= PSD_Normed(k+1)))
                    if (k > 2 && k < 63)
                        j = [-2 2];
                    elseif (k >= 63 && k < 127)
                        j = [-3 -2 2 3];
                    elseif (k >= 127 && k <= 250)
                        j = [-6 -5 -4 -3 -2 2 3 4 5 6];
                    end
    
                    is_tone = 1;
                    for i = j
                        if ~(PSD_Normed(k) - PSD_Normed(k+i) >= 7)
                            is_tone = 0;
                        end
                    end
    
                    if (is_tone)
                        tones(counter,1) = k;
                        tones(counter,2) = 10*log10(10^(PSD_Normed(k-1)/10) + 10^(PSD_Normed(k)/10) + 10^(PSD_Normed(k+1)/10));
                        counter = counter + 1;
                    end
                end
            end
        end

        %% Finding the masking theshold of each tone

        % For every tone a masking threshold is found. Go through the table
        % from ISO 11172-3, the distance from tone is calculated and based
        % on this a the masking threshold is calculated
        tones_func = zeros(length(tones), length(Table_D1c)) -200;

        for i = 1:length(Table_D1c)

            zi = Table_D1c(i, 2);  
           
            if not(isempty(tones))
                for k = 1:length(tones(:,1))
                    j  = tones(k, 1);
	                zj = Table_D1c(Map(j), 2);
	                dz = zi - zj;
	              
                    if (dz >= -3 && dz < 8)

                        avtm = -1.525 - 0.275 * zj - 4.5;
	                 
	                    if (-3 <= dz && dz < -1)
	                        vf = 17 * (dz + 1) - (0.4 * PSD_Normed(j) + 6);
	                    elseif (-1 <= dz && dz < 0)
	                        vf = (0.4 * PSD_Normed(j) + 6) * dz;
	                    elseif (0 <= dz && dz < 1)
	                        vf = -17 * dz;
	                    elseif (1 <= dz && dz < 8)
	                        vf = - (dz - 1) * (17 - 0.15 * PSD_Normed(j)) - 17;
	                    end
                 
	                    tones_func(k, i) = tones(k, 2) + avtm + vf;

                    end
                end
            end
        end


        %% Global threshold
        
        % Just takes the values of absolute threshold and sum the masking
        % thresholds of tones
        N = length(Table_D1c);
        m = length(tones_func(:, 1));
        global_thresh = zeros(N,1);

        for i = 1:N
           
            temp = 10^(Table_D1c(i,3) / 10);        % Transfer to magnitude to be able to sum

            if not(isempty(tones_func))
                for j = 1:m
                    temp = temp + 10^(tones_func(j, i) / 10);
                end
            end
           
            global_thresh(i) = 10 * log10(temp);    % Transfer back to dB
        end
    
        %% Quantization Error DFT
        
        % Quantize the signal to 4, 6, 8 and 12 bits. Calculate the error
        % and transfer to PSD.
        PSD_error_CB_all = zeros(length(bins), length(Table_D1c));
        kvanted = zeros(length(bins), length(s));

        for i = 1:length(bins)
      
            s_kvantovane = interp1(bins{i}, bins{i}, s, 'nearest');
            kvanted(i,:) = s_kvantovane;

            kvantovani_error = s - s_kvantovane;
            error_F = abs(fft(kvantovani_error .* h, FFTLength));
            error_F = error_F(1:(FFTLength/2));
            error_PSD_dB = 20 * log10(error_F / FFTLength);
    
            error_PSD_Normed = error_PSD_dB + 96;

            % Map the error to length of table from ISO 11172-3
            for j = 1:length(Table_D1c)
                bin_indices = find(Map == j);
                if ~isempty(bin_indices)
                    PSD_error_CB_all(i, j) = max(error_PSD_Normed(bin_indices));
                else
                    PSD_error_CB_all(i, j) = -200;
                end
            end

        end
        
        current_bit = 0;
        for i = 1:length(bins)
            flag = 1;
            for j = 1:length(PSD_error_CB_all(1,:))
                if PSD_error_CB_all(i,j) > global_thresh(j)
                    flag = 0;
                end
            end
            if (flag)
                current_bit = i;
                break;
            end
        end
    
        switch current_bit
            case 1
                filtered = kvanted(1,:);
                c_4bit = c_4bit + 1;
            case 2
                filtered = kvanted(2,:);
                c_6bit = c_6bit + 1;
            case 3
                filtered = kvanted(3,:);
                c_8bit = c_8bit + 1;
            case 4
                filtered = kvanted(4,:);
                c_12bit = c_12bit + 1;
            otherwise
                filtered = kvanted(4,:);
        end
        c_all = c_all + 1;
        
    else
        filtered = s;
    end
    
    xr(start:start + FFTLength - 1) = xr(start:start + FFTLength - 1) + filtered;
    start = start + FFTLength - FFTLength/16;               % Overllap 

end


uni_puvod = length(unique(x));
uni_new = length(unique(xr));

fprintf('Counter: %d\n', c_all);
fprintf('Counter 4 bits: %d\n', c_4bit);
fprintf('Counter 6 bits: %d\n', c_6bit);
fprintf('Counter 8 bits: %d\n', c_8bit);
fprintf('Counter 12 bits: %d\n', c_12bit);
fprintf('Number of unique characters in original data: %d\n', uni_puvod);
fprintf('Number of unique characters in filtered: %d\n', uni_new);



%% Huffmann coding

% This takes a while

[root, encoded_data] = Huffman_encode(xr);
save('encoded_data.mat', 'encoded_data');
decoded_data = Huffman_decode(encoded_data, root);
error = max(abs(encoded_data - decoded_data));




%% Plotting

s = x(53761:53761 + FFTLength - 1);

h = (sqrt(8/3) * hanning(512, 'periodic'))';

x_F = abs(fft(s.*h, FFTLength));
x_F = x_F(1:(FFTLength/2));
PSD_dB = 20 * log10(x_F / FFTLength);

PSD_Normed = PSD_dB - max(PSD_dB) + 96;

tones = [];
counter = 1;

for k = 1:FFTLength/2
    if (k > 2 && k < 250)
        if ((PSD_Normed(k) > PSD_Normed(k-1)) && (PSD_Normed(k) >= PSD_Normed(k+1)))
            if (k > 2 && k < 63)
                j = [-2 2];
            elseif (k >= 63 && k < 127)
                j = [-3 -2 2 3];
            elseif (k >= 127 && k <= 250)
                j = [-6 -5 -4 -3 -2 2 3 4 5 6];
            end

            is_tone = 1;
            for i = j
                if ~(PSD_Normed(k) - PSD_Normed(k+i) >= 7)
                    is_tone = 0;
                end
            end

            if (is_tone)
                tones(counter,1) = k;
                tones(counter,2) = 10*log10(10^(PSD_Normed(k-1)/10) + 10^(PSD_Normed(k)/10) + 10^(PSD_Normed(k+1)/10));
                counter = counter + 1;
            end
        end
    end
end

tones_func = zeros(length(tones), length(Table_D1c)) -200;

for i = 1:length(Table_D1c)

    zi = Table_D1c(i, 2);  
   
    if not(isempty(tones))
        for k = 1:length(tones(:,1))
            j  = tones(k, 1);
            zj = Table_D1c(Map(j), 2);
            dz = zi - zj;
          
            if (dz >= -3 && dz < 8)

                avtm = -1.525 - 0.275 * zj - 4.5;
             
                if (-3 <= dz && dz < -1)
                    vf = 17 * (dz + 1) - (0.4 * PSD_Normed(j) + 6);
                elseif (-1 <= dz && dz < 0)
                    vf = (0.4 * PSD_Normed(j) + 6) * dz;
                elseif (0 <= dz && dz < 1)
                    vf = -17 * dz;
                elseif (1 <= dz && dz < 8)
                    vf = - (dz - 1) * (17 - 0.15 * PSD_Normed(j)) - 17;
                end
         
                tones_func(k, i) = tones(k, 2) + avtm + vf;

            end
        end
    end
end

N = length(Table_D1c);
m = length(tones_func(:, 1));
global_thresh = zeros(N,1);

for i = 1:N
   
    temp = 10^(Table_D1c(i,3) / 10);

    if not(isempty(tones_func))
        for j = 1:m
            temp = temp + 10^(tones_func(j, i) / 10);
        end
    end
   
    global_thresh(i) = 10 * log10(temp);
end

PSD_error_CB_all = zeros(length(bins), length(Table_D1c));
kvanted = zeros(length(bins), length(s));

for i = 1:length(bins)

    s_kvantovane = interp1(bins{i}, bins{i}, s, 'nearest');
    kvanted(i,:) = s_kvantovane;

    kvantovani_error = s - s_kvantovane;
    error_F = abs(fft(kvantovani_error .* h, FFTLength));
    error_F = error_F(1:(FFTLength/2));
    error_PSD_dB = 20 * log10(error_F / FFTLength);

    error_PSD_Normed = error_PSD_dB + 96;

    for j = 1:length(Table_D1c)
        bin_indices = find(Map == j);
        if ~isempty(bin_indices)
            PSD_error_CB_all(i, j) = max(error_PSD_Normed(bin_indices));
        else
            PSD_error_CB_all(i, j) = -200;
        end
    end

end


%% Plot
figure;
hold on;
plot(PSD_Normed, 'k');
plot(tones(:,1), tones(:,2), 'ro');
grid minor
xlabel('Index [-]');
ylabel('Power spectrum density [dB]');
xlim([0 256]);
ylim([-40 120]);
pause;

figure;
hold on
grid minor
ylim([-40 120]);
xlim([Table_D1c(1,2) Table_D1c(end,2)]);
xlabel('Frequency [bark]');
ylabel('Power [dB]');
plot(Table_D1c(:,2), Table_D1c(:,3), 'k');
pause;
plot(Table_D1c(:,2), global_thresh, 'k--', 'LineWidth', 1.5);
pause;
for j = 1:length(tones)
    plot(Table_D1c(:,2), tones_func(j, :), 'r:');
end
pause;

plot(Table_D1c(:,2), PSD_error_CB_all(1,:), 'LineWidth', 2, 'Color', "#003f5c");
pause;
plot(Table_D1c(:,2), PSD_error_CB_all(2,:), 'LineWidth', 2, 'Color', "#58508d");
pause;
plot(Table_D1c(:,2), PSD_error_CB_all(3,:), 'LineWidth', 2, 'Color', "#bc5090");
pause;
plot(Table_D1c(:,2), PSD_error_CB_all(4,:), 'LineWidth', 2, 'Color', "#ff6361");
pause;

close all;






