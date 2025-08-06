num_dots =100;
square_length = 1000;
mass = 10;
[positions, dot_pos] = trajectoryGenerator(num_dots, square_length, mass);
figure;
Colors = parula(num_dots);
subplot(1, 2, 1); hold on
for i = 1:size(dot_pos, 1)-1
    plot(dot_pos([i i+1], 1), dot_pos([i i+1], 2), 'Color', Colors(i, :));
end

subplot(1, 2, 2); hold on;
for i = 1:size(positions, 1)-1
    plot(positions([i i+1], 1), positions([i i+1], 2), 'Color', Colors(i, :));
end

%%
timestamp = datestr(now, 'yyyymmddHHMM');
trajectory_name = sprintf('trajectory_%s.mat', timestamp);
moviefolder = './Simulation/Trajectories';
save(fullfile(moviefolder, trajectory_name), 'num_dots', 'square_length', 'positions', 'dot_pos');

%%
fv = linspace(0, 1, 20);                                % Normalised Frequencies
a = 1./(1 + fv*2);                                      % Amplitudes Of ‘1/f’
b = firls(42, fv, a);                                   % Filter Numerator Coefficients
figure(1)
freqz(b, 1, 2^17)                                       % Filter Bode Plot
N = 1E+5;
ns = rand(1, N);
invfn = filtfilt(b, 1, ns);                             % Create ‘1/f’ Noise
figure(2)
plot([0:N-1], invfn)                                    % Plot Noise In Time Domain
grid
FTn = fft(invfn-mean(invfn))/N;                         % Fourier Transform
figure(3)
plot([0:N/2], abs(FTn(1:N/2+1))*2)                      % Plot Fourier Transform Of Noise
grid
%%
fv = linspace(0, 1, 20);                                % Normalised Frequencies
a = zeros(1,20);
a(1:10) = 1./(1 + fv(1:10)*2);                          % Amplitudes Of  1/fv until 0.5
a(11:20) = a(10);                                       % after 0.5 it gets flat
b = firls(42, fv, a);                                   % Filter Numerator Coefficients
figure(1)
freqz(b, 1, 2^17)                                       % Filter Bode Plot
N = 1E+6;
ns = randn(1, N);
invfn = filtfilt(b, 1, ns);                             % Create ‘1/f’ Noise
figure(2)
plot([0:N-1], invfn)                                    % Plot Noise In Time Domain
grid
FTn = fft(invfn-mean(invfn))/N;                         % Fourier Transform
figure(3)
plot([0:N/2], abs(FTn(1:N/2+1))*2)                      % Plot Fourier Transform Of Noise
grid