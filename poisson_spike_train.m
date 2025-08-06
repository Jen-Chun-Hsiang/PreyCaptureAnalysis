% ---------------------- Poisson Spike Generation ----------------------
% Generate spike times from continuous rate trace using Poisson process
% rate: 1×T vector of firing rates (Hz), dt: bin duration in seconds
function spike_times = poisson_spike_train(rate, dt)
    % Ensure row vector
    rate = rate(:)';
    T = numel(rate);
    spikes = [];
    for t = 1:T
        lambda = rate(t) * dt;           % expected count in bin
        n_spikes = poissrnd(lambda);     % draw from Poisson
        bin_start = (t-1) * dt;
        % uniform spike times within bin
        times = bin_start + rand(1, n_spikes) * dt;
        spikes = [spikes, times];
    end
    spike_times = sort(spikes);
end