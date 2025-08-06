
deviations = STAmat - mean_spike_triggered_stimulus;
stc = (deviations' * deviations) / (total_spikes - 1);