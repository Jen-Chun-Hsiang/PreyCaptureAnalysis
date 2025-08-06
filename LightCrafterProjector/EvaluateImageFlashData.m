channel_id = 1;
Sampling_rate = 5.92;
Fz = size(ImageData, 1)*size(ImageData, 1)*Sampling_rate;
figure; 
a = [];
for i = 1:size(ImageData, 3)
    a = [a; reshape(squeeze(ImageData(:, :, i, channel_id))', [], 1)];
end
% a = smoothdata(a, "gaussian", round(Fz/1000));
t = (0:length(a)-1)/Fz;
plot(t, a);

