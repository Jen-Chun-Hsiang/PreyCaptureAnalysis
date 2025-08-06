function [ComputerIntensity, MinContrast] = ProjIntCalib(SigmParams, TargetContrast)

fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
flogit = @(param,xval) param(3)-log10(complex(((param(2)-param(1))./(xval-param(1))-1)))/param(4);

MinContrast = fsigm(SigmParams, 0);
TargetContrast(TargetContrast<MinContrast) = MinContrast;

try
    G = gpuArray(TargetContrast);
    ComputerIntensity = abs(gather(flogit(SigmParams, G)));
catch
    ComputerIntensity = abs(flogit(SigmParams, TargetContrast));
end
end