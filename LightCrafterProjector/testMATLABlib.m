addpath('./matlab-lcr');
lightCrafter = MockLightCrafter4500(60, 'uv');
lightCrafter.connect();
%%
stageClient = stage.core.network.StageClient();
stageClient.connect('localhost', '5678');
stageClient.setMonitorGamma(1);