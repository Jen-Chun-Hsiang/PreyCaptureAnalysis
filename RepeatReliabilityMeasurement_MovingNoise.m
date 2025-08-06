clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
recording_id = 115;
Is_Calculate_STC = 0;
switch recording_id
    case 1
        recordingname = 'a120723';
        load([ephys_data_folder recordingname '_0010_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 2
        recordingname = 'b120723';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 3
        recordingname = 'c120723';
        load([ephys_data_folder recordingname '_0005_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 4
        recordingname = 'd010224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_2_MovingNoise_1.mat']);
    case 5
        recordingname = 'a030124';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 100;
    case 6
        recordingname = 'c030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 7
        recordingname = 'd030124';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 8
        recordingname = 'e030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 9
        recordingname = 'f030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 150;
    case 10
        recordingname = 'a030624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 11
        recordingname = 'a030924';
        load([ephys_data_folder recordingname '_0000_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 120;
    case 12
        recordingname = 'b030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 13
        recordingname = 'c030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 14
        recordingname = 'd030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 15
        recordingname = 'a032924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 16
        recordingname = 'b032924';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 2000;
    case 17
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 18
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 19
        recordingname = 'a040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 25
        recordingname = 'b040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 26
        recordingname = 'c040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 40
        recordingname = 'd040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 41
        recordingname = 'e040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1500;
    case 42
        recordingname = 'f040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1100;
    case 43
        recordingname = 'g040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;

    case 20
        recordingname = 'a040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 21
        recordingname = 'b040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 22
        recordingname = 'c040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 350;
    case 23
        recordingname = 'd040224';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 24
        recordingname = 'e040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 27
        recordingname = 'a040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 150;
    case 28
        recordingname = 'b040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 29
        recordingname = 'c040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 30
        recordingname = 'd040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 31
        recordingname = 'e040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 32
        recordingname = 'a040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 33
        recordingname = 'b040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 34
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 35
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;
    case 36
        recordingname = 'd040624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
    case 37
        recordingname = 'e040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
    case 38
        recordingname = 'f040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 39
        recordingname = 'g040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 44
        recordingname = 'a042424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 45
        recordingname = 'b042424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 46
        recordingname = 'd042424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 47
        recordingname = 'e042424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 48
        recordingname = 'e042424';
        save_recording_name = 'e042424B';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
    case 49
        recordingname = 'd042424';
        save_recording_name = 'd042424B';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 50
        recordingname = 'b042424';
        save_recording_name = 'b042424B';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 51
        recordingname = 'f042424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 52
        recordingname = 'g042424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 53
        recordingname = 'h042424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 54
        recordingname = 'ih042424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 55
        recordingname = 'h042424';
        save_recording_name = 'h042424B';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 56
        recordingname = 'g042424';
        save_recording_name = 'g042424B';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 57
        recordingname = 'f042424';
        save_recording_name = 'f042424B';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1200;
    case 58
        recordingname = 'a081024';
        save_recording_name = 'a04242401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 59
        recordingname = 'b081024';
        save_recording_name = 'b08102401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 270;
    case 60
        recordingname = 'c081024';
        save_recording_name = 'c08102401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 61
        recordingname = 'd081024';
        save_recording_name = 'd08102401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 62
        recordingname = 'a081224';
        save_recording_name = 'a08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 63
        recordingname = 'b081224';
        save_recording_name = 'b08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 64
        recordingname = 'c081224';
        save_recording_name = 'c08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 65
        recordingname = 'd081224';
        save_recording_name = 'd08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 66
        recordingname = 'e081224';
        save_recording_name = 'e08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 67
        recordingname = 'f081224';
        save_recording_name = 'f08122401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 68
        recordingname = 'a090724';
        save_recording_name = 'a09072401';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 69
        recordingname = 'b090724';
        save_recording_name = 'b09072401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 70
        recordingname = 'c090724';
        save_recording_name = 'c09072401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 71
        recordingname = 'd090724';
        save_recording_name = 'd09072401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 72
        recordingname = 'a091024';
        save_recording_name = 'a09102401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 73
        recordingname = 'b091024';
        save_recording_name = 'b09102401';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
     case 74
        recordingname = 'c091024';
        save_recording_name = 'c09102401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;
     case 75
        recordingname = 'd091024';
        save_recording_name = 'd09102401';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
     case 76
        recordingname = 'a092324';
        save_recording_name = 'a09232401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
     case 77
        recordingname = 'a092324';
        save_recording_name = 'a09232402';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 250;
     case 78
        recordingname = 'b092324';
        save_recording_name = 'b09232401';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
     case 79
        recordingname = 'b092324';
        save_recording_name = 'b09232402';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 80
        recordingname = 'e100724';
        save_recording_name = 'e100724';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
     case 81
        recordingname = 'f100724';
        save_recording_name = 'f100724';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
     case 82
        recordingname = 'a101224';
        save_recording_name = 'a101224';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 350;
      case 83
        recordingname = 'b101224';
        save_recording_name = 'b101224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 350;
      case 84
        recordingname = 'c101224';
        save_recording_name = 'c101224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 85
        recordingname = 'd101224';
        save_recording_name = 'd101224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;
      case 86
        recordingname = 'e101224';
        save_recording_name = 'e101224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 87
        recordingname = 'a101424';
        save_recording_name = 'a101424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 2000;
     case 88
        recordingname = 'b101424';
        save_recording_name = 'b101424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight =500;
     case 89
        recordingname = 'c101424';
        save_recording_name = 'c101424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight =1000;
      case 90
        recordingname = 'd101424';
        save_recording_name = 'd101424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 91
        recordingname = 'e101424';
        save_recording_name = 'e101424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 92
        recordingname = 'f101424';
        save_recording_name = 'f101424';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
      case 93
        recordingname = 'a101624';
        save_recording_name = 'a101624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 250;
      case 94
        recordingname = 'b101624';
        save_recording_name = 'b101624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 95
        recordingname = 'c101624';
        save_recording_name = 'c101624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 96
        recordingname = 'd101624';
        save_recording_name = 'd101624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
      case 97
        recordingname = 'e101624';
        save_recording_name = 'e101624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 98
        recordingname = 'b101924';
        save_recording_name = 'b101924';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
      case 99
        recordingname = 'c101924';
        save_recording_name = 'c101924';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1200;
      case 100
        recordingname = 'd101924';
        save_recording_name = 'd101924';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
      case 101
        recordingname = 'e101924';
        save_recording_name = 'e101924';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
      case 102
        recordingname = 'b103124';
        save_recording_name = 'b103124';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
      case 103
        recordingname = 'c103124';
        save_recording_name = 'c103124';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
      case 104
        recordingname = 'e103124';
        save_recording_name = 'e103124';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 150;
     case 105
        recordingname = 'f103124';
        save_recording_name = 'f103124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 120;
     case 106
        recordingname = 'a110424';
        save_recording_name = 'a110424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 107
        recordingname = 'b110424';
        save_recording_name = 'b110424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
     case 108
        recordingname = 'c110424';
        save_recording_name = 'c110424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
     case 109
        recordingname = 'd110424';
        save_recording_name = 'd110424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 100;
     case 110
        recordingname = 'e110424';
        save_recording_name = 'e110424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 120;
     case 111
        recordingname = 'f110424';
        save_recording_name = 'f110424';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
     case 112
        recordingname = 'g110424';
        save_recording_name = 'g110424';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 113
        recordingname = 'a110924';
        save_recording_name = 'a110924';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 114
        recordingname = 'b110924';
        save_recording_name = 'b110924';
        load([ephys_data_folder recordingname '_0008_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_008_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
     case 115
        recordingname = 'c110924';
        save_recording_name = 'c110924';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
     case 116
        recordingname = 'a111224';
        save_recording_name = 'a111224';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
     case 117
        recordingname = 'b111224';
        save_recording_name = 'b111224';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
end
%%
OUT = MNS1_OUT;
IN = MNS1_IN;
% OLED.height = 600;
% OLED.width = 800;
% OLED.pixelSize = 2.5; %size of each pixel in micron on retina
clear MNS1_OUT MNS1_IN

assert(abs(length(sectionData)/(OUT.endT-OUT.startT) -10000)<1);
%%
if ~exist('save_recording_name', 'var')
    save_recording_name = recordingname;
end
%% Detrend
close all
figure;
subplot(2, 1, 1); hold on
plot(sectionData, 'k');
% plot(smoothdata(sectionData, 'movmedian', 1*10000), 'r');
lpass = lowpass(sectionData, 0.05, 10000, 'Steepness', 0.99);
plot(lpass, 'r');

subplot(2, 1, 2); hold on
plot(sectionData, 'k');
plot(sectionData-lpass, 'r');
%%
% keyboard;
%% Use the figure plotted above to determine MinPeakHeight
% perform manual checking and input
%%
sectionData = sectionData-lpass;
%%
Fz = 100;
WinT = [-0.5 0];
% convert trace to spike or firing rate
[pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 30);
figure; hold on
plot(sectionData(1:locs(end)), 'k');
plot(locs(1:end), -pks(1:end), 'rx');

% keyboard;
% downsample to Fz
%%
fs = 10000;
signals = zeros(size(sectionData));
signals(locs) = 1;
% t = (0:length(sectionData)-1)/fs;
downsample_size = round(fs/Fz);
add_more = mod(length(signals), downsample_size);
if add_more ~= 0
    add_more = downsample_size - add_more;
    sig = [signals(:)' zeros(1, add_more)];
end
assert(mod(length(sig), downsample_size) == 0);
sig = reshape(sig(:), downsample_size, []);
sig = sum(sig, 1);

%%
UniBlc = unique(OUT.FrmTable(:, 2));
nBlock = length(UniBlc);
bTab = nan(nBlock, 8);
bcount = 1;
for b = 1:length(UniBlc)
    cTab = OUT.FrmTable(OUT.FrmTable(:, 2) == UniBlc(b) , :);
    dT = diff(cTab(:, 1));
    bTab(bcount, :) = [cTab(1, 1), cTab(end, 1), cTab(1, 3),...
        cTab(1, 2), 1/mean(dT),1/min(dT), 1/max(dT),...
        size(cTab, 1)];
    bcount = bcount + 1;
end
bTab(:, 9) = bTab(:, 2)-bTab(:, 1);

%-------------------------------------------------------------------
% Compute binned spike‐count traces for each repeated stimulus trial
%-------------------------------------------------------------------
% Assumptions:
%   • OUT            – struct containing FrmTable (NxM) and startT
%   • bTab           – [nBlocks×8] matrix, with:
%                       • bTab(:,3): block type
%                       • bTab(:,4): block ID (matches FrmTable(:,2))
%   • sig            – binary spike vector (1 at spike times, 0 elsewhere)
%   • Fz             – sampling frequency (Hz), so bin width = 1/Fz
%
% Output:
%   repeatTraces    – [nRepeats × nBinsMax] matrix; each row is one repeat,
%                     columns are spike counts in successive 1/Fz‐s bins.
%-------------------------------------------------------------------

%% 1) Identify the repeats
repeatIdx   = find(bTab(:,3) == 11);    % indices of repeating‐stimulus blocks
nRepeats    = numel(repeatIdx);

%% 2) Compute number of bins per repeat, find the maximum
binsPerRep = zeros(nRepeats,1);
for i = 1:nRepeats
    blockID   = bTab(repeatIdx(i),4);
    frmRows   = OUT.FrmTable(:,2) == blockID;
    tAll      = OUT.FrmTable(frmRows,1) - OUT.startT;       % times (s) relative to start
    duration  = tAll(end) - tAll(1);                        % stimulus duration
    binsPerRep(i) = floor(duration * Fz);                   % full bins per trial
end
nBinsMax = max(binsPerRep);

%% 3) Pre‐allocate output matrix (pad shorter repeats with NaN)
repeatTraces = nan(nRepeats, nBinsMax);

%% 4) Precompute spike times from 'sig'
cData = (find(sig > 0) - 1) / Fz;   % vector of spike times in seconds

%% 5) Fill in each repeat’s binned counts
kernel_size = 100;
sigma = 0.08;
is_smoothed = 1;

for i = 1:nRepeats
    blockID   = bTab(repeatIdx(i),4);
    frmRows   = OUT.FrmTable(:,2) == blockID;
    tAll      = OUT.FrmTable(frmRows,1) - OUT.startT;
    t0        = tAll(1);             % start time of this repeat
    nBins     = binsPerRep(i);
    
    % count spikes in each 1/Fz bin
    counts    = zeros(1, nBins);
    for b = 1:nBins
        sT = t0 + (b-1)/Fz;          % bin start
        eT = sT + 1/Fz;              % bin end
        counts(b) = sum(cData >= sT & cData < eT);
    end
    
    % store, leaving trailing columns as NaN if this repeat is shorter
    if is_smoothed
        counts = gaussian_smooth_1d(counts, kernel_size, Fz, sigma);
    end

    repeatTraces(i,1:nBins) = counts;
end
%
% 'repeatTraces' is now your [nRepeats × nBinsMax] output matrix.
figure; plot(repeatTraces')
trace = repeatTraces;
trace(:, any(isnan(repeatTraces), 1)) = [];
R = corr(trace');

mask = triu(true(size(R)), 1);    % create mask of entries above the diagonal
vals = R(mask);                   % extract those entries into a vector
r = mean(vals);             % compute their mean
%%
keyboard;
% (optional) 
if is_smoothed
    % 1) Convert sigma to samples
    sigma_samp = sigma * Fz;
    
    % 2) Build symmetric Gaussian kernel
    half_k = (kernel_size - 1) / 2;
    x      = -half_k : half_k;                         % sample offsets
    h      = exp(- (x.^2) / (2 * sigma_samp^2));       % unnormalized
    h      = h / sum(h);                               % normalize so sum(h)=1
    
    % 3) Compute beta = sum of squares
    beta = sum(h.^2);
    
else
    beta = 1;
end
fprintf('Beta (sum of h.^2) = %.6f\n', beta);
%
rgc_time_std = 0.061;

equivalent_noise_level = rgc_time_std*(1-r)/(beta*r)

%%
empirical_fac = 0.2;
base_trace = rgc_time_std*randn(1, size(trace, 2));
traces = repmat(base_trace, nRepeats*10, 1);
traces = traces+equivalent_noise_level*empirical_fac*randn(size(traces));
R = corr(traces');

mask = triu(true(size(R)), 1);    % create mask of entries above the diagonal
vals = R(mask);                   % extract those entries into a vector
r_s = mean(vals) 

%%
T    = size(trace, 2);            % number of bins
r_target = r;
initial_fac            = 0.7;   % starting guess
max_bin_count = max(trace, [], 'all');

load('betaParams_52808.mat','alphaHat','betaHat')
pd2 = makedist('Beta','a',alphaHat,'b',betaHat);
% r2 = random(pd2, 5000, 1);
% Run optimization
[fac_opt, err_min] = optimize_empirical_fac_minimize(...
    rgc_time_std, equivalent_noise_level, nRepeats, r_target, T, max_bin_count, pd2, initial_fac...
);


fprintf('Optimized fac = %.4f, minimum error = %.6f\n', fac_opt, err_min);
fprintf('Optimal equivalent_noise_level: %.4f\n', fac_opt*equivalent_noise_level);

% Verify
% base = rgc_time_std * randn(1, T);
base = random(pd2, 1, T);
base = rgc_time_std*base/std(base);
traces = repmat(base, nRepeats*10, 1) + equivalent_noise_level*fac_opt*randn(nRepeats*10, T);
traces = max(traces, 0);
traces = double(traces > mean(traces, 'all'));
% traces = traces*max_bin_count*40;
% traces = poissrnd(traces);
R = corr(traces'); mask = triu(true(size(R)), 1);
r_sim = mean(R(mask));
fprintf('Simulated reliability = %.3f (target %.3f)\n', r_sim, r_target);

%%
keyboard;
%% ---------------------- rectification optimization ----------------------
% This script demonstrates how to call optimize_empirical_params.m
% Define parameters based on your data and simulation settings:
% rgc_time_std          = 0.10;    % example std of smoothed, rectified rgc_time
% equivalent_noise_level = 0.05;    % pre-smoothing noise STD 'x'
% nRepeats             = 20;      % number of repeats to simulate (per trial)
% r_target              = 0.3;     % measured repeat reliability from neuro data
mean_rt               = 0.041;     % mean of smoothed, rectified rgc_time
rect_thr              = 0;     % rectification threshold used in simulation
% T                     = 100;     % length of each time-series (number of bins)
initial_guess         = [0.85, 0.0001];  % initial guesses [fac0, offset0]

% Call the optimizer
defaults = struct(); % ensure path
% [fac_opt, offset_opt] = optimize_empirical_params( ...
%     rgc_time_std, ...
%     equivalent_noise_level, ...
%     nRepeats, ...
%     r_target, ...
%     mean_rt, ...
%     rect_thr, ...
%     T, ...
%     initial_guess ...
% );
fac_range   = linspace(0, 1, 41);         % from 0 to 2 in steps of 0.05
offset_range = linspace(-0.05, 0.05, 41); 
[fac_opt, offset_opt] = optimize_empirical_params_grid(...
    rgc_time_std, equivalent_noise_level, nRepeats, r_target, mean_rt,...
    rect_thr, T, fac_range, offset_range...
);

% Display results
fprintf('Optimized fac = %.4f\n', fac_opt);
fprintf('Optimized offset = %.4f\n', offset_opt);

% Verify with one simulation
nTr = nRepeats*10;
raw = offset_opt + rgc_time_std * randn(nTr, T);
noisy = raw + equivalent_noise_level * fac_opt * randn(size(raw));
rectified = max(noisy, rect_thr);
R = corr(rectified'); mask = triu(true(size(R)), 1);
r_sim = mean(R(mask)); mean_sim = mean(rectified(:));
fprintf('Simulated reliability = %.3f (target %.3f)\n', r_sim, r_target);
fprintf('Simulated mean post-rect = %.3f (target %.3f)\n', mean_sim, mean_rt);


