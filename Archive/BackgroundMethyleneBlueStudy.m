close all; clear
%%
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))

msotBackgroundMB = load('/media/gehrun01/work-io/flow-phantom/data/5_nigrosin_backgrounds/msot_MB_20runs_17wav_4bgd_40uM200uM_sos57opt.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
waterCorrectionSpectum = load('/media/gehrun01/work-io/flow-phantom/data/spectra/WaterCoeffs.mat');
load('/media/gehrun01/work-io/flow-phantom/data/spectra/NdFilter-Spectrum.mat');

dataCollection = {
    'spectrometer_MB_0_005','pO2_MB_0_005';
    'spectrometer_MB_0025_01','pO2_MB_0025_01';
    'spectrometer_MB_0_005','pO2_MB_0_005';
    'spectrometer_MB_0025_01','pO2_MB_0025_01';
    }
%%
wavelengths=msotBackgroundMB.msot.wavelengths(:);

%figure; hold on;
msotSpectra = [];
msotDefaultMeans=[];
clarioStarMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

wavelength660nm = []

for k=1:size(dataCollection,1)
    flowSpectrometer = load(['/media/gehrun01/work-io/flow-phantom/data/5_nigrosin_backgrounds/',dataCollection{k,1},'.mat']);
    for i=1:2
        msotData = squeeze(mean(msotBackgroundMB.msot.recon(:,:,:,:,k,i),3));
        load(['BackgroundMethyleneBlueStudy-ROIs/tubeROI_' num2str(k) '_' num2str(i) '.mat']);
        
        % Spectral unmixing with MSOT Default spectra
        [tf locDefault] = ismember(wavelengths,msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1));
        msotDefaultUnmixing = msotDefaultSpectra.MSOT.MB_in_water_10uM(locDefault,2);
        msotDefaultUnmixing = msotDefaultUnmixing/max(msotDefaultUnmixing);
        [msotDefaultMSP defR2 aR2 fi] = msp_pinv2(msotData,[msotDefaultUnmixing]);
        msotDefaultMeans(k,i)=mean(msotDefaultMSP(tubeRoi));
        
        at660 = msotData(:,:,1);
        
        wavelength660nm(k,i) = mean(at660(tubeRoi));
        
        %         figure;imshow(msotData(:,:,8),[])
        %         h = imellipse;
        %         position = wait(h);
        %         tubeRoi = h.createMask()
        %         save(['BackgroundMethyleneBlueStudy-ROIs/tubeROI_' num2str(k) '_' num2str(i) '.mat'],'tubeRoi')
    end
end



%%
