close all; clear;
%%
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))

%msotBackgroundMB = load('/media/gehrun01/work-io/flow-phantom/data/5_nigrosin_backgrounds/msot_MB_20runs_17wav_4bgd_40uM200uM_sos57opt.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
waterCorrectionSpectum = load('/media/gehrun01/work-io/flow-phantom/data/spectra/WaterCoeffs.mat');
load('/media/gehrun01/work-io/flow-phantom/data/spectra/NdFilter-Spectrum.mat');
nigrosinSpectra = csvread('/media/gehrun01/work-io/flow-phantom/data/spectra/nigrosin.csv');
interpolatedNigrosinSpectra = [round(interp(nigrosinSpectra(:,1),2)), interp(nigrosinSpectra(:,2),2)];


dataCollection = {
    'msot_blood_0_sos55opt','spectrometer_blood_0_005','pO2_blood_0_005',0;
    'msot_blood_0025_sos55opt','spectrometer_blood_0025_01','pO2_blood_0025_01',0.025;
    'msot_blood_005_sos55opt','spectrometer_blood_0_005','pO2_blood_0_005',0.05;
    'msot_blood_01_sos55opt','spectrometer_blood_0025_01','pO2_blood_0025_01',0.1;
    }
%% ND filter spectrum
interpolatedSpectrum = [round(interp(ndFilterTmSpectrum(:,1),5)), interp(ndFilterTmSpectrum(:,2),5)];
%figure;plot(interpolatedSpectrum(:,1),interpolatedSpectrum(:,2))
%%
%

%figure; hold on;
msotSpectra = [];
msotDefaultMeans=[];
clarioStarMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

wavelength660nm = [];

so2Vals = {};

for k=1:size(dataCollection,1)
    flowSpectrometer = load(['/media/gehrun01/work-io/flow-phantom/data/5_nigrosin_backgrounds/',dataCollection{k,2},'.mat']);
    msotData = load(['/media/gehrun01/work-io/flow-phantom/data/5_nigrosin_backgrounds/',dataCollection{k,1},'.mat']);
    %msotData = squeeze(mean(msotBackgroundMB.msot.recon(:,:,:,:,k,i),3));
    load(['BackgroundBloodStudy-ROIs/tubeROI_' num2str(k) '.mat']);
    wavelengths=msotData.msot.wavelengths(:);
    [tf locTm] = ismember(round(flowSpectrometer.spec.wavs(:,1)),interpolatedSpectrum(:,1));
    sameRange = interpolatedSpectrum(locTm,2);
    
    flowHbO2 = [round(flowSpectrometer.spec.wavs(:,1)),mean(flowSpectrometer.spec.absorbance(:,1:5),2)./sameRange];
    flowHb = [round(flowSpectrometer.spec.wavs(:,1)),mean(flowSpectrometer.spec.absorbance(:,end-5:end),2)./sameRange];
    
    [tf locHb] = ismember(wavelengths,flowHb(:,1));
    hbUnmixingFlow = flowHb(locHb,2);

    [tf locHbO2] = ismember(wavelengths,flowHbO2(:,1));
    hbO2UnmixingFlow = flowHbO2(locHbO2,2);
    
    msotDefaultMeansSo2 = [];
    msotDefaultMeansCorrSo2 = [];
    
    correctionSpectrum = interpolatedNigrosinSpectra(:,2)./max(interpolatedNigrosinSpectra(:,2))*dataCollection{k,4};
    
    
    
    [tf locNigrosin] = ismember(wavelengths,interpolatedNigrosinSpectra(:,1));
    nigrosinCorrection = correctionSpectrum(locNigrosin,1);
    for j=1:size(msotData.msot.recon,3)
       cFrameMsot = squeeze(msotData.msot.recon(:,:,j,:));
       cFrameMsotCorr = squeeze(msotData.msot.recon(:,:,j,:));
       for y=1:length(wavelengths)
        cFrameMsotCorr(:,:,y) = cFrameMsotCorr(:,:,y)./exp(-8.5*nigrosinCorrection(y));
       end
       
       [tf locHb] = ismember(wavelengths,msotDefaultSpectra.MSOT.Hb(:,1));
       hbUnmixing = msotDefaultSpectra.MSOT.Hb(locHb,2);
      
       [tf locHbO2] = ismember(wavelengths,msotDefaultSpectra.MSOT.HbO2(:,1));
       hbO2Unmixing = msotDefaultSpectra.MSOT.HbO2(locHbO2,2);
       
       msotDefaultMSP = msp_pinv(cFrameMsot,[hbUnmixing, hbO2Unmixing]);
       msotDefaultMSPSo2 = msotDefaultMSP(:,:,2)./(msotDefaultMSP(:,:,1)+msotDefaultMSP(:,:,2));
       msotDefaultMeansSo2(j) = mean(msotDefaultMSPSo2(tubeRoi));   
       
       msotFlowMSP = msp_pinv(cFrameMsot,[hbUnmixingFlow, hbO2UnmixingFlow]);
       msotFlowMSPSo2 = msotFlowMSP(:,:,2)./(msotFlowMSP(:,:,1)+msotFlowMSP(:,:,2));
       msotFlowMeansSo2(j) = mean(msotFlowMSPSo2(tubeRoi));   
       
       
       msotDefaultMSPCorr = msp_pinv(cFrameMsotCorr,[hbUnmixing, hbO2Unmixing]);
       msotDefaultMSPCorrSo2 = msotDefaultMSPCorr(:,:,2)./(msotDefaultMSPCorr(:,:,1)+msotDefaultMSPCorr(:,:,2));
       msotDefaultMeansCorrSo2(j) = mean(msotDefaultMSPCorrSo2(tubeRoi));    
       
       msotFlowMSPCorr = msp_pinv(cFrameMsotCorr,[hbUnmixingFlow, hbO2UnmixingFlow]);
       msotFlowMSPCorrSo2 = msotFlowMSPCorr(:,:,2)./(msotFlowMSPCorr(:,:,1)+msotFlowMSPCorr(:,:,2));
       msotFlowMeansCorrSo2(j) = mean(msotFlowMSPCorrSo2(tubeRoi));   
       %msotDefaultMSPHb = msotDefaultMSP(:,:,1);
       %msotDefaultMSPHbO2 = msotDefaultMSP(:,:,2);

       %msotDefaultMeansMb(i)=mean(msotDefaultMSPMb(tubeRoi));
       
       
    end
    so2Vals{k,1} = msotDefaultMeansSo2;
    so2Vals{k,2} = msotDefaultMeansCorrSo2;
    so2Vals{k,3} = msotFlowMeansSo2;
    so2Vals{k,4} = msotFlowMeansCorrSo2;
%     figure;imshow(msotData.msot.recon(:,:,1,8),[])
%     h = imellipse;
%     position = wait(h);
%     tubeRoi = h.createMask()
%     save(['BackgroundBloodStudy-ROIs/tubeROI_' num2str(k) '.mat'],'tubeRoi')
end

%%
cmap = lines(4);
%cmap = repmat(cmap,2,1)
figure;hold on;
plot(1:50,so2Vals{1,1}(1:50).*100,'r')
plot(1:50,so2Vals{2,1}(1:50).*100,'c')
plot(1:50,so2Vals{3,1}(1:50).*100,'b')
plot(1:50,so2Vals{4,1}(1:50).*100,'g')
hold off;


plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Flow spectra'
plt.Legend = {'0 cm^{-1}','0.025 cm^{-1}','0.05 cm^{-1}','0.1 cm^{-1}'};
plt.LegendLoc = 'Best';
plt.LegendBox = true;
plt.BoxDim = [8 6];
plt.YLim = [50 100];
plt.XLim = [1 50];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-','-','-'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/BackgroundBlood-UnCorrected.pdf');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
%%
cmap = lines(4);
%cmap = repmat(cmap,2,1)
figure;hold on;
plot(1:50,so2Vals{1,2}(1:50).*100,'r')
plot(1:50,so2Vals{2,2}(1:50).*100,'c')
plot(1:50,so2Vals{3,2}(1:50).*100,'b')
plot(1:50,so2Vals{4,2}(1:50).*100,'g')

hold off;


plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Flow spectra'
%plt.Legend = {'MSOT','Online flow spectrometer','pO₂ probe'};
plt.Legend = {'0 cm^{-1}','0.025 cm^{-1}','0.05 cm^{-1}','0.1 cm^{-1}'};
plt.LegendLoc = 'Best';
plt.LegendBox = true;
plt.BoxDim = [8 6];
plt.YLim = [50 100];
plt.XLim = [1 50];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-.','-.','-.','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/BackgroundBlood-Corrected.pdf');
%%
figure;hold on;
ourbar = bar([1 1;2 2;3 3],[mean(so2Vals{1,1}(1:50).*100),mean(so2Vals{1,2}(1:50).*100);mean(so2Vals{3,1}(1:50).*100),mean(so2Vals{3,2}(1:50).*100);mean(so2Vals{4,1}(1:50).*100),mean(so2Vals{4,2}(1:50).*100)],'BarWidth',0.9)
ourerror = errorbar(0.86,mean(so2Vals{1,1}(1:50).*100),std(so2Vals{1,1}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(1.14,mean(so2Vals{1,2}(1:50).*100),std(so2Vals{1,2}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(1.86,mean(so2Vals{3,1}(1:50).*100),std(so2Vals{3,1}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(2.14,mean(so2Vals{3,2}(1:50).*100),std(so2Vals{3,2}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(2.86,mean(so2Vals{4,1}(1:50).*100),std(so2Vals{4,1}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(3.14,mean(so2Vals{4,2}(1:50).*100),std(so2Vals{4,2}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
hold off;
ourbar(1).FaceColor = [.5 .5 .5]
%ourerror.Color = [0 0 0]

xtickangle(45)
XTick=[0.86,1.14,1.86,2.14,2.86,3.14];
set(gca, 'XTick',XTick);
%set(gca,'TickLabelInterpreter','none')
set(gca, 'XTickLabel', {'0 cm^{-1}','(corrected)','0.05 cm^{-1}','(corrected)','0.1 cm^{-1}','(corrected)'});
line('XData', [0 5], 'YData', [81.7 81.7],'Color','r','LineWidth',2,'LineStyle','--','HandleVisibility','off');


plt = Plot(); % create a Plot object and grab the current figure
%plt.XLabel = 'Tubes (I.D. / O.D. in µm)'; % xlabel
%plt.YLabel = 'Signal-to-background ratio [a.u.]'; %ylabel
plt.XLim = [0.5 3.5]
plt.YLim = [70 95]
plt.BoxDim = [5 5];
plt.LineWidth = 1.5;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.XMinorTick = false;
plt.ShowBox = false;
plt.XLabel = 'Phantoms with different background absorptions'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%set(gca, 'Layer', 'Top');
%rectangle('Position',[7.5 0 1 1.1], 'FaceColor', [0 1 0 0.1],'LineStyle','none');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/BackgroundBlood-Attenuation-Comparison-Default.pdf');

%%
figure;hold on;
ourbar = bar([1 1;2 2;3 3],[mean(so2Vals{1,3}(1:50).*100),mean(so2Vals{1,4}(1:50).*100);mean(so2Vals{3,3}(1:50).*100),mean(so2Vals{3,4}(1:50).*100);mean(so2Vals{4,3}(1:50).*100),mean(so2Vals{4,4}(1:50).*100)],'BarWidth',0.9)
ourerror = errorbar(0.86,mean(so2Vals{1,3}(1:50).*100),std(so2Vals{1,3}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(1.14,mean(so2Vals{1,4}(1:50).*100),std(so2Vals{1,4}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(1.86,mean(so2Vals{3,3}(1:50).*100),std(so2Vals{3,3}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(2.14,mean(so2Vals{3,4}(1:50).*100),std(so2Vals{3,4}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(2.86,mean(so2Vals{4,3}(1:50).*100),std(so2Vals{4,3}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
errorbar(3.14,mean(so2Vals{4,4}(1:50).*100),std(so2Vals{4,4}(1:50).*100),'black','LineWidth',2, 'linestyle', 'none');
hold off;
ourbar(1).FaceColor = [.5 .5 .5]
%ourerror.Color = [0 0 0]

xtickangle(45)
XTick=[0.86,1.14,1.86,2.14,2.86,3.14];
set(gca, 'XTick',XTick);
%set(gca,'TickLabelInterpreter','none')
set(gca, 'XTickLabel', {'0 cm^{-1}','(corrected)','0.05 cm^{-1}','(corrected)','0.1 cm^{-1}','(corrected)'});
line('XData', [0 5], 'YData', [81.7 81.7],'Color','r','LineWidth',2,'LineStyle','--','HandleVisibility','off');


plt = Plot(); % create a Plot object and grab the current figure
%plt.XLabel = 'Tubes (I.D. / O.D. in µm)'; % xlabel
%plt.YLabel = 'Signal-to-background ratio [a.u.]'; %ylabel
plt.XLim = [0.5 3.5]
plt.YLim = [70 140]
plt.BoxDim = [5 5];
plt.LineWidth = 1.5;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.XMinorTick = false;
plt.ShowBox = false;
plt.XLabel = 'Phantoms with different background absorptions'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%set(gca, 'Layer', 'Top');
%rectangle('Position',[7.5 0 1 1.1], 'FaceColor', [0 1 0 0.1],'LineStyle','none');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/BackgroundBlood-Attenuation-Comparison-Online.pdf');
