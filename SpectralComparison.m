close all; clear
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))

msotMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/msot_MB_sos58opt.mat');
msotMetyleneBlueReconMean = squeeze(mean(msotMethyleneBlue.msot.recon,3));

msotIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/msot_ICG_sos67opt.mat');
msotIndocyanineGreenReconMean = squeeze(mean(msotIndocyanineGreen.msot.recon,3));

%clarioStarMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/CLARIOstar_MB.mat');
flowSpectrometerMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/spectrometer_MB.mat');
flowSpectrometerIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/spectrometer_ICG_5uMcorrected.mat');

clarioStarIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/CLARIOstar_ICG.mat');
clarioStarMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/CLARIOstar_MB.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');

%%
lineStyle = {'-'};
figure;
hold on;
plot(msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1),msotDefaultSpectra.MSOT.MB_in_water_10uM(:,2),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorption coefficient [cm^{-1}]'; %ylabel

plt.Legend = {'MB (10μM) literature spectrum','ICG (6.5μM) literature spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 200000];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[133/255 148/255 205/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 200000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-MB-Literature-Spectrum.pdf');
%%
lineStyle = {'-'};
figure;
hold on;
plot(flowSpectrometerMethyleneBlue.spec.wavs(:,5),flowSpectrometerMethyleneBlue.spec.mean(:,5),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'MB (10μM) online spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 0.1];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[133/255 148/255 205/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-MB-DS-Online-Spectrum.pdf');
%%
lineStyle = {'-'};
figure;
hold on;
plot(clarioStarMethyleneBlue.clario.wavs(:),clarioStarMethyleneBlue.clario.mean(:,4),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'MB (10μM) offline spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 1.5];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[133/255 148/255 205/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-MB-Offline-Spectrum.pdf');
%%
allWavelengths = squeeze(msotMetyleneBlueReconMean(:,:,:,5));
wavelengths=[660,664,680,684,694,700,708,715,730,735,760,770,775,779,800,850,950];
errorData = [];
load(['MethyleneBlue-Series-ROIs/tubeROI_MB_DS_5.mat'])
for k=1:length(wavelengths)
    meanAlongReps =allWavelengths(:,:,k);
    wavelengthsDataSlice(k) = mean(meanAlongReps(tubeRoi));
    errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
end
lineStyle = {'-'};
figure;
hold on;
plot(wavelengths,wavelengthsDataSlice,'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'MB (10μM) PAT spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 1000];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[133/255 148/255 205/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-MB-MSOT-Spectrum.pdf');
%%
lineStyle = {'-'};
figure;
hold on;
plot(msotDefaultSpectra.MSOT.ICG(:,1),msotDefaultSpectra.MSOT.ICG(:,2),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorption coefficient [cm^{-1}]'; %ylabel

plt.Legend = {'ICG (6.5μM) literature spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 230000];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[52/255 113/255 71/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-ICG-Literature-Spectrum.pdf');
%%
lineStyle = {'-'};
figure;
hold on;
plot(flowSpectrometerIndocyanineGreen.spec.wavs(:,4),flowSpectrometerIndocyanineGreen.spec.mean(:,4),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'ICG (5μM) online spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 0.06];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[52/255 113/255 71/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-ICG-DS-Online-Spectrum.pdf');
%%
lineStyle = {'-'};
figure;
hold on;
plot(clarioStarIndocyanineGreen.clario.wavs(:),clarioStarIndocyanineGreen.clario.mean(:,4),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'ICG (5μM) offline spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 1.1];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[52/255 113/255 71/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-ICG-Offline-Spectrum.pdf');
%%
allWavelengths = squeeze(msotIndocyanineGreenReconMean(:,:,:,4));
wavelengths=[660,664,680,684,694,700,708,715,730,735,760,770,775,779,800,850,950];
errorData = [];
load(['ICG-Series-ROIs/tubeROI_ICG_DS_4.mat'])
for k=1:length(wavelengths)
    meanAlongReps =allWavelengths(:,:,k);
    wavelengthsDataSlice(k) = mean(meanAlongReps(tubeRoi));
    errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
end
lineStyle = {'-'};
figure;
hold on;
plot(wavelengths,wavelengthsDataSlice,'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel

plt.Legend = {'ICG (5μM) PAT spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 1000];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[52/255 113/255 71/255]};
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 230000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
%text(890,1080,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/SpectralComparison-ICG-MSOT-Spectrum.pdf');