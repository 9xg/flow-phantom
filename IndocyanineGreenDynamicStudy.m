close all; clear
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
msotIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch2_same-as-dilution-series/msot_ICG_dynamic_sos67opt.mat');

%clarioStarIndocyanineGreen = load('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/Flow/2018-06-05,06,07_Dilution_series/CLARIOstar_ICG.mat');
flowSpectrometerIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch2_same-as-dilution-series/spectrometer_ICG_dynamic.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
%%
wavelengths=msotIndocyanineGreen.msot.wavelengths;
msotIndocyanineGreen.msot.averagedTimes = mean(msotIndocyanineGreen.msot.times,2);
msotDefaultMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

radialProfileDefault = {};
radialProfileFlow = {};
%%
for i=1:size(msotIndocyanineGreen.msot.recon,3)
    allWavelengths = squeeze(msotIndocyanineGreen.msot.recon(:,:,i,:));
    load('IndocyanineGreen-Dynamic-ROI/tubeROI.mat')

    errorData = [];
    for k=1:length(wavelengths)
        meanAlongReps =allWavelengths(:,:,k);
        wavelengthsDataSlice(k,i) = mean(meanAlongReps(tubeRoi));
        errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
    end
    
    
    
    %figure;
    %errorbar(wavelengths,wavelengthsDataSlice,errorData);
    [tf locWater] = ismember(wavelengths,msotDefaultSpectra.MSOT.H2O(:,1));
    waterUnmixing = msotDefaultSpectra.MSOT.H2O(locWater,2);
    waterUnmixing = waterUnmixing/max(waterUnmixing);
    
    [tf locUniform] = ismember(wavelengths,msotDefaultSpectra.MSOT.uniform(:,1));
    uniformUnmixing = msotDefaultSpectra.MSOT.uniform(locUniform,2);
    uniformUnmixing = uniformUnmixing/max(uniformUnmixing);
    
    % Spectral unmixing with MSOT Default spectra
    [tf locDefault] = ismember(wavelengths,msotDefaultSpectra.MSOT.ICG(:,1));
    msotDefaultUnmixing = msotDefaultSpectra.MSOT.ICG(locDefault,2);
    msotDefaultUnmixing = msotDefaultUnmixing/max(msotDefaultUnmixing);
    [msotDefaultMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[msotDefaultUnmixing]);
    msotDefaultMeans(i)=mean(msotDefaultMSP(tubeRoi));
    
    %find corresponding flowSpectrum
    [c index] = min(abs(flowSpectrometerIndocyanineGreen.spec.times-msotIndocyanineGreen.msot.averagedTimes(i)));
    
    [tf locFlow] = ismember(wavelengths,round(flowSpectrometerIndocyanineGreen.spec.wavs(:,index)));
    flowSpectrometerUnmixing = flowSpectrometerIndocyanineGreen.spec.absorbance(locFlow,index);
    flowSpectrometerUnmixing = flowSpectrometerUnmixing/max(flowSpectrometerUnmixing);
    [flowSpectrometerMSP fsR2 aR2 fi] = msp_pinv2(allWavelengths,[flowSpectrometerUnmixing]);
    flowSpectrometerMeans(i)=mean(flowSpectrometerMSP(tubeRoi));
    

    msotTubeSpectra = [];
    for h=1:length(wavelengths)
        tmpMaskedTube = allWavelengths(:,:,h);
        msotTubeSpectra(h) = mean(tmpMaskedTube(tubeRoi));
    end
    
    tubeProps = regionprops(tubeRoi,'BoundingBox');
    
    topLeftX = round(tubeProps(1).BoundingBox(1));
    topLeftY = round(tubeProps(1).BoundingBox(2));
    edgeLength = max(round(tubeProps(1).BoundingBox(3:4)));
    [ZrD, R] = radialavg(squeeze(msotDefaultMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileDefault{i} = ZrD;
    [ZrF, R] = radialavg(squeeze(flowSpectrometerMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileFlow{i} = ZrF;

end
%%
timeVector = [zeros(90,1);[1:1395]'];

times = timeVector; % seconds
flowRateIn = 100/60; % uL/min converted to uL/s
circuitVol = 5500; % uL
flowRateOut = 100/60; % uL/min converted to uL/s
concInjected = 100; % uM (umol per L *10^-6 to convert to umol per uL)

estConcDefault = @(I) (I-385.71)/34.814;
estConcFlow = @(I) (I-573.27)/33.121;
estConcClario = @(I) (I-640.17)/33.196;

%intercept 0
estConcDefault = @(I) I/43.299;
estConcFlow = @(I) I/45.731;
estConcClario = @(I) I/47.278;

C3 = (circuitVol*concInjected*10^-6 - circuitVol*concInjected*10^-6*exp(-times.*flowRateIn/circuitVol))./(circuitVol*10^-6);
%%
localCmap = lines(3)
localCmap(1,:) = [52/255 113/255 71/255]
%figure;plot(msotIndocyanineGreen.msot.averagedTimes(1:end-1),estConcDefault(msotDefaultMeans))
figure;hold on;
plot([1:1485]',C3,msotIndocyanineGreen.msot.averagedTimes,estConcDefault(msotDefaultMeans)-6.6346)
plot(msotIndocyanineGreen.msot.averagedTimes,estConcFlow(flowSpectrometerMeans)-1.0069)
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Indocyanine green concentration [µM]'; %ylabel
%plt.Title = 'MSOT Default spectra'
plt.Legend = {'Modelled concentration','Measured concentration (literature spectrum)','Measured concentration (online spectrum)'};
plt.LegendLoc = 'NorthWest';
plt.LegendBox = true;
plt.BoxDim = [5 5];
plt.YLim = [-5 70];
plt.XLim = [0 1483];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(localCmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.','--'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
line('XData', [79 79], 'YData', [-20 15],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
text(77,18,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-Dynamic-ActualEstimated.pdf');
%%
figure;plot(msotIndocyanineGreen.msot.averagedTimes(1:end),estConcFlow(flowSpectrometerMeans),[1:1485]',C3)

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Indocyanine green concentration [µM]'; %ylabel
plt.Title = 'Flow spectrometer spectra'
plt.Legend = {'Actual concentration','Estimated concentration'};
plt.LegendLoc = 'Best';
plt.BoxDim = [5 5];
plt.YLim = [-20 70];
plt.XLim = [0 1500];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
%%plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
line('XData', [90 90], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
text(77,55,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-Dynamic-ActualEstimated-FlowSpec.pdf');
%% Plot radial profiles
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],106);
figure; hold on;
for j=1:length(radialProfileDefault)
    if mod(j,10)==0
    plot(1:length(radialProfileDefault{j}),estConcDefault(radialProfileDefault{j}));
    end
end
hold off;
colormap(cmap)
colorbar;
caxis([1 1000]) % change here

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Radial tube profile [px]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAT intensity [a.u.]'; %ylabel
plt.Title = 'MSOT default radial profile';
plt.BoxDim = [5 5];
plt.YLim = [-20 350];
plt.XLim = [1 14];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 0.5;
plt.export('figure-raw-drafts/IndocyanineGreen-Dynamic-Radial-Default.pdf');

figure; hold on;
for j=1:length(radialProfileFlow)
    if mod(j,10)==0
    plot(1:length(radialProfileFlow{j}),estConcFlow(radialProfileFlow{j}));
    end
end
hold off;
colormap(cmap)
colorbar;
caxis([1 1000]) % change here

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Radial tube profile [px]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAT intensity [a.u.]'; %ylabel
plt.Title = 'Flow spec radial profile';
plt.BoxDim = [5 5];
plt.YLim = [-20 100];
plt.XLim = [1 14];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
%plt.LineStyle = lineStyle;
plt.LineWidth = 0.5;
plt.export('figure-raw-drafts/IndocyanineGreen-Dynamic-Radial-flowSpec.pdf');


