close all; clear
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
msotMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch1/msot_MB_dynamic_60.mat');

%clarioStarMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch2_same-as-dilution-series/CLARIOstar_MB.mat');
flowSpectrometerMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch1/spectrometer_MB_dynamic.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
%%
wavelengths=msotMethyleneBlue.msot.wavelengths;
msotMethyleneBlue.msot.averagedTimes = mean(msotMethyleneBlue.msot.times,2);
msotDefaultMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

radialProfileDefault = {};
radialProfileFlow = {};
%%
f = waitbar(0,'Please wait...');
for i=1:size(msotMethyleneBlue.msot.recon,3)
    waitbar(i/size(msotMethyleneBlue.msot.recon,3),f,'Spectral unmixing and tube profile extraction');
    allWavelengths = squeeze(msotMethyleneBlue.msot.recon(:,:,i,:));
    load('MethyleneBlue-Dynamic-ROI/tubeROI_b1.mat')

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
    [tf locDefault] = ismember(wavelengths,msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1));
    msotDefaultUnmixing = msotDefaultSpectra.MSOT.MB_in_water_10uM(locDefault,2);
    msotDefaultUnmixing = msotDefaultUnmixing/max(msotDefaultUnmixing);
    [msotDefaultMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[msotDefaultUnmixing]);
    msotDefaultMeans(i)=mean(msotDefaultMSP(tubeRoi));
    
    %find corresponding flowSpectrum
    [c index] = min(abs(flowSpectrometerMethyleneBlue.spec.times-msotMethyleneBlue.msot.averagedTimes(i)));
    
    [tf locFlow] = ismember(wavelengths,round(flowSpectrometerMethyleneBlue.spec.wavs(:,index)));
    flowSpectrometerUnmixing = flowSpectrometerMethyleneBlue.spec.absorbance(locFlow,index);
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
close(f)
%%
timeVector = [zeros(79,1);[1:2000]'];

times = timeVector; % seconds
flowRateIn = 500/60; % uL/min converted to uL/s
circuitVol = 5500; % uL
flowRateOut = 500/60; % uL/min converted to uL/s
concInjected = 429; % uM (umol per L *10^-6 to convert to umol per uL)

estConcDefault = @(I) (I-270.55)/8.0979;
estConcFlow = @(I) (I-236.98)/8.5319;
estConcClario = @(I) (I-262.81)/8.2236;

%Intercep 0
estConcDefault = @(I) I/9.8287;
estConcFlow = @(I) I/10.048;
estConcClario = @(I) I/9.9054;

C3 = (circuitVol*concInjected*10^-6 - circuitVol*concInjected*10^-6*exp(-times.*flowRateIn/circuitVol))./(circuitVol*10^-6);
%%
figure;plot(msotMethyleneBlue.msot.averagedTimes,estConcDefault(msotDefaultMeans),[1:2079]',C3)

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Methylene blue concentration [µM]'; %ylabel
plt.Title = 'MSOT Default spectra'
plt.Legend = {'Actual concentration','Estimated concentration'};
plt.LegendLoc = 'Best';
plt.BoxDim = [5 5];
plt.YLim = [-10 220];
plt.XLim = [0 1800];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
%%plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/MethyleneBlue-Dynamic-ActualEstimated-Default.pdf');
%%
figure;plot(msotMethyleneBlue.msot.averagedTimes,estConcFlow(flowSpectrometerMeans),[1:2079]',C3)

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Methylene blue concentration [µM]'; %ylabel
plt.Title = 'Flow spectrometer spectra'
plt.Legend = {'Actual concentration','Estimated concentration'};
plt.LegendLoc = 'Best';
plt.BoxDim = [5 5];
plt.YLim = [-10 220];
plt.XLim = [0 1800];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
%%plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/MethyleneBlue-Dynamic-ActualEstimated-FlowSpec.pdf');
%% Plot radial profiles
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],106);
figure; hold on;
for j=1:length(radialProfileDefault)
    if mod(j,50)==0
    plot(1:length(radialProfileDefault{j}),estConcDefault(radialProfileDefault{j}));
    end
end
hold off;
colormap(cmap)
colorbar;
caxis([1 1000]) % change here

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Radial tube profile [px]'; % xlabel
plt.YLabel = 'Spectrally unmixed MSOT intensity [a.u.]'; %ylabel
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
plt.export('figure-raw-drafts/MethyleneBlue-Dynamic-Radial-Default.pdf');

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
plt.YLabel = 'Spectrally unmixed MSOT intensity [a.u.]'; %ylabel
plt.Title = 'Flow spec radial profile';
plt.BoxDim = [5 5];
plt.YLim = [-20 350];
plt.XLim = [1 14];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
%plt.LineStyle = lineStyle;
plt.LineWidth = 0.5;
plt.export('figure-raw-drafts/MethyleneBlue-Dynamic-Radial-flowSpec.pdf');


