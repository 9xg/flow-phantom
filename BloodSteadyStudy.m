close all; clear
%%
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
msotBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/msot_steady_sos58opt.mat');

%clarioStarBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch2_same-as-dilution-series/CLARIOstar_MB.mat');
flowSpectrometerBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/spectrometer_steady.mat');

partialOxygen = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/pO2_steady.mat')

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
waterCorrectionSpectum = load('/media/gehrun01/work-io/flow-phantom/data/spectra/WaterCoeffs.mat');
load('/media/gehrun01/work-io/flow-phantom/data/spectra/NdFilter-Spectrum.mat');

%crop data to first 250 seconds
msotBloodMB.msot.recon = msotBloodMB.msot.recon(:,:,1:148,:)
msotBloodMB.msot.times = msotBloodMB.msot.times(1:148,:)

%%
wavelengths=msotBloodMB.msot.wavelengths(:);
msotBloodMB.msot.averagedTimes = mean(msotBloodMB.msot.times,2);
msotDefaultMeansHb=[];
msotDefaultMeansHbO2=[];

radialProfileDefault={};
radialProfileDefaultMatrix=[];
%% ND filter spectrum
interpolatedSpectrum = [round(interp(ndFilterTmSpectrum(:,1),5)), interp(ndFilterTmSpectrum(:,2),5)];
% %figure;plot(interpolatedSpectrum(:,1),interpolatedSpectrum(:,2))
% 
[tf locTm] = ismember(round(flowSpectrometerBloodMB.spec.wavs(:,1)),interpolatedSpectrum(:,1));
% 
sameRange = interpolatedSpectrum(locTm,2);

%%
flowHbO2 = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbance(:,1:5),2)./sameRange];
flowHb = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbance(:,end-5:end),2)./sameRange];
% 
figure;plot(flowHb(:,1),flowHb(:,2))
figure;plot(flowHbO2(:,1),flowHbO2(:,2))
%figure;plot(round(flowSpectrometerBloodMB.spec.wavs(:,1)),flowHbO2)
%%
%flowHbO2 = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbanceCorrected(:,1:5),2)];
%flowHb = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbanceCorrected(:,end-5:end),2)];

%figure;plot(flowHb(:,1),flowHb(:,2))
%figure;plot(flowHbO2(:,1),flowHbO2(:,2))
%figure;plot(round(flowSpectrometerBloodMB.spec.wavs(:,1)),flowHbO2)
%%
[tf locHb] = ismember(wavelengths,flowHb(:,1));
hbUnmixing = flowHb(locHb,2);
figure;plot(hbUnmixing)

[tf locHbO2] = ismember(wavelengths,flowHbO2(:,1));
hbO2Unmixing = flowHbO2(locHbO2,2);
figure;plot(hbO2Unmixing)
%%


f = waitbar(0,'Please wait...');
    load('Blood-Dynamic-ROI/tubeROI_e1_steady.mat')
    %tubeRoi = imerode(tubeRoi,strel('sphere',6));
for i=1:size(msotBloodMB.msot.recon,3)
    waitbar(i/size(msotBloodMB.msot.recon,3),f,'Spectral unmixing and tube profile extraction');
    allWavelengths = squeeze(msotBloodMB.msot.recon(:,:,i,:));

%         figure;imshow(allWavelengths(:,:,17),[])
%     h = imellipse;
%      position = wait(h);
%      tubeRoi = h.createMask()
%    save(['Blood-Dynamic-ROI/tubeROI_e1_steady.mat'],'tubeRoi')
%     
    
%     figure;imshowpair(tubeRoi,allWavelengths(:,:,1))

    errorData = [];
    for k=1:length(wavelengths)
        meanAlongReps = allWavelengths(:,:,k);
        wavelengthsDataSlice(k,i) = mean(meanAlongReps(tubeRoi));
        errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
    end
    
    
    
    %figure;
    %errorbar(wavelengths,wavelengthsDataSlice,errorData);
     [tf locHb] = ismember(wavelengths,msotDefaultSpectra.MSOT.Hb(:,1));
     hbUnmixing = msotDefaultSpectra.MSOT.Hb(locHb,2);

    
     [tf locHbO2] = ismember(wavelengths,msotDefaultSpectra.MSOT.HbO2(:,1));
     hbO2Unmixing = msotDefaultSpectra.MSOT.HbO2(locHbO2,2);

    
    [tf locWater] = ismember(wavelengths,msotDefaultSpectra.MSOT.H2O(:,1));
    waterUnmixing = msotDefaultSpectra.MSOT.H2O(locWater,2);
    %waterUnmixing = waterUnmixing/max(waterUnmixing);
    
    [tf locMb] = ismember(wavelengths,msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1));
    mbUnmixing = msotDefaultSpectra.MSOT.MB_in_water_10uM(locMb,2);
    %mbUnmixing = mbUnmixing/max(mbUnmixing)
    
    [tf locUni] = ismember(wavelengths,msotDefaultSpectra.MSOT.uniform(:,1));
    uniUnmixing = msotDefaultSpectra.MSOT.uniform(locUni,2);
    uniUnmixing = uniUnmixing/max(uniUnmixing);
    
    
    % Spectral unmixing with MSOT Default spectra
    %[msotDefaultMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[hbUnmixing, hbO2Unmixing]);
    msotDefaultMSP = msp_pinv(allWavelengths,[hbUnmixing, hbO2Unmixing,waterUnmixing]);
    %,waterUnmixing.*100000
    
    msotDefaultMSPHb = msotDefaultMSP(:,:,1);
    msotDefaultMSPHbO2 = msotDefaultMSP(:,:,2);
    msotDefaultMSPSo2 = msotDefaultMSP(:,:,2)./(msotDefaultMSP(:,:,1)+msotDefaultMSP(:,:,2));
    %msotDefaultMSPMb = msotDefaultMSP(:,:,3);
    msotDefaultMeansHb(i)=mean(msotDefaultMSPHb(tubeRoi));
    msotDefaultMeansHbO2(i)=mean(msotDefaultMSPHbO2(tubeRoi));
    %msotDefaultMeansMb(i)=mean(msotDefaultMSPMb(tubeRoi));
    
    msotTubeSpectra = [];
    for h=1:length(wavelengths)
        tmpMaskedTube = allWavelengths(:,:,h);
        msotTubeSpectra(h) = mean(tmpMaskedTube(tubeRoi));
    end
    
    tubeProps = regionprops(tubeRoi,'BoundingBox');
    
    topLeftX = round(tubeProps(1).BoundingBox(1))-7;
    topLeftY = round(tubeProps(1).BoundingBox(2))-7;
    edgeLength = max(round(tubeProps(1).BoundingBox(3:4)))+14;
    [ZrD, R] = radialavg(squeeze(msotDefaultMSPSo2(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileDefault{i} = ZrD;
    radialProfileDefaultMatrix(i,:) = ZrD;
    
    
end
close(f)
%% 1:325 pO2 curve
severinghausPreSO2 = calcSeveringhaus(partialOxygen.pO2.prePhantomPO2(1:325));
severinghausPostSO2 = calcSeveringhaus(partialOxygen.pO2.postPhantomPO2(1:325));

kelmanPreSO2 = calcKelman(partialOxygen.pO2.prePhantomPO2(1:325),partialOxygen.pO2.prePhantomTemp(1:325));
kelmanPostSO2 = calcKelman(partialOxygen.pO2.postPhantomPO2(1:325),partialOxygen.pO2.postPhantomTemp(1:325));

figure;
plot(partialOxygen.pO2.times(1:325),partialOxygen.pO2.prePhantomPO2(1:325),'r',partialOxygen.pO2.times(1:325),partialOxygen.pO2.postPhantomPO2(1:325),'b')

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Partial pressure of pO₂ [mmHg]'; %ylabel
%plt.Title = 'Flow spectra'
%plt.Legend = {'Blood oxygenation',['LM blood oxygenation / R²=' num2str(round(dlm3.Rsquared.Adjusted,3))]};
plt.LegendLoc = 'Best';
plt.BoxDim = [3 2];
plt.YLim = [0 200];
plt.XLim = [0 200];
plt.FontName = 'Arial';
plt.FontSize = 11;
plt.ShowBox = false;
%%plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Steady-pO2-Probe.pdf');
%%
cmap = hot(50);

radialProfileDefaultMean = mean(radialProfileDefaultMatrix,2);

figure; hold on;
for j=1:length(radialProfileDefault)
    if mod(j,10)==0
    plot(1:length(radialProfileDefault{j}),radialProfileDefault{j}.*100);
    end
end
hold off;
plots = flipud(get(gca, 'children'));
%legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
%newOrder = flip(1:12,2);
%legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Radial tube profile [px]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
plt.BoxDim = [4 3];
plt.YLim = [75 101];
plt.XLim = [2 16];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;

rectangle('Position',[0 0 10 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[10 0 4 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(14.1,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
text(3,2700,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(11,2700,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Steady-Radial-Default.pdf');

%%
cmap = [.6,.1,.1;.9,.2,.2;.8,.6,.6];

firstFlowSpectrum = flowSpectrometerBloodMB.spec.absorbance(:,1)./sameRange
unmixings = []
for n=1:196
nthFlowSpectrum = flowSpectrometerBloodMB.spec.absorbance(:,n)./sameRange;
unmixings(n) = msp_pinv(nthFlowSpectrum(556:1067)',[firstFlowSpectrum(556:1067)]);
end
yVals = msotDefaultMeansHbO2./(msotDefaultMeansHbO2+msotDefaultMeansHb);
dlm3 = fitlm(msotBloodMB.msot.averagedTimes,yVals);

figure;hold on;
plot(msotBloodMB.msot.averagedTimes,msotDefaultMeansHbO2./(msotDefaultMeansHbO2+msotDefaultMeansHb).*100);
plot(partialOxygen.pO2.times(1:325),mean([kelmanPreSO2,kelmanPostSO2],2).*100)
plot(flowSpectrometerBloodMB.spec.times(1:196),unmixings.*100);
%plot(dlm3.predict([0:200]').*100,'-');
%plot(partialOxygen.pO2.times(1:325),partialOxygen.pO2.prePhantomPO2(1:325),'r',partialOxygen.pO2.times(1:325),partialOxygen.pO2.postPhantomPO2(1:325),'b')
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = {'MSOT','Online flow spectrometer','pO₂ probe'};
newOrder = flip(1:3,2);
legend(plots(newOrder),legendLabels(newOrder),'Location','best');

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Flow spectra'
%plt.Legend = {'MSOT','Online flow spectrometer','pO₂ probe'};
plt.LegendLoc = 'East';
plt.BoxDim = [4 3];
plt.YLim = [75 101];
plt.XLim = [0 200];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Steady-MSOT-Flow-pO2-Comp.pdf');