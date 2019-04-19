close all; clear
%%
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
addpath(genpath('./imtool3D-master'))
msotBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/msot_dynamic_sos58opt.mat');

%clarioStarBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/3_MB_ICG_dynamic/batch2_same-as-dilution-series/CLARIOstar_MB.mat');
flowSpectrometerBloodMB = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/spectrometer_dynamic.mat');

partialOxygen = load('/media/gehrun01/work-io/flow-phantom/data/4_blood_dynamic/Experiment1/pO2_dynamic.mat')

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
waterCorrectionSpectum = load('/media/gehrun01/work-io/flow-phantom/data/spectra/WaterCoeffs.mat');
load('/media/gehrun01/work-io/flow-phantom/data/spectra/NdFilter-Spectrum.mat');
 

%%
wavelengths=msotBloodMB.msot.wavelengths(1:16);
msotBloodMB.msot.averagedTimes = mean(msotBloodMB.msot.times,2);
msotDefaultMeansHb=[];
msotDefaultMeansHbO2=[];
msotDefaultMSPs = [];
msotFlowMSPs = [];
msotFlowMeansHb=[];
msotFlowMeansHbO2=[];
radialProfileDefault = {};
radialProfileFlow = {};
rSquaredTube = [];
%%
%% ND filter spectrum
interpolatedSpectrum = [round(interp(ndFilterTmSpectrum(:,1),5)), interp(ndFilterTmSpectrum(:,2),5)];
%figure;plot(interpolatedSpectrum(:,1),interpolatedSpectrum(:,2))

[tf locTm] = ismember(round(flowSpectrometerBloodMB.spec.wavs(:,1)),interpolatedSpectrum(:,1));

sameRange = interpolatedSpectrum(locTm,2);

%%
flowHbO2 = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbance(:,1:5),2)./sameRange];
flowHb = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbance(:,end-5:end),2)./sameRange];
% 
% figure;plot(flowHb(:,1),flowHb(:,2))
% figure;plot(flowHbO2(:,1),flowHbO2(:,2))
%figure;plot(round(flowSpectrometerBloodMB.spec.wavs(:,1)),flowHbO2)
%%
%flowHbO2 = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbanceCorrected(:,1:5),2)];
%flowHb = [round(flowSpectrometerBloodMB.spec.wavs(:,1)),mean(flowSpectrometerBloodMB.spec.absorbanceCorrected(:,end-5:end),2)];

%figure;plot(flowHb(:,1),flowHb(:,2))
%figure;plot(flowHbO2(:,1),flowHbO2(:,2))
%figure;plot(round(flowSpectrometerBloodMB.spec.wavs(:,1)),flowHbO2)
%%
[tf locHb] = ismember(wavelengths,flowHb(:,1));
hbUnmixingFlow = flowHb(locHb,2);
%figure;plot(hbUnmixing)

[tf locHbO2] = ismember(wavelengths,flowHbO2(:,1));
hbO2UnmixingFlow = flowHbO2(locHbO2,2);
%figure;plot(hbO2Unmixing)

[tf locHb] = ismember(wavelengths,msotDefaultSpectra.MSOT.Hb(:,1));
hbUnmixing = msotDefaultSpectra.MSOT.Hb(locHb,2);
% 
%     
[tf locHbO2] = ismember(wavelengths,msotDefaultSpectra.MSOT.HbO2(:,1));
hbO2Unmixing = msotDefaultSpectra.MSOT.HbO2(locHbO2,2);

figure;
plot(wavelengths,hbO2Unmixing./500,'r',wavelengths,hbUnmixing./500,'b');
plots = flipud(get(gca, 'children'));
legendLabels = {'Oxyhemoglobin','Deoxyhemoglobin'};
newOrder = flip(1:2,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorption coefficient [a.u.]'; %yl
plt.BoxDim = [3 2.5];
plt.YLim = [0 7];
plt.XLim = [660 850];
plt.FontName = 'Arial';
plt.FontSize = 12;
plt.ShowBox = false;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
plt.export('figure-raw-drafts/Blood-Dynamic-HbO2Hb-Default.pdf');


figure;
plot(wavelengths,hbO2UnmixingFlow,'r',wavelengths,hbUnmixingFlow,'b');
plots = flipud(get(gca, 'children'));
legendLabels = {'Oxyhemoglobin','Deoxyhemoglobin'};
newOrder = flip(1:2,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorption coefficient [a.u.]'; %ylabel
plt.BoxDim = [3 2.5];
plt.YLim = [0 7];
plt.XLim = [660 850];
plt.FontName = 'Arial';
plt.FontSize = 12;
plt.ShowBox = false;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
plt.export('figure-raw-drafts/Blood-Dynamic-HbO2Hb-Flow.pdf');

%%


f = waitbar(0,'Please wait...');
    load('Blood-Dynamic-ROI/tubeROI_e1_steady.mat')
    %tubeRoi = imerode(tubeRoi,strel('sphere',6));
for i=1:size(msotBloodMB.msot.recon,3)
    waitbar(i/size(msotBloodMB.msot.recon,3),f,'Spectral unmixing and tube profile extraction');
    allWavelengths = squeeze(msotBloodMB.msot.recon(:,:,i,1:16));
    %allWavelengthts = allWavelengths./exp(-0.085*waterCorrectionSpectum.water(1:16,2))

%         figure;imshow(allWavelengths(:,:,17),[])
%     h = imellipse;
%      position = wait(h);
%      tubeRoi = h.createMask()
%    save(['Blood-Dynamic-ROI/tubeROI_e1_steady.mat'],'tubeRoi')
%     
    
%     figure;imshowpair(tubeRoi,allWavelengths(:,:,1))

    errorData = [];
    for k=1:length(wavelengths)
        %allWavelengths(:,:,k) = allWavelengths(:,:,k)./exp(-0.085*waterCorrectionSpectum.water(k,2));
        meanAlongReps = allWavelengths(:,:,k);
        wavelengthsDataSlice(k,i) = mean(meanAlongReps(tubeRoi));
        errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
    end
    
    %figure;
    %errorbar(wavelengths,wavelengthsDataSlice,errorData);
    
    [tf locWater] = ismember(wavelengths,msotDefaultSpectra.MSOT.H2O(:,1));
    waterUnmixing = msotDefaultSpectra.MSOT.H2O(locWater,2);
    %waterUnmixing = waterUnmixing/max(waterUnmixing);
    
    [tf locMb] = ismember(wavelengths,msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1));
    mbUnmixing = msotDefaultSpectra.MSOT.MB_in_water_10uM(locMb,2);
    %mbUnmixing = mbUnmixing/max(mbUnmixing)
    
    [tf locUni] = ismember(wavelengths,msotDefaultSpectra.MSOT.uniform(:,1));
    uniUnmixing = msotDefaultSpectra.MSOT.uniform(locUni,2);
    %uniUnmixing = uniUnmixing/max(uniUnmixing);
    
    
    % Spectral unmixing with MSOT Default spectra

    [msotDefaultMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[hbUnmixing, hbO2Unmixing]);
    [msotFlowMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[hbUnmixingFlow, hbO2UnmixingFlow]);
    %msotDefaultMSP = msp_pinv(allWavelengths,[hbUnmixing, hbO2Unmixing]);
    %,waterUnmixing.*100000
    %msotDefaultMSPs(:,:,:,i) = msotDefaultMSP;
    msotFlowMSPs(:,:,:,i) = msotFlowMSP;
    %rSquaredTube(i,:) = [mean(defR2(tubeRoi))]; 
    
    msotDefaultMSPHb = msotDefaultMSP(:,:,1);
    msotDefaultMSPHbO2 = msotDefaultMSP(:,:,2);
    msotDefaultMSPSo2 = msotDefaultMSP(:,:,2)./(msotDefaultMSP(:,:,1)+msotDefaultMSP(:,:,2));
    msotDefaultMeansHb(i)=mean(msotDefaultMSPHb(tubeRoi));
    msotDefaultMeansHbO2(i)=mean(msotDefaultMSPHbO2(tubeRoi));

    msotFlowMSPHb = msotFlowMSP(:,:,1);
    msotFlowMSPHbO2 = msotFlowMSP(:,:,2);
    msotFlowMSPSo2 = msotFlowMSP(:,:,2)./(msotFlowMSP(:,:,1)+msotFlowMSP(:,:,2));
    msotFlowMeansHb(i)=mean(msotFlowMSPHb(tubeRoi));
    msotFlowMeansHbO2(i)=mean(msotFlowMSPHbO2(tubeRoi));
    
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
    
    [ZrD, R] = radialavg(squeeze(msotFlowMSPSo2(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileFlow{i} = ZrD;
    
end
close(f)
%%



figure;hold on;
plot(msotBloodMB.msot.averagedTimes,msotDefaultMeansHb,msotBloodMB.msot.averagedTimes,msotDefaultMeansHbO2)
plot(msotBloodMB.msot.averagedTimes,msotFlowMeansHb,msotBloodMB.msot.averagedTimes,msotFlowMeansHbO2)
hold off;

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Intensity [a.u.]'; %ylabel
plt.Title = 'Unmix Default spectra'
%plt.Legend = {'Actual concentration','Estimated concentration'};
plt.LegendLoc = 'Best';
plt.BoxDim = [5 5];
%plt.YLim = [-10 220];
%plt.XLim = [0 1800];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
%%plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-.'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
%plt.export('figure-raw-drafts/BloodMB-Dynamic-ActualEstimated-Default.pdf');
%%
unmixings = [];
for n=1:535
nthFlowSpectrum = flowSpectrometerBloodMB.spec.absorbance(:,n)./sameRange;
unmixings(n,:) = msp_pinv(nthFlowSpectrum(556:1067)',[flowHbO2(556:1067,2),flowHb(556:1067,2)]);
end
flowSO2 = unmixings(:,1)./(unmixings(:,1)+unmixings(:,2));
%%
cmap = [145,40,40;
254,204,92;
253,141,60;
240,59,32;
189,0,38]/255;
cmap = lines(4)
cmap(1,:) = [145,40,40]/255
severinghausPreSO2 = calcSeveringhaus(partialOxygen.pO2.prePhantomPO2);
severinghausPostSO2 = calcSeveringhaus(partialOxygen.pO2.postPhantomPO2);

kelmanPreSO2 = calcKelman(partialOxygen.pO2.prePhantomPO2,partialOxygen.pO2.prePhantomTemp);
kelmanPostSO2 = calcKelman(partialOxygen.pO2.postPhantomPO2,partialOxygen.pO2.postPhantomTemp);


figure;hold on;
plot(partialOxygen.pO2.times,mean([severinghausPreSO2,severinghausPostSO2],2).*100)
plot(flowSpectrometerBloodMB.spec.times,flowSO2.*100)
plot(msotBloodMB.msot.averagedTimes,(msotFlowMeansHbO2./(msotFlowMeansHbO2+msotFlowMeansHb)).*100)
plot(msotBloodMB.msot.averagedTimes,(msotDefaultMeansHbO2./(msotDefaultMeansHbO2+msotDefaultMeansHb)).*100)


%plot(partialOxygen.pO2.times,mean([kelmanPreSO2,kelmanPostSO2],2).*100)


hold off;

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Flow spectra'
%plt.Legend = {'sO₂ (MSOT / Literature spectra)','sO₂ (MSOT / Online spectra)','sO₂ (Online spectrometer)','sO₂ (Severinghaus model)'};
plt.Legend = {'sO₂ (Severinghaus model)','sO₂ (Online spectrometer)','sO₂ (MSOT / Online spectra)','sO₂ (MSOT / Literature spectra)',};
plt.LegendLoc = 'Best';
plt.LegendBox = 'on';
plt.BoxDim = [5.5 4.5];
plt.YLim = [0 100];
plt.XLim = [0 560];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-',':','--','--'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Dynamic-MSOT-Flow-pO2-Comp.pdf');

%%
cmap = [145,40,40;
254,204,92;
253,141,60;
240,59,32;
189,0,38]/255;
cmap = lines(4)
cmap(1,:) = [145,40,40]/255
severinghausPreSO2 = calcSeveringhaus(partialOxygen.pO2.prePhantomPO2);
severinghausPostSO2 = calcSeveringhaus(partialOxygen.pO2.postPhantomPO2);

kelmanPreSO2 = calcKelman(partialOxygen.pO2.prePhantomPO2,partialOxygen.pO2.prePhantomTemp);
kelmanPostSO2 = calcKelman(partialOxygen.pO2.postPhantomPO2,partialOxygen.pO2.postPhantomTemp);


figure;hold on;
%plot(partialOxygen.pO2.times,mean([severinghausPreSO2,severinghausPostSO2],2).*100)
plot(flowSpectrometerBloodMB.spec.times,flowSO2.*100)
plot(msotBloodMB.msot.averagedTimes,(msotFlowMeansHbO2./(msotFlowMeansHbO2+msotFlowMeansHb)).*100)
plot(msotBloodMB.msot.averagedTimes,(msotDefaultMeansHbO2./(msotDefaultMeansHbO2+msotDefaultMeansHb)).*100)


%plot(partialOxygen.pO2.times,mean([kelmanPreSO2,kelmanPostSO2],2).*100)


hold off;

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time [s]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Flow spectra'
%plt.Legend = {'sO₂ (MSOT / Literature spectra)','sO₂ (MSOT / Online spectra)','sO₂ (Online spectrometer)','sO₂ (Severinghaus model)'};
plt.Legend = {'sO₂ (Online spectrometer)','sO₂ (MSOT / Online spectra)','sO₂ (MSOT / Literature spectra)',};
plt.LegendLoc = 'SouthWest';
plt.LegendBox = 'on';
plt.BoxDim = [5.5 4.5];
plt.YLim = [0 105];
plt.XLim = [0 560];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = {'-','-',':','--','--'};
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
%line('XData', [79 79], 'YData', [-20 50],'Color',[1 0 0 0.5],'LineWidth',2,'HandleVisibility','off');
%text(77,60,{'Injection', 'Start'},'Color',[1 0.4 0.4],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Dynamic-MSOT-Flow-pO2-Comp-Talk-Version.pdf');

%%
other = msotFlowMeansHbO2./(msotFlowMeansHbO2+msotFlowMeansHb);
other = msotDefaultMeansHbO2./(msotDefaultMeansHbO2+msotDefaultMeansHb);
figure;
scatter(other(1:301),downsample(mean([kelmanPreSO2,kelmanPostSO2],2),3))
%%
cmap = colorGradient([255/255 0/255 0/255],[0/255 0/255 255/255],33);

figure; hold on;
for j=1:length(radialProfileDefault)
    reversedValues = radialProfileDefault{j};
    if mod(j,10)==0
    plot([1:length(radialProfileDefault{j})].*0.0753,radialProfileDefault{j}.*100);
    end
end
hold off;
colormap(cmap)
cba = colorbar;
caxis([1 1000]) % change here


plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Radial tube profile [mm]'; % xlabel
plt.YLabel = 'Blood oxygenation sO₂ [%]'; %ylabel
%plt.Title = 'Clario Star radial profile';
plt.BoxDim = [5 5];
plt.YLim = [45 90];
plt.XLim = [2 14].*0.0753;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
%plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;

rectangle('Position',[0 0 10*0.0753 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[10*0.0753 0 4*0.0753 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
rectangle('Position',[11*0.0753 0 4*0.0753 3300], 'FaceColor', [1 1 1 0.9],'LineStyle','none');
text(5*0.0753,3050,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(11*0.0753,3050,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
text(14.1*0.0753,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/Blood-Dynamic-Radial-Default.pdf');

