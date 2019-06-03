close all; clear
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
msotMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/msot_MB_sos58opt.mat');
msotMetyleneBlueReconMean = squeeze(mean(msotMethyleneBlue.msot.recon,3));

clarioStarMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/CLARIOstar_MB.mat');
flowSpectrometerMethyleneBlue = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/spectrometer_MB.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
%%
dropoutIndices = [2,7,8,9]; % 0µM repeat, all 25µM
msotMethyleneBlue.msot.recon(:,:,:,:,dropoutIndices) = []; 
flowSpectrometerMethyleneBlue.spec.conc(dropoutIndices) = []; 
flowSpectrometerMethyleneBlue.spec.wavs(:,dropoutIndices) = []; 
flowSpectrometerMethyleneBlue.spec.mean(:,dropoutIndices) = []; 
flowSpectrometerMethyleneBlue.spec.std(:,dropoutIndices) = []; 
flowSpectrometerMethyleneBlue.spec.reps(dropoutIndices) = []; 
msotMetyleneBlueReconMean(:,:,:,dropoutIndices) = [];

clarioStarDropout = 6; % 25µM
clarioStarMethyleneBlue.clario.mean(:,clarioStarDropout) = [];
clarioStarMethyleneBlue.clario.std(:,clarioStarDropout) = [];
clarioStarMethyleneBlue.clario.reps(:,clarioStarDropout) = [];
clarioStarMethyleneBlue.clario.conc(:,clarioStarDropout) = [];
%% Adapt baseline by correcting for offset
flowSpectrometerMethyleneBlue.spec.mean(:,2) = flowSpectrometerMethyleneBlue.spec.mean(:,2) - mean(flowSpectrometerMethyleneBlue.spec.mean(1068:1157,2)); 
flowSpectrometerMethyleneBlue.spec.mean(:,10) = flowSpectrometerMethyleneBlue.spec.mean(:,10) - mean(flowSpectrometerMethyleneBlue.spec.mean(1068:1157,10)); 


%%
wavelengths=[660,664,680,684,694,700,708,715,730,735,760,770,775,779,800,850,950];
mbConcentrationList = [0,2.5,5,10,18.75,37.5,50,75,100];
%mbConcentrationList = [0,2.5,5,10,18.75,37.5,50,75,100,175,250,500];
%figure; hold on;
msotSpectra = [];
msotDefaultMeans=[];
clarioStarMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

radialProfileDefault = {};
radialProfileFlow = {};
radialProfileClario = {};
% figure; hold on;
% subplotCounter = 1;
% subplotM = 13;
% subplotN = 3;
%figure; hold on;
for i=1:length(mbConcentrationList)
    allWavelengths = squeeze(msotMetyleneBlueReconMean(:,:,:,i));
    %figure;imshow(allWavelengths(:,:,17),[])
    %h = imellipse;
    %position = wait(h);
    %tubeRoi = h.createMask()
    %save(['MethyleneBlue-Series-ROIs/tubeROI_MB_DS_' num2str(i) '.mat'],'tubeRoi')

    load(['MethyleneBlue-Series-ROIs/tubeROI_MB_DS_' num2str(i) '.mat'])
    %figure;imshow(tubeRoi)
    %tubeRoi = imerode(tubeRoi,strel('sphere',6));
    %figure;imshow(tubeRoi)
    %return;
    %close all
    %figure;imshowpair(allWavelengths(:,:,14),tubeRoi);
    %continue;
    %wavelengthsData = [];
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
    
    %figure;imshow(R2,[0 1]);
    
    % Spectral unmixing with Clario Star
    [tf locClario] = ismember(wavelengths,clarioStarMethyleneBlue.clario.wavs);
    clarioStarUnmixing = clarioStarMethyleneBlue.clario.mean(locClario,i);
    clarioStarUnmixing = clarioStarUnmixing/max(clarioStarUnmixing);
    [clarioStarMSP csR2 aR2 fi] = msp_pinv2(allWavelengths,[clarioStarUnmixing]);
    clarioStarMeans(i)=mean(clarioStarMSP(tubeRoi));
    
   % figure;imshow(csR2,[0 1]);
   % title(mbConcentrationList(i))
    % Spectral unmixing with Spectrometer
    [tf locFlow] = ismember(wavelengths,round(flowSpectrometerMethyleneBlue.spec.wavs(:,1)));
    flowSpectrometerUnmixing = flowSpectrometerMethyleneBlue.spec.mean(locFlow,i);
    flowSpectrometerUnmixing = flowSpectrometerUnmixing/max(flowSpectrometerUnmixing);
    [flowSpectrometerMSP fsR2 aR2 fi] = msp_pinv2(allWavelengths,[flowSpectrometerUnmixing]);
    flowSpectrometerMeans(i)=mean(flowSpectrometerMSP(tubeRoi));
    
    %figure;imshow(R2,[0 1]);
    rSquaredTube(i,:) = [mean(defR2(tubeRoi)),mean(csR2(tubeRoi)),mean(fsR2(tubeRoi))]; 

    msotTubeSpectra = [];
    for h=1:length(wavelengths)
        tmpMaskedTube = allWavelengths(:,:,h);
        msotTubeSpectra(h) = mean(tmpMaskedTube(tubeRoi));
    end
    
    tubeProps = regionprops(tubeRoi,'BoundingBox');
    
    topLeftX = round(tubeProps(1).BoundingBox(1))-7;
    topLeftY = round(tubeProps(1).BoundingBox(2))-7;
    edgeLength = max(round(tubeProps(1).BoundingBox(3:4)))+14;
    [ZrD, R] = radialavg(squeeze(msotDefaultMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileDefault{i} = ZrD;
    [ZrF, R] = radialavg(squeeze(flowSpectrometerMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileFlow{i} = ZrF;
    [ZrC, R] = radialavg(squeeze(clarioStarMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength)),round(edgeLength/2));
    radialProfileClario{i} = ZrC;
    
    
    
    %plot(wavelengths,msotTubeSpectra);
%     [jj,ii] = ind2sub([subplotN,subplotM],subplotCounter);
%     i_rowwise1 = sub2ind([subplotM,subplotN],ii,jj);
%     
%     [jj,ii] = ind2sub([subplotN,subplotM],subplotCounter+1);
%     i_rowwise2 = sub2ind([subplotM,subplotN],ii,jj);
%     
%     [jj,ii] = ind2sub([subplotN,subplotM],subplotCounter+2);
%     i_rowwise3 = sub2ind([subplotM,subplotN],ii,jj);
%     
%     subplot(3,13,i_rowwise1);plot(wavelengths,msotTubeSpectra);
%     title(num2str(mbConcentrationList(i)))
%     subplot(3,13,i_rowwise2);plot(wavelengths,clarioStarMethyleneBlue.clario.mean(locClario,i));
%     subplot(3,13,i_rowwise3);plot(wavelengths,flowSpectrometerMethyleneBlue.spec.mean(locFlow,i));
%     subplotCounter = subplotCounter + 3;
    
%     
%     errorData = [];
%     %extract mean from wavelengths
%     for k=1:length(wavelengths)
%         meanAlongReps = mean(Recon(:,:,:,k),3);
%         wavelengthsData(k) = mean(meanAlongReps(tubeRoi));
%         errorData(k) = std(meanAlongReps(tubeRoi))/sqrt(length(meanAlongReps(tubeRoi)));
%     end
%     %averageInRoi(i) = mean2(mean(Recon(:,:,10,1),3))
%     %figure;errorbar(wavelengths,wavelengthsData,errorData);
%     
%     
%     wavelengthsI = [round(interp(AvgofReps(:,1),2)),interp(AvgofReps(:,i),2)];
%     [tf loc] = ismember(wavelengths,wavelengthsI(:,1));
%     icgUnmixing = wavelengthsI(loc,2)
%     icgUnmixing = icgUnmixing/max(icgUnmixing);
    % 
%     figure;imshow(MSP,[]);colormap(gca,jet)
%     mspVal(i)=mean(MSP(tubeRoi));

end
%hold off;
%%

dlm1 = fitlm(mbConcentrationList(1:9),msotDefaultMeans(1:9));
dlm2 = fitlm(mbConcentrationList(1:9),flowSpectrometerMeans(1:9));
dlm3 = fitlm(mbConcentrationList(1:9),clarioStarMeans(1:9));

dlm1i = fitlm(mbConcentrationList(1:9),msotDefaultMeans(1:9),'Intercept',false);
dlm2i = fitlm(mbConcentrationList(1:9),flowSpectrometerMeans(1:9),'Intercept',false);
dlm3i = fitlm(mbConcentrationList(1:9),clarioStarMeans(1:9),'Intercept',false);



cmap = lines(3);
cmap = [cmap;cmap]
figure;hold on;
plot(mbConcentrationList,msotDefaultMeans,'*',mbConcentrationList,flowSpectrometerMeans,'*',mbConcentrationList,clarioStarMeans,'*');
plot(dlm1.predict([1:100]'),'-')
plot(dlm2.predict([1:100]'),'-')
plot(dlm3.predict([1:100]'),'-')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'MB concentration [µM]'; % xlabel
plt.YLabel = 'Spectrally unmixed MB PA signal intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 1300];
%plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.Legend = {'Literature spectrum','Online flow spectrum','Offline spectrum',['LM literature spectrum, R²=' num2str(round(dlm1.Rsquared.Adjusted,3))],['LM online spectrum, R²=' num2str(round(dlm2.Rsquared.Adjusted,3))],['LM offline spectrum, R²=' num2str(round(dlm3.Rsquared.Adjusted,3))]}
plt.LegendLoc = 'SouthEast'
plt.LegendBox = true;
plt.TickLength = [0.01 0.01];
plt.LineWidth = [2.5 2.5 2.5 1.5];
plt.LineStyle = {'none','none','none','-'};
set(gca, 'Layer', 'Top');
%rectangle('Position',[250 2000 250 1200], 'FaceColor', [0.7 0 0 0.075],'LineStyle','none');
%text(330,2150,'Non-linear','Color',[0.7 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-Unmixing-Comparison.pdf');

%%

dlm1 = fitlm(mbConcentrationList(1:9),msotDefaultMeans(1:9));
dlm2 = fitlm(mbConcentrationList(1:9),flowSpectrometerMeans(1:9));
dlm3 = fitlm(mbConcentrationList(1:9),clarioStarMeans(1:9));

dlm1i = fitlm(mbConcentrationList(1:9),msotDefaultMeans(1:9),'Intercept',false);
dlm2i = fitlm(mbConcentrationList(1:9),flowSpectrometerMeans(1:9),'Intercept',false);
dlm3i = fitlm(mbConcentrationList(1:9),clarioStarMeans(1:9),'Intercept',false);



cmap = lines(3);
cmap = [cmap;cmap]
figure;hold on;
plot(mbConcentrationList,msotDefaultMeans,'*',mbConcentrationList,flowSpectrometerMeans,'*');

hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'MB concentration [µM]'; % xlabel
plt.YLabel = 'Spectrally unmixed MB PA signal intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 1300];
%plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.Legend = {'MSOT (literature) spectrum','Online flow spectrum','Offline spectrum',['LM vendor spectrum, R²=' num2str(round(dlm1.Rsquared.Adjusted,3))],['LM flow spectrum, R²=' num2str(round(dlm2.Rsquared.Adjusted,3))],['LM offline spectrum, R²=' num2str(round(dlm3.Rsquared.Adjusted,3))]}
plt.LegendLoc = 'SouthEast'
plt.LegendBox = true;
plt.TickLength = [0.01 0.01];
plt.LineWidth = [2.5 2.5 2.5 1.5];
plt.LineStyle = {'-','-','-','-.'};
set(gca, 'Layer', 'Top');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-Unmixing-Comparison-Talk-Version.pdf');

%%
lineStyle = {'-',':','-.','-',':','-.','-',':','-.'}
figure;hold on;
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],12);
for k=1:length(mbConcentrationList)
    plot(flowSpectrometerMethyleneBlue.spec.wavs(282:1066,1),flowSpectrometerMethyleneBlue.spec.mean(282:1066,k),'Color','r')
end

hold off;
%title('Flow Spectrometer')
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Online spectrometer absorbance [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 0.4];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 1.2], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
text(890,0.385,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(813,0.3,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-FlowSpectrometer.pdf');
%plt.Title = 'Voltage as a function of time'; % plot title
%%
lineStyle = {'-',':','-.','-',':','-.','-',':','-.'}
figure;hold on;
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],12);
for k=1:length(mbConcentrationList)
    plot(clarioStarMethyleneBlue.clario.wavs(101:551,1),clarioStarMethyleneBlue.clario.mean(101:551,k),'Color','r')
end


hold off;
%title('Flow Spectrometer')
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 7.1];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01];
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 7.1], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
text(890,3.0,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(820,2.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
line('XData', [550 700], 'YData', [7.1 7.1],'Color','r','LineWidth',2,'HandleVisibility','off');
text(540,7.3,'Spectrometer saturation','Color','r','FontSize',11,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-ClarioStar.pdf');
%% MSOT plot 
lineStyle = {'none','-','none',':','none','-.','none','-','none',':','none','-.','none','-','none',':','none','-.','none','-','none',':','none','-.'};
figure;hold on;
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],12);
cmap = reshape(repmat(cmap(:)',2,[]),[],3);
shownHandles = [];
for k=1:length(mbConcentrationList)
    f = fit(wavelengths',wavelengthsDataSlice(:,k),'smoothingspline');
    plot(wavelengths,wavelengthsDataSlice(:,k),'*');
    shownHandles(end+1) = plot(660:950,f(660:950));
end

hold off;
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(shownHandles(newOrder),legendLabels(newOrder))

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'PA signal intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 2500];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 660 3000], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
text(890,1650,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(840,700.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',12,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-MSOT.pdf');
%plt.Title = 'Voltage as a function of time'; % plot title
%%
lineStyle = {'-',':','-.','-',':','-.','-',':','-.','-',':','-.'};
cmap = colorGradient([183/255 198/255 255/255],[23/255 39/255 115/255],12);
figure;
hold on;
plot(msotDefaultSpectra.MSOT.MB_in_water_10uM(:,1),msotDefaultSpectra.MSOT.MB_in_water_10uM(:,2),'Color','b')
plot(msotDefaultSpectra.MSOT.ICG(:,1),msotDefaultSpectra.MSOT.ICG(:,2),'Color','b')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorption coefficient [cm^{-1}]'; %ylabel

plt.Legend = {'MB (10μM) vendor spectrum','ICG (6.5μM) vendor spectrum'}
plt.LegendLoc = 'NorthWest'
plt.LegendBox = true;

plt.BoxDim = [5 5];
plt.YLim = [0 200000];
plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = {[133/255 148/255 205/255],[52/255 113/255 71/255]};
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
plt.export('figure-raw-drafts/MB-ICG-Default-Spectrum.pdf');

%% Plot radial profiles
figure; hold on;
for j=1:length(radialProfileDefault)
    plot([1:length(radialProfileDefault{j})].*0.0753,radialProfileDefault{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [µm]'; % xlabel
plt.YLabel = 'Spectrally unmixed MB PA signal intensity [a.u.]'; %ylabel
%plt.Title = 'MSOT default radial profile';
plt.BoxDim = [5 5];
plt.YLim = [0 3300];
plt.XLim = [0.1506 1.2801];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 0.750 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(1.128,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
text(0.350,3150,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,3150,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-Radial-MSOTDefault.pdf');


figure; hold on;
for j=1:length(radialProfileFlow)
    plot([1:length(radialProfileFlow{j})].*0.0753,radialProfileFlow{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [mm]'; % xlabel
plt.YLabel = 'Spectrally unmixed MB PA signal intensity [a.u.]'; %ylabel
%plt.Title = 'Flow spec radial profile';
plt.BoxDim = [5 5];
plt.YLim = [0 3300];
plt.XLim = [0.1506 1.2801];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 0.750 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(0.350,2850,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,2850,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
text(1.128,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-Radial-flowSpec.pdf');

figure; hold on;
for j=1:length(radialProfileClario)
    plot([1:length(radialProfileClario{j})].*0.0753,radialProfileClario{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(mbConcentrationList)));
newOrder = flip(1:9,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [mm]'; % xlabel
plt.YLabel = 'Spectrally unmixed MB PA signal intensity [a.u.]'; %ylabel
%plt.Title = 'Clario Star radial profile';
plt.BoxDim = [5 5];
plt.YLim = [0 3300];
plt.XLim = [0.1506 1.280];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;
set(gca, 'Layer', 'Top');
rectangle('Position',[0 0 0.750 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(0.350,3050,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,3050,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
text(1.128,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
set(get(gca,'XLabel'),'FontWeight','bold','FontSize',15)
set(get(gca,'YLabel'),'FontWeight','bold','FontSize',15)
plt.export('figure-raw-drafts/MethyleneBlue-DS-Radial-clarioStar.pdf');