close all; clear
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
msotIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/msot_ICG_sos67opt.mat');
msotIndocyanineGreenReconMean = squeeze(mean(msotIndocyanineGreen.msot.recon,3));

clarioStarIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/CLARIOstar_ICG.mat');
flowSpectrometerIndocyanineGreen = load('/media/gehrun01/work-io/flow-phantom/data/2_MB_ICG_dilution_series/spectrometer_ICG.mat');

msotDefaultSpectra = load('/media/gehrun01/work-io/flow-phantom/data/spectra/MSOT_default_spectra.mat');
%%
wavelengths=[660,664,680,684,694,700,708,715,730,735,760,770,775,779,800,850,950];
icgConcentrationList = [0,1.25,2.5,5,10,12.5,18.75,25,37.5,50,75,100];

%figure; hold on;
msotSpectra = [];
msotDefaultMeans=[];
clarioStarMeans=[];
flowSpectrometerMeans=[];
rSquaredTube=[];

radialProfileDefault = {};
radialProfileFlow = {};
radialProfileClario = {};
for i=1:length(icgConcentrationList)
    allWavelengths = squeeze(msotIndocyanineGreenReconMean(:,:,:,i));
    
    
     %figure;imshow(allWavelengths(:,:,17),[])
     %h = imellipse;
     %position = wait(h);
     %tubeRoi = h.createMask()
     %save(['ICG-Series-ROIs/tubeROI_ICG_DS_' num2str(i) '.mat'],'tubeRoi')
 
    load(['ICG-Series-ROIs/tubeROI_ICG_DS_' num2str(i) '.mat'])
    %tubeRoi = imerode(tubeRoi,strel('sphere',6));
%figure;imshowpair(allWavelengths(:,:,15),tubeRoi)

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
    [msotDefaultMSP defR2 aR2 fi] = msp_pinv2(allWavelengths,[msotDefaultUnmixing,waterUnmixing]);
    msotDefaultMeans(i)=mean(msotDefaultMSP(tubeRoi));
    
    %figure;imshow(R2,[0 1]);
    
    % Spectral unmixing with Clario Star
    [tf locClario] = ismember(wavelengths,clarioStarIndocyanineGreen.clario.wavs);
    clarioStarUnmixing = clarioStarIndocyanineGreen.clario.mean(locClario,i);
    clarioStarUnmixing = clarioStarUnmixing/max(clarioStarUnmixing);
    [clarioStarMSP csR2 aR2 fi] = msp_pinv2(allWavelengths,[clarioStarUnmixing,waterUnmixing]);
    clarioStarMeans(i)=mean(clarioStarMSP(tubeRoi));
    
   % figure;imshow(csR2,[0 1]);
   % title(icgConcentrationList(i))
    % Spectral unmixing with Spectrometer
    [tf locFlow] = ismember(wavelengths,round(flowSpectrometerIndocyanineGreen.spec.wavs(:,1)));
    flowSpectrometerUnmixing = flowSpectrometerIndocyanineGreen.spec.mean(locFlow,i);
    flowSpectrometerUnmixing = flowSpectrometerUnmixing/max(flowSpectrometerUnmixing);
    [flowSpectrometerMSP fsR2 aR2 fi] = msp_pinv2(allWavelengths,[flowSpectrometerUnmixing,waterUnmixing]);
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
    [ZrD, R] = radialavg(squeeze(msotDefaultMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength,1)),round(edgeLength/2));
    radialProfileDefault{i} = ZrD;
    [ZrF, R] = radialavg(squeeze(flowSpectrometerMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength,1)),round(edgeLength/2));
    radialProfileFlow{i} = ZrF;
    [ZrC, R] = radialavg(squeeze(clarioStarMSP(topLeftY:topLeftY+edgeLength,topLeftX:topLeftX+edgeLength,1)),round(edgeLength/2));
    radialProfileClario{i} = ZrC;

end
%hold off;
%%
dlm1 = fitlm(icgConcentrationList(1:9),msotDefaultMeans(1:9));
dlm2 = fitlm(icgConcentrationList(1:9),flowSpectrometerMeans(1:9));
dlm3 = fitlm(icgConcentrationList(1:9),clarioStarMeans(1:9));

dlm1i = fitlm(icgConcentrationList(1:9),msotDefaultMeans(1:9),'Intercept',false);
dlm2i = fitlm(icgConcentrationList(1:9),flowSpectrometerMeans(1:9),'Intercept',false);
dlm3i = fitlm(icgConcentrationList(1:9),clarioStarMeans(1:9),'Intercept',false);

cmap = repmat(lines(3),2,1);
figure;hold on;
plot(icgConcentrationList,msotDefaultMeans,'*',icgConcentrationList,flowSpectrometerMeans,'*',icgConcentrationList,clarioStarMeans,'*');
plot(dlm1.predict([0:100]'),'-')
plot(dlm2.predict([0:100]'),'-')
plot(dlm3.predict([0:100]'),'-')
hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Indocyanine green concentration [µM]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAI intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 3300];
%plt.XLim = [500 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.Legend = {'MSOT vendor spectrum','Online flow spectrum','Offline spectrum',['LM vendor spectrum, R²=' num2str(round(dlm1.Rsquared.Adjusted,3))],['LM flow spectrum, R²=' num2str(round(dlm2.Rsquared.Adjusted,3))],['LM offline spectrum, R²=' num2str(round(dlm3.Rsquared.Adjusted,3))]}
plt.LegendLoc = 'SouthEast'
plt.LegendBox = true;
plt.TickLength = [0.01 0.01];
plt.LineWidth = [2.5 2.5 2.5 1.5];
plt.LineStyle = {'-','-','-','-.'};
set(gca, 'Layer', 'Top');
rectangle('Position',[37.5 1700 62.5 1600], 'FaceColor', [0.7 0 0 0.075],'LineStyle','none');
text(65,2000,'Non-linear','Color',[0.7 0.5 0.5],'FontSize',12,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-Unmixing-Comparison.pdf');

%%

%%
dlm1 = fitlm(icgConcentrationList(1:9),msotDefaultMeans(1:9));
dlm2 = fitlm(icgConcentrationList(1:9),flowSpectrometerMeans(1:9));

dlm1i = fitlm(icgConcentrationList(1:9),msotDefaultMeans(1:9),'Intercept',false);
dlm2i = fitlm(icgConcentrationList(1:9),flowSpectrometerMeans(1:9),'Intercept',false);

cmap = repmat(lines(3),2,1);
figure;hold on;
plot(icgConcentrationList,msotDefaultMeans,'*',icgConcentrationList,flowSpectrometerMeans,'*');

hold off;
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Indocyanine green concentration [µM]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAI intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 3300];
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
%rectangle('Position',[37.5 1700 62.5 1600], 'FaceColor', [0.7 0 0 0.075],'LineStyle','none');
%text(65,2000,'Non-linear','Color',[0.7 0.5 0.5],'FontSize',12,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-Unmixing-Comparison-Talk-Version.pdf');



%%
lineStyle = {'-',':','-.','-',':','-.','-',':','-.','-',':','-.'}
figure;hold on;
cmap = colorGradient([160/255 213/255 171/255],[12/255 73/255 31/255],12);
for k=1:length(icgConcentrationList)
    plot(flowSpectrometerIndocyanineGreen.spec.wavs(282:1066,1),flowSpectrometerIndocyanineGreen.spec.mean(282:1066,k),'Color','r')
end

hold off;
%title('Flow Spectrometer')
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 0.65];
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
text(890,0.28,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
%text(863,0.2,'MSR','Color',[0.5 0.5 0.5],'FontSize',14,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-FlowSpectrometer.pdf');
%plt.Title = 'Voltage as a function of time'; % plot title
%%
lineStyle = {'-',':','-.','-',':','-.','-',':','-.','-',':','-.'}
figure;hold on;
cmap = colorGradient([160/255 213/255 171/255],[12/255 73/255 31/255],12);
for k=1:length(icgConcentrationList)
    plot(clarioStarIndocyanineGreen.clario.wavs(101:551,1),clarioStarIndocyanineGreen.clario.mean(101:551,k),'Color','r')
end

hold off;
%title('Flow Spectrometer')
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
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
%text(870,2.0,'MSR','Color',[0.5 0.5 0.5],'FontSize',14,'FontWeight','bold');
line('XData', [670 790], 'YData', [7.1 7.1],'Color','r','LineWidth',2,'HandleVisibility','off');
text(650,7.3,'Spectrometer saturation','Color','r','FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-ClarioStar.pdf');
%% MSOT plot 
lineStyle = {'-',':','-.','-',':','-.','-',':','-.','-',':','-.'}
figure;hold on;
cmap = colorGradient([160/255 213/255 171/255],[12/255 73/255 31/255],12);
for k=1:length(icgConcentrationList)
    plot(wavelengths,wavelengthsDataSlice(:,k),'Color','r');
end

hold off;
%title('Flow Spectrometer')
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'PAI intensity [a.u.]'; %ylabel
plt.BoxDim = [5 5];
plt.YLim = [0 3000];
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
text(890,1285,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-MSOT.pdf');
%plt.Title = 'Voltage as a function of time'; % plot title
%% Plot radial profiles
figure; hold on;
for j=1:length(radialProfileDefault)
    plot([1:length(radialProfileClario{j})].*0.0753,radialProfileDefault{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [mm]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAI intensity [a.u.]'; %ylabel
%plt.Title = 'MSOT default radial profile';
plt.BoxDim = [5 5];
plt.YLim = [0 3700];
plt.XLim = [0.1506 1.2801];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;

rectangle('Position',[0 0 0.750 3700], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3700], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(1.128,1580,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
text(0.350,3550,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,3550,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-Radial-MSOTDefault.pdf');


figure; hold on;
for j=1:length(radialProfileFlow)
    plot([1:length(radialProfileClario{j})].*0.0753,radialProfileFlow{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [mm]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAI intensity [a.u.]'; %ylabel
%plt.Title = 'Flow spec radial profile';
plt.BoxDim = [5 5];
plt.YLim = [0 3400];
plt.XLim = [0.1506 1.2801];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.LineStyle = lineStyle;
plt.LineWidth = 2.5;

rectangle('Position',[0 0 0.750 3400], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3400], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(0.350,3275,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,3275,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
text(1.128,1470,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-Radial-flowSpec.pdf');

figure; hold on;
for j=1:length(radialProfileClario)
    plot([1:length(radialProfileClario{j})].*0.0753,radialProfileClario{j});
end
hold off;
plots = flipud(get(gca, 'children'));
legendLabels = strcat(strsplit(num2str(icgConcentrationList)));
newOrder = flip(1:12,2);
legend(plots(newOrder),legendLabels(newOrder));

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Distance from tube centre [mm]'; % xlabel
plt.YLabel = 'Spectrally unmixed PAI intensity [a.u.]'; %ylabel
%plt.Title = 'Clario Star radial profile';
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

rectangle('Position',[0 0 0.750 3300], 'FaceColor', [0 0 0 0.075],'LineStyle','none');
rectangle('Position',[0.750 0 0.300 3300], 'FaceColor', [0 0 0 0.15],'LineStyle','none');
text(0.350,3150,'Tube lumen','Color',[0.5 0.5 0.5],'FontSize',11,'FontWeight','bold');
text(0.850,3150,'Wall','Color',[1 1 1],'FontSize',11,'FontWeight','bold');
text(1.128,1420,'in µM','Color',[0 0 0],'FontSize',11,'FontWeight','bold');
plt.export('figure-raw-drafts/IndocyanineGreen-DS-Radial-clarioStar.pdf');