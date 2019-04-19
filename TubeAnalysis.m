close all; clear
%%
opengl('save', 'software');
addpath(genpath('/media/gehrun01/work-io/cruk-phd-data/flow-phantom/com.itheramedical.msotlib_beta_rev629'))
addpath(genpath('/home/gehrun01/Documents/MATLAB'))
cmap = jet(11)

%excluded '562_750',
tubes = {'300_630','375_500','432_865','500_600','580_960','630_1190','667_1000','1500_2100','1570_2410','2660_2800','2800_3150'};
%tubes = {'630_1190','667_1000','1500_2100','1570_2410','2660_2800','2800_3150'};
%tubes = {'375_500'};
[optimizer, metric] = imregconfig('monomodal');
%figure; hold on;
allTubeSpectra = [];
allTubeSBR = [];
allDemoImgs = []
for k=1:length(tubes)
    msotTube = load(['/media/gehrun01/work-io/flow-phantom/data/1_tubes/Experiment_2/msot_' tubes{k}]);
    wavelengths=msotTube.msot.wavelengths(:);
    firstSlice = squeeze(msotTube.msot.recon(:,:,:,1,:));
    for l=2:size(firstSlice,3)
        tform = imregtform(firstSlice(:,:,l,10), firstSlice(:,:,1,10), 'affine', optimizer, metric);
        for n=1:size(firstSlice,4)
            firstSlice(:,:,l,n) = imwarp(firstSlice(:,:,l,n),tform,'OutputView',imref2d(size(firstSlice(:,:,1,10))));
        end
    end
    
    newSlice = squeeze(mean(firstSlice,3));
    load(['TubeAnalysis-ROIs/tubeROI_' tubes{k} '.mat'])  
    tubeRoiBG = load(['TubeAnalysis-ROIs/tubeROI_' tubes{k} '_BG.mat'])  
    
    allDemoImgs(:,:,k) = squeeze(newSlice(:,:,13))
    
    msotTubeSpectra = [];
    msotTubeSBR = []
    for h=1:length(wavelengths)
        tmpMaskedTube = newSlice(:,:,h);
        msotTubeSpectra(h) = mean(tmpMaskedTube(tubeRoi));
        onlyPhantom = xor(tubeRoiBG.tubeRoi,tubeRoi);
        msotTubeSBR(h) = mean(tmpMaskedTube(tubeRoi))/mean(tmpMaskedTube(onlyPhantom));
    end
    allTubeSpectra(k,:)=msotTubeSpectra;
    allTubeSBR(k,:)=msotTubeSBR;

%      figure;imshow(mean(newSlice(:,:,17),3),[])
%      title('Mean');
%      h = imellipse;
%      position = wait(h);
%     tubeRoi = h.createMask()
%      save(['TubeAnalysis-ROIs/tubeROI_' tubes{k} '_BG.mat'],'tubeRoi')   
end 
%% Clariostar stuff
figure;
plot(wavelengths,allTubeSBR)
colormap(gca,jet)
legend(tubes,'Interpreter','none')

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Signal-to-Background ratios [a.u.]'; %ylabel
plt.BoxDim = [5 5];
%plt.YLim = [0 7.1];
plt.XLim = [660 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.LineWidth = 2.5;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.XMinorTick = false;
plt.ShowBox = false;


set(gca, 'Layer', 'Top');
rectangle('Position',[7.5 0 1 1.1], 'FaceColor', [0 1 0 0.1],'LineStyle','none');
plt.export('figure-raw-drafts/TubeAnalysis-MSOT-SBR-Wavelengths.pdf');
%%
cmap = jet(10)

figure;
plot(wavelengths,allTubeSpectra)
colormap(gca,jet)
legend(tubes,'Interpreter','none')

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Wavelength [nm]'; % xlabel
plt.YLabel = 'Absorbance [a.u.]'; %ylabel
plt.BoxDim = [5 5];
%plt.YLim = [0 7.1];
plt.XLim = [660 950];
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.ShowBox = false;
plt.Colors = num2cell(cmap,2);
plt.LineWidth = 2.5;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.XMinorTick = false;
plt.ShowBox = false;


set(gca, 'Layer', 'Top');
rectangle('Position',[7.5 0 1 1.1], 'FaceColor', [0 1 0 0.1],'LineStyle','none');
plt.export('figure-raw-drafts/TubeAnalysis-MSOT.pdf');
%plot(wavelengths,msotTubeSpectra,'Color',cmap(k,:))
%% weird rank
winnerWhere = []
for m=1:size(allTubeSBR,2)
    winnerWhere(m,:) = allTubeSBR(:,m)/max(allTubeSBR(:,m))
end
figure;
errorbar([1:11]',mean(winnerWhere,1),std(winnerWhere,1));
%%
figure;hold on;
ourbar = bar([1:11]',mean(winnerWhere,1),'c')
ourerror = errorbar(mean(winnerWhere,1),std(winnerWhere,1), 'linestyle', 'none');
hold off;
ourbar(:).FaceColor = [.5 .5 .5]
ourerror.Color = [0 0 0]
XTickLabel=strrep(tubes,'_',' / ');
XTick=1:11;
set(gca, 'XTick',XTick);
set(gca, 'XTickLabel', XTickLabel); 
set(gca,'TickLabelInterpreter','none')
xtickangle(45)

plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Tubes (I.D. / O.D. in µm)'; % xlabel
plt.YLabel = 'Signal-to-background ratio [a.u.]'; %ylabel
plt.YLim = [0 1.1]
plt.BoxDim = [7 5];
plt.LineWidth = 1.5;
plt.FontName = 'Arial';
plt.FontSize = 14;
plt.TickDir = 'out';
plt.TickLength = [0.01 0.01]
plt.XMinorTick = false;
plt.ShowBox = false;
set(gca, 'Layer', 'Top');
%rectangle('Position',[7.5 0 1 1.1], 'FaceColor', [0 1 0 0.1],'LineStyle','none');
plt.export('figure-raw-drafts/TubeAnalysis-Signal-To-Background.pdf');
%%
figure;
subplot(2,2,1);imshow(allDemoImgs(:,:,6),[0 1500])
xlim([164 194]);ylim([168 199])
title(tubes(6),'Interpreter','none');
subplot(2,2,2);imshow(allDemoImgs(:,:,7),[0 1500])
xlim([145 175]);ylim([146 176])
title(tubes(7),'Interpreter','none');
subplot(2,2,3);imshow(allDemoImgs(:,:,8),[0 1500])
xlim([152 182]);ylim([151 182])
title(tubes(8),'Interpreter','none');
subplot(2,2,4);imshow(allDemoImgs(:,:,9),[0 1500])
xlim([145 175]);ylim([156 186])
title(tubes(9),'Interpreter','none');
%savefig(gcf,'figure-raw-drafts/TubeAnalysis-exampleSlices.pdf')
print('figure-raw-drafts/TubeAnalysis-exampleSlices','-dpdf')