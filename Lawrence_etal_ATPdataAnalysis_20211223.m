%%% Created by GA
%%% last edited by GA on 20211223
%%% Use this to analyze Spastin severing on SSNA1-coated and not coated
%%% microtubules in the presence of ATP.

close all

alpha = 0.5; 


columnno=1;
protein=sprintf('tubulin');


SSNAx=5;
NoSSNAx=15;
    
        
% first day
path=sprintf('~/pathtofile1/');
name=sprintf('Stabilized-20210811-Exp1-4_25nMSpastin-ATP');
Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,0,10)
SSNADataDay1=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name));
SSNAmeansDay1=SSNADataDay1(:,3);
SSNAmeansDay1(any(isnan(SSNAmeansDay1), 2), :) = [];
SSNAsemsDay1=SSNADataDay1(:,4);
SSNAsemsDay1(any(isnan(SSNAsemsDay1), 2), :) = [];
SSNANpixelDay1=SSNADataDay1(:,5);
SSNANpixelDay1(any(isnan(SSNANpixelDay1), 2), :) = [];
NoSSNADataDay1=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name));
NoSSNAmeansDay1=NoSSNADataDay1(:,3);
NoSSNAmeansDay1(any(isnan(NoSSNAmeansDay1), 2), :) = [];
NoSSNAsemsDay1=NoSSNADataDay1(:,4);
NoSSNAsemsDay1(any(isnan(NoSSNAsemsDay1), 2), :) = [];
NoSSNANpixelDay1=NoSSNADataDay1(:,5);
NoSSNANpixelDay1(any(isnan(NoSSNANpixelDay1), 2), :) = [];
NSSNA11=size(SSNAmeansDay1,1);
NNoSSNA11=size(NoSSNAmeansDay1,1);
SSNAColor=ones(NSSNA11,1);
NoSSNAColor=ones(NNoSSNA11,1);

SSNAoutlieridxDay1=isoutlier(SSNAmeansDay1(:,columnno));
NoSSNAoutlieridxDay1=isoutlier(NoSSNAmeansDay1(:,columnno));

% second day, first FOV
path=sprintf('~/pathtofile2/');
name=sprintf('Stabilized-FOV1-20210812-Exp1-4_25nMSpastin-ATP');
Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,0,10)
SSNADataDay2one=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name));
SSNAmeansDay2one=SSNADataDay2one(:,3);
SSNAmeansDay2one(any(isnan(SSNAmeansDay2one), 2), :) = [];
SSNAsemsDay2one=SSNADataDay2one(:,4);
SSNAsemsDay2one(any(isnan(SSNAsemsDay2one), 2), :) = [];
NoSSNADataDay2one=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name));
NoSSNAmeansDay2one=NoSSNADataDay2one(:,3);
NoSSNAmeansDay2one(any(isnan(NoSSNAmeansDay2one), 2), :) = [];
NoSSNAsemsDay2one=NoSSNADataDay2one(:,4);
NoSSNAsemsDay2one(any(isnan(NoSSNAsemsDay2one), 2), :) = [];
NSSNA12one=size(SSNAmeansDay2one,1);
NNoSSNA12one=size(NoSSNAmeansDay2one,1);
SSNAColor12one=2*ones(NSSNA12one,1);
NoSSNAColor12one=2*ones(NNoSSNA12one,1);
% second day, second FOV
name=sprintf('~/pathtofile3/');
Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,0,10)
SSNADataDay2two=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name));
SSNAmeansDay2two=SSNADataDay2two(:,3);
SSNAmeansDay2two(any(isnan(SSNAmeansDay2two), 2), :) = [];
SSNAsemsDay2two=SSNADataDay2two(:,4);
SSNAsemsDay2two(any(isnan(SSNAsemsDay2two), 2), :) = [];
NoSSNADataDay2two=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name));
NoSSNAmeansDay2two=NoSSNADataDay2two(:,3);
NoSSNAmeansDay2two(any(isnan(NoSSNAmeansDay2two), 2), :) = [];
NoSSNAsemsDay2two=NoSSNADataDay2two(:,4);
NoSSNAsemsDay2two(any(isnan(NoSSNAsemsDay2two), 2), :) = [];
NSSNA12two=size(SSNAmeansDay2two,1);
NNoSSNA12two=size(NoSSNAmeansDay2two,1);
SSNAColor12two=3*ones(NSSNA12two,1);
NoSSNAColor12two=3*ones(NNoSSNA12two,1);

SSNAmeansDay2=[SSNAmeansDay2one;SSNAmeansDay2two];
SSNAsemsDay2=[SSNAsemsDay2one;SSNAsemsDay2two];
NoSSNAmeansDay2=[NoSSNAmeansDay2one;NoSSNAmeansDay2two];
NoSSNAsemsDay2=[NoSSNAsemsDay2one;NoSSNAsemsDay2two];
SSNAColor=[SSNAColor;SSNAColor12one;SSNAColor12two];
NoSSNAColor=[NoSSNAColor;NoSSNAColor12one;NoSSNAColor12two];
SSNAoutlieridxDay2=isoutlier(SSNAmeansDay2(:,columnno));
NoSSNAoutlieridxDay2=isoutlier(NoSSNAmeansDay2(:,columnno));

% third day, first exp, first FOV
path=sprintf('~/pathtofile4/');
name=sprintf('Stabilized-FOV1-20210817-Exp1-4_25nMSpastin-ATP');
Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,0,15)
SSNADataDay3one=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name));
SSNAmeansDay3one=SSNADataDay3one(:,3);
SSNAmeansDay3one(any(isnan(SSNAmeansDay3one), 2), :) = [];
SSNAsemsDay3one=SSNADataDay3one(:,4);
SSNAsemsDay3one(any(isnan(SSNAsemsDay3one), 2), :) = [];
NoSSNADataDay3one=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name));
NoSSNAmeansDay3one=NoSSNADataDay3one(:,3);
NoSSNAmeansDay3one(any(isnan(NoSSNAmeansDay3one), 2), :) = [];
NoSSNAsemsDay3one=NoSSNADataDay3one(:,4);
NoSSNAsemsDay3one(any(isnan(NoSSNAsemsDay3one), 2), :) = [];
NSSNA17one=size(SSNAmeansDay3one,1);
NNoSSNA17one=size(NoSSNAmeansDay3one,1);
SSNAColor17one=4*ones(NSSNA17one,1);
NoSSNAColor17one=4*ones(NNoSSNA17one,1);
% third day, first exp, second FOV
name=sprintf('~/pathtofile5/');
Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,0,15)
SSNADataDay3two=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name));
SSNAmeansDay3two=SSNADataDay3two(:,3);
SSNAmeansDay3two(any(isnan(SSNAmeansDay3two), 2), :) = [];
SSNAsemsDay3two=SSNADataDay3two(:,4);
SSNAsemsDay3two(any(isnan(SSNAsemsDay3two), 2), :) = [];
NoSSNADataDay3two=load(sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name));
NoSSNAmeansDay3two=NoSSNADataDay3two(:,3);
NoSSNAmeansDay3two(any(isnan(NoSSNAmeansDay3two), 2), :) = [];
NoSSNAsemsDay3two=NoSSNADataDay3two(:,4);
NoSSNAsemsDay3two(any(isnan(NoSSNAsemsDay3two), 2), :) = [];
NSSNA17two=size(SSNAmeansDay3two,1);
NNoSSNA17two=size(NoSSNAmeansDay3two,1);
SSNAColor17two=5*ones(NSSNA17two,1);
NoSSNAColor17two=5*ones(NNoSSNA17two,1);


SSNAmeansDay3=[SSNAmeansDay3one;SSNAmeansDay3two];
SSNAsemsDay3=[SSNAsemsDay3one;SSNAsemsDay3two];
NoSSNAmeansDay3=[NoSSNAmeansDay3one;NoSSNAmeansDay3two];
NoSSNAsemsDay3=[NoSSNAsemsDay3one;NoSSNAsemsDay3two];
SSNAColor=[SSNAColor;SSNAColor17one;SSNAColor17two];
NoSSNAColor=[NoSSNAColor;NoSSNAColor17one;NoSSNAColor17two];
SSNAoutlieridxDay3=isoutlier(SSNAmeansDay3(:,columnno));
NoSSNAoutlieridxDay3=isoutlier(NoSSNAmeansDay3(:,columnno));



close all


SSNAmeansAll=[SSNAmeansDay1;SSNAmeansDay2;SSNAmeansDay3];
SSNAsemsAll=[SSNAsemsDay1;SSNAsemsDay2;SSNAsemsDay3];
NoSSNAmeansAll=[NoSSNAmeansDay1;NoSSNAmeansDay2;NoSSNAmeansDay3];
NoSSNAsemsAll=[NoSSNAsemsDay1;NoSSNAsemsDay2;NoSSNAsemsDay3];

size(SSNAmeansAll)
size(NoSSNAmeansAll)
SSNAoutlieridxAll=isoutlier(SSNAmeansAll,'mean');
NoSSNAoutlieridxAll=isoutlier(NoSSNAmeansAll,'mean');

[XSSNA,YSSNA]=BeeHive_20210107(SSNAx,SSNAmeansAll(:,columnno),20,5);
SSNAplot=sortrows([XSSNA,YSSNA],2);
SSNAsort=sortrows([SSNAmeansAll(:,columnno) SSNAColor SSNAsemsAll(:,columnno) SSNAoutlieridxAll]);
if (~isequal(SSNAsort(:,1),SSNAplot(:,2)))
    fprintf('SSNA colors are mismatching\n')
end

[XNoSSNA,YNoSSNA]=BeeHive_20210107(NoSSNAx,NoSSNAmeansAll(:,columnno),20,5);
NoSSNAplot=sortrows([XNoSSNA,YNoSSNA],2);
NoSSNAsort=sortrows([NoSSNAmeansAll(:,columnno) NoSSNAColor NoSSNAsemsAll(:,columnno) NoSSNAoutlieridxAll]);
if (~isequal(NoSSNAsort(:,1),NoSSNAplot(:,2)))
    fprintf('NoSSNA colors are mismatching\n')
end



    SSNAoutlieridx=SSNAsort(:,4);
    NoSSNAoutlieridx=NoSSNAsort(:,4);


for i=1:size(SSNAplot,1)
    if (SSNAsort(i,2)==1) % Day 1
        color=[117,107,177]/255;
    elseif (SSNAsort(i,2)==2) % Day 2, FOV1
        color=[222,45,38]/255;%[222,45,38]/255;
    elseif (SSNAsort(i,2)==3) % Day 2, FOV2
        color=[222,45,38]/255;%[165,15,21]/255;
    elseif (SSNAsort(i,2)==4) % Day 3, FOV1
        color=[49,163,84]/255;%[116,196,118]/255;   
    elseif (SSNAsort(i,2)==5) % Day 3, FOV2
        color=[49,163,84]/255;%[49,163,84]/255;   
    else
        fprintf('something is wrong\n')
    end
    
    if (~SSNAoutlieridx(i,1))
        figure(columnno)
        h=errorbar(SSNAplot(i,1),SSNAplot(i,2),SSNAsort(i,3),'Color',color);
        hold on
        set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
        set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
        hold on
        scatter(SSNAplot(i,1),SSNAplot(i,2),'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
    end
end

for i=1:size(NoSSNAplot,1)
    if (NoSSNAsort(i,2)==1) % Day 1
        color=[117,107,177]/255;
    elseif (NoSSNAsort(i,2)==2) % Day 2, FOV1
        color=[222,45,38]/255;%[222,45,38]/255;
    elseif (NoSSNAsort(i,2)==3) % Day 2, FOV2
        color=[222,45,38]/255;%[165,15,21]/255;
    elseif (NoSSNAsort(i,2)==4) % Day 3, FOV1
        color=[49,163,84]/255;%[116,196,118]/255;   
    elseif (NoSSNAsort(i,2)==5) % Day 3, FOV2
        color=[49,163,84]/255;%[49,163,84]/255;   
    else
        fprintf('something is wrong\n')
    end
            
    if (~NoSSNAoutlieridx(i,1))
        figure(columnno)
        h=errorbar(NoSSNAplot(i,1),NoSSNAplot(i,2),SSNAsort(i,3),'Color',color);
        hold on
        set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
        set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
        hold on
        scatter(NoSSNAplot(i,1),NoSSNAplot(i,2),'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%%% Each Day Means %%%
%%%%%%%%%%%%%%%%%%%%%%
figure(columnno)
fidmeans=fopen(sprintf('~/Desktop/TIRF/August2021/ATP-%sWeightedMeans_20211223.dat',protein),'w');
SSNAWeightsDay1=1./(SSNAsemsDay1(:,columnno).*SSNAsemsDay1(:,columnno));
NSSNADay1=sum(~SSNAoutlieridxDay1);
WeightedMeanSSNADay1=sum(SSNAmeansDay1(~SSNAoutlieridxDay1,columnno).*SSNAWeightsDay1(~SSNAoutlieridxDay1,1))/sum(SSNAWeightsDay1(~SSNAoutlieridxDay1,1));
WeightedSTDSSNADay1=sqrt(var(SSNAmeansDay1(~SSNAoutlieridxDay1,columnno),SSNAWeightsDay1(~SSNAoutlieridxDay1,1)))*sqrt(NSSNADay1/(NSSNADay1-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx-1,WeightedMeanSSNADay1,WeightedSTDSSNADay1/sqrt(NSSNADay1),'o','Color','k','MarkerFaceColor',[117,107,177]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day1 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay1,WeightedSTDSSNADay1/sqrt(NSSNADay1),NSSNADay1);

SSNAWeightsDay2=1./(SSNAsemsDay2(:,columnno).*SSNAsemsDay2(:,columnno));
NSSNADay2=sum(~SSNAoutlieridxDay2);
WeightedMeanSSNADay2=sum(SSNAmeansDay2(~SSNAoutlieridxDay2,columnno).*SSNAWeightsDay2(~SSNAoutlieridxDay2,1))/sum(SSNAWeightsDay2(~SSNAoutlieridxDay2,1));
WeightedSTDSSNADay2=sqrt(var(SSNAmeansDay2(~SSNAoutlieridxDay2,columnno),SSNAWeightsDay2(~SSNAoutlieridxDay2,1)))*sqrt(NSSNADay2/(NSSNADay2-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx,WeightedMeanSSNADay2,WeightedSTDSSNADay2/sqrt(NSSNADay2),'o','Color','k','MarkerFaceColor',[222,45,38]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day2 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay2,WeightedSTDSSNADay2/sqrt(NSSNADay2),NSSNADay2);

SSNAWeightsDay3=1./(SSNAsemsDay3(:,columnno).*SSNAsemsDay3(:,columnno));
NSSNADay3=sum(~SSNAoutlieridxDay3);
WeightedMeanSSNADay3=sum(SSNAmeansDay3(~SSNAoutlieridxDay3,columnno).*SSNAWeightsDay3(~SSNAoutlieridxDay3,1))/sum(SSNAWeightsDay3(~SSNAoutlieridxDay3,1));
WeightedSTDSSNADay3=sqrt(var(SSNAmeansDay3(~SSNAoutlieridxDay3,columnno),SSNAWeightsDay3(~SSNAoutlieridxDay3,1)))*sqrt(NSSNADay3/(NSSNADay3-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx+1,WeightedMeanSSNADay3,WeightedSTDSSNADay3/sqrt(NSSNADay3),'o','Color','k','MarkerFaceColor',[49,163,84]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day3 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay3,WeightedSTDSSNADay3/sqrt(NSSNADay3),NSSNADay3);


NoSSNAWeightsDay1=1./(NoSSNAsemsDay1(:,columnno).*NoSSNAsemsDay1(:,columnno));
NNoSSNADay1=sum(~NoSSNAoutlieridxDay1);
WeightedMeanNoSSNADay1=sum(NoSSNAmeansDay1(~NoSSNAoutlieridxDay1,columnno).*NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1))/sum(NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1));
WeightedSTDNoSSNADay1=sqrt(var(NoSSNAmeansDay1(~NoSSNAoutlieridxDay1,columnno),NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1)))*sqrt(NNoSSNADay1/(NNoSSNADay1-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx-1,WeightedMeanNoSSNADay1,WeightedSTDNoSSNADay1,'o','Color','k','MarkerFaceColor',[117,107,177]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day1 NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNADay1,WeightedSTDNoSSNADay1/sqrt(NNoSSNADay1),NNoSSNADay1);

NoSSNAWeightsDay2=1./(NoSSNAsemsDay2(:,columnno).*NoSSNAsemsDay2(:,columnno));
NNoSSNADay2=sum(~NoSSNAoutlieridxDay2);
WeightedMeanNoSSNADay2=sum(NoSSNAmeansDay2(~NoSSNAoutlieridxDay2,columnno).*NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1))/sum(NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1));
WeightedSTDNoSSNADay2=sqrt(var(NoSSNAmeansDay2(~NoSSNAoutlieridxDay2,columnno),NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1)))*sqrt(NNoSSNADay2/(NNoSSNADay2-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx,WeightedMeanNoSSNADay2,WeightedSTDNoSSNADay2/sqrt(NNoSSNADay2),'o','Color','k','MarkerFaceColor',[222,45,38]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day2 NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNADay2,WeightedSTDNoSSNADay2/sqrt(NNoSSNADay2),NNoSSNADay2);

NoSSNAWeightsDay3=1./(NoSSNAsemsDay3(:,columnno).*NoSSNAsemsDay3(:,columnno));
NNoSSNADay3=sum(~NoSSNAoutlieridxDay3);
WeightedMeanNoSSNADay3=sum(NoSSNAmeansDay3(~NoSSNAoutlieridxDay3,columnno).*NoSSNAWeightsDay3(~NoSSNAoutlieridxDay3,1))/sum(NoSSNAWeightsDay3(~NoSSNAoutlieridxDay3,1));
WeightedSTDNoSSNADay3=sqrt(var(NoSSNAmeansDay3(~NoSSNAoutlieridxDay3,columnno),NoSSNAWeightsDay3(~NoSSNAoutlieridxDay3,1)))*sqrt(NNoSSNADay3/(NNoSSNADay3-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx+1,WeightedMeanNoSSNADay3,WeightedSTDNoSSNADay3/sqrt(NNoSSNADay3),'o','Color','k','MarkerFaceColor',[49,163,84]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day3 NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNADay3,WeightedSTDNoSSNADay3/sqrt(NNoSSNADay3),NNoSSNADay3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% means of "each day means" %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fidmeans,'\n');
SSNAWeightsAll=1./(SSNAsemsAll(:,columnno).*SSNAsemsAll(:,columnno));
NSSNAAll=sum(~SSNAoutlieridxAll);
WeightedMeanSSNAAll=sum(SSNAmeansAll(~SSNAoutlieridxAll,columnno).*SSNAWeightsAll(~SSNAoutlieridxAll,1))/sum(SSNAWeightsAll(~SSNAoutlieridxAll,1));
WeightedSTDSSNAAll=sqrt(var(SSNAmeansAll(~SSNAoutlieridxAll,columnno),SSNAWeightsAll(~SSNAoutlieridxAll,1)))*sqrt(NSSNAAll/(NSSNAAll-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx,WeightedMeanSSNAAll,WeightedSTDSSNAAll/sqrt(NSSNAAll),'x','Color','k','MarkerFaceColor','k','MarkerSize',12,'Linewidth',1)
hold on
fprintf(fidmeans,'AllDays SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNAAll,WeightedSTDSSNAAll/sqrt(NSSNAAll),NSSNAAll);

NoSSNAWeightsAll=1./(NoSSNAsemsAll(:,columnno).*NoSSNAsemsAll(:,columnno));
NNoSSNAAll=sum(~NoSSNAoutlieridxAll);
WeightedMeanNoSSNAAll=sum(NoSSNAmeansAll(~NoSSNAoutlieridxAll,columnno).*NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1))/sum(NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1));
WeightedSTDNoSSNAAll=sqrt(var(NoSSNAmeansAll(~NoSSNAoutlieridxAll,columnno),NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1)))*sqrt(NNoSSNAAll/(NNoSSNAAll-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx,WeightedMeanNoSSNAAll,WeightedSTDNoSSNAAll/sqrt(NNoSSNAAll),'x','Color','k','MarkerFaceColor','k','MarkerSize',12,'Linewidth',1)
hold on
fprintf(fidmeans,'AllDays NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNAAll,WeightedSTDNoSSNAAll/sqrt(NNoSSNAAll),NNoSSNAAll);



fig=figure(columnno);
ylim([-0.2 1.6])
xlim([0 20])
xticks([5 15])
xticklabels({'SSNA coated','Not coated'})
ylabel(sprintf('%s intensity (a.u.)',protein))
pbaspect([1 1 1])
exportgraphics(fig,sprintf('~/outputpath/%sIntensity_ATP_20211223-PointsWithErrorbars.png',protein));
saveas(fig,sprintf('~/outputpath/%sIntensity_ATP_20211223-PointsWithErrorbars.fig',protein));



