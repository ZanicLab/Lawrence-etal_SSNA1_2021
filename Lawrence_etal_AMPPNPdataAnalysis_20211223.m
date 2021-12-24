%%% Created by GA
%%% last edited by GA on 20211223
%%% Use this to analyze Spastin binding on SSNA1-coated and not coated
%%% microtubules in the presence of AMPPNP.

close all

alpha = 0.5; 

for columnno=3:4 % 2:tubulin, 3:spastin, 4:SSNA1
    if (columnno==2)
        protein=sprintf('tubulin');
        ylimmax=9000;
    elseif (columnno==3)
        protein=sprintf('spastin');
        ylimmax=1400;
        ylimmin=0;
    elseif (columnno==4)
        protein=sprintf('SSNA1');
        ylimmax=4000;
        ylimmin=-250;
    end
    
    
tstart=3;
tend=5;
SSNAx=5;
NoSSNAx=15;
fprintf('\nAnalyzing between %d to %d mins\n',tstart,tend)
    
        
% first day
path=sprintf('~/pathtofile1/');
name=sprintf('Stabilized-20210811-Exp1-3_25nMSpastin-AMPPNP');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeansDay1=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeansDay1(any(isnan(SSNAmeansDay1), 2), :) = [];
SSNAstdsDay1=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstdsDay1(any(isnan(SSNAstdsDay1), 2), :) = [];
NoSSNAmeansDay1=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeansDay1(any(isnan(NoSSNAmeansDay1), 2), :) = [];
NoSSNAstdsDay1=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstdsDay1(any(isnan(NoSSNAstdsDay1), 2), :) = [];
NSSNA11=size(SSNAmeansDay1,1);
NNoSSNA11=size(NoSSNAmeansDay1,1);
SSNAColor=ones(NSSNA11,1);
NoSSNAColor=ones(NNoSSNA11,1);

SSNAoutlieridxDay1=isoutlier(SSNAmeansDay1(:,columnno));
NoSSNAoutlieridxDay1=isoutlier(NoSSNAmeansDay1(:,columnno));

% second day, first FOV
path=sprintf('~/pathtofile2/');
name=sprintf('Stabilized-FOV1-20210812-Exp1-3_25nMSpastin-AMPPNP');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeansDay2=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeansDay2(any(isnan(SSNAmeansDay2), 2), :) = [];
SSNAstdsDay2=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstdsDay2(any(isnan(SSNAstdsDay2), 2), :) = [];
NoSSNAmeansDay2=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeansDay2(any(isnan(NoSSNAmeansDay2), 2), :) = [];
NoSSNAstdsDay2=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstdsDay2(any(isnan(NoSSNAstdsDay2), 2), :) = [];
NSSNA12one=size(SSNAmeansDay2,1);
NNoSSNA12one=size(NoSSNAmeansDay2,1);
SSNAColor=[SSNAColor;2*ones(NSSNA12one,1)];
NoSSNAColor=[NoSSNAColor;2*ones(NNoSSNA12one,1)];
% second day, second FOV
name=sprintf('~/pathtofile3/');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeans=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeans(any(isnan(SSNAmeans), 2), :) = [];
SSNAstds=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstds(any(isnan(SSNAstds), 2), :) = [];
NoSSNAmeans=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeans(any(isnan(NoSSNAmeans), 2), :) = [];
NoSSNAstds=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstds(any(isnan(NoSSNAstds), 2), :) = [];
SSNAmeansDay2=[SSNAmeansDay2;SSNAmeans];
SSNAstdsDay2=[SSNAstdsDay2;SSNAstds];
NoSSNAmeansDay2=[NoSSNAmeansDay2;NoSSNAmeans];
NoSSNAstdsDay2=[NoSSNAstdsDay2;NoSSNAstds];
NSSNA12two=size(SSNAmeans,1);
NNoSSNA12two=size(NoSSNAmeans,1);
SSNAColor=[SSNAColor;3*ones(NSSNA12two,1)];
NoSSNAColor=[NoSSNAColor;3*ones(NNoSSNA12two,1)];

SSNAoutlieridxDay2=isoutlier(SSNAmeansDay2(:,columnno));
NoSSNAoutlieridxDay2=isoutlier(NoSSNAmeansDay2(:,columnno));

% third day, first exp, first FOV
path=sprintf('~/pathtofile4/');
name=sprintf('Stabilized-FOV1-20210817-Exp1-3_25nMSpastin-AMPPNP');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeansDay3=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeansDay3(any(isnan(SSNAmeansDay3), 2), :) = [];
SSNAstdsDay3=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstdsDay3(any(isnan(SSNAstdsDay3), 2), :) = [];
NoSSNAmeansDay3=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeansDay3(any(isnan(NoSSNAmeansDay3), 2), :) = [];
NoSSNAstdsDay3=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstdsDay3(any(isnan(NoSSNAstdsDay3), 2), :) = [];
NSSNA17one=size(SSNAmeansDay3,1);
NNoSSNA17one=size(NoSSNAmeansDay3,1);
SSNAColor=[SSNAColor;4*ones(NSSNA17one,1)];
NoSSNAColor=[NoSSNAColor;4*ones(NNoSSNA17one,1)];
% third day, first exp, second FOV
path=sprintf('~/pathtofile5/');
name=sprintf('Stabilized-FOV2-20210817-Exp1-3_25nMSpastin-AMPPNP');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeans=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeans(any(isnan(SSNAmeans), 2), :) = [];
SSNAstds=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstds(any(isnan(SSNAstds), 2), :) = [];
NoSSNAmeans=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeans(any(isnan(NoSSNAmeans), 2), :) = [];
NoSSNAstds=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstds(any(isnan(NoSSNAstds), 2), :) = [];
SSNAmeansDay3=[SSNAmeansDay3;SSNAmeans];
SSNAstdsDay3=[SSNAstdsDay3;SSNAstds];
NoSSNAmeansDay3=[NoSSNAmeansDay3;NoSSNAmeans];
NoSSNAstdsDay3=[NoSSNAstdsDay3;NoSSNAstds];
NSSNA17two=size(SSNAmeans,1);
NNoSSNA17two=size(NoSSNAmeans,1);
SSNAColor=[SSNAColor;5*ones(NSSNA17two,1)];
NoSSNAColor=[NoSSNAColor;5*ones(NNoSSNA17two,1)];
% third day, first exp, third FOV
path=sprintf('~/pathtofile6/');
name=sprintf('Stabilized-FOV3-20210817-Exp1-3_25nMSpastin-AMPPNP');
Lawrence_etal_SpastinSSNA_Intensity_20211223(path,sprintf('BGsub-%s.tif',name),name,tstart,tend);
SSNAmeans=load(sprintf('%s/WeightedMeans_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAmeans(any(isnan(SSNAmeans), 2), :) = [];
SSNAstds=load(sprintf('%s/WeightedSTDs_SSNACoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
SSNAstds(any(isnan(SSNAstds), 2), :) = [];
NoSSNAmeans=load(sprintf('%s/WeightedMeans_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAmeans(any(isnan(NoSSNAmeans), 2), :) = [];
NoSSNAstds=load(sprintf('%s/WeightedSTDs_NotCoated_%dto%dmin_BGsub-%s.dat',path,tstart,tend,name));
NoSSNAstds(any(isnan(NoSSNAstds), 2), :) = [];
SSNAmeansDay3=[SSNAmeansDay3;SSNAmeans];
SSNAstdsDay3=[SSNAstdsDay3;SSNAstds];
NoSSNAmeansDay3=[NoSSNAmeansDay3;NoSSNAmeans];
NoSSNAstdsDay3=[NoSSNAstdsDay3;NoSSNAstds];
NSSNA17three=size(SSNAmeans,1);
NNoSSNA17three=size(NoSSNAmeans,1);
SSNAColor=[SSNAColor;6*ones(NSSNA17three,1)];
NoSSNAColor=[NoSSNAColor;6*ones(NNoSSNA17three,1)];

SSNAoutlieridxDay3=isoutlier(SSNAmeansDay3(:,columnno));
NoSSNAoutlieridxDay3=isoutlier(NoSSNAmeansDay3(:,columnno));

close all

SSNAmeansAll=[SSNAmeansDay1;SSNAmeansDay2;SSNAmeansDay3];
SSNAstdsAll=[SSNAstdsDay1;SSNAstdsDay2;SSNAstdsDay3];
NoSSNAmeansAll=[NoSSNAmeansDay1;NoSSNAmeansDay2;NoSSNAmeansDay3];
NoSSNAstdsAll=[NoSSNAstdsDay1;NoSSNAstdsDay2;NoSSNAstdsDay3];

size(SSNAmeansAll)
size(NoSSNAmeansAll)

if (columnno==3)
    SSNAoutlieridxAll=( isoutlier(SSNAmeansAll(:,3),'mean') | isoutlier(SSNAmeansAll(:,4),'mean') );
    NoSSNAoutlieridxAll=( isoutlier(NoSSNAmeansAll(:,3),'mean') | isoutlier(NoSSNAmeansAll(:,4),'mean') );
end



[XSSNA,YSSNA]=BeeHive_20210107(SSNAx,SSNAmeansAll(:,columnno),50,5);
SSNAplot=sortrows([XSSNA,YSSNA],2);
SSNAsort=sortrows([SSNAmeansAll(:,columnno) SSNAColor SSNAstdsAll(:,columnno) SSNAoutlieridxAll]);
if (~isequal(SSNAsort(:,1),SSNAplot(:,2)))
    fprintf('SSNA colors are mismatching\n')
end

[XNoSSNA,YNoSSNA]=BeeHive_20210107(NoSSNAx,NoSSNAmeansAll(:,columnno),50,5);
NoSSNAplot=sortrows([XNoSSNA,YNoSSNA],2);
NoSSNAsort=sortrows([NoSSNAmeansAll(:,columnno) NoSSNAColor NoSSNAstdsAll(:,columnno) NoSSNAoutlieridxAll]);
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
    elseif (SSNAsort(i,2)==6) % Day 3, FOV3
        color=[49,163,84]/255;%[0,109,44]/255;     
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
    elseif (NoSSNAsort(i,2)==6) % Day 3, FOV3
        color=[49,163,84]/255;%[0,109,44]/255;     
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
fidmeans=fopen(sprintf('~/Desktop/TIRF/August2021/AMPPNP-%sWeightedMeans_20211223.dat',protein),'w');
SSNAWeightsDay1=1./(SSNAstdsDay1(:,columnno).*SSNAstdsDay1(:,columnno));
NSSNADay1=sum(~SSNAoutlieridxDay1);
WeightedMeanSSNADay1=sum(SSNAmeansDay1(~SSNAoutlieridxDay1,columnno).*SSNAWeightsDay1(~SSNAoutlieridxDay1,1))/sum(SSNAWeightsDay1(~SSNAoutlieridxDay1,1));
WeightedSTDSSNADay1=sqrt(var(SSNAmeansDay1(~SSNAoutlieridxDay1,columnno),SSNAWeightsDay1(~SSNAoutlieridxDay1,1)))*sqrt(NSSNADay1/(NSSNADay1-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx-1,WeightedMeanSSNADay1,WeightedSTDSSNADay1/sqrt(NSSNADay1),'o','Color','k','MarkerFaceColor',[117,107,177]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day1 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay1,WeightedSTDSSNADay1/sqrt(NSSNADay1),NSSNADay1);

SSNAWeightsDay2=1./(SSNAstdsDay2(:,columnno).*SSNAstdsDay2(:,columnno));
NSSNADay2=sum(~SSNAoutlieridxDay2);
WeightedMeanSSNADay2=sum(SSNAmeansDay2(~SSNAoutlieridxDay2,columnno).*SSNAWeightsDay2(~SSNAoutlieridxDay2,1))/sum(SSNAWeightsDay2(~SSNAoutlieridxDay2,1));
WeightedSTDSSNADay2=sqrt(var(SSNAmeansDay2(~SSNAoutlieridxDay2,columnno),SSNAWeightsDay2(~SSNAoutlieridxDay2,1)))*sqrt(NSSNADay2/(NSSNADay2-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx,WeightedMeanSSNADay2,WeightedSTDSSNADay2/sqrt(NSSNADay2),'o','Color','k','MarkerFaceColor',[222,45,38]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day2 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay2,WeightedSTDSSNADay2/sqrt(NSSNADay2),NSSNADay2);

SSNAWeightsDay3=1./(SSNAstdsDay3(:,columnno).*SSNAstdsDay3(:,columnno));
NSSNADay3=sum(~SSNAoutlieridxDay3);
WeightedMeanSSNADay3=sum(SSNAmeansDay3(~SSNAoutlieridxDay3,columnno).*SSNAWeightsDay3(~SSNAoutlieridxDay3,1))/sum(SSNAWeightsDay3(~SSNAoutlieridxDay3,1));
WeightedSTDSSNADay3=sqrt(var(SSNAmeansDay3(~SSNAoutlieridxDay3,columnno),SSNAWeightsDay3(~SSNAoutlieridxDay3,1)))*sqrt(NSSNADay3/(NSSNADay3-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx+1,WeightedMeanSSNADay3,WeightedSTDSSNADay3/sqrt(NSSNADay3),'o','Color','k','MarkerFaceColor',[49,163,84]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day3 SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNADay3,WeightedSTDSSNADay3/sqrt(NSSNADay3),NSSNADay3);


NoSSNAWeightsDay1=1./(NoSSNAstdsDay1(:,columnno).*NoSSNAstdsDay1(:,columnno));
NNoSSNADay1=sum(~NoSSNAoutlieridxDay1);
WeightedMeanNoSSNADay1=sum(NoSSNAmeansDay1(~NoSSNAoutlieridxDay1,columnno).*NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1))/sum(NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1));
WeightedSTDNoSSNADay1=sqrt(var(NoSSNAmeansDay1(~NoSSNAoutlieridxDay1,columnno),NoSSNAWeightsDay1(~NoSSNAoutlieridxDay1,1)))*sqrt(NNoSSNADay1/(NNoSSNADay1-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx-1,WeightedMeanNoSSNADay1,WeightedSTDNoSSNADay1,'o','Color','k','MarkerFaceColor',[117,107,177]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day1 NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNADay1,WeightedSTDNoSSNADay1/sqrt(NNoSSNADay1),NNoSSNADay1);

NoSSNAWeightsDay2=1./(NoSSNAstdsDay2(:,columnno).*NoSSNAstdsDay2(:,columnno));
NNoSSNADay2=sum(~NoSSNAoutlieridxDay2);
WeightedMeanNoSSNADay2=sum(NoSSNAmeansDay2(~NoSSNAoutlieridxDay2,columnno).*NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1))/sum(NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1));
WeightedSTDNoSSNADay2=sqrt(var(NoSSNAmeansDay2(~NoSSNAoutlieridxDay2,columnno),NoSSNAWeightsDay2(~NoSSNAoutlieridxDay2,1)))*sqrt(NNoSSNADay2/(NNoSSNADay2-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx,WeightedMeanNoSSNADay2,WeightedSTDNoSSNADay2/sqrt(NNoSSNADay2),'o','Color','k','MarkerFaceColor',[222,45,38]/255,'MarkerSize',8,'Linewidth',1)
hold on
fprintf(fidmeans,'Day2 NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNADay2,WeightedSTDNoSSNADay2/sqrt(NNoSSNADay2),NNoSSNADay2);

NoSSNAWeightsDay3=1./(NoSSNAstdsDay3(:,columnno).*NoSSNAstdsDay3(:,columnno));
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
SSNAWeightsAll=1./(SSNAstdsAll(:,columnno).*SSNAstdsAll(:,columnno));
NSSNAAll=sum(~SSNAoutlieridxAll);
WeightedMeanSSNAAll=sum(SSNAmeansAll(~SSNAoutlieridxAll,columnno).*SSNAWeightsAll(~SSNAoutlieridxAll,1))/sum(SSNAWeightsAll(~SSNAoutlieridxAll,1));
WeightedSTDSSNAAll=sqrt(var(SSNAmeansAll(~SSNAoutlieridxAll,columnno),SSNAWeightsAll(~SSNAoutlieridxAll,1)))*sqrt(NSSNAAll/(NSSNAAll-1)); % to correct for normalizing with N-1, instead of N
errorbar(SSNAx,WeightedMeanSSNAAll,WeightedSTDSSNAAll/sqrt(NSSNAAll),'x','Color','k','MarkerFaceColor','k','MarkerSize',12,'Linewidth',1)
hold on
fprintf(fidmeans,'AllDays SSNA:\t%f\t%f\t%d\n',WeightedMeanSSNAAll,WeightedSTDSSNAAll/sqrt(NSSNAAll),NSSNAAll);

NoSSNAWeightsAll=1./(NoSSNAstdsAll(:,columnno).*NoSSNAstdsAll(:,columnno));
NNoSSNAAll=sum(~NoSSNAoutlieridxAll);
WeightedMeanNoSSNAAll=sum(NoSSNAmeansAll(~NoSSNAoutlieridxAll,columnno).*NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1))/sum(NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1));
WeightedSTDNoSSNAAll=sqrt(var(NoSSNAmeansAll(~NoSSNAoutlieridxAll,columnno),NoSSNAWeightsAll(~NoSSNAoutlieridxAll,1)))*sqrt(NNoSSNAAll/(NNoSSNAAll-1)); % to correct for normalizing with N-1, instead of N
errorbar(NoSSNAx,WeightedMeanNoSSNAAll,WeightedSTDNoSSNAAll/sqrt(NNoSSNAAll),'x','Color','k','MarkerFaceColor','k','MarkerSize',12,'Linewidth',1)
hold on
fprintf(fidmeans,'AllDays NoSSNA:\t%f\t%f\t%d\n',WeightedMeanNoSSNAAll,WeightedSTDNoSSNAAll/sqrt(NNoSSNAAll),NNoSSNAAll);



fig=figure(columnno);
ylim([ylimmin inf])
xlim([0 20])
xticks([5 15])
xticklabels({'SSNA coated','Not coated'})
ylabel(sprintf('%s intensity (a.u.)',protein))
pbaspect([1 1 1])
exportgraphics(fig,sprintf('~/outputpath/%sIntensity_AMPPNP_20211223-PointsWithErrorbars.png',protein));
saveas(fig,sprintf('~/outputpath/%sIntensity_AMPPNP_20211223-PointsWithErrorbars.fig',protein));

end

