%%% Created by GA
%%% last edited by GA on 20211222
%%% Use this to postprocess curvature data to plot intensity as a function
%%% of curvature.
%%%
%%% USAGE:
%%% for taxol
%%% CurvatureIntensity_PostProcess_20211222GA('./taxol/','Curves-SUM_Stabilized_5uM488-SSNA1_561TaxolSeeds_30s_001_pos','SSNAchannel-SUM_Stabilized_5uM488-SSNA1_561TaxolSeeds_30s_001_pos')
%%% for dynamic ones
%%% CurvatureIntensity_PostProcess_20211222GA('./.','20180316-curves','MAX_SSNA1Channel_20180316_Stabilized 5uM488SSNA1_8uMTub647_5s_001-1')

function Lawrence_etal_CurvatureIntensity_PostProcess_20211222(path,excelname,maxprojname)
close all
clearvars -except path excelname maxprojname
totalNUniqueCurve=0;
totalnoofpoints=0;
curvaturebin=0.1;

writetable(array2table([[],[]]),sprintf('%s/Curvature-Intensity_AllPoints_%s.dat',path,datetime('today')));
fidaverages=fopen(sprintf('%s/Curvature-Intensity_Averages_%s.dat',path,datetime('today')),'w');
fprintf(fidaverages,'AVEcurvature\tAVEintensity\tSTDcurvature\tSTDintensity\n');

%for pos=7:10 % use this line for taxol
% read the max projection image
%[MaxProj,map]=imread(sprintf('%s/%s%d.tif',path,maxprojname,pos)); % use this line for taxol
[MaxProj,map]=imread(sprintf('%s/%s.tif',path,maxprojname)); % use this line for dynamic
Maskall=zeros(size(MaxProj));
% read excel data file
%T=readtable(sprintf('%s/%s%d.csv',path,excelname,pos),'ReadVariableNames',false); % use this line for taxol
T=readtable(sprintf('%s/%s.csv',path,excelname),'ReadVariableNames',false); % use this line for dynamic
UniqueCurve = table2cell(unique(T(:,1)));
NUniqueCurve=size(UniqueCurve,1);
totalNUniqueCurve=totalNUniqueCurve+NUniqueCurve;

Ncolors=ceil(NUniqueCurve/2);
colormatrixone=autumn(Ncolors);
colormatrixtwo=winter(Ncolors);

colormatrix=NaN(2*Ncolors,3);
for i=1:Ncolors
    colormatrix(i*2-1,:)=colormatrixone(i,:);
    colormatrix(i*2,:)=colormatrixtwo(end-i+1,:);
end
rotateangle=360/NUniqueCurve;


for c=1:NUniqueCurve
    idx=find(strcmp(table2cell(T(:,1)),UniqueCurve(c,1)));
    
    Xcoord=table2array(T(idx,5))./0.16;
    Ycoord=table2array(T(idx,6))./0.16;
    
    Maskindividual=zeros(size(MaxProj));
    Curvature=table2array(T(idx,7));
    Intensity=NaN(size(Xcoord,1),1);
    length=0;
    for s=1:size(Xcoord,1)
        if (s>1)
            length=length+0.16*sqrt( (Xcoord(s,1)-Xcoord(s-1,1))^2 + (Ycoord(s,1)-Ycoord(s-1,1))^2 );
        end
        Maskall(floor(Ycoord(s,1)),floor(Xcoord(s,1)))=1;
        Maskindividual(floor(Ycoord(s,1)),floor(Xcoord(s,1)))=1;
        Intensity(s,1)=MaxProj(floor(Ycoord(s,1)),floor(Xcoord(s,1)));
    end
    totalnoofpoints=totalnoofpoints+size(Curvature,1);
    writetable(array2table([Curvature,Intensity]),sprintf('%s/Curvature-Intensity_AllPoints_%s.dat',path,datetime('today')),'WriteMode','Append');
    
    
    % plot selected curves
    alpha=(atan2(Ycoord(1,1)-Ycoord(end,1),Xcoord(end,1)-Xcoord(1,1)))*180/pi;
    theta=((c-1)*rotateangle)-alpha;
    Rotation=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    RotatedCoords=[Xcoord, Ycoord]*Rotation;
    %fprintf('%d\t%.1f\t%.1f\t%.1f\n',c,((c-1)*rotateangle),alpha,theta)
    figure(2)
    plot(0.16*(RotatedCoords(:,1)-RotatedCoords(1,1)),0.16*(RotatedCoords(1,2)-RotatedCoords(:,2)),'LineWidth',2,'color',colormatrix(c,:))
    hold on
   
    [maxcurvature,idxmaxcurvature]=max(Curvature(:,1)); 
    plot(0.16*(RotatedCoords(idxmaxcurvature,1)-RotatedCoords(1,1)),0.16*(RotatedCoords(1,2)-RotatedCoords(idxmaxcurvature,2)),'k*')
    hold on
end

%end % use this line for taxol

fig=figure(2);
pbaspect([1 1 1])

xlabel('x-coordinate (µm)')
ylabel('y-coordinate (µm)')
exportgraphics(fig,sprintf('%s/Curls_%s.pdf',path,datetime('today')));
saveas(fig,sprintf('%s/Curls_%s.fig',path,datetime('today')));


CurvInt=readtable(sprintf('%s/Curvature-Intensity_AllPoints_%s.dat',path,datetime('today')));
if (totalnoofpoints==size(CurvInt,1))
    AllCurvature=table2array(CurvInt(:,1));
    AllIntensity=table2array(CurvInt(:,2));
    
    maxcurvature=max(AllCurvature);
    lastcurvaturebin=maxcurvature-mod(maxcurvature,curvaturebin)+curvaturebin;
    Nbins=floor(lastcurvaturebin/curvaturebin);
    bincounter=0;
    Averages=NaN(Nbins,2);
    for leftbin=0:curvaturebin:lastcurvaturebin-curvaturebin
        bincounter=bincounter+1;
        rightbin=leftbin+curvaturebin;
        idx=(AllCurvature>=leftbin) & (AllCurvature<rightbin);
        BinCurvature=AllCurvature(idx,1);
        BinIntensity=AllIntensity(idx,1);
        outlieridx=isoutlier(BinIntensity);
        BinCurvature(outlieridx)=[];
        BinIntensity(outlieridx)=[];
        NpointsInBin=size(BinCurvature,1);
        if (NpointsInBin>10)
            AVEcurvature=mean(BinCurvature);
            AVEintensity=mean(BinIntensity);
            STDcurvature=std(BinCurvature);
            STDintensity=std(BinIntensity);
            figure(1)
            plot(BinCurvature,BinIntensity,'.','color',[0.8 0.8 0.8])
            hold on
            errorbar(AVEcurvature,AVEintensity,STDintensity,STDintensity,STDcurvature,STDcurvature,'ko','MarkerFaceColor','k','LineWidth',2)
            hold on
            Averages(bincounter,1)=AVEcurvature;
            Averages(bincounter,2)=AVEintensity;
            fprintf(fidaverages,'%f\t%f\t%f\t%f\n',AVEcurvature,AVEintensity,STDcurvature,STDintensity);
            fprintf('%d\t%d\t%d\n',sum(idx),size(BinCurvature,1),size(BinIntensity,1)); 
        end
    end
    
    % remove all NaN rows that were skipped due to number of points restriction
    Averages(all(isnan(Averages),2),:) = [];
    
    
    ft='A + (x*B)';
    [fo, gof] = fit(Averages(:,1),Averages(:,2), ft, 'StartPoint',[1 1]);
    CI=confint(fo);
    fprintf('intercept= %f [%f - %f]\n',fo.A,CI(1,1),CI(2,1))
    fprintf('slope= %f [%f - %f]\n',fo.B,CI(1,2),CI(2,2))
        
    fig=figure(1);
    p1=plot(fo);
    set(p1,'lineWidth',2,'color','k');
    hold on
    ylim([0 inf]);
    pbaspect([1 1 1]);
    xlabel('curvature (µm^-^1)')
    ylabel('intensity (a.u.)')
    legend('hide')
    set(gca, 'box', 'off')
    exportgraphics(fig,sprintf('%s/Curvature-Intensity_SD_%s.pdf',path,datetime('today')));
    saveas(fig,sprintf('%s/Curvature-Intensity_SD_%s.fig',path,datetime('today')));
    
    
    
else
    fprintf('Total number of points is wrong for Curvature-vs-Intensity plot\n');
end

fclose(fidaverages);

fprintf('Total number of curves used: %d\n',totalNUniqueCurve)



end

