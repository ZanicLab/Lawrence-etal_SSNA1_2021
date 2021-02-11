% Created by GA on 19  Nov 2020
% Last edited by GA on 10 Feb 2021 to add input descriptions and comments
% Read kymographdirect tracking results on tubulin channel from excel file.
% Perform a sliding window velocity analysis with a given window size
% Optional: Perform intensity analysis on SSNA1 channel
% Inputs:
% 1) path to excel file,
% 2) excel file name,
% 3) intensity analysis true/false,
% 4) window size for smoothening tracking results. 0 was used for this analysis,
% 5) window size for velocity calculation
% 6) interval for shifting velocity window
% (Window sizes are in number of frames.)

% USAGE: Lawrence_etal_SSNA1_SlidingWindowVelocity('~/ExamplePath','ExampleExcelFile.xlsx',false,0,2,2)

function Lawrence_etal_SSNA1_SlidingWindowVelocity(path,excelfilename,intensityanalysis,positionsmoothenwindowsizeinminutes,velocitywindowsizeinminutes,shiftintervalinminutes)
close all
clearvars -except path excelfilename intensityanalysis positionsmoothenwindowsizeinminutes velocitywindowsizeinminutes shiftintervalinminutes

plotflag=true;

% Beehive plot parameters
Nbins=10;
xwidth=1.5;

% Time conversion factor to seconds (assuming excel file includes ms).
timeconversionfactor=1/1000; % 1/1000

% Position conversion factor to nm (assuming excel file includes um).
positionconversionfactor=1000;

% read tracking data
DATA=readtable(sprintf('%s/%s',path,excelfilename),'NumHeaderLines',3);

% read kymograph names (if intensity analysis is also done)
if (intensityanalysis)
    KYMOindexes = table2array(readtable(sprintf('%s/%s',path,excelfilename), 'Range', '1:1'));
    sumwidth=5; 
end

% linear fit formula
ft = fittype('(a*x)+b');

velocitywindowsizeinseconds=velocitywindowsizeinminutes*60;

maxtimeinseconds=max(max(table2array(DATA(:,2:3:size(DATA,2)))))*timeconversionfactor;

% find midpoint intervals for velocity fit
midpoints=transpose([velocitywindowsizeinseconds/2:(shiftintervalinminutes*60):maxtimeinseconds]);


TrackVelocity=NaN(floor(size(DATA,2)/3),size(midpoints,1));
TrackVelocityError=NaN(floor(size(DATA,2)/3),size(midpoints,1));
if (intensityanalysis)
    TrackIntensity=NaN(floor(size(DATA,2)/3),size(midpoints,1));
end

%colormap=jet(7);
% 1)dark green, 2)orange, 3)purple, 4)pink, 5)light green, 6)yellow, 7)light blue
%colormapmatrix=[ 27,158,119 ; 217,95,2 ; 117,112,179 ; 231,41,138 ; 102,166,30 ; 230,171,2 ; 166,206,227 ]/255;
% 1)black, 2)dark purple, 3)purple, 4)blue, 5)green, 6)yellow, 7)orange, 8)red
colormapmatrix=[ 0,0,0; 84,39,143 ; 117,112,179 ; 0,55,200 ; 27,158,119 ; 230,171,2 ; 217,95,2 ; 180,0,0]/255;

fidNsTimeVelocity=fopen(sprintf('%s/Ns_TimeVelocity.dat',path),'w');
fidNsTimeIntensity=fopen(sprintf('%s/Ns_TimeIntensity.dat',path),'w');

% loop through individual kymographs
% NOTE: each kymograph had 3 columns: 1) pos, 2) time, 3) normalized pos
for c=1:3:size(DATA,2)
    fprintf('%s\t(%d of %d)\n',DATA.Properties.VariableNames{c},(c+2)/3,size(DATA,2)/3)
    % flag to check first positive velocity
    velocityflag=false;
    
    % convert to numerical array
    M=table2array(DATA(:,c:c+2));
    
    if (M(1,3)~=0)
        fprintf('WARNING: Check number of headers\n')
    end
    
    % remove rows that contain all NaNs
    M(all(isnan(M),2),:) = [];
    
    % Convert time to seconds
    M(:,2)=M(:,2).*timeconversionfactor; % convert to seconds
    timeinterval=((M(2,2)-M(1,2)));
    
    % Convert normalized position to nm
    M(:,3)=M(:,3).*positionconversionfactor;
    
    % Smoothen normalized position (actual position, i.e. M(:,1), can be smoothened as well)
    MM=M;
    if ( floor( positionsmoothenwindowsizeinminutes * 60 / ((M(2,2)-M(1,2))) ) > 0 )
        MM(:,3)=smoothdata(M(:,3),'movmedian',positionsmoothenwindowsizeinminutes);
    end
    
    % read kymographs
    if (intensityanalysis)
        [TB,mapTB]=imread(sprintf('%s/TUB/kymograph/kymograph_%d/kymograph%d.tif',path,KYMOindexes(1,c+0),KYMOindexes(1,c+0)));
        [SSNA1,mapSSNA1]=imread(sprintf('%s/SSNA1/kymograph/kymograph_%d/kymograph%d.tif',path,KYMOindexes(1,c+1),KYMOindexes(1,c+1)));
        TBoverlay=zeros(size(TB));
        SSNA1intensityoverlaylattice=zeros(size(SSNA1));
        SSNA1intensityoverlaysolution=zeros(size(SSNA1));
        for t=1:size(M,1)
            TBoverlay(floor(M(t,2)/timeinterval),floor(M(t,1)*positionconversionfactor/160))=1;
        end
        imwrite(imfuse(TB,TBoverlay,'ColorChannels',[1 2 0]),sprintf('%s/UsedKymos/TUBoverlay-%d.tif',path,KYMOindexes(1,c+0)))
        
    end
    
    % find midpoint intervals for velocity fit
    startpoints=midpoints-(velocitywindowsizeinseconds/2)+timeinterval;
    endpoints=midpoints+(velocitywindowsizeinseconds/2);
    
    % calculate velocity with fit using sliding window
    %fidvelocityfit=fopen(sprintf('%s/Time_versus_%dFrameVelocityFit_K%d.dat',path,velocitywindowsizeinminutes,KYMOindexes(1,c+0)),'w');
    Vfit=NaN(size(midpoints,1),4);
    colorindex=NaN(size(midpoints,1),3);
    for t=1:size(midpoints,1)
        startframe=find(M(:,2)==startpoints(t,1));
        midpoint=find(M(:,2)==midpoints(t,1));
        endframe=find(M(:,2)==endpoints(t,1));
        
        if ( ~isempty(startframe) && ~isempty(midpoint) && ~isempty(endframe) )
            
            [fo,gof] = fit(MM(startframe:endframe,2), MM(startframe:endframe,3), ft, 'StartPoint', [1 1]);
            Coeffs = coeffvalues(fo);
            growthv=Coeffs(1);
            ci = confint(fo);
            growthvci=ci(:,1);
            growthverror=(growthvci(2)-growthvci(1))/2; % (Upper95%CI - Lower95%CI)/2
            
            Vfit(t,1)=MM(midpoint,2);
            Vfit(t,2)=growthv;
            Vfit(t,3)=growthverror;
            Vfit(t,4)=(MM(endframe,3)-MM(startframe,3))/(MM(endframe,2)-MM(startframe,2)); % empirically calculated velocity, DeltaX/DeltaT
            
            if ( (~velocityflag) && (Vfit(t,2)>0) )
                velocityflag=true;
            end
            
            if (velocityflag)
                
                % write data: 1) time in seconds, 2) velocity in nm/s, 3) velocity error in nm/s 4) empirical velocity in nm/s
                %fprintf(fidvelocityfit,'%f\t%f\t%f\t%f\n',Vfit(t,1),Vfit(t,2),Vfit(t,3),Vfit(t,4));
                
                % assign colors
                if (Vfit(t,2)<=0)
                    colorindex(t,:)=colormapmatrix(1,:);
                elseif ( (Vfit(t,2)>0) && (Vfit(t,2)<=2.5) )
                    colorindex(t,:)=colormapmatrix(2,:);
                elseif ( (Vfit(t,2)>2.5) && (Vfit(t,2)<=5) )
                    colorindex(t,:)=colormapmatrix(3,:);
                elseif ( (Vfit(t,2)>5) && (Vfit(t,2)<=7.5) )
                    colorindex(t,:)=colormapmatrix(4,:);
                elseif ( (Vfit(t,2)>7.5) && (Vfit(t,2)<=10) )
                    colorindex(t,:)=colormapmatrix(5,:);
                elseif ( (Vfit(t,2)>10) && (Vfit(t,2)<=12.5) )
                    colorindex(t,:)=colormapmatrix(6,:);
                elseif ( (Vfit(t,2)>12.5) && (Vfit(t,2)<=15) )
                    colorindex(t,:)=colormapmatrix(7,:);
                elseif (Vfit(t,2)>15)
                    colorindex(t,:)=colormapmatrix(8,:);
                end
                
                if (plotflag)
                    figure(1)
                    plot(MM(startframe:endframe,2)/60,MM(startframe:endframe,3)/1000,'-','Color',colorindex(t,:))
                    hold on
                end
                
                % filter out windows including shrinkage
                if (Vfit(t,2)>-2.5)
                    TrackVelocity(((c+2)/3),t)=Vfit(t,2);
                    TrackVelocityError(((c+2)/3),t)=Vfit(t,3);
                    
                    
                    % perform intensity analysis
                    if (intensityanalysis)
                        sumI=0;
                        counter=0;
                        for tt=startframe:endframe
                            sumflag=false;
                            row=floor(M(tt,2)/timeinterval);
                            col=floor(M(tt,1)*positionconversionfactor/160);
                            if ((M(end,1)-M(1,1))<0) % plus end grows towards left
                                cols=col:col+sumwidth-1;
                                bgcols=cols-(sumwidth+3);
                                if (min(bgcols)>0)
                                    sumflag=true;
                                end
                            else
                                cols=col-(sumwidth-1):col;
                                bgcols=cols+(sumwidth+3);
                                if (max(bgcols)<=size(SSNA1,2))
                                    sumflag=true;
                                end
                            end
                            %SSNA1copy(row,cols)=max(max(SSNA1));
                            if (sumflag)
                                aveBG=mean(SSNA1(row,bgcols));
                                Ibgsub=mean(SSNA1(row,cols))-aveBG;
                                sumI=sumI+Ibgsub;
                                counter=counter+1;
                                SSNA1intensityoverlaylattice(row,cols)=1;
                                SSNA1intensityoverlaysolution(row,bgcols)=1;
                            end
                        end
                        
                        if (counter>0)
                            TrackIntensity(((c+2)/3),t)=sumI/counter;
                            if (plotflag)
                                figure(3)
                                plot(Vfit(t,2),TrackIntensity(((c+2)/3),t),'.','Color',colorindex(t,:))
                                hold on
                            end
                        end
                    end
                end
            end
        end
    end
    %fclose(fidvelocityfit);
    if (intensityanalysis)
        overlay=NaN;
        overlay=imfuse(SSNA1intensityoverlaylattice,SSNA1intensityoverlaysolution,'ColorChannels',[2 1 0]);
        imwrite(overlay,sprintf('%s/UsedKymos/SSNA1overlay-%d.tif',path,KYMOindexes(1,c+0)))
    end
end

% Beahive plot
for t=1:size(midpoints,1)
    if (~all(isnan(TrackVelocity(:,t))))
        % The function that is called produces a beahive plot for a dataset with
        % the same x-axis values. y-axis values are unaltered, but x-axis values
        % are slightly altered to produce plot similar to a beehive. The output of
        % this function is x-y pairs of input dataset.
        [Tout,Vout]=BeeHive_20210107(midpoints(t,1)/60,TrackVelocity(:,t),Nbins,xwidth);
        if ( (intensityanalysis) && (~all(isnan(TrackIntensity(:,t)))) )
            [ToutTversusI,IoutTversusI]=BeeHive_20210107(midpoints(t,1)/60,TrackIntensity(:,t),Nbins,xwidth);
        end
        % assign colors and plot
        for c=1:size(Vout,1)
            
            if (Vout(c,1)<=0)
                colorindex=colormapmatrix(1,:);
            elseif ( (Vout(c,1)>0) && (Vout(c,1)<=2.5) )
                colorindex=colormapmatrix(2,:);
            elseif ( (Vout(c,1)>2.5) && (Vout(c,1)<=5) )
                colorindex=colormapmatrix(3,:);
            elseif ( (Vout(c,1)>5) && (Vout(c,1)<=7.5) )
                colorindex=colormapmatrix(4,:);
            elseif ( (Vout(c,1)>7.5) && (Vout(c,1)<=10) )
                colorindex=colormapmatrix(5,:);
            elseif ( (Vout(c,1)>10) && (Vout(c,1)<=12.5) )
                colorindex=colormapmatrix(6,:);
            elseif ( (Vout(c,1)>12.5) && (Vout(c,1)<=15) )
                colorindex=colormapmatrix(7,:);
            elseif (Vout(c,1)>15)
                colorindex=colormapmatrix(8,:);
            end
            
            if (plotflag)
            figure(2)
            plot(Tout(c,1),Vout(c,1),'.','Color',colorindex)
            hold on
            end
            
            if ( (intensityanalysis) && (~all(isnan(TrackIntensity(:,t)))) )
                if (plotflag)
                figure(4)
                plot(ToutTversusI(c,1),IoutTversusI(c,1),'.','Color',colorindex)
                hold on
                end
            end
            
        end
        fprintf(fidNsTimeVelocity,'%d, ',sum(~isnan(Vout)));
        if (intensityanalysis)
            fprintf(fidNsTimeIntensity,'%d, ',sum(~isnan(IoutTversusI),1));
        end
        medianV=nanmedian(Vout);
        if (plotflag)
        figure(2)
        plot([(midpoints(t,1)/60)-(xwidth)/2 (midpoints(t,1)/60)+(xwidth)/2],[medianV medianV],'r-','LineWidth',1)
        hold on
        end
        
        fidcolumns=fopen(sprintf('%s/VelocityColumnsAt-%dmin.dat',path,(midpoints(t,1)/60)),'w');
        for row=1:size(Vout,1)
            fprintf(fidcolumns,'%f\n',Vout(row,1));
        end
        fclose(fidcolumns);
    end
end 

fprintf(fidNsTimeVelocity,'\n');
fclose(fidNsTimeVelocity);
if (intensityanalysis)
    fprintf(fidNsTimeIntensity,'\n');
    fclose(fidNsTimeIntensity);
end

% this is both a double check and also to plot median as dots
medianV=nanmedian(TrackVelocity);

if (plotflag)
fig=figure(1);
xlabel('time (min)')
ylabel(sprintf('smooth position (um) (%d min smoothening)',positionsmoothenwindowsizeinminutes))
xlim([0 maxtimeinseconds/60]);
xticks(1:2:(maxtimeinseconds/60)+1);
ylim([0 22]);
yticks(0:5:22);
colormap(colormapmatrix);
colorbar('Ticks',[0:1/8:1],'TickLabels',{'<0','0','2.5','5','7.5','10','12.5','15','>15nm/s'})
hold off
saveas(fig,sprintf('%s/Time-Position.png',path));

fig=figure(2);
plot(midpoints/60,medianV,'ro','MarkerFace','r')
xlabel('time (min)')
ylabel('fit velocity (nm/s)')
xlim([0 maxtimeinseconds/60]);
xticks(1:2:(maxtimeinseconds/60)+1);
ylim([-3 21]);
yticks(-5:5:25);
title(sprintf('median'))
colormap(colormapmatrix);
colorbar('Ticks',[0:1/8:1],'TickLabels',{'<0','0','2.5','5','7.5','10','12.5','15','>15nm/s'})
hold off
saveas(fig,sprintf('%s/Time-Velocity.png',path));

if (intensityanalysis)
    fig=figure(3);
    xlabel('fit velocity (nm/s)')
    ylabel('average intensity (a.u.)')
    xlim([-3 21]);
    xticks(-5:5:25);
    ylim([-500 8500]);
    yticks(-500:500:8500);
    colormap(colormapmatrix);
    colorbar('Ticks',[0:1/8:1],'TickLabels',{'<0','0','2.5','5','7.5','10','12.5','15','>15nm/s'})
    hold off
    saveas(fig,sprintf('%s/Velocity-Intensity.png',path));
    
    fig=figure(4);
    xlabel('time (min)')
    ylabel('average intensity (a.u.)')
    xlim([0 maxtimeinseconds/60]);
    xticks(1:2:(maxtimeinseconds/60)+1);
    ylim([-500 8500]);
    yticks(-500:500:8500);
    colormap(colormapmatrix);
    colorbar('Ticks',[0:1/8:1],'TickLabels',{'<0','0','2.5','5','7.5','10','12.5','15','>15nm/s'})
    hold off
    saveas(fig,sprintf('%s/Time-Intensity_Beehive.png',path));
    
end
end


end




