
%%% Created by GA
%%% last edited by GA on 20211223
%%% Use this to quantify intensity for Spastin severing on SSNA1-coated and not coated
%%% microtubules in the presence of ATP.


function Lawrence_etal_SpastinSSNA_Intensity_SingleFrame_20211223(path,name,tstart,tend)
close all
cutoffflag=0;
fstart=(tstart*60/15)+1;
fend=(tend*60/15)+1;
Nframe=fend-fstart+1;

lengthcutoff=2;

Nchannel=3;

info=imfinfo(sprintf('%s/BGsub-%s.tif',path,name));
Ntotalframes=size(info,1)/Nchannel;
if (Ntotalframes<fend)
    fprintf('Frame number is wrong\n')
    return;
end

% read data file for coordinates
TSSNACoated=load(sprintf('%s/Points-SSNACoatedMTs-%s.txt',path,name));
UniqueCurveSSNACoated = (unique(TSSNACoated(:,1)));
NUniqueCurveSSNACoated=size(UniqueCurveSSNACoated,1);

TNotCoated=load(sprintf('%s/Points-NotCoatedMTs-%s.txt',path,name));
UniqueCurveNotCoated = (unique(TNotCoated(:,1)));
NUniqueCurveNotCoated=size(UniqueCurveNotCoated,1);

IndividualMTdataSSNACoated=NaN(Nframe*NUniqueCurveSSNACoated,6);
IndividualMTdataNotCoated=NaN(Nframe*NUniqueCurveNotCoated,6);

for condition=1:2
    if (condition==1)
        T=TSSNACoated;
        NUniqueCurve=NUniqueCurveSSNACoated;
        tag=sprintf('SSNACoated');
        plotcolor=[1 0 0];
        SSNAcoatedCounter=0;
    elseif (condition==2)
        T=TNotCoated;
        NUniqueCurve=NUniqueCurveNotCoated;
        tag=sprintf('NotCoated');
        plotcolor=[0 0 1];
        NotcoatedCounter=0;
    end
    
    Time=NaN(Nframe,1);
    WeightedDataTB=NaN(Nframe,3);
    WeightedDataSpastin=NaN(Nframe,3);
    WeightedDataSSNA=NaN(Nframe,3);
    
    MaskAllCurves=zeros(info(1).Height,info(1).Width);
    
    fidweighteddataTB=fopen(sprintf('%s/WeightedMeansTBChannel-%s_%s_20211223.dat',path,tag,name),'w');
    fidweighteddataSpastin=fopen(sprintf('%s/WeightedMeansSpastinChannel-%s_%s_20211223.dat',path,tag,name),'w');
    fidweighteddataSSNA=fopen(sprintf('%s/WeightedMeansSSNAChannel-%s_%s_20211223.dat',path,tag,name),'w');
    for frame=fstart:fend
        CurveLength=NaN(NUniqueCurve,1);
        FrameMean=NaN(NUniqueCurve,3);
        FrameSEM=NaN(NUniqueCurve,3);
        FrameWeight=NaN(NUniqueCurve,3);
        % read frames from each channel
        frameR=(frame*Nchannel)-(Nchannel-1);
        [MR,mapR]=imread(sprintf('%s/BGsub-%s.tif',path,name),frameR);
        
        frameG=(frame*Nchannel)-(Nchannel-2);
        [MG,mapG]=imread(sprintf('%s/BGsub-%s.tif',path,name),frameG);
        
        frameM=(frame*Nchannel)-(Nchannel-3);
        [MM,mapM]=imread(sprintf('%s/BGsub-%s.tif',path,name),frameM);
        
        curvecounter=0;
        for c=1:NUniqueCurve
            idx=find(((T(:,1))==c));
            Xcoord=round(T(idx,2));
            Ycoord=round(T(idx,3));
            Mask=zeros(info(1).Height,info(1).Width);
            length=0.16;
            lcount=0;
            prevx=Xcoord(1,1);
            prevy=Ycoord(1,1);
            for s=1:size(Xcoord,1)
                if ( (prevx==Xcoord(s,1)) && (prevy==Ycoord(s,1)) && (s>1) )
                    continue
                end
                if  ( (Xcoord(s,1)>(info(1).Width/2)) && (Xcoord(s,1)<=info(1).Width) && (Ycoord(s,1)<=info(1).Height) ) % analyze only the right side of the FOV
                    Mask(Ycoord(s,1),Xcoord(s,1))=1;
                    lcount=lcount+1;
                    deltaL=0.16*sqrt( ((Xcoord(s,1)-prevx)^2) + ((Ycoord(s,1)-prevy)^2)   );
                    length=length+deltaL;
                    prevx=Xcoord(s,1);
                    prevy=Ycoord(s,1);
                    
                else
                    prevx=Xcoord(s,1);
                    prevy=Ycoord(s,1);
                end
                
                
            end
            
            MaskAllCurves=MaskAllCurves+Mask;
            CurveLength(c,1)=length;
            if (length>lengthcutoff)
                imwrite(Mask,sprintf('%s/Mask_%d.tif',path,c));
                idxmask=(Mask==1);
                
                RI=double(MR(idxmask));
                GI=double(MG(idxmask));
                MI=double(MM(idxmask));
                
                
                if ( (~isnan(1/(std(MI)^2))) && (~isnan(1/(std(MI)^2))) && (~isnan(1/(std(MI)^2))) )
                    curvecounter=curvecounter+1;
                    FrameMean(c,1)=mean(RI);
                    FrameMean(c,2)=mean(GI);
                    FrameMean(c,3)=mean(MI);
                    FrameSEM(c,1)=std(RI)/sqrt(size(RI,1));
                    FrameSEM(c,2)=std(GI)/sqrt(size(GI,1));
                    FrameSEM(c,3)=std(MI)/sqrt(size(MI,1));
                    FrameWeight(c,1)=1/(FrameSEM(c,1)^2);
                    FrameWeight(c,2)=1/(FrameSEM(c,2)^2);
                    FrameWeight(c,3)=1/(FrameSEM(c,3)^2);

                    if (condition==1)
                        SSNAcoatedCounter=SSNAcoatedCounter+1;
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,1)=frame;
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,2)=c;
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,3)=FrameMean(c,1);
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,4)=FrameSEM(c,1);
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,5)=size(RI,1);
                        IndividualMTdataSSNACoated(SSNAcoatedCounter,6)=length;
                    else
                        NotcoatedCounter=NotcoatedCounter+1;
                        IndividualMTdataNotCoated(NotcoatedCounter,1)=frame;
                        IndividualMTdataNotCoated(NotcoatedCounter,2)=c;
                        IndividualMTdataNotCoated(NotcoatedCounter,3)=FrameMean(c,1);
                        IndividualMTdataNotCoated(NotcoatedCounter,4)=FrameSEM(c,1);
                        IndividualMTdataNotCoated(NotcoatedCounter,5)=size(RI,1);
                        IndividualMTdataNotCoated(NotcoatedCounter,6)=length;
                    end
                end

            end
        end
            % calculate weighted mean and STD using averages from different frames
            FrameMean(all(isnan(FrameMean),2),:) = [];
            FrameSEM(all(isnan(FrameSEM),2),:) = [];
            FrameWeight(all(isnan(FrameWeight),2),:) = [];

            WeightedMean=NaN(1,3);
            WeightedSTD=NaN(1,3);
            for j=1:3
                WeightedMean(1,j)=sum(FrameMean(:,j).*FrameWeight(:,j))/sum(FrameWeight(:,j));
                WeightedSTD(1,j)=sqrt(var(FrameMean(:,j),FrameWeight(:,j)))*sqrt(curvecounter/(curvecounter-1)); % to correct for normalizing with N-1, instead of N
            end
            Time(frame,1)=(frame-2)*15/60;
            fprintf(fidweighteddataTB,'%f\t%f\t%f\t%d\n',Time(frame,1),WeightedMean(1,1),WeightedSTD(1,1),curvecounter);
            fprintf(fidweighteddataSpastin,'%f\t%f\t%f\t%d\n',Time(frame,1),WeightedMean(1,2),WeightedSTD(1,2),curvecounter);
            fprintf(fidweighteddataSSNA,'%f\t%f\t%f\t%d\n',Time(frame,1),WeightedMean(1,3),WeightedSTD(1,3),curvecounter);
        
        WeightedDataTB(frame,1)=WeightedMean(1,1);
        WeightedDataTB(frame,2)=WeightedSTD(1,1);
        WeightedDataTB(frame,3)=curvecounter;
        WeightedDataSpastin(frame,1)=WeightedMean(1,2);
        WeightedDataSpastin(frame,2)=WeightedSTD(1,2);
        WeightedDataSpastin(frame,3)=curvecounter;
        WeightedDataSSNA(frame,1)=WeightedMean(1,3);
        WeightedDataSSNA(frame,2)=WeightedSTD(1,3);
        WeightedDataSSNA(frame,3)=curvecounter;
    end
    
    tbnormalize=mean(WeightedDataTB(2,1));
    WeightedDataTBNormalized=NaN(Nframe,2);
    WeightedDataTBNormalized(:,1)=WeightedDataTB(:,1)/tbnormalize;
    WeightedDataTBNormalized(:,2)=WeightedDataTB(:,2)/tbnormalize;
    WeightedDataTBNormalized(:,3)=WeightedDataTB(:,3);
    
    if (condition==1)
        IndividualMTdataSSNACoated(:,3)=IndividualMTdataSSNACoated(:,3)/tbnormalize;
        IndividualMTdataSSNACoated(:,4)=IndividualMTdataSSNACoated(:,4)/tbnormalize;
    else
        IndividualMTdataNotCoated(:,3)=IndividualMTdataNotCoated(:,3)/tbnormalize;
        IndividualMTdataNotCoated(:,4)=IndividualMTdataNotCoated(:,4)/tbnormalize;
    end
    
    figure(1)
    errorbar(Time(2:end,1),WeightedDataTBNormalized(2:end,1),WeightedDataTBNormalized(2:end,2),'s','color',plotcolor,'MarkerFaceColor',plotcolor)
    hold on
    figure(2)
    errorbar(Time(2:end,1),WeightedDataTBNormalized(2:end,1),WeightedDataTBNormalized(2:end,2)./sqrt(WeightedDataTBNormalized(2:end,3)),'s','color',plotcolor,'MarkerFaceColor',plotcolor)
    hold on
    
    
    imwrite(MaskAllCurves,sprintf('%s/Mask_%s_%s.tif',path,tag,name));

    fclose(fidweighteddataTB);
    fclose(fidweighteddataSpastin);
    fclose(fidweighteddataSSNA);
    
    if (condition==2) 
       for frame=2:size(WeightedDataTBNormalized,1)-1
          if ( (cutoffflag==0) && (WeightedDataTBNormalized(frame,1)<0.2) && (WeightedDataTBNormalized(frame+1,1)<0.2) ) 
              fprintf('')
              cutoffflag=1;
              SSNAcoatedidx=IndividualMTdataSSNACoated(:,1)==frame;
              Notcoatedidx=IndividualMTdataNotCoated(:,1)==frame;
              writematrix(IndividualMTdataSSNACoated(SSNAcoatedidx,:),sprintf('%s/IndividualTBIntensitiesAtCuroff-SSNAcoated-%s.dat',path,name),'Delimiter','tab');
              writematrix(IndividualMTdataNotCoated(Notcoatedidx,:),sprintf('%s/IndividualTBIntensitiesAtCuroff-Notcoated-%s.dat',path,name),'Delimiter','tab');
          end
       end
        
    end
    
end

fig=figure(1);
pbaspect([1 1 1])
ylim([ 0 inf])
xlabel('time (min)')
ylabel('normalized TB intensity ± SD (a.u.)')
legend('SSNA1 coated','not coated')
exportgraphics(fig,sprintf('%s/TBintensities_WeightedMean-WeightedSD-%s_20211223.png',path,name));
saveas(fig,sprintf('%s/TBintensities_WeightedMean-WeightedSD-%s_20211223.fig',path,name));

fig=figure(2);
pbaspect([1 1 1])
ylim([ 0 inf])
xlabel('time (min)')
ylabel('normalized TB intensity ± SEM (a.u.)')
legend('SSNA1 coated','not coated')
exportgraphics(fig,sprintf('%s/TBintensities_WeightedMean-WeightedSEM-%s_20211223.png',path,name));
saveas(fig,sprintf('%s/TBintensities_WeightedMean-WeightedSEM-%s_20211223.fig',path,name));






end