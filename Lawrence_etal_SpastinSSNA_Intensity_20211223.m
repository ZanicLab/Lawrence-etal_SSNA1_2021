%%% Created by GA
%%% last edited by GA on 20211223
%%% Use this to quantify intensity for Spastin binding on SSNA1-coated and not coated
%%% microtubules in the presence of AMPPNP.

function Lawrence_etal_SpastinSSNA_Intensity_20211223(path,filename,pointsname,tstart,tend)
close all

fstart=(tstart*60/15)+1;
fend=(tend*60/15)+1;
Nframe=fend-fstart+1;

lengthcutoff=2;

Nchannel=3;

[filepath,name,ext] = fileparts(filename);

info=imfinfo(sprintf('%s/%s',path,filename));
Ntotalframes=size(info,1)/Nchannel;
if (Ntotalframes<fend)
    fprintf('Frame number is wrong\n')
    return;
end

% read data file for coordinates
TSSNACoated=load(sprintf('%s/Points-SSNACoatedMTs-%s.txt',path,pointsname));
UniqueCurveSSNACoated = (unique(TSSNACoated(:,1)));
NUniqueCurveSSNACoated=size(UniqueCurveSSNACoated,1);

TNotCoated=load(sprintf('%s/Points-NotCoatedMTs-%s.txt',path,pointsname));
UniqueCurveNotCoated = (unique(TNotCoated(:,1)));
NUniqueCurveNotCoated=size(UniqueCurveNotCoated,1);




for condition=1:2
    if (condition==1)
        T=TSSNACoated;
        NUniqueCurve=NUniqueCurveSSNACoated;
        tag=sprintf('SSNACoated');
    elseif (condition==2)
        T=TNotCoated;
        NUniqueCurve=NUniqueCurveNotCoated;
        tag=sprintf('NotCoated');
    end
    
    CurveLength=NaN(NUniqueCurve,1);
    WeightedMean=NaN(NUniqueCurve,3);
    WeightedSTD=NaN(NUniqueCurve,3);
    MaskAllCurves=zeros(info(1).Height,info(1).Width);
    
    fidweightedmean=fopen(sprintf('%s/WeightedMeans_%s_%dto%dmin_%s.dat',path,tag,tstart,tend,name),'w');
    fidweightedSTD=fopen(sprintf('%s/WeightedSTDs_%s_%dto%dmin_%s.dat',path,tag,tstart,tend,name),'w');
    
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
            
            idxmask=(Mask==1);
            
            FrameMean=NaN(Nframe,3);
            FrameWeight=NaN(Nframe,3);
            
            for frame=fstart:fend
                % read frames from each channel
                frameR=(frame*Nchannel)-(Nchannel-1);
                [MR,mapR]=imread(sprintf('%s/%s',path,filename),frameR);
                
                frameG=(frame*Nchannel)-(Nchannel-2);
                [MG,mapG]=imread(sprintf('%s/%s',path,filename),frameG);
                
                frameM=(frame*Nchannel)-(Nchannel-3);
                [MM,mapM]=imread(sprintf('%s/%s',path,filename),frameM);
                
                RI=double(MR(idxmask));
                GI=double(MG(idxmask));
                MI=double(MM(idxmask));
                
                
                FrameMean(frame-fstart+1,1)=mean(RI);
                FrameMean(frame-fstart+1,2)=mean(GI);
                FrameMean(frame-fstart+1,3)=mean(MI);
                FrameWeight(frame-fstart+1,1)=1/((std(RI)/sqrt(size(RI,1)))^2);
                FrameWeight(frame-fstart+1,2)=1/((std(GI)/sqrt(size(GI,1)))^2);
                FrameWeight(frame-fstart+1,3)=1/((std(MI)/sqrt(size(MI,1)))^2);
            end
            
            % calculate weighted mean and STD using averages from different frames
            for j=1:3
                WeightedMean(c,j)=sum(FrameMean(:,j).*FrameWeight(:,j))/sum(FrameWeight(:,j));
                WeightedSTD(c,j)=sqrt(var(FrameMean(:,j),FrameWeight(:,j)))*sqrt(Nframe/(Nframe-1)); % to correct for normalizing with N-1, instead of N
            end
            fprintf(fidweightedmean,'%d\t%f\t%f\t%f\n',c,WeightedMean(c,1),WeightedMean(c,2),WeightedMean(c,3));
            fprintf(fidweightedSTD,'%d\t%f\t%f\t%f\n',c,WeightedSTD(c,1),WeightedSTD(c,2),WeightedSTD(c,3));
        end
    end
    
    imwrite(MaskAllCurves,sprintf('%s/Mask_%s_%s.tif',path,tag,name));
    
    [Lcdf,L]=ecdf(CurveLength(:,1));
    fig=figure(900+c);
    plot(L,Lcdf)
    xlabel('length (um)')
    ylabel('CDF')
    pbaspect([1 1 1])
    exportgraphics(fig,sprintf('%s/LengthDistribution_%s_%s.png',path,tag,name));
    close(900+c);
    
    fclose(fidweightedmean);
    fclose(fidweightedSTD);
    
end



end