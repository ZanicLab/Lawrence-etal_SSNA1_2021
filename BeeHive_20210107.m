% Created by Goker Arpag on 7 January 2021
% Last edited by GA on 10 Feb 2021 to add description
% This function takes y-axis values for a given x-axis value to produce a
% plot similar to a beehive.
% New x-axis values are produced by binning and sorting y-axis values which
% are not altered.
% Inputs:
% 1) single x-axis value.
% 2) array containing y-axis values
% 3) number of bins to group y-axis values
% 4) the width that the x-axis values will be altered to produce a plot
% similar to a beehive
% The outputs are y-axis values with corresponding new x-axis values. 
% To confirm that the y-values are not altered, do the following example in
% matlab command line and compare mean and standard deviation:
% >> Yvals=100*rand(50,1);
% >> Xval=5;
% >> [Xout,Yout]=BeeHive_20210107(Xval,Yvals,10,2);
% >> mean(Yvals)
% >> mean(Yout)
% >> std(Yvals)
% >> std(Yout)

function [Xout,Yout]=BeeHive_20210107(xval,Yarray,Nbins,Xwidth)
clearvars -except xval Yarray Nbins Xwidth

[N,edges] = histcounts(Yarray,Nbins);
Xout=NaN(size(Yarray,1),1);
Yout=NaN(size(Yarray,1),1);

maxN=max(N);
xspread=Xwidth/maxN;

%colorindex=jet(Nbins);
counter=0;
for b=1:Nbins
    Ysmall=sort(Yarray( (Yarray>=edges(b))&(Yarray<edges(b+1)) ));
    if (N(b)>2)
        if (mod(N(b),2)==1)
            midpoint=Ysmall(1,1);
            Y2D=reshape(Ysmall(2:end),2,[]);
        else
            midpoint=[];
            Y2D=reshape(Ysmall(1:end),2,[]);
        end
        leftarray=sort(Y2D(1,:),'descend');
        rightarray=sort(Y2D(2,:),'ascend');
        
        Ysmallarranged=transpose([leftarray midpoint rightarray]);
    else
        Ysmallarranged=Ysmall;
    end
    
    for i=1:N(b)
        offset=xspread*(N(b)+1)/2;
        Xnew=xval+(xspread*i)-offset;
        counter=counter+1;
        Xout(counter)=Xnew;
        Yout(counter)=Ysmallarranged(i);
%         plot(Xnew,Ysmallarranged(i),'o','Color',colorindex(b,:),'MarkerFaceColor',colorindex(b,:))
%         hold on
    end
    
end






end