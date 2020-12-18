function [kgraph2, signalstoplot,firstframe,endframe] = kymograph2a(clist,cell,varargin)
% Syntax:
% kymograph(cellList,cell);
% kymograph(cellList,cell,[nbins],[frame2time],[pix2mu],[framerange]);
% kymograph(cellList,cell,...,[signal]);
% kymograph(cellList,cell,...,'area');
% kymograph(cellList,cell,...,'volume');
% kymograph(cellList,cell,...,'normalize');
% kymograph(cellList,cell,...,'nodisp');
% kgraph = kymograph(cellList,cell,...);
%
% This function produces a kymograph of the signal profiles for a single
% cell in a list of obtained by MicrobeTracker.
%
% <cellList> - input cell list (the other possible input, which can be used
%     instead of the curvlist). If neither curvlist nor cellList is
%     supplied, the program will request to open a file with the cellList.
% <cell> - the cell to produce a kymorgaph for.
% <nbins> (the arguments in square brackets are optional) - number of bins
%     in the kymograph (the profile is resampled into a fixed number of
%     bins and is therefore produced in the relative coordinates along the
%     cell to compensate for cell growth; default: 50).
% <frame2time> - conversion factor from frames to time. If not supplied,
%     the kymograph is plotted vs. frame. To supply frame2time but not
%     nbins, supply an empty array ([]) for nbins.
% <pix2mu> - conversion factor from pixels to microns. Only works when
%     'absolute' parameter is provided. To supply pix2mu but not nbins or
%     frame2time, supply empty arrays ([]) for nbins and frame2time.
% <framerange> - two-element array, indicating the starting and the final
%     frames of a range to use for the kymograph. If not supplied, the
%     program starts from the beginning and proceeds until the cell is
%     divided or the end of the list is reached.
% <signal> - the signal field which will be processed. Must be in single
%     quotes and must start with the word 'signal'. Default: 'signal1'.
% <framerange> - two-element array, indicating the starting and the final
%     frames of a range to use for the kymograph. If not supplied, the
%     program starts from the beginning and proceeds until the cell is
%     divided or the end of the list is reached.
% 'area' - normalize the profile by the area of each segment producing the
%     average image intensity in each segment.
% 'volume' - normalize the profile by the estimated volume of each segment
%     producing the estimate of signal concentration in each segment.
% 'normalize' - normalize each profile by its sum to compensate for
%     photobleaching of the signal (if known to be constant).
% 'nodisp' - suppress displaying the results as a figure. Note, the
%     function does not produce a new figure by itself. Therefore, you
%     should use the figure command to avoid replotting an existing figure.
% 'integrate' - obtain the values in the kymograpg by integration instead
%     of interpolation of the raw data values.
% <kgraph> - the output kymograph matrix with the first dimension being
%     the relative position in the cell and the second one being the frame
%     number in the range.

if ~iscell(clist)
    disp('Incorrect cell list format')
    return
end
if ~isnumeric(cell)
    disp('Incorrect cell number')
    return
end
minframe = 1;
maxframe = Inf;
sfield = 'signal1';
spotsfield='spots';
%klength = 50;
normalize = false;
normarea = false;
normvolume = false;
numcount = 1;
framefactor = 1;
framefactorset = false;
display = true;
integrate = false;
scale= false;
pix2mu = [];
for i=1:length(varargin)
    if ischar(varargin{i}) && isequal(varargin{i},'normalize')
        normalize = true;
    elseif ischar(varargin{i}) && isequal(varargin{i},'area')
        if ~normvolume, normarea = true; end
    elseif ischar(varargin{i}) && isequal(varargin{i},'volume')
        if ~normarea, normvolume = true; end
    elseif ischar(varargin{i}) && isequal(varargin{i},'nodisp')
        display = false;
    elseif ischar(varargin{i}) && isequal(varargin{i},'integrate')
        integrate = true;
    elseif ischar(varargin{i}) && isequal(varargin{i},'scale')
        scale = true;
    elseif ischar(varargin{i}) && length(varargin{i})>6 && (isequal(varargin{i}(1:6),'signal') || isequal(varargin{i}(1:8),'nucleoid') )
        sfield = varargin{i};
    elseif ischar(varargin{i}) && length(varargin{i})>4 && isequal(varargin{i}(1:5),'spots')
        spotsfield = varargin{i};
    elseif isnumeric(varargin{i}) && length(varargin{i})<=1 && numcount==1
        if length(varargin{i})==1, klength = varargin{i}; end
        numcount = numcount+1;
    elseif isnumeric(varargin{i}) && length(varargin{i})<=1 && numcount==2
        if length(varargin{i})==1, framefactor = varargin{i}; framefactorset = true; end
        numcount = numcount+1;
    elseif isnumeric(varargin{i}) && length(varargin{i})<=1 && numcount==3
        if length(varargin{i})==1, pix2mu = varargin{i}; end
    elseif isnumeric(varargin{i}) && length(varargin{i})==2
        minframe = varargin{i}(1);
        maxframe = varargin{i}(2);
    end
end

kgraph = [];
kgraph2 = [];
signalstoplot=[];
frange = [];


flag=0;
divisions=[];
firstframe=[];
endframe=[];
for frame=max(1,minframe):min(length(clist),maxframe)
    if ~isempty(clist{frame}) && length(clist{frame})>=cell && ~isempty(clist{frame}{cell}) && length(clist{frame}{cell}.mesh)>4 && isfield(clist{frame}{cell},sfield)
        if flag==0
            firstframe=frame;
            flag=1;
        end
        divisions=clist{frame}{cell}.divisions;
        endframe=frame;
    end
end

if isempty(firstframe)
    disp('Cell not found')
    return
end

d=divisions(divisions>firstframe);
if isempty(d)
    endframe;
elseif length(d)==1
    endframe;
else
    endframe=d(2)-1; %the endframe or before the second divisions in the range, whichever is first
end




xmaxmax = 0;
xmaxmaxmother=0;
for frame=firstframe:endframe
    L=clist{frame}{cell}.lengthvector(end);
    if ~isempty(clist{frame}{cell}.divisions) && clist{frame}{cell}.divisions(end)>firstframe
        L=L+clist{frame}{clist{frame}{cell}.descendants(end)}.lengthvector(end)+2;
    end
    xmaxmax = max(xmaxmax,L);
    if xmaxmax==L
        xmaxmaxmother=clist{frame}{cell}.lengthvector(end);
        klength=round(xmaxmax);
    end
end

allspots=[];
spotsframes=[];

i=1;
for frame=firstframe:endframe
    cells=cell;
    if ~isempty(clist{frame}{cell}.divisions) && clist{frame}{cell}.divisions(end)>firstframe
        cells=[cells clist{frame}{cell}.descendants(end)];
    end
    isignal=0;
    %signalslist=[];
    j=1;
    for Cell=cells
        
        eval(['signal=clist{frame}{Cell}.' sfield ';'])
        if isempty(signal), break; end
        x = clist{frame}{Cell}.lengthvector;
        xmin = clist{frame}{Cell}.lengthvector(1);
        xmax = clist{frame}{Cell}.lengthvector(end);
        if isfield(clist{frame}{Cell},spotsfield)
            spots= clist{frame}{Cell}.(spotsfield).l;
        else
            spots=[];
        end
        %x = [0;x;xmax];
        delta = 1E-10;
        signal = signal+delta;
        if length(cells)==1
            x = x+(xmaxmax-xmax)/2;
            spots=spots+(xmaxmax-xmax)/2;
            xmax = xmaxmax;
        else if Cell==cells(1) %mother
                x = x+xmaxmaxmother-xmax;
                spots =spots + xmaxmaxmother-xmax;
                xmax = xmaxmax;
            else %daughter
                x=xmax-x;%invert order
                spots=xmax-spots;
                xmax=max(x);
                x = x+xmaxmaxmother-xmin+3;
                spots=spots+xmaxmaxmother-xmin+3;
                xmax = xmaxmax;
            end
        end
        
        
        if integrate
            isignal = integinterp(x/xmax,signal,klength);
            if normalize
                isignal = isignal./sum(isignal(~isnan(isignal)));
            end
            if normarea
                isignalcell = isignal./(delta+integinterp(x/xmax,[0;clist{frame}{Cell}.steparea;0],klength));
            elseif normvolume
                isignalcell = isignal./integinterp(x/xmax,[delta;clist{frame}{Cell}.stepvolume;delta],klength);
            end
        else
            if normarea
                %signal = signal./[delta;clist{frame}{cell}.steparea;delta];%bad
                %line
                totalsignal=sum(signal);
                signal = signal./clist{frame}{Cell}.steparea;
                signal=smooth(x,signal,7,'sgolay',2);
                signal=signal/(totalsignal/clist{frame}{Cell}.area);
                
            elseif normvolume
                %signal =
                %signal./[delta;clist{frame}{cell}.stepvolume;delta];%bad
                %line
                signal = signal./clist{frame}{Cell}.stepvolume;
            end
            ix = linspace(0,1,klength);
            spotpos=[];
            if ~isempty(spots)
                for k=1:length(spots)
                    [~,spotpos(k)]=min(abs(ix-spots(k)/xmax));
                end
            end
            allspots=[allspots spotpos];
            spotsframes=[spotsframes frame*ones(1,length(spotpos))];
            signalslist{1,j}=x';
            signalslist{2,j}=signal';
            j=j+1;
            isignalcell = interp1(x/xmax,signal,ix,'linear');
            
            if normalize
                isignalcell = isignalcell./sum(isignalcell(~isnan(isignalcell)));
            elseif scale
                isignalcell=isignalcell/max(isignalcell);
            end
        end

        
        isignal=isignal+isignalcell;
    end
    signalstoplot{i}=signalslist;
    i=i+1;
    kgraph = [kgraph isignal'];
    frange = [frange frame];
    
end

if display
    if ~isempty(frange)
        %kgraph2 = 63*(kgraph-min(min(kgraph(kgraph~=0))))/(max(max(kgraph))-min(min(kgraph(kgraph~=0))))+2;
        %kgraph2 = 63*(kgraph)/(max(max(kgraph)))+2;
        %kgraph2(kgraph==0) = 1;
        kgraph2=kgraph;
        %kgraph2(kgraph2==0) = NaN;
        %kgraph2 = kgraph/max(max(kgraph))*62+2;
        %kgraph2(kgraph2==2) = 1;
        
        
        if framefactorset && ~isempty(pix2mu)
            image((frange-1)*framefactor,linspace(0,1,klength)*xmaxmax*pix2mu,kgraph2);
            xlabel('Time','FontSize',14)
            ylabel('Position, um','FontSize',14)
            xlim([frange(1)-1.5 frange(end)-0.5]*framefactor)
        elseif ~framefactorset && ~isempty(pix2mu)
            image(frange*framefactor,linspace(0,1,klength)*xmaxmax*pix2mu,kgraph2);
            xlabel('Frame','FontSize',14)
            ylabel('Position, um','FontSize',14)
            xlim([frange(1)-0.5 frange(end)+0.5])
        elseif framefactorset && isempty(pix2mu)
            image((frange-1)*framefactor,linspace(0,1,klength)*xmaxmax,kgraph2);
            xlabel('Time','FontSize',14)
            ylabel('Position, pixels','FontSize',14)
            xlim([frange(1)-1.5 frange(end)-0.5]*framefactor)
        elseif ~framefactorset && isempty(pix2mu)
            h=imagesc(frange*framefactor,linspace(0,1,klength)*xmaxmax,kgraph2);%,'CDataMapping','scaled');
            %freezeColors;
            set(h,'AlphaData',~isnan(kgraph2));
            xlabel('Frame','FontSize',14)
            ylabel('Position, pixels','FontSize',14)
            xlim([frange(1)-0.5 frange(end)+0.5])
        end
        
        ax=gca;
        %colormap(ax,[[0.8 0.8 0.8];jet(256)]);
        %colormap(ax,[[1 1 1];parula(256)]);
        colorbar;
        %set(h,'alphadata',~isnan(kgraph2));
        
            hold on;
            %plot(spotsframes, allspots, 'x');
            hold off;
        
    end
    set(gca,'FontSize',12,'YDir','normal')
else
    disp('No data for this cell')
end
end

