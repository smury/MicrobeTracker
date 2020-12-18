function [kgraph] = kymograph2(clist,cell,varargin)
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
% 'absolute' - plot the kymograph in absolute spatial coordinated (pixels
% or microns, if provided)
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
%klength = 50;
normalize = false;
normarea = false;
normvolume = false;
numcount = 1;
framefactor = 1;
framefactorset = false;
display = true;
integrate = false;
absolute = false;
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
    elseif ischar(varargin{i}) && isequal(varargin{i},'absolute')
        absolute = true;
    elseif ischar(varargin{i}) && length(varargin{i})>6 && isequal(varargin{i}(1:6),'signal')
        sfield = varargin{i};
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
frange = [];


flag=0;
for frame=max(1,minframe):min(length(clist),maxframe)
    if length(clist{frame})>=cell && ~isempty(clist{frame}) && isfield(clist{frame}{cell},sfield)
        if flag==0 firstframe=frame; flag=1;end
        divisions=clist{frame}{cell}.divisions;
        descendants=clist{frame}{cell}.descendants;
        endframe=frame;
    end
end
if ~isempty(divisions)
endframe=min(endframe,divisions(divisions>firstframe)-1);%the endframe or the first divisions in the range, whichever is first
end

if absolute
    xmaxmax = 0;
    for frame=firstframe:endframe
            xmaxmax = max(xmaxmax,clist{frame}{cell}.lengthvector(end));
            if xmaxmax==clist{frame}{cell}.lengthvector(end)
            xlastmin=clist{frame}{cell}.lengthvector(1);
            klength=length(clist{frame}{cell}.lengthvector);
            end
    end
end

for frame=firstframe:endframe
        eval(['signal=clist{frame}{cell}.' sfield ';'])
        if isempty(signal), break; end
        x = clist{frame}{cell}.lengthvector;
        xmin = clist{frame}{cell}.lengthvector(1);
        xmax = clist{frame}{cell}.lengthvector(end);
        %x = [0;x;xmax];
        delta = 1E-10;
        signal = signal+delta;
        if absolute
            x = x+(xmaxmax-xmax)/2;
            xmax = xmaxmax;
            xmin=xlastmin;
        end
        
        if integrate
            isignal = integinterp(x/xmax,signal,klength);
            if normalize
                isignal = isignal./sum(isignal(~isnan(isignal)));
            end
            if normarea
                isignal = isignal./(delta+integinterp(x/xmax,[0;clist{frame}{cell}.steparea;0],klength));
            elseif normvolume
                isignal = isignal./integinterp(x/xmax,[delta;clist{frame}{cell}.stepvolume;delta],klength);
            end
        else
            if normalize
                signal = signal./sum(signal(~isnan(signal)));
            end
            if normarea
                %signal = signal./[delta;clist{frame}{cell}.steparea;delta];%bad
                %line
                signal = [signal(1:end)./clist{frame}{cell}.steparea];
            elseif normvolume
                %signal =
                %signal./[delta;clist{frame}{cell}.stepvolume;delta];%bad
                %line
                signal = [signal(1:end)./clist{frame}{cell}.stepvolume];
            end
            ix = linspace(xmin/xmax,1,klength);
      
            isignal = interp1(x/xmax,signal,ix,'linear',0);%value outside the range are 0 - needed for absolute option

        end
%         if normalize %moved earlier
%             isignal = isignal./sum(isignal(~isnan(isignal)));
%         end
        kgraph = [kgraph isignal'];
        frange = [frange frame];
end

if display
    if ~isempty(frange)
        if absolute
            kgraph2 = 63*(kgraph-min(min(kgraph(kgraph~=0))))/(max(max(kgraph))-min(min(kgraph(kgraph~=0))))+2;
            kgraph2(kgraph==0) = 1;
            %kgraph2 = kgraph/max(max(kgraph))*62+2;
            %kgraph2(kgraph2==2) = 1;
            colormap([[0.8 0.8 0.8];jet])
        else
            kgraph2 = 63*(kgraph-min(min(kgraph)))/(max(max(kgraph))-min(min(kgraph)))+1;%scale to fit whole 64 value colormap
            %kgraph2 = kgraph/max(max(kgraph))*63+1;
            colormap jet;
        end
        if framefactorset && absolute && ~isempty(pix2mu)
            image((frange-1)*framefactor,linspace(0,1,klength)*xmaxmax*pix2mu,kgraph2);
            xlabel('Time','FontSize',14)
            ylabel('Position, um','FontSize',14)
            xlim([frange(1)-1.5 frange(end)-0.5]*framefactor)
        elseif ~framefactorset && absolute && ~isempty(pix2mu)
            image(frange*framefactor,linspace(0,1,klength)*xmaxmax*pix2mu,kgraph2);
            xlabel('Frame','FontSize',14)
            ylabel('Position, um','FontSize',14)
            xlim([frange(1)-0.5 frange(end)+0.5])
        elseif framefactorset && absolute && isempty(pix2mu)
            image((frange-1)*framefactor,linspace(0,1,klength)*xmaxmax,kgraph2);
            xlabel('Time','FontSize',14)
            ylabel('Position, pixels','FontSize',14)
            xlim([frange(1)-1.5 frange(end)-0.5]*framefactor)
        elseif ~framefactorset && absolute && isempty(pix2mu)
            image(frange*framefactor,linspace(0,1,klength)*xmaxmax,kgraph2);
            xlabel('Frame','FontSize',14)
            ylabel('Position, pixels','FontSize',14)
            xlim([frange(1)-0.5 frange(end)+0.5])
        elseif framefactorset && ~absolute
            image((frange-1)*framefactor,linspace(0,1,klength),kgraph2);
            xlabel('Time','FontSize',14)
            ylabel('Relative position','FontSize',14)
            xlim([frange(1)-1.5 frange(end)-0.5]*framefactor)
        elseif ~framefactorset && ~absolute
            image(frange*framefactor,linspace(0,1,klength),kgraph2);
            xlabel('Frame','FontSize',14)
            ylabel('Relative position','FontSize',14)
            xlim([frange(1)-0.5 frange(end)+0.5])
        end
        set(gca,'FontSize',12,'YDir','normal')
    else
        disp('No data for this cell')
    end
end

