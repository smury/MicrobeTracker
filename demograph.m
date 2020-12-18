

function demograph(varargin) %cellList,maxCellNum,maxCellLength,numPixelsMovingAverage,signal,frameNum,descriptor
%-----------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
%function demograph(varargin)
%
%@author:  Jason Hocking
%@date:    October 27, 2011
%@modified:  Ahmad Paintdakhi -- August 1, 2013
%@copyright 2011-2013 Yale University
%====================================================================================================
%**********output********:
%No output arguements required, the function oly plots information.
%**********Input********:
%cellList:    cellList structure
%maxCellNum:    maximum number of cells to be included for final demograph
%maxCellLength:    maximum length of cell
%numPixelsMovingAverage:    number of pixels to be used for the moving average.  This
%                            routine finds the segment where max intensity of the signal
%                            is located.
%signal:    an array for signal information.  For example, to use only signal 1 the
%            array should be [1,0], for signal 2 --> [0,1] and both signal 1 and
%           random signal 2 ---> [1,1].
%frameNum:  frame # to be used for analysis or [] vector (use all frames in
%           a dataset). 
%descriptor:    The descriptor value is a key for the type of demograph to be drawn.
%                The different keys are 'randomN','randomNasorted', 'constriction_noNormalization',
%                'sort_by_constriction','constriction'.
%Purpose:  script was designed to provide a colormap of relative segment intensities for every
%            cell in an asynchorous cellList, sorted by cell length in ascending order.
%====================================================================================================

if length(varargin) < 7 || length(varargin) > 7
    disp('Only 7 arguments are accepted');
    return;
end
if ~isstruct(varargin{1}) && ~iscell(varargin{1})
    disp('cellList must be a struct or cell array')
    return;
end
 
if length(varargin{5}) ~=2
   disp('signal should be a vector of length 2, such as [1,0] or [0,1] or [1,1]')
   return;
end

if ~ischar(varargin{7})
    disp(['descriptor needs to be a string such as ' 'randomN'])
    return;
end
cellList = varargin{1};
maxCellNum = varargin{2};
maxCellLength = varargin{3};
numPixelsMovingAverage = varargin{4};
signal = varargin{5};
frameNum = varargin{6};
descriptor = varargin{7};

signal1='signal1';
signal2='signal2';


warning('off','MATLAB:colon:nonIntegerIndex');
descriptorValues = {'randomN'
                    'randomNasorted'
                    'constriction_noNormalization'
                    'sort_by_constriction'
                    'constriction'};
%---------------------------------------------------------------------------------


if isempty(frameNum)
    frameList = 1:length(cellList);
else
    frameList = frameNum;
end

switch descriptor

    case 'randomNasorted'
           
            replacement=false;

            %%finds the maximum number of stepareas inside of a cell from the cellList

            maxsizelarray=[];
            n=0;
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    if isempty(cellList{frame}{cellNum}) || ...
                            length(cellList{frame}{cellNum}.mesh)<4 ...
                            ||~isfield(cellList{frame}{cellNum},signal1) ...
                            || isempty(cellList{frame}{cellNum}.(signal1)) ...
                            || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                        n=n+1;
                end
            end
            if n<=maxCellNum
                maxCellNum=n;
            end
            rand=randsample(n,maxCellNum,replacement);
            n=0;
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    if isempty(cellList{frame}{cellNum}) ...
                            || length(cellList{frame}{cellNum}.mesh)<4 ...
                            || ~isfield(cellList{frame}{cellNum},signal1) ...
                            || isempty(cellList{frame}{cellNum}.(signal1)) ...
                            || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                    n = n+1;
                    b=rand==n;
                   
                    if sum(b)~=1
                        continue
                    end
                 
                    maxsizelarray=[maxsizelarray length(cellList{frame}{cellNum}.lengthvector)];%#ok<AGROW>  
                end
            end

            %using the maxima from above, a matrix consiting of zeros is created to be
            %filled in by mesh intensities

            relintarray1=zeros(maxCellLength,maxCellNum);
            maxsizel=maxCellLength;

            % if maxCellLength > maxsizel;
            %     maxsizel=maxCellLength;
            % end
            maxsizel2 = ceil(maxsizel); if mod(maxsizel2,2)==0, maxsizel2=maxsizel2+1; end
            maxsizel2a = maxsizel2/2+0.5;

            n=0;
            passed=0;
            cellLength=[];

            %zeroarray is replaced with relative segment intensity data from the cell
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    place=1;%#ok

                    if isempty(cellList{frame}{cellNum}) ...
                        || length(cellList{frame}{cellNum}.mesh)<4 ...
                        ||~isfield(cellList{frame}{cellNum},signal1) ...
                        || isempty(cellList{frame}{cellNum}.(signal1)) ...
                        || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                    n = n+1;
                    b=rand==n;
                   
                    if sum(b)~=1
                        continue
                    end
                   
                    passed=passed+1;
                    if signal(1) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                    elseif signal(2) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    elseif sum(signal) == 2
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                       
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    else
                        disp('provide information in signal variable')
                        return;
                    end
                       
                    %%% A MOVING AVERAGE IS CALCULATED FOR EACH OF THE SEGMENTS TO FIND THE SINGLE BRIGHTEST SEGMENT AREA
                    cellList{frame}{cellNum}.meshavg=[];
                    if signal(1) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint1)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint1(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1;%#ok
                        end
                    elseif signal(2) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint2)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint2(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1;%#ok
                        end
                    end
                    %%WITH THE BRIGHTEST SEGMENT CALCULATED ABOVE WE CAN ORIENT THE
                    %%CELL SO THAT THE BRIGHTEST SEGMENT IS ON THE RIGHTS (i.e. WITH FtsZ BEING POLAR ON RIGHT(NEW POLE)
                    %%AND LARGER STALK CELL BIAS LETTING THE FtsZ RING BE ON THE RIGHT
                    %%AS WELL)
                    [~,maxavg]=max(cellList{frame}{cellNum}.meshavg);
                    if  maxavg<=length(cellList{frame}{cellNum}.meshavg)/2+1;
                        if signal(1) ==1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                        elseif signal(2) == 1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        elseif sum(signal) == 2
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        end
                    end
                   
                    k = floor(cellList{frame}{cellNum}.length/2);
                    temp = cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2;
                    if signal(1) ==1 && sum(signal) == 1
                        interpint1 = interp1(temp(1:length(cellList{frame}{cellNum}.relint1)),cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                        relintarray1(maxsizel2a-k:maxsizel2a+k,passed)=interpint1;
                    elseif signal(2) == 1 && sum(signal) == 1
                        interpint2 = interp1(temp(1:length(cellList{frame}{cellNum}.relint2)),cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                        relintarray2(maxsizel2a-k:maxsizel2a+k,passed)=interpint2;%#ok
                    elseif sum(signal) == 2
                        interpint1 = interp1(temp(1:length(cellList{frame}{cellNum}.relint1)),cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                        relintarray1(maxsizel2a-k:maxsizel2a+k,passed)=interpint1;
                        interpint2 = interp1(temp(1:length(cellList{frame}{cellNum}.relint2)),cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                        relintarray2(maxsizel2a-k:maxsizel2a+k,passed)=interpint2;%#ok
                    end
                   
                    cellLength=[cellLength cellList{frame}{cellNum}.length];%#ok
                end
            end

            % % cells length array is concatonated with the fluorescence matrix. This matrix is then sorted by length in ascending order
            numlist = [1:1:maxCellNum];%#ok
            lvint0=cat(2,numlist',cellLength');
            if signal(1) ==1 && sum(signal) == 1
                lvint1=cat(2,lvint0,relintarray1');
                lnumsort1=sortrows(lvint1,[2]);%#ok
            elseif signal(2) == 1 && sum(signal) == 1
                lvint2=cat(2,lvint0,relintarray2');
                lnumsort2=sortrows(lvint2,[2]);%#ok
            elseif sum(signal) == 2
                lvint1=cat(2,lvint0,relintarray1');
                lvint2=cat(2,lvint0,relintarray2');
                lnumsort1=sortrows(lvint1,[2]);%#ok
                lnumsort2=sortrows(lvint2,[2]);%#ok
            end
           
            if signal(1) ==1 && sum(signal) == 1
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            elseif signal(2) == 1 && sum(signal) == 1
               %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            elseif sum(signal) == 2
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
               
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            end
    case 'randomN'
     
            replacement=false;

            %%finds the maximum number of stepareas inside of a cell from the cellList
            maxsizelarray=[];
            n=0;
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    if isempty(cellList{frame}{cellNum}) || length(cellList{frame}{cellNum}.mesh)<4 ...
                            ||~isfield(cellList{frame}{cellNum},signal1) ...
                            || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                        n=n+1;
                end
            end
            if n<=maxCellNum
                maxCellNum=n;
            end
            rand=randsample(n,maxCellNum,replacement);
            n=0;
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    if isempty(cellList{frame}{cellNum}) ...
                            || length(cellList{frame}{cellNum}.mesh)<4 ...
                            || ~isfield(cellList{frame}{cellNum},signal1) ...
                            || isempty(cellList{frame}{cellNum}.(signal1)) ...
                            || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                    n = n+1;
                    b=rand==n;
                   
                    if sum(b)~=1
                        continue
                    end
                        maxsizelarray=[maxsizelarray length(cellList{frame}{cellNum}.lengthvector)];%#ok<AGROW>
                end
            end

            %using the maxima from above, a matrix consiting of zeros is created to be
            %filled in by mesh intensities
            relintarray1=zeros(maxCellLength,maxCellNum);
            maxsizel=maxCellLength;
            % if maxCellLength > maxsizel;
            %     maxsizel=maxCellLength;
            % end
            maxsizel2 = ceil(maxsizel); if mod(maxsizel2,2)==0, maxsizel2=maxsizel2+1; end
            maxsizel2a = maxsizel2/2+0.5;
            n=0;
            passed=0;
            cellLength=[];

            %zeroarray is replaced with relative segment intensity data from the cell
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    place=1;%#ok
                    if isempty(cellList{frame}{cellNum}) || length(cellList{frame}{cellNum}.mesh)<4 ...
                        ||~isfield(cellList{frame}{cellNum},signal1) ...
                        || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                    n = n+1;
                    b=rand==n;
                   
                    if sum(b)~=1
                        continue
                    end
                    passed=passed+1;
                    if signal(1) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                    elseif signal(2) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    elseif sum(signal) == 2
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                       
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    else
                        disp('provide information in signal variable')
                        return;
                    end
                       
                    %%% A MOVING AVERAGE IS CALCULATED FOR EACH OF THE SEGMENTS TO FIND THE SINGLE BRIGHTEST SEGMENT AREA
                    cellList{frame}{cellNum}.meshavg=[];
                    if signal(1) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint1)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint1(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1; %#ok
                        end
                    elseif signal(2) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint2)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint2(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1; %#ok
                        end
                    end
                    %%WITH THE BRIGHTEST SEGMENT CALCULATED ABOVE WE CAN ORIENT THE
                    %%CELL SO THAT THE BRIGHTEST SEGMENT IS ON THE RIGHTS (i.e. WITH FtsZ BEING POLAR ON RIGHT(NEW POLE)
                    %%AND LARGER STALK CELL BIAS LETTING THE FtsZ RING BE ON THE RIGHT
                    %%AS WELL)
                    [~,maxavg]=max(cellList{frame}{cellNum}.meshavg);
                    if maxavg<=length(cellList{frame}{cellNum}.meshavg)/2+1; %#ok
                        if signal(1) ==1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                        elseif signal(2) == 1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        elseif sum(signal) == 2
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        end
                           
                    end
                    k = floor(cellList{frame}{cellNum}.length/2);
                    temp = cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2;
                    if signal(1) ==1 && sum(signal) == 1
                        interpint1 = interp1(temp(1:length(cellList{frame}{cellNum}.relint1)),cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                        relintarray1(maxsizel2a-k:maxsizel2a+k,passed)=interpint1;
                    elseif signal(2) == 1 && sum(signal) == 1
                        interpint2 = interp1(temp(1:length(cellList{frame}{cellNum}.relint2)),cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                        relintarray2(maxsizel2a-k:maxsizel2a+k,passed)=interpint2; %#ok<AGROW>
                    elseif sum(signal) == 2
                        interpint1 = interp1(temp(1:length(cellList{frame}{cellNum}.relint1)),cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                        relintarray1(maxsizel2a-k:maxsizel2a+k,passed)=interpint1;
                        interpint2 = interp1(temp(1:length(cellList{frame}{cellNum}.relint2)),cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                        relintarray2(maxsizel2a-k:maxsizel2a+k,passed)=interpint2; %#ok<AGROW>
                    end
                   
                    cellLength=[cellLength cellList{frame}{cellNum}.length]; %#ok<AGROW>
                end
            end

            % % cells length array is concatonated with the fluorescence matrix. This matrix is then sorted by length in ascending order
            numlist=[1:1:maxCellNum]; %#ok
            lvint0=cat(2,numlist',cellLength');
            if signal(1) ==1 && sum(signal) == 1
                lvint1=cat(2,lvint0,relintarray1');
                lnumsort1=sortrows(lvint1,[2]); %#ok
            elseif signal(2) == 1 && sum(signal) == 1
                lvint2=cat(2,lvint0,relintarray2');
                lnumsort2=sortrows(lvint2,[2]); %#ok
            elseif sum(signal) == 2
                lvint1=cat(2,lvint0,relintarray1');
                lvint2=cat(2,lvint0,relintarray2');
                lnumsort1=sortrows(lvint1,[2]); %#ok
                lnumsort2=sortrows(lvint2,[2]); %#ok
            end
           
            if signal(1) ==1 && sum(signal) == 1
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            elseif signal(2) == 1 && sum(signal) == 1
               %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            elseif sum(signal) == 2
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
               
                %relative intensities are plotted accoring to a colormap
                figure
                x=[.0643*-maxCellLength./2 .0643*maxCellLength./2];
                y=[1 length(maxsizel+2)];
                xlabel('Distance From Midcell','FontSize',18)
                ylabel('Number of Cells','FontSize',18)
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            end
    case 'constriction'
            DC = [];
            replacement=false; %#ok
            %%finds the maximum number of stepareas inside of a cell from the cellList
            sizel=[];
            sizelarray=[];
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                        sizelarray= [sizelarray sizel];%#ok<AGROW>
                end
            end
            %using the maxima from above, a (max X cell number)matrix consiting of zeros is created
            % rand=randsample(length(sizelarray),maxCellNum,replacement);
            maxCellLength2 = ceil(maxCellLength); if mod(maxCellLength2,2)==0, maxCellLength2=maxCellLength2+1; end
            maxCellLength2a = maxCellLength2/2+0.5;
            relintarray1=zeros(maxCellLength2,length(sizelarray));
            lengthvectorarray=zeros(maxCellLength,length(sizelarray));
            n=0;
            cellLength=[];
            cellArea=[];
            %zeroarray is replaced with relative segment intensity data from the cell
            for frame = frameList
                for cellNum = 1:length(cellList{frame})
                    place=1; %#ok
                    if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                        continue
                    end
                    n = n+1;
                    if signal(1) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                    elseif signal(2) == 1 && sum(signal) == 1
                        %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    elseif sum(signal) == 2
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                       
                         %%calculates the fluorescent intensities in each segment normalized
                        %%by the area of that segment
                        if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                        else
                            cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                        end
                        %%segments are then normalized to the brightest segment so that
                        %%this sigment is represented as 1.
                        cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                    else
                        disp('provide information in signal variable')
                        return;
                    end
                    cellList{frame}{cellNum}.meshavg=[];
                    if signal(1) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint1)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint1(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1; %#ok
                        end
                    elseif signal(2) == 1 && sum(signal) == 1
                        for place = 1:(length(cellList{frame}{cellNum}.relint2)-(numPixelsMovingAverage-1));
                            cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint2(place:(place+(numPixelsMovingAverage-1))))];
                            place=place+1; %#ok
                        end
                    end

            %         [qwert,maxavg]=max(cellList{f}{c}.meshavg);
            %         if maxavg<=length(cellList{f}{c}.meshavg)/2+1;
            %             %             cellList{f}{c}.relint2=flipud(cellList{f}{c}.relint2);
            %             cellList{f}{c}.relint1=flipud(cellList{f}{c}.relint1);
            %         end
                    lngvector = cellList{frame}{cellNum}.lengthvector;
                    lng = cellList{frame}{cellNum}.length;
                    k = floor(lng/2);
                    if signal(1) == 1 && sum(signal) == 1
                        relint1=cellList{frame}{cellNum}.relint1;
                        interpint1 = interp1(lngvector-lng/2,relint1,-k:k,'linear','extrap');
                        relintarray1(maxCellLength2a-k:maxCellLength2a+k,n)=interpint1;
                        ind1 = length(relint1);
                    elseif signal(2) == 1 && sum(signal) == 1                       
                        relint2=cellList{frame}{cellNum}.relint2;
                        interpint2 = interp1(lngvector-lng/2,relint2,-k:k,'linear','extrap');
                        relintarray2(maxCellLength2a-k:maxCellLength2a+k,n)=interpint2; %#ok<AGROW>
                        ind1 = length(relint2);
                    elseif sum(signal) == 2
                        relint1=cellList{frame}{cellNum}.relint1;
                        interpint1 = interp1(lngvector-lng/2,relint1,-k:k,'linear','extrap');
                        relintarray1(maxCellLength2a-k:maxCellLength2a+k,n)=interpint1;
                        relint2=cellList{frame}{cellNum}.relint2;
                        interpint2 = interp1(lngvector-lng/2,relint2,-k:k,'linear','extrap');
                        relintarray2(maxCellLength2a-k:maxCellLength2a+k,n)=interpint2;%#ok<AGROW>
                    end
                   
                    ind2 = round(maxCellLength/2-ind1/2);
                    lengthvectorarray(ind2+1:ind2+ind1,n)=cellList{frame}{cellNum}.lengthvector;
                    cellL=cellList{frame}{cellNum}.length;
                    cellLength=[cellLength cellL];%#ok<AGROW>
                    cellArea=[cellArea cellList{frame}{cellNum}.area];%#ok<AGROW> 

                    prf = cellList{frame}{cellNum}.signal0;
                            if isempty(prf),break; end
                            for i=1:2
                                prf = 0.5*prf + 0.25*(prf([1 1:end-1])+prf([2:end end]));
                            end
                            minima = [false reshape((prf(2:end-1)<prf(1:end-2))&(prf(2:end-1)<prf(3:end)),1,[]) false];
                            if isempty(minima) || sum(prf)==0
                                minsize=0;
                                ctpos = []; %#ok
                            else
                                im = find(minima);
                                minsize = 0; %#ok
                                ctpos = 0; %#ok
                                dh = [];
                                dhi = [];
                                hgt = [];
                                for k=1:length(im)
                                    i=im(k);
                                    half1 = prf(1:i-1);
                                    half2 = prf(i+1:end);
                                    dh1 = max(half1)-prf(i);
                                    dh2 = max(half2)-prf(i);
                                    dh(k) = min(dh1,dh2); %#ok
                                    dhi(k) = mean([dh1 dh2]); %#ok
                                    hgt(k) = prf(i)+dhi(k); %#ok
                                end
                                [~,i] = max(dh);
                                minsizeabs = dhi(i);
                                minsize = minsizeabs/hgt(i);
                                ctpos = im(i); %#ok
                                if isempty(minsize), minsize=0; end
                            end
                        DC = [DC minsize]; %#ok<AGROW>
                end
            end
            % % % cells are sorted by length in ascending order
            numlist=[1:1:n]; %#ok
            reverseNumlist=[n:-1:1]; %#ok
            plotsizel=2+maxCellLength;
            if signal(1) ==1 && sum(signal) == 1
                lvint1=cat(2,cellLength',DC');
                lvint1=cat(2,lvint1,relintarray1');
                lnumsort1=sortrows(lvint1,[1]); %#ok
            elseif signal(2) == 1 && sum(signal) == 1
                lvint2=cat(2,cellLength',DC');
                lvint2=cat(2,lvint2,relintarray2');
                lnumsort2=sortrows(lvint2,[1]); %#ok
            elseif sum(signal) == 2
                lvint1=cat(2,cellLength',DC');
                lvint1=cat(2,lvint1,relintarray1');
                lnumsort1=sortrows(lvint1,[1]); %#ok

                lvint2=cat(2,cellLength',DC');
                lvint2=cat(2,lvint2,relintarray2');
                lnumsort2=sortrows(lvint2,[1]); %#ok
            end
   
            % stepwidth=cat(2,cellLength',numlist');
            % stepwidth=cat(2,stepwidth,relwidtharray');
            % lwidthsort=sortrows(stepwidth,[1]);
            % plotsizel=2+maxCellLength;
% % %             lengthvectorarray=cat(2,cellLength',lengthvectorarray');
% % %             lengthvectorarray=sortrows(lengthvectorarray,[1]);
% % %             length_int=cat(2,lengthvectorarray,lnumsort1);
% % %             sizelengthint=length(length_int);
% % %             %relative intensities are plotted accoring to a colormap
% % %             rand=randsample(length(lnumsort1),maxCellNum,replacement);
            if signal(1) ==1 && sum(signal) == 1
                figure
                x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
                y=[1 length(plotsizel)];
                image(x,y,64*lnumsort1(1:end,3:plotsizel));colormap jet;
                figure,scatter(lnumsort1(1:end,2),reverseNumlist)
            elseif signal(2) == 1 && sum(signal) == 1               
                figure
                x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
                y=[1 length(plotsizel)];
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
                figure,scatter(lnumsort2(1:end,2),reverseNumlist)
            elseif sum(signal) == 2
                %signal 1
                figure
                x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
                y=[1 length(plotsizel)];
                image(x,y,64*lnumsort1(1:end,3:plotsizel));colormap jet;
                figure,scatter(lnumsort1(1:end,2),reverseNumlist)
                %signal 2
                figure
                x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
                y=[1 length(plotsizel)];
                image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
                figure,scatter(lnumsort2(1:end,2),reverseNumlist)
            end
           
    case 'sort_by_constriction'
        %%finds the maximum number of stepareas inside of a cell from the cellList
        maxsizelarray=[];
        sizel=[];
        sizelarray=[];
        for frame = frameList
            for cellNum = 1:length(cellList{frame})
                if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                    continue
                end
            maxsizelarray=[maxsizelarray length(cellList{frame}{cellNum}.lengthvector)]; %#ok<AGROW>
            sizelarray= [sizelarray sizel]; %#ok<AGROW>
            end
        end
        %using the maxima from above, a matrix consiting of zeros is created to be
        %filled in by mesh intensities
        maxsizel=max(maxsizelarray);
        maxsizel2 = ceil(maxsizel); if mod(maxsizel2,2)==0, maxsizel2=maxsizel2+1; end
        maxsizel2a = maxsizel2/2+0.5;
        relintarray1 = [ ];
        relintarray2 = [ ];
        lengthvectorarray=zeros(maxsizel,length(sizelarray)); %#ok
        n=0;
        cellLength=[];
        DC=[];
        %zeroarray is replaced with relative segment intensity data from the cell
        for frame = frameList
            for cellNum = 1:length(cellList{frame})
                place=1; %#ok
                if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                    continue
                end

                n = n+1;
                %%calculates the fluorescent intensities in each segment normalized
                %%by the area of that segment
                if signal(1) == 1 && sum(signal) == 1
                    if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                elseif signal(2) == 1 && sum(signal) == 1
                    %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                elseif sum(signal) == 2
                     %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));

                     %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                else
                    disp('provide information in signal variable')
                    return;
                end

                %%% A MOVING AVERAGE IS CALCULATED FOR EACH OF THE SEGMENTS TO FIND THE SINGLE BRIGHTEST SEGMENT AREA
                cellList{frame}{cellNum}.meshavg=[];
                if signal(1) == 1 && sum(signal) == 1
                    for place = 1:(length(cellList{frame}{cellNum}.relint1)-(numPixelsMovingAverage-1));
                        cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint1(place:(place+(numPixelsMovingAverage-1))))];
                        place=place+1; %#ok
                    end
                elseif signal(2) == 1 && sum(signal) == 1
                     for place = 1:(length(cellList{frame}{cellNum}.relint2)-(numPixelsMovingAverage-1));
                        cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint2(place:(place+(numPixelsMovingAverage-1))))];
                        place=place+1;%#ok
                     end
                end
                %%WITH THE BRIGHTEST SEGMENT CALCULATED ABOVE WE CAN ORIENT THE
                %%CELL SO THAT THE BRIGHTEST SEGMENT IS ON THE RIGHTS (i.e. WITH FtsZ BEING POLAR ON RIGHT(NEW POLE)
                %%AND LARGER STALK CELL BIAS LETTING THE FtsZ RING BE ON THE RIGHT
                %%AS WELL)
                [~,maxavg]=max(cellList{frame}{cellNum}.meshavg);
                if maxavg<=length(cellList{frame}{cellNum}.meshavg)/2+1;
                        if signal(1) ==1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                        elseif signal(2) == 1 && sum(signal) == 1
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        elseif sum(signal) == 2
                            cellList{frame}{cellNum}.relint1=flipud(cellList{frame}{cellNum}.relint1);
                            cellList{frame}{cellNum}.relint2=flipud(cellList{frame}{cellNum}.relint2);
                        end
                end

                k = floor(cellList{frame}{cellNum}.length/2);
                v=(1/(2*maxsizel)):(1/maxsizel):1;%#ok
                if signal(1) == 1 && sum(signal) == 1
                    interpint1 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                    relintarray1(maxsizel2a-k:maxsizel2a+k,n)=interpint1;%#ok<AGROW>
                elseif signal(2) == 1 && sum(signal) == 1
                    interpint2 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                    relintarray2(maxsizel2a-k:maxsizel2a+k,n)=interpint2; %#ok<AGROW>
                elseif  sum(signal) == 2
                    interpint1 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                    relintarray1(maxsizel2a-k:maxsizel2a+k,n)=interpint1;%#ok<AGROW>
                   
                    interpint2 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                    relintarray2(maxsizel2a-k:maxsizel2a+k,n)=interpint2; %#ok<AGROW>
                end
                cellLength=[cellLength cellList{frame}{cellNum}.length]; %#ok<AGROW>

                prf = cellList{frame}{cellNum}.signal0;
                if isempty(prf),  break; end
                for i=1:2
                    prf = 0.5*prf + 0.25*(prf([1 1:end-1])+prf([2:end end]));
                end
                minima = [false reshape((prf(2:end-1)<prf(1:end-2))&(prf(2:end-1)<prf(3:end)),1,[]) false];
                if isempty(minima) || sum(prf)==0
                    minsize=0;
                    ctpos = [];%#ok
                else
                    im = find(minima);
                    minsize = 0;%#ok
                    ctpos = 0;%#ok
                    dh = [];
                    dhi = [];
                    hgt = [];
                    for k=1:length(im)
                        i=im(k);
                        half1 = prf(1:i-1);
                        half2 = prf(i+1:end);
                        dh1 = max(half1)-prf(i);
                        dh2 = max(half2)-prf(i);
                        dh(k) = min(dh1,dh2);%#ok
                        dhi(k) = mean([dh1 dh2]);%#ok
                        hgt(k) = prf(i)+dhi(k);%#ok
                    end
                    [~,i] = max(dh);
                    minsizeabs = dhi(i);
                    minsize = minsizeabs/hgt(i);
                    ctpos = im(i);%#ok
                    if isempty(minsize), minsize=0; end
                end
                DC = [DC minsize]; %#ok<AGROW>
            end
        end

        % % cells length array is concatonated with the fluorescence matrix. This matrix is then sorted by length in ascending order
        numlist=[1:1:n]; %#ok
        reverseNumlist=[n:-1:1];%#ok

        if signal(1) ==1 && sum(signal) == 1
            lvint1=cat(2,cellLength',DC');
            lvint1=cat(2,lvint1,relintarray1');
            lnumsort1=sortrows(lvint1,[2]); %#ok
        elseif signal(2) == 1 && sum(signal) == 1
            lvint2=cat(2,cellLength',DC');
            lvint2=cat(2,lvint2,relintarray2');
            lnumsort2=sortrows(lvint2,[2]); %#ok
        elseif sum(signal) == 2
            lvint1=cat(2,cellLength',DC');
            lvint1=cat(2,lvint1,relintarray1');
            lnumsort1=sortrows(lvint1,[2]); %#ok

            lvint2=cat(2,cellLength',DC');
            lvint2=cat(2,lvint2,relintarray2');
            lnumsort2=sortrows(lvint2,[2]); %#ok
        end
        %relative intensities are plotted accoring to a colormap
        if signal(1) ==1 && sum(signal) == 1
            figure
            x=[-.5 .5];
            y=[1 length(maxsizel+2)];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            figure,scatter(lnumsort1(1:end,2),reverseNumlist)
        elseif signal(2) == 1 && sum(signal) == 1
            figure
            x=[-.5 .5];
            y=[1 length(maxsizel+2)];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            figure,scatter(lnumsort2(1:end,2),reverseNumlist)
        elseif sum(signal) == 2
            figure
            x=[-.5 .5];
            y=[1 length(maxsizel+2)];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            figure,scatter(lnumsort1(1:end,2),reverseNumlist)
           
            figure
            x=[-.5 .5];
            y=[1 length(maxsizel+2)];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            figure,scatter(lnumsort2(1:end,2),reverseNumlist)
        end
    case 'constriction_noNormalization'
        %%finds the maximum number of stepareas inside of a cell from the cellList
        maxsizelarray=[];
        sizel=[];
        sizelarray=[];
         for frame = frameList
            for cellNum = 1:length(cellList{frame})
                if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                    continue
                end
            maxsizelarray=[maxsizelarray length(cellList{frame}{cellNum}.lengthvector)]; %#ok<AGROW>
            sizelarray= [sizelarray sizel]; %#ok<AGROW>
            end
        end
        %using the maxima from above, a matrix consiting of zeros is created to be
        %filled in by mesh intensities

        maxsizel=max(maxsizelarray);
        plotsizel=length(maxsizelarray);
        maxsizel2 = ceil(maxsizel); if mod(maxsizel2,2)==0, maxsizel2=maxsizel2+1; end
        maxsizel2a = maxsizel2/2+0.5;
        relintarray1=zeros(maxsizel2,length(sizelarray));
        lengthvectorarray=zeros(maxsizel,length(sizelarray));%#ok
        n=0;
        cellLength=[];
        lengthvectorarray=[];%#ok
        DC=[];
        %zeroarray is replaced with relative segment intensity data from the cell
        for frame = frameList
            for cellNum = 1:length(cellList{frame})
                place=1; %#ok
                if ~isfield(cellList{frame}{cellNum},signal1) || isempty(cellList{frame}{cellNum}.(signal1)) || cellList{frame}{cellNum}.length>maxCellLength
                    continue
                end
                n = n+1;
                %%calculates the fluorescent intensities in each segment normalized
                %%by the area of that segment
                if signal(1) == 1 && sum(signal) == 1
                    if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));
                elseif signal(2) == 1 && sum(signal) == 1
                    %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                elseif sum(signal) == 2
                     %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal1)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal1 = (cellList{frame}{cellNum}.(signal1)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal1))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint1 = (cellList{frame}{cellNum}.relsignal1./max(cellList{frame}{cellNum}.relsignal1));

                     %%calculates the fluorescent intensities in each segment normalized
                    %%by the area of that segment
                    if length(cellList{frame}{cellNum}.(signal2)) > length(cellList{frame}{cellNum}.steparea)
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)(1:length(cellList{frame}{cellNum}.steparea))./cellList{frame}{cellNum}.steparea);
                    else
                        cellList{frame}{cellNum}.relsignal2 = (cellList{frame}{cellNum}.(signal2)./cellList{frame}{cellNum}.steparea(1:length(cellList{frame}{cellNum}.(signal2))));
                    end
                    %%segments are then normalized to the brightest segment so that
                    %%this sigment is represented as 1.
                    cellList{frame}{cellNum}.relint2 = (cellList{frame}{cellNum}.relsignal2./max(cellList{frame}{cellNum}.relsignal2));
                else
                    disp('provide information in signal variable')
                    return;
                end

                %%% A MOVING AVERAGE IS CALCULATED FOR EACH OF THE SEGMENTS TO FIND THE SINGLE BRIGHTEST SEGMENT AREA
                cellList{frame}{cellNum}.meshavg=[];
                if signal(1) == 1 && sum(signal) == 1
                    for place = 1:(length(cellList{frame}{cellNum}.relint1)-(numPixelsMovingAverage-1));
                        cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint1(place:(place+(numPixelsMovingAverage-1))))];
                        place=place+1;%#ok
                    end
                elseif signal(2) == 1 && sum(signal) == 1
                    for place = 1:(length(cellList{frame}{cellNum}.relint2)-(numPixelsMovingAverage-1));
                        cellList{frame}{cellNum}.meshavg=[cellList{frame}{cellNum}.meshavg mean(cellList{frame}{cellNum}.relint2(place:(place+(numPixelsMovingAverage-1))))];
                        place=place+1;%#ok
                    end
                end
                %%WITH THE BRIGHTEST SEGMENT CALCULATED ABOVE WE CAN ORIENT THE
                %%CELL SO THAT THE BRIGHTEST SEGMENT IS ON THE RIGHTS (i.e. WITH FtsZ BEING POLAR ON RIGHT(NEW POLE)
                %%AND LARGER STALK CELL BIAS LETTING THE FtsZ RING BE ON THE RIGHT
                %%AS WELL)
                k = floor(cellList{frame}{cellNum}.length/2);
                v=1/(2*maxsizel):1/maxsizel:1;%#ok
                if signal(1) == 1 && sum(signal) == 1
                    interpint1 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                    relintarray1(maxsizel2a-k:maxsizel2a+k,n)=interpint1;
                elseif signal(2) == 1 && sum(signal) == 1
                    interpint2 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                    relintarray2(maxsizel2a-k:maxsizel2a+k,n)=interpint2; %#ok<AGROW>
                elseif  sum(signal) == 2
                    interpint1 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint1,-k:k,'linear','extrap');
                    relintarray1(maxsizel2a-k:maxsizel2a+k,n)=interpint1;
                   
                    interpint2 = interp1(cellList{frame}{cellNum}.lengthvector-cellList{frame}{cellNum}.length/2,cellList{frame}{cellNum}.relint2,-k:k,'linear','extrap');
                    relintarray2(maxsizel2a-k:maxsizel2a+k,n)=interpint2; %#ok<AGROW>
                end
                cellLength=[cellLength cellList{frame}{cellNum}.length]; %#ok<AGROW>

                        prf = cellList{frame}{cellNum}.signal0;
                        if isempty(prf),break; end
                        for i=1:2
                            prf = 0.5*prf + 0.25*(prf([1 1:end-1])+prf([2:end end]));
                        end
                        minima = [false reshape((prf(2:end-1)<prf(1:end-2))&(prf(2:end-1)<prf(3:end)),1,[]) false];
                        if isempty(minima) || sum(prf)==0
                            minsize=0;
                            ctpos = [];%#ok
                        else
                            im = find(minima);
                            minsize = 0;%#ok
                            ctpos = 0;%#ok
                            dh = [];
                            dhi = [];
                            hgt = [];
                            for k=1:length(im)
                                i=im(k);
                                half1 = prf(1:i-1);
                                half2 = prf(i+1:end);
                                dh1 = max(half1)-prf(i);
                                dh2 = max(half2)-prf(i);
                                dh(k) = min(dh1,dh2);%#ok
                                dhi(k) = mean([dh1 dh2]);%#ok
                                hgt(k) = prf(i)+dhi(k);%#ok
                            end
                            [~,i] = max(dh);
                            minsizeabs = dhi(i);
                            minsize = minsizeabs/hgt(i);
                            ctpos = im(i);%#ok
                            if isempty(minsize), minsize=0; end
                        end
                    DC = [DC minsize];%#ok<AGROW>
            end
        end

        % % cells length array is concatonated with the fluorescence matrix. This matrix is then sorted by length in ascending order
        numlist=[1:1:n];%#ok
        reverseNumlist=[n:-1:1];%#ok
        if signal(1) == 1 && sum(signal) == 1
            lvint1=cat(2,cellLength',DC');
            lvint1=cat(2,lvint1,relintarray1');
            lnumsort1=sortrows(lvint1,[2]);%#ok
        elseif signal(2) == 1 && sum(signal) == 1
            lvint2=cat(2,cellLength',DC');
            lvint2=cat(2,lvint2,relintarray2');
            lnumsort2=sortrows(lvint2,[2]);%#ok
        elseif sum(signal) == 2
            lvint1=cat(2,cellLength',DC');
            lvint1=cat(2,lvint1,relintarray1');
            lnumsort1=sortrows(lvint1,[2]);%#ok
           
            lvint2=cat(2,cellLength',DC');
            lvint2=cat(2,lvint2,relintarray2');
            lnumsort2=sortrows(lvint2,[2]);%#ok
        end

        if signal(1) == 1 && sum(signal) == 1
            %relative intensities are plotted accoring to a colormap
            figure
            x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
            y=[plotsizel 1];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            % x2=[0 max(DC)];
            figure,scatter(lnumsort1(1:end,2),reverseNumlist)
        elseif signal(2) == 1 && sum(signal) == 1
            %relative intensities are plotted accoring to a colormap
            figure
            x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
            y=[plotsizel 1];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            % x2=[0 max(DC)];
            figure,scatter(lnumsort2(1:end,2),reverseNumlist)
           
        elseif sum(signal) == 2
            %relative intensities are plotted accoring to a colormap
            figure
            x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
            y=[plotsizel 1];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort1(1:end,3:end));colormap jet;
            % x2=[0 max(DC)];
            figure,scatter(lnumsort1(1:end,2),reverseNumlist)
           
            %relative intensities are plotted accoring to a colormap
            figure
            x=[-.0643*max(cellLength)./2 .0643*max(cellLength)./2];
            y=[plotsizel 1];
            xlabel('Distance From Midcell','FontSize',18)
            ylabel('Number of Cells','FontSize',18)
            image(x,y,64*lnumsort2(1:end,3:end));colormap jet;
            % x2=[0 max(DC)];
            figure,scatter(lnumsort2(1:end,2),reverseNumlist)
           
        end   
       
    otherwise
        disp('descriptor variable must contain one of the following values')
        descriptorValues %#ok
end

end


