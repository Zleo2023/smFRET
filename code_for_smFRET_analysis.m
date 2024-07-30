

clear all
close all

directory = 'D:\xxxxx\';
filename = 'xxxxx.tif';

infilename = strcat(directory,filename);

% channel = 'Cy5';
% channel = 'Cy3';
channel = 'whole';

% overwrite = 'yes';
overwrite = 'yes';
% 
% stdfactor = 2; % Threshold for particle finding
stdfactor = 3;

%Frames to use in averaging
startframe = 1;
endframe = 30;

plotcolor = 'r';

selectionbyflucts = 0;

analyze_all = 'y'; %Analyze all traces without filtering based on model results
% analyze_all = 'n';
analysis = 'manual';
% analysis = 'auto';

% mincounts = 1000; % For Cy5 in 5/25/2012
% mincounts = 1000; % For Cy3 in 5/25/2012
% mincounts = 200; % For expt 5/22/2012
mincounts = 5000;

driftcorr = 'no';
% driftcorr = 'no';
driftbins = 5;

maxL = 500; % Maximum intensity value allowed for the "low" (no molecule) state


outfilename = strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_avg_bgsub.tif');
tracesname = strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_traces.dat');
coordsname = strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_coords.dat');

% Check for files containing the molecule numbers, coordinates, and
% intervals for fitting.  If they exist; add to them.  If not, start at
% beginning.
if exist(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molintervals.mat'))==2 && strcmp(overwrite,'yes')~=1
    molnum = load(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molnum.mat'));
    molnum = molnum.molnum;
    molcoords = load(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molcoords.mat'));
    molcoords = molcoords.molcoords;
    molintervals = load(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molintervals.mat'));
    molintervals = molintervals.molintervals;
    mm = size(molnum,1)+1;
elseif exist(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molnum.mat'))==2
    molnum = load(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molnum.mat'));
    molnum = molnum.molnum;
    molcoords = load(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molcoords.mat'));
    molcoords = molcoords.molcoords;
    molintervals = cell(1,0);
    overwrite = 'no';
    mm=1;
else
    molnum = zeros(0,1);
    molcoords = zeros(0,2);
    molintervals = cell(1,0);
    mm = 1;
end


binsize = endframe-startframe+1;
% binsize = 300; % Number of frames to include in average


% Read TIF info

mov1info = imfinfo(infilename,'tif');

frames = size(mov1info,1);
xdim = mov1info(1,1).Width;
ydim = mov1info(1,1).Height;

if endframe > frames
    endframe = frames;
end

if endframe <= frames && startframe >= 1
    
    if strcmp(driftcorr,'yes')==1
        [driftxy]=drift_corr_tif(infilename,1,driftbins,5,5);
        for dd = 1:size(driftxy,1)
            for ff = 1:size(driftxy,2)
                driftorigin = driftxy{dd,ff}(round((endframe+startframe)/2),1:2);
                for ii = 1:size(driftxy{dd,ff})
                    driftxy{dd,ff}(ii,:)=driftxy{dd,ff}(ii,:)-driftorigin;
                end
            end
        end
    else
        driftxy{1,1}=zeros(frames,2);
    end
    
    
    
    if exist(outfilename) ~= 2 || strcmp(overwrite,'yes') ==1  % Only create averaged image if doesn't already exist
        
        tic;
        imagesum = zeros(ydim,xdim);
         if selectionbyflucts == 1
            fmap = zeros(ydim,xdim);
            for frame = startframe:endframe
                adata = imread(infilename,'tif','Index',frame);
                adata = im2double(adata);
                
                if frame > startframe
                    adiff = abs(adata-lastframe);
                    adiff = adiff - imopen(adiff,strel('disk',15)); % This step is one speed bottleneck
                    astd = std(reshape(adiff, size(adiff,1)*size(adiff,2),1));
                    adiff = imhmax(adiff,stdfactor*astd,4); % This step is another speed bottleneck

                    fmap = fmap + adiff;
                end
                
%                 imagesum = (imagesum + adata);
                lastframe = adata;
                %     disp(strcat('Working: frame ',num2str(frame)));
            end
            imageavg = fmap/binsize;
            imageavg = imageavg./max(max(imageavg));
            imshow(imageavg);
%             imageavg = im2uint16(imageavg);
        else
        for frame = startframe:endframe
            adata = imread(infilename,'tif','Index',frame);
            adata = im2double(adata);
            
            imagesum = (imagesum + adata);         
        end
        imageavg = imagesum/binsize;
        end
        imageavg = imageavg-imopen(imageavg,strel('disk',15));
        imageavg = imageavg-min(min(imageavg));
        
        imageavg2 = imageavg/max(max(imageavg));
        
        tfinish1 = toc;
        
        disp(strcat('Finished averaging',{' '},num2str(endframe-startframe+1),{' '},'frames in',{' '},num2str(tfinish1),{' '},'seconds.'));
        
        imshow(imageavg2);
        
        
        pause(0.01);
       
        
        % imageavg = im2uint16(imageavg);
        
        imwrite(im2uint8(imageavg2),outfilename,'tif');
        
        clear adata;
        
        clear imageavg2;
        
    else
        
              
        disp('Averaged image exists; skipping...');
        
        % Read averaged TIF
        
        adata = imread(outfilename, 'tif');

        imageavg = im2double(adata);
        
        imageavg2 = (imageavg-min(min(imageavg)));
        imageavg2 = imageavg2/max(max(imageavg2));
        
        imshow(imageavg2);
        clear imageavg2;
    end
    
    
    
    %*********************************************
    
    % Find maxima in image
    
    if strcmp(channel,'whole')==1
        aIm = imageavg;
        
    else
    
    if strcmp(channel,'Cy5')==1
%         Rect = [257 0 256 512];
          Rect = [xdim/2+1 0 xdim/2 ydim];
    elseif strcmp(channel,'Cy3')==1
%         Rect = [0 0 256 512];
          Rect = [0 0 xdim/2 ydim];
    else
          Rect = [0 0 xdim ydim];
    end
    
    aIm = imcrop(imageavg, Rect);
    end
    
    aImlinear = reshape(aIm,1,size(aIm,1)*size(aIm,2));
    aImstd = std(aImlinear);
    aImmedian = median(aImlinear);
    %     end
    %   Hard intensity cutoff -- ignore pixels below cutoff intensity
    Icutoff = aImmedian; %+ stdfactor*aImstd; % In my experience, histogramming the Cy5 channel intensities showed that the median is a reasonable cutoff.
    % aImcut = aIm-Icutoff;
    aImcut = aIm;
    
    aImcut2 = aImcut/max(max(aImcut));
    
    imshow(aImcut2);
    
    pause(0.01);
    
    hold on
    
    % clear aImcut2;
    
    for i = 1:size(aImcut,1)
        for j = 1:size(aImcut,2)
            if aImcut(i,j)<0
                aImcut(i,j)=0;
            end
        end
    end
    
    aImSupp = imhmax(aImcut,aImstd*stdfactor, 8); %Suppress local maxima
    % imshow(aImSupp);
    aImrMax = imregionalmax(aImSupp);
    % imshow(aImrMax);
    regions = bwlabel(aImrMax);
    %     figure(2)
    %     imshow(regions);
    centroids = regionprops(regions, 'centroid');
    RegArr = struct2cell(centroids);
    RegMat = cell2mat(RegArr.');
    
    if exist(coordsname) == 2 && strcmp(overwrite,'yes')~=1 % Record molecule coordinates if coordinates file doesn't already exist
        
        coordscell = load(coordsname,'-mat');
        X = round(coordscell.coords(:,1));
        Y = round(coordscell.coords(:,2));
        molecules = coordscell.coords;
        if min(molecules(:,1)) > xdim/2-1
%             molecules(:,1)=molecules(:,1)-256;
              molecules(:,1)=molecules(:,1)-xdim/2;
        end
        
        nummols=size(molecules,1);
        
    else
        
        nummols=size(RegMat,1);
        molecules = zeros(0,2);
        for n = 1:nummols
            if RegMat(n,1)>30 && RegMat(n,1)<size(aIm,2)-30 && RegMat(n,2)>30 && RegMat(n,2)<size(aIm,1)-30
                molecules=cat(1,molecules,RegMat(n,:));
                %         if n>190
                %             disp(RegMat(n,:));
                
                %         end
            end
        end
        % for m = 1:nummols
        %     molnum(m,1)=m;
        % end
        % molecules = cat(2,molnum,RegMat);
        nummols=size(molecules,1);
        
        % clear RegMat RegArr
        
    end
    
    plot(molecules(:,1),molecules(:,2), 'wo');
    
    hold off
    
    startframe=1;
    endframe=frames;
    
    % **********************************************
    
    % Create intensity traces of molecules
    
    if exist(tracesname) ~= 2 || exist(coordsname) ~= 2 || strcmp(overwrite,'yes')==1
        
        tic;
        
        X = zeros(nummols,1);
        Y = zeros(nummols,1);
        
        for m = 1:nummols
            
            if strcmp(channel,'Cy5')==1
                X(m) = molecules(m,1)+xdim/2;
            else
                X(m) = molecules(m,1);
            end
            Y(m) = molecules(m,2);
            
        end
        
             
        if endframe > frames
            endframe=frames;
        end
        
        
        traces=zeros(nummols,frames);
        
        for i = startframe:endframe

            frame = i;
            
            b = imread(infilename,'tif','Index',frame);
            
            bnew = b;
            
            for m = 1:nummols
                %         disp('Molecule:');
                %         disp(num2str(m));
                %         disp(num2str(X(m)));
                %         disp(num2str(Y(m)));
                Xr(m)=round(X(m)+driftxy{1,1}(i,1));
                Yr(m)=round(Y(m)+driftxy{1,1}(i,2));
                region = bnew(Yr(m)-10:Yr(m)+10,Xr(m)-10:Xr(m)+10);
                region_bg = double(region);
                region_bg(9:13,9:13) = nan(5,5);
                %         if i == startframe
                background = nanmedian(reshape(region_bg,1,size(region_bg,1)*size(region_bg,2)));
                %         end
                particle = double(region(9:13,9:13))-background;
                I = sum(sum(particle));
                % traceframe = cat(2,i,I);
                traces(m,i) = I;
            end
            %         figure(1)
            %         imshow(particle/max(max(particle)), 'InitialMagnification', 1000);
            %         pause(0.02)
            % imsum = imsum+bnew;
            if mod(i,10)==0
                disp('Done with frame:');
                disp(frame);
                % figure(1)
                % region = region/max(max(bnew));
                % imshow(region, 'InitialMagnification', 1000);
                % figure(2)
                % plot(trace(:,1),trace(:,2))
                % xlim([0 max(trace(:,1))]);
                % ylim([0 max(trace(:,2))]);
            end
        end
        % imsum = imsum/(endframe-startframe);
        % imshow(imsum);
        coords = cat(2,round(X),round(Y));
        tfinish2 = toc;
        disp(strcat('Finished averaging',{' '},num2str(endframe-startframe+1),{' '},'frames in',{' '},num2str(tfinish2),{' '},'seconds.'));;
        save(coordsname, 'coords');
        save(tracesname, 'traces');
    else
        coordscell = load(coordsname,'-mat');
        X = coordscell.coords(:,1);
        Y = coordscell.coords(:,2);
        tracescell = load(tracesname,'-mat');
        traces = tracescell.traces;
    end
    
    
    m = 1;
    % mm = 1;
    tprob_list = zeros(0,1);
    firstmol = 1;
    while m <= nummols
        if m < 1
            m = 1;
        end
        if size(molnum,1)>size(molintervals,2) % If the molecule list is longer than the list of intervals, go to the next molecule in the list that has no corresponding intervals
            mm = size(molintervals,2)+1;
            m = molnum(mm,1);
        end
        figure(1)
        imshow(aImcut2);
        hold on
        plot(molecules(:,1),molecules(:,2), 'wo');
        plot(molecules(m,1),molecules(m,2), 'ro');
        hold off
        %     figure(2)
        figure(2)
        %     [Ihist,Ibins]=hist(traces(m,:)-min(min(traces(m,:))),15);
        [Ihist,Ibins]=hist(traces(m,:),30);
        bar(Ibins,Ihist);
        

        figure(3)
        
%         if strcmp(channel,'Cy3')==1
            plot(traces(m,:), plotcolor, 'LineWidth', 2);
%         else
%             plot(traces(m,:), 'r', 'LineWidth', 2);
%         end
        xlim([startframe-1 endframe+1]);
        ylim([min(traces(m,:))-1000 max(traces(m,:))+1000]);
        xlabel('Frame');
        ylabel('Fluorescence (A.U.)');
        text(0.45, 0.9, sprintf('Candidate %d',m),'units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
        text(0.45, 0.8, strcat('(',num2str(X(m)),',',num2str(Y(m)),')'),'units','normalized','FontSize', 16, 'FontWeight', 'bold');
        %     hold on
        
        
        disp(sprintf('Candidate %d',m));
        disp(strcat('(',num2str(X(m)),',',num2str(Y(m)),')'));
        
        if strcmp(analysis,'manual')==1
            user1 = input('Command? n = next, p = previous, s = skip to trace, a = analyze (n)', 's');
        else
        user1 = 'a';
        end
        
        if strcmp(user1,'p') == 1
            m = m-1;
            hold off
        elseif strcmp(user1,'s') == 1
            m = input('What trace would you like to skip to?');
            if m < 1
                m = 1;
            elseif m > nummols
                m = nummols;
            end
        elseif strcmp(user1,'a') == 1
            
            FRET = mat2cell(traces(m,:));
%             FRET = mat2cell(traces(m,1:500));
            vbFRET_no_gui;
            hold on
            plot(x_hat{1,3});
            %         plot(abs(x_hat{1,3}-traces(m,:)'),'g');
            hold off
            figure(4)
            [pop,state]=hist(x_hat{1,3});
            bar(state,pop);
            [C,I] = max(bestOut{1,3}.m);
            Hind = I;
            [C2,I2] = min(bestOut{1,3}.m);
            Lind = I2;
            Mind = 6-I-I2;
            Htmat = bestOut{1,3}.Wa(Hind,:);
            Htmat(Hind)=0;                   % Remove the high-to-high transition from consideration
            [dummy,H_fastest] = max(Htmat);  % Find the fastest transition to/from the high state
            ycenter = 0;
            if size(molnum,1)>size(molintervals,2) || strcmp(analysis,'manual')==1  % If this is not the first set to be analyzed, analyze all molecules from molnum list
                    discard = 'n';
                    if bestOut{1,3}.m(H_fastest)<bestOut{1,3}.m(Hind)/3   % If the fastest transition from the high state is to a state less than 1/3 its own intensity, take the high state as the true single-molecule signature
                        ycenter = bestOut{1,3}.m(Hind);
                        tprob = sum(bestOut{1,3}.Wa(Hind,:))-bestOut{1,3}.Wa(Hind,Hind);
                        tprob = tprob/frames;
                    else                                                % Otherwise, take mid state as the true single-molecule signature
                        ycenter = bestOut{1,3}.m(Mind);
                        tprob = sum(bestOut{1,3}.Wa(Mind,:))-bestOut{1,3}.Wa(Mind,Mind);
                        tprob = tprob/frames;
                    end
            elseif C2 < maxL                        % Otherwise, require that the low state is below maxL counts
                if bestOut{1,3}.m(Hind) > mincounts
                    if bestOut{1,3}.m(H_fastest)<bestOut{1,3}.m(Hind)/3   % If the fastest transition from the high state is to a state less than 1/3 its own intensity, take the high state as the true single-molecule signature
                        ycenter = bestOut{1,3}.m(Hind);
                        tprob = sum(bestOut{1,3}.Wa(Hind,:))-bestOut{1,3}.Wa(Hind,Hind);
                        tprob = tprob/frames;
                        discard = 'n';
                    elseif bestOut{1,3}.m(Mind)>mincounts    % Otherwise, if the mid state is at least mincounts in intensity, take this as the true single-molecule signature
                        ycenter = bestOut{1,3}.m(Mind);
                        tprob = sum(bestOut{1,3}.Wa(Mind,:))-bestOut{1,3}.Wa(Mind,Mind);
                        tprob = tprob/frames;
                        discard = 'n';
                    else
                        discard = 'y';
                    end
                else
                    discard = 'y';              % Otherwise, discard this molecule
                    ycenter = 0;
                end
            else
                discard = 'y';
                ycenter = 0;
            end
%             ycenter = mincounts;
            %         for bb = 1:size(state,2)
            %             if state(1,bb)<max(bestOut{1,3}.m)/3 || state(1,bb) < 1500; %1500 for Cy5, 2000 for Cy3
            %                 pop(1,bb)=0;
            %             end
            %         end
            %         [C,I]=max(pop);
            %         [C2,I2]=min(abs(bestOut{1,3}.m-state(1,I)));
            %         ycenter = bestOut{1,3}.m(1,I2);
            ytop = 1000000;
            ybot = ycenter-0.25*ycenter;
            
            %         disp('Please select the intensity window for analysis.')
            
            disp(sprintf('Candidate %d',m));
            disp(strcat('(',num2str(X(m)),',',num2str(Y(m)),')'));
            disp(sprintf(' '));
            
            %         k = waitforbuttonpress;              % Hold program until user selects region
            %         point1 = get(gca,'CurrentPoint');    % button down detected
            %         finalRect = rbbox;                   % return figure units
            %         point2 = get(gca,'CurrentPoint');    % button up detected
            %         point1 = point1(1,1:2);             % extract x and y
            %         point2 = point2(1,1:2);
            %         Iwind = [min(point1(1,2),point2(1,2)) max(point1(1,2),point2(1,2))]+min(min(traces(m,:)));
            
            x = 1:1:frames;
            y1 = x.*0.+ybot;
            y2 = x.*0.+ytop;
            
            figure(3)
            hold on
            plot(x,y1, 'k', 'LineWidth', 2);
            plot(x,y2, 'k', 'LineWidth', 2);
            hold off
            
            %         discard
            if firstmol == 1
                    disp('Press any key to continue with analysis.');
                    firstmol = 0;
                    pause;
            end
            
            goodframes = zeros(0,1);
            for n = startframe:endframe
%                 if traces(m,n) > ybot && traces(m,n) < ytop
                if x_hat{1,3}(n) > ybot && x_hat{1,3}(n) < ytop
                    goodframes = cat(1,goodframes,n);
                end
            end
            if size(goodframes,1)~=0
                intervals = zeros(0,2);
                k=1;
                intervals(k,1)=goodframes(k,1);
                for p = 1:size(goodframes,1)-1;
                    if goodframes(p+1,1)==goodframes(p,1)+1
                    else
                        intervals(k,2)=goodframes(p,1);
                        k = k+1;
                        intervals(k,1)=goodframes(p+1,1);
                    end
                end
                intervals(k,2)=goodframes(size(goodframes,1),1);
                dwellt_on = zeros(size(intervals,1),1);
                dwellt_off = zeros(size(intervals,1)-1,1);
                for mol = 1:size(intervals,1)
                    dwellt_on(mol) = intervals(mol,2)-intervals(mol,1) + 1;
                    if mol > 1
                        dwellt_off(mol-1) = intervals(mol,1)-intervals(mol-1,2) - 1;
                    end
                end
                %     intervals
                %             disp(sprintf('[%d %d;', intervals(1,1), intervals(1,2)));
                %             for s = 2:size(intervals,1)-1
                %                 disp(sprintf('%d %d;', intervals(s,1),intervals(s,2)));
                %             end
                %             disp(sprintf('%d %d];', intervals(size(intervals,1),1), intervals(size(intervals,1),2)));
                %             Wanorm = bestOut{1,3}.Wa/max(max(bestOut{1,3}.Wa))
                %             Wa = bestOut{1,3}.Wa
                
                int_string =  sprintf('[%d %d;', intervals(1,1), intervals(1,2));
                for s = 2:size(intervals,1)-1
                    string1 = sprintf('%d %d;', intervals(s,1),intervals(s,2));
                    int_string = cat(2,int_string,string1);
                end
                int_string = cat(2,int_string,sprintf('%d %d];', intervals(size(intervals,1),1), intervals(size(intervals,1),2)));
                disp('intervals:');
                disp(num2str(size(intervals,1)))             
                disp(int_string);
%                 pause; 
%                 offdiagsum = sum(sum(bestOut{1,3}.Wa))-(bestOut{1,3}.Wa(1,1)+bestOut{1,3}.Wa(2,2)+bestOut{1,3}.Wa(3,3)); %Sum of off-diagonal transition probabilities

                %             offdiagsum = sum(sum(bestOut{1,2}.Wa))-(bestOut{1,2}.Wa(1,1)+bestOut{1,2}.Wa(2,2));
                %             keepyn = input('Keep this molecule? (y)','s');
                %             meddeviation = median(abs(x_hat{1,2}-traces(m,:)'))
                if strcmp(discard,'n')==1
                    tprob_list = cat(1,tprob_list,tprob);
                    if strcmp(analyze_all,'y')==1
                        keepyn='y';
                    elseif tprob < 0.025 && size(molnum,1)==size(molintervals,2) % Check that transitions are frequent
                        keepyn = 'n';
                        disp('Do not keep');
%                         pause(0.5);
                    else
                        keepyn = 'y';
                        disp('Keep');
                    end
                else
                    keepyn = 'n';
                    disp('Do not keep');
                end
                %             keepyn = 'n';
                if strcmp(keepyn,'n')==0
                    
                    if size(molnum,1)==size(molintervals,2)
                        molnum(mm,1) = m;
                        molcoords(mm,1:2)=[X(m) Y(m)];
                    end
                    molintervals{mm}=intervals;
                    mm = mm + 1;
                    figure(2)
                    pause(0.1);
                    set(gcf, 'PaperPositionMode', 'auto');%Make sure on-screen dimensions are preserved in output file
                    print('-f2', '-r600', '-djpeg', strcat(infilename(1,1:(size(infilename,2)-4)),'_cand_',num2str(m),'_',channel,'_hist.jpg'));
                    figure(3)
                    pause(0.1);
                    set(gcf, 'PaperPositionMode', 'auto');%Make sure on-screen dimensions are preserved in output file
                    print('-f3', '-r300', '-djpeg', strcat(infilename(1,1:(size(infilename,2)-4)),'_cand_',num2str(m),'_',channel,'_trace.jpg'));
                    save(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molintervals.mat'),'molintervals');
                    save(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molnum.mat'),'molnum');
                    save(strcat(infilename(1,1:(size(infilename,2)-4)),'_',channel,'_molcoords.mat'),'molcoords');
                    disp('Recorded molecule.');
                else
                    disp('Discarded.');
                end
                %             hold off
                %             user1 = input('Command? n = next, p = previous (n)', 's');
                pause(0.5);
                user1 = 'n';
                if strcmp(user1,'p')==1
                    m = m-1;
                    %                 hold off
                else
                    m = m+1;
                    %                 hold off
                end
            else
                %             hold off
                m = m+1;
            end
        else
            m = m+1;
            %         hold off
        end
        
    end
    
    
    
else
    
    disp('Error: Check averaging interval.  Movie was NOT analyzed.');
    
end