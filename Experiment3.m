clear;
close all;
% Load the tiff file
filename = 'MbNB_Pos4004_z2-14_t1-341 cropped.tif';
[ I, n_f, n_s, n_c, row, col ] = loadimage4D( filename );
se = strel('disk',5);
%Initiat variables
imsum = zeros(row,col);
im_mean = zeros(row,col);
LastCellCount = 0;
timeLapse = 0;
t = 0;
oldTime=1;
% Find the Big Neuroblast center
for time = 1:n_f
    for Z = 1:n_s
        image = I(:,:,time,Z,2);
        imsum = imsum+image;
    end
    
    im_mean = im_mean + imsum/Z;
    imsum = zeros(row,col);
end
im_mean = im_mean/time;
bigcell = (im_mean>4 & im_mean<8);
bigcell = medfilt2(bigcell);
bigcell = imopen(bigcell,se);
bigcell = imfill(bigcell,'holes');
bigcell_stat = regionprops(bigcell,'Centroid');
big_cell_center = bigcell_stat.Centroid;
clear bigcell im_mean;
%Process each time frame
for time = 1:n_f
    diameters_array = [];
    centers_array=[];
    diameters_1 = zeros(1);
    centers_1 = zeros(1,2);
    diameters_2 = zeros(1);
    centers_2 = zeros(1,2);
    z_array = [];
    %Process each slice
    for Z = 1:n_s
        z_val = [];
        im = I(:,:,time,Z,2);
        %filtering and morphing
        im_morph = medfilt2(im);
        im_gray  = im_morph;
        thresh   = graythresh(im_gray);
        bwImage  = (im_gray>=20 & im_gray<=200);
        im_morph = imfill(bwImage,'holes');
        im_morph = imopen(im_morph,se);
        %Find the centers and diameters of each blob
        stats = regionprops('table',...
            im_morph,'Centroid','MajorAxisLength','MinorAxisLength');
        centers = stats.Centroid;
        diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
        %Get Z value of each blob
        if ~isempty(diameters)
            z_val(1:length(diameters))=Z;
        end
        %Collect Centers, Diameters and Z values of all the blobs in all
        %slices to be processed
        diameters_array = cat(1,diameters_array, diameters);
        centers_array = cat(1,centers_array,centers);
        if ~isempty(z_val)
            z_array = cat(2,z_array,z_val);
        end
        %Sum up all slices for better view
        imsum=imsum +im;
    end
    %Search for cell around the big Neuroblast only- Not near the center
    %nor far away from it and creat a new array
    ds=1;
    for d = 1:length(diameters_array)
        if diameters_array(d) < 30 && diameters_array(d) > 10 ...
                &&(centers_array(d,1) > (big_cell_center(1,1) - 50) ...
                && centers_array(d,2) < (big_cell_center(1,2) + 50) ...
                && centers_array(d,1) > (big_cell_center(1,1) - 50) ...
                && centers_array(d,2) < (big_cell_center(1,2) + 50) ...
                && ~(centers_array(d,1) > (big_cell_center(1,1) - 10) ...
                &&   centers_array(d,2) < (big_cell_center(1,2) + 10) ...
                &&   centers_array(d,1) < (big_cell_center(1,1) + 10) ...
                &&   centers_array(d,2) > (big_cell_center(1,2) - 10)))
            diameters_1(ds) = diameters_array(d)/2;
            centers_1(ds,:) = centers_array(d,:);
            z_1(ds)=z_array(d);
            ds= ds + 1;
        end
    end
    % Create a new array that only has the largest diameters of each blob
    % to get the where it is located in the Z dimension
    d=1;s=1;
    big_D = diameters_1(d);
    big_C = centers_1(d,:);
    big_Z = z_1(d);
    z_2=[0];
    while d <= length(diameters_1)
        while s <= length(diameters_1)
            if big_D < diameters_1(s)  ...
                    && (big_C(1) > centers_1(s,1) - 5 ...
                    &&  big_C(2) < centers_1(s,2) + 5 ...
                    &&  big_C(1) < centers_1(s,1) + 5 ...
                    &&  big_C(2) > centers_1(s,2) - 5)
                big_D = diameters_1(s);
                big_C = centers_1(s,:);
                big_Z = z_1(s);
                
                
            end
            s=s+1;
        end
        if big_D ~= diameters_2(length(diameters_2))
            centers_2=[centers_2;big_C];
            diameters_2 = [diameters_2;big_D];
            z_2 = [z_2;big_Z];
        end
        s=1;
        big_D = diameters_1(d);
        big_C = centers_1(d,:);
        big_Z = z_1(d);
        d=d+1;
    end
    %View the image and only keep the unique values of the diameters,
    %centers and their corresponding Z values
    figure(1);imagesc(imsum);
    centers_final=unique(centers_2,'rows');
    diameters_final = unique(diameters_2);
    centers_final(1,:)=[];
    diameters_final(1) = [];
    X=[];
    for x=1:length(diameters_final)
        X=cat(1,X,find(diameters_array==diameters_final(x)*2));
    end
    z_final=z_array(X);
    imsum = zeros(row,col);
    % Calcuate the time lapse between consecutive births
    if length(diameters_final)>LastCellCount
        t=t+1;
        if t>2
            t = 0;
            timeLapse=time - oldTime;
            LastCellCount = length(diameters_final);
            oldTime = time;
        end
    end
    % Mark the duaghter cells
    hold on
    viscircles(centers_final,diameters_final);
    for tex = 1:length(z_final)
        text(centers_final(tex,1), centers_final(tex,2), sprintf('%d', z_final(tex)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
    hold off;
    % <=======Set a Breakpoint here====
    clear diameters_1 ; clear centers_1 ;
    clear diameters_2; clear centers_2;
    clear diameters_final ; clear centers_final ;
    clear diameters_array ; clear centers_array ;
end

