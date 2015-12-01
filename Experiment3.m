%clear;
close all;
filename = 'MbNB_Pos4004_z2-14_t1-341 cropped.tif';
%[ I, n_f, n_s, n_c, row, col ] = loadimage4D( filename );
se = strel('disk',5);
imsum = zeros(row,col);
im_mean = zeros(row,col);
clear diameters_array;
clear centers_array;
clear diameters_1;
clear centers_1;
clear diameters_2;
clear centers_2;
diameters_1 = zeros(1);
centers_1 = zeros(1,2);
diameters_2 = zeros(1);
centers_2 = zeros(1,2);
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
%imagesc(bigcell);
bigcell_stat = regionprops(bigcell,'Centroid');
big_cell_center = bigcell_stat.Centroid;
clear bigcell im_mean;
for time = 88:n_f
    diameters_array = [];
    centers_array=[];
    for Z = 1:n_s
        im = I(:,:,time,Z,2);
        im_morph = medfilt2(im);
        im_gray  = im_morph;
        thresh   = graythresh(im_gray);
        bwImage  = (im_gray>=20 & im_gray<=200);
        %(im_gray>=thresh*max(max(im_gray)));%
        im_morph = imfill(bwImage,'holes');
        im_morph = imopen(im_morph,se);
        figure(1);imagesc(im_morph);
        stats = regionprops('table',...
            im_morph,'Centroid','MajorAxisLength','MinorAxisLength');
        centers = stats.Centroid;
        diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
        diameters_array = cat(1,diameters_array, diameters);
        centers_array = cat(1,centers_array,centers);
        imsum=imsum +im;
    hold on
    viscircles(centers,diameters/2);
    hold off;
    end

    ds=1;
    for d = 1:length(diameters_array)
        if diameters_array(d) < 30 && diameters_array(d) > 10 ...
                &&(centers_array(d,1) > (big_cell_center(1,1) - 40) ...
                && centers_array(d,2) < (big_cell_center(1,2) + 40) ...
                && centers_array(d,1) > (big_cell_center(1,1) - 40) ...
                && centers_array(d,2) < (big_cell_center(1,2) + 40) ...
                && ~(centers_array(d,1) > (big_cell_center(1,1) - 5) ...
                && centers_array(d,2) < (big_cell_center(1,2) + 5) ...
                && centers_array(d,1) < (big_cell_center(1,1) + 5) ...
                && centers_array(d,2) > (big_cell_center(1,2) - 5)))
            diameters_1(ds) = diameters_array(d)/2;
            centers_1(ds,:) = centers_array(d,:);
            ds= ds + 1;
        end
    end

    d=1;s=1;
    centers_2=[0,0];
    diameters_2 = [0];
    big_D = diameters_1(d);
    big_C = centers_1(d,:);
    while d <= length(diameters_1)
        while s <= length(diameters_1)
            if big_D < diameters_1(s)  ...
                    && (big_C(1) > centers_1(s,1) - 5 ...
                    &&  big_C(2) < centers_1(s,2) + 5 ...
                    &&  big_C(1) < centers_1(s,1) + 5 ...
                    &&  big_C(2) > centers_1(s,2) - 5)
                big_D = diameters_1(s);
                big_C = centers_1(s,:);
                
            end
            s=s+1;
        end
        if big_D ~= diameters_2(length(diameters_2))
            centers_2=[centers_2;big_C];
            diameters_2 = [diameters_2;big_D];
        end
        s=1;
        big_D = diameters_1(d);
        big_C = centers_1(d,:);
        d=d+1;
        
    end
    
    figure(3);imagesc(imsum);
    centers_final=unique(centers_2,'rows');
    diameters_final = unique(diameters_2);
    imsum = zeros(row,col);
    hold on
    viscircles(centers_final,diameters_final);
    hold off;
    clear diameters_1 ; clear centers_1 ;
    clear diameters_2; clear centers_2;
    clear diameters_final ; clear centers_final ;
    clear diameters_array ; clear centers_array ;
end

