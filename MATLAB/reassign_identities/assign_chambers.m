%Annika Rings July 2018
%function for reassigning chamber numbers and identities of flies that were
%tracked by Flytracker

%depends on the nestStructSort package - this must be added to the path or
%located in the folder from which the function is executed.

%USAGE: in matlab command window, type: assign_chambers(inputdir)
%inputdir is the absolute or relative path to the input directory

%reads in data from a file called 'inputdir-calibration.mat' located in
%inputdir
%reads additional data from a file called 'trx.mat' located in the folder
%called 'inputdir_JAABA' in inputdir
%saves a structure 'ids' in a file called 'ids.mat' in the inputdir

%the reassign identities function automatically calls this function

function assign_chambers(inputdir)
startdir = pwd;
cd(inputdir);
% calibfile=strcat(inputdir,'-calibration.mat');
calibfile = '../calibration.mat';
JAABAfolder = strcat(inputdir, '_JAABA');
trxfile = fullfile(JAABAfolder, 'trx.mat');
trxfile_new = strrep(trxfile, '.mat', '_id_corrected.mat');
%load calibration file
load(calibfile);
%load trx file
load(trxfile);
%determine old ids assigned by the tracker
id_old = arrayfun(@(f) f.id, trx);
%determine mean x and y values of each fly
meanx = arrayfun(@(f) mean(rmmissing(f.x)), trx);
meany = arrayfun(@(f) mean(rmmissing(f.y)), trx);
meanx_c = num2cell(meanx);
%add meanx and meany to trx file
[trx(:).meanx] = deal(meanx_c{:});
meany_c = num2cell(meany);
[trx(:).meany] = deal(meany_c{:});
%determine number of rows and columns of chambers
centroid = sortrows(calib.centroids, 1);
numchambers = size(centroid, 1);
diff_y_centroid = diff(centroid(:, 1));
av_increase_y = mean(diff_y_centroid);
y_jumps = find(diff_y_centroid > av_increase_y);
numrows = size(y_jumps, 1) + 1;
newrow_starts = 1 + [0; y_jumps; numchambers];


numcolumns = ceil(numchambers/numrows);
%sort the chambers by row first, then within the row sort by column
for c = 1:numrows
    first_in_row = newrow_starts(c);
    last_in_row = newrow_starts(c+1) -1 ;
    
    centroid(first_in_row:last_in_row, :) = sortrows(centroid(first_in_row:last_in_row, :), 2);
end


maxx_chambers = centroid(:, 2) + round((calib.rois{1}(4))/2);
minx_chambers = centroid(:, 2) - round((calib.rois{1}(4))/2);
maxy_chambers = centroid(:, 1) + round((calib.rois{1}(4))/2);
miny_chambers = centroid(:, 1) - round((calib.rois{1}(4))/2);
x = arrayfun(@(f) rmmissing(f.x), trx, 'UniformOutput', false);
firstx = cellfun(@(f) f(1, 1), x);

chambers = arrayfun(@(f) assign_one(f, maxy_chambers, miny_chambers, maxx_chambers, minx_chambers, numchambers), trx);
startpos = left_right_assign(chambers, firstx, numchambers);


id_new = assign_identity(chambers, calib.n_flies, numchambers);
id_new_c = num2cell(id_new);
[trx(:).id] = deal(id_new_c{:});
fields = {'meanx', 'meany'};
trx = rmfield(trx, fields);
trx = nestedSortStruct(trx, 'id');
[trx(:).startpos] = deal(startpos{:});

ids = struct('id_old', id_old, 'id_new', id_new, 'chambers', chambers);
save('ids', 'ids');
save(trxfile_new, 'trx');

cd(startdir);
clear all;
end


function chambernumber = assign_one(f, maxy_chambers, miny_chambers, maxx_chambers, minx_chambers, numchamb)
%assigns the chamber number according to the x value to the flies
for i = 1:numchamb
    if f.meany < maxy_chambers(i) && f.meany > miny_chambers(i) && f.meanx < maxx_chambers(i) && f.meanx > minx_chambers(i)
        chambernumber = i;
    end
end
end

function id_new = assign_identity(chambers, nflies, nchambers)
id_new = zeros(1, numel(chambers));
occupied_chambers = unique(chambers);

nums = [1:nflies * nchambers];
individuals = transpose(reshape(nums, [nflies, nchambers]));
for i = 1:numel(occupied_chambers)
    
    id_new(chambers == occupied_chambers(i)) = individuals(occupied_chambers(i), :);
    
end
end


function startpos = left_right_assign(chambers, firstx, numchambs)
%assigns startposition 'l' for left or 'r' for right to flies
%according to their position in the first frame of the tracking
startpos = {};
for i = 1:numchambs
    firstxi = firstx(chambers == i);
    if numel(firstxi) == 2
        if firstxi(1) < firstxi(2)
            
            startpos{(size(startpos, 1) + 1), 1} = {'l'};
            startpos{(size(startpos, 1) + 1), 1} = {'r'};
        else
            startpos{(size(startpos, 1) + 1), 1} = {'r'};
            startpos{(size(startpos, 1) + 1), 1} = {'l'};
        end
    elseif numel(firstxi) == 1
        startpos{(size(startpos, 1) + 1), 1} = {'l'};
    end
    
end
end
