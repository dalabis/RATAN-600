%%
% This script creates a .avi video file from a series of .fits files
% located at ..\200110_21_23\MG1_200110_21_23
directory = '200110_21_23\MG1_200110_21_23';
dirData = dir(directory);
for i = 3:size(dirData, 1)
    filename = dirData(i).name;
    frame = fitsread([directory, '\', filename]);
    frames(1:size(frame, 1), 1:size(frame, 2), i-2) = frame;
    
    info = fitsinfo([directory, '\', filename]);
    date(i-2) = info.PrimaryData.Keywords(9, 2);
    CRPIX1(i-2) = cell2mat(info.PrimaryData.Keywords(21, 2));
    CRPIX2(i-2) = cell2mat(info.PrimaryData.Keywords(22, 2));
end
telescope = info.PrimaryData.Keywords(10, 2);
instrument = info.PrimaryData.Keywords(11, 2);
SOLAR_R = cell2mat(info.PrimaryData.Keywords(44, 2));
SC_ROLL = cell2mat(info.PrimaryData.Keywords(48, 2));

%%
v = VideoWriter('peaks.avi');
open(v);

framesCentered(1:4*SOLAR_R, 1:4*SOLAR_R, 1:size(frames, 3)) = 0;
for i = 1:size(frames, 3)
    framesCentered((2*SOLAR_R - CRPIX1(i) + 1):(2*SOLAR_R - CRPIX1(i) + size(frames, 1)), (2*SOLAR_R - CRPIX2(i) + 1):(2*SOLAR_R - CRPIX2(i) + size(frames, 2)), i) = frames(:, :, i);
    framesCentered(:, :, i) = flipud(framesCentered(:, :, i));
    framesCentered(:, :, i) = imrotate(framesCentered(:, :, i), SC_ROLL, 'crop');
end

for i = 1:size(frames, 3)
    imshow(1000 - framesCentered(fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), i), [0 1000])
    text(1, 1, date(i), 'VerticalAlignment', 'top');
    text(1, size(framesCentered(fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), 2), i), telescope, 'VerticalAlignment', 'bottom');
    text(size(framesCentered(fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), fix(0.7*SOLAR_R):fix(3.3*SOLAR_R), 1), i), 1, instrument, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    frame = getframe(gcf);
    
    for j = 1:10
        writeVideo(v,frame);
    end
end

close(v);