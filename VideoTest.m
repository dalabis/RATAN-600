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
end
telescope = info.PrimaryData.Keywords(10, 2);
instrument = info.PrimaryData.Keywords(11, 2);
%%
%Z = peaks;
%surf(Z); 
%axis tight manual 
%set(gca,'nextplot','replacechildren');

v = VideoWriter('peaks.avi');
open(v);

for i = 1:size(frames, 3)
   %surf(sin(2*pi*k/20)*Z,Z)
   imshow(1000-frames(:, :, i), [0 1000])
   text(1, 1, date(i), 'VerticalAlignment', 'top');
   text(1, size(frames(:, :, i), 1), telescope, 'VerticalAlignment', 'bottom');
   text(size(frames(:, :, i), 2), 1, instrument, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
   frame = getframe(gcf);
   
   for j = 1:10
        writeVideo(v,frame);
   end
end

close(v);