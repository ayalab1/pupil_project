
figure(1);
[X,map]=imread('trees.tif');% built in image
gry = double(rgb2gray(ind2rgb(X,map)));

subplot(2,2,1);
imshow(uint8(gry*255));
title('Original gray');

subplot(2,2,2);
col = gray2pseudocolor(gry,'traditional');
imshow(uint8(col*255));
title('traditional');

subplot(2,2,3);
col = gray2pseudocolor(gry,'ocean');
imshow(uint8(col*255));
title('ocean');

subplot(2,2,4);
col = gray2pseudocolor(gry,'rainbow');
imshow(uint8(col*255));
title('rainbow');


figure(2);

ax1=subplot(2,2,1);
surf(peaks, 'EdgeColor', 'none');
colormap(ax1,colormapformulae('jet'));
colorbar;
title('jet');

ax2=subplot(2,2,2);
surf(peaks, 'EdgeColor', 'none');
colormap(ax2,colormapformulae([7,5,15])); % it is same as 'traditional'
colorbar;
title('traditional');

ax3=subplot(2,2,3);
surf(peaks, 'EdgeColor', 'none');
colormap(ax3,colormapformulae([23,28,3])); % it is same as 'ocean'
colorbar;
title('ocean');

ax4=subplot(2,2,4);
surf(peaks, 'EdgeColor', 'none');
colormap(ax4,colormapformulae([34,35,36])); % it is same as 'black-red-yellow-white'
colorbar;
title('black-red-yellow-white');



