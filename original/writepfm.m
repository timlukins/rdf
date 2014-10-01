function [] = writepfm(im,str)
%writepfm - write color PFM image
%
%writepfm(im,str)

[h,w,c] = size(im);

disp(sprintf('Writing Image as %s, [H,W] = [%d,%d]',str,h,w));

fp = fopen(str,'wb');
fprintf(fp, 'PF\n%d %d\n%f\n',  w , h, -1.0);

fdata = zeros(h*w*3,1);
if(c==1)
    im = flipud(im)';
    im = im(:);
    fdata(1:3:end) = im;
    fdata(2:3:end) = im;
    fdata(3:3:end) = im;
elseif(c==3)
    imr = flipud(im(:,:,1))';
    imr = imr(:);
    img = flipud(im(:,:,2))';
    img = img(:);
    imb = flipud(im(:,:,3))';
    imb = imb(:);
    fdata(1:3:end) = imr;
    fdata(2:3:end) = img;
    fdata(3:3:end) = imb;
else
    disp('Not a valid image');
end
fwrite(fp,fdata,'float');
fclose(fp);
