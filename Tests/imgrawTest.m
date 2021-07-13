function [img,raw,rawimg] = imgrawTest(imgpath,rawpath,ispacked)
%IMGRAWTEST Get Image Array from image file and raw file saved by Imager from the same image buffer

imgpath = 'C:\Users\fff00\Pictures\Test\HalfMono12Packed-Epoch0-Frame0.TIFF';
rawpath = 'C:\Users\fff00\Pictures\Test\HalfMono12Packed-Epoch0-Frame0.Raw';
ispacked = true;

img = imread(imgpath);
fid=fopen(rawpath,'r');
raw=fread(fid,'*uint8');
fclose(fid);

if ispacked
    rawimg = zeros(numel(img),1,'like',img);
    j=0;
    for i = 1:3:length(raw)
        rawimg((1:2)+2*j)=unpackb(raw(i:i+2));
        j=j+1;
    end
    rawimg=reshape(rawimg,size(img))';
else
    rawimg=reshape(raw,size(img))';
end

%% unpack Mono12packed 3 bytes to 2 uint16 ( the two 12bits values | A | B | are packed in | A11 A10 A9 A8 A7 A6 A5 A4 | B3 B2 B1 B0 A3 A2 A1 A0 | B11 B10 B9 B8 B7 B6 B5 B4 | 3 bytes )
    function [v1,v2]=unpack(b3)
        b = de2bi(b3,8);
        v1=bi2de([b(2,1:4) b(1,:) 0 0 0 0]);
        v2=bi2de([b(2,5:8) b(3,:) 0 0 0 0]);
    end

    function [v1,v2]=unpackb(b3)
        v1= bitshift( typecast([bitshift( b3(2),4), b3(1)],'uint16'),-4);
        v2= bitshift( typecast([b3(2), b3(3)],'uint16'),-4);
    end
end

