function [img] = readraw(filepath,width,height,pixfmt)
%READRAW Read Raw File saved by Imager
%   Detailed explanation goes here

if nargin ==3
    pixfmt='Mono8';
end

fid=fopen(filepath,'r');
raw=fread(fid,'*uint8');
fclose(fid);

switch pixfmt
    case 'Mono8'
        img = reshape(raw,width,height)';
    case 'Mono12Packed'
        % Unpack Mono12packed 3 bytes to 2 uint16
        % The two 12bits pixels | A | B | are packed in | A11 A10 A9 A8 A7 A6 A5 A4 | B3 B2 B1 B0 A3 A2 A1 A0 | B11 B10 B9 B8 B7 B6 B5 B4 | 3 bytes
        
        % unpacked bytes
        img = zeros(length(raw)/3*4,1,'uint8');
        % second bytes
        r2 = raw(2:3:end);
        % second bytes shift 4 bits left | A3 A2 A1 A0 0 0 0 0 |
        r2 = bitshift(r2,4);
        img(1:4:end)=r2;
        img(2:4:end)=raw(1:3:end);
        img(3:4:end)=raw(2:3:end);
        img(4:4:end)=raw(3:3:end);
        % shift uint16s 4 bits right
        img = bitshift(typecast(img,'uint16'),-4);
        img = reshape(img,width,height)';
    otherwise
        error(['Pixel Format: ',pixfmt,' is not supported.']);
end

end

