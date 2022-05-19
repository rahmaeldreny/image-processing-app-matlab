function [outimg]=logic(img,s,e)
[rows cols ]=size(img);

%gimg=rgb2gray(img);
for i=1:rows
   for j=1:cols
    if (img(i,j)>=s&&img(i,j)<=e)
        outimg(i,j)=img(i,j)&&1;
    else
        outimg(i,j)=img(i,j);
    end
end
end
end
% function [outimg]=threshold(img);
% img=getimage(handles.axes2);
% k=inputdlg('value of k=');
% k=str2num(k{1});
% [rows cols matreciesNo]=size(img);
% doubleImg=double(img);
% for i=1:rows
%     for j=1:cols
%         if doubleImg(i,j)>=k
%             doubleImg(i,j)=255;
%         else
%             doubleImg(i,j)=0;
%         end
%     end
% end
% s=uint8(doubleImg);
% axes(handles.axes2);
% imshow(s);
% 
% function [outimg]=and_zero(image,s,e)
% [rows cols]=size(img);
% for i=1:rows
%    for j=1:cols
%     if (image(i,j)>=s&&image(i,j)<=e)
%         outimg(i,j)=image(i,j)&&0;
%     else
%         outimg(i,j)=image(i,j);
%     end
% end
% end
% end
% 
% function [outimg]=or_one(image,s,e)
% [rows cols]=size(img);
% for i=1:rows
%    for j=1:cols
%     if (image(i,j)>=s&&image(i,j)<=e)
%         outimg(i,j)=image(i,j)||1;
%     else
%         outimg(i,j)=image(i,j);
%     end
% end
% end
% end
% 
% 
% function [outimg]=or_zero(image,s,e)
% [rows cols]=size(img);
% for i=1:rows
%    for j=1:cols
%     if (image(i,j)>=s&&image(i,j)<=e)
%         outimg(i,j)=image(i,j)||0;
%     else
%         outimg(i,j)=image(i,j);
%     end
% end
% end
% end