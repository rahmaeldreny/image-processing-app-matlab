function [outimg]=or_one(img,s,e)
[rows cols]=size(img);
for i=1:rows
   for j=1:cols
    if (img(i,j)>=s&&img(i,j)<=e)
        outimg(i,j)=img(i,j)||1;
    else
        outimg(i,j)=img(i,j);
    end
end
end
end
