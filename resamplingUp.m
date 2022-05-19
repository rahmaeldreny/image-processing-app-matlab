% function [resizedImg]=resamplingUp(img,f);
% [m n matrixno]=size(img);
% for metricesIndex=1:1:matrixno
%     resizedImg(:,:,metricesIndex)=resamplingUp(img(:,:,metricesIndex),f);
% end
% function [img2]=resamplingUp(image,f)
% [m n]=size(image);
% x2=1;
% i2=1;
% for i=1:m
%     for l=1:f
%         for x=1:n
%             for l2=1:f
%                 img2(i2,x2)=image(i,x);
%                 x2=x2+l;
%             end
%         end
%         i2=i2+l;
%         x2=l;
%     end
% end
% end
% end
function [outImage] = upSampling(image, upSamplingFactor)
[rows cols] = size(image);
newRows = rows*upSamplingFactor;
newCols = cols*upSamplingFactor;
rowStart = 1;
for rowsIndex=1:upSamplingFactor:newRows
    colStart = 1;
    for columnIndex=1:upSamplingFactor:newCols
        outImage(rowsIndex:rowsIndex+upSamplingFactor-1,columnIndex:columnIndex+upSamplingFactor-1) = image(rowStart,colStart);
    colStart = colStart + 1;
    end
    rowStart = rowStart + 1;
end
