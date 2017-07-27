function s = sum_square(y1,x1,y2,x2,r,im1,im2)
    patch1 = im1(y1-r:y1+r,x1-r:x1+r);
    patch2 = im2(y2-r:y2+r,x2-r:x2+r);
    dif = patch1-patch2;
    s = sum(sum(dif.^2));
end