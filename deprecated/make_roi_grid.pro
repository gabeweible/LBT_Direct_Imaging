function make_roi, xc0,yc0,half
    xc = fltarr(csize,csize)
    yc = fltarr(csize,csize)
    for i=0,csize-1 do begin xc[i,*] = xc0 + i; yc[i,*] = yc0 + indgen(csize) endfor
    return, {xc:xc, yc:yc}