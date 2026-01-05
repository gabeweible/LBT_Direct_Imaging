function getzone, cube, annulus, angles

; AD
; Dec 2015 - Undisclosed Location
;
;
; The .indices returns an array of the (linear) indices of the pixels
; that are extracted in the zone. This can be used to copy these back,
; after processing, to the final image.
; =======================

 sc = size(cube)

 mask  = fltarr(sc[1],sc[2])
 
 ; get the x- and y-dimensions of the image (number of px, not max. px index)
mask_sz_x = sc[1]
mask_sz_y = sc[2]

; star center
x_center = mask_sz_x / 2.
y_center = mask_sz_y / 2.

inrad = annulus[0] & outrad = annulus[1]

; Loop through all pixels in the image
for x = 0, mask_sz_x - 1 do begin
	for y = 0, mask_sz_y - 1 do begin
	
		; distance of pixel at (x, y) from center of the image
		distance = sqrt( ((x-x_center)^2.) + ((y-y_center)^2.) )
	
		;Pythagorean theorem baby! (everything within the radius of pixels is turned into NaN)
		if (distance lt inrad) or (distance) gt outrad then mask[x,y] = !values.f_nan
		
	endfor; y-pixel loop
endfor; x-pixel loop
 
 ang = mask
 ang = angmask(mask, sc[1]/2.-1, sc[2]/2.-1)

 good = where(ang ge angles[0]*!DTOR and ang le angles[1]*!DTOR and finite(mask) eq 1)

 indices=good
 data = fltarr(n_elements(good), sc[3])

                                ; Go through the images, layer by
                                ; layer, and extract the data points
                                ; required here
 for ii=0, sc[3]-1 do begin
    img= cube[*,*,ii]
    data[*,ii] = img[good]
 endfor

return, {data:data, good:good, num:n_elements(good), indices:indices}
end

