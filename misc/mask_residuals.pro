pro mask_residuals, image, radius, output_file
;+
; NAME:
;       MASK_RESIDUALS

; PURPOSE:
;       Reduce PSF-subtracted residuals in an the center of an image to zero,
;       within a radius

; CALLING SEQUENCE:
;       MASK_RESIDUALS, image, radius, output_file

; INPUTS:
;       image = a single 2D FITS image
; 
;       radius = radius in pixels to mask in a circle

;       output_file = string output FITS to write the masked image to
;       
; RESTRICTIONS:
;       (1) image must be a string path to the FITS image, radius may be a float
;           or an integer.
;
; EXAMPLE:
;       Mask the center of an image within a radius of 10 pixels:
;       MASK_RESIDUALS, '/path/to/image.fits', 10, '/path/to/masked/image.fits'
;
; EXTERNAL PROCEDURES USED:
;       READFITS, WRITEFITS
;
; MODIFICATION HISTORY:
;       WRITTEN, 2024 Gabriel Weible
;-



; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

; read-in the unmasked image
image = readfits(image)

; get the x- and y-dimensions of the image (number of px, not max. px index)
image_sz_x = (size(image))[1]
image_sz_y = (size(image))[2]

; Loop through all pixels in the image
for x = 0, image_sz_x - 1 do begin
	for y = 0, image_sz_y - 1 do begin
	
		;Pythagorean theorem baby! (everything within the radius of pixels is divided into oblivion)
		if sqrt( ((x-(image_sz_x-1)/2)^2.) + ((y-(image_sz_y-1)/2)^2.) ) lt radius then image[x,y] *= 0.00000000000000000000000001
		
	endfor; y-pixel loop
endfor; x-pixel loop

; write the masked image
writefits, output_file, image

end