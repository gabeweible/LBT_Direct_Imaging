pro bad_pixels, output_folder, obj_name, stripe
;'~/OneDrive/Research/HII1348/testing/', 'HII1348' for current run
compile_opt idl2

;+
; NAME:
;       BAD_PIXELS
;       
; PURPOSE:
;       Finds and removes bad pixels from a data cube of FITS files,
;       writing a mask and a new cube (both as FITS files as well)
;       
; CALLING SEQUENCE:
;       BAD_PIXELS, output_folder, obj_name
;
; INPUTS:
;       OUTPUT_FOLDER = String folder to write the bad-pixel-free data
;       cube FITS into
;       
;       OBJ_NAME = String object name to look for in FITS header
;
; OPTIONAL INPUT:
;       None
; OPTIONAL INPUT KEYWORD:
;       None
; OUTPUTS:
;       None
; RESTRICTIONS:
;       (1) Make sure that output_folder and obj_name are string literals
;
; EXAMPLE:
;       Remove bad pixels from a data cube for Vega in ~/Desktop/Cube:
;       
;       BAD_PIXELS, '~/Desktop/Cube', 'Vega'
;
; PROCEDURES USED:
;       READFITS, WRITEFITS
;
; MODIFICATION HISTORY:
;       WRITTEN, Kevin Wagner ????
;       Moved into this procedure and edited by Gabriel Weible 2021
;-

;Read in our object cube created above from our individual raw FITS files
obj_cube = readfits(output_folder + obj_name + '_cube.fits')
restore, output_folder + obj_name + '_parang.sav'
print, 'Object cube size:', size(obj_cube)

if (size(obj_cube))[3] gt 250 then badpix = mean(obj_cube[*,*,0:250], dim=3); mean image for the first 251(?) frames
if (size(obj_cube))[3] lt 251 then badpix = mean(obj_cube[*,*,0:(size(obj_cube))[3]-1], dim=3); mean image for the first 251(?) frames

 ; Initialize the mask with all 1s
if stripe eq 'center1024' then badpix_mask = replicate(1, 2048, 1024)
if stripe eq 'second512' then badpix_mask = replicate(1, 2048, 512)

boxhsize = 8 & threshold = 5 & n_bad_pixels = 0. & n_tests = 0. ; 5 factor threshold to conclude a pixel is bad
;Trying a factor of 4...
for xx = 1.*boxhsize, (size(obj_cube))[1]-1-(1.*boxhsize) do begin; (size(obj_cube))[2] is our image width (x dimension)
   for yy = 1.*boxhsize,(size(obj_cube))[2]-1-(1.*boxhsize) do begin; (size(obj_cube))[1] is our image height (y dimension)
      test = badpix
      test[xx,yy] = !values.f_nan
      box = test[xx-boxhsize:xx+boxhsize-1, yy-boxhsize:yy+boxhsize-1]
      boxstdv = stddev(box, /nan)
      this = badpix[xx,yy] - mean(box, /nan)

      ;Do this stuff to find a bad pixel
      if abs(this) / boxstdv gt threshold then begin
         badpix_mask[xx, yy] = 0; 0 means a bad pixel!
         ;print, 'Bad pixel found!'
         n_bad_pixels += 1.
      endif; bad pixel locating if

      n_tests += 1.

      ;Print our progress every 100 frames (every frame is just too jittery) p.s., I should probably make a real widget for this sort of thing :P
         print, 'Testing pixel at (x, y): ', fix(xx), fix(yy)
         print, 'n_bad_pixels: ', fix(n_bad_pixels), ' n_tests: ', long(n_tests)
         print, 'Progress: ', 100. * n_tests / ((size(obj_cube))[1] * (size(obj_cube))[2]), '%'
         print, 'Bad Pixel Percent =', 100. * n_bad_pixels / n_tests, '%' + string(10B)
   endfor; yy for

   ;Write the bad pixel mask based on the given y-value that we're on in our loop (why?)
   ;writefits, strcompress(output_folder + obj_name + '_new_badpix_mask.fits', /rem), badpix_mask
endfor; xx for

;badpix_mask[where(badpix_mask eq 0.)]=2.
;badpix_mask[where(badpix_mask ne 2.)]=0.
;badpix_mask[where(badpix_mask eq 2.)]=1.
;badpix_mask=fix(badpix_mask)

; Write the bad pixel mask for a final time since we've looped through every pixel now
print, 'Writing bad pixel mask...'
writefits, strcompress(output_folder + obj_name + '_new_badpix_mask.fits', /rem), badpix_mask

n_frames = (size(obj_cube))[3]
for ii=0,n_frames-1 do begin; go through each frame in the data cube
   print, 'Working on frame index', ii, ' / ', n_frames-1
   frame = obj_cube[*,*,ii]; get our frame to fix
   fixed_frame = frame; initialize a frame that will be our fixed frame soon!
   ;frame[where(badpix_mask gt 0)]=!values.f_nan
   fixpix, frame, badpix_mask, fixed_frame, npix=20;, /silent; Fix the frame (why 20?)
   obj_cube[*,*,ii] = fixed_frame; reassign our frame in the obj cube to the fixed frame
endfor; fix frame for

; write our new FITS with the bad pixels removed
print, 'Writing cube with bad pixels removed...'
writefits, output_folder + obj_name + '_cube_no_bad_pixels.fits', obj_cube
print, 'Cube written!'

end