pro bad_pixels, output_folder, obj_name, stripe, sigma_threshold=sigma_threshold,$
	boxhsize=boxhsize,pre_sky_sub=pre_sky_sub, bad_px_arr=bad_px_arr,$
	n_groups=n_groups, fix_bad_cols=fix_bad_cols;, do_destripe=do_destripe; removed
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

		;Stripe = 'center1024'(HII 1348) or 'second512' (Alcor)
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
if not keyword_set(pre_sky_sub) then pre_sky_sub = 0

; don't split the bad-pixel correction into multiple groups by default
if not keyword_set(n_groups) then n_groups=1

if pre_sky_sub eq 1 then obj_cube = readfits(output_folder + obj_name + '_cube.fits')
if pre_sky_sub eq 0 then obj_cube = readfits(output_folder + obj_name + '_cube_skysub.fits')

restore, output_folder + obj_name + '_parang.sav'
print, 'Object cube size:', size(obj_cube)

;if (size(obj_cube))[3] gt 51 then badpix = mean(obj_cube[*,*,1:51], dim=3, /double); mean image for the first 250 frames

; calculate # of frames per group based on n_groups
total_frames = (size(obj_cube))[3]
frames_per_group = total_frames / n_groups

; how many stddevs from mean to mask a bad pixel?
if not keyword_set(sigma_threshold) then threshold = 5
if keyword_set(sigma_threshold) then threshold = sigma_threshold

; how may pixels for half a box length?
if not keyword_set(boxhsize) then boxhsize = 4

for group=0,n_groups-1 do begin

	; Initialize the mask with all 1s (assume pixels are good by default)
	if stripe eq 'center1024' then badpix_mask = replicate(1, 2048, 1024)
	if stripe eq 'second512' then badpix_mask = replicate(1, 2048, 512)

	n_bad_pixels = 0. & n_tests = 0.

	start_i = group*frames_per_group
	
	; decide whether we get a full or a partial group size (only relevant at the end)
	if (group+1)*frames_per_group le total_frames then begin; full group
	
		end_i = (group+1)*frames_per_group-1
		; mean image for all the frames in this group
		badpix = mean(obj_cube[*,*,start_i:end_i], dim=3, /double, /nan)
		
	endif else begin; partial group (only at the very end, sometimes)
		
		end_i = total_frames-1
		; mean image for all the frames in this group
		badpix = mean(obj_cube[*,*,start_i:end_i], dim=3, /double, /nan)
		
	endelse
	
	; loop through pixels to test (given by boxhsize and the 2-D shape of the frames)
	for xx = 1.*boxhsize, (size(obj_cube))[1]-1-(1.*boxhsize) do begin; (size(obj_cube))[2] is our image width (x dimension)
	   	for yy = 1.*boxhsize,(size(obj_cube))[2]-1-(1.*boxhsize) do begin; (size(obj_cube))[1] is our image height (y dimension)
	   		test = badpix; replicate frame
		  	test[xx,yy] = !values.f_nan; set value @ test pixel location to NAN
		  	
		  	; grab a box out to boxhsize in all 4 directions to test
		  	box = test[xx-boxhsize:xx+boxhsize-1, yy-boxhsize:yy+boxhsize-1]
		  
		  	; standard deviation of all the pixels besides the test location
		  	boxstdv = stddev(box, /nan, /double)
		  
		  	; distance of actual value at (xx, yy) from the surrounding mean
		  	this = badpix[xx,yy] - mean(box, /nan, /double)
		  	;this = badpix[xx,yy] - median(box, /even, /double)
	
		  	;Do this stuff if a bad pixel is found (more stdevs away from mean than threshold)
		  	if (abs(this) / boxstdv) gt threshold then begin
				badpix_mask[xx, yy] = 0; 0 means a bad pixel!
			 	;print, 'Bad pixel found!'
			 	n_bad_pixels += 1.
			 	;Print our progress every 100 frames (every frame is just too jittery) p.s., I should probably make a real widget for this sort of thing :P
				print, 'Group #: ', group+1
				;print, 'Testing pixel at (x, y): ', fix(xx), fix(yy)
				;print, 'n_bad_pixels: ', fix(n_bad_pixels), ' n_tests: ', long(n_tests)
				print, 'Progress: ', 100. * n_tests / ((size(obj_cube))[1] * (size(obj_cube))[2]), '%'
				print, 'Bad Pixel Percent =', 100. * n_bad_pixels / n_tests, '%' + string(10B)
		  	endif; bad pixel locating if
		  
		  	if not finite(this) then begin
				badpix_mask[xx, yy] = 0; 0 means a bad pixel!
				;print, 'Bad pixel found!'
				n_bad_pixels += 1.
			endif
		  
		  	if (this lt -1000) or (this gt 1000) then begin
				badpix_mask[xx, yy] = 0; 0 means a bad pixel!
				;print, 'Bad pixel found!'
				n_bad_pixels += 1.
		  	endif; lt -1000
	
		  	n_tests += 1.
   		endfor; yy for
	endfor; xx for
	
	; START FIXING PIXELS IN THIS GROUP!!!
	; hard-coded bad pixels
	if keyword_set(bad_px_arr) then begin
		foreach bad_px,bad_px_arr do begin
			badpix_mask[bad_px[0],bad_px[1]]=0
		endforeach
	endif
	
	; known bad columns, etc.
	if keyword_set(fix_bad_cols) then begin
		badpix_mask[504,*]=0
		badpix_mask[512,*]=0
		badpix_mask[568,*]=0
		badpix_mask[576,*]=0
		badpix_mask[384,*]=0
		badpix_mask[960,*]=0
		badpix_mask[1024,*]=0
		badpix_mask[1080,*]=0
		badpix_mask[1272,*]=0
		badpix_mask[1280:1282,*]=0
		badpix_mask[1344,*]=0
		badpix_mask[1336,*]=0
		badpix_mask[1408,*]=0
		badpix_mask[1472,*]=0
	endif
	
	; Write the bad pixel mask for a final time since we've looped through every pixel now
	print, 'Writing bad pixel mask...'
	writefits, strcompress(output_folder + obj_name + '_new_badpix_mask_'+string(group)+'.fits', /rem), badpix_mask
	
	n_frames = end_i - start_i + 1; in this group
	for ii=start_i,start_i+n_frames-1 do begin; go through each frame in the data cube
	   print, 'Working on frame index', ii, ' / ', start_i + n_frames-1
	   frame = obj_cube[*,*,ii]; get our frame to fix
	   fixed_frame = frame; initialize a frame that will be our fixed frame soon!
	   ;frame[where(badpix_mask gt 0)]=!values.f_nan
	   fixpix, frame, badpix_mask, fixed_frame, npix=35;, /silent; Fix the frame
	   obj_cube[*,*,ii] = fixed_frame; reassign our frame in the obj cube to the fixed frame
	endfor; fix frame for
	
endfor; n_groups foreach


; write our new FITS with the bad pixels removed
print, 'Writing cube with bad pixels removed...'
writefits, output_folder + obj_name + '_cube_no_bad_pixels.fits', obj_cube
print, 'Cube written!'

end