pro clean, obj_name, stripe, cube_folder, corr_thresh, half_cropped_sz
;'HII1348', '~/OneDrive/Research/HII1348/testing'
ct = corr_thresh 
compile_opt idl2

;HII 1348 (SX AND DX)
if stripe eq 'center1024' then for runs=1,4 do begin

	; Do this for runs eq 1 and runs eq 3
	if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'

	;Do this for runs eq 1 and runs eq 2
	if runs lt 3 then output_folder = cube_folder + '/processed_left/' else output_folder = cube_folder + '/processed_right/'

	print, 'Finding bad frames...'
	cube = readfits(output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', hdr)
	restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'

	;192,195

	bad_pixels = []
	goods = findgen((size(cube))[3])

	;initialize a new cube to clean
	cube2 = cube
	;remove saturated center
	;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.

	for xx = 0, (size(cube))[2] - 1 do begin
		for yy = 0, (size(cube))[2] - 1 do begin
			;Pythagorean theorem baby! (everything outside of r = 20px is divided into oblivion)
			if sqrt( ((xx - half_cropped_sz-1)^2.) + ((yy - half_cropped_sz-1)^2.) ) gt 20 then cube2[xx,yy,*] *= 0.0000001
		endfor
	endfor

	;compute pupil median
	lmed = median(cube2, dim=3)

	writefits, output_folder + dither_folder + obj_name + '_test-med_' + string(ct) +  '_corrthresh.fits', lmed

	;First frame special case
	maxcor = max(crosscorr(cube2[*,*,0], lmed))
	print, 'Frame ', 0, ' has max correlation value of ', maxcor, ' with pupil median.'
	;First frame vs. the rest
	corrs = maxcor

	for jj=1, (size(cube))[3] - 1 do begin
		maxcor = max(crosscorr(cube2[*,*,jj], lmed))
		print, 'Frame ', jj, ' has max correlation value of ', maxcor, ' with pupil median.'
		;First frame vs. the rest
		corrs = [corrs, maxcor]
	endfor; max correlation for

	bad_pixels = where(corrs lt corr_thresh)

	; Commenting out this stuff so that I don't have to worry about the plots...
	;plot, corrs - median(corrs)
	;level = corrs
	;level[*] = corr_thresh
	;print, level
	!p.linestyle = 2
	;oplot, level - median(corrs)
	!p.linestyle = 0

	print, bad_pixels
	print, 'Found ', n_elements(bad_pixels), ' bad frames.'
	print, 'Average correlation = ', mean(corrs)
	print, 'Median correlation = ', median(corrs)
	print, 'Standard deviation = ', stdev(corrs)
	;hak
	print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'

	if n_elements(bad_pixels) gt 1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + string(ct) +  '_corrthresh.fits', left_bad_pixels_cube
	endif else if n_elements(bad_pixels) eq 1 and bad_pixels ne -1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + string(ct) +  '_corrthresh.fits', left_bad_pixels_cube
	endif

	cube = cube[*,*,goods]; Reassign our cube to only have our good pixels
	medarr, cube, medframe; Get a median frame for our new bad-pixel-free cube

	;Write our clean cube, our bad pixels, and our median frame
	writefits, output_folder + dither_folder + obj_name + string(ct) +  '_cube_skysub_cen_clean.fits', cube, hdr
	writefits, output_folder + dither_folder + obj_name + string(ct) +  '_cube_skysub_cen_bad_pixels.fits', left_bad_pixels_cube, hdr
	writefits, output_folder + dither_folder + obj_name + string(ct) + '_pupil.fits', medframe
	save, filename = output_folder + dither_folder + obj_name + string(ct) +  '_parang_clean.sav', angles

	cgcleanup; Found online, closes any graphics or widget windows from our plots
	print, 'Run: ' + string(runs) + ' complete
endfor	; runs for

; ALCOR (SX ONLY)
if stripe eq 'second512' then for runs=1,2 do begin
; Do this for runs eq 1 and runs eq 3
	if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'

	output_folder = cube_folder + '/processed_left/'

	print, 'Finding bad frames...'
	cube = readfits(output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', hdr)
	restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'

	;192,195

	bad_pixels = []
	goods = findgen((size(cube))[3])

	;initialize a new cube to clean
	cube2 = cube
	;remove saturated center
	;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.

	for xx = 0, (size(cube))[2] - 1 do begin
		for yy = 0, (size(cube))[2] - 1 do begin
			;Pythagorean theorem baby! (everything within r = 20px is divided into oblivion)
			if sqrt( ((xx - half_cropped_sz)^2.) + ((yy - half_cropped_sz)^2.) ) gt 20 then cube2[xx,yy,*] *= 0.0000001
		endfor
	endfor

	;compute pupil median
	lmed = median(cube2, dim=3)

	writefits, output_folder + dither_folder + obj_name + '_test-med_' + string(ct) +  '_corrthresh.fits', lmed

	;First frame special case
	maxcor = max(crosscorr(cube2[*,*,0], lmed))
	print, 'Frame ', 0, ' has max correlation value of ', maxcor, ' with pupil median.'
	;First frame vs. the rest
	corrs = maxcor

	for jj=1, (size(cube))[3] - 1 do begin
		maxcor = max(crosscorr(cube2[*,*,jj], lmed))
		print, 'Frame ', jj, ' has max correlation value of ', maxcor, ' with pupil median.'
		;First frame vs. the rest
		corrs = [corrs, maxcor]
	endfor; max correlation for

	bad_pixels = where(corrs lt corr_thresh)

	; Commenting out this stuff so that I don't have to worry about the plots...
	;plot, corrs - median(corrs)
	;level = corrs
	;level[*] = corr_thresh
	;print, level
	!p.linestyle = 2
	;oplot, level - median(corrs)
	!p.linestyle = 0

	print, bad_pixels
	print, 'Found ', n_elements(bad_pixels), ' bad frames.'
	print, 'Average correlation = ', mean(corrs)
	print, 'Median correlation = ', median(corrs)
	print, 'Standard deviation = ', stdev(corrs)
	;hak
	print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'

	if n_elements(bad_pixels) gt 1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + string(ct) +  '_corrthresh.fits', left_bad_pixels_cube
	endif else if n_elements(bad_pixels) eq 1 and bad_pixels ne -1 then begin
		remove, bad_pixels, goods, angles
		left_bad_pixels_cube = cube[*,*,bad_pixels]
		writefits, output_folder + dither_folder + obj_name + '_test-bad_pixels_' + string(ct) +  '_corrthresh.fits', left_bad_pixels_cube
	endif

	cube = cube[*,*,goods]; Reassign our cube to only have our good pixels
	medarr, cube, medframe; Get a median frame for our new bad-pixel-free cube

	;Write our clean cube, our bad pixels, and our median frame
	writefits, output_folder + dither_folder + obj_name + string(ct) +  '_cube_skysub_cen_clean.fits', cube, hdr
	writefits, output_folder + dither_folder + obj_name + string(ct) +  '_cube_skysub_cen_bad_pixels.fits', left_bad_pixels_cube, hdr
	writefits, output_folder + dither_folder + obj_name + string(ct) + '_pupil.fits', medframe
	save, filename = output_folder + dither_folder + obj_name + string(ct) +  '_parang_clean.sav', angles

	cgcleanup; Found online, closes any graphics or widget windows from our plots
	print, 'Run: ' + string(runs) + ' complete
endfor

print, 'Finished!'

end
