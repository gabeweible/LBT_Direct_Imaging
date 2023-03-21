FUNCTION DO_ADI, obj_name, cube_folder, use_injection, do_destripe, filter, suffix, ct,$
	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, nod=nod, runs=runs
	
compile_opt idl2
newline = string(10B)
	
	; Do this for runs eq 1 and runs eq 3
if runs mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
;Do this for runs eq 1 and runs eq 2
if runs lt 3 then begin
	output_folder = cube_folder + 'processed_left'
	truenorth = truenorth_sx
endif else begin
	truenorth = truenorth_dx
	output_folder = cube_folder + 'processed_right'
endelse


if use_injection then begin
   obj_cube = readfits(output_folder + dither_folder + obj_name + string(ct) +$
   	'_cube_skysub_cen_clean_inj.fits')
   inj_string = '_inj'
endif else begin
   obj_cube = readfits(output_folder + dither_folder + obj_name + string(ct) +$
   	'_cube_skysub_cen_clean.fits')
   inj_string = ''
endelse

n_frames = (size(obj_cube))[3]
if do_destripe then begin
   print, 'destriping 90 degrees...'
   
   for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 90.,$
   	clip_level=0.0, /nodisp)
   	
   print, 'destriping 0 degrees...', newline
   
   for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 0.,$
   	clip_level=0.0, /nodisp)
   	
endif; destripe if

; Subtract a smoothed version of the image
if filter gt 1 then for iii=0, n_frames - 1 do obj_cube[*,*,iii] -=$
	smooth(obj_cube[*,*,iii], filter)
	
; Normalize all of the frames so that they each have a max value of 1 (should help with residuals around the star)
if keyword_set(normal) then begin
   if normal eq 1 then for iv=0, n_frames - 1 do obj_cube[*,*,iv] /= max(obj_cube[*,*,iv])
endif

restore, filename = output_folder + dither_folder + obj_name + string(ct) +$
	'_parang_clean.sav'

medarr, obj_cube, medframe; medframe is the output here of calling medarr on our obj_cube

;Initialize our cube for Angular Differential Imaging
adi_cube = obj_cube

for ii=0, n_frames - 1 do begin
   if not keyword_set(silent) or silent eq 0 then print, 'Rotating by ', angles[ii],$
   	 newline
   ; Rotate obj_cube frame to truenorth
   ; (negatives because these are E of N, a.k.a. CCW, ROT does CW)
   obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
   ; Subtract the median frame (star PSF) from the new adi_cube frame
   adi_cube[*,*,ii] -= medframe
   ; Rotate the (PSF-subtracted) adi_cube frame to truenorth
   ; (negatives because these are E of N, a.k.a. CCW, ROT does CW)
   adi_cube[*,*,ii] = rot(adi_cube[*,*,ii], -angles[ii] - truenorth, /interp)
endfor; rotate if

suffix = cube_folder + 'combined/' + obj_name + 'ct_' + string(ct)$
	+ 'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj) + '_uncert_' +$
	string(uncert)

; Write the de-rotated obj_cube (No ADI)
writefits, strcompress(suffix + '_cube_skysub_cen_filt_derot.fits', /rem),$
	obj_cube

combine_type = 'nwadi'

if combine_type eq 'median' then medarr, obj_cube, medframe
if combine_type eq 'mean' then medframe = mean(obj_cube, dim=3)
if combine_type eq 'nwadi' then medframe = nw_ang_comb(obj_cube, angles)

writefits, strcompress(suffix + '_median_derot.fits', /rem), medframe

if combine_type eq 'median' then medarr, adi_cube, adiframe
if combine_type eq 'mean' then adiframe = mean(adi_cube, dim=3)
if combine_type eq 'nwadi' then adiframe = nw_ang_comb(adi_cube, angles)

adiframe[where(finite(adiframe) ne 1)] = 0.

writefits, strcompress(suffix + '_median_derot_adi.fits', /rem), adiframe

size = 500.
width = 8.72059;(3.8*1E-6) / (8.4) * 206265. / 0.0107; Where does this come from?
print, 'PSF Width: ', width
PSF = psf_Gaussian(npixel=size, FWHM=[width, width])
PSFN = PSF / MAX(PSF); N for normalized
adiframe_c = convolve(adiframe, PSFN)

print, 'Writing FITS for run: ', runs, ' ...'
writefits, strcompress(suffix + '_median_derot_adi_conv.fits', /rem), adiframe_c
print, 'FITS for run: ', runs, ' written!', newline

return, adiframe
END; do_adi function end