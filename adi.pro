pro adi, obj_name, cube_folder, use_injection, do_destripe, filter, suffix, ct,$
	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, norm=norm, uncert=uncert,$
	silent=silent, ssh=ssh
; silent will supress all the "Rotating by ..." prints.
compile_opt idl2
newline = string(10B)

if ssh eq 0 then path='/Users/gabeweible/OneDrive/research/HII1348/macbook_'
if ssh eq 1 then path='/Users/gabe/research/macbook_'

for run=1,4 do begin
   
; Do this for run eq 1 and run eq 3
if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
;Do this for run eq 1 and run eq 2
if run lt 3 then begin
   output_folder = cube_folder + 'processed_left'
   truenorth = -1.39
endif else begin
   truenorth = 0.59
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
if filter gt 1 then for iii=0, n_frames - 1 do obj_cube[*,*,iii] -= $
	smooth(obj_cube[*,*,iii], filter)
; Normalize all of the frames so that they each have a max value of 1 (should help with residuals around the star)
if keyword_set(norm) then begin
   if norm eq 1 then for iv=0, n_frames - 1 do obj_cube[*,*,iv] /= $
   	max(obj_cube[*,*,iv])
endif

restore, filename = output_folder + dither_folder + obj_name + string(ct) +$
	'_parang_clean.sav'

medarr, obj_cube, medframe; medframe is the output here of calling medarr on our obj_cube

;Initialize our cube for Angular Differential Imaging
adi_cube = obj_cube

for ii=0, n_frames - 1 do begin
   if not keyword_set(silent) or silent eq 0 then print, 'Rotating by ', angles[ii], newline
   obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
   adi_cube[*,*,ii] -= medframe
   adi_cube[*,*,ii] = rot(adi_cube[*,*,ii], -angles[ii] - truenorth, /interp)
endfor; rotate if

writefits, output_folder + dither_folder + obj_name + inj_string + 'ct_' + $
	string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj) + $
	'_uncert_' + string(uncert) +  '_cube_skysub_cen_filt_derot.fits', obj_cube

combine_type = 'nwadi'

if combine_type eq 'median' then medarr, obj_cube, medframe
if combine_type eq 'mean' then medframe = mean(obj_cube, dim=3)
if combine_type eq 'nwadi' then medframe = nw_ang_comb(obj_cube, angles)

writefits, strcompress(output_folder + dither_folder + obj_name + inj_string + $
	'ct_' + string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj)$
 	+ '_uncert_' + string(uncert) +  '_median_derot.fits', /rem), medframe

if combine_type eq 'median' then medarr, adi_cube, adiframe
if combine_type eq 'mean' then adiframe = mean(adi_cube, dim=3)
if combine_type eq 'nwadi' then adiframe = nw_ang_comb(adi_cube, angles)

adiframe[where(finite(adiframe) ne 1)] = 0.

writefits, strcompress(output_folder + dither_folder + obj_name + inj_string + $
	'ct_' + string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + $
   string(neg_inj) + '_uncert_' + string(uncert) +  '_median_derot_adi.fits', /rem), adiframe

size = 500.
width = 8.72059;(3.8*1E-6) / (8.4) * 206265. / 0.0107; Where does this come from?
print, 'PSF Width: ', width
PSF = psf_Gaussian(npixel=size, FWHM=[width, width])
PSFN = PSF / MAX(PSF); N for normalized
adiframe_c = convolve(adiframe, PSFN)

print, 'Writing FITS for run: ', run, ' ...'
writefits, strcompress(output_folder + dither_folder + obj_name + inj_string + $
	'ct_' + string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj)$
 	+ '_uncert_' + string(uncert) +  '_median_derot_adi_conv.fits', /rem), adiframe_c
print, 'FITS for run: ', run, ' written!', newline

if run eq 1 then adinods = adiframe else adinods = [[[adinods]], [[adiframe]]]

; I'm having trouble with trying to *not* manually type in the folder here (and in the lines below)
writefits, strcompress(path + string(coadd) + '/combined/' + obj_name + '_nods_adi'$
	 + suffix + 'ct_' + string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + $
	 string(neg_inj) + '_uncert_' + string(uncert) +  '.fits', /rem), adinods

if run eq 4 then begin
   e = mean(adinods, dim=3)
   writefits, strcompress(path + string(coadd) + '/combined/' + obj_name + $
   'ct_' + string(ct) +   'filt_'  +  string(filter) + '_neg_inj_' + $
   string(neg_inj) + '_uncert_' + string(uncert) +  '_total_adi.fits', /rem), e
   ;What is the correction factor for?
   ;
   ; This line is giving me problems so I'm going to try and comment it out...Something is going wrong so that tests is undefined
   ref_file=output_folder+dither_folder+obj_name+'_pupil.fits' ;unsaturated data can just use the pupil image
   
   if keyword_set(fs) then begin
   find_sources,strcompress(path + string(coadd) +$
       '/combined/' +obj_name + 'ct_' + string(ct) +   'filt_'  +  string(filter)$
        		+ '_neg_inj_' +$
        string(neg_inj) + '_total_adi.fits',/rem),reference=output_folder+$
        		dither_folder+obj_name+$
        string(ct) +  '_pupil.fits',platescale=0.0107,$
        correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
        fwhm=8.7,ct=ct,$
        filter=filter
   endif; fs set if

   e = mean(adinods[*,*,0:1], dim=3)
   writefits, strcompress(path + string(coadd) + '/combined/' + obj_name + $
   '_left_adi' + suffix + 'ct_' + string(ct) +   'filt_'  +  string(filter) + $
   '_neg_inj_' + string(neg_inj) + '_uncert_' + string(uncert) +  '.fits', /rem), e

   e = mean(adinods[*,*,2:3], dim=3)
   writefits, strcompress(path + string(coadd) + '/combined/' + obj_name + $
   '_right_adi' + suffix + 'ct_' + string(ct) +   'filt_'  +  string(filter) +$
    '_neg_inj_' + string(neg_inj) + '_uncert_' + string(uncert) +  '.fits', /rem), e
endif; last run (run eq 4)

endfor; run for

print, 'Finished.'

end
