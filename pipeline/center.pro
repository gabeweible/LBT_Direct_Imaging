pro center, obj_name, stripe, cube_folder, do_block_right, do_block_left,$
	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub, half_cropped_sz
; 'HII1348', '~/OneDrive/Research/HII1348/testing', 1, 1, 1, 1, 1, 0 for current run (let's try cen_median_sub later)
compile_opt idl2

; I think this is the second time that the runs matter (we are now post-split)
; HII 1348 CENTER 1024 PIXEL ROWS (DOUBLE-SIDED)
if stripe eq 'center1024' then for runs=1,4 do begin
   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   
   ;Do this for runs eq 1 and runs eq 2
   if runs lt 3 then output_folder = cube_folder + '/processed_left/' else output_folder = cube_folder + '/processed_right/'

   print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
   print, 'Reading in a sky_sub cube post-split...'
   obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
   print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
   restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
   print, 'Post-split flags and angles read! Changing some stuff in the cube...'
   obj_cube[where(finite(obj_cube) eq 0)] = 0.
   print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'
   ;correcting bad frames
   ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-left/' then begin
   ;if dither_folder eq '/dith1/' and night eq 2 then obj_cube[*, *, 36] = !values.f_nan
   ;if dither_folder eq '/dith2/' and night eq 2 then obj_cube[*, *, 24:38] = !values.f_nan
   ;endif
   
   ;bad_pixels = [37, 90, 91, 92, 93, 94] - 1
   ;goods = indgen(n_frames)
   
   ;remove, bad_pixels, goods
   ;print, angles
   ;obj_cube = obj_cube[*, *, goods]
   ;angles = angles[goods]
   ;print, angles
   
   ;endif
   
   ;if do_median_subtraction then for ii=0, n_frames-1 do obj_cube[*,*,ii] = obj_cube[*,*,ii] - medframe
   
   ;initialize a new cube to modify with the same data
   obj_cube2 = obj_cube
   print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   ; Divide into oblivion if we're blocking a side out
   if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*cropped_half_sz-1, *, *] *= 0.000001
   if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] *= 0.000001
   ;if do_block_left then obj_cube2[700:1023, *, *] *= 0.000001
   ;if do_block_right then obj_cube2[0:450, *, *] *= 0.000001
   if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] *= 0.000001
   if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*cropped_half_sz-1, *] *= 0.000001
   print, 'Requested blocking complete! Changing some stuff in the second cube...'
   
   ;obj_cube2[*, 0:first_center[0]-25., *] /= oblivion
   ;if do_block_top then obj_cube2[*, first_center[0]+25.:1023., *] /= oblivion
   
   ;obj_cube[*,*,0] = fshift(obj_cube[*,*,0], 300. - first_center[0], 450. - first_center[1])
   ;obj_cube2[*,*,0] = fshift(obj_cube2[*,*,0], 300. -first_center[0], 450. - first_center[1])
   
   obj_cube2[where(finite(obj_cube) eq 0)] = 0.

   obj_cube2[where(obj_cube2 lt 0)] *= 0.000000001 ;reduce negative values so CC doesn't get confused on subtracted PSF.
      print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
   half_sizex = half_cropped_sz
   half_sizey = half_cropped_sz
   
   n_frames = (size(obj_cube))[3]
   if do_cen_filter then begin
      for ii=0, n_frames - 1 do begin
         print, 'Cen_filter on frame:' + string(ii)
         print, 'Run:' + string(runs)
         obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 20., /NAN)
         obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
      endfor
   endif; do_cen_filter if
   
   print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
   ; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
   ;Commented out for now, let's try without this first.
   ;if do_cen_median_sub then medarr, obj_cube2, medframe
   
   for ii=0, n_frames-1 do begin
      
      ;Commented out for now, let's try without this first.
      ;if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
      
      ; Changed from 130(arbitrary frame?)
      if (size(obj_cube2))[3] gt 129 then corr = crosscorr(obj_cube2[*,*,130], obj_cube2[*,*,ii], pmax, /cut)
      if (size(obj_cube2))[3] lt 130 then corr = crosscorr(obj_cube2[*,*,(size(obj_cube2))[3]-1], obj_cube2[*,*,ii], pmax, /cut)
      if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
      if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
   
      print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      print, 'Run:' + string(runs)
      obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], -dx, -dy); fshift shifts frames by non-int values, which we're using to center
   endfor; ii for
   print, 'First frame shifting complete! Starting the second shift...'
   ;Basically the same as the medarr used above with obj_cube2? These don't seem to be used in this step so I'll comment them out.
   med = median(obj_cube, dim=3)
   peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
   xx = A[4] & yy = A[5]
   
   ; What is the difference between this fshift and the one above with dx and dy?
   for ii=0, n_frames-1 do begin
      obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], half_cropped_sz-xx, half_cropped_sz-yy)
      print, 'Second shift on frame: ' + string(ii)
      print, 'Run: ' + string(runs)
   endfor
   print, 'Second shifting complete! Writing centered FITS file...'
   
   ;obj_cube=obj_cube[*,fix(0.3*2*half_cropped_sz):749,*]
   ;
   ;Write our centered FITS file
   writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
   print, 'Centered FITS written for run:' + string(runs)
  
endfor; runs = 1,4 for (HII 1348)


;ALCOR SECOND SET OF 512 PIXEL ROWS (SINGLE-SIDED SX OBSERVATIONS)
if stripe eq 'second512' then for runs=1,2 do begin
   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   
   output_folder = cube_folder + '/processed_left/'

   print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
   print, 'Reading in a sky_sub cube post-split...'
   obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
   print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
   restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
   print, 'Post-split flags and angles read! Changing some stuff in the cube...'
   obj_cube[where(finite(obj_cube) eq 0)] = 0.
   print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'
   ;correcting bad frames
   ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-left/' then begin
   ;if dither_folder eq '/dith1/' and night eq 2 then obj_cube[*, *, 36] = !values.f_nan
   ;if dither_folder eq '/dith2/' and night eq 2 then obj_cube[*, *, 24:38] = !values.f_nan
   ;endif
   
   ;bad_pixels = [37, 90, 91, 92, 93, 94] - 1
   ;goods = indgen(n_frames)
   
   ;remove, bad_pixels, goods
   ;print, angles
   ;obj_cube = obj_cube[*, *, goods]
   ;angles = angles[goods]
   ;print, angles
   
   ;endif
   
   ;if do_median_subtraction then for ii=0, n_frames-1 do obj_cube[*,*,ii] = obj_cube[*,*,ii] - medframe
   
   ;initialize a new cube to modify with the same data
   obj_cube2 = obj_cube
   print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   ; Divide into oblivion if we're blocking a side out
   if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *, *] *= 0.000001
   if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] *= 0.000001
   ;if do_block_left then obj_cube2[700:1023, *, *] *= 0.000001
   ;if do_block_right then obj_cube2[0:450, *, *] *= 0.000001
   if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] *= 0.000001
   if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *] *= 0.000001
   print, 'Requested blocking complete! Changing some stuff in the second cube...'
   
   ;obj_cube2[*, 0:first_center[0]-25., *] /= oblivion
   ;if do_block_top then obj_cube2[*, first_center[0]+25.:1023., *] /= oblivion
   
   ;obj_cube[*,*,0] = fshift(obj_cube[*,*,0], 300. - first_center[0], 450. - first_center[1])
   ;obj_cube2[*,*,0] = fshift(obj_cube2[*,*,0], 300. -first_center[0], 450. - first_center[1])
   
   obj_cube2[where(finite(obj_cube) eq 0)] = 0.

   obj_cube2[where(obj_cube2 lt 0)] *= 0.000000001 ;reduce negative values so CC doesn't get confused on subtracted PSF.
      print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
   half_sizex = half_cropped_sz
   half_sizey = half_cropped_sz
   
   n_frames = (size(obj_cube))[3]
   if do_cen_filter then begin
      for ii=0, n_frames - 1 do begin
         print, 'Cen_filter on frame:' + string(ii)
         print, 'Run:' + string(runs)
         obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 20., /NAN)
         obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
      endfor
   endif; do_cen_filter if
   
   print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
   ; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
   ;Commented out for now, let's try without this first.
   ;if do_cen_median_sub then medarr, obj_cube2, medframe
   
   for ii=0, n_frames-1 do begin
      
      ;Commented out for now, let's try without this first.
      ;if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
      
      ; Changed from 130
      if (size(obj_cube2))[3] gt 129 then corr = crosscorr(obj_cube2[*,*,130], obj_cube2[*,*,ii], pmax, /cut)
      if (size(obj_cube2))[3] lt 130 then corr = crosscorr(obj_cube2[*,*,(size(obj_cube2))[3]-1], obj_cube2[*,*,ii], pmax, /cut)
      if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
      if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
   
      print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      print, 'Run:' + string(runs)
      obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], -dx, -dy); fshift shifts frames by non-int values, which we're using to center
   endfor; ii for
   print, 'First frame shifting complete! Starting the second shift...'
   ;Basically the same as the medarr used above with obj_cube2? These don't seem to be used in this step so I'll comment them out.
   med = median(obj_cube, dim=3)
   peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
   xx = A[4] & yy = A[5]
   
   ; What is the difference between this fshift and the one above with dx and dy?
   for ii=0, n_frames-1 do begin
      obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], half_cropped_sz-xx, half_cropped_sz-yy)
      print, 'Second shift on frame: ' + string(ii)
      print, 'Run: ' + string(runs)
   endfor
   print, 'Second shifting complete! Writing centered FITS file...'
   
   ;obj_cube=obj_cube[*,fix(0.3*2*half_cropped_sz):749,*]
   ;
   ;Write our centered FITS file
   writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
   print, 'Centered FITS written for run:' + string(runs)
  
endfor; runs = 1,2 for (Alcor)
print, 'Finished!'

end