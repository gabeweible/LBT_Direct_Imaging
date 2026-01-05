pro center_old, obj_name, stripe, cube_folder, do_block_right, do_block_left,$
	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub, half_cropped_sz,$
	do_smooth=do_smooth, post_vip=post_vip, align_radius_thresh=align_radius_thresh,$
	align_inner_radius=align_inner_radius, cenframewidth=cenframewidth,$
	inner_rad_cen=inner_rad_cen, outer_rad_cen=outer_rad_cen, rot_center=rot_center
; 'HII1348', '~/OneDrive/Research/HII1348/testing', 1, 1, 1, 1, 1, 0 for current run (let's try cen_median_sub later)
compile_opt idl2

new_filter=27.0
censteps=21	;steps in x and y - i.e. total steps will be square of this
censquare=0.05	;test square size in pixels (please enter a float value)
; censquare=2. means search +- 1 pixel square
radius = half_cropped_sz

sz = 2*half_cropped_sz

; I think this is the second time that the runs matter (we are now post-split)
; HII 1348 CENTER 1024 PIXEL ROWS (DOUBLE-SIDED)
if stripe eq 'center1024' then for runs=1,4 do begin
   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   
   ;Do this for runs eq 1 and runs eq 2
   if runs lt 3 then output_folder = cube_folder + '/processed_left/' else output_folder = cube_folder + '/processed_right/'

   print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
   print, 'Reading in a sky_sub cube post-split...'
   
   
    if keyword_set(post_vip) then begin
        obj_cube = readfits_fast(output_folder + dither_folder + obj_name + '_cube_skysub_cen.fits')
    endif else begin
         obj_cube = readfits_fast(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
    endelse
   
   print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
    if keyword_set(post_vip) then begin
        angles = readfits_fast(output_folder + dither_folder + obj_name + '_angles.fits')
    endif else begin
        restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
    endelse
   print, 'Post-split angles read! nans to 0...'
   
   obj_cube[where(finite(obj_cube) eq 0)] = 0.
   print, 'Nans changed to zero! Creating the obj_cube copy...'
   
   ;initialize a new cube to modify with the same data
   obj_cube2 = obj_cube
   print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   ; Divide into oblivion if we're blocking a side out
   if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *, *] = 0.
   if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] *= 0.
   if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] = 0.
   if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *] = 0.
   print, 'Requested blocking complete! Changing some stuff in the second cube...'
   
   obj_cube2[where(finite(obj_cube) eq 0)] = 0.

   obj_cube2[where(obj_cube2 lt 0)] = 0.;reduce negative values so CC doesn't get confused on subtracted PSF.
      print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
   half_sizex = half_cropped_sz
   half_sizey = half_cropped_sz
   
   n_frames = (size(obj_cube))[3]
   if do_cen_filter then begin
      for ii=0, n_frames - 1 do begin
         print, 'Cen_filter on frame:' + string(ii)
         print, 'Run:' + string(runs)
         
        ;; high- and low-pass filter the frame so we can better find its peak
        ; obj_cube2[*,*,ii] -= filter_image(temporary(obj_cube2[*,*,ii]), fwhm=new_filter,$
        ; 	PSF=PSFnew, /ALL_PIXELS)
         
         obj_cube2[*,*,ii] = filter_image(temporary(obj_cube2[*,*,ii]), fwhm=2.0, PSF=PSF1,$
         	/ALL_PIXELS)
      endfor
   endif; do_cen_filter if
   
   ;focus only on inner part
    for xx=0,2.0*half_cropped_sz-1 do for yy=0,2.0*half_cropped_sz-1 do if sqrt((xx-half_cropped_sz)^2.0+(yy-half_cropped_sz)^2.0) gt align_radius_thresh then obj_cube2[xx,yy,*]=0.


    ;block saturated part
    for xx=0,2.0*half_cropped_sz-1 do for yy=0,2.0*half_cropped_sz-1 do if sqrt((xx-half_cropped_sz)^2.0+(yy-half_cropped_sz)^2.0) lt align_inner_radius then obj_cube2[xx,yy,*]=0.

   
   print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
   
   obj_cube2 = obj_cube2[0:2*half_cropped_sz-1, 0:2*half_cropped_sz-1, *]
    refframe = median(obj_cube2, dim=3, /double, /even)
    for ii=0, n_frames-1 do begin
      
        ; Changed from 130(arbitrary frame?)
        corr = crosscorr(refframe, obj_cube2[*,*,ii], pmax, /cut)
      
        ; sub-pixel shifts required
        if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
        if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
   
        print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
         print, 'Run:' + string(runs)
      
        obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), -dx, -dy, cubic=-1.0); shift_sub shifts frames by non-int values, which we're using to center
      
        ; low-pass filter. Should help cubic convolution.
	  	;    if do_smooth eq 1 then begin
	  		   ; frame_ii = obj_cube[*,*,ii]
	  		  ;  nan_indices = where(finite(frame_ii) ne 1, nan_count)
			 ;   if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
			;    obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		  ;  endif
      
   endfor; ii for
   
   print, 'First frame shifting complete! Starting the second shift...'
   
   
    if not keyword_set(rot_center) then begin
        ; absolute registration based on 2D Gaussian fit to median frame (before was co-registration.)
       med = median(obj_cube, dim=3, /even, /double)
       peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
       xx = A[4] & yy = A[5]
       
       for ii=0, n_frames-1 do begin
          
          obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), half_cropped_sz-xx,$
            half_cropped_sz-yy, cubic=-1.0)
            
             ; low-pass filter. Should help cubic convolution.
           ; if do_smooth eq 1 then begin
              ;  frame_ii = obj_cube[*,*,ii]
             ;   nan_indices = where(finite(frame_ii) ne 1, nan_count)
            ;    if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
           ;     obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
          ;  endif
            
          print, 'Second shift on frame: ' + string(ii)
          print, 'Run: ' + string(runs)
       endfor
       
        
       print, 'Second shifting complete!'
       
       ;Write our centered FITS file
       writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen_IDL.fits', obj_cube
       print, 'Centered FITS written for run:' + string(runs)
    endif else begin; rotation centering

        obj_cubei=obj_cube
        obj_cubei=obj_cubei[radius-cenframewidth/2.:radius+cenframewidth/2.,radius-cenframewidth/2.:radius+cenframewidth/2.,*]
        writefits, '~/Desktop/obj_cubi.fits', obj_cubei
        
        bestresiduals=99e99
        for kk=0, censteps-1 do begin
            for jj=0, censteps-1 do begin
                print, 'On position ', (kk)*censteps + jj+1 , ' out of ', censteps*censteps
                
                obj_cube2=obj_cubei;obj_cubes
                
                ;shift cube frames
                shiftx=float(kk)*censquare/float(censteps-1) - censquare/2.
                shifty=float(jj)*censquare/float(censteps-1) - censquare/2.
                
                print, 'Testing shift x', shiftx, ' shift y', shifty
                
                for ll=0, (size(obj_cube2))[3]-1 do obj_cube2[*,*,ll]=shift_sub(obj_cube2[*,*,ll],shiftx,shifty,$
                                                                        cubic=-1.0)
                
                ;take median
                medarr, obj_cube2, medframe
                
                ;derotate images by 180 degrees, then subtract the pupil-stabilized median
                for ii=0, (size(obj_cube2))[3]-1 do begin
                    obj_cube2[*,*,ii]=rot(obj_cube2[*,*,ii],180.,$
   		                1.0, cenframewidth/2., cenframewidth/2.,missing=0, /pivot,$
   		                cubic=-1.0) - medframe
                endfor
                
                ;block outer regions 
                for xx=0,(size(obj_cube2))[1]-1 do begin
                    for yy=0,(size(obj_cube2))[2]-1 do begin
                        xd=abs(xx-cenframewidth/2.0)
                        yd=abs(yy-cenframewidth/2.0)
                        if sqrt(xd*xd + yd*yd) gt outer_rad_cen then obj_cube2[xx,yy,*]=0.
                    endfor
                endfor
                
                ;block inner regions 
                for xx=0,(size(obj_cube2))[1]-1  do begin
                    for yy=0,(size(obj_cube2))[2]-1  do begin
                        xd=abs(xx-cenframewidth/2.0)
                        yd=abs(yy-cenframewidth/2.0)
                        if sqrt(xd*xd + yd*yd) lt inner_rad_cen then obj_cube2[xx,yy,*]=0.
                    endfor
                endfor
                
                obj_cube2[where(finite(obj_cube2) eq 0)]=0.
                residuals=total(abs(obj_cube2))
                print, residuals
                
                if residuals le bestresiduals then begin
                    
                    print, 'Found better residual!'
                
                    bestresiduals=residuals
                    bestimage=median(obj_cube2,dim=3)
                    bestshiftx=shiftx
                    bestshifty=shifty
                    
                endif
            endfor	;kk censteps for
        endfor	;jj censteps for
        
        print, 'best x shift: ', bestshiftx
        print, 'best y shift: ', bestshifty
        
        if bestshiftx eq censquare/2. or bestshifty eq censquare/2. then begin
            print, 'Precise centering step maxed out! Increase box size and re-run!!!'
        endif
        
        for ll=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ll]=shift_sub(obj_cube[*,*,ll],bestshiftx,bestshifty,cubic=-1.0)
        
         
       print, 'Second shifting complete!'
       
       ;Write our centered FITS file
       writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_rot_cen_IDL.fits', obj_cube
       print, 'Centered FITS written for run:' + string(runs)
    endelse
   
endfor; runs = 1,4 for (HII 1348, TYC 5709)


;ALCOR SECOND SET OF 512 PIXEL ROWS (SINGLE-SIDED SX OBSERVATIONS)
if stripe eq 'second512' then for runs=3,4 do begin
   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   
   output_folder = cube_folder + '/processed_right/'

   print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
   print, 'Reading in a sky_sub cube post-split...'
   obj_cube = readfits_fast(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
   print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
   restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
   print, 'Post-split flags and angles read! Changing some stuff in the cube...'
   obj_cube[where(finite(obj_cube) eq 0)] = 0.
   print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'

   ;initialize a new cube to modify with the same data
   obj_cube2 = obj_cube
   print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   ; set to zero if we're blocking a side out
   if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *, *] /= 1000000000000000000.
   if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] /= 1000000000000000000.
   if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] /= 1000000000000000000.
   if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *] /= 1000000000000000000.
   print, 'Requested blocking complete! Changing some stuff in the second cube...'

   ;obj_cube2[where(obj_cube2 lt 0)] =0. ;reduce negative values so CC doesn't get confused on subtracted PSF.
      print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
   half_sizex = half_cropped_sz
   half_sizey = half_cropped_sz
  
	obj_cube2[where(obj_cube2 lt 0)] /= 1000000000000000000. ;reduce negative values so CC doesn't get confused on subtracted PSF.
   
   n_frames = (size(obj_cube))[3]
   if do_cen_filter then begin
   		print, 'starting cen_filter...'
      for ii=0, n_frames - 1 do begin
         print, 'Cen_filter on frame:' + string(ii)
         print, 'Run:' + string(runs)
		; high- and low-pass filter the frame so we can better find its peak
         obj_cube2[*,*,ii] -= filter_image(temporary(obj_cube2[*,*,ii]), fwhm=new_filter,$
         	PSF=PSFnew, /ALL_PIXELS)
         
         obj_cube2[*,*,ii] = filter_image(temporary(obj_cube2[*,*,ii]), fwhm=1.,$
         	PSF=PSF1, /ALL_PIXELS)
      endfor
   endif; do_cen_filter if
   
   ; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
    if do_cen_median_sub then medarr, obj_cube2, medframe
	
   print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
   
    ref_frame = obj_cube2[*,*,130]
    
	; first shift: use a single reference frame to co-register images (i.e., to each other,
	; referenced to frame index 130, in this case)
   for ii=0, n_frames-1 do begin
      
      frame_ii = obj_cube2[*,*,ii]
      
      ; make sure that this frame (index 131 in DS9) looks good for both nodding positions...
      ; this is crosscorrelating the filtered frame with another filtered frame, but this won't actually be shifted
      corr = crosscorr(ref_frame, frame_ii, pmax, /cut)
      
		if pmax[0] ge 0. then dx = half_sizex-pmax[0] else dx = -half_sizex+abs(pmax[0])
		if pmax[1] ge 0. then dy = half_sizey-pmax[1] else dy = -half_sizey+abs(pmax[1])
   
      print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      print, 'Run:' + string(runs)
      
      ; shift to peak cross-correlation location with ref_frame
      obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), -dx, -dy, cubic=-1.0); shift_sub shifts frames by non-int values, which we're using to center
      
       ; low-pass filter. Should help cubic convolution.
	  	if do_smooth eq 1 then begin
	  		frame_ii = obj_cube[*,*,ii]
	  		nan_indices = where(finite(frame_ii) ne 1, nan_count)
			if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
			obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		endif
   endfor; ii for
   
   print, 'First frame shifting complete! Starting the second shift...'
   
   ; take the median of the co-registered frames (still stacked at the native ref_frame location)
   med = median(obj_cube, dim=3, /even, /double)
   
   ; find the peak of where the co-registration leaves the median
   peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
   xx = A[4] & yy = A[5]
   
   ; Second shift: move each frame by a common shift (xx, yy) in order to bring them
   ; to the center of the image (centered at [half_cropped_sz, half_cropped_sz])
   	for ii=0, n_frames-1 do begin
      	obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), half_cropped_sz-xx,$
      		half_cropped_sz-yy, cubic=-1.0)
      		
      	 ; low-pass filter. Should help cubic convolution.
	  	if do_smooth eq 1 then begin
	  		frame_ii = obj_cube[*,*,ii]
	  		nan_indices = where(finite(frame_ii) ne 1, nan_count)
			if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
			obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		endif
		
      	print, 'Second shift on frame: ' + string(ii)
      	print, 'Run: ' + string(runs)
   	endfor
   	print, 'Second shifting complete!'
   
   	; third shift - now that the images are reasonably co-registered and common-centered
   	; to the frame center, we can use a new 'med' to further co-register
   	
   	; remove NaNs (again)
   	obj_cube[where(finite(obj_cube) eq 0)] = 0.
   	;initialize a new cube to modify with the same data
   	obj_cube2 = obj_cube
   	print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   	; set to zero if we're blocking a side out
   	if do_block_right then obj_cube2[fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *, *] /= 1000000000000000000.
   	if do_block_left then obj_cube2[0:fix(0.3*2*half_cropped_sz), *, *] /= 1000000000000000000.
   	if do_block_bottom then obj_cube2[*, 0:fix(0.3*2*half_cropped_sz), *] /= 1000000000000000000.
   	if do_block_top then obj_cube2[*,fix(0.7*2*half_cropped_sz):2*half_cropped_sz-1, *] /= 1000000000000000000.
  	print, 'Requested blocking complete! Changing some stuff in the second cube...'

	obj_cube2[where(obj_cube2 lt 0)] /= 1000000000000000000. ;reduce negative values so CC doesn't get confused on subtracted PSF.

   	if do_cen_filter then begin
   		print, 'starting cen_filter for third shift...'
      	for ii=0, n_frames - 1 do begin
         	print, 'Cen_filter on frame:' + string(ii)
         	print, 'Run:' + string(runs)
			; high- and low-pass filter the frame so we can better find its peak
         	obj_cube2[*,*,ii] -= filter_image(temporary(obj_cube2[*,*,ii]), fwhm=new_filter,$
         		PSF=PSFnew, /ALL_PIXELS)
         
         	obj_cube2[*,*,ii] = filter_image(temporary(obj_cube2[*,*,ii]), fwhm=1., PSF=PSF1,$
         		/ALL_PIXELS)
      	endfor
   	endif; do_cen_filter if
   	
   	; take the filtered reference to co-register (again)
   	ref_frame = median(obj_cube2, dim=3, /even, /double)
   	for ii=0, n_frames-1 do begin
      
		frame_ii = obj_cube2[*,*,ii]
      
      	corr = crosscorr(ref_frame, frame_ii, pmax, /cut)
     	if pmax[0] ge 0. then dx = half_sizex-pmax[0] else dx = -half_sizex+abs(pmax[0])
		if pmax[1] ge 0. then dy = half_sizey-pmax[1] else dy = -half_sizey+abs(pmax[1])
   
      	print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      	print, 'Run:' + string(runs)
      
      	; shift to peak cross-correlation location with new median-ref_frame
      	obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), -dx, -dy, cubic=-1.0); shift_sub shifts frames by non-int values, which we're using to center
      	
      	 ; low-pass filter. Should help cubic convolution.
	  	if do_smooth eq 1 then begin
	  		frame_ii = obj_cube[*,*,ii]
	  		nan_indices = where(finite(frame_ii) ne 1, nan_count)
			if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
			obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		endif
   endfor; ii for
   print, 'Third frame shifting complete! Starting the fourth shift...'
   
   ; take the median of the co-registered frames (should be stacked very close to zero
   ; from the third shift to a centered, filtered medframe (2nd co-reg))
   med = median(obj_cube, dim=3, /even, /double)
   
   ; find the peak of where the median now is, after all other shifts
   peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
   xx = A[4] & yy = A[5]
   
   ; Fourth shift: move each frame by a common shift (xx, yy) in order to bring them
   ; to the center of the image (centered at [half_cropped_sz, half_cropped_sz])
   	for ii=0, n_frames-1 do begin
      	obj_cube[*,*,ii] = shift_sub(temporary(obj_cube[*,*,ii]), half_cropped_sz-xx,$
      		half_cropped_sz-yy, cubic=-1.0)
      		
      	 ; low-pass filter. Should help cubic convolution.
	  	if do_smooth eq 1 then begin
	  		frame_ii = obj_cube[*,*,ii]
	  		nan_indices = where(finite(frame_ii) ne 1, nan_count)
			if nan_count gt 0 then frame_ii[nan_indices] = median(frame_ii, /double, /even)
			obj_cube[*,*,ii] = convolve(frame_ii, lp_PSF)
		endif
      	print, 'Fourth shift on frame: ' + string(ii)
      	print, 'Run: ' + string(runs)
   	endfor
   	print, 'Fourth shifting complete!'
   
   ;Write our centered FITS file
   writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
   print, 'Centered FITS written for run:' + string(runs)
  
endfor; runs = 1,2 for (Alcor)
print, 'Finished!'

end