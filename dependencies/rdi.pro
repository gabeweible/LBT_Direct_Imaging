pro rdi, obj, half_cropped_sz, nod, output_path, use_injection, filter, bin, bin_type,$
    do_annmode, combine_type, k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout_sx, annmode_inout_dx,$
    suffix, ct, do_cen_filter, coadd, trial=trial, fs=fs, neg_inj=neg_inj,$
    truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx,$
    pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl, refcube=refcube, ref_ang=ref_ang,$
    do_create_basis=do_create_basis, do_new_binning=do_new_binning, save_basis=save_basis,$
    rdi_extra=rdi_extra,mask_pt=mask_pt, dxklip=dxklip,$
    normal=normal, do_smooth=do_smooth, keep_number=keep_number,$
    peak_thresh=peak_thresh, stddev_thresh=stddev_thresh

    ; Optimized RDI processing procedure
    compile_opt idl2, logical_predicate
    newline = string(10B)
    
    print, newline, 'Starting RDI...'
    
    ; Default settings
    if ~keyword_set(magnify) then magnify = 0
    if ~keyword_set(do_create_basis) then do_create_basis = 0
    if ~keyword_set(do_smooth) then do_smooth=0
    
	sz = 2*half_cropped_sz
	xs = sz & ys = sz
	xhs = half_cropped_sz & yhs = half_cropped_sz
	
	; Pre-compute coordinate grid once
	xx_grid = rebin(indgen(fix(xs)), fix(xs), fix(ys))
	yy_grid = rebin(reform(indgen(fix(ys)), 1, fix(ys)), fix(xs), fix(ys))
	
	boxehs = fwhm/2.
	
    obj = strcompress(obj, /rem)
    if ~keyword_set(do_new_binning) then do_new_binning = 0
    
    ; Generate a unique suffix for output files
    if keyword_set(rdi_extra) then begin
        short_suffix = rdi_extra
    endif else begin
        short_suffix = strcompress(string(angsep, anglemax, format='("_",F0.4,"as_",F0.2,"am")'), /rem)
    endelse
    
    ; Efficiently handle RDI processing for DX aperture
    if nod eq 'dx_only' then begin
        ; Preallocate arrays for efficiency
        nods = fltarr(sz, sz, 2, /nozero)
        nods_c = fltarr(sz, sz, 2, /nozero)
        
        for runs = 3,4 do begin
        
            ; Determine dither folder
            if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
            
            ; Set parameters for DX processing
            output_folder = output_path + 'processed_right/'
            print, dither_folder
            truenorth = truenorth_dx
            annmode_inout = annmode_inout_dx
            pxscale = pxscale_dx
            
            if do_new_binning eq 1 then begin
            	
            	; Read object cube with efficient file handling
				if use_injection then begin
					obj_cube = readfits_fast(output_folder + dither_folder + obj + string(ct) + '_cube_skysub_cen_clean_inj.fits')
				endif else begin
					cube_str = strcompress(output_folder + dither_folder + obj+ '_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
						'_stddev_thresh_'+string(stddev_thresh) + '_cube_skysub_cen_clean.fits', /r)
						
					print, 'reading: ' + cube_str
					obj_cube = readfits_fast(cube_str)
				endelse
				
				; Restore parallactic angles
				restore, filename=strcompress(output_folder + dither_folder + obj + '_' +$
					string(keep_number) + '_peak_thresh_' + string(peak_thresh) + $
					'_stddev_thresh_' + string(stddev_thresh) +  '_parang_clean.sav', /r)
				
				; Efficient destriping and filtering
				processed_obj_cube = obj_cube
            
				if bin ge 2 then begin
					; Optimized binning using array reshaping
					print, 'Binning science frames...'
					n_frames = (size(processed_obj_cube, /dim))[2]
					n_bins = floor(n_frames / bin)
					
					; Reshape and reduce reference cube
					reshaped_cube = reform(processed_obj_cube[*,*,0:n_bins*bin-1], $
										   (size(processed_obj_cube, /dim))[0], $
										   (size(processed_obj_cube, /dim))[1], $
										   bin, n_bins)
					
					; Choose binning method
					if bin_type eq 'median' then begin
						processed_obj_cube = median(reshaped_cube, dimension=3, /double, /even)
					endif else begin
						processed_obj_cube = mean(reshaped_cube, dimension=3, /nan, /double)
					endelse
					
					; Average angles for binned frames
					processed_obj_angles = total(reform(angles[0:n_bins*bin-1], bin, n_bins), 1) / bin
				endif else begin
					processed_obj_angles = angles ; no binning
				endelse
				
				; by turning bad pixels to 0, then filtering, I effectively replace them with
				; nearby values.
				
				; High-pass filtering AND/OR low-pass filtering
				if filter gt 0 and do_smooth gt 0 then begin
					
					lp_npix = fix(11*do_smooth)
					if ~ODD(lp_npix) then lp_npix += 1
					lp_PSF = psf_Gaussian(NPIX=lp_npix, FWHM=[do_smooth,do_smooth], /normalize)
					
					hp_npix = fix(11*filter)
					if ~ODD(hp_npix) then hp_npix += 1
		 			hp_PSF = psf_Gaussian(NPIX=hp_npix, FWHM=[filter, filter], /normalize)
					
					print, 'Filtering science frames...'
					for iii = 0, (size(processed_obj_cube, /dim))[2]-1 do begin
					
						frame_iii = processed_obj_cube[*,*,iii]
						; Replace NaNs with zeros
						nan_indices = where(finite(frame_iii) ne 1, nan_count)
						if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
					
						; high-pass
						processed_obj_cube[*,*,iii] -= convolve(frame_iii, hp_PSF)
						
						frame_iii = processed_obj_cube[*,*,iii]
						; Replace NaNs with zeros
						nan_indices = where(finite(frame_iii) ne 1, nan_count)
						if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
						
						; low-pass
						processed_obj_cube[*,*,iii] = convolve(frame_iii, lp_PSF)
						
					endfor
				endif else begin; only one filtering or the other.
				
					if filter gt 0 then begin; just high-pass filtering
						hp_npix = fix(11*filter)
						if ~ODD(hp_npix) then hp_npix += 1
		 				hp_PSF = psf_Gaussian(NPIX=hp_npix, FWHM=[filter, filter], /normalize)
						
						print, 'Filtering science frames...'
						for iii = 0, (size(processed_obj_cube, /dim))[2]-1 do begin
							frame_iii = processed_obj_cube[*,*,iii]
							
							; Replace NaNs with zeros
							nan_indices = where(finite(frame_iii) ne 1, nan_count)
							if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
							
							; high-pass
							processed_obj_cube[*,*,iii] -= convolve(frame_iii, hp_PSF)
						endfor
					endif
					
					if do_smooth gt 0 then begin; just low-pass filtering
						lp_npix = fix(11*do_smooth)
						if ~ODD(lp_npix) then lp_npix += 1
						lp_PSF = psf_Gaussian(NPIX=lp_npix, FWHM=[do_smooth,do_smooth], /normalize)
					
						print, 'Filtering science frames...'
						for iii = 0, (size(processed_obj_cube, /dim))[2]-1 do begin
						frame_iii = processed_obj_cube[*,*,iii]
						
						; Replace NaNs with zeros
						nan_indices = where(finite(frame_iii) ne 1, nan_count)
						if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
						
						; low-pass
						processed_obj_cube[*,*,iii] = convolve(frame_iii, lp_PSF)
						endfor
					endif
				endelse
				
				; Normalize all of the frames so that they each have a max value of 1 (should help with residuals around the star)
				if keyword_set(normal) then begin
					print, 'Normalizing frames...'
				   if normal eq 1 then for iv=0, n_frames - 1 do begin; only normalize to the smoothed peak		
						; try a more robust normalization?
						
						frame_iv = processed_obj_cube[*,*,iv]
						
						; grab pixel values within a radius of fwhm/2.
						center_pixels = []
						for dx=-fwhm*(0.75),fwhm*(0.75) do begin
						for dy=-fwhm*(0.75),fwhm*(0.75) do begin
							if sqrt( dx^2. + dy^2. ) lt fwhm*0.75 then center_pixels = [[center_pixels],$
								frame_iv[xhs+dx, yhs+dy]]
						endfor
						endfor
						
						processed_obj_cube[*,*,iv] /= (median(center_pixels, /even, /double) / 100.)
							
				   endfor; double-check on normalization
				endif
				
				; Save preprocessed obj cube
				print, 'deleting old objs and saving new objs...'
				objs_filename = strcompress(output_path+obj+short_suffix+'_processed_objs_'+'run'+string(runs)+'.sav',/r)
				file_delete, objs_filename, /allow_nonexistent
				save, processed_obj_cube, processed_obj_angles, $
					  filename=objs_filename  
					  
			endif else begin; read-in object cube (no new binning)
			
				; Restore existing processed obj cube
				processed_objs_name = strcompress(output_path+obj+short_suffix+'_processed_objs_'+'run'+string(runs)+'.sav',/r)
				print, 'restoring: ' + processed_objs_name
				restore, filename=processed_objs_name

			endelse
			
			; Optimized reference cube preprocessing
			if ~keyword_set(dxklip) then begin; bin and/or read references
			
				if do_new_binning eq 1 then begin
					processed_refcube = bin_refcube(refcube, bin, ref_ang, filter,$
						do_smooth, output_path, obj, short_suffix, processed_ref_angles,$
						bin_type=bin_type,half_cropped_sz=half_cropped_sz, normal=normal, fwhm=fwhm)
					; processed_ref_angles is set after running
				endif else begin
					restore_filename = output_path+obj+short_suffix+'_processed_refs.sav'
					print, 'restoring ' + restore_filename
					restore, filename=restore_filename
				endelse
				
			endif else begin
			
				; if doing dxklip, then the objects are the references
				if dxklip eq 1 then begin
				 	print, 'Doing true ADI-KLIP w/obj as references'
					processed_refcube = processed_obj_cube
					processed_ref_angles = processed_obj_angles
				endif
			endelse
            
            ; Efficient KLIP processing
            klip_cube = processed_obj_cube ; Copy for KLIP subtraction
            numbers_of_frames = (size(processed_obj_cube))[3]
            
            ; get numbers of frames for each nod
			if runs eq 3 then begin
				numbers_of_frames = (size(processed_obj_cube))[3]
			endif else begin
				numbers_of_frames = [[numbers_of_frames], (size(processed_obj_cube))[3]]
			endelse
            
            replicate_inplace, klip_cube, 0.; will be replaced, anyway

            for ii = 0, (size(processed_obj_cube))[3]-1 do begin
            
                ; Use existing rdiklip2 function with optimized parameters
                klip_frame_ii = rdiklip2(processed_obj_cube, processed_refcube,$
                	k_klip, target=ii, posang=processed_obj_angles,$
                	wl=wl, diam=8.4,pixelscale=pxscale, angsep=angsep,$
                	anglemax=anglemax,obj=obj, nrings=nrings, wr=wr, n_ang=n_ang,$
                	annmode_inout=annmode_inout, create_basis=do_create_basis,$
                	ref_angles=processed_ref_angles, filter=filter, bin=bin, $
                    FWHM=fwhm, save_basis=save_basis, short_suffix=short_suffix, $
                    output_path=output_path)
                    
                if ii eq 0 then writefits, '~/Desktop/rot_klip_frame_ii_test.fits', klip_frame_ii
                
                ; Rotate frame to North-up
                klip_cube[*,*,ii] = rot(klip_frame_ii,$
                	-processed_obj_angles[ii]-truenorth, 1.0, half_cropped_sz,$
                	half_cropped_sz,cubic=-1.0, missing=median(klip_frame_ii, /double, /even),$
                	/pivot)
            endfor
            
            ;writefits, '~/Desktop/processed_obj_cube.fits', processed_obj_cube
            ;writefits, '~/Desktop/klip_cube_test.fits', klip_cube
            
            ; Combine frames based on specified method
            if combine_type eq 'median' then medarr, klip_cube, medframe
            if combine_type eq 'mean' then medframe = mean(klip_cube, dimension=3, /double, /nan)
            if combine_type eq 'nwadi' then medframe = nw_ang_comb(klip_cube, processed_obj_angles)
            
			print, 'Filtering klipframe...'
		
			; High-pass filtering AND/OR low-pass filtering
			if filter gt 0 and do_smooth gt 0 then begin
			
				frame_iii = medframe
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
			
				; high-pass
				medframe -= convolve(frame_iii, hp_PSF)
				
				frame_iii = medframe
				; Replace NaNs with zeros
				nan_indices = where(finite(frame_iii) ne 1, nan_count)
				if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
				; low-pass
				medframe = convolve(frame_iii, lp_PSF)
		
			endif else begin; only one filtering or the other.
			
				if filter gt 0 then begin; just high-pass filtering
					frame_iii = medframe
					
					; Replace NaNs with zeros
					nan_indices = where(finite(frame_iii) ne 1, nan_count)
					if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
					
					; high-pass
					medframe -= convolve(frame_iii, hp_PSF)
				endif
				
				if do_smooth gt 0 then begin; just low-pass filtering
					frame_iii = medframe
					
					; Replace NaNs with zeros
					nan_indices = where(finite(frame_iii) ne 1, nan_count)
					if nan_count gt 0 then frame_iii[nan_indices] = median(frame_iii, /double, /even)
					
					; low-pass
					medframe = convolve(frame_iii, lp_PSF)
				endif
			endelse
		
			; Store nods
			nods[*,*,runs-3] = medframe
        endfor; nods for
        
        ; FINAL STUFF AFTER EACH NOD
        
        ; right combination only
		final = (numbers_of_frames[0]*nods[*,*,0] + numbers_of_frames[1]*nods[*,*,1])/total(numbers_of_frames)
	
        ; Write output files
        super_suffix = strcompress(output_path + 'combined/' + obj +'_k_klip_'+string(sigfig(k_klip,3)) +$
	   		'_comb_type_'+combine_type+'_angsep_'+ string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	   		string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_filter_' + string(filter) + '_bin_' + string(sigfig(bin,2)) + '_type_' + bin_type +  '_neg_inj_' + string(neg_inj), /r)
        
        if keyword_set(trial) then super_suffix += '_trial_' + string(trial, format='(F0.4)')
        
        ; Write combined and individual nod results
        writefits, strcompress(super_suffix + '_combined.fits', /rem), final
        
        writefits, strcompress(super_suffix + '_nod3.fits', /rem), nods[*,*,0]
        writefits, strcompress(super_suffix + '_nod4.fits', /rem), nods[*,*,1]
        
        ; also, convolve the combined-over-nods frame after the fact and write this.
        ; convolving by fwhm/2. should be the same as what snr_map.pro does.
        PSF = psf_Gaussian(NPIX=73, FWHM=[fwhm/2.,fwhm/2.], /normalize)
	    
	    final[~finite(final)] = median(final, /double, /even); get rid of NaN before convolving
		right_klip_conv = convolve(final, PSF)
		
		; Write combined and individual nod results
        writefits, strcompress(super_suffix + '_combined_conv.fits', /rem), right_klip_conv

    endif ; 'dx_only' (only thing implemented right now - for Alcor)
end