pro klip, obj, half_cropped_sz, nod, cube_folder, use_injection, do_destripe, filter, bin, bin_type,$
	do_hyper, do_annmode, combine_type, klip_fraction, start_frame, end_frame, fill,$
	k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout_sx, annmode_inout_dx,$
	suffix, ct, do_cen_filter, coadd, trial=trial, fs=fs, neg_inj=neg_inj,$
	truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, pxscale_sx=pxscale_sx,$
	pxscale_dx=pxscale_dx, magnify=magnify, fwhm=fwhm, wl=wl,rdi=rdi,refcube=refcube,ref_angles=ref_angles,$
	basis=basis,new_binning=new_binning
   
newline = string(10B)
compile_opt IDL2
filter=filter
bin=bin
basis=basis
new_binning=new_binning
refcube=refcube
ref_angles=ref_angles

; default to not doing rdi and not passing a basis
if not keyword_set(basis) then basis=0
if not keyword_set(rdi) then rdi=0

sz=2*half_cropped_sz
width = fwhm

if rdi eq 1 and new_binning eq 1 then begin
	;if do_destripe eq 1 then begin
	;	print, 'destriping 90 degrees...'
	;	for ii=0, (size(refcube))[3]-1 do refcube[*,*,ii]=destripe(refcube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;	print, 'destriping 0 degrees...'
	;	for ii=0, (size(refcube))[3]-1 do refcube[*,*,ii]=destripe(refcube[*,*,ii],0.,clip_level=0.0,/nodisp)
	;endif; destripe w/rdi
	
	if filter gt 1 then for iii=0,(size(refcube))[3]-1 do refcube[*,*,iii]=refcube[*,*,iii]-$
		smooth(refcube[*,*,iii],filter)
		
	if bin gt 1 and bin_type eq 'median' then begin
	   st=1
	   for ii=0.,(size(refcube))[3]-1. do begin
		  print, fix(ii+1.) mod fix(bin)
		  
		  if fix(ii+1.) mod fix(bin) eq 1 then begin; start over on 1
				binned = refcube[*,*,ii]
				binned_angle = ref_angles[ii]
		  endif else begin; concatenate, or add
				binned = [ [[binned]], [[ refcube[*,*,ii] ]] ]
				binned_angle = binned_angle + ref_angles[ii]
		  
		  
		  if fix(ii+1.) mod fix(bin) eq 0 then begin; bin on 0
		  
			 print,'Binning RDI library frames...'
			 
			  if st eq 1 then begin ; first binning (set at top of loop)
				binned_cube = median(binned, dim=3, /even, /double)
				binned_angles = [binned_angle/bin]
				st=0 
			  endif else begin; st eq 1 if, else is for future binning
				new_binned = median(binned,dim=3,/even,/double)
				new_binned_angle = binned_angle / bin
				binned_cube=[[[binned_cube]],[[new_binned]]]; add to cube
				binned_angles=[[binned_angles],[new_binned_angle]]; add to angles
			  endelse; st neq 1
			  
		  endif; mod 0 if
		  endelse; not starting over if
	   endfor   ; cube loop for

	   print, size(binned_cube)
	   print, size(binned_angles)

	   refcube_rdi = binned_cube
	   ref_angles_rdi = binned_angles
	   ;print, newline, ref_angles_rdi, newline
	endif; bin gt 1

		;if bin gt 1 and bin_type eq 'mean' then begin
		 ;  st=1
		  ; for ii=0.,(size(refcube))[3]-1. do begin
			;  print, fix(ii+1.) mod fix(bin)
			 ; if fix(ii+1.) mod fix(bin) eq 1 then binned=refcube[*,*,ii] else binned=binned +refcube[*,*,ii] 
			;	if fix(ii+1.) mod fix(bin) eq 0 then begin
			;		 print,'Binning RDI library frames...'
			;		 
			;		if st eq 1 then begin 
			;			 binned_cube = binned / float(bin)
			;			 st=0 
			;		endif else begin
			;			binned_cube=[[[binned_cube]],[[binned/float(bin)]]]
			;		endelse
			;			  
			;		print, size(binned_cube)
			;	endif; fix mod 0 if  
		   ;endfor; cube loop for
		   
		   ;print, size(binned_cube)
		   ;st=1
		   
		   ;for ii=0,(size(refcube))[3]-1. do begin
			;  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=ref_angles[ii] else binned_angle=binned_angle+$
			;	 ref_angles[ii]
			 ; if fix(ii+1.) mod fix(bin) eq 0 then begin
			;	 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
			;		[[[binned_angles]],[[binned_angle/bin]]]
			 ; endif; fix mod eq 0 if
		   ;endfor; cube loop for
   
		  ; refcube_rdi=binned_cube
		  ; ref_angles_rdi=binned_angles
		;endif; bin gt 1 if
		
   lbc=0
	bads=fltarr( (size(refcube_rdi))[3] )
	bads[*]=0.

	refcube_rdi[where( finite( refcube_rdi ) eq 0 )]=0.

	for ii=0, (size(refcube_rdi))[3]-1 do begin
	   if total(refcube_rdi[*,*,ii]) eq 0 then begin
		  print, ii
		remove, ii-lbc, ref_angles_rdi
		  bads[ii]=1.
		  lbc=lbc+1
	   endif
	endfor
	refcube_rdi=refcube_rdi[*,*,where(bads eq 0.)]
	
	save, refcube_rdi, ref_angles_rdi, filename='~/Desktop/Alcor_rdi_refs.sav'
		
endif else if rdi eq 1 and new_binning eq 0 then restore, '~/Desktop/Alcor_rdi_refs.sav'

print, size(binned_cube)
print, size(binned_angles)

if not keyword_set(magnify) then magnify = 0; Default to no magnification

if nod eq 'total' then begin
	min_pxscale = min([pxscale_sx, pxscale_dx])
	for runs=1,4 do begin
		if rdi eq 1 then refcube_rdi = refcube
	   ; Do this for runs eq 1 and runs eq 3
	   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
	   ;Do this for runs eq 1 and runs eq 2
	   if runs lt 3 then begin
		  output_folder = cube_folder + 'processed_left/'
		  truenorth = truenorth_sx
		  annmode_inout = annmode_inout_sx
		  pxscale = pxscale_sx
	   endif else begin; runs lt 3 if
	   		output_folder = cube_folder + 'processed_right/'
	   		truenorth = truenorth_dx
	   		annmode_inout = annmode_inout_dx
	   		pxscale = pxscale_dx
	   endelse; runs gt 2

		obj=strcompress(obj,/rem)
		
		if use_injection then obj_cube = readfits(output_folder + dither_folder + obj + string(ct) +$
		   '_cube_skysub_cen_clean_inj.fits') else obj_cube = readfits(output_folder + dither_folder +$
		   obj + string(ct) +  '_cube_skysub_cen_clean.fits')

		; destriping
		;if do_destripe eq 1 then begin
		;	print, 'destriping 90 degrees...'
		;	for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
		;	print, 'destriping 0 degrees...'
		;	for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],0.,clip_level=0.0,/nodisp)
		;endif
		
		; high-pass filter
		if filter gt 1 then for iii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,iii]=obj_cube[*,*,iii]-$
		   smooth(obj_cube[*,*,iii],filter)

		restore,filename=output_folder + dither_folder + obj + string(ct) +  '_parang_clean.sav'
		
		; restrict (unbinned cube)
		if klip_fraction eq 1 then begin
		   obj_cube=obj_cube[*,*,start_frame:end_frame]
		   angles=angles[start_frame:end_frame]
		endif
		
		if bin gt 1 and bin_type eq 'median' then begin
		   st=1
		   for ii=0.,(size(obj_cube))[3]-1. do begin
			  print, fix(ii+1.) mod fix(bin)
			  if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=[ [[binned]],$
				 [[ obj_cube[*,*,ii] ]] ]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
			  
				 print,'Binning left frames...'
				 
				 if st eq 1 then begin 
				 	binned=median(binned,dim=3,/even,/double)
				 	binned_cube = binned 
				 	st=0 
				 endif else begin; st eq 1 if
					binned=median(binned,dim=3,/even,/double)
				  	binned_cube=[[[binned_cube]],[[binned]]]
				  endelse; st neq 1
				 print, size(binned_cube)
			  endif; mod 0 if
		   endfor   ; cube loop for
   
		   print, size(binned_cube)
		   st=1
		   
		   for ii=0,(size(obj_cube))[3]-1. do begin
			  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+$
				 angles[ii]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
				 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
					[[[binned_angles]],[[binned_angle/bin]]]
			  endif; fix mod eq 0 if
		   endfor; cube loop for
   
		   obj_cube=binned_cube
		   angles=binned_angles
		endif; bin gt 1

		if bin gt 1 and bin_type eq 'mean' then begin
		   st=1
		   for ii=0.,(size(obj_cube))[3]-1. do begin
			  print, fix(ii+1.) mod fix(bin)
			  if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=binned + obj_cube[*,*,ii] 
				if fix(ii+1.) mod fix(bin) eq 0 then begin
					 print,'Binning left frames...'
					 
					if st eq 1 then begin 
						 binned_cube = binned / float(bin)
						 st=0 
					endif else begin
						binned_cube=[[[binned_cube]],[[binned/float(bin)]]]
					endelse
						  
					print, size(binned_cube)
				endif; fix mod 0 if  
		   endfor; cube loop for
		   
		   print, size(binned_cube)
		   st=1
		   
		   for ii=0,(size(obj_cube))[3]-1. do begin
			  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+$
				 angles[ii]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
				 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
					[[[binned_angles]],[[binned_angle/bin]]]
			  endif
		   endfor; cube loop for
   
		   obj_cube=binned_cube
		   angles=binned_angles
		endif; bin gt 1 if

		lbc=0
		bads=fltarr( (size(obj_cube))[3] )
		bads[*]=0.
		
		obj_cube[where( finite( obj_cube ) eq 0 )]=0.
		
		for ii=0, (size(obj_cube))[3]-1 do begin
		   if total(obj_cube[*,*,ii]) eq 0 then begin
			  print, ii
			  remove, ii-lbc, angles
			  bads[ii]=1.
			  lbc=lbc+1
		   endif; total of frame is 0
		endfor; cube loop for
		
		obj_cube=obj_cube[*,*,where(bads eq 0.)]

		klip_cube=obj_cube
		
		num=1 ;vestigal
		
		if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
			klip_cube_basis_arr = rdiklip(obj_cube, refcube_rdi, (size(refcube_rdi))[3], target=0,$
					posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,angsep=angsep,$
					anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,annmode_inout=annmode_inout,$
					basis=basis,ref_angles=ref_angles,filter=filter,bin=bin)
		endif
		
		if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
			restore,filename='~/Desktop/'+obj+'_filt_'+strcompress(filter,/r)+'_bin_'+strcompress(bin,/r)+'RDI_basis.sav'
			klip_cube_basis_arr = kl_basis_arr
		endif
		
		for ii=0, (size(obj_cube))[3]-1 do begin
			Print, '-------[ KLIPing image ', ii, ' out of ', (size(obj_cube))[3]-1, ' ]----------'
			print, 'Run:', runs
			
			if do_hyper ne 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
				trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
				pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr,$
				n_ang =n_ang, num=758) 
			
			if do_hyper eq 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
			   trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
			   pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang,$
			   num=758, /hyper) 
			
			if do_annmode eq 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask,$
			   distmask=distmask, trial=trial, posang=angles, wl=wl,diam=8.4, pixelscale=pxscale, angsep=angsep,$
			   anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758, annmode_inout=annmode_inout) 
			   
			if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
				klip_cube[*,*,ii] = rdiklip(obj_cube, refcube_rdi, k_klip,$
				target=ii,posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,$
				angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,$
				annmode_inout=annmode_inout, basis=0, kl_basis=klip_cube_basis_arr,$
				ref_angles=ref_angles,filter=filter,bin=bin)
			endif;annmode with rdi
			
			if do_annmode eq 1 and rdi eq 1 and basis eq 1 then begin
				klip_cube[*,*,ii] = rdiklip(obj_cube, refcube_rdi, k_klip,$
				target=ii,posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,$
				angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,$
				annmode_inout=annmode_inout, basis=1,$
				ref_angles=ref_angles,filter=filter,bin=bin)
			endif;annmode with rdi
			
		endfor; cube loop for

		for ii=0, (size(obj_cube))[3]-1 do begin
		   print, 'Rotating by ', angles[ii]
		   klip_cube[*,*,ii]=rot(klip_cube[*,*,ii],-angles[ii]-truenorth,/interp)
		   obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],-angles[ii]-truenorth,/interp)
		
		   framei=klip_cube[*,*,ii]
		   frameifull=obj_cube[*,*,ii]
		   
		   ;fill in the rest of the image
		   if fill eq 1 then framei[where(finite(framei) eq 0)]=frameifull[where(finite(framei) eq 0)]
		   klip_cube[*,*,ii]=framei
		endfor; cube loop for

	   ;writefits,output_folder+dither_folder+obj+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type +$
		;  '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+$
		 ; string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		  ;string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +$
		  ;'_cube_klip.fits', klip_cube
	   
		if combine_type eq 'median' then medarr, klip_cube, medframe
		if combine_type eq 'mean' then medframe = mean(klip_cube, dim=3, /double)
		if combine_type eq 'nwadi' then medframe = nw_ang_comb(klip_cube, angles)
	   
	   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip'+suffix+ '_bin_' +$
		  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
		  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe
   
	   print, 'PSF Width: ',width
	   PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	   PSFN = PSF/MAX(PSF)
	
	   medframe[where(finite(medframe) ne 1)]=0.
	   medframe_c = convolve(medframe, PSFN)
			 
	   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip_conv'+suffix+$
		  '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+$
		  '_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe_c
	
		klipframe=medframe
	   
	   if runs eq 1 then nods = klipframe else nods = [[[nods]], [[klipframe]]]
	   writefits, strcompress(cube_folder+ 'combined/' + obj + '_nods_klip' + suffix + '_bin_' +$
		  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
		  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits', /rem), nods
	endfor; Runs=1,4 for
	
	; Do this stuff at the end to combine
	left_klip = nods[*,*,0:1]
	right_klip = nods[*,*,2:3]
	
	; Get the factor we need to magnify the side with the greater pxscale by
	; Bigger pxscale means less pixels per arcsec so we need to scale it up to get
	; more pixels! Both are then at the smaller plate scale.
	
	if magnify eq 1 then begin
		pxscale = min_pxscale
		mag_factor = max([pxscale_sx, pxscale_dx]) / min([pxscale_sx, pxscale_dx])
		old_dim = (size(right_klip[*,*,0]))[1]
		new_dim = old_dim * mag_factor
		dim_diff = new_dim - old_dim
		start_i = dim_diff / 2
		end_i = (new_dim-1) - start_i
		
		if pxscale_sx gt pxscale_dx then begin
			new_left_klip = CONGRID(left_klip[*,*,0], new_dim, new_dim, /INTERP)
			
			new_left_klip = [[[new_left_klip]], [[CONGRID(left_klip[*,*,1], new_dim,$
				new_dim, /INTERP)]]]
			
			left_klip = new_left_klip[start_i:end_i, start_i:end_i, *]
		endif else begin
			new_right_klip = CONGRID(right_klip[*,*,0], new_dim, new_dim, /INTERP)
			
			new_right_klip = [[[new_right_klip]], [[CONGRID(right_klip[*,*,1], new_dim,$
				new_dim, /INTERP)]]]
			
			right_klip = new_right_klip[start_i:end_i, start_i:end_i, *]
		endelse
		
	endif; magnify if
	
	if combine_type eq 'median' then begin
		medarr, left_klip, left_klip_median
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_left_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), left_klip_median
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		left_klip_mean=mean(left_klip,dim=3,/double)
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_left_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), left_klip_mean
	endif; mean combine
	
	if combine_type eq 'median' then begin
		medarr, right_klip, right_klip_median
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), right_klip_median
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		right_klip_mean=mean(right_klip,dim=3,/double)
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), right_klip_mean
	endif; mean combine
	
	; write the combined left + right sides
	if combine_type eq 'median' then begin
		medarr, [[[left_klip_median]], [[right_klip_median]]], total_klip_median
		
		super_suffix = cube_folder + 'combined/' + obj + '_bin_' + string(sigfig(bin,1)) + '_type_' +$
		   bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+$
		   string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		   string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj)
		   
		if keyword_set(trial) then super_suffix += '_trial_' + string(sigfig(trial, 4))
		
		writefits, strcompress(super_suffix + '_total_klip.fits', /rem), total_klip_median
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		total_klip_mean=mean([[[left_klip_mean]], [[right_klip_mean]]],dim=3,/double)
		
		super_suffix = cube_folder + 'combined/' + obj + '_bin_' + string(sigfig(bin,1)) + '_type_' +$
		   bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+$
		   string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		   string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj)
		   
		if keyword_set(trial) then super_suffix += '_trial_' + string(sigfig(trial, 4))
		
		writefits, strcompress(super_suffix + '_total_klip.fits', /rem), total_klip_mean
	endif; mean combine
	
	writefits, strcompress(super_suffix + '_klip_nod1.fits', /rem), left_klip[*,*,0]
	writefits, strcompress(super_suffix + '_klip_nod2.fits', /rem), left_klip[*,*,1]
	writefits, strcompress(super_suffix + '_klip_nod3.fits', /rem), right_klip[*,*,0]
	writefits, strcompress(super_suffix + '_klip_nod4.fits', /rem), right_klip[*,*,1]
	;Where does this correction factor come from?
	; I'm having trouble with trying to *not* manually type in the folder here.
	if keyword_set(fs) then begin
		print, 'FS = ', fs
		if fs eq 1 then begin
			
			ref_file = output_folder+dither_folder+obj+string(ct) + '_pupil.fits'
			
			find_sources,strcompress(super_suffix + '_total_klip.fits',/rem),$
				reference=ref_file, platescale=pxscale,$
				correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
				fwhm=fwhm,ct=ct,filter=filter
				
		endif; find_sources if
	endif; keyword_set(fs) if
	
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(npixel=sz, FWHM=[width, width])
	PSFN = PSF / MAX(PSF)
	
	if combine_type eq 'median' then begin
		total_klip_convolve = convolve(total_klip_median, PSFN)
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		total_klip_convolve = convolve(total_klip_mean, PSFN)
	endif; mean combine
	
	writefits, strcompress(super_suffix +  '_total_klip_conv.fits', /rem), total_klip_convolve
endif; nod eq 'total' if

; NOT ALCOR
if nod eq 'sx_only' then begin
	for runs=1,2 do begin
		if rdi eq 1 then refcube_rdi = refcube
	   ; Do this for runs eq 1 and runs eq 3
	   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
	   ;Do this for runs eq 1 and runs eq 2
		  output_folder = cube_folder + 'processed_left/'
		  truenorth = truenorth_sx
		  annmode_inout = annmode_inout_sx
		  pxscale = pxscale_sx
	
		obj=strcompress(obj,/rem)
		
		if use_injection then obj_cube = readfits(output_folder + dither_folder + obj + string(ct) +$
		   '_cube_skysub_cen_clean_inj.fits') else obj_cube = readfits(output_folder + dither_folder +$
		   obj + string(ct) +  '_cube_skysub_cen_clean.fits')
		
		; destriping
		;if do_destripe eq 1 then begin
		;	print, 'destriping 90 degrees...'
		;	for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
		;	print, 'destriping 0 degrees...'
		;	for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],0.,clip_level=0.0,/nodisp)
		;endif
		
		; high-pass filter
		if filter gt 1 then for iii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,iii]=obj_cube[*,*,iii]-$
		   smooth(obj_cube[*,*,iii],filter)
		   
		restore,filename=output_folder + dither_folder + obj + string(ct) +  '_parang_clean.sav'
		
		; restrict (unbinned cube)
		if klip_fraction eq 1 then begin
		   obj_cube=obj_cube[*,*,start_frame:end_frame]
		   angles=angles[start_frame:end_frame]
		endif
		
		if bin gt 1 and bin_type eq 'median' then begin
		   st=1
		   for ii=0.,(size(obj_cube))[3]-1. do begin
			  print, fix(ii+1.) mod fix(bin)
			  if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=[ [[binned]],$
				 [[ obj_cube[*,*,ii] ]] ]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
			  
				 print,'Binning left frames...'
				 
				 if st eq 1 then begin 
				 
				 ;medarr,binned,binned
		
					binned=median(binned,dim=3,/even,/double)
				 binned_cube = binned 
				 st=0 
				 
				  endif else begin
				  ;medarr,binned,binned
		
					binned=median(binned,dim=3,/even,/double)
		
				  binned_cube=[[[binned_cube]],[[binned]]]
				  
				  endelse
				 print, size(binned_cube)
			  endif
		   endfor   
		   print, size(binned_cube)
		   st=1
		   for ii=0,(size(obj_cube))[3]-1. do begin
			  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+$
				 angles[ii]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
				 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
					[[[binned_angles]],[[binned_angle/bin]]]
			  endif
		   endfor   
		   
		   obj_cube=binned_cube
		   angles=binned_angles
		endif
		
		if bin gt 1 and bin_type eq 'mean' then begin
		   st=1
		   for ii=0.,(size(obj_cube))[3]-1. do begin
			  print, fix(ii+1.) mod fix(bin)
			  if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=binned + obj_cube[*,*,ii] 
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
			  
				 print,'Binning left frames...'
				 
				 if st eq 1 then begin 
				 
				 ;medarr,binned,binned
				 binned_cube = binned / float(bin)
				 st=0 
				 
				  endif else begin
				 ; medarr,binned,binned
		
				  binned_cube=[[[binned_cube]],[[binned/float(bin)]]]
				  
				  endelse
				 print, size(binned_cube)
			  endif
		   endfor   
		   print, size(binned_cube)
		   st=1
		   for ii=0,(size(obj_cube))[3]-1. do begin
			  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+$
				 angles[ii]
			  if fix(ii+1.) mod fix(bin) eq 0 then begin
				 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
					[[[binned_angles]],[[binned_angle/bin]]]
			  endif
		   endfor   
		   
		   obj_cube=binned_cube
		   angles=binned_angles
		endif
		
		lbc=0
		bads=fltarr( (size(obj_cube))[3] )
		bads[*]=0.
		
		obj_cube[where( finite( obj_cube ) eq 0 )]=0.
		
		for ii=0, (size(obj_cube))[3]-1 do begin
		   
		   if total(obj_cube[*,*,ii]) eq 0 then begin
			  print, ii
			  remove, ii-lbc, angles
			  bads[ii]=1.
			  lbc=lbc+1
		   endif
		   
		endfor
		obj_cube=obj_cube[*,*,where(bads eq 0.)]
		
		klip_cube=obj_cube
		
		num=1 ;vestigal
		
		if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
			klip_cube_basis_arr = rdiklip(obj_cube, refcube_rdi, (size(refcube_rdi))[3], target=0,$
					posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,angsep=angsep,$
					anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,annmode_inout=annmode_inout,$
					basis=basis,ref_angles=ref_angles,filter=filter,bin=bin)
		endif
		
		if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
			restore,filename='~/Desktop/'+obj+'_filt_'+strcompress(filter,/r)+'_bin_'+strcompress(bin,/r)+'RDI_basis.sav'
			klip_cube_basis_arr = kl_basis_arr
		endif
		
		for ii=0, (size(obj_cube))[3]-1 do begin
			Print, '-------[ KLIPing image ', ii, ' out of ', (size(obj_cube))[3]-1, ' ]----------'
			print, 'Run:', runs
		
			if do_hyper ne 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
				trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
				pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr,$
				n_ang =n_ang, num=758) 
		
			if do_hyper eq 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
				trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
				pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang,$
				num=758, /hyper) 
		
			if do_annmode eq 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask,$
				distmask=distmask, trial=trial, posang=angles, wl=wl,diam=8.4, pixelscale=pxscale, angsep=angsep,$
				anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758, annmode_inout=annmode_inout) 
				
			if do_annmode eq 1 and rdi eq 1 and basis eq 1 then begin
				klip_cube_basis_arr = rdiklip(obj_cube, refcube_rdi, k_klip, target=ii,$
						posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,angsep=angsep,$
						anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,annmode_inout=annmode_inout,$
						basis=1, ref_angles=ref_angles_rdi, filter=filter,bin=bin)
			endif
				
			if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
				klip_cube[*,*,ii] = rdiklip(obj_cube, refcube_rdi, k_klip,$
				target=ii,posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,$
				angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,$
				annmode_inout=annmode_inout, basis=0, kl_basis=klip_cube_basis_arr,$
				 ref_angles=ref_angles,filter=filter,bin=bin)
			endif;annmode with rdi
			
		endfor; adiklip loop over frames
		
		for ii=0, (size(obj_cube))[3]-1 do begin
		   print, 'Rotating by ', angles[ii]
		   klip_cube[*,*,ii]=rot(klip_cube[*,*,ii],-angles[ii]-truenorth,/interp)
		   obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],-angles[ii]-truenorth,/interp)
		
		   framei=klip_cube[*,*,ii]
		   frameifull=obj_cube[*,*,ii]
		   
		   ;fill in the rest of the image
		   if fill eq 1 then framei[where(finite(framei) eq 0)]=frameifull[where(finite(framei) eq 0)]
		   klip_cube[*,*,ii]=framei
		endfor; de-rotation to north-up
		
		 ;  writefits,output_folder+dither_folder+obj+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type +$
		;	  '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+$
		;	  string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		;	  string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +$
		;	  '_cube_klip.fits', klip_cube
		   
			if combine_type eq 'median' then medarr, klip_cube, medframe
			if combine_type eq 'mean' then medframe = mean(klip_cube, dim=3, /double)
			if combine_type eq 'nwadi' then medframe = nw_ang_comb(klip_cube, angles)
		   
		   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip'+suffix+ '_bin_' +$
			  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
			  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
			  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
			  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe
		
		   print, 'PSF Width: ',width
		   PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
		   PSFN = PSF/MAX(PSF)
		
		   medframe[where(finite(medframe) ne 1)]=0.
		   medframe_c = convolve(medframe, PSFN)
				 
		   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip_conv'+suffix+$
			  '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+$
			  '_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
			  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
			  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe_c
		
		klipframe=medframe
		   
		   if runs eq 1 then nods = klipframe else nods = [[[nods]], [[klipframe]]]
		   writefits, strcompress(cube_folder+ 'combined/' + obj + '_nods_klip' + suffix + '_bin_' +$
			  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
			  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
			  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
			  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits', /rem), nods
	endfor; Runs=1,2 for
	
	; Do this stuff at the end to combine
	left_klip = nods
	
	if combine_type eq 'median' then begin
		medarr, left_klip, left_klip_median
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_left_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), left_klip_median
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		left_klip_mean=mean(left_klip,dim=3,/double)
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_left_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), left_klip_mean
	endif; mean combine
	
	super_suffix = cube_folder + 'combined/' + obj + '_bin_' + string(sigfig(bin,1)) + '_type_' +$
	   bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,2))+'_angsep_'+$
	   string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	   string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj)
	   
	if keyword_set(trial) then super_suffix += '_trial_' + string(sigfig(trial, 4))
	
	writefits, strcompress(super_suffix + '_klip_nod1.fits', /rem), left_klip[*,*,0]
	writefits, strcompress(super_suffix + '_klip_nod2.fits', /rem), left_klip[*,*,1]
	;Where does this correction factor come from?
	; I'm having trouble with trying to *not* manually type in the folder here.
	if keyword_set(fs) then begin
		print, 'FS = ', fs
		if fs eq 1 then begin
			
			ref_file = output_folder+dither_folder+obj+string(ct) + '_pupil.fits'
			
			find_sources,strcompress(super_suffix + '_left_klip.fits',/rem),$
				reference=ref_file, platescale=pxscale,$
				correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
				fwhm=fwhm,ct=ct,filter=filter
				
		endif; find_sources if
	endif; keyword_set(fs) if
	
	;print, 'PSF Width: ',width
	PSF = psf_Gaussian(npixel=sz, FWHM=[width, width])
	PSFN = PSF / MAX(PSF)
	
	if combine_type eq 'median' then begin
		left_klip_convolve = convolve(left_klip_median, PSFN)
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		left_klip_convolve = convolve(left_klip_mean, PSFN)
	endif; mean combine
	
	writefits, strcompress(super_suffix +  '_left_klip_conv.fits', /rem), left_klip_convolve
endif; nod eq 'sx_only if'

if nod eq 'dx_only' then begin
	for runs=3,4 do begin
	   ; Do this for runs eq 1 and runs eq 3
	   if runs mod 2 then begin
	   		dither_folder = 'dith1/'
	  		not_dither_folder = 'dith2/'
	   endif else begin
	   		dither_folder = 'dith2/'
	   		not_dither_folder='dith1/'
	   	endelse
	
		  output_folder = cube_folder + 'processed_right/'
		  truenorth = truenorth_dx
		  annmode_inout = annmode_inout_dx
		  pxscale = pxscale_dx
	
		obj=strcompress(obj,/rem)
		
		if use_injection then obj_cube = readfits(output_folder + dither_folder + obj + string(ct) +$
		   '_cube_skysub_cen_clean_inj.fits') else obj_cube = readfits(output_folder + dither_folder +$
		   obj + string(ct) +  '_cube_skysub_cen_clean.fits')
		
		; destriping
	;	if do_destripe eq 1 then begin
	;		print, 'destriping 90 degrees...'
	;		for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
	;		print, 'destriping 0 degrees...'
	;		for ii=0, (size(obj_cube))[3]-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],0.,clip_level=0.0,/nodisp)
	;	endif
		
		; high-pass filter
		if filter gt 1 then for iii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,iii]=obj_cube[*,*,iii]-$
		   smooth(obj_cube[*,*,iii],filter)
		   
		restore,filename=output_folder + dither_folder + obj + string(ct) +  '_parang_clean.sav'
		print,newline, 'obj_cube restored angles size:', size(angles)
		
		; restrict (unbinned cube)
		;if (klip_fraction eq 1) then begin
		 ;  obj_cube=obj_cube[*,*,start_frame:end_frame]
		  ; angles=angles[start_frame:end_frame]
		;endif
		
		if bin gt 1 and bin_type eq 'median' then begin
		   st=1
		   for ii=0.,(size(obj_cube))[3]-1. do begin
			  print, fix(ii+1.) mod fix(bin)
			  
			  if fix(ii+1.) mod fix(bin) eq 1 then begin; start over on 1
			  		binned = obj_cube[*,*,ii]
			  		binned_angle = angles[ii]
			  endif else begin; concatenate, or add
			  		binned = [ [[binned]], [[ obj_cube[*,*,ii] ]] ]
			  		binned_angle = binned_angle + angles[ii]
			  
			  
			  if fix(ii+1.) mod fix(bin) eq 0 then begin; bin on 0
			  
				 print,'Binning Science frames...'
				 
				  if st eq 1 then begin ; first binning (set at top of loop)
				 	binned_cube = median(binned, dim=3, /even, /double)
				 	binned_angles = [binned_angle/bin]
				 	st=0 
				  endif else begin; st eq 1 if, else is for future binning
					new_binned = median(binned,dim=3,/even,/double)
					new_binned_angle = binned_angle / bin
				  	binned_cube=[[[binned_cube]],[[new_binned]]]; add to cube
				  	binned_angles=[[binned_angles],[new_binned_angle]]; add to angles
				  endelse; st neq 1
				  
			  endif; mod 0 if
			  
			  endelse; not starting over if
		   endfor   ; cube loop for
   
		   ;print, size(binned_cube)
		   ;print, size(binned_angles)
   
		   obj_cube = binned_cube
		   angles = binned_angles
		   ;print, newline, ref_angles_rdi, newline
		endif; bin gt 1
		
		lbc=0
		bads=fltarr( (size(obj_cube))[3] )
		bads[*]=0.

		obj_cube[where( finite( obj_cube ) eq 0 )]=0.

		for ii=0, (size(obj_cube))[3]-1 do begin
		   if total(obj_cube[*,*,ii]) eq 0 then begin
			  print, ii
			remove, ii-lbc, angles
			  bads[ii]=1.
			  lbc=lbc+1
		   endif
		endfor
		
		obj_cube=obj_cube[*,*,where(bads eq 0.)]
		
	;	if bin gt 1 and bin_type eq 'mean' then begin
	;	   st=1
	;	   for ii=0.,(size(obj_cube))[3]-1. do begin
	;		  print, fix(ii+1.) mod fix(bin)
	;		  if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=binned + obj_cube[*,*,ii] 
	;		  if fix(ii+1.) mod fix(bin) eq 0 then begin
	;		  
	;			 print,'Binning right frames...'
	;			 
	;			 if st eq 1 then begin 
	;			 	binned_cube = binned / float(bin)
	;			 	st=0
	;			endif else begin
	;				binned_cube=[[[binned_cube]],[[binned/float(bin)]]]
	;			endelse
	;			
	;			 print, size(binned_cube)
	;		  endif
	;	   endfor   
	;	   print, size(binned_cube)
	;	   st=1
	;	   for ii=0,(size(obj_cube))[3]-1. do begin
	;		  if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+$
	;			 angles[ii]
	;		  if fix(ii+1.) mod fix(bin) eq 0 then begin
	;			 if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=$
	;				[[[binned_angles]],[[binned_angle/bin]]]
	;		  endif
	;	   endfor   
	;	   
	;	   obj_cube=binned_cube
	;	   angles=binned_angles
	;	endif
	;	
	;;	
	;	lbc=0
	;	bads=fltarr( (size(obj_cube))[3] )
	;	bads[*]=0.
	;	
	;	obj_cube[where( finite( obj_cube ) eq 0 )]=0.
		
	;	for ii=0, (size(obj_cube))[3]-1 do begin
	;	   
	;	   if total(obj_cube[*,*,ii]) eq 0 then begin
	;		  print, ii
	;		  remove, ii-lbc, angles
	;		  bads[ii]=1.
	;		  lbc=lbc+1
	;	   endif
;
;		endfor
;		obj_cube=obj_cube[*,*,where(bads eq 0.)]
		
		klip_cube=obj_cube
	
		num=1 ;vestigal
		
		if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
			print, 'Restoring basis....'
			restore,filename='~/Desktop/'+obj+'_filt_'+strcompress(filter,/r)+'_bin_'+strcompress(bin,/r)+'RDI_basis.sav'
			klip_cube_basis_arr = kl_basis_arr
		endif
		
		for ii=0, (size(obj_cube))[3]-1 do begin
			Print, '-------[ KLIPing image ', ii, ' out of ', (size(obj_cube))[3]-1, ' ]----------'
			print, 'Run:', runs
		
			if do_hyper ne 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
				trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
				pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr,$
				n_ang =n_ang, num=758) 
		
			if do_hyper eq 1 and do_annmode ne 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii,$
				trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=wl,diam=8.4,$
				pixelscale=pxscale, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang,$
				num=758, /hyper) 
		
			if do_annmode eq 1 and rdi eq 0 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask,$
				distmask=distmask, trial=trial, posang=angles, wl=wl,diam=8.4, pixelscale=pxscale, angsep=angsep,$
				anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758, annmode_inout=annmode_inout)
								
			if do_annmode eq 1 and rdi eq 1 and basis eq 1 then begin
				klip_cube_basis_arr = rdiklip(obj_cube, refcube_rdi, k_klip, target=ii,$
						posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,angsep=angsep,$
						anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,n_ang=n_ang,annmode_inout=annmode_inout,$
						basis=1, ref_angles=ref_angles_rdi, filter=filter,bin=bin)
			endif
								
			if do_annmode eq 1 and rdi eq 1 and basis eq 0 then begin
				print, 'size(obj_cube):', size(obj_cube)
				klip_cube[*,*,ii] = rdiklip(obj_cube, refcube_rdi, k_klip,$
				target=ii,posang=angles,wl=wl,diam=8.4,pixelscale=pxscale,$
				angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,$
				n_ang=n_ang,annmode_inout=annmode_inout, basis=0,$
				kl_basis_arr=klip_cube_basis_arr, ref_angles=ref_angles_rdi,$
				filter=filter,bin=bin)
			endif;annmode with rdi
			 
		endfor; adiklip loop over frames
		
		for ii=0, (size(obj_cube))[3]-1 do begin
		   print, 'Rotating by ', angles[ii]
		   klip_cube[*,*,ii]=rot(klip_cube[*,*,ii],-angles[ii]-truenorth,/interp)
		endfor; de-rotation to north-up
		
	 ;  writefits,output_folder+dither_folder+obj+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type +$
	;	  '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+$
	;	  string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	;	  string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +$
	;	  '_cube_klip.fits', klip_cube
	  
		if combine_type eq 'median' then medarr, klip_cube, medframe
		if combine_type eq 'mean' then medframe = mean(klip_cube, dim=3, /double)
		if combine_type eq 'nwadi' then medframe = nw_ang_comb(klip_cube, angles)
	   
	   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip'+suffix+ '_bin_' +$
		  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
		  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe
	
	   print, 'PSF Width: ',width
	   PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	   PSFN = PSF/MAX(PSF)
	
	   medframe[where(finite(medframe) ne 1)]=0.
	   medframe_c = convolve(medframe, PSFN)
			 
	   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip_conv'+suffix+$
		  '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+$
		  '_k_klip_'+string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe_c
		
		klipframe=medframe
		   
	   ; first vs 2nd nod
	   if runs eq 3 then nods = klipframe else nods = [[[nods]], [[klipframe]]]
	   writefits, strcompress(cube_folder+ 'combined/' + obj + '_nods_klip' + suffix + '_bin_' +$
		  string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
		  string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
		  string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
		  string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits', /rem), nods
	endfor; Runs=3,4 for
	
	; Do this stuff at the end to combine
	right_klip = nods
	
	if combine_type eq 'median' then begin
		medarr, right_klip, right_klip_median
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), right_klip_median
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		right_klip_mean=mean(right_klip,dim=3,/double)
		
		writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), right_klip_mean
	endif; mean combine
	
	super_suffix = cube_folder + 'combined/' + obj + '_bin_' + string(sigfig(bin,1)) + '_type_' +$
	   bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,2))+'_angsep_'+$
	   string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
	   string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj)
	   
	if keyword_set(trial) then super_suffix += '_trial_' + string(sigfig(trial, 4))
	
	writefits, strcompress(super_suffix + '_klip_nod3.fits', /rem), right_klip[*,*,0]
	writefits, strcompress(super_suffix + '_klip_nod4.fits', /rem), right_klip[*,*,1]
	;Where does this correction factor come from?
	; I'm having trouble with trying to *not* manually type in the folder here.
	if keyword_set(fs) then begin
		print, 'FS = ', fs
		if fs eq 1 then begin
			
			ref_file = output_folder+dither_folder+obj+string(ct) + '_pupil.fits'
			
			find_sources,strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem),$
				reference=ref_file, platescale=pxscale,$
				correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
				fwhm=fwhm,ct=ct,filter=filter
				
		endif; find_sources if
	endif; keyword_set(fs) if
	
	PSF = psf_Gaussian(npixel=sz, FWHM=[width, width])
	PSFN = PSF / MAX(PSF)
	
	if combine_type eq 'median' then begin
		right_klip_convolve = convolve(right_klip_median, PSFN)
	endif; median combine
	if (combine_type eq 'mean') or (combine_type eq 'nwadi') then begin
		right_klip_convolve = convolve(right_klip_mean, PSFN)
	endif; mean combine
	
	writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip_conv' + suffix + '_bin_' +$
	   string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
	   string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
	   string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
	   string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), right_klip_convolve
	   
endif; nod eq 'dx_only if'
   
end; end procedure
