; aper does aperture photometry, cntrd gets the "center of light" for the sources in the apertures

; I probably need to go full grid-search and loop over negative contrasts at all positions as well.
; I think it makes sense to pass the xx and yy values straight into the planet injection as well,
; as something could go slightly wrong with either conversion to and from radii and angles, but
; I've been having problems with this so for right now the injection is passed through a radius and
; an angle.
; The centroiding should get me within a pixel of the true best negative injection center

; stdev of differences in astrometry (picked up & injected) are uncertainties in astrometry and photometry
; make sure that there isn't a big difference in one axis and not the other, but otherwise just worrying
; about the absolute differences in position. Final results in PA and sep. More uncertainty in PA than in
; sep because of the orientation of LBTI

; another way to do it is to just take stdev of recovered fluxes and their mean being
; as close to the actual companion to get the photometric uncertainty (divide by the mean to get a
; percent error)

pro uncert, coadd=coadd, type=type, planet_spots=planet_spots, use_gauss=use_gauss,$
	write=write, nod=nod
   ; Type is 'ADI' or 'KLIP'
   ; Nod is 'total', '1', '2', '3', or '4'
   ; planet_spots is meant to be a multiple of 16 (the number that fit naturally)

   COMPILE_OPT IDL2
   newline = string(10B)
   start_time = systime(/JULIAN); Get the current time in a Julian date format

   ;------------------------------[ Start User Input ]---------------------------------

   ; Where to find our files and put the results
   
	output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_'+$
		strcompress(coadd,/r)+'_'+strcompress(nod,/r)+'/'
   	
   obj = 'HII1348'

   ; Parameters needed to read in our total_klip or total_adi file
   bin = 3 & bin_type = 'mean' & combine_type = 'nwadi'
   k_klip = 7 & angsep= 1. & anglemax = 360. & nrings = 4.
   n_ang = 2 & do_cen_filter = 1 & filter = 17. & ct = 0.994

   ; Inject_planets parameters
   if not keyword_set(planet_spots) then planet_spots = 16; Includes HII 1348 b!
   if not keyword_set(use_gauss) then use_gauss = 0
	pxscale_sx = 0.010648 ; +0.000039 or -0.000050 arcsec/pixel (from Steve and Jared)
	pxscale_dx = 0.010700 ; +0.000042 or -0.000051 arcsec/pixel (from Steve and Jared)

	min_pxscale = min([pxscale_sx, pxscale_dx]); Both are matched to this in KLIP/ADI
	
   contrast_total = 0.00940411; from photometry.pro and made positive
   fwhm = 8.72059 ; px ``width'' in reduce_lbti_HII1348.pro
   real_theta_total = 1.31301; rad, from photometry.pro results
   planet_r_total = 1.14026; arcsec, from photometry.pro results
   x_total = 277.300 & y_total = 353.549
   
   ; four separate nod contrasts
   contrast_1 = 0.00958910
   contrast_2 = 0.00914963
   contrast_3 = 0.00969777
   contrast_4 = 0.00924111
   
   ; four separate nod radii in arcsec
   planet_r_1 = 1.14132
   planet_r_2 = 1.13583
   planet_r_3 = 1.14147
   planet_r_4 = 1.14188
   
   ; four separate nod x-positions
   x_1 = 277.133
   x_2 = 277.272
   x_3 = 277.419
   x_4 = 277.377
   
   ; four separate nod y-positions
   y_1 = 353.696
   y_2 = 353.125
   y_3 = 353.635
   y_4 = 353.685
   
   ; four separate nod thetas in radians
   real_theta_1 = 1.31488
   real_theta_2 = 1.31226
   real_theta_3 = 1.31215
   real_theta_4 = 1.31265
   
   ; aperture photometry parameters
   aper_rad = fwhm/2.
   sky_rad=[5,7] ;does not do anything, just needs to be set

   ;------------------------------[ End User Input ]----------------------------
   ; Make a folder to put our stuff
   file_mkdir, strcompress(output_path + 'photometry/uncert_' + nod, /r)
   
   ; Set values for our specific nod, if necessary
	if nod eq 'total' then begin
		contrast = contrast_total
		real_theta = real_theta_total
		planet_r = planet_r_total
		x = x_total
		y = y_total
	endif
	if nod eq '1' then begin
		contrast = contrast_1
		real_theta = real_theta_1
		planet_r = planet_r_1
		x = x_1
		y = y_1
	endif
	if nod eq '2' then begin
		contrast = contrast_2
		real_theta = real_theta_2
		planet_r = planet_r_2
		x = x_2
		y = y_2
	endif
	if nod eq '3' then begin
		contrast = contrast_3
		real_theta = real_theta_3
		planet_r = planet_r_3
		x = x_3
		y = y_3
	endif
	if nod eq '4' then begin
		contrast = contrast_4
		real_theta = real_theta_4
		planet_r = planet_r_4
		x = x_4
		y = y_4
	endif

   ; Initialize arrays for our results
   rec_xxs=[] & rec_yys=[] & rec_rhos=[] & rec_thetas=[] & rec_fluxes=[]

   ; Default to ADI for now
   if not keyword_set(type) then type = 'ADI'
   
   if not keyword_set(write) then write = 0
   
   if planet_spots mod 16 then begin
      print, "Error enter an integer multiple of 16 for planet_spots"
      stop
   endif

   print, 'Reading in total ' + type +  ' image...'
   ;if type eq 'KLIP' then begin
   ;
   ;   og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
   ;      +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
   ;      string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
   ;      string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
   ;      '_total_klip.fits', /rem))
   ;
   ;endif
   
	if nod eq 'total' then begin

		print, 'Reading in total ' + type +  ' image...'
		;if type eq 'KLIP' then begin
		;  
		;   og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
		;      +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
		;      string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		;      string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
		;      '_total_klip.fits', /rem))
		;    
		;endif

		if type eq 'ADI' then begin
			og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+$
			'filt_'+string(filter)+'_neg_inj_'+string(0)+'_uncert_0_total_adi.fits',/rem))
		endif
	
	endif else begin; nod eq 'total' if
	
		print, 'Reading in nod ' + type +  nod + ' image...'
		;if type eq 'KLIP' then begin
		;  
		;   og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
		;      +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
		;      string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		;      string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
		;      '_total_klip.fits', /rem))
		;    
		;endif

		if type eq 'ADI' then begin
			og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+$
			'filt_'+string(filter)+'_neg_inj_'+string(0)+'_uncert_0_adi_nod' + nod +$
			'.fits',/rem))
		endif

	endelse; nod ne 'total' else
	
   x_size = (size(og_image))[1] & y_size = (size(og_image))[2]
	;half-sizes
	x_hsize = x_size / 2.
	y_hsize = y_size / 2.
	print, x_hsize
	print, y_hsize

   print, 'Original image read', newline
   print, 'Starting loop over thetas'
   
   ; Create an array of planet_spots which start at the position of HII 1348 b
   thetas = (2*!DPI * findgen(planet_spots) / planet_spots) + real_theta
   
   ; remove those on and too close to HII 1348 b
   n_to_remove = planet_spots / 16
   remove, 0, thetas; Get rid of the one on top of HII 1348 b
   ; Remove those that are too close (always the first or the last element)
   for j = 0,((n_to_remove-1)/2) do begin
      remove, 0, thetas
      remove, (size(thetas))[1]-1, thetas
   endfor; j removal for
   
   thetas[where(thetas gt 2*!DPI)] -= 2*!DPI ; Keep angles on [0, 2pi) rad
   
   ; Calculate injected x and y values
   xxs = (planet_r * (1./min_pxscale) * Cos(thetas)) + x_hsize
   yys = (planet_r * (1./min_pxscale) * Sin(thetas)) + y_hsize
   
   trial = 0
   foreach theta, thetas do begin
      ; Do our injections
      hii1348_pipeline, rho=planet_r, theta=theta, contrast=contrast, pre_inj=0,$
         neg_inj=0, uncert=1, trial=trial, coadd=coadd, use_gauss=use_gauss, klip=0,$
         fs=0, extra=nod, nod=nod; Inject and run ADI

      ; Read in the total KLIP or ADI file after the negative injection
      print, 'Reading in neg-injected file'
      
      ;           if type eq 'KLIP' then begin
      ;
      ;              image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+$
      ;                 '_type_'+bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+$
      ;                 '_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'$
      ;                 +string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'_neg_inj_'+string(1)+$
      ;                 '_total_klip.fits', /rem))
      ;
      ;           endif
      ;           if type eq 'ADI' then begin (only real option for now)
      ;           endif
      
      if nod eq 'total' then begin
				
			uncert_image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+$
				string(ct)+'filt_'+string(filter)+'_neg_inj_0'+$
				'_uncert_1_total_adi.fits',/rem))
						
		endif else begin; nod eq 'total' if
				
			uncert_image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+$
				string(ct)+'filt_'+string(filter)+'_neg_inj_0'+$
				'_uncert_1_adi_nod'+nod+'.fits',/rem))
						
		endelse; nod neq 'total' else
         
      print, 'Read in image after negative injection'
      print, 'Finding injected planets...', newline
      ; Calculate the centroid of the companion
      ; CNTRD was giving me issues, so I'm just going to use GCNTRD here.
      gcntrd, uncert_image, xxs[trial], yys[trial], XCEN, YCEN, fwhm
      print, 'GCNTRD xcen, ycen:', XCEN, YCEN
      print, 'Injected xx, yy:', xxs[trial], yys[trial]
      cen_x = XCEN & cen_y = YCEN
      
      ; Do some photometry on our recovered source
      aper,uncert_image,cen_x,cen_y,rec_flux,rec_fluxerr,sky,skyerr,1.75,aper_rad,$
      	sky_rad,[-99E99,99E99],/flux,SETSKYVAL=0,/exact,/nan
      
      ; Calculate rho, theta from cen_x, cen_y
      cen_rho = min_pxscale * SQRT( ((cen_x-x_hsize)^2.)  +((cen_y-y_hsize)^2.) )
      ; radians (ccw from top at 0 rad)
      cen_theta = ATAN( (cen_y-y_hsize) / (cen_x-x_hsize) )
      if (cen_x - x_hsize) lt 0 then cen_theta += !DPI; Accounting for limited range of arctan
      if cen_theta lt 0 then cen_theta += 2*!DPI; only use positive angles for consistency
      print, 'cen_rho, cen_theta:', cen_rho, cen_theta
      print, 'Injected rho, theta:', planet_r, theta, newline
      
      ; Save our results  to arrays
      
      rec_xxs=[rec_xxs,cen_x] & rec_yys=[rec_yys,cen_y] & rec_fluxes=[rec_fluxes,$
      	rec_flux]
      	
      rec_rhos=[rec_rhos,cen_rho] & rec_thetas=[rec_thetas,cen_theta]
      
      if write eq 1 then begin
      	print, 'Done.'+newline+'Writing FITS...'
      	writefits, strcompress(output_path+'photometry/uncert_'+nod+'/'+obj+'_trial_'+$
      		string(sigfig(trial,4))+'.fits', /rem), uncert_image
      endif
      
      print, 'Done.'+newline+'Incrementing trial...'
      trial += 1
   endforeach; thetas foreach

   ; Save our results
   print, 'Saving...'
   
   save,filename=output_path+'photometry/uncert_'+nod+'/'+obj+'_inj_data.sav',xxs,yys,$
   	thetas,rec_xxs,rec_yys,rec_rhos,rec_thetas,rec_fluxes
   	
   print, 'Done.'
   
   ; Combine all of the trials into a cube and write it to the same folder
   
   if write eq 1 then begin
   	print, newline, 'Saving the trials into one FITS cube'
   	
   	folder_cube, output_path+'photometry/uncert_'+nod+'/',$
   		output_path+'photometry/uncert_'+nod+'/', type='list'
   		
   endif
   
   print, 'FITS cube created! Starting analysis', newline
   
   ; Create an array of HII 1348 b radii
   rhos = fltarr((size(thetas))[1]) + planet_r
   
   ; Get vectors of differences between injected and recovered values
   x_diff=xxs-rec_xxs & y_diff=yys-rec_yys & rho_diff=rhos-rec_rhos & theta_diff=thetas-rec_thetas
   x_uncert=STDDEV(x_diff) & y_uncert=STDDEV(y_diff)
   
   ; Get SEP and PA uncertainties with error propagation
   rho_uncert = 0.5*sqrt( (x_uncert^2. / (2.*(x-x_hsize))) + (y_uncert^2. /$
   	(2.*(y-y_hsize))) )
   	
   theta_uncert = sqrt( ((y-y_hsize)^2. * x_uncert^2.) + ((x-x_hsize)^2. *$
   	y_uncert^2.) ) / (y^2. - 2.*y_hsize*y + x^2. + 2.*y_hsize*x_hsize - 2.*x_hsize*x)
   
   ; Some average differences to make sure that their isn't an obvious systematic error
   avg_xdiff = mean(x_diff) & avg_ydiff = mean(y_diff)
   avg_rho_diff = mean(rho_diff) & avg_theta_diff = mean(theta_diff)
   
   ; Get contrast uncertainty
   flux_err = stddev(rec_fluxes) & rel_flux_err = flux_err / mean(rec_fluxes)
   print, 'Injected and recovered: ', (size(thetas))[1], ' artificial planets'
   print, 'Flux error: ', flux_err
   print, 'Relative flux error: ', rel_flux_err * 100, ' %', newline
   
   ; Print our results to the terminal:
   print, 'x (RA) uncertainty:', x_uncert, ' (px) = ', x_uncert * min_pxscale, ' (arcsec)'
   print, 'mean x (RA) difference:', avg_xdiff, ' (px) = ', avg_xdiff * min_pxscale, '(arcsec)', newline
   
   print, 'y uncertainty:', y_uncert, ' (px) = ', y_uncert * min_pxscale, ' (arcsec)'
   print, 'mean y (DEC) difference:', avg_ydiff, ' (px) = ', avg_ydiff * min_pxscale, '(arcsec)', newline
   
   print, 'rho uncertainty:', rho_uncert, ' (arcsec)'
   print, 'mean rho (SEP) difference:', avg_rho_diff, ' (arcsec)', newline
    
   ; Converted from (rad) to (deg)
   print, 'theta uncertainty:', !RADEG * theta_uncert, ' (deg)'
   print, 'mean theta (PA) difference:', !RADEG * avg_theta_diff, ' (deg)', newline
   
   save, filename=output_path+'photometry/uncert_'+nod+'/'+obj+'_inj_calc_results.sav',$
   	flux_err, rel_flux_err, x_uncert, avg_xdiff, y_uncert, avg_ydiff, rho_uncert,$
   	avg_rho_diff, theta_uncert, avg_theta_diff
   
   print, 'Completed photometric and astrometric uncertainties in ',$
		(systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!