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

pro stdev_photometry, coadd=coadd, type=type
   ; Type is 'ADI' or 'KLIP'

   COMPILE_OPT IDL2
   newline = string(10B)
   start_time = systime(/JULIAN); Get the current time in a Julian date format

   ;------------------------------[ Start User Input ]---------------------------------

   ; Where to find our files and put the results
   output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_'+strcompress(coadd,/r)+'/'
   obj = 'HII1348'

   ; Parameters needed to read in our total_klip or total_adi file
   bin = 3 & bin_type = 'mean' & combine_type = 'nwadi'
   k_klip = 7 & angsep= 1. & anglemax = 360. & nrings = 4.
   n_ang = 2 & do_cen_filter = 1 & filter = 17. & ct = 0.994

   ; Inject_planets parameters
   use_gauss = 0; Fit a gaussian to the pupil median PSF to try and get rid of the 'lobes'
   n_planets = 48; Includes HII 1348 b!
   pxscale = 0.0107 ; arcsec/pixel
   contrast = 0.008914755; Average of last two runs with photometry.pro and made positive
   fwhm = 8.72059 ; px ``width'' in reduce_lbti_HII1348.pro
   real_theta = 1.312585; rad, from photometry.pro results
   planet_r = 1.14512; arcsec, from photometry.pro results

   ;------------------------------[ End User Input ]---------------------------------

   ; Initialize arrays for our results
   rec_xxs=[] & rec_yys=[] & rec_rhos=[] & rec_thetas=[] & rec_contrasts=[]
   
   dither_folder = 'processed_left/dith1/'; Just choose dith 1 for simplicity

   obj_cube = readfits(output_path + dither_folder + obj + string(ct) +  '_cube_skysub_cen_clean.fits')
   ;unsaturated data can just use the pupil image
   ref = readfits(output_path + dither_folder + obj + string(ct) +  '_pupil.fits')
   ;if use_gauss then ref = gauss2dfit(ref, /tilt); I don't think that I really need this here.

   big_ref=obj_cube[*,*,0]
   big_ref[*]=0

   ;build bigger reference image
   for xx=0, (size(ref))[1] - 1 do begin
      for yy=0, (size(ref))[1] - 1 do begin

         big_ref[250.-(size(ref))[1]/2.+xx,250.-(size(ref))[1]/2.+yy]=ref[xx,yy]

      endfor
   endfor
   
   writefits, strcompress(output_path+'stdev_photometry/ref/'+obj+'_reference.fits', /rem), big_ref
   ; flux of what we recover dividied by the flux of this reference gives us a recovered contrast
   
   aper,big_ref,0,0,ref_flux,ref_fluxerr,sky,skyerr,1.75,aper_rad,sky_rad,[-99E99,999E99],/flux,SETSKYVAL=0,/exact
   ; ref_flux is now saved

   ; Default to ADI for now
   if not keyword_set(type) then type = 'ADI'

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
      og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+'filt_'+$
         string(filter)+'_neg_inj_'+string(0)+'_total_adi.fits',/rem))
   endif
   print, 'Original image read, finding companion centroid...', newline

   print, 'Starting loop over thetas'
   ; Create an array of n_planets which start at the position of HII 1348 b
   thetas = (2*!DPI * findgen(n_planets) / n_planets) + real_theta
   ; remove the three closest to HII 1348 b
   remove, [0, 1, (size(thetas))[1]-1], thetas
   thetas[where(thetas ge 2*!DPI)] -= 2*!DPI ; Keep angles on [0, 2pi) rad
   
   ; Calculate injected x and y values
   xxs = planet_r * (1./pxscale) * Cos(thetas) + 250.
   yys = planet_r * (1./pxscale) * Sin(thetas) + 250.
   
   trial = 0
   foreach theta, thetas do begin
      
      ; Do our injections
      hii1348_pipeline, rho=planet_r, theta=theta, contrast=contrast, pre_inj_stuff=0,$
         neg_inj=0, uncert=1, trial=trial, coadd=coadd, use_gauss=use_gauss; Inject and run ADI

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
      uncert_image=readfits(strcompress(output_path+'combined/'+obj+'ct_'+string(ct)+'filt_'+ $
         string(filter)+'_neg_inj_'+string(0)+'_uncert_'+string(1)+'_total_adi.fits',/rem))
      ;           endif
      print, 'Read in image after negative injection'
      print, 'Finding injected planets...', newline
      ; Calculate the centroid of the companion
      ; CNTRD was giving me issues, so I'm just going to use GCNTRD here.
      gcntrd, uncert_image, xxs[trial], yys[trial], XCEN, YCEN, fwhm
      print, 'GCNTRD xcen, ycen:', XCEN, YCEN
      print, 'Injected xx, yy:', xxs[trial], yys[trial]
      cen_x = XCEN & cen_y = YCEN
      
      ; Calculate rho, theta from cen_x, cen_y
      cen_rho = pxscale * SQRT(((cen_x-250.)^2.)+((cen_y-250.)^2.))
      ; radians (ccw from top at 0 rad)
      cen_theta = ATAN((cen_y-250.)/(cen_x-250.))
      if (cen_x-250.) lt 0 then cen_theta += !DPI; Accounting for limited range of arctan
      if cen_theta lt 0 then cen_theta += 2*!DPI; only use positive angles for consistency
      print, 'cen_rho, cen_theta:', cen_rho, cen_theta
      print, 'Injected rho, theta:', planet_r, theta, newline
      
      
      ; Save our results  to arrays
      rec_xxs=[rec_xxs,cen_x] & rec_yys=[rec_yys,cen_y]
      rec_rhos=[rec_rhos,cen_rho] & rec_thetas=[rec_thetas,cen_theta]
      
      print, 'Done.'+newline+'Writing FITS...'
      writefits, strcompress(output_path+'stdev_photometry/'+obj+'_trial_'+string(sigfig(trial,4))+'.fits', /rem), uncert_image
      print, 'Done.'+newline+'Incrementing trial...'
      trial += 1
   endforeach

   ; Save our results
   print, 'Saving...'
   save,filename=output_path+'stdev_photometry/'+obj+'_negative_inj_data.sav',xxs,yys,thetas,rec_xxs,rec_yys,rec_rhos,rec_thetas
   print, 'Done.'
   
   ; Combine all of the trials into a cube and write it to the same folder
   print, newline, 'Saving the trials into one FITS cube'
   folder_cube, output_path+'stdev_photometry/'
   print, 'FITS cube created! Starting analysis', newline
   rhos = fltarr(n_planets-1) + planet_r; Create an array of HII 1348 b radii
   ; Get vectors of differences between injected and recovered values
   x_diff=xxs-rec_xxs & y_diff=yys-rec_yys & rho_diff=rhos-rec_rhos & theta_diff=thetas-rec_thetas
   x_uncert=STDDEV(x_diff) & y_uncert=STDDEV(y_diff)
   
   x = 277.328 & y = 353.4725
   rho_uncert = 0.5*sqrt( (x_uncert^2. / (2.*(x-250.))) + (y_uncert^2. / (2.*(y-250.))) )
   theta_uncert = sqrt( ((y-250.)^2. * x_uncert^2.) + ((x-250.)^2. * y_uncert^2.) ) / (y^2. - 500.*y + x^2. + 125000. - 500.*x)
   
   
   avg_xdiff = mean(x_diff) & avg_ydiff = mean(y_diff) & avg_rho_diff = mean(rho_diff) & avg_theta_diff = mean(theta_diff)
   
   print, 'x (RA) uncertainty:', x_uncert, ' (px) = ', x_uncert * pxscale, ' (arcsec)'
   print, 'mean x (RA) difference:', avg_xdiff, ' (px) = ', avg_xdiff * pxscale, '(arcsec)', newline
   
   print, 'y uncertainty:', y_uncert, ' (px) = ', y_uncert * pxscale, ' (arcsec)'
   print, 'mean y (DEC) difference:', avg_ydiff, ' (px) = ', avg_ydiff * pxscale, '(arcsec)', newline
   
   ; I need to do actual error propagation for these, I think?!?!?!?!!?
   print, 'rho uncertainty:', rho_uncert, ' (arcsec)'
   print, 'mean rho (SEP) difference:', avg_rho_diff, ' (arcsec)', newline
    
   ; Converted from (rad) to (arcsec)
   print, 'theta uncertainty:', !RADEG * 3600. * theta_uncert, ' (arcsec)'
   print, 'mean theta (PA) difference:', !RADEG * 3600. * avg_theta_diff, ' (arcsec)', newline
   
   print, 'Completed photometry uncertainties in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!