pro contrastcurve_BD45,display=display

;USER INPUTS
pas=( findgen(72)*5.0 ) ;+150.

;filter='F410M' 

d=156 ;pc
obj_id='BD45'
star_name='BD45'

filter='Lp'

threshold=5.		;SNR recovery threshold
precision=0.1 		;how close to the threshold should each retrieved value be?
					;NOTE: smaller values require much longer run time
pxscale=0.0107

contrast_guess=5e-3
limit_to_guess=1	;will limit all trials to ≤ the above value

style_1='klip' 		;reduction style for inner region
style_2='adi'		;reduction style for outer region
style_thresh=10. ;arcsec cutoff to switch to reduction style 2


;include blackbody planet models with some extra heat
internal_heating_1=0.5
internal_heating_2=0.1
albedo=0.3 ;(between 0 and 1; E,N,U,S are about 0.3, jupiter is about 0.5)
Rpl=[0.,0.09,1.7*0.09,0.4,1.,2.0]	;Jupiter radii (Saturn is 0.9 RJ, Neptune is 0.3 RJ, 
									;Earth is 0.09 RJ) 

;synthetic box filter for planet models
startwave=9.8
endwave=12.4
wavestep=0.1 

 y_range=[0.5e-7,5e-4]	;only affects plotting




	separations=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,3.0]

root_dir='/Volumes/RAID36TB/LMIRCam/BD45598/'
reduced_dir='reduced/'


if pxscale gt 0.05 then x_range=[0.5,1.1*max(separations)] else  x_range=[0.25,1.1*max(separations)]



;-------------------- END USER INPUT -------------------------------


cgLoadCT, 1

;save start time
contstarttime=systime(/JULIAN)




;constants to be used later (for planet models)
au2sun=215.032 ;au in units of Solar radii
RJ2sun=0.10049 ;RJ in units of Solar radii
mic2ang=1E4 ;1 micron in units of angstroms

;set up arrays
results=separations & results[*]=0.
contrasts=separations & contrasts[*]=contrast_guess
curves=[]
testn=0

;start a loop over position angle
for m=0,n_elements(pas)-1 do begin

PA=separations & PA[*]=pas[m]*!DTOR ;set all same PA for now

if m eq 1 then contrasts=curves > 0
if m gt 1 then contrasts=mean(curves,dim=2,/nan) > 0
;start a loop over separation
for n=0,n_elements(separations)-1 do begin
	testn=testn+1
	plot_ind=1
	;uncomment this line to end up with a plot with a single curve at the end
	;if m eq n_elements(pas)-1 and n eq n_elements(separations)-1 then plot_ind=0 else plot_ind=1
	print, 'Working on test separation',n,'/',n_elements(separations)-1
	result=99E99 ;start significance off at infinity
	runs=0
	while result gt 1.1*threshold or result lt 0.9*threshold do begin
		
		if limit_to_guess then contrasts=contrasts < contrast_guess $
			else contrasts=contrasts < 1.0 ; just so that everything makes sense
		contrasts=contrasts>1e-8

		rho=separations[n] & theta=PA[n] & contrast=contrasts[n]
		if n ge 1 and finite(contrast) eq 0 then contrast=contrasts[n-1]
		if n le 1 and finite(contrast) eq 0 then contrast=contrasts_guess
		print, 'Testing rho, theta, contrast = ',rho,theta/!DTOR,contrast
		runs=runs+1
		;set boundaries to remove the line (NEAR SPECIFIC)
		
		;aa=0. & bb=0.  & ab=0. & ba=0. ;fix to zero
		;perform reduction
		
		transmission=1.0 ;no coronagraph, yet
		contrasti=contrast
		contrast=contrast*transmission ;correct for throughput (transmission)
		
		;limit to a certain range
		
		style=style_1
		if rho gt style_thresh then style=style_2
		print, 'Theta = ',theta
		if style eq 'klip' then reduce_lmircam_BD45,rho=rho,theta=theta,contrast=contrast,/klip;,/debug

		if style eq 'adi' then reduce_lmircam_BD45,rho=rho,theta=theta,contrast=contrast,/adi


		contrast=contrasti

			if style eq 'klip' then image=readfits(root_dir+reduced_dir+'BD45_left_klip.fits',hd)
			if style eq 'adi' then image=readfits(root_dir+reduced_dir+'BD45_left_adi.fits',hd)
		
			xs=(size(image))(1) & ys=(size(image))(2)
			xhs=float(xs)/2. & yhs=float(ys)/2. ;half-sizes
			xs=(size(image))(1) & ys=(size(image))(2)
			xhs=float(xs)/2. & yhs=float(ys)/2. ;half-sizes
		
			for xx=0,xhs*2.0-1 do for yy=0,yhs*2.0-1 do if sqrt((xx-xhs)^2. + (yy-yhs)^2.) lt rho/pxscale - 12 then image[xx,yy]=!values.f_nan
			for xx=0,xhs*2.0-1 do for yy=0,yhs*2.0-1 do if sqrt((xx-xhs)^2. + (yy-yhs)^2.) gt rho/pxscale + 13 then image[xx,yy]=!values.f_nan
			;truncate wild values
			;zmi=where(finite(image) eq 0)
			;image=image>(mean(image,/nan)-3.0*stddev(image,/nan) ) ;exclude pixels that are lt 3x the standard deviation of other pixels
			;image[zmi]=!values.f_nan ;mid step to replace all the close to zero values
			;image=image<(mean(image,/nan)+5.0*stddev(image,/nan) )
			;image[zmi]=!values.f_nan
		
	
			;old combination methods
			;image=mean([[[image1*sqrt(tt1/tt)]],[[image2*sqrt(tt2/tt)]],[[image3*sqrt(tt3/tt)]],[[image4*sqrt(tt4/tt)]]],dim=3,/nan) 
			;image=mean([[[image1*(tt1/tt)]],[[image2*(tt2/tt)]],[[image3*(tt3/tt)]],[[image4*(tt4/tt)]]],dim=3,/nan) 
			;image=image1+image2+image3+image4

			;weighted mean image combination for dealing with NaNs correctly
			;read in an image just to set the size
			;image=readfits(root_dir+reduced_dir+filter+'_'+style+'_inj.fits',hd)
			
			
			;imagecube=[[[image1]],[[image2]],[[image3]],[[image4]]]
			
			;store things in the first folder	

		xs=(size(image))(1) & ys=(size(image))(2)
		xhs=float(xs)/2. & yhs=float(ys)/2. ;half-sizes
		
		xp=xhs + ( (rho/pxscale) * Cos(theta)  ) & yp=yhs + ( (rho/pxscale) * Sin(theta)  )
		
		;re-position on the peak
		boxehs=5.
		;maxxpyp=max(image[xp-boxehs:xp+boxehs-1,yp-boxehs:yp+boxehs-1],xpyp )
		;xp=xpyp[0] & yp=xpyp[1]
		;xp = (xpyp MOD (2.*boxehs)) + xp - boxehs
		;yp = ((xpyp)/(2.*boxehs)) + yp - boxehs
		
		;cimage=image[xp-boxehs:xp+boxehs-1,yp-boxehs:yp+boxehs-1]
		;cntrd,cimage,boxehs,boxehs,xps,yps,5.
		;xp=xps+xp-boxehs & yp=yps+yp-boxehs
		
		;if rho le 0.4 then cntrd,image,xp,yp,xp,yp,5. ;shift to find max when the inner part of the PSF might be washed out


		;xshift= planet_r[ii] * (1./pxscale) * Cos(planet_theta[ii])
		;yshift= planet_r[ii] * (1./pxscale) * Sin(planet_theta[ii])
		print, 'Xp,Yp = ',xp,yp
		;hak
		xc=xhs & yc=yhs  & xsize=xs & ysize=ys
		print, xhs, yhs

		;measure noise
		nimage=image ;make a copy
		if pxscale gt 0.05 then boxehs=8 else boxehs = 8
		nimage[xp-boxehs:xp+boxehs-1,yp-boxehs:yp+boxehs-1]=!values.f_nan
		
		;BD45 background star
		if pxscale gt 0.05 then begin & xp1=183 & yp1=196 & endif else begin & xp1=218 & yp1=223 & endelse
		nimage[xp1-boxehs:xp1+boxehs-1,yp1-boxehs:yp1+boxehs-1]=!values.f_nan
		if pxscale gt 0.05 then begin & xp1=188 & yp1=191 & endif else begin & xp1=226 & yp1=215 & endelse
		nimage[xp1-boxehs:xp1+boxehs-1,yp1-boxehs:yp1+boxehs-1]=!values.f_nan
		if pxscale gt 0.05 then begin & xp1=193 & yp1=186 & endif else begin & xp1=210 & yp1=231 & endelse
		nimage[xp1-boxehs:xp1+boxehs-1,yp1-boxehs:yp1+boxehs-1]=!values.f_nan
		
		

		;block center (more of a test) ;not needed anymore since this is done radially above
	;	boxehs=12. ;source exclusion half-size
	;	nimage[xhs-boxehs:xhs+boxehs-1,yhs-boxehs:yhs+boxehs-1]=!values.f_nan
		;bkg=mean(nimage,/nan)	
		;noise=stddev(nimage-bkg,/nan)
		
		

		;measure signal as the peak... just a rough first guess
		signal=max(image[xp-1:xp+1,yp-1:yp+1])
		;nimage= image > 0
		;signal=total(nimage[xp-1:xp+1,yp-1:yp+1],/nan)
		;signal=signal/sqrt(9.)

	
		
		
		
			
			nantr=10
			nanmask=where(nimage lt (mean(nimage,/nan)-nantr*stddev(nimage,/nan) ) or nimage gt (mean(nimage,/nan)+nantr*stddev(nimage,/nan) ))
			image[nanmask]=!values.f_nan
			nimage[nanmask]=!values.f_nan
		
		noise=stddev(nimage,/nan)

		
		
		result=signal/noise

	;if rho gt 0.5 then begin ;use apertures
	if rho gt 0. then begin ;use apertures



	fwhm=(4.0E-6)/6.5 *206265 / pxscale
	print, fwhm; & hak

		;measure signal in an aperture
		aper_rad=fwhm/2.;1.      ;fwhm ;& skyradi=4. & skyrado=4.
		zimage=image
		zimage[where(finite(zimage) ne 1)]=0.
		
		aper, zimage, xp,yp, signal,err,sky,skyerr,1,aper_rad,[0,0],[-99E99,99E99],/flux,/exact,/silent,SETSKYVAL = 0,/nan;bkg
		


	;copied from find_sources
	;start at source, progress by each pixel (to 2PiR) and calculate flux in each aperture
	ntheta=float( round( ( 2.*!PI*float(rho/pxscale)/(2.0*aper_rad) )-1 ) ) ;about one angular bin per resolution element




		nzimage=nimage
		nzimage[where(finite(nzimage) ne 1)]=0.
		;writefits,'/Users/kevinwagner/Desktop/nzimage.fits',nzimage
		;writefits,'/Users/kevinwagner/Desktop/zimage.fits',zimage
;ensure at least 4 apertures are used for SNR measurements, if available, otherwise reject those adjacent to the source

;if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta & endif else begin & stitheta=3 & stendtheta=ntheta-2 & endelse ;used this one for 0.3 l/d apertures to skip the source and ADI wings 
;if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta & endif else begin & stitheta=2 & stendtheta=ntheta-1 & endelse ;using this one for l/d apertures to not reject too many

;shouldn't this be:
if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta-1 & endif else begin & stitheta=2 & stendtheta=ntheta-2 & endelse

	for itheta=stitheta,stendtheta do begin
		ttheta=theta+(float(itheta)/float(ntheta))*2.*!PI;+(!PI/31.)
		xtest=xhs+(float(rho/pxscale)*cos(ttheta)) & ytest=yhs+(float(rho/pxscale)*sin(ttheta))
		;print, size(frame)
		aper, nzimage, xtest, ytest, flux,fluxerr,sky,skyerr,1.75,aper_rad,[0,0],[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact,/nan

		if itheta eq stitheta then fluxarr=flux else fluxarr=[fluxarr,flux]
		if itheta eq stitheta then xarr=xtest else xarr=[xarr,xtest]
		if itheta eq stitheta then yarr=ytest else yarr=[yarr,ytest]
	endfor ;;itheta for
		print, fluxarr
		fluxarr[(where(fluxarr eq 0.))]=!values.f_nan
		noise=stddev(fluxarr,/nan)

		signal=signal-mean(fluxarr,/nan)
		noise=noise*sqrt(1.+(1./n_elements(fluxarr[where(finite(fluxarr) ne 0)]))) ;small sample statistics correction: Mawet+2014
			result=signal/noise
		print, 'Source recovered with significance = ',result






	endif;aperture if

		

		if result lt 0. then result=0.5 ;if no signal, double next the contrast guess next time
		;if result eq 0. then begin & result=threshold & contrasts[n] = !values.f_nan & endif ;in blocked area? move on!
		
		;if the brightness has already been modulated 5 times and is still higher, then assume this would be detected
		if runs gt 5 and result ge threshold then result=threshold 
		if runs gt 7 and result lt 0.9*threshold then begin & result=threshold & contrasts[n] = !values.f_nan & endif ;in blocked area? move on!
		print, 'Source recovered with significance = ',result

		;if result is not within 10% of the threshold, revise the injected contrast
		if abs(1.-(result/threshold) ) gt 0.1 then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold))	
		;contrasts=contrasts < contrast_guess
		
		;if result lt 0.95*threshold then contrasts[n:n_elements(separations)-1]=contrast/((result/threshold)) 
		results[n]=result ;store results for diagnostics
		print, contrasts
		print, results
		print, separations
		;hak

		;start plotting things
		cgps_open,root_dir+'contrast/'+obj_id+'_'+filter+'_contrast_curve.ps';,xsize=9.5, ysize=6.5, /inches,/landscape
		;y_range=[1E-8,max([1E-3,2.*max(contrasts,/nan)])]
			;open a new plot if it is the first position angle
			if m eq 0 then cgplot,separations,contrasts,xtitle='Angular Separation (arcsec)',xstyle=8,ytitle='Contrast',/ylog,xrange=x_range,yrange=y_range,charsize=1	
		
			;open a new plot for the first PA and then oplot for the remaining
			if m gt 0 then begin
				cgplot,[0,0],[0,0],xtitle='Angular Separation (arcsec)',ytitle='Contrast',/ylog,xstyle=8,xrange=x_range,yrange=y_range,charsize=1,color='gray'
				;if m gt 1 then for l=1,m-1 do  if plot_ind then cgoplot,separations,curves[*,l],color='gray'

				;if plot_ind then cgoplot, separations, contrasts,color='gray' ;now plot the current contrasts

				if m gt 0 then for l=0,m-1 do  if plot_ind then cgoplot,separations,curves[*,l],color='gray'

				if m eq 1 and plot_ind then cgoplot,separations,curves,color='gray'

				if plot_ind then cgoplot, separations, contrasts,color='gray' ;now plot the current contrasts
			endif

;Procyon is 79 Jy, as estimated via Vizier (mostly constrained by IRAS 11.6 µm point at 79 Jy)
;cgaxis,/YAxis,charsize=1.,yrange=y_range*starJy*1000.,ystyle=1,ytitle='mJy'

cgaxis,/XAxis,charsize=1.,xrange=x_range*1.0*d,xstyle=1,xtitle=strcompress('Projected Separation (au)')
;cgtext,max(x_range)*0.55,max(y_range)*0.48,'SNR = '+string(sigfig(result,3)),charsize=1;+' / '+'Floor = '+string(sigfig(min(contrasts,/nan),3)),charsize=1

cgtext,max(x_range)*0.6,max(y_range)*0.42,strcompress('Test '+string(testn)+' / '+string(n_elements(separations)*n_elements(pas))+' -- SNR = '+string(sigfig(result,3))),charsize=1;,charsize=1


;define some things
;au2sun=215.032 ;au in units of Solar radii
;;RJ2sun=0.10049 ;RJ in units of Solar radii
;mic2ang=1E4 ;1 micron in units of angstroms
;startwave=9.8
;endwave=12.4
;wavestep=0.1 ;making a synthetic box filter 



;albedo=0.3 ;albedo (between 0 and 1; E,N,U,S are about 0.3, jupiter is about 0.5)
;Rpl=[0.,0.09,1.7*0.09,0.4,1.,1.7];0.09;1. ;Jupiter radii (NB: Saturn is 0.9 RJ, Neptune is 0.3 RJ, Earth is 0.09 RJ)
;Apl=[0.,0.31,0.31,0.3,0.5,0.3]

;print, 'At this distance, FoV = 0.2" to 5" (AU) = ',[0.2,5.0]*d


		;plot mean contrast curve from the data reductions

		if m gt 1 then begin
			medcurve=mean(curves,dim=2,/nan)
			cgoplot,separations,medcurve,thick=4
		endif




		cgps_close,/pdf
		
		
		writefits,strcompress(root_dir+'contrast/'+obj_id+'_'+filter+'_contrastcurve_testimage.fits',/rem),image,hd
		writefits,strcompress(root_dir+'contrast/'+obj_id+'_'+filter+'_contrastcurve_testnoise.fits',/rem),nimage,hd
	
	;	hak
	
		hsz=((size(image))(1))/2.0
		
		if keyword_set(display) then cgimage,[[image[hsz-(rho/pxscale)-20:hsz+(rho/pxscale)+19,hsz-(rho/pxscale)-20:hsz+(rho/pxscale)+19],nimage[hsz-(rho/pxscale)-20:hsz+(rho/pxscale)+19,hsz-(rho/pxscale)-20:hsz+(rho/pxscale)+19]]],/keep,stretch=1,clip=4,color=2;,missing_color='white',missing_value=!values.f_nan;,/display
		;if keyword_set(display) and m gt 1 then cgoplot,separations,medcurve
		;if keyword_set(display) and m eq 0 then cgoplot,separations,curves

		if keyword_set(display) then cgtext,0.2,0.05,textoidl(strcompress('SNR(\rho='+string(sigfig(rho,2))+'",\theta='+string(sigfig(theta/!DTOR,3))+'^o,cont.='+string(sigfig(contrast,2))+')='+string(sigfig(result,3)),/rem)),color='white',charsize=1.5,charthick=1.25
		if keyword_set(display) then cgtext,0.42,0.85,obj_id+' - '+filter,color='white',charsize=1.5,charthick=1.25

		;if keyword_set(display) then cgcleanup
		;output images once the threshold has been found
		;if result ge 0.95*threshold and result le 1.05*threshold then $
			;if n eq 0 and m eq 0 then imgs=image else 
			if n eq 0 and m eq 0 then imgs=image
		if result ge (1.0-precision)*threshold and result le (1.0+precision)*threshold then $			
			if n gt 0 or m gt 0 then imgs=[ [[imgs]], [[image]] ]
			;print, n, m & hak			
			;writefits,'/Users/kevinwagner/Desktop/test.fits',imgs		
			if n gt 0 or m gt 0 then writefits,root_dir+'/contrast/'+obj_id+'_'+filter+'_contrastcurve_images.fits',imgs
	endwhile
endfor ;n

;placeholder
curves=[ [curves], [contrasts] ]
print, 'Saving...'
save,filename=strcompress(root_dir+'contrast/'+obj_id+'_'+filter+'_contrast_output.sav',/rem),curves,pas,separations,threshold

endfor ;m



print, 'Completed contrast curve generation in ',(systime(/JULIAN)-contstarttime)*86400./60.,' minutes.'

end
