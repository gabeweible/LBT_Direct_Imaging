pro cc_hii1348_v2, klip=klip, k_klip=k_klip, coadd=coadd, extra=extra, half_pas=half_pas, third_pas=third_pas
; CURRENT CODE (01/29/2025)
COMPILE_OPT IDL2

; Trying this out...
SET_PLOT, 'Z'

newline = string(10B)

;USER INPUTS
; Going being 2.0 and with more increments to see if that improves our results
; at all (though this will take longer to run, to be sure)
separations=[0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35,$
				 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]; Arcsec

; Position angles (sorta), 18 of them all separated by 20 deg
; starting at 5 deg
pas = findgen(18)*20. + 5.

if half_pas eq 1 then remove, [0, 2, 4, 6, 8, 10, 12, 14, 16], pas; this should do every other?
if third_pas eq 1 then remove, [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17], pas; this should do every 2/third?
print, 'Number of PAs to test:', n_elements(pas)

threshold = 5.; Sigma/Significant/SNR
pxscale = 0.010648; Plate scale for our observations [arcsec/px]
contrast_guess = 0.0035 ;guess at the starting point for the contrast curve
pos = [277.3, 353.6] ; approx. position of HII 1348 b

; Where to find our files and put the results
;output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_25/';'/Users/gabe/reduction/macbook_25/'
output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/HII1348/macbook_'+$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)
output_path += '/'

; Rotate/KLIP parameters (needed to read in our total_klip file)
if keyword_set(klip) then klip=[klip]
bin = 3
bin_type = 'mean'
combine_type = 'nwadi'
k_klip = k_klip
angsep= 1.
anglemax = 360.
nrings = 4.
n_ang = 2

;-------------------- END USER INPUT -------------------------------

;save start time
contstarttime=systime(/JULIAN)

;set up arrays
results=separations & replicate_inplace, results, 0.
contrasts=separations & replicate_inplace, contrasts, contrast_guess
curves=[]

;start a loop over position angle
runs = 0
for m = 0,n_elements(pas)-1 do begin

PA=separations & replicate_inplace, PA, pas[m]*!DTOR ;set all same PA for now

;start a loop over separation
for n = 0,n_elements(separations)-1 do begin
	print, 'Working on test separation',n,'/',n_elements(separations)-1
	result=99E99 ;start significance off at infinity
	
	; Loop until we get within 5% of 5-sigma
	while result gt 1.05*threshold or result lt 0.95*threshold do begin
	
		rho = separations[n] & theta = PA[n] & contrast = contrasts[n]
		if n gt 0 and finite(contrast) eq 0 then contrast=contrasts[n-1]
		print, 'Testing rho, theta, contrast = ', rho, theta/!DTOR, contrast
		print, 'Run:', runs
		runs += 1

		print, 'running HII 1348B pipeline for k_klip = ', k_klip
		hii1348_pipeline, rho=rho, theta=theta, contrast=contrast, pre_inj=0, neg_inj=0, coadd=25, uncert=0, klip=klip,$
			use_gauss=0, extra='kklip_'+strcompress(k_klip, /r), k_klip=k_klip,$
			nod='total', adi=0
		
		;read image
		if klip eq 1 then begin ; KLIP	
			
			image=readfits(strcompress(output_path + 'combined/' + 'HII1348' +$
		 		'_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+$
		 		combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
		 		string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+$
		 		'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) +$
		  		'_neg_inj_0' + '_total_klip.fits', /rem))
		  		
		endif else begin; ADI
		
			image = readfits(strcompress(output_path + 'combined/' +$
				'HII1348ct_0.994000filt_0.00000_neg_inj_0_uncert_0_total_adi.fits', /rem))

		endelse
		; remove zeros
		image[where(image eq 0)] = !values.f_nan

		;what is rho, theta on the image?
		xs = (size(image))[1] & ys = (size(image))[2]
		xhs = xs/2. & yhs = ys/2. ;half-sizes
		print, xhs, yhs
		print, rho
		xp=xhs + ( (rho/pxscale) * Cos(theta)  ) & yp=yhs + ( (rho/pxscale) * Sin(theta)  )
		print, 'Xp,Yp = ',xp,yp
		
		xc=xhs & yc=yhs  & xsize=xs & ysize=ys
		print, xhs, yhs
		
		; NOISE IMAGE SETUP
		nimage = image ;make a copy
		boxehs = 10. ;source exclusion half-size
		; Exclude injected source for the 'noise' image
		nimage[xp-boxehs:xp+boxehs-1, yp-boxehs:yp+boxehs-1] = !values.f_nan

		;exclude known companion for the 'noise' image
		boxehs=20. ;source exclusion half-size
		nimage[pos[0]-boxehs:pos[0]+boxehs-1, pos[1]-boxehs:pos[1]+boxehs-1] = !values.f_nan
		
		; get rid of really high or low values (outside 10 sigma from the mean)
		nantr=10
		nanmask=where(nimage lt (mean(nimage,/nan)-nantr*stddev(nimage,/nan) ) or nimage gt (mean(nimage,/nan)+nantr*stddev(nimage,/nan) ))
		image[nanmask]=!values.f_nan
		nimage[nanmask]=!values.f_nan
		
		; full width at half maximum in px
		fwhm=9.6; hard-coded, could fit a gaussian to the median frame
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

		if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta-1 & endif else begin & stitheta=2 & stendtheta=ntheta-2 & endelse

		for itheta=stitheta,stendtheta do begin
			ttheta=theta+(float(itheta)/float(ntheta))*2.*!PI;+(!PI/31.)
			xtest=xhs+(float(rho/pxscale)*cos(ttheta)) & ytest=yhs+(float(rho/pxscale)*sin(ttheta))

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

		; Changing the gt .1 to .2 to see if that fixes my nan contrast problem...
		; Spoiler! It doesn't...but I'm debugging this...
		; Maybe gt 1.0? idk
		; 2.0?
		; Putting it to 5.0 just so that I (probably) won't have to worry about it anymore...
		;if runs gt 5 and abs(1.-(result/threshold) ) gt 5.0 then begin & print, 'abs:',$
		;	abs(1.-(result/threshold) ) & result=threshold & contrasts[n] = !values.f_nan &$
		;	print, "Hey dude, I'm setting the contrast to NaN now..." & nan_problem = 1 &$
		;	stop & endif

		if runs gt 5 and abs(1.-(result/threshold) ) le 0.1 then result=threshold 
		
		;if result is not within 5% of the threshold, revise the injected contrast
		if abs(1.-(result/threshold) ) gt 0.05 then contrasts[n:n_elements(separations)-1]=$
			contrast/((result/threshold))	

		results[n]=result ;store results for diagnostics
		print, contrasts
		print, results
		print, separations

		writefits,strcompress(output_path + 'HII1348_contrastcurve_testimage.fits',/rem),image
		writefits,strcompress(output_path + 'HII1348_contrastcurve_testnoise.fits',/rem),nimage
		
		;output images once the threshold has been found
		if result ge 0.95*threshold and result le 1.05*threshold then $
			if n eq 0 and m eq 0 then imgs=image else imgs=[ [[imgs]], [[image]] ]
			
		if n gt 0 or m gt 0 then writefits,strcompress(output_path +$
			'HII1348_contrastcurve_images.fits',/rem),imgs
			
	endwhile
endfor ;n/separations

curves=[ [curves], [contrasts] ]

save,filename=strcompress(output_path + 'HII1348_contrast_output.sav',/rem),$
	curves,pas,separations,threshold
; Write curves to a CSV.
WRITE_CSV, strcompress(output_path + 'HII1348_contrast_output.csv',/rem), curves

endfor ;m/PAs

print, 'Completed contrast curve generation in ',(systime(/JULIAN)-contstarttime)*$
	86400./60.,' minutes.'

end