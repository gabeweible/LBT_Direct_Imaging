pro coadd_ct_snr
COMPILE_OPT IDL2

; Trying this out...
SET_PLOT, 'Z'

newline = string(10B)

;------------------------------[ Start User Input ]---------------------------------

half_cropped_sz = 190.;250. ; px (used to be 250, i.e., 500 px by 500 px cropped images)

; General/Combine Parameters
obj = 'Alcor'
band = 'M'
raw_path = '/Volumes/T7/alcor_raw'
obj_path = '/Users/gabeweible/OneDrive/research/Alcor/'

pos = [189, 196] ; approx. position of Alcor B

fwhm = 13.0773; Average of dith1 and dith2 fwhm in x and in y (from pupil medians)

filter = 0; high-pass filter kernel size for ADI and KLIP; or 0 for no filter
neg_inj=0
uncert=0

; 20 APPEARS TO BE THE SMALLEST COADD THAT WORKS ON MY MACHINE
coadds = reverse(indgen(1, start=20, increment=5)); reverse to do smaller datacube to larger
corr_threshes = findgen(20, start=0.9900, increment=0.0005)
SNR_list = list()

foreach coadd, coadds do begin; loop through coadds

	; get the output path for this coadd
	output_path = obj_path + 'macbook_' + strcompress(coadd, /r)
	if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)
	output_path += '/'

	ct_i = 0
	foreach corr_thresh, corr_threshes do begin

		; helps find the right left_adi FITS file
		adi_suffix = output_path + 'combined/' + obj_name + 'ct_' + string(corr_thresh)$
			+ 'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj) + '_uncert_' +$
			string(uncert)

		; Read in the left_adi FITS file for this ct and coadd
		image = readfits(strcompress(adi_suffix + '_left_adi.fits',/r))

		; remove zeros
		image[where(image eq 0)] = !values.f_nan

		;x- and y-sizes and half-sizes
		xs = (size(image))[1] & ys = (size(image))[2]
		xhs = xs/2. & yhs = ys/2. ;half-sizes
		print, xhs, yhs

		xc=xhs & yc=yhs  & xsize=xs & ysize=ys

		; NOISE IMAGE SETUP
		nimage = image ;make a copy

		;exclude known companion for the 'noise' image
		boxehs=13. ;source exclusion half-size (approx. fwhm)
		nimage[pos[0]-boxehs:pos[0]+boxehs-1, pos[1]-boxehs:pos[1]+boxehs-1] = !values.f_nan

		; get rid of really high or low values (outside 10 sigma from the mean)
		nantr=10
		nanmask=where(nimage lt (mean(nimage,/nan)-nantr*stddev(nimage,/nan) ) or nimage gt (mean(nimage,/nan)+nantr*stddev(nimage,/nan) ))
		image[nanmask]=!values.f_nan
		nimage[nanmask]=!values.f_nan

		; full width at half maximum in px
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
		SNR=signal/noise
		print, 'Source recovered with SNR = ',SNR
		print, 'for coadd = '+string(coadd)+' and ct = '+string(corr)
		SNR_list.Add, SNR

	endforeach; corr_threshes
endforeach; coadds

SNR_array = SNR_list.toArray(/NO_COPY)
save, filename=strcompress(obj_path+'coadd_ct_snr_test1.sav',/rem), SNR_array

end