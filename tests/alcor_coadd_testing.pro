pro alcor_coadd_testing

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

; Allowed coadds that, I think, actually can run on my machine
; descending so that the coadds less likely process (and longer to process)
; are run last.

; 20 APPEARS TO BE THE SMALLEST COADD THAT WORKS ON MY MACHINE
coadds = reverse(indgen(1, start=20, increment=5)); reverse to do smaller datacube to larger
corr_threshes = findgen(20, start=0.9900, increment=0.0005)

foreach coadd, coadds do begin

ct_i = 0
foreach corr_thresh, corr_threshes do begin

	if ct_i eq 0 then begin

		alcor_pipeline, neg_inj=0, pre_inj=1, coadd=coadd,$
			uncert=0, klip=0, fs=1, nod='sx_only', cds='endcds',$
			corr_thresh=corr_thresh, pre_clean=1, outpath='/Volumes/T7/macbook_'+$
			strcompress(coadd, /r)+'/'
			
		ct_i += 1
		
	endif else begin; first correlation threshold at this coadd
		
		alcor_pipeline, neg_inj=0, pre_inj=1, coadd=coadd,$
			uncert=0, klip=0, fs=1, nod='sx_only', cds='endcds',$
			corr_thresh=corr_thresh, pre_clean=0, outpath='/Volumes/T7/macbook_'+$
			strcompress(coadd, /r)+'/'
		
		ct_i += 1
		
	endelse; other cleaning correlation thresholds at this coadd.
		
endforeach; corr_threshes
endforeach; coadds
  
end