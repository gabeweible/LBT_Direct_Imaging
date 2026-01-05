; Function to extract nod number for proper sorting
function extract_nod_number, filename
	; Extract the nod number from filenames like '_nod05_' or '_nod10_'
	pos = strpos(filename, '_nod')
	if pos eq -1 then return, -1
	
	; Find the end of the nod number (next underscore)
	start_pos = pos + 4  ; Skip '_nod'
	end_pos = strpos(filename, '_', start_pos)
	if end_pos eq -1 then return, -1
	
	nod_str = strmid(filename, start_pos, end_pos - start_pos)
	return, fix(nod_str)
end