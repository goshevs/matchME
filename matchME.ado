********************************************************************************
*** matchME
***
***
**  v.0.01
**
** Simo Goshev
**



*** dist   : type of distance measure
*** cutoff : cut-off distance
*** smp    : size of matching pool
*** force  : match beyond reciprocal matches


capture program drop matchME
program define matchME

	syntax varlist, dist(string asis) [cutoff(real 0.0) smp(integer 3) force]

	qui {
		tempname mymat res
		
		*** Compute the distance matrix
		
		if "`dist'" == "Euclidean" {
			matrix dis `mymat' = `varlist', `dist'	
		}
		else if "`dist'" == "nEuclidean" {
			*** ## Normalized Euclidean distance (non-scale preserving)
			mahascores `varlist', genmat(`mymat') euclidean compute_invcovarmat
		}
		else if "`dist'" == "Mahalanobis" {	
			*** ## Mahalanobis distance
			mahascores `varlist', genmat(`mymat') compute_invcovarmat
		}
		else {
			di in r "Distance measure is not supported"
			exit 489
		}
		
		*** Find first, second and third order matches
		if `smp' == 0 {
			di in r "The size of the matching pool must be strictly positive"
			exit 489
		}
		mata: findMatch("`mymat'", "`res'", `smp')
		
		* mat list `res'
		
		*** Merge in the results from the matching

		tempfile mymatches

		preserve
		clear
		
		svmat `res'
		
		ren `res'1 id
		
		if `smp' > 1 {
			forval i = 1/`smp' {
				ren `res'`=`i'*2' `dist'`i'm
				ren `res'`=`i'*2 + 1' `dist'`i'd
			}
		}
		save `mymatches', replace
		
		restore
		
		merge 1:1 id using `mymatches', nogen
		
		sort id
		gen `dist'Matched = ""
		
		*** Apply conditional cutoff
		if `cutoff' > 0.0 {
			qui forval i = 1/`=_N' {
				local check = `dist'Matched[`i']
				if "`check'" == "" {
					if `dist'1d[`i'] > `cutoff' {
						replace `dist'Matched = "OR" if _n == `i'
					}
				}
			}
		}
			
		gen eligible = ""
		*** Find the Best Match
		
		*** Create the initial set of eligible matches
		forval i =1/`=_N' {
			local check = `dist'Matched[`i']
			if "`check'" == "" {
				local mypool ""
				forval j = 1/`smp' {
					if `cutoff' > 0.0 {
						local mypool "`mypool' `=cond(`dist'`j'd[`i'] <= `cutoff',`dist'`j'm[`i'],.)'"
					}
					else {
						local mypool "`mypool' `=`dist'`j'm[`i']'"
					}
				}
				replace eligible = itrim(trim("`mypool'"))  if _n == `i'
			}
		}
		*** Remove non-eligible (coded as missing) from eligible string
		replace eligible = trim(itrim(subinstr(eligible, ".", "", .)))
		
		
		*** Make all RECIPROCAL MATCHES by order
		*** ------->>> Reciprocal matching loop
		
		forval l = 1/`smp' {
			forval i = 1/`=_N' {
				local check = `dist'Matched[`i']
				if "`check'" == "" {
					local elig "`=eligible[`i']'"
					local bestMatch: word 1 of `elig' 
					local soM "`i' `bestMatch'"
					local soM: list sort soM
					replace `dist'Matched = "`soM'" if _n == `i'
				}
			}
			*** Find the pairs
			bys `dist'Matched: gen ind = _N if `dist'Matched ~= "" 
			replace `dist'Matched = "" if ind ~= 2 & `dist'Matched ~= "OR"
			drop ind
			replace eligible  = "" if `dist'Matched ~= ""
			
			*** Update eligible

			*** Not matched ids:
			* levelsof id if `dist'Matched == "", local(`dist'NM)
			*** Matched ids:
			levelsof id if `dist'Matched != "" & `dist'Matched != "OR", local(`dist'M)
			
			qui forval i = 1/`=_N' {
				local check = `dist'Matched[`i']
				if "`check'" == "" {
					local matches "`=eligible[`i']'"
					local intersection: list matches & `dist'M
					local matches: list matches - intersection
					replace eligible = "`matches'" if _n == `i'
				}
			}
			
			replace `dist'Matched = "PNM" if `dist'Matched == "" & eligible == ""	
			sort id
		}
		
		*** -------<<< Reciprocal matching loop
		

		*** Forcing the remaining matches
		if "`force'" ~= "" { 
		
			preserve
			
			*** Expand the dataset
			gen nexpand = 1
			forval i = 1/`=_N' {
				local check = `dist'Matched[`i']
				if "`check'" == "" {
					local elig "`=eligible[`i']'"
					loc nElig: word count `elig'
					replace nexpand = `nElig' if _n ==`i'
				}
			}
			
			expand nexpand
			replace nexpand =0 if nexpand == 1
			sort `dist'Matched id

			bys id: replace nexpand = _n if _n ~= _N
			forval i = 1/`=_N'  {
				if nexpand[`i'] > 0 {
					local cElig "`=eligible[`i']'"
					* di "`cElig'"
					local myVal: word `=nexpand[`i']' of `cElig'
					replace eligible = "`myVal'" if _n == `i'
				}
			}
			
			gen dist =.
			forval i = 1/`=_N' {
				local check = `dist'Matched[`i']
				if "`check'" == "" {
					local elig "`=eligible[`i']'"
					local bestMatch: word 1 of `elig'
					forval j = 1/`smp' {
						if `dist'`j'm[`i'] == `bestMatch' {
							replace dist = `dist'`j'd[`i'] if _n == `i'
						}
					}
				}
			}
			
			drop if `dist'Matched ~= "" & `dist'Matched ~= "PNM"
			
			levelsof id if `dist'Matched == "", local(idLevs)
			
			clonevar elig1 = eligible
			destring elig1, replace
			levelsof elig1, local(elLevs)
			drop elig1
			
			*local nIter: list idLevs & elLevs
			*local nIter: word count `nIter'
			local myComb "`idLevs' `elLevs'"
			local nIter: list uniq myComb
			local nIter: word count `nIter'
			
			sort dist
			
			di "Number of interations: `nIter'"
			
			forval i = 1/`nIter' { 
				*** Match
				local omatch "`=eligible[1]'"
				* levelsof id if id == `omatch'
				local match "`=id[1]' `omatch'"
				local match: list sort match
				replace `dist'Matched = "`match'" if _n==1 | id == `omatch'
				
				duplicates drop id `dist'Matched if id == `omatch', force

				*** Update eligible and id
				replace eligible = regexr(eligible, "^(`=subinstr("`match'", " ", "|", .)')$", "") if _n ~= 1
				gen idStr = string(id)
				gen flag = 1 if regexm(idStr, "^(`=subinstr("`match'", " ", "|",.)')$") & _n ~= 1 & `dist'Matched == ""
				
				drop if flag == 1
				drop idStr flag
			
				*** If id apprears in id once --> label as PNM; if appears multiple times --> keep all with non-empty eligible
				levelsof id if eligible ~= "" , local(checkDups)
				levelsof id if `dist'Matched == "" & eligible == "" , local(iter)
				
				foreach val of local iter {
					local inList: list val & checkDups
					if "`inList'" ~= "" {
						* noi di "`inList'"
						drop if id == `val' & eligible == ""
					}
					else {
						replace `dist'Matched = "PNM" if `dist'Matched == "" & eligible == "" & id == `val'
					}
				}
				replace dist = . if _n == 1 | eligible == ""
				
				duplicates drop id `dist'Matched if eligible == "", force
				
				sort dist
				
				*** Check whether need to continue matching
				count if eligible == ""
				if "`r(N)'" == "`=_N'" {
					continue, break
				}
			
			}
			sort id

			drop eligible nexpand dist
			
			keep id `dist'Matched
			
			tempfile mergebackin
			
			save "`mergebackin'", replace
			restore
			merge 1:1 id using "`mergebackin'", nogen update replace
		}

		replace `dist'Matched = "NM" if `dist'Matched == "PNM" | `dist'Matched == ""
		drop eligible
	}
	noi di in g "Matching completed."

end


cap mata: mata drop findMatch()
mata:
void function findMatch(string scalar mymat, string scalar outmat, real scalar smp) {
	
	real matrix res
	real vector ordered
	
	scores = st_matrix(mymat)
	vlen = rows(scores)
	_diag(scores, J(1, vlen, 999))
	
	res = J(vlen, smp*2 + 1,.) 
	
	for(i=1;i<=vlen;i++) {
		ordered = sort(scores[i,1...]',1)'
		minvals = ordered[1,1..smp]
		res[i,1] = i
		for (j=1;j<=vlen;j++) {
			for (k=1;k<=smp;k++) {
				if (scores[i,j] == minvals[1,k]) {
					res[i,k*2..k*2+1] = (j, minvals[1,k])
				}	
			}
		}
	}
	st_matrix(outmat, res)
}
end
