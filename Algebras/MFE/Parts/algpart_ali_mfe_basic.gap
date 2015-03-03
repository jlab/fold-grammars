  mfecovar sadd(Subsequence lb, mfecovar x) {
    mfecovar res;
	
	int sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
	res.mfe = x.mfe + (sbase_sum / float(rows(lb)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar cadd(mfecovar x, mfecovar y) {
	  mfecovar res;
	  res.mfe = x.mfe + y.mfe;
	  res.covar = x.covar + y.covar;
	  return res;
  }
  mfecovar edl(Subsequence ldangle, mfecovar x, Subsequence rb) {
    mfecovar res;
    
	Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar edr(Subsequence lb, mfecovar x, Subsequence rdangle) {
    mfecovar res;
    
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar edlr(Subsequence ldangle, mfecovar x, Subsequence rdangle) {
    mfecovar res;
    
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar drem(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res;
    
    res.mfe = x.mfe + (termau_energy(lb, rb) / float(rows(lb)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar dall(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar;
	return res;
  }
  mfecovar sr(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + (sr_energy(lb, rb) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar hl(Subsequence lb, Subsequence r, Subsequence rb) {
    mfecovar res;
    
	res.mfe = (hl_energy(r) / float(rows(r)));
	res.covar = covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar bl(Subsequence lb, Subsequence lr, mfecovar x, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + (bl_energy(lr, rb) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar br(Subsequence lb, mfecovar x, Subsequence rr, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + (br_energy(lb, rr) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar il(Subsequence lb, Subsequence lr, mfecovar x, Subsequence rr, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + (il_energy(lr, rr) / float(rows(lr)));
	//~ if (lr.j-lr.i + rr.j-rr.i > 30) { res.mfe = 99999; } // ugly hack to realize a filter that rejects internal loops whose combined unpaired loop regions exeed 30 bases. Grammar filter causes errors with --kbacktrace. Georg and I don't know why.
    res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar mldl(Subsequence lb, Subsequence dl, mfecovar x, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar mldr(Subsequence lb, mfecovar x, Subsequence dr, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar mldlr(Subsequence lb, Subsequence dl, mfecovar x, Subsequence dr, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar ml(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res;
    
	res.mfe = x.mfe + ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb)));
	res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovar mlall(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar + covscore(lb, lb.i, rb.i);
    return res;
  }
  mfecovar incl(mfecovar x) {
    mfecovar res;
    
	res.mfe = x.mfe + ul_energy();
	res.covar = x.covar;
	
    return res;
  }
  mfecovar addss(mfecovar x, Subsequence r) {
    mfecovar res;
    
	res.mfe = x.mfe + (ss_energy(r) / float(rows(r)));
	res.covar = x.covar;
	
    return res;
  }
  mfecovar nil(Subsequence n) {
	mfecovar res;
	res.mfe = 0;
	res.covar = 0;
    return res;
  }
  choice [mfecovar] h([mfecovar] i) {
    return list(minimum(i));
  }
  
