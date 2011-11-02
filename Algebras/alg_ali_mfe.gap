algebra alg_ali_mfe implements sig_foldrna(alphabet = M_Char, answer = mfecovar) {
  mfecovar sadd(Subsequence lb, mfecovar x) {
    mfecovar res = x;
	
	int sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
	res.mfe = res.mfe + (sbase_sum / float(rows(lb)));
	
    return res;
  }
  mfecovar cadd(mfecovar x, mfecovar y) {
    return (x + y);
  }
  mfecovar edl(Subsequence ldangle, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
	Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	res.mfe = res.mfe + ((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));

    return res;
  }
  mfecovar edr(Subsequence lb, mfecovar x, Subsequence rdangle) {
    mfecovar res = x;
    
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.mfe = res.mfe + ((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));

    return res;
  }
  mfecovar edlr(Subsequence ldangle, mfecovar x, Subsequence rdangle) {
    mfecovar res = x;
    
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.mfe = res.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));

    return res;
  }
  mfecovar drem(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
    res.mfe = res.mfe + (termau_energy(lb, rb) / float(rows(lb)));

    return res;
  }
  mfecovar sr(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + (sr_energy(lb, rb) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar hl(Subsequence lb, Subsequence r, Subsequence rb) {
    mfecovar res;
    
	res.mfe = (hl_energy(r) / float(rows(r)));
	res.covar = covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar bl(Subsequence lb, Subsequence lr, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + (bl_energy(lr, rb) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar br(Subsequence lb, mfecovar x, Subsequence rr, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + (br_energy(lb, rr) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar il(Subsequence lb, Subsequence lr, mfecovar x, Subsequence rr, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + (il_energy(lr, rr) / float(rows(lr)));
	if (lr.j-lr.i + rr.j-rr.i > 30) { res.mfe = 99999; } // ugly hack to realize a filter that rejects internal loops whose combined unpaired loop regions exeed 30 bases. Grammar filter causes errors with --kbacktrace. Georg and I don't know why.
    res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar mldl(Subsequence lb, Subsequence dl, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar mldr(Subsequence lb, mfecovar x, Subsequence dr, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar mldlr(Subsequence lb, Subsequence dl, mfecovar x, Subsequence dr, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar ml(Subsequence lb, mfecovar x, Subsequence rb) {
    mfecovar res = x;
    
	res.mfe = res.mfe + ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);

    return res;
  }
  mfecovar incl(mfecovar x) {
    mfecovar res = x;
    
	res.mfe = res.mfe + ul_energy();

    return res;
  }
  mfecovar addss(mfecovar x, Subsequence r) {
    mfecovar res = x;
    
	res.mfe = res.mfe + (ss_energy(r) / float(rows(r)));

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
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  mfecovar acomb(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar combine(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar trafo(mfecovar e) {mfecovar x; return x;}
  mfecovar ssadd(Subsequence lb,mfecovar e) {mfecovar x; return x;}
  mfecovar mladl(Subsequence lb,Subsequence dl,mfecovar e,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladldr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mldladr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladlr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladr(Subsequence lb,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar ambd_Pr(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar ambd(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
}
