algebra alg_ali_pfunc implements sig_foldrna(alphabet = M_Char, answer = answer_ali_pfunc) {
  answer_ali_pfunc sadd(Subsequence lb, answer_ali_pfunc x) {
    answer_ali_pfunc res;
	
	float sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
	res.pfunc = scale(1) * x.pfunc * mk_pf(sbase_sum / float(rows(lb)));
	res.covar = scale(1) * x.covar;
	
    return res;
  }
  answer_ali_pfunc cadd(answer_ali_pfunc x, answer_ali_pfunc y) {
	answer_ali_pfunc res;
	  
	res.pfunc = x.pfunc * y.pfunc;
	res.covar = x.covar * y.covar;
	  
    return res;
  }
  answer_ali_pfunc edl(Subsequence ldangle, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
	  
	Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	res.pfunc = x.pfunc * scale(1) * mk_pf((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));
	res.covar = x.covar * scale(1);
	  
    return res;
  }
  answer_ali_pfunc edr(Subsequence lb, answer_ali_pfunc x, Subsequence rdangle) {
    answer_ali_pfunc res;
	  
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.pfunc = x.pfunc * scale(1) * mk_pf((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));
	res.covar = x.covar * scale(1);
	  
    return res;
  }
  answer_ali_pfunc edlr(Subsequence ldangle, answer_ali_pfunc x, Subsequence rdangle) {
    answer_ali_pfunc res;
    
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	res.pfunc = x.pfunc * scale(2) * mk_pf((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));
	res.covar = x.covar * scale(2);
	  
    return res;
  }
  answer_ali_pfunc drem(Subsequence lb, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
    
    res.pfunc = x.pfunc * mk_pf(termau_energy(lb, rb) / float(rows(lb)));
	res.covar = x.covar;
	  
    return res;
  }
  answer_ali_pfunc sr(Subsequence lb, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(2) * mk_pf(int(sr_energy(lb, rb) / float(rows(lb))));
	res.covar = x.covar * scale(2) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc hl(Subsequence lb, Subsequence r, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = scale(2+r.j-r.i) * mk_pf(int(hl_energy(r) / float(rows(r))));
	res.covar = scale(2+r.j-r.i) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc bl(Subsequence lb, Subsequence lr, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(2+lr.j-lr.i) * mk_pf(bl_energy(lr, rb) / float(rows(lb)));
	res.covar = x.covar * scale(2+lr.j-lr.i) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc br(Subsequence lb, answer_ali_pfunc x, Subsequence rr, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(2+rr.j-rr.i) * mk_pf(br_energy(lb, rr) / float(rows(lb)));
	res.covar = x.covar * scale(2+rr.j-rr.i) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc il(Subsequence lb, Subsequence lr, answer_ali_pfunc x, Subsequence rr, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(2+lr.j-lr.i+rr.j-rr.i) * mk_pf(il_energy(lr, rr) / float(rows(lr)));
	//~ if (lr.j-lr.i + rr.j-rr.i > 30) { res.mfe = 99999; } // ugly hack to realize a filter that rejects internal loops whose combined unpaired loop regions exeed 30 bases. Grammar filter causes errors with --kbacktrace. Georg and I don't know why.
    res.covar = x.covar * scale(2+lr.j-lr.i+rr.j-rr.i) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc mldl(Subsequence lb, Subsequence dl, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb))));
	res.covar = x.covar * scale(3) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc mldr(Subsequence lb, answer_ali_pfunc x, Subsequence dr, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb))));
	res.covar = x.covar * scale(3) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc mldlr(Subsequence lb, Subsequence dl, answer_ali_pfunc x, Subsequence dr, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(4) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb))));
	res.covar = x.covar * scale(4) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc ml(Subsequence lb, answer_ali_pfunc x, Subsequence rb) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(2) * mk_pf(ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb))));
	res.covar = x.covar * scale(2) * mk_pf(covscore(lb, lb.i, rb.i));

    return res;
  }
  answer_ali_pfunc incl(answer_ali_pfunc x) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * mk_pf(ul_energy());
	res.covar = x.covar;
	
    return res;
  }
  answer_ali_pfunc addss(answer_ali_pfunc x, Subsequence r) {
    answer_ali_pfunc res;
    
	res.pfunc = x.pfunc * scale(r.j-r.i) * mk_pf(ss_energy(r) / float(rows(r)));
	res.covar = x.covar * scale(r.j-r.i);
	
    return res;
  }
  answer_ali_pfunc nil(Subsequence n) {
	answer_ali_pfunc res;

	res.pfunc = 1;
	res.covar = 1;

    return res;
  }
  choice [answer_ali_pfunc] h([answer_ali_pfunc] i) {
    return list(sum(i));
    //~ return i;
  }
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  answer_ali_pfunc acomb(answer_ali_pfunc le,Subsequence b,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc combine(answer_ali_pfunc le,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc trafo(answer_ali_pfunc e) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc ssadd(Subsequence lb,answer_ali_pfunc e) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc mladl(Subsequence lb,Subsequence dl,answer_ali_pfunc e,Subsequence rb) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc mladldr(Subsequence lb,Subsequence dl,answer_ali_pfunc e,Subsequence dr,Subsequence rb) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc mldladr(Subsequence lb,Subsequence dl,answer_ali_pfunc e,Subsequence dr,Subsequence rb) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc mladlr(Subsequence lb,Subsequence dl,answer_ali_pfunc e,Subsequence dr,Subsequence rb) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc mladr(Subsequence lb,answer_ali_pfunc e,Subsequence dr,Subsequence rb) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc ambd_Pr(answer_ali_pfunc le,Subsequence b,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc ambd(answer_ali_pfunc le,Subsequence b,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc cadd_Pr_Pr_Pr(answer_ali_pfunc le,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc cadd_Pr_Pr(answer_ali_pfunc le,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
  answer_ali_pfunc cadd_Pr(answer_ali_pfunc le,answer_ali_pfunc re) {answer_ali_pfunc x; return x;}
}

algebra alg_ali_pfunc_id extends alg_ali_pfunc {
  choice [answer_ali_pfunc] h([answer_ali_pfunc] i) {
    return i;
  }
}

