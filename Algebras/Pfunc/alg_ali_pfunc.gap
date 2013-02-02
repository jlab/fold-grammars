//we can't use a two part partition function (one component for energy, the other for covariation), because the scorings for covariation and the Boltzman function have a very unintuitive behaviour if combined --> very strange results for stochastical backtracing. Better directly fuse energy and covariation into one double!
algebra alg_ali_pfunc implements sig_foldrna(alphabet = M_Char, answer = double) {
  double sadd(Subsequence lb, double x) {
    double res;
	
	float sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
	return x * scale(1) * mk_pf(sbase_sum / float(rows(lb)));
  }
  double cadd(double x, double y) {
    return x * y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
	Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	
	return x * scale(1) * mk_pf((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	  
	return x * scale(1) * mk_pf((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	  
	return x * scale(2) * mk_pf((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x * mk_pf(termau_energy(lb, rb) / float(rows(lb)));
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
	return x * scale(2) * mk_pf(int(sr_energy(lb, rb) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
	return scale(2+r.j-r.i) * mk_pf(int(hl_energy(r) / float(rows(r))) + covscore(lb, lb.i, rb.i));
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
	return x * scale(2+lr.j-lr.i) * mk_pf(bl_energy(lr, rb) / float(rows(lb)) + covscore(lb, lb.i, rb.i));
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x * scale(2+rr.j-rr.i) * mk_pf(br_energy(lb, rr) / float(rows(lb)) + covscore(lb, lb.i, rb.i));
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
	return x * scale(2+lr.j-lr.i+rr.j-rr.i) * mk_pf(il_energy(lr, rr) / float(rows(lr)) + covscore(lb, lb.i, rb.i));
	//~ if (lr.j-lr.i + rr.j-rr.i > 30) { res.mfe = 99999; } // ugly hack to realize a filter that rejects internal loops whose combined unpaired loop regions exeed 30 bases. Grammar filter causes errors with --kbacktrace. Georg and I don't know why.
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
    return x * scale(4) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x * scale(2) * mk_pf(ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
  }
  double incl(double x) {
    return x * mk_pf(ul_energy());
  }
  double addss(double x, Subsequence r) {
    return x * scale(r.j-r.i) * mk_pf(ss_energy(r) / float(rows(r)));
  }
  double nil(Subsequence n) {
    return 1.0;
  }
  choice [double] h([double] i) {
    return list(sum(i));
    //~ return i;
  }
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {double x; return x;}
  double combine(double le,double re) {double x; return x;}
  double trafo(double e) {double x; return x;}
  double ssadd(Subsequence lb,double e) {double x; return x;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {double x; return x;}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double ambd_Pr(double le,Subsequence b,double re) {double x; return x;}
  double ambd(double le,Subsequence b,double re) {double x; return x;}
  double cadd_Pr_Pr_Pr(double le,double re) {double x; return x;}
  double cadd_Pr_Pr(double le,double re) {double x; return x;}
  double cadd_Pr(double le,double re) {double x; return x;}
}

algebra alg_ali_pfunc_id extends alg_ali_pfunc {
  choice [double] h([double] i) {
    return i;
  }
}

