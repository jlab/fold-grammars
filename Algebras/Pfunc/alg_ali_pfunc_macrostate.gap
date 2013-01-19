algebra alg_ali_pfunc_macrostate implements sig_foldrna(alphabet = M_Char, answer = answer_ali_pfunc_macrostate) {
  answer_ali_pfunc_macrostate sadd(Subsequence lb,answer_ali_pfunc_macrostate e) {
	answer_ali_pfunc_macrostate res;
	  
	res.firststem = e.firststem;
	
	float sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
    
    res.pfunc.q1 = scale(1) * e.pfunc.q1 * mk_pf(sbase_sum / float(rows(lb)));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
	
	res.covar.q1 = scale(1) * e.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
	
    return res;
  }

  answer_ali_pfunc_macrostate cadd(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pfunc.q1 = le.pfunc.q1 * re.pfunc.q1;
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = le.covar.q1 * re.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pfunc.q1 = le.pfunc.q1 * sum_elems(re.pfunc);
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = le.covar.q1 * sum_elems(re.covar);
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    	
	res.firststem = le.firststem;
    res.pfunc = mk_tuple(le.firststem, le.pfunc.q1 * re.pfunc.q1);
    res.covar = mk_tuple(le.firststem, le.covar.q1 * re.covar.q1);
    
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr_Pr_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;
    res.pfunc = mk_tuple(le.firststem, le.pfunc.q1 * sum_elems(re.pfunc));
    res.covar = mk_tuple(le.firststem, le.covar.q1 * sum_elems(re.covar));
    
    return res;
  }

  answer_ali_pfunc_macrostate ambd(answer_ali_pfunc_macrostate le,Subsequence b,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pfunc.q1 = scale(1) * check_tuple(le.pfunc.q1, le.firststem, re.firststem, b, re.pfunc);
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(1) * check_tuple(le.covar.q1, le.firststem, re.firststem, b, re.covar);
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate ambd_Pr(answer_ali_pfunc_macrostate le,Subsequence b,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;
    res.pfunc = mk_tuple(le.firststem, scale(1) * check_tuple(le.pfunc.q1, le.firststem, re.firststem, b, re.pfunc));
    res.covar = mk_tuple(le.firststem, scale(1) * check_tuple(le.covar.q1, le.firststem, re.firststem, b, re.covar));
    
    return res;
  }

  answer_ali_pfunc_macrostate nil(Subsequence loc) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = seq_size(loc);
    res.firststem.j = seq_size(loc);
    res.firststem.seq = loc.seq;

    res.pfunc.q1 = 1.0;
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;

    res.covar.q1 = 1.0;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
	  
    return res;
  }

  answer_ali_pfunc_macrostate edl(Subsequence lb,answer_ali_pfunc_macrostate e, Subsequence rloc) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;

    res.pfunc.q1 = scale(1) * e.pfunc.q1 * mk_pf(int((dl_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem)) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
	
	res.covar.q1 = scale(1) * e.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate edr(Subsequence lloc, answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    
	res.pfunc.q1 = scale(1) * e.pfunc.q1 * mk_pf(int((dr_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem)) / float(rows(rb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
	res.covar.q1 = scale(1) * e.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate edlr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;

    //this minimization is necessary since Turner2004 parameters introduced the ext_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
	int help = min(min(ext_mismatch_energy(e.firststem, e.firststem), dl_energy(e.firststem, e.firststem)), dr_energy(e.firststem, e.firststem));
    
	res.firststem = e.firststem;
    res.pfunc.q1 = scale(2) * e.pfunc.q1 * mk_pf(int((help + termau_energy(e.firststem, e.firststem)) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(2) * e.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate drem(Subsequence lloc, answer_ali_pfunc_macrostate e, Subsequence rloc) {
    answer_ali_pfunc_macrostate res;

	res.firststem = e.firststem;
	
    res.pfunc.q1 = e.pfunc.q1 * mk_pf(int(termau_energy(e.firststem, e.firststem)) / float(rows(lloc)));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;

    res.covar.q1 = e.covar.q1;
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;

    return res;
  }

  answer_ali_pfunc_macrostate sr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pfunc.q1 = scale(2) * e.pfunc.q1 * mk_pf(int(sr_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(2) * e.covar.q1 * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;

    res.pfunc.q1 = scale(region.j - region.i + 2) * mk_pf(int(hl_energy(region) / float(rows(region))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(region.j - region.i + 2) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }


  answer_ali_pfunc_macrostate bl(Subsequence lb,Subsequence lregion,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pfunc.q1 = scale(lregion.j - lregion.i + 2) * e.pfunc.q1 * mk_pf(int(bl_energy(lregion,rb) / float(rows(lregion))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
	res.covar.q1 = scale(lregion.j - lregion.i + 2) * e.covar.q1 * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate br(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rregion,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pfunc.q1 = scale(rregion.j - rregion.i + 2) * e.pfunc.q1 * mk_pf(int(br_energy(lb, rregion) / float(rows(rregion))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;

    res.covar.q1 = scale(rregion.j - rregion.i + 2) * e.covar.q1 * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;

    return res;
  }

  answer_ali_pfunc_macrostate il(Subsequence lb,Subsequence lregion,answer_ali_pfunc_macrostate e,Subsequence rregion,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pfunc.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 2) * e.pfunc.q1 * mk_pf(int(il_energy(lregion, rregion) / float(rows(lregion))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 2) * e.covar.q1 * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate ml(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pfunc.q1 = scale(2) * sum_elems(e.pfunc) * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(2) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }
  
  answer_ali_pfunc_macrostate mldr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    res.pfunc.q1 = scale(3) * sum_elems(e.pfunc) * mk_pf(ml_energy() + ul_energy() + int((dri_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem)) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(3) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mladr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	
	int dangleInternal_wc  = 0.0;
	int dangleInternal_wob = 0.0;
	int dangleClosing = 0.0;
	for (int k = 0; k < int(rows(lb)); k=k+1) {
		base_t rightdanglingBase = base_t(column(seq_char(dr, dr.i), k));
		base_t rightmostBaselastStem = base_t(column(seq_char(e.firststem, dr.i-1), k));
		dangleInternal_wc  = dangleInternal_wc  + dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
		dangleInternal_wob = dangleInternal_wob + dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
		dangleClosing = dangleClosing + getEnergyAtRow(res.firststem, k, 2);
	}
	dangleInternal_wc  = int(dangleInternal_wc  / float(rows(lb)));
	dangleInternal_wob = int(dangleInternal_wob / float(rows(lb)));
	dangleClosing      = int(dangleClosing      / float(rows(lb)));
		
    float amdangle = (e.pfunc.q1 + e.pfunc.q3) * mk_pf(min(dangleInternal_wc,  dangleClosing)) +
                     (e.pfunc.q2 + e.pfunc.q4) * mk_pf(min(dangleInternal_wob, dangleClosing));
    
    res.pfunc.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(3) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mldlr(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    //this minimization is necessary since Turner2004 parameters introduced the ml_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
	int help = min(min(ml_mismatch_energy(res.firststem,res.firststem), dli_energy(res.firststem,res.firststem)), dri_energy(res.firststem,res.firststem));
    
    res.pfunc.q1 = scale(4) * sum_elems(e.pfunc) * mk_pf(ml_energy() + ul_energy() + int((help + termau_energy(res.firststem,res.firststem)) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(4) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mladlr(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	
	int dangleInternal_first_wc  = 0;
	int dangleInternal_first_wob = 0;
	int dangleInternal_last_wc  = 0;
	int dangleInternal_last_wob = 0;
	int dangle_l = 0;
	int dangle_r = 0;
	for (int k = 0; k < int(rows(lb)); k=k+1) {
		base_t leftdanglingBase = base_t(column(seq_char(dl, dl.i),k));
		base_t rightdanglingBase = base_t(column(seq_char(dr, dr.i), k));
		base_t leftmostBasefirstStem = base_t(column(seq_char(e.firststem, dl.i+1), k));
		base_t rightmostBaselastStem = base_t(column(seq_char(e.firststem, dr.i-1), k));
		dangleInternal_first_wc  = dangleInternal_first_wc  + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem));
		dangleInternal_first_wob = dangleInternal_first_wob + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem));
		dangle_l = dangle_l + getEnergyAtRow(res.firststem, k, 1);
		dangleInternal_last_wc  = dangleInternal_last_wc  + dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
		dangleInternal_last_wob = dangleInternal_last_wob + dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
		dangle_r = dangle_r + getEnergyAtRow(res.firststem, k, 2);
	}
	dangleInternal_first_wc  = int(dangleInternal_first_wc  / float(rows(lb)));
	dangleInternal_first_wob = int(dangleInternal_first_wob / float(rows(lb)));
	dangleInternal_last_wc   = int(dangleInternal_last_wc   / float(rows(lb)));
	dangleInternal_last_wob  = int(dangleInternal_last_wob  / float(rows(lb)));
	dangle_l                 = int(dangle_l                 / float(rows(lb)));
	dangle_r                 = int(dangle_r                 / float(rows(lb)));
	
	float amdangle = (e.pfunc.q1 * mk_pf(min(dangleInternal_first_wc,  dangle_l) + min(dangleInternal_last_wc,  dangle_r)) +
				      e.pfunc.q2 * mk_pf(min(dangleInternal_first_wc,  dangle_l) + min(dangleInternal_last_wob, dangle_r)) +
				      e.pfunc.q3 * mk_pf(min(dangleInternal_first_wob, dangle_l) + min(dangleInternal_last_wc,  dangle_r)) +
				      e.pfunc.q4 * mk_pf(min(dangleInternal_first_wob, dangle_l) + min(dangleInternal_last_wob, dangle_r)));
	
    res.pfunc.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(4) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mldladr(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	
	int dangle_l = 0;
	int dangle_r = 0;
	int dangleInternal_wc = 0;
	int dangleInternal_wob = 0;
	for (int k = 0; k < int(rows(lb)); k=k+1) {
		base_t rightdanglingBase = base_t(column(seq_char(dr, dr.i), k));
		base_t rightmostBaselastStem = base_t(column(seq_char(e.firststem, dr.i-1), k));
		dangle_l = dangle_l + getEnergyAtRow(res.firststem, k, 1);
		dangle_r = dangle_r + getEnergyAtRow(res.firststem, k, 2);
		dangleInternal_wc  = dangleInternal_wc  + dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
		dangleInternal_wob = dangleInternal_wob + dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase);
	}
	dangle_l           = int(dangle_l           / float(rows(lb)));
	dangle_r           = int(dangle_r           / float(rows(lb)));
	dangleInternal_wc  = int(dangleInternal_wc  / float(rows(lb)));
	dangleInternal_wob = int(dangleInternal_wob / float(rows(lb)));
	double amdangle = ((e.pfunc.q1 * mk_pf(dangle_l) + e.pfunc.q3) * mk_pf(min(dangleInternal_wc, dangle_r)) +
					   (e.pfunc.q2 + e.pfunc.q4) * mk_pf(min(dangleInternal_wob, dangle_r)));
	
    res.pfunc.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(4) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mladldr(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
	int dangleInternal_wc = 0;
	int dangleInternal_wob = 0;
	int dangle_l = 0;
	int dangle_r = 0;
	for (int k = 0; k < int(rows(lb)); k=k+1) {
		base_t leftdanglingBase = base_t(column(seq_char(dl, dl.i), k));
		base_t leftmostBasefirstStem = base_t(column(seq_char(e.firststem, dl.i+1), k));
		dangleInternal_wc  = dangleInternal_wc  + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem));
		dangleInternal_wob = dangleInternal_wob + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem));
		dangle_l = dangle_l + getEnergyAtRow(res.firststem, k, 1);
		dangle_r = dangle_r + getEnergyAtRow(res.firststem, k, 2);
	}
	dangleInternal_wc  = int(dangleInternal_wc  / float(rows(lb)));
	dangleInternal_wob = int(dangleInternal_wob / float(rows(lb)));
	dangle_l           = int(dangle_l           / float(rows(lb)));
	dangle_r           = int(dangle_r           / float(rows(lb)));
	
	float amdangle = ((e.pfunc.q1 + e.pfunc.q2) * mk_pf(min(dangleInternal_wc, dangle_l)) +
				      (e.pfunc.q3 + e.pfunc.q4  * mk_pf(dangle_r)) * mk_pf(min(dangleInternal_wob, dangle_l)));
	
    res.pfunc.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(4) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mldl(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    res.pfunc.q1 = scale(3) * sum_elems(e.pfunc) * mk_pf(ml_energy() + ul_energy() + int((dli_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem)) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(3) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mladl(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
	int dangleInternal_wc = 0.0;
	int dangleInternal_wob = 0.0;
	int dangleClosing = 0.0;
	for (int k = 0; k < int(rows(lb)); k=k+1) {
		base_t leftdanglingBase = base_t(column(seq_char(dl, dl.i), k));
		base_t leftmostBasefirstStem = base_t(column(seq_char(e.firststem, dl.i+1), k));
		dangleInternal_wc  = dangleInternal_wc + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem));
		dangleInternal_wob = dangleInternal_wc + dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem));
		dangleClosing = dangleClosing + getEnergyAtRow(res.firststem, k, 1);
	}
	dangleInternal_wc  = int(dangleInternal_wc  / float(rows(lb)));
	dangleInternal_wob = int(dangleInternal_wob / float(rows(lb)));
	dangleClosing      = int(dangleClosing      / float(rows(lb)));
	
	float amdangle = (e.pfunc.q1 + e.pfunc.q2) * mk_pf(min(dangleInternal_wc,  dangleClosing)) +
                     (e.pfunc.q3 + e.pfunc.q4) * mk_pf(min(dangleInternal_wob, dangleClosing));

    res.pfunc.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))));
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = scale(3) * sum_elems(e.covar) * mk_pf(covscore(lb, lb.i, rb.i));
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate addss(answer_ali_pfunc_macrostate e,Subsequence rregion) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    res.pfunc = mult_tup(scale(rregion.j - rregion.i) * mk_pf(int(ss_energy(rregion) / float(rows(rregion)))), e.pfunc);
    res.covar = mult_tup(scale(rregion.j - rregion.i), e.covar);

    return res;
  }

  answer_ali_pfunc_macrostate ssadd(Subsequence lregion,answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    res.pfunc = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pfunc.q1 * mk_pf(ul_energy() + int(ss_energy(lregion)) / float(rows(lregion))));
    res.covar = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.covar.q1);
    
    return res;
  }

  answer_ali_pfunc_macrostate trafo(answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
	
    res.pfunc.q1 = sum_elems(e.pfunc);
    res.pfunc.q2 = 0.0;
    res.pfunc.q3 = 0.0;
    res.pfunc.q4 = 0.0;
    
    res.covar.q1 = sum_elems(e.covar);
    res.covar.q2 = 0.0;
    res.covar.q3 = 0.0;
    res.covar.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate incl(answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
	
    res.pfunc = mk_tuple(e.firststem, e.pfunc.q1 * mk_pf(ul_energy()));
    res.covar = mk_tuple(e.firststem, e.covar.q1);

    return res;
  }

  answer_ali_pfunc_macrostate combine(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem = le.firststem;
    
    res.pfunc.q1 = (le.pfunc.q1 + le.pfunc.q2) * (re.pfunc.q1 + re.pfunc.q3);
    res.pfunc.q2 = (le.pfunc.q1 + le.pfunc.q2) * (re.pfunc.q2 + re.pfunc.q4);
    res.pfunc.q3 = (le.pfunc.q3 + le.pfunc.q4) * (re.pfunc.q3 + re.pfunc.q1);
    res.pfunc.q4 = (le.pfunc.q4 + le.pfunc.q3) * (re.pfunc.q4 + re.pfunc.q2);
    
	res.covar.q1 = (le.covar.q1 + le.covar.q2) * (re.covar.q1 + re.covar.q3);
    res.covar.q2 = (le.covar.q1 + le.covar.q2) * (re.covar.q2 + re.covar.q4);
    res.covar.q3 = (le.covar.q3 + le.covar.q4) * (re.covar.q3 + re.covar.q1);
    res.covar.q4 = (le.covar.q4 + le.covar.q3) * (re.covar.q4 + re.covar.q2);
    
    return res;
  }

  answer_ali_pfunc_macrostate acomb(answer_ali_pfunc_macrostate le,Subsequence b,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
   
    int  wcDr = 0;
    int wobDr = 0;
    int  wcDl = 0;
    int wobDl = 0;
	for (int k = 0; k < int(rows(b)); k=k+1) {
		base_t baseLeftStem = base_t(column(seq_char(le.firststem, b.i-1), k));
		base_t baseRightStem = base_t(column(seq_char(re.firststem, b.i+1), k));
		base_t baseAmbigious = base_t(column(seq_char(b, b.i), k));
		 wcDr =  wcDr + dr_dangle_dg(  wc_comp(baseLeftStem), baseLeftStem, baseAmbigious);
		wobDr = wobDr + dr_dangle_dg( wob_comp(baseLeftStem), baseLeftStem, baseAmbigious);
		 wcDl =  wcDl + dl_dangle_dg(baseAmbigious, baseRightStem,  wc_comp(baseRightStem));
		wobDl = wobDl + dl_dangle_dg(baseAmbigious, baseRightStem, wob_comp(baseRightStem));
	}
	 wcDr = int( wcDr / float(rows(b)));
    wobDr = int(wobDr / float(rows(b)));
	 wcDl = int( wcDl / float(rows(b)));
	wobDl = int(wobDl / float(rows(b)));
	
    res.pfunc.q1 = le.pfunc.q1 * (re.pfunc.q1 * mk_pf(min( wcDr, wcDl)) + re.pfunc.q3 * mk_pf(min( wcDr,wobDl))) + 
                   le.pfunc.q2 * (re.pfunc.q1 * mk_pf(min(wobDr, wcDl)) + re.pfunc.q3 * mk_pf(min(wobDr,wobDl)));
    res.pfunc.q2 = le.pfunc.q2 * (re.pfunc.q2 * mk_pf(min(wobDr, wcDl)) + re.pfunc.q4 * mk_pf(min(wobDr,wobDl))) + 
                   le.pfunc.q1 * (re.pfunc.q2 * mk_pf(min( wcDr, wcDl)) + re.pfunc.q4 * mk_pf(min( wcDr,wobDl)));
    res.pfunc.q3 = le.pfunc.q3 * (re.pfunc.q3 * mk_pf(min( wcDr,wobDl)) + re.pfunc.q1 * mk_pf(min( wcDr, wcDl))) +
                   le.pfunc.q4 * (re.pfunc.q3 * mk_pf(min(wobDr,wobDl)) + re.pfunc.q1 * mk_pf(min(wobDr, wcDl)));
    res.pfunc.q4 = le.pfunc.q4 * (re.pfunc.q4 * mk_pf(min(wobDr,wobDl)) + re.pfunc.q2 * mk_pf(min(wobDr, wcDl))) +
                   le.pfunc.q3 * (re.pfunc.q4 * mk_pf(min( wcDr,wobDl)) + re.pfunc.q2 * mk_pf(min( wcDr, wcDl)));

	res.pfunc.q1 = res.pfunc.q1 * scale(1);
	res.pfunc.q2 = res.pfunc.q2 * scale(1);
	res.pfunc.q3 = res.pfunc.q3 * scale(1);
	res.pfunc.q4 = res.pfunc.q4 * scale(1);

	res.covar.q1 = le.covar.q1 * (re.covar.q1 + re.covar.q3) + 
				   le.covar.q2 * (re.covar.q1 + re.covar.q3);
	res.covar.q2 = le.covar.q2 * (re.covar.q2 + re.covar.q4) + 
				   le.covar.q1 * (re.covar.q2 + re.covar.q4);
	res.covar.q3 = le.covar.q3 * (re.covar.q3 + re.covar.q1) +
				   le.covar.q4 * (re.covar.q3 + re.covar.q1);
	res.covar.q4 = le.covar.q4 * (re.covar.q4 + re.covar.q2) +
				   le.covar.q3 * (re.covar.q4 + re.covar.q2);

	res.covar.q1 = res.covar.q1 * scale(1);
	res.covar.q2 = res.covar.q2 * scale(1);
	res.covar.q3 = res.covar.q3 * scale(1);
	res.covar.q4 = res.covar.q4 * scale(1);

	res.firststem = le.firststem;
	
    return res;
  }

  choice [answer_ali_pfunc_macrostate] h([answer_ali_pfunc_macrostate] i) {
    return list(sum(i));
    //~ return i;
  }
}


algebra alg_ali_pfunc_macrostate_id extends alg_ali_pfunc_macrostate {
  choice [answer_ali_pfunc_macrostate] h([answer_ali_pfunc_macrostate] l) {
    return l;
  }
}


