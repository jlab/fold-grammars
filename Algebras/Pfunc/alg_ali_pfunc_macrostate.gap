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
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(sbase_sum / float(rows(lb)));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
	
    return res;
  }

  answer_ali_pfunc_macrostate cadd(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pf.q1 = le.pf.q1 * re.pf.q1;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
        
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pf.q1 = le.pf.q1 * sum_elems(re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    	
	res.firststem = le.firststem;
    res.pf = mk_tuple(le.firststem, le.pf.q1 * re.pf.q1);
    
    return res;
  }

  answer_ali_pfunc_macrostate cadd_Pr_Pr_Pr(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;
    res.pf = mk_tuple(le.firststem, le.pf.q1 * sum_elems(re.pf));
    
    return res;
  }

  answer_ali_pfunc_macrostate ambd(answer_ali_pfunc_macrostate le,Subsequence b,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;

    res.pf.q1 = scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate ambd_Pr(answer_ali_pfunc_macrostate le,Subsequence b,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = le.firststem;
    res.pf = mk_tuple(le.firststem, scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf));
    
    return res;
  }

  answer_ali_pfunc_macrostate nil(Subsequence loc) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = seq_size(loc);
    res.firststem.j = seq_size(loc);
    res.firststem.seq = loc.seq;

    res.pf.q1 = 1.0;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  answer_ali_pfunc_macrostate edl(Subsequence lb,answer_ali_pfunc_macrostate e, Subsequence rloc) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;

    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(int((dl_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem)) / float(rows(lb))));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
	
    return res;
  }

  answer_ali_pfunc_macrostate edr(Subsequence lloc, answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    
	res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(int((dr_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem)) / float(rows(rb))));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate edlr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;

    //this minimization is necessary since Turner2004 parameters introduced the ext_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
	int help = min(min(ext_mismatch_energy(e.firststem, e.firststem), dl_energy(e.firststem, e.firststem)), dr_energy(e.firststem, e.firststem));
    
	res.firststem = e.firststem;
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(int((help + termau_energy(e.firststem, e.firststem)) / float(rows(lb))));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate drem(Subsequence lloc, answer_ali_pfunc_macrostate e, Subsequence rloc) {
    answer_ali_pfunc_macrostate res;

	res.firststem = e.firststem;
	
    res.pf.q1 = e.pf.q1 * mk_pf(int(termau_energy(e.firststem, e.firststem)) / float(rows(lloc)));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  answer_ali_pfunc_macrostate sr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(int(sr_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;

    res.pf.q1 = scale(region.j - region.i + 2) * mk_pf(int(hl_energy(region) / float(rows(region))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }


  answer_ali_pfunc_macrostate bl(Subsequence lb,Subsequence lregion,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pf.q1 = scale(lregion.j - lregion.i + 2) * e.pf.q1 * mk_pf(int(bl_energy(lregion,rb) / float(rows(lregion))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate br(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rregion,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pf.q1 = scale(rregion.j - rregion.i + 2) * e.pf.q1 * mk_pf(int(br_energy(lb, rregion) / float(rows(rregion))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  answer_ali_pfunc_macrostate il(Subsequence lb,Subsequence lregion,answer_ali_pfunc_macrostate e,Subsequence rregion,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pf.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 2) * e.pf.q1 * mk_pf(int(il_energy(lregion, rregion) / float(rows(lregion))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate ml(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
	res.firststem.seq = lb.seq;
    
    res.pf.q1 = scale(2) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }
  
  answer_ali_pfunc_macrostate mldr(Subsequence lb,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + int((dri_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
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
		
    float amdangle = (e.pf.q1 + e.pf.q3) * mk_pf(min(dangleInternal_wc,  dangleClosing)) +
                     (e.pf.q2 + e.pf.q4) * mk_pf(min(dangleInternal_wob, dangleClosing));
    
    res.pf.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mldlr(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence dr,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    //this minimization is necessary since Turner2004 parameters introduced the ml_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
	int help = min(min(ml_mismatch_energy(res.firststem,res.firststem), dli_energy(res.firststem,res.firststem)), dri_energy(res.firststem,res.firststem));
    
    res.pf.q1 = scale(4) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + int((help + termau_energy(res.firststem,res.firststem)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
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
	
	float amdangle = (e.pf.q1 * mk_pf(min(dangleInternal_first_wc,  dangle_l) + min(dangleInternal_last_wc,  dangle_r)) +
				      e.pf.q2 * mk_pf(min(dangleInternal_first_wc,  dangle_l) + min(dangleInternal_last_wob, dangle_r)) +
				      e.pf.q3 * mk_pf(min(dangleInternal_first_wob, dangle_l) + min(dangleInternal_last_wc,  dangle_r)) +
				      e.pf.q4 * mk_pf(min(dangleInternal_first_wob, dangle_l) + min(dangleInternal_last_wob, dangle_r)));
	
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
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
	double amdangle = ((e.pf.q1 * mk_pf(dangle_l) + e.pf.q3) * mk_pf(min(dangleInternal_wc, dangle_r)) +
					   (e.pf.q2 + e.pf.q4) * mk_pf(min(dangleInternal_wob, dangle_r)));
	
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
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
	
	float amdangle = ((e.pf.q1 + e.pf.q2) * mk_pf(min(dangleInternal_wc, dangle_l)) +
				      (e.pf.q3 + e.pf.q4  * mk_pf(dangle_r)) * mk_pf(min(dangleInternal_wob, dangle_l)));
	
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate mldl(Subsequence lb,Subsequence dl,answer_ali_pfunc_macrostate e,Subsequence rb) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    res.firststem.seq = lb.seq;
	  
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + int((dli_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
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
	
	float amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dangleInternal_wc,  dangleClosing)) +
                     (e.pf.q3 + e.pf.q4) * mk_pf(min(dangleInternal_wob, dangleClosing));

    res.pf.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + int(termau_energy(res.firststem,res.firststem) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate addss(answer_ali_pfunc_macrostate e,Subsequence rregion) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    res.pf = mult_tup(scale(rregion.j - rregion.i) * mk_pf(int(ss_energy(rregion) / float(rows(rregion)))), e.pf);

    return res;
  }

  answer_ali_pfunc_macrostate ssadd(Subsequence lregion,answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
    res.pf = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(ul_energy() + int(ss_energy(lregion)) / float(rows(lregion))));
    
    return res;
  }

  answer_ali_pfunc_macrostate trafo(answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
	
    res.pf.q1 = sum_elems(e.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_ali_pfunc_macrostate incl(answer_ali_pfunc_macrostate e) {
    answer_ali_pfunc_macrostate res;
    
	res.firststem = e.firststem;
	
    res.pf = mk_tuple(e.firststem, e.pf.q1 * mk_pf(ul_energy()));

    return res;
  }

  answer_ali_pfunc_macrostate combine(answer_ali_pfunc_macrostate le,answer_ali_pfunc_macrostate re) {
    answer_ali_pfunc_macrostate res;
    
    res.firststem = le.firststem;
    
    res.pf.q1 = (le.pf.q1 + le.pf.q2) * (re.pf.q1 + re.pf.q3);
    res.pf.q2 = (le.pf.q1 + le.pf.q2) * (re.pf.q2 + re.pf.q4);
    res.pf.q3 = (le.pf.q3 + le.pf.q4) * (re.pf.q3 + re.pf.q1);
    res.pf.q4 = (le.pf.q4 + le.pf.q3) * (re.pf.q4 + re.pf.q2);
    
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
	
    res.pf.q1 = le.pf.q1 * (re.pf.q1 * mk_pf(min( wcDr, wcDl)) + re.pf.q3 * mk_pf(min( wcDr,wobDl))) + 
                   le.pf.q2 * (re.pf.q1 * mk_pf(min(wobDr, wcDl)) + re.pf.q3 * mk_pf(min(wobDr,wobDl)));
    res.pf.q2 = le.pf.q2 * (re.pf.q2 * mk_pf(min(wobDr, wcDl)) + re.pf.q4 * mk_pf(min(wobDr,wobDl))) + 
                   le.pf.q1 * (re.pf.q2 * mk_pf(min( wcDr, wcDl)) + re.pf.q4 * mk_pf(min( wcDr,wobDl)));
    res.pf.q3 = le.pf.q3 * (re.pf.q3 * mk_pf(min( wcDr,wobDl)) + re.pf.q1 * mk_pf(min( wcDr, wcDl))) +
                   le.pf.q4 * (re.pf.q3 * mk_pf(min(wobDr,wobDl)) + re.pf.q1 * mk_pf(min(wobDr, wcDl)));
    res.pf.q4 = le.pf.q4 * (re.pf.q4 * mk_pf(min(wobDr,wobDl)) + re.pf.q2 * mk_pf(min(wobDr, wcDl))) +
                   le.pf.q3 * (re.pf.q4 * mk_pf(min( wcDr,wobDl)) + re.pf.q2 * mk_pf(min( wcDr, wcDl)));

	res.pf.q1 = res.pf.q1 * scale(1);
	res.pf.q2 = res.pf.q2 * scale(1);
	res.pf.q3 = res.pf.q3 * scale(1);
	res.pf.q4 = res.pf.q4 * scale(1);

	res.firststem = le.firststem;
	
    return res;
  }

  choice [answer_ali_pfunc_macrostate] h([answer_ali_pfunc_macrostate] i) {
    return list(sum(i));
    //~ return i;
  }
}


algebra alg_ali_pfunc_id_macrostate extends alg_ali_pfunc_macrostate {
  choice [answer_ali_pfunc_macrostate] h([answer_ali_pfunc_macrostate] l) {
    return l;
  }
}


