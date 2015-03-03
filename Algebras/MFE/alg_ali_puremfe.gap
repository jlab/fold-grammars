algebra alg_ali_puremfe implements sig_foldrna(alphabet = M_Char, answer = int) {
  int sadd(Subsequence lb, int x) {
	int sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
    return x + (sbase_sum / float(rows(lb)));
  }
  int cadd(int x, int y) {
	  return x+y;
  }
  int edl(Subsequence ldangle, int x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	return x + ((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));
  }
  int edr(Subsequence lb, int x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	return x + ((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));
  }
  int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	return x + ((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));
  }
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x + (termau_energy(lb, rb) / float(rows(lb)));
  }
  int dall(Subsequence lb, int x, Subsequence rb) {
	return x + ((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
  }
  int sr(Subsequence lb, int x, Subsequence rb) {
    return x + (sr_energy(lb, rb) / float(rows(lb)));
  }
  int hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return (hl_energy(r) / float(rows(r)));
  }
  int bl(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
	return x + (bl_energy(lr, rb) / float(rows(lb)));
  }
  int br(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
    return x + (br_energy(lb, rr) / float(rows(lb)));
  }
  int il(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
    return x + (il_energy(lr, rr) / float(rows(lr)));
  }
  int mldl(Subsequence lb, Subsequence dl, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb)));
  }
  int mldr(Subsequence lb, int x, Subsequence dr, Subsequence rb) {
    return x + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb)));
  }
  int mldlr(Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb) {
    return x + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb)));
  }
  int mlall(Subsequence lb, int x, Subsequence rb) {
	return x + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
  }
  int incl(int x) {
    return x + ul_energy();
  }
  int addss(int x, Subsequence r) {
    return x + (ss_energy(r) / float(rows(r)));
  }
  int nil(Subsequence n) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
	
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  int acomb(int le,Subsequence b,int re) {return le+re;}
  int combine(int le,int re) {return le+re;}
  int trafo(int e) {return e;}
  int ssadd(Subsequence lb,int e) {return e;}
  int mladl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e;}
  int mladldr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mldladr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e;}
  int ambd_Pr(int le,Subsequence b,int re) {return le+re;}
  int ambd(int le,Subsequence b,int re) {return le+re;}
  int cadd_Pr_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr(int le,int re) {return le+re;}
}
