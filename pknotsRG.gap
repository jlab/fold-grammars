import rna

import stacklen

input rna

type shape_t = shape
// type base_t = extern // XXX
type Rope = extern

type string_t = string

signature Algebra(alphabet, comp) {
  comp sumend(Subsequence, Subsequence, Subsequence);
  comp sumss(Subsequence);
	comp sadd(Subsequence, comp);
	comp cadd(comp, comp);
	comp nil(void);
	comp is(Subsequence, comp, Subsequence);
	comp edl(Subsequence, comp, Subsequence);
	comp edr(Subsequence, comp, Subsequence);
	comp edlr(Subsequence, comp, Subsequence);
	comp pk(comp);
	comp pknot(Subsequence, comp, Subsequence, comp, Subsequence,
                comp, Subsequence,
                comp, comp, comp, comp);
	comp kndl(Subsequence, comp);
	comp kndr(comp, Subsequence);
	comp kndlr(Subsequence, comp, Subsequence);
	comp sr(Subsequence, comp, Subsequence);
	comp hl(Subsequence, Subsequence, Subsequence, Subsequence, Subsequence);
	comp bl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp br(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp il(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp ml(Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp mldl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp mldr(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp mldlr(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp addss(comp, Subsequence);
	comp mlstem(comp);
	comp pkml(comp);
	comp frd(comp, Subsequence); //frd j
	comp ul(comp);
	comp emptymid(void); //emptymid k l
	comp midbase(void); //midbase k l
	comp middlro(void); //middlro k l
	comp midregion(comp ;  int, int);
	comp middl(Subsequence, comp); //middl k
	comp middr(comp, Subsequence); //middr   l
	comp middlr(Subsequence, comp, Subsequence); //middlr k l
	comp bkd(Subsequence, comp); //bkd i
	comp pss(Subsequence);
	choice [comp] h([comp]);
}

algebra pretty implements Algebra(alphabet = char, comp = string_t) {
	string_t sadd(Subsequence b, string_t x) {
		string_t res;
		append(res, '.');
		append(res, x);
		return res;
	}

	string_t cadd(string_t x, string_t y) {
		string_t res;
		append(res, x);
		append(res, y);
		return res;
	}

	string_t nil(void) {
		string_t res;
		return res;
	}

	string_t is(Subsequence ld, string_t x, Subsequence rd) {
		return x;
	}

	string_t edl(Subsequence ld, string_t x, Subsequence rd) {
		string_t res;
		append(res, '.');
		append(res, x);
		return res;
	}
 
	string_t edr(Subsequence ld, string_t x, Subsequence rd) {
		string_t res;
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t edlr(Subsequence ld, string_t x, Subsequence rd) {
		string_t res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t pk(string_t x) {
		return x;
	}

	string_t pknot(Subsequence a, string_t frt, Subsequence b,
                string_t mid, Subsequence at,
                string_t bck, Subsequence bt,
                string_t e1, string_t e2, string_t e3, string_t e4) {
          string_t res;
          append(res, '(', size(a));
          append(res, '.');
          append(res, frt);
          append(res, '{', size(b));
          append(res, mid);
          append(res, ')', size(at));
          append(res, bck);
          append(res, '.', 2);
          append(res, '}', size(bt));
          return res;
	}

	string_t kndl(Subsequence ld, string_t x) {
		string_t res;
		append(res, '.');
		append(res, x);
		return res;
	}

	string_t kndr(string_t x, Subsequence rd) {
		string_t res;
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t kndlr(Subsequence ld, string_t x, Subsequence rd) {
		string_t res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

  string_t sumend(Subsequence lb, Subsequence r, Subsequence rb)
  {
string_t res;
return res;
}
  string_t sumss(Subsequence r)
  {
string_t res;
return res;
}

	string_t sr(Subsequence lb, string_t x, Subsequence rb) {
		string_t res;
		append(res, '(');
		append(res, x);
		append(res, ')');
		return res;
	}

	string_t hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, '.', size(r));
		append(res, "))", 2);
		return res;
	}

	string_t bl(Subsequence llb, Subsequence lb, Subsequence lr, string_t x, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	string_t br(Subsequence llb, Subsequence lb, string_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, x);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}

	string_t il(Subsequence llb, Subsequence lb, Subsequence lr, string_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, x);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}

	string_t ml(Subsequence llb, Subsequence lb, string_t x, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	string_t mldl(Subsequence llb, Subsequence lb, Subsequence ld, string_t x, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, '.');
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	string_t mldr(Subsequence llb, Subsequence lb, string_t x, Subsequence rd, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, x);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string_t mldlr(Subsequence llb, Subsequence lb, Subsequence ld, string_t x, Subsequence rd, Subsequence rb, Subsequence rrb) {
		string_t res;
		append(res, "((", 2);
		append(res, '.');
		append(res, x);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string_t addss(string_t x, Subsequence r) {
		string_t res;
		append(res, x);
		append(res, '.', size(r));
		return res;
	}

	string_t mlstem(string_t x) {
		return x;
	}

	string_t pkml(string_t x) {
		return x;
	}


	string_t frd(string_t x, Subsequence ld) { //frd j
		string_t res;
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t ul(string_t x) {
		return x;
	}

	string_t emptymid(void) { //emptymid k l
		string_t res;
		return res;
	}

	string_t midbase(void) { //midbase k l
		string_t res;
		append(res, '.'); //if k+1==l
		return res;
	}

	string_t middlro(void) { //middlro k l
		string_t res;
		append(res, "..", 2); //if k+2==l
		return res;
	}

	string_t midregion(string_t x ; int k, int l) {
		return x;
	}

	string_t middl(Subsequence ld, string_t x) { //middl k
		string_t res;
		append(res, '.');
		append(res, x);
		return res;
	}

	string_t middr(string_t x, Subsequence rd) { //middr   l
		string_t res;
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t middlr(Subsequence ld, string_t x, Subsequence rd) { //middlr k l
		string_t res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

	string_t bkd(Subsequence rd, string_t x) { //bkd i
		string_t res;
		append(res, '.');
		append(res, x);
		return res;
	}
 
	string_t pss(Subsequence r) {
		string_t res;
		append(res, '.', size(r));
		return res;
	}

	choice [string_t] h([string_t] i) {
		return i;
	}
}


grammar pknotsRG uses Algebra(axiom = struct) {
    struct       = sadd(BASE,      struct) |
                   cadd(dangle_Pr, struct) |
                   nil (EMPTY) 
                   # h;

    dangle_Pr    = dangle | 
                   dangleknot
                   # h;
    
    dangle       = is   (LOC,  closed, LOC ) |
                   edl  (BASE, closed, LOC ) |
                   edr  (LOC,  closed, BASE) |
                   edlr (BASE, closed, BASE) 
                   # h;
    
    dangleknot   = pk   (      knot        ) |
                   kndl (BASE, knot        ) |
                   kndr (      knot,   BASE) |
                   kndlr(BASE, knot,   BASE) 
                   # h;

    closed       ={stack   | 
                   hairpin |
                   leftB   | 
                   rightB  | 
                   iloop   | 
                   multiloop} with stackpairing 
                   # h;

    stack        = sr   (      BASE,                          closed,                                            BASE      ) # h;
    hairpin      = hl   (BASE, BASE,                          {REGION with minsize(3)},                          BASE, BASE) # h;
    leftB        = bl   (BASE, BASE, REGION with maxsize(30), closed,                                            BASE, BASE) # h;
    rightB       = br   (BASE, BASE,                          closed,                   REGION with maxsize(30), BASE, BASE) # h;
    iloop        = il   (BASE, BASE, REGION with maxsize(30), closed,                   REGION with maxsize(30), BASE, BASE) # h;
    
    multiloop    ={ml   (BASE, BASE,                          ml_comps1,                                         BASE, BASE) |
                   mldl (BASE, BASE, BASE,                    ml_comps1,                                         BASE, BASE) |
                   mldr (BASE, BASE,                          ml_comps1,                BASE,                    BASE, BASE) |
                   mldlr(BASE, BASE, BASE,                    ml_comps1,                BASE,                    BASE, BASE) } with stackpairing
                   # h;

    ml_comps1    = sadd (BASE,             ml_comps1) |
                   cadd (mldangle,         ml_comps)  |
                   addss(pkml(dangleknot), REGION0)
                   # h ;

                     
    ml_comps     = sadd (BASE,             ml_comps) |
                   cadd (mldangle,         ml_comps) |
                   addss(mldangle,         REGION0)
                   # h ;

    mldangle     = mlstem(dangle)     |
                   pkml  (dangleknot)
                   # h;
                     
    knot         = 

      .[
         int i = t_0_i;
         int j = t_0_j;
         int k = t_0_k_0;
         int l = t_0_k_1;
         if (j-i<11)
           continue;
         if (k-i < 3 || j-l < 4)
           continue;
         if (l-k < 4)
           continue;
         int alphamaxlen = stacklen(t_0_seq, i, l);
         if (alphamaxlen < 2)
           continue;
         int alphareallen = min(alphamaxlen, k-i-1);
         if (alphareallen < 2)
           continue;
         int betamaxlen = stacklen(t_0_seq, k, j);
         if (betamaxlen < 2)
           continue;
         int betatemplen = min(betamaxlen, j-l-2);
         int betareallen = min(betatemplen, l-k-alphareallen);
         if (betareallen < 2)
           continue;
       ].
      {
         pknot(REGION, REGION, REGION) .{
           pknot(REGION[i, i+alphareallen],
              front[i+alphareallen+1, k],
              REGION[k, k+betareallen],
              middle[k+betareallen, l-alphareallen] .(k, l).,
              REGION[l-alphareallen, l],
              back[l, j-betareallen-2],
              REGION[j-betareallen, j],
              stacknrg[k, j],
              stacknrg[i, l],
              stacknrg[i+alphareallen, l-alphareallen],
              stacknrg[k+betareallen, j-betareallen] ) 
         }.
      } # h;

/*

>	    pknot       (i,j) = [pk' energy a u b v a' w b' (0,0) | i+11<=j, k <- [i+7 .. j-4],
>				                 (alphanrg, alphalen) <- stacklen (i,k),
>                                alphalen >= 2,
>                                l <- [i+3 .. k-4],
>                                -- don't let a-a' run into b-b' from behind
>                                let h = min alphalen (l-i-1),
>                                h>=2,
>					             (betanrg, betalen) <- stacklen (l,j),
>				                 betalen >=2, 
>                                -- don't let b-b' run into a-a' from behind
>                                let tmph' = min betalen (j-k-2),
>                                -- don't let a-a' and b-b'  collide in the middle
>					             let h' = min tmph'  (k-l-h),
>					             h' >= 2,

>					             a <- region   (i     , i+h  ),
>					             u <- front  j (i+h+1,    l  ),
>					             b <- region   (l    , l+h'  ),
>				                 v <- middle (j-h') (i+h) (l+h', k-h  ),
>					             a'<- region   (k-h  , k     ),
>					             w <- back   i (k    , j-h'-2),
>					             b'<- region   (j-h' , j     ),

>                                -- recalculate the energy of shrinked helices
>                                (acorrectionterm, _) <- stacklen (i+h -1,k-h +1),
>					             (bcorrectionterm, _) <- stacklen (l+h'-1,j-h'+1),
>					             let energy = alphanrg + betanrg - acorrectionterm - bcorrectionterm
>    					    ]	                     

*/
    
                     
    front       = front_Pr               |
                   frd  (front_Pr, BASE)
                   # h;
              
    front_Pr     = ul(emptystrand) |
                   pk_comps
				   # h;
               
    middle(int k, int l)     = emptymid(EMPTY)                     |
                   midbase(EMPTY)                      |
                   middlro(EMPTY)                      |
                   midregion     (      mid    ;  k, l) |
                   middl        (BASE, mid      ) |
                   middr        (      mid, BASE) |
                   middlr      (BASE, mid, BASE) 
                   # h;
    
    mid          = ul(singlestrand) |
                   pk_comps
                   # h;
          
    back        = back_Pr               |
                   bkd  (BASE, back_Pr) 
                   # h;
             
    back_Pr      = ul(emptystrand) |
                   pk_comps
                   # h;
              
    pk_comps     = cadd(singlestrand, pk_comps)    |
                   cadd(mldangle, pk_comps)        |
                   cadd(mldangle, ul(emptystrand)) 
                   # h;
    
    singlestrand = pss(REGION) # h;
    
    emptystrand  = pss(REGION0) # h ;

/*
    stacknrg = sr(BASE, stacknrg, BASE) with stackpairing |
               hairpin # h ;
*/

    stacknrg = sr(BASE, stacknrg, BASE) with basepairing |
               sumend(BASE, REGION with minsize(3), BASE) with basepairing |
               sumss(REGION0) # h ;

/*

>   stacklen = tabulated(
>              (sum    <<<   base +~~ stacklen                        ~~+ base)  `with` basepair  |||
>              (sumend <<<   base +~~ (region `with` (minloopsize 3)) ~~+ base)  `with` basepair  ...hmin)
>           where    sum      lb (c,k) rb   = (c + sr_energy inp (lb, rb), k+1)
>                    sumend   lb _ rb   = (0,1)
>   hmin []  = []
>   hmin xs = [minimum xs]

*/

}


instance pretty = pknotsRG(pretty) ;


