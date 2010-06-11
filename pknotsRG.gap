import rna

import stacklen

input rna

type shape_t = shape
type base_t = extern
type Rope = extern

signature Algebra(alphabet, comp) {
	comp sadd(Subsequence, comp);
	comp cadd(comp, comp);
	comp nil(void);
	comp is(Subsequence, comp, Subsequence);
	comp edl(Subsequence, comp, Subsequence);
	comp edr(Subsequence, comp, Subsequence);
	comp edlr(Subsequence, comp, Subsequence);
	comp pk(comp);
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
	comp knot(Subsequence, Subsequence, Subsequence);
	comp frd(comp, Subsequence); //frd j
	comp ul(comp);
	comp emptymid(void); //emptymid k l
	comp midbase(void); //midbase k l
	comp middlro(void); //middlro k l
	comp midregion(comp);
	comp middl(Subsequence, comp); //middl k
	comp middr(comp, Subsequence); //middr   l
	comp middlr(Subsequence, comp, Subsequence); //middlr k l
	comp bkd(Subsequence, comp); //bkd i
	comp pss(Subsequence);
	comp sum(Subsequence, comp, Subsequence);
	comp sumend(Subsequence, Subsequence, Subsequence);
	choice [comp] h([comp]);
}

algebra pretty implements jensAlgebra(alphabet = char, comp = Rope) {
	Rope sadd(Subsequence b, Rope x) {
		Rope res;
		append(res, '.');
		append(res, e);
		return res;
	}

	Rope cadd(Rope x, Rope y) {
		Rope res;
		append(res, x);
		append(res, y);
		return res;
	}

	Rope nil(void) {
		Rope res;
		return res;
	}

	Rope is(Subsequence ld, Rope x, Subsequence rd) {
		return x;
	}

	Rope edl(Subsequence ld, Rope x, Subsequence rd) {
		Rope res;
		append(res, '.');
		append(res, x);
		return res;
	}
 
	Rope edr(Subsequence ld, Rope x, Subsequence rd) {
		Rope res;
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope edlr(Subsequence ld, Rope x, Subsequence rd) {
		Rope res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope pk(Rope x) {
		return x;
	}

	Rope kndl(Subsequence ld, Rope x) {
		Rope res;
		append(res, '.');
		append(res, x);
		return res;
	}

	Rope kndr(Rope x, Subsequence rd) {
		Rope res;
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope kndlr(Subsequence ld, Rope x, Subsequence rd) {
		Rope res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope sr(Subsequence lb, Rope x, Subsequence rb) {
		Rope res;
		append(res, '(');
		append(res, x);
		append(res, ')');
		return res;
	}

	Rope hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(r));
		append(res, "))", 2);
		return res;
	}

	Rope bl(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	Rope br(Subsequence llb, Subsequence lb, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, x);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}

	Rope il(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, x);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}

	Rope ml(Subsequence llb, Subsequence lb, Rope x, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	Rope mldl(Subsequence llb, Subsequence lb, Subsequence ld, Rope x, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.');
		append(res, x);
		append(res, "))", 2);
		return res;
	}

	Rope mldr(Subsequence llb, Subsequence lb, Rope x, Subsequence rd, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, x);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	Rope mldlr(Subsequence llb, Subsequence lb, Subsequence ld, Rope x, Subsequence rd, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.');
		append(res, x);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	Rope addss(Rope x, Subsequence r) {
		Rope res;
		append(res, x);
		append(res, '.', size(r));
		return res;
	}

	Rope mlstem(Rope x) {
		return x;
	}

	Rope pkml(Rope x) {
		return x;
	}

	Rope knot(Subsequence a, Subsequence b, Subsequence c) {
		Rope res;
		return res;
	}

	Rope frd(Rope x, Subsequence ld) { //frd j
		Rope res;
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope ul(Rope x) {
		return x;
	}

	Rope emptymid(void) { //emptymid k l
		Rope res;
		return res;
	}

	Rope midbase(void) { //midbase k l
		Rope res;
		append(res, '.'); //if k+1==l
		return res;
	}

	Rope middlro(void) { //middlro k l
		Rope res;
		append(res, "..", 2); //if k+2==l
		return res;
	}

	Rope midregion(Rope x) {
		return x;
	}

	Rope middl(Subsequence ld, Rope x) { //middl k
		Rope res;
		append(res, '.');
		append(res, x);
		return res;
	}

	Rope middr(Rope x, Subsequence rd) { //middr   l
		Rope res;
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope middlr(Subsequence ld, Rope x, Subsequence rd) { //middlr k l
		Rope res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
	}

	Rope bkd(Subsequence rd, Rope x) { //bkd i
		Rope res;
		append(res, '.');
		append(res, x);
		return res;
	}
 
	Rope pss(Subsequence r) {
		Rope res;
		append(res, '.', size(r));
		return res;
	}

	Rope sum(Subsequence lb, Rope x, Subsequence rb) {
		Rope res;
		append(res, '(');
		append(res, r);
		append(res, ')');
		return res;
	}

	Rope sumend(Subsequence lb, Subsequence r, Subsequence rb) {
		Rope res;
		append(res, '(');
		append(res, '.', size(r));
		append(res, ')');
		return res;
	}

	choice [Rope] h([Rope] i) {
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
         unsigned i = t_0_i;
         unsigned j = t_0_j;
         unsigned k = t_0_k_0;
         unsigned l = t_0_k_1;
         if (i+11>j)
           continue;
         if (k-i < 3 || j-l < 4)
           continue;
         if (l-k < 4)
           continue;
         unsigned alphamaxlen = stacklen(i, l);
         if (alphamaxlen < 2)
           continue;
         unsigned alphareallen = min(alphamaxlen, k-i-1);
         if (alphareallen < 2)
           continue;
         unsigned betamaxlen = stacklen(k, j);
         if (betamaxlen < 2)
           continue;
         unsigned betatemplen = min(betamaxlen, j-l-2);
         unsigned betareallen = min(betatemplen, l-k-alphareallen);
         if (betareallen < 2)
           continue;
       ].
      {
         pk(REGION, REGION, REGION) .{
           pk(REGION[i, i+alphareallen],
              front[i+alphareallen+1, k],
              REGION[k, k+betareallen],
              middle[k+betareallen, l-alphareallen],
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
    
                     
    front j      = front_Pr               |
                   frd j (front_Pr, BASE)
                   # h;
              
    front_Pr     = ul(emptystrand) |
                   pk_comps
				   # h;
               
    middle k l   = emptymid  k l                   |
                   midbase   k l                   |
                   middlro   k l                   |
                   midregion     (      mid      ) |
                   middl     k   (BASE, mid      ) |
                   middr       l (      mid, BASE) |
                   middlr    k l (BASE, mid, BASE) 
                   # h;
    
    mid          = ul(singlestrand) |
                   pk_comps
                   # h;
          
    back i       = back_Pr               |
                   bkd i (BASE, back_Pr) 
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

    stacknrg = sr(BASE, stacknrg, BASE) with stackpairing |
               sr(BASE, REGION with minsize(3), BASE) # h ;

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
