//WARNING: this grammar is of experimental state (Stefan Janssen, 25.07.2011)

//This is the grammar, developed by Bjoern Voss, for the probablistic shape analysis of RNAshapes 2006 release. It is also known as "canonicals_nonamb" in the Haskell version of RNAshapes, or "adpf_nonamb"

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  Non-Terminals:
//    left_dg = left_dangle
//    noleft_dg = noleft_dangle
//    edgl = edanglel
//    edgr = edangler
//    edglr = edanglelr
//    nodg = nodangle
//    mc1 = ml_comps1
//    mc2 = ml_comps2
//    mc3 = ml_comps3
//    mc4 = ml_comps4
//    mcadd1 = no_dl_no_ss_end
//    mcadd2 = dl_or_ss_left_no_ss_end
//    mcadd3 = dl_or_ss_left_ss_end
//    mcadd4 = no_dl_ss_end
grammar gra_macrostate_centers uses sig_foldrna(axiom = struct) {
	struct = left_danglePre | trafo(noleft_danglePre) | left_unpairedPre | left_unpairedEnd# h;

	left_unpairedEnd  = sadd(BASE, left_unpairedEnd)  | sadd(BASE, nil(LOC))        # h;
	left_unpairedPre  = sadd(BASE, left_unpairedPre)  | sadd(BASE, left_danglePre)  # h;
	left_unpairedPost = sadd(BASE, left_unpairedPost) | sadd(BASE, left_danglePost) # h;

	left_danglePre  = ambd(edanglel, BASE, noleft_danglePre)  | cadd_Pr(edanglel, noleft_danglePre)               | cadd(edanglelr, {left_danglePre  | left_unpairedPre})   | ambd(edanglelTag, BASE, noleft_danglePost) | cadd_Pr(edanglelTag, {noleft_danglePost | nil(LOC)}) | cadd(edanglelrTag, {left_danglePost | left_unpairedPost}) # h;
	left_danglePost = ambd(edanglel, BASE, noleft_danglePost) | cadd_Pr(edanglel, {noleft_danglePost | nil(LOC)}) | cadd(edanglelr, {left_danglePost | left_unpairedPost}) | nil(LOC) # h;

	noleft_danglePre  = cadd_Pr_Pr(edangler, {left_danglePre | left_unpairedPre})   | cadd_Pr_Pr_Pr(nodangle, noleft_danglePre)               | ambd_Pr(nodangle, BASE, noleft_danglePre)  | cadd_Pr_Pr(edanglerTag, {left_danglePost | left_unpairedPost}) | cadd_Pr_Pr_Pr(nodangleTag, {noleft_danglePost | nil(LOC)}) | ambd_Pr(nodangleTag, BASE, noleft_danglePost) # h;
	noleft_danglePost = cadd_Pr_Pr(edangler, {left_danglePost | left_unpairedPost}) | cadd_Pr_Pr_Pr(nodangle, {noleft_danglePost | nil(LOC)}) | ambd_Pr(nodangle, BASE, noleft_danglePost) # h;

	edanglelTag = edl(BASE, closedTag, LOC) # h;
	edanglel    = edl(BASE, closed, LOC)    # h;

	edanglerTag = edr(LOC, closedTag, BASE) # h;
	edangler    = edr(LOC, closed, BASE)    # h;

	edanglelrTag = edlr(BASE, closedTag, BASE) # h;
	edanglelr    = edlr(BASE, closed, BASE)    # h;

	nodangleTag = drem(LOC, closedTag, LOC) # h;
	nodangle    = drem(LOC, closed, LOC)    # h;

	closedTag = stackTag | hairpinTag | multiloopPre | leftBTag | rightBTag | iloopTag # h;
	closed    = stack    | hairpin    | multiloop    | leftB    | rightB    | iloop    # h;

	multiloopPre = {mldl(BASE, BASE, BASE, ml_comps1Pre, BASE, BASE) | 
		mladl(BASE, BASE, BASE, ml_comps2Pre, BASE, BASE) | 
		mldr(BASE, BASE, ml_comps3Pre, BASE, BASE, BASE) | 
		mladr(BASE, BASE, ml_comps2Pre, BASE, BASE, BASE) | 
		mldlr(BASE, BASE, BASE, ml_comps4Pre, BASE, BASE, BASE) | 
		mladlr(BASE, BASE, BASE, ml_comps2Pre, BASE, BASE, BASE) | 
		mldladr(BASE, BASE, BASE, ml_comps1Pre, BASE, BASE, BASE) | 
		mladldr(BASE, BASE, BASE, ml_comps3Pre, BASE, BASE, BASE)  | 
		ml(BASE, BASE, ml_comps2Pre, BASE, BASE)} with stackpairing # h;

	multiloop = {mldl(BASE, BASE, BASE, ml_comps1Post, BASE, BASE) | 
		mladl(BASE, BASE, BASE, ml_comps2Post, BASE, BASE) | 
		mldr(BASE, BASE, ml_comps3Post, BASE, BASE, BASE) | 
		mladr(BASE, BASE, ml_comps2Post, BASE, BASE, BASE) | 
		mldlr(BASE, BASE, BASE, ml_comps4Post, BASE, BASE, BASE) | 
		mladlr(BASE, BASE, BASE, ml_comps2Post, BASE, BASE, BASE) | 
		mldladr(BASE, BASE, BASE, ml_comps1Post, BASE, BASE, BASE) | 
		mladldr(BASE, BASE, BASE, ml_comps3Post, BASE, BASE, BASE)  | 
		ml(BASE, BASE, ml_comps2Post, BASE, BASE)} with stackpairing # h;

	ml_comps1Pre  = combine(block_dl, no_dl_no_ss_endPre)  | combine(block_dlr, dl_or_ss_left_no_ss_endPre)  | acomb(block_dl, BASE, no_dl_no_ss_endPre)  | combine(block_dlTag, no_dl_no_ss_endPost) | combine(block_dlrTag, dl_or_ss_left_no_ss_endPost) | acomb(block_dlTag, BASE, no_dl_no_ss_endPost) # h;
	ml_comps1Post = combine(block_dl, no_dl_no_ss_endPost) | combine(block_dlr, dl_or_ss_left_no_ss_endPost) | acomb(block_dl, BASE, no_dl_no_ss_endPost) # h;

	ml_comps2Pre  = combine(incl(nodangle), no_dl_no_ss_endPre)  | combine(incl(edangler), dl_or_ss_left_no_ss_endPre)  | acomb(incl(nodangle), BASE, no_dl_no_ss_endPre)  | combine(incl(nodangleTag), no_dl_no_ss_endPost) | combine(incl(edanglerTag), dl_or_ss_left_no_ss_endPost) | acomb(incl(nodangleTag), BASE, no_dl_no_ss_endPost) # h;
	ml_comps2Post = combine(incl(nodangle), no_dl_no_ss_endPost) | combine(incl(edangler), dl_or_ss_left_no_ss_endPost) | acomb(incl(nodangle), BASE, no_dl_no_ss_endPost) # h;

	ml_comps3Pre  = combine(incl(edangler), dl_or_ss_left_ss_endPre)  | combine(incl(nodangle), no_dl_ss_endPre)  | acomb(incl(nodangle), BASE, no_dl_ss_endPre)  | combine(incl(edanglerTag), dl_or_ss_left_ss_endPost) | combine(incl(nodangleTag), no_dl_ss_endPost) | acomb(incl(nodangleTag), BASE, no_dl_ss_endPost) # h;
	ml_comps3Post = combine(incl(edangler), dl_or_ss_left_ss_endPost) | combine(incl(nodangle), no_dl_ss_endPost) | acomb(incl(nodangle), BASE, no_dl_ss_endPost) # h;

	ml_comps4Pre  = combine(block_dl, no_dl_ss_endPre)  | combine(block_dlr, dl_or_ss_left_ss_endPre)  | acomb(block_dl, BASE, no_dl_ss_endPre)  | combine(block_dlTag, no_dl_ss_endPost) | combine(block_dlrTag, dl_or_ss_left_ss_endPost) | acomb(block_dlTag, BASE, no_dl_ss_endPost) # h;
	ml_comps4Post = combine(block_dl, no_dl_ss_endPost) | combine(block_dlr, dl_or_ss_left_ss_endPost) | acomb(block_dl, BASE, no_dl_ss_endPost) # h;

	block_dlTag = ssadd(REGION, edanglelTag) | incl(edanglelTag) # h;
	block_dl    = ssadd(REGION, edanglel)    | incl(edanglel)    # h;

	block_dlrTag = ssadd(REGION, edanglelrTag) | incl(edanglelrTag) # h;
	block_dlr    = ssadd(REGION, edanglelr)    | incl(edanglelr)    # h;

	no_dl_no_ss_endPre  = ml_comps2Pre  | incl(nodangleTag) # h;
	no_dl_no_ss_endPost = ml_comps2Post | incl(nodangle) # h;

	dl_or_ss_left_no_ss_endPre  = ml_comps1Pre  | block_dlTag # h;
	dl_or_ss_left_no_ss_endPost = ml_comps1Post | block_dl # h;

	no_dl_ss_endPre  = ml_comps3Pre  | incl(edanglerTag) | addss(incl(edanglerTag), REGION) # h;
	no_dl_ss_endPost = ml_comps3Post | incl(edangler)    | addss(incl(edangler),    REGION) # h;

	dl_or_ss_left_ss_endPre  = ml_comps4Pre  | block_dlrTag | addss(block_dlrTag, REGION) # h;
	dl_or_ss_left_ss_endPost = ml_comps4Post | block_dlr    | addss(block_dlr,    REGION) # h;

	stackTag = sr(BASE, closedTag, BASE) with basepairing # h;
	stack    = sr(BASE, closed,    BASE) with basepairing # h;

	hairpinTag = hlTag(BASE, BASE, REGION with minsize(3), BASE, BASE) with stackpairing # h;
	hairpin    = hl(BASE, BASE, REGION with minsize(3), BASE, BASE) with stackpairing # h;

	leftBTag = bl(BASE, BASE, REGION, closedTag, BASE, BASE) with stackpairing # h;
	leftB    = bl(BASE, BASE, REGION, closed,    BASE, BASE) with stackpairing # h;

	rightBTag = br(BASE, BASE, closedTag, REGION, BASE, BASE) with stackpairing # h;
	rightB    = br(BASE, BASE, closed,    REGION, BASE, BASE) with stackpairing # h;

	iloopTag = il(BASE, BASE, REGION with maxsize(30), closedTag, REGION with maxsize(30), BASE, BASE) with stackpairing # h;
	iloop    = il(BASE, BASE, REGION with maxsize(30), closed,    REGION with maxsize(30), BASE, BASE) with stackpairing # h;

}