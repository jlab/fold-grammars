//WARNING: this grammar is of experimental state (Stefan Janssen, 25.07.2011)

grammar gra_macrostate_centers uses sig_foldrna(axiom = struct) {
  struct = left_danglePre | trafo(noleft_danglePre) | left_unpairedPre | left_unpairedEnd# h;

  left_unpairedEnd  = sadd(BASE, left_unpairedEnd)  | sadd(BASE, nil(LOC))        # h;
  left_unpairedPre  = sadd(BASE, left_unpairedPre)  | sadd(BASE, left_danglePre)  # h;
  left_unpairedPost = sadd(BASE, left_unpairedPost) | sadd(BASE, left_danglePost) # h;

  left_danglePre  = ambd(edanglel, BASE, noleft_danglePre)  | cadd_Pr(edanglel, noleft_danglePre)               | cadd(edanglelr, {left_danglePre  | left_unpairedPre})   | ambd(edanglelTag, BASE, noleft_danglePost) | cadd_Pr(edanglelTag, {noleft_danglePost | nil(LOC)}) | cadd(edanglelrTag, {left_danglePost | left_unpairedPost}) # h;
  left_danglePost = ambd(edanglel, BASE, noleft_danglePost) | cadd_Pr(edanglel, {noleft_danglePost | nil(LOC)}) | cadd(edanglelr, {left_danglePost | left_unpairedPost}) | nil(LOC) # h;

  noleft_danglePre  = cadd_Pr_Pr(edangler, {left_danglePre | left_unpairedPre})   | cadd_Pr_Pr_Pr(nodangle, noleft_danglePre)               | ambd_Pr(nodangle, BASE, noleft_danglePre)  | cadd_Pr_Pr(edanglerTag, {left_danglePost | left_unpairedPost}) | cadd_Pr_Pr_Pr(nodangleTag, {noleft_danglePost | nil(LOC)}) | ambd_Pr(nodangleTag, BASE, noleft_danglePost) # h;
  noleft_danglePost = cadd_Pr_Pr(edangler, {left_danglePost | left_unpairedPost}) | cadd_Pr_Pr_Pr(nodangle, {noleft_danglePost | nil(LOC)}) | ambd_Pr(nodangle, BASE, noleft_danglePost) # h;

  edanglelTag = edl(BASE, strongTag, LOC) # h;
  edanglel    = edl(BASE, strong, LOC)    # h;

  edanglerTag = edr(LOC, strongTag, BASE) # h;
  edangler    = edr(LOC, strong, BASE)    # h;

  edanglelrTag = edlr(BASE, strongTag, BASE) # h;
  edanglelr    = edlr(BASE, strong, BASE)    # h;

  nodangleTag = drem(LOC, strongTag, LOC) # h;
  nodangle    = drem(LOC, strong, LOC)    # h;

  strongTAG = {sr(BASE, weakTag, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weakTag                     } with allowLonelyBasepairs(true)  # h;
  strong    = {sr(BASE, weak,    BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weak                        } with allowLonelyBasepairs(true)  # h;

  weakTag = stackTag | hairpinTag | multiloopPre | leftBTag | rightBTag | iloopTag # h;
  weak    = stack    | hairpin    | multiloop    | leftB    | rightB    | iloop    # h;

  multiloopPre = {mldl   (BASE, BASE, ml_comps1Pre,        BASE) with basepair | 
                  mladl  (BASE, BASE, ml_comps2Pre,        BASE) with basepair | 
                  mldr   (BASE,       ml_comps3Pre,  BASE, BASE) with basepair | 
                  mladr  (BASE,       ml_comps2Pre,  BASE, BASE) with basepair | 
                  mldlr  (BASE, BASE, ml_comps4Pre,  BASE, BASE) with basepair | 
                  mladlr (BASE, BASE, ml_comps2Pre,  BASE, BASE) with basepair | 
                  mldladr(BASE, BASE, ml_comps1Pre,  BASE, BASE) with basepair | 
                  mladldr(BASE, BASE, ml_comps3Pre,  BASE, BASE) with basepair | 
                  ml     (BASE,       ml_comps2Pre,        BASE) with basepair} with basepair # h;

  multiloop =    {mldl   (BASE, BASE, ml_comps1Post,       BASE) with basepair | 
                  mladl  (BASE, BASE, ml_comps2Post,       BASE) with basepair | 
                  mldr   (BASE,       ml_comps3Post, BASE, BASE) with basepair | 
                  mladr  (BASE,       ml_comps2Post, BASE, BASE) with basepair | 
                  mldlr  (BASE, BASE, ml_comps4Post, BASE, BASE) with basepair | 
                  mladlr (BASE, BASE, ml_comps2Post, BASE, BASE) with basepair | 
                  mldladr(BASE, BASE, ml_comps1Post, BASE, BASE) with basepair | 
                  mladldr(BASE, BASE, ml_comps3Post, BASE, BASE) with basepair | 
                  ml     (BASE,       ml_comps2Post,       BASE) with basepair} with basepair # h;

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

  stack      = sr(BASE,                          weak,                               BASE) with basepair # h;
  stackTag   = sr(BASE,                          weakTag,                            BASE) with basepair # h;

  hairpin    = hl   (BASE,                          REGION with minsize(3),             BASE) with basepair # h;
  hairpinTag = hlTag(BASE,                          REGION with minsize(3),             BASE) with basepair # h;

  leftB      = bl(BASE, REGION,                  strong,                             BASE) with basepair # h;
  leftBTag   = bl(BASE, REGION,                  strongTag,                          BASE) with basepair # h;

  rightB     = br(BASE,                          strong, REGION,                     BASE) with basepair # h;
  rightBTag  = br(BASE,                          strongTag, REGION,                  BASE) with basepair # h;

  iloop      = il(BASE, REGION with maxsize(30), strong, REGION with maxsize(30),    BASE) with basepair # h;
  iloopTag   = il(BASE, REGION with maxsize(30), strongTag, REGION with maxsize(30), BASE) with basepair # h;
}