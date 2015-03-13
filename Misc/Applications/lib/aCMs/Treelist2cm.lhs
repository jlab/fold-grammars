> module Treelist2cm where

> import Data.List
> import Datatypes
> import Probs
> import Data.Map (Map)
> import qualified Data.Map as Map

> get_train_code :: [TreeList Int] -> String
> get_train_code consensi = "import isntimes\n\n" ++
>                           "type Rope = extern\n" ++
>                           signature ++ "\n" ++
>                           "algebra alg_count auto count;\n" ++
>                           "algebra alg_enum auto enum;\n\n" ++
>                           gen_train_alg_parse ++
>                           gen_train_alg_fake ++
>                           grammar ++ 
>                           "instance count = gra_train(alg_count);\n" ++
>                           "instance train = gra_train(alg_enum);\n"
>         where signature = gen_train_sig
>               grammar   = "grammar gra_train uses sig_cm(axiom = start) {\n" ++
>                           "  start = " ++ (foldr1 (\x -> \y -> x ++ " | " ++ y) [ jump | x <- indexedConsensi, 
>                                                                                           let subNT = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                           let jump = if ((fst (getTreeListIndex x)) == 1) then subNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-1)) ++ "), " ++ subNT ++ ")"
>                                                                                 ]) ++ " # h;\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_train_grammar indexedConsensi)) ++ "\n" ++
>                           "}\n"
>               indexedConsensi = indexTreeList consensi 1

> get_build_code :: [TreeList Int] -> String
> get_build_code consensi = "import probabilities\n" ++
>                           "import isntimes\n\n" ++
>                           "type Rope = extern\n" ++
>                           "type ali = (Rope model, Rope seq, Rope cons)\n\n" ++
>                           signature ++
>                           "algebra alg_count auto count;\n\n" ++
>                           "algebra alg_enum auto enum;\n\n" ++
>                           gen_build_alg_cyk ++ "\n" ++
>                           gen_build_alg_align ++ "\n" ++
>                           grammar ++ "\n" ++
>                           "instance count = gra_build(alg_count);\n" ++
>                           "instance train = gra_build(alg_enum);\n" ++
>                           "instance cyk = gra_build(alg_cyk);\n" ++
>                           "instance cykali = gra_build(alg_cyk * alg_align);\n"
>         where signature = gen_build_sig
>               grammar   = "grammar gra_build uses sig_cm(axiom = start) {\n" ++
>                           "  start = " ++ (foldr1 (\x -> \y -> x ++ " | " ++ y) [ "a_" ++ (show (snd (getTreeListIndex x))) | x <- indexedConsensi ]) ++ " # h;\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_build_grammar indexedConsensi)) ++ "\n" ++
>                           "}\n"
>               indexedConsensi = indexTreeList consensi 1


> -- algebra functions are shared between all alternative consensus, thus they might appear several times in the fusion tree but they have to be implemented just once and without knowledge about their alternative offspring. collectAlgfkts returns the uniq list of algebra function names to be implemented.
> collectAlgfkts :: [TreeList (Int,Int)] -> [(String, Int)]
> collectAlgfkts [El (pos,i)]       = nub (("INS", pos):[
>                                                         ("NIL", pos)
>                                                       ])
> collectAlgfkts [Pl (pos,i) xs ys] = nub (("INS", pos):[
>                                                         ("PK", pos),
>                                                         ("Lr", pos),
>                                                         ("lR", pos),
>                                                         ("bg", pos)
>                                                       ] ++ (collectAlgfkts xs) ++ (collectAlgfkts ys))
> collectAlgfkts [Ol (pos,i) xs]    = nub (("INS", pos):[
>                                                         ("MAT", pos),
>                                                         ("DEL", pos)
>                                                       ] ++ (collectAlgfkts xs))
> collectAlgfkts (a:as)             = nub ((collectAlgfkts [a]) ++ (collectAlgfkts as))
> collectAlgfkts _                  = error ("collectAlgfkts: this case should never ever been reached")


============== SIGNATURES ==============

> -- generates the constant Bellman's GAP code signature for training with given structure-consensi
> gen_train_sig :: String
> gen_train_sig = "signature sig_cm(alphabet, answer) {\n"
>              ++ "  answer INS(alphabet, alphabet, answer; int);\n"
>              ++ "  answer NIL(void; int);\n"
>              ++ "  answer MAT(alphabet, alphabet, answer; int);\n"
>              ++ "  answer DEL(alphabet, alphabet, answer; int);\n"
>              ++ "  answer PK(alphabet, alphabet, answer, alphabet, alphabet, answer; int);\n"
>              ++ "  answer Lr(alphabet, alphabet, answer, alphabet, alphabet, answer; int);\n"
>              ++ "  answer lR(alphabet, alphabet, answer, alphabet, alphabet, answer; int);\n"
>              ++ "  answer bg(alphabet, alphabet, answer, alphabet, alphabet, answer; int);\n"
>              ++ "  answer jump(Subsequence, answer);\n"
>              ++ "  choice [answer] h([answer]);\n"
>              ++ "}\n"

> -- generates the Bellman's GAP code signature for searching with given structure-consensi as a [TreeList Int]
> gen_build_sig :: String
> gen_build_sig = "signature sig_cm(alphabet, answer) {\n"
>              ++ "  answer INS(alphabet, answer; int);\n"
>              ++ "  answer NIL(void; int);\n"
>              ++ "  answer MAT(alphabet, answer; int);\n"
>              ++ "  answer DEL(answer; int);\n"
>              ++ "  answer PK(alphabet, answer, alphabet, answer; int);\n"
>              ++ "  answer Lr(alphabet, answer,           answer; int);\n"
>              ++ "  answer lR(          answer, alphabet, answer; int);\n"
>              ++ "  answer bg(          answer,           answer; int);\n"
>              ++ "  choice [answer] h([answer]);\n"
>              ++ "}\n"


============== ALGEBRAS ==============

> -- generates the constant Bellman's GAP code algebra to print hints about the parse. Needed because enum algebra has no hints to positions if constant.
> gen_train_alg_parse :: String
> gen_train_alg_parse = "algebra alg_parse implements sig_cm(alphabet = char, answer = Rope) {\n"
>                    ++ "  Rope INS(char a, char p, Rope x; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, \"INS \");\n    append(res, pos);\n    append(res, ' ');\n    append(res, a);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope NIL(void; int pos) {\n    Rope res;\n    append(res, \"NIL \");\n    append(res, pos);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope MAT(char a, char p, Rope x; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, \"MAT \");\n    append(res, pos);\n    append(res, ' ');\n    append(res, a);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope DEL(char a, char p, Rope x; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, \"DEL \");\n    append(res, pos);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope PK(char a, char p, Rope x, char b, char q, Rope y; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, y);\n    append(res, \"PK \");\n    append(res, pos);\n    append(res, ' ');\n    append(res, a);\n    append(res, ' ');\n    append(res, b);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope Lr(char a, char p, Rope x, char b, char q, Rope y; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, y);\n    append(res, \"Lr \");\n    append(res, pos);\n    append(res, ' ');\n    append(res, a);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope lR(char a, char p, Rope x, char b, char q, Rope y; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, y);\n    append(res, \"lR \");\n    append(res, pos);\n    append(res, ' ');\n    append(res, b);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope bg(char a, char p, Rope x, char b, char q, Rope y; int pos) {\n    Rope res;\n    append(res, x);\n    append(res, y);\n    append(res, \"bg \");\n    append(res, pos);\n    append(res, ';');\n    return res;\n  }\n"
>                    ++ "  Rope jump(Subsequence a, Rope x) {\n    return x;\n  }\n"
>                    ++ "  choice [Rope] h([Rope] i) {\n    return i;\n  }\n"
>                    ++ "}\n"

> -- generates the Bellman's GAP code algebra for searching unknown sequences with the CYK algorithm with given structure-consensi as a [TreeList Int]
> gen_build_alg_cyk :: String
> gen_build_alg_cyk = "algebra alg_cyk implements sig_cm(alphabet = char, answer = float) {\n"
>                  ++ "  float NIL(void; int pos) {\n    return getTransition_NIL(pos);\n  }\n"
>                  ++ "  float INS(char a, float x; int pos) {\n    return x + getTransition_INS(pos) + getEmission_INS(pos, a);\n  }\n"
>                  ++ "  float MAT(char a, float x; int pos) {\n    return x + getTransition_MAT(pos) + getEmission_MAT(pos, a);\n  }\n"
>                  ++ "  float DEL(float x; int pos) {\n    return x + getTransition_DEL(pos);\n  }\n"
>                  ++ "  float PK (char a, float x, char b, float y; int pos) {\n    return x + y + getTransition_PK(pos) + getEmission_PK(pos, a,b);\n  }\n"
>                  ++ "  float Lr (char a, float x,         float y; int pos) {\n    return x + y + getTransition_Lr(pos) + getEmission_Lr(pos, a);\n  }\n"
>                  ++ "  float lR (        float x, char b, float y; int pos) {\n    return x + y + getTransition_lR(pos) + getEmission_lR(pos, b);\n  }\n"
>                  ++ "  float bg (        float x,         float y; int pos) {\n    return x + y + getTransition_bg(pos);\n  }\n"
>                  ++ "  choice [float] h([float] i) {\n    return list(maximum(i));\n  }\n"
>                  ++ "}\n"

> -- this is an algebra computing a nonsense score to enable backtracking capabilities for Bellman's GAP, which significantly improves speed of parsing a sequence-structure pair
> gen_train_alg_fake :: String
> gen_train_alg_fake = "algebra alg_fake implements sig_cm(alphabet = char, answer = int) {\n"
>                   ++ "  int INS(char a, char p, int x; int pos) {\n    return x;\n  }\n"
>                   ++ "  int NIL(void; int pos) {\n    return 0;\n  }\n"
>                   ++ "  int MAT(char a, char p, int x; int pos) {\n    return x;\n  }\n"
>                   ++ "  int DEL(char a, char p, int x; int pos) {\n    return x;\n  }\n"
>                   ++ "  int PK(char a, char p, int x, char b, char q, int y; int pos) {\n    return x+y;\n  }\n"
>                   ++ "  int Lr(char a, char p, int x, char b, char q, int y; int pos) {\n    return x+y;\n  }\n"
>                   ++ "  int lR(char a, char p, int x, char b, char q, int y; int pos) {\n    return x+y;\n  }\n"
>                   ++ "  int bg(char a, char p, int x, char b, char q, int y; int pos) {\n    return x+y;\n  }\n"
>                   ++ "  int jump(Subsequence a, int x) {\n    return x;\n  }\n"
>                   ++ "  choice [int] h([int] i) {\n    return list(maximum(i));\n  }\n"
>                   ++ "}\n"

> -- generates the Bellman's GAP code algebra for searching unknown sequences with the CYK algorithm with given structure-consensi as a [TreeList Int]
> gen_build_alg_align = "algebra alg_align implements sig_cm(alphabet = char, answer = ali) {\n"
>                    ++ "  ali NIL(void; int pos) {\n    ali res;\n    return res;  \n}\n"
>                    ++ "  ali INS(char a, ali x; int pos) {\n    ali res;\n    append(res.model, '-'); append(res.seq, a); append(res.cons, getConsensus_INS(pos));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    return res;\n  }\n"
>                    ++ "  ali MAT(char a, ali x; int pos) {\n    ali res;\n    append(res.model, '*'); append(res.seq, a); append(res.cons, getConsensus_MAT(pos));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    return res;\n  }\n"
>                    ++ "  ali DEL(ali x; int pos) {\n    ali res;\n    append(res.model, '*'); append(res.seq, '-'); append(res.cons, getConsensus_DEL(pos));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    return res;\n  }\n"
>                    ++ "  ali PK (char a, ali x, char b, ali y; int pos) {\n    ali res;\n    append(res.model, '<'); append(res.seq, a  ); append(res.cons, getConsensus_PK(pos,0));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    append(res.model, '>'); append(res.seq, b  ); append(res.cons, getConsensus_PK(pos,1));\n    append(res.model, y.model); append(res.seq, y.seq); append(res.cons, y.cons);\n    return res;\n  }\n"
>                    ++ "  ali Lr (char a, ali x,         ali y; int pos) {\n    ali res;\n    append(res.model, '<'); append(res.seq, a  ); append(res.cons, getConsensus_Lr(pos,0));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    append(res.model, '>'); append(res.seq, '-'); append(res.cons, getConsensus_Lr(pos,1));\n    append(res.model, y.model); append(res.seq, y.seq); append(res.cons, y.cons);\n    return res;\n  }\n"
>                    ++ "  ali lR (        ali x, char b, ali y; int pos) {\n    ali res;\n    append(res.model, '<'); append(res.seq, '-'); append(res.cons, getConsensus_lR(pos,0));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    append(res.model, '>'); append(res.seq, b  ); append(res.cons, getConsensus_lR(pos,1));\n    append(res.model, y.model); append(res.seq, y.seq); append(res.cons, y.cons);\n    return res;\n  }\n"
>                    ++ "  ali bg (        ali x,         ali y; int pos) {\n    ali res;\n    append(res.model, '<'); append(res.seq, '-'); append(res.cons, getConsensus_bg(pos,0));\n    append(res.model, x.model); append(res.seq, x.seq); append(res.cons, x.cons);\n    append(res.model, '>'); append(res.seq, '-'); append(res.cons, getConsensus_bg(pos,1));\n    append(res.model, y.model); append(res.seq, y.seq); append(res.cons, y.cons);\n    return res;\n  }\n"
>                    ++ "  choice [ali] h([ali] i) {\n    return i;\n  }\n"
>                    ++ "}\n"


============== GRAMMARS ==============

> -- generates the Bellman's GAP code grammar for training with given structure-consensi as a [TreeList Int]
> gen_train_grammar :: [TreeList (Int,Int)] -> [String]
> gen_train_grammar [El (pos,i)]       = [(gen_train_grammarINS (pos,i)) ++ " | NIL(EMPTY; " ++ (show pos) ++ ") # h;"]
> gen_train_grammar [Pl (pos,i) xs ys] = [(gen_train_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs, y <- ys,
>                                                                                                  let leftNT   = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rightNT  = "a_" ++ (show (snd (getTreeListIndex y))),
>                                                                                                  let leftJump = if ((fst (getTreeListIndex x))-pos == 1) then leftNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-pos-1)) ++ "), " ++ leftNT ++ ")",
>                                                                                                  let rightJump = if ((fst (getTreeListIndex y))-(fst (getTreeListLastIndex x)) == 1) then rightNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex y))-(fst (getTreeListLastIndex x))-1)) ++ "), " ++ rightNT ++ ")",
>                                                                                                  let rules = " | \n\tPK(CHAR with isAnyBase, CHAR('<'), " ++ leftJump ++ ", CHAR with isAnyBase, CHAR('>'), " ++ rightJump ++ "; " ++ me ++ ") | Lr(CHAR with isAnyBase, CHAR('<'), " ++ leftJump ++ ", CHAR('.'), CHAR('>'), " ++ rightJump ++ "; " ++ me ++ ") | lR(CHAR('.'), CHAR('<'), " ++ leftJump ++ ", CHAR with isAnyBase, CHAR('>'), " ++ rightJump ++ "; " ++ me ++ ") | bg(CHAR('.'), CHAR('<'), " ++ leftJump ++ ", CHAR('.'), CHAR('>'), " ++ rightJump ++ "; " ++ me ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_train_grammar xs) ++ (gen_train_grammar ys)
>   where me    = show pos
> gen_train_grammar [Ol (pos,i) xs]    = [(gen_train_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs,
>                                                                                                  let subNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let jump = if ((fst (getTreeListIndex x))-pos == 1) then subNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-pos-1)) ++ "), " ++ subNT ++ ")",
>                                                                                                  let rules = " | MAT(CHAR with isAnyBase, CHAR('*'), " ++ jump ++ "; " ++ me ++ ") | DEL(CHAR('.'), CHAR('*'), " ++ jump ++ "; " ++ me ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_train_grammar xs)
>   where me    = show pos
> gen_train_grammar (a:as)             = (gen_train_grammar [a]) ++ (gen_train_grammar as)
> gen_train_grammar _                  = error ("gen_train_grammar: this case should never ever been reached")
> gen_train_grammarINS (pos,i)         = "a_" ++ (show i) ++ " = INS(CHAR with isAnyBase, CHAR('-'), a_" ++ (show i) ++ "; " ++ (show pos) ++ ")"


> -- generates the Bellman's GAP code grammar for training with given structure-consensi as a [TreeList Int]
> gen_build_grammar :: [TreeList (Int,Int)] -> [String]
> gen_build_grammar [El (pos,i)]       = [(gen_build_grammarINS (pos,i)) ++ " | NIL(EMPTY; " ++ (show pos) ++ ") # h;"]
> gen_build_grammar [Pl (pos,i) xs ys] = [(gen_build_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs, y <- ys,
>                                                                                                  let leftNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rightNT = "a_" ++ (show (snd (getTreeListIndex y))),
>                                                                                                  let rules = " | PK(CHAR, " ++ leftNT ++ ", CHAR, " ++ rightNT ++ "; " ++ me ++ ") | Lr(CHAR, " ++ leftNT ++ ", " ++ rightNT ++ "; " ++ me ++ ") | lR(" ++ leftNT ++ ", CHAR, " ++ rightNT ++ "; " ++ me ++ ") | bg(" ++ leftNT ++ ", " ++ rightNT ++ "; " ++ me ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_build_grammar xs) ++ (gen_build_grammar ys)
>   where me = show pos
> gen_build_grammar [Ol (pos,i) xs]    = [(gen_build_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs,
>                                                                                                  let subNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rules = " | MAT(CHAR, " ++ subNT ++ "; " ++ me ++ ") | DEL(" ++ subNT ++ "; " ++ me ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_build_grammar xs)
>   where me = show pos
> gen_build_grammar (a:as)             = (gen_build_grammar [a]) ++ (gen_build_grammar as)
> gen_build_grammar _                  = error ("gen_build_grammar: this case should never ever been reached")
> gen_build_grammarINS (pos,i)         = "a_" ++ (show i) ++ " = INS(CHAR, a_" ++ (show i) ++ "; " ++ (show pos) ++ ")"








> -- just my toy example
> ex :: [TreeList Int] -- "<<<<<***------------------------***>>>>>" "<<<<**<<<<<<<<<****>>>>>*>>>>**-****>>>>"
> ex = [Pl 1 [Pl 2 [Pl 3 [Pl 4 [Pl 5 [Ol 6 [Ol 7 [Ol 8 [Ol 33 [Ol 34 [Ol 35 [El 36]]]]]]] [El 37],Ol 5 [Ol 6 [Pl 7 [Pl 8 [Pl 9 [Pl 10 [Pl 11 [Pl 12 [Pl 13 [Pl 14 [Pl 15 [Ol 16 [Ol 17 [Ol 18 [Ol 19 [El 20]]]]] [El 21]] [El 22]] [El 23]] [El 24]] [Ol 25 [El 26]]] [El 27]] [El 28]] [El 29]] [Ol 30 [Ol 31 [Ol 33 [Ol 34 [Ol 35 [Ol 36 [El 37]]]]]]]]]] [El 38]] [El 39]] [El 40]] [El 41]]
