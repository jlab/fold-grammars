> module Treelist2cm where

> import Data.List
> import Datatypes
> import Probs
> import Data.Map (Map)
> import qualified Data.Map as Map

> get_train_code :: [TreeList Int] -> String
> get_train_code consensi = "import isntimes\n\n" ++
>                           signature ++ "\n" ++
>                           "algebra alg_count auto count;\n" ++
>                           "algebra alg_enum auto enum;\n\n" ++
>                           grammar ++ 
>                           "instance count = gra_train(alg_count);\n" ++
>                           "instance train = gra_train(alg_enum);\n"
>         where signature = "signature sig_cm(alphabet, answer) {\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_train_sig (collectAlgfkts indexedConsensi))) ++ "\n" ++
>                           "  answer jump(Subsequence, answer);\n" ++ 
>                           "  choice [answer] h([answer]);\n" ++
>                           "}\n"
>               grammar   = "grammar gra_train uses sig_cm(axiom = start) {\n" ++
>                           "  base = CHAR('A') | CHAR('C') | CHAR('G') | CHAR('U') | CHAR('R') | CHAR('Y') | CHAR('M') | CHAR('K') | CHAR('W') | CHAR('S') | CHAR('B') | CHAR('D') | CHAR('H') | CHAR('V') | CHAR('N');\n" ++
>                           "  start = " ++ (foldr1 (\x -> \y -> x ++ " | " ++ y) [ jump | x <- indexedConsensi, 
>                                                                                           let subNT = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                           let jump = if ((fst (getTreeListIndex x)) == 1) then subNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-1)) ++ "), " ++ subNT ++ ")"
>                                                                                 ]) ++ " # h;\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_train_grammar indexedConsensi)) ++ "\n" ++
>                           "}\n"
>               indexedConsensi = indexTreeList consensi 1

> get_build_code :: [TreeList Int] -> String
> get_build_code consensi = "type Rope = extern\n" ++
>                           "type ali = (Rope model, Rope seq)\n\n" ++
>                           signature ++
>                           "algebra alg_count auto count;\n\n" ++
>                           "algebra alg_enum auto enum;\n\n" ++
>                           alg_cyk ++ "\n" ++
>                           alg_align ++ "\n" ++
>                           grammar ++ "\n" ++
>                           "instance count = gra_build(alg_count);\n" ++
>                           "instance train = gra_build(alg_enum);\n" ++
>                           "instance cyk = gra_build(alg_cyk);\n" ++
>                           "instance cykali = gra_build(alg_cyk * alg_align);\n"
>         where signature = "signature sig_cm(alphabet, answer) {\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_build_sig (collectAlgfkts indexedConsensi))) ++ "\n" ++
>                           "  choice [answer] h([answer]);\n" ++
>                           "}\n"
>               grammar   = "grammar gra_build uses sig_cm(axiom = start) {\n" ++
>                           "  start = " ++ (foldr1 (\x -> \y -> x ++ " | " ++ y) [ "a_" ++ (show (snd (getTreeListIndex x))) | x <- indexedConsensi ]) ++ " # h;\n" ++
>                           "  " ++ (foldl1 (\x y -> x ++ "\n  " ++ y) (gen_build_grammar indexedConsensi)) ++ "\n" ++
>                           "}\n"
>               indexedConsensi = indexTreeList consensi 1
>               alg_cyk   = "algebra alg_cyk implements sig_cm(alphabet = char, answer = float) {\n" ++
>                            (foldl1 (\x y -> x ++ y) (gen_build_alg_cyk (collectAlgfkts indexedConsensi))) ++ "\n" ++
>                           "  choice [float] h([float] i) {\n    return list(maximum(i));\n  }\n" ++
>                           "}\n"
>               alg_align = "algebra alg_align implements sig_cm(alphabet = char, answer = ali) {\n" ++
>                            (foldl1 (\x y -> x ++ y) (gen_build_alg_align (collectAlgfkts indexedConsensi))) ++ "\n" ++
>                           "  choice [ali] h([ali] i) {\n    return i;\n  }\n" ++
>                           "}\n"


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

> -- generates the Bellman's GAP code signature for training with given structure-consensi as a [TreeList Int]
> gen_train_sig :: [(String, Int)] -> [String]
> gen_train_sig []           = []
> gen_train_sig [("NIL", i)] = ["answer NIL_" ++ (show i) ++ "(void);"]
> gen_train_sig [("INS", i)] = ["answer INS_" ++ (show i) ++ "(alphabet, alphabet, answer);"]
> gen_train_sig [("MAT", i)] = ["answer MAT_" ++ (show i) ++ "(alphabet, alphabet, answer);"]
> gen_train_sig [("DEL", i)] = ["answer DEL_" ++ (show i) ++ "(alphabet, alphabet, answer);"]
> gen_train_sig [("PK",  i)] = ["answer PK_" ++ (show i) ++ "(alphabet, alphabet, answer, alphabet, alphabet, answer);"]
> gen_train_sig [("Lr",  i)] = ["answer Lr_" ++ (show i) ++ "(alphabet, alphabet, answer, alphabet, alphabet, answer);"]
> gen_train_sig [("lR",  i)] = ["answer lR_" ++ (show i) ++ "(alphabet, alphabet, answer, alphabet, alphabet, answer);" ]
> gen_train_sig [("bg",  i)] = ["answer bg_" ++ (show i) ++ "(alphabet, alphabet, answer, alphabet, alphabet, answer);"]
> gen_train_sig (a:as)       = (gen_train_sig [a]) ++ (gen_train_sig as)

> -- generates the Bellman's GAP code signature for searching with given structure-consensi as a [TreeList Int]
> gen_build_sig :: [(String, Int)] -> [String]
> gen_build_sig []           = []
> gen_build_sig [("NIL", i)] = ["answer NIL_" ++ (show i) ++ "(void);"]
> gen_build_sig [("INS", i)] = ["answer INS_" ++ (show i) ++ "(alphabet, answer);"]
> gen_build_sig [("MAT", i)] = ["answer MAT_" ++ (show i) ++ "(alphabet, answer);"]
> gen_build_sig [("DEL", i)] = ["answer DEL_" ++ (show i) ++ "(answer);"]
> gen_build_sig [("PK",  i)] = ["answer PK_" ++ (show i) ++ "(alphabet, answer, alphabet, answer);"]
> gen_build_sig [("Lr",  i)] = ["answer Lr_" ++ (show i) ++ "(alphabet, answer,           answer);"]
> gen_build_sig [("lR",  i)] = ["answer lR_" ++ (show i) ++ "(          answer, alphabet, answer);"]
> gen_build_sig [("bg",  i)] = ["answer bg_" ++ (show i) ++ "(          answer,           answer);"]
> gen_build_sig (a:as)       = (gen_build_sig [a]) ++ (gen_build_sig as)


============== ALGEBRAS ==============

> -- generates the Bellman's GAP code algebra for searching unknown sequences with the CYK algorithm with given structure-consensi as a [TreeList Int]
> gen_build_alg_cyk :: [(String, Int)] -> [String]
> gen_build_alg_cyk []           = []
> gen_build_alg_cyk [("NIL", i)] = ["  float NIL_" ++ (show i) ++ "(void) {\n    return " ++ (getTransition ("NIL_" ++ (show i))) ++ "+0.0;\n  }\n"]
> gen_build_alg_cyk [("INS", i)] = ["  float INS_" ++ (show i) ++ "(char a, float x) {\n" ++ (printUnpair ("INS_" ++ (show i)) "a") ++ (getTransition ("INS_" ++ (show i))) ++ "+x;\n  }\n"]
> gen_build_alg_cyk [("MAT", i)] = ["  float MAT_" ++ (show i) ++ "(char a, float x) {\n" ++ (printUnpair ("MAT_" ++ (show i)) "a") ++ (getTransition ("MAT_" ++ (show i))) ++ "+x;\n  }\n"]
> gen_build_alg_cyk [("DEL", i)] = ["  float DEL_" ++ (show i) ++ "(float x) {\n    return " ++ (getTransition ("DEL_" ++ (show i))) ++ "+x;\n  }\n"]
> gen_build_alg_cyk [("PK",  i)] = ["  float PK_"  ++ (show i) ++ "(char a, float x, char b, float y) {\n" ++ (printPair ("PK_" ++ (show i))) ++ "    return " ++ (getTransition ("PK_" ++ (show i))) ++ "+e+x+y;\n  }\n"]
> gen_build_alg_cyk [("Lr",  i)] = ["  float Lr_"  ++ (show i) ++ "(char a, float x,         float y) {\n" ++ (printUnpair ("Lr_" ++ (show i)) "a") ++ (getTransition ("Lr_" ++ (show i)))  ++ "+x+y;\n  }\n"]
> gen_build_alg_cyk [("lR",  i)] = ["  float lR_"  ++ (show i) ++ "(        float x, char b, float y) {\n" ++ (printUnpair ("lR_" ++ (show i)) "b") ++ (getTransition ("lR_" ++ (show i)))  ++ "+x+y;\n  }\n"]
> gen_build_alg_cyk [("bg",  i)] = ["  float bg_"  ++ (show i) ++ "(        float x,         float y) {\n    return " ++ (getTransition ("bg_" ++ (show i)))  ++ "+x+y;\n  }\n"]
> gen_build_alg_cyk (a:as)       = (gen_build_alg_cyk [a]) ++ (gen_build_alg_cyk as)

> printUnpair :: String -> String -> String
> printUnpair id baseName = if (test) then "    return " else "    float e = 0.0;\n  " ++ (mergeLines ["  if (" ++ baseName ++ " == '" ++ x ++ "') e = " ++ (getEmission id x) ++ ";" | x <- alphabet]) ++ "\n    return e+"
>   where test = ["0","0","0","0"] == [ getEmission id a | a <- alphabet ]

> printPair :: String -> String
> printPair id = "    float e = 0.0;\n  " ++ (mergeLines ["  if ((a == '" ++ x ++ "') && (b == '" ++ y ++ "')) e = " ++ (getEmission id (x++y)) ++ ";" | x <- alphabet, y <- alphabet]) ++ "\n"

> mergeLines :: [String] -> String
> mergeLines xs = foldl1 (\x y -> x ++ "\n  " ++ y) xs

> alphabet :: [String] 
> alphabet = ["A","C","G","U"]

> pairs :: [String]
> pairs = [ x++y | x <- alphabet, y <- alphabet]


> -- generates the Bellman's GAP code algebra for searching unknown sequences with the CYK algorithm with given structure-consensi as a [TreeList Int]
> gen_build_alg_align :: [(String, Int)] -> [String]
> gen_build_alg_align []           = []
> gen_build_alg_align [("NIL", i)] = ["  ali NIL_" ++ (show i) ++ "(void) {    ali res;\n    return res;\n}\n"]
> gen_build_alg_align [("INS", i)] = ["  ali INS_" ++ (show i) ++ "(char a, ali x) {    ali res;\n    append(res.model, '-'); append(res.seq, a);\n    append(res.model, x.model); append(res.seq, x.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("MAT", i)] = ["  ali MAT_" ++ (show i) ++ "(char a, ali x) {    ali res;\n    append(res.model, '*'); append(res.seq, a);\n    append(res.model, x.model); append(res.seq, x.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("DEL", i)] = ["  ali DEL_" ++ (show i) ++ "(ali x) {    ali res;\n    append(res.model, '*'); append(res.seq, '-');\n    append(res.model, x.model); append(res.seq, x.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("PK",  i)] = ["  ali PK_"  ++ (show i) ++ "(char a, ali x, char b, ali y) {    ali res;\n    append(res.model, '<'); append(res.seq, a  );\n    append(res.model, x.model); append(res.seq, x.seq);\n    append(res.model, '>'); append(res.seq, b  );\n    append(res.model, y.model); append(res.seq, y.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("Lr",  i)] = ["  ali Lr_"  ++ (show i) ++ "(char a, ali x,         ali y) {    ali res;\n    append(res.model, '<'); append(res.seq, a  );\n    append(res.model, x.model); append(res.seq, x.seq);\n    append(res.model, '>'); append(res.seq, '-');\n    append(res.model, y.model); append(res.seq, y.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("lR",  i)] = ["  ali lR_"  ++ (show i) ++ "(        ali x, char b, ali y) {    ali res;\n    append(res.model, '<'); append(res.seq, '-');\n    append(res.model, x.model); append(res.seq, x.seq);\n    append(res.model, '>'); append(res.seq, b  );\n    append(res.model, y.model); append(res.seq, y.seq);\n    return res;\n  }\n"]
> gen_build_alg_align [("bg",  i)] = ["  ali bg_"  ++ (show i) ++ "(        ali x,         ali y) {    ali res;\n    append(res.model, '<'); append(res.seq, '-');\n    append(res.model, x.model); append(res.seq, x.seq);\n    append(res.model, '>'); append(res.seq, '-');\n    append(res.model, y.model); append(res.seq, y.seq);\n    return res;\n  }\n"]
> gen_build_alg_align (a:as)       = (gen_build_alg_align [a]) ++ (gen_build_alg_align as)


============== GRAMMARS ==============

> -- generates the Bellman's GAP code grammar for training with given structure-consensi as a [TreeList Int]
> gen_train_grammar :: [TreeList (Int,Int)] -> [String]
> gen_train_grammar [El (pos,i)]       = [(gen_train_grammarINS (pos,i)) ++ " | NIL_" ++ (show pos) ++ "(EMPTY) # h;"]
> gen_train_grammar [Pl (pos,i) xs ys] = [(gen_train_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs, y <- ys,
>                                                                                                  let leftNT   = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rightNT  = "a_" ++ (show (snd (getTreeListIndex y))),
>                                                                                                  let leftJump = if ((fst (getTreeListIndex x))-pos == 1) then leftNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-pos-1)) ++ "), " ++ leftNT ++ ")",
>                                                                                                  let rightJump = if ((fst (getTreeListIndex y))-(fst (getTreeListLastIndex x)) == 1) then rightNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex y))-(fst (getTreeListLastIndex x))-1)) ++ "), " ++ rightNT ++ ")",
>                                                                                                  let rules = " | \n\tPK_" ++ me ++ "(base, CHAR('<'), " ++ leftJump ++ ", base, CHAR('>'), " ++ rightJump ++ ") | Lr_" ++ me ++ "(base, CHAR('<'), " ++ leftJump ++ ", CHAR('.'), CHAR('>'), " ++ rightJump ++ ") | lR_" ++ me ++ "(CHAR('.'), CHAR('<'), " ++ leftJump ++ ", base, CHAR('>'), " ++ rightJump ++ ") | bg_" ++ me ++ "(CHAR('.'), CHAR('<'), " ++ leftJump ++ ", CHAR('.'), CHAR('>'), " ++ rightJump ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_train_grammar xs) ++ (gen_train_grammar ys)
>   where me    = show pos
> gen_train_grammar [Ol (pos,i) xs]    = [(gen_train_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs,
>                                                                                                  let subNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let jump = if ((fst (getTreeListIndex x))-pos == 1) then subNT else "jump(REGION0 with isntimes('_'," ++ (show ((fst (getTreeListIndex x))-pos-1)) ++ "), " ++ subNT ++ ")",
>                                                                                                  let rules = " | MAT_" ++ me ++ "(base, CHAR('*'), " ++ jump ++ ") | DEL_" ++ me ++ "(CHAR('.'), CHAR('*'), " ++ jump ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_train_grammar xs)
>   where me    = show pos
> gen_train_grammar (a:as)             = (gen_train_grammar [a]) ++ (gen_train_grammar as)
> gen_train_grammar _                  = error ("gen_train_grammar: this case should never ever been reached")
> gen_train_grammarINS (pos,i)         = "a_" ++ (show i) ++ " = INS_" ++ (show pos) ++ "(base, CHAR('-'), a_" ++ (show i) ++ ")"


> -- generates the Bellman's GAP code grammar for training with given structure-consensi as a [TreeList Int]
> gen_build_grammar :: [TreeList (Int,Int)] -> [String]
> gen_build_grammar [El (pos,i)]       = [(gen_build_grammarINS (pos,i)) ++ " | NIL_" ++ (show pos) ++ "(EMPTY) # h;"]
> gen_build_grammar [Pl (pos,i) xs ys] = [(gen_build_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs, y <- ys,
>                                                                                                  let leftNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rightNT = "a_" ++ (show (snd (getTreeListIndex y))),
>                                                                                                  let rules = " | PK_" ++ me ++ "(CHAR, " ++ leftNT ++ ", CHAR, " ++ rightNT ++ ") | Lr_" ++ me ++ "(CHAR, " ++ leftNT ++ ", " ++ rightNT ++ ") | lR_" ++ me ++ "(" ++ leftNT ++ ", CHAR, " ++ rightNT ++ ") | bg_" ++ me ++ "(" ++ leftNT ++ ", " ++ rightNT ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_build_grammar xs) ++ (gen_build_grammar ys)
>   where me = show pos
> gen_build_grammar [Ol (pos,i) xs]    = [(gen_build_grammarINS (pos,i)) ++ (foldr1 (++) [ rules | x <- xs,
>                                                                                                  let subNT  = "a_" ++ (show (snd (getTreeListIndex x))),
>                                                                                                  let rules = " | MAT_" ++ me ++ "(CHAR, " ++ subNT ++ ") | DEL_" ++ me ++ "(" ++ subNT ++ ")"
>                                                                                        ]) ++ " # h;"] ++ (gen_build_grammar xs)
>   where me = show pos
> gen_build_grammar (a:as)             = (gen_build_grammar [a]) ++ (gen_build_grammar as)
> gen_build_grammar _                  = error ("gen_build_grammar: this case should never ever been reached")
> gen_build_grammarINS (pos,i)         = "a_" ++ (show i) ++ " = INS_" ++ (show pos) ++ "(CHAR, a_" ++ (show i) ++ ")"








> -- just my toy example
> ex :: [TreeList Int] -- "<<<<<***------------------------***>>>>>" "<<<<**<<<<<<<<<****>>>>>*>>>>**-****>>>>"
> ex = [Pl 1 [Pl 2 [Pl 3 [Pl 4 [Pl 5 [Ol 6 [Ol 7 [Ol 8 [Ol 33 [Ol 34 [Ol 35 [El 36]]]]]]] [El 37],Ol 5 [Ol 6 [Pl 7 [Pl 8 [Pl 9 [Pl 10 [Pl 11 [Pl 12 [Pl 13 [Pl 14 [Pl 15 [Ol 16 [Ol 17 [Ol 18 [Ol 19 [El 20]]]]] [El 21]] [El 22]] [El 23]] [El 24]] [Ol 25 [El 26]]] [El 27]] [El 28]] [El 29]] [Ol 30 [Ol 31 [Ol 33 [Ol 34 [Ol 35 [Ol 36 [El 37]]]]]]]]]] [El 38]] [El 39]] [El 40]] [El 41]]
