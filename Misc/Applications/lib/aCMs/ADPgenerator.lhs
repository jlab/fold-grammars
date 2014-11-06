> module ADPgenerator where

> import ADPCombinators
> import Data.Array
> import Data.List
> import Datatypes

The signature:

> data X 
>   =  KeepPair  Int X X
>   |  DelPair   Int X X
>   |  InsPair   Int X X
>   |  MakePair  Int X X
>   |  BreakPair Int X X
>   |  DelLeft   Int X X
>   |  DelRight  Int X X
>   |  FindLeft  Int X X
>   |  FindRight Int X X
>   |  KeepOpen  Int X
>   |  DelOpen   Int X
>   |  InsOpen   Int X
>   |  KeepSkip  Int X
>   |  Nil       Int
>                        deriving (Eq, Show)



Algebra type:

> type Algebra alphabet answer = (
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> answer,
>   alphabet -> alphabet -> answer -> answer,
>   Int -> answer,
>   [answer] -> [answer]
>   )

Counting algebra:

> count :: Algebra Int Int
> count = (keeppair, delpair, inspair, makepair, breakpair, delleft, delright, findleft, findright, keepopen, delopen, insopen, keepskip, nil, h) where
>   keeppair  a b c d e f = c * f
>   delpair   a b c d e f = c * f
>   inspair   a b c d e f = c * f
>   makepair  a b c d e f = c * f
>   breakpair a b c d e f = c * f
>   delleft   a b c d e f = c * f
>   delright  a b c d e f = c * f
>   findleft  a b c d e f = c * f
>   findright a b c d e f = c * f
>   keepopen  a b c       = c
>   delopen   a b c       = c
>   insopen   a b c       = c
>   keepskip  a b c       = c
>   nil i                 = 1
>   h []                  = []
>   h xs                  = [sum xs]

Enumeration algebra:

> enum :: Algebra Int X
> enum = (keeppair, delpair, inspair, makepair, breakpair, delleft, delright, findleft, findright, keepopen, delopen, insopen, keepskip, nil, h) where
>   keeppair  a b c d e f = KeepPair  (div b 2) c f
>   delpair   a b c d e f = DelPair   (div b 2) c f
>   inspair   a b c d e f = InsPair   (div b 2) c f
>   makepair  a b c d e f = MakePair  (div b 2) c f
>   breakpair a b c d e f = BreakPair (div b 2) c f
>   delleft   a b c d e f = DelLeft   (div b 2) c f
>   delright  a b c d e f = DelRight  (div b 2) c f
>   findleft  a b c d e f = FindLeft  (div b 2) c f
>   findright a b c d e f = FindRight (div b 2) c f
>   keepopen  a b c       = KeepOpen  (div b 2) c
>   delopen   a b c       = DelOpen   (div b 2) c
>   insopen   a b c       = InsOpen   (div b 2) c
>   keepskip  a b c       = KeepSkip  (div b 2) c
>   nil i                 = Nil      ((div i 2)+1)
>   h                     = id

> gen :: Algebra Int (Tree Int, Tree Int)
> gen = (keeppair, delpair, inspair, makepair, breakpair, delleft, delright, findleft, findright, keepopen, delopen, insopen, keepskip, nil, h) where
>   keeppair  _ i (x1,x2) _ _ (y1,y2) = (P (div i 2) x1 y1, P (div i 2) x2 y2)
>   delpair   _ i (x1,x2) _ _ (y1,y2) = (P (div i 2) x1 y1, x2 +:+ y2        )
>   inspair   _ i (x1,x2) _ _ (y1,y2) = (x1 +:+ y1        , P (div i 2) x2 y2)
>   makepair  _ i (x1,x2) _ _ (y1,y2) = (O (div i 2) (x1 +:+ (O i' (E (i'+1))) +:+ y1), P (div i 2) x2 y2)
>                            where i' = (getIndex(getEnd x1))
>   breakpair _ i (x1,x2) _ _ (y1,y2) = (P (div i 2) x1 y1, O (div i 2) (x2 +:+ (O i' (E (i'+1))) +:+ y2))
>                            where i' = (getIndex(getEnd x2))
>   delleft   _ i (x1,x2) _ _ (y1,y2) = (P (div i 2) x1 y1, x2 +:+ (O i' y2))
>                            where i' = (getIndex(getEnd x2))
>   delright  _ i (x1,x2) _ _ (y1,y2) = (P (div i 2) x1 y1, O (div i 2) (x2 +:+ y2))
>   findleft  _ i (x1,x2) _ _ (y1,y2) = (x1 +:+ (O i' y1), P (div i 2) x2 y2)
>                            where i' = (getIndex(getEnd x1))
>   findright _ i (x1,x2) _ _ (y1,y2) = (O (div i 2) (x1 +:+ y1), P (div i 2) x2 y2)
>   keepopen  _ i (x1,x2)             = (O (div i 2) x1   , O (div i 2) x2   )
>   delopen   _ i (x1,x2)             = (O (div i 2) x1   , x2               )
>   insopen   _ i (x1,x2)             = (x1               , O (div i 2) x2   )
>   keepskip  _ i (x1,x2)             = (x1               , x2               )
>   nil i                             = (E ((div i 2)+1)  , E ((div i 2)+1)  )
>   h                                 = id

Algebra product operation:

> (***) :: Eq answer1 => 
>    (Algebra Int answer1) -> 
>    (Algebra Int answer2) -> 
>    Algebra Int (answer1,answer2)
> infix ***
> alg1 *** alg2 = (keeppair, delpair, inspair, makepair, breakpair, delleft, delright, findleft, findright, keepopen, delopen, insopen, keepskip, nil, h) where
>   (keeppair1, delpair1, inspair1, makepair1, breakpair1, delleft1, delright1, findleft1, findright1, keepopen1, delopen1, insopen1, keepskip1, nil1, h1) = alg1
>   (keeppair2, delpair2, inspair2, makepair2, breakpair2, delleft2, delright2, findleft2, findright2, keepopen2, delopen2, insopen2, keepskip2, nil2, h2) = alg2
> 
>   keeppair  a b (c1,c2) d e (f1,f2) = (keeppair1  a b c1 d e f1, keeppair2  a b c2 d e f2)
>   delpair   a b (c1,c2) d e (f1,f2) = (delpair1   a b c1 d e f1, delpair2   a b c2 d e f2)
>   inspair   a b (c1,c2) d e (f1,f2) = (inspair1   a b c1 d e f1, inspair2   a b c2 d e f2)
>   makepair  a b (c1,c2) d e (f1,f2) = (makepair1  a b c1 d e f1, makepair2  a b c2 d e f2)
>   breakpair a b (c1,c2) d e (f1,f2) = (breakpair1 a b c1 d e f1, breakpair2 a b c2 d e f2)
>   delleft   a b (c1,c2) d e (f1,f2) = (delleft1   a b c1 d e f1, delleft2   a b c2 d e f2)
>   delright  a b (c1,c2) d e (f1,f2) = (delright1  a b c1 d e f1, delright2  a b c2 d e f2)
>   findleft  a b (c1,c2) d e (f1,f2) = (findleft1  a b c1 d e f1, findleft2  a b c2 d e f2)
>   findright a b (c1,c2) d e (f1,f2) = (findright1 a b c1 d e f1, findright2 a b c2 d e f2)
>   keepopen  a b (c1,c2)             = (keepopen1  a b c1       , keepopen2  a b c2       )
>   delopen   a b (c1,c2)             = (delopen1   a b c1       , delopen2   a b c2       )
>   insopen   a b (c1,c2)             = (insopen1   a b c1       , insopen2   a b c2       )
>   keepskip  a b (c1,c2)             = (keepskip1  a b c1       , keepskip2  a b c2       )
>   nil i                             = (nil1       i            , nil2       i            )
> 
>   h xs = [(x1,x2)| x1 <- nub $ h1 [ y1 | (y1,y2) <- xs],
>                    x2 <-       h2 [ y2 | (y1,y2) <- xs, y1 == x1]]

The yield grammar:

> -- input for "grammar" is a tuple of "Consensus", i.e. ((String, String), (String, String)), which is transformed to one string holding an interlace of both consensus structures. Btw. it checks if both structures have the same length.
> grammar alg (a,b) = axiom x where
>   (keeppair, delpair, inspair, makepair, breakpair, delleft, delright, findleft, findright, keepopen, delopen, insopen, keepskip, nil, h) = alg
> 
>   x = tabulated(
>       keeppair  <<< char '<' -~~ char '<' ~~! x ~~- char '>' ~~- char '>' ~~!! x ||| 
>       delpair   <<< char '<' -~~ char '-' ~~! x ~~- char '>' ~~- char '-' ~~!! x ||| 
>       inspair   <<< char '-' -~~ char '<' ~~! x ~~- char '-' ~~- char '>' ~~!! x ||| 
>       makepair  <<< char '*' -~~ char '<' ~~! x ~~- char '*' ~~- char '>' ~~!! x ||| 
>       breakpair <<< char '<' -~~ char '*' ~~! x ~~- char '>' ~~- char '*' ~~!! x ||| 
>       delleft   <<< char '<' -~~ char '-' ~~! x ~~- char '>' ~~- char '*' ~~!! x ||| 
>       delright  <<< char '<' -~~ char '*' ~~! x ~~- char '>' ~~- char '-' ~~!! x ||| 
>       findleft  <<< char '-' -~~ char '<' ~~! x ~~- char '*' ~~- char '>' ~~!! x ||| 
>       findright <<< char '*' -~~ char '<' ~~! x ~~- char '-' ~~- char '>' ~~!! x ||| 
>       keepopen  <<< char '*' -~~ char '*' ~~! x                                  ||| 
>       delopen   <<< char '*' -~~ char '-' ~~! x                                  ||| 
>       insopen   <<< char '-' -~~ char '*' ~~! x                                  ||| 
>       keepskip  <<< char '-' -~~ char '-' ~~! x                                  ||| 
>       nil       <<< empti                                                        ... h)
>          where
>          infixl 7 ~~!
>          (~~!)  = (~~*) (2,2) 0
>          infixl 7 ~~!!
>          (~~!!) = (*~*) 4 0


Bind input:

>   z         = if (length (snd a) == length (snd b)) then 
>                  mk (foldr1 (++) (map (\(x,y) -> x:[y]) (zip (snd a) (snd b)))) 
>               else 
>                  error("Consensus structures '" ++ (fst a) ++ "' and '" ++ (fst b) ++ "' have different lengths. They are not part of the same alignment and thus cannot build a multiple covariance model!")
>   (_,n)     = bounds z
>   achar     = achar' z
>   axiom     = axiom' n
>   tabulated = table n
>   char      = charp' z

> empti        :: Parser Int
> empti  (i,j) =  [i | i == j]
> charp'           ::  Eq a => Array Int a -> a -> Parser Int
> charp' z c (i,j) =  [j | i+1 == j, z!j == c]

> -- given a list of n elements, the function getPairs returns a list of all (n over 2) pairs of elements
> getPairs :: [a] -> [(a,a)]
> getPairs [] = []
> getPairs [x] = []
> getPairs (x:xs) = [(x,y) | y <- xs] ++ getPairs xs



> checkCompatibility structures = warning
>    where warning = if (length namesIncPairs > 0) then
>                       error ("Structurs '"++(fst (head namesIncPairs))++"' and '"++(snd (head namesIncPairs))++"', which you provided in the input alignment file, are incompatible to each other. Please correct them and try again.")
>                    else
>                       ""
>          namesIncPairs = map (\(((a,_),(b,_)),_) -> (a,b)) incompatiblePairs     -- reduce information of an invalid pair to the pair of its structure names
>          incompatiblePairs = filter (\(_,c) -> c == []) (zip pairs checkResults) -- filter out all pairs of structures, where the parsing results in no valid candidates
>          checkResults = map (grammar count) pairs                                -- check if a pair of structures is compatible as the result of parsing both structures with "grammar" and count the number of valid candidates. Do this for all pairs via map.
>          pairs = getPairs structures                                             -- create all n over 2 possible pairs of structures
