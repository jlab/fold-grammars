Haskell header:

> module SingleTree where

> import ADPCombinators
> import Data.Array
> import Data.List
> import Datatypes

The signature:

> data S
>   =  Nil  Int 
>   |  Pair Int S S
>   |  Open Int S
>   |  Skip Int S
>                        deriving (Eq, Show)



Algebra type:

> type Algebra alphabet answer = (
>   Int -> answer,
>   alphabet -> answer -> alphabet -> answer -> answer,
>   alphabet -> answer -> answer,
>   alphabet -> answer -> answer,
>   [answer] -> [answer]
>   )

Enumeration algebra:

> enum :: Algebra Int S
> enum = (nil, pair, open, skip, h) where
>   nil i         = Nil (i+1)
>   pair a x a' y = Pair a x y
>   open a x      = Open a x
>   skip a x      = Skip a x
>   h             = id

> alg_tree :: Algebra Int (Tree Int)
> alg_tree = (nil, pair, open, skip, h) where
>   nil i         = E (i+1)
>   pair a x a' y = P a x y
>   open a x      = O a x
>   skip a x      = x
>   h             = id

Counting algebra:

> count :: Algebra Int Int
> count = (nil, pair, open, skip, h) where
>   nil a        = 1
>   pair a b c d = b * d
>   open a b     = b
>   skip a b     = b
>   h []         = []
>   h xs         = [sum xs]

Algebra product operation:

> (***) :: Eq answer1 => 
>    (Algebra Int answer1) -> 
>    (Algebra Int answer2) -> 
>    Algebra Int (answer1,answer2)
> infix ***
> alg1 *** alg2 = (nil, pair, open, skip, h) where
>   (nil1, pair1, open1, skip1, h1) = alg1
>   (nil2, pair2, open2, skip2, h2) = alg2
> 
>   nil a                    = (nil1 a, nil2 a)
>   pair a (b1,b2) c (d1,d2) = (pair1 a b1 c d1, pair2 a b2 c d2)
>   open a (b1,b2)           = (open1 a b1, open2 a b2)
>   skip a (b1,b2)           = (skip1 a b1, skip2 a b2)
> 
>   h xs = [(x1,x2)| x1 <- nub $ h1 [ y1 | (y1,y2) <- xs],
>                    x2 <-       h2 [ y2 | (y1,y2) <- xs, y1 == x1]]

The yield grammar:

> gra_g5skip alg f = axiom s where
>   (nil, pair, open, skip, h) = alg
> 
>   s = tabulated(
>       nil  <<< empti                             ||| 
>       pair <<< char '<' -~~ s ~~- char '>' ~~~ s ||| 
>       open <<< char '*' -~~ s                    ||| 
>       skip <<< char '-' -~~ s                    ... h)


Bind input:

>   z         = mk f
>   (_,n)     = bounds z
>   achar     = achar' z
>   axiom     = axiom' n
>   tabulated = table n
>   char      = charp' z


> empti        :: Parser Int
> empti  (i,j) =  [i | i == j]
> charp'           ::  Eq a => Array Int a -> a -> Parser Int
> charp' z c (i,j) =  [j | i+1 == j, z!j == c]

