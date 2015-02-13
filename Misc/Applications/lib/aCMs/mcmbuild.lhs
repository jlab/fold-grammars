> module Main where

> import System.Environment(getArgs,getProgName)

Haskell header:

> import Datatypes
> import Treelist2cm
> import SingleTree

> main = do   
>	  args <- getArgs
>	  me <- getProgName
>	  case args of
>          []     -> putStrLn("Usage: " ++ me ++ " [shortcut] <one or more consensus structures>")
>          ("test":consensi) -> putStrLn(gen_build_alg_cyk)
>          consensi -> putStrLn(
>                              "/*\n" ++ 
>                              "(multiple) covariance model for the following " ++ (show (length (a:as))) ++ " consensus structure(s):\n  " ++ 
>                              (foldr1 (\x -> \y -> x ++ "\n  " ++ y) (map (\(c, t) -> (fst c) ++ ":  " ++ (replicate (longestName - (length (fst c))) ' ')  ++ (snd c) ++ "    (as tree: '" ++ (show t) ++ "')") (zip (a:as) allTrees)) ) ++ "\n" ++
>                              "forming the unifying tree:\n" ++ 
>                              "  " ++ (show (commonTree)) ++ "\n" ++
>                              "*/\n\n" ++
>                              (get_build_code commonTree)
>                              )
>                  where commonTree = foldr (\x -> \y -> fuseTreelists y x) [] allTrees -- fuse all consensus trees into one large tree holding all consensus information with maximal overlaps between all single trees 
>                        allTrees = [ head (gra_g5skip alg_tree structure) | (name,structure) <- (a:as) ] 
>                        longestName = getLongestName (a:as)                            -- get the length of the longes consensus name for pretty printing consensus information as a comment in the generated code file
>                        (a:as) = map (\x -> (read x)::Consensus) (consensi)            -- parse input as a list of tupels of strings, i.e. a list of "Consensus"


> t1 = P 1 (P 2 (E 3) (E 4)) (O 5 (E 6))
> t2 = P 1 (P 2 (E 3) (E 4)) (E 5)

