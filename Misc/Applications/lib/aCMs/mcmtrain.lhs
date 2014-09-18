> module Main where

> import System.Environment(getArgs,getProgName)

Haskell header:

> import Datatypes
> import Treelist2cm
> import ADPgenerator
> import SingleTree

> main = do   
>	  args <- getArgs
>	  me <- getProgName
>	  case args of
>          []     -> putStrLn("Usage: " ++ me ++ " [shortcut] <one or more consensus structures>")
>          ("shortcut":ft) -> putStrLn(
>                              "/*\n" ++ 
>                              "(multiple) covariance model for an unknown number of consensus structure(s):\n" ++ 
>                              "given is just the unifying tree:\n" ++
>                              "  " ++ (show t) ++ "\n" ++
>                              "*/\n\n" ++
>                               get_train_code t
>                             )
>                  where t :: [TreeList Int]
>                        t = read (head ft)
>          ("trees":(consensi)) -> putStrLn(
>                                        foldr1 (\x -> \y -> x ++ "\n" ++ y) (map (\(name, tree) -> name ++ ":  " ++ (replicate (longestName - (length name)) ' ') ++ (show tree)) (zip (map fst (a:as)) allTrees))
>                                      )
>                  where allTrees = map snd [head (grammar gen (a, b)) | b <- (a:as)]
>                        longestName = getLongestName (a:as)
>                        (a:as) = map (\x -> (read x)::Consensus) (consensi)
>          consensi -> if (warning == "") then putStrLn(
>                              "/*\n" ++ 
>                              "(multiple) covariance model for the following " ++ (show (length (a:as))) ++ " consensus structure(s):\n  " ++ 
>                              (foldr1 (\x -> \y -> x ++ "\n  " ++ y) (map (\(c, t) -> (fst c) ++ ":  " ++ (replicate (longestName - (length (fst c))) ' ')  ++ (snd c) ++ "    (as tree: '" ++ (show t) ++ "')") (zip (a:as) allTrees)) ) ++ "\n" ++
>                              "forming the unifying tree:\n" ++ 
>                              "  " ++ (show (commonTree)) ++ "\n" ++
>                              "*/\n\n" ++
>                              (get_train_code commonTree)
>                            )
>                       else putStrLn(warning)
>                  where commonTree = foldr (\x -> \y -> fuseTreelists y x) [] allTrees -- fuse all consensus trees into one large tree holding all consensus information with maximal overlaps between all single trees 
> --                       allTrees = map snd [head (grammar gen (a, b)) | b <- (a:as)]   -- each given consensus is forced to form a pair with the first given consensus (yes, also with itsself) which is then transformed into a pair of aligned trees via ADPgenerator. Thus, it is ensured that both trees use the same alignment positions. The first consensus is just for those positions, the second is the consensus tree of interest.
>                        allTrees = [ head (gra_g5skip alg_tree structure) | (name,structure) <- (a:as) ] 
>                        warning = if (sum (foldr (++) [0] [(grammar ADPgenerator.count (a,b)) | b <- (a:as)]) < length (a:as)) then "your alignment contains incompatible consensus structrues!" else ""
>                        longestName = getLongestName (a:as)                            -- get the length of the longes consensus name for pretty printing consensus information as a comment in the generated code file
>                        (a:as) = map (\x -> (read x)::Consensus) (consensi)            -- parse input as a list of tupels of strings, i.e. a list of "Consensus"
						
						



> t1 = P 1 (P 2 (E 3) (E 4)) (O 5 (E 6))
> t2 = P 1 (P 2 (E 3) (E 4)) (E 5)

