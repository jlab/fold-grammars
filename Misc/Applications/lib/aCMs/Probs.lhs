> module Probs where

> import Data.Map (Map)
> import qualified Data.Map as Map

learned transition and emission probabilities for 'test'.

> getTransition :: String -> String
> getTransition nt = case (Map.lookup nt transitions) of
>                          Nothing -> "0"
>                          Just a -> a

> getEmission :: String -> String -> String
> getEmission nt sym = case (Map.lookup nt emissions) of
>                            Nothing -> "0"
>                            Just a -> case (Map.lookup sym a) of
>                                            Nothing -> "0"
>                                            Just a -> a

> getSeqCons :: String -> String
> getSeqCons nt = case (Map.lookup nt seqCons) of
>                       Nothing -> ".."
>                       Just a -> a

> transitions :: Map.Map String String
> transitions = Map.fromList [
>    ("MAT_1", "-0.321928094887362")
>  ]

> emissions :: Map.Map String (Map.Map String String)
> emissions = Map.fromList [
>    ("MAT_1", Map.fromList [
>                 ("a", "-1.73696559416621"),
>                 ("c", "-1.73696559416621"),
>                 ("g", "-1.73696559416621"),
>                 ("u", "-3.32192809488736")
>              ]
>    )
>  ]
> seqCons :: Map.Map String String
> seqCons = Map.fromList [
>    ("MAT_1", "..")
>  ]
