> module Datatypes where

> import Data.List

> data Tree a = E a | O a (Tree a) | P a (Tree a) (Tree a) deriving (Eq, Show)

> -- concatenate two trees, neccessary if a pair is deleted
> (+:+) :: Tree a -> Tree a -> Tree a
> (+:+) a b = case (a,b) of
>                  (E i    , _        ) -> b
>                  (O i x  , _        ) -> O i (x +:+ b)
>                  (P i x y, _        ) -> P i x (y +:+ b)

> getIndex :: Tree a -> a
> getIndex (E i) = i
> getIndex (O i _) = i
> getIndex (P i _ _) = i

> getEnd :: Tree a -> Tree a
> getEnd (E i) = E i
> getEnd (O _ x) = getEnd x
> getEnd (P _ x y) = getEnd y 



> data TreeList a = El a | Ol a [TreeList a] | Pl a [TreeList a] [TreeList a] deriving (Eq, Show, Read)

> -- applies depth first numbers to all nodes, i.e. nodes from the tree structure and alternatives
> indexTreeList :: [TreeList a] -> Int -> [TreeList (a,Int)]
> indexTreeList [El a] i = [El (a, i)]
> indexTreeList [Ol a xs] i = [Ol (a,i) (indexTreeList xs (i+1))]
> indexTreeList [Pl a xs ys] i = [Pl (a,i) left right]
>   where left = (indexTreeList xs (i+1))
>         right = (indexTreeList ys (i+1+(getTreeListSize left)))
> indexTreeList (a:as) i = first ++ (indexTreeList as (i+(getTreeListSize first)))
>   where first = indexTreeList [a] i

> -- gets the size of a TreeList, i.e. number of all nodes may they be from the tree structure or alternatives
> getTreeListSize :: [TreeList a] -> Int
> getTreeListSize [El _] = 1
> getTreeListSize [Ol _ xs] = 1 + getTreeListSize xs
> getTreeListSize [Pl _ xs ys] = 1 + (getTreeListSize xs) + (getTreeListSize ys)
> getTreeListSize (a:as) = (getTreeListSize [a]) + (getTreeListSize as)



> getTreeListIndex :: TreeList a -> a
> getTreeListIndex (El a) = a
> getTreeListIndex (Ol a _) = a
> getTreeListIndex (Pl a _ _) = a

> getTreeListLastIndex :: TreeList a -> a
> getTreeListLastIndex (El a) = a
> getTreeListLastIndex (Ol _ (x:xs)) = getTreeListLastIndex x
> getTreeListLastIndex (Pl _ _ (y:ys)) = getTreeListLastIndex y

> -- converts a Tree into a TreeList
> tree2treelist :: Tree a -> [TreeList a]
> tree2treelist (E i)     = [El i]
> tree2treelist (O i x)   = [Ol i (tree2treelist x)]
> tree2treelist (P i x y) = [Pl i (tree2treelist x) (tree2treelist y)]

> -- checks if a Tree is contained inside of a TreeList
> elemTree :: (Eq a) => Tree a -> [TreeList a] -> Bool
> elemTree _         []               = False
> elemTree (E i)     (El i':ts)       = (i == i') || elemTree (E i) ts
> elemTree (E i)     (_    :ts)       = elemTree (E i) ts
> elemTree (O i x)   (Ol i' x':ts)    = ((i == i') && (elemTree x x')) || elemTree (O i x) ts
> elemTree (O i x)   (_       :ts)    = elemTree (O i x) ts
> elemTree (P i x y) (Pl i' x' y':ts) = ((i == i') && (elemTree x x') && (elemTree y y')) || elemTree (P i x y) ts
> elemTree (P i x y) (_          :ts) = elemTree (P i x y) ts

> -- returns true if a root Node of a TreeList is of type a and has index i
> isNodeNr :: (Eq a) => Char -> a -> TreeList a -> Bool
> isNodeNr 'E' i (El a) = (i == a)
> isNodeNr 'O' i (Ol a _) = (i == a)
> isNodeNr 'P' i (Pl a _ _) = (i == a)
> isNodeNr _ _ _ = False

> fuseTreelists :: (Eq a) => [TreeList a] -> Tree a -> [TreeList a]
> fuseTreelists ts (E i)     = if (elemTree (E i) ts) then ts else ts ++ [El i]
> fuseTreelists [] (O i x)   = tree2treelist (O i x)
> fuseTreelists ts (O i x)   = if (elemTree (O i x) ts) then ts else nonos ++ rest
>          where rest        = case os of
>                                   [] -> fuseTreelists [] (O i x)
>                                   [Ol i' x'] -> [Ol i (fuseTreelists x' x)]
>                (os, nonos) = partition (\x -> isNodeNr 'O' i x) ts
> fuseTreelists [] (P i x y) = tree2treelist (P i x y)
> fuseTreelists ts (P i x y) = if (elemTree (P i x y) ts) then ts else nonps ++ rest
>          where rest        = case ps of
>                                   [] -> fuseTreelists [] (P i x y)
>                                   [Pl i' x' y'] -> [Pl i (fuseTreelists x' x) (fuseTreelists y' y)]
>                (ps, nonps) = partition (\x -> isNodeNr 'P' i x) ts


> type Consensus = (String, String)
> getLongestName :: [Consensus] -> Int
> getLongestName (a:as) = maximum (map (\x -> length(fst x)) (a:as))
