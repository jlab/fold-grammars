> module Datatypes where

> import Data.List

> -- data structure for a G5-Tree, consisting of _E_nd of structures, _O_pen aka unpaired bases and _P_aired bases
> data Tree a = E a | O a (Tree a) | P a (Tree a) (Tree a) deriving (Eq, Show)

> -- data structure for G5-Tree-Lists. This is neccessary if two or more G5-Tree shall be merged into one combined tree. Used for alternative Sub-Structures for RNA families.
> data TreeList a = El a | Ol a [TreeList a] | Pl a [TreeList a] [TreeList a] deriving (Eq, Show, Read)


> -- Function to fuse a potentially empty list of G5-Tree-Lists (ts) with one more G5-Tree (tree).
> -- If tree is already included in ts then ts can be returned unchanged, otherwise we simultaniously traverse ts and tree until they begin to differ in their indexed topology.
> -- ts is partitioned into those G5-Tree-Lists that start with the same indexed root node as tree (which by construction can one at most) = sametype and all other nodes. The first case is the lucky one, because that does not increase the G5-Tree-List at a close to the root node.
> fuseTreelists :: (Eq a) => [TreeList a] -> Tree a -> [TreeList a]
> fuseTreelists ts tree = if (elemTree ts tree) then ts else difftype ++ rest
>          where rest   = case (sametype,      tree   ) of
>                              ([],            tree   ) -> tree2treelist tree
>                              ([El i'],       E _    ) -> [El i']
>                              ([Ol i' x'],    O _ x  ) -> [Ol i' (fuseTreelists x' x)]
>                              ([Pl i' x' y'], P _ x y) -> [Pl i' (fuseTreelists x' x) (fuseTreelists y' y)]
>                (sametype, difftype) = partition (\x -> isIDroot x tree) ts


> -- Checks if a G5-Tree-Lists already containes the given G5-Tree.
> -- Recursion is two-fold. One direction is down the first tree of the list if the type and index of the root node is identical to the tree
> elemTree :: (Eq a) => [TreeList a] -> Tree a -> Bool
> elemTree [] _         = False
> elemTree (t:ts) tree  = ((isIDroot t tree) && (haveSameChildren)) || elemTree ts tree 
>           where haveSameChildren = case (t,          tree   ) of
>                                         (El _,       E _    ) -> True
>                                         (Ol _ xs,    O _ x  ) -> elemTree xs x
>                                         (Pl _ xs ys, P _ x y) -> elemTree xs x && elemTree ys y


> -- converts a Tree into a TreeList
> tree2treelist :: Tree a -> [TreeList a]
> tree2treelist (E i)     = [El i]
> tree2treelist (O i x)   = [Ol i (tree2treelist x)]
> tree2treelist (P i x y) = [Pl i (tree2treelist x) (tree2treelist y)]


> -- returns true if the root Node type of a Tree and a TreeList as well as their labels are equal
> isIDroot :: (Eq a) => TreeList a -> Tree a -> Bool
> isIDroot (El i')     (E i)     = (i == i')
> isIDroot (Ol i' _)   (O i _)   = (i == i')
> isIDroot (Pl i' _ _) (P i _ _) = (i == i')
> isIDroot _           _         = False










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






> type Consensus = (String, String)
> getLongestName :: [Consensus] -> Int
> getLongestName (a:as) = maximum (map (\x -> length(fst x)) (a:as))
