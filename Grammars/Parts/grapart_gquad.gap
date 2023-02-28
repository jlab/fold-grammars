/* for layout definition see
 * - ViennaRNA 2.5.1 src/ViennaRNA/params/basic.h, lines 45ff
 *   #define   VRNA_GQUAD_MAX_STACK_SIZE     7
 *   #define   VRNA_GQUAD_MIN_STACK_SIZE     2
 *   #define   VRNA_GQUAD_MAX_LINKER_LENGTH  15
 *   #define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
 *   #define   VRNA_GQUAD_MIN_BOX_SIZE       ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + \
 *                                            (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
 *   #define   VRNA_GQUAD_MAX_BOX_SIZE       ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + \
 *                                            (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))
 * - ViennaRNA 2.5.1 src/ViennaRNA/gquad.c, lines 1848 ff
 */
  gquadruplex     = gquad(REGION with minsize(2) with maxsize(7) with onlychar(G_BASE),
                          REGION with minsize(1) with maxsize(15) with unpaired,
                          REGION with minsize(2) with maxsize(7) with onlychar(G_BASE),
                          REGION with minsize(1) with maxsize(15) with unpaired,
                          REGION with minsize(2) with maxsize(7) with onlychar(G_BASE),
                          REGION with minsize(1) with maxsize(15) with unpaired,
                          REGION with minsize(2) with maxsize(7) with onlychar(G_BASE)
                         ) with_overlay gquad_same_quarted_sizes # h;
