from pylib.gapc import *

#import pylib.gapc
#// A dynamic programming evaluator generated by GAP-C.
#//
#//   GAP-C version:
#//     2022.02.01
#//
#//   GAP-C call:
#//     gapc nodangle.gap -p alg_pfunc
#//
#//


#define GAPC_MOD_TRANSLATION_UNIT
#include "out.hh"

#include <rtlib/subopt.hh>
#include "rna.hh"
#include "Extensions/singlefold.hh"
#include "Extensions/mfesubopt.hh"
#include "Extensions/probabilities.hh"
#include "Extensions/shapes.hh"
#include "Extensions/mea.hh"
#include "Extensions/probing.hh"

#include <rtlib/generic_opts.hh>
#include "rtlib/pareto_dom_sort.hh"
#include "rtlib/pareto_yukish_ref.hh"

global weak_table
global strong_table
global iloop_table
global ml_comps_table
global ml_comps1_table
global dangle_table
global struct_table

global t_0_seq

def init(inputsequence):
    gapcrna.librna_read_param_file(None)

    global t_0_seq
    t_0_seq = inputsequence.upper().replace('A','\1').replace('C','\2').replace('G','\3').replace('U','\4')

    global weak_table
    weak_table = DPtable(len(t_0_seq))

    global strong_table
    strong_table = DPtable(len(t_0_seq))

    global iloop_table
    iloop_table = DPtable(len(t_0_seq))

    global ml_comps_table
    ml_comps_table = DPtable(len(t_0_seq))

    global ml_comps1_table
    ml_comps1_table = DPtable(len(t_0_seq))

    global dangle_table
    dangle_table = DPtable(len(t_0_seq))

    global struct_table
    struct_table = DPtable(len(t_0_seq))

def nt_dangle(t_0_i:int, t_0_j:int) -> float:
    if (dangle_table.is_tabulated(t_0_i, t_0_j)):
#     {
        return dangle_table.get(t_0_i, t_0_j)
#     }
#
    answers = []
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 5)):
#   {
        ret_3 = LOC(t_0_seq, t_0_j, t_0_j)
#     TUSubsequence a_2 = ret_3;
        if (is_not_empty(ret_3)):
#     {
            ret_1 = LOC(t_0_seq, t_0_i, t_0_i)
#       TUSubsequence a_0 = ret_1;
            if (is_not_empty(ret_1)):
#       {
                ret_2 = nt_strong(t_0_i, t_0_j)
#         double a_1 = ret_2;
                if (is_not_empty(ret_2)):
#           {
                    ret_0 = drem(ret_1, ret_2, ret_3)
#           }
#
#         else
#           {
#             empty( ret_0);
#           }
#
#         erase( a_1);
#       }
#
#       else
#         {
#           empty( ret_0);
#         }
#
#       erase( a_0);
#     }
#
#     else
#       {
#         empty( ret_0);
#       }
#
#     erase( a_2);
#   }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
#     {
        answers.append(ret_0)
#       push_back_sum( answers, ret_0);
#     }
#
    eval = h(answers)
#   erase( answers);
    dangle_table.set( t_0_i, t_0_j, eval)
    return dangle_table.get(t_0_i, t_0_j)
# }

def nt_hairpin(t_0_i:int, t_0_j:int) -> float:
#{
  if (((t_0_j - t_0_i) < 5)):
#     {
    return float_zero
#     }
#
#   double answers;
#   empty(answers);
#   empty( answers);
  answers = []
#   double ret_0;
  ret_0 = np.nan
  if (((t_0_j - t_0_i) >= 5)):
#     {
    if (basepair(t_0_seq, t_0_i, t_0_j)):
#       {
      ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#         TUSubsequence a_2 = ret_3;
      if (is_not_empty(ret_3)):
#         {
        ret_2 = np.nan;
        if ((minsize(t_0_seq, (t_0_i + 1), (t_0_j - 1), 3) and unpaired(t_0_seq, (t_0_i + 1), (t_0_j - 1)))):
#             {
          ret_2 = REGION(t_0_seq, (t_0_i + 1), (t_0_j - 1))
#             }
#
#           else
#             {
#               empty( ret_2);
#             }
#
#           TUSubsequence a_1 = ret_2;
          if (is_not_empty(ret_2)):
#           {
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#             TUSubsequence a_0 = ret_1;
            if (is_not_empty(ret_1)):
#               {
              ret_0 = hl(ret_1, ret_2, ret_3)
#               }
#
#             else
#               {
#                 empty( ret_0);
#               }
#
#             erase( a_0);
#           }
#
#           else
#             {
#               empty( ret_0);
#             }
#
#           erase( a_1);
#         }
#
#         else
#           {
#             empty( ret_0);
#           }
#
#         erase( a_2);
#       }
#
#       else
#       {
#         empty( ret_0);
#         empty( ret_0);
#       }
#
#     }
#
#   else
#     {
#       empty( ret_0);
#     }
#
  if (is_not_empty(ret_0)):
#     {
#       push_back_sum( answers, ret_0);
    answers.append(ret_0)
#     }
#
#   double eval = h(answers);
#   erase( answers);
  return h(answers)
# }
#
def nt_iloop(t_0_i:int, t_0_j:int) -> float:
# double &  out::nt_iloop(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (iloop_table.is_tabulated(t_0_i, t_0_j)):
#     {
        return iloop_table.get(t_0_i, t_0_j)
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
#
    if (((t_0_j - t_0_i) >= 9)):
#     {
        if (basepair(t_0_seq, t_0_i, t_0_j)):
#         {
            t_0_k_0 = (t_0_i + 2)
            while ((t_0_k_0 <= (t_0_j - 7)) and (t_0_k_0 <= (t_0_i + 31))):
#           for(          unsigned int t_0_k_0 = (t_0_i + 2); ((t_0_k_0 <= (t_0_j - 7)) && (t_0_k_0 <= (t_0_i + 31))); ++t_0_k_0)
#           {
                t_0_k_1 = ((t_0_j - 31)) if (((t_0_j - (t_0_k_0 + 5)) >= 31)) else ((t_0_k_0 + 5))
                while (t_0_k_1 <= (t_0_j - 2)):
#             for(            unsigned int t_0_k_1 = (((t_0_j - (t_0_k_0 + 5)) >= 31)) ? ((t_0_j - 31)) : ((t_0_k_0 + 5)); (t_0_k_1 <= (t_0_j - 2)); ++t_0_k_1)
#             {
                    ret_5 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#               TUSubsequence a_4 = ret_5;
                    if (is_not_empty(ret_5)):
#               {
#                 TUSubsequence ret_4;
                        ret_4 = np.nan;
                        if ((maxsize(t_0_seq, t_0_k_1, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_1, (t_0_j - 1)))):
#                   {
                            ret_4 = REGION(t_0_seq, t_0_k_1, (t_0_j - 1))
#                   }
#
#                 else
#                   {
#                     empty( ret_4);
#                   }
#
#                 TUSubsequence a_3 = ret_4;
                        if (is_not_empty(ret_4)):
#                 {
                            ret_2 = np.nan
                            if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
#                     {
                                ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
#                     }
#
#                   else
#                     {
#                       empty( ret_2);
#                     }
#
#                   TUSubsequence a_1 = ret_2;
                            if (is_not_empty(ret_2)):
#                   {
                                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#                     TUSubsequence a_0 = ret_1;
                                if (is_not_empty(ret_1)):
#                     {
                                    ret_3 = nt_strong(t_0_k_0, t_0_k_1)
#                       double a_2 = ret_3;
                                    if (is_not_empty(ret_3)):
#                       {
                                        answers.append(il(ret_1, ret_2, ret_3, ret_4, ret_5))
#                         push_back_sum( answers, ans);
#                       }
#
#                       erase( a_2);
#                     }
#
#                     erase( a_0);
#                   }
#
#                   erase( a_1);
#                 }
#
#                 erase( a_3);
#               }
#
#               erase( a_4);
#             }
                    t_0_k_1 += 1
#
#           }
                t_0_k_0 += 1
#
#         }
#
#       else
#         {
#           empty( answers);
#         }
#
#     }
#
    eval = h(answers)
#   erase( answers);
    iloop_table.set( t_0_i, t_0_j, eval)
    return iloop_table.get(t_0_i, t_0_j)
# }
#
def nt_leftB(t_0_i:int, t_0_j:int) -> float:
# double out::nt_leftB(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (((t_0_j - t_0_i) < 8)):
#     {
       return float_zero;
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
#
    if (((t_0_j - t_0_i) >= 8)):
#     {
        if (basepair(t_0_seq, t_0_i, t_0_j)):
#         {
            t_0_k_0 = t_0_i + 2
            while ((t_0_k_0 <= (t_0_j - 6)) and (t_0_k_0 <= (t_0_i + 31))):
#           for(          unsigned int t_0_k_0 = (t_0_i + 2); ((t_0_k_0 <= (t_0_j - 6)) && (t_0_k_0 <= (t_0_i + 31))); ++t_0_k_0)
#           {
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#             TUSubsequence a_3 = ret_4;
                if (is_not_empty(ret_4)):
#             {
                    ret_2 = np.nan
                    if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
#                 {
                        ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
#                 }
#
#               else
#                 {
#                   empty( ret_2);
#                 }
#
#               TUSubsequence a_1 = ret_2;
                    if (is_not_empty(ret_2)):
#               {
                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#                 TUSubsequence a_0 = ret_1;
                        if (is_not_empty(ret_1)):
#                 {
                            ret_3 = nt_strong(t_0_k_0, (t_0_j - 1))
#                   double a_2 = ret_3;
                            if (is_not_empty(ret_3)):
#                   {
                                ans = bl(ret_1, ret_2, ret_3, ret_4)
                                answers.append(ans)
                t_0_k_0 += 1
#                   }
#
#                   erase( a_2);
#                 }
#
#                 erase( a_0);
#               }
#
#               erase( a_1);
#             }
#
#             erase( a_3);
#           }
#
#         }
#
#       else
#         {
#           empty( answers);
#         }
#
#     }
#
#   double eval = h(answers);
#   erase( answers);
    return h(answers)
# }
#
def nt_ml_comps(t_0_i:int, t_0_j:int) -> float:
# double &  out::nt_ml_comps(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (ml_comps_table.is_tabulated(t_0_i, t_0_j)):
#     {
        return ml_comps_table.get(t_0_i, t_0_j)
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 11)):
#   {
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
#       {
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#       }
#
#     else
#       {
#         empty( ret_1);
#       }
#
#     TUSubsequence a_0 = ret_1;
        if (is_not_empty(ret_1)):
#     {
            ret_2 = nt_ml_comps((t_0_i + 1), t_0_j)
#       double a_1 = ret_2;
            if (is_not_empty(ret_2)):
#         {
                ret_0 = sadd(ret_1, ret_2)
#         }
#
#       else
#         {
#           empty( ret_0);
#         }
#
#       erase( a_1);
#     }
#
#     else
#       {
#         empty( ret_0);
#       }
#
#     erase( a_0);
#   }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
#     {
        answers.append(ret_0)
#     }
#
#
    if (((t_0_j - t_0_i) >= 10)):
#     {
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= (t_0_j - 5)):
#       for(      unsigned int t_0_k_0 = (t_0_i + 5); (t_0_k_0 <= (t_0_j - 5)); ++t_0_k_0)
#       {
            ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
#         double a_4 = ret_6;
            if (is_not_empty(ret_6)):
#         {
                ret_4 = np.nan
                if (((t_0_k_0 - t_0_i) >= 5)):
#           {
                    ret_5 = nt_dangle(t_0_i, t_0_k_0)
#             double a_3 = ret_5;
                    if (is_not_empty(ret_5)):
#               {
                        ret_4 = incl(ret_5)
#               }
#
#             else
#               {
#                 empty( ret_4);
#               }
#
#             erase( a_3);
#           }
#
#           else
#             {
#               empty( ret_4);
#             }
#
#           double a_2 = ret_4;
                if (is_not_empty(ret_4)):
#           {
                    answers.append(cadd(ret_4, ret_6))
#             push_back_sum( answers, ans);
#           }
#
#           erase( a_2);
#         }
#
#         erase( a_4);
#       }
            t_0_k_0 += 1
#
#     }
#
    eval = h(answers)
#   erase( answers);
    ml_comps_table.set( t_0_i, t_0_j, eval)
    return ml_comps_table.get(t_0_i, t_0_j)
# }
#
def nt_ml_comps1(t_0_i:int, t_0_j:int) -> float:
# double &  out::nt_ml_comps1(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (ml_comps1_table.is_tabulated(t_0_i, t_0_j)):
#     {
        return ml_comps1_table.get(t_0_i, t_0_j)
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 6)):
#   {
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
#       {
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#       }
#
#     else
#       {
#         empty( ret_1);
#       }
#
#     TUSubsequence a_0 = ret_1;
        if (is_not_empty(ret_1)):
#     {
            ret_2 = nt_ml_comps1((t_0_i + 1), t_0_j)
#       double a_1 = ret_2;
            if (is_not_empty(ret_2)):
#         {
                ret_0 = sadd(ret_1, ret_2)
#         }
#
#       else
#         {
#           empty( ret_0);
#         }
#
#       erase( a_1);
#     }
#
#     else
#       {
#         empty( ret_0);
#       }
#
#     erase( a_0);
#   }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
#     {
        answers.append(ret_0)
#     }
#
#
    if (((t_0_j - t_0_i) >= 10)):
#     {
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= (t_0_j - 5)):
#       for(      unsigned int t_0_k_0 = (t_0_i + 5); (t_0_k_0 <= (t_0_j - 5)); ++t_0_k_0)
#       {
            ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
#         double a_4 = ret_6;
            if (is_not_empty(ret_6)):
#         {
                ret_4 = np.nan
                if (((t_0_k_0 - t_0_i) >= 5)):
#           {
                    ret_5 = nt_dangle(t_0_i, t_0_k_0)
#             double a_3 = ret_5;
                    if (is_not_empty(ret_5)):
#               {
                        ret_4 = incl(ret_5)
#               }
#
#             else
#               {
#                 empty( ret_4);
#               }
#
#             erase( a_3);
#           }
#
#           else
#             {
#               empty( ret_4);
#             }
#
#           double a_2 = ret_4;
                if (is_not_empty(ret_4)):
#           {
#             double ans = cadd(a_2, a_4);
                    answers.append(cadd(ret_4, ret_6))
#           }
#
#           erase( a_2);
#         }
#
#         erase( a_4);
#       }
            t_0_k_0 += 1
#
#     }
#
    ret_7 = np.nan
    if (((t_0_j - t_0_i) >= 5)):
#   {
        ret_8 = nt_dangle(t_0_i, t_0_j)
#     double a_5 = ret_8;
        if (is_not_empty(ret_8)):
#       {
            ret_7 = incl(ret_8)
#       }
#
#     else
#       {
#         empty( ret_7);
#       }
#
#     erase( a_5);
#   }
#
#   else
#     {
#       empty( ret_7);
#     }
#
    if (is_not_empty(ret_7)):
#     {
        answers.append(ret_7)
#     }
#
#
    if (((t_0_j - t_0_i) >= 6)):
#     {
        t_0_k_1 = (t_0_i + 5)
        while (t_0_k_1 <= (t_0_j - 1)):
#        for(      unsigned int t_0_k_1 = (t_0_i + 5); (t_0_k_1 <= (t_0_j - 1)); ++t_0_k_1)
#       {
            ret_12 = np.nan
            if (unpaired(t_0_seq, t_0_k_1, t_0_j)):
#           {
                ret_12 = REGION(t_0_seq, t_0_k_1, t_0_j)
#           }
#
#         else
#           {
#             empty( ret_12);
#           }
#
#         TUSubsequence a_8 = ret_12;
            if (is_not_empty(ret_12)):
#         {
                ret_10 = np.nan
                if (((t_0_k_1 - t_0_i) >= 5)):
#           {
                    ret_11 = nt_dangle(t_0_i, t_0_k_1)
#             double a_7 = ret_11;
                    if (is_not_empty(ret_11)):
#               {
                        ret_10 = incl(ret_11)
#               }
#
#             else
#               {
#                 empty( ret_10);
#               }
#
#             erase( a_7);
#           }
#
#           else
#             {
#               empty( ret_10);
#             }
#
#           double a_6 = ret_10;
                if (is_not_empty(ret_10)):
#           {
                    answers.append(addss(ret_10, ret_12))
#             push_back_sum( answers, ans);
#           }
#
#           erase( a_6);
#         }
#
#         erase( a_8);
#       }
            t_0_k_1 += 1
#
#     }
#
    eval = h(answers)
#   erase( answers);
    ml_comps1_table.set( t_0_i, t_0_j, eval)
    return ml_comps1_table.get(t_0_i, t_0_j)
# }
#
def nt_multiloop(t_0_i:int, t_0_j:int) -> float:
# double out::nt_multiloop(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (((t_0_j - t_0_i) < 12)):
#     {
        return float_zero;
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 12)):
#     {
        if (basepair(t_0_seq, t_0_i, t_0_j)):
#       {
            ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#         TUSubsequence a_2 = ret_3;
            if (is_not_empty(ret_3)):
#         {
                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#           TUSubsequence a_0 = ret_1;
                if (is_not_empty(ret_1)):
#           {
                    ret_2 = nt_ml_comps((t_0_i + 1), (t_0_j - 1))
#             double a_1 = ret_2;
                    if (is_not_empty(ret_2)):
#               {
                        ret_0 = ml(ret_1, ret_2, ret_3)
#               }
#
#             else
#               {
#                 empty( ret_0);
#               }
#
#             erase( a_1);
#           }
#
#           else
#             {
#               empty( ret_0);
#             }
#
#           erase( a_0);
#         }
#
#         else
#           {
#             empty( ret_0);
#           }
#
#         erase( a_2);
#       }
#
#       else
#       {
#         empty( ret_0);
#         empty( ret_0);
#       }
#
#     }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
#     {
        answers.append(ret_0)
#     }
#
#   double eval = h(answers);
#   erase( answers);
    return h(answers)
# }
#
def nt_rightB(t_0_i:int, t_0_j:int) -> float:
# double out::nt_rightB(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (((t_0_j - t_0_i) < 8)):
#     {
        return float_zero
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
#

    if (((t_0_j - t_0_i) >= 8)):
#     {
        if (basepair(t_0_seq, t_0_i, t_0_j)):
#         {
            t_0_k_0 = ((t_0_j - 31)) if (((t_0_j - (t_0_i + 6)) >= 31)) else ((t_0_i + 6))
            while (t_0_k_0 <= (t_0_j - 2)):
#           for(          unsigned int t_0_k_0 = (((t_0_j - (t_0_i + 6)) >= 31)) ? ((t_0_j - 31)) : ((t_0_i + 6)); (t_0_k_0 <= (t_0_j - 2)); ++t_0_k_0)
#           {
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#             TUSubsequence a_3 = ret_4;
                if (is_not_empty(ret_4)):
#             {
                    ret_3 = np.nan
                    if ((maxsize(t_0_seq, t_0_k_0, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_0, (t_0_j - 1)))):
#                 {
                        ret_3 = REGION(t_0_seq, t_0_k_0, (t_0_j - 1))
#                 }
#
#               else
#                 {
#                   empty( ret_3);
#                 }
#
#               TUSubsequence a_2 = ret_3;
                    if (is_not_empty(ret_3)):
#               {
                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#                 TUSubsequence a_0 = ret_1;
                        if (is_not_empty(ret_1)):
#                 {
                            ret_2 = nt_strong((t_0_i + 1), t_0_k_0)
#                   double a_1 = ret_2;
                            if (is_not_empty(ret_2)):
#                   {
                                answers.append(br(ret_1, ret_2, ret_3, ret_4))
                t_0_k_0 += 1

#                     push_back_sum( answers, ans);
#                   }
#
#                   erase( a_1);
#                 }
#
#                 erase( a_0);
#               }
#
#               erase( a_2);
#             }
#
#             erase( a_3);
#           }
#
#         }
#
#       else
#         {
#           empty( answers);
#         }
#
#     }
#
#   double eval = h(answers);
#   erase( answers);
#   return eval;
    return h(answers)
# }
#
def nt_stack(t_0_i:int, t_0_j:int) -> float:
# double out::nt_stack(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (((t_0_j - t_0_i) < 7)):
#     {
       return float_zero;
#     }
#
#   double answers;
    answers = []
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 7)):
#     {
       if (basepair(t_0_seq, t_0_i, t_0_j)):
#       {
         ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#         TUSubsequence a_2 = ret_3;
         if (is_not_empty(ret_3)):
#         {
           ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#           TUSubsequence a_0 = ret_1;
           if (is_not_empty(ret_1)):
#           {
             ret_2 = nt_weak((t_0_i + 1), (t_0_j - 1))
#             double a_1 = ret_2;
             if (is_not_empty(ret_2)):
#               {
                 ret_0 = sr(ret_1, ret_2, ret_3)
#               }
#
#             else
#               {
#                 empty( ret_0);
#               }
#
#             erase( a_1);
#           }
#
#           else
#             {
#               empty( ret_0);
#             }
#
#           erase( a_0);
#         }
#
#         else
#           {
#             empty( ret_0);
#           }
#
#         erase( a_2);
#       }
#
#       else
#       {
#         empty( ret_0);
#         empty( ret_0);
#       }
#
#     }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
#     {
#       push_back_sum( answers, ret_0);
        answers.append(ret_0)
#     }
#
#   double eval = h(answers);
#   erase( answers);
#   return eval;
    return h(answers)
# }
#
def nt_strong(t_0_i:int, t_0_j:int) -> float:
# double &  out::nt_strong(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (strong_table.is_tabulated(t_0_i, t_0_j)):
#     {
       return strong_table.get(t_0_i, t_0_j)
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
#   double ret_0;
#   empty( ret_0);
    ret_0 = np.nan
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, False)):
#   {
        ret_1 = np.nan
        if (((t_0_j - t_0_i) >= 7)):
#       {
            if (basepair(t_0_seq, t_0_i, t_0_j)):
#         {
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
#           TUSubsequence a_2 = ret_4;
                if (is_not_empty(ret_4)):
#           {
                    ret_2 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#             TUSubsequence a_0 = ret_2;
                    if (is_not_empty(ret_2)):
#             {
                        ret_3 = nt_weak((t_0_i + 1), (t_0_j - 1))
#               double a_1 = ret_3;
                        if (is_not_empty(ret_3)):
#                 {
                            ret_0 = sr(ret_2, ret_3, ret_4)
#                 }
#
#               else
#                 {
#                   empty( ret_1);
#                 }
#
#               erase( a_1);
#             }
#
#             else
#               {
#                 empty( ret_1);
#               }
#
#             erase( a_0);
#           }
#
#           else
#             {
#               empty( ret_1);
#             }
#
#           erase( a_2);
#         }
#
#         else
#         {
#           empty( ret_1);
#           empty( ret_1);
#         }
#
#       }
#
#     else
#       {
#         empty( ret_1);
#       }
#
#     ret_0 = ret_1;
#   }
#
    if (is_not_empty(ret_0)):
        answers.append(ret_0)
#     {
#       push_back_sum( answers, ret_0);
#     }
#
    ret_5 = np.nan
#   empty( ret_5);
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
#   {
        ret_5 = nt_weak(t_0_i, t_0_j)
#     ret_5 = ret_6;
#   }
#
    if (is_not_empty(ret_5)):
#     {
       answers.append(ret_5)
#     }
#
    eval = h(answers)
#   erase( answers);
    strong_table.set( t_0_i, t_0_j, eval)
    return strong_table.get(t_0_i, t_0_j)
# }
#
def nt_struct(t_0_i:int) -> float:
    t_0_right_most = len(t_0_seq)
# double &  out::nt_struct(unsigned int t_0_i)
# {
    if (struct_table.is_tabulated(t_0_i, 0)):
#     {
        return struct_table.get(t_0_i, 0)
#     }
#
    answers = []
#   double answers;
#   empty(answers);
#   empty( answers);
    ret_0 = np.nan
    if (((t_0_right_most - t_0_i) >= 1)):
#   {
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
#       {
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
#       }
#
#     else
#       {
#         empty( ret_1);
#       }
#
#     TUSubsequence a_0 = ret_1;
        if (is_not_empty(ret_1)):
#     {
            ret_2 = nt_struct((t_0_i + 1))
#       double a_1 = ret_2;
            if (is_not_empty(ret_2)):
#         {
                ret_0 = sadd(ret_1, ret_2)
#         }
#
#       else
#         {
#           empty( ret_0);
#         }
#
#       erase( a_1);
#     }
#
#     else
#       {
#         empty( ret_0);
#       }
#
#     erase( a_0);
#   }
#
#   else
#     {
#       empty( ret_0);
#     }
#
    if (is_not_empty(ret_0)):
        answers.append(ret_0)
#     {
#       push_back_sum( answers, ret_0);
#     }
#
#
    if (((t_0_right_most - t_0_i) >= 5)):
#     {
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= t_0_right_most):
#        for(      unsigned int t_0_k_0 = (t_0_i + 5); (t_0_k_0 <= t_0_right_most); ++t_0_k_0)
#       {
            ret_5 = nt_struct(t_0_k_0)
#         double a_3 = ret_5;
            if (is_not_empty(ret_5)):
#         {
                ret_4 = nt_dangle(t_0_i, t_0_k_0)
#           double a_2 = ret_4;
                if (is_not_empty(ret_4)):
#           {
                    answers.append(cadd(ret_4, ret_5))
#             push_back_sum( answers, ans);
#           }
#
#           erase( a_2);
#         }
#
#         erase( a_3);
#       }
            t_0_k_0 += 1
#
#     }
#
    ret_6 = np.nan
    if ((((t_0_right_most - t_0_i) >= 0) and ((t_0_right_most - t_0_i) <= 0))):
#   {
        ret_7 = LOC(t_0_seq, t_0_i, t_0_i)
#     TUSubsequence a_4 = ret_7;
        if (is_not_empty(ret_7)):
#       {
            ret_6 = nil(ret_7);
#       }
#
#     else
#       {
#         empty( ret_6);
#       }
#
#     erase( a_4);
#   }
#
#   else
#     {
#       empty( ret_6);
#     }
#
    if (is_not_empty(ret_6)):
#     {
        answers.append(ret_6)
#       push_back_sum( answers, ret_6);
#     }
#
    eval = h(answers)
#   erase( answers);
    struct_table.set( t_0_i, 0, eval)
    return struct_table.get(t_0_i, 0)
# }
#
def nt_weak(t_0_i:int, t_0_j:int) -> float:
    global weak_table
# double &  out::nt_weak(unsigned int t_0_i, unsigned int t_0_j)
# {
    if (weak_table.is_tabulated(t_0_i, t_0_j)):
#     {
       return weak_table.get(t_0_i, t_0_j)
#     }
#
#   double answers;
    answers = []
#   empty(answers);
#   empty( answers);
#
#
    ret_1 = nt_stack(t_0_i, t_0_j)
    if (is_not_empty(ret_1)):
        answers.append(ret_1)
#     {
#       push_back_sum( answers, ret_1);
#     }
#
    ret_2 = nt_hairpin(t_0_i, t_0_j)
    if (is_not_empty(ret_2)):
        answers.append(ret_2)
#     {
#       push_back_sum( answers, ret_2);
#     }
#
    ret_3 = nt_leftB(t_0_i, t_0_j)
    if (is_not_empty(ret_3)):
        answers.append(ret_3)
#     {
#       push_back_sum( answers, ret_3);
#     }
#
    ret_4 = nt_rightB(t_0_i, t_0_j)
    if (is_not_empty(ret_4)):
        answers.append(ret_4)
#     {
#       push_back_sum( answers, ret_4);
#     }
#
    ret_5 = nt_iloop(t_0_i, t_0_j)
    if (is_not_empty(ret_5)):
        answers.append(ret_5)
#     {
#       push_back_sum( answers, ret_5);
#     }
#
    ret_6 = nt_multiloop(t_0_i, t_0_j)
    if (is_not_empty(ret_6)):
        answers.append(ret_6)
#     {
#       push_back_sum( answers, ret_6);
#     }
#
    eval = h(answers)
#   erase( answers);
    weak_table.set( t_0_i, t_0_j, eval)
    return weak_table.get(t_0_i, t_0_j)
# }
#
#

if True:
    # pfunc algebra
    def addss(x:float, r:Basic_Subsequence):
        return ((scale((r.j - r.i)) * x) * mk_pf(ss_energy(r)))

    def bl(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return ((scale(((2 + lr.j) - lr.i)) * x) * mk_pf(bl_energy(lr, rb)))

    def br(lb:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return ((scale(((2 + rr.j) - rr.i)) * x) * mk_pf(br_energy(lb, rr)))

    def cadd(x:float, y:float):
        return x*y

    def drem(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return (x * mk_pf(termau_energy(lb, rb)))

    def hl(lb:Basic_Subsequence, r:Basic_Subsequence, rb:Basic_Subsequence):
        return (scale(((2 + r.j) - r.i)) * mk_pf(hl_energy(r)))

    def il(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return ((scale(((((2 + lr.j) - lr.i) + rr.j) - rr.i)) * x) * mk_pf(il_energy(lr, rr)))

    def incl(x:float):
        return (x * mk_pf(ul_energy()))

    def ml(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return ((scale(2) * x) * mk_pf(((ml_energy() + ul_energy()) + termau_energy(lb, rb))))

    def nil(n:Basic_Subsequence):
        return 1

    def sadd(lb:Basic_Subsequence, x:float):
        return ((scale(1) * x) * mk_pf(sbase_energy()))

    def sr(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
       return ((scale(2) * x) * mk_pf(sr_energy(lb, rb)))

    def h(i:[float]):
        if len(i) > 0:
            return np.sum(i)
        else:
            return np.nan

else:
    # mfe algebra
    def addss(x:float, r:Basic_Subsequence):
        return x + ss_energy(r)

    def bl(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return x + bl_energy(lr, rb)

    def br(lb:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return x + br_energy(lb, rr)

    def cadd(x:float, y:float):
        return x + y

    def drem(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return x + termau_energy(lb, rb)

    def hl(lb:Basic_Subsequence, r:Basic_Subsequence, rb:Basic_Subsequence):
        return hl_energy(r)

    def il(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return x + il_energy(lr, rr)

    def incl(x:float):
        return x + ul_energy()

    def ml(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return x + ml_energy() + ul_energy() + termau_energy(lb, rb)

    def nil(n:Basic_Subsequence):
        return 0

    def sadd(lb:Basic_Subsequence, x:float):
        return x + sbase_energy()

    def sr(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
       return x + sr_energy(lb, rb)

    def h(i:[float]) -> [float]:
        if len(i) > 0:
            return np.min(i)
        else:
            return np.nan
