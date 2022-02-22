from pylib.gapc import *

design_original = ["dangle","iloop","ml_comps","ml_comps1","strong","struct","weak"]
design_all = ["hairpin","leftB","multiloop","rightB","stack","dangle","iloop","ml_comps","ml_comps1","strong","struct","weak"]

def init(inputsequence, algebra='pfunc', printstack=False, printBTstack=False, tabulateNTs=design_all):
    gapcrna.librna_read_param_file(None)
    global PRINTSTACK
    PRINTSTACK = printstack
    global PRINTBTSTACK
    PRINTBTSTACK = printBTstack
    global INDENT
    INDENT = ""
    global ALGEBRA
    ALGEBRA = algebra
    global COMPUTATIONALSTEPS
    COMPUTATIONALSTEPS = 0
    global STORAGE
    STORAGE = 0

    global t_0_seq
    t_0_seq = inputsequence.upper().replace('A','\1').replace('C','\2').replace('G','\3').replace('U','\4')

    global tables
    tables = dict()
    for nt in tabulateNTs:
        tables[nt] = DPtable(len(t_0_seq), nt)

    if 'struct' in tables:
        base = np.nan
        if algebra in ['pfunc', 'count']:
            base = 1.0
        elif algebra == 'mfe':
            base = 0.0
            #tables['struct'].bt_set(0,0,0.0)
        tables['struct'].bt_set(0,0,base)
        tables['struct'].bt_set_v2(0,0,base)
        #tables['struct'].bt_set(len(t_0_seq),0,base)

def incr():
    global INDENT
    INDENT = INDENT + " "
def decr():
    global INDENT
    INDENT = INDENT[:-1]
def computed():
    global COMPUTATIONALSTEPS
    COMPUTATIONALSTEPS += 1
def stored():
    global STORAGE
    STORAGE += 1

def nt_dangle(t_0_i:int, t_0_j:int, name="dangle", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # drem(LOC, strong, LOC)
    if (((t_0_j - t_0_i) >= 5)):
        ret_3 = LOC(t_0_seq, t_0_j, t_0_j)
        if (is_not_empty(ret_3)):
            ret_1 = LOC(t_0_seq, t_0_i, t_0_i)
            if (is_not_empty(ret_1)):
                ret_2 = nt_strong(t_0_i, t_0_j)
                if (is_not_empty(ret_2)):
                    ret_0 = drem(ret_1, ret_2, ret_3)
                    if (is_not_empty(ret_0)):
                        answers.append(ret_0)
                        add_trace(tables, 'strong', t_0_i,t_0_j, 'dangle', t_0_i, t_0_j, algfct=drem, algparams=[ret_1, 'x', ret_3], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_hairpin(t_0_i:int, t_0_j:int, name="hairpin", bwdpass=False) -> float:
    computed()
    if (((t_0_j - t_0_i) < 5)):
        return float_zero

    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall_nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # hl(BASE, REGION, BASE)
    if (((t_0_j - t_0_i) >= 5)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
            if (is_not_empty(ret_3)):
                ret_2 = np.nan;
                if ((minsize(t_0_seq, (t_0_i + 1), (t_0_j - 1), 3) and unpaired(t_0_seq, (t_0_i + 1), (t_0_j - 1)))):
                    ret_2 = REGION(t_0_seq, (t_0_i + 1), (t_0_j - 1))

                    if (is_not_empty(ret_2)):
                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                        if (is_not_empty(ret_1)):
                            ret_0 = hl(ret_1, ret_2, ret_3)
                            if (is_not_empty(ret_0)):
                                answers.append(ret_0)
                                add_trace(tables, 'hairpin', None, None, None, t_0_i, t_0_j, algfct=hl, algparams=[ret_1, ret_2, ret_3], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_iloop(t_0_i:int, t_0_j:int, name="iloop", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # il(BASE, REGION, strong, REGION, BASE)
    if (((t_0_j - t_0_i) >= 9)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = (t_0_i + 2)
            while ((t_0_k_0 <= (t_0_j - 7)) and (t_0_k_0 <= (t_0_i + 31))):
                t_0_k_1 = ((t_0_j - 31)) if (((t_0_j - (t_0_k_0 + 5)) >= 31)) else ((t_0_k_0 + 5))
                while (t_0_k_1 <= (t_0_j - 2)):
                    ret_5 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                    if (is_not_empty(ret_5)):
                        if ((maxsize(t_0_seq, t_0_k_1, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_1, (t_0_j - 1)))):
                            ret_4 = REGION(t_0_seq, t_0_k_1, (t_0_j - 1))
                            if (is_not_empty(ret_4)):
                                if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
                                    ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
                                    if (is_not_empty(ret_2)):
                                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                                        if (is_not_empty(ret_1)):
                                            ret_3 = nt_strong(t_0_k_0, t_0_k_1)
                                            if (is_not_empty(ret_3)):
                                                res = il(ret_1, ret_2, ret_3, ret_4, ret_5)
                                                answers.append(res)
                                                add_trace(tables, 'strong', t_0_i,t_0_j, 'iloop', t_0_k_0, t_0_k_1, algfct=il, algparams=[ret_1, ret_2, 'x', ret_4, ret_5], bwdpass=bwdpass)
                    t_0_k_1 += 1
                t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_leftB(t_0_i:int, t_0_j:int, name="leftB", bwdpass=False) -> float:
    computed()
    if (((t_0_j - t_0_i) < 8)):
       return float_zero;

    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # bl(BASE, REGION, strong, BASE)
    if (((t_0_j - t_0_i) >= 8)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = t_0_i + 2
            while ((t_0_k_0 <= (t_0_j - 6)) and (t_0_k_0 <= (t_0_i + 31))):
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                if (is_not_empty(ret_4)):
                    if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
                        ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
                        if (is_not_empty(ret_2)):
                            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                            if (is_not_empty(ret_1)):
                                ret_3 = nt_strong(t_0_k_0, (t_0_j - 1))
                                if (is_not_empty(ret_3)):
                                    res = bl(ret_1, ret_2, ret_3, ret_4)
                                    answers.append(res)
                                    add_trace(tables, 'strong', t_0_i,t_0_j, 'leftB', t_0_k_0, t_0_j-1, algfct=bl, algparams=[ret_1,ret_2,'x', ret_4], bwdpass=bwdpass)
                t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_ml_comps(t_0_i:int, t_0_j:int, name="ml_comps", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    if (t_0_j-t_0_i >= 5):
        ret_0 = nt_dangle(t_0_i,t_0_j)
        if is_not_empty(ret_0):
            answers.append(ret_0)
            add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps', t_0_i,t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # # sadd(BASE, ml_comps)
    # if (((t_0_j - t_0_i) >= 11)):
    #     if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
    #         ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
    #         if (is_not_empty(ret_1)):
    #             ret_2 = nt_ml_comps((t_0_i + 1), t_0_j)
    #             if (is_not_empty(ret_2)):
    #                 ret_0 = sadd(ret_1, ret_2)
    #                 if (is_not_empty(ret_0)):
    #                     answers.append(ret_0)
    #                     add_trace(tables, 'ml_comps', t_0_i,t_0_j, 'ml_comps', t_0_i+1,t_0_j, algfct=sadd, algparams=[ret_1, 'x'], bwdpass=bwdpass)
    #
    # # cadd(incl(dangle), ml_comps1)
    # if (((t_0_j - t_0_i) >= 10)):
    #     t_0_k_0 = (t_0_i + 5)
    #     while (t_0_k_0 <= (t_0_j - 5)):
    #         ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
    #         if (is_not_empty(ret_6)):
    #             if (((t_0_k_0 - t_0_i) >= 5)):
    #                 ret_5 = nt_dangle(t_0_i, t_0_k_0)
    #                 if (is_not_empty(ret_5)):
    #                     ret_4 = incl(ret_5)
    #                     if (is_not_empty(ret_4)):
    #                         res = cadd(ret_4, ret_6)
    #                         answers.append(res)
    #                         add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'], bwdpass=bwdpass)
    #                         add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x), y), algparams=['x', ret_6], bwdpass=bwdpass)
    #         t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_ml_comps1(t_0_i:int, t_0_j:int, name="ml_comps1", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()

    answers = []

    # sadd(BASE, ml_comps1)
    if (((t_0_j - t_0_i) >= 6)):
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
            if (is_not_empty(ret_1)):
                ret_2 = nt_ml_comps1((t_0_i + 1), t_0_j)
                if (is_not_empty(ret_2)):
                    ret_0 = sadd(ret_1, ret_2)
                    if (is_not_empty(ret_0)):
                        answers.append(ret_0)
                        add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps1', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'], bwdpass=bwdpass)

    # cadd(incl(dangle), ml_comps1)
    if (((t_0_j - t_0_i) >= 10)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= (t_0_j - 5)):
            ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
            if (is_not_empty(ret_6)):
                if (((t_0_k_0 - t_0_i) >= 5)):
                    ret_5 = nt_dangle(t_0_i, t_0_k_0)
                    if (is_not_empty(ret_5)):
                        ret_4 = incl(ret_5)
                        if (is_not_empty(ret_4)):
                            res = cadd(ret_4, ret_6)
                            answers.append(res)
                            add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps1', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'], bwdpass=bwdpass)
                            add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x),y), algparams=['x', ret_6], bwdpass=bwdpass)
            t_0_k_0 += 1

    # incl(dangle)
    if (((t_0_j - t_0_i) >= 5)):
        ret_8 = nt_dangle(t_0_i, t_0_j)
        if (is_not_empty(ret_8)):
            ret_7 = incl(ret_8)
            if (is_not_empty(ret_7)):
                answers.append(ret_7)
                add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_j, algfct=incl, algparams=['x'], bwdpass=bwdpass)

    # addss(incl(dangle), REGION)
    if (((t_0_j - t_0_i) >= 6)):
        t_0_k_1 = (t_0_i + 5)
        while (t_0_k_1 <= (t_0_j - 1)):
            if (unpaired(t_0_seq, t_0_k_1, t_0_j)):
                ret_12 = REGION(t_0_seq, t_0_k_1, t_0_j)
                if (is_not_empty(ret_12)):
                    if (((t_0_k_1 - t_0_i) >= 5)):
                        ret_11 = nt_dangle(t_0_i, t_0_k_1)
                        if (is_not_empty(ret_11)):
                            ret_10 = incl(ret_11)
                            if (is_not_empty(ret_10)):
                                res = addss(ret_10, ret_12)
                                answers.append(res)
                                add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_1, algfct=lambda x,y: addss(incl(x),y), algparams=['x', ret_12], bwdpass=bwdpass)
            t_0_k_1 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_multiloop(t_0_i:int, t_0_j:int, name='multiloop', bwdpass=False) -> float:
    computed()
    if (((t_0_j - t_0_i) < 12)):
        return float_zero;

    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    # ml(BASE, ml_comps, BASE)
    if (((t_0_j - t_0_i) >= 12)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
            if (is_not_empty(ret_3)):
                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                if (is_not_empty(ret_1)):
                    ret_2 = nt_ml_comps((t_0_i + 1), (t_0_j - 1))
                    if (is_not_empty(ret_2)):
                        ret_0 = ml(ret_1, ret_2, ret_3)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)
                            add_trace(tables, 'ml_comps', t_0_i,t_0_j, 'multiloop', t_0_i+1, t_0_j-1, algfct=ml, algparams=[ret_1, 'x', ret_3], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_rightB(t_0_i:int, t_0_j:int, name="rightB", bwdpass=False) -> float:
    computed()
    if (((t_0_j - t_0_i) < 8)):
        return float_zero

    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # br(BASE, strong, REGION, BASE)
    if (((t_0_j - t_0_i) >= 8)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = ((t_0_j - 31)) if (((t_0_j - (t_0_i + 6)) >= 31)) else ((t_0_i + 6))
            while (t_0_k_0 <= (t_0_j - 2)):
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                if (is_not_empty(ret_4)):
                    if ((maxsize(t_0_seq, t_0_k_0, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_0, (t_0_j - 1)))):
                        ret_3 = REGION(t_0_seq, t_0_k_0, (t_0_j - 1))
                        if (is_not_empty(ret_3)):
                            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                            if (is_not_empty(ret_1)):
                                ret_2 = nt_strong((t_0_i + 1), t_0_k_0)
                                if (is_not_empty(ret_2)):
                                    res = br(ret_1, ret_2, ret_3, ret_4)
                                    answers.append(res)
                                    add_trace(tables, 'strong', t_0_i,t_0_j, 'rightB', t_0_i+1, t_0_k_0, algfct=br, algparams=[ret_1, 'x', ret_3, ret_4], bwdpass=bwdpass)
                t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_stack(t_0_i:int, t_0_j:int, name="stack", bwdpass=False) -> float:
    computed()
    if (((t_0_j - t_0_i) < 7)):
       return float_zero;

    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # sr(BASE, weak, BASE)
    if (((t_0_j - t_0_i) >= 7)):
       if (basepair(t_0_seq, t_0_i, t_0_j)):
         ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
         if (is_not_empty(ret_3)):
           ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
           if (is_not_empty(ret_1)):
             ret_2 = nt_weak((t_0_i + 1), (t_0_j - 1))
             if (is_not_empty(ret_2)):
                 ret_0 = sr(ret_1, ret_2, ret_3)
                 if (is_not_empty(ret_0)):
                     answers.append(ret_0)
                     add_trace(tables, 'weak', t_0_i,t_0_j, 'stack', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_1, 'x', ret_3], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_strong(t_0_i:int, t_0_j:int, name="strong", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # sr(BASE, weak, BASE)
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, False)):
        if (((t_0_j - t_0_i) >= 7)):
            if (basepair(t_0_seq, t_0_i, t_0_j)):
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                if (is_not_empty(ret_4)):
                    ret_2 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                    if (is_not_empty(ret_2)):
                        ret_3 = nt_weak((t_0_i + 1), (t_0_j - 1))
                        if (is_not_empty(ret_3)):
                            ret_0 = sr(ret_2, ret_3, ret_4)
                            if (is_not_empty(ret_0)):
                                answers.append(ret_0)
                                add_trace(tables, 'weak', t_0_i,t_0_j, 'strong', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_2, 'x', ret_4], bwdpass=bwdpass)

    # weak
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
        ret_5 = nt_weak(t_0_i, t_0_j)
        if (is_not_empty(ret_5)):
           answers.append(ret_5)
           add_trace(tables, 'weak', t_0_i,t_0_j, 'strong', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_struct(t_0_i:int, name="struct", bwdpass=False) -> float:
    computed()
    t_0_j = 0
    t_0_right_most = len(t_0_seq)
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_right_most))
        incr()
    answers = []

    # sadd(BASE, struct)
    if (((t_0_right_most - t_0_i) >= 1)):
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
            if (is_not_empty(ret_1)):
                ret_2 = nt_struct((t_0_i + 1))
                if (is_not_empty(ret_2)):
                    ret_0 = sadd(ret_1, ret_2)
                    if (is_not_empty(ret_0)):
                        answers.append(ret_0)
                        add_trace(tables, 'struct', t_0_i,t_0_j, 'struct', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'], bwdpass=bwdpass)

    # cadd(dangle, struct)
    if (((t_0_right_most - t_0_i) >= 5)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= t_0_right_most):
            ret_5 = nt_struct(t_0_k_0)
            if (is_not_empty(ret_5)):
                ret_4 = nt_dangle(t_0_i, t_0_k_0)
                if (is_not_empty(ret_4)):
                    ret_0 = cadd(ret_4, ret_5)
                    answers.append(ret_0)
                    add_trace(tables, 'dangle', t_0_i,t_0_j, 'struct', t_0_i, t_0_k_0, algfct=cadd, algparams=['x', ret_5], bwdpass=bwdpass)
                    add_trace(tables, 'struct', t_0_i,t_0_j, 'struct', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'], bwdpass=bwdpass)
            t_0_k_0 += 1

    # nil(LOC)
    if ((((t_0_right_most - t_0_i) >= 0) and ((t_0_right_most - t_0_i) <= 0))):
        ret_7 = LOC(t_0_seq, t_0_i, t_0_i)
        if (is_not_empty(ret_7)):
            ret_6 = nil(ret_7);
            if (is_not_empty(ret_6)):
                answers.append(ret_6)
                add_trace(tables, 'struct', None, None, None, t_0_i, t_0_j, algfct=nil, algparams=[ret_7], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_right_most, eval))
    if name in tables:
        tables[name].set( t_0_i, 0, eval)
        stored()
        return tables[name].get(t_0_i, 0)
    else:
        return eval
def nt_weak(t_0_i:int, t_0_j:int, name="weak", bwdpass=False) -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
           return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []

    # stack
    ret_1 = nt_stack(t_0_i, t_0_j)
    if (is_not_empty(ret_1)):
        answers.append(ret_1)
        if 'stack' in tables:
            add_trace(tables, 'stack', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # hairpin
    ret_2 = nt_hairpin(t_0_i, t_0_j)
    if (is_not_empty(ret_2)):
        answers.append(ret_2)
        add_trace(tables, 'hairpin', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # leftB
    ret_3 = nt_leftB(t_0_i, t_0_j)
    if (is_not_empty(ret_3)):
        answers.append(ret_3)
        add_trace(tables, 'leftB', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # rightB
    ret_4 = nt_rightB(t_0_i, t_0_j)
    if (is_not_empty(ret_4)):
        answers.append(ret_4)
        add_trace(tables, 'rightB', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # iloop
    ret_5 = nt_iloop(t_0_i, t_0_j)
    if (is_not_empty(ret_5)):
        answers.append(ret_5)
        add_trace(tables, 'iloop', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    # multiloop
    ret_6 = nt_multiloop(t_0_i, t_0_j)
    if (is_not_empty(ret_6)):
        answers.append(ret_6)
        add_trace(tables, 'multiloop', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'], bwdpass=bwdpass)

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if name in tables:
        tables[name].set( t_0_i, t_0_j, eval)
        stored()
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval

if True:
    def bt_struct(t_0_i:int, t_0_j, name="struct") -> float:
        t_0_j = 0
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # struct = sadd(BASE, struct) --> struct = sadd(BASE, *struct*)
        if (t_0_i - 1 >= 0):
            if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
                ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                if (is_not_empty(ret_1)):
                    #ret_2 = bt_struct(t_0_i+1, t_0_j) # nt_struct((t_0_i + 1))
                    ret_2 = bt_struct(t_0_i-1, t_0_j) # nt_struct((t_0_i + 1))
                    if (is_not_empty(ret_2)):
                        ret_0 = sadd(ret_1, ret_2)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)

        # struct = cadd(dangle, struct) --> struct = cadd(dangle, *struct*)
        t_0_k_0 = t_0_i - 5
        while (t_0_k_0 >= 0):
            ret_4 = nt_dangle(t_0_k_0, t_0_i)
            if (is_not_empty(ret_4)):
                ret_5 = bt_struct(t_0_k_0, t_0_j)
                if (is_not_empty(ret_5)):
                    ret_0 = cadd(ret_4, ret_5)
                    answers.append(ret_0)
            t_0_k_0 -= 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_hairpin(t_0_i:int, t_0_j, name="hairpin") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = hairpin --> hairpin = *weak*
        if (((t_0_j - t_0_i) >= 5)):
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_weak(t_0_i:int, t_0_j, name="weak") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # strong = weak --> weak = *strong*
        if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
            ret_5 = bt_strong(t_0_i, t_0_j)
            if (is_not_empty(ret_5)):
               answers.append(ret_5)

        # strong = sr(BASE, weak, BASE) --> sr(BASE, *strong*, BASE)
        if (t_0_i-1 >= 0) and (t_0_j <= len(t_0_seq)):
            if (allowLonelyBasepairs(t_0_seq, t_0_i-1, t_0_j+1, False)):
                if (basepair(t_0_seq, t_0_i-1, t_0_j+1)):
                    ret_4 = BASE(t_0_seq, (t_0_j), t_0_j+1)
                    if (is_not_empty(ret_4)):
                        ret_2 = BASE(t_0_seq, t_0_i-1, (t_0_i))
                        if (is_not_empty(ret_2)):
                            ret_3 = bt_strong((t_0_i - 1), (t_0_j + 1))
                            if (is_not_empty(ret_3)):
                                ret_0 = sr(ret_2, ret_3, ret_4)
                                if (is_not_empty(ret_0)):
                                    answers.append(ret_0)

        # stack = sr(BASE, weak, BASE) --> weak = sr(BASE, stack, BASE)
        if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):
           if (basepair(t_0_seq, t_0_i-1, t_0_j+1)):
             ret_3 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
             if (is_not_empty(ret_3)):
               ret_1 = BASE(t_0_seq, t_0_i-1, (t_0_i ))
               if (is_not_empty(ret_1)):
                 ret_2 = bt_stack((t_0_i - 1), (t_0_j + 1))
                 if (is_not_empty(ret_2)):
                     ret_0 = sr(ret_1, ret_2, ret_3)
                     if (is_not_empty(ret_0)):
                         answers.append(ret_0)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_strong(t_0_i:int, t_0_j, name="strong") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # dangle = drem(LOC, strong, LOC) --> strong = drem(LOC, *dangle*, LOC)
        ret_3 = LOC(t_0_seq, t_0_j, t_0_j)
        if (is_not_empty(ret_3)):
            ret_1 = LOC(t_0_seq, t_0_i, t_0_i)
            if (is_not_empty(ret_1)):
                ret_2 = bt_dangle(t_0_i, t_0_j)
                if (is_not_empty(ret_2)):
                    ret_0 = drem(ret_1, ret_2, ret_3)
                    if (is_not_empty(ret_0)):
                        answers.append(ret_0)

        # leftB = bl(BASE, REGION, strong, BASE) --> strong = bl(BASE, REGION, leftB, BASE)
        if t_0_j+1 <= len(t_0_seq):
            t_0_k_0 = t_0_i - 2
            while (t_0_k_0 >= 0):
                if (basepair(t_0_seq, t_0_k_0, t_0_j+1)):
                    ret_4 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
                    if (is_not_empty(ret_4)):
                        if ((maxsize(t_0_seq, (t_0_k_0 + 1), t_0_i, 30) and unpaired(t_0_seq, (t_0_k_0 + 1), t_0_i))):
                            ret_2 = REGION(t_0_seq, (t_0_k_0 + 1), t_0_i)
                            if (is_not_empty(ret_2)):
                                ret_1 = BASE(t_0_seq, t_0_k_0, (t_0_k_0 + 1))
                                if (is_not_empty(ret_1)):
                                    ret_3 = bt_leftB(t_0_k_0, t_0_j + 1)
                                    if (is_not_empty(ret_3)):
                                        res = bl(ret_1, ret_2, ret_3, ret_4)
                                        answers.append(res)
                t_0_k_0 -= 1

        # br(BASE, strong, REGION, BASE) --> strong = br(BASE, rightB, REGION, BASE)
        if (t_0_i-1 >= 0):
            t_0_k_0 = t_0_j + 2
            while (t_0_k_0 <= len(t_0_seq)):
                if (basepair(t_0_seq, t_0_i-1, t_0_k_0)):
                    ret_4 = BASE(t_0_seq, (t_0_k_0 - 1), t_0_k_0)
                    if (is_not_empty(ret_4)):
                        if ((maxsize(t_0_seq, t_0_j, (t_0_k_0 - 1), 30) and unpaired(t_0_seq, t_0_j, (t_0_k_0 - 1)))):
                            ret_3 = REGION(t_0_seq, t_0_j, (t_0_k_0 - 1))
                            if (is_not_empty(ret_3)):
                                ret_1 = BASE(t_0_seq, t_0_i-1, (t_0_i))
                                if (is_not_empty(ret_1)):
                                    ret_2 = bt_rightB((t_0_i - 1), t_0_k_0)
                                    if (is_not_empty(ret_2)):
                                        res = br(ret_1, ret_2, ret_3, ret_4)
                                        answers.append(res)
                t_0_k_0 += 1

        # iloop -> il(BASE, REGION, strong, REGION, BASE) --> strong = il(BASE, REGION, iloop, REGION, BASE)
        t_0_k_0 = t_0_i - 2
        while (t_0_k_0 >= 0):
            t_0_k_1 = t_0_j + 2
            while (t_0_k_1 <= len(t_0_seq)):
                if (basepair(t_0_seq, t_0_k_0, t_0_k_1)):
                    ret_5 = BASE(t_0_seq, (t_0_k_1 - 1), t_0_k_1)
                    if (is_not_empty(ret_5)):
                        if ((maxsize(t_0_seq, t_0_j, (t_0_k_1 - 1), 30) and unpaired(t_0_seq, t_0_j, (t_0_k_1 - 1)))):
                            ret_4 = REGION(t_0_seq, t_0_j, (t_0_k_1 - 1))
                            if (is_not_empty(ret_4)):
                                if ((maxsize(t_0_seq, (t_0_k_0 + 1), t_0_i, 30) and unpaired(t_0_seq, (t_0_k_0 + 1), t_0_i))):
                                    ret_2 = REGION(t_0_seq, (t_0_k_0+1), t_0_i)
                                    if (is_not_empty(ret_2)):
                                        ret_1 = BASE(t_0_seq, t_0_k_0, (t_0_k_0 + 1))
                                        if (is_not_empty(ret_1)):
                                            ret_3 = bt_iloop(t_0_k_0, t_0_k_1)
                                            if (is_not_empty(ret_3)):
                                                res = il(ret_1, ret_2, ret_3, ret_4, ret_5)
                                                answers.append(res)
                t_0_k_1 += 1
            t_0_k_0 -= 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_dangle(t_0_i:int, t_0_j, name="dangle") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # struct = cadd(dangle, struct) --> dangle = cadd(*struct*, struct)
        if (t_0_j - t_0_i) >= 5:
            ret_5 = nt_struct(t_0_j)
            if (is_not_empty(ret_5)):
                ret_4 = bt_struct(t_0_i, t_0_j)
                if (is_not_empty(ret_4)):
                    ret_0 = cadd(ret_4, ret_5)
                    answers.append(ret_0)

        # tmp: ml_comps = dangle --> dangle = ml_comps
        if (t_0_j-t_0_i >= 5):
            ret_0 = bt_ml_comps(t_0_i,t_0_j)
            if is_not_empty(ret_0):
                answers.append(ret_0)


        # # ml_comps = cadd(incl(dangle), ml_comps1) --> dangle = cadd(incl(*ml_comps*), ml_comps1)
        # if (t_0_j - t_0_i) >= 5:
        #     ret_5 = nt_ml_comps1(t_0_i, t_0_j)
        #     if (is_not_empty(ret_5)):
        #         ret_4 = bt_ml_comps(t_0_i, t_0_j)
        #         if (is_not_empty(ret_4)):
        #             ret_6 = incl(ret_4)
        #             if (is_not_empty(ret_6)):
        #                 ret_0 = cadd(ret_6, ret_5)
        #                 answers.append(ret_0)
        #
        # # ml_comps1 = cadd(incl(dangle), ml_comps1) --> cadd(incl(*ml_comps1*), ml_comps1)
        # if (t_0_j - t_0_i) >= 5:
        #     ret_5 = nt_ml_comps1(t_0_i, t_0_j)
        #     if (is_not_empty(ret_5)):
        #         ret_4 = bt_ml_comps1(t_0_i, t_0_j)
        #         if (is_not_empty(ret_4)):
        #             ret_6 = incl(ret_4)
        #             if (is_not_empty(ret_6)):
        #                 ret_0 = cadd(ret_6, ret_5)
        #                 answers.append(ret_0)
        #
        # # ml_comps1 = incl(dangle) --> dangle = incl(*ml_comps1*)
        # ret_0 = bt_ml_comps1(t_0_i, t_0_j)
        # if (is_not_empty(ret_0)):
        #     ret_1 = incl(ret_0)
        #     if (is_not_empty(ret_1)):
        #         answers.append(ret_1)
        #
        # # ml_comps1 = addss(incl(dangle), REGION) --> dangle = addss(incl(*ml_comp1*), REGION)
        # t_0_k_1 = t_0_j + 1
        # while (t_0_k_1 <= len(t_0_seq)):
        #     if (unpaired(t_0_seq, t_0_j, t_0_k_1)):
        #         ret_12 = REGION(t_0_seq, t_0_j, t_0_k_1)
        #         if (is_not_empty(ret_12)):
        #             ret_11 = bt_ml_comps1(t_0_i, t_0_k_1)
        #             if (is_not_empty(ret_11)):
        #                 ret_10 = incl(ret_11)
        #                 if (is_not_empty(ret_10)):
        #                     res = addss(ret_10, ret_12)
        #                     answers.append(res)
        #     t_0_k_1 += 1


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_stack(t_0_i:int, t_0_j, name="stack") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = hairpin --> hairpin = *weak*
        if (((t_0_j - t_0_i) >= 7)):
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_leftB(t_0_i:int, t_0_j, name="leftB") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = leftB --> leftB = *weak*
        if (((t_0_j - t_0_i) >= 8)):
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_rightB(t_0_i:int, t_0_j, name="rightB") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = rightB --> rightB = *weak*
        if (((t_0_j - t_0_i) >= 8)):
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_iloop(t_0_i:int, t_0_j, name="iloop") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = iloop --> iloop = *weak*
        if (((t_0_j - t_0_i) >= 9)):
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_multiloop(t_0_i:int, t_0_j, name="multiloop") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # weak = multiloop --> multiloop = *weak*
        if t_0_j - t_0_i >= 12:
            ret_2 = bt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_2)):
                answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    def bt_ml_comps(t_0_i:int, t_0_j, name="ml_comps") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
                return tables[name].bt_get_v2(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # multiloop = ml(BASE, ml_comps, BASE) --> ml_comps = ml(BASE, multiloop, BASE)
        if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):
            if (basepair(t_0_seq, t_0_i-1, t_0_j+1)):
                ret_3 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
                if (is_not_empty(ret_3)):
                    ret_1 = BASE(t_0_seq, t_0_i-1, (t_0_i ))
                    if (is_not_empty(ret_1)):
                        ret_2 = bt_multiloop((t_0_i - 1), (t_0_j + 1))
                        if (is_not_empty(ret_2)):
                            ret_0 = ml(ret_1, ret_2, ret_3)
                            if (is_not_empty(ret_0)):
                                answers.append(ret_0)

        # # multiloop = ml(BASE, ml_comps, BASE) --> ml_comps = ml(BASE, multiloop, BASE)
        # if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):
        #     if (basepair(t_0_seq, t_0_i-1, t_0_j+1)):
        #         ret_3 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
        #         if (is_not_empty(ret_3)):
        #             ret_1 = BASE(t_0_seq, t_0_i-1, (t_0_i ))
        #             if (is_not_empty(ret_1)):
        #                 ret_2 = bt_multiloop((t_0_i - 1), (t_0_j + 1))
        #                 if (is_not_empty(ret_2)):
        #                     ret_0 = ml(ret_1, ret_2, ret_3)
        #                     if (is_not_empty(ret_0)):
        #                         answers.append(ret_0)
        #
        # # productions:
        # # ml_comps = sadd(BASE, ml_comps) --> ml_comps = sadd(BASE, *ml_comps*)
        # if (t_0_i - 1 >= 0):
        #     if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
        #         ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
        #         if (is_not_empty(ret_1)):
        #             ret_2 = bt_ml_comps(t_0_i-1, t_0_j) # nt_struct((t_0_i + 1))
        #             if (is_not_empty(ret_2)):
        #                 ret_0 = sadd(ret_1, ret_2)
        #                 if (is_not_empty(ret_0)):
        #                     answers.append(ret_0)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set_v2( t_0_i, t_0_j, eval)
            return tables[name].bt_get_v2(t_0_i, t_0_j)
        else:
            return eval
    # def bt_ml_comps1(t_0_i:int, t_0_j, name="ml_comps1") -> float:
    #     if name in tables:
    #         if (tables[name].bt_is_tabulated_v2(t_0_i, t_0_j)):
    #             if PRINTBTSTACK:
    #                 print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get_v2(t_0_i, t_0_j)))
    #             return tables[name].bt_get_v2(t_0_i, t_0_j)
    #     if PRINTBTSTACK:
    #         print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
    #         incr()
    #
    #     answers = []
    #
    #     # productions:
    #     # ml_comps = cadd(incl(dangle), ml_comps1) --> ml_comps1 = cadd(incl(dangle), *ml_comps*)
    #     t_0_k_0 = t_0_i - 5
    #     while (t_0_k_0 >= 0):
    #         ret_4 = nt_dangle(t_0_k_0, t_0_i)
    #         if (is_not_empty(ret_4)):
    #             ret_6 = incl(ret_4)
    #             if (is_not_empty(ret_6)):
    #                 ret_5 = bt_ml_comps(t_0_k_0, t_0_j)
    #                 if (is_not_empty(ret_5)):
    #                     ret_0 = cadd(ret_6, ret_5)
    #                     answers.append(ret_0)
    #         t_0_k_0 -= 1
    #
    #     # ml_comps1 = cadd(incl(dangle), ml_comps1) --> ml_comps1 = cadd(incl(dangle), *ml_comps1*)
    #     t_0_k_0 = t_0_i - 5
    #     while (t_0_k_0 >= 0):
    #         ret_4 = nt_dangle(t_0_k_0, t_0_i)
    #         if (is_not_empty(ret_4)):
    #             ret_6 = incl(ret_4)
    #             if (is_not_empty(ret_6)):
    #                 ret_5 = bt_ml_comps1(t_0_k_0, t_0_j)
    #                 if (is_not_empty(ret_5)):
    #                     ret_0 = cadd(ret_6, ret_5)
    #                     answers.append(ret_0)
    #         t_0_k_0 -= 1
    #
    #     # ml_comps1 = sadd(BASE, ml_comps1) --> ml_comps1 = sadd(BASE, *ml_comps1*)
    #     if (t_0_i - 1 >= 0):
    #         if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
    #             ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
    #             if (is_not_empty(ret_1)):
    #                 ret_2 = bt_ml_comps1(t_0_i-1, t_0_j) # nt_struct((t_0_i + 1))
    #                 if (is_not_empty(ret_2)):
    #                     ret_0 = sadd(ret_1, ret_2)
    #                     if (is_not_empty(ret_0)):
    #                         answers.append(ret_0)
    #
    #     eval = h(answers)
    #     if PRINTBTSTACK:
    #         decr()
    #         print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    #     if name in tables:
    #         tables[name].bt_set_v2( t_0_i, t_0_j, eval)
    #         return tables[name].bt_get_v2(t_0_i, t_0_j)
    #     else:
    #         return eval

if False:
    def bt_dangle(t_0_i:int, t_0_j:int, name="dangle") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # Productions
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
            #print(answers, h(answers))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_hairpin(t_0_i:int, t_0_j:int, name="hairpin") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_iloop(t_0_i:int, t_0_j:int, name="iloop") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_leftB(t_0_i:int, t_0_j:int, name="leftB") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # bl(BASE, REGION, strong, BASE) -> weak
        if (((t_0_j - t_0_i) >= 8)):
            if (basepair(t_0_seq, t_0_i, t_0_j)):
                t_0_k_0 = t_0_i + 2
                while ((t_0_k_0 <= (t_0_j - 6)) and (t_0_k_0 <= (t_0_i + 31))):
                    ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                    if (is_not_empty(ret_4)):
                        if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
                            ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
                            if (is_not_empty(ret_2)):
                                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                                if (is_not_empty(ret_1)):
                                    ret_3 = 0# bt_weak(t_0_i, t_0_j)
                                    print("calling weak(%i,%i) region(%i,%i)" % (t_0_i, t_0_j, t_0_i+1, t_0_k_0))
                                    if (is_not_empty(ret_3)):
                                        res = bl(ret_1, ret_2, ret_3, ret_4)
                                        answers.append(res)
                    t_0_k_0 += 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_ml_comps(t_0_i:int, t_0_j:int, name="ml_comps") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_ml_comps1(t_0_i:int, t_0_j:int, name="ml_comps1") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        sys.exit()

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_multiloop(t_0_i:int, t_0_j:int, name="multiloop") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_rightB(t_0_i:int, t_0_j:int, name="rightB") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_stack(t_0_i:int, t_0_j:int, name="stack") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        sys.exit()



        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_strong(t_0_i:int, t_0_j:int, name="strong") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        sys.exit()

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_struct(t_0_i:int, t_0_j:int, name="struct") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        sys.exit()


        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval
    def bt_weak(t_0_i:int, t_0_j:int, name="weak") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        sys.exit()

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, 0, eval)
            return tables[name].bt_get(t_0_i, 0)
        else:
            return eval

if False:
    def bt_dangle(t_0_i:int, t_0_j:int, name="dangle") -> float:
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # cadd(incl(*dangle*), ml_comps1)
        if (((t_0_j + 5) <= len(t_0_seq))):
            t_0_k_0 = (t_0_j + 5)
            while (t_0_k_0 <= len(t_0_seq)):
                ret_6 = nt_ml_comps1(t_0_j, t_0_k_0)
                if (is_not_empty(ret_6)):
                    if (((t_0_k_0 - t_0_i) >= 10)):
                        ret_5 = bt_ml_comps1(t_0_i, t_0_k_0) # nt_dangle(t_0_i, t_0_k_0)
                        if (is_not_empty(ret_5)):
                            ret_4 = incl(ret_5)
                            if (is_not_empty(ret_4)):
                                res = cadd(ret_4, ret_6)
                                answers.append(res)
                t_0_k_0 += 1

        # incl(*dangle*)
        #if (((t_0_j - t_0_i) >= 5)):
        ret_8 = bt_ml_comps1(t_0_i, t_0_j) # nt_dangle(t_0_i, t_0_j)
        if (is_not_empty(ret_8)):
            ret_7 = incl(ret_8)
            if (is_not_empty(ret_7)):
                answers.append(ret_7)

        # addss(incl(*dangle*), REGION)
        if (((t_0_j + 1) <= len(t_0_seq))):
            t_0_k_1 = (t_0_j + 1)
            while (t_0_k_1 <= len(t_0_seq)):
                if (unpaired(t_0_seq, t_0_j, t_0_k_1)):
                    ret_12 = REGION(t_0_seq, t_0_j, t_0_k_1)
                    if (is_not_empty(ret_12)):
                        if (((t_0_k_1 - t_0_i) >= 6)):
                            ret_11 = bt_ml_comps1(t_0_i, t_0_k_1) # nt_dangle(t_0_i, t_0_k_1)
                            if (is_not_empty(ret_11)):
                                ret_10 = incl(ret_11)
                                if (is_not_empty(ret_10)):
                                    res = addss(ret_10, ret_12)
                                    answers.append(res)
                t_0_k_1 += 1

        # cadd(incl(*dangle*), ml_comps1)
        if (((t_0_j + 5) <= len(t_0_seq))):
            t_0_k_0 = (t_0_j + 5)
            while (t_0_k_0 <= len(t_0_seq)):
                ret_6 = nt_ml_comps1(t_0_j, t_0_k_0)
                if (is_not_empty(ret_6)):
                    if (((t_0_k_0 - t_0_i) >= 10)):
                        ret_5 = bt_ml_comps(t_0_i, t_0_k_0) # nt_dangle(t_0_i, t_0_k_0)
                        if (is_not_empty(ret_5)):
                            ret_4 = incl(ret_5)
                            if (is_not_empty(ret_4)):
                                res = cadd(ret_4, ret_6)
                                answers.append(res)
                t_0_k_0 += 1

        # cadd(*dangle*, struct)
        # if True: # (t_0_i <= 0) and (t_0_j >= len(t_0_seq)):
        #     ret_5 = nt_struct(t_0_j)
        #     if (is_not_empty(ret_5)):
        #         ret_4 = bt_struct(t_0_i)
        #         if (is_not_empty(ret_4)):
        #             ret_0 = cadd(ret_4, ret_5)
        #             answers.append(ret_0)
        # if (((t_0_j + 0 <= len(t_0_seq)))):
            t_0_k_0 = t_0_j
            while (t_0_k_0 <= len(t_0_seq)):
                #print("nt_dangle(%i,%i) = " % (t_0_i,t_0_j))
                ret_5 = nt_struct(t_0_k_0)
                if (is_not_empty(ret_5)):
                    ret_4 = bt_struct(t_0_k_0) # nt_dangle(t_0_i, t_0_k_0)
                    if (is_not_empty(ret_4)):
                        ret_0 = cadd(ret_4, ret_5)
                        answers.append(ret_0)
                t_0_k_0 += 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
            #print(answers, h(answers))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_hairpin(t_0_i:int, t_0_j:int, name="hairpin") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # hl(BASE, REGION, BASE)
        if False:
            if (((t_0_j - t_0_i) >= 5)):
                if (basepair(t_0_seq, t_0_i, t_0_j)):
                    ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                    if (is_not_empty(ret_3)):
                        ret_2 = np.nan;
                        if ((minsize(t_0_seq, (t_0_i + 1), (t_0_j - 1), 3) and unpaired(t_0_seq, (t_0_i + 1), (t_0_j - 1)))):
                            ret_2 = REGION(t_0_seq, (t_0_i + 1), (t_0_j - 1))
                            if (is_not_empty(ret_2)):
                                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                                if (is_not_empty(ret_1)):
                                    ret_0 = hl(ret_1, ret_2, ret_3)
                                    if (is_not_empty(ret_0)):
                                        answers.append(ret_0)

        # hairpin
        ret_2 = bt_weak(t_0_i, t_0_j) #nt_hairpin(t_0_i, t_0_j)
        if (is_not_empty(ret_2)):
            answers.append(ret_2)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_iloop(t_0_i:int, t_0_j:int, name="iloop") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # iloop
        ret_5 = bt_weak(t_0_i, t_0_j) #nt_iloop(t_0_i, t_0_j)
        if (is_not_empty(ret_5)):
            answers.append(ret_5)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_leftB(t_0_i:int, t_0_j:int, name="leftB") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # leftB
        ret_3 = bt_weak(t_0_i, t_0_j) #nt_leftB(t_0_i, t_0_j)
        if (is_not_empty(ret_3)):
            answers.append(ret_3)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_ml_comps(t_0_i:int, t_0_j:int, name="ml_comps") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        # sadd(BASE, *ml_comps*)
        #if (((t_0_j - t_0_i-1) >= 11)):
        if (t_0_i-1 >= 0):
            if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
                ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                if (is_not_empty(ret_1)):
                    ret_2 = bt_ml_comps(t_0_i-1,t_0_j) #nt_ml_comps((t_0_i + 1), t_0_j)
                    if (is_not_empty(ret_2)):
                        ret_0 = sadd(ret_1, ret_2)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)

        # ml(BASE, ml_comps, BASE)
        #if (((t_0_j - t_0_i-2) >= 12)):
        if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):
            if (basepair(t_0_seq, t_0_i-1, t_0_j+1)):
                ret_3 = BASE(t_0_seq, t_0_j, t_0_j+1)
                if (is_not_empty(ret_3)):
                    ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                    if (is_not_empty(ret_1)):
                        #ret_2 = bt_multiloop(t_0_i+1, t_0_j-1) # nt_ml_comps((t_0_i + 1), (t_0_j - 1))
                        ret_2 = bt_multiloop(t_0_i-1, t_0_j+1) # nt_ml_comps((t_0_i + 1), (t_0_j - 1))
                        if (is_not_empty(ret_2)):
                            ret_0 = ml(ret_1, ret_2, ret_3)
                            if (is_not_empty(ret_0)):
                                answers.append(ret_0)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_ml_comps1(t_0_i:int, t_0_j:int, name="ml_comps1") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        # cadd(incl(dangle), *ml_comps1*)
        if (((t_0_j - t_0_i) >= 5)):
            t_0_k_0 = t_0_i - 5
            while (t_0_k_0 >= 0):
                ret_6 = bt_ml_comps(t_0_k_0, t_0_j) # nt_ml_comps1(t_0_k_0, t_0_j)
                if (is_not_empty(ret_6)):
                    if (((t_0_k_0 - t_0_j) >= 10)):
                        ret_5 = nt_dangle(t_0_k_0, t_0_i)
                        if (is_not_empty(ret_5)):
                            ret_4 = incl(ret_5)
                            if (is_not_empty(ret_4)):
                                res = cadd(ret_4, ret_6)
                                answers.append(res)
                t_0_k_0 -= 1

        # sadd(BASE, *ml_comps1*)
        #if (((t_0_j - t_0_i-1) >= 6)):
        if (t_0_i-1 >= 0):
            if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
                ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                if (is_not_empty(ret_1)):
                    #ret_2 = bt_ml_comps1(t_0_i+1, t_0_j) # nt_ml_comps1((t_0_i + 1), t_0_j)
                    ret_2 = bt_ml_comps1(t_0_i-1, t_0_j) # nt_ml_comps1((t_0_i + 1), t_0_j)
                    if (is_not_empty(ret_2)):
                        ret_0 = sadd(ret_1, ret_2)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)

        # cadd(incl(dangle), *ml_comps1*)
        if (((t_0_j - t_0_i) >= 10)):
            t_0_k_0 = (t_0_i - 5)
            while (t_0_k_0 >= 0):
                ret_6 = bt_ml_comps1(t_0_k_0, t_0_j) # nt_ml_comps1(t_0_k_0, t_0_j)
                if (is_not_empty(ret_6)):
                    if (((t_0_k_0 - t_0_i) >= 5)):
                        ret_5 = nt_dangle(t_0_k_0, t_0_i)
                        if (is_not_empty(ret_5)):
                            ret_4 = incl(ret_5)
                            if (is_not_empty(ret_4)):
                                res = cadd(ret_4, ret_6)
                                answers.append(res)
                t_0_k_0 -= 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_multiloop(t_0_i:int, t_0_j:int, name="multiloop") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        # multiloop
        ret_6 = bt_weak(t_0_i, t_0_j) # nt_multiloop(t_0_i, t_0_j)
        if (is_not_empty(ret_6)):
            answers.append(ret_6)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_rightB(t_0_i:int, t_0_j:int, name="rightB") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        # rightB
        ret_4 = bt_weak(t_0_i, t_0_j) # nt_rightB(t_0_i, t_0_j)
        if (is_not_empty(ret_4)):
            answers.append(ret_4)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_stack(t_0_i:int, t_0_j:int, name="stack") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)

        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()
        answers = []

        # productions:
        # stack
        ret_1 = bt_weak(t_0_i, t_0_j) # nt_stack(t_0_i, t_0_j)
        if (is_not_empty(ret_1)):
            answers.append(ret_1)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_strong(t_0_i:int, t_0_j:int, name="strong") -> float:
        #if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
        #    return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # drem(LOC, *strong*, LOC)
        if (t_0_i >= 0) and (t_0_j <= len(t_0_seq)):
            ret_3 = LOC(t_0_seq, t_0_j, t_0_j)
            if (is_not_empty(ret_3)):
                ret_1 = LOC(t_0_seq, t_0_i, t_0_i)
                if (is_not_empty(ret_1)):
                    ret_2 = bt_dangle(t_0_i, t_0_j) # nt_strong(t_0_i, t_0_j)
                    if (is_not_empty(ret_2)):
                        ret_0 = drem(ret_1, ret_2, ret_3)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)

        # br(BASE, *strong*, REGION, BASE)
        #if (((t_0_j - t_0_i) >= 8)):
        if (t_0_i-1 >= 0):
            t_0_k_0 = t_0_j+2 #((t_0_j - 31)) if (((t_0_j - (t_0_i + 6)) >= 31)) else ((t_0_i + 6))
            while (t_0_k_0 <= len(t_0_seq)):
                if (basepair(t_0_seq, t_0_i-1, t_0_k_0)):
                    ret_4 = BASE(t_0_seq, (t_0_k_0 - 1), t_0_k_0)
                    if (is_not_empty(ret_4)):
                        if (maxsize(t_0_seq, t_0_j, t_0_k_0 - 1, 30) and unpaired(t_0_seq, t_0_j, t_0_k_0-1)):
                            ret_3 = REGION(t_0_seq, t_0_j, t_0_k_0 - 1)
                            if (is_not_empty(ret_3)):
                                ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                                if (is_not_empty(ret_1)):
                                    #ret_2 = bt_rightB(t_0_i+1, t_0_k_0) # nt_strong((t_0_i + 1), t_0_k_0)
                                    ret_2 = bt_rightB(t_0_i-1, t_0_k_0) # nt_strong((t_0_i + 1), t_0_k_0)
                                    if (is_not_empty(ret_2)):
                                        res = br(ret_1, ret_2, ret_3, ret_4)
                                        answers.append(res)
                t_0_k_0 += 1

        # il(BASE, REGION, *strong*, REGION, BASE)
        if (((t_0_j - t_0_i) >= 9)):
            t_0_k_0 = t_0_i - 2
            while ((t_0_k_0 >= 0) and (t_0_i - t_0_k_0 <= 31)):
                t_0_k_1 = t_0_j+2#((t_0_j - 31)) if (((t_0_j - (t_0_k_0 + 5)) >= 31)) else ((t_0_k_0 + 5))
                while (t_0_k_1 <= len(t_0_seq)):
                    if (basepair(t_0_seq, t_0_k_0, t_0_k_1)):
                        ret_5 = BASE(t_0_seq, t_0_k_1-1, t_0_k_1)
                        if (is_not_empty(ret_5)):
                            if ((maxsize(t_0_seq, t_0_j, t_0_k_1 - 1, 30) and unpaired(t_0_seq, t_0_j, t_0_k_1 - 1))):
                                ret_4 = REGION(t_0_seq, t_0_j, t_0_k_1 - 1)
                                if (is_not_empty(ret_4)):
                                    if ((maxsize(t_0_seq, t_0_k_0+1, t_0_i, 30) and unpaired(t_0_seq, (t_0_k_0 + 1), t_0_i))):
                                        ret_2 = REGION(t_0_seq, t_0_k_0 + 1, t_0_i)
                                        if (is_not_empty(ret_2)):
                                            ret_1 = BASE(t_0_seq, t_0_k_0, t_0_k_0 + 1)
                                            if (is_not_empty(ret_1)):
                                                ret_3 = bt_iloop(t_0_k_0, t_0_k_1) # nt_strong(t_0_k_0, t_0_k_1)
                                                if (is_not_empty(ret_3)):
                                                    res = il(ret_1, ret_2, ret_3, ret_4, ret_5)
                                                    answers.append(res)
                    t_0_k_1 += 1
                t_0_k_0 -= 1

        # bl(BASE, REGION, *strong*, BASE)
        if (t_0_i-2 >= 0) and (t_0_j+1 <= len(t_0_seq)):# (((t_0_j - t_0_i) >= 8)):
            t_0_k_0 = t_0_i - 2
            while (t_0_k_0 >= 0):#((t_0_k_0 <= (t_0_j - 6)) and (t_0_k_0 <= (t_0_i + 31))):
                if (basepair(t_0_seq, t_0_k_0, t_0_j+1)):
                    ret_4 = BASE(t_0_seq, t_0_j, t_0_j+1)
                    if (is_not_empty(ret_4)):
                        if ((maxsize(t_0_seq, (t_0_k_0 + 1), t_0_i, 30) and unpaired(t_0_seq, (t_0_k_0 + 1), t_0_i))):
                            ret_2 = REGION(t_0_seq, (t_0_k_0 + 1), t_0_i)
                            if (is_not_empty(ret_2)):
                                ret_1 = BASE(t_0_seq, t_0_k_0, t_0_k_0 + 1)
                                if (is_not_empty(ret_1)):
                                    #ret_3 = bt_leftB(t_0_k_0, t_0_j-1) # nt_strong(t_0_k_0, (t_0_j - 1))
                                    ret_3 = bt_leftB(t_0_k_0, t_0_j+1) # nt_strong(t_0_k_0, (t_0_j - 1))
                                    if (is_not_empty(ret_3)):
                                        res = bl(ret_1, ret_2, ret_3, ret_4)
                                        answers.append(res)
                t_0_k_0 -= 1

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_struct(t_0_i:int, name="struct") -> float:
        t_0_j = 0
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # sadd(BASE, *struct*)
        if (t_0_i - 1 >= 0):
            if (unpaired(t_0_seq, t_0_i-1, t_0_i)):
                ret_1 = BASE(t_0_seq, t_0_i-1, t_0_i)
                if (is_not_empty(ret_1)):
                    #ret_2 = bt_struct(t_0_i+1, t_0_j) # nt_struct((t_0_i + 1))
                    ret_2 = bt_struct(t_0_i-1) # nt_struct((t_0_i + 1))
                    if (is_not_empty(ret_2)):
                        ret_0 = sadd(ret_1, ret_2)
                        if (is_not_empty(ret_0)):
                            answers.append(ret_0)

        # # cadd(dangle, *struct*)
        # ret_5 = nt_struct(t_0_j)
        # if (is_not_empty(ret_5)):
        #     ret_4 = bt_struct(t_0_i)
        #     if (is_not_empty(ret_4)):
        #         ret_0 = cadd(ret_4, ret_5)
        #         answers.append(ret_0)
        if (((t_0_i-5) >= 0)):
            t_0_k_0 = (t_0_i - 5)
            while (t_0_k_0 >= 0):
                ret_5 = bt_struct(t_0_k_0) # nt_struct(t_0_k_0)
                if (is_not_empty(ret_5)):
                    ret_4 = nt_dangle(t_0_k_0, t_0_i)
                    if (is_not_empty(ret_4)):
                        ret_0 = cadd(ret_4, ret_5)
                        answers.append(ret_0)
                t_0_k_0 -= 1

        # # # nil(LOC)
        # if (t_0_i == 0):
        #     ret_7 = LOC(t_0_seq, t_0_i, t_0_i)
        #     if (is_not_empty(ret_7)):
        #         ret_6 = nil(ret_7);
        #         if (is_not_empty(ret_6)):
        #             answers.append(ret_6)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval
    def bt_weak(t_0_i:int, t_0_j:int, name="weak") -> float:
        if (t_0_i < 0) or (t_0_j > len(t_0_seq)):
            return float_zero
        if name in tables:
            if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
                if PRINTBTSTACK:
                    print("%sretrieved bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, tables[name].bt_get(t_0_i, t_0_j)))
                return tables[name].bt_get(t_0_i, t_0_j)
        if PRINTBTSTACK:
            print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
            incr()

        answers = []

        # productions:
        # sr(BASE, weak, BASE)
        if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, False)):
            if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):#if (((t_0_j - t_0_i) >= 7)):
                if (basepair(t_0_seq, t_0_i, t_0_j)):
                    ret_4 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
                    if (is_not_empty(ret_4)):
                        ret_2 = BASE(t_0_seq, t_0_i-1, (t_0_i ))
                        if (is_not_empty(ret_2)):
                            #ret_3 = bt_strong(t_0_i+1, t_0_j-1) # nt_weak((t_0_i + 1), (t_0_j - 1))
                            ret_3 = bt_strong(t_0_i-1, t_0_j+1) # nt_weak((t_0_i + 1), (t_0_j - 1))
                            if (is_not_empty(ret_3)):
                                ret_0 = sr(ret_2, ret_3, ret_4)
                                if (is_not_empty(ret_0)):
                                    answers.append(ret_0)

        # weak
        if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
            ret_5 = bt_strong(t_0_i, t_0_j) # nt_weak(t_0_i, t_0_j)
            if (is_not_empty(ret_5)):
               answers.append(ret_5)

        # sr(BASE, *weak*, BASE)
        if (t_0_i-1 >= 0) and (t_0_j+1 <= len(t_0_seq)):#(((t_0_j - t_0_i) >= 7)):
           if (basepair(t_0_seq, t_0_i, t_0_j)):
             ret_3 = BASE(t_0_seq, (t_0_j ), t_0_j+1)
             if (is_not_empty(ret_3)):
               ret_1 = BASE(t_0_seq, t_0_i-1, (t_0_i))
               if (is_not_empty(ret_1)):
                 #ret_2 = bt_stack(t_0_i+1, t_0_j-1) # nt_weak((t_0_i + 1), (t_0_j - 1))
                 ret_2 = bt_stack(t_0_i-1, t_0_j+1) # nt_weak((t_0_i + 1), (t_0_j - 1))
                 if (is_not_empty(ret_2)):
                     ret_0 = sr(ret_1, ret_2, ret_3)
                     if (is_not_empty(ret_0)):
                         answers.append(ret_0)

        eval = h(answers)
        if PRINTBTSTACK:
            decr()
            print("%s} set bt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
        if name in tables:
            tables[name].bt_set( t_0_i, t_0_j, eval)
            return tables[name].bt_get(t_0_i, t_0_j)
        else:
            return eval


msg = "Function '%s' for algebra '%s' is not implemented (yet?)!"
def addss(x:float, r:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale((r.j - r.i)) * x) * mk_pf(ss_energy(r)))
    elif ALGEBRA == 'mfe':
        return x + ss_energy(r)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('addss', ALGEBRA))
def bl(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale(((2 + lr.j) - lr.i)) * x) * mk_pf(bl_energy(lr, rb)))
    elif ALGEBRA == 'mfe':
        return x + bl_energy(lr, rb)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('bl', ALGEBRA))
def br(lb:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale(((2 + rr.j) - rr.i)) * x) * mk_pf(br_energy(lb, rr)))
    elif ALGEBRA == 'mfe':
        return x + br_energy(lb, rr)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('br', ALGEBRA))
def cadd(x:float, y:float):
    if ALGEBRA == 'pfunc':
        return x*y
    elif ALGEBRA == 'mfe':
        return x + y
    elif ALGEBRA == 'count':
        return x * y
    else:
        raise ValueError(msg % ('cadd', ALGEBRA))
def drem(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return (x * mk_pf(termau_energy(lb, rb)))
    elif ALGEBRA == 'mfe':
        return x + termau_energy(lb, rb)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('drem', ALGEBRA))
def hl(lb:Basic_Subsequence, r:Basic_Subsequence, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return (scale(((2 + r.j) - r.i)) * mk_pf(hl_energy(r)))
    elif ALGEBRA == 'mfe':
        return hl_energy(r)
    elif ALGEBRA == 'count':
        return 1
    else:
        raise ValueError(msg % ('hl', ALGEBRA))
def il(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale(((((2 + lr.j) - lr.i) + rr.j) - rr.i)) * x) * mk_pf(il_energy(lr, rr)))
    elif ALGEBRA == 'mfe':
        return x + il_energy(lr, rr)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('il', ALGEBRA))
def incl(x:float):
    if ALGEBRA == 'pfunc':
        return (x * mk_pf(ul_energy()))
    elif ALGEBRA == 'mfe':
        return x + ul_energy()
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('incl', ALGEBRA))
def ml(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale(2) * x) * mk_pf(((ml_energy() + ul_energy()) + termau_energy(lb, rb))))
    elif ALGEBRA == 'mfe':
        return x + ml_energy() + ul_energy() + termau_energy(lb, rb)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('ml', ALGEBRA))
def nil(n:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return 1
    elif ALGEBRA == 'mfe':
        return 0
    elif ALGEBRA == 'count':
        return 1
    else:
        raise ValueError(msg % ('nil', ALGEBRA))
def sadd(lb:Basic_Subsequence, x:float):
    if ALGEBRA == 'pfunc':
        return ((scale(1) * x) * mk_pf(sbase_energy()))
    elif ALGEBRA == 'mfe':
        return x + sbase_energy()
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('sadd', ALGEBRA))
def sr(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
    if ALGEBRA == 'pfunc':
        return ((scale(2) * x) * mk_pf(sr_energy(lb, rb)))
    elif ALGEBRA == 'mfe':
        return x + sr_energy(lb, rb)
    elif ALGEBRA == 'count':
        return x
    else:
        raise ValueError(msg % ('sr', ALGEBRA))
def h(i:[float]) -> [float]:
    if ALGEBRA == 'pfunc':
        if len(i) > 0:
            return np.sum(i)
        else:
            return np.nan
    elif ALGEBRA == 'mfe':
        if len(i) > 0:
            return np.min(i)
        else:
            return np.nan
    elif ALGEBRA == 'count':
        if len(i) > 0:
            return np.sum(i)
        else:
            return 0
    else:
        raise ValueError(msg % ('h', ALGEBRA))


def backtrace(t_0_i:int, t_0_j:int, name:str) -> float:
    if name not in tables:
        print("Do something to reconstruct %s(%i,%i)" % (name, t_0_i,t_0_j))
        globals()["nt_%s" % name](t_0_i,t_0_j,bwdpass=True)
        sys.exit()

    if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
        return tables[name].bt_get(t_0_i, t_0_j)

    if PRINTBTSTACK:
        print("%scall bt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()

    answers = []
    edges = tables[name].backtrace.loc[t_0_i, t_0_j]
    rep = ""
    if is_not_empty(edges):
        for edge in edges:
            if edge['nt'] is not None:
                params = [backtrace(edge['i'],edge['j'],edge['nt']) if p == 'x' else p for p in edge['algparams']]
                # there can arise situations in the fwd pass where input is split into two ore more non-terminals (like cadd(x,y))
                # x(i,k) evaluates to a valid parse, but y(k,j) don't. Then, x(i,k) tables are filled AND backtrace information,
                # which will somewhere end into a nan edge. We have to prune params here for this case
                if all([is_not_empty(p) for p in params]):
                    algfct = edge['algfct']
                    if algfct is None:
                        res = params[0]
                    else:
                        res = algfct(*params)
                    answers.append(res)
                    #rep += ', ' + str(params)

    eval = h(answers)
    if PRINTBTSTACK:
        decr()
        print("%s} set bt_%s(%i,%i) = %s %s" % (INDENT, name, t_0_i, t_0_j, eval, rep))
    tables[name].bt_set( t_0_i, t_0_j, eval)
    return tables[name].bt_get(t_0_i, t_0_j)
