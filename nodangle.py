from pylib.gapc import *

def init(inputsequence, printstack=False, printBTstack=False, taball=False):
    gapcrna.librna_read_param_file(None)
    global PRINTSTACK
    PRINTSTACK = printstack
    global PRINTBTSTACK
    PRINTBTSTACK = printBTstack
    global TABALL
    TABALL = taball
    global INDENT
    INDENT = ""

    global t_0_seq
    t_0_seq = inputsequence.upper().replace('A','\1').replace('C','\2').replace('G','\3').replace('U','\4')

    global tables
    tables = dict()
    for nt in [
        "hairpin","leftB","multiloop","rightB","stack",
        "dangle","iloop","ml_comps","ml_comps1","strong","struct","weak",
        ]:
        tables[nt] = DPtable(len(t_0_seq), nt)

def incr():
    global INDENT
    INDENT = INDENT + " "
def decr():
    global INDENT
    INDENT = INDENT[:-1]

def nt_dangle(t_0_i:int, t_0_j:int, name="dangle") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
        return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
    ret_2 = np.nan
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
        tables['strong'].add_trace2(t_0_i,t_0_j, 'dangle', t_0_i, t_0_j, algfct=drem, algparams=[ret_1, 'x', ret_3])

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)
def nt_hairpin(t_0_i:int, t_0_j:int, name="hairpin") -> float:
    if (((t_0_j - t_0_i) < 5)):
        return float_zero

    if TABALL:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall_nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
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
        tables['hairpin'].add_trace2(None, None, None, t_0_i, t_0_j, algfct=hl, algparams=[ret_1, ret_2, ret_3])

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if TABALL:
        tables[name].set( t_0_i, t_0_j, eval)
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_iloop(t_0_i:int, t_0_j:int, name="iloop") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
        return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    if (((t_0_j - t_0_i) >= 9)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = (t_0_i + 2)
            while ((t_0_k_0 <= (t_0_j - 7)) and (t_0_k_0 <= (t_0_i + 31))):
                t_0_k_1 = ((t_0_j - 31)) if (((t_0_j - (t_0_k_0 + 5)) >= 31)) else ((t_0_k_0 + 5))
                while (t_0_k_1 <= (t_0_j - 2)):
                    ret_5 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                    if (is_not_empty(ret_5)):
                        ret_4 = np.nan;
                        if ((maxsize(t_0_seq, t_0_k_1, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_1, (t_0_j - 1)))):
                            ret_4 = REGION(t_0_seq, t_0_k_1, (t_0_j - 1))
                        if (is_not_empty(ret_4)):
                            ret_2 = np.nan
                            if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
                                ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)

                            if (is_not_empty(ret_2)):
                                ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                                if (is_not_empty(ret_1)):
                                    ret_3 = nt_strong(t_0_k_0, t_0_k_1)
                                    if (is_not_empty(ret_3)):
                                        res = il(ret_1, ret_2, ret_3, ret_4, ret_5)
                                        answers.append(res)
                                        tables['strong'].add_trace2(t_0_i,t_0_j, 'iloop', t_0_k_0, t_0_k_1, algfct=il, algparams=[ret_1, ret_2, 'x', ret_4, ret_5])

                    t_0_k_1 += 1
                t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)
def nt_leftB(t_0_i:int, t_0_j:int, name="leftB") -> float:
    if (((t_0_j - t_0_i) < 8)):
       return float_zero;

    if TABALL:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    if (((t_0_j - t_0_i) >= 8)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = t_0_i + 2
            while ((t_0_k_0 <= (t_0_j - 6)) and (t_0_k_0 <= (t_0_i + 31))):
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                if (is_not_empty(ret_4)):
                    ret_2 = np.nan
                    if ((maxsize(t_0_seq, (t_0_i + 1), t_0_k_0, 30) and unpaired(t_0_seq, (t_0_i + 1), t_0_k_0))):
                        ret_2 = REGION(t_0_seq, (t_0_i + 1), t_0_k_0)
                    if (is_not_empty(ret_2)):
                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                        if (is_not_empty(ret_1)):
                            ret_3 = nt_strong(t_0_k_0, (t_0_j - 1))
                            if (is_not_empty(ret_3)):
                                res = bl(ret_1, ret_2, ret_3, ret_4)
                                answers.append(res)
                                tables['strong'].add_trace2(t_0_i,t_0_j, 'leftB', t_0_k_0, t_0_j-1, algfct=bl, algparams=[ret_1,ret_2,'x', ret_4])

                t_0_k_0 += 1
    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if TABALL:
        tables[name].set( t_0_i, t_0_j, eval)
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_ml_comps(t_0_i:int, t_0_j:int, name="ml_comps") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
        return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 11)):
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
        if (is_not_empty(ret_1)):
            ret_2 = nt_ml_comps((t_0_i + 1), t_0_j)
            if (is_not_empty(ret_2)):
                ret_0 = sadd(ret_1, ret_2)

    if (is_not_empty(ret_0)):
        answers.append(ret_0)
        tables['ml_comps'].add_trace2(t_0_i,t_0_j, 'ml_comps', t_0_i+1,t_0_j, algfct=sadd, algparams=[ret_1, 'x'])

    if (((t_0_j - t_0_i) >= 10)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= (t_0_j - 5)):
            ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
            if (is_not_empty(ret_6)):
                ret_4 = np.nan
                if (((t_0_k_0 - t_0_i) >= 5)):
                    ret_5 = nt_dangle(t_0_i, t_0_k_0)
                    if (is_not_empty(ret_5)):
                        ret_4 = incl(ret_5)
                if (is_not_empty(ret_4)):
                    res = cadd(ret_4, ret_6)
                    answers.append(res)
                    tables['ml_comps1'].add_trace2(t_0_i,t_0_j, 'ml_comps', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])
                    tables['dangle'].add_trace2(t_0_i,t_0_j, 'ml_comps', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x), y), algparams=['x', ret_6])


            t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)
def nt_ml_comps1(t_0_i:int, t_0_j:int, name="ml_comps1") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
        return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()

    answers = []
    ret_0 = np.nan
    if (((t_0_j - t_0_i) >= 6)):
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
        if (is_not_empty(ret_1)):
            ret_2 = nt_ml_comps1((t_0_i + 1), t_0_j)
            if (is_not_empty(ret_2)):
                ret_0 = sadd(ret_1, ret_2)

    if (is_not_empty(ret_0)):
        answers.append(ret_0)
        tables['ml_comps1'].add_trace2(t_0_i,t_0_j, 'ml_comps1', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'])


    if (((t_0_j - t_0_i) >= 10)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= (t_0_j - 5)):
            ret_6 = nt_ml_comps1(t_0_k_0, t_0_j)
            if (is_not_empty(ret_6)):
                ret_4 = np.nan
                if (((t_0_k_0 - t_0_i) >= 5)):
                    ret_5 = nt_dangle(t_0_i, t_0_k_0)
                    if (is_not_empty(ret_5)):
                        ret_4 = incl(ret_5)

                if (is_not_empty(ret_4)):
                    res = cadd(ret_4, ret_6)
                    answers.append(res)
                    tables['ml_comps1'].add_trace2(t_0_i,t_0_j, 'ml_comps1', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])
                    tables['dangle'].add_trace2(t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x),y), algparams=['x', ret_6])

            t_0_k_0 += 1

    ret_7 = np.nan
    if (((t_0_j - t_0_i) >= 5)):
        ret_8 = nt_dangle(t_0_i, t_0_j)
        if (is_not_empty(ret_8)):
            ret_7 = incl(ret_8)

    if (is_not_empty(ret_7)):
        answers.append(ret_7)
        tables['dangle'].add_trace2(t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_j, algfct=incl, algparams=['x'])


    if (((t_0_j - t_0_i) >= 6)):
        t_0_k_1 = (t_0_i + 5)
        while (t_0_k_1 <= (t_0_j - 1)):
            ret_12 = np.nan
            if (unpaired(t_0_seq, t_0_k_1, t_0_j)):
                ret_12 = REGION(t_0_seq, t_0_k_1, t_0_j)
            if (is_not_empty(ret_12)):
                ret_10 = np.nan
                if (((t_0_k_1 - t_0_i) >= 5)):
                    ret_11 = nt_dangle(t_0_i, t_0_k_1)
                    if (is_not_empty(ret_11)):
                        ret_10 = incl(ret_11)
                if (is_not_empty(ret_10)):
                    res = addss(ret_10, ret_12)
                    answers.append(res)
                    tables['dangle'].add_trace2(t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_1, algfct=lambda x,y: addss(incl(x),y), algparams=['x', ret_12])

            t_0_k_1 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)
def nt_multiloop(t_0_i:int, t_0_j:int, name='multiloop') -> float:
    if (((t_0_j - t_0_i) < 12)):
        return float_zero;

    if TABALL:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
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
        tables['ml_comps'].add_trace2(t_0_i,t_0_j, 'multiloop', t_0_i+1, t_0_j-1, algfct=ml, algparams=[ret_1, 'x', ret_3])

    #    print("multiloop(%i,%i) = " % (t_0_i,t_0_j), answers, ret_0, ret_1, ret_2, ret_3)
    #print(" set nt_multiloop(%i,%i)" % (t_0_i, t_0_j))
    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if TABALL:
        tables[name].set( t_0_i, t_0_j, eval)
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_rightB(t_0_i:int, t_0_j:int, name="rightB") -> float:
    if (((t_0_j - t_0_i) < 8)):
        return float_zero

    if TABALL:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    if (((t_0_j - t_0_i) >= 8)):
        if (basepair(t_0_seq, t_0_i, t_0_j)):
            t_0_k_0 = ((t_0_j - 31)) if (((t_0_j - (t_0_i + 6)) >= 31)) else ((t_0_i + 6))
            while (t_0_k_0 <= (t_0_j - 2)):
                ret_4 = BASE(t_0_seq, (t_0_j - 1), t_0_j)
                if (is_not_empty(ret_4)):
                    ret_3 = np.nan
                    if ((maxsize(t_0_seq, t_0_k_0, (t_0_j - 1), 30) and unpaired(t_0_seq, t_0_k_0, (t_0_j - 1)))):
                        ret_3 = REGION(t_0_seq, t_0_k_0, (t_0_j - 1))

                    if (is_not_empty(ret_3)):
                        ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
                        if (is_not_empty(ret_1)):
                            ret_2 = nt_strong((t_0_i + 1), t_0_k_0)
                            if (is_not_empty(ret_2)):
                                res = br(ret_1, ret_2, ret_3, ret_4)
                                answers.append(res)
                                tables['strong'].add_trace2(t_0_i,t_0_j, 'rightB', t_0_i+1, t_0_k_0, algfct=br, algparams=[ret_1, 'x', ret_3, ret_4])

                t_0_k_0 += 1

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if TABALL:
        tables[name].set( t_0_i, t_0_j, eval)
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_stack(t_0_i:int, t_0_j:int, name="stack") -> float:
    if (((t_0_j - t_0_i) < 7)):
       return float_zero;

    if TABALL:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
            return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
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
        tables['weak'].add_trace2(t_0_i,t_0_j, 'stack', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_1, 'x', ret_3])

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    if TABALL:
        tables[name].set( t_0_i, t_0_j, eval)
        return tables[name].get(t_0_i, t_0_j)
    else:
        return eval
def nt_strong(t_0_i:int, t_0_j:int, name="strong") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
       return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_0 = np.nan
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, False)):
        ret_1 = np.nan
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
        tables['weak'].add_trace2(t_0_i,t_0_j, 'strong', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_2, 'x', ret_4])


    ret_5 = np.nan
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
        ret_5 = nt_weak(t_0_i, t_0_j)

    if (is_not_empty(ret_5)):
       answers.append(ret_5)
       tables['weak'].add_trace2(t_0_i,t_0_j, 'strong', t_0_i, t_0_j, algfct=None, algparams=['x'])


    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)
def nt_struct(t_0_i:int, name="struct") -> float:
    t_0_j = 0
    t_0_right_most = len(t_0_seq)
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
        return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_right_most))
        incr()
    answers = []
    ret_0 = np.nan
    if (((t_0_right_most - t_0_i) >= 1)):
        ret_1 = np.nan
        if (unpaired(t_0_seq, t_0_i, (t_0_i + 1))):
            ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1))
        if (is_not_empty(ret_1)):
            ret_2 = nt_struct((t_0_i + 1))
            if (is_not_empty(ret_2)):
                ret_0 = sadd(ret_1, ret_2)
    if (is_not_empty(ret_0)):
        answers.append(ret_0)
        tables['struct'].add_trace2(t_0_i,t_0_j, 'struct', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'])


    if (((t_0_right_most - t_0_i) >= 5)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= t_0_right_most):
            ret_5 = nt_struct(t_0_k_0)
            if (is_not_empty(ret_5)):
                ret_4 = nt_dangle(t_0_i, t_0_k_0)
                if (is_not_empty(ret_4)):
                    ret_0 = cadd(ret_4, ret_5)
                    answers.append(ret_0)
                    tables['dangle'].add_trace2(t_0_i,t_0_j, 'struct', t_0_i, t_0_k_0, algfct=cadd, algparams=['x', ret_5])
                    tables['struct'].add_trace2(t_0_i,t_0_j, 'struct', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])


            t_0_k_0 += 1

    ret_6 = np.nan
    if ((((t_0_right_most - t_0_i) >= 0) and ((t_0_right_most - t_0_i) <= 0))):
        ret_7 = LOC(t_0_seq, t_0_i, t_0_i)
        if (is_not_empty(ret_7)):
            ret_6 = nil(ret_7);

    if (is_not_empty(ret_6)):
        answers.append(ret_6)
        tables['struct'].add_trace2(None, None, None, t_0_i, t_0_j, algfct=nil, algparams=[ret_7])


    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_right_most, eval))
    tables[name].set( t_0_i, 0, eval)
    return tables[name].get(t_0_i, 0)
def nt_weak(t_0_i:int, t_0_j:int, name="weak") -> float:
    if (tables[name].is_tabulated(t_0_i, t_0_j)):
       return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_1 = nt_stack(t_0_i, t_0_j)
    if (is_not_empty(ret_1)):
        answers.append(ret_1)
        tables['stack'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_2 = nt_hairpin(t_0_i, t_0_j)
    if (is_not_empty(ret_2)):
        answers.append(ret_2)
        tables['hairpin'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_3 = nt_leftB(t_0_i, t_0_j)
    if (is_not_empty(ret_3)):
        answers.append(ret_3)
        tables['leftB'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_4 = nt_rightB(t_0_i, t_0_j)
    if (is_not_empty(ret_4)):
        answers.append(ret_4)
        tables['rightB'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_5 = nt_iloop(t_0_i, t_0_j)
    if (is_not_empty(ret_5)):
        answers.append(ret_5)
        tables['iloop'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_6 = nt_multiloop(t_0_i, t_0_j)
    if (is_not_empty(ret_6)):
        answers.append(ret_6)
        tables['multiloop'].add_trace2(t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    eval = h(answers)
    if PRINTSTACK:
        decr()
        print("%s} set nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].set( t_0_i, t_0_j, eval)
    return tables[name].get(t_0_i, t_0_j)

global ALGEBRA
ALGEBRA = 'pfunc'
if ALGEBRA == 'pfunc':
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

elif ALGEBRA == 'mfe':
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
elif ALGEBRA == 'count':
    # mfe algebra
    def addss(x:float, r:Basic_Subsequence):
        return x
    def bl(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        #print("calling bl(<%i,%i>, <%i,%i>, %f, <%i,%i>)" % (lb.i,lb.j,lr.i,lr.j,x,rb.i,rb.j))
        return x
    def br(lb:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return x
    def cadd(x:float, y:float):
        #print("Calling cadd(x=%s,y=%s)" % (x,y))
        return x * y
    def drem(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return x
    def hl(lb:Basic_Subsequence, r:Basic_Subsequence, rb:Basic_Subsequence):
        return 1
    def il(lb:Basic_Subsequence, lr:Basic_Subsequence, x:float, rr:Basic_Subsequence, rb:Basic_Subsequence):
        return x
    def incl(x:float):
        return x
    def ml(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
        return x
    def nil(n:Basic_Subsequence):
        return 1
    def sadd(lb:Basic_Subsequence, x:float):
        return x
    def sr(lb:Basic_Subsequence, x:float, rb:Basic_Subsequence):
       return x
    def h(i:[float]) -> [float]:
        if len(i) > 0:
            return np.sum(i)
        else:
            return 0

def backtrace(t_0_i:int, t_0_j:int, name:str) -> float:
    if (tables[name].bt_is_tabulated(t_0_i, t_0_j)):
        return tables[name].bt_get(t_0_i, t_0_j)

    if PRINTBTSTACK:
        print("%scall bt_nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()

    answers = []
    edges = tables[name].backtrace.loc[t_0_i, t_0_j]
    if is_not_empty(edges):
        for edge in edges:
            if edge['nt'] is not None:
                params = [backtrace(edge['i'],edge['j'],edge['nt']) if p == 'x' else p for p in edge['algparams']]
                # there can arise situations in the fwd pass where input is split into two ore more non-terminals (like cadd(x,y))
                # x(i,k) evaluates to a valid parse, but y(k,j) don't. Then, x(i,k) tables are filled AND backtrace information,
                # which will somewhere end into a nan edge. We have to prune params here for this case
                params = [p for p in params if is_not_empty(p)]
                if params != []:
                    algfct = edge['algfct']
                    if algfct is None:
                        res = params[0]
                    else:
                        res = algfct(*params)
                    answers.append(res)

    eval = h(answers)
    if PRINTBTSTACK:
        decr()
        print("%s} set bt_nt_%s(%i,%i) = %s" % (INDENT, name, t_0_i, t_0_j, eval))
    tables[name].bt_set( t_0_i, t_0_j, eval)
    return tables[name].bt_get(t_0_i, t_0_j)
