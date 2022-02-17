from pylib.gapc import *

def init(inputsequence, algebra='pfunc', printstack=False, printBTstack=False,
         tabulateNTs=["hairpin","leftB","multiloop","rightB","stack",
                      "dangle","iloop","ml_comps","ml_comps1","strong","struct","weak"]):
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
        if algebra in ['pfunc', 'count']:
            tables['struct'].bt_set(0,0,1.0)
        elif algebra == 'mfe':
            tables['struct'].bt_set(0,0,0.0)

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

def nt_dangle(t_0_i:int, t_0_j:int, name="dangle") -> float:
    computed()
    if name in tables:
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
        add_trace(tables, 'strong', t_0_i,t_0_j, 'dangle', t_0_i, t_0_j, algfct=drem, algparams=[ret_1, 'x', ret_3])

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
def nt_hairpin(t_0_i:int, t_0_j:int, name="hairpin") -> float:
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
        add_trace(tables, 'hairpin', None, None, None, t_0_i, t_0_j, algfct=hl, algparams=[ret_1, ret_2, ret_3])

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
def nt_iloop(t_0_i:int, t_0_j:int, name="iloop") -> float:
    computed()
    if name in tables:
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
                                        add_trace(tables, 'strong', t_0_i,t_0_j, 'iloop', t_0_k_0, t_0_k_1, algfct=il, algparams=[ret_1, ret_2, 'x', ret_4, ret_5])

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
def nt_leftB(t_0_i:int, t_0_j:int, name="leftB") -> float:
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
                                add_trace(tables, 'strong', t_0_i,t_0_j, 'leftB', t_0_k_0, t_0_j-1, algfct=bl, algparams=[ret_1,ret_2,'x', ret_4])

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
def nt_ml_comps(t_0_i:int, t_0_j:int, name="ml_comps") -> float:
    computed()
    if name in tables:
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
        add_trace(tables, 'ml_comps', t_0_i,t_0_j, 'ml_comps', t_0_i+1,t_0_j, algfct=sadd, algparams=[ret_1, 'x'])

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
                    add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])
                    add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x), y), algparams=['x', ret_6])


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
def nt_ml_comps1(t_0_i:int, t_0_j:int, name="ml_comps1") -> float:
    computed()
    if name in tables:
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
        add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps1', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'])


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
                    add_trace(tables, 'ml_comps1', t_0_i,t_0_j, 'ml_comps1', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])
                    add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_0, algfct=lambda x,y: cadd(incl(x),y), algparams=['x', ret_6])

            t_0_k_0 += 1

    ret_7 = np.nan
    if (((t_0_j - t_0_i) >= 5)):
        ret_8 = nt_dangle(t_0_i, t_0_j)
        if (is_not_empty(ret_8)):
            ret_7 = incl(ret_8)

    if (is_not_empty(ret_7)):
        answers.append(ret_7)
        add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_j, algfct=incl, algparams=['x'])


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
                    add_trace(tables, 'dangle', t_0_i,t_0_j, 'ml_comps1', t_0_i, t_0_k_1, algfct=lambda x,y: addss(incl(x),y), algparams=['x', ret_12])

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
def nt_multiloop(t_0_i:int, t_0_j:int, name='multiloop') -> float:
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
        add_trace(tables, 'ml_comps', t_0_i,t_0_j, 'multiloop', t_0_i+1, t_0_j-1, algfct=ml, algparams=[ret_1, 'x', ret_3])

    #    print("multiloop(%i,%i) = " % (t_0_i,t_0_j), answers, ret_0, ret_1, ret_2, ret_3)
    #print(" set nt_multiloop(%i,%i)" % (t_0_i, t_0_j))
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
def nt_rightB(t_0_i:int, t_0_j:int, name="rightB") -> float:
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
                                add_trace(tables, 'strong', t_0_i,t_0_j, 'rightB', t_0_i+1, t_0_k_0, algfct=br, algparams=[ret_1, 'x', ret_3, ret_4])

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
def nt_stack(t_0_i:int, t_0_j:int, name="stack") -> float:
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
        add_trace(tables, 'weak', t_0_i,t_0_j, 'stack', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_1, 'x', ret_3])

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
def nt_strong(t_0_i:int, t_0_j:int, name="strong") -> float:
    computed()
    if name in tables:
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
        add_trace(tables, 'weak', t_0_i,t_0_j, 'strong', t_0_i+1, t_0_j-1, algfct=sr, algparams=[ret_2, 'x', ret_4])


    ret_5 = np.nan
    if (allowLonelyBasepairs(t_0_seq, t_0_i, t_0_j, True)):
        ret_5 = nt_weak(t_0_i, t_0_j)

    if (is_not_empty(ret_5)):
       answers.append(ret_5)
       add_trace(tables, 'weak', t_0_i,t_0_j, 'strong', t_0_i, t_0_j, algfct=None, algparams=['x'])


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
def nt_struct(t_0_i:int, name="struct") -> float:
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
        add_trace(tables, 'struct', t_0_i,t_0_j, 'struct', t_0_i+1, t_0_j, algfct=sadd, algparams=[ret_1, 'x'])


    if (((t_0_right_most - t_0_i) >= 5)):
        t_0_k_0 = (t_0_i + 5)
        while (t_0_k_0 <= t_0_right_most):
            ret_5 = nt_struct(t_0_k_0)
            if (is_not_empty(ret_5)):
                ret_4 = nt_dangle(t_0_i, t_0_k_0)
                if (is_not_empty(ret_4)):
                    ret_0 = cadd(ret_4, ret_5)
                    answers.append(ret_0)
                    add_trace(tables, 'dangle', t_0_i,t_0_j, 'struct', t_0_i, t_0_k_0, algfct=cadd, algparams=['x', ret_5])
                    add_trace(tables, 'struct', t_0_i,t_0_j, 'struct', t_0_k_0, t_0_j, algfct=cadd, algparams=[ret_4, 'x'])


            t_0_k_0 += 1

    ret_6 = np.nan
    if ((((t_0_right_most - t_0_i) >= 0) and ((t_0_right_most - t_0_i) <= 0))):
        ret_7 = LOC(t_0_seq, t_0_i, t_0_i)
        if (is_not_empty(ret_7)):
            ret_6 = nil(ret_7);

    if (is_not_empty(ret_6)):
        answers.append(ret_6)
        add_trace(tables, 'struct', None, None, None, t_0_i, t_0_j, algfct=nil, algparams=[ret_7])


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
def nt_weak(t_0_i:int, t_0_j:int, name="weak") -> float:
    computed()
    if name in tables:
        if (tables[name].is_tabulated(t_0_i, t_0_j)):
           return tables[name].get(t_0_i, t_0_j)

    if PRINTSTACK:
        print("%scall nt_%s(%i,%i) {" % (INDENT, name, t_0_i, t_0_j))
        incr()
    answers = []
    ret_1 = nt_stack(t_0_i, t_0_j)
    if (is_not_empty(ret_1)):
        answers.append(ret_1)
        if 'stack' in tables:
            add_trace(tables, 'stack', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_2 = nt_hairpin(t_0_i, t_0_j)
    if (is_not_empty(ret_2)):
        answers.append(ret_2)
        add_trace(tables, 'hairpin', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_3 = nt_leftB(t_0_i, t_0_j)
    if (is_not_empty(ret_3)):
        answers.append(ret_3)
        add_trace(tables, 'leftB', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_4 = nt_rightB(t_0_i, t_0_j)
    if (is_not_empty(ret_4)):
        answers.append(ret_4)
        add_trace(tables, 'rightB', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_5 = nt_iloop(t_0_i, t_0_j)
    if (is_not_empty(ret_5)):
        answers.append(ret_5)
        add_trace(tables, 'iloop', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

    ret_6 = nt_multiloop(t_0_i, t_0_j)
    if (is_not_empty(ret_6)):
        answers.append(ret_6)
        add_trace(tables, 'multiloop', t_0_i,t_0_j, 'weak', t_0_i, t_0_j, algfct=None, algparams=['x'])

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
                if all([is_not_empty(p) for p in params]):
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
