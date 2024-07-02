import re

class Type:
    regex = r''
    dtype = None
    name = "TAnswer"

    def __init__(self):
        raise NotImplementedError()

    def getType(self):
        return self.dtype

    def getRegex(self):
        return self.regex

    def __str__(self):
        return "GAPC-Type[%s]" % self.name

    def setName(self, name: str):
        self.name = name

    def getName(self):
        return self.name


class TypeInt(Type):
    def __init__(self, name: str):
        self.regex = r'(-?\d+)'
        self.dtype = int
        self.name = name

class TypeFloat(Type):
    def __init__(self, name: str):
        self.regex = r'(-?\d+\.?\d*)'
        self.dtype = float
        self.name = name

class TypeMFE(TypeFloat):
    def __init__(self):
        super().__init__('mfe')


class Product:
    __left = None
    __right = None
    __regex = ""
    __dtypes = []
    __dnames = []

    def isSingle(self):
        return self.__right is None

    def chain(product):
        #print("CHAIN", type(product))
        if isinstance(product, Type):
            return {'regexs': product.getRegex(),
                    'dtypes': [product.getType()],
                    'dnames': [product.getName()]}
        elif product.isSingle():
            return {'regexs': product.__left.getRegex(),
                    'dtypes': [product.__left.getType()],
                    'dnames': [product.__left.getName()]}
        else:
            res_left = Product.chain(product.__left)
            res_right = Product.chain(product.__right)
            return {'regexs': '\( %s , %s \)' % (res_left['regexs'], res_right['regexs']),
                    'dtypes': res_left['dtypes'] + res_right['dtypes'],
                    'dnames': res_left['dnames'] + res_right['dnames']}

    def __init__(self, left, right=None):
        if left is None:
            raise ValueError("left component cannot be None!")
        assert isinstance(left, Product) or isinstance(left, Type)

        if right is None:
            if isinstance(left, Product):
                raise ValueError("Cannot make a Product of a single Product!")

        if right is not None:
            assert isinstance(right, Product) or isinstance(right, Type)

        self.__left = left
        self.__right = right

        res = Product.chain(self)
        self.__regex = res['regexs']
        self.__dtypes = res['dtypes']
        self.__dnames = res['dnames']

    def getRegex(self):
        return self.__regex

    def getNames(self):
        return self.__dnames

    def getTypes(self):
        return self.__dtypes

    def parse_lines(self, lines:[str]):
        pattern = re.compile(self.getRegex())

        results = []
        nonhitlines = 0
        for i, answer in enumerate(lines):
            hit = pattern.fullmatch(answer)
            if hit:
                results.append({
                    name: dtype(value)
                    for (name, dtype, value)
                    in zip(self.getNames(), self.getTypes(), hit.groups())})
            else:
                nonhitlines += 1

        return results

# GAPC_TYPES = {
#     'int': (int, '(-?\d+)'),
#     'energy': (float, r'(-?\d+\.?\d*)'),
#     'rna': (str, r'([ACGU -]+)'),
#     'pairs': (str, r'([\|: ]+)'),
#     #'pp_hybrid_sophie': re.compile(r''), # (25682, target 5' U     C       GA      U U 3',            UGCUG UGGCCUU  AGCCCC U ,            ||||| |:|||:|  |||||| : ,            ACGAC AUCGGGA--UCGGGG G , miRNA 3'        A               C U 5')
#
# }

# def product_to_pattern(product: tuple) -> str:
#     if isinstance(product, tuple):
#         assert len(product) == 2
#
#         res = []
#         for part in product:
#             if isinstance(part, tuple):
#                 res.append(product_to_pattern(part[1]))
#             else:
#                 res.append(part[1])
#         return "\( %s \)" % ' , '.join(res)
#     else:
#         return product[1]

# def parse_gapc(answers: [str], product: tuple):
#     pattern = re.compile(product_to_pattern(product))
#
#     results = []
#     for i, answer in enumerate(answers):
#         p = pattern.fullmatch(answer)
#         if p:
#             results.append(list(p.groups()))
#     return results

def mainP():
    p = Product(TypeInt('a'), TypeInt('b'))
    obs = p.parse_lines(["( 23 , -10 )", "( -44 , 10 )"])
    print(obs)
    #t2 = Product(TypeInt(), TypeInt())  # t2 = (GAPC_TYPES['int'], GAPC_TYPES['int'])
    # t3 = GAPProduct((GAPInt, (GAPInt, GAPInt))) # t3 = (GAPC_TYPES['int'], (GAPC_TYPES['int'], GAPC_TYPES['int']))
    # t4 = GAPProduct(((GAPInt, GAPInt), GAPInt)) # t4 = ((GAPC_TYPES['int'], GAPC_TYPES['int']), GAPC_TYPES['int'])
    # #
    # # print(parse_gapc(['23'], product=t1))
    # #
    # # print("fertig")
    # x = TypeInt()
    # print(x)
    #print(t1.getRegex())
    #print(t2.getRegex())
    # #print(t3)
    # #print(t4)
    # print(t1.chain())

if __name__ == '__main__':
    mainP()
