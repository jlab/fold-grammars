import re


class Type:
    regex = r''
    dtypes = []
    name = "TAnswer"
    component_names = []

    def __init__(self, name: str, regex: str, dtypes: [str], component_names: [str]=[]):
        if (name is None) or (name == ""):
            raise ValueError("You need to give the new Type a name!")
        self.name = name

        if ('(' not in regex) or (')' not in regex):
            raise ValueError("Regex: looks like you don't capture any pattern!")
        self.regex = regex

        if not isinstance(dtypes, list):
            dtypes = [dtypes]
        if (dtypes is None) or (len(dtypes) <= 0) or (not isinstance(dtypes, list)):
            raise ValueError("You need to specify at least one dtype as a one element list!")
        if (component_names is None) or (len(component_names) <= 0) or (not isinstance(component_names, list)):
            component_names = [name]
        if (len(component_names) != len(dtypes)):
            raise ValueError("Number of component names and dtypes must be identical!")

        self.dtypes = dtypes
        self.component_names = component_names

    def getTypes(self) -> []:
        return self.dtypes

    def getRegex(self) -> str:
        return self.regex

    def __str__(self) -> str:
        return "GAPC-Type[%s]" % self.name

    def getName(self) -> str:
        return self.name

    def getComponentNames(self) -> [str]:
        return self.component_names


class TypeInt(Type):
    def __init__(self, name: str):
        super().__init__(name,
            r'(-?\d+)',
            int)


class TypeFloat(Type):
    def __init__(self, name: str):
        super().__init__(name,
            r'(-?\d+\.?\d*)',
            float)


class TypeMFE(TypeFloat):
    def __init__(self):
        super().__init__('mfe')


class TypeHybrid(Type):
    def __init__(self):
        super().__init__('hybrid',
            r'\((\d+), target 5\' ([ACGU\- ]+) 3\',           ([ACGU\- ]+),           ([|: ]+),           ([ACGU\- ]+), miRNA 3\'  ([ACGU\- ]+) 5\'\)',
            [int, str, str, str, str, str],
            ['target_position', 'target_unpaired', 'target_stacked', 'pairs', 'mirna_stacked', 'mirna_unpaired'])


class Product:
    __left = None
    __right = None
    __regex = ""
    __dtypes = []
    __component_names = []

    def isSingle(self):
        return self.__right is None

    def chain(product):
        if isinstance(product, Type):
            return {'regexs': product.getRegex(),
                    'dtypes': product.getTypes(),
                    'component_names': product.getComponentNames()}
        elif product.isSingle():
            return {'regexs': product.__left.getRegex(),
                    'dtypes': product.__left.getTypes(),
                    'component_names': product.__left.getComponentNames()}
        else:
            res_left = Product.chain(product.__left)
            res_right = Product.chain(product.__right)
            return {'regexs': '\( %s , %s \)' % (res_left['regexs'], res_right['regexs']),
                    'dtypes': res_left['dtypes'] + res_right['dtypes'],
                    'component_names': res_left['component_names'] + res_right['component_names']}

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
        self.__component_names = res['component_names']

    def getRegex(self):
        return self.__regex

    def getComponentNames(self) -> [str]:
        return self.__component_names

    def getTypes(self) -> []:
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
                    in zip(self.getComponentNames(), self.getTypes(), hit.groups())})
            else:
                nonhitlines += 1

        return results


def mainP():
    p = Product(TypeHybrid())
    obs = p.parse_lines(["(25682, target 5' U     C       GA      U U 3',            UGCUG UGGCCUU  AGCCCC U ,            ||||| |:|||:|  |||||| : ,            ACGAC AUCGGGA--UCGGGG G , miRNA 3'        A               C U 5')"])
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
