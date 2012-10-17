#ifndef PKISS_TABLES_HH
#define PKISS_TABLES_HH

class ThreeD_private;
class ThreeD {
  private:
    ThreeD_private *d;

    ThreeD(const ThreeD &);
    ThreeD &operator=(const ThreeD &);
  public:
    ThreeD();
    ~ThreeD();

    void set(int i, int j, int l, int k, int mfe, int n);

    const std::pair<int, int> &get(int i, int j, int l, int n);

};

#endif
