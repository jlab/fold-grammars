#ifndef SHAPES_HH
#define SHAPES_HH

static const char openParen = '[';
static const char closeParen = ']';

// If Rope size is >= 2, returns contents from second to secondlast character
// used to obtain "y" from "[y]"
template<typename X>
inline Rope inner(const rope::Ref<X> &str) {
	rope::Ref<X> &x = const_cast<rope::Ref<X>&>(str);
	typename rope::Ref<X>::iterator it = x.begin();
	Rope res;
	if (str.size() <= 2) {
		res.empty();
	} else {
		++it;
		for (unsigned int i = 1; i+1 < str.size(); ++i) {
			append(res, (char) *it);
			++it;
		}
	}
	return res;
}

#endif
