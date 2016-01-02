#include <sstream>
#include "segment.h"

Segment::Segment(int id, int frame_num) {
    _id = id;
    stringstream s;
    s << _id;
    symbol sym(s.str());
    _id_symbol = sym;
    _frame_num = frame_num;
}

void Segment::push_back(Bound* bound) {
    vector<float*> data = bound -> data();
    // need to fix the copy thing
    int ptr = _data.size();
    _data.resize(ptr + data.size());
    copy(data.begin(), data.end(), _data.begin() + ptr); 
    _frame_num = (int) _data.size();
}

Segment::~Segment() {
}
