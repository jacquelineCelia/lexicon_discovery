#include <iostream>
#include "bound.h"

using namespace std;

Bound::Bound(int dim, int start_frame, int end_frame) {
    _dim = dim;
    _start_frame = start_frame;
    _end_frame = end_frame;
    _data = NULL;
    _data_ptr.resize(_end_frame - _start_frame + 1, NULL);
}

void Bound::set_data(float** data, int frame_num) {
    _data = data;
    for (int i = 0; i < frame_num; ++i) {
        _data_ptr[i] = _data[i];
    }
}

Bound::~Bound() {
    int frame_num = get_frame_num(); 
    for (int i = 0; i < frame_num; ++i) {
        delete[] _data[i];
    }
    delete[] _data;
}
