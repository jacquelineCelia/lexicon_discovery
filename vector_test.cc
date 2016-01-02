#include <vector>
#include <iostream>

using namespace std;

int main() {
    vector<int> v(10, 0);
    vector<int>::iterator iter = v.begin();
    int counter = 7;
    for (; iter != v.end(); ++iter) {
        if (counter % 5 == 3) {
            iter = v.insert(iter, -1);
        }
        ++counter;
    }

    cout << v.size() << endl;

    for (size_t i = 0; i < v.size(); ++i) {
        cout << v[i] << " ";
    }
    cout << endl;

    return 0;
}
