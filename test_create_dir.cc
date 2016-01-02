#include <string>
#include <iostream>
#include <sys/stat.h>

using namespace std;

int main() {
    string output_dir = "dir";
    mkdir(output_dir.c_str(), S_IRWXU | S_IRWXO | S_IRWXG);
    return 0;
}
