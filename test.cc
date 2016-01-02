#include <iostream>
#include <vector>
#include <map>

#include "interface.h"

using namespace std;

int main() {
    Interface interface;
    if (!interface.load_config("configuration")) {
        cout << "Can't load config file" << endl;
    }
    else {
        cout << "successfully loaded config file!" << endl;
    }
    if (!interface.load_speech_data("/usr/users/chiaying/AG/extend_pycfg/1.input")) {
        cout << "failed" << endl;
    }
    else {
        cout << "success! Loaded speech data!" << endl;
    }

    if (!interface.load_clusters("/usr/users/chiaying/AG/extend_pycfg/dphmm.snapshot", "dphmm", "dphmm.id")) {
        cout << "Can't load dphmm model" << endl;
    }
    else {
        cout << "dphmm model loaded successfully!" << endl;
    }

    /*
    if (!interface.load_clusters("/usr/users/chiaying/AG/extend_pycfg/hdphmm.snapshot", "hdphmm", "dphmm.id")) {
        cout << "Can't load hdphmm model" << endl;
    }
    else {
        cout << "hdphmm model loaded successfully!" << endl;
    }
    */
    // cout << "Showing clusters" << endl;
    // interface.show_clusters();

    // cout << "Showing counters" << endl;
    // interface.show_counters();

    vector<symbol> train = interface.retrieve_symbols(0);
    for (size_t i = 0; i < train.size(); ++i) {
        cout << train[i] << " ";
    }
    cout << endl;

    cout << "Initializing state and mixture seq" << endl;
    interface.initialize_state_mixture_seq();

    // cout << "Showing updated counters" << endl;
    // interface.show_counters();

    cout << "Inserting empty strings" << endl;
    interface.insert_empty_strings();

    vector<symbol> new_train = interface.retrieve_symbols(0);
    for (size_t i = 0; i < new_train.size(); ++i) {
        cout << new_train[i] << " ";
    }
    cout << endl;


    /*
    map<char, vector<int> > rules;

    rules['a'].push_back(5);
    rules['a'].push_back(5);
    rules['b'].push_back(10);

    map<char, vector<int> >::iterator iter = rules.begin();
    for (; iter != rules.end(); ++iter) {
        vector<int>::iterator i_iter = (iter -> second).begin();
        for (; i_iter != (iter -> second).end(); ++i_iter) {
            cout << iter -> first << " " << *i_iter << " ";
        }
    }
    */

    return 0;
}
