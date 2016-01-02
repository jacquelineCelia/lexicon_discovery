#include "gmm.h"
#include "toolkit.h"

GMM::GMM(Config* config) {
    _config = config;
}

GMM::GMM(Config* config, int mix_num) {
    _config = config;
    _mix_num = mix_num;
    _weight.resize(_mix_num, 0);
    for (int i = 0; i < _mix_num; ++i) {
        Mixture mixture(_config);
        _mixtures.push_back(mixture);
    }
}

void GMM::load(ifstream& fin) {
    fin.read(reinterpret_cast<char*> (&_mix_num), sizeof(int));
    _weight.resize(_mix_num, 0);
    for (int i = 0; i < _mix_num; ++i) {
        Mixture mixture(_config);
        _mixtures.push_back(mixture);
    }
    fin.read(reinterpret_cast<char*> (&_weight[0]), sizeof(float) * \
            _mix_num);
    for (int m = 0; m < _mix_num; ++m) {
        float det;
        vector<float> mean(_config -> dim(), 0);  
        vector<float> pre(_config -> dim(), 0);
        fin.read(reinterpret_cast<char*> (&det), sizeof(float));
        fin.read(reinterpret_cast<char*> (&mean[0]), sizeof(float) * \
                _config -> dim());
        fin.read(reinterpret_cast<char*> (&pre[0]), sizeof(float) * \
                _config -> dim());
        _mixtures[m].set_det(det);
        _mixtures[m].set_mean(mean);
        _mixtures[m].set_pre(pre);
    }
}

void GMM::set_mixture_det(int s, float det) {
    _mixtures[s].set_det(det);
}

void GMM::set_mixture_mean(int s, vector<float>& mean) {
    _mixtures[s].set_mean(mean);
}

void GMM::set_mixture_pre(int s, vector<float>& pre) {
    _mixtures[s].set_pre(pre);
}


float GMM::ComputeLikehood(float* data) {
    vector<float> likelihood;
    for (int i = 0; i < _mix_num; ++i) {
        likelihood.push_back(_weight[i] + _mixtures[i].likelihood(data));
    } 
    return ToolKit::SumLogs(likelihood);
}

void GMM::ComputeLikehood(vector<float*> data, float* likelihood) {
    for (int i = 0; i < (int) data.size(); ++i) {
        likelihood[i] = ComputeLikehood(data[i]);
    }
}

vector<float> GMM::ComponentLikelihood(float* data) {
    vector<float> likelihood;
    for (int i = 0; i < _mix_num; ++i) {
        likelihood.push_back(_weight[i] + _mixtures[i].likelihood(data));
    } 
    return likelihood;
}

void GMM::Minus(float* data, int index) {
    --_weight[index];
    _mixtures[index].Minus(data);
}

void GMM::Plus(float* data, int index) {
    ++_weight[index];
    _mixtures[index].Plus(data);
}

void GMM::print() {
    for (int i = 0; i < _mix_num; ++i) {
        cout << "weight[" << i << "]: " << _weight[i] << endl;
        _mixtures[i].print();
    }
}

void GMM::Save(ofstream& fout) {
    fout.write(reinterpret_cast<char*> (&_mix_num), sizeof(int));
    fout.write(reinterpret_cast<char*> (&_weight[0]), sizeof(float) * _mix_num);
    for (int m = 0; m < _mix_num; ++m) {
       float det = mixture(m).det();
       vector<float> mean = mixture(m).mean();
       vector<float> pre = mixture(m).pre();
       fout.write(reinterpret_cast<char*> (&det), sizeof(float));
       fout.write(reinterpret_cast<char*> (&mean[0]), sizeof(float) * mean.size());
       fout.write(reinterpret_cast<char*> (&pre[0]), sizeof(float) * pre.size()); 
    }
}

GMM::~GMM() {
}
