#ifndef BPNN_H
#define BPNN_H

#include <vector>
#include <cmath>

inline double sigmoid(double x)
{
    return 1.0/(1.0+exp(-x));
}

class BP_neural_networks
{
public:
    BP_neural_networks(std::vector<int>& numNodes);
    void train(std::vector<double>& input, std::vector<double>& output);
    void getResult(std::vector<double>& input, std::vector<double>& output);
    void initialWeights();

    std::vector<std::vector<std::vector<double>>> mWeights;

private:
    int mNumLayers;
    std::vector<int> mNumNodes;
    std::vector<std::vector<double>> mNodeValue;
    std::vector<std::vector<double>> mErr;

};

#endif // BPNN_H

