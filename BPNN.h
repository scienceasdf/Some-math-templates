#ifndef BPNN_H
#define BPNN_H

#include <vector>
#include <cmath>
#include <boost/serialization/vector.hpp>


inline double sigmoid(double x)
{
    return 1.0/(1.0+exp(-x));
}

class BP_neural_networks
{
public:

    BP_neural_networks() {}

    BP_neural_networks(std::vector<int>& numNodes);
    void train(std::vector<double>& input, std::vector<double>& output);
    void getResult(std::vector<double>& input, std::vector<double>& output);
    void initialWeights();

    std::vector<std::vector<std::vector<double>>> mWeights;
    double learningRate = 1.0;
private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // for boost::archive::binary_iarchive, it may not work
        // I don't know the reasons
        // And not suite for xml archive
        
        ar & learningRate;
        ar & mNumLayers;
        ar & mWeights;
        ar & mNumNodes;
        ar & mNodeValue;
        ar & mErr;
        
    }

    int mNumLayers;
    std::vector<int> mNumNodes;
    std::vector<std::vector<double>> mNodeValue;
    std::vector<std::vector<double>> mErr;

};

#endif // BPNN_H

