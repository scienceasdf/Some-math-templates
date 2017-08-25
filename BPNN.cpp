#include "BPNN.h"
#include <random>

BP_neural_networks::BP_neural_networks(std::vector<int> &numNodes)
{
    mNumLayers=numNodes.size();
    mNumNodes=numNodes;
    mNodeValue.resize(mNumLayers);
    mWeights.resize(mNumLayers);
    mErr.resize(mNumLayers);
    for(int i=0; i<mNumLayers; ++i){
        mNodeValue[i].resize(numNodes[i]);
        mWeights[i].resize(numNodes[i]);
        mErr[i].resize(numNodes[i]);
        for(int j=0; j<numNodes[i]; ++j){
            if(i < mNumLayers - 1) {
                mWeights[i][j].resize(numNodes[i+1]);
                std::fill(mWeights[i][j].begin(),mWeights[i][j].end(),.3);
            }
        }
    }
}

void BP_neural_networks::train(std::vector<double> &input, std::vector<double> &output)
{
    //std::copy(input.begin(),input.end(),mNodeValue[0]);
    for(int i=0; i<mNumNodes[0]; ++i){
        mNodeValue[0][i] = input[i];
    }

    for(int i=0; i<mNumLayers-1; ++i){
        for(int j=0; j<mNumNodes[i+1]; ++j){
            // Calculate the nodes' values on the next layer, following the i'th layer
            double sum = .0;
            for(int k=0; k<mNumNodes[i]; ++k){
                sum += mNodeValue[i][k]*mWeights[i][k][j];
            }
            mNodeValue[i+1][j] = sigmoid(sum);
        }
    }

    for(int j=0; j<mNumNodes[mNumLayers-1]; ++j){
        double e = output[j] - mNodeValue[mNumLayers-1][j];
        mErr[mNumLayers-1][j] =
        mNodeValue[mNumLayers-1][j] * (1.0 - mNodeValue[mNumLayers-1][j]) * e;
        for(int k=0; k<mNumNodes[mNumLayers-2]; ++k){
            mWeights[mNumLayers-2][k][j] +=
            1.0 * mErr[mNumLayers-1][j] * mNodeValue[mNumLayers-2][k];
        }
    }

    for(int i=mNumLayers-1; i>1; --i){
        for(int j=0; j<mNumNodes[i-1]; ++j){
            double g=0;
            for(int k=0; k<mNumNodes[i]; ++k){
                g += mWeights[i-1][j][k] * mErr[i][k];
            }
            mErr[i-1][j] = mNodeValue[i-1][j] * (1.0 - mNodeValue[i-1][j]) * g;

            for(int k=0; k<mNumNodes[i-2]; ++k){
                mWeights[i-2][k][j] += 1.0 * mErr[i-1][j] * mNodeValue[i-2][k];
            }
        }
    }
}

void BP_neural_networks::getResult(std::vector<double>& input, std::vector<double>& output)
{
    //std::copy(input.begin(),input.end(),mNodeValue[0]);
    for(int i=0; i<mNumNodes[0]; ++i){
        mNodeValue[0][i] = input[i];
    }

    for(int i=0; i<mNumLayers-1; ++i){
        for(int j=0; j<mNumNodes[i+1]; ++j){
            // Calculate the nodes' values on the next layer, following the i'th layer
            double sum=.0;
            for(int k=0; k<mNumNodes[i]; ++k){
                sum += mNodeValue[i][k]*mWeights[i][k][j];
            }
            mNodeValue[i+1][j] = sigmoid(sum);
        }
    }

    //std::copy(mNodeValue[mNumLayers-1].begin(), mNodeValue[mNumLayers-1].end(),output);
    for(int i=0; i<mNumNodes[mNumLayers-1]; ++i){
        output[i] = mNodeValue[mNumLayers-1][i];
    }
}

void BP_neural_networks::initialWeights()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    for(int i=0; i<mNumLayers-1; ++i){
        for(int j=0; j<mNumNodes[i]; j++){
            for(int k=0; k<mNumNodes[i+1]; ++k){
                double ub = 1.0 / std::pow(mNumNodes[i],.5);
                std::uniform_real_distribution<double> dis(-ub,ub);
                mWeights[i][j][k]=dis(gen);
            }
        }
    }

}

