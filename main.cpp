#define GNUPLOT             "gnuplot"
#define GNUPLOT_WIN_WIDTH   1280
#define GNUPLOT_WIN_HEIGHT  720

#include <iostream>
#include <cmath>
#include <vector>
#include "MultiplicativePRNG.h"

using namespace std;


double calcDiscrepancy(const vector<double> &real, const vector<double> &result)
{
    double maxVal = -1.0;
    int size = real.size();
    for (int i = 0; i < size; i++)
    {
        maxVal = max(abs(real[i] - result[i]), maxVal);
    }
    return maxVal;
}


double calcNext(const vector<vector<double>> &coefsMtr, const vector<vector<double>> &markovMtr,
                const vector<double> &vecF, const vector<double> &markovVec, const vector<double> &h,
                const vector<int> &trajectory)
{
    int indexPrev;
    int indexCur;
    int trajLength = trajectory.size();
    double result;
    double q;

    indexCur = trajectory[0];
    q = (markovVec[indexCur] > 0.0) ? h[indexCur] / markovVec[indexCur] : 0.0;
    result = q * vecF[indexCur];
    for (int i = 1; i < trajLength; i++)
    {
        indexPrev = trajectory[i - 1];
        indexCur = trajectory[i];
        q *= (markovMtr[indexPrev][indexCur] > 0.0) ? coefsMtr[indexPrev][indexCur] / markovMtr[indexPrev][indexCur] : 0.0;

        result += q * vecF[indexCur];
    }

    return result;
}

void calcSystem(const vector<vector<double>> &coefsMtr, const vector<double> &vecF, int chainsNum, int chainLength,
                vector<double> &result)
{

    int size = vecF.size();
    PRNG* prng = new MultiplicativePRNG(pow(2, 31), 262147, 262147);
    vector<vector<double>> markovMtr(size, vector<double>(size));
    vector<double> markovVec(size);
    vector<double> h(size);
    vector<vector<int>> trajectory(chainsNum, vector<int>(chainLength));

    for (int i = 0; i < size; i++)
    {
        markovVec[i] = 1.0 / size;

        for (int j = 0; j < size; j++)
        {
            markovMtr[i][j] = 1.0 / size;
        }
    }

    for (int i = 0; i < chainsNum; i++)
    {
        for (int j = 0; j < chainLength; j++)
        {
            trajectory[i][j] = prng->nextInt(0, size - 1);
        }
    }

    for (int i = 0; i < chainsNum; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fill(h.begin(), h.end(), 0.0);
            h[j] = 1.0;
            result[j] += calcNext(coefsMtr, markovMtr, vecF, markovVec, h, trajectory[i]);
        }
    }

    for (int i = 0; i < size; i++)
    {
        result[i] /= chainsNum;
    }

    delete prng;
}

bool checkCoefs(vector<vector<double>> coefs)
{
    double sum;
    int size = coefs.size();

    for (int i = 0; i < size; i++)
    {
        sum = 0.0;
        for (int j = 0; j < size; j++)
        {
            sum += abs(coefs[i][j]);
        }

        if (sum >= 1.0)
        {
            return false;
        }
    }

    return true;
}


int main() {
    FILE* pipe;

#ifdef WIN32
    pipe = _popen(GNUPLOT, "w");
#else
    pipe = popen(GNUPLOT, "w");
#endif

    const double coefsMtrIn[3][3] = {
            {1.2, -0.3, 0.4},
            {0.4, 0.7, -0.2},
            {0.2, -0.3, 0.9}
    };
    const double vecFIn[3] = {-4.0, 2.0, 0.0};
    const double realIn[3] = {-2.8285714286, 5.1428571429, 2.3428571429};
    const int size = 3;

    vector<vector<double>> coefsMtr(size, vector<double>(size));
    vector<double> vec(size);
    vector<double> real(size);
    vector<double> result(size);
    const int chainLengthMin = 5;
    const int chainLengthStep = 1;
    const int chainLengthStepNum = 45;
    const int chainNumMin = 1000;
    const int chainNumStep = 1000;
    const int chainNumStepNum = 20;

    for (int i = 0; i < size; i++)
    {
        vec[i] = vecFIn[i];
        real[i] = realIn[i];
        for (int j = 0; j < size; j++)
        {
            coefsMtr[i][j] = (i == j) ? -(coefsMtrIn[i][j] - 1.0) : -coefsMtrIn[i][j];
        }
    }

    if (!checkCoefs(coefsMtr))
    {
        printf("Incorrect matrix");
        return -1;
    }

    if (pipe != nullptr)
    {
        fprintf(pipe, "set term wxt size %d, %d\n", GNUPLOT_WIN_WIDTH, GNUPLOT_WIN_HEIGHT);
        fprintf(pipe, "set title \"Error of system solution\"\n");
        fprintf(pipe, "set dgrid3d\n");
        fprintf(pipe, "set xlabel \"Markov chain length\"\n");
        fprintf(pipe, "set ylabel \"Number of Markov chains\"\n");
        fprintf(pipe, "set zlabel \"Error\"\n");
        fprintf(pipe, "splot '-' with lines\n");

        int chainLengthCur = chainLengthMin;
        int chainNumCur;
        for (int i = 0; i < chainLengthStepNum; i++)
        {
            chainNumCur = chainNumMin;
            for (int j = 0; j < chainNumStepNum; j++)
            {
                calcSystem(coefsMtr, vec, chainNumCur, chainLengthCur, result);
                fprintf(pipe, "%d %d %lf\n", chainLengthCur, chainNumCur, calcDiscrepancy(real, result));
                chainNumCur += chainNumStep;
            }
            chainLengthCur += chainLengthStep;

            printf("Iteration: %d/%d\n", i + 1, chainLengthStepNum);
        }

        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        cin.clear();
        cin.get();

#ifdef WIN32
        _pclose(pipe);
#else
        pclose(pipe);
#endif
    }
    else
    {
        printf("Could not open pipe");
    }

    return 0;
}