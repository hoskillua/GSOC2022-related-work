#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>
#include <chrono>

#define N 100000

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;

inline double sampledRatio(Epic::Vector_3 p, double r, Epic::Vector_3 a, Epic::Vector_3 b, Epic::Vector_3 c, int res = 5) {
    double R = r * r;
    double in = 0; 
    double acc = 1.0 / res;
    for (double i = 0.334 * acc; i < 1; i += acc)
        for (double j = 0.334 * acc; j < 1 - i; j += acc)
        {
            if ((i * a + j * b + (1 - i - j) * c - p).squared_length() < R)
                in++;
        }
    return in / (res * (res + 1) / 2);
}

inline double faceInclusionRatio(Epic::Vector_3 p, double r, Epic::Vector_3 a, Epic::Vector_3 b, Epic::Vector_3 c)
{
    Epic::Vector_3 m = (a+b+c)/3;
    double d_min = (m - p).squared_length();
    double d_max = d_min;
    double da = (a - p).squared_length();
    double db = (b - p).squared_length();
    double dc = (c - p).squared_length();
    d_max = sqrt(std::max(std::max(da, d_max), std::max(db, dc)));
    d_min = sqrt(std::min(std::min(da, d_min), std::min(db, dc)));

    if (d_max <= r) return 1.0;
    else if (r <= d_min) return 0.0;
    return (r - d_min) / (d_max - d_min);
}


int main(int argc, char* argv[])
{
    srand(time(NULL));
    std::vector <Epic::Vector_3> TP(N * 4);
    std::vector <double> R(N);
    std::vector <double> Ref(N);
    std::vector <double> Sampled(N);
    std::vector <double> MinMax(N);

    for (int i = 0; i < N * 4; i++)
        TP[i] = {rand(), rand(), rand()};

    for (int i = 0; i < N; i++)
        R[i] = rand();

    double err_sample1x1 = 0;
    double err_sample2x2 = 0;
    double err_sample3x3 = 0;
    double err_sample5x5 = 0;
    double err_sample9x9 = 0;
    double err_minmax = 0;
    double time_sample1x1 = 0;
    double time_sample2x2 = 0;
    double time_sample3x3 = 0;
    double time_sample5x5 = 0;
    double time_sample9x9 = 0;
    double time_minmax = 0;

    std::chrono::steady_clock::time_point begin, end;
    for (int i = 0; i < N; i++)
        Ref[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 100);


    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        MinMax[i] = faceInclusionRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3]);
    end = std::chrono::steady_clock::now();
    time_minmax += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_minmax += abs(Ref[i] - MinMax[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 1);
    end = std::chrono::steady_clock::now();
    time_sample1x1 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample1x1 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 2);
    end = std::chrono::steady_clock::now();
    time_sample2x2 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample2x2 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 3);
    end = std::chrono::steady_clock::now();
    time_sample3x3 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample3x3 += abs(Ref[i] - Sampled[i]);


    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 5);
    end = std::chrono::steady_clock::now();
    time_sample5x5 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample5x5 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 9);
    end = std::chrono::steady_clock::now();
    time_sample9x9 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample9x9 += abs(Ref[i] - Sampled[i]);

    std::cout << "N = " << N << std::endl;
    std::cout << "minmax approx:        err = " << err_minmax << ", time = " << time_minmax << std::endl;
    std::cout << "sampled res=1x1(=1):  err = " << err_sample1x1 << ", time = " << time_sample1x1 << std::endl;
    std::cout << "sampled res=2x2(=3):  err = " << err_sample2x2 << ", time = " << time_sample2x2 << std::endl;
    std::cout << "sampled res=3x3(=6):  err = " << err_sample3x3 << ", time = " << time_sample3x3 << std::endl;
    std::cout << "sampled res=5x5(=15): err = " << err_sample5x5 << ", time = " << time_sample5x5 << std::endl;
    std::cout << "sampled res=9x9(=45): err = " << err_sample9x9 << ", time = " << time_sample9x9 << std::endl;

}
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>
#include <chrono>

#define N 100000

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;

inline double sampledRatio(Epic::Vector_3 p, double r, Epic::Vector_3 a, Epic::Vector_3 b, Epic::Vector_3 c, int res = 5) {
    double R = r * r;
    double in = 0; 
    double acc = 1.0 / res;
    for (double i = 0.334 * acc; i < 1; i += acc)
        for (double j = 0.334 * acc; j < 1 - i; j += acc)
        {
            if ((i * a + j * b + (1 - i - j) * c - p).squared_length() < R)
                in++;
        }
    return in / (res * (res + 1) / 2);
}

inline double faceInclusionRatio(Epic::Vector_3 p, double r, Epic::Vector_3 a, Epic::Vector_3 b, Epic::Vector_3 c)
{
    Epic::Vector_3 m = (a+b+c)/3;
    double d_min = (m - p).squared_length();
    double d_max = d_min;
    double da = (a - p).squared_length();
    double db = (b - p).squared_length();
    double dc = (c - p).squared_length();
    d_max = sqrt(std::max(std::max(da, d_max), std::max(db, dc)));
    d_min = sqrt(std::min(std::min(da, d_min), std::min(db, dc)));

    if (d_max <= r) return 1.0;
    else if (r <= d_min) return 0.0;
    return (r - d_min) / (d_max - d_min);
}


int main(int argc, char* argv[])
{
    srand(time(NULL));
    std::vector <Epic::Vector_3> TP(N * 4);
    std::vector <double> R(N);
    std::vector <double> Ref(N);
    std::vector <double> Sampled(N);
    std::vector <double> MinMax(N);

    for (int i = 0; i < N * 4; i++)
        TP[i] = {rand(), rand(), rand()};

    for (int i = 0; i < N; i++)
        R[i] = rand();

    double err_sample1x1 = 0;
    double err_sample2x2 = 0;
    double err_sample3x3 = 0;
    double err_sample5x5 = 0;
    double err_sample9x9 = 0;
    double err_minmax = 0;
    double time_sample1x1 = 0;
    double time_sample2x2 = 0;
    double time_sample3x3 = 0;
    double time_sample5x5 = 0;
    double time_sample9x9 = 0;
    double time_minmax = 0;

    std::chrono::steady_clock::time_point begin, end;
    for (int i = 0; i < N; i++)
        Ref[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 100);


    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        MinMax[i] = faceInclusionRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3]);
    end = std::chrono::steady_clock::now();
    time_minmax += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_minmax += abs(Ref[i] - MinMax[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 1);
    end = std::chrono::steady_clock::now();
    time_sample1x1 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample1x1 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 2);
    end = std::chrono::steady_clock::now();
    time_sample2x2 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample2x2 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 3);
    end = std::chrono::steady_clock::now();
    time_sample3x3 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample3x3 += abs(Ref[i] - Sampled[i]);


    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 5);
    end = std::chrono::steady_clock::now();
    time_sample5x5 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample5x5 += abs(Ref[i] - Sampled[i]);

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++)
        Sampled[i] = sampledRatio(TP[i * 4], R[i], TP[i * 4 + 1], TP[i * 4 + 2], TP[i * 4 + 3], 9);
    end = std::chrono::steady_clock::now();
    time_sample9x9 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

    for (int i = 0; i < N; i++)
        err_sample9x9 += abs(Ref[i] - Sampled[i]);

    std::cout << "N = " << N << std::endl;
    std::cout << "minmax approx:        err = " << err_minmax << ", time = " << time_minmax << std::endl;
    std::cout << "sampled res=1x1(=1):  err = " << err_sample1x1 << ", time = " << time_sample1x1 << std::endl;
    std::cout << "sampled res=2x2(=3):  err = " << err_sample2x2 << ", time = " << time_sample2x2 << std::endl;
    std::cout << "sampled res=3x3(=6):  err = " << err_sample3x3 << ", time = " << time_sample3x3 << std::endl;
    std::cout << "sampled res=5x5(=15): err = " << err_sample5x5 << ", time = " << time_sample5x5 << std::endl;
    std::cout << "sampled res=9x9(=45): err = " << err_sample9x9 << ", time = " << time_sample9x9 << std::endl;

}
