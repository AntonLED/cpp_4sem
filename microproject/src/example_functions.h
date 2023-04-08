#pragma once

namespace  exampleFunctions {
    double mccormick(const std::vector<double> & x) {
        auto a = x[0];
        auto b = x[1];
        return sin(a + b) + (a - b) * (a - b) + 1.0 + 2.5 * b - 1.5 * a;
    }

    double test1(const std::vector<double> & x) {
        return x[1]*x[1] + std::sin(x[0]); 
    }

    double test(const std::vector<double> & x) {
        auto a = x[0];
        auto b = x[1];
        return (a+2)*(a+2) + (b-2)*(b-2);
    }

    double michalewicz(const std::vector<double> & x) {
        auto m = 10;
        auto d = x.size();
        auto sum = 0.0;
        for (int i = 1; i < d; ++i) {
            auto j = x[i - 1];
            auto k = sin(i * j * j / std::acos(-1.0));
            sum += sin(j) * pow(k, (2.0 * m));
        }
        return -sum;
    }
}
