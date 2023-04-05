#include "Swarm.h"



double PSO::uniformAB(double A, double B) {
    static std::random_device rd;  
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(A, B);
    return distribution(gen);
}

PSO::PSO() { }

PSO::~PSO() { }

void PSO::init(
    double w,
    double phi_p, 
    double phi_g, 
    const std::vector<std::pair<double, double>> bounds,
    unsigned swarmSize
) {
    // params
    PSO::swarmParams = { w, phi_p, phi_g }; 
    PSO::swarmState.swarmParameters = PSO::swarmParams; 
    // dims calc    
    unsigned dimensions = bounds.size(); 
    PSO::swarmState.dimention = dimensions; 
    // poss init 
    PSO::swarmState.swarmPoses.resize(swarmSize);
    for (int i = 0; i < swarmSize; ++i) {
        PSO::swarmState.swarmPoses[i].resize(dimensions); 
        for (int d = 0; d < dimensions; ++d)
            PSO::swarmState.swarmPoses[i][d] = uniformAB(bounds[d].first, bounds[d].second); 
    }
    // best_poss init
    PSO::swarmState.swarmBestPoses.resize(swarmSize);
    for (int i = 0; i < swarmSize; ++i) {
        PSO::swarmState.swarmBestPoses[i].resize(dimensions); 
        for (int d = 0; d < dimensions; ++d)
            PSO::swarmState.swarmBestPoses[i][d] = PSO::swarmState.swarmPoses[i][d]; 
    }
    // velsinit 
    PSO::swarmState.swarmVelos.resize(swarmSize);
    for (int i = 0; i < swarmSize; ++i) {
        PSO::swarmState.swarmVelos[i].resize(dimensions);
        for (int d = 0; d < dimensions; ++d) 
            PSO::swarmState.swarmVelos[i][d] = uniformAB(
                -1 * std::abs( bounds[d].first - bounds[d].second ),
                1 * std::abs( bounds[d].first - bounds[d].second )
            );
    }
    
    PSO::swarmState.bounds.resize(dimensions); 
    for (int d = 0; d < dimensions; ++d) {
        PSO::swarmState.bounds[d].first = bounds[d].first;
        PSO::swarmState.bounds[d].second = bounds[d].second;
    }
    // swarm size init 
    PSO::swarmState.swarmSize = swarmSize; 

    PSO::swarmState.globalBestPos.resize(dimensions); 

    PSO::swarmState.swarmBestVals.resize(swarmSize); 

    PSO::swarmState.curNumIterations = 0; 

    PSO::swarmState.globalBestVal =  HUGE_VAL;
}

void PSO::makeStep(
    const std::function<double(const std::vector<double>&)> &targFunc,
    std::ofstream& outfile
) {
    double r_g = uniformAB(0.0, 1.0);

    for (int i = 0; i < PSO::swarmState.swarmSize; ++i) { 
        if (targFunc(PSO::swarmState.swarmPoses[i]) < targFunc(PSO::swarmState.swarmBestPoses[i])) {
            for (int d = 0; d < PSO::swarmState.dimention; ++d)
                 PSO::swarmState.swarmBestPoses[i][d] = PSO::swarmState.swarmPoses[i][d]; 
            if (targFunc(PSO::swarmState.swarmBestPoses[i]) < targFunc(PSO::swarmState.globalBestPos)) 
                for (int d = 0; d < PSO::swarmState.dimention; ++d)
                    PSO::swarmState.globalBestPos[d] = PSO::swarmState.swarmBestPoses[i][d];
        }

        for (int d = 0; d < PSO::swarmState.dimention; ++d) 
            outfile << PSO::swarmState.swarmPoses[i][d] << ' '; 
        outfile << '\n'; 

        double r_p = uniformAB(0.0, 1.0); 
        bool ok = true;
        for (int d = 0; d < PSO::swarmState.dimention; ++d) {
            PSO::swarmState.swarmVelos[i][d] = PSO::swarmState.swarmParameters.w * PSO::swarmState.swarmVelos[i][d] +
                PSO::swarmState.swarmParameters.phi_p * r_p * (PSO::swarmState.swarmBestPoses[i][d] - PSO::swarmState.swarmPoses[i][d]) + 
                PSO::swarmState.swarmParameters.phi_g * r_g * (PSO::swarmState.globalBestPos[d] - PSO::swarmState.swarmPoses[i][d]); 
            PSO::swarmState.swarmPoses[i][d] += PSO::swarmState.swarmVelos[i][d]; 
            ok = ok && PSO::swarmState.bounds[d].first < PSO::swarmState.swarmPoses[i][d] && PSO::swarmState.bounds[d].second > PSO::swarmState.swarmPoses[i][d];
        }
        if (!ok)
            for (int d = 0; d < PSO::swarmState.dimention; ++d) 
                PSO::swarmState.swarmPoses[i][d] = PSO::swarmState.bounds[d].first + (PSO::swarmState.bounds[d].second - PSO::swarmState.bounds[d].first) * uniformAB(0.0, 1.0); 
    }
    PSO::swarmState.curNumIterations++; 
}

 void PSO::run(
    const std::function<double(const std::vector<double>&)> &targFunc,
    unsigned numOfIterations,
    std::ofstream& outfile
) {
    // std::ofstream outfile("./_dump.txt");  
    while (PSO::swarmState.curNumIterations < numOfIterations)
        PSO::makeStep(targFunc, outfile);
    outfile.close();

    PSO::swarmState.globalBestVal = targFunc(PSO::swarmState.globalBestPos); 
    // file.open(); 
    outfile << "bounds:\n"; 
    for (int d = 0; d < PSO::swarmState.dimention; ++d)
        outfile << "    " << PSO::swarmState.bounds[d].first << ' ' << PSO::swarmState.bounds[d].second << '\n'; 
    outfile << "num of particles: " << std::to_string(PSO::swarmState.swarmSize) << '\n'; 
    outfile << "num of iterations: " << std::to_string(PSO::swarmState.curNumIterations) << '\n'; 
    outfile << "global minimum in ["; 
    for (int d = 0; d < PSO::swarmState.dimention; ++d) 
        outfile << PSO::swarmState.globalBestPos[d] << ' ';
    outfile << "]\n";
    outfile << "MIN = "; 
    outfile << PSO::swarmState.globalBestVal << '\n'; 
    outfile.close();
 }

void PSO::makeVTKsnapshot(
    const std::function<double(const std::vector<double>&)> &targFunc,
    unsigned snapshotNum,
    std::string folderName
) {
    auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto smth = vtkSmartPointer<vtkDoubleArray>::New();
    smth->SetName("function");

    std::vector tmpVec = {0.0, 0.0, 0.0};
    tmpVec[0] = PSO::swarmState.bounds[0].first; 
    tmpVec[1] = PSO::swarmState.bounds[1].first; 
    double stepX = (PSO::swarmState.bounds[0].second - PSO::swarmState.bounds[0].first) / 1000; 
    double stepY = (PSO::swarmState.bounds[1].second - PSO::swarmState.bounds[1].first) / 1000; 
    for (unsigned i = 0; i < 1000; ++i) {
        dumpPoints->InsertNextPoint(tmpVec[0], tmpVec[1], tmpVec[2]);
        smth->InsertNextValue(targFunc(tmpVec)); 
        tmpVec[0] += stepX; 
        tmpVec[1] += stepY;
    }

    unstructuredGrid->SetPoints(dumpPoints); 
    unstructuredGrid->GetPointData()->AddArray(smth);

    // for (unsigned i = 0; i < PSO::swarmState.swarmSize; ++i) {
    //     for (auto position : PSO::swarmState.swarmPoses) {
    //         auto particlePoint = vtkSmartPointer<vtkPoints>::New();
    //         unstructuredGrid->InsertNextCell(position[0], position[1], position[2]);
    //     }
    // }

    std::string fileName = folderName + "task-step-" + std::to_string(snapshotNum) + ".vtu";
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void PSO::visualize(
    const std::function<double(const std::vector<double>&)> &targFunc,
    unsigned numOfIterations
) {
    for (unsigned snapshotNum = 0; snapshotNum < numOfIterations; ++snapshotNum) 
        PSO::makeVTKsnapshot(targFunc, snapshotNum, "./vtu_frames/task-step-");
}

void PSO::run_and_visualize(
    const std::function<double(const std::vector<double>&)> &targFunc,
    unsigned numOfIterations,
    std::ofstream& outfile
) {
    PSO::run(targFunc, numOfIterations, outfile);
    PSO::visualize(targFunc, numOfIterations); 
}

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
    return (a+2)*(a+2) + (b-3)*(b-3);
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

int main(int argc, char* argv[]) {
    const double PI = std::acos(-1.0); 

    std::ofstream oufile("./run/test.txt"); 

    PSO pso; 
    pso.init(
        0.6, 0.1, 0.3,
        {{-1.5, 4.0}, {-3.0, 4.0}},
        3000
    );
    pso.run_and_visualize(mccormick, 50, oufile); 
// 
    // pso.init(
    //     0.6, 0.5, 0.6,
    //     {{-4, 4}, {-4, 4}},
    //     1000
    // );
    // pso.run_and_visualize(test, 100); 

    // pso.init(
    //     0.3, 0.7, 0.5,
    //     {{-PI, PI}, {-PI, PI}},
    //     1000
    // );
    // // pso.run(test1, 40); 
    // pso.run_and_visualize(test1, 40); 

// //
//     pso.init(
//         0.7, 0.4, 0.6,
//         {{0, std::acos(-1.0)}, {0, std::acos(-1.0)}},
//         500
//     );
//     pso.run_and_visualize(michalewicz, 100); 
//
    return 0; 
}

