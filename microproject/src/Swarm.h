#pragma once

// STD lib
#include <algorithm>
#include <functional>
#include <random>
#include <fstream>
#include <vector>
#include <filesystem>
// VTK 
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay2D.h>
#include <vtkVertexGlyphFilter.h>



class PSO final {
private:
    double uniformAB(double A, double B); 
    struct pso_params { double w, phi_p, phi_g; };
    struct pso_state {
        unsigned curNumIterations;
        std::vector<double> globalBestPos;
        double globalBestVal; 
        std::vector<std::pair<double, double>> bounds; 
        pso_params swarmParameters; 
        std::vector<std::vector<double>> swarmPoses; 
        std::vector<std::vector<double>> swarmVelos;
        std::vector<std::vector<double>> swarmBestPoses; 
        std::vector<double> swarmBestVals; 
        unsigned swarmSize; 
        unsigned dimention; 
    };
    void makeStep(
        const std::function<double(const std::vector<double>&)>& targFunc,
        std::string visualDirName = ""
    ); 
    void makeVTKsnapshot(
        unsigned snapshotNum, 
        std::string visualDirName
    );     
    
    pso_params swarmParams; 
    pso_state swarmState;

public:
    PSO(); 
    ~PSO(); 
    void init(
        double w,
        double phi_p, 
        double phi_g, 
        std::vector<std::pair<double, double>> bounds,
        unsigned swarmSize
    ); 
    void run(
        const std::function<double(const std::vector<double>&)>& targFunc,
        unsigned numOfIterations,
        std::string resultFilePath,
        std::string visualDirName = ""
    ); 
    void run_and_visualize(
        const std::function<double(const std::vector<double>&)>& targFunc,
        unsigned numOfIterations,
        std::string resultFileName,
        std::string visualDirName
    );
    void visualize(
        const std::function<double(const std::vector<double>&)>& targFunc,
        unsigned numOfIterations,
        std::string visualDirName
    ); 
}; 


