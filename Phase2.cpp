//Good Programming Practices: https://www.browserstack.com/guide/coding-standards-best-practices
//Sources:
//removing file extension from file name: https://stackoverflow.com/questions/6417817/easy-way-to-remove-extension-from-a-filename
//converting strings to ints and doubles: https://www.geeksforgeeks.org/cpp-program-for-char-to-int-conversion/
//making random unique numbers: https://www.youtube.com/watch?v=wnjYD_euSJQ
//extracting numbers from string: https://codescracker.com/cpp/program/cpp-extract-numbers-from-string.htm
//dynamic 2d arrays: https://www.geeksforgeeks.org/how-to-declare-a-2d-array-dynamically-in-c-using-new-operator/
//remove the space from nextLine: https://stackoverflow.com/questions/23834624/remove-first-and-last-character-c
//new random number algorithm: https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution

#include <iostream>
#include <fstream>
#include <string>
#include <stdbool.h>
#include <random>
#include <vector>

int dimensionality;
int numberOfPoints;
int iterations;
int numberOfClusters;
double convergenceThreshold;
int runs;

class Point {
private:
    double* pointData = new double[dimensionality];
    void stringToDoubleStar(std::string string) {
        double stringDouble = 0;
        for (int i = 0; i < dimensionality; i++) {
            stringDouble = std::stod(string);
            pointData[i] = std::stod(string);
            //have to remove entered data point from nextLine
            if (string.find(" ") != std::string::npos) {
                string = string.substr(string.find(" "));
                string.erase(0, 1);
            }
        }
    }
public:
    Point(std::string nextLine) {
        stringToDoubleStar(nextLine);
    }

    double getValue(int dimension) {
        return pointData[dimension];
    }
};

class Cluster {
private:
    std::vector<double> centroid;
    std::vector<Point> clusteredPoints;
public:
    Cluster(std::vector<double> newCentroid) {
        centroid = newCentroid;
    }
    void addPoint(Point point) {
        clusteredPoints.push_back(point);
    }
    void removeAllPoints() {
        clusteredPoints.clear();
    }
    Point getPoint(int pointCount) {
        return clusteredPoints[pointCount];
    }
    double getPointValue(int pointCount, int dimension) {
        return clusteredPoints[pointCount].getValue(dimension);
    }
    double getCentroidValue(int dimension) {
        return centroid[dimension];
    }
    int getNumberOfPointsInCluster() {
        return clusteredPoints.size();
    }
    void setCentroid(std::vector<double> newCentroid) {
        centroid = newCentroid;
    }
};

class KMeans {
private:
    std::vector<Point> allPoints;
    std::vector<Cluster> clusters;
    std::string outFileName;
public:
    KMeans(std::vector<Point> mainPoints, std::string outFileName) {
        allPoints = mainPoints;
        this->outFileName = outFileName;
        run();
    }
    void run() {
        int bestRun = 0;
        double bestSSE = 0;
        std::vector<double> SSE;
        std::ofstream outfile(outFileName);
        for (int i = 0; i < runs; i++) {
            SSE.clear();
            clusters.clear();
            initClusters();
            outfile << "Run " + std::to_string(i + 1) + "\n-----\n";
            int iteration = 0;
            while (iteration < 2) {
                for (int j = 0; j < numberOfPoints; j++) {
                    int nearestCluster = getNearestCluster(allPoints[j]);
                    clusters[nearestCluster].addPoint(allPoints[j]);
                }
                resetCentroid();
                outfile << "Iteration " + std::to_string(iteration + 1) + ": SSE = ";
                SSE.push_back(getSSE());
                outfile << std::to_string(SSE[iteration]) + "\n";
                iteration++;
            }
            while (iteration < iterations && (abs(SSE[iteration - 2] - SSE[iteration - 1]) / SSE[iteration - 2]) > convergenceThreshold) {
                for (int j = 0; j < numberOfPoints; j++) {
                    int nearestCluster = getNearestCluster(allPoints[j]);
                    clusters[nearestCluster].addPoint(allPoints[j]);
                }
                resetCentroid();
                outfile << "Iteration " + std::to_string(iteration + 1) + ": SSE = ";
                SSE.push_back(getSSE());
                outfile << std::to_string(SSE[iteration]) + "\n";
                for (int j = 0; j < numberOfClusters; j++) {
                    clusters[j].removeAllPoints();
                }
                iteration++;
            }
            if (i == 0) {
                bestSSE = SSE[iteration - 1];
                bestRun = 1;
            }
            else if (SSE[iteration - 1] < bestSSE) {
                bestSSE = SSE[iteration - 1];
                bestRun = iteration;
            }
        }
        outfile << "\nBest Run: " + std::to_string(bestRun) + ": SSE = " + std::to_string(bestSSE);
    }
    void initClusters() {
        int* randomCenters = new int[numberOfClusters];
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, numberOfPoints - 1);
        bool unique;
        int newRandom;
        for (int i = 0; i < numberOfClusters; i++) {
            do {
                newRandom = distrib(gen);
                unique = true;
                for (int j = 0; j < i; j++)
                    if (randomCenters[j] == newRandom) {
                        unique = false;
                    }
            } while (!unique);
            randomCenters[i] = newRandom;
            std::vector<double> centroid;
            for (int j = 0; j < dimensionality; j++) {
                centroid.push_back(allPoints[newRandom].getValue(j));
            }
            clusters.push_back(centroid);
        }
    }

    int getNearestCluster(Point point) {
        double distance;
        double minimumDistance;
        double sum = 0;
        double x;
        int nearestCluster = 0;
        for (int i = 0; i < dimensionality; i++) {
            x = clusters[0].getCentroidValue(i) - point.getValue(i);
            x = x * x;
            sum += x;
        }
        minimumDistance = sum;
        for (int i = 1; i < numberOfClusters; i++) {
            sum = 0;
            for (int j = 0; j < dimensionality; j++) {
                x = clusters[i].getCentroidValue(j) - point.getValue(j);
                x = x * x;
                sum += x;
            }
            if (minimumDistance > sum) {
                minimumDistance = sum;
                nearestCluster = i;
            }
        }
        return nearestCluster;
    }

    void resetCentroid() {
        for (int i = 0; i < numberOfClusters; i++) {
            std::vector<double> newCentroid;
            int numberOfPointsInCluster = clusters[i].getNumberOfPointsInCluster();
            for (int j = 0; j < dimensionality; j++) {
                double sum = 0;
                if (numberOfPointsInCluster > 0) {
                    for (int k = 0; k < numberOfPointsInCluster; k++) {
                        sum += clusters[i].getPointValue(k, j);
                    }
                    newCentroid.push_back(sum / numberOfPointsInCluster);
                }
            }
            if (numberOfPointsInCluster > 0) {
                clusters[i].setCentroid(newCentroid);
            }
        }
    }

    double getSSE() {
        double SSE = 0;
        double difference = 0;
        for (int i = 0; i < numberOfClusters; i++) {
            for (int k = 0; k < clusters[i].getNumberOfPointsInCluster(); k++) {
                for (int j = 0; j < dimensionality; j++) {
                    difference = clusters[i].getCentroidValue(j) - clusters[i].getPointValue(k, j);
                    SSE += difference * difference;
                }
            }
        }
        return SSE;
    }
};

int main(int argc, char* argv[])
{
    //initialize arguments
    std::string fileName = argv[1];
    numberOfClusters = std::stoi(argv[2]);
    iterations = std::stoi(argv[3]);
    convergenceThreshold = std::atof(argv[4]);
    runs = std::stoi(argv[5]);

    //open file and get point info
    std::ifstream infile(fileName);
    std::string firstLine;
    getline(infile, firstLine);
    size_t lastindex = fileName.find_last_of(".");
    std::string outFileName = fileName.substr(0, lastindex) + "_output.txt";

    //get numberOfPoints and dimensionality
    numberOfPoints = std::stoi(firstLine);
    dimensionality = std::stoi(firstLine.substr(firstLine.find(" ")));

    //get individual data points
    std::string nextLine;
    std::vector<Point> points;
    for (int i = 0; i < numberOfPoints; i++) {
        std::getline(infile, nextLine);
        Point point(nextLine);
        points.push_back(point);
    }
    infile.close();
    KMeans kmeans(points, outFileName);
}